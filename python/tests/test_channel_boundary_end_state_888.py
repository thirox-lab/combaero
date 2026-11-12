"""A channel's end state at a node at rest is the duct's own (#888).

A node at rest (Pt = P: a plenum, a combustor, a PressureBoundary or
MassFlowBoundary) holds a state that is not the duct's: a reservoir where
the duct leaves it (P = Pt, T = Tt), a static pressure where it enters. Read
as the duct's end state it put ``mach_in`` at the reservoir density (0.58 in
a 1 mm duct carrying Mach 0.89) and ``mach_out`` above 1 in an unchoked duct.
At such an end the diagnostics now read the Fanno station the duct's own
residual solves (compressible) or the residual's own density state
(incompressible). A moving node (a momentum chamber) keeps its state, and no
residual changes: the heat-transfer correlation was already at the duct's
inlet static state (#636).
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero import _solver_tools
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)
from combaero.network.components import SmoothModel

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
X_AIR = list(cb.mass_to_mole(Y_AIR))
# The case of the issue: air from a reservoir through a short duct to a static
# exit pressure.
TT, PT, D = 808.0, 243.45e3, 0.047
AREA = math.pi * D * D / 4.0


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=60)
    assert r["__success__"], r["__message__"]
    return s, r


def _duct(pr: float, L: float = 1e-3, regime: str = "compressible") -> dict:
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=PT, Tt=TT, Y=Y_AIR))
    g.add_node(PressureBoundary("out", Pt=pr * PT, Tt=TT, Y=Y_AIR))
    g.add_element(
        ChannelElement(
            "c",
            "in",
            "out",
            length=L,
            diameter=D,
            roughness=0.0,
            regime=regime,
            surface=ConvectiveSurface(model=SmoothModel()),
        )
    )
    return _solve(g)[1]["__element_diag__"]["c"]


@pytest.mark.parametrize("pr", [0.600, 0.545, 0.537])
def test_a_short_duct_reports_one_mach_at_both_ends(pr: float) -> None:
    """In a 1 mm duct the entry Mach is the exit Mach; it was 0.58 against
    0.96-1.09 at the ends, from the reservoir state at the inflow and the
    boundary's (Tt, static P) at the outflow."""
    d = _duct(pr)
    assert d["choked"] == 0.0
    assert d["mach_in"] == pytest.approx(d["M_in_duct"], rel=1e-6)
    assert d["mach_out"] == pytest.approx(d["M_exit_duct"], rel=1e-6)
    assert d["mach_out"] < 1.0
    assert abs(d["mach_in"] - d["mach_out"]) < 0.015
    if pr == 0.600:
        assert abs(d["mach_in"] - d["mach_out"]) < 2e-3


def test_a_choked_duct_reports_sonic_at_its_exit() -> None:
    d = _duct(0.500)
    assert d["choked"] == 1.0
    assert d["mach_out"] == pytest.approx(1.0, rel=1e-6)


def test_the_inlet_reference_is_the_ducts_static_state() -> None:
    """The ref block (velocity, rho_ref, mu_ref, ...) is the duct's inlet
    Fanno station, and T_aw the recovery temperature there:
    T + Pr^(1/3) u^2 / (2 cp) at the duct's static state."""
    d = _duct(0.600)
    G = d["m_dot"] / AREA
    st = cb.fanno_state_at_mach(G, TT, d["M_in_duct"], X_AIR)
    assert d["rho_ref"] == pytest.approx(st.rho, rel=1e-6)
    assert d["velocity"] == pytest.approx(st.u, rel=1e-6)
    cs = cb.complete_state(st.T, st.P, X_AIR)
    assert d["mu_ref"] == pytest.approx(cs.transport.mu, rel=1e-9)
    r = cs.transport.Pr ** (1.0 / 3.0)
    T_rec = st.T + r * st.u * st.u / (2.0 * cs.thermo.cp)
    assert d["T_aw"] == pytest.approx(T_rec, abs=0.05)
    assert d["T_aw"] < TT


def test_an_incompressible_duct_reads_its_residuals_density_state() -> None:
    """Incompressible: the inflow end at the approach static P of the
    residual (#651) and the stagnation T (the incompressible convention); the
    outflow end, which takes the flow at its static pressure, is unchanged."""
    d = _duct(0.98, L=0.1, regime="incompressible")
    G = d["m_dot"] / AREA
    a = _solver_tools.approach_static_pressure(PT, TT, G, Y_AIR, True)
    assert d["rho_ref"] == pytest.approx(cb.density(TT, a.P, X_AIR), rel=1e-9)
    assert d["P_in"] == PT  # the node's, as before
    rho_out = cb.density(TT, 0.98 * PT, X_AIR)
    v_out = d["m_dot"] / (rho_out * AREA)
    assert d["mach_out"] == pytest.approx(v_out / cb.speed_of_sound(TT, X_AIR), rel=1e-9)


def _two_ducts(mid: object) -> tuple[NetworkSolver, dict]:
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=PT, Tt=TT, Y=Y_AIR))
    g.add_node(mid)
    g.add_node(PressureBoundary("out", Pt=0.8 * PT, Tt=TT, Y=Y_AIR))
    for name, a, b in (("c1", "in", "p"), ("c2", "p", "out")):
        g.add_element(
            ChannelElement(
                name, a, b, length=0.05, diameter=D, roughness=0.0, regime="compressible"
            )
        )
    return _solve(g)


def test_a_plenum_end_reports_the_ducts_state() -> None:
    """A plenum holds its gas at rest (Pt = P): both ducts at it report
    their own Fanno station there, the case the C3X segments hit."""
    _, r = _two_ducts(PlenumNode("p"))
    d = r["__element_diag__"]
    assert PlenumNode.state_at_rest
    for name in ("c1", "c2"):
        assert d[name]["mach_in"] == pytest.approx(d[name]["M_in_duct"], rel=1e-6)
        assert d[name]["mach_out"] == pytest.approx(d[name]["M_exit_duct"], rel=1e-6)
    assert d["c2"]["mach_in"] > 0.3


def test_a_moving_node_end_keeps_its_node_state() -> None:
    """A momentum chamber carries its own moving static state (Pt > P), here
    at twice the duct's area: the duct leaving it still reads mach_in at the
    chamber's (T, P), not its own inlet station."""
    s, r = _two_ducts(MomentumChamberNode("p", area=2.0 * AREA))
    assert not MomentumChamberNode.state_at_rest
    d = r["__element_diag__"]
    x = np.array(r["__x_solution__"])
    st = s._report_state(s.network.nodes["p"], x)
    assert st.Pt > st.P
    cs = cb.complete_state(st.T, st.P, list(st.X))
    v = d["c2"]["m_dot"] / (cs.thermo.rho * AREA)
    assert d["c2"]["mach_in"] == pytest.approx(v / cs.thermo.a, rel=1e-12)
    assert abs(d["c2"]["mach_in"] - d["c2"]["M_in_duct"]) > 0.01  # not the duct's
    assert d["c1"]["mach_in"] == pytest.approx(d["c1"]["M_in_duct"], rel=1e-6)


def _heated_pair() -> FlowNetwork:
    """A hot duct from a reservoir to a static exit, cooled through a wall by
    a duct fed from a MassFlowBoundary: every end is a flow boundary, the
    heat is in both ducts' momentum (#672)."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("hot_in", Pt=3.0e5, Tt=1200.0, Y=Y_AIR))
    g.add_node(PressureBoundary("hot_out", Pt=2.0e5, Tt=1200.0, Y=Y_AIR))
    g.add_node(MassFlowBoundary("cold_in", m_dot=0.01, Tt=400.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cold_out", Pt=2.0e5, Tt=400.0, Y=Y_AIR))
    for name in ("hot", "cold"):
        g.add_element(
            ChannelElement(
                name,
                f"{name}_in",
                f"{name}_out",
                length=0.4,
                diameter=0.01,
                roughness=0.0,
                regime="compressible",
                surface=ConvectiveSurface(area=math.pi * 0.01 * 0.4),
            )
        )
    g.add_wall(
        ThermalWall(id="w", element_a="hot", element_b="cold", layers=[WallLayer(0.001, 25.0)])
    )
    return g


def test_heated_ducts_report_their_kernel_mach_at_boundary_ends() -> None:
    """Both ends, both ducts: the Fanno-Rayleigh kernel's Mach (the exit at
    the duct's own heated exit Tt), the MassFlowBoundary inlet included."""
    _, r = _solve(_heated_pair())
    d = r["__element_diag__"]
    for name in ("hot", "cold"):
        assert d[name]["rayleigh_in_momentum"] == 1.0
        assert d[name]["mach_in"] == pytest.approx(d[name]["M_in_duct"], rel=1e-6)
        assert d[name]["mach_out"] == pytest.approx(d[name]["M_exit_duct"], rel=1e-6)


def test_the_heated_rows_are_exact() -> None:
    """Every Jacobian row of the heated boundary-end network against central
    differences at the root, columns scaled by |x|: the wall heat rows read
    the ducts' heat-transfer state, which this change leaves alone."""
    s, r = _solve(_heated_pair())
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()
    fd = fd_jacobian(s, x)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    assert float(np.max(np.abs(Js - fds) / scale)) < 1e-5
