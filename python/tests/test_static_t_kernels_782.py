"""Kernels that take a STATIC temperature are given one (#782, Tier 2 physics).

A network node carries its STAGNATION temperature (``state.T`` is ``Tt`` at
every interior node). Four C++ inputs are documented as static and form a
density from the static P beside them:

* ``BranchInput.T`` (the compressible tee, ``solver_interface.h``);
* ``station_half_drop`` and ``channel_entry_drop`` (``CrossflowSegmentElement``,
  used by ``EffusionLiner`` and the impingement arrays);
* ``crossflow_velocity_ratio`` (a duct-fed ``EffusionPlateElement``).

Fed Tt, each density was low by Tt/T - 1: about 2 % at M = 0.3, 6-10 % at
M = 0.5-0.6. Each test below solves a network at M ~ 0.5-0.6 and recomputes
the kernel's quantity from the solved state with a static T found HERE, from
h(Tt) = h(T) + (G R T / P)^2 / 2 by bisection -- independently of the C++
closure the elements now use.
"""

from __future__ import annotations

import copy
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    TeeJunctionElement,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_bleed_station_bassett as tb  # noqa: E402
import test_effusion_liner as tl  # noqa: E402

X_AIR = list(cb.species.dry_air())
Y_AIR = list(cb.mole_to_mass(X_AIR))


def _static_T(Tt: float, P: float, G: float, X: list[float]) -> float:
    """Static T of gas at stagnation Tt, static P, mass flux G (bisection)."""
    R = cb.specific_gas_constant(X)
    ht = cb.h_mass(Tt, X)
    lo, hi = 0.3 * Tt, Tt
    for _ in range(200):
        T = 0.5 * (lo + hi)
        if ht - cb.h_mass(T, X) - 0.5 * (G * R * T / P) ** 2 > 0.0:
            lo = T
        else:
            hi = T
    return 0.5 * (lo + hi)


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r.get("__message__")
    return s, r


def _jac_err(s: NetworkSolver, x: np.ndarray) -> float:
    """RESIDUAL-ROW-SCALING: compare J_ij |x_j| (see _build_residual_scales)."""
    _, J = s._residuals_and_jacobian(x)
    rel = np.array([1e-5 if n.endswith("m_dot") else 1e-7 for n in s.unknown_names])
    fd = fd_jacobian(s, x, rel=rel, floor=1e-10)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    return float(np.max(np.abs(Js - fds) / scale))


# ---------------------------------------------------------------------------
# Bleed station: station_half_drop + channel_entry_drop
# ---------------------------------------------------------------------------


def _bleed() -> tuple[NetworkSolver, dict]:
    """Reservoir -> entry -> bleed station b -> outlet, M ~ 0.52 at b."""
    return _solve(tb._bleed_net("compressible", drive=2.0e4, length=0.05))


def test_the_bleed_station_drops_use_the_static_density() -> None:
    """Both segments' momentum terms at b, by hand at rho(P_b, T_static).

    s0 carries the entry acceleration (K_in = 0) and its half of b's station;
    s1 the other half. The station's static T is the one of its mean flow,
    the same for both halves. Fed Tt, both came out 5.5 % low."""
    _, r = _bleed()
    A, kappa = 1e-3, cb.STATION_KAPPA_BLEED_BASSETT
    m_a, m_b = r["s0.m_dot"], r["s1.m_dot"]
    P, Tt = r["b.P"], r["b.Tt"]
    T = _static_T(Tt, P, 0.5 * (m_a + m_b) / A, X_AIR)
    assert Tt / T - 1.0 > 0.05  # the case is a real test
    rho = cb.density(T, P, X_AIR)
    half = (m_b * abs(m_b) - m_a * abs(m_a) + kappa * m_a * (m_a - m_b)) / (2.0 * rho * A * A)
    entry = m_a * m_a / (2.0 * rho * A * A)
    d = r["__element_diag__"]
    assert d["s0"]["dP_momentum"] == pytest.approx(entry + half, rel=1e-9)
    assert d["s1"]["dP_momentum"] == pytest.approx(half, rel=1e-9)


def test_the_bleed_station_jacobian_matches_finite_differences() -> None:
    s, r = _bleed()
    assert _jac_err(s, np.array(r["__x_solution__"], dtype=float)) < 1e-5


# ---------------------------------------------------------------------------
# Duct-fed effusion panel: crossflow_velocity_ratio
# ---------------------------------------------------------------------------


def _liner() -> tuple[FlowNetwork, NetworkSolver, dict]:
    """Four stations, duct Mach ~ 0.58."""
    g = tl._net(tl._liner(4), p_in=1.5e5, p_out=1.1e5)
    s, r = _solve(g)
    return g, s, r


@pytest.mark.slow
def test_the_panel_crossflow_ratio_uses_the_static_duct_state() -> None:
    """McGreehan-Schotsch's U1/Vi: U1 = m_mean / (rho A) and Vi the isentropic
    jet velocity from the duct's STATIC (P, T) to the gas (their extraction
    item 22). Both now see the station's static T."""
    g, _, r = _liner()
    d = r["__element_diag__"]
    gaps = []
    for i in (1, 2, 3, 4):
        plate = g.elements[f"ln__p{i}"]
        node = f"ln__b{i}"
        A = plate.crossflow_area
        m_mean = 0.5 * (r[f"ln__s{i - 1}.m_dot"] + r[f"ln__s{i}.m_dot"])
        P, Tt = r[f"{node}.P"], r[f"{node}.Tt"]
        T = _static_T(Tt, P, abs(m_mean) / A, X_AIR)
        gaps.append(Tt / T - 1.0)
        U1 = abs(m_mean) / (cb.density(T, P, X_AIR) * A)
        w = cb.jet_impulse(1.0, P, T, 1.0e5, X_AIR).J
        Vi = math.sqrt(w * w + cb.CROSSFLOW_VI_FLOOR**2)
        assert d[f"ln__p{i}"]["Cd_inlet_crossflow_U1_over_Vi"] == pytest.approx(U1 / Vi, rel=1e-9)
    assert min(gaps) > 0.05


@pytest.mark.slow
def test_the_panel_crossflow_ratio_partials_match_finite_differences() -> None:
    """d(U1/Vi) in the station's T0 (``{node}.T``), static P and both segment
    flows, now through the static T as well. Checked on the panel's own
    term: the panel row's Pt column is off by ~3 % at this Mach on main as
    well, in the compressible orifice, not in this term."""
    g, s, r = _liner()
    x = np.array(r["__x_solution__"], dtype=float)
    plate = g.elements["ln__p2"]
    node = g.nodes["ln__b2"]
    st_in = s._get_node_state(node, x)
    st_out = s._get_node_state(g.nodes["gas"], x)
    flows = {"ln__s1": r["ln__s1.m_dot"], "ln__s2": r["ln__s2.m_dot"]}
    _, d = plate._hole_cd.supply_crossflow(st_in, st_out, flows)

    def ratio(dP: float = 0.0, dT: float = 0.0, ds1: float = 0.0, ds2: float = 0.0) -> float:
        a = copy.copy(st_in)
        a.__dict__.pop("_X", None)
        a.P += dP
        a.Tt += dT
        f = {"ln__s1": flows["ln__s1"] + ds1, "ln__s2": flows["ln__s2"] + ds2}
        return plate._hole_cd.supply_crossflow(a, st_out, f)[0]

    hP, hT, hm = 1.0, 1e-3, 1e-6
    fd = {
        "ln__b2.P": (ratio(dP=hP) - ratio(dP=-hP)) / (2 * hP),
        "ln__b2.Tt": (ratio(dT=hT) - ratio(dT=-hT)) / (2 * hT),
        "ln__s1.m_dot": (ratio(ds1=hm) - ratio(ds1=-hm)) / (2 * hm),
        "ln__s2.m_dot": (ratio(ds2=hm) - ratio(ds2=-hm)) / (2 * hm),
    }
    for key, v in fd.items():
        assert d[key] == pytest.approx(v, rel=1e-5), key


# ---------------------------------------------------------------------------
# Compressible tee: BranchInput.T
# ---------------------------------------------------------------------------

_A = 2e-3
_D = math.sqrt(4.0 * _A / math.pi)


def _tee(kind: str, port_com: type = MomentumChamberNode) -> FlowNetwork:
    """A tee at M ~ 0.5-0.6 in its common arm, MomentumChamberNode ports
    (the common one optionally a plenum)."""
    g = FlowNetwork()
    if kind == "branching":
        g.add_node(MassFlowBoundary("in", m_dot=0.5, Tt=500.0, Y=Y_AIR))
        g.add_node(PressureBoundary("o_s", Pt=1.5e5, Tt=500.0, Y=Y_AIR))
        g.add_node(PressureBoundary("o_b", Pt=1.5e5, Tt=500.0, Y=Y_AIR))
        edges = (("ch_c", "in", "c"), ("ch_s", "s", "o_s"), ("ch_b", "b", "o_b"))
    else:
        g.add_node(MassFlowBoundary("i_s", m_dot=0.3, Tt=500.0, Y=Y_AIR))
        g.add_node(MassFlowBoundary("i_b", m_dot=0.3, Tt=300.0, Y=Y_AIR))
        g.add_node(PressureBoundary("out", Pt=1.5e5, Tt=500.0, Y=Y_AIR))
        edges = (("ch_s", "i_s", "s"), ("ch_b", "i_b", "b"), ("ch_c", "c", "out"))
    g.add_node(port_com("c", area=_A) if port_com is MomentumChamberNode else port_com("c"))
    g.add_node(MomentumChamberNode("s", area=_A))
    g.add_node(MomentumChamberNode("b", area=_A))
    for eid, a, b in edges:
        g.add_element(ChannelElement(eid, a, b, length=0.05, diameter=_D))
    g.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="c",
            straight_node="s",
            branch_node="b",
            theta=math.pi / 2.0,
            F_C=_A,
            psi=1.0,
            tee_type=kind,
        )
    )
    return g


def _kernel_inputs(
    monkeypatch: pytest.MonkeyPatch, s: NetworkSolver, x: np.ndarray, kind: str
) -> tuple:
    """The (com, str, bra) BranchInputs the tee hands its C++ kernel at x."""
    name = f"compressible_{kind}_tee_rj"
    real = getattr(cb._core, name)
    seen: list[tuple] = []

    def spy(com, st, bra):  # noqa: ANN001, ANN202
        seen.append(tuple((b.P_static, b.T, b.m_dot, b.A) for b in (com, st, bra)))
        return real(com, st, bra)

    monkeypatch.setattr(cb._core, name, spy)
    s._residuals_and_jacobian(x, compute_jacobian=False)
    monkeypatch.setattr(cb._core, name, real)
    assert seen
    return seen[-1]


@pytest.mark.parametrize("kind", ["branching", "merging"])
def test_the_tee_kernel_reads_the_static_temperature_of_a_chamber_port(
    kind: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """At a MomentumChamberNode port the kernel's T is the static T of the
    arm's own flow, so its rho = P / (R T) and u = m / (rho A) are the
    port's. Fed Tt, the common arm's rho was 6 % (branching) and 6 %
    (merging) low."""
    s, r = _solve(_tee(kind))
    x = np.array(r["__x_solution__"], dtype=float)
    inputs = _kernel_inputs(monkeypatch, s, x, kind)
    gap = 0.0
    for node, (P, T, m, A) in zip("csb", inputs, strict=True):
        Tt = r[f"{node}.Tt"]
        assert pytest.approx(r[f"{node}.P"], rel=1e-12) == P
        T_ref = _static_T(Tt, P, abs(m) / A, X_AIR)
        assert pytest.approx(T_ref, rel=1e-9) == T, node
        gap = max(gap, Tt / T - 1.0)
    assert gap > 0.05


def test_a_plenum_port_is_at_rest(monkeypatch: pytest.MonkeyPatch) -> None:
    """A plenum's static state is its total: P = Pt, gas at rest. Its T is
    handed over unchanged."""
    s, r = _solve(_tee("merging", port_com=PlenumNode))
    x = np.array(r["__x_solution__"], dtype=float)
    (_, T_c, _, _), _, _ = _kernel_inputs(monkeypatch, s, x, "merging")
    assert T_c == r["c.Tt"]


@pytest.mark.parametrize("kind", ["branching", "merging"])
def test_the_tee_jacobian_matches_finite_differences(kind: str) -> None:
    s, r = _solve(_tee(kind))
    assert _jac_err(s, np.array(r["__x_solution__"], dtype=float)) < 1e-5
