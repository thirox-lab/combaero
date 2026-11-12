"""exit_head_recovered: the dynamic head a node at rest silently recovers.

Information only (drive-convention audit, 2026-11-12). A total-receiving
element (channel, pressure loss, area change, lossless link, splitter
reference, vortex) matches the receiving node's Pt; into a node at rest the
exit head 0.5 rho V^2 is recovered with no dump loss. The key reports it, 0
where the receiver carries its own head (a MomentumChamberNode) or the
element loses it (a boundary declaring exit_head_lost).

The networks are the audit's probes S1 (channel -> node -> wall hole) and S4
(last element into an outflow PressureBoundary).
"""

from __future__ import annotations

import math
import warnings

import pytest

import combaero as cb
from combaero.network import (
    AreaChangeElement,
    ChannelElement,
    FlowNetwork,
    LosslessConnectionElement,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
    VortexElement,
)
from combaero.network.pressure_loss import ConstantFractionLoss, ConstantHeadLoss

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
X_AIR = list(cb.species.dry_air())
D = 0.05
A = math.pi * D * D / 4.0
KEY = "exit_head_recovered"


def _solve(net: FlowNetwork) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(net).solve(timeout=30)
    assert r["__success__"], r.get("__message__")
    return r


def _q_isentropic(m: float, area: float, Pt: float, Tt: float) -> float:
    """0.5 rho V^2 of air at stagnation (Pt, Tt) and flux m / area: the
    calorically perfect isentrope at the stagnation gamma, solved by fixed
    point. An independent reference, within O(dgamma) of the C++ state."""
    th = cb.complete_state(Tt, Pt, X_AIR).thermo
    g, cp = th.gamma, th.cp
    R = cp * (g - 1.0) / g
    G = m / area
    M = 0.0
    for _ in range(400):
        T = Tt / (1.0 + 0.5 * (g - 1.0) * M * M)
        rho = Pt * (T / Tt) ** (g / (g - 1.0)) / (R * T)
        V = G / rho
        M = V / math.sqrt(g * R * T)
    return 0.5 * rho * V * V


def _s1(node_kind: str) -> dict:
    """S1: PB 1.3e5 -> channel -> X -> 40 mm wall hole -> PB 1.0e5."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=1.30e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PlenumNode("x") if node_kind == "plenum" else MomentumChamberNode("x", area=A))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(ChannelElement("ch", "in", "x", length=0.5, diameter=D))
    g.add_element(
        OrificeElement("or", "x", "out", diameter=0.04, correlation="fixed", Cd=0.7, in_pipe=False)
    )
    return _solve(g)


def _s4(kind: str) -> dict:
    """S4: PB 1.05e5 -> channel -> a (chamber, area A) -> X -> outflow PB 1.0e5."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=1.05e5, Tt=300.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("a", area=A))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(ChannelElement("c", "in", "a", length=0.5, diameter=D))
    if kind == "area_change":
        g.add_element(AreaChangeElement("x", "a", "out", F0=A, F1=4.0 * A))
    elif kind == "loss":
        g.add_element(
            PressureLossElement("x", "a", "out", ConstantHeadLoss(zeta=0.0, area=A), area=A)
        )
    else:
        g.add_element(ChannelElement("x", "a", "out", length=1e-4, diameter=D))
    return _solve(g)


def test_s1_channel_into_plenum_reports_its_exit_head() -> None:
    r = _s1("plenum")
    m = r["ch.m_dot"]
    q = r["__element_diag__"]["ch"][KEY]
    assert q == pytest.approx(_q_isentropic(m, A, r["x.Pt"], 300.0), rel=2e-3)
    # The audit's 4448 Pa is the same head at the STAGNATION density.
    assert 4400.0 < q < 4700.0


def test_s1_channel_into_momentum_chamber_of_its_area_is_zero() -> None:
    r = _s1("chamber")
    assert r["__element_diag__"]["ch"][KEY] == 0.0


@pytest.mark.parametrize(("kind", "area"), [("loss", A), ("area_change", 4.0 * A)])
def test_s4_total_receiving_element_into_outflow_boundary(kind: str, area: float) -> None:
    """PressureLoss and AreaChange do not read exit_head_lost: the outflow
    boundary recovers their head, at the exit area (F1 for the change)."""
    r = _s4(kind)
    m = r["x.m_dot"]
    q = r["__element_diag__"]["x"][KEY]
    assert q > 0.0
    assert q == pytest.approx(_q_isentropic(m, area, 1.0e5, 300.0), rel=2e-3)


def test_s4_channel_into_boundary_with_exit_head_lost_is_zero() -> None:
    r = _s4("short_channel")
    diag = r["__element_diag__"]
    assert diag["x"][KEY] == 0.0
    # The feeding channel discharges into a momentum chamber: 0 as well.
    assert diag["c"][KEY] == 0.0


def test_channel_into_total_coupled_boundary_reports_the_head() -> None:
    """coupling="total" declares the head recovered: the channel reports it."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=1.05e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR, coupling="total"))
    g.add_element(ChannelElement("c", "in", "out", length=0.5, diameter=D))
    r = _solve(g)
    q = r["__element_diag__"]["c"][KEY]
    assert q == pytest.approx(_q_isentropic(r["c.m_dot"], A, 1.0e5, 300.0), rel=2e-3)


def test_fraction_loss_with_no_velocity_reference_is_zero() -> None:
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=1.05e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PlenumNode("p"))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(
        OrificeElement("o1", "in", "p", diameter=0.02, correlation="fixed", Cd=0.7, in_pipe=False)
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g.add_element(PressureLossElement("pl", "p", "out", ConstantFractionLoss(0.01)))
    r = _solve(g)
    assert r["__element_diag__"]["pl"]["area_source"] == "default"
    assert r["__element_diag__"]["pl"][KEY] == 0.0


def test_lossless_carries_the_head_of_the_chamber_it_leaves() -> None:
    """A lossless link has no area: into a plenum it carries the head of the
    momentum chamber it draws from, 0 from a plenum."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=1.05e5, Tt=300.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("a", area=A))
    g.add_node(PlenumNode("p"))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(ChannelElement("c", "in", "a", length=0.5, diameter=D))
    g.add_element(LosslessConnectionElement("l1", "a", "p"))
    g.add_element(LosslessConnectionElement("l2", "p", "out"))
    r = _solve(g)
    diag = r["__element_diag__"]
    m = r["c.m_dot"]
    assert diag["l1"][KEY] == pytest.approx(_q_isentropic(m, A, r["p.Pt"], 300.0), rel=2e-3)
    assert diag["l2"][KEY] == 0.0


def test_vortex_into_plenum_reports_the_swirl_head() -> None:
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=1.05e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PlenumNode("p"))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(VortexElement("v", "in", "p", r_in=0.05, r_out=0.1, v_theta_in=40.0))
    g.add_element(
        OrificeElement("o", "p", "out", diameter=0.02, correlation="fixed", Cd=0.7, in_pipe=False)
    )
    r = _solve(g)
    d = r["__element_diag__"]["v"]
    v = d["v_theta_out"]
    rho = cb.density(d["T_out"], d["P_swirl_out"], X_AIR)
    # At 20 m/s the static T is within 0.2 K of Tt: rho at Tt is close.
    assert d[KEY] > 0.0
    assert d[KEY] == pytest.approx(0.5 * rho * v * v, rel=2e-3)


def test_incompressible_channel_into_plenum_uses_its_constant_density() -> None:
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=1.02e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PlenumNode("x"))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(ChannelElement("ch", "in", "x", length=0.5, diameter=D, regime="incompressible"))
    g.add_element(
        OrificeElement("or", "x", "out", diameter=0.04, correlation="fixed", Cd=0.7, in_pipe=False)
    )
    r = _solve(g)
    m = r["ch.m_dot"]
    rho0 = cb.density(300.0, r["x.Pt"], X_AIR)
    q = r["__element_diag__"]["ch"][KEY]
    assert q > 0.0
    assert q == pytest.approx(m * m / (2.0 * rho0 * A * A), rel=1e-12)
