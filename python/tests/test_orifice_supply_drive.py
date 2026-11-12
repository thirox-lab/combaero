"""The plain orifice and the connection elements drive Pt -> static (#1087).

The owner's rule (2026-11-12): a hole or orifice drives from the upstream
total to the downstream static. #1122 gave it to the film, effusion and slot
holes; here OrificeElement, EffectiveAreaConnectionElement and
DiameterDischargeCoefficientConnectionElement follow. Into a
MomentumChamberNode their direction comes from Pt_from - P_to: a forward
flow discharges against the chamber's static, a reversed flow draws from
that static and discharges against Pt_from. Inside the band
P_to < Pt_from < Pt_to the orifice blows, where the Pt - Pt direction
reversed it.

Excluded, Pt - Pt exactly as before (OrificeElement._discharges_to_exit_static):
an in-pipe orifice with a velocity of approach in either direction
(beta > 0), one that recovers pressure in a pipe beyond (#645), a
pass-through bore (#788, #1016) and a closed restriction (#1024).

Jacobian checks are central differences with pressure steps sized on the
local drop (#1107, the shared _fd_jacobian helper).
"""

from __future__ import annotations

import copy
import warnings
from collections.abc import Callable

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian, scaled_error

import combaero as cb
from combaero.network import (
    ChannelElement,
    DiameterDischargeCoefficientConnectionElement,
    EffectiveAreaConnectionElement,
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
P_GAS = 5e5


def _orifice() -> OrificeElement:
    return OrificeElement("orf", "cool", "ch", diameter=4e-3)


def _eff_area() -> EffectiveAreaConnectionElement:
    return EffectiveAreaConnectionElement("orf", "cool", "ch", effective_area=8e-6)


def _d_cd() -> DiameterDischargeCoefficientConnectionElement:
    return DiameterDischargeCoefficientConnectionElement("orf", "cool", "ch", 4e-3, Cd=0.7)


MAKERS: list[Callable[[], OrificeElement]] = [_orifice, _eff_area, _d_cd]


def _band_net(element: OrificeElement, dp: float) -> FlowNetwork:
    """Hot gas through a merge chamber at ~5 bar (172 m/s at 2 kg/s, a
    dynamic head of ~3.4 % of P); the orifice's supply plenum at
    5e5 (1 + dp) discharges into it."""
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("gas_in", m_dot=2.0, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("exit", Pt=P_GAS, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cool", Pt=P_GAS * (1.0 + dp), Tt=600.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("ch", area=0.01, main_inlet="feed"))
    g.add_element(LosslessConnectionElement("feed", "gas_in", "ch"))
    g.add_element(LosslessConnectionElement("tail", "ch", "exit"))
    g.add_element(element)
    return g


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r.get("__message__")
    return s, r


def _states(s: NetworkSolver, r: dict) -> tuple:
    x = np.array([r[n] for n in s.unknown_names])
    si = s._get_node_state(s.network.nodes["cool"], x)
    so = s._get_node_state(s.network.nodes["ch"], x)
    si.m_dot = so.m_dot = r["orf.m_dot"]
    return si, so


def _fd_err(s: NetworkSolver, x: np.ndarray) -> float:
    J = s._residuals_and_jacobian(x, compute_jacobian=True)[1].toarray()
    return float(scaled_error(J, fd_jacobian(s, x), x).max())


# --- the band: forward Pt -> static -------------------------------------------


@pytest.mark.parametrize("make", MAKERS)
def test_an_orifice_between_the_chambers_static_and_total_blows(make) -> None:
    s, r = _solve(_band_net(make(), dp=-0.01))
    pt_c = P_GAS * 0.99
    assert r["ch.P"] < pt_c < r["ch.Pt"]  # inside the band
    el = s.network.elements["orf"]
    assert el._hole_cd.exit_static is True
    assert r["orf.m_dot"] > 0.0
    assert r["__element_diag__"]["orf"].get("is_ingesting", 0.0) == 0.0
    # Falsification: the Pt - Pt direction reverses the same states.
    si, so = _states(s, r)
    cd = el._effective_Cd(si, so)
    assert el._flow_and_partials(si, so, cd)[0] > 0.0
    el._hole_cd.exit_static = False
    try:
        assert el._flow_and_partials(si, so, cd)[0] < 0.0
    finally:
        el._hole_cd.exit_static = True


@pytest.mark.parametrize("make", MAKERS)
def test_an_orifice_ingests_only_below_the_chambers_static(make) -> None:
    s, r = _solve(_band_net(make(), dp=-0.05))
    assert r["ch.P"] > P_GAS * 0.95
    assert r["orf.m_dot"] < 0.0


@pytest.mark.parametrize("make", MAKERS)
def test_the_flow_passes_zero_smoothly_at_the_chambers_static(make) -> None:
    # At fixed node states m_calc(Pt_supply) is monotone and continuous
    # through the chamber static: both forms vanish there. The Pt - Pt
    # direction jumped by most of a jet across Pt_chamber.
    s, r = _solve(_band_net(make(), dp=-0.01))
    el = s.network.elements["orf"]
    si, so = _states(s, r)

    def m_at(pt: float) -> float:
        st = copy.copy(si)
        st.P = st.Pt = pt
        return el._flow_and_partials(st, so, 0.7)[0]

    P_m = float(so.P)
    m_ref = abs(m_at(P_m + 2e3))
    steps = (-300.0, -30.0, -3.0, -1.0, -0.3, 0.0, 0.3, 1.0, 3.0, 30.0, 300.0)
    ms = [m_at(P_m + d) for d in steps]
    assert all(b > a for a, b in zip(ms, ms[1:], strict=False))  # monotone through zero
    assert abs(m_at(P_m)) < 1e-12 + 1e-3 * m_ref
    assert max(abs(m) for m in ms[2:-2]) < 0.06 * m_ref  # ~sqrt(3 / 2000)
    el._hole_cd.exit_static = False
    try:
        Pt_m = float(so.Pt)
        assert m_at(Pt_m + 3.0) - m_at(Pt_m - 3.0) > 0.3 * m_ref
    finally:
        el._hole_cd.exit_static = True


# --- the Jacobian ----------------------------------------------------------------


@pytest.mark.parametrize("make", MAKERS)
@pytest.mark.parametrize("dp", [-0.01, -0.05, 0.02])  # band, reverse, plain forward
def test_the_jacobian_is_analytic(make, dp: float) -> None:
    s, r = _solve(_band_net(make(), dp=dp))
    x = np.array(r["__x_solution__"])
    assert _fd_err(s, x) < 1e-5


@pytest.mark.parametrize("make", MAKERS)
@pytest.mark.parametrize("offset", [-0.5, -0.05, 0.05, 0.5])
def test_the_jacobian_is_analytic_across_the_switch(make, offset: float) -> None:
    # The supply moved onto the chamber static, either side of it inside
    # the 0.1 Pa logistic blend, at a state, not a root. The central
    # pressure step is |offset| / 400: m ~ sqrt(dp) there, so a step sized
    # on |P| or on the blend width is truncation (the error falls 4x per
    # halving, second order: the analytic column is right). Exactly AT the
    # static dm/dp of the sqrt law is unbounded and the kernels report 0, a
    # measure-zero point the film holes of #1122 share.
    s, r = _solve(_band_net(make(), dp=-0.01))
    x = np.array(r["__x_solution__"])
    s.network.nodes["cool"].Pt = float(r["ch.P"]) + offset
    J = s._residuals_and_jacobian(x, compute_jacobian=True)[1].toarray()
    fd = np.zeros_like(J)
    for j, name in enumerate(s.unknown_names):
        h = abs(offset) / 400.0 if name.endswith((".P", ".Pt")) else 1e-6 * max(abs(x[j]), 1e-6)
        e = np.zeros_like(x)
        e[j] = h
        fd[:, j] = (np.asarray(s._residuals(x + e)) - np.asarray(s._residuals(x - e))) / (2 * h)
    assert float(scaled_error(J, fd, x).max()) < 1e-5


# --- the second law ---------------------------------------------------------------


@pytest.mark.parametrize("make", MAKERS)
@pytest.mark.parametrize("dp", [-0.01, -0.05])  # band (forward), reverse leg
def test_the_second_law_holds_on_both_legs(make, dp: float) -> None:
    # The control volume runs between the drive pair: forward from the
    # supply's Pt to the chamber static, reversed from that static to the
    # supply's Pt. Closed at the chamber's Pt, the blowing band orifice
    # would destroy entropy.
    _, r = _solve(_band_net(make(), dp=dp))
    v = r["__second_law__"]["orf"]
    assert v["verdict"] == "pass"
    assert v["sigma_norm"] > 0.0
    assert r["__second_law_failures__"] == []


# --- the exclusions: Pt - Pt exactly as before ------------------------------------


def _pipe_net(recover: bool, chamber: bool = True) -> FlowNetwork:
    """A metering orifice between pipes at momentum chambers (beta = 0.5),
    recovering into the downstream pipe (#645) or discharging into a
    chamber with no pipe beyond."""
    g = FlowNetwork()
    node = MomentumChamberNode if chamber else None
    g.add_node(PressureBoundary("in", Pt=1.2e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(node("a", area=7.85e-3))
    g.add_node(node("b", area=7.85e-3))
    g.add_element(ChannelElement("p1", "in", "a", length=0.5, diameter=0.1))
    g.add_element(OrificeElement("orf", "a", "b", diameter=0.05))
    if recover:
        g.add_element(ChannelElement("p2", "b", "out", length=0.5, diameter=0.1))
    else:
        g.add_element(LosslessConnectionElement("p2", "b", "out"))
    return g


@pytest.mark.parametrize("recover", [True, False])
def test_an_in_pipe_orifice_keeps_the_pt_pt_drive(recover: bool) -> None:
    g = _pipe_net(recover)
    g.resolve_all_topology()
    el = g.elements["orf"]
    assert el.beta > 0.0
    assert el._recovery_fwd is recover
    assert el._discharges_to_exit_static() is False
    assert el._hole_cd.exit_static is False
    s = NetworkSolver(g)
    x = np.array(s._build_x0(), dtype=float)
    si = s._get_node_state(g.nodes["a"], x)
    so = s._get_node_state(g.nodes["b"], x)
    assert el._hole_cd.drive_pair(si, so, None) is None


def test_a_reverse_approach_pipe_keeps_the_pt_pt_drive() -> None:
    # Plenum-fed (beta = 0) but a single pipe beyond: reversed flow
    # approaches through it (beta_reverse > 0), excluded.
    g = FlowNetwork()
    g.add_node(PressureBoundary("cool", Pt=1.2e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("ch", area=7.85e-3))
    g.add_element(OrificeElement("orf", "cool", "ch", diameter=0.05))
    g.add_element(ChannelElement("p2", "ch", "out", length=0.5, diameter=0.1))
    g.resolve_all_topology()
    el = g.elements["orf"]
    assert el.beta == 0.0 and el.beta_reverse > 0.0
    assert el._hole_cd.exit_static is False


def test_a_pass_through_and_a_closed_orifice_keep_the_pt_pt_drive() -> None:
    g = _band_net(OrificeElement("orf", "cool", "ch", area=0.01), dp=-0.01)  # = chamber
    g.resolve_all_topology()
    assert g.elements["orf"]._pass_through is True
    assert g.elements["orf"]._hole_cd.exit_static is False
    g = _band_net(EffectiveAreaConnectionElement("orf", "cool", "ch", effective_area=0.0), -0.01)
    g.resolve_all_topology()
    assert g.elements["orf"].closed is True
    assert g.elements["orf"]._hole_cd.exit_static is False


@pytest.mark.parametrize("make", MAKERS)
def test_a_still_receiver_is_unchanged(make) -> None:
    # Into a plenum (P = Pt) there is no drive pair: the original path.
    g = FlowNetwork()
    g.add_node(PressureBoundary("cool", Pt=1.3e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("ch", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(make())
    g.resolve_all_topology()
    el = g.elements["orf"]
    assert el._discharges_to_exit_static() is True
    assert el._hole_cd.exit_static is False
