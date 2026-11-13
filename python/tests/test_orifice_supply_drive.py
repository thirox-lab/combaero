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

The gate (OrificeElement._discharges_to_exit_static) is the #1087 band's
own structure: a SIDE JET into a merge chamber (a MomentumChamberNode with
main_inlet declared, this element not its main inlet). Everything else
keeps Pt - Pt exactly as before: the orifice as the chamber's main inlet, a
chamber without main_inlet, an orifice in line with a pipe, an in-pipe
orifice that recovers pressure in a pipe beyond (ISO permanent loss, #645),
a pass-through bore (#788, #1016) and a closed restriction (#1024). A side
jet's reversed leg draws from the chamber static with no velocity of
approach (_approach_along).

Jacobian checks are central differences with pressure steps sized on the
local drop (#1107, the shared _fd_jacobian helper).
"""

from __future__ import annotations

import copy
import math
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


def _pipe_net(recover: bool, merge: bool = False) -> FlowNetwork:
    """A metering orifice between pipes at momentum chambers (beta = 0.5),
    recovering into the downstream pipe (#645) or discharging into a
    chamber with no pipe beyond. ``merge``: 'b' is a merge chamber whose
    main inlet is another stream, so the orifice is a side jet into it."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=1.2e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("a", area=7.85e-3))
    g.add_node(MomentumChamberNode("b", area=7.85e-3, main_inlet="feed" if merge else None))
    if merge:
        g.add_node(MassFlowBoundary("gas_in", m_dot=0.5, Tt=300.0, Y=Y_AIR))
        g.add_element(LosslessConnectionElement("feed", "gas_in", "b"))
    g.add_element(ChannelElement("p1", "in", "a", length=0.5, diameter=0.1))
    g.add_element(OrificeElement("orf", "a", "b", diameter=0.05))
    if recover:
        g.add_element(ChannelElement("p2", "b", "out", length=0.5, diameter=0.1))
    else:
        g.add_element(LosslessConnectionElement("p2", "b", "out"))
    return g


def test_a_recovering_in_pipe_orifice_keeps_the_pt_pt_drive() -> None:
    # ISO's permanent-loss row (#645) reads both nodes' Pt: excluded.
    g = _pipe_net(recover=True, merge=True)  # a side jet but for recovery
    g.resolve_all_topology()
    el = g.elements["orf"]
    assert el.beta > 0.0 and el._recovery_fwd is True
    assert el._discharges_to_exit_static(g) is False
    assert el._hole_cd.exit_static is False
    s = NetworkSolver(g)
    x = np.array(s._build_x0(), dtype=float)
    si = s._get_node_state(g.nodes["a"], x)
    so = s._get_node_state(g.nodes["b"], x)
    assert el._hole_cd.drive_pair(si, so, None) is None


def test_an_approach_pipe_alone_does_not_exclude_a_side_jet() -> None:
    # beta > 0 without a recovering pipe beyond, into a merge chamber.
    g = _pipe_net(recover=False, merge=True)
    g.resolve_all_topology()
    el = g.elements["orf"]
    assert el.beta > 0.0 and not (el._recovery_fwd or el._recovery_rev)
    assert el._hole_cd.exit_static is True


def _in_line_net(p_exit: float, main_inlet: str | None, area: float | None) -> FlowNetwork:
    """The reviewer's in-line networks (re-review of de3ad53d): the supply
    boundary feeds the orifice (d 50 mm) into a chamber whose only other
    element is a single D 0.1 pipe to the outlet; ``extra`` side elements
    are added by the caller."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=1.2e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("out", Pt=p_exit, Tt=300.0, Y=Y_AIR))
    kw: dict = {"main_inlet": main_inlet}
    if area is not None:
        kw["area"] = area
    g.add_node(MomentumChamberNode("b", **kw))
    g.add_element(OrificeElement("orf", "in", "b", diameter=0.05))
    g.add_element(ChannelElement("p2", "b", "out", length=0.5, diameter=0.1))
    return g


@pytest.mark.parametrize("p_exit", [1.3e5, 1.0e5])  # reversed, forward
@pytest.mark.parametrize(
    ("main_inlet", "area"),
    [(None, None), (None, 7.85e-3), (None, 1.9635e-3), ("orf", None)],
)
def test_an_in_line_or_main_inlet_orifice_is_unchanged(p_exit, main_inlet, area) -> None:
    # Not a side jet into a merge chamber: the Pt - Pt form exactly (no
    # drive pair, the pre-#1087 path; bit-identical to main by F, J and
    # solved-x hashes, tmp/orifdrive/bitid4.py). The in-line patch of
    # de3ad53d moved these by -23 % (chamber area = bore area) and -16 %
    # (main_inlet="orf", inherited area).
    g = _in_line_net(p_exit, main_inlet, area)
    g.resolve_all_topology()
    el = g.elements["orf"]
    assert el._discharges_to_exit_static(g) is False
    assert el._hole_cd.exit_static is False
    s = NetworkSolver(g)
    x = np.array(s._build_x0(), dtype=float)
    si = s._get_node_state(g.nodes["in"], x)
    so = s._get_node_state(g.nodes["b"], x)
    assert el._hole_cd.drive_pair(si, so, None) is None


@pytest.mark.parametrize("extra", ["closed", "tiny"])
def test_an_inert_element_at_the_chamber_does_not_gate_the_main_inlet(extra: str) -> None:
    # A closed orifice or a 0.1 mm side hole into the main_inlet chamber:
    # the main inlet stays Pt - Pt (de3ad53d's predicate flipped on it,
    # +6.8 %). The tiny hole is itself a side jet and is gated; the closed
    # one is not.
    g = _in_line_net(1.3e5, "orf", None)
    g.add_node(PressureBoundary("x", Pt=1.25e5, Tt=300.0, Y=Y_AIR))
    side = (
        OrificeElement("side", "x", "b", area=0.0)
        if extra == "closed"
        else OrificeElement("side", "x", "b", diameter=1e-4, in_pipe=False)
    )
    g.add_element(side)
    g.resolve_all_topology()
    assert g.elements["orf"]._hole_cd.exit_static is False
    assert g.elements["side"]._hole_cd.exit_static is (extra == "tiny")


# --- an approach pipe: beta_reverse > 0 (side injection), beta > 0 -----------


def _side_net(dp: float, tail_channel: bool, in_pipe: bool = True) -> FlowNetwork:
    """A side-injection hole into a duct-fed chamber (the reviewer's case):
    a single channel at the chamber gives the default in_pipe orifice a
    reverse approach (beta_reverse > 0) unless the outlet is a second
    channel. The supply at Pt = 5e5 (1 + dp)."""
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("gas_in", m_dot=2.0, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("exit", Pt=P_GAS, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cool", Pt=P_GAS * (1.0 + dp), Tt=600.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("ch", main_inlet="duct"))
    g.add_element(ChannelElement("duct", "gas_in", "ch", length=0.05, diameter=0.1128))
    if tail_channel:
        g.add_element(ChannelElement("tail", "ch", "exit", length=0.05, diameter=0.1128))
    else:
        g.add_element(LosslessConnectionElement("tail", "ch", "exit"))
    g.add_element(OrificeElement("orf", "cool", "ch", diameter=4e-3, in_pipe=in_pipe))
    return g


def _band_dp(tail_channel: bool) -> float:
    """The supply dp at the middle of the chamber's band (P, Pt)."""
    _, r = _solve(_side_net(-0.01, tail_channel, in_pipe=False))
    return 0.5 * (r["ch.P"] + r["ch.Pt"]) / P_GAS - 1.0


def test_a_side_injection_hole_in_the_band_blows() -> None:
    # Reviewer's repro: in the band the default (in_pipe) hole blew only
    # without an approach pipe; it read -8.06e-4 kg/s with one.
    s, r = _solve(_side_net(-0.01, tail_channel=False))
    el = s.network.elements["orf"]
    assert el.beta_reverse > 0.0 and el._hole_cd.exit_static is True
    assert r["ch.P"] < P_GAS * 0.99 < r["ch.Pt"]
    assert r["orf.m_dot"] > 0.0


@pytest.mark.parametrize("leg", ["band", "reverse"])
@pytest.mark.parametrize("tail_channel", [False, True])
def test_the_outcome_does_not_depend_on_the_chambers_outlet(tail_channel: bool, leg: str) -> None:
    # The outlet element sets beta_reverse (one channel at the chamber) or
    # not (two). Neither leg reads it: forward there is no approach pipe
    # (beta = 0), and the reversed leg draws from the chamber static, which
    # carries no approach head (_approach_along). So in_pipe does not
    # change the root; inside the band the hole blows, below the static it
    # ingests.
    dp = _band_dp(tail_channel)
    if leg == "reverse":
        _, r0 = _solve(_side_net(dp, tail_channel, in_pipe=False))
        dp = (r0["ch.P"] - 0.02 * P_GAS) / P_GAS - 1.0
    _, r_pipe = _solve(_side_net(dp, tail_channel, in_pipe=True))
    _, r_wall = _solve(_side_net(dp, tail_channel, in_pipe=False))
    if leg == "band":
        assert r_pipe["ch.P"] < P_GAS * (1.0 + dp) < r_pipe["ch.Pt"]
        assert r_pipe["orf.m_dot"] > 0.0
    else:
        assert r_pipe["orf.m_dot"] < 0.0
    assert r_pipe["orf.m_dot"] == pytest.approx(r_wall["orf.m_dot"], rel=1e-6)


def _side_pipe_net(dp: float, in_pipe: bool = True) -> FlowNetwork:
    """Side injection with a large reverse approach ratio (beta_reverse 0.4):
    a flowing chamber fed by a lossless link and drained by one 10 mm pipe;
    the orifice (d 4 mm) is NOT in line, so its reversed leg reads no
    approach."""
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("gas_in", m_dot=0.02, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("exit", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cool", Pt=1.0e5 * (1.0 + dp), Tt=300.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("ch", area=math.pi * 0.01**2 / 4.0, main_inlet="feed"))
    g.add_element(LosslessConnectionElement("feed", "gas_in", "ch"))
    g.add_element(ChannelElement("p2", "ch", "exit", length=0.05, diameter=0.01))
    g.add_element(OrificeElement("orf", "cool", "ch", diameter=4e-3, in_pipe=in_pipe))
    return g


def test_a_side_injection_reversed_leg_reads_no_approach_head() -> None:
    # Not in line: the chamber's through-flow is not an approach to the
    # bore, so beta_reverse must not change the root (E and the #651
    # correction would count the pipe's head twice). Reversed below the
    # chamber static.
    _, r0 = _solve(_side_pipe_net(0.0, in_pipe=False))
    dp = (r0["ch.P"] - 2e3) / 1.0e5 - 1.0
    s, r_pipe = _solve(_side_pipe_net(dp, in_pipe=True))
    el = s.network.elements["orf"]
    assert el.beta_reverse == pytest.approx(0.4)
    assert el._hole_cd.exit_static is True
    _, r_wall = _solve(_side_pipe_net(dp, in_pipe=False))
    assert r_pipe["orf.m_dot"] < 0.0
    assert r_pipe["orf.m_dot"] == pytest.approx(r_wall["orf.m_dot"], rel=1e-9)


def _approach_net(dp: float, **orifice_kw) -> FlowNetwork:
    """beta > 0: a supply pipe into a chamber at the orifice's inlet (in_pipe,
    beta = 0.4), discharging into the flowing chamber of _band_net."""
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("gas_in", m_dot=2.0, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("exit", Pt=P_GAS, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cool", Pt=P_GAS * (1.0 + dp), Tt=600.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("ch", area=0.01, main_inlet="feed"))
    g.add_node(MomentumChamberNode("a", area=7.85e-5))
    g.add_element(LosslessConnectionElement("feed", "gas_in", "ch"))
    g.add_element(LosslessConnectionElement("tail", "ch", "exit"))
    g.add_element(ChannelElement("sup", "cool", "a", length=0.05, diameter=0.01))
    g.add_element(OrificeElement("orf", "a", "ch", diameter=4e-3, **orifice_kw))
    return g


@pytest.mark.parametrize(
    ("build", "dp", "sign"),
    [
        (lambda dp: _side_net(dp, tail_channel=False), -0.01, 1.0),  # beta_reverse > 0: band
        (lambda dp: _side_net(dp, tail_channel=False), -0.06, -1.0),  # reverse
        (lambda dp: _side_net(dp, tail_channel=False), 0.02, 1.0),  # forward
        (lambda dp: _side_pipe_net(dp), -0.2, -1.0),  # side, beta_rev 0.4: reverse
        (lambda dp: _side_pipe_net(dp), 0.2, 1.0),
        (_approach_net, -0.01, 1.0),  # beta > 0: band
        (_approach_net, -0.06, -1.0),  # reverse
        (_approach_net, 0.02, 1.0),  # forward
    ],
)
def test_the_jacobian_is_analytic_with_an_approach_pipe(build, dp: float, sign: float) -> None:
    s, r = _solve(build(dp))
    el = s.network.elements["orf"]
    assert el._hole_cd.exit_static is True
    assert el.beta > 0.0 or el.beta_reverse > 0.0
    assert sign * r["orf.m_dot"] > 0.0
    x = np.array(r["__x_solution__"])
    assert _fd_err(s, x) < 1e-5


def test_dp_tap_is_the_reversed_kernels_drive() -> None:
    # Two chambers (the supply pipe's 'a', P < Pt, and the flowing 'ch'),
    # the orifice reversed: the kernel draws from ch's static and discharges
    # against a's Pt, so dP_tap = P_ch - Pt_a (a fixed Cd: no tap rescale).
    # It used to read a's own static, off by a's dynamic head.
    s, r = _solve(_approach_net(-0.06, Cd=0.7, correlation="fixed"))
    assert r["orf.m_dot"] < 0.0
    assert r["a.P"] < r["a.Pt"] < r["ch.P"]
    d = r["__element_diag__"]["orf"]
    assert d["dP_tap"] == pytest.approx(r["ch.P"] - r["a.Pt"], rel=1e-12)
    assert abs(d["dP_tap"] - (r["ch.P"] - r["a.P"])) > 1.0  # not the old reading
    assert d["beta_flow"] == 0.0  # the reversed leg reads no approach


# --- the two-chamber band and the second-law pair ---------------------------------


def _two_chambers(p_a: float) -> FlowNetwork:
    """The reviewer's two-chamber case: supply chamber 'chA' (its main
    inlet a Cd connection from a boundary at p_a, 500 K) jets through a fixed
    Cd connection into the merge chamber 'chB' (2 kg/s at 1400 K, exit 5e5)."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("supA", Pt=p_a, Tt=500.0, Y=Y_AIR))
    g.add_node(MassFlowBoundary("inB", m_dot=2.0, Tt=1400.0, Y=Y_AIR))
    g.add_node(PressureBoundary("exB", Pt=P_GAS, Tt=1400.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("chA", main_inlet="fA"))
    g.add_node(MomentumChamberNode("chB", area=0.01, main_inlet="fB"))
    g.add_element(DiameterDischargeCoefficientConnectionElement("fA", "supA", "chA", 0.02, Cd=0.8))
    g.add_element(LosslessConnectionElement("fB", "inB", "chB"))
    g.add_element(LosslessConnectionElement("tB", "chB", "exB"))
    g.add_element(
        DiameterDischargeCoefficientConnectionElement("orf", "chA", "chB", 6e-3, zeta=1.2)
    )
    return g


@pytest.mark.parametrize(("p_a", "sign"), [(5.1e5, 1.0), (P_GAS * 0.995, 1.0), (4.7e5, -1.0)])
def test_the_two_chamber_side_jet(p_a: float, sign: float) -> None:
    # Forward, inside the band (P_B < Pt_A < Pt_B: it blows), and reversed.
    s, r = _solve(_two_chambers(p_a))
    el = s.network.elements["orf"]
    assert el._hole_cd.exit_static is True
    assert s.network.elements["fA"]._hole_cd.exit_static is False  # chA's main inlet
    assert sign * r["orf.m_dot"] > 0.0
    if p_a == P_GAS * 0.995:
        assert r["chB.P"] < r["chA.Pt"] < r["chB.Pt"]
    x = np.array(r["__x_solution__"])
    assert _fd_err(s, x) < 1e-5


@pytest.mark.parametrize(
    "build",
    [
        lambda: _two_chambers(5.1e5),
        lambda: _two_chambers(4.7e5),
        lambda: _band_net(_orifice(), -0.01),
        lambda: _band_net(_orifice(), -0.05),
        lambda: _approach_net(-0.06, Cd=0.7, correlation="fixed"),
        lambda: _side_pipe_net(-0.2),
    ],
)
def test_the_second_law_pair_is_the_kernels_pair(build) -> None:
    # On every gated leg the control volume runs between the pressures the
    # kernel expands between: forward from Pt_from to the chamber static
    # P_to; reversed from P_to to Pt_from. Read independently from the node
    # states, not from the drive pair.
    from combaero.network.second_law import SecondLawReport

    s, r = _solve(build())
    el = s.network.elements["orf"]
    assert el._hole_cd.exit_static is True
    x = np.array(r["__x_solution__"])
    cv = SecondLawReport(s, x, r)._two_port("orf", el)
    (inflow,), (outflow,) = cv.inflows, cv.outflows
    pt_from = float(r[f"{el.from_node}.Pt"]) if f"{el.from_node}.Pt" in r else None
    if pt_from is None:  # a PressureBoundary supply: P = Pt
        pt_from = float(s.network.nodes[el.from_node].Pt)
    p_to = float(r[f"{el.to_node}.P"])
    if r["orf.m_dot"] >= 0.0:
        assert (inflow.Pt, outflow.Pt) == (pt_from, p_to)
    else:
        assert (inflow.Pt, outflow.Pt) == (p_to, pt_from)
    assert r["__second_law__"]["orf"]["verdict"] == "pass"
