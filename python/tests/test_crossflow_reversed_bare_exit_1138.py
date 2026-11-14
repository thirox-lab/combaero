"""A bare crossflow segment's exit follows ChannelElement's (#1138).

A ``CrossflowSegmentElement`` with no station at either end and no
``entry_K`` has no convention of its own: ``_duct_ends`` runs it Pt to Pt, as
a ``ChannelElement``. Its ``_exit_head_lost`` opted out with every other
segment (#588, whose reason, station nodes at static continuity, holds only
for a segment with a station or an entry), so it never paid the exit dynamic
head into a PressureBoundary, in either direction, in either regime: it
flowed 2.9x a ChannelElement between the same nodes (3.72e-2 against
1.28e-2 kg/s at 2 kPa) and the regime gap was +5.6 % (+8.0 % boundary to
boundary). It now loses the head where a ChannelElement does: into a
PressureBoundary, forward or reversed. Into a PlenumNode or a
MomentumChamberNode a ChannelElement recovers the head, and so does the bare
segment, as before.

The probe: two boundaries 2 kPa apart. ``direct``: the element between them;
``first``: the element from ``sup`` into a mid node, a channel on to ``exit``;
``last``: a channel into the mid node, the element on to ``exit``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian, scaled_error

import combaero as cb
from combaero.network import (
    ChannelElement,
    CrossflowSegmentElement,
    FlowNetwork,
    MomentumChamberNode,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)

Y = list(cb.mole_to_mass(cb.species.dry_air()))
A_DUCT = 2.0e-4
D_DUCT = math.sqrt(4.0 * A_DUCT / math.pi)
P0 = 1.0e5
DP = 2000.0


def _element(seg: bool, regime: str, a: str, b: str) -> ChannelElement:
    if seg:
        return CrossflowSegmentElement("s1x", a, b, length=0.05, area=A_DUCT, regime=regime)
    # The segment's roughness (0): the channel's default is 1e-5 m.
    return ChannelElement("s1x", a, b, length=0.05, diameter=D_DUCT, roughness=0.0, regime=regime)


def _net(layout: str, regime: str, reversed_: bool, seg: bool = True) -> FlowNetwork:
    p_sup, p_exit = (P0, P0 + DP) if reversed_ else (P0 + DP, P0)
    g = FlowNetwork()
    g.add_node(PressureBoundary("sup", Pt=p_sup, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("exit", Pt=p_exit, Tt=300.0, Y=Y))
    if layout == "direct":
        g.add_element(_element(seg, regime, "sup", "exit"))
        return g
    mid, pos = layout.split("_")
    g.add_node(
        MomentumChamberNode("mc", area=A_DUCT, regime=regime)
        if mid == "chamber"
        else PlenumNode("mc")
    )
    tail = {"length": 0.05, "diameter": D_DUCT, "roughness": 0.0, "regime": regime}
    if pos == "first":
        g.add_element(_element(seg, regime, "sup", "mc"))
        g.add_element(ChannelElement("tail", "mc", "exit", **tail))
    else:
        g.add_element(ChannelElement("tail", "sup", "mc", **tail))
        g.add_element(_element(seg, regime, "mc", "exit"))
    return g


def _flow(g: FlowNetwork) -> float:
    r = NetworkSolver(g).solve()
    assert r["__success__"], r.get("__message__")
    assert not r["__second_law_failures__"], r["__second_law__"]
    return float(r["s1x.m_dot"])


# (layout, reversed): where the element discharges.
INTO_BOUNDARY = [
    ("plenum_last", False),  # forward into the exit boundary
    ("chamber_last", False),
    ("plenum_first", True),  # reversed into its from-side boundary
    ("chamber_first", True),
    ("direct", False),  # boundary to boundary, both signs
    ("direct", True),
]
INTO_MID_NODE = [
    ("plenum_first", False),
    ("chamber_first", False),
    ("plenum_last", True),
    ("chamber_last", True),
]
REGIMES = ["incompressible", "compressible"]


@pytest.mark.parametrize(("layout", "reversed_"), INTO_BOUNDARY + INTO_MID_NODE)
@pytest.mark.parametrize("regime", REGIMES)
def test_bare_segment_flows_as_a_channel(layout: str, reversed_: bool, regime: str) -> None:
    """Every discharge, both regimes: the bare segment solves the channel's
    flow. On main the boundary cases ran 2.9x it (3.6x boundary to
    boundary); the mid-node ones already matched."""
    ms = _flow(_net(layout, regime, reversed_))
    mc = _flow(_net(layout, regime, reversed_, seg=False))
    assert (ms < 0.0) == reversed_
    assert abs(ms - mc) <= 1e-6 * abs(mc), (ms, mc)


@pytest.mark.parametrize(("layout", "reversed_"), INTO_BOUNDARY)
def test_bare_segment_regime_gap_is_order_mach_squared(layout: str, reversed_: bool) -> None:
    """The regime gap into a boundary was +5.6 % (+8.0 % boundary to
    boundary), the head missing in both regimes but not equally; it is now
    the channel's, below 0.2 % at exit Mach ~0.16."""
    inc = _flow(_net(layout, "incompressible", reversed_))
    com = _flow(_net(layout, "compressible", reversed_))
    assert abs((inc - com) / com) < 3e-3, (inc, com)


@pytest.mark.parametrize("reversed_", [False, True])
def test_bare_segment_loses_the_head_into_a_boundary(reversed_: bool) -> None:
    """Pins the head itself, not only the channel equivalence: boundary to
    boundary, incompressible, the 2 kPa drive is friction plus one exit head,
    (f L / D + 1) m^2 / (2 rho A^2). Friction alone is ~6 % of it (on main
    it was all of it, at 3.6x the flow). The reported rho_ref is the inlet's,
    not the residual's own density: 2 % slack."""
    r = NetworkSolver(_net("direct", "incompressible", reversed_)).solve()
    assert r["__success__"], r.get("__message__")
    m = float(r["s1x.m_dot"])
    q = m * m / (2.0 * float(r["s1x.rho_ref"]) * A_DUCT**2)
    drive = (float(r["s1x.f"]) * 0.05 / D_DUCT + 1.0) * q
    assert abs(drive - DP) < 0.025 * DP, (m, drive)


def _fresh_jacobian_error(g: FlowNetwork, x: np.ndarray) -> float:
    # A fresh solver: the residual depends on solver state after a solve (#1142).
    s = NetworkSolver(g)
    s._build_x0()
    J = s._residuals_and_jacobian(x)[1]
    fd = fd_jacobian(s, x, floor=1e-10)
    rows = [i for i, (owner, _kind) in enumerate(s._row_kinds(x)) if "'s1x'" in owner]
    return float(np.max(scaled_error(J, fd, x)[rows]))


@pytest.mark.parametrize("point", ["root", "perturbed", "zero_flow", "plus", "minus"])
@pytest.mark.parametrize(("layout", "reversed_"), [c for c in INTO_BOUNDARY if c[0] != "direct"])
@pytest.mark.parametrize("regime", REGIMES)
def test_bare_segment_jacobian(layout: str, reversed_: bool, regime: str, point: str) -> None:
    """The segment row's analytic Jacobian against central differences in
    scaled variables, at the root, a perturbation and across m = 0 (the head
    m|m| / (2 rho A^2) has a slope |m| / (rho A^2): C1 at m = 0). Measured
    below 1e-8 everywhere. The boundary-to-boundary network has the flow as
    its only unknown: at m = 0 exactly its incompressible column is the
    channel's laminar friction slope, which the analytic form reads as 0 on
    main and for a ChannelElement alike, so it is left out."""
    s = NetworkSolver(_net(layout, regime, reversed_))
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    x = np.array(r["__x_solution__"], dtype=float)
    im = s.unknown_names.index("s1x.m_dot")
    if point == "perturbed":
        x = x * (1.0 + 1e-3 * np.random.default_rng(7).uniform(-1.0, 1.0, x.size))
    elif point != "root":
        x[im] = {"zero_flow": 0.0, "plus": 1e-7, "minus": -1e-7}[point]
    assert _fresh_jacobian_error(_net(layout, regime, reversed_), x) < 1e-7


@pytest.mark.parametrize("regime", REGIMES)
def test_stationed_segment_still_opts_out(regime: str) -> None:
    """A segment with a station or an entry keeps #588's opt-out: its exit
    already discharges at static continuity, and ChannelElement's head would
    charge it twice."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("sup", Pt=P0 + DP, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("exit", Pt=P0, Tt=300.0, Y=Y))
    g.add_element(
        CrossflowSegmentElement(
            "s1x", "sup", "exit", length=0.05, area=A_DUCT, entry_K=0.5, regime=regime
        )
    )
    g.resolve_all_topology()
    seg = g.elements["s1x"]
    # Not vacuous: ChannelElement's rule would charge it, both ways.
    assert ChannelElement._exit_head_lost(seg, 1.0)
    assert ChannelElement._exit_head_lost(seg, -1.0)
    assert not seg._exit_head_lost(1.0)
    assert not seg._exit_head_lost(-1.0)


@pytest.mark.parametrize("regime", REGIMES)
def test_a_segment_leaving_a_to_only_station_still_opts_out(regime: str) -> None:
    """A segment with no ``from_kappa`` leaving a station the previous
    segment registers (its ``to_kappa``, B1) is not bare: it is fed at the
    station's static pressure, the gas already moving, so charging its exit
    head into the boundary would make it pay a head it never had to make.
    Its residual is main's (bit-identical, #1138's hash check)."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("sup", Pt=P0 + DP, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("exit", Pt=P0, Tt=300.0, Y=Y))
    g.add_node(PlenumNode("n1"))
    g.add_element(
        CrossflowSegmentElement(
            "s0",
            "sup",
            "n1",
            length=0.05,
            area=A_DUCT,
            entry_K=0.0,
            to_kappa=cb.STATION_KAPPA_BLEED_BASSETT,
            next_seg="s1x",
            regime=regime,
        )
    )
    g.add_element(
        CrossflowSegmentElement("s1x", "n1", "exit", length=0.05, area=A_DUCT, regime=regime)
    )
    g.resolve_all_topology()
    seg = g.elements["s1x"]
    assert g.nodes["n1"]._crossflow_station.segment_ids == ("s0", "s1x")
    assert ChannelElement._exit_head_lost(seg, 1.0)
    assert not seg._exit_head_lost(1.0)
    assert NetworkSolver(g).solve()["__success__"]
