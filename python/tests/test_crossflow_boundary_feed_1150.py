"""A stationed crossflow segment fed straight from a node at rest (#1150).

A ``CrossflowSegmentElement`` end fed from a node at rest (a plenum, a
combustor, a boundary: ``state_at_rest``) with no ``entry_K`` and no station
of its own is a lossless reservoir entry, the same as ``entry_K=0`` (owner
decision 2026-11-14). It applies to every segment with a station, its own or
a neighbour's, under one "has a station" predicate (``_is_bare``) that every
gate reads. A bare segment is a ChannelElement, its feed already in Pt.

Before, the incompressible row had no term for such a feed, and the gates
disagreed on what a station is: the reservoir ends and the duct ends asked
only for the segment's own ``from_kappa`` / ``entry_K``. Two boundaries 2 kPa
apart, two 0.05 m segments of 2e-4 m^2:

* self-stationed (x has its own to-station, no ``entry_K``): 3.8798e-2 kg/s
  forward against -1.2840e-2 reversed incompressible; compressible
  1.2824e-2 both ways (3.0x regime gap);
* B1 (x leaves the to-only station s0 registers) reversed and B2 (x arrives
  at s1's from-station) forward: 3.8798e-2 incompressible, 3.7111e-2
  compressible, against 1.284e-2 the other way.

Now each layout flows 1.28403e-2 incompressible and 1.28241e-2 compressible
in both directions (regime gap +0.126 %). A segment with no ``entry_K`` flows
the one with ``entry_K=0`` bit for bit, and a segment whose only station is a
neighbour's flows the one whose station is its own (to 1e-9).
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
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

Y = list(cb.mole_to_mass(cb.species.dry_air()))
A = 2.0e-4
D = math.sqrt(4.0 * A / math.pi)
P0 = 1.0e5
DP = 2000.0
KAPPA = cb.STATION_KAPPA_BLEED_BASSETT
REGIMES = ["incompressible", "compressible"]


def _seg(eid: str, a: str, b: str, regime: str, **kw: object) -> CrossflowSegmentElement:
    return CrossflowSegmentElement(eid, a, b, length=0.05, area=A, regime=regime, **kw)


def _ends(reversed_: bool) -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(PressureBoundary("sup", Pt=P0 if reversed_ else P0 + DP, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("exit", Pt=P0 + DP if reversed_ else P0, Tt=300.0, Y=Y))
    g.add_node(PlenumNode("n1"))
    return g


def _selfst(
    regime: str, reversed_: bool, entry_K: float | None = None, feed: str = "boundary",
    bleed: bool = False,
) -> FlowNetwork:  # fmt: skip
    """x src -> n1 with its own to-station (next s1); s1 n1 -> exit owns the
    from-station. ``feed``: x straight from ``sup``, or from a node ``p0``
    (a PlenumNode or a MomentumChamberNode) that a channel feeds."""
    g = _ends(reversed_)
    src = "sup"
    if feed != "boundary":
        g.add_node(
            PlenumNode("p0")
            if feed == "plenum"
            else MomentumChamberNode("p0", area=A, regime=regime)
        )
        g.add_element(
            ChannelElement(
                "feed", "sup", "p0", length=0.05, diameter=D, roughness=0.0, regime=regime
            )
        )
        src = "p0"
    kw = {} if entry_K is None else {"entry_K": entry_K}
    g.add_element(_seg("x", src, "n1", regime, to_kappa=KAPPA, next_seg="s1", **kw))
    g.add_element(_seg("s1", "n1", "exit", regime, from_kappa=KAPPA, prev_seg="x"))
    if bleed:
        g.add_node(PressureBoundary("sink", Pt=P0, Tt=300.0, Y=Y))
        g.add_element(
            OrificeElement("bl", "n1", "sink", diameter=0.006, Cd=0.6, correlation="fixed")
        )
    return g


def _b1(
    regime: str, reversed_: bool, own: bool = False, chamber: bool = False,
    area_to: float | None = None,
) -> FlowNetwork:  # fmt: skip
    """s0 (entry_K = 0, to_kappa, next x) sup -> n1; x leaves n1 with no
    station of its own (``own``: its own from_kappa), into ``exit`` or
    (``chamber``) into a MomentumChamberNode a channel drains; ``area_to``
    tapers x."""
    g = _ends(reversed_)
    g.add_element(_seg("s0", "sup", "n1", regime, entry_K=0.0, to_kappa=KAPPA, next_seg="x"))
    kw: dict = {"from_kappa": KAPPA, "prev_seg": "s0"} if own else {}
    if area_to is not None:
        kw["area_to"] = area_to
    if chamber:
        g.add_node(MomentumChamberNode("mc", area=A, regime=regime))
        g.add_element(_seg("x", "n1", "mc", regime, **kw))
        g.add_element(
            ChannelElement(
                "tail", "mc", "exit", length=0.05, diameter=D, roughness=0.0, regime=regime
            )
        )
    else:
        g.add_element(_seg("x", "n1", "exit", regime, **kw))
    return g


def _b2(regime: str, reversed_: bool, own: bool = False) -> FlowNetwork:
    """x sup -> n1 with no station of its own (``own``: its own to-station);
    s1 n1 -> exit owns the from-station at n1, x arriving."""
    g = _ends(reversed_)
    kw = {"to_kappa": KAPPA, "next_seg": "s1"} if own else {}
    g.add_element(_seg("x", "sup", "n1", regime, **kw))
    g.add_element(_seg("s1", "n1", "exit", regime, from_kappa=KAPPA, prev_seg="x"))
    return g


def _solve(g: FlowNetwork) -> tuple[dict, NetworkSolver]:
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    assert not r["__second_law_failures__"], r["__second_law__"]
    return r, s


def _x(g: FlowNetwork) -> float:
    return float(_solve(g)[0]["x.m_dot"])


LAYOUTS = {
    "selfst": _selfst,
    "b1": _b1,
    "b2": _b2,
    "b1_chamber": lambda r, v: _b1(r, v, chamber=True),
}


@pytest.mark.parametrize("layout", list(LAYOUTS))
@pytest.mark.parametrize("regime", REGIMES)
def test_boundary_fed_segment_is_direction_symmetric(layout: str, regime: str) -> None:
    """Each layout is its own mirror image. Before, the boundary-fed leg ran
    3.0x the other incompressible (3.8798e-2 against 1.2840e-2), and 2.9x
    compressible in B1 / B2 (3.7111e-2); B1 into a chamber ran 9.20e-3
    forward against -3.24e-2 reversed. Now 4e-15 incompressible, 5e-9
    compressible; into the chamber 1.0e-6 incompressible, the same as with
    x's own station on main (the chamber's own closure)."""
    fwd = _x(LAYOUTS[layout](regime, False))
    rev = _x(LAYOUTS[layout](regime, True))
    assert abs(abs(rev) / fwd - 1.0) < 1e-5, (fwd, rev)


@pytest.mark.parametrize("layout", list(LAYOUTS))
@pytest.mark.parametrize("reversed_", [False, True])
def test_boundary_fed_segment_regime_gap_is_order_mach_squared(
    layout: str, reversed_: bool
) -> None:
    """Regime gap at exit Mach ~0.16: +0.126 % (+0.191 % into the chamber).
    The self-stationed forward leg was +203 % (3.8798e-2 against 1.2824e-2)."""
    inc = _x(LAYOUTS[layout]("incompressible", reversed_))
    com = _x(LAYOUTS[layout]("compressible", reversed_))
    assert abs((inc - com) / com) < 0.003, (inc, com)


def _FJ(g: FlowNetwork, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    g.resolve_all_topology()
    s = NetworkSolver(g)
    s._build_x0()
    F, J = s._residuals_and_jacobian(x)
    J = J.toarray() if hasattr(J, "toarray") else np.asarray(J, dtype=float)
    return np.asarray(F, dtype=float), np.asarray(J, dtype=float)


@pytest.mark.parametrize(
    ("build", "kw"),
    [(_selfst, {}), (_selfst, {"feed": "plenum"}), (_selfst, {"bleed": True})],
    ids=["boundary", "plenum", "bleed"],
)
@pytest.mark.parametrize("reversed_", [False, True])
@pytest.mark.parametrize("regime", REGIMES)
def test_no_entry_K_is_entry_K_zero_bit_for_bit(
    build: object, kw: dict, reversed_: bool, regime: str
) -> None:
    """The owner's rule itself: no entry_K from a node at rest IS entry_K=0.
    Residuals and Jacobian on fresh, resolved solvers, at the root and a
    perturbation, bit for bit in both regimes (compressible the entry adds
    K/(1+K) = 0 of its drop, so the term is left out)."""
    r, _ = _solve(build(regime, reversed_, 0.0, **kw))
    x0 = np.array(r["__x_solution__"], dtype=float)
    rng = np.random.default_rng(5)
    for x in (x0, x0 * (1.0 + 1e-3 * rng.standard_normal(x0.size))):
        F0, J0 = _FJ(build(regime, reversed_, None, **kw), x)
        F1, J1 = _FJ(build(regime, reversed_, 0.0, **kw), x)
        assert F0.tobytes() == F1.tobytes()
        assert J0.tobytes() == J1.tobytes()


NEIGHBOUR_VS_OWN = {
    "b1": (lambda r, v: _b1(r, v), lambda r, v: _b1(r, v, own=True)),
    "b2": (lambda r, v: _b2(r, v), lambda r, v: _b2(r, v, own=True)),
    "b1_chamber": (
        lambda r, v: _b1(r, v, chamber=True),
        lambda r, v: _b1(r, v, own=True, chamber=True),
    ),
    # x tapers to 0.6 A: the area change between the station's static and
    # the exit's is Bernoulli's, at a neighbour's station as at its own
    # (B1 taper was 1.258e-2 / -3.285e-2 incompressible, now 7.959e-3 /
    # -1.2408e-2, the own station's).
    "b1_taper": (
        lambda r, v: _b1(r, v, area_to=0.6 * A),
        lambda r, v: _b1(r, v, own=True, area_to=0.6 * A),
    ),
}


@pytest.mark.parametrize("layout", list(NEIGHBOUR_VS_OWN))
@pytest.mark.parametrize("reversed_", [False, True])
@pytest.mark.parametrize("regime", REGIMES)
def test_a_neighbours_station_is_a_station(layout: str, reversed_: bool, regime: str) -> None:
    """A segment whose only station is a neighbour's flows the one whose
    station is its own (no side stream, so no half-drop either way): to
    4e-15 incompressible and 1.4e-9 compressible. Its end at that station
    takes the node's P as the duct static, closed to its moving Pt. With
    only the gate flipped, that end still read the node's P as a reservoir
    Pt and B1 forward charged the head twice: 9.34e-3 against 1.28e-2
    (-27 %)."""
    neighbour, own = NEIGHBOUR_VS_OWN[layout]
    m = _x(neighbour(regime, reversed_))
    m_own = _x(own(regime, reversed_))
    assert abs(m - m_own) <= 1e-8 * abs(m_own), (m, m_own)


# Main's flows (ad5b4159) for a segment that sets entry_K = 0.5: unchanged.
ENTRY_HALF = {
    ("incompressible", False): 0.010671447204308378,
    ("incompressible", True): -0.012840326341502627,
    ("compressible", False): 0.010670618594559488,
    ("compressible", True): -0.012824138575944104,
}


@pytest.mark.parametrize("reversed_", [False, True])
@pytest.mark.parametrize("regime", REGIMES)
def test_a_set_entry_K_is_unchanged(regime: str, reversed_: bool) -> None:
    """A segment that sets entry_K keeps its own entry: no second head on
    top (its residuals are main's bit for bit, the PR's hash check)."""
    m = _x(_selfst(regime, reversed_, 0.5))
    assert abs(m - ENTRY_HALF[(regime, reversed_)]) <= 1e-9 * abs(m), m


def test_the_one_predicate() -> None:
    """Which ends are terminal and which feed is a lossless entry."""
    for g, eid, feed in (
        (_selfst("incompressible", False), "x", 0.0),
        (_selfst("incompressible", False, 0.5), "x", 0.5),
        (_selfst("incompressible", False, feed="plenum"), "x", 0.0),
        # A moving node is not a node at rest: outside the rule (#1150).
        (_selfst("incompressible", False, feed="chamber"), "x", None),
        (_b2("incompressible", False), "x", 0.0),
        # x's from-end is the station s0 registers: no feed there.
        (_b1("incompressible", False), "x", None),
    ):
        g.resolve_all_topology()
        assert g.elements[eid]._feed_K(True) == feed, (eid, feed)
    g = _b1("incompressible", False)
    g.resolve_all_topology()
    x = g.elements["x"]
    assert not x._is_bare()
    assert x._station_end("from") and not x._station_end("to")
    assert x._terminal_end("to") and not x._terminal_end("from")
    # Compressible: the duct end already takes the node's Pt, no term.
    g = _selfst("compressible", False)
    g.resolve_all_topology()
    assert g.elements["x"]._feed_K(False) is None
    # A bare segment has no terminal ends: a ChannelElement.
    g = FlowNetwork()
    g.add_node(PressureBoundary("sup", Pt=P0 + DP, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("exit", Pt=P0, Tt=300.0, Y=Y))
    g.add_element(_seg("x", "sup", "exit", "incompressible"))
    g.resolve_all_topology()
    x = g.elements["x"]
    assert x._is_bare() and not x._terminal_end("from") and x._feed_K(True) is None


POINTS = ["root", "perturbed", "m=0", "m=+1e-7", "m=-1e-7"]


@pytest.mark.parametrize("point", POINTS)
@pytest.mark.parametrize("layout", list(LAYOUTS))
@pytest.mark.parametrize("reversed_", [False, True])
@pytest.mark.parametrize("regime", REGIMES)
def test_jacobian(layout: str, reversed_: bool, regime: str, point: str) -> None:
    """Analytic Jacobian against central differences in scaled variables on
    a fresh solver with its topology resolved, at the root, a perturbation,
    and with every segment flow at 0 and +-1e-7 (the entry head is m|m|-like
    one side and 0 the other: C1 at m = 0). Measured: 1.4e-8 at worst."""
    build = LAYOUTS[layout]
    r, s = _solve(build(regime, reversed_))
    x = np.array(r["__x_solution__"], dtype=float)
    names = list(s.unknown_names)
    segs = [
        i for i, n in enumerate(names)
        if n.endswith(".m_dot") and isinstance(s.network.elements.get(n[:-6]), CrossflowSegmentElement)
    ]  # fmt: skip
    if point == "perturbed":
        x = x * (1.0 + 1e-3 * np.random.default_rng(7).uniform(-1.0, 1.0, x.size))
    elif point != "root":
        x[segs] = float(point.split("=")[1])
    g = build(regime, reversed_)
    g.resolve_all_topology()
    f = NetworkSolver(g)
    f._build_x0()
    J = f._residuals_and_jacobian(x)[1]
    J = J.toarray() if hasattr(J, "toarray") else np.asarray(J, dtype=float)
    err = float(np.max(scaled_error(J, fd_jacobian(f, x, floor=1e-10), x)))
    assert err < 1e-7, err


def _bleed(g: FlowNetwork, node: str) -> FlowNetwork:
    """A bleed orifice from ``node`` to a sink 1 kPa below ``P0``: a side
    stream at the station."""
    if "sink" not in g.nodes:
        g.add_node(PressureBoundary("sink", Pt=P0 - 1000.0, Tt=300.0, Y=Y))
    g.add_element(
        OrificeElement(f"bl_{node}", node, "sink", diameter=0.006, Cd=0.6, correlation="fixed")
    )
    return g


def _interior(regime: str, reversed_: bool) -> FlowNetwork:
    """s0 (entry_K = 0, to_kappa, next x) sup -> n1; x n1 -> n2 with no
    station of its own, both ends its neighbours'; s2 n2 -> exit owns the
    from-station at n2, x arriving. A bleed at each station."""
    g = _ends(reversed_)
    g.add_node(PlenumNode("n2"))
    g.add_element(_seg("s0", "sup", "n1", regime, entry_K=0.0, to_kappa=KAPPA, next_seg="x"))
    g.add_element(_seg("x", "n1", "n2", regime))
    g.add_element(_seg("s2", "n2", "exit", regime, from_kappa=KAPPA, prev_seg="x"))
    return _bleed(_bleed(g, "n1"), "n2")


SIDE_STREAM = {
    "b1_bleed": (lambda r, v: _bleed(_b1(r, v), "n1"), "s0+x"),
    "b2_bleed": (lambda r, v: _bleed(_b2(r, v), "n1"), "s1+x"),
    "interior_bleed": (_interior, "s0+s2+x"),
}


@pytest.mark.parametrize("layout", list(SIDE_STREAM))
@pytest.mark.parametrize("reversed_", [False, True])
@pytest.mark.parametrize("regime", REGIMES)
def test_a_side_stream_at_a_neighbours_station_passes_the_second_law(
    layout: str, reversed_: bool, regime: str
) -> None:
    """A segment that names no neighbour (its only station is the one its
    neighbour registers) is in that neighbour's chain: the second law's
    control volume spans both. Before the chains were walked both ways, x
    was a chain of its own whenever it came first, its control volume
    booked the station's static P as a Pt, and B2 with a bleed at n1
    forward was refused (x = 1.331e-2 kg/s, s1 sigma_norm -2.6e-4); the
    interior chain with bleeds was refused on main too (s2 -2.9e-4)."""
    build, group = SIDE_STREAM[layout]
    r, _ = _solve(build(regime, reversed_))
    sl = r["__second_law__"]
    assert sl["x"]["verdict"] == "pass", sl["x"]
    assert sl["x"].get("group") == group, sl["x"]


def _hot(regime: str, entry_K: float | None = None) -> FlowNetwork:
    """The self-stationed x fed from a 900 K boundary; a 300 K side stream
    joins at n1, so the station's gas (~780 K) is not the reservoir's."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("sup", Pt=P0 + DP, Tt=900.0, Y=Y))
    g.add_node(PressureBoundary("exit", Pt=P0, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("cold", Pt=P0 + 2.0 * DP, Tt=300.0, Y=Y))
    g.add_node(PlenumNode("n1"))
    kw = {} if entry_K is None else {"entry_K": entry_K}
    g.add_element(_seg("x", "sup", "n1", regime, to_kappa=KAPPA, next_seg="s1", **kw))
    g.add_element(_seg("s1", "n1", "exit", regime, from_kappa=KAPPA, prev_seg="x"))
    g.add_element(OrificeElement("c", "cold", "n1", diameter=0.006, Cd=0.6, correlation="fixed"))
    return g


def test_the_lossless_entry_accelerates_the_reservoirs_gas() -> None:
    """The implicit entry is entry_K = 0's, so its head is at the density
    of the gas it carries, the reservoir's (#1137), not the station's mixed
    gas: regime gap -0.075 % (5.68978e-3 incompressible, 5.69403e-3
    compressible). At the station's gas it was +4.6 %."""
    inc = _x(_hot("incompressible"))
    com = _x(_hot("compressible"))
    assert abs((inc - com) / com) < 0.005, (inc, com)
