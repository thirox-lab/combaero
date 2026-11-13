"""Crossflow station nodes must be PlenumNodes.

A station of a crossflow chain (a segment's from-station, a to-station with a
next segment, the node an impingement row discharges into) holds the
channel's static pressure, the pressure a side stream there sees. A
MomentumChamberNode there modelled the same stream twice, lumped and
distributed, and crashed inside ``_station_jets`` (the chamber's inflows are
triples). It is refused at validate, naming the rule. A chain's exit segment
(no to-station) into a chamber is a terminal node, not a station: it solves.
"""

import warnings

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    ImpingementCrossflowElement,
    ImpingementPlateElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)

D = 0.00254
SPAN = 12 * 4.0 * D
Y = list(cb.mole_to_mass(cb.species.dry_air()))
KW = {"height": 2.0 * D, "span": SPAN}


def _network(rows: int, bypass: bool, station: str, pt: float = 1.05e5) -> FlowNetwork:
    """A rows-row plate over a crossflow gap; ``station`` picks the node type
    at each row: 'plenum', or a MomentumChamberNode whose main inlet is the
    arriving 'segment', the row's 'plate', or 'none'."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("sup", Pt=pt, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("exit", Pt=1.0e5, Tt=300.0, Y=Y))
    for i in range(1, rows + 1):
        arriving = f"x{i - 1}" if (i > 1 or bypass) else None
        if station == "plenum":
            g.add_node(PlenumNode(f"c{i}"))
        else:
            main = {"segment": arriving, "plate": f"p{i}", "none": None}[station]
            g.add_node(MomentumChamberNode(f"c{i}", main_inlet=main))
    if bypass:
        g.add_node(MassFlowBoundary("feed_in", m_dot=0.01, Tt=300.0, Y=Y))
        g.add_node(PlenumNode("byp"))
        g.add_element(ChannelElement("feed", "feed_in", "byp", length=0.2, diameter=0.02))
        g.add_element(
            ImpingementCrossflowElement("x0", "byp", "c1", length=2.5 * D, bypass_inlet=True, **KW)
        )
    for i in range(1, rows + 1):
        g.add_element(
            ImpingementPlateElement(
                f"p{i}",
                "sup",
                f"c{i}",
                d_jet=D,
                xn_d=5.0,
                yn_d=4.0,
                z_d=2.0,
                span=SPAN,
                plate_thickness=D,
                row=i,
            )
        )
        to = f"c{i + 1}" if i < rows else "exit"
        g.add_element(ImpingementCrossflowElement(f"x{i}", f"c{i}", to, length=5.0 * D, **KW))
    return g


def _solve(g: FlowNetwork) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return NetworkSolver(g).solve()


_REFUSAL = r"must be (a PlenumNode|PlenumNodes).*either lumped \(a MomentumChamberNode\)"


@pytest.mark.parametrize("station", ["segment", "plate", "none"])
@pytest.mark.parametrize("rows", [1, 2])
@pytest.mark.parametrize("bypass", [False, True])
def test_momentum_chamber_station_is_refused(station: str, rows: int, bypass: bool) -> None:
    if station == "segment" and rows == 1 and not bypass:
        pytest.skip("no segment arrives at row 1 without a bypass")
    with pytest.raises(ValueError, match=_REFUSAL):
        _solve(_network(rows, bypass, station))


@pytest.mark.parametrize("rows", [1, 2])
@pytest.mark.parametrize("bypass", [False, True])
def test_plenum_stations_solve(rows: int, bypass: bool) -> None:
    r = _solve(_network(rows, bypass, "plenum"))
    assert r["__success__"], r.get("__message__")
    for i in range(1, rows + 1):
        assert r[f"p{i}.m_dot"] > 0.0


def test_generic_segment_station_is_refused() -> None:
    """The rule sits on CrossflowSegmentElement, not only the impingement
    subclass: a jet station on a MomentumChamberNode is refused too."""
    from combaero.network import CrossflowSegmentElement, OrificeElement

    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=1.05e5, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("jet", Pt=1.06e5, Tt=300.0, Y=Y))
    g.add_node(PlenumNode("a"))
    g.add_node(MomentumChamberNode("s", main_inlet="s1"))
    A = 1.0e-4
    g.add_element(ChannelElement("feed", "in", "a", length=0.05, diameter=0.011))
    g.add_element(
        CrossflowSegmentElement(
            "s1",
            "a",
            "s",
            length=0.05,
            area=A,
            to_kappa=cb.STATION_KAPPA_MERGE_NORMAL,
            next_seg="s2",
        )
    )
    g.add_element(
        CrossflowSegmentElement(
            "s2",
            "s",
            "out",
            length=0.05,
            area=A,
            from_kappa=cb.STATION_KAPPA_MERGE_NORMAL,
            prev_seg="s1",
        )
    )
    g.add_element(OrificeElement("o", "jet", "s", diameter=0.004, Cd=0.6, correlation="fixed"))
    with pytest.raises(ValueError, match=r"CrossflowSegmentElement 's1'.*station node 's'"):
        _solve(g)


@pytest.mark.parametrize("main_inlet", [None, "x2"])
def test_chain_exit_into_momentum_chamber_solves(main_inlet: str | None) -> None:
    """An exit segment (no to-station) into a chamber is a terminal node: the
    duct's exit Pt meets the chamber's Pt (``_duct_ends``), not refused."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("sup", Pt=1.05e5, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("exit", Pt=1.0e5, Tt=300.0, Y=Y))
    for i in (1, 2):
        g.add_node(PlenumNode(f"c{i}"))
        g.add_element(
            ImpingementPlateElement(
                f"p{i}",
                "sup",
                f"c{i}",
                d_jet=D,
                xn_d=5.0,
                yn_d=4.0,
                z_d=2.0,
                span=SPAN,
                plate_thickness=D,
                row=i,
            )
        )
    g.add_node(MomentumChamberNode("mc", main_inlet=main_inlet))
    g.add_element(ImpingementCrossflowElement("x1", "c1", "c2", length=5.0 * D, **KW))
    g.add_element(ImpingementCrossflowElement("x2", "c2", "mc", length=5.0 * D, **KW))
    g.add_element(ChannelElement("tail", "mc", "exit", length=0.05, diameter=0.012))
    r = _solve(g)
    assert r["__success__"], r.get("__message__")
    # The chamber carries the exit's moving state: Pt above P.
    assert r["mc.Pt"] > r["mc.P"]
