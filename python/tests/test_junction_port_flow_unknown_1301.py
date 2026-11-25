"""A junction port must be fed by an element with its own m_dot (#1301).

A MultiPortChamberBase reads each port's throughflow from the outer element's
``<id>.m_dot`` unknown. A legacy TeeJunctionElement leg (m_dot_com /
m_dot_branch) or another junction (P_jct) has none, and the port read no
flow: on main the issue's network converged at |F| 3.6e-5 while the ledger
reported mass -3.69e-1 and energy -2.26e-1, and a reversed outlet handed its
port node a leaving stream. Such a port is now refused at build time, with
the fix named: a LosslessConnectionElement through one more chamber node,
which gives the port a flow unknown of its own.
"""

import math
import warnings

import pytest

from combaero.network import FlowNetwork, NetworkSolver
from combaero.network import components as cbn
from combaero.network.ideal_tee_element import IdealMomentumTeeElement

_F = 2e-3


def _tee_into_junction(via_lossless: bool) -> FlowNetwork:
    """#1301's network: a branching legacy tee whose straight leg feeds an
    ideal tee's inlet port, directly or through a lossless connection."""
    g = FlowNetwork()
    nodes = ["pc", "ps", "pb", "o1", "o2"] + (["pt"] if via_lossless else [])
    for nid in nodes:
        g.add_node(cbn.MomentumChamberNode(nid, area=_F))
    g.add_node(cbn.PressureBoundary("c", Pt=3e5, Tt=400.0))
    g.add_node(cbn.PressureBoundary("b", Pt=2.9e5, Tt=300.0))
    g.add_node(cbn.PressureBoundary("B1", Pt=2.8e5, Tt=300.0))
    g.add_node(cbn.PressureBoundary("B2", Pt=2.85e5, Tt=300.0))
    g.add_element(cbn.LosslessConnectionElement("lc", "c", "pc"))
    g.add_element(cbn.LosslessConnectionElement("lb", "pb", "b"))
    g.add_element(
        cbn.TeeJunctionElement(
            id="tee",
            common_node="pc",
            straight_node="pt" if via_lossless else "ps",
            branch_node="pb",
            theta=math.radians(90.0),
            F_C=_F,
            psi=1.0,
            tee_type="branching",
        )
    )
    if via_lossless:
        g.add_element(cbn.LosslessConnectionElement("lt", "pt", "ps"))
    g.add_element(
        IdealMomentumTeeElement(
            id="jct", inlet_nodes=["ps"], outlet_nodes=["o1", "o2"], port_areas=[_F] * 3
        )
    )
    d = math.sqrt(4.0 * _F / math.pi) * 0.6
    g.add_element(cbn.OrificeElement("e1", "o1", "B1", Cd=0.8, diameter=d, correlation="fixed"))
    g.add_element(cbn.OrificeElement("e2", "o2", "B2", Cd=0.8, diameter=d, correlation="fixed"))
    return g


def test_a_legacy_tee_leg_into_a_junction_port_is_refused_with_the_fix() -> None:
    with pytest.raises(
        ValueError, match=r"port 'ps' is fed by 'tee' \(TeeJunctionElement\)"
    ) as exc:
        NetworkSolver(_tee_into_junction(via_lossless=False)).solve()
    assert "LosslessConnectionElement" in str(exc.value)


def test_a_junction_port_fed_by_another_junction_is_refused() -> None:
    g = FlowNetwork()
    for nid in ("a", "m", "o1", "o2", "o3"):
        g.add_node(cbn.MomentumChamberNode(nid, area=_F))
    g.add_node(cbn.PressureBoundary("c", Pt=3e5, Tt=300.0))
    for k, pt in (("B1", 2.8e5), ("B2", 2.8e5), ("B3", 2.8e5)):
        g.add_node(cbn.PressureBoundary(k, Pt=pt, Tt=300.0))
    g.add_element(cbn.LosslessConnectionElement("l0", "c", "a"))
    g.add_element(
        IdealMomentumTeeElement(
            id="j1", inlet_nodes=["a"], outlet_nodes=["m", "o1"], port_areas=[_F] * 3
        )
    )
    g.add_element(
        IdealMomentumTeeElement(
            id="j2", inlet_nodes=["m"], outlet_nodes=["o2", "o3"], port_areas=[_F] * 3
        )
    )
    for e, n, b in (("e1", "o1", "B1"), ("e2", "o2", "B2"), ("e3", "o3", "B3")):
        g.add_element(cbn.OrificeElement(e, n, b, Cd=0.8, diameter=0.02, correlation="fixed"))
    with pytest.raises(ValueError, match=r"has no flow unknown of its own"):
        NetworkSolver(g).solve()


def test_the_named_fix_solves_with_a_closed_mass_balance() -> None:
    """The lossless connection the message names gives the port its own
    flow; the network solves and every node's mass closes."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol = NetworkSolver(_tee_into_junction(via_lossless=True)).solve()
    assert sol["__success__"]
    com, branch, straight = sol["tee.m_dot_com"], sol["tee.m_dot_branch"], sol["lt.m_dot"]
    assert straight > 0.0
    assert abs(com - branch - straight) <= 1e-9 * com
    assert abs(straight - sol["e1.m_dot"] - sol["e2.m_dot"]) <= 1e-9 * straight
