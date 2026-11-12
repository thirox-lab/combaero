"""GUI 'impingement_array' node: one configuration, expanded into rows (#465).

The node is the channel-exit configuration only. The backend expands it into
ImpingementArray's row chain, reports the rows back under the node's own id,
and spreads a thermal edge to it over the rows.
"""

from __future__ import annotations

import pytest

from combaero.network import ImpingementCrossflowElement, ImpingementPlateElement, PlenumNode
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import ImpingementArrayData, NetworkGraphSchema, OrificeData

# Over 20 s on one worker: runs in CI's parallel slow job (see
# validation-tests.yml), still before every merge.
pytestmark = pytest.mark.slow


def _node(nid, ntype, data):
    return {"id": nid, "type": ntype, "position": {"x": 0, "y": 0}, "data": data}


def _edge(eid, src, tgt, data=None):
    return {"id": eid, "source": src, "target": tgt, "data": data or {}}


def _schema(with_wall: bool = False, upstream_channel: bool = False) -> dict:
    nodes = [
        _node("supply", "pressure_boundary", {"Pt": 1.03e5, "Tt": 300.0}),
        _node("exit", "pressure_boundary", {"Pt": 1.0e5, "Tt": 300.0}),
        _node("arr", "impingement_array", {}),  # every default: the working example
    ]
    edges = [_edge("e2", "arr", "exit")]
    if upstream_channel:
        nodes.append(_node("feed", "channel", {"L": 0.1, "D": 0.1}))
        edges += [_edge("e0", "supply", "feed"), _edge("e1", "feed", "arr")]
    else:
        edges.append(_edge("e1", "supply", "arr"))
    if with_wall:
        nodes += [
            _node("hot_in", "pressure_boundary", {"Pt": 1.02e5, "Tt": 1200.0}),
            _node("hot_out", "pressure_boundary", {"Pt": 1.0e5, "Tt": 1200.0}),
            _node("hot", "channel", {"L": 0.127, "D": 0.05}),
        ]
        edges += [
            _edge("h1", "hot_in", "hot"),
            _edge("h2", "hot", "hot_out"),
            _edge(
                "wall",
                "hot",
                "arr",
                {"type": "thermal", "thickness": 0.001, "conductivity": 20.0},
            ),
        ]
    return {"nodes": nodes, "edges": edges}


def test_the_node_expands_into_the_row_chain() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_schema()))
    n = ImpingementArrayData().n_rows
    for i in range(1, n + 1):
        assert isinstance(net.elements[f"arr__p{i}"], ImpingementPlateElement)
        assert isinstance(net.elements[f"arr__x{i}"], ImpingementCrossflowElement)
        assert net.elements[f"arr__p{i}"].from_node == "supply"
    assert net.elements[f"arr__x{n}"].to_node == "exit"


def test_the_default_is_a_working_example() -> None:
    """Defaults solve and land inside Florschuetz's range on every row."""
    result = NetworkRunner.from_dict(_schema()).solve(timeout=120.0)
    r = result._element_results["arr"]
    assert r.success
    d = r.model_dump()
    assert d["Re_j_min"] > 2.5e3 and d["Re_j_max"] < 7.0e4
    assert d["Gc_Gj_max"] < 0.8
    assert d["surface_extrapolated"] == 0.0
    assert len(d["rows_Nu"]) == ImpingementArrayData().n_rows


def test_an_upstream_element_feeds_through_a_plenum() -> None:
    """A junction MomentumChamberNode cannot split one stream into N plates."""
    net = build_network_from_schema(NetworkGraphSchema(**_schema(upstream_channel=True)))
    supply = net.elements["arr__p1"].from_node
    assert isinstance(net.nodes[supply], PlenumNode)
    result = NetworkRunner.from_dict(_schema(upstream_channel=True)).solve(timeout=120.0)
    assert result._element_results["arr"].success


def test_one_thermal_edge_heats_every_row_and_reports_as_one() -> None:
    schema = _schema(with_wall=True)
    net = build_network_from_schema(NetworkGraphSchema(**schema))
    n = ImpingementArrayData().n_rows
    assert sorted(net.walls) == sorted(f"wall__r{i}" for i in range(1, n + 1))
    result = NetworkRunner.from_dict(schema).solve(timeout=120.0)
    edge = result._edge_results["wall"]
    assert len(edge["rows_Q"]) == n
    assert edge["Q"] == pytest.approx(sum(edge["rows_Q"]), rel=1e-12)
    assert edge["Q"] > 0.0
    assert edge["T_hot"] == max(edge["rows_T_hot"])
    # The per-row walls are internal ids: kept beside the aggregate, a sum of
    # the export's wall rows counted the array's heat twice.
    assert sorted(result._edge_results) == ["wall"]
    walls = result.to_dataframe().query("type == 'thermal_wall'")
    assert walls["Q [W]"].sum() == pytest.approx(edge["Q"], rel=1e-12)


def test_orifice_offers_lichtarowicz() -> None:
    assert (
        OrificeData(hole_cd={"correlation": "Lichtarowicz"}).hole_cd.correlation == "Lichtarowicz"
    )


def _bypass_schema(feed: str) -> dict:
    """The default array with its "b" handle wired. ``feed``: 'channel' (a
    duct from a pressure boundary), 'boundary' (wired straight in) or
    'array' (an upstream array's exit, the gap continuing)."""
    s = _schema()
    if feed == "array":
        s["nodes"].append(_node("up", "impingement_array", {"n_rows": 4}))
        s["edges"] += [
            _edge("u1", "supply", "up"),
            _edge("u2", "up", "arr"),
        ]
        s["edges"][-1]["targetHandle"] = "port-bypass-target"
        return s
    s["nodes"].append(_node("bsrc", "pressure_boundary", {"Pt": 1.032e5, "Tt": 300.0}))
    if feed == "channel":
        s["nodes"].append(_node("bduct", "channel", {"L": 0.1, "D": 0.02}))
        s["edges"] += [_edge("b0", "bsrc", "bduct"), _edge("b1", "bduct", "arr")]
    else:
        s["edges"].append(_edge("b1", "bsrc", "arr"))
    s["edges"][-1]["targetHandle"] = "port-bypass-target"
    return s


def test_an_unwired_bypass_handle_is_the_closed_end() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_schema()))
    assert "arr__x0" not in net.elements


@pytest.mark.parametrize("feed", ["channel", "boundary", "array"])
def test_the_bypass_handle_adds_the_lead_in(feed: str) -> None:
    """#467: the "b" handle is a second inflow that does not count as the
    supply. The lead-in joins whatever feeds it to row 1's node."""
    net = build_network_from_schema(NetworkGraphSchema(**_bypass_schema(feed)))
    x0 = net.elements["arr__x0"]
    assert isinstance(x0, ImpingementCrossflowElement) and x0.bypass_inlet
    assert x0.to_node == "arr__c1"
    assert net.elements["arr__p1"].from_node == "supply"
    if feed == "boundary":
        assert x0.from_node == "bsrc"
    else:
        assert isinstance(net.nodes[x0.from_node], PlenumNode)
    if feed == "array":
        assert net.elements["up__x4"].to_node == x0.from_node


@pytest.mark.parametrize("feed", ["channel", "array"])
def test_a_wired_bypass_solves_and_reports(feed: str) -> None:
    result = NetworkRunner.from_dict(_bypass_schema(feed)).solve(timeout=120.0)
    d = result._element_results["arr"].model_dump()
    assert result._element_results["arr"].success
    assert d["m_dot_bypass"] > 0.0
    assert d["rows_Gc_Gj"][0] > 0.0  # row 1 now sees a crossflow
    assert d["heat_transfer_outside_source"] == 0.0


def test_two_bypass_links_are_refused() -> None:
    s = _bypass_schema("boundary")
    s["nodes"].append(_node("bsrc2", "pressure_boundary", {"Pt": 1.032e5, "Tt": 300.0}))
    s["edges"].append(_edge("b2", "bsrc2", "arr"))
    s["edges"][-1]["targetHandle"] = "port-bypass-target"
    with pytest.raises(ValueError, match="two bypass links"):
        build_network_from_schema(NetworkGraphSchema(**s))
