"""GUI 'slot_film' node: a 2-D film slot as a thermal wall's film (#714).

The node is the slot; a thermal wall edge on the gas duct downstream names
it (``film``), and the edge comes back with the film the wall saw. The
default node is Marek & Tacina's D-7958 rig, on which the film is scored.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from combaero.network import SlotFilmElement
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import NetworkGraphSchema, SlotFilmData

REPO = Path(__file__).resolve().parents[2]


def _node(nid: str, ntype: str, data: dict) -> dict:
    return {"id": nid, "type": ntype, "position": {"x": 0, "y": 0}, "data": data}


def _edge(eid: str, src: str, tgt: str, data: dict | None = None) -> dict:
    return {"id": eid, "source": src, "target": tgt, "data": data or {}}


def _schema(film: dict | None = None, slot: dict | None = None) -> dict:
    """Hot gas into a plenum the slot discharges into, then a duct walled
    against a cold duct (the film_cooling GUI test's network)."""
    wall = {"type": "thermal", "thickness": 0.001, "conductivity": 20.0, "area": 0.05}
    wall.update(film or {})
    return {
        "nodes": [
            _node("hi", "pressure_boundary", {"Pt": 1.03e5, "Tt": 1600.0}),
            _node("ho", "pressure_boundary", {"Pt": 1.0e5, "Tt": 1600.0}),
            _node("cool", "pressure_boundary", {"Pt": 1.05e5, "Tt": 700.0}),
            _node("ci", "pressure_boundary", {"Pt": 1.02e5, "Tt": 400.0}),
            _node("co", "pressure_boundary", {"Pt": 1.0e5, "Tt": 400.0}),
            _node("hin", "orifice", {"diameter": 0.03, "correlation": "fixed", "Cd": 0.8}),
            _node("g", "plenum", {}),
            _node(
                "slot",
                "slot_film",
                slot if slot is not None else {"slot_height": 1.0e-3, "span": 0.05},
            ),
            _node("hch", "channel", {"L": 0.2, "D": 0.04}),
            _node("cch", "channel", {"L": 0.2, "D": 0.02}),
        ],
        "edges": [
            _edge("h0", "hi", "hin"),
            _edge("h1", "hin", "g"),
            _edge("f1", "cool", "slot"),
            _edge("f2", "slot", "g"),
            _edge("h2", "g", "hch"),
            _edge("h3", "hch", "ho"),
            _edge("c1", "ci", "cch"),
            _edge("c2", "cch", "co"),
            _edge("w", "hch", "cch", wall),
        ],
    }


FILM = {"film": "slot", "film_side": "a", "film_start": 0.0, "film_end": 0.03}


def test_the_default_node_is_the_d7958_rig() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_schema(FILM, slot={})))
    e = net.elements["slot"]
    assert isinstance(e, SlotFilmElement)
    assert (e.slot_height, e.span, e.turbulence_intensity) == (6.3e-3, 0.37, 0.14)
    assert e.C_m_multiplier == 1.0
    assert net.walls["w"].film == "slot"


def test_the_canvas_drops_the_same_defaults_as_the_schema() -> None:
    """NetworkCanvas's new-node data and SlotFilmData describe one rig."""
    src = (REPO / "gui/frontend/src/components/NetworkCanvas.tsx").read_text()
    block = src[src.index('type === "slot_film"') :]
    block = block[: block.index("} as any;")]
    d = SlotFilmData()
    for key in (
        "slot_height",
        "span",
        "lip_thickness",
        "turbulence_intensity",
        "C_m_multiplier",
    ):
        m = re.search(rf"{key}: ([0-9.e-]+),", block)
        assert m, key
        assert float(m.group(1)) == pytest.approx(getattr(d, key)), key
    # The hole Cd (#1098): the schema's None is the element's IdelchikThick.
    assert d.hole_cd.correlation is None
    assert 'hole_cd: { correlation: "IdelchikThick" }' in block


def test_the_wall_edge_comes_back_with_the_slot_film() -> None:
    bare = NetworkRunner.from_dict(_schema()).solve(timeout=120.0)
    filmed = NetworkRunner.from_dict(_schema(FILM)).solve(timeout=120.0)
    assert bare.success, bare.message
    assert filmed.success, filmed.message
    r = filmed._edge_results["w"]
    assert 0.0 < r["eta_film"] < 1.0
    assert r["T_aw_a"] == pytest.approx(
        r["T_aw_unfilmed"] - r["eta_film"] * (r["T_aw_unfilmed"] - r["T_coolant_film"])
    )
    assert r["Q"] < bare._edge_results["w"]["Q"]
    assert filmed._element_results["slot"].model_dump()["m_dot"] > 0.0
