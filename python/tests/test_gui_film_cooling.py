"""GUI 'film_cooling' node and a thermal wall's film (#386).

The node is the holes; a thermal wall edge on the gas duct downstream names
it (``film``), and the edge comes back with the film the wall saw. The
default node is Colban's Eq. (19) at Gritsch configuration B, the geometry it
is scored on.
"""

from __future__ import annotations

import pytest

from combaero.network import FilmCoolingElement
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import FilmCoolingData, NetworkGraphSchema


def _node(nid: str, ntype: str, data: dict) -> dict:
    return {"id": nid, "type": ntype, "position": {"x": 0, "y": 0}, "data": data}


def _edge(eid: str, src: str, tgt: str, data: dict | None = None) -> dict:
    return {"id": eid, "source": src, "target": tgt, "data": data or {}}


def _schema(film: dict | None = None, holes: dict | None = None) -> dict:
    """Hot gas into a plenum the holes discharge into, then a duct walled
    against a cold duct."""
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
            _node("film", "film_cooling", holes or {}),
            _node("hch", "channel", {"L": 0.2, "D": 0.04}),
            _node("cch", "channel", {"L": 0.2, "D": 0.02}),
        ],
        "edges": [
            _edge("h0", "hi", "hin"),
            _edge("h1", "hin", "g"),
            _edge("f1", "cool", "film"),
            _edge("f2", "film", "g"),
            _edge("h2", "g", "hch"),
            _edge("h3", "hch", "ho"),
            _edge("c1", "ci", "cch"),
            _edge("c2", "cch", "co"),
            _edge("w", "hch", "cch", wall),
        ],
    }


FILM = {"film": "film", "film_side": "a", "film_start": 0.0, "film_end": 0.03}


def test_the_node_builds_colban_configuration_b_by_default() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_schema(FILM)))
    e = net.elements["film"]
    assert isinstance(e, FilmCoolingElement)
    d = FilmCoolingData()
    assert e.hole_type == "shaped"
    assert e.lateral_pitch / e.hole_diameter == pytest.approx(6.0)
    assert e.breakout_width / e.lateral_pitch == pytest.approx(0.49)
    assert e.area_ratio == d.area_ratio == 3.5
    assert e.angle_deg == 30.0
    w = net.walls["w"]
    assert (w.film, w.film_side, w.film_start, w.film_end) == ("film", "a", 0.0, 0.03)


def test_a_saved_wall_without_a_film_builds_without_one() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_schema()))
    assert net.walls["w"].film is None


def test_the_wall_edge_comes_back_with_the_film_it_saw() -> None:
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
    assert "eta_film" not in bare._edge_results["w"]
    assert filmed._element_results["film"].model_dump()["m_dot"] > 0.0


def test_cylindrical_holes_drop_the_shaped_geometry() -> None:
    net = build_network_from_schema(
        NetworkGraphSchema(**_schema(FILM, {"hole_type": "cylindrical", "lateral_pitch": 3e-3}))
    )
    e = net.elements["film"]
    assert e.hole_type == "cylindrical"
    assert e.area_ratio is None and e.breakout_width is None


def test_the_hole_angle_cd_term_threads_through_with_idelchik_thick_only() -> None:
    assert FilmCoolingData().cd_inclination == "none"
    on = build_network_from_schema(
        NetworkGraphSchema(**_schema(FILM, {"cd_inclination": "idelchik_3_2"}))
    ).elements["film"]
    assert on.cd_inclination == "idelchik_3_2"
    kept = build_network_from_schema(
        NetworkGraphSchema(
            **_schema(FILM, {"cd_inclination": "idelchik_3_2", "correlation": "fixed"})
        )
    ).elements["film"]
    assert kept.cd_inclination == "none"


def test_a_film_naming_something_that_is_not_a_provider_is_refused() -> None:
    with pytest.raises(ValueError, match="not a film provider"):
        build_network_from_schema(NetworkGraphSchema(**_schema({**FILM, "film": "hin"})))


def test_the_exit_crossflow_cd_term_threads_through() -> None:
    assert FilmCoolingData().cd_exit_crossflow == "none"
    off = build_network_from_schema(NetworkGraphSchema(**_schema())).elements["film"]
    assert off.cd_exit_crossflow == "none" and off._exit_xf is None
    on = build_network_from_schema(
        NetworkGraphSchema(
            **_schema(
                holes={
                    "cd_exit_crossflow": "rowbury_2001",
                    "exit_crossflow_region": "pressure_surface",
                    "exit_crossflow_orientation_deg": 5.0,
                }
            )
        )
    ).elements["film"]
    assert on.cd_exit_crossflow == "rowbury_2001"
    assert on.exit_crossflow_region == "pressure_surface"
    assert on.exit_crossflow_orientation_deg == 5.0
