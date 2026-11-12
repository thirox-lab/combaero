"""GUI 'film_cooling' node and a thermal wall's film (#386).

The node is the holes; a thermal wall edge on the gas duct downstream names
it (``film``), and the edge comes back with the film the wall saw. The
default node is Colban's Eq. (19) at Gritsch configuration B, the geometry it
is scored on.
"""

from __future__ import annotations

import warnings

import pytest
from _asgi_client import request

from combaero.network import FilmCoolingElement
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.main import app
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


def test_the_hole_angle_cd_term_threads_through() -> None:
    assert FilmCoolingData().hole_cd.inclination == "auto"
    # The node's default hole is shaped, which "auto" leaves without the term;
    # a cylindrical hole gets it.
    shaped = build_network_from_schema(NetworkGraphSchema(**_schema(FILM))).elements["film"]
    assert shaped.hole_type == "shaped" and shaped.cd_inclination == "none"
    cyl = build_network_from_schema(
        NetworkGraphSchema(**_schema(FILM, {"hole_type": "cylindrical", "lateral_pitch": 3e-3}))
    ).elements["film"]
    assert cyl.cd_inclination == "idelchik_3_2"
    on = build_network_from_schema(
        NetworkGraphSchema(**_schema(FILM, {"hole_cd": {"inclination": "idelchik_3_2"}}))
    ).elements["film"]
    assert on.cd_inclination == "idelchik_3_2"


CYL = {"hole_type": "cylindrical", "lateral_pitch": 3e-3}


def test_auto_applies_the_hole_angle_term_to_an_inclined_cylindrical_hole() -> None:
    # The builder hands the node's "auto" to the element, which resolves it.
    e = build_network_from_schema(NetworkGraphSchema(**_schema(FILM, CYL))).elements["film"]
    assert e.angle_deg == 30.0 and e.correlation == "IdelchikThick"
    assert e.cd_inclination == "idelchik_3_2" and e._hole_cd.inclination_deg == 30.0
    assert e._hole_cd.inclination_dzeta() > 0.0


def test_auto_leaves_a_shaped_hole_without_the_hole_angle_term() -> None:
    e = build_network_from_schema(NetworkGraphSchema(**_schema(FILM))).elements["film"]
    assert e.hole_type == "shaped" and e.correlation == "IdelchikThick"
    assert e.cd_inclination == "none" and e._hole_cd.inclination_deg == 90.0
    assert e._hole_cd.inclination_dzeta() == 0.0


@pytest.mark.parametrize("hole_type", ["shaped", "cylindrical"])
@pytest.mark.parametrize("correlation", ["fixed", "Lichtarowicz", "McGreehanSchotsch"])
def test_an_explicit_hole_angle_term_the_correlation_cannot_take_is_refused(
    hole_type: str, correlation: str
) -> None:
    # The builder used to turn this into "none" silently; the element refuses it.
    holes = {
        **(CYL if hole_type == "cylindrical" else {}),
        "hole_cd": {
            "correlation": correlation,
            **({"Cd": 0.6} if correlation == "fixed" else {}),
            "inclination": "idelchik_3_2",
        },
    }
    with pytest.raises(ValueError, match="does not compose with correlation"):
        build_network_from_schema(NetworkGraphSchema(**_schema(FILM, holes)))


def test_the_refused_hole_angle_term_reaches_the_user_as_a_readable_400() -> None:
    holes = {"hole_cd": {"correlation": "fixed", "Cd": 0.6, "inclination": "idelchik_3_2"}}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = request(app, "POST", "/solve", json_body=_schema(FILM, holes))
    assert r.status_code == 400, r.text[:300]
    detail = r.json()["detail"]
    assert "'film' [" in detail, detail
    assert "does not compose with correlation='fixed'" in detail, detail
    assert "Hint: the Hole Cd section's Hole Angle Term" in detail, detail


def test_a_film_naming_something_that_is_not_a_provider_is_refused() -> None:
    with pytest.raises(ValueError, match="not a film provider"):
        build_network_from_schema(NetworkGraphSchema(**_schema({**FILM, "film": "hin"})))


def test_the_exit_crossflow_cd_term_threads_through() -> None:
    assert FilmCoolingData().hole_cd.exit_crossflow == "auto"
    off = build_network_from_schema(NetworkGraphSchema(**_schema())).elements["film"]
    assert off.cd_exit_crossflow == "none" and off._hole_cd.exit_xf is None
    on = build_network_from_schema(
        NetworkGraphSchema(
            **_schema(
                holes={
                    "hole_cd": {"exit_crossflow": "rowbury_2001"},
                    "exit_crossflow_region": "pressure_surface",
                    "exit_crossflow_orientation_deg": 5.0,
                }
            )
        )
    ).elements["film"]
    assert on.cd_exit_crossflow == "rowbury_2001"
    assert on.exit_crossflow_region == "pressure_surface"
    assert on.exit_crossflow_orientation_deg == 5.0
