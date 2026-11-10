"""GUI radiation inputs (#1071): a thermal wall's two faces, an effusion
plate's and an effusion liner's hot face.

Kind "none", or no field at all (a file saved before it existed), builds no
radiation; a chosen kind builds the library's GasRadiation / CasingRadiation
on the face the edge names: side A is the edge's source, also when one end is
an impingement array.
"""

from __future__ import annotations

import pytest
import test_gui_effusion_liner as liner_t
import test_gui_effusion_plate as plate_t
import test_gui_impingement_array as imp_t
import test_gui_wall_orientation as wall_t
from pydantic import ValidationError

from combaero.network import CasingRadiation, GasRadiation
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import (
    EffusionLinerData,
    EffusionPlateData,
    NetworkGraphSchema,
    ThermalWallData,
)

_GAS = {
    "kind": "gas",
    "beam_length": 0.04,
    "emissivity_wall": 0.8,
    "luminosity": 2.0,
    "T_gas": 1700.0,
    "fuel_air_ratio": 0.03,
    "pressure": 2.0e5,
}
_CASING = {"kind": "casing", "area_ratio": 0.7, "emissivity_casing": 0.5, "T_casing": 500.0}


def _build(graph: dict):
    return build_network_from_schema(NetworkGraphSchema(**graph))


def _wall_graph(**radiation: dict) -> dict:
    graph = wall_t._schema(None)
    graph["edges"][-1]["data"].update(radiation)
    return graph


def _array_graph(array_is_source: bool, **radiation: dict) -> dict:
    graph = imp_t._schema(with_wall=True)
    wall = next(e for e in graph["edges"] if e["id"] == "wall")
    if array_is_source:
        wall["source"], wall["target"] = "arr", "hot"
    wall["data"].update(radiation)
    return graph


def test_old_payloads_without_the_fields_load_with_radiation_off() -> None:
    assert ThermalWallData().radiation_a is None and ThermalWallData().radiation_b is None
    assert EffusionPlateData().radiation is None
    assert EffusionLinerData().radiation is None
    net = _build(wall_t._schema(None))
    assert net.walls["w"].radiation_a is None and net.walls["w"].radiation_b is None
    assert _build(plate_t._liner()).elements["eff"].radiation is None


def test_kind_none_builds_the_same_objects_as_no_field() -> None:
    """A node or edge that carries the field with kind "none" (a new one, or
    one switched back) builds exactly what a file without it builds."""
    off = {"kind": "none", "beam_length": 0.1, "area_ratio": 0.5}
    plain = _build(wall_t._schema(None)).walls["w"]
    none = _build(_wall_graph(radiation_a=off, radiation_b=off)).walls["w"]
    assert vars(none).keys() == vars(plain).keys()
    for name in ("radiation_a", "radiation_b", "film", "orientation", "contact_area"):
        assert getattr(none, name) == getattr(plain, name)
    p_plain = _build(plate_t._liner()).elements["eff"]
    p_none = _build(plate_t._liner({"radiation": {"kind": "none"}})).elements["eff"]
    assert p_none.radiation is None
    assert vars(p_none).keys() == vars(p_plain).keys()


def test_a_wall_side_builds_the_chosen_spec_on_its_own_face() -> None:
    w = _build(_wall_graph(radiation_a=_GAS, radiation_b=_CASING)).walls["w"]
    assert w.radiation_a == GasRadiation(
        beam_length=0.04,
        emissivity_wall=0.8,
        luminosity=2.0,
        T_gas=1700.0,
        fuel_air_ratio=0.03,
        pressure=2.0e5,
    )
    assert w.radiation_b == CasingRadiation(
        area_ratio=0.7, emissivity_wall=0.7, emissivity_casing=0.5, T_casing=500.0
    )
    swapped = _build(_wall_graph(radiation_a=_CASING, radiation_b=_GAS)).walls["w"]
    assert isinstance(swapped.radiation_a, CasingRadiation)
    assert isinstance(swapped.radiation_b, GasRadiation)


def test_the_bordbar_emissivity_builds_without_the_lefebvre_fields() -> None:
    """emissivity_model "bordbar2021" (#1071): the node keeps its luminosity
    and fuel-air ratio, which the WSGG does not read, and builds without
    them; the default is Eq. 8.11."""
    gas = {**_GAS, "emissivity_model": "bordbar2021"}
    w = _build(_wall_graph(radiation_a=gas)).walls["w"]
    assert w.radiation_a == GasRadiation(
        beam_length=0.04,
        emissivity_wall=0.8,
        T_gas=1700.0,
        pressure=2.0e5,
        emissivity_model="bordbar2021",
    )
    assert (
        _build(_wall_graph(radiation_a=_GAS)).walls["w"].radiation_a.emissivity_model == "lefebvre"
    )
    with pytest.raises(ValueError, match="emissivity_model"):
        _build(_wall_graph(radiation_a={**_GAS, "emissivity_model": "wsgg"}))


def test_the_gui_defaults_are_the_libraries_with_working_geometry() -> None:
    w = _build(_wall_graph(radiation_a={"kind": "gas"}, radiation_b={"kind": "casing"})).walls["w"]
    assert w.radiation_a == GasRadiation(beam_length=0.05)
    assert w.radiation_b == CasingRadiation(area_ratio=0.8)


@pytest.mark.parametrize("array_is_source", [False, True])
def test_an_array_wall_keeps_side_a_on_the_edges_source(array_is_source: bool) -> None:
    """The gas spec on the hot channel's face, the casing spec on the array's
    target face, whichever end of the edge the array is."""
    hot_side, arr_side = (
        ("radiation_b", "radiation_a") if array_is_source else ("radiation_a", "radiation_b")
    )
    net = _build(_array_graph(array_is_source, **{hot_side: _GAS, arr_side: _CASING}))
    rows = [w for wid, w in net.walls.items() if wid.startswith("wall__r")]
    assert rows
    for w in rows:
        faces = {w.element_a: w.radiation_a, w.element_b: w.radiation_b}
        assert isinstance(faces["hot"], GasRadiation)
        plate = next(k for k in faces if k != "hot")
        assert plate.startswith("arr__p")
        assert isinstance(faces[plate], CasingRadiation)


def test_the_plate_takes_gas_radiation_and_refuses_casing() -> None:
    e = _build(plate_t._liner({"radiation": {"kind": "gas", "beam_length": 0.03}})).elements["eff"]
    assert e.radiation == GasRadiation(beam_length=0.03)
    # A plate face sees the gas: its schema has no casing kind.
    with pytest.raises(ValidationError):
        EffusionPlateData(radiation={"kind": "casing"})


def test_the_liner_passes_its_radiation_to_every_panel() -> None:
    net = _build(liner_t._schema({"radiation": {"kind": "gas", "beam_length": 0.06}}))
    panels = [e for e in net.elements.values() if type(e).__name__ == "EffusionPlateElement"]
    assert panels
    assert all(p.radiation == GasRadiation(beam_length=0.06) for p in panels)
    off = _build(liner_t._schema())
    assert all(
        e.radiation is None
        for e in off.elements.values()
        if type(e).__name__ == "EffusionPlateElement"
    )


def test_a_library_refusal_is_the_normal_build_error() -> None:
    with pytest.raises(ValueError, match="beam_length"):
        _build(_wall_graph(radiation_a={"kind": "gas", "beam_length": 0.0}))


def test_the_results_carry_the_radiation_keys() -> None:
    result = NetworkRunner.from_dict(
        _wall_graph(radiation_a={"kind": "gas"}, radiation_b={"kind": "casing"})
    ).solve(timeout=120.0)
    assert result.success, result.message
    edge = result._edge_results["w"]
    for key in ("q_rad_a", "Q_rad_a", "T_surface_a", "eps_g_a", "q_rad_b", "Q_rad_b"):
        assert key in edge, key
    assert "eps_g_b" not in edge  # a casing face has no gas emissivity
    plate = NetworkRunner.from_dict(plate_t._liner({"radiation": {"kind": "gas"}})).solve(
        timeout=120.0
    )
    assert plate.success, plate.message
    r = plate._element_results["eff"].model_dump()
    assert r["q_rad_hot"] > 0.0 and 0.0 < r["eps_g"] < 1.0
