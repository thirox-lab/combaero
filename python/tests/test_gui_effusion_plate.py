"""GUI 'effusion_plate' node (#471): the default is a working example.

Discharged into a momentum chamber's side handle the plate takes its gas side
from the chamber; discharged into a plenum it takes the imposed heat flux.
"""

from __future__ import annotations

import pytest

from combaero.network import EffusionPlateElement
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import EffusionPlateData, NetworkGraphSchema


def _node(nid, ntype, data):
    return {"id": nid, "type": ntype, "position": {"x": 0, "y": 0}, "data": data}


def _edge(eid, src, tgt, tgt_handle=None):
    e = {"id": eid, "source": src, "target": tgt, "data": {}}
    if tgt_handle:
        e["targetHandle"] = tgt_handle
    return e


def _liner(plate: dict | None = None) -> dict:
    return {
        "nodes": [
            _node("gas_in", "pressure_boundary", {"Pt": 1.003e5, "Tt": 1500.0}),
            _node("exit", "pressure_boundary", {"Pt": 1.0e5, "Tt": 1500.0}),
            _node("cool", "pressure_boundary", {"Pt": 1.03e5, "Tt": 600.0}),
            _node("duct", "channel", {"L": 0.2, "D": 0.08}),
            _node("tail", "channel", {"L": 0.2, "D": 0.08}),
            _node("ch", "momentum_chamber", {}),
            _node("eff", "effusion_plate", plate or {}),
        ],
        "edges": [
            _edge("e1", "gas_in", "duct"),
            _edge("e2", "duct", "ch", "flow-target"),
            _edge("e3", "ch", "tail"),
            _edge("e4", "tail", "exit"),
            _edge("e5", "cool", "eff"),
            _edge("e6", "eff", "ch", "side-target"),
        ],
    }


def test_the_node_builds_andrews_plate_c_by_default() -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_liner()))
    e = net.elements["eff"]
    assert isinstance(e, EffusionPlateElement)
    d = EffusionPlateData()
    assert e.hole_diameter == d.hole_diameter == 3.27e-3
    assert e.pitch_x == e.pitch_y == 15.24e-3
    assert e.n_holes == round(0.152 * 0.152 / 15.24e-3**2)
    assert net.nodes["ch"].main_inlet == "duct"


def _liner_in(regime: str) -> dict:
    g = _liner()
    g["solver_settings"] = {"global_regime": regime}
    return g


def test_the_default_liner_solves_with_a_chamber_gas_side() -> None:
    # Incompressible: the compressible solve of this liner stalls (#585).
    result = NetworkRunner.from_dict(_liner_in("incompressible")).solve(timeout=120.0)
    r = result._element_results["eff"].model_dump()
    assert r["success"]
    assert r["gas_side_chamber"] == 1.0
    assert 600.0 < r["T_wall_cold"] <= r["T_wall_hot"] < r["T_gas"]
    assert 0.0 < r["eta_overall"] < 1.0


def test_the_default_liner_solves_in_the_gui_default_regime() -> None:
    """Compressible is the default to come (#585). At M ~ 0.04 it agrees with
    the incompressible solve: measured 4e-6 apart in the plate flow. It used to
    stall; its causes were fixed in #582, #587 and #554 (orifices read the
    upstream side along the flow)."""
    comp = NetworkRunner.from_dict(_liner_in("compressible")).solve(timeout=120.0)
    inc = NetworkRunner.from_dict(_liner_in("incompressible")).solve(timeout=120.0)
    r_c = comp._element_results["eff"].model_dump()
    r_i = inc._element_results["eff"].model_dump()
    assert r_c["success"] and r_i["success"]
    assert r_c["m_dot"] == pytest.approx(r_i["m_dot"], rel=1e-4)


def test_a_plenum_discharge_uses_the_imposed_flux() -> None:
    schema = {
        "nodes": [
            _node("cool", "pressure_boundary", {"Pt": 1.03e5, "Tt": 600.0}),
            _node("burner", "pressure_boundary", {"Pt": 1.0e5, "Tt": 1500.0}),
            _node("eff", "effusion_plate", {"gas_heat_flux": 1.0e5}),
        ],
        "edges": [_edge("e1", "cool", "eff"), _edge("e2", "eff", "burner")],
    }
    result = NetworkRunner.from_dict(schema).solve(timeout=120.0)
    r = result._element_results["eff"].model_dump()
    assert r["gas_side_chamber"] == 0.0
    assert r["q_wall"] == pytest.approx(1.0e5)
    assert r["T_wall_hot"] > r["T_wall_cold"] > 600.0


def test_the_flat_plate_gas_baseline_is_exposed_and_the_chamber_stays_default() -> None:
    """The node carries the element's gas_baseline; its length is read only
    for the flat plate, so a stored length with the chamber baseline builds."""
    assert EffusionPlateData().gas_baseline == "chamber"
    e = build_network_from_schema(NetworkGraphSchema(**_liner())).elements["eff"]
    assert e.gas_baseline == "chamber" and e.gas_start_length is None
    kept = build_network_from_schema(
        NetworkGraphSchema(**_liner({"gas_baseline": "chamber", "gas_start_length": 0.02}))
    ).elements["eff"]
    assert kept.gas_baseline == "chamber" and kept.gas_start_length is None
    g = _liner_in("incompressible")
    g["nodes"][-1]["data"] = {"gas_baseline": "flat_plate", "gas_start_length": 0.02}
    r = NetworkRunner.from_dict(g).solve(timeout=120.0)._element_results["eff"].model_dump()
    assert r["success"]
    assert r["h_gas_unblown_x_start"] == pytest.approx(0.02)
    assert r["h_gas_unblown_x_end"] == pytest.approx(0.02 + 0.152)


def test_the_hole_angle_cd_term_is_on_by_default_and_threads_through() -> None:
    auto = build_network_from_schema(NetworkGraphSchema(**_liner({"angle_deg": 30.0})))
    assert EffusionPlateData().cd_inclination == "auto"
    assert auto.elements["eff"].cd_inclination == "idelchik_3_2"
    assert auto.elements["eff"]._hole_cd.inclination_deg == 30.0
    off = build_network_from_schema(
        NetworkGraphSchema(**_liner({"angle_deg": 30.0, "cd_inclination": "none"}))
    )
    assert off.elements["eff"].cd_inclination == "none"
    on = build_network_from_schema(
        NetworkGraphSchema(**_liner({"angle_deg": 30.0, "cd_inclination": "idelchik_3_2"}))
    )
    assert on.elements["eff"].cd_inclination == "idelchik_3_2"
    assert on.elements["eff"]._hole_cd.inclination_deg == 30.0


def test_auto_applies_the_hole_angle_term_on_idelchik_thick_only() -> None:
    on = build_network_from_schema(NetworkGraphSchema(**_liner({"angle_deg": 30.0})))
    assert on.elements["eff"]._hole_cd.inclination_diagnostics()["Cd_inclination_dzeta"] > 0.0
    off = build_network_from_schema(
        NetworkGraphSchema(**_liner({"angle_deg": 30.0, "correlation": "Lichtarowicz"}))
    )
    assert off.elements["eff"].cd_inclination == "none"
    assert off.elements["eff"]._hole_cd.inclination_diagnostics()["Cd_inclination_dzeta"] == 0.0


def test_an_explicit_hole_angle_term_off_idelchik_thick_is_refused() -> None:
    # The builder used to turn this into "none" silently; the element refuses it.
    with pytest.raises(ValueError, match="does not compose with correlation='Lichtarowicz'"):
        build_network_from_schema(
            NetworkGraphSchema(
                **_liner({"correlation": "Lichtarowicz", "cd_inclination": "idelchik_3_2"})
            )
        )


def test_the_exit_crossflow_cd_term_threads_through() -> None:
    off = build_network_from_schema(NetworkGraphSchema(**_liner({"angle_deg": 30.0})))
    assert off.elements["eff"].cd_exit_crossflow == "none"
    on = build_network_from_schema(
        NetworkGraphSchema(**_liner({"angle_deg": 30.0, "cd_exit_crossflow": "rowbury_2001"}))
    )
    assert on.elements["eff"].cd_exit_crossflow == "rowbury_2001"
    assert on.elements["eff"].exit_crossflow_region == "suction_surface"
