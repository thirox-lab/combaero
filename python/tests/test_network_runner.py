"""Tests for NetworkRunner and NetworkResult.

Covers: from_file, from_dict, solve (no overrides), BC overrides, sweep,
NetworkResult.get, node_state, to_dataframe, and the multi_wall_tee.json
integration case (python/tests/data, #783). The sweep's results table
itself is tested in test_results_table.py (#784).
"""

import math
import warnings
from pathlib import Path

import pandas as pd
import pytest

from gui.backend.runner import NetworkResult, NetworkRunner

_FIXTURE_DIR = Path(__file__).parent / "data"
_MULTI_WALL_TEE = _FIXTURE_DIR / "multi_wall_tee.json"

# ---------------------------------------------------------------------------
# Minimal inline fixture -- one mass boundary, one channel, one pressure boundary.
# Uses the incompressible (default) solver for speed.
# ---------------------------------------------------------------------------

_SIMPLE: dict = {
    "nodes": [
        {
            "id": "inlet",
            "type": "mass_boundary",
            "position": {"x": 0, "y": 0},
            "data": {
                "label": "inlet",
                "m_dot": 1.0,
                "Tt": 300.0,
                "composition": {"source": "dry_air"},
            },
        },
        {
            "id": "chan",
            "type": "channel",
            "position": {"x": 200, "y": 0},
            "data": {"label": "chan", "L": 1.0, "D": 0.1, "roughness": 1e-5},
        },
        {
            "id": "outlet",
            "type": "pressure_boundary",
            "position": {"x": 400, "y": 0},
            "data": {
                "label": "outlet",
                "Pt": 101325.0,
                "Tt": 300.0,
                "composition": {"source": "dry_air"},
            },
        },
    ],
    "edges": [
        {"id": "e1", "source": "inlet", "target": "chan"},
        {"id": "e2", "source": "chan", "target": "outlet"},
    ],
}


# Fixture with two pressure boundaries so one can be swapped to mass without
# leaving the network without a pressure reference (which the solver requires).
_DUAL_PRESSURE: dict = {
    "nodes": [
        {
            "id": "p_inlet",
            "type": "pressure_boundary",
            "position": {"x": 0, "y": 0},
            "data": {
                "label": "p_inlet",
                "Pt": 200_000.0,
                "Tt": 300.0,
                "composition": {"source": "dry_air"},
            },
        },
        {
            "id": "chan",
            "type": "channel",
            "position": {"x": 200, "y": 0},
            "data": {"label": "chan", "L": 1.0, "D": 0.1, "roughness": 1e-5},
        },
        {
            "id": "p_outlet",
            "type": "pressure_boundary",
            "position": {"x": 400, "y": 0},
            "data": {
                "label": "p_outlet",
                "Pt": 101_325.0,
                "Tt": 300.0,
                "composition": {"source": "dry_air"},
            },
        },
    ],
    "edges": [
        {"id": "e1", "source": "p_inlet", "target": "chan"},
        {"id": "e2", "source": "chan", "target": "p_outlet"},
    ],
}


# Fixture with UUID-style IDs and separate labels (tests label-based resolution)
_UUID_LABELED: dict = {
    "nodes": [
        {
            "id": "node_abc123",
            "type": "mass_boundary",
            "position": {"x": 0, "y": 0},
            "data": {
                "label": "air_inlet",
                "m_dot": 1.0,
                "Tt": 300.0,
                "composition": {"source": "dry_air"},
            },
        },
        {
            "id": "node_def456",
            "type": "channel",
            "position": {"x": 200, "y": 0},
            "data": {"label": "hot_channel", "L": 1.0, "D": 0.1, "roughness": 1e-5},
        },
        {
            "id": "node_ghi789",
            "type": "pressure_boundary",
            "position": {"x": 400, "y": 0},
            "data": {
                "label": "exhaust",
                "Pt": 101325.0,
                "Tt": 300.0,
                "composition": {"source": "dry_air"},
            },
        },
    ],
    "edges": [
        {"id": "e1", "source": "node_abc123", "target": "node_def456"},
        {"id": "e2", "source": "node_def456", "target": "node_ghi789"},
    ],
}


# ---------------------------------------------------------------------------
# Construction tests
# ---------------------------------------------------------------------------


def test_from_dict():
    runner = NetworkRunner.from_dict(_SIMPLE)
    assert isinstance(runner, NetworkRunner)


def test_from_file(tmp_path):
    import json

    p = tmp_path / "simple.json"
    p.write_text(json.dumps(_SIMPLE))
    runner = NetworkRunner.from_file(p)
    assert isinstance(runner, NetworkRunner)


def test_from_dict_normalises_solver_settings_camelcase():
    """GUI exports use camelCase 'solverSettings'; from_dict must translate it."""
    import copy

    d = copy.deepcopy(_SIMPLE)
    d["solverSettings"] = {"global_regime": "compressible", "method": "hybr"}
    runner = NetworkRunner.from_dict(d)
    assert runner._schema.solver_settings.global_regime == "compressible"


def test_from_dict_snake_case_solver_settings_unchanged():
    """If snake_case solver_settings is already present it takes precedence."""
    import copy

    d = copy.deepcopy(_SIMPLE)
    d["solver_settings"] = {"global_regime": "incompressible"}
    d["solverSettings"] = {"global_regime": "compressible"}
    runner = NetworkRunner.from_dict(d)
    assert runner._schema.solver_settings.global_regime == "incompressible"


# ---------------------------------------------------------------------------
# Solve -- basic
# ---------------------------------------------------------------------------


def test_solve_success():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    assert result.success
    assert result.final_norm is not None
    assert result.final_norm < 1e-3


def test_solve_returns_network_result():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    assert isinstance(result, NetworkResult)


def test_solve_mass_conservation():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    m_dot = result.get("chan.m_dot")
    assert math.isclose(m_dot, 1.0, rel_tol=1e-3)


# ---------------------------------------------------------------------------
# Solve -- statelessness (repeated calls must not interfere)
# ---------------------------------------------------------------------------


def test_solve_stateless():
    runner = NetworkRunner.from_dict(_SIMPLE)
    r1 = runner.solve({"inlet.m_dot": 0.5})
    r2 = runner.solve({"inlet.m_dot": 1.5})
    r3 = runner.solve()
    assert math.isclose(result_m_dot(r1), 0.5, rel_tol=1e-3)
    assert math.isclose(result_m_dot(r2), 1.5, rel_tol=1e-3)
    assert math.isclose(result_m_dot(r3), 1.0, rel_tol=1e-3)


def result_m_dot(result: NetworkResult) -> float:
    return result.get("chan.m_dot")


# ---------------------------------------------------------------------------
# Overrides
# ---------------------------------------------------------------------------


def test_override_m_dot():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve({"inlet.m_dot": 2.0})
    assert result.success
    assert math.isclose(result.get("chan.m_dot"), 2.0, rel_tol=1e-3)


def test_override_by_node_id_fallback():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve({"inlet.m_dot": 0.7})
    assert math.isclose(result.get("chan.m_dot"), 0.7, rel_tol=1e-3)


def test_override_missing_label_raises():
    runner = NetworkRunner.from_dict(_SIMPLE)
    with pytest.raises(KeyError, match="no_such_node"):
        runner.solve({"no_such_node.m_dot": 1.0})


def test_override_bad_key_format_raises():
    runner = NetworkRunner.from_dict(_SIMPLE)
    with pytest.raises(ValueError, match="<label>.<attribute>"):
        runner.solve({"m_dot": 1.0})


# ---------------------------------------------------------------------------
# NetworkResult API
# ---------------------------------------------------------------------------


def test_result_get_known_key():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    v = result.get("chan.m_dot")
    assert isinstance(v, float)


def test_result_get_missing_key_raises():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    with pytest.raises(KeyError):
        result.get("nonexistent.key")


def test_result_get_by_label_when_id_differs():
    """get() must resolve '<label>.<quantity>' when label != node id."""
    runner = NetworkRunner.from_dict(_UUID_LABELED)
    result = runner.solve()
    # 'hot_channel' is the label for node id 'node_def456'
    v = result.get("hot_channel.m_dot")
    assert isinstance(v, float)
    assert math.isclose(v, 1.0, rel_tol=1e-3)


def test_result_get_label_fallback_does_not_shadow_raw_key():
    """A key present verbatim in raw must be returned without label resolution."""
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    v_direct = result.get("chan.m_dot")
    v_label = result.get("chan.m_dot")
    assert v_direct == v_label


def test_node_state_by_label():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    state = result.node_state("outlet")
    assert "T" in state
    assert "P" in state
    assert "Y" in state


def test_node_state_missing_label_raises():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    with pytest.raises(KeyError):
        result.node_state("no_such")


# ---------------------------------------------------------------------------
# to_dataframe
# ---------------------------------------------------------------------------


def test_to_dataframe_returns_dataframe():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    df = result.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0


def test_to_dataframe_has_type_column():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    df = result.to_dataframe()
    assert "type" in df.columns


def test_to_dataframe_has_element_row():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    df = result.to_dataframe()
    assert "element" in df["type"].values


def test_to_dataframe_unit_annotated_columns():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    df = result.to_dataframe()
    assert "T [K]" in df.columns
    assert "P [Pa]" in df.columns
    assert "m_dot [kg/s]" in df.columns


# ---------------------------------------------------------------------------
# sweep
# ---------------------------------------------------------------------------


def test_sweep_metrics_shape():
    # One row per case, (kind, entity, quantity) columns (#784).
    runner = NetworkRunner.from_dict(_SIMPLE)
    params = pd.DataFrame({"inlet.m_dot": [0.5, 1.0, 1.5]})
    df = runner.sweep(params, metrics=["chan.m_dot"])
    assert len(df) == 3
    assert ("outputs", "chan", "m_dot") in df.columns
    assert ("inputs", "inlet", "m_dot") in df.columns
    assert ("outputs", "_solve", "success") in df.columns
    # metrics keeps only the named outputs (and the solve status)
    assert set(df["outputs"].columns.get_level_values(0)) == {"_solve", "chan"}


def test_sweep_metrics_values():
    runner = NetworkRunner.from_dict(_SIMPLE)
    params = pd.DataFrame({"inlet.m_dot": [0.5, 1.0, 2.0]})
    df = runner.sweep(params, metrics=["chan.m_dot"])
    for i, expected in enumerate([0.5, 1.0, 2.0]):
        assert math.isclose(df["outputs"].loc[i, ("chan", "m_dot")], expected, rel_tol=1e-3)


def test_sweep_no_metrics_returns_full_detail():
    runner = NetworkRunner.from_dict(_SIMPLE)
    params = pd.DataFrame({"inlet.m_dot": [0.5, 1.0]})
    df = runner.sweep(params)
    assert len(df) == 2
    entities = set(df["outputs"].columns.get_level_values(0))
    assert {"inlet", "chan", "outlet", "_solve"} <= entities
    assert ("outputs", "outlet", "T") in df.columns


def test_sweep_empty_params():
    runner = NetworkRunner.from_dict(_SIMPLE)
    params = pd.DataFrame({"inlet.m_dot": []})
    df = runner.sweep(params, metrics=["chan.m_dot"])
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 0


def test_sweep_solves_cold_when_the_warm_seed_fails():
    # #1054: the previous row's solution carries absolute pressures, so a jump
    # from 1.5 to 6 bar seeds the 6 bar row far from its root, and the warm
    # solve stalls (no_progress on every platform at this jump). The row is
    # solved cold instead, and the row after it continues from it.
    import gui

    example = Path(gui.__file__).parent / "examples" / "iso5167_orifice_metering.json"
    runner = NetworkRunner.from_file(example)
    params = pd.DataFrame(
        {"supply.Pt": [1.5e5, 6.0e5, 6.0e5], "sink.Pt": [1.4e5, 5.99e5, 5.975e5]},
        index=["low", "jump", "next"],
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        df = runner.sweep(params, metrics=["plate.m_dot"])
    assert df["outputs"]["_solve", "success"].all()
    assert df.attrs["warm_start"] == {"low": None, "jump": None, "next": "jump"}
    cold = runner.solve({"supply.Pt": 6.0e5, "sink.Pt": 5.99e5})
    assert df["outputs"].loc["jump"].xs("m_dot", level="quantity").iloc[0] == pytest.approx(
        cold.get("plate.m_dot"), rel=1e-6
    )


# ---------------------------------------------------------------------------
# swap_boundary
# ---------------------------------------------------------------------------


def test_swap_mass_to_pressure_returns_runner():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    p_runner = result.swap_boundary("inlet", "pressure_boundary")
    assert isinstance(p_runner, NetworkRunner)


def test_swap_mass_to_pressure_pt_in_schema():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    solved_Pt = result.node_state("inlet")["Pt"]
    p_runner = result.swap_boundary("inlet", "pressure_boundary")
    inlet_node = next(n for n in p_runner._schema.nodes if n.id == "inlet")
    assert math.isclose(inlet_node.data["Pt"], solved_Pt, rel_tol=1e-6)


def test_swap_mass_to_pressure_resolves_same():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    p_runner = result.swap_boundary("inlet", "pressure_boundary")
    result2 = p_runner.solve()
    assert math.isclose(result.get("chan.m_dot"), result2.get("chan.m_dot"), rel_tol=1e-3)
    assert math.isclose(
        result.node_state("outlet")["T"], result2.node_state("outlet")["T"], rel_tol=1e-3
    )


def test_swap_pressure_to_mass_returns_runner():
    # _DUAL_PRESSURE has two pressure BCs; swapping the inlet keeps a pressure reference.
    runner = NetworkRunner.from_dict(_DUAL_PRESSURE)
    result = runner.solve()
    m_runner = result.swap_boundary("p_inlet", "mass_boundary")
    assert isinstance(m_runner, NetworkRunner)


def test_swap_pressure_to_mass_mdot_in_schema():
    runner = NetworkRunner.from_dict(_DUAL_PRESSURE)
    result = runner.solve()
    solved_m_dot = result.get("chan.m_dot")
    m_runner = result.swap_boundary("p_inlet", "mass_boundary")
    inlet_node = next(n for n in m_runner._schema.nodes if n.id == "p_inlet")
    assert math.isclose(inlet_node.data["m_dot"], solved_m_dot, rel_tol=1e-3)


def test_swap_pressure_to_mass_resolves_same():
    runner = NetworkRunner.from_dict(_DUAL_PRESSURE)
    result = runner.solve()
    m_runner = result.swap_boundary("p_inlet", "mass_boundary")
    result2 = m_runner.solve()
    assert math.isclose(result.get("chan.m_dot"), result2.get("chan.m_dot"), rel_tol=1e-3)
    assert math.isclose(
        result.node_state("p_outlet")["T"], result2.node_state("p_outlet")["T"], rel_tol=1e-3
    )


def test_swap_same_type_raises():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    with pytest.raises(ValueError, match="already"):
        result.swap_boundary("inlet", "mass_boundary")


def test_swap_unsupported_new_type_raises():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    with pytest.raises(ValueError):
        result.swap_boundary("inlet", "plenum")


def test_swap_unsupported_old_type_raises():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    with pytest.raises(ValueError, match="cannot be swapped"):
        result.swap_boundary("chan", "pressure_boundary")


def test_swap_missing_label_raises():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    with pytest.raises(KeyError):
        result.swap_boundary("no_such_node", "pressure_boundary")


# ---------------------------------------------------------------------------
# Integration -- multi_wall_tee.json
# ---------------------------------------------------------------------------


def test_multi_wall_tee_solve():
    # init_strategy="default": the fixture's stored results are from a
    # FAILED solve, so since the warm-start success gate they no longer
    # seed the solver -- and the incompressible_warmstart (removed in 0.8.0) this
    # test used to pin only converged FROM that poisoned seed. The plain
    # default converges cold.
    runner = NetworkRunner.from_file(_MULTI_WALL_TEE)
    result = runner.solve(
        method="hybr",
        init_strategy="default",
        timeout=60.0,
    )
    assert result.success, f"Solver failed: {result.message}"
    assert result.final_norm < 1e-3


def test_multi_wall_tee_to_dataframe():
    runner = NetworkRunner.from_file(_MULTI_WALL_TEE)
    result = runner.solve(
        method="hybr",
        init_strategy="default",
        timeout=60.0,
    )
    df = result.to_dataframe()
    assert len(df) > 0
    assert "thermal_wall" in df["type"].values
