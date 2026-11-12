"""Every name deprecated in 0.7.0 is gone in 0.8.0; every 0.8.0 deprecation is
gone in 0.9.0.

0.7.0 made each deprecated name warn and name its removal release; 0.8.0
removes them. This is the inventory: each removed function, alias, element,
keyword and strategy must fail the way a missing name fails (AttributeError,
ImportError, TypeError, ValueError), not keep working silently, and the
replacements must stay.

The second part is the inventory for the names deprecated in 0.8.0 (the
stagnation spellings, #782, Stolz #876, Al 6061 #761) and removed in 0.9.0:
each old spelling fails, and its replacement works without a warning.
"""

from __future__ import annotations

import dataclasses
import warnings

import pytest

import combaero as cb
import combaero.network as net
from combaero import _core, _solver_tools
from combaero.network import (
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
    VortexElement,
    pressure_loss,
)


@pytest.mark.parametrize(
    "name",
    ["standard_dry_air_composition", "smr_wgs_equilibrium", "smr_wgs_equilibrium_adiabatic"],
)
def test_the_removed_functions_are_gone(name: str) -> None:
    assert name not in cb.__all__
    with pytest.raises(AttributeError):
        getattr(cb, name)
    with pytest.raises(AttributeError):
        getattr(_core, name)


def test_the_replacement_functions_remain() -> None:
    assert len(cb.species.dry_air()) == cb.num_species()
    assert callable(cb.reforming_equilibrium)
    assert callable(cb.reforming_equilibrium_adiabatic)


@pytest.mark.parametrize("alias", ["ConstantLoss", "LinearThetaLoss"])
@pytest.mark.parametrize("module", [net, pressure_loss], ids=["network", "pressure_loss"])
def test_the_pressure_loss_aliases_are_gone(module: object, alias: str) -> None:
    with pytest.raises(AttributeError):
        getattr(module, alias)
    assert net.ConstantFractionLoss is pressure_loss.ConstantFractionLoss
    assert net.LinearThetaFractionLoss is pressure_loss.LinearThetaFractionLoss


def test_the_border_carnot_loss_element_is_gone() -> None:
    assert "BorderCarnotLossElement" not in net.__all__
    with pytest.raises(AttributeError):
        net.BorderCarnotLossElement  # noqa: B018
    with pytest.raises(ImportError):
        from combaero.network.components import BorderCarnotLossElement  # noqa: F401
    for module in (_core, _solver_tools):
        assert not hasattr(module, "border_carnot_loss_residual_and_jacobian")
    assert not hasattr(_core, "BorderCarnotLossResult")


@pytest.mark.parametrize(
    "kwargs",
    [{"omega_rpm": 3000.0}, {"r_c": 0.02}, {"n": 2.0}],
    ids=["omega_rpm", "r_c", "n"],
)
def test_the_removed_vortex_keywords_raise(kwargs: dict[str, float]) -> None:
    with pytest.raises(TypeError):
        VortexElement("v", "a", "b", **kwargs)


def test_the_vortex_omega_rpm_attribute_is_gone() -> None:
    el = VortexElement("v", "a", "b", shaft_speed_rpm=3000.0)
    assert el.shaft_speed_rpm == 3000.0
    with pytest.raises(AttributeError):
        el.omega_rpm  # noqa: B018


def test_the_vatistas_functions_remain() -> None:
    """Only the element's r_c/n went; the free functions are unchanged."""
    assert callable(cb.vatistas_v0_bar)


def test_incompressible_warmstart_raises_and_names_the_valid_strategies() -> None:
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=200000.0, Tt=300.0))
    g.add_node(PressureBoundary("out", Pt=100000.0, Tt=300.0))
    g.add_element(OrificeElement("orf", "in", "out", 0.6, 0.05, correlation="fixed"))
    with pytest.raises(ValueError, match="outletref_warmstart") as info:
        NetworkSolver(g).solve(init_strategy="incompressible_warmstart", timeout=10.0)
    assert "'incompressible_warmstart'" in str(info.value)  # the rejected value
    for valid in ("default", "analytical_pt_prop", "homotopy", "continuation"):
        assert valid in str(info.value)


def test_the_warmstart_maxfev_keyword_went_with_it() -> None:
    with pytest.raises(TypeError):
        NetworkSolver(FlowNetwork()).solve(warmstart_maxfev=150)  # type: ignore[call-arg]


def test_an_unknown_network_attribute_still_raises() -> None:
    with pytest.raises(AttributeError):
        net.NoSuchElement  # noqa: B018


def test_a_gui_file_saved_with_incompressible_warmstart_still_loads() -> None:
    """The library refuses the strategy; a saved GUI file (or a headless
    payload) gets its outlet-referenced successor, as the frontend maps it."""
    from gui.backend.schemas import SolverSettings

    s = SolverSettings(init_strategy="incompressible_warmstart")  # type: ignore[arg-type]
    assert s.init_strategy == "outletref_warmstart"


# ---------------------------------------------------------------------------
# Deprecated in 0.8.0, removed in 0.9.0: the stagnation spellings (#782).
#
# P/T are static and Pt/Tt total. Each old spelling (a keyword, a function
# name, a struct field) went in 0.9.0 and must fail; its replacement must
# work without a warning. The inventory below is the one combaero._renamed
# carried in 0.8.
# ---------------------------------------------------------------------------


_TP = {"T0": "Tt", "P0": "Pt"}
_EJ = {"p0": "Pt", "t0": "Tt"}

#: Keyword renames: function -> {old keyword: new keyword}.
KWARGS: dict[str, dict[str, str]] = {
    "nozzle_flow": _TP,
    "solve_A_eff_from_mdot": _TP,
    "solve_P_back_from_mdot": _TP,
    "solve_Pt_from_mdot": {"T0": "Tt"},
    "critical_pressure_ratio": _TP,
    "mach_from_pressure_ratio": _TP,
    "mass_flux_isentropic": _TP,
    "nozzle_cd": _TP,
    "nozzle_thrust_cd": _TP,
    "duct_static_state": _TP,
    "duct_htc_static_state": _TP,
    "duct_htc_static_state_dY": _TP,
    "sonic_static_pressure": {"T0": "Tt"},
    "orifice_compressible_mdot_and_jacobian": _TP,
    "T_from_stagnation": {"T0": "Tt"},
    "P_from_stagnation": _TP,
    "approach_static_pressure": {"P_total": "Pt", "T_total": "Tt"},
    "orifice_residuals_and_jacobian": {"T_up": "Tt_up"},
    "channel_residuals_and_jacobian": {"T_up": "Tt_up"},
    "orifice_compressible_residuals_and_jacobian": {"T_up": "Tt_up"},
    "heat_entropy_flux": {"T0_in": "Tt_in", "T0_out": "Tt_out"},
    "combustor_head_loss": {"T_in": "Tt_in", "P_in": "Pt_in"},
    "chamber_merge_face_state": {"T_main": "Tt_main", "T": "Tt"},
    "momentum_chamber_residual_and_jacobian": {"P_total": "Pt", "T": "Tt"},
    "orifice_tap_pressure_and_jacobian": {"P_total_up": "Pt_up", "P_total_down": "Pt_down"},
    "ejector_choked_mass_flow_and_jacobian": _EJ,
    "ejector_cd_nozzle_mass_flow_and_jacobian": _EJ,
    "ejector_cd_nozzle_exit_static_and_jacobian": _EJ,
}

#: Function renames: old name -> new name.
FUNCTIONS: dict[str, str] = {
    "T0_from_static": "Tt_from_static",
    "T0_from_static_v": "Tt_from_static_v",
    "P0_from_static": "Pt_from_static",
    "solve_P0_from_mdot": "solve_Pt_from_mdot",
    "T0_from_static_and_jacobian_M": "Tt_from_static_and_jacobian_M",
    "P0_from_static_and_jacobian_M": "Pt_from_static_and_jacobian_M",
}

#: Field renames: class -> {old field: new field}.
FIELDS: dict[str, dict[str, str]] = {
    "MixtureState": {"P_total": "Pt", "T_total": "Tt"},
    "MassStream": {"T": "Tt", "P_total": "Pt"},
    "StreamJacobian": {"d_T": "d_Tt", "d_P_total": "d_Pt"},
    "NozzleSolution": {"T0": "Tt", "P0": "Pt"},
    "DuctStaticState": {
        "dT_dT0": "dT_dTt",
        "dT_dP0": "dT_dPt",
        "dP_dT0": "dP_dTt",
        "dP_dP0": "dP_dPt",
        "drho_dT0": "drho_dTt",
        "drho_dP0": "drho_dPt",
        "dGeff_dT0": "dGeff_dTt",
        "dGeff_dP0": "dGeff_dPt",
    },
    "SonicStatic": {"dP_dT0": "dP_dTt"},
    "OrificeResult": {"d_mdot_dT_up": "d_mdot_dTt_up"},
    "ChannelSolverResult": {"d_dP_dT_up": "d_dP_dTt_up"},
    "CombustorLossResult": {"dxi_dT_in": "dxi_dTt_in", "dxi_dP_in": "dxi_dPt_in"},
    "MixerResult": {
        "T_mix": "Tt_mix",
        "P_total_mix": "Pt_mix",
        "dT_mix_d_stream": "dTt_mix_d_stream",
        "dP_total_mix_d_stream": "dPt_mix_d_stream",
        "dT_mix_d_delta_h": "dTt_mix_d_delta_h",
        "dT_mix_dQ": "dTt_mix_dQ",
    },
}


_X = list(cb.species.dry_air())
_Y = list(cb.mole_to_mass(_X))
_TT, _PT = 640.0, 3.1e5
_NOZZLE_GEOM = {
    "A_inlet": 3e-3,
    "A_throat": 1e-3,
    "A_exit": 2.5e-3,
    "x_throat": 0.05,
    "x_exit": 0.15,
}
_UP = {"m_dot": 0.4, "P_total_up": _PT, "P_static_up": 3.0e5, "Tt_up": _TT, "Y_up": _Y}
_EJ_TP = {"Pt": _PT, "Tt": _TT}
_EJ_GEOM = {
    "area_throat": 3.14e-5,
    "area_exit": 1.0e-4,
    "gamma": 1.4,
    "r_gas": 287.0,
    "eta": 0.95,
    "eps_frac": 0.01,
}


def _calls() -> dict[str, dict[str, object]]:
    """A valid all-keyword call, in the NEW names, for every renamed-kwarg function."""
    tp = {"Tt": _TT, "Pt": _PT}
    return {
        "nozzle_flow": {**tp, "P_back": 2.2e5, "A_eff": 1.3e-3, "X": _X},
        "solve_A_eff_from_mdot": {**tp, "P_back": 2.2e5, "mdot_target": 0.7, "X": _X},
        "solve_P_back_from_mdot": {**tp, "A_eff": 1.3e-3, "mdot_target": 0.5, "X": _X},
        "solve_Pt_from_mdot": {
            "Tt": _TT,
            "P_back": 2.2e5,
            "A_eff": 1.3e-3,
            "mdot_target": 0.7,
            "X": _X,
        },
        "critical_pressure_ratio": {**tp, "X": _X},
        "mach_from_pressure_ratio": {**tp, "P": 2.4e5, "X": _X},
        "mass_flux_isentropic": {**tp, "P": 2.4e5, "X": _X},
        "nozzle_cd": {**tp, "P_exit": 0.4e5, **_NOZZLE_GEOM, "X": _X},
        "nozzle_thrust_cd": {**tp, "P_design": 0.4e5, "P_amb": 1e5, **_NOZZLE_GEOM, "X": _X},
        "duct_static_state": {**tp, "G": 300.0, "X": _X},
        "duct_htc_static_state": {**tp, "G": 300.0, "X": _X},
        "duct_htc_static_state_dY": {
            **tp,
            "G": 300.0,
            "X": _X,
            "state": cb.duct_htc_static_state(_TT, _PT, 300.0, _X),
        },
        "sonic_static_pressure": {"Tt": _TT, "G": 700.0, "Y": _Y},
        "orifice_compressible_mdot_and_jacobian": {
            **tp,
            "P_back": 2.0e5,
            "X": _X,
            "Cd": 0.62,
            "area": 1e-3,
        },
        "T_from_stagnation": {"Tt": _TT, "M": 0.6, "X": _X},
        "P_from_stagnation": {**tp, "M": 0.6, "X": _X},
        "approach_static_pressure": {**tp, "G": 250.0, "Y": _Y, "incompressible": False},
        "orifice_residuals_and_jacobian": {**_UP, "P_static_down": 2.5e5, "Cd": 0.62, "area": 1e-3},
        "orifice_compressible_residuals_and_jacobian": {
            **_UP,
            "P_static_down": 2.5e5,
            "Cd": 0.62,
            "area": 1e-3,
        },
        "channel_residuals_and_jacobian": {
            **_UP,
            "P_static_down": 2.5e5,
            "L": 1.2,
            "D": 0.03,
            "roughness": 1e-5,
            "friction_model": "haaland",
        },
        "heat_entropy_flux": {
            "m_dot": 0.8,
            "Tt_in": 600.0,
            "Tt_out": 900.0,
            "Y": _Y,
            "M_in": 0.1,
            "M_out": 0.3,
        },
        "combustor_head_loss": {
            "theta": 0.2,
            "mdot": 0.4,
            "Tt_in": 700.0,
            "Pt_in": 2.0e5,
            "Y_in": _Y,
            "k": 1.0,
            "zeta0": 3.0,
            "area": 0.05,
        },
        "chamber_merge_face_state": {
            "m_main": 0.2,
            "Tt_main": 700.0,
            "X_main": _X,
            "m_out": 0.23,
            "P": 1.5e5,
            "Tt": 650.0,
            "X": _X,
            "side_momentum": 8.0,
            "area": 3.14e-4,
        },
        "momentum_chamber_residual_and_jacobian": {
            "P": 2.5e5,
            "Pt": 2.6e5,
            "m_dot": 0.7,
            "Tt": _TT,
            "Y": _Y,
            "area": 4e-3,
        },
        "orifice_tap_pressure_and_jacobian": {
            "P_static_up": 3.0e5,
            "Pt_up": _PT,
            "Pt_down": 2.9e5,
            "beta": 0.5,
            "C": 0.61,
            "compressible": True,
        },
        "ejector_choked_mass_flow_and_jacobian": {
            **_EJ_TP,
            "area_throat": 3.14e-5,
            "gamma": 1.4,
            "r_gas": 287.0,
            "eta": 0.95,
        },
        "ejector_cd_nozzle_mass_flow_and_jacobian": {
            **_EJ_TP,
            "p_static_down": 2.0e5,
            **_EJ_GEOM,
        },
        "ejector_cd_nozzle_exit_static_and_jacobian": {
            "m_dot": 0.005,
            **_EJ_TP,
            **_EJ_GEOM,
        },
    }


def _function(name: str) -> object:
    return getattr(cb, name, None) or getattr(_solver_tools, name)


def test_the_782_inventory_is_covered() -> None:
    assert set(_calls()) == set(KWARGS)


@pytest.mark.parametrize("fn_name", sorted(KWARGS))
def test_an_old_keyword_is_refused(fn_name: str) -> None:
    fn = _function(fn_name)
    calls = _calls()[fn_name]
    for old, new in KWARGS[fn_name].items():
        old_kwargs = {(old if k == new else k): v for k, v in calls.items()}
        with pytest.raises(TypeError):
            fn(**old_kwargs)


@pytest.mark.parametrize("fn_name", sorted(KWARGS))
def test_the_new_keywords_work_without_a_warning(fn_name: str) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        _function(fn_name)(**_calls()[fn_name])


def test_no_function_is_wrapped_any_more() -> None:
    for name in KWARGS:
        core_fn = getattr(_core, name, None)
        if core_fn is not None:
            assert not hasattr(core_fn, "__renamed_kwargs__")
    with pytest.raises(ImportError):
        from combaero import _renamed  # noqa: F401


_FUNCTION_ARGS: dict[str, tuple[object, ...]] = {
    "Tt_from_static": (560.0, 0.6, _X),
    "Tt_from_static_v": (560.0, 280.0, _X),
    "Pt_from_static": (2.5e5, 560.0, 0.6, _X),
    "solve_Pt_from_mdot": (_TT, 2.2e5, 1.3e-3, 0.7, _X),
    "Tt_from_static_and_jacobian_M": (560.0, 0.6, _X),
    "Pt_from_static_and_jacobian_M": (2.5e5, 560.0, 0.6, _X),
}


@pytest.mark.parametrize(("old", "new"), sorted(FUNCTIONS.items()))
def test_an_old_function_name_is_gone(old: str, new: str) -> None:
    assert old not in cb.__all__
    for module in (cb, _core, _solver_tools):
        assert not hasattr(module, old)
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        _function(new)(*_FUNCTION_ARGS[new])


def test_the_new_function_names_stay_exported() -> None:
    for new in ("Tt_from_static", "Tt_from_static_v", "Pt_from_static", "solve_Pt_from_mdot"):
        assert new in cb.__all__


def _instances() -> dict[str, object]:
    s1 = cb.MassStream(0.8, 600.0, 3.0e5, _Y)
    s2 = cb.MassStream(0.5, 400.0, 3.1e5, _Y)
    mix = _solver_tools.mixer_from_streams_and_jacobians([s1, s2], 1500.0, 0.0)
    calls = _calls()
    return {
        "MixerResult": mix,
        "MixtureState": cb.MixtureState(2.5e5, 3.0e5, 560.0, 640.0, 0.7, _Y),
        "MassStream": s1,
        "StreamJacobian": mix.dTt_mix_d_stream[0],
        "NozzleSolution": cb.nozzle_cd(**calls["nozzle_cd"]),
        "DuctStaticState": cb.duct_static_state(_TT, _PT, 300.0, _X),
        "SonicStatic": _solver_tools.sonic_static_pressure(_TT, 700.0, _Y),
        "OrificeResult": _solver_tools.orifice_residuals_and_jacobian(
            **calls["orifice_residuals_and_jacobian"]
        ),
        "ChannelSolverResult": _solver_tools.channel_residuals_and_jacobian(
            **calls["channel_residuals_and_jacobian"]
        ),
        "CombustorLossResult": cb.combustor_head_loss(**calls["combustor_head_loss"]),
    }


def test_the_782_field_inventory_is_covered() -> None:
    assert set(_instances()) == set(FIELDS)


_FIELD_CASES = [
    (cls, old, new)
    for cls, mapping in sorted(FIELDS.items())
    for old, new in sorted(mapping.items())
]


@pytest.mark.parametrize(
    ("cls", "old", "new"), _FIELD_CASES, ids=[f"{c}.{o}" for c, o, _ in _FIELD_CASES]
)
def test_an_old_field_is_gone(cls: str, old: str, new: str) -> None:
    obj = _instances()[cls]
    assert type(obj).__name__ == cls
    with pytest.raises(AttributeError):
        getattr(obj, old)
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        getattr(obj, new)


def test_mixture_state_old_constructor_keywords_are_refused() -> None:
    cb.MixtureState(P=2.5e5, Pt=3.0e5, T=560.0, Tt=640.0, m_dot=0.7, Y=_Y)
    with pytest.raises(TypeError):
        cb.MixtureState(P=2.5e5, P_total=3.0e5, T=560.0, T_total=640.0, m_dot=0.7, Y=_Y)


def test_compressible_nozzle_flow_old_keywords_are_refused() -> None:
    from combaero import compressible

    compressible.nozzle_flow(Tt=_TT, Pt=_PT, X=_X, P_back=2.2e5, A_eff=1.3e-3)
    with pytest.raises(TypeError):
        compressible.nozzle_flow(T=_TT, P=_PT, X=_X, P_back=2.2e5, A_eff=1.3e-3)


# ---------------------------------------------------------------------------
# Deprecated in 0.8.0, removed in 0.9.0: Al 6061 (#761) and Stolz (#876).
# ---------------------------------------------------------------------------


def test_aluminum_6061_is_gone() -> None:
    assert "k_aluminum_6061" not in cb.__all__
    assert not hasattr(cb, "k_aluminum_6061")
    assert not hasattr(_core, "k_aluminum_6061")
    assert "al6061" not in cb.list_materials()
    with pytest.raises(ValueError):
        cb.get_material_conductivity("al6061", 300.0)


def test_stolz_is_gone_from_the_metering_selector() -> None:
    assert not hasattr(cb.MeteringCdCorrelation, "Stolz")
    assert hasattr(cb.MeteringCdCorrelation, "ReaderHarrisGallagher")
    for internal in ("_metering_cd_and_dRe", "_metering_cd_in_range"):
        assert not hasattr(_core, internal)


def test_a_stolz_orifice_element_is_refused_naming_the_replacement() -> None:
    with pytest.raises(ValueError, match=r"(?s)Stolz.*removed in 0\.9\.0.*ReaderHarrisGallagher"):
        OrificeElement("o", "a", "b", diameter=0.05, correlation="Stolz")


def test_a_saved_gui_orifice_with_stolz_is_refused_naming_the_replacement() -> None:
    """A GUI file saved with the Stolz Cd fails to load with a message that
    names Reader-Harris/Gallagher; it is not silently migrated (#876)."""
    from pydantic import ValidationError

    from gui.backend.schemas import OrificeData

    # A pre-hole_cd file (the loose field) and a current one alike.
    for saved in ({"correlation": "Stolz"}, {"hole_cd": {"correlation": "Stolz"}}):
        with pytest.raises(ValidationError, match=r"(?s)Stolz.*0\.9\.0.*ReaderHarrisGallagher"):
            OrificeData(**saved)
    kept = OrificeData(correlation="ReaderHarrisGallagher")
    assert kept.hole_cd.correlation == "ReaderHarrisGallagher"


def test_the_gui_solve_response_has_no_deprecations_field() -> None:
    from gui.backend.schemas import SolveResponse

    assert "deprecations" not in SolveResponse.model_fields


def test_a_gui_solve_of_a_saved_stolz_orifice_is_a_400_naming_the_replacement() -> None:
    import asyncio

    from fastapi import HTTPException

    from gui.backend.main import solve
    from gui.backend.schemas import NetworkGraphSchema

    def node(i: str, t: str, data: dict, x: float) -> dict:
        return {"id": i, "type": t, "position": {"x": x, "y": 0}, "data": data}

    def edge(i: str, s: str, t: str) -> dict:
        return {
            "id": i,
            "source": s,
            "sourceHandle": "flow-source",
            "target": t,
            "targetHandle": "flow-target",
            "data": {},
        }

    graph = {
        "nodes": [
            node("pb_in", "pressure_boundary", {"Pt": 2e5, "Tt": 300}, 0),
            node("ori", "orifice", {"diameter": 0.01, "correlation": "Stolz"}, 100),
            node("pb_out", "pressure_boundary", {"Pt": 1e5, "Tt": 300}, 200),
        ],
        "edges": [edge("e1", "pb_in", "ori"), edge("e2", "ori", "pb_out")],
    }
    with pytest.raises(HTTPException) as info:
        asyncio.run(solve(NetworkGraphSchema(**graph)))
    assert info.value.status_code == 400
    assert "Stolz" in info.value.detail and "ReaderHarrisGallagher" in info.value.detail


# --- Deprecated in 0.12.0, removed in 0.13.0: the per-element hole-Cd keywords.
# HoleCd (combaero.network) carries the correlation, the fixed Cd and each
# term's model; the keywords it replaces still build the same element (bit for
# bit, test_hole_cd_spec.py), with ONE DeprecationWarning per construction
# naming every keyword given, its HoleCd field and the removal release.

_D_JET = 0.00254

_HOLE_ELEMENTS = {
    "EffusionPlateElement": lambda **kw: net.EffusionPlateElement(
        "h", "a", "b", hole_diameter=1e-3, wall_thickness=2e-3, pitch=6e-3, panel_area=0.01,
        angle_deg=30.0, **kw,
    ),
    "FilmCoolingElement": lambda **kw: net.FilmCoolingElement(
        "h", "a", "b", hole_diameter=1e-3, holes_per_row=10, lateral_pitch=3e-3,
        wall_thickness=3e-3, **kw,
    ),
    "SlotFilmElement": lambda **kw: net.SlotFilmElement(
        "h", "a", "b", slot_height=1e-3, span=0.05, lip_thickness=4e-3, **kw
    ),
    "ImpingementPlateElement": lambda **kw: net.ImpingementPlateElement(
        "h", "a", "b", d_jet=_D_JET, xn_d=5.0, yn_d=4.0, z_d=1.0, span=48 * _D_JET,
        plate_thickness=_D_JET, **kw,
    ),
    "ImpingementArray": lambda **kw: net.ImpingementArray(
        "h", n_rows=3, d_jet=_D_JET, xn_d=5.0, yn_d=4.0, z_d=1.0, span=48 * _D_JET,
        plate_thickness=_D_JET, **kw,
    ),
}  # fmt: skip

# Each keyword, as given, and the HoleCd field the warning names for it.
_KW = {
    "correlation": ({"correlation": "McGreehanSchotsch"}, "correlation='McGreehanSchotsch'"),
    "Cd": ({"correlation": "fixed", "Cd": 0.7}, "Cd=0.7"),
    "cd_crossflow": ({"cd_crossflow": "none"}, "exit_crossflow='none'"),
}
_COOLING = ("correlation", "Cd")
_JET = ("correlation", "Cd", "cd_crossflow")
_HOLE_CD_KEYWORDS = [
    (name, *_KW[k])
    for name, kws in (
        ("EffusionPlateElement", _COOLING),
        ("FilmCoolingElement", _COOLING),
        ("SlotFilmElement", ("correlation", "Cd")),
        ("ImpingementPlateElement", _JET),
        ("ImpingementArray", _JET),
    )
    for k in kws
]


@pytest.mark.parametrize(("name", "kw", "field"), _HOLE_CD_KEYWORDS)
def test_a_deprecated_hole_cd_keyword_warns_once_and_names_its_replacement(
    name: str, kw: dict, field: str
) -> None:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        el = _HOLE_ELEMENTS[name](**kw)
    dep = [w for w in caught if issubclass(w.category, DeprecationWarning)]
    assert len(dep) == 1, [str(w.message) for w in dep]
    msg = str(dep[0].message)
    assert all(k in msg for k in kw), msg
    assert "hole_cd=HoleCd(" in msg and field in msg, msg
    assert "deprecated since 0.12.0" in msg and "removed in 0.13.0" in msg, msg
    # Attributed to the caller's line, not to the package.
    assert dep[0].filename == __file__
    # The same element as its HoleCd builds, without a warning.
    fields = {("exit_crossflow" if k == "cd_crossflow" else k): v for k, v in kw.items()}
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        new = _HOLE_ELEMENTS[name](hole_cd=net.HoleCd(**fields))
    if isinstance(el, net.ImpingementArray):  # its fields keep what was passed
        assert el._rows_hole_cd == new._rows_hole_cd
    else:
        assert el.hole_cd == new.hole_cd


def test_an_array_built_from_old_fields_warns_once_not_per_row() -> None:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        arr = _HOLE_ELEMENTS["ImpingementArray"](correlation="fixed", Cd=0.8, cd_crossflow="none")
        g = FlowNetwork()
        g.add_node(PressureBoundary("sup", Pt=1.03e5, Tt=300.0))
        g.add_node(PressureBoundary("ex", Pt=1.0e5, Tt=300.0))
        arr.add_to(g, "sup", "ex")
    dep = [w for w in caught if issubclass(w.category, DeprecationWarning)]
    assert len(dep) == 1 and "correlation, Cd, cd_crossflow" in str(dep[0].message)
    # The fields keep what was passed; dataclasses.replace rebuilds from them.
    assert (arr.correlation, arr.Cd, arr.cd_crossflow, arr.hole_cd) == ("fixed", 0.8, "none", None)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        again = dataclasses.replace(arr, Cd=0.7)
    assert again._rows_hole_cd == net.HoleCd(correlation="fixed", Cd=0.7, exit_crossflow="none")
    assert all(g.elements[f"h__p{i}"].cd_crossflow == "none" for i in (1, 2, 3))


def test_the_plain_orifice_keeps_its_correlation_and_cd_keywords() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        OrificeElement("o", "a", "b", Cd=0.7, diameter=0.01, correlation="fixed")
        OrificeElement("o", "a", "b", diameter=0.01, correlation="IdelchikThick")
