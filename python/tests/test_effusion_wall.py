"""EffusionPlateElement owns its wall; the discharge node decides the gas side (#471).

* Into a MomentumChamberNode (flow): the chamber's own surface correlation
  is the unblown gas coefficient, its velocity sets the blowing ratio.
* Into a plenum (state, no flow): an imposed heat flux.

The wall is the plate itself, solved through C++; its heat goes from the gas
into the effusing coolant and back into the same node, so it is an output.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    EffusionPlateElement,
    FlowNetwork,
    MomentumChamberNode,
    NetworkMixtureState,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)

Y_AIR = cb.mole_to_mass(cb.species.dry_air())


def _pb(name: str, Pt: float, Tt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y_AIR
    return b


def _plate(**kw) -> EffusionPlateElement:
    args = {
        "hole_diameter": 0.8e-3,
        "wall_thickness": 2e-3,
        "pitch": 6e-3,
        "panel_length": 0.1,
        "panel_width": 0.1,
        "angle_deg": 30.0,
    }
    args.update(kw)
    return EffusionPlateElement("eff", "cool", "ch", **args)


def _liner(**kw) -> tuple[NetworkSolver, dict]:
    g = FlowNetwork()
    g.add_node(_pb("gas_in", 1.003e5, 1500.0))
    g.add_node(_pb("exit", 1.0e5, 1500.0))
    g.add_node(_pb("cool", 1.03e5, 600.0))
    g.add_node(MomentumChamberNode("ch", main_inlet="duct"))
    g.add_element(ChannelElement("duct", "gas_in", "ch", length=0.2, diameter=0.08, roughness=0.0))
    g.add_element(ChannelElement("tail", "ch", "exit", length=0.2, diameter=0.08, roughness=0.0))
    g.add_element(_plate(**kw))
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    return s, r


def _plenum(**kw) -> dict:
    g = FlowNetwork()
    g.add_node(_pb("cool", 1.03e5, 600.0))
    g.add_node(_pb("ch", 1.0e5, 1500.0))
    g.add_element(_plate(**kw))
    r = NetworkSolver(g).solve()
    assert r["__success__"], r.get("__message__")
    return r["__element_diag__"]["eff"]


def test_a_chamber_discharge_takes_its_gas_side_from_the_chamber() -> None:
    s, r = _liner()
    d = r["__element_diag__"]["eff"]
    ch = s.network.nodes["ch"]
    T, P = s._derived_states["ch"][0], r["ch.P"]
    state = NetworkMixtureState(P=P, Pt=r["ch.Pt"], Tt=T, m_dot=0.0, Y=Y_AIR)
    # The plate owns the wall; the chamber supplies only its h per unit area,
    # with or without a wetted area of its own (#1022: this chamber has none,
    # so its own htc_and_T, gated on surface.area, returns None).
    assert ch.surface.area == 0.0
    gas = ch._gas_side_htc_and_T(state)
    assert d["gas_side_chamber"] == 1.0
    assert d["h_gas_unblown"] == pytest.approx(gas.h, rel=1e-6)
    assert 600.0 < d["T_wall_cold"] < d["T_wall_hot"] < d["T_gas"]
    # What the augmentation multiplies is reported with it (#476).
    model = ch.surface.model
    assert d["h_gas_unblown_model"].startswith(type(model).__name__)
    assert d["h_gas_unblown_Re"] == pytest.approx(gas.Re, rel=1e-6)
    assert d["h_gas_unblown_Pr"] == pytest.approx(gas.Pr, rel=1e-6)
    assert d["h_gas_unblown_Nu"] == pytest.approx(gas.Nu, rel=1e-6)


def test_the_gas_temperature_is_the_approaching_main_stream() -> None:
    """A merge chamber's state is the MIXED outlet, diluted by the effused
    coolant; the wall sees the approaching main stream."""
    s, r = _liner()
    d = r["__element_diag__"]["eff"]
    assert d["T_gas_from_main_inlet"] == 1.0
    assert d["T_gas"] == pytest.approx(1500.0)
    assert s._derived_states["ch"][0] < d["T_gas"] - 5.0


def test_the_wall_is_the_validated_closure_when_conduction_vanishes() -> None:
    """With no conduction resistance the wall reduces to overall_effectiveness,
    the two-resistance closure scored on Andrews 88-GT-290."""
    s, r = _liner(wall_conductivity=1e9, gas_augmentation=1.7)
    d = r["__element_diag__"]["eff"]
    eff = s.network.elements["eff"]
    supply = NetworkMixtureState(P=1.03e5, Pt=1.03e5, Tt=600.0, m_dot=d["m_dot"], Y=Y_AIR)
    ref = eff.overall_effectiveness(supply, d["h_gas_unblown"], d["T_gas"], gas_augmentation=1.7)
    assert d["eta_overall"] == pytest.approx(ref["eta_overall"], rel=1e-6)
    assert d["T_wall_hot"] == pytest.approx(ref["T_wall"], rel=1e-6)


def test_conduction_opens_a_temperature_drop_across_the_wall() -> None:
    _, r = _liner(wall_conductivity=20.0)
    d = r["__element_diag__"]["eff"]
    assert d["T_wall_hot"] - d["T_wall_cold"] == pytest.approx(d["q_wall"] * 2e-3 / 20.0, rel=1e-6)


def test_augmentation_heats_the_wall() -> None:
    """The caller's knob, never fitted: more gas-side h, hotter wall."""
    cold = _liner(gas_augmentation=1.0)[1]["__element_diag__"]["eff"]["T_wall_hot"]
    hot = _liner(gas_augmentation=2.0)[1]["__element_diag__"]["eff"]["T_wall_hot"]
    assert hot > cold


def test_the_film_is_off_by_default_and_reported_when_selected() -> None:
    off = _liner()[1]["__element_diag__"]["eff"]
    on = _liner(gas_film="baldauf_sellers")[1]["__element_diag__"]["eff"]
    assert off["eta_film"] == 0.0
    assert 0.0 < on["eta_film"] < 1.0
    assert on["T_wall_hot"] < off["T_wall_hot"]
    # 6 mm pitch on 0.8 mm holes: s/D 7.5 is outside Baldauf's 2-5.
    assert on["film_extrapolated"] == 1.0


def test_the_film_matches_the_validation_runners_panel_average() -> None:
    """One definition: the element and effusion_overall_runner share C++'s
    effusion_panel_film_effectiveness."""
    _, r = _liner(gas_film="baldauf_sellers")
    d = r["__element_diag__"]["eff"]
    rows = round(0.1 / 6e-3)
    film = cb.effusion_panel_film_effectiveness(
        rows, 6e-3 / 0.8e-3, 6e-3 / 0.8e-3, d["blowing_ratio"], d["density_ratio"], 30.0, 0.05
    )
    assert d["eta_film"] == pytest.approx(film.eta, rel=1e-12)


def test_a_plenum_discharge_takes_the_imposed_heat_flux() -> None:
    d = _plenum(gas_heat_flux=2.0e5, wall_conductivity=20.0)
    assert d["gas_side_chamber"] == 0.0
    assert d["q_wall"] == 2.0e5
    assert d["Q_wall"] == pytest.approx(2.0e5 * 0.01)
    assert d["T_wall_cold"] == pytest.approx(600.0 + 2.0e5 / d["h_internal_plate_area"], rel=1e-9)
    assert d["T_wall_hot"] - d["T_wall_cold"] == pytest.approx(2.0e5 * 2e-3 / 20.0, rel=1e-9)


def test_a_plenum_discharge_defaults_to_an_adiabatic_gas_side() -> None:
    d = _plenum()
    assert d["q_wall"] == 0.0
    assert d["T_wall_hot"] == pytest.approx(600.0)


def test_a_channel_fed_supply_is_flagged() -> None:
    """The 2-port plate is plenum-fed; a supply node a channel runs through
    feeds it a crossflow it does not model (U1/Vi stays 0)."""
    g = FlowNetwork()
    g.add_node(_pb("cool_in", 1.04e5, 600.0))
    g.add_node(_pb("cool_out", 1.035e5, 600.0))
    g.add_node(PlenumNode("cool"))
    g.add_node(_pb("ch", 1.0e5, 1500.0))
    g.add_element(ChannelElement("c1", "cool_in", "cool", length=0.1, diameter=0.02))
    g.add_element(ChannelElement("c2", "cool", "cool_out", length=0.1, diameter=0.02))
    g.add_element(_plate())
    r = NetworkSolver(g).solve()
    assert r["__success__"], r.get("__message__")
    assert r["__element_diag__"]["eff"]["Cd_inlet_crossflow_ignored"] == 1.0
    assert _plenum()["Cd_inlet_crossflow_ignored"] == 0.0


def test_film_selection_needs_the_panel_length() -> None:
    with pytest.raises(ValueError, match="panel_length"):
        EffusionPlateElement(
            "eff",
            "a",
            "b",
            hole_diameter=1e-3,
            wall_thickness=2e-3,
            pitch=6e-3,
            panel_area=0.01,
            gas_film="baldauf_sellers",
        )


def test_the_energy_balance_is_untouched() -> None:
    """The wall heat leaves the gas and returns with the coolant into the same
    node: the chamber outlet temperature is the adiabatic mix."""
    s, r = _liner()
    d = r["__element_diag__"]
    m_g, m_c = d["duct"]["m_dot"], d["eff"]["m_dot"]
    X = cb.species.dry_air()
    h_mix = (m_g * cb.h_mass(1500.0, X) + m_c * cb.h_mass(600.0, X)) / (m_g + m_c)
    T_ch = s._derived_states["ch"][0]
    assert cb.h_mass(T_ch, X) == pytest.approx(h_mix, rel=1e-6)
    assert math.isfinite(d["eff"]["Q_wall"]) and d["eff"]["Q_wall"] > 0.0
