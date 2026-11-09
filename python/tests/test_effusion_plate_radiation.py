"""EffusionPlateElement's opt-in gas radiation on the hot face (#1071).

`radiation=GasRadiation(...)` is the ThermalWall face spec (#793): Lefebvre &
Ballal (2010) Eq. 8.7 with Eq. 8.11's emissivity, solved with the wall in
C++'s `wall_radiation_sides`. Checked here against the printed equations
written out by hand and against the wall's own energy balance; off, the
plate reports what it always did.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network import (
    CasingRadiation,
    ChannelElement,
    EffusionPlateElement,
    FlowNetwork,
    GasRadiation,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
)

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
T_GAS = 1600.0
GEOM = {
    "hole_diameter": 0.8e-3,
    "wall_thickness": 2e-3,
    "pitch": 6e-3,
    "panel_length": 0.1,
    "panel_width": 0.1,
    "angle_deg": 30.0,
    "wall_conductivity": 20.0,
}


def _products() -> list[float]:
    """Lean methane/air complete-combustion products: H2O and CO2 to radiate."""
    x = list(cb.species.dry_air())
    i = cb.species.species_index("CH4")
    n_f = 0.5 * 0.2095 / 2.0  # phi 0.5: CH4 + 2 O2
    x = [v / (1.0 + n_f) for v in x]
    x[i] += n_f / (1.0 + n_f)
    return list(cb.mole_to_mass(list(cb.complete_combustion(700.0, x, 1.0e6).X)))


Y_GAS = _products()


def _pb(name: str, Pt: float, Tt: float, Y: list[float]) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y
    return b


def _diag(chamber: bool = True, **kw) -> dict:
    g = FlowNetwork()
    g.add_node(_pb("cool", 10.3e5, 700.0, Y_AIR))
    if chamber:
        g.add_node(_pb("gas_in", 10.03e5, T_GAS, Y_GAS))
        g.add_node(_pb("exit", 10.0e5, T_GAS, Y_GAS))
        g.add_node(MomentumChamberNode("ch", main_inlet="duct"))
        g.add_element(
            ChannelElement("duct", "gas_in", "ch", length=0.2, diameter=0.08, roughness=0.0)
        )
        g.add_element(
            ChannelElement("tail", "ch", "exit", length=0.2, diameter=0.08, roughness=0.0)
        )
    else:
        g.add_node(_pb("ch", 10.0e5, T_GAS, Y_GAS))
    g.add_element(EffusionPlateElement("eff", "cool", "ch", **{**GEOM, **kw}))
    r = NetworkSolver(g).solve()
    assert r["__success__"], r.get("__message__")
    return r["__element_diag__"]["eff"]


@pytest.fixture(scope="module")
def off() -> dict:
    return _diag()


@pytest.fixture(scope="module")
def on() -> dict:
    return _diag(radiation=GasRadiation(beam_length=0.05, emissivity_wall=0.8))


def test_off_reports_no_radiation(off) -> None:
    assert "q_rad_hot" not in off and "eps_g" not in off


def test_the_flux_is_eq_8_7_with_eq_8_11_written_out(on) -> None:
    Tg, Tw = on["T_gas"], on["T_wall_hot"]
    assert on["T_gas_from_main_inlet"] == 1.0 and Tg == pytest.approx(T_GAS, rel=1e-12)
    q = cb.elemental_fuel_air_ratio(Y_GAS)[0]
    assert 0.02 < q < 0.04
    P_kPa = 10.0e5 * 1e-3  # the chamber's static P, within its tiny dynamic head
    eps_g = 1.0 - math.exp(-290.0 * P_kPa * 1.0 * (q * 0.05) ** 0.5 * Tg**-1.5)
    assert on["eps_g"] == pytest.approx(eps_g, rel=2e-3)
    r1 = 0.5 * cb.STEFAN_BOLTZMANN * 1.8 * on["eps_g"] * Tg**1.5 * (Tg**2.5 - Tw**2.5)
    assert on["q_rad_hot"] == pytest.approx(r1, rel=1e-9)
    assert on["q_rad_hot"] > 1.0e4


def test_the_wall_balances_convection_radiation_and_conduction(on) -> None:
    q = on["q_wall"]
    gas = on["h_gas"] * (on["T_adiabatic_wall"] - on["T_wall_hot"]) + on["q_rad_hot"]
    cond = (on["T_wall_hot"] - on["T_wall_cold"]) / (GEOM["wall_thickness"] / 20.0)
    assert gas == pytest.approx(q, rel=1e-9)
    assert cond == pytest.approx(q, rel=1e-9)
    assert on["Q_wall"] == pytest.approx(q * 0.1 * 0.1, rel=1e-12)


def test_radiation_heats_the_wall_and_lowers_eta(off, on) -> None:
    assert on["T_wall_hot"] > off["T_wall_hot"] + 5.0
    assert on["q_wall"] > off["q_wall"]
    assert on["eta_overall"] < off["eta_overall"]
    # The convective side and the flow are untouched: it is the wall's output.
    for key in ("h_gas", "T_gas", "T_adiabatic_wall", "blowing_ratio", "h_internal_plate_area"):
        assert on[key] == off[key], key


def test_a_longer_beam_radiates_more() -> None:
    """Falsifying the metric: the emissivity, hence the flux, must follow l_b."""
    short = _diag(radiation=GasRadiation(beam_length=0.01))
    long = _diag(radiation=GasRadiation(beam_length=0.2))
    assert long["eps_g"] > short["eps_g"] and long["q_rad_hot"] > short["q_rad_hot"]
    assert long["T_wall_hot"] > short["T_wall_hot"]


def test_a_fixed_spec_overrides_the_stream() -> None:
    d = _diag(radiation=GasRadiation(beam_length=0.05, fuel_air_ratio=0.0))
    assert d["eps_g"] == 0.0 and d["q_rad_hot"] == 0.0


def test_a_plenum_keeps_its_imposed_flux() -> None:
    d = _diag(chamber=False, gas_heat_flux=1.0e4, radiation=GasRadiation(beam_length=0.05))
    assert d["q_wall"] == 1.0e4
    assert math.isnan(d["q_rad_hot"]) and math.isnan(d["eps_g"])


def test_only_gas_radiation_faces_the_hot_side() -> None:
    with pytest.raises(ValueError, match="GasRadiation"):
        EffusionPlateElement("e", "a", "b", **GEOM, radiation=CasingRadiation(area_ratio=0.8))
