"""Hole-Cd envelope flags as INFORMATION, by variable (hole-Cd audit stage 4).

Extrapolation is inevitable and the user can do nothing better about it, so
the results state which input of a hole's Cd lies outside its source's range
(``Cd_extrapolated_Re`` / ``_l_d`` / ``_edge``, and Rogers & Hersh's
``Cd_crossflow_extrapolated_*``), never a warning. The flags change no Cd:
the values are bit-identical to before (checked over the stress set by F, J
hashes, and here on the C++ calls directly).
"""

from __future__ import annotations

import math
import warnings

import pytest

import combaero as cb
from combaero.network import (
    EffusionPlateElement,
    FilmCoolingElement,
    FlowNetwork,
    NetworkSolver,
    PressureBoundary,
)

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
C = cb.DischargeCdCorrelation
KEYS = ("Cd_extrapolated_Re", "Cd_extrapolated_l_d", "Cd_extrapolated_edge")


def _solve(element, p_in: float = 1.3e5) -> dict:
    g = FlowNetwork()
    g.add_node(PressureBoundary("cool", Pt=p_in, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("gas", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(element)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(g).solve()
    assert r["__success__"], r.get("__message__")
    return r["__element_diag__"][element.id]


def _film(l_over_d: float, correlation: str, d: float = 2e-3) -> FilmCoolingElement:
    # Normal holes: the drilled length is the wall thickness.
    return FilmCoolingElement(
        "film", "cool", "gas", hole_diameter=d, holes_per_row=4, lateral_pitch=4 * d,
        wall_thickness=l_over_d * d, angle_deg=90.0, hole_type="cylindrical",
        correlation=correlation,
    )  # fmt: skip


def test_a_thick_hole_beyond_its_table_names_l_d_only() -> None:
    # Idelchik 4-18a's table ends at l/d 4; beyond it zeta' is held and only
    # bore friction grows (a gentle fall). The flag names l/d, not Re.
    d = _solve(_film(6.0, "IdelchikThick"))
    assert d["Cd_extrapolated_l_d"] == 1.0
    assert d["Cd_extrapolated_Re"] == 0.0 and d["Cd_extrapolated_edge"] == 0.0
    assert d["Cd_in_range"] == 0.0
    inside = _solve(_film(3.0, "IdelchikThick"))
    assert all(inside[k] == 0.0 for k in KEYS) and inside["Cd_in_range"] == 1.0


def test_lichtarowicz_at_its_high_re_constant_is_in_range() -> None:
    # Above Eq. (12)'s fitted Re 2e4 the source takes Cd constant at Eq. (7)
    # (design recommendation (2)), in range to the end of its data (~1e5).
    d = _solve(_film(4.0, "Lichtarowicz", d=5e-3), p_in=1.6e5)
    re = 4.0 * abs(d["m_dot"]) / 4 / (math.pi * 5e-3 * 1.85e-5)
    assert 2.0e4 < re < 1.0e5
    assert d["Cd_in_range"] == 1.0
    assert all(d[k] == 0.0 for k in KEYS)


def test_mcgreehan_schotsch_beyond_l_d_10_is_flagged_by_variable() -> None:
    # Eq. (14)'s g keeps falling beyond the Lichtarowicz data it was fitted
    # to (Fig. 3, L/d to 10): Cd 0.49 at l/d 40, the steepest extrapolation
    # in the family. Reported, not refused, and the Cd is unchanged.
    d = _solve(_film(20.0, "McGreehanSchotsch", d=5e-3), p_in=1.6e5)
    assert d["Cd_extrapolated_l_d"] == 1.0
    assert d["Cd_extrapolated_Re"] == 0.0
    assert d["Cd_in_range"] == 0.0


def test_mcgreehan_schotsch_below_its_re_floor_is_info_not_stderr(capfd) -> None:
    # The held value below Re 1e4 used to print a C++ warning on every Cd
    # call (thousands per solve). Now the result says it, once, by name.
    capfd.readouterr()
    plate = EffusionPlateElement(
        "eff", "cool", "gas", hole_diameter=0.6e-3, wall_thickness=1e-3, pitch=3e-3,
        panel_length=0.02, panel_width=0.02, correlation="McGreehanSchotsch",
    )  # fmt: skip
    d = _solve(plate, p_in=1.02e5)
    assert d["Cd_extrapolated_Re"] == 1.0
    assert "below the correlation's floor" not in capfd.readouterr().err


def test_the_flags_are_the_c_plus_plus_envelope() -> None:
    h = cb.DischargeHoleGeometry(d=1e-3, L=40e-3)
    e = cb.discharge_cd_envelope(C.McGreehanSchotsch1988, h, cb.DischargeHoleState(Re=5e4))
    assert (e.re, e.l_over_d, e.edge, e.in_range()) == (False, True, False, False)
    assert (
        cb.discharge_cd_in_range(C.McGreehanSchotsch1988, h, cb.DischargeHoleState(Re=5e4)) is False
    )
    h = cb.DischargeHoleGeometry(d=1e-3, L=4e-3)
    e = cb.discharge_cd_envelope(C.Lichtarowicz1965, h, cb.DischargeHoleState(Re=5e4))
    assert e.in_range()
    # Beyond the collected data (about 1e5) the held Eq. (7) is extrapolated.
    e = cb.discharge_cd_envelope(C.Lichtarowicz1965, h, cb.DischargeHoleState(Re=1e6))
    assert (e.re, e.l_over_d) == (True, False)


def test_a_fixed_cd_carries_no_envelope_flags() -> None:
    e = FilmCoolingElement(
        "film", "cool", "gas", hole_diameter=1e-3, holes_per_row=4, lateral_pitch=4e-3,
        wall_thickness=20e-3, angle_deg=90.0, hole_type="cylindrical",
        correlation="fixed", Cd=0.7,
    )  # fmt: skip
    d = _solve(e)
    assert all(d[k] == 0.0 for k in KEYS) and d["Cd_in_range"] == 1.0


@pytest.mark.parametrize(
    "args,named",
    [
        ((0.1, 1.0, 50.0, False), ()),
        ((0.1, 1.0, 50.0, True), ()),
        ((0.001, 1.0, 50.0, False), ("velocity_ratio",)),
        ((0.7, 1.0, 50.0, True), ("velocity_ratio",)),
        ((0.1, 6.0, 50.0, False), ("l_over_d",)),
        # Injection keeps L/D 0.12-2.0; suction only Fig. 11's band 0.6-1.0.
        ((0.1, 2.0, 50.0, False), ()),
        ((0.1, 2.0, 50.0, True), ("l_over_d",)),
        ((0.1, 0.4, 50.0, True), ("l_over_d",)),
        # Above Vi/V_inf 0.35 the injection L/D independence is unconfirmed.
        ((0.45, 2.0, 50.0, False), ("l_over_d",)),
        ((0.1, 1.0, 120.0, False), ("grazing_velocity",)),
    ],
)
def test_rogers_hersh_envelope_names_each_variable(args: tuple, named: tuple) -> None:
    env = cb.rogers_hersh_1976_envelope(*args)
    for v in ("velocity_ratio", "l_over_d", "grazing_velocity"):
        assert getattr(env, v) == (v in named), v
    assert env.in_range() == (not named)
