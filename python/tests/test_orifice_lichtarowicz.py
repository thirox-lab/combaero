"""Lichtarowicz (1965) long-orifice Cd through OrificeElement, and the range flag (#465).

The correlation existed in C++ but no element could select it. A jet plate's
holes are long (t/d 1-3) and run at Re 1e3-7e4, which is where it belongs.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network import (
    EffusionPlateElement,
    FlowNetwork,
    HoleCd,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)

D = 2.0e-3


def _boundary(name: str, Pt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt = Pt
    b.Tt = 300.0
    b.Y = cb.species.dry_air_mass()
    return b


def _solve(el: OrificeElement) -> dict:
    g = FlowNetwork()
    g.add_node(_boundary("a", 1.05e5))
    g.add_node(_boundary("b", 1.00e5))
    g.add_element(el)
    res = NetworkSolver(g).solve()
    assert res["__success__"], res.get("__message__")
    return res["__element_diag__"][el.id]


def test_the_element_returns_the_c_plus_plus_value() -> None:
    o = OrificeElement("o", "a", "b", diameter=D, correlation="Lichtarowicz", plate_thickness=3 * D)
    diag = _solve(o)
    mu = cb.transport_state(300.0, 1.05e5, cb.species.dry_air()).mu
    Re = 4.0 * diag["m_dot"] / (math.pi * D * mu)
    ref = cb.discharge_cd(
        cb.DischargeCdCorrelation.Lichtarowicz1965,
        cb.DischargeHoleGeometry(d=D, L=3 * D, r=0.0),
        cb.DischargeHoleState(Re=Re),
    )
    assert diag["Cd"] == pytest.approx(ref, rel=1e-3)
    assert diag["Cd_in_range"] == 1.0


def test_a_short_hole_is_refused_at_set_up_naming_the_element() -> None:
    """Below l/d = 1.5 the source reports hysteresis; refuse before solving."""
    o = OrificeElement("short", "a", "b", diameter=D, correlation="Lichtarowicz", plate_thickness=D)
    with pytest.raises(ValueError, match="'short'.*1.5"):
        _solve(o)


def test_the_flag_reports_an_out_of_range_hole() -> None:
    """l/d 1.8 is computed (the source's flat 0.81 band) but outside 2-10."""
    o = OrificeElement(
        "o", "a", "b", diameter=D, correlation="Lichtarowicz", plate_thickness=1.8 * D
    )
    assert _solve(o)["Cd_in_range"] == 0.0


def test_fixed_cd_is_always_in_range() -> None:
    assert (
        _solve(OrificeElement("o", "a", "b", Cd=0.6, diameter=D, correlation="fixed"))[
            "Cd_in_range"
        ]
        == 1.0
    )


def test_an_effusion_plate_can_select_it() -> None:
    e = EffusionPlateElement(
        "p",
        "a",
        "b",
        hole_diameter=D,
        wall_thickness=4 * D,
        pitch=8 * D,
        panel_area=1e-3,
        hole_cd=HoleCd(correlation="Lichtarowicz"),
    )
    diag = _solve(e)
    assert 0.6 < diag["Cd"] < 0.85
