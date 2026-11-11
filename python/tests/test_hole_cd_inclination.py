"""The hole-inclination Cd term (Idelchik diagram 3-2 in Eq. 4-20).

Default "auto": on for a cylindrical hole on IdelchikThick, off for a shaped
hole and for every other correlation; a normal hole is unchanged bit for
bit. Forced on, it composes with IdelchikThick only; it reports
Cd_inclination_dzeta. The accuracy evidence (Hay 1983, Gritsch 1998,
cross-source) is validation/cooling/hole_cd_inclination_runner.py; the
ratio test below pins the property the term exists for.
"""

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
THICK = cb.DischargeCdCorrelation.Idelchik1966Thick


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
    return EffusionPlateElement("eff", "cool", "gas", **args)


def _film(**kw) -> FilmCoolingElement:
    args = {
        "hole_diameter": 1e-3,
        "holes_per_row": 10,
        "lateral_pitch": 3e-3,
        "wall_thickness": 3e-3,
        "angle_deg": 30.0,
        "hole_type": "cylindrical",
    }
    args.update(kw)
    return FilmCoolingElement("film", "cool", "gas", **args)


def _solve(element) -> dict:
    g = FlowNetwork()
    g.add_node(PressureBoundary("cool", Pt=1.3e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("gas", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(element)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(g).solve()
    assert r["__success__"], r.get("__message__")
    return r


# --- the C++ term ---------------------------------------------------------


@pytest.mark.parametrize(
    "delta,zeta_in",
    [(20, 0.96), (30, 0.91), (45, 0.81), (60, 0.70), (70, 0.63), (80, 0.56), (90, 0.50)],
)
def test_the_increment_is_diagram_3_2s_table(delta: float, zeta_in: float) -> None:
    assert 0.5 + cb.discharge_inclination_dzeta(delta) == pytest.approx(zeta_in, abs=0.005)


def test_a_normal_hole_has_exactly_no_increment() -> None:
    assert cb.discharge_inclination_dzeta(90.0) == 0.0
    with pytest.raises(ValueError):
        cb.discharge_inclination_dzeta(0.0)


def test_the_geometry_carries_the_angle_and_defaults_to_normal() -> None:
    h = cb.DischargeHoleGeometry(d=1e-3, L=4e-3)
    assert h.inclination_deg == 90.0 and not h.is_inclined()
    h.inclination_deg = 30.0
    assert h.is_inclined()


def test_the_derivative_contract_holds_with_the_term_on() -> None:
    h = cb.DischargeHoleGeometry(d=1e-3, L=4e-3, inclination_deg=30.0)
    for re in (500.0, 5e3, 5e4):
        cd, dcd, du = cb.discharge_cd_and_derivatives(THICK, h, cb.DischargeHoleState(Re=re))
        step = 1e-6 * re
        fd = (
            cb.discharge_cd(THICK, h, cb.DischargeHoleState(Re=re + step))
            - cb.discharge_cd(THICK, h, cb.DischargeHoleState(Re=re - step))
        ) / (2.0 * step)
        assert dcd == pytest.approx(fd, rel=1e-5, abs=1e-12)
        assert du == 0.0


@pytest.mark.parametrize(
    "corr",
    ["McGreehanSchotsch1988", "Lichtarowicz1965", "Idelchik1966Sharp", "Idelchik1966Rounded"],
)
def test_other_correlations_refuse_an_inclined_hole(corr: str) -> None:
    h = cb.DischargeHoleGeometry(d=1e-3, L=4e-3, inclination_deg=30.0)
    with pytest.raises(ValueError, match="Idelchik1966Thick"):
        cb.discharge_cd(getattr(cb.DischargeCdCorrelation, corr), h, cb.DischargeHoleState(Re=2e4))


# --- the elements -----------------------------------------------------------


@pytest.mark.parametrize("make", [_plate, _film])
def test_switched_off_it_is_reported_as_zero(make) -> None:
    e = make(cd_inclination="none")
    assert e.cd_inclination == "none"
    d = _solve(e)["__element_diag__"][e.id]
    assert d["Cd_inclination_dzeta"] == 0.0


@pytest.mark.parametrize("make", [_plate, _film])
def test_the_default_is_the_term_on_a_cylindrical_idelchik_thick_hole(make) -> None:
    # "auto" (the default since the 2026-11-10 hole-Cd audit) resolves to the
    # term for a cylindrical hole on IdelchikThick, and is the explicit
    # opt-in exactly, not a near copy of it.
    e = make()
    assert e.cd_inclination == "idelchik_3_2"
    auto = _solve(e)
    forced = _solve(make(cd_inclination="idelchik_3_2"))
    assert auto[f"{e.id}.m_dot"] == forced[f"{e.id}.m_dot"]
    assert auto["__element_diag__"][e.id] == forced["__element_diag__"][e.id]


def test_the_default_leaves_a_shaped_hole_without_the_term() -> None:
    # No shape term offsets it, and alone it moves Gritsch 1998's fan and
    # laidback holes away from the data (audit 2026-11-10), so a shaped
    # hole keeps the length-only Cd by default.
    shaped = {"hole_type": "shaped", "area_ratio": 3.0, "breakout_width": 2.5e-3}
    e = _film(**shaped)
    assert e.cd_inclination == "none"
    auto = _solve(e)
    off = _solve(_film(cd_inclination="none", **shaped))
    assert auto["film.m_dot"] == off["film.m_dot"]
    assert auto["__element_diag__"]["film"]["Cd_inclination_dzeta"] == 0.0
    # Forcing it on a shaped hole stays the user's choice.
    assert _film(cd_inclination="idelchik_3_2", **shaped).cd_inclination == "idelchik_3_2"


@pytest.mark.parametrize("make", [_plate, _film])
@pytest.mark.parametrize("corr", ["McGreehanSchotsch", "Lichtarowicz", "fixed"])
def test_the_default_on_a_non_composing_correlation_is_off_not_refused(make, corr) -> None:
    kw = {"Cd": 0.7} if corr == "fixed" else {}
    e = make(correlation=corr, **kw)
    assert e.cd_inclination == "none"
    assert e._hole_cd.inclination_deg == 90.0


@pytest.mark.parametrize("make", [_plate, _film])
def test_the_term_lowers_the_cd_by_the_diagram_3_2_increment(make) -> None:
    off = _solve(make(cd_inclination="none"))["__element_diag__"]
    on = _solve(make(cd_inclination="idelchik_3_2"))["__element_diag__"]
    eid = next(iter(off))
    dz = cb.discharge_inclination_dzeta(30.0)
    assert on[eid]["Cd_inclination_dzeta"] == pytest.approx(dz)
    assert on[eid]["Cd"] < off[eid]["Cd"]
    assert on[eid]["m_dot"] < off[eid]["m_dot"]


@pytest.mark.parametrize("make", [_plate, _film])
@pytest.mark.parametrize("inc", ["auto", "idelchik_3_2"])
def test_a_normal_hole_is_unchanged_by_the_term(make, inc: str) -> None:
    # Bit for bit: every solved unknown and every diagnostic.
    off = _solve(make(angle_deg=90.0, cd_inclination="none"))
    on = _solve(make(angle_deg=90.0, cd_inclination=inc))
    # "auto" on a normal hole names what it applies: no inclined inlet.
    assert make(angle_deg=90.0, cd_inclination=inc).cd_inclination == (
        "none" if inc == "auto" else "idelchik_3_2"
    )
    assert {k: v for k, v in on.items() if not k.startswith("__")} == {
        k: v for k, v in off.items() if not k.startswith("__")
    }
    assert on["__element_diag__"] == off["__element_diag__"]
    assert on["__element_diag__"][next(iter(on["__element_diag__"]))]["Cd_inclination_dzeta"] == 0.0


@pytest.mark.parametrize("make", [_plate, _film])
@pytest.mark.parametrize("corr", ["McGreehanSchotsch", "Lichtarowicz", "fixed"])
def test_an_opt_in_on_a_non_composing_correlation_is_refused(make, corr: str) -> None:
    kw = {"Cd": 0.7} if corr == "fixed" else {}
    with pytest.raises(ValueError, match="does not compose"):
        make(correlation=corr, cd_inclination="idelchik_3_2", **kw)


def test_an_unknown_inclination_model_is_refused() -> None:
    with pytest.raises(ValueError, match="cd_inclination"):
        _plate(cd_inclination="weisbach")


def test_the_range_flag_reports_a_hole_too_short_to_reattach() -> None:
    # t/d 0.5 at 30 deg: L/d 1.0 < the l/d 2 floor -> flagged, still solved.
    short = _solve(_plate(wall_thickness=0.4e-3, cd_inclination="idelchik_3_2"))
    assert short["__element_diag__"]["eff"]["Cd_in_range"] == 0.0


def test_the_30_over_90_ratio_is_hays_order() -> None:
    # Hay et al. (1983) measured Cd(30)/Cd(90) 0.85-0.93 (no coolant
    # crossflow, L/D 6), carpets 0.90-0.91. Off, the angle enters only
    # through the length, and at equal L/D the ratio is exactly 1; on, it
    # must fall into Hay's range. Cross-source accuracy, not a fit: the
    # coefficients are Idelchik's.
    def cd(alpha: float, inc: str) -> float:
        h = cb.DischargeHoleGeometry(d=3.5e-3, L=6 * 3.5e-3)
        if inc == "on":
            h.inclination_deg = alpha
        return cb.discharge_cd(THICK, h, cb.DischargeHoleState(Re=5e4))

    assert cd(30.0, "off") / cd(90.0, "off") == 1.0
    ratio = cd(30.0, "on") / cd(90.0, "on")
    assert 0.85 <= ratio <= 0.93
