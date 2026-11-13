"""HoleCd: the hole discharge-coefficient family's options as one spec.

Every hole element (OrificeElement, the effusion, film, slot and impingement
elements) and both builders (EffusionLiner, ImpingementArray) take
``hole_cd=HoleCd(...)``: the correlation, a fixed Cd and each term's model.
A field left at None/"auto" is the element's own default, which differs per
element and is exactly what the element did before the spec. The
per-element keywords it replaces that shipped (``correlation``, ``Cd``,
``cd_crossflow``) still work, deprecated: their warning is checked in
test_deprecations.py, and here that they build the same element, residuals
and Jacobian, bit for bit. The pre-release ``cd_inclination`` and
``cd_exit_crossflow`` never shipped and are gone (TypeError).
"""

from __future__ import annotations

import dataclasses
import warnings

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    EffusionLiner,
    EffusionPlateElement,
    FilmCoolingElement,
    FlowNetwork,
    HoleCd,
    ImpingementArray,
    ImpingementPlateElement,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    SlotFilmElement,
)

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
D_JET = 0.00254


def _plate(**kw) -> EffusionPlateElement:
    a = {
        "hole_diameter": 0.8e-3,
        "wall_thickness": 2e-3,
        "pitch": 6e-3,
        "panel_length": 0.1,
        "panel_width": 0.1,
        "angle_deg": 30.0,
    }
    a.update(kw)
    return EffusionPlateElement("h", "cool", "ch", **a)


def _film(**kw) -> FilmCoolingElement:
    a = {"hole_diameter": 1e-3, "holes_per_row": 10, "lateral_pitch": 3e-3, "wall_thickness": 3e-3}
    a.update(kw)
    return FilmCoolingElement("h", "cool", "ch", **a)


def _shaped(**kw) -> FilmCoolingElement:
    return _film(hole_type="shaped", area_ratio=3.0, breakout_width=2.5e-3, **kw)


def _slot(**kw) -> SlotFilmElement:
    return SlotFilmElement("h", "cool", "ch", slot_height=1e-3, span=0.05, lip_thickness=4e-3, **kw)


def _row(**kw) -> ImpingementPlateElement:
    return ImpingementPlateElement(
        "p", "s", "c", d_jet=D_JET, xn_d=5.0, yn_d=4.0, z_d=1.0, span=48 * D_JET,
        plate_thickness=D_JET, **kw,
    )  # fmt: skip


def _array(**kw) -> ImpingementArray:
    return ImpingementArray(
        "ia", n_rows=4, d_jet=D_JET, xn_d=5.0, yn_d=4.0, z_d=1.0, span=48 * D_JET,
        plate_thickness=D_JET, **kw,
    )  # fmt: skip


def _chamber(element, dp: float = 0.03) -> FlowNetwork:
    """Hot gas through a merge chamber the holes discharge into."""
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("gas_in", m_dot=2.0, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("exit", Pt=5e5, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cool", Pt=5e5 * (1.0 + dp), Tt=600.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("ch", area=0.01, main_inlet="feed"))
    g.add_element(LosslessConnectionElement("feed", "gas_in", "ch"))
    g.add_element(LosslessConnectionElement("tail", "ch", "exit"))
    g.add_element(element)
    return g


def _impingement(array: ImpingementArray) -> FlowNetwork:
    """The array under a bypass as strong as its supply, so the crossflow
    law is not just the hole's own Cd (CR-3630 Table 4.1 then moves F)."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("sup", Pt=3e5 * 1.03, Tt=600.0, Y=Y_AIR))
    g.add_node(PressureBoundary("ex", Pt=3e5, Tt=600.0, Y=Y_AIR))
    g.add_node(PressureBoundary("byp", Pt=3e5 * 1.03, Tt=600.0, Y=Y_AIR))
    array.add_to(g, "sup", "ex", bypass="byp")
    return g


def _liner(**kw) -> tuple[FlowNetwork, EffusionLiner]:
    g = FlowNetwork()
    g.add_node(PressureBoundary("cin", Pt=1.5e5, Tt=600.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cout", Pt=1.1e5, Tt=600.0, Y=Y_AIR))
    g.add_node(PressureBoundary("gas", Pt=1.0e5, Tt=1500.0, Y=Y_AIR))
    liner = EffusionLiner(
        "ln", n_segments=3, length=0.152, width=0.152, duct_height=0.01, hole_diameter=3.27e-3,
        wall_thickness=6.3e-3, pitch_x=15.24e-3, pitch_y=15.24e-3, angle_deg=30.0, **kw,
    )  # fmt: skip
    liner.add_to(g, "cin", "cout", "gas")
    return g, liner


def _fj(g: FlowNetwork) -> list[np.ndarray]:
    """F and dense J at the start guess and one perturbation of it, and
    every hole element's diagnostics at the start guess (as arrays)."""
    g.resolve_all_topology()
    s = NetworkSolver(g)
    x0 = np.array(s._build_x0(), dtype=float)
    rng = np.random.default_rng(7)
    out = []
    for x in (x0, x0 * (1.0 + 0.02 * rng.standard_normal(x0.size))):
        F, J = s._residuals_and_jacobian(x)
        out += [np.asarray(F, dtype=float), np.asarray(J.toarray(), dtype=float)]
    return out


def _legacy(build, **kw):
    """Build through the deprecated keywords, their warning silenced."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        return build(**kw)


# --- the spec ---------------------------------------------------------------


def test_the_spec_is_frozen_and_defaults_to_every_elements_own_default() -> None:
    spec = HoleCd()
    assert (spec.correlation, spec.Cd) == (None, None)
    assert (spec.inclination, spec.inlet_crossflow, spec.exit_crossflow) == ("auto",) * 3
    with pytest.raises(dataclasses.FrozenInstanceError):
        spec.Cd = 0.7  # type: ignore[misc]


@pytest.mark.parametrize(
    ("kw", "match"),
    [
        ({"correlation": "Miller"}, "correlation 'Miller' is not one of"),
        ({"correlation": "Stolz"}, "correlation 'Stolz' is not one of"),
        ({"correlation": "IdelchikThick", "Cd": 0.8}, "Cd=0.8 with correlation='IdelchikThick'"),
        ({"correlation": "fixed", "Cd": float("nan")}, "must be a finite number"),
        ({"Cd": float("inf")}, "must be a finite number"),
        ({"inclination": "weisbach"}, "inclination 'weisbach' is not one of"),
        ({"inlet_crossflow": "duct"}, "inlet_crossflow 'duct' is not one of"),
        ({"exit_crossflow": "bogus"}, "exit_crossflow 'bogus' is not one of"),
    ],
)
def test_the_spec_refuses_what_no_element_can_build(kw: dict, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        HoleCd(**kw)


@pytest.mark.parametrize(
    ("build", "spec", "match"),
    [
        # No inclined hole: the inclination term does not exist there.
        (_row, HoleCd(inclination="idelchik_3_2"), "does not apply"),
        (_slot, HoleCd(inclination="idelchik_3_2"), "does not apply"),
        # The exit laws belong to their own hole families.
        (_row, HoleCd(exit_crossflow="rowbury_2001"), "does not apply"),
        (_slot, HoleCd(exit_crossflow="rowbury_2001"), "does not apply"),
        (_plate, HoleCd(exit_crossflow="rogers_hersh_1976"), "does not apply"),
        (_film, HoleCd(exit_crossflow="florschuetz_metzger_su_1982"), "does not apply"),
        # McGreehan-Schotsch's supply crossflow needs a duct and its own Cd.
        (_plate, HoleCd(correlation="McGreehanSchotsch", inlet_crossflow="mcgreehan_schotsch"),
         "plenum-fed"),
        (_film, HoleCd(correlation="McGreehanSchotsch", inlet_crossflow="mcgreehan_schotsch"),
         "plenum-fed"),
        (_row, HoleCd(inlet_crossflow="mcgreehan_schotsch"), "plenum-fed"),
        (lambda hole_cd: _plate(hole_cd=hole_cd, crossflow_segments=(None, None),
                                crossflow_area=1e-3),
         HoleCd(inlet_crossflow="mcgreehan_schotsch"), "does not compose"),
        # CR-3630 Table 4.1 covers the plates its authors tested only.
        (lambda hole_cd: ImpingementPlateElement(
            "p", "s", "c", d_jet=D_JET, xn_d=5.0, yn_d=8.0, z_d=2.0, span=0.1,
            plate_thickness=D_JET, hole_cd=hole_cd),
         HoleCd(exit_crossflow="florschuetz_metzger_su_1982"), "Table 4.1"),
    ],
)  # fmt: skip
def test_an_element_refuses_a_term_it_has_not(build, spec: HoleCd, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        build(hole_cd=spec)


def test_a_plain_orifice_refuses_the_cooling_hole_terms() -> None:
    for spec in (HoleCd(inclination="idelchik_3_2"), HoleCd(exit_crossflow="rowbury_2001")):
        with pytest.raises(ValueError, match="does not apply"):
            OrificeElement("o", "a", "b", diameter=0.01, hole_cd=spec)


@pytest.mark.parametrize("build", [_plate, _film, _slot, _row, _array])
def test_hole_cd_must_be_a_hole_cd(build) -> None:
    with pytest.raises(ValueError, match="hole_cd is a HoleCd"):
        build(hole_cd={"correlation": "fixed"})


def test_the_element_holds_its_resolved_spec() -> None:
    # The correlation and Cd as used; the terms as given.
    assert _plate().hole_cd == HoleCd(correlation="IdelchikThick")
    assert _row().hole_cd == HoleCd(correlation="fixed", Cd=cb.FLORSCHUETZ_1981_DEFAULT_CD)
    assert _slot(hole_cd=HoleCd(correlation="fixed", Cd=0.7)).hole_cd.Cd == 0.7
    o = OrificeElement("o", "a", "b", diameter=0.01)
    assert o.hole_cd == HoleCd(correlation="ReaderHarrisGallagher")


# --- the per-element defaults ("auto") ---------------------------------------


def test_auto_is_each_elements_own_default() -> None:
    # Inclination: on for an inclined cylindrical IdelchikThick hole, off for
    # a normal or shaped one and for any other correlation.
    assert _plate().cd_inclination == "idelchik_3_2"
    assert _plate(angle_deg=90.0).cd_inclination == "none"
    assert _shaped().cd_inclination == "none"
    assert _film(hole_cd=HoleCd(correlation="McGreehanSchotsch")).cd_inclination == "none"
    # Exit crossflow: off on film and effusion holes, Rogers & Hersh on a row.
    assert _plate().cd_exit_crossflow == "none" and _plate()._hole_cd.exit_xf is None
    assert _row().cd_crossflow == "rogers_hersh_1976"
    # Supply crossflow: read on a duct-fed panel only.
    assert _plate()._hole_cd.supply_duct is False
    duct = _plate(
        hole_cd=HoleCd(correlation="McGreehanSchotsch"),
        crossflow_segments=(None, None),
        crossflow_area=1e-3,
    )
    assert duct._hole_cd.supply_duct is True


def test_inlet_crossflow_none_takes_a_duct_fed_panel_as_plenum_fed() -> None:
    # The caller's explicit choice: U1/Vi held at 0, and the panel is no
    # longer refused for a correlation without the crossflow term.
    g, off = _liner(hole_cd=HoleCd(inlet_crossflow="none"))
    for el in (g.elements[off.panel_id(i)] for i in range(1, 4)):
        assert el._hole_cd.supply_duct is False
    g, _ = _liner(hole_cd=HoleCd(correlation="IdelchikThick", inlet_crossflow="none"))
    g.resolve_all_topology()
    g.validate()
    g, _ = _liner(hole_cd=HoleCd(correlation="IdelchikThick"))
    g.resolve_all_topology()
    with pytest.raises(ValueError, match="duct-fed panel needs the supply-side crossflow"):
        g.validate()


# --- builders pass it through ---------------------------------------------------


def test_the_liner_passes_its_spec_to_every_panel_with_its_own_default() -> None:
    g, liner = _liner()
    for i in range(1, 4):
        p = g.elements[liner.panel_id(i)]
        assert p.correlation == "McGreehanSchotsch" and p._hole_cd.supply_duct
    g, liner = _liner(hole_cd=HoleCd(exit_crossflow="none", inclination="none"))
    for i in range(1, 4):
        p = g.elements[liner.panel_id(i)]
        assert p.hole_cd == HoleCd(
            correlation="McGreehanSchotsch", inclination="none", exit_crossflow="none"
        )
    # The default and the spelled-out default are one network.
    a, b = _fj(_liner()[0]), _fj(_liner(hole_cd=HoleCd(correlation="McGreehanSchotsch"))[0])
    assert all(np.array_equal(x, y) for x, y in zip(a, b, strict=True))


def test_the_array_passes_its_spec_to_every_row() -> None:
    g = _impingement(_array())
    rows = [e for e in g.elements.values() if isinstance(e, ImpingementPlateElement)]
    assert len(rows) == 4 and all(r.cd_crossflow == "rogers_hersh_1976" for r in rows)
    spec = HoleCd(correlation="IdelchikThick", exit_crossflow="none")
    g = _impingement(_array(hole_cd=spec))
    rows = [e for e in g.elements.values() if isinstance(e, ImpingementPlateElement)]
    assert all(r.hole_cd == spec and r.cd_crossflow == "none" for r in rows)


# --- the deprecated keywords build the same element -----------------------------


_CHAMBER_CASES = {
    "plate|default": (_plate, {}, HoleCd()),
    "plate|fixed": (_plate, {"correlation": "fixed", "Cd": 0.7}, HoleCd(correlation="fixed", Cd=0.7)),
    "plate|mgs": (_plate, {"correlation": "McGreehanSchotsch"},
                  HoleCd(correlation="McGreehanSchotsch")),
    "film|default": (_film, {}, HoleCd()),
    "film|lich": (_film, {"correlation": "Lichtarowicz"}, HoleCd(correlation="Lichtarowicz")),
    "film|fixed": (_film, {"correlation": "fixed", "Cd": 0.65},
                   HoleCd(correlation="fixed", Cd=0.65)),
    "shaped|mgs": (_shaped, {"correlation": "McGreehanSchotsch"},
                   HoleCd(correlation="McGreehanSchotsch")),
    "slot|default": (_slot, {}, HoleCd()),
    "slot|fixed": (_slot, {"correlation": "fixed", "Cd": 0.6}, HoleCd(correlation="fixed", Cd=0.6)),
    "slot|mgs": (_slot, {"correlation": "McGreehanSchotsch"},
                 HoleCd(correlation="McGreehanSchotsch")),
}  # fmt: skip


@pytest.mark.parametrize("case", sorted(_CHAMBER_CASES))
def test_the_old_keywords_and_the_spec_are_one_element(case: str) -> None:
    build, old, spec = _CHAMBER_CASES[case]
    a = _legacy(build, **old)
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        b = build(hole_cd=spec)
    assert a.hole_cd == b.hole_cd
    assert (a.correlation, a.Cd) == (b.correlation, b.Cd)
    fa, fb = _fj(_chamber(a)), _fj(_chamber(b))
    assert all(np.array_equal(x, y) for x, y in zip(fa, fb, strict=True))


@pytest.mark.parametrize(
    ("old", "spec"),
    [
        ({}, HoleCd()),
        ({"cd_crossflow": "none"}, HoleCd(exit_crossflow="none")),
        ({"cd_crossflow": "florschuetz_metzger_su_1982", "Cd": 0.85},
         HoleCd(Cd=0.85, exit_crossflow="florschuetz_metzger_su_1982")),
        ({"correlation": "IdelchikThick"}, HoleCd(correlation="IdelchikThick")),
        ({"correlation": "McGreehanSchotsch", "cd_crossflow": "rogers_hersh_1976"},
         HoleCd(correlation="McGreehanSchotsch", exit_crossflow="rogers_hersh_1976")),
    ],
)  # fmt: skip
def test_the_old_array_fields_and_the_spec_are_one_network(old: dict, spec: HoleCd) -> None:
    fa = _fj(_impingement(_legacy(_array, **old)))
    fb = _fj(_impingement(_array(hole_cd=spec)))
    assert all(np.array_equal(x, y) for x, y in zip(fa, fb, strict=True))
    row_a = _legacy(_row, **old)
    row_b = _row(hole_cd=spec)
    assert row_a.hole_cd == row_b.hole_cd and row_a.cd_crossflow == row_b.cd_crossflow


@pytest.mark.parametrize(
    ("build", "old"),
    [
        (_plate, {"correlation": "IdelchikThick"}),
        (_plate, {"Cd": 0.7}),
        (_film, {"Cd": 0.7}),
        (_slot, {"correlation": "fixed"}),
        (_row, {"cd_crossflow": "none"}),
        (_array, {"correlation": "fixed"}),
    ],
)
def test_the_spec_and_an_old_keyword_together_are_refused(build, old: dict) -> None:
    with pytest.raises(ValueError, match="deprecated keyword"):
        build(hole_cd=HoleCd(), **old)


# --- the plain orifice ------------------------------------------------------------


def _orifice_net(el: OrificeElement) -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=1.2e5, Tt=600.0))
    g.add_node(PressureBoundary("b", Pt=1.0e5, Tt=300.0))
    g.add_element(el)
    return g


@pytest.mark.parametrize("corr", ["IdelchikThick", "McGreehanSchotsch", "ReaderHarrisGallagher"])
def test_a_plain_orifice_takes_its_correlation_either_way(corr: str) -> None:
    kw = {"diameter": 0.01, "plate_thickness": 0.02}
    a = OrificeElement("o", "a", "b", correlation=corr, **kw)
    b = OrificeElement("o", "a", "b", hole_cd=HoleCd(correlation=corr), **kw)
    c = OrificeElement("o", "a", "b", correlation=corr, hole_cd=HoleCd(correlation=corr), **kw)
    assert a.hole_cd == b.hole_cd == c.hole_cd
    fa, fb = _fj(_orifice_net(a)), _fj(_orifice_net(b))
    assert all(np.array_equal(x, y) for x, y in zip(fa, fb, strict=True))
    # Its own keyword is not deprecated: no warning.
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        OrificeElement("o", "a", "b", correlation="fixed", Cd=0.7, **kw)


def test_a_plain_orifice_refuses_a_keyword_that_disagrees_with_its_spec() -> None:
    with pytest.raises(ValueError, match="disagree"):
        OrificeElement(
            "o", "a", "b", diameter=0.01, correlation="ReaderHarrisGallagher",
            hole_cd=HoleCd(correlation="IdelchikThick"),
        )  # fmt: skip
    with pytest.raises(ValueError, match="disagree"):
        OrificeElement(
            "o", "a", "b", diameter=0.01, Cd=0.6, hole_cd=HoleCd(correlation="fixed", Cd=0.7)
        )


# --- serialization -----------------------------------------------------------


def test_a_network_round_trips_through_its_spec_without_a_warning() -> None:
    g = FlowNetwork()
    g.add_node(PressureBoundary("cool", Pt=1.3e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PlenumNode("ch"))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(_film(hole_cd=HoleCd(correlation="fixed", Cd=0.7, inclination="none")))
    g.add_element(OrificeElement("o", "ch", "out", diameter=0.01, correlation="IdelchikThick"))
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        again = FlowNetwork.from_dict(g.to_dict())
    for eid in ("h", "o"):
        assert again.elements[eid].hole_cd == g.elements[eid].hole_cd
    assert "correlation" not in g.to_dict()["elements"]["h"]["kwargs"]


@pytest.mark.parametrize("build", [_plate, _film])
@pytest.mark.parametrize("kw", ["cd_inclination", "cd_exit_crossflow"])
def test_the_pre_release_keywords_are_gone(build, kw: str) -> None:
    # Never in a release, so no deprecation: an ordinary unknown keyword.
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        build(**{kw: "none"})


def test_inlet_crossflow_none_on_a_duct_fed_panel_is_reported() -> None:
    # The term off is the caller's choice, not hidden: the duct's U1/Vi is
    # still reported, the Cd does not read it, and Cd_inlet_crossflow_ignored
    # says the panel is taken as plenum-fed. Panel 2 of a 3-station liner.
    def panel(spec: HoleCd | None) -> dict:
        g, liner = _liner(**({} if spec is None else {"hole_cd": spec}))
        r = NetworkSolver(g).solve()
        assert r["__success__"], r.get("__message__")
        return r["__element_diag__"][liner.panel_id(2)]

    # The pre-0.12 drive (supply_drive="static": station static, node T0),
    # measured on this liner (solved values, 3-digit pins): the duct's
    # crossflow is strong (U1/Vi ~2) and McGreehan-Schotsch's term halves
    # the Cd; switched off, the Cd is the plenum-fed hole's.
    on = panel(HoleCd(supply_drive="static"))
    off = panel(HoleCd(inlet_crossflow="none", supply_drive="static"))
    assert on["Cd_inlet_crossflow_ignored"] == 0.0 and on["Cd_inlet_crossflow_degraded"] == 1.0
    assert on["Cd_inlet_crossflow_U1_over_Vi"] == pytest.approx(2.067, abs=0.001)
    assert on["Cd"] == pytest.approx(0.559, abs=0.001)
    assert off["Cd_inlet_crossflow_ignored"] == 1.0
    assert off["Cd_inlet_crossflow_U1_over_Vi"] == pytest.approx(
        2.349, abs=0.001
    )  # reported, read by nothing
    assert (
        off["Cd_inlet_crossflow_beyond_8pct"] == 0.0 and off["Cd_inlet_crossflow_degraded"] == 0.0
    )
    assert off["Cd"] == pytest.approx(0.807, abs=0.001)
    # The default, every hole from the duct total to the exit static (#1087).
    # With the term: Eq. 17 converted exactly, so the flow is the static
    # drive's up to the static-temperature correction (+1.6 % here) and the
    # reported Cd is on the total basis (83 % of this liner's drive is the
    # duct's dynamic head, Cd_supply_drive_head_fraction).
    on_t, off_t = panel(None), panel(HoleCd(inlet_crossflow="none"))
    assert on_t["Cd_supply_drive_head_fraction"] == pytest.approx(0.829, abs=0.001)
    assert on_t["Cd_inlet_crossflow_U1_over_Vi"] == pytest.approx(2.089, abs=0.001)
    assert on_t["Cd"] == pytest.approx(0.240, abs=0.001)
    assert on_t["m_dot"] / on["m_dot"] - 1.0 == pytest.approx(0.016, abs=0.002)
    # Without the term, the default f = 0 (owner, 2026-11-13; the held
    # normal-hole data): the plain Cd referenced to the station static and
    # converted onto the total exactly, so again the static drive's flow up to
    # the static-temperature correction (+1.8 % here).
    assert off_t["Cd_supply_recovery_fraction"] == 0.0
    assert off_t["Cd_supply_drive_head_fraction"] == pytest.approx(0.864, abs=0.001)
    assert off_t["Cd_supply_recovery_q1"] == pytest.approx(28406.0, rel=1e-3)
    assert off_t["Cd_supply_recovery_ratio"] == pytest.approx(0.386, abs=0.001)
    assert off_t["Cd"] == pytest.approx(0.312, abs=0.001)
    assert off_t["m_dot"] / off["m_dot"] - 1.0 == pytest.approx(0.018, abs=0.002)
    # supply_recovery=1.0: the plain Cd total -> static, the full head
    # recovered (the co-aligned behaviour): this panel bleeds 2.3x the static
    # drive's, and the reported U1/Vi (on the station static's jet velocity,
    # McGreehan-Schotsch's definition) rises with the duct's speed.
    off_t = panel(HoleCd(inlet_crossflow="none", supply_recovery=1.0))
    assert off_t["Cd_supply_recovery_fraction"] == 1.0
    assert off_t["Cd_supply_recovery_ratio"] == 1.0
    assert off_t["Cd_supply_drive_head_fraction"] == pytest.approx(0.991, abs=0.001)
    assert off_t["Cd_inlet_crossflow_U1_over_Vi"] == pytest.approx(10.265, abs=0.001)
    assert off_t["Cd"] == pytest.approx(0.801, abs=0.001)
    assert off_t["m_dot"] / off["m_dot"] == pytest.approx(2.29, abs=0.02)


@pytest.mark.parametrize(
    ("build", "match"),
    [
        # The element's default correlation, not the user's, refuses the Cd.
        (lambda: _plate(hole_cd=HoleCd(Cd=0.7)),
         r"EffusionPlateElement 'h': HoleCd: Cd=0.7 with correlation='IdelchikThick'.*"
         r"this element's default"),
        (lambda: _film(hole_cd=HoleCd(Cd=0.7)), r"FilmCoolingElement 'h': .*element's default"),
        (lambda: _slot(hole_cd=HoleCd(Cd=0.7)), r"SlotFilmElement 'h': .*element's default"),
        (lambda: _liner(hole_cd=HoleCd(Cd=0.7)), r"EffusionLiner 'ln': .*element's default"),
        (lambda: OrificeElement("o", "a", "b", diameter=0.01, correlation="Foo"),
         r"OrificeElement 'o': HoleCd: correlation 'Foo' is not one of"),
        (lambda: ImpingementPlateElement(
            "p", "s", "c", d_jet=D_JET, xn_d=5.0, yn_d=8.0, z_d=2.0, span=0.1,
            plate_thickness=D_JET, hole_cd=HoleCd(exit_crossflow="florschuetz_metzger_su_1982")),
         r"ImpingementPlateElement 'p': HoleCd\(exit_crossflow='florschuetz_metzger_su_1982'\)"),
    ],
)  # fmt: skip
def test_every_refusal_names_the_element(build, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        build()


def test_an_array_with_a_spec_survives_dataclasses_replace() -> None:
    arr = _array(hole_cd=HoleCd(exit_crossflow="none"))
    again = dataclasses.replace(arr, n_rows=2)
    assert again.hole_cd == arr.hole_cd and again._rows_hole_cd == HoleCd(exit_crossflow="none")
