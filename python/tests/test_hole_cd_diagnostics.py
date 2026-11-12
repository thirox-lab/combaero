"""The hole-Cd diagnostic group: one stable shape on every hole element.

Every OrificeElement subclass reports every key of HOLE_CD_DIAGNOSTIC_KEYS,
whether or not a term is on; a term that is off reads model "none",
applied 0 and neutral values. The 0.11 names that shipped are still emitted
beside the group for one release (deprecated since 0.12.0, removed in
0.13.0), equal to their new keys and on the same elements as before; this
file is the only reader of those aliases, as the deprecation contract.
Bit-identity of the values against the pre-group keys was checked over every
stress family, term probes and the GUI examples (PR description).
"""

from __future__ import annotations

import math
import warnings

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    EffusionLiner,
    EffusionPlateElement,
    FilmCoolingElement,
    FlowNetwork,
    HoleCd,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    SlotFilmElement,
)
from combaero.network.hole_cd import (
    HOLE_CD_DIAGNOSTIC_ALIASES,
    HOLE_CD_DIAGNOSTIC_KEYS,
    HOLE_CD_DIAGNOSTIC_TEXT_KEYS,
    LINER_CD_DIAGNOSTIC_ALIASES,
)
from combaero.network.impingement_array import ImpingementArray

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _solve(g: FlowNetwork) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(g).solve()
    assert r["__success__"], r.get("__message__")
    return r


def _chamber(element) -> FlowNetwork:
    """A hole from a coolant plenum into a flowing chamber."""
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("gas_in", m_dot=2.0, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("exit", Pt=5e5, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cool", Pt=5.15e5, Tt=600.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("ch", area=0.01, main_inlet="feed"))
    g.add_element(LosslessConnectionElement("feed", "gas_in", "ch"))
    g.add_element(LosslessConnectionElement("tail", "ch", "exit"))
    g.add_element(element)
    return g


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


def _slot() -> SlotFilmElement:
    return SlotFilmElement("h", "cool", "ch", slot_height=1e-3, span=0.05, lip_thickness=1e-3)


def _metering() -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=3e5, Tt=500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("b", Pt=2.6e5, Tt=500.0, Y=Y_AIR))
    g.add_node(PlenumNode("m1"))
    g.add_node(PlenumNode("m2"))
    g.add_element(ChannelElement("p1", "a", "m1", length=0.5, diameter=0.05))
    g.add_element(OrificeElement("h", "m1", "m2", diameter=0.03))
    g.add_element(ChannelElement("p2", "m2", "b", length=0.5, diameter=0.05))
    return g


def _impingement(law: str) -> tuple[FlowNetwork, ImpingementArray]:
    g = FlowNetwork()
    D = 0.00254
    g.add_node(PressureBoundary("sup", Pt=3e5 * 1.5, Tt=600.0, Y=Y_AIR))
    g.add_node(PressureBoundary("ex", Pt=3e5, Tt=600.0, Y=Y_AIR))
    arr = ImpingementArray(
        "ia", n_rows=6, d_jet=D, xn_d=5.0, yn_d=4.0, z_d=1.0, span=48.0 * D, plate_thickness=D,
        hole_cd=HoleCd(exit_crossflow=law),
    )  # fmt: skip
    arr.add_to(g, "sup", "ex")
    return g, arr


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


def _assert_full_group(d: dict) -> None:
    missing = [k for k in HOLE_CD_DIAGNOSTIC_KEYS if k not in d]
    assert not missing, missing
    for k in HOLE_CD_DIAGNOSTIC_KEYS:
        assert isinstance(d[k], str) == (k in HOLE_CD_DIAGNOSTIC_TEXT_KEYS), k


@pytest.mark.parametrize(
    "build",
    [
        lambda: _chamber(_plate()),
        lambda: _chamber(_film()),
        lambda: _chamber(_film(hole_cd=HoleCd(exit_crossflow="rowbury_2001"))),
        lambda: _chamber(_slot()),
        _metering,
    ],
    ids=["effusion", "film", "film_rowbury", "slot", "metering"],
)
def test_every_hole_element_reports_the_whole_group(build) -> None:
    d = _solve(build())["__element_diag__"]["h"]
    _assert_full_group(d)
    assert d["Cd_exit_crossflow_ratio"] == pytest.approx(d["Cd"] / d["Cd_base"], rel=1e-12)


def test_a_term_that_is_off_reads_neutral_values() -> None:
    d = _solve(_metering())["__element_diag__"]["h"]
    assert d["Cd_correlation"] == "ReaderHarrisGallagher"
    assert d["Cd_base"] == d["Cd"]
    for term in ("inclination", "inlet_crossflow", "exit_crossflow"):
        assert d[f"Cd_{term}_model"] == "none"
        assert d[f"Cd_{term}_applied"] == 0.0
    assert d["Cd_inclination_dzeta"] == 0.0
    assert d["Cd_inlet_crossflow_U1_over_Vi"] == 0.0
    assert d["Cd_exit_crossflow_ratio"] == 1.0 and d["Cd_exit_crossflow_rho_u2"] == 0.0
    assert all(math.isnan(d[f"Cd_exit_crossflow_{k}"]) for k in ("k_perp", "k_used", "k_floor"))
    assert d["Cd_exit_crossflow_in_range"] == 1.0
    for k in ("ignored", "beyond_8pct", "degraded"):
        assert d[f"Cd_inlet_crossflow_{k}"] == 0.0
    # No pre-0.12 name appears where 0.11 had none.
    assert not set(HOLE_CD_DIAGNOSTIC_ALIASES) & set(d)


def test_the_terms_that_are_on_say_so() -> None:
    film = _solve(_chamber(_film(hole_cd=HoleCd(exit_crossflow="rowbury_2001"))))
    d = film["__element_diag__"]["h"]
    assert d["Cd_inclination_model"] == "idelchik_3_2" and d["Cd_inclination_applied"] == 1.0
    assert d["Cd_inclination_dzeta"] > 0.0
    assert d["Cd_exit_crossflow_model"] == "rowbury_2001"
    assert d["Cd_exit_crossflow_applied"] == 1.0 and d["Cd_exit_crossflow_k_perp"] > 0.0
    assert d["Cd"] != d["Cd_base"]


@pytest.mark.parametrize("law", ["rogers_hersh_1976", "florschuetz_metzger_su_1982"])
def test_impingement_rows_report_the_jet_law_as_the_exit_term_and_alias_cd_own(law: str) -> None:
    g, arr = _impingement(law)
    diag = _solve(g)["__element_diag__"]
    rows = [diag[arr.plate_id(i)] for i in range(1, 7)]
    for p in rows:
        _assert_full_group(p)
        assert p["Cd_exit_crossflow_model"] == law
        assert p["Cd_correlation"] == "fixed" and p["Cd_base"] == 0.79
        assert p["Cd_own"] == p["Cd_base"]  # deprecated alias, until 0.13.0
        acts = p["Cd_exit_crossflow_applied"] == 1.0
        assert (p["Cd"] != p["Cd_base"]) == acts
        assert p["Cd_exit_crossflow_ratio"] == (p["Cd"] / p["Cd_base"] if acts else 1.0)
    # Rogers & Hersh is capped at the hole's own Cd on this mild crossflow;
    # Florschuetz's law replaces it on every row.
    applied = [p["Cd_exit_crossflow_applied"] for p in rows]
    assert applied == [float(law == "florschuetz_metzger_su_1982")] * 6


def test_duct_fed_panels_and_the_liner_summary_keep_the_011_names_as_aliases() -> None:
    for spec in (None, HoleCd(inlet_crossflow="none")):
        g, liner = _liner(**({} if spec is None else {"hole_cd": spec}))
        diag = _solve(g)["__element_diag__"]
        for i in (1, 2, 3):
            p = diag[liner.panel_id(i)]
            _assert_full_group(p)
            assert p["Cd_inlet_crossflow_U1_over_Vi"] > 0.0
            assert p["Cd_inlet_crossflow_applied"] == float(spec is None)
            assert p["Cd_inlet_crossflow_ignored"] == float(spec is not None)
            for old, new in HOLE_CD_DIAGNOSTIC_ALIASES.items():
                if old != "Cd_own":
                    assert p[old] == p[new], old
            assert "Cd_own" not in p
        s = liner.summarize(diag)
        for old, new in LINER_CD_DIAGNOSTIC_ALIASES.items():
            assert s[old] == s[new], old
        assert s["Cd_inlet_crossflow_U1_over_Vi_max"] == max(
            s["rows_Cd_inlet_crossflow_U1_over_Vi"]
        )


def test_a_plenum_fed_plate_emits_only_the_ignored_alias() -> None:
    # 0.11 gave coolant_crossflow_ignored on every effusion plate, and U1/Vi
    # with its bounds on a duct-fed one only: the aliases follow that.
    d = _solve(_chamber(_plate()))["__element_diag__"]["h"]
    assert d["coolant_crossflow_ignored"] == d["Cd_inlet_crossflow_ignored"] == 0.0
    assert not {"U1_over_Vi", "crossflow_cd_beyond_8pct", "crossflow_cd_degraded"} & set(d)


def test_the_model_is_the_configuration_even_where_the_term_reads_zero() -> None:
    # An explicit idelchik_3_2 on a normal (90 deg) hole: configured, but the
    # term is exactly zero there.
    d = _solve(_chamber(_plate(angle_deg=90.0, hole_cd=HoleCd(inclination="idelchik_3_2"))))
    p = d["__element_diag__"]["h"]
    assert p["Cd_inclination_model"] == "idelchik_3_2"
    assert p["Cd_inclination_applied"] == 0.0 and p["Cd_inclination_dzeta"] == 0.0
    auto = _solve(_chamber(_plate(angle_deg=90.0)))["__element_diag__"]["h"]
    assert auto["Cd_inclination_model"] == "none"


def test_wall_heat_transfer_keeps_the_011_ignored_key() -> None:
    # 0.11's public wall_heat_transfer() returned coolant_crossflow_ignored;
    # it still does until 0.13.0, and the diagnostics carry it from there.
    g = _chamber(_plate())
    g.resolve_all_topology()
    s = NetworkSolver(g)
    x0 = s._build_x0()
    s._propagate_states(x0)
    el = g.elements["h"]
    si = s._report_state(g.nodes[el.from_node], x0)
    so = s._report_state(g.nodes[el.to_node], x0)
    assert el.wall_heat_transfer(si, so)["coolant_crossflow_ignored"] == 0.0
    assert el.diagnostics(si, so)["coolant_crossflow_ignored"] == 0.0
