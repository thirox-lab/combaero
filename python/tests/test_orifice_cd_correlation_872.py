"""Cd belongs to correlation='fixed' alone (#872, 0.8.0).

Until 0.8.0 an orifice-family element took ``Cd=0.6`` by default and, under
any correlation but 'fixed' (the default was ReaderHarrisGallagher), used it
only as a start guess: ``OrificeElement(..., Cd=0.8)`` ran the correlation's
Cd, about 0.6, without a word. Now Cd defaults to None, a Cd given with a
correlation raises, and 'fixed' without a Cd raises (the impingement plate
keeps Florschuetz's documented 0.79). The start guess keeps the old 0.6
(0.79 for the plate), so a network that never passed a Cd with a
correlation is unchanged to the bit.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

import combaero as cb
from combaero.network import (
    EffusionPlateElement,
    FilmCoolingElement,
    FlowNetwork,
    HoleCd,
    ImpingementArray,
    ImpingementPlateElement,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
    SlotFilmElement,
)
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import NetworkGraphSchema

_REPO = Path(__file__).resolve().parents[2]
_ISO_GUI = _REPO / "gui" / "examples" / "iso5167_orifice_metering.json"

D = 0.01

# One constructor per orifice-family element, taking (correlation, **cd): the
# plain orifice through its own keywords, every other through its HoleCd.
_FAMILY = {
    "OrificeElement": lambda corr, **kw: OrificeElement(
        "o", "a", "b", diameter=D, correlation=corr, plate_thickness=3 * D, **kw
    ),
    "EffusionPlateElement": lambda corr, **kw: EffusionPlateElement(
        "o",
        "a",
        "b",
        hole_diameter=1e-3,
        wall_thickness=3e-3,
        pitch=6e-3,
        panel_area=1e-2,
        hole_cd=HoleCd(correlation=corr, **kw),
    ),
    "FilmCoolingElement": lambda corr, **kw: FilmCoolingElement(
        "o",
        "a",
        "b",
        hole_diameter=1e-3,
        holes_per_row=10,
        lateral_pitch=4e-3,
        wall_thickness=3e-3,
        hole_cd=HoleCd(correlation=corr, **kw),
    ),
    "SlotFilmElement": lambda corr, **kw: SlotFilmElement(
        "o",
        "a",
        "b",
        slot_height=2e-3,
        span=0.1,
        lip_thickness=6e-3,
        hole_cd=HoleCd(correlation=corr, **kw),
    ),
    "ImpingementPlateElement": lambda corr, **kw: ImpingementPlateElement(
        "o",
        "a",
        "b",
        d_jet=2.5e-3,
        xn_d=5.0,
        yn_d=4.0,
        z_d=2.0,
        span=0.1,
        plate_thickness=2.5e-3,
        hole_cd=HoleCd(correlation=corr, **kw),
    ),
    "ImpingementArray": lambda corr, **kw: ImpingementArray(
        "o",
        n_rows=3,
        d_jet=2.5e-3,
        xn_d=5.0,
        yn_d=4.0,
        z_d=2.0,
        span=0.1,
        plate_thickness=2.5e-3,
        hole_cd=HoleCd(correlation=corr, **kw),
    ),
}


@pytest.mark.parametrize("name", sorted(_FAMILY))
def test_a_cd_with_a_correlation_is_refused(name: str) -> None:
    """The old silent ignore. Fails on the pre-#872 code, which built it."""
    with pytest.raises(ValueError, match="correlation='fixed'"):
        _FAMILY[name]("IdelchikThick", Cd=0.8)


@pytest.mark.parametrize(
    "corr", ["ReaderHarrisGallagher", "IdelchikThick", "Lichtarowicz", "McGreehanSchotsch"]
)
def test_every_orifice_correlation_refuses_a_cd(corr: str) -> None:
    with pytest.raises(ValueError, match="'o'.*Cd=0.8.*correlation='fixed'"):
        _FAMILY["OrificeElement"](corr, Cd=0.8)


def test_the_default_correlation_refuses_a_cd() -> None:
    """The issue's own case: ReaderHarrisGallagher by default, Cd=0.8 given."""
    with pytest.raises(ValueError, match="ReaderHarrisGallagher"):
        OrificeElement("nozzle", "a", "b", Cd=0.8, diameter=0.195441)


@pytest.mark.parametrize(
    "name", ["OrificeElement", "EffusionPlateElement", "FilmCoolingElement", "SlotFilmElement"]
)
def test_fixed_without_a_cd_is_refused(name: str) -> None:
    with pytest.raises(ValueError, match="correlation='fixed' needs a Cd"):
        _FAMILY[name]("fixed")


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_a_non_finite_fixed_cd_is_refused(bad: float) -> None:
    with pytest.raises(ValueError, match="must be a finite number"):
        OrificeElement("o", "a", "b", Cd=bad, diameter=0.01, correlation="fixed")


def test_the_impingement_plate_keeps_florschuetz_as_its_fixed_default() -> None:
    """Its documented default: 'fixed' at Florschuetz's measured 0.79."""
    p = _FAMILY["ImpingementPlateElement"]("fixed")
    assert p.Cd == cb.FLORSCHUETZ_1981_DEFAULT_CD
    assert _FAMILY["ImpingementPlateElement"]("fixed", Cd=0.7).Cd == 0.7
    assert _FAMILY["ImpingementPlateElement"]("IdelchikThick").Cd is None


@pytest.mark.parametrize("name", sorted(_FAMILY))
def test_a_correlation_without_a_cd_builds_with_cd_none(name: str) -> None:
    assert _FAMILY[name]("IdelchikThick").Cd is None


def _boundary(name: str, Pt: float) -> PressureBoundary:
    return PressureBoundary(name, Pt=Pt, Tt=300.0, Y=cb.species.dry_air_mass())


def _solve(el: OrificeElement) -> dict:
    g = FlowNetwork()
    g.add_node(_boundary("a", 1.05e5))
    g.add_node(_boundary("b", 1.00e5))
    g.add_element(el)
    res = NetworkSolver(g).solve()
    assert res["__success__"], res.get("__message__")
    return res["__element_diag__"][el.id]


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_the_fixed_path_uses_the_given_cd(regime: str) -> None:
    fixed = _solve(
        OrificeElement("o", "a", "b", Cd=0.8, diameter=D, correlation="fixed", regime=regime)
    )
    assert fixed["Cd"] == 0.8
    rhg = _solve(OrificeElement("o", "a", "b", diameter=D, regime=regime))
    assert rhg["Cd"] == pytest.approx(0.6, abs=0.02)  # the correlation's, not 0.8
    assert fixed["m_dot"] > 1.25 * rhg["m_dot"]


def test_the_seed_assumes_the_old_default_cd_under_a_correlation() -> None:
    """The start guess's Cd: the fixed one, else 0.6 standing in for the
    correlation's (0.79 on an impingement plate), the value the unset Cd
    had before 0.8.0. Before resolve (no geometry) _effective_Cd falls back
    to it too."""
    rho = 1.2
    o = OrificeElement("o", "a", "b", diameter=D)
    A = o.area
    assert o._seed_resistance(rho, 1e5, 1.0) == 1.0 / (2.0 * rho * (0.6 * A) ** 2)
    assert o._effective_Cd(None, None) == 0.6  # type: ignore[arg-type]
    f = OrificeElement("o", "a", "b", Cd=0.8, diameter=D, correlation="fixed")
    assert f._seed_resistance(rho, 1e5, 1.0) == 1.0 / (2.0 * rho * (0.8 * A) ** 2)
    p = _FAMILY["ImpingementPlateElement"]("IdelchikThick")
    assert p._seed_Cd == cb.FLORSCHUETZ_1981_DEFAULT_CD


def test_the_seed_does_not_move_the_solved_flow() -> None:
    """A different seed Cd converges to the same root."""

    class _OtherSeed(OrificeElement):
        _seed_Cd_default = 0.9

    base = _solve(OrificeElement("o", "a", "b", diameter=D, plate_thickness=0.0))
    other = _solve(_OtherSeed("o", "a", "b", diameter=D, plate_thickness=0.0))
    assert other["m_dot"] == pytest.approx(base["m_dot"], rel=1e-8)


def _iso_with_cd(cd: float | None) -> dict:
    """The ISO example as a pre-hole_cd (v0.11.1) file: the correlation and
    the Cd loose on the node, lifted into hole_cd on load."""
    data = json.loads(_ISO_GUI.read_text())
    plate = next(n for n in data["nodes"] if n["id"] == "plate")
    plate["data"]["correlation"] = plate["data"].pop("hole_cd")["correlation"]
    if cd is not None:
        plate["data"]["Cd"] = cd
    return data


def test_an_old_gui_file_with_a_cd_and_a_correlation_still_loads() -> None:
    """The shipped ISO example carried Cd 0.6 beside ReaderHarrisGallagher
    until #872 (the Inspector shows the field for 'fixed' only, so every
    saved orifice has one). It builds with Cd None and solves exactly as the
    file without it: the Cd was never used and still is not."""
    old = _iso_with_cd(0.8)
    net = build_network_from_schema(NetworkGraphSchema(**copy.deepcopy(old)))
    plate = net.elements["plate"]
    assert plate.correlation == "ReaderHarrisGallagher" and plate.Cd is None
    r_old = NetworkRunner.from_dict(old).solve(timeout=120.0)
    r_new = NetworkRunner.from_dict(_iso_with_cd(None)).solve(timeout=120.0)
    assert r_old.success and r_new.success
    m_old = r_old._element_results["plate"].model_dump()["m_dot"]
    m_new = r_new._element_results["plate"].model_dump()["m_dot"]
    assert m_old == m_new


def test_a_gui_fixed_orifice_passes_its_cd() -> None:
    data = _iso_with_cd(0.7)
    plate = next(n for n in data["nodes"] if n["id"] == "plate")
    plate["data"]["correlation"] = "fixed"
    net = build_network_from_schema(NetworkGraphSchema(**data))
    assert net.elements["plate"].Cd == 0.7
