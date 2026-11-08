"""Measurement conventions: are the two sides measuring the same quantity?

The last structural gap #389 names, and the one no range check can close:
every coordinate is in range, the quantity is simply not the same one.
This project has hit it three times -- Rohde's total-referenced `Cd`
against McGreehan-Schotsch's static, Andrews' `h_m` over the hole length
against his own `h` on the plate area (worth -60%), and a gas-side `h`
quoted without its unblown baseline (worth 1.75x).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.cooling import conventions as C  # noqa: E402, N812
from validation.cooling.schema import load_dataset  # noqa: E402


def _fake_series(source_name: str, y_axis: str):
    """A stand-in for a SeriesMetadata that does not exist in the dataset.

    The refusal paths must be exercised, and the dataset is deliberately
    free of them -- so the mismatch is constructed rather than waited for.
    """
    from types import SimpleNamespace

    return SimpleNamespace(
        label=f"synthetic/{y_axis}",
        y_axis=y_axis,
        source=SimpleNamespace(name=source_name),
    )


@pytest.fixture(scope="module")
def dataset():
    return load_dataset()


def test_every_scored_series_declares_a_convention(dataset) -> None:
    """A new source cannot arrive and be scored silently.

    The declaration is a table keyed (source, y_axis) rather than a field
    repeated on 103 series, because the convention is a property of how a
    rig was instrumented and not of the individual curve. This is what
    stops the table going stale: add a scored source and this fails until
    its convention is written down.
    """
    missing = sorted(
        {
            (s.source.name, s.y_axis)
            for s in dataset
            if s.scores is not None and C.series_convention(s) is None
        }
    )
    assert not missing, (
        f"these scored (source, quantity) pairs declare no convention: "
        f"{missing}. Add them to conventions.SOURCE_PUBLISHES -- scoring a "
        "quantity whose basis is unrecorded is what #389 exists to stop."
    )

    undeclared_sets = sorted(
        {s.scores for s in dataset if s.scores is not None and s.scores not in C.SET_PRODUCES}
    )
    assert not undeclared_sets, (
        f"these sets score series but declare no output convention: {undeclared_sets}"
    )


def test_the_whole_dataset_resolves(dataset) -> None:
    """No scored series may be unregistered or incompatible today.

    Not a tautology: `resolve` returns those statuses and the scorecard
    refuses on them. This asserts the current dataset is clean, so a
    future mismatch shows up here rather than as a large error someone
    reads as a model limitation.
    """
    statuses: dict[str, list[str]] = {}
    for s in dataset:
        if s.scores is None:
            continue
        r = C.resolve(s, s.scores)
        statuses.setdefault(r.status, []).append(s.label)

    assert set(statuses) <= {C.STATUS_DIRECT, C.STATUS_CONVERTED}, {
        k: v[:3] for k, v in statuses.items()
    }
    # +4 #402, +1 #435, -2 fig4.51 Gbar 45par/45crs resolved as duplicates of
    # G, +3 Taslim & Spring square-channel friction (#403 item 2), +1 fig4.53
    # 60 deg V on its own shape (#434), +7 Taslim & Spring Nu series on
    # their ratio sets (#444), +2 Colban 2011 Figs. 6-7 shaped-hole film (#386)
    # +5 Taslim & Spring AR 0.5 / 3.5 two-side friction, every eligible set.
    # +1 Perry 1949 air via Jobson Fig. 5a, K on the reservoir state (#854).
    # +1 Schiller 1933 steam via Jobson Fig. 5b, same convention (#861).
    # +9 Behrendt 2008 eta_tot through EffusionPlateElement (eta_overall, direct).
    assert len(statuses[C.STATUS_DIRECT]) == 132
    # +12 CR-4396 Figs. 6.1/6.2, Kays-Perkins to Dittus-Boelter (#452).
    cr4396 = [
        f"johnson1993/fig{fig}_{wall}_seg{seg}"
        for fig in ("6.1", "6.2")
        for wall in ("leading", "trailing")
        for seg in "BCD"
    ]
    assert sorted(statuses[C.STATUS_CONVERTED]) == sorted(
        [
            *cr4396,
            "rohde1969/fig10_rd0_triangles",
            "rohde1969/fig10_rd0p195_squares",
            "rohde1969/fig10_rd0p488_circles",
        ]
    )


def test_rohde_is_the_acceptance_case(dataset) -> None:
    """#389 names it as the acceptance case because the answer is known.

    Rohde reports `Cd` on duct TOTAL pressure; McGreehan and Schotsch's is
    on STATIC. `orifice_runner` has always applied
    `sqrt(VHR/(VHR-1))` correctly -- the gap was that nothing DECLARED a
    conversion was needed, so nothing would notice a series arriving the
    other way round.
    """
    rohde = next(s for s in dataset if s.label == "rohde1969/fig10_rd0_triangles")
    assert C.series_convention(rohde) == "Cd_total"
    assert C.SET_PRODUCES[rohde.scores] == "Cd_static"

    r = C.resolve(rohde, rohde.scores)
    assert r.status == C.STATUS_CONVERTED
    assert r.conversion is not None
    assert "sqrt(VHR/(VHR-1))" in r.conversion.derivation
    assert "orifice_runner" in r.conversion.applied_by

    # The registered factor is the one the runner uses, not a copy of it.
    from validation.cooling.orifice_runner import vhr_to_static_cd

    for vhr in (2.0, 10.0, 50.0):
        assert r.conversion.factor(vhr) == pytest.approx(math.sqrt(vhr / (vhr - 1.0)))
        assert vhr_to_static_cd(vhr, 0.6) == pytest.approx(0.6 * r.conversion.factor(vhr))
    # The size that makes this a category error rather than a disagreement.
    assert r.conversion.factor(2.0) == pytest.approx(1.414, abs=0.01)
    assert r.conversion.factor(50.0) == pytest.approx(1.010, abs=0.01)


def test_a_mislabelled_rohde_series_is_caught_not_converted_twice(dataset) -> None:
    """The guard the registry buys, exercised directly.

    Before this, the conversion fired on the AXIS NAME. A Rohde-shaped
    series already on the static basis would have been converted a second
    time -- 41% high at VHR 2 -- and nothing would have said so.
    """

    from validation.cooling import orifice_runner

    rohde = next(s for s in dataset if s.label == "rohde1969/fig10_rd0_triangles")
    assert orifice_runner.run_series(rohde), "the series scores today"

    try:
        C.SERIES_OVERRIDE[rohde.label] = "Cd_static"
        assert C.resolve(rohde, rohde.scores).status == C.STATUS_DIRECT
        with pytest.raises(ValueError, match="applied twice"):
            orifice_runner.run_series(rohde)
    finally:
        C.SERIES_OVERRIDE.pop(rohde.label, None)

    # And restored, so the suite order cannot matter.
    assert C.resolve(rohde, rohde.scores).status == C.STATUS_CONVERTED


def test_adiabatic_and_overall_effectiveness_are_declared_incompatible() -> None:
    """The refusal that matters most, and it is not hypothetical.

    An adiabatic wall passes no heat, so its effectiveness carries no
    coefficient; an overall effectiveness contains the internal convection
    AND the gas side. Converting one to the other needs the whole
    resistance network -- a model, not a conversion -- and an overall-eta
    correlation used as a closure would double-count the convection the
    network computes (#387).
    """
    why = C.INCOMPATIBLE[("eta_adiabatic", "eta_overall")]
    assert "double-count" in why

    r = C.resolve(_fake_series("andrews1988", "eta_overall"), "baldauf_2002_sellers")
    assert r.status == C.STATUS_INCOMPATIBLE
    assert not r.scorable


def test_an_unregistered_mismatch_is_refused_not_scored() -> None:
    """The default is refusal, which is the whole design.

    A pair nobody thought about must not score. Registering a conversion
    is a deliberate act with a derivation attached.
    """

    r = C.resolve(_fake_series("rohde1969", "Cd"), "florschuetz_1981_inline")
    assert r.status == C.STATUS_UNREGISTERED
    assert not r.scorable
    assert "category error" in r.detail


def test_a_registered_conversion_is_definition_not_a_fit() -> None:
    """Every conversion must be geometry or algebra, never a correction.

    A conversion that needed fitting would be a model and would belong in
    the library with its own evidence. Han's `G_bar = 1.2 G` is the case
    that tests this: it looks like a conversion and is NOT registered as
    one, because the measured ratio spans 1.096 to 1.413 with 69% of the
    variance between rigs (#401).
    """
    assert ("G_ribbed_wall", "G_four_wall") not in C.CONVERSIONS
    assert ("G_ribbed_wall", "G_four_wall") in C.INCOMPATIBLE
    assert "CONSTANT" in C.INCOMPATIBLE[("G_ribbed_wall", "G_four_wall")]

    for (src, dst), conv in C.CONVERSIONS.items():
        assert conv.derivation, (src, dst)
        assert conv.applied_by, (src, dst)
        # The derivation must point at where a reader can check it.
        assert "extractions/" in conv.derivation, (src, dst)


def test_the_andrews_h_pair_is_named_even_though_the_axis_already_split_it(
    dataset,
) -> None:
    """Honest about what this does and does not guard.

    `h_hole_length` and `h_internal` are distinct `y_axis` values, so the
    effusion runner already dispatches on them and `resolve` cannot catch
    a mislabel -- the axis IS the declaration. What the registry adds here
    is that both are NAMED in one place with the -60% consequence written
    down, and that the conversion between them is recorded as pure
    geometry rather than living only in a runner.
    """
    pair = C.CONVERSIONS[("h_hole_internal_area", "h_plate_approach_area")]
    assert "pi D L" in pair.derivation and "X^2" in pair.derivation
    assert pair.factor is None, (
        "the factor needs the series geometry, so the runner owns it; a "
        "factor here would be a second implementation to drift"
    )

    a86 = next(s for s in dataset if s.label.startswith("andrews1986/fig8_hm_plate_a"))
    a88 = next(s for s in dataset if s.label == "andrews1988/fig8_h_effusionC")
    assert C.series_convention(a86) == "h_hole_internal_area"
    assert C.series_convention(a88) == "h_plate_approach_area"
    # Both are produced by the one set, so both resolve direct.
    assert C.resolve(a86, a86.scores).status == C.STATUS_DIRECT
    assert C.resolve(a88, a88.scores).status == C.STATUS_DIRECT


def test_every_named_convention_is_documented_and_used() -> None:
    """A convention nobody can read is a string, not a declaration."""
    used = set(C.SOURCE_PUBLISHES.values()) | set(C.SET_PRODUCES.values())
    for extra in C.SET_ALSO_PRODUCES.values():
        used |= set(extra)
    for src, dst in list(C.CONVERSIONS) + list(C.INCOMPATIBLE):
        used |= {src, dst}

    undocumented = sorted(used - set(C.CONVENTION_NOTES))
    assert not undocumented, f"no note explains: {undocumented}"

    unused = sorted(set(C.CONVENTION_NOTES) - used)
    assert not unused, f"documented but referenced by nothing: {unused}"


def test_the_undeclared_count_does_not_grow(dataset) -> None:
    """`undeclared` scores, so it needs a ratchet or it becomes permanent.

    The design lets an undeclared pair through rather than refusing, so
    that incomplete bookkeeping cannot hide working results. That is only
    honest if the gap is counted and cannot quietly widen -- the
    `conventions.py` docstring promised this test and it was written
    after a falsification run showed nothing held the promise.

    Today the count is ZERO. If it ever rises, the right fix is to
    declare the convention, not to raise the bound.
    """
    from validation.cooling.scorecard import build, run_dataset

    cells = build(run_dataset(dataset), dataset)
    undeclared = sorted(c.label for c in cells if c.convention == C.STATUS_UNDECLARED)
    assert not undeclared, (
        f"{len(undeclared)} scored rows declare no convention: "
        f"{undeclared[:5]}. Add them to conventions.SOURCE_PUBLISHES."
    )

    # And the status really is reachable -- a ratchet pinned at zero that
    # could never be non-zero would assert nothing.
    probe = C.resolve(_fake_series("not_a_source", "not_an_axis"), "han_1988_orthogonal")
    assert probe.status == C.STATUS_UNDECLARED
    assert probe.scorable, "undeclared must still score; refusing would hide results"


def test_a_refusal_replaces_the_numbers_rather_than_sitting_beside_them(
    dataset,
) -> None:
    """How a refusal reaches the scorecard.

    A category error reported as a large error reads as a model
    limitation, so `build` zeroes the scored count and writes the reason
    instead. Exercised by pointing a real series at a set producing an
    incompatible quantity.
    """
    import dataclasses

    from validation.cooling.scorecard import build, run_dataset

    murray = next(s for s in dataset if s.label == "murray2018/fig6_eta_exp_M0p19")
    assert C.resolve(murray, murray.scores).status == C.STATUS_DIRECT

    crossed = dataclasses.replace(murray, scores="andrews_1986_effusion_internal")
    r = C.resolve(crossed, crossed.scores)
    assert r.status == C.STATUS_UNREGISTERED, r.status
    assert not r.scorable

    cell = next(
        c for c in build(run_dataset([crossed]), [crossed]) if c.label.startswith("murray2018/")
    )
    assert cell.n_scored == 0, "a refused series must report no score"
    assert cell.reason and r.status in cell.reason
    assert math.isnan(cell.mae) and math.isnan(cell.bias)
