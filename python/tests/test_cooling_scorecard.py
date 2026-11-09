"""The aggregate scorecard: one runner per series, one view of all of them.

Three runners exist because three physics families score differently. Before
this, `scorecard.py` knew only the rib runner, so the jet-array and orifice
series read as "not scored by any set" and there was no single place to see
whether cooling as a whole was healthy.
"""

from __future__ import annotations

import math

import pytest

from validation.cooling import jet_array_runner, orifice_runner
from validation.cooling.schema import load_dataset
from validation.cooling.scorecard import (
    OUT_OF_DOMAIN_MARK,
    Cell,
    build,
    rollup,
    run_dataset,
)


@pytest.fixture(scope="module")
def dataset():
    return load_dataset()


@pytest.fixture(scope="module")
def cells(dataset):
    return build(run_dataset(dataset), dataset)


def test_every_series_is_claimed_exactly_once(dataset) -> None:
    """The bug this guards against is double counting.

    Every runner returns reason-carrying records for series it cannot score,
    so "the runner returned something" is not an ownership test -- the first
    version of run_dataset used it and the jet-array runner swallowed the
    orifice series, reporting them unscored.
    """
    records = run_dataset(dataset)
    per_series: dict[str, set[str]] = {}
    for r in records:
        per_series.setdefault(r.series.label, set()).add(type(r).__module__)

    assert len(per_series) == len(dataset), "a series was dropped or duplicated"
    for label, modules in per_series.items():
        assert len(modules) == 1, f"{label} scored by {modules}"


def test_ownership_predicates_are_disjoint(dataset) -> None:
    """No series may be owned by two specialists."""
    for s in dataset:
        assert not (jet_array_runner.owns(s) and orifice_runner.owns(s)), s.label


def test_the_specialist_series_are_now_scored(cells) -> None:
    """Florschuetz and the orifice sources used to report 0 scored points."""
    scored = {c.label: c for c in cells if c.n_scored > 0}

    assert any(c.startswith("florschuetz1981/") for c in scored)
    assert any(c.startswith("rohde1969/") for c in scored)
    assert any(c.startswith("mcgreehan_schotsch1988/") for c in scored)


def test_rohde_is_reported_per_velocity_head_band(cells) -> None:
    """Pooling Rohde's bands would mix a near-exact comparison with one
    dominated by the total-to-static conversion (factor 1.01 vs 1.33), so the
    runner splits them and the scorecard must keep them split."""
    rohde = [c.label for c in cells if c.label.startswith("rohde1969/fig10_rd0_")]
    assert len(rohde) == 3, rohde
    assert all("[VHR " in lbl for lbl in rohde)

    # And the split is not cosmetic: the bands disagree by a wide margin.
    by_band = {
        c.label.split("[")[1]: c for c in cells if c.label.startswith("rohde1969/fig10_rd0_")
    }
    low = by_band["VHR 1-10]"].bias
    high = by_band["VHR 25-60]"].bias
    assert low > 3.0 * high, f"bands no longer separate: {low:.3f} vs {high:.3f}"


def test_rollup_covers_every_implemented_set(cells) -> None:
    sets = {c.scored_by for c in rollup(cells)}
    for expected in (
        "han_1988_orthogonal",
        "han_park_1988_angled",
        "rallabandi_2009_high_re",
        "florschuetz_1981_inline",
        "mcgreehan_schotsch_1988_cd",
        "mcgreehan_schotsch_1988_crossflow_cd",
    ):
        assert expected in sets, f"{expected} missing from the rollup"


def test_rollup_matches_the_runners_own_numbers(dataset) -> None:
    """The aggregate must not quietly disagree with the runner it aggregates.

    Scored independently here, through each runner's own API, and compared
    against what the rollup reports.
    """
    rows = rollup(build(run_dataset(dataset), dataset))

    def only(set_name: str):
        """The one rollup row for a set.

        Rollup rows are segregated by sampling completeness, so a set with
        both completely- and partially-sampled series has several. The two
        cross-checked here have one apiece; assert that rather than assume
        it, so this stops being a valid comparison loudly if it changes.
        """
        matching = [c for c in rows if c.scored_by == set_name]
        assert len(matching) == 1, (
            f"{set_name} now spans {len(matching)} sampling classes; comparing "
            "it against a single whole-set number from the runner is no "
            "longer like-for-like"
        )
        return matching[0]

    flor = [s for s in dataset if s.source.name == "florschuetz1981"]
    recs = [r for r in jet_array_runner.run_all(flor) if r.predicted is not None]
    errs = [r.rel_error for r in recs]
    direct_bias = sum(errs) / len(errs)
    assert only("florschuetz_1981_inline").bias == pytest.approx(direct_bias, abs=1e-9)
    assert only("florschuetz_1981_inline").n_scored == len(errs)

    ms = [s for s in dataset if s.source.name == "mcgreehan_schotsch1988"]
    sc = orifice_runner.score(orifice_runner.run_all(ms))
    assert only("mcgreehan_schotsch_1988_crossflow_cd").bias == pytest.approx(sc.bias, abs=1e-9)


def test_a_set_that_answered_nothing_cannot_look_perfect(cells) -> None:
    """The shape of mistake this harness exists to catch: an unscored series
    reporting 0.0% error rather than '-'."""
    for c in cells:
        if c.n_scored == 0:
            assert math.isnan(c.mae) and math.isnan(c.rmse) and math.isnan(c.bias)


def test_within_is_blank_where_no_band_is_stated(cells) -> None:
    """`uncertainty` is the band agreement is judged against. Han and
    Florschuetz state one (3-8%); Rohde and McGreehan-Schotsch state none, so
    their `within` must read '-' rather than a structural 0%.
    """
    for c in cells:
        if c.label.startswith(("rohde1969/", "mcgreehan_schotsch1988/")):
            assert math.isnan(c.within), f"{c.label} reports a within fraction"
        if c.label.startswith("han2012/fig4.46_R_eD") and c.n_scored:
            assert not math.isnan(c.within)


# ---------------------------------------------------------------------------
# Lau 1990 -- the first rib source outside Han (#385)
# ---------------------------------------------------------------------------


def _by_series(cells, label):
    """Every cell for one series, across the in/out-of-domain split.

    `build` splits a series whose points straddle the correlation set's own
    validity box (#389), so a lookup by exact label now finds only the
    in-domain half. Tests that mean "the whole series" say so.
    """
    return [c for c in cells if c.label.split(f"  {OUT_OF_DOMAIN_MARK}")[0] == label]


def _pooled(cells, label):
    """The whole series as one cell-like tuple: (n, n_scored)."""
    cs = _by_series(cells, label)
    assert cs, f"no cells for {label}"
    return sum(c.n for c in cs), sum(c.n_scored for c in cs)


def test_lau_independently_confirms_han_G(cells) -> None:
    """A different lab, rig and decade reproducing han_1988_orthogonal's
    heat-transfer roughness function.

    Han prints G = 3.7 (e+)^0.28, Lau 4.218 (e+)^0.257 -- different
    coefficient AND exponent, agreeing because the forms cross near e+ = 300.
    """
    n, n_scored = _pooled(cells, "lau1990/fig_table2_90deg_G")
    assert n_scored == 9

    # All nine are inside Han's declared box. The top point, e+ = 636, maps
    # to Re 60024 at this geometry: 0.04% over Han's stated 60000 ceiling,
    # which is that ceiling quoted to its own precision rather than a higher
    # Reynolds number. It was split out as out-of-domain until bounds were
    # read with RibRange::rel_tol.
    in_domain = next(c for c in cells if c.label == "lau1990/fig_table2_90deg_G")
    assert in_domain.n_scored == 9
    assert in_domain.domain == "in-domain"
    assert not any(c.label == f"lau1990/fig_table2_90deg_G  {OUT_OF_DOMAIN_MARK}" for c in cells)

    assert in_domain.rmse < 0.03, f"agreement degraded to {in_domain.rmse:.4f}"
    assert abs(in_domain.bias) < 0.02
    # Every point inside Lau's own stated +/-5.8% Stanton uncertainty
    # -- a MEASUREMENT band from the source, not a model-derived one.
    assert in_domain.within == 1.0


def test_lau_R_is_scored_now_that_an_absolute_R_branch_exists(cells) -> None:
    """This asserted the opposite until the e+ path gained an absolute-R
    branch: the series was committed `scores: null` so the ~13%
    disagreement stayed visible rather than absent, because scoring it
    would have compared R against G.

    The disagreement itself is pinned in
    test_lau_R_disagreement_is_now_scored_not_just_narrated."""
    n, n_scored = _pooled(cells, "lau1990/fig_table2_90deg_R")
    assert n == 9
    assert n_scored == 9
    assert all(not math.isnan(c.rmse) for c in _by_series(cells, "lau1990/fig_table2_90deg_R"))


def test_lau_Gbar_independently_confirms_hans_1p2_factor(cells) -> None:
    """Han's `G_bar = 1.2 G` rests on a single printed label on figure 4.47.
    Lau's Eq. (8) Gbar is the same construction -- St_avg over two ribbed and
    two smooth walls -- so it is a like-for-like check by another lab.

    That Han's Nu(AV) really is that average is verified, not assumed: NASA
    CR-3837's appendix prints Nu(R), Nu(S) and Nu(AV) per run, and Nu(AV) is
    their two-way mean to 0.084% over 33 rows.

    An earlier version of this test asserted the opposite -- that tagging
    Lau's Gbar would apply a "Prandtl factor" -- by quoting the WITHDRAWN
    form of han_ribbed.md item 10. See #392.
    """
    n, n_scored = _pooled(cells, "lau1990/fig_table2_90deg_Gbar")
    assert n == 9
    assert n_scored == 9

    # Two labs' wall-averaging ratios differ by 3-4% (Lau 1.235-1.251 against
    # Han's constant 1.200), so this confirms the relationship without
    # reproducing it exactly. Pinned as a band, not a point. All nine are
    # in-domain: the top one sits on Han's Reynolds ceiling (Re 60024 against
    # 60000), see test_lau_independently_confirms_han_G.
    c = next(c for c in cells if c.label == "lau1990/fig_table2_90deg_Gbar")
    assert c.n_scored == 9
    assert 0.02 < c.rmse < 0.06, f"RMSE {c.rmse:.3f} outside the recorded band"
    assert c.bias < 0.0, "Han's 1.2 should sit below Lau's ratio, not above"
    assert c.within >= 0.85


def test_gbar_series_score_at_every_angle_and_are_labelled_accuracy(dataset, cells) -> None:
    """`G_bar = 1.2 G` applied as published, with the cost reported.

    This asserted the opposite until 2026-09-26 -- that off-90 G_bar
    series must come back unscored. Refusing withheld a number Han
    publishes; the ratio's real variability is an accuracy limit, and
    accuracy limits get measured, not hidden.

    `fig4.51_Gbar_60par` is the case that moved: 60 deg parallel ribs,
    scored against han_park_1988_angled at MAE 4.9%, which is the size of
    error the 1.2 approximation actually costs there.
    """
    by_label = {c.label.split("  [")[0]: c for c in cells}
    scored_off_90 = [
        s
        for s in dataset
        if s.y_axis == "G_bar"
        and s.alpha_deg not in (None, 90.0)
        and (by_label.get(s.label) or Cell("", "", None)).n_scored
    ]
    assert scored_off_90, (
        "no off-90 G_bar series scores; the published constant is being withheld again"
    )
    cell = by_label["han2012/fig4.51_Gbar_60par"]
    assert cell.n_scored == cell.n
    assert cell.basis == "accuracy", (
        "figure 4.51 is Han, Zhang and Lee (1991), a different paper from "
        "the set scoring it, so this is an accuracy result"
    )
    assert 0.02 < cell.mae < 0.10


def test_rollup_never_pools_across_sampling_completeness(dataset) -> None:
    """A tabulated series and an overplotted one do not measure the same
    thing, so the rollup must not average them into one row.

    Where a figure overplots symbol classes, only the spatially isolated
    marks can be digitised, and those are the ones furthest from the
    cluster centre. So a `partial` row's error is an upper bound, and a
    `complete` row's is an estimate. One number covering both is neither.
    """
    summary = rollup(build(run_dataset(dataset), dataset))
    han = {c.sampling: c for c in summary if c.scored_by == "han_1988_orthogonal"}
    assert "complete" in han and "partial" in han, (
        "han_1988_orthogonal must report completely- and partially-sampled "
        f"series separately; got {sorted(han)}"
    )
    assert han["complete"].label != han["partial"].label


def test_segregation_recovers_the_sources_own_stated_agreement(dataset) -> None:
    """The reason segregation matters, pinned as a number.

    Han prints that 95% of his data lies within 8% on G and 6% on R.
    Pooled, `han_1988_orthogonal` reported 73% within, which reads as a
    deficiency in the correlation. Segregated, the completely-sampled
    series reach about 91% -- close to the source's own claim -- and the
    partially-sampled ones about 65%.

    The gap is the pick, not the model. If this test fails because the
    complete subset fell, something real changed in the correlation; if
    it fails because the two converged, the sampling classification has
    stopped discriminating and is no longer worth its complexity.
    """
    # Restricted to G. The comparison has to be like-for-like in QUANTITY:
    # pooling R in drags the completely-sampled figure from 84% to 76%, not
    # because sampling got worse but because lau1990's R carries a genuine
    # ~11% disagreement between two labs. A real model disagreement in one
    # quantity would otherwise masquerade as a sampling effect -- the same
    # comparability trap #389 is about.
    cells_by_label = {c.label.split("  [")[0]: c for c in build(run_dataset(dataset), dataset)}
    buckets: dict[str, list[tuple[float, float, int]]] = {"complete": [], "partial": []}
    for series in dataset:
        cell = cells_by_label.get(series.label)
        if cell is None or cell.scored_by != "han_1988_orthogonal":
            continue
        if series.y_axis != "G" or not cell.n_scored:
            continue
        if cell.sampling in buckets:
            buckets[cell.sampling].append((cell.within, cell.mae, cell.n_scored))

    def weighted(rows, index):
        n = sum(k for *_, k in rows)
        return sum(row[index] * row[2] for row in rows) / n

    assert buckets["complete"] and buckets["partial"]
    complete_within = weighted(buckets["complete"], 0)
    partial_within = weighted(buckets["partial"], 0)

    assert complete_within > 0.80, (
        f"completely-sampled G series now agree only {complete_within:.0%} of "
        "the time, against the source's stated 95%"
    )
    assert partial_within < 0.70
    assert complete_within - partial_within > 0.15, (
        "sampling completeness no longer separates the two populations "
        f"({complete_within:.0%} against {partial_within:.0%})"
    )
    # And the upper-bound reading must hold: partial error is the larger.
    assert weighted(buckets["partial"], 1) > weighted(buckets["complete"], 1)


def test_recovered_bounds_are_model_free_and_discarded_when_loose(dataset) -> None:
    """Recovery bounds a missing mark by the OTHER marks near that
    abscissa -- never by the correlation plus its stated band, which
    would bound the observation with the model being scored.

    The price is that the envelope is only sometimes tight enough to
    constrain, and that is a property of the figure. Figure 4.46's
    classes hug one correlation, so its envelopes are usable; figure
    4.51's classes genuinely separate, so most of its are not and are
    discarded rather than counted as weak evidence.
    """
    from validation.cooling.recovery import MAX_WIDTH, recover

    by_figure: dict[str, int] = {}
    for series in dataset:
        bounds = recover(series, dataset)
        for bound in bounds:
            assert bound.lo < bound.hi
            assert bound.width <= MAX_WIDTH, (
                f"{series.label} kept a {bound.width:.0%} envelope; anything "
                f"over {MAX_WIDTH:.0%} admits nearly any prediction"
            )
            assert bound.n_marks >= 2
        if bounds:
            by_figure[series.figure] = by_figure.get(series.figure, 0) + len(bounds)

    assert by_figure.get("4.46", 0) >= 20, (
        f"figure 4.46 should recover a useful number of runs; got {by_figure.get('4.46', 0)}"
    )
    assert by_figure.get("4.51", 0) <= 3, (
        "figure 4.51's classes separate, so its envelopes should almost all "
        f"be discarded; kept {by_figure.get('4.51', 0)}"
    )


def test_absolute_R_series_are_scored_against_R_not_G(cells) -> None:
    """Sources that TABULATE the roughness functions print absolute `R`,
    not the normalised ordinate figure 4.46 plots. The e+ path had no
    branch for it, so six series sat committed with `scores: null` purely
    for want of one -- the measurements were never in doubt.

    Pinned as a band rather than a point: these are measurements, and the
    numbers should move if the correlation does.
    """
    for label in (
        "han_park_lei1984/table_a90_R",
        "han_park_lei1984/table_a45_R",
        "lau1990/fig_table2_90deg_R",
    ):
        cell = next(c for c in cells if c.label == label)
        assert cell.n_scored == cell.n > 0, f"{label} is still unscored"
        assert not math.isnan(cell.rmse)


def test_eq_4_17_is_confirmed_out_of_sample_on_R(cells) -> None:
    """CR-3837 is e/D = 0.063, a blockage Han and Park never measured, so
    scoring its R against Eq. 4.17 is a genuine out-of-sample test of the
    friction roughness function.

    It passes: inside the source's own 6.6% friction uncertainty at 90 and
    45 deg, and the 90 deg set has every point within band.
    """
    ninety = next(c for c in cells if c.label == "han_park_lei1984/table_a90_R")
    assert ninety.mae < 0.066, (
        f"Eq. 4.17 now misses CR-3837's 90 deg R by {ninety.mae:.1%}, outside "
        "the source's stated 6.6% friction uncertainty"
    )
    assert ninety.within == 1.0


def test_lau_R_disagreement_is_now_scored_not_just_narrated(cells) -> None:
    """Reviewer item L2. Lau's R runs well above `han_1988_orthogonal`'s
    constant 3.2, which converts to two labs ~12% apart in ribbed-wall
    friction for nominally the same 90 deg configuration.

    Until the absolute-R branch existed this lived only in prose. It is now
    a harness number, and the point of pinning it is that it must stay
    VISIBLE -- a silent drift to agreement would mean the correlation or the
    conversion moved, not that two labs reconciled.
    """
    cell = next(c for c in cells if c.label == "lau1990/fig_table2_90deg_R")
    assert cell.bias < -0.08, (
        f"Lau's R bias is now {cell.bias:.1%}; the recorded disagreement is "
        "about -11%, so something moved"
    )
    # Lau states +/-10.9% on friction, and the gap is larger than his own band.
    assert cell.within == 0.0


def test_unknown_y_axis_on_the_eplus_path_is_refused(dataset) -> None:
    """The branch used to end in `else: predicted = g`, so any y_axis that
    was not G_bar or R_normalised -- absolute R among them -- was silently
    scored against G. Nothing hit it only because every such series carried
    `scores: null`.

    A quantity the path cannot predict must come back unscored with a
    reason, never compared against a different quantity.
    """
    import dataclasses

    from validation.cooling.runner import run_series

    series = next(s for s in dataset if s.label == "han_park_lei1984/table_a90_R" and s.scores)
    impostor = dataclasses.replace(series, y_axis="Nu_ratio")
    records = run_series(impostor)
    assert records
    for record in records:
        assert record.predicted is None, (
            "an unpredictable y_axis was given a number; the catch-all is back"
        )
        assert "Nu_ratio" in (record.reason or "")


def test_fidelity_and_accuracy_are_never_pooled(dataset) -> None:
    """Fidelity (the correlation's own paper) and accuracy (cross-source)
    answer different questions and fail differently: a fidelity miss is our
    transcription, an accuracy miss is the model's limitation. One number
    covering both is neither.

    See docs/VALIDATION_POLICY.md.
    """
    rows = rollup(build(run_dataset(dataset), dataset))
    for set_name in ("han_1988_orthogonal", "han_park_1988_angled"):
        bases = {r.basis for r in rows if r.scored_by == set_name}
        assert {"fidelity", "accuracy"} <= bases, (
            f"{set_name} reports {sorted(bases)}; both must appear as separate rows"
        )
    for row in rows:
        assert row.basis in ("fidelity", "accuracy", "unknown")


def test_every_scored_set_declares_which_paper_it_is(dataset) -> None:
    """`basis_of` returns "unknown" for a set with no recorded origin rather
    than defaulting to fidelity, so a new set cannot quietly claim the
    stronger of the two. This asserts none is currently unknown -- i.e. the
    map has kept up with the sets."""
    from validation.cooling.scorecard import SET_ORIGIN

    scored = {c.scored_by for c in build(run_dataset(dataset), dataset) if c.n_scored}
    missing = {s for s in scored if s not in SET_ORIGIN}
    assert not missing, (
        f"these sets score series but declare no origin paper: {sorted(missing)}. "
        "Add them to SET_ORIGIN or their results cannot be read as fidelity "
        "or accuracy."
    )


def test_han_1988_fidelity_rests_on_very_few_clean_points(dataset) -> None:
    """The interaction the pooled number hid, pinned so it stays visible.

    `han_1988_orthogonal`'s own data survives only as an overplotted
    scatter, so almost all its FIDELITY evidence is sampling-biased and its
    cross-source ACCURACY looks better than its fidelity purely because the
    cross-source data is tabulated.

    If this starts failing because the clean fidelity set grew, that is
    good news and the number should be updated.
    """
    rows = rollup(build(run_dataset(dataset), dataset))
    clean = [
        r
        for r in rows
        if r.scored_by == "han_1988_orthogonal"
        and r.basis == "fidelity"
        and r.sampling == "complete"
    ]
    assert clean, "han_1988_orthogonal reports no cleanly-sampled own-data row"
    assert sum(r.n_scored for r in clean) < 30, (
        "the clean fidelity set has grown; update docs/VALIDATION_POLICY.md's "
        "worked example, which says it rests on eight points"
    )


def test_every_scored_series_basis_is_what_its_provenance_says(dataset) -> None:
    """`basis_of` matches a substring, and the haystack includes the source
    CITATION -- which is wide enough that a citation naming another paper
    could claim a false fidelity. So the label is pinned per source rather
    than trusted.

    The citation is in the haystack because a folder name cannot always
    carry the key: `andrews1986` holds 86-GT-225's Fig. 8 AND Fig. 10, and
    the set is keyed "86-GT-225" to stop it matching the 1988 study. That
    string appears only in the citation, so without it the paper's own data
    reported "accuracy" -- understating a fidelity result, the mirror of
    the error the key exists to prevent.

    The pair of Andrews figures is the case worth stating outright: same
    author, same group, same figure number, and the basis differs.
    """
    from validation.cooling.scorecard import SET_ORIGIN, basis_of

    expected = {
        # (source folder, correlation set): basis
        ("andrews1986", "andrews_1986_effusion_internal"): "fidelity",
        ("andrews1988", "andrews_1986_effusion_internal"): "accuracy",
        ("florschuetz1981", "florschuetz_1981_inline"): "fidelity",
        ("mcgreehan_schotsch1988", "mcgreehan_schotsch_1988_crossflow_cd"): "fidelity",
        # Perry's data in Jobson's figure: the theory has no data of its own.
        ("jobson1955", "jobson_1955_mass_flow"): "accuracy",
        ("rohde1969", "mcgreehan_schotsch_1988_cd"): "accuracy",
        ("murray2018", "baldauf_2002_sellers"): "accuracy",
        ("andrei2014", "baldauf_2002_sellers"): "accuracy",
        ("lau1990", "han_1988_orthogonal"): "accuracy",
        ("han_park_lei1984", "han_park_1988_angled"): "accuracy",
        # Alstom wedge passage, W/H 6.64: another lab, another geometry.
        ("naik2014", "han_park_1988_angled"): "accuracy",
    }
    seen = set()
    for series in dataset:
        if series.scores is None or series.scores not in SET_ORIGIN:
            continue
        key = (series.source.name, series.scores)
        got = basis_of(series, series.scores)
        seen.add(key)
        if key in expected:
            assert got == expected[key], f"{series.label} reports {got}"
        else:
            # Not enumerated above, but it must still resolve to one of the
            # two -- never "unknown", and never silently something new.
            assert got in ("fidelity", "accuracy"), f"{series.label}: {got}"

    missing = set(expected) - seen
    assert not missing, f"these source/set pairs no longer exist: {sorted(missing)}"


def test_the_andrews_figures_report_as_separate_rows(dataset) -> None:
    """One correlation set, four questions, and they must not pool.

    All four rows come from `andrews_1986_effusion_internal`, and every
    pair of them would be misleading to average:

      86-GT-225 Fig. 8   the correlations' OWN paper -- fidelity
      88-GT-290 Fig. 8   the same lab's later rig -- accuracy
      88-GT-290 Fig. 10  overall eta, plate C -- a different QUANTITY
      88-GT-290 Fig. 10  overall eta, plate B -- a different REGIME

    Pooling the first two would hide that the transcription is sound
    while the correlation under-predicts. Pooling the h rows with the eta
    rows would average unlike quantities. Pooling the two plates would
    give +13%, describing neither the 3% fit nor the 24% miss.
    """
    # Behrendt 2008's rows ride on the same set (the element's coolant side)
    # and are pinned by test_behrendt_rows_report_each_gas_baseline_apart.
    rows = [
        r
        for r in rollup(build(run_dataset(dataset), dataset))
        if r.scored_by == "andrews_1986_effusion_internal" and "eta_tot," not in r.label
    ]
    keyed = {
        (r.basis, r.label.split("  [")[1].rstrip("]") if r.label.count("  [") > 1 else ""): r
        for r in rows
    }
    assert len(keyed) == len(rows) == 4, f"rows: {[r.label for r in rows]}"

    h_fidelity = keyed[("fidelity", "")]
    h_accuracy = keyed[("accuracy", "")]
    eta_b = keyed[("accuracy", "overall eta, plate B")]
    eta_c = keyed[("accuracy", "overall eta, plate C")]

    assert h_fidelity.n_scored == 41
    assert h_accuracy.n_scored == 10
    assert -0.16 < h_fidelity.bias < -0.05
    assert -0.20 < h_accuracy.bias < -0.05

    assert eta_c.n_scored == 72 and eta_b.n_scored == 66
    assert 0.0 < eta_c.bias < 0.08
    assert 0.17 < eta_b.bias < 0.30
    # The two plates must stay far apart; a pooled row would be +13%.
    assert eta_b.bias - eta_c.bias > 0.15


def test_behrendt_rows_report_each_gas_baseline_apart(dataset) -> None:
    """Behrendt 2008 through EffusionPlateElement: one rollup row per gas-side
    baseline, both ACCURACY, never pooled, and the flat plate's out of its
    kernel's Re_x span (so out-of-domain) -- reported, not dropped.

    The numbers are the element at gas_augmentation 1.0; the per-curve
    table and trends are pinned in test_effusion_network_runner.py."""
    rows = {
        r.label: r for r in rollup(build(run_dataset(dataset), dataset)) if "eta_tot," in r.label
    }
    chamber = rows["andrews_1986_effusion_internal  [eta_tot, chamber gas side]  [acc]"]
    flat = rows[
        "andrews_1986_effusion_internal  [eta_tot, flat_plate gas side]  [acc] "
        f"{OUT_OF_DOMAIN_MARK}"
    ]
    chamber_rad = rows[
        "andrews_1986_effusion_internal  [eta_tot, chamber gas side + gas radiation]  [acc]"
    ]
    flat_rad = rows[
        "andrews_1986_effusion_internal  [eta_tot, flat_plate gas side + gas radiation]  [acc] "
        f"{OUT_OF_DOMAIN_MARK}"
    ]
    assert len(rows) == 4, sorted(rows)
    for r in (chamber, flat, chamber_rad, flat_rad):
        assert r.basis == "accuracy" and r.sampling == "complete"
        # Nine curves, the x/d 18-50 vertices of each.
        assert r.n_scored == 821
    assert chamber.domain == "in-domain" and flat.domain == "out-of-domain"
    # Both over-predict on every curve; the restarted boundary layer roughly
    # halves the miss, and does not close it.
    assert 0.20 < chamber.bias < 0.24
    assert 0.11 < flat.bias < 0.14
    assert chamber.within == 0.0 and flat.within == 0.0
    # Gas radiation (#1071) lowers both by 3-5 points and closes neither.
    assert 0.15 < chamber_rad.bias < chamber.bias - 0.03
    assert 0.07 < flat_rad.bias < flat.bias - 0.03


# ---------------------------------------------------------------------------
# In-domain / out-of-domain (#389): a set must not be judged on conditions
# its authors never claimed.
# ---------------------------------------------------------------------------


def test_scored_error_is_split_by_the_sets_own_validity_box(cells) -> None:
    """`extrapolated` was a column; now it is a partition.

    Every point already carried the flag -- `evaluate_rib` sets it when any
    of `Re`, `e/D`, `p/e`, `W/H`, `alpha` or `e+` falls outside the set's
    declared range -- and the scorecard totalled it into a number beside an
    error that pooled both regimes. Pooling is what this removes.
    """
    split = [c for c in cells if c.domain == "out-of-domain"]
    assert split, "nothing is out of domain; the partition scores nothing"

    for c in cells:
        assert c.domain in ("in-domain", "out-of-domain")
        # The flag and the partition must agree, or the column is lying.
        if c.domain == "out-of-domain":
            assert c.n_extrapolated == c.n, c.label
            assert c.label.endswith(OUT_OF_DOMAIN_MARK)
        else:
            assert c.n_extrapolated == 0, c.label
            assert not c.label.endswith(OUT_OF_DOMAIN_MARK)


def test_the_rollup_never_pools_across_the_validity_box(dataset) -> None:
    """The row that would otherwise mislead, named.

    `baldauf_2002_sellers` scores 138 points and every one of them is
    outside its envelope -- Andrei's s/D of 7.37 against a stated maximum
    of 5, and Murray's 5.75. A row reading "48% MAE" without saying that
    invites reading it as the model's accuracy rather than as what
    extrapolating it costs.
    """
    rows = rollup(build(run_dataset(dataset), dataset))
    for r in rows:
        assert r.domain in ("in-domain", "out-of-domain")
        assert (r.domain == "out-of-domain") == r.label.endswith(OUT_OF_DOMAIN_MARK)
        # THE INVARIANT THAT ACTUALLY BINDS THE BUCKETING, and it was
        # missing. Every row must be PURE. Without this the rollup can
        # stop keying on domain, mix both regimes into one row and label
        # it from whichever cell happened to be first -- which is exactly
        # what a falsification run did while this test stayed green.
        assert r.n_extrapolated in (0, r.n), (
            f"{r.label} mixes {r.n_extrapolated} extrapolated points into "
            f"{r.n}; the row's error describes neither regime"
        )
        assert (r.n_extrapolated == r.n) == (r.domain == "out-of-domain")

    baldauf = [r for r in rows if r.scored_by == "baldauf_2002_sellers"]
    assert len(baldauf) == 1, "baldauf now spans both domains"
    assert baldauf[0].domain == "out-of-domain"
    assert baldauf[0].n_scored == 138

    # And a set that straddles must report BOTH rows, never one.
    han = [r for r in rows if r.scored_by == "han_1988_orthogonal"]
    assert {r.domain for r in han} == {"in-domain", "out-of-domain"}


def test_out_of_domain_is_not_assumed_to_be_worse(dataset) -> None:
    """The reason the split is a REPORT and not a filter.

    `han_1988_orthogonal`'s largest fidelity row scores BETTER outside its
    declared box than inside it -- roughly 4.6% against 6.1% MAE. Han's
    stated `e+ >= 50` floor and 10,000-60,000 Reynolds ceiling are
    conservative there, and dropping or down-weighting extrapolated points
    would have thrown away the better half.

    So the harness separates them and judges neither. If this ever flips it
    is a finding, not a failure -- but it should be noticed, which is why
    it is pinned.
    """
    rows = rollup(build(run_dataset(dataset), dataset))
    pairs: dict[tuple, dict[str, float]] = {}
    for r in rows:
        key = (r.scored_by, r.sampling, r.basis)
        if r.n_scored:
            pairs.setdefault(key, {})[r.domain] = r.mae

    both = {k: v for k, v in pairs.items() if len(v) == 2}
    assert both, "no set straddles its own validity box any more"

    better_outside = [k for k, v in both.items() if v["out-of-domain"] < v["in-domain"]]
    assert better_outside, (
        "every straddling set now scores worse outside its box; the claim "
        "that extrapolation is not automatically worse needs remeasuring"
    )
    assert ("han_1988_orthogonal", "partial", "fidelity") in better_outside


def test_a_measured_accuracy_is_never_used_as_a_band(dataset) -> None:
    """The guard, and the reason `usable_as_band()` exists.

    Two of the three rib sets carry THIS PROJECT's own measured error in
    `accuracy_R`/`accuracy_G`. Judging a model against its own error is
    what #415 removed from the dataset's `uncertainty`; it must not come
    back through the correlation set.

    `within` is computed from the SERIES' `uncertainty` and from nothing
    else, which is asserted here directly rather than trusted.
    """
    import combaero as cb
    from validation.cooling.scorecard import render_set_accuracy

    measured = [
        (name, field)
        for name in ("han_park_1988_angled", "rallabandi_2009_high_re")
        for field in ("accuracy_R", "accuracy_G")
        if getattr(getattr(cb, name)(), field).provenance == cb.AccuracyProvenance.Measured
    ]
    assert measured, "the measured figures this guards have gone"

    # None of those values may appear as a band anywhere in the scoring.
    cells = build(run_dataset(dataset), dataset)
    banded = {c.label: c.within for c in cells if not math.isnan(c.within)}
    for label, within in banded.items():
        series = next(s for s in dataset if s.label == label.split("  ")[0])
        assert series.uncertainty is not None, (
            f"{label} reports a within fraction with no series band; it must "
            "be coming from somewhere it should not"
        )
        assert 0.0 <= within <= 1.0

    # And the report states the provenance rather than printing a bare
    # number a reader would take for an author claim.
    text = "\n".join(render_set_accuracy(["han_1988_orthogonal", "han_park_1988_angled"]))
    assert "STATED" in text and "measured" in text
    assert "NOT a band" in text


def test_naik2014_rows_are_accuracy_out_of_domain_and_an_upper_bound(cells) -> None:
    """Naik et al. 2014 (Alstom wedge TE passage, W/H 6.64) as the scorecard
    prints it: one row per eligible set, every one ACCURACY and OUT OF DOMAIN,
    and sampling `unknown` -- Fig. 16 overplots the angled circle with the
    chevron triangle at 60k and 145k, so 1 of the 3 marks is clean and the
    rows are an upper bound, not an estimate. `within` is judged against the
    paper's stated HTC uncertainty (9.5%), never against any model's error.
    The unscored series report 0 scored points rather than a perfect score."""
    rows = {c.label: c for c in cells if c.series_label == "naik2014/fig16_angled_tapered"}
    own = rows[f"naik2014/fig16_angled_tapered  {OUT_OF_DOMAIN_MARK}"]
    assert own.scored_by == "han_park_1988_angled"
    assert len(rows) == 4, sorted(rows)
    for c in rows.values():
        assert c.basis == "accuracy" and c.domain == "out-of-domain"
        assert c.sampling == "unknown" and c.n_scored == 3
        assert c.convention == "direct"
    # Recorded 2026-11-09: han_park_1988_angled reads LOW (only the 60k point
    # inside 9.5%); the other three read 25-45% high.
    assert own.bias == pytest.approx(-0.148, abs=0.005)
    assert own.within == pytest.approx(1.0 / 3.0)
    others = [c for c in rows.values() if c is not own]
    assert all(0.20 < c.bias < 0.50 and c.within == 0.0 for c in others)

    unscored = [c for c in cells if c.series_label.startswith("naik2014/") and c.scored_by is None]
    assert len(unscored) == 20
    assert all(c.n_scored == 0 for c in unscored)
