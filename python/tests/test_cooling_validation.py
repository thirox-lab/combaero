"""Scored validation of the cooling correlations against digitised data.

These are not self-consistency checks. Every expectation here is a number
read off a published figure, and the bands come from the source's own
stated accuracy -- not from what the implementation currently returns.
See validation/cooling/ and issue #333.
"""

from __future__ import annotations

import dataclasses
import math

import pytest

from validation.cooling.runner import Record, run_all, run_series
from validation.cooling.schema import load_dataset, load_points


@pytest.fixture(scope="module")
def dataset():
    return load_dataset()


@pytest.fixture(scope="module")
def records(dataset) -> list[Record]:
    return run_all(dataset)


def _series(dataset, stem: str):
    for s in dataset:
        if s.path.stem == stem:
            return s
    raise AssertionError(f"no series {stem}")


def test_every_series_declares_a_cross_check(dataset) -> None:
    """The harness rule: a set without a check that could fail is not data."""
    for s in dataset:
        assert s.cross_check.strip(), f"{s.label} has no cross_check"
        assert s.item and s.series, f"{s.label} is missing provenance"


def test_han_1988_reproduces_its_own_figure_454_line(dataset) -> None:
    """The implemented set against the drawn 90 deg correlation line.

    Han and Zhang (1992) state the 90 deg reference line IS Han (1988), so
    this is the set meeting a redrawing of itself in a later paper. The 6%
    band is Han's own stated accuracy for 95% of the data, not a tolerance
    fitted to what the code returns.
    """
    recs = run_series(_series(dataset, "fig4.54_90deg_continuous_Gbar_curve"))
    errs = [r.rel_error for r in recs if r.rel_error is not None]
    assert len(errs) == 7
    rmse = math.sqrt(sum(e * e for e in errs) / len(errs))
    assert rmse < 0.06, f"RMSE {rmse:.1%} exceeds Han's stated 6%"


def test_han_1988_against_figure_454_measured_points(dataset) -> None:
    """Measured scatter must be worse than the drawn line, and still in band.

    If the points scored BETTER than the line they were drawn through, the
    two series would have been swapped in the metadata.
    """
    curve = run_series(_series(dataset, "fig4.54_90deg_continuous_Gbar_curve"))
    data = run_series(_series(dataset, "fig4.54_90deg_continuous_Gbar_data"))

    def rmse(recs: list[Record]) -> float:
        errs = [r.rel_error for r in recs if r.rel_error is not None]
        return math.sqrt(sum(e * e for e in errs) / len(errs))

    assert rmse(data) > rmse(curve)
    assert rmse(data) < 0.08


def test_45deg_series_is_refused_not_scored(dataset) -> None:
    """A 90 deg-only set must decline 45 deg data rather than score it.

    Scoring it produced a 26.7% bias that reads as model error when it is
    really the wrong correlation entirely. The refusal is the feature.

    Built from a synthetic series (dataclasses.replace on a real one) rather
    than pinned to a specific filed CSV: which series happen to carry
    scores=han_1988_orthogonal at alpha=45 is metadata, and changes as
    understanding improves -- fig4.193c_old_correlation_G_curve was scored
    against han_1988_orthogonal via this exact refusal path once, then
    re-labelled scores=null once a better reason (drawn line, no
    angled-correlation set yet) was found. The mechanism under test should
    not depend on any one series still holding that combination.
    """

    series_45 = dataclasses.replace(
        next(s for s in dataset if s.y_axis == "G"),
        alpha_deg=45.0,
        scores="han_1988_orthogonal",
        geometry=None,
    )
    recs = run_series(series_45)
    assert all(r.predicted is None for r in recs)
    assert all("90-90 deg only" in (r.reason or "") for r in recs)


def test_figure_454_x_scale_is_the_only_consistent_decade(dataset) -> None:
    """Falsify the recorded x_scale: neighbouring decades must be far worse.

    The digitiser dropped the axis's `x 10^-2` multiplier. This pins the
    recovery so a future re-digitisation cannot silently reintroduce it.
    """
    series = _series(dataset, "fig4.54_45deg_vbroken_Gbar_data")
    assert series.x_scale == 1.0e4
    raw = [(p.x / series.x_scale, p.y) for p in load_points(series)]

    def rms_against_fit(scale: float) -> float:
        # The extraction's free fit for this series, G_bar = 3.242 e+^0.2882.
        errs = [3.242 * (x * scale) ** 0.2882 / y - 1.0 for x, y in raw]
        return math.sqrt(sum(e * e for e in errs) / len(errs))

    assert rms_against_fit(1.0e4) < 0.02
    for wrong in (1.0, 1.0e3, 1.0e5):
        assert rms_against_fit(wrong) > 0.40


def test_45deg_vbroken_outperforms_90deg_line(dataset) -> None:
    """Lower G is better. The text calls broken V ribs the best performers.

    Reads the two series against each other at a shared e+, which a sign or
    axis error in either digitisation would break.
    """
    v_broken = load_points(_series(dataset, "fig4.54_45deg_vbroken_Gbar_data"))
    line = load_points(_series(dataset, "fig4.54_90deg_continuous_Gbar_curve"))

    def fit(points):
        n = len(points)
        sx = sum(math.log(p.x) for p in points)
        sy = sum(math.log(p.y) for p in points)
        sxx = sum(math.log(p.x) ** 2 for p in points)
        sxy = sum(math.log(p.x) * math.log(p.y) for p in points)
        b = (n * sxy - sx * sy) / (n * sxx - sx * sx)
        return math.exp((sy - b * sx) / n), b

    c_v, n_v = fit(v_broken)
    c_l, n_l = fit(line)
    at = 500.0
    ratio = (c_v * at**n_v) / (c_l * at**n_l)
    assert 0.70 < ratio < 0.78, f"expected ~0.739 at e+=500, got {ratio:.3f}"


def test_every_digitised_series_passes_its_figure_card() -> None:
    """The figure card and the points are two independent channels.

    The card records what the printed axes and equations say, read off the
    page without reference to where the digitiser put anything. A
    disagreement means one of the two is wrong, and which one is not
    decided here -- it is raised for a human.

    Verified against injected bugs: a dropped axis multiplier (the real
    figure 4.54 defect) fails x-span, two curves swapped fails
    printed-exponent, a double-picked mark fails count and distinct.
    """
    from validation.cooling.verify import check_all

    failures = [f for f in check_all() if not f.ok]
    assert not failures, "\n".join(f"{f.series}: {f.check}: {f.detail}" for f in failures)


def test_pooled_figure_446_validates_the_implemented_set(records) -> None:
    """han_1988_orthogonal against the figure it was extracted from.

    Pooled rather than per-class: symbols on figure 4.46 overlap and a
    mark cannot always be assigned to its class, but the set carries zero
    geometry exponents in both R and G, so every class must land on one
    curve and pooling removes a dependence on labels that cannot be fully
    trusted. See the note at the head of the han2012 metadata.

    The bands are Han's own stated accuracy, not fitted to the code.
    """
    from validation.cooling.scorecard import pool

    n_r, mae_r, _, bias_r, _ = pool(records, "han2012/fig4.46_R_eD")
    assert n_r >= 60, f"only {n_r} R points pooled"
    assert abs(bias_r) < 0.02, f"R bias {bias_r:.1%} against the printed 3.2"
    assert mae_r < 0.06

    n_g, mae_g, _, bias_g, _ = pool(records, "han2012/fig4.46_G_eD")
    assert n_g >= 40, f"only {n_g} G points pooled"
    assert abs(bias_g) < 0.06, f"G bias {bias_g:.1%}"


def test_stated_6pct_behaves_as_one_sigma_not_a_95pct_bound(records) -> None:
    """Han states 'within 6% for 95% of the data'. The data says otherwise.

    About 70% of the lower-panel cloud falls inside 6%, with a measured
    standard deviation near 6%. That is 1 sigma, not a 95% bound -- which
    is what harness tolerances here are set from. Pinned because the
    alternative reading would justify a band roughly twice as wide, and a
    band twice as wide is how a real error hides.
    """
    from validation.cooling.scorecard import pool

    _, _, rmse, _, within = pool(records, "han2012/fig4.46_R_eD")
    assert 0.04 < rmse < 0.08, f"spread {rmse:.1%} is not near Han's 6%"
    assert 0.55 < within < 0.85, (
        f"{within:.0%} inside 6% -- consistent with 1 sigma (~68%), not with the stated 95%"
    )


def test_disputed_class_labels_are_declared_and_surfaced(dataset) -> None:
    """A disputed class label must stay visible, not decay into a comment.

    The four fig4.51 45-deg series returned G_bar/G below 1, which is not
    physically attainable. Resolved 2026-10-01 against the primary paper,
    Han, Zhang and Lee (1991) JHT 113, 590, Table 2: both panels carry the
    G marks, so the two G series are confirmed and the two G_bar series
    stay disputed, unscored, as a record of the duplicate
    (extractions/han_zhang_lee_1991_jht.md).

    fig4.46_G_eD0.047_pe10_wh2 was the fifth and is now resolved -- see
    test_fig446_wh2_class_is_pinned_to_its_source_runs.

    This pins three things: the series stays marked disputed, its
    cross_check still carries the specific question to ask of the book,
    and the verifier reports it on every run. Deleting any of those makes
    the uncertainty invisible, which is worse than the uncertainty.
    """
    from validation.cooling.verify import check_all

    disputed = [s for s in dataset if s.class_confidence == "disputed"]
    assert disputed, "the known disputed series has lost its marking"

    for s in disputed:
        # The route to resolution differs -- a figure may settle one, a
        # primary paper another -- but a dispute with no route recorded is
        # just an unexplained flag. A RESOLVED dispute keeps its label when
        # the resolution is "this is not the quantity it is labelled as", so
        # nothing pools it; it must then cite the evidence and stay unscored.
        if "RESOLVED" in s.cross_check:
            assert "extractions/" in s.cross_check, f"{s.label}: resolution cites no record"
            assert s.scores is None, f"{s.label}: resolved as mislabelled but still scored"
        else:
            assert "NEEDS THE" in s.cross_check, (
                f"{s.label} is disputed but records no route to resolve it"
            )

    reported = {f.series for f in check_all() if f.check == "class-label"}
    for s in disputed:
        assert s.label in reported, f"{s.label} is disputed but not surfaced"


def test_fig446_wh2_class_is_pinned_to_its_source_runs(dataset) -> None:
    """The e/D 0.047, P/e 10, W/H 2 class was disputed on the grounds that
    its two panels pair on only 2 of 4 marks. Resolved 2026-09-26 by going
    behind the figure to NASA CR-4015 appendix 7.3 (pages 141-145), which
    tabulates the five source runs with self-identifying headers, so the
    class is found by label and the overplotted cluster is never resolved.

    This pins the finding rather than the prose: every digitised mark in
    both panels must land on one of the five documented runs. If someone
    re-picks these series and a mark drifts off the run grid, that is a
    digitisation error and should fail here.
    """
    # CR-4015 p.141-145, E/D=0.047 P/E=10.0 ALPHA=90 HYD DIA=2.667 IN.
    # e+ from the printed R = 3.2 via f; the paper states R holds 95% of
    # its data within 6%, which is about +/-3% on e+, so the tolerance is
    # that plus digitisation error at the crowded right-hand edge.
    RUN_EPLUS = [80.7, 150.5, 256.2, 485.8, 512.1]
    TOL = 0.10

    seen = 0
    for series in dataset:
        if "fig4.46" not in series.label or "eD0.047_pe10_wh2" not in series.label:
            continue
        assert series.class_confidence == "confirmed", (
            f"{series.label} lost its resolution against CR-4015"
        )
        assert "CR-4015" in series.cross_check, (
            f"{series.label} no longer records how the dispute was closed"
        )
        seen += 1
        for pt in load_points(series):
            x = pt.x
            nearest = min(RUN_EPLUS, key=lambda r: abs(x / r - 1.0))
            assert abs(x / nearest - 1.0) <= TOL, (
                f"{series.label}: mark at e+={x:.1f} is "
                f"{100 * abs(x / nearest - 1):.1f}% from the nearest "
                f"CR-4015 run ({nearest}); no source run explains it"
            )
    assert seen == 2, f"expected the G and R pair, found {seen}"


def test_disputed_labels_are_refused_by_geometry_binding_sets(dataset) -> None:
    """Harmless today, refused the moment it would matter.

    han_1988_orthogonal carries zero geometry exponents, so a wrong class
    label changes nothing and the series scores normally. A set that binds
    geometry uses the label AS an input, and the runner must decline
    rather than feed a disputed one in. This guards the 4.47/4.48 work
    before it exists.
    """
    import combaero as cb
    from validation.cooling.runner import _binds_geometry

    assert not _binds_geometry(cb.han_1988_orthogonal()), (
        "han_1988_orthogonal now binds geometry; the disputed series must "
        "be resolved against the page before it can be scored again"
    )

    # A dispute must not by itself stop a series scoring. Check one that
    # is otherwise scorable: 90 degrees, so the set's valid_alpha admits
    # it, and G rather than G_bar so no ratio conversion is involved.
    from validation.cooling.runner import run_series

    scorable = [
        s
        for s in dataset
        if s.class_confidence == "disputed"
        and s.scores
        and s.alpha_deg in (None, 90.0)
        and s.y_axis == "G"
    ]
    for s in scorable:
        recs = run_series(s)
        assert any(r.predicted is not None for r in recs), (
            f"{s.label} is disputed but otherwise scorable, and was refused"
        )


def test_every_figure_redraws(tmp_path, dataset) -> None:
    """The plotter must run for every figure in the dataset.

    A visual check only helps if it still works, and it is the kind of
    tool that rots silently -- nobody notices until they need it. This
    also pins that every series carries the axis spec the plot needs, so
    a new series cannot be added without one.
    """
    pytest.importorskip("matplotlib", reason="needs the 'examples' extra")

    from validation.cooling.plot import plot_figure

    figures = {s.figure for s in dataset if s.figure}
    assert figures, "no series declares a figure"

    for figure in sorted(figures):
        out = plot_figure(figure, tmp_path / f"{figure}.png", dataset)
        assert out.exists() and out.stat().st_size > 5000, f"{figure} drew nothing"


def test_every_figure_has_a_committed_plot(dataset) -> None:
    """A figure in the dataset must have its redrawn image committed.

    The images are committed so a reviewer sees the cloud's shape in the
    diff. That only works if they exist, and the way it fails is filing a
    new figure and forgetting to regenerate -- so that case is pinned.

    This does not check the image is CURRENT; PNG bytes are not stable
    across matplotlib versions, so comparing them would fail for reasons
    that have nothing to do with the data. Regenerate with
    `uv run python -m validation.cooling.plot --all` after changing data
    or a card.
    """
    from validation.cooling.plot import plot_path

    pairs = sorted({(s.source.name, s.figure) for s in dataset if s.figure})
    for source, figure in pairs:
        img = plot_path(source, figure)
        assert img.exists(), (
            f"{source}/{figure} has no committed plot; run "
            "`uv run python -m validation.cooling.plot --all`"
        )
        assert img.stat().st_size > 5000, f"{img.name} looks empty"


def test_figure_references_are_namespaced_by_source(dataset) -> None:
    """A figure number alone is not a unique identifier.

    Two books can each have a figure 4.46. Data, cards and plots are all
    namespaced by source; a bare number is accepted only while it happens
    to be unique, and becomes an error naming the alternatives rather than
    a silent pick once it is not.
    """
    from validation.cooling.plot import resolve_figure

    source, figure = resolve_figure(dataset, "4.46")
    assert (source, figure) == ("han2012", "4.46")
    assert resolve_figure(dataset, "han2012/4.46") == ("han2012", "4.46")

    for bad in ("nope/4.46", "han2012/9.99", "9.99"):
        with pytest.raises(SystemExit):
            resolve_figure(dataset, bad)

    # Ambiguity must be refused, not resolved by luck of ordering.

    other = dataclasses.replace(
        next(s for s in dataset if s.figure == "4.46"),
        source=dataclasses.replace(
            next(s for s in dataset if s.figure == "4.46").source,
            name="someoneelse2030",
        ),
    )
    with pytest.raises(SystemExit, match="ambiguous"):
        resolve_figure([*dataset, other], "4.46")


def test_axis_specs_are_declared_for_plotting(dataset) -> None:
    """Every series needs a log/linear declaration and plot limits.

    Defaulting these would silently draw a log figure on linear axes,
    which looks plausible and is wrong -- figure 4.53 is linear while the
    rest are log-log.
    """
    for s in dataset:
        card = s.verification or {}
        assert card.get("x_axis_type") in ("log", "linear"), f"{s.label} declares no x_axis_type"
        assert card.get("y_axis_type") in ("log", "linear"), f"{s.label} declares no y_axis_type"
        assert s.figure and s.panel, f"{s.label} is not assigned to a panel"


def test_gbar_ratio_is_hans_published_constant_at_every_angle(dataset) -> None:
    """`G_bar = 1.2 G` is what Han publishes, so it is applied everywhere.

    This test asserted the OPPOSITE until 2026-09-26: that 1.2 is a 90
    degree result and applying it to an angled rib "would manufacture a
    number that looks like a measurement", so the runner refused.

    What changed is not the evidence but what to do about it. The ratio is
    indeed not constant -- pooling every closed-form measurement available
    (91-GT-3 table 3, Lau table 2, CR-3837's per-run Nu(R)/Nu(AV) split,
    16 configurations across three rigs) gives 1.096 to 1.413 at
    e+ = 300, a 29% spread. But refusing withheld a number Han does
    publish, and turned a known ACCURACY limit into a missing answer.

    Faithful implementation uses the published constant; what it costs is
    measured and labelled. See docs/VALIDATION_POLICY.md and G_BAR_OVER_G.
    """

    from validation.cooling.runner import _gbar_reason

    template = next(s for s in dataset if s.y_axis == "G_bar")
    for alpha in (30.0, 45.0, 60.0, 90.0):
        series = dataclasses.replace(template, alpha_deg=alpha)
        assert _gbar_reason(series) is None, (
            f"G_bar at {alpha} deg was refused; Han publishes 1.2 for it"
        )


def test_keeping_hans_constant_is_justified_by_the_pooled_measurements() -> None:
    """Why no "improved" G_bar/G correlation was written.

    A generalised replacement was considered and rejected on evidence, and
    this pins the evidence so the decision can be revisited rather than
    re-argued. 16 configurations, three rigs:

        91-GT-3  mean 1.160, internal spread 12%
        CR-3837  mean 1.323, internal spread  3%
        lau1990  mean 1.327, internal spread 14%

    **69% of the variance is BETWEEN rigs, not within them.** No angle or
    shape term can reach a rig offset, so a "better" correlation would be
    fitting rig identity. Against the pooled population Han's 1.2 scores
    MAE 8.3%, a defensible one-number answer to a quantity spanning 29%.

    ALL THREE RIGS ARE TEXAS A&M, which "three rigs" does not convey.
    91-GT-3 and CR-3837 are Han's own laboratory; lau1990 is Lau's, in the
    same department. That cuts both ways and both are asserted below:

      * the conclusion is STRONGER -- three rigs inside one laboratory
        tradition separating this far makes the offset rig-level, not an
        artefact of comparing institutions, so there is even less for an
        angle or shape term to reach;
      * the MAE is OPTIMISTIC -- this is not a sample of the world's rigs,
        so 29% is a LOWER BOUND on the true spread and 1.2 may do worse
        against a non-TAMU rig.

    Same monoculture as #403 item 2, where it blocks the Lau friction
    disagreement outright.
    """
    import statistics

    def ratio(g, gb, e_plus=300.0):
        return (gb[0] * e_plus ** gb[1]) / (g[0] * e_plus ** g[1])

    gt3 = [
        ratio(g, gb)
        for g, gb in (
            ((1.61, 0.42), (2.94, 0.35)),
            ((1.82, 0.41), (1.90, 0.42)),
            ((1.04, 0.48), (1.14, 0.48)),
            ((1.01, 0.47), (1.63, 0.42)),
        )
    ]
    lau = [
        ratio(g, gb)
        for g, gb in (
            ((4.218, 0.257), (5.450, 0.250)),
            ((1.819, 0.355), (2.656, 0.335)),
            ((1.299, 0.399), (2.163, 0.360)),
            ((1.685, 0.376), (2.719, 0.336)),
            ((1.983, 0.352), (3.460, 0.315)),
            ((2.739, 0.324), (4.351, 0.291)),
            ((1.992, 0.362), (3.685, 0.313)),
        )
    ]
    cr = [1.340, 1.328, 1.325, 1.301, 1.320]

    groups = (gt3, lau, cr)
    pooled = gt3 + lau + cr
    grand = statistics.mean(pooled)
    between = sum(len(g) * (statistics.mean(g) - grand) ** 2 for g in groups)
    within = sum(sum((x - statistics.mean(g)) ** 2 for x in g) for g in groups)

    assert max(pooled) / min(pooled) > 1.25, "the spread that rules out one exact constant"
    assert between / (between + within) > 0.6, (
        "the between-rig share has dropped; if most variance is now WITHIN "
        "rigs, a generalised correlation may be worth revisiting"
    )
    mae = statistics.mean(abs(1.2 / v - 1) for v in pooled)
    assert mae < 0.12, f"Han's 1.2 now misses the pooled population by {mae:.1%}"

    # The pool is three rigs from ONE institution, so this MAE is a
    # best case. If a non-TAMU source is ever added here, the assertion
    # above should be re-derived rather than widened -- the number
    # getting worse would be the finding, not a regression.
    assert len(groups) == 3, (
        "a fourth source joined the pool; if it is not Texas A&M, the "
        "'lower bound on the spread' caveat in this docstring is now "
        "testable and should be tested"
    )


def test_90_degree_gbar_still_scores(dataset) -> None:
    """The guard must not block where the ratio was actually established.

    Figure 4.46 prints both lines, 3.7 and 4.5, a ratio of 1.216, and
    figure 4.47 labels its dashed line G_bar = 1.2 G. At 90 degrees the
    conversion is evidenced and must keep working.
    """
    from validation.cooling.runner import run_series

    scored = 0
    for s in dataset:
        if s.y_axis != "G_bar" or s.alpha_deg != 90.0 or not s.scores:
            continue
        recs = run_series(s)
        if any(r.predicted is not None for r in recs):
            scored += 1

    assert scored, "the 90 degree G_bar series stopped scoring"


def test_candidate_mode_runs(tmp_path, dataset) -> None:
    """The digitising path must keep working, not just the filed path.

    check_candidate builds a SeriesMetadata by hand, so every field added
    to that dataclass has to be added here too. It broke silently when
    `figure`, `panel` and `class_confidence` arrived, because 262 checks
    covered the filed data and none covered the path actually used while
    picking points.
    """
    import yaml

    from validation.cooling.verify import check_candidate

    csv_path = tmp_path / "candidate.csv"
    csv_path.write_text("x, y\n100, 12.0\n300, 17.0\n900, 22.0\n")
    card_path = tmp_path / "card.yaml"
    card_path.write_text(
        yaml.safe_dump(
            {
                "verification": {
                    "x_axis_type": "log",
                    "y_axis_type": "log",
                    "x_ticks": [50, 100, 1000],
                    "x_multiplier": 1.0,
                    "y_ticks": [8, 40],
                    "printed_curve": {"power_law": {"C": 3.7, "n": 0.28}},
                    "curve_tolerance": 0.2,
                }
            }
        )
    )
    findings = check_candidate(csv_path, card_path)
    assert findings, "candidate mode produced no findings at all"
    assert all(f.ok for f in findings), [f"{f.check}: {f.detail}" for f in findings if not f.ok]


def test_candidate_mode_catches_a_dropped_multiplier(tmp_path) -> None:
    """The bug candidate mode exists to catch, end to end."""
    import yaml

    from validation.cooling.verify import check_candidate

    csv_path = tmp_path / "candidate.csv"
    # Two decades low: the real figure 4.54 and 4.51 defect.
    csv_path.write_text("x, y\n1.0, 12.0\n3.0, 17.0\n9.0, 22.0\n")
    card_path = tmp_path / "card.yaml"
    card_path.write_text(
        yaml.safe_dump(
            {
                "verification": {
                    "x_axis_type": "log",
                    "y_axis_type": "log",
                    "x_ticks": [1, 10],
                    "x_multiplier": 1.0e2,
                    "y_ticks": [8, 40],
                }
            }
        )
    )
    findings = check_candidate(csv_path, card_path)
    failed = [f.check for f in findings if not f.ok]
    assert "x-span" in failed, f"dropped multiplier not caught; got {failed}"


def test_figure_448_resolves_extraction_item_41(dataset) -> None:
    """The W/H=1/2 branch boundary: which side does Han's own figure use?

    Extraction item 41 found Eq. 4.19's two branches disagreeing by ~10%
    at W/H=0.5, unresolved. Figure 4.48's own drawn correlation at that
    exact W/H settles it: its exponent must sit far closer to the wide
    branch (n=0.35) than the narrow one (n=0.258).
    """
    series = next(s for s in dataset if s.path.name == "fig4.48_G_correlation_WH0.5.csv")
    pts = load_points(series)
    fitted = _power_slope_test(pts)
    assert abs(fitted - 0.35) < abs(fitted - 0.258), (
        f"fitted exponent {fitted:.4f} is closer to the narrow branch "
        "(0.258) than the wide one (0.35) -- item 41's resolution direction "
        "has reversed"
    )
    assert 0.30 < fitted < 0.36


def _power_slope_test(points) -> float:
    xs = [math.log(p.x) for p in points]
    ys = [math.log(p.y) for p in points]
    n = len(points)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(v * v for v in xs)
    sxy = sum(a * b for a, b in zip(xs, ys, strict=False))
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


def test_figure_448_R_axis_is_linear_not_log(dataset) -> None:
    """Panel (a)'s y-axis was initially misread as log; it is linear.

    The tell is the digitised frame's own picks: 1.5, 2.5, 3.5, 4.5 are
    half-integer minor ticks, which a log axis cannot produce -- its
    conventional minor ticks sit at the same mantissa (2..9) in every
    decade, never at a half-integer. Confirmed on the printed page too:
    one minor tick sits exactly centred between each major label.
    """
    for s in dataset:
        if s.figure == "4.48" and s.panel == "top":
            card = s.verification or {}
            assert card.get("y_axis_type") == "linear", (
                f"{s.label} declares y_axis_type={card.get('y_axis_type')}, "
                "but panel (a)'s R axis is linear"
            )


# Both sets that declare the passage-average friction (f_channel_two_ribbed)
# and bind 90 deg transverse ribs, so both are eligible by rule for every
# two-wall Taslim friction series.
TASLIM_FRICTION_SETS = ("han_1988_orthogonal", "han_park_1988_angled")

TASLIM_SQUARE = {
    "taslim_spring1987/fig11_ar_1.0_solid_circles_2_sides_eDh_0.083": 0.083,
    "taslim_spring1987/fig11_ar_1.0_solid_tri_2_sides_eDh_0.167": 0.167,
    "taslim_spring1987/fig11_ar_1.0_solid_half_circles_2_sides_eDh_0.250": 0.250,
}


def test_taslim_square_channel_friction_is_scored_against_han(records) -> None:
    """#403 item 2: the independent-lab friction that sits on the other side.

    Taslim & Spring (1987) -- Northeastern / GE -- measure Han's own
    square-channel configuration (two opposite walls, 90 deg in-line, P/e
    10). Han under-predicts their passage friction; in R terms their 2.77
    sits ~13% below Han's 3.2 while Lau's sits ~12% above it, so the three
    labs bracket Han rather than one of them arbitrating. Every e/D here is
    above Han's band, so every row is out of domain by construction.

    The recovered Re range also guards figure 11's x units (Re x 10^-4):
    read raw, Re would be 2-10 and the smooth-friction term would swamp f.
    """
    rows = {}
    for r in records:
        if r.series.label in TASLIM_SQUARE:
            rows.setdefault((r.series.label, r.scored_by), []).append(r)
    assert set(rows) == {(label, name) for label in TASLIM_SQUARE for name in TASLIM_FRICTION_SETS}
    for (label, name), rs in rows.items():
        assert all(r.predicted is not None for r in rs), label
        assert all(r.extrapolated for r in rs), f"{label}: e/D above Han's band"
        assert all(20_000 <= r.re_used <= 110_000 for r in rs), label
        bias = sum(r.rel_error for r in rs) / len(rs)
        assert bias < 0.0, f"{label} @{name}: Han should under-predict, bias {bias:+.1%}"

    near = rows[
        ("taslim_spring1987/fig11_ar_1.0_solid_circles_2_sides_eDh_0.083", "han_1988_orthogonal")
    ]
    bias = sum(r.rel_error for r in near) / len(near)
    assert -0.20 < bias < -0.10, f"recorded at -14.9% when #403 item 2 closed, now {bias:+.1%}"


TASLIM_OTHER_AR = (
    "taslim_spring1987/fig11_ar_0.5_solid_hex_2_sides_eDh_0.250",
    "taslim_spring1987/fig11_ar_0.5_solid_square_2_sides_eDh_0.125",
    "taslim_spring1987/fig11_ar_3.5_solid_tri_down_2_sides_eDh_0.053",
    "taslim_spring1987/fig11_ar_3.5_solid_tri_left_2_sides_eDh_0.107",
    "taslim_spring1987/fig11_ar_3.5_solid_tri_right_2_sides_eDh_0.161",
)


def test_taslim_friction_is_scored_at_every_aspect_ratio(records) -> None:
    """Not only the square channel that reads best. AR 3.5 (Han W/H 0.286,
    below his 1-4) is the only independent check of Han's friction in a
    narrow channel, and Han reads 34-44% low there, as recorded when these
    were first scored."""
    rows = {}
    for r in records:
        if r.series.label in TASLIM_OTHER_AR:
            rows.setdefault((r.series.label, r.scored_by), []).append(r)
    assert set(rows) == {
        (label, name) for label in TASLIM_OTHER_AR for name in TASLIM_FRICTION_SETS
    }
    for (label, name), rs in rows.items():
        assert all(r.predicted is not None for r in rs), label
        assert all(r.extrapolated for r in rs), f"{label} @{name}: outside the set's domain"
    for name in TASLIM_FRICTION_SETS:
        narrow = [
            r for (label, n), rs in rows.items() if n == name and "_ar_3.5_" in label for r in rs
        ]
        bias = sum(r.rel_error for r in narrow) / len(narrow)
        assert -0.45 < bias < -0.30, f"{name}: AR 3.5 recorded at about -37%, now {bias:+.1%}"


# Per-(series, set) friction bias over every two-wall Taslim series, recorded
# 2026-11-10 when han_park_1988_angled declared the passage friction. Every
# row is out of both sets' boxes (e/D above 0.078, or W/H 0.286 below 1) and
# cross-source ACCURACY (Northeastern / GE, not Texas A&M). Pinned so a change
# to either set's friction, or to _channel_fbar, shows up as a number moving.
TASLIM_FRICTION_BIAS = {
    ("fig11_ar_1.0_solid_circles_2_sides_eDh_0.083", "han_1988_orthogonal"): -0.149,
    ("fig11_ar_1.0_solid_circles_2_sides_eDh_0.083", "han_park_1988_angled"): -0.118,
    ("fig11_ar_1.0_solid_tri_2_sides_eDh_0.167", "han_1988_orthogonal"): -0.249,
    ("fig11_ar_1.0_solid_tri_2_sides_eDh_0.167", "han_park_1988_angled"): -0.205,
    ("fig11_ar_1.0_solid_half_circles_2_sides_eDh_0.250", "han_1988_orthogonal"): -0.491,
    ("fig11_ar_1.0_solid_half_circles_2_sides_eDh_0.250", "han_park_1988_angled"): -0.447,
    ("fig11_ar_0.5_solid_hex_2_sides_eDh_0.250", "han_1988_orthogonal"): -0.117,
    ("fig11_ar_0.5_solid_hex_2_sides_eDh_0.250", "han_park_1988_angled"): -0.005,
    ("fig11_ar_0.5_solid_square_2_sides_eDh_0.125", "han_1988_orthogonal"): 0.126,
    ("fig11_ar_0.5_solid_square_2_sides_eDh_0.125", "han_park_1988_angled"): 0.193,
    ("fig11_ar_3.5_solid_tri_down_2_sides_eDh_0.053", "han_1988_orthogonal"): -0.346,
    ("fig11_ar_3.5_solid_tri_down_2_sides_eDh_0.053", "han_park_1988_angled"): -0.336,
    ("fig11_ar_3.5_solid_tri_left_2_sides_eDh_0.107", "han_1988_orthogonal"): -0.343,
    ("fig11_ar_3.5_solid_tri_left_2_sides_eDh_0.107", "han_park_1988_angled"): -0.328,
    ("fig11_ar_3.5_solid_tri_right_2_sides_eDh_0.161", "han_1988_orthogonal"): -0.437,
    ("fig11_ar_3.5_solid_tri_right_2_sides_eDh_0.161", "han_park_1988_angled"): -0.420,
}


def test_taslim_friction_bias_per_set_is_pinned(records) -> None:
    """Each set keeps its own row; nothing is pooled or picked between them."""
    got: dict[tuple[str, str], list[float]] = {}
    for r in records:
        if (
            r.series.source.name == "taslim_spring1987"
            and r.series.y_axis == "f_fanning_passage"
            and r.rel_error is not None
        ):
            got.setdefault((r.series.path.stem, r.scored_by), []).append(r.rel_error)
    assert set(got) == set(TASLIM_FRICTION_BIAS)
    for key, bias in TASLIM_FRICTION_BIAS.items():
        mean = sum(got[key]) / len(got[key])
        assert mean == pytest.approx(bias, abs=0.005), f"{key}: {mean:+.3f}"


def test_taslim_friction_is_compared_as_fbar_not_four_sided_f() -> None:
    """The convention trap, falsified rather than asserted.

    Taslim's f is the passage average of a channel with two ribbed walls --
    Han's measured fbar. Read as the four-sided f_r that R is built on, Han's
    prediction would sit ~60% HIGH at e/D 0.083 instead of ~15% low.
    """
    import csv

    import combaero as cb
    from validation.cooling.runner import SMOOTH_F_COEF, SMOOTH_F_EXP
    from validation.cooling.schema import DATA_ROOT

    # DATA_ROOT, not a cwd-relative path: CI runs pytest from python/tests.
    path = DATA_ROOT / "taslim_spring1987" / "fig11_ar_1.0_solid_circles_2_sides_eDh_0.083.csv"
    pts = [(float(x) * 1.0e4, float(y)) for x, y in list(csv.reader(path.open()))[1:]]
    g = cb.RibGeometry(e_D=0.083, p_e=10.0, W_H=1.0, alpha_deg=90.0)
    s = cb.han_1988_orthogonal()
    four_sided = [cb.evaluate_rib(s, g, re).f / f - 1.0 for re, f in pts]
    channel = [
        (cb.evaluate_rib(s, g, re).f + SMOOTH_F_COEF * re**SMOOTH_F_EXP) / 2.0 / f - 1.0
        for re, f in pts
    ]
    assert sum(four_sided) / len(four_sided) > 0.40
    assert sum(channel) / len(channel) < 0.0


def test_taslim_friction_path_refuses_what_it_cannot_predict(dataset) -> None:
    import dataclasses

    from validation.cooling.runner import run_series

    one_side = next(
        s
        for s in dataset
        if s.label == "taslim_spring1987/fig11_ar_1.0_open_circles_1_side_eDh_0.083"
    )
    base = next(s for s in dataset if s.label in TASLIM_SQUARE)
    for impostor, word in (
        (dataclasses.replace(one_side, scores="han_1988_orthogonal"), "one-side"),
        (dataclasses.replace(base, y_axis="Nu"), "only predicts f_fanning_passage"),
        (dataclasses.replace(base, geometry=None), "no rig geometry"),
    ):
        recs = run_series(impostor)
        assert recs and all(r.predicted is None for r in recs)
        assert all(word in (r.reason or "") for r in recs), recs[0].reason


def test_figure_448_nothing_is_scored_against_the_wrong_paper(dataset) -> None:
    """Overlapping data must not silently score against a cross-paper gap.

    The '[Ref. 3]' W/H=1 curve on figure 4.48 sits inside
    han_1988_orthogonal's valid range, but it reproduces Eq. 4.18
    (C=2.24, n=0.35), a DIFFERENT published form from the one
    han_1988_orthogonal encodes (Eq. 4.15/16, n=0.28). Scoring it would
    report that ~15% cross-paper gap as model error.

    #402 added han_1989_narrow_channel, Eq. 4.19's own set, and wired the
    four W/H<1 raw-scatter series to it -- that is the RIGHT paper, not the
    wrong one this test guards against, so those four are now expected to
    score. Every OTHER figure 4.48 series (the drawn correlation lines and
    the W/H=1 [Ref. 3] curve/error bars) must still be unscored, regardless
    of geometry.
    """
    fig448 = [s for s in dataset if s.figure == "4.48"]
    assert fig448, "figure 4.48 is not filed"
    narrow_channel_series = {
        "fig4.48_R_WH0.25",
        "fig4.48_R_WH0.5",
        "fig4.48_G_WH0.5",
        "fig4.48_G_WH0.25",
    }
    for s in fig448:
        if s.path.stem in narrow_channel_series:
            assert s.scores == "han_1989_narrow_channel", s.label
        else:
            assert s.scores is None, s.label


def test_figure_453_90deg_scores_on_the_performance_curve(records) -> None:
    """#435: Nu_ratio against f_ratio, with Re recovered from f_ratio.

    Han and Zhang state Re 15,000-80,000 and nothing else about the rig.
    The chain's recovered Re landing inside that range is a check the
    borrowed geometry and the fbar convention could have failed.
    """
    rows = [
        r
        for r in records
        if r.series.label == "han2012/fig4.53_90deg_continuous_Nu_ratio"
        and r.scored_by == r.series.scores
    ]
    assert len(rows) == 5
    assert all(r.predicted is not None for r in rows), [r.reason for r in rows]
    assert all(15_000 <= r.re_used <= 80_000 for r in rows), [r.re_used for r in rows]
    mae = sum(abs(r.rel_error) for r in rows) / len(rows)
    assert mae < 0.10, f"MAE {mae:.1%}; recorded at 7.9% when #435 landed"


def test_figure_453_reading_f_ratio_as_four_sided_f_is_caught() -> None:
    """The convention trap, falsified rather than asserted.

    figure 4.53's f_ratio is the channel-average fbar over f_s. Read as
    evaluate_rib's four-sided f instead, the model's ratio is ~1.8x the
    measured one across the whole stated Re range, so no Re reaches it.
    """
    import combaero as cb
    from validation.cooling.runner import SMOOTH_F_COEF, SMOOTH_F_EXP

    s = cb.han_1988_orthogonal()
    g = cb.RibGeometry(e_D=0.0625, p_e=10.0, W_H=1.0, alpha_deg=90.0)
    for re in (15_000, 80_000):
        f4 = cb.evaluate_rib(s, g, re).f / (SMOOTH_F_COEF * re**SMOOTH_F_EXP)
        assert f4 > 8.0, "four-sided ratio overlaps the measured 5.4-6.3"


def test_figure_453_60deg_vee_scores_against_its_own_shape(dataset, records) -> None:
    """#434: the 60 deg V series scores against the 60 deg forward-V row of
    Han, Zhang and Lee (1991) Table 2 -- never a parallel set, which would
    alias a V onto a parallel rib. Through the full f -> Nu chain on the
    performance-curve path it lands at ~2.5% MAE."""
    series = next(s for s in dataset if s.label == "han2012/fig4.53_60deg_vcontinuous_Nu_ratio")
    assert series.scores == "han_zhang_lee_1991_60vee"
    assert series.rib_shape == "v"
    rows = [r for r in records if r.series.label == series.label]
    assert rows and all(r.predicted is not None for r in rows)
    mae = sum(abs(r.rel_error) for r in rows) / len(rows)
    assert mae < 0.05, f"recorded at 2.5% when #434 landed, now {mae:.1%}"


def test_rib_shape_is_bound_not_aliased(dataset) -> None:
    """A set binds its rib shape. A V series offered to a parallel set, a
    series with no recorded shape offered to a crossed set, and a broken-V
    series offered to a continuous-V set are each refused for the SHAPE --
    the reason that used to surface only incidentally through rib angle."""
    import dataclasses

    from validation.cooling.runner import run_series

    vee = next(s for s in dataset if s.label == "han2012/fig4.51_G_60vee")
    crs = next(s for s in dataset if s.label == "han2012/fig4.51_G_60crs")
    broken = next(s for s in dataset if s.label == "han2012/fig4.54_45deg_vbroken_Gbar_data")
    cases = (
        (dataclasses.replace(vee, scores="han_park_1988_angled"), "binds Parallel"),
        (dataclasses.replace(crs, rib_shape=None), "records no rib shape"),
        (
            dataclasses.replace(
                broken,
                scores="han_zhang_lee_1991_45vee",
                geometry={"e_D": 0.0625, "p_e": 10.0, "W_H": 1.0},
            ),
            "v_broken ribs",
        ),
    )
    for impostor, word in cases:
        recs = run_series(impostor)
        assert recs and all(r.predicted is None for r in recs), impostor.label
        assert all(word in (r.reason or "") for r in recs), recs[0].reason


def test_printed_g_bar_is_used_not_the_1_2_proxy(records) -> None:
    """Table 2 prints G_bar per configuration (G_bar/G 1.07-1.26). Where a set
    carries it, a G_bar series is scored against THAT, not Han's 1.2 G:
    falsifiable, since at 60 deg V the printed ratio is ~1.16, not 1.2."""
    import combaero as cb

    s = cb.han_zhang_lee_1991(cb.RibShape.V, 60.0)
    rows = [r for r in records if r.series.label == "han2012/fig4.51_Gbar_60vee"]
    assert rows
    for r in rows:
        g = cb.RibGeometry(e_D=0.0625, p_e=10.0, W_H=1.0, alpha_deg=60.0)
        res = cb.evaluate_rib(s, g, r.re_used)
        assert res.has_G_bar
        assert abs(r.predicted - res.G_bar) < 1e-9 * res.G_bar
        assert abs(r.predicted - 1.2 * res.G) > 0.01 * res.G


def test_f_ratio_path_refuses_what_it_cannot_predict(dataset) -> None:
    import dataclasses

    from validation.cooling.runner import run_series

    series = next(s for s in dataset if s.label == "han2012/fig4.53_90deg_continuous_Nu_ratio")
    for impostor, word in (
        (dataclasses.replace(series, y_axis="G"), "only predicts Nu_ratio"),
        (dataclasses.replace(series, geometry=None), "no rig geometry"),
    ):
        recs = run_series(impostor)
        assert all(r.predicted is None for r in recs)
        assert all(word in (r.reason or "") for r in recs)


def test_han_park_1988_angled_scores_its_own_source(records) -> None:
    """han_park_1988_angled must reproduce figure 4.47, its own source.

    Both R and G panels scored through evaluate_rib with a representative
    geometry (the cloud carries no legend). Bands are what this project
    measured, not an author claim -- see the C++ set's accuracy_R/accuracy_G
    comment for why that distinction matters here specifically.

    THE COUNTS ARE POOLED ACROSS THE VALIDITY SPLIT, and that is itself
    worth knowing: `build` now separates in- and out-of-domain points
    (#389), and part of each panel falls outside han_park_1988_angled's
    own declared box. So the `measured` accuracy figures the C++ set
    carries -- 10.5% and 8.8%, taken from these very panels -- mix the two
    regimes. One more reason they are reported rather than used as a band.
    """
    from validation.cooling.scorecard import build

    cells = build(records)

    def panel(label):
        cs = [c for c in cells if c.label.split("  {out}")[0] == label]
        assert cs, label
        n_scored = sum(c.n_scored for c in cs)
        worst = max(c.rmse for c in cs)
        return n_scored, worst, len(cs)

    r_scored, r_rmse, r_parts = panel("han2012/fig4.47_R_vs_alpha")
    g_scored, g_rmse, g_parts = panel("han2012/fig4.47_G_vs_eplus")

    assert r_scored == 39
    assert r_rmse < 0.13
    assert g_scored == 115
    assert g_rmse < 0.11
    # Both panels really do straddle the box; if one stopped, the note
    # above about the measured figures mixing regimes would need redoing.
    assert r_parts == 2 and g_parts == 2


def test_figure_451_carries_its_stated_rig(dataset) -> None:
    """Han, Zhang and Lee (1991) p. 590 states the rig: square, e/D 0.0625,
    P/e 10. Without it the runner probes the set's validity midpoint (P/e 15
    for han_park_1988_angled), which reads Eq. 4.18 ~4% HIGH through its
    (P/e/10)^0.1 term and flatters the model -- a regression that no score
    threshold would catch, since the score gets better."""
    stated = {"e_D": 0.0625, "p_e": 10.0, "W_H": 1.0}
    measured = [s for s in dataset if s.figure == "4.51" and s.kind == "measured"]
    assert len(measured) == 18
    for s in measured:
        assert s.geometry == stated, f"{s.label}: {s.geometry}"


def test_han_park_1988_angled_confirmed_independently_by_figure_451(
    records,
) -> None:
    """A different paper's parallel-rib data, scored against this equation.

    Figure 4.51 (Han et al. 1991) is not han_park_1988_angled's own source
    (Han and Park 1988) -- an independent cross-paper check, on the two
    classes (45/60 deg PARALLEL) that match the rib shape Eq. 4.17/4.18
    was fitted to. The other seven shapes on figure 4.51 (crossed, V,
    lambda) are deliberately not scored against this set: it has no way to
    represent rib shape beyond angle, and scoring a shape it cannot model
    would repeat the mistake extraction item 23 already recorded once.
    """
    from validation.cooling.scorecard import build

    cells = {c.label: c for c in build(records)}
    scored_parallel = [
        c
        for name, c in cells.items()
        if "par" in name
        and "4.51" in name
        and not c.unsupported
        and c.scored_by == "han_park_1988_angled"
    ]
    assert scored_parallel, "no figure 4.51 parallel class scored"

    # Both parallel classes read LOW -- the recorded Eq. 4.18 accuracy limit
    # (CR-3837: -0.8% to -19.9%, han_ribbed.md "Recorded as a known
    # accuracy limit"), now in a second independent dataset. Scored at the
    # rig Han, Zhang and Lee (1991) states (e/D 0.0625, P/e 10); the earlier
    # probe P/e 15 read the model ~4% high and hid part of it. An ACCURACY
    # result, pinned as one rather than asserted away:
    #   * every parallel cell, in or out of domain, keeps the low sign;
    #   * the IN-DOMAIN cells -- what the set claims -- carry a ratchet:
    #     60 deg ~10% RMS against 15%, 45 deg ~17% against 20%.
    # Out-of-domain cells are single extrapolated points (Re above 60k);
    # they are reported, and only their sign is held.
    assert any("fig4.51_G_45par" in c.label for c in scored_parallel), (
        "45 deg parallel G should score now that its label is resolved"
    )
    bounds = {"60par": 0.15, "45par": 0.20}
    for c in scored_parallel:
        assert c.bias < 0.0, f"{c.label}: bias {c.bias:+.1%}, expected the known low reading"
        if c.domain == "out-of-domain":
            continue
        bound = next(b for k, b in bounds.items() if k in c.label)
        assert c.rmse < bound, f"{c.label}: RMSE {c.rmse:.1%}"

    # The other shapes must NEVER be scored against this set, which cannot
    # represent rib shape. Since #434 they score against Han, Zhang and Lee
    # (1991)'s own row for their shape -- and only that.
    other_shapes = [
        c
        for name, c in cells.items()
        if "4.51" in name and any(t in name for t in ("crs", "vee", "lam")) and c.n > 0
    ]
    assert other_shapes
    for c in other_shapes:
        if c.unsupported:
            continue  # the Gbar_45crs duplicate (#403 item 1)
        cls = c.label.split("_")[-1].split()[0]
        assert c.scored_by == f"han_zhang_lee_1991_{cls}", (c.label, c.scored_by)


def test_han_park_r_alpha_normalisation_is_not_raw_R(records) -> None:
    """Figure 4.47's R panel plots R divided by its own normalisers.

    Falsifiable directly: if the runner scored raw R instead of the
    normalised quantity, the bias would be enormous (R_pe and the W/H^m
    term are both > 1 across the validity box), not a few percent.
    """
    from validation.cooling.scorecard import build

    cells = {c.label: c for c in build(records)}
    r_cell = cells["han2012/fig4.47_R_vs_alpha"]
    assert abs(r_cell.bias) < 0.15, (
        f"bias {r_cell.bias:.1%} -- looks like raw R is being compared "
        "against the normalised digitised values"
    )


def test_fig446_span_check_is_one_sided_and_its_exceedances_are_known(dataset) -> None:
    """The fig4.46 cross_checks record that their span argument is ONE-SIDED:
    only isolated marks can be picked, so the picked span is a subset of the
    true span and fitting inside the predicted band is close to vacuous.
    Exceeding it is the informative direction.

    This pins both halves of that claim so neither can rot:

      1. the check really is weak where it passes -- mean coverage of the
         predicted range is well under 100%;
      2. the series that DO exceed the band are exactly the four recorded.

    A re-pick that moves a series across the boundary changes this list and
    should make someone look, rather than silently joining a set of series
    whose prose says they agree.
    """
    import csv as _csv
    import math

    # Figure 4.46's own legend: "This study" 10,000-60,000 (e/D 0.047, 0.078);
    # "Han (1984)" 8,000-80,000 (e/D 0.063, 0.042, 0.021).
    RE_RANGE = {0.047: (10000, 60000), 0.078: (10000, 60000)}
    DEFAULT_RE = (8000, 80000)
    MARGIN = 0.15  # f comes from the printed R = 3.2, worth about +/-3% on e+
    KNOWN_EXCEEDANCES = {
        "fig4.46_G_eD0.063_pe10_wh1",
        "fig4.46_R_eD0.042_pe10_wh1",
        "fig4.46_R_eD0.063_pe10_wh1",
        "fig4.46_R_eD0.063_pe20_wh1",
    }

    def e_plus(re_d: float, e_over_d: float, w_over_h: float) -> float:
        geom = 2 * e_over_d * (2 * w_over_h / (w_over_h + 1.0))
        f = 2.0 / (3.2 - 2.5 * math.log(geom) - 2.5) ** 2
        return e_over_d * re_d * math.sqrt(f / 2.0)

    coverage: list[float] = []
    exceed: set[str] = set()
    for series in dataset:
        stem = series.path.stem
        if not stem.startswith("fig4.46_") or series.kind != "measured":
            continue
        geom = series.geometry or {}
        e_over_d, w_over_h = geom.get("e_D"), geom.get("W_H")
        if e_over_d is None or w_over_h is None:
            continue
        with open(series.path) as fh:
            xs = [float(r[0]) for r in list(_csv.reader(fh))[1:] if r]
        if not xs:
            continue
        lo_re, hi_re = RE_RANGE.get(e_over_d, DEFAULT_RE)
        lo, hi = e_plus(lo_re, e_over_d, w_over_h), e_plus(hi_re, e_over_d, w_over_h)
        coverage.append((max(xs) / min(xs)) / (hi / lo))
        if any(x < lo * (1 - MARGIN) or x > hi * (1 + MARGIN) for x in xs):
            exceed.add(stem)

    assert coverage, "no fig4.46 measured series found"
    mean_cov = sum(coverage) / len(coverage)
    assert mean_cov < 0.95, (
        f"picked marks now cover {mean_cov:.0%} of the predicted range; the "
        "cross_checks describe this check as weak because coverage is partial"
    )
    assert exceed == KNOWN_EXCEEDANCES, (
        f"the set of series exceeding the predicted band changed: "
        f"added {sorted(exceed - KNOWN_EXCEEDANCES)}, "
        f"gone {sorted(KNOWN_EXCEEDANCES - exceed)}"
    )


def test_rib_angle_reaches_the_correlation_on_the_eplus_path() -> None:
    """A series' rib angle must reach `evaluate_rib`, not be replaced by the
    probe geometry's default of 90 degrees.

    It was. `run_series` copied `e_D`, `p_e` and `W_H` out of
    `series.geometry` but never `alpha_deg`, which is a top-level field
    rather than part of that mapping. Every angled series on the e+ path was
    therefore scored as though its ribs were transverse -- against
    `han_park_1988_angled`, the one set whose entire subject is rib angle.

    Checked against the library directly rather than against a stored
    number, so the test cannot drift with the data.
    """
    import combaero as cb
    from validation.cooling.runner import run_series

    dataset = load_dataset()
    series = next(s for s in dataset if s.path.stem == "fig4.51_G_60par")
    assert series.alpha_deg == 60.0

    record = next(r for r in run_series(series) if r.predicted is not None)
    geometry = cb.RibGeometry(e_D=0.0625, p_e=10.0, W_H=1.0, alpha_deg=series.alpha_deg)
    at_declared = cb.evaluate_rib(cb.han_park_1988_angled(), geometry, record.re_used)
    geometry.alpha_deg = 90.0
    at_ninety = cb.evaluate_rib(cb.han_park_1988_angled(), geometry, record.re_used)
    # The two must be far enough apart for this to be a real check.
    assert abs(at_declared.G / at_ninety.G - 1.0) > 0.05

    assert abs(record.predicted / at_declared.G - 1.0) < 0.05, (
        f"prediction {record.predicted:.3f} is not the declared "
        f"alpha={series.alpha_deg} value ({at_declared.G:.3f}); it matches "
        f"alpha=90 ({at_ninety.G:.3f}), so the rib angle is being dropped"
    )


def _cr4396(records, figure: str) -> list[Record]:
    return [
        r for r in records if r.series.source.name == "johnson1993" and r.series.figure == figure
    ]


def test_cr4396_scores_the_han_sets_cross_source_on_a_common_basis(dataset, records) -> None:
    """#452: an independent lab (Pratt & Whitney), so ACCURACY, never fidelity.

    Both Han sets over-predict CR-4396's cast-style, staggered trips by a
    flat level: +12% on the normal trips, +18% on the 45 deg skewed ones,
    with no Re trend. e/d 0.1 is above both sets' 0.047-0.078, so every
    point is out of domain. Bounds are the values recorded when #452 was
    scored, wide enough for a re-digitisation, far too narrow to absorb
    the 14.6% the basis conversion moves them.
    """
    from validation.cooling.conventions import CONVERSIONS
    from validation.cooling.scorecard import basis_of

    to_db = CONVERSIONS[("Nu_ratio_ribbed_KP", "Nu_ratio_ribbed_DB")].factor()
    assert to_db == pytest.approx(0.0176 / (0.023 * 0.72**0.4))
    for figure, set_name, lo, hi in (
        ("6.1", "han_1988_orthogonal", 0.09, 0.15),
        ("6.2", "han_park_1988_angled", 0.15, 0.22),
    ):
        rs = [r for r in _cr4396(records, figure) if r.scored_by == set_name]
        assert len(rs) == 24
        assert {r.series.scores for r in rs} == {set_name}
        assert basis_of(rs[0].series, set_name) == "accuracy"
        assert all(r.predicted is not None and r.extrapolated for r in rs)
        for r in rs:
            digitised = next(p.y for p in load_points(r.series) if p.x == r.x)
            assert r.measured == pytest.approx(digitised * to_db)
        bias = sum(r.rel_error for r in rs) / len(rs)
        assert lo < bias < hi, f"{figure}: recorded +12.1% / +18.4% in #452, now {bias:+.1%}"


def test_cr4396_on_the_raw_kays_perkins_basis_would_read_as_agreement(dataset) -> None:
    """The conversion falsified, not asserted.

    Declared as though its ordinate were already Dittus-Boelter, a CR-4396
    series must be REFUSED, not scored -- and the number it would have
    produced is the trap: Han within a few percent, where the common basis
    says 12-18% high.
    """
    from validation.cooling import conventions

    series = _series(dataset, "fig6.1_leading_segB")
    honest = run_series(series)
    try:
        conventions.SERIES_OVERRIDE[series.label] = "Nu_ratio_ribbed_DB"
        refused = run_series(series)
    finally:
        conventions.SERIES_OVERRIDE.pop(series.label, None)
    assert all(r.predicted is None for r in refused)
    assert "registered conversion" in (refused[0].reason or "")

    to_db = conventions.CONVERSIONS[("Nu_ratio_ribbed_KP", "Nu_ratio_ribbed_DB")].factor()
    raw = [r.predicted / (r.measured / to_db) - 1.0 for r in honest]
    fair = [r.rel_error for r in honest]
    assert abs(sum(raw) / len(raw)) < 0.05
    assert sum(fair) / len(fair) > 0.10


# Every set the harness can score CR-4396 with, and what each reads. Recorded
# when multi-set scoring landed; a change here is a change in what the report
# shows a designer choosing between them.
CR4396_SETS = {
    "6.1": {
        "han_1988_orthogonal": (24, 0, 0.121),
        "han_zhang_lee_1991_90": (24, 0, 0.052),
    },
    "6.2": {
        "han_1989_narrow_channel": (24, 0, 0.071),
        "han_park_1988_angled": (24, 0, 0.184),
        "han_zhang_lee_1991_45par": (24, 0, -0.016),
        "rallabandi_2009_high_re": (24, 12, 0.050),
    },
}


def test_cr4396_is_scored_by_every_eligible_set_and_no_fewer(dataset, records) -> None:
    """Several sets on one series is honest only if the list was not chosen.

    `eligible_sets` derives it by rule -- shape and angle bind, the quantity
    resolves, a prediction comes back -- and the declared `scores` +
    `also_scores` must equal it. Leaving off han_park_1988_angled, the set
    that reads worst here (+18%), fails this test; so would adding a set the
    rule rejects.
    """
    from validation.cooling.runner import eligible_sets, sets_for

    johnson = [s for s in dataset if s.source.name == "johnson1993"]
    assert len(johnson) == 12
    for s in johnson:
        assert set(sets_for(s)) == set(eligible_sets(s)), s.label
        assert set(sets_for(s)) == set(CR4396_SETS[s.figure]), s.label

    # Falsified: the same series with the worst set dropped is caught.
    s = next(s for s in johnson if s.figure == "6.2")
    worst_dropped = dataclasses.replace(
        s,
        scores="rallabandi_2009_high_re",
        also_scores=("han_1989_narrow_channel", "han_zhang_lee_1991_45par"),
    )
    assert set(sets_for(worst_dropped)) != set(eligible_sets(worst_dropped))

    for figure, sets in CR4396_SETS.items():
        for set_name, (n, n_in, bias) in sets.items():
            rs = [r for r in _cr4396(records, figure) if r.scored_by == set_name]
            assert len(rs) == n, set_name
            assert sum(not r.extrapolated for r in rs) == n_in, set_name
            got = sum(r.rel_error for r in rs) / len(rs)
            assert got == pytest.approx(bias, abs=0.005), f"{set_name}: {got:+.3f}"


def test_several_sets_are_reported_side_by_side_never_ranked(dataset, records) -> None:
    """One row per set; none picked, none pooled with another."""
    from validation.cooling.runner import SETS
    from validation.cooling.scorecard import (
        SET_LAB,
        build,
        multi_set_rows,
        pool,
        render_multi_set,
    )

    johnson = [s for s in dataset if s.source.name == "johnson1993"]
    cells = build([r for r in records if r.series.source.name == "johnson1993"], dataset)
    rows = multi_set_rows(cells)
    assert len(rows) == 12
    for series, sets in rows.items():
        names = [name for name, _, _ in sets]
        assert names == sorted(names), "ordered by name, never by score"
        figure = next(s.figure for s in johnson if s.label == series)
        assert set(names) == set(CR4396_SETS[figure])

    text = "\n".join(render_multi_set(cells))
    assert "han_park_1988_angled" in text  # the worst-reading set is shown
    assert "best" not in text.lower()
    assert "NOT independent corroboration" in text

    # A pool takes each series' own set and never mixes two: one series
    # pools to its 4 own-set points, not 16; Figs. 6.1 and 6.2 together
    # (two own sets) refuse.
    assert pool(records, "johnson1993/fig6.2_leading_segB")[0] == 4
    with pytest.raises(ValueError, match="spans sets"):
        pool(records, "johnson1993/fig6")

    # Every rib set states its lab, so the independence note cannot go
    # silent for a set nobody labelled.
    assert set(SETS) <= set(SET_LAB)


def test_only_the_rib_runner_scores_one_series_with_several_sets(dataset) -> None:
    """A specialist runner would silently use `scores` alone and report it as
    though it were the only eligible set."""
    from validation.cooling.scorecard import run_dataset

    film = next(s for s in dataset if s.scores == "baldauf_2002_sellers")
    with pytest.raises(ValueError, match="only supported by the rib runner"):
        run_dataset([dataclasses.replace(film, also_scores=("colban_2011_shaped",))])


def test_every_scored_rib_series_lists_every_eligible_set(dataset) -> None:
    """The CR-4396 rule, everywhere the rib runner scores. A series whose
    `scores` is set must be scored by every set `eligible_sets` derives --
    so no series shows only the set that flatters it. `scores: null` stays a
    recorded decision (a printed correlation curve, a resolved mislabel),
    not an omission this test could judge."""
    from validation.cooling.runner import eligible_sets, sets_for
    from validation.cooling.scorecard import rib_owned

    checked = 0
    for s in dataset:
        if s.scores is None or not rib_owned(s):
            continue
        checked += 1
        assert set(sets_for(s)) == set(eligible_sets(s)), (
            f"{s.label}: declared {sorted(sets_for(s))}, eligible {sorted(eligible_sets(s))}"
        )
    assert checked >= 75


# Measured series a set COULD score but that are deliberately not scored,
# each with its recorded reason. Anything else with an eligible set must be
# scored: an unscored series is otherwise a choice made after the fact.
RECORDED_UNSCORED = {
    "han2012/fig4.51_Gbar_45par": "carries the G marks a second time (#403 item 1)",
    "han2012/fig4.51_Gbar_45crs": "carries the G marks a second time (#403 item 1)",
    "han2012/fig4.48_R_errorbar_WH1_ref3": "a published correlation's own scatter band",
    "han2012/fig4.48_G_errorbar_WH1_ref3": "a published correlation's own scatter band",
}


def test_no_measured_rib_series_is_left_unscored_without_a_reason(dataset) -> None:
    """The other half of the CR-4396 rule. Taslim's two-wall AR 0.5 and 3.5
    friction sat unscored under a note #443 had already made stale, while
    the AR 1.0 series beside them were scored -- the omitted ones read
    25-45% worse. A printed correlation curve is not a measurement and is
    exempt by kind."""
    from validation.cooling.runner import eligible_sets
    from validation.cooling.scorecard import rib_owned

    found = set()
    for s in dataset:
        if s.scores is not None or s.kind != "measured" or not rib_owned(s):
            continue
        if not eligible_sets(s):
            continue
        found.add(s.label)
        assert s.label in RECORDED_UNSCORED, (
            f"{s.label}: unscored, but {list(eligible_sets(s))} can score it"
        )
    assert found == set(RECORDED_UNSCORED), "a recorded exemption no longer applies"


def test_eligibility_needs_a_known_angle_and_a_declared_quantity(dataset, monkeypatch) -> None:
    """Two ways a set can LOOK scorable and not be.

    A single-angle set does not qualify for a cloud that does not state its
    angle: figure 4.47's G pools 30-90 deg, and han_1988_orthogonal (90 deg
    only) would score it as all-transverse. And an undeclared measurement
    convention lets a series' own set score it, but qualifies no other.
    """
    from validation.cooling.runner import eligible_sets

    pooled = next(s for s in dataset if s.label == "han2012/fig4.47_G_vs_eplus")
    assert pooled.alpha_deg is None
    got = eligible_sets(pooled)
    assert "han_park_1988_angled" in got
    for single_angle in ("han_1988_orthogonal", "rallabandi_2009_high_re", "han_zhang_lee_1991_90"):
        assert single_angle not in got
    # Falsified: the same series with its angle stated admits the 90 deg set.
    assert "han_1988_orthogonal" in eligible_sets(dataclasses.replace(pooled, alpha_deg=90.0))

    from validation.cooling import conventions

    monkeypatch.setattr(conventions, "series_convention", lambda series: None)
    assert conventions.resolve(pooled, "han_1989_narrow_channel").status == (
        conventions.STATUS_UNDECLARED
    )
    assert eligible_sets(pooled) == ("han_park_1988_angled",)


def test_a_set_scoring_another_sets_figure_is_accuracy(dataset) -> None:
    """Figures 4.47 and 4.48 share one book citation but are two sets' own
    data. Once every eligible set scores both, the shared key would call
    han_park_1988_angled on the narrow-channel figure "fidelity" -- a miss
    claimed as our transcription bug on data that set was never fitted to.
    Reprints of a set's own paper stay fidelity: figure 4.51 is Han, Zhang
    and Lee (1991) Fig. 12."""
    from validation.cooling.scorecard import basis_of

    def first(figure: str):
        return next(
            s for s in dataset if s.source.name == "han2012" and s.figure == figure and s.scores
        )

    assert basis_of(first("4.47"), "han_park_1988_angled") == "fidelity"
    assert basis_of(first("4.48"), "han_1989_narrow_channel") == "fidelity"
    assert basis_of(first("4.48"), "han_park_1988_angled") == "accuracy"
    assert basis_of(first("4.47"), "han_1989_narrow_channel") == "accuracy"
    g45 = next(s for s in dataset if s.label == "han2012/fig4.51_G_45par")
    assert basis_of(g45, "han_zhang_lee_1991_45par") == "fidelity"


# ---------------------------------------------------------------------------
# Naik, Retzko, Gritsch and Sedlov (2014), GT2014-25841: Alstom wedge
# trailing-edge passage, W/H 6.64. Record:
# validation/cooling/extractions/naik_2014_te_wedge_ribs.md.
# ---------------------------------------------------------------------------


def _naik(dataset, stem: str):
    return next(s for s in dataset if s.source.name == "naik2014" and s.path.stem == stem)


def _naik_points(dataset, fig: str, series: str) -> dict[float, float]:
    return {p.x: p.y for p in load_points(_naik(dataset, f"fig{fig}_{series}"))}


# Every set the harness can score the angled-tapered Nu/Nu0 with: (points,
# points inside the set's own box, bias). Recorded when the series was
# committed (2026-11-09); every point is outside every box (W/H 6.64, e/Dh
# 0.089, a wedge, tapered ribs, Re to 145k).
NAIK_NU_SETS = {
    "han_1989_narrow_channel": (3, 0, 0.451),
    "han_park_1988_angled": (3, 0, -0.148),
    "han_zhang_lee_1991_45par": (3, 0, 0.342),
    "rallabandi_2009_high_re": (3, 0, 0.246),
}


def test_naik2014_nu_is_scored_by_every_eligible_set_as_accuracy(dataset, records) -> None:
    """Alstom's rig, not Texas A&M: cross-source ACCURACY on every row, and
    out of domain for every set -- reported, not dropped. The list of sets is
    derived by rule (eligible_sets), so han_park_1988_angled, the only one
    that reads LOW, cannot be the only one shown, nor the three that read
    25-45% high be left off.

    The ordinate is already on the sets' Dittus-Boelter basis (Eq. 3), so it
    resolves DIRECT: the measured value scored is the digitised one,
    unconverted.
    """
    from validation.cooling import conventions
    from validation.cooling.runner import eligible_sets, sets_for
    from validation.cooling.scorecard import basis_of

    s = _naik(dataset, "fig16_angled_tapered")
    assert set(sets_for(s)) == set(eligible_sets(s)) == set(NAIK_NU_SETS)
    digitised = _naik_points(dataset, "16", "angled_tapered")
    assert sorted(digitised) == [22000.0, 60000.0, 145000.0]
    for set_name, (n, n_in, bias) in NAIK_NU_SETS.items():
        assert conventions.resolve(s, set_name).status == conventions.STATUS_DIRECT
        assert basis_of(s, set_name) == "accuracy"
        rs = [r for r in records if r.series is s and r.scored_by == set_name]
        assert len(rs) == n, set_name
        assert sum(not r.extrapolated for r in rs) == n_in, set_name
        for r in rs:
            assert r.measured == digitised[r.x]
        got = sum(r.rel_error for r in rs) / len(rs)
        assert got == pytest.approx(bias, abs=0.005), f"{set_name}: {got:+.3f}"


def test_naik2014_nu_axis_and_convention_must_agree(dataset) -> None:
    """A Dittus-Boelter axis that resolves through a conversion contradicts
    itself, as a Kays-Perkins axis that resolves direct does (the CR-4396
    trap). Either is refused -- so no conversion is skipped or applied twice."""
    from validation.cooling import conventions

    s = _naik(dataset, "fig16_angled_tapered")
    assert all(r.predicted is not None for r in run_series(s))
    try:
        conventions.SERIES_OVERRIDE[s.label] = "Nu_ratio_ribbed_KP"
        refused = run_series(s)
    finally:
        conventions.SERIES_OVERRIDE.pop(s.label, None)
    assert all(r.predicted is None for r in refused)
    assert "Dittus-Boelter basis, unconverted" in (refused[0].reason or "")


def test_naik2014_friction_is_scored_by_rule_and_the_conversion_is_definition(
    dataset, records
) -> None:
    """f/f0 is scored by han_park_1988_angled, the one set eligible by rule:
    the only 45 deg set that declares the passage-average friction
    (f_channel_two_ribbed). han_1988_orthogonal declares it too but binds
    transverse 90 deg ribs, so it refuses on angle, not on quantity.

    The measured side is exactly (f/f0) * Eq. 5 -- the registered conversion,
    definition, not a fit. Every row is out of the set's box (W/H 6.64 > 4,
    e/D 0.089 > 0.078) and cross-source ACCURACY. Recorded 2026-11-10: the
    set reads +42/+64/+44% HIGH, none inside the paper's 10% friction band."""
    from validation.cooling import conventions
    from validation.cooling.runner import eligible_sets, sets_for
    from validation.cooling.scorecard import basis_of

    s = _naik(dataset, "fig17_angled_tapered")
    assert sets_for(s) == eligible_sets(s) == ("han_park_1988_angled",)
    resolution = conventions.resolve(s, "han_park_1988_angled")
    assert resolution.status == conventions.STATUS_CONVERTED
    assert resolution.detail == "f_ratio_passage_petukhov -> f_channel_two_ribbed"
    assert basis_of(s, "han_park_1988_angled") == "accuracy"
    orthogonal = run_series(s, set_name="han_1988_orthogonal")
    assert all(r.predicted is None for r in orthogonal)
    assert "binds Transverse ribs only" in (orthogonal[0].reason or "")

    rs = [r for r in records if r.series is s]
    digitised = _naik_points(dataset, "17", "angled_tapered")
    assert len(rs) == 3
    for r in rs:
        assert r.scored_by == "han_park_1988_angled"
        assert r.predicted is not None and r.extrapolated
        assert r.measured == pytest.approx(
            digitised[r.x] * 2.0 / (2.236 * math.log(r.x) - 4.639) ** 2
        )
        assert not r.within_uncertainty
    bias = sum(r.rel_error for r in rs) / len(rs)
    assert bias == pytest.approx(0.499, abs=0.005)


def test_naik_petukhov_is_fanning() -> None:
    """Eq. 5 is Petukhov's Darcy (0.790 ln Re - 1.64)^-2 over 4 -- Fanning,
    like Eq. 4's f -- to 0.2%. A Darcy reading would be 4x off."""
    import combaero as cb
    from validation.cooling.conventions import naik_petukhov_fanning

    for re in (22000.0, 60000.0, 145000.0):
        assert naik_petukhov_fanning(re) == pytest.approx(cb.friction_petukhov(re) / 4.0, rel=0.003)


NAIK_OWN = ("angled_tapered", "broken", "chevron")
NAIK_REPRINTS = (
    "reprint_rhee2003_45degV",
    "reprint_rhee2003_60degV",
    "reprint_wright2004_parallel",
    "reprint_wright2004_discreteW",
)


def _naik_pi_triples(dataset, series: str) -> list[tuple[float, float, float]]:
    """(Nu/Nu0, f/f0, PI) at every Re all three figures carry."""
    nu = _naik_points(dataset, "16", series)
    f = _naik_points(dataset, "17", series)
    pi = _naik_points(dataset, "18", series)
    common = sorted(set(nu) & set(f) & set(pi))
    return [(nu[x], f[x], pi[x]) for x in common]


def test_naik_fig18_own_points_are_3_nu_over_f(dataset) -> None:
    """The paper's own Fig. 18 does not follow its own Eq. 6.

    PI = (Nu/Nu0)/(f/f0)^(1/3) is printed; the 8 points of the paper's own
    configurations are instead 3 (Nu/Nu0)/(f/f0) to 0.3% -- A/(B^1/3)
    evaluated with ^ before /. They sit 24-172% below Eq. 6, so Fig. 18's
    own curves (and the abstract's PI ranking) are not Eq. 6. Guards the
    digitisation of Figs. 16-18 too: three figures read independently agree
    to 0.3% under one identity.
    """
    n = 0
    for series in NAIK_OWN:
        for nu, f, pi in _naik_pi_triples(dataset, series):
            n += 1
            assert pi == pytest.approx(3.0 * nu / f, rel=0.004), series
            assert nu / f ** (1.0 / 3.0) / pi - 1.0 > 0.20, series
    assert n == 8  # chevron at 22k has f/f0 only


def test_naik_fig18_reprints_follow_eq6(dataset) -> None:
    """The 12 reprinted points DO follow Eq. 6, to 0.13% -- so the error is
    in the paper's own column, not in our reading of the figure."""
    n = 0
    for series in NAIK_REPRINTS:
        for nu, f, pi in _naik_pi_triples(dataset, series):
            n += 1
            assert pi == pytest.approx(nu / f ** (1.0 / 3.0), rel=0.002), series
            assert abs(3.0 * nu / f / pi - 1.0) > 0.10, series
    assert n == 12


def test_naik2014_unscored_series_have_no_eligible_set(dataset) -> None:
    """Broken and chevron ribs (two rows of discrete ribs) and the reprints
    (no geometry in Naik) have no set by rule, and Fig. 18 is not a scored
    quantity: nothing here is unscored that a set could score."""
    from validation.cooling.runner import eligible_sets

    naik = [s for s in dataset if s.source.name == "naik2014"]
    assert len(naik) == 21
    scored = [s.path.stem for s in naik if s.scores]
    assert scored == ["fig16_angled_tapered", "fig17_angled_tapered"]
    for s in naik:
        if s.scores is None:
            assert eligible_sets(s) == (), s.label
        if "reprint" in s.path.stem:
            assert s.after and s.geometry is None and s.scores is None, s.label
