"""Aggregate every cooling runner into one scorecard.

Three runners exist because three physics families score differently -- ribs
need a Reynolds-number bisection, Nu/Nu1 cancels it, and the orifice chain
works in the source's own pressure coordinates. Each owns its series and
declines the rest by returning nothing, so ``run_dataset`` can give every
series to exactly one of them. Without that, a series declared in metadata but
scored by a specialist runner appeared here as "not scored by any set", and
there was no single place to see whether cooling as a whole was healthy.

Metrics per series:
    N             digitised points
    N_scored      points the set had a binding for
    MAE           mean absolute relative error
    RMSE          root mean square relative error
    bias          signed mean relative error (model - source)
    within        fraction inside the source's stated uncertainty
    extrap        points outside the set's advisory validity

    sampling      how completely the marks could be picked
    bnd           runs recovered as interval observations
    held          fraction of those the prediction lands inside

A series the set has no binding for reports '-' rather than 0. A zero
there would read as a perfect score for a model that answered nothing,
which is the shape of mistake this harness exists to catch.

**Rows are segregated by sampling completeness and never pooled across
it.** A figure that overplots several symbol classes yields only its
spatially isolated marks to a digitisation, and those are the ones
furthest from the cluster centre -- so a `partial` row's MAE/RMSE/within
are an UPPER BOUND on the model's error, not an estimate of it. The
difference is not academic: pooled, `han_1988_orthogonal` reported 73.0%
within against Han's own printed claim of 95%, which reads as a model
deficiency. Segregated, its completely-sampled series reach 91.1% and
its partially-sampled ones 65.3%. The gap was the pick, not the
correlation. See #393 and `recovery.py`.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass

from validation.cooling import conventions, fidelity
from validation.cooling.runner import Record
from validation.cooling.schema import Point, load_points


# The marker for a row scored outside its correlation set's declared
# validity box. Defined once: it is written by `build` onto a series label
# and by `rollup` onto a set label, and the two drifted apart the first
# time they were written separately.
OUT_OF_DOMAIN_MARK = "{out}"


@dataclass
class Cell:
    label: str
    kind: str
    scored_by: str | None
    n: int = 0
    n_scored: int = 0
    mae: float = float("nan")
    rmse: float = float("nan")
    bias: float = float("nan")
    within: float = float("nan")
    n_extrapolated: int = 0
    reason: str | None = None  # why nothing was scored
    # Sampling completeness -- see `sampling_of`. "partial" means the
    # digitisation could not pick every mark, so MAE/RMSE/within are an
    # UPPER BOUND on the model's error rather than an estimate of it: the
    # marks that could be picked are the ones furthest from the cluster
    # centre (#393).
    sampling: str = "unknown"
    # "fidelity" (the correlation's own paper) or "accuracy" (cross-source).
    # Never pooled: they answer different questions and fail differently.
    basis: str = "unknown"
    # How the series' measurement convention met the set's: "direct",
    # "converted", "undeclared", "incompatible" or "unregistered". The last
    # two are REFUSALS -- scoring two different quantities is a category
    # error, not a disagreement, and no range check can see it (#389).
    convention: str = "undeclared"
    convention_detail: str = ""
    # "in-domain" or "out-of-domain" against the correlation set's OWN
    # declared validity box. Never pooled either: a set must not be judged
    # on conditions its authors never claimed, and today several are --
    # `baldauf_2002_sellers` scores 138 points of which 138 are outside its
    # envelope, and `han_1988_orthogonal`'s largest row is 52 of 98. One
    # number over both says nothing about either (#389).
    domain: str = "in-domain"
    # The series this row belongs to, and whether several sets score it --
    # which is how the several-sets section finds its rows, without
    # parsing labels.
    series_label: str = ""
    multi_set: bool = False
    n_bounded: int = 0  # runs recovered as interval observations
    held: float = float("nan")  # fraction of those the prediction lands inside

    @property
    def unsupported(self) -> bool:
        return self.n_scored == 0


def _specialists() -> tuple:
    from validation.cooling import (
        effusion_internal_runner,
        effusion_network_runner,
        effusion_overall_runner,
        film_runner,
        jet_array_runner,
        orifice_runner,
        ratio_runner,
        shaped_film_runner,
    )

    return (
        (jet_array_runner.owns, jet_array_runner.run_series),
        (ratio_runner.owns, ratio_runner.run_series),
        (orifice_runner.owns, orifice_runner.run_series),
        # Before film_runner, which claims every x_over_D / eta_adiabatic
        # series as cylindrical-hole Baldauf + Sellers.
        (shaped_film_runner.owns, shaped_film_runner.run_series),
        (film_runner.owns, film_runner.run_series),
        (effusion_internal_runner.owns, effusion_internal_runner.run_series),
        (effusion_overall_runner.owns, effusion_overall_runner.run_series),
        # Behrendt 2008 through the element's network (both gas baselines).
        (effusion_network_runner.owns, effusion_network_runner.run_series),
    )


def rib_owned(series) -> bool:
    """True when no specialist claims the series, so the rib runner scores it."""
    return not any(owns(series) for owns, _ in _specialists())


def run_dataset(dataset) -> list[Record]:
    """Give every series to exactly one runner.

    Ownership is an explicit predicate, not "the runner returned something".
    Every runner returns reason-carrying records for series it cannot score --
    that is the point of them, so a silent zero never masquerades as a good
    score -- which makes truthiness useless as a dispatch test and was the
    first version's bug: the jet-array runner claimed the orifice series and
    reported them unscored.
    """
    from validation.cooling.runner import run_series as rib_run
    from validation.cooling.runner import sets_for

    specialists = _specialists()
    out: list[Record] = []
    for series in dataset:
        for owns, run in specialists:
            if owns(series):
                if series.also_scores:
                    # Only the rib runner scores one series with several
                    # sets. Dropping the others silently would report one
                    # set as though it were the only one eligible.
                    raise ValueError(
                        f"{series.label}: also_scores is only supported by "
                        "the rib runner"
                    )
                out.extend(run(series))
                break
        else:
            # The rib runner is the fallback and reports why, rather than
            # dropping a series no specialist wanted.
            for name in sets_for(series) or (None,):
                out.extend(rib_run(series, set_name=name))
    return out


# Which paper each correlation set IS. Fidelity can only be judged on data
# from that paper; everything else is accuracy, and must be labelled so --
# see docs/VALIDATION_POLICY.md.
#
# This cannot be read off the C++ set: `RibCorrelationSet.provenance` is an
# Extracted/Fitted/User enum, not a citation. It is matched against a
# series' `after` field (the original author of a reprinted figure), its
# source name, and its source CITATION, so a textbook reprint of the set's
# own paper still counts as the author's own data.
#
# The CITATION is in the haystack because a folder name cannot always carry
# the identifier a set is keyed on. `andrews1986` holds BOTH 86-GT-225's
# Fig. 8 and Fig. 10; the set is keyed "86-GT-225" to keep it from matching
# the 1988 study, and that string appears only in the citation. Without it
# the paper's own data reported "accuracy", understating a fidelity result
# -- the mirror of the error the key was introduced to prevent.
SET_ORIGIN: dict[str, str] = {
    "han_1988_orthogonal": "ASME J. Heat Transfer 110, 321",
    "han_park_1988_angled": "IJHMT 31(1), 183",
    # Fig. 4.48 (where all of this set's scored series live) is credited
    # "Han, J.C. et al. (1988). IJHMT 31(1), 183" -- the same journal
    # citation as han_park_1988_angled's Fig. 4.47, even though Eq. 4.19's G
    # is attributed in the running text to a companion Han et al. (1989).
    # Same figure, same research effort: Fig. 4.48 is this set's own data.
    "han_1989_narrow_channel": "IJHMT 31(1), 183",
    "rallabandi_2009_high_re": "rallabandi",
    "florschuetz_1981_inline": "florschuetz1981",
    "mcgreehan_schotsch_1988_cd": "mcgreehan_schotsch1988",
    "mcgreehan_schotsch_1988_crossflow_cd": "mcgreehan_schotsch1988",
    # Jobson (1955) is theory with no data of its own: the Perry points his
    # Fig. 5a plots are a DIFFERENT study, so they are cross-source accuracy.
    # SET_OWN_FIGURES below pins that (no figure of his is his own data).
    "jobson_1955_mass_flow": "jobson1955",
    # Baldauf's own data is not in this repo, so nothing scored by this set
    # can ever read "fidelity" -- which is correct. Naming the origin makes
    # the verdict "accuracy" rather than "unknown", so the report states the
    # cross-source basis instead of declining to say.
    "baldauf_2002_sellers": "baldauf",
    # 86-GT-225, NOT "andrews": the acceptance data is 88-GT-290, a
    # different study from the same lab, which the validation policy counts
    # as CROSS-SOURCE. Matching on the bare author name would have matched
    # the 1988 source and reported fidelity, overstating the claim.
    "andrews_1986_effusion_internal": "86-GT-225",
    # The scored points are Gritsch et al. (2005), reprinted in Colban's Figs.
    # 6-7 and among the data Eq. (19) was fitted on -- the correlation's own
    # data, so fidelity (on its fitting set, which the record says).
    "colban_2011_shaped": "colban2011",
}

# Taslim & Spring (1987) ratio sets (#444): their own paper's data.
for _ar, _ed in ((0.5, 0.125), (0.5, 0.250), (1.0, 0.083), (1.0, 0.167),
                 (3.5, 0.053), (3.5, 0.107), (3.5, 0.161)):
    SET_ORIGIN[f"taslim_spring_1987_ar{_ar:.1f}_eD{_ed:.3f}"] = "taslim"

# Han, Zhang and Lee (1991) Table 2 -- every configuration's origin is the
# same paper, whose Fig. 12 figure 4.51 reprints.
for _tag in ("90", "60par", "60crs", "60vee", "60lam", "45par", "45crs", "45vee", "45lam"):
    SET_ORIGIN[f"han_zhang_lee_1991_{_tag}"] = "ASME J. Heat Transfer 113, 590"

# Two sets keyed to ONE citation: the book credits figures 4.47 and 4.48 both
# to IJHMT 31(1), 183, but han_park_1988_angled is fitted to figure 4.47
# (square channel, angled ribs) and han_1989_narrow_channel to figure 4.48
# (narrow channels). While each scored only its own figure the shared key
# was harmless. Scored by every eligible set, each would read the OTHER's
# figure as "fidelity" -- claiming a miss on data the correlation was never
# fitted to is our transcription bug. So within that source, fidelity is
# the set's own figure only.
SET_OWN_FIGURES: dict[str, tuple[str, frozenset[str]]] = {
    "han_park_1988_angled": ("han2012", frozenset({"4.47"})),
    "han_1989_narrow_channel": ("han2012", frozenset({"4.48"})),
    # Every point under jobson1955 is someone else's measurement (Perry 1949),
    # so none is the theory's own data: accuracy, never fidelity (#854).
    "jobson_1955_mass_flow": ("jobson1955", frozenset()),
}

# Which LAB each rib set comes from. Two sets from one lab agreeing is not
# two labs agreeing: they share a rig family, a data-reduction convention
# and often the data itself. Every shipped rib set is Han's laboratory at
# Texas A&M (#403 item 2), so where several score one series the report
# says their agreement is not independent corroboration -- only the
# measured series, from another lab, is. Required for every set in the
# rib runner's registry (test_every_rib_set_names_its_lab).
SET_LAB: dict[str, str] = {
    "han_1988_orthogonal": "Texas A&M (Han)",
    "han_park_1988_angled": "Texas A&M (Han)",
    "han_1989_narrow_channel": "Texas A&M (Han)",
    "rallabandi_2009_high_re": "Texas A&M (Han)",
}
for _tag in ("90", "60par", "60crs", "60vee", "60lam", "45par", "45crs", "45vee", "45lam"):
    SET_LAB[f"han_zhang_lee_1991_{_tag}"] = "Texas A&M (Han)"

# Matching on the citation widens the haystack, so a citation that NAMES
# another paper could claim a false fidelity. Bounded by
# `test_every_scored_series_basis_is_what_its_provenance_says`, which pins
# the label of every scored series rather than trusting the substring.


def basis_of(series, set_name: str | None) -> str:
    """"fidelity" when the data is the correlation's own paper, else "accuracy".

    A set with no recorded origin returns "unknown" rather than defaulting
    to fidelity, so adding a set without declaring its paper cannot quietly
    claim the stronger of the two.
    """
    origin = SET_ORIGIN.get(set_name or "")
    if origin is None:
        return "unknown"
    own = SET_OWN_FIGURES.get(set_name or "")
    if own is not None and series.source.name == own[0] and series.figure not in own[1]:
        return "accuracy"
    haystack = (
        f"{series.after or ''} {series.source.name} {series.source.citation}"
    ).lower()
    return "fidelity" if origin.lower() in haystack else "accuracy"


def sampling_of(series, dataset) -> str:
    """How completely this series' marks could be picked.

    "complete"  a table, a drawn curve, or a panel it does not share --
                nothing could have occluded a mark.
    "partial"   MEASURED shortfall: the figure's other panel proves runs
                exist that this panel has no mark for.
    "unknown"   shares a panel with other classes, but nothing independent
                says how many runs it should have.

    Derived, never declared, so it cannot go stale against the data.
    """
    from validation.cooling.recovery import PAIR_TOL, _partner

    # "pdf-vector-exact" joins "tabulated" for the same reason: the marks are
    # separate path objects in the PDF, each with its own stroke colour, so
    # overplotting cannot occlude one. There is no pixel measurement that a
    # neighbouring curve could have hidden -- the failure mode `partial`
    # exists to flag simply cannot occur.
    if series.extraction in ("tabulated", "pdf-vector-exact") or series.kind in (
        "correlation",
        "frame",
    ):
        return "complete"
    partner = _partner(series, dataset)
    if partner is not None:
        mine = [p.x for p in load_points(series)]
        missing = sum(
            1
            for q in load_points(partner)
            if not any(abs(q.x / x - 1.0) <= PAIR_TOL for x in mine)
        )
        return "partial" if missing else "complete"
    shares = sum(
        1
        for other in dataset
        if other.kind == "measured"
        and other.source.name == series.source.name
        and other.figure == series.figure
        and other.panel == series.panel
        and other.path != series.path
    )
    return "unknown" if shares else "complete"


def score_recovered(series, dataset, set_name: str | None = None) -> tuple[int, float]:
    """Check the set against runs the digitisation could not pick.

    Returns (count, fraction the prediction lands inside the bound).

    Deliberately NOT folded into MAE/RMSE. A point observation gives a
    signed error; an interval gives containment. Averaging the two into
    one number would be a category error, and would let a loose bound
    flatter a set that a tight point disagrees with.

    Read `held` with the envelope in mind: it is the spread of OTHER
    classes near that abscissa, so it measures whether the prediction
    lands in the local cloud -- not whether it matches this class.
    """
    from validation.cooling.recovery import recover
    from validation.cooling.runner import run_series

    bounds = recover(series, dataset)
    if not bounds:
        return 0, float("nan")
    probes = [Point(x=b.x, y=(b.lo + b.hi) / 2.0) for b in bounds]
    records = run_series(series, probes, set_name)
    checked = [
        (b, r.predicted)
        for b, r in zip(bounds, records, strict=True)
        if r.predicted is not None
    ]
    if not checked:
        return 0, float("nan")
    inside = sum(1 for b, pred in checked if b.contains(pred))
    return len(checked), inside / len(checked)


def build(records: list[Record], dataset=None) -> list[Cell]:
    """One row per series.

    Pass `dataset` to fill in sampling completeness and the recovered
    interval observations. Without it those columns stay 'unknown' and
    empty -- the metrics are unchanged either way.
    """
    grouped: dict[tuple[str, str], list[Record]] = defaultdict(list)
    for r in records:
        # A runner may split one series into reported groups -- Rohde's are
        # per velocity-head-ratio band, because pooling them would mix a
        # near-exact comparison with one dominated by a conversion factor.
        group = getattr(r, "group", None)
        label = r.series.label
        scored_by = getattr(r, "scored_by", None)
        if scored_by is not None and scored_by != r.series.scores:
            # One row per set. An ADDED set rides on the label so two sets'
            # rows never merge, and it goes before any group so `rollup`
            # still finds the group after "  [". The series' own `scores`
            # set keeps the bare label: it is the pairing by provenance
            # (the set's own paper, or the one its issue scoped), not one
            # picked by score, and its rows stay as they were.
            label += f"  @{scored_by}"
        label += f"  [{group}]" if group else ""
        # And the harness splits every series by the set's own validity box.
        # `extrapolated` already flags each point; until now it was only
        # counted in a column, which made it a number rather than a
        # decision (#389).
        domain = (
            "out-of-domain" if getattr(r, "extrapolated", False) else "in-domain"
        )
        grouped[(label, domain)].append(r)

    cells: list[Cell] = []
    for label, domain in sorted(grouped):
        rs = grouped[(label, domain)]
        series = rs[0].series
        errs = [r.rel_error for r in rs if r.rel_error is not None]
        cell = Cell(
            # The marker rides on the label so a split series is legible in
            # the per-series listing without a column nothing else needs.
            label=label
            + (f"  {OUT_OF_DOMAIN_MARK}" if domain == "out-of-domain" else ""),
            series_label=series.label,
            multi_set=bool(series.also_scores),
            domain=domain,
            kind=series.kind,
            scored_by=getattr(rs[0], "scored_by", None) or series.scores,
            n=len(rs),
            n_scored=len(errs),
            n_extrapolated=sum(1 for r in rs if getattr(r, "extrapolated", False)),
            reason=next(
                (r.reason for r in rs if getattr(r, "reason", None)), None
            ),
        )
        resolution = conventions.resolve(series, cell.scored_by)
        cell.convention = resolution.status
        cell.convention_detail = resolution.detail
        if not resolution.scorable:
            # Refused, and the refusal replaces the numbers rather than
            # sitting beside them. A category error reported as a large
            # error would be read as a model limitation.
            cell.reason = f"{resolution.status}: {resolution.detail}"
            cell.n_scored = 0
            cells.append(cell)
            continue

        if dataset is not None:
            cell.basis = basis_of(series, cell.scored_by)
            cell.sampling = sampling_of(series, dataset)
            cell.n_bounded, cell.held = score_recovered(
                series, dataset, cell.scored_by
            )
        if errs:
            cell.mae = sum(abs(e) for e in errs) / len(errs)
            cell.rmse = math.sqrt(sum(e * e for e in errs) / len(errs))
            cell.bias = sum(errs) / len(errs)
            # A series with no stated band gets '-', not 0%. Zero would read
            # as "nothing agreed" when the truth is "the source states no
            # band to agree within".
            cell.within = (
                sum(1 for r in rs if r.within_uncertainty) / len(errs)
                if series.uncertainty is not None
                else float("nan")
            )
        cells.append(cell)
    return cells


def pool(records: list[Record], prefix: str) -> tuple:
    """Aggregate every series whose label starts with prefix.

    Where class attribution is provisional -- overlapping symbols, marks
    that appear in only one panel -- the pooled cloud is what the scoring
    should rest on. For a set with no geometry exponents every class must
    land on one curve anyway, so pooling costs nothing and removes a
    dependence on labels that cannot be fully trusted.
    """
    # Each series' own `scores` set only. The other sets' records for the
    # same series are a different question, answered on their own rows.
    records = [
        r for r in records
        if r.series.label.startswith(prefix)
        and r.rel_error is not None
        and getattr(r, "scored_by", r.series.scores) == r.series.scores
    ]
    sets = {r.series.scores for r in records}
    if len(sets) > 1:
        # A pool across sets is a best-of in disguise: it averages one
        # set's misses with another's hits into a number neither earns.
        raise ValueError(f"pool {prefix!r} spans sets {sorted(map(str, sets))}")
    errs = [r.rel_error for r in records]
    if not errs:
        return (0, float("nan"), float("nan"), float("nan"), float("nan"))
    within = [r for r in records if r.within_uncertainty]
    mean = sum(errs) / len(errs)
    return (
        len(errs),
        sum(abs(e) for e in errs) / len(errs),
        math.sqrt(sum(e * e for e in errs) / len(errs)),
        mean,
        len(within) / len(errs),
    )


def rollup(cells: list[Cell]) -> list[Cell]:
    """One row per correlation set -- the whole-of-cooling health view.

    Grouped by (set, reporting group) rather than by set alone. Where a runner
    split a series into groups it did so because pooling them would be
    misleading, and a rollup that ignored that would put the misleading number
    back at the top of the page.

    Unscored series contribute their point count but no error, so a set that
    answered nothing cannot look perfect.
    """
    buckets: dict[tuple[str, str, str, str], list[Cell]] = defaultdict(list)
    for c in cells:
        if c.scored_by is None:
            continue
        # The out-of-domain mark trails the group's bracket, so it comes off
        # first: stripping "]" before it left "group]" on a grouped row that
        # is out of domain (first hit by behrendt2008's flat-plate rows).
        group = (
            c.label.split("  [")[1]
            .replace(f"  {OUT_OF_DOMAIN_MARK}", "")
            .rstrip("]")
            if "  [" in c.label
            else ""
        )
        # Segregated by sampling completeness as well as by set and group.
        # A tabulated series and an overplotted one do not measure the same
        # thing: the first gives the model's error, the second an upper
        # bound on it, because only the marks furthest from the cluster
        # centre could be picked. Averaging them produces a number that is
        # neither (#393).
        buckets[(c.scored_by, group, c.sampling, c.basis, c.domain)].append(c)

    out: list[Cell] = []
    for (set_name, group, sampling, basis, domain), cs in sorted(buckets.items()):
        n = sum(c.n for c in cs)
        n_scored = sum(c.n_scored for c in cs)
        tag = (
            f"  [{basis[:3]}]"
            + ("" if sampling == "complete" else f" <{sampling[:4]}>")
            + ("" if domain == "in-domain" else f" {OUT_OF_DOMAIN_MARK}")
        )
        agg = Cell(
            label=set_name + (f"  [{group}]" if group else "") + tag,
            kind=f"{len(cs)} series",
            scored_by=set_name,
            n=n,
            n_scored=n_scored,
            n_extrapolated=sum(c.n_extrapolated for c in cs),
            sampling=sampling,
            basis=basis,
            domain=domain,
            n_bounded=sum(c.n_bounded for c in cs),
        )
        held = [(c.held, c.n_bounded) for c in cs if c.n_bounded]
        if held:
            agg.held = sum(h * k for h, k in held) / sum(k for _, k in held)
        if n_scored:
            # Weight each series by the points it actually scored, so a
            # nine-point series does not count the same as a one-point one.
            w = [(c, c.n_scored) for c in cs if c.n_scored]
            agg.mae = sum(c.mae * k for c, k in w) / n_scored
            agg.bias = sum(c.bias * k for c, k in w) / n_scored
            agg.rmse = math.sqrt(sum(c.rmse**2 * k for c, k in w) / n_scored)
            agg.within = sum(c.within * k for c, k in w) / n_scored
        out.append(agg)
    return out


def _pct(v: float) -> str:
    return f"{'-':>7}" if math.isnan(v) else f"{v * 100:6.1f}%"


def render_conventions(cells: list[Cell]) -> list[str]:
    """How each scored series' measurement convention met its set's.

    A range check cannot see a convention mismatch -- every coordinate is
    in range, the quantity is simply not the same one. So it is resolved
    explicitly and reported, including the refusals (#389).
    """
    from collections import Counter

    counts = Counter(
        c.convention for c in cells if c.convention != conventions.STATUS_NOT_SCORED
    )
    if not counts:
        return []
    out = ["Measurement conventions, series against the set scoring it:"]
    for status in (
        conventions.STATUS_DIRECT,
        conventions.STATUS_CONVERTED,
        conventions.STATUS_UNDECLARED,
        conventions.STATUS_UNREGISTERED,
        conventions.STATUS_INCOMPATIBLE,
    ):
        if counts.get(status):
            out.append(f"  {status:<14} {counts[status]:>4}")
    seen: set[tuple[str, str]] = set()
    for c in cells:
        key = (c.label.split("  ")[0], c.convention_detail)
        if c.convention == conventions.STATUS_CONVERTED and key not in seen:
            seen.add(key)
            out.append(f"    converted: {key[0]:<40} {c.convention_detail}")
        elif c.convention in (
            conventions.STATUS_UNREGISTERED,
            conventions.STATUS_INCOMPATIBLE,
        ):
            out.append(f"    REFUSED:   {key[0]:<40} {c.convention_detail[:90]}")
    out.append(
        "  A refusal is a CATEGORY error, not a disagreement: no range check "
        "can see it, so it is declared."
    )
    return out


# Correlation sets whose accuracy is carried in C++. The others (film,
# orifice, effusion) are not `RibCorrelationSet`s and state no band here.
def _rib_set(name: str):
    # The rib runner's own registry, so a set added there is reported here
    # without a second list to keep in step.
    from validation.cooling.runner import SETS

    factory = SETS.get(name)
    return factory() if factory else None


def render_set_accuracy(set_names) -> list[str]:
    """What each set claims about its OWN accuracy, and whether that is a
    band anything may be judged against.

    #389 asked for the sets' `accuracy_R`/`accuracy_G` to be used as the
    source's own band. Reading them showed that cannot be done as stated:
    only `han_1988_orthogonal`'s are the author's claim. The other two are
    THIS PROJECT's measurements -- `han_park_1988_angled`'s taken through
    `evaluate_rib`, the very code path they would be judging. Judging a
    model against its own error is the circularity #415 removed from the
    dataset's `uncertainty`; it was still here, one level down.

    So they are REPORTED with their provenance and used as a band by
    nothing. `usable_as_band()` is the C++ predicate that says which.
    """
    import combaero as cb

    rows = []
    for name in set_names:
        s = _rib_set(name or "")
        if s is None:
            continue
        for label, acc in (("R", s.accuracy_R), ("G", s.accuracy_G)):
            if acc.provenance == cb.AccuracyProvenance.Unstated:
                rows.append((name, label, "  --  ", "unstated", ""))
                continue
            rows.append((
                name,
                label,
                f"{acc.value:.1%}",
                "STATED" if acc.usable_as_band() else "measured",
                "usable as a band" if acc.usable_as_band()
                else "this project's own measurement -- NOT a band",
            ))
    if not rows:
        return []
    out = ["Correlation-set accuracy, as each set declares it:"]
    for name, label, value, prov, note in rows:
        out.append(f"  {name:<26} accuracy_{label}  {value:>6}  "
                   f"{prov:<8} {note}")
    out.append(
        "  A `measured` figure is the error of the model it describes, so "
        "nothing is judged against it (#389)."
    )
    return out


def multi_set_rows(cells: list[Cell]) -> dict[str, list[tuple[str, Cell | None, Cell | None]]]:
    """Series scored by more than one set: per series, per set, the
    in-domain and out-of-domain cells. Sets in name order -- never by score."""
    out: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(lambda: [None, None]))
    for c in cells:
        if not c.multi_set or c.scored_by is None:
            continue
        out[c.series_label][c.scored_by][0 if c.domain == "in-domain" else 1] = c
    return {
        series: [(name, *by_set[name]) for name in sorted(by_set)]
        for series, by_set in sorted(out.items())
    }


def _merged(cs: list[Cell | None]) -> tuple[int, float]:
    """Scored points and bias over the cells given, weighted by points."""
    cs = [c for c in cs if c is not None and c.n_scored]
    n = sum(c.n_scored for c in cs)
    return n, (sum(c.bias * c.n_scored for c in cs) / n if n else float("nan"))


def render_multi_set(cells: list[Cell]) -> list[str]:
    """Every eligible set on every series that has more than one.

    The sets are listed in name order and none is marked best; the spread
    between them is reported as what it is -- the model-form uncertainty a
    designer inherits by choosing one. Eligibility is derived by rule
    (runner.eligible_sets), so this list cannot have been trimmed to the
    sets that score well.
    """
    rows = multi_set_rows(cells)
    if not rows:
        return []
    out = [
        "Series scored by several sets (every eligible set, in name order; "
        "none is ranked, none is picked):",
        f"  {'series / set':<52} {'all':>4} {'bias':>7}   {'in':>3} {'bias':>7}  lab",
    ]
    one_lab: set[str] = set()
    for series, sets in rows.items():
        out.append(f"  {series}")
        biases, labs = [], set()
        for name, inside, outside in sets:
            n_all, b_all = _merged([inside, outside])
            n_in, b_in = _merged([inside])
            lab = SET_LAB.get(name, "unknown")
            labs.add(lab)
            if n_all:
                biases.append(b_all)
            out.append(
                f"    {name:<50} {n_all:>4} {_pct(b_all)}   {n_in:>3} {_pct(b_in)}  {lab}"
            )
        if len(biases) > 1:
            out.append(
                f"    spread of bias across sets: {_pct(min(biases)).strip()} to "
                f"{_pct(max(biases)).strip()} -- the cost of choosing one"
            )
        if len(labs) == 1:
            one_lab |= labs
    out.append(
        "  all = every scored point; in = inside that set's own validity box."
    )
    for lab in sorted(one_lab - {"unknown"}):
        out.append(
            f"  Where every set is {lab}, the sets agreeing with each other is "
            "NOT independent corroboration: only the measured series is."
        )
    return out


def render(
    cells: list[Cell], pools: dict | None = None, dataset=None
) -> str:
    pools = pools or {}
    head = (
        f"{'series':<44} {'kind':<12} {'N':>3} {'scored':>6} "
        f"{'MAE':>7} {'RMSE':>7} {'bias':>7} {'within':>7} {'extrap':>6} "
        f"{'basis':<9} {'sampling':<9} {'bnd':>3} {'held':>7}"
    )
    lines = [head, "-" * len(head)]
    for c in cells:
        lines.append(
            f"{c.label:<44} {c.kind:<12} {c.n:>3} {c.n_scored:>6} "
            f"{_pct(c.mae)} {_pct(c.rmse)} {_pct(c.bias)} {_pct(c.within)} "
            f"{c.n_extrapolated:>6} {c.basis:<9} {c.sampling:<9} "
            f"{c.n_bounded if c.n_bounded else '':>3} {_pct(c.held)}"
        )
    summary = rollup(cells)
    if summary:
        lines.append("")
        lines.append("By correlation set:")
        lines.append("-" * len(head))
        for c in summary:
            lines.append(
                f"{c.label:<44} {c.kind:<12} {c.n:>3} {c.n_scored:>6} "
                f"{_pct(c.mae)} {_pct(c.rmse)} {_pct(c.bias)} {_pct(c.within)} "
                f"{c.n_extrapolated:>6} {c.basis:<9} {c.sampling:<9} "
                f"{c.n_bounded if c.n_bounded else '':>3} {_pct(c.held)}"
            )
        if any(c.sampling == "partial" for c in summary):
            lines.append("")
            lines.append(
                "  <partial>: only the marks a digitisation could separate are in "
                "these rows, and those are the ones"
            )
            lines.append(
                "  furthest from the cluster centre, so MAE/RMSE/within are an "
                "UPPER BOUND on the error, not an estimate."
            )
            lines.append(
                "  <unknown>: shares a panel with other classes, with nothing "
                "independent to say how many runs it should have."
            )
        lines.append(
            "  [fidelity]: the correlation's OWN paper -- a miss is our "
            "transcription. [accuracy]: cross-source, so a miss is the"
        )
        lines.append(
            "  model's limitation, not ours. Never pooled. See "
            "docs/VALIDATION_POLICY.md."
        )
        if any(c.domain == "out-of-domain" for c in summary):
            lines.append(
                "  {out}: scored OUTSIDE the correlation set's own declared "
                "validity box. Reported separately because a set"
            )
            lines.append(
                "  must not be judged on conditions its authors never "
                "claimed -- not because extrapolating is always worse."
            )
        multi = render_multi_set(cells)
        if multi:
            lines.append("")
            lines.extend(multi)
        lines.append("")
        lines.extend(fidelity.render(dataset))
        lines.append("")
        lines.extend(render_conventions(cells))
        lines.append("")
        lines.extend(render_set_accuracy(sorted({c.scored_by for c in summary})))
    unsupported = [c for c in cells if c.unsupported]
    if pools:
        lines.append("")
        lines.append("Pooled (the scored path where class attribution is provisional):")
        for name, (n, mae, rmse, bias, within) in sorted(pools.items()):
            lines.append(
                f"  {name:<44} {'':<12} {n:>3} {n:>6} "
                f"{_pct(mae)} {_pct(rmse)} {_pct(bias)} {_pct(within)}"
            )
    if unsupported:
        lines.append("")
        lines.append("Not scored:")
        for c in unsupported:
            lines.append(f"  {c.label:<44} {c.reason}")
    return "\n".join(lines)


def main() -> None:
    from validation.cooling.schema import load_dataset

    dataset = load_dataset()
    records = run_dataset(dataset)
    pools = {
        "han2012/fig4.46_R_*  (lower panel)": pool(records, "han2012/fig4.46_R_eD"),
        "han2012/fig4.46_G_*  (upper panel)": pool(records, "han2012/fig4.46_G_eD"),
    }
    print(render(build(records, dataset), pools, dataset))


if __name__ == "__main__":
    main()
