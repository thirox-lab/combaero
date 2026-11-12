"""Score the ribbed-channel friction against the number of ribbed walls.

Data: Han, Dutta and Ekkad (2013) Figs. 4.88 and 4.90, digitised -- after
Chandra and Cook (1994), a square channel with 1-4 ribbed walls, and Chandra,
Niland and Han (1997), a 1:2 rectangle with ribs on wide and narrow walls.
Each case's friction ratio f/f0 is read in each panel it is plotted in; those
READINGS are the points here (validation/cooling/data/han2012/metadata.yaml).

EVERYTHING HERE IS CROSS-SOURCE ACCURACY. The rule under test -- Han's area
weighting, fbar = phi f + (1 - phi) f0 with phi the ribbed share of the
perimeter (`cb.rib_channel_friction`, #1106) -- comes from Han and Park
(1988). Chandra et al. are the same lab, different studies
(docs/VALIDATION_POLICY.md), so no row is fidelity.

THE BOOK DOES NOT PRINT THE RIB. It gives the wall count, which walls are
ribbed and the 1:2 wall widths, and "transverse ribs" for Fig. 4.88 case C;
no e/D, P/e, other rib angle or Re range, and neither paper is held. So:

A. MODEL-FREE (scored). The same study's two-wall case C fixes the four-sided
   ribbed friction at each Re rank, F_i = 1 + (C_i - 1)/phi_C, and the rule
   then predicts every other case with no correlation and no geometry:
   pred = 1 + phi_n (F_i - 1). Pairing by Re rank assumes every case was run
   at the same five Re, which the book does not print. Falsified by the
   n-independent rule the element used before #1106 (four-sided f for any
   wall count).

B. THROUGH THE ELEMENT (conditional). The ribbed `ChannelElement` with
   han_1988_orthogonal, n_ribbed_walls and W_H from the printed walls. Its
   absolute level needs the unprinted e/D, P/e and Re, so it is NOT scored;
   `level_envelope` reports only whether the measured range falls inside what
   the element gives across the set's own validity box. What is scored is the
   ratio to the same study's two-wall case at the same rib and Re,
   el(n)/el(C) against meas(n)/meas(C), in which the unknown rib largely
   cancels. It is reported as an INTERVAL over that box (e/D, P/e and Re
   at the set's lo/mid/hi), never a value at a chosen point. It is
   CONDITIONAL on transverse ribs: stated for Fig. 4.88 case C only (B and E
   assumed to share its rib), not printed for Fig. 4.90.
   Three-wall case D is refused (the element accepts 1, 2 or 4 walls).

Flagged readings (stacked, overlapped, hidden, split_merge; the owner judged
their position error well below scoring significance) are pooled apart from
clean ones. Nothing is fitted, and no multiplier is applied except in
`falsify`, which perturbs f_multiplier on the scored configuration to show
the element metric moves.

    uv run python -m validation.cooling.rib_wall_count_runner
"""

from __future__ import annotations

import csv
import itertools
import statistics
from dataclasses import dataclass
from functools import lru_cache

from validation.cooling import runner
from validation.cooling.rib_element_friction_runner import element_fbar
from validation.cooling.schema import DATA_ROOT, SeriesMetadata, load_dataset, load_points

FIGURES = ("4.88", "4.90")
FLAGS_CSV = DATA_ROOT / "han2012" / "fig4.88_fig4.90_point_flags.csv"
ANCHOR_CASE = "C"
SET_NAME = "han_1988_orthogonal"
CLEAN = "clean"
# A flag row and a CSV point are the same reading when both coordinates
# agree to the flags file's 4 decimals.
MATCH_TOL = 6e-5
# Fig. 4.90 case F has all four walls ribbed, so no ribbed wall defines W/H.
# Declared before scoring: the set's envelope is W/H 1-4, so the pooled row
# uses the wide-wall orientation (2.0); 0.5 is reported beside it, unpooled.
CASE_F_W_H = 2.0
CASE_F_W_H_SENSITIVITY = 0.5


@dataclass(frozen=True)
class Reading:
    figure: str
    case: str
    series: str
    rank: int  # 0 = lowest friction ratio (lowest Re) in its series
    f_ratio: float
    flag: str


def _series() -> list[SeriesMetadata]:
    return [
        s
        for s in load_dataset()
        if s.source.name == "han2012" and s.figure in FIGURES and s.rig is not None
    ]


def load_flags() -> list[dict]:
    with open(FLAGS_CSV, newline="") as fh:
        return list(csv.DictReader(fh))


def load_readings() -> list[Reading]:
    """Every friction-ratio reading, each with its own confidence flag."""
    flags = load_flags()
    out: list[Reading] = []
    for s in _series():
        pts = sorted(load_points(s), key=lambda p: p.x)
        mine = [f for f in flags if f["series"] == s.path.name]
        for rank, p in enumerate(pts):
            hit = [
                f["flag"]
                for f in mine
                if abs(float(f["x"]) - p.x) < MATCH_TOL and abs(float(f["y"]) - p.y) < MATCH_TOL
            ]
            if len(hit) > 1:
                raise ValueError(f"{s.path.name}: {len(hit)} flag rows match one reading")
            out.append(
                Reading(s.figure or "", s.rig["case"], s.path.name, rank, p.x, hit[0] if hit else CLEAN)
            )
    return out


def case_rig(figure: str) -> dict[str, dict]:
    """Per case, the printed rig (identical across that case's series)."""
    out: dict[str, dict] = {}
    for s in _series():
        if s.figure != figure:
            continue
        rig = {k: v for k, v in s.rig.items() if k != "walls_plotted"}
        prev = out.setdefault(rig["case"], rig)
        if prev != rig:
            raise ValueError(f"Fig. {figure} case {rig['case']}: series disagree on the rig")
    return out


def anchor(readings: list[Reading], figure: str) -> list[float]:
    """Case C's friction ratio per Re rank, the mean over its readings."""
    by_rank: dict[int, list[float]] = {}
    for r in readings:
        if r.figure == figure and r.case == ANCHOR_CASE:
            if r.flag != CLEAN:
                raise ValueError("the two-wall anchor must rest on clean readings")
            by_rank.setdefault(r.rank, []).append(r.f_ratio)
    return [statistics.fmean(by_rank[k]) for k in sorted(by_rank)]


def panel_spread(readings: list[Reading] | None = None) -> dict[str, tuple[float, int]]:
    """Per figure: the largest relative spread between readings of one run.

    Every run is plotted in each panel at the same friction ratio, so the
    readings of one (case, Re rank) are independent picks of one value.
    Returns (max spread relative to the mean, runs compared). Hidden marks
    whose x was taken from the covering series are not independent; they
    agree by construction and only make this check easier to pass.
    """
    readings = readings if readings is not None else load_readings()
    out: dict[str, tuple[float, int]] = {}
    for figure in FIGURES:
        runs: dict[tuple[str, int], list[float]] = {}
        for r in readings:
            if r.figure == figure:
                runs.setdefault((r.case, r.rank), []).append(r.f_ratio)
        multi = [v for v in runs.values() if len(v) > 1]
        worst = max((max(v) - min(v)) / statistics.fmean(v) for v in multi)
        out[figure] = (worst, len(multi))
    return out


def _stats(errs: list[float]) -> tuple[int, float, float]:
    if not errs:
        return 0, float("nan"), float("nan")
    return len(errs), statistics.fmean(errs), statistics.fmean(abs(e) for e in errs)


# --------------------------------------------------------------------------
# A. Model-free area weighting
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class AreaRow:
    reading: Reading
    phi: float
    predicted: float  # area weighting from case C
    alternative: float  # n-independent rule (pre-#1106 element)


def area_weighting(readings: list[Reading] | None = None) -> list[AreaRow]:
    readings = readings if readings is not None else load_readings()
    rows: list[AreaRow] = []
    for figure in FIGURES:
        rigs = case_rig(figure)
        c = anchor(readings, figure)
        phi_c = float(rigs[ANCHOR_CASE]["phi"])
        for r in readings:
            if r.figure != figure or r.case == ANCHOR_CASE:
                continue
            phi = float(rigs[r.case]["phi"])
            if phi == 0.0:
                pred = alt = 1.0  # all smooth: the rule says f = f0
            else:
                four_sided = 1.0 + (c[r.rank] - 1.0) / phi_c
                pred = 1.0 + phi * (four_sided - 1.0)
                alt = four_sided
            rows.append(AreaRow(r, phi, pred, alt))
    return rows


# --------------------------------------------------------------------------
# B. Through the element, conditional on transverse ribs
# --------------------------------------------------------------------------


def _box() -> list[tuple[float, float, float]]:
    """han_1988_orthogonal's own validity box at lo / mid / hi, nothing chosen."""
    s = runner.SETS[SET_NAME]()

    def three(rng, geometric: bool = False) -> tuple[float, float, float]:
        mid = (rng.lo * rng.hi) ** 0.5 if geometric else 0.5 * (rng.lo + rng.hi)
        return (rng.lo, mid, rng.hi)

    return list(itertools.product(three(s.valid_eD), three(s.valid_pe), three(s.valid_Re, True)))


@lru_cache(maxsize=None)
def _el(n: int, w_h: float, e_d: float, p_e: float, re: float, f_mult: float) -> float:
    rib_set = runner.SETS[SET_NAME]()
    geom = runner._probe_geometry(rib_set)
    geom.e_D, geom.p_e, geom.alpha_deg, geom.W_H = e_d, p_e, 90.0, w_h
    return element_fbar(rib_set, geom, n, re, f_mult)


def _config(rig: dict, w_h_override: float | None = None) -> tuple[int, float] | None:
    n = int(rig["n_ribbed_walls"])
    if n not in (1, 2, 4):
        return None
    w_h = w_h_override if w_h_override is not None else rig["W_H"]
    if w_h is None:
        w_h = CASE_F_W_H
    return n, float(w_h)


@dataclass(frozen=True)
class ElementRow:
    figure: str
    case: str
    label: str
    pooled: bool
    condition: str
    n_ribbed: int | None
    w_h: float | None
    refused: str | None
    # pool -> (n readings, bias per box point, MAE per box point), the box
    # points in `_box()` order so rows pool point by point
    pools: dict


def _condition(figure: str, rig: dict) -> str:
    if figure == "4.88" and rig.get("alpha_deg") is not None:
        return "alpha 90 stated (p. 414)"
    if figure == "4.88":
        return "alpha 90 assumed: case C's stated rib"
    return "alpha 90 assumed: not printed"


def element_ratio(
    readings: list[Reading] | None = None, f_mult: float = 1.0
) -> list[ElementRow]:
    """el(n)/el(C) against meas(n)/meas(C), as an interval over the box.

    `f_mult` multiplies the SCORED configuration's friction only, never the
    two-wall reference, so a perturbation shows in the ratio.
    """
    readings = readings if readings is not None else load_readings()
    box = _box()
    rows: list[ElementRow] = []
    for figure in FIGURES:
        rigs = case_rig(figure)
        c = anchor(readings, figure)
        ref = _config(rigs[ANCHOR_CASE])
        assert ref is not None
        for case in sorted(rigs):
            if case in (ANCHOR_CASE, "A"):
                continue
            rig = rigs[case]
            variants: list[tuple[str, bool, float | None]] = [("", True, None)]
            if rig["n_ribbed_walls"] == 4 and rig["W_H"] is None:
                variants = [
                    (f"W/H {CASE_F_W_H:g}", True, CASE_F_W_H),
                    (f"W/H {CASE_F_W_H_SENSITIVITY:g} (sensitivity)", False, CASE_F_W_H_SENSITIVITY),
                ]
            for suffix, pooled, w_override in variants:
                cfg = _config(rig, w_override)
                label = f"case {case}" + (f", {suffix}" if suffix else "")
                if cfg is None:
                    rows.append(
                        ElementRow(figure, case, label, False, "", None, None,
                                   f"{rig['n_ribbed_walls']} ribbed walls: the element accepts 1, 2 or 4", {})
                    )
                    continue
                mine = [r for r in readings if r.figure == figure and r.case == case]
                per_pool: dict[str, list[list[float]]] = {}
                for e_d, p_e, re in box:
                    el_ratio = _el(*cfg, e_d, p_e, re, f_mult) / _el(*ref, e_d, p_e, re, 1.0)
                    for pool in (CLEAN, "flagged"):
                        sub = [r for r in mine if (r.flag == CLEAN) == (pool == CLEAN)]
                        errs = [el_ratio / (r.f_ratio / c[r.rank]) - 1.0 for r in sub]
                        per_pool.setdefault(pool, []).append(errs)
                pools = {}
                for pool, err_sets in per_pool.items():
                    k = len(err_sets[0])
                    if k == 0:
                        continue
                    biases = [statistics.fmean(e) for e in err_sets]
                    maes = [statistics.fmean(abs(x) for x in e) for e in err_sets]
                    pools[pool] = (k, biases, maes)
                rows.append(
                    ElementRow(figure, case, label, pooled, _condition(figure, rig), cfg[0], cfg[1], None, pools)
                )
    return rows


def pooled_element(rows: list[ElementRow], figure: str, pool: str) -> tuple[int, float, float] | None:
    """Pooled bias over the box for one figure and pool: (n, lo, hi).

    Readings-weighted at each box point, then the extremes over the box.
    """
    sub = [r.pools[pool] for r in rows if r.figure == figure and r.pooled and pool in r.pools]
    if not sub:
        return None
    k = sum(s[0] for s in sub)
    at_point = [sum(s[0] * s[1][g] for s in sub) / k for g in range(len(sub[0][1]))]
    return k, min(at_point), max(at_point)


def level_envelope(readings: list[Reading] | None = None) -> list[tuple[str, str, float, float, float, float]]:
    """Information, NOT a score: element f/f0 over the box vs the measured range.

    f0 is taken as 0.046 Re^-0.2 (Han and Park's smooth reference); the book
    does not define Chandra's baseline, which is one more reason this is not
    scored.
    """
    readings = readings if readings is not None else load_readings()
    out = []
    for figure in FIGURES:
        for case, rig in sorted(case_rig(figure).items()):
            cfg = _config(rig)
            if cfg is None or case == "A":
                continue
            vals = [
                _el(*cfg, e_d, p_e, re, 1.0) / (runner.SMOOTH_F_COEF * re**runner.SMOOTH_F_EXP)
                for e_d, p_e, re in _box()
            ]
            meas = [r.f_ratio for r in readings if r.figure == figure and r.case == case]
            out.append((figure, case, min(vals), max(vals), min(meas), max(meas)))
    return out


def falsify(delta: float = 0.10) -> dict[str, tuple[float, float]]:
    """Pooled clean element bias (box midpoint of the interval) at f_multiplier
    1 and 1 + delta on the scored configuration; the shift should be ~delta."""
    readings = load_readings()
    base = element_ratio(readings, 1.0)
    pert = element_ratio(readings, 1.0 + delta)
    out = {}
    for figure in FIGURES:
        b = pooled_element(base, figure, CLEAN)
        p = pooled_element(pert, figure, CLEAN)
        if b and p:
            out[figure] = (0.5 * (b[1] + b[2]), 0.5 * (p[1] + p[2]))
    return out


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------


def _pool_line(label: str, errs: list[float]) -> str:
    k, bias, mae = _stats(errs)
    if k == 0:
        return f"    {label:<34} n=0"
    return f"    {label:<34} n={k:<3} bias {bias:+6.1%}  MAE {mae:5.1%}"


def render() -> list[str]:
    readings = load_readings()
    n_flag = sum(1 for r in readings if r.flag != CLEAN)
    lines = [
        "Ribbed-channel friction vs number of ribbed walls -- ACCURACY (CROSS-SOURCE)",
        "Han 2013 Figs. 4.88 (after Chandra and Cook 1994) and 4.90 (after Chandra, Niland and Han 1997);",
        f"{len(readings)} friction readings, {len(readings) - n_flag} clean, {n_flag} flagged (pooled apart).",
        "Cross-panel check, one run read in each panel: "
        + "; ".join(
            f"Fig. {fig} max spread {w:.1%} over {k} runs" for fig, (w, k) in panel_spread(readings).items()
        ),
        "",
        "A. Model-free area weighting from the same study's two-wall case C (by Re rank)",
    ]
    rows = area_weighting(readings)
    for figure in FIGURES:
        lines.append(f"  Fig. {figure}")
        lines.append(f"    {'case':<5} {'phi':>6} {'pool':<8} {'n':>3} {'bias':>7} {'MAE':>6}   n-independent bias")
        for case in sorted({r.reading.case for r in rows if r.reading.figure == figure}):
            for pool in (CLEAN, "flagged"):
                sub = [
                    r for r in rows
                    if r.reading.figure == figure and r.reading.case == case
                    and (r.reading.flag == CLEAN) == (pool == CLEAN)
                ]
                if not sub:
                    continue
                k, bias, mae = _stats([r.predicted / r.reading.f_ratio - 1.0 for r in sub])
                _, alt, _ = _stats([r.alternative / r.reading.f_ratio - 1.0 for r in sub])
                tag = "  (all-smooth baseline, not pooled)" if case == "A" else ""
                lines.append(
                    f"    {case:<5} {sub[0].phi:6.3f} {pool:<8} {k:>3} {bias:+7.1%} {mae:6.1%}   {alt:+7.1%}{tag}"
                )
        for pool in (CLEAN, "flagged"):
            errs = [
                r.predicted / r.reading.f_ratio - 1.0
                for r in rows
                if r.reading.figure == figure and r.reading.case != "A"
                and (r.reading.flag == CLEAN) == (pool == CLEAN)
            ]
            lines.append(_pool_line(f"POOLED accuracy (cross-source) {pool}", errs))
    lines.append("")
    lines.append(
        "B. Through the element (han_1988_orthogonal), ratio to case C, CONDITIONAL on transverse ribs;"
    )
    lines.append(
        "   bias/MAE are intervals over the set's validity box (e/D, P/e, Re at lo/mid/hi), not a fitted point"
    )
    erows = element_ratio(readings)
    for figure in FIGURES:
        lines.append(f"  Fig. {figure}")
        for r in [r for r in erows if r.figure == figure]:
            if r.refused:
                lines.append(f"    {r.label:<30} refused: {r.refused}")
                continue
            for pool, (k, biases, maes) in r.pools.items():
                blo, bhi, mlo, mhi = min(biases), max(biases), min(maes), max(maes)
                lines.append(
                    f"    {r.label:<30} n={r.n_ribbed} W/H {r.w_h:<4g} {pool:<8} {k:>3}  "
                    f"bias {blo:+6.1%}..{bhi:+6.1%}  MAE {mlo:5.1%}..{mhi:5.1%}  [{r.condition}]"
                    + ("" if r.pooled else "  (not pooled)")
                )
        for pool in (CLEAN, "flagged"):
            p = pooled_element(erows, figure, pool)
            if p:
                lines.append(
                    f"    POOLED accuracy (cross-source, conditional) {pool:<8} n={p[0]:<3} "
                    f"bias {p[1]:+6.1%}..{p[2]:+6.1%}"
                )
    lines.append("")
    lines.append("   Level, element f/f0 over the box vs measured (INFORMATION, not a score: e/D, P/e, Re, f0 unprinted)")
    for figure, case, lo, hi, mlo, mhi in level_envelope(readings):
        inside = "inside" if lo <= mlo and mhi <= hi else "NOT inside"
        lines.append(
            f"    Fig. {figure} case {case}: element {lo:5.2f}..{hi:5.2f}  measured {mlo:5.2f}..{mhi:5.2f}  ({inside})"
        )
    lines.append("")
    lines.append("Falsification: f_multiplier 1.10 on the scored configuration (pooled clean, interval midpoint)")
    for figure, (b, p) in falsify().items():
        lines.append(f"    Fig. {figure}: bias {b:+.1%} -> {p:+.1%} (shift {p - b:+.1%})")
    return lines


def main() -> None:
    print("\n".join(render()))


if __name__ == "__main__":
    main()
