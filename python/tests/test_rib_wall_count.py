"""Friction against the number of ribbed walls: Han 2013 Figs. 4.88 / 4.90.

Chandra and Cook (1994, square, 1-4 ribbed walls) and Chandra, Niland and Han
(1997, 1:2 rectangle, ribs on wide and narrow walls), digitised from the
book and owner-approved 2026-11-12. Everything is cross-source accuracy for
Han's area weighting (#1106). The book prints no e/D, P/e or Re for these
rigs, so the element is scored only as a ratio to the same study's two-wall
case, as an interval over the set's validity box, conditional on transverse
ribs; see validation/cooling/rib_wall_count_runner.py.

Data are read through DATA_ROOT (CI runs pytest from python/tests).
"""

from __future__ import annotations

import csv
import statistics

import pytest

from validation.cooling import rib_wall_count_runner as wc
from validation.cooling.schema import DATA_ROOT, load_dataset, load_points

HAN = DATA_ROOT / "han2012"


@pytest.fixture(scope="module")
def readings() -> list[wc.Reading]:
    return wc.load_readings()


@pytest.fixture(scope="module")
def series():
    return [s for s in load_dataset() if s.source.name == "han2012" and s.figure in wc.FIGURES]


@pytest.fixture(scope="module")
def element_rows(readings) -> list[wc.ElementRow]:
    return wc.element_ratio(readings)


def _xs(name: str) -> list[float]:
    with open(HAN / name, newline="") as fh:
        rows = [r for r in csv.reader(fh) if r and r[0].strip() != "x"]
    return sorted(float(r[0]) for r in rows)


def _ys(name: str) -> list[float]:
    with open(HAN / name, newline="") as fh:
        rows = [r for r in csv.reader(fh) if r and r[0].strip() != "x"]
    return sorted(float(r[1]) for r in rows)


def test_every_flag_row_names_exactly_one_reading(series, readings) -> None:
    """41 flagged of 102 readings; each flag row matches one point by x and y
    (the index column is the pick order, the CSVs are sorted by x)."""
    flags = wc.load_flags()
    assert len(series) == 21
    assert len(readings) == 102
    assert len(flags) == 41
    assert {f["flag"] for f in flags} == {"stacked", "overlapped", "hidden", "split_merge"}
    names = {s.path.name: s for s in series}
    for f in flags:
        pts = load_points(names[f["series"]])
        hits = [
            p
            for p in pts
            if abs(p.x - float(f["x"])) < wc.MATCH_TOL and abs(p.y - float(f["y"])) < wc.MATCH_TOL
        ]
        assert len(hits) == 1, f
    assert sum(1 for r in readings if r.flag != wc.CLEAN) == 41


def test_reads_reproduce_the_ranges_printed_in_the_text() -> None:
    """p. 414 prints Chandra and Cook's case C and E ranges: the third channel."""
    c_f, e_f = _xs("fig4.88_lower_caseC_ribbed.csv"), _xs("fig4.88_lower_caseE_ribbed.csv")
    c_nu, e_nu = _ys("fig4.88_lower_caseC_ribbed.csv"), _ys("fig4.88_lower_caseE_ribbed.csv")
    for got, printed, tol in (
        (c_f[0], 4.35, 0.005),
        (c_f[-1], 6.29, 0.005),
        (e_f[0], 7.96, 0.005),
        (e_f[-1], 11.45, 0.005),
        (c_nu[-1], 2.64, 0.01),
        (c_nu[0], 1.92, 0.01),
        (e_nu[-1], 2.99, 0.01),
        (e_nu[0], 2.12, 0.01),
    ):
        assert abs(got / printed - 1.0) < tol, (got, printed)


def test_one_run_reads_the_same_friction_in_every_panel(readings) -> None:
    """Each run is plotted in each panel at one abscissa: 1.1% / 0.8% worst."""
    spread = wc.panel_spread(readings)
    assert spread["4.88"][1] == 15 and spread["4.90"][1] == 20
    assert spread["4.88"][0] < 0.015
    assert spread["4.90"][0] < 0.012


def test_no_rig_parameter_is_invented(series) -> None:
    """The book prints no e/D, P/e or Re: null, and the rib runner refuses."""
    from validation.cooling.runner import eligible_sets

    for s in series:
        assert s.geometry is None, s.label
        assert s.scores is None, s.label
        assert s.uncertainty is None, s.label
        assert s.rig["e_D"] is None and s.rig["p_e"] is None and s.rig["Re_range"] is None, s.label
        stated = s.figure == "4.88" and s.rig["case"] == "C"
        assert (s.rig["alpha_deg"] == 90.0) if stated else (s.rig["alpha_deg"] is None), s.label
        assert eligible_sets(s) == (), s.label


def test_printed_y_label_misprint_is_recorded(series) -> None:
    """Fig. 4.88 upper prints '5' at the tick read as 3; kept, flagged."""
    upper = [s for s in series if s.figure == "4.88" and s.panel == "upper"]
    assert upper
    for s in upper:
        assert s.verification["y_label_misprint"] == {"printed": 5, "read_as": 3}
        assert s.verification["y_ticks"] == [0, 1, 2, 3]


def test_area_weighting_predicts_every_wall_count_from_two_walls(readings) -> None:
    """Model-free, cross-source: case C predicts the others within 4% bias in
    both pools; the n-independent rule misses n < 4 by +32% to +322%."""
    rows = wc.area_weighting(readings)
    for figure in wc.FIGURES:
        for case in {r.reading.case for r in rows if r.reading.figure == figure} - {"A"}:
            for clean in (True, False):
                sub = [
                    r
                    for r in rows
                    if r.reading.figure == figure
                    and r.reading.case == case
                    and (r.reading.flag == wc.CLEAN) == clean
                ]
                if not sub:
                    continue
                bias = statistics.fmean(r.predicted / r.reading.f_ratio - 1.0 for r in sub)
                alt = statistics.fmean(r.alternative / r.reading.f_ratio - 1.0 for r in sub)
                assert abs(bias) < 0.04, (figure, case, clean, bias)
                if sub[0].phi < 1.0:
                    assert alt > 0.30, (figure, case, alt)


def test_element_ratio_rows_are_conditional_and_case_d_square_is_refused(element_rows) -> None:
    d = [r for r in element_rows if r.figure == "4.88" and r.case == "D"]
    assert len(d) == 1 and d[0].refused and "1, 2 or 4" in d[0].refused
    for r in element_rows:
        if r.refused is None:
            assert "alpha 90" in r.condition
    unpooled = [r for r in element_rows if not r.pooled and r.refused is None]
    assert [(r.figure, r.case, r.w_h) for r in unpooled] == [("4.90", "F", 0.5)]


def test_element_misses_narrow_ribbed_walls_across_the_whole_box(element_rows) -> None:
    """Conditional on transverse ribs: with ribs on the 3.81 cm walls (W/H 0.5,
    outside the set's 1-4) the element's ratio to case C is low at every point
    of the validity box, case E by 20-45%; the wide-wall cases bracket zero."""
    by = {(r.figure, r.case, r.w_h): r for r in element_rows if r.refused is None}
    e_clean = by[("4.90", "E", 0.5)].pools[wc.CLEAN]
    assert max(e_clean[1]) < -0.15
    for key in (("4.90", "B", 2.0), ("4.90", "F", 2.0), ("4.88", "E", 1.0)):
        biases = by[key].pools[wc.CLEAN][1]
        assert min(biases) < 0.0 < max(biases), key


def test_level_is_information_and_narrow_wall_cases_fall_outside(readings) -> None:
    out = {(f, c): (lo, hi, mlo, mhi) for f, c, lo, hi, mlo, mhi in wc.level_envelope(readings)}
    for key in (
        ("4.88", "B"),
        ("4.88", "C"),
        ("4.88", "E"),
        ("4.90", "B"),
        ("4.90", "C"),
        ("4.90", "F"),
    ):
        lo, hi, mlo, mhi = out[key]
        assert lo <= mlo and mhi <= hi, key
    for key in (("4.90", "D"), ("4.90", "E")):
        lo, hi, mlo, mhi = out[key]
        assert mhi > hi, key


def test_element_metric_moves_with_f_multiplier() -> None:
    """Falsify the metric: f_multiplier 1.1 on the scored configuration moves
    the pooled bias by 0.1 (1 + bias), +8.7 and +10.5 points here."""
    for figure, (base, pert) in wc.falsify(0.10).items():
        assert pert - base == pytest.approx(0.10 * (1.0 + base), abs=0.005), figure


def test_report_is_labelled_cross_source() -> None:
    text = "\n".join(wc.render())
    assert "ACCURACY (CROSS-SOURCE)" in text
    assert "fidelity" not in text.lower()
    assert "CONDITIONAL on transverse ribs" in text
    assert "INFORMATION, not a score" in text
    assert "clean" in text and "flagged" in text
