"""Kaddah (1964) liner-hole tables: loader, read-confidence flags, pinned scores.

The dataset is the thesis's printed Tables 11-30 (validation/cooling/data/
kaddah1964). Scores are ACCURACY, CROSS-SOURCE (docs/VALIDATION_POLICY.md):
the pins below guard against a silent change in the data, the conversion or
the model, not a target the model should meet.
"""

from __future__ import annotations

import csv
import inspect
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.cooling import conventions as C  # noqa: E402, N812
from validation.cooling import kaddah_runner as K  # noqa: E402, N812

COLUMNS = [
    "p1/pj",
    "alpha_%",
    "mu",
    "L1_lb/hr",
    "Lj_lb/hr",
    "rho1_lb/cu.ft",
    "V1_ft/sec",
    "rhoj_lb/cu.ft",
    "Vj_ft/sec",
    "Lth_lb/hr",
    "C",
    "K",
    "phi_deg",
    "delta",
]


@pytest.fixture(scope="module")
def tables():
    return K.load_tables()


def test_all_twenty_tables_load_as_printed(tables) -> None:
    assert [t.meta["table"] for t in tables] == list(range(11, 31))
    codes = [t.code for t in tables]
    assert codes == [f"C{i}" for i in range(1, 6)] + [f"O{i}" for i in range(1, 8)] + [
        f"R{i}" for i in range(1, 9)
    ]
    assert sum(len(t.rows) for t in tables) == 269
    for t in tables:
        assert {r.block for r in t.rows} == {"noxflow", "80", "130"}
        for r in t.rows:
            for c in COLUMNS:
                float(r.printed[c])  # every printed cell is a number


def test_flags_point_at_real_cells_with_the_committed_value(tables) -> None:
    by_key = {(t.meta["table"], r.block, r.row): r for t in tables for r in t.rows}
    with open(K.DATA / "cell_flags.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert rows
    for f in rows:
        r = by_key[(int(f["table"]), f["block"], int(f["row"]))]
        assert r.printed[f["column"]] == f["value_as_transcribed"], f
        assert f["flag"] in ("sequence_judged", "arithmetic_confirmed", "printed_inconsistent")
        assert f["alternative"], f


def test_c_is_lj_over_lth_except_where_flagged(tables) -> None:
    """C = Lj/Lth (App. 2, A2.1.4) on every row, to the printing precision,
    unless the row carries a printed_inconsistent flag -- which is how the
    transcription is guarded and how a new misprint gets noticed."""
    misses = []
    for t in tables:
        for r in t.rows:
            ratio = r.f("C") / (r.f("Lj_lb/hr") / r.f("Lth_lb/hr"))
            if abs(ratio - 1.0) > 0.02:
                misses.append((t.meta["table"], r.block, r.row, r.excluded))
    assert misses and all(excluded for *_, excluded in misses), misses


def test_kaddah_cd_is_declared_total_and_converted_not_compared_raw() -> None:
    assert C.SOURCE_PUBLISHES[("kaddah1964", "Cd")] == "Cd_total"
    assert K.static_factor(2.0) == pytest.approx(math.sqrt(2.0))


def test_both_u1_definitions_must_be_named() -> None:
    """Never one silently: run_all has no default U1 definition."""
    sig = inspect.signature(K.run_all)
    assert sig.parameters["u1"].default is inspect.Parameter.empty
    assert K.U1_DEFS == ("approach", "mean")
    row = K.load_tables()[0].rows[5]  # C1, 80 ft/s, row 1
    assert K.crossflow_ratio(row, "mean") == pytest.approx(
        (1.0 - row.alpha / 2.0) * K.crossflow_ratio(row, "approach")
    )


def test_exclusions_follow_the_flags(tables) -> None:
    excluded = [(r.table, r.block, r.row) for t in tables for r in t.rows if r.excluded]
    assert sorted(excluded) == [
        (13, "80", 1),
        (14, "80", 2),
        (16, "130", 4),
        (17, "130", 1),
        (24, "80", 3),
        (28, "130", 3),
    ]
    assert [t.code for t in tables if t.excluded] == ["R5"]
    assert all(rec.row.code != "R5" for rec in K.run_all("approach", "eq17"))


# Pinned (n, bias, MAE) for the liner-envelope configurations, clean inputs.
PINS = {
    ("eq17", "approach"): (54, -0.0533, 0.0715),
    ("eq17", "mean"): (54, -0.0299, 0.0582),
    ("chain", "approach"): (26, -0.0729, 0.0788),
    ("chain", "mean"): (26, -0.0674, 0.0758),
}


@pytest.mark.parametrize("mode,u1", list(PINS))
def test_liner_envelope_scores_are_pinned(mode: str, u1: str) -> None:
    recs = [
        r
        for r in K.run_all(u1, mode)
        if r.envelope == "inside" and not r.row.sequence_judged and r.predicted is not None
    ]
    s = K.score(recs)
    n, bias, mae = PINS[(mode, u1)]
    assert s.n == n
    assert s.bias == pytest.approx(bias, abs=5e-4)
    assert s.mae == pytest.approx(mae, abs=5e-4)


def test_eq17_reads_low_under_strong_crossflow() -> None:
    """The finding, pinned in sign: at K < 2 (strong crossflow) Eq. (17)
    under-predicts Kaddah's thin-wall Cd by 30-45% under either U1
    definition -- the TN 3663 direction (#807)."""
    for u1 in K.U1_DEFS:
        low = [r for r in K.run_all(u1, "eq17") if r.k < 2.0 and r.predicted is not None]
        assert low and K.score(low).bias < -0.25


def test_slots_are_refused_by_the_chain_not_scored() -> None:
    recs = K.run_all("approach", "chain")
    slots = [r for r in recs if not r.row.code.startswith("C")]
    assert slots and all(r.predicted is None and r.reason for r in slots)
