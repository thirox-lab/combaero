"""The hole-inclination Cd term against Hay 1983 and Gritsch 1998.

ACCURACY, cross-source: the term is Idelchik's, neither paper is. These pins
hold the numbers the runner reports (and the PR body quotes); nothing is
fitted, and the series stay unscored in the scorecard.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.cooling import hole_cd_inclination_runner as hr  # noqa: E402
from validation.cooling.schema import load_dataset  # noqa: E402


@pytest.fixture(scope="module")
def ds():
    return load_dataset()


def _bias(rows: list, which: str) -> float:
    return sum(r.err(which) for r in rows) / len(rows)


def test_only_no_coolant_crossflow_series_are_in_scope(ds) -> None:
    labels = {s.label for s in hr.in_scope(ds)}
    assert labels == {
        "gritsch1998/fig3_Mac0p0_Mam0p0",
        "hay1983/fig6_open_circles",
        "hay1983/fig8_open_circles",
        "hay1983/fig12_open_circles",
    }
    whys = dict(hr.out_of_scope(ds))
    assert whys["hay1983/fig5_open_circles"] == "coolant (inlet) crossflow"
    assert whys["gritsch1998/fig7_Mac0p0_Mam0p0"] == "shaped hole"


def test_gritsch_30_deg_no_crossflow(ds) -> None:
    (s,) = [s for s in hr.in_scope(ds) if s.source.name == "gritsch1998"]
    rows = hr.rows_for(s)
    assert _bias(rows, "none") == pytest.approx(0.064, abs=0.005)
    assert _bias(rows, "idelchik_3_2") == pytest.approx(-0.043, abs=0.005)


def test_hay_carpets_the_term_moves_only_the_30_deg_holes() -> None:
    rows = hr.carpet_rows()
    by = {r.label: [x for x in rows if x.label == r.label] for r in rows}
    for rs in by.values():
        if rs[0].alpha == 90.0:
            assert _bias(rs, "none") == _bias(rs, "idelchik_3_2")
        else:
            assert abs(_bias(rs, "idelchik_3_2")) < abs(_bias(rs, "none"))
    assert _bias(by["hay1983/carpet_fig10_parallel"], "idelchik_3_2") == pytest.approx(
        -0.006, abs=0.005
    )


def test_the_30_over_90_ratio(ds) -> None:
    (row,) = hr.angle_ratio(ds, (1.5,))
    assert row["none"] == pytest.approx(1.0, abs=1e-3)
    assert row["idelchik_3_2"] == pytest.approx(0.90, abs=0.01)
    assert row["measured"] == pytest.approx(0.904, abs=0.01)
