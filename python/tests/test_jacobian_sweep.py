"""Jacobian coverage sweep, fast tier: fixed points of ``_jacobian_sweep``.

Each point is one the full sweep (test_jacobian_sweep_full.py, ``slow``)
reaches. The guards each caught a past bug on the commit before its fix and
are clean on main; the open patterns (``_jacobian_sweep.OPEN``) are strict
xfails, one reproducer each. The sweep's docstring has the method and its
exclusion rules.

* #1233 (its first push, 24067756): a branching legacy tee with its branch
  leg reversed (pattern 10 at ratio 3: the straight run reverses too) and the
  boundaries' gases redrawn: 1.2e-1.
* #1235 (main before #1233, 39d7451f): a forward legacy tee with a side
  stream into a port: 1.3e-5 on the side stream's column at the start
  vector; a tee with port chambers, two side streams, a leg reversed: 1.3e-4.
  Both points also fail on 24067756.

The metric is falsified on its own: a dropped or halved entry, a spurious
entry where the true partial is 0, a bug only an equilibrium combustor has,
a wrong entry at an exact zero, and one confined to tiny flows (round-off
bound at 1e-9 kg/s, caught at the always-kept 3e-5 kg/s point) are each a
mismatch no open issue covers.

The open reproducers on a Fanno network (over a second each) carry the
``slow`` mark; the rest of the module runs in about two seconds.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pytest
from _jacobian_sweep import OPEN, Case, Report, evaluate, open_issue, point

from combaero.network import NetworkSolver


def _check(case: Case, desc: str, **kw: object) -> Report:
    x0, pt = point(case, desc)
    rep = Report()
    evaluate(case, pt, x0, rep, **kw)
    assert rep.compared == 1, (case.label(), desc, rep.excluded, rep.nonfinite)
    return rep


#: Points that caught a past bug before its fix, or reach a pattern no
#: hand-written test does; clean on main.
GUARDS = [
    pytest.param(
        Case("legacy_tee", 2, 1, 0), "tee pattern 10 ratio 3", id="1233-reversed-leg-gases"
    ),
    pytest.param(Case("legacy_tee", 2, 2, 1), "start vector", id="1235-side-stream-forward-tee"),
    pytest.param(
        Case("legacy_tee_ports", 3, 3, 2), "tee pattern 10 ratio 1", id="1235-tee-ports-reversed"
    ),
    pytest.param(Case("legacy_tee_ports", 0, 3, 2), "pc flows +-1e-09", id="1238-port-near-zero"),
    pytest.param(Case("tee3", 0, 1, 0), "jct pattern 001 ratio 1", id="mpce-arm-reversed-gases"),
    pytest.param(Case("film", 0, 1, 0), "h pattern 1 ratio 1", id="film-reversed-gases"),
]


@pytest.mark.parametrize(("case", "desc"), GUARDS)
def test_a_guard_point_matches_its_finite_difference(case: Case, desc: str) -> None:
    rep = _check(case, desc)
    assert not rep.mismatches(), "\n".join(f.message for f in rep.mismatches())


def _scale_entry(
    row_owner: str, column: str, factor: float
) -> Callable[[NetworkSolver, np.ndarray], None]:
    def perturb(s: NetworkSolver, J: np.ndarray) -> None:
        j = s.unknown_names.index(column)
        rows = [k for k, (o, _) in enumerate(s._row_kinds(s._build_x0())) if o == row_owner]
        i = max(rows, key=lambda k: abs(J[k, j]))  # the owner's row that reads the column
        J[i, j] *= factor

    return perturb


@pytest.mark.parametrize("factor", [0.0, 0.5])
def test_the_sweep_flags_a_wrong_entry(factor: float) -> None:
    """Falsify the metric: the #1235 composition partial dropped, or halved,
    in the analytic J is a mismatch named by row, row kind and column. The
    entry is 4.1e-6 of its row (scaled), so half of it is 2.0e-6, twice TOL."""
    rep = _check(
        Case("legacy_tee", 2, 2, 1),
        "start vector",
        perturb_J=_scale_entry("TeeJunctionElement 'tee'", "ss0.m_dot", factor),
    )
    [f] = rep.mismatches()
    assert (f.element, f.row_kind, f.column) == ("TeeJunctionElement 'tee'", "p", "ss0.m_dot")
    assert "TeeJunctionElement 'tee' [p]: scaled" in f.message and "ss0.m_dot step" in f.message


@pytest.mark.parametrize("size", [1e-5, 0.1])
def test_a_spurious_entry_is_a_mismatch(size: float) -> None:
    """A wrong entry where the true partial is 0 (the tee row on ps.Pt; the
    FD reads exact zeros) is a mismatch, judged against the row scale, not
    against the FD value: an estimate relative to |C| called it unresolved."""

    def spurious(s: NetworkSolver, J: np.ndarray) -> None:
        x = s._build_x0()
        j = s.unknown_names.index("ps.Pt")
        sc = np.maximum(np.abs(x), 1e-12)
        i = next(
            k
            for k, (o, _) in enumerate(s._row_kinds(x))
            if o == "TeeJunctionElement 'tee'" and J[k, j] == 0.0
        )
        J[i, j] = size * np.max(np.abs(J[i] * sc)) / sc[j]

    rep = _check(Case("legacy_tee", 2, 2, 1), "start vector", perturb_J=spurious)
    [f] = rep.mismatches()
    assert (f.element, f.column) == ("TeeJunctionElement 'tee'", "ps.Pt")
    assert f.error == pytest.approx(size, rel=1e-3) and open_issue(f) is None


@pytest.mark.parametrize("factor", [2.0, -1.0])
def test_a_bug_confined_to_tiny_flows_is_caught_at_the_kept_small_flow(factor: float) -> None:
    """At 1e-9 kg/s the tee's own column is round-off bound, so an error
    there hides; the always-kept point at 3e-5 kg/s resolves it."""
    rep = _check(
        Case("legacy_tee", 0, 3, 2),
        "tee flows 3e-05",
        perturb_J=_scale_entry("TeeJunctionElement 'tee'", "tee.m_dot_com", factor),
    )
    new = [f for f in rep.mismatches() if open_issue(f) is None]
    assert [(f.element, f.column) for f in new] == [("TeeJunctionElement 'tee'", "tee.m_dot_com")]


def _halve_orifice_rows(s: NetworkSolver, J: np.ndarray) -> None:
    """Every OrificeElement row's partials halved, in a network with an
    equilibrium combustor only: a bug no other network has."""
    if not any(getattr(n, "method", None) == "equilibrium" for n in s.network.nodes.values()):
        return
    for i, (owner, _) in enumerate(s._row_kinds(s._build_x0())):
        if owner.startswith("OrificeElement"):
            J[i] *= 0.5


@pytest.mark.parametrize(
    ("case", "desc"),
    [
        (Case("combustor", 1, 0, 0), "start vector"),
        # A point #1244 covers: the bug must still be new there.
        (Case("combustor", 1, 0, 0), "a_in flows 3e-05"),
    ],
)
def test_a_bug_only_in_an_equilibrium_combustor_is_new(case: Case, desc: str) -> None:
    """The OPEN hits are narrow (row, column, point kind, a cap): halving
    every orifice row in an equilibrium-combustor network is new at a clean
    point and at a point an open issue already covers."""
    rep = _check(case, desc, perturb_J=_halve_orifice_rows)
    new = [f for f in rep.mismatches() if open_issue(f) is None]
    assert new, [f.message for f in rep.mismatches()]


def test_a_wrong_entry_at_an_exact_zero_is_new() -> None:
    """An exact zero is compared like any other point (#1245 fixed its own
    misses there): a 0.5 error on a compressor-network row there is new."""
    rep = _check(
        Case("machine", 0, 3, 2),
        "oi.m_dot = 0",
        perturb_J=_scale_entry("OrificeElement 'oo'", "oo.m_dot", 0.5),
    )
    new = [f for f in rep.mismatches() if open_issue(f) is None]
    assert [(f.element, f.column) for f in new] == [("OrificeElement 'oo'", "oo.m_dot")]


@pytest.mark.parametrize(
    ("case", "desc"),
    [
        pytest.param(
            o.case,
            o.point,
            marks=[
                pytest.mark.xfail(
                    strict=o.strict,
                    reason=f"{o.issue}: {o.what}"
                    + ("" if o.strict else " (borderline, platform-dependent)"),
                )
            ]
            + ([pytest.mark.slow] if o.slow else []),
            id=o.issue.lstrip("#"),
        )
        for o in OPEN
    ],
)
def test_an_open_mismatch_still_reproduces(case: Case, desc: str) -> None:
    rep = _check(case, desc)
    assert not rep.mismatches(), "\n".join(f.message for f in rep.mismatches())
