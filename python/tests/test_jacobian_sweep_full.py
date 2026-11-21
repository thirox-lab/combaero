"""Jacobian coverage sweep, slow tier: ``_jacobian_sweep.ci_cases`` against
finite differences, and the coverage gate.

A mismatch fails its case unless an open pattern covers it
(``_jacobian_sweep.OPEN``, each with a strict-xfail reproducer in
test_jacobian_sweep.py); the failure names the point, then ``worst``'s row,
row kind, column, scaled error and step. The cases are split by family into
classes so ``--dist loadscope`` spreads them over the workers. The whole sweep
(every seed and variant, 40 points a case) is the module's script:
``uv run python python/tests/_jacobian_sweep.py``.
"""

from __future__ import annotations

import pytest
from _jacobian_sweep import (
    BUILDERS,
    MAX_POINTS_CI,
    MAX_UNKNOWNS,
    UNRESOLVED_LIMIT,
    Case,
    Report,
    _small,
    ci_cases,
    network_classes,
    open_issue,
    run_case,
    sweep,
)

pytestmark = pytest.mark.slow

GROUPS = {
    "junctions": ["legacy_tee", "legacy_tee_ports", "tee3", "ideal_nport", "ejector", "splitter"],
    "ducts": ["random", "pin_fin", "series", "port_discharge", "machine"],
    "walls": ["hotcold", "hotcold_mild"],
    "cooling": ["impingement", "impingement_bypass", "effusion", "bleed_station", "jet_chain"],
    "combustion": ["fuel_controller", "combustor", "walled_combustor", "film"],
}


def _check(case: Case) -> None:
    if not _small(case):
        pytest.skip(f"over {MAX_UNKNOWNS} unknowns")
    rep = Report()
    run_case(case, rep, max_points=MAX_POINTS_CI)
    new = [f for f in rep.mismatches() if open_issue(f) is None]
    assert not new, "\n".join(f"{f.case} {f.point}: {f.message}" for f in new)
    # An unresolved row FD could have resolved (not round-off bound) fails
    # above UNRESOLVED_LIMIT: unresolved never passes silently.
    bad = rep.unresolved(UNRESOLVED_LIMIT)
    assert not bad, "\n".join(
        f"unresolved {f.error:.2e} (est {f.fd_estimate:.1e}) {f.element} [{f.row_kind}] col "
        f"{f.column} | {f.case} | {f.point}"
        for f in bad
    )


def _params(group: str) -> list[Case]:
    # Not filtered by size here: collection would build every network.
    return ci_cases(GROUPS[group], small_only=False)


def test_the_groups_cover_every_builder() -> None:
    assert sorted(f for g in GROUPS.values() for f in g) == sorted(BUILDERS)


class TestJunctions:
    @pytest.mark.parametrize("case", _params("junctions"), ids=Case.label)
    def test_matches_fd(self, case: Case) -> None:
        _check(case)


class TestDucts:
    @pytest.mark.parametrize("case", _params("ducts"), ids=Case.label)
    def test_matches_fd(self, case: Case) -> None:
        _check(case)


class TestWalls:
    @pytest.mark.parametrize("case", _params("walls"), ids=Case.label)
    def test_matches_fd(self, case: Case) -> None:
        _check(case)


class TestCooling:
    @pytest.mark.parametrize("case", _params("cooling"), ids=Case.label)
    def test_matches_fd(self, case: Case) -> None:
        _check(case)


class TestCombustion:
    @pytest.mark.parametrize("case", _params("combustion"), ids=Case.label)
    def test_matches_fd(self, case: Case) -> None:
        _check(case)


#: Exported classes never instantiated in a network.
ABSTRACT = {
    "NetworkElement": "abstract base of every element",
    "NetworkNode": "abstract base of every node",
    "MultiPortChamberBase": "base of the junctions (MPCE, ConstantKTee, IdealMomentumTee, ejector)",
    "TurbomachineElement": "base of CompressorElement and TurbineElement",
}
#: Classes whose rows never carry a composition partial, with the reason. A
#: listed class the sweep does reach that way fails the gate: the reason is
#: then wrong. (A FuelBoundary is not listed: its Tt_target row reads its
#: combustor's mix, and is counted there.)
NO_COMPOSITION = {
    "PressureBoundary": "a fixed state: nothing it receives mixes into it",
    "MassFlowBoundary": "a fixed state at the network's edge; one feeding a node mid-network "
    "is that node's stream, counted at the node",
    "WallNode": "a dead end on one element: at most one stream arrives, nothing mixes",
}
#: Classes that cannot run reversed. None: every class has a leg that can.
NO_REVERSAL: dict[str, str] = {}


class TestCoverageGate:
    def test_every_class_is_reached_reversed_and_mixing(self) -> None:
        """Over the slow tier's cases and points, every element and node
        class exported from combaero.network (and the ThermalWall) is
        reached at a compared point, with a leg reversed against the start,
        and where its composition partials are non-zero (``_record_coverage``:
        a node mixing differing streams, an element or wall touching one)."""
        rep = sweep(ci_cases(), max_points=MAX_POINTS_CI, coverage_only=True)
        missing = []
        for name in sorted(network_classes()):
            c = rep.coverage.get(name)
            if name in ABSTRACT:
                assert c is None, f"{name} is reached: not abstract"
                continue
            if c is None or not c.reached:
                missing.append(f"{name}: never reached")
                continue
            for what, n, allow in (
                ("reversed", c.reversed, NO_REVERSAL),
                ("with differing composition", c.composition, NO_COMPOSITION),
            ):
                if name in allow:
                    assert not n, f"{name} is reached {what} ({n}x): {allow[name]!r} is wrong"
                elif not n:
                    missing.append(f"{name}: never {what}")
        assert not missing, "\n".join(missing)
