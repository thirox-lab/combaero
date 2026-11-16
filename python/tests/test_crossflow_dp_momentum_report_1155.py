"""``dP_momentum`` reports the momentum terms the row applies (#1155).

#1150 made a stationed segment fed from a node at rest with no ``entry_K``
the one with ``entry_K = 0``, and #1188 made their compressible rows the
same dict. The compressible report still summed the entry term's full
(1 + K) q, which the mass-flow row never adds (its dynamic head is the
outlet closure's, only the loss K/(1 + K) of it goes on the outlet Pt):
1759.4 Pa for ``entry_K = 0`` against 3.6e-14 for no ``entry_K`` (boundary
feed), 1663.9 Pa against -3.3e-13 (plenum feed). The report now applies the
row's own share of each term (``_applied_momentum_drop``), so the two report
the same in every feed at rest, both regimes, both directions; a set
``entry_K`` reports its loss K q.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

from combaero.network import NetworkSolver

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_crossflow_boundary_feed_1150 as t  # noqa: E402

pytestmark = pytest.mark.slow


def _solve(g: object) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(g).solve()
    assert r["__success__"], r.get("__message__")
    return r


def _report(r: dict) -> dict[str, float | str]:
    """Every scalar the solve reports, element diagnostics included."""
    out: dict[str, float | str] = {}
    for k, v in r.items():
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            out[k] = float(v)
    for eid, d in r["__element_diag__"].items():
        for k, v in d.items():
            out[f"{eid}.{k}"] = float(v) if isinstance(v, (int, float)) else repr(v)
    return out


# Feeds at rest (#1150's rule) and, compressible, a moving chamber: there
# the terminal from-end takes the chamber's Pt with or without entry_K, so
# the rows are the same too. Incompressible a chamber feed is outside the
# rule (#1152).
CASES = [
    (regime, feed, bleed)
    for regime in t.REGIMES
    for feed in ("boundary", "plenum", "chamber")
    for bleed in (False, True)
    if not (feed == "chamber" and regime == "incompressible")
]


@pytest.mark.parametrize("reversed_", [False, True])
@pytest.mark.parametrize(("regime", "feed", "bleed"), CASES)
def test_no_entry_K_reports_entry_K_zero(
    regime: str, feed: str, bleed: bool, reversed_: bool
) -> None:
    """Every reported scalar, bit for bit. Compressible forward x reported
    1759.4 Pa (boundary) and 1663.9 Pa (plenum, chamber) with entry_K = 0
    against ~1e-13 with none; now both read the latter."""
    a = _report(_solve(t._selfst(regime, reversed_, None, feed=feed, bleed=bleed)))
    b = _report(_solve(t._selfst(regime, reversed_, 0.0, feed=feed, bleed=bleed)))
    assert a == b


@pytest.mark.parametrize("regime", t.REGIMES)
def test_no_entry_K_reports_entry_K_zero_with_a_hot_feed(regime: str) -> None:
    """The 900 K reservoir with a 300 K side stream at the station: the
    compressible x reported 1405.7 Pa against 367.5 Pa."""
    assert _report(_solve(t._hot(regime))) == _report(_solve(t._hot(regime, 0.0)))


def _states(g: object, r: dict, eid: str) -> tuple:
    g.resolve_all_topology()
    s = NetworkSolver(g)
    s._build_x0()
    x = np.array(r["__x_solution__"], dtype=float)
    el = g.elements[eid]
    a = s._get_node_state(g.nodes[el.from_node], x)
    b = s._get_node_state(g.nodes[el.to_node], x)
    a.m_dot = b.m_dot = float(r[f"{eid}.m_dot"])
    return el, a, b, s._network_flows(el, x)


@pytest.mark.parametrize("regime", t.REGIMES)
def test_a_set_entry_K_reports_what_its_row_applies(regime: str) -> None:
    """entry_K = 0.5: incompressible the row subtracts the full (1 + K) q,
    compressible it adds K/(1 + K) of it (the loss K q) to the outlet Pt.
    The report follows the row in both."""
    K = 0.5
    g = t._selfst(regime, False, K)
    r = _solve(g)
    el, a, b, flows = _states(t._selfst(regime, False, K), r, "x")
    terms = {k[5]: k[0].dP for k in el._momentum_terms(a, b, flows)}
    assert terms["entry"] > 1000.0  # the case is a real test
    share = 1.0 if regime == "incompressible" else K / (1.0 + K)
    expected = share * terms["entry"] + terms["to"]
    assert r["__element_diag__"]["x"]["dP_momentum"] == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("regime", t.REGIMES)
def test_the_report_is_continuous_as_entry_K_goes_to_zero(regime: str) -> None:
    """entry_K = 1e-9 reports entry_K = 0 to O(K q). Compressible it was
    (1 + K) q at any K > 0 and 0 at K = 0 with no entry_K, a step of q."""
    small = _solve(t._selfst(regime, False, 1e-9))["__element_diag__"]["x"]["dP_momentum"]
    zero = _solve(t._selfst(regime, False, None))["__element_diag__"]["x"]["dP_momentum"]
    assert small == pytest.approx(zero, abs=1e-5)
