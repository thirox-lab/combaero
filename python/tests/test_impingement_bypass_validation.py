"""An impingement array with an initial crossflow, against the authors' own
flow data: Florschuetz, Metzger and Su (1982), NASA CR-3630, Table D.1 (#467).

42 tests: 11 inline geometries (xn/d 5, 10; yn/d 4, 8; z/d 1, 2, 3), each at
initial-to-jet flow ratios mc/mj of 0 to about 1, ten rows each. Per row the
table gives P/P0 (channel static over jet plenum), Gj/Gj_mean and Gc/Gj. The
jet flows are the authors' own, inferred from the measured P with a constant
plate Cd (their Eq. 4.1): that is the quantity this network computes, so the
comparison is like for like. Fidelity: the source's own data
(docs/VALIDATION_POLICY.md).

The rig is reproduced as built: jet plenum at the measured P0 and T0; the
total jet flow from Re_j; the initial crossflow a fixed flow delivered along
the gap's upstream end (CR-3630 ran it 16 to 95 hydraulic diameters before
the array, so it arrives with its momentum, which the lead-in carries); the
plate's own mean Cd from Table 3.1; channel width 18.29 cm (18 holes at
yn/d 4, d 0.254 cm).

What the network leaves out, and where it shows: the authors' own analysis
says the discharge coefficient falls with crossflow and wall shear matters
at z/d = 1 with large crossflow, and that no adequate friction model exists
for these flows. Those cases are reported, with a looser bound, not dropped.
"""

from __future__ import annotations

import csv
import math
import statistics
from collections import defaultdict
from functools import cache

import pytest

import combaero as cb
from combaero.network import (
    FlowNetwork,
    HoleCd,
    ImpingementArray,
    MassFlowBoundary,
    NetworkSolver,
    PressureBoundary,
)
from validation.cooling.schema import DATA_ROOT

TABLE = DATA_ROOT / "florschuetz_metzger_su1982" / "table_d1_flow.csv"
N_ROWS = 10
SPAN = 0.1829  # channel width [m]: 18 holes at yn/d 4, d 0.254 cm
# Table 3.1: hole diameter by streamwise pitch, mean plate Cd by (xn/d, yn/d).
D_JET = {5: 0.00254, 10: 0.00127}
PLATE_CD = {(5, 4): 0.85, (5, 8): 0.80, (10, 4): 0.76, (10, 8): 0.76}
Y_AIR = cb.mole_to_mass(cb.species.dry_air())
X_AIR = cb.species.dry_air()

# Over 20 s on one worker: runs in CI's parallel slow job.
pytestmark = pytest.mark.slow

Key = tuple[int, int, int, float]

# Tests the network does not solve from a cold start, each with its reason
# (a strict xfail below). None since #802: B(5,4,1) at mc/mj 0.99 stalled on a
# false Fanno choke at its exit segment, a station's static P fed to the
# kernel as P0 (99.7% of that choke; 81% of the true one). It now solves cold.
# Its stagnant rows 1-4 (data Gj/Gj_mean 0.00-0.16) came out with
# alternating signs under a constant Cd (#799); under the holes' crossflow
# Cd (#800) they are monotone, as the data are.
STALLED: dict[Key, str] = {}


@cache
def _table() -> dict[Key, list[dict[str, float]]]:
    out: dict[Key, list[dict[str, float]]] = defaultdict(list)
    with TABLE.open() as fh:
        for r in csv.DictReader(fh):
            key = (int(r["xn_d"]), int(r["yn_d"]), int(r["z_d"]), float(r["mc_mj"]))
            out[key].append({k: (float(v) if v else math.nan) for k, v in r.items()})
    return dict(out)


def _solved() -> list[Key]:
    """Every test but the STALLED ones, in table order."""
    return [k for k in sorted(_table()) if k not in STALLED]


@cache
def _solve(key: Key, with_bypass: bool = True) -> tuple[list[float], list[float], list[float]]:
    """Per row (Gj/Gj_mean, Gc/Gj, P/P0). ``with_bypass=False`` drops the
    initial crossflow, keeping the jets: the falsification."""
    xn, yn, z, mc = key
    rows = _table()[key]
    d = D_JET[xn]
    T0, P0, Re_j = rows[0]["T0_K"], rows[0]["P0_kPa"] * 1.0e3, rows[0]["Re_j"]
    mu = cb.complete_state(T0, P0, X_AIR).transport.mu
    n_holes = round(SPAN / (yn * d))
    m_j = Re_j * mu * N_ROWS * n_holes * math.pi * d / 4.0
    m_c = mc * m_j if with_bypass else 0.0

    net = FlowNetwork()
    supply = PressureBoundary("supply")
    supply.Pt, supply.Tt, supply.Y = P0, T0, Y_AIR
    net.add_node(supply)
    net.add_node(MassFlowBoundary("exit", m_dot=m_j + m_c, Tt=T0, Y=Y_AIR))
    arr = ImpingementArray(
        "a",
        n_rows=N_ROWS,
        d_jet=d,
        xn_d=xn,
        yn_d=yn,
        z_d=z,
        span=SPAN,
        plate_thickness=d,  # Table 3.1: b = d
        hole_cd=HoleCd(Cd=PLATE_CD[(xn, yn)]),
    )
    bypass = None
    if m_c > 0.0:
        # The metered initial crossflow, straight into the lead-in. The rig's
        # 24 cm development run only sets the feed pressure, which floats
        # here: at a fixed bypass flow it cannot move the split (see
        # test_impingement_bypass.py::test_the_feed_node_does_not_matter).
        net.add_node(MassFlowBoundary("feed", m_dot=m_c, Tt=T0, Y=Y_AIR))
        bypass = "feed"
    arr.add_to(net, "supply", "exit", bypass=bypass)
    # The default start is warm for impingement networks (#706).
    res = NetworkSolver(net).solve()
    assert res["__success__"], f"{key}: {res.get('__message__')}"
    diag = res["__element_diag__"]
    m = [diag[arr.plate_id(i)]["m_dot"] for i in range(1, N_ROWS + 1)]
    mean = sum(m) / N_ROWS
    return (
        [x / mean for x in m],
        [diag[arr.plate_id(i)]["Gc_Gj"] for i in range(1, N_ROWS + 1)],
        [res[f"{arr.node_id(i)}.P"] / P0 for i in range(1, N_ROWS + 1)],
    )


def _gj_err(key: Key, with_bypass: bool = True) -> float:
    """Worst row's |Gj/Gj_mean - measured| for one test."""
    gj, _, _ = _solve(key, with_bypass)
    return max(abs(g - r["Gj_Gjbar"]) for g, r in zip(gj, _table()[key], strict=True))


def _closed_form(xn: int, yn: int, z: int, mc: float) -> list[float]:
    """The authors' own model, CR-3630 Eq. 4.10 (constant Cd, no friction)."""
    B = N_ROWS * math.sqrt(2.0) * PLATE_CD[(xn, yn)] * (math.pi / 4.0) / (yn * z)
    return [
        B * ((1 + mc) * math.cosh(B * x) - mc * math.cosh(B * (1 - x))) / math.sinh(B)
        for x in ((i + 0.5) / N_ROWS for i in range(N_ROWS))
    ]


def test_the_transcription_closes_on_itself() -> None:
    """Table D.1 is parsed from the report's text layer. Two of its own
    identities catch a misread digit: the jet flows average to 1, and each
    Gc/Gj follows from the initial crossflow and the jets upstream of it,

        Gc_i/Gj_i = (pi/4)/((yn/d)(z/d)) (N (mc/mj)/MC + sum_{k<i} g_k) / g_i

    MC is the table's mass closure (inferred jet flows over the metered
    total, 0.91 to 1.02): the authors' Gc counts the metered initial
    crossflow against the inferred jets. Without it B(5,8,1) at mc/mj 0.52,
    MC 0.91, misses by 11% -- the printed 1.18, confirmed on the page."""
    table = _table()
    assert len(table) == 42
    for (xn, yn, z, mc), rows in table.items():
        g = [r["Gj_Gjbar"] for r in rows]
        assert len(g) == N_ROWS
        assert sum(g) / N_ROWS == pytest.approx(1.0, abs=0.003), (xn, yn, z, mc)
        closure = rows[0]["mass_closure"]
        for i, r in enumerate(rows):
            # 0.30: below that the printed two decimals dominate the ratio;
            # a few rows have Gj 0.0-0.02, where nothing can be checked.
            if r["Gc_Gj"] < 0.30 or g[i] < 0.10:
                continue
            want = (math.pi / 4) / (yn * z) * (N_ROWS * mc / closure + sum(g[:i])) / g[i]
            # Measured worst 4.4%, two-decimal rounding.
            assert r["Gc_Gj"] == pytest.approx(want, rel=0.06), (xn, yn, z, mc, i + 1)


@pytest.mark.parametrize(
    "key",
    [
        pytest.param(k, marks=pytest.mark.xfail(reason=STALLED[k], strict=True))
        if k in STALLED
        else k
        for k in sorted(_table())
    ],
    ids=str,
)
def test_each_test_solves(key: Key) -> None:
    """Every Table D.1 test reaches a root; the STALLED ones are expected not
    to, strictly, so a solve there is looked at rather than scored."""
    _solve(key)


def test_the_flow_distribution_matches_the_authors_data() -> None:
    """Row-by-row Gj/Gj_mean over all 42 tests (STALLED aside: none).

    Measured 2026-10-28 (#800, the holes' crossflow Cd, Rogers & Hersh
    1976 by default), worst row per test: z/d >= 2 median 0.021, max 0.085;
    z/d = 1 up to 0.231. P/P0 mean absolute error 0.0073.

    NOT strictly like for like: the authors reduced Table D.1's Gj/Gj_mean
    with their own Table 4.1 Cd(Gc/Gj) above each plate's break point
    (CR-3630 Sec. 4; the extraction record), the network uses CR-2681's
    generic law. The measured quantity is P/P0. ``cd_crossflow=
    'florschuetz_metzger_su_1982'`` is the authors' law.

    History: #802 (constant Cd) z/d >= 2 max 0.111 (B(5,4,2) mc/mj 1.01,
    row 10), z/d = 1 up to 0.269 (B(5,4,1) mc/mj 0.99, row 4: the stagnant
    rows alternating in sign, #799), P/P0 0.0075. Before #782 (kernels fed
    T0): median 0.021, max 0.100.
    """
    wide, narrow, p_err = [], [], []
    for key in _solved():
        rows = _table()[key]
        err = _gj_err(key)
        (narrow if key[2] == 1 else wide).append(err)
        _, _, p = _solve(key)
        p_err += [abs(a - r["P_P0"]) for a, r in zip(p, rows, strict=True)]
    print(
        f"\nGj/Gj_mean worst row: z/d>=2 median {statistics.median(wide):.3f} max {max(wide):.3f};"
        f" z/d=1 max {max(narrow):.3f}; P/P0 mean |err| {statistics.mean(p_err):.4f}"
    )
    assert len(wide) + len(narrow) == 42
    assert statistics.median(wide) < 0.03
    # From the measurement band, not from the score: CR-3630 Appendix B
    # (Section 4.4.1) puts Gj/Gj_mean at about +/-2%. The worst z/d >= 2 row
    # is B(5,4,2) mc/mj 1.01 row 10, measured Gj/Gj_mean 1.98, where the
    # network missed by 0.100 with the kernels fed T0 and a constant Cd
    # (the authors' Table 4.1 has CD falling with Gc/Gj; 0.085 under the
    # crossflow Cd, #800). A change smaller than that row's band,
    # 0.02 * 1.98 = 0.04, cannot be told from the data: the guard is
    # 0.100 + 0.04.
    assert max(wide) < 0.14
    assert max(narrow) < 0.30
    assert statistics.mean(p_err) < 0.010


def test_row_ones_crossflow_is_the_bypass() -> None:
    """Gc/Gj at row 1 is the initial crossflow alone: where it is large the
    table pins the hand-off directly. Rows with z/d >= 2 and measured Gc/Gj
    above 0.3; measured worst 2026-10-22: 11%, which includes the table's own
    mass closure (down to 0.91) between the metered bypass and the jets."""
    errs = []
    for key, rows in _table().items():
        if key[2] == 1 or key[3] == 0.0 or rows[0]["Gc_Gj"] < 0.3:
            continue
        _, gc, _ = _solve(key)
        errs.append(abs(gc[0] / rows[0]["Gc_Gj"] - 1.0))
    print(f"\nrow-1 Gc/Gj: n {len(errs)}, worst {max(errs):.1%}")
    assert len(errs) >= 8
    assert max(errs) < 0.15


def _jet_compressibility_gap(key: Key, pp0: list[float]) -> float:
    """Worst row's |Gj/Gj_mean| difference between CR-3630 Eq. 4.1 (the
    isentropic orifice, constant Cd) and its incompressible form
    sqrt(1 - P/P0), both applied to the same P/P0 and normalised by their
    mean: how far compressibility alone moves the split on this test."""
    xn, yn, z, mc = key
    T0, P0 = _table()[key][0]["T0_K"], _table()[key][0]["P0_kPa"] * 1.0e3
    g = cb.complete_state(T0, P0, X_AIR).thermo.gamma
    comp = [p ** (1.0 / g) * math.sqrt(max(1.0 - p ** ((g - 1.0) / g), 0.0)) for p in pp0]
    inc = [math.sqrt(max(1.0 - p, 0.0)) for p in pp0]
    cm, im = sum(comp) / len(comp), sum(inc) / len(inc)
    return max(abs(a / cm - b / im) for a, b in zip(comp, inc, strict=True))


def test_the_network_discretises_the_authors_model() -> None:
    """At xn/d 5, z/d >= 2 -- low pressure ratio, little wall shear -- the
    network is the authors' Eq. 4.10 discretised. At xn/d 10 (P0 up to 134
    kPa through 1.27 mm holes) the compressible jets make the supply
    measurably more uniform than the incompressible closed form; that is
    physics the closed form omits.

    MODEL VS MODEL, not fidelity: Eq. 4.10 is the authors' incompressible
    closed form (constant Cd, no friction, constant density in the jets and
    the channel); the network is compressible throughout. The bound is
    therefore the incompressible/compressible gap itself, per test: the
    network may sit no further from Eq. 4.10 on any row than Eq. 4.1 sits
    from its own incompressible form on the network's solved P/P0. That gap
    runs from 0.002 (B(5,8,3) at mc/mj 0) to 0.062 (B(5,4,2) at mc/mj 1.01). It is
    taken on the solved P/P0 because Table D.1 prints P/P0 to two decimals,
    which flattens the z/d = 3 tests to one value on all ten rows.
    Measured 2026-10-27 (#782): the network uses 59-84% of each test's gap,
    worst row 0.052 against a gap of 0.062 (B(5,4,2) at mc/mj 1.01); with
    the kernels fed T0 the worst row was 0.041."""
    worst, used = 0.0, 0.0
    for key in _solved():
        xn, yn, z, mc = key
        if xn != 5 or z == 1:
            continue
        gj, _, p = _solve(key)
        cf = _closed_form(xn, yn, z, mc)
        err = max(abs(a - b) for a, b in zip(gj, cf, strict=True))
        gap = _jet_compressibility_gap(key, p)
        worst, used = max(worst, err), max(used, err / gap)
        assert err < gap, (key, err, gap)
    print(f"\nnetwork vs Eq. 4.10, xn/d 5, z/d >= 2: worst row {worst:.3f}, max {used:.0%} of gap")


def test_dropping_the_bypass_scores_worse() -> None:
    """Falsify the metric: the same jets with no initial crossflow must miss
    the tests that had one. Every z/d >= 2 test at mc/mj near 1."""
    for key in _table():
        if key[2] == 1 or key[3] < 0.95:
            continue
        assert _gj_err(key, with_bypass=False) > 2.0 * _gj_err(key), key
