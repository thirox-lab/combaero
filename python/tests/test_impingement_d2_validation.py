"""An impingement array with an initial crossflow, against the heat transfer
of Florschuetz, Metzger and Su (1982), NASA CR-3630, Table D.2 (#467).

CROSS-SOURCE: accuracy, not fidelity (docs/VALIDATION_POLICY.md). Same lab,
different studies: the bypass reference temperature's eta_r closure is
fitted to CR-3936 Appendix A, the jet h to Florschuetz et al. (1981); neither
saw these points. The extraction record is
validation/cooling/extractions/florschuetz_metzger_su_1982_cr3630_table_d2.md
(signed off 2026-10-24).

The rig is solved as built, at each D.2 test's own Re_j and mc/mj (the
network of test_impingement_bypass_validation.py), with the initial
crossflow 20 K above the jets. Per row the plate reports Nu and T_aw, so

    eta = (T_aw - T_j) / (T_0 - T_j)

is D.2's ETA (Eq. 5.3: referenced to the entrance crossflow temperature).
The jet plenum's P0 and T0 come from the D.1 test of the same geometry and
nearest mc/mj; D.2 prints neither. B(10,4,2) has no D.1 test, so it is
left out rather than run at an invented pressure: 10 inline geometries, 40
tests, 400 in-array rows (300 with an initial crossflow).

Staggered B(5,4,3)S is outside the eta_r closure and is not scored here.
"""

from __future__ import annotations

import csv
import math
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

DATA = DATA_ROOT / "florschuetz_metzger_su1982"
N_ROWS = 10
SPAN = 0.1829  # channel width [m]: 18 holes at yn/d 4, d 0.254 cm
D_JET = {5: 0.00254, 10: 0.00127}  # Table 3.1
PLATE_CD = {(5, 4): 0.85, (5, 8): 0.80, (10, 4): 0.76, (10, 8): 0.76}  # Table 3.1
Y_AIR = cb.mole_to_mass(cb.species.dry_air())
X_AIR = cb.species.dry_air()
# The bypass above the jets. The score does not depend on it: 5 K gives
# the same eta RMSE (0.078) as 20 K.
BYPASS_DT = 20.0  # [K]

# About 20 s on one worker: runs in CI's parallel slow job.
pytestmark = pytest.mark.slow

Geom = tuple[int, int, int]
Test = tuple[int, int, int, float, float]  # geometry, mc/mj, Re_j / 1000

# Tests the network does not solve from a cold start, each with its reason
# (a strict xfail below). None since #802: B(5,4,1) at mc/mj 1.01 stalled on
# the false Fanno choke of its exit segment (a station's static P fed to the
# kernel as P0: 99.7% of that choke, 81% of the true one). It solves cold
# again, its rows 1-5 starved (jets reversed or Gc/Gj past 4.44), as before
# #782.
STALLED: dict[Test, str] = {}


@cache
def _d1() -> dict[tuple[int, int, int, float], list[dict[str, str]]]:
    out: dict[tuple[int, int, int, float], list[dict[str, str]]] = defaultdict(list)
    with (DATA / "table_d1_flow.csv").open() as f:
        for r in csv.DictReader(f):
            out[(int(r["xn_d"]), int(r["yn_d"]), int(r["z_d"]), float(r["mc_mj"]))].append(r)
    return dict(out)


@cache
def _d2() -> dict[Test, list[dict[str, float]]]:
    """Inline tests, the ten in-array rows each."""
    out: dict[Test, list[dict[str, float]]] = defaultdict(list)
    with (DATA / "table_d2_heat_transfer.csv").open() as f:
        for r in csv.DictReader(f):
            if r["pattern"] != "I" or not 0.0 < float(r["x_L"]) < 1.0:
                continue
            key = (
                int(r["xn_d"]),
                int(r["yn_d"]),
                int(r["z_d"]),
                float(r["mc_mj"]),
                float(r["Re_mean_k"]),
            )
            out[key].append(
                {"Nu": float(r["Nu"]), "eta": float(r["eta"]) if r["eta"] else math.nan}
            )
    return dict(out)


def _flow_test(key: Test) -> list[dict[str, str]] | None:
    """The D.1 test of the same geometry with the nearest mc/mj."""
    cands = [k for k in _d1() if k[:3] == key[:3]]
    if not cands:
        return None
    return _d1()[min(cands, key=lambda k: abs(k[3] - key[3]))]


@cache
def _solve(key: Test) -> list[dict[str, float]]:
    """Per row the plate diagnostics, plus ``eta`` from T_aw."""
    xn, yn, z, mc, re_k = key
    flow = _flow_test(key)
    assert flow is not None
    T0, P0 = float(flow[0]["T0_K"]), float(flow[0]["P0_kPa"]) * 1.0e3
    d = D_JET[xn]
    mu = cb.complete_state(T0, P0, X_AIR).transport.mu
    n_holes = round(SPAN / (yn * d))
    m_j = re_k * 1.0e3 * mu * N_ROWS * n_holes * math.pi * d / 4.0
    m_c = mc * m_j

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
        net.add_node(MassFlowBoundary("feed", m_dot=m_c, Tt=T0 + BYPASS_DT, Y=Y_AIR))
        bypass = "feed"
    arr.add_to(net, "supply", "exit", bypass=bypass)
    res = NetworkSolver(net).solve()
    assert res["__success__"], f"{key}: {res.get('__message__')}"
    diag = res["__element_diag__"]
    unresolved = arr.summarize(diag)["rows_unresolved"]
    rows = []
    for i in range(1, N_ROWS + 1):
        p = dict(diag[arr.plate_id(i)])
        p["unresolved"] = unresolved[i - 1]
        p["eta"] = (p["T_aw"] - T0) / BYPASS_DT if m_c > 0.0 else math.nan
        rows.append(p)
    return rows


def _scored() -> list[Test]:
    return [k for k in sorted(_d2()) if _flow_test(k) is not None]


def _solving() -> list[Test]:
    return [k for k in _scored() if k not in STALLED]


def _bias_rms(v: list[float]) -> tuple[float, float]:
    return sum(v) / len(v), math.sqrt(sum(x * x for x in v) / len(v))


def test_the_scored_set() -> None:
    """What the score rests on: B(10,4,2) is the one geometry dropped."""
    scored = _scored()
    assert {k[:3] for k in _d2()} - {k[:3] for k in scored} == {(10, 4, 2)}
    assert len(scored) == 40
    assert sum(1 for k in scored if k[3] > 0.0) == 30


@pytest.mark.parametrize(
    "key",
    [pytest.param(k, marks=pytest.mark.xfail(reason=r, strict=True)) for k, r in STALLED.items()],
    ids=str,
)
def test_the_stalled_tests_still_stall(key: Test) -> None:
    """Strict: a STALLED test that starts solving fails here, so it is looked
    at and returned to the score rather than left out silently."""
    _solve(key)


def test_the_bypass_reference_temperature() -> None:
    """D.2's ETA from the network's T_aw = T_j + eta_r phi_0 (T_0 - T_j).

    Measured 2026-10-27 (#802) over all 300 rows: bias +0.013, RMSE 0.078.
    Without B(5,4,1) mc/mj 1.01, which #782 to #802 left out (STALLED), 290
    rows: bias +0.010, RMSE 0.075, the same before and after #802.
    Measured 2026-10-24 over 300 rows: bias +0.013, RMSE 0.078; by row the
    bias runs from -0.09 at row 1 to +0.06 at row 10. Contrasts on the same
    rows, from the same solves: the mixing share alone (eta_r = 1, T_aw the
    crossflow's mixed-mean temperature) RMSE 0.21, nothing (eta = 0) 0.52.
    Before #728 a reversed jet row dropped to eta_r = 0: RMSE 0.163.
    """
    err, mix, none = [], [], []
    for key in _solving():
        if key[3] == 0.0:
            continue
        for p, r in zip(_solve(key), _d2()[key], strict=True):
            err.append(p["eta"] - r["eta"])
            mix.append(p["bypass_fraction"] - r["eta"])
            none.append(-r["eta"])
    bias, rms = _bias_rms(err)
    print(
        f"\nD.2 eta: n {len(err)} bias {bias:+.3f} RMSE {rms:.3f};"
        f" eta_r=1 RMSE {_bias_rms(mix)[1]:.3f}; eta=0 RMSE {_bias_rms(none)[1]:.3f}"
    )
    assert len(err) == 300
    assert rms < 0.085
    assert abs(bias) < 0.025
    # The closure must earn its place over the mixing share alone.
    assert rms < 0.5 * _bias_rms(mix)[1]


def test_the_bypass_closure_at_the_measured_flows() -> None:
    """eta_r phi_0 at D.1's measured Gc/Gj and jet split, no network: the
    closure alone, without the flow model's error. Rests on matching each
    D.2 test to a D.1 run (mc/mj within 0.05; the report does not say they
    are the same runs): 28 tests, 279 rows.

    Measured 2026-10-24: bias +0.024, RMSE 0.077 -- the same as end to end.
    """
    err = []
    for key in _scored():
        if key[3] == 0.0:
            continue
        cands = [k for k in _d1() if k[:3] == key[:3] and k[3] > 0.0]
        k1 = min(cands, key=lambda k: abs(k[3] - key[3]))
        if abs(k1[3] - key[3]) > 0.05:
            continue
        flow = _d1()[k1]
        closure = float(flow[0]["mass_closure"])
        g = [float(r["Gj_Gjbar"]) for r in flow]
        share = k1[3] / closure
        for i, r in enumerate(_d2()[key]):
            if not flow[i]["Gc_Gj"]:
                continue
            phi_0 = share / (share + sum(g[:i]) / N_ROWS)
            eta_r = cb.impingement_eta_r(float(flow[i]["Gc_Gj"]), *key[:3]).eta
            err.append(eta_r * phi_0 - r["eta"])
    bias, rms = _bias_rms(err)
    print(f"\nD.2 eta at D.1 flows: n {len(err)} bias {bias:+.3f} RMSE {rms:.3f}")
    assert len(err) == 279
    assert rms < 0.085
    assert abs(bias) < 0.030


def test_the_jet_heat_transfer() -> None:
    """Florschuetz (1981) Nu per row, at the network's own Re_j and Gc/Gj.

    Measured 2026-10-28 (#800, the holes' crossflow Cd, Rogers & Hersh 1976
    by default; all 40 tests, Nu_model / Nu_D2 - 1):

    - no initial crossflow, 100 rows: bias +5.8%, RMS 8.4%;
    - with it, inside the correlation's range, 197 rows: -2.7%, 21.5%;
    - with it, flagged outside the range (surface_extrapolated), 103 rows:
      -17.9%, 35.4%. Among them the 11 starved rows (network jet reversed or
      Gc/Gj past 4.44), at -81%: the surface is cooled by the crossflow and
      the jet correlation has no channel-flow floor (CR-3936 conclusion 6;
      the open #467 item). Pinned so a floor shows here. The crossflow Cd
      weakens the starved jets, so more rows pass Gc/Gj 4.44 (8 -> 11).

    #802 to #800 (constant Cd): inside 196, outside 104 at -16.9% / 35.1%,
    starved 8 at -95%.

    #782 to #802 (that test STALLED): outside 94 rows at -13.8% / 29.0%,
    starved 3; the other groups as at #802. Before #782 (2026-10-24): inside
    195 rows, outside 105 at -16.9% / 34.9%, starved 8.

    The array summary marks only starved rows unresolved (#729): 2 of the 11
    under the crossflow Cd (all 8 under a constant Cd).
    """
    groups: dict[str, list[float]] = defaultdict(list)
    n_unresolved = 0
    for key in _solving():
        for p, r in zip(_solve(key), _d2()[key], strict=True):
            e = p["Nu"] / r["Nu"] - 1.0
            starved = p["Gc_Gj"] < 0.0 or p["Gc_Gj"] > 4.44
            # Only a starved row can be unresolved (#729). Under the crossflow
            # Cd (#800) most starved rows' small flows take a drop above the
            # chain's resolution: counted below.
            assert not p["unresolved"] or starved, (key, p["Gc_Gj"])
            n_unresolved += bool(p["unresolved"])
            if key[3] == 0.0:
                groups["none"].append(e)
            elif p["surface_extrapolated"]:
                groups["outside"].append(e)
                if starved:
                    groups["starved"].append(e)
            else:
                groups["inside"].append(e)
    stats = {k: _bias_rms(v) for k, v in groups.items()}
    print(f"\nD.2 unresolved rows: {n_unresolved}")
    print(
        "\nD.2 Nu: "
        + "; ".join(
            f"{k} n {len(groups[k])} bias {b:+.1%} RMS {s:.1%}" for k, (b, s) in stats.items()
        )
    )
    assert {k: len(v) for k, v in groups.items()} == {
        "none": 100,
        "inside": 197,
        "outside": 103,
        "starved": 11,
    }
    assert n_unresolved == 2
    assert abs(stats["none"][0]) < 0.07 and stats["none"][1] < 0.09
    assert abs(stats["inside"][0]) < 0.04 and stats["inside"][1] < 0.23
    assert stats["outside"][0] < -0.10 and stats["outside"][1] < 0.37
    assert stats["starved"][0] < -0.75
