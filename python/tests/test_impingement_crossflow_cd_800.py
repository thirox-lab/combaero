"""The jet holes' Cd under the gap's crossflow (#800).

A crossflow passing the jets' exits lowers their Cd, steeply at large Gc/Gj:
Rogers & Hersh (1976, NASA CR-2681) on single holes, and Florschuetz,
Metzger & Su (1982, NASA CR-3630 Table 4.1) on the impingement plates
themselves. ``HoleCd(exit_crossflow=...)`` selects the law; the default is CR-2681's.

The stagnant upstream rows of CR-3630 B(5,4,1) at mc/mj 0.99 (Gc/Gj ~125)
used to come out with alternating-sign jet flows (#799): a constant Cd makes
the row flow ~ sqrt(dP), whose infinite slope at dP = 0 lets the centred
station scheme carry a period-2 mode. Under the crossflow law the flow goes
as dP^1.07 there and the rows are monotone, as the data are.
"""

from __future__ import annotations

import csv
import math

import numpy as np
import pytest
from scipy.optimize._numdiff import approx_derivative

import combaero as cb
from combaero.network import (
    FlowNetwork,
    HoleCd,
    ImpingementPlateElement,
    MassFlowBoundary,
    NetworkSolver,
    PressureBoundary,
)
from combaero.network.impingement_array import ImpingementArray
from validation.cooling.schema import DATA_ROOT

Y_AIR = cb.mole_to_mass(cb.species.dry_air())
X_AIR = cb.species.dry_air()

# CR-3630 B(5,4,1): 18 holes of 0.254 cm per row, Table 3.1 Cd 0.85; each
# test's own T0, P0 and Re_j from Table D.1.
D_JET, SPAN, N_ROWS = 0.00254, 0.1829, 10
TABLE = DATA_ROOT / "florschuetz_metzger_su1982" / "table_d1_flow.csv"


def _test_row(mc: float) -> dict[str, float]:
    with TABLE.open() as fh:
        for r in csv.DictReader(fh):
            key = (int(r["xn_d"]), int(r["yn_d"]), int(r["z_d"]), float(r["mc_mj"]))
            if key == (5, 4, 1, mc):
                return {k: float(v) for k, v in r.items() if v}
    raise KeyError(mc)


def _b541(
    mc: float, law: str, x0: np.ndarray | None = None
) -> tuple[NetworkSolver, dict, ImpingementArray]:
    row = _test_row(mc)
    T0, P0 = row["T0_K"], row["P0_kPa"] * 1.0e3
    mu = cb.complete_state(T0, P0, X_AIR).transport.mu
    m_j = row["Re_j"] * mu * N_ROWS * 18 * math.pi * D_JET / 4.0
    net = FlowNetwork()
    s = PressureBoundary("supply")
    s.Pt, s.Tt, s.Y = P0, T0, Y_AIR
    net.add_node(s)
    net.add_node(MassFlowBoundary("exit", m_dot=m_j * (1.0 + mc), Tt=T0, Y=Y_AIR))
    net.add_node(MassFlowBoundary("feed", m_dot=mc * m_j, Tt=T0, Y=Y_AIR))
    arr = ImpingementArray(
        "a", n_rows=N_ROWS, d_jet=D_JET, xn_d=5, yn_d=4, z_d=1, span=SPAN,
        plate_thickness=D_JET, hole_cd=HoleCd(Cd=0.85, exit_crossflow=law),
    )  # fmt: skip
    arr.add_to(net, "supply", "exit", bypass="feed")
    solver = NetworkSolver(net)
    res = solver.solve() if x0 is None else solver.solve(x0=x0)
    assert res["__success__"], res.get("__message__")
    return solver, res, arr


def _solved(mc: float, law: str) -> tuple[NetworkSolver, dict, ImpingementArray]:
    """Every case cold. Until #1021 Table 4.1 at mc/mj 0.99 started from the
    default law's root here: cold, it failed the primary solve (|F| 4.5e-1)
    and the outlet-referenced retry (1.7e-2) before the constant-Cd retry
    reached the root, 16.9-25 s in all. The solver now runs that stage first
    for Table 4.1 plates."""
    return _b541(mc, law)


def test_table_4_1_solves_cold_from_the_constant_cd_stage() -> None:
    """#1021: CR-3630 B(5,4,1) at mc/mj 0.99 under Table 4.1. The cold x0
    seeds the stagnant rows backwards (p1-p3 at -60x, -21x, -7x the root) in
    a basin too narrow for the primary attempt; the constant-Cd root is a
    start inside it. The solver runs that stage FIRST and lands on the root
    a warm start from the default law's root finds (measured 7e-8 relative,
    5.5 s against 25.5 s through the retries on main 4843598f)."""
    _, cold, _ = _b541(0.99, "florschuetz_metzger_su_1982")
    assert "run first for the Table 4.1 Cd" in cold["__message__"]
    _, rh, _ = _b541(0.99, "rogers_hersh_1976")
    _, warm, _ = _b541(
        0.99, "florschuetz_metzger_su_1982", np.array(rh["__x_solution__"], dtype=float)
    )
    x, w = np.array(cold["__x_solution__"]), np.array(warm["__x_solution__"])
    assert np.max(np.abs(x - w) / np.maximum(np.abs(w), 1e-12)) < 1e-6
    # The default law keeps its cold path: no constant-Cd stage.
    assert "constant-Cd" not in rh["__message__"]


def _jets(res: dict, arr: ImpingementArray) -> list[float]:
    return [float(res[f"{arr.plate_id(i)}.m_dot"]) for i in range(1, N_ROWS + 1)]


@pytest.mark.parametrize("law", ["rogers_hersh_1976", "florschuetz_metzger_su_1982"])
def test_the_stagnant_rows_are_monotone(law: str) -> None:
    """#799: the constant Cd gives rows 2 and 4 backwards; the data
    (Gj/Gj_mean 0.00, 0.02, 0.07, 0.16, 0.35, ...) rise monotonically."""
    _, res, arr = _solved(0.99, law)
    m = _jets(res, arr)
    mean = sum(m) / N_ROWS
    assert all(b > a for a, b in zip(m[1:], m[2:], strict=False))
    assert min(m) / mean > -0.01  # row 1 at -0.0005 of the mean, not -0.07


def test_the_constant_cd_still_alternates() -> None:
    """The falsification: with ``'none'`` the #799 mode is back."""
    _, res, arr = _b541(0.99, "none")
    m = _jets(res, arr)
    assert any(x < 0.0 for x in m[1:5])


@pytest.mark.parametrize("law", ["rogers_hersh_1976", "florschuetz_metzger_su_1982"])
def test_the_jacobian_matches_finite_differences(law: str) -> None:
    """Every row's Cd reads its own drive and the approaching crossflow; at
    the root and nearby, the analytic Jacobian is the FD one. Step 1e-7: the
    1e-6 central difference is truncation (it falls 100x for a 10x smaller
    step). With a constant Cd the same check reads 0.23 at the root: the
    sqrt(dP) curvature of the near-stagnant rows that the law removes."""
    solver, res, _ = _b541(0.51, law)
    x0 = np.array(res["__x_solution__"], dtype=float)
    for x in (x0, x0 * (1.0 + 0.002 * np.sin(np.arange(len(x0))))):
        J = solver._residuals_and_jacobian(x)[1].toarray()
        fd = approx_derivative(
            solver._residuals, x, method="3-point", abs_step=np.abs(x) * 1e-7 + 1e-12
        )
        cols = np.maximum(np.abs(x), 1e-12)[None, :]
        scale = np.maximum(np.max(np.abs(fd * cols), axis=1, keepdims=True), 1e-300)
        assert float(np.max(np.abs((J - fd) * cols) / scale)) < 1e-5


def test_a_weak_crossflow_keeps_the_holes_own_cd() -> None:
    """Rogers & Hersh is capped at the hole's own Cd below Gc/Gj ~1.5 (at
    0.85): rows there report exactly that Cd."""
    _, res, arr = _b541(0.0, "rogers_hersh_1976")
    diag = res["__element_diag__"]
    held = 0
    for i in range(1, N_ROWS + 1):
        d = diag[arr.plate_id(i)]
        if d["Gc_Gj"] < 1.0:
            assert d["Cd"] == d["Cd_own"] == 0.85
            held += 1
    assert held >= 5


def test_the_reported_cd_is_the_law_at_the_solved_state() -> None:
    """Diagnostics report the Cd the residual used: CR-2681's law at the
    row's ideal and crossflow mass velocities."""
    solver, res, arr = _b541(0.99, "rogers_hersh_1976")
    x = np.array(res["__x_solution__"], dtype=float)
    solver._residuals(x)
    row = 3
    pid = arr.plate_id(row)
    el = solver.network.elements[pid]
    a = solver._get_node_state(solver.network.nodes["supply"], x)
    b = solver._get_node_state(solver.network.nodes[arr.node_id(row)], x)
    rho = cb.density(a.Tt, a.P, a.X)
    G_i = math.sqrt(2.0 * rho * (a.Pt - b.P))
    m_c = float(res[f"{arr.segment_id(row - 1)}.m_dot"])
    expect = cb.rogers_hersh_1976_cd(0.85, G_i, m_c / el.channel_area, False).Cd
    assert res["__element_diag__"][pid]["Cd"] == pytest.approx(expect, rel=1e-9)
    assert expect < 0.2


def test_table_4_1_is_refused_for_a_plate_it_does_not_cover() -> None:
    with pytest.raises(ValueError, match="Table 4.1"):
        ImpingementPlateElement(
            "p", "s", "c", d_jet=D_JET, xn_d=5, yn_d=8, z_d=2, span=SPAN,
            plate_thickness=D_JET,
            hole_cd=HoleCd(exit_crossflow="florschuetz_metzger_su_1982"),
        )  # fmt: skip


def test_an_unknown_law_is_refused() -> None:
    with pytest.raises(ValueError, match="exit_crossflow 'mcgreehan' is not one of"):
        HoleCd(exit_crossflow="mcgreehan")
    with pytest.raises(ValueError, match="does not apply to this element"):
        ImpingementPlateElement(
            "p", "s", "c", d_jet=D_JET, xn_d=5, yn_d=4, z_d=1, span=SPAN,
            plate_thickness=D_JET, hole_cd=HoleCd(exit_crossflow="rowbury_2001"),
        )  # fmt: skip


def test_the_envelope_is_reported_where_the_law_acts() -> None:
    """Rogers & Hersh's envelope (CR-2681: Vi/V_inf 0.004-0.5, branch-aware
    L/D, V_inf 24-85 m/s), reported per variable as information, never
    a warning, and only on rows where the law acts. At mc/mj 0.99 the
    first row's jet barely moves (Vi/V_inf below 0.004) and the rows the
    cap holds (the hole's own Cd) carry no flag."""
    _, res, arr = _solved(0.99, "rogers_hersh_1976")
    diag = res["__element_diag__"]
    rows = [diag[arr.plate_id(i)] for i in range(1, N_ROWS + 1)]
    assert rows[0]["Cd_crossflow_extrapolated_velocity_ratio"] == 1.0
    for p in rows:
        assert p["Cd_crossflow_extrapolated_l_d"] == 0.0  # L/D 1, Fig. 11's own
        if p["Cd"] == p["Cd_own"]:  # capped: the law does not act
            assert p["Cd_crossflow_extrapolated_velocity_ratio"] == 0.0
            assert p["Cd_crossflow_extrapolated_grazing_velocity"] == 0.0
    summary = arr.summarize(diag)
    assert summary["Cd_crossflow_extrapolated_velocity_ratio"] == 1.0
    _, res, arr = _solved(0.99, "none")
    for i in range(1, N_ROWS + 1):
        p = res["__element_diag__"][arr.plate_id(i)]
        assert p["Cd_crossflow_extrapolated_velocity_ratio"] == 0.0
