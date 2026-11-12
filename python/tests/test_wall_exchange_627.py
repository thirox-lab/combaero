"""A ThermalWall passes the two streams' exchange, epsilon-NTU (#627).

The linear U A (T_aw,a - T_aw,b) at the inlets was its small-NTU limit: it
moved more heat than the streams hold at NTU of order 1, and kept a finite Q
as one stream stopped, so the stress hot-cold networks wandered to a stalled
hot side at -2e4 K. The arrangement (parallel or counter) follows from the
flow signs and the wall's declared orientation.
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian
from test_second_law import Y_AIR, _hot_cold_wall

import combaero as cb
from combaero.network import NetworkSolver, ThermalWall

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
from networks import hot_cold  # noqa: E402

X_AIR = list(cb.mass_to_mole(Y_AIR))


def _solve(net, **kw) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return s, s.solve(**kw)


def _jac_err(s: NetworkSolver, x: np.ndarray) -> float:
    """RESIDUAL-ROW-SCALING: compare J_ij |x_j|, relative to each row's largest."""
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()
    fd = fd_jacobian(s, x)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    a, b = J * cols, fd * cols
    scale = np.maximum(np.abs(b).max(axis=1, keepdims=True), 1e-300)
    return float((np.abs(a - b) / scale).max())


def _reorient(g, orientation: str):
    w = g.walls["w"]
    g.walls["w"] = ThermalWall(
        id="w",
        element_a=w.element_a,
        element_b=w.element_b,
        layers=w.layers,
        orientation=orientation,
    )
    return g


@pytest.mark.parametrize("seed", [11, 12, 16])
def test_the_stress_hot_cold_networks_that_stalled_now_solve_with_bounded_outlets(
    seed: int,
) -> None:
    """Seeds 11, 12, 16 of the hot-cold family: the linear wall pulled 1.7-28 kW
    out of a hot side LM had all but stopped, which went below 0 K. Now each
    converges, and every outlet lies between the two inlets."""
    g, _ = hot_cold(seed)
    _, r = _solve(g)
    assert r["__success__"], r["__message__"]
    T_hi, T_ci = g.nodes["hi"].Tt, g.nodes["ci"].Tt
    assert T_ci < r["hm.Tt"] < T_hi and T_ci < r["cm.Tt"] < T_hi
    assert r["__second_law_walls__"]["w"]["verdict"] == "pass"


def test_the_heat_is_at_most_what_the_weaker_stream_can_take() -> None:
    """Q <= Q_max = min over the streams of m (h(T_other,in) - h(T_own,in)),
    the variable-property bound: the capacity rates are secant-cp ones over
    the inlet span (cp_secant_mass_and_jacobian)."""
    g = _hot_cold_wall(0.006, 0.2, 5.0, 1600.0, 400.0, 0.7, (0.0005, 200.0))
    _, r = _solve(g)
    assert r["__success__"]
    q_max = min(
        abs(r[f"{s}ch.m_dot"]) * abs(cb.h_mass(1600.0, X_AIR) - cb.h_mass(400.0, X_AIR))
        for s in ("h", "c")
    )
    assert 0.0 < r["w.Q"] <= q_max


def test_an_opposed_wall_runs_counterflow_and_moves_more_heat() -> None:
    """The same channels, declared to lie opposite ways along the wall: at
    NTU ~ 1 counterflow moves more heat than parallel, as epsilon says."""
    par = _solve(_hot_cold_wall(0.006, 0.2, 5.0, 1600.0, 400.0, 0.7, (0.0005, 200.0)))[1]
    g = _reorient(_hot_cold_wall(0.006, 0.2, 5.0, 1600.0, 400.0, 0.7, (0.0005, 200.0)), "opposed")
    cnt = _solve(g)[1]
    assert par["__success__"] and cnt["__success__"]
    assert par["w.counterflow"] == 0.0 and cnt["w.counterflow"] == 1.0
    assert cnt["w.Q"] > 1.02 * par["w.Q"]
    assert cnt["hm.Tt"] > 400.0 and cnt["cm.Tt"] < 1600.0


def test_a_reversed_flow_is_the_same_wall_as_an_opposed_declaration() -> None:
    """The arrangement is geometry plus flow signs: a cold channel declared
    backwards (its flow negative) on an aligned wall is the counterflow wall
    of a forward channel on an opposed one -- same heat."""
    g_opp = _reorient(
        _hot_cold_wall(0.006, 0.2, 5.0, 1600.0, 400.0, 0.7, (0.0005, 200.0)), "opposed"
    )
    g_rev = _hot_cold_wall(
        0.006, 0.2, 5.0, 1600.0, 400.0, 0.7, (0.0005, 200.0), cold_backwards=True
    )
    r_opp, r_rev = _solve(g_opp)[1], _solve(g_rev)[1]
    assert r_opp["__success__"] and r_rev["__success__"]
    assert r_rev["cch.m_dot"] < 0.0 < r_opp["cch.m_dot"]
    assert r_rev["w.counterflow"] == 1.0
    assert r_rev["w.Q"] == pytest.approx(r_opp["w.Q"], rel=1e-6)


@pytest.mark.parametrize("orientation", ["aligned", "opposed"])
def test_the_wall_jacobian_is_complete(orientation: str) -> None:
    """The capacity rates' flow and temperature columns (dQ/dC dC/dm, the
    secant cp's dependence on BOTH inlets) are in the analytic Jacobian."""
    g = _reorient(_hot_cold_wall(0.006, 0.2, 5.0, 1600.0, 400.0, 0.7, (0.0005, 200.0)), orientation)
    s, r = _solve(g)
    assert r["__success__"]
    assert _jac_err(s, np.array(r["__x_solution__"])) < 1e-5


def test_the_secant_cp_is_the_enthalpy_slope_and_its_partials_are_exact() -> None:
    T1, T2 = 1600.0, 400.0
    cp, d1, d2 = cb.cp_secant_mass_and_jacobian(T1, T2, X_AIR)
    assert cp == pytest.approx((cb.h_mass(T1, X_AIR) - cb.h_mass(T2, X_AIR)) / (T1 - T2), rel=1e-14)
    h = 1e-3
    for T, d, at in ((T1, d1, lambda v: (v, T2)), (T2, d2, lambda v: (T1, v))):
        fd = (
            cb.cp_secant_mass_and_jacobian(*at(T + h), X_AIR)[0]
            - cb.cp_secant_mass_and_jacobian(*at(T - h), X_AIR)[0]
        ) / (2 * h)
        assert d == pytest.approx(fd, rel=1e-5)
    # the limit T2 -> T1 is cp, continuous across the switch
    near = cb.cp_secant_mass_and_jacobian(800.0, 800.0 * (1 + 2e-6), X_AIR)[0]
    at = cb.cp_secant_mass_and_jacobian(800.0, 800.0 * (1 + 5e-7), X_AIR)[0]
    assert near == pytest.approx(cb.cp_mass(800.0, X_AIR), rel=1e-5)
    assert at == pytest.approx(near, rel=1e-6)


def test_an_unknown_orientation_is_refused() -> None:
    with pytest.raises(ValueError, match="orientation"):
        ThermalWall(id="w", element_a="a", element_b="b", orientation="counter")
    assert math.isfinite(1.0)
