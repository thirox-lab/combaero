"""A film's undiluted gas takes its node's Tt_target heat (#1181).

``_film_unmix`` mixes a film discharge node's inflows, without the film's
own coolant, into the T_inf the film's effectiveness is defined against
(#694), with the node's energy boundaries: Q-mode and fraction-mode heat
and, since #1169, the wall heat. A Tt_target boundary (#792) entered with
its ``Q`` (0 in target mode), not the heat it applies (``heat``, its
``{id}.Q`` unknown at x): T_inf left out the solved heat (1494.1 K against
1357.8 K on this network). The owner's decision (2026-11-15): include it,
whatever mode set the heat. The mix now takes each boundary's ``heat``, set
from x before it runs, and its relay takes the unknown's column.

The network is #1169's (``test_film_inf_wall_heat_1169._net``), its plenum
"g" (the film's discharge node) given a Tt_target boundary.
"""

from __future__ import annotations

import warnings

import numpy as np
from _fd_jacobian import fd_jacobian, scaled_error
from test_film_inf_wall_heat_1169 import KEY, _cold, _net

from combaero.network import EnergyBoundary, FlowNetwork, NetworkSolver

TT_TARGET = 1350.0
#: The analytic Jacobian against FD, scaled (``scaled_error``): measured
#: <= 2.5e-7; without the heat unknown's column in T_inf's relay 4.5e-5
#: at the root and 6.9e-3 at the perturbed point.
JAC_TOL = 1e-6

#: The solved root (repr round-trips each float exactly).
X_ROOT = {
    "c0_in.P": 201066.59408178055,
    "c0_in.Pt": 201066.59408178055,
    "g.P": 219132.56099795087,
    "g.Pt": 219132.56099795087,
    "heat.Q": -79469.42597223246,
    "hp.P": 218961.37439959578,
    "hp.Pt": 218961.37439959578,
    "cp.P": 200654.7359772041,
    "cp.Pt": 200654.7359772041,
    "c0p.P": 200725.67399531679,
    "c0p.Pt": 200725.67399531679,
    "h0.m_dot": 0.48458359915527355,
    "c0.m_dot": 0.05,
    "c0o.m_dot": 0.026295235420508243,
    "film.m_dot": 0.006112143220199574,
    "hot.m_dot": 0.4906957423754732,
    "ho.m_dot": 0.4906957423754732,
    "cold.m_dot": 0.023704764579491756,
    "co.m_dot": 0.023704764579491756,
}


def _net_with(eb: EnergyBoundary) -> FlowNetwork:
    g = _net()
    g.nodes["g"].add_energy_boundary(eb)
    return g


def _fresh(eb: EnergyBoundary | None = None) -> NetworkSolver:
    s = NetworkSolver(_net_with(eb or EnergyBoundary("heat", Tt_target=TT_TARGET)))
    s.network.resolve_all_topology()
    s._build_x0()
    return s


def _x(s: NetworkSolver) -> np.ndarray:
    return np.array([X_ROOT[n] for n in s.unknown_names], dtype=float)


def _F(s: NetworkSolver, x: np.ndarray) -> np.ndarray:
    _cold(s)
    return np.asarray(s._residuals(x), dtype=float)


def test_the_undiluted_gas_takes_the_target_heat_as_a_fixed_q() -> None:
    """At the root, T_inf equals that of the same network with a fixed Q
    equal to the solved heat, evaluated at the same state: the target
    boundary's heat enters as any other heat. It is far below the 1494 K
    the mix gave without it."""
    s = _fresh()
    x = _x(s)
    _F(s, x)
    T_target = s._film_inf_states[KEY][0]
    assert T_target < 1400.0

    Q = X_ROOT["heat.Q"]
    q = _fresh(EnergyBoundary("heat", Q=Q))
    xq = np.array([X_ROOT[n] for n in q.unknown_names], dtype=float)
    assert len(xq) == len(x) - 1
    _F(q, xq)
    T_fixed = q._film_inf_states[KEY][0]
    assert abs(T_target - T_fixed) <= q._PROPAGATION_SETTLE_TOL * T_fixed


def test_one_pass_reads_the_target_heat_at_this_x() -> None:
    """A single pass (cap 1) on a fresh solver gives the settled T_inf: the
    heat unknown's value is set before the mix reads it. Set after, the
    first pass mixed the boundary's initial heat, 0 (1494.1 K); only a
    repeat pass picked x's up."""
    s = _fresh()
    s._MAX_PROPAGATION_PASSES = 1
    x = _x(s)
    with warnings.catch_warnings():
        # One pass does not settle the outlets' lagged reads: warned.
        warnings.simplefilter("ignore", RuntimeWarning)
        _F(s, x)
    T_one_pass = s._film_inf_states[KEY][0]
    full = _fresh()
    _F(full, x)
    T_settled = full._film_inf_states[KEY][0]
    assert T_one_pass == T_settled
    assert T_one_pass < 1400.0


def test_the_target_solve_and_the_fixed_q_solve_agree() -> None:
    """The two solves reach the same root: the fixed-Q network solved with
    the target's solved heat holds the node at Tt_target and gives the same
    T_inf, to the solve's tolerance."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = NetworkSolver(_net_with(EnergyBoundary("heat", Tt_target=TT_TARGET)))
        r = s.solve()
        assert r["__success__"]
        T_target = s._film_inf_states[KEY][0]
        q = NetworkSolver(_net_with(EnergyBoundary("heat", Q=float(r["heat.Q"]))))
        rq = q.solve()
        assert rq["__success__"]
    assert abs(q._derived_states["g"][0] - TT_TARGET) < 1e-6
    assert abs(q._film_inf_states[KEY][0] - T_target) < 1e-6


def _jac_err(s: NetworkSolver, x: np.ndarray) -> float:
    _cold(s)
    J = s._residuals_and_jacobian(x)[1]
    orig = s._residuals

    def F(z: np.ndarray) -> np.ndarray:
        _cold(s)
        return orig(z)

    s._residuals = F
    try:
        fd = fd_jacobian(s, x)
    finally:
        del s._residuals
    return float(scaled_error(J, fd, x).max())


def test_the_jacobian_reaches_the_target_heat_in_the_undiluted_gas() -> None:
    """FD over every row, at the root and a perturbed point: d(T_inf)/dQ
    reaches the film's T_aw and the rows the film reaches. Dropping the
    heat unknown's column from T_inf's relay (its value kept) fails it."""
    x = _x(_fresh())
    xp = x * (1.0 + 0.03 * np.sin(np.arange(len(x)) + 1.0))
    for point in (x, xp):
        assert _jac_err(_fresh(), point) < JAC_TOL

    s = _fresh()
    unmix = s._film_unmix
    i_q = s.unknown_names.index("heat.Q")

    def no_target_column(nid, pid, stream_info, walls, x, relay, up, with_relay):  # noqa: ANN001, ANN202
        done = unmix(nid, pid, stream_info, walls, x, relay, up, with_relay)
        if done and with_relay:
            relay[f"__film_inf__{pid}"].pop(i_q, None)
        return done

    s._film_unmix = no_target_column
    for point in (x, xp):
        assert _jac_err(s, point) > 10.0 * JAC_TOL
