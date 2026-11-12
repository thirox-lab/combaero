"""The wall relay's pressure column (#496).

A wall's heat moves with the channel's static pressure at fixed flow through
the adiabatic wall temperature only: T_aw = T + r v^2/(2 cp), v = mdot/(rho A)
and rho ~ P, so dT_aw/dP = -2 (T_aw - T)/P; h does not move (Re, Pr, k do
not). The relay had no pressure column: 1e-3 of the scaled Jacobian in a
two-wall network, its largest entry once #463 had landed.
"""

from __future__ import annotations

import numpy as np
import pytest
import test_solver_robustness_481 as tsr
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import NetworkSolver

AIR = np.array(cb.species.dry_air())


@pytest.mark.parametrize("corr", ["gnielinski", "petukhov"])
def test_the_adiabatic_wall_temperature_moves_with_pressure_at_fixed_flow(corr: str) -> None:
    T, D, L, mdot = 500.0, 0.02, 0.5, 0.05
    A = np.pi * D * D / 4

    def at(P: float) -> cb.ChannelResult:
        v = mdot / (cb.density(T, P, list(AIR)) * A)
        return cb.channel_smooth(T, P, AIR, v, D, L, 300.0, corr, True, 1.0, 0.0, 1.0, 1.0)

    P, h = 1.5e5, 10.0
    r = at(P)
    assert r.dT_aw_dP == pytest.approx((at(P + h).T_aw - at(P - h).T_aw) / (2 * h), rel=1e-6)
    assert r.dT_aw_dP == pytest.approx(-2.0 * (r.T_aw - T) / P, rel=1e-12)
    assert at(P + h).h == pytest.approx(at(P - h).h, rel=1e-12)  # h does not move


@pytest.mark.parametrize("Dh", [None, 0.012], ids=["circular", "flat"])
def test_the_wall_coupled_network_jacobian_is_exact(Dh: float | None) -> None:
    g = tsr._series_walls()
    if Dh is not None:
        for eid in ("e1", "e2", "cc"):
            g.elements[eid].Dh = Dh
    s = NetworkSolver(g)
    x = np.array(s.solve()["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    # RESIDUAL-ROW-SCALING: compare J_ij |x_j| (see _build_residual_scales).
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    assert float(np.max(np.abs(Js - fds) / scale)) < 1e-7  # was 1e-3
