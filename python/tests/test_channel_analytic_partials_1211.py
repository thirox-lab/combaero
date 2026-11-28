"""channel_residuals_and_jacobian: the static-pressure, Tt and composition
partials are the chain rule through rho and mu, not differences (#1211).

The references are Richardson-extrapolated central differences of the drop
itself, resolved to well below the tolerances here. Until #1211 the kernel
differenced the drop internally: forward in Y (up to 6e-6 off), central in P
and Tt.
"""

from __future__ import annotations

import pytest

import combaero as cb
from combaero import _solver_tools

L, D = 1.0, 0.01
P_UP, T_UP = 1e5, 300.0


def _air_with_fuel() -> list[float]:
    """Air with a little CH4 and H2O, so heavy and light species both move mu."""
    Y = [0.0] * cb.num_species()
    for name, y in (("N2", 0.74), ("O2", 0.22), ("H2O", 0.02), ("CH4", 0.015), ("AR", 0.005)):
        Y[cb.species_index_from_name(name)] = y
    return Y


def _drop(m: float, P: float, T: float, Y: list[float], model: str, head: bool) -> float:
    r = _solver_tools.channel_residuals_and_jacobian(
        m, P, P, T, Y, P - 1.0, L, D, 1e-5, model, 1.0, head
    )
    return r.dP_calc


def _richardson(g, h: float) -> float:
    d1 = (g(h) - g(-h)) / (2.0 * h)
    d2 = (g(0.5 * h) - g(-0.5 * h)) / h
    return (4.0 * d2 - d1) / 3.0


def _one_sided(g, h: float) -> float:
    """Second-order forward difference, Richardson-extrapolated: for a species
    at Y = 0, whose central step would go negative (normalize_fractions clamps
    it back to 0, halving a central difference)."""

    def fwd(t: float) -> float:
        return (-3.0 * g(0.0) + 4.0 * g(t) - g(2.0 * t)) / (2.0 * t)

    return (4.0 * fwd(0.5 * h) - fwd(h)) / 3.0


def _y_step(Y: list[float], k: int, t: float) -> list[float]:
    """Y + t e_k renormalised: the step mole_gradient_to_mass's gradient is along."""
    Ys = list(Y)
    Ys[k] += t
    s = sum(Ys)
    return [y / s for y in Ys]


CASES = [
    # m [kg/s]: laminar (Re ~ 70), transitional, turbulent (Re ~ 1.4e5).
    (1e-5, "haaland"),
    (3e-3, "haaland"),
    (2e-2, "haaland"),
    (2e-2, "colebrook"),
    (-3e-3, "haaland"),  # reversed flow: the drop is odd in m
]


@pytest.mark.parametrize("head", [False, True])
@pytest.mark.parametrize(("m", "model"), CASES)
def test_the_drop_partials_match_resolved_differences(m: float, model: str, head: bool) -> None:
    Y = _air_with_fuel()
    r = _solver_tools.channel_residuals_and_jacobian(
        m, P_UP, P_UP, T_UP, Y, P_UP - 1.0, L, D, 1e-5, model, 1.0, head
    )
    scale = abs(r.dP_calc)
    fd_P = _richardson(lambda t: _drop(m, P_UP + t, T_UP, Y, model, head), 50.0)
    fd_T = _richardson(lambda t: _drop(m, P_UP, T_UP + t, Y, model, head), 0.05)
    assert r.d_dP_dP_static_up == pytest.approx(fd_P, rel=1e-8, abs=1e-12 * scale / P_UP)
    assert r.d_dP_dTt_up == pytest.approx(fd_T, rel=1e-8, abs=1e-12 * scale / T_UP)
    # Every species, present or absent: a row's composition columns are all
    # of them (an absent species arrives with an inflow at zero flow).
    for k in range(cb.num_species()):
        g = lambda t, k=k: _drop(m, P_UP, T_UP, _y_step(Y, k, t), model, head)  # noqa: E731
        fd_Y = _richardson(g, 1e-4) if Y[k] > 1e-3 else _one_sided(g, 1e-4)
        assert r.d_dP_dY_up[k] == pytest.approx(fd_Y, rel=1e-7, abs=1e-9 * scale), cb.species_name(
            k
        )


def test_the_partials_at_exactly_zero_flow_are_zero() -> None:
    """dP = 0 at m = 0 whatever the state, so its state partials are 0, not 0 * inf."""
    r = _solver_tools.channel_residuals_and_jacobian(
        0.0, P_UP, P_UP, T_UP, _air_with_fuel(), P_UP, L, D, 1e-5, "haaland", 1.0, True
    )
    assert r.d_dP_dP_static_up == 0.0
    assert r.d_dP_dTt_up == 0.0
    assert all(v == 0.0 for v in r.d_dP_dY_up)
