"""A channel's laminar drop is Hagen-Poiseuille, linear in its flow.

The network friction kernels floored the laminar branch at Re 64 with the
turbulent fits' guard, so below that f saturated at 1 and the drop went as
m^2, with zero slope at zero flow; Colebrook had no laminar branch at all.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero import _solver_tools

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
L, D, P, TT = 0.5, 0.05, 2e5, 300.0


def _drop(m: float, model: str) -> tuple[float, float]:
    r = _solver_tools.channel_residuals_and_jacobian(
        m, P, P, TT, Y_AIR, P - 1e3, L, D, 0.0, model, 1.0, False
    )
    return float(r.dP_calc), float(r.d_dP_d_mdot)


@pytest.mark.parametrize("model", ["haaland", "serghides", "colebrook"])
def test_the_low_flow_drop_is_hagen_poiseuille(model: str) -> None:
    X = cb.mass_to_mole(Y_AIR)
    rho = cb.density(TT, P, X)
    mu = cb.transport_state(TT, P, X).mu
    A = math.pi * D * D / 4.0
    slope_hp = 32.0 * mu * L / (rho * A * D * D)  # dP/dm, f Re = 64
    for m in (1e-8, 1e-6, 1e-4):
        dP, dP_dm = _drop(m, model)
        Re = m * D / (mu * A)
        assert Re < 300.0
        assert dP == pytest.approx(slope_hp * m, rel=2e-3), (model, m)
        assert dP_dm == pytest.approx(slope_hp, rel=2e-3), (model, m)


@pytest.mark.parametrize("model", ["haaland", "serghides", "colebrook"])
def test_the_slope_at_exactly_zero_flow_is_the_laminar_limit(model: str) -> None:
    """At m = 0, |m| f is 0 * inf: the kernel's slope read 0 there (#1156)."""
    X = cb.mass_to_mole(Y_AIR)
    rho = cb.density(TT, P, X)
    mu = cb.transport_state(TT, P, X).mu
    A = math.pi * D * D / 4.0
    # Hagen-Poiseuille times the blend's laminar share at Re = 0,
    # 1 - w(0) = 1 - (1 + tanh(-2300 / 400)) / 2: tight enough to tell 64 from
    # 63.99935 (w(0) = 1.0e-5), which a bare Hagen-Poiseuille bound cannot.
    w0 = 0.5 * (1.0 + math.tanh(-2300.0 / 400.0))
    slope_hp = 32.0 * mu * L / (rho * A * D * D) * (1.0 - w0)
    dP0, dP_dm0 = _drop(0.0, model)
    assert dP0 == 0.0
    assert dP_dm0 == pytest.approx(slope_hp, rel=1e-9), model
    # Continuous through zero: the one-sided slopes just off it agree.
    for m in (1e-8, -1e-8):
        assert _drop(m, model)[1] == pytest.approx(dP_dm0, rel=1e-7), (model, m)


def test_petukhov_slope_at_zero_flow_is_finite() -> None:
    """Petukhov has no laminar branch (f -> 0 as the flow stops); its slope at
    m = 0 is 0, where 0 * inf read NaN (#1156)."""
    assert _drop(0.0, "petukhov") == (0.0, 0.0)
