"""A lone incompressible duct's Jacobian at exactly zero flow (#1156).

A ChannelElement, or a bare CrossflowSegmentElement, alone between two
PressureBoundaries, incompressible: at m = 0 the analytic d(row)/dm was 0,
while a difference read the laminar friction slope (-6.5e4 Pa s/kg on the 1 m,
1 cm duct here). The C++ kernel's slope ``|m| (f + |m| f' dRe/dm / 2)`` is
0 * inf at m = 0 (f -> 64 / Re); it now takes the limit, f Re at Re = 0
(``friction_f_re_at_zero``). F is unchanged. Below the laminar branch's
re_floor (Re 1e-6, |m| ~ 1e-13 kg/s here) the drop turns quadratic; a
difference step of 1e-9 kg/s is far above it, and a central difference
converges on the limit as h (1.1e-3, 1.1e-4, 1.1e-5, 1.1e-6 relative at h =
1e-6 .. 1e-9 kg/s, the m|m| terms of the flow blend) until the floor, below
h ~ 1e-11.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian, scaled_error

import combaero as cb
from combaero.network import (
    ChannelElement,
    CrossflowSegmentElement,
    FlowNetwork,
    NetworkSolver,
    PressureBoundary,
)

TOL = 1e-5


def _boundary(name: str, Pt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, 300.0, cb.mole_to_mass(cb.species.dry_air())
    return b


def _fresh(kind: str, dp: float) -> NetworkSolver:
    net = FlowNetwork()
    net.add_node(_boundary("a", 1.0e5 + dp))
    net.add_node(_boundary("b", 1.0e5))
    if kind == "channel":
        el = ChannelElement("c", "a", "b", length=1.0, diameter=0.01, regime="incompressible")
    else:
        el = CrossflowSegmentElement(
            "c", "a", "b", length=1.0, area=7.85e-5, regime="incompressible"
        )
    net.add_element(el)
    s = NetworkSolver(net)
    s._build_x0()
    el.resolve_topology(net)
    return s


def _jacobian(s: NetworkSolver, x: np.ndarray) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return s._residuals_and_jacobian(x)[1].toarray()


@pytest.mark.parametrize("dp", [1000.0, 0.0])
@pytest.mark.parametrize("kind", ["channel", "segment"])
@pytest.mark.parametrize("m", [0.0, 1e-7, -1e-7])
def test_the_flow_column_meets_fd_at_and_near_zero_flow(kind: str, dp: float, m: float) -> None:
    s = _fresh(kind, dp)
    x = np.array([m])
    J = _jacobian(s, x)
    fd = fd_jacobian(s, x)
    assert J[0, 0] < 0.0
    assert float(np.max(scaled_error(J, fd, x))) < TOL, (J, np.asarray(fd))


@pytest.mark.parametrize("kind", ["channel", "segment"])
def test_one_newton_step_from_zero_flow_moves_to_the_laminar_estimate(kind: str) -> None:
    """J = 0 made the zero start a singular Newton step. From m = 0 the step
    is now -F / J, the laminar flow: on a 1 Pa drive (Re ~ 100) within 3% of
    the root, the exit head the residual adds (2% of the drive here)."""
    s = _fresh(kind, 1.0)
    x0 = np.array([0.0])
    F0 = np.asarray(s._residuals(x0))
    assert F0[0] == 1.0  # Pt_a - Pt_b - dP(0), dP(0) = 0
    m1 = float(-F0[0] / _jacobian(s, x0)[0, 0])
    assert 0.0 < m1 < 1e-4
    assert abs(float(np.asarray(s._residuals(np.array([m1])))[0])) < 3e-2


@pytest.mark.parametrize("dp", [1.0, 100.0, 1e4])
@pytest.mark.parametrize("kind", ["channel", "segment"])
def test_a_solve_from_the_zero_flow_seed_converges(kind: str, dp: float) -> None:
    """From the all-zero-flow seed (``_zero_flow_seed``, #558) hybr stopped
    on its progress test near m ~ 1e-8 with J = 0 there: 5 of these 6 failed
    before #1156 (segment at 1e4 Pa got through). The cold start, and equal
    boundary pressures, converged before and still do."""
    s = _fresh(kind, dp)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=20, x0=s._zero_flow_seed())
        cold = _fresh(kind, dp).solve(timeout=20)
    assert r["__success__"], r.get("__message__")
    assert cold["__success__"], cold.get("__message__")
    assert r["c.m_dot"] == pytest.approx(cold["c.m_dot"], rel=1e-8)
