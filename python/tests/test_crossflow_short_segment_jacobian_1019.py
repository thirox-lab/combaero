"""Short crossflow segments: the Fanno kernel's choked derivatives (#1019).

A 1 um ``CrossflowSegmentElement`` (the stress family ``bleed_station``)
chokes at its sonic flux less ~1e-6: the choked flux of a duct of length L
sits 1 - M_in ~ sqrt(f L / D) below sonic. The choked branch differentiated
L*(G; Pt0, Tt0, Y) = L through a re-solved inlet Mach with 1e-6 steps, which
then crossed sonic (M_in clamps to 1, L* to 0): dG/dTt0 came out 17% high at
L = 1 um, and the station rows' columns that move the station's Tt (the
segments' and the hot injection's flows, through the mixing relay) 7-17% off.

Below the 1 Pa drive floor the blend's exponent sigma stays held fixed, a
measured approximation (see fanno_channel_flow): its bound is recorded here.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    CrossflowSegmentElement,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

X_AIR = cb.species.dry_air()
Y_AIR = cb.mole_to_mass(X_AIR)
D = math.sqrt(4e-3 / math.pi)  # area 1e-3, as the stress family's segments
DIRS = [[float(j == k) for j in range(len(X_AIR))] for k, x in enumerate(X_AIR) if x > 0.0]


def _flow(L: float, Pt0: float, Tt0: float, P: float, X=X_AIR, derivatives: bool = True):
    return cb.fanno_channel_flow(
        Pt0,
        Tt0,
        X,
        P,
        True,
        L,
        D,
        0.0,
        "haaland",
        1.0,
        -1.0,
        DIRS if derivatives else [],
        derivatives,
    )


@pytest.mark.parametrize("L", [0.5, 0.05, 1e-3, 1e-5, 1e-6, 1e-7])
def test_choked_flux_derivatives_hold_down_to_a_micron(L: float) -> None:
    """J_ij |x_j| / G against central differences, to 1e-6 at every length.
    Before the fix: 1.2e-5 at 10 um, 8.6e-2 (dG/dTt0) at 1 um."""
    Pt0, Tt0, P = 517338.6, 468.26, 1.0e5  # the issue's reversed station duct
    f = _flow(L, Pt0, Tt0, P)
    assert f.choked
    h = 1e-7
    fd_T = (
        _flow(L, Pt0, Tt0 * (1 + h), P, derivatives=False).G
        - _flow(L, Pt0, Tt0 * (1 - h), P, derivatives=False).G
    ) / (2 * h * Tt0)
    fd_P = (
        _flow(L, Pt0 * (1 + h), Tt0, P, derivatives=False).G
        - _flow(L, Pt0 * (1 - h), Tt0, P, derivatives=False).G
    ) / (2 * h * Pt0)
    assert abs(f.dG_dTt0 - fd_T) * Tt0 / f.G < 1e-6
    assert abs(f.dG_dPt0 - fd_P) * Pt0 / f.G < 1e-6
    hy = 1e-6
    for i, v in enumerate(DIRS):

        def g(s: float, v=v) -> float:
            Ys = [y + s * vk for y, vk in zip(Y_AIR, v, strict=True)]
            return _flow(L, Pt0, Tt0, P, cb.mass_to_mole(Ys), derivatives=False).G

        fd_v = (g(hy) - g(-hy)) / (2 * hy)
        assert abs(f.dG_ddir[i] - fd_v) / f.G < 1e-6


def test_the_drive_floor_holds_sigma_fixed_within_its_measured_bound() -> None:
    """Half the floor of drive on a 5 cm duct, where f moves with Re and so
    sigma with Tt0: the frozen sigma leaves dG/dTt0 Tt0 / G 4e-3 off (the
    same for K steps from 0.1 to 1, so it is not the oracle's noise).
    Differencing sigma tripled a short-segment network's Jacobian cost and
    slowed its solves 2-3x, so it stays an approximation, bounded here."""
    L, Pt0, Tt0 = 0.05, 5.0e5, 300.0
    P = Pt0 - 0.5 * cb.FANNO_FLOW_DRIVE_FLOOR
    f = _flow(L, Pt0, Tt0, P)
    for hT in (1.0, 0.1):
        fd = (
            _flow(L, Pt0, Tt0 + hT, P, derivatives=False).G
            - _flow(L, Pt0, Tt0 - hT, P, derivatives=False).G
        ) / (2 * hT)
        assert abs(f.dG_dTt0 - fd) * Tt0 / f.G < 1e-2


def _station_net(length: float) -> FlowNetwork:
    """The stress family's bleed_station seed 20: a hot (600 K) normal
    injection between two segments, 5.1 / 5.0 / 5.2 bar. The station's Tt
    mixes the cold through-flow with the jet, so it moves with both flows."""
    g = FlowNetwork()
    for n, pt, tt in (("in", 5.1e5, 300.0), ("out", 5.0e5, 300.0), ("lat", 5.2e5, 600.0)):
        b = PressureBoundary(n)
        b.Pt, b.Tt, b.Y = pt, tt, Y_AIR
        g.add_node(b)
    g.add_node(PlenumNode("b"))
    kappa = cb.STATION_KAPPA_MERGE_NORMAL
    g.add_element(
        CrossflowSegmentElement(
            "s0",
            "in",
            "b",
            length=length,
            area=1e-3,
            entry_K=0.0,
            to_kappa=kappa,
            next_seg="s1",
            regime="compressible",
        )
    )
    g.add_element(
        CrossflowSegmentElement(
            "s1",
            "b",
            "out",
            length=length,
            area=1e-3,
            from_kappa=kappa,
            prev_seg="s0",
            regime="compressible",
        )
    )
    g.add_element(
        OrificeElement(
            "lat_o", "lat", "b", diameter=0.03, Cd=0.6, correlation="fixed", regime="compressible"
        )
    )
    return g


# The solver's best iterate on seed 20 (issue #1019): s0's kernel runs
# reversed and choked there.
ISSUE_X = {
    "b.P": 508091.4272960788,
    "b.Pt": 508091.4378957084,
    "s0.m_dot": 0.11780735172365611,
    "s1.m_dot": 0.27424118894773436,
    "lat_o.m_dot": 0.14693563628770168,
}


@pytest.mark.parametrize("length", [0.5, 0.05, 1e-3, 1e-5, 1e-6])
def test_station_rows_match_finite_differences_over_segment_length(length: float) -> None:
    """J_ij |x_j| against central differences, scaled by each row's largest
    entry. Before the fix: 1.4e-2 at 1 um (s0's row in s0.m_dot and
    lat_o.m_dot, J 1.727 / -0.600 against FD 1.617 / -0.511).

    The oracle is Richardson-extrapolated central differences, relative step
    1e-7 and 2e-7: s1's drive sits inside the 1 Pa floor, whose cubic has
    curvature on the Pa scale, so a plain 1e-6 step (0.5 Pa of b.P) reads
    0.7% truncation error there, falling as the step squared, and a step
    small enough to drop it (1e-8) sits on the residual's rounding noise
    (~1.5e-6)."""
    s = NetworkSolver(_station_net(length))
    s.network.resolve_all_topology()
    s._build_x0()
    x = np.array([ISSUE_X[n] for n in s.unknown_names], dtype=float)
    J = s._residuals_and_jacobian(x)[1].toarray()
    d1 = fd_jacobian(s, x, rel=1e-7, floor=0.0)
    d2 = fd_jacobian(s, x, rel=2e-7, floor=0.0)
    fd = (4.0 * d1 - d2) / 3.0
    cols = np.abs(x)[None, :]
    scale = np.maximum(np.max(np.abs(fd * cols), axis=1, keepdims=True), 1e-300)
    assert float(np.max(np.abs((J - fd) * cols) / scale)) < 1e-6
