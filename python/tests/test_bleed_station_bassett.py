"""Bleed through a duct wall: the station is Bassett's separating straight run (#471).

Bassett, Winterbone & Pearson (2001), separating flow, straight-through loss
K2/K5 (Eq. 15): K = q^2 - 1.5 q + 0.5, q = m_out/m_in, theta- and
psi-independent. The station momentum balance with the bleed carrying
kappa = 0.75 of the arriving axial velocity reproduces it for every split, so
the closure is Bassett's and is scored here against his own measured K5
(Fig. 7c, theta 45/60/90, psi 1/3).
"""

from __future__ import annotations

import math
from pathlib import Path

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

DATA = Path(__file__).resolve().parents[2] / "validation" / "junction" / "data" / "bassett2001"
X_AIR = cb.species.dry_air()
Y_AIR = cb.mole_to_mass(X_AIR)


def _k5(q: float) -> float:
    """Bassett Eq. (15), written out here so the check is independent of C++."""
    return q * q - 1.5 * q + 0.5


def _station_K(q: float, kappa: float) -> float:
    """Total-pressure loss coefficient of one whole station from C++, on the
    arriving dynamic head, at constant density."""
    P, T, A, m_a = 1.2e5, 500.0, 4e-3, 0.4
    rho = cb.density(T, P, X_AIR)
    half = cb.station_half_drop(m_a, q * m_a, P, T, X_AIR, A, kappa)
    ua, ub = m_a / (rho * A), q * m_a / (rho * A)
    return (2.0 * half.dP + 0.5 * rho * (ua**2 - ub**2)) / (0.5 * rho * ua**2)


@pytest.mark.parametrize("q", [0.0, 0.3, 0.6, 0.9, 0.97, 1.0])
def test_the_bleed_station_is_bassett_k5(q: float) -> None:
    assert _station_K(q, cb.STATION_KAPPA_BLEED_BASSETT) == pytest.approx(_k5(q), abs=1e-12)


def test_normal_injection_is_not_bassett() -> None:
    """kappa = 0 (the impingement merge) is a different closure: it would
    over-state the bleed's regain."""
    assert abs(_station_K(0.9, cb.STATION_KAPPA_MERGE_NORMAL) - _k5(0.9)) > 0.1


def _fig7c() -> list[tuple[float, float]]:
    pts = []
    for f in sorted(DATA.glob("bassett_fig07c_K5_*_measured.csv")):
        for line in f.read_text().splitlines()[1:]:
            q, k = map(float, line.split(","))
            pts.append((q, k))
    return pts


def test_scored_against_bassett_fig7c() -> None:
    """Fidelity on the author's own data, all four series (theta 45/60/90,
    psi 1/3: theta- and psi-independent as Bassett states). Absolute error in
    K, since K crosses zero. The bleed-relevant end is q >= 0.6 (an effusion
    station bleeds a few percent)."""
    pts = _fig7c()
    assert len(pts) == 43
    err = np.array([_station_K(q, cb.STATION_KAPPA_BLEED_BASSETT) - k for q, k in pts])
    hi = np.array([_station_K(q, cb.STATION_KAPPA_BLEED_BASSETT) - k for q, k in pts if q >= 0.6])
    print(
        f"\nall {len(err)}: bias {err.mean():+.4f}, MAE {np.abs(err).mean():.4f}, "
        f"max {np.abs(err).max():.4f}"
        f"\nq>=0.6 {len(hi)}: bias {hi.mean():+.4f}, MAE {np.abs(hi).mean():.4f}, "
        f"max {np.abs(hi).max():.4f}"
    )
    # Measured 2026-10-08: all 43 bias +0.017, MAE 0.039, max 0.158 (q ~ 0,
    # full diversion, formula 0.49 vs 0.35 measured); q >= 0.6: bias -0.019,
    # MAE 0.021, max 0.046. Reported, not tuned. Near q -> 1 the misses sit at
    # q -> 1, where the data carry a +0.04 offset the formula (K5(1) = 0)
    # cannot have: a measurement datum, most likely the rig's own friction.
    assert np.abs(err).mean() < 0.045
    assert np.abs(hi).mean() < 0.025


def _bleed_net(
    regime: str, drive: float = 1.0e3, length: float = 1e-6, A: float = 1e-3
) -> FlowNetwork:
    """Reservoir -> duct -> bleed station -> duct -> outlet, the bleed to a
    sink as far below the outlet as the reservoir is above it."""
    g = FlowNetwork()
    for n, pt in (("in", 1.0e5 + drive), ("out", 1.0e5), ("lat", 1.0e5 - drive)):
        b = PressureBoundary(n)
        b.Pt, b.Tt, b.Y = pt, 300.0, Y_AIR
        g.add_node(b)
    g.add_node(PlenumNode("b"))
    kappa = cb.STATION_KAPPA_BLEED_BASSETT
    g.add_element(
        CrossflowSegmentElement(
            "s0",
            "in",
            "b",
            length=length,
            area=A,
            entry_K=0.0,
            to_kappa=kappa,
            next_seg="s1",
            regime=regime,
        )
    )
    g.add_element(
        CrossflowSegmentElement(
            "s1", "b", "out", length=length, area=A, from_kappa=kappa, prev_seg="s0", regime=regime
        )
    )
    g.add_element(
        OrificeElement(
            "bleed", "b", "lat", diameter=0.01, Cd=0.6, correlation="fixed", regime=regime
        )
    )
    return g


def test_a_solved_bleed_network_satisfies_bassett() -> None:
    """Near-frictionless and nearly isothermal: the outlet static must sit
    below the reservoir by exactly the entry dynamic head plus Bassett's
    station drop at the SOLVED split -- computed here from K5 alone,
    independently of the C++ station. Incompressible: the check below holds
    the density constant."""
    A = 1e-3
    r = NetworkSolver(_bleed_net("incompressible", A=A)).solve()
    assert r["__success__"], r.get("__message__")
    d = r["__element_diag__"]
    m_a, m_b = d["s0"]["m_dot"], d["s1"]["m_dot"]
    q = m_b / m_a
    assert 0.5 < q < 0.99
    rho = cb.density(300.0, r["b.P"], X_AIR)
    qa, qb = 0.5 * m_a**2 / (rho * A * A), 0.5 * m_b**2 / (rho * A * A)
    # Static: P_a = Pt_in - qa (entry, K_in = 0); station Pt loss K5 qa.
    P_a = 1.010e5 - qa
    P_b = P_a - (_k5(q) * qa - (qa - qb))
    # The outlet dumps (static continuity).
    assert P_b == pytest.approx(1.000e5, abs=2e-3 * (qa + qb))
    assert math.isfinite(r["b.P"])


def _split(regime: str, drive: float) -> tuple[float, float]:
    r = NetworkSolver(_bleed_net(regime, drive=drive, length=0.05)).solve()
    assert r["__success__"], r.get("__message__")
    return r["s0.m_dot"], r["s1.m_dot"]


def test_the_compressible_station_is_the_incompressible_one_at_low_mach() -> None:
    """The compressible segment row is a MASS FLOW (#481); the station drop is
    a pressure, and was once subtracted from that row: a 1 um segment then
    held 1 kPa at zero flow, the bleed reversed, and the flows came out 50x
    too small. Now the drop shifts the Fanno target, and the two regimes part
    as Ma^2: quartering the drive (Ma 0.12 -> 0.06) quarters the gap.

    A real segment length: with a near-frictionless one the friction drive
    sits inside the Fanno flux's 1 Pa regularisation, an absolute offset that
    does not shrink with Mach."""
    gaps = []
    for drive in (1.0e3, 250.0):
        (ci, di), (cc, dc) = _split("incompressible", drive), _split("compressible", drive)
        assert cc > 0.0 and dc > 0.0 and dc < cc  # forward, and the bleed takes its share
        gaps.append(max(abs(cc - ci) / ci, abs(dc - di) / di))
    assert gaps[0] < 0.12**2
    assert 3.0 < gaps[0] / gaps[1] < 5.0


def test_the_compressible_station_jacobian_matches_finite_differences() -> None:
    """Own flow, the neighbour segment's flow (through the station), and the
    station node's pressure, all through the shifted Fanno target."""

    s = NetworkSolver(_bleed_net("compressible", length=0.05))
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    x = np.array(r["__x_solution__"], dtype=float)
    J = s._residuals_and_jacobian(x)[1].toarray()
    # Flows step 1e-5, pressures 1e-7: the Fanno iteration is warm-started,
    # so a residual repeats only to ~1e-12.
    rel = np.array([1e-5 if n.endswith(".m_dot") else 1e-7 for n in s.unknown_names])
    fd = fd_jacobian(s, x, rel=rel, floor=1e-10)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    scale = np.maximum(np.max(np.abs(fd * cols), axis=1, keepdims=True), 1e-300)
    assert float(np.max(np.abs((J - fd) * cols) / scale)) < 1e-5


def _entry_net(entry_K: float, p_res: float, p_duct: float, regime: str) -> FlowNetwork:
    """Reservoir --[segment, entry_K]--> plenum boundary."""
    g = FlowNetwork()
    for n, pt in (("res", p_res), ("duct", p_duct)):
        b = PressureBoundary(n)
        b.Pt, b.Tt, b.Y = pt, 300.0, Y_AIR
        g.add_node(b)
    g.add_element(
        CrossflowSegmentElement(
            "s", "res", "duct", length=0.05, area=1e-3, entry_K=entry_K, regime=regime
        )
    )
    return g


def _entry_segment(entry_K: float, p_res: float, p_duct: float, regime: str) -> float:
    """The solved flow of ``_entry_net``."""
    res = NetworkSolver(_entry_net(entry_K, p_res, p_duct, regime)).solve()
    assert res["__success__"], res.get("__message__")
    return float(res["s.m_dot"])


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_the_entry_loss_applies_to_inflow_only(regime: str) -> None:
    """Into the duct, the reservoir pays the acceleration and the entry
    loss, so K throttles the flow. Reversed, the duct discharges into the
    reservoir at static continuity -- an exit, the segment's own convention
    -- so K must not matter. The mirrored m|m| form charged the outflow
    (1 + K) dynamic heads instead (#705)."""
    fwd = [_entry_segment(k, 1.01e5, 1.0e5, regime) for k in (0.0, 0.5)]
    rev = [_entry_segment(k, 1.0e5, 1.01e5, regime) for k in (0.0, 0.5)]
    assert fwd[0] > fwd[1] > 0.0
    assert rev[0] < 0.0
    assert rev[1] == pytest.approx(rev[0], rel=1e-9)
    # Reversed, the plenum boundary at rest is the reservoir: its gas pays the
    # same acceleration (#802). Before, the reversed flow paid only friction
    # and ran ~5x the forward one. With K = 0 the two directions are mirror
    # images, in both regimes: the incompressible head takes rho at the
    # duct's far-end static P either way (#1136; at the reservoir's own P it
    # was 0.49 % off the mirror).
    assert abs(rev[0]) == pytest.approx(fwd[0], rel=1e-12)


def test_reversed_out_of_an_exit_the_regimes_agree() -> None:
    """Low Mach, the compressible duct ends and the incompressible momentum
    terms charge the same head where the gas leaves an exit node at rest. Before
    the incompressible reversed flow paid friction only, ~5x the compressible
    one. They agreed to 0.53 % with the head at the reservoir's (stagnation)
    density, and agree to 0.04 % with it at the duct's far-end static one
    (#1136), the forward entry's 0.04 %."""
    inc = _entry_segment(0.5, 1.0e5, 1.01e5, "incompressible")
    com = _entry_segment(0.5, 1.0e5, 1.01e5, "compressible")
    assert inc < 0.0
    assert com == pytest.approx(inc, rel=1e-3)


def test_the_reversed_exit_inflow_jacobian_matches_finite_differences() -> None:
    """The incompressible reversed-exit head, through the solver's assembly."""

    s = NetworkSolver(_entry_net(0.5, 1.0e5, 1.01e5, "incompressible"))
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    x = np.array(r["__x_solution__"], dtype=float)
    assert x[s.unknown_names.index("s.m_dot")] < 0.0
    J = s._residuals_and_jacobian(x)[1].toarray()
    fd = fd_jacobian(s, x, floor=1e-10)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    scale = np.maximum(np.max(np.abs(fd * cols), axis=1, keepdims=True), 1e-300)
    assert float(np.max(np.abs((J - fd) * cols) / scale)) < 1e-6
