"""A duct's heat transfer is evaluated at its own static state (#636).

ChannelElement and PressureLossElement handed the correlation the feed
node's T, which for a plenum or a chamber is the stagnation T, as the duct's
static T. T_aw = T + r v^2 / 2cp then came out ABOVE the stagnation
temperature, by about (1 + r)(Tt - T): 1602.3 K for a 1600 K stream at
Mach 0.1 (30 K in a GUI case). The static state now comes from the feed's
(Tt, Pt) and the mass flux by the isentropic acceleration (duct_static_state).
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FlowNetwork,
    NetworkSolver,
    PressureBoundary,
)
from combaero.network.components import NetworkMixtureState, SmoothModel

X_AIR = cb.species.dry_air()
Y_AIR = list(cb.mole_to_mass(X_AIR))
T0, P0 = 1600.0, 2.0e5


def _flux_at_mach(M: float) -> float:
    """G of the isentropic duct state at Mach M (bisection on the subsonic branch)."""
    lo, hi = 1e-6, 0.999
    G_hi = cb.duct_static_state(T0, P0, 1e6, X_AIR)
    G_star = G_hi.rho * G_hi.M * cb.speed_of_sound(G_hi.T, X_AIR)
    g_lo, g_hi = 0.0, G_star
    for _ in range(80):
        g = 0.5 * (g_lo + g_hi)
        if cb.duct_static_state(T0, P0, g, X_AIR).M < M:
            g_lo = g
        else:
            g_hi = g
    assert lo < M < hi
    return 0.5 * (g_lo + g_hi)


@pytest.mark.parametrize("G", [5.0, 50.0, 150.0])
def test_the_duct_state_closes_energy_and_its_stagnation_pressure(G: float) -> None:
    s = cb.duct_static_state(T0, P0, G, X_AIR)
    assert not s.choked
    v = G / s.rho
    assert cb.h_mass(T0, X_AIR) - cb.h_mass(s.T, X_AIR) == pytest.approx(0.5 * v * v, rel=1e-9)
    assert cb.Pt_from_static(s.P, s.T, s.M, X_AIR) == pytest.approx(P0, rel=1e-10)


@pytest.mark.parametrize("G", [5.0, 150.0])
def test_the_duct_state_partials_match_central_differences(G: float) -> None:
    s = cb.duct_static_state(T0, P0, G, X_AIR)
    for arg, step in (("Tt", 1e-2), ("Pt", 1.0), ("G", 1e-4 * G)):

        def at(d: float, arg: str = arg) -> cb.DuctStaticState:
            kw = {"Tt": T0, "Pt": P0, "G": G}
            kw[arg] += d
            return cb.duct_static_state(kw["Tt"], kw["Pt"], kw["G"], X_AIR)

        sp, sm = at(step), at(-step)
        for q in ("T", "P", "rho"):
            fd = (getattr(sp, q) - getattr(sm, q)) / (2 * step)
            assert getattr(s, f"d{q}_d{arg}") == pytest.approx(fd, rel=1e-5, abs=1e-12), (q, arg)


def test_beyond_choke_the_duct_state_is_sonic_and_finite() -> None:
    s = cb.duct_static_state(T0, P0, 1e4, X_AIR)
    assert s.choked and pytest.approx(1.0, abs=1e-6) == s.M
    assert s.dT_dG == 0.0 and s.dP_dG == 0.0
    bad = cb.duct_static_state(T0, -1.5e6, 40.0, X_AIR)  # an unphysical iterate
    assert all(math.isfinite(v) for v in (bad.T, bad.P, bad.rho))


@pytest.mark.parametrize("M", [0.1, 0.3, 0.5, 0.7])
def test_a_channels_t_aw_lies_below_the_stagnation_temperature(M: float) -> None:
    """At every Mach: T_static < T_aw < Tt, and T_aw is the recovery value at
    the exact isentropic static state."""
    D = 0.02
    A = math.pi * D * D / 4.0
    G = _flux_at_mach(M)
    el = ChannelElement(
        "c", "a", "b", length=0.3, diameter=D, surface=ConvectiveSurface(model=SmoothModel())
    )
    el._fill_convective_area()
    state = NetworkMixtureState(P=P0, Pt=P0, Tt=T0, m_dot=G * A, Y=Y_AIR)
    r = el.htc_and_T(state)
    s = cb.duct_static_state(T0, P0, G, X_AIR)
    assert s.T < r.T_aw < T0, (s.T, r.T_aw)
    v = G / s.rho
    expected = cb.T_adiabatic_wall(s.T, v, s.T, s.P, X_AIR)
    assert r.T_aw == pytest.approx(expected, rel=1e-12)


def _plenum_fed_wall_network() -> FlowNetwork:
    from combaero.network import PlenumNode, ThermalWall, WallLayer

    g = FlowNetwork()
    for nid, Pt, Tt in (("hi", 1.6e5, 1500.0), ("ci", 1.6e5, 400.0)):
        g.add_node(PressureBoundary(nid, Pt=Pt, Tt=Tt))
    for nid in ("ho", "co"):
        g.add_node(PressureBoundary(nid, Pt=1.0e5, Tt=300.0))
    for nid in ("hp", "cp"):
        g.add_node(PlenumNode(nid))
    from combaero.network import LosslessConnectionElement

    # Lossless feeds: the compressible orifice row carries its own, unrelated
    # 2e-3 Jacobian error at its m_dot column (identical before this change).
    g.add_element(LosslessConnectionElement("ho1", "hi", "hp"))
    g.add_element(LosslessConnectionElement("co1", "ci", "cp"))
    for eid, a, b in (("hch", "hp", "ho"), ("cch", "cp", "co")):
        g.add_element(
            ChannelElement(
                eid,
                a,
                b,
                length=0.5,
                diameter=0.015,
                regime="compressible",
                surface=ConvectiveSurface(model=SmoothModel()),
            )
        )
    g.add_wall(ThermalWall("w", "hch", "cch", layers=[WallLayer(0.001, 20.0)]))
    return g


def test_a_plenum_fed_wall_network_has_an_exact_jacobian() -> None:
    """A plenum-fed pair of heated ducts keeps an exact Jacobian with the
    duct-state chain. (The Pt relay columns, dh_dPt and dT_aw_dPt, are pinned
    by test_wall_relay_pressure_496, which fails without them.)"""
    s = NetworkSolver(_plenum_fed_wall_network())
    r = s.solve()
    assert r["__success__"], r["__message__"]
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    # RESIDUAL-ROW-SCALING: compare J_ij |x_j|.
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    assert float(np.max(np.abs(Js - fds) / scale)) < 1e-6
    assert r["hp.Pt"] > r["hp.P"] or r["hp.Pt"] == pytest.approx(r["hp.P"])


def test_past_choke_the_duct_heat_transfer_stays_bounded() -> None:
    """A solver iterate beyond the duct's capacity: the sonic state and the
    sonic velocity, so T_aw stays below Tt and the heat bounded; it used to
    take G / rho at the iterate's G (53 caught mixer failures in the stress
    sweep)."""
    D = 0.005
    A = math.pi * D * D / 4.0
    el = ChannelElement(
        "c", "a", "b", length=0.3, diameter=D, surface=ConvectiveSurface(model=SmoothModel())
    )
    el._fill_convective_area()
    for G in (1e3, 1e4, 1e5):
        state = NetworkMixtureState(P=P0, Pt=P0, Tt=T0, m_dot=G * A, Y=Y_AIR)
        r = el.htc_and_T(state)
        assert cb.duct_static_state(T0, P0, G, X_AIR).choked
        assert r.T_aw < T0
        assert r.dh_dmdot == 0.0 and r.dT_aw_dmdot == 0.0
