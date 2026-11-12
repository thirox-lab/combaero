"""A crossflow station takes the axial momentum its joining jets declare (#813).

A merge chamber already read each side stream's ``injection_momentum``
(m w cos(theta), #471); a ``CrossflowSegmentElement`` station took the
side stream's axial velocity as ``kappa`` times the receiving duct's own,
fixed when the chain was built, and never read the hole's angle. Both now
read the same declaration. A normal jet (90 deg: impingement rows, plain
holes) declares none, so every station that existed before is untouched --
the Florschuetz and CR-3630 scorecards are the guard.

The second-law ledger groups a momentum-bearing hole with the chain: closed
at the station's Pt, the hole's jet momentum would be booked as pressure the
chain created.
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
from combaero.network.second_law import SecondLawReport, _sigma

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
AREA = 2e-3


def _net(theta: float, regime: str, p_res: float = 1.03e5, p_sup: float = 1.06e5) -> FlowNetwork:
    """Reservoir --s0--> station --s1--> outlet; one hole from a supply
    plenum joins at the station at ``theta`` to the duct axis."""
    g = FlowNetwork()
    for n, pt in (("res", p_res), ("sup", p_sup), ("out", 1.0e5)):
        b = PressureBoundary(n)
        b.Pt, b.Tt, b.Y = pt, 300.0, Y_AIR
        g.add_node(b)
    g.add_node(PlenumNode("st"))
    g.add_element(
        CrossflowSegmentElement(
            "s0", "res", "st", length=0.02, area=AREA, entry_K=0.0, to_kappa=0.0,
            next_seg="s1", regime=regime,
        )
    )  # fmt: skip
    g.add_element(
        CrossflowSegmentElement(
            "s1", "st", "out", length=0.02, area=AREA, from_kappa=0.0, prev_seg="s0",
            regime=regime,
        )
    )  # fmt: skip
    g.add_element(
        OrificeElement(
            "h", "sup", "st", Cd=0.8, diameter=0.02, regime=regime, correlation="fixed",
            injection_angle_deg=theta,
        )
    )  # fmt: skip
    return g


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    return s, r


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_the_station_balances_the_jets_momentum(regime: str) -> None:
    """Over the station, (P_a - P_b) A = m_b u_b - m_a u_a - J, with J the
    hole's m w cos(theta) at its supply state and the station's static P.
    The two half-terms are rebuilt here from the solved state; their sum is
    the station drop the two segments carry."""
    s, r = _solve(_net(30.0, regime))
    m_a, m_b, m_h = (float(r[k]) for k in ("s0.m_dot", "s1.m_dot", "h.m_dot"))
    P = float(r["st.P"])
    x = np.array(r["__x_solution__"], dtype=float)
    s._residuals(x)
    st = s._get_node_state(s.network.nodes["st"], x)
    X = list(cb.mass_to_mole(st.Y))
    jet = cb.jet_impulse(m_h, 1.06e5, 300.0, P, X)
    J = jet.J * math.cos(math.radians(30.0))
    assert J > 0.0
    # The segments' two halves, each at the station's moving static state.
    halves = 0.0
    for eid in ("s0", "s1"):
        el = s.network.elements[eid]
        a = s._get_node_state(s.network.nodes[el.from_node], x)
        b = s._get_node_state(s.network.nodes[el.to_node], x)
        a.m_dot = b.m_dot = float(r[f"{eid}.m_dot"])
        for t in el._momentum_terms(a, b, s._network_flows(el, x)):
            if t[1] == "st":
                halves += t[0].dP
    # The same halves without the jet: the station's kappa = 0 drop.
    rho_terms = 0.0
    for eid in ("s0", "s1"):
        el = s.network.elements[eid]
        a = s._get_node_state(s.network.nodes[el.from_node], x)
        b = s._get_node_state(s.network.nodes[el.to_node], x)
        a.m_dot = b.m_dot = float(r[f"{eid}.m_dot"])
        saved = el._station_jets
        el._station_jets = lambda node, station: (0.0, {})  # type: ignore[method-assign]
        try:
            for t in el._momentum_terms(a, b, s._network_flows(el, x)):
                if t[1] == "st":
                    rho_terms += t[0].dP
        finally:
            el._station_jets = saved  # type: ignore[method-assign]
    assert m_b == pytest.approx(m_a + m_h, rel=1e-9)
    assert halves == pytest.approx(rho_terms - J / AREA, rel=1e-9, abs=1e-9)


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_a_more_axial_jet_pumps_the_crossflow(regime: str) -> None:
    """The jets' axial momentum lowers the station pressure and draws more
    from the reservoir; 90 deg brings none."""
    runs = [_solve(_net(th, regime))[1] for th in (90.0, 60.0, 30.0, 0.0)]
    m0 = [float(r["s0.m_dot"]) for r in runs]
    P = [float(r["st.P"]) for r in runs]
    assert all(b > a for a, b in zip(m0, m0[1:], strict=False))
    assert all(b < a for a, b in zip(P, P[1:], strict=False))


def test_a_normal_jet_declares_nothing() -> None:
    """90 deg: the station never sees the hole, the term carries no extra
    partials -- the stations of an impingement array are unchanged."""
    s, r = _solve(_net(90.0, "compressible"))
    x = np.array(r["__x_solution__"], dtype=float)
    s._residuals(x)
    el = s.network.elements["s0"]
    st = s._get_node_state(s.network.nodes["st"], x)
    assert el._station_jets(s.network.nodes["st"], st) == (0.0, {})


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.parametrize("theta", [30.0, 0.0])
def test_the_jet_jacobian_matches_finite_differences(regime: str, theta: float) -> None:
    """At the root and nearby. Not near a duct's flow reversal: there m ~
    sqrt(dPt) and a central difference of 1e-6 relative is truncation."""
    s, r = _solve(_net(theta, regime))
    x0 = np.array(r["__x_solution__"], dtype=float)
    for x in (x0, x0 * (1.0 + 0.002 * np.sin(np.arange(len(x0))))):
        J = s._residuals_and_jacobian(x)[1].toarray()
        fd = fd_jacobian(s, x, floor=1e-10)
        cols = np.maximum(np.abs(x), 1e-12)[None, :]
        scale = np.maximum(np.max(np.abs(fd * cols), axis=1, keepdims=True), 1e-300)
        assert float(np.max(np.abs((J - fd) * cols) / scale)) < 1e-4


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_the_ledger_credits_the_jet_to_the_chain(regime: str) -> None:
    """A jet that alone pumps the chain (reservoir at the outlet's Pt). The
    hole is a member of the chain's control volume and the root passes.
    Closing the hole at the station's Pt instead, as for a normal jet, books
    the jet's momentum as pressure the chain created: negative."""
    s, r = _solve(_net(0.0, regime, p_res=1.0e5, p_sup=1.3e5))
    assert float(r["s0.m_dot"]) > 0.0  # drawn in by the jet
    sl = r["__second_law__"]
    assert sl["h"]["group"] == "h+s0+s1"
    assert sl["h"]["verdict"] == "pass"
    report = SecondLawReport(s, np.array(r["__x_solution__"], dtype=float), r)
    assert _sigma(report._chain(["s0", "s1"], []))["sigma_norm_high"] < 0.0


def test_a_choked_jets_jacobian_is_exact() -> None:
    """A jet from an orifice-fed plenum at ~4 bar into a 1 bar chain is
    choked; the plenum's Pt is an unknown, so its momentum's dJ/dPt (C++
    jet_impulse) enters the station rows. Analytic since #816: on main those
    entries read 2% off FD (central differences inside C++ around a sonic
    point noisy to ~1e-8 of Pt, golden-section), here 1e-7."""
    g = FlowNetwork()
    for n, pt in (("res", 1.0e5), ("hp", 4.0e5), ("out", 1.0e5)):
        b = PressureBoundary(n)
        b.Pt, b.Tt, b.Y = pt, 300.0, Y_AIR
        g.add_node(b)
    g.add_node(PlenumNode("sup"))
    g.add_node(PlenumNode("st"))
    g.add_element(OrificeElement("feed", "hp", "sup", Cd=0.8, diameter=0.04, correlation="fixed"))
    g.add_element(
        CrossflowSegmentElement(
            "s0", "res", "st", length=0.02, area=AREA, entry_K=0.0, to_kappa=0.0, next_seg="s1"
        )
    )
    g.add_element(
        CrossflowSegmentElement(
            "s1", "st", "out", length=0.02, area=AREA, from_kappa=0.0, prev_seg="s0"
        )
    )
    g.add_element(
        OrificeElement(
            "h", "sup", "st", Cd=0.8, diameter=0.02, correlation="fixed", injection_angle_deg=30.0
        )
    )
    s, r = _solve(g)
    x = np.array(r["__x_solution__"], dtype=float)
    st = s._get_node_state(s.network.nodes["st"], x)
    jet = cb.jet_impulse(
        float(r["h.m_dot"]), float(r["sup.P"]), 300.0, float(r["st.P"]), list(st.X)
    )
    assert jet.choked
    J = s._residuals_and_jacobian(x)[1].toarray()
    fd = fd_jacobian(s, x, floor=1e-12)
    col = s.unknown_names.index("sup.Pt")
    big = np.abs(fd[:, col]) > 1e-9
    assert big.sum() >= 3  # the plenum's own row and both segments' station rows
    rel = np.abs(J[big, col] - fd[big, col]) / np.abs(fd[big, col])
    assert float(rel.max()) < 1e-5
