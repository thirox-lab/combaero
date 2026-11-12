"""A ribbed or pin-fin drop at a reservoir end reads the duct's state (#926).

Ribbed and pin-fin surfaces own the channel's pressure drop and read the from
node's (Tt, P) for it. On a node at rest (a plenum, a combustor, a
PressureBoundary or MassFlowBoundary) that is the reservoir, P = Pt, whose
density rho(Tt, Pt) is about M^2 / 2 above the duct's static density (13% at
Mach 0.5), so the drop came out low by as much. At such an end the drop now
reads the duct's own inlet static state, the isentropic acceleration from the
node's (Tt, Pt) to the duct's mass flux, with its analytic partials. A moving
node (a momentum chamber) keeps its own static state, unchanged.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)
from combaero.network.components import NetworkMixtureState, PinFinModel, RibbedModel

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
X_AIR = list(cb.mass_to_mole(Y_AIR))
D, L = 0.01, 0.1
AREA = math.pi * D * D / 4.0


def _surface(kind: str) -> ConvectiveSurface:
    if kind == "ribbed":
        return ConvectiveSurface(area=0.01, model=RibbedModel(e_D=0.06, p_e=10.0, W_H=1.0))
    return ConvectiveSurface(area=0.01, model=PinFinModel())


def _channel(eid: str, a: str, b: str, kind: str, regime: str) -> ChannelElement:
    return ChannelElement(eid, a, b, length=L, diameter=D, regime=regime, surface=_surface(kind))


def _reservoir_net(kind: str, regime: str, feed: str) -> FlowNetwork:
    g = FlowNetwork()
    if feed == "pressure":
        g.add_node(PressureBoundary("a", Pt=2.6e5, Tt=600.0, Y=Y_AIR))
    elif feed == "mass_flow":
        g.add_node(MassFlowBoundary("a", m_dot=0.02, Tt=400.0, Y=Y_AIR))
    else:
        # A plenum whose Pt, P and Tt are solver unknowns / relays, fed by an
        # incompressible duct (its rows are FD-exact, so every row is checked).
        g.add_node(PressureBoundary("s", Pt=2.6e5, Tt=600.0, Y=Y_AIR))
        g.add_node(PlenumNode("a"))
        g.add_element(
            ChannelElement("feed", "s", "a", length=0.05, diameter=0.03, regime="incompressible")
        )
    g.add_node(PressureBoundary("b", Pt=2.0e5, Tt=600.0, Y=Y_AIR))
    g.add_element(_channel("ch", "a", "b", kind, regime))
    return g


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=60)
    assert r["__success__"], r["__message__"]
    return s, r


def _reservoir_states(m_dot: float) -> tuple[NetworkMixtureState, NetworkMixtureState]:
    a = NetworkMixtureState(P=2.6e5, Pt=2.6e5, Tt=300.0, m_dot=m_dot, Y=Y_AIR)
    b = NetworkMixtureState(P=2.0e5, Pt=2.0e5, Tt=300.0, m_dot=m_dot, Y=Y_AIR)
    return a, b


def _ribbed_on(node: object) -> ChannelElement:
    el = _channel("ch", "a", "b", "ribbed", "incompressible")
    el.area = AREA
    el._upstream_node = node
    el._downstream_node = PressureBoundary("b", Pt=2.0e5, Tt=300.0, Y=Y_AIR)
    return el


def test_the_ribbed_drop_reads_the_duct_static_density_at_a_reservoir() -> None:
    """Ribbed f has no Re dependence, so at a fixed flow the drop is K / rho:
    from a reservoir it is the node-density drop times rho_node / rho_static,
    13% at Mach 0.5 here."""
    m = 0.036
    a, b = _reservoir_states(m)
    reservoir = PressureBoundary("a", Pt=2.6e5, Tt=300.0, Y=Y_AIR)
    moving = MomentumChamberNode("a", area=AREA)  # state_at_rest False
    dP_new = a.Pt - b.Pt - _ribbed_on(reservoir).residuals(a, b)[0][0]
    dP_node = a.Pt - b.Pt - _ribbed_on(moving).residuals(a, b)[0][0]
    s = cb.duct_htc_static_state(300.0, 2.6e5, m / AREA, X_AIR)
    rho_node = cb.complete_state(300.0, 2.6e5, X_AIR).thermo.rho
    assert pytest.approx(0.5, abs=0.03) == s.M
    assert dP_new / dP_node == pytest.approx(rho_node / s.rho, rel=1e-9)
    assert dP_new / dP_node - 1.0 > 0.12


@pytest.mark.parametrize("kind", ["ribbed", "pin_fin"])
@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
@pytest.mark.parametrize("feed", ["pressure", "mass_flow", "plenum"])
def test_every_row_of_a_reservoir_fed_surface_network_matches_fd(
    kind: str, regime: str, feed: str
) -> None:
    """Every row and column, J_ij |x_j| against central differences, at the
    solution and at a perturbed point."""
    g = _reservoir_net(kind, regime, feed)
    s, r = _solve(g)
    x_sol = np.array([float(r[k]) for k in s.unknown_names])
    rng = np.random.default_rng(3)
    for x in (x_sol, x_sol * (1.0 + 0.02 * rng.standard_normal(x_sol.shape))):
        F, J = s._residuals_and_jacobian(x)
        J = J.toarray() if hasattr(J, "toarray") else np.asarray(J)
        FD = fd_jacobian(s, x)
        sc = np.maximum(np.abs(x), 1e-3)
        for i in range(len(F)):
            a, b = J[i] * sc, FD[i] * sc
            den = max(np.max(np.abs(a)), np.max(np.abs(b)))
            assert np.max(np.abs(a - b)) <= 1e-6 * den, (s.unknown_names, i, a, b)


def test_diagnostics_report_the_state_the_drop_reads() -> None:
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=2.6e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("b", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(_channel("ch", "a", "b", "ribbed", "incompressible"))
    _, r = _solve(g)
    s = cb.duct_htc_static_state(300.0, 2.6e5, float(r["ch.m_dot"]) / AREA, X_AIR)
    assert s.M > 0.45
    assert r["ch.mach_in"] == pytest.approx(s.M, rel=1e-6)


@pytest.mark.parametrize("kind", ["ribbed", "pin_fin"])
def test_a_moving_node_keeps_its_own_state(kind: str) -> None:
    """A momentum chamber carries its own static state (Pt > P): the drop
    keeps reading it, bit-identical to before."""
    el = _channel("ch", "c", "b", kind, "incompressible")
    el.area = AREA
    el._upstream_node = MomentumChamberNode("c", area=AREA)
    st = NetworkMixtureState(P=2.2e5, Pt=2.3e5, Tt=500.0, m_dot=0.02, Y=Y_AIR)
    assert el._surface_rest_feed(st, 0.02, AREA) is None
    el._upstream_node = PlenumNode("c")
    assert el._surface_rest_feed(st, 0.02, AREA) is not None
