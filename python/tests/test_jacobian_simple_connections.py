"""Analytic Jacobian vs central differences for the elements no other test
checked (pre-0.7.0 sweep): EffectiveAreaConnectionElement,
DiameterDischargeCoefficientConnectionElement (Cd and zeta forms) and
VortexElement (free vortex and shaft speed), each between two channels, fed
by a pressure or a mass-flow boundary.

Metric (RESIDUAL-ROW-SCALING): J_ij |x_j| relative to each row's largest.
Measured at the solution: 2e-8 to 1e-7 at a 1e-7 relative step.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    DiameterDischargeCoefficientConnectionElement,
    EffectiveAreaConnectionElement,
    FlowNetwork,
    MassFlowBoundary,
    NetworkElement,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    VortexElement,
)

Y_AIR = cb.mole_to_mass(cb.species.dry_air())

ELEMENTS: dict[str, Callable[[str, str], NetworkElement]] = {
    "effective-area": lambda a, b: EffectiveAreaConnectionElement("x", a, b, 8e-4),
    "diameter-Cd": lambda a, b: DiameterDischargeCoefficientConnectionElement(
        "x", a, b, diameter=0.03, Cd=0.7
    ),
    "diameter-zeta": lambda a, b: DiameterDischargeCoefficientConnectionElement(
        "x", a, b, diameter=0.03, zeta=2.5
    ),
    "vortex-free": lambda a, b: VortexElement("x", a, b, r_out=0.1, r_in=0.05, v_theta_in=30.0),
    "vortex-shaft": lambda a, b: VortexElement(
        "x", a, b, r_out=0.1, r_in=0.05, shaft_speed_rpm=3000.0
    ),
}


def _network(make: Callable[[str, str], NetworkElement], mass_inlet: bool) -> FlowNetwork:
    g = FlowNetwork()
    if mass_inlet:
        g.add_node(MassFlowBoundary("in", m_dot=0.3, Tt=500.0, Y=Y_AIR))
    else:
        g.add_node(PressureBoundary("in", Pt=2.0e5, Tt=500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("out", Pt=1.6e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PlenumNode("a"))
    g.add_node(PlenumNode("b"))
    g.add_element(ChannelElement("c1", "in", "a", length=0.5, diameter=0.05, roughness=1e-5))
    g.add_element(make("a", "b"))
    g.add_element(ChannelElement("c2", "b", "out", length=0.5, diameter=0.06, roughness=1e-5))
    return g


def _jac_err(s: NetworkSolver, x: np.ndarray) -> float:
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x, rel=1e-7, floor=1e-12)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    a, b = J.toarray() * cols, fd * cols
    scale = np.maximum(np.abs(b).max(axis=1, keepdims=True), 1e-300)
    return float((np.abs(a - b) / scale).max())


@pytest.mark.parametrize("mass_inlet", [False, True], ids=["pressure-fed", "mass-fed"])
@pytest.mark.parametrize("kind", list(ELEMENTS))
def test_the_jacobian_matches_central_differences(kind: str, mass_inlet: bool) -> None:
    s = NetworkSolver(_network(ELEMENTS[kind], mass_inlet))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r.get("__message__")
    assert _jac_err(s, np.asarray(r["__x_solution__"])) < 1e-5
