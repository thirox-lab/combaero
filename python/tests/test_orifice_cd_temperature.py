"""A Reynolds-dependent orifice Cd moves with the upstream temperature.

Re = 4 |m| / (pi D mu(T)): ISO metering (Re_D) and discharge-hole (Re_hole)
correlations alike. The row's Jacobian had the flow column (#639) but not
the temperature one, so an orifice behind a heated plenum was 1e-3 off in
the column that moves the plenum's T (here the feed orifice's flow, through
Q / m).
"""

from __future__ import annotations

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

from combaero.network import (
    EnergyBoundary,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)


def _network(corr: str, regime: str) -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=2.0e5, Tt=600.0))
    g.add_node(PressureBoundary("b", Pt=1.0e5, Tt=300.0))
    p = PlenumNode("p")
    p.add_energy_boundary(EnergyBoundary("q", Q=2.0e4))
    g.add_node(p)
    g.add_element(
        OrificeElement("o1", "a", "p", Cd=0.7, diameter=0.02, correlation="fixed", regime=regime)
    )
    g.add_element(OrificeElement("o2", "p", "b", diameter=0.02, correlation=corr, regime=regime))
    return g


@pytest.mark.parametrize("corr", ["ReaderHarrisGallagher", "IdelchikThick", "McGreehanSchotsch"])
@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
def test_a_reynolds_dependent_cd_has_its_temperature_column(corr: str, regime: str) -> None:
    s = NetworkSolver(_network(corr, regime))
    r = s.solve()
    assert r["__success__"], r["__message__"]
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    # RESIDUAL-ROW-SCALING: compare J_ij |x_j|. Before: 6e-4 to 3.5e-3.
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    assert float(np.max(np.abs(Js - fds) / scale)) < 1e-6
