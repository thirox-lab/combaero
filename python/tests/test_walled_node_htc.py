"""A thermal wall on a node: the combustor's and momentum chamber's htc at the
node's static state, and an exact Jacobian for both.

Three defects, each found on the same walled-node network:
- the wall relay put the node's own-flow htc derivative on the node's P
  column, not on its inflows (the node's first unknown is P);
- a combustor handed its stagnation T to the correlation as static, so T_aw
  sat above Tt (+46 K at Mach 0.49), the duct defect of #636;
- the combustor's mixing Jacobian dropped d(Q / m)/dm: with heat on the node,
  dT/dm was 0 where the difference gives 1658 K per kg/s.
"""

from __future__ import annotations

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    CombustorNode,
    ConvectiveSurface,
    EnergyBoundary,
    FlowNetwork,
    LosslessConnectionElement,
    MomentumChamberNode,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)
from combaero.network.components import SmoothModel


def _node(kind: str, area: float, surface: ConvectiveSurface | None) -> object:
    if kind == "combustor":
        return CombustorNode("hn", method="complete", area=area, Dh=0.02, surface=surface)
    if kind == "equilibrium":
        return CombustorNode("hn", method="equilibrium", area=area, Dh=0.02, surface=surface)
    return MomentumChamberNode("hn", area=area, Dh=0.02, surface=surface)


def _network(kind: str, area: float, walled: bool = True, fuel: bool = False) -> FlowNetwork:
    g = FlowNetwork()
    Y = None
    if fuel:
        X = list(cb.species.dry_air())
        X[cb.species_index_from_name("CH4")] = 0.02
        Y = list(cb.mole_to_mass(cb.normalize_fractions(X)))
    g.add_node(PressureBoundary("hi", Pt=1.6e5, Tt=1500.0 if not fuel else 700.0, Y=Y))
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
    surface = ConvectiveSurface(model=SmoothModel(), area=0.05) if walled else None
    node = _node(kind, area, surface)
    g.add_node(node)
    for eid, a, b in (("hin", "hi", "hn"), ("hout", "hn", "ho")):
        g.add_element(ChannelElement(eid, a, b, length=0.3, diameter=0.02, regime="compressible"))
    if walled:
        g.add_node(PressureBoundary("ci", Pt=1.6e5, Tt=400.0))
        g.add_node(PressureBoundary("co", Pt=1.0e5, Tt=300.0))
        g.add_node(PlenumNode("cp"))
        g.add_element(LosslessConnectionElement("co1", "ci", "cp"))
        g.add_element(
            ChannelElement(
                "cch",
                "cp",
                "co",
                length=0.5,
                diameter=0.015,
                regime="compressible",
                surface=ConvectiveSurface(model=SmoothModel()),
            )
        )
        g.add_wall(ThermalWall("w", "hn", "cch", layers=[WallLayer(0.001, 20.0)]))
    else:
        node.add_energy_boundary(EnergyBoundary("heat", Q=-5.0e3))
    return g


def _jacobian_error(s: NetworkSolver, x: np.ndarray) -> float:
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    # RESIDUAL-ROW-SCALING: compare J_ij |x_j|.
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    return float(np.max(np.abs(Js - fds) / scale))


@pytest.mark.parametrize("kind", ["combustor", "chamber"])
@pytest.mark.parametrize("area", [3.14e-4, 1e-3])
def test_a_walled_node_has_an_exact_jacobian(kind: str, area: float) -> None:
    s = NetworkSolver(_network(kind, area))
    r = s.solve()
    assert r["__success__"], r["__message__"]
    # Before: 3e4 at the node's P column, the inflow column missing.
    assert _jacobian_error(s, np.array(r["__x_solution__"])) < 1e-6


@pytest.mark.parametrize("kind", ["combustor", "chamber"])
def test_a_walled_nodes_t_aw_lies_below_its_stagnation_temperature(kind: str) -> None:
    s = NetworkSolver(_network(kind, 3.14e-4))
    r = s.solve()
    assert r["__success__"], r["__message__"]
    x = np.array(r["__x_solution__"])
    node = s.network.nodes["hn"]
    d = node.diagnostics(s._get_node_state(node, x))
    # Mach about 0.5: the recovery factor below one puts T_aw a few K under
    # Tt. The combustor read 1430.6 K against Tt = 1384.2 K before.
    assert d["T_aw"] < d["Tt"]
    assert d["T_aw"] > d["Tt"] - 0.05 * d["Tt"]


@pytest.mark.parametrize("kind", ["combustor", "equilibrium"])
@pytest.mark.parametrize("fuel", [False, True])
def test_a_heated_combustor_has_an_exact_jacobian(kind: str, fuel: bool) -> None:
    """Heat on a combustor (an EnergyBoundary, no wall) is spread over its
    throughput: T moves with the inflow through Q / m."""
    s = NetworkSolver(_network(kind, 1e-3, walled=False, fuel=fuel))
    r = s.solve()
    assert r["__success__"], r["__message__"]
    assert _jacobian_error(s, np.array(r["__x_solution__"])) < 1e-5
