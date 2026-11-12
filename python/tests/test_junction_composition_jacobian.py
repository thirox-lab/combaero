"""Junction Jacobians on merges of different gases and temperatures (#517).

A junction port node's T and Y are not unknowns: the solver derives them by
mixing and relays any ``{node}.T`` / ``{node}.Y[k]`` Jacobian entry through
that chain. Three residuals depended on composition without saying so:

- ``MomentumChamberNode``'s Pt closure (rho, a, s), 7.6e-5 on a CO2 / air merge;
- ``MultiPortChamberElement``'s rows through the port density, 3.9e-6;
- ``ConstantKTeeElement``'s common-port head (rho, gamma), 7.5e-5.

Measured as J_ij |x_j| against the row's largest such term, at a converged
state perturbed by 1e-3; the central-difference floor is ~2e-7.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
)
from combaero.network.mpce_element import ConstantKTeeElement, MultiPortChamberElement

_D = 0.05
_A = math.pi * (_D / 2.0) ** 2
_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
_CO2 = list(cb.mole_to_mass(list(cb.species.pure_species("CO2"))))


def _merge(kind: str, Y_bra: list[float], Tt_bra: float) -> FlowNetwork:
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("mfb_str", m_dot=0.03, Tt=300.0, Y=_AIR))
    net.add_node(MassFlowBoundary("mfb_bra", m_dot=0.02, Tt=Tt_bra, Y=Y_bra))
    net.add_node(PressureBoundary("pb_out", Pt=1.9e5, Tt=300.0, Y=_AIR))
    for nid in ("mc_com", "mc_str", "mc_bra"):
        net.add_node(MomentumChamberNode(nid, area=_A))
    net.add_element(ChannelElement("ch_str", "mfb_str", "mc_str", length=0.5, diameter=_D))
    net.add_element(ChannelElement("ch_bra", "mfb_bra", "mc_bra", length=0.5, diameter=_D))
    net.add_element(ChannelElement("ch_out", "mc_com", "pb_out", length=0.5, diameter=_D))
    kw = {
        "id": "jct",
        "inlet_nodes": ["mc_str", "mc_bra"],
        "outlet_nodes": ["mc_com"],
        "inlet_angles_deg": [0.0, 90.0],
        "outlet_angles_deg": [0.0],
        "flow_direction": "merge",
    }
    if kind == "mynard":
        net.add_element(MultiPortChamberElement(**kw, strict=False))
    else:
        net.add_element(ConstantKTeeElement(**kw, K_ports={0: 0.4, 1: 1.0}))
    return net


@pytest.mark.parametrize("kind", ["mynard", "constant_k"])
@pytest.mark.parametrize(
    ("Y_bra", "Tt_bra"),
    [(_AIR, 450.0), (_CO2, 300.0), (_CO2, 450.0)],
    ids=["hot-air", "cold-CO2", "hot-CO2"],
)
def test_every_row_matches_central_differences(
    kind: str, Y_bra: list[float], Tt_bra: float
) -> None:
    solver = NetworkSolver(_merge(kind, Y_bra, Tt_bra))
    sol = solver.solve(timeout=60.0, use_jac=False)
    assert sol["__success__"], sol.get("__message__")
    x = np.array([sol[n] for n in solver.unknown_names], dtype=float)
    x = x * (1.0 + 1e-3 * np.random.default_rng(3).standard_normal(x.size))
    _, J = solver._residuals_and_jacobian(x)
    J = J.toarray()
    J_fd = fd_jacobian(solver, x, rel=1e-7, floor=1e-10)
    scaled_err = np.abs(J - J_fd) * np.abs(x)
    row_ref = np.max(np.abs(J_fd) * np.abs(x), axis=1) + 1e-30
    assert float(np.max(scaled_err.max(axis=1) / row_ref)) < 1e-6


def test_momentum_chamber_reports_its_composition_derivative() -> None:
    """The kernel's d(res)/dY_k against a difference on normalised fractions,
    zero for species absent at the node."""
    Y = [0.5 * a + 0.5 * c for a, c in zip(_AIR, _CO2, strict=True)]
    P, T, m, A = 1.9e5, 350.0, 0.05, _A
    r = cb._core.momentum_chamber_residual_and_jacobian(P, 1.95e5, m, T, Y, A)
    for k, y in enumerate(Y):
        if y <= 0.0:
            assert r.d_res_dY[k] == 0.0
            continue
        h = 1e-4
        res = []
        for dy in (h, -h):
            Y_p = list(Y)
            Y_p[k] += dy
            res.append(
                cb._core.momentum_chamber_residual_and_jacobian(P, 1.95e5, m, T, Y_p, A).residual
            )
        assert r.d_res_dY[k] == pytest.approx((res[0] - res[1]) / (2.0 * h), rel=1e-4, abs=1e-6)
