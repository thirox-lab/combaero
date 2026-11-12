"""A discrete pressure loss follows the flow (#550).

``PressureLossElement`` applied ``Pt_to = (1 - xi) Pt_from`` in its declared
direction whatever the flow did. Reversed, a constant fraction RAISED total
pressure by 1/(1 - xi) -- found by the second-law check (#545), in solves
reported as successful -- and a head loss, whose correlation returns 0 for
m_dot <= 0, lost nothing. The loss now runs from the flow's upstream node at
|m_dot|, blended through zero flow as the area change is.

A constant fraction holds a pressure difference at zero flow, so a reversing
flow had to jump through it, which the solver cannot cross. A fraction loss
takes a design flow ``m_ref``: at and above it the loss is exactly xi, below
it the loss scales with the dynamic head and vanishes at zero flow.
"""

from __future__ import annotations

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    FlowNetwork,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
)
from combaero.network.components import FLOW_REVERSAL_M_SCALE
from combaero.network.pressure_loss import (
    ConstantFractionLoss,
    ConstantHeadLoss,
    LinearThetaFractionLoss,
)

Y = list(cb.species.dry_air_mass())


def _net(correlation, P_a: float, P_b: float, mirrored: bool = False) -> FlowNetwork:
    """a -> loss -> p -> partner -> b; ``mirrored`` declares both elements the
    other way, b -> partner -> p -> loss -> a.

    The partner is itself a head loss, which follows the flow: an orifice
    would not do, because a reversed orifice still reads its density at its
    declared inlet (#554) and mirrors only approximately.
    """
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=P_a, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("b", Pt=P_b, Tt=300.0, Y=Y))
    g.add_node(PlenumNode("p"))
    loss, part = (("p", "a"), ("b", "p")) if mirrored else (("a", "p"), ("p", "b"))
    g.add_element(PressureLossElement("loss", *loss, correlation=correlation, area=1e-3))
    g.add_element(
        PressureLossElement(
            "part", *part, correlation=ConstantHeadLoss(zeta=2.0, area=1e-3), area=1e-3
        )
    )
    return g


def _solve(g: FlowNetwork) -> dict:
    r = NetworkSolver(g).solve()
    assert r["__success__"], r["__message__"]
    return r


def test_a_reversed_fraction_loss_loses_total_pressure_along_the_flow() -> None:
    """Flow b -> a through a loss declared a -> p: Pt_a = (1 - xi) Pt_p."""
    r = _solve(_net(ConstantFractionLoss(xi=0.2, m_ref=0.05), 1.0e5, 2.0e5))
    assert r["loss.m_dot"] < 0.0
    assert r["p.Pt"] == pytest.approx(1.0e5 / 0.8, rel=1e-9)
    assert r["loss.xi"] == pytest.approx(0.2)


@pytest.mark.parametrize(
    "correlation",
    [ConstantFractionLoss(xi=0.2, m_ref=0.05), ConstantHeadLoss(zeta=2.0, area=1e-3)],
    ids=["fraction", "head"],
)
def test_a_reversed_loss_is_its_mirror_image(correlation) -> None:
    """Declared against the flow or with it: the same physical loss, so the
    same flow and the same intermediate pressure."""
    against = _solve(_net(correlation, 1.0e5, 2.0e5))
    along = _solve(_net(correlation, 1.0e5, 2.0e5, mirrored=True))
    assert against["loss.m_dot"] == pytest.approx(-along["loss.m_dot"], rel=1e-8)
    assert against["p.Pt"] == pytest.approx(along["p.Pt"], rel=1e-8)
    assert along["p.Pt"] > 1.0e5  # the head loss is no longer zero reversed


def _jacobian_error(g: FlowNetwork, x: np.ndarray) -> float:
    s = NetworkSolver(g)
    s._residuals(x)  # build the solver's maps
    s._build_x0()  # name the unknowns: the FD step reads which are pressures
    _, jac = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x, rel=1e-7, floor=1e-11)
    scale = np.abs(fd * x).max(axis=1, keepdims=True)
    scale[scale == 0.0] = 1.0
    return float((np.abs(jac.toarray() - fd) * np.abs(x) / scale).max())


@pytest.mark.parametrize(
    "correlation",
    [ConstantHeadLoss(zeta=2.0, area=1e-3), ConstantFractionLoss(xi=0.2, m_ref=0.05)],
    ids=["head", "fraction"],
)
@pytest.mark.parametrize("m_over_scale", [-900.0, -40.0, -1.5, 0.3, 40.0, 900.0])
def test_the_jacobian_holds_through_the_reversal(correlation, m_over_scale: float) -> None:
    """Analytic against central differences, row-scaled, at flows in and out
    of the blend band and either side of the design flow (900 scales)."""
    g = _net(correlation, 1.0e5, 1.0e5)
    s = NetworkSolver(g)
    r = s.solve()
    x = np.array(r["__x_solution__"], dtype=float)
    x[s.unknown_names.index("loss.m_dot")] = m_over_scale * FLOW_REVERSAL_M_SCALE
    assert _jacobian_error(g, x) < 1e-6


def _ctx(m: float):
    from types import SimpleNamespace

    return SimpleNamespace(mdot_total=m, theta=0.5)


def _low_flow_factor(ctx: object, m_ref: float | None) -> float:
    """g(|m|/m_ref), read through the C++ fraction loss with xi0 = 1."""
    return ConstantFractionLoss(xi=1.0, m_ref=m_ref)(ctx)[0]


def test_the_low_flow_factor_is_one_at_design_and_c2_there() -> None:
    m_ref = 0.05
    g = [_low_flow_factor(_ctx(m), m_ref) for m in (0.0, 0.025, 0.05, 0.08, -0.08)]
    assert g == pytest.approx([0.0, 1.0 - 0.75**3, 1.0, 1.0, 1.0])
    # zero slope and curvature at r = 1, from below
    h = 1e-4 * m_ref
    f = [_low_flow_factor(_ctx(m_ref - k * h), m_ref) for k in (0, 1, 2)]
    assert (f[0] - f[1]) / h == pytest.approx(0.0, abs=1e-6 / m_ref)
    assert (f[0] - 2 * f[1] + f[2]) / h**2 == pytest.approx(0.0, abs=1e-2 / m_ref**2)
    # head-loss scaling at low flow: 3 r^2
    assert _low_flow_factor(_ctx(1e-4), m_ref) == pytest.approx(3 * (1e-4 / m_ref) ** 2, rel=1e-3)
    assert _low_flow_factor(_ctx(0.01), None) == 1.0


def test_the_design_flow_leaves_the_fraction_exact_above_it() -> None:
    """At and above m_ref the loss is exactly xi, so forward results above the
    design flow do not move; theta-linear fractions scale the same way."""
    assert ConstantFractionLoss(xi=0.2, m_ref=0.05)(_ctx(0.07)) == (0.2, 0.0)
    xi, dxi = LinearThetaFractionLoss(k=0.1, xi0=0.02, m_ref=0.05)(_ctx(0.025))
    g = 1.0 - 0.75**3
    assert (xi, dxi) == pytest.approx(((0.1 * 0.5 + 0.02) * g, 0.1 * g))


def test_the_gui_passes_the_design_flow_through() -> None:
    from gui.backend.graph_builder import _build_discrete_loss_correlation

    c = _build_discrete_loss_correlation("constant_fraction", 0.2, 0, 0, 0, 0, 0.1, m_ref=0.05)
    assert c == ConstantFractionLoss(xi=0.2, m_ref=0.05)
    t = _build_discrete_loss_correlation("linear_theta_fraction", 0, 0.1, 0.02, 0, 0, 0.1)
    assert t.m_ref is None
