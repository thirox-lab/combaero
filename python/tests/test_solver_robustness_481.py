"""Robustness core of the #481 audit: conservation under reversed flow, a
residual that depends on x alone, Jacobians that hold through reversal, and
a solver that survives non-physical probes.
"""

from __future__ import annotations

import io
from contextlib import redirect_stdout

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero import _solver_tools
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)

with redirect_stdout(io.StringIO()):
    import test_impingement_array as array_cases

from _closure_check import closure

X_AIR = cb.species.dry_air()
Y_AIR = cb.mole_to_mass(X_AIR)


def _pb(name: str, Pt: float, Tt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y_AIR
    return b


def _orifice(eid: str, a: str, b: str, d: float = 0.01, in_pipe: bool = True) -> OrificeElement:
    # in_pipe=False where the bore meets its pipe: a lumped Cd*A
    # restriction, not the no-orifice of a bore equal to its pipe (#788).
    return OrificeElement(eid, a, b, Cd=0.7, diameter=d, correlation="fixed", in_pipe=in_pipe)


def _jac_vs_fd(s: NetworkSolver, x: np.ndarray) -> float:
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x, rel=1e-5, floor=1e-7)
    # RESIDUAL-ROW-SCALING (NetworkSolver._build_residual_scales). In scaled
    # variables (J_ij |x_j|): rows mix kg/s and Pa columns, and raw
    # row scaling hides a wrong entry in the small-unit columns.
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    return float(np.max(np.abs(Js - fds) / scale))


# --- reversed flow: the plenum takes the stream that actually feeds it ------


def test_a_reversed_line_feeds_the_plenum_from_its_real_source() -> None:
    """Both orifices declared A -> p -> B, the pressures drive B -> p -> A:
    the plenum is fed at B's 600 K, not A's 400 K, and energy closes."""
    g = FlowNetwork()
    g.add_node(_pb("A", 1.02e5, 400.0))
    g.add_node(_pb("B", 1.03e5, 600.0))
    g.add_node(PlenumNode("p"))
    g.add_element(_orifice("e1", "A", "p"))
    g.add_element(_orifice("e2", "p", "B"))
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"]
    assert r["e1.m_dot"] < 0.0 and r["e2.m_dot"] < 0.0
    assert s._derived_states["p"][0] == pytest.approx(600.0, abs=1e-6)
    c = closure(s, r)
    assert abs(c["energy"]) < c["energy_tol"]


def test_a_reversed_branch_brings_its_heat_in() -> None:
    """A third line declared n -> c delivers from the hot c INTO n (the
    audit's case 1a): n mixes it in and the network balance closes."""
    g = FlowNetwork()
    g.add_node(_pb("a", 2.0e5, 300.0))
    g.add_node(_pb("c", 1.6e5, 800.0))
    g.add_node(_pb("out", 1.0e5, 300.0))
    g.add_node(PlenumNode("n"))
    g.add_element(_orifice("e1", "a", "n"))
    g.add_element(_orifice("e2", "n", "out", 0.014))
    g.add_element(_orifice("e3", "n", "c"))
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"]
    m1, m3 = r["e1.m_dot"], -r["e3.m_dot"]
    assert m3 > 0.0
    h_mix = (m1 * cb.h_mass(300.0, X_AIR) + m3 * cb.h_mass(800.0, X_AIR)) / (m1 + m3)
    assert cb.h_mass(s._derived_states["n"][0], X_AIR) == pytest.approx(h_mix, rel=1e-9)
    c = closure(s, r)
    assert abs(c["energy"]) < c["energy_tol"]


def test_the_relay_carries_the_reversed_stream() -> None:
    """The global Jacobian through a reversed, non-isothermal plenum."""
    g = FlowNetwork()
    g.add_node(_pb("A", 1.02e5, 400.0))
    g.add_node(_pb("B", 1.03e5, 600.0))
    g.add_node(PlenumNode("p"))
    g.add_element(_orifice("e1", "A", "p", in_pipe=False))
    g.add_element(ChannelElement("e2", "p", "B", length=0.5, diameter=0.01))
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"] and r["e2.m_dot"] < 0.0
    assert _jac_vs_fd(s, np.array(r["__x_solution__"])) < 1e-4


def test_a_reversed_channel_jacobian_has_the_right_sign() -> None:
    """dP = f m|m| ... is odd in m, so its slope is even: the C++ term used
    signed m and gave reversed flow the wrong sign (+9.7e4 vs -8.0e4)."""
    g = FlowNetwork()
    g.add_node(_pb("a", 1.0e5, 300.0))
    g.add_node(_pb("b", 1.01e5, 300.0))
    g.add_node(PlenumNode("p"))
    g.add_element(ChannelElement("c", "a", "p", length=1.0, diameter=0.02, roughness=1e-5))
    g.add_element(_orifice("o", "p", "b", 0.02, in_pipe=False))
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"] and r["c.m_dot"] < 0.0
    assert _jac_vs_fd(s, np.array(r["__x_solution__"])) < 1e-4


# --- a residual that depends on x alone --------------------------------------


def _series_walls() -> FlowNetwork:
    """Two channels in series coupled by a wall: the downstream one feeds
    heat back upstream -- a back-edge, read ahead of the propagation order."""
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("a", m_dot=0.05, Tt=900.0, Y=Y_AIR))
    g.add_node(PlenumNode("n2"))
    g.add_node(PlenumNode("n3"))
    g.add_node(_pb("out", 1.0e5, 300.0))
    area = 3.14159 * 0.03
    g.add_element(
        ChannelElement(
            "e1", "a", "n2", length=1.0, diameter=0.03, surface=ConvectiveSurface(area=area)
        )
    )
    g.add_element(
        ChannelElement(
            "e2", "n2", "n3", length=1.0, diameter=0.03, surface=ConvectiveSurface(area=area)
        )
    )
    g.add_element(_orifice("ex", "n3", "out", 0.02))
    g.add_node(MassFlowBoundary("cool", m_dot=0.05, Tt=300.0, Y=Y_AIR))
    g.add_node(_pb("cool_out", 1.0e5, 300.0))
    g.add_node(PlenumNode("cp"))
    g.add_element(
        ChannelElement(
            "cc", "cool", "cp", length=1.0, diameter=0.03, surface=ConvectiveSurface(area=area)
        )
    )
    g.add_element(_orifice("cx", "cp", "cool_out", 0.02))
    g.add_wall(
        ThermalWall(id="w1", element_a="e2", element_b="e1", layers=[WallLayer(0.002, 20.0)])
    )
    g.add_wall(
        ThermalWall(id="w2", element_a="e1", element_b="cc", layers=[WallLayer(0.002, 20.0)])
    )
    return g


def test_the_residual_does_not_depend_on_the_previous_evaluation() -> None:
    s = NetworkSolver(_series_walls())
    r = s.solve()
    assert r["__success__"]
    x = np.array(r["__x_solution__"])
    xp = x * (1.0 + 0.01 * np.sign(np.arange(len(x)) % 2 - 0.5))
    f1, _ = s._residuals_and_jacobian(xp)
    s._residuals_and_jacobian(x)
    f2, _ = s._residuals_and_jacobian(xp)
    f3, _ = s._residuals_and_jacobian(xp)
    assert np.max(np.abs(f1 - f2)) <= 1e-9 * max(np.max(np.abs(f1)), 1.0)
    assert np.max(np.abs(f2 - f3)) <= 1e-9 * max(np.max(np.abs(f1)), 1.0)


def test_a_temperature_dependent_wall_does_not_remember_trial_points() -> None:
    """k(T) used to be lagged from the previous call -- whatever trial point
    the solver had probed. Now the same inputs give the same heat."""
    wall = ThermalWall(
        id="w", element_a="a", element_b="b", layers=[WallLayer(0.003, 15.0, material="inconel718")]
    )
    q1 = wall.compute_coupling(500.0, 1500.0, 0.1, 800.0, 600.0, 0.1).Q
    wall.compute_coupling(50.0, 400.0, 0.1, 80.0, 300.0, 0.1)
    q2 = wall.compute_coupling(500.0, 1500.0, 0.1, 800.0, 600.0, 0.1).Q
    assert q2 == pytest.approx(q1, rel=1e-10)


def test_wall_coupling_floors_a_non_positive_htc_smoothly() -> None:
    """h <= 0 (a correlation past its range at an iterate) gives a finite,
    positive coupling, and dQ/dh is the slope of the floored value."""
    for h in (-50.0, -1e-3, 0.0, 5e-3, 2e-2, 10.0):
        r = _solver_tools.wall_coupling_and_jacobian_multilayer(
            h, 1500.0, 200.0, 500.0, [1e-4], 0.1, 0.0
        )
        assert np.isfinite(r.Q) and r.Q > 0.0
        step = max(abs(h) * 1e-6, 1e-8)
        qp = _solver_tools.wall_coupling_and_jacobian_multilayer(
            h + step, 1500.0, 200.0, 500.0, [1e-4], 0.1, 0.0
        ).Q
        qm = _solver_tools.wall_coupling_and_jacobian_multilayer(
            h - step, 1500.0, 200.0, 500.0, [1e-4], 0.1, 0.0
        ).Q
        assert r.dQ_dh_a == pytest.approx((qp - qm) / (2.0 * step), rel=1e-5, abs=1e-12)
    assert cb.WALL_HTC_KNEE == 1e-2


def test_with_coupling_off_the_walls_report_no_heat() -> None:
    g = _series_walls()
    g.thermal_coupling_enabled = False
    r = NetworkSolver(g).solve()
    assert r["__success__"]
    assert not any(k.startswith(("w1.", "w2.")) or k.endswith(".Q_wall_out") for k in r)


# --- the solver survives non-physical probes ---------------------------------


_HOT_ARRAYS = [(5, 1.01e5), (5, 1.02e5), (8, 1.05e5), (8, 1.1e5)]


def _hot_array(n_rows: int, p_hot: float, regime: str) -> FlowNetwork:
    g = array_cases._net(array_cases._array(n_rows), hot=True)
    g.nodes["hot_in"].Pt = p_hot
    for item in list(g.elements.values()) + list(g.nodes.values()):
        if hasattr(item, "regime"):
            item.regime = regime
    return g


@pytest.mark.parametrize(("n_rows", "p_hot"), _HOT_ARRAYS)
def test_a_hot_impingement_array_converges_cold(n_rows: int, p_hot: float) -> None:
    """Without any retry. These four failed cold before #481: the penalty
    accepted steps into the unphysical region, and reversed jets mixed as
    negative mass parked the chambers at 1500-3100 K. Pinned incompressible:
    its pressure-drop rows are linear in the chamber pressures, which is what
    this cold guarantee rests on (see the compressible test below)."""
    g = _hot_array(n_rows, p_hot, "incompressible")
    s = NetworkSolver(g)
    r = s.solve(auto_retry=False)
    assert r["__success__"], r["__message__"]
    assert all(s._derived_states[f"ia__c{i}"][0] < 1200.0 for i in range(1, n_rows + 1))


@pytest.mark.parametrize(("n_rows", "p_hot"), _HOT_ARRAYS)
def test_a_hot_impingement_array_converges_compressible(n_rows: int, p_hot: float) -> None:
    """The same arrays compressible, on the default solve path (#605). Their
    cold attempt alone is not guaranteed: the crossflow segments drive their
    friction with a few Pa (the station momentum takes the rest of the drop),
    where the Fanno mass-flow row goes as sqrt(drive), and the wall heat
    changes node as a jet reverses. 5 rows at 1.01 bar stalls cold at |F| ~
    0.13, and which cases do is not ordered by the drive. The incompressible
    warm-start retry (#603) is the designed recovery; it must reach the
    physical root, the chambers well below the hot gas."""
    g = _hot_array(n_rows, p_hot, "compressible")
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"], r["__message__"]
    assert all(s._derived_states[f"ia__c{i}"][0] < 1200.0 for i in range(1, n_rows + 1))


def test_a_probe_the_model_refuses_is_rejected_not_accepted(monkeypatch) -> None:
    """A residual that throws beyond a flow: the penalty must read as WORSE
    than where the solver stands, so it steps back and still converges."""
    g = FlowNetwork()
    g.add_node(_pb("a", 2.0e5, 300.0))
    g.add_node(_pb("b", 1.0e5, 300.0))
    g.add_element(_orifice("o", "a", "b", 0.02))
    s = NetworkSolver(g)
    original = s._residuals_and_jacobian
    limit = 0.2

    def refusing(x, **kw):
        if abs(x[s.unknown_names.index("o.m_dot")]) > limit:
            raise ValueError("refused")
        return original(x, **kw)

    s._residuals_and_jacobian = refusing
    r = s.solve()
    assert r["__success__"]
    assert 0.0 < r["o.m_dot"] < limit


def test_a_wall_on_a_reversed_channel_relays_the_right_sign() -> None:
    """Surface correlations take |m_dot|, so their dh/dm is in |m|: for a
    reversed element the relay must flip it (+0.096 vs -0.096 before)."""
    g = FlowNetwork()
    g.add_node(_pb("A", 1.00e5, 400.0))
    g.add_node(_pb("B", 1.02e5, 900.0))
    g.add_node(PlenumNode("p"))
    area = 3.14159 * 0.02
    surface = ConvectiveSurface(area=area)
    g.add_element(ChannelElement("hot", "A", "p", length=1.0, diameter=0.02, surface=surface))
    g.add_element(_orifice("o", "p", "B", 0.02, in_pipe=False))
    g.add_node(MassFlowBoundary("cin", m_dot=0.01, Tt=300.0, Y=Y_AIR))
    g.add_node(PlenumNode("cp"))
    g.add_node(_pb("cout", 1.0e5, 300.0))
    g.add_element(
        ChannelElement(
            "cold", "cin", "cp", length=1.0, diameter=0.02, surface=ConvectiveSurface(area=area)
        )
    )
    g.add_element(_orifice("co", "cp", "cout", 0.02, in_pipe=False))
    g.add_wall(
        ThermalWall(id="w", element_a="hot", element_b="cold", layers=[WallLayer(0.002, 20.0)])
    )
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"] and r["hot.m_dot"] < 0.0 and r["w.Q"] > 0.0
    assert _jac_vs_fd(s, np.array(r["__x_solution__"])) < 1e-4
