"""A solver's first residual is F(x), on networks with a back-edge (#1162).

When a propagation pass reads a state ahead of its order (a wall back-edge,
a recirculation loop) the pass repeats until what it reads settles. The
stopping test compared only the temperatures: on stress ``random`` seed 37
the first evaluation of a fresh solver stopped after two passes, T settled
to 1e-16 while a lagged composition still moved by 0.025, and combustor n4
sat at 610 K where the fixed point is 1992 K (F off by 0.24). The contraction
there is ~0.55 per pass and needs ~82 passes, past the old cap of 50.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
import test_solver_robustness_481 as tsr
from _fd_jacobian import fd_jacobian, scaled_error

import combaero as cb
from combaero.network import (
    CombustorNode,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    PlenumNode,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
from networks import FAMILIES  # noqa: E402

# The final x of the failing random-37 solve (#1158, tmp/1158/r37_x.npy),
# where a fresh solver's first evaluation was 0.24 off its second.
X_RANDOM_37 = np.array(
    [
        float.fromhex(h)
        for h in (
            "0x1.e5f6dc3381ecep+16",
            "0x1.e5f6dc3381e5ap+16",
            "0x1.0ac9204845be9p+21",
            "0x1.0ac9204845e9bp+21",
            "0x1.0b4df700308a8p+21",
            "0x1.0b4df700308c4p+21",
            "0x1.e5f7e1c51060dp+16",
            "0x1.e5f7e1c5105e9p+16",
            "0x1.9dc4f2afcd85dp+20",
            "0x1.9dc4f2afcd860p+20",
            "0x1.0b51ca86ff012p+21",
            "0x1.0b51ca86fefffp+21",
            "0x1.0af9a6408c1bbp+21",
            "0x1.0af9a6408c1bdp+21",
            "0x1.8ea18766ea716p+16",
            "0x1.8ea18766ea714p+16",
            "0x1.0624dd2f1a8f0p-10",
            "0x1.e75aa33335f98p-10",
            "0x1.cf1d12b7efea0p-5",
            "-0x1.d8c2123a9903fp-3",
            "-0x1.f9507d6e2db9cp-4",
            "-0x1.89de3f53062c7p-5",
            "-0x1.6f60cab1c966dp-2",
            "0x1.4b3e3f1e3ec57p-5",
            "-0x1.855c4b9ce7da6p-7",
            "0x1.d3ce586e099abp-6",
            "-0x1.4abab253c1ac2p+1",
            "-0x1.1982d1c762971p-4",
            "0x1.572e6592d8d43p+1",
            "0x0.0p+0",
            "0x1.7fa1c6d7ba674p-3",
            "0x1.8c9db038aa008p-2",
        )
    ]
)


def _fresh(net: FlowNetwork) -> NetworkSolver:
    s = NetworkSolver(net)
    s._build_x0()
    return s


def _first_and_second(net: FlowNetwork, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    s = _fresh(net)
    assert len(s.unknown_names) == x.size
    return np.asarray(s._residuals(x), float), np.asarray(s._residuals(x), float)


def _recirculating_combustor() -> FlowNetwork:
    """A combustor whose exit plenum feeds part of its flow back (c -> p ->
    c): at a trial point with both loop flows positive, c reads p ahead of
    the order; 0.18 of the 0.2 kg/s goes round again (0.9 per pass)."""
    y_ch4 = list(cb.mole_to_mass(cb.species.pure_species("CH4")))
    g = FlowNetwork()
    g.add_node(tsr._pb("air", 2.0e5, 600.0))
    g.add_node(MassFlowBoundary("fuel", m_dot=0.002, Tt=300.0, Y=y_ch4))
    g.add_node(CombustorNode("c", method="complete"))
    g.add_node(PlenumNode("p"))
    g.add_node(tsr._pb("out", 1.0e5, 300.0))
    g.add_element(tsr._orifice("ea", "air", "c", 0.02))
    g.add_element(tsr._orifice("ef", "fuel", "c", 0.005))
    g.add_element(tsr._orifice("cp", "c", "p", 0.03))
    g.add_element(tsr._orifice("pc", "p", "c", 0.03))
    g.add_element(tsr._orifice("po", "p", "out", 0.02))
    return g


def _recirculating_x() -> np.ndarray:
    s = _fresh(_recirculating_combustor())
    x = np.array(s._build_x0(), float)
    flows = {"ea": 0.05, "ef": 0.002, "cp": 0.2, "pc": 0.18, "po": 0.05}
    for eid, m in flows.items():
        x[s.unknown_names.index(f"{eid}.m_dot")] = m
    return x


def test_random_37_first_evaluation_is_the_fixed_point() -> None:
    """Was 0.24 (T-only test) and 7.8e-9 (the 50-pass cap)."""
    net, _ = FAMILIES["random"](37)
    s = _fresh(net)
    F1 = np.asarray(s._residuals(X_RANDOM_37), float)
    # Burned: the fixed point, not the 610 K of the two-pass stop.
    assert s._derived_states["n4"][0] == pytest.approx(1992.29, abs=0.01)
    F2 = np.asarray(s._residuals(X_RANDOM_37), float)
    assert np.max(np.abs(F1 - F2)) < 1e-12


def test_two_wall_series_first_evaluation_is_the_fixed_point() -> None:
    """The wall back-edge of #481, at a point off its root."""
    s = NetworkSolver(tsr._series_walls())
    x = np.array(s.solve()["__x_solution__"])
    xp = x * (1.0 + 0.01 * np.sign(np.arange(len(x)) % 2 - 0.5))
    F1, F2 = _first_and_second(tsr._series_walls(), xp)
    assert np.max(np.abs(F1 - F2)) < 1e-12


def test_recirculating_combustor_first_evaluation_is_the_fixed_point() -> None:
    """Was 3.5e-8 at the 50-pass cap."""
    F1, F2 = _first_and_second(_recirculating_combustor(), _recirculating_x())
    assert np.max(np.abs(F1 - F2)) < 1e-12


def test_random_37_jacobian_matches_fd_at_the_fixed_point() -> None:
    """The relay iterates the borrowed sensitivities with the states, so J
    is the derivative of the settled fixed point (not of lagged constants):
    9.3e-8 scaled against a central FD (2e-5 before, FD of a residual that
    depended on the previous evaluation)."""
    net, _ = FAMILIES["random"](37)
    net.resolve_all_topology()
    s = _fresh(net)
    _, J = s._residuals_and_jacobian(X_RANDOM_37, compute_jacobian=True)
    err = scaled_error(J, fd_jacobian(s, X_RANDOM_37), X_RANDOM_37)
    assert float(np.max(err)) < 1e-6


def test_an_unsettled_propagation_warns_once() -> None:
    """At the cap the residual depends on the previous evaluation: reported
    as a RuntimeWarning, once per solver, like the penalty residual."""
    s = _fresh(_recirculating_combustor())
    s._MAX_PROPAGATION_PASSES = 3
    x = _recirculating_x()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        s._residuals(x)
        s._residuals(x * (1.0 + 1e-3))
    hits = [m for m in w if "did not settle" in str(m.message)]
    assert len(hits) == 1
    assert issubclass(hits[0].category, RuntimeWarning)
