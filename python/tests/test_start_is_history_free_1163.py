"""The cold start is a function of the user's inputs alone (#1163).

A FuelBoundary's ``m_dot`` holds the solver's current iterate of its unknown
(#892). The start estimates read it, so a start built after any evaluation
took the fuel flow that evaluation left behind: on the Brayton cycle a 6e-17
change left by a failed stage decided whether the next stage converged.
"""

from __future__ import annotations

import contextlib
import sys
import warnings
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest

from combaero.network import FlowNetwork, NetworkSolver

sys.path.insert(0, str(Path(__file__).parent))
from test_fuel_boundary_892 import _net, _routed  # noqa: E402
from test_turbomachine_elements_893 import _brayton  # noqa: E402

BUILDS: dict[str, Callable[[], FlowNetwork]] = {
    "brayton": _brayton,
    "direct": lambda: _routed("direct"),
    "premix": lambda: _routed("premix"),
    "line": lambda: _routed("line"),
    "pressure_driven": lambda: _net(1500.0),
    "rich": lambda: _net(2000.0, branch="rich"),
    "with_guess": lambda: _net(1500.0, m_dot_guess=0.01),
}


def _quiet(f: Callable[[], object]) -> object:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return f()


@pytest.mark.parametrize("build", list(BUILDS))
def test_the_start_is_the_same_after_an_evaluation_anywhere(build: str) -> None:
    """``_build_x0`` gives the fresh start, bit for bit, after ``_residuals``
    at any x and after a whole solve."""
    s = NetworkSolver(BUILDS[build]())
    x0 = np.array(_quiet(s._build_x0))
    rng = np.random.default_rng(1163)
    for scale in (1e-12, 0.05, 0.5):
        x = x0 * (1.0 + scale * rng.standard_normal(x0.size))
        # An evaluation that fails part-way still leaves its iterate.
        with contextlib.suppress(ValueError, RuntimeError, ArithmeticError):
            _quiet(lambda x=x: s._residuals(x))
        assert np.array_equal(np.array(_quiet(s._build_x0)), x0), scale
    _quiet(s.solve)
    assert np.array_equal(np.array(_quiet(s._build_x0)), x0)


def test_an_evaluation_between_stages_does_not_change_the_brayton_solve(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """One more evaluation at a nested (stage) solve's returned x moves the
    fuel iterate by ~1e-16; the cycle's solve went from converged to |F| 52.6
    on it. It now returns the same root either way."""
    ref = _quiet(NetworkSolver(_brayton()).solve)
    assert ref["__success__"], ref["__message__"]

    orig = NetworkSolver.solve
    depth = [0]
    evaluated = [0]

    def nested_then_evaluate(self: NetworkSolver, *a: object, **k: object) -> dict:
        depth[0] += 1
        try:
            r = orig(self, *a, **k)
        finally:
            depth[0] -= 1
        x = r.get("__x_solution__")
        if depth[0] >= 1 and x is not None:
            # A point next to the returned one, as a stage's last probe is.
            self._residuals(np.asarray(x, dtype=float) * (1.0 + 1e-15))
            evaluated[0] += 1
        return r

    monkeypatch.setattr(NetworkSolver, "solve", nested_then_evaluate)
    r = _quiet(NetworkSolver(_brayton()).solve)
    assert evaluated[0] >= 1
    assert r["__success__"], r["__message__"]
    assert np.array_equal(np.asarray(r["__x_solution__"]), np.asarray(ref["__x_solution__"]))


def test_a_leftover_fuel_flow_is_not_the_start() -> None:
    """Whatever a FuelBoundary's m_dot holds, the start is the fresh one; a
    given ``m_dot_guess`` is what it starts from. The given flow every start
    estimate reads (the ejector's closed-form seed too) ignores the iterate."""
    for build in ("brayton", "with_guess"):
        net = BUILDS[build]()
        s = NetworkSolver(net)
        x0 = np.array(_quiet(s._build_x0))
        for left in (1e-6, 0.0027, 0.05, 3.0):
            net.nodes["fuel"].m_dot = left
            assert np.array_equal(np.array(_quiet(s._build_x0)), x0), (build, left)
            guess = net.nodes["fuel"].m_dot_guess
            assert net.nodes["fuel"]._given_m_dot() == (0.0 if guess is None else guess)


def test_a_homotopy_rung_that_raises_restores_the_boundary_flows(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The homotopy ladder scales every MassFlowBoundary's m_dot; a rung that
    raised left them at a fifth of the user's flow."""
    net = _routed("direct")
    s = NetworkSolver(net)
    orig = NetworkSolver._solve_impl
    calls = [0]

    def raising(self: NetworkSolver, *a: object, **k: object) -> dict:
        calls[0] += 1
        if calls[0] > 1:  # the first rung, inside the ladder
            raise RuntimeError("rung failed")
        return orig(self, *a, **k)

    monkeypatch.setattr(NetworkSolver, "_solve_impl", raising)
    with pytest.raises(RuntimeError, match="rung failed"):
        _quiet(lambda: s.solve(init_strategy="homotopy"))
    assert calls[0] == 2
    assert net.nodes["air"].m_dot == 0.3
