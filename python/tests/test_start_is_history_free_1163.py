"""The cold start is a function of the user's inputs alone (#1163).

A FuelBoundary's ``m_dot`` holds the solver's current iterate of its unknown
(#892). The start estimates read it, so a start built after any evaluation
took the fuel flow that evaluation left behind: on the Brayton cycle a 6e-17
change left by a failed stage decided whether the next stage converged.
"""

from __future__ import annotations

import contextlib
import math
import sys
import warnings
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest

from combaero.network import FlowNetwork, NetworkSolver

sys.path.insert(0, str(Path(__file__).parent))
from test_fuel_boundary_892 import CH4, _net, _routed  # noqa: E402
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


@pytest.mark.parametrize(
    ("build", "air", "fuel_el", "out"),
    [("brayton", "comp", "inj", "turb"), ("pressure_driven", "inlet", "inj", "exit")],
)
def test_the_fuel_fixed_stage_starts_with_the_fuel_it_holds(
    build: str, air: str, fuel_el: str, out: str
) -> None:
    """The fuel-fixed stage holds each controller at its fresh start's flow,
    and its own start seeds the flows past the combustor with that flow, not
    with none: their start closes the combustor's mass balance at the held
    fuel. Held only inside the stage."""
    s = NetworkSolver(BUILDS[build]())
    fuel = s.network.nodes["fuel"]
    orig = s._build_x0
    held: list[tuple[float, float]] = []

    def recording() -> np.ndarray:
        x = orig()
        if fuel._frozen_ln_m is not None:
            m = {e: float(x[s.unknown_names.index(f"{e}.m_dot")]) for e in (air, out)}
            held.append((m[out] - m[air], math.exp(fuel._frozen_ln_m)))
        return x

    s._build_x0 = recording  # type: ignore[method-assign]
    r = _quiet(s.solve)
    assert r["__success__"], r["__message__"]
    assert held, "the stage built no start while holding the fuel"
    for seeded, m_fuel in held:
        assert seeded == pytest.approx(m_fuel, rel=1e-9)
    assert fuel._frozen_ln_m is None
    assert fuel._given_m_dot() == 0.0


def test_a_stage_that_raises_holds_no_fuel_afterwards(monkeypatch: pytest.MonkeyPatch) -> None:
    """The held flow is cleared on every exit of the stage, an exception
    included, so no later start reads it."""
    s = NetworkSolver(_net(1500.0))
    fuel = s.network.nodes["fuel"]
    orig = NetworkSolver.solve

    def raising_inside(self: NetworkSolver, *a: object, **k: object) -> dict:
        if getattr(self, "_in_fuel_fixed_stage", False):
            assert fuel._frozen_ln_m is not None
            raise RuntimeError("stage failed")
        return orig(self, *a, **k)

    monkeypatch.setattr(NetworkSolver, "solve", raising_inside)
    with pytest.raises(RuntimeError, match="stage failed"):
        _quiet(s.solve)
    assert fuel._frozen_ln_m is None
    assert fuel._given_m_dot() == 0.0


def test_an_ejector_fed_by_a_fuel_controller_is_seeded_from_its_given_flow() -> None:
    """The ejector re-derives its closed-form seeds on every resolve, which
    every start runs; a fuel controller as its primary seeded them from the
    last iterate."""
    from test_ejector_element import (
        _CHANNEL_D,
        _DRY_AIR_Y,
        _MIXING_AREA,
        _NOZZLE_EXIT_AREA,
        _THROAT_AREA,
    )

    from combaero.network import (
        CombustorNode,
        FuelBoundary,
        MomentumChamberNode,
        PressureBoundary,
    )
    from combaero.network.components import ChannelElement, OrificeElement
    from combaero.network.ejector_element import EjectorElement

    def build(guess: float | None) -> FlowNetwork:
        n = FlowNetwork()
        n.add_node(
            FuelBoundary(
                "fuel", combustor="burner", Tt_target=1200.0, Tt=440.0, Y=CH4, m_dot_guess=guess
            )
        )
        n.add_node(PressureBoundary("air", Pt=1.0e5, Tt=300.0, Y=_DRY_AIR_Y))
        n.add_node(PressureBoundary("out", Pt=0.9e5, Tt=300.0, Y=_DRY_AIR_Y))
        area = math.pi * (_CHANNEL_D / 2.0) ** 2
        for k in ("mc_p", "mc_s", "mc_o"):
            n.add_node(MomentumChamberNode(k, area=area))
        n.add_node(CombustorNode("burner"))
        n.add_element(ChannelElement("ch_p", "fuel", "mc_p", length=0.2, diameter=_CHANNEL_D))
        n.add_element(ChannelElement("ch_s", "air", "mc_s", length=0.2, diameter=_CHANNEL_D))
        n.add_element(ChannelElement("ch_o", "mc_o", "burner", length=0.2, diameter=_CHANNEL_D))
        n.add_element(
            OrificeElement(
                "exit", "burner", "out", diameter=0.03, Cd=0.8, correlation="fixed", in_pipe=False
            )
        )
        n.add_element(
            EjectorElement(
                id="ej",
                primary_node="mc_p",
                secondary_node="mc_s",
                outlet_node="mc_o",
                throat_area=_THROAT_AREA,
                nozzle_exit_area=_NOZZLE_EXIT_AREA,
                mixing_area=_MIXING_AREA,
                regime="compressible",
                anchor_closure="kracik_dvorak",
            )
        )
        return n

    for guess in (None, 0.01):
        n = build(guess)
        n.resolve_all_topology()
        seeds = dict(n.elements["ej"]._auto_seeds)
        assert seeds
        n.nodes["fuel"].m_dot = 0.05  # what an evaluation leaves behind
        n.resolve_all_topology()
        assert n.elements["ej"]._auto_seeds == seeds, guess
