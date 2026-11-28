"""After the Newton homotopy, a failed cold solve gets continuation in its
boundary conditions (#642).

Every pressure boundary moves from the mean of their Pt to its own, every
flow and heat boundary from zero to its own; the gentle end is solved cold,
then each step warm-starts from the last root. Opt-in
(``NetworkSolver.bc_ramp_retry``); the GUI, NetworkRunner and the stress
harness turn it on. On the stress random draws the homotopy left unsolved it
rescued 3 of 21 (31, 34, 35).
"""

from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path
from typing import Any

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    EnergyBoundary,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
from networks import FAMILIES  # noqa: E402

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


@pytest.mark.slow
def test_a_draw_the_homotopy_left_is_rescued_by_the_ramp() -> None:
    """Stress random 48: the cold solve, its retries and the homotopy all
    fail; the ramp solves it. (Random 34 was the draw until the
    cold start learnt its compressible elements' kernels (#971), random 35
    until it seeded a fraction loss as its own ratio (#989), random 19
    until #1002 replaced its loop-closing fraction loss with a head loss,
    and random 111 until #1128 routed static-discharge and unanchored
    fraction losses through the LM fallback: all four now converge before
    the ramp.)"""
    net, _ = FAMILIES["random"](48)
    s = NetworkSolver(net)
    s.homotopy_retry = True
    s.bc_ramp_retry = True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=90.0)
    assert r["__success__"], r["__message__"]
    assert "boundary-condition ramp" in r["__message__"]


def _small() -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=3.0e5, Tt=400.0, Y=Y_AIR))
    g.add_node(PressureBoundary("b", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(MassFlowBoundary("m", m_dot=0.02, Tt=500.0, Y=Y_AIR))
    p = PlenumNode("p")
    p.add_energy_boundary(EnergyBoundary("q", Q=1.0e3, fraction=0.01))
    g.add_node(p)
    g.add_element(ChannelElement("c1", "a", "p", length=1.0, diameter=0.02))
    g.add_element(ChannelElement("c2", "m", "p", length=1.0, diameter=0.02))
    g.add_element(ChannelElement("c3", "p", "b", length=1.0, diameter=0.02))
    return g


def test_a_failed_ramp_restores_every_boundary_exactly(monkeypatch: pytest.MonkeyPatch) -> None:
    """Each attempt sees the ramped boundaries; afterwards every one is back
    bit for bit, and the failed result is the one returned."""
    g = _small()
    s = NetworkSolver(g)
    seen: list[tuple[float, float, float, float, float]] = []

    def fail(self: NetworkSolver, *a: Any, **kw: Any) -> dict[str, Any]:
        n = g.nodes
        eb = n["p"].energy_boundaries[0]
        seen.append((n["a"].Pt, n["b"].Pt, n["m"].m_dot, eb.Q, eb.fraction))
        return {"__success__": False}

    monkeypatch.setattr(NetworkSolver, "solve", fail)
    failed = {"__success__": False, "__message__": "cold failed"}
    assert s._bc_ramp_retry(failed, "hybr", 30.0, None, True, time.time()) is failed
    lam = 0.05
    assert seen[0] == pytest.approx(
        (2.0e5 + lam * 1.0e5, 2.0e5 - lam * 1.0e5, lam * 0.02, lam * 1.0e3, lam * 0.01)
    )
    assert len(seen) >= 2 and seen[1][2] == pytest.approx(0.5 * seen[0][2])  # gentler
    n = g.nodes
    eb = n["p"].energy_boundaries[0]
    assert (n["a"].Pt, n["b"].Pt, n["m"].m_dot, eb.Q, eb.fraction) == (
        3.0e5,
        1.0e5,
        0.02,
        1.0e3,
        0.01,
    )


def test_with_the_ramp_on_an_uncapped_homotopy_takes_half_the_time_left(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The GUI's homotopy has no cap of its own; with the ramp after it, it
    gets half of the time left, and its setting is restored."""
    g = _small()
    s = NetworkSolver(g)
    s.homotopy_retry = True
    s.homotopy_budget = None
    s.bc_ramp_retry = True
    budgets: list[float | None] = []
    calls: list[str] = []
    original = NetworkSolver.solve

    def cold_fails(self: NetworkSolver, *a: Any, **kw: Any) -> dict[str, Any]:
        if getattr(self, "_in_homotopy", False):
            return {"__success__": False, "__converged__": False, "__message__": "cold"}
        return original(self, *a, **kw)

    def homotopy(self: NetworkSolver, failed: dict, *a: Any) -> dict:
        budgets.append(self.homotopy_budget)
        calls.append("homotopy")
        return failed

    def ramp(self: NetworkSolver, failed: dict, *a: Any) -> dict:
        calls.append("ramp")
        return failed

    monkeypatch.setattr(NetworkSolver, "solve", cold_fails)
    monkeypatch.setattr(NetworkSolver, "_newton_homotopy_retry", homotopy)
    monkeypatch.setattr(NetworkSolver, "_bc_ramp_retry", ramp)
    s.solve(timeout=40.0)
    assert calls == ["homotopy", "ramp"]
    assert budgets[0] is not None and 15.0 < budgets[0] <= 20.0
    assert s.homotopy_budget is None
