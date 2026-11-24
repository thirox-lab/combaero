"""A node's streams all arrive (#1239).

The propagation points every stream it hands a node by the sign of its flow:
a reversed two-port element delivers -m_dot from its to_node (#481), and a
junction with a reversed port delivers by the sign of its port flows (#517;
the legacy tee, #1231). So no node takes a leaving stream, and a node's
signed inflow sum -- the throughput ``_total_m_dot`` its momentum closure,
merge face, coefficient and wall relays read -- is the total of the streams
that arrive, the one the C++ mixer and combustor kernels weight (#1232).
Two exceptions are known. An inline MassFlowBoundary's own stream (an
extraction) is negative, and only its own node takes it; no swept case has
one. A junction port fed by a legacy tee has no flow index to reverse by, so
it can be handed a leaving stream (#1301); no builder draws that topology,
and the network should refuse it.

#1239 asked whether the combustor's signed sums could shrink or change sign
while gas still arrives. Outside #1301 they cannot; this pins why. The points
are the Jacobian sweep's (``_jacobian_sweep.points``): every sign pattern of
each element's legs at two magnitude ratios, a flow at exactly zero, flows
near zero, and side streams of another gas into junction ports and interior
nodes.

A direct caller of ``combustor_residuals_and_jacobians`` can still pass
streams of both signs; the kernel's pressure-loss context then carries the
arriving total too, the flow the kernel mixes and burns.
"""

from __future__ import annotations

from collections import Counter
from typing import Any

import numpy as np
import pytest
from _jacobian_sweep import _quiet, case_points, ci_cases, resolved_solver

import combaero as cb
from combaero import _solver_tools as _st
from combaero.network import components as cbn

_NODES = (cbn.PlenumNode, cbn.MomentumChamberNode, cbn.CombustorNode, cbn.WallNode)


def _record(monkeypatch: pytest.MonkeyPatch, seen: list) -> None:
    """Wrap each node class's compute_derived_state to log its streams."""
    for cls in (*_NODES, cbn.MassFlowBoundary):
        real = cls.compute_derived_state

        def wrapped(self: Any, upstream_states: list, _real: Any = real) -> Any:
            seen.append(
                (self, [(getattr(s, "_element_id", None), float(s.m_dot)) for s in upstream_states])
            )
            return _real(self, upstream_states)

        monkeypatch.setattr(cls, "compute_derived_state", wrapped)


def _reversed_feeds(s: Any, x: np.ndarray, node: Any, streams: list) -> set[str]:
    """How a stream into ``node`` reached it against its element's declared
    direction: a reversed two-port element, or a reversed junction port."""
    out = set()
    rev2 = s._reversed_elements(x)
    revj = getattr(s, "_rev_junctions", {})
    for eid, _ in streams:
        if eid in rev2:
            out.add("two_port")
        elif eid in revj:
            out.add("junction")
    return out


@pytest.mark.slow
def test_no_node_takes_a_leaving_stream(monkeypatch: pytest.MonkeyPatch) -> None:
    seen: list = []
    _record(monkeypatch, seen)
    bad = []
    covered: Counter = Counter()
    n_points = 0
    for case in ci_cases(small_only=False):
        try:
            _, pts = case_points(case, None)
        except Exception:  # noqa: BLE001 -- a variant the network refuses to build
            continue
        for pt in pts:
            with _quiet():
                s = resolved_solver(case)
                seen.clear()
                try:
                    s._residuals(pt.x)
                except Exception:  # noqa: BLE001 -- a point the model refuses
                    continue
            n_points += 1
            for node, streams in seen:
                kind = type(node).__name__
                neg = [(e, m) for e, m in streams if m < 0.0]
                if isinstance(node, cbn.MassFlowBoundary):
                    # Its own injection: no element id, its own m_dot.
                    neg = [(e, m) for e, m in neg if not (e is None and m == float(node.m_dot))]
                if neg:
                    bad.append(f"{case.label()} {pt.desc}: {kind} {node.id} {neg}")
                for how in _reversed_feeds(s, pt.x, node, streams):
                    covered[(kind, how)] += 1
    assert not bad, "\n".join(bad[:20])
    # The points reach each way a stream gets turned round (2626 points on
    # this commit): a combustor fed through a reversed element (66), a
    # chamber through one (299) and through a reversed junction port (621).
    # No builder puts a combustor on a junction port; the propagation that
    # points the streams does not depend on the node's class.
    assert n_points > 1000
    for key in (
        ("CombustorNode", "two_port"),
        ("MomentumChamberNode", "two_port"),
        ("MomentumChamberNode", "junction"),
        ("PlenumNode", "junction"),
    ):
        assert covered[key] > 0, key


AIR = list(cb.mole_to_mass(cb.species.dry_air()))
FUEL = [0.0] * len(AIR)
FUEL[cb.species.species_index("CH4")] = 1.0


def _loss_total(streams: list, method: str = "complete") -> float:
    """The mdot_total the kernel hands its pressure-loss callable."""
    got: list[float] = []

    def loss(ctx: Any) -> tuple[float, float]:
        got.append(float(ctx.mdot_total))
        return (0.0, 0.0)

    _st.combustor_residuals_and_jacobians(
        streams, pressure_loss=loss, use_equilibrium=(method == "equilibrium")
    )
    assert len(got) == 1
    return got[0]


@pytest.mark.parametrize("method", ["complete", "equilibrium"])
def test_the_loss_sees_the_arriving_total_beside_a_leaving_stream(method: str) -> None:
    # 0.1 kg/s air and 0.004 kg/s methane arrive; 0.103 kg/s leaves through
    # a third stream. The signed sum, 0.001 kg/s, is what the loss saw.
    s = [
        cb.MassStream(0.1, 700.0, 2.0e5, AIR),
        cb.MassStream(0.004, 300.0, 2.1e5, FUEL),
        cb.MassStream(-0.103, 1500.0, 1.9e5, AIR),
    ]
    assert _loss_total(s, method) == pytest.approx(0.104, rel=1e-15)


def test_one_sign_keeps_the_signed_sum() -> None:
    # Streams of one sign: the arriving total is their sum, as before. With
    # none arriving the kernel mixes every stream (the #1224 reverse mean),
    # so the loss sees their (negative) sum.
    fwd = [cb.MassStream(0.1, 700.0, 2.0e5, AIR), cb.MassStream(0.004, 300.0, 2.1e5, FUEL)]
    assert _loss_total(fwd) == 0.1 + 0.004
    rev = [cb.MassStream(-0.1, 700.0, 2.0e5, AIR), cb.MassStream(-0.004, 300.0, 2.1e5, FUEL)]
    assert _loss_total(rev) == -0.1 - 0.004
