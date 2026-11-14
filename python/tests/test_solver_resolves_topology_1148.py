"""A solver evaluates the model a solve uses, built or solved (#1148, #1142).

``NetworkSolver.solve`` resolves the network's topology (each component reads
its neighbours: a merge chamber its inflow elements, an orifice its upstream
pipe); a solver only built (``_build_x0``) and evaluated did not. On a merge
chamber whose main stream is relayed (its T0 and Y move with the flows), the
unresolved chamber had no inflow elements: its main face lost the side jets'
momentum from F, and the main stream's T0 and Y columns from J. F of a fresh
solver differed from the solved one's by up to 2.8e3 at the same x
(incompressible); the segment row's FD error was 3.4e-3 to 3.1e-2. Every
FD audit on a fresh solver (#1142) checked that unresolved model.

``_build_x0`` -- which every residual and Jacobian evaluation needs -- now
resolves the topology. Mutation-tested: without that call all 8 cases of the
first test fail, and the 12 forward root and perturbed cases of the second
(at m = 0 and +-1e-7 the main stream's T0 column vanishes with its flow).
"""

from __future__ import annotations

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian, scaled_error
from test_crossflow_entry_gas_1137 import A_DUCT, CH4, K_IN, P0, Y, _merge

from combaero.network import (
    ChannelElement,
    CombustorNode,
    CrossflowSegmentElement,
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)


def _burning(regime: str, rev: bool = False, dp: float = 4000.0) -> FlowNetwork:
    """A burning combustor (air orifice + fixed CH4) as the segment's
    reservoir, the segment the main inlet of a merge chamber that a 900 K jet
    joins: the main stream's T0 AND composition are relayed."""
    p_sup, p_exit = (P0, P0 + dp) if rev else (P0 + dp, P0)
    g = FlowNetwork()
    g.add_node(PressureBoundary("exit", Pt=p_exit, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("air", Pt=p_sup, Tt=600.0, Y=Y))
    g.add_node(CombustorNode("res"))
    g.add_element(
        OrificeElement("o", "air", "res", diameter=0.03, Cd=0.8, correlation="fixed", in_pipe=False)
    )
    g.add_node(MassFlowBoundary("fuel", m_dot=2e-4, Tt=300.0, Y=CH4))
    g.add_element(LosslessConnectionElement("inj", "fuel", "res"))
    g.add_node(MomentumChamberNode("mc", area=A_DUCT, regime=regime, main_inlet="s1x"))
    g.add_element(
        CrossflowSegmentElement(
            "s1x", "res", "mc", length=0.05, area=A_DUCT, entry_K=K_IN, regime=regime
        )
    )
    g.add_node(PressureBoundary("hot", Pt=p_sup, Tt=900.0, Y=Y))
    g.add_element(OrificeElement("h", "hot", "mc", diameter=0.006, Cd=0.6, correlation="fixed"))
    g.add_element(ChannelElement("tail", "mc", "exit", length=0.05, diameter=0.02, regime=regime))
    return g


def _dense(J) -> np.ndarray:
    return J.toarray() if hasattr(J, "toarray") else np.asarray(J)


BUILDERS = {
    "mix": lambda regime, rev: _merge(regime, 2000.0, "chamber", rev=rev, reservoir="mix"),
    "burning": _burning,
}


@pytest.mark.parametrize("rev", [False, True])
@pytest.mark.parametrize("kind", ["mix", "burning"])
@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_fresh_solver_is_the_solved_model(regime: str, kind: str, rev: bool) -> None:
    """At the root of a solve, a freshly built solver's F and J are the solved
    solver's. Measured: identical incompressible; compressible within 2e-14
    (F) and 1e-12 of max|J| -- the round-off one solver shows between two
    evaluations at the same x. Unresolved, F differed by 2e-2 to 2.8e3."""

    def build() -> FlowNetwork:
        return BUILDERS[kind](regime, rev)

    s = NetworkSolver(build())
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    x = np.array(r["__x_solution__"], dtype=float)
    F1, J1 = s._residuals_and_jacobian(x)
    fresh = NetworkSolver(build())
    fresh._build_x0()
    assert fresh.unknown_names == s.unknown_names
    F2, J2 = fresh._residuals_and_jacobian(x)
    F1, F2 = np.asarray(F1, dtype=float), np.asarray(F2, dtype=float)
    J1, J2 = _dense(J1), _dense(J2)
    assert np.max(np.abs(F1 - F2)) < 1e-9
    assert np.max(np.abs(J1 - J2)) < 1e-12 * np.max(np.abs(J1))


def _segment_row_fd_error(build, x: np.ndarray) -> float:
    s = NetworkSolver(build())
    s._build_x0()
    J = s._residuals_and_jacobian(x)[1]
    fd = fd_jacobian(s, x, floor=1e-10)
    rows = [i for i, (label, _k) in enumerate(s._row_kinds(x)) if label.endswith("'s1x'")]
    return float(scaled_error(J, fd, x)[rows[0]])


@pytest.mark.parametrize("point", ["root", "perturbed", "m0", "m+", "m-"])
@pytest.mark.parametrize("rev", [False, True])
@pytest.mark.parametrize("res", ["mix", "chain", "mcres"])
@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_merge_face_relays_the_main_streams_t0(
    regime: str, res: str, rev: bool, point: str
) -> None:
    """The merge chamber's main face carries the main stream's relayed T0
    (dry reservoirs whose T0 moves with the flows, #1137's): the segment row
    on a fresh solver to 1e-7 scaled, at the root, perturbed, at m = 0 and
    either side of it, both directions. Measured <= 1.2e-8; unresolved it was
    1.2e-3 to 7.9e-3 forward. A reservoir whose composition moves (burning,
    wet) still lacks the entry drop's d/dY (#1146), so it is checked above
    for F and J identity only."""

    def build() -> FlowNetwork:
        return _merge(regime, 2000.0, "chamber", rev=rev, reservoir=res)

    s = NetworkSolver(build())
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    x = np.array(r["__x_solution__"], dtype=float)
    if point == "perturbed":
        x = x * (1.0 + 1e-3 * np.random.default_rng(7).uniform(-1.0, 1.0, x.size))
    elif point in ("m0", "m+", "m-"):
        x[s.unknown_names.index("s1x.m_dot")] = {"m0": 0.0, "m+": 1e-7, "m-": -1e-7}[point]
    assert _segment_row_fd_error(build, x) < 1e-7
