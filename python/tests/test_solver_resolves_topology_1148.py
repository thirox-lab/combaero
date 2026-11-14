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
resolves the topology, and ``resolve_all_topology`` runs to its fixed point:
a node that inherits from an element whose own geometry is inherited (a
chamber's area from a loss element's, a GUI junction's from a defaulted
channel) settles only on a later pass, so one pass still left a built solver
with another model than ``solve()``, which resolves several times.

Mutation-tested: no resolve in ``_build_x0`` (or none on a solver's first
call): 39 of 91 cases fail; the resolve after the unknowns are gathered: the
13 start cases fail; a single resolve pass (no fixed point): 8 fail (the 4
loss-chain and the GUI-chain identity cases, the fixed-point case, and the
C3X and random-30 idempotency cases);
the inferred-area round trip restored in the chamber and combustor nodes:
the C3X idempotency case fails.
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
    PlenumNode,
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


def _loss_chain(regime: str, rev: bool) -> FlowNetwork:
    """A MomentumChamberNode between two pressure-loss elements: the chamber's
    area comes from the loss element's, which is inherited from the channel
    upstream -- settled on the second resolve pass, not the first."""
    from combaero.network.components import ConstantFractionLoss, PressureLossElement

    p_in, p_out = (1.0e5, 1.05e5) if rev else (1.05e5, 1.0e5)
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=p_in, Tt=300.0))
    g.add_node(PressureBoundary("out", Pt=p_out, Tt=300.0))
    g.add_node(PlenumNode("p"))
    g.add_node(MomentumChamberNode("m", regime=regime))
    g.add_element(ChannelElement("c1", "in", "p", length=0.2, diameter=0.01, regime=regime))
    g.add_element(PressureLossElement("pl", "p", "m", correlation=ConstantFractionLoss(xi=0.05)))
    g.add_element(PressureLossElement("pl2", "m", "out", correlation=ConstantFractionLoss(xi=0.05)))
    return g


def _gui_channel_chain(regime: str, rev: bool) -> FlowNetwork:
    """The GUI's two-channel chain with no diameter given: the first channel
    takes the marked default, the junction between them its area from it --
    on the second pass (test_resolve_preview)."""
    from test_resolve_preview import _two_channel_graph

    from gui.backend.graph_builder import build_network_from_schema
    from gui.backend.schemas import NetworkGraphSchema

    return build_network_from_schema(NetworkGraphSchema(**_two_channel_graph(None, None)))


def _dense(J) -> np.ndarray:
    return J.toarray() if hasattr(J, "toarray") else np.asarray(J)


BUILDERS = {
    "mix": lambda regime, rev: _merge(regime, 2000.0, "chamber", rev=rev, reservoir="mix"),
    "burning": _burning,
}
#: Networks whose geometry settles only on a later resolve pass.
CHAINS = {"loss_chain": _loss_chain, "gui_channel_chain": _gui_channel_chain}


@pytest.mark.parametrize("rev", [False, True])
@pytest.mark.parametrize("kind", ["mix", "burning", "loss_chain"])
@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_fresh_solver_is_the_solved_model(regime: str, kind: str, rev: bool) -> None:
    """At the root of a solve, a freshly built solver's F and J are the solved
    solver's. Measured: identical incompressible; compressible within 2e-14
    (F) and 1e-12 of max|J| -- the round-off one solver shows between two
    evaluations at the same x. Unresolved, F differed by 2e-2 to 2.8e3."""

    _assert_fresh_is_solved(lambda: {**BUILDERS, **CHAINS}[kind](regime, rev))


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_a_gui_junction_chain_fresh_solver_is_the_solved_model() -> None:
    """The GUI chain's junction takes its area on the second pass: one pass
    left it at 0.1 m^2 against 7.85e-3 after the solve, F off by 5.9e4."""
    _assert_fresh_is_solved(lambda: _gui_channel_chain("incompressible", False))


def _assert_fresh_is_solved(build) -> None:
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


START_CASES = [
    (regime, kind, rev)
    for kind in ("mix", "burning", "loss_chain")
    for regime in ("incompressible", "compressible")
    for rev in (False, True)
] + [("incompressible", "gui_channel_chain", False)]


@pytest.mark.parametrize(("regime", "kind", "rev"), START_CASES)
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_the_start_is_the_resolved_networks(regime: str, kind: str, rev: bool) -> None:
    """A built solver's start and unknowns are those of a network resolved
    beforehand: the resolve comes before the unknowns are gathered, not after."""

    def build() -> FlowNetwork:
        return {**BUILDERS, **CHAINS}[kind](regime, rev)

    pre = build()
    pre.resolve_all_topology()
    a, b = NetworkSolver(build()), NetworkSolver(pre)
    x_a, x_b = np.asarray(a._build_x0()), np.asarray(b._build_x0())
    assert a.unknown_names == b.unknown_names
    assert np.array_equal(x_a, x_b)


def test_resolve_reaches_its_fixed_point() -> None:
    """One call settles the geometry: a second changes nothing."""
    g = _loss_chain("incompressible", False)
    g.resolve_all_topology()
    once = g._resolved_signature()
    g.resolve_all_topology()
    assert g._resolved_signature() == once
    assert g.nodes["m"].area == pytest.approx(np.pi * 0.01**2 / 4.0, rel=1e-15)


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


def _c3x_gui_network() -> FlowNetwork:
    import json
    from pathlib import Path

    from gui.backend.graph_builder import build_network_from_schema
    from gui.backend.schemas import NetworkGraphSchema

    path = Path(__file__).resolve().parents[2] / "gui" / "examples" / "hylton1983_c3x_vane.json"
    return build_network_from_schema(NetworkGraphSchema(**json.loads(path.read_text())))


def _stress_network(family: str, seed: int) -> FlowNetwork:
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
    from networks import FAMILIES

    return FAMILIES[family](seed)[0]


IDEMPOTENCY_NETS = {
    "c3x_gui": _c3x_gui_network,
    "random_30": lambda: _stress_network("random", 30),
    "ejector_3": lambda: _stress_network("ejector", 3),
    "ideal_nport_5": lambda: _stress_network("ideal_nport", 5),
}


@pytest.mark.parametrize("name", sorted(IDEMPOTENCY_NETS))
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_a_second_resolve_changes_nothing(name: str) -> None:
    """Resolving a resolved network leaves every scalar attribute of every
    node and element as it was, and an area a node inferred is the inferred
    area itself: the area -> Dh -> area round trip used to move it by an ulp
    on every resolve after the first (the C3X burner's 0.1 m^2 reported as
    0.09999999999999998 on some solves and 0.1 on others)."""
    from combaero.network.components import _infer_flow_area

    g = IDEMPOTENCY_NETS[name]()
    g.resolve_all_topology()
    once = g._resolved_signature()
    g.resolve_all_topology()
    assert g._resolved_signature() == once
    for nid, node in g.nodes.items():
        if not getattr(node, "_Dh_from_area", False):
            continue
        area, source = _infer_flow_area(g, [nid])
        if source == node._area_source:
            assert node.area == area, nid
    if name == "c3x_gui":
        assert g.nodes["burner"].area == 0.1
