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

Mutation-tested (97 cases): no resolve in ``_build_x0`` (or none on a
solver's first call): 39 fail; the resolve after the unknowns are gathered:
the 13 start cases fail; a single resolve pass (no fixed point): 8 fail (the
4 loss-chain and the GUI-chain identity cases, the fixed-point case, the C3X
and random-30 idempotency cases); the inferred-area round trip restored: the
C3X and random-30 idempotency and the 2 area-after-resolve cases fail; the
derived Dh kept even after the user changes it: the 2 Dh-after-resolve
cases fail; the GUI resolve preview back to one pass: both preview cases
fail.
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


def _deep(o, depth: int = 0, seen: frozenset = frozenset()):  # noqa: ANN001, ANN202
    """A comparable deep copy of a component's state: scalars by repr (NaN
    equals itself), arrays by bytes, other components by id, cycles cut."""
    if o is None or isinstance(o, (bool, int, str)):
        return o
    if isinstance(o, float):
        return repr(o)
    if isinstance(o, np.ndarray):
        return ("nd", o.shape, np.ascontiguousarray(o).tobytes())
    if isinstance(o, (np.floating, np.integer)):
        return repr(o.item())
    if depth > 6:
        return "<deep>"
    if isinstance(o, (list, tuple)):
        return [type(o).__name__] + [_deep(x, depth + 1, seen) for x in o]
    if isinstance(o, (set, frozenset)):
        return ["set"] + sorted(repr(_deep(x, depth + 1, seen)) for x in o)
    if isinstance(o, dict):
        return {repr(k): _deep(v, depth + 1, seen) for k, v in o.items()}
    if isinstance(o, FlowNetwork):
        return "<net>"
    d = getattr(o, "__dict__", None)
    if d is None:
        return repr(o) if not callable(o) else "<callable>"
    if depth > 0 and getattr(o, "id", None) is not None:
        return f"<ref {type(o).__name__} {o.id}>"
    if id(o) in seen:
        return "<cycle>"
    return {k: _deep(v, depth + 1, seen | {id(o)}) for k, v in d.items()}


def _deep_snapshot(g: FlowNetwork) -> dict:
    return {
        f"{kind}:{cid}": _deep(comp)
        for kind, comps in (("node", g.nodes), ("elem", g.elements))
        for cid, comp in comps.items()
    }


@pytest.mark.parametrize("name", sorted(IDEMPOTENCY_NETS))
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_a_second_resolve_changes_nothing(name: str) -> None:
    """Resolving a resolved network leaves the state of every node and
    element as it was (a deep snapshot, not just the scalars the fixed point
    compares), and an area a node inferred is the inferred
    area itself: the area -> Dh -> area round trip used to move it by an ulp
    on every resolve after the first (the C3X burner's 0.1 m^2 reported as
    0.09999999999999998 on some solves and 0.1 on others)."""
    from combaero.network.components import _infer_flow_area

    g = IDEMPOTENCY_NETS[name]()
    g.resolve_all_topology()
    once = _deep_snapshot(g)
    g.resolve_all_topology()
    assert _deep_snapshot(g) == once
    for nid, node in g.nodes.items():
        if getattr(node, "_Dh_derived", None) is None or node.Dh != node._Dh_derived:
            continue
        area, source = _infer_flow_area(g, [nid])
        if source == node._area_source:
            assert node.area == area, nid
    if name == "c3x_gui":
        assert g.nodes["burner"].area == 0.1


def _inferring_node_network(kind: str) -> FlowNetwork:
    """A chamber or combustor that infers its area from the channel after it
    (no channel feeds it): its Dh is derived from that area."""
    from combaero.network.components import ConstantFractionLoss, PressureLossElement

    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=2.0e5, Tt=300.0))
    g.add_node(PressureBoundary("out", Pt=1.9e5, Tt=300.0))
    g.add_node(PlenumNode("p"))
    g.add_node(MomentumChamberNode("m") if kind == "chamber" else CombustorNode("m"))
    g.add_element(ChannelElement("c1", "in", "p", length=0.2, diameter=0.01))
    g.add_element(PressureLossElement("pl", "p", "m", correlation=ConstantFractionLoss(xi=0.05)))
    g.add_element(ChannelElement("c2", "m", "out", length=0.2, diameter=0.05))
    return g


@pytest.mark.parametrize("kind", ["chamber", "combustor"])
def test_a_dh_set_after_a_resolve_sets_the_area(kind: str) -> None:
    """The inferred area stands only while Dh is still the one derived from
    it: a Dh the user sets after a resolve sets the area on the next."""
    g = _inferring_node_network(kind)
    g.resolve_all_topology()
    m = g.nodes["m"]
    assert m.area == np.pi * 0.05**2 / 4.0
    m.Dh = 0.02
    g.resolve_all_topology()
    assert m.area == pytest.approx(np.pi * 0.02**2 / 4.0, rel=1e-15)


@pytest.mark.parametrize("kind", ["chamber", "combustor"])
def test_an_area_set_after_a_resolve_stands(kind: str) -> None:
    """An area attribute set after a resolve, with the derived Dh untouched,
    stands on the next resolve (the area is the primary quantity here; it is
    no longer derived back from the stale Dh)."""
    g = _inferring_node_network(kind)
    g.resolve_all_topology()
    m = g.nodes["m"]
    m.area = 0.003
    g.resolve_all_topology()
    assert m.area == 0.003


@pytest.mark.parametrize("name", ["gui_channel_chain", "c3x_gui"])
@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_the_resolve_preview_reports_the_geometry_a_solve_uses(name: str) -> None:
    """The GUI's resolve preview runs the same fixed point: every element's and
    node's geometry and source are the solved network's. On the junction chain
    one node pass reported the junction at 0.1 m^2 "default" while the solve
    used 7.854e-3 m^2 from the defaulted channel."""
    from gui.backend.resolve_preview import _extract, resolve_preview

    build = {"gui_channel_chain": lambda: _gui_channel_chain("incompressible", False)}.get(
        name, _c3x_gui_network
    )
    preview = resolve_preview(build())
    solved = NetworkSolver(build())
    solved.solve()
    net = solved.network
    assert set(preview) == set(net.nodes) | set(net.elements)
    for item_id, item in {**net.nodes, **net.elements}.items():
        assert preview[item_id].fields == _extract(item), item_id
    if name == "gui_channel_chain":
        jct = preview["__junction__chan__chan2"].fields
        assert jct["area"].value == pytest.approx(np.pi * 0.1**2 / 4.0, rel=1e-15)
        assert jct["area"].source != "default"
