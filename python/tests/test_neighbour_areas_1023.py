"""An element's explicit area against the pipe or chamber it adjoins (#1023).

The design rule (agreed 2026-11-06):

* an area equal to its neighbour is a pass-through, as a bore equal to its
  pipe is (#788): an AreaChangeElement with F0 == F1 is no area change, zero
  loss, the row Pt_in - Pt_out with the constant Jacobian of a lossless
  connection;
* an explicit AreaChange F0/F1 or tee F_C/F_branch that contradicts the pipe
  attached at its face is refused, naming both;
* an ejector mixing area larger than its outlet pipe is refused; a wider
  outlet pipe stays accepted;
* a momentum chamber smaller than an adjoining pipe is warned about, never
  refused, and never for the DEFAULT_CHAMBER_AREA placeholder;
* a tee branch wider than the main (psi < 1) stays accepted.

Equal means to NEIGHBOUR_AREA_RTOL on the equivalent diameter.
"""

from __future__ import annotations

import math
import os
import sys
import warnings

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    AreaChangeElement,
    ChannelElement,
    FlowNetwork,
    LosslessConnectionElement,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    TeeJunctionElement,
)
from combaero.network.components import (
    DEFAULT_CHAMBER_AREA,
    NEIGHBOUR_AREA_RTOL,
    NetworkMixtureState,
)
from combaero.network.ejector_element import EjectorElement
from gui.backend.resolve_preview import resolve_preview

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
D = 0.05
A = math.pi * D * D / 4.0


def _area_net(
    F0: float | None,
    F1: float | None,
    *,
    regime: str = "compressible",
    D_up: float = D,
    D_down: float = D,
    lossless: bool = False,
) -> FlowNetwork:
    """in -> [channel D_up] -> a -> area change -> b -> [channel D_down] -> out.
    lossless: a LosslessConnectionElement in the area change's place."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    n.add_node(PlenumNode("a"))
    n.add_node(PlenumNode("b"))
    n.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    n.add_element(ChannelElement("c1", "in", "a", length=1.0, diameter=D_up, regime=regime))
    n.add_element(ChannelElement("c2", "b", "out", length=1.0, diameter=D_down, regime=regime))
    if lossless:
        n.add_element(LosslessConnectionElement("x", "a", "b"))
    else:
        n.add_element(AreaChangeElement("x", "a", "b", F0=F0, F1=F1, regime=regime))
    return n


def _fd_error(s: NetworkSolver, x: np.ndarray) -> float:
    F, J = s._residuals_and_jacobian(x)
    J = J.toarray() if hasattr(J, "toarray") else np.asarray(J)
    Jfd = fd_jacobian(s, x)
    return float(np.max(np.abs(J - Jfd)) / np.max(np.abs(J)))


# --- rule 1: equal faces are a pass-through ----------------------------------
# Falsified by removing the _pass_through branch of AreaChangeElement.residuals
# (the closure then runs at n = 1: no loss to set the flow, the solve fails or
# the row is no Pa row), and by restoring the F0 == F1 refusal in validate.


@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
@pytest.mark.parametrize("F", [A, A * (1.0 + NEIGHBOUR_AREA_RTOL)])
def test_equal_faces_are_no_area_change(regime: str, F: float) -> None:
    """Same flow and pressures as a lossless connection in its place."""
    r = NetworkSolver(_area_net(A, F, regime=regime)).solve()
    ref = NetworkSolver(_area_net(None, None, regime=regime, lossless=True)).solve()
    assert r["__success__"] and ref["__success__"]
    assert r["x.m_dot"] == pytest.approx(ref["x.m_dot"], rel=1e-9)
    assert r["a.Pt"] == pytest.approx(r["b.Pt"], rel=1e-12)
    assert r["a.Pt"] == pytest.approx(ref["a.Pt"], rel=1e-9)


def test_the_pass_through_row_is_a_pressure_row_with_an_exact_jacobian() -> None:
    net = _area_net(None, None)  # both faces inherited from the D pipes: equal
    s = NetworkSolver(net)
    net.resolve_all_topology()
    x_el = net.elements["x"]
    assert x_el._pass_through and x_el._F0_source == "channel 'c1'"
    assert x_el.residual_scale_kind == "p"
    assert x_el._seed_resistance(1.2, 1e5, 0.1) == 0.0
    r = s.solve()
    assert r["__success__"]
    x = np.array(r["__x_solution__"])
    assert _fd_error(s, x) < 1e-6
    rng = np.random.default_rng(1023)
    assert _fd_error(s, x * (1.0 + 0.02 * rng.standard_normal(x.size))) < 1e-6
    st = NetworkMixtureState(P=1.5e5, Pt=1.6e5, Tt=300.0, m_dot=0.3, Y=Y_AIR)
    st2 = NetworkMixtureState(P=1.2e5, Pt=1.25e5, Tt=300.0, m_dot=0.3, Y=Y_AIR)
    res, jac = x_el.residuals(st, st2)
    assert res == [pytest.approx(0.35e5)]
    assert dict(jac) == {0: {"a.Pt": 1.0, "b.Pt": -1.0}}
    d = x_el.diagnostics(st, st2)
    assert d["dPt_loss"] == 0.0 and d["zeta"] == 0.0 and d["choked"] == 0.0
    assert d["area_ratio"] == 1.0 and d["mach_loss_ref"] > 0.0


def test_faces_just_outside_the_band_are_still_an_area_change() -> None:
    F1 = A * (1.0 + 10.0 * NEIGHBOUR_AREA_RTOL)
    net = _area_net(None, F1, D_down=math.sqrt(4.0 * F1 / math.pi))
    net.resolve_all_topology()
    assert not net.elements["x"]._pass_through
    assert net.elements["x"].residual_scale_kind == "mdot"


def test_a_loop_of_pass_throughs_is_a_flow_independent_loop() -> None:
    """Pt_out = Pt_in at every flow: two in parallel leave the split
    undetermined, as two lossless connections do (#988). Falsified by
    removing the pass-through case of graph._flow_independent_dpt."""
    n = _area_net(A, A)
    n.add_element(AreaChangeElement("x2", "a", "b", F0=A, F1=A))
    n.resolve_all_topology()
    with pytest.raises(ValueError, match="closed loop"):
        n.validate()


# --- rule 2: an explicit face that contradicts its pipe is refused -----------
# Falsified by removing the _check_face_against_pipe calls (AreaChange and tee
# resolve_topology): every case below is then accepted.


@pytest.mark.parametrize(
    ("F0", "F1", "name", "pipe"), [(4.0 * A, None, "F0", "c1"), (None, 0.25 * A, "F1", "c2")]
)
def test_an_area_change_face_contradicting_its_pipe_is_refused(
    F0: float | None, F1: float | None, name: str, pipe: str
) -> None:
    net = _area_net(F0, F1)
    with pytest.raises(ValueError, match=rf"'x': {name} = .* contradicts the channel '{pipe}'"):
        net.resolve_all_topology()


def test_an_area_change_face_equal_to_its_pipe_is_consistent() -> None:
    F1 = 2.0 * A
    net = _area_net(A * (1.0 + NEIGHBOUR_AREA_RTOL), F1, D_down=math.sqrt(4.0 * F1 / math.pi))
    r = NetworkSolver(net).solve()
    assert r["__success__"]
    assert not net.elements["x"]._pass_through


def test_a_face_meeting_a_chamber_of_its_own_area_is_not_compared_with_the_pipe() -> None:
    """The face meets the chamber, not the pipe beyond it: the chamber against
    the pipe is rule 4's warning, not a refusal (the stress random family
    draws such faces)."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    n.add_node(MomentumChamberNode("a", area=4.0 * A))
    n.add_node(MomentumChamberNode("b", area=8.0 * A))
    n.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    # The pipes are A, the faces 4 A and 8 A: their chambers' own areas.
    n.add_element(ChannelElement("c1", "in", "a", length=1.0, diameter=D))
    n.add_element(AreaChangeElement("x", "a", "b", F0=4.0 * A, F1=8.0 * A))
    n.add_element(ChannelElement("c2", "b", "out", length=1.0, diameter=D))
    n.resolve_all_topology()  # no refusal
    assert n.elements["x"].F0 == 4.0 * A


def _tee_net(
    F_C: float | None, F_branch: float | None, *, D_com: float = D, D_bra: float = D
) -> FlowNetwork:
    """A branching tee with a pipe on each arm."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("bc", Pt=2.1e5, Tt=300.0, Y=Y_AIR))
    n.add_node(PressureBoundary("bs", Pt=2.098e5, Tt=300.0, Y=Y_AIR))
    n.add_node(PressureBoundary("bb", Pt=2.06e5, Tt=300.0, Y=Y_AIR))
    for nid in ("pc", "ps", "pb"):
        n.add_node(MomentumChamberNode(nid))
    n.add_element(ChannelElement("cc", "bc", "pc", length=0.5, diameter=D_com))
    n.add_element(ChannelElement("cs", "ps", "bs", length=0.5, diameter=D_com))
    n.add_element(ChannelElement("cb", "pb", "bb", length=0.5, diameter=D_bra))
    n.add_element(
        TeeJunctionElement(
            "tee",
            common_node="pc",
            straight_node="ps",
            branch_node="pb",
            theta=math.pi / 2.0,
            F_C=F_C,
            F_branch=F_branch,
            tee_type="branching",
        )
    )
    return n


@pytest.mark.parametrize(
    ("F_C", "F_branch", "match"),
    [
        (10.0 * A, None, r"'tee': F_C = .* contradicts the channel 'c[cs]'"),
        (A, 10.0 * A, r"'tee': F_branch = .* contradicts the channel 'cb'"),
    ],
)
def test_a_tee_area_contradicting_its_arm_pipe_is_refused(
    F_C: float, F_branch: float | None, match: str
) -> None:
    with pytest.raises(ValueError, match=match):
        _tee_net(F_C, F_branch).resolve_all_topology()


def test_a_tee_branch_wider_than_the_main_is_accepted() -> None:
    """psi < 1 is buildable (rule 5); consistent with its pipes, nothing
    refuses it, given or inherited."""
    for F_C, F_b in ((A, 4.0 * A), (None, None)):
        net = _tee_net(F_C, F_b, D_bra=2.0 * D)
        net.resolve_all_topology()
        assert net.elements["tee"].psi == pytest.approx(0.25)


# --- rule 3: an ejector mixing area above its outlet pipe is refused ---------
# Falsified by removing _check_mixing_against_outlet_pipe's call.


def _ejector_net(mixing_area: float | None, D_out: float = D) -> FlowNetwork:
    g = FlowNetwork()
    for nid, p in (("bp", 7e5), ("bs", 0.5e5), ("bo", 1.0e5)):
        g.add_node(PressureBoundary(nid, Pt=p, Tt=300.0, Y=Y_AIR))
    for nid in ("mp", "ms", "mo"):
        g.add_node(MomentumChamberNode(nid, area=A))
    g.add_element(ChannelElement("cp", "bp", "mp", length=0.2, diameter=D))
    g.add_element(ChannelElement("cs", "bs", "ms", length=0.2, diameter=D))
    g.add_element(ChannelElement("co", "mo", "bo", length=0.2, diameter=D_out))
    g.add_element(
        EjectorElement(
            id="ej",
            primary_node="mp",
            secondary_node="ms",
            outlet_node="mo",
            throat_area=3.14e-5,
            nozzle_exit_area=4.675e-5,
            mixing_area=mixing_area,
        )
    )
    return g


def test_a_mixing_area_larger_than_the_outlet_pipe_is_refused() -> None:
    with pytest.raises(ValueError, match=r"'ej': mixing_area .* larger than its outlet pipe"):
        _ejector_net(1.2 * A).resolve_all_topology()


@pytest.mark.parametrize("A3", [3.069e-4, A, A * (1.0 + NEIGHBOUR_AREA_RTOL)])
def test_a_mixing_area_at_or_inside_the_outlet_pipe_is_accepted(A3: float) -> None:
    """A wider outlet pipe is a diffuser exit; equal is the duct continuing."""
    net = _ejector_net(A3)
    net.resolve_all_topology()
    assert net.elements["ej"].mixing_area == A3


def test_the_stress_ejector_family_fits_its_mixing_duct_into_the_pipe() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(here, "..", "..", "scripts", "stress"))
    from networks import FAMILIES

    for seed in (1, 11, 68, 93):
        net, _ = FAMILIES["ejector"](seed)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            net.resolve_all_topology()
        assert net.elements["ej"].mixing_area == pytest.approx(A, rel=1e-15)


# --- rule 4: a chamber narrower than an adjoining pipe warns -----------------
# Falsified by removing the warning loop in FlowNetwork.validate (no warning)
# or the comparison in MomentumChamberNode.pipe_area_warnings.


def _chamber_net(area: float | None, D_up: float = D, D_down: float = D) -> FlowNetwork:
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    n.add_node(MomentumChamberNode("mc", area=area))
    n.add_node(PressureBoundary("out", Pt=1.9e5, Tt=300.0, Y=Y_AIR))
    n.add_element(ChannelElement("c1", "in", "mc", length=1.0, diameter=D_up))
    n.add_element(ChannelElement("c2", "mc", "out", length=1.0, diameter=D_down))
    return n


def _chamber_warnings(net: FlowNetwork) -> list[str]:
    net.resolve_all_topology()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        net.validate()
    return [str(x.message) for x in w if "MomentumChamberNode" in str(x.message)]


def test_a_chamber_narrower_than_its_pipes_warns_naming_both() -> None:
    msgs = _chamber_warnings(_chamber_net(1e-4))
    assert len(msgs) == 2
    for msg, pipe in zip(msgs, ("c1", "c2"), strict=True):
        assert "'mc'" in msg and f"'{pipe}'" in msg
        assert "1e-04 m^2" in msg or "0.0001 m^2" in msg
        assert f"{A:.6g} m^2" in msg and "AreaChangeElement" in msg
    # A warning, not a refusal: a chamber at half the pipe's area solves.
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        assert NetworkSolver(_chamber_net(0.5 * A)).solve()["__success__"]
    assert any("smaller than the adjoining pipe" in str(x.message) for x in w)


def test_an_inherited_chamber_area_narrower_than_the_far_pipe_warns() -> None:
    """Inherited from the upstream pipe, the chamber is narrower than the
    wider pipe downstream: an expansion with no loss."""
    msgs = _chamber_warnings(_chamber_net(None, D_up=D, D_down=2.0 * D))
    assert len(msgs) == 1 and "'c2'" in msgs[0] and "expands" in msgs[0]


@pytest.mark.parametrize("area", [A, A * (1.0 - NEIGHBOUR_AREA_RTOL), 4.0 * A])
def test_a_chamber_at_least_its_pipes_area_does_not_warn(area: float) -> None:
    assert _chamber_warnings(_chamber_net(area)) == []


def test_a_chamber_sized_by_a_pipe_sized_later_does_not_warn() -> None:
    """The pipe into the chamber takes its diameter from an area change
    upstream only in the element pass, after the chamber's first look found
    nothing: one resolve pass left the chamber at the default, the fixed point
    (what every solve sees, #1148) gives it the pipe's 0.5 m^2. Equal to the
    pipe, so nothing warns."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    n.add_node(PlenumNode("n1"))
    n.add_node(MomentumChamberNode("mc"))
    n.add_node(PressureBoundary("out", Pt=1.9e5, Tt=300.0, Y=Y_AIR))
    n.add_element(AreaChangeElement("x", "in", "n1", F0=A, F1=0.5))
    n.add_element(ChannelElement("c", "n1", "mc", length=1.0))
    n.add_element(LosslessConnectionElement("l", "mc", "out"))
    assert _chamber_warnings(n) == []
    mc, c = n.nodes["mc"], n.elements["c"]
    assert mc._area_source == "channel 'c'"
    assert math.pi * c.diameter**2 / 4.0 == pytest.approx(0.5)
    assert mc.area == pytest.approx(0.5)


def test_the_default_chamber_area_never_warns() -> None:
    """A chamber nothing sizes holds DEFAULT_CHAMBER_AREA, a placeholder: fed
    through a lumped restriction (an orifice with in_pipe=False) and left by
    a lossless connection, it has no neighbouring geometry at the fixed
    point. The placeholder is not compared with anything, so nothing warns."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    n.add_node(MomentumChamberNode("mc"))
    n.add_node(PressureBoundary("out", Pt=1.9e5, Tt=300.0, Y=Y_AIR))
    n.add_element(
        OrificeElement("o", "in", "mc", diameter=0.01, Cd=0.6, correlation="fixed", in_pipe=False)
    )
    n.add_element(LosslessConnectionElement("l", "mc", "out"))
    assert _chamber_warnings(n) == []
    mc = n.nodes["mc"]
    assert mc._area_source == "default"
    assert mc.area == DEFAULT_CHAMBER_AREA


def test_the_resolve_preview_carries_the_chamber_warning() -> None:
    rep = resolve_preview(_chamber_net(1e-4))
    assert rep["mc"].status == "resolved"
    assert "smaller than the adjoining pipe 'c1'" in rep["mc"].warning
    assert resolve_preview(_chamber_net(A))["mc"].warning == ""


# --- the GUI preview ---------------------------------------------------------


def test_the_resolve_preview_reports_the_pass_through_and_the_refusal() -> None:
    rep = resolve_preview(_area_net(None, None))
    assert rep["x"].status == "resolved"
    assert rep["x"].fields["pass_through"].value == 1.0
    rep = resolve_preview(_area_net(4.0 * A, None))
    assert rep["x"].status == "unresolved"
    assert "contradicts the channel 'c1'" in rep["x"].reason
    rep = resolve_preview(_ejector_net(1.2 * A))
    assert rep["ej"].status == "unresolved"
    assert "outlet pipe" in rep["ej"].reason


def _gui(**ac: float) -> dict:
    def node(i: str, t: str, **data: object) -> dict:
        return {"id": i, "type": t, "position": {"x": 0.0, "y": 0.0}, "data": data}

    def edge(i: str, s: str, t: str) -> dict:
        return {"id": i, "source": s, "target": t, "data": {}}

    return {
        "nodes": [
            node("in", "pressure_boundary", Pt=2.0e5, Tt=300.0),
            node("c1", "channel", L=1.0, D=D),
            node("ac", "area_change", **ac),
            node("c2", "channel", L=1.0, D=D),
            node("out", "pressure_boundary", Pt=1.0e5, Tt=300.0),
        ],
        "edges": [
            edge("e1", "in", "c1"),
            edge("e2", "c1", "ac"),
            edge("e3", "ac", "c2"),
            edge("e4", "c2", "out"),
        ],
    }


@pytest.mark.parametrize(
    ("ac", "status", "pass_through"),
    [({}, "resolved", 1.0), ({"F0": A}, "resolved", 1.0), ({"F0": 4.0 * A}, "unresolved", 0.0)],
)
def test_the_gui_previews_an_area_change_between_equal_pipes(
    ac: dict, status: str, pass_through: float
) -> None:
    """Drawn between two D pipes (through the builder's junction chambers,
    which take the pipes' areas): inherited or equal faces are no area
    change; a face 4x the pipe is refused."""
    from gui.backend.graph_builder import build_network_from_schema
    from gui.backend.schemas import NetworkGraphSchema

    rep = resolve_preview(build_network_from_schema(NetworkGraphSchema(**_gui(**ac))))
    assert rep["ac"].status == status
    assert rep["ac"].fields["pass_through"].value == pass_through
    if status == "unresolved":
        assert "contradicts the channel 'c1'" in rep["ac"].reason
