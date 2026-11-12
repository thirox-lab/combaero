"""The orifice bore against the pipe it sits in (#788, 0.8.0).

* bore == pipe (beta = 1, to ORIFICE_BORE_PIPE_RTOL) is NO orifice: zero
  loss, the flow passes unchanged, the row is Pt_in - Pt_out with the
  constant Jacobian of a lossless connection;
* bore > pipe is refused;
* a lumped Cd*A restriction is written explicitly, in_pipe=False: no pipe
  diameter, beta = 0, the wall hole the old silent fallback produced.
"""

from __future__ import annotations

import os
import sys
import warnings

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    LosslessConnectionElement,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)
from combaero.network.components import ORIFICE_BORE_PIPE_RTOL, NetworkMixtureState
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.resolve_preview import resolve_preview
from gui.backend.schemas import NetworkGraphSchema

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
D = 0.05


def _net(
    bore: float | None,
    *,
    in_pipe: bool = True,
    D_up: float = D,
    D_down: float = D,
    regime: str = "compressible",
    correlation: str = "fixed",
    pipes: bool = True,
) -> FlowNetwork:
    """in -> [channel D_up] -> a -> orifice -> b -> [channel D_down] -> out.
    bore None: a lossless connection in the orifice's place."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    n.add_node(PlenumNode("a"))
    n.add_node(PlenumNode("b"))
    n.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    if pipes:
        n.add_element(ChannelElement("c1", "in", "a", length=1.0, diameter=D_up, regime=regime))
        n.add_element(ChannelElement("c2", "b", "out", length=1.0, diameter=D_down, regime=regime))
    else:
        n.add_element(LosslessConnectionElement("c1", "in", "a"))
        n.add_element(LosslessConnectionElement("c2", "b", "out"))
    if bore is None:
        n.add_element(LosslessConnectionElement("o", "a", "b"))
    else:
        n.add_element(
            OrificeElement(
                "o",
                "a",
                "b",
                Cd=0.6 if correlation == "fixed" else None,
                diameter=bore,
                correlation=correlation,
                regime=regime,
                in_pipe=in_pipe,
            )
        )
    return n


def _fd_error(s: NetworkSolver, x: np.ndarray) -> float:
    F, J = s._residuals_and_jacobian(x)
    J = J.toarray() if hasattr(J, "toarray") else np.asarray(J)
    F = np.asarray(F)
    Jfd = fd_jacobian(s, x)
    return float(np.max(np.abs(J - Jfd)) / np.max(np.abs(J)))


# --- bore == pipe: no orifice ------------------------------------------------


@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
@pytest.mark.parametrize("bore", [D, D * (1.0 + 0.5 * ORIFICE_BORE_PIPE_RTOL)])
def test_a_bore_equal_to_its_pipe_is_no_orifice(regime: str, bore: float) -> None:
    """Same flow and pressures as a lossless connection in its place."""
    r = NetworkSolver(_net(bore, regime=regime)).solve()
    ref = NetworkSolver(_net(None, regime=regime)).solve()
    assert r["__success__"] and ref["__success__"]
    assert r["o.m_dot"] == pytest.approx(ref["o.m_dot"], rel=1e-9)
    assert r["a.Pt"] == pytest.approx(r["b.Pt"], rel=1e-12)
    assert r["a.Pt"] == pytest.approx(ref["a.Pt"], rel=1e-9)


def test_the_pass_through_row_is_a_pressure_row_with_an_exact_jacobian() -> None:
    net = _net(D)
    s = NetworkSolver(net)
    net.resolve_all_topology()
    o = net.elements["o"]
    assert o._pass_through and o.beta == 1.0
    assert o.residual_scale_kind == "p"
    assert o._seed_resistance(1.2, 1e5, 0.1) == 0.0
    r = s.solve()
    assert r["__success__"]
    x = np.array(r["__x_solution__"])
    assert _fd_error(s, x) < 1e-6
    # Off the root too: the row is linear, its Jacobian finite and constant.
    rng = np.random.default_rng(788)
    assert _fd_error(s, x * (1.0 + 0.02 * rng.standard_normal(x.size))) < 1e-6
    st = NetworkMixtureState(P=1.5e5, Pt=1.6e5, Tt=300.0, m_dot=0.3, Y=Y_AIR)
    st2 = NetworkMixtureState(P=1.2e5, Pt=1.25e5, Tt=300.0, m_dot=0.3, Y=Y_AIR)
    res, jac = o.residuals(st, st2)
    assert res == [pytest.approx(0.35e5)]
    assert jac == {0: {"a.Pt": 1.0, "b.Pt": -1.0}}


def test_a_bore_equal_to_the_reverse_pipe_is_no_orifice_too() -> None:
    """The pipe on the to-side sets beta in reverse: equal there is equal."""
    net = _net(0.03, D_up=0.04, D_down=0.03)
    net.resolve_all_topology()
    assert net.elements["o"]._pass_through


def test_a_bore_just_inside_the_tolerance_band_below_is_still_an_orifice() -> None:
    net = _net(D * (1.0 - 10.0 * ORIFICE_BORE_PIPE_RTOL))
    net.resolve_all_topology()
    o = net.elements["o"]
    assert not o._pass_through and 0.99 < o.beta < 1.0
    assert o.residual_scale_kind == "mdot"


# --- bore > pipe: refused ----------------------------------------------------


@pytest.mark.parametrize(("D_up", "D_down"), [(0.04, 0.1), (0.1, 0.04)])
def test_a_bore_larger_than_its_pipe_is_refused(D_up: float, D_down: float) -> None:
    net = _net(0.05, D_up=D_up, D_down=D_down)
    with pytest.raises(ValueError, match="larger than the pipe.*in_pipe=False"):
        NetworkSolver(net).solve()


# --- in_pipe=False: the lumped restriction, as before ------------------------


@pytest.mark.parametrize("correlation", ["fixed", "IdelchikThick", "ReaderHarrisGallagher"])
@pytest.mark.parametrize("bore", [0.03, D, 0.08])
def test_no_pipe_diameter_is_the_plenum_fed_wall_hole(correlation: str, bore: float) -> None:
    """in_pipe=False between pipes evaluates exactly as the same orifice with
    no pipe anywhere near it: beta = 0 both ways, no recovery -- the old
    beta = 0 wall hole, now said rather than inferred."""
    kw = {"correlation": correlation, "regime": "compressible"}
    lumped = _net(bore, in_pipe=False, D_up=0.04, D_down=0.1, **kw)
    plenum = _net(bore, pipes=False, **kw)
    if correlation == "IdelchikThick":
        for n in (lumped, plenum):
            n.elements["o"].plate_thickness = 0.5 * bore
    lumped.resolve_all_topology()
    plenum.resolve_all_topology()
    a, b = lumped.elements["o"], plenum.elements["o"]
    assert a.beta == a.beta_reverse == 0.0 and a.upstream_diameter is None
    assert not (a._recovery_fwd or a._recovery_rev or a._pass_through)
    for m in (0.2, -0.05):
        si = NetworkMixtureState(P=1.6e5, Pt=1.7e5, Tt=320.0, m_dot=m, Y=Y_AIR)
        so = NetworkMixtureState(P=1.2e5, Pt=1.25e5, Tt=300.0, m_dot=m, Y=Y_AIR)
        assert a.residuals(si, so) == b.residuals(si, so)


def test_in_pipe_false_keeps_the_old_wall_hole_numbers() -> None:
    """The pre-0.8 silent fallback (bore 0.06 between 0.05 pipes, 'fixed',
    compressible), recorded on origin/main 5f962b16 before #788. The explicit
    form reproduced it bit for bit there (macOS); rel 1e-12 here leaves room
    for another platform's libm."""
    net = _net(0.06, in_pipe=False, D_up=0.05, D_down=0.05)
    net.resolve_all_topology()
    o = net.elements["o"]
    si = NetworkMixtureState(P=1.6e5, Pt=1.7e5, Tt=320.0, m_dot=0.2, Y=Y_AIR)
    so = NetworkMixtureState(P=1.2e5, Pt=1.25e5, Tt=300.0, m_dot=0.2, Y=Y_AIR)
    res, jac = o.residuals(si, so)
    assert res[0] == pytest.approx(-0.40417964873021467, rel=1e-12)
    assert jac[0]["o.m_dot"] == 1.0
    assert jac[0]["a.Pt"] == pytest.approx(-5.867947674649058e-06, rel=1e-10)
    assert jac[0]["a.Tt"] == pytest.approx(0.0009462317464889918, rel=1e-10)
    assert jac[0]["b.P"] == pytest.approx(3.2780954663343763e-06, rel=1e-10)


# --- the migrated fixtures ---------------------------------------------------


def test_the_stress_random_family_declares_its_lumped_restrictions() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(here, "..", "..", "scripts", "stress"))
    from networks import FAMILIES

    lumped = 0
    for seed in range(40):
        net, _ = FAMILIES["random"](seed)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            net.resolve_all_topology()
        for e in net.elements.values():
            if type(e) is OrificeElement:
                assert not e._pass_through
                lumped += not e.in_pipe
    assert lumped > 0


def _gui(in_pipe: bool | None) -> dict:
    orf = {"Cd": 0.6, "diameter": 0.05}
    if in_pipe is not None:
        orf["in_pipe"] = in_pipe

    def node(i: str, t: str, **data: object) -> dict:
        return {"id": i, "type": t, "position": {"x": 0.0, "y": 0.0}, "data": data}

    def edge(i: str, s: str, t: str) -> dict:
        return {"id": i, "source": s, "target": t, "data": {}}

    return {
        "nodes": [
            node("in", "pressure_boundary", Pt=2.0e5, Tt=300.0),
            node("c1", "channel", L=1.0, D=0.05),
            node("orf", "orifice", **orf),
            node("c2", "channel", L=1.0, D=0.05),
            node("out", "pressure_boundary", Pt=1.0e5, Tt=300.0),
        ],
        "edges": [
            edge("e1", "in", "c1"),
            edge("e2", "c1", "orf"),
            edge("e3", "orf", "c2"),
            edge("e4", "c2", "out"),
        ],
    }


@pytest.mark.parametrize(("in_pipe", "expect"), [(None, 1.0), (True, 1.0), (False, 0.0)])
def test_the_gui_carries_in_pipe_and_previews_the_pass_through(
    in_pipe: bool | None, expect: float
) -> None:
    net = build_network_from_schema(NetworkGraphSchema(**_gui(in_pipe)))
    assert net.elements["orf"].in_pipe is (in_pipe is not False)
    rep = resolve_preview(net)
    assert rep["orf"].status == "resolved"
    assert rep["orf"].fields["pass_through"].value == expect
    assert rep["orf"].fields["beta"].value == expect


def test_the_gui_previews_a_bore_above_its_pipe_as_unresolved() -> None:
    g = _gui(None)
    g["nodes"][2]["data"]["diameter"] = 0.08
    rep = resolve_preview(build_network_from_schema(NetworkGraphSchema(**g)))
    assert rep["orf"].status == "unresolved"
    assert "in_pipe=False" in rep["orf"].reason
