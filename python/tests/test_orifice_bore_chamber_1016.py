"""The orifice bore against a momentum chamber on its flow path (#1016).

#788's rule, with areas: the bore's own geometric area (a plate's total hole
area, never Cd A) against the area of a MomentumChamberNode it opens into or
draws from, compared on the equivalent diameter to ORIFICE_BORE_PIPE_RTOL:

* larger is refused at resolve_topology, naming the orifice, the chamber and
  both areas;
* equal is NO orifice: zero loss, the row is Pt_in - Pt_out;
* smaller is unchanged, even with Cd A close to the chamber area.

A chamber area that says nothing about the geometry is not compared: the
DEFAULT_CHAMBER_AREA placeholder, an area a merge chamber inherits from an
orifice as its main inlet, and (with in_pipe=False, #788's waiver) an area
inherited from a channel.
"""

from __future__ import annotations

import math
import os
import sys
import warnings

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    EffusionPlateElement,
    FlowNetwork,
    LosslessConnectionElement,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)
from combaero.network.components import (
    DEFAULT_CHAMBER_AREA,
    ORIFICE_BORE_PIPE_RTOL,
    NetworkMixtureState,
)
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.resolve_preview import resolve_preview
from gui.backend.schemas import NetworkGraphSchema

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
A_CH = 2.0e-3  # chamber area [m^2]
D_CH = math.sqrt(4.0 * A_CH / math.pi)


def _orifice(bore: float, **kw: object) -> OrificeElement:
    return OrificeElement(
        "o", kw.pop("a", "in"), kw.pop("b", "ch"), Cd=0.6, diameter=bore, correlation="fixed", **kw
    )


def _net(bore: float | None, *, chamber_upstream: bool = False, **kw: object) -> FlowNetwork:
    """in -> orifice -> ch (A_CH) -> r -> out, or, chamber_upstream,
    in -> r -> ch -> orifice -> out, with r a small fixed restriction that
    sets the flow. bore None: a lossless connection in the orifice's
    place."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=1.2e5, Tt=300.0, Y=Y_AIR))
    n.add_node(MomentumChamberNode("ch", area=A_CH))
    n.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    a, b = ("ch", "out") if chamber_upstream else ("in", "ch")
    ra, rb = ("in", "ch") if chamber_upstream else ("ch", "out")
    n.add_element(OrificeElement("r", ra, rb, Cd=0.6, diameter=0.01, correlation="fixed"))
    if bore is None:
        n.add_element(LosslessConnectionElement("o", a, b))
    else:
        n.add_element(_orifice(bore, a=a, b=b, **kw))
    return n


# --- larger: refused -----------------------------------------------------------


@pytest.mark.parametrize("chamber_upstream", [False, True])
@pytest.mark.parametrize("in_pipe", [True, False])
def test_a_bore_wider_than_its_chamber_is_refused(chamber_upstream: bool, in_pipe: bool) -> None:
    """Falsified: with _bore_against_chambers returning "" both solve
    (bore area 1.27x the chamber's) instead of raising."""
    bore = 1.13 * D_CH
    net = _net(bore, chamber_upstream=chamber_upstream, in_pipe=in_pipe)
    side = "draws from" if chamber_upstream else "opens into"
    A_bore = math.pi * bore**2 / 4.0
    with pytest.raises(ValueError) as exc:
        net.resolve_all_topology()
    msg = str(exc.value)
    assert f"OrificeElement 'o': bore {A_bore:.6g} m^2 is larger" in msg
    assert f"momentum chamber it {side} ('ch', A = {A_CH:.6g} m^2)" in msg


def test_cd_a_below_the_chamber_does_not_excuse_a_wider_bore() -> None:
    """The comparison is geometric (#1016's decided scope): Cd A = 0.73 A_ch
    but the bore 1.21 A_ch."""
    net = _net(1.1 * D_CH)
    with pytest.raises(ValueError, match="larger than the momentum chamber"):
        net.resolve_all_topology()


# --- equal: no orifice ---------------------------------------------------------


@pytest.mark.parametrize("bore", [D_CH, D_CH * (1.0 + 0.5 * ORIFICE_BORE_PIPE_RTOL)])
@pytest.mark.parametrize("chamber_upstream", [False, True])
def test_a_bore_equal_to_its_chamber_is_no_orifice(bore: float, chamber_upstream: bool) -> None:
    """Same flow and pressures as a lossless connection in its place.
    Falsified: with _bore_against_chambers returning "" the orifice is an
    orifice (row kind "mdot", Pt_in != Pt_out) and the flows differ."""
    net = _net(bore, chamber_upstream=chamber_upstream)
    net.resolve_all_topology()
    o = net.elements["o"]
    assert o._pass_through and o.residual_scale_kind == "p"
    assert o._pass_through_reason == "bore = chamber 'ch' area: no orifice"
    r = NetworkSolver(_net(bore, chamber_upstream=chamber_upstream)).solve()
    ref = NetworkSolver(_net(None, chamber_upstream=chamber_upstream)).solve()
    assert r["__success__"] and ref["__success__"]
    assert r["o.m_dot"] == pytest.approx(ref["o.m_dot"], rel=1e-9)
    assert r["ch.Pt"] == pytest.approx(ref["ch.Pt"], rel=1e-12)
    if chamber_upstream:
        # Outflow into PressureBoundary(coupling="auto") loses the chamber
        # exit dynamic head identically to LosslessConnectionElement (#1128).
        assert r["ch.P"] == pytest.approx(r["out.Pt"], rel=1e-12)
    else:
        assert r["in.Pt"] == pytest.approx(r["ch.Pt"], rel=1e-12)


# --- smaller: unchanged ----------------------------------------------------------


def test_a_bore_just_inside_the_band_below_is_still_an_orifice() -> None:
    net = _net(D_CH * (1.0 - 10.0 * ORIFICE_BORE_PIPE_RTOL))
    net.resolve_all_topology()
    o = net.elements["o"]
    assert not o._pass_through and o.residual_scale_kind == "mdot"
    si = NetworkMixtureState(P=1.2e5, Pt=1.2e5, Tt=300.0, m_dot=0.1, Y=Y_AIR)
    so = NetworkMixtureState(P=1.1e5, Pt=1.15e5, Tt=300.0, m_dot=0.1, Y=Y_AIR)
    _, jac = o.residuals(si, so)
    assert jac[0]["o.m_dot"] == 1.0


# --- chamber areas that are not a statement ---------------------------------------


def test_a_defaulted_chamber_area_is_not_compared() -> None:
    """Nothing sizes the chamber (orifices on both sides; a bore is never an
    area to inherit): its DEFAULT_CHAMBER_AREA is a placeholder."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=1.2e5, Tt=300.0, Y=Y_AIR))
    n.add_node(MomentumChamberNode("ch"))
    n.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    bore = 2.0 * math.sqrt(4.0 * DEFAULT_CHAMBER_AREA / math.pi)
    n.add_element(_orifice(bore, a="in", b="ch"))
    n.add_element(OrificeElement("o2", "ch", "out", Cd=0.6, diameter=bore, correlation="fixed"))
    n.resolve_all_topology()
    assert n.nodes["ch"]._area_source == "default"
    assert not n.elements["o"]._pass_through


def test_a_merge_chambers_area_inherited_from_its_orifice_is_not_compared() -> None:
    """A merge chamber's area is its main inlet's; an orifice main inlet would
    otherwise read as equal to itself and vanish."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=1.2e5, Tt=300.0, Y=Y_AIR))
    n.add_node(PressureBoundary("side", Pt=1.2e5, Tt=300.0, Y=Y_AIR))
    n.add_node(MomentumChamberNode("ch", main_inlet="o"))
    n.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    n.add_element(_orifice(0.05, a="in", b="ch"))
    n.add_element(OrificeElement("s", "side", "ch", Cd=0.6, diameter=0.01, correlation="fixed"))
    n.add_element(ChannelElement("c", "ch", "out", length=1.0, diameter=0.1))
    n.resolve_all_topology()
    assert n.nodes["ch"]._area_source == "inherited from o"
    assert not n.elements["o"]._pass_through


def test_a_lumped_restriction_is_not_compared_with_a_pipes_end() -> None:
    """in_pipe=False is #788's waiver: a chamber that took a channel's area is
    that pipe's end, so a lumped Cd*A at or above it stands."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=1.2e5, Tt=300.0, Y=Y_AIR))
    n.add_node(MomentumChamberNode("ch"))
    n.add_node(PlenumNode("p"))
    n.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    n.add_element(ChannelElement("c", "in", "ch", length=1.0, diameter=0.05))
    n.add_element(_orifice(0.06, a="ch", b="p", in_pipe=False))
    n.add_element(ChannelElement("c2", "p", "out", length=1.0, diameter=0.1))
    n.resolve_all_topology()
    assert n.nodes["ch"]._area_source == "channel 'c'"
    assert not n.elements["o"]._pass_through
    # In a pipe, the same bore is #788's refusal (the pipe), as before.
    n.elements["o"].in_pipe = True
    with pytest.raises(ValueError, match="larger than the pipe"):
        n.resolve_all_topology()


# --- the plates: total hole area ------------------------------------------------


def _panel_net(chamber_area_over_holes: float) -> tuple[FlowNetwork, float]:
    plate = EffusionPlateElement(
        "eff",
        "cool",
        "ch",
        hole_diameter=0.8e-3,
        wall_thickness=2e-3,
        pitch=6e-3,
        panel_length=0.1,
        panel_width=0.1,
        angle_deg=30.0,
    )
    A = chamber_area_over_holes * plate.total_hole_area
    n = FlowNetwork()
    n.add_node(PressureBoundary("cool", Pt=1.05e5, Tt=600.0, Y=Y_AIR))
    n.add_node(MomentumChamberNode("ch", area=A))
    n.add_node(PressureBoundary("exit", Pt=1.0e5, Tt=1500.0, Y=Y_AIR))
    n.add_element(plate)
    n.add_element(
        ChannelElement("tail", "ch", "exit", length=0.2, diameter=math.sqrt(4 * A / math.pi))
    )
    return n, plate.total_hole_area


def test_a_panel_whose_holes_outsize_its_chamber_is_refused() -> None:
    """Falsified: with _bore_against_chambers returning "" it resolves."""
    n, holes = _panel_net(0.5)
    with pytest.raises(ValueError) as exc:
        n.resolve_all_topology()
    assert f"EffusionPlateElement 'eff': total hole area {holes:.6g} m^2" in str(exc.value)
    assert "momentum chamber it opens into ('ch'" in str(exc.value)


def test_a_panel_whose_holes_equal_its_chamber_is_no_plate() -> None:
    """Falsified: with _bore_against_chambers returning "" the row is the
    plate's m - m_calc."""
    n, _ = _panel_net(1.0)
    n.resolve_all_topology()
    e = n.elements["eff"]
    assert e._pass_through and e.residual_scale_kind == "p"
    si = NetworkMixtureState(P=1.05e5, Pt=1.05e5, Tt=600.0, m_dot=0.01, Y=Y_AIR)
    so = NetworkMixtureState(P=1.0e5, Pt=1.02e5, Tt=600.0, m_dot=0.01, Y=Y_AIR)
    res, jac = e.residuals(si, so)
    assert res == [pytest.approx(0.03e5)]
    assert jac == {0: {"cool.Pt": 1.0, "ch.Pt": -1.0}}


def test_a_panel_smaller_than_its_chamber_is_unchanged() -> None:
    n, _ = _panel_net(1.5)
    n.resolve_all_topology()
    assert not n.elements["eff"]._pass_through


# --- the generator and the GUI ------------------------------------------------------


def test_the_stress_random_family_fits_its_bores_into_their_chambers() -> None:
    """No random draw: random 82 still draws its d = 0.3 m bore into the
    0.02 m^2 chamber n4 (the issue's repro), now fitted to Cd A = 0.5 A."""
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(here, "..", "..", "scripts", "stress"))
    from networks import FAMILIES

    net, _ = FAMILIES["random"](82)
    e = net.elements["e11"]
    assert e.Cd * e.area == pytest.approx(0.5 * net.nodes["n4"].area, rel=1e-12)
    for seed in range(60):
        net, _ = FAMILIES["random"](seed)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                net.resolve_all_topology()
            except ValueError as exc:
                assert "momentum chamber" not in str(exc), seed
                continue
        assert not any(getattr(e, "_pass_through", False) for e in net.elements.values())


def _gui(chamber_area: float) -> dict:
    def node(i: str, t: str, **data: object) -> dict:
        return {"id": i, "type": t, "position": {"x": 0.0, "y": 0.0}, "data": data}

    def edge(i: str, s: str, t: str) -> dict:
        return {"id": i, "source": s, "target": t, "data": {}}

    return {
        "nodes": [
            node("in", "pressure_boundary", Pt=1.2e5, Tt=300.0),
            node("orf", "orifice", Cd=0.6, diameter=D_CH),
            node("ch", "momentum_chamber", area=chamber_area),
            node("c", "channel", L=1.0, D=0.2),
            node("out", "pressure_boundary", Pt=1.0e5, Tt=300.0),
        ],
        "edges": [
            edge("e1", "in", "orf"),
            edge("e2", "orf", "ch"),
            edge("e3", "ch", "c"),
            edge("e4", "c", "out"),
        ],
    }


def test_the_gui_previews_the_chamber_verdict() -> None:
    rep = resolve_preview(build_network_from_schema(NetworkGraphSchema(**_gui(A_CH))))
    assert rep["orf"].status == "resolved"
    pt = rep["orf"].fields["pass_through"]
    assert pt.value == 1.0 and pt.source == "bore = chamber 'ch' area: no orifice"
    rep = resolve_preview(build_network_from_schema(NetworkGraphSchema(**_gui(0.5 * A_CH))))
    assert rep["orf"].status == "unresolved"
    assert "larger than the momentum chamber it opens into ('ch'" in rep["orf"].reason
