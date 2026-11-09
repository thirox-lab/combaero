"""#902: a FlowNetwork exported as a GUI file solves like the network itself."""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import pytest

import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, RibbedModel, SmoothModel
from combaero.network import (
    CasingRadiation,
    ChannelElement,
    CombustorNode,
    FlowNetwork,
    FuelBoundary,
    GasRadiation,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
    ThermalWall,
    WallLayer,
)
from combaero.network.components import OrificeElement
from combaero.network.pressure_loss import (
    ConstantFractionLoss,
    ConstantHeadLoss,
    LinearThetaFractionLoss,
)
from gui.backend.export import FLOW_HANDLES, UnsupportedForGUIError, network_to_gui
from gui.backend.runner import NetworkRunner

_REPO = Path(__file__).resolve().parents[2]
AIR = list(cb.mole_to_mass(cb.species.dry_air()))
CH4 = list(cb.species.to_mass(cb.species.pure_species("CH4")))


def _both(net: FlowNetwork) -> tuple[dict, dict, dict]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        doc = network_to_gui(net)
        direct = NetworkSolver(net).solve()
        res = NetworkRunner.from_dict(json.loads(json.dumps(doc))).solve(timeout=120.0)
    assert direct["__success__"], direct.get("__message__")
    assert res.success, res.message
    return doc, direct, res._raw


def _same(direct: dict, gui: dict, keys: list[str], rel: float = 1e-6) -> None:
    for k in keys:
        assert gui[k] == pytest.approx(direct[k], rel=rel, abs=1e-9), k


def test_handles_match_the_generated_table() -> None:
    table = json.loads(
        (_REPO / "gui" / "frontend" / "src" / "utils" / "nodeHandles.generated.json").read_text()
    )
    for kind, (out, inp) in FLOW_HANDLES.items():
        assert out in table[kind]["outlets"], kind
        # None on a kind with no inlets at all: it only supplies (a fuel boundary).
        assert inp in table[kind]["inlets"] or (inp is None and not table[kind]["inlets"]), kind
    assert "side-target" in table["momentum_chamber"]["inlets"]


def test_channel_with_surface_and_dh() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("a", Pt=2.0e5, Tt=600.0, Y=AIR))
    net.add_node(PressureBoundary("b", Pt=1.8e5, Tt=600.0, Y=AIR))
    net.add_element(
        ChannelElement(
            "c",
            "a",
            "b",
            length=0.4,
            diameter=0.03,
            Dh=0.02,
            roughness=2e-6,
            regime="compressible",
            surface=ConvectiveSurface(model=SmoothModel(), Nu_multiplier=1.3, f_multiplier=0.9),
        )
    )
    doc, d, g = _both(net)
    assert doc["nodes"][2]["data"]["Dh"] == 0.02
    _same(d, g, ["c.m_dot", "c.htc", "c.Nu"])


def test_burner_plenum_lossless_and_losses() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("air", Pt=2.5e5, Tt=390.0, Y=AIR))
    net.add_node(MassFlowBoundary("fuel", m_dot=0.002, Tt=300.0, Y=CH4))
    net.add_node(CombustorNode("burner", method="complete"))
    net.add_node(PlenumNode("plenum"))
    net.add_node(PressureBoundary("out", Pt=2.0e5, Tt=800.0, Y=AIR))
    net.add_element(
        PressureLossElement("loss", "air", "burner", correlation=ConstantFractionLoss(0.04))
    )
    net.add_element(LosslessConnectionElement("inj", "fuel", "burner"))
    net.add_element(
        PressureLossElement(
            "duct",
            "burner",
            "plenum",
            correlation=LinearThetaFractionLoss(k=0.3, xi0=0.01),
            theta_source="burner",
        )
    )
    net.add_element(
        PressureLossElement(
            "exit", "plenum", "out", correlation=ConstantHeadLoss(zeta=2.0, area=0.002), area=0.002
        )
    )
    doc, d, g = _both(net)
    kinds = {n["id"]: n["type"] for n in doc["nodes"]}
    assert kinds["inj"] == "lossless_connection" and kinds["plenum"] == "plenum"
    by_id = {e["id"]: e for e in doc["edges"]}
    assert by_id["e_plenum_exit"]["sourceHandle"] == "flow-right"
    assert by_id["e_inj_burner"]["sourceHandle"] is None
    _same(d, g, ["loss.m_dot", "exit.m_dot", "plenum.Tt", "plenum.Pt"])


def test_orifice_between_pipes() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("a", Pt=2.0e5, Tt=300.0, Y=AIR))
    net.add_node(PlenumNode("n1"))
    net.add_node(PlenumNode("n2"))
    net.add_node(PressureBoundary("b", Pt=1.95e5, Tt=300.0, Y=AIR))
    net.add_element(ChannelElement("up", "a", "n1", length=2.0, diameter=0.1))
    net.add_element(
        OrificeElement("plate", "n1", "n2", diameter=0.05, correlation="ReaderHarrisGallagher")
    )
    net.add_element(ChannelElement("down", "n2", "b", length=1.0, diameter=0.1))
    _, d, g = _both(net)
    _same(d, g, ["plate.m_dot", "plate.Cd", "plate.dP_tap"])


def test_merge_chamber_side_stream_uses_the_side_handle() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("main", Pt=1.5e5, Tt=500.0, Y=AIR))
    net.add_node(MassFlowBoundary("side", m_dot=0.05, Tt=320.0, Y=AIR))
    net.add_node(MomentumChamberNode("mix", main_inlet="m_in"))
    net.add_node(PressureBoundary("out", Pt=1.4e5, Tt=450.0, Y=AIR))
    net.add_element(ChannelElement("m_in", "main", "mix", length=0.2, diameter=0.05))
    net.add_element(ChannelElement("s_in", "side", "mix", length=0.2, diameter=0.02))
    net.add_element(ChannelElement("x", "mix", "out", length=0.2, diameter=0.05))
    doc, d, g = _both(net)
    by_id = {e["id"]: e for e in doc["edges"]}
    assert by_id["e_s_in_mix"]["targetHandle"] == "side-target"
    assert by_id["e_m_in_mix"]["targetHandle"] == "flow-target"
    _same(d, g, ["x.m_dot", "mix.Tt"], rel=1e-5)


def test_thermal_wall_with_a_material_layer() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("h_in", Pt=1.2e5, Tt=900.0, Y=AIR))
    net.add_node(PressureBoundary("h_out", Pt=1.1e5, Tt=900.0, Y=AIR))
    net.add_node(MassFlowBoundary("c_in", m_dot=0.02, Tt=350.0, Y=AIR))
    net.add_node(PressureBoundary("c_out", Pt=1.5e5, Tt=350.0, Y=AIR))
    net.add_element(
        ChannelElement(
            "hot",
            "h_in",
            "h_out",
            length=0.3,
            diameter=0.04,
            surface=ConvectiveSurface(model=SmoothModel()),
        )
    )
    net.add_element(
        ChannelElement(
            "cold",
            "c_in",
            "c_out",
            length=0.3,
            diameter=0.01,
            surface=ConvectiveSurface(model=SmoothModel()),
        )
    )
    net.add_wall(
        ThermalWall(
            "w",
            "hot",
            "cold",
            layers=[WallLayer(0.002, 20.0, material="ss310")],
            contact_area=0.01,
            orientation="opposed",
        )
    )
    doc, d, g = _both(net)
    wall = next(e for e in doc["edges"] if e["id"] == "w")
    assert wall["type"] == "thermal" and wall["data"]["layers"][0]["material"] == "ss310"
    _same(d, g, ["w.Q", "w.T_hot", "cold.m_dot"], rel=1e-5)


def test_what_has_no_gui_form_is_refused_and_named() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("a", Pt=2.0e5, Tt=300.0, Y=AIR))
    net.add_node(PressureBoundary("b", Pt=1.9e5, Tt=300.0, Y=AIR, coupling="static"))
    net.add_element(
        ChannelElement(
            "c",
            "a",
            "b",
            length=1.0,
            diameter=0.02,
            surface=ConvectiveSurface(area=0.1, model=RibbedModel()),
        )
    )
    with pytest.raises(UnsupportedForGUIError) as err:
        network_to_gui(net)
    text = str(err.value)
    assert "b: PressureBoundary coupling='static'" in text
    assert "c: RibbedModel surface" in text
    assert "c: an explicit convective area" in text


def test_a_wall_without_a_side_b_stream_is_refused() -> None:
    """#970: the GUI has no field for an adiabatic or heat-loss side yet."""
    net = FlowNetwork()
    net.add_node(PressureBoundary("a", Pt=2.0e5, Tt=600.0, Y=AIR))
    net.add_node(PressureBoundary("b", Pt=1.9e5, Tt=600.0, Y=AIR))
    net.add_element(
        ChannelElement(
            "c", "a", "b", length=1.0, diameter=0.02, surface=ConvectiveSurface(model=SmoothModel())
        )
    )
    net.add_wall(ThermalWall("w", "c", None, side_b="adiabatic", contact_area=0.01))
    with pytest.raises(UnsupportedForGUIError, match="w: a wall with no stream on side B"):
        network_to_gui(net)


def test_a_head_loss_area_must_be_the_element_area() -> None:
    """The GUI's loss element has one area for both (#902): refuse, don't approximate."""
    net = FlowNetwork()
    net.add_node(PressureBoundary("a", Pt=2.0e5, Tt=300.0, Y=AIR))
    net.add_node(PressureBoundary("b", Pt=1.9e5, Tt=300.0, Y=AIR))
    net.add_element(
        PressureLossElement("x", "a", "b", correlation=ConstantHeadLoss(zeta=2.0, area=0.002))
    )
    with pytest.raises(UnsupportedForGUIError, match="give the element area=0.002"):
        network_to_gui(net)


def test_layouts() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("a", Pt=2.0e5, Tt=300.0, Y=AIR))
    net.add_node(PressureBoundary("b", Pt=1.9e5, Tt=300.0, Y=AIR))
    net.add_element(ChannelElement("c", "a", "b", length=1.0, diameter=0.02))
    rows = {n["id"]: n["position"]["x"] for n in network_to_gui(net)["nodes"]}
    assert rows["a"] < rows["c"] < rows["b"]
    given = network_to_gui(
        net, layout={"a": (0, 0), "b": (9, 9), "c": (5, 5)}, labels={"c": "Pipe"}
    )
    assert {n["id"]: n["position"]["x"] for n in given["nodes"]} == {"a": 0, "b": 9, "c": 5}
    assert next(n for n in given["nodes"] if n["id"] == "c")["data"]["label"] == "Pipe"
    with pytest.raises(ValueError, match="no position"):
        network_to_gui(net, layout={"a": (0, 0)})


@pytest.mark.parametrize(("target", "branch"), [(1500.0, "lean"), (2000.0, "rich")])
def test_a_fuel_controller_round_trips(target: float, branch: str) -> None:
    """#892: exported as a fuel boundary naming its combustor -- not as the
    MassFlowBoundary it subclasses, which would fix the fuel flow -- and the
    GUI solve finds the same fuel flow on the same branch."""
    net = FlowNetwork()
    net.add_node(PressureBoundary("air", Pt=2.5e5, Tt=600.0, Y=AIR))
    net.add_node(
        FuelBoundary("fuel", combustor="burner", Tt_target=target, branch=branch, Tt=300.0, Y=CH4)
    )
    net.add_node(CombustorNode("burner"))
    net.add_node(PressureBoundary("out", Pt=2.3e5, Tt=target, Y=AIR))
    net.add_element(
        OrificeElement(
            "inlet", "air", "burner", diameter=0.06, Cd=0.8, correlation="fixed", in_pipe=False
        )
    )
    net.add_element(LosslessConnectionElement("inj", "fuel", "burner"))
    net.add_element(
        OrificeElement(
            "exit", "burner", "out", diameter=0.05, Cd=0.8, correlation="fixed", in_pipe=False
        )
    )
    doc, direct, gui = _both(net)
    fuel = next(n for n in doc["nodes"] if n["id"] == "fuel")
    assert fuel["type"] == "fuel_boundary"
    assert fuel["data"]["combustor_id"] == "burner"
    assert fuel["data"]["branch"] == branch
    _same(direct, gui, ["fuel.m_dot", "fuel.phi", "burner.Tt", "inlet.m_dot"])
    assert gui["burner.Tt"] == pytest.approx(target, abs=1e-3)


@pytest.mark.parametrize("eta", [{"eta_poly": 0.88}, {"eta_is": 0.85}])
def test_a_compressor_and_turbine_round_trip(eta: dict) -> None:
    """#893: exported as compressor / turbine nodes with their pressure ratio
    and efficiency; the GUI solve gives the direct solve's flow, Tt and power."""
    from combaero.network import CompressorElement, TurbineElement

    net = FlowNetwork()
    net.add_node(PressureBoundary("in", Pt=1.0e5, Tt=300.0, Y=AIR))
    net.add_node(PlenumNode("a"))
    net.add_node(PlenumNode("b"))
    net.add_node(PlenumNode("c"))
    net.add_node(PressureBoundary("out", Pt=1.6e5, Tt=300.0, Y=AIR))
    net.add_element(
        OrificeElement("o1", "in", "a", diameter=0.05, Cd=0.8, correlation="fixed", in_pipe=False)
    )
    net.add_element(CompressorElement("comp", "a", "b", 4.0, **eta))
    net.add_element(TurbineElement("turb", "b", "c", 2.0, **eta))
    net.add_element(
        OrificeElement("o2", "c", "out", diameter=0.05, Cd=0.8, correlation="fixed", in_pipe=False)
    )
    doc, direct, gui = _both(net)
    kinds = {n["id"]: n["type"] for n in doc["nodes"]}
    assert kinds["comp"] == "compressor" and kinds["turb"] == "turbine"
    comp = next(n for n in doc["nodes"] if n["id"] == "comp")["data"]
    assert comp["pressure_ratio"] == 4.0
    assert comp["efficiency_type"] == ("polytropic" if "eta_poly" in eta else "isentropic")
    _same(direct, gui, ["comp.m_dot", "b.Tt", "c.Tt", "comp.power", "turb.power"])


def test_wall_radiation_round_trips() -> None:
    """#1071: each face's radiation has a GUI form, so it is exported, not dropped."""
    net = FlowNetwork()
    net.add_node(PressureBoundary("h_in", Pt=1.2e5, Tt=1600.0, Y=AIR))
    net.add_node(PressureBoundary("h_out", Pt=1.1e5, Tt=1600.0, Y=AIR))
    net.add_node(MassFlowBoundary("c_in", m_dot=0.02, Tt=600.0, Y=AIR))
    net.add_node(PressureBoundary("c_out", Pt=1.5e5, Tt=600.0, Y=AIR))
    for eid, a, b, D in (("hot", "h_in", "h_out", 0.04), ("cold", "c_in", "c_out", 0.01)):
        net.add_element(
            ChannelElement(
                eid, a, b, length=0.3, diameter=D, surface=ConvectiveSurface(model=SmoothModel())
            )
        )
    net.add_wall(
        ThermalWall(
            "w",
            "hot",
            "cold",
            layers=[WallLayer(0.002, 20.0)],
            contact_area=0.01,
            radiation_a=GasRadiation(beam_length=0.08, fuel_air_ratio=0.03),
            radiation_b=CasingRadiation(area_ratio=0.7),
        )
    )
    doc, d, g = _both(net)
    wall = next(e for e in doc["edges"] if e["id"] == "w")["data"]
    assert wall["radiation_a"]["kind"] == "gas" and wall["radiation_a"]["beam_length"] == 0.08
    assert wall["radiation_b"] == {
        "kind": "casing",
        "area_ratio": 0.7,
        "emissivity_wall": 0.7,
        "emissivity_casing": 0.4,
        "T_casing": None,
    }
    _same(d, g, ["w.Q", "w.q_rad_a", "w.q_rad_b", "w.eps_g_a"], rel=1e-5)
