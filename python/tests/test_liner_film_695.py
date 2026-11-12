"""An EffusionLiner as the film a downstream wall sees (#695).

The liner's panels are ONE film provider, registered under the liner's id:
each panel a group of Baldauf rows at its own M, DR and coolant temperature,
layered row by row along the gas path. One panel must be the plate's own
film exactly; equal panels must be Sellers over all their rows; the wall's
Jacobian must reach every panel's flow and supply temperature.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    EffusionLiner,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)
from combaero.network.effusion_liner import LinerFilm

Y_AIR = cb.mole_to_mass(cb.species.dry_air())

LINER = {
    "length": 0.06,
    "width": 0.15,
    "duct_height": 0.02,
    "hole_diameter": 0.8e-3,
    "wall_thickness": 2e-3,
    "pitch_x": 6e-3,
    "pitch_y": 6e-3,
    "angle_deg": 30.0,
}


def _liner(n: int = 3, **kw) -> EffusionLiner:
    return EffusionLiner("ln", n_segments=n, **{**LINER, **kw})


def _net(
    liner: EffusionLiner, film: bool = True, hot_mix: bool = False, mixed_supply: bool = False
) -> FlowNetwork:
    """Hot gas into a plenum the liner's panels discharge into, then a duct
    walled against a cold duct; the liner's backside duct runs from "cin" to
    "cout". ``hot_mix`` adds a cooler second hot inflow (the undiluted gas is
    a mix of unknown flows); ``mixed_supply`` feeds the liner from a plenum
    mixing two coolant streams, so the stations' temperatures are unknowns."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("hot_in", Pt=2.2e5, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("out", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cin", Pt=2.4e5, Tt=700.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cout", Pt=2.35e5, Tt=700.0, Y=Y_AIR))
    g.add_node(MassFlowBoundary("cold_in", m_dot=0.05, Tt=600.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cold_out", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    for n in ("g", "hp", "cp"):
        g.add_node(PlenumNode(n))
    g.add_element(OrificeElement("hin", "hot_in", "g", Cd=0.8, diameter=0.06, correlation="fixed"))
    if hot_mix:
        g.add_node(PressureBoundary("hot2", Pt=2.15e5, Tt=1100.0, Y=Y_AIR))
        g.add_element(
            OrificeElement("hin2", "hot2", "g", Cd=0.8, diameter=0.03, correlation="fixed")
        )
    supply = "cin"
    if mixed_supply:
        g.add_node(PlenumNode("csup"))
        g.add_node(MassFlowBoundary("cin2", m_dot=0.04, Tt=400.0, Y=Y_AIR))
        g.add_element(
            OrificeElement("cs1", "cin", "csup", Cd=0.8, diameter=0.02, correlation="fixed")
        )
        g.add_element(
            OrificeElement("cs2", "cin2", "csup", Cd=0.8, diameter=0.02, correlation="fixed")
        )
        supply = "csup"
    liner.add_to(g, supply, "cout", "g")
    area = np.pi * 0.08 * 0.1
    g.add_element(
        ChannelElement(
            "hot", "g", "hp", length=0.1, diameter=0.08, surface=ConvectiveSurface(area=area)
        )
    )
    g.add_element(OrificeElement("ho", "hp", "out", Cd=0.8, diameter=0.07, correlation="fixed"))
    g.add_element(
        ChannelElement(
            "cold", "cold_in", "cp", length=0.1, diameter=0.03, surface=ConvectiveSurface(area=area)
        )
    )
    g.add_element(
        OrificeElement(
            "co", "cp", "cold_out", Cd=0.8, diameter=0.03, correlation="fixed", in_pipe=False
        )
    )
    kw = {"film": "ln", "film_start": 0.0, "film_end": 0.05} if film else {}
    g.add_wall(
        ThermalWall(
            id="w", element_a="hot", element_b="cold", layers=[WallLayer(0.002, 25.0)], **kw
        )
    )
    return g


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r.get("__message__")
    return s, r


def _jac_err(s: NetworkSolver, x: np.ndarray) -> float:
    """RESIDUAL-ROW-SCALING: compare J_ij |x_j|, relative to each row's largest,
    over the rows outside the liner. The film's terms land there (the wall's
    heat reaches the ducts' downstream nodes); the liner's own rows carry its
    McGreehan-Schotsch and station chain, checked in test_effusion_liner.py
    to that file's own tolerance."""
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()
    fd = fd_jacobian(s, x)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    a, b = J * cols, fd * cols
    scale = np.maximum(np.abs(b).max(axis=1, keepdims=True), 1e-300)
    rows = [i for i, n in enumerate(s.unknown_names) if not n.startswith("ln__")]
    return float((np.abs(a - b) / scale)[rows].max())


def test_the_liner_registers_one_provider_with_a_source_per_panel() -> None:
    g = _net(_liner(3))
    p = g.film_provider("ln")
    assert isinstance(p, LinerFilm)
    assert [s.element_id for s in p.film_sources()] == ["ln__p1", "ln__p2", "ln__p3"]
    assert [s.supply_node for s in p.film_sources()] == ["ln__b1", "ln__b2", "ln__b3"]
    assert p.film_discharge_node() == "g"
    assert p.rows_per_panel() == round(0.02 / 6e-3)


def test_one_panel_is_the_plates_own_film() -> None:
    g = _net(_liner(1))
    p = g.film_provider("ln")
    plate = g.elements["ln__p1"]
    a = p.film_span_average_groups(0.0, 0.05, [2.1], [1.9])
    b = plate.film_span_average(0.0, 0.05, 2.1, 1.9)
    assert a.eta == pytest.approx(b.eta, rel=1e-14)
    assert a.weight[0] == pytest.approx(b.eta, rel=1e-14)
    assert a.deta_dM[0] == pytest.approx(b.deta_dM, rel=1e-12)


def test_equal_panels_are_one_plate_with_all_the_rows() -> None:
    p = _net(_liner(3)).film_provider("ln")
    a = p.film_span_average_groups(0.0, 0.05, [2.1] * 3, [1.9] * 3)
    D = LINER["hole_diameter"]
    b = cb.film_span_average_baldauf_2002(
        3 * p.rows_per_panel(), 6e-3 / D, 0.0, 0.05 / D, 2.1, 1.9, 30.0, 6e-3 / D, 0.05
    )
    assert a.eta == pytest.approx(b.eta, rel=1e-13)


def test_the_gas_direction_orders_the_panels() -> None:
    w = _net(_liner(3)).film_provider("ln")
    a = _net(_liner(3, film_direction="against_coolant")).film_provider("ln")
    assert [s.element_id for s in a.film_sources()] == ["ln__p3", "ln__p2", "ln__p1"]
    M = [1.0, 2.0, 3.0]
    fw = w.film_span_average_groups(0.0, 0.05, M, [1.9] * 3)
    fa = a.film_span_average_groups(0.0, 0.05, M[::-1], [1.9] * 3)
    # Panel by panel the same physics, reversed along the gas path.
    assert fa.eta != pytest.approx(fw.eta, rel=1e-6)


def test_the_wall_downstream_sees_the_liners_film() -> None:
    _, bare = _solve(_net(_liner(3), film=False))
    _, r = _solve(_net(_liner(3)))
    eta, T_c, T_ref = r["w.eta_film"], r["w.T_coolant_film"], r["w.T_aw_unfilmed"]
    assert 0.0 < eta < 1.0
    assert r["w.film_gas_undiluted"] == 1.0
    # sum_g w_g T_c,g = eta T_coolant_film, by its definition.
    assert r["w.T_aw_a"] == pytest.approx(T_ref - eta * (T_ref - T_c), rel=1e-12)
    assert 600.0 < T_c < 800.0
    assert r["w.T_aw_b"] == pytest.approx(bare["w.T_aw_b"], rel=1e-6)


def test_the_jacobian_reaches_every_panels_flow_and_supply_temperature() -> None:
    """The panels' flows and their station temperatures are unknowns outside
    the wall; dropping any of the film's relay terms fails the check."""
    # 2 mm holes: above McGreehan-Schotsch's Re floor (1e4), where the panel's
    # own derivative is exact (test_effusion_liner.py).
    g = _net(
        _liner(3, hole_diameter=2e-3, pitch_x=10e-3, pitch_y=10e-3), hot_mix=True, mixed_supply=True
    )
    s, r = _solve(g)
    assert r["w.film_gas_undiluted"] == 1.0
    x = np.array(r["__x_solution__"])
    assert _jac_err(s, x) < 1e-5
    assert _jac_err(s, x * (1.0 + 0.03 * np.sin(np.arange(len(x)) + 1.0))) < 1e-5

    full = s._relay_film

    def without(attr: str):
        def relay(out, wall, ch, factor, relay_, up):
            if not hasattr(ch, "film"):
                return
            stub = type("Stub", (), {})()
            for a in ("dT_aw_dfilm_mdot", "dT_aw_dT_coolant", "dT_aw_dT_inf"):
                v = getattr(ch, a)
                if a == attr:
                    v = [0.0] * len(v) if isinstance(v, list) else 0.0
                setattr(stub, a, v)
            full(out, wall, stub, factor, relay_, up)

        return relay

    for attr in ("dT_aw_dfilm_mdot", "dT_aw_dT_coolant", "dT_aw_dT_inf"):
        s._relay_film = without(attr)
        assert _jac_err(s, x) > 1e-5, attr
    s._relay_film = full


def test_a_wall_cannot_film_the_liners_own_backside() -> None:
    g = _net(_liner(3), film=False)
    g.add_element(
        ChannelElement(
            "bk", "ln__b2", "cp", length=0.05, diameter=0.02, surface=ConvectiveSurface(area=0.01)
        )
    )
    with pytest.raises(ValueError, match="coolant supply"):
        g.add_wall(ThermalWall(id="w2", element_a="ln__b2", element_b="cold", film="ln"))
