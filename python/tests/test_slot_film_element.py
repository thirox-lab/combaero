"""SlotFilmElement: a 2-D film slot as a film provider (#714).

The kernel is Marek & Tacina (1975) D-7958 Eqs. (4)-(5) (test_slot_film.py
scores it on the report's Table I); here the span average a wall sees, and
the element in a network: its blowing ratio is the slot's mass flux over the
gas duct's -- for D-7958's rig that is the report's (m_s/m_H)(H/S) -- and the
wall Jacobian reaches the slot's flow and coolant temperature.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FlowNetwork,
    HoleCd,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    SlotFilmElement,
    ThermalWall,
    WallLayer,
)
from combaero.network.components import SLOT_FILM_DEFAULT_TU

X_AIR = cb.species.dry_air()
Y_AIR = cb.mole_to_mass(X_AIR)
S, W = 2.0e-3, 0.25  # slot height and span [m]


# --- the C++ span average -------------------------------------------------


def test_a_point_span_is_the_kernel() -> None:
    f = cb.slot_film_span_average_marek_tacina_1975(12.0, 12.0, 0.8, 0.14)
    assert f.eta == pytest.approx(
        cb.slot_film_effectiveness_marek_tacina_1975(12.0, 0.8, 0.14), rel=1e-14
    )
    assert f.deta_dP == 0.0


@pytest.mark.parametrize(("x0", "x1", "M", "Tu"), [(0.0, 30.0, 1.0, 0.14), (5.0, 33.0, 0.5, 0.35)])
def test_the_span_average_is_the_integral(x0: float, x1: float, M: float, Tu: float) -> None:
    """Eq. (4) integrates in closed form: ln(1 + a x) / a over the span."""
    a = 0.25 * Tu / M  # Tu a fraction (#773)
    exact = (math.log1p(a * x1) - math.log1p(a * x0)) / (a * (x1 - x0))
    f = cb.slot_film_span_average_marek_tacina_1975(x0, x1, M, Tu)
    assert f.eta == pytest.approx(exact, rel=1e-9)


def test_the_span_average_derivative_is_exact() -> None:
    def eta(M: float) -> float:
        return cb.slot_film_span_average_marek_tacina_1975(0.0, 30.0, M, 0.14).eta

    h = 1e-6
    d = cb.slot_film_span_average_marek_tacina_1975(0.0, 30.0, 0.9, 0.14).deta_dM
    assert d == pytest.approx((eta(0.9 + h) - eta(0.9 - h)) / (2 * h), rel=1e-7)


# --- the element ----------------------------------------------------------


def test_the_element_refuses_what_is_not_a_slot() -> None:
    for kw in (
        {"slot_height": 0.0},
        {"span": 1e-3},
        {"lip_thickness": 0.0},
        {"turbulence_intensity": -1.0},
        {"turbulence_intensity": 14.0},  # a percent passed as a fraction (#773)
    ):
        args = {"slot_height": S, "span": W, "lip_thickness": 1e-3, **kw}
        with pytest.raises(ValueError):
            SlotFilmElement("s", "a", "b", **args)
    with pytest.raises(ValueError, match="NORMED metering"):
        SlotFilmElement(
            "s", "a", "b", S, W, 1e-3, hole_cd=HoleCd(correlation="ReaderHarrisGallagher")
        ).validate()


def test_the_hole_reynolds_number_is_the_slots() -> None:
    """The discharge correlations take 4 m/(pi d mu) per hole; with d the
    hydraulic diameter and the equivalent hole count, that is m Dh/(A mu)."""
    e = SlotFilmElement("s", "a", "b", S, W, 1e-3)
    m, mu = 0.02, 1.8e-5
    per_hole = 4.0 * (m / e._hole_count()) / (math.pi * e.hydraulic_diameter * mu)
    assert per_hole == pytest.approx(m * e.hydraulic_diameter / (S * W * mu), rel=1e-12)


def test_the_default_turbulence_is_the_documented_example() -> None:
    # A fraction, as on every element (#773): 0.14 = D-7958's 14%.
    assert SLOT_FILM_DEFAULT_TU == 0.14
    assert SlotFilmElement("s", "a", "b", S, W, 1e-3).turbulence_intensity == 0.14


def _net(span: tuple[float, float] = (0.0, 0.06), mixed_supply: bool = False) -> FlowNetwork:
    """Hot gas through a plenum the slot discharges into, then a walled duct
    whose backside is a cold duct (the #386 film-wall network)."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("hot_in", Pt=2.2e5, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("out", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(MassFlowBoundary("cold_in", m_dot=0.05, Tt=600.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cold_out", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    for n in ("g", "hp", "cp"):
        g.add_node(PlenumNode(n))
    if mixed_supply:
        g.add_node(PressureBoundary("cool", Pt=2.4e5, Tt=800.0, Y=Y_AIR))
        g.add_node(MassFlowBoundary("cool2", m_dot=0.004, Tt=450.0, Y=Y_AIR))
        g.add_node(PlenumNode("supply"))
        g.add_element(
            OrificeElement("cin", "cool", "supply", Cd=0.8, diameter=0.006, correlation="fixed")
        )
        g.add_element(
            OrificeElement("c2", "cool2", "supply", Cd=0.8, diameter=0.01, correlation="fixed")
        )
        supply = "supply"
    else:
        g.add_node(PressureBoundary("cool", Pt=2.3e5, Tt=700.0, Y=Y_AIR))
        supply = "cool"
    g.add_element(OrificeElement("hin", "hot_in", "g", Cd=0.8, diameter=0.06, correlation="fixed"))
    g.add_element(SlotFilmElement("film", supply, "g", slot_height=S, span=W, lip_thickness=1e-3))
    area = np.pi * 0.08 * 0.1
    g.add_element(
        ChannelElement(
            "hot",
            "g",
            "hp",
            length=0.1,
            diameter=0.08,
            roughness=0.0,
            surface=ConvectiveSurface(area=area),
        )
    )
    g.add_element(OrificeElement("ho", "hp", "out", Cd=0.8, diameter=0.07, correlation="fixed"))
    g.add_element(
        ChannelElement(
            "cold",
            "cold_in",
            "cp",
            length=0.1,
            diameter=0.03,
            roughness=0.0,
            surface=ConvectiveSurface(area=area),
        )
    )
    g.add_element(
        OrificeElement(
            "co", "cp", "cold_out", Cd=0.8, diameter=0.03, correlation="fixed", in_pipe=False
        )
    )
    g.add_wall(
        ThermalWall(
            id="w",
            element_a="hot",
            element_b="cold",
            layers=[WallLayer(0.002, 25.0)],
            film="film",
            film_start=span[0],
            film_end=span[1],
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


def test_the_wall_sees_the_slot_film_at_the_slots_blowing_ratio() -> None:
    """M = (m_slot / (S W)) / (m_gas / A_duct) -- for D-7958's rig, with the
    duct H x W and the slot S x W, the report's (m_s/m_H)(H/S)."""
    g = _net()
    _, r = _solve(g)
    m_slot = r["__element_diag__"]["film"]["m_dot"]
    m_gas = r["__element_diag__"]["hot"]["m_dot"]
    M = (m_slot / (S * W)) / (m_gas / (math.pi * 0.08**2 / 4.0))
    expected = cb.slot_film_span_average_marek_tacina_1975(0.0, 0.06 / S, M, 0.14).eta
    assert 0.0 < r["w.eta_film"] < 1.0
    assert r["w.eta_film"] == pytest.approx(expected, rel=1e-9)
    assert r["w.T_aw_a"] == pytest.approx(
        r["w.T_aw_unfilmed"] - r["w.eta_film"] * (r["w.T_aw_unfilmed"] - r["w.T_coolant_film"]),
        rel=1e-12,
    )


def test_the_film_decays_along_the_wall() -> None:
    near = _solve(_net(span=(0.0, 0.01)))[1]["w.eta_film"]
    far = _solve(_net(span=(0.05, 0.1)))[1]["w.eta_film"]
    assert near > far > 0.0


def _jac_err(s: NetworkSolver, x: np.ndarray) -> float:
    """RESIDUAL-ROW-SCALING: compare J_ij |x_j|, relative to each row's largest."""
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    a, b = J.toarray() * cols, fd * cols
    scale = np.maximum(np.abs(b).max(axis=1, keepdims=True), 1e-300)
    return float((np.abs(a - b) / scale).max())


def test_the_jacobian_reaches_the_slot_flow_and_coolant_temperature() -> None:
    s, r = _solve(_net(mixed_supply=True))
    x = np.array(r["__x_solution__"])
    assert _jac_err(s, x) < 1e-5
