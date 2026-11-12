"""#962: the film-side heat transfer coefficient near a slot, opt-in on
SlotFilmElement (``film_htc``) for the ThermalWall side its film covers.

The kernels (C++, tests/test_slot_film_htc.cpp) are Seban (1960) Fig. 12 as
drawn and the Ballal & Lefebvre (1973) wall jet. Here: the wall side's h in
a network is the kernel's at the solved state, the default changes nothing,
the hand-back to the duct h, the diagnostics, and the Jacobian through the
slot's flow, its coolant temperature and the gas flow (central differences
of the residual).

Scoring on Seban's own points and the cross-source SME4 set:
test_slot_film_htc_validation_962.py.
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
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    SlotFilmElement,
    ThermalWall,
    WallLayer,
)

X_AIR = cb.species.dry_air()
Y_AIR = cb.mole_to_mass(X_AIR)
S, W = 2.0e-3, 0.25  # slot height and span [m]
D_HOT = 0.08


def _net(
    model: str = "seban1960",
    span: tuple[float, float] = (0.0, 0.05),
    cool_Pt: float = 2.3e5,
    supply_d: float = 0.03,
    cool2_mdot: float = 0.01,
) -> FlowNetwork:
    """Hot gas into a plenum the slot discharges into, then a walled duct
    whose backside is a cold duct. The slot's supply plenum mixes two
    coolant streams, so its temperature moves with the unknowns."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("hot_in", Pt=2.2e5, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("out", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(MassFlowBoundary("cold_in", m_dot=0.05, Tt=600.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cold_out", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cool", Pt=cool_Pt, Tt=800.0, Y=Y_AIR))
    g.add_node(MassFlowBoundary("cool2", m_dot=cool2_mdot, Tt=450.0, Y=Y_AIR))
    for n in ("g", "hp", "cp", "supply"):
        g.add_node(PlenumNode(n))
    g.add_element(
        OrificeElement("cin", "cool", "supply", Cd=0.8, diameter=supply_d, correlation="fixed")
    )
    g.add_element(
        OrificeElement("c2", "cool2", "supply", Cd=0.8, diameter=0.03, correlation="fixed")
    )
    g.add_element(OrificeElement("hin", "hot_in", "g", Cd=0.8, diameter=0.06, correlation="fixed"))
    g.add_element(
        SlotFilmElement(
            "film", "supply", "g", slot_height=S, span=W, lip_thickness=1e-3, film_htc=model
        )
    )
    area = np.pi * D_HOT * 0.1
    g.add_element(
        ChannelElement(
            "hot",
            "g",
            "hp",
            length=0.1,
            diameter=D_HOT,
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


def _jac_err(s: NetworkSolver, x: np.ndarray) -> float:
    """RESIDUAL-ROW-SCALING: compare J_ij |x_j|, relative to each row's largest."""
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    a, b = J.toarray() * cols, fd * cols
    scale = np.maximum(np.abs(b).max(axis=1, keepdims=True), 1e-300)
    return float((np.abs(a - b) / scale).max())


def _kernel(model: str) -> object:
    return {
        "seban1960": cb.slot_film_htc_seban_1960,
        "ballal_lefebvre1973": cb.slot_film_htc_ballal_lefebvre_1973,
    }[model]


@pytest.mark.parametrize("model", ["seban1960", "ballal_lefebvre1973"])
def test_the_side_h_is_the_kernels_at_the_solved_state(model: str) -> None:
    s, r = _solve(_net(model))
    diag = r["__element_diag__"]
    m_slot, m_gas = diag["film"]["m_dot"], diag["hot"]["m_dot"]
    M = (m_slot / (S * W)) / (m_gas / (math.pi * D_HOT**2 / 4.0))
    assert r["w.film_blowing_ratio"] == pytest.approx(M, rel=1e-9)
    sup = s._get_node_state(s.network.nodes["supply"], np.array(r["__x_solution__"]))
    k = _kernel(model)(
        0.0, 0.05, S, m_slot / (S * W), M, r["w.film_h_duct"], sup.Tt, sup.P, list(sup.X)
    )
    assert r["w.h_a"] == pytest.approx(k.h, rel=1e-9)
    assert r["w.film_h"] == pytest.approx(k.h, rel=1e-9)
    assert r["w.film_h_film"] == pytest.approx(k.h_film, rel=1e-9)
    # m > 1.1 here: fully the film's h, well above the duct's.
    assert M > 1.5
    assert r["w.film_h_weight_M"] == 1.0
    assert r["w.h_a"] > 2.0 * r["w.film_h_duct"]
    expected_branch = (
        cb.SLOT_FILM_BRANCH_NEAR if model == "seban1960" else cb.SLOT_FILM_BRANCH_WALL_JET
    )
    assert r["w.film_h_branch"] == expected_branch
    for key in ("film_h_Re_s", "film_h_in_box", "film_h_x_in_box", "film_h_M_in_box"):
        assert f"w.{key}" in r


def test_the_default_changes_nothing() -> None:
    """film_htc="none" (the default) reports no film-side h and keeps the
    duct's h: the same solution as before #962 (bit-identity over the stress
    families is in the PR's evidence)."""
    _, r = _solve(_net("none"))
    assert "w.film_h" not in r
    _, r1 = _solve(_net("seban1960"))
    assert r1["w.film_h_duct"] != r1["w.h_a"]
    assert r["w.Q"] != pytest.approx(r1["w.Q"], rel=1e-3)


@pytest.mark.parametrize(
    ("model", "kw"),
    [
        ("seban1960", {}),
        ("ballal_lefebvre1973", {}),
        # The span crossing the B/C blend, and past the wall jet's x/s 50.
        ("seban1960", {"span": (0.05, 0.1)}),
        ("ballal_lefebvre1973", {"span": (0.05, 0.15)}),
    ],
)
def test_the_jacobian_reaches_the_slot_flow_and_coolant_temperature(model: str, kw: dict) -> None:
    s, r = _solve(_net(model, **kw))
    x = np.array(r["__x_solution__"])
    assert _jac_err(s, x) < 1e-5


def _band_net(model: str, M_target: float) -> tuple[NetworkSolver, np.ndarray]:
    """The solved network, then x moved so m sits inside the hand-back band
    (the residual need not be zero for a Jacobian check)."""
    s, r = _solve(_net(model))
    x = np.array(r["__x_solution__"])
    i = s.unknown_names.index("film.m_dot")
    x[i] *= M_target / r["w.film_blowing_ratio"]
    return s, x


@pytest.mark.parametrize(("model", "M"), [("seban1960", 1.04), ("ballal_lefebvre1973", 1.35)])
def test_the_jacobian_inside_the_hand_back_band(model: str, M: float) -> None:
    """dh/dM enters through both the slot flow and the gas flow."""
    s, x = _band_net(model, M)
    wall = s.network.walls["w"]
    obj = s.network.elements["hot"]
    state = s._get_node_state(s.network.nodes["g"], x)
    state.m_dot = x[s.unknown_names.index("hot.m_dot")]
    s._propagate_states(x, with_relay=False)
    ch = s._wall_side_htc(wall, obj, state, x, lagged=False)
    assert 0.0 < ch.film["film_h_weight_M"] < 1.0
    assert _jac_err(s, x) < 1e-5


def test_below_unit_blowing_the_side_hands_back_to_the_duct() -> None:
    _, r = _solve(_net("seban1960", cool_Pt=2.05e5, supply_d=0.006, cool2_mdot=0.002))
    assert r["w.film_blowing_ratio"] < 1.0
    assert r["w.film_h_weight_M"] == 0.0
    assert r["w.h_a"] == r["w.film_h_duct"]
    assert r["w.film_h_M_in_box"] == 0.0


def test_refusals() -> None:
    with pytest.raises(ValueError, match="film_htc"):
        SlotFilmElement("s", "a", "b", S, W, 1e-3, film_htc="seban")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="film_end > 0"):
        _net("seban1960", span=(0.0, 0.0))
    # The point value away from the lip is allowed.
    _net("seban1960", span=(0.03, 0.03))


def _relay_err(s: NetworkSolver, x: np.ndarray) -> float:
    """Largest |relay - central difference| x |x_j| [K] of the two
    wall-heated nodes' temperatures, over every unknown."""
    relay = s._propagate_states(x, with_relay=True)

    def T(v: np.ndarray, nid: str) -> float:
        s._propagate_states(v, with_relay=False)
        return float(s._derived_states[nid][0])

    worst = 0.0
    for nid in ("hp", "cp"):
        for j in range(len(x)):
            h = 1e-6 * max(abs(x[j]), 1e-3)
            e = np.zeros_like(x)
            e[j] = h
            fd = (T(x + e, nid) - T(x - e, nid)) / (2 * h)
            an = relay[nid].get(j, {}).get("T", 0.0)
            worst = max(worst, abs(an - fd) * abs(x[j]))
    return worst


@pytest.mark.parametrize(
    ("model", "M"),
    [("seban1960", None), ("ballal_lefebvre1973", None), ("seban1960", 1.04)],
)
def test_the_relay_is_exact(model: str, M: float | None) -> None:
    """The wall-heated nodes' temperature sensitivities against central
    differences, in every unknown: the film h's slot flow, coolant
    temperature (through the slot fluid's mu, cp, k) and gas flow terms.
    Dropping the coolant-temperature term moves them by 0.05-0.1 K, the
    duct-h scaling (1 - w f_x) by 0.01-0.03 K; this holds them to 1e-5 K."""
    if M is None:
        s, r = _solve(_net(model))
        x = np.array(r["__x_solution__"])
    else:
        s, x = _band_net(model, M)
    assert _relay_err(s, x) < 1e-5
