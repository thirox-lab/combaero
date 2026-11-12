"""A film covers a downstream wall: ThermalWall(film=...) (#386).

The film provider (FilmCoolingElement) owns no wall. A ThermalWall on the gas
duct downstream of its holes names it, and that side's adiabatic wall
temperature becomes T_aw - eta_bar (T_aw - T_coolant); h, the stream and the
heat booking stay the duct's. eta_bar is C++'s span average: the per-row
kernel (Baldauf or Colban), Sellers over the rows, averaged over the wall.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian
from scipy.integrate import trapezoid

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    EffusionPlateElement,
    FilmCoolingElement,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)

X_AIR = cb.species.dry_air()
Y_AIR = cb.mole_to_mass(X_AIR)

# Colban configuration B geometry (AR 3.5, P/D 6) and a Baldauf-range array.
SHAPED = {"hole_type": "shaped", "lateral_pitch": 6e-3, "area_ratio": 3.5, "breakout_width": 3e-3}
CYLINDRICAL = {"hole_type": "cylindrical", "lateral_pitch": 3e-3}


# --- the C++ span average -------------------------------------------------


def _sellers_baldauf(x: float, n: int, pitch: float, M: float, DR: float) -> float:
    rows = [
        cb.film_effectiveness_baldauf_2002(x + j * pitch, M, DR, 30.0, 3.0, 0.05) for j in range(n)
    ]
    return cb.film_superposition_sellers(rows)


def _sellers_colban(x: float, n: int, pitch: float, M: float) -> float:
    rows = [cb.film_effectiveness_colban_2011(x + j * pitch, M, 3.5, 6.0, 0.5) for j in range(n)]
    return cb.film_superposition_sellers(rows)


def test_a_point_span_is_sellers_over_the_rows() -> None:
    b = cb.film_span_average_baldauf_2002(3, 10.0, 4.0, 4.0, 0.8, 1.8, 30.0, 3.0, 0.05)
    c = cb.film_span_average_colban_2011(3, 10.0, 4.0, 4.0, 0.8, 3.5, 6.0, 0.5)
    assert b.eta == pytest.approx(_sellers_baldauf(4.0, 3, 10.0, 0.8, 1.8), rel=1e-14)
    assert c.eta == pytest.approx(_sellers_colban(4.0, 3, 10.0, 0.8), rel=1e-14)


def test_one_shaped_row_at_its_breakout_is_t_over_p() -> None:
    """Colban's built-in boundary condition survives the composition."""
    c = cb.film_span_average_colban_2011(1, 0.0, 0.0, 0.0, 1.0, 3.5, 6.0, 0.5)
    assert c.eta == pytest.approx(0.5, rel=1e-14)


@pytest.mark.parametrize("n, x0, x1", [(1, 0.0, 100.0), (3, 0.0, 40.0), (5, 5.0, 60.0)])
def test_the_span_average_is_the_integral(n: int, x0: float, x1: float) -> None:
    """Against a dense trapezoid, graded towards the steep near-hole end."""
    xs = np.concatenate([[x0], x0 + (x1 - x0) * np.geomspace(1e-7, 1.0, 200001)])
    xb = np.maximum(xs, 1e-9)  # Baldauf is defined for x > 0
    ref_b = trapezoid([_sellers_baldauf(x, n, 10.0, 0.8, 1.8) for x in xb], xs) / (x1 - x0)
    ref_c = trapezoid([_sellers_colban(x, n, 10.0, 0.8) for x in xs], xs) / (x1 - x0)
    b = cb.film_span_average_baldauf_2002(n, 10.0, x0, x1, 0.8, 1.8, 30.0, 3.0, 0.05)
    c = cb.film_span_average_colban_2011(n, 10.0, x0, x1, 0.8, 3.5, 6.0, 0.5)
    assert b.eta == pytest.approx(ref_b, rel=1e-7)
    assert c.eta == pytest.approx(ref_c, rel=1e-7)


def test_the_span_average_derivatives_are_exact() -> None:
    h = 1e-6

    def b(M: float, DR: float) -> float:
        return cb.film_span_average_baldauf_2002(3, 10.0, 2.0, 50.0, M, DR, 30.0, 3.0, 0.05).eta

    def c(M: float) -> float:
        return cb.film_span_average_colban_2011(3, 10.0, 2.0, 50.0, M, 3.5, 6.0, 0.5).eta

    rb = cb.film_span_average_baldauf_2002(3, 10.0, 2.0, 50.0, 0.8, 1.8, 30.0, 3.0, 0.05)
    rc = cb.film_span_average_colban_2011(3, 10.0, 2.0, 50.0, 0.8, 3.5, 6.0, 0.5)
    assert rb.deta_dM == pytest.approx((b(0.8 + h, 1.8) - b(0.8 - h, 1.8)) / (2 * h), rel=1e-7)
    assert rb.deta_dP == pytest.approx((b(0.8, 1.8 + h) - b(0.8, 1.8 - h)) / (2 * h), rel=1e-6)
    assert rc.deta_dM == pytest.approx((c(0.8 + h) - c(0.8 - h)) / (2 * h), rel=1e-7)
    assert rc.deta_dP == 0.0


def test_the_span_must_lie_downstream_of_the_last_row() -> None:
    with pytest.raises(ValueError, match="x_start_over_D <= x_end_over_D"):
        cb.film_span_average_colban_2011(1, 0.0, 5.0, 2.0, 1.0, 3.5, 6.0, 0.5)
    with pytest.raises(ValueError, match="row_pitch_over_D > 0"):
        cb.film_span_average_colban_2011(2, 0.0, 0.0, 2.0, 1.0, 3.5, 6.0, 0.5)


# --- the network ----------------------------------------------------------


def _net(
    film: bool = True,
    holes: dict | None = None,
    cool_pt: float = 2.3e5,
    span: tuple[float, float] = (0.0, 0.05),
    n_rows: int = 2,
    mixed_supply: bool = False,
    panel: dict | None = None,
    hot_mix: bool = False,
    discharge: str = "g",
) -> FlowNetwork:
    """Hot gas through a plenum the film discharges into, then a walled duct
    whose backside is a cold duct. ``mixed_supply`` feeds the holes from a
    plenum mixing two coolant streams, so the coolant temperature is an
    unknown the wall's Jacobian must reach. ``panel`` makes the provider an
    EffusionPlateElement with those keywords instead of film holes.
    ``hot_mix`` adds a second, cooler hot-gas inflow into "g", so the
    undiluted gas the film sees is itself a mix of unknown flows (#694);
    ``discharge`` is the node the film's holes discharge into."""
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
            OrificeElement("cin", "cool", "supply", Cd=0.8, diameter=0.004, correlation="fixed")
        )
        g.add_element(
            OrificeElement("c2", "cool2", "supply", Cd=0.8, diameter=0.01, correlation="fixed")
        )
        supply = "supply"
    else:
        g.add_node(PressureBoundary("cool", Pt=cool_pt, Tt=700.0, Y=Y_AIR))
        supply = "cool"
    g.add_element(OrificeElement("hin", "hot_in", "g", Cd=0.8, diameter=0.06, correlation="fixed"))
    if hot_mix:
        g.add_node(PressureBoundary("hot2", Pt=2.15e5, Tt=1100.0, Y=Y_AIR))
        g.add_element(
            OrificeElement("hin2", "hot2", "g", Cd=0.8, diameter=0.03, correlation="fixed")
        )
    if panel is not None:
        g.add_element(EffusionPlateElement("film", supply, discharge, **panel))
    else:
        g.add_element(
            FilmCoolingElement(
                "film",
                supply,
                discharge,
                hole_diameter=1e-3,
                holes_per_row=30,
                wall_thickness=2e-3,
                n_rows=n_rows,
                row_pitch=0.01 if n_rows > 1 else None,
                **(holes or SHAPED),
            )
        )
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
    kw = {"film": "film", "film_start": span[0], "film_end": span[1]} if film else {}
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
    """RESIDUAL-ROW-SCALING: compare J_ij |x_j|, relative to each row's largest."""
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()
    fd = fd_jacobian(s, x)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    a, b = J * cols, fd * cols
    scale = np.maximum(np.abs(b).max(axis=1, keepdims=True), 1e-300)
    return float((np.abs(a - b) / scale).max())


@pytest.mark.parametrize("holes", [SHAPED, CYLINDRICAL], ids=["shaped", "cylindrical"])
def test_the_film_lowers_only_the_covered_sides_driving_temperature(holes: dict) -> None:
    _, bare = _solve(_net(film=False, holes=holes))
    _, r = _solve(_net(holes=holes))
    eta, T_c, T_aw = r["w.eta_film"], r["w.T_coolant_film"], r["w.T_aw_unfilmed"]
    assert 0.0 < eta < 1.0
    assert r["w.T_aw_a"] == pytest.approx(T_aw - eta * (T_aw - T_c), rel=1e-12)
    # Nothing else on either side moves: same h, same backside, less heat.
    assert r["w.h_a"] == pytest.approx(bare["w.h_a"], rel=1e-3)
    assert r["w.T_aw_b"] == pytest.approx(bare["w.T_aw_b"], rel=1e-6)
    assert r["w.Q"] < bare["w.Q"]
    assert "w.eta_film" not in bare


def test_a_fully_effective_film_puts_the_wall_at_the_coolant_temperature() -> None:
    """One shaped row whose breakout spans its pitch, read at the breakout:
    eta = t/P = 1, so T_aw,f = T_coolant exactly."""
    holes = dict(SHAPED, breakout_width=SHAPED["lateral_pitch"])
    _, r = _solve(_net(holes=holes, span=(0.0, 0.0), n_rows=1))
    assert r["w.eta_film"] == pytest.approx(1.0, rel=1e-14)
    assert r["w.T_aw_a"] == pytest.approx(r["w.T_coolant_film"], rel=1e-12)


def test_the_film_decays_along_the_wall() -> None:
    near = _solve(_net(span=(0.0, 0.01)))[1]["w.eta_film"]
    far = _solve(_net(span=(0.1, 0.2)))[1]["w.eta_film"]
    assert near > far > 0.0


def test_no_forward_coolant_flow_no_film() -> None:
    """film_on_wall_side hands back the side's own result when the holes do
    not discharge."""
    g = _net()
    s, r = _solve(g)
    hot = g.elements["hot"]
    state = s._get_node_state(g.nodes["g"], np.array(r["__x_solution__"]))
    state.m_dot = 0.25
    base = hot.htc_and_T(state)
    from combaero.network.components import film_on_wall_side

    for m_c in (0.0, -1e-3):
        out = film_on_wall_side(
            g.elements["film"], hot, base, state, [(state, m_c)], 0.25, 0.0, 0.05
        )
        assert out is base


def test_the_energy_balance_still_closes() -> None:
    """The film changes T_aw only; heat stays between the two streams."""
    g = _net(mixed_supply=True)
    s, r = _solve(g)
    bnd = {n for n, o in g.nodes.items() if o.is_flow_boundary}

    def h(node: str) -> float:
        T = g.nodes[node].Tt if node in bnd else s._derived_states[node][0]
        return cb.h_mass(T, X_AIR)

    e_in = e_out = 0.0
    for eid, e in g.elements.items():
        m = r["__element_diag__"][eid]["m_dot"]
        f, t = (e.from_node, e.to_node) if m >= 0.0 else (e.to_node, e.from_node)
        m = abs(m)
        if f in bnd:
            e_in += m * h(f)
        if t in bnd:
            e_out += m * h(f)
    e_out += sum(r.get(f"{b}.Q_wall_out", 0.0) for b in bnd)
    assert abs(e_out - e_in) < 1e-9 * e_in


def test_the_jacobian_reaches_the_coolant_flow_and_temperature() -> None:
    """With the holes fed from a mixing plenum, T_aw,f depends on the film's
    flow (through M) and on the coolant temperature, both outside the wall.
    The check discriminates: dropping either relay term fails it."""
    g = _net(mixed_supply=True)
    s, r = _solve(g)
    x = np.array(r["__x_solution__"])
    assert _jac_err(s, x) < 1e-5
    assert _jac_err(s, x * (1.0 + 0.03 * np.sin(np.arange(len(x)) + 1.0))) < 1e-5

    full = s._relay_film

    def without(attr: str):
        def relay(out, wall, ch, factor, relay_, up):
            if not hasattr(ch, "film"):
                return
            stub = type("Stub", (), {})()
            stub.dT_aw_dfilm_mdot = ch.dT_aw_dfilm_mdot
            stub.dT_aw_dT_coolant = ch.dT_aw_dT_coolant
            stub.dT_aw_dT_inf = ch.dT_aw_dT_inf
            old = getattr(stub, attr)
            setattr(stub, attr, [0.0] * len(old) if isinstance(old, list) else 0.0)
            full(out, wall, stub, factor, relay_, up)

        return relay

    for attr in ("dT_aw_dfilm_mdot", "dT_aw_dT_coolant"):
        s._relay_film = without(attr)
        assert _jac_err(s, x) > 1e-5, attr
    s._relay_film = full


def test_a_wall_refuses_a_film_that_is_not_one() -> None:
    g = _net(film=False)
    with pytest.raises(ValueError, match="not a film provider"):
        g.add_wall(ThermalWall(id="w2", element_a="hot", element_b="cold", film="ho"))
    with pytest.raises(ValueError, match="GAS side downstream"):
        g.add_wall(ThermalWall(id="w3", element_a="film", element_b="cold", film="film"))
    with pytest.raises(ValueError, match="film_start <= film_end"):
        ThermalWall(
            id="w4", element_a="hot", element_b="cold", film="film", film_start=0.2, film_end=0.1
        )


def test_the_network_round_trips_with_its_film() -> None:
    g = _net()
    again = FlowNetwork.from_dict(g.to_dict())
    w = again.walls["w"]
    assert (w.film, w.film_side, w.film_start, w.film_end) == ("film", "a", 0.0, 0.05)
    _, r0 = _solve(g)
    _, r1 = _solve(again)
    assert r1["w.T_aw_a"] == pytest.approx(r0["w.T_aw_a"], rel=1e-9)


def test_shaped_holes_need_their_exit_geometry() -> None:
    with pytest.raises(ValueError, match="area_ratio"):
        FilmCoolingElement(
            "f",
            "a",
            "b",
            hole_diameter=1e-3,
            holes_per_row=10,
            lateral_pitch=6e-3,
            wall_thickness=2e-3,
            hole_type="shaped",
        )


# --- an effusion panel as the provider (#386 B2) ---------------------------

# 3 rows of 0.8 mm holes at 30 deg on a 6 mm square pitch.
PANEL = {
    "hole_diameter": 0.8e-3,
    "wall_thickness": 2e-3,
    "pitch": 6e-3,
    "panel_length": 0.018,
    "panel_width": 0.18,
    "angle_deg": 30.0,
}


def test_an_effusion_panel_covers_the_wall_downstream_with_its_film() -> None:
    """The panel's rows, Baldauf each and Sellers over them, read past its
    last row: the same C++ span average at the panel's own geometry."""
    _, bare = _solve(_net(film=False, panel=PANEL))
    _, r = _solve(_net(panel=PANEL))
    D = PANEL["hole_diameter"]
    ref = cb.film_span_average_baldauf_2002(
        3,
        6e-3 / D,
        0.0,
        0.05 / D,
        r["w.film_blowing_ratio"],
        r["w.film_density_ratio"],
        30.0,
        6e-3 / D,
        0.05,
    )
    assert r["w.eta_film"] == pytest.approx(ref.eta, rel=1e-12)
    T_aw = r["w.T_aw_unfilmed"]
    assert r["w.T_aw_a"] == pytest.approx(T_aw - ref.eta * (T_aw - r["w.T_coolant_film"]))
    # The film is read against the gas before the panel's coolant mixed in
    # (#694), so it cools the wall below THAT gas. At this blowing ratio
    # (M ~ 3.4, jets lifted off) eta ~ 0.016 is below the coolant's mixed-in
    # share (~0.031), so the wall ends HOTTER than the unfilmed run, which
    # assumes the coolant fully mixed: the conservative side of #694.
    assert r["w.film_gas_undiluted"] == 1.0
    assert bare["w.T_aw_a"] < r["w.T_aw_a"] < T_aw


def test_the_panel_film_jacobian_reaches_the_coolant_flow_and_temperature() -> None:
    g = _net(mixed_supply=True, panel=PANEL)
    s, r = _solve(g)
    x = np.array(r["__x_solution__"])
    assert _jac_err(s, x) < 1e-5
    s._relay_film = lambda *a, **k: None
    assert _jac_err(s, x) > 1e-5


def test_a_panel_without_its_length_cannot_cover_a_wall() -> None:
    no_length = {k: v for k, v in PANEL.items() if k not in ("panel_length", "panel_width")}
    g = _net(film=False, panel={**no_length, "panel_area": 0.018 * 0.18})
    with pytest.raises(ValueError, match="needs panel_length"):
        g.add_wall(ThermalWall(id="w2", element_a="hot", element_b="cold", film="film"))


# --- the undiluted gas the film is read against (#694) ---------------------


def _own_T_aw(s: NetworkSolver, g: FlowNetwork, r: dict) -> tuple[float, float]:
    """The covered duct's own T_aw at the solution (the mixed gas), and its
    feed node's mixed T."""
    x = np.array(r["__x_solution__"])
    state = s._get_node_state(g.nodes["g"], x)
    state.m_dot = r["__element_diag__"]["hot"]["m_dot"]
    return float(g.elements["hot"].htc_and_T(state).T_aw), float(state.Tt)


def test_the_film_reads_the_gas_before_its_coolant_mixed_in() -> None:
    """The holes discharge into "g", whose other inflow is the 1500 K gas:
    T_aw is shifted from the mixed gas to it, and DR takes it."""
    g = _net()
    s, r = _solve(g)
    own, T_mix = _own_T_aw(s, g, r)
    assert r["w.film_gas_undiluted"] == 1.0
    assert T_mix < 1500.0 - 10.0
    assert r["w.T_aw_unfilmed"] == pytest.approx(own + 1500.0 - T_mix, rel=1e-12)
    T_c = r["w.T_coolant_film"]
    # DR takes the undiluted gas at its STATIC T where the jets meet it: the
    # duct closure at (T_inf, Pt, G_gas) (#782 audit), below T_inf / T_c.
    state = s._get_node_state(g.nodes["g"], np.array(r["__x_solution__"]))
    G = abs(r["__element_diag__"]["hot"]["m_dot"]) / g.elements["hot"].area
    T_gs = float(cb.duct_htc_static_state(1500.0, state.Pt, G, state.X).T)
    assert r["w.film_density_ratio"] == pytest.approx(T_gs / T_c, rel=1e-9)
    assert r["w.film_density_ratio"] < 1500.0 / T_c


def test_a_side_not_fed_by_the_discharge_node_reads_its_own_state() -> None:
    """Holes discharging elsewhere ("hp", downstream of the duct): the duct's
    gas does not hold their coolant, and the film reads its own state."""
    g = _net(discharge="hp")
    s, r = _solve(g)
    own, _ = _own_T_aw(s, g, r)
    assert r["w.film_gas_undiluted"] == 0.0
    assert r["w.T_aw_unfilmed"] == pytest.approx(own, rel=1e-12)


def test_the_jacobian_reaches_the_undiluted_gas() -> None:
    """With two hot inflows into "g", T_inf is a mix of unknown flows; its
    relay is part of the Jacobian, and dropping it fails the check."""
    g = _net(mixed_supply=True, hot_mix=True)
    s, r = _solve(g)
    assert r["w.film_gas_undiluted"] == 1.0
    x = np.array(r["__x_solution__"])
    assert _jac_err(s, x) < 1e-5
    assert _jac_err(s, x * (1.0 + 0.03 * np.sin(np.arange(len(x)) + 1.0))) < 1e-5

    full = s._relay_film

    def without_inf(out, wall, ch, factor, relay_, up):
        if not hasattr(ch, "film"):
            return
        stub = type("Stub", (), {})()
        stub.dT_aw_dfilm_mdot = ch.dT_aw_dfilm_mdot
        stub.dT_aw_dT_coolant = ch.dT_aw_dT_coolant
        stub.dT_aw_dT_inf = 0.0
        full(out, wall, stub, factor, relay_, up)

    s._relay_film = without_inf
    assert _jac_err(s, x) > 1e-5
    s._relay_film = full
