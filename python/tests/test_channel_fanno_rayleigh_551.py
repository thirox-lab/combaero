"""A wall-heated compressible channel carries its heat in the momentum (#551).

The compressible channel was adiabatic Fanno flow: a wall's heat only raised
the downstream node's Tt, so the duct omitted the Rayleigh change (#660
flagged it). Where the duct's receiving node is fed by this duct alone, that
node's Tt IS the duct's exit Tt, and the duct now uses the Fanno-Rayleigh
kernel with it. Into a node that mixes other streams, or a boundary, the
solver forms the duct's own exit state from its stream and its own wall heat
(#672), and the duct uses that.

Heat moves a duct's choke: cooling can now carry it to the feed's full sonic
flux, choked at its inlet, where the wall's htc state (the duct's inlet
static state, #636) had an infinite slope. duct_htc_static_state bends the
flux C1 just below sonic there.
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
from networks import FAMILIES  # noqa: E402

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
X_AIR = list(cb.species.dry_air())


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=60)
    assert r["__success__"], r["__message__"]
    return s, r


def _pair(
    mix_into_cold: bool = False, cold_exit_plenum: bool = False, warm_bleed: bool = False
) -> FlowNetwork:
    """A hot and a cold compressible duct on one wall, each into its own
    plenum; with mix_into_cold a second stream also feeds the cold plenum.
    With cold_exit_plenum the cold duct discharges into a plenum of its own,
    joined losslessly to the cold plenum. With warm_bleed a 600 K stream
    drawn by pressure also feeds the cold plenum, so its T moves with x."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("hot_in", Pt=3.0e5, Tt=1200.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cold_in", Pt=3.0e5, Tt=400.0, Y=Y_AIR))
    for n in ("hot_out", "cold_out"):
        g.add_node(PressureBoundary(n, Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    for n in ("hot_p", "cold_p"):
        g.add_node(PlenumNode(n))
    if cold_exit_plenum:
        g.add_node(PlenumNode("cold_x"))
        g.add_element(LosslessConnectionElement("cold_link", "cold_x", "cold_p"))
    for name in ("hot", "cold"):
        g.add_element(
            ChannelElement(
                name,
                f"{name}_in",
                "cold_x" if (name == "cold" and cold_exit_plenum) else f"{name}_p",
                length=0.4,
                diameter=0.01,
                roughness=0.0,
                regime="compressible",
                surface=ConvectiveSurface(area=math.pi * 0.01 * 0.4),
            )
        )
        g.add_element(
            OrificeElement(
                f"{name}_o",
                f"{name}_p",
                f"{name}_out",
                Cd=0.8,
                diameter=0.012,
                correlation="fixed",
                in_pipe=False,
            )
        )
    if warm_bleed:
        g.add_node(PressureBoundary("bleed2", Pt=2.8e5, Tt=600.0, Y=Y_AIR))
        g.add_element(
            OrificeElement("b2_o", "bleed2", "cold_p", Cd=0.8, diameter=0.004, correlation="fixed")
        )
    if mix_into_cold:
        g.add_node(MassFlowBoundary("bleed", m_dot=0.005, Tt=300.0, Y=Y_AIR))
        g.add_element(
            OrificeElement("b_o", "bleed", "cold_p", Cd=0.8, diameter=0.004, correlation="fixed")
        )
    g.add_wall(
        ThermalWall(id="w", element_a="hot", element_b="cold", layers=[WallLayer(0.001, 25.0)])
    )
    return g


def test_the_duct_flow_is_the_fanno_rayleigh_flow_of_its_end_states() -> None:
    s, r = _solve(_pair())
    d = r["__element_diag__"]
    for name in ("hot", "cold"):
        el = s.network.elements[name]
        assert d[name]["rayleigh_in_momentum"] == 1.0
        assert "rayleigh_omitted" not in d[name]
        G = cb.fanno_rayleigh_channel_flow(
            r[f"{name}_in.Pt"],
            r[f"{name}_in.Tt"],
            r[f"{name}_p.Tt"],
            X_AIR,
            r[f"{name}_p.Pt"],
            True,
            el.length,
            el.diameter,
            0.0,
            el.friction_model,
        ).G
        assert r[f"{name}.m_dot"] == pytest.approx(G * el.area, rel=1e-9)
    # The heat really went across: hot cooled, cold heated.
    assert r["hot_p.Tt"] < 1200.0 and r["cold_p.Tt"] > 400.0
    assert r["__second_law_failures__"] == []


def test_the_global_jacobian_matches_finite_differences() -> None:
    s, r = _solve(_pair())
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    # RESIDUAL-ROW-SCALING: compare J_ij |x_j| per row.
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    assert float(np.max(np.abs(Js - fds) / scale)) < 1e-4


def test_into_a_mixing_node_the_duct_derives_its_own_exit_state() -> None:
    """#672: the cold plenum also takes a bleed, so its Tt is not the cold
    duct's exit Tt. The solver now forms that exit state from the duct's own
    stream and wall heat. The answer is the one the explicit outlet plenum
    gives, where the exact path of #551 applies by construction."""
    s, r = _solve(_pair(mix_into_cold=True))
    d = r["__element_diag__"]
    assert d["hot"]["rayleigh_in_momentum"] == 1.0
    assert d["cold"]["rayleigh_in_momentum"] == 1.0
    assert "rayleigh_omitted" not in d["cold"]
    exit_T = s.network.elements["cold"]._rayleigh_exit[1]
    _, ref = _solve(_pair(mix_into_cold=True, cold_exit_plenum=True))
    assert exit_T == pytest.approx(ref["cold_x.Tt"], rel=1e-6)
    for key in ("cold.m_dot", "hot.m_dot", "cold_p.Tt", "hot_p.Tt"):
        assert r[key] == pytest.approx(ref[key], rel=1e-6), key
    # The heat is in the momentum: without it the flow differs.
    assert exit_T > r["cold_in.Tt"] + 1.0


def test_into_a_mixing_node_the_global_jacobian_matches_finite_differences() -> None:
    """The exit state's derivatives chain through the wall's heat (both
    streams' htc, the capacity rates) by the same relay as a node's."""
    s, r = _solve(_pair(mix_into_cold=True))
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    assert float(np.max(np.abs(Js - fds) / scale)) < 1e-4


@pytest.mark.parametrize("seed", [10, 18, 26, 29])
def test_the_stress_pairs_that_failed_the_second_law_now_pass(seed: int) -> None:
    """Their heated duct omitted a 4-9% Pt loss (#598 audit). Seed 10's cooled
    duct now chokes at its inlet -- the case the htc state's bend is for."""
    g, _ = FAMILIES["hotcold"](seed)
    _, r = _solve(g)
    d = r["__element_diag__"]
    assert d["cch"]["rayleigh_in_momentum"] == 1.0 and d["hch"]["rayleigh_in_momentum"] == 1.0
    assert r["__second_law_failures__"] == []
    assert abs(r["__ledger__"]["mass_closure_rel"]) < 1e-9


# --- the htc state's bend below sonic ------------------------------------------

T0, P0 = 900.0, 3.0e5
GS = cb.fanno_sonic_mass_flux(P0, T0, X_AIR)
W = cb.DUCT_HTC_SONIC_MARGIN


def test_below_the_band_it_is_the_plain_duct_state() -> None:
    for g in (0.1, 0.5, 0.99, 1.0 - 2.1 * W):
        a = cb.duct_static_state(T0, P0, g * GS, X_AIR)
        b = cb.duct_htc_static_state(T0, P0, g * GS, X_AIR)
        assert (b.T, b.P, b.M, b.dT_dG, b.G_eff) == (a.T, a.P, a.M, a.dT_dG, g * GS)


def test_through_choke_it_stays_subsonic_and_c1() -> None:
    Ms = [cb.duct_htc_static_state(T0, P0, g * GS, X_AIR).M for g in (0.999, 1.0, 1.5, 10.0)]
    assert max(Ms) < 0.99
    assert Ms[1] == Ms[2] == Ms[3]  # capped: no flow dependence beyond G*

    # C1 at both band edges: the slope dT/dG on either side of each edge
    # agrees, against its size where the band starts.
    def slope(g: float) -> float:
        return cb.duct_htc_static_state(T0, P0, g * GS, X_AIR).dT_dG

    ref = abs(slope(1.0 - 2.0 * W))
    for edge in (1.0 - 2.0 * W, 1.0):
        assert abs(slope(edge * (1 + 1e-7)) - slope(edge * (1 - 1e-7))) < 1e-3 * ref


@pytest.mark.parametrize("g", [0.9985, 0.9995, 1.2])
def test_its_partials_match_differences(g: float) -> None:
    s = cb.duct_htc_static_state(T0, P0, g * GS, X_AIR)

    def f(T: float, P: float, G: float) -> np.ndarray:
        t = cb.duct_htc_static_state(T, P, G, X_AIR)
        return np.array([t.T, t.P, t.rho, t.G_eff])

    G = g * GS
    for name, i, h in (("Tt", 0, 1e-5 * T0), ("Pt", 1, 1e-5 * P0), ("G", 2, 1e-7 * G)):
        args_p, args_m = [T0, P0, G], [T0, P0, G]
        args_p[i] += h
        args_m[i] -= h
        fd = (f(*args_p) - f(*args_m)) / (2 * h)
        got = np.array(
            [getattr(s, f"d{q}_d{name}") for q in ("T", "P", "rho")]
            + [getattr(s, f"dGeff_d{name}")]
        )
        np.testing.assert_allclose(got, fd, rtol=2e-5, atol=1e-9 * np.abs(fd).max(), err_msg=name)


# --- heat area beyond the friction-wetted area ---------------------------------


def test_a_smooth_wall_with_heat_area_beyond_its_wetted_area_is_flagged() -> None:
    """Heat over surface the friction never sees breaks Reynolds analogy; the
    stress harness's crossflow segments carried 3.9x and became cooling-driven
    pumps. At its own wetted area a duct is not flagged."""
    g = _pair()
    _, r = _solve(g)
    d = r["__element_diag__"]
    assert d["hot"]["heat_area_ratio"] == pytest.approx(1.0, rel=1e-12)
    assert d["hot"]["heat_area_exceeds_wetted"] == 0.0
    g = _pair()
    g.elements["cold"].surface.area *= 3.9
    _, r = _solve(g)
    assert r["__element_diag__"]["cold"]["heat_area_ratio"] == pytest.approx(3.9, rel=1e-12)
    assert r["__element_diag__"]["cold"]["heat_area_exceeds_wetted"] == 1.0


def test_the_impingement_crossflow_segments_solve_at_their_wetted_area() -> None:
    """Impingement seed 29: six short wall-coupled crossflow segments in
    series, each carrying its heat in the momentum."""
    g, _ = FAMILIES["impingement"](29)
    _, r = _solve(g)
    d = r["__element_diag__"]
    for i in range(6):
        assert d[f"g{i}"]["heat_area_exceeds_wetted"] == 0.0
    assert r["__second_law_failures__"] == []


def test_a_cold_solve_goes_through_the_adiabatic_momentum_solution() -> None:
    """Cold, the heated channels are first solved adiabatic and the exact
    problem then from that root; it is the root an unstaged exact solve
    finds."""
    s, staged = _solve(_pair())
    assert "adiabatic-momentum" in staged["__message__"]
    plain = NetworkSolver(_pair())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = plain.solve(timeout=60, auto_retry=False)
    assert r["__success__"], r["__message__"]
    assert "adiabatic-momentum" not in r["__message__"]
    for name in ("hot", "cold"):
        assert staged[f"{name}.m_dot"] == pytest.approx(r[f"{name}.m_dot"], rel=1e-9)
    # Nothing is left suspended.
    assert not any(getattr(el, "_rayleigh_suspended", False) for el in s.network.elements.values())


def test_the_adiabatic_seed_stage_is_not_refused_by_the_second_law() -> None:
    """Its root leaves the heat out of the momentum on purpose -- the sigma < 0
    of #551 -- and only seeds the exact stage. Refused (#598), that stage ran
    the whole retry chain: a grid example went 3.5 -> 14 s. The exact stage
    is judged in full. Hotcold seed 18: its adiabatic heated duct failed the
    second law in the #598 audit."""
    s = NetworkSolver(FAMILIES["hotcold"](18)[0])
    seeds: list[dict] = []
    inner = s.solve

    def spy(*args, **kw):
        r = inner(*args, **kw)
        if s._in_rayleigh_stage and "x0" not in kw:
            seeds.append(r)
        return r

    s.solve = spy
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=60)
    assert r["__success__"] and "adiabatic-momentum" in r["__message__"]
    assert seeds and seeds[0]["__success__"]
    assert seeds[0]["__second_law_failures__"]  # reported, not applied
    assert r["__second_law_failures__"] == []


def _pair_into_boundaries() -> FlowNetwork:
    """A hot and a cold compressible duct on one wall, each straight into a
    pressure boundary: the heat lands on the boundary (Q_wall_out)."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("hot_in", Pt=3.0e5, Tt=1200.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cold_in", Pt=3.0e5, Tt=400.0, Y=Y_AIR))
    for n in ("hot_out", "cold_out"):
        g.add_node(PressureBoundary(n, Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    for name in ("hot", "cold"):
        g.add_element(
            ChannelElement(
                name,
                f"{name}_in",
                f"{name}_out",
                length=0.4,
                diameter=0.01,
                roughness=0.0,
                regime="compressible",
                surface=ConvectiveSurface(area=math.pi * 0.01 * 0.4),
            )
        )
    g.add_wall(
        ThermalWall(id="w", element_a="hot", element_b="cold", layers=[WallLayer(0.001, 25.0)])
    )
    return g


def test_into_a_boundary_the_duct_carries_its_heat_and_its_jacobian_holds() -> None:
    """#672: a boundary's Tt is fixed, never the duct's exit Tt; the duct's
    own exit state carries the heat instead."""
    s, r = _solve(_pair_into_boundaries())
    d = r["__element_diag__"]
    for name in ("hot", "cold"):
        assert d[name]["rayleigh_in_momentum"] == 1.0
    T_hot, T_cold = (s.network.elements[n]._rayleigh_exit[1] for n in ("hot", "cold"))
    assert T_hot < 1200.0 - 1.0 and T_cold > 400.0 + 1.0
    assert r["__second_law_failures__"] == []
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    assert float(np.max(np.abs(Js - fds) / scale)) < 1e-4


def test_across_a_flow_sign_and_pt_mismatch_the_rise_is_carried_and_differentiated() -> None:
    """An iterate whose cold duct flow is reversed while its Pt still falls
    forward: the exit state sits at the duct's inlet by flow sign, and the
    kernel takes the duct's own temperature rise instead of dropping to
    adiabatic. Its Jacobian still matches finite differences there."""
    s, r = _solve(_pair(mix_into_cold=True, warm_bleed=True))
    x = np.array(r["__x_solution__"])
    i = s.unknown_names.index("cold.m_dot")
    x[i] = -0.3 * x[i]
    s._residuals(x)
    cold = s.network.elements["cold"]
    assert cold._rayleigh_exit[0] == "cold_in"  # by flow sign
    st_in = s._get_node_state(s.network.nodes["cold_in"], x)
    st_out = s._get_node_state(s.network.nodes["cold_p"], x)
    state = cold._rayleigh_exit_state("cold_p", st_out, "cold_in", st_in)
    assert state is not None and len(state[1]) == 3
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    assert float(np.max(np.abs(Js - fds) / scale)) < 1e-4


def test_a_staged_solve_restores_each_ducts_own_suspension() -> None:
    """A solve nested inside stage 1 (the wall-free seed solves the same
    elements) used to clear every duct's suspension on its way out, so the
    outer adiabatic stage carried the heat again (#672). The flag a duct had
    before a staged solve is the flag it has after."""
    g = _pair()
    s = NetworkSolver(g)
    heated = s._rayleigh_heated_elements()
    assert heated
    for el in heated:
        el._rayleigh_suspended = True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s._solve_rayleigh_staged(heated, {"timeout": 30.0})
    assert all(el._rayleigh_suspended for el in heated)


def test_with_thermal_coupling_off_nothing_is_staged() -> None:
    g = _pair()
    g.thermal_coupling_enabled = False
    assert NetworkSolver(g)._rayleigh_heated_elements() == []
