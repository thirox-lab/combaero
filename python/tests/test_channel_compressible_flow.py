"""Compressible ChannelElement in mass-flow form (#481).

The residual is m - m_calc with m_calc the Fanno flow the end states drive
(fanno_channel_flow). It saturates at the choked flow, so a choked duct
converges ON its choked flow -- the former drop-given-m form patched a drop in
past choke with a barrier and converged choked ducts 5-20% above it.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

X = cb.species.dry_air()
Y = list(cb.mole_to_mass(X))
D, L, ROUGH = 0.02, 1.0, 1e-5
A = math.pi * D * D / 4
PT, TT = 2.0e5, 300.0
M_CHOKE = cb.fanno_choked_mass_flux(PT, TT, X, L, D, ROUGH, "haaland") * A


def _pb(name: str, Pt: float, Tt: float = TT, coupling: str | None = None) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y
    if coupling:
        b.coupling = coupling
    return b


def _duct(P_down: float, coupling: str | None = None) -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(_pb("in", PT))
    g.add_node(_pb("out", P_down, coupling=coupling))
    g.add_element(
        ChannelElement(
            "c", "in", "out", length=L, diameter=D, roughness=ROUGH, regime="compressible"
        )
    )
    return g


def _jac_err(s: NetworkSolver, x: np.ndarray) -> float:
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    # RESIDUAL-ROW-SCALING (NetworkSolver._build_residual_scales). In the
    # solver's own scaled variables (J_ij |x_j|): a row mixes kg/s and
    # Pa columns, and raw row scaling would hide a wrong Pa entry entirely
    # (A dG/dPt0 ~ 5e-7 next to d/dm = 1).
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    return float(np.max(np.abs(Js - fds) / scale))


@pytest.mark.parametrize("P_down", [6.0e4, 3.0e4, 1.0e4])
def test_a_choked_duct_passes_exactly_its_choked_flow(P_down: float) -> None:
    """Exit head lost (bare exit). Was 1.05-1.07 x m_choke, and 60 kPa failed."""
    s = NetworkSolver(_duct(P_down))
    r = s.solve()
    assert r["__success__"], r["__message__"]
    assert r["c.m_dot"] == pytest.approx(M_CHOKE, rel=1e-8)
    assert r["__element_diag__"]["c"]["choked"] == 1.0
    assert r["__element_diag__"]["c"]["M_exit_duct"] == 1.0


@pytest.mark.parametrize("P_down", [1.3e5, 9.0e4, 3.0e4])
def test_a_recovered_exit_chokes_at_the_same_flow(P_down: float) -> None:
    """Total coupling caps the exit at its stagnation pressure; choked, the
    flow is still the duct's choked flow (was up to 1.2 x)."""
    r = NetworkSolver(_duct(P_down, coupling="total")).solve()
    assert r["__success__"], r["__message__"]
    assert r["c.m_dot"] == pytest.approx(M_CHOKE, rel=1e-8)


@pytest.mark.parametrize(
    ("P_down", "m_ref"),
    [(1.9e5, 0.04717), (1.6e5, 0.087336), (1.3e5, 0.104801), (1.1e5, 0.110505)],
)
def test_unchoked_flows_are_the_fanno_flows(P_down: float, m_ref: float) -> None:
    """Below choke the physics is unchanged: the same Fanno flows as before
    (values after air's viscosity correction and the exact collision
    integrals, #485: 0.1-0.3% lower than before, f(Re))."""
    r = NetworkSolver(_duct(P_down)).solve()
    assert r["__success__"]
    assert r["c.m_dot"] == pytest.approx(m_ref, abs=1e-6)
    assert r["__element_diag__"]["c"]["choked"] == 0.0


def test_an_imposed_flow_above_choke_raises_the_supply_pressure() -> None:
    """1.2 x the choked flow at 2 bar: the supply must rise until the duct
    can carry it (choked flux scales with Pt), not converge past choke."""
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("in", m_dot=1.2 * M_CHOKE, Tt=TT, Y=Y))
    g.add_node(PlenumNode("p"))
    g.add_node(_pb("out", 1.0e4))
    g.add_element(
        OrificeElement("o", "in", "p", Cd=0.95, diameter=0.06, correlation="fixed", in_pipe=False)
    )
    g.add_element(
        ChannelElement(
            "c", "p", "out", length=L, diameter=D, roughness=ROUGH, regime="compressible"
        )
    )
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"], r["__message__"]
    Pt_p = r["p.Pt"]
    T_p = s._derived_states["p"][0]
    assert cb.fanno_choked_mass_flux(Pt_p, T_p, X, L, D, ROUGH, "haaland") * A == pytest.approx(
        1.2 * M_CHOKE, rel=1e-7
    )
    assert Pt_p > PT


def test_reversed_flow_is_fed_from_the_downstream_node() -> None:
    """Declared a -> b with b the hot, high-pressure end: the duct is fed from
    b's stagnation state -- the mirror of the duct declared the other way."""

    def run(rev: bool) -> float:
        g = FlowNetwork()
        g.add_node(_pb("a", 1.5e5, 300.0))
        g.add_node(_pb("b", 2.0e5, 900.0))
        if rev:
            g.add_element(
                ChannelElement(
                    "c", "a", "b", length=L, diameter=D, roughness=ROUGH, regime="compressible"
                )
            )
        else:
            g.add_element(
                ChannelElement(
                    "c", "b", "a", length=L, diameter=D, roughness=ROUGH, regime="compressible"
                )
            )
        r = NetworkSolver(g).solve()
        assert r["__success__"]
        return r["c.m_dot"]

    m_rev, m_fwd = run(True), run(False)
    assert m_rev < 0.0
    assert -m_rev == pytest.approx(m_fwd, rel=1e-10)
    G = cb.fanno_channel_flow(2.0e5, 900.0, X, 1.5e5, False, L, D, ROUGH, "haaland").G
    assert m_fwd == pytest.approx(G * A, rel=1e-9)


@pytest.mark.parametrize(
    ("P_out", "choked"), [(1.75e5, False), (1.65e5, False), (1.61e5, True), (3.0e4, True)]
)
def test_the_global_jacobian_matches_finite_differences(P_out: float, choked: bool) -> None:
    """Between two plenums, so both end pressures are unknowns (a boundary's
    pressure is fixed and would hide d/dP_target); hot supply."""
    g = FlowNetwork()
    g.add_node(_pb("in", 2.2e5, 700.0))
    g.add_node(PlenumNode("p1"))
    g.add_node(PlenumNode("p2"))
    g.add_node(_pb("out", P_out))
    g.add_element(
        OrificeElement("o1", "in", "p1", Cd=0.8, diameter=0.06, correlation="fixed", in_pipe=False)
    )
    g.add_element(
        ChannelElement(
            "c", "p1", "p2", length=L, diameter=D, roughness=ROUGH, regime="compressible"
        )
    )
    g.add_element(
        OrificeElement("o2", "p2", "out", Cd=0.8, diameter=0.04, correlation="fixed", in_pipe=False)
    )
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"], r["__message__"]
    assert r["__element_diag__"]["c"]["choked"] == float(choked)
    assert _jac_err(s, np.array(r["__x_solution__"])) < 1e-4


def test_the_global_jacobian_through_a_reversed_duct() -> None:
    g = FlowNetwork()
    g.add_node(_pb("a", 1.5e5, 300.0))
    g.add_node(PlenumNode("p"))
    g.add_node(_pb("b", 2.0e5, 900.0))
    g.add_element(
        ChannelElement("c", "a", "p", length=L, diameter=D, roughness=ROUGH, regime="compressible")
    )
    g.add_element(
        OrificeElement("o", "p", "b", Cd=0.8, diameter=0.04, correlation="fixed", in_pipe=False)
    )
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"] and r["c.m_dot"] < 0.0
    assert _jac_err(s, np.array(r["__x_solution__"])) < 1e-4


def _mixing_net(gas: dict[str, float], P_gas: float, reversed_duct: bool) -> FlowNetwork:
    """Air and a second gas mix in a plenum that feeds the duct, so the duct's
    feeding composition moves with the two inflows. reversed_duct: every
    element is declared against its flow (out -> p, p -> air, p -> gas), so the
    duct is fed from its to_node and the mixing relay runs on reversed
    orifices."""
    Xg = [0.0] * cb.num_species()
    for k, v in gas.items():
        Xg[cb.species_index_from_name(k)] = v
    g = FlowNetwork()
    g.add_node(_pb("air", 2.0e5, 300.0))
    gas_b = PressureBoundary("gas")
    gas_b.Pt, gas_b.Tt, gas_b.Y = P_gas, 500.0, list(cb.mole_to_mass(Xg))
    g.add_node(gas_b)
    g.add_node(_pb("out", 1.0e5))
    g.add_node(PlenumNode("p"))
    for eid, src, d in (("oa", "air", 0.015), ("og", "gas", 0.012)):
        ends_o = ("p", src) if reversed_duct else (src, "p")
        g.add_element(OrificeElement(eid, *ends_o, Cd=0.8, diameter=d, correlation="fixed"))
    ends = ("out", "p") if reversed_duct else ("p", "out")
    g.add_element(
        ChannelElement("c", *ends, length=L, diameter=D, roughness=ROUGH, regime="compressible")
    )
    return g


@pytest.mark.parametrize("reversed_duct", [False, True])
@pytest.mark.parametrize(
    ("gas", "P_gas"),
    [({"CO2": 0.8, "H2O": 0.2}, 2.05e5), ({"H2": 1.0}, 2.1e5), ({"CH4": 1.0}, 1.95e5)],
    ids=["co2", "h2", "ch4"],
)
def test_the_global_jacobian_carries_the_feeding_composition(
    gas: dict[str, float], P_gas: float, reversed_duct: bool
) -> None:
    """Without d/dY of the mixing plenum the inflow columns were 1.6-2.9%
    off (the inflows set the composition, the composition sets G)."""
    s = NetworkSolver(_mixing_net(gas, P_gas, reversed_duct))
    r = s.solve()
    assert r["__success__"], r["__message__"]
    assert (r["c.m_dot"] < 0.0) == reversed_duct
    assert _jac_err(s, np.array(r["__x_solution__"])) < 1e-6


def test_only_a_moving_composition_is_differentiated() -> None:
    """The column costs two length integrals per direction, so it is computed
    only where the feeding node's composition moves, and only along the
    directions it moves in."""
    mixing = NetworkSolver(_mixing_net({"CO2": 1.0}, 2.05e5, False))
    air = {cb.species_name(i): x for i, x in enumerate(X) if x > 0.0}
    single = NetworkSolver(_mixing_net(air, 2.05e5, False))
    for s, expected in ((mixing, {"p"}), (single, set())):
        x = np.array(s.solve()["__x_solution__"])
        dirs = s._composition_directions(s._propagate_states(x), x)
        assert set(dirs) == expected
        # Two gases mix: one direction, Y_air - Y_co2 normalised.
        assert all(len(v) == 1 for v in dirs.values())
        s._residuals_and_jacobian(x)
        jac = s.network.elements["c"].residuals(
            s._get_node_state(s.network.nodes["p"], x), s._get_node_state(s.network.nodes["out"], x)
        )[1]
        has_y = any(".Y[" in k for k in jac[0])
        assert has_y == bool(expected)


def test_the_residual_row_is_scaled_as_a_mass_flow() -> None:
    """Mis-scaled by ref_p it stalled a GUI tee network (697 evaluations)."""
    c = ChannelElement("c", "a", "b", length=L, diameter=D, regime="compressible")
    assert c.residual_scale_kind == "mdot"
    i = ChannelElement("i", "a", "b", length=L, diameter=D, regime="incompressible")
    assert i.residual_scale_kind == "p"


def test_the_spatial_profile_starts_at_the_ducts_inlet_static_state() -> None:
    """It raised TypeError for every compressible channel."""
    g = _duct(1.6e5)
    s = NetworkSolver(g)
    r = s.solve()
    st = s._get_node_state(g.nodes["in"], np.array(r["__x_solution__"]))
    st.m_dot = r["c.m_dot"]
    prof = g.elements["c"].get_spatial_profile(st)
    duct = cb.fanno_duct(PT, TT, r["c.m_dot"] / A, X, L, D, ROUGH, "haaland")
    assert pytest.approx(duct.inlet.P, rel=1e-10) == prof[0].P
    assert pytest.approx(duct.exit.P, rel=1e-6) == prof[-1].P
