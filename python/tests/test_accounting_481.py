"""Energy accounting semantics (#481): an injecting MassFlowBoundary mixes its
own stream in, 'fraction' scales sensible enthalpy, heat is spread over the
real flow, and heat a node cannot take is reported, not lost.
"""

from __future__ import annotations

import pytest
from _closure_check import closure

import combaero as cb
from combaero import _solver_tools
from combaero.network import (
    EnergyBoundary,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    WallNode,
)

X_AIR = cb.species.dry_air()
Y_AIR = list(cb.mole_to_mass(X_AIR))
I_CH4 = cb.species_index_from_name("CH4")
I_O2 = cb.species_index_from_name("O2")


def _pb(name: str, Pt: float, Tt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y_AIR
    return b


def _orifice(eid: str, a: str, b: str, d: float) -> OrificeElement:
    return OrificeElement(eid, a, b, Cd=0.7, diameter=d, correlation="fixed")


def _assert_closed(s: NetworkSolver, r: dict) -> None:
    c = closure(s, r)
    assert c is not None and abs(c["energy"]) < c["energy_tol"], c


def _sensible(T: float, X) -> float:
    return cb.h_mass(T, X) - cb.h_mass(cb.SENSIBLE_ENTHALPY_REF_T, X)


# --- an injecting MassFlowBoundary -------------------------------------------


def test_an_injection_mixes_its_own_stream_in() -> None:
    """50 g/s of hot methane-laden gas injected between two orifices: the
    node takes it at its own Tt and Y (it stayed at the inflow's 300 K and
    air composition)."""
    Y_inj = list(Y_AIR)
    Y_inj[I_CH4] += 0.2
    Y_inj = [y / sum(Y_inj) for y in Y_inj]
    g = FlowNetwork()
    g.add_node(_pb("a", 2e5, 300.0))
    g.add_node(_pb("out", 1e5, 300.0))
    g.add_node(MassFlowBoundary("inj", m_dot=0.05, Tt=800.0, Y=Y_inj))
    g.add_element(_orifice("e1", "a", "inj", 0.02))
    g.add_element(_orifice("e2", "inj", "out", 0.03))
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"]
    m1 = r["e1.m_dot"]
    assert r["e2.m_dot"] == pytest.approx(m1 + 0.05, rel=1e-9)
    T, Y, _ = s._derived_states["inj"]
    H = m1 * cb.h_mass(300.0, X_AIR) + 0.05 * cb.h_mass(800.0, cb.mass_to_mole(Y_inj))
    assert (m1 + 0.05) * cb.h_mass(T, cb.mass_to_mole(list(Y))) == pytest.approx(H, rel=1e-9)
    assert Y[I_CH4] == pytest.approx(0.05 * Y_inj[I_CH4] / (m1 + 0.05), rel=1e-9)
    _assert_closed(s, r)


# --- 'fraction' is a fraction of the sensible enthalpy ------------------------


def _fraction_net(T_in: float, fraction: float) -> tuple[NetworkSolver, dict]:
    g = FlowNetwork()
    g.add_node(_pb("a", 2e5, T_in))
    g.add_node(_pb("out", 1e5, 300.0))
    p = PlenumNode("p")
    p.add_energy_boundary(EnergyBoundary("loss", fraction=fraction))
    g.add_node(p)
    g.add_element(_orifice("e1", "a", "p", 0.02))
    g.add_element(_orifice("e2", "p", "out", 0.02))
    s = NetworkSolver(g)
    return s, s.solve()


@pytest.mark.parametrize("T_in", [600.0, 1500.0])
def test_a_five_percent_loss_removes_five_percent_of_the_sensible_heat(T_in: float) -> None:
    s, r = _fraction_net(T_in, -0.05)
    assert r["__success__"]
    T = s._derived_states["p"][0]
    assert _sensible(T, X_AIR) == pytest.approx(0.95 * _sensible(T_in, X_AIR), rel=1e-9)
    m = r["e1.m_dot"]
    assert r["p.Q_fraction"] == pytest.approx(-0.05 * m * _sensible(T_in, X_AIR), rel=1e-9)
    _assert_closed(s, r)


def test_a_loss_does_nothing_to_gas_at_the_reference_temperature() -> None:
    """On absolute h a '-5%' HEATED air at 298 K (formation enthalpy)."""
    s, r = _fraction_net(cb.SENSIBLE_ENTHALPY_REF_T, -0.05)
    assert s._derived_states["p"][0] == pytest.approx(cb.SENSIBLE_ENTHALPY_REF_T, abs=1e-9)


def test_a_combustor_loss_is_a_fraction_of_its_products_sensible_heat() -> None:
    Yf = [0.0] * len(Y_AIR)
    Yf[I_CH4] = 1.0

    def burn(f: float):
        st = [cb.MassStream(1.0, 600.0, 1e5, Y_AIR), cb.MassStream(0.03, 300.0, 1e5, Yf)]
        return _solver_tools.adiabatic_T_complete_and_jacobian_T_from_streams(st, 1e5, 0.0, f)

    r0, r1 = burn(0.0), burn(-0.05)
    Xb = cb.mass_to_mole(list(r1.Y_mix))
    assert _sensible(r1.Tt_mix, Xb) == pytest.approx(0.95 * _sensible(r0.Tt_mix, Xb), rel=1e-9)


@pytest.mark.parametrize("method", ["complete", "equilibrium"])
@pytest.mark.parametrize(("Q", "fraction"), [(-2e5, 0.0), (2e5, 0.0), (0.0, -0.05)])
def test_the_combustor_heat_shift_jacobian_holds(method: str, Q: float, fraction: float) -> None:
    """A post-combustion heat shift (Q or 'fraction') moves T_out through the
    burnt composition at fixed h_out as well as through T_ad: dT_out = [(1+f)
    cp_ad dT_ad + sum_j w_j dYb_j]/cp_out. Without the second term dT/dm_fuel
    was off by 110 K/(kg/s) at Q = -2e5 W and 21 at f = -0.05 (complete)."""
    Yf = [0.0] * len(Y_AIR)
    Yf[I_CH4] = 1.0
    fn = (
        _solver_tools.adiabatic_T_complete_and_jacobian_T_from_streams
        if method == "complete"
        else _solver_tools.adiabatic_T_equilibrium_and_jacobians_from_streams
    )

    def burn(m_air: float, m_fuel: float):
        st = [cb.MassStream(m_air, 700.0, 2e5, Y_AIR), cb.MassStream(m_fuel, 300.0, 2e5, Yf)]
        return fn(st, 2e5, Q, fraction)

    m_air, m_fuel, e = 1.0, 0.03, 1e-6
    r = burn(m_air, m_fuel)
    fd_air = (burn(m_air + e, m_fuel).Tt_mix - burn(m_air - e, m_fuel).Tt_mix) / (2 * e)
    fd_fuel = (burn(m_air, m_fuel + e).Tt_mix - burn(m_air, m_fuel - e).Tt_mix) / (2 * e)
    # The composition derivatives are analytic since #716 (they were forward
    # differences inside: 7e-6 complete, 2e-3 equilibrium). The rest is the
    # central difference's own error at e = 1e-6.
    tol = 2e-6
    assert r.dTt_mix_d_stream[0].d_mdot == pytest.approx(fd_air, rel=tol)
    assert r.dTt_mix_d_stream[1].d_mdot == pytest.approx(fd_fuel, rel=tol)


def _burner(fraction: float) -> tuple[FlowNetwork, str]:
    from combaero.network import CombustorNode

    net = FlowNetwork()
    net.add_node(_pb("air", 2.0e5, 700.0))
    fuel = PressureBoundary("fuel")
    Yf = [0.0] * len(Y_AIR)
    Yf[I_CH4] = 1.0
    fuel.Pt, fuel.Tt, fuel.Y = 2.05e5, 300.0, Yf
    net.add_node(fuel)
    net.add_node(_pb("out", 1.0e5, 300.0))
    comb = CombustorNode("comb", method="complete")
    comb.add_energy_boundary(EnergyBoundary("eb", fraction=fraction))
    net.add_node(comb)
    net.add_element(_orifice("oa", "air", "comb", 0.04))
    net.add_element(_orifice("of", "fuel", "comb", 0.004))
    net.add_element(_orifice("on", "comb", "out", 0.06))
    return net, "comb"


def test_a_combustor_fraction_loss_closes_the_ledger() -> None:
    """The heat a 'fraction' boundary takes from a combustor is reported as
    {node}.Q_fraction, the products' sensible heat at T_ad times f. It was
    never written for a combustor, so the ledger saw an unexplained 6.4%
    energy imbalance and demoted a correct root."""
    net, nid = _burner(-0.05)
    s = NetworkSolver(net)
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    assert abs(r["__ledger__"]["energy_closure_rel"]) < 1e-12
    comb = net.nodes[nid]
    # The combustor's solved state is its mixer output (Tt_mix, Y_mix).
    n_sp = len(Y_AIR)
    X_b = cb.mass_to_mole([r[f"{nid}.Y[{i}]"] for i in range(n_sp)])
    T_out = float(r[f"{nid}.Tt"])
    m = float(comb._total_m_dot)
    # h_out - h_ref = (1 + f)(h_ad - h_ref): the heat is f m (h_ad - h_ref)
    expected = -0.05 * m * _sensible(T_out, X_b) / 0.95
    assert r[f"{nid}.Q_fraction"] == pytest.approx(expected, rel=1e-9)
    assert expected < -1e3


# --- the mixer's Jacobian with heat and 'fraction' ---------------------------


@pytest.mark.parametrize(
    ("Q", "fraction"), [(0.0, 0.0), (3000.0, 0.0), (0.0, -0.05), (3000.0, -0.05)]
)
def test_the_mixer_jacobian_holds_with_heat_and_fraction(Q: float, fraction: float) -> None:
    """d/dm carried (1+f) on the composition term, and d/dY used the base
    enthalpy where the final one belongs: 46% off with Q, 280% with f."""
    Yb = list(Y_AIR)
    Yb[I_CH4] += 0.05
    Yb = [y / sum(Yb) for y in Yb]
    ms, Ts, Ys = [0.1, 0.05], [400.0, 900.0], [Y_AIR, Yb]

    def mix(ms_, Ts_, Ys_):
        st = [cb.MassStream(m, T, 1e5, Y) for m, T, Y in zip(ms_, Ts_, Ys_, strict=True)]
        return _solver_tools.mixer_from_streams_and_jacobians(st, Q=Q, fraction=fraction)

    r = mix(ms, Ts, Ys)
    for i in range(2):
        mp, mm = list(ms), list(ms)
        mp[i] += 1e-7
        mm[i] -= 1e-7
        fd = (mix(mp, Ts, Ys).Tt_mix - mix(mm, Ts, Ys).Tt_mix) / 2e-7
        assert r.dTt_mix_d_stream[i].d_mdot == pytest.approx(fd, rel=1e-6)
        tp, tm = list(Ts), list(Ts)
        tp[i] += 1e-3
        tm[i] -= 1e-3
        fd = (mix(ms, tp, Ys).Tt_mix - mix(ms, tm, Ys).Tt_mix) / 2e-3
        assert r.dTt_mix_d_stream[i].d_Tt == pytest.approx(fd, rel=1e-6)
        for k in (I_CH4, I_O2):
            yp, ym = [list(y) for y in Ys], [list(y) for y in Ys]
            yp[i][k] += 1e-6
            ym[i][k] -= 1e-6
            fd = (mix(ms, Ts, yp).Tt_mix - mix(ms, Ts, ym).Tt_mix) / 2e-6
            assert r.dTt_mix_d_stream[i].d_Y[k] == pytest.approx(fd, rel=1e-6)


# --- heat over the real flow; what a stagnant node cannot take ---------------


@pytest.mark.parametrize("m", [1e-6, 1e-4, 1e-3])
def test_heat_is_spread_over_the_real_flow_down_to_a_milligram_per_second(m: float) -> None:
    """The floor was 2 g/s (20% of Q withheld at 1 g/s); it is 1 mg/s now."""
    assert cb.MIXER_HEAT_MDOT_FLOOR == 1e-6
    q = 0.5 * m * 1000.0  # a 500 K rise at cp ~ 1 kJ/kg/K
    mix = _solver_tools.mixer_from_streams_and_jacobians(
        [cb.MassStream(m, 400.0, 1e5, Y_AIR)], Q=q, fraction=0.0
    )
    dh = cb.h_mass(mix.Tt_mix, X_AIR) - cb.h_mass(400.0, X_AIR)
    assert m * dh == pytest.approx(q, rel=1e-9)


def test_heat_given_to_a_stagnant_node_is_reported_as_withheld() -> None:
    """A plenum in a dead-end branch (closed by a WallNode) with a 100 W source: no flow carries it anywhere,
    so the solver reports it and the network balance still closes."""
    g = FlowNetwork()
    g.add_node(_pb("a", 2e5, 300.0))
    g.add_node(_pb("out", 1e5, 300.0))
    g.add_node(PlenumNode("p"))
    dead = PlenumNode("dead")
    dead.add_energy_boundary(EnergyBoundary("heater", Q=100.0))
    g.add_node(dead)
    g.add_element(_orifice("e1", "a", "p", 0.02))
    g.add_element(_orifice("e2", "p", "out", 0.02))
    g.add_node(WallNode("end"))
    g.add_element(_orifice("e3", "p", "dead", 0.01))
    g.add_element(_orifice("e4", "dead", "end", 0.01))
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"]
    assert abs(r["e3.m_dot"]) < 1e-9
    assert r["dead.Q_withheld"] == pytest.approx(100.0, rel=1e-6)
    # It holds the temperature of what reaches it: the C1 floor this replaced
    # spread the 100 W over 0.5 mg/s, 1e5 K, which a [200, 5000] K clamp hid.
    assert r["dead.Tt"] == pytest.approx(r["p.Tt"], abs=1e-3)
    _assert_closed(s, r)


@pytest.mark.parametrize("Q", [0.0, 2.0e5])
@pytest.mark.parametrize("m", [0.09, 0.36, 1.0])
def test_a_flowing_heated_node_reports_exactly_none_withheld(m: float, Q: float) -> None:
    """The key is there whatever the heat (a 0 W source too), and it is 0.0:
    m * (1/m) rounds to 1 - 1 ulp for about one flow in eight (0.09 and 0.36
    kg/s here, not 1.0), which surfaced as Q_withheld ~ 1e-16 Q on one
    platform and not another -- a result key that came and went with the
    last bit of a solved flow."""
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("src", m_dot=m, Tt=400.0, Y=Y_AIR))
    g.add_node(_pb("out", 1e5, 300.0))
    p = PlenumNode("p")
    p.add_energy_boundary(EnergyBoundary("heater", Q=Q))
    g.add_node(p)
    g.add_element(_orifice("e1", "src", "p", 0.3))
    g.add_element(_orifice("e2", "p", "out", 0.3))
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"]
    assert r["e1.m_dot"] == pytest.approx(m, rel=1e-12)
    assert r["p.Q_withheld"] == 0.0
    assert r["p.Q_fraction"] == 0.0  # no fraction boundary: none applied, still reported
    _assert_closed(s, r)


def test_a_fraction_boundary_with_no_flow_reports_none_applied() -> None:
    """No flow for the fraction to be a fraction of: Q_fraction is 0.0, not
    absent. (A dead branch off a flowing node is no test: it carries ~1e-22
    kg/s, not none, so the flowing branch is taken.) The orifices' rows have a
    finite slope at zero drop since #1256, so the linear solve lands on zero
    flow to round-off (~1e-30 kg/s), not on an exact 0."""
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("src", m_dot=0.0, Tt=600.0, Y=Y_AIR))
    g.add_node(_pb("out", 1e5, 300.0))
    p = PlenumNode("p")
    p.add_energy_boundary(EnergyBoundary("loss", fraction=0.05))
    g.add_node(p)
    g.add_element(_orifice("e1", "src", "p", 0.02))
    g.add_element(_orifice("e2", "p", "out", 0.02))
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"]
    assert abs(r["e1.m_dot"]) < 1e-25
    assert r["p.Q_fraction"] is not None
    assert abs(r["p.Q_fraction"]) < 1e-20


@pytest.mark.parametrize("m", [0.0, 1e-9, 3e-7, 8e-7, 1e-6, 1e-4, -5e-7, -1e-4])
def test_the_heat_taken_up_fades_c1_from_the_floor_to_none_at_zero_flow(m: float) -> None:
    """phi = 1/|m| at or above the floor; below it the share taken up,
    |m| phi = 4 u^3 - 3 u^4 (u = |m| / m0), runs 0 -> 1 with zero slope at
    both ends, and dphi/dm is phi's derivative (analytic vs central FD)."""
    m0 = cb.MIXER_HEAT_MDOT_FLOOR
    phi, dphi = cb.mixer_heat_per_unit_flow(m)
    u = abs(m) / m0
    if u >= 1.0:
        assert phi == pytest.approx(1.0 / abs(m), rel=1e-15)
    else:
        assert abs(m) * phi == pytest.approx(4 * u**3 - 3 * u**4, abs=1e-15)
        assert phi <= 1.06 / m0
    h = 1e-12
    fd = (cb.mixer_heat_per_unit_flow(m + h)[0] - cb.mixer_heat_per_unit_flow(m - h)[0]) / (2 * h)
    # (rel 1e-5: at m0 itself the difference straddles a curvature step)
    assert dphi == pytest.approx(fd, rel=1e-5, abs=1e-6 / m0**2)
    # C1 at the floor itself
    lo = cb.mixer_heat_per_unit_flow(m0 * (1 - 1e-9))
    hi = cb.mixer_heat_per_unit_flow(m0 * (1 + 1e-9))
    assert lo[0] == pytest.approx(hi[0], rel=1e-8) and lo[1] == pytest.approx(hi[1], rel=1e-6)


def test_a_root_inside_the_heat_fade_is_refused_and_the_retry_finds_the_real_one() -> None:
    """100 kW into a driven chain: the real root flows 0.25 g/s at ~2e5 K.
    From the cold start LM lands at 0.86 mg/s instead, inside the fade below
    the floor, where the node takes up 90% of its heat -- a root the fade
    makes. It is refused by name; the warm-start retry finds the real one.
    The resistance-network seed (the default) starts on the real one."""
    g = FlowNetwork()
    g.add_node(_pb("a", 1.1e5, 300.0))
    g.add_node(_pb("b", 1.0e5, 300.0))
    p = PlenumNode("p")
    p.add_energy_boundary(EnergyBoundary("heat", Q=1e5))
    g.add_node(p)
    g.add_node(PlenumNode("q"))
    for eid, a, b in (("o1", "a", "p"), ("o2", "p", "q"), ("o3", "q", "b")):
        g.add_element(OrificeElement(eid, a, b, Cd=0.8, diameter=0.01, correlation="fixed"))
    s = NetworkSolver(g)
    s.resistance_seed = False
    with pytest.warns(UserWarning, match="taking up only part"):
        r = s.solve(auto_retry=False)
    assert r["__converged__"] and not r["__success__"]
    assert r["__outcome__"] == "inconsistent"
    assert r["__heat_faded__"] == ["p"]
    r = NetworkSolver(g).solve()
    assert r["__success__"] and r["__heat_faded__"] == []
    assert r["o1.m_dot"] > 1e-4
    r = NetworkSolver(g).solve(auto_retry=False)
    assert r["__success__"] and r["o1.m_dot"] > 1e-4


def test_the_dead_fuel_boundary_api_is_gone() -> None:
    """It stored a boundary nothing read: fuel counted only if wired through
    an element. An injecting MassFlowBoundary is the way now."""
    from combaero.network import CombustorNode

    assert not hasattr(CombustorNode("c"), "set_fuel_boundary")
