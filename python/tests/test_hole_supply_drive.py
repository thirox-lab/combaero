"""A hole's drive pair (#1087): the supply drive of a duct-fed hole
(``HoleCd.supply_drive``) and the direction of a cooling hole into a flowing
chamber.

Supply side. Every hole drives from the upstream TOTAL to the downstream
static (the owner's rule, 2026-11-12); its Cd model owns how much supply
head it recovers. A duct-fed EffusionPlateElement's supply node is a duct
station holding the duct's STATIC pressure, so ``supply_drive="auto"``
rebuilds the duct's stagnation pressure there (C++ hole_supply_drive, the
momentum-chamber closure). Every base Cd is referenced to the partly
recovered drive P1 + f q1 (``HoleCd.supply_recovery``, default f = 0 at every
angle and for every base) and converted onto it exactly (C++
hole_supply_recovery_cd_ratio): at f = 0 the flow is the static-referenced
model's own, at the station's static temperature (McGreehan & Schotsch's
Eq. 17 is static-referenced by its authors); supply_recovery=1.0 recovers the
whole head. "static" is the v0.11.1 drive.
Design record: docs/archive/hole_drive_convention.md.

Exit side. A film, effusion or slot hole into a MomentumChamberNode takes its
direction from its supply drive against the chamber's STATIC pressure, the
pressure the jet discharges against, and a reversed jet draws from that
static. Before, the direction came from Pt - Pt, which reversed a hole that
blows wherever P_static,chamber < Pt_coolant < Pt_chamber.

Every Jacobian check here is central differences with a pressure step sized
on the local pressure DROP, not on |P| (#1107): m ~ sqrt(dp), so a step of
1e-6 |P| on a small drop is truncation, not a Jacobian error.
"""

from __future__ import annotations

import copy
import math
import warnings

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian, scaled_error

import combaero as cb
from combaero.network import (
    ChannelElement,
    CrossflowSegmentElement,
    EffusionLiner,
    EffusionPlateElement,
    FilmCoolingElement,
    FlowNetwork,
    HoleCd,
    ImpingementPlateElement,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    SlotFilmElement,
)
from combaero.network.components import _moving_static_T, _moving_total_P
from combaero.network.hole_cd import HOLE_CD_DIAGNOSTIC_KEYS

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
X_AIR = list(cb.species.dry_air())
P_ATM = 101325.0
A_DUCT = 1.2e-3
D_HOLE = 0.01


# --- the kernel --------------------------------------------------------------


@pytest.mark.parametrize("incompressible", [False, True])
def test_the_total_drive_is_the_momentum_chamber_closures_total(incompressible: bool) -> None:
    # Static + the closure's dynamic head: what a MomentumChamberNode at the
    # station's static state and flow would hold as its Pt.
    P, Tt, m = 1.2e5, 450.0, 0.35
    tot = cb.hole_supply_drive(cb.HoleSupplyDrive.Total, P, Tt, Y_AIR, m, A_DUCT, incompressible)
    sta = cb.hole_supply_drive(cb.HoleSupplyDrive.Static, P, Tt, Y_AIR, m, A_DUCT, incompressible)
    c = cb._core.momentum_chamber_residual_and_jacobian(
        P, 0.0, m, Tt, Y_AIR, A_DUCT, incompressible
    )
    assert tot.q_dyn == sta.q_dyn == pytest.approx(c.q_dyn, rel=1e-12) and c.q_dyn > 0.0
    assert tot.P_drive == -c.residual  # the closure's value, bounded branch included
    assert tot.P_drive == pytest.approx(P + c.q_dyn, rel=1e-15)
    assert sta.P_drive == P and sta.dP_dP == 1.0 and sta.dP_dTt == 0.0 and sta.dP_dm_dot == 0.0
    if not incompressible:
        assert tot.P_drive == pytest.approx(_moving_total_P(P, Tt, Y_AIR, m, A_DUCT).Pt, rel=1e-15)
    else:
        rho = P / (cb.specific_gas_constant(X_AIR) * Tt)
        assert tot.q_dyn == pytest.approx(0.5 * (m / A_DUCT) ** 2 / rho, rel=1e-12)
    # Even in the flow.
    rev = cb.hole_supply_drive(cb.HoleSupplyDrive.Total, P, Tt, Y_AIR, -m, A_DUCT, incompressible)
    assert rev.P_drive == tot.P_drive and rev.dP_dm_dot == -tot.dP_dm_dot


@pytest.mark.parametrize("incompressible", [False, True])
def test_the_drives_partials_are_analytic(incompressible: bool) -> None:
    # Central differences, each step small against the head it moves.
    P, Tt, m = 1.2e5, 450.0, 0.35
    r = cb.hole_supply_drive(cb.HoleSupplyDrive.Total, P, Tt, Y_AIR, m, A_DUCT, incompressible)

    def f(P=P, Tt=Tt, m=m, Y=Y_AIR) -> float:
        return cb.hole_supply_drive(
            cb.HoleSupplyDrive.Total, P, Tt, Y, m, A_DUCT, incompressible
        ).P_drive

    hP, hT, hm = 1e-4 * r.q_dyn, 1e-3, 1e-6 * m
    assert (f(P=P + hP) - f(P=P - hP)) / (2 * hP) == pytest.approx(r.dP_dP, rel=1e-7)
    assert (f(Tt=Tt + hT) - f(Tt=Tt - hT)) / (2 * hT) == pytest.approx(r.dP_dTt, rel=1e-6)
    assert (f(m=m + hm) - f(m=m - hm)) / (2 * hm) == pytest.approx(r.dP_dm_dot, rel=1e-7)
    k = cb.species_index_from_name("O2")
    hY = 1e-6
    Yp, Ym = list(Y_AIR), list(Y_AIR)
    Yp[k] += hY
    Ym[k] -= hY
    fd = (f(Y=Yp) - f(Y=Ym)) / (2 * hY)
    assert fd == pytest.approx(r.dP_dY[k], rel=1e-4, abs=1e-6 * r.q_dyn)


# --- a duct-fed panel --------------------------------------------------------


def _duct_net(
    spec: HoleCd,
    angle: float = 30.0,
    pt_in: float = 1.30,
    p_out: float = 1.20,
    p_gas: float = 1.0,
    crossflow: dict | None = None,
):
    """One station of a backside duct (the EffusionLiner wiring: a PlenumNode
    between two crossflow segments), a single-hole panel into 1 atm.
    ``crossflow``: the panel's crossflow declaration (default: the liner's,
    both segments and the duct's area)."""
    if crossflow is None:
        crossflow = {"crossflow_segments": ("s0", "s1"), "crossflow_area": A_DUCT}
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=pt_in * P_ATM, Tt=290.0, Y=Y_AIR))
    g.add_node(PressureBoundary("out", Pt=p_out * P_ATM, Tt=290.0, Y=Y_AIR))
    g.add_node(PressureBoundary("gas", Pt=p_gas * P_ATM, Tt=290.0, Y=Y_AIR))
    g.add_node(PlenumNode("n1"))
    kappa = cb.STATION_KAPPA_BLEED_BASSETT
    g.add_element(
        CrossflowSegmentElement(
            "s0", "in", "n1", length=0.03, area=A_DUCT, entry_K=0.0, to_kappa=kappa, next_seg="s1"
        )
    )
    g.add_element(
        CrossflowSegmentElement(
            "s1", "n1", "out", length=0.03, area=A_DUCT, from_kappa=kappa, prev_seg="s0"
        )
    )
    t = 6 * D_HOLE * math.sin(math.radians(angle))
    g.add_element(
        EffusionPlateElement(
            "h", "n1", "gas", hole_diameter=D_HOLE, wall_thickness=t, pitch_x=0.06,
            pitch_y=0.06, panel_length=0.06, panel_width=0.06, angle_deg=angle,
            hole_cd=spec, **crossflow,
        )
    )  # fmt: skip
    return g


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r.get("__message__")
    return s, r


def _ideal_flux(pt: float, Tt: float, p_exit: float) -> float:
    """Isentropic mass flux [kg/(m^2 s)] from (pt, Tt) to p_exit, an
    independent closed form (calorically perfect at the supply's gamma)."""
    st = cb.complete_state(Tt, pt, X_AIR).thermo
    k, R = st.gamma, st.cp - st.cv
    r = max(p_exit, pt * (2.0 / (k + 1.0)) ** (k / (k - 1.0))) / pt
    return (
        pt
        * r ** ((k + 1.0) / (2.0 * k))
        * math.sqrt(2.0 * k / ((k - 1.0) * R * Tt) * (r ** (-(k - 1.0) / k) - 1.0))
    )


@pytest.mark.parametrize("ms", [True, False])
def test_the_total_drive_flows_the_cd_from_the_ducts_total(ms: bool) -> None:
    # The flow over the ideal flow from the duct's total at the station
    # (static + the closure's head at the mean of the two segment flows) is
    # the element's Cd; from the static it is not (that is the old wiring,
    # where the same ratio from the STATIC was the Cd).
    # Bases without the supply-crossflow term: with it, the drive converts
    # Eq. 17's static-referenced Cd exactly (the conversion tests), so the
    # flow over the total-drive ideal flow is not the correlation's Cd.
    spec = HoleCd(
        correlation="McGreehanSchotsch" if ms else "IdelchikThick", inlet_crossflow="none"
    )
    a_h = math.pi * D_HOLE**2 / 4.0
    for drive in ("recovery", "static"):
        # The full recovery (f = 1) is the Cd from the total; the v0.11.1
        # static drive is the Cd from the static.
        kw = {"supply_recovery": 1.0} if drive == "recovery" else {"supply_drive": "static"}
        _, r = _solve(_duct_net(dataclass_replace(spec, **kw)))
        d = r["__element_diag__"]["h"]
        P = r["n1.P"]
        m_d = 0.5 * (r["s0.m_dot"] + r["s1.m_dot"])
        pt_duct = _moving_total_P(P, 290.0, Y_AIR, m_d, A_DUCT).Pt
        from_total = r["h.m_dot"] / (a_h * _ideal_flux(pt_duct, 290.0, P_ATM))
        from_static = r["h.m_dot"] / (a_h * _ideal_flux(P, 290.0, P_ATM))
        assert d["Cd_supply_drive_model"] == drive
        if drive == "recovery":
            assert from_total == pytest.approx(d["Cd"], rel=2e-4)  # the closed form's gamma
            assert d["Cd_supply_drive_applied"] == 1.0
            assert d["dP_drive"] == pytest.approx(pt_duct - P_ATM, rel=1e-9)
        else:
            assert from_static == pytest.approx(d["Cd"], rel=2e-4)
            assert d["Cd_supply_drive_applied"] == 0.0
        # How much of the drive is the duct's head (information).
        assert d["Cd_supply_drive_head_fraction"] == pytest.approx(
            (pt_duct - P) / (pt_duct - P_ATM), rel=1e-6
        )


def dataclass_replace(spec: HoleCd, **kw) -> HoleCd:
    import dataclasses

    return dataclasses.replace(spec, **kw)


def test_a_hole_whose_station_static_is_below_the_exit_blows() -> None:
    # The supply-side band: station static < exit pressure < duct total. The
    # hole leans with the duct flow and blows (Gritsch 1998, Ma_c 0.6 at low
    # PR); from the static drive it reverses.
    spec = HoleCd(correlation="IdelchikThick", inlet_crossflow="none", supply_recovery=1.0)
    _, r = _solve(_duct_net(spec, pt_in=1.10, p_out=0.90))
    d = r["__element_diag__"]["h"]
    assert r["n1.P"] < P_ATM  # the station static is below the exit
    assert r["h.m_dot"] > 0.0 and d["is_ingesting"] == 0.0 and d["dP_drive"] > 0.0
    static = dataclass_replace(spec, supply_recovery=None, supply_drive="static")
    _, r = _solve(_duct_net(static, pt_in=1.10, p_out=0.90))
    assert r["n1.P"] < P_ATM and r["h.m_dot"] < 0.0
    assert r["__element_diag__"]["h"]["is_ingesting"] == 1.0


def test_static_on_a_station_reads_it_as_before() -> None:
    # supply_drive="static" on a PlenumNode station is the pre-0.12 wiring:
    # no drive pair, the original evaluation path, bit for bit.
    g = _duct_net(HoleCd(correlation="McGreehanSchotsch", supply_drive="static"))
    g.resolve_all_topology()
    s = NetworkSolver(g)
    s._propagate_states(np.array(s._build_x0(), dtype=float))
    el = g.elements["h"]
    assert el._hole_cd.supply_drive == "static"
    si = s._get_node_state(g.nodes["n1"], np.array(s._build_x0(), dtype=float))
    so = s._get_node_state(g.nodes["gas"], np.array(s._build_x0(), dtype=float))
    assert el._hole_cd.drive_pair(si, so, {"s0": 0.3, "s1": 0.2}) is None


@pytest.mark.parametrize("make", ["plate", "film", "orifice", "slot"])
def test_static_is_refused_where_the_hole_is_not_duct_fed(make: str) -> None:
    spec = HoleCd(supply_drive="static")
    build = {
        "plate": lambda: EffusionPlateElement(
            "e", "a", "b", hole_diameter=1e-3, wall_thickness=2e-3, pitch=6e-3, panel_area=0.01,
            hole_cd=spec,
        ),
        "film": lambda: FilmCoolingElement(
            "e", "a", "b", hole_diameter=1e-3, holes_per_row=10, lateral_pitch=3e-3,
            wall_thickness=3e-3, hole_cd=spec,
        ),
        "orifice": lambda: OrificeElement("e", "a", "b", diameter=0.01, hole_cd=spec),
        "slot": lambda: SlotFilmElement(
            "e", "a", "b", slot_height=1e-3, span=0.05, lip_thickness=1e-3,
            hole_cd=HoleCd(correlation="IdelchikThick", supply_drive="static"),
        ),
    }[make]  # fmt: skip
    with pytest.raises(ValueError, match="supply_drive='static'"):
        build()


def test_eq_17_on_the_total_drive_is_the_static_referenced_flow() -> None:
    # McGreehan-Schotsch's Eq. 17 Cd is static-referenced; on the duct-total
    # drive it is converted exactly, so the hole's flow is that Cd times the
    # ideal flow from the station's STATIC state (p1 and its static T) -- the
    # authors' own model -- and the reported Cd is that over the ideal flow
    # from the duct's total. Independent closed form, not the kernel.
    s, r = _solve(_duct_net(HoleCd(correlation="McGreehanSchotsch")))
    d = r["__element_diag__"]["h"]
    el = s.network.elements["h"]
    x = np.array(r["__x_solution__"])
    si = s._get_node_state(s.network.nodes["n1"], x)
    so = s._get_node_state(s.network.nodes["gas"], x)
    si.m_dot = so.m_dot = r["h.m_dot"]
    cd_s = el._effective_Cd(si, so, s._network_flows(el, x))
    m_d = 0.5 * (r["s0.m_dot"] + r["s1.m_dot"])
    T1 = _moving_static_T(si, m_d, A_DUCT).T
    a_h = math.pi * D_HOLE**2 / 4.0
    P1 = r["n1.P"]
    assert r["h.m_dot"] == pytest.approx(cd_s * a_h * _ideal_flux(P1, T1, P_ATM), rel=2e-4)
    pt = _moving_total_P(P1, 290.0, Y_AIR, m_d, A_DUCT).Pt
    assert d["Cd"] == pytest.approx(r["h.m_dot"] / (a_h * _ideal_flux(pt, 290.0, P_ATM)), rel=2e-4)
    assert d["Cd_supply_drive_model"] == "recovery" and d["Cd_supply_drive_applied"] == 1.0
    # It passes zero at the station static (f = 0), so the direction, and
    # dP_drive / is_ingesting with it, are taken from there.
    assert d["dP_drive"] == pytest.approx(P1 - P_ATM, rel=1e-12)
    # The pre-0.12 drive gives the same model at the node's TOTAL temperature:
    # the only difference is the static-temperature correction.
    _, r0 = _solve(_duct_net(HoleCd(correlation="McGreehanSchotsch", supply_drive="static")))
    corr = r["h.m_dot"] / r0["h.m_dot"] - 1.0
    assert 0.0 < corr < 0.03


def test_the_converted_eq_17_hole_passes_zero_at_the_station_static() -> None:
    # f = 0: the converted Cd's forward flow vanishes where the exit meets the
    # station static, and above it the hole does not blow (its reversed jet
    # discharges against that static; Eq. 17's Cd underflows there, as on the
    # static drive). Evaluated at fixed states, with drop-sized FD.
    s, r = _solve(_duct_net(HoleCd(correlation="McGreehanSchotsch")))
    x = np.array(r["__x_solution__"])
    p1 = r["n1.P"]
    i = s.unknown_names.index("h.m_dot")
    m = []
    for dpg in (-200.0, -5.0, -0.05, 0.05, 5.0, 200.0):
        s.network.nodes["gas"].Pt = p1 + dpg
        F = np.asarray(s._residuals(x))
        m.append(x[i] - F[i])
        assert _fd_check(s, x, rows=[i]) < 1e-5
    assert m[0] > m[1] > m[2] >= 0.0 and m[3] <= 0.0
    assert abs(m[2]) < 1e-2 * m[0]


def test_eq_17_on_the_total_drive_has_an_analytic_jacobian() -> None:
    s, r = _solve(_duct_net(HoleCd(correlation="McGreehanSchotsch")))
    x = np.array(r["__x_solution__"])
    assert _fd_check(s, x) < 1e-5


@pytest.mark.parametrize("area", [6e-5, 2e-5, 1.2e-5])
def test_a_hole_fed_from_a_flowing_chamber_passes_zero_continuously(area: float) -> None:
    # Review of e471aced: a hole fed FROM a MomentumChamberNode (not duct-fed)
    # into a flowing chamber. Its drive is the supply's Pt; the reversed jet
    # must discharge against that same Pt, not the supply's static, or the
    # flow jumps by the supply's dynamic head at the switch (it was -0.02 to
    # -0.16 of a full jet at the old switch). The supply chamber's area sets
    # its head: about 100 / 1000 / 3000 Pa here.
    g = FlowNetwork()
    g.add_node(PressureBoundary("cin", Pt=5.15e5, Tt=600.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("cool", area=area))
    g.add_element(LosslessConnectionElement("cfeed", "cin", "cool"))
    g.add_node(MassFlowBoundary("gas_in", m_dot=2.0, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("exit", Pt=5e5, Tt=1500.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("ch", area=0.01, main_inlet="feed"))
    g.add_element(LosslessConnectionElement("feed", "gas_in", "ch"))
    g.add_element(LosslessConnectionElement("tail", "ch", "exit"))
    g.add_element(_film())
    s, r = _solve(g)
    el = s.network.elements["eff"]
    x = np.array(r["__x_solution__"])
    si = s._get_node_state(s.network.nodes["cool"], x)
    so = s._get_node_state(s.network.nodes["ch"], x)
    si.m_dot = so.m_dot = r["eff.m_dot"]
    q = float(si.Pt) - float(si.P)
    assert q > 50.0
    flows = s._network_flows(el, x)

    def m_at(dp: float) -> tuple[float, float]:
        # The chamber static dp below the supply's total, at fixed states.
        st = copy.copy(so)
        st.P = float(si.Pt) - dp
        m = el._flow_and_partials(si, st, el._effective_Cd(si, st, flows), flows)[0]
        return m, el._hole_cd.drive_pressures(si, st, flows)[0] - st.P

    m_ref = m_at(2e4)[0]
    for dp in (-3.0, -0.3, 0.0, 0.3, 3.0):
        m, drive = m_at(dp)
        assert abs(m) < 0.02 * m_ref, (q, dp, m / m_ref)
        if abs(dp) >= 0.3:
            assert (m > 0.0) == (drive > 0.0), (dp, m, drive)  # is_ingesting agrees
    # Falsification: the supply's static as the reversed jet's back pressure
    # leaves the old jump.
    pair = el._hole_cd.drive_pair(si, so, flows)
    assert pair.p_rev == float(si.Pt)


def test_the_drive_is_smooth_through_a_duct_flow_reversal() -> None:
    # The station's mean duct flow enters only through the closure's |m| /
    # A, whose head and partials are even and smooth through zero: the hole
    # row's Jacobian holds at states where the mean flow crosses zero.
    s, r = _solve(_duct_net(HoleCd(correlation="McGreehanSchotsch")))
    x0 = np.array(r["__x_solution__"])
    i_h = s.unknown_names.index("h.m_dot")
    i_a, i_b = s.unknown_names.index("s0.m_dot"), s.unknown_names.index("s1.m_dot")
    for mean in (-1e-3, -1e-6, 0.0, 1e-6, 1e-3):
        x = x0.copy()
        x[i_a] = mean + 0.5 * x0[i_h]
        x[i_b] = mean - 0.5 * x0[i_h]
        assert _fd_check(s, x, rows=[i_h]) < 1e-5, mean
    for spec in (
        HoleCd(correlation="IdelchikThick", inlet_crossflow="none"),  # f = 0
        HoleCd(correlation="IdelchikThick", inlet_crossflow="none", supply_recovery=0.5),
        HoleCd(correlation="IdelchikThick", inlet_crossflow="none", supply_recovery=1.0),
    ):
        s2, r2 = _solve(_duct_net(spec))
        x = np.array(r2["__x_solution__"])
        x[i_a], x[i_b] = 0.5 * x[i_h], -0.5 * x[i_h]
        assert _fd_check(s2, x, rows=[i_h]) < 1e-5


def test_an_unknown_supply_drive_is_refused() -> None:
    with pytest.raises(ValueError, match="supply_drive 'bogus' is not one of"):
        HoleCd(supply_drive="bogus")


def test_the_group_carries_the_drive_keys_everywhere() -> None:
    # Fixed shape: a plenum-fed hole reports the drive as its node's total.
    _, r = _solve(_chamber_net(_film(), dp=0.03))
    d = r["__element_diag__"]["eff"]
    assert all(k in d for k in HOLE_CD_DIAGNOSTIC_KEYS)
    assert d["Cd_supply_drive_model"] == "recovery"
    assert d["Cd_supply_drive_applied"] == 0.0 and d["Cd_supply_drive_head_fraction"] == 0.0
    # f does not apply to a hole that is not duct-fed.
    assert d["Cd_supply_recovery_q1"] == 0.0 and d["Cd_supply_recovery_ratio"] == 1.0
    assert math.isnan(d["Cd_supply_recovery_fraction"])


# --- the exit side: a hole into a flowing chamber (#1087) --------------------


def _film(**kw) -> FilmCoolingElement:
    a: dict = {"hole_diameter": 1e-3, "holes_per_row": 10, "lateral_pitch": 3e-3}
    a.update({"wall_thickness": 3e-3, **kw})
    return FilmCoolingElement("eff", "cool", "ch", angle_deg=30.0, **a)


def _plate(**kw) -> EffusionPlateElement:
    a: dict = {"hole_diameter": 0.8e-3, "wall_thickness": 2e-3, "pitch": 6e-3}
    a.update({"panel_area": 0.01, **kw})
    return EffusionPlateElement("eff", "cool", "ch", angle_deg=30.0, **a)


def _slot(**kw) -> SlotFilmElement:
    a: dict = {"slot_height": 1e-3, "span": 0.05, "lip_thickness": 1e-3}
    a.update(kw)
    return SlotFilmElement("eff", "cool", "ch", **a)


def _chamber_net(element, dp: float, m_gas: float = 2.0, area: float = 0.01) -> FlowNetwork:
    """Hot gas through a merge chamber at ~5 bar (172 m/s at 2 kg/s: a
    dynamic head of ~3.4 % of P); the holes' coolant plenum at
    5e5 (1 + dp)."""
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("gas_in", m_dot=m_gas, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("exit", Pt=5e5, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cool", Pt=5e5 * (1.0 + dp), Tt=600.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("ch", area=area, main_inlet="feed"))
    g.add_element(LosslessConnectionElement("feed", "gas_in", "ch"))
    g.add_element(LosslessConnectionElement("tail", "ch", "exit"))
    g.add_element(element)
    return g


@pytest.mark.parametrize("make", [_film, _plate, _slot])
def test_a_hole_between_the_chambers_static_and_total_blows(make) -> None:
    # The #1087 band: P_static,chamber < Pt_coolant < Pt_chamber. The jet
    # discharges against the static, so it blows; is_ingesting agrees with
    # the sign of the solved flow.
    s, r = _solve(_chamber_net(make(), dp=-0.01))
    d = r["__element_diag__"]["eff"]
    pt_c = 5e5 * 0.99
    assert r["ch.P"] < pt_c < r["ch.Pt"]  # inside the band
    assert r["eff.m_dot"] > 0.0
    assert d["dP_drive"] == pytest.approx(pt_c - r["ch.P"], rel=1e-9)
    assert d.get("is_ingesting", 0.0) == 0.0
    # Falsification: with the Pt - Pt direction the same root reverses.
    el = s.network.elements["eff"]
    x = np.array([r[n] for n in s.unknown_names])
    si = s._get_node_state(s.network.nodes["cool"], x)
    so = s._get_node_state(s.network.nodes["ch"], x)
    si.m_dot = so.m_dot = r["eff.m_dot"]
    cd = el._effective_Cd(si, so)
    assert el._flow_and_partials(si, so, cd)[0] > 0.0
    el._hole_cd.exit_static = False
    try:
        assert el._flow_and_partials(si, so, cd)[0] < 0.0
    finally:
        el._hole_cd.exit_static = True


def test_a_hole_ingests_only_below_the_chambers_static() -> None:
    s, r = _solve(_chamber_net(_film(), dp=-0.05))
    d = r["__element_diag__"]["eff"]
    assert r["ch.P"] > 5e5 * 0.95
    assert r["eff.m_dot"] < 0.0 and d["is_ingesting"] == 1.0


def test_the_flow_passes_zero_where_the_coolant_meets_the_static() -> None:
    # At fixed node states, m_calc(Pt_coolant) is continuous through the
    # chamber's static: both forms vanish there. With the Pt - Pt direction
    # the forward form still gave a full jet at Pt_coolant = Pt_chamber and
    # the reverse form ~0, a jump inside the 3 Pa blend.
    s, r = _solve(_chamber_net(_film(), dp=-0.01))
    el = s.network.elements["eff"]
    x = np.array([r[n] for n in s.unknown_names])
    so = s._get_node_state(s.network.nodes["ch"], x)
    si = s._get_node_state(s.network.nodes["cool"], x)

    def m_at(pt: float) -> float:
        st = copy.copy(si)
        st.P = st.Pt = pt
        return el._flow_and_partials(st, so, 0.7)[0]

    P_m = float(so.P)
    m_ref = abs(m_at(P_m + 2e3))
    ms = [m_at(P_m + dp) for dp in (-3.0, -1.0, -0.3, 0.0, 0.3, 1.0, 3.0)]
    assert all(b > a for a, b in zip(ms, ms[1:], strict=False))  # monotone through zero
    assert max(abs(m) for m in ms) < 0.06 * m_ref  # ~sqrt(3 / 2000) of it
    assert abs(m_at(P_m)) < 1e-12 + 1e-3 * m_ref
    # The old switch: across Pt_chamber the flow jumps by most of a jet.
    el._hole_cd.exit_static = False
    try:
        Pt_m = float(so.Pt)
        jump = m_at(Pt_m + 3.0) - m_at(Pt_m - 3.0)
        assert jump > 0.3 * m_ref
    finally:
        el._hole_cd.exit_static = True


def test_a_still_chamber_is_unchanged() -> None:
    # A plenum (P = Pt) discharge: no drive pair, the original path.
    g = FlowNetwork()
    g.add_node(PressureBoundary("cool", Pt=1.3e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("ch", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(_film())
    g.resolve_all_topology()
    el = g.elements["eff"]
    assert el._hole_cd.exit_static is False
    s = NetworkSolver(g)
    x = np.array(s._build_x0(), dtype=float)
    si = s._get_node_state(g.nodes["cool"], x)
    so = s._get_node_state(g.nodes["ch"], x)
    assert el._hole_cd.drive_pair(si, so, None) is None


def test_second_law_holds_inside_the_band_and_on_a_liner() -> None:
    # The hole's control volume runs between its drive pair: from the
    # coolant's total to the chamber static it generates entropy; closed at
    # the chamber's Pt (the old way) a blowing hole in the band destroyed it.
    s, r = _solve(_chamber_net(_film(), dp=-0.01))
    assert r["__second_law__"]["eff"]["verdict"] == "pass"
    assert r["__second_law__"]["eff"]["sigma_norm"] > 0.0
    for chamber in (False, True):
        for spec in (
            HoleCd(),  # the default: Eq. 17 converted onto the duct total
            HoleCd(inlet_crossflow="none", supply_recovery=1.0),  # f = 1
            HoleCd(inlet_crossflow="none"),  # f = 0, the default
            HoleCd(supply_drive="static"),
        ):
            _, r = _solve(_liner_net(spec, chamber=chamber))
            assert r["__second_law_failures__"] == []
            verdicts = {eid: v["verdict"] for eid, v in r["__second_law__"].items()}
            assert all(v == "pass" for v in verdicts.values()), verdicts


# --- the Jacobian, FD with drop-sized pressure steps (#1107) ----------------


def _liner_net(spec: HoleCd | None = None, chamber: bool = False) -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(PressureBoundary("cin", Pt=1.5e5, Tt=600.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cout", Pt=1.1e5, Tt=600.0, Y=Y_AIR))
    liner = EffusionLiner(
        "ln", n_segments=3, length=0.152, width=0.152, duct_height=0.01,
        hole_diameter=3.27e-3, wall_thickness=6.3e-3, pitch_x=15.24e-3, pitch_y=15.24e-3,
        angle_deg=30.0, **({} if spec is None else {"hole_cd": spec}),
    )  # fmt: skip
    if chamber:
        g.add_node(MassFlowBoundary("gas_in", m_dot=0.4, Tt=1500.0, Y=Y_AIR))
        g.add_node(PressureBoundary("exit", Pt=1.0e5, Tt=1500.0, Y=Y_AIR))
        g.add_node(MomentumChamberNode("ch", area=0.02, main_inlet="feed"))
        g.add_element(LosslessConnectionElement("feed", "gas_in", "ch"))
        g.add_element(LosslessConnectionElement("tail", "ch", "exit"))
        liner.add_to(g, "cin", "cout", "ch")
    else:
        g.add_node(PressureBoundary("gas", Pt=1.0e5, Tt=1500.0, Y=Y_AIR))
        liner.add_to(g, "cin", "cout", "gas")
    return g


def _fd_check(s: NetworkSolver, x: np.ndarray, rows: list[int] | None = None) -> float:
    """Worst scaled |J - J_fd| (J_ij |x_j|) over the rows (all by default):
    the shared helper's central differences, pressure steps capped at a
    fraction of the smallest pressure gap, i.e. sized on the drops (#1107)."""
    J = s._residuals_and_jacobian(x, compute_jacobian=True)[1].toarray()
    fd = fd_jacobian(s, x)
    err = scaled_error(J, fd, x)
    return float((err if rows is None else err[rows]).max())


@pytest.mark.parametrize("chamber", [False, True])
@pytest.mark.parametrize(
    "variant", ["default", "term_off", "term_off_f0", "term_off_half", "static"]
)
def test_the_liner_jacobian_is_analytic(variant: str, chamber: bool) -> None:
    spec = {
        "default": HoleCd(),
        "term_off": HoleCd(inlet_crossflow="none", supply_recovery=1.0),
        "term_off_f0": HoleCd(inlet_crossflow="none"),
        "term_off_half": HoleCd(inlet_crossflow="none", supply_recovery=0.5),
        "static": HoleCd(supply_drive="static"),
    }[variant]
    s, r = _solve(_liner_net(spec, chamber=chamber))
    x = np.array(r["__x_solution__"])
    assert _fd_check(s, x) < 1e-5


@pytest.mark.parametrize("make", [_film, _plate])
def test_the_band_jacobian_is_analytic(make) -> None:
    s, r = _solve(_chamber_net(make(), dp=-0.01))
    x = np.array([r[n] for n in s.unknown_names])
    assert _fd_check(s, x) < 1e-5


def test_the_jacobian_is_analytic_inside_the_blend() -> None:
    # Coolant within 0.05 Pa of the chamber static: inside the 0.1 Pa
    # logistic blend of the two forms, at a state (not a root).
    s, r = _solve(_chamber_net(_film(), dp=-0.01))
    x = np.array([r[n] for n in s.unknown_names])
    # Move the coolant boundary onto the chamber's static.
    s.network.nodes["cool"].Pt = float(r["ch.P"]) + 0.05
    # Inside the 0.1 Pa logistic the shared helper's pressure-step floor
    # (1e-8 |P| = 5e-3 Pa here) is too coarse: a 5e-4 Pa central step.
    J = s._residuals_and_jacobian(x, compute_jacobian=True)[1].toarray()
    fd = np.zeros_like(J)
    for j, name in enumerate(s.unknown_names):
        h = 5e-4 if name.endswith((".P", ".Pt")) else 1e-6 * max(abs(x[j]), 1e-6)
        e = np.zeros_like(x)
        e[j] = h
        fd[:, j] = (np.asarray(s._residuals(x + e)) - np.asarray(s._residuals(x - e))) / (2 * h)
    assert float(scaled_error(J, fd, x).max()) < 1e-5


def test_a_duct_fed_hole_that_ingests_has_an_analytic_jacobian() -> None:
    # Gas above the duct's total: the panel ingests, and the reversed jet
    # discharges against the duct total (the drive), whose derivatives by
    # the station's P, T and the segment flows carry the row. The hole's own
    # row only: the segment rows' d/d(hole flow) of an INGESTING inclined
    # panel has the wrong sign before this change too (the joining jet's
    # momentum at the station, #813 path; reproduced on the static drive).
    spec = HoleCd(correlation="IdelchikThick", inlet_crossflow="none", supply_recovery=1.0)
    s, r = _solve(_duct_net(spec, pt_in=1.30, p_out=1.20, p_gas=1.32))
    d = r["__element_diag__"]["h"]
    assert r["h.m_dot"] < 0.0 and d["is_ingesting"] == 1.0
    x = np.array([r[n] for n in s.unknown_names])
    row = [s.unknown_names.index("h.m_dot")]
    assert _fd_check(s, x, rows=row) < 1e-5


# --- a duct-fed panel at a MomentumChamberNode station (review C1) -----------


def _chamber_station_net(spec: HoleCd, area: float = 3e-4, angle: float = 30.0) -> FlowNetwork:
    """A dead-end duct into a MomentumChamberNode station (it carries its
    own Pt), a single-hole duct-fed panel into 1 atm: crossflow_segments
    ("s0", None)."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=1.30 * P_ATM, Tt=290.0, Y=Y_AIR))
    g.add_node(PressureBoundary("gas", Pt=P_ATM, Tt=290.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("n1", area=area))
    g.add_element(CrossflowSegmentElement("s0", "in", "n1", length=0.03, area=A_DUCT, entry_K=0.0))
    t = 6 * D_HOLE * math.sin(math.radians(angle))
    g.add_element(
        EffusionPlateElement(
            "h", "n1", "gas", hole_diameter=D_HOLE, wall_thickness=t, pitch_x=0.06,
            pitch_y=0.06, panel_length=0.06, panel_width=0.06, angle_deg=angle,
            hole_cd=spec, crossflow_segments=("s0", None), crossflow_area=A_DUCT,
        )
    )  # fmt: skip
    return g


def _chamber_static_T(Tt: float, P: float, G: float) -> float:
    """Static T at stagnation Tt, static P and mass flux G, from h(Tt) = h(T)
    + (G R T / P)^2 / 2 by bisection: independent of the C++ closure."""
    R = cb.specific_gas_constant(X_AIR)
    ht = cb.h_mass(Tt, X_AIR)
    lo, hi = 0.3 * Tt, Tt
    for _ in range(200):
        T = 0.5 * (lo + hi)
        if ht - cb.h_mass(T, X_AIR) - 0.5 * (G * R * T / P) ** 2 > 0.0:
            lo = T
        else:
            hi = T
    return 0.5 * (lo + hi)


def test_eq_17_at_a_chamber_station_is_the_static_referenced_flow() -> None:
    # A MomentumChamberNode station carries its own Pt, the drive already.
    # Eq. 17's static-referenced Cd is converted onto it exactly, as at a
    # PlenumNode station: the flow is Cd_S times the ideal flow from the
    # chamber's STATIC P and static T (its own closure at its through-flow).
    # Unconverted (Cd_S on the chamber's Pt, the pre-fix wiring) it is 1.5 %
    # more here. Independent closed form, not the kernel.
    area = 3e-4
    s, r = _solve(_chamber_station_net(HoleCd(correlation="McGreehanSchotsch"), area=area))
    el = s.network.elements["h"]
    x = np.array(r["__x_solution__"])
    si = s._get_node_state(s.network.nodes["n1"], x)
    so = s._get_node_state(s.network.nodes["gas"], x)
    si.m_dot = so.m_dot = r["h.m_dot"]
    cd_s = el._effective_Cd(si, so, s._network_flows(el, x))
    P1, Pt1, m = r["n1.P"], r["n1.Pt"], r["h.m_dot"]
    assert Pt1 - P1 > 500.0  # a chamber with a head worth converting
    T1 = _chamber_static_T(290.0, P1, abs(m) / area)
    a_h = math.pi * D_HOLE**2 / 4.0
    assert m == pytest.approx(cd_s * a_h * _ideal_flux(P1, T1, P_ATM), rel=2e-4)
    unconverted = cd_s * a_h * _ideal_flux(Pt1, 290.0, P_ATM)
    assert 0.005 < unconverted / m - 1.0 < 0.03
    d = r["__element_diag__"]["h"]
    assert d["Cd"] == pytest.approx(m / (a_h * _ideal_flux(Pt1, 290.0, P_ATM)), rel=2e-4)
    # It passes zero at the chamber static (f = 0).
    assert d["dP_drive"] == pytest.approx(P1 - P_ATM, rel=1e-12)


def test_eq_17_at_a_chamber_station_has_an_analytic_jacobian() -> None:
    for area in (1.2e-3, 3e-4):
        s, r = _solve(_chamber_station_net(HoleCd(correlation="McGreehanSchotsch"), area=area))
        x = np.array(r["__x_solution__"])
        assert _fd_check(s, x) < 1e-5


# --- duct-flow reversal: a refused proxy root is named (review B1, #1125) -----


def test_a_proxy_root_refused_by_the_second_law_is_named() -> None:
    # A 3-segment liner on a 5 mm duct, 1.3e5 into 1.25e5 Pa with the gas
    # side a flowing chamber: the duct flow reverses past the first station
    # (+---), and booked at the duct total the station's reversed-flow model
    # manufactures total pressure (#1125). The solve fails; its message must
    # say the warm-start proxy root was refused by the second law, name the
    # element and its sigma_norm and point at #1125, not just "did not
    # converge".
    g = FlowNetwork()
    g.add_node(PressureBoundary("cin", Pt=1.3e5, Tt=600.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cout", Pt=1.25e5, Tt=600.0, Y=Y_AIR))
    g.add_node(PressureBoundary("gas_in", Pt=1.003e5, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("exit", Pt=1.0e5, Tt=1500.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("ch", main_inlet="duct"))
    g.add_element(ChannelElement("duct", "gas_in", "ch", length=0.2, diameter=0.08))
    g.add_element(ChannelElement("tail", "ch", "exit", length=0.2, diameter=0.08))
    EffusionLiner(
        "ln", n_segments=3, length=0.152, width=0.152, duct_height=0.005,
        hole_diameter=3.27e-3, wall_thickness=6.3e-3, pitch_x=15.24e-3,
        pitch_y=15.24e-3, angle_deg=30.0,
    ).add_to(g, "cin", "cout", "ch")  # fmt: skip
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(g).solve(timeout=60)
    assert not r["__success__"]
    msg = r["__message__"]
    assert "proxy root was refused by the second law at ln__s" in msg
    assert "(sigma_norm -" in msg
    assert "#1125" in msg


# --- the station's area from the chain (B1) ----------------------------------


def _duct_net_declared(**crossflow) -> FlowNetwork:
    """_duct_net with the panel's crossflow declaration given by the test."""
    return _duct_net(HoleCd(correlation="McGreehanSchotsch"), crossflow=crossflow)


def test_an_omitted_crossflow_area_is_the_stations_own() -> None:
    """A duct-fed panel on a chain station takes the station's section when
    crossflow_area is omitted, and flows what the declared area flows (the
    chain's area is re-derived from a diameter: equal to 1e-12, not bit for
    bit)."""
    _, given = _solve(_duct_net(HoleCd(correlation="McGreehanSchotsch")))
    g = _duct_net_declared(crossflow_segments=("s0", "s1"))
    g.resolve_all_topology()
    assert g.elements["h"]._crossflow_station.area == g.elements["s1"].area_from
    _, derived = _solve(g)
    for key in ("h.m_dot", "s0.m_dot", "n1.P", "n1.Pt", "n1.T"):
        assert derived[key] == pytest.approx(given[key], rel=1e-10)


@pytest.mark.parametrize(
    ("kw", "match"),
    [
        ({"crossflow_segments": ("s1", "s0"), "crossflow_area": A_DUCT}, "not the segments"),
        ({"crossflow_segments": ("s0", None), "crossflow_area": A_DUCT}, "not the segments"),
        ({"crossflow_segments": ("s0", "s1"), "crossflow_area": 2.0 * A_DUCT}, "not the area"),
    ],
)
def test_a_declaration_the_station_contradicts_is_refused(kw: dict, match: str) -> None:
    """One station has one moving state: declared segments or an area the
    chain's station contradicts are refused at validate."""
    g = _duct_net_declared(**kw)
    g.resolve_all_topology()
    with pytest.raises(ValueError, match=match):
        g.validate()


def test_an_omitted_crossflow_area_needs_a_chain_station() -> None:
    g = FlowNetwork()
    g.add_node(PlenumNode("b"))
    g.add_node(PressureBoundary("gas", Pt=1.0e5, Tt=1500.0, Y=Y_AIR))
    g.add_element(
        EffusionPlateElement(
            "p", "b", "gas", hole_diameter=1e-3, wall_thickness=2e-3, pitch=6e-3,
            panel_area=0.01, crossflow_segments=(None, None),
        )
    )  # fmt: skip
    g.resolve_all_topology()
    with pytest.raises(ValueError, match="needs the backside duct's crossflow_area"):
        g.validate()


# --- the recovered fraction f (HoleCd.supply_recovery, owner 2026-11-13) ----


@pytest.mark.parametrize("f", [None, 0.0, 0.4, 1.0])
@pytest.mark.parametrize("angle", [90.0, 30.0])
def test_any_base_is_referenced_to_the_partly_recovered_drive(
    f: float | None, angle: float
) -> None:
    # A base without the supply-crossflow term on the duct-total drive: the
    # flow is the base Cd times the ideal flow from P1 + f q1 at the
    # isentropic T there; f = None is 0 at every angle. Independent closed
    # form (_ideal_flux), not the kernel.
    spec = HoleCd(correlation="IdelchikThick", inlet_crossflow="none", supply_recovery=f)
    s, r = _solve(_duct_net(spec, angle=angle))
    d = r["__element_diag__"]["h"]
    el = s.network.elements["h"]
    x = np.array(r["__x_solution__"])
    si = s._get_node_state(s.network.nodes["n1"], x)
    so = s._get_node_state(s.network.nodes["gas"], x)
    si.m_dot = so.m_dot = r["h.m_dot"]
    cd_b = el._effective_Cd(si, so, s._network_flows(el, x))
    m_d = 0.5 * (r["s0.m_dot"] + r["s1.m_dot"])
    T1 = _moving_static_T(si, m_d, A_DUCT).T
    P1 = r["n1.P"]
    pt = _moving_total_P(P1, 290.0, Y_AIR, m_d, A_DUCT).Pt
    fv = 0.0 if f is None else f
    pf = P1 + fv * (pt - P1)
    Tf = T1 * (290.0 / T1) ** (math.log(pf / P1) / math.log(pt / P1))
    a_h = math.pi * D_HOLE**2 / 4.0
    assert r["h.m_dot"] == pytest.approx(cd_b * a_h * _ideal_flux(pf, Tf, P_ATM), rel=2e-4)
    assert d["Cd_supply_recovery_fraction"] == fv
    assert d["Cd_supply_recovery_q1"] == pytest.approx(pt - P1, rel=1e-9)
    assert d["Cd_supply_recovery_ratio"] == pytest.approx(
        _ideal_flux(pf, Tf, P_ATM) / _ideal_flux(pt, 290.0, P_ATM), rel=2e-4
    )
    # It passes zero at P1 + f q1: dP_drive is taken from there.
    assert d["dP_drive"] == pytest.approx(pf - P_ATM, rel=1e-9)
    assert _fd_check(s, x) < 1e-5


def test_f_zero_is_the_static_referenced_conversion_bit_for_bit() -> None:
    # A plenum-fed base with f = 0 given and with the default: one model, the
    # static-referenced conversion (hole_static_to_total_cd_ratio is the
    # kernel's f = 0 case).
    off = HoleCd(correlation="IdelchikThick", inlet_crossflow="none")
    _, a = _solve(_duct_net(off))
    _, b = _solve(_duct_net(dataclass_replace(off, supply_recovery=0.0)))
    assert a["h.m_dot"] == b["h.m_dot"] and a["n1.P"] == b["n1.P"]
    Y = Y_AIR
    for inc in (False, True):
        k0 = cb.hole_supply_recovery_cd_ratio(1.2e5, 430.0, 1.32e5, 450.0, 1.05e5, Y, 0.0, inc)
        ks = cb.hole_static_to_total_cd_ratio(1.2e5, 430.0, 1.32e5, 450.0, 1.05e5, Y, inc)
        assert (k0.ratio, k0.d_dP1, k0.d_dT1, k0.d_dPt1, k0.d_dTt1, k0.d_dP2) == (
            ks.ratio, ks.d_dP1, ks.d_dT1, ks.d_dPt1, ks.d_dTt1, ks.d_dP2,
        )  # fmt: skip


def test_f_one_takes_the_total_drive_unconverted() -> None:
    # supply_recovery=1.0: no conversion at all, the Cd from the duct total
    # (the drive pair carries no conv and passes zero at the duct total).
    s, r = _solve(_duct_net(HoleCd(correlation="IdelchikThick", inlet_crossflow="none",
                                   supply_recovery=1.0)))  # fmt: skip
    el = s.network.elements["h"]
    x = np.array(r["__x_solution__"])
    si = s._get_node_state(s.network.nodes["n1"], x)
    so = s._get_node_state(s.network.nodes["gas"], x)
    pair = el._hole_cd.drive_pair(si, so, s._network_flows(el, x))
    assert pair.conv is None and pair.p_rev == pair.p_up
    assert r["__element_diag__"]["h"]["Cd_supply_recovery_ratio"] == 1.0


@pytest.mark.parametrize("make", ["plate", "film", "orifice", "slot", "impingement"])
def test_supply_recovery_is_refused_where_the_hole_is_not_duct_fed(make: str) -> None:
    spec = HoleCd(supply_recovery=0.0)
    build = {
        "plate": lambda: EffusionPlateElement(
            "e", "a", "b", hole_diameter=1e-3, wall_thickness=2e-3, pitch=6e-3, panel_area=0.01,
            hole_cd=spec,
        ),
        "film": lambda: FilmCoolingElement(
            "e", "a", "b", hole_diameter=1e-3, holes_per_row=10, lateral_pitch=3e-3,
            wall_thickness=3e-3, hole_cd=spec,
        ),
        "orifice": lambda: OrificeElement("e", "a", "b", diameter=0.01, hole_cd=spec),
        "slot": lambda: SlotFilmElement(
            "e", "a", "b", slot_height=1e-3, span=0.05, lip_thickness=1e-3,
            hole_cd=HoleCd(correlation="IdelchikThick", supply_recovery=1.0),
        ),
        "impingement": lambda: ImpingementPlateElement(
            "e", "a", "b", d_jet=2.54e-3, xn_d=5.0, yn_d=4.0, z_d=1.0, span=0.12,
            plate_thickness=2.54e-3,
            hole_cd=HoleCd(correlation="fixed", Cd=0.8, supply_recovery=0.0),
        ),
    }[make]  # fmt: skip
    with pytest.raises(ValueError, match="applies only to holes fed by a duct crossflow"):
        build()


@pytest.mark.parametrize("bad", [-0.1, 1.5, math.nan, math.inf])
def test_a_supply_recovery_outside_zero_one_is_refused(bad: float) -> None:
    with pytest.raises(ValueError, match="supply_recovery"):
        HoleCd(supply_recovery=bad)


def test_supply_recovery_with_the_static_drive_is_refused() -> None:
    with pytest.raises(ValueError, match="supply_drive='static'"):
        HoleCd(supply_drive="static", supply_recovery=0.0)


def test_the_removed_total_drive_is_refused() -> None:
    # supply_drive="total" never shipped: supply_recovery=1.0 replaces it.
    with pytest.raises(ValueError, match="supply_drive 'total' is not one of"):
        HoleCd(supply_drive="total")


@pytest.mark.parametrize("f", [0.0, 0.5, 1.0])
def test_supply_recovery_is_refused_where_the_correlation_owns_the_crossflow(f: float) -> None:
    # Review of C1: McGreehan-Schotsch's U1/Vi term models the supply
    # crossflow on the station static (f = 0 by definition); a recovery on
    # top rebuilds the unconverted duct-total model and counts the crossflow
    # twice. Any value is refused, 0.0 included, so the rule is one sentence.
    with pytest.raises(ValueError, match="already models the supply crossflow"):
        _duct_net(HoleCd(correlation="McGreehanSchotsch", supply_recovery=f))
    with pytest.raises(ValueError, match="already models the supply crossflow"):
        _liner_net(HoleCd(supply_recovery=f))
    # With the term off the same base is plenum-fed, and f applies.
    _duct_net(HoleCd(correlation="McGreehanSchotsch", inlet_crossflow="none", supply_recovery=f))


@pytest.mark.parametrize("f", [None, 0.5])
def test_the_incompressible_liner_jacobian_is_analytic(f: float | None) -> None:
    # Review of C1: the incompressible kernel's sqrt(1 - (1 - f) X) through a
    # term-off duct-fed liner, every element in the incompressible regime.
    g = _liner_net(HoleCd(inlet_crossflow="none", supply_recovery=f))
    for el in g.elements.values():
        if hasattr(el, "regime"):
            el.regime = "incompressible"
    s, r = _solve(g)
    h = [e for e in g.elements.values() if isinstance(e, EffusionPlateElement)]
    assert h and all(e.regime == "incompressible" for e in h)
    d = r["__element_diag__"][h[1].id]
    assert d["Cd_supply_recovery_fraction"] == (0.0 if f is None else f)
    assert d["Cd_supply_recovery_ratio"] < 1.0
    x = np.array(r["__x_solution__"])
    assert _fd_check(s, x) < 1e-5
