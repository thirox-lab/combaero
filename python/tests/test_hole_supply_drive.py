"""A hole's drive pair (#1087): the supply drive of a duct-fed hole
(``HoleCd.supply_drive``) and the direction of a cooling hole into a flowing
chamber.

Supply side. A duct-fed EffusionPlateElement's supply node is a duct station
holding the duct's STATIC pressure. ``supply_drive="total"`` takes
the hole's ideal flow from the duct's stagnation pressure there, static plus
the dynamic head of its mean flow (C++ hole_supply_drive, the momentum
chamber closure): the reference the plenum-fed Cd data are taken on, the
hole-Cd audit's option (a). "static" ("auto", the default and the wiring
before) drives from the station static, McGreehan & Schotsch's convention.

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

import combaero as cb
from combaero.network import (
    CrossflowSegmentElement,
    EffusionLiner,
    EffusionPlateElement,
    FilmCoolingElement,
    FlowNetwork,
    HoleCd,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    SlotFilmElement,
)
from combaero.network.components import _moving_total_P
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
    assert tot.q_dyn == sta.q_dyn == c.q_dyn > 0.0
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
    spec: HoleCd, angle: float = 30.0, pt_in: float = 1.30, p_out: float = 1.20, p_gas: float = 1.0
):
    """One station of a backside duct (the EffusionLiner wiring: a PlenumNode
    between two crossflow segments), a single-hole panel into 1 atm."""
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
            hole_cd=spec, crossflow_segments=("s0", "s1"), crossflow_area=A_DUCT,
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
    spec = (
        HoleCd(correlation="McGreehanSchotsch")
        if ms
        else HoleCd(correlation="IdelchikThick", inlet_crossflow="none")
    )
    a_h = math.pi * D_HOLE**2 / 4.0
    for drive in ("total", "static"):
        _, r = _solve(_duct_net(dataclass_replace(spec, supply_drive=drive)))
        d = r["__element_diag__"]["h"]
        P = r["n1.P"]
        m_d = 0.5 * (r["s0.m_dot"] + r["s1.m_dot"])
        pt_duct = _moving_total_P(P, 290.0, Y_AIR, m_d, A_DUCT).Pt
        from_total = r["h.m_dot"] / (a_h * _ideal_flux(pt_duct, 290.0, P_ATM))
        from_static = r["h.m_dot"] / (a_h * _ideal_flux(P, 290.0, P_ATM))
        assert d["Cd_supply_drive_model"] == drive
        if drive == "total":
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
    spec = HoleCd(correlation="IdelchikThick", inlet_crossflow="none", supply_drive="total")
    _, r = _solve(_duct_net(spec, pt_in=1.10, p_out=0.90))
    d = r["__element_diag__"]["h"]
    assert r["n1.P"] < P_ATM  # the station static is below the exit
    assert r["h.m_dot"] > 0.0 and d["is_ingesting"] == 0.0 and d["dP_drive"] > 0.0
    _, r = _solve(_duct_net(dataclass_replace(spec, supply_drive="static"), pt_in=1.10, p_out=0.90))
    assert r["n1.P"] < P_ATM and r["h.m_dot"] < 0.0
    assert r["__element_diag__"]["h"]["is_ingesting"] == 1.0


@pytest.mark.parametrize("drive", ["auto", "static"])
def test_static_on_a_station_reads_it_as_before(drive: str) -> None:
    # supply_drive="static" -- the default on a duct-fed hole, "auto" -- on a
    # PlenumNode station is the wiring as before: no drive pair, the original
    # evaluation path, bit for bit.
    g = _duct_net(HoleCd(correlation="McGreehanSchotsch", supply_drive=drive))
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


def test_an_unknown_supply_drive_is_refused() -> None:
    with pytest.raises(ValueError, match="supply_drive 'bogus' is not one of"):
        HoleCd(supply_drive="bogus")


def test_the_group_carries_the_drive_keys_everywhere() -> None:
    # Fixed shape: a plenum-fed hole reports the drive as its node's total.
    _, r = _solve(_chamber_net(_film(), dp=0.03))
    d = r["__element_diag__"]["eff"]
    assert all(k in d for k in HOLE_CD_DIAGNOSTIC_KEYS)
    assert d["Cd_supply_drive_model"] == "total"
    assert d["Cd_supply_drive_applied"] == 0.0 and d["Cd_supply_drive_head_fraction"] == 0.0


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
    for chamber, drive in ((False, "auto"), (True, "auto"), (False, "total"), (True, "total")):
        _, r = _solve(_liner_net(HoleCd(supply_drive=drive), chamber=chamber))
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


def _fd_check(
    s: NetworkSolver, x: np.ndarray, dp_min: float, rows: list[int] | None = None
) -> float:
    """Worst |J - J_fd| |x_j| over its row's largest |J_fd| |x_j|. Pressure
    columns step min(1e-6 |P|, 1e-3 dp_min) (#1107), others 1e-6 |x|."""
    _, J = s._residuals_and_jacobian(x, compute_jacobian=True)
    J = J.toarray()
    Jfd = np.zeros_like(J)
    for j, name in enumerate(s.unknown_names):
        if name.endswith((".P", ".Pt")):
            h = min(1e-6 * abs(x[j]), 1e-3 * dp_min)
        else:
            h = 1e-6 * max(abs(x[j]), 1e-6)
        e = np.zeros_like(x)
        e[j] = h
        Jfd[:, j] = (np.asarray(s._residuals(x + e)) - np.asarray(s._residuals(x - e))) / (2 * h)
    if rows is not None:
        J, Jfd = J[rows], Jfd[rows]
    scale = np.abs(Jfd * x[None, :]).max(axis=1, keepdims=True)
    scale[scale == 0.0] = 1.0
    return float((np.abs((J - Jfd) * x[None, :]) / scale).max())


@pytest.mark.parametrize("chamber", [False, True])
@pytest.mark.parametrize("drive", ["total", "static"])
def test_the_liner_jacobian_is_analytic(drive: str, chamber: bool) -> None:
    s, r = _solve(_liner_net(HoleCd(supply_drive=drive), chamber=chamber))
    x = np.array([r[n] for n in s.unknown_names])
    d = r["__element_diag__"]
    dp_min = min(abs(d[e]["dP_drive"]) for e in ("ln__p1", "ln__p2", "ln__p3"))
    assert _fd_check(s, x, dp_min) < 1e-5


@pytest.mark.parametrize("make", [_film, _plate])
def test_the_band_jacobian_is_analytic(make) -> None:
    s, r = _solve(_chamber_net(make(), dp=-0.01))
    x = np.array([r[n] for n in s.unknown_names])
    dp = abs(r["__element_diag__"]["eff"]["dP_drive"])
    assert _fd_check(s, x, dp) < 1e-5


def test_the_jacobian_is_analytic_inside_the_blend() -> None:
    # Coolant within 0.05 Pa of the chamber static: inside the 0.1 Pa
    # logistic blend of the two forms, at a state (not a root).
    s, r = _solve(_chamber_net(_film(), dp=-0.01))
    x = np.array([r[n] for n in s.unknown_names])
    # Move the coolant boundary onto the chamber's static.
    s.network.nodes["cool"].Pt = float(r["ch.P"]) + 0.05
    assert _fd_check(s, x, 0.5) < 1e-5


def test_a_duct_fed_hole_that_ingests_has_an_analytic_jacobian() -> None:
    # Gas above the duct's total: the panel ingests, and the reversed jet
    # discharges against the duct total (the drive), whose derivatives by
    # the station's P, T and the segment flows carry the row. The hole's own
    # row only: the segment rows' d/d(hole flow) of an INGESTING inclined
    # panel has the wrong sign before this change too (the joining jet's
    # momentum at the station, #813 path; reproduced on the static drive).
    spec = HoleCd(correlation="IdelchikThick", inlet_crossflow="none", supply_drive="total")
    s, r = _solve(_duct_net(spec, pt_in=1.30, p_out=1.20, p_gas=1.32))
    d = r["__element_diag__"]["h"]
    assert r["h.m_dot"] < 0.0 and d["is_ingesting"] == 1.0
    x = np.array([r[n] for n in s.unknown_names])
    row = [s.unknown_names.index("h.m_dot")]
    assert _fd_check(s, x, abs(d["dP_drive"]), rows=row) < 1e-5
