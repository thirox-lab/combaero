"""The opt-in exit-crossflow Cd term (Rowbury, Oldfield & Lock 2001, #1074).

Default off: no diagnostics keys, and the stress networks' residuals and
Jacobians hash identically to main (the PR's bit-identity check). On, the
film and effusion elements put Eqs. 7-8 on their no-crossflow Cd with the
discharge chamber as the mainstream; the row's Jacobian is analytic and is
checked against central differences here. The fidelity and accuracy
evidence is validation/cooling/exit_crossflow_cd_runner.py, pinned below.
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
)

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.cooling import exit_crossflow_cd_runner as xr  # noqa: E402

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
SUCTION = cb.Rowbury2001Set.CylindricalSuctionSurface


# The exit-crossflow term opted in, every other choice the element's default.
ROWBURY = HoleCd(exit_crossflow="rowbury_2001")


def _plate(**kw) -> EffusionPlateElement:
    args = {
        "hole_diameter": 0.8e-3,
        "wall_thickness": 2e-3,
        "pitch": 6e-3,
        "panel_length": 0.1,
        "panel_width": 0.1,
        "angle_deg": 30.0,
    }
    args.update(kw)
    return EffusionPlateElement("eff", "cool", "ch", **args)


def _film(**kw) -> FilmCoolingElement:
    args = {
        "hole_diameter": 1e-3,
        "holes_per_row": 10,
        "lateral_pitch": 3e-3,
        "wall_thickness": 3e-3,
        "angle_deg": 30.0,
    }
    args.update(kw)
    return FilmCoolingElement("eff", "cool", "ch", **args)


def _net(
    element,
    m_gas: float = 2.0,
    dp: float = 0.03,
    area: float = 0.01,
    mixed_supply: bool = False,
    mixed_main: bool = False,
) -> FlowNetwork:
    """Hot gas through a merge chamber the holes discharge into, as the
    Behrendt 2008 network runner builds it. ``mixed_supply``: the coolant
    plenum is fed by two streams of different T and composition, so its Tt
    and Y move with the feed orifices' flows (unknowns)."""
    g = FlowNetwork()
    if mixed_main:
        # The main stream mixes two gases first: its Tt and Y move with
        # their flows, through the plenum the feed draws from.
        i_co2 = cb.species_index_from_name("CO2")
        y_co2 = [0.6 * y + (0.4 if i == i_co2 else 0.0) for i, y in enumerate(Y_AIR)]
        g.add_node(PressureBoundary("ga", Pt=5.3e5, Tt=1300.0, Y=Y_AIR))
        g.add_node(PressureBoundary("gb", Pt=5.3e5, Tt=1700.0, Y=y_co2))
        g.add_node(PlenumNode("gas_in"))
        for eid, src, d in (("ga_o", "ga", 0.06), ("gb_o", "gb", 0.05)):
            g.add_element(
                OrificeElement(eid, src, "gas_in", Cd=0.7, diameter=d, correlation="fixed")
            )
    else:
        g.add_node(MassFlowBoundary("gas_in", m_dot=m_gas, Tt=1500.0, Y=Y_AIR))
    ex = PressureBoundary("exit")
    ex.Pt, ex.Tt, ex.Y = 5e5, 1500.0, Y_AIR
    g.add_node(ex)
    if mixed_supply:
        i_n2 = cb.species_index_from_name("N2")
        y_n2 = [0.5 * y + (0.5 if i == i_n2 else 0.0) for i, y in enumerate(Y_AIR)]
        g.add_node(PressureBoundary("sa", Pt=5.4e5, Tt=450.0, Y=Y_AIR))
        g.add_node(PressureBoundary("sb", Pt=5.4e5, Tt=750.0, Y=y_n2))
        g.add_node(PlenumNode("cool"))
        for eid, src, d in (("oa", "sa", 4e-3), ("ob", "sb", 3e-3)):
            g.add_element(OrificeElement(eid, src, "cool", Cd=0.7, diameter=d, correlation="fixed"))
    else:
        g.add_node(PressureBoundary("cool", Pt=5e5 * (1.0 + dp), Tt=600.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("ch", area=area, main_inlet="feed"))
    g.add_element(LosslessConnectionElement("feed", "gas_in", "ch"))
    g.add_element(LosslessConnectionElement("tail", "ch", "exit"))
    g.add_element(element)
    return g


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r.get("__message__")
    return s, r


# --- defaults and refusals ------------------------------------------------


@pytest.mark.parametrize("make", [_plate, _film])
def test_off_by_default_reports_the_term_as_off(make) -> None:
    # The hole-Cd group keeps its shape: the term's keys are there, neutral.
    e = make()
    assert e.cd_exit_crossflow == "none" and e._hole_cd.exit_xf is None
    _, r = _solve(_net(e))
    d = r["__element_diag__"]["eff"]
    assert d["Cd_exit_crossflow_model"] == "none" and d["Cd_exit_crossflow_applied"] == 0.0
    assert d["Cd_exit_crossflow_ratio"] == 1.0 and d["Cd_base"] == d["Cd"]
    assert d["Cd_exit_crossflow_rho_u2"] == 0.0 and math.isnan(d["Cd_exit_crossflow_k_perp"])


@pytest.mark.parametrize("make", [_plate, _film])
def test_opting_in_without_a_flowing_gas_side_is_refused(make) -> None:
    g = FlowNetwork()
    g.add_node(PressureBoundary("cool", Pt=1.3e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("ch", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(make(hole_cd=ROWBURY))
    with pytest.raises(ValueError, match="MomentumChamberNode"):
        NetworkSolver(g).solve()


def test_a_hole_shallower_than_the_sets_zero_crossing_is_refused() -> None:
    # Suction cylindrical A = 0.682 sin(theta) - 0.24 < 0 below 20.6 deg.
    with pytest.raises(ValueError, match="A > 0 and B > 0"):
        _plate(angle_deg=20.0, hole_cd=ROWBURY)
    # The pressure-surface set reaches it (Table 1 rows 1-2 are 20-21 deg).
    _plate(angle_deg=20.0, hole_cd=ROWBURY, exit_crossflow_region="pressure_surface")


def test_fan_shaped_has_no_leading_edge_set() -> None:
    shaped = {"hole_type": "shaped", "area_ratio": 3.0, "breakout_width": 2.5e-3}
    with pytest.raises(ValueError, match="no fan-shaped leading-edge"):
        _film(hole_cd=ROWBURY, exit_crossflow_region="leading_edge", **shaped)
    e = _film(hole_cd=ROWBURY, **shaped)
    assert e._hole_cd.exit_xf.coefficient_set == cb.Rowbury2001Set.FanShapedSuctionSurface
    assert e._hole_cd.exit_xf.exit_area_ratio == 3.0


@pytest.mark.parametrize(
    "kw",
    [
        {"exit_crossflow_region": "hub"},
        {"exit_crossflow_orientation_deg": 95.0},
        {"exit_crossflow_orientation_deg": 90.0},
    ],
)
def test_bad_options_are_refused(kw: dict) -> None:
    with pytest.raises(ValueError):
        _plate(hole_cd=ROWBURY, **kw)


def test_an_unknown_exit_crossflow_model_is_refused() -> None:
    with pytest.raises(ValueError, match="exit_crossflow 'bogus' is not one of"):
        HoleCd(exit_crossflow="bogus")


# --- the element ------------------------------------------------------------


def _fd_rel_error(s: NetworkSolver, r: dict) -> float:
    """Worst |J - J_fd| |x_j| over the row's largest |J_fd| |x_j|, across
    every row and column (central differences)."""
    x0 = np.array([r.get(n, 0.0) for n in s.unknown_names])
    _, J = s._residuals_and_jacobian(x0, compute_jacobian=True)
    J = J.toarray()
    Jfd = fd_jacobian(s, x0)
    w = np.maximum(np.abs(x0), 1e-3)
    row = np.max(np.abs(Jfd) * w, axis=1) + 1e-30
    return float(np.max(np.abs(J - Jfd) * w / row[:, None]))


@pytest.mark.parametrize("make", [_plate, _film])
@pytest.mark.parametrize(
    ("spec", "kw"),
    [
        ({}, {}),
        ({"inclination": "idelchik_3_2"}, {}),
        ({"correlation": "fixed", "Cd": 0.7}, {}),
        ({}, {"angle_deg": 90.0}),
    ],
)
def test_the_jacobian_matches_central_differences(make, spec: dict, kw: dict) -> None:
    hole_cd = HoleCd(exit_crossflow="rowbury_2001", **spec)
    s, r = _solve(_net(make(hole_cd=hole_cd, **kw)))
    assert r["__element_diag__"]["eff"]["Cd_exit_crossflow_k_perp"] > 0.0
    assert _fd_rel_error(s, r) < 1e-6


def test_a_chamber_fed_only_by_holes_is_refused() -> None:
    # The holes' own discharge is not the crossflow they blow into.
    g = FlowNetwork()
    g.add_node(PressureBoundary("cool", Pt=1.3e5, Tt=300.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("ch", area=0.01))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(_film(hole_cd=ROWBURY))
    g.add_element(LosslessConnectionElement("tail", "ch", "out"))
    with pytest.raises(ValueError, match="approaching mainstream"):
        NetworkSolver(g).solve()


def test_the_mainstream_is_the_approaching_flow_not_the_holes() -> None:
    # The chamber's main inlet is the mainstream (a chamber with several
    # inflows must declare one): rho u^2 must not grow with the holes' own
    # discharge, which the chamber's total flow includes.
    def net(dp: float) -> FlowNetwork:
        g = FlowNetwork()
        g.add_node(MassFlowBoundary("gas_in", m_dot=2.0, Tt=1500.0, Y=Y_AIR))
        ex = PressureBoundary("exit")
        ex.Pt, ex.Tt, ex.Y = 5e5, 1500.0, Y_AIR
        g.add_node(ex)
        g.add_node(PressureBoundary("cool", Pt=5e5 * (1.0 + dp), Tt=600.0, Y=Y_AIR))
        g.add_node(MomentumChamberNode("ch", area=0.01, main_inlet="feed"))
        g.add_element(LosslessConnectionElement("feed", "gas_in", "ch"))
        g.add_element(LosslessConnectionElement("tail", "ch", "exit"))
        g.add_element(_plate(hole_cd=ROWBURY))
        return g

    s, a = _solve(net(0.03))
    assert s.network.elements["eff"]._hole_cd.exit_mainstream == ("feed",)
    _, b = _solve(net(0.08))
    qa = a["__element_diag__"]["eff"]["Cd_exit_crossflow_rho_u2"]
    qb = b["__element_diag__"]["eff"]["Cd_exit_crossflow_rho_u2"]
    assert b["eff.m_dot"] > 1.3 * a["eff.m_dot"]
    # The main stream's face density barely moves (the side streams'
    # impulse shifts the face pressure slightly): 0.06% for 1.3x the flow.
    assert qb == pytest.approx(qa, rel=2e-3)
    assert _fd_rel_error(s, a) < 1e-6


def test_the_jacobian_reads_a_mixed_coolant_supply() -> None:
    # rho_h, the isentropic exit density, reads the supply's Pt, Tt and composition:
    # two coolant streams of different T and Y mix in a plenum, so both
    # depend on the feed orifices' flows (unknowns).
    g = _net(_film(hole_cd=ROWBURY), mixed_supply=True)
    s, r = _solve(g)
    assert r["__element_diag__"]["eff"]["Cd_exit_crossflow_k_perp"] > 0.0
    assert _fd_rel_error(s, r) < 1e-6


def test_the_jacobian_reads_a_mixed_main_stream() -> None:
    # rho_m, the main stream's face density, reads its Tt and composition:
    # the main gas mixes two streams first, so both move with their flows.
    s, r = _solve(_net(_plate(hole_cd=ROWBURY), mixed_main=True))
    assert r["__element_diag__"]["eff"]["Cd_exit_crossflow_k_perp"] > 0.0
    assert _fd_rel_error(s, r) < 1e-6


def test_the_jacobian_holds_where_the_term_is_held() -> None:
    # A normal hole on the suction set (B 1.18 > 1, k_floor 0.23) weakly
    # INGESTING from a fast chamber: k below the floor, the Cd held (slope
    # zero), the row exact. (Forward, a hole into a moving chamber needs
    # Pt_c > Pt_m (#1087), so k >~ Cd^2 here; a set or angle whose floor is
    # higher -- the pressure set at theta 90, alpha 60: 2.9 -- is held while
    # blowing too.)
    s, r = _solve(_net(_plate(angle_deg=90.0, hole_cd=ROWBURY), dp=-0.029))
    d = r["__element_diag__"]["eff"]
    assert r["eff.m_dot"] < 0.0 and d["Cd_exit_crossflow_held"] == 1.0
    assert _fd_rel_error(s, r) < 1e-6


@pytest.mark.parametrize("make", [_plate, _film])
def test_the_elements_cd_is_eqs_7_8_on_its_own_cd0(make) -> None:
    _, r = _solve(_net(make(hole_cd=ROWBURY)))
    d = r["__element_diag__"]["eff"]
    k = d["Cd_exit_crossflow_k_perp"]
    assert d["Cd_exit_crossflow_k_used"] == k and d["Cd_exit_crossflow_held"] == 0.0
    expect = cb.rowbury_2001_exit_crossflow_cd(d["Cd_base"], k, SUCTION, 30.0, 0.0)
    assert d["Cd"] == pytest.approx(expect.Cd, rel=1e-12)
    # k from the solved state, by hand: the hole's actual flux, the coolant's
    # isentropic exit density, the approaching gas's G^2 / rho_chamber.
    e_area = np.pi * (1e-3 if make is _film else 0.8e-3) ** 2 / 4.0
    n = 10 if make is _film else d["n_holes"]
    G_h = abs(r["eff.m_dot"]) / (n * e_area)
    X = list(cb.species.dry_air())
    rho_h = cb.rowbury_2001_exit_density(5e5 * 1.03, 600.0, r["ch.P"], X).rho
    T_c = cb.rowbury_2001_exit_density(5e5 * 1.03, 600.0, r["ch.P"], X).T
    assert cb.s_mass(T_c, X, r["ch.P"]) == pytest.approx(cb.s_mass(600.0, X, 5e5 * 1.03))
    k_hand = G_h**2 * 0.5 / (rho_h * d["Cd_exit_crossflow_rho_u2"])
    assert k == pytest.approx(k_hand, rel=1e-9)


def test_crossflow_lowers_cd_at_low_k_and_raises_it_at_high_k() -> None:
    # Eq. 8's crossover: delta > 0 while A k^-B > C, < 0 beyond, inside the
    # set's fitted k (the leading-edge rows reach k_perp 4.3).
    lo = _solve(_net(_film(hole_cd=ROWBURY), m_gas=4.0, dp=0.01))[1]
    le = _film(
        angle_deg=88.0,
        hole_cd=ROWBURY,
        exit_crossflow_region="leading_edge",
        exit_crossflow_orientation_deg=60.0,
    )
    hi = _solve(_net(le, m_gas=1.5, dp=0.10))[1]
    dlo, dhi = lo["__element_diag__"]["eff"], hi["__element_diag__"]["eff"]
    assert dlo["Cd_exit_crossflow_delta"] > 0.0 and dlo["Cd_exit_crossflow_ratio"] < 1.0
    assert dhi["Cd_exit_crossflow_faded"] == 0.0
    assert dhi["Cd_exit_crossflow_delta"] < 0.0 and dhi["Cd_exit_crossflow_ratio"] > 1.0


@pytest.mark.parametrize("make", [_plate, _film])
def test_a_vanishing_mainstream_returns_cd0_continuously(make) -> None:
    # Beyond the fitted k the -C offset fades: no jump to Cd0 at u_m = 0,
    # and the faded region's row is exact.
    ratios = []
    for m_gas in (0.3, 0.03, 0.003):
        s, r = _solve(_net(make(hole_cd=ROWBURY), m_gas=m_gas))
        d = r["__element_diag__"]["eff"]
        assert d["Cd_exit_crossflow_faded"] == 1.0
        ratios.append(abs(d["Cd_exit_crossflow_ratio"] - 1.0))
        assert _fd_rel_error(s, r) < 1e-6
    assert ratios[0] > ratios[1] > ratios[2]
    assert ratios[2] < 0.005


def test_leading_edge_compound_angle_row_is_exact() -> None:
    # alpha 60 (the leading-edge rows): k on the surface-normal component,
    # u_perp = u / sqrt(cot^2 theta + 1 + tan^2 alpha); row checked by FD.
    e = _film(
        angle_deg=88.0,
        hole_cd=ROWBURY,
        exit_crossflow_region="leading_edge",
        exit_crossflow_orientation_deg=60.0,
    )
    s, r = _solve(_net(e))
    d = r["__element_diag__"]["eff"]
    assert d["Cd_exit_crossflow_in_range"] == 1.0
    assert _fd_rel_error(s, r) < 1e-6


def test_shaped_holes_take_the_exit_plane_flux() -> None:
    shaped = {"hole_type": "shaped", "area_ratio": 3.0, "breakout_width": 2.5e-3}
    _, r = _solve(_net(_film(hole_cd=ROWBURY, **shaped)))
    d = r["__element_diag__"]["eff"]
    G_m = abs(r["eff.m_dot"]) / (10 * np.pi * 1e-3**2 / 4.0)
    rho_h = cb.rowbury_2001_exit_density(5e5 * 1.03, 600.0, r["ch.P"], list(cb.species.dry_air()))
    k_exit = (G_m / 3.0) ** 2 * 0.5 / (rho_h.rho * d["Cd_exit_crossflow_rho_u2"])
    assert d["Cd_exit_crossflow_k_perp"] == pytest.approx(k_exit, rel=1e-9)


# --- the runner's numbers (cross-checked in the PR body) --------------------


def test_fidelity_reproduces_the_authors_fig9_predictions() -> None:
    # Ma_m 0.6 / 1.2 within the eye-reading (+-0.005): max 0.0027 both.
    # Ma_m 0.3's two lowest-PR markers (PR 1.016: -0.036, 1.028: -0.014) sit
    # on a near-vertical part of the curve where a 0.005 PR reading error is
    # 0.03 in Cd; the other six are within 0.008.
    rows = xr.fidelity_rows(as_printed=True)
    for mach in (0.6, 1.2):
        d = [cd - read for m, _, read, cd in rows if m == mach]
        assert max(abs(x) for x in d) < 0.003
    d03 = [cd - read for m, pr, read, cd in rows if m == 0.3 and pr > 1.05]
    assert max(abs(x) for x in d03) < 0.008
    # As shipped the suction set fades only beyond the authors' own Fig. 9
    # application (k_fit 3.8, tier 2): every marker is Eq. 8 as printed.
    shipped = xr.fidelity_rows(as_printed=False)
    for (_, _, _, c1), (_, _, _, c2) in zip(rows, shipped, strict=True):
        assert c2 == pytest.approx(c1, abs=1e-12)
    # The first version's rho(Tt_c, p_m) missed Ma_m 1.2 by +0.008 (the bug).
    old = xr.fidelity_rows(rho_h_supply_tt=True)
    d12 = [cd - read for m, _, read, cd in old if m == 1.2]
    assert sum(d12) / len(d12) > 0.006


def test_the_suction_sets_tier_2_is_the_authors_fig9_reach() -> None:
    # Fig. 9's Ma_m 0.3 prediction runs to PR 2.0; there k_perp is 3.8 on
    # the authors' own curve (Cd,no 0.776, isentropic exit density).
    _, k, _ = xr.solve_point(
        xr.GRITSCH_RIG, 2.0, 0.3, True, False, cd_no=xr.LICHTAROWICZ_CDU_L_OVER_D_6, as_printed=True
    )
    assert k == pytest.approx(cb.rowbury_2001_k_fit_max(SUCTION), abs=0.05)
    assert cb.rowbury_2001_k_data_max(SUCTION) == 0.58


def test_accuracy_before_and_after() -> None:
    got = dict(xr.accuracy())
    expect = {
        # (bias off, bias on), inclination term on: cross-source accuracy.
        "gritsch1998/fig3_Mac0p0_Mam0p3": (0.008, -0.028),
        "gritsch1998/fig3_Mac0p0_Mam0p6": (0.159, 0.032),
        "gritsch1998/fig3_Mac0p0_Mam1p2": (0.468, 0.086),
        # 90 deg: outside the suction rows (30-78 deg); depends on the
        # k_unique margin (M 0.25 / 0.5 / 0.75 / 1.0: -0.2 / -4.6 / -7.4 / -9.3%).
        "hay1983/fig6_open_circles": (0.112, -0.046),
        "hay1983/fig8_open_circles": (0.114, 0.068),
        "hay1983/fig12_open_circles": (0.123, 0.088),
    }
    for label, (off, on) in expect.items():
        assert xr.summary(got[label], True, False)[0] == pytest.approx(off, abs=0.003)
        assert xr.summary(got[label], True, True)[0] == pytest.approx(on, abs=0.003)


def test_no_mainstream_momentum_returns_the_no_crossflow_cd() -> None:
    # The authors' no-crossflow prediction is Cd,no by definition (Fig. 9,
    # "Predict Mm = 0"), and the faded Eq. 8 tends to it.
    s, r = _solve(_net(_film(hole_cd=ROWBURY)))
    e = s.network.elements["eff"]
    e._hole_cd.exit_node._inflows = [
        (eid, _Probe(0.0, 1500.0, Y_AIR), {}) for eid, _, _ in e._hole_cd.exit_node._inflows
    ]
    r, _ = e._hole_cd._evaluate(
        _Probe(1e-3, 600.0, Y_AIR), _Probe(0.0, 1500.0, Y_AIR), {}, jet=False
    )
    assert not r.exit_active
    assert r.cd_hole == r.cd_base
    assert r.dCd_dG_h == r.dCd_drho_h == r.dCd_drho_m_u_m2 == 0.0


class _Probe:
    """The few state fields _exit_crossflow reads."""

    def __init__(self, m_dot: float, Tt: float, Y: list[float]) -> None:
        self.m_dot, self.Tt, self.Y = m_dot, Tt, Y
        self.P = self.Pt = 5e5
        self.X = list(cb.mass_to_mole(Y))


def test_the_faded_flag_reports_the_k_tier() -> None:
    # 0 inside the suction rows' data (k_perp <= 0.58), 0.5 inside only the
    # authors' Fig. 9 application (<= 3.8, Eq. 8 unfaded), 1 beyond (faded).
    tiers = {}
    for m_gas, dp in ((4.0, 0.01), (1.5, 0.10), (0.03, 0.03)):
        _, r = _solve(_net(_film(hole_cd=ROWBURY), m_gas=m_gas, dp=dp))
        d = r["__element_diag__"]["eff"]
        tiers[d["Cd_exit_crossflow_faded"]] = d["Cd_exit_crossflow_k_used"]
    assert set(tiers) == {0.0, 0.5, 1.0}
    assert tiers[0.0] <= 0.58 < tiers[0.5] <= 3.8 < tiers[1.0]
