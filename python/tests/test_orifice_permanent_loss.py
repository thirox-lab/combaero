"""An orifice plate between pipes charges ISO 5167-2's permanent loss (#645).

Before, pipe -> orifice -> pipe booked the whole tapping differential as the
orifice's loss: about 2x too much at beta 0.7. The row now relates the two
nodes' STAGNATION difference dw to the tapping differential dp through
ISO 5167-2:2003 5.4.1,

    dw = k dp,  k = [s - C b^2] / [s + C b^2],  s = sqrt(1 - b^4 (1 - C^2)),

and C stays referenced to dp. Discharging into a node with no pipe beyond
it, the jet's head is lost there and the drop is still the tapping
differential. The tap type (#645, owner's ruling) is measurement: it only
picks the differential the diagnostics report.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

D = 0.1


def _k_by_hand(beta: float, C: float) -> float:
    s = math.sqrt(1.0 - beta**4 * (1.0 - C * C))
    return (s - C * beta * beta) / (s + C * beta * beta)


# --- the C++ relation --------------------------------------------------------


@pytest.mark.parametrize("beta", [0.1, 0.3, 0.5, 0.7, 0.75])
@pytest.mark.parametrize("C", [0.6, 0.62, 0.8])
def test_the_loss_ratio_is_isos_and_its_derivative_exact(beta: float, C: float) -> None:
    k, dk = cb.permanent_loss_ratio_and_dC(beta, C)
    assert k == pytest.approx(_k_by_hand(beta, C), rel=1e-14)
    h = 1e-6
    fd = (_k_by_hand(beta, C + h) - _k_by_hand(beta, C - h)) / (2 * h)
    assert dk == pytest.approx(fd, rel=1e-7)


@pytest.mark.parametrize("beta", [0.3, 0.5, 0.7])
def test_the_loss_ratio_agrees_with_isos_short_form(beta: float) -> None:
    """ISO 5167-2 5.4.1 also gives dw/dp ~ 1 - beta^1.9. A cross-check of the
    long form at C = 0.6, not a second source: they agree to 0.2% at beta
    0.3-0.5 and part by 4.8% at beta 0.7 (0.516 vs 0.492)."""
    k, _ = cb.permanent_loss_ratio_and_dC(beta, 0.6)
    assert k == pytest.approx(1.0 - beta**1.9, rel=0.05)


@pytest.mark.parametrize(("beta", "C"), [(0.0, 0.6), (1.0, 0.6), (0.5, 0.0), (0.5, -0.1)])
def test_the_loss_ratio_refuses_outside_its_domain(beta: float, C: float) -> None:
    with pytest.raises(ValueError, match="0 < beta < 1 and C > 0"):
        cb.permanent_loss_ratio_and_dC(beta, C)


@pytest.mark.parametrize("compressible", [False, True])
def test_the_tap_pressure_derivatives_are_exact(compressible: bool) -> None:
    from combaero import _solver_tools

    f = _solver_tools.orifice_tap_pressure_and_jacobian
    x0 = [1.15e5, 1.17e5, 1.05e5, 0.6, 0.61]
    t = f(*x0, compressible)
    assert not t.floored
    k = _k_by_hand(0.6, 0.61)
    assert t.p_tap == pytest.approx(x0[0] - (x0[1] - x0[2]) / k, rel=1e-14)
    names = ["d_dP_static_up", "d_dP_total_up", "d_dP_total_down", None, "d_dC"]
    for i, name in enumerate(names):
        if name is None:
            continue
        h = 1e-6 * abs(x0[i])
        up, dn = list(x0), list(x0)
        up[i] += h
        dn[i] -= h
        fd = (f(*up, compressible).p_tap - f(*dn, compressible).p_tap) / (2 * h)
        assert getattr(t, name) == pytest.approx(fd, rel=1e-6, abs=1e-12), name


def test_the_compressible_tap_pressure_is_floored_below_choke_only() -> None:
    """A recovered drop the incompressible relation turns into p_tap <= 0 is
    held at 0.1 P1 compressible (choked there, so the nozzle cannot see it);
    incompressible has no floor."""
    from combaero import _solver_tools

    args = (1.0e5, 1.0e5, 0.4e5, 0.7, 0.61)
    t_c = _solver_tools.orifice_tap_pressure_and_jacobian(*args, True)
    t_i = _solver_tools.orifice_tap_pressure_and_jacobian(*args, False)
    assert t_c.floored and t_c.p_tap == pytest.approx(0.1e5)
    assert t_c.d_dP_total_up == 0.0 and t_c.d_dC == 0.0
    assert not t_i.floored and t_i.p_tap < 0.0
    assert not t_c.in_range and not t_i.in_range


def test_the_range_flag_is_isos_pressure_ratio() -> None:
    from combaero import _solver_tools

    f = _solver_tools.orifice_tap_pressure_and_jacobian
    k = _k_by_hand(0.5, 0.6)
    # p_tap / p1 = 1 - dw / (k p1): 0.76 in range, 0.74 out.
    assert f(1e5, 1e5, 1e5 - 0.24e5 * k, 0.5, 0.6, False).in_range
    assert not f(1e5, 1e5, 1e5 - 0.26e5 * k, 0.5, 0.6, False).in_range


# --- the network row ---------------------------------------------------------


def _node(node_cls: type, nid: str, regime: str) -> object:
    """A node in the network's regime: a MomentumChamberNode closes Pt(P)
    compressibly by default, which an incompressible network would not."""
    return node_cls(nid, regime=regime) if node_cls is MomentumChamberNode else node_cls(nid)


def _between_pipes(
    node_cls: type,
    beta: float,
    regime: str,
    correlation: str,
    tap: str = "Flange",
    downstream_pipe: bool = True,
) -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=1.2e5, Tt=300.0))
    g.add_node(_node(node_cls, "a", regime))
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
    g.add_element(ChannelElement("p1", "hi", "a", length=0.5, diameter=D, regime=regime))
    to = "b" if downstream_pipe else "ho"
    if downstream_pipe:
        g.add_node(_node(node_cls, "b", regime))
    g.add_element(
        OrificeElement(
            "o",
            "a",
            to,
            Cd=0.6 if correlation == "fixed" else None,
            diameter=beta * D,
            regime=regime,
            correlation=correlation,
            tap=tap,
        )
    )
    if downstream_pipe:
        g.add_element(ChannelElement("p2", "b", "ho", length=0.5, diameter=D, regime=regime))
    return g


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r["__message__"]
    return s, r


def _iso_dp(r: dict, beta: float, C: float) -> float:
    """ISO's tapping differential for the solved flow: (m / (C E A))^2 / (2 rho1),
    with rho1 at the approach pipe's static pressure (#651)."""
    rho = cb.density(300.0, r["o.P_approach"], cb.species.dry_air())
    E = 1.0 / math.sqrt(1.0 - beta**4)
    A = math.pi * (beta * D) ** 2 / 4.0
    return (r["o.m_dot"] / (C * E * A)) ** 2 / (2.0 * rho)


@pytest.mark.parametrize("node_cls", [PlenumNode, MomentumChamberNode])
@pytest.mark.parametrize("beta", [0.3, 0.5, 0.7])
def test_between_pipes_the_orifice_charges_isos_permanent_loss(node_cls: type, beta: float) -> None:
    """Pt_a - Pt_b = k dp, with dp ISO's tapping differential of the solved
    flow. Before #645 it was dp itself: k = 0.90 / 0.73 / 0.52 here."""
    _, r = _solve(_between_pipes(node_cls, beta, "incompressible", "fixed"))
    dp = _iso_dp(r, beta, 0.6)
    k = _k_by_hand(beta, 0.6)
    # 1e-4: the network's air vs species.dry_air in rho; the change is 10-48%.
    assert r["a.Pt"] - r["b.Pt"] == pytest.approx(k * dp, rel=1e-4)
    assert r["o.pressure_recovery"] == 1.0
    assert r["o.dP_permanent"] == pytest.approx(r["a.Pt"] - r["b.Pt"], rel=1e-12)
    assert r["o.loss_ratio"] == pytest.approx(k, rel=1e-12)


@pytest.mark.parametrize("beta", [0.5, 0.7])
def test_into_a_node_with_no_pipe_the_tapping_differential_is_lost(beta: float) -> None:
    """pipe -> orifice -> boundary: the jet's head is lost in the sink, so the
    row is ISO's flow equation on p_a - p_sink, as before #645."""
    _, r = _solve(
        _between_pipes(PlenumNode, beta, "incompressible", "fixed", downstream_pipe=False)
    )
    assert r["o.pressure_recovery"] == 0.0
    assert _iso_dp(r, beta, 0.6) == pytest.approx(r["o.P_approach"] - 1.0e5, rel=1e-4)


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.parametrize("beta", [0.5, 0.7])
def test_the_node_type_does_not_move_the_flow(regime: str, beta: float) -> None:
    """The user draws pipe -> orifice -> pipe; the node type in between must
    not decide the flow (#651). The orifice takes p1 and rho1 from the
    approach pipe, not from the node: a PlenumNode reports P = Pt by
    definition, which read them a dynamic head too high, 1.5% flow at beta
    0.7 (pipe Mach 0.18). Its momentum-chamber twin has the pipe's area, so
    both nodes now reduce to the same pipe state."""
    _, rp = _solve(_between_pipes(PlenumNode, beta, regime, "fixed"))
    _, rc = _solve(_between_pipes(MomentumChamberNode, beta, regime, "fixed"))
    # The solves' own tolerance; the defect was 1.5e-2.
    assert rc["o.m_dot"] == pytest.approx(rp["o.m_dot"], rel=1e-6)
    assert rc["o.P_approach"] == pytest.approx(rp["o.P_approach"], rel=1e-6)
    # And the chamber's own static pressure is that pipe state.
    assert rc["a.P"] == pytest.approx(rc["o.P_approach"], rel=1e-9)


@pytest.mark.parametrize("correlation", ["ReaderHarrisGallagher"])
def test_the_tap_never_moves_the_flow_only_the_reported_differential(correlation: str) -> None:
    """Owner's ruling on #645: the tap is measurement. Every tap gives the
    same solved flow, to the last bit; dP_tap moves as 1/C_tap^2."""
    runs = {
        tap: _solve(_between_pipes(PlenumNode, 0.7, "incompressible", correlation, tap=tap))[1]
        for tap in ("Corner", "D_D2", "Flange")
    }
    m = {tap: r["o.m_dot"] for tap, r in runs.items()}
    assert m["Corner"] == m["D_D2"] == m["Flange"]
    dp = {tap: r["o.dP_tap"] for tap, r in runs.items()}
    assert len(set(dp.values())) == 3
    r = runs["Flange"]
    g = cb.OrificeGeometry()
    g.d, g.D = 0.07, D
    mu = cb.transport_state(300.0, r["a.P"], cb.species.dry_air()).mu
    Re_D = 4.0 * r["o.m_dot"] / (math.pi * D * mu)
    c = {}
    for tap in dp:
        g.tap = getattr(cb.OrificeTap, tap)
        c[tap] = cb.metering_cd_and_dRe(getattr(cb.MeteringCdCorrelation, correlation), g, Re_D)[0]
    for tap in dp:
        assert dp[tap] * c[tap] ** 2 == pytest.approx(dp["D_D2"] * c["D_D2"] ** 2, rel=1e-4)


@pytest.mark.parametrize("correlation", ["ReaderHarrisGallagher"])
def test_the_network_metering_cd_is_cpps_at_dd2_taps(correlation: str) -> None:
    """No metering formula in Python: the row's Cd is metering_cd_and_dRe at
    the pipe Reynolds number and D-D/2 taps, the taps ISO's loss relation is
    derived for."""
    g = _between_pipes(PlenumNode, 0.5, "incompressible", correlation)
    _, r = _solve(g)
    geom = g.elements["o"]._orifice_geom
    assert geom.tap == cb.OrificeTap.D_D2
    mu = cb.transport_state(300.0, r["a.P"], cb.species.dry_air()).mu
    Re_D = 4.0 * r["o.m_dot"] / (math.pi * D * mu)
    C = cb.metering_cd_and_dRe(getattr(cb.MeteringCdCorrelation, correlation), geom, Re_D)[0]
    assert r["o.Cd"] == pytest.approx(C, rel=1e-12)
    k = _k_by_hand(0.5, C)
    assert r["a.Pt"] - r["b.Pt"] == pytest.approx(k * _iso_dp(r, 0.5, C), rel=1e-4)


def test_a_reversed_orifice_recovers_into_the_pipe_it_flows_into() -> None:
    """Declared a -> b, the flow runs b -> a: the approach pipe is p2, the
    recovering pipe p1, and the loss is still k dp (#554's beta_reverse)."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=1.0e5, Tt=300.0))
    g.add_node(PlenumNode("a"))
    g.add_node(PlenumNode("b"))
    g.add_node(PressureBoundary("ho", Pt=1.2e5, Tt=300.0))
    g.add_element(ChannelElement("p1", "hi", "a", length=0.5, diameter=D, regime="incompressible"))
    g.add_element(
        OrificeElement(
            "o", "a", "b", Cd=0.6, diameter=0.06, regime="incompressible", correlation="fixed"
        )
    )
    g.add_element(ChannelElement("p2", "b", "ho", length=0.5, diameter=D, regime="incompressible"))
    _, r = _solve(g)
    assert r["o.m_dot"] < 0.0
    rho = cb.density(300.0, r["o.P_approach"], cb.species.dry_air())
    E = 1.0 / math.sqrt(1.0 - 0.6**4)
    A = math.pi * 0.06**2 / 4.0
    dp = (r["o.m_dot"] / (0.6 * E * A)) ** 2 / (2.0 * rho)
    assert r["b.Pt"] - r["a.Pt"] == pytest.approx(_k_by_hand(0.6, 0.6) * dp, rel=1e-4)


_HOLES = cb.DischargeCdCorrelation
# (correlation, its C++ selector, geometry as a fraction of the bore)
_HOLE_CASES = [
    ("IdelchikBeveled", _HOLES.Idelchik1966Beveled, {"bevel_depth": 0.04}),
    ("IdelchikRounded", _HOLES.Idelchik1966Rounded, {"edge_radius": 0.04}),
    ("IdelchikRounded", _HOLES.Idelchik1966Rounded, {"edge_radius": 0.16}),
]


def _hole_between_pipes(
    correlation: str, geom: dict[str, float], beta: float, D_pipe: float = D
) -> FlowNetwork:
    """pipe -> plenum -> hole -> plenum -> pipe, incompressible like
    Idelchik's Section IV."""
    d = beta * D_pipe
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=1.05e5, Tt=300.0))
    g.add_node(PlenumNode("a"))
    g.add_node(PlenumNode("b"))
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
    reg = "incompressible"
    g.add_element(ChannelElement("p1", "hi", "a", length=0.05, diameter=D_pipe, regime=reg))
    g.add_element(
        OrificeElement(
            "o",
            "a",
            "b",
            diameter=d,
            correlation=correlation,
            regime=reg,
            **{k: v * d for k, v in geom.items()},
        )
    )
    g.add_element(ChannelElement("p2", "b", "ho", length=0.05, diameter=D_pipe, regime=reg))
    return g


def _hole(d: float, geom: dict[str, float]) -> cb.DischargeHoleGeometry:
    h = cb.DischargeHoleGeometry(
        d=d, L=geom.get("plate_thickness", 0.0) * d, r=geom.get("edge_radius", 0.0) * d
    )
    h.bevel = geom.get("bevel_depth", 0.0) * d
    return h


def _zeta_w1(r: dict, D_pipe: float) -> float:
    """The solved permanent loss on the approach pipe's velocity head, at its
    own static density (#651)."""
    rho = cb.density(300.0, r["o.P_approach"], cb.species.dry_air())
    w1 = r["o.m_dot"] / (rho * math.pi * D_pipe * D_pipe / 4.0)
    return (r["a.Pt"] - r["b.Pt"]) / (0.5 * rho * w1 * w1)


@pytest.mark.parametrize(("correlation", "selector", "geom"), _HOLE_CASES)
@pytest.mark.parametrize("beta", [0.3, 0.5, 0.7])
def test_a_hole_between_pipes_charges_idelchiks_in_pipe_loss(
    correlation: str, selector: object, geom: dict[str, float], beta: float
) -> None:
    """Between pipes a discharge hole loses Idelchik's in-pipe zeta1 (4-15,
    4-16), not its whole jet head: the plenum Cd carried into the pipe through
    the recovery row (#645). Beveled and rounded have no Re term, so the
    comparison is direct. It used to be 16-20% low in flow at beta 0.5."""
    _, r = _solve(_hole_between_pipes(correlation, geom, beta))
    expected = cb.discharge_pipe_zeta(selector, _hole(beta * D, geom), beta)
    # The old behaviour was 50-160% off in zeta.
    assert _zeta_w1(r, D) == pytest.approx(expected, rel=2e-3)


@pytest.mark.parametrize("ld", [0.0, 1.0])
@pytest.mark.parametrize("beta", [0.3, 0.5, 0.7])
def test_a_thick_hole_between_pipes_is_idelchiks_at_high_re(ld: float, beta: float) -> None:
    """The thick hole's wall Cd keeps its Re dependence; at Re >= 1e6, where
    diagram 4-17/4-18a reach their high-Re value, the in-pipe loss is 4-14's."""
    D_pipe = 1.0
    _, r = _solve(_hole_between_pipes("IdelchikThick", {"plate_thickness": ld}, beta, D_pipe))
    hole = _hole(beta * D_pipe, {"plate_thickness": ld})
    expected = cb.discharge_pipe_zeta(_HOLES.Idelchik1966Thick, hole, beta)
    # Exact for the short hole. With length, the wall Cd's friction is the
    # smooth-pipe lam(Re) (0.012 here) where 4-14's table is computed at
    # lam = 0.02: 0.5-0.8% less loss at l/d 1, and only that.
    assert _zeta_w1(r, D_pipe) == pytest.approx(expected, rel=1e-4 if ld == 0.0 else 1e-2)
    if ld > 0.0:
        assert _zeta_w1(r, D_pipe) < expected


@pytest.mark.parametrize("correlation", ["IdelchikThick", "McGreehanSchotsch", "Lichtarowicz"])
def test_a_hole_recovers_and_reports_its_in_pipe_cd_factor(correlation: str) -> None:
    """Every discharge-hole correlation recovers into the far pipe, with the
    in-pipe factor of its own approach pipe in each direction."""
    g = _between_pipes(PlenumNode, 0.5, "incompressible", correlation)
    o = g.elements["o"]
    o.plate_thickness = 0.02
    g.resolve_all_topology()
    assert o._recovery_fwd and o._recovery_rev
    sel = {
        "IdelchikThick": _HOLES.Idelchik1966Thick,
        "McGreehanSchotsch": _HOLES.McGreehanSchotsch1988,
        "Lichtarowicz": _HOLES.Lichtarowicz1965,
    }[correlation]
    R = cb.discharge_pipe_cd_ratio(sel, o._hole_cd.discharge_hole(), 0.5)
    assert R != 1.0
    assert o._pipe_cd_ratio_fwd == R
    assert o._pipe_cd_ratio_rev == R


def test_a_reversed_hole_takes_the_factor_of_the_pipe_it_flows_from() -> None:
    """Reversed, the approach is the to-node's pipe (#554): its beta sets the
    hole's in-pipe factor, which the diagnostics report along the flow."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=1.0e5, Tt=300.0))
    g.add_node(PlenumNode("a"))
    g.add_node(PlenumNode("b"))
    g.add_node(PressureBoundary("ho", Pt=1.05e5, Tt=300.0))
    reg = "incompressible"
    g.add_element(ChannelElement("p1", "hi", "a", length=0.05, diameter=0.1, regime=reg))
    g.add_element(
        OrificeElement(
            "o",
            "a",
            "b",
            diameter=0.05,
            correlation="IdelchikRounded",
            edge_radius=0.002,
            regime=reg,
        )
    )
    g.add_element(ChannelElement("p2", "b", "ho", length=0.05, diameter=0.07, regime=reg))
    _, r = _solve(g)
    assert r["o.m_dot"] < 0.0
    hole = g.elements["o"]._hole_cd.discharge_hole()
    R_rev = cb.discharge_pipe_cd_ratio(_HOLES.Idelchik1966Rounded, hole, 0.05 / 0.07)
    assert r["o.Cd_pipe_ratio"] == pytest.approx(R_rev, rel=1e-14)
    assert R_rev != pytest.approx(
        cb.discharge_pipe_cd_ratio(_HOLES.Idelchik1966Rounded, hole, 0.5), rel=1e-3
    )
    # And the row charges it: the loss on the reverse approach pipe's head is
    # 4-16's zeta1 at that pipe's beta.
    rho = cb.density(300.0, r["o.P_approach"], cb.species.dry_air())
    w1 = abs(r["o.m_dot"]) / (rho * math.pi * 0.07**2 / 4.0)
    zeta = (r["b.Pt"] - r["a.Pt"]) / (0.5 * rho * w1 * w1)
    expected = cb.discharge_pipe_zeta(_HOLES.Idelchik1966Rounded, hole, 0.05 / 0.07)
    assert zeta == pytest.approx(expected, rel=2e-3)


def test_a_plenum_fed_hole_keeps_its_wall_cd() -> None:
    """No approach pipe: no factor, no recovery -- the plenum Cd as before."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=1.05e5, Tt=300.0))
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
    g.add_element(OrificeElement("o", "hi", "ho", diameter=0.01, correlation="IdelchikThick"))
    g.resolve_all_topology()
    o = g.elements["o"]
    assert o._pipe_cd_ratio_fwd == 1.0 and o._pipe_cd_ratio_rev == 1.0
    assert not o._recovery_fwd and not o._recovery_rev


def test_tap_must_be_named() -> None:
    with pytest.raises(ValueError, match="tap 'Vena' is not one of"):
        OrificeElement("o", "a", "b", diameter=0.02, tap="Vena")


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.parametrize("node_cls", [PlenumNode, MomentumChamberNode])
@pytest.mark.parametrize(
    "correlation", ["fixed", "ReaderHarrisGallagher", "IdelchikThick", "IdelchikRounded"]
)
def test_the_recovering_rows_jacobian_is_exact(
    regime: str, node_cls: type, correlation: str
) -> None:
    """Every column of the orifice row against central differences, in the
    J_ij |x_j| metric: the tap pressure's Pt and P columns and, with
    Reader-Harris/Gallagher, dk/dC through the flow column. A discharge hole
    adds its Re-dependent Cd times the in-pipe factor (#645)."""
    s, r = _solve(_between_pipes(node_cls, 0.7, regime, correlation))
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()

    def F(v: np.ndarray) -> np.ndarray:
        return np.asarray(s._residuals_and_jacobian(v, compute_jacobian=False)[0])

    rows = s._row_kinds(x)
    i = next(k for k, (owner, _) in enumerate(rows) if "'o'" in owner)
    scale = np.abs(J[i]) * np.abs(x)
    checked = 0
    for j, name in enumerate(s.unknown_names):
        h = 1e-6 * max(abs(x[j]), 1e-3)
        e = np.zeros_like(x)
        e[j] = h
        fd = (F(x + e)[i] - F(x - e)[i]) / (2 * h)
        assert J[i, j] * abs(x[j]) == pytest.approx(
            fd * abs(x[j]), rel=1e-5, abs=1e-7 * scale.max()
        ), name
        checked += abs(J[i, j]) > 0.0
    assert checked >= 3
