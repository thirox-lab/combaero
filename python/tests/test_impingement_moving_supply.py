"""Impingement rows fed by a moving supply (C2, #1087).

An ImpingementPlateElement whose supply node is a crossflow-chain station
(a PlenumNode carrying the chain's _CrossflowStation) is duct-fed: its jets
drive from the station's rebuilt duct total, the composed Cd (base, jet law,
cap) referenced to P1 + f q1 and converted onto it, f =
HoleCd.supply_recovery (None = 0). A row fed by a plenum is unchanged, a row
fed by a MomentumChamberNode drives from the chamber's Pt and refuses f.

No fidelity set exists for this configuration (Florschuetz 1981 and CR-3630
are plenum-fed); its accuracy is unscored. These tests pin the wiring, the
analytic Jacobian and the diagnostics only.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian, scaled_error

import combaero as cb
from combaero.network import (
    CrossflowSegmentElement,
    FlowNetwork,
    HoleCd,
    ImpingementCrossflowElement,
    ImpingementPlateElement,
    LosslessConnectionElement,
    MomentumChamberNode,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)

D = 0.00254
SPAN = 0.122
Z_D = 2.0
N = 3
Y_AIR = cb.mole_to_mass(cb.species.dry_air())


def _pb(name: str, Pt: float, Tt: float = 300.0) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y_AIR
    return b


def _row(i: int, supply: str, spec: HoleCd | None) -> ImpingementPlateElement:
    kw = {} if spec is None else {"hole_cd": spec}
    return ImpingementPlateElement(
        f"p{i}", supply, f"c{i}", d_jet=D, xn_d=5.0, yn_d=4.0, z_d=Z_D, span=SPAN,
        plate_thickness=D, row=i, **kw,
    )  # fmt: skip


def _gap(net: FlowNetwork, regime: str) -> None:
    dx = 5.0 * D
    for i in range(1, N + 1):
        net.add_element(
            ImpingementCrossflowElement(
                f"x{i}", f"c{i}", f"c{i + 1}" if i < N else "exit", length=dx,
                height=Z_D * D, span=SPAN, regime=regime,
            )
        )  # fmt: skip


def _channel_net(
    spec: HoleCd | None = None,
    regime: str = "compressible",
    area: float = 1e-3,
    p_sink: float = 1.02e5,
    bypass_P: float | None = None,
) -> FlowNetwork:
    """N rows fed by a supply channel (CrossflowSegmentElement chain) that
    runs along the jet plate from ``src`` to ``sink``; their gap drains to
    ``exit``. A small ``area`` with the sink above the stations feeds the
    channel from both ends. ``bypass_P``: an initial (bypass) crossflow
    into the gap ahead of row 1 from a boundary at that pressure, strong
    enough for Rogers & Hersh's law to act (below its cap)."""
    net = FlowNetwork()
    net.add_node(_pb("src", 1.03e5))
    net.add_node(_pb("sink", p_sink))
    net.add_node(_pb("exit", 1.0e5))
    bleed = cb.STATION_KAPPA_BLEED_BASSETT
    dx = 5.0 * D
    for i in range(1, N + 1):
        net.add_node(PlenumNode(f"s{i}"))
        net.add_node(PlenumNode(f"c{i}"))
    net.add_element(
        CrossflowSegmentElement(
            "s0x", "src", "s1", length=0.5 * dx, area=area, entry_K=0.0, to_kappa=bleed,
            next_seg="s1x", regime=regime,
        )
    )  # fmt: skip
    for i in range(1, N + 1):
        last = i == N
        net.add_element(
            CrossflowSegmentElement(
                f"s{i}x", f"s{i}", "sink" if last else f"s{i + 1}",
                length=(0.5 if last else 1.0) * dx, area=area, from_kappa=bleed,
                prev_seg=f"s{i - 1}x", to_kappa=None if last else bleed,
                next_seg=None if last else f"s{i + 1}x", regime=regime,
            )
        )  # fmt: skip
        net.add_element(_row(i, f"s{i}", spec))
    _gap(net, regime)
    if bypass_P is not None:
        net.add_node(_pb("byp", bypass_P))
        net.add_element(
            ImpingementCrossflowElement(
                "x0", "byp", "c1", length=dx, height=Z_D * D, span=SPAN, regime=regime,
                bypass_inlet=True,
            )
        )  # fmt: skip
    for e in net.elements.values():
        if hasattr(e, "regime"):
            e.regime = regime
    return net


def _plenum_net(spec: HoleCd | None = None) -> FlowNetwork:
    net = FlowNetwork()
    net.add_node(_pb("sup", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    for i in range(1, N + 1):
        net.add_node(PlenumNode(f"c{i}"))
        net.add_element(_row(i, "sup", spec))
    _gap(net, "compressible")
    return net


def _chamber_net(spec: HoleCd | None = None) -> FlowNetwork:
    """One row fed by a MomentumChamberNode (a chamber feeds one element:
    its through-flow is the row's)."""
    net = FlowNetwork()
    net.add_node(_pb("src", 1.03e5))
    net.add_node(MomentumChamberNode("ch", area=1e-4))
    net.add_node(_pb("exit", 1.0e5))
    net.add_node(PlenumNode("c1"))
    net.add_element(LosslessConnectionElement("feed", "src", "ch"))
    net.add_element(_row(1, "ch", spec))
    net.add_element(
        ImpingementCrossflowElement("x1", "c1", "exit", length=5.0 * D, height=Z_D * D, span=SPAN)
    )
    return net


def _solve(net: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    return s, r


def _rows(s: NetworkSolver) -> list[int]:
    return [s.unknown_names.index(f"p{i}.m_dot") for i in range(1, N + 1)]


def _fd_rows(s: NetworkSolver, x: np.ndarray) -> float:
    J = s._residuals_and_jacobian(x, compute_jacobian=True)[1].toarray()
    return float(scaled_error(J, fd_jacobian(s, x), x)[_rows(s)].max())


def _perturbed(s: NetworkSolver, x: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    xp = x.copy()
    for j, name in enumerate(s.unknown_names):
        if name.endswith((".P", ".Pt")):
            xp[j] += 30.0 * rng.normal()
        else:
            xp[j] *= 1.0 + 0.03 * rng.normal()
    return xp


# --- the wiring ---------------------------------------------------------------


def test_a_row_fed_by_a_crossflow_station_is_duct_fed() -> None:
    s, r = _solve(_channel_net())
    assert r["__success__"]
    for i in range(1, N + 1):
        d = r["__element_diag__"][f"p{i}"]
        assert d["Cd_supply_drive_model"] == "recovery"
        assert d["Cd_supply_drive_applied"] == 1.0
        assert d["Cd_supply_recovery_fraction"] == 0.0
        assert d["Cd_supply_recovery_q1"] > 0.0
        assert 0.0 < d["Cd_supply_recovery_ratio"] < 1.0
        # The reported Cd is on the total basis: the composed Cd times it.
        assert d["Cd"] == pytest.approx(d["Cd_base"])
        assert d["Cd"] == pytest.approx(0.79 * d["Cd_supply_recovery_ratio"], rel=1e-12)
        assert d["Cd_inlet_crossflow_ignored"] == 1.0
        assert d["surface_extrapolated_supply_crossflow"] == 1.0
        assert d["surface_extrapolated"] == 1.0


def test_the_drive_is_the_stations_duct_total() -> None:
    s, r = _solve(_channel_net())
    x = np.array(r["__x_solution__"])
    el = s.network.elements["p2"]
    si = s._get_node_state(s.network.nodes["s2"], x)
    so = s._get_node_state(s.network.nodes["c2"], x)
    flows = s._network_flows(el, x)
    pair = el._hole_cd.drive_pair(si, so, flows)
    station = s.network.nodes["s2"]._crossflow_station
    assert el._crossflow_station is station
    ms = station.moving_state(si, flows, incompressible=False)
    assert pair.p_up == pytest.approx(ms.Pt, rel=1e-12)
    # f = 0: the flow passes zero at the station static.
    assert pair.p_rev == float(si.P)
    assert {"s1x", "s2x"} <= set(flows)


@pytest.mark.parametrize("f", [None, 0.5])
def test_rogers_hersh_reads_the_kernels_referenced_state(f: float | None) -> None:
    # A bypass crossflow makes Rogers & Hersh's law act (below its cap) on
    # every row. Its G_ideal is the ideal jet at the conversion kernel's own
    # referenced state (Pf = P1 + f q1, Tf) to the gap static: at f = 0 the
    # station's static P and static T, not the supply node's (Pt, Tt).
    spec = None if f is None else HoleCd(supply_recovery=f)
    s, r = _solve(_channel_net(spec, bypass_P=1.03e5))
    assert r["__success__"]
    x = np.array(r["__x_solution__"])
    acts = [
        i for i in range(1, N + 1) if r["__element_diag__"][f"p{i}"]["Cd_exit_crossflow_applied"]
    ]
    # Every row at f = 0, row 1 at f = 0.5 (the law caps downstream there).
    assert acts == ([1, 2, 3] if f is None else [1])
    for i in acts:
        el = s.network.elements[f"p{i}"]
        si = s._get_node_state(s.network.nodes[f"s{i}"], x)
        so = s._get_node_state(s.network.nodes[f"c{i}"], x)
        flows = s._network_flows(el, x)
        pair = el._hole_cd.drive_pair(si, so, flows)
        T1 = el._crossflow_station.static_T(si, flows).T
        c = cb.hole_supply_recovery_cd_ratio(
            float(si.P), T1, pair.p_up, float(si.Tt), float(so.P), list(si.Y), f or 0.0
        )
        assert c.Pf == pytest.approx(pair.p_rev, rel=1e-14)
        if f is None:
            assert c.Pf == float(si.P) and c.Tf == T1
        rho = float(cb.density(c.Tf, c.Pf, list(si.X)))
        G = math.sqrt(2.0 * rho * (c.Pf - float(so.P)))
        jc = el._hole_cd._jet_drive(si, so, flows)
        assert jc.forward and jc.G_i == pytest.approx(G, rel=1e-12)
        # Not the plenum frame (the node's Pt and Tt).
        rho_pl = float(cb.density(float(si.Tt), float(si.Pt), list(si.X)))
        G_pl = math.sqrt(2.0 * rho_pl * (float(si.Pt) - float(so.P)))
        assert abs(jc.G_i / G_pl - 1.0) > 1e-4


def test_the_referenced_temperature_has_analytic_partials() -> None:
    Y = list(Y_AIR)
    args = {"P1": 1.2e5, "T1": 430.0, "Pt1": 1.32e5, "Tt1": 450.0, "P2": 1.05e5}
    for f in (0.0, 0.3, 1.0):
        c = cb.hole_supply_recovery_cd_ratio(**args, Y=Y, f=f)
        for key, dname in (("P1", "dTf_dP1"), ("T1", "dTf_dT1"), ("Pt1", "dTf_dPt1"),
                           ("Tt1", "dTf_dTt1")):  # fmt: skip
            h = 1e-6 * args[key]
            up = cb.hole_supply_recovery_cd_ratio(**{**args, key: args[key] + h}, Y=Y, f=f)
            dn = cb.hole_supply_recovery_cd_ratio(**{**args, key: args[key] - h}, Y=Y, f=f)
            fd = (up.Tf - dn.Tf) / (2.0 * h)
            assert getattr(c, dname) == pytest.approx(fd, rel=1e-6, abs=1e-12)
        if 0.0 < f < 1.0:
            up = cb.hole_supply_recovery_cd_ratio(**args, Y=Y, f=f + 1e-6)
            dn = cb.hole_supply_recovery_cd_ratio(**args, Y=Y, f=f - 1e-6)
            assert c.dTf_df == pytest.approx((up.Tf - dn.Tf) / 2e-6, rel=1e-6)
    assert cb.hole_supply_recovery_cd_ratio(**args, Y=Y, f=0.0).Tf == 430.0
    assert cb.hole_supply_recovery_cd_ratio(**args, Y=Y, f=1.0).Tf == 450.0
    # Set where the ratio is 0 (Pf at or below P2) too.
    low = cb.hole_supply_recovery_cd_ratio(**{**args, "P2": 1.25e5}, Y=Y, f=0.0)
    assert low.ratio == 0.0 and low.Tf == 430.0


def test_f_one_takes_the_duct_total_unconverted() -> None:
    _, r0 = _solve(_channel_net())
    _, r1 = _solve(_channel_net(HoleCd(supply_recovery=1.0)))
    for i in range(1, N + 1):
        d = r1["__element_diag__"][f"p{i}"]
        assert d["Cd_supply_recovery_ratio"] == 1.0
        assert d["Cd"] == pytest.approx(0.79, rel=1e-12)
        assert r1[f"p{i}.m_dot"] > r0[f"p{i}.m_dot"]


def test_a_plenum_fed_row_is_unchanged_and_refuses_f() -> None:
    _, r = _solve(_plenum_net())
    d = r["__element_diag__"]["p1"]
    assert d["Cd_supply_drive_applied"] == 0.0
    assert d["Cd_supply_recovery_q1"] == 0.0
    assert math.isnan(d["Cd_supply_recovery_fraction"])
    assert d["surface_extrapolated_supply_crossflow"] == 0.0
    with pytest.raises(ValueError, match="is a PressureBoundary that no chain registers"):
        _solve(_plenum_net(HoleCd(supply_recovery=0.0)))


def test_a_chamber_fed_row_drives_from_the_chambers_total_and_refuses_f() -> None:
    s, r = _solve(_chamber_net())
    assert r["__success__"]
    x = np.array(r["__x_solution__"])
    el = s.network.elements["p1"]
    assert el._crossflow_station is None and not el._hole_cd.duct_fed
    si = s._get_node_state(s.network.nodes["ch"], x)
    so = s._get_node_state(s.network.nodes["c1"], x)
    assert el._hole_cd.drive_pair(si, so, s._network_flows(el, x)) is None
    d = r["__element_diag__"]["p1"]
    assert d["Cd_supply_drive_applied"] == 0.0
    assert d["surface_extrapolated_supply_crossflow"] == 0.0
    for f in (0.0, 1.0):
        with pytest.raises(ValueError, match="MomentumChamberNode"):
            _solve(_chamber_net(HoleCd(supply_recovery=f)))


def test_the_static_drive_is_refused_on_a_row() -> None:
    with pytest.raises(ValueError, match="supply_drive='static'"):
        _row(1, "s1", HoleCd(supply_drive="static"))


# --- the analytic Jacobian ----------------------------------------------------


@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
@pytest.mark.parametrize(
    "spec",
    [
        None,
        HoleCd(supply_recovery=0.5),
        HoleCd(exit_crossflow="florschuetz_metzger_su_1982"),
        HoleCd(exit_crossflow="florschuetz_metzger_su_1982", supply_recovery=0.5),
    ],
    ids=["f0", "f05", "cr3630_f0", "cr3630_f05"],
)
def test_the_rows_jacobian_is_analytic(regime: str, spec: HoleCd | None) -> None:
    s, r = _solve(_channel_net(spec, regime))
    x = np.array(r["__x_solution__"])
    assert _fd_rows(s, x) < 1e-5
    assert _fd_rows(s, _perturbed(s, x, 5)) < 1e-5


@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
@pytest.mark.parametrize("f", [None, 0.5])
def test_the_jacobian_is_analytic_where_rogers_hersh_acts(regime: str, f: float | None) -> None:
    # G_ideal through the drive pair and the referenced (Pf, Tf): the jet
    # law's slope reaches the station's P, Tt, Y and segment flows.
    spec = None if f is None else HoleCd(supply_recovery=f)
    s, r = _solve(_channel_net(spec, regime, bypass_P=1.03e5))
    assert r["__element_diag__"]["p1"]["Cd_exit_crossflow_applied"] == 1.0
    x = np.array(r["__x_solution__"])
    assert _fd_rows(s, x) < 1e-5
    assert _fd_rows(s, _perturbed(s, x, 5)) < 1e-5


def test_the_jacobian_is_analytic_where_the_supply_reverses() -> None:
    # Fed from both ends: the channel flow reverses between two stations.
    s, r = _solve(_channel_net(area=1.2e-4))
    assert r["s3x.m_dot"] < 0.0 < r["s0x.m_dot"]
    x = np.array(r["__x_solution__"])
    assert _fd_rows(s, x) < 1e-5
    assert _fd_rows(s, _perturbed(s, x, 9)) < 1e-5
