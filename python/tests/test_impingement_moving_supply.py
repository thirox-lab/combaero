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
) -> FlowNetwork:
    """N rows fed by a supply channel (CrossflowSegmentElement chain) that
    runs along the jet plate from ``src`` to ``sink``; their gap drains to
    ``exit``. A small ``area`` with the sink above the stations feeds the
    channel from both ends."""
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
    with pytest.raises(ValueError, match="crossflow-chain station"):
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


def test_the_jacobian_is_analytic_where_the_supply_reverses() -> None:
    # Fed from both ends: the channel flow reverses between two stations.
    s, r = _solve(_channel_net(area=1.2e-4))
    assert r["s3x.m_dot"] < 0.0 < r["s0x.m_dot"]
    x = np.array(r["__x_solution__"])
    assert _fd_rows(s, x) < 1e-5
    assert _fd_rows(s, _perturbed(s, x, 9)) < 1e-5
