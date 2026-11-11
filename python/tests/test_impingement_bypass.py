"""An impingement array with an initial (bypass) crossflow (#467).

The bypass enters the gap's upstream end through a lead-in segment declared
``bypass_inlet``. The claims here are about the wiring and the momentum
hand-off; the flow distribution against Florschuetz, Metzger and Su (1982,
NASA CR-3630) Table D.1 is in test_impingement_bypass_validation.py.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize._numdiff import approx_derivative

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    HoleCd,
    ImpingementArray,
    ImpingementCrossflowElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)

D = 0.00254
SPAN = 12 * 4.0 * D
Y_AIR = cb.mole_to_mass(cb.species.dry_air())


def _pb(name: str, Pt: float, Tt: float = 300.0) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y_AIR
    return b


def _array(id: str = "ia", n_rows: int = 10, **kw) -> ImpingementArray:
    return ImpingementArray(
        id, n_rows=n_rows, d_jet=D, xn_d=5.0, yn_d=4.0, z_d=2.0, span=SPAN, plate_thickness=D, **kw
    )


def _solve(net: FlowNetwork) -> dict:
    res = NetworkSolver(net).solve()
    assert res["__success__"], res.get("__message__")
    return res


def _rows(res: dict, arr: ImpingementArray, key: str = "m_dot") -> list[float]:
    d = res["__element_diag__"]
    return [d[arr.plate_id(i)][key] for i in range(1, int(arr.n_rows) + 1)]


def _with_bypass(m_c: float, node: str = "chamber", Tt: float = 300.0, **kw):
    """Supply and exit at fixed pressure; the bypass a fixed flow through a
    feed duct into ``node`` (a momentum chamber or a plenum)."""
    net = FlowNetwork()
    net.add_node(_pb("supply", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    net.add_node(MassFlowBoundary("feed_in", m_dot=m_c, Tt=Tt, Y=Y_AIR))
    net.add_node(MomentumChamberNode("byp") if node == "chamber" else PlenumNode("byp"))
    net.add_element(ChannelElement("feed", "feed_in", "byp", length=0.2, diameter=0.02))
    arr = _array(**kw)
    arr.add_to(net, "supply", "exit", bypass="byp")
    return net, arr


def _closed() -> tuple[dict, ImpingementArray]:
    net = FlowNetwork()
    net.add_node(_pb("supply", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    arr = _array()
    arr.add_to(net, "supply", "exit")
    return _solve(net), arr


def test_the_lead_in_is_declared_and_half_a_pitch_long() -> None:
    net, arr = _with_bypass(0.001)
    x0 = net.elements["ia__x0"]
    assert isinstance(x0, ImpingementCrossflowElement) and x0.bypass_inlet
    assert (x0.from_node, x0.to_node) == ("byp", "ia__c1")
    assert x0.length == pytest.approx(0.5 * 5.0 * D)  # CR-3630: L starts xn/2 ahead of row 1
    assert x0.area == pytest.approx(net.elements["ia__x1"].area)


def test_a_zero_bypass_is_the_closed_end() -> None:
    """The regression that keeps #465 intact: an inlet carrying nothing is
    Florschuetz et al.'s (1981) closed end, row by row."""
    res_c, arr_c = _closed()
    net, arr = _with_bypass(0.0)
    res = _solve(net)
    for a, b in zip(_rows(res, arr), _rows(res_c, arr_c), strict=True):
        assert a == pytest.approx(b, rel=1e-6)


def test_row_one_sees_the_bypass_as_its_crossflow() -> None:
    net, arr = _with_bypass(0.004)
    res = _solve(net)
    d = res["__element_diag__"]
    m1 = d[arr.plate_id(1)]["m_dot"]
    gc = 0.004 / m1 * (3.141592653589793 / 4.0) / (4.0 * 2.0)
    assert d[arr.plate_id(1)]["Gc_Gj"] == pytest.approx(gc, rel=1e-6)
    s = arr.summarize(d)
    assert s["m_dot_bypass"] == pytest.approx(0.004, rel=1e-6)
    # The bypass starves the upstream rows: they discharge against it.
    res_c, arr_c = _closed()
    assert _rows(res, arr)[0] < _rows(res_c, arr_c)[0]


def test_starved_rows_are_reported_unresolved() -> None:
    """A row whose plate drop is below 1% of the gap's dynamic head is below
    the crossflow chain's resolution (#729). Closed, no row is; a bypass that
    starves the upstream rows marks them, and only them. With the holes' Cd
    held constant (``HoleCd(exit_crossflow='none')``): under the crossflow Cd (#800)
    the same starved rows need a drop well above that (next test)."""
    res_c, arr_c = _closed()
    s = arr_c.summarize(res_c["__element_diag__"])
    assert s["n_rows_unresolved"] == 0.0
    net, arr = _with_bypass(0.05, hole_cd=HoleCd(exit_crossflow="none"))
    d = _solve(net)["__element_diag__"]
    s = arr.summarize(d)
    flagged = [i for i, f in enumerate(s["rows_unresolved"], 1) if f]
    # Measured: rows 1-5 at 0.02-0.41% of the head, row 6 at 5.8%.
    assert flagged == [1, 2, 3, 4, 5]
    m = s["rows_m_dot"]
    assert max(abs(m[i - 1]) for i in flagged) < 0.1 * max(m)


def test_under_the_crossflow_cd_the_starved_rows_are_resolved() -> None:
    """The crossflow Cd (#800, default) makes a starved row's flow go as
    dP^1.07: the small flow now takes a drop above the chain's resolution, so
    no row is flagged, and the rows' flows rise monotonically from the
    upstream end (none reversed)."""
    net, arr = _with_bypass(0.05)
    d = _solve(net)["__element_diag__"]
    s = arr.summarize(d)
    assert s["n_rows_unresolved"] == 0.0
    m = s["rows_m_dot"]
    assert all(b > a for a, b in zip(m, m[1:], strict=False))
    assert m[0] > 0.0


def test_the_feed_node_does_not_matter() -> None:
    """The lead-in accelerates the bypass from the feed node's TOTAL
    pressure, and a channel keeps its total pressure into either node. So a
    momentum chamber and a plenum feed the gap alike: no dynamic head is
    lost at the port, whichever the user (or the GUI) puts there."""
    res_ch = _solve(_with_bypass(0.004, "chamber")[0])
    res_pl = _solve(_with_bypass(0.004, "plenum")[0])
    arr = _array()
    for a, b in zip(_rows(res_ch, arr), _rows(res_pl, arr), strict=True):
        assert a == pytest.approx(b, rel=1e-6)
    assert res_pl["feed_in.P"] == pytest.approx(res_ch["feed_in.P"], rel=1e-9)
    assert res_pl["byp.Pt"] == pytest.approx(res_ch["byp.Pt"], rel=1e-9)


def test_an_entry_loss_costs_feed_pressure_not_distribution() -> None:
    """``bypass_entry_K`` (an abrupt entry from a large plenum) is paid by
    the feed, at the same bypass flow."""
    res_0 = _solve(_with_bypass(0.004)[0])
    res_k = _solve(_with_bypass(0.004, bypass_entry_K=0.5)[0])
    arr = _array()
    for a, b in zip(_rows(res_0, arr), _rows(res_k, arr), strict=True):
        assert a == pytest.approx(b, rel=1e-6)
    assert res_k["feed_in.P"] > res_0["feed_in.P"]


def test_arrays_in_series_continue_the_gap() -> None:
    """A's exit into B's bypass: the stream keeps its momentum. With a
    vanishing lead-in, two five-row arrays are one ten-row array."""
    net = FlowNetwork()
    net.add_node(_pb("supply", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    net.add_node(PlenumNode("joint"))
    a = _array("a", n_rows=5)
    b = _array("b", n_rows=5, bypass_length=1.0e-9)
    a.add_to(net, "supply", "joint")
    b.add_to(net, "supply", "exit", bypass="joint")
    res = _solve(net)
    x0 = net.elements["b__x0"]
    assert x0.entry_K is None and x0.from_kappa is not None and x0.prev_seg == "a__x5"
    res_c, arr_c = _closed()
    chained = _rows(res, a) + _rows(res, b)
    # Measured 1.6e-5: the joint is one more node with its own density.
    for got, want in zip(chained, _rows(res_c, arr_c), strict=True):
        assert got == pytest.approx(want, rel=1e-4)


def test_gaps_in_series_must_match() -> None:
    net = FlowNetwork()
    net.add_node(_pb("supply", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    net.add_node(PlenumNode("joint"))
    _array("a", n_rows=3).add_to(net, "supply", "joint")
    ImpingementArray(
        "b", n_rows=3, d_jet=D, xn_d=5.0, yn_d=4.0, z_d=3.0, span=SPAN, plate_thickness=D
    ).add_to(net, "supply", "exit", bypass="joint")
    with pytest.raises(ValueError, match="cross-section differs"):
        NetworkSolver(net).solve()


def test_an_undeclared_inflow_is_still_refused() -> None:
    """Only a declared inlet admits a bypass: a stream wired straight into a
    row's crossflow node is still refused."""
    net = FlowNetwork()
    net.add_node(_pb("supply", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    net.add_node(_pb("side", 1.02e5))
    arr = _array(n_rows=3)
    arr.add_to(net, "supply", "exit")
    net.add_element(ChannelElement("side_in", "side", "ia__c2", length=0.01, diameter=0.01))
    with pytest.raises(ValueError, match="bypass_inlet"):
        NetworkSolver(net).solve()


def test_a_bypass_temperature_inside_the_data_is_not_flagged() -> None:
    """A hotter bypass moves each row's reference temperature through
    CR-3936's eta_r (test_impingement_eta_r.py); inside the inline data that
    is in-source, so nothing is flagged."""
    net, arr = _with_bypass(0.004)
    s = arr.summarize(_solve(net)["__element_diag__"])
    assert s["heat_transfer_outside_source"] == 0.0
    net, arr = _with_bypass(0.004, Tt=400.0)
    s = arr.summarize(_solve(net)["__element_diag__"])
    assert s["bypass_Tt_minus_supply"] == pytest.approx(100.0, abs=1.0)
    assert s["heat_transfer_outside_source"] == 0.0


def test_bad_bypass_parameters_are_refused() -> None:
    with pytest.raises(ValueError, match="bypass_length"):
        _array(bypass_length=0.0)
    with pytest.raises(ValueError, match="bypass_entry_K"):
        _array(bypass_entry_K=-0.1)


def test_a_reversed_bypass_is_flagged() -> None:
    """Feed pressure below the gap: the gap drains back out of the inlet.
    Outside CR-3630: the lead-in then discharges at static continuity (the
    entry term is inflow only, C1 through zero); the summary says it happened."""
    net = FlowNetwork()
    net.add_node(_pb("supply", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    net.add_node(_pb("bsrc", 1.0005e5))
    net.add_node(PlenumNode("byp"))
    net.add_element(ChannelElement("feed", "bsrc", "byp", length=0.1, diameter=0.02))
    arr = _array()
    arr.add_to(net, "supply", "exit", bypass="byp")
    s = arr.summarize(_solve(net)["__element_diag__"])
    assert s["m_dot_bypass"] < 0.0
    assert s["bypass_reversed"] == 1.0
    net, arr = _with_bypass(0.004)
    assert arr.summarize(_solve(net)["__element_diag__"])["bypass_reversed"] == 0.0


def _set_regime(net: FlowNetwork, regime: str) -> FlowNetwork:
    for item in (*net.elements.values(), *net.nodes.values()):
        if hasattr(item, "regime"):
            item.regime = regime
    return net


def _pressure_fed(p_feed: float, wall: bool = False) -> FlowNetwork:
    net = FlowNetwork()
    net.add_node(_pb("supply", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    net.add_node(_pb("bsrc", p_feed))
    net.add_node(PlenumNode("byp"))
    net.add_element(ChannelElement("feed", "bsrc", "byp", length=0.1, diameter=0.02))
    arr = _array(n_rows=5)
    arr.add_to(net, "supply", "exit", bypass="byp")
    if wall:
        from combaero.network import ConvectiveSurface, WallLayer

        net.add_node(_pb("hot_in", 1.02e5, 1200.0))
        net.add_node(_pb("hot_out", 1.0e5, 1200.0))
        L = 5 * 5.0 * D
        net.add_element(
            ChannelElement(
                "hot",
                "hot_in",
                "hot_out",
                length=L,
                diameter=0.05,
                surface=ConvectiveSurface(area=L * SPAN),
            )
        )
        arr.add_wall(net, "w", "hot", [WallLayer(thickness=0.001, conductivity=20.0)])
    return net


def _series() -> FlowNetwork:
    net = FlowNetwork()
    net.add_node(_pb("supply", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    net.add_node(PlenumNode("joint"))
    _array("a", n_rows=4).add_to(net, "supply", "joint")
    _array("b", n_rows=4).add_to(net, "supply", "exit", bypass="joint")
    return net


def _jacobian_error(
    net: FlowNetwork,
) -> tuple[float, int, NetworkSolver, np.ndarray, np.ndarray, np.ndarray]:
    """Worst relative error of the analytic Jacobian against central
    differences at the solved point, over each row's significant entries,
    and how many significant entries the analytic one leaves at zero."""
    s = NetworkSolver(net)
    res = s.solve(init_strategy="outletref_warmstart")
    assert res["__success__"], res.get("__message__")
    x = np.array(res["__x_solution__"], dtype=float)
    # Flows 1e-5, pressures and temperatures 1e-7: the compressible channel's
    # Fanno iteration repeats to ~1e-12 (see test_impingement_plate_element).
    rel = np.array([1e-5 if n.endswith(".m_dot") else 1e-7 for n in s.unknown_names])
    fd = approx_derivative(
        lambda v: s._residuals(v), x, method="3-point", abs_step=np.maximum(np.abs(x) * rel, 1e-10)
    )
    jac = s._residuals_and_jacobian(x)[1].toarray()
    sig = np.abs(fd) > 1e-3 * np.abs(fd).max(axis=1, keepdims=True)
    err = np.abs(jac - fd)[sig] / np.abs(fd[sig])
    return float(err.max()), int(((jac == 0.0) & sig).sum()), s, x, jac, fd


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.parametrize(
    "case",
    [
        "chamber",  # entry from a node's total pressure, bypass flow fixed
        "pressure",  # entry, the bypass flow an unknown
        "series",  # the continuation station
        "reversed",  # the gap draining back out: the entry term off
    ],
)
def test_the_jacobian_is_analytic_and_complete(case: str, regime: str) -> None:
    """Measured 2026-10-22: incompressible <= 1.6e-9, compressible <= 2.7e-5
    (the Fanno iteration's repeatability), no significant entry missing."""
    net = {
        "chamber": lambda: _with_bypass(0.004)[0],
        "pressure": lambda: _pressure_fed(1.032e5),
        "series": _series,
        "reversed": lambda: _pressure_fed(1.0005e5),
    }[case]()
    worst, missing, *_ = _jacobian_error(_set_regime(net, regime))
    assert missing == 0
    assert worst < (1e-6 if regime == "incompressible" else 1e-4), worst


def test_the_heat_transfer_reads_the_bypass_in_the_jacobian() -> None:
    """With a wall, row 1's h depends on the bypass flow (its Gc/Gj): the
    lead-in's flow column must carry dh/dm_c through the wall relay, exactly
    where nothing else contributes to the entry."""
    _, _, s, x, jac, fd = _jacobian_error(_pressure_fed(1.032e5, wall=True))
    j = s.unknown_names.index("ia__x0.m_dot")
    s._relay_flow_inputs = lambda *a, **k: None
    without = s._residuals_and_jacobian(x)[1].toarray()
    moved = np.where(np.abs(jac[:, j] - without[:, j]) > 1e-9 * np.abs(fd[:, j]).max())[0]
    assert moved.size > 0, "the bypass flow never reaches the heat transfer"
    for i in moved:
        assert abs(jac[i, j] - fd[i, j]) < abs(without[i, j] - fd[i, j]), s.unknown_names[i]
        assert abs(jac[i, j] - fd[i, j]) <= 0.02 * abs(fd[i, j]), (i, jac[i, j], fd[i, j])
