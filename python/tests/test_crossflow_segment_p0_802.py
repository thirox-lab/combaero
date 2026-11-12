"""A crossflow segment's Fanno kernel runs from the duct's TRUE stagnation state (#802).

The station nodes of an impingement array or an effusion liner are plenums
holding the channel's STATIC pressure at the middle of each station. Before
#802 ``CrossflowSegmentElement`` fed that static pressure to the Fanno kernel
as the duct's stagnation pressure, a dynamic head too low, so the kernel's
choke came far too early: CR-3630 D.2 B(5,4,1) at mc/mj 1.01 ran its exit
segment at 99.7% of that false choke, 81% of the true one.

Each check below builds the duct's end states independently, HERE: the static
T by bisection on h(T0) = h(T) + (G R T / P)^2 / 2, the half-station drop by
hand at that density, and P0 from the static state by isentropic entropy
balance. The Fanno kernel itself is not under test; what the element hands it
is.
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
    FlowNetwork,
    ImpingementArray,
    MassFlowBoundary,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)
from combaero.network.components import ImpingementCrossflowElement, NetworkMixtureState

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_effusion_liner as tl  # noqa: E402

X_AIR = list(cb.species.dry_air())
Y_AIR = list(cb.mole_to_mass(X_AIR))
R_AIR = cb.specific_gas_constant(X_AIR)
T0 = 300.0
D_JET = 2.54e-3


def _static_T(Tt: float, P: float, G: float) -> float:
    """Static T of gas at stagnation Tt, static P, mass flux G (bisection)."""
    ht = cb.h_mass(Tt, X_AIR)
    lo, hi = 0.3 * Tt, Tt
    for _ in range(200):
        T = 0.5 * (lo + hi)
        if ht - cb.h_mass(T, X_AIR) - 0.5 * (G * R_AIR * T / P) ** 2 > 0.0:
            lo = T
        else:
            hi = T
    return 0.5 * (lo + hi)


def _P0(P: float, G: float) -> float:
    """Isentropic stagnation pressure of the static state (P, T(T0, P, G))."""
    T = _static_T(T0, P, G)
    return P * math.exp((cb.s_mass(T0, X_AIR, P) - cb.s_mass(T, X_AIR, P)) / R_AIR)


def _half_drop_normal(m_a: float, m_b: float, P: float, A: float) -> float:
    """Half a normal-injection station's static drop at its own static density."""
    rho = P / (R_AIR * _static_T(T0, P, 0.5 * (m_a + m_b) / A))
    return 0.5 * (m_b * m_b - m_a * m_a) / (rho * A * A)


def _fanno(seg: ImpingementCrossflowElement, P0: float, P_exit_static: float) -> float:
    """Mass flow the duct passes from P0 (at T0) to an exit STATIC pressure."""
    flow = cb.fanno_channel_flow(
        P0,
        T0,
        X_AIR,
        P_exit_static,
        False,
        seg.length,
        seg.diameter,
        seg.roughness,
        seg.friction_model,
    )
    return seg.area * flow.G


def _state(P: float, m: float) -> NetworkMixtureState:
    return NetworkMixtureState(P=P, Pt=P, Tt=T0, m_dot=m, Y=Y_AIR)


def _exit_pair() -> tuple[ImpingementCrossflowElement, ImpingementCrossflowElement]:
    """x1: c1 -> c2 and x2: c2 -> exit, normal-injection stations, resolved."""
    g = FlowNetwork()
    g.add_node(PlenumNode("c1"))
    g.add_node(PlenumNode("c2"))
    g.add_node(PressureBoundary("exit", Pt=1.0e5, Tt=T0, Y=Y_AIR))
    kw = {"length": 5 * D_JET, "height": D_JET, "span": 0.05}
    g.add_element(ImpingementCrossflowElement("x1", "c1", "c2", **kw))
    g.add_element(ImpingementCrossflowElement("x2", "c2", "exit", **kw))
    for e in g.elements.values():
        e.resolve_topology(g)
    return g.elements["x1"], g.elements["x2"]


def test_the_exit_segment_passes_what_its_true_stagnation_state_drives() -> None:
    """An exit segment after a long array: 90% of its flow arrives along the
    channel. Its duct inlet static is P_c2 less half the station, and its P0
    is that static state's, a dynamic head above P_c2. Set the exit static so
    that the duct passes m from that true P0 at 85% of its choked flow: the
    element's row must then vanish.

    Fed P_c2 as P0 (before #802), the kernel's choked flow is BELOW m: no
    exit pressure at all passes it, the false choke of #802."""
    _, seg = _exit_pair()
    A, P, frac = seg.area, 1.3e5, 0.9
    m = 0.02
    for _ in range(60):  # m at 85% of the true choked flow
        e_in = P - _half_drop_normal(frac * m, m, P, A)
        P0 = _P0(e_in, m / A)
        m = 0.5 * m + 0.5 * 0.85 * _fanno(seg, P0, 1.0)
    e_in = P - _half_drop_normal(frac * m, m, P, A)
    P0 = _P0(e_in, m / A)
    assert P0 / P - 1.0 > 0.15  # a real test: the dynamic head is large
    assert _fanno(seg, P, 1.0) < m  # fed P as P0, the duct chokes below m

    lo, hi = 1.0, e_in  # the exit static the true state needs
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _fanno(seg, P0, mid) > m:
            lo = mid
        else:
            hi = mid
    P_exit = 0.5 * (lo + hi)
    assert _fanno(seg, P0, P_exit) == pytest.approx(m, rel=1e-9)

    res, _ = seg.residuals(_state(P, m), _state(P_exit, m), {"x1": frac * m})
    assert abs(res[0]) < 1e-5 * m
    d = seg.diagnostics(_state(P, m), _state(P_exit, m), {"x1": frac * m})
    assert d["choked"] == 0.0


def test_a_mid_array_segment_runs_between_its_end_statics() -> None:
    """Station to station: inlet static P_c1 - half of c1, outlet static
    P_c2 + half of c2, each end's P0 its own static state's. At the outlet
    static the true duct state needs, the row vanishes."""
    x1, _ = _exit_pair()
    A, P = x1.area, 1.3e5
    m_a, m, m_b = 0.02, 0.03, 0.04  # arriving at c1, x1, leaving c2
    e_in = P - _half_drop_normal(m_a, m, P, A)
    P0 = _P0(e_in, m / A)
    lo, hi = 1.0, e_in
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _fanno(x1, P0, mid) > m:
            lo = mid
        else:
            hi = mid
    e_out = 0.5 * (lo + hi)
    # c2's own static: the outlet static less half of c2 (solved by fixed point).
    P2 = e_out
    for _ in range(50):
        P2 = e_out - _half_drop_normal(m, m_b, P2, A)
    assert e_in / P0 < 0.9  # M ~ 0.4 at the duct inlet

    # x1's own prev segment is None (c1 is a closed end in _exit_pair), so
    # give it one: the arriving flow m_a.
    g = FlowNetwork()
    for nid in ("c0", "c1", "c2"):
        g.add_node(PlenumNode(nid))
    g.add_node(PressureBoundary("exit", Pt=1.0e5, Tt=T0, Y=Y_AIR))
    kw = {"length": 5 * D_JET, "height": D_JET, "span": 0.05}
    for eid, a, b in (("x0", "c0", "c1"), ("x1", "c1", "c2"), ("x2", "c2", "exit")):
        g.add_element(ImpingementCrossflowElement(eid, a, b, **kw))
    for e in g.elements.values():
        e.resolve_topology(g)
    seg = g.elements["x1"]
    res, _ = seg.residuals(_state(P, m), _state(P2, m), {"x0": m_a, "x2": m_b})
    assert abs(res[0]) < 1e-5 * m


def _array(m_exit: float, n: int = 3, regime: str = "compressible") -> FlowNetwork:
    net = FlowNetwork()
    net.add_node(PressureBoundary("supply", Pt=2.0e5, Tt=T0, Y=Y_AIR))
    net.add_node(MassFlowBoundary("exit", m_dot=m_exit, Tt=T0, Y=Y_AIR))
    ImpingementArray(
        "a", n_rows=n, d_jet=D_JET, xn_d=5, yn_d=4, z_d=1, span=0.05, plate_thickness=D_JET
    ).add_to(net, "supply", "exit")
    if regime != "compressible":
        for e in net.elements.values():
            e.regime = regime
    return net


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r.get("__message__")
    return s, r


def test_low_mach_limit_is_the_incompressible_balance() -> None:
    """No dynamic head counted twice: at M ~ 0.03 the compressible array's
    jet split is the incompressible one to O(M^2)."""
    _, rc = _solve(_array(0.0025))
    _, ri = _solve(_array(0.0025, regime="incompressible"))
    for i in (1, 2, 3):
        assert rc[f"a__p{i}.m_dot"] == pytest.approx(ri[f"a__p{i}.m_dot"], rel=2e-3)


def _row_jac_err(s: NetworkSolver, x: np.ndarray, rows: list[int]) -> float:
    """RESIDUAL-ROW-SCALING: compare J_ij |x_j| on the given rows."""
    _, J = s._residuals_and_jacobian(x)
    # A 1e-6 flow step: the reservoir entry row is curved enough in its own
    # flow that 1e-5 leaves a 1e-5 truncation error (it falls as the step).
    rel = np.array([1e-6 if n.endswith("m_dot") else 1e-7 for n in s.unknown_names])
    fd = fd_jacobian(s, x, rel=rel, floor=1e-10)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = (J.toarray() * cols)[rows], (fd * cols)[rows]
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    return float(np.max(np.abs(Js - fds) / scale))


def _segment_rows(s: NetworkSolver, x: np.ndarray, marker: str) -> list[int]:
    return [i for i, (owner, _) in enumerate(s._row_kinds(x)) if marker in owner]


def test_an_array_passes_what_the_false_choke_refused() -> None:
    """Six rows into a fixed exit flow, the exit segment at duct M ~ 0.7.
    Fed static P as P0 (before #802), the kernel choked below that flow:
    the solve ran onto the choke's plateau, the free exit pressure went
    negative and the diagnostics raised out of solve() (from 0.041 kg/s).
    It also reported every exit segment choked (M_exit_duct 1.0) at flows
    where the true duct runs at M ~ 0.35."""
    _, r = _solve(_array(0.043, n=6))
    d = r["__element_diag__"]["a__x6"]
    assert d["choked"] == 0.0
    assert 0.65 < d["M_in_duct"] < d["M_exit_duct"] < 0.75
    _, r = _solve(_array(0.03, n=6))
    assert r["__element_diag__"]["a__x6"]["M_exit_duct"] == pytest.approx(0.35, abs=0.02)


def test_the_array_jacobian_matches_finite_differences() -> None:
    """Six rows, the exit segment at duct M ~ 0.58."""
    s, r = _solve(_array(0.04, n=6))
    x = np.array(r["__x_solution__"], dtype=float)
    assert r["__element_diag__"]["a__x6"]["M_in_duct"] > 0.55
    rows = _segment_rows(s, x, "__x")
    assert len(rows) == 6
    assert _row_jac_err(s, x, rows) < 1e-5


@pytest.mark.slow
def test_the_liner_segment_jacobian_matches_finite_differences() -> None:
    """Four bleed stations at duct M ~ 0.6: the reservoir entry, three
    station-to-station segments and the exit. Segment rows only: the panel
    rows carry #801's separate orifice-column error."""
    s, r = _solve(tl._net(tl._liner(4), p_in=1.5e5, p_out=1.1e5))
    x = np.array(r["__x_solution__"], dtype=float)
    rows = _segment_rows(s, x, "__s")
    assert len(rows) == 5
    assert _row_jac_err(s, x, rows) < 1e-5


def _chain3() -> FlowNetwork:
    """x0: c0 -> c1, x1: c1 -> c2, x2: c2 -> exit (a plenum), resolved."""
    g = FlowNetwork()
    for nid in ("c0", "c1", "c2", "exit"):
        g.add_node(PlenumNode(nid))
    g.add_node(PressureBoundary("ref", Pt=1.0e5, Tt=T0, Y=Y_AIR))
    kw = {"length": 5 * D_JET, "height": D_JET, "span": 0.05}
    for eid, a, b in (("x0", "c0", "c1"), ("x1", "c1", "c2"), ("x2", "c2", "exit")):
        g.add_element(ImpingementCrossflowElement(eid, a, b, **kw))
    for e in g.elements.values():
        e.resolve_topology(g)
    return g


@pytest.mark.parametrize(
    ("seg", "P_in", "P_out", "m", "flows"),
    [
        ("x1", 1.3e5, 1.2e5, 0.03, {"x0": 0.02, "x2": 0.04}),  # subsonic
        ("x1", 1.3e5, 0.4e5, 0.03, {"x0": 0.02, "x2": 0.04}),  # outlet below P*
        ("x1", 0.5e5, 0.4e5, 0.03, {"x0": 0.02, "x2": 0.04}),  # both ends below P*
        ("x1", 1.2e5, 1.3e5, -0.03, {"x0": -0.02, "x2": -0.04}),  # reversed
        ("x2", 1.3e5, 1.0e5, 0.04, {"x1": 0.03}),  # exit, static continuity
        ("x2", 1.3e5, 0.3e5, 0.04, {"x1": 0.03}),  # exit below P*
        ("x2", 1.0e5, 1.3e5, -0.04, {"x1": -0.03}),  # the exit plenum feeds
    ],
)
def test_the_duct_end_partials_match_finite_differences(
    seg: str, P_in: float, P_out: float, m: float, flows: dict[str, float]
) -> None:
    """Every partial of both duct-end Pt the row chains: node P and T0, own
    flow, neighbour flows, on each branch of the closure (subsonic, the C1
    continuation below the sonic static, reversal, a terminal end)."""
    g = _chain3()
    e = g.elements[seg]
    a_id, b_id = e.from_node, e.to_node

    def ends(dx: dict[str, float]) -> tuple[float, float]:
        sa = _state(P_in + dx.get(f"{a_id}.P", 0.0), m + dx.get(f"{seg}.m_dot", 0.0))
        sb = _state(P_out + dx.get(f"{b_id}.P", 0.0), m + dx.get(f"{seg}.m_dot", 0.0))
        for st, nid in ((sa, a_id), (sb, b_id)):
            st.Tt += dx.get(f"{nid}.Tt", 0.0)
            st.Pt = st.P
        f = {k: v + dx.get(f"{k}.m_dot", 0.0) for k, v in flows.items()}
        src, dst, _, _ = e._duct_ends(sa, sb, e._momentum_terms(sa, sb, f))
        return src.Pt, dst.Pt

    sa, sb = _state(P_in, m), _state(P_out, m)
    _, _, d_in, d_out = e._duct_ends(sa, sb, e._momentum_terms(sa, sb, flows))
    keys = [f"{a_id}.P", f"{b_id}.P", f"{a_id}.Tt", f"{b_id}.Tt", f"{seg}.m_dot"]
    keys += [f"{k}.m_dot" for k in flows]
    steps = {"P": 1.0, "Tt": 1e-3, "m_dot": 1e-7}
    for key in keys:
        h = steps[key.rsplit(".", 1)[1]]
        hi, lo = ends({key: h}), ends({key: -h})
        for i, d in enumerate((d_in, d_out)):
            fd = (hi[i] - lo[i]) / (2.0 * h)
            an = d.get(key, 0.0) + (
                d.get(key.replace(".P", ".Pt"), 0.0) if key.endswith(".P") else 0.0
            )
            assert an == pytest.approx(fd, rel=1e-5, abs=1e-6 * max(1.0, abs(fd))), (key, i)
