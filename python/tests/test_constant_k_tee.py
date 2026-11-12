"""ConstantKTeeElement: the "simplest model" junction tier.

Fixed handbook loss coefficients referenced to the common-leg dynamic
head (Idelchik convention) instead of the Mynard closure:

    Pt_i - Pt_jct + sign * K_i * q_dyn_com = 0     (non-common ports)

Anchored on the analytical closure itself (calibration methodology:
anchor on analytical, never fit measurements directly): at convergence
the port total-pressure drops must reproduce K_i * q_com exactly.
"""

import math

import numpy as np
import pytest
from _chamber_static import chamber_isentropic_rise
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
)
from combaero.network.mpce_element import ConstantKTeeElement

_D = 0.05
_A = math.pi * (_D / 2.0) ** 2
_DRY_AIR_Y = list(cb.mole_to_mass(cb.species.dry_air()))


def _branch_net(K_straight: float, K_branch: float, m_in: float = 0.3) -> FlowNetwork:
    """MFB-driven separating tee, free split between two equal PB sinks.

    MFB drive pins the common-arm direction (PB-driven variants of this
    fixture kept converging onto reversed-arm roots: with all losses
    referenced to q_com the split is hypersensitive to the sink spread).
    The long arm channels add resistance so the K asymmetry produces an
    interior forward split: with (K_b - K_s) * q_com = 2 * (q_str -
    q_bra) the analytical split at K = (0.5, 1.5), m_in = 0.3 is
    m_str = 0.225 / m_bra = 0.075.
    """
    Y = _DRY_AIR_Y
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("mfb_in", m_dot=m_in, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("pb_str", Pt=1.9e5, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("pb_bra", Pt=1.9e5, Tt=300.0, Y=Y))
    for nid in ("mc_com", "mc_str", "mc_bra"):
        net.add_node(MomentumChamberNode(nid, area=_A, regime="incompressible"))
    net.add_element(
        ChannelElement(
            "ch_in", "mfb_in", "mc_com", length=0.5, diameter=_D, regime="incompressible"
        )
    )
    net.add_element(
        ChannelElement(
            "ch_str", "mc_str", "pb_str", length=5.0, diameter=_D, regime="incompressible"
        )
    )
    net.add_element(
        ChannelElement(
            "ch_bra", "mc_bra", "pb_bra", length=5.0, diameter=_D, regime="incompressible"
        )
    )
    net.add_element(
        ConstantKTeeElement(
            id="jct",
            inlet_nodes=["mc_com"],
            outlet_nodes=["mc_str", "mc_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0],
            flow_direction="branch",
            K_ports={1: K_straight, 2: K_branch},
        )
    )
    return net


def _merge_net(K_straight: float, K_branch: float) -> FlowNetwork:
    """MFB-driven joining tee: two imposed suppliers into one PB sink."""
    Y = _DRY_AIR_Y
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("mfb_str", m_dot=0.18, Tt=300.0, Y=Y))
    net.add_node(MassFlowBoundary("mfb_bra", m_dot=0.12, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("pb_out", Pt=1.9e5, Tt=300.0, Y=Y))
    for nid in ("mc_com", "mc_str", "mc_bra"):
        net.add_node(MomentumChamberNode(nid, area=_A, regime="incompressible"))
    net.add_element(
        ChannelElement(
            "ch_str", "mfb_str", "mc_str", length=0.5, diameter=_D, regime="incompressible"
        )
    )
    net.add_element(
        ChannelElement(
            "ch_bra", "mfb_bra", "mc_bra", length=0.5, diameter=_D, regime="incompressible"
        )
    )
    net.add_element(
        ChannelElement(
            "ch_out", "mc_com", "pb_out", length=0.5, diameter=_D, regime="incompressible"
        )
    )
    net.add_element(
        ConstantKTeeElement(
            id="jct",
            inlet_nodes=["mc_str", "mc_bra"],
            outlet_nodes=["mc_com"],
            inlet_angles_deg=[0.0, 90.0],
            outlet_angles_deg=[0.0],
            flow_direction="merge",
            K_ports={0: K_straight, 1: K_branch},
        )
    )
    return net


def _q_com(sol: dict, m_key: str, p_key: str) -> float:
    """The common port's reference head, in the ports' regime.

    ConstantKTeeElement references K to whatever its port MomentumChamberNodes
    call dynamic head. These networks are incompressible throughout, and the
    ports follow (#588): q = m^2 / (2 rho A^2) at rho(P, T0). The compressible
    ports close on the isentropic stagnation rise instead (#357, #576), which
    this reduces to as M -> 0 -- kept below as the cross-check.
    """
    X_air = list(cb.mass_to_mole(_DRY_AIR_Y))
    T0 = 300.0
    P = sol[p_key]
    m = sol[m_key]
    rho = float(cb.density(T0, P, X_air))
    q = m * m / (2.0 * rho * _A * _A)
    q_isentropic = chamber_isentropic_rise(P, T0, m, _A, X_air)
    assert q_isentropic > q, "the compressible rise must exceed its low-Mach limit at M > 0"
    assert q_isentropic == pytest.approx(q, rel=0.02)
    return q


def test_branch_pt_drops_match_K_analytically():
    """The converged port Pt drops ARE the closure: Pt_com - Pt_arm =
    K_arm * q_com (common port has Pt continuity to P_jct... via
    Pt_com = Pt_jct)."""
    K_s, K_b = 0.5, 1.5
    net = _branch_net(K_s, K_b)
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=30.0)
    assert sol["__success__"], sol.get("__message__")
    assert sol["ch_in.m_dot"] > 0.0

    q = _q_com(sol, "ch_in.m_dot", "mc_com.P")
    drop_str = sol["mc_com.Pt"] - sol["mc_str.Pt"]
    drop_bra = sol["mc_com.Pt"] - sol["mc_bra.Pt"]
    assert drop_str == pytest.approx(K_s * q, rel=1e-6)
    assert drop_bra == pytest.approx(K_b * q, rel=1e-6)


def test_branch_split_responds_to_K_asymmetry():
    """Raising K_branch must reduce the branch share of the flow.

    Equal K gives an exactly even split between the equal sinks;
    K_b = 1.4 shifts it while keeping the branch arm forward
    ((K_b - K_s) * q_com stays below the arm-channel resistance scale).
    """
    lo = NetworkSolver(_branch_net(0.5, 0.5)).solve(timeout=30.0)
    hi = NetworkSolver(_branch_net(0.5, 1.4)).solve(timeout=30.0)
    assert lo["__success__"] and hi["__success__"]
    share_lo = lo["ch_bra.m_dot"] / lo["ch_in.m_dot"]
    share_hi = hi["ch_bra.m_dot"] / hi["ch_in.m_dot"]
    assert share_lo == pytest.approx(0.5, abs=0.02)
    assert 0.0 < share_hi < share_lo


def test_merge_converges_and_collector_below_suppliers():
    """Joining mode: Pt_arm = Pt_jct + K * q_com => the collector's Pt
    sits BELOW both supplier port Pts (dissipative junction)."""
    net = _merge_net(0.4, 1.0)
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=30.0)
    assert sol["__success__"], sol.get("__message__")
    assert sol["ch_str.m_dot"] > 0.0
    assert sol["ch_bra.m_dot"] > 0.0
    assert sol["ch_out.m_dot"] > 0.0
    assert sol["mc_com.Pt"] < sol["mc_str.Pt"]
    assert sol["mc_com.Pt"] < sol["mc_bra.Pt"]

    q = _q_com(sol, "ch_out.m_dot", "mc_com.P")
    assert sol["mc_str.Pt"] - sol["mc_com.Pt"] == pytest.approx(0.4 * q, rel=1e-6)
    assert sol["mc_bra.Pt"] - sol["mc_com.Pt"] == pytest.approx(1.0 * q, rel=1e-6)


@pytest.mark.parametrize("make_net", [_branch_net, _merge_net])
def test_jacobian_matches_fd(make_net):
    """Analytical junction rows vs central differences at x0."""
    net = make_net(0.5, 1.2)
    solver = NetworkSolver(net)
    net.resolve_all_topology()
    net.validate()
    x0 = np.asarray(solver._build_x0(), dtype=float)

    _, jac_analytical = solver._residuals_and_jacobian(x0)
    jac_analytical = jac_analytical.toarray()

    jac_fd = fd_jacobian(solver, x0, rel=1e-7, floor=1e-10)
    scale = np.maximum(np.abs(jac_fd), np.abs(jac_analytical)) + 1.0
    max_rel = float(np.max(np.abs(jac_analytical - jac_fd) / scale))
    assert max_rel < 5e-4, f"jacobian mismatch: {max_rel:.3e}"


def test_direction_guard_inherited_from_v2():
    """The MultiPortChamberElement direction check applies: a wrong-direction solution
    dict is rejected (memory rule: constant-K merges always pair with
    the post-solve guard)."""
    net = _merge_net(0.4, 1.0)
    net.resolve_all_topology()
    jct = net.elements["jct"]
    good = {"ch_str.m_dot": 0.3, "ch_bra.m_dot": 0.2, "ch_out.m_dot": 0.5}
    bad = {"ch_str.m_dot": 0.3, "ch_bra.m_dot": -0.2, "ch_out.m_dot": 0.1}
    assert jct.verify_solution_consistent(good) is True
    assert jct.verify_solution_consistent(bad) is False


def test_diagnostics_report_fixed_K():
    net = _branch_net(0.5, 1.5)
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=30.0)
    assert sol["__success__"]
    net.resolve_all_topology()
    jct = net.elements["jct"]
    states = [solver._get_node_state(net.nodes[n], solver._last_x) for n in jct.port_nodes]
    port_mdots = [
        float(jct._port_signs[i]) * sol[f"{jct._port_element_ids[i]}.m_dot"] for i in range(jct.N)
    ]
    diag = jct.diagnostics(states, sol["jct.Pt_jct"], port_mdots)
    assert diag["K_straight"] == 0.5
    assert diag["K_branch"] == 1.5
    assert "mass_flow_ratio" in diag
