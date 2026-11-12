"""IdealMomentumTeeElement: the ideal momentum junction kernel in a network (#517).

The kernel's closed forms are held in test_ideal_junction.py/.cpp. What is new
here is the wiring: the solver's unknowns reach the kernel in the right slots,
the angle convention maps the tee's port angles onto outward axes, ``P_jct`` is
a static pressure, and the solve converges -- including through a reversed
port, which the closure models and the Mynard tee refuses.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    IdealMomentumTeeElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
)
from validation.junction import joining_variants as jv

_D = 0.05
_A = math.pi * (_D / 2.0) ** 2
_Y = list(cb.mole_to_mass(cb.species.dry_air()))


def _branch_net(
    m_in: float = 0.05,
    theta: float = 90.0,
    regime: str = "incompressible",
    pt_bra: float = 1.9e5,
    tt_bra: float = 300.0,
    loss_factor: float | list[float] = 1.0,
) -> FlowNetwork:
    """MFB-fed dividing tee into two PB sinks through long arm channels."""
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("src", m_dot=m_in, Tt=300.0, Y=_Y))
    net.add_node(PressureBoundary("pb_str", Pt=1.9e5, Tt=300.0, Y=_Y))
    net.add_node(PressureBoundary("pb_bra", Pt=pt_bra, Tt=tt_bra, Y=_Y))
    for nid in ("mc_com", "mc_str", "mc_bra"):
        net.add_node(MomentumChamberNode(nid, area=_A))
    net.add_element(
        ChannelElement("ch_in", "src", "mc_com", length=0.5, diameter=_D, regime=regime)
    )
    net.add_element(
        ChannelElement("ch_str", "mc_str", "pb_str", length=5.0, diameter=_D, regime=regime)
    )
    net.add_element(
        ChannelElement("ch_bra", "mc_bra", "pb_bra", length=5.0, diameter=_D, regime=regime)
    )
    net.add_element(
        IdealMomentumTeeElement(
            id="jct",
            inlet_nodes=["mc_com"],
            outlet_nodes=["mc_str", "mc_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, theta],
            loss_factor=loss_factor,
        )
    )
    return net


def _merge_net(
    m: tuple[float, float] = (0.03, 0.02),
    theta: float = 90.0,
    regime: str = "incompressible",
    loss_factor: float | list[float] = 1.0,
) -> FlowNetwork:
    """Two MFB suppliers joining into one PB sink."""
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("mfb_str", m_dot=m[0], Tt=300.0, Y=_Y))
    net.add_node(MassFlowBoundary("mfb_bra", m_dot=m[1], Tt=300.0, Y=_Y))
    net.add_node(PressureBoundary("pb_out", Pt=1.9e5, Tt=300.0, Y=_Y))
    for nid in ("mc_com", "mc_str", "mc_bra"):
        net.add_node(MomentumChamberNode(nid, area=_A))
    net.add_element(
        ChannelElement("ch_str", "mfb_str", "mc_str", length=0.5, diameter=_D, regime=regime)
    )
    net.add_element(
        ChannelElement("ch_bra", "mfb_bra", "mc_bra", length=0.5, diameter=_D, regime=regime)
    )
    net.add_element(
        ChannelElement("ch_out", "mc_com", "pb_out", length=0.5, diameter=_D, regime=regime)
    )
    net.add_element(
        IdealMomentumTeeElement(
            id="jct",
            inlet_nodes=["mc_str", "mc_bra"],
            outlet_nodes=["mc_com"],
            inlet_angles_deg=[0.0, theta],
            outlet_angles_deg=[0.0],
            loss_factor=loss_factor,
        )
    )
    return net


def _solve(net: FlowNetwork) -> tuple[dict, NetworkSolver]:
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=60.0)
    assert sol["__success__"], sol.get("__message__")
    return sol, solver


def _rho(sol: dict, node: str) -> float:
    return float(sol[f"{node}.rho"])


def _receiver_identity(sol: dict, node: str, m_out: float, V_n: float) -> float:
    """P_j - p_jn + rho_j u_j (u_j - V.n_j), in Pa: the receiver row."""
    u = m_out / (_rho(sol, node) * _A)
    return sol[f"{node}.P"] - sol["jct.P_jct"] + _rho(sol, node) * u * (u - V_n)


# --- the rows hold on the converged network ---------------------------------


@pytest.mark.parametrize("theta", [45.0, 90.0, 135.0])
def test_dividing_rows_hold_at_the_converged_state(theta: float) -> None:
    """Supplier at p_jn; each receiver's momentum balance with V = v_supplier.

    Exact up to the role-smoothing bias, which is tiny next to the head.
    """
    sol, _ = _solve(_branch_net(theta=theta))
    m_c, m_s, m_b = sol["ch_in.m_dot"], sol["ch_str.m_dot"], sol["ch_bra.m_dot"]
    assert m_c > 0.0 and m_s > 0.0 and m_b > 0.0
    head = sol["mc_com.Pt"] - sol["mc_com.P"]
    assert sol["mc_com.P"] == pytest.approx(sol["jct.P_jct"], abs=1e-4 * head)
    # Supplier velocity along the main axis: the straight receiver sees all of
    # it, the branch receiver its projection on the branch axis.
    u_c = m_c / (_rho(sol, "mc_com") * _A)
    assert _receiver_identity(sol, "mc_str", m_s, u_c) == pytest.approx(0.0, abs=1e-4 * head)
    V_n = u_c * math.cos(math.radians(theta))
    assert _receiver_identity(sol, "mc_bra", m_b, V_n) == pytest.approx(0.0, abs=1e-4 * head)


def test_dividing_K_is_vector_borda_at_low_mach() -> None:
    """K_straight = (1 - q_s)^2, K_branch = 1 + q_b^2 - 2 q_b cos(theta)
    (equal areas), against the port nodes' own Pt. M ~ 0.03, so the
    isentropic port heads differ from rho u^2/2 by ~2e-4: the tolerance."""
    theta = 60.0
    sol, _ = _solve(_branch_net(theta=theta))
    q_c = sol["mc_com.Pt"] - sol["mc_com.P"]
    q_s = sol["ch_str.m_dot"] / sol["ch_in.m_dot"]
    q_b = sol["ch_bra.m_dot"] / sol["ch_in.m_dot"]
    K_s = (sol["mc_com.Pt"] - sol["mc_str.Pt"]) / q_c
    K_b = (sol["mc_com.Pt"] - sol["mc_bra.Pt"]) / q_c
    assert K_s == pytest.approx((1.0 - q_s) ** 2, abs=2e-3)
    assert K_b == pytest.approx(1.0 + q_b**2 - 2.0 * q_b * math.cos(math.radians(theta)), abs=2e-3)


@pytest.mark.parametrize("theta", [45.0, 90.0, 120.0])
def test_merge_K_is_blaisdell_manson_theory_at_low_mach(theta: float) -> None:
    """The converged network's K11/K12 are the scored ``BM theory`` variant."""
    sol, _ = _solve(_merge_net(theta=theta))
    q_c = sol["mc_com.Pt"] - sol["mc_com.P"]
    q_lat = sol["ch_bra.m_dot"] / sol["ch_out.m_dot"]
    U, A, th = jv.tee(q_lat, 1.0, theta)
    ref = jv.MODELS["BM theory"](U, A, th)
    assert (sol["mc_str.Pt"] - sol["mc_com.Pt"]) / q_c == pytest.approx(ref[1], abs=2e-3)
    assert (sol["mc_bra.Pt"] - sol["mc_com.Pt"]) / q_c == pytest.approx(ref[2], abs=2e-3)


# --- solve behaviour ----------------------------------------------------------


def _reversed_branch_solution() -> dict:
    """Declared dividing, but the branch sink sits ABOVE the junction, so the
    branch flows back in and the tee runs as a merge."""
    sol = NetworkSolver(_branch_net(pt_bra=1.95e5)).solve(timeout=60.0)
    assert sol["ch_bra.m_dot"] < 0.0  # back into the junction
    return sol


def test_a_reversed_port_is_a_state_not_an_error() -> None:
    """The branch flows back in and the declared dividing tee runs as a merge.
    The Mynard tee raises here (strict) or barriers it away (soft); this
    closure converges to the merge: both suppliers at p_jn, the one receiver
    on its momentum balance with V the mass-weighted supplied velocity."""
    sol = _reversed_branch_solution()
    assert sol["__success__"], sol.get("__message__")
    assert sol["__consistent__"] is True
    m_c, m_b, m_s = sol["ch_in.m_dot"], -sol["ch_bra.m_dot"], sol["ch_str.m_dot"]
    assert m_s == pytest.approx(m_c + m_b, rel=1e-9)
    head = sol["mc_str.Pt"] - sol["mc_str.P"]
    assert sol["mc_com.P"] == pytest.approx(sol["jct.P_jct"], abs=1e-4 * head)
    assert sol["mc_bra.P"] == pytest.approx(sol["jct.P_jct"], abs=1e-4 * head)
    # Common inflow arrives along the main axis, the branch inflow normal to
    # it, so only the common stream carries momentum along the straight run.
    u_c = m_c / (_rho(sol, "mc_com") * _A)
    V_x = m_c * u_c / (m_c + m_b)
    assert _receiver_identity(sol, "mc_str", m_s, V_x) == pytest.approx(0.0, abs=1e-4 * head)


def test_a_reversed_port_node_keeps_its_dynamic_head() -> None:
    """The solver feeds a junction port by the sign of its flow (#517). It
    used to follow the declared direction: the reversed port node saw the
    junction's outflow as a negative inflow cancelling the real one, and its
    Pt fell to P (0 Pa against a 711 Pa head)."""
    sol = _reversed_branch_solution()
    m = sol["ch_bra.m_dot"]
    rho = sol["mc_bra.rho"]
    head = 0.5 * rho * (m / (rho * _A)) ** 2
    # Isentropic port head vs rho u^2/2 at M ~ 0.07: ~1e-3 apart.
    assert sol["mc_bra.Pt"] - sol["mc_bra.P"] == pytest.approx(head, rel=5e-3)


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_network_jacobian_matches_central_differences_at_a_reversed_state(regime: str) -> None:
    """The reversed-port streams carry their own d(mass)/dx; checked at the
    reversed root perturbed by 1e-3, so every port keeps its sign. The
    reversed branch draws from a HOT sink, so the receiver's mixed T -- and
    with it its density -- depends on each stream's share."""
    net = _branch_net(pt_bra=1.95e5, regime=regime, tt_bra=400.0)
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=60.0, use_jac=False)
    assert sol["__success__"] and sol["ch_bra.m_dot"] < 0.0
    x = np.array([sol[n] for n in solver.unknown_names], dtype=float)
    x = x * (1.0 + 1e-3 * np.random.default_rng(1).standard_normal(x.size))
    _, J = solver._residuals_and_jacobian(x)
    J = J.toarray()
    J_fd = fd_jacobian(solver, x, rel=1e-7, floor=1e-10)
    scaled_err = np.abs(J - J_fd) * np.abs(x)
    row_ref = np.max(np.abs(J_fd) * np.abs(x), axis=1) + 1e-30
    assert float(np.max(scaled_err.max(axis=1) / row_ref)) < 1e-5


@pytest.mark.parametrize("builder", [_branch_net, _merge_net])
@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_network_jacobian_matches_central_differences(builder, regime: str) -> None:
    """At a converged state perturbed by 1e-3, NOT at x0.

    x0 seeds every pressure equal, which sits a compressible channel exactly
    on its direction switch, where its own derivative is one-sided; that is
    the channel's, not this element's, and a check there measures the kink.
    Compared as J_ij |x_j| against the row's largest such term.
    """
    net = builder(regime=regime)
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=60.0, use_jac=False)
    assert sol["__success__"], sol.get("__message__")
    x = np.array([sol[n] for n in solver.unknown_names], dtype=float)
    x = x * (1.0 + 1e-3 * np.random.default_rng(0).standard_normal(x.size))
    _, J = solver._residuals_and_jacobian(x)
    J = J.toarray()
    J_fd = fd_jacobian(solver, x, rel=1e-7, floor=1e-10)
    scaled_err = np.abs(J - J_fd) * np.abs(x)
    row_ref = np.max(np.abs(J_fd) * np.abs(x), axis=1) + 1e-30
    assert float(np.max(scaled_err.max(axis=1) / row_ref)) < 1e-5


def test_compressible_merge_converges_with_the_analytic_jacobian_from_a_warm_start() -> None:
    """From the FD solve's root perturbed by 1%, the analytic Jacobian walks
    back: the element's (f, J) is good away from the x0 tie."""
    net = _merge_net(m=(0.18, 0.12), theta=45.0, regime="compressible")
    solver = NetworkSolver(net)
    ref = solver.solve(timeout=60.0, use_jac=False)
    assert ref["__success__"]
    x0 = np.array([ref[n] for n in solver.unknown_names], dtype=float) * 1.01
    sol = NetworkSolver(_merge_net(m=(0.18, 0.12), theta=45.0, regime="compressible")).solve(
        timeout=60.0, x0=x0
    )
    assert sol["__success__"], sol.get("__message__")
    assert sol["jct.P_jct"] == pytest.approx(ref["jct.P_jct"], rel=1e-8)


def test_static_seed_moves_P_jct_below_the_common_Pt() -> None:
    """``analytical_pt_prop`` seeds P_jct from the common port's Pt; this
    closure's P_jct is static, so its seed sits one common head lower."""
    net = _branch_net(m_in=0.3)
    solver = NetworkSolver(net)
    net.resolve_all_topology()
    net.validate()
    ref = solver._infer_reference_state()
    overrides = solver._propagate_analytical_pt_prop(ref)
    pt = solver._propagate_pressure_guess(ref)["mc_com"]
    m = overrides["ch_in.m_dot"]
    rho = max(pt / (287.0 * float(ref.get("T", 300.0))), 1e-3)
    assert overrides["jct.P_jct"] == pytest.approx(pt - m * m / (2.0 * rho * _A * _A), rel=1e-12)
    assert overrides["jct.P_jct"] < pt


# --- element-level contracts --------------------------------------------------


def test_axes_follow_the_tee_angle_convention() -> None:
    """Outlet: the angle. Inlet: pi - angle (its flow arrives along the angle)."""
    el = IdealMomentumTeeElement(
        "j",
        inlet_nodes=["a", "b"],
        outlet_nodes=["c"],
        inlet_angles_deg=[0.0, 45.0],
        outlet_angles_deg=[30.0],
        port_areas=[1.0, 1.0, 1.0],
    )
    assert el.axis_rad() == pytest.approx(
        [math.pi, math.pi - math.radians(45.0), math.radians(30.0)]
    )


def test_more_ports_than_the_kernel_takes_is_refused() -> None:
    with pytest.raises(ValueError, match="at most 6"):
        IdealMomentumTeeElement("j", inlet_nodes=["a", "b", "c"], outlet_nodes=["d", "e", "f", "g"])


def test_four_port_cross_converges_with_both_side_arms_receiving() -> None:
    """N > 3 is in scope (the Mynard tee refuses it): a cross with one
    supply, a straight run and two opposed 90-deg side arms splits the side
    flow evenly by symmetry."""
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("src", m_dot=0.05, Tt=300.0, Y=_Y))
    for nid in ("pb_s", "pb_l", "pb_r"):
        net.add_node(PressureBoundary(nid, Pt=1.9e5, Tt=300.0, Y=_Y))
    for nid in ("mc_c", "mc_s", "mc_l", "mc_r"):
        net.add_node(MomentumChamberNode(nid, area=_A))
    net.add_element(ChannelElement("ch_in", "src", "mc_c", length=0.5, diameter=_D))
    for arm in ("s", "l", "r"):
        net.add_element(
            ChannelElement(f"ch_{arm}", f"mc_{arm}", f"pb_{arm}", length=5.0, diameter=_D)
        )
    net.add_element(
        IdealMomentumTeeElement(
            "jct",
            inlet_nodes=["mc_c"],
            outlet_nodes=["mc_s", "mc_l", "mc_r"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0, -90.0],
        )
    )
    sol, _ = _solve(net)
    assert sol["ch_l.m_dot"] == pytest.approx(sol["ch_r.m_dot"], rel=1e-8)
    assert sol["ch_s.m_dot"] > sol["ch_l.m_dot"] > 0.0


def test_energy_verifier_rejects_a_junction_that_creates_flow_work() -> None:
    net = _merge_net()
    net.resolve_all_topology()
    jct = net.elements["jct"]
    flows = {"ch_str.m_dot": 0.03, "ch_bra.m_dot": 0.02, "ch_out.m_dot": 0.05}
    good = {**flows, "mc_str.Pt": 1.9010e5, "mc_bra.Pt": 1.9012e5, "mc_com.Pt": 1.9005e5}
    bad = {**flows, "mc_str.Pt": 1.9010e5, "mc_bra.Pt": 1.9012e5, "mc_com.Pt": 1.9020e5}
    assert jct.verify_solution_consistent(good) is True
    assert jct.verify_solution_consistent(bad) is False
    # A reversed port is judged on energy alone: here the reversed branch
    # leaves at a Pt above the supplier's, paid for by the outlet.
    reversed_ = {**good, "ch_bra.m_dot": -0.01, "ch_str.m_dot": 0.06}
    assert jct.verify_solution_consistent(reversed_) is True


def test_diagnostics_name_the_regime_from_the_converged_flows() -> None:
    sol, solver = _solve(_branch_net(theta=90.0))
    net = solver.network
    jct = net.elements["jct"]
    states = [solver._get_node_state(net.nodes[n], solver._last_x) for n in jct.port_nodes]
    mdots = [jct._port_signs[i] * sol[f"{jct._port_element_ids[i]}.m_dot"] for i in range(3)]
    diag = jct.diagnostics(states, sol["jct.P_jct"], mdots)
    assert diag["P_jct"] == sol["jct.P_jct"]
    assert "p_jn" not in diag  # the duplicate went in #773
    assert diag["mass_flow_ratio"] == pytest.approx(sol["ch_bra.m_dot"] / sol["ch_in.m_dot"])
    q_b = diag["mass_flow_ratio"]
    assert diag["K_branch"] == pytest.approx(1.0 + q_b * q_b, abs=2e-3)
    # A merge coefficient on a split has no value: present, NaN.
    assert math.isnan(diag["K11"]) and math.isnan(diag["K12"])


def _cross_net(pt_left: float, tt_left: float = 300.0) -> FlowNetwork:
    """One supply, a straight run and two opposed 90-deg side arms."""
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("src", m_dot=0.05, Tt=300.0, Y=_Y))
    net.add_node(PressureBoundary("pb_s", Pt=1.9e5, Tt=300.0, Y=_Y))
    net.add_node(PressureBoundary("pb_l", Pt=pt_left, Tt=tt_left, Y=_Y))
    net.add_node(PressureBoundary("pb_r", Pt=1.9e5, Tt=300.0, Y=_Y))
    for nid in ("mc_c", "mc_s", "mc_l", "mc_r"):
        net.add_node(MomentumChamberNode(nid, area=_A))
    net.add_element(ChannelElement("ch_in", "src", "mc_c", length=0.5, diameter=_D))
    for arm in ("s", "l", "r"):
        net.add_element(
            ChannelElement(f"ch_{arm}", f"mc_{arm}", f"pb_{arm}", length=5.0, diameter=_D)
        )
    net.add_element(
        IdealMomentumTeeElement(
            "jct",
            inlet_nodes=["mc_c"],
            outlet_nodes=["mc_s", "mc_l", "mc_r"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0, -90.0],
        )
    )
    return net


def test_two_suppliers_into_two_receivers() -> None:
    """A reversed side arm makes two suppliers (common, left) feed two
    receivers (straight, right): each receiver carries the supplied mixture
    at its own flow, so every port node holds its own head, and the
    Jacobian of those split streams matches central differences. The left
    supply is hot, so the receivers' T depends on the split."""
    net = _cross_net(pt_left=1.93e5, tt_left=400.0)
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=60.0)
    assert sol["__success__"], sol.get("__message__")
    assert sol["ch_l.m_dot"] < 0.0 < sol["ch_r.m_dot"] and sol["ch_s.m_dot"] > 0.0
    for node, elem in (("mc_c", "ch_in"), ("mc_s", "ch_s"), ("mc_l", "ch_l"), ("mc_r", "ch_r")):
        m, rho = sol[f"{elem}.m_dot"], sol[f"{node}.rho"]
        head = 0.5 * rho * (m / (rho * _A)) ** 2
        assert sol[f"{node}.Pt"] - sol[f"{node}.P"] == pytest.approx(head, rel=5e-3), node
    # Both receivers carry the same supplied mixture.
    assert sol["mc_s.Tt"] == pytest.approx(sol["mc_r.Tt"], rel=1e-9)
    assert 300.0 < sol["mc_s.Tt"] < 400.0
    x = np.array([sol[n] for n in solver.unknown_names], dtype=float)
    x = x * (1.0 + 1e-3 * np.random.default_rng(2).standard_normal(x.size))
    _, J = solver._residuals_and_jacobian(x)
    J = J.toarray()
    J_fd = fd_jacobian(solver, x, rel=1e-7, floor=1e-10)
    scaled_err = np.abs(J - J_fd) * np.abs(x)
    row_ref = np.max(np.abs(J_fd) * np.abs(x), axis=1) + 1e-30
    assert float(np.max(scaled_err.max(axis=1) / row_ref)) < 1e-5


def test_jacobian_carries_port_temperature_and_composition() -> None:
    """A forward merge of air with HOT CO2: each port's density follows the
    T and Y the solver mixes into it, so the rows need d/d(port T) and
    d/d(port Y[k]). Without them the hot merge alone was off by 1.2e-4
    (scaled) in the flow columns."""
    net = _merge_net()
    co2 = list(cb.species.pure_species("CO2"))
    net.nodes["mfb_bra"].Tt = 450.0
    net.nodes["mfb_bra"].Y = list(cb.mole_to_mass(co2))
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=60.0, use_jac=False)
    assert sol["__success__"], sol.get("__message__")
    x = np.array([sol[n] for n in solver.unknown_names], dtype=float)
    x = x * (1.0 + 1e-3 * np.random.default_rng(3).standard_normal(x.size))
    _, J = solver._residuals_and_jacobian(x)
    J = J.toarray()
    J_fd = fd_jacobian(solver, x, rel=1e-7, floor=1e-10)
    scaled_err = np.abs(J - J_fd) * np.abs(x)
    row_ref = np.max(np.abs(J_fd) * np.abs(x), axis=1) + 1e-30
    rel = scaled_err.max(axis=1) / row_ref
    # The junction's N + 1 rows come last: element rows follow node rows,
    # and it is the last element added. Only those are judged. The port
    # node's own Pt closure (MomentumChamberNode) relays d/dT but not
    # d/dY[k], and is off by 7.6e-5 here with either junction -- a
    # separate, pre-existing gap.
    assert len(rel) == len(solver.unknown_names) == 14
    assert float(np.max(rel[-4:])) < 1e-5


# --- receiver loss factor (#526) -----------------------------------------------


def _fd_scaled_error(solver: NetworkSolver, sol: dict, seed: int) -> float:
    x = np.array([sol[n] for n in solver.unknown_names], dtype=float)
    x = x * (1.0 + 1e-3 * np.random.default_rng(seed).standard_normal(x.size))
    _, J = solver._residuals_and_jacobian(x)
    J = J.toarray()
    J_fd = fd_jacobian(solver, x, rel=1e-7, floor=1e-10)
    scaled_err = np.abs(J - J_fd) * np.abs(x)
    row_ref = np.max(np.abs(J_fd) * np.abs(x), axis=1) + 1e-30
    return float(np.max(scaled_err.max(axis=1) / row_ref))


@pytest.mark.parametrize("builder", [_branch_net, _merge_net])
def test_unit_loss_factor_is_the_closure_bit_for_bit(builder) -> None:
    """f = 1, scalar or per port, is the parameter-free closure: same
    residuals and Jacobian, exactly, at the same state."""
    ref = NetworkSolver(builder())
    sol = ref.solve(timeout=60.0)
    x = np.array([sol[n] for n in ref.unknown_names], dtype=float)
    r0, J0 = ref._residuals_and_jacobian(x)
    for f in (1.0, [1.0, 1.0, 1.0]):
        other = NetworkSolver(builder(loss_factor=f))
        other.solve(timeout=60.0)  # resolves topology (inherited port areas)
        assert other.unknown_names == ref.unknown_names
        r1, J1 = other._residuals_and_jacobian(x)
        assert np.array_equal(r0, r1)
        assert np.array_equal(J0.toarray(), J1.toarray())


@pytest.mark.parametrize("f", [0.05, 0.5, 1.5])
def test_dividing_loss_is_f_times_vector_borda(f: float) -> None:
    """One supplier, no mixing: each receiver loses f times its vector Borda
    head at the split the network converged to (low Mach, as the f = 1 test)."""
    theta = 60.0
    sol, _ = _solve(_branch_net(theta=theta, loss_factor=f))
    q_c = sol["mc_com.Pt"] - sol["mc_com.P"]
    q_s = sol["ch_str.m_dot"] / sol["ch_in.m_dot"]
    q_b = sol["ch_bra.m_dot"] / sol["ch_in.m_dot"]
    K_s = (sol["mc_com.Pt"] - sol["mc_str.Pt"]) / q_c
    K_b = (sol["mc_com.Pt"] - sol["mc_bra.Pt"]) / q_c
    assert K_s == pytest.approx(f * (1.0 - q_s) ** 2, abs=2e-3)
    borda = 1.0 + q_b**2 - 2.0 * q_b * math.cos(math.radians(theta))
    assert K_b == pytest.approx(f * borda, abs=2e-3)


@pytest.mark.parametrize("f", [0.0, 0.5, 1.5])
def test_merge_scales_only_the_collector_loss(f: float) -> None:
    """Two suppliers: f scales the collector's loss against the MIXED stream,
    rho/2 |v_c - V|^2, and leaves the suppliers' mixing alone, so the K pair
    moves from Blaisdell & Manson's by exactly (1 - f) L_c / q_c."""
    theta = 90.0
    sol, _ = _solve(_merge_net(theta=theta, loss_factor=f))
    q_c = sol["mc_com.Pt"] - sol["mc_com.P"]
    m_s, m_b, m_c = sol["ch_str.m_dot"], sol["ch_bra.m_dot"], sol["ch_out.m_dot"]
    q_lat = m_b / m_c
    U, A, th = jv.tee(q_lat, 1.0, theta)
    ref = jv.MODELS["BM theory"](U, A, th)
    # Mixed supplied velocity: the straight stream along +x, the lateral normal.
    u_s = m_s / (_rho(sol, "mc_str") * _A)
    u_b = m_b / (_rho(sol, "mc_bra") * _A)
    u_c = m_c / (_rho(sol, "mc_com") * _A)
    V_x, V_y = m_s * u_s / m_c, m_b * u_b / m_c
    L_c = 0.5 * _rho(sol, "mc_com") * ((u_c - V_x) ** 2 + V_y**2)
    shift = (1.0 - f) * L_c / q_c
    assert (sol["mc_str.Pt"] - sol["mc_com.Pt"]) / q_c == pytest.approx(ref[1] - shift, abs=2e-3)
    assert (sol["mc_bra.Pt"] - sol["mc_com.Pt"]) / q_c == pytest.approx(ref[2] - shift, abs=2e-3)
    assert sol["__consistent__"] is True  # dissipation >= 0 down to f = 0


@pytest.mark.parametrize("builder", [_branch_net, _merge_net])
@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_network_jacobian_with_unequal_loss_factors(builder, regime: str) -> None:
    net = builder(regime=regime, loss_factor=[0.3, 0.6, 1.4])
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=60.0, use_jac=False)
    assert sol["__success__"], sol.get("__message__")
    assert _fd_scaled_error(solver, sol, seed=2) < 1e-5


def test_network_jacobian_with_loss_factors_at_a_reversed_state() -> None:
    """The reversed branch supplies while its factor (0.3) is set: the factor
    must act only through the receiver weight, and the hot sink's density
    columns must still carry it."""
    net = _branch_net(
        pt_bra=1.95e5, regime="compressible", tt_bra=400.0, loss_factor=[0.5, 1.4, 0.3]
    )
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=60.0, use_jac=False)
    assert sol["__success__"] and sol["ch_bra.m_dot"] < 0.0
    assert _fd_scaled_error(solver, sol, seed=3) < 1e-5


@pytest.mark.parametrize("f", [-0.1, float("nan"), float("inf"), [1.0, 1.0]])
def test_invalid_loss_factor_is_refused(f) -> None:
    with pytest.raises(ValueError, match="loss_factor"):
        _branch_net(loss_factor=f)


def test_lossless_split_is_dissipative_only_at_order_mach_squared() -> None:
    """f = 0 removes a split's only loss as written at low Mach (rho u^2/2):
    it conserves the INCOMPRESSIBLE Pt. The port nodes' Pt is isentropic,
    larger by ~M^2/4 of each port's head, and the supplier's head is the
    largest, so scored on the ports' Pt the split loses a little: measured
    1.4e-4 of the supplied flow's head at port Mach 0.033 (M^2/4 = 2.7e-4).

    Before #576 the kernel read each port's density at the chamber's
    STAGNATION T, inflating the receivers' heads; the same split then created
    flow work (-3.2e-4) and the energy check refused it. With the port's
    static density on both sides it converges and passes.

    The ports are pinned compressible: the M^2/4 excess IS the subject."""
    net = _branch_net(theta=60.0, loss_factor=0.0)
    for node in net.nodes.values():
        if isinstance(node, MomentumChamberNode):
            node.regime = "compressible"
    sol = NetworkSolver(net).solve(timeout=60.0)
    assert sol["__success__"], sol.get("__message__")
    el = net.elements["jct"]
    out = sum(
        float(el._port_signs[i])
        * sol[f"{el._port_element_ids[i]}.m_dot"]
        * sol[f"{el.port_nodes[i]}.Pt"]
        for i in range(el.N)
    )
    head = sol["mc_com.Pt"] - sol["mc_com.P"]
    loss = -out / (sol["ch_in.m_dot"] * head)
    assert 0.0 < loss < 3e-4


def test_lossless_split_with_incompressible_ports_is_consistent() -> None:
    """The same f = 0 split on incompressible ports (#588). Its ports differ in
    density by ~1e-3, so mass-weighted Pt is not what the closure conserves:
    scored that way the split gained 1.1e-7 and the energy check refused an
    exact, dissipative root. Scored volumetrically, as the closure's guarantee
    is stated, it dissipates (1.5e-4 of the heads) and passes."""
    net = _branch_net(theta=60.0, loss_factor=0.0)
    for node in net.nodes.values():
        if isinstance(node, MomentumChamberNode):
            node.regime = "incompressible"
    sol = NetworkSolver(net).solve(timeout=60.0, auto_retry=False)
    assert sol["__converged__"], sol.get("__message__")
    assert sol["__consistent__"] is True, sol.get("__message__")
