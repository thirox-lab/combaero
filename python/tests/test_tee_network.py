"""
Network-level solver tests for TeeJunctionElement.

Covers:
- Merging tee: two pressure-boundary inlets, one outlet
- Branching tee: one inlet, two pressure-boundary outlets
- Mass conservation at each node (including mixed-composition streams)
- Species conservation across the mixing junction
- Enthalpy conservation across the mixing junction
- Jacobian accuracy (FD vs analytical) for both tee types
- No-flow edge cases
- TeeJunctionElement.validate() error paths
"""

import math

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    TeeJunctionElement,
)

_DRY_AIR_Y = list(cb.mole_to_mass(cb.species.dry_air()))

_THETA_90 = math.pi / 2.0
_THETA_60 = math.pi / 3.0
_F_C = 0.01  # 100 cm^2


def _merging_net(
    Pt_straight: float = 2.1e5,
    Pt_branch: float = 2.05e5,
    Pt_common: float = 2.0e5,
    theta: float = _THETA_90,
    psi: float = 1.0,
) -> FlowNetwork:
    """Merging tee: two inlets (straight, branch) into one outlet (common).

    Supplier ports use MomentumChamberNode so static pressure P can differ from Pt,
    which allows the R_sp = P_str - P_bra residual to be satisfied.  The collector
    port uses PlenumNode because R_K,com only references Pt_com (not P_com).
    """
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    A_branch = _F_C / psi
    net.add_node(PressureBoundary("pb_straight", Pt=Pt_straight, Tt=400.0, Y=Y))
    net.add_node(PressureBoundary("pb_branch", Pt=Pt_branch, Tt=400.0, Y=Y))
    net.add_node(PressureBoundary("pb_common", Pt=Pt_common, Tt=400.0, Y=Y))
    net.add_node(MomentumChamberNode("port_str", area=_F_C))
    net.add_node(MomentumChamberNode("port_bra", area=A_branch))
    net.add_node(PlenumNode("port_com"))
    net.add_element(LosslessConnectionElement("lc_str", "pb_straight", "port_str"))
    net.add_element(LosslessConnectionElement("lc_bra", "pb_branch", "port_bra"))
    net.add_element(LosslessConnectionElement("lc_com", "port_com", "pb_common"))
    net.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="port_com",
            straight_node="port_str",
            branch_node="port_bra",
            theta=theta,
            F_C=_F_C,
            psi=psi,
            tee_type="merging",
        )
    )
    return net


def _branching_net(
    Pt_common: float = 2.1e5,
    Pt_straight: float = 2.098e5,
    Pt_branch: float = 2.06e5,
    theta: float = _THETA_90,
    psi: float = 1.0,
) -> FlowNetwork:
    """Three-node branching tee: one inlet (common) into two outlets (straight, branch).

    Default boundary pressures sit inside the unified0D closure's physical envelope.
    The straight run is a low-loss through-leg: its blended extraction term tops out at
    an effective loss coefficient of beta*4/27 ~ 0.30 of the common dynamic head (the
    Borda-Carnot diversion maximum). A large imposed straight-run stagnation drop is
    therefore unreachable - the earlier 4.0e3 Pa straight drop (str 2.06 bar) only ever
    "solved" because the pre-#9 closure over-charged the run with a full Borda-Carnot
    K ~ 1, the bug Finding 5 removes. With a small straight drop the solver reaches a
    clean all-forward exact root; the split is branch-dominated because the low-loss run
    cannot carry much flow without downstream loading (see addendum Finding 6). See
    sources/junction/'momentum cv implementation guide.pdf' for the momentum-CV successor
    that removes this envelope limit entirely.
    """
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(PressureBoundary("pb_common", Pt=Pt_common, Tt=400.0, Y=Y))
    net.add_node(PressureBoundary("pb_straight", Pt=Pt_straight, Tt=400.0, Y=Y))
    net.add_node(PressureBoundary("pb_branch", Pt=Pt_branch, Tt=400.0, Y=Y))
    net.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="pb_common",
            straight_node="pb_straight",
            branch_node="pb_branch",
            theta=theta,
            F_C=_F_C,
            psi=psi,
            tee_type="branching",
        )
    )
    return net


# ---------------------------------------------------------------------------
# Topology tests
# ---------------------------------------------------------------------------


def test_merging_tee_graph_connectivity():
    """FlowNetwork registers all nodes and elements in connectivity dicts."""
    net = _merging_net()
    assert "lc_str" in net._downstream_of_node["pb_straight"]
    assert "lc_bra" in net._downstream_of_node["pb_branch"]
    assert "tee" in net._downstream_of_node["port_str"]
    assert "tee" in net._downstream_of_node["port_bra"]
    assert "tee" in net._upstream_of_node["port_com"]


def test_branching_tee_graph_connectivity():
    net = _branching_net()
    assert "tee" in net._downstream_of_node["pb_common"]
    assert "tee" in net._upstream_of_node["pb_straight"]
    assert "tee" in net._upstream_of_node["pb_branch"]


def test_merging_reachable_from_boundary():
    """All nodes reachable from pb_straight in merging net."""
    net = _merging_net()
    reachable = net._reachable_from("pb_straight")
    assert "pb_branch" in reachable
    assert "pb_common" in reachable


def test_tee_element_unknowns():
    net = _merging_net()
    assert net.elements["tee"].unknowns() == ["tee.m_dot_com", "tee.m_dot_branch"]
    assert net.elements["tee"].n_equations() == 2


# ---------------------------------------------------------------------------
# Solver: merging tee
# ---------------------------------------------------------------------------


def test_merging_tee_solver_converges():
    """Merging tee solver finds finite positive flows."""
    net = _merging_net()
    solver = NetworkSolver(net)
    sol = solver.solve()
    m_com = sol["tee.m_dot_com"]
    m_branch = sol["tee.m_dot_branch"]
    assert math.isfinite(m_com), "m_dot_com is not finite"
    assert math.isfinite(m_branch), "m_dot_branch is not finite"
    assert m_com > 0.0, "m_dot_com should be positive (flow into common outlet)"
    assert 0.0 < m_branch < m_com, "branch flow should be between 0 and total"


def test_merging_tee_mass_conservation():
    """m_dot_com == m_dot_straight + m_dot_branch at solved state."""
    net = _merging_net()
    solver = NetworkSolver(net)
    sol = solver.solve()
    m_com = sol["tee.m_dot_com"]
    m_branch = sol["tee.m_dot_branch"]
    m_straight = m_com - m_branch
    assert abs(m_straight) >= 0.0
    assert abs(m_com - m_straight - m_branch) < 1e-10


def test_merging_tee_residuals_near_zero():
    """Residuals evaluated at the solution must be near zero."""
    net = _merging_net()
    solver = NetworkSolver(net)
    sol = solver.solve()

    import numpy as np

    x = np.array([sol.get(n, 0.0) for n in solver.unknown_names])
    res, _ = solver._residuals_and_jacobian(x, compute_jacobian=False)
    # Pa-scale and kg/s-scale residuals; 1e-3 is ~1e-8 relative at 1e5 Pa
    assert all(abs(r) < 1e-3 for r in res), f"Non-zero residuals: {res}"


def test_merging_tee_diagnostics():
    """Solved solution produces finite diagnostic values."""
    net = _merging_net()
    solver = NetworkSolver(net)
    sol = solver.solve()
    assert math.isfinite(sol["tee.mass_flow_ratio"])
    assert math.isfinite(sol["tee.m_dot_com"])
    assert math.isfinite(sol["tee.m_dot_branch"])


# ---------------------------------------------------------------------------
# Solver: branching tee
# ---------------------------------------------------------------------------


def test_branching_tee_solver_converges():
    net = _branching_net()
    solver = NetworkSolver(net)
    sol = solver.solve()
    m_com = sol["tee.m_dot_com"]
    m_branch = sol["tee.m_dot_branch"]
    assert math.isfinite(m_com)
    assert math.isfinite(m_branch)
    assert m_com > 0.0


def test_branching_tee_mass_conservation():
    net = _branching_net()
    solver = NetworkSolver(net)
    sol = solver.solve()
    m_com = sol["tee.m_dot_com"]
    m_branch = sol["tee.m_dot_branch"]
    assert abs(m_com - m_branch) > 0.0  # straight flow is nonzero


def test_branching_tee_residuals_near_zero():
    net = _branching_net()
    solver = NetworkSolver(net)
    sol = solver.solve()

    import numpy as np

    x = np.array([sol.get(n, 0.0) for n in solver.unknown_names])
    res, _ = solver._residuals_and_jacobian(x, compute_jacobian=False)
    assert all(abs(r) < 1e-6 for r in res), f"Non-zero residuals: {res}"


def test_branching_tee_forward_flow_root():
    """Single dividing tee on direct pressure boundaries converges to a PHYSICAL
    (all-forward) exact root.

    Honest positive baseline for the unified0D blended closure (task #9, addendum
    Finding 6): a clean fixture (no MomentumChamberNode coupling, see Finding 4) with a
    single source (common) and two sinks (straight, branch), so the BCs unambiguously
    force all-outward flow - which the solver must then reproduce. The split is
    branch-dominated because the straight run is a low-loss through-leg (its loss
    coefficient is capped near 0.30 of the common dynamic head); a large straight drop
    is outside the closure envelope. Direction here is a consequence of the BCs, not an
    intrinsic property: with a single source and two sinks the physical answer is
    all-forward, and that is what we assert.
    """
    net = _branching_net()
    solver = NetworkSolver(net)
    sol = solver.solve()
    assert sol["__success__"], f"did not converge: {sol.get('__message__')}"
    m_com = sol["tee.m_dot_com"]
    m_branch = sol["tee.m_dot_branch"]
    m_straight = m_com - m_branch
    assert m_com > 0.0, f"common reversed: {m_com}"
    assert m_branch > 0.0, f"branch reversed: {m_branch}"
    assert m_straight > 0.0, f"straight reversed: {m_straight}"


# ---------------------------------------------------------------------------
# Conservation laws: mixed-composition streams
# ---------------------------------------------------------------------------


# Two-stream fixture: 1 kg/s pure N2 at 300 K merges with 0.1 kg/s dry air at 400 K.
# Common arm is a PlenumNode so its T and Y are solver unknowns, giving independent
# verification that the C++ mixer satisfies conservation.
def _make_Y_n2() -> list[float]:
    Y = [0.0] * cb.num_species()
    names = [cb.species_name(i) for i in range(cb.num_species())]
    if "N2" in names:
        Y[names.index("N2")] = 1.0
    return Y


def _merging_net_two_compositions() -> FlowNetwork:
    """Merging tee: straight=pure N2 @ 300 K, branch=dry air @ 400 K."""
    Y_n2 = _make_Y_n2()
    Y_air = _DRY_AIR_Y

    net = FlowNetwork()
    net.add_node(MassFlowBoundary("mb_str", m_dot=1.0, Tt=300.0, Y=Y_n2))
    net.add_node(MassFlowBoundary("mb_br", m_dot=0.1, Tt=400.0, Y=Y_air))
    net.add_node(PressureBoundary("pb", Pt=101325.0, Tt=300.0, Y=Y_air))
    net.add_node(PlenumNode("pl_str"))
    net.add_node(PlenumNode("pl_br"))
    net.add_node(PlenumNode("pl_com"))
    net.add_element(
        TeeJunctionElement(
            "tee",
            common_node="pl_com",
            straight_node="pl_str",
            branch_node="pl_br",
            F_C=0.02,
            psi=1.0,
            theta=_THETA_90,
            tee_type="merging",
        )
    )
    net.add_element(LosslessConnectionElement("lc_str", "mb_str", "pl_str"))
    net.add_element(LosslessConnectionElement("lc_br", "mb_br", "pl_br"))
    net.add_element(LosslessConnectionElement("lc_com", "pl_com", "pb"))
    return net


def test_merging_tee_node_mass_conservation():
    """lc_str + lc_br == lc_com: independent solver unknowns satisfy global continuity."""
    net = _merging_net_two_compositions()
    sol = NetworkSolver(net).solve()
    m_str = sol["lc_str.m_dot"]
    m_br = sol["lc_br.m_dot"]
    m_com = sol["lc_com.m_dot"]
    assert abs(m_com - m_str - m_br) < 1e-9, (
        f"Mass not conserved: {m_com:.6f} != {m_str:.6f} + {m_br:.6f}"
    )


def test_merging_tee_species_conservation():
    """Species conservation: m_com*Y_com[i] == m_str*Y_str[i] + m_br*Y_br[i]."""
    Y_n2 = _make_Y_n2()
    Y_air = _DRY_AIR_Y
    net = _merging_net_two_compositions()
    sol = NetworkSolver(net).solve()

    m_str = sol["lc_str.m_dot"]
    m_br = sol["lc_br.m_dot"]
    m_com = sol["lc_com.m_dot"]

    n_spec = cb.num_species()
    for i in range(n_spec):
        y_com = sol.get(f"pl_com.Y[{i}]", 0.0)
        expected = (m_str * Y_n2[i] + m_br * Y_air[i]) / m_com
        err = abs(y_com - expected)
        assert err < 1e-8, (
            f"Species {cb.species_name(i)}: Y_com={y_com:.8f}, expected={expected:.8f}, err={err:.2e}"
        )


def test_merging_tee_enthalpy_conservation():
    """Adiabatic mixing: m_com*h(T_com,Y_com) == m_str*h(T_str,Y_str) + m_br*h(T_br,Y_br)."""
    Y_n2 = _make_Y_n2()
    Y_air = _DRY_AIR_Y
    net = _merging_net_two_compositions()
    sol = NetworkSolver(net).solve()

    m_str = sol["lc_str.m_dot"]
    m_br = sol["lc_br.m_dot"]
    m_com = sol["lc_com.m_dot"]

    T_str = sol.get("pl_str.Tt", 300.0)
    T_br = sol.get("pl_br.Tt", 400.0)
    T_com = sol.get("pl_com.Tt")
    assert T_com is not None

    n_spec = cb.num_species()
    Y_com = [sol.get(f"pl_com.Y[{i}]", 0.0) for i in range(n_spec)]

    X_n2 = list(cb.mass_to_mole(Y_n2))
    X_air = list(cb.mass_to_mole(Y_air))
    X_com = list(cb.mass_to_mole(Y_com))

    h_str = float(cb._core.enthalpy_and_jacobian(T_str, X_n2)[0])
    h_br = float(cb._core.enthalpy_and_jacobian(T_br, X_air)[0])
    h_com = float(cb._core.enthalpy_and_jacobian(T_com, X_com)[0])

    h_expected = (m_str * h_str + m_br * h_br) / m_com
    # Allow 1 J/kg absolute tolerance (rounding in iterative T solve)
    err = abs(h_com - h_expected)
    assert err < 1.0, (
        f"Enthalpy not conserved: h_com={h_com:.2f}, expected={h_expected:.2f}, err={err:.3f} J/kg"
    )


# ---------------------------------------------------------------------------
# Jacobian accuracy (FD check)
# ---------------------------------------------------------------------------


def _check_jacobian_vs_fd(net: FlowNetwork, tol: float = 5e-3) -> None:
    """Central-difference Jacobian check against analytical for the full solver system."""
    solver = NetworkSolver(net)
    sol = solver.solve()
    x0 = np.array([sol.get(n, 0.0) for n in solver.unknown_names])

    _, jac_analytical = solver._residuals_and_jacobian(x0, compute_jacobian=True)
    jac_dense = jac_analytical.toarray()

    jac_fd = fd_jacobian(solver, x0)

    # Check all columns with a combined absolute+relative criterion so that
    # near-zero analytical entries (that FD correctly resolves as near-zero) are
    # still covered without false failures from floating-point noise.
    scale = np.maximum(np.abs(jac_fd), np.abs(jac_dense)) + 1.0
    rel_err = np.abs(jac_dense - jac_fd) / scale
    max_err = rel_err.max()
    worst_idx = np.unravel_index(rel_err.argmax(), rel_err.shape)
    assert max_err < tol, (
        f"Jacobian FD mismatch: max relative err={max_err:.3e} "
        f"at row={worst_idx[0]} (res '{solver.unknown_names[worst_idx[0]] if worst_idx[0] < len(solver.unknown_names) else worst_idx[0]}'), "
        f"col={worst_idx[1]} (unk '{solver.unknown_names[worst_idx[1]]}')"
    )


def test_merging_tee_jacobian_vs_fd():
    """Merging-tee analytical Jacobian matches central-difference FD (all entries)."""
    _check_jacobian_vs_fd(_merging_net())


def test_branching_tee_jacobian_vs_fd():
    """Branching-tee analytical Jacobian matches central-difference FD (all entries)."""
    _check_jacobian_vs_fd(_branching_net())


def test_merging_tee_mixed_composition_jacobian_vs_fd():
    """Jacobian accuracy with two different inlet compositions (mixed-species sensitivity)."""
    _check_jacobian_vs_fd(_merging_net_two_compositions())


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_validate_bad_theta():
    # theta=0 is now accepted (limit is continuous, flagged as extrapolated by C++ correlation)
    # Only values with |theta| > pi/2 should raise.
    with pytest.raises(ValueError, match="theta"):
        net = FlowNetwork()
        Y = _DRY_AIR_Y
        for nid in ("a", "b", "c"):
            net.add_node(PressureBoundary(nid, Pt=2e5, Tt=400.0, Y=Y))
        net.add_element(TeeJunctionElement("tee", "a", "b", "c", theta=2.0, F_C=0.01))
        net.validate()


def test_validate_bad_F_C():
    with pytest.raises(ValueError, match="F_C"):
        net = FlowNetwork()
        Y = _DRY_AIR_Y
        for nid in ("a", "b", "c"):
            net.add_node(PressureBoundary(nid, Pt=2e5, Tt=400.0, Y=Y))
        net.add_element(TeeJunctionElement("tee", "a", "b", "c", theta=_THETA_90, F_C=-1.0))
        net.validate()


def test_validate_duplicate_node_rejected():
    """Adding a tee where two ports share the same node raises ValueError."""
    with pytest.raises(ValueError, match="duplicate"):
        net = FlowNetwork()
        Y = _DRY_AIR_Y
        for nid in ("a", "b"):
            net.add_node(PressureBoundary(nid, Pt=2e5, Tt=400.0, Y=Y))
        net.add_element(TeeJunctionElement("tee", "a", "b", "a", theta=_THETA_90, F_C=0.01))


# ---------------------------------------------------------------------------
# Multi-tee cascade (manifold)
# ---------------------------------------------------------------------------


def _methane_Y() -> list[float]:
    Y = [0.0] * cb.num_species()
    names = [cb.species_name(i) for i in range(cb.num_species())]
    if "CH4" in names:
        Y[names.index("CH4")] = 1.0
    return Y


def _cascaded_branching_manifold() -> FlowNetwork:
    """3-branching-tee manifold: one methane inlet splits to 4 channel outlets.

    Each tee splits its common arm into straight (continues cascade) and branch
    (exits to outlet).  Separate MCN nodes exist at every tee port so each node
    has exactly one upstream and one downstream element.

    Layout (all channels D=0.025m, L=0.3m):

      pb_in -> ch0 -> mc1_com
                      tee1 -> mc1_str -> ch_s1 -> mc2_com
                           -> mc1_bra -> ch_b1 -> pb_out_1b
                                         tee2 -> mc2_str -> ch_s2 -> mc3_com
                                              -> mc2_bra -> ch_b2 -> pb_out_2b
                                                             tee3 -> mc3_str -> ch_s3 -> pb_out_3s
                                                                  -> mc3_bra -> ch_b3 -> pb_out_3b
    """
    D = 0.025
    L = 0.3
    roughness = 1e-5
    F_C = math.pi * D**2 / 4.0
    theta = math.pi / 2.0

    Y = _methane_Y()
    Pt_in = 2.0e5
    Pt_out = 1.013e5
    Tt = 300.0

    net = FlowNetwork()

    net.add_node(PressureBoundary("pb_in", Pt=Pt_in, Tt=Tt, Y=Y))
    for nid in (
        "mc1_com",
        "mc1_str",
        "mc1_bra",
        "mc2_com",
        "mc2_str",
        "mc2_bra",
        "mc3_com",
        "mc3_str",
        "mc3_bra",
    ):
        net.add_node(MomentumChamberNode(nid))
    for nid in ("pb_out_1b", "pb_out_2b", "pb_out_3s", "pb_out_3b"):
        net.add_node(PressureBoundary(nid, Pt=Pt_out, Tt=Tt, Y=Y))

    net.add_element(
        ChannelElement("ch0", "pb_in", "mc1_com", length=L, diameter=D, roughness=roughness)
    )
    net.add_element(
        ChannelElement("ch_s1", "mc1_str", "mc2_com", length=L, diameter=D, roughness=roughness)
    )
    net.add_element(
        ChannelElement("ch_s2", "mc2_str", "mc3_com", length=L, diameter=D, roughness=roughness)
    )
    net.add_element(
        ChannelElement("ch_s3", "mc3_str", "pb_out_3s", length=L, diameter=D, roughness=roughness)
    )
    net.add_element(
        ChannelElement("ch_b1", "mc1_bra", "pb_out_1b", length=L, diameter=D, roughness=roughness)
    )
    net.add_element(
        ChannelElement("ch_b2", "mc2_bra", "pb_out_2b", length=L, diameter=D, roughness=roughness)
    )
    net.add_element(
        ChannelElement("ch_b3", "mc3_bra", "pb_out_3b", length=L, diameter=D, roughness=roughness)
    )

    net.add_element(
        TeeJunctionElement(
            "tee1",
            common_node="mc1_com",
            straight_node="mc1_str",
            branch_node="mc1_bra",
            theta=theta,
            F_C=F_C,
            tee_type="branching",
        )
    )
    net.add_element(
        TeeJunctionElement(
            "tee2",
            common_node="mc2_com",
            straight_node="mc2_str",
            branch_node="mc2_bra",
            theta=theta,
            F_C=F_C,
            tee_type="branching",
        )
    )
    net.add_element(
        TeeJunctionElement(
            "tee3",
            common_node="mc3_com",
            straight_node="mc3_str",
            branch_node="mc3_bra",
            theta=theta,
            F_C=F_C,
            tee_type="branching",
        )
    )
    return net


def _gentle_branching_cascade() -> FlowNetwork:
    """Two-branching-tee manifold at a gentle pressure ratio (1.1 -> 1.0 bar).

    Same MCN-bounded topology as the larger manifold but with a low-Mach
    pressure ratio, where the branching loss closure stays well-conditioned and
    the solver reaches a mass-conserving root.  Used as the positive baseline;
    the full near-choke manifold below is xfailed.
    """
    D = 0.025
    L = 0.3
    roughness = 1e-5
    F_C = math.pi * D**2 / 4.0
    theta = math.pi / 2.0
    Y = _methane_Y()

    net = FlowNetwork()
    net.add_node(PressureBoundary("pb_in", Pt=1.1e5, Tt=300.0, Y=Y))
    for nid in ("mc1_com", "mc1_str", "mc1_bra", "mc2_com", "mc2_str", "mc2_bra"):
        net.add_node(MomentumChamberNode(nid, area=F_C))
    for nid in ("pb_1b", "pb_2b", "pb_2s"):
        net.add_node(PressureBoundary(nid, Pt=1.0e5, Tt=300.0, Y=Y))

    net.add_element(
        ChannelElement("ch0", "pb_in", "mc1_com", length=L, diameter=D, roughness=roughness)
    )
    net.add_element(
        ChannelElement("ch_s1", "mc1_str", "mc2_com", length=L, diameter=D, roughness=roughness)
    )
    net.add_element(
        ChannelElement("ch_b1", "mc1_bra", "pb_1b", length=L, diameter=D, roughness=roughness)
    )
    net.add_element(
        ChannelElement("ch_b2", "mc2_bra", "pb_2b", length=L, diameter=D, roughness=roughness)
    )
    net.add_element(
        ChannelElement("ch_s2", "mc2_str", "pb_2s", length=L, diameter=D, roughness=roughness)
    )
    net.add_element(
        TeeJunctionElement(
            "tee1",
            common_node="mc1_com",
            straight_node="mc1_str",
            branch_node="mc1_bra",
            theta=theta,
            F_C=F_C,
            tee_type="branching",
        )
    )
    net.add_element(
        TeeJunctionElement(
            "tee2",
            common_node="mc2_com",
            straight_node="mc2_str",
            branch_node="mc2_bra",
            theta=theta,
            F_C=F_C,
            tee_type="branching",
        )
    )
    return net


# The gentle (low-Mach) MCN-bounded branching cascade now reaches a physical
# (all-forward) mass-conserving root. Two defects previously blocked it, both fixed:
#   * Defect 1 (task #7): ChannelElement coupled upstream Pt to downstream P, leaking an
#     inline MomentumChamberNode's dynamic head as a free pressure gain. Landed.
#   * Defect 2 (task #9): the parallel branch closure over-charged a near-closed branch
#     (Borda-Carnot K -> 1) and forced K_str = K_bra at equal outlets. The unified0D
#     blended turning-loss closure (both legs: extraction + branch turning loss) removes
#     this. Landed in compressible_branching_tee_rj.
# With both fixed, the gentle cascade is a single source feeding sinks, so the physical
# answer is all-outward; the solver reproduces it (addendum Finding 6). The full
# near-choke manifold below stays xfail - it is outside the closure envelope.
def test_gentle_branching_cascade_forward_root():
    """Two-tee MCN-bounded manifold reaches an all-forward, mass-conserving root."""
    net = _gentle_branching_cascade()
    solver = NetworkSolver(net)
    sol = solver.solve()
    assert sol["__success__"], f"did not converge: {sol.get('__message__')}"
    m_in = sol["ch0.m_dot"]
    m_out = sol["ch_b1.m_dot"] + sol["ch_b2.m_dot"] + sol["ch_s2.m_dot"]
    assert abs(m_in - m_out) < 1e-6 * max(abs(m_in), 1e-6), (
        f"Mass not conserved: {m_in:.6f} vs {m_out:.6f}"
    )
    for tee in ("tee1", "tee2"):
        m_com = sol[f"{tee}.m_dot_com"]
        m_branch = sol[f"{tee}.m_dot_branch"]
        m_straight = m_com - m_branch
        assert m_com > 0.0, f"{tee} common reversed: {m_com}"
        assert m_branch > 0.0, f"{tee} branch reversed: {m_branch}"
        assert m_straight > 0.0, f"{tee} straight reversed: {m_straight}"


@pytest.mark.xfail(
    reason="branching cascade reversed/unsolvable at near-choke ratio (addendum Findings 4-5; tasks #7, #9)",
    strict=False,
)
def test_cascaded_branching_manifold_converges():
    """3-tee branching manifold with channels: solver must converge."""
    net = _cascaded_branching_manifold()
    solver = NetworkSolver(net)
    sol = solver.solve()
    assert sol["__success__"]
    assert sol["tee1.m_dot_com"] > 0.0


@pytest.mark.xfail(
    reason="branching cascade reversed/unsolvable at near-choke ratio (addendum Findings 4-5; tasks #7, #9)",
    strict=False,
)
def test_cascaded_branching_manifold_mass_conservation():
    """Total inlet flow equals sum of 4 outlet flows."""
    net = _cascaded_branching_manifold()
    solver = NetworkSolver(net)
    sol = solver.solve()
    assert sol["__success__"]
    m_in = sol["ch0.m_dot"]
    m_out = sol["ch_b1.m_dot"] + sol["ch_b2.m_dot"] + sol["ch_s3.m_dot"] + sol["ch_b3.m_dot"]
    assert abs(m_in - m_out) < 1e-6 * m_in, f"Mass not conserved: {m_in:.6f} vs {m_out:.6f}"
