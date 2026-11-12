"""EnergyBoundary(Tt_target=...) (#792): a heat source/sink with a set outlet
TOTAL temperature. The heat flow is a solver unknown ``{eb}.Q`` closed by one
temperature row, Tt_node - Tt_target = 0; its Jacobian column is C++'s
dT_mix/dQ relayed through the node's propagation chain.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest
import test_combustor_regression as tcr
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    CombustorNode,
    EnergyBoundary,
    FlowNetwork,
    LinearThetaFractionLoss,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    WallNode,
)
from combaero.network.components import validate_energy_boundaries
from combaero.network.second_law import LEDGER_CLOSURE_TOL

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _chain(
    eb: EnergyBoundary | None,
    node: str = "plenum",
    regime: str = "compressible",
    Pa: float = 1.3e5,
    Pb: float = 1.0e5,
) -> FlowNetwork:
    """a -> o1 -> p (the heated node) -> o2 -> q -> o3 -> b, all at 400 K."""
    g = FlowNetwork()
    for name, Pt in (("a", Pa), ("b", Pb)):
        g.add_node(PressureBoundary(name, Pt=Pt, Tt=400.0, Y=Y_AIR))
    p = PlenumNode("p") if node == "plenum" else MomentumChamberNode("p", area=0.002, regime=regime)
    g.add_node(p)
    g.add_node(PlenumNode("q"))
    for eid, a, b in (("o1", "a", "p"), ("o2", "p", "q"), ("o3", "q", "b")):
        g.add_element(
            OrificeElement(eid, a, b, Cd=0.8, diameter=0.03, correlation="fixed", regime=regime)
        )
    if eb is not None:
        p.add_energy_boundary(eb)
    return g


def _jac_err(s: NetworkSolver, x: np.ndarray) -> float:
    """RESIDUAL-ROW-SCALING: compare J_ij |x_j| (see _build_residual_scales)."""
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    return float(np.max(np.abs(Js - fds) / scale))


# --- construction and the unknown --------------------------------------------


def test_the_target_mode_adds_one_unknown_and_one_temperature_row() -> None:
    eb = EnergyBoundary("h", Tt_target=600.0)
    assert eb.is_target and eb.unknown_name == "h.Q" and eb.heat == 0.0
    g = _chain(eb)
    s = NetworkSolver(g)
    s._build_x0()
    assert "h.Q" in s.unknown_names
    assert g.nodes["p"].unknowns() == ["p.P", "p.Pt", "h.Q"]
    x0 = s._build_x0()
    kinds = s._row_kinds(x0)
    assert ("node 'p' Tt_target", "T") in kinds
    # the temperature row is weighted by a temperature scale [K]
    D_f = s._build_residual_scales(x0)
    i_T = kinds.index(("node 'p' Tt_target", "T"))
    assert D_f[i_T] == pytest.approx(s._reference_T_scale())
    assert s._reference_T_scale() == pytest.approx(200.0)  # targets 600, boundaries 400


def test_the_start_heat_is_the_heat_to_the_target_at_the_start_flow() -> None:
    s = NetworkSolver(_chain(EnergyBoundary("h", Tt_target=600.0)))
    x0 = s._build_x0()
    q0 = x0[s._name_to_index["h.Q"]]
    m0 = x0[s._name_to_index["o1.m_dot"]]
    X = cb.mass_to_mole(Y_AIR)
    assert q0 == pytest.approx(m0 * (cb.h_mass(600.0, X) - cb.h_mass(400.0, X)), rel=1e-9)


# --- round trips --------------------------------------------------------------


@pytest.mark.parametrize("node", ["plenum", "chamber"])
@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
@pytest.mark.parametrize("Q1", [5.0e4, -2.0e4], ids=["heating", "cooling"])
def test_round_trip_recovers_the_heat(node: str, regime: str, Q1: float) -> None:
    """Solve with Q = Q1 for Tt1; solve with Tt_target = Tt1 and get Q1 back."""
    r1 = NetworkSolver(_chain(EnergyBoundary("h", Q=Q1), node, regime)).solve()
    assert r1["__success__"]
    Tt1 = r1["p.Tt"]
    assert (Tt1 > 400.0) == (Q1 > 0.0)
    assert r1["h.Q"] == Q1 and r1["p.Q"] == Q1  # reported in fixed mode too

    r2 = NetworkSolver(_chain(EnergyBoundary("h", Tt_target=Tt1), node, regime)).solve()
    assert r2["__success__"], r2["__message__"]
    assert r2["h.Q"] == pytest.approx(Q1, rel=1e-8)
    assert r2["p.Q"] == pytest.approx(Q1, rel=1e-8)
    assert r2["p.Tt"] == pytest.approx(Tt1, abs=1e-6)
    for k in ("o1.m_dot", "p.P", "q.Tt"):
        assert r2[k] == pytest.approx(r1[k], rel=1e-7)


def test_round_trip_on_a_combustor_is_the_heat_after_the_reaction() -> None:
    """On a CombustorNode the target is the EXIT Tt: Q is the heat taken out
    after the reaction to reach it (here a loss), the fuel flow untouched."""
    net, comb = tcr._make_effective_area_inlet_network()
    comb.add_energy_boundary(EnergyBoundary("loss", Q=-1.5e5))
    r1 = NetworkSolver(net).solve()
    assert r1["__success__"]

    net2, comb2 = tcr._make_effective_area_inlet_network()
    comb2.add_energy_boundary(EnergyBoundary("loss", Tt_target=r1["comb.Tt"]))
    r2 = NetworkSolver(net2).solve()
    assert r2["__success__"], r2["__message__"]
    assert r2["loss.Q"] == pytest.approx(-1.5e5, rel=1e-8)
    assert r2["comb.Q"] == pytest.approx(-1.5e5, rel=1e-8)
    assert r2["comb.phi"] == pytest.approx(r1["comb.phi"], rel=1e-10)


def test_a_target_beside_a_fixed_heat_solves_for_the_rest() -> None:
    """A fixed Q boundary and a Tt_target one on the same node: the target's Q
    is what the fixed one leaves; the node reports the sum."""
    r1 = NetworkSolver(_chain(EnergyBoundary("h", Q=3.0e4))).solve()
    g = _chain(EnergyBoundary("fixed", Q=1.0e4))
    g.nodes["p"].add_energy_boundary(EnergyBoundary("trim", Tt_target=r1["p.Tt"]))
    r2 = NetworkSolver(g).solve()
    assert r2["__success__"]
    assert r2["fixed.Q"] == 1.0e4
    assert r2["trim.Q"] == pytest.approx(2.0e4, rel=1e-8)
    assert r2["p.Q"] == pytest.approx(3.0e4, rel=1e-8)


def test_reverse_flow_through_the_node() -> None:
    """Flow from b to a: the node mixes what now flows INTO it (from q), and
    the target still fixes Q, the same heat the fixed-Q twin applies."""
    r1 = NetworkSolver(_chain(EnergyBoundary("h", Q=2.0e4), Pa=1.0e5, Pb=1.3e5)).solve()
    assert r1["__success__"] and r1["o1.m_dot"] < 0.0 and r1["o2.m_dot"] < 0.0
    g = _chain(EnergyBoundary("h", Tt_target=r1["p.Tt"]), Pa=1.0e5, Pb=1.3e5)
    s = NetworkSolver(g)
    r2 = s.solve()
    assert r2["__success__"], r2["__message__"]
    assert r2["o1.m_dot"] < 0.0
    assert r2["h.Q"] == pytest.approx(2.0e4, rel=1e-8)
    assert _jac_err(s, np.array(r2["__x_solution__"])) < 1e-6


# --- refusals -----------------------------------------------------------------


@pytest.mark.parametrize("kw", [{"Q": 100.0}, {"fraction": -0.01}])
def test_a_target_with_a_heat_or_fraction_is_refused(kw: dict) -> None:
    with pytest.raises(ValueError, match="over-specified"):
        EnergyBoundary("h", Tt_target=600.0, **kw)


@pytest.mark.parametrize("T", [0.0, -5.0, math.nan, math.inf])
def test_a_non_physical_target_is_refused(T: float) -> None:
    with pytest.raises(ValueError, match="positive, finite"):
        EnergyBoundary("h", Tt_target=T)


def test_a_heat_set_after_construction_is_refused_at_network_build() -> None:
    eb = EnergyBoundary("h", Tt_target=600.0)
    g = _chain(eb)
    eb.Q = 500.0
    with pytest.raises(ValueError, match="over-specified"):
        NetworkSolver(g).solve()


def test_two_targets_on_one_node_are_refused() -> None:
    p = PlenumNode("p")
    p.add_energy_boundary(EnergyBoundary("t1", Tt_target=600.0))
    with pytest.raises(ValueError, match="second time"):
        p.add_energy_boundary(EnergyBoundary("t2", Tt_target=700.0))
    # and at network build, should the list be filled directly
    p.energy_boundaries.append(EnergyBoundary("t2", Tt_target=700.0))
    with pytest.raises(ValueError, match="at most one"):
        validate_energy_boundaries(p)


def test_a_target_on_a_node_that_cannot_honour_it_is_refused() -> None:
    """A WallNode (a closed end) has no energy balance a heat enters."""
    w = WallNode("w")
    w.energy_boundaries = [EnergyBoundary("h", Tt_target=600.0)]
    with pytest.raises(ValueError, match="cannot honour it"):
        validate_energy_boundaries(w)


def test_every_target_node_type_is_accepted() -> None:
    for node in (PlenumNode("p"), MomentumChamberNode("m", area=0.01), CombustorNode("c")):
        node.add_energy_boundary(EnergyBoundary("h", Tt_target=900.0))
        validate_energy_boundaries(node)


# --- zero flow ----------------------------------------------------------------


@pytest.mark.no_closure_check  # nothing converges: there is no balance to close
@pytest.mark.parametrize("dP", [0.0, 1e-6], ids=["no_drive", "micro_drive"])
def test_no_through_flow_leaves_the_heat_undetermined_and_says_so(dP: float) -> None:
    """At zero flow Q does not reach the node temperature (dT/dQ = phi(m)/cp
    = 0): any Q, or none, meets the target. The solve does not succeed, names
    the boundary, and reports no Q (NaN) -- never the iterate's Q."""
    s = NetworkSolver(_chain(EnergyBoundary("h", Tt_target=600.0), Pa=1.0e5 + dP, Pb=1.0e5))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=30.0)
    assert not r["__success__"]
    assert r["__heat_undetermined__"] == ["h"]
    assert "no through-flow" in r["__message__"]
    assert math.isnan(r["h.Q"]) and math.isnan(r["p.Q"])


# --- Jacobian -----------------------------------------------------------------


@pytest.mark.parametrize("node", ["plenum", "chamber"])
@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
def test_the_jacobian_matches_differences(node: str, regime: str) -> None:
    """The Q column reaches the target row AND, through the relay, every row
    downstream that reads p's temperature (o2, q, o3)."""
    s = NetworkSolver(_chain(EnergyBoundary("h", Tt_target=700.0), node, regime))
    r = s.solve()
    assert r["__success__"]
    x = np.array(r["__x_solution__"])
    assert _jac_err(s, x) < 1e-6
    # and away from the root, where the target row is not yet zero
    x_off = x.copy()
    x_off[s._name_to_index["h.Q"]] *= 0.7
    x_off[s._name_to_index["o1.m_dot"]] *= 1.1
    assert _jac_err(s, x_off) < 1e-6
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()
    col = J[:, s._name_to_index["h.Q"]]
    assert np.count_nonzero(col) > 1  # not only the target row


def test_the_jacobian_matches_differences_on_a_combustor_with_a_theta_loss() -> None:
    """The burned temperature feeds a theta-sourced loss downstream."""
    net, comb = tcr._make_effective_area_inlet_network(
        loss_correlation=LinearThetaFractionLoss(k=0.01, xi0=0.02)
    )
    comb.add_energy_boundary(EnergyBoundary("loss_q", Tt_target=1500.0))
    s = NetworkSolver(net)
    r = s.solve()
    assert r["__success__"], r["__message__"]
    assert r["comb.Tt"] == pytest.approx(1500.0, abs=1e-6)
    assert _jac_err(s, np.array(r["__x_solution__"])) < 2e-5


# --- second law and the balance ledger -----------------------------------------


@pytest.mark.parametrize("Q1", [5.0e4, -2.0e4], ids=["heating", "cooling"])
def test_the_ledgers_use_the_solved_heat(Q1: float) -> None:
    """The ledger counts the SOLVED Q (the declared Q of a target is 0, which
    would leave the balance open by Q1 -- far above its bound), and the second
    law sees the same heat as the fixed-Q twin."""
    r1 = NetworkSolver(_chain(EnergyBoundary("h", Q=Q1))).solve()
    r2 = NetworkSolver(_chain(EnergyBoundary("h", Tt_target=r1["p.Tt"]))).solve()
    assert r2["__success__"]
    led = r2["__ledger__"]
    assert led["Q_given"] == pytest.approx(r2["h.Q"], rel=1e-12)
    assert led["Q_given"] == pytest.approx(Q1, rel=1e-8)
    assert abs(led["energy_closure_rel"]) < 1e-9
    h_scale = max(abs(led["H_in"]), abs(led["H_out"]))
    assert abs(Q1) / h_scale > 1e3 * LEDGER_CLOSURE_TOL  # declared 0 would fail
    v1, v2 = r1["__second_law__"], r2["__second_law__"]
    assert v1.keys() == v2.keys()
    for eid in v1:
        assert v2[eid]["verdict"] == v1[eid]["verdict"]
        if "sigma_norm" in v1[eid]:
            assert v2[eid]["sigma_norm"] == pytest.approx(v1[eid]["sigma_norm"], rel=1e-6)


def test_a_fraction_boundary_reports_the_heat_it_applied() -> None:
    """{eb}.Q in fraction mode is the fraction heat; {node}.Q sums every mode."""
    g = _chain(EnergyBoundary("f", fraction=-0.02))
    g.nodes["p"].add_energy_boundary(EnergyBoundary("h", Q=500.0))
    r = NetworkSolver(g).solve()
    assert r["__success__"]
    assert r["f.Q"] == pytest.approx(r["p.Q_fraction"], rel=1e-12)
    assert r["f.Q"] < 0.0
    assert r["p.Q"] == pytest.approx(r["f.Q"] + 500.0, rel=1e-12)
