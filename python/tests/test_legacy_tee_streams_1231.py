"""The legacy TeeJunctionElement feeds its port nodes by the sign of its legs (#1231).

The tee used to hand its sink nodes signed leg flows. A merging tee whose
branch leg drew flow out of the common node gave that node a negative
stream, and the branch node, which the leg now flows into, received nothing
and kept the 300 K fallback. At a reversed root the energy ledger closed only
because the signed mean subtracted that fallback state again, and the common
node reached up to 17,800 K from 300-900 K inputs.

A tee with a leg against its tee_type is now a reversed junction, fed as the
#517 junctions are (``NetworkSolver._junction_streams``): each port node
takes the streams arriving at it. The tee has a second-law control volume
like the other junctions (#545), and a root with a leg against its tee_type
is outside Bassett's coefficients for that pattern, so it is demoted with
the legs named (``verify_solution_consistent``).
"""

from __future__ import annotations

import math
import warnings

import numpy as np
from _fd_jacobian import fd_jacobian, worst

import combaero as cb
from combaero import _solver_tools
from combaero.network import FlowNetwork, NetworkSolver
from combaero.network.components import (
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    OrificeElement,
    PressureBoundary,
    TeeJunctionElement,
)

AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _fuel_air() -> list[float]:
    Y = [0.98 * y for y in AIR]
    Y[cb.species.species_index("CH4")] += 0.02
    return Y


def _merging(
    Pis: float,
    Pib: float,
    Po: float,
    Tis: float,
    Tib: float,
    To: float,
    F_C: float,
    psi: float,
    Y_is: list[float] | None = None,
) -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(PressureBoundary("is", Pt=Pis, Tt=Tis, Y=Y_is or AIR))
    g.add_node(PressureBoundary("ib", Pt=Pib, Tt=Tib, Y=AIR))
    g.add_node(PressureBoundary("o", Pt=Po, Tt=To, Y=AIR))
    g.add_node(MomentumChamberNode("ps", area=F_C))
    g.add_node(MomentumChamberNode("pb", area=F_C / psi))
    g.add_node(MomentumChamberNode("pc", area=F_C))
    g.add_element(LosslessConnectionElement("ls", "is", "ps"))
    g.add_element(LosslessConnectionElement("lb", "ib", "pb"))
    g.add_element(LosslessConnectionElement("lc", "pc", "o"))
    g.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="pc",
            straight_node="ps",
            branch_node="pb",
            theta=math.radians(90.0),
            F_C=F_C,
            psi=psi,
            tee_type="merging",
        )
    )
    return g


def _solve(g: FlowNetwork) -> dict:
    s = NetworkSolver(g)
    s.homotopy_retry = True
    s.bc_ramp_retry = True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return s.solve(timeout=20)


def _at(g: FlowNetwork, m_com: float, m_b: float) -> tuple[NetworkSolver, np.ndarray]:
    """A resolved solver and its start vector with the tee's legs set (the
    connections carry the same flows, so the mass rows hold)."""
    s = NetworkSolver(g)
    g.resolve_all_topology()
    x = np.asarray(s._build_x0(), dtype=float)
    idx = s._name_to_index
    for name, m in (
        ("tee.m_dot_com", m_com),
        ("tee.m_dot_branch", m_b),
        ("lc.m_dot", m_com),
        ("lb.m_dot", m_b),
        ("ls.m_dot", m_com - m_b),
    ):
        x[idx[name]] = m
    return s, x


def _net() -> FlowNetwork:
    return _merging(2.0e5, 1.98e5, 1.97e5, 500.0, 900.0, 700.0, 1e-2, 2.0, Y_is=_fuel_air())


def test_a_reversed_branch_delivers_the_straight_stream_to_its_node() -> None:
    """Branch leg reversed: the straight leg is the tee's only supplier, so
    the common node and the branch node both take its gas. The branch node
    used to receive nothing (300 K, no fuel) and the common node a signed
    mean of +1.3 kg/s at 500 K and -0.3 kg/s at that 300 K."""
    s, x = _at(_net(), 1.0, -0.3)
    s._propagate_states(x)
    T_s, Y_s, _ = s._derived_states["ps"]
    ich4 = cb.species.species_index("CH4")
    assert math.isclose(T_s, 500.0, rel_tol=1e-12)
    assert math.isclose(Y_s[ich4], 0.02, rel_tol=1e-12)
    for nid in ("pb", "pc"):
        T, Y, _ = s._derived_states[nid]
        assert math.isclose(T, T_s, rel_tol=1e-12), nid
        assert np.allclose(Y, Y_s, rtol=1e-12, atol=1e-15), nid
    # And the node inflows the heat checks read: one stream each, at the leg's flow.
    assert [m for m, _, _ in s._node_inflows("pb", x)] == [0.3]
    assert [m for m, _, _ in s._node_inflows("pc", x)] == [1.0]


def test_a_forward_tee_still_mixes_its_legs_at_the_common_node() -> None:
    s, x = _at(_net(), 1.0, 0.3)
    s._propagate_states(x)
    T_c = s._derived_states["pc"][0]
    mixed = _solver_tools.mixer_from_streams_and_jacobians(
        [
            cb.MassStream(0.7, 500.0, 2.0e5, _fuel_air()),
            cb.MassStream(0.3, 900.0, 1.98e5, AIR),
        ],
        Q=0.0,
        fraction=0.0,
    )
    assert math.isclose(T_c, mixed.Tt_mix, rel_tol=1e-12)


def _fd_error(m_com: float, m_b: float) -> tuple[float, str]:
    s, x = _at(_net(), m_com, m_b)
    _, J = s._residuals_and_jacobian(x)
    return worst(J, fd_jacobian(s, x), x, list(s.unknown_names))


def test_the_jacobian_matches_a_difference_with_a_leg_reversed() -> None:
    """The streams' weights are the legs' flows, linear in m_com and m_b
    (the straight leg in both), so the relay is exact there too."""
    for m_com, m_b in ((1.0, -0.3), (-0.4, 0.1), (-0.5, -0.2)):
        e, where = _fd_error(m_com, m_b)
        assert e < 1e-6, (m_com, m_b, where)


def test_case_108_is_refused_and_its_ledger_closes() -> None:
    """The #1230 review's case 108: the branch leg reversed and the straight
    inlet's air leaving at 203.5 kPa from 200 kPa. The tee's control volume
    destroys entropy, so the second law fails it, and the leg against the
    tee_type puts it outside the tee's model. Base reported it converged
    with its second law not checked."""
    r = _solve(
        _merging(200000.0, 197704.39100378475, 203503.56229127114, 300.0, 900.0, 700.0, 0.03, 2.0)
    )
    assert r["__converged__"]
    assert r["tee.m_dot_branch"] < 0.0 < r["tee.m_dot_com"]
    sl = r["__second_law__"]["tee"]
    assert sl["verdict"] == "fail" and sl["sigma_norm_high"] < -1e-3
    assert r["__second_law_failures__"] == ["tee"]
    assert r["__outcome__"] == "inconsistent"
    assert r["__inconsistent_elements__"] == ["tee"]
    assert "runs leg(s) ['pb'] against its tee_type" in r["__message__"]
    assert r["__element_diag__"]["tee"]["n_reversed_ports"] == 1.0
    assert r["__element_diag__"]["tee"]["correlation_extrapolated"] == 1.0
    # The branch node carries the straight leg's gas, which it delivers.
    assert math.isclose(r["pb.Tt"], r["ps.Tt"], rel_tol=1e-12)
    L = r["__ledger__"]
    assert abs(L["mass_closure_rel"]) < 1e-12 and abs(L["energy_closure_rel"]) < 1e-12


def test_a_leg_against_the_tee_type_is_outside_the_model() -> None:
    """verify_solution_consistent, both tee types, at the tolerance."""
    merging = TeeJunctionElement("t", "c", "s", "b", theta=0.5, F_C=1e-2, tee_type="merging")
    branching = TeeJunctionElement("t", "c", "s", "b", theta=0.5, F_C=1e-2, tee_type="branching")

    def ok(e: TeeJunctionElement, m_com: float, m_b: float) -> bool:
        return e.verify_solution_consistent({"t.m_dot_com": m_com, "t.m_dot_branch": m_b})

    for e in (merging, branching):
        assert ok(e, 1.0, 0.4)
        assert ok(e, 1.0, 0.0) and ok(e, 1.0, 1.0)  # a leg at rest is in the pattern
        assert ok(e, 1.0, -1e-7)  # within TEE_REVERSE_FLOW_TOL of the largest leg
        assert not ok(e, 1.0, -1e-3)  # branch reversed
        assert not ok(e, 1.0, 1.2)  # straight reversed
        assert not ok(e, -1.0, -0.4)  # every leg reversed
    assert merging.verify_solution_consistent({}) is True  # missing keys pass


def test_a_branching_tee_at_boundaries_books_the_mixed_stream_it_delivers() -> None:
    """A branching tee wired straight to boundaries, its straight boundary
    above the common one: both supply the branch, which receives their mix.
    The ledger books the branch's outflow at that mix; it used to take the
    common boundary's enthalpy alone, and the 300 K / 800 K energy balance
    stayed open."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("c", Pt=2.00e5, Tt=800.0, Y=AIR))
    g.add_node(PressureBoundary("s", Pt=2.01e5, Tt=300.0, Y=AIR))
    g.add_node(PressureBoundary("b", Pt=1.90e5, Tt=300.0, Y=AIR))
    g.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="c",
            straight_node="s",
            branch_node="b",
            theta=math.radians(90.0),
            F_C=1e-2,
            psi=1.0,
            tee_type="branching",
        )
    )
    r = _solve(g)
    assert r["__converged__"], r["__message__"]
    m_com, m_b = r["tee.m_dot_com"], r["tee.m_dot_branch"]
    assert 0.0 < m_com < m_b  # the straight leg runs into the tee
    assert r["__outcome__"] == "inconsistent"
    assert r["__inconsistent_elements__"] == ["tee"]
    L = r["__ledger__"]
    assert abs(L["mass_closure_rel"]) < 1e-12 and abs(L["energy_closure_rel"]) < 1e-12


def test_a_forward_branching_tee_into_colder_boundaries_passes_the_second_law() -> None:
    """A boundary holds its set state, not the gas the tee delivers into it:
    each outflow leg leaves at the state of what supplies it. Taken at the
    boundaries' own 300 K, the 800 K supply read as sigma_norm -3.5 (stress
    legacy_tee seed 3)."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("c", Pt=1.5e5, Tt=800.0, Y=AIR))
    g.add_node(PressureBoundary("s", Pt=1.49925e5, Tt=300.0, Y=AIR))
    g.add_node(PressureBoundary("b", Pt=1.485e5, Tt=300.0, Y=AIR))
    g.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="c",
            straight_node="s",
            branch_node="b",
            theta=math.radians(60.0),
            F_C=3e-2,
            psi=1.0,
            tee_type="branching",
        )
    )
    r = _solve(g)
    assert r["__outcome__"] == "converged", r["__message__"]
    assert 0.0 < r["tee.m_dot_branch"] < r["tee.m_dot_com"]
    v = r["__second_law__"]["tee"]
    assert v["verdict"] == "pass" and v["sigma_norm"] > 0.0


def _co2_air(f: float) -> list[float]:
    Y = [(1.0 - f) * y for y in AIR]
    Y[cb.species.species_index("CO2")] += f
    return Y


def _branching_mc(side: str | None = None) -> FlowNetwork:
    """A branching tee between three chambers fed by lossless connections,
    each boundary its own gas; ``side`` injects a hot CO2-laden stream into
    that port chamber through an orifice."""
    g = FlowNetwork()
    gases = {
        "c": (2.0e5, 650.0, _co2_air(0.05)),
        "s": (2.02e5, 420.0, _fuel_air()),
        "b": (1.99e5, 880.0, AIR),
    }
    for k, (P, T, Y) in gases.items():
        g.add_node(PressureBoundary("B" + k, Pt=P, Tt=T, Y=Y))
        g.add_node(MomentumChamberNode("p" + k, area=1e-2 if k != "b" else 5e-3))
        a, b = ("B" + k, "p" + k) if k == "c" else ("p" + k, "B" + k)
        g.add_element(LosslessConnectionElement("l" + k, a, b))
    if side is not None:
        g.add_node(MassFlowBoundary("Mx", m_dot=0.05, Tt=1200.0, Y=_co2_air(0.1)))
        g.add_element(
            OrificeElement("ox", "Mx", "p" + side, Cd=0.7, diameter=0.03, correlation="fixed")
        )
    g.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="pc",
            straight_node="ps",
            branch_node="pb",
            theta=math.radians(60.0),
            F_C=1e-2,
            psi=2.0,
            tee_type="branching",
        )
    )
    return g


def _fd_at(g: FlowNetwork, legs: dict[str, float]) -> tuple[float, str]:
    s = NetworkSolver(g)
    g.resolve_all_topology()
    x = np.asarray(s._build_x0(), dtype=float)
    for name, m in legs.items():
        x[s._name_to_index[name]] = m
    _, J = s._residuals_and_jacobian(x)
    return worst(J, fd_jacobian(s, x), x, list(s.unknown_names))


def test_the_jacobian_carries_the_port_nodes_gas_constant_and_gamma() -> None:
    """The kernel takes each port node's R and gamma as numbers. With the
    straight leg of a branching tee reversed its node mixes gases that move
    with m_com and m_branch (0.2-0.7 % off on those columns without the
    chain), and a side stream of another gas into the common chamber moves
    its R with the side stream's flow (2 % off on ox.m_dot, #1235). A
    merging tee reads its suppliers' gamma too."""
    reversed_straight = {"tee.m_dot_com": 1.0, "tee.m_dot_branch": 1.3, "lc.m_dot": 1.0}
    reversed_straight |= {"ls.m_dot": -0.3, "lb.m_dot": 1.3}
    e, where = _fd_at(_branching_mc(), reversed_straight)
    assert e < 1e-6, where
    forward = {"tee.m_dot_com": 1.0, "tee.m_dot_branch": 0.4, "lc.m_dot": 0.95, "ox.m_dot": 0.05}
    forward |= {"ls.m_dot": 0.6, "lb.m_dot": 0.4}
    e, where = _fd_at(_branching_mc(side="c"), forward)
    assert e < 1e-6, where
    # Merging: a side stream of another gas and T into the straight chamber
    # moves that supplier's R and gamma.
    g = _net()
    g.add_node(MassFlowBoundary("Mx", m_dot=0.05, Tt=1200.0, Y=_co2_air(0.1)))
    g.add_element(OrificeElement("ox", "Mx", "ps", Cd=0.7, diameter=0.03, correlation="fixed"))
    legs = {"tee.m_dot_com": 1.0, "tee.m_dot_branch": 0.3, "lc.m_dot": 1.0, "lb.m_dot": 0.3}
    legs |= {"ls.m_dot": 0.65, "ox.m_dot": 0.05}
    e, where = _fd_at(g, legs)
    assert e < 1e-6, where


def test_a_side_stream_into_a_port_chamber_is_not_the_tees_outflow() -> None:
    """The tee's control volume books each leg it delivers at the mix it
    delivers, not at the receiving chamber's state, which also holds a side
    stream. A forward branching tee whose straight chamber takes 0.05 kg/s
    of 300 K air read its own outflow at that chamber's mix and failed the
    second law (the #1233 review: 8 forward roots base accepts, sigma_norm
    -0.16 to -2.1)."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("Bc", Pt=2.0e5, Tt=800.0, Y=AIR))
    g.add_node(PressureBoundary("Bs", Pt=1.999e5, Tt=300.0, Y=AIR))
    g.add_node(PressureBoundary("Bb", Pt=1.98e5, Tt=300.0, Y=AIR))
    g.add_node(MomentumChamberNode("pc", area=1e-2))
    g.add_node(MomentumChamberNode("ps", area=1e-2, main_inlet="tee"))
    g.add_node(MomentumChamberNode("pb", area=1e-2))
    g.add_element(LosslessConnectionElement("lc", "Bc", "pc"))
    g.add_element(LosslessConnectionElement("ls", "ps", "Bs"))
    g.add_element(LosslessConnectionElement("lb", "pb", "Bb"))
    g.add_node(MassFlowBoundary("Mx", m_dot=0.05, Tt=300.0, Y=AIR))
    g.add_element(OrificeElement("ox", "Mx", "ps", Cd=0.7, diameter=0.02, correlation="fixed"))
    g.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="pc",
            straight_node="ps",
            branch_node="pb",
            theta=math.radians(90.0),
            F_C=1e-2,
            psi=1.0,
            tee_type="branching",
        )
    )
    r = _solve(g)
    assert r["__outcome__"] == "converged", r["__message__"]
    assert 0.0 < r["tee.m_dot_branch"] < r["tee.m_dot_com"]
    # The chamber mixes the 800 K leg with the 300 K side stream; the tee's
    # outflow there is the 800 K leg.
    assert r["ps.Tt"] < 790.0
    v = r["__second_law__"]["tee"]
    assert v["verdict"] == "pass" and v["sigma_norm"] > 0.0


def test_the_gas_constant_and_gamma_gradients_match_a_difference() -> None:
    """cb.specific_gas_constant_dY and cb.isentropic_expansion_coefficient_dT_dY
    against central differences in the mass fractions (renormalised, as a
    node's X is) and in T."""
    Y = np.asarray(_co2_air(0.07), dtype=float)
    Y[cb.species.species_index("CH4")] += 0.02
    Y /= Y.sum()
    T = 700.0

    def X(Yv: np.ndarray) -> list[float]:
        return list(cb.mass_to_mole(list(Yv / Yv.sum())))

    dR = np.asarray(cb.specific_gas_constant_dY(X(Y)))
    dg_dT, dg = cb.isentropic_expansion_coefficient_dT_dY(T, X(Y))
    h = 1e-6
    for k in range(Y.size):
        Yp, Ym = Y.copy(), Y.copy()
        Yp[k] += h
        Ym[k] -= h
        if Ym[k] < 0.0:
            continue
        fd_R = (cb.specific_gas_constant(X(Yp)) - cb.specific_gas_constant(X(Ym))) / (2 * h)
        g = cb.isentropic_expansion_coefficient
        fd_g = (g(T, X(Yp)) - g(T, X(Ym))) / (2 * h)
        assert math.isclose(dR[k], fd_R, rel_tol=1e-6, abs_tol=1e-6), k
        assert math.isclose(dg[k], fd_g, rel_tol=1e-6, abs_tol=1e-9), k
    g = cb.isentropic_expansion_coefficient
    fd_T = (g(T + 1e-3, X(Y)) - g(T - 1e-3, X(Y))) / 2e-3
    assert math.isclose(dg_dT, fd_T, rel_tol=1e-6)


def _stress_legacy_tee(seed: int) -> FlowNetwork:
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
    from networks import FAMILIES

    return FAMILIES["legacy_tee"](seed)[0]


def test_a_failed_forward_retry_leaves_the_refused_solve_on_the_solver() -> None:
    """Stress legacy_tee seed 30: the cold solve's root is refused (a leg
    against the tee_type) and the forward retry does not converge. The solver's
    last iterate and diagnostics are the returned solve's, not the retry's
    (the #1233 review: 3.3-4.6 relative off on four seeds)."""
    s = NetworkSolver(_stress_legacy_tee(30))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=20)
    assert r["__outcome__"] == "inconsistent"
    assert "did not converge" in r["__message__"]
    x = np.asarray(r["__x_solution__"], dtype=float)
    assert np.array_equal(np.asarray(s._last_x, dtype=float), x)
    assert np.array_equal(np.asarray(s._diagnostic_data["solution"], dtype=float), x)


def test_the_forward_retry_honours_auto_retry() -> None:
    s = NetworkSolver(_stress_legacy_tee(30))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=20, auto_retry=False)
    assert "forward start" not in r["__message__"]
