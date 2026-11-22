"""Tests for MomentumChamberNode proper implementation."""

import pytest

import combaero as cb
from combaero.network import (
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)


def test_momentum_chamber_dynamic_pressure():
    """Verify MomentumChamberNode enforces Pt = P + 0.5*rho*v^2."""
    graph = FlowNetwork()

    Y_air = cb.species.dry_air_mass()
    inlet = MassFlowBoundary("inlet", m_dot=0.2, Tt=600.0, Y=Y_air)
    chamber = MomentumChamberNode(
        "chamber", area=0.02
    )  # 0.02 m^2, larger area for easier convergence
    outlet = PressureBoundary("outlet", Pt=101325.0)

    orf1 = OrificeElement(
        "orf1", "inlet", "chamber", Cd=0.6, diameter=0.035682, correlation="fixed"
    )
    orf2 = OrificeElement(
        "orf2", "chamber", "outlet", Cd=0.6, diameter=0.035682, correlation="fixed"
    )

    graph.add_node(inlet)
    graph.add_node(chamber)
    graph.add_node(outlet)
    graph.add_element(orf1)
    graph.add_element(orf2)

    solver = NetworkSolver(graph)
    sol = solver.solve()

    assert sol["__success__"]

    # Get chamber state
    P = sol["chamber.P"]
    Pt = sol["chamber.Pt"]
    m_dot = 0.2  # From inlet
    T = solver._derived_states["chamber"][0]
    Y = solver._derived_states["chamber"][1]

    # Verify dynamic pressure relation: Pt = P + 0.5*rho*v^2
    X = cb.mass_to_mole(cb.normalize_fractions(Y))
    rho = cb.density(T, P, X)
    v = m_dot / (rho * chamber.area)
    q_dynamic = 0.5 * rho * v * v

    # Check that Pt = P + q_dynamic
    assert Pt == pytest.approx(P + q_dynamic, rel=1e-4)

    # Verify dynamic pressure is non-zero (unlike plenum)
    assert q_dynamic > 1.0  # Should be several Pa

    print("Momentum chamber test passed:")
    print(f"  P = {P:.1f} Pa")
    print(f"  Pt = {Pt:.1f} Pa")
    print(f"  q_dynamic = {q_dynamic:.1f} Pa")
    print(f"  v = {v:.2f} m/s")


def test_momentum_chamber_vs_plenum():
    """Compare MomentumChamberNode to PlenumNode - should have different Pt."""
    from combaero.network import PlenumNode

    Y_air = cb.species.dry_air_mass()
    m_dot = 1.0  # Higher flow for larger dynamic pressure

    # Test with momentum chamber
    graph_mom = FlowNetwork()
    inlet_mom = MassFlowBoundary("inlet", m_dot=m_dot, Tt=600.0, Y=Y_air)
    chamber_mom = MomentumChamberNode("chamber", area=0.01)
    outlet_mom = PressureBoundary("outlet", Pt=101325.0)

    graph_mom.add_node(inlet_mom)
    graph_mom.add_node(chamber_mom)
    graph_mom.add_node(outlet_mom)
    graph_mom.add_element(
        OrificeElement("orf1", "inlet", "chamber", Cd=0.6, diameter=0.035682, correlation="fixed")
    )
    graph_mom.add_element(
        OrificeElement("orf2", "chamber", "outlet", Cd=0.6, diameter=0.035682, correlation="fixed")
    )

    solver_mom = NetworkSolver(graph_mom)
    sol_mom = solver_mom.solve()

    # Test with plenum (Pt = P)
    graph_ple = FlowNetwork()
    inlet_ple = MassFlowBoundary("inlet", m_dot=m_dot, Tt=600.0, Y=Y_air)
    chamber_ple = PlenumNode("chamber")
    outlet_ple = PressureBoundary("outlet", Pt=101325.0)

    graph_ple.add_node(inlet_ple)
    graph_ple.add_node(chamber_ple)
    graph_ple.add_node(outlet_ple)
    graph_ple.add_element(
        OrificeElement("orf1", "inlet", "chamber", Cd=0.6, diameter=0.035682, correlation="fixed")
    )
    graph_ple.add_element(
        OrificeElement("orf2", "chamber", "outlet", Cd=0.6, diameter=0.035682, correlation="fixed")
    )

    solver_ple = NetworkSolver(graph_ple)
    sol_ple = solver_ple.solve()

    assert sol_mom["__success__"]
    assert sol_ple["__success__"]

    # Plenum: Pt = P
    P_ple = sol_ple["chamber.P"]
    Pt_ple = sol_ple["chamber.Pt"]
    assert Pt_ple == pytest.approx(P_ple, rel=1e-6)

    # Momentum chamber: Pt > P (due to dynamic pressure)
    P_mom = sol_mom["chamber.P"]
    Pt_mom = sol_mom["chamber.Pt"]
    assert Pt_mom > P_mom

    # Dynamic pressure should be significant
    q_dynamic = Pt_mom - P_mom
    assert q_dynamic > 10.0  # At least 10 Pa

    print("Plenum vs Momentum Chamber:")
    print(f"  Plenum: P={P_ple:.1f} Pa, Pt={Pt_ple:.1f} Pa, Delta_P={Pt_ple - P_ple:.1f} Pa")
    print(f"  Momentum: P={P_mom:.1f} Pa, Pt={Pt_mom:.1f} Pa, Delta_P={q_dynamic:.1f} Pa")


def test_momentum_chamber_with_heat_exchange():
    """Verify MomentumChamberNode supports energy boundaries (heat exchange)."""
    from combaero.network import EnergyBoundary

    graph = FlowNetwork()

    inlet = MassFlowBoundary("inlet", m_dot=0.5, Tt=600.0)
    inlet.Y = cb.species.dry_air_mass()
    chamber = MomentumChamberNode("chamber", area=0.01)
    outlet = PressureBoundary("outlet")
    outlet.Pt = 101325.0
    outlet.Tt = 300.0
    outlet.Y = cb.species.dry_air_mass()

    # Add heat exchange to chamber
    chamber.add_energy_boundary(EnergyBoundary("heater", Q=5000.0))

    orf1 = OrificeElement(
        "orf1", "inlet", "chamber", Cd=0.6, diameter=0.035682, correlation="fixed"
    )
    orf2 = OrificeElement(
        "orf2", "chamber", "outlet", Cd=0.6, diameter=0.035682, correlation="fixed"
    )

    graph.add_node(inlet)
    graph.add_node(chamber)
    graph.add_node(outlet)
    graph.add_element(orf1)
    graph.add_element(orf2)

    solver = NetworkSolver(graph)
    sol = solver.solve()

    assert sol["__success__"]

    # Verify temperature increased due to heat addition
    T_chamber = solver._derived_states["chamber"][0]
    T_inlet = 600.0
    assert T_chamber > T_inlet

    # Verify dynamic pressure relation still holds
    P = sol["chamber.P"]
    Pt = sol["chamber.Pt"]
    m_dot = 0.5
    Y = solver._derived_states["chamber"][1]

    X = cb.mass_to_mole(cb.normalize_fractions(Y))
    rho = cb.density(T_chamber, P, X)
    v = m_dot / (rho * chamber.area)
    q_dynamic = 0.5 * rho * v * v

    assert Pt == pytest.approx(P + q_dynamic, rel=1e-4)

    print("Momentum chamber with heat exchange:")
    print(
        f"  T_in = {T_inlet:.1f} K -> T_chamber = {T_chamber:.1f} K (Delta_T = {T_chamber - T_inlet:.1f} K)"
    )
    print(f"  Pt = P + q_dynamic: {Pt:.1f} = {P:.1f} + {q_dynamic:.1f} Pa")


# ---------------------------------------------------------------------------
# The node's T is the STAGNATION temperature (#576)
# ---------------------------------------------------------------------------

_X_AIR = list(cb.species.dry_air())
_Y_AIR = list(cb.mole_to_mass(_X_AIR))


def _closure(P: float, Pt: float, m: float, T0: float, Y: list[float], A: float):
    from combaero import _solver_tools

    return _solver_tools.momentum_chamber_residual_and_jacobian(P, Pt, m, T0, Y, A)


@pytest.mark.parametrize("M", [0.3, 0.8, 1.0, 1.5])
def test_closure_holds_at_an_exact_stagnation_state(M: float) -> None:
    """Build a flow at Mach M from (Pt, T0) with the library's own isentropic
    relations; the closure, fed T0 as the node carries it, must close on it.
    With T0 read as static (before #576) the residual was -0.1% of Pt at M 0.3
    and -12% at M 1.0."""
    T0, Pt, A = 600.0, 3.0e5, 2.0e-3
    T = cb.T_from_stagnation(T0, M, _X_AIR)
    P = cb.P_from_stagnation(Pt, T0, M, _X_AIR)
    m = cb.density(T, P, _X_AIR) * M * cb.speed_of_sound(T, _X_AIR) * A
    r = _closure(P, Pt, m, T0, _Y_AIR, A)
    assert r.T_static == pytest.approx(T, rel=1e-9)
    mach = r.M
    assert mach == pytest.approx(M, rel=1e-8)
    assert r.rho_static == pytest.approx(cb.density(T, P, _X_AIR), rel=1e-9)
    if M <= 1.0:
        assert r.residual == pytest.approx(0.0, abs=1e-9 * Pt)
        # The dynamic head, returned without cancellation against P.
        assert r.q_dyn == pytest.approx(Pt - r.residual - P, rel=1e-12)
        return
    # Supersonic: no root (#788). The closure is continued past the sonic
    # point by its mirror image 2 Pt* - Pt_cl(P), and at this state Pt_cl = Pt,
    # so the residual is Pt - (2 Pt* - Pt) = 2 (Pt - Pt*) > 0. The static
    # state fields stay the closure's at P.
    assert r.bounded
    # The sonic state of the same flux G = m / A and T0: P* = G R T* / a*,
    # on the isentrope Pt* / P* = Pt / P(M = 1).
    T_s = cb.T_from_stagnation(T0, 1.0, _X_AIR)
    R = cb.specific_gas_constant(_X_AIR)
    P_star = (m / A) * R * T_s / cb.speed_of_sound(T_s, _X_AIR)
    Pt_star = P_star * Pt / cb.P_from_stagnation(Pt, T0, 1.0, _X_AIR)
    assert r.P_sonic == pytest.approx(P_star, rel=1e-8)
    assert r.Pt_sonic == pytest.approx(Pt_star, rel=1e-8)
    assert r.Pt_sonic < Pt
    assert r.residual == pytest.approx(2.0 * (Pt - r.Pt_sonic), rel=1e-8)
    assert r.q_dyn == pytest.approx(Pt - P, rel=1e-9)


def test_closure_passes_the_true_choked_flux() -> None:
    """At fixed (Pt, T0) the largest mass flux the closure admits subsonically, over the
    static pressure, is the choked flux rho* a* A. T0 read as static capped
    it at 92.0% (peak at M 0.913), and Newton stalled on any chamber asked to
    carry more -- an ejector outlet pipe of about A_3 (#567)."""
    from scipy.optimize import brentq

    T0, Pt, A = 300.0, 2.0e5, 1.0e-3

    T_star = cb.T_from_stagnation(T0, 1.0, _X_AIR)
    P_star = cb.P_from_stagnation(Pt, T0, 1.0, _X_AIR)
    m_star = cb.density(T_star, P_star, _X_AIR) * cb.speed_of_sound(T_star, _X_AIR) * A

    def m_at(P: float) -> float:
        # Bracketed below the flux where the continued closure turns back
        # (#788: past the sonic point it is 2 Pt* - Pt_cl, not monotone in m
        # at fixed P).
        f = lambda m: _closure(P, Pt, m, T0, _Y_AIR, A).residual  # noqa: E731
        return brentq(f, 0.0, 1.1 * m_star, xtol=1e-13)

    peak = max(m_at(P_star * s) for s in (1.0, 1.01, 1.03))
    assert peak == pytest.approx(m_star, rel=1e-6)
    assert m_at(P_star) == pytest.approx(m_star, rel=1e-9)
    # Below P* the closure is continued past the sonic point (#788): its root
    # there carries MORE than the choked flux, and says so (bounded) -- the
    # flux a choked duct forces into a chamber that cannot carry it at Pt.
    for s in (0.97, 0.99):
        f = lambda m, s=s: _closure(P_star * s, Pt, m, T0, _Y_AIR, A).residual  # noqa: E731
        m = brentq(f, 0.5 * m_star, 1.1 * m_star, xtol=1e-13)
        assert m > m_star
        assert _closure(P_star * s, Pt, m, T0, _Y_AIR, A).bounded


@pytest.mark.parametrize("M", [0.05, 0.5, 0.95, 1.4])
@pytest.mark.parametrize("sign", [1.0, -1.0])
def test_closure_jacobian_matches_central_differences(M: float, sign: float) -> None:
    """Residual and static-density sensitivities in P, m_dot, T0 and Y, on a
    CO2 / air mix (Y moves R, h and s)."""
    T0, Pt, A = 900.0, 3.0e5, 2.0e-3
    X = [0.0] * len(_X_AIR)
    for k, x in enumerate(_X_AIR):
        X[k] = 0.9 * x
    X[cb.species_index_from_name("CO2")] += 0.1
    Y = list(cb.mole_to_mass(X))
    T = cb.T_from_stagnation(T0, M, X)
    P = cb.P_from_stagnation(Pt, T0, M, X) * 0.98
    m = sign * cb.density(T, P, X) * M * cb.speed_of_sound(T, X) * A
    r = _closure(P, Pt, m, T0, Y, A)

    def cd(f, x: float, h: float) -> float:
        return (f(x + h) - f(x - h)) / (2.0 * h)

    def both(**kw):
        a = {"P": P, "m": m, "T0": T0, "Y": Y} | kw
        rr = _closure(a["P"], Pt, a["m"], a["T0"], a["Y"], A)
        return rr.residual, rr.rho_static

    for k, (an_r, an_rho, f, x, h) in {
        "P": (r.d_res_dP, r.drho_dP, lambda v: both(P=v), P, P * 1e-6),
        "m": (r.d_res_dmdot, r.drho_dmdot, lambda v: both(m=v), m, abs(m) * 1e-6),
        "T0": (r.d_res_dT, r.drho_dT, lambda v: both(T0=v), T0, T0 * 1e-6),
    }.items():
        fd_r = cd(lambda v, f=f: f(v)[0], x, h)
        fd_rho = cd(lambda v, f=f: f(v)[1], x, h)
        assert an_r == pytest.approx(fd_r, rel=2e-5, abs=1e-9 * Pt / abs(x)), k
        assert an_rho == pytest.approx(fd_rho, rel=2e-5, abs=1e-12), k

    k_co2 = cb.species_index_from_name("CO2")
    h_Y = 1e-4
    Yp, Ym = list(Y), list(Y)
    Yp[k_co2] += h_Y
    Ym[k_co2] -= h_Y
    fd_r = (both(Y=Yp)[0] - both(Y=Ym)[0]) / (2 * h_Y)
    fd_rho = (both(Y=Yp)[1] - both(Y=Ym)[1]) / (2 * h_Y)
    assert r.d_res_dY[k_co2] == pytest.approx(fd_r, rel=1e-6)
    assert r.drho_dY[k_co2] == pytest.approx(fd_rho, rel=1e-6)


def _chamber_closure(
    P: float, Pt: float, m: float, T0: float, Y: list[float], A: float, *, incompressible=False
):
    from combaero import _solver_tools

    return _solver_tools.momentum_chamber_residual_and_jacobian(P, Pt, m, T0, Y, A, incompressible)


def test_static_state_below_the_data_has_an_exact_pressure_partial_1255() -> None:
    """#1255: a chamber whose static T (168 K) is below the thermo data's
    200 K and whose T0 (317 K) is above it. One Gauss-Legendre rule across
    the low tail's joint (C2: the third derivative of cp jumps) missed the
    enthalpy and entropy integrals by enough that d(res)/dP, which uses the
    integrands, was 2.2e-6 off the residual's own slope (the ejector's
    suction chamber with its leg reversed). The quadrature now splits there."""
    A = 1.9634954084936207e-3
    P, m, T0 = 28857.756662553693, 0.6445370338698646, 317.2289549347265
    r = _chamber_closure(P, P, m, T0, _Y_AIR, A)
    assert r.T_static < 200.0 < T0
    for h in (1e-2, 1e-1):
        fd = (
            _chamber_closure(P + h, P, m, T0, _Y_AIR, A).residual
            - _chamber_closure(P - h, P, m, T0, _Y_AIR, A).residual
        ) / (2 * h)
        assert r.d_res_dP == pytest.approx(fd, rel=2e-8)


def test_static_state_across_a_c2_bridge_edge_has_an_exact_pressure_partial() -> None:
    """The same quadrature flaw at a C2 bridge's window edge (950 K, where the
    re-fit's third derivative of cp jumps): T_static 935 K, T0 1100 K. One
    rule across it left d(res)/dP 4.7e-8 off the residual's own slope; the
    quadrature splits at the active thermo model's breakpoints."""
    A = 1.9635e-3
    P, m, T0 = 2.0e5, 0.9, 1100.0
    r = _chamber_closure(P, P, m, T0, _Y_AIR, A)
    assert r.T_static < 950.0
    h = 1.0
    fd = (
        _chamber_closure(P + h, P, m, T0, _Y_AIR, A).residual
        - _chamber_closure(P - h, P, m, T0, _Y_AIR, A).residual
    ) / (2 * h)
    assert r.d_res_dP == pytest.approx(fd, rel=5e-9)


@pytest.mark.filterwarnings("ignore:.*STATIC temperature:FutureWarning")
def test_node_reports_the_static_state_its_closure_used() -> None:
    """mach, rho and T_static of a solved chamber come from the static state,
    not from its stagnation T."""
    graph = FlowNetwork()
    graph.add_node(MassFlowBoundary("inlet", m_dot=0.5, Tt=600.0, Y=cb.species.dry_air_mass()))
    chamber = MomentumChamberNode("chamber", area=3e-3, regime="compressible")
    graph.add_node(chamber)
    graph.add_node(PressureBoundary("outlet", Pt=101325.0))
    # Incompressible orifices put the chamber above Mach 0.2, where the
    # static and stagnation states are told apart; compressible ones pass
    # the same flow at a higher chamber pressure, Mach 0.15.
    for eid, a, b in (("orf1", "inlet", "chamber"), ("orf2", "chamber", "outlet")):
        graph.add_element(
            OrificeElement(
                eid, a, b, Cd=0.6, diameter=0.04, correlation="fixed", regime="incompressible"
            )
        )
    solver = NetworkSolver(graph)
    sol = solver.solve()
    assert sol["__success__"]
    T0 = sol["chamber.Tt"]
    P = sol["chamber.P"]
    T_static = sol["chamber.T"]
    v = 0.5 / (sol["chamber.rho"] * chamber.area)
    # Energy: h(T0) = h(T) + v^2 / 2, at the reported static state.
    assert cb.h_mass(T0, _X_AIR) - cb.h_mass(T_static, _X_AIR) == pytest.approx(
        0.5 * v * v, rel=1e-8
    )
    assert sol["chamber.rho"] == pytest.approx(cb.density(T_static, P, _X_AIR), rel=1e-9)
    assert sol["chamber.mach"] == pytest.approx(v / cb.speed_of_sound(T_static, _X_AIR), rel=1e-8)
    assert sol["chamber.mach"] > 0.2


def test_an_incompressible_chamber_reports_its_low_mach_state() -> None:
    """The incompressible closure (#588): the static state IS the stagnation
    one and Pt = P + v^2 rho / 2 at rho(P, T0) -- no static-T drop, no sonic
    limit. The node reports exactly the state its residual used."""
    graph = FlowNetwork()
    graph.add_node(MassFlowBoundary("inlet", m_dot=0.5, Tt=600.0, Y=cb.species.dry_air_mass()))
    chamber = MomentumChamberNode("chamber", area=3e-3, regime="incompressible")
    graph.add_node(chamber)
    graph.add_node(PressureBoundary("outlet", Pt=101325.0))
    graph.add_element(
        OrificeElement("orf1", "inlet", "chamber", Cd=0.6, diameter=0.04, correlation="fixed")
    )
    graph.add_element(
        OrificeElement("orf2", "chamber", "outlet", Cd=0.6, diameter=0.04, correlation="fixed")
    )
    sol = NetworkSolver(graph).solve()
    assert sol["__success__"]
    T0, P, Pt = sol["chamber.Tt"], sol["chamber.P"], sol["chamber.Pt"]
    rho = float(cb.density(T0, P, _X_AIR))
    v = 0.5 / (rho * chamber.area)
    assert sol["chamber.T"] == pytest.approx(T0, rel=1e-12)
    assert sol["chamber.rho"] == pytest.approx(rho, rel=1e-9)
    assert Pt - P == pytest.approx(0.5 * rho * v * v, rel=1e-9)
