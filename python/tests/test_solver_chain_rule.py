import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network.components import (
    CombustorNode,
    MassFlowBoundary,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)
from combaero.network.graph import FlowNetwork
from combaero.network.solver import NetworkSolver


def test_chain_rule_orifice_through_combustor():
    """
    Verify chain-rule relay: d(Residual_E2)/d(E1.m_dot)
    E1 (Fuel) -> Combustor -> E2 (Orifice)
    Changing E1.m_dot changes T_ad in Combustor, which changes E2.mdot_calc, affecting E2 residual.
    """
    net = FlowNetwork()

    # Air stream
    net.add_node(MassFlowBoundary("AIR_IN", m_dot=10.0, Tt=300.0))
    # Fuel stream (CH4)
    # Using mole fractions: 100% CH4
    X_fuel = np.zeros(cb.num_species())
    X_fuel[cb.species_index_from_name("CH4")] = (
        1.0  # CH4 (NASA order might differ, but mixer uses Y)
    )
    Y_fuel = list(cb.mole_to_mass(X_fuel))
    net.add_node(MassFlowBoundary("FUEL_IN", m_dot=0.5, Tt=300.0, Y=Y_fuel))

    net.add_node(CombustorNode("BURN", method="complete"))
    net.add_node(PressureBoundary("OUT", Pt=100000.0, Tt=300.0))

    net.add_element(
        OrificeElement("E_AIR", "AIR_IN", "BURN", Cd=0.6, diameter=0.356825, correlation="fixed")
    )
    net.add_element(
        OrificeElement("E_FUEL", "FUEL_IN", "BURN", Cd=0.6, diameter=0.112838, correlation="fixed")
    )
    net.add_element(
        OrificeElement("E_OUT", "BURN", "OUT", Cd=0.6, diameter=0.356825, correlation="fixed")
    )

    solver = NetworkSolver(net)

    # Solve to get a physical state
    sol = solver.solve()
    assert sol["__success__"]

    # Extract solution vector
    x_sol = np.array([sol[name] for name in solver.unknown_names])

    # Compute analytical and FD Jacobians
    res, jac_analytical_sparse = solver._residuals_and_jacobian(x_sol)
    jac_analytical = jac_analytical_sparse.toarray()
    jac_fd = fd_jacobian(solver, x_sol)

    # Compare
    max_err = np.max(np.abs(jac_analytical - jac_fd))
    print(f"Max Jacobian error: {max_err:.2e}")

    # Standard tolerance for FD vs Analytical comparison
    assert max_err < 1e-2, f"Jacobian mismatch: max error {max_err:.2e}"


def test_chain_rule_single_stream_passthrough():
    """
    Verify dT/dT_in = 1.0 through a series of nodes.
    """
    net = FlowNetwork()
    net.add_node(PressureBoundary("IN", Pt=200000.0, Tt=500.0))
    net.add_node(PlenumNode("N1"))
    net.add_node(PlenumNode("N2"))
    net.add_node(PressureBoundary("OUT", Pt=100000.0, Tt=300.0))

    net.add_element(
        OrificeElement("E1", "IN", "N1", Cd=0.6, diameter=0.112838, correlation="fixed")
    )
    net.add_element(
        OrificeElement("E2", "N1", "N2", Cd=0.6, diameter=0.112838, correlation="fixed")
    )
    net.add_element(
        OrificeElement("E3", "N2", "OUT", Cd=0.6, diameter=0.112838, correlation="fixed")
    )

    solver = NetworkSolver(net)
    sol = solver.solve()
    assert sol["__success__"]

    x_sol = np.array([sol[name] for name in solver.unknown_names])
    res, jac_analytical_sparse = solver._residuals_and_jacobian(x_sol)
    jac_analytical = jac_analytical_sparse.toarray()
    jac_fd = fd_jacobian(solver, x_sol)

    max_err = np.max(np.abs(jac_analytical - jac_fd))
    assert max_err < 1e-4


if __name__ == "__main__":
    pytest.main([__file__])
