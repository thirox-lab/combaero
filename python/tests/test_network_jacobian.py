import numpy as np
from _fd_jacobian import fd_jacobian

from combaero.network import (
    ChannelElement,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)


def _assert_jacobian_matches_fd(solver: NetworkSolver, max_rel_tol: float) -> None:
    x0 = solver._build_x0()

    _, jac_analytical = solver._residuals_and_jacobian(x0)
    jac_analytical = jac_analytical.toarray()

    jac_numerical = fd_jacobian(solver, x0, rel=1e-7, floor=1e-10)

    diff = np.abs(jac_analytical - jac_numerical)
    mask = np.abs(jac_numerical) > 1e-5
    if np.any(mask):
        rel_diff = diff[mask] / np.abs(jac_numerical[mask])
        max_rel_diff = np.max(rel_diff)
        assert max_rel_diff < max_rel_tol, f"Jacobian mismatch! Max relative diff: {max_rel_diff}"
    else:
        max_diff = np.max(diff)
        assert max_diff < 1e-5


def test_jacobian_accuracy():
    """Global Jacobian FD check for a representative incompressible network."""
    net = FlowNetwork()

    # Nodes
    n_in = PressureBoundary("in", Pt=200000.0, Tt=300.0)
    n1 = PlenumNode("n1")
    n_out = PressureBoundary("out", Pt=101325.0, Tt=300.0)

    net.add_node(n_in)
    net.add_node(n1)
    net.add_node(n_out)

    # Elements. The bore sits in the 0.1 m channel's line (beta 0.5 against
    # it, its reverse approach pipe); it was 0.112838 m, larger than that
    # pipe, a geometry a solve refuses and an unresolved fresh solver did not
    # check (#1148).
    orf = OrificeElement(
        "orf", from_node="in", to_node="n1", Cd=0.6, diameter=0.05, correlation="fixed"
    )
    channel = ChannelElement(
        "channel", from_node="n1", to_node="out", length=1.0, diameter=0.1, roughness=1e-5
    )

    net.add_element(orf)
    net.add_element(channel)

    solver = NetworkSolver(net)
    _assert_jacobian_matches_fd(solver, max_rel_tol=0.01)


def test_jacobian_accuracy_compressible_network() -> None:
    """Global Jacobian FD check for a compressible channel-orifice network."""
    net = FlowNetwork()

    n_in = PressureBoundary("in", Pt=4.8e5, Tt=540.0)
    n_mid = PlenumNode("mid")
    n_out = PressureBoundary("out", Pt=1.9e5, Tt=520.0)

    net.add_node(n_in)
    net.add_node(n_mid)
    net.add_node(n_out)

    net.add_element(
        ChannelElement(
            "channel",
            from_node="in",
            to_node="mid",
            length=3.2,
            diameter=0.022,
            roughness=8.0e-5,
            regime="compressible",
            friction_model="haaland",
        )
    )
    net.add_element(
        OrificeElement(
            "orf",
            from_node="mid",
            to_node="out",
            Cd=0.69,
            diameter=0.012361,
            regime="compressible",
            correlation="fixed",
        )
    )

    solver = NetworkSolver(net)
    _assert_jacobian_matches_fd(solver, max_rel_tol=0.05)


if __name__ == "__main__":
    test_jacobian_accuracy()
