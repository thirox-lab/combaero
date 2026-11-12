"""The state relay of a single-source, multi-sink element (branching junctions).

A node's mixed temperature depends on the flows INTO it when it carries heat
(EnergyBoundary): T = T_in + Q/(m cp). For a single-source element the stream
into a sink node is ``flow_at_node(sink)``, so its derivative must be taken
at the sink too. The relay took it at the SOURCE node, so a branching tee's
straight sink saw d/d(m_com) only and lost the -1 on m_branch: that Jacobian
entry was 0 against a finite difference of -1.9.
"""

from __future__ import annotations

import math

import numpy as np
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    EnergyBoundary,
    FlowNetwork,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    TeeJunctionElement,
)

Y = cb.mole_to_mass(cb.species.dry_air())


def _max_rel_jacobian_error(net: FlowNetwork) -> float:
    s = NetworkSolver(net)
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    x = np.array(r["__x_solution__"], dtype=float)
    _, jac = s._residuals_and_jacobian(x)
    # Step 1e-6: at 1e-7 the smallest entries (~6e-3) are round-off bound --
    # the FD drifts from J as the step shrinks (4.7e-6 at 1e-5, 2.9e-4 at 1e-7
    # locally, 1.1e-3 on Linux, 1e-2 at 1e-8). At 1e-6 the clean case measures
    # 6e-5 and a 0.5% error in its smallest entry still reads 5e-3 (#612).
    fd = fd_jacobian(s, x, floor=1e-10)
    mask = np.abs(fd) > 1e-6
    return float((np.abs(jac.toarray() - fd)[mask] / np.abs(fd[mask])).max())


def test_branching_tee_into_a_heated_node() -> None:
    net = FlowNetwork()
    net.add_node(PressureBoundary("pb_common", Pt=2.1e5, Tt=400.0, Y=Y))
    net.add_node(PlenumNode("s"))
    net.add_node(PressureBoundary("pb_out", Pt=2.06e5, Tt=400.0, Y=Y))
    net.add_node(PressureBoundary("pb_branch", Pt=2.06e5, Tt=400.0, Y=Y))
    net.nodes["s"].add_energy_boundary(EnergyBoundary("q", Q=5000.0))
    net.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="pb_common",
            straight_node="s",
            branch_node="pb_branch",
            theta=math.pi / 2,
            F_C=0.01,
            psi=1.0,
            tee_type="branching",
        )
    )
    # The straight pipe is F_C's own: an explicit F_C that contradicts the
    # pipe at its arm is refused (#1023); 0.1128 m was F_C to 3 digits.
    D_C = math.sqrt(4.0 * 0.01 / math.pi)
    net.add_element(ChannelElement("c", "s", "pb_out", length=0.5, diameter=D_C))
    assert _max_rel_jacobian_error(net) < 1e-3
