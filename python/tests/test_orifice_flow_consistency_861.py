"""The own-flow consistency solve of an orifice in a pipe returns a root (#861).

OrificeElement solves m* = m_calc(m*, x) itself when its drive reads its own
flow through the approach pipe (#651), and the row is m - m*. Where
m_calc(m) carries a kink, plain Newton 2-cycles across it to the iteration
cap and the m* it stopped at is no root. Found on the stress family
bleed_station (seeds 6 and 39) with compressible_model="jobson1955" at
beta 0.84, past the model's stated beta <= 0.75, where Jobson's held
discriminant puts the kink in m_calc: the row was 6% off and its Pt column
matched neither side, and the solve needed 685 LM evaluations, not 437.
The state below is an LM iterate from that run.
"""

from __future__ import annotations

import copy
import warnings

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    CrossflowSegmentElement,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    components,
)

Y_AIR = list(cb.species.dry_air_mass())
# b.P, b.Pt, s0.m_dot, s1.m_dot, lat_o.m_dot
X_ITERATE = {
    "b.P": 499711.50207837514,
    "b.Pt": 499711.0416041446,
    "s0.m_dot": 0.13523234095077175,
    "s1.m_dot": 0.11229202713690327,
    "lat_o.m_dot": 0.03342006698406128,
}


def _bleed_station(model: str) -> FlowNetwork:
    """scripts/stress/networks.py bleed_station(6), written out."""
    g = FlowNetwork()
    for nid, Pt, Tt in (
        ("in", 501000.0, 300.0),
        ("out", 500000.0, 300.0),
        ("lat", 499000.0, 600.0),
    ):
        g.add_node(PressureBoundary(nid, Pt=Pt, Tt=Tt, Y=Y_AIR))
    g.add_node(PlenumNode("b"))
    kappa = cb.STATION_KAPPA_BLEED_BASSETT
    g.add_element(
        CrossflowSegmentElement(
            "s0",
            "in",
            "b",
            length=1e-6,
            area=1e-3,
            entry_K=0.0,
            to_kappa=kappa,
            next_seg="s1",
            regime="compressible",
        )
    )
    g.add_element(
        CrossflowSegmentElement(
            "s1",
            "b",
            "out",
            length=1e-6,
            area=1e-3,
            from_kappa=kappa,
            prev_seg="s0",
            regime="compressible",
        )
    )
    g.add_element(
        OrificeElement(
            "lat_o",
            "b",
            "lat",
            diameter=0.03,
            Cd=0.6,
            correlation="fixed",
            regime="compressible",
            compressible_model=model,
        )
    )
    return g


def _solver_at_iterate(model: str) -> tuple[NetworkSolver, np.ndarray]:
    s = NetworkSolver(_bleed_station(model))
    s.network.resolve_all_topology()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s._build_x0()  # the unknowns and their order; no solve
    x = np.array([X_ITERATE[name] for name in s.unknown_names])
    return s, x


@pytest.mark.parametrize("model", ["nozzle", "jobson1955"])
def test_the_row_is_m_minus_a_root(model: str) -> None:
    s, x = _solver_at_iterate(model)
    el = s.network.elements["lat_o"]
    assert el._approach_reads_flow()  # the path under test
    seen = []
    orig = components.OrificeElement._flow_consistent_residuals

    def spy(self, si, so, flows, cd):
        res, jac = orig(self, si, so, flows, cd)
        seen.append((copy.copy(si), copy.copy(so), flows, cd, res[0]))
        return res, jac

    components.OrificeElement._flow_consistent_residuals = spy
    try:
        s._residuals_and_jacobian(x)
    finally:
        components.OrificeElement._flow_consistent_residuals = orig
    si, so, flows, cd, res = seen[-1]
    m_star = float(si.m_dot) - res
    si.m_dot = so.m_dot = m_star
    m_calc = el._flow_and_partials(si, so, cd)[0]
    assert m_calc == pytest.approx(m_star, rel=1e-9)


@pytest.mark.parametrize("model", ["nozzle", "jobson1955"])
def test_the_rows_jacobian_matches_central_differences(model: str) -> None:
    s, x = _solver_at_iterate(model)
    kinds = s._row_kinds(x)
    row = next(i for i, (owner, _) in enumerate(kinds) if "lat_o" in owner)
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray() if hasattr(J, "toarray") else np.asarray(J)
    fd_row = fd_jacobian(s, x)[row]
    for j, name in enumerate(s.unknown_names):
        fd = fd_row[j]
        scale = max(abs(fd), 1e-12) if fd != 0.0 else 1.0
        assert abs(J[row, j] - fd) <= 1e-5 * scale + 1e-14, (name, J[row, j], fd)
