"""Rogers & Hersh's cap band carries the hole's own Cd slope (#1248).

The impingement crossflow law is joined C1 to the hole's own Cd0 over
+/-ROGERS_HERSH_1976_CAP_BAND (C++ ``rogers_hersh_1976_cd``): inside the
band Cd = Cd0 f(cd / Cd0), so it reads Cd0 -- for an IdelchikThick hole its
Reynolds number, and so the row's own flow. ``hole_cd_compose`` zeroed the
hole's partials wherever the law acted, so in the band the own-flow slope
dCd/dm was dropped. The row solves m* = Cd(m*) K(x) (#582) and divides every
column by 1 - K dCd/dm: the whole row was ~0.8% off, the sweep's chamber
pressure column first (scaled 7.7e-3, Case("impingement_bypass", 0, 1, 0)
at its start vector). The law's dCd/dCd0 now scales the hole's partials.
"""

from __future__ import annotations

import pytest
from _jacobian_sweep import Case, point, resolved_solver

import combaero as cb

CASE = Case("impingement_bypass", 0, 1, 0)


def _row_inputs(s, el, x, dP: float = 0.0, m: float | None = None):
    """The plate's (state_in, state_out, flows) at x, the out static moved
    by dP and the own flow set to m."""
    s._propagate_states(x, with_relay=False)
    si = s._get_node_state(s.network.nodes[el.from_node], x)
    so = s._get_node_state(s.network.nodes[el.to_node], x)
    own = float(x[s._unknown_indices[el.id][0]]) if m is None else m
    si.m_dot = so.m_dot = own
    so.P += dP
    return si, so, s._network_flows(el, x)


@pytest.fixture(scope="module")
def solver_at_start():
    x0, _ = point(CASE, "start vector")
    s = resolved_solver(CASE)
    s._residuals_and_jacobian(x0)
    return s, x0


def test_the_cap_band_cd_carries_the_own_flow_slope(solver_at_start) -> None:
    """Row 1 sits in the band at the start vector (row 3 is capped there
    and enters it at the sweep's perturbed points)."""
    s, x0 = solver_at_start
    el = s.network.elements["ia__p1"]
    si, so, flows = _row_inputs(s, el, x0)
    r, _ = el._hole_cd._evaluate(si, so, flows, jet=True)
    # The point is inside the band: the law acts and reads the hole's Cd.
    assert r.jet_acts
    jc = el._hole_cd._jet_drive(si, so, flows)
    law = cb.rogers_hersh_1976_cd(
        float(r.cd_hole), jc.G_i, abs(jc.m_c) / el.channel_area, not jc.forward
    )
    assert 0.0 < law.dCd_dCd0 < 1.0
    m = float(si.m_dot)
    h = 1e-6 * abs(m)
    fd = (
        el._effective_Cd(*_row_inputs(s, el, x0, m=m + h))
        - el._effective_Cd(*_row_inputs(s, el, x0, m=m - h))
    ) / (2.0 * h)
    analytic = el._dCd_dnames(si, so, flows).get(f"{el.id}.m_dot", 0.0)
    assert abs(fd) > 1e-2  # the band's slope is real here
    assert analytic == pytest.approx(fd, rel=1e-5)


@pytest.mark.parametrize("row", ["ia__p1", "ia__p3"])
def test_the_rows_chamber_pressure_column_matches_fd(solver_at_start, row: str) -> None:
    s, x0 = solver_at_start
    el = s.network.elements[row]
    res, jac = el.residuals(*_row_inputs(s, el, x0))
    h = 1.0  # Pa, on ~3e5 Pa; the row is smooth in P here
    fp = el.residuals(*_row_inputs(s, el, x0, dP=h))[0][0]
    fm = el.residuals(*_row_inputs(s, el, x0, dP=-h))[0][0]
    fd = (fp - fm) / (2.0 * h)
    assert jac[0][f"{el.to_node}.P"] == pytest.approx(fd, rel=1e-6)
