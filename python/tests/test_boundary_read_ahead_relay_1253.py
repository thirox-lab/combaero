"""A flow boundary read ahead of the propagation order relays its own state (#1253).

The splitter's mass-flow inlet ``in`` feeds two outlets. With outlet o0's
feed reversed, ``in`` receives o0's gas and derives a mix of it, but it comes
last in the propagation order (the reversed legs close a cycle through the
outlet pressure boundary). Outlet o1 reads ``in`` before that: the residual
takes the boundary's own fixed T and Y (``_get_node_state_with_prev`` does not
lag a boundary), yet the relay borrowed the previous pass's MIX sensitivities
of ``in``. Channel a1's row then carried d/d(a0.m_dot) = 1.85e12 where the
residual has no dependence at all (FD exactly 0).
"""

from __future__ import annotations

import numpy as np
from _jacobian_sweep import Case, _dense, point, resolved_solver


def test_read_ahead_boundary_carries_no_borrowed_mix_sensitivity() -> None:
    case = Case("splitter", 0, 3, 2)
    _, pt = point(case, "o0 flows +-1e-09")
    s = resolved_solver(case)
    names = s.unknown_names
    j = names.index("a0.m_dot")
    F, J = s._residuals_and_jacobian(pt.x)
    J = _dense(J)
    F = np.asarray(F, dtype=float)
    # One-sided away from zero: a0.m_dot is -1e-9 and stays negative.
    h = 1e-11
    x = pt.x.copy()
    x[j] -= h
    F2, _ = s._residuals_and_jacobian(x)
    fd = (np.asarray(F2, dtype=float) - F) / (-h)
    # Scaled by |x_j| over the row's scale, as the sweep compares.
    err = np.abs(J[:, j] - fd) * abs(pt.x[j]) / np.maximum(np.abs(F), 1.0)
    assert err.max() < 1e-4, (names[int(np.argmax(err))], err.max())
    # Every row with no dependence in F (the sibling channel a1's among
    # them) has none in J either.
    assert np.all(J[fd == 0.0, j] == 0.0)
