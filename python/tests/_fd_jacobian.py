"""One finite-difference Jacobian for the network tests (#1107).

A central difference has a truncation error of about ``(h / s)^2``, where
``s`` is the pressure scale the residual curves on, not the pressure itself.
On ``m ~ sqrt(dp)`` that scale is the DROP: a step of ``1e-6 |P|`` is 3.1 Pa
at 31 bar, and on a 42.7 Pa duct it read as a 6e-4 "Jacobian error" that was
the step (#1107). The same rule as the sweep warm start: a step must be small
against dp, not against P (#1054).

So on a pressure column:

* the step is capped at ``DROP_FRACTION`` of the smallest gap between its
  value and any other pressure in the network (pressure unknowns and
  boundary Pt, its own node's skipped): a lower bound on the drops of the
  elements it touches, with no topology needed;
* never below ``P_STEP_FLOOR |P|``, where round-off takes over (residuals
  with an inner solve repeat only to ~1e-12);
* and the difference is Richardson-extrapolated, ``(4 D(h) - D(2h)) / 3``,
  which cancels the h^2 term, so the cap can be wide and a smooth network
  keeps the caller's step: a tight check (1e-7) read 1.5e-7 of round-off at a
  fixed 1e-7 |P|, 1.6e-8 at its own 1e-6.

On #1107's network the Jacobian reads 0.4-2.6e-6 against 6e-4 at the old
step. A drop the unknowns do not show is not covered: rows whose drive sits
near the 1 Pa flow-drive floor (``FANNO_FLOW_DRIVE_FLOOR``) curve on the Pa
scale between unknowns 70 Pa apart: a crossflow segment's static-pressure
column reads 5e-2 at 0.07 Pa and 1e-5 at 0.007 Pa, falling as the step
squared (test_liner_film_695; its tests carry their own tolerance).
The cap uses the gaps between PRESSURE UNKNOWNS, which are not the true drive
of a duct-fed hole with a recovered fraction f > 0 (it runs from P1 + f q1,
a derived pressure): a flag within ~1 Pa of that hole's zero flow is a step
artefact; refining the step converges as h^2 (to 1e-6 on the C2 rows).

The other columns keep a plain central difference at ``max(rel |x_j|,
floor)``. Callers pass ``rel`` per column where the residual's round-off
needs it (a warm-started Fanno iteration repeats only to ~1e-12, so flows
want 1e-5 there).

Compare in scaled variables, ``J_ij |x_j|`` (RESIDUAL-ROW-SCALING, see
``NetworkSolver._build_residual_scales``): ``scaled_error``.
"""

from __future__ import annotations

import numpy as np

from combaero.network import NetworkSolver, PressureBoundary

#: A pressure step as a fraction of the smallest pressure gap. With the
#: extrapolation the error is ~ (h / dp)^4.
DROP_FRACTION = 1e-2
#: Gaps below this fraction of |P| count as the same level (lossless links).
SAME_LEVEL = 1e-9
#: The cap never takes a pressure step below this fraction of |P|: a
#: film-cooling row capped at 1.5e-10 |P| read 3.5e-6 of round-off, against
#: 2e-7 at 1.5e-9 |P|.
P_STEP_FLOOR = 1e-8


def is_pressure(name: str) -> bool:
    """Pressure unknowns are ``.P``, ``.Pt``, ``.P_*`` or ``.Pt_*``."""
    last = name.rsplit(".", 1)[-1]
    return last in ("P", "Pt") or last.startswith(("P_", "Pt_"))


def fd_steps(
    s: NetworkSolver,
    x: np.ndarray,
    rel: float | np.ndarray = 1e-6,
    floor: float = 1e-9,
) -> np.ndarray:
    """Per-column absolute steps: ``max(rel |x_j|, floor)``, with ``rel`` one
    value or one per column; a pressure column's is then capped at
    ``DROP_FRACTION`` of its smallest pressure gap (see the module docstring)."""
    x = np.asarray(x, dtype=float)
    names = s.unknown_names
    if len(names) != x.size:
        raise ValueError(
            f"{len(names)} unknown names for {x.size} unknowns: build them first (s._build_x0())"
        )
    h = np.maximum(np.abs(x) * np.asarray(rel, dtype=float), floor)

    # Every pressure in the network with its owner (node or element id). A
    # column's gaps skip its own owner: a node's P against its Pt is its
    # dynamic head, which its rows carry smoothly, not a drop.
    p_cols = [j for j, n in enumerate(names) if is_pressure(n)]
    bounds = [(nid, n.Pt) for nid, n in s.network.nodes.items() if isinstance(n, PressureBoundary)]
    owner = np.array([names[j].rsplit(".", 1)[0] for j in p_cols] + [b[0] for b in bounds])
    levels = np.array([x[j] for j in p_cols] + [b[1] for b in bounds], dtype=float)
    for j in p_cols:
        p = abs(x[j])
        gaps = np.abs(levels - x[j])[owner != names[j].rsplit(".", 1)[0]]
        gaps = gaps[gaps > SAME_LEVEL * p]
        if gaps.size:
            h[j] = min(h[j], max(DROP_FRACTION * gaps.min(), P_STEP_FLOOR * p))
    return h


def fd_jacobian(
    s: NetworkSolver,
    x: np.ndarray,
    rel: float | np.ndarray = 1e-6,
    floor: float = 1e-9,
) -> np.ndarray:
    """Finite-difference Jacobian of ``s._residuals`` at ``x``: central, and
    Richardson-extrapolated on the pressure columns."""
    x = np.asarray(x, dtype=float)
    h = fd_steps(s, x, rel, floor)

    def central(j: int, step: float) -> np.ndarray:
        xp, xm = x.copy(), x.copy()
        xp[j] += step
        xm[j] -= step
        return (np.asarray(s._residuals(xp)) - np.asarray(s._residuals(xm))) / (xp[j] - xm[j])

    cols = []
    for j, name in enumerate(s.unknown_names):
        d = central(j, h[j])
        if is_pressure(name):
            d = (4.0 * d - central(j, 2.0 * h[j])) / 3.0
        cols.append(d)
    return np.column_stack(cols)


def scaled_error(J: np.ndarray, fd: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Per row: max_j |J_ij - fd_ij| |x_j| over max_j |fd_ij| |x_j|."""
    J = J.toarray() if hasattr(J, "toarray") else np.asarray(J)
    cols = np.maximum(np.abs(np.asarray(x, dtype=float)), 1e-12)[None, :]
    scale = np.maximum(np.max(np.abs(fd * cols), axis=1), 1e-300)
    return np.max(np.abs((J - fd) * cols), axis=1) / scale
