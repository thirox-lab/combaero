"""One finite-difference Jacobian for the network tests (#1107).

A central difference has a truncation error of about ``(h / s)^2``, where
``s`` is the pressure scale the residual curves on, not the pressure itself.
On ``m ~ sqrt(dp)`` that scale is the DROP: a step of ``1e-6 |P|`` is 3.1 Pa
at 31 bar, and on a 42.7 Pa duct it read as a 6e-4 "Jacobian error" that was
the step (#1107). The same rule as the sweep warm start: a step must be small
against dp, not against P (#1054).

Pressure columns (``fd_jacobian``, the default)
-----------------------------------------------
Each pressure column is differenced on a LADDER of steps halving from twice
the caller's ``max(rel |x_j|, floor)``, every level Richardson-extrapolated
(``R(h) = (4 D(h) - D(2h)) / 3``, which cancels the h^2 term). Each ROW then
takes its own level, between two bounds:

* above, for a duct's row, its drive: every call the row's element makes to
  the Fanno kernels is recorded with its drive (``Pt0 - P_target``, signed
  by the duct's direction), and a level is allowed only while each drive
  stays on its own side of the 1 Pa flow-drive floor
  (``FANNO_FLOW_DRIVE_FLOOR``; the flux is C1 there, so G'' jumps), its
  stencil moves the drive by at most ``FLOOR_MARGIN`` of its distance to the
  floor and by at most 1 / ``DRIVE_DIVISOR`` of the drive itself (of the
  floor, inside it). Above the floor, the row's curvature scale is its drive,
  so it takes the largest allowed level (#1157: s2x's drive 1.06 Pa, the old
  gap-capped step straddled the joint and read 3.8e-3). Below the floor
  (#1216), a creeping duct's real laminar-to-turbulent curvature scale can be
  orders of magnitude below 1 Pa (e.g. ~1e-3 Pa), so sub-floor duct rows use
  adaptive ``_pick`` like non-duct rows, may extend down to
  ``P_SUBFLOOR_STEP_FLOOR |P|`` (``5e-10 |P|``), when the 4th-order
  Richardson spread remains above ``UNRESOLVED`` may take the 6th-order
  Neville-Richardson column ``R2`` on the same ``D`` ladder (zero extra
  residual evaluations) only if its 3-window spread is strictly tighter, and
  if still above ``UNRESOLVED`` at the bottom of the 2x halving ladder
  inserts geometric midpoints (``sqrt(2)`` step ratio) around the unresolved
  step window so a 3-window fits inside the narrow sweet spot between a
  ~1e-3 Pa zero-crossing / transition knee and IEEE-754 pressure round-off;
* below, round-off: the ladder never steps below ``P_STEP_FLOOR |P|`` (or
  ``P_SUBFLOOR_STEP_FLOOR |P|`` for a sub-floor duct row).

Any other row takes the tightest window of three consecutive levels
(``_pick``): its spread estimates the error (truncation shows as a trend,
falling ~16x per level, round-off as scatter), the largest-step window within
``PREFER_LARGER`` of the tightest is taken, and its median is the value. The
ladder starts with ``FIRST_LEVELS`` levels and extends down, two at a time,
while some row's tightest window is its smallest and above ``ACCURACY`` of its
scale, or its estimate is above ``UNRESOLVED`` (the levels so far may all lie
outside the asymptotic range). So the step follows the row's own curvature (a
step large against the row's drop moves R by ~(h / dp)^4) and its own
round-off (a step pinned near the floor reads platform-dependent noise:
#1187's c2.P column read 7.9e-8 on macOS and 6e-7 on CI x86_64 at 1e-8 |P|),
instead of the gaps between pressure unknowns standing in for both. The
choice reads only F, never the analytic Jacobian, so it cannot hide a wrong
entry: a wrong partial misses at every level.

A row with no allowed level (a drive within a few round-off steps of the
floor) falls back to a one-sided Richardson difference, ``2 D+(h) - D+(2h)``,
on the side the drive moves AWAY from the floor, on the same ladder.

Each entry's base step and what bound it are kept on the result
(``FDJacobian.why(i, j)``; ``worst`` puts it in a failure message): "drive
cap of <element>, drive <d> Pa" (a duct's row at its largest allowed level;
"top level, under the drive cap" where the cap cut nothing), "fit" (the
window the spread chose), "roundoff" (the ladder reached ``P_STEP_FLOOR |P|``
with its tightest window the smallest), "one-sided <element>", "plain" for a
non-pressure column.

Kernel warm starts (``_last_G``, ``_last_M_exit``) are cleared before every
evaluation: with them, one x read twice differs by ~1e-11, the noise in a
structurally zero column (#1157), and a warm-started flow column needed a
1e-5 step. F is then a function of x up to the state other fixed points
keep: a ThermalWall with k(T) layers starts its conductivity iteration from
the previous call's (``components.py``, converged to ~1e-12 relative).

Cost: a pressure column takes 8 residual evaluations, 4 more per extension
(at most 18 from the default step), against 4 for the old rule.

Other columns, and the old rule
-------------------------------
Non-pressure columns keep a plain central difference at ``max(rel |x_j|,
floor)``; callers may pass ``rel`` per column.

``fd_jacobian(..., adaptive=False)`` is the previous rule (``fd_steps``): one
step per pressure column, capped at ``DROP_FRACTION`` of the smallest gap
between its value and any other pressure in the network, never below
``P_STEP_FLOOR |P|``. That cap pinned columns at the round-off floor wherever
two levels sat close (#1187: a 0.04 Pa gap; #1177: 0.14 Pa, a 1.4e-3 Pa step).

Compare in scaled variables, ``J_ij |x_j|`` (RESIDUAL-ROW-SCALING, see
``NetworkSolver._build_residual_scales``): ``scaled_error``, or ``worst`` for
the worst row with the step that produced it.
"""

from __future__ import annotations

import sys
from collections.abc import Callable, Iterator
from contextlib import contextmanager

import numpy as np

import combaero as cb
from combaero.network import NetworkSolver, PressureBoundary

#: A pressure step as a fraction of the smallest pressure gap (the old rule,
#: ``adaptive=False``). With the extrapolation the error is ~ (h / dp)^4.
DROP_FRACTION = 1e-2
#: Gaps below this fraction of |P| count as the same level (lossless links).
SAME_LEVEL = 1e-9
#: No pressure step below this fraction of |P|: a film-cooling row capped at
#: 1.5e-10 |P| read 3.5e-6 of round-off, against 2e-7 at 1.5e-9 |P|. The
#: ladder's bottom level.
P_STEP_FLOOR = 1e-8
#: Sub-floor ducts (#1216) solve the low-Mach balance directly in drive (no
#: Pt0 - P cancellation inside the kernel) and can sit at drives ~1e-3 Pa;
#: allowing their ladder to extend down to 1e-9 |P| keeps them on central
#: Richardson levels instead of falling back to one-sided O(h^2) steps.
P_SUBFLOOR_STEP_FLOOR = 5e-10
#: The Fanno kernels' joints [Pa]: the flow-drive floor either way, where
#: the flux is C1 (``fanno_mach.h``). Zero drive is not one: the flux is an
#: odd cubic through it, and a duct whose ends hold the same gas reverses
#: smoothly (at #1157's reversed root x1 sits at -0.002 Pa; differenced
#: across zero it meets the 1e-7 check, one-sided away from it it read
#: 3.4e-6).
DRIVE_JOINTS = (cb.FANNO_FLOW_DRIVE_FLOOR, -cb.FANNO_FLOW_DRIVE_FLOOR)
#: A level's stencil (+-2h) may move a drive by at most this fraction of its
#: distance to the nearest joint, and by at most 1 / DRIVE_DIVISOR of the
#: drive (of the floor below it): h at most 1/5 of the distance and 1/50 of
#: the drive, the cap #1157 measured by hand (the flux ~ sqrt(drive), whose
#: Richardson error is ~0.2 (h / drive)^4).
FLOOR_MARGIN = 0.4
DRIVE_DIVISOR = 25.0
#: The ladder's first levels (D(h) values, halving from twice the caller's
#: step): three extrapolations, the caller's own (base rel |x|) the first.
FIRST_LEVELS = 4
#: A row extends its column's ladder down while its tightest window is the
#: smallest one and its estimate is above this fraction of its scale
#: (scaled, as ``scaled_error``): 1e-3 of the tightest tolerance the tests
#: use (1e-7).
ACCURACY = 1e-10
#: ... and whatever its tightest window, while its estimate is above this
#: (scaled): the levels so far may all lie outside the asymptotic range
#: (#695's ln__s3, ``_pick``).
UNRESOLVED = 1e-7

#: The window taken is the largest-step one whose spread is within this
#: factor of the smallest (``_pick``): round-off falls with a larger step,
#: and a window of round-off can agree by luck at a small one.
PREFER_LARGER = 2.0

#: The kernels the network's ducts call, and the position of P_target in
#: their arguments (``Pt0`` is the first).
_KERNELS = {"fanno_channel_flow": 3, "fanno_rayleigh_channel_flow": 4}
#: The kernels' warm starts (``solver._KERNEL_WARM_STARTS``).
_WARM_STARTS = ("_last_G", "_last_M_exit")


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
    """The old rule's per-column steps (``adaptive=False``): ``max(rel |x_j|,
    floor)``, with ``rel`` one value or one per column; a pressure column's is
    then capped at ``DROP_FRACTION`` of its smallest pressure gap."""
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


class FDJacobian(np.ndarray):
    """The FD Jacobian, with each entry's step (``steps``), what bound it
    (``bounds``) and its estimated error (``est``, nan where none), scaled as
    ``scaled_error``: ``|J_ij| |x_j|`` units over the row's largest
    ``|fd_ij x_j|``."""

    steps: np.ndarray
    bounds: np.ndarray
    est: np.ndarray
    names: list[str]

    def __array_finalize__(self, obj: object) -> None:
        for k in ("steps", "bounds", "est", "names"):
            setattr(self, k, getattr(obj, k, None))

    def why(self, i: int, j: int) -> str:
        if self.steps is None or self.steps.shape != self.shape:
            return f"column {j}"
        e = self.est[i, j]
        est = "" if not np.isfinite(e) else f", est. scaled error {e:.1e}"
        return f"{self.names[j]} step {self.steps[i, j]:.3g} ({self.bounds[i, j]}{est})"


def _wrap(name: str, fn: Callable, record: list) -> Callable:
    pos = _KERNELS[name]

    def recorder(*args, **kwargs):  # noqa: ANN002, ANN003, ANN202
        Pt0 = kwargs.get("Pt0", args[0] if args else np.nan)
        P_t = kwargs.get("P_target", args[pos] if len(args) > pos else np.nan)
        frame = sys._getframe(1)
        el, fwd = None, True
        for _ in range(4):
            if frame is None:
                break
            loc = frame.f_locals
            if el is None and hasattr(loc.get("self"), "id"):
                el = loc["self"].id
                fwd = bool(loc.get("forward", True))
            frame = frame.f_back
        drive = float(Pt0) - float(P_t)
        record.append((el, drive if fwd else -drive))
        return fn(*args, **kwargs)

    return recorder


@contextmanager
def _recording(record: list) -> Iterator[None]:
    """Record (element id, signed drive) for every Fanno kernel call."""
    saved = {k: getattr(cb, k) for k in _KERNELS}
    try:
        for k, fn in saved.items():
            setattr(cb, k, _wrap(k, fn, record))
        yield
    finally:
        for k, fn in saved.items():
            setattr(cb, k, fn)


def _evaluate(s: NetworkSolver, x: np.ndarray) -> tuple[np.ndarray, dict[str, list[float]]]:
    """F at x with the kernels started cold, and each element's drives in
    call order."""
    for el in getattr(getattr(s, "network", None), "elements", {}).values():
        for key in _WARM_STARTS:
            if key in vars(el):
                delattr(el, key)
    record: list = []
    with _recording(record):
        F = np.asarray(s._residuals(x), dtype=float)
    drives: dict[str, list[float]] = {}
    for el, d in record:
        drives.setdefault(el, []).append(d)
    return F, drives


def _margin(d: float) -> float:
    return min(abs(d - jt) for jt in DRIVE_JOINTS)


def _side(d: float) -> tuple[bool, ...]:
    return tuple(d > jt for jt in DRIVE_JOINTS)


def _safe(d0: list[float], dp: list[float] | None, dm: list[float] | None) -> bool:
    """Whether both sides of a stencil keep every drive on its side of the
    joints, within FLOOR_MARGIN of its distance to them and within
    1 / DRIVE_DIVISOR of its own curvature scale (the drive, or the floor
    inside it, where ``_pick`` resolves sub-floor curvature adaptively)."""
    for d in (dp, dm):
        if d is None or len(d) != len(d0):
            return False
        for a, b in zip(d0, d, strict=True):
            if _side(a) != _side(b) or abs(b - a) > FLOOR_MARGIN * _margin(a):
                return False
            if abs(b - a) > max(abs(a), cb.FANNO_FLOW_DRIVE_FLOOR) / DRIVE_DIVISOR:
                return False
    return True


def _row_elements(s: NetworkSolver, x: np.ndarray, n_rows: int) -> list[str | None]:
    """The element that owns each residual row (None: a node's row)."""
    try:
        kinds = s._row_kinds(x)
    except AttributeError:
        return [None] * n_rows
    # Owners read "<Type> '<id>'": the id runs from the first quote to the
    # last, so an id holding a quote parses whole.
    out = [
        None if o.startswith("node") or o.count("'") < 2 else o[o.index("'") + 1 : o.rindex("'")]
        for o, _ in kinds
    ]
    return out if len(out) == n_rows else [None] * n_rows


def _pick(R: np.ndarray, h: np.ndarray, allowed: np.ndarray) -> tuple[float, int, float, bool]:
    """One row's extrapolation from its allowed levels: (value, index of the
    level whose base step it reports, estimate, whether smaller steps might
    do better). ``R[k]`` is the extrapolation at base step ``h[k]``
    (decreasing).

    The estimate of a window of three consecutive allowed levels is its
    spread (max - min): truncation shows there as a trend (it falls ~16x per
    level), round-off as scatter, and a window agrees by luck only if two of
    its pairs do. The window taken is the largest-step one whose spread is
    within PREFER_LARGER of the smallest, its value the median of its three.
    A pair alone misled: an entry that is round-off at every step (#1157's
    s4x row in s5.P, an off-diagonal J of 3e-5) had two neighbouring levels
    agree by luck at 3e-3 Pa, scaled 1.1e-7, against 2e-9 at 0.05 Pa; two
    levels outside the asymptotic range agree too (#695's liner segment
    ln__s3: 0.47 and 0.235 Pa, 8e-2 off J, which the FD meets below 0.015
    Pa). With two allowed levels the smaller step's value is taken (estimate
    their difference), with one that level's (nan)."""
    idx = np.flatnonzero(allowed)
    if idx.size == 1:
        return float(R[idx[0]]), int(idx[0]), np.nan, True
    if idx.size == 2:
        return float(R[idx[1]]), int(idx[1]), float(abs(R[idx[1]] - R[idx[0]])), True
    win = np.array([R[idx[n : n + 3]] for n in range(idx.size - 2)])
    spread = win.max(axis=1) - win.min(axis=1)
    best = int(np.argmin(spread))
    n = int(np.flatnonzero(spread <= PREFER_LARGER * spread[best])[0])
    # Smaller steps may help while the tightest window is the smallest one.
    more = best == spread.size - 1
    return float(np.median(win[n])), int(idx[n + 1]), float(spread[n]), more


class _Column:
    """One pressure column's ladder: per level its step, central D and the
    drives at both sides of the stencil."""

    def __init__(self, s: NetworkSolver, x: np.ndarray, j: int, top: float) -> None:
        self.s, self.x, self.j = s, x, j
        self.lo = P_STEP_FLOOR * abs(x[j])
        self.h: list[float] = []
        self.D: list[np.ndarray] = []
        self.dp: list[dict] = []
        self.dm: list[dict] = []
        self._add(top)
        while len(self.h) < FIRST_LEVELS:
            self._add(self.h[-1] / 2.0)

    def _add(self, h: float) -> None:
        xp, xm = self.x.copy(), self.x.copy()
        xp[self.j] += h
        xm[self.j] -= h
        (Fp, dp), (Fm, dm) = _evaluate(self.s, xp), _evaluate(self.s, xm)
        self.h.append(h)
        self.D.append((Fp - Fm) / (xp[self.j] - xm[self.j]))
        self.dp.append(dp)
        self.dm.append(dm)

    def can_extend(self, floor_frac: float = P_STEP_FLOOR) -> bool:
        lo = floor_frac * abs(self.x[self.j])
        return self.h[-1] > lo * (1.0 + 1e-12)

    def extend(self, floor_frac: float = P_STEP_FLOOR) -> None:
        """Two more levels, halving, the last clipped at the round-off floor."""
        lo = floor_frac * abs(self.x[self.j])
        for _ in range(2):
            if self.can_extend(floor_frac):
                self._add(max(self.h[-1] / 2.0, lo))

    def refine_around(self, target_steps: list[float]) -> bool:
        """Insert geometric midpoints (sqrt(2) step ratio) for intervals
        bracketing ``target_steps``, keeping levels sorted by decreasing h.
        Returns True if any new level was added."""
        mids: set[float] = set()
        h_arr = np.array(self.h)
        for h_t in target_steps:
            k_near = int(np.argmin(np.abs(h_arr - h_t)))
            for k in range(max(0, k_near - 2), min(len(self.h) - 1, k_near + 2)):
                if self.h[k] / self.h[k + 1] > 1.6:
                    mids.add(float(np.sqrt(self.h[k] * self.h[k + 1])))
        if not mids:
            return False
        for hm in sorted(mids, reverse=True):
            self._add(hm)
        order = sorted(range(len(self.h)), key=lambda idx: self.h[idx], reverse=True)
        self.h = [self.h[idx] for idx in order]
        self.D = [self.D[idx] for idx in order]
        self.dp = [self.dp[idx] for idx in order]
        self.dm = [self.dm[idx] for idx in order]
        return True

    def select(
        self,
        d0: dict[str, list[float]],
        row_el: list[str | None],
        scale: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], np.ndarray]:
        """Per row: value, base step, estimate, bound, and whether the row
        wants more levels (estimate above ACCURACY of its scale, or no allowed
        level)."""
        h = np.array(self.h)
        D = np.array(self.D)
        r2 = (h[:-1] / h[1:])[:, None] ** 2
        R = (r2 * D[1:] - D[:-1]) / (r2 - 1.0)
        # Second Neville-Richardson column (O(h^6) truncation, zero extra F
        # evaluations). Between R[k] (error ~ h[k]^2 h[k+1]^2) and R[k+1]
        # (error ~ h[k+1]^2 h[k+2]^2) the exact cancellation ratio is
        # (h[k] / h[k+2])^2 (= 16 on halving steps, exact even when the last
        # step is clipped at the floor).
        r4 = (h[:-2] / h[2:])[:, None] ** 2
        R2 = (r4 * R[1:] - R[:-1]) / (r4 - 1.0)
        safe = {
            el: np.array([_safe(dd, self.dp[k].get(el), self.dm[k].get(el)) for k in range(h.size)])
            for el, dd in d0.items()
        }
        m = D.shape[1]
        val, step, est = np.zeros(m), np.zeros(m), np.full(m, np.nan)
        why: list[str] = [""] * m
        more = np.zeros(m, dtype=bool)
        xj = abs(self.x[self.j])
        for i in range(m):
            el = row_el[i]
            ok = np.ones(h.size - 1, dtype=bool)
            if el in safe:
                ok = safe[el][:-1] & safe[el][1:]
            if not ok.any():
                more[i] = True
                why[i] = "none"
                continue
            allowed = np.flatnonzero(ok)
            # Above the drive floor a duct's curvature scale is its drive, which
            # the cap reads, so it takes the largest allowed step (round-off
            # falls with it). Below the floor (#1216) the low-Mach inertial
            # crossover can be orders of magnitude smaller than 1 Pa, so let
            # _pick resolve the row's curvature adaptively.
            sub_floor = el in safe and any(abs(d) < cb.FANNO_FLOW_DRIVE_FLOOR for d in d0[el])
            if el in safe and not sub_floor:
                k = int(allowed[0])
                e = float(abs(R[k, i] - R[allowed[1], i])) if allowed.size > 1 else np.nan
                v, smaller = float(R[k, i]), False
            else:
                v, k, e, smaller = _pick(R[:, i], h[1:], ok)
                # When a sub-floor duct's O(h^4) truncation meets IEEE-754
                # pressure round-off above UNRESOLVED (curving on a ~1e-3 Pa
                # scale against 1e5 Pa pressures), check the O(h^6) column R2
                # on the same D ladder and adopt it only if its 3-window spread
                # is strictly tighter (never peeking at J).
                ok2 = ok[:-1] & ok[1:]
                rel_e0 = e * xj / scale[i] if scale[i] > 0.0 else 0.0
                if sub_floor and np.count_nonzero(ok2) >= 3 and not (rel_e0 <= UNRESOLVED):
                    v2, k2, e2, smaller2 = _pick(R2[:, i], h[2:], ok2)
                    if np.isfinite(e2) and (not np.isfinite(e) or e2 < e):
                        v, k, e, smaller = v2, k2 + 1, e2, smaller2
            val[i], step[i], est[i] = v, h[k + 1], e
            if el in safe and not sub_floor:
                cut = "drive cap" if not ok.all() else "top level, under the drive cap"
                reason = f"{cut} of {el}, drive {min(d0[el], key=_margin):.4g} Pa"
            else:
                reason = "roundoff" if smaller and not self.can_extend() else "fit"
            why[i] = reason
            # Smaller steps may help while the tightest window is the
            # smallest one and above ACCURACY, and whenever the estimate is
            # above UNRESOLVED: then the asymptotic range may lie below.
            rel_e = e * xj / scale[i] if scale[i] > 0.0 else 0.0
            more[i] = (el not in safe or sub_floor) and (
                not (rel_e <= ACCURACY) and (smaller or not (rel_e <= UNRESOLVED))
            )
        return val, step, est, why, more

    def one_sided(
        self, F0: np.ndarray, d0: list[float], el: str
    ) -> tuple[np.ndarray, np.ndarray, float]:
        """Forward or backward Richardson ``(r D+(h) - D+(r h)) / (r - 1)`` on
        the side that moves the element's nearest-joint drive away from its
        joint, on the ladder's steps: (R per level, base steps, the side)."""
        q = min(range(len(d0)), key=lambda n: _margin(d0[n]))
        a = d0[q]
        joint = min(DRIVE_JOINTS, key=lambda jt: abs(a - jt))
        probe = self.dp[-1].get(el, [])
        slope = probe[q] - a if len(probe) == len(d0) else 0.0
        away = 1.0 if (a - joint) * slope >= 0.0 else -1.0
        h = np.array(self.h)
        Dp = []
        for hk in h:
            xs = self.x.copy()
            xs[self.j] += away * hk
            Dp.append((_evaluate(self.s, xs)[0] - F0) / (xs[self.j] - self.x[self.j]))
        Dp_a = np.array(Dp)
        r = (h[:-1] / h[1:])[:, None]
        return (r * Dp_a[1:] - Dp_a[:-1]) / (r - 1.0), h[1:], away


def _old_rule(s: NetworkSolver, x: np.ndarray, rel: float | np.ndarray, floor: float) -> np.ndarray:
    """The previous fd_jacobian: one step per column (``fd_steps``), central,
    Richardson on the pressure columns, warm starts as the solver left them."""
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


def fd_jacobian(
    s: NetworkSolver,
    x: np.ndarray,
    rel: float | np.ndarray = 1e-6,
    floor: float = 1e-9,
    *,
    adaptive: bool = True,
) -> FDJacobian:
    """Finite-difference Jacobian of ``s._residuals`` at ``x``: central, and
    on the pressure columns a floor-aware Richardson ladder per row (module
    docstring); ``adaptive=False`` is the old one-step rule (``fd_steps``)."""
    x = np.asarray(x, dtype=float)
    names = list(s.unknown_names)
    n = x.size
    if len(names) != n:
        raise ValueError(
            f"{len(names)} unknown names for {n} unknowns: build them first (s._build_x0())"
        )
    if not adaptive:
        cols = _old_rule(s, x, rel, floor)
        out = cols.view(FDJacobian)
        out.steps = np.broadcast_to(fd_steps(s, x, rel, floor), cols.shape).copy()
        out.bounds = np.full(cols.shape, "old rule", dtype=object)
        out.est, out.names = np.full(cols.shape, np.nan), names
        return out

    rel_a = np.broadcast_to(np.asarray(rel, dtype=float), (n,))
    F0, d0 = _evaluate(s, x)
    m = F0.size
    cols = np.zeros((m, n))
    steps = np.zeros((m, n))
    bounds = np.empty((m, n), dtype=object)
    est = np.full((m, n), np.nan)
    row_el = _row_elements(s, x, m)

    ladders: dict[int, _Column] = {}
    for j, name in enumerate(names):
        h = max(abs(x[j]) * rel_a[j], floor)
        if not is_pressure(name):
            xp, xm = x.copy(), x.copy()
            xp[j] += h
            xm[j] -= h
            cols[:, j] = (_evaluate(s, xp)[0] - _evaluate(s, xm)[0]) / (xp[j] - xm[j])
            steps[:, j], bounds[:, j] = h, "plain"
            continue
        # The top level is twice the caller's step: the first extrapolation
        # is then the caller's own (base rel |x|, stencil 2 rel |x|).
        ladders[j] = _Column(s, x, j, 2.0 * h)
        cols[:, j] = ladders[j].select(d0, row_el, np.full(m, np.inf))[0]

    # Each row's scale (RESIDUAL-ROW-SCALING): its largest |fd_ij x_j|, from
    # the first levels. A row whose estimate is above ACCURACY of it extends
    # its column's ladder down, to the round-off floor at most.
    scale = np.max(np.abs(cols) * np.maximum(np.abs(x), 1e-12)[None, :], axis=1)
    for j, col in ladders.items():
        while True:
            val, step, e, why, more = col.select(d0, row_el, scale)
            if not more.any():
                break
            need_subfloor = any(
                more[i]
                and row_el[i] in d0
                and any(abs(d) < cb.FANNO_FLOW_DRIVE_FLOOR for d in d0[row_el[i]])
                for i in range(m)
            )
            floor_frac = P_SUBFLOOR_STEP_FLOOR if need_subfloor else P_STEP_FLOOR
            if not col.can_extend(floor_frac):
                # Sub-floor duct rows sitting near ~1e-3 Pa can have a narrow
                # sweet spot between a zero-drive curvature knee and IEEE-754
                # pressure round-off where a 2x halving 3-window spans 8x-16x.
                # Refine unresolved sub-floor intervals once with sqrt(2)
                # geometric midpoints (reading only F's own spread, never J).
                xj_abs = abs(x[j])
                unresolved_steps = [
                    float(step[i])
                    for i in range(m)
                    if more[i]
                    and row_el[i] in d0
                    and any(abs(d) < cb.FANNO_FLOW_DRIVE_FLOOR for d in d0[row_el[i]])
                    and scale[i] > 0.0
                    and not (e[i] * xj_abs / scale[i] <= UNRESOLVED)
                    and step[i] > 0.0
                ]
                if unresolved_steps and col.refine_around(unresolved_steps):
                    val, step, e, why, more = col.select(d0, row_el, scale)
                break
            col.extend(floor_frac)
        cols[:, j], steps[:, j], est[:, j] = val, step, e
        bounds[:, j] = why
        for i in np.flatnonzero(np.array(why) == "none"):
            el = row_el[i]
            Rs, hs, away = col.one_sided(F0, d0[el], el)
            v, k, ek, _ = _pick(Rs[:, i], hs, np.ones(hs.size, dtype=bool))
            cols[i, j], steps[i, j], est[i, j] = v, away * hs[k], ek
            bounds[i, j] = f"one-sided {el}, drive {min(d0[el], key=_margin):.4g} Pa"

    # The estimates in scaled units, as ``scaled_error`` reports.
    xs = np.maximum(np.abs(x), 1e-12)[None, :]
    row = np.maximum(np.max(np.abs(cols) * xs, axis=1), 1e-300)[:, None]
    out = cols.view(FDJacobian)
    out.steps, out.bounds, out.est, out.names = steps, bounds, est * xs / row, names
    return out


def scaled_error(J: np.ndarray, fd: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Per row: max_j |J_ij - fd_ij| |x_j| over max_j |fd_ij| |x_j|."""
    J = J.toarray() if hasattr(J, "toarray") else np.asarray(J)
    fd_a = np.asarray(fd)
    cols = np.maximum(np.abs(np.asarray(x, dtype=float)), 1e-12)[None, :]
    scale = np.maximum(np.max(np.abs(fd_a * cols), axis=1), 1e-300)
    return np.max(np.abs((J - fd_a) * cols), axis=1) / scale


def worst(
    J: np.ndarray, fd: np.ndarray, x: np.ndarray, rows: list[str] | None = None
) -> tuple[float, str]:
    """The worst row's scaled error and a message naming the row, the column
    that sets it and (from an ``FDJacobian``) that column's step and bound."""
    J = J.toarray() if hasattr(J, "toarray") else np.asarray(J)
    fd_a = np.asarray(fd)
    e = scaled_error(J, fd_a, x)
    i = int(np.argmax(e))
    cols = np.maximum(np.abs(np.asarray(x, dtype=float)), 1e-12)
    j = int(np.argmax(np.abs((J[i] - fd_a[i]) * cols)))
    row = rows[i] if rows is not None else f"row {i}"
    why = fd.why(i, j) if isinstance(fd, FDJacobian) else f"column {j}"
    return float(e[i]), f"{row}: scaled {e[i]:.2e} at {why}"
