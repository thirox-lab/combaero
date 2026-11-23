"""hybr on a Jacobian that is singular but for underflow-level entries (#1027).

The hybrj twin of #1020. On stress effusion 105 (an EffusionLiner, 5
segments) the plate rows are the only ones fixing each chamber's pressure
level, and their McGreehan-Schotsch crossflow Cd = exp(-Rv^n) underflows at a
far iterate: the row of chamber ln__b1 carries d/dP = 7.6e-309. In the
Jacobian hybrj requested at its 48th evaluation that entry became R's pivot.
hybrj's dogleg guards only an exactly zero pivot, so the Gauss-Newton step
overflowed, the dogleg step and the trust radius went NaN, and evaluations
49-56 were all NaN until hybr gave up ("not making good progress").

The test runs effusion 103, the other seed #1027 found. A far path like
this is round-off sensitive: by #1198, 105 had left it, and when #1223 moved
the test to 103, 103 without the flush went NaN from its 25th evaluation (8
of 32). At 2aa860fd it no longer does: without the flush, 103's primary hybr
runs 46 evaluations, none NaN (with it, 53), and no effusion seed 0-199 goes
NaN (#1298). So this test does not currently bind to the flush; it guards
the finite-iterate property until a binding case replaces it.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import combaero.network.solver as solver_mod
from combaero.network import NetworkSolver

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
from networks import FAMILIES  # noqa: E402


class _PrimaryHybrDone(BaseException):
    """Ends the solve once the primary hybr attempt returns. A BaseException,
    so no ``except Exception`` in the solver catches it."""


def test_primary_hybr_never_proposes_a_non_finite_iterate(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Effusion 103 cold: the primary hybr attempt sees only finite iterates
    (when #1223 chose the seed, without the trace flush 8 of its 32
    evaluations were NaN; see the module docstring and #1298).

    No wall-clock dependence (#1118): the solver's budgets and stall
    handover read the clock, and a loaded machine ended the primary attempt
    early (with a clock advancing 0.5 s per read it stopped after 12 of its
    53 evaluations), so the precondition tripped under parallel pytest. The
    solver's clock is frozen here, so every budget stays whole and no stall
    is seen: the path of an unloaded machine, ended by iteration limits
    alone. The solve stops once the primary attempt returns; the retries
    after it are not what this pins.

    On main at 2aa860fd the case no longer reaches the underflow pivot: with
    the flush removed it still passes, as does every effusion seed 0-199
    (#1298). It stays as a guard on the finite-iterate property until a
    binding case replaces it."""
    hybr_x: list[np.ndarray] = []
    phase = {"primary": True}
    real_root = solver_mod.root

    def traced_root(fun: Any, x0: np.ndarray, method: str = "hybr", **kw: Any) -> Any:
        if method != "hybr" or not phase["primary"]:
            phase["primary"] = False
            return real_root(fun, x0, method=method, **kw)
        phase["primary"] = False

        def f(x: np.ndarray) -> Any:
            hybr_x.append(np.array(x, copy=True))
            return fun(x)

        # From a finally, so an exception inside hybr also ends the solve
        # (chained) instead of leaving it to run unbudgeted on the frozen clock.
        try:
            real_root(f, x0, method=method, **kw)
        finally:
            raise _PrimaryHybrDone

    monkeypatch.setattr(solver_mod, "root", traced_root)
    t_frozen = solver_mod.time.time()
    monkeypatch.setattr(solver_mod.time, "time", lambda: t_frozen)
    monkeypatch.setattr(solver_mod.time, "perf_counter", lambda: t_frozen)
    net, _ = FAMILIES["effusion"](103)
    solver = NetworkSolver(net)
    solver.homotopy_retry = False
    solver.bc_ramp_retry = False
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(_PrimaryHybrDone):
            solver.solve(timeout=60.0)

    # Enough evaluations to pass the point where the NaN step was taken.
    assert len(hybr_x) > 25, f"hybr ran only {len(hybr_x)} evaluations: case changed"
    n_bad = sum(not np.all(np.isfinite(x)) for x in hybr_x)
    assert n_bad == 0, f"{n_bad} of {len(hybr_x)} hybr iterates non-finite"
