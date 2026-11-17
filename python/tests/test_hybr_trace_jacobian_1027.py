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
this is round-off sensitive: by #1198, 105 had left it (with the trace flush
removed, hybr proposed no NaN iterate there, before and after #1198), and
#1198's 1e-16 Jacobian change cut its run below the old 48-evaluation
precondition. On 103, without the flush, hybr's 25th evaluation is the first
NaN, 8 of 32 in all, before and after #1198; with it, 0 of 53. Re-check the
seed (flush removed, it must still go NaN) whenever this test's precondition
trips.
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


def test_primary_hybr_never_proposes_a_non_finite_iterate(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Effusion 103 cold, minus the retries: the primary hybr attempt sees
    only finite iterates (without the trace flush: 8 of its 32 evaluations
    NaN, the first at evaluation 25)."""
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

        return real_root(f, x0, method=method, **kw)

    monkeypatch.setattr(solver_mod, "root", traced_root)
    net, _ = FAMILIES["effusion"](103)
    solver = NetworkSolver(net)
    solver.homotopy_retry = False
    solver.bc_ramp_retry = False
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = solver.solve(timeout=60.0)

    # Enough evaluations to pass the point where the NaN step was taken.
    assert len(hybr_x) > 25, f"hybr ran only {len(hybr_x)} evaluations: case changed"
    n_bad = sum(not np.all(np.isfinite(x)) for x in hybr_x)
    assert n_bad == 0, f"{n_bad} of {len(hybr_x)} hybr iterates non-finite"
    assert res["__converged__"], res.get("__message__")
