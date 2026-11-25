"""Property tests of the combustion equilibrium kernel (Hypothesis pilot, nightly).

Each recent equilibrium bug was found by hand on a thin slice of composition
space: heptane with a trace of air (#1287), H2 or CO fuel beside a trace
hydrocarbon (#1288), and a lean burn with neither steam nor CO2 (#1293). This
module draws whole mixtures instead: an oxidizer, a fuel at a log-uniform
fraction (or the oxidizer as the trace, at the fuel-rich edge), and
log-uniform traces of the product and fuel species, at a drawn inlet
temperature and pressure. For each it checks what must hold whatever the
mixture:

* the solve converges to a finite temperature, and conserves every element;
* every column of dT/dX that a central difference resolves matches it.

``hypothesis.target`` steers the search toward the largest error. A case the
search finds becomes an issue and, once fixed, a plain regression test beside
the others (``test_combustion_jacobians.py`` holds #1287, #1288 and #1293's).

Measured when this module was written: run unchanged against the commit
before each fix, it fails at #1287's case (heptane with 2e-13 O2), #1288's
(CO with a 1e-6 trace of CH4, not converged) and #1293's (lean CO in dry air
and H2 in O2/N2, Jacobian refused), in every run tried, in under a minute.

Nightly only (``-m nightly``; ``COMBAERO_NIGHTLY=1`` in the scheduled job):
Hypothesis's search budget is spent on the nightly run, not on PR CI.
"""

from __future__ import annotations

import math
import os

import numpy as np
import pytest

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import HealthCheck, given, settings, target  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

import combaero as cb  # noqa: E402

pytestmark = pytest.mark.nightly

NAMES: list[str] = list(cb.species.names)
N = len(NAMES)
IDX = {n: i for i, n in enumerate(NAMES)}

# Atoms per molecule, for the conservation check (test-side data: the names
# are the formulas). AR is argon.
_ATOMS: dict[str, dict[str, int]] = {
    "N2": {"N": 2}, "O2": {"O": 2}, "AR": {"Ar": 1}, "CO2": {"C": 1, "O": 2},
    "H2O": {"H": 2, "O": 1}, "CH4": {"C": 1, "H": 4}, "C2H6": {"C": 2, "H": 6},
    "C3H8": {"C": 3, "H": 8}, "IC4H10": {"C": 4, "H": 10}, "NC5H12": {"C": 5, "H": 12},
    "NC6H14": {"C": 6, "H": 14}, "NC7H16": {"C": 7, "H": 16}, "CO": {"C": 1, "O": 1},
    "H2": {"H": 2}, "NH3": {"N": 1, "H": 3},
}  # fmt: skip
ELEMENTS = ("C", "H", "O", "N", "Ar")
A = np.array([[_ATOMS[n].get(e, 0) for n in NAMES] for e in ELEMENTS], dtype=float)

FUELS = ("CH4", "C2H6", "C3H8", "IC4H10", "NC5H12", "NC6H14", "NC7H16", "CO", "H2")
TRACE_SPECIES = ("H2O", "CO2", "CO", "H2", "CH4", "O2")


def _oxidizer(kind: str) -> np.ndarray:
    X = np.zeros(N)
    if kind == "dry_air":
        X += np.array(cb.species.dry_air())
    elif kind == "humid_air":
        X += np.array(cb.species.humid_air(300.0, 101325.0, 0.6))
    else:  # O2/N2, no argon or CO2
        X[IDX["O2"]], X[IDX["N2"]] = 0.21, 0.79
    return X / X.sum()


@st.composite
def mixtures(
    draw: st.DrawFn, fuel: str, oxidizer: str, edge: str
) -> tuple[float, list[float], float]:
    """(T_in, X, P) for one fuel in one oxidizer: the fuel at a log-uniform
    mole fraction (or, at the fuel-rich edge, the oxidizer at a log-uniform
    fraction down to 1e-12) and up to three log-uniform traces, normalised.

    The fuel and the oxidizer are enumerated by the tests' parametrisation,
    not drawn: Hypothesis's search does not sample uniformly, and a drawn
    categorical axis left whole (fuel, oxidizer) cells unvisited in some runs
    (#1293's lean CO in dry air was found in 3 of 5 runs at 5000 drawn
    examples). Each cell gets its own search over the continuous axes."""
    ox = _oxidizer(oxidizer)
    if edge == "oxidizer":
        x_f = 10.0 ** draw(st.floats(-8.0, math.log10(0.6)))
    else:  # the fuel-rich edge of the simplex: the oxidizer is the trace
        x_f = 1.0 - 10.0 ** draw(st.floats(-12.0, math.log10(0.4)))
    X = (1.0 - x_f) * ox
    X[IDX[fuel]] += x_f
    for name in draw(st.lists(st.sampled_from(TRACE_SPECIES), max_size=3, unique=True)):
        X[IDX[name]] += 10.0 ** draw(st.floats(-20.0, -3.0))
    X = X / X.sum()
    T = draw(st.floats(250.0, 1500.0))
    P = 10.0 ** draw(st.floats(4.0, math.log10(5.0e6)))
    return T, [float(v) for v in X], P


OXIDIZERS = ("dry_air", "humid_air", "o2n2")

# The fuel-rich edge (the oxidizer a trace, down to 1e-12) holds open
# failures: a heavy alkane with trace air and a CO trace does not converge at
# 250 K (#1312), and dT/dX columns there miss their resolved difference. Its
# cells keep searching and report, but do not fail the nightly run until
# #1312 is fixed; then drop the mark so the edge gates too.
_RICH_EDGE_OPEN = pytest.mark.xfail(
    reason="fuel-rich edge: open equilibrium failures (#1312)", strict=False
)
CELLS = pytest.mark.parametrize(
    ("fuel", "oxidizer", "edge"),
    [(f, o, "oxidizer") for f in FUELS for o in OXIDIZERS]
    + [pytest.param(f, o, "fuel", marks=_RICH_EDGE_OPEN) for f in FUELS for o in OXIDIZERS],
)


def _atoms(X: list[float] | np.ndarray) -> np.ndarray:
    return A @ np.asarray(X, dtype=float)


SETTINGS = settings(
    # Per (fuel, oxidizer, edge) cell; the nightly job raises it.
    max_examples=int(os.environ.get("COMBAERO_HYPOTHESIS_EXAMPLES", "200")),
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)


@CELLS
@SETTINGS
@given(data=st.data())
def test_converges_and_conserves_elements(
    fuel: str, oxidizer: str, edge: str, data: st.DataObject
) -> None:
    _conserves(data.draw(mixtures(fuel, oxidizer, edge), label="T_in, X, P"))


def _conserves(case: tuple[float, list[float], float]) -> None:
    T, X, P = case
    r = cb._core.combustion_equilibrium(T, X, P)
    assert r.converged, case
    assert math.isfinite(r.state.T) and r.state.T > 0.0, (case, r.state.T)
    a_in, a_out = _atoms(X), _atoms(r.state.X)
    # Per mole of mixture: total atoms change with the mole count, so compare
    # element ratios to the total, which the reaction must keep.
    frac_in, frac_out = a_in / a_in.sum(), a_out / a_out.sum()
    err = float(np.max(np.abs(frac_out - frac_in)))
    target(err, label="element fraction error")
    assert err <= 1e-10, (case, err)


@CELLS
@SETTINGS
@given(data=st.data())
def test_dT_dX_matches_resolved_central_differences(
    fuel: str, oxidizer: str, edge: str, data: st.DataObject
) -> None:
    _dT_dX_matches(data.draw(mixtures(fuel, oxidizer, edge), label="T_in, X, P"))


def _central(T: float, X: list[float], P: float, k: int, h: float) -> float:
    Xp, Xm = list(X), list(X)
    Xp[k] += h
    Xm[k] -= h
    return (
        cb._core.combustion_equilibrium(T, Xp, P).state.T
        - cb._core.combustion_equilibrium(T, Xm, P).state.T
    ) / (2.0 * h)


def _resolved_fd(T: float, X: list[float], P: float, k: int, scale: float) -> float | None:
    """dT/dX_k from the most converged pair of a Richardson ladder, or None.

    Richardson on central differences at h and h/2 is O(h^4), but the
    curvature it must resolve varies by orders of magnitude: a plain central
    difference at 1e-4 X_k is truncation-limited at about 1e-4 of the row on
    rich heptane, and near stoichiometric methane (where dT/dX turns at the
    temperature peak) even Richardson at 1e-4 X_k is off by 3e-5 and at
    3e-4 X_k by 4e-3, while 3e-5 X_k agrees to 2e-7. So three steps are tried
    and the two adjacent estimates that agree best are used; a column whose
    best pair still differs by more than a tenth of the tolerance is not
    resolved and is skipped."""
    est = []
    for rel in (1e-4, 1e-5, 1e-6):
        h = rel * X[k]
        est.append((4.0 * _central(T, X, P, k, 0.5 * h) - _central(T, X, P, k, h)) / 3.0)
    spread, best = min((abs(a - b), b) for a, b in zip(est[:-1], est[1:], strict=True))
    return best if spread <= 0.1 * TOL * scale else None


TOL = 1e-4  # dT/dX error over the row's largest entry


def _dT_dX_matches(case: tuple[float, list[float], float]) -> None:
    T, X, P = case
    jac = cb._core.combustion_equilibrium_and_jacobian(T, X, P)
    scale = max(abs(v) for v in jac.dT_dX)
    assert math.isfinite(scale), (case, list(jac.dT_dX))
    if scale == 0.0:
        return  # T does not depend on X here (nothing burns): nothing to compare
    worst = 0.0
    for k in range(N):
        if X[k] < 1e-6:
            continue  # only columns a central difference resolves
        fd = _resolved_fd(T, X, P, k, scale)
        if fd is None:
            continue  # no step resolves it: not a test of J
        worst = max(worst, abs(jac.dT_dX[k] - fd) / scale)
    target(worst, label="dT/dX error over row scale")
    assert worst <= TOL, (case, worst)
