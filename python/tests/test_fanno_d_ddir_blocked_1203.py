"""dG/dv along a composition direction that no stencil can step (#1203).

`d_ddir` (src/fanno_mach.cpp) differences a duct's flux along a mass-fraction
direction v. Where v has opposite-sign components on two species that are
absent or below h|v_k|, neither side of the step keeps every fraction >= 0,
and it returned 0: the duct's J entry was 0 where the true dG/dv is 58-66.
dG/dv is linear in v, so the direction is now split into its blocking species
(each a unit vector, feasible on its positive side) and the rest, and the
parts are summed. The references: that linearity (sum_k v_k dG/dY_k from
per-species unit directions -- for an all-blocked direction this is the same
arithmetic, so it only shows the entry is no longer 0), and, independent of
d_ddir, the limit of plain central differences along v from feasible traces.
"""

import subprocess
import sys

import numpy as np
import pytest

import combaero as cb

_N = [cb.species_name(k) for k in range(cb.num_species())]
_I = {n: k for k, n in enumerate(_N)}


def _gas(trace: dict[str, float]) -> list[float]:
    """Dry air with its CO2 moved into N2 (as in #1203), plus mass-fraction traces."""
    Y = np.zeros(len(_N))
    for name, y in {"N2": 0.7557, "O2": 0.2314, "AR": 0.0129, **trace}.items():
        Y[_I[name]] = y
    Y /= Y.sum()
    return list(cb.mass_to_mole(list(Y)))


def _vec(d: dict[str, float]) -> list[float]:
    v = [0.0] * len(_N)
    for name, x in d.items():
        v[_I[name]] = x
    return v


def _dG_ddir(X: list[float], dirs: list[list[float]]) -> list[float]:
    r = cb._core.fanno_channel_flow(
        2e5, 600.0, X, 1.5e5, False, 0.5, 0.02, 0.0, "haaland", 1.0, -1.0, dirs
    )
    return list(r.dG_ddir)


@pytest.mark.parametrize(
    ("trace", "direction"),
    [
        ({}, {"CO2": 0.5, "H2O": -0.5, "O2": 0.5, "N2": -0.5}),
        ({"CO2": 1e-7, "H2O": 1e-7}, {"CO2": 0.5, "H2O": -0.5, "O2": 0.5, "N2": -0.5}),
        ({"CO2": 1e-7, "H2O": 1e-7}, {"CO2": 0.5, "H2O": -0.5}),
        ({"H2O": 0.9e-5}, {"CO2": 0.5, "H2O": -0.5}),
    ],
    ids=["both-absent", "both-trace", "pure-trace-trade", "co2-absent-h2o-trace"],
)
def test_a_direction_no_stencil_can_step_takes_its_linear_value(
    trace: dict[str, float], direction: dict[str, float]
) -> None:
    X = _gas(trace)
    got = _dG_ddir(X, [_vec(direction)])[0]
    parts = _dG_ddir(X, [_vec({name: 1.0}) for name in direction])
    ref = sum(direction[name] * p for name, p in zip(direction, parts, strict=True))
    assert abs(ref) > 10.0  # an order-one entry, not a vanishing one
    assert got == pytest.approx(ref, rel=1e-7)


def test_no_step_where_a_trace_crosses_the_stencil() -> None:
    """dG/dv is continuous as the H2O trace rises past the one-sided step
    (0.9e-5 split, 1.2e-5 one-sided): on main it went 0 -> 58."""
    v = _vec({"CO2": 0.5, "H2O": -0.5})
    below = _dG_ddir(_gas({"H2O": 0.9e-5}), [v])[0]
    above = _dG_ddir(_gas({"H2O": 1.2e-5}), [v])[0]
    assert below == pytest.approx(above, rel=1e-5)


def _feasible_central(trace: float, direction: dict[str, float], h: float) -> float:
    """dG along v by a plain central difference at a state where both traced
    species sit at `trace`, so every step of size h stays feasible."""
    X0 = _gas({name: trace for name in direction if name in ("CO2", "H2O")})
    Y0 = np.asarray(cb.mole_to_mass(X0))
    v = np.asarray(_vec(direction))

    def G(s: float) -> float:
        X = list(cb.mass_to_mole(list(Y0 + s * v)))
        return cb._core.fanno_channel_flow(
            2e5, 600.0, X, 1.5e5, False, 0.5, 0.02, 0.0, "haaland", 1.0, -1.0
        ).G

    return (G(h) - G(-h)) / (2.0 * h)


@pytest.mark.parametrize(
    "direction",
    [{"CO2": 0.5, "H2O": -0.5}, {"CO2": 0.5, "H2O": -0.5, "O2": 0.5, "N2": -0.5}],
    ids=["pure-trade", "four-species"],
)
def test_the_split_matches_the_limit_of_feasible_differences(direction: dict[str, float]) -> None:
    """An independent reference, not d_ddir: plain central differences along v
    itself at traces where the step is feasible (1e-3, 5e-4, 2.5e-4), taken to
    trace 0 by quadratic extrapolation. The split at both species absent must
    meet that one-sided (Y >= 0) limit."""
    ts = np.array([1e-3, 5e-4, 2.5e-4])
    vals = np.array([_feasible_central(t, direction, 0.25 * t) for t in ts])
    limit = float(np.polyval(np.polyfit(ts, vals, 2), 0.0))
    got = _dG_ddir(_gas({}), [_vec(direction)])[0]
    assert got == pytest.approx(limit, rel=1e-6)


def test_a_negative_fraction_does_not_recurse() -> None:
    """A state off the simplex (a -1e-20 fraction, only from a direct kernel
    call: the network clamps Y to [0, 1]) blocks every direction, so no split
    shrinks; it keeps the 0 it had instead of recursing until the process
    dies. Run in a subprocess so a regression fails here rather than killing
    the test session."""
    code = (
        "import combaero as cb\n"
        "N=[cb.species_name(i) for i in range(cb.num_species())]\n"
        "X=[0.0]*len(N)\n"
        "for k,x in {'N2':0.78084,'O2':0.20946,'AR':0.00934}.items(): X[N.index(k)]=x\n"
        "X[N.index('CO2')]=-1e-20\n"
        "v=[0.0]*len(N); v[N.index('O2')]=0.5; v[N.index('N2')]=-0.5\n"
        "r=cb._core.fanno_channel_flow(2e5,600.0,X,1.5e5,False,0.5,0.02,0.0,'haaland',1.0,-1.0,[v])\n"
        "print(r.dG_ddir[0])\n"
    )
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    assert done.returncode == 0, (done.returncode, done.stderr[-500:])
    assert float(done.stdout.split()[-1]) == 0.0
