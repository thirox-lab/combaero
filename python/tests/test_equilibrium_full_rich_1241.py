"""Full equilibrium of rich products (#1241): methanation and NH3 synthesis always.

Up to 0.12.0 the reforming stage ran only where the burned state held more
than 1e-10 of a hydrocarbon, and NH3 decomposition only above 1e-10 of NH3.
Below the cuts the traces were frozen, and no CH4 or NH3 formed from CO, H2
and N2 at all: very rich syngas at pressure came out 320 K cold. Both
reactions now run wherever the burned state holds a reducing species (H2, CO,
a hydrocarbon or NH3) and the elements they need; a state without one (a
complete lean burn, a fuel-free mix) takes the frozen path as before, bit for
bit.

The reference temperatures are Cantera 3.2 HP equilibria on combaero's own
species set, NH3 included, and its NASA-9 data (2026-11-27). Cantera is not
imported here.
"""

from __future__ import annotations

import numpy as np
import pytest

import combaero as cb

_core = cb._core
_NAMES = [cb.species_name(i) for i in range(cb.num_species())]
_DRY_AIR = {"N2": 0.78084, "O2": 0.20946, "AR": 0.00934, "CO2": 0.00036}
_O2_N2 = {"O2": 0.21, "N2": 0.79}


def _vec(d: dict[str, float]) -> np.ndarray:
    x = np.zeros(len(_NAMES))
    for k, v in d.items():
        x[_NAMES.index(k)] += v
    return x / x.sum()


# (T_in [K], P [Pa], fuel, phi in dry air): Cantera T [K], X_CH4, X_NH3.
# Main (0.12.0) gave 903.39, 1622.17, 1561.18, 1565.62, 703.77 and 534.90 K.
_CANTERA = {
    "syngas_phi14": ((354.0, 32.7e5, {"H2": 0.5, "CO": 0.5}, 14.0), 1223.8548, 5.08474e-02, 3.56441e-04),
    "syngas_phi4": ((354.0, 32.7e5, {"H2": 0.5, "CO": 0.5}, 4.0), 1622.4924, 3.08759e-05, 8.30951e-05),
    "h2_phi4_100bar": ((300.0, 1e7, {"H2": 1.0}, 4.0), 1562.8290, 2.09934e-06, 1.04527e-03),
    "ch4_phi2_100bar": ((300.0, 1e7, {"CH4": 1.0}, 2.0), 1566.0592, 1.09165e-04, 2.76299e-04),
    "nh3_co_phi8": ((300.0, 2e6, {"NH3": 0.5, "CO": 0.5}, 8.0), 1070.7784, 6.78351e-02, 5.26253e-04),
    "nh3_co_phi16": ((300.0, 2e6, {"NH3": 0.5, "CO": 0.5}, 16.0), 1032.3665, 1.24427e-01, 5.43781e-04),
}  # fmt: skip


@pytest.mark.parametrize("case", list(_CANTERA))
def test_rich_products_match_cantera(case: str) -> None:
    """Methanation and NH3 synthesis at their equilibrium: T to 1e-3 K and
    CH4, NH3 to 1e-4 of themselves (4e-5 K and 2e-6 measured). The #1241
    issue's case, H2/CO at phi 14 and 32.7 bar, was 320 K cold (5 % CH4 not
    formed); NH3 with CO 367 to 497 K cold; H2 at phi 4 and 100 bar 1.65 K."""
    (T_in, P, fuel, phi), T_ct, x_ch4, x_nh3 = _CANTERA[case]
    X = np.asarray(cb.set_equivalence_ratio_mole(phi, _vec(fuel), _vec(_DRY_AIR)))
    r = _core.combustion_equilibrium(T_in, X, P)
    assert r.converged
    assert pytest.approx(T_ct, abs=1e-3) == r.state.T
    assert r.state.X[_NAMES.index("CH4")] == pytest.approx(x_ch4, rel=1e-4)
    assert r.state.X[_NAMES.index("NH3")] == pytest.approx(x_nh3, rel=1e-4)
    _core.combustion_equilibrium_and_jacobian(T_in, X, P)


@pytest.mark.parametrize("fuel", ["CH4", "C2H6", "C3H8", "NC7H16", "H2"])
@pytest.mark.parametrize("T", [650.0, 900.0])
@pytest.mark.parametrize("P", [2e5, 2e6])
@pytest.mark.parametrize("trace", [1e-20, 1e-16, 1e-13])
def test_a_warm_fuel_trace_in_air_converges(fuel: str, T: float, P: float, trace: float) -> None:
    """A fuel trace in warm air: the phi = 0 softmin leaves half of it
    unburnt, which now reforms. Newton stalled there (151 refused states,
    1e-20 to 1e-13 of fuel, 500-900 K) because the heavy alkanes'
    equilibrium amounts are no double (n-heptane about 1e-391 beside 1e-15 of
    CH4 at 700 K): the log step underflowed and no step was finite. A key
    whose equilibrium amount is below the normal range is boxed at 0, its
    reaction run to completion away from it (joint_key_unrepresentable)."""
    X = np.array(cb.species.dry_air())
    X[_NAMES.index(fuel)] += trace
    X /= X.sum()
    r = _core.combustion_equilibrium(T, X, P)
    assert r.converged
    _core.combustion_equilibrium_and_jacobian(T, X, P)


# (state, T_in [K], P [Pa], the directions a reactant arrives along): states
# that lack a species, or a whole element, the arriving reactant feeds.
_DIRECTIONAL = [
    ("dry air, fuel-free", _DRY_AIR, 1100.0, 1.2e4, ("CH4", "CO", "NH3", "NC7H16", "H2")),
    ("O2/N2, carbon-free", _O2_N2, 800.0, 7.5e4, ("NC7H16", "CH4", "CO")),
    ("rich H2", {"H2": 0.84, "O2": 0.034, "N2": 0.126}, 769.0, 6.5e5, ("CH4", "CO", "C3H8")),
    ("rich NH3", {"NH3": 0.83, "O2": 0.036, "N2": 0.134}, 476.0, 4.8e6, ("CO", "CH4")),
    ("CO, hydrogen-free", {"CO": 0.8, "O2": 0.043, "N2": 0.157}, 503.0, 1.47e4, ("NC7H16", "H2", "CH4")),
]  # fmt: skip


@pytest.mark.parametrize(
    ("state", "T", "P", "species"),
    [(s, T, P, sp) for lab, s, T, P, sps in _DIRECTIONAL for sp in sps],
    ids=[f"{lab}-{sp}" for lab, s, T, P, sps in _DIRECTIONAL for sp in sps],
)
def test_the_slope_along_an_arriving_reactant_matches_a_difference(
    state: dict[str, float], T: float, P: float, species: str
) -> None:
    """Where the burned state lacks a species the arriving reactant feeds,
    T has a one-sided slope along each direction, and in dry air no dT/dX
    matrix holds them (steam is a trace too, so the slopes along two fuels do
    not add: CH4 with CO at 335 K, 11 %). combustion_equilibrium_and_jacobian
    takes the direction (dX_dir), seeds its trace states along it, and J dX
    is that direction's slope: to 2e-5 of the best forward difference over
    1e-8 to 1e-11 (6e-6 at worst measured). Without dX_dir the matrix was off
    by up to 4x on the rich and hydrogen-free burns, and of the wrong sign for
    CO into rich NH3 (-18.7 against 816 K)."""
    x = _vec(state)
    e = np.zeros(len(_NAMES))
    e[_NAMES.index(species)] = 1.0
    d = e - x
    j = _core.combustion_equilibrium_and_jacobian(T, x, P, dX_dir=list(d))
    an = float(np.asarray(j.dT_dX) @ d)
    T0 = _core.combustion_equilibrium(T, x, P).state.T

    def fd(h: float) -> float:
        xx = x + h * d
        return (_core.combustion_equilibrium(T, xx / xx.sum(), P).state.T - T0) / h

    best = min(abs(fd(h) / an - 1.0) for h in (1e-8, 1e-9, 1e-10, 1e-11))
    assert best < 2e-5, (an, best)


def test_a_direction_that_feeds_nothing_absent_changes_nothing() -> None:
    """dX_dir only changes the derivatives where the direction feeds a
    species the burned state lacks: at an ordinary rich state the Jacobian
    is the matrix's, bit for bit."""
    x = np.asarray(cb.set_equivalence_ratio_mole(1.5, _vec({"CH4": 1.0}), _vec(_DRY_AIR)))
    d = _vec({"CH4": 1.0}) - x
    a = _core.combustion_equilibrium_and_jacobian(600.0, x, 5e5)
    b = _core.combustion_equilibrium_and_jacobian(600.0, x, 5e5, dX_dir=list(d))
    assert list(a.dT_dX) == list(b.dT_dX)
    assert [list(r) for r in a.dXout_dX] == [list(r) for r in b.dXout_dX]


@pytest.mark.parametrize("x_o2", [1e-12, 1e-11, 3e-11, 1e-10, 1.17808e-10, 3e-10, 1e-9, 1e-8])
def test_a_cold_rich_alkane_with_an_oxygen_trace_has_the_right_jacobian(x_o2: float) -> None:
    """#1339, found by the 0.12.0 release gate: rich isobutane at 250 K with
    an O2 trace of 1e-10 to 1e-9 left a burnt hydrocarbon of about 1e-10, at
    the reforming cut, and J took a branch F did not have (dT/dX_CO2 10.15
    against a central difference of 10.67, 5 %). Without the cut, the CO2
    column matches a central difference of 1e-3 X_k to 1e-4 over the issue's
    O2 sweep (9e-5 at worst, 3e-10; 1.6 % to 2.6 % on 0.12.0)."""
    X = _vec({"IC4H10": 0.99998222, "CO2": 1.7782478e-5, "N2": 4.3907e-10, "AR": 5.23e-12})
    X[_NAMES.index("CO")] = 5.62e-16
    X[_NAMES.index("O2")] = x_o2
    T, P = 250.0, 177827.94
    j = _core.combustion_equilibrium_and_jacobian(T, X, P)
    k = _NAMES.index("CO2")
    h = 1e-3 * X[k]

    def T_out(dx: float) -> float:
        xx = X.copy()
        xx[k] += dx
        return _core.combustion_equilibrium(T, xx, P).state.T

    fd = (T_out(h) - T_out(-h)) / (2.0 * h)
    assert j.dT_dX[k] == pytest.approx(fd, rel=1e-4)
