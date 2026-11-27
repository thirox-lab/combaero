"""NH3 decomposition, 2 NH3 = N2 + 3 H2, in the combustion equilibrium (#1314).

Rich, complete combustion leaves the unburnt NH3, and the reforming + water-gas
shift equilibrium passed it through: NH3 at 0.217 where equilibrium holds
1e-5, 273 K hot (NH3 0.40, O2 0.13, N2 0.47 at 300 K, 1 bar). The
decomposition now joins the joint solve, with N in the element balances:
wherever NH3 is in the burned state (#1314), and since #1241 also forming NH3
from N2 and H2 (no threshold).

The reference values are Cantera 3.2 HP (TP for the isothermal cases)
equilibria on combaero's own species set, NH3 included, and its NASA-9 data
(tmp/ct_ref.py's gas(exclude=()) on 2026-11-25). Cantera is not imported here.
Where Cantera keeps a 1e-6 trace of O2 (H2 + O2 = H2O equilibrium, outside
the model: no reaction consumes the softmin's O2 trace, and none forms O2 from
steam), T differs by up to 0.04 K; everywhere else by 2e-5 K or less.
"""

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


def _phi_mix(phi: float, fuel: dict[str, float], air: dict[str, float]) -> np.ndarray:
    return np.asarray(cb.set_equivalence_ratio_mole(phi, _vec(fuel), _vec(air)))


# (feed, T_in [K], P [Pa])
_HP_CASES = {
    # The issue's table rows (NH3/O2/N2 rich; CH4 + NH3 rich).
    "rich_1bar": (lambda: _vec({"NH3": 0.40, "O2": 0.13, "N2": 0.47}), 300.0, 1e5),
    "rich_20bar": (lambda: _vec({"NH3": 0.40, "O2": 0.13, "N2": 0.47}), 300.0, 2e6),
    "phi1.2_1bar": (lambda: _phi_mix(1.2, {"NH3": 1.0}, _DRY_AIR), 300.0, 1e5),
    "phi1.2_700K_20bar": (lambda: _phi_mix(1.2, {"NH3": 1.0}, _DRY_AIR), 700.0, 2e6),
    "ch4_nh3_1bar": (
        lambda: _vec({"CH4": 0.15, "NH3": 0.05, "O2": 0.17, "N2": 0.63}),
        300.0,
        1e5,
    ),
    "ch4_nh3_phi1.5_20bar": (
        lambda: _phi_mix(1.5, {"CH4": 0.5, "NH3": 0.5}, _DRY_AIR),
        300.0,
        2e6,
    ),
    # Very rich, cold, high pressure: the cracking cools to 850 and 600 K and
    # equilibrium keeps 2 % and 14 % NH3 -- the reaction is equilibrated, not
    # just cracked. O2/N2 air: no carbon, so no methanation (a hydrocarbon
    # absent from the burned state is not formed; Cantera forms 1e-4 CH4 from
    # dry air's CO2).
    "crack_phi6_700K_50bar": (lambda: _phi_mix(6.0, {"NH3": 1.0}, _O2_N2), 700.0, 5e6),
    "crack_phi10_600K_20bar": (lambda: _phi_mix(10.0, {"NH3": 1.0}, _O2_N2), 600.0, 2e6),
}

# Cantera: T [K] and mole fractions.
_HP_CANTERA = {
    "rich_1bar": (
        1369.1923342109114,
        {"NH3": 9.513179226e-06, "H2": 0.2677048125, "N2": 0.5275593173, "H2O": 0.204726357},
    ),
    "rich_20bar": (
        1369.4735838167196,
        {"NH3": 0.0001898349984, "H2": 0.2674826049, "N2": 0.5275642868, "H2O": 0.2047632733},
    ),
    "phi1.2_1bar": (
        1954.7622489525618,
        {
            "NH3": 2.382418409e-07,
            "H2": 0.05724069861,
            "N2": 0.649232794,
            "H2O": 0.2868850821,
            "CO": 0.0001154000295,
            "CO2": 0.0001310380528,
        },
    ),
    "phi1.2_700K_20bar": (
        2251.16945034672,
        {
            "NH3": 3.04487401e-06,
            "H2": 0.05722769891,
            "N2": 0.6492318896,
            "H2O": 0.2868941363,
            "CO": 0.0001283100732,
            "CO2": 0.0001181281981,
        },
    ),
    "ch4_nh3_1bar": (
        1572.766889354246,
        {
            "NH3": 2.901057528e-06,
            "H2": 0.182668975,
            "N2": 0.5550849151,
            "H2O": 0.1351241938,
            "CO": 0.1012257725,
            "CO2": 0.02589323388,
            "CH4": 8.633311529e-09,
        },
    ),
    "ch4_nh3_phi1.5_20bar": (
        1858.3165462805591,
        {
            "NH3": 1.233292106e-05,
            "H2": 0.09755933979,
            "N2": 0.6153405574,
            "H2O": 0.1960663718,
            "CO": 0.05623575549,
            "CO2": 0.02792694457,
            "CH4": 1.457704021e-08,
        },
    ),
    "crack_phi6_700K_50bar": (
        847.453643153442,
        {"NH3": 0.0199300939, "H2": 0.4862140633, "N2": 0.390634002, "H2O": 0.1032218408},
    ),
    "crack_phi10_600K_20bar": (
        600.7606626774398,
        {"NH3": 0.1438999397, "H2": 0.4608232942, "N2": 0.3200908546, "H2O": 0.07518591153},
    ),
}

# Cases where Cantera's O2 trace (outside the model) moves T and H2/H2O.
_O2_TRACE = {"phi1.2_1bar", "phi1.2_700K_20bar"}


@pytest.mark.parametrize("case", list(_HP_CASES))
def test_rich_nh3_equilibrium_matches_cantera(case: str) -> None:
    """On main every case kept the unburnt NH3 (0.02-0.65 of the moles) and
    came out 34-544 K hot. Measured: within 2e-5 K and 4e-9 in mole fraction
    of Cantera (NH3 within 1e-7 of itself), and within 0.04 K, 7e-6 and 2e-4
    where Cantera keeps its O2 trace."""
    make, T_in, P = _HP_CASES[case]
    T_ct, X_ct = _HP_CANTERA[case]
    r = _core.combustion_equilibrium(T_in, make(), P)
    assert r.converged
    loose = case in _O2_TRACE
    assert abs(r.state.T - T_ct) <= (0.05 if loose else 1e-4), (r.state.T, T_ct)
    for name, x_ct in X_ct.items():
        x = r.state.X[_NAMES.index(name)]
        assert abs(x - x_ct) <= (1e-5 if loose else 1e-8), (name, x, x_ct)
    nh3 = r.state.X[_NAMES.index("NH3")]
    assert nh3 == pytest.approx(X_ct["NH3"], rel=1e-3 if loose else 1e-6)


# Isothermal: reforming_equilibrium(T, X, P) against Cantera TP.
_TP_CASES = {
    "nh3_700K_50bar": (
        {"NH3": 0.6, "N2": 0.2, "H2O": 0.2},
        700.0,
        5e6,
        {"NH3": 0.08481549322, "H2": 0.4829854751, "N2": 0.296597095, "H2O": 0.1356019367},
    ),
    "nh3_ch4_800K_20bar": (
        {"NH3": 0.3, "CH4": 0.1, "H2O": 0.3, "N2": 0.3},
        800.0,
        2e6,
        {
            "NH3": 0.006856446639,
            "H2": 0.3386870915,
            "N2": 0.3449982014,
            "H2O": 0.2320305191,
            "CO": 3.745543199e-05,
            "CO2": 0.0001081543391,
            "CH4": 0.07728177858,
        },
    ),
}


@pytest.mark.parametrize("case", list(_TP_CASES))
def test_nh3_isothermal_equilibrium_matches_cantera(case: str) -> None:
    """A cracker's range (700-800 K, 20-50 bar), where equilibrium keeps NH3
    at 0.7-8 %: on main NH3 stayed at its feed. Measured: within 1e-12."""
    feed, T, P, X_ct = _TP_CASES[case]
    r = _core.reforming_equilibrium(T, _vec(feed), P)
    assert r.converged
    for name, x_ct in X_ct.items():
        assert r.state.X[_NAMES.index(name)] == pytest.approx(x_ct, rel=1e-8, abs=1e-12), name


def test_lean_nh3_burns_completely_as_before() -> None:
    """Lean, complete combustion burns all the NH3 (the phi = 1 softmin's
    trace underflows to 0 this far from phi 1), and the products are N2, O2,
    H2O, no NH3 and no H2 (Cantera: 1248.5972 K): no reducing species, so the
    equilibrium stage is skipped, bit for bit as before #1241."""
    X = _vec({"NH3": 0.10, "O2": 0.19, "N2": 0.71})
    r = _core.combustion_equilibrium(300.0, X, 1e5)
    assert r.converged
    assert abs(r.state.T - 1248.5973) <= 1e-4
    assert r.state.X[_NAMES.index("NH3")] == 0.0
    assert r.state.X[_NAMES.index("H2")] == 0.0


def _g_over_RT(species: str, T: float) -> float:
    e = np.zeros(len(_NAMES))
    e[_NAMES.index(species)] = 1.0
    R = 8.31446261815324
    return (_core.h(T, e) - T * _core.s(T, e, 101325.0, 101325.0)) / (R * T)


@pytest.mark.parametrize("phi", [0.9993, 0.9996])
def test_the_lean_softmin_trace_of_nh3_is_decomposed_to_equilibrium(phi: float) -> None:
    """At phi 0.9993 the phi = 1 softmin leaves 8.5e-12 of the NH3 unburnt.
    Up to 0.12.0 it passed through below EQUILIBRIUM_NH3_MIN_FRACTION (1e-10)
    and was decomposed above it, a step in T and dX_out/dX at the cut. The
    decomposition is now always equilibrated (#1241): 6.4e-22 is left, at
    ln Q - ln K = 0 to 1e-9 (1e-14 measured; K from h and s of the pure
    species, independent of the solve), and the N and H it releases are
    conserved."""
    X = _phi_mix(phi, {"NH3": 1.0}, _O2_N2)
    burned = np.array(_core.complete_combustion_products_and_jacobian(X).X_b)
    i_nh3, i_n2, i_h2 = (_NAMES.index(k) for k in ("NH3", "N2", "H2"))
    assert 0.0 < burned[i_nh3] < 1e-8
    r = _core.combustion_equilibrium(300.0, X, 1e5)
    assert r.converged
    x = np.array(r.state.X)
    T = r.state.T
    assert x[i_nh3] < 1e-6 * burned[i_nh3]
    lnK = -(_g_over_RT("N2", T) + 3.0 * _g_over_RT("H2", T) - 2.0 * _g_over_RT("NH3", T))
    lnQ = (
        np.log(x[i_n2])
        + 3.0 * np.log(x[i_h2])
        - 2.0 * np.log(x[i_nh3])
        + 2.0 * np.log(1e5 / 101325.0)
    )
    assert abs(lnQ - lnK) < 1e-9
    # H released as H2 (no other H carrier moves: no carbon, O2 inert here).
    assert x[i_h2] == pytest.approx(
        1.5 * (burned[i_nh3] - x[i_nh3]) * x[i_n2] / burned[i_n2], rel=1e-6
    )


_ATOMS = {
    "N2": (0, 0, 0, 2),
    "O2": (0, 0, 2, 0),
    "AR": (0, 0, 0, 0),
    "CO2": (1, 0, 2, 0),
    "H2O": (0, 2, 1, 0),
    "CH4": (1, 4, 0, 0),
    "C2H6": (2, 6, 0, 0),
    "C3H8": (3, 8, 0, 0),
    "IC4H10": (4, 10, 0, 0),
    "NC5H12": (5, 12, 0, 0),
    "NC6H14": (6, 14, 0, 0),
    "NC7H16": (7, 16, 0, 0),
    "CO": (1, 0, 1, 0),
    "H2": (0, 2, 0, 0),
    "NH3": (0, 3, 0, 1),
}


@pytest.mark.parametrize("case", list(_HP_CASES))
def test_nh3_equilibrium_conserves_n_and_the_other_elements(case: str) -> None:
    """C, H, O and N per mole of burned gas, before and after the
    equilibrium: N moves between NH3 and N2 now."""
    make, T_in, P = _HP_CASES[case]
    X = make()
    burned = np.asarray(_core.complete_combustion(T_in, X, P).X)
    r = _core.combustion_equilibrium(T_in, X, P)
    A = np.array([_ATOMS[n] for n in _NAMES], dtype=float)
    # The moles change; the elements' shares of all atoms do not.
    b0 = burned @ A
    b1 = np.asarray(r.state.X) @ A
    assert np.allclose(b1 / b1.sum(), b0 / b0.sum(), rtol=1e-12, atol=0.0), (b0, b1)


# Per-column relative steps: the fuel, O2 and N2 are majors; the dry air's CO2
# (3e-4) and CH4 need larger ones (round-off at 1e-5 of CO2 is 1e-7).
_FD_COLUMNS = {
    "rich_1bar": {"NH3": 1e-4, "O2": 1e-4, "N2": 1e-4},
    "phi1.2_700K_20bar": {"NH3": 1e-4, "O2": 1e-4, "N2": 1e-4, "CO2": 1e-3},
    "ch4_nh3_1bar": {"NH3": 1e-4, "CH4": 1e-4, "O2": 1e-4, "N2": 1e-4},
    "ch4_nh3_phi1.5_20bar": {"NH3": 1e-4, "CH4": 1e-4, "O2": 1e-4, "CO2": 1e-3},
    "crack_phi6_700K_50bar": {"NH3": 1e-4, "O2": 1e-4, "N2": 1e-4},
}
_FD_PRODUCTS = ("NH3", "H2", "N2", "H2O", "CO", "CO2")


@pytest.mark.parametrize("case", list(_FD_COLUMNS))
def test_nh3_equilibrium_jacobian_matches_central_differences(case: str) -> None:
    """dT/dX and dXout/dX of the joint solve with NH3 decomposition, by the
    implicit function (the reaction's exact dF/dn0, dF/dT, dF/dP), against
    central differences; and dT/dP, dT/dT_in, dXout/dP. Measured: dT/dX
    within 1e-8 of itself, dXout/dX within 6e-7 of itself, dT/dP within 1e-7
    (a step of 1e-3 P: dT/dP is 2e-9 K/Pa at phi 1.2, and smaller steps are
    round-off)."""
    make, T, P = _HP_CASES[case]
    X = make()
    jac = _core.combustion_equilibrium_and_jacobian(T, X, P)
    dXo = np.asarray(jac.dXout_dX)

    def run(Xq: np.ndarray, Tq: float = T, Pq: float = P) -> tuple[float, np.ndarray]:
        r = _core.combustion_equilibrium(Tq, Xq, Pq)
        assert r.converged
        return r.state.T, np.asarray(r.state.X)

    for name, rel in _FD_COLUMNS[case].items():
        k = _NAMES.index(name)
        h = rel * X[k]
        Xp, Xm = X.copy(), X.copy()
        Xp[k] += h
        Xm[k] -= h
        (Tp, xp), (Tm, xm) = run(Xp), run(Xm)
        fd = (Tp - Tm) / (2.0 * h)
        assert abs(jac.dT_dX[k] - fd) <= 1e-6 * abs(fd), (name, jac.dT_dX[k], fd)
        fdx = (xp - xm) / (2.0 * h)
        scale = np.max(np.abs(fdx))
        for prod in _FD_PRODUCTS:
            j = _NAMES.index(prod)
            assert abs(dXo[j, k] - fdx[j]) <= 1e-5 * abs(fdx[j]) + 1e-7 * scale, (
                name,
                prod,
                dXo[j, k],
                fdx[j],
            )
    hP = 1e-3 * P
    (TPp, xPp), (TPm, xPm) = run(X, Pq=P + hP), run(X, Pq=P - hP)
    assert jac.dT_dP == pytest.approx((TPp - TPm) / (2.0 * hP), rel=1e-6)
    fdxP = (xPp - xPm) / (2.0 * hP)
    for prod in ("NH3", "H2"):
        j = _NAMES.index(prod)
        assert jac.dXout_dP[j] == pytest.approx(fdxP[j], rel=1e-5, abs=1e-7 * np.max(np.abs(fdxP)))
    hT = 1e-3
    fd_T = (run(X, Tq=T + hT)[0] - run(X, Tq=T - hT)[0]) / (2.0 * hT)
    assert jac.dT_dT_in == pytest.approx(fd_T, rel=1e-7)


# Very rich NH3 with a trace of carbon (#1314 review). With no hydrocarbon, the
# round cracks NH3 (making H2) and then runs the shift; the adiabatic search
# probes down to the 10 K floor, where CO2 at 1e-11 beside H2 leaves CO at
# 1e-185 from its bound, out of reach of solve_wgs_extent's 100 bisections:
# 10 of these 420 states threw "solve_wgs_extent: no convergence" (none on
# main, which did not crack the NH3; the shift solve alone threw there too).
@pytest.mark.parametrize("co2", [1e-12, 1e-9, 1e-6])
@pytest.mark.parametrize("phi", [3.0, 30.0, 60.0, 100.0])
def test_very_rich_nh3_with_a_carbon_trace_converges(phi: float, co2: float) -> None:
    X = _phi_mix(phi, {"NH3": 1.0}, {"O2": 0.21, "N2": 0.79, "CO2": co2})
    for T in (250.0, 550.0, 850.0):
        for P in (1e4, 1e5, 1e6, 1e7):
            r = _core.combustion_equilibrium(T, X, P)
            assert r.converged, (T, P)
            jac = _core.combustion_equilibrium_and_jacobian(T, X, P)
            assert np.isfinite(jac.dT_dX).all(), (T, P)


@pytest.mark.parametrize("T", [10.0, 20.0, 30.0])
def test_a_cold_shift_with_a_trace_of_co2_converges(T: float) -> None:
    """The shift alone at the thermo floor: H2-rich, no CO, CO2 1e-11. The
    root, CO at 1e-217 to 1e-79, sits that close to its bound; the second
    pass of solve_wgs_extent steps in the log of that distance and reaches
    it: CO = CO2 H2 / (H2O K), so a tenth of the CO2 gives a tenth of the CO
    (H2 and steam do not move), which a point that is not the root misses."""
    i_co = _NAMES.index("CO")
    co = []
    for x_co2 in (1e-11, 1e-12):
        X = _vec({"H2": 0.6, "H2O": 0.1, "CO2": x_co2, "N2": 0.3})
        r = _core.wgs_equilibrium(T, X, 1e5)
        assert r.state.X[_NAMES.index("CO2")] == pytest.approx(X[_NAMES.index("CO2")], rel=1e-9)
        co.append(r.state.X[i_co])
    assert 0.0 < co[1] < co[0] < 1e-70
    assert co[1] / co[0] == pytest.approx(0.1, rel=1e-9)
