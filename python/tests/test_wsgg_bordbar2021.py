"""#1071: the Bordbar et al. (2021) pressure-dependent WSGG gas emissivity.

FIDELITY on the paper's own data: its Figs. 3-5 slab results (1, 30 and
0.5 atm; 18 WSGG curves, 3609 points extracted from the PDF's vector paths)
re-solved with cb's gray gases (validation/radiation/bordbar2021_slab.py).
The other seven pressure sets are not tested by any figure; their constants
are pinned here by hand-summed printed rows instead.

Everything else here is the implementation's own contract: the analytic
partials against central differences, the ln P interpolation between the
sets (continuous at the nodes, the node itself exactly the set), the clamps.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

import combaero as cb
from combaero.network import GasRadiation

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.radiation import bordbar2021_slab as slab  # noqa: E402

SETS = (0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0, 80.0)
ATM = 101325.0


# --- the constants ------------------------------------------------------------


@pytest.mark.parametrize(
    ("p_atm", "i", "printed"),
    [
        # Table A.1 (0.1 atm), d_1k and d_4k; at Mr = 1 kappa_i is the row sum.
        (0.1, 0, (-1.017554e-3, 1.040575e-2, -3.644496e-2, 5.170700e-2, 3.165477e-1)),
        (0.1, 3, (5.193055, -5.536241e1, 2.238764e2, -4.350202e2, 8.549649e2)),
        # Table A.2 (0.5 atm), d_4k: its gas 4 is the WEAKEST band, so gray
        # gas i is not the same band from set to set.
        (0.5, 3, (-7.166366e-4, 7.658088e-3, -3.043043e-2, 5.603451e-2, 7.690453e-2)),
    ],
)
def test_printed_kappa_rows(p_atm: float, i: int, printed: tuple) -> None:
    g = cb.wsgg_bordbar2021_gray_gases(p_atm, 1500.0, 1.0)
    assert g.kappa[i] == pytest.approx(sum(printed), rel=1e-12)
    # At Mr = 2 the powers are Mr^(4-k): k = 0 multiplies Mr^4.
    g2 = cb.wsgg_bordbar2021_gray_gases(p_atm, 1500.0, 2.0)
    assert g2.kappa[i] == pytest.approx(sum(c * 2.0 ** (4 - k) for k, c in enumerate(printed)))


def test_gray_gases_refuse_outside_the_fit() -> None:
    for args in ((3.0, 1500.0, 1.0), (1.0, 299.0, 1.0), (1.0, 3001.0, 1.0), (1.0, 1500.0, 4.1)):
        with pytest.raises(ValueError):
            cb.wsgg_bordbar2021_gray_gases(*args)


# --- fidelity: the paper's own slab results ---------------------------------


@pytest.fixture(scope="module")
def slab_scores() -> list:
    return slab.score()


def test_reproduces_the_papers_slab_curves(slab_scores: list) -> None:
    """Max 1.25 %, mean 0.28 % of each curve's peak over 3609 points (no
    curve above 1.3 %). The residual is the extraction (axis calibration up
    to 0.37 units on spans of 100-2000) and the discretisation. What the
    bound catches: a misread exponent, reversed Mr or T powers, a wrong
    T_ref (each moves it by orders of magnitude) and a 5 % change in a
    dominant gray gas's weight (to 2-4 %); it does NOT catch a few-% change
    of a weak gas. The constants' transcription is checked separately: the
    parsed Appendix A equals the supplement's own coefficient files to
    2.5e-16 (validation/sources.yaml)."""
    assert len(slab_scores) == 18
    n = sum(s.n for s in slab_scores)
    assert n == 3609
    assert max(s.max_norm for s in slab_scores) < 0.0130
    assert sum(s.mean_norm * s.n for s in slab_scores) / n < 0.0030


def test_model_error_against_lbl_matches_table_2(slab_scores: list) -> None:
    """Not fidelity of ours: the paper's model-vs-LBL error, seen through our
    solve. Table 2's p = 30 atm Case 2 qr is 10.0 % max, 2.82 % mean."""
    s = next(s for s in slab_scores if (s.p_atm, s.case, s.quantity) == (30.0, 2, "qr"))
    assert s.lbl_max_norm == pytest.approx(0.100, abs=0.01)
    assert s.lbl_mean_norm == pytest.approx(0.0282, abs=0.005)


# --- across the sets --------------------------------------------------------


def _eps_at_set(p_atm: float, T: float, xh: float, xc: float, L: float, P_path: float) -> float:
    g = cb.wsgg_bordbar2021_gray_gases(p_atm, T, xh / xc)
    S = (xh + xc) * P_path / ATM * L
    return sum(a * (1.0 - math.exp(-k * S)) for a, k in zip(g.a, g.kappa, strict=True))


@pytest.mark.parametrize("p_atm", SETS)
def test_a_node_is_its_set(p_atm: float) -> None:
    """Eq. 7 at the tabulated pressure, exactly."""
    T, xh, xc, L = 1600.0, 0.08, 0.06, 0.2
    e = cb.gas_emissivity_bordbar2021(T, p_atm * ATM, xh, xc, L)
    assert e.eps == pytest.approx(_eps_at_set(p_atm, T, xh, xc, L, p_atm * ATM), rel=1e-12)


def test_between_sets_each_set_sees_the_actual_path() -> None:
    """eps interpolates linearly in ln P between the bracketing sets, both
    evaluated at the actual pressure path (never the coefficients)."""
    T, xh, xc, L = 1567.0, 0.0697, 0.0613, 0.126
    P = 5.0e5
    lo, hi = 2.0, 5.0
    w = math.log(P / ATM / lo) / math.log(hi / lo)
    want = (1 - w) * _eps_at_set(lo, T, xh, xc, L, P) + w * _eps_at_set(hi, T, xh, xc, L, P)
    assert cb.gas_emissivity_bordbar2021(T, P, xh, xc, L).eps == pytest.approx(want, rel=1e-12)


@pytest.mark.parametrize("p_atm", SETS[1:-1])
def test_continuous_at_the_nodes(p_atm: float) -> None:
    """eps is continuous across a node; d eps/dP may kink there."""
    args = (1700.0, None, 0.09, 0.07, 0.3)
    P = p_atm * ATM

    def eps(Pv: float) -> float:
        return cb.gas_emissivity_bordbar2021(args[0], Pv, *args[2:]).eps

    assert eps(P * (1 - 1e-9)) == pytest.approx(eps(P * (1 + 1e-9)), abs=1e-9)


def test_analytic_partials() -> None:
    rng = np.random.default_rng(1071)
    for _ in range(200):
        T = rng.uniform(320.0, 2950.0)
        xc = rng.uniform(0.01, 0.19)
        xh = xc * rng.uniform(0.11, 3.9)
        L = rng.uniform(0.01, 2.0)
        P = math.exp(rng.uniform(math.log(0.11), math.log(79.0))) * ATM
        if min(abs(math.log(P / ATM / p)) for p in SETS) < 1e-3:
            continue  # a kink
        x = [T, P, xh, xc, L]
        e = cb.gas_emissivity_bordbar2021(*x)
        for i, name in enumerate(("deps_dT", "deps_dP", "deps_dx_h2o", "deps_dx_co2", "deps_dl_b")):
            h = 1e-6 * x[i]
            hi, lo = list(x), list(x)
            hi[i] += h
            lo[i] -= h
            fd = (
                cb.gas_emissivity_bordbar2021(*hi).eps - cb.gas_emissivity_bordbar2021(*lo).eps
            ) / (2 * h)
            assert getattr(e, name) == pytest.approx(fd, rel=1e-5, abs=1e-9 / x[i]), name


def test_mixture_chains_through_the_mole_fractions() -> None:
    names = list(cb.species.names)
    Y = np.zeros(len(names))
    for sp, v in (("N2", 0.72), ("O2", 0.12), ("H2O", 0.07), ("CO2", 0.09)):
        Y[names.index(sp)] = v
    m = cb.gas_emissivity_bordbar2021_mixture(1700.0, 12e5, Y, 0.1)
    X = cb.mass_to_mole(list(Y))
    e = cb.gas_emissivity_bordbar2021(
        1700.0, 12e5, X[names.index("H2O")], X[names.index("CO2")], 0.1
    )
    assert m.eps == pytest.approx(e.eps, rel=1e-14)
    assert (m.deps_dT, m.deps_dP, m.deps_dl_b) == pytest.approx((e.deps_dT, e.deps_dP, e.deps_dl_b))
    for k in range(len(Y)):
        # Y >= 0 is required, so an absent species gets a forward difference
        # at a smaller step (its truncation error is first order).
        h = 1e-6 if Y[k] > 0.0 else 1e-9
        hi, lo = Y.copy(), Y.copy()
        hi[k] += h
        if Y[k] > 0.0:
            lo[k] -= h
        fd = (
            cb.gas_emissivity_bordbar2021_mixture(1700.0, 12e5, hi, 0.1).eps
            - cb.gas_emissivity_bordbar2021_mixture(1700.0, 12e5, lo, 0.1).eps
        ) / (hi[k] - lo[k])
        assert m.deps_dY[k] == pytest.approx(fd, rel=1e-5, abs=1e-7), names[k]


# --- clamps -------------------------------------------------------------------


def test_clamps() -> None:
    f = cb.gas_emissivity_bordbar2021
    # Temperature: the fitted 300-3000 K, flat outside.
    hot = f(3400.0, 10 * ATM, 0.1, 0.08, 0.3)
    assert hot.eps == f(3000.0, 10 * ATM, 0.1, 0.08, 0.3).eps
    assert hot.deps_dT == 0.0
    # Mr: hydrogen products (no CO2) are Mr = 4, dry CO2 is Mr = 0.1.
    assert f(1800.0, 10 * ATM, 0.2, 0.0, 0.3).eps == pytest.approx(
        f(1800.0, 10 * ATM, 0.16, 0.04, 0.3).eps
    )
    assert f(1800.0, 10 * ATM, 0.0, 0.2, 0.3).eps == pytest.approx(
        f(1800.0, 10 * ATM, 0.2 / 11.0, 2.0 / 11.0, 0.3).eps
    )
    # No absorber, no pressure or no path: no emission.
    assert f(1800.0, 10 * ATM, 0.0, 0.0, 0.3).eps == 0.0
    assert f(1800.0, 0.0, 0.1, 0.1, 0.3).eps == 0.0
    assert f(1800.0, 10 * ATM, 0.1, 0.1, 0.0).eps == 0.0
    # Optically thick: bounded by 1 though sum a_i reaches 1.039.
    for p in SETS:
        for T in (300.0, 800.0, 1500.0, 3000.0):
            e = f(T, p * ATM, 0.5, 0.5, 1000.0).eps
            assert 0.0 <= e <= 1.0
    with pytest.raises(ValueError):
        f(1800.0, 10 * ATM, 0.7, 0.4, 0.3)


# --- the network spec -------------------------------------------------------


def test_gas_radiation_spec() -> None:
    assert GasRadiation(beam_length=0.1).emissivity_model == "lefebvre"
    GasRadiation(beam_length=0.1, emissivity_model="bordbar2021")
    with pytest.raises(ValueError, match="emissivity_model"):
        GasRadiation(beam_length=0.1, emissivity_model="wsgg")
    with pytest.raises(ValueError, match="fuel_air_ratio"):
        GasRadiation(beam_length=0.1, emissivity_model="bordbar2021", fuel_air_ratio=0.03)
    with pytest.raises(ValueError, match="luminosity"):
        GasRadiation(beam_length=0.1, emissivity_model="bordbar2021", luminosity=1.7)


def test_dlr_operating_points() -> None:
    """ACCURACY context, cross-source (not fidelity): Bordbar against Eq.
    8.11 at the Behrendt 2008 rig's gas, beam 0.126 m. Eq. 8.11 is 3 %
    above at 5 bar and 32 % above at 10 bar."""
    from validation.cooling.effusion_network_runner import hot_gas

    Y = list(hot_gas()[0])
    q = cb.elemental_fuel_air_ratio(Y)[0]
    for P, T, want in ((5e5, 1567.0, 0.1325), (10e5, 1532.0, 0.1981)):
        e = cb.gas_emissivity_bordbar2021_mixture(T, P, Y, 0.126).eps
        assert e == pytest.approx(want, abs=5e-4)
        assert cb.gas_emissivity_luminous(T, P, q, 0.126).eps > e
