"""Area change closures and element (#731): Idelchik (1966) diagrams 4-1, 3-9 and
5-2, carried to Mach by their control volume (include/area_change.h).

Fidelity is judged on Idelchik's own printed tables (validation/area_change,
signed off 2026-10-25); the solver properties on the closures themselves.
"""

import csv
import math
from collections import defaultdict

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from validation.area_change import DATA_ROOT

IDELCHIK = DATA_ROOT / "idelchik1966"
Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _table(name: str, row: str, col: str, val: str) -> dict[float, dict[float, float]]:
    out: dict[float, dict[float, float]] = defaultdict(dict)
    with (IDELCHIK / name).open() as f:
        for r in csv.DictReader(f):
            out[float(r[row])][float(r[col])] = float(r[val])
    return dict(out)


# ---------------------------------------------------------------------------
# Fidelity: the closures against the printed tables
# ---------------------------------------------------------------------------


def test_sudden_expansion_matches_diagram_4_1() -> None:
    """Every printed cell of 4-1 part 2 (10 < Re <= 3500), 78 cells, but three:

    - n 0.6 and 0.5 at Re 2000, lowered so the loss rises with the flow (the
      printed n 0.6 column falls between Re 2000 and 3000; see area_change.h):
      0.41 for 0.50, 0.61 for 0.65;
    - n 0.3 at Re 3500, printed 0.50 for Borda-Carnot's 0.49.
    Part 1 (Re >= 3500): Borda-Carnot exactly. Below Re 8: 26/Re.
    """
    tab = _table("diag4-1_sudden_expansion_low_re.csv", "F0_F2", "Re", "zeta")
    moved = {}
    for n, row in tab.items():
        for Re, zeta in row.items():
            z, _ = cb.sudden_expansion_zeta(Re, n)
            if abs(z - zeta) > 1e-12:
                moved[(n, Re)] = round(z, 3)
    assert moved == {(0.6, 2000.0): 0.407, (0.5, 2000.0): 0.61, (0.3, 3500.0): 0.49}
    for n in (0.0001, 0.1, 0.33, 0.6, 0.9):
        for Re in (3500.0, 1e4, 1e7):
            assert cb.sudden_expansion_zeta(Re, n)[0] == pytest.approx((1 - n) ** 2, rel=1e-14)
    for Re in (0.01, 1.0, 7.9):
        for n in (0.1, 0.4, 0.6):
            assert cb.sudden_expansion_zeta(Re, n)[0] == pytest.approx(26.0 / Re, rel=1e-14)


def test_conical_diffuser_matches_diagram_5_2() -> None:
    """Every printed cell of 5-2 (11 area ratios by 14 full angles), exactly."""
    tab = _table("diag5-2_conical_diffuser.csv", "F0_F1", "alpha_deg", "zeta_exp")
    n_cells = 0
    for n, row in tab.items():
        for alpha, zeta in row.items():
            assert cb.conical_diffuser_zeta(alpha, n) == pytest.approx(zeta, abs=1e-12), (n, alpha)
            n_cells += 1
    assert n_cells == 154


def test_sudden_contraction_matches_diagrams_3_10_and_3_9() -> None:
    """Every printed cell of 3-10 (Karev), 6 ratios by 13 Reynolds numbers,
    exactly: the last printed column is 3-9's 0.5 (1 - n). Of the two columns
    headed 10^4 the other is not used; it differs at n 0.1 only (0.50 against
    0.45). Above Re 1e4: 0.5 (1 - n). Below Re 8: 27/Re."""
    rows = list(csv.DictReader((IDELCHIK / "diag3-10_sudden_contraction_low_re.csv").open()))
    n_cells = 0
    for r in rows:
        n, Re, zeta = float(r["F0_F1"]), float(r["Re"]), float(r["zeta"])
        if r["column"] == "13":
            continue
        assert cb.sudden_contraction_zeta(Re, n)[0] == pytest.approx(zeta, abs=1e-12), (n, Re)
        n_cells += 1
    assert n_cells == 78
    for n in (0.0001, 0.3, 0.9):
        for Re in (1e4, 1e6):
            assert cb.sudden_contraction_zeta(Re, n)[0] == pytest.approx(0.5 * (1 - n), rel=1e-14)
    for Re in (0.5, 7.9):
        for n in (0.1, 0.6):
            assert cb.sudden_contraction_zeta(Re, n)[0] == pytest.approx(27.0 / Re, rel=1e-14)


def test_conical_nozzle_matches_diagrams_3_7_and_3_6() -> None:
    """zeta' (1 - n) with zeta' at every printed cell of 3-6 exactly."""
    rows = list(csv.DictReader((IDELCHIK / "diag3-6_conical_bellmouth_end_wall.csv").open()))
    assert len(rows) == 54
    for r in rows:
        a, lD, zp = float(r["alpha_deg"]), float(r["l_Dh"]), float(r["zeta_prime"])
        for n in (0.0, 0.25):
            assert cb.conical_nozzle_zeta(a, lD, n) == pytest.approx(zp * (1 - n), abs=1e-12)


# ---------------------------------------------------------------------------
# Solver properties of the 4-1 interpolation
# ---------------------------------------------------------------------------


def _exponent(Re: float, n: float, f=cb.sudden_expansion_zeta) -> float:
    """d ln(zeta Re^2) / d ln Re: the loss's exponent in the flow."""
    z, dz = f(Re, n)
    return 2.0 + Re * dz / z


def test_the_contraction_loss_rises_with_the_flow_everywhere() -> None:
    """3-10 needs no regularisation: every printed secant exponent is 0.95 or
    more. Measured 2026-10-25: the smallest local exponent of the
    interpolation is 0.857."""
    worst = min(
        _exponent(Re, n, cb.sudden_contraction_zeta)
        for n in (0.01, 0.1, 0.15, 0.3, 0.5, 0.6, 0.8, 0.99)
        for Re in np.geomspace(1.0, 2e4, 6000)
    )
    assert 0.8 < worst < 0.9, worst


@pytest.mark.parametrize("n", [0.1, 0.35, 0.6, 0.8])
def test_the_contraction_zeta_is_c1_in_re(n: float) -> None:
    for Re0 in (8.0, 10.0, 20.0, 100.0, 2000.0, 4000.0, 5000.0, 1e4):
        lo = cb.sudden_contraction_zeta(Re0 * (1 - 1e-10), n)
        hi = cb.sudden_contraction_zeta(Re0 * (1 + 1e-10), n)
        assert hi[0] == pytest.approx(lo[0], rel=1e-8), Re0
        assert hi[1] == pytest.approx(lo[1], rel=1e-6, abs=1e-12), Re0


def test_the_expansion_loss_rises_with_the_flow_everywhere() -> None:
    """zeta Re^2 must rise with Re, or one drop has two flows. Measured
    2026-10-25: the smallest exponent is 0.058 (n >= 0.6 near Re 3130),
    where the steepest printed drop turns onto Borda-Carnot's slope 2; a
    cubic cannot do better than ~0.067 there."""
    worst = min(
        _exponent(Re, n)
        for n in (0.01, 0.1, 0.15, 0.2, 0.3, 0.45, 0.5, 0.55, 0.6, 0.8, 0.99)
        for Re in np.geomspace(1.0, 5000.0, 6000)
    )
    assert 0.05 < worst < 0.07, worst


@pytest.mark.parametrize("n", [0.1, 0.25, 0.5, 0.6, 0.75])
def test_the_expansion_zeta_is_c1_in_re(n: float) -> None:
    """Value and Re-slope continuous at the laminar join (Re 8), every knot
    and the Borda-Carnot join (Re 3500)."""
    for Re0 in (8.0, 10.0, 15.0, 100.0, 1000.0, 2000.0, 3000.0, 3500.0):
        lo = cb.sudden_expansion_zeta(Re0 * (1 - 1e-10), n)
        hi = cb.sudden_expansion_zeta(Re0 * (1 + 1e-10), n)
        assert hi[0] == pytest.approx(lo[0], rel=1e-8), Re0
        assert hi[1] == pytest.approx(lo[1], rel=1e-6, abs=1e-12), Re0


def test_the_expansion_dzeta_dre_is_its_slope() -> None:
    for n in (0.2, 0.6):
        for Re in (9.0, 37.0, 640.0, 2500.0, 3300.0):
            h = 1e-6 * Re
            fd = (
                cb.sudden_expansion_zeta(Re + h, n)[0] - cb.sudden_expansion_zeta(Re - h, n)[0]
            ) / (2 * h)
            assert cb.sudden_expansion_zeta(Re, n)[1] == pytest.approx(fd, rel=1e-6)


# ---------------------------------------------------------------------------
# Compressible Borda-Carnot
# ---------------------------------------------------------------------------


def test_compressible_borda_carnot_low_mach_limit() -> None:
    """1 - Pt2/Pt1 -> (1 - n)^2 gamma M1^2 / 2."""
    for n in (0.1, 0.5, 0.9):
        s = cb.compressible_borda_carnot(1e-4, n, 1.4)
        assert 1.0 - s.pt_ratio == pytest.approx((1 - n) ** 2 * 0.7e-8, rel=1e-6)
        assert pytest.approx(n * 1e-4, rel=1e-6) == s.M2


def test_compressible_borda_carnot_partials() -> None:
    for M1, n, g in ((0.2, 0.3, 1.4), (0.7, 0.6, 1.33), (0.95, 0.2, 1.3)):
        s = cb.compressible_borda_carnot(M1, n, g)
        h = 1e-7
        fd_M = (
            cb.compressible_borda_carnot(M1 + h, n, g).pt_ratio
            - cb.compressible_borda_carnot(M1 - h, n, g).pt_ratio
        ) / (2 * h)
        fd_g = (
            cb.compressible_borda_carnot(M1, n, g + h).pt_ratio
            - cb.compressible_borda_carnot(M1, n, g - h).pt_ratio
        ) / (2 * h)
        assert s.dpt_dM1 == pytest.approx(fd_M, rel=1e-6)
        assert s.dpt_dgamma == pytest.approx(fd_g, rel=1e-5, abs=1e-10)


def _thermally_perfect_step(Pt: float, Tt: float, G1: float, n: float) -> tuple[float, float]:
    """(M1, 1 - Pt2/Pt1) from the control volume with cp(T): isentropic
    from (Pt, Tt) to flux G1, then mass, momentum and h0 to F1 / n."""
    from scipy.optimize import brentq

    X = cb.species.dry_air()
    R = cb.specific_gas_constant(X)
    h0 = cb.h_mass(Tt, X)

    def state(T: float) -> tuple[float, float]:
        P = Pt * math.exp((cb.s_mass(T, X, 101325.0) - cb.s_mass(Tt, X, 101325.0)) / R)
        return P, math.sqrt(max(2.0 * (h0 - cb.h_mass(T, X)), 0.0))

    def flux(T: float) -> float:
        P, u = state(T)
        return P / (R * T) * u

    T1 = brentq(lambda T: flux(T) - G1, 0.84 * Tt, Tt * (1 - 1e-12))
    P1, u1 = state(T1)
    G2 = n * G1
    impulse = P1 + n * G1 * u1

    def T_of(u2: float) -> float:
        return cb.calc_T_from_h_mass(h0 - 0.5 * u2 * u2, X, Tt)

    u2 = brentq(lambda u: G2 * R * T_of(u) / u + G2 * u - impulse, 1e-6 * u1, u1)
    T2 = T_of(u2)
    Pt2 = G2 * R * T2 / u2 * math.exp((cb.s_mass(Tt, X, 101325.0) - cb.s_mass(T2, X, 101325.0)) / R)
    a1 = math.sqrt(cb.isentropic_expansion_coefficient(T1, X) * R * T1)
    return u1 / a1, 1.0 - Pt2 / Pt


@pytest.mark.parametrize("Tt", [300.0, 1500.0])
def test_frozen_gamma_against_the_thermally_perfect_control_volume(Tt: float) -> None:
    """The kernel freezes gamma at Tt across the step. Measured 2026-10-25,
    relative loss error up to M1 0.65: 3.3e-4 at 300 K, 2.4e-3 at 1500 K."""
    worst = 0.0
    for n in (0.2, 0.5, 0.8):
        for G in (100.0, 200.0, 300.0, 400.0):
            G1 = G * math.sqrt(300.0 / Tt)
            M1, loss_tp = _thermally_perfect_step(2e5, Tt, G1, n)
            F_s = 1e-3
            k = cb.area_change_loss(G1 * F_s, 2e5, Tt, Y_AIR, F_s, F_s / n)
            assert k.mach == pytest.approx(M1, rel=2e-3)
            worst = max(worst, abs(k.dPt / 2e5 / loss_tp - 1.0))
    assert worst < (5e-4 if Tt < 1000 else 3e-3), worst


def test_mach_raises_the_expansion_loss() -> None:
    """At the same flow the compressible loss exceeds the incompressible one,
    and the gap closes as M^2. Measured 2026-10-25: 1.62x Borda-Carnot at
    M1 0.62, n 0.8 (the removed k1 = 1 + 0.35 M^2 gave 1.13x)."""
    F_s, n = 1e-3, 0.8
    gaps = []
    for G in (20.0, 40.0, 80.0):
        c = cb.area_change_loss(G * F_s, 2e5, 300.0, Y_AIR, F_s, F_s / n)
        i = cb.area_change_loss(G * F_s, 2e5, 300.0, Y_AIR, F_s, F_s / n, incompressible=True)
        gaps.append((c.dPt / i.dPt - 1.0) / c.mach**2)
    assert all(g > 0 for g in gaps)
    assert gaps[0] == pytest.approx(gaps[-1], rel=0.05)
    hi = cb.area_change_loss(380.0 * F_s, 2e5, 300.0, Y_AIR, F_s, F_s / n)
    assert hi.zeta / (1 - n) ** 2 > 1.5


# ---------------------------------------------------------------------------
# The element's flow: what the solver relies on
# ---------------------------------------------------------------------------

_CASES = {
    "sharp expansion": {"F_in": 1e-3, "F_out": 4e-3, "model": "sharp"},
    "sharp contraction": {"F_in": 4e-3, "F_out": 1e-3, "model": "sharp"},
    "conical diffuser": {"F_in": 1e-3, "F_out": 4e-3, "model": "conical", "length": 0.05},
    "conical nozzle": {"F_in": 4e-3, "F_out": 1e-3, "model": "conical", "length": 0.05},
    "small expansion (low Re)": {"F_in": 7.85e-7, "F_out": 3.14e-6, "model": "sharp"},
    "small contraction (low Re)": {"F_in": 3.14e-6, "F_out": 7.85e-7, "model": "sharp"},
}


@pytest.mark.parametrize("case", list(_CASES))
@pytest.mark.parametrize("incompressible", [False, True])
def test_the_flow_rises_with_the_drive_and_saturates(case: str, incompressible: bool) -> None:
    """m(drive) is odd, rises strictly, its slope is the analytic one, and
    (compressible) it saturates at the sonic flow with zero slope there."""
    kw = _CASES[case]
    drives = np.concatenate([np.geomspace(1e-6, 1.2e5, 400)])
    m = [
        cb.area_change_flow(2e5, 300.0, Y_AIR, 2e5 - d, incompressible=incompressible, **kw)
        for d in drives
    ]
    mdot = np.array([f.m_dot for f in m])
    slope = np.array([-f.dm_dPt_target for f in m])
    choked = np.array([f.choked for f in m])
    assert np.all(np.diff(mdot[~choked]) > 0.0)
    assert np.all(slope[~choked] > 0.0)
    if incompressible:
        assert not choked.any()
    elif "low Re" not in case:
        assert choked.any()
        assert np.all(np.diff(mdot[choked]) == 0.0)
        # C1 into the saturation: the slope just below choke is near zero.
        last = np.flatnonzero(~choked)[-1]
        assert slope[last] * drives[last] / mdot[last] < 0.05
    neg = cb.area_change_flow(2e5, 300.0, Y_AIR, 2e5 + 0.3, incompressible=incompressible, **kw)
    pos = cb.area_change_flow(2e5, 300.0, Y_AIR, 2e5 - 0.3, incompressible=incompressible, **kw)
    assert neg.m_dot == pytest.approx(-pos.m_dot, rel=1e-12)


def _floor(kw: dict) -> float:
    """A branch's drive floor: its incompressible zeta times 2 Pa of head."""
    z = cb.area_change_loss(1e-9, 2e5, 300.0, Y_AIR, incompressible=True, **kw).zeta
    if kw["F_in"] < kw["F_out"] and kw["model"] == "sharp":
        z = (1.0 - kw["F_in"] / kw["F_out"]) ** 2  # Borda-Carnot, not the low-Re zeta
    return 2.0 * z


@pytest.mark.parametrize("case", list(_CASES))
def test_the_element_flow_is_c1_through_zero_drive(case: str) -> None:
    """The two directions are different closures (an expansion one way, a
    contraction the other) fed from different nodes; their blend is still C1,
    and so is each branch at its drive floor, up to the O(M^2) gap between
    the floor's incompressible exponent and the compressible law's: a slope
    step of at most 1.1e-4 (measured 2026-10-25 over 0.2 mm - 0.3 m ducts,
    n 0.05 - 0.99, 1 - 20 bar)."""
    kw = _CASES[case]
    F_a, F_b = kw["F_in"], kw["F_out"]
    extra = {k: v for k, v in kw.items() if k not in ("F_in", "F_out")}
    rev = dict(kw, F_in=F_b, F_out=F_a)

    def flow(drive: float) -> tuple[float, float]:
        f = cb.area_change_element_flow(
            2e5, 300.0, Y_AIR, 2e5 - drive, 400.0, Y_AIR, F_a, F_b, **extra
        )
        return f.m_dot, -f.dm_dPt_b  # the drive moves Pt_b

    width = min(_floor(kw), _floor(rev))
    # The slope on either side of d0 differs by its curvature (2h over a
    # floor-wide scale) plus any step: measured 2026-10-25 the jump shrinks
    # with h everywhere (C1), and is at most ~7e-5 at a floor at this h.
    h = 3e-5 * width  # above the rounding of 2e5 Pa
    for d0, jump_tol, fd_tol in (
        (0.0, 1e-4, 2e-5),
        (_floor(kw), 2e-4, 1e-4),
        (-_floor(rev), 2e-4, 1e-4),
        (30.0 * width, 1e-4, 2e-5),
    ):
        lo, hi = flow(d0 - h), flow(d0 + h)
        assert hi[1] == pytest.approx(lo[1], rel=jump_tol), d0
        fd = (hi[0] - lo[0]) / (2 * h)
        assert flow(d0)[1] == pytest.approx(fd, rel=fd_tol), d0
    m = np.array([flow(d)[0] for d in np.linspace(-40.0 * width, 40.0 * width, 801)])
    assert np.all(np.diff(m) > 0.0)


def test_equal_areas_are_refused_by_the_kernel_and_a_pass_through_in_the_element() -> None:
    """The closure has no loss to set a flow at n = 1; the element never runs
    it there: equal faces are no area change, a pass-through (#1023,
    test_neighbour_areas_1023.py)."""
    with pytest.raises(ValueError, match="no area change"):
        cb.area_change_flow(2e5, 300.0, Y_AIR, 1.9e5, 1e-3, 1e-3)
    el = cb.network.AreaChangeElement("ac", "a", "b", F0=0.01, F1=0.01)
    el.validate()
    assert el._pass_through and el.residual_scale_kind == "p"


# ---------------------------------------------------------------------------
# Incompressible standalone kernels (given rho, mu)
# ---------------------------------------------------------------------------


def test_expansion_forward_positive_dp() -> None:
    res = cb.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04)
    assert res.dP == pytest.approx(0.5625 * 0.1**2 / (2 * 1.2 * 0.01**2), rel=1e-6)


def test_contraction_less_than_expansion() -> None:
    res_exp = cb.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04)
    res_con = cb.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.04, F1=0.01)
    assert 0.0 < res_con.dP < res_exp.dP


def test_reverse_flow_negative_dp() -> None:
    res = cb.sharp_area_change(m_dot=-0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04)
    assert res.dP < 0


def test_equal_areas_no_loss() -> None:
    """Borda-Carnot at n = 1 is zero (the old Re factor left 0.07 heads)."""
    res = cb.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.02, F1=0.02)
    assert res.dP == 0.0


@pytest.mark.parametrize(
    "m_dot, F0, F1, label",
    [
        (0.1, 0.01, 0.04, "expansion_fwd"),
        (0.1, 0.04, 0.01, "contraction_fwd"),
        (-0.1, 0.01, 0.04, "expansion_rev"),
        (-0.1, 0.04, 0.01, "contraction_rev"),
        (3e-4, 0.01, 0.04, "expansion_low_re"),
    ],
)
def test_jacobian_fd_accuracy(m_dot: float, F0: float, F1: float, label: str) -> None:
    rho, mu = 1.2, 1.8e-5
    res = cb.sharp_area_change(m_dot=m_dot, rho=rho, mu=mu, F0=F0, F1=F1)
    for name, kw, h in (
        ("dS_dm", "m_dot", 1e-7 * max(abs(m_dot), 1e-3)),
        ("dS_drho", "rho", 1e-6),
        ("dS_dmu", "mu", 1e-11),
    ):
        base = {"m_dot": m_dot, "rho": rho, "mu": mu, "F0": F0, "F1": F1}
        p, m = dict(base), dict(base)
        p[kw] += h
        m[kw] -= h
        fd = (cb.sharp_area_change(**p).dP - cb.sharp_area_change(**m).dP) / (2 * h)
        np.testing.assert_allclose(
            getattr(res, name), fd, rtol=1e-4, atol=1e-12, err_msg=f"{label} {name}"
        )


@pytest.mark.parametrize("F0, F1", [(0.01, 0.04), (0.04, 0.01)])
def test_zero_crossing_solver_properties(F0: float, F1: float) -> None:
    rho, mu, m_scale = 1.2, 1.8e-5, 1e-4
    ms = np.linspace(-5 * m_scale, 5 * m_scale, 201)
    r = [cb.sharp_area_change(m, rho, mu, F0, F1, m_scale=m_scale) for m in ms]
    dPs = np.array([x.dP for x in r])
    dS = np.array([x.dS_dm for x in r])
    assert np.all(np.isfinite(dPs)) and np.all(np.isfinite(dS))
    assert np.all(ms * dPs >= -1e-15)
    assert np.all(dS >= 0)


def test_dh_sets_the_reynolds_number() -> None:
    """At low Re the loss depends on Re, so D_h matters; 0 is circular."""
    F0, F1, m = 0.01, 0.04, 3e-4
    auto = cb.sharp_area_change(m_dot=m, rho=1.2, mu=1.8e-5, F0=F0, F1=F1, D_h=0.0)
    manual = cb.sharp_area_change(
        m_dot=m, rho=1.2, mu=1.8e-5, F0=F0, F1=F1, D_h=float(np.sqrt(4 * F0 / np.pi))
    )
    other = cb.sharp_area_change(m_dot=m, rho=1.2, mu=1.8e-5, F0=F0, F1=F1, D_h=0.05)
    assert auto.dP == pytest.approx(manual.dP, rel=1e-12)
    assert other.dP != pytest.approx(auto.dP, rel=1e-3)


def test_zero_area_is_refused() -> None:
    with pytest.raises(ValueError):
        cb.sharp_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.0, F1=0.04)


def test_conical_diffuser_gt_nozzle() -> None:
    d = cb.conical_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04, length=0.1)
    z = cb.conical_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.04, F1=0.01, length=0.1)
    assert d.dP > z.dP > 0


def test_conical_longer_less_loss() -> None:
    short = cb.conical_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04, length=0.05)
    long_ = cb.conical_area_change(m_dot=0.1, rho=1.2, mu=1.8e-5, F0=0.01, F1=0.04, length=0.5)
    assert long_.dP < short.dP


def test_area_change_exports() -> None:
    for name in (
        "sharp_area_change",
        "conical_area_change",
        "sudden_expansion_zeta",
        "sudden_contraction_zeta",
        "conical_diffuser_zeta",
        "conical_nozzle_zeta",
        "compressible_borda_carnot",
        "area_change_loss",
        "area_change_flow",
        "area_change_element_flow",
        "AreaChangeResult",
        "AreaChangeLoss",
        "AreaChangeFlow",
        "AreaChangeElementFlow",
        "BordaCarnotStep",
    ):
        assert getattr(cb, name) is getattr(cb._core, name), name
    for gone in (
        "area_change_residuals_and_jacobian",
        "conical_area_change_residuals_and_jacobian",
    ):
        assert not hasattr(cb._solver_tools, gone)


def test_area_change_validation() -> None:
    cb.network.AreaChangeElement("ac1", "n1", "n2", F0=0.01, F1=0.02).validate()
    with pytest.raises(ValueError, match="invalid Upstream Area"):
        cb.network.AreaChangeElement("ac1", "n1", "n2", F0=0.0, F1=0.02).validate()
    with pytest.raises(ValueError, match="invalid Hydraulic Diameter"):
        cb.network.AreaChangeElement("ac1", "n1", "n2", F0=0.01, F1=0.02, D_h=-1.0).validate()
    with pytest.raises(ValueError, match="needs length"):
        cb.network.AreaChangeElement(
            "ac1", "n1", "n2", F0=0.01, F1=0.02, model_type="conical"
        ).validate()


# ---------------------------------------------------------------------------
# Network
# ---------------------------------------------------------------------------
_F_BIG = 0.25 * np.pi * 0.1**2
_F_SMALL = 0.25 * np.pi * 0.05**2


def _step_into_chamber(
    F0: float, F1: float, model_type: str = "sharp", regime: str = "compressible"
):
    """mfb -> chamber(F0) -> area change -> chamber(F1) -> channel -> boundary.

    Both chambers keep their dynamic head (Pt = P + q), as the node the GUI
    inserts between two elements does.
    """
    from combaero.network import (
        AreaChangeElement,
        ChannelElement,
        FlowNetwork,
        LosslessConnectionElement,
        MassFlowBoundary,
        MomentumChamberNode,
        NetworkSolver,
        PressureBoundary,
    )

    net = FlowNetwork()
    net.add_node(MassFlowBoundary("mfb", m_dot=0.06, Tt=300.0, Y=Y_AIR))
    net.add_node(MomentumChamberNode("up", area=F0))
    net.add_node(MomentumChamberNode("down", area=F1))
    net.add_node(PressureBoundary("pb", Pt=101325.0, Tt=300.0, Y=Y_AIR))
    net.add_element(LosslessConnectionElement("feed", "mfb", "up"))
    net.add_element(
        AreaChangeElement(
            "step", "up", "down", F0=F0, F1=F1, model_type=model_type, length=0.05, regime=regime
        )
    )
    D1 = float(np.sqrt(4.0 * F1 / np.pi))
    net.add_element(ChannelElement("duct", "down", "pb", length=0.5, diameter=D1))
    return NetworkSolver(net)


@pytest.mark.parametrize(
    ("F0", "F1", "model_type"),
    [
        (_F_BIG, _F_SMALL, "sharp"),
        (_F_SMALL, _F_BIG, "sharp"),
        (_F_BIG, _F_SMALL, "conical"),
        (_F_SMALL, _F_BIG, "conical"),
    ],
    ids=["sharp-contraction", "sharp-expansion", "conical-contraction", "conical-expansion"],
)
def test_step_into_a_chamber_loses_exactly_its_loss(F0: float, F1: float, model_type: str) -> None:
    """Pt drops across the step by the element's own loss, never rises (#541)."""
    sol = _step_into_chamber(F0, F1, model_type).solve()
    assert sol["__success__"], sol["__message__"]
    q_down = sol["down.Pt"] - sol["down.P"]
    assert q_down > 10.0  # the chamber really carries a dynamic head
    step = sol["__element_diag__"]["step"]
    loss = sol["up.Pt"] - sol["down.Pt"]
    assert loss > 0.0
    assert loss == pytest.approx(step["dPt_loss"], rel=1e-8)


def test_step_into_a_chamber_jacobian_matches_fd() -> None:
    """Global analytic Jacobian vs FD, row-scaled (J_ij |x_j| against the row's
    largest term), 1e-6 relative step."""

    solver = _step_into_chamber(_F_BIG, _F_SMALL)
    sol = solver.solve()
    assert sol["__success__"]
    x = np.array([sol[n] for n in solver.unknown_names], dtype=float)
    _, jac = solver._residuals_and_jacobian(x)
    jac = jac.toarray()
    fd = fd_jacobian(solver, x, floor=1e-10)
    row_scale = np.abs(fd * x).max(axis=1, keepdims=True) + 1e-300
    err = np.abs(jac - fd) * np.abs(x) / row_scale
    assert err.max() < 1e-6, err.max()


@pytest.mark.parametrize(
    ("F0", "F1"), [(_F_BIG, _F_SMALL), (_F_SMALL, _F_BIG)], ids=["contraction", "expansion"]
)
def test_the_diagnostics_report_the_loss_that_was_solved(F0: float, F1: float) -> None:
    """zeta and dPt_loss are the solved Pt drop, at M_small ~0.4, and the
    reference block is the small section's (#732)."""
    solver = _step_into_chamber(F0, F1)
    solver.network.nodes["mfb"].m_dot = 0.3
    sol = solver.solve()
    assert sol["__success__"], sol["__message__"]
    step = sol["__element_diag__"]["step"]
    # The closure's Mach: the small section's (expansion), the vena
    # contracta's (contraction, ~0.66 here).
    assert 0.3 < step["mach_loss_ref"] < 0.7
    loss = sol["up.Pt"] - sol["down.Pt"]
    assert loss == pytest.approx(step["dPt_loss"], rel=1e-8)
    assert step["velocity"] == pytest.approx(0.3 / (step["rho_ref"] * _F_SMALL), rel=1e-12)


@pytest.mark.parametrize(("F0", "F1"), [(_F_BIG, _F_SMALL), (_F_SMALL, _F_BIG)])
def test_the_incompressible_regime_is_the_low_mach_limit(F0: float, F1: float) -> None:
    """Same network, same flow: the step's loss in the two regimes, and the
    gap closes as M^2 (#588). Measured 2026-10-25 at 0.06 kg/s (M ~0.07)."""
    gaps = []
    for m in (0.03, 0.06):
        losses = {}
        for regime in ("compressible", "incompressible"):
            s = _step_into_chamber(F0, F1, regime=regime)
            s.network.nodes["mfb"].m_dot = m
            sol = s.solve()
            assert sol["__success__"], sol["__message__"]
            losses[regime] = sol["up.Pt"] - sol["down.Pt"]
            mach = sol["__element_diag__"]["step"]["mach_loss_ref"]
        gaps.append((losses["compressible"] / losses["incompressible"] - 1.0) / mach**2)
    assert gaps[0] == pytest.approx(gaps[1], rel=0.1)
    assert abs(gaps[1]) < 2.0


def test_a_contraction_chokes_and_the_network_still_solves() -> None:
    """Supply at 3 bar into the atmosphere through a 4:1 sharp contraction:
    the vena contracta chokes, the flow is the sonic flow and independent of
    the back pressure."""
    from combaero.network import AreaChangeElement, FlowNetwork, NetworkSolver, PressureBoundary

    flows = []
    for p_back in (1.0e5, 0.8e5):
        net = FlowNetwork()
        net.add_node(PressureBoundary("hi", Pt=3e5, Tt=300.0, Y=Y_AIR))
        net.add_node(PressureBoundary("lo", Pt=p_back, Tt=300.0, Y=Y_AIR))
        net.add_element(AreaChangeElement("c", "hi", "lo", F0=4e-3, F1=1e-3))
        sol = NetworkSolver(net).solve()
        assert sol["__success__"], sol["__message__"]
        assert sol["__element_diag__"]["c"]["choked"] == 1.0
        flows.append(sol["c.m_dot"])
    assert flows[0] == pytest.approx(flows[1], rel=1e-10)
    sonic = cb.area_change_flow(3e5, 300.0, Y_AIR, 0.0 + 1.0, 4e-3, 1e-3)
    assert flows[0] == pytest.approx(sonic.m_dot, rel=1e-10)


def test_reversed_flow_reads_the_downstream_node() -> None:
    """Flow from to_node to from_node is fed from to_node's state (#731 item 8):
    a hot to_node changes the reversed flow, a hot from_node does not."""
    from combaero.network import AreaChangeElement, FlowNetwork, NetworkSolver, PressureBoundary

    def m(T_from: float, T_to: float) -> float:
        net = FlowNetwork()
        net.add_node(PressureBoundary("a", Pt=1.00e5, Tt=T_from, Y=Y_AIR))
        net.add_node(PressureBoundary("b", Pt=1.02e5, Tt=T_to, Y=Y_AIR))
        net.add_element(AreaChangeElement("s", "a", "b", F0=1e-3, F1=4e-3))
        sol = NetworkSolver(net).solve()
        assert sol["__success__"], sol["__message__"]
        return sol["s.m_dot"]

    base = m(300.0, 300.0)
    assert base < 0.0
    assert m(600.0, 300.0) == pytest.approx(base, rel=1e-12)
    assert abs(m(300.0, 600.0) / base - 1.0) > 0.1


@pytest.mark.parametrize("drive", [-80.0, -0.5, 0.0, 0.2, 5.0, 2e3])
def test_the_element_flow_partials_inside_and_outside_the_blend(drive: float) -> None:
    """Every partial of area_change_element_flow against central differences,
    both nodes at different temperatures and compositions (a fuel-air
    mixture on one side), inside the forward/reverse blend and past it."""
    X_fuel = [0.0] * cb.num_species()
    for name, x in (("N2", 0.74), ("O2", 0.2), ("CH4", 0.06)):
        X_fuel[cb.species_index_from_name(name)] = x
    Y_a, Y_b = list(cb.mole_to_mass(X_fuel)), Y_AIR
    args = {"Pt_a": 2e5, "Tt_a": 320.0, "Y_a": Y_a, "Pt_b": 2e5 - drive, "Tt_b": 600.0, "Y_b": Y_b}
    geo = {"F_a": 1e-3, "F_b": 4e-3, "model": "sharp"}
    f = cb.area_change_element_flow(**args, **geo)

    def m(**kw: object) -> float:
        return cb.area_change_element_flow(**{**args, **kw}, **geo, with_derivatives=False).m_dot

    scale = abs(f.m_dot) + abs(f.dm_dPt_a) * 1.0
    for key, attr, h in (
        ("Pt_a", "dm_dPt_a", 1e-3),
        ("Pt_b", "dm_dPt_b", 1e-3),
        ("Tt_a", "dm_dTt_a", 1e-3),
        ("Tt_b", "dm_dTt_b", 1e-3),
    ):
        fd = (m(**{key: args[key] + h}) - m(**{key: args[key] - h})) / (2 * h)
        assert getattr(f, attr) == pytest.approx(fd, rel=1e-5, abs=1e-9 * scale), key
    for node, grad in (("Y_a", f.dm_dY_a), ("Y_b", f.dm_dY_b)):
        Y = args[node]
        for k in range(len(Y)):
            h = 1e-6
            yp = list(Y)
            yp[k] += h
            if Y[k] > h:
                ym = list(Y)
                ym[k] -= h
                fd = (m(**{node: yp}) - m(**{node: ym})) / (2 * h)
            else:
                fd = (m(**{node: yp}) - f.m_dot) / h
            assert grad[k] == pytest.approx(fd, rel=1e-4, abs=1e-6 * scale), (node, k)


@pytest.mark.parametrize("dp", [0.0, 0.5, -0.5, 5.0])
def test_a_dead_or_creeping_area_change_solves(dp: float) -> None:
    """Zero and near-zero drive between a hot and a cold node, through a
    plenum the step shares with a feed: the root sits inside the C1
    forward/reverse blend and the drive floor, where the flow is linear in
    the drive with a finite slope (no square-root cusp for Newton)."""
    from combaero.network import (
        AreaChangeElement,
        FlowNetwork,
        NetworkSolver,
        OrificeElement,
        PlenumNode,
        PressureBoundary,
    )

    net = FlowNetwork()
    net.add_node(PressureBoundary("hi", Pt=1.2e5, Tt=300.0, Y=Y_AIR))
    net.add_node(PressureBoundary("lo", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    net.add_node(PlenumNode("p"))
    net.add_element(OrificeElement("o1", "hi", "p", Cd=0.6, diameter=0.01, correlation="fixed"))
    net.add_element(OrificeElement("o2", "p", "lo", Cd=0.6, diameter=0.01, correlation="fixed"))
    sol = NetworkSolver(net).solve()
    assert sol["__success__"]
    p_mid = sol["p.P"]
    # A hot side boundary at the plenum's pressure, plus dp.
    net.add_node(PressureBoundary("side", Pt=p_mid + dp, Tt=600.0, Y=Y_AIR))
    net.add_element(AreaChangeElement("ac", "p", "side", F0=1e-3, F1=4e-3))
    sol = NetworkSolver(net).solve()
    assert sol["__success__"], sol["__message__"]
    m = sol["ac.m_dot"]
    if dp == 0.0:
        # Zero to the solve's pressure tolerance: 1e-10 of 1.1e5 Pa is 1e-5 Pa
        # of drive, 2e-8 kg/s through this step (measured 1.8e-8).
        assert abs(m) < 1e-5 * sol["o1.m_dot"]
    else:
        assert np.sign(m) == -np.sign(dp)


def test_conical_loss_re_uses_the_station_static_viscosity() -> None:
    """area_change_loss's conical Re takes mu at the closure station's static T.

    It took mu(Tt0) before (#782 audit), 4% low at M 0.62. Diagnostic only:
    the loss itself does not read Re, so dPt is unchanged.
    """
    X = list(cb.species.dry_air())
    Pt0, Tt0, F_in, F_out = 3e5, 900.0, 1e-3, 3e-3
    fl = cb.area_change_flow(Pt0, Tt0, Y_AIR, 2.8e5, F_in, F_out, "conical", 0.05, 0.0, False)
    loss = cb.area_change_loss(fl.m_dot, Pt0, Tt0, Y_AIR, F_in, F_out, "conical", 0.05, 0.0, False)
    assert loss.mach > 0.5
    assert loss.dPt == pytest.approx(Pt0 - 2.8e5, rel=1e-9)
    g = cb.isentropic_expansion_coefficient(Tt0, X)
    T_static = Tt0 / (1.0 + 0.5 * (g - 1.0) * loss.mach**2)
    G_D = fl.m_dot / F_in * math.sqrt(4.0 * F_in / math.pi)
    assert loss.Re == pytest.approx(G_D / cb.viscosity(T_static, Pt0, X), rel=1e-12)
    assert loss.Re > 1.03 * G_D / cb.viscosity(Tt0, Pt0, X)
