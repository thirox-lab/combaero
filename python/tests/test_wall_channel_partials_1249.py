"""Partials of wall-coupled channels (#1249, #1254, #1264, #1265).

* #1249: the Fanno-Rayleigh kernel's implicit derivatives were central
  differences of the march at a relative step of 1e-4, O(h^2): 1.9e-6 off in
  dG/dP_target (1.3e-6 in dG/dTt_exit) at an exit Mach of 0.83 on the
  hotcold_mild reproducer's 30 cm duct. Richardson's extrapolation of the
  steps h and h / 2 takes them to O(h^4) (1e-10 there), with less rounding
  at low drive than before; near the choke the plain difference stays.
* #1254: a pin-fin surface's dh_dT was the smooth correlation's, scaled to
  h_eff: 1.398 against 1.771 W/(m^2 K^2) on a 300 K air array (21% low),
  which reached the wall-coupled flow partial through the duct's static T
  (0.9%, 3.5e-5 scaled). It is now the array's own, through mu, k and Pr.
* #1264: a ribbed surface's dh_dmdot divided by the flow, 0 at exactly zero
  flow where h_r is linear in |m| with a finite slope.
* #1275: a ribbed surface's dh_dT left out the ribbed side (72-80% low).
* #1265: a one-stream mix (a wall-heated duct's own exit state) at exactly
  zero flow sat in the mixer's no-flow band and lost its partials, dropping
  the exit Tt's 1:1 move with the feed. A single stream has no 1 / m pole and
  keeps them.
"""

from __future__ import annotations

import math

import pytest
from _jacobian_sweep import Case, Report, evaluate, point
from test_mixer_net_reverse_1224 import _Y, AIR, _analytic, _mix

import combaero as cb
from combaero.network.components import ConvectiveSurface, PinFinModel, RibbedModel

# The #1249 reproducer's cold duct (hotcold_mild:0/v1, start vector).
X_1249 = [
    0.7261700834721987,
    0.194841998575084,
    0.008649310676602776,
    0.000372013362434528,
    0.0,
    0.06996659391367997,
    *([0.0] * 9),
]
ARGS_1249 = (1e6, 260.0, 260.00625083968873, X_1249, 907299.0788979029)


def _rayleigh(Pt0: float, Tt0: float, Tte: float, Pd: float, derivs: bool = False) -> object:
    return cb.fanno_rayleigh_channel_flow(
        Pt0, Tt0, Tte, X_1249, Pd, True, 0.3, 0.012, 0.0, "haaland", 1.0, -1.0, [], derivs
    )


def test_1249_rayleigh_partials_match_a_converged_difference() -> None:
    Pt0, Tt0, Tte, _, Pd = ARGS_1249
    base = [Pt0, Tt0, Tte, Pd]
    r = _rayleigh(*base, derivs=True)
    assert 0.8 < r.M_exit < 0.9 and not r.choked
    for k, name in enumerate(("dG_dPt0", "dG_dTt0", "dG_dTt_exit", "dG_dP_target")):
        h = 1e-6 * abs(base[k])
        up, dn = list(base), list(base)
        up[k] += h
        dn[k] -= h
        fd = (_rayleigh(*up).G - _rayleigh(*dn).G) / (2.0 * h)
        # The O(h^2) difference was 1.9e-6 off in dG/dP_target.
        assert getattr(r, name) == pytest.approx(fd, rel=1e-7), name


def test_1249_rayleigh_partials_near_sonic_band() -> None:
    # A 3 m duct at exit Mach 0.977 (zeta 2.1e-3): the zeta check cuts the
    # step, and a plain central difference at the cut step carried ~1.2e-5
    # of truncation in dG/dP_target by construction. The reference is a
    # 4-level Richardson table of the march itself at a step well inside
    # the zeta-limited one.
    air = list(cb.species.dry_air())
    Pt0, Tt0, Te, P, L, D = 8e5, 400.0, 410.0, 232242.448, 3.0, 0.01
    r = cb.fanno_rayleigh_channel_flow(
        Pt0, Tt0, Te, air, P, False, L, D, 0.0, "haaland", 1.0, -1.0, [], True
    )
    G = r.G

    def E(g: float, p0: float = Pt0, t0: float = Tt0, te: float = Te) -> float:
        return cb.fanno_rayleigh_march(g, p0, t0, te, air, L, D, 0.0, "haaland", 1.0).P_exit

    zeta = cb.fanno_rayleigh_march(G, Pt0, Tt0, Te, air, L, D, 0.0, "haaland", 1.0).zeta_min
    assert 1e-3 < zeta < 1e-2 and 0.97 < r.M_exit < 0.98

    def rich(f: object, x: float, h: float) -> float:
        d = [(f(x + h / 2**i) - f(x - h / 2**i)) / (2.0 * h / 2**i) for i in range(4)]
        r1 = [(4.0 * d[i + 1] - d[i]) / 3.0 for i in range(3)]
        r2 = [(16.0 * r1[i + 1] - r1[i]) / 15.0 for i in range(2)]
        return (64.0 * r2[1] - r2[0]) / 63.0

    h = 0.005 * zeta
    RG = rich(E, G, h * G)
    ref = {
        "dG_dPt0": -rich(lambda v: E(G, p0=v), Pt0, h * Pt0) / RG,
        "dG_dTt0": -rich(lambda v: E(G, t0=v), Tt0, h * Tt0) / RG,
        "dG_dTt_exit": -rich(lambda v: E(G, te=v), Te, h * Te) / RG,
        "dG_dP_target": 1.0 / RG,
    }
    for name, v in ref.items():
        # The plain cut difference was 1.2e-5 off in dG/dP_target.
        assert getattr(r, name) == pytest.approx(v, rel=3e-7), name


def _pin_fin_surface(nu_set: object = None) -> ConvectiveSurface:
    model = PinFinModel(
        pin_diameter=0.005,
        S_D=2.5,
        X_D=2.5,
        H_D=1.0,
        N_rows=10,
        arrangement=cb.PinArrangement.Inline,
        k_pin=20.0,
        nu_set=nu_set,
    )
    return ConvectiveSurface(area=0.01, model=model)


AIR_X = [0.0] * len(AIR)
AIR_X[0], AIR_X[1], AIR_X[2] = 0.7808, 0.2095, 0.0097


@pytest.mark.parametrize(
    "nu_set",
    [None, cb.chyu_1998_nu(cb.PinArrangement.Inline, cb.PinNuSurface.Total)],
    ids=["metzger-no-pr", "chyu-pr-0.4"],
)
@pytest.mark.parametrize("T", [300.0, 900.0])
def test_1254_pin_fin_dh_dT_is_the_arrays_own(nu_set: object, T: float) -> None:
    # At fixed flow and P, the convention of every surface's dh_dT.
    surf = _pin_fin_surface(nu_set)
    P, D, m = 2.2e5, 0.009, 0.0244
    area = math.pi / 4.0 * D * D

    def h(Tx: float) -> object:
        rho = cb.density(Tx, P, AIR_X)
        return surf.htc_and_T(Tx, P, AIR_X, m / (area * rho), D, 0.02, flow_area=area)

    r = h(T)
    d = 1e-4 * T
    fd = (h(T + d).h - h(T - d).h) / (2.0 * d)
    assert r.dh_dT == pytest.approx(fd, rel=1e-6)


def test_pin_fin_nu_pr_slope() -> None:
    g = cb.PinFinGeometry(S_D=2.5, X_D=2.5, H_D=1.0, N_rows=10)
    s = cb.chyu_1998_nu(cb.PinArrangement.Staggered, cb.PinNuSurface.Total)
    r = cb.evaluate_pin_fin_nu(s, g, 2e4, 0.7)
    fd = (
        cb.evaluate_pin_fin_nu(s, g, 2e4, 0.7 + 1e-6).Nu
        - cb.evaluate_pin_fin_nu(s, g, 2e4, 0.7 - 1e-6).Nu
    ) / 2e-6
    assert r.dNu_dPr == pytest.approx(fd, rel=1e-7)
    assert r.dNu_dPr == pytest.approx(0.4 * r.Nu / 0.7, rel=1e-6)
    assert cb.evaluate_pin_fin_nu(cb.metzger_1986_staggered_nu(), g, 2e4, 0.7).dNu_dPr == 0.0


@pytest.mark.parametrize("v", [0.0, 1e-13, 1.0])
def test_1264_ribbed_dh_dmdot_is_the_one_sided_slope(v: float) -> None:
    # A real rib (the defaults e_D = p_e = 0 are degenerate: h ~ 7e7).
    surf = ConvectiveSurface(area=0.01, model=RibbedModel(e_D=0.0625, p_e=10.0))
    T, P, D = 400.0, 2e5, 0.02
    area = math.pi / 4.0 * D * D
    rho = cb.density(T, P, AIR_X)

    def r(vel: float) -> object:
        return surf._ribbed_result(T, P, AIR_X, vel, D, 0.2, math.nan, True)

    # One-sided (the forward side) at and next to zero flow, central above.
    if v < 1e-6:
        dm = 1e-9
        fd = (r(v + dm / (rho * area)).h - r(v).h) / dm
    else:
        dm = 1e-5 * rho * v * area
        fd = (r(v + dm / (rho * area)).h - r(v - dm / (rho * area)).h) / (2.0 * dm)
    assert r(v).dh_dmdot == pytest.approx(fd, rel=1e-5)


@pytest.mark.parametrize("m", [0.003, 0.03])
@pytest.mark.parametrize("T", [300.0, 900.0])
def test_1275_ribbed_dh_dT_includes_the_ribbed_side(T: float, m: float) -> None:
    # At fixed flow and P. Smooth side only, it was 72-80% low.
    surf = ConvectiveSurface(area=0.01, model=RibbedModel(e_D=0.0625, p_e=10.0))
    P, D = 2e5, 0.02
    area = math.pi / 4.0 * D * D

    def h(Tx: float) -> object:
        rho = cb.density(Tx, P, AIR_X)
        return surf.htc_and_T(Tx, P, AIR_X, m / (area * rho), D, 0.2, flow_area=area)

    d = 1e-4 * T
    fd = (h(T + d).h - h(T - d).h) / (2.0 * d)
    assert h(T).dh_dT == pytest.approx(fd, rel=1e-6)


@pytest.mark.parametrize("m", [0.0, 1e-13, -1e-13])
def test_1265_one_stream_keeps_its_partials_in_the_band(m: float) -> None:
    Y = _Y({"N2": 0.7, "O2": 0.2, "CO2": 0.1})
    r = _mix([cb.MassStream(m, 700.0, 2e5, Y)], Q=50.0)
    assert r.Tt_mix == pytest.approx(700.0, abs=1e-6)
    for var, h in (("Tt", 1e-3), ("Pt", 1.0)):
        up = _mix([cb.MassStream(m, 700.0 + h * (var == "Tt"), 2e5 + h * (var == "Pt"), Y)], Q=50.0)
        dn = _mix([cb.MassStream(m, 700.0 - h * (var == "Tt"), 2e5 - h * (var == "Pt"), Y)], Q=50.0)
        fd = [
            (a - b) / (2.0 * h)
            for a, b in zip([up.Tt_mix, up.Pt_mix], [dn.Tt_mix, dn.Pt_mix], strict=True)
        ]
        assert _analytic(r, 0, var)[:2] == pytest.approx(fd, rel=1e-6, abs=1e-9), var
    # Two streams in the band keep zero partials (#1238).
    two = _mix([cb.MassStream(m, 700.0, 2e5, Y), cb.MassStream(m, 300.0, 2e5, AIR)])
    assert all(v == 0.0 for v in _analytic(two, 0, "Tt"))


@pytest.mark.slow
def test_1265_heated_duct_into_a_boundary_at_zero_flow() -> None:
    # impingement:2, g1 (a wall-heated channel into a PressureBoundary) at
    # exactly zero flow: its row on g0.m_dot was J -1.412e-4 against FD
    # +1.420e-4. The crossflow rows' misses on g1.m_dot at this point are not
    # this pattern; only g1's own row is checked.
    case = Case("impingement", 2, 0, 0)
    x0, pt = point(case, "g1.m_dot = 0")
    rep = Report()
    evaluate(case, pt, x0, rep)
    assert rep.compared == 1
    assert not [f for f in rep.findings if f.element == "ChannelElement 'g1'"]
