"""fanno_rayleigh_channel_flow's cold path below the choke (#1161).

The cold scan (no warm start) gave up on its last bracket, whose upper end
the march reads as choked: every root above ~0.93 Gc was missed, and the
kernel returned G = Gc unchoked, with a nonzero dG/dP_target on a flat
plateau. A warm start found the real root, so the flux depended on the
element's _last_G (found under #1158, stress hotcold 10/11/12/16/18/26).
"""

from __future__ import annotations

import pytest

import combaero as cb

X = list(cb.species.dry_air())
# The cold duct of stress hotcold seed 10.
PT, TT, TTE, L, D = 1.0e6, 300.0, 697.9029088757375, 0.2, 0.008
DIR = [1.0, -1.0] + [0.0] * (len(X) - 2)


def _flow(
    drive: float, guess: float = -1.0, et: bool = True, der: bool = False, **kw: float
) -> cb.FannoChannelFlow:
    a = {"Pt0": PT, "Tt0": TT, "Tte": TTE, "X": X}
    a.update(kw)
    return cb.fanno_rayleigh_channel_flow(
        a["Pt0"],
        a["Tt0"],
        a["Tte"],
        a["X"],
        PT - drive,
        et,
        L,
        D,
        0.0,
        "haaland",
        1.0,
        guess,
        [DIR] if der else [],
        der,
    )


def _choking_drive(et: bool) -> float:
    a, b = 1.0, 0.9 * PT
    for _ in range(45):
        m = 0.5 * (a + b)
        a, b = (a, m) if _flow(m, et=et).choked else (m, b)
    return b


@pytest.mark.parametrize("et", [True, False])
def test_cold_and_warm_agree_and_the_flux_rises_to_the_choke(et: bool) -> None:
    Gc = cb.fanno_rayleigh_choked_mass_flux(PT, TT, TTE, X, L, D, 0.0, "haaland", 1.0)
    dc = _choking_drive(et)
    G_prev, G_warm = 0.0, -1.0
    for k in range(1, 41):
        drive = dc * (0.80 + 0.2 * k / 40) * (1.0 - 1e-6)
        c = _flow(drive, et=et)
        w = _flow(drive, G_warm, et=et)
        assert not c.choked and not w.choked, drive
        assert pytest.approx(w.G, rel=1e-11) == c.G, drive
        assert G_prev < c.G, drive
        assert Gc * (1.0 - 1e-7) > c.G, drive  # off the plateau
        for g in (0.3 * c.G, 0.999 * Gc, 1.2 * c.G):
            assert pytest.approx(c.G, rel=1e-11) == _flow(drive, g, et=et).G, (drive, g)
        G_prev, G_warm = c.G, w.G
    assert Gc == _flow(dc * (1.0 + 1e-6), et=et).G


@pytest.mark.parametrize("et", [True, False])
@pytest.mark.parametrize("frac", [0.9, 0.98, 0.999])
def test_derivatives_match_differences_on_the_former_plateau(et: bool, frac: float) -> None:
    """Central differences of the kernel's own G, steps of 1e-4 of the drive
    (>= 20 Pa): clear of the C1 blend below the 1 Pa drive floor. The
    kernel's own steps carry ~4e-5 of truncation where the exit is near sonic
    (kRayleighZetaStepFrac); before #1161 they straddled the sonic point within
    ~1e-4 G of the choke, 4-96% off."""
    drive = frac * _choking_drive(et)
    f = _flow(drive, et=et, der=True)
    assert not f.choked
    h = 1e-4 * drive
    fd_P = (_flow(drive - h, et=et).G - _flow(drive + h, et=et).G) / (2 * h)
    assert f.dG_dP_target == pytest.approx(fd_P, rel=1e-4)
    fd_Pt = (_flow(drive, et=et, Pt0=PT + h).G - _flow(drive, et=et, Pt0=PT - h).G) / (2 * h)
    assert f.dG_dPt0 == pytest.approx(fd_Pt, rel=1e-4)
    hT = 1e-5 * TT
    fd_T0 = (_flow(drive, et=et, Tt0=TT + hT).G - _flow(drive, et=et, Tt0=TT - hT).G) / (2 * hT)
    assert f.dG_dTt0 * TT / f.G == pytest.approx(fd_T0 * TT / f.G, abs=1e-4)
    hT = 1e-5 * TTE
    fd_Te = (_flow(drive, et=et, Tte=TTE + hT).G - _flow(drive, et=et, Tte=TTE - hT).G) / (2 * hT)
    assert f.dG_dTt_exit * TTE / f.G == pytest.approx(fd_Te * TTE / f.G, abs=1e-4)
    Y = list(cb.mole_to_mass(X))
    hy = 1e-6
    Xp = list(cb.mass_to_mole([y + hy * v for y, v in zip(Y, DIR, strict=True)]))
    Xm = list(cb.mass_to_mole([y - hy * v for y, v in zip(Y, DIR, strict=True)]))
    fd_Y = (_flow(drive, et=et, X=Xp).G - _flow(drive, et=et, X=Xm).G) / (2 * hy)
    assert f.dG_ddir[0] / f.G == pytest.approx(fd_Y / f.G, abs=1e-4)
