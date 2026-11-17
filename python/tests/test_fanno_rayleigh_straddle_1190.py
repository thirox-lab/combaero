"""fanno_rayleigh_channel_flow: a straddling cold-scan bracket with three crossings (#1190).

The cold scan took the first ln-2 bracket of flux whose excess (exit pressure
minus the target) went from positive to non-positive and Brent-solved it at
once. On a cooled duct whose needed drive folds in its laminar range, that
bracket can hold a fall, a rise and a later fall, and Brent returned the
third crossing: the flux was not the first fall from zero flow, and a warm
start from below it gave a different one. The C++ tests
(test_fanno_rayleigh_straddle_1190.cpp) sweep more ducts; this pins the
issue's duct at the binding.
"""

from __future__ import annotations

import pytest

import combaero as cb

X = list(cb.species.dry_air())

# The issue's duct: dry air, static-pressure exit. Its needed drive peaks at
# 2.15199 Pa near 0.0161 Gc and dips to 2.12681 Pa near 0.0197 Gc.
PT0, TT0, TTE, L, D = 199111.5, 868.17, 348.98, 0.09662, 9.1945e-3
# The roots are defined only to the march's rounding: needed - drive changes
# sign within up to 3.3e-7 of them in ln G on this duct (measured). 30x that,
# far below the 41% jump to the third crossing.
ROOT_BAND = 1e-5


def _gc() -> float:
    return cb.fanno_rayleigh_choked_mass_flux(PT0, TT0, TTE, X, L, D, 0.0, "haaland", 1.0)


def _flow(drive: float, guess: float = -1.0) -> cb.FannoChannelFlow:
    return cb.fanno_rayleigh_channel_flow(
        PT0, TT0, TTE, X, PT0 - drive, False, L, D, 0.0, "haaland", 1.0, guess, [], False
    )


def _needed(G: float) -> float:
    return PT0 - cb.fanno_rayleigh_march(G, PT0, TT0, TTE, X, L, D, 0.0, "haaland").P_exit


def _first_fall(drive: float, gc: float) -> float:
    """First flux on a 0.0001 Gc grid whose needed drive reaches drive, bisected."""
    a = 0.005 * gc
    assert _needed(a) < drive
    b = a
    for i in range(1, 301):
        b = gc * (0.005 + 0.0001 * i)
        if _needed(b) >= drive:
            break
        a = b
    for _ in range(60):
        m = 0.5 * (a + b)
        if _needed(m) >= drive:
            b = m
        else:
            a = m
    return b


def test_issue_drive_returns_the_first_fall() -> None:
    gc = _gc()
    c = _flow(2.145)
    assert not c.choked
    # The issue's table: the first fall at 0.0150235 Gc, the third crossing
    # (what the cold scan returned) at 0.021177 Gc.
    assert c.G / gc == pytest.approx(0.0150235, abs=1e-6)
    assert pytest.approx(_first_fall(2.145, gc), rel=ROOT_BAND) == c.G
    # A warm start from below the first fall agrees.
    assert pytest.approx(c.G, rel=ROOT_BAND) == _flow(2.145, 0.0145 * gc).G


@pytest.mark.parametrize("drive", [2.1280, 2.1340, 2.1380, 2.1420, 2.1460, 2.1490, 2.1510, 2.1519])
def test_fold_band_returns_the_first_fall(drive: float) -> None:
    # Between the valley and the peak: the old scan returned the third
    # crossing from 2.1378 Pa up to the peak.
    gc = _gc()
    c = _flow(drive)
    assert not c.choked
    assert 0.0162 * gc > c.G  # below the peak's flux
    assert pytest.approx(_first_fall(drive, gc), rel=ROOT_BAND) == c.G
