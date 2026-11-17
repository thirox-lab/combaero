"""fanno_rayleigh_channel_flow: a warm start above a later fall (#1191).

A warm start (G_guess, the channel element's _last_G) brackets the root
around the guess. On a cooled duct whose needed drive rises again after a
fold's valley, the excess falls, rises and falls again, and a guess within
6.4 % above the later fall returned that later crossing while the cold scan
returned the first fall: the flux depended on the solver's history. The C++
tests (test_fanno_rayleigh_warm_root_1191.cpp) cover more guesses and the
derivatives; this pins the issue's duct at the binding.
"""

from __future__ import annotations

import pytest

import combaero as cb

X = list(cb.species.dry_air())

# The issue's duct: total-pressure exit, first fall at 0.666323 Gc, a later
# one at 0.957556 Gc.
PT0, TT0, TTE, L, D, DRIVE = 1.28643e6, 1984.90, 1248.21, 0.15866, 5.8745e-3, 6797.70
# needed - drive changes sign within -5.3e-12..+6.3e-12 of the root in ln G
# (43 sign changes, measured); 1e-9 leaves 80x that, far below the 44% jump.
ROOT_BAND = 1e-9


def _flow(guess: float, tte: float = TTE) -> cb.FannoChannelFlow:
    return cb.fanno_rayleigh_channel_flow(
        PT0, TT0, tte, X, PT0 - DRIVE, True, L, D, 0.0, "haaland", 1.0, guess, [], False
    )


def _gc(tte: float = TTE) -> float:
    return cb.fanno_rayleigh_choked_mass_flux(PT0, TT0, tte, X, L, D, 0.0, "haaland", 1.0)


@pytest.mark.parametrize("s", [0.9, 0.95, 0.99, 0.999])
def test_warm_start_returns_the_first_fall(s: float) -> None:
    gc = _gc()
    cold = _flow(-1.0)
    assert cold.G / gc == pytest.approx(0.666323, abs=1e-6)
    warm = _flow(s * gc)
    assert not warm.choked
    # 0.95, 0.99 and 0.999 Gc returned 0.957556 Gc before.
    assert pytest.approx(cold.G, rel=ROOT_BAND) == warm.G


@pytest.mark.parametrize("tte", [TT0, 2400.0])
def test_heated_and_adiabatic_warm_starts_agree_with_cold(tte: float) -> None:
    gc = _gc(tte)
    cold = _flow(-1.0, tte)
    for s in (0.05, 0.5, 0.9, 0.99):
        assert pytest.approx(cold.G, rel=ROOT_BAND) == _flow(s * gc, tte).G
