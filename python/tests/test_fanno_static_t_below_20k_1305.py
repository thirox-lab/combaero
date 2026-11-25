"""Fanno-Rayleigh with a stagnation temperature below 20 K (#1305).

`static_T` (src/fanno_mach.cpp) bracketed the static temperature on
[max(20 K, 0.05 Tt), Tt]. For Tt at or below 20 K the bracket was inverted,
and the safeguarded Newton could return the 20 K bound: a static T above Tt,
an exit static pressure above the exit total, and an E(G) that flipped
between two branches with the iteration's start. On a duct cooled from 450 K
to 9.7 K the march's exit total jumped between about 9.6e5 and 1.1e5 Pa with
G, and a warm-started solve returned a bracket edge as a root (residual
2.4e4 Pa, no flag). The floor is now 0.05 Tt where 20 K is not below Tt.
"""

import numpy as np
import pytest

import combaero as cb

_N = [cb.species_name(k) for k in range(cb.num_species())]
_X = [0.0] * len(_N)
for _name, _x in {"N2": 0.78084, "O2": 0.20946, "AR": 0.00934, "CO2": 0.00036}.items():
    _X[_N.index(_name)] = _x

# #1305's case: Pt0 1 MPa, Tt 450 K -> 9.7155 K, exit total, L 0.5 m, D 4 mm.
PT0, TT0, TTE, PT_TARGET, L, D = 1e6, 450.0, 9.7155, 9.4463e5, 0.5, 0.004


def _march(G: float) -> cb._core.FannoRayleighMarch:
    return cb._core.fanno_rayleigh_march(G, PT0, TT0, TTE, _X, L, D, 0.0, "haaland")


def test_the_exit_state_is_physical_below_20k() -> None:
    """Static T at or below the stagnation T, static P at or below the total,
    and the exit total continuous in G, across the range where it flipped."""
    Gs = np.geomspace(0.004, 0.0095, 40)
    pts = []
    for G in Gs:
        s = _march(G)
        assert s.T_exit <= TTE * (1.0 + 1e-12), (G, s.T_exit)
        assert s.P_exit <= s.Pt_exit * (1.0 + 1e-12), (G, s.P_exit, s.Pt_exit)
        pts.append(s.Pt_exit)
    # On main it jumped by about 8.6e5 Pa between neighbours; it varies by
    # well under 1 Pa over this whole range.
    assert np.max(np.abs(np.diff(pts))) < 1.0


@pytest.mark.parametrize("G_guess", [-1.0, 1.761])
def test_cold_and_warm_starts_agree_and_choke(G_guess: float) -> None:
    """Below choke the drive never reaches the target (it peaks near 8.7e3
    Pa against 5.5e4), so the duct chokes. Both starts say so; on main the
    warm start returned G = 0.0044, not a root and not flagged."""
    r = cb._core.fanno_rayleigh_channel_flow(
        PT0, TT0, TTE, _X, PT_TARGET, True, L, D, 0.0, "haaland", 1.0, G_guess, [], False
    )
    Gc = cb._core.fanno_rayleigh_choked_mass_flux(PT0, TT0, TTE, _X, L, D, 0.0, "haaland")
    assert r.choked
    assert pytest.approx(Gc, rel=1e-9) == r.G
