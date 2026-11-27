"""fanno_rayleigh_channel_flow on a cooled duct whose needed drive folds (#1167).

The cold scan (no warm start) tested only the ends of each ln-2 bracket of
flux. On a cooled duct the exit-pressure excess can fall through zero and
rise above it again inside one bracket: the drive the duct needs peaks and
turns back down. The scan then saw no fall and reported the duct choked,
while a warm start found its first crossing unchoked, so the flow depended
on the element's _last_G. The C++ tests (test_fanno_rayleigh_fold_1167.cpp)
sweep the folds; these pin the issue's cases at the binding.
"""

from __future__ import annotations

import pytest

import combaero as cb

X = list(cb.species.dry_air())


def _flow(
    duct: tuple[float, float, float, float, float], drive: float, guess: float, et: bool
) -> cb.FannoChannelFlow:
    Pt0, Tt0, Tte, L, D = duct
    return cb.fanno_rayleigh_channel_flow(
        Pt0, Tt0, Tte, X, Pt0 - drive, et, L, D, 0.0, "haaland", 1.0, guess, [], False
    )


LAST = (3.0e6, 1400.0, 1000.0, 0.03, 0.01)
INTERIOR = (2.0e5, 1800.0, 400.0, 0.3, 0.01)
# Needed drive peaks at 354.198 Pa; the excess spans ~20x the drive.
NEAR_PEAK = (204095.106, 627.788, 327.196, 0.0475957, 0.0329094)


@pytest.mark.parametrize(
    ("duct", "drive", "et", "root"),
    [
        (LAST, 260294.0, False, 0.97533),  # the fold in the last bracket
        (LAST, 260957.0, False, 0.97963),
        (LAST, 261621.0, False, 0.98531),
        (INTERIOR, 100.0, False, 0.12232),  # the fold in an interior bracket
        (INTERIOR, 13.7, True, None),  # the laminar fold, total-pressure exit
        (NEAR_PEAK, 353.9, False, 0.47434),  # 0.3 Pa below a fold's peak
    ],
)
def test_cold_finds_the_first_crossing_and_warm_agrees(
    duct: tuple[float, float, float, float, float], drive: float, et: bool, root: float | None
) -> None:
    Pt0, Tt0, Tte, L, D = duct
    Gc = cb.fanno_rayleigh_choked_mass_flux(Pt0, Tt0, Tte, X, L, D, 0.0, "haaland", 1.0)
    c = _flow(duct, drive, -1.0, et)
    assert not c.choked
    assert Gc * (1.0 - 1e-4) > c.G
    if root is not None:
        assert c.G / Gc == pytest.approx(root, abs=1e-5)
    for s in (0.05, 0.3, 0.9, 0.99, 1.01, 1.2, 1.5):
        w = _flow(duct, drive, s * c.G, et)
        assert not w.choked, s
        assert pytest.approx(c.G, rel=1e-8) == w.G, s
