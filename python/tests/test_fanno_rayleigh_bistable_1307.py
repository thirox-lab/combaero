"""A cooled duct can be bistable: a subsonic root and choked flow (#1307).

Cooling raises a duct's Pt (the Rayleigh gain), so its needed drive
Pt0 - Pt_exit(G) is not monotone in G: it rises from zero flow to a peak and
falls back. At a drive below that peak the duct has a subsonic root on the
rising branch, and choked flow holds too wherever the exit pressure at the
choked flux is at or above the target. As Tt_exit or L moves the peak below
the drive, the subsonic roots merge and vanish (a saddle-node), and the flux
jumps to the choke: no continuous root choice removes that.

The kernel keeps the start-up root (the first stable one, which a slow ramp
from rest reaches). fanno_rayleigh_choked_alternative says whether choked
flow also holds for the boundary values given; whether a network has a
choked solution is #1338. Case from the issue: Pt0 1e5 Pa, Tt0 400 K, L = D = 1 cm, exit
total, air.
"""

from __future__ import annotations

import math
from collections.abc import Callable

import numpy as np
import pytest

import combaero as cb

X = list(cb.species.dry_air())
Y = list(cb.mole_to_mass(cb.species.dry_air()))
PT0, TT0, D = 1e5, 400.0, 0.01


def _flow(Te: float, drive: float, L: float = 0.01) -> cb.FannoChannelFlow:
    return cb.fanno_rayleigh_channel_flow(
        PT0, TT0, Te, X, PT0 - drive, True, L, D, 0.0, "haaland", 1.0, -1.0, [], False
    )


def _alt(Te: float, drive: float, L: float = 0.01) -> cb.FannoChokedAlternative:
    return cb.fanno_rayleigh_choked_alternative(
        PT0, TT0, Te, X, PT0 - drive, True, L, D, 0.0, "haaland"
    )


@pytest.mark.parametrize("drive", [1.5, 3.0, 10.0])
def test_the_jump_to_choke_is_flagged_from_its_subsonic_side(drive: float) -> None:
    """Bisect Tt_exit to the fold, to 1e-11 K: below it the duct chokes, above
    it a subsonic root an order of magnitude smaller in G. That root is
    flagged, and its alternative is the choked flux the other side takes, so
    a solve sitting on either side names the other state."""
    lo, hi = 380.0, 400.0
    assert _flow(lo, drive).choked and not _flow(hi, drive).choked
    for _ in range(41):
        mid = 0.5 * (lo + hi)
        if _flow(mid, drive).choked:
            lo = mid
        else:
            hi = mid
    choked, sub = _flow(lo, drive), _flow(hi, drive)
    alt = _alt(hi, drive)
    # Measured: 7.9x at 1.5 Pa, 5.2x at 3 Pa, 2.5x at 10 Pa.
    assert choked.G / sub.G > 2.0
    assert alt.admissible and alt.exit_margin > 0.0
    assert pytest.approx(choked.G, rel=1e-12) == alt.G


@pytest.mark.parametrize(
    ("label", "values", "call"),
    [
        ("Tt_exit", np.linspace(388.0, 394.0, 25), lambda v: (v, 1.5, 0.01)),
        ("L", np.linspace(0.112, 0.128, 17), lambda v: (250.0, 1.5, float(v))),
    ],
)
def test_every_branch_flip_in_a_scan_is_flagged(
    label: str, values: np.ndarray, call: Callable[[float], tuple[float, float, float]]
) -> None:
    """Across the issue's Tt_exit and L folds (above the 1 Pa drive floor),
    each step between a choked and an unchoked result has the unchoked one
    flagged bistable. The scan must cross a fold, or it tests nothing."""
    flows = [_flow(*call(v)) for v in values]
    flips = [i for i in range(len(flows) - 1) if flows[i].choked != flows[i + 1].choked]
    assert flips, f"the {label} scan crosses no fold"
    for i in flips:
        j = i if not flows[i].choked else i + 1
        assert _alt(*call(values[j])).admissible, (label, values[j])


@pytest.mark.parametrize("drive", [1.0, 5.0])
@pytest.mark.parametrize("Te", [400.0, 420.0, 396.0])
def test_choked_flow_is_not_admissible_where_the_drive_is_monotone(Te: float, drive: float) -> None:
    """Adiabatic, heated, and cooled too little for the choke's exit Pt to
    reach the target: the exit pressure at choke is below it, no flag."""
    alt = _alt(Te, drive)
    assert not alt.admissible
    assert alt.exit_margin < 0.0


def test_a_static_exit_duct_can_be_bistable_too() -> None:
    """A static exit loses the exit head, which narrows the band but does not
    close it: strongly cooled, choked flow holds beside a subsonic root an
    order of magnitude smaller. Its alternative is the flux the same duct
    chokes at under a large drive."""
    args = (2.64e5, 899.6, 302.0, X)
    geom = (False, 0.738, 0.0243, 0.0, "haaland")
    sub = cb.fanno_rayleigh_channel_flow(*args, 2.64e5 - 19.6, *geom, 1.0, -1.0, [], False)
    alt = cb.fanno_rayleigh_choked_alternative(*args, 2.64e5 - 19.6, *geom)
    hard = cb.fanno_rayleigh_channel_flow(*args, 1e5, *geom, 1.0, -1.0, [], False)
    assert not sub.choked and hard.choked
    assert alt.admissible and alt.G > 10.0 * sub.G
    assert pytest.approx(hard.G, rel=1e-10) == alt.G


@pytest.mark.parametrize("drive", [0.0, -5.0])
def test_no_forward_drive_has_no_choked_alternative(drive: float) -> None:
    """With no forward drive the channel flow is no flow, and there is no
    forward choked state to compare: not admissible, nothing evaluated. Te
    391 K sits inside the band, where a forward drive would admit it."""
    alt = _alt(391.0, drive)
    assert not alt.admissible
    assert alt.G == 0.0 and math.isnan(alt.exit_margin)
    assert _flow(391.0, drive).G == 0.0
    assert _alt(391.0, 1.5).admissible
