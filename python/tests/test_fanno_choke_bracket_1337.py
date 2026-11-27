"""Regression tests for #1337: long laminar Fanno and Fanno-Rayleigh choke brackets."""

from __future__ import annotations

import math

import pytest

import combaero as cb

X = list(cb.species.dry_air())

LONG_LAMINAR_CASES = [
    ("issue_1337_repro", 1.0e5, 1000.0, 1000.0, 0.10, 1.0e-4),
    ("pr_1336_review", 35416.0, 787.0, 787.0, 4.28, 1.01e-3),
    ("table_row1_T1000_1bar_D0p1mm", 1.0e5, 1000.0, 1000.0, 0.072, 1.0e-4),
    ("table_row2_T1000_1bar_D1mm", 1.0e5, 1000.0, 1000.0, 7.2, 1.0e-3),
    ("table_row3_T1000_1bar_D10mm", 1.0e5, 1000.0, 1000.0, 720.0, 1.0e-2),
    ("table_row4_T1000_10bar_D0p1mm", 1.0e6, 1000.0, 1000.0, 0.72, 1.0e-4),
    ("table_row5_T300_1bar_D0p1mm", 1.0e5, 300.0, 300.0, 0.30, 1.0e-4),
    ("table_row6_T300_10bar_D0p1mm", 1.0e6, 300.0, 300.0, 3.0, 1.0e-4),
    ("extreme_LD_1e5", 1.0e5, 300.0, 300.0, 100.0, 1.0e-3),
]


@pytest.mark.parametrize(("name", "pt0", "tt0", "tte", "L", "D"), LONG_LAMINAR_CASES)
def test_long_laminar_choked_flux_and_flow_match_adiabatic(
    name: str, pt0: float, tt0: float, tte: float, L: float, D: float
) -> None:
    gc_ad = cb.fanno_choked_mass_flux(pt0, tt0, X, L, D, 0.0, "haaland", 1.0)
    gc_ray = cb.fanno_rayleigh_choked_mass_flux(pt0, tt0, tte, X, L, D, 0.0, "haaland", 1.0)
    assert math.isfinite(gc_ad) and gc_ad > 0.0, name
    assert math.isfinite(gc_ray) and gc_ray > 0.0, name
    assert gc_ray == pytest.approx(gc_ad, rel=1e-7), name

    f_ad = cb.fanno_channel_flow(pt0, tt0, X, 1.0, False, L, D, 0.0, "haaland", 1.0)
    f_ray = cb.fanno_rayleigh_channel_flow(pt0, tt0, tte, X, 1.0, False, L, D, 0.0, "haaland", 1.0)
    assert f_ad.choked and f_ray.choked, name
    assert pytest.approx(f_ad.G, rel=1e-7) == f_ray.G, name

    f_sub = cb.fanno_rayleigh_channel_flow(
        pt0, tt0, tte, X, 0.5 * pt0, False, L, D, 0.0, "haaland", 1.0, -1.0, [], True
    )
    assert not f_sub.choked, name
    assert 0.0 < f_sub.G < gc_ray, name
    assert f_sub.dG_dP_target < 0.0, name


@pytest.mark.parametrize(
    ("name", "tt0", "tte", "L", "D"),
    [
        ("heated_LD1000", 300.0, 600.0, 1.0, 1.0e-3),
        ("cooled_LD1000", 600.0, 300.0, 1.0, 1.0e-3),
        ("heated_LD1e5", 300.0, 600.0, 100.0, 1.0e-3),
        ("cooled_LD1e5", 600.0, 300.0, 100.0, 1.0e-3),
    ],
)
def test_nonadiabatic_long_laminar_choked_and_unchoked_jacobians(
    name: str, tt0: float, tte: float, L: float, D: float
) -> None:
    pt0 = 1.0e5
    fc = cb.fanno_rayleigh_channel_flow(
        pt0, tt0, tte, X, 1.0, False, L, D, 0.0, "haaland", 1.0, -1.0, [], True
    )
    assert fc.choked, name
    assert fc.dG_dP_target == 0.0, name

    h = 1e-4 * pt0
    fd_pt = (
        cb.fanno_rayleigh_choked_mass_flux(pt0 + h, tt0, tte, X, L, D, 0.0, "haaland", 1.0)
        - cb.fanno_rayleigh_choked_mass_flux(pt0 - h, tt0, tte, X, L, D, 0.0, "haaland", 1.0)
    ) / (2.0 * h)
    assert fc.dG_dPt0 == pytest.approx(fd_pt, rel=1e-5), name

    fu = cb.fanno_rayleigh_channel_flow(
        pt0, tt0, tte, X, 0.5 * pt0, False, L, D, 0.0, "haaland", 1.0, -1.0, [], True
    )
    assert not fu.choked, name
    hp = 1e-4 * 0.5 * pt0
    fd_p = (
        cb.fanno_rayleigh_channel_flow(
            pt0, tt0, tte, X, 0.5 * pt0 + hp, False, L, D, 0.0, "haaland", 1.0
        ).G
        - cb.fanno_rayleigh_channel_flow(
            pt0, tt0, tte, X, 0.5 * pt0 - hp, False, L, D, 0.0, "haaland", 1.0
        ).G
    ) / (2.0 * hp)
    assert fu.dG_dP_target == pytest.approx(fd_p, rel=1e-5), name
