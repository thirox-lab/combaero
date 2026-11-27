"""Regression tests for #1276: exact discrete tangents in fanno_rayleigh_channel_flow
and second temperature derivatives d2cp_R_dT2 / d2cp_dT2."""

from __future__ import annotations

from collections.abc import Callable

import pytest

import combaero as cb


def _rfd(fn: Callable[[float], float], x0: float, h: float) -> float:
    d1 = (fn(x0 + h) - fn(x0 - h)) / (2.0 * h)
    d2 = (fn(x0 + 0.5 * h) - fn(x0 - 0.5 * h)) / h
    return (4.0 * d2 - d1) / 3.0


def test_d2cp_python_bindings_match_fd_across_regions() -> None:
    air = cb.species.dry_air()
    for T in (120.0, 350.0, 975.0, 1000.0, 1025.0, 1800.0, 6200.0):
        h = 1e-4 * T
        ana = cb.d2cp_dT2(T, air)
        fd = _rfd(lambda t: cb.dcp_dT(t, air), T, h)
        assert ana == pytest.approx(fd, rel=1e-8, abs=1e-12)

    n2_idx = cb.species_index_from_name("N2")
    for T in (80.0, 400.0, 990.0, 1010.0, 2500.0, 6100.0, 25000.0):
        h = 1e-4 * T
        ana = cb.d2cp_R_dT2(n2_idx, T)
        fd = _rfd(lambda t: cb.dcp_R_dT(n2_idx, t), T, h)
        assert ana == pytest.approx(fd, rel=1e-8, abs=1e-12)
    assert cb.d2cp_R_dT2(n2_idx, 30000.0) == 0.0


@pytest.mark.parametrize(
    ("label", "Pt0", "Tt0", "Tt_ex", "P_t", "exit_total", "L", "D", "rough", "model", "fmult"),
    [
        ("heated_haaland", 5e5, 800.0, 1200.0, 3.0e5, False, 0.30, 0.02, 1e-5, "haaland", 1.0),
        ("cooled_colebrook", 5e5, 900.0, 700.0, 4.6e5, True, 2.0, 0.02, 5e-5, "colebrook", 1.0),
        ("fixed_static", 4e5, 500.0, 900.0, 2.5e5, False, 0.40, 0.02, 0.0, "fixed", 0.025),
    ],
)
def test_unchoked_discrete_tangents_match_march_fd(
    label: str,
    Pt0: float,
    Tt0: float,
    Tt_ex: float,
    P_t: float,
    exit_total: bool,
    L: float,
    D: float,
    rough: float,
    model: str,
    fmult: float,
) -> None:
    air = cb.species.dry_air()
    r = cb.fanno_rayleigh_channel_flow(
        Pt0, Tt0, Tt_ex, air, P_t, exit_total, L, D, rough, model, fmult, -1.0, [], True
    )
    assert not r.choked, label

    def e_of(Gx: float, Ptx: float, T0x: float, Tex: float) -> float:
        m = cb.fanno_rayleigh_march(Gx, Ptx, T0x, Tex, air, L, D, rough, model, fmult)
        return m.Pt_exit if exit_total else m.P_exit

    e_g = _rfd(lambda v: e_of(v, Pt0, Tt0, Tt_ex), r.G, 1e-4 * r.G)
    e_pt = _rfd(lambda v: e_of(r.G, v, Tt0, Tt_ex), Pt0, 1e-4 * Pt0)
    e_t0 = _rfd(lambda v: e_of(r.G, Pt0, v, Tt_ex), Tt0, 1e-4 * Tt0)
    e_te = _rfd(lambda v: e_of(r.G, Pt0, Tt0, v), Tt_ex, 1e-4 * Tt_ex)

    assert r.dG_dP_target == pytest.approx(1.0 / e_g, rel=1e-8)
    assert r.dG_dPt0 == pytest.approx(-e_pt / e_g, rel=1e-8)
    assert r.dG_dTt0 == pytest.approx(-e_t0 / e_g, rel=1e-8)
    assert r.dG_dTt_exit == pytest.approx(-e_te / e_g, rel=1e-8)


@pytest.mark.parametrize(
    ("label", "Pt0", "Tt0", "Tt_ex", "L", "D", "model", "fmult", "tol"),
    [
        ("kInlet", 5e5, 1000.0, 300.0, 0.02, 0.05, "haaland", 1.0, 1e-8),
        ("kExit_heated", 5e5, 600.0, 1100.0, 0.50, 0.02, "serghides", 1.0, 1e-7),
        ("kExit_cooled", 5e5, 800.0, 650.0, 2.00, 0.02, "petukhov", 1.0, 1e-7),
        ("kThroat", 5e5, 600.0, 320.0, 1.00, 0.02, "haaland", 1.0, 1e-4),
    ],
)
def test_choked_tangents_match_choked_flux_fd(
    label: str,
    Pt0: float,
    Tt0: float,
    Tt_ex: float,
    L: float,
    D: float,
    model: str,
    fmult: float,
    tol: float,
) -> None:
    air = cb.species.dry_air()
    r = cb.fanno_rayleigh_channel_flow(
        Pt0, Tt0, Tt_ex, air, 1e4, False, L, D, 1e-5, model, fmult, -1.0, [], True
    )
    assert r.choked, label
    assert r.dG_dP_target == 0.0

    def gc(Ptx: float, T0x: float, Tex: float) -> float:
        return cb.fanno_rayleigh_choked_mass_flux(Ptx, T0x, Tex, air, L, D, 1e-5, model, fmult)

    fd_pt = _rfd(lambda v: gc(v, Tt0, Tt_ex), Pt0, 1e-4 * Pt0)
    fd_t0 = _rfd(lambda v: gc(Pt0, v, Tt_ex), Tt0, 1e-4 * Tt0)
    fd_te = _rfd(lambda v: gc(Pt0, Tt0, v), Tt_ex, 1e-4 * Tt_ex)

    assert r.dG_dPt0 == pytest.approx(fd_pt, rel=tol, abs=1e-12)
    assert r.dG_dTt0 == pytest.approx(fd_t0, rel=tol, abs=1e-12)
    assert r.dG_dTt_exit == pytest.approx(fd_te, rel=tol, abs=1e-12)
