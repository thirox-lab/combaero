"""NASA TN D-7156 (Marek and Juhasz 1973) as a network case (#764).

A plate in a combustor's exhaust with a slot film on its hot face, a
parallel-flow convection passage under it, or both; 62 runs by 10 hot-face
thermocouples (Table II, transcribed as printed and signed off 2026-11-03,
validation/cooling/extractions/marek_juhasz_1973_d7156.md). The network and
every labelled assumption are in validation/cooling/marek_juhasz_network_runner.py.

CROSS-SOURCE ACCURACY, every number (same lab as D-7958's kernel, different
study; docs/VALIDATION_POLICY.md). Nothing is tuned to these data: the film is
D-7958's kernel at the element's default Tu 0.14, the hot-side h without a
film is combaero's duct h, the film-side h Seban's. The pins record what the
network does; the +-11 K band is the report's own measurement accuracy (p. 13).
"""

from __future__ import annotations

import csv
import math

import numpy as np
import pytest

import combaero as cb
from validation.cooling import marek_juhasz_network_runner as mj

pytestmark = pytest.mark.slow

TOL = 1.0  # K on a pinned bias or RMS (platform-to-platform solver noise)


# --- the transcription ----------------------------------------------------


def _rows(name: str) -> list[dict[str, str]]:
    with (mj.DATA / name).open() as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def test_the_run_classes_partition_table_ii() -> None:
    runs = mj.runs()
    assert sorted(runs) == list(range(1, 63))
    ids = [i for ids in mj.STAGES.values() for i in ids]
    assert sorted(ids) == list(range(1, 63))
    for r in runs.values():
        assert (r.m_film > 0.0) == (r.stage in ("film", "combined")), r.run
        assert (r.m_conv > 0.0) == (r.stage in ("convection", "combined")), r.run


def test_the_si_table_is_the_us_table() -> None:
    """Every SI cell against its US-customary counterpart, to the print's
    rounding (the transcription's cross-unit check)."""
    a = {r["run"]: r for r in _rows("table2a_experimental_data_SI.csv")}
    b = {r["run"]: r for r in _rows("table2b_experimental_data_US.csv")}
    lb = 0.45359237
    for run, ra in a.items():
        rb = b[run]
        assert float(ra["m_hot_kg_s"]) == pytest.approx(float(rb["m_hot_lb_s"]) * lb, abs=0.011)
        assert float(ra["m_film_kg_s"]) == pytest.approx(float(rb["m_film_lb_s"]) * lb, abs=0.0011)
        for ka, kb in zip(list(ra)[10:], list(rb)[10:], strict=True):
            T_F = float(rb[kb])
            assert float(ra[ka]) == pytest.approx((T_F - 32.0) / 1.8 + 273.15, abs=1.6), (run, ka)


def test_the_hot_gas_area_reproduces_table_i() -> None:
    """The reconstructed 0.0370 m^2 (38.1 x 10.6 cm less the 3.3 x 10.2 cm
    test section) returns Table I's printed velocities within 1 %."""
    for mh, Th, u in ((1.59, 589.0, 72.0), (1.59, 1033.0, 125.0), (2.27, 1033.0, 180.0)):
        X = cb.mass_to_mole(list(mj.vitiated(Th)[0]))
        assert mh / (cb.density(Th, mj.P_ATM, X) * mj.HOT_AREA) == pytest.approx(u, rel=0.01)


def test_the_hot_gas_is_vitiated_air() -> None:
    """Kerosene stand-in n-heptane, complete combustion: 290 -> 589 K takes
    f/a 0.0069, on to 1033 K 0.0185 overall; O2 falls, CO2 and H2O rise."""
    _, fa1 = mj.vitiated(589.0)
    Y, fa = mj.vitiated(1033.0)
    assert fa1 == pytest.approx(0.0069, abs=2e-4)
    assert fa == pytest.approx(0.0185, abs=3e-4)
    X = dict(zip(cb.species.names, cb.mass_to_mole(list(Y)), strict=True))
    assert X["O2"] == pytest.approx(0.1475, abs=2e-3)
    assert X["H2O"] > X["CO2"] > 0.03
    assert X[mj.SURROGATE] == 0.0


# --- the network: staged scores -------------------------------------------


def _pin(sc: mj.Score, bias: float, rms: float) -> None:
    assert sc.bias == pytest.approx(bias, abs=TOL)
    assert sc.rms == pytest.approx(rms, abs=TOL)


class TestConvectionOnly:
    """15 runs. The duct h is low for a plate whose boundary layer starts at
    its leading edge 6 cm upstream: the wall reads 49 K cold, TC 4-10 58 K."""

    def test_scores(self) -> None:
        res = mj.solve_all(mj.CONVECTION_ONLY)
        s = mj.stage_scores(res)
        _pin(s["all"], -49.2, 59.1)
        _pin(s["TC1-3"], -28.7, 51.7)
        _pin(s["TC4-10"], -58.0, 62.0)
        exits = np.array([r.T_conv_exit - r.run.T_conv_exit for r in res])
        _pin(mj.score(exits), -12.4, 15.2)


class TestFilmOnly:
    """28 runs, adiabatic backside: the wall IS the film's adiabatic wall
    temperature, so this stage scores the film effectiveness alone."""

    def test_scores(self) -> None:
        res = mj.solve_all(mj.FILM_ONLY)
        s = mj.stage_scores(res)
        _pin(s["all"], 24.8, 34.8)
        _pin(s["TC1-3"], 13.6, 22.7)
        _pin(s["TC4-10"], 29.6, 38.9)
        # An adiabatic side B (#970): no heat, the wall sits at T_aw exactly.
        for r in res:
            assert all(q == 0.0 for q in r.extras["Q"]), r.run.run
            assert list(r.T_wall) == list(r.extras["T_aw_a"]), r.run.run


class TestCombined:
    """17 runs, film and convection together, film_htc='seban1960'."""

    def test_scores(self) -> None:
        res = mj.solve_all(mj.COMBINED)
        s = mj.stage_scores(res)
        _pin(s["all"], 2.4, 43.9)
        _pin(s["TC1-3"], 23.9, 32.5)
        _pin(s["TC4-10"], -6.9, 48.0)
        hi, lo = mj._split(res)
        assert all(r.extras["film_blowing_ratio"][0] >= 1.0 for r in hi)
        assert all(r.extras["film_blowing_ratio"][0] < 1.0 for r in lo)
        _pin(mj.stage_scores(hi)["all"], 36.6, 38.8)
        _pin(mj.stage_scores(lo)["all"], -21.5, 47.1)
        exits = np.array([r.T_conv_exit - r.run.T_conv_exit for r in res])
        _pin(mj.score(exits), -5.8, 8.1)


class TestFilmSideHtcComparison:
    """METHOD COMPARISON, combined stage. Falsification: the film-side h must
    act where Seban's box applies (m >= 1) and nowhere else. Its higher h near
    the slot drives the wall towards the film's T_aw, hotter here, so the
    m >= 1 bias must RISE from 'none' to 'seban1960'; at m < 1 the side hands
    back to the duct's h and the three must agree exactly."""

    def test_none(self) -> None:
        res = mj.solve_all(mj.COMBINED, mj.Options(film_htc="none"))
        hi, lo = mj._split(res)
        _pin(mj.stage_scores(res)["all"], -7.8, 37.3)
        _pin(mj.stage_scores(hi)["all"], 11.8, 14.4)
        _pin(mj.stage_scores(lo)["all"], -21.5, 47.1)
        seban = {r.run.run: r for r in mj.solve_all(mj.COMBINED_M_GE_1[:2])}
        for r in hi[:2]:
            assert np.mean(seban[r.run.run].errors) > np.mean(r.errors) + 10.0
        assert all(w == 0.0 for r in lo for w in r.extras["film_h_weight_M"] if w is not None)

    def test_ballal_lefebvre(self) -> None:
        res = mj.solve_all(mj.COMBINED, mj.Options(film_htc="ballal_lefebvre1973"))
        hi, lo = mj._split(res)
        _pin(mj.stage_scores(res)["all"], 0.7, 44.2)
        _pin(mj.stage_scores(hi)["all"], 32.5, 39.6)
        _pin(mj.stage_scores(lo)["all"], -21.5, 47.1)


def test_the_wall_k_is_immaterial() -> None:
    """0.159 cm of stainless: k from 13.8 to 22.5 W/m K moves one combined
    run's wall by under 1 K, so the wall's k(T) barely matters."""
    run = mj.runs()[30]
    lo = mj.solve(run, mj.Options(wall_k=13.8)).T_wall
    hi = mj.solve(run, mj.Options(wall_k=22.5)).T_wall
    assert max(abs(a - b) for a, b in zip(lo, hi, strict=True)) < 1.0
    assert not math.isclose(lo[-1], hi[-1], abs_tol=1e-6)


def test_run_19_starts_at_its_root_and_converges_cold() -> None:
    """Run 19's coolant passage runs at up to Mach 0.7 at 297 K. Seeded at its
    thermal group's mean 446 K, its pressures started 19% high (172.8 against
    a 145.2 kPa root) and the default cold solve took the warm-start retry;
    the resistance seed now runs at each stream's own Tt (#990)."""
    import warnings

    from combaero.network import NetworkSolver

    net = mj.build(mj.runs()[19])
    s = NetworkSolver(net)
    net.resolve_all_topology()
    seed = s._build_x0()
    x0 = dict(zip(s.unknown_names, seed, strict=True))
    assert x0["k1.P"] == pytest.approx(145.2e3, rel=0.02)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        # stall_handover=False: the cold attempt's stall detector measures
        # progress per 5 s of wall clock, so on a loaded runner (the slow tier
        # on 4 xdist workers) it handed this ~6 s solve to the retry although
        # the solve converges cold; this test is about the seed (#990).
        r = NetworkSolver(mj.build(mj.runs()[19])).solve(stall_handover=False)
    assert r["__success__"], r["__message__"]
    assert not [w for w in caught if "auto-retrying" in str(w.message)]
