"""EffusionPlateElement against Behrendt et al. (2008), through the network.

ACCURACY, cross-source: no closure of the element comes from this paper or
lab. Both gas-side baselines are scored at gas_augmentation 1.0 and never
pooled; nothing is fitted. These pins hold the table in
validation/cooling/extractions/behrendt_2008_effusion_eta_tot.md.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.cooling import effusion_network_runner as nr  # noqa: E402
from validation.cooling.schema import load_dataset  # noqa: E402

# Plateau (x/d 18-50) predictions at gas_augmentation 1.0: (chamber, flat plate).
PINNED = {
    "fig9_eta_tot_E20_dp2": (0.926, 0.860),
    "fig9_eta_tot_E20_dp3": (0.934, 0.874),
    "fig9_eta_tot_E20_dp4": (0.939, 0.883),
    "fig6_eta_tot_E30_Tc470": (0.924, 0.856),
    "fig6_eta_tot_E30_Tc540": (0.925, 0.857),
    "fig6_eta_tot_E30_Tc580": (0.924, 0.856),
    "fig10_eta_tot_E30_5bar_dp2": (0.915, 0.840),
    "fig10_eta_tot_E30_10bar_dp3": (0.914, 0.839),
    "fig10_eta_tot_E30_10bar_dp2": (0.906, 0.824),
}


@pytest.fixture(scope="module")
def dataset():
    return load_dataset()


@pytest.fixture(scope="module")
def series(dataset):
    return {s.path.stem: s for s in nr.scored(dataset)}


def test_owns_every_behrendt_series_and_scores_only_eta(dataset) -> None:
    owned = [s for s in dataset if nr.owns(s)]
    assert len(owned) == 12
    assert {s.path.stem for s in nr.scored(dataset)} == set(PINNED)
    for s in owned:
        if s.y_axis != "eta_overall":
            assert all(r.predicted is None for r in nr.run_series(s)), s.label


def test_the_plateau_predictions_are_pinned(series) -> None:
    for stem, (chamber, flat) in PINNED.items():
        s = series[stem]
        assert nr.predict(s, "chamber")["eta_overall"] == pytest.approx(chamber, abs=1e-3), stem
        assert nr.predict(s, "flat_plate")["eta_overall"] == pytest.approx(flat, abs=1e-3), stem


def test_no_tuner_is_applied(series) -> None:
    for s in series.values():
        for b in nr.BASELINES:
            for rad in nr.RADIATION:
                d = nr.predict(s, b, radiation=rad)
                assert d["gas_augmentation"] == 1.0
                assert d["internal_Nu_multiplier"] == 1.0
                assert d["eta_film"] == 0.0


def test_both_baselines_over_predict_and_the_flat_plate_less(series) -> None:
    """The restarted boundary layer roughly halves the miss and does not close
    it: every curve stays above the authors' 3% band on both baselines."""
    for stem, s in series.items():
        m = nr.plateau(s)
        ch = nr.predict(s, "chamber")["eta_overall"] / m - 1.0
        fp = nr.predict(s, "flat_plate")["eta_overall"] / m - 1.0
        assert 0.16 < ch < 0.35, stem
        assert 0.08 < fp < 0.25, stem
        assert fp < ch - 0.07, stem


def test_the_flat_plate_rows_are_flagged_outside_the_kernels_span(series) -> None:
    """Re_x 1.6-5.9e4 against the Stanford data's 1.5e5 floor: evaluated,
    reported out-of-domain, never dropped."""
    for s in series.values():
        d = nr.predict(s, "flat_plate")
        assert d["h_gas_unblown_extrapolated"] == 1.0
        assert 1.5e4 < d["h_gas_unblown_Re_x_start"] < d["h_gas_unblown_Re_x_end"] < 6.0e4
        recs = [r for r in nr.run_series(s) if r.group == "eta_tot, flat_plate gas side"]
        assert recs and all(r.extrapolated for r in recs)
        recs = [r for r in nr.run_series(s) if r.group == "eta_tot, chamber gas side"]
        assert recs and not any(r.extrapolated for r in recs)


def test_only_the_window_is_scored(series) -> None:
    for s in series.values():
        for r in nr.run_series(s):
            assert (r.predicted is not None) == (18.0 <= r.x <= 50.0)


def test_the_trends(dataset) -> None:
    """Sign right on angle, dp and pressure for both baselines; the coolant
    temperature effect missed entirely by both (the closure is T_c-invariant
    to first order)."""
    angle, tc, dp, p3, p2 = nr.trends(dataset)  # in nr.TRENDS order
    assert angle["measured"] == pytest.approx(-0.073, abs=1e-3)
    assert tc["measured"] == pytest.approx(-0.090, abs=1e-3)
    assert dp["measured"] == pytest.approx(0.015, abs=1e-3)
    for b in (*nr.BASELINES, *(f"{x}+rad" for x in nr.BASELINES)):
        assert angle[b] < 0.0 and angle[b] / angle["measured"] < 0.3
        assert abs(tc[b]) < 2e-3
        assert dp[b] > 0.0
        assert p3[b] < 0.0 and p2[b] < 0.0
    # The flat plate overshoots the dp trend (1.5x) where the chamber
    # undershoots it (0.86x), and closes half the pressure trend's gap.
    assert dp["chamber"] / dp["measured"] == pytest.approx(0.86, abs=0.03)
    assert dp["flat_plate"] / dp["measured"] == pytest.approx(1.54, abs=0.05)
    assert p3["flat_plate"] / p3["measured"] == pytest.approx(0.68, abs=0.05)
    # Gas radiation (#1071) grows with pressure (eps_g 0.14 -> 0.26), so it
    # closes most of the flat plate's pressure trend; it also steepens the
    # dp trend, whose overshoot grows to 1.8x.
    assert p3["flat_plate+rad"] / p3["measured"] == pytest.approx(0.84, abs=0.05)
    assert dp["flat_plate+rad"] / dp["measured"] == pytest.approx(1.82, abs=0.06)


def test_the_start_length_moves_the_score_and_the_duct_area_does_not(series) -> None:
    """Falsifying the metric: a farther boundary-layer start lowers h and
    raises eta, monotonically. The 80e3 vs 8e3 mm^2 duct reading moves the
    chamber baseline only -- the flat plate never sees Dh."""
    for stem, s in series.items():
        v = nr.sensitivity(s)
        assert (
            v["flat_plate_from_sample_le"]
            < v["flat_plate_bleed_half"]
            < v["flat_plate"]
            < v["flat_plate_bleed_double"]
        ), stem
        assert v["flat_plate"] - v["flat_plate_from_sample_le"] > 0.01, stem
        assert v["chamber"] - v["chamber_dh_alt"] > 0.015, stem
        assert v["flat_plate_dh_alt"] == pytest.approx(v["flat_plate"], abs=1e-9), stem
        # The whole-sample panel averages a thicker boundary layer.
        assert v["flat_plate_whole_sample"] > v["flat_plate"], stem
        # Radiation: a shorter beam (80 x 100 mm duct) or a duller wall
        # radiates less, so eta rises toward the no-radiation value.
        assert v["flat_plate_rad"] < v["flat_plate_rad_eps_alt"] < v["flat_plate"], stem
        assert v["flat_plate_rad"] < v["flat_plate_rad_beam_alt"] < v["flat_plate"], stem


# Plateau predictions with gas radiation (rig beam length, eps_w 0.9):
# (chamber, flat plate), and the plate's hot-face radiative flux [kW/m^2].
PINNED_RAD = {
    "fig9_eta_tot_E20_dp2": (0.896, 0.832, 38.2),
    "fig9_eta_tot_E20_dp3": (0.907, 0.849, 38.5),
    "fig9_eta_tot_E20_dp4": (0.914, 0.860, 38.8),
    "fig6_eta_tot_E30_Tc470": (0.894, 0.829, 39.3),
    "fig6_eta_tot_E30_Tc540": (0.894, 0.829, 38.1),
    "fig6_eta_tot_E30_Tc580": (0.892, 0.827, 37.3),
    "fig10_eta_tot_E30_5bar_dp2": (0.882, 0.811, 38.9),
    "fig10_eta_tot_E30_10bar_dp3": (0.879, 0.808, 67.4),
    "fig10_eta_tot_E30_10bar_dp2": (0.868, 0.791, 66.8),
}


def test_the_radiation_rows_are_pinned_and_score_apart(series) -> None:
    """ACCURACY, cross-source, on their own rows: gas radiation lowers every
    curve by 2.5-3.5 points and leaves all of them over-predicted."""
    for stem, (chamber, flat, q_rad) in PINNED_RAD.items():
        s = series[stem]
        dc = nr.predict(s, "chamber", radiation=True)
        df = nr.predict(s, "flat_plate", radiation=True)
        assert dc["eta_overall"] == pytest.approx(chamber, abs=1e-3), stem
        assert df["eta_overall"] == pytest.approx(flat, abs=1e-3), stem
        assert df["q_rad_hot"] * 1e-3 == pytest.approx(q_rad, abs=0.2), stem
        assert df["eta_overall"] > nr.plateau(s), stem
        groups = {r.group for r in nr.run_series(s)}
        assert groups == {nr._group(b, rad) for b in nr.BASELINES for rad in nr.RADIATION}


# Plateau predictions with the Bordbar 2021 WSGG emissivity (#1071) in place
# of Eq. 8.11, same beam length and wall: (chamber, flat plate, q_rad hot
# face [kW/m^2], eps_g). ACCURACY, cross-source, on their own rows.
PINNED_RAD_BORDBAR = {
    "fig9_eta_tot_E20_dp2": (0.897, 0.833, 37.0, 0.132),
    "fig9_eta_tot_E20_dp3": (0.908, 0.850, 37.4, 0.132),
    "fig9_eta_tot_E20_dp4": (0.915, 0.861, 37.6, 0.132),
    "fig6_eta_tot_E30_Tc470": (0.895, 0.830, 38.1, 0.132),
    "fig6_eta_tot_E30_Tc540": (0.895, 0.830, 36.9, 0.132),
    "fig6_eta_tot_E30_Tc580": (0.893, 0.828, 36.2, 0.132),
    "fig10_eta_tot_E30_5bar_dp2": (0.883, 0.812, 37.7, 0.132),
    "fig10_eta_tot_E30_10bar_dp3": (0.888, 0.815, 51.2, 0.198),
    "fig10_eta_tot_E30_10bar_dp2": (0.877, 0.799, 50.7, 0.198),
}


def test_the_bordbar_radiation_rows_are_pinned(series, dataset) -> None:
    """At 5 bar the WSGG eps_g (0.132) sits 4 % under Eq. 8.11's (0.137) and
    eta barely moves (+0.1 point); at 10 bar it is 0.198 against 0.262, the
    radiative flux falls 24 % and eta rises 0.7-0.9 point. Every curve stays
    over-predicted, and the 5 -> 10 bar trend, which Eq. 8.11 closed to 0.84
    of measured, falls back to 0.58 on the flat plate."""
    for stem, (chamber, flat, q_rad, eps_g) in PINNED_RAD_BORDBAR.items():
        s = series[stem]
        dc = nr.predict(s, "chamber", radiation="bordbar2021")
        df = nr.predict(s, "flat_plate", radiation="bordbar2021")
        assert dc["eta_overall"] == pytest.approx(chamber, abs=1e-3), stem
        assert df["eta_overall"] == pytest.approx(flat, abs=1e-3), stem
        assert df["q_rad_hot"] * 1e-3 == pytest.approx(q_rad, abs=0.2), stem
        assert df["eps_g"] == pytest.approx(eps_g, abs=1e-3), stem
        assert df["eta_overall"] > nr.plateau(s), stem
        # Radiation still lowers eta, by less than Eq. 8.11's.
        lef = nr.predict(s, "flat_plate", radiation=True)["eta_overall"]
        assert lef < df["eta_overall"] < nr.predict(s, "flat_plate")["eta_overall"], stem
    p3 = nr.trends(dataset)[3]
    assert p3["flat_plate+rad_bordbar"] / p3["measured"] == pytest.approx(0.58, abs=0.05)
