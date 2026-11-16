#include <gtest/gtest.h>

#include "../include/fanno_mach.h"
#include "../include/humidair.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

using namespace combaero;

// ---------------------------------------------------------------------------
// The cold scan of fanno_rayleigh_channel_flow (no warm start) walks log G up
// in ln-2 brackets and takes the first bracket whose excess (exit pressure
// minus the target) falls through zero. On a cooled duct the drive the duct
// needs can peak and turn back down, so the excess falls below zero and rises
// above it again inside one bracket. Both ends then read positive, the scan
// saw no fall, and the duct was reported choked while a warm start found its
// first crossing, unchoked (#1167; the last-bracket form of #1161's defect).
// ---------------------------------------------------------------------------

namespace {

const std::string kModel = "haaland";

struct Duct {
  double Pt0, Tt0, Tt_exit, L, D;
};

// The issue's reproductions, static-pressure exit: the fold in the last
// bracket (and its milder twins at Tt_exit 1050-1150 K), and the fold in an
// interior bracket (whose needed drive then falls below zero).
const Duct kLast{3.0e6, 1400.0, 1000.0, 0.03, 0.01};
const Duct kLast1050{3.0e6, 1400.0, 1050.0, 0.03, 0.01};
const Duct kLast1150{3.0e6, 1400.0, 1150.0, 0.03, 0.01};
const Duct kInterior{2.0e5, 1800.0, 400.0, 0.3, 0.01};

FannoChannelFlow flow(const Duct& d, double drive, double G_guess = -1.0,
                      bool with_derivatives = false, double Pt0_shift = 0.0,
                      double Tt0_shift = 0.0, double Tte_shift = 0.0, bool exit_total = false) {
  const std::vector<double> X = dry_air();
  return fanno_rayleigh_channel_flow(d.Pt0 + Pt0_shift, d.Tt0 + Tt0_shift, d.Tt_exit + Tte_shift,
                                     X, d.Pt0 - drive, exit_total, d.L, d.D, 0.0, kModel, 1.0,
                                     G_guess, {}, with_derivatives);
}

double choked_flux(const Duct& d) {
  return fanno_rayleigh_choked_mass_flux(d.Pt0, d.Tt0, d.Tt_exit, dry_air(), d.L, d.D, 0.0,
                                         kModel);
}

// The drive the duct needs at flux G: Pt0 - P_exit(G) of the public march
// (Pt_exit for a total-pressure exit).
double needed(const Duct& d, double G, bool exit_total = false) {
  const FannoRayleighMarch m =
      fanno_rayleigh_march(G, d.Pt0, d.Tt0, d.Tt_exit, dry_air(), d.L, d.D, 0.0, kModel);
  return d.Pt0 - (exit_total ? m.Pt_exit : m.P_exit);
}

// The fold: the flux in [a, b] Gc at which the needed drive peaks, and the
// peak, by golden section after a coarse sample.
std::pair<double, double> fold(const Duct& d, double a, double b) {
  const double Gc = choked_flux(d);
  auto f = [&](double s) { return -needed(d, s * Gc); };
  double best = a;
  for (int i = 0; i <= 200; ++i) {
    const double s = a + (b - a) * i / 200.0;
    if (f(s) < f(best)) {
      best = s;
    }
  }
  double lo = std::max(a, best - (b - a) / 200.0), hi = std::min(b, best + (b - a) / 200.0);
  const double g = 0.5 * (std::sqrt(5.0) - 1.0);
  for (int i = 0; i < 80; ++i) {
    const double c = hi - g * (hi - lo), e = lo + g * (hi - lo);
    (f(c) < f(e) ? hi : lo) = (f(c) < f(e) ? e : c);
  }
  const double s = 0.5 * (lo + hi);
  return {s * Gc, -f(s)};
}

}  // namespace

TEST(FannoRayleighFold1167, IssueCasesFindTheFirstCrossing) {
  const double Gc = choked_flux(kLast);
  // The issue's table: cold read Gc choked, a warm start from 0.9915 Gc the
  // root. Both must give the root now.
  const double drives[] = {260294.0, 260957.0, 261621.0};
  const double roots[] = {0.97537, 0.97971, 0.98554};
  for (int i = 0; i < 3; ++i) {
    const FannoChannelFlow c = flow(kLast, drives[i]);
    const FannoChannelFlow w = flow(kLast, drives[i], 0.9915 * Gc);
    EXPECT_FALSE(c.choked) << drives[i];
    EXPECT_FALSE(w.choked) << drives[i];
    EXPECT_NEAR(c.G / w.G, 1.0, 1e-11) << drives[i];
    EXPECT_NEAR(c.G / Gc, roots[i], 1e-5) << drives[i];
    EXPECT_NEAR(needed(kLast, c.G), drives[i], 1e-6 * drives[i]) << drives[i];
  }
  // The interior fold: at 100 Pa the issue's warm start from 0.12 Gc found
  // 0.12215 Gc, the cold scan Gc choked.
  const double Gi = choked_flux(kInterior);
  const FannoChannelFlow c = flow(kInterior, 100.0);
  EXPECT_FALSE(c.choked);
  EXPECT_NEAR(c.G / Gi, 0.12215, 1e-5);
  EXPECT_NEAR(c.G / flow(kInterior, 100.0, 0.12 * Gi).G, 1.0, 1e-11);
}

TEST(FannoRayleighFold1167, ColdAndWarmAgreeThroughEachFold) {
  // Guesses over 0.05-1.5 G at drives from below each fold's band up past its
  // peak. The interior duct's sweep starts at 25 Pa, above its laminar fold
  // (see LaminarFoldColdScanFindsTheFirstFall), and its fluxes agree to
  // 1e-9: at 2 bar and a few Pa of drive the root carries ~2e-10 of rounding
  // whichever way it is bracketed (the same before #1167).
  struct Case {
    Duct d;
    double a, b;  // where the fold's peak lies, in Gc
    double from;  // lowest drive [Pa]; <= 0: just below the drive at the choke
    double tol;
  };
  for (const Case& k :
       {Case{kLast, 0.95, 0.9999, 0.0, 1e-11}, Case{kLast1050, 0.99, 0.99999, 0.0, 1e-11},
        Case{kLast1150, 0.999, 0.9999999, 0.0, 1e-11}, Case{kInterior, 0.12, 0.2, 25.0, 1e-9}}) {
    const double Gc = choked_flux(k.d);
    const auto [G_peak, D_peak] = fold(k.d, k.a, k.b);
    const double D_end = needed(k.d, Gc * (1.0 - 1e-6));
    ASSERT_LT(D_end, D_peak);
    const double lo = (k.from > 0.0) ? k.from : D_end - 0.1 * (D_peak - D_end);
    for (int i = 0; i <= 16; ++i) {
      const double drive = lo + (1.02 * D_peak - lo) * i / 16.0;
      const FannoChannelFlow c = flow(k.d, drive);
      if (drive < D_peak) {
        // The first fall from zero flow: below the fold's peak, unchoked.
        EXPECT_FALSE(c.choked) << "Tt_exit " << k.d.Tt_exit << " drive " << drive;
        EXPECT_LT(c.G, G_peak) << "Tt_exit " << k.d.Tt_exit << " drive " << drive;
      } else {
        EXPECT_TRUE(c.choked) << "Tt_exit " << k.d.Tt_exit << " drive " << drive;
      }
      for (const double s : {0.05, 0.2, 0.5, 0.9, 0.99, 1.0, 1.01, 1.1, 1.5}) {
        const FannoChannelFlow w = flow(k.d, drive, s * c.G);
        EXPECT_NEAR(w.G / c.G, 1.0, k.tol)
            << "Tt_exit " << k.d.Tt_exit << " drive " << drive << " guess " << s;
        EXPECT_EQ(w.choked, c.choked)
            << "Tt_exit " << k.d.Tt_exit << " drive " << drive << " guess " << s;
      }
    }
  }
}

TEST(FannoRayleighFold1167, LaminarFoldColdScanFindsTheFirstFall) {
  // The interior duct's needed drive also peaks in its laminar range and
  // dips before the transition raises it again (static exit: 19.7 Pa near
  // 0.03 Gc, then 18.4 Pa near 0.047 Gc; total-pressure exit: 14.0 Pa near
  // 0.02 Gc, then below zero). Between those drives the excess falls, rises
  // and falls within a factor of 1.6 of flux, and the cold scan must take the
  // first fall. On the total-pressure exit both crossings lay in one bracket
  // from 13.02 Pa up to the peak, and the old scan read the duct choked. Reference: the
  // first flux on a fine grid whose needed drive reaches the target, refined
  // by bisection.
  const double Gc = choked_flux(kInterior);
  struct Case {
    bool exit_total;
    double drive;
  };
  for (const Case& k : {Case{false, 18.6}, Case{false, 18.9}, Case{false, 19.2},
                        Case{false, 19.5}, Case{false, 19.65}, Case{true, 13.0},
                        Case{true, 13.5}, Case{true, 13.7}, Case{true, 13.9}}) {
    double a = 0.002 * Gc, b = a;
    ASSERT_LT(needed(kInterior, a, k.exit_total), k.drive);
    for (int i = 1; i <= 4000; ++i) {
      b = Gc * (0.002 + 0.06 * i / 4000.0);
      if (needed(kInterior, b, k.exit_total) >= k.drive) {
        break;
      }
      a = b;
    }
    for (int i = 0; i < 60; ++i) {
      const double m = 0.5 * (a + b);
      (needed(kInterior, m, k.exit_total) >= k.drive ? b : a) = m;
    }
    const FannoChannelFlow c = flow(kInterior, k.drive, -1.0, false, 0.0, 0.0, 0.0, k.exit_total);
    EXPECT_FALSE(c.choked) << k.exit_total << " " << k.drive;
    EXPECT_NEAR(c.G / b, 1.0, 1e-8) << k.exit_total << " " << k.drive;
    // Warm starts agree, from either side of the first fall and from the
    // rise after it: a cooled duct's warm bracket widens to at most 6.4% of
    // its guess, and the cold scan decides the rest.
    for (const double s : {0.05, 0.2, 0.5, 0.9, 0.99, 1.0, 1.01, 1.1, 1.3, 1.5}) {
      const FannoChannelFlow w =
          flow(kInterior, k.drive, s * c.G, false, 0.0, 0.0, 0.0, k.exit_total);
      // 1e-8: near a fold's peak the excess is flat in G, and the root's
      // rounding (~2e-10 at these drives) is amplified.
      EXPECT_NEAR(w.G / c.G, 1.0, 1e-8) << k.exit_total << " " << k.drive << " guess " << s;
    }
  }
}

TEST(FannoRayleighFold1167, FluxRisesContinuouslyToTheFoldThenJumpsToTheChoke) {
  // The interior duct has a second, smaller fold at its laminar-turbulent
  // transition (the needed drive peaks at 19.7 Pa near 0.03 Gc and dips to
  // 18.4 Pa near 0.047 Gc): its sweep starts above it, at 25 Pa.
  for (const Duct& d : {kLast, kInterior}) {
    const double Gc = choked_flux(d);
    const auto [G_peak, D_peak] = fold(d, d.Pt0 > 1e6 ? 0.95 : 0.12, d.Pt0 > 1e6 ? 0.9999 : 0.2);
    const double D_from = d.Pt0 > 1e6 ? 0.0 : 25.0;
    // Up to the peak G rises with the drive, without a jump, toward the
    // peak's flux; there the needed drive is flat, so G goes as sqrt.
    double G_prev = 0.0;
    double step_prev = 0.0;
    const int n = 60;
    for (int i = 1; i <= n; ++i) {
      const double drive =
          D_from + (D_peak - D_from) *
                       (1.0 - std::pow(1.0 - static_cast<double>(i) / n, 2.0) - 1e-9);
      const FannoChannelFlow f = flow(d, drive);
      ASSERT_FALSE(f.choked) << "Tt_exit " << d.Tt_exit << " drive " << drive;
      EXPECT_GT(f.G, G_prev) << "Tt_exit " << d.Tt_exit << " drive " << drive;
      EXPECT_LT(f.G, G_peak * (1.0 + 1e-9)) << "drive " << drive;
      if (i > 3) {
        // The drive steps shrink linearly toward the peak while G's slope
        // grows as 1 / sqrt(D_peak - drive): a step at most a few times the
        // last is continuous; a jump to another branch is not.
        EXPECT_LT(f.G - G_prev, 4.0 * step_prev + 1e-6 * Gc) << "jump at drive " << drive;
      }
      step_prev = f.G - G_prev;
      G_prev = f.G;
    }
    EXPECT_GT(G_prev / G_peak, 1.0 - 1e-3);
    // Just past the peak no flux below the choke matches the drive: choked.
    const FannoChannelFlow past = flow(d, D_peak * (1.0 + 1e-6));
    EXPECT_TRUE(past.choked);
    EXPECT_DOUBLE_EQ(past.G, Gc);
  }
}

TEST(FannoRayleighFold1167, DerivativesAtTheNewRootsMatchDifferences) {
  // The roots the old scan skipped. Steps of 1e-3 of the distance to the
  // fold's peak (G goes as sqrt there), so a step never crosses it.
  for (const Duct& d : {kLast, kInterior}) {
    const auto [G_peak, D_peak] = fold(d, d.Pt0 > 1e6 ? 0.95 : 0.12, d.Pt0 > 1e6 ? 0.9999 : 0.2);
    const double band = d.Pt0 > 1e6 ? 2000.0 : 10.0;
    for (const double frac : {0.2, 0.5, 0.9}) {
      const double drive = D_peak - band * (1.0 - frac);
      const FannoChannelFlow f = flow(d, drive, -1.0, true);
      ASSERT_FALSE(f.choked) << drive;
      ASSERT_GT(f.G, 0.5 * G_peak);
      const double h = 1e-3 * (D_peak - drive);
      const double fdP = (flow(d, drive - h).G - flow(d, drive + h).G) / (2.0 * h);
      EXPECT_NEAR(f.dG_dP_target / fdP, 1.0, 1e-4) << "drive " << drive;
      // flow() holds the target at d.Pt0 - drive while it shifts Pt0.
      const double fdPt =
          (flow(d, drive, -1.0, false, h).G - flow(d, drive, -1.0, false, -h).G) / (2.0 * h);
      EXPECT_NEAR(f.dG_dPt0 / fdPt, 1.0, 1e-4) << "drive " << drive;
      const double hT = 1e-7 * d.Tt0;
      const double fdT0 = (flow(d, drive, -1.0, false, 0.0, hT).G -
                           flow(d, drive, -1.0, false, 0.0, -hT).G) /
                          (2.0 * hT);
      EXPECT_NEAR(f.dG_dTt0 / fdT0, 1.0, 1e-4) << "drive " << drive;
      const double hTe = 1e-7 * d.Tt_exit;
      const double fdTe = (flow(d, drive, -1.0, false, 0.0, 0.0, hTe).G -
                           flow(d, drive, -1.0, false, 0.0, 0.0, -hTe).G) /
                          (2.0 * hTe);
      EXPECT_NEAR(f.dG_dTt_exit / fdTe, 1.0, 1e-4) << "drive " << drive;
    }
  }
}
