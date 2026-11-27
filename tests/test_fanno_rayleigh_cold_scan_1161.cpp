#include <gtest/gtest.h>

#include "../include/fanno_mach.h"
#include "../include/humidair.h"

#include <cmath>
#include <string>
#include <vector>

using namespace combaero;

// ---------------------------------------------------------------------------
// The cold scan of fanno_rayleigh_channel_flow (no warm start) steps up in
// factors of 2 to the choke. The march sonics ~3e-7 below the choke
// relation's flux, so the last bracket's upper end always read as choked, and
// the scan gave up without looking inside it: every root above ~0.93 Gc was
// missed, and the kernel returned G = Gc unchoked with a nonzero dG/dP_target
// on a flat plateau. A warm start found the real root, so the flux depended on
// the solver's history (#1161, found under #1158).
// ---------------------------------------------------------------------------

namespace {

const std::string kModel = "haaland";

struct Duct {
  double Pt0, Tt0, Tt_exit, L, D;
  bool exit_total;
};

// The cold duct of stress hotcold seed 10 (the issue's reproduction), its
// static-exit twin, and the hot cooled ducts of seeds 11 and 12.
const Duct kIssue{1.0e6, 300.0, 697.9029088757375, 0.2, 0.008, true};
const Duct kIssueStatic{1.0e6, 300.0, 697.9029088757375, 0.2, 0.008, false};
const Duct kCooled{1.0e6, 1600.0, 937.4701685805791, 1.0, 0.006, true};
// The hot cooled duct of seed 12 on a static-pressure exit.
const Duct kCooledStatic{1.0e6, 1200.0, 824.0167669773003, 0.5, 0.006, false};

FannoChannelFlow flow(const Duct& d, double drive, double G_guess = -1.0,
                      bool with_derivatives = false, double Pt0_shift = 0.0,
                      double Tt0_shift = 0.0, double Tte_shift = 0.0) {
  const std::vector<double> X = dry_air();
  return fanno_rayleigh_channel_flow(d.Pt0 + Pt0_shift, d.Tt0 + Tt0_shift, d.Tt_exit + Tte_shift,
                                     X, d.Pt0 - drive, d.exit_total, d.L, d.D, 0.0, kModel, 1.0,
                                     G_guess, {}, with_derivatives);
}

double choked_flux(const Duct& d) {
  return fanno_rayleigh_choked_mass_flux(d.Pt0, d.Tt0, d.Tt_exit, dry_air(), d.L, d.D, 0.0,
                                         kModel);
}

// The drive at which the cold path first reports the duct choked.
double choking_drive(const Duct& d) {
  double a = 1.0, b = 0.9 * d.Pt0;
  for (int i = 0; i < 45; ++i) {
    const double m = 0.5 * (a + b);
    (flow(d, m).choked ? b : a) = m;
  }
  return b;
}

}  // namespace

TEST(FannoRayleighColdScan1161, IssueDuctFindsTheRootBelowTheChoke) {
  const double Gc = choked_flux(kIssue);
  // Continue the root up in drive from 140 kPa (0.91 Gc, which the old scan
  // found), each step warm-started from the last: the issue's warm column.
  // Every cold call must land on the same root; the old one returned Gc,
  // unchoked, from 160 to 195 kPa.
  double G_warm = flow(kIssue, 140e3).G;
  ASSERT_LT(G_warm, 0.93 * Gc);
  for (double drive = 142e3; drive <= 194e3; drive += 2e3) {
    const FannoChannelFlow w = flow(kIssue, drive, G_warm);
    const FannoChannelFlow c = flow(kIssue, drive);
    ASSERT_FALSE(w.choked) << drive;
    EXPECT_FALSE(c.choked) << drive;
    EXPECT_NEAR(c.G / w.G, 1.0, 1e-11) << drive;
    // Not the plateau: the root lies below the choked flux, the exit subsonic.
    EXPECT_LT(c.G, Gc * (1.0 - 1e-5)) << drive;
    EXPECT_LT(c.M_exit, 1.0) << drive;
    G_warm = w.G;
  }
}

TEST(FannoRayleighColdScan1161, FluxIsMonotoneAndContinuousToTheChoke) {
  for (const Duct& d : {kIssue, kIssueStatic, kCooled}) {
    const double Gc = choked_flux(d);
    const double dc = choking_drive(d);
    double G_prev = 0.0;
    double slope_prev = 0.0;
    const int n = 80;
    for (int i = 1; i <= n; ++i) {
      const double drive = 1.2 * dc * i / n;
      const FannoChannelFlow f = flow(d, drive);
      EXPECT_GE(f.G, G_prev) << "drive " << drive;
      if (!f.choked) {
        EXPECT_LT(f.G, Gc * (1.0 - 1e-7)) << "unchoked at the choked flux, drive " << drive;
      }
      // Continuity: no step much larger than the previous step's.
      const double step = f.G - G_prev;
      if (i > 2 && !f.choked) {
        EXPECT_LT(step, 3.0 * slope_prev + 1e-9 * Gc) << "jump at drive " << drive;
      }
      slope_prev = step;
      G_prev = f.G;
    }
    // Just below and above the choking drive the flux meets Gc.
    EXPECT_NEAR(flow(d, dc * (1.0 - 1e-6)).G / Gc, 1.0, 1e-4);
    EXPECT_DOUBLE_EQ(flow(d, dc * (1.0 + 1e-6)).G, Gc);
  }
}

TEST(FannoRayleighColdScan1161, ColdAndWarmStartsAgree) {
  for (const Duct& d : {kIssue, kIssueStatic, kCooled}) {
    const double dc = choking_drive(d);
    const double Gc = choked_flux(d);
    for (const double frac : {0.5, 0.8, 0.9, 0.95, 0.99, 0.999, 1.001, 1.1}) {
      const FannoChannelFlow c = flow(d, frac * dc);
      for (const double g : {0.05 * c.G, 0.5 * c.G, 0.999 * c.G, c.G, 1.001 * c.G, 0.999 * Gc,
                             1.5 * c.G}) {
        const FannoChannelFlow w = flow(d, frac * dc, g);
        EXPECT_NEAR(w.G / c.G, 1.0, 1e-12) << "drive " << frac << " dc, guess " << g;
        EXPECT_EQ(w.choked, c.choked) << "drive " << frac << " dc, guess " << g;
      }
    }
  }
}

TEST(FannoRayleighColdScan1161, DerivativesMatchDifferencesOnTheFormerPlateau) {
  // Steps of 1e-4 of the drive (>= 10 Pa here): clear of the C1 blend below
  // the 1 Pa kFannoFlowDriveFloor.
  for (const Duct& d : {kIssue, kIssueStatic, kCooled}) {
    const double dc = choking_drive(d);
    for (const double frac : {0.95, 0.99, 0.999}) {
      const double drive = frac * dc;
      const FannoChannelFlow f = flow(d, drive, -1.0, true);
      ASSERT_FALSE(f.choked);
      const double hP = 1e-4 * drive;
      const double fdP = (flow(d, drive - hP).G - flow(d, drive + hP).G) / (2.0 * hP);
      EXPECT_NEAR(f.dG_dP_target / fdP, 1.0, 1e-4) << "drive " << frac << " dc";
      const double fdPt =
          (flow(d, drive, -1.0, false, hP).G - flow(d, drive, -1.0, false, -hP).G) /
          (2.0 * hP);
      EXPECT_NEAR(f.dG_dPt0 / fdPt, 1.0, 1e-4) << "drive " << frac << " dc";
      const double hT = 1e-5 * d.Tt0;
      const double fdT0 = (flow(d, drive, -1.0, false, 0.0, hT).G -
                           flow(d, drive, -1.0, false, 0.0, -hT).G) /
                          (2.0 * hT);
      EXPECT_NEAR(f.dG_dTt0 * d.Tt0 / f.G, fdT0 * d.Tt0 / f.G, 1e-4) << "drive " << frac << " dc";
      const double hTe = 1e-5 * d.Tt_exit;
      const double fdTe = (flow(d, drive, -1.0, false, 0.0, 0.0, hTe).G -
                           flow(d, drive, -1.0, false, 0.0, 0.0, -hTe).G) /
                          (2.0 * hTe);
      EXPECT_NEAR(f.dG_dTt_exit * d.Tt_exit / f.G, fdTe * d.Tt_exit / f.G, 1e-4)
          << "drive " << frac << " dc";
    }
  }
}

TEST(FannoRayleighColdScan1161, TheChokeEdgeIsSharpWithFiniteDerivatives) {
  // Up to the choking drive the duct carries a subsonic root right up to the
  // march's own sonic flux (~3e-7 below Gc), with finite derivatives; past
  // it, Gc choked. A static exit's E(G) is so steep there that the fluxes
  // closest to the march's sonic point (zeta < kRayleighZetaFloor) span ~2 Pa
  // of drive: reported choked, not as Gc unchoked, and never as a root whose
  // cut derivative steps all round to the same exit pressure (dG/dP = inf).
  for (const Duct& d : {kIssue, kIssueStatic, kCooledStatic}) {
    const double Gc = choked_flux(d);
    const double dc = choking_drive(d);
    ASSERT_LT(dc, 0.7 * d.Pt0);
    const FannoChannelFlow far = flow(d, dc - 1000.0, -1.0, true);
    double G_prev = far.G;
    for (const double below : {100.0, 10.0, 1.0, 0.1, 0.01, 0.001}) {
      const FannoChannelFlow f = flow(d, dc - below, -1.0, true);
      EXPECT_FALSE(f.choked) << below << " Pa below";
      EXPECT_GE(f.G, G_prev) << below << " Pa below";
      EXPECT_LT(f.G, Gc) << below << " Pa below";  // the march sonics 7e-10-3e-7 below Gc
      EXPECT_LT(f.M_exit, 1.0) << below << " Pa below";
      EXPECT_TRUE(std::isfinite(f.dG_dP_target)) << below << " Pa below";
      EXPECT_LE(f.dG_dP_target, 0.0) << below << " Pa below";
      EXPECT_LE(std::abs(f.dG_dP_target), 2.0 * std::abs(far.dG_dP_target)) << below;
      EXPECT_TRUE(std::isfinite(f.dG_dPt0) && std::isfinite(f.dG_dTt0) &&
                  std::isfinite(f.dG_dTt_exit))
          << below << " Pa below";
      G_prev = f.G;
    }
    // The root reaches the march's sonic flux: within 1e-6 of Gc.
    EXPECT_GT(G_prev / Gc, 1.0 - 1e-6);
    for (const double above : {0.01, 1.0, 100.0}) {
      const FannoChannelFlow f = flow(d, dc + above, -1.0, true);
      EXPECT_TRUE(f.choked) << above << " Pa above";
      EXPECT_DOUBLE_EQ(f.G, Gc) << above << " Pa above";
      EXPECT_EQ(f.dG_dP_target, 0.0) << above << " Pa above";
    }
  }
}

TEST(FannoRayleighColdScan1161, StaticChokingDriveIsTheMarchsSubsonicEnd) {
  // A static exit chokes at the drive Pt0 - P_exit(G_s), G_s the largest flux
  // the march carries with zeta_min above the floor (1e-10). Reference: plain
  // bisection on the public march, to the last bit of log G. The kernel's
  // Illinois search must land within 0.1 Pa of it; a looser subsonic-end
  // tolerance (1e-3 in log G) moved the edge by 72-518 Pa, a floor of 1e-6
  // by 248-785 Pa. And within 10 Pa below the sonic-exit drive: the floor's
  // band is ~2 Pa on these ducts.
  const std::vector<double> X = dry_air();
  for (const Duct& d : {kIssueStatic, kCooledStatic}) {
    const double Gc = choked_flux(d);
    auto march = [&](double logG) {
      return fanno_rayleigh_march(std::exp(logG), d.Pt0, d.Tt0, d.Tt_exit, X, d.L, d.D, 0.0,
                                  kModel);
    };
    double lo = std::log(0.9 * Gc), hi = std::log(Gc);
    ASSERT_GT(march(lo).zeta_min, 1e-10);
    for (int i = 0; i < 200; ++i) {
      const double mid = 0.5 * (lo + hi);
      if (mid <= lo || mid >= hi) {
        break;
      }
      (march(mid).zeta_min > 1e-10 ? lo : hi) = mid;
    }
    const FannoRayleighMarch m_gc = march(std::log(Gc));
    const double ref = d.Pt0 - march(lo).P_exit;
    const double sonic = d.Pt0 - m_gc.P_exit;
    const double dc = choking_drive(d);
    EXPECT_NEAR(dc, ref, 0.1) << "Tt0 " << d.Tt0;
    if (m_gc.zeta_min <= 1e-10) {
      EXPECT_LT(dc, sonic) << "Tt0 " << d.Tt0;
      EXPECT_GT(dc, sonic - 10.0) << "Tt0 " << d.Tt0;
    } else {
      // When the forward march stays above the zeta floor all the way to Gc,
      // ref equals sonic (0 Pa floor cutback) and choking_drive's upper
      // bisection bracket (width ~2.6e-8 Pa) sits right above sonic.
      EXPECT_NEAR(dc, sonic, 1e-6) << "Tt0 " << d.Tt0;
    }
  }
}
