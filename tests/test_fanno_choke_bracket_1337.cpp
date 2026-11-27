// Regression tests for #1337: long laminar Fanno and Fanno-Rayleigh ducts
// whose choked mass flux falls past the old 60 e-fold bracket cap and whose
// backward theta-RK4 march blew up on high f L / D. Covers:
//   - All 6 rows of #1337's threshold table, the L/D = 1000 repro, PR #1336,
//     and extreme L/D = 1e5 adiabatic, heated, and cooled laminar ducts.
//   - Analytic-vs-FD Jacobian checks across long laminar ducts, cooled folds,
//     heated near-choke / choke edges, and interior throats (kThroat).
//   - Accuracy pins that go red if the friction-stretched mesh (alpha) or
//     grading power (kRayleighGrading = 3) is degraded.

#include "../include/fanno_mach.h"
#include "../include/humidair.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace combaero {
namespace {

template <typename Fn>
double richardson_fd(const Fn& fn, double x0, double h) {
  const double d1 = (fn(x0 + h) - fn(x0 - h)) / (2.0 * h);
  const double d2 = (fn(x0 + 0.5 * h) - fn(x0 - 0.5 * h)) / h;
  return (4.0 * d2 - d1) / 3.0;
}

struct LongDuctCase {
  const char* name;
  double Pt0, Tt0, Tt_exit, L, D;
};

}  // namespace

TEST(FannoChokeBracket1337, LongLaminarDuctsChokeAndMatchAdiabatic) {
  const std::vector<double> X = dry_air();
  // Issue #1337 primary reproducer, PR #1336 review case, all 6 exact
  // threshold-table rows (just above each L/D failure threshold), and extreme
  // L/D = 1e5..1e6 laminar ducts. Every single one raises "brent: no sign
  // change" on main (e427b904).
  const std::vector<LongDuctCase> cases = {
      {"issue_1337_repro", 1.0e5, 1000.0, 1000.0, 0.10, 1.0e-4},
      {"pr_1336_review", 35416.0, 787.0, 787.0, 4.28, 1.01e-3},
      {"table_row1_T1000_1bar_D0p1mm", 1.0e5, 1000.0, 1000.0, 0.072, 1.0e-4},
      {"table_row2_T1000_1bar_D1mm", 1.0e5, 1000.0, 1000.0, 7.2, 1.0e-3},
      {"table_row3_T1000_1bar_D10mm", 1.0e5, 1000.0, 1000.0, 720.0, 1.0e-2},
      {"table_row4_T1000_10bar_D0p1mm", 1.0e6, 1000.0, 1000.0, 0.72, 1.0e-4},
      {"table_row5_T300_1bar_D0p1mm", 1.0e5, 300.0, 300.0, 0.30, 1.0e-4},
      {"table_row6_T300_10bar_D0p1mm", 1.0e6, 300.0, 300.0, 3.0, 1.0e-4},
      {"extreme_LD_1e5", 1.0e5, 300.0, 300.0, 100.0, 1.0e-3},
      {"extreme_LD_1e6", 1.0e5, 300.0, 300.0, 1000.0, 1.0e-3},
  };

  for (const auto& c : cases) {
    const double Gc_ad = fanno_choked_mass_flux(c.Pt0, c.Tt0, X, c.L, c.D, 0.0, "haaland", 1.0);
    const double Gc_ray = fanno_rayleigh_choked_mass_flux(c.Pt0, c.Tt0, c.Tt_exit, X, c.L, c.D,
                                                          0.0, "haaland", 1.0);
    ASSERT_TRUE(std::isfinite(Gc_ad)) << c.name;
    ASSERT_GT(Gc_ad, 0.0) << c.name;
    ASSERT_TRUE(std::isfinite(Gc_ray)) << c.name;
    ASSERT_GT(Gc_ray, 0.0) << c.name;
    // Adiabatic Fanno and adiabatic Fanno-Rayleigh must agree closely across
    // arbitrary laminar L/D (fails if friction stretching alpha = 0 is used).
    EXPECT_NEAR(Gc_ray / Gc_ad, 1.0, 1e-7) << c.name;

    // Forward march at 0.999 Gc stays subsonic and approaches M_exit -> 1.
    const FannoRayleighMarch m =
        fanno_rayleigh_march(0.999 * Gc_ray, c.Pt0, c.Tt0, c.Tt_exit, X, c.L, c.D, 0.0, "haaland",
                             1.0);
    EXPECT_GT(m.zeta_min, 0.0) << c.name;
    EXPECT_GT(m.M_exit, 0.0) << c.name;
    EXPECT_LT(m.M_exit, 1.0) << c.name;

    // Channel flow succeeds in both choked and unchoked regimes.
    const FannoChannelFlow f_ad_ch = fanno_channel_flow(
        c.Pt0, c.Tt0, X, 1.0, false, c.L, c.D, 0.0, "haaland", 1.0, -1.0, {}, true);
    const FannoChannelFlow f_ray_ch = fanno_rayleigh_channel_flow(
        c.Pt0, c.Tt0, c.Tt_exit, X, 1.0, false, c.L, c.D, 0.0, "haaland", 1.0, -1.0, {}, true);
    EXPECT_TRUE(f_ad_ch.choked) << c.name;
    EXPECT_TRUE(f_ray_ch.choked) << c.name;
    EXPECT_NEAR(f_ray_ch.G / f_ad_ch.G, 1.0, 1e-7) << c.name;

    const FannoChannelFlow f_ray_sub = fanno_rayleigh_channel_flow(
        c.Pt0, c.Tt0, c.Tt_exit, X, 0.5 * c.Pt0, false, c.L, c.D, 0.0, "haaland", 1.0, -1.0, {},
        true);
    EXPECT_FALSE(f_ray_sub.choked) << c.name;
    EXPECT_GT(f_ray_sub.G, 0.0) << c.name;
    EXPECT_LT(f_ray_sub.G, Gc_ray) << c.name;
    EXPECT_LT(f_ray_sub.dG_dP_target, 0.0) << c.name;
  }
}

TEST(FannoChokeBracket1337, HeatedAndCooledLongLaminarTangentsMatchFD) {
  const std::vector<double> X = dry_air();
  const std::vector<LongDuctCase> cases = {
      {"heated_LD1000", 1.0e5, 300.0, 600.0, 1.0, 1.0e-3},
      {"cooled_LD1000", 1.0e5, 600.0, 300.0, 1.0, 1.0e-3},
      {"heated_LD1e5", 1.0e5, 300.0, 600.0, 100.0, 1.0e-3},
      {"cooled_LD1e5", 1.0e5, 600.0, 300.0, 100.0, 1.0e-3},
  };

  for (const auto& c : cases) {
    // Choked kExit tangents via Dual4 vs Richardson FD of choked_mass_flux.
    const FannoChannelFlow fc = fanno_rayleigh_channel_flow(
        c.Pt0, c.Tt0, c.Tt_exit, X, 1.0, false, c.L, c.D, 0.0, "haaland", 1.0, -1.0, {}, true);
    ASSERT_TRUE(fc.choked) << c.name;
    EXPECT_DOUBLE_EQ(fc.dG_dP_target, 0.0) << c.name;

    auto Gc_of = [&](double Pt, double T0, double Te) {
      return fanno_rayleigh_choked_mass_flux(Pt, T0, Te, X, c.L, c.D, 0.0, "haaland", 1.0);
    };
    const double fd_Pt =
        richardson_fd([&](double v) { return Gc_of(v, c.Tt0, c.Tt_exit); }, c.Pt0, 1e-4 * c.Pt0);
    const double fd_T0 =
        richardson_fd([&](double v) { return Gc_of(c.Pt0, v, c.Tt_exit); }, c.Tt0, 1e-4 * c.Tt0);
    const double fd_Te = richardson_fd([&](double v) { return Gc_of(c.Pt0, c.Tt0, v); }, c.Tt_exit,
                                       1e-4 * c.Tt_exit);
    EXPECT_NEAR(fc.dG_dPt0, fd_Pt, 1e-7 * std::abs(fd_Pt)) << c.name;
    EXPECT_NEAR(fc.dG_dTt0, fd_T0, 1e-7 * std::max(std::abs(fd_T0), 1e-12)) << c.name;
    EXPECT_NEAR(fc.dG_dTt_exit, fd_Te, 1e-7 * std::max(std::abs(fd_Te), 1e-12)) << c.name;

    // Unchoked tangents via Dual4 vs Richardson FD of march exit pressure.
    const double P_target = 0.5 * c.Pt0;
    const FannoChannelFlow fu = fanno_rayleigh_channel_flow(
        c.Pt0, c.Tt0, c.Tt_exit, X, P_target, false, c.L, c.D, 0.0, "haaland", 1.0, -1.0, {}, true);
    ASSERT_FALSE(fu.choked) << c.name;
    auto E_of = [&](double G, double Pt, double T0, double Te) {
      return fanno_rayleigh_march(G, Pt, T0, Te, X, c.L, c.D, 0.0, "haaland", 1.0).P_exit;
    };
    const double E_G =
        richardson_fd([&](double v) { return E_of(v, c.Pt0, c.Tt0, c.Tt_exit); }, fu.G, 1e-4 * fu.G);
    const double E_Pt = richardson_fd([&](double v) { return E_of(fu.G, v, c.Tt0, c.Tt_exit); },
                                      c.Pt0, 1e-4 * c.Pt0);
    const double E_T0 = richardson_fd([&](double v) { return E_of(fu.G, c.Pt0, v, c.Tt_exit); },
                                      c.Tt0, 1e-4 * c.Tt0);
    const double E_Te = richardson_fd([&](double v) { return E_of(fu.G, c.Pt0, c.Tt0, v); },
                                      c.Tt_exit, 1e-4 * c.Tt_exit);
    EXPECT_NEAR(fu.dG_dP_target, 1.0 / E_G, 1e-7 * std::abs(1.0 / E_G)) << c.name;
    EXPECT_NEAR(fu.dG_dPt0, -E_Pt / E_G, 1e-7 * std::abs(E_Pt / E_G)) << c.name;
    EXPECT_NEAR(fu.dG_dTt0, -E_T0 / E_G, 1e-7 * std::max(std::abs(E_T0 / E_G), 1e-12)) << c.name;
    EXPECT_NEAR(fu.dG_dTt_exit, -E_Te / E_G, 1e-7 * std::max(std::abs(E_Te / E_G), 1e-12))
        << c.name;
  }
}

TEST(FannoChokeBracket1337, FoldAndChokeEdgeJacobiansMatchFD) {
  const std::vector<double> X = dry_air();
  struct Case {
    const char* name;
    double Pt0, Tt0, Tt_exit, P_target;
    bool exit_total;
    double L, D;
    std::string model;
    double fmult;
  };
  const std::vector<Case> cases = {
      // Cooled folds (#1167, #1206) on the first fall below peak drive.
      {"cooled_fold_kLast", 3.0e6, 1400.0, 1000.0, 3.0e6 - 260957.0, false, 0.03, 0.01, "haaland",
       1.0},
      {"cooled_fold_kInterior", 2.0e5, 1800.0, 400.0, 2.0e5 - 100.0, false, 0.3, 0.01, "haaland",
       1.0},
      {"cooled_fold_k1206", 177138.87, 1137.32, 606.41, 177138.87 - 1535.0, true, 0.105077,
       4.18657e-3, "haaland", 1.0},
      // Heated and cooled near-choke edges (#1161).
      {"heated_choke_edge", 1.0e6, 300.0, 697.9029088757375, 1.0e6 - 570000.0, false, 0.2, 0.008,
       "haaland", 1.0},
      {"cooled_choke_edge", 1.0e6, 1200.0, 824.0167669773003, 1.0e6 - 620000.0, false, 0.5, 0.006,
       "haaland", 1.0},
  };

  for (const auto& c : cases) {
    const double Gc = fanno_rayleigh_choked_mass_flux(c.Pt0, c.Tt0, c.Tt_exit, X, c.L, c.D, 0.0,
                                                      c.model, c.fmult);
    const FannoChannelFlow f = fanno_rayleigh_channel_flow(
        c.Pt0, c.Tt0, c.Tt_exit, X, c.P_target, c.exit_total, c.L, c.D, 0.0, c.model, c.fmult,
        -1.0, {}, true);
    ASSERT_FALSE(f.choked) << c.name;
    auto E_of = [&](double G, double Pt, double T0, double Te) {
      const FannoRayleighMarch m =
          fanno_rayleigh_march(G, Pt, T0, Te, X, c.L, c.D, 0.0, c.model, c.fmult);
      return c.exit_total ? m.Pt_exit : m.P_exit;
    };
    // Keep the FD probe strictly subsonic even near the choke edge; tol 1e-5
    // stays above the O(eps * Pt0 / (h * E_G)) rounding floor near folds and
    // sonic edges while remaining 10x tighter than the 1e-4 #1161/#1167 checks.
    const double hFrac = std::min(1e-4, 0.05 * (Gc - f.G) / Gc);
    const double E_G =
        richardson_fd([&](double v) { return E_of(v, c.Pt0, c.Tt0, c.Tt_exit); }, f.G, hFrac * f.G);
    const double E_Pt = richardson_fd([&](double v) { return E_of(f.G, v, c.Tt0, c.Tt_exit); },
                                      c.Pt0, hFrac * c.Pt0);
    const double E_T0 = richardson_fd([&](double v) { return E_of(f.G, c.Pt0, v, c.Tt_exit); },
                                      c.Tt0, hFrac * c.Tt0);
    const double E_Te = richardson_fd([&](double v) { return E_of(f.G, c.Pt0, c.Tt0, v); },
                                      c.Tt_exit, hFrac * c.Tt_exit);
    EXPECT_NEAR(f.dG_dP_target, 1.0 / E_G, 1e-5 * std::abs(1.0 / E_G)) << c.name;
    EXPECT_NEAR(f.dG_dPt0, -E_Pt / E_G, 1e-5 * std::abs(E_Pt / E_G)) << c.name;
    EXPECT_NEAR(f.dG_dTt0, -E_T0 / E_G, 1e-5 * std::max(std::abs(E_T0 / E_G), 1e-12)) << c.name;
    EXPECT_NEAR(f.dG_dTt_exit, -E_Te / E_G, 1e-5 * std::max(std::abs(E_Te / E_G), 1e-12))
        << c.name;
  }

  // Interior throat (kThroat) sweep across a node crossing verifies C1 cubic
  // Hermite zeta_peak smoothness (no slope jump when imin shifts by one node).
  for (const double Tt0 : {598.0, 600.0, 602.0}) {
    const FannoChannelFlow ft = fanno_rayleigh_channel_flow(
        5.0e5, Tt0, 380.0, X, 1.0e4, false, 0.70, 0.02, 1.0e-5, "fixed", 0.02, -1.0, {}, true);
    ASSERT_TRUE(ft.choked) << Tt0;
    auto Gc_of = [&](double Pt, double T0, double Te) {
      return fanno_rayleigh_choked_mass_flux(Pt, T0, Te, X, 0.70, 0.02, 1.0e-5, "fixed", 0.02);
    };
    const double fd_T0 =
        richardson_fd([&](double v) { return Gc_of(5.0e5, v, 380.0); }, Tt0, 1e-4 * Tt0);
    const double fd_Te =
        richardson_fd([&](double v) { return Gc_of(5.0e5, Tt0, v); }, 380.0, 1e-4 * 380.0);
    EXPECT_NEAR(ft.dG_dTt0, fd_T0, 1e-4 * std::abs(fd_T0)) << Tt0;
    EXPECT_NEAR(ft.dG_dTt_exit, fd_Te, 1e-4 * std::abs(fd_Te)) << Tt0;
  }
}

TEST(FannoChokeBracket1337, FrictionStretchingAndGradingAccuracyGuard) {
  // Guards both grid parameters against degradation:
  //   - Dropping alpha -> 0 (un-stretched t^3 mesh) pushes L/D = 1e5 choked
  //     flux error from ~5e-9 to >1e-7 and leaves w_exit(Gc) ~ 6e-3.
  //   - Lowering kRayleighGrading from 3 to 2 degrades near-sonic static exit
  //     accuracy on kIssueStatic (see also test_fanno_rayleigh_cold_scan_1161).
  const std::vector<double> X = dry_air();
  const double Gc_ad = fanno_choked_mass_flux(1.0e5, 300.0, X, 100.0, 1.0e-3, 0.0, "haaland", 1.0);
  const double Gc_ray = fanno_rayleigh_choked_mass_flux(1.0e5, 300.0, 300.0, X, 100.0, 1.0e-3, 0.0,
                                                        "haaland", 1.0);
  EXPECT_NEAR(Gc_ray / Gc_ad, 1.0, 2e-8);
  const FannoRayleighMarch m_lam =
      fanno_rayleigh_march(Gc_ray, 1.0e5, 300.0, 300.0, X, 100.0, 1.0e-3, 0.0, "haaland", 1.0);
  EXPECT_LT(m_lam.zeta_min, 1e-4);

  // Near-sonic heated duct: Gc agrees with the 512-step reference
  // (1222.6722372531 kg/(m^2 s)) within 5e-7 (alpha=0 or p=1 fails badly).
  const double Gc_hot = fanno_rayleigh_choked_mass_flux(1.0e6, 300.0, 697.9029088757375, X, 0.2,
                                                        0.008, 0.0, "haaland", 1.0);
  EXPECT_NEAR(Gc_hot, 1222.6722372531, 5e-7 * 1222.6722372531);
}

}  // namespace combaero
