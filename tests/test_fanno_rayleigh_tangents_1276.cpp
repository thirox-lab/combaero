// Regression tests for #1276: exact forward-mode discrete tangents (DualN<4>)
// through the 32-step RK4 Fanno-Rayleigh march, plus d2cp_R_dT2 / d2cp_dT2.

#include "../include/composition.h"
#include "../include/fanno_mach.h"
#include "../include/thermo.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace combaero {
namespace {

std::vector<double> dry_air() {
  std::vector<double> X(num_species(), 0.0);
  X[species_index_from_name("N2")] = 0.78084;
  X[species_index_from_name("O2")] = 0.20946;
  X[species_index_from_name("AR")] = 0.00934;
  X[species_index_from_name("CO2")] = 0.00036;
  return X;
}

std::vector<double> hot_products() {
  std::vector<double> X(num_species(), 0.0);
  X[species_index_from_name("N2")] = 0.72;
  X[species_index_from_name("CO2")] = 0.08;
  X[species_index_from_name("H2O")] = 0.16;
  X[species_index_from_name("O2")] = 0.04;
  return X;
}

template <typename Fn>
double richardson_fd(const Fn& fn, double x0, double h) {
  const double d1 = (fn(x0 + h) - fn(x0 - h)) / (2.0 * h);
  const double d2 = (fn(x0 + 0.5 * h) - fn(x0 - 0.5 * h)) / h;
  return (4.0 * d2 - d1) / 3.0;
}

// Symmetric degree-5 odd-polynomial least-squares derivative at x0: fits
// c1 * u + c3 * u^3 + c5 * u^5 to 0.5 * (fn(x0 + u*span) - fn(x0 - u*span))
// on u_k = k / M (k = 1..M) and returns c1 / span. Gives O(span^6) truncation
// while averaging down low-drive march root jitter (~2.5e-9 Pa on 2e5 Pa).
template <typename Fn>
double odd_poly5_deriv(const Fn& fn, double x0, double span, int M = 20) {
  double a[3][3] = {};
  double b[3] = {};
  for (int k = 1; k <= M; ++k) {
    const double u = static_cast<double>(k) / static_cast<double>(M);
    const double p[3] = {u, u * u * u, u * u * u * u * u};
    const double y = 0.5 * (fn(x0 + u * span) - fn(x0 - u * span));
    for (int i = 0; i < 3; ++i) {
      b[i] += p[i] * y;
      for (int j = 0; j < 3; ++j) {
        a[i][j] += p[i] * p[j];
      }
    }
  }
  // 3x3 Gaussian elimination with partial pivoting.
  for (int col = 0; col < 3; ++col) {
    int piv = col;
    for (int r = col + 1; r < 3; ++r) {
      if (std::abs(a[r][col]) > std::abs(a[piv][col])) {
        piv = r;
      }
    }
    if (piv != col) {
      for (int c = col; c < 3; ++c) {
        std::swap(a[col][c], a[piv][c]);
      }
      std::swap(b[col], b[piv]);
    }
    for (int r = col + 1; r < 3; ++r) {
      const double f = a[r][col] / a[col][col];
      for (int c = col; c < 3; ++c) {
        a[r][c] -= f * a[col][c];
      }
      b[r] -= f * b[col];
    }
  }
  double c[3] = {};
  for (int r = 2; r >= 0; --r) {
    double s = b[r];
    for (int col = r + 1; col < 3; ++col) {
      s -= a[r][col] * c[col];
    }
    c[r] = s / a[r][r];
  }
  return c[0] / span;
}

TEST(FannoRayleighTangents1276, SecondDerivativeOfCpMatchesFDAllRegions) {
  // Check every species across low tail, NASA-9 low, C2 bridge, NASA-9 high,
  // high tail fade, and constant high tail.
  const std::vector<double> temps = {60.0, 150.0, 350.0, 975.0, 1000.0, 1025.0,
                                     1800.0, 5950.0, 6200.0, 6800.0, 21000.0, 25000.0};
  for (std::size_t s = 0; s < num_species(); ++s) {
    for (double T : temps) {
      const double h = 1e-4 * T;
      const double ana = d2cp_R_dT2(s, T);
      const double fd = richardson_fd([&](double t) { return dcp_R_dT(s, t); }, T, h);
      const double tol = 1e-8 * std::max(std::abs(fd), 1e-6);
      EXPECT_NEAR(ana, fd, tol) << "species=" << species_name(s) << " T=" << T;
    }
    // Beyond the high-tail fade (T_hi * 1.1 <= 22000 K), cp/R is constant.
    EXPECT_DOUBLE_EQ(d2cp_R_dT2(s, 30000.0), 0.0) << "species=" << species_name(s);
  }

  // Mixture d2cp_dT2 on dry air and combustion products.
  for (const auto& X : {dry_air(), hot_products()}) {
    for (double T : {120.0, 300.0, 980.0, 1000.0, 1500.0, 3500.0}) {
      const double h = 1e-4 * T;
      const double ana = d2cp_dT2(T, X);
      const double fd = richardson_fd([&](double t) { return dcp_dT(t, X); }, T, h);
      EXPECT_NEAR(ana, fd, 1e-9 * std::max(std::abs(fd), 1e-6)) << "T=" << T;
    }
  }

  // Out-of-bounds species throws.
  EXPECT_THROW(d2cp_R_dT2(num_species(), 300.0), std::out_of_range);
}

TEST(FannoRayleighTangents1276, UnchokedTangentsMatchRelationFiniteDifferences) {
  const auto air = dry_air();
  const auto prod = hot_products();
  struct Case {
    const char* name;
    double Pt0, Tt0, Tt_exit, P_target;
    bool exit_total;
    double L, D, rough;
    std::string model;
    double fmult;
    bool use_prod;
  };
  const std::vector<Case> cases = {
      {"heated_haaland_static", 5e5, 800.0, 1200.0, 3.0e5, false, 0.30, 0.02, 1e-5, "haaland", 1.0,
       false},
      {"heated_serghides_total", 5e5, 700.0, 1100.0, 4.4e5, true, 0.40, 0.02, 2e-5, "serghides", 1.0,
       true},
      {"cooled_colebrook_static", 5e5, 1100.0, 650.0, 4.6e5, false, 0.50, 0.02, 5e-5, "colebrook",
       1.0, false},
      {"adiabatic_petukhov_total", 4e5, 600.0, 600.0, 3.8e5, true, 0.60, 0.025, 0.0, "petukhov", 1.0,
       false},
      {"heated_fixed_static", 4e5, 500.0, 900.0, 2.5e5, false, 0.40, 0.02, 0.0, "fixed", 0.025,
       false},
      {"c2_bridge_crossing", 3e5, 920.0, 1080.0, 2.2e5, false, 0.25, 0.02, 1e-5, "haaland", 1.0,
       true},
  };

  for (const auto& c : cases) {
    const auto& X = c.use_prod ? prod : air;
    const FannoChannelFlow r = fanno_rayleigh_channel_flow(
        c.Pt0, c.Tt0, c.Tt_exit, X, c.P_target, c.exit_total, c.L, c.D, c.rough, c.model, c.fmult,
        -1.0, {}, true);
    ASSERT_FALSE(r.choked) << c.name;

    // Direct 4th-order Richardson FD of the march exit pressure E(G, Pt0, Tt0, Tt_exit)
    // avoids Brent root-tolerance noise and verifies the exact discrete tangent.
    auto E_of = [&](double Gx, double Ptx, double T0x, double Tex) {
      const FannoRayleighMarch m =
          fanno_rayleigh_march(Gx, Ptx, T0x, Tex, X, c.L, c.D, c.rough, c.model, c.fmult);
      return c.exit_total ? m.Pt_exit : m.P_exit;
    };
    const double E_G =
        richardson_fd([&](double v) { return E_of(v, c.Pt0, c.Tt0, c.Tt_exit); }, r.G, 1e-4 * r.G);
    const double E_Pt = richardson_fd(
        [&](double v) { return E_of(r.G, v, c.Tt0, c.Tt_exit); }, c.Pt0, 1e-4 * c.Pt0);
    const double E_T0 = richardson_fd(
        [&](double v) { return E_of(r.G, c.Pt0, v, c.Tt_exit); }, c.Tt0, 1e-4 * c.Tt0);
    const double E_Te = richardson_fd(
        [&](double v) { return E_of(r.G, c.Pt0, c.Tt0, v); }, c.Tt_exit, 1e-4 * c.Tt_exit);

    EXPECT_NEAR(r.dG_dP_target, 1.0 / E_G, 1e-8 * std::abs(1.0 / E_G)) << c.name;
    EXPECT_NEAR(r.dG_dPt0, -E_Pt / E_G, 1e-8 * std::abs(E_Pt / E_G)) << c.name;
    EXPECT_NEAR(r.dG_dTt0, -E_T0 / E_G, 1e-8 * std::abs(E_T0 / E_G)) << c.name;
    EXPECT_NEAR(r.dG_dTt_exit, -E_Te / E_G, 1e-8 * std::abs(E_Te / E_G)) << c.name;

    // Primal values must be bit-identical whether derivatives are requested or not.
    const FannoChannelFlow r_no = fanno_rayleigh_channel_flow(
        c.Pt0, c.Tt0, c.Tt_exit, X, c.P_target, c.exit_total, c.L, c.D, c.rough, c.model, c.fmult,
        -1.0, {}, false);
    EXPECT_DOUBLE_EQ(r.G, r_no.G) << c.name;
    EXPECT_DOUBLE_EQ(r.M_in, r_no.M_in) << c.name;
    EXPECT_DOUBLE_EQ(r.M_exit, r_no.M_exit) << c.name;
    EXPECT_EQ(r.choked, r_no.choked) << c.name;
  }
}

TEST(FannoRayleighTangents1276, ChokedTangentsMatchFiniteDifferencesAllModes) {
  const auto air = dry_air();
  struct Case {
    const char* name;
    double Pt0, Tt0, Tt_exit, P_target;
    bool exit_total;
    double L, D, rough;
    std::string model;
    double fmult;
    double rel_tol;
  };
  const std::vector<Case> cases = {
      // kInlet: strong cooling keeps sonic inlet subsonic downstream.
      {"kInlet_cooled", 5e5, 1000.0, 300.0, 1.0e5, false, 0.02, 0.05, 1e-6, "haaland", 1.0, 1e-8},
      // kExit: heated & cooled sonic exit.
      {"kExit_heated", 5e5, 600.0, 1100.0, 1.0e5, false, 0.50, 0.02, 1e-5, "serghides", 1.0, 1e-7},
      {"kExit_cooled", 5e5, 800.0, 650.0, 1.0e5, true, 2.00, 0.02, 1e-5, "petukhov", 1.0, 1e-7},
      {"kExit_fixed", 4e5, 500.0, 900.0, 1.0e5, false, 0.40, 0.02, 0.0, "fixed", 0.02, 1e-7},
      // kThroat: cooled duct with interior sonic throat (re-solve branch).
      {"kThroat_haaland", 5e5, 600.0, 320.0, 1.0e4, false, 1.00, 0.02, 1e-5, "haaland", 1.0, 1e-4},
      {"kThroat_fixed", 5e5, 600.0, 380.0, 1.0e4, false, 0.70, 0.02, 1e-5, "fixed", 0.02, 1e-4},
  };

  for (const auto& c : cases) {
    const FannoChannelFlow r = fanno_rayleigh_channel_flow(
        c.Pt0, c.Tt0, c.Tt_exit, air, c.P_target, c.exit_total, c.L, c.D, c.rough, c.model, c.fmult,
        -1.0, {}, true);
    ASSERT_TRUE(r.choked) << c.name;
    EXPECT_DOUBLE_EQ(r.dG_dP_target, 0.0) << c.name;

    auto Gc = [&](double Ptx, double T0x, double Tex) {
      return fanno_rayleigh_choked_mass_flux(Ptx, T0x, Tex, air, c.L, c.D, c.rough, c.model,
                                             c.fmult);
    };
    const double fd_Pt = richardson_fd([&](double v) { return Gc(v, c.Tt0, c.Tt_exit); }, c.Pt0,
                                       1e-4 * c.Pt0);
    const double fd_T0 = richardson_fd([&](double v) { return Gc(c.Pt0, v, c.Tt_exit); }, c.Tt0,
                                       1e-4 * c.Tt0);
    const double fd_Te = richardson_fd([&](double v) { return Gc(c.Pt0, c.Tt0, v); }, c.Tt_exit,
                                       1e-4 * c.Tt_exit);

    EXPECT_NEAR(r.dG_dPt0, fd_Pt, c.rel_tol * std::max(std::abs(fd_Pt), 1e-12)) << c.name;
    EXPECT_NEAR(r.dG_dTt0, fd_T0, c.rel_tol * std::max(std::abs(fd_T0), 1e-12)) << c.name;
    EXPECT_NEAR(r.dG_dTt_exit, fd_Te, c.rel_tol * std::max(std::abs(fd_Te), 1e-12)) << c.name;
  }
}

TEST(FannoRayleighTangents1276, LowDriveAndNearChokeStayNoiseFree) {
  const auto air = dry_air();
  // Low drive (1.5 Pa on 2 bar): eager (cold) and warm-started tangents must
  // match cleanly with no FD cancellation scatter, and agree with resolved
  // march FD of E(G, Pt0, Tt0, Tt_exit).
  const double Pt0 = 2e5;
  const double P_target = Pt0 - 1.5;
  const FannoChannelFlow cold = fanno_rayleigh_channel_flow(
      Pt0, 600.0, 750.0, air, P_target, false, 0.05, 0.01, 1e-5, "haaland", 1.0, -1.0, {}, true);
  const FannoChannelFlow warm = fanno_rayleigh_channel_flow(
      Pt0, 600.0, 750.0, air, P_target, false, 0.05, 0.01, 1e-5, "haaland", 1.0, cold.G * 1.001,
      {}, true);
  ASSERT_FALSE(cold.choked);
  EXPECT_NEAR(cold.dG_dP_target, warm.dG_dP_target, 1e-8 * std::abs(cold.dG_dP_target));
  EXPECT_NEAR(cold.dG_dPt0, warm.dG_dPt0, 1e-8 * std::abs(cold.dG_dPt0));
  EXPECT_NEAR(cold.dG_dTt0, warm.dG_dTt0, 1e-8 * std::abs(cold.dG_dTt0));
  EXPECT_NEAR(cold.dG_dTt_exit, warm.dG_dTt_exit, 1e-8 * std::abs(cold.dG_dTt_exit));

  auto E_ld = [&](double Gx, double Ptx, double T0x, double Tex) {
    return fanno_rayleigh_march(Gx, Ptx, T0x, Tex, air, 0.05, 0.01, 1e-5, "haaland", 1.0).P_exit;
  };
  // At 1.5 Pa on 2e5 Pa, fanno_rayleigh_march.P_exit carries ~2.5e-9 Pa root
  // jitter from inlet_theta; odd_poly5_deriv on measured plateau spans resolves
  // dG/dP, dG/dPt0, dG/dTt_exit to ~1e-8 and dG/dTt0 to ~1e-7.
  const double E_ld_G =
      odd_poly5_deriv([&](double v) { return E_ld(v, Pt0, 600.0, 750.0); }, cold.G, 0.08 * cold.G);
  const double E_ld_Pt =
      odd_poly5_deriv([&](double v) { return E_ld(cold.G, v, 600.0, 750.0); }, Pt0, 0.01 * Pt0);
  const double E_ld_T0 =
      odd_poly5_deriv([&](double v) { return E_ld(cold.G, Pt0, v, 750.0); }, 600.0, 70.0);
  const double E_ld_Te =
      odd_poly5_deriv([&](double v) { return E_ld(cold.G, Pt0, 600.0, v); }, 750.0, 40.0);
  EXPECT_NEAR(cold.dG_dP_target, 1.0 / E_ld_G, 1e-7 * std::abs(1.0 / E_ld_G));
  EXPECT_NEAR(cold.dG_dPt0, -E_ld_Pt / E_ld_G, 1e-7 * std::abs(E_ld_Pt / E_ld_G));
  EXPECT_NEAR(cold.dG_dTt0, -E_ld_T0 / E_ld_G, 1e-6 * std::abs(E_ld_T0 / E_ld_G));
  EXPECT_NEAR(cold.dG_dTt_exit, -E_ld_Te / E_ld_G, 1e-7 * std::abs(E_ld_Te / E_ld_G));

  // Near-choke: check a cleanly resolved point (1 - G/Gc = 1e-4, step 1e-6 on
  // a 1e-9..4e-8 plateau) against Richardson march FD, and an extreme point
  // (1 - G/Gc = 1e-6, where FD sits on the noise floor) for finiteness, sign,
  // and physical sqrt(1 - G/Gc) scaling.
  const double Gc =
      fanno_rayleigh_choked_mass_flux(5e5, 800.0, 1200.0, air, 0.30, 0.02, 1e-5, "haaland", 1.0);
  auto E_nc = [&](double Gx, double Ptx, double T0x, double Tex) {
    return fanno_rayleigh_march(Gx, Ptx, T0x, Tex, air, 0.30, 0.02, 1e-5, "haaland", 1.0).P_exit;
  };
  double dG_dP_1e4 = 0.0;
  double dG_dPt0_1e4 = 0.0;
  {
    const FannoRayleighMarch m_nc = fanno_rayleigh_march(Gc * (1.0 - 1e-4), 5e5, 800.0, 1200.0, air,
                                                         0.30, 0.02, 1e-5, "haaland", 1.0);
    const FannoChannelFlow r_nc = fanno_rayleigh_channel_flow(
        5e5, 800.0, 1200.0, air, m_nc.P_exit, false, 0.30, 0.02, 1e-5, "haaland", 1.0, -1.0, {},
        true);
    ASSERT_FALSE(r_nc.choked);
    dG_dP_1e4 = r_nc.dG_dP_target;
    dG_dPt0_1e4 = r_nc.dG_dPt0;
    const double E_G =
        richardson_fd([&](double v) { return E_nc(v, 5e5, 800.0, 1200.0); }, r_nc.G, 1e-6 * r_nc.G);
    const double E_Pt =
        richardson_fd([&](double v) { return E_nc(r_nc.G, v, 800.0, 1200.0); }, 5e5, 1e-6 * 5e5);
    EXPECT_NEAR(r_nc.dG_dP_target, 1.0 / E_G, 1e-7 * std::abs(1.0 / E_G));
    EXPECT_NEAR(r_nc.dG_dPt0, -E_Pt / E_G, 1e-7 * std::abs(E_Pt / E_G));
  }
  {
    const FannoRayleighMarch m_nc = fanno_rayleigh_march(Gc * (1.0 - 1e-6), 5e5, 800.0, 1200.0, air,
                                                         0.30, 0.02, 1e-5, "haaland", 1.0);
    const FannoChannelFlow r_nc = fanno_rayleigh_channel_flow(
        5e5, 800.0, 1200.0, air, m_nc.P_exit, false, 0.30, 0.02, 1e-5, "haaland", 1.0, -1.0, {},
        true);
    ASSERT_FALSE(r_nc.choked);
    EXPECT_TRUE(std::isfinite(r_nc.dG_dP_target));
    EXPECT_TRUE(std::isfinite(r_nc.dG_dPt0));
    EXPECT_LT(r_nc.dG_dP_target, 0.0);
    EXPECT_GT(r_nc.dG_dPt0, 0.0);
    // Near a quadratic sonic fold, dG/dP scales as sqrt(1 - G/Gc) (~0.1x from
    // 1e-4 to 1e-6) while dG/dPt0 stays close to dGc/dPt0 (~0.99x).
    const double ratio_P = r_nc.dG_dP_target / dG_dP_1e4;
    const double ratio_Pt = r_nc.dG_dPt0 / dG_dPt0_1e4;
    EXPECT_GT(ratio_P, 0.05);
    EXPECT_LT(ratio_P, 0.2);
    EXPECT_GT(ratio_Pt, 0.9);
    EXPECT_LT(ratio_Pt, 1.1);
  }
}

TEST(FannoRayleighTangents1276, CompositionDirectionsMatchFD) {
  const auto X = hot_products();
  const std::vector<double> Y = mole_to_mass(X);
  std::vector<double> dir(Y.size(), 0.0);
  dir[species_index_from_name("N2")] = -0.5;
  dir[species_index_from_name("H2O")] = 0.5;
  const std::vector<std::vector<double>> dirs{dir};

  auto check_case = [&](double Pt0, double Tt0, double Tt_ex, double P_t, double L, double D) {
    const FannoChannelFlow r = fanno_rayleigh_channel_flow(
        Pt0, Tt0, Tt_ex, X, P_t, false, L, D, 1e-5, "haaland", 1.0, -1.0, dirs, true);
    ASSERT_EQ(r.dG_ddir.size(), 1u);
    auto G_at_s = [&](double s) {
      std::vector<double> Ys = Y;
      for (std::size_t k = 0; k < Ys.size(); ++k) {
        Ys[k] += s * dir[k];
      }
      return fanno_rayleigh_channel_flow(Pt0, Tt0, Tt_ex, mass_to_mole(Ys), P_t, false, L, D, 1e-5,
                                         "haaland", 1.0, r.G, {}, false)
          .G;
    };
    const double fd = richardson_fd(G_at_s, 0.0, 1e-4);
    EXPECT_NEAR(r.dG_ddir[0], fd, 1e-6 * std::max(std::abs(fd), 1e-12));
  };

  // Unchoked, kExit, and kInlet.
  check_case(5e5, 800.0, 1200.0, 3.0e5, 0.30, 0.02);
  check_case(5e5, 600.0, 1100.0, 1.0e5, 0.50, 0.02);
  check_case(5e5, 1000.0, 300.0, 1.0e5, 0.02, 0.05);
}

}  // namespace
}  // namespace combaero
