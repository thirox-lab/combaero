#include "area_change.h"
#include "composition.h"
#include "math_constants.h"
#include "solver_interface.h"
#include "stagnation.h"
#include "thermo.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace combaero;
using namespace combaero::solver;

// Global flag to enable detailed difference reporting
bool g_print_jacobian_differences = false;

// File output settings (empty string means no file output configured)
std::string g_jacobian_output_file =
    ""; // Will be set by JACOBIAN_OUTPUT_FILE env var
bool g_write_jacobian_to_file = false;

// Helper to collect and report differences
struct JacobianDifference {
  std::string test_name;
  std::string derivative;
  double analytical;
  double finite_diff;
  double absolute_diff;
  double relative_diff;
};

std::vector<JacobianDifference> g_jacobian_differences;

void report_jacobian_difference(const std::string &test_name,
                                const std::string &derivative,
                                double analytical, double finite_diff) {
  double abs_diff = std::abs(analytical - finite_diff);
  double rel_diff =
      (std::abs(finite_diff) > 1e-15) ? abs_diff / std::abs(finite_diff) : 0.0;

  // Console output
  if (g_print_jacobian_differences) {
    std::cout << std::fixed << std::setprecision(6)
              << "JACOBIAN_DIFF: " << test_name << " | " << derivative
              << " | analytical=" << analytical << " | fd=" << finite_diff
              << " | abs_diff=" << abs_diff << " | rel_diff=" << rel_diff
              << std::endl;
  }

  // File output
  if (g_write_jacobian_to_file && !g_jacobian_output_file.empty()) {
    std::ofstream file(g_jacobian_output_file, std::ios::app);
    if (file.is_open()) {
      file << std::fixed << std::setprecision(12) << test_name << ","
           << derivative << "," << analytical << "," << finite_diff << ","
           << abs_diff << "," << rel_diff << std::endl;
      file.close();
    }
  }

  g_jacobian_differences.push_back(
      {test_name, derivative, analytical, finite_diff, abs_diff, rel_diff});
}

// Print summary of all collected differences
void print_jacobian_difference_summary() {
  if (g_jacobian_differences.empty())
    return;

  std::cout << "\n=== JACOBIAN DIFFERENCE SUMMARY ===" << std::endl;
  std::cout << std::setw(25) << "Test" << std::setw(20) << "Derivative"
            << std::setw(15) << "Abs Diff" << std::setw(15) << "Rel Diff"
            << std::endl;
  std::cout << std::string(75, '-') << std::endl;

  for (const auto &diff : g_jacobian_differences) {
    std::cout << std::fixed << std::setprecision(3) << std::setw(25)
              << diff.test_name << std::setw(20) << diff.derivative
              << std::setw(15) << diff.absolute_diff << std::setw(15)
              << diff.relative_diff << std::endl;
  }

  // Statistics
  double max_abs = 0, max_rel = 0;
  for (const auto &diff : g_jacobian_differences) {
    max_abs = std::max(max_abs, diff.absolute_diff);
    max_rel = std::max(max_rel, diff.relative_diff);
  }

  std::cout << std::string(75, '-') << std::endl;
  std::cout << "Maximum absolute difference: " << std::scientific << max_abs
            << std::endl;
  std::cout << "Maximum relative difference: " << std::scientific << max_rel
            << std::endl;
  std::cout << "Total comparisons: " << g_jacobian_differences.size()
            << std::endl;
  std::cout << "=== END SUMMARY ===" << std::endl;
}

// Check environment variables once at startup
bool check_jacobian_diff_flag() {
  const char *env = std::getenv("PRINT_JACOBIAN_DIFFERENCES");
  return env && (std::string(env) == "1" || std::string(env) == "true");
}

void initialize_jacobian_file_output() {
  const char *file_env = std::getenv("JACOBIAN_OUTPUT_FILE");
  if (file_env) {
    g_jacobian_output_file = std::string(file_env);
    g_write_jacobian_to_file = true;

    // Write CSV header
    std::ofstream file(g_jacobian_output_file, std::ios::app);
    if (file.is_open()) {
      // Check if file is empty, write header if needed
      file.seekp(0, std::ios::end);
      if (file.tellp() == 0) {
        file << "test_name,derivative,analytical,finite_diff,abs_diff,rel_diff"
             << std::endl;
      }
      file.close();
    }
  }
}

// Initialize flags
static bool jacobian_diff_flag_initialized = []() {
  g_print_jacobian_differences = check_jacobian_diff_flag();
  initialize_jacobian_file_output();
  return true;
}();

TEST(SolverJacobianTest, OrificeDerivatives) {
  // Boilerplate for species
  std::size_t ns = num_species();
  std::vector<double> Y(ns, 0.0);
  // Assume air-like for testing
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      Y[i] = 0.77;
    else if (species_name(i) == "O2")
      Y[i] = 0.23;
  }

  double m_dot = 0.1;
  double P_total_up = 200000.0;
  double P_static_up = 190000.0;
  double T_up = 300.0;
  double P_static_down = 180000.0;
  double Cd = 0.62;
  double area = 0.001;

  double beta = 0.5;
  OrificeResult res = orifice_residuals_and_jacobian(
      m_dot, P_total_up, P_static_up, T_up, Y, P_static_down, Cd, area, beta);

  // FD checks
  double eps = 1.0;

  // d_mdot_dP_total_up
  auto res_p =
      orifice_residuals_and_jacobian(m_dot, P_total_up + eps, P_static_up, T_up,
                                     Y, P_static_down, Cd, area, beta);
  double fd_dP_total_up = (res_p.m_dot_calc - res.m_dot_calc) / eps;
  report_jacobian_difference("OrificeDerivatives", "d_mdot_dP_total_up",
                             res.d_mdot_dP_total_up, fd_dP_total_up);
  EXPECT_NEAR(res.d_mdot_dP_total_up, fd_dP_total_up, 1e-10);

  // d_mdot_dP_static_down
  auto res_down_p =
      orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up, Y,
                                     P_static_down + eps, Cd, area, beta);
  auto res_down_m =
      orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up, Y,
                                     P_static_down - eps, Cd, area, beta);
  double fd_dP_static_down = (res_down_p.m_dot_calc - res_down_m.m_dot_calc) / (2.0 * eps);
  report_jacobian_difference("OrificeDerivatives", "d_mdot_dP_static_down",
                             res.d_mdot_dP_static_down, fd_dP_static_down);
  EXPECT_NEAR(res.d_mdot_dP_static_down, fd_dP_static_down, 1e-10);

  // d_mdot_dP_static_up: with beta > 0 the drive (approach_from_static) plus
  // the density, so the stagnation column above is exactly zero.
  auto res_up_p =
      orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up + eps, T_up,
                                     Y, P_static_down, Cd, area, beta);
  auto res_up_m =
      orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up - eps, T_up,
                                     Y, P_static_down, Cd, area, beta);
  double fd_dP_static_up = (res_up_p.m_dot_calc - res_up_m.m_dot_calc) / (2.0 * eps);
  report_jacobian_difference("OrificeDerivatives", "d_mdot_dP_static_up",
                             res.d_mdot_dP_static_up, fd_dP_static_up);
  EXPECT_NEAR(res.d_mdot_dP_static_up, fd_dP_static_up,
              std::abs(fd_dP_static_up) * 1e-7 + 1e-12);
  EXPECT_EQ(res.d_mdot_dP_total_up, 0.0);

  // A named density reference takes the density column with it.
  const double P_rho = 185000.0;
  const OrificeResult rr = orifice_residuals_and_jacobian(
      m_dot, P_total_up, P_static_up, T_up, Y, P_static_down, Cd, area, beta, P_rho);
  const double fd_dP_rho =
      (orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up, Y,
                                      P_static_down, Cd, area, beta, P_rho + eps)
           .m_dot_calc -
       orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up, Y,
                                      P_static_down, Cd, area, beta, P_rho - eps)
           .m_dot_calc) /
      (2.0 * eps);
  EXPECT_NEAR(rr.d_mdot_dP_rho_ref, fd_dP_rho, std::abs(fd_dP_rho) * 1e-7 + 1e-14);
  EXPECT_GT(std::abs(rr.d_mdot_dP_rho_ref), 0.0);

  // d_mdot_dTt_up
  double eps_T = 1e-4;
  auto res_T_p = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up,
                                                T_up + eps_T, Y, P_static_down,
                                                Cd, area, beta);
  auto res_T_m = orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up,
                                                T_up - eps_T, Y, P_static_down,
                                                Cd, area, beta);
  double fd_dT_up = (res_T_p.m_dot_calc - res_T_m.m_dot_calc) / (2.0 * eps_T);
  report_jacobian_difference("OrificeDerivatives", "d_mdot_dTt_up",
                             res.d_mdot_dTt_up, fd_dT_up);
  EXPECT_NEAR(res.d_mdot_dTt_up, fd_dT_up, std::abs(fd_dT_up) * 1e-6);

  // Y derivatives
  double eps_Y = 1e-6;
  for (size_t i = 0; i < ns; ++i) {
    std::vector<double> Y_p = Y;
    Y_p[i] += eps_Y;
    auto res_Y_p =
        orifice_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up,
                                       Y_p, P_static_down, Cd, area, beta);
    double fd_dYi = (res_Y_p.m_dot_calc - res.m_dot_calc) / eps_Y;
    EXPECT_NEAR(res.d_mdot_dY_up[i], fd_dYi, std::abs(fd_dYi) * 1e-3 + 1e-8);
  }
}

TEST(SolverJacobianTest, PipeDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> Y(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      Y[i] = 0.77;
    else if (species_name(i) == "O2")
      Y[i] = 0.23;
  }

  double m_dot = 0.5;
  double P_total_up = 200000.0;
  double P_static_up = 190000.0;
  double T_up = 300.0;
  double P_static_down = 180000.0;
  double L = 1.0;
  double D = 0.05;
  double roughness = 1e-5;
  std::string model = "haaland";

  PipeResult res =
      channel_residuals_and_jacobian(m_dot, P_total_up, P_static_up, T_up, Y,
                                     P_static_down, L, D, roughness, model);

  double eps_m = 1e-4;
  auto res_m_p = channel_residuals_and_jacobian(
      m_dot + eps_m, P_total_up, P_static_up, T_up, Y, P_static_down, L, D,
      roughness, model);
  auto res_m_m = channel_residuals_and_jacobian(
      m_dot - eps_m, P_total_up, P_static_up, T_up, Y, P_static_down, L, D,
      roughness, model);
  double fd_dmdot = (res_m_p.dP_calc - res_m_m.dP_calc) / (2.0 * eps_m);
  report_jacobian_difference("PipeDerivatives", "d_dP_d_mdot", res.d_dP_d_mdot,
                             fd_dmdot);
  EXPECT_NEAR(res.d_dP_d_mdot, fd_dmdot, std::abs(fd_dmdot) * 1e-8);

  double eps_P = 1.0;
  auto res_P_p = channel_residuals_and_jacobian(
      m_dot, P_total_up, P_static_up + eps_P, T_up, Y, P_static_down, L, D,
      roughness, model);
  auto res_P_m = channel_residuals_and_jacobian(
      m_dot, P_total_up, P_static_up - eps_P, T_up, Y, P_static_down, L, D,
      roughness, model);
  double fd_dP_static = (res_P_p.dP_calc - res_P_m.dP_calc) / (2.0 * eps_P);
  report_jacobian_difference("PipeDerivatives", "d_dP_dP_static_up",
                             res.d_dP_dP_static_up, fd_dP_static);
  EXPECT_NEAR(res.d_dP_dP_static_up, fd_dP_static,
              std::abs(fd_dP_static) * 1e-8);

  double eps_T = 0.1;
  auto res_T_p = channel_residuals_and_jacobian(m_dot, P_total_up, P_static_up,
                                                T_up + eps_T, Y, P_static_down,
                                                L, D, roughness, model);
  auto res_T_m = channel_residuals_and_jacobian(m_dot, P_total_up, P_static_up,
                                                T_up - eps_T, Y, P_static_down,
                                                L, D, roughness, model);
  double fd_dT = (res_T_p.dP_calc - res_T_m.dP_calc) / (2.0 * eps_T);
  report_jacobian_difference("PipeDerivatives", "d_dP_dTt_up", res.d_dP_dTt_up,
                             fd_dT);
  EXPECT_NEAR(res.d_dP_dTt_up, fd_dT, std::abs(fd_dT) * 1e-8);

  double eps_Y = 1e-6;
  for (size_t i = 0; i < ns; ++i) {
    std::vector<double> Y_p = Y;
    Y_p[i] += eps_Y;
    auto res_Y_p = channel_residuals_and_jacobian(
        m_dot, P_total_up, P_static_up, T_up, Y_p, P_static_down, L, D,
        roughness, model);
    double fd_dYi = (res_Y_p.dP_calc - res.dP_calc) / eps_Y;
    EXPECT_NEAR(res.d_dP_dY_up[i], fd_dYi, std::abs(fd_dYi) * 1e-2 + 1e-8);
  }
}

// exit_head_lost (#588): the drop gains m|m| / (2 rho A^2), rho the upstream
// static density, and every derivative carries the term. Both flow directions:
// the head is odd in m, its slope even.
TEST(SolverJacobianTest, PipeDerivativesExitHeadLost) {
  std::size_t ns = num_species();
  std::vector<double> Y(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      Y[i] = 0.77;
    else if (species_name(i) == "O2")
      Y[i] = 0.23;
  }
  const double P_total_up = 200000.0;
  const double P_static_up = 190000.0;
  const double T_up = 300.0;
  const double P_static_down = 180000.0;
  const double L = 1.0;
  const double D = 0.05;
  const double roughness = 1e-5;
  const std::string model = "haaland";
  auto eval = [&](double m, double P, double T, const std::vector<double> &Yv,
                  bool lost) {
    return channel_residuals_and_jacobian(m, P_total_up, P, T, Yv,
                                          P_static_down, L, D, roughness,
                                          model, 1.0, lost);
  };

  for (double m_dot : {0.5, -0.5}) {
    PipeResult base = eval(m_dot, P_static_up, T_up, Y, false);
    PipeResult res = eval(m_dot, P_static_up, T_up, Y, true);

    const std::vector<double> X =
        combaero::mass_to_mole(combaero::normalize_fractions(Y));
    const double rho = combaero::density(T_up, P_static_up, X);
    const double area = 0.25 * M_PI * D * D;
    const double q = 0.5 * m_dot * std::abs(m_dot) / (rho * area * area);
    EXPECT_NEAR(res.dP_calc - base.dP_calc, q, std::abs(q) * 1e-12);

    const double eps_m = 1e-4;
    const double fd_m = (eval(m_dot + eps_m, P_static_up, T_up, Y, true).dP_calc -
                         eval(m_dot - eps_m, P_static_up, T_up, Y, true).dP_calc) /
                        (2.0 * eps_m);
    EXPECT_NEAR(res.d_dP_d_mdot, fd_m, std::abs(fd_m) * 1e-7);

    const double eps_P = 1.0;
    const double fd_P = (eval(m_dot, P_static_up + eps_P, T_up, Y, true).dP_calc -
                         eval(m_dot, P_static_up - eps_P, T_up, Y, true).dP_calc) /
                        (2.0 * eps_P);
    EXPECT_NEAR(res.d_dP_dP_static_up, fd_P, std::abs(fd_P) * 1e-6);

    const double eps_T = 0.1;
    const double fd_T = (eval(m_dot, P_static_up, T_up + eps_T, Y, true).dP_calc -
                         eval(m_dot, P_static_up, T_up - eps_T, Y, true).dP_calc) /
                        (2.0 * eps_T);
    EXPECT_NEAR(res.d_dP_dTt_up, fd_T, std::abs(fd_T) * 1e-6);

    const double eps_Y = 1e-6;
    for (size_t i = 0; i < ns; ++i) {
      std::vector<double> Y_p = Y;
      Y_p[i] += eps_Y;
      const double fd_Y =
          (eval(m_dot, P_static_up, T_up, Y_p, true).dP_calc - res.dP_calc) /
          eps_Y;
      EXPECT_NEAR(res.d_dP_dY_up[i], fd_Y, std::abs(fd_Y) * 1e-2 + 1e-8);
    }
  }
}

TEST(SolverJacobianTest, OrificeCompressibleDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> Y(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      Y[i] = 0.77;
    else if (species_name(i) == "O2")
      Y[i] = 0.23;
  }

  const double m_dot = 0.30;
  const double P_total_up = 200000.0;
  const double P_static_up = 190000.0;  // a flowing node: Pt > P
  const double T_up = 300.0;
  const double P_static_down = 150000.0;
  const double Cd = 0.65;
  const double area = 1e-4;

  // beta = 0 drives from the stagnation pressure, beta = 0.5 from the static
  // one (approach_from_static); each column must match its difference, and
  // the column the rule does not use must be exactly zero.
  for (const double beta : {0.0, 0.5}) {
    auto f = [&](double pt, double ps, double t, double pb) {
      return orifice_compressible_residuals_and_jacobian(m_dot, pt, ps, t, Y, pb, Cd,
                                                         area, beta);
    };
    const OrificeResult res = f(P_total_up, P_static_up, T_up, P_static_down);

    const double eps_P = 1.0;
    const double fd_dPt = (f(P_total_up + eps_P, P_static_up, T_up, P_static_down).m_dot_calc -
                           f(P_total_up - eps_P, P_static_up, T_up, P_static_down).m_dot_calc) /
                          (2.0 * eps_P);
    const double fd_dPs = (f(P_total_up, P_static_up + eps_P, T_up, P_static_down).m_dot_calc -
                           f(P_total_up, P_static_up - eps_P, T_up, P_static_down).m_dot_calc) /
                          (2.0 * eps_P);
    report_jacobian_difference("OrificeCompressibleDerivatives", "d_mdot_dP_total_up",
                               res.d_mdot_dP_total_up, fd_dPt);
    EXPECT_NEAR(res.d_mdot_dP_total_up, fd_dPt, std::abs(fd_dPt) * 1e-8 + 1e-12);
    EXPECT_NEAR(res.d_mdot_dP_static_up, fd_dPs, std::abs(fd_dPs) * 1e-8 + 1e-12);
    if (beta > 0.0) {
      // The drive is the static pressure; Pt only sets the static
      // temperature on the isentrope (#645): a small, nonzero column (a
      // higher Pt at fixed P means a colder, denser jet).
      EXPECT_GT(res.d_mdot_dP_total_up, 0.0);
      EXPECT_LT(res.d_mdot_dP_total_up, 0.1 * std::abs(res.d_mdot_dP_static_up));
      // Behind a node with no velocity (P = Pt) the static state is the
      // stagnation state: T1 = T_up to the bit.
      const OrificeResult plenum = f(P_static_up, P_static_up, T_up, P_static_down);
      const auto X = combaero::mass_to_mole(combaero::normalize_fractions(Y));
      const double m_nozzle = std::get<0>(orifice_compressible_mdot_and_jacobian(
          T_up, P_static_up, P_static_down, X, Cd, area, beta));
      EXPECT_EQ(plenum.m_dot_calc, m_nozzle);
      // Its Pt column is not zero: the isentrope still moves T1 with Pt.
      EXPECT_GT(plenum.d_mdot_dP_total_up, 0.0);
      const double fd_pl_Pt = (f(P_static_up + eps_P, P_static_up, T_up, P_static_down).m_dot_calc -
                               f(P_static_up - eps_P, P_static_up, T_up, P_static_down).m_dot_calc) /
                              (2.0 * eps_P);
      const double fd_pl_Ps = (f(P_static_up, P_static_up + eps_P, T_up, P_static_down).m_dot_calc -
                               f(P_static_up, P_static_up - eps_P, T_up, P_static_down).m_dot_calc) /
                              (2.0 * eps_P);
      EXPECT_NEAR(plenum.d_mdot_dP_total_up, fd_pl_Pt, std::abs(fd_pl_Pt) * 1e-6 + 1e-12);
      EXPECT_NEAR(plenum.d_mdot_dP_static_up, fd_pl_Ps, std::abs(fd_pl_Ps) * 1e-6 + 1e-12);
      // With velocity the nozzle starts from the static T: the flow is the
      // nozzle's at T1 < T_up, so higher than at T_up.
      EXPECT_GT(res.m_dot_calc, std::get<0>(orifice_compressible_mdot_and_jacobian(
                                    T_up, P_static_up, P_static_down, X, Cd, area, beta)));
    } else {
      EXPECT_EQ(res.d_mdot_dP_static_up, 0.0);
    }

    const double fd_dPb = (f(P_total_up, P_static_up, T_up, P_static_down + eps_P).m_dot_calc -
                           f(P_total_up, P_static_up, T_up, P_static_down - eps_P).m_dot_calc) /
                          (2.0 * eps_P);
    report_jacobian_difference("OrificeCompressibleDerivatives", "d_mdot_dP_static_down",
                               res.d_mdot_dP_static_down, fd_dPb);
    EXPECT_NEAR(res.d_mdot_dP_static_down, fd_dPb, std::abs(fd_dPb) * 1e-7 + 1e-12);

    const double eps_T = 1e-3;
    const double fd_dT = (f(P_total_up, P_static_up, T_up + eps_T, P_static_down).m_dot_calc -
                          f(P_total_up, P_static_up, T_up - eps_T, P_static_down).m_dot_calc) /
                         (2.0 * eps_T);
    report_jacobian_difference("OrificeCompressibleDerivatives", "d_mdot_dTt_up",
                               res.d_mdot_dTt_up, fd_dT);
    EXPECT_NEAR(res.d_mdot_dTt_up, fd_dT, std::abs(fd_dT) * 1e-8 + 1e-12);
  }
}

// The incompressible closure (#588): T_static = T0, Pt = P + G^2 / (2 rho),
// every derivative against a central difference, both flow directions; and
// the low-Mach limit of the compressible closure, gap ~ M^2.
TEST(SolverJacobianTest, MomentumChamberIncompressibleDerivativesAndLimit) {
  std::size_t ns = num_species();
  std::vector<double> Y(ns, 0.0);
  Y[species_index_from_name("N2")] = 0.767;
  Y[species_index_from_name("O2")] = 0.233;
  const double P = 101325.0, Pt = 105000.0, T = 600.0, area = 0.01;
  auto f = [&](double P_, double Pt_, double m_, double T_, const std::vector<double> &Y_) {
    return momentum_chamber_residual_and_jacobian(P_, Pt_, m_, T_, Y_, area, true);
  };
  for (double m : {0.5, -0.5}) {
    const auto r = f(P, Pt, m, T, Y);
    EXPECT_DOUBLE_EQ(r.T_static, T);
    const double fdP = (f(P + 1, Pt, m, T, Y).residual - f(P - 1, Pt, m, T, Y).residual) / 2;
    const double fdPt = (f(P, Pt + 1, m, T, Y).residual - f(P, Pt - 1, m, T, Y).residual) / 2;
    const double hm = 1e-5;
    const double fdm =
        (f(P, Pt, m + hm, T, Y).residual - f(P, Pt, m - hm, T, Y).residual) / (2 * hm);
    const double fdT =
        (f(P, Pt, m, T + 1e-3, Y).residual - f(P, Pt, m, T - 1e-3, Y).residual) / 2e-3;
    EXPECT_NEAR(r.d_res_dP, fdP, 1e-8 * std::abs(fdP));
    EXPECT_NEAR(r.d_res_dP_total, fdPt, 1e-12);
    EXPECT_NEAR(r.d_res_dmdot, fdm, 1e-7 * std::abs(fdm));
    EXPECT_NEAR(r.d_res_dT, fdT, 1e-7 * std::abs(fdT));
    const double fdrP =
        (f(P + 1, Pt, m, T, Y).rho_static - f(P - 1, Pt, m, T, Y).rho_static) / 2;
    EXPECT_NEAR(r.drho_dP, fdrP, 1e-8 * std::abs(fdrP));
    for (std::size_t k = 0; k < ns; ++k) {
      std::vector<double> Yp = Y, Ym = Y;
      double fdY = 0.0;
      if (!(Y[k] > 0.0)) {
        // An absent species keeps its partial, one-sided as it enters
        // (#1245): a stream arriving at zero flow brings it.
        Yp[k] += 1e-7;
        fdY = (f(P, Pt, m, T, Yp).residual - r.residual) / 1e-7;
      } else {
        Yp[k] += 1e-5;
        Ym[k] -= 1e-5;
        fdY = (f(P, Pt, m, T, Yp).residual - f(P, Pt, m, T, Ym).residual) / 2e-5;
      }
      EXPECT_NEAR(r.d_res_dY[k], fdY, 1e-4 * std::abs(fdY) + 1e-9) << k;
    }
  }
  // Low-Mach limit of the compressible closure: the dynamic heads agree to
  // O(M^2), the gap quartering as M halves.
  const std::vector<double> X = combaero::mass_to_mole(combaero::normalize_fractions(Y));
  const double rho = combaero::density(T, P, X);
  const double a = combaero::speed_of_sound(T, X);
  double prev = 0.0;
  for (double M : {0.1, 0.05, 0.025}) {
    const double m = M * a * rho * area;
    const auto ri = momentum_chamber_residual_and_jacobian(P, P, m, T, Y, area, true);
    const auto rc = momentum_chamber_residual_and_jacobian(P, P, m, T, Y, area);
    const double rel = (ri.q_dyn - rc.q_dyn) / rc.q_dyn;
    if (prev != 0.0) {
      EXPECT_NEAR(prev / rel, 4.0, 0.1) << M;
    }
    prev = rel;
  }
  EXPECT_LT(std::abs(prev), 1e-4);
}

TEST(SolverJacobianTest, MomentumChamberDerivatives) {
  // Test momentum chamber residual: P_total = P + 0.5 * rho * v^2
  std::size_t ns = num_species();
  std::vector<double> Y(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    Y[i] = (i == species_index_from_name("N2"))   ? 0.767
           : (i == species_index_from_name("O2")) ? 0.233
                                                  : 0.0;
  }

  double P = 101325.0;
  double P_total = 105000.0;
  double m_dot = 0.5;
  double T = 600.0;
  double area = 0.01; // 0.01 m^2

  MomentumChamberResult res =
      momentum_chamber_residual_and_jacobian(P, P_total, m_dot, T, Y, area);

  // Finite difference checks
  double eps_P = 1.0;
  auto res_P_p = momentum_chamber_residual_and_jacobian(P + eps_P, P_total,
                                                        m_dot, T, Y, area);
  double fd_dP = (res_P_p.residual - res.residual) / eps_P;
  report_jacobian_difference("MomentumChamberDerivatives", "d_res_dP",
                             res.d_res_dP, fd_dP);
  EXPECT_NEAR(res.d_res_dP, fd_dP, std::abs(fd_dP) * 1e-6 + 1e-10);

  auto res_Pt_p = momentum_chamber_residual_and_jacobian(P, P_total + eps_P,
                                                         m_dot, T, Y, area);
  double fd_dP_total = (res_Pt_p.residual - res.residual) / eps_P;
  report_jacobian_difference("MomentumChamberDerivatives", "d_res_dP_total",
                             res.d_res_dP_total, fd_dP_total);
  EXPECT_NEAR(res.d_res_dP_total, fd_dP_total,
              std::abs(fd_dP_total) * 1e-8 + 1e-12);

  double eps_mdot = 1e-4;
  auto res_mdot_p = momentum_chamber_residual_and_jacobian(
      P, P_total, m_dot + eps_mdot, T, Y, area);
  auto res_mdot_m = momentum_chamber_residual_and_jacobian(
      P, P_total, m_dot - eps_mdot, T, Y, area);
  double fd_dmdot =
      (res_mdot_p.residual - res_mdot_m.residual) / (2.0 * eps_mdot);
  report_jacobian_difference("MomentumChamberDerivatives", "d_res_dmdot",
                             res.d_res_dmdot, fd_dmdot);
  EXPECT_NEAR(res.d_res_dmdot, fd_dmdot, std::abs(fd_dmdot) * 1e-8 + 1e-12);

  // Verify the residual is the COMPRESSIBLE closure: P_total - P0(P, Ts, M),
  // the isentropic stagnation state from entropy conservation with variable
  // cp. It used to be P_total - P - 0.5*rho*v^2, which under-predicts the
  // stagnation rise -- 0.2% at M = 0.33 but 11% at M = 1.0 and 34% at M = 1.4
  // -- and does so without bound past sonic, so a station driven transonic got
  // a closure that was both wrong and unable to say so. See issue #357.
  //
  // T is the chamber's STAGNATION temperature (#576): the static Ts solves
  // h(T) = h(Ts) + v^2 / 2 with v = m_dot / (rho(Ts, P) A).
  std::vector<double> X = mass_to_mole(normalize_fractions(Y));
  const double Ts = res.T_static;
  EXPECT_LT(Ts, T);
  double rho = density(Ts, P, X);
  double v = m_dot / (rho * area);
  EXPECT_NEAR(h_mass(T, X) - h_mass(Ts, X), 0.5 * v * v, 1e-8 * 0.5 * v * v);
  double M = std::abs(v) / speed_of_sound(Ts, X);
  EXPECT_NEAR(res.M, M, 1e-10);
  EXPECT_NEAR(res.rho_static, rho, 1e-12 * rho);
  double expected_residual = P_total - Pt_from_static(P, Ts, M, X);
  EXPECT_NEAR(res.residual, expected_residual, 1e-6);

  // The incompressible form is the limit it reduces to, not an equal: the
  // stagnation rise must exceed 0.5*rho*v^2 at any finite Mach.
  double q_incompressible = 0.5 * rho * v * v;
  EXPECT_GT(Pt_from_static(P, Ts, M, X) - P, q_incompressible);

  // Static-density sensitivities against central differences.
  const double hP = P * 1e-6;
  const double fd_rho_P =
      (momentum_chamber_residual_and_jacobian(P + hP, P_total, m_dot, T, Y, area).rho_static -
       momentum_chamber_residual_and_jacobian(P - hP, P_total, m_dot, T, Y, area).rho_static) /
      (2.0 * hP);
  EXPECT_NEAR(res.drho_dP, fd_rho_P, std::abs(fd_rho_P) * 1e-6);
  const double hT = T * 1e-6;
  const double fd_rho_T =
      (momentum_chamber_residual_and_jacobian(P, P_total, m_dot, T + hT, Y, area).rho_static -
       momentum_chamber_residual_and_jacobian(P, P_total, m_dot, T - hT, Y, area).rho_static) /
      (2.0 * hT);
  EXPECT_NEAR(res.drho_dT, fd_rho_T, std::abs(fd_rho_T) * 1e-6);
  const double fd_rho_m = (res_mdot_p.rho_static - res_mdot_m.rho_static) / (2.0 * eps_mdot);
  EXPECT_NEAR(res.drho_dmdot, fd_rho_m, std::abs(fd_rho_m) * 1e-5 + 1e-12);
}

TEST(SolverJacobianTest, AdiabaticTCompleteDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> X(ns, 0.0);

  // Fuel-rich mixture (CH4 + air, phi > 1)
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "CH4")
      X[i] = 0.12;
    else if (species_name(i) == "O2")
      X[i] = 0.18;
    else if (species_name(i) == "N2")
      X[i] = 0.70;
  }

  double T_in = 300.0;
  double P = 101325.0;

  auto [T_ad, dT_ad_dT_in, X_products] =
      adiabatic_T_complete_and_jacobian_T(T_in, P, X);

  // Central FD for dT_ad/dT_in
  double eps_T = 1e-3;
  auto [T_ad_plus, dummy1, dummy2] =
      adiabatic_T_complete_and_jacobian_T(T_in + eps_T, P, X);
  auto [T_ad_minus, dummy3, dummy4] =
      adiabatic_T_complete_and_jacobian_T(T_in - eps_T, P, X);
  double fd_dT_ad_dT_in = (T_ad_plus - T_ad_minus) / (2.0 * eps_T);

  report_jacobian_difference("AdiabaticTCompleteDerivatives", "dT_ad_dT_in",
                             dT_ad_dT_in, fd_dT_ad_dT_in);
  EXPECT_NEAR(dT_ad_dT_in, fd_dT_ad_dT_in,
              std::abs(fd_dT_ad_dT_in) * 1e-8 + 1e-12);
  EXPECT_GT(T_ad, T_in); // Sanity: adiabatic T should be elevated
}

TEST(SolverJacobianTest, AdiabaticTEquilibriumDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> X(ns, 0.0);

  // Stoichiometric mixture (CH4 + air)
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "CH4")
      X[i] = 0.095;
    else if (species_name(i) == "O2")
      X[i] = 0.19;
    else if (species_name(i) == "N2")
      X[i] = 0.715;
  }

  double T_in = 300.0;
  double P = 101325.0;

  auto [T_ad, dT_ad_dT_in, dT_ad_dP, X_products] =
      adiabatic_T_equilibrium_and_jacobians(T_in, P, X);

  // Central FD for dT_ad/dT_in
  double eps_T = 1e-3;
  auto [T_ad_T_plus, dummy1, dummy2, dummy3] =
      adiabatic_T_equilibrium_and_jacobians(T_in + eps_T, P, X);
  auto [T_ad_T_minus, dummy4, dummy5, dummy6] =
      adiabatic_T_equilibrium_and_jacobians(T_in - eps_T, P, X);
  double fd_dT_ad_dT_in = (T_ad_T_plus - T_ad_T_minus) / (2.0 * eps_T);

  report_jacobian_difference("AdiabaticTEquilibriumDerivatives", "dT_ad_dT_in",
                             dT_ad_dT_in, fd_dT_ad_dT_in);
  EXPECT_NEAR(dT_ad_dT_in, fd_dT_ad_dT_in,
              std::abs(fd_dT_ad_dT_in) * 1e-8 + 1e-12);

  // Central FD for dT_ad/dP
  double eps_P = 10.0;
  auto [T_ad_P_plus, dummy7, dummy8, dummy9] =
      adiabatic_T_equilibrium_and_jacobians(T_in, P + eps_P, X);
  auto [T_ad_P_minus, dummy10, dummy11, dummy12] =
      adiabatic_T_equilibrium_and_jacobians(T_in, P - eps_P, X);
  double fd_dT_ad_dP = (T_ad_P_plus - T_ad_P_minus) / (2.0 * eps_P);

  report_jacobian_difference("AdiabaticTEquilibriumDerivatives", "dT_ad_dP",
                             dT_ad_dP, fd_dT_ad_dP);
  EXPECT_NEAR(dT_ad_dP, fd_dT_ad_dP, std::abs(fd_dT_ad_dP) * 1e-8 + 1e-12);
  EXPECT_GT(T_ad, T_in); // Sanity: adiabatic T should be elevated
}

TEST(SolverJacobianTest, T0FromStaticDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> X(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      X[i] = 0.79;
    else if (species_name(i) == "O2")
      X[i] = 0.21;
  }

  double T = 300.0;
  double M = 0.5;

  auto [T0, dT0_dM] = Tt_from_static_and_jacobian_M(T, M, X);

  // Central FD for dT0/dM
  double eps_M = 1e-5;
  auto [T0_plus, dummy1] = Tt_from_static_and_jacobian_M(T, M + eps_M, X);
  auto [T0_minus, dummy2] = Tt_from_static_and_jacobian_M(T, M - eps_M, X);
  double fd_dT0_dM = (T0_plus - T0_minus) / (2.0 * eps_M);

  report_jacobian_difference("T0FromStaticDerivatives", "dT0_dM", dT0_dM,
                             fd_dT0_dM);
  EXPECT_NEAR(dT0_dM, fd_dT0_dM, std::abs(fd_dT0_dM) * 1e-6 + 1e-10);
  EXPECT_GT(T0, T); // Sanity: T0 > T for M > 0

  // Test at higher Mach
  M = 0.8;
  auto [T0_high, dT0_dM_high] = Tt_from_static_and_jacobian_M(T, M, X);
  auto [T0_high_plus, dummy3] = Tt_from_static_and_jacobian_M(T, M + eps_M, X);
  auto [T0_high_minus, dummy4] = Tt_from_static_and_jacobian_M(T, M - eps_M, X);
  double fd_dT0_dM_high = (T0_high_plus - T0_high_minus) / (2.0 * eps_M);

  report_jacobian_difference("T0FromStaticDerivatives", "dT0_dM_high",
                             dT0_dM_high, fd_dT0_dM_high);

  EXPECT_NEAR(dT0_dM_high, fd_dT0_dM_high,
              std::abs(fd_dT0_dM_high) * 1e-6 + 1e-10);
}

TEST(SolverJacobianTest, P0FromStaticDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> X(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      X[i] = 0.79;
    else if (species_name(i) == "O2")
      X[i] = 0.21;
  }

  double P = 101325.0;
  double T = 300.0;
  double M = 0.5;

  auto [P0, dP0_dM] = Pt_from_static_and_jacobian_M(P, T, M, X);

  // Central FD for dP0/dM
  double eps_M = 1e-5;
  auto [P0_plus, dummy1] = Pt_from_static_and_jacobian_M(P, T, M + eps_M, X);
  auto [P0_minus, dummy2] = Pt_from_static_and_jacobian_M(P, T, M - eps_M, X);
  double fd_dP0_dM = (P0_plus - P0_minus) / (2.0 * eps_M);

  report_jacobian_difference("P0FromStaticDerivatives", "dP0_dM", dP0_dM,
                             fd_dP0_dM);
  EXPECT_NEAR(dP0_dM, fd_dP0_dM, std::abs(fd_dP0_dM) * 1e-6 + 1e-8);
  EXPECT_GT(P0, P); // Sanity: P0 > P for M > 0

  // Test at higher Mach
  M = 0.8;
  auto [P0_high, dP0_dM_high] = Pt_from_static_and_jacobian_M(P, T, M, X);
  auto [P0_high_plus, dummy3] =
      Pt_from_static_and_jacobian_M(P, T, M + eps_M, X);
  auto [P0_high_minus, dummy4] =
      Pt_from_static_and_jacobian_M(P, T, M - eps_M, X);
  double fd_dP0_dM_high = (P0_high_plus - P0_high_minus) / (2.0 * eps_M);

  report_jacobian_difference("P0FromStaticDerivatives", "dP0_dM_high",
                             dP0_dM_high, fd_dP0_dM_high);
  EXPECT_NEAR(dP0_dM_high, fd_dP0_dM_high,
              std::abs(fd_dP0_dM_high) * 1e-6 + 1e-8);
}

TEST(SolverJacobianTest, OrificeCompressibleMdotDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> X(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      X[i] = 0.79;
    else if (species_name(i) == "O2")
      X[i] = 0.21;
  }

  double T0 = 350.0;
  double P0 = 200000.0;
  double P_back = 150000.0;
  double Cd = 0.65;
  double area = 1e-4;
  double beta = 0.5;

  auto [mdot, dmdot_dP0, dmdot_dP_back, dmdot_dT0] =
      orifice_compressible_mdot_and_jacobian(T0, P0, P_back, X, Cd, area, beta);

  // Central FD for dmdot/dT0
  double eps_T = 1e-3;
  auto [mdot_T_plus, dummy1, dummy2, dummy3] =
      orifice_compressible_mdot_and_jacobian(T0 + eps_T, P0, P_back, X, Cd,
                                             area, beta);
  auto [mdot_T_minus, dummy4, dummy5, dummy6] =
      orifice_compressible_mdot_and_jacobian(T0 - eps_T, P0, P_back, X, Cd,
                                             area, beta);
  double fd_dT0 = (mdot_T_plus - mdot_T_minus) / (2.0 * eps_T);
  report_jacobian_difference("OrificeCompressibleMdotDerivatives", "dmdot_dT0",
                             dmdot_dT0, fd_dT0);
  EXPECT_NEAR(dmdot_dT0, fd_dT0, std::abs(fd_dT0) * 1e-6 + 1e-12);

  // Central FD for dmdot/dP0
  double eps_P = 1.0;
  auto [mdot_P_plus, dummy7, dummy8, dummy9] =
      orifice_compressible_mdot_and_jacobian(T0, P0 + eps_P, P_back, X, Cd,
                                             area, beta);
  auto [mdot_P_minus, dummy10, dummy11, dummy12] =
      orifice_compressible_mdot_and_jacobian(T0, P0 - eps_P, P_back, X, Cd,
                                             area, beta);
  double fd_dP0 = (mdot_P_plus - mdot_P_minus) / (2.0 * eps_P);
  report_jacobian_difference("OrificeCompressibleMdotDerivatives", "dmdot_dP0",
                             dmdot_dP0, fd_dP0);
  EXPECT_NEAR(dmdot_dP0, fd_dP0, std::abs(fd_dP0) * 1e-6 + 1e-12);

  // Central FD for dmdot/dP_back
  auto [mdot_Pb_plus, dummy13, dummy14, dummy15] =
      orifice_compressible_mdot_and_jacobian(T0, P0, P_back + eps_P, X, Cd,
                                             area, beta);
  auto [mdot_Pb_minus, dummy16, dummy17, dummy18] =
      orifice_compressible_mdot_and_jacobian(T0, P0, P_back - eps_P, X, Cd,
                                             area, beta);
  double fd_dP_back = (mdot_Pb_plus - mdot_Pb_minus) / (2.0 * eps_P);
  report_jacobian_difference("OrificeCompressibleMdotDerivatives",
                             "dmdot_dP_back", dmdot_dP_back, fd_dP_back);
  EXPECT_NEAR(dmdot_dP_back, fd_dP_back, std::abs(fd_dP_back) * 1e-6 + 1e-12);
}

TEST(SolverJacobianTest, MachNumberDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> X(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      X[i] = 0.79;
    else if (species_name(i) == "O2")
      X[i] = 0.21;
  }

  double v = 200.0;
  double T = 300.0;

  auto [M, dM_dv] = mach_number_and_jacobian_v(v, T, X);

  // Central FD for dM/dv
  double eps_v = 1e-3;
  auto [M_plus, dummy1] = mach_number_and_jacobian_v(v + eps_v, T, X);
  auto [M_minus, dummy2] = mach_number_and_jacobian_v(v - eps_v, T, X);
  double fd_dM_dv = (M_plus - M_minus) / (2.0 * eps_v);
  report_jacobian_difference("MachNumberDerivatives", "dM_dv", dM_dv, fd_dM_dv);
  EXPECT_NEAR(dM_dv, fd_dM_dv, std::abs(fd_dM_dv) * 1e-6 + 1e-14);
}

TEST(SolverJacobianTest, AdiabaticWallDerivatives) {
  std::size_t ns = num_species();
  std::vector<double> X(ns, 0.0);
  for (size_t i = 0; i < ns; ++i) {
    if (species_name(i) == "N2")
      X[i] = 0.79;
    else if (species_name(i) == "O2")
      X[i] = 0.21;
  }

  double T_static = 350.0;
  double v = 150.0;
  double T = 300.0;
  double P = 101325.0;

  auto [T_aw, dT_aw_dv] = T_adiabatic_wall_and_jacobian_v(T_static, v, T, P, X);

  // Central FD for dT_aw/dv
  double eps_v = 1e-3;
  auto [T_aw_plus, dummy1] =
      T_adiabatic_wall_and_jacobian_v(T_static, v + eps_v, T, P, X);
  auto [T_aw_minus, dummy2] =
      T_adiabatic_wall_and_jacobian_v(T_static, v - eps_v, T, P, X);
  double fd_dT_aw_dv = (T_aw_plus - T_aw_minus) / (2.0 * eps_v);
  report_jacobian_difference("AdiabaticWallDerivatives", "dT_aw_dv", dT_aw_dv,
                             fd_dT_aw_dv);
  EXPECT_NEAR(dT_aw_dv, fd_dT_aw_dv, std::abs(fd_dT_aw_dv) * 1e-6 + 1e-11);
}

// ---------------------------------------------------------------------------
// Area Change Element Jacobians
// ---------------------------------------------------------------------------

TEST(SolverJacobianTest, SharpAreaChangeExpansionDerivatives) {
  double m_dot = 0.1;
  double rho = 1.2;
  double mu = 1.8e-5;
  double F0 = 0.01;   // small area
  double F1 = 0.04;   // large area (expansion for positive flow)
  double m_scale = 1e-4;

  auto res = combaero::sharp_area_change(m_dot, rho, mu, F0, F1, m_scale);

  // Central FD for dS/dm
  double eps_m = 1e-7;
  auto res_mp = combaero::sharp_area_change(m_dot + eps_m, rho, mu, F0, F1, m_scale);
  auto res_mm = combaero::sharp_area_change(m_dot - eps_m, rho, mu, F0, F1, m_scale);
  double fd_dS_dm = (res_mp.dP - res_mm.dP) / (2.0 * eps_m);
  report_jacobian_difference("SharpAreaChangeExpansion", "dS_dm", res.dS_dm, fd_dS_dm);
  EXPECT_NEAR(res.dS_dm, fd_dS_dm, std::abs(fd_dS_dm) * 1e-4 + 1e-10);

  // Central FD for dS/drho
  double eps_rho = 1e-6;
  auto res_rp = combaero::sharp_area_change(m_dot, rho + eps_rho, mu, F0, F1, m_scale);
  auto res_rm = combaero::sharp_area_change(m_dot, rho - eps_rho, mu, F0, F1, m_scale);
  double fd_dS_drho = (res_rp.dP - res_rm.dP) / (2.0 * eps_rho);
  report_jacobian_difference("SharpAreaChangeExpansion", "dS_drho", res.dS_drho, fd_dS_drho);
  EXPECT_NEAR(res.dS_drho, fd_dS_drho, std::abs(fd_dS_drho) * 1e-4 + 1e-10);
}

TEST(SolverJacobianTest, SharpAreaChangeContractionDerivatives) {
  double m_dot = 0.1;
  double rho = 1.2;
  double mu = 1.8e-5;
  double F0 = 0.04;   // large area
  double F1 = 0.01;   // small area (contraction for positive flow)
  double m_scale = 1e-4;

  auto res = combaero::sharp_area_change(m_dot, rho, mu, F0, F1, m_scale);

  // Central FD for dS/dm
  double eps_m = 1e-7;
  auto res_mp = combaero::sharp_area_change(m_dot + eps_m, rho, mu, F0, F1, m_scale);
  auto res_mm = combaero::sharp_area_change(m_dot - eps_m, rho, mu, F0, F1, m_scale);
  double fd_dS_dm = (res_mp.dP - res_mm.dP) / (2.0 * eps_m);
  report_jacobian_difference("SharpAreaChangeContraction", "dS_dm", res.dS_dm, fd_dS_dm);
  EXPECT_NEAR(res.dS_dm, fd_dS_dm, std::abs(fd_dS_dm) * 1e-4 + 1e-10);

  // Central FD for dS/drho
  double eps_rho = 1e-6;
  auto res_rp = combaero::sharp_area_change(m_dot, rho + eps_rho, mu, F0, F1, m_scale);
  auto res_rm = combaero::sharp_area_change(m_dot, rho - eps_rho, mu, F0, F1, m_scale);
  double fd_dS_drho = (res_rp.dP - res_rm.dP) / (2.0 * eps_rho);
  report_jacobian_difference("SharpAreaChangeContraction", "dS_drho", res.dS_drho, fd_dS_drho);
  EXPECT_NEAR(res.dS_drho, fd_dS_drho, std::abs(fd_dS_drho) * 1e-4 + 1e-10);
}

TEST(SolverJacobianTest, SharpAreaChangeReverseFlowDerivatives) {
  // Reverse flow through expansion geometry
  double m_dot = -0.1;
  double rho = 1.2;
  double mu = 1.8e-5;
  double F0 = 0.01;
  double F1 = 0.04;
  double m_scale = 1e-4;

  auto res = combaero::sharp_area_change(m_dot, rho, mu, F0, F1, m_scale);

  double eps_m = 1e-7;
  auto res_mp = combaero::sharp_area_change(m_dot + eps_m, rho, mu, F0, F1, m_scale);
  auto res_mm = combaero::sharp_area_change(m_dot - eps_m, rho, mu, F0, F1, m_scale);
  double fd_dS_dm = (res_mp.dP - res_mm.dP) / (2.0 * eps_m);
  report_jacobian_difference("SharpAreaChangeReverse", "dS_dm", res.dS_dm, fd_dS_dm);
  EXPECT_NEAR(res.dS_dm, fd_dS_dm, std::abs(fd_dS_dm) * 1e-4 + 1e-10);

  double eps_rho = 1e-6;
  auto res_rp = combaero::sharp_area_change(m_dot, rho + eps_rho, mu, F0, F1, m_scale);
  auto res_rm = combaero::sharp_area_change(m_dot, rho - eps_rho, mu, F0, F1, m_scale);
  double fd_dS_drho = (res_rp.dP - res_rm.dP) / (2.0 * eps_rho);
  report_jacobian_difference("SharpAreaChangeReverse", "dS_drho", res.dS_drho, fd_dS_drho);
  EXPECT_NEAR(res.dS_drho, fd_dS_drho, std::abs(fd_dS_drho) * 1e-4 + 1e-10);
}

TEST(SolverJacobianTest, SharpAreaChangeZeroFlowSmooth) {
  // At m_dot=0: verify smoothness and finite values.
  // Exact FD agreement is not expected here because smooth_abs(0) = sqrt(eps)
  // creates inherent O(sqrt(eps)) bias at the regularization scale.
  double rho = 1.2;
  double mu = 1.8e-5;
  double F0 = 0.01;
  double F1 = 0.04;
  double m_scale = 1e-4;

  auto res = combaero::sharp_area_change(0.0, rho, mu, F0, F1, m_scale);

  // No NaN or inf
  EXPECT_FALSE(std::isnan(res.dP));
  EXPECT_FALSE(std::isnan(res.dS_dm));
  EXPECT_FALSE(std::isnan(res.dS_drho));
  EXPECT_FALSE(std::isinf(res.dS_dm));

  // dP should be ~0 at zero flow
  EXPECT_NEAR(res.dP, 0.0, 1e-6);

  // dS/dm should be non-negative (more flow -> more loss)
  EXPECT_GE(res.dS_dm, 0.0);

  // Verify continuity: nearby points produce similar dP
  double m_tiny = 1e-8;
  auto res_p = combaero::sharp_area_change(+m_tiny, rho, mu, F0, F1, m_scale);
  auto res_m = combaero::sharp_area_change(-m_tiny, rho, mu, F0, F1, m_scale);
  EXPECT_NEAR(res_p.dP, 0.0, 1e-3);
  EXPECT_NEAR(res_m.dP, 0.0, 1e-3);
}

// ---------------------------------------------------------------------------
// Area change element kernel (area_change_flow, #731): every analytic partial
// against central differences, across the closure's regimes.
// ---------------------------------------------------------------------------

namespace {

std::vector<double> air_fuel_Y() {
  std::vector<double> X(combaero::num_species(), 0.0);
  for (std::size_t i = 0; i < X.size(); ++i) {
    if (combaero::species_name(i) == "N2") X[i] = 0.75;
    if (combaero::species_name(i) == "O2") X[i] = 0.20;
    if (combaero::species_name(i) == "CH4") X[i] = 0.05;
  }
  return combaero::mole_to_mass(X);
}

struct AreaChangeFlowCase {
  const char *name;
  double Pt0, Tt0, drive, F_in, F_out;
  const char *model;
  double length;
  bool incompressible;
  double hp = 0.0;  // pressure step [Pa]; 0: min(1e-3 drive, 1e-3)
};

void check_area_change_flow(const AreaChangeFlowCase &c) {
  const std::vector<double> Y = air_fuel_Y();
  const double Pt_t = c.Pt0 - c.drive;
  auto m = [&](double Pt0, double Tt0, const std::vector<double> &Yv, double Ptt) {
    return combaero::area_change_flow(Pt0, Tt0, Yv, Ptt, c.F_in, c.F_out, c.model,
                                      c.length, 0.0, c.incompressible, false)
        .m_dot;
  };
  const auto r = combaero::area_change_flow(c.Pt0, c.Tt0, Y, Pt_t, c.F_in, c.F_out,
                                            c.model, c.length, 0.0, c.incompressible, true);
  ASSERT_GT(std::abs(r.m_dot), 0.0) << c.name;
  const double tol = 2e-5;
  // Pressure steps well inside the drive (and the 1 Pa floor), T in K.
  const double hp = (c.hp > 0.0) ? c.hp : std::min(1e-3 * std::abs(c.drive), 1e-3);
  const double fd_Pt0 = (m(c.Pt0 + hp, c.Tt0, Y, Pt_t) - m(c.Pt0 - hp, c.Tt0, Y, Pt_t)) / (2 * hp);
  const double fd_Ptt = (m(c.Pt0, c.Tt0, Y, Pt_t + hp) - m(c.Pt0, c.Tt0, Y, Pt_t - hp)) / (2 * hp);
  const double hT = 1e-3;
  const double fd_Tt = (m(c.Pt0, c.Tt0 + hT, Y, Pt_t) - m(c.Pt0, c.Tt0 - hT, Y, Pt_t)) / (2 * hT);
  // Scales: the flow over the drive, the pressure and the temperature.
  const double sP = std::abs(r.m_dot) / std::max(std::abs(c.drive), 1e-12);
  const double sT = std::abs(r.m_dot) / c.Tt0;
  report_jacobian_difference(c.name, "dm_dPt0", r.dm_dPt0, fd_Pt0);
  report_jacobian_difference(c.name, "dm_dPt_target", r.dm_dPt_target, fd_Ptt);
  report_jacobian_difference(c.name, "dm_dTt0", r.dm_dTt0, fd_Tt);
  EXPECT_NEAR(r.dm_dPt0, fd_Pt0, tol * sP) << c.name;
  EXPECT_NEAR(r.dm_dPt_target, fd_Ptt, tol * sP) << c.name;
  EXPECT_NEAR(r.dm_dTt0, fd_Tt, tol * sT) << c.name;
  for (std::size_t k = 0; k < Y.size(); ++k) {
    const double hY = 1e-6;
    std::vector<double> Yp = Y, Ym = Y;
    Yp[k] += hY;
    double fd = 0.0;
    if (Y[k] > hY) {
      Ym[k] -= hY;
      fd = (m(c.Pt0, c.Tt0, Yp, Pt_t) - m(c.Pt0, c.Tt0, Ym, Pt_t)) / (2 * hY);
    } else {
      fd = (m(c.Pt0, c.Tt0, Yp, Pt_t) - r.m_dot) / hY;
    }
    EXPECT_NEAR(r.dm_dY[k], fd, 1e-4 * std::max(std::abs(r.m_dot), std::abs(fd)))
        << c.name << " Y[" << combaero::species_name(k) << "]";
  }
}

}  // namespace

TEST(AreaChangeFlowTest, DerivativesMatchCentralDifferences) {
  const double F_s = 7.853981633974483e-7;  // 1 mm
  const AreaChangeFlowCase cases[] = {
      {"expansion M 0.3", 2e5, 300.0, 7e3, 1e-3, 4e-3, "sharp", 0.0, false},
      {"expansion high Mach", 2e5, 600.0, 3e4, 1e-3, 2e-3, "sharp", 0.0, false},
      {"expansion Re ~100", 2e5, 300.0, 1.0, F_s, 4 * F_s, "sharp", 0.0, false},
      {"expansion Re ~2500", 2e5, 300.0, 300.0, F_s, 4 * F_s, "sharp", 0.0, false},
      // Laminar (Re ~ 0.4): linear in the drive, so a step of 40% of it is exact,
      // and well above the rounding of 1e5 Pa.
      {"expansion near zero drive", 1e5, 300.0, 1e-5, 1e-4, 4e-4, "sharp", 0.0, false, 4e-6},
      {"expansion reversed (odd)", 2e5, 300.0, -0.5, F_s, 4 * F_s, "sharp", 0.0, false},
      {"contraction", 2e5, 300.0, 5e3, 4e-3, 1e-3, "sharp", 0.0, false},
      {"contraction Re ~100", 2e5, 300.0, 1.0, 4 * F_s, F_s, "sharp", 0.0, false},
      {"contraction Re ~5000", 2e5, 300.0, 300.0, 4 * F_s, F_s, "sharp", 0.0, false},
      {"contraction below floor", 2e5, 300.0, 0.4, 4e-3, 1e-3, "sharp", 0.0, false},
      {"contraction choked", 2e5, 300.0, 1e5, 4e-3, 1e-3, "sharp", 0.0, false},
      {"conical diffuser", 2e5, 300.0, 3e3, 1e-3, 4e-3, "conical", 0.05, false},
      {"conical nozzle", 2e5, 300.0, 3e3, 4e-3, 1e-3, "conical", 0.05, false},
      {"expansion incompressible", 2e5, 300.0, 5e3, 1e-3, 4e-3, "sharp", 0.0, true},
      {"expansion incompressible Re ~2500", 2e5, 300.0, 300.0, F_s, 4 * F_s, "sharp", 0.0,
       true},
  };
  for (const auto &c : cases) {
    check_area_change_flow(c);
  }
}

TEST(AreaChangeFlowTest, ChokedFlowIsIndependentOfTheTarget) {
  const std::vector<double> Y = air_fuel_Y();
  const auto a = combaero::area_change_flow(2e5, 300.0, Y, 1.0e5, 4e-3, 1e-3, "sharp", 0.0, 0.0, false);
  const auto b = combaero::area_change_flow(2e5, 300.0, Y, 0.5e5, 4e-3, 1e-3, "sharp", 0.0, 0.0, false);
  EXPECT_TRUE(a.choked);
  EXPECT_DOUBLE_EQ(a.m_dot, b.m_dot);
  EXPECT_EQ(a.dm_dPt_target, 0.0);
}

// ---------------------------------------------------------------------------
// Conical Area Change -- low-level Jacobians
// ---------------------------------------------------------------------------

TEST(SolverJacobianTest, ConicalAreaChangeDiffuserDerivatives) {
  double m_dot = 0.1, rho = 1.2, mu = 1.8e-5;
  double F0 = 0.01, F1 = 0.04, L = 0.1, m_scale = 1e-4;

  auto res = combaero::conical_area_change(m_dot, rho, mu, F0, F1, L, m_scale);

  double eps_m = 1e-7;
  auto rp = combaero::conical_area_change(m_dot + eps_m, rho, mu, F0, F1, L, m_scale);
  auto rm = combaero::conical_area_change(m_dot - eps_m, rho, mu, F0, F1, L, m_scale);
  double fd_dm = (rp.dP - rm.dP) / (2.0 * eps_m);
  report_jacobian_difference("ConicalDiffuser", "dS_dm", res.dS_dm, fd_dm);
  EXPECT_NEAR(res.dS_dm, fd_dm, std::abs(fd_dm) * 1e-4 + 1e-10);

  double eps_rho = 1e-6;
  rp = combaero::conical_area_change(m_dot, rho + eps_rho, mu, F0, F1, L, m_scale);
  rm = combaero::conical_area_change(m_dot, rho - eps_rho, mu, F0, F1, L, m_scale);
  double fd_drho = (rp.dP - rm.dP) / (2.0 * eps_rho);
  report_jacobian_difference("ConicalDiffuser", "dS_drho", res.dS_drho, fd_drho);
  EXPECT_NEAR(res.dS_drho, fd_drho, std::abs(fd_drho) * 1e-4 + 1e-10);

  // dS_dmu must be exactly 0
  EXPECT_EQ(res.dS_dmu, 0.0);
}

TEST(SolverJacobianTest, ConicalAreaChangeNozzleDerivatives) {
  double m_dot = 0.1, rho = 1.2, mu = 1.8e-5;
  double F0 = 0.04, F1 = 0.01, L = 0.1, m_scale = 1e-4;

  auto res = combaero::conical_area_change(m_dot, rho, mu, F0, F1, L, m_scale);

  double eps_m = 1e-7;
  auto rp = combaero::conical_area_change(m_dot + eps_m, rho, mu, F0, F1, L, m_scale);
  auto rm = combaero::conical_area_change(m_dot - eps_m, rho, mu, F0, F1, L, m_scale);
  double fd_dm = (rp.dP - rm.dP) / (2.0 * eps_m);
  report_jacobian_difference("ConicalNozzle", "dS_dm", res.dS_dm, fd_dm);
  EXPECT_NEAR(res.dS_dm, fd_dm, std::abs(fd_dm) * 1e-4 + 1e-10);

  double eps_rho = 1e-6;
  rp = combaero::conical_area_change(m_dot, rho + eps_rho, mu, F0, F1, L, m_scale);
  rm = combaero::conical_area_change(m_dot, rho - eps_rho, mu, F0, F1, L, m_scale);
  double fd_drho = (rp.dP - rm.dP) / (2.0 * eps_rho);
  report_jacobian_difference("ConicalNozzle", "dS_drho", res.dS_drho, fd_drho);
  EXPECT_NEAR(res.dS_drho, fd_drho, std::abs(fd_drho) * 1e-4 + 1e-10);
}

TEST(SolverJacobianTest, ConicalAreaChangeReverseDerivatives) {
  double m_dot = -0.1, rho = 1.2, mu = 1.8e-5;
  double F0 = 0.01, F1 = 0.04, L = 0.1, m_scale = 1e-4;

  auto res = combaero::conical_area_change(m_dot, rho, mu, F0, F1, L, m_scale);

  double eps_m = 1e-7;
  auto rp = combaero::conical_area_change(m_dot + eps_m, rho, mu, F0, F1, L, m_scale);
  auto rm = combaero::conical_area_change(m_dot - eps_m, rho, mu, F0, F1, L, m_scale);
  double fd_dm = (rp.dP - rm.dP) / (2.0 * eps_m);
  report_jacobian_difference("ConicalReverse", "dS_dm", res.dS_dm, fd_dm);
  EXPECT_NEAR(res.dS_dm, fd_dm, std::abs(fd_dm) * 1e-4 + 1e-10);

  double eps_rho = 1e-6;
  rp = combaero::conical_area_change(m_dot, rho + eps_rho, mu, F0, F1, L, m_scale);
  rm = combaero::conical_area_change(m_dot, rho - eps_rho, mu, F0, F1, L, m_scale);
  double fd_drho = (rp.dP - rm.dP) / (2.0 * eps_rho);
  report_jacobian_difference("ConicalReverse", "dS_drho", res.dS_drho, fd_drho);
  EXPECT_NEAR(res.dS_drho, fd_drho, std::abs(fd_drho) * 1e-4 + 1e-10);
}

// ---------------------------------------------------------------------------
// Conical -> Sharp limit consistency
// ---------------------------------------------------------------------------

TEST(SolverJacobianTest, ConicalSharpLimitConsistency) {
  // At length -> 0 the cone is 180 deg: Idelchik 5-2 prints 0.58 at n 0.25,
  // Borda-Carnot (the sharp expansion at Re 6e4) is 0.5625.
  double m_dot = 0.1, rho = 1.2, mu = 1.8e-5;
  double F0 = 0.01, F1 = 0.04;

  auto sharp = combaero::sharp_area_change(m_dot, rho, mu, F0, F1);
  auto conical = combaero::conical_area_change(m_dot, rho, mu, F0, F1, 1e-6);

  double rel = std::abs(conical.dP - sharp.dP) / std::abs(sharp.dP);
  EXPECT_NEAR(rel, 0.58 / 0.5625 - 1.0, 1e-5) << "Conical(L->0) vs Sharp: " << conical.dP << " vs " << sharp.dP;
}

// Test fixture to print summary at the end
class JacobianTestEnvironment : public ::testing::Environment {
public:
  void TearDown() override { print_jacobian_difference_summary(); }
};

// Register the environment
::testing::Environment *const jacobian_env =
    ::testing::AddGlobalTestEnvironment(new JacobianTestEnvironment());

// An orifice between pipes (#645): the tap pressure from the recovered drop,
// every derivative against central differences, and the compressible floor.
TEST(OrificeTapPressure, DerivativesAreExactAndTheFloorIsCompressibleOnly) {
  const double x0[5] = {1.15e5, 1.17e5, 1.05e5, 0.6, 0.61};
  auto eval = [](const double *x, bool c) {
    return orifice_tap_pressure_and_jacobian(x[0], x[1], x[2], x[3], x[4], c);
  };
  for (bool comp : {false, true}) {
    const auto t = eval(x0, comp);
    EXPECT_FALSE(t.floored);
    const double analytic[5] = {t.d_dP_static_up, t.d_dP_total_up, t.d_dP_total_down,
                                0.0, t.d_dC};
    for (int i : {0, 1, 2, 4}) {
      double up[5], dn[5];
      std::copy(x0, x0 + 5, up);
      std::copy(x0, x0 + 5, dn);
      const double h = 1e-6 * std::abs(x0[i]);
      up[i] += h;
      dn[i] -= h;
      const double fd = (eval(up, comp).p_tap - eval(dn, comp).p_tap) / (2.0 * h);
      EXPECT_NEAR(analytic[i], fd, 1e-6 * std::max(1.0, std::abs(fd))) << i;
    }
  }
  // A drop the incompressible relation maps below zero.
  const auto c = orifice_tap_pressure_and_jacobian(1e5, 1e5, 0.4e5, 0.7, 0.61, true);
  const auto i = orifice_tap_pressure_and_jacobian(1e5, 1e5, 0.4e5, 0.7, 0.61, false);
  EXPECT_TRUE(c.floored);
  EXPECT_DOUBLE_EQ(c.p_tap, ORIFICE_TAP_PRESSURE_FLOOR * 1e5);
  EXPECT_EQ(c.d_dC, 0.0);
  EXPECT_FALSE(i.floored);
  EXPECT_LT(i.p_tap, 0.0);
}

// The static state of an orifice's approach pipe (#651): both closures'
// derivatives against central differences, and the momentum chamber's own
// residual zero at the returned pressure.
TEST(ApproachStatic, DerivativesAreExactAndItIsTheChambersState) {
  std::vector<double> Y(combaero::num_species(), 0.0);
  Y[combaero::species_index_from_name("N2")] = 0.767;
  Y[combaero::species_index_from_name("O2")] = 0.233;
  const double x0[3] = {1.5e5, 350.0, 120.0};
  for (bool inc : {true, false}) {
    auto P = [&](const double *x) { return approach_static_pressure(x[0], x[1], x[2], Y, inc).P; };
    const auto a = approach_static_pressure(x0[0], x0[1], x0[2], Y, inc);
    EXPECT_FALSE(a.choked);
    const double analytic[3] = {a.dP_dPt, a.dP_dTt, a.dP_dG};
    for (int i = 0; i < 3; ++i) {
      double up[3], dn[3];
      std::copy(x0, x0 + 3, up);
      std::copy(x0, x0 + 3, dn);
      const double h = 1e-6 * x0[i];
      up[i] += h;
      dn[i] -= h;
      const double fd = (P(up) - P(dn)) / (2.0 * h);
      EXPECT_NEAR(analytic[i], fd, 1e-6 * std::abs(fd)) << inc << " " << i;
    }
    const double A = 0.01;
    const auto r = momentum_chamber_residual_and_jacobian(a.P, x0[0], x0[2] * A, x0[1], Y, A, inc);
    EXPECT_NEAR(r.residual, 0.0, 1e-7 * x0[0]) << inc;
  }
}

// Past the approach pipe's sonic limit (#958): the compressible state is the
// sonic one, P* = P*(Pt, Tt), flat in G. Its Tt column was 0 (cp varies, so
// P*/Pt moves with Tt) and its value carried the golden-section noise of
// critical_pressure_ratio, which central differences read as a derivative.
TEST(ApproachStatic, BeyondTheSonicLimitTheDerivativesMatchTheValue) {
  std::vector<double> Y(combaero::num_species(), 0.0);
  Y[combaero::species_index_from_name("N2")] = 0.767;
  Y[combaero::species_index_from_name("O2")] = 0.233;
  for (double Tt : {400.0, 1900.0}) {
    const double x0[3] = {40e5, Tt, 1e5};
    auto P = [&](const double *x) {
      return approach_static_pressure(x[0], x[1], x[2], Y, false).P;
    };
    const auto a = approach_static_pressure(x0[0], x0[1], x0[2], Y, false);
    ASSERT_TRUE(a.choked);
    EXPECT_EQ(a.dP_dG, 0.0);
    EXPECT_GT(std::abs(a.dP_dTt), 0.0);
    const double analytic[2] = {a.dP_dPt, a.dP_dTt};
    for (int i = 0; i < 2; ++i) {
      double up[3], dn[3];
      std::copy(x0, x0 + 3, up);
      std::copy(x0, x0 + 3, dn);
      const double h = 1e-6 * x0[i];
      up[i] += h;
      dn[i] -= h;
      const double fd = (P(up) - P(dn)) / (2.0 * h);
      EXPECT_NEAR(analytic[i], fd, 1e-6 * std::abs(fd)) << Tt << " " << i;
    }
  }
}

// #964: the orifice kernels' and the approach state's composition columns are
// analytic; they were forward differences, off a central one by up to 7e-4 of
// the row. Checked against a central difference in each mass fraction (a
// one-sided second-order one for an absent species, which cannot go below 0).
namespace {

std::vector<double> burnt_mass_fractions() {
  std::vector<double> X(combaero::num_species(), 0.0);
  X[combaero::species_index_from_name("N2")] = 0.72;
  X[combaero::species_index_from_name("O2")] = 0.10;
  X[combaero::species_index_from_name("CO2")] = 0.06;
  X[combaero::species_index_from_name("H2O")] = 0.11;
  X[combaero::species_index_from_name("AR")] = 0.009;
  X[combaero::species_index_from_name("CH4")] = 0.001;
  return combaero::mole_to_mass(combaero::normalize_fractions(X));
}

template <typename F>
double y_difference(const F &f, const std::vector<double> &Y, std::size_t k) {
  const double h = 1e-6;
  std::vector<double> a = Y, b = Y;
  if (Y[k] > 2.0 * h) {
    a[k] += h;
    b[k] -= h;
    return (f(a) - f(b)) / (2.0 * h);
  }
  a[k] += h;
  b[k] += 2.0 * h;
  return (-3.0 * f(Y) + 4.0 * f(a) - f(b)) / (2.0 * h);
}

// Every column within tol of the largest central difference.
template <typename F>
void expect_y_columns(const F &f, const std::vector<double> &Y,
                      const std::vector<double> &analytic, double tol,
                      const std::string &what) {
  ASSERT_EQ(analytic.size(), Y.size()) << what;
  std::vector<double> fd(Y.size());
  double scale = 0.0;
  for (std::size_t k = 0; k < Y.size(); ++k) {
    fd[k] = y_difference(f, Y, k);
    scale = std::max(scale, std::abs(fd[k]));
  }
  ASSERT_GT(scale, 0.0) << what;
  for (std::size_t k = 0; k < Y.size(); ++k) {
    EXPECT_NEAR(analytic[k], fd[k], tol * scale) << what << " Y[" << k << "]";
  }
}

} // namespace

TEST(OrificeKernels, CompositionColumnsMatchCentralDifferences) {
  const std::vector<double> Y = burnt_mass_fractions();
  const double Pt = 10e5;
  const double Tt = 900.0;
  const double A = 7.85e-5;
  const double tol = 1e-6;
  {
    const double P1 = 9.5e5, Pb = 7e5, beta = 0.5;
    auto m = [&](const std::vector<double> &y) {
      return orifice_residuals_and_jacobian(0.0, Pt, P1, Tt, y, Pb, 0.62, A, beta)
          .m_dot_calc;
    };
    expect_y_columns(
        m, Y, orifice_residuals_and_jacobian(0.0, Pt, P1, Tt, Y, Pb, 0.62, A, beta).d_mdot_dY_up,
        tol, "incompressible");
  }
  // Compressible nozzle: plenum-fed and from the approach static state,
  // subsonic, choked and reversed.
  for (double beta : {0.0, 0.6}) {
    for (double Pb : {7e5, 3e5, 11e5}) {
      const double P1 = beta > 0.0 ? 9.5e5 : Pt;
      auto m = [&](const std::vector<double> &y) {
        return orifice_compressible_residuals_and_jacobian(0.0, Pt, P1, Tt, y, Pb, 0.62, A, beta)
            .m_dot_calc;
      };
      const auto r = orifice_compressible_residuals_and_jacobian(0.0, Pt, P1, Tt, Y, Pb, 0.62, A,
                                                                 beta);
      expect_y_columns(m, Y, r.d_mdot_dY_up, tol,
                       "nozzle beta " + std::to_string(beta) + " Pb " + std::to_string(Pb));
    }
  }
  // Jobson and ISO 5167, composition through gamma(Tt) and R.
  for (double Pb : {7e5, 3e5}) {
    const double P1 = 9.5e5, beta = 0.5;
    auto mj = [&](const std::vector<double> &y) {
      return orifice_jobson_residuals_and_jacobian(0.0, Pt, P1, Tt, y, Pb, 0.62, A, beta)
          .m_dot_calc;
    };
    expect_y_columns(
        mj, Y,
        orifice_jobson_residuals_and_jacobian(0.0, Pt, P1, Tt, Y, Pb, 0.62, A, beta).d_mdot_dY_up,
        tol, "jobson Pb " + std::to_string(Pb));
    auto mi = [&](const std::vector<double> &y) {
      return orifice_iso5167_residuals_and_jacobian(0.0, Pt, P1, Tt, y, Pb, 0.62, A, beta)
          .m_dot_calc;
    };
    expect_y_columns(
        mi, Y,
        orifice_iso5167_residuals_and_jacobian(0.0, Pt, P1, Tt, Y, Pb, 0.62, A, beta).d_mdot_dY_up,
        tol, "iso5167 Pb " + std::to_string(Pb));
  }
  // The approach pipe's static pressure: low-Mach closure, subsonic, choked.
  for (bool incompressible : {true, false}) {
    for (double G : {300.0, 1000.0, 1e5}) {
      if (incompressible && G > 2000.0) {
        continue;
      }
      auto p = [&](const std::vector<double> &y) {
        return approach_static_pressure(Pt, Tt, G, y, incompressible).P;
      };
      const auto a = approach_static_pressure(Pt, Tt, G, Y, incompressible);
      EXPECT_EQ(a.choked, G > 2000.0) << G;
      expect_y_columns(p, Y, a.dP_dY, tol,
                       "approach G " + std::to_string(G) + " inc " + std::to_string(incompressible));
    }
  }
}
