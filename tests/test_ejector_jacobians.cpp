// Jacobian-accuracy check for the ejector C++ port (include/ejector.h):
// compares the analytic d/dp_g, d/dt_g, d/dp_e, d/dt_e against the
// pre-baked central-difference targets (h = 1e-6*|x|) in
// validation/ejector/data/huang1999_reference_data.h, for all 31 rows.
// Mirrors tests/test_solver_jacobians.cpp's OrificeCompressibleDerivatives
// EXPECT_NEAR idiom.

#include "ejector.h"
#include "validation/ejector/data/huang1999_reference_data.h"
#include <cmath>
#include <gtest/gtest.h>

using namespace combaero::solver;
using combaero::validation::ejector::kHuangCases;

namespace {

// Relative tolerance for analytic-vs-CD comparison. The reference file's
// own CD step (h = 1e-6*|x|) has truncation/round-off noise around
// 1e-8-1e-7 relative on these expressions (confirmed in
// scripts/derive_ejector_jacobian.py's numeric cross-check); 1e-5 leaves
// comfortable margin while still catching a real derivation error.
constexpr double kRelTol = 1.0e-5;
constexpr double kAbsFloor = 1.0e-9;

void expect_jacobian_close(double analytic, double target, const char* what,
                           const char* case_id) {
  double tol = kRelTol * std::abs(target) + kAbsFloor;
  EXPECT_NEAR(analytic, target, tol) << what << " mismatch for case " << case_id
                                     << " (analytic=" << analytic << ", target=" << target << ")";
}

} // namespace

TEST(EjectorJacobians, EntrainmentRatioMatchesCentralDifferenceForAllCases) {
  for (const auto& c : kHuangCases) {
    EjectorGeometry geom{c.area_ratio_nozzle, c.area_ratio_mix};
    EjectorEntrainmentJacobian jac =
        ejector_entrainment_ratio_and_jacobian(c.p_g_pa, c.t_g_k, c.p_e_pa, c.t_e_k, geom, c.gamma);

    expect_jacobian_close(jac.domega_dp_g, c.domega_dp_g, "domega_dp_g", c.id);
    expect_jacobian_close(jac.domega_dt_g, c.domega_dt_g, "domega_dt_g", c.id);
    expect_jacobian_close(jac.domega_dp_e, c.domega_dp_e, "domega_dp_e", c.id);
    expect_jacobian_close(jac.domega_dt_e, c.domega_dt_e, "domega_dt_e", c.id);
  }
}

TEST(EjectorJacobians, CriticalBackPressureMatchesCentralDifferenceForAllCases) {
  for (const auto& c : kHuangCases) {
    EjectorGeometry geom{c.area_ratio_nozzle, c.area_ratio_mix};
    EjectorCriticalPressureJacobian jac = ejector_critical_back_pressure_and_jacobian(
        c.p_g_pa, c.t_g_k, c.p_e_pa, c.t_e_k, geom, c.gamma, c.r_gas, c.recovery_efficiency);

    expect_jacobian_close(jac.dpc_dp_g, c.dpc_dp_g, "dpc_dp_g", c.id);
    expect_jacobian_close(jac.dpc_dt_g, c.dpc_dt_g, "dpc_dt_g", c.id);
    expect_jacobian_close(jac.dpc_dp_e, c.dpc_dp_e, "dpc_dp_e", c.id);
    expect_jacobian_close(jac.dpc_dt_e, c.dpc_dt_e, "dpc_dt_e", c.id);
  }
}

TEST(EjectorJacobians, DeadHeadBackPressureMatchesCentralDifferenceForAllCases) {
  for (const auto& c : kHuangCases) {
    EjectorGeometry geom{c.area_ratio_nozzle, c.area_ratio_mix};
    EjectorCriticalPressureJacobian jac = ejector_dead_head_back_pressure_and_jacobian(
        c.p_g_pa, c.t_g_k, c.p_e_pa, c.t_e_k, geom, c.gamma, c.r_gas, c.recovery_efficiency);

    expect_jacobian_close(jac.dpc_dp_g, c.dpb0_dp_g, "dpb0_dp_g", c.id);
    // At zero entrainment p03 = p_g q(lambda1)/q(lambda3) and lambda1 = v/c*,
    // both ~sqrt(t_g): dP_b0/dt_g is exactly 0. The golden central difference
    // is its noise (~5e-7 Pa/K against a natural scale P_b0/t_g ~ 300 Pa/K), so
    // compare on that scale.
    const double scale_t = c.p_b0_pa / c.t_g_k;
    EXPECT_NEAR(jac.dpc_dt_g, c.dpb0_dt_g, 1e-8 * scale_t) << "dpb0_dt_g for case " << c.id;
    EXPECT_NEAR(jac.dpc_dt_g, 0.0, 1e-10 * scale_t) << "dpb0_dt_g not ~0 for case " << c.id;
    expect_jacobian_close(jac.dpc_dp_e, c.dpb0_dp_e, "dpb0_dp_e", c.id);
    expect_jacobian_close(jac.dpc_dt_e, c.dpb0_dt_e, "dpb0_dt_e", c.id);
  }
}

TEST(EjectorJacobians, ChokedMassFlowMatchesCentralDifference) {
  // huang1999_reference_data.h has no column for choked_mass_flow (the
  // reference model never needed the absolute mass flow, only the ratio
  // form) -- CD spot-check computed directly here instead.
  double gamma = 1.3;
  double r_gas = 287.0;
  double area_throat = 5.0e-5;
  double eta = 0.95;

  for (auto [p0, t0] : {std::pair{4.0e5, 400.0}, std::pair{1.2e5, 320.0}}) {
    auto [mdot, d_dp0, d_dt0] =
        ejector_choked_mass_flow_and_jacobian(p0, t0, area_throat, gamma, r_gas, eta);

    double eps_p = 1.0e-6 * p0;
    double fd_dp0 = (ejector_choked_mass_flow(p0 + eps_p, t0, area_throat, gamma, r_gas, eta) -
                     ejector_choked_mass_flow(p0 - eps_p, t0, area_throat, gamma, r_gas, eta)) /
                    (2.0 * eps_p);
    double eps_t = 1.0e-6 * t0;
    double fd_dt0 = (ejector_choked_mass_flow(p0, t0 + eps_t, area_throat, gamma, r_gas, eta) -
                     ejector_choked_mass_flow(p0, t0 - eps_t, area_throat, gamma, r_gas, eta)) /
                    (2.0 * eps_t);

    EXPECT_NEAR(d_dp0, fd_dp0, std::abs(fd_dp0) * 1e-8 + 1e-12);
    EXPECT_NEAR(d_dt0, fd_dt0, std::abs(fd_dt0) * 1e-8 + 1e-12);
    EXPECT_GT(mdot, 0.0);
  }
}

TEST(EjectorJacobians, CDNozzleAnalyticMatchesCentralDifference) {
  // C-D nozzle R0 residual: analytic (p0, t0, p_py) partials vs central
  // difference of the value function. Covers unchoked, near-threshold, and
  // choked (flat) points.
  const double A_t = 3.14e-5, A_e = 1.0e-4, g = 1.40, R = 287.0, eta = 0.95, eps = 1e-3;
  struct Pt {
    double p0, t0, p_py;
  };
  const Pt pts[] = {
      {101325.0, 300.0, 100207.0}, // unchoked (reported root)
      {101325.0, 300.0, 98700.0},  // near the choke threshold
      {101325.0, 300.0, 60000.0},  // choked (flat region)
      {604000.0, 368.0, 300000.0}, // higher motive, unchoked
  };
  auto val = [&](double p0, double t0, double p_py) {
    return ejector_cd_nozzle_mass_flow(p0, t0, p_py, A_t, A_e, g, R, eta, eps);
  };
  for (const auto& p : pts) {
    EjectorCDNozzleJacobian j =
        ejector_cd_nozzle_mass_flow_and_jacobian(p.p0, p.t0, p.p_py, A_t, A_e, g, R, eta, eps);
    double hp0 = 1e-6 * p.p0, ht0 = 1e-6 * p.t0, hpy = 1e-6 * p.p_py;
    double cd_p0 = (val(p.p0 + hp0, p.t0, p.p_py) - val(p.p0 - hp0, p.t0, p.p_py)) / (2 * hp0);
    double cd_t0 = (val(p.p0, p.t0 + ht0, p.p_py) - val(p.p0, p.t0 - ht0, p.p_py)) / (2 * ht0);
    double cd_py = (val(p.p0, p.t0, p.p_py + hpy) - val(p.p0, p.t0, p.p_py - hpy)) / (2 * hpy);
    EXPECT_NEAR(j.dmdot_dp0, cd_p0, 1e-6 * std::abs(cd_p0) + 1e-12) << "dp0 at p_py=" << p.p_py;
    EXPECT_NEAR(j.dmdot_dt0, cd_t0, 1e-6 * std::abs(cd_t0) + 1e-12) << "dt0 at p_py=" << p.p_py;
    EXPECT_NEAR(j.dmdot_dp_py, cd_py, 1e-5 * std::abs(cd_py) + 1e-12) << "dp_py at p_py=" << p.p_py;
  }
}

TEST(EjectorJacobians, JetPumpDischargeValueAndJacobian) {
  // R3 building block: value parity vs the Python mixing closure, and analytic
  // partials (p_g, t_g, p_e, t_e, p_py, omega) vs central difference.
  struct Case {
    double p_g, t_g, p_e, t_e, p_py, omega, gamma, rec, expected;
  };
  const Case cs[] = {
      {101325.0, 300.0, 101325.0, 300.0, 100207.0, 7.0, 1.4, 1.0, 101325.00000000015},
      {102435.0, 300.0, 101325.0, 288.15, 100000.0, 3.5, 1.4, 1.0, 101544.85636197214},
      {110000.0, 320.0, 100000.0, 300.0, 95000.0, 1.5, 1.4, 0.9, 92957.65455779647},
  };
  auto val = [](double pg, double tg, double pe, double te, double ppy, double om, double gm,
                double rec) {
    return ejector_jetpump_discharge_and_jacobian(pg, tg, pe, te, ppy, om, gm, rec).p03;
  };
  for (const auto& c : cs) {
    auto j = ejector_jetpump_discharge_and_jacobian(c.p_g, c.t_g, c.p_e, c.t_e, c.p_py, c.omega,
                                                    c.gamma, c.rec);
    EXPECT_NEAR(j.p03, c.expected, 1e-7 * c.expected) << "value at p_py=" << c.p_py;
    double h[6] = {1e-6 * c.p_g, 1e-6 * c.t_g, 1e-6 * c.p_e,
                   1e-6 * c.t_e, 1e-6 * c.p_py, 1e-6 * c.omega};
    double an[6] = {j.dp03_dp_g, j.dp03_dt_g, j.dp03_dp_e, j.dp03_dt_e, j.dp03_dp_py, j.dp03_domega};
    for (int k = 0; k < 6; ++k) {
      double a[8] = {c.p_g, c.t_g, c.p_e, c.t_e, c.p_py, c.omega, c.gamma, c.rec};
      double b[8];
      for (int m = 0; m < 8; ++m) b[m] = a[m];
      a[k] += h[k];
      b[k] -= h[k];
      double cd = (val(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]) -
                   val(b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7])) /
                  (2 * h[k]);
      // Abs floor 1e-5 tolerates central-difference round-off on partials that
      // are genuinely ~0 at the symmetric degenerate point (p_g=p_e, t_g=t_e);
      // real partials here are O(0.01)-O(1000), caught by the relative term.
      EXPECT_NEAR(an[k], cd, 1e-5 * std::abs(cd) + 1e-5) << "seed " << k << " at p_py=" << c.p_py;
    }
  }
}

// The nozzles' efficiencies in the mixing (#588): each stream mixes at the
// velocity its lossy nozzle delivers. At zero drive (both supplies at the same
// Pt and T) the isentropic mixing returned the supply Pt exactly -- the nozzle
// losses handed back -- while the lossy streams must arrive below it. And the
// partials stay analytic.
TEST(EjectorJacobians, JetPumpDischargeChargesTheNozzleLosses) {
  const double P = 101325.0, T = 300.0, ppy = 100207.0, om = 7.0, g = 1.4;
  const double lossless = ejector_jetpump_discharge_and_jacobian(P, T, P, T, ppy, om, g, 1.0).p03;
  EXPECT_NEAR(lossless, P, 1e-6);
  const double lossy =
      ejector_jetpump_discharge_and_jacobian(P, T, P, T, ppy, om, g, 1.0, 0.95, 0.85).p03;
  EXPECT_LT(lossy, P - 50.0);  // ~1.1 kPa of head behind it: a real loss

  auto val = [&](const double* a) {
    return ejector_jetpump_discharge_and_jacobian(a[0], a[1], a[2], a[3], a[4], a[5], g, 0.9,
                                                  0.95, 0.85)
        .p03;
  };
  const double c[6] = {110000.0, 320.0, 100000.0, 300.0, 95000.0, 1.5};
  const auto j = ejector_jetpump_discharge_and_jacobian(c[0], c[1], c[2], c[3], c[4], c[5], g,
                                                        0.9, 0.95, 0.85);
  const double an[6] = {j.dp03_dp_g, j.dp03_dt_g, j.dp03_dp_e,
                        j.dp03_dt_e, j.dp03_dp_py, j.dp03_domega};
  for (int k = 0; k < 6; ++k) {
    double a[6], b[6];
    for (int m = 0; m < 6; ++m) a[m] = b[m] = c[m];
    const double h = 1e-6 * c[k];
    a[k] += h;
    b[k] -= h;
    const double cd = (val(a) - val(b)) / (2 * h);
    EXPECT_NEAR(an[k], cd, 1e-5 * std::abs(cd) + 1e-5) << "seed " << k;
  }
}

TEST(EjectorJacobians, CDNozzleExitStaticImplicitDerivatives) {
  // Inverse of the C-D nozzle: analytic implicit-function partials
  // (m_dot, p0, t0) vs central difference of the bisection value. Also checks
  // the round-trip cd_nozzle(exit_static(m_dot)) == m_dot.
  const double A_t = 3.14e-5, A_e = 1.0e-4, g = 1.40, R = 287.0, eta = 0.95, eps = 1e-3;
  struct Pt {
    double m_dot, p0, t0;
  };
  const Pt pts[] = {
      {0.0051, 101325.0, 300.0}, // reported
      {0.0030, 101325.0, 300.0}, // deeper unchoked
      {0.0065, 101325.0, 300.0}, // near choke
      {0.0300, 604000.0, 368.0}, // higher motive
  };
  auto val = [&](double md, double p0, double t0) {
    return ejector_cd_nozzle_exit_static_and_jacobian(md, p0, t0, A_t, A_e, g, R, eta, eps).p_py;
  };
  for (const auto& p : pts) {
    auto j = ejector_cd_nozzle_exit_static_and_jacobian(p.m_dot, p.p0, p.t0, A_t, A_e, g, R, eta,
                                                        eps);
    // round-trip
    double back = ejector_cd_nozzle_mass_flow(p.p0, p.t0, j.p_py, A_t, A_e, g, R, eta, eps);
    EXPECT_NEAR(back, p.m_dot, 1e-9 * p.m_dot) << "round-trip at m_dot=" << p.m_dot;
    // implicit derivatives vs FD
    double hm = 1e-7 * p.m_dot, hp = 1e-6 * p.p0, ht = 1e-6 * p.t0;
    double cd_m = (val(p.m_dot + hm, p.p0, p.t0) - val(p.m_dot - hm, p.p0, p.t0)) / (2 * hm);
    double cd_p = (val(p.m_dot, p.p0 + hp, p.t0) - val(p.m_dot, p.p0 - hp, p.t0)) / (2 * hp);
    double cd_t = (val(p.m_dot, p.p0, p.t0 + ht) - val(p.m_dot, p.p0, p.t0 - ht)) / (2 * ht);
    EXPECT_NEAR(j.dp_py_dm_dot, cd_m, 1e-4 * std::abs(cd_m) + 1e-6) << "dm at " << p.m_dot;
    EXPECT_NEAR(j.dp_py_dp0, cd_p, 1e-4 * std::abs(cd_p) + 1e-9) << "dp0 at " << p.m_dot;
    EXPECT_NEAR(j.dp_py_dt0, cd_t, 1e-4 * std::abs(cd_t) + 1e-9) << "dt0 at " << p.m_dot;
  }
}

TEST(EjectorJacobians, ElementResidualsMatchCentralDifference) {
  // Whole-element (f, J): every entry of the 4x9 analytic Jacobian vs a
  // central difference of the residual value, at operating points inside each
  // regime. This directly guards the DualN<9> assembly in
  // ejector_element_residuals_and_jacobian (the C++ home of the element
  // assembly; 1:1 with EjectorElement.residuals()). Points are chosen
  // ASYMMETRIC (p_g != p_e or t_g != t_e) and away from the s_choke/s_sub
  // seams so no partial is genuinely ~0 -- the same reason the Python FD test
  // avoids the symmetric degenerate root.
  const double A_t = 3.14e-5, A_e = 1.0e-4, A_mix = 8.0e-4, A_s = A_mix - A_e;
  const EjectorGeometry geom{A_e / A_t, A_mix / A_t};
  const double R = 287.0, ep = 0.95, es = 0.85, rec = 1.0, eps = 1e-3;
  const double sc_lo = 0.90, sc_hi = 0.999, ss_lo = 0.98, ss_hi = 1.02;
  // gamma is a frozen coefficient here (matches the element's frozen-property
  // Jacobian); use a fixed representative air value so the FD holds it fixed.
  const double g = 1.4;

  struct Case {
    const char* id;
    double u[9]; // mp, ms, mdot_out, p_g, t_g, p_e, t_e, p_out, p_py
  };
  const Case cases[] = {
      // Unchoked jet-pump regime (s_choke -> 0), asymmetric streams.
      {"jetpump", {0.0051, 0.030, 0.0351, 105000.0, 305.0, 100000.0, 295.0, 100500.0, 99000.0}},
      // Deeper unchoked, different values.
      {"jetpump2", {0.0040, 0.022, 0.0260, 130000.0, 310.0, 101325.0, 292.0, 101000.0, 98500.0}},
      // Double-choked critical regime (s_choke -> 1, outlet < P_c* -> s_sub 0).
      {"critical", {0.44, 0.13, 0.57, 700000.0, 440.0, 22850.0, 280.0, 40000.0, 100000.0}},
      // Choked subcritical (#619): outlet between P_c* = 58.0 kPa and the
      // dead-head P_b0 = 175.9 kPa -> the omega_sub droop, outlet free.
      {"subcritical", {0.44, 0.20, 0.64, 700000.0, 440.0, 22850.0, 280.0, 100000.0, 60000.0}},
      // Choked, outlet beyond P_b0 -> the outlet pin at recovery * P_b0.
      {"beyond_pb0", {0.44, 0.05, 0.49, 700000.0, 440.0, 22850.0, 280.0, 200000.0, 60000.0}},
  };

  auto resid = [&](const double u[9], int row) {
    return ejector_element_residuals_and_jacobian(u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7],
                                                  u[8], geom, A_t, A_e, A_s, g, R, ep, es, rec, eps,
                                                  sc_lo, sc_hi, ss_lo, ss_hi)
        .residuals[row];
  };

  for (const auto& c : cases) {
    auto rj = ejector_element_residuals_and_jacobian(c.u[0], c.u[1], c.u[2], c.u[3], c.u[4], c.u[5],
                                                     c.u[6], c.u[7], c.u[8], geom, A_t, A_e, A_s, g,
                                                     R, ep, es, rec, eps, sc_lo, sc_hi, ss_lo, ss_hi);
    for (int row = 0; row < 4; ++row) {
      for (int k = 0; k < 9; ++k) {
        double a[9], b[9];
        for (int m = 0; m < 9; ++m) a[m] = b[m] = c.u[m];
        double h = 1e-6 * std::abs(c.u[k]);
        if (h == 0.0) h = 1e-6;
        a[k] += h;
        b[k] -= h;
        double cd = (resid(a, row) - resid(b, row)) / (2 * h);
        double an = rj.jacobian[row][k];
        // Residuals span mass rows (O(0.01-1)) and pressure rows (O(1e4-1e5)),
        // so the abs floor is scaled to the residual magnitude to tolerate CD
        // round-off on the pressure rows without masking mass-row errors.
        double floor = 1e-6 * std::abs(rj.residuals[row]) + 1e-7;
        EXPECT_NEAR(an, cd, 1e-5 * std::abs(cd) + floor)
            << "case " << c.id << " row " << row << " seed " << k;
      }
    }
  }
}


// -----------------------------------------------------------------------------
// The s_choke blend band (issue #306)
// -----------------------------------------------------------------------------
//
// The element test above chooses its points "away from the s_choke/s_sub seams
// so no partial is genuinely ~0". That is right for its purpose and is exactly
// why `ejector_smootherstep`'s INTERIOR was executed by nothing: every case sat
// on a saturated end, where the function returns a constant 0 or 1 and its
// derivative term is never evaluated. Measured with llvm-cov, the three lines
// computing `val`, `dval` and the chained result had zero coverage, so a
// shipping analytic derivative had no cross-check at all -- which the repo's
// solver rule forbids.
//
// The seam-avoidance does not apply here, and inverting it is the point.
// Inside the band `s` is strictly between 0 and 1 and its slope is nonzero, so
// the partials are LIVE rather than degenerate; it is the saturated ends that
// zero them. A band point is therefore a better finite-difference test than
// the ones above, not a worse one.
//
// The thresholds are arguments, so rather than hunting for a state that
// happens to land in a fixed band, the band is placed around the state's own
// mp/cap ratio -- 0.6857 for the "jetpump" point.
TEST(EjectorJacobians, SmootherstepBandIsExercisedAndItsDerivativeIsCorrect) {
  const double A_t = 3.14e-5, A_e = 1.0e-4, A_mix = 8.0e-4, A_s = A_mix - A_e;
  const EjectorGeometry geom{A_e / A_t, A_mix / A_t};
  const double R = 287.0, ep = 0.95, es = 0.85, rec = 1.0, eps = 1e-3, g = 1.4;
  const double ss_lo = 0.98, ss_hi = 1.02;

  // mp/cap = 0.6857 at this state.
  const double u[9] = {0.0051, 0.030, 0.0351, 105000.0, 305.0, 100000.0, 295.0, 100500.0, 99000.0};

  auto evaluate = [&](const double v[9], double sc_lo, double sc_hi) {
    return ejector_element_residuals_and_jacobian(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7],
                                                  v[8], geom, A_t, A_e, A_s, g, R, ep, es, rec, eps,
                                                  sc_lo, sc_hi, ss_lo, ss_hi);
  };

  // Three bands: one that saturates s to 0, one that saturates it to 1, and
  // one straddling the state's ratio.
  const double kBandLo = 0.60, kBandHi = 0.75;
  auto below = evaluate(u, 0.90, 0.999);  // ratio 0.686 < lo  -> s = 0
  auto above = evaluate(u, 0.30, 0.50);   // ratio 0.686 > hi  -> s = 1
  auto band = evaluate(u, kBandLo, kBandHi);

  // The band result must sit strictly between the two saturated ones on at
  // least one row, or the interior is not being reached and everything below
  // would be testing a saturated end again.
  bool strictly_between = false;
  for (int row = 0; row < 4; ++row) {
    double a = below.residuals[row];
    double b = above.residuals[row];
    double m = band.residuals[row];
    if (std::abs(a - b) < 1e-9) continue;
    double lo = std::min(a, b);
    double hi = std::max(a, b);
    if (m > lo + 1e-9 * std::abs(hi - lo) && m < hi - 1e-9 * std::abs(hi - lo)) {
      strictly_between = true;
    }
  }
  EXPECT_TRUE(strictly_between)
      << "the blend is saturated, so the smootherstep interior is still unreached";

  // And the analytic Jacobian must be right there, which is the cross-check
  // the derivative never had.
  auto resid = [&](const double v[9], int row) { return evaluate(v, kBandLo, kBandHi).residuals[row]; };
  for (int row = 0; row < 4; ++row) {
    for (int k = 0; k < 9; ++k) {
      double a[9], b[9];
      for (int m = 0; m < 9; ++m) a[m] = b[m] = u[m];
      double h = 1e-6 * std::abs(u[k]);
      if (h == 0.0) h = 1e-6;
      a[k] += h;
      b[k] -= h;
      double cd = (resid(a, row) - resid(b, row)) / (2 * h);
      double an = band.jacobian[row][k];
      double floor = 1e-6 * std::abs(band.residuals[row]) + 1e-7;
      EXPECT_NEAR(an, cd, 1e-5 * std::abs(cd) + floor) << "row " << row << " seed " << k;
    }
  }
}

// -----------------------------------------------------------------------------
// The "huang1999" anchor closures: Huang's printed P_c* chain and the
// sudden-expansion dead head (validation/ejector/README.md). No golden central
// differences are baked for them, so the difference is taken here.
// -----------------------------------------------------------------------------
namespace {

using AnchorFn = EjectorCriticalPressureJacobian (*)(double, double, double, double,
                                                     const EjectorGeometry&, double, double,
                                                     double);

void expect_anchor_jacobian(AnchorFn fn, const char* name) {
  for (const auto& c : kHuangCases) {
    const EjectorGeometry geom{c.area_ratio_nozzle, c.area_ratio_mix};
    const double x0[4] = {c.p_g_pa, c.t_g_k, c.p_e_pa, c.t_e_k};
    auto at = [&](const double x[4]) {
      return fn(x[0], x[1], x[2], x[3], geom, c.gamma, c.r_gas, 1.0).value.p_c_pa;
    };
    const auto jac = fn(x0[0], x0[1], x0[2], x0[3], geom, c.gamma, c.r_gas, 1.0);
    const double analytic[4] = {jac.dpc_dp_g, jac.dpc_dt_g, jac.dpc_dp_e, jac.dpc_dt_e};
    for (int i = 0; i < 4; ++i) {
      double a[4], b[4];
      for (int m = 0; m < 4; ++m) a[m] = b[m] = x0[m];
      const double h = 1e-6 * std::abs(x0[i]);
      a[i] += h;
      b[i] -= h;
      const double cd = (at(a) - at(b)) / (2.0 * h);
      // J_i |x_i| against the value (the residual-row-scaling metric).
      EXPECT_NEAR(analytic[i] * x0[i], cd * x0[i], 1e-6 * jac.value.p_c_pa)
          << name << " d/dx[" << i << "] for case " << c.id;
    }
  }
}

} // namespace

TEST(EjectorJacobians, HuangCriticalBackPressureMatchesCentralDifferenceForAllCases) {
  expect_anchor_jacobian(&ejector_huang_critical_back_pressure_and_jacobian, "huang P_c*");
}

TEST(EjectorJacobians, SuddenExpansionDeadHeadMatchesCentralDifferenceForAllCases) {
  expect_anchor_jacobian(&ejector_sudden_expansion_dead_head_and_jacobian, "sudden-expansion P_b0");
}

TEST(EjectorPhysics, SuddenExpansionDeadHeadIsANormalShockWithNoAnnulus) {
  // A_3 = A_p1: no suction annulus, so the sudden expansion is a normal shock
  // at the nozzle-exit Mach and P_b0 the stagnation pressure behind it.
  const double g = 1.4, arn = 1.86, p_g = 5.0e5;
  const auto res = ejector_sudden_expansion_dead_head_and_jacobian(
      p_g, 293.15, 101325.0, 293.15, EjectorGeometry{arn, arn}, g, 287.05, 1.0);
  const double m1 = ejector_mach_from_area_ratio_supersonic(arn, g);
  const double m1s = m1 * m1;
  const double p02_p01 = std::pow((g + 1) * m1s / ((g - 1) * m1s + 2), g / (g - 1)) *
                         std::pow((g + 1) / (2 * g * m1s - (g - 1)), 1 / (g - 1));
  EXPECT_NEAR(res.value.p_c_pa, p_g * p02_p01, 1e-9 * p_g);
}

TEST(EjectorJacobians, ElementResidualsMatchCentralDifferenceWithHuangAnchors) {
  // The element above, with the (P_c*, P_b0) pair switched: plateau, droop and
  // the pin above P_b0 (P_c* ~ 37 kPa, P_b0 ~ 59 kPa at this supply).
  const double A_t = 3.14e-5, A_e = 1.0e-4, A_mix = 8.0e-4, A_s = A_mix - A_e;
  const EjectorGeometry geom{A_e / A_t, A_mix / A_t};
  const double R = 287.0, ep = 0.95, es = 0.85, rec = 1.0, eps = 1e-3, g = 1.4;
  const double sc_lo = 0.90, sc_hi = 0.999, ss_lo = 0.98, ss_hi = 1.02;
  const auto huang = EjectorAnchorClosure::huang1999;
  struct Case {
    const char* id;
    double u[9];
  };
  const Case cases[] = {
      {"critical", {0.44, 0.13, 0.57, 700000.0, 440.0, 22850.0, 280.0, 30000.0, 100000.0}},
      {"subcritical", {0.44, 0.20, 0.64, 700000.0, 440.0, 22850.0, 280.0, 52000.0, 60000.0}},
      {"beyond_pb0", {0.44, 0.05, 0.49, 700000.0, 440.0, 22850.0, 280.0, 100000.0, 60000.0}},
  };
  auto rj_at = [&](const double u[9]) {
    return ejector_element_residuals_and_jacobian(u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7],
                                                  u[8], geom, A_t, A_e, A_s, g, R, ep, es, rec,
                                                  eps, sc_lo, sc_hi, ss_lo, ss_hi, huang);
  };
  for (const auto& c : cases) {
    const auto rj = rj_at(c.u);
    for (int row = 0; row < 4; ++row) {
      for (int k = 0; k < 9; ++k) {
        double a[9], b[9];
        for (int m = 0; m < 9; ++m) a[m] = b[m] = c.u[m];
        const double h = 1e-6 * std::abs(c.u[k]);
        a[k] += h;
        b[k] -= h;
        const double cd = (rj_at(a).residuals[row] - rj_at(b).residuals[row]) / (2 * h);
        const double floor = 1e-6 * std::abs(rj.residuals[row]) + 1e-7;
        EXPECT_NEAR(rj.jacobian[row][k], cd, 1e-5 * std::abs(cd) + floor)
            << "case " << c.id << " row " << row << " seed " << k;
      }
    }
  }
}

TEST(EjectorJacobians, ElementGasPartialsMatchCentralDifference) {
  // d_dgamma and d_dr_gas (#1246): gamma and R are the secondary port gas's,
  // which moves with the flows when a leg reverses, so the element carries
  // its rows' partials in both. A central difference in gamma moves the
  // nozzle-exit Mach M_p1 too (the bisection's root), so this also checks
  // its implicit derivative. Every regime, both anchor families, and a
  // reversed secondary (ms < 0: omega's C1 floor, whose level moves with the
  // choked flow, hence with gamma and R).
  const double A_t = 3.14e-5, A_e = 1.0e-4, A_mix = 8.0e-4, A_s = A_mix - A_e;
  const EjectorGeometry geom{A_e / A_t, A_mix / A_t};
  const double ep = 0.95, es = 0.85, rec = 0.9, eps = 1e-3;
  const double sc_lo = 0.90, sc_hi = 0.999, ss_lo = 0.98, ss_hi = 1.02;
  struct Case {
    const char* id;
    double u[9];
  };
  const Case cases[] = {
      {"jetpump", {0.0051, 0.030, 0.0351, 105000.0, 305.0, 100000.0, 295.0, 100500.0, 99000.0}},
      {"jetpump_reversed", {0.0051, -0.010, -0.0049, 105000.0, 305.0, 100000.0, 295.0, 100500.0,
                            99000.0}},
      {"critical", {0.44, 0.13, 0.57, 700000.0, 440.0, 22850.0, 280.0, 30000.0, 100000.0}},
      {"subcritical", {0.44, 0.20, 0.64, 700000.0, 440.0, 22850.0, 280.0, 52000.0, 60000.0}},
      {"beyond_pb0", {0.44, 0.05, 0.49, 700000.0, 440.0, 22850.0, 280.0, 200000.0, 60000.0}},
      {"choked_reversed", {0.44, -0.07, 0.37, 700000.0, 440.0, 22850.0, 280.0, 40000.0,
                           100000.0}},
  };
  for (const auto closure : {EjectorAnchorClosure::kracik_dvorak, EjectorAnchorClosure::huang1999}) {
    for (const auto& c : cases) {
      auto rj_at = [&](double gamma, double r_gas) {
        return ejector_element_residuals_and_jacobian(
            c.u[0], c.u[1], c.u[2], c.u[3], c.u[4], c.u[5], c.u[6], c.u[7], c.u[8], geom, A_t,
            A_e, A_s, gamma, r_gas, ep, es, rec, eps, sc_lo, sc_hi, ss_lo, ss_hi, closure);
      };
      const double g = 1.37, R = 291.0;
      const auto rj = rj_at(g, R);
      const double hg = 1e-6 * g, hr = 1e-6 * R;
      const auto gp = rj_at(g + hg, R), gm = rj_at(g - hg, R);
      const auto rp = rj_at(g, R + hr), rm = rj_at(g, R - hr);
      for (int row = 0; row < 4; ++row) {
        const double cd_g = (gp.residuals[row] - gm.residuals[row]) / (2 * hg);
        const double cd_r = (rp.residuals[row] - rm.residuals[row]) / (2 * hr);
        // On the scale of the row's own terms (|R_i| and the gamma / R
        // columns times their values), as the element's other checks.
        const double scale = std::abs(rj.residuals[row]) + std::abs(cd_g) * g +
                             std::abs(cd_r) * R + 1e-9;
        EXPECT_NEAR(rj.d_dgamma[row] * g, cd_g * g, 1e-6 * scale)
            << c.id << " closure " << static_cast<int>(closure) << " row " << row << " gamma";
        EXPECT_NEAR(rj.d_dr_gas[row] * R, cd_r * R, 1e-6 * scale)
            << c.id << " closure " << static_cast<int>(closure) << " row " << row << " r_gas";
      }
    }
  }
}
