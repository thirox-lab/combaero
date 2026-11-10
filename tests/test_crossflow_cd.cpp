#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <vector>

#include "composition.h"
#include "crossflow_cd.h"
#include "math_constants.h"
#include "thermo.h"

using combaero::CrossflowCd;
using combaero::florschuetz_metzger_su_1982_cd;
using combaero::florschuetz_metzger_su_1982_has_plate;
using combaero::rogers_hersh_1976_cd;
using combaero::rogers_hersh_1976_envelope;

namespace {

constexpr double kGc = 40.0;  // crossflow mass velocity [kg/(m^2 s)]

// Central differences of a law in (G_ideal, G_c).
template <class F>
void expect_partials(F law, double Gi, double Gc, double rel) {
  const CrossflowCd r = law(Gi, Gc);
  const double hi = Gi * 1e-6;
  const double hc = Gc * 1e-6;
  const double fd_i = (law(Gi + hi, Gc).Cd - law(Gi - hi, Gc).Cd) / (2.0 * hi);
  const double fd_c = (law(Gi, Gc + hc).Cd - law(Gi, Gc - hc).Cd) / (2.0 * hc);
  EXPECT_NEAR(r.dCd_dG_ideal, fd_i, rel * std::abs(fd_i) + 1e-14);
  EXPECT_NEAR(r.dCd_dG_c, fd_c, rel * std::abs(fd_c) + 1e-14);
}

double zeta(double x) { return 1.06 / std::pow(x + 0.806, 0.602); }

}  // namespace

// ---- Rogers & Hersh (1976) --------------------------------------------------

TEST(RogersHersh1976, ResolvesThePapersFitInTheActualJetVelocity) {
  // The paper: Cd = a (Vi/V_inf)^b with Vi the ACTUAL jet velocity, Gj = Cd G_ideal.
  for (const bool suction : {false, true}) {
    const double a = suction ? 0.80 : 0.97;
    const double b = suction ? 0.47 : 0.534;
    for (const double Gi : {0.05, 0.5, 5.0, 20.0}) {
      const CrossflowCd r = rogers_hersh_1976_cd(0.85, Gi, kGc, suction);
      ASSERT_FALSE(r.capped);
      ASSERT_LT(r.Cd, 0.9 * 0.85);  // below the cap's blend band
      const double ratio = r.Cd * Gi / kGc;  // Gj / Gc
      EXPECT_NEAR(r.Cd, a * std::pow(ratio, b), 1e-12);
    }
  }
}

TEST(RogersHersh1976, HoldsTheHolesOwnCdUnderAWeakCrossflow) {
  EXPECT_TRUE(rogers_hersh_1976_cd(0.79, 200.0, kGc, false).capped);
  EXPECT_DOUBLE_EQ(rogers_hersh_1976_cd(0.79, 200.0, kGc, false).Cd, 0.79);
  EXPECT_EQ(rogers_hersh_1976_cd(0.79, 200.0, kGc, false).dCd_dG_c, 0.0);
  // No crossflow at all: the hole's own Cd, whatever the drive.
  EXPECT_TRUE(rogers_hersh_1976_cd(0.79, 1e-6, 0.0, false).capped);
  // Injection reaches the 0.79 cap at Gc/Gj = 1.47: Florschuetz's arrays
  // (Gc/Gj <= 0.77) never leave it.
  const double q = std::pow(0.79 / std::pow(0.97, 1.0 / 0.466), 0.466 / 0.534);
  EXPECT_NEAR(1.0 / (0.79 * q), 1.469, 2e-3);
}

TEST(RogersHersh1976, JoinsTheCapC1) {
  // Across the blend band the Cd and its slope are continuous, and it never
  // exceeds Cd0 or the fit.
  double prev = 0.0;
  double prev_slope = 0.0;
  for (double Gi = 25.0; Gi <= 60.0; Gi += 0.01) {
    const CrossflowCd r = rogers_hersh_1976_cd(0.85, Gi, kGc, false);
    const double fit = std::pow(0.97, 1.0 / 0.466) * std::pow(Gi / kGc, 0.534 / 0.466);
    EXPECT_LE(r.Cd, 0.85 + 1e-15);
    EXPECT_LE(r.Cd, fit + 1e-15);
    if (prev > 0.0) {
      // The secant over the step is the mean of the end slopes (C1, no jump).
      const double secant = (r.Cd - prev) / 0.01;
      EXPECT_NEAR(secant, 0.5 * (r.dCd_dG_ideal + prev_slope), 1e-5) << Gi;
    }
    prev = r.Cd;
    prev_slope = r.dCd_dG_ideal;
  }
}

TEST(RogersHersh1976, NoDriveNoFlow) {
  EXPECT_EQ(rogers_hersh_1976_cd(0.85, 0.0, kGc, false).Cd, 0.0);
  EXPECT_EQ(rogers_hersh_1976_cd(0.85, -1.0, kGc, true).Cd, 0.0);
}

TEST(RogersHersh1976, AtTheStagnantRowsOfCr3630) {
  // Gc/Gj ~ 125 (B(5,4,1) mc/mj 0.99, rows 2-4): Cd ~ 0.07-0.08 either way.
  for (const bool suction : {false, true}) {
    double Gi = 1.0;
    for (int it = 0; it < 60; ++it) {  // find Gi with Gj/Gc = 1/125
      const double cd = rogers_hersh_1976_cd(0.85, Gi, kGc, suction).Cd;
      Gi *= std::sqrt((kGc / 125.0) / (cd * Gi));
    }
    const double cd = rogers_hersh_1976_cd(0.85, Gi, kGc, suction).Cd;
    EXPECT_NEAR(cd, suction ? 0.8 * std::pow(0.008, 0.47) : 0.97 * std::pow(0.008, 0.534), 1e-6);
    EXPECT_LT(cd, 0.09);
  }
}

TEST(RogersHersh1976, AnalyticPartials) {
  for (const bool suction : {false, true}) {
    auto law = [suction](double Gi, double Gc) {
      return rogers_hersh_1976_cd(0.85, Gi, Gc, suction);
    };
    for (const double Gi : {0.02, 0.4, 8.0, 30.0, 38.0, 41.0}) {
      expect_partials(law, Gi, kGc, 1e-6);
    }
  }
}

TEST(RogersHersh1976, EnvelopeNamesEachVariable) {
  // Fig. 11's own conditions: in range.
  EXPECT_TRUE(rogers_hersh_1976_envelope(0.1, 1.0, 50.0).in_range());
  const auto lo = rogers_hersh_1976_envelope(0.001, 1.0, 50.0);
  EXPECT_TRUE(lo.velocity_ratio);
  EXPECT_FALSE(lo.l_over_d);
  EXPECT_FALSE(lo.grazing_velocity);
  const auto thick = rogers_hersh_1976_envelope(0.1, 6.0, 50.0);
  EXPECT_TRUE(thick.l_over_d);
  EXPECT_FALSE(thick.velocity_ratio);
  const auto slow = rogers_hersh_1976_envelope(0.1, 1.0, 10.0);
  EXPECT_TRUE(slow.grazing_velocity);
  EXPECT_FALSE(slow.in_range());
  // The edges are inside.
  EXPECT_TRUE(rogers_hersh_1976_envelope(0.004, 0.12, 24.0).in_range());
  EXPECT_TRUE(rogers_hersh_1976_envelope(0.9, 2.0, 85.0).in_range());
}

TEST(RogersHersh1976, RefusesANonPositiveCd0) {
  EXPECT_THROW(rogers_hersh_1976_cd(0.0, 1.0, kGc, false), std::invalid_argument);
}

// ---- Florschuetz, Metzger & Su (1982), CR-3630 Table 4.1 --------------------

TEST(FlorschuetzMetzgerSu1982, CoversTheTestedPlatesOnly) {
  EXPECT_TRUE(florschuetz_metzger_su_1982_has_plate(5, 4, 1));
  EXPECT_TRUE(florschuetz_metzger_su_1982_has_plate(10, 8, 1));
  EXPECT_FALSE(florschuetz_metzger_su_1982_has_plate(5, 8, 2));
  EXPECT_FALSE(florschuetz_metzger_su_1982_has_plate(15, 4, 1));
  EXPECT_THROW(florschuetz_metzger_su_1982_cd(5, 8, 2, 1.0, kGc), std::invalid_argument);
}

// The law in x = Gc/Gj as implemented: the printed pieces joined
// continuously (see crossflow_cd.h).
double table41(double xn, double yn, double z, double x) {
  if (xn == 5 && yn == 4) {
    const double mult = z == 1 ? 1.0 : (z == 2 ? 1.08 : 1.11);
    return std::min(0.85, mult * zeta(x));
  }
  const bool b581 = xn == 5 && yn == 8;
  const double A = b581 ? -0.169 : -0.128;
  const double end = b581 ? 1.5 : 1.8;
  const double cd0 = b581 ? 0.80 : 0.76;
  const double B = zeta(end) - A * end;
  return x <= end ? std::min(cd0, A * x + B) : zeta(x);
}

TEST(FlorschuetzMetzgerSu1982, IsTheTableInGcOverGj) {
  struct Case {
    double xn, yn, z;
  };
  for (const Case c : {Case{5, 4, 1}, Case{5, 4, 2}, Case{5, 4, 3}, Case{5, 8, 1}, Case{10, 4, 1}}) {
    for (const double Gi : {0.05, 0.3, 2.0, 8.0, 30.0, 45.0, 52.0, 56.0, 60.0, 80.0, 400.0}) {
      const CrossflowCd r = florschuetz_metzger_su_1982_cd(c.xn, c.yn, c.z, Gi, kGc);
      const double x = kGc / (r.Cd * Gi);  // Gc / Gj at the solved Cd
      const double cd0 = c.xn == 10 ? 0.76 : (c.yn == 8 ? 0.80 : 0.85);
      if (r.capped) {
        EXPECT_DOUBLE_EQ(r.Cd, cd0);
      } else if (r.Cd <= 0.9 * cd0) {
        // Below the join band: the table's falling piece exactly.
        EXPECT_NEAR(r.Cd, table41(c.xn, c.yn, c.z, x), 1e-12)
            << c.xn << "," << c.yn << "," << c.z << " Gi " << Gi;
      } else {
        // In the band: between the piece and cd0, within 2.5% of either.
        EXPECT_LE(r.Cd, cd0);
        EXPECT_NEAR(r.Cd, table41(c.xn, c.yn, c.z, x), 0.025 * cd0);
      }
    }
  }
}

TEST(FlorschuetzMetzgerSu1982, StaysWithinThePrintedDigitsOfThePrintedPieces) {
  // The printed Table 4.1, piece by piece, against the continuous law.
  for (double x = 0.0; x <= 6.0; x += 0.01) {
    const double p541 = x <= 0.63 ? 0.85 : zeta(x);
    const double p581 = x <= 0.54 ? 0.80 : (x <= 1.5 ? -0.169 * x + 0.893 : zeta(x));
    const double p1041 = x <= 0.54 ? 0.76 : (x <= 1.8 ? -0.128 * x + 0.825 : zeta(x));
    EXPECT_NEAR(table41(5, 4, 1, x), p541, 3e-3) << x;
    EXPECT_NEAR(table41(5, 8, 1, x), p581, 3e-3) << x;
    EXPECT_NEAR(table41(10, 4, 1, x), p1041, 4e-3) << x;
  }
}

TEST(FlorschuetzMetzgerSu1982, AnalyticPartials) {
  auto law41 = [](double Gi, double Gc) { return florschuetz_metzger_su_1982_cd(5, 4, 1, Gi, Gc); };
  auto law81 = [](double Gi, double Gc) { return florschuetz_metzger_su_1982_cd(5, 8, 1, Gi, Gc); };
  for (const double Gi : {0.05, 1.0, 10.0, 40.0, 52.0, 56.0}) {
    expect_partials(law41, Gi, kGc, 1e-6);
  }
  // On B(5,8,1)'s linear piece (x 0.54-1.5) and on zeta beyond it.
  for (const double Gi : {35.0, 60.0, 5.0}) {
    expect_partials(law81, Gi, kGc, 1e-6);
  }
}

// ---- Rowbury, Oldfield & Lock (2001): exit crossflow -----------------------

using combaero::Rowbury2001Coefficients;
using combaero::Rowbury2001Set;
using combaero::rowbury_2001_coefficients;
using combaero::rowbury_2001_exit_crossflow_cd;
using combaero::rowbury_2001_in_range;
using combaero::rowbury_2001_k_floor;
using combaero::rowbury_2001_k_perp;

namespace {

double sind(double deg) { return std::sin(deg * M_PI / 180.0); }

}  // namespace

TEST(Rowbury2001, CoefficientsAreTable2) {
  // Table 2, p. 263, evaluated by hand at theta 40, alpha 6 (a suction row).
  const double st = sind(40.0);
  const double sa = sind(6.0);
  struct Case {
    Rowbury2001Set set;
    double A, B, C;
  };
  const Case cases[] = {
      {Rowbury2001Set::CylindricalPressureSurface, (2.95 * st - 0.88) * (1 + sa),
       (0.714 * st + 0.393) * (1 + sa), 0.10},
      {Rowbury2001Set::CylindricalSuctionSurface, (0.682 * st - 0.24) * (1 + sa),
       (0.777 * st + 0.405) * (1 + sa), 0.10},
      {Rowbury2001Set::FanShapedPressureSurface, (0.435 * st - 0.140) * (1 + sa),
       (0.350 * st + 0.570) * (1 + sa), 0.02},
      {Rowbury2001Set::FanShapedSuctionSurface, (0.225 * st - 0.110) * (1 + sa),
       (0.350 * st + 0.570) * (1 + 2.8 * sa), 0.05},
  };
  for (const auto& c : cases) {
    const auto k = rowbury_2001_coefficients(c.set, 40.0, 6.0);
    EXPECT_NEAR(k.A, c.A, 1e-14);
    EXPECT_NEAR(k.B, c.B, 1e-14);
    EXPECT_EQ(k.C, c.C);
  }
  // The leading edge at its own rows' angles (theta 84, alpha 60).
  const auto le = rowbury_2001_coefficients(Rowbury2001Set::CylindricalLeadingEdge, 84.0, 60.0);
  EXPECT_NEAR(le.A, 0.079 * sind(84.0) * (1 + sind(60.0)), 1e-14);
  EXPECT_NEAR(le.B, (4.795 * sind(84.0) - 4.487) * (1 + sind(60.0)), 1e-14);
  EXPECT_EQ(le.C, 0.28);
}

TEST(Rowbury2001, RefusesAnglesWhereEq8HasNoMeaning) {
  // Suction cylindrical A crosses zero at sin(theta) = 0.24/0.682 (20.6 deg),
  // the leading edge's B at 4.487/4.795 (69.3 deg).
  EXPECT_THROW(rowbury_2001_coefficients(Rowbury2001Set::CylindricalSuctionSurface, 20.0, 0.0),
               std::invalid_argument);
  EXPECT_NO_THROW(rowbury_2001_coefficients(Rowbury2001Set::CylindricalSuctionSurface, 21.0, 0.0));
  EXPECT_THROW(rowbury_2001_coefficients(Rowbury2001Set::CylindricalLeadingEdge, 60.0, 60.0),
               std::invalid_argument);
  EXPECT_THROW(rowbury_2001_coefficients(Rowbury2001Set::CylindricalSuctionSurface, 0.0, 0.0),
               std::invalid_argument);
  EXPECT_THROW(rowbury_2001_coefficients(Rowbury2001Set::CylindricalSuctionSurface, 30.0, -1.0),
               std::invalid_argument);
}

TEST(Rowbury2001, Eq7AndEq8) {
  const auto c = rowbury_2001_coefficients(Rowbury2001Set::CylindricalSuctionSurface, 30.0, 0.0);
  for (const double k : {0.05, 0.3, combaero::rowbury_2001_k_fit_max(
                                        Rowbury2001Set::CylindricalSuctionSurface)}) {
    const auto r =
        rowbury_2001_exit_crossflow_cd(0.776, k, Rowbury2001Set::CylindricalSuctionSurface, 30.0, 0.0);
    const double delta = c.A * std::pow(k, -c.B) - c.C;
    EXPECT_NEAR(r.delta, delta, 1e-14);
    EXPECT_NEAR(1.0 / (r.Cd * r.Cd) - 1.0 / (0.776 * 0.776), delta, 1e-12);
    EXPECT_FALSE(r.held);
  }
}

TEST(Rowbury2001, CrossoverInsideTheFitAndCd0AtInfiniteK) {
  // Inside the fitted k the -C offset raises Cd above Cd0 (crossover; the
  // leading-edge rows reach it, k_fit 4.3); beyond it the offset fades, and
  // k = inf (no mainstream) is Cd0 exactly.
  const auto le = Rowbury2001Set::CylindricalLeadingEdge;
  const double cd0 = 0.8;
  EXPECT_GT(rowbury_2001_exit_crossflow_cd(cd0, combaero::rowbury_2001_k_fit_max(le), le, 88.0, 60.0).Cd,
            cd0);
  const auto set = Rowbury2001Set::CylindricalSuctionSurface;
  const auto inf =
      rowbury_2001_exit_crossflow_cd(cd0, std::numeric_limits<double>::infinity(), set, 30.0, 0.0);
  EXPECT_EQ(inf.Cd, cd0);
  EXPECT_TRUE(inf.faded);
  // Continuous towards it: |Cd - Cd0| falls with k far out.
  double prev = 1.0;
  for (const double k : {1e2, 1e4, 1e8, 1e16}) {
    const double d = std::abs(rowbury_2001_exit_crossflow_cd(cd0, k, set, 30.0, 0.0).Cd - cd0);
    EXPECT_LT(d, prev);
    prev = d;
  }
  EXPECT_LT(prev, 2e-4);
}

TEST(Rowbury2001, TheFadeJoinsEq8C1AtTheFittedKMax) {
  const auto set = Rowbury2001Set::CylindricalLeadingEdge;
  const double kmax = combaero::rowbury_2001_k_fit_max(set);
  const auto c = rowbury_2001_coefficients(set, 88.0, 60.0);
  const double e = 1e-9 * kmax;
  const auto lo = rowbury_2001_exit_crossflow_cd(0.7, kmax - e, set, 88.0, 60.0);
  const auto hi = rowbury_2001_exit_crossflow_cd(0.7, kmax + e, set, 88.0, 60.0);
  EXPECT_FALSE(lo.faded);
  EXPECT_TRUE(hi.faded);
  EXPECT_NEAR(lo.delta, c.A * std::pow(kmax, -c.B) - c.C, 1e-8);
  EXPECT_NEAR(hi.delta, lo.delta, 1e-8);
  EXPECT_NEAR(hi.dCd_dk, lo.dCd_dk, 1e-6 * std::abs(lo.dCd_dk));
  // Low slope beyond: |d delta / d ln k| <= B A k^-B + 0.65 C / width.
  for (const double k : {4.0, 10.0, 100.0, 1e4}) {
    const auto r = rowbury_2001_exit_crossflow_cd(0.7, k, set, 88.0, 60.0);
    const double ddelta_dlnk = -2.0 * r.dCd_dk / std::pow(r.Cd, 3) * k;
    EXPECT_LE(std::abs(ddelta_dlnk), c.B * c.A * std::pow(k, -c.B) +
                                         0.65 * c.C / combaero::ROWBURY_2001_FADE_WIDTH + 1e-12);
  }
}

TEST(Rowbury2001, ExitDensityIsTheIsentropeAndItsPartialsMatchFD) {
  const std::vector<double> X = {0.7809, 0.2095, 0.0093, 0.0003};
  std::vector<double> Xf(combaero::num_species(), 0.0);
  // dry air in the species order N2, O2, AR, CO2 is not assumed: map by name
  const char* names[] = {"N2", "O2", "AR", "CO2"};
  for (int i = 0; i < 4; ++i) {
    Xf[combaero::species_index_from_name(names[i])] = X[i];
  }
  const double Pt = 2.0e5, Tt = 400.0, P = 1.2e5;
  const auto r = combaero::rowbury_2001_exit_density(Pt, Tt, P, Xf);
  EXPECT_NEAR(combaero::s_mass(r.T, Xf, P), combaero::s_mass(Tt, Xf, Pt), 1e-8);
  auto rho = [&](double pt, double tt, double p) {
    return combaero::rowbury_2001_exit_density(pt, tt, p, Xf).rho;
  };
  EXPECT_NEAR(r.drho_dPt, (rho(Pt + 1, Tt, P) - rho(Pt - 1, Tt, P)) / 2.0, 1e-7 * std::abs(r.drho_dPt));
  EXPECT_NEAR(r.drho_dP, (rho(Pt, Tt, P + 1) - rho(Pt, Tt, P - 1)) / 2.0, 1e-7 * std::abs(r.drho_dP));
  EXPECT_NEAR(r.drho_dTt, (rho(Pt, Tt + 1e-3, P) - rho(Pt, Tt - 1e-3, P)) / 2e-3,
              1e-6 * std::abs(r.drho_dTt));
  // Mass-fraction partials against differences in Y (renormalised to X).
  const std::vector<double> Y = combaero::mole_to_mass(Xf);
  for (std::size_t k = 0; k < Y.size(); ++k) {
    if (Y[k] <= 0.0) continue;
    const double h = 1e-7;
    std::vector<double> Yp = Y, Ym = Y;
    Yp[k] += h;
    Ym[k] -= h;
    const double fd = (combaero::rowbury_2001_exit_density(Pt, Tt, P, combaero::mass_to_mole(Yp)).rho -
                       combaero::rowbury_2001_exit_density(Pt, Tt, P, combaero::mass_to_mole(Ym)).rho) /
                      (2.0 * h);
    EXPECT_NEAR(r.drho_dY[k], fd, 1e-5 * std::abs(fd) + 1e-9) << k;
  }
}

TEST(Rowbury2001, PartialsMatchCentralDifferences) {
  for (const auto set : {Rowbury2001Set::CylindricalSuctionSurface,
                         Rowbury2001Set::CylindricalPressureSurface,
                         Rowbury2001Set::FanShapedSuctionSurface}) {
    for (const double k : {0.4, 1.5, 9.0}) {
      const double cd0 = 0.75;
      const auto r = rowbury_2001_exit_crossflow_cd(cd0, k, set, 40.0, 3.0);
      const double hk = 1e-6 * k;
      const double hc = 1e-6;
      const double fd_k = (rowbury_2001_exit_crossflow_cd(cd0, k + hk, set, 40.0, 3.0).Cd -
                           rowbury_2001_exit_crossflow_cd(cd0, k - hk, set, 40.0, 3.0).Cd) /
                          (2.0 * hk);
      const double fd_c = (rowbury_2001_exit_crossflow_cd(cd0 + hc, k, set, 40.0, 3.0).Cd -
                           rowbury_2001_exit_crossflow_cd(cd0 - hc, k, set, 40.0, 3.0).Cd) /
                          (2.0 * hc);
      EXPECT_NEAR(r.dCd_dk, fd_k, 1e-7 * std::abs(fd_k) + 1e-12);
      EXPECT_NEAR(r.dCd_dCd0, fd_c, 1e-7 * std::abs(fd_c) + 1e-12);
    }
  }
}

TEST(Rowbury2001, HeldBelowTheFloorValueKeptSlopeZero) {
  const auto set = Rowbury2001Set::CylindricalSuctionSurface;
  // B < 1 at 30 deg: the floor is the authors' lowest evaluated k.
  EXPECT_EQ(rowbury_2001_k_floor(set, 30.0, 0.0), combaero::ROWBURY_2001_K_PERP_MIN);
  const double kf = rowbury_2001_k_floor(set, 30.0, 0.0);
  const auto at = rowbury_2001_exit_crossflow_cd(0.8, kf, set, 30.0, 0.0);
  const auto below = rowbury_2001_exit_crossflow_cd(0.8, 0.0, set, 30.0, 0.0);
  EXPECT_FALSE(at.held);
  EXPECT_TRUE(below.held);
  EXPECT_EQ(below.Cd, at.Cd);
  EXPECT_EQ(below.dCd_dk, 0.0);
  EXPECT_GT(below.Cd, 0.0);
  // B > 1 (a normal hole): (B - 1) A k^-B = M (1 - C) at the floor.
  const auto c = rowbury_2001_coefficients(set, 90.0, 0.0);
  ASSERT_GT(c.B, 1.0);
  const double ku = rowbury_2001_k_floor(set, 90.0, 0.0);
  EXPECT_NEAR((c.B - 1.0) * c.A * std::pow(ku, -c.B),
              combaero::ROWBURY_2001_UNIQUE_MARGIN * (1.0 - c.C), 1e-12);
}

TEST(Rowbury2001, FlowRootIsUniqueAboveTheFloor) {
  // m = Cd(m) K with k ~ m^2: d ln Cd / d ln m = B Cd^2 A k^-B must stay
  // below 1 on the unheld branch, for the worst B of every set at Table 1's
  // angles and any Cd0 <= 1.
  for (const auto set :
       {Rowbury2001Set::CylindricalPressureSurface, Rowbury2001Set::CylindricalLeadingEdge,
        Rowbury2001Set::CylindricalSuctionSurface, Rowbury2001Set::FanShapedPressureSurface,
        Rowbury2001Set::FanShapedSuctionSurface}) {
    for (const double th : {35.0, 60.0, 85.0, 90.0}) {
      for (const double al : {0.0, 60.0}) {
        Rowbury2001Coefficients c;
        try {
          c = rowbury_2001_coefficients(set, th, al);
        } catch (const std::invalid_argument&) {
          continue;
        }
        const double kf = rowbury_2001_k_floor(set, th, al);
        for (const double cd0 : {0.5, 0.8, 1.0}) {
          const auto r = rowbury_2001_exit_crossflow_cd(cd0, kf, set, th, al);
          const double slope = c.B * r.Cd * r.Cd * c.A * std::pow(kf, -c.B);
          EXPECT_LT(slope, 1.0) << "theta " << th << " alpha " << al << " Cd0 " << cd0;
        }
      }
    }
  }
}

TEST(Rowbury2001, FitRangeIsTheLargestKTheRowsReach) {
  // k_perp = Cd^2 Ma_j^2 / Ma_m^2 f / AR^2 at the row's design PR (Cd 0.70)
  // or printed sweep, gamma 1.4 (see crossflow_cd.h, "THE TWO ENDS").
  auto k = [](double th, double al, double mm, double pr, double cd, double ar) {
    const double mj2 = 2.0 / 0.4 * (std::pow(pr, 0.4 / 1.4) - 1.0);
    const double ct = th < 90.0 ? 1.0 / std::tan(th * M_PI / 180.0) : 0.0;
    const double ta = std::tan(al * M_PI / 180.0);
    return cd * cd * mj2 / (mm * mm) / std::sqrt(ct * ct + 1.0 + ta * ta) / (ar * ar);
  };
  using combaero::rowbury_2001_k_fit_max;
  EXPECT_NEAR(rowbury_2001_k_fit_max(Rowbury2001Set::CylindricalPressureSurface),
              k(58, 0, 0.05, 1.02, 0.70, 1.0), 0.05);  // row 6
  EXPECT_NEAR(rowbury_2001_k_fit_max(Rowbury2001Set::CylindricalLeadingEdge),
              k(90, 60, 0.04, 1.02, 0.70, 1.0), 0.05);  // row 7
  EXPECT_NEAR(combaero::rowbury_2001_k_data_max(Rowbury2001Set::CylindricalSuctionSurface),
              k(78, 0, 0.32, 1.09, 0.70, 1.0), 0.005);  // row 11
  // Tier 2: the authors' own application, Fig. 9 Ma_m 0.3 to PR 2.0 (the
  // runner's evaluation of their curve; python test_exit_crossflow_cd).
  EXPECT_EQ(rowbury_2001_k_fit_max(Rowbury2001Set::CylindricalSuctionSurface), 3.8);
  for (const auto set : {Rowbury2001Set::CylindricalPressureSurface,
                         Rowbury2001Set::CylindricalLeadingEdge,
                         Rowbury2001Set::FanShapedPressureSurface,
                         Rowbury2001Set::FanShapedSuctionSurface}) {
    EXPECT_EQ(rowbury_2001_k_fit_max(set), combaero::rowbury_2001_k_data_max(set));
  }
  EXPECT_NEAR(rowbury_2001_k_fit_max(Rowbury2001Set::FanShapedPressureSurface),
              k(30, 7, 0.08, 1.02, 0.70, 2.9), 0.005);  // row 4
  EXPECT_NEAR(rowbury_2001_k_fit_max(Rowbury2001Set::FanShapedSuctionSurface),
              k(30, 6, 0.96, 2.5, 0.75, 2.9), 0.0005);  // row 14, Fig. 8
}

TEST(Rowbury2001, InRangeIsTable1) {
  EXPECT_TRUE(rowbury_2001_in_range(Rowbury2001Set::CylindricalSuctionSurface, 30.0, 0.0));
  EXPECT_TRUE(rowbury_2001_in_range(Rowbury2001Set::CylindricalSuctionSurface, 78.0, 6.0));
  EXPECT_FALSE(rowbury_2001_in_range(Rowbury2001Set::CylindricalSuctionSurface, 90.0, 0.0));
  EXPECT_TRUE(rowbury_2001_in_range(Rowbury2001Set::CylindricalLeadingEdge, 89.0, 60.0));
  EXPECT_FALSE(rowbury_2001_in_range(Rowbury2001Set::CylindricalLeadingEdge, 89.0, 0.0));
  EXPECT_FALSE(rowbury_2001_in_range(Rowbury2001Set::FanShapedSuctionSurface, 50.0, 0.0));
}

TEST(Rowbury2001, KPerpAndItsPartials) {
  const double G = 120.0, rho = 2.0, q = 3.0e4, th = 35.0;
  const auto r = rowbury_2001_k_perp(G, rho, q, th, 0.0);
  EXPECT_NEAR(r.k, G * G * sind(th) / (rho * q), 1e-14);
  const double h = 1e-6;
  EXPECT_NEAR(r.dk_dG,
              (rowbury_2001_k_perp(G + h, rho, q, th, 0.0).k -
               rowbury_2001_k_perp(G - h, rho, q, th, 0.0).k) /
                  (2 * h),
              1e-7 * std::abs(r.dk_dG));
  EXPECT_NEAR(r.dk_drho_h, -r.k / rho, 1e-14);
  EXPECT_NEAR(r.dk_dq_m, -r.k / q, 1e-18);
  EXPECT_THROW(rowbury_2001_k_perp(G, rho, 0.0, th, 0.0), std::invalid_argument);
  EXPECT_THROW(rowbury_2001_k_perp(G, rho, q, th, 90.0), std::invalid_argument);
}

TEST(Rowbury2001, KPerpTakesTheSurfaceNormalComponent) {
  // Fig. 3: axis along (cot theta, 1, tan alpha), alpha from the normal.
  // A normal-ish showerhead hole (theta 90, alpha 60) has u_perp = u / 2.
  const double G = 100.0, rho = 1.0, q = 1.0e4;
  EXPECT_NEAR(rowbury_2001_k_perp(G, rho, q, 90.0, 60.0).k,
              0.5 * rowbury_2001_k_perp(G, rho, q, 90.0, 0.0).k, 1e-12);
  const double th = 40.0, al = 6.0;
  const double ct = 1.0 / std::tan(th * M_PI / 180.0), ta = std::tan(al * M_PI / 180.0);
  EXPECT_NEAR(rowbury_2001_k_perp(G, rho, q, th, al).k,
              G * G / std::sqrt(ct * ct + 1.0 + ta * ta) / (rho * q), 1e-12);
}
