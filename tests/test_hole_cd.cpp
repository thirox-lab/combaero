#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <tuple>
#include <vector>

#include "composition.h"
#include "crossflow_cd.h"
#include "hole_cd.h"
#include "humidair.h"
#include "orifice.h"
#include "solver_interface.h"

// hole_cd_compose: each term against the kernel it composes (exact equality:
// the chain must not re-derive anything), and every dCd/d(input) against
// central differences of the composed Cd.

using combaero::HOLE_CD_FIXED_MIN;
using combaero::hole_cd_compose;
using combaero::HoleCdBase;
using combaero::HoleCdInputs;
using combaero::HoleCdJetCrossflow;
using combaero::HoleCdModel;
using combaero::HoleCdResult;
using combaero::Rowbury2001Set;

namespace {

DischargeHoleGeometry hole(double inclination_deg = 90.0) {
  DischargeHoleGeometry h;
  h.d = 1e-3;
  h.L = 2e-3;
  h.r = 0.0;
  h.inclination_deg = inclination_deg;
  return h;
}

OrificeGeometry pipe_orifice() {
  OrificeGeometry g;
  g.d = 0.05;
  g.D = 0.1;
  return g;
}

HoleCdModel discharge(DischargeCdCorrelation id) {
  HoleCdModel m;
  m.base = HoleCdBase::Discharge;
  m.discharge = id;
  return m;
}

HoleCdModel with_exit(HoleCdModel m) {
  m.exit_crossflow = true;
  m.exit_set = Rowbury2001Set::CylindricalSuctionSurface;
  m.exit_theta_deg = 30.0;
  m.exit_alpha_deg = 0.0;
  return m;
}

HoleCdInputs state() {
  HoleCdInputs in;
  in.Re = 2e4;
  in.U1_over_Vi = 0.2;
  in.Re_D = 1e5;
  in.G_h = 200.0;      // k_perp = G^2 sin(30) / (rho q) = 0.25: not held
  in.rho_h = 4.0;
  in.rho_m_u_m2 = 2e4;
  in.G_ideal = 10.0;   // q = 0.25: Rogers & Hersh acts
  in.G_c = 40.0;
  return in;
}

// Central differences of the composed Cd in each input, against the
// analytic partials. Relative step 1e-6; tolerance 1e-6 relative.
void expect_partials(const HoleCdModel& model, const DischargeHoleGeometry& h,
                     const HoleCdInputs& in) {
  const HoleCdResult r = hole_cd_compose(model, h, pipe_orifice(), in);
  struct Arg {
    double HoleCdInputs::*field;
    double analytic;
    const char* name;
  };
  const Arg args[] = {
      {&HoleCdInputs::Re, r.dCd_dRe, "Re"},
      {&HoleCdInputs::U1_over_Vi, r.dCd_dU1_over_Vi, "U1_over_Vi"},
      {&HoleCdInputs::Re_D, r.dCd_dRe_D, "Re_D"},
      {&HoleCdInputs::G_h, r.dCd_dG_h, "G_h"},
      {&HoleCdInputs::rho_h, r.dCd_drho_h, "rho_h"},
      {&HoleCdInputs::rho_m_u_m2, r.dCd_drho_m_u_m2, "rho_m_u_m2"},
      {&HoleCdInputs::G_ideal, r.dCd_dG_ideal, "G_ideal"},
      {&HoleCdInputs::G_c, r.dCd_dG_c, "G_c"},
  };
  for (const Arg& a : args) {
    const double x = in.*(a.field);
    const double step = 1e-6 * (std::abs(x) > 0.0 ? std::abs(x) : 1.0);
    HoleCdInputs hi = in;
    HoleCdInputs lo = in;
    hi.*(a.field) = x + step;
    lo.*(a.field) = x - step;
    const double fd = (hole_cd_compose(model, h, pipe_orifice(), hi).Cd -
                       hole_cd_compose(model, h, pipe_orifice(), lo).Cd) /
                      (2.0 * step);
    EXPECT_NEAR(a.analytic, fd, 1e-6 * std::abs(fd) + 1e-12) << a.name;
  }
}

}  // namespace

// ---- term 1: the base Cd ------------------------------------------------------

TEST(HoleCdCompose, FixedIsTheCallersValueClampedWithNoPartials) {
  HoleCdModel m;
  m.fixed_cd = 0.7;
  EXPECT_EQ(hole_cd_compose(m, hole(), pipe_orifice(), state()).Cd, 0.7);
  m.fixed_cd = 1.3;
  EXPECT_EQ(hole_cd_compose(m, hole(), pipe_orifice(), state()).Cd, 1.0);
  m.fixed_cd = 0.0;
  const HoleCdResult r = hole_cd_compose(m, hole(), pipe_orifice(), state());
  EXPECT_EQ(r.Cd, HOLE_CD_FIXED_MIN);
  EXPECT_EQ(r.dCd_dRe, 0.0);
  EXPECT_EQ(r.dCd_dU1_over_Vi, 0.0);
  EXPECT_EQ(r.dCd_dRe_D, 0.0);
}

TEST(HoleCdCompose, MeteringIsTheIsoCorrelationAtThePipeReynoldsNumber) {
  HoleCdModel m;
  m.base = HoleCdBase::Metering;
  const auto [cd, dRe] = metering_cd_and_dRe(m.metering, pipe_orifice(), 1e5);
  const HoleCdResult r = hole_cd_compose(m, hole(), pipe_orifice(), state());
  EXPECT_EQ(r.Cd, cd);
  EXPECT_EQ(r.dCd_dRe_D, dRe);
  EXPECT_EQ(r.dCd_dRe, 0.0);
  expect_partials(m, hole(), state());
}

TEST(HoleCdCompose, DischargeIsTheCorrelationsValueAndPartials) {
  for (const DischargeCdCorrelation id :
       {DischargeCdCorrelation::McGreehanSchotsch1988, DischargeCdCorrelation::Idelchik1966Thick,
        DischargeCdCorrelation::Lichtarowicz1965}) {
    const auto c = make_discharge_correlation(id);
    DischargeHoleState flow;
    flow.Re = 2e4;
    flow.U1_over_Vi = 0.2;
    const auto [cd, dRe, dU] = c->Cd_and_derivatives(hole(), flow);
    (void)cd;
    const HoleCdResult r = hole_cd_compose(discharge(id), hole(), pipe_orifice(), state());
    EXPECT_EQ(r.Cd, c->Cd(hole(), flow));
    EXPECT_EQ(r.cd_base, r.Cd);
    EXPECT_EQ(r.dCd_dRe, dRe);
    EXPECT_EQ(r.dCd_dU1_over_Vi, dU);
    EXPECT_FALSE(r.re_held);
    expect_partials(discharge(id), hole(), state());
  }
}

TEST(HoleCdCompose, InclinationIsInsideTheHoleGeometry) {
  const HoleCdModel m = discharge(DischargeCdCorrelation::Idelchik1966Thick);
  const DischargeHoleGeometry h = hole(30.0);
  const auto c = make_discharge_correlation(m.discharge);
  DischargeHoleState flow;
  flow.Re = 2e4;
  flow.U1_over_Vi = 0.2;
  const HoleCdResult r = hole_cd_compose(m, h, pipe_orifice(), state());
  EXPECT_EQ(r.Cd, c->Cd(h, flow));
  EXPECT_EQ(r.inclination_dzeta, discharge_inclination_dzeta(30.0));
  EXPECT_LT(r.Cd, hole_cd_compose(m, hole(), pipe_orifice(), state()).Cd);
  EXPECT_EQ(hole_cd_compose(m, hole(), pipe_orifice(), state()).inclination_dzeta, 0.0);
  expect_partials(m, h, state());
}

TEST(HoleCdCompose, BelowTheReFloorTheHeldValueHasNoSlopeUnlessContinued) {
  HoleCdModel m = discharge(DischargeCdCorrelation::McGreehanSchotsch1988);
  HoleCdInputs in = state();
  in.Re = 0.5 * discharge_cd_re_floor(m.discharge);
  const auto c = make_discharge_correlation(m.discharge);
  DischargeHoleState flow;
  flow.Re = in.Re;
  flow.U1_over_Vi = in.U1_over_Vi;
  const double continued = std::get<1>(c->Cd_and_derivatives(hole(), flow));
  ASSERT_NE(continued, 0.0);
  const HoleCdResult held = hole_cd_compose(m, hole(), pipe_orifice(), in);
  EXPECT_TRUE(held.re_held);
  EXPECT_EQ(held.dCd_dRe, 0.0);
  m.continue_re_floor = true;
  EXPECT_EQ(hole_cd_compose(m, hole(), pipe_orifice(), in).dCd_dRe, continued);
}

// ---- term 2: the exit crossflow (Rowbury et al. 2001) ---------------------------

TEST(HoleCdCompose, ExitTermIsRowburyOnTheBaseCd) {
  const HoleCdModel m = with_exit(discharge(DischargeCdCorrelation::McGreehanSchotsch1988));
  const HoleCdInputs in = state();
  const HoleCdResult r = hole_cd_compose(m, hole(), pipe_orifice(), in);
  const auto kp = combaero::rowbury_2001_k_perp(in.G_h, in.rho_h, in.rho_m_u_m2, 30.0, 0.0);
  const auto ex =
      combaero::rowbury_2001_exit_crossflow_cd(r.cd_base, kp.k, m.exit_set, 30.0, 0.0);
  ASSERT_TRUE(r.exit_active);
  ASSERT_FALSE(ex.held);
  EXPECT_EQ(r.k_perp.k, kp.k);
  EXPECT_EQ(r.Cd, ex.Cd);
  EXPECT_EQ(r.cd_hole, ex.Cd);
  EXPECT_LT(r.Cd, r.cd_base);
  EXPECT_EQ(r.exit_ratio, ex.Cd / r.cd_base);
  EXPECT_EQ(r.exit_k_floor, combaero::rowbury_2001_k_floor(m.exit_set, 30.0, 0.0));
  EXPECT_EQ(r.exit_in_range, combaero::rowbury_2001_in_range(m.exit_set, 30.0, 0.0));
  EXPECT_EQ(r.dCd_dRe, ex.dCd_dCd0 * r.dcd_base_dRe);
  EXPECT_EQ(r.dCd_dG_h, ex.dCd_dk * kp.dk_dG);
  expect_partials(m, hole(), in);
}

TEST(HoleCdCompose, ExitTermOnAFixedAndAMeteringBase) {
  HoleCdModel fixed = with_exit(HoleCdModel{});
  fixed.fixed_cd = 0.8;
  expect_partials(fixed, hole(), state());
  HoleCdModel metering = with_exit(HoleCdModel{});
  metering.base = HoleCdBase::Metering;
  expect_partials(metering, hole(), state());
}

TEST(HoleCdCompose, ExitTermIsOffWithoutMainstreamMomentum) {
  const HoleCdModel m = with_exit(discharge(DischargeCdCorrelation::McGreehanSchotsch1988));
  HoleCdInputs in = state();
  in.rho_m_u_m2 = 0.0;
  const HoleCdResult r = hole_cd_compose(m, hole(), pipe_orifice(), in);
  EXPECT_FALSE(r.exit_active);
  EXPECT_EQ(r.Cd, r.cd_base);
  EXPECT_EQ(r.exit_ratio, 1.0);
  EXPECT_EQ(r.dCd_dG_h, 0.0);
  EXPECT_EQ(r.exit_k_tier, 0.0);
}

TEST(HoleCdCompose, HeldExitTermHasNoKPartials) {
  const HoleCdModel m = with_exit(discharge(DischargeCdCorrelation::McGreehanSchotsch1988));
  HoleCdInputs in = state();
  in.G_h = 1.0;  // k_perp << floor
  const HoleCdResult r = hole_cd_compose(m, hole(), pipe_orifice(), in);
  ASSERT_TRUE(r.exit.held);
  EXPECT_EQ(r.dCd_dG_h, 0.0);
  EXPECT_EQ(r.dCd_drho_h, 0.0);
  EXPECT_EQ(r.dCd_drho_m_u_m2, 0.0);
  EXPECT_EQ(r.dCd_dRe, r.exit.dCd_dCd0 * r.dcd_base_dRe);
}

// ---- term 3: the impingement jet-crossflow law ----------------------------------

TEST(HoleCdCompose, RogersHershReplacesTheHolesCdWithItsOwnPartials) {
  HoleCdModel m = discharge(DischargeCdCorrelation::Idelchik1966Thick);
  m.jet = HoleCdJetCrossflow::RogersHersh1976;
  const HoleCdInputs in = state();
  const HoleCdResult r = hole_cd_compose(m, hole(), pipe_orifice(), in);
  const auto law = combaero::rogers_hersh_1976_cd(r.cd_hole, in.G_ideal, in.G_c, false);
  ASSERT_FALSE(law.capped);
  EXPECT_TRUE(r.jet_acts);
  EXPECT_TRUE(r.jet_partials);
  EXPECT_EQ(r.Cd, law.Cd);
  EXPECT_EQ(r.dCd_dG_ideal, law.dCd_dG_ideal);
  EXPECT_EQ(r.dCd_dG_c, law.dCd_dG_c);
  // The hole's Cd is the law's cap only: no partials through it.
  EXPECT_EQ(r.dCd_dRe, 0.0);
  EXPECT_EQ(r.dCd_dU1_over_Vi, 0.0);
  HoleCdInputs suction = in;
  suction.suction = true;
  EXPECT_EQ(hole_cd_compose(m, hole(), pipe_orifice(), suction).Cd,
            combaero::rogers_hersh_1976_cd(r.cd_hole, in.G_ideal, in.G_c, true).Cd);
}

TEST(HoleCdCompose, CappedRogersHershHoldsTheHolesCdAndPartials) {
  HoleCdModel m = discharge(DischargeCdCorrelation::McGreehanSchotsch1988);
  m.jet = HoleCdJetCrossflow::RogersHersh1976;
  HoleCdInputs in = state();
  in.G_ideal = 200.0;
  const HoleCdResult r = hole_cd_compose(m, hole(), pipe_orifice(), in);
  ASSERT_TRUE(r.jet_capped);
  EXPECT_FALSE(r.jet_acts);
  EXPECT_EQ(r.Cd, r.cd_hole);
  EXPECT_EQ(r.dCd_dRe, r.dcd_base_dRe);
  EXPECT_EQ(r.dCd_dG_ideal, 0.0);
}

TEST(HoleCdCompose, NoDriveLeavesTheLawWithoutPartials) {
  HoleCdModel m = discharge(DischargeCdCorrelation::Idelchik1966Thick);
  m.jet = HoleCdJetCrossflow::RogersHersh1976;
  HoleCdInputs in = state();
  in.G_ideal = 0.0;
  const HoleCdResult r = hole_cd_compose(m, hole(), pipe_orifice(), in);
  EXPECT_TRUE(r.jet_acts);
  EXPECT_FALSE(r.jet_partials);
  EXPECT_EQ(r.dCd_dG_c, 0.0);
  EXPECT_EQ(r.dCd_dRe, 0.0);
}

TEST(HoleCdCompose, FlorschuetzMetzgerSuIsTheTablesOwnCd) {
  HoleCdModel m = discharge(DischargeCdCorrelation::Idelchik1966Thick);
  m.jet = HoleCdJetCrossflow::FlorschuetzMetzgerSu1982;
  m.xn_d = 5.0;
  m.yn_d = 4.0;
  m.z_d = 1.0;
  const HoleCdInputs in = state();
  const auto law =
      combaero::florschuetz_metzger_su_1982_cd(5.0, 4.0, 1.0, in.G_ideal, in.G_c);
  const HoleCdResult r = hole_cd_compose(m, hole(), pipe_orifice(), in);
  EXPECT_TRUE(r.jet_acts);
  EXPECT_EQ(r.Cd, law.Cd);
  EXPECT_EQ(r.jet_capped, law.capped);
  expect_partials(m, hole(), in);
}

TEST(HoleCdCompose, EveryPartialMatchesCentralDifferences) {
  for (const bool exit : {false, true}) {
    for (const HoleCdJetCrossflow jet :
         {HoleCdJetCrossflow::Off, HoleCdJetCrossflow::RogersHersh1976}) {
      HoleCdModel m = discharge(DischargeCdCorrelation::McGreehanSchotsch1988);
      if (exit) {
        m = with_exit(m);
      }
      m.jet = jet;
      expect_partials(m, hole(), state());
    }
  }
}

// hole_supply_drive (#1087): Static is the station static; Total adds the
// momentum-chamber closure's dynamic head, with analytic partials checked
// against central differences (steps small against the head, #1107).
TEST(HoleSupplyDrive, StaticIsTheStationStaticAndTotalAddsTheClosuresHead) {
  using combaero::hole_supply_drive;
  using combaero::HoleSupplyDrive;
  const std::vector<double> Y = combaero::mole_to_mass(combaero::dry_air());
  const double P = 1.2e5, Tt = 450.0, m = 0.35, A = 1.2e-3;
  for (const bool inc : {false, true}) {
    const auto s = hole_supply_drive(HoleSupplyDrive::Static, P, Tt, Y, m, A, inc);
    const auto t = hole_supply_drive(HoleSupplyDrive::Total, P, Tt, Y, m, A, inc);
    const auto c =
        combaero::solver::momentum_chamber_residual_and_jacobian(P, 0.0, m, Tt, Y, A, inc);
    EXPECT_EQ(s.P_drive, P);
    EXPECT_EQ(s.dP_dP, 1.0);
    EXPECT_NEAR(s.q_dyn, c.q_dyn, 1e-12 * c.q_dyn);  // -residual - P: round-off apart
    EXPECT_GT(t.q_dyn, 0.0);
    EXPECT_EQ(t.P_drive, -c.residual);
    EXPECT_NEAR(t.P_drive, P + c.q_dyn, 1e-12 * t.P_drive);
    auto f = [&](double P_, double T_, double m_) {
      return hole_supply_drive(HoleSupplyDrive::Total, P_, T_, Y, m_, A, inc).P_drive;
    };
    const double hP = 1e-4 * t.q_dyn, hT = 1e-3, hm = 1e-6 * m;
    EXPECT_NEAR((f(P + hP, Tt, m) - f(P - hP, Tt, m)) / (2 * hP), t.dP_dP, 1e-7 * std::abs(t.dP_dP));
    EXPECT_NEAR((f(P, Tt + hT, m) - f(P, Tt - hT, m)) / (2 * hT), t.dP_dTt,
                1e-6 * std::abs(t.dP_dTt));
    EXPECT_NEAR((f(P, Tt, m + hm) - f(P, Tt, m - hm)) / (2 * hm), t.dP_dm_dot,
                1e-7 * std::abs(t.dP_dm_dot));
    // Even in the flow.
    const auto r = hole_supply_drive(HoleSupplyDrive::Total, P, Tt, Y, -m, A, inc);
    EXPECT_EQ(r.P_drive, t.P_drive);
    EXPECT_EQ(r.dP_dm_dot, -t.dP_dm_dot);
  }
}

// Past the station's sonic static the closure is its bounded continuation
// (#788): the drive's value and its partials must still be one function.
TEST(HoleSupplyDrive, BoundedBranchPartialsMatchTheValue) {
  using combaero::hole_supply_drive;
  using combaero::HoleSupplyDrive;
  const std::vector<double> Y = combaero::mole_to_mass(combaero::dry_air());
  const double P = 1.2e5, Tt = 450.0, A = 1.2e-3;
  for (const double m : {0.35, 0.8, 1.2}) {
    const auto t = hole_supply_drive(HoleSupplyDrive::Total, P, Tt, Y, m, A, false);
    const auto c = combaero::solver::momentum_chamber_residual_and_jacobian(P, 0.0, m, Tt, Y, A,
                                                                            false);
    EXPECT_EQ(t.P_drive, -c.residual);
    EXPECT_EQ(t.q_dyn, t.P_drive - P);
    auto f = [&](double P_, double T_, double m_) {
      return hole_supply_drive(HoleSupplyDrive::Total, P_, T_, Y, m_, A, false).P_drive;
    };
    const double hP = 1e-3 * std::abs(t.q_dyn), hT = 1e-4, hm = 1e-7 * m;
    EXPECT_NEAR((f(P + hP, Tt, m) - f(P - hP, Tt, m)) / (2 * hP), t.dP_dP,
                1e-5 * std::abs(t.dP_dP) + 1e-9);
    EXPECT_NEAR((f(P, Tt + hT, m) - f(P, Tt - hT, m)) / (2 * hT), t.dP_dTt,
                1e-5 * std::abs(t.dP_dTt) + 1e-6);
    EXPECT_NEAR((f(P, Tt, m + hm) - f(P, Tt, m - hm)) / (2 * hm), t.dP_dm_dot,
                1e-5 * std::abs(t.dP_dm_dot) + 1e-3);
  }
}

// The static-to-total Cd conversion: the ratio of the kernel's own ideal
// flows, 0 where P1 <= P2, partials against central differences.
TEST(HoleStaticToTotalCdRatio, IsTheKernelsIdealFlowRatio) {
  using combaero::hole_static_to_total_cd_ratio;
  const std::vector<double> Y = combaero::mole_to_mass(combaero::dry_air());
  const double P1 = 1.2e5, T1 = 430.0, Pt1 = 1.32e5, Tt1 = 450.0, P2 = 1.05e5;
  const auto r = hole_static_to_total_cd_ratio(P1, T1, Pt1, Tt1, P2, Y);
  const auto s = combaero::solver::orifice_compressible_residuals_and_jacobian(0.0, P1, P1, T1, Y,
                                                                               P2, 1.0, 1.0, 0.0);
  const auto t = combaero::solver::orifice_compressible_residuals_and_jacobian(
      0.0, Pt1, Pt1, Tt1, Y, P2, 1.0, 1.0, 0.0);
  EXPECT_DOUBLE_EQ(r.ratio, s.m_dot_calc / t.m_dot_calc);
  EXPECT_EQ(hole_static_to_total_cd_ratio(P2, T1, Pt1, Tt1, P2, Y).ratio, 0.0);
  auto f = [&](double a, double b, double c, double d, double e) {
    return hole_static_to_total_cd_ratio(a, b, c, d, e, Y).ratio;
  };
  const double hp = 1.0, ht = 1e-3;
  EXPECT_NEAR((f(P1 + hp, T1, Pt1, Tt1, P2) - f(P1 - hp, T1, Pt1, Tt1, P2)) / (2 * hp), r.d_dP1,
              1e-6 * std::abs(r.d_dP1));
  EXPECT_NEAR((f(P1, T1 + ht, Pt1, Tt1, P2) - f(P1, T1 - ht, Pt1, Tt1, P2)) / (2 * ht), r.d_dT1,
              1e-6 * std::abs(r.d_dT1));
  EXPECT_NEAR((f(P1, T1, Pt1 + hp, Tt1, P2) - f(P1, T1, Pt1 - hp, Tt1, P2)) / (2 * hp), r.d_dPt1,
              1e-6 * std::abs(r.d_dPt1));
  EXPECT_NEAR((f(P1, T1, Pt1, Tt1 + ht, P2) - f(P1, T1, Pt1, Tt1 - ht, P2)) / (2 * ht), r.d_dTt1,
              1e-6 * std::abs(r.d_dTt1));
  EXPECT_NEAR((f(P1, T1, Pt1, Tt1, P2 + hp) - f(P1, T1, Pt1, Tt1, P2 - hp)) / (2 * hp), r.d_dP2,
              1e-6 * std::abs(r.d_dP2));
}
