#include <gtest/gtest.h>

#include "../include/composition.h"
#include "../include/fanno_mach.h"
#include "../include/humidair.h"
#include "../include/solver_interface.h"
#include "../include/thermo.h"
#include "../include/transport.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <tuple>
#include <vector>

using namespace combaero;

// ---------------------------------------------------------------------------
// #1216: Compressible ducts at creeping flow recover Hagen-Poiseuille (to the
// laminar blend's 1 - w(0)) whenever the incompressible floor-drop fraction
// chi = dp_inc(at.G) / kFannoFlowDriveFloor >= kFannoLowMachDropHi (satisfied
// by every unchoked real-gas duct), with continuous value and Jacobian handoff
// across kFannoFlowDriveFloor away from pre-existing compressible branch flips
// (#1307).
// ---------------------------------------------------------------------------

namespace {

double hp_flux(double Pt0, double Tt0, const std::vector<double>& X, double dp, double L, double D,
               const std::string& model) {
  const double rho = density(Tt0, Pt0, X);
  const double mu = viscosity(Tt0, Pt0, X);
  const double f_re_0 = solver::friction_f_re_at_zero(model);
  return 2.0 * rho * D * D * dp / (f_re_0 * mu * L);
}

}  // namespace

TEST(FannoCreepingFlow1216, FrictionFReAndJacobianExactAndSmooth) {
  const std::vector<std::string> models = {"haaland", "serghides", "colebrook", "petukhov"};
  for (const std::string& m : models) {
    const double e_D = (m == "petukhov") ? 0.0 : 1e-4;
    const auto [f0, df0] = solver::friction_f_re_and_jacobian(m, 0.0, e_D);
    EXPECT_DOUBLE_EQ(f0, solver::friction_f_re_at_zero(m)) << m;
    EXPECT_TRUE(std::isfinite(df0)) << m;

    // One-sided FD at Re = 0
    const double h0 = 1e-4;
    const double f_p = std::get<0>(solver::friction_f_re_and_jacobian(m, h0, e_D));
    EXPECT_NEAR(df0, (f_p - f0) / h0, 1e-6 * std::max(1.0, std::abs(df0))) << m;

    // Central FD across laminar, transitional and turbulent regimes
    // (step >= 1e-4 so ulp(64) / (2 h) stays below 1e-10).
    for (const double Re : {1e-2, 1.0, 10.0, 2300.0, 1e5}) {
      const auto [val, slope] = solver::friction_f_re_and_jacobian(m, Re, e_D);
      const double h = std::max(1e-5 * Re, 1e-4);
      const double vp = std::get<0>(solver::friction_f_re_and_jacobian(m, Re + h, e_D));
      const double vm = std::get<0>(solver::friction_f_re_and_jacobian(m, Re - h, e_D));
      const double fd = (vp - vm) / (2.0 * h);
      EXPECT_NEAR(slope, fd, 1e-6 * std::max(1.0, std::abs(fd))) << m << " Re=" << Re;
      if (Re >= 10.0) {
        const double f_kern = std::get<0>(solver::friction_and_jacobian(m, Re, e_D).result);
        EXPECT_NEAR(val / (Re * f_kern), 1.0, 1e-11) << m << " Re=" << Re;
      }
    }
  }
}

TEST(FannoCreepingFlow1216, HagenPoiseuilleLimitAcrossLengthRatiosAndModels) {
  const std::vector<double> X = dry_air();
  const double Pt0 = 1e5;
  const double Tt0 = 300.0;
  const double D = 0.01;
  const double P_target = Pt0 - 1e-9;
  const double dp_actual = Pt0 - P_target;

  const std::vector<std::string> models = {"haaland", "serghides", "colebrook", "petukhov"};
  for (const std::string& model : models) {
    for (const double L_D : {1.0, 10.0, 100.0, 1000.0}) {
      const double L = L_D * D;
      const double g_hp = hp_flux(Pt0, Tt0, X, dp_actual, L, D, model);
      const FannoChannelFlow f_stat =
          fanno_channel_flow(Pt0, Tt0, X, P_target, false, L, D, 0.0, model);
      const FannoChannelFlow f_tot =
          fanno_channel_flow(Pt0, Tt0, X, P_target, true, L, D, 0.0, model);
      const FannoChannelFlow r_stat =
          fanno_rayleigh_channel_flow(Pt0, Tt0, Tt0, X, P_target, false, L, D, 0.0, model);
      const FannoChannelFlow r_tot =
          fanno_rayleigh_channel_flow(Pt0, Tt0, Tt0, X, P_target, true, L, D, 0.0, model);

      EXPECT_NEAR(f_stat.G / g_hp, 1.0, 1e-5) << model << " L/D=" << L_D;
      EXPECT_NEAR(f_tot.G / g_hp, 1.0, 1e-5) << model << " L/D=" << L_D;
      EXPECT_NEAR(r_stat.G / g_hp, 1.0, 1e-5) << model << " L/D=" << L_D;
      EXPECT_NEAR(r_tot.G / g_hp, 1.0, 1e-5) << model << " L/D=" << L_D;
      EXPECT_NEAR(f_stat.G / f_tot.G, 1.0, 1e-5) << model << " L/D=" << L_D;
    }
  }
}

TEST(FannoCreepingFlow1216, HagenPoiseuilleLimitAcrossPressuresAndTemperatures) {
  const std::vector<double> X = dry_air();
  const double D = 0.01;
  for (const double Pt0 : {1e4, 1e5, 1e6, 1e7}) {
    for (const double Tt0 : {300.0, 1000.0, 2000.0}) {
      for (const double L_D : {1.0, 100.0}) {
        const double L = L_D * D;
        // Scale dp small enough that dynamic head q = G_hp^2 / (2 rho) << dp,
        // while staying well above ulp(Pt0) so Pt0 - P_target is non-zero.
        const double rho = density(Tt0, Pt0, X);
        const double mu = viscosity(Tt0, Pt0, X);
        const double dp_inert =
            (64.0 * mu * L) * (64.0 * mu * L) / (2.0 * rho * D * D * D * D);
        const double dp_req = std::max(std::min(1e-6, 1e-5 * dp_inert), 1e-12 * Pt0);
        const double P_target = Pt0 - dp_req;
        const double dp_actual = Pt0 - P_target;
        const double g_hp = hp_flux(Pt0, Tt0, X, dp_actual, L, D, "haaland");
        const FannoChannelFlow fs =
            fanno_channel_flow(Pt0, Tt0, X, P_target, false, L, D, 0.0, "haaland");
        const FannoChannelFlow ft =
            fanno_channel_flow(Pt0, Tt0, X, P_target, true, L, D, 0.0, "haaland");
        // For exit_total=false, physical exit head reduces G by O(dp / dp_inert).
        const double tol_s = std::max(1e-4, 2.0 * dp_actual / dp_inert);
        EXPECT_NEAR(fs.G / g_hp, 1.0, tol_s)
            << "Pt0=" << Pt0 << " Tt0=" << Tt0 << " L/D=" << L_D;
        EXPECT_NEAR(ft.G / g_hp, 1.0, 1e-4)
            << "Pt0=" << Pt0 << " Tt0=" << Tt0 << " L/D=" << L_D;
      }
    }
  }
}

TEST(FannoCreepingFlow1216, DiabaticCreepingFlowMatchesVariablePropertyHagenPoiseuille) {
  const std::vector<double> X = dry_air();
  const double Pt0 = 1e5;
  const double D = 0.01;
  const double L = 0.1;
  const double P_target = Pt0 - 1e-8;
  const double dp_actual = Pt0 - P_target;
  const double R = specific_gas_constant(X);
  const double f_re_0 = solver::friction_f_re_at_zero("haaland");

  for (const auto [Tt0, Tt_exit] :
       std::vector<std::pair<double, double>>{{300.0, 900.0}, {900.0, 400.0}}) {
    // Midpoint reference integral of mu(Tt(xi)) * Tt(xi) on [0, 1]
    constexpr int N = 2000;
    double integral = 0.0;
    for (int i = 0; i < N; ++i) {
      const double xi = (i + 0.5) / N;
      const double T = (1.0 - xi) * Tt0 + xi * Tt_exit;
      integral += viscosity(T, Pt0, X) * T;
    }
    integral /= N;
    const double g_ref = 2.0 * Pt0 * D * D * dp_actual / (f_re_0 * R * L * integral);
    for (const bool exit_total : {false, true}) {
      const FannoChannelFlow r =
          fanno_rayleigh_channel_flow(Pt0, Tt0, Tt_exit, X, P_target, exit_total, L, D, 0.0,
                                      "haaland");
      EXPECT_NEAR(r.G / g_ref, 1.0, 1e-5)
          << "Tt0=" << Tt0 << " Tt_exit=" << Tt_exit << " exit_total=" << exit_total;
    }
  }
}

TEST(FannoCreepingFlow1216, ValueAndAllJacobiansContinuousAcrossDriveFloor) {
  const std::vector<double> X = dry_air();
  std::vector<double> v_o2(num_species(), 0.0);
  v_o2[species_index_from_name("O2")] = 1.0;
  const std::vector<std::vector<double>> dirs = {v_o2};

  const double Pt0 = 2.0e5;
  const double D = 0.01;
  const double floor = kFannoFlowDriveFloor;
  const double eps = 1e-7 * floor;

  // Adiabatic fanno_channel_flow across L/D = 1, 10, 100 and exit_total
  for (const double L_D : {1.0, 10.0, 100.0}) {
    const double L = L_D * D;
    for (const bool exit_total : {false, true}) {
      const FannoChannelFlow lo =
          fanno_channel_flow(Pt0, 400.0, X, Pt0 - (floor - eps), exit_total, L, D, 0.0, "haaland",
                             1.0, -1.0, dirs, true);
      const FannoChannelFlow at =
          fanno_channel_flow(Pt0, 400.0, X, Pt0 - floor, exit_total, L, D, 0.0, "haaland", 1.0,
                             -1.0, dirs, true);
      const FannoChannelFlow hi =
          fanno_channel_flow(Pt0, 400.0, X, Pt0 - (floor + eps), exit_total, L, D, 0.0, "haaland",
                             1.0, -1.0, dirs, true);
      EXPECT_NEAR(lo.G / hi.G, 1.0, 1e-6) << "adiabatic L/D=" << L_D << " et=" << exit_total;
      EXPECT_NEAR(lo.dG_dP_target / hi.dG_dP_target, 1.0, 1e-5)
          << "adiabatic L/D=" << L_D << " et=" << exit_total;
      EXPECT_NEAR(lo.dG_dPt0 / hi.dG_dPt0, 1.0, 1e-5)
          << "adiabatic L/D=" << L_D << " et=" << exit_total;
      EXPECT_NEAR(lo.dG_dTt0 / hi.dG_dTt0, 1.0, 1e-5)
          << "adiabatic L/D=" << L_D << " et=" << exit_total;
      // Compare sub-floor limit (floor - eps) directly against the compressible
      // branch at drive == floor (raw_flow's d_ddir at 1 Pa carries ~1e-2
      // root-rounding jitter between 1.0 and 1.0 + 1e-7 Pa, #1019).
      EXPECT_NEAR(lo.dG_ddir[0] / at.dG_ddir[0], 1.0, 1e-6)
          << "adiabatic L/D=" << L_D << " et=" << exit_total;
    }
  }

  // Diabatic fanno_rayleigh_channel_flow (heated, adiabatic, cooled)
  for (const auto [Tt0, Tt_exit] : std::vector<std::pair<double, double>>{
           {400.0, 400.0}, {400.0, 600.0}, {600.0, 450.0}}) {
    for (const bool exit_total : {false, true}) {
      const double L = 0.2;
      const FannoChannelFlow lo =
          fanno_rayleigh_channel_flow(Pt0, Tt0, Tt_exit, X, Pt0 - (floor - eps), exit_total, L, D,
                                      0.0, "haaland", 1.0, -1.0, dirs, true);
      const FannoChannelFlow at =
          fanno_rayleigh_channel_flow(Pt0, Tt0, Tt_exit, X, Pt0 - floor, exit_total, L, D, 0.0,
                                      "haaland", 1.0, -1.0, dirs, true);
      const FannoChannelFlow hi =
          fanno_rayleigh_channel_flow(Pt0, Tt0, Tt_exit, X, Pt0 - (floor + eps), exit_total, L, D,
                                      0.0, "haaland", 1.0, -1.0, dirs, true);
      EXPECT_NEAR(lo.G / hi.G, 1.0, 1e-6)
          << "rayleigh " << Tt0 << "->" << Tt_exit << " et=" << exit_total;
      EXPECT_NEAR(lo.dG_dP_target / hi.dG_dP_target, 1.0, 1e-5)
          << "rayleigh " << Tt0 << "->" << Tt_exit << " et=" << exit_total;
      EXPECT_NEAR(lo.dG_dPt0 / hi.dG_dPt0, 1.0, 1e-5)
          << "rayleigh " << Tt0 << "->" << Tt_exit << " et=" << exit_total;
      EXPECT_NEAR(lo.dG_dTt0 / hi.dG_dTt0, 1.0, 1e-4)
          << "rayleigh " << Tt0 << "->" << Tt_exit << " et=" << exit_total;
      EXPECT_NEAR(lo.dG_dTt_exit / hi.dG_dTt_exit, 1.0, 1e-4)
          << "rayleigh " << Tt0 << "->" << Tt_exit << " et=" << exit_total;
      EXPECT_NEAR(lo.dG_ddir[0] / at.dG_ddir[0], 1.0, 1e-6)
          << "rayleigh " << Tt0 << "->" << Tt_exit << " et=" << exit_total;
    }
  }
}

TEST(FannoCreepingFlow1216, SubFloorPartialsMatchFiniteDifferences) {
  const std::vector<double> X = dry_air();
  const std::vector<double> Y = mole_to_mass(X);
  const std::size_t k_o2 = species_index_from_name("O2");
  std::vector<double> v_o2(num_species(), 0.0);
  v_o2[k_o2] = 1.0;
  const std::vector<std::vector<double>> dirs = {v_o2};

  const double Pt0 = 2.0e4;
  const double Tt0 = 400.0;
  const double Tt_exit = 550.0;
  const double D = 0.01;

  for (const double L : {0.01, 0.5}) {
    for (const bool heated : {false, true}) {
      for (const bool exit_total : {false, true}) {
        auto eval = [&](double p0, double t0, double te, const std::vector<double>& xk,
                        double pb, bool derivs) {
          if (heated) {
            return fanno_rayleigh_channel_flow(p0, t0, te, xk, pb, exit_total, L, D, 0.0,
                                               "haaland", 1.0, -1.0, derivs ? dirs :
                                               std::vector<std::vector<double>>{}, derivs);
          }
          return fanno_channel_flow(p0, t0, xk, pb, exit_total, L, D, 0.0, "haaland", 1.0, -1.0,
                                    derivs ? dirs : std::vector<std::vector<double>>{}, derivs);
        };
        for (const double drive : {0.0, 1e-4, 0.05, 0.5, 0.9}) {
          const double Pb = Pt0 - drive;
          const FannoChannelFlow r = eval(Pt0, Tt0, Tt_exit, X, Pb, true);
          const double hp = (drive == 0.0) ? 1e-5 : std::min(1e-4 * drive, 1e-4 * (1.0 - drive));
          if (drive == 0.0) {
            // Second-order one-sided FD so O(hp / dp_inert) curvature cancels.
            const double g1_P = eval(Pt0, Tt0, Tt_exit, X, Pt0 - hp, false).G;
            const double g2_P = eval(Pt0, Tt0, Tt_exit, X, Pt0 - 2.0 * hp, false).G;
            const double fd_P = (4.0 * g1_P - g2_P) / (2.0 * hp);
            const double g1_Pt = eval(Pt0 + hp, Tt0, Tt_exit, X, Pb, false).G;
            const double g2_Pt = eval(Pt0 + 2.0 * hp, Tt0, Tt_exit, X, Pb, false).G;
            const double fd_Pt = (4.0 * g1_Pt - g2_Pt) / (2.0 * hp);
            EXPECT_NEAR(-r.dG_dP_target / fd_P, 1.0, 1e-6)
                << "L=" << L << " heated=" << heated << " et=" << exit_total;
            EXPECT_NEAR(r.dG_dPt0 / fd_Pt, 1.0, 1e-6)
                << "L=" << L << " heated=" << heated << " et=" << exit_total;
            EXPECT_EQ(r.dG_dTt0, 0.0);
            EXPECT_EQ(r.dG_dTt_exit, 0.0);
            EXPECT_EQ(r.dG_ddir[0], 0.0);
          } else {
            const double fd_P = (eval(Pt0, Tt0, Tt_exit, X, Pb + hp, false).G -
                                 eval(Pt0, Tt0, Tt_exit, X, Pb - hp, false).G) /
                                (2.0 * hp);
            const double fd_Pt = (eval(Pt0 + hp, Tt0, Tt_exit, X, Pb, false).G -
                                  eval(Pt0 - hp, Tt0, Tt_exit, X, Pb, false).G) /
                                 (2.0 * hp);
            const double hT = 1e-4 * Tt0;
            const double fd_T0 = (eval(Pt0, Tt0 + hT, Tt_exit, X, Pb, false).G -
                                  eval(Pt0, Tt0 - hT, Tt_exit, X, Pb, false).G) /
                                 (2.0 * hT);
            EXPECT_NEAR(r.dG_dP_target / fd_P, 1.0, 1e-4)
                << "L=" << L << " drive=" << drive << " heated=" << heated << " et=" << exit_total;
            EXPECT_NEAR(r.dG_dPt0 / fd_Pt, 1.0, 5e-4)
                << "L=" << L << " drive=" << drive << " heated=" << heated << " et=" << exit_total;
            EXPECT_NEAR(r.dG_dTt0 / fd_T0, 1.0, 5e-4)
                << "L=" << L << " drive=" << drive << " heated=" << heated << " et=" << exit_total;
            if (heated) {
              const double fd_Te = (eval(Pt0, Tt0, Tt_exit + hT, X, Pb, false).G -
                                    eval(Pt0, Tt0, Tt_exit - hT, X, Pb, false).G) /
                                   (2.0 * hT);
              EXPECT_NEAR(r.dG_dTt_exit / fd_Te, 1.0, 5e-4)
                  << "L=" << L << " drive=" << drive << " et=" << exit_total;
            }
            // Step 1e-3 sits cleanly on the FD plateau (1e-4 sits in round-off
            // on heated exit_total ducts under MSVC where ratio dipped to 0.987).
            const double hy = 1e-3;
            std::vector<double> Yp = Y, Ym = Y;
            Yp[k_o2] += hy;
            Ym[k_o2] -= hy;
            const double fd_Y = (eval(Pt0, Tt0, Tt_exit, mass_to_mole(Yp), Pb, false).G -
                                 eval(Pt0, Tt0, Tt_exit, mass_to_mole(Ym), Pb, false).G) /
                                (2.0 * hy);
            EXPECT_NEAR(r.dG_ddir[0] / fd_Y, 1.0, 1e-2)
                << "L=" << L << " drive=" << drive << " heated=" << heated << " et=" << exit_total;
          }
        }
      }
    }
  }
}

TEST(FannoCreepingFlow1216, FixedFrictionKeepsOddCubicFallbackAndContinuity) {
  const std::vector<double> X = dry_air();
  const double Pt0 = 2.0e5;
  const double floor = kFannoFlowDriveFloor;
  const double eps = 1e-7 * floor;
  const FannoChannelFlow lo =
      fanno_channel_flow(Pt0, 400.0, X, Pt0 - (floor - eps), false, 0.5, 0.01, 0.0, "fixed", 0.02);
  const FannoChannelFlow hi =
      fanno_channel_flow(Pt0, 400.0, X, Pt0 - (floor + eps), false, 0.5, 0.01, 0.0, "fixed", 0.02);
  EXPECT_NEAR(lo.G / hi.G, 1.0, 1e-6);
  EXPECT_NEAR(lo.dG_dP_target / hi.dG_dP_target, 1.0, 1e-5);
  const FannoChannelFlow z =
      fanno_channel_flow(Pt0, 400.0, X, Pt0, false, 0.5, 0.01, 0.0, "fixed", 0.02);
  EXPECT_EQ(z.G, 0.0);
  EXPECT_GT(-z.dG_dP_target, 0.0);
  EXPECT_TRUE(std::isfinite(z.dG_dP_target));
}

TEST(FannoCreepingFlow1216, RayleighSubFloorPartialsAtEqualTt) {
  // Independent Tt0 and Tt_exit partials of fanno_rayleigh_channel_flow must
  // match finite differences and remain continuous at Tt_exit == Tt0 (where
  // fanno_channel_flow ties them together).
  const std::vector<double> X = dry_air();
  const double Pt0 = 2.0e5;
  const double Tt0 = 400.0;
  const double L = 0.05;
  const double D = 0.01;
  const double drive = 0.3;
  const double P_target = Pt0 - drive;
  const double hT = 1e-4 * Tt0;

  for (const bool exit_total : {false, true}) {
    const FannoChannelFlow r0 = fanno_rayleigh_channel_flow(
        Pt0, Tt0, Tt0, X, P_target, exit_total, L, D, 0.0, "haaland");
    const FannoChannelFlow f0 =
        fanno_channel_flow(Pt0, Tt0, X, P_target, exit_total, L, D, 0.0, "haaland");
    EXPECT_NEAR(r0.dG_dTt0 + r0.dG_dTt_exit, f0.dG_dTt0,
                1e-5 * std::max(std::abs(f0.dG_dTt0), 1e-9))
        << "exit_total=" << exit_total;
    const double fd_T0 =
        (fanno_rayleigh_channel_flow(Pt0, Tt0 + hT, Tt0, X, P_target, exit_total, L, D, 0.0,
                                     "haaland", 1.0, -1.0, {}, false)
             .G -
         fanno_rayleigh_channel_flow(Pt0, Tt0 - hT, Tt0, X, P_target, exit_total, L, D, 0.0,
                                     "haaland", 1.0, -1.0, {}, false)
             .G) /
        (2.0 * hT);
    const double fd_Te =
        (fanno_rayleigh_channel_flow(Pt0, Tt0, Tt0 + hT, X, P_target, exit_total, L, D, 0.0,
                                     "haaland", 1.0, -1.0, {}, false)
             .G -
         fanno_rayleigh_channel_flow(Pt0, Tt0, Tt0 - hT, X, P_target, exit_total, L, D, 0.0,
                                     "haaland", 1.0, -1.0, {}, false)
             .G) /
        (2.0 * hT);
    const double dT_scale =
        std::max({std::abs(r0.dG_dTt0), std::abs(r0.dG_dTt_exit), 1e-9});
    EXPECT_NEAR(r0.dG_dTt0, fd_T0, 1e-3 * dT_scale) << "exit_total=" << exit_total;
    EXPECT_NEAR(r0.dG_dTt_exit, fd_Te, 1e-3 * dT_scale) << "exit_total=" << exit_total;

    // Continuity of both temperature partials across Tt_exit = Tt0 +- eps
    // (including across c_acc = 0 when exit_total = true). Scale by the
    // channel's temperature-derivative scale since dG_dTt0 sits near a zero
    // crossing at exit_total = false.
    const double eps = 1e-5;
    const FannoChannelFlow r_lo = fanno_rayleigh_channel_flow(
        Pt0, Tt0, Tt0 - eps, X, P_target, exit_total, L, D, 0.0, "haaland");
    const FannoChannelFlow r_hi = fanno_rayleigh_channel_flow(
        Pt0, Tt0, Tt0 + eps, X, P_target, exit_total, L, D, 0.0, "haaland");
    EXPECT_NEAR(r0.dG_dTt0, r_lo.dG_dTt0, 1e-5 * dT_scale) << "exit_total=" << exit_total;
    EXPECT_NEAR(r0.dG_dTt0, r_hi.dG_dTt0, 1e-5 * dT_scale) << "exit_total=" << exit_total;
    EXPECT_NEAR(r0.dG_dTt_exit, r_lo.dG_dTt_exit, 1e-5 * dT_scale)
        << "exit_total=" << exit_total;
    EXPECT_NEAR(r0.dG_dTt_exit, r_hi.dG_dTt_exit, 1e-5 * dT_scale)
        << "exit_total=" << exit_total;
  }
}

TEST(FannoCreepingFlow1216, NoDiscontinuityAcrossCooledOrLowPressureStates) {
  const std::vector<double> X = dry_air();

  // Reproducer A: cooled exit_total duct (L/D = 1) stepping across 389.148 K
  // must stay smooth (no guard flip discontinuity).
  {
    const double Pt0 = 1.0e5;
    const double Tt0 = 400.0;
    const double L = 0.01;
    const double D = 0.01;
    const double P_target = Pt0 - 0.5;
    const double Tte_lo = 389.148;
    const double Tte_hi = 389.149;
    const FannoChannelFlow r_lo =
        fanno_rayleigh_channel_flow(Pt0, Tt0, Tte_lo, X, P_target, true, L, D, 0.0, "haaland");
    const FannoChannelFlow r_hi =
        fanno_rayleigh_channel_flow(Pt0, Tt0, Tte_hi, X, P_target, true, L, D, 0.0, "haaland");
    EXPECT_GT(r_lo.G, 0.0);
    EXPECT_GT(r_hi.G, 0.0);
    EXPECT_FALSE(r_lo.choked);
    EXPECT_FALSE(r_hi.choked);
    EXPECT_NEAR(r_lo.G / r_hi.G, 1.0, 5e-3);
  }

  // Reproducer B: low-pressure adiabatic duct (Pt0 ~ 10 Pa) must stay smooth
  // with no threshold jump where M_ex(1 Pa) ~ 0.05.
  {
    const double Tt0 = 300.0;
    const double L = 0.01;
    const double D = 0.01;
    const double Pt_lo = 10.000;
    const double Pt_hi = 10.006;
    const FannoChannelFlow r_lo =
        fanno_channel_flow(Pt_lo, Tt0, X, Pt_lo - 0.5, false, L, D, 0.0, "haaland");
    const FannoChannelFlow r_hi =
        fanno_channel_flow(Pt_hi, Tt0, X, Pt_hi - 0.5, false, L, D, 0.0, "haaland");
    EXPECT_NEAR(r_lo.G / r_hi.G, 1.0, 2e-3);
  }

  // Reproducer C: cooled exit_total duct (L/D = 10) near 271.9 K must stay
  // strictly positive and smooth.
  {
    const double Pt0 = 1.0e5;
    const double Tt0 = 400.0;
    const double L = 0.10;
    const double D = 0.01;
    const double P_target = Pt0 - 0.5;
    const FannoChannelFlow r_a =
        fanno_rayleigh_channel_flow(Pt0, Tt0, 272.0, X, P_target, true, L, D, 0.0, "haaland");
    const FannoChannelFlow r_b =
        fanno_rayleigh_channel_flow(Pt0, Tt0, 271.9, X, P_target, true, L, D, 0.0, "haaland");
    EXPECT_GT(r_a.G, 0.0);
    EXPECT_GT(r_b.G, 0.0);
    EXPECT_FALSE(r_a.choked);
    EXPECT_FALSE(r_b.choked);
  }
}

TEST(FannoCreepingFlow1216, CooledExitTotalCreepingFlowMatchesHagenPoiseuille) {
  // Strongly cooled exit_total ducts (1500 -> 900 K, L/D = 1 and 10) whose
  // 1 Pa floor state chokes still recover exact variable-property
  // Hagen-Poiseuille as drive -> 0, with choked = false below the floor.
  const std::vector<double> X = dry_air();
  const double Pt0 = 1.0e5;
  const double Tt0 = 1500.0;
  const double Tt_exit = 900.0;
  const double D = 0.01;
  const double P_target = Pt0 - 1e-9;
  const double dp_actual = Pt0 - P_target;
  const double R = specific_gas_constant(X);
  const double f_re_0 = solver::friction_f_re_at_zero("haaland");

  // Midpoint rule (1000 panels) for exact integral_0^1 mu(Tt(xi)) * Tt(xi) dxi.
  double integral = 0.0;
  constexpr int kN = 1000;
  for (int i = 0; i < kN; ++i) {
    const double xi = (static_cast<double>(i) + 0.5) / kN;
    const double Tt = (1.0 - xi) * Tt0 + xi * Tt_exit;
    integral += viscosity(Tt, Pt0, X) * Tt / kN;
  }

  for (const double L_D : {1.0, 10.0}) {
    const double L = L_D * D;
    const double g_hp = (2.0 * Pt0 * D * D * dp_actual) / (f_re_0 * L * R * integral);
    const FannoChannelFlow r = fanno_rayleigh_channel_flow(
        Pt0, Tt0, Tt_exit, X, P_target, true, L, D, 0.0, "haaland");
    EXPECT_FALSE(r.choked) << "L/D=" << L_D;
    EXPECT_NEAR(r.G / g_hp, 1.0, 1e-5) << "L/D=" << L_D;
  }
}

TEST(FannoCreepingFlow1216, StrictlyMonotoneAndBoundedForShortAndChokedDucts) {
  // Vanishingly short ducts (L = 1e-9 m, exit_total = true) and strongly
  // cooled exit_total ducts (1500 -> 900 K) where low_1.G >> at.G must remain
  // strictly increasing and bounded by at.G across [0, floor].
  const std::vector<double> X = dry_air();
  const double Pt0 = 1.0e5;
  const double floor = kFannoFlowDriveFloor;
  const double drives[] = {0.0, 1e-9, 1e-6, 1e-4, 1e-2, 0.1, 0.3, 0.5, 0.8, 0.99, floor};

  // Case 1: vanishingly short adiabatic exit_total duct (L = 1e-9 m, D = 0.028 m)
  {
    const FannoChannelFlow at =
        fanno_channel_flow(Pt0, 300.0, X, Pt0 - floor, true, 1e-9, 0.028, 0.0, "haaland");
    double prev_G = -1.0;
    for (const double d : drives) {
      const FannoChannelFlow r =
          fanno_channel_flow(Pt0, 300.0, X, Pt0 - d, true, 1e-9, 0.028, 0.0, "haaland");
      EXPECT_GT(r.G, prev_G) << "short duct d=" << d;
      EXPECT_LE(r.G, at.G * (1.0 + 1e-12)) << "short duct d=" << d;
      if (d < floor) {
        EXPECT_GT(-r.dG_dP_target, 0.0) << "short duct d=" << d;
      }
      prev_G = r.G;
    }
  }

  // Case 2: strongly cooled exit_total duct (1500 -> 900 K, L/D = 1)
  {
    const FannoChannelFlow at = fanno_rayleigh_channel_flow(
        Pt0, 1500.0, 900.0, X, Pt0 - floor, true, 0.01, 0.01, 0.0, "haaland");
    double prev_G = -1.0;
    for (const double d : drives) {
      const FannoChannelFlow r = fanno_rayleigh_channel_flow(
          Pt0, 1500.0, 900.0, X, Pt0 - d, true, 0.01, 0.01, 0.0, "haaland");
      EXPECT_GT(r.G, prev_G) << "cooled duct d=" << d;
      EXPECT_LE(r.G, at.G * (1.0 + 1e-12)) << "cooled duct d=" << d;
      if (d < floor) {
        EXPECT_GT(-r.dG_dP_target, 0.0) << "cooled duct d=" << d;
      }
      prev_G = r.G;
    }
  }

  // Case 3: smoothstep window chi = dp_inc(at.G) / kFannoFlowDriveFloor in
  // [kFannoLowMachDropLo = 0.125, kFannoLowMachDropHi = 0.25].
  // Covered in full by SonicDropSmoothstepWindowBoundAndContinuous below.
}

TEST(FannoCreepingFlow1216, CooledSubFloorPartialsMatchFiniteDifferences) {
  // Cooled exit_total Rayleigh ducts have large d(sigma)/dtheta at the 1 Pa
  // floor (e.g. Pt0 = 1e5, 400 -> 330 K, L/D = 10 has sigma ~ 1.52 and
  // dsigma/dTt0 ~ +1.55e-2, dsigma/dTt_exit ~ -3.09e-2). Dropping ds/dtheta
  // causes an 89% wrong-sign error in dG/dTt0 at drive = 0.5 Pa.
  const std::vector<double> X = dry_air();
  const std::vector<double> Y = mole_to_mass(X);
  const std::size_t k_o2 = species_index_from_name("O2");
  std::vector<double> v_o2(num_species(), 0.0);
  v_o2[k_o2] = 1.0;
  const std::vector<std::vector<double>> dirs = {v_o2};

  struct Case {
    double Pt0;
    double Tt0;
    double Tt_exit;
    double L;
    double D;
    double drive;
  };
  const std::vector<Case> cases = {
      {1.0e5, 400.0, 330.0, 0.10, 0.01, 0.3},
      {1.0e5, 400.0, 330.0, 0.10, 0.01, 0.5},
      {1.0e5, 400.0, 330.0, 0.10, 0.01, 0.7},
      {1.0e5, 400.0, 300.0, 0.10, 0.01, 0.5},
      {1.0e5, 1500.0, 900.0, 0.01, 0.01, 0.5},
  };

  for (const Case& c : cases) {
    const double Pb = c.Pt0 - c.drive;
    auto eval_G = [&](double p0, double t0, double te, const std::vector<double>& xk, double pb) {
      return fanno_rayleigh_channel_flow(p0, t0, te, xk, pb, true, c.L, c.D, 0.0, "haaland", 1.0,
                                         -1.0, {}, false)
          .G;
    };
    const FannoChannelFlow r = fanno_rayleigh_channel_flow(
        c.Pt0, c.Tt0, c.Tt_exit, X, Pb, true, c.L, c.D, 0.0, "haaland", 1.0, -1.0, dirs, true);

    // Resolved Richardson-extrapolated central FDs (4*D(h/2) - D(h)) / 3
    auto richardson = [](auto&& f, double h) {
      const double d1 = (f(h) - f(-h)) / (2.0 * h);
      const double d2 = (f(0.5 * h) - f(-0.5 * h)) / h;
      return (4.0 * d2 - d1) / 3.0;
    };

    const double hp = 2e-4;
    const double fd_P =
        richardson([&](double d) { return eval_G(c.Pt0, c.Tt0, c.Tt_exit, X, Pb + d); }, hp);
    const double fd_Pt =
        richardson([&](double d) { return eval_G(c.Pt0 + d, c.Tt0, c.Tt_exit, X, Pb); }, hp);
    // Use hT = 0.2 K so both hT and hT/2 = 0.1 K sit on the FD plateau
    // (0.05 K enters rayleigh_raw_flow root round-off near dG/dTt0 ~ 0).
    const double hT = 0.2;
    const double fd_T0 =
        richardson([&](double d) { return eval_G(c.Pt0, c.Tt0 + d, c.Tt_exit, X, Pb); }, hT);
    const double fd_Te =
        richardson([&](double d) { return eval_G(c.Pt0, c.Tt0, c.Tt_exit + d, X, Pb); }, hT);
    const double hy = 2e-3;
    const double fd_Y = richardson(
        [&](double d) {
          std::vector<double> Yk = Y;
          Yk[k_o2] += d;
          return eval_G(c.Pt0, c.Tt0, c.Tt_exit, mass_to_mole(Yk), Pb);
        },
        hy);

    EXPECT_NEAR(r.dG_dP_target / fd_P, 1.0, 1e-4)
        << "T=" << c.Tt0 << "->" << c.Tt_exit << " drive=" << c.drive;
    EXPECT_NEAR(r.dG_dPt0 / fd_Pt, 1.0, 5e-4)
        << "T=" << c.Tt0 << "->" << c.Tt_exit << " drive=" << c.drive;
    EXPECT_NEAR(r.dG_dTt0 / fd_T0, 1.0, 1e-2)
        << "T=" << c.Tt0 << "->" << c.Tt_exit << " drive=" << c.drive;
    EXPECT_NEAR(r.dG_dTt_exit / fd_Te, 1.0, 5e-3)
        << "T=" << c.Tt0 << "->" << c.Tt_exit << " drive=" << c.drive;
    EXPECT_NEAR(r.dG_ddir[0] / fd_Y, 1.0, 1e-2)
        << "T=" << c.Tt0 << "->" << c.Tt_exit << " drive=" << c.drive;
  }
}

TEST(FannoCreepingFlow1216, SonicDropSmoothstepWindowBoundAndContinuous) {
  // Verifies the gas-dynamics-derived sub-floor choking smoothstep window
  //   chi = dp_inc(at.G) / kFannoFlowDriveFloor in [kFannoLowMachDropLo=0.125, kFannoLowMachDropHi=0.25]
  // (below the exact Fanno sonic limit (2/(gamma+1))^((gamma+1)/(gamma-1)) = 0.3349
  // for gamma = 1.4) fails if the window is removed, moved lower/higher, or replaced
  // by a hard switch.
  const std::vector<double> X = dry_air();
  const double Pt0 = 1.0e5;
  const double Tt0 = 300.0;
  const double D = 0.028;
  const double floor = kFannoFlowDriveFloor;
  const double drive = 0.4;
  const double t = drive / floor;

  // Helper: pure odd cubic value and zero-drive slope at (Pt0, Tt0, L, D).
  auto odd_cubic_ref = [&](double L, double d_eval) {
    const FannoChannelFlow at =
        fanno_channel_flow(Pt0, Tt0, X, Pt0 - floor, true, L, D, 0.0, "haaland");
    const double sigma = std::clamp(-at.dG_dP_target * floor / at.G, 0.0, kDriveFloorSigmaMax);
    const DriveFloorCubic c = drive_floor_cubic(d_eval / floor, sigma);
    return std::make_pair(at.G * c.h, at.G * c.dh / floor);
  };

  // 1. Below chi = 0.125 (L = 1e-9 m has chi ~ 1.2e-5; L = 1e-5 m has chi ~ 0.118):
  //    must equal the pure odd cubic (fails if window is removed or moved lower).
  {
    const double P_target = Pt0 - drive;
    const double dp_actual = Pt0 - P_target;
    for (const double L_lo : {1.0e-9, 1.0e-5}) {
      const auto [G_cub, dG0_cub] = odd_cubic_ref(L_lo, dp_actual);
      const auto [G_cub0, slope0_cub] = odd_cubic_ref(L_lo, 0.0);
      (void)G_cub0;
      const FannoChannelFlow r =
          fanno_channel_flow(Pt0, Tt0, X, P_target, true, L_lo, D, 0.0, "haaland");
      const FannoChannelFlow r0 =
          fanno_channel_flow(Pt0, Tt0, X, Pt0, true, L_lo, D, 0.0, "haaland");
      EXPECT_NEAR(r.G / G_cub, 1.0, 1e-12) << "L=" << L_lo;
      EXPECT_NEAR(-r0.dG_dP_target / slope0_cub, 1.0, 1e-12) << "L=" << L_lo;
      EXPECT_GT(dG0_cub, 0.0);
    }
  }

  // 2. Inside chi in (0.125, 0.25) (L = 1.5e-5 m has chi ~ 0.177):
  //    G lies strictly between pure odd cubic and pure viscous rational cubic
  //    (fails if replaced by a hard switch or moved out of [0.125, 0.25]),
  //    and all partials including dw/dchi * dchi/dtheta match FD to 1e-4.
  {
    const double L_mid = 1.5e-5;
    const double Pb = Pt0 - drive;
    const FannoChannelFlow r =
        fanno_channel_flow(Pt0, Tt0, X, Pb, true, L_mid, D, 0.0, "haaland");
    const auto [G_cub, unused_s] = odd_cubic_ref(L_mid, Pt0 - Pb);
    (void)unused_s;
    const FannoChannelFlow at =
        fanno_channel_flow(Pt0, Tt0, X, Pt0 - floor, true, L_mid, D, 0.0, "haaland");
    EXPECT_GT(r.G, G_cub + 0.05 * (at.G - G_cub));
    EXPECT_LT(r.G, at.G - 0.05 * (at.G - G_cub));

    const double hp = 1e-4;
    const double fd_P =
        (fanno_channel_flow(Pt0, Tt0, X, Pb + hp, true, L_mid, D, 0.0, "haaland", 1.0, -1.0, {},
                            false)
             .G -
         fanno_channel_flow(Pt0, Tt0, X, Pb - hp, true, L_mid, D, 0.0, "haaland", 1.0, -1.0, {},
                            false)
             .G) /
        (2.0 * hp);
    const double fd_Pt =
        (fanno_channel_flow(Pt0 + hp, Tt0, X, Pb, true, L_mid, D, 0.0, "haaland", 1.0, -1.0, {},
                            false)
             .G -
         fanno_channel_flow(Pt0 - hp, Tt0, X, Pb, true, L_mid, D, 0.0, "haaland", 1.0, -1.0, {},
                            false)
             .G) /
        (2.0 * hp);
    const double hT = 1e-3 * Tt0;
    const double fd_T0 =
        (fanno_channel_flow(Pt0, Tt0 + hT, X, Pb, true, L_mid, D, 0.0, "haaland", 1.0, -1.0, {},
                            false)
             .G -
         fanno_channel_flow(Pt0, Tt0 - hT, X, Pb, true, L_mid, D, 0.0, "haaland", 1.0, -1.0, {},
                            false)
             .G) /
        (2.0 * hT);
    EXPECT_NEAR(r.dG_dP_target / fd_P, 1.0, 1e-4);
    EXPECT_NEAR(r.dG_dPt0 / fd_Pt, 1.0, 5e-4);
    EXPECT_NEAR(r.dG_dTt0 / fd_T0, 1.0, 1e-4);
  }

  // 3. Continuity across chi = 0.125 and chi = 0.25: stepping L in 2% increments
  //    from L = 1e-5 m (chi ~ 0.118) to 3e-5 m (chi ~ 0.35) produces no jump
  //    larger than 5% of at.G (a hard switch jumps by ~35% of at.G in one step).
  {
    const FannoChannelFlow at_ref =
        fanno_channel_flow(Pt0, Tt0, X, Pt0 - floor, true, 1.5e-5, D, 0.0, "haaland");
    double prev_G =
        fanno_channel_flow(Pt0, Tt0, X, Pt0 - drive, true, 1.0e-5, D, 0.0, "haaland").G;
    for (int k = 1; k <= 56; ++k) {
      const double L_k = 1.0e-5 * std::pow(1.02, k);
      const double G_k =
          fanno_channel_flow(Pt0, Tt0, X, Pt0 - drive, true, L_k, D, 0.0, "haaland").G;
      EXPECT_LT(std::abs(G_k - prev_G), 0.05 * at_ref.G) << "jump at L=" << L_k;
      prev_G = G_k;
    }
  }

  // 4. Above chi = 0.25: unchoked short exit_total ducts (L = 3.0e-5 m at D = 0.028 m
  //    with chi ~ 0.35; 1 bar L/D = 1e-3 with chi ~ 0.40; 1 bar L/D = 1e-2 with
  //    chi ~ 0.92; 1 kPa L/D = 1 with chi ~ 1.0) must recover exact Hagen-Poiseuille
  //    as drive -> 0 (fails if window is moved higher above chi = 0.35).
  struct ShortCase {
    double P0;
    double L;
    double D_c;
  };
  const std::vector<ShortCase> short_cases = {
      {1.0e5, 3.0e-5, 0.028},
      {1.0e5, 1.0e-5, 0.01},
      {1.0e5, 1.0e-4, 0.01},
      {1.0e3, 1.0e-2, 0.01},
  };
  for (const ShortCase& sc : short_cases) {
    // Use tiny drive well below dp_inert so dynamic/turbulent head is negligible,
    // and evaluate hp_flux at dp_actual = P0 - (P0 - dp) to avoid float64 ULP bias.
    const double dp = 1e-9;
    const double P_target = sc.P0 - dp;
    const double dp_actual = sc.P0 - P_target;
    const double g_hp = hp_flux(sc.P0, Tt0, X, dp_actual, sc.L, sc.D_c, "haaland");
    const FannoChannelFlow r =
        fanno_channel_flow(sc.P0, Tt0, X, P_target, true, sc.L, sc.D_c, 0.0, "haaland");
    EXPECT_NEAR(r.G / g_hp, 1.0, 1e-4)
        << "Pt0=" << sc.P0 << " L=" << sc.L << " D=" << sc.D_c << " t=" << t;
  }
}
