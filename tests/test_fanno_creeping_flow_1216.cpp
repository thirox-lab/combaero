#include <gtest/gtest.h>

#include "../include/composition.h"
#include "../include/fanno_mach.h"
#include "../include/humidair.h"
#include "../include/solver_interface.h"
#include "../include/thermo.h"
#include "../include/transport.h"

#include <cmath>
#include <string>
#include <tuple>
#include <vector>

using namespace combaero;

// ---------------------------------------------------------------------------
// #1216: Compressible ducts at creeping flow must recover Hagen-Poiseuille
// (to the laminar blend's 1 - w(0)) across all L/D, pressures and exit
// coupling modes, with C1 continuity of G and all Jacobian partials across
// kFannoFlowDriveFloor.
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

  for (const std::string& model : {"haaland", "serghides", "colebrook", "petukhov"}) {
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
            const double hy = 1e-4;
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
