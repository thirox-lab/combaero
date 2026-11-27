// Regression tests for #1242: adiabatic fanno_channel_flow Mach bracket below
// the former hardcoded 1e-6 floor (raw_flow in src/fanno_mach.cpp).
//
// On main, halving from Me = 0.5 stopped at lo = 2^-20 = 9.536743e-7 <= 1e-6,
// throwing "fanno (Mach) brent: no sign change" whenever the unchoked exit
// Mach at the 1 Pa drive floor lay below 1e-6 (both high-Tt near-zero-start
// solver probes from #1242's net_deadend and normal-Tt narrow laminar ducts).

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

#include "composition.h"
#include "fanno_mach.h"
#include "friction.h"
#include "thermo.h"
#include "transport.h"

using namespace combaero;

namespace {

std::vector<double> air_X() {
  std::vector<double> X(num_species(), 0.0);
  X[species_index_from_name("N2")] = 0.78084;
  X[species_index_from_name("O2")] = 0.20946;
  X[species_index_from_name("AR")] = 0.00934;
  X[species_index_from_name("CO2")] = 0.00036;
  return X;
}

}  // namespace

TEST(FannoMachBracket1242, UltraLowMachStatesSolveCleanlyColdAndWarm) {
  const std::vector<double> X = air_X();
  struct Case {
    const char* name;
    double Pt0, Tt0, drive;
    bool exit_total;
    double L, D, rough;
    std::string model;
    double fmult;
  };
  const Case cases[] = {
      // Exact failing intermediate probes captured from #1242 net_deadend (sc=1e-11).
      {"deadend_probe_5.95e5K", 1.110425e5, 5.953859e5, 1.0, true, 0.5, 0.02, 1e-5, "haaland", 1.0},
      {"deadend_probe_3.08e6K", 1.157294e5, 3.076902e6, 1.0, true, 0.5, 0.02, 1e-5, "haaland", 1.0},
      {"deadend_probe_8.25e6K", 1.157294e5, 8.25e6, 1.0, true, 0.5, 0.02, 1e-5, "haaland", 1.0},
      // Normal-temperature narrow laminar duct: Hagen-Poiseuille Me ~ 1.17e-7 < 1e-6.
      {"narrow_laminar_1000K", 1.0e5, 1000.0, 1.0, true, 0.1, 1.0e-4, 0.0, "haaland", 1.0},
      // Static-pressure exit twin of the narrow laminar duct.
      {"narrow_laminar_static", 1.0e5, 1000.0, 1.0, false, 0.1, 1.0e-4, 0.0, "haaland", 1.0},
      // Sub-floor drive (calls raw_flow at kFannoFlowDriveFloor = 1 Pa).
      {"narrow_laminar_subfloor", 1.0e5, 1000.0, 0.25, true, 0.1, 1.0e-4, 0.0, "haaland", 1.0},
      // Ultra-narrow and capillary ducts pinning the machine-range guard below
      // 1e-11 (measured Me = 1.53e-9, 4.84e-11, 4.84e-12), plus a high-fL/D
      // fixed-f duct (measured Me = 1.2e-8).
      {"ultra_narrow_D1e-5", 1.0e5, 1000.0, 1.0, true, 1.0, 1.0e-5, 0.0, "haaland", 1.0},
      {"capillary_D1e-6_L100", 1.0e5, 1000.0, 1.0, true, 100.0, 1.0e-6, 0.0, "haaland", 1.0},
      {"capillary_D1e-7_L1000", 1.0e5, 1000.0, 1.0, true, 1000.0, 1.0e-7, 0.0, "haaland", 1.0},
      {"fixed_huge_fLD", 1.0e5, 300.0, 1.0, true, 1.0e6, 1.0e-3, 0.0, "fixed", 100.0},
  };

  for (const auto& c : cases) {
    const double P_target = c.Pt0 - c.drive;
    const FannoChannelFlow cold = fanno_channel_flow(
        c.Pt0, c.Tt0, X, P_target, c.exit_total, c.L, c.D, c.rough, c.model, c.fmult, -1.0, {},
        true);
    EXPECT_FALSE(cold.choked) << c.name;
    EXPECT_GT(cold.G, 0.0) << c.name;
    EXPECT_LT(cold.M_exit_root, 1.0e-6) << c.name;
    EXPECT_TRUE(std::isfinite(cold.dG_dP_target) && cold.dG_dP_target < 0.0) << c.name;
    EXPECT_TRUE(std::isfinite(cold.dG_dPt0) && cold.dG_dPt0 > 0.0) << c.name;
    EXPECT_TRUE(std::isfinite(cold.dG_dTt0)) << c.name;

    // Warm starts below, near, and far above the ultra-low-Mach root must agree
    // with the cold root to the inlet-solve rounding floor eps / (drive / Pt0)
    // (~2e-11 at 1 Pa on 1 bar; on main, Me_guess < 1e-6 clamped lo up to 1e-6).
    for (const double s : {0.5, 0.99, 1.01, 20.0}) {
      const FannoChannelFlow warm = fanno_channel_flow(
          c.Pt0, c.Tt0, X, P_target, c.exit_total, c.L, c.D, c.rough, c.model, c.fmult,
          s * cold.M_exit_root, {}, false);
      EXPECT_FALSE(warm.choked) << c.name << " scale=" << s;
      EXPECT_NEAR(warm.G / cold.G, 1.0, 1e-9) << c.name << " scale=" << s;
    }

    // Above the 1 Pa floor, verify the solved (G, M_in, M_exit) satisfies
    // fanno_length_between(G, Tt0, M_in, M_exit) == L.
    if (c.drive >= kFannoFlowDriveFloor) {
      const double L_check = fanno_length_between(
          cold.G, c.Tt0, cold.M_in, cold.M_exit, X, c.D, c.rough, c.model, c.fmult);
      EXPECT_NEAR(L_check, c.L, 1e-8 * c.L) << c.name;
    }
  }
}

TEST(FannoMachBracket1242, NarrowLaminarMatchesHagenPoiseuilleAndDerivativesMatchFD) {
  const std::vector<double> X = air_X();
  const double Pt0 = 1.0e5;
  const double Tt0 = 1000.0;
  const double drive = 1.0;
  const double P_target = Pt0 - drive;
  const double L = 0.1;
  const double D = 1.0e-4;

  const FannoChannelFlow r =
      fanno_channel_flow(Pt0, Tt0, X, P_target, true, L, D, 0.0, "haaland", 1.0, -1.0, {}, true);
  // Hagen-Poiseuille check: G = rho * u = (Pt0 / (R Tt0)) * (2 drive D^2 / (64 mu L)).
  const double R = specific_gas_constant(X);
  const double mu = viscosity(Tt0, Pt0, X);
  const double G_hp = (Pt0 / (R * Tt0)) * (2.0 * drive * D * D / (laminar_friction::f_re * mu * L));
  EXPECT_NEAR(r.G / G_hp, 1.0, 1e-3);

  // Analytic partials vs central differences at drive = 2.0 Pa (clear of the 1 Pa floor).
  const double d2 = 2.0;
  const FannoChannelFlow r2 =
      fanno_channel_flow(Pt0, Tt0, X, Pt0 - d2, true, L, D, 0.0, "haaland", 1.0, -1.0, {}, true);
  const double hP = 1e-3;
  const double fd_P =
      (fanno_channel_flow(Pt0, Tt0, X, Pt0 - d2 + hP, true, L, D, 0.0, "haaland", 1.0, -1.0, {},
                          false)
           .G -
       fanno_channel_flow(Pt0, Tt0, X, Pt0 - d2 - hP, true, L, D, 0.0, "haaland", 1.0, -1.0, {},
                          false)
           .G) /
      (2.0 * hP);
  const double fd_Pt =
      (fanno_channel_flow(Pt0 + hP, Tt0, X, Pt0 - d2, true, L, D, 0.0, "haaland", 1.0, -1.0, {},
                          false)
           .G -
       fanno_channel_flow(Pt0 - hP, Tt0, X, Pt0 - d2, true, L, D, 0.0, "haaland", 1.0, -1.0, {},
                          false)
           .G) /
      (2.0 * hP);
  const double hT = 1e-2;
  const double fd_T =
      (fanno_channel_flow(Pt0, Tt0 + hT, X, Pt0 - d2, true, L, D, 0.0, "haaland", 1.0, -1.0, {},
                          false)
           .G -
       fanno_channel_flow(Pt0, Tt0 - hT, X, Pt0 - d2, true, L, D, 0.0, "haaland", 1.0, -1.0, {},
                          false)
           .G) /
      (2.0 * hT);

  EXPECT_NEAR(r2.dG_dP_target / fd_P, 1.0, 1e-5);
  EXPECT_NEAR(r2.dG_dPt0 / fd_Pt, 1.0, 1e-5);
  EXPECT_NEAR(r2.dG_dTt0 / fd_T, 1.0, 1e-5);
}
