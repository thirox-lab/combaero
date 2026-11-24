#include <gtest/gtest.h>

#include "../include/fanno_mach.h"
#include "../include/humidair.h"
#include "../include/thermo.h"

#include <cmath>
#include <string>
#include <vector>

using namespace combaero;

// ---------------------------------------------------------------------------
// Below kFannoFlowDriveFloor the flux is blended to h(t) G0, but the end
// states were returned at the floor: M_in at a 0.035 Pa drive on 20.8 bar was
// 23x the flow's, which made the second-law heat band refuse a converged root
// (stress random 99, #1017). The Mach numbers now follow G; G and its
// partials are unchanged.
// ---------------------------------------------------------------------------

namespace {

constexpr double kPt0 = 20.8e5;
constexpr double kTt0 = 300.0;
constexpr double kTtExit = 400.0;
constexpr double kL = 0.3;
constexpr double kD = 0.1;
const std::string kModel = "haaland";

FannoChannelFlow flow(double drive, bool heated, double Pt0 = kPt0) {
  const std::vector<double> X = dry_air();
  if (heated) {
    return fanno_rayleigh_channel_flow(Pt0, kTt0, kTtExit, X, Pt0 - drive, false, kL, kD, 0.0,
                                       kModel);
  }
  return fanno_channel_flow(Pt0, kTt0, X, Pt0 - drive, false, kL, kD, 0.0, kModel);
}

// Mach of flux G at the inlet, G sqrt(R T) / (P sqrt(gamma)) at the stagnation
// state: exact to O(M^2), ~1e-9 here.
double mach_from_flux(double G, double Pt0) {
  const std::vector<double> X = dry_air();
  const double R = specific_gas_constant(X);
  const double g = isentropic_expansion_coefficient(kTt0, X);
  return G * std::sqrt(R * kTt0) / (Pt0 * std::sqrt(g));
}

}  // namespace

TEST(FannoDriveFloorMach1017, InletMachFollowsTheBlendedFlux) {
  for (const bool heated : {false, true}) {
    for (const double drive : {0.035, 0.35, 0.9}) {
      const FannoChannelFlow f = flow(drive, heated);
      ASSERT_GT(f.G, 0.0);
      const double M_ref = mach_from_flux(f.G, kPt0);
      EXPECT_NEAR(f.M_in / M_ref, 1.0, 0.01) << "heated " << heated << " drive " << drive;
      // The exit Mach keeps the floor's exit-to-inlet ratio: heating raises it,
      // the adiabatic duct barely moves it.
      const FannoChannelFlow at = flow(kFannoFlowDriveFloor, heated);
      EXPECT_NEAR(f.M_exit / f.M_in, at.M_exit / at.M_in, 1e-9);
      // The warm start stays the root actually solved, the floor's.
      EXPECT_DOUBLE_EQ(f.M_exit_root, at.M_exit);
    }
  }
  // 1 bar, 0.01 Pa: the floor's Mach was 80x the flow's.
  const FannoChannelFlow f = flow(0.01, false, 1e5);
  EXPECT_NEAR(f.M_in / mach_from_flux(f.G, 1e5), 1.0, 0.01);
}

TEST(FannoDriveFloorMach1017, MachIsContinuousAcrossTheFloor) {
  for (const bool heated : {false, true}) {
    const FannoChannelFlow below = flow(kFannoFlowDriveFloor * (1.0 - 1e-9), heated);
    const FannoChannelFlow above = flow(kFannoFlowDriveFloor * (1.0 + 1e-9), heated);
    EXPECT_NEAR(below.M_in / above.M_in, 1.0, 1e-7) << "heated " << heated;
    EXPECT_NEAR(below.M_exit / above.M_exit, 1.0, 1e-7) << "heated " << heated;
    EXPECT_NEAR(below.G / above.G, 1.0, 1e-7) << "heated " << heated;
    EXPECT_DOUBLE_EQ(above.M_exit_root, above.M_exit);
  }
}

TEST(FannoDriveFloorMach1017, MachGoesToZeroWithTheDrive) {
  // With #1216, below the floor the duct follows Hagen-Poiseuille as drive -> 0
  // rather than extrapolating the inertia-dominated 1 Pa flux linearly: on this
  // wide 20.8 bar duct (L/D = 3) G(1e-6 Pa) / G(1 Pa) ~ 2e-4 (vs 1.25e-6 under
  // the old cubic), and scales linearly to zero as drive falls further.
  for (const bool heated : {false, true}) {
    const FannoChannelFlow f6 = flow(1e-6, heated);
    const FannoChannelFlow f9 = flow(1e-9, heated);
    const FannoChannelFlow at = flow(kFannoFlowDriveFloor, heated);
    EXPECT_GT(f6.M_in, 0.0);
    EXPECT_LT(f6.M_in, 1e-3 * at.M_in) << "heated " << heated;
    EXPECT_GT(f9.M_in, 0.0);
    EXPECT_LT(f9.M_in, 1e-6 * at.M_in) << "heated " << heated;
    EXPECT_NEAR(f9.M_in / f6.M_in, 1e-3, 1e-4) << "heated " << heated;
  }
}
