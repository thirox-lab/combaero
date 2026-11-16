#include <gtest/gtest.h>

#include "../include/composition.h"
#include "../include/fanno_mach.h"
#include "../include/humidair.h"
#include "../include/thermo.h"

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

using namespace combaero;

// ---------------------------------------------------------------------------
// dG_ddir differences the duct along a mass-fraction direction v. Where the
// central step would take a fraction below zero it is one-sided. A component
// of -2.2e-20 on an absent species (round-off in the solver's orthonormal
// basis) switched the stencil, and a first-order one at that: J moved by
// 2.3e-5 relative (stress random 30, #1199). A round-off component now no
// longer decides the stencil, and the one-sided stencil is second order.
// ---------------------------------------------------------------------------

namespace {

constexpr double kPt0 = 2.0e5;
constexpr double kTt0 = 600.0;
constexpr double kTtExit = 800.0;
constexpr double kL = 0.5;
constexpr double kD = 0.02;
const std::string kModel = "haaland";

// Dry air with its CO2 moved into N2: CO2 absent (Y = 0).
std::vector<double> air_without_co2() {
  std::vector<double> Y = mole_to_mass(dry_air());
  const std::size_t i_co2 = species_index_from_name("CO2");
  Y[species_index_from_name("N2")] += Y[i_co2];
  Y[i_co2] = 0.0;
  return Y;
}

// A unit direction trading species a for species b.
std::vector<double> trade(const std::string& a, const std::string& b) {
  std::vector<double> v(num_species(), 0.0);
  v[species_index_from_name(a)] = 1.0 / std::sqrt(2.0);
  v[species_index_from_name(b)] = -1.0 / std::sqrt(2.0);
  return v;
}

double ddir(const std::vector<double>& Y, const std::vector<double>& v, double P_target,
            bool heated) {
  const std::vector<double> X = mass_to_mole(Y);
  const std::vector<std::vector<double>> dirs{v};
  const FannoChannelFlow f =
      heated ? fanno_rayleigh_channel_flow(kPt0, kTt0, kTtExit, X, P_target, false, kL, kD, 0.0,
                                           kModel, 1.0, -1.0, dirs)
             : fanno_channel_flow(kPt0, kTt0, X, P_target, false, kL, kD, 0.0, kModel, 1.0, -1.0,
                                  dirs);
  EXPECT_EQ(f.dG_ddir.size(), 1u);
  return f.dG_ddir.empty() ? 0.0 : f.dG_ddir[0];
}

}  // namespace

TEST(FannoDdirRoundoff1199, RoundoffComponentOnAnAbsentSpeciesDoesNotMoveTheResult) {
  const std::vector<double> Y = air_without_co2();
  const std::size_t i_co2 = species_index_from_name("CO2");
  for (const bool heated : {false, true}) {
    for (const double P : {1.5e5, 0.5e5}) {
      std::vector<double> v = trade("O2", "N2");
      const double exact = ddir(Y, v, P, heated);
      ASSERT_TRUE(std::isfinite(exact));
      ASSERT_NE(exact, 0.0);
      for (const double c : {-2.2e-20, 2.2e-20, -1e-14, 1e-300}) {
        v[i_co2] = c;
        // Bit for bit: the component is dropped from the step.
        EXPECT_EQ(ddir(Y, v, P, heated), exact) << "heated " << heated << " P " << P << " c " << c;
      }
    }
  }
}

TEST(FannoDdirRoundoff1199, OneSidedStencilIsSecondOrder) {
  // A genuine component on the absent species: only s > 0 is feasible. The
  // reference extrapolates central derivatives at s = 2h, 4h, 6h (feasible,
  // O(h^2)) back to s = 0, quadratically: weights 3, -3, 1. Its own scatter is
  // ~5e-8; the first-order stencil this replaced was off by 3.3e-6.
  const std::vector<double> Y = air_without_co2();
  const std::vector<double> v = trade("CO2", "N2");
  const double h = kFannoDirStep * std::sqrt(2.0);
  auto at = [&](double s, double P, bool heated) {
    std::vector<double> Ys = Y;
    for (std::size_t k = 0; k < Ys.size(); ++k) {
      Ys[k] += s * v[k];
    }
    return ddir(Ys, v, P, heated);
  };
  for (const bool heated : {false, true}) {
    for (const double P : {1.5e5, 1.0e5}) {
      const double d0 = at(0.0, P, heated);
      const double ref = 3.0 * at(2.0 * h, P, heated) - 3.0 * at(4.0 * h, P, heated) +
                         at(6.0 * h, P, heated);
      ASSERT_TRUE(std::isfinite(ref));
      EXPECT_NEAR(d0 / ref, 1.0, 3e-7) << "heated " << heated << " P " << P;
    }
  }
}

TEST(FannoDdirRoundoff1199, BackwardOneSidedIsTheForwardNegated) {
  // The SVD sign of a solver direction is arbitrary, so -v (one-sided towards
  // s = -h) is as likely as v. dG/d(-v) = -dG/dv exactly: the backward stencil
  // runs the forward one's arithmetic with s negated.
  const std::vector<double> Y = air_without_co2();
  const std::vector<double> v = trade("CO2", "N2");
  std::vector<double> minus_v = v;
  for (double& c : minus_v) {
    c = -c;
  }
  for (const bool heated : {false, true}) {
    const double fwd = ddir(Y, v, 1.5e5, heated);
    ASSERT_NE(fwd, 0.0);
    EXPECT_EQ(ddir(Y, minus_v, 1.5e5, heated), -fwd) << "heated " << heated;
  }
}

TEST(FannoDdirRoundoff1199, HalvedStepWhereTheDoubledStepIsInfeasible) {
  // CO2 absent blocks s < 0; H2O at 1.5e-5 with h|v_H2O| = 1e-5 lets s = h
  // through but blocks 2h. The stencil then runs at half the step, still
  // second order. The first-order fallback it replaced was off by 1.2e-5 here,
  // a jump in dG/dv as Y_H2O crossed 2e-5. Reference: the full-step values at
  // Y_H2O = 3, 4, 5e-5, extrapolated quadratically to 1.5e-5 (weights 4.375,
  // -5.25, 1.875; their scatter amplified, ~2e-7).
  const std::vector<double> v = trade("CO2", "H2O");
  auto with_h2o = [](double y) {
    std::vector<double> Y = air_without_co2();
    Y[species_index_from_name("H2O")] = y;
    Y[species_index_from_name("N2")] -= y;
    return Y;
  };
  std::vector<double> minus_v = v;
  for (double& c : minus_v) {
    c = -c;
  }
  for (const bool heated : {false, true}) {
    for (const double P : {1.5e5, 1.0e5}) {
      const double d = ddir(with_h2o(1.5e-5), v, P, heated);
      const double ref = 4.375 * ddir(with_h2o(3e-5), v, P, heated) -
                         5.25 * ddir(with_h2o(4e-5), v, P, heated) +
                         1.875 * ddir(with_h2o(5e-5), v, P, heated);
      ASSERT_GT(std::abs(ref), 1.0);
      EXPECT_NEAR(d / ref, 1.0, 1e-6) << "heated " << heated << " P " << P;
      EXPECT_EQ(ddir(with_h2o(1.5e-5), minus_v, P, heated), -d);
    }
  }
}

TEST(FannoDdirRoundoff1199, SmallComponentOnAPresentSpeciesIsKept) {
  // Only a component that would take its fraction below zero is dropped. On
  // CO2 present (dry air) a component at the threshold stays in the step, so
  // the result moves (here by FD rounding, ~1e-8; physically by ~4e-12).
  const std::vector<double> Y = mole_to_mass(dry_air());
  std::vector<double> v = trade("O2", "N2");
  const double base = ddir(Y, v, 1.5e5, false);
  v[species_index_from_name("CO2")] = kFannoDirRoundoff / std::sqrt(2.0);
  EXPECT_NE(ddir(Y, v, 1.5e5, false), base);
}
