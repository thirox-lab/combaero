// The network friction kernels are laminar below transition: f Re = 64
// (Hagen-Poiseuille) down to essentially zero flow, for every model, and the
// turbulent fits are untouched above it. The laminar branch used to share the
// turbulent fits' Re floor of 64, so f saturated at 1 and a channel's drop
// went as m^2 at low flow; Colebrook had no laminar branch at all.

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "composition.h"
#include "friction.h"
#include "math_constants.h"
#include "solver_interface.h"
#include "thermo.h"
#include "transport.h"

using combaero::solver::friction_and_jacobian;

namespace {

const char* const kModels[] = {"haaland", "serghides", "colebrook", "petukhov"};

double f_of(const std::string& m, double Re, double eD) {
  return std::get<0>(friction_and_jacobian(m, Re, eD).result);
}

double df_of(const std::string& m, double Re, double eD) {
  return std::get<1>(friction_and_jacobian(m, Re, eD).result);
}

}  // namespace

TEST(LaminarFriction, IsHagenPoiseuilleBelowTransition) {
  for (const char* m : kModels) {
    for (double Re : {0.01, 0.1, 1.0, 10.0, 64.0, 200.0, 500.0, 1000.0}) {
      for (double eD : {0.0, 1e-3}) {
        // The blend's turbulent share at Re 1000 is 1e-3: the bound covers it.
        EXPECT_NEAR(f_of(m, Re, eD) * Re, 64.0, 64.0 * 2e-3) << m << " Re " << Re;
      }
    }
  }
}

TEST(LaminarFriction, LeavesTheTurbulentFitsUntouched) {
  for (double Re : {1e4, 1e5, 1e6}) {
    for (double eD : {0.0, 1e-4, 1e-2}) {
      EXPECT_NEAR(f_of("haaland", Re, eD), friction_haaland(Re, eD),
                  1e-12 * friction_haaland(Re, eD));
      EXPECT_NEAR(f_of("serghides", Re, eD), friction_serghides(Re, eD),
                  1e-12 * friction_serghides(Re, eD));
      EXPECT_NEAR(f_of("colebrook", Re, eD), friction_colebrook(Re, eD),
                  1e-9 * friction_colebrook(Re, eD));
      // Petukhov at the turbulent fits' regularized Re, sqrt(Re^2 + 64^2):
      // that shifts f by 5.8e-6 relative at Re 1e4 (#1210).
      const double re_eff = std::hypot(Re, laminar_friction::turbulent_re_floor);
      EXPECT_NEAR(f_of("petukhov", Re, eD), friction_petukhov(re_eff),
                  1e-12 * friction_petukhov(re_eff));
      EXPECT_NEAR(f_of("petukhov", Re, eD), friction_petukhov(Re), 1e-5 * friction_petukhov(Re));
    }
  }
}

// A channel's drop, f Re |m| ~ m: finite and linear as the flow stops.
TEST(LaminarFriction, KeepsTheDropLinearAsTheFlowStops) {
  for (const char* m : kModels) {
    const double ratio = (f_of(m, 1e-2, 0.0) * 1e-2) / (f_of(m, 1.0, 0.0) * 1.0);
    EXPECT_NEAR(ratio, 1.0, 1e-6) << m;
    EXPECT_TRUE(std::isfinite(f_of(m, 0.0, 0.0))) << m;
  }
}

// The drop goes as f Re^2: it must rise with the flow at every Re. The raw
// Petukhov fit has a pole at Re 7.97 and fell 700x across it (#1210).
TEST(LaminarFriction, DropRisesWithFlowAtEveryRe) {
  for (const char* m : kModels) {
    double prev = 0.0;
    for (double Re = 1e-3; Re < 1e6; Re *= 1.02) {
      const double f = f_of(m, Re, 0.0);
      ASSERT_TRUE(std::isfinite(f) && f > 0.0) << m << " Re " << Re;
      const double drop = f * Re * Re;
      ASSERT_GT(drop, prev) << m << " Re " << Re;
      ASSERT_GT(df_of(m, Re, 0.0) * Re * Re + 2.0 * f * Re, 0.0) << m << " Re " << Re;
      prev = drop;
    }
  }
}

TEST(LaminarFriction, DerivativesMatchCentralDifferences) {
  for (const char* m : kModels) {
    for (double Re : {0.002, 1.0, 50.0, 800.0, 1900.0, 2300.0, 2700.0, 3500.0, 2e4}) {
      for (double eD : {0.0, 1e-3}) {
        const double h = 1e-6 * Re;
        const double fd = (f_of(m, Re + h, eD) - f_of(m, Re - h, eD)) / (2.0 * h);
        EXPECT_NEAR(df_of(m, Re, eD), fd, 1e-6 * std::abs(fd) + 1e-12)
            << m << " Re " << Re << " eD " << eD;
      }
    }
  }
}

// The drop's slope at EXACTLY zero flow (#1156). |m| f is 0 * inf there (f ->
// 64 / Re), and the kernel's slope formula read 0 where the drop is laminar,
// linear in m: Newton stalled at an exact zero start. The slope is the limit
// f_re_0 mu L / (2 rho A D^2), Hagen-Poiseuille's 32 mu L / (rho A D^2) times
// the blend's laminar share at Re = 0; the drop itself is unchanged.
namespace {

std::vector<double> air_Y() {
  std::vector<double> Y(combaero::num_species(), 0.0);
  for (std::size_t i = 0; i < Y.size(); ++i) {
    if (combaero::species_name(i) == "N2") {
      Y[i] = 0.767;
    } else if (combaero::species_name(i) == "O2") {
      Y[i] = 0.233;
    }
  }
  return Y;
}

combaero::solver::ChannelResult channel_at(double m, const std::string& model, double f_mult,
                                           bool exit_head) {
  return combaero::solver::channel_residuals_and_jacobian(m, 1e5, 1e5, 300.0, air_Y(), 0.99e5,
                                                          1.0, 0.01, 1e-5, model, f_mult,
                                                          exit_head);
}

}  // namespace

TEST(LaminarFriction, ZeroReLimitOfFRe) {
  using namespace laminar_friction;
  const double w0 = 0.5 * (1.0 + std::tanh(-blend_center / blend_width));
  for (const char* m : kModels) {
    EXPECT_DOUBLE_EQ(combaero::solver::friction_f_re_at_zero(m), f_re * (1.0 - w0)) << m;
    // The kernel's own f Re approaches it once Re is well above the floor:
    // at Re 0.1 the floor's share is (re_floor / Re)^2 / 2 = 5e-11 and the
    // blend weight's own slope moves f Re by ~1e-7.
    EXPECT_NEAR(f_of(m, 0.1, 0.0) * 0.1, combaero::solver::friction_f_re_at_zero(m), 1e-6)
        << m;
  }
  EXPECT_THROW(combaero::solver::friction_f_re_at_zero("nope"), std::invalid_argument);
}

TEST(LaminarFriction, ChannelSlopeAtZeroFlowIsTheLaminarLimit) {
  const std::vector<double> X = combaero::mass_to_mole(air_Y());
  const double rho = combaero::density(300.0, 1e5, X);
  const double mu = combaero::viscosity(300.0, 1e5, X);
  const double D = 0.01;
  const double A = 0.25 * M_PI * D * D;
  for (const char* m : kModels) {
    for (bool exit_head : {false, true}) {
      const auto r0 = channel_at(0.0, m, 1.3, exit_head);
      EXPECT_EQ(r0.dP_calc, 0.0) << m;
      // Built from the blend's constants, not from friction_f_re_at_zero, so
      // the check is independent of the function under test: Hagen-Poiseuille
      // 32 mu L / (rho A D^2) times the laminar share 1 - w(0), w(0) = 1.013e-5.
      const double w0 = 0.5 * (1.0 + std::tanh(-laminar_friction::blend_center /
                                               laminar_friction::blend_width));
      const double hp = 1.3 * 32.0 * mu * 1.0 / (rho * A * D * D) * (1.0 - w0);
      EXPECT_NEAR(r0.d_dP_d_mdot, hp, 1e-12 * hp) << m;
      // The one-sided limits from either side, well above the laminar
      // re_floor (Re 1e-6, |m| ~ 1e-13 kg/s here; Re 0.07 at 1e-8): the
      // slope is continuous through zero. The exit head, m|m| / (2 rho A^2),
      // adds |m| / (rho A^2) to the slope and half that to dP / m.
      for (double mm : {1e-8, -1e-8}) {
        const auto r = channel_at(mm, m, 1.3, exit_head);
        const double q = exit_head ? std::abs(mm) / (rho * A * A) : 0.0;
        EXPECT_NEAR(r.d_dP_d_mdot, r0.d_dP_d_mdot + q, 1e-7 * hp) << m << " m " << mm;
        EXPECT_NEAR(r.dP_calc / mm, r0.d_dP_d_mdot + 0.5 * q, 1e-7 * hp) << m << " m " << mm;
      }
    }
  }
}
