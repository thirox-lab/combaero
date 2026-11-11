// Combustor pressure-loss correlations (#718): every partial against central
// differences of the value, and the head family against its definition
// written out through the density (q = 0.5 rho v^2 / P, v = mdot/(rho A)).

#include "combustor_loss.h"
#include "composition.h"
#include "thermo.h"
#include "thermo_transport_data.h"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cb = combaero;

namespace {

std::vector<double> air_with_fuel() {
  std::vector<double> Y(cb::molar_masses.size(), 0.0);
  Y[cb::species_index_from_name("N2")] = 0.73;
  Y[cb::species_index_from_name("O2")] = 0.20;
  Y[cb::species_index_from_name("CO2")] = 0.03;
  Y[cb::species_index_from_name("H2O")] = 0.02;
  Y[cb::species_index_from_name("CH4")] = 0.02;
  return Y;
}

double head_xi(double th, double m, double T, double P, const std::vector<double> &Y) {
  return cb::combustor_head_loss(th, m, T, P, Y, 1.0, 3.0, 0.05).xi;
}

} // namespace

TEST(CombustorLoss, HeadLossIsTheDynamicHeadThroughTheDensity) {
  const std::vector<double> Y = air_with_fuel();
  const double T = 650.0, P = 1.6e5, m = 0.8, A = 0.05, theta = 1.7;
  const double rho = cb::density(T, P, cb::mass_to_mole(Y));
  const double v = m / (rho * A);
  const double expected = (1.0 * theta + 3.0) * 0.5 * rho * v * v / P;
  EXPECT_NEAR(head_xi(theta, m, T, P, Y) / expected, 1.0, 1e-12);
}

TEST(CombustorLoss, HeadLossPartialsMatchDifferences) {
  const std::vector<double> Y = air_with_fuel();
  const double x0[4] = {1.7, 0.8, 650.0, 1.6e5};
  const auto r = cb::combustor_head_loss(x0[0], x0[1], x0[2], x0[3], Y, 1.0, 3.0, 0.05);
  const double an[4] = {r.dxi_dtheta, r.dxi_dmdot, r.dxi_dTt_in, r.dxi_dPt_in};
  for (int i = 0; i < 4; ++i) {
    double a[4] = {x0[0], x0[1], x0[2], x0[3]};
    double b[4] = {x0[0], x0[1], x0[2], x0[3]};
    const double h = std::abs(x0[i]) * 1e-6;
    a[i] += h;
    b[i] -= h;
    const double fd = (head_xi(a[0], a[1], a[2], a[3], Y) - head_xi(b[0], b[1], b[2], b[3], Y)) / (2.0 * h);
    EXPECT_NEAR(an[i] / fd, 1.0, 1e-7) << "partial " << i;
  }
  ASSERT_EQ(r.dxi_dY.size(), Y.size());
  for (std::size_t k = 0; k < Y.size(); ++k) {
    std::vector<double> a = Y, b = Y;
    a[k] += 1e-6;
    b[k] -= 1e-6;
    const double fd = (head_xi(x0[0], x0[1], x0[2], x0[3], a) - head_xi(x0[0], x0[1], x0[2], x0[3], b)) / 2e-6;
    EXPECT_NEAR(r.dxi_dY[k], fd, 1e-8 * r.xi) << "species " << k;
  }
}

TEST(CombustorLoss, HeadLossVanishesWithoutFlow) {
  const auto r = cb::combustor_head_loss(1.0, 0.0, 600.0, 1e5, air_with_fuel(), 1.0, 3.0, 0.05);
  EXPECT_EQ(r.xi, 0.0);
  EXPECT_EQ(r.dxi_dmdot, 0.0);
}

TEST(CombustorLoss, FractionLossFadesBelowTheDesignFlow) {
  const double m_ref = 0.5;
  // At and above the design flow: exactly k theta + xi0, flat in mdot.
  const auto at = cb::combustor_fraction_loss(1.5, 0.7, 0.01, 0.02, m_ref);
  EXPECT_DOUBLE_EQ(at.xi, 0.035);
  EXPECT_DOUBLE_EQ(at.dxi_dtheta, 0.01);
  EXPECT_EQ(at.dxi_dmdot, 0.0);
  EXPECT_TRUE(at.dxi_dY.empty());
  // Below it: g = 1 - (1 - r^2)^3, with its mdot partial, for either sign.
  for (double m : {0.3, -0.3}) {
    const auto r = cb::combustor_fraction_loss(1.5, m, 0.01, 0.02, m_ref);
    const double rr = std::abs(m) / m_ref;
    EXPECT_NEAR(r.xi, 0.035 * (1.0 - std::pow(1.0 - rr * rr, 3)), 1e-15);
    const double h = 1e-7;
    const double fd = (cb::combustor_fraction_loss(1.5, m + h, 0.01, 0.02, m_ref).xi -
                       cb::combustor_fraction_loss(1.5, m - h, 0.01, 0.02, m_ref).xi) / (2.0 * h);
    EXPECT_NEAR(r.dxi_dmdot / fd, 1.0, 1e-7);
  }
  // No design flow: a constant fraction at every flow.
  EXPECT_DOUBLE_EQ(cb::combustor_fraction_loss(0.0, 1e-6, 0.0, 0.03, 0.0).xi, 0.03);
}

// ---------------------------------------------------------------------------
// pressure_loss_compose (#725): each kind is its kernel, bit for bit, and
// every partial matches central differences of the composed value.

namespace {

cb::PressureLossModel fraction_model(double k, double xi0, double m_ref) {
  cb::PressureLossModel m;
  m.kind = cb::PressureLossKind::Fraction;
  m.k = k;
  m.base = xi0;
  m.m_ref = m_ref;
  return m;
}

cb::PressureLossModel head_model(double k, double zeta0, double area) {
  cb::PressureLossModel m;
  m.kind = cb::PressureLossKind::Head;
  m.k = k;
  m.base = zeta0;
  m.area = area;
  return m;
}

cb::PressureLossInputs inputs_at(double theta, double mdot, double Tt, double Pt,
                                 const std::vector<double> &Y) {
  cb::PressureLossInputs in;
  in.mdot = mdot;
  in.Tt_in = Tt;
  in.Pt_in = Pt;
  in.Y_in = Y;
  in.x = {theta};
  return in;
}

} // namespace

TEST(PressureLossCompose, EachKindIsItsKernel) {
  const std::vector<double> Y = air_with_fuel();
  for (double m_ref : {0.0, 0.5, 2.0}) {
    const auto k = cb::combustor_fraction_loss(1.7, 0.8, 0.01, 0.02, m_ref);
    const auto c = cb::pressure_loss_compose(fraction_model(0.01, 0.02, m_ref),
                                             inputs_at(1.7, 0.8, 650.0, 1.6e5, Y));
    EXPECT_EQ(c.xi, k.xi);
    EXPECT_EQ(c.dxi_dmdot, k.dxi_dmdot);
    ASSERT_EQ(c.dxi_dx.size(), 1u);
    EXPECT_EQ(c.dxi_dx[cb::PRESSURE_LOSS_INPUT_THETA], k.dxi_dtheta);
    EXPECT_TRUE(c.dxi_dY.empty());
    EXPECT_EQ(c.dxi_dTt_in, 0.0);
    EXPECT_EQ(c.dxi_dPt_in, 0.0);
  }
  const auto k = cb::combustor_head_loss(1.7, 0.8, 650.0, 1.6e5, Y, 1.0, 3.0, 0.05);
  const auto c = cb::pressure_loss_compose(head_model(1.0, 3.0, 0.05),
                                           inputs_at(1.7, 0.8, 650.0, 1.6e5, Y));
  EXPECT_EQ(c.xi, k.xi);
  EXPECT_EQ(c.dxi_dmdot, k.dxi_dmdot);
  EXPECT_EQ(c.dxi_dTt_in, k.dxi_dTt_in);
  EXPECT_EQ(c.dxi_dPt_in, k.dxi_dPt_in);
  EXPECT_EQ(c.dxi_dx[cb::PRESSURE_LOSS_INPUT_THETA], k.dxi_dtheta);
  ASSERT_EQ(c.dxi_dY.size(), Y.size());
  for (std::size_t i = 0; i < Y.size(); ++i) {
    EXPECT_EQ(c.dxi_dY[i], k.dxi_dY[i]);
  }
}

TEST(PressureLossCompose, PartialsMatchDifferences) {
  const std::vector<double> Y = air_with_fuel();
  const cb::PressureLossModel models[] = {
      fraction_model(0.01, 0.02, 2.0), // below the design flow: live mdot partial
      fraction_model(0.01, 0.02, 0.0), head_model(1.0, 3.0, 0.05), head_model(0.0, 4.0, 0.05)};
  const double x0[4] = {1.7, 0.8, 650.0, 1.6e5}; // theta, mdot, Tt, Pt
  for (const auto &model : models) {
    auto xi_at = [&](const double *v, const std::vector<double> &Yv) {
      return cb::pressure_loss_compose(model, inputs_at(v[0], v[1], v[2], v[3], Yv)).xi;
    };
    const auto r = cb::pressure_loss_compose(model, inputs_at(x0[0], x0[1], x0[2], x0[3], Y));
    const double an[4] = {r.dxi_dx[cb::PRESSURE_LOSS_INPUT_THETA], r.dxi_dmdot, r.dxi_dTt_in,
                          r.dxi_dPt_in};
    for (int i = 0; i < 4; ++i) {
      double a[4] = {x0[0], x0[1], x0[2], x0[3]};
      double b[4] = {x0[0], x0[1], x0[2], x0[3]};
      const double h = std::abs(x0[i]) * 1e-6;
      a[i] += h;
      b[i] -= h;
      const double fd = (xi_at(a, Y) - xi_at(b, Y)) / (2.0 * h);
      EXPECT_NEAR(an[i], fd, 1e-7 * std::abs(fd) + 1e-14) << "partial " << i;
    }
    for (std::size_t k = 0; k < r.dxi_dY.size(); ++k) {
      std::vector<double> a = Y, b = Y;
      a[k] += 1e-6;
      b[k] -= 1e-6;
      const double fd = (xi_at(x0, a) - xi_at(x0, b)) / 2e-6;
      EXPECT_NEAR(r.dxi_dY[k], fd, 1e-8 * r.xi) << "species " << k;
    }
  }
}

TEST(PressureLossCompose, TheValueReadsXThePartialsY) {
  // An unnormalised Y (a perturbed Newton iterate): the value and its theta
  // partial are the kernel's at mole_to_mass(X), the state partials the
  // kernel's at the raw Y.
  std::vector<double> Y = air_with_fuel();
  for (double &y : Y) {
    y *= 1.03;
  }
  std::vector<double> Yn = air_with_fuel();
  auto in = inputs_at(1.7, 0.8, 650.0, 1.6e5, Y);
  in.X_in = cb::mass_to_mole(Yn);
  const auto c = cb::pressure_loss_compose(head_model(1.0, 3.0, 0.05), in);
  const auto at_x = cb::combustor_head_loss(1.7, 0.8, 650.0, 1.6e5, cb::mole_to_mass(in.X_in),
                                            1.0, 3.0, 0.05);
  const auto at_y = cb::combustor_head_loss(1.7, 0.8, 650.0, 1.6e5, Y, 1.0, 3.0, 0.05);
  EXPECT_EQ(c.xi, at_x.xi);
  EXPECT_EQ(c.dxi_dx[cb::PRESSURE_LOSS_INPUT_THETA], at_x.dxi_dtheta);
  EXPECT_EQ(c.dxi_dmdot, at_y.dxi_dmdot);
  EXPECT_EQ(c.dxi_dTt_in, at_y.dxi_dTt_in);
  EXPECT_EQ(c.dxi_dPt_in, at_y.dxi_dPt_in);
  for (std::size_t i = 0; i < Y.size(); ++i) {
    EXPECT_EQ(c.dxi_dY[i], at_y.dxi_dY[i]);
  }
}

TEST(PressureLossCompose, InputsAModelDoesNotReadAreIgnored) {
  const std::vector<double> Y = air_with_fuel();
  // k = 0: theta is not read, x may be empty or carry anything.
  cb::PressureLossInputs in = inputs_at(0.0, 0.8, 650.0, 1.6e5, Y);
  in.x.clear();
  const auto c0 = cb::pressure_loss_compose(head_model(0.0, 4.0, 0.05), in);
  in.x = {std::nan(""), 3.0, 4.0};
  const auto c1 = cb::pressure_loss_compose(head_model(0.0, 4.0, 0.05), in);
  EXPECT_EQ(c0.xi, c1.xi);
  EXPECT_TRUE(c0.dxi_dx.empty());
  ASSERT_EQ(c1.dxi_dx.size(), 3u);
  for (double d : c1.dxi_dx) {
    EXPECT_EQ(d, 0.0);
  }
  // A fraction without a design flow does not read the flow.
  in.mdot = std::nan("");
  EXPECT_EQ(cb::pressure_loss_compose(fraction_model(0.0, 0.03, 0.0), in).xi, 0.03);
  EXPECT_THROW(cb::pressure_loss_compose(fraction_model(0.0, 0.03, 0.5), in),
               std::invalid_argument);
  // k != 0 needs theta.
  in.mdot = 0.8;
  in.x.clear();
  EXPECT_THROW(cb::pressure_loss_compose(fraction_model(0.01, 0.02, 0.0), in),
               std::invalid_argument);
}
