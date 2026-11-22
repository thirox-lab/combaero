// The wall's two-stream exchange (#627): epsilon-NTU as the textbook writes
// it (Incropera, Table 11.3) for both arrangements, symmetric in the sides,
// finite at a reservoir, at C_a = C_b and as a stream stops, the linear
// coupling at NTU -> 0, and analytic derivatives against central differences.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "heat_transfer.h"

using combaero::WallCouplingResult;
using combaero::wall_coupling_and_jacobian;
using combaero::wall_exchange_and_jacobian;

namespace {

const double kInf = std::numeric_limits<double>::infinity();
const std::vector<double> kLayers = {0.0005 / 200.0};
constexpr double kHa = 500.0;
constexpr double kHb = 300.0;
constexpr double kTa = 1500.0;
constexpr double kTb = 400.0;
constexpr double kA = 0.01;

WallCouplingResult exchange(double C_a, double C_b, bool counterflow, double h_a = kHa,
                            double h_b = kHb, double T_a = kTa, double T_b = kTb) {
  return wall_exchange_and_jacobian(h_a, T_a, C_a, h_b, T_b, C_b, kLayers, kA, 0.0,
                                    counterflow);
}

double eps_textbook(double ntu, double cr, bool counterflow) {
  if (!counterflow) {
    return (1.0 - std::exp(-ntu * (1.0 + cr))) / (1.0 + cr);
  }
  if (std::abs(cr - 1.0) < 1e-12) {
    return ntu / (1.0 + ntu);
  }
  const double e = std::exp(-ntu * (1.0 - cr));
  return (1.0 - e) / (1.0 - cr * e);
}

struct Case {
  double C_a;
  double C_b;
};

const Case kCases[] = {{1.0, 2.0}, {2.0, 1.0}, {1.5, 1.5},  {0.3, 30.0}, {1.0, kInf},
                       {kInf, 0.7}, {1e-3, 5.0}, {50.0, 80.0}, {0.05, 0.05}};

}  // namespace

TEST(WallExchange, MatchesTheTextbookEffectivenessForBothArrangements) {
  for (bool cf : {false, true}) {
    for (const Case& c : kCases) {
      const WallCouplingResult r = exchange(c.C_a, c.C_b, cf);
      const double C_min = std::min(c.C_a, c.C_b);
      const double C_max = std::max(c.C_a, c.C_b);
      const double cr = std::isinf(C_max) ? 0.0 : C_min / C_max;
      const double ref = eps_textbook(r.UA / C_min, cr, cf) * C_min * (kTa - kTb);
      EXPECT_NEAR(r.Q, ref, 1e-12 * std::abs(ref)) << cf << " " << c.C_a << " " << c.C_b;
      EXPECT_NEAR(r.effectiveness, r.Q / (C_min * (kTa - kTb)), 1e-14) << cf;
    }
  }
}

TEST(WallExchange, IsSymmetricInItsSides) {
  for (bool cf : {false, true}) {
    for (const Case& c : kCases) {
      const WallCouplingResult ab = exchange(c.C_a, c.C_b, cf, kHa, kHb, kTa, kTb);
      const WallCouplingResult ba = exchange(c.C_b, c.C_a, cf, kHb, kHa, kTb, kTa);
      EXPECT_NEAR(ab.Q, -ba.Q, 1e-12 * std::abs(ab.Q)) << cf << " " << c.C_a;
      EXPECT_NEAR(ab.dQ_dC_a, -ba.dQ_dC_b, 1e-12 * (std::abs(ab.dQ_dC_a) + 1e-300)) << cf;
    }
  }
}

// Q <= C_min dT always, Q -> C dT as a stream stops (it reaches the other's
// inlet temperature), and the two arrangements meet there: a reversing flow
// passes C = 0, so the arrangement flips with no jump in Q.
TEST(WallExchange, IsBoundedAndTheArrangementsMeetAsAStreamStops) {
  const double dT = kTa - kTb;
  for (bool cf : {false, true}) {
    for (const Case& c : kCases) {
      const WallCouplingResult r = exchange(c.C_a, c.C_b, cf);
      EXPECT_LE(r.Q, std::min(c.C_a, c.C_b) * dT * (1.0 + 1e-14)) << cf << " " << c.C_a;
      EXPECT_LE(r.Q, r.UA * dT * (1.0 + 1e-14)) << cf << " " << c.C_a;
    }
    EXPECT_EQ(exchange(0.0, 2.0, cf).Q, 0.0);
    EXPECT_DOUBLE_EQ(exchange(0.0, 2.0, cf).dQ_dC_a, dT);
    EXPECT_NEAR(exchange(1e-9, 2.0, cf).Q / 1e-9, dT, 1e-6 * dT);
  }
  for (double C : {1e-6, 1e-3, 1e-2}) {
    const double par = exchange(C, 2.0, false).Q;
    const double cnt = exchange(C, 2.0, true).Q;
    EXPECT_NEAR(par, cnt, 2.0 * C / 2.0 * C * dT + 1e-12) << C;
  }
}

// #1245: a stopped stream's partial is the forward slope K -> C as C -> 0+
// beside a flowing one, and 0 when both are stopped (K(C, 0) = 0 for any C).
TEST(WallExchange, StoppedStreamsTakeTheirOneSidedPartials) {
  for (bool cf : {false, true}) {
    const double h = 1e-9;
    const WallCouplingResult one = exchange(0.0, 2.0, cf);
    EXPECT_NEAR(one.dQ_dC_a, exchange(h, 2.0, cf).Q / h, 1e-6 * std::abs(one.dQ_dC_a)) << cf;
    EXPECT_EQ(one.dQ_dC_b, 0.0) << cf;
    const WallCouplingResult both = exchange(0.0, 0.0, cf);
    EXPECT_EQ(both.Q, 0.0) << cf;
    EXPECT_EQ(both.dQ_dC_a, exchange(h, 0.0, cf).Q / h) << cf;
    EXPECT_EQ(both.dQ_dC_b, exchange(0.0, h, cf).Q / h) << cf;
    EXPECT_EQ(both.dQ_dC_a, 0.0) << cf;
    EXPECT_EQ(both.dQ_dC_b, 0.0) << cf;
  }
}

// At NTU -> 0 the exchange is the linear wall: same Q and h, T_aw slopes.
TEST(WallExchange, IsTheLinearCouplingAtSmallNtu) {
  const WallCouplingResult lin = wall_coupling_and_jacobian(kHa, kTa, kHb, kTb, kLayers, kA, 0.0);
  for (bool cf : {false, true}) {
    const WallCouplingResult r = exchange(1e9, 2e9, cf);
    EXPECT_NEAR(r.Q, lin.Q, 1e-8 * std::abs(lin.Q));
    EXPECT_NEAR(r.dQ_dh_a, lin.dQ_dh_a, 1e-8 * std::abs(lin.dQ_dh_a));
    EXPECT_NEAR(r.dQ_dT_aw_a, lin.dQ_dT_aw_a, 1e-8 * lin.dQ_dT_aw_a);
    const WallCouplingResult res = exchange(kInf, kInf, cf);
    EXPECT_DOUBLE_EQ(res.Q, lin.Q);
  }
}

TEST(WallExchange, DerivativesMatchCentralDifferences) {
  for (bool cf : {false, true}) {
    for (const Case& c : kCases) {
      const WallCouplingResult r = exchange(c.C_a, c.C_b, cf);
      auto check = [&](double analytic, double value, auto&& at) {
        const double h = 1e-6 * std::abs(value);
        const double fd = (at(value + h).Q - at(value - h).Q) / (2.0 * h);
        EXPECT_NEAR(analytic, fd, 1e-6 * std::abs(fd) + 1e-9 * std::abs(r.Q) / std::abs(value))
            << cf << " " << c.C_a << " " << c.C_b;
      };
      if (std::isfinite(c.C_a)) {
        check(r.dQ_dC_a, c.C_a, [&](double v) { return exchange(v, c.C_b, cf); });
      }
      if (std::isfinite(c.C_b)) {
        check(r.dQ_dC_b, c.C_b, [&](double v) { return exchange(c.C_a, v, cf); });
      }
      check(r.dQ_dh_a, kHa, [&](double v) { return exchange(c.C_a, c.C_b, cf, v); });
      check(r.dQ_dh_b, kHb, [&](double v) { return exchange(c.C_a, c.C_b, cf, kHa, v); });
      check(r.dQ_dT_aw_a, kTa,
            [&](double v) { return exchange(c.C_a, c.C_b, cf, kHa, kHb, v); });
      check(r.dQ_dT_aw_b, kTb,
            [&](double v) { return exchange(c.C_a, c.C_b, cf, kHa, kHb, kTa, v); });
      // #932: dQ/d(t/k) of the layer.
      check(r.dQ_dR_wall, kLayers[0], [&](double v) {
        return wall_exchange_and_jacobian(kHa, kTa, c.C_a, kHb, kTb, c.C_b, {v}, kA, 0.0, cf);
      });
    }
  }
}
