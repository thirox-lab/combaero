// Bordbar et al. (2021) pressure-dependent WSGG emissivity (#1071): a node
// is its own set (Eq. 7), the analytic partials against central
// differences, and the clamps. Fidelity against the paper's slab results is
// python/tests/test_wsgg_bordbar2021.py.

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <stdexcept>

#include "wsgg_bordbar2021.h"

using combaero::gas_emissivity_bordbar2021;
using combaero::WSGG_PA_PER_ATM;
using combaero::wsgg_bordbar2021_gray_gases;
using combaero::WsggEmissivityResult;

TEST(WsggBordbar2021, NodeIsItsSet) {
  const double T = 1600.0, xh = 0.08, xc = 0.06, L = 0.2;
  for (double p : {0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0, 80.0}) {
    const auto g = wsgg_bordbar2021_gray_gases(p, T, xh / xc);
    const double S = (xh + xc) * p * L;
    double want = 0.0;
    for (int i = 0; i < 4; ++i) {
      want += g.a[i] * (1.0 - std::exp(-g.kappa[i] * S));
    }
    const auto e = gas_emissivity_bordbar2021(T, p * WSGG_PA_PER_ATM, xh, xc, L);
    EXPECT_NEAR(e.eps, want, 1e-12 * want) << p;
  }
}

TEST(WsggBordbar2021, PartialsMatchCentralDifferences) {
  const std::array<std::array<double, 5>, 4> points{{
      {1567.0, 5.0e5, 0.0697, 0.0613, 0.126},
      {1532.0, 10.0e5, 0.0697, 0.0613, 0.076},
      {2200.0, 3.0e6, 0.12, 0.05, 0.4},
      {900.0, 0.3e5, 0.03, 0.09, 1.5},
  }};
  for (const auto &x : points) {
    const WsggEmissivityResult e =
        gas_emissivity_bordbar2021(x[0], x[1], x[2], x[3], x[4]);
    const std::array<double, 5> an{e.deps_dT, e.deps_dP, e.deps_dx_h2o,
                                   e.deps_dx_co2, e.deps_dl_b};
    for (int i = 0; i < 5; ++i) {
      auto hi = x;
      auto lo = x;
      const double h = 1e-6 * x[i];
      hi[i] += h;
      lo[i] -= h;
      const double fd =
          (gas_emissivity_bordbar2021(hi[0], hi[1], hi[2], hi[3], hi[4]).eps -
           gas_emissivity_bordbar2021(lo[0], lo[1], lo[2], lo[3], lo[4]).eps) /
          (2.0 * h);
      EXPECT_NEAR(an[i], fd, 1e-5 * std::abs(fd) + 1e-9 / x[i]) << i;
    }
  }
}

TEST(WsggBordbar2021, Clamps) {
  const double P = 10.0 * WSGG_PA_PER_ATM;
  const auto hot = gas_emissivity_bordbar2021(3400.0, P, 0.1, 0.08, 0.3);
  EXPECT_EQ(hot.eps, gas_emissivity_bordbar2021(3000.0, P, 0.1, 0.08, 0.3).eps);
  EXPECT_EQ(hot.deps_dT, 0.0);
  EXPECT_NEAR(gas_emissivity_bordbar2021(1800.0, P, 0.2, 0.0, 0.3).eps,
              gas_emissivity_bordbar2021(1800.0, P, 0.16, 0.04, 0.3).eps,
              1e-14);
  EXPECT_EQ(gas_emissivity_bordbar2021(1800.0, P, 0.0, 0.0, 0.3).eps, 0.0);
  // 2 atm, 300 K, Mr 1: sum a_i = 1.0285, so the thick limit is clamped.
  const auto thick = gas_emissivity_bordbar2021(300.0, 2.0 * WSGG_PA_PER_ATM,
                                                0.5, 0.5, 1000.0);
  EXPECT_EQ(thick.eps, 1.0);
  EXPECT_EQ(thick.deps_dl_b, 0.0);
  EXPECT_THROW(gas_emissivity_bordbar2021(1800.0, P, 0.7, 0.4, 0.3),
               std::invalid_argument);
  EXPECT_THROW(wsgg_bordbar2021_gray_gases(3.0, 1500.0, 1.0),
               std::invalid_argument);
}
