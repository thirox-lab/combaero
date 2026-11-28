// Regression and analytic-vs-FD Jacobian tests for #1128:
//   - solver::discharge_total_pressure (subsonic, choked, C1 across sonic,
//     zero flow, incompressible, multi-species Y)
//   - area_change_flow / area_change_element_flow with exit_head_lost = true
//     (sharp/conical, expansion/contraction, compressible/incompressible,
//     low-Re exit-face choke, and bit-identity when exit_head_lost = false).

#include "../include/area_change.h"
#include "../include/composition.h"
#include "../include/humidair.h"
#include "../include/solver_interface.h"
#include "../include/thermo.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace combaero {
namespace {

std::vector<double> fuel_mix_mass() {
  std::vector<double> X(num_species(), 0.0);
  X[species_index_from_name("N2")] = 0.72;
  X[species_index_from_name("O2")] = 0.18;
  X[species_index_from_name("CH4")] = 0.07;
  X[species_index_from_name("CO2")] = 0.03;
  return mole_to_mass(X);
}

TEST(ExitHeadLost1128, DischargeTotalPressurePartialsMatchFD) {
  const std::vector<double> Y = fuel_mix_mass();
  const double P = 1.0e5;
  const double Tt = 850.0;

  for (const bool inc : {false, true}) {
    for (const double G : {0.0, 80.0, 220.0, 350.0}) {
      const solver::DischargeTotal r = solver::discharge_total_pressure(P, Tt, G, Y, inc);
      auto pt_of = [&](double p, double t, double g, const std::vector<double>& y) {
        return solver::discharge_total_pressure(p, t, g, y, inc).Pt;
      };

      const double hP = 1.0;
      const double fd_P = (pt_of(P + hP, Tt, G, Y) - pt_of(P - hP, Tt, G, Y)) / (2.0 * hP);
      EXPECT_NEAR(r.dPt_dP, fd_P, 1e-6 * std::abs(fd_P) + 1e-9) << "inc=" << inc << " G=" << G;

      const double hT = 1e-2;
      const double fd_T = (pt_of(P, Tt + hT, G, Y) - pt_of(P, Tt - hT, G, Y)) / (2.0 * hT);
      EXPECT_NEAR(r.dPt_dTt, fd_T, 1e-6 * std::abs(fd_T) + 1e-9) << "inc=" << inc << " G=" << G;

      const double hG = 1e-3;
      const double fd_G = (pt_of(P, Tt, G + hG, Y) - pt_of(P, Tt, G - hG, Y)) / (2.0 * hG);
      EXPECT_NEAR(r.dPt_dG, fd_G, 1e-6 * std::abs(fd_G) + 1e-6) << "inc=" << inc << " G=" << G;

      const double hY = 1e-6;
      for (std::size_t k = 0; k < Y.size(); ++k) {
        if (!(Y[k] > 0.0)) {
          continue;
        }
        std::vector<double> yp = Y;
        std::vector<double> ym = Y;
        yp[k] += hY;
        ym[k] -= hY;
        const double fd_Y = (pt_of(P, Tt, G, yp) - pt_of(P, Tt, G, ym)) / (2.0 * hY);
        EXPECT_NEAR(r.dPt_dY[k], fd_Y, 1e-5 * std::abs(fd_Y) + 1e-6 * r.Pt)
            << "inc=" << inc << " G=" << G << " k=" << k;
      }
    }
  }
}

TEST(ExitHeadLost1128, DischargeTotalPressureIsC1AcrossSonicChoke) {
  const std::vector<double> Y = fuel_mix_mass();
  const double Tt = 600.0;
  const double G = 200.0;
  const double P_star = solver::sonic_static_pressure(Tt, G, Y).P;
  const double eps = 1e-6 * P_star;
  const solver::DischargeTotal sub = solver::discharge_total_pressure(P_star + eps, Tt, G, Y, false);
  const solver::DischargeTotal chk = solver::discharge_total_pressure(P_star - eps, Tt, G, Y, false);
  EXPECT_FALSE(sub.choked);
  EXPECT_TRUE(chk.choked);
  EXPECT_NEAR(sub.Pt, chk.Pt, 1e-10 * chk.Pt);
  EXPECT_LT(std::abs(sub.dPt_dP), 1e-5);
  EXPECT_DOUBLE_EQ(chk.dPt_dP, 0.0);
  EXPECT_NEAR(sub.dPt_dG, chk.dPt_dG, 1e-5 * chk.dPt_dG);
  EXPECT_NEAR(sub.dPt_dTt, chk.dPt_dTt, 1e-5 * std::abs(chk.dPt_dTt));
}

TEST(ExitHeadLost1128, AreaChangeFlowExitHeadLostPartialsAndLowReExitChoke) {
  const std::vector<double> Y_a = fuel_mix_mass();
  const std::vector<double> Y_b = mole_to_mass(dry_air());
  const double Pt_a = 3.0e5;
  const double Tt_a = 420.0;
  const double Tt_b = 600.0;

  for (const char* model : {"sharp", "conical"}) {
    for (const bool inc : {false, true}) {
      for (const double drive : {-50.0, -0.5, 0.0, 0.5, 50.0, 2.0e4}) {
        const double Pt_b = Pt_a - drive;
        const AreaChangeElementFlow f = area_change_element_flow(
            Pt_a, Tt_a, Y_a, Pt_b, Tt_b, Y_b, 1e-3, 4e-3, model, 0.08, 0.0, inc, true, true, false);
        auto m_of = [&](double pa, double ta, double pb, double tb) {
          return area_change_element_flow(pa, ta, Y_a, pb, tb, Y_b, 1e-3, 4e-3, model, 0.08, 0.0,
                                          inc, false, true, false)
              .m_dot;
        };
        const double hp = (std::abs(drive) < 1.0) ? 1e-4 : 1e-2;
        const double ht = 1e-2;
        const double scale = std::abs(f.m_dot) + std::abs(f.dm_dPt_a) + 1e-12;
        const double fd_pa = (m_of(Pt_a + hp, Tt_a, Pt_b, Tt_b) - m_of(Pt_a - hp, Tt_a, Pt_b, Tt_b)) / (2.0 * hp);
        const double fd_pb = (m_of(Pt_a, Tt_a, Pt_b + hp, Tt_b) - m_of(Pt_a, Tt_a, Pt_b - hp, Tt_b)) / (2.0 * hp);
        const double fd_ta = (m_of(Pt_a, Tt_a + ht, Pt_b, Tt_b) - m_of(Pt_a, Tt_a - ht, Pt_b, Tt_b)) / (2.0 * ht);
        const double fd_tb = (m_of(Pt_a, Tt_a, Pt_b, Tt_b + ht) - m_of(Pt_a, Tt_a, Pt_b, Tt_b - ht)) / (2.0 * ht);
        EXPECT_NEAR(f.dm_dPt_a, fd_pa, 1e-4 * std::abs(fd_pa) + 1e-8 * scale);
        EXPECT_NEAR(f.dm_dPt_b, fd_pb, 1e-4 * std::abs(fd_pb) + 1e-8 * scale);
        EXPECT_NEAR(f.dm_dTt_a, fd_ta, 1e-4 * std::abs(fd_ta) + 1e-8 * scale);
        EXPECT_NEAR(f.dm_dTt_b, fd_tb, 1e-4 * std::abs(fd_tb) + 1e-8 * scale);
      }
    }
  }

  // Low-Re exit-face choke (M_exit = 1 reached before M_in = 1).
  const double Pt0 = 1.0e4;
  const double Tt0 = 1000.0;
  const AreaChangeFlow fc = area_change_flow(Pt0, Tt0, Y_b, 1.0e2, 1e-8, 4e-8, "sharp", 0.0, 0.0,
                                             false, true, true);
  EXPECT_TRUE(fc.choked);
  EXPECT_LT(fc.mach, 1.0);
  EXPECT_DOUBLE_EQ(fc.dm_dPt_target, 0.0);
  const double fd_pt =
      (area_change_flow(Pt0 + 1.0, Tt0, Y_b, 1.0e2, 1e-8, 4e-8, "sharp", 0.0, 0.0, false, false, true).m_dot -
       area_change_flow(Pt0 - 1.0, Tt0, Y_b, 1.0e2, 1e-8, 4e-8, "sharp", 0.0, 0.0, false, false, true).m_dot) /
      2.0;
  EXPECT_NEAR(fc.dm_dPt0, fd_pt, 1e-5 * std::abs(fd_pt));
}

}  // namespace
}  // namespace combaero
