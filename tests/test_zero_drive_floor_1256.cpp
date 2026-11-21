#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "composition.h"
#include "fanno_mach.h"
#include "flow_drive_floor.h"
#include "solver_interface.h"
#include "thermo.h"

// The zero-drive convention (#1193, #1209, #1256; flow_drive_floor.h): below
// kJetDriveFloor of drive the jet kernels are the floor's flow times the odd
// cubic h(t) with the jet exponent 1/2: a finite slope at zero drive, exact
// partials, no round-off staircase below the floor; above it, the law as
// before (including its own round-off: the nozzle's staircase reaches ~8.7e-7
// Pa at 1 bar, ~9.7e-4 Pa at 100 bar and past the floor at 1000 bar, #1261).
// The Fanno kernels return their cubic's slope at exactly zero drive.

using namespace combaero;
using combaero::solver::orifice_compressible_mdot_and_jacobian;
using combaero::solver::orifice_compressible_residuals_and_jacobian;
using combaero::solver::orifice_jobson_residuals_and_jacobian;

namespace {

// The partials are differenced on a 2 kPa reservoir: a step on 1e5 Pa cannot
// resolve a fraction of a 1e-3 Pa floor (its ulp is 1.5e-11 Pa). The
// staircase test runs at 1 bar: main's staircase there, over 0 .. 4e-9 Pa.
constexpr double kPt = 2.0e3;
constexpr double kPtBar = 1.01e5;
constexpr double kTt = 300.0;
constexpr double kCd = 0.6;
constexpr double kArea = 2e-5;
constexpr double kFloor = kJetDriveFloor;

std::vector<double> air_Y() {
  std::vector<double> X(num_species(), 0.0);
  X[species_index_from_name("N2")] = 0.7808;
  X[species_index_from_name("O2")] = 0.2095;
  X[species_index_from_name("AR")] = 0.0097;
  return mole_to_mass(X);
}

std::vector<double> air_X() { return mass_to_mole(air_Y()); }

double central(const std::function<double(double)>& f, double x, double h) {
  return (f(x + h) - f(x - h)) / (2.0 * h);
}

// d f / d(drive) for f a function of the drive: one-sided forward over 1e-3
// of the floor within 2e-3 floors of zero (the jet's reverse side is another
// state's law), else central over 1e-3 floors or a quarter of the distance to
// the floor, whichever is less (the cubic is C1 but not C2 there).
double d_ddrive(const std::function<double(double)>& f, double d) {
  if (d < 2e-3 * kFloor) {
    const double h = 1e-3 * kFloor;
    return (f(d + h) - f(d)) / h;
  }
  const double room = std::abs(d - kFloor);
  return central(f, d, std::min({1e-3 * kFloor, 0.25 * room, 1e-3}));
}

// Drives at which the partials are checked: zero, deep and just inside the
// floor, and on the exact law above it.
const double kDrives[] = {0.0,           1e-7 * kFloor, 1e-3 * kFloor, 0.3 * kFloor,
                          0.9 * kFloor,  1.1 * kFloor,  5.0 * kFloor,  1.0,
                          300.0};

// The flow's temperature (and, from the static state, stagnation pressure)
// partials are differenced over 0.1 K (5 Pa): nozzle_flow's own root carries
// noise of ~1e-4 of a 1e-3 K difference at a 1 Pa drive, and ~2e-6 of a 0.1 K
// one. Every column is held to kTol: the raw law's own differences (the
// nozzle's root, above the floor and unchanged from main) scatter to ~1e-5
// on a 2 kPa reservoir.
constexpr double kTol = 5e-5;
constexpr double kStepT = 0.1;
constexpr double kStepPt = 5.0;

}  // namespace

TEST(ZeroDriveFloor1256, NozzlePartialsMatchDifferences) {
  const std::vector<double> X = air_X();
  auto m = [&](double Pt, double Pb, double Tt) {
    return std::get<0>(orifice_compressible_mdot_and_jacobian(Tt, Pt, Pb, X, kCd, kArea, 0.0));
  };
  for (const double d : kDrives) {
    const double Pb = kPt - d;
    const auto [mdot, dPt, dPb, dTt] =
        orifice_compressible_mdot_and_jacobian(kTt, kPt, Pb, X, kCd, kArea, 0.0);
    // The back pressure lowers the drive; the stagnation pressure raises it
    // and the floor's flow with it.
    const double fd_Pb = -d_ddrive([&](double dd) { return m(kPt, kPt - dd, kTt); }, d);
    const double fd_Pt = d_ddrive([&](double dd) { return m(kPt + (dd - d), Pb, kTt); }, d);
    const double fd_Tt = central([&](double t) { return m(kPt, Pb, t); }, kTt, kStepT);
    EXPECT_TRUE(std::isfinite(dPb) && dPb < 0.0) << d;
    EXPECT_NEAR(dPb / fd_Pb, 1.0, kTol) << "drive " << d;
    EXPECT_NEAR(dPt / fd_Pt, 1.0, kTol) << "drive " << d;
    if (d > 0.0) {
      EXPECT_NEAR(dTt / fd_Tt, 1.0, kTol) << "drive " << d;
    } else {
      EXPECT_EQ(mdot, 0.0);
      EXPECT_EQ(dTt, 0.0);
    }
  }
}

TEST(ZeroDriveFloor1256, NozzleIsContinuousAtTheFloorAndOddThroughZero) {
  const std::vector<double> X = air_X();
  auto m = [&](double drive) {
    return std::get<0>(
        orifice_compressible_mdot_and_jacobian(kTt, kPt, kPt - drive, X, kCd, kArea, 0.0));
  };
  auto slope = [&](double drive) {
    return -std::get<2>(
        orifice_compressible_mdot_and_jacobian(kTt, kPt, kPt - drive, X, kCd, kArea, 0.0));
  };
  const double e = 1e-7 * kFloor;
  EXPECT_NEAR(m(kFloor - e) / m(kFloor + e), 1.0, 1e-6);
  EXPECT_NEAR(slope(kFloor - e) / slope(kFloor + e), 1.0, 1e-5);
  // (Each side's slope against differences of the value: the partials test,
  // at 0.9 and 1.1 floors. A one-sided difference straddling the floor over
  // 1e-8 Pa reads the nozzle root's own scatter, 4e-4 on MinGW.)
  // Odd through zero to the floor's flow's own variation with its reservoir
  // (P0 is the upstream end either way, a drive apart: ~drive / P0).
  for (const double d : {1e-9 * kFloor, 0.01 * kFloor, 0.5 * kFloor}) {
    EXPECT_NEAR(m(-d) / -m(d), 1.0, 1e-5) << d;
  }
  // The slope at zero is the cubic's (3 - 1/2)/2 G0 / floor.
  EXPECT_NEAR(slope(0.0) / (1.25 * m(kFloor) / kFloor), 1.0, 1e-4);
}

TEST(ZeroDriveFloor1256, NozzleHasNoRoundOffStaircaseNearZero) {
  // Main: over 0 .. 4e-9 Pa of drive the flow took 28 distinct values with 14
  // reversals of its increments, and its slope read 0 below ~1.5e-10 Pa.
  const std::vector<double> X = air_X();
  double prev = -1.0;
  for (int k = 0; k <= 40; ++k) {
    const double d = 1e-10 * k;
    const double mdot = std::get<0>(
        orifice_compressible_mdot_and_jacobian(kTt, kPtBar, kPtBar - d, X, kCd, kArea, 0.0));
    EXPECT_GE(mdot, prev) << d;
    prev = mdot;
  }
}

TEST(ZeroDriveFloor1256, ResidualKernelPartialsFromTheStaticStateToo) {
  // With E (beta 0.5) the nozzle expands from the approach static state.
  const std::vector<double> Y = air_Y();
  const double beta = 0.5;
  const double P1 = kPt - 1e-2 * kFloor;
  auto r = [&](double Pt, double p1, double Tt, double Pb) {
    return orifice_compressible_residuals_and_jacobian(0.0, Pt, p1, Tt, Y, Pb, kCd, kArea, beta);
  };
  for (const double d : kDrives) {
    const double Pb = P1 - d;
    const auto res = r(kPt, P1, kTt, Pb);
    const double fd_Pb = -d_ddrive([&](double dd) { return r(kPt, P1, kTt, P1 - dd).m_dot_calc; }, d);
    const double fd_P1 =
        d_ddrive([&](double dd) { return r(kPt, P1 + (dd - d), kTt, Pb).m_dot_calc; }, d);
    EXPECT_NEAR(res.d_mdot_dP_static_down / fd_Pb, 1.0, kTol) << "drive " << d;
    EXPECT_NEAR(res.d_mdot_dP_static_up / fd_P1, 1.0, kTol) << "drive " << d;
    if (d > 0.0) {
      const double fd_Pt =
          central([&](double p) { return r(p, P1, kTt, Pb).m_dot_calc; }, kPt, kStepPt);
      const double fd_Tt =
          central([&](double t) { return r(kPt, P1, t, Pb).m_dot_calc; }, kTt, kStepT);
      EXPECT_NEAR(res.d_mdot_dP_total_up / fd_Pt, 1.0, 1e-4) << "drive " << d;
      EXPECT_NEAR(res.d_mdot_dTt_up / fd_Tt, 1.0, kTol) << "drive " << d;
    }
  }
}

TEST(ZeroDriveFloor1256, JobsonPartialsMatchDifferences) {
  // Every column, with and without an approach (beta 0.5, 0.7: the drive is
  // the approach static P1, and an approach head q = Pt - P1 of 0, 1e-2 and 1
  // floor moves Jobson's own exponent at the floor from 0.5 to 0.8). The
  // cubic's exponent is the constant 1/2, so the band's partials are exact;
  // a matched, held exponent left dm/dPt 22-27% and dm/dP1 7-15% off there.
  const std::vector<double> Y = air_Y();
  for (const double beta : {0.0, 0.5, 0.7}) {
    for (const double qf : {0.0, 1e-2, 1.0}) {
      if (beta == 0.0 && qf > 0.0) {
        continue;  // without E the drive is Pt itself: no approach head
      }
      const double P1 = kPt - qf * kFloor;
      const double pd = beta > 0.0 ? P1 : kPt;
      auto r = [&](double Pt, double p1, double Tt, double Pb, double Cd) {
        return orifice_jobson_residuals_and_jacobian(0.0, Pt, p1, Tt, Y, Pb, Cd, kArea, beta);
      };
      auto m = [&](double Pt, double p1, double Pb) { return r(Pt, p1, kTt, Pb, kCd).m_dot_calc; };
      for (const double d : kDrives) {
        for (const double sgn : {1.0, -1.0}) {
          if (d == 0.0 && sgn < 0.0) {
            continue;
          }
          const double Pb = pd - sgn * d;
          const auto res = r(kPt, P1, kTt, Pb, kCd);
          const std::string at = "beta " + std::to_string(beta) + " q/floor " +
                                 std::to_string(qf) + " drive " + std::to_string(sgn * d);
          // In the drive |pd - Pb|: forward, Pb lowers it and pd raises it;
          // reversed, the other way round.
          const double fd_Pb = -sgn * d_ddrive([&](double dd) { return m(kPt, P1, pd - sgn * dd); }, d);
          EXPECT_TRUE(std::isfinite(res.d_mdot_dP_static_down) && res.d_mdot_dP_static_down < 0.0)
              << at;
          EXPECT_NEAR(res.d_mdot_dP_static_down / fd_Pb, 1.0, kTol) << at;
          if (beta > 0.0) {
            // Pt and P1 also move q, and the flow curves in q on the floor's
            // own scale: differenced over 1e-5 floors (forward in the drive
            // near zero), to 2e-4 (the 2 kPa reservoir's ulp over that step).
            const double hq = 1e-5 * kFloor;
            const double fd_P1 =
                (d < 2e-3 * kFloor)
                    ? (m(kPt, P1 + sgn * hq, Pb) - res.m_dot_calc) / (sgn * hq)
                    : central([&](double p) { return m(kPt, p, Pb); }, P1, hq);
            const double fd_Pt = central([&](double p) { return m(p, P1, Pb); }, kPt, hq);
            EXPECT_NEAR(res.d_mdot_dP_static_up / fd_P1, 1.0, 2e-4) << at;
            EXPECT_NEAR(res.d_mdot_dP_total_up, fd_Pt, 2e-4 * std::abs(fd_P1)) << at;
          } else {
            const double fd_Pt =
                sgn * d_ddrive([&](double dd) { return m(kPt + sgn * (dd - d), P1, Pb); }, d);
            EXPECT_NEAR(res.d_mdot_dP_total_up / fd_Pt, 1.0, kTol) << at;
          }
          if (d > 0.0) {
            const double fd_Tt =
                central([&](double t) { return r(kPt, P1, t, Pb, kCd).m_dot_calc; }, kTt, kStepT);
            const double fd_Cd =
                central([&](double c) { return r(kPt, P1, kTt, Pb, c).m_dot_calc; }, kCd, 1e-5);
            EXPECT_NEAR(res.d_mdot_dTt_up / fd_Tt, 1.0, 1e-4) << at;
            EXPECT_NEAR(res.d_mdot_dCd / fd_Cd, 1.0, kTol) << at;
          }
        }
      }
    }
  }
}

TEST(ZeroDriveFloor1209, FannoSlopeAtExactlyZeroDriveIsTheCubics) {
  const std::vector<double> X = air_X();
  const std::string model = "haaland";
  for (const bool heated : {false, true}) {
    for (const bool exit_total : {true, false}) {
      auto flow = [&](double Pt0, double P_target) {
        if (heated) {
          return fanno_rayleigh_channel_flow(Pt0, kTt, kTt + 50.0, X, P_target, exit_total, 1.0,
                                             0.01, 0.0, model);
        }
        return fanno_channel_flow(Pt0, kTt, X, P_target, exit_total, 1.0, 0.01, 0.0, model);
      };
      const FannoChannelFlow f0 = flow(kPt, kPt);
      EXPECT_EQ(f0.G, 0.0);
      EXPECT_EQ(f0.M_in, 0.0);
      EXPECT_EQ(f0.M_exit, 0.0);
      EXPECT_EQ(f0.M_exit_root, 0.0);
      // One-sided: the flux is 0 at and below zero drive.
      const double h = 1e-4;
      const double fd_P = (flow(kPt, kPt - h).G - f0.G) / h;  // = -dG/dP_target
      const double fd_Pt = (flow(kPt + h, kPt).G - f0.G) / h;
      EXPECT_GT(-f0.dG_dP_target, 0.0);
      EXPECT_NEAR(-f0.dG_dP_target / fd_P, 1.0, 1e-4) << heated << exit_total;
      EXPECT_NEAR(f0.dG_dPt0 / fd_Pt, 1.0, 1e-4) << heated << exit_total;
      EXPECT_EQ(f0.dG_dTt0, 0.0);
      // Without derivatives, nothing is solved: the value alone, 0.
      const FannoChannelFlow v = heated ? fanno_rayleigh_channel_flow(
                                              kPt, kTt, kTt + 50.0, X, kPt, exit_total, 1.0, 0.01,
                                              0.0, model, 1.0, -1.0, {}, false)
                                        : fanno_channel_flow(kPt, kTt, X, kPt, exit_total, 1.0,
                                                             0.01, 0.0, model, 1.0, -1.0, {},
                                                             false);
      EXPECT_EQ(v.G, 0.0);
      EXPECT_EQ(v.dG_dP_target, 0.0);
    }
  }
}
