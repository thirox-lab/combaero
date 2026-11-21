#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "composition.h"
#include "fanno_mach.h"
#include "flow_drive_floor.h"
#include "solver_interface.h"
#include "thermo.h"

// The zero-drive convention (#1193, #1209, #1256; flow_drive_floor.h): below
// kJetDriveFloor of drive the jet kernels are the floor's flow times the C1
// odd cubic h(t), with a finite slope at zero drive and no round-off
// staircase; above it, the law as before. The Fanno kernels return the
// cubic's slope at exactly zero drive.

using namespace combaero;
using combaero::solver::orifice_compressible_mdot_and_jacobian;
using combaero::solver::orifice_compressible_residuals_and_jacobian;
using combaero::solver::orifice_jobson_residuals_and_jacobian;

namespace {

// The partials are differenced on a 2 kPa reservoir: a step on 1e5 Pa cannot
// resolve a fraction of a 1e-3 Pa floor (its ulp is 1.5e-11 Pa). The
// staircase test runs at 1 bar, where main's staircase was measured.
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

TEST(ZeroDriveFloor1256, NozzleIsC1AtTheFloorAndOddThroughZero) {
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
  // One-sided differences on each side of the floor agree with J there.
  const double h = 1e-5 * kFloor;
  EXPECT_NEAR((m(kFloor) - m(kFloor - h)) / h / slope(kFloor), 1.0, 3e-4);
  EXPECT_NEAR((m(kFloor + h) - m(kFloor)) / h / slope(kFloor), 1.0, 3e-4);
  // Odd through zero to the floor's flow's own variation with its reservoir
  // (P0 is the upstream end either way, a drive apart: ~drive / P0).
  for (const double d : {1e-9 * kFloor, 0.01 * kFloor, 0.5 * kFloor}) {
    EXPECT_NEAR(m(-d) / -m(d), 1.0, 1e-5) << d;
  }
  // The slope at zero is the cubic's (3 - sigma)/2 G0 / floor, sigma ~ 1/2.
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
  // sigma is held fixed below the floor. Jobson's contraction moves the
  // exponent with Cd and, with an approach, with the static drive: measured,
  // d(m)/d(Cd) is then up to 4% off below the floor at beta 0.5 (1e-5 at
  // beta 0), d(m)/d(P1) 2e-5. The pressures that set the drive are exact.
  const std::vector<double> Y = air_Y();
  for (const double beta : {0.0, 0.5}) {
    const double P1 = beta > 0.0 ? kPt - 1e-2 * kFloor : kPt;
    const double pd = beta > 0.0 ? P1 : kPt;
    const double tol_cd = beta > 0.0 ? 5e-2 : 1e-4;
    auto r = [&](double Pt, double p1, double Tt, double Pb, double Cd) {
      return orifice_jobson_residuals_and_jacobian(0.0, Pt, p1, Tt, Y, Pb, Cd, kArea, beta);
    };
    for (const double d : kDrives) {
      for (const double sgn : {1.0, -1.0}) {
        if (d == 0.0 && sgn < 0.0) {
          continue;
        }
        const double Pb = pd - sgn * d;
        const auto res = r(kPt, P1, kTt, Pb, kCd);
        // In the drive |pd - Pb|: forward, Pb lowers it; reversed, raises it.
        const double fd_Pb =
            -sgn * d_ddrive([&](double dd) { return r(kPt, P1, kTt, pd - sgn * dd, kCd).m_dot_calc; },
                            d);
        EXPECT_TRUE(std::isfinite(res.d_mdot_dP_static_down) && res.d_mdot_dP_static_down < 0.0)
            << "beta " << beta << " drive " << sgn * d;
        EXPECT_NEAR(res.d_mdot_dP_static_down / fd_Pb, 1.0, 1e-6)
            << "beta " << beta << " drive " << sgn * d;
        if (d > 0.0) {
          const double fd_Tt =
              central([&](double t) { return r(kPt, P1, t, Pb, kCd).m_dot_calc; }, kTt, kStepT);
          const double fd_Cd =
              central([&](double c) { return r(kPt, P1, kTt, Pb, c).m_dot_calc; }, kCd, 1e-5);
          EXPECT_NEAR(res.d_mdot_dTt_up / fd_Tt, 1.0, 1e-4)
              << "beta " << beta << " drive " << sgn * d;
          EXPECT_NEAR(res.d_mdot_dCd / fd_Cd, 1.0, tol_cd)
              << "beta " << beta << " drive " << sgn * d;
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
