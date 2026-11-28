// Regression tests for #1269: shared zero-drive sub-floor law
// (include/flow_drive_floor.h; #1259) ported to the ejector C-D nozzle.
//
// Before #1269, the exit flux grew as sqrt(1 - r) (r = P_py / Pt), so
// dmdot/dP_py diverged to -infinity as r -> 1 (-2.34 kg/(s Pa) at 1 - r = 1e-14)
// and jumped to 0 at r = 1 exactly. Below kJetDriveFloor (1e-3 Pa) of drive,
// the flow now takes m(floor) * h(drive / floor) with kJetFloorExponent = 1/2:
// finite odd-symmetric slope through r = 1, exact DualN partials across the
// sub-floor band and handover, and bit-identical forward flow above the floor.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>

#include "ejector.h"
#include "flow_drive_floor.h"

using combaero::kJetDriveFloor;
using combaero::kJetFloorExponent;
using namespace combaero::solver;

namespace {

constexpr double kPt = 2.0e3;      // [Pa] resolves sub-floor FD steps cleanly (see #1256)
constexpr double kPtBar = 1.0e5;   // [Pa] 1 bar case for the 1 - r = 1e-14..1e-3 check
constexpr double kTt = 300.0;      // [K]
constexpr double kAt = 3.14e-5;    // [m^2]
constexpr double kAe = 1.0e-4;     // [m^2]
constexpr double kGamma = 1.4;
constexpr double kRGas = 287.0;
constexpr double kEta = 0.95;
constexpr double kEpsFrac = 1.0e-3;
constexpr double kFloor = kJetDriveFloor;

double central(const std::function<double(double)>& f, double x, double h) {
  return (f(x + h) - f(x - h)) / (2.0 * h);
}

// Derivative w.r.t. signed drive d = p0 - p_py, staying clear of the handover
// edges at +-kFloor (the cubic is C0 with a small slope step there).
double d_ddrive(const std::function<double(double)>& f, double d) {
  const double room = std::min(std::abs(d - kFloor), std::abs(d + kFloor));
  const double h = std::min({1e-3 * kFloor, 0.25 * room, 1e-3});
  return central(f, d, h);
}

}  // namespace

TEST(EjectorZeroDriveFloor1269, PartialsMatchCentralDifferencesAcrossBandAndHandover) {
  // Signed drives d = p0 - p_py spanning reverse above/inside floor, exact
  // r = 1 (d = 0), forward sub-floor, handover, and exact law above the floor.
  const double drives[] = {-5.0 * kFloor, -1.1 * kFloor, -0.9 * kFloor, -0.3 * kFloor,
                           -1e-3 * kFloor, 0.0,           1e-7 * kFloor, 1e-3 * kFloor,
                           0.3 * kFloor,  0.9 * kFloor,  1.1 * kFloor,  5.0 * kFloor,
                           1.0,           10.0};
  auto m_of = [&](double p0, double t0, double p_py) {
    return ejector_cd_nozzle_mass_flow(p0, t0, p_py, kAt, kAe, kGamma, kRGas, kEta, kEpsFrac);
  };

  constexpr double kRelTol = 1e-6;
  constexpr double kStepT = 1e-3;
  for (const double d : drives) {
    const double p_py = kPt - d;
    const EjectorCDNozzleJacobian j = ejector_cd_nozzle_mass_flow_and_jacobian(
        kPt, kTt, p_py, kAt, kAe, kGamma, kRGas, kEta, kEpsFrac);
    EXPECT_NEAR(j.mdot, m_of(kPt, kTt, p_py), 1e-13 * std::abs(j.mdot) + 1e-18) << "drive=" << d;

    const double fd_p_py = -d_ddrive([&](double dd) { return m_of(kPt, kTt, kPt - dd); }, d);
    const double fd_p0 = d_ddrive([&](double dd) { return m_of(kPt + (dd - d), kTt, p_py); }, d);
    const double fd_t0 = central([&](double t) { return m_of(kPt, t, p_py); }, kTt, kStepT);

    EXPECT_TRUE(std::isfinite(j.dmdot_dp_py) && j.dmdot_dp_py < 0.0) << "drive=" << d;
    EXPECT_TRUE(std::isfinite(j.dmdot_dp0) && j.dmdot_dp0 > 0.0) << "drive=" << d;
    EXPECT_NEAR(j.dmdot_dp_py / fd_p_py, 1.0, kRelTol) << "drive=" << d;
    EXPECT_NEAR(j.dmdot_dp0 / fd_p0, 1.0, kRelTol) << "drive=" << d;
    if (d == 0.0) {
      EXPECT_DOUBLE_EQ(j.mdot, 0.0);
      EXPECT_DOUBLE_EQ(j.dmdot_dt0, 0.0);
    } else {
      EXPECT_NEAR(j.dmdot_dt0 / fd_t0, 1.0, kRelTol) << "drive=" << d;
    }
  }
}

TEST(EjectorZeroDriveFloor1269, FiniteOddSymmetricSlopeAtZeroDriveAndNoBlowUp) {
  // Gaps 1 - r from 1e-3 down to 1e-14 and r = 1 exactly at 1 bar (#1269):
  // on main dmdot/dp_py blew up to -2.34 kg/(s Pa) at 1 - r = 1e-14 and jumped
  // to 0 at r = 1. With the sub-floor cubic, the slope is bounded by the zero-
  // drive cubic slope 1.25 * m(floor) / floor and matches both sides at r = 1.
  const double m_floor = ejector_cd_nozzle_mass_flow(kPtBar, kTt, kPtBar - kFloor, kAt, kAe, kGamma,
                                                     kRGas, kEta, kEpsFrac);
  const double slope_zero_expected = -0.5 * (3.0 - kJetFloorExponent) * m_floor / kFloor;

  const EjectorCDNozzleJacobian at_zero = ejector_cd_nozzle_mass_flow_and_jacobian(
      kPtBar, kTt, kPtBar, kAt, kAe, kGamma, kRGas, kEta, kEpsFrac);
  EXPECT_DOUBLE_EQ(at_zero.mdot, 0.0);
  EXPECT_NEAR(at_zero.dmdot_dp_py / slope_zero_expected, 1.0, 1e-12);
  EXPECT_NEAR(at_zero.dmdot_dp0 / (-slope_zero_expected), 1.0, 1e-12);

  // Gaps 1 - r = 1e-3 down to 1e-14 (forward) and r = 1 + 1e-14 (reverse).
  for (const double gap : {1e-3, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14}) {
    const double p_fwd = kPtBar * (1.0 - gap);
    const EjectorCDNozzleJacobian jf = ejector_cd_nozzle_mass_flow_and_jacobian(
        kPtBar, kTt, p_fwd, kAt, kAe, kGamma, kRGas, kEta, kEpsFrac);
    EXPECT_TRUE(std::isfinite(jf.dmdot_dp_py) && jf.dmdot_dp_py < 0.0) << "gap=" << gap;
    EXPECT_LE(std::abs(jf.dmdot_dp_py), 1.01 * std::abs(slope_zero_expected)) << "gap=" << gap;
  }
  const EjectorCDNozzleJacobian rev_near = ejector_cd_nozzle_mass_flow_and_jacobian(
      kPtBar, kTt, kPtBar * (1.0 + 1e-14), kAt, kAe, kGamma, kRGas, kEta, kEpsFrac);
  EXPECT_NEAR(rev_near.dmdot_dp_py / at_zero.dmdot_dp_py, 1.0, 1e-6);

  // Continuous value across the handover at +-kFloor, and small slope step.
  const double e = 1e-7 * kFloor;
  auto eval = [&](double d) {
    return ejector_cd_nozzle_mass_flow_and_jacobian(kPt, kTt, kPt - d, kAt, kAe, kGamma, kRGas,
                                                    kEta, kEpsFrac);
  };
  const auto lo = eval(kFloor - e);
  const auto hi = eval(kFloor + e);
  EXPECT_NEAR(lo.mdot / hi.mdot, 1.0, 1e-6);
  EXPECT_NEAR(lo.dmdot_dp_py / hi.dmdot_dp_py, 1.0, 1e-3);

  // Odd symmetry through zero inside the sub-floor band.
  for (const double d : {1e-9 * kFloor, 1e-3 * kFloor, 0.5 * kFloor}) {
    const double mf = eval(d).mdot;
    const double mr = eval(-d).mdot;
    EXPECT_NEAR(mr / -mf, 1.0, 1e-5) << "d=" << d;
  }
}

TEST(EjectorZeroDriveFloor1269, ForwardStatesAboveFloorAreBitIdenticalToRawFormula) {
  // Independent evaluation of the raw C-D nozzle formula (no sub-floor branch)
  // confirms bit-identity for every forward state at drive >= kJetDriveFloor.
  auto raw_ref = [&](double p0, double t0, double p_py) {
    const double gm1 = kGamma - 1.0;
    const double gp1 = kGamma + 1.0;
    const double cap = ejector_choked_mass_flow(p0, t0, kAt, kGamma, kRGas, kEta);
    const double x = (p_py - p0) / p0;
    double mach = std::sqrt(2.0 / gm1 * std::expm1(-gm1 / kGamma * std::log1p(x)));
    mach = std::min(mach, 1.0);
    const double flux = mach * std::pow(1.0 + 0.5 * gm1 * mach * mach, -gp1 / (2.0 * gm1));
    const double exit_flux =
        kAe * p0 / std::sqrt(t0) * std::sqrt(kGamma / kRGas) * flux * std::sqrt(kEta);
    const double eps = kEpsFrac * cap;
    const double diff = exit_flux - cap;
    return 0.5 * (exit_flux + cap - std::sqrt(diff * diff + eps * eps));
  };

  for (const double d : {kFloor, 2.0 * kFloor, 1.0, 100.0, 50000.0}) {
    const double p_py = kPtBar - d;
    const double actual =
        ejector_cd_nozzle_mass_flow(kPtBar, kTt, p_py, kAt, kAe, kGamma, kRGas, kEta, kEpsFrac);
    EXPECT_DOUBLE_EQ(actual, raw_ref(kPtBar, kTt, p_py)) << "drive=" << d;
  }
}

TEST(EjectorZeroDriveFloor1269, WholeElementGammaAndRGasPartialsMatchFDInSubFloorBand) {
  // Verify the DualN gamma and r_gas propagation through the sub-floor cubic
  // inside ejector_element_residuals_and_jacobian (p_py within 0.4 * kFloor of p_g).
  const double p_g = kPt;
  const double t_g = kTt;
  const double p_e = kPt - 10.0;
  const double t_e = kTt;
  const double p_out = kPt - 5.0;
  const double p_py = p_g - 0.4 * kFloor;
  const double mp = 1e-6, ms = 1e-6, mdot_out = 2e-6;
  const EjectorGeometry geom{kAe / kAt, 6.44};
  const double A_s = (geom.area_ratio_mix - geom.area_ratio_nozzle) * kAt;
  const double rec = 0.85;

  auto eval = [&](double gam, double rg) {
    return ejector_element_residuals_and_jacobian(
        mp, ms, mdot_out, p_g, t_g, p_e, t_e, p_out, p_py, geom, kAt, kAe, A_s, gam, rg, kEta,
        0.85, rec, kEpsFrac, 0.90, 0.98, 0.90, 0.98);
  };
  const EjectorElementResidualJacobian sys = eval(kGamma, kRGas);
  const double hg = 1e-6 * kGamma;
  const double hr = 1e-6 * kRGas;
  const double fd_g =
      (eval(kGamma + hg, kRGas).residuals[0] - eval(kGamma - hg, kRGas).residuals[0]) / (2.0 * hg);
  const double fd_r =
      (eval(kGamma, kRGas + hr).residuals[0] - eval(kGamma, kRGas - hr).residuals[0]) / (2.0 * hr);
  EXPECT_NEAR(sys.d_dgamma[0], fd_g, 1e-5 * std::abs(fd_g) + 1e-12);
  EXPECT_NEAR(sys.d_dr_gas[0], fd_r, 1e-5 * std::abs(fd_r) + 1e-12);
}
