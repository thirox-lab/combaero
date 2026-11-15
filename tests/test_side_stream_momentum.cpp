#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "composition.h"
#include "compressible.h"
#include "side_stream_momentum.h"
#include "stagnation.h"
#include "tee_junction.h"
#include "thermo.h"

using combaero::chamber_merge_face_state;
using combaero::jet_impulse;

namespace {
std::vector<double> air() {
  std::vector<double> X(combaero::num_species(), 0.0);
  X[combaero::species_index_from_name("N2")] = 0.7808;
  X[combaero::species_index_from_name("O2")] = 0.2095;
  X[combaero::species_index_from_name("AR")] = 0.0097;
  return X;
}
}  // namespace

// ---- merge chamber face ----------------------------------------------------

TEST(ChamberMergeFace, NoSideStreamsIsTheChamberItself) {
  const auto X = air();
  const double m = 0.5, P = 2.0e5, T = 900.0, A = 0.005;
  const auto r = chamber_merge_face_state(m, T, X, m, P, T, X, 0.0, A);
  EXPECT_NEAR(r.P_face, P, 1e-9 * P);
  // T is a stagnation temperature (#578): the face's static state is the
  // chamber's, from the same energy solve.
  const auto st = combaero::chamber_static_state(T, P, m / A, X);
  const double rho = combaero::density(st.T, P, X);
  const double M = m / (rho * A) / combaero::speed_of_sound(st.T, X);
  EXPECT_NEAR(r.Pt_face, combaero::Pt_from_static(P, st.T, M, X), 1e-6 * P);
  EXPECT_NEAR(r.M_face, M, 1e-9);
  EXPECT_FALSE(r.choked);
}

TEST(ChamberMergeFace, TheImpulseBalanceHoldsAtTheFaceDensity) {
  // Mach ~0.4 at the outlet: compressible, and exact by construction.
  const auto X = air();
  const double mm = 1.0, mo = 1.15, J = 25.0, P = 1.5e5, T = 1100.0, Tm = 1200.0;
  const double A = 0.01;
  const auto r = chamber_merge_face_state(mm, Tm, X, mo, P, T, X, J, A);
  // Both temperatures are stagnation (#578): densities at the static states.
  const double Tf = combaero::chamber_static_state(Tm, r.P_face, mm / A, X).T;
  const double Ts = combaero::chamber_static_state(T, P, mo / A, X).T;
  const double rho_f = combaero::density(Tf, r.P_face, X);
  const double rho = combaero::density(Ts, P, X);
  // Energy at the face: h(T0_main) = h(T_f) + u_f^2 / 2.
  const double u_f = mm / (rho_f * A);
  EXPECT_NEAR(combaero::h_mass(Tm, X) - combaero::h_mass(Tf, X), 0.5 * u_f * u_f,
              1e-8 * 0.5 * u_f * u_f);
  const double lhs = r.P_face * A + mm * mm / (rho_f * A) + J;
  const double rhs = P * A + mo * mo / (rho * A);
  EXPECT_NEAR(lhs, rhs, 1e-9 * rhs);
  EXPECT_GT(r.M_face, 0.2);
}

TEST(ChamberMergeFace, LowMachIsTheIncompressibleMerge) {
  const auto X = air();
  const double mm = 0.03, mo = 0.036, J = 0.05, P = 1.0e5, T = 300.0, A = 0.01;
  const auto r = chamber_merge_face_state(mm, T, X, mo, P, T, X, J, A);
  const double rho = combaero::density(T, P, X);
  const double incompressible = (mo * mo - mm * mm) / (rho * A * A) - J / A;
  EXPECT_NEAR(r.P_face - P, incompressible, 1e-3 * std::abs(incompressible));
}

TEST(ChamberMergeFace, DerivativesMatchCentralDifferences) {
  const auto X = air();
  // T = 1000 K would straddle the NASA-polynomial breakpoint, where h(T)
  // jumps slightly (#549) and a central difference picks the jump up.
  const double mm = 0.8, mo = 0.95, J = 12.0, P = 1.6e5, T = 1040.0, Tm = 1150.0;
  const double A = 0.01;
  const auto r = chamber_merge_face_state(mm, Tm, X, mo, P, T, X, J, A);
  auto f = [&](double a, double b, double c, double d, double e2, double g) {
    return chamber_merge_face_state(a, b, X, c, d, e2, X, g, A);
  };
  struct Probe {
    double dPf, dPtf, h;
    int which;
  };
  const Probe probes[] = {
      {r.dPf_dm_main, r.dPtf_dm_main, 1e-6, 0}, {r.dPf_dT_main, r.dPtf_dT_main, 1e-3, 1},
      {r.dPf_dm_out, r.dPtf_dm_out, 1e-6, 2},   {r.dPf_dP, r.dPtf_dP, 1.0, 3},
      {r.dPf_dT, r.dPtf_dT, 1e-3, 4},           {r.dPf_dJ, r.dPtf_dJ, 1e-4, 5}};
  for (const auto& p : probes) {
    double up[6] = {mm, Tm, mo, P, T, J}, dn[6] = {mm, Tm, mo, P, T, J};
    up[p.which] += p.h;
    dn[p.which] -= p.h;
    const auto a = f(up[0], up[1], up[2], up[3], up[4], up[5]);
    const auto b = f(dn[0], dn[1], dn[2], dn[3], dn[4], dn[5]);
    const double fdP = (a.P_face - b.P_face) / (2 * p.h);
    const double fdPt = (a.Pt_face - b.Pt_face) / (2 * p.h);
    EXPECT_NEAR(p.dPf, fdP, 1e-4 * std::abs(fdP) + 1e-6) << p.which;
    EXPECT_NEAR(p.dPtf, fdPt, 1e-4 * std::abs(fdPt) + 1e-6) << p.which;
  }
}

// ---- incompressible face (#588) ---------------------------------------------

TEST(ChamberMergeFaceIncompressible, IsTheImpulseAtChamberPressureDensities) {
  const auto X = air();
  const double mm = 0.8, mo = 0.95, J = 12.0, P = 1.6e5, T = 1040.0, Tm = 1150.0;
  const double A = 0.01;
  const auto r = chamber_merge_face_state(mm, Tm, X, mo, P, T, X, J, A, true);
  const double rho = combaero::density(T, P, X);
  const double rho_m = combaero::density(Tm, P, X);
  const double Pf = P + mo * mo / (rho * A * A) - J / A - mm * mm / (rho_m * A * A);
  EXPECT_NEAR(r.P_face, Pf, 1e-9 * P);
  EXPECT_NEAR(r.Pt_face - r.P_face, 0.5 * mm * mm / (rho_m * A * A), 1e-9 * P);
  EXPECT_FALSE(r.choked);
}

TEST(ChamberMergeFaceIncompressible, ApproachesTheCompressibleFaceAsMachSquared) {
  const auto X = air();
  const double P = 1.0e5, T = 600.0, A = 0.01;
  const double rho = combaero::density(T, P, X);
  const double a = combaero::speed_of_sound(T, X);
  double prev = 0.0;
  for (double M : {0.1, 0.05, 0.025}) {
    const double mo = M * a * rho * A, mm = 0.8 * mo;
    const auto ci = chamber_merge_face_state(mm, T, X, mo, P, T, X, 0.0, A, true);
    const auto cc = chamber_merge_face_state(mm, T, X, mo, P, T, X, 0.0, A);
    const double rel = (ci.Pt_face - cc.Pt_face) / (cc.Pt_face - P);
    if (prev != 0.0) {
      EXPECT_NEAR(prev / rel, 4.0, 0.2) << M;
    }
    prev = rel;
  }
  EXPECT_LT(std::abs(prev), 1e-4);
}

TEST(ChamberMergeFaceIncompressible, DerivativesMatchCentralDifferences) {
  const auto X = air();
  const double mm = 0.8, mo = 0.95, J = 12.0, P = 1.6e5, T = 1040.0, Tm = 1150.0;
  const double A = 0.01;
  const auto r = chamber_merge_face_state(mm, Tm, X, mo, P, T, X, J, A, true);
  auto f = [&](const double* v) {
    return chamber_merge_face_state(v[0], v[1], X, v[2], v[3], v[4], X, v[5], A, true);
  };
  const double dPf[6] = {r.dPf_dm_main, r.dPf_dT_main, r.dPf_dm_out, r.dPf_dP, r.dPf_dT, r.dPf_dJ};
  const double dPtf[6] = {r.dPtf_dm_main, r.dPtf_dT_main, r.dPtf_dm_out,
                          r.dPtf_dP,      r.dPtf_dT,      r.dPtf_dJ};
  const double h[6] = {1e-6, 1e-3, 1e-6, 1.0, 1e-3, 1e-4};
  for (int k = 0; k < 6; ++k) {
    double up[6] = {mm, Tm, mo, P, T, J}, dn[6] = {mm, Tm, mo, P, T, J};
    up[k] += h[k];
    dn[k] -= h[k];
    const auto a = f(up);
    const auto b = f(dn);
    const double fdP = (a.P_face - b.P_face) / (2 * h[k]);
    const double fdPt = (a.Pt_face - b.Pt_face) / (2 * h[k]);
    EXPECT_NEAR(dPf[k], fdP, 1e-6 * std::abs(fdP) + 1e-8) << k;
    EXPECT_NEAR(dPtf[k], fdPt, 1e-6 * std::abs(fdPt) + 1e-8) << k;
  }
}

TEST(ChamberMergeFace, AnImpulseThatCannotBeCarriedSubsonicallyIsFlagged) {
  const auto X = air();
  // Side jets pushing harder than the main stream's impulse can absorb: the
  // quadratic P_f^2 - Pi P_f + c has no real root (Pi = 0.84e5 against
  // 2 sqrt(c) = 1.17e5 here), so the face would have to choke.
  const auto r = chamber_merge_face_state(1.0, 1200.0, X, 1.0, 1.0e5, 1200.0, X, 500.0, 0.01);
  EXPECT_TRUE(r.choked);
  EXPECT_TRUE(std::isfinite(r.P_face));
  EXPECT_TRUE(std::isfinite(r.Pt_face));
}

// ---- jet impulse -------------------------------------------------------------

TEST(JetImpulse, LowPressureRatioIsBernoulli) {
  const auto X = air();
  const double Pt = 1.001e5, Tt = 600.0, P = 1.0e5, m = 0.01;
  const auto r = jet_impulse(m, Pt, Tt, P, X);
  const double rho = combaero::density(Tt, P, X);
  EXPECT_NEAR(r.dJ_dm, std::sqrt(2.0 * (Pt - P) / rho), 1e-3 * r.dJ_dm);
  EXPECT_NEAR(r.J, m * r.dJ_dm, 1e-12 * r.J);
  EXPECT_FALSE(r.choked);
}

TEST(JetImpulse, PressureDerivativeIsEuler) {
  // Isentropic expansion: w dw = -dP / rho_static.
  const auto X = air();
  const double Pt = 1.4e5, Tt = 700.0, P = 1.0e5, m = 0.02;
  const auto r = jet_impulse(m, Pt, Tt, P, X);
  const double w = r.dJ_dm;
  const double M = combaero::mach_from_pressure_ratio(Tt, Pt, P, X);
  const double Ts = combaero::T_from_stagnation(Tt, M, X);
  const double rho_s = combaero::density(Ts, P, X);
  EXPECT_NEAR(r.dJ_dP / m, -1.0 / (rho_s * w), 1e-4 / (rho_s * w));
}

TEST(JetImpulse, IsContinuousAndSmoothAcrossChoking) {
  const auto X = air();
  const double Pt = 3.0e5, Tt = 600.0, m = 0.02;
  // jet_impulse's own sonic point (the isentrope at M = 1); the golden-
  // section critical_pressure_ratio agrees to ~1e-8.
  const double Tstar = combaero::T_from_stagnation(Tt, 1.0, X);
  const double R = combaero::specific_gas_constant(X);
  const double Pstar =
      Pt * std::exp((combaero::s_mass(Tstar, X, Pt) - combaero::s_mass(Tt, X, Pt)) / R);
  EXPECT_NEAR(Pstar / Pt, combaero::critical_pressure_ratio(Tt, Pt, X), 1e-6);
  const auto a = jet_impulse(m, Pt, Tt, Pstar * (1.0 + 1e-6), X);
  const auto b = jet_impulse(m, Pt, Tt, Pstar * (1.0 - 1e-6), X);
  EXPECT_FALSE(a.choked);
  EXPECT_TRUE(b.choked);
  EXPECT_NEAR(a.J, b.J, 1e-4 * a.J);
  // Analytic on both sides (#816): every partial joins to O(1e-6).
  EXPECT_NEAR(a.dJ_dP, b.dJ_dP, 1e-4 * std::abs(a.dJ_dP));
  EXPECT_NEAR(a.dJ_dPt, b.dJ_dPt, 1e-4 * std::abs(a.dJ_dPt));
  EXPECT_NEAR(a.dJ_dTt, b.dJ_dTt, 1e-4 * std::abs(a.dJ_dTt));
  // Below choking the jet still gains thrust as the back pressure falls.
  EXPECT_GT(jet_impulse(m, Pt, Tt, 0.5 * Pstar, X).J, b.J);
}

TEST(JetImpulse, AnalyticPartialsMatchFiniteDifferences) {
  // #816: the partials were central differences inside C++ and read ~1%
  // off once choked (an inner solve's tolerance inside a Pt*1e-7 step).
  // Analytic now; checked here against differences of J itself, unchoked
  // and choked, over the supply temperature.
  const auto X = air();
  const double m = 0.05;
  for (const double Tt : {300.0, 600.0, 1200.0}) {
    for (const double ratio : {1.05, 1.6, 2.5, 4.0}) {  // Pt / P
      const double P = 1.0e5;
      const double Pt = ratio * P;
      const auto r = jet_impulse(m, Pt, Tt, P, X);
      const double hP = 1e-5 * P;
      const double hPt = 1e-5 * Pt;
      const double hT = 1e-5 * Tt;
      const double fd_P = (jet_impulse(m, Pt, Tt, P + hP, X).J - jet_impulse(m, Pt, Tt, P - hP, X).J) / (2 * hP);
      const double fd_Pt =
          (jet_impulse(m, Pt + hPt, Tt, P, X).J - jet_impulse(m, Pt - hPt, Tt, P, X).J) / (2 * hPt);
      const double fd_T =
          (jet_impulse(m, Pt, Tt + hT, P, X).J - jet_impulse(m, Pt, Tt - hT, P, X).J) / (2 * hT);
      EXPECT_NEAR(r.dJ_dP, fd_P, 1e-4 * std::abs(fd_P)) << Tt << " " << ratio;
      EXPECT_NEAR(r.dJ_dPt, fd_Pt, 1e-4 * std::abs(fd_Pt)) << Tt << " " << ratio;
      EXPECT_NEAR(r.dJ_dTt, fd_T, 1e-4 * std::abs(fd_T)) << Tt << " " << ratio;
    }
  }
}

TEST(JetImpulse, NoOutflowNoMomentum) {
  const auto X = air();
  const auto r = jet_impulse(0.01, 1.0e5, 600.0, 1.0e5, X);
  EXPECT_EQ(r.J, 0.0);
  EXPECT_EQ(r.dJ_dm, 0.0);
}

// ---- area change over a segment (#812) --------------------------------------

TEST(SegmentAreaChangeDrop, IsBernoulliEvenInFlowAndZeroAtEqualAreas) {
  const auto X = air();
  const double P = 1.0e5, T = 300.0, m = 0.1;
  const auto r = combaero::segment_area_change_drop(m, P, T, X, 2e-3, 1e-3);
  const double rho = combaero::density(T, P, X);
  EXPECT_NEAR(r.dP, m * m / (2.0 * rho) * (1.0 / 1e-6 - 1.0 / 4e-6), 1e-9 * r.dP);
  EXPECT_DOUBLE_EQ(combaero::segment_area_change_drop(-m, P, T, X, 2e-3, 1e-3).dP, r.dP);
  const auto z = combaero::segment_area_change_drop(m, P, T, X, 2e-3, 2e-3);
  EXPECT_EQ(z.dP, 0.0);
  EXPECT_EQ(z.d_dm_a, 0.0);
  // An expansion recovers static pressure: a negative drop.
  EXPECT_LT(combaero::segment_area_change_drop(m, P, T, X, 1e-3, 2e-3).dP, 0.0);
}

TEST(SegmentAreaChangeDrop, AnalyticPartials) {
  const auto X = air();
  const double P = 1.2e5, T = 450.0, m = 0.07;
  auto f = [&](double mm, double PP, double TT) {
    return combaero::segment_area_change_drop(mm, PP, TT, X, 3e-3, 1.2e-3).dP;
  };
  const auto r = combaero::segment_area_change_drop(m, P, T, X, 3e-3, 1.2e-3);
  EXPECT_NEAR(r.d_dm_a, (f(m + 1e-6, P, T) - f(m - 1e-6, P, T)) / 2e-6, 1e-7 * std::abs(r.d_dm_a));
  EXPECT_NEAR(r.d_dP, (f(m, P + 1.0, T) - f(m, P - 1.0, T)) / 2.0, 1e-6 * std::abs(r.d_dP));
  EXPECT_NEAR(r.d_dT, (f(m, P, T + 1e-3) - f(m, P, T - 1e-3)) / 2e-3, 1e-6 * std::abs(r.d_dT));
  EXPECT_THROW(combaero::segment_area_change_drop(m, P, T, X, 0.0, 1e-3), std::invalid_argument);
}

// ---- one station, half at a time (#471) ------------------------------------

TEST(StationHalfDrop, BleedAtBassettKappaIsBassettK5ForEverySplit) {
  // Constant density: both halves at one state make the whole station.
  const auto X = air();
  const double P = 1.2e5, T = 500.0, A = 4e-3, m_a = 0.4;
  const double rho = combaero::density(T, P, X);
  for (double q : {0.0, 0.2, 0.5, 0.8, 0.95, 1.0}) {
    const double m_b = q * m_a;
    const auto h = combaero::station_half_drop(m_a, m_b, P, T, X, A,
                                               combaero::STATION_KAPPA_BLEED_BASSETT);
    const double drop = 2.0 * h.dP;  // P_a - P_b
    const double ua = m_a / (rho * A), ub = m_b / (rho * A);
    const double K = (drop + 0.5 * rho * (ua * ua - ub * ub)) / (0.5 * rho * ua * ua);
    EXPECT_NEAR(K, combaero::K5(q), 1e-12) << q;
  }
}

TEST(StationHalfDrop, NormalMergeIsTheImpingementTerm) {
  const auto X = air();
  const double P = 1.0e5, T = 300.0, A = 5e-4;
  const auto h = combaero::station_half_drop(0.015, 0.021, P, T, X, A,
                                             combaero::STATION_KAPPA_MERGE_NORMAL);
  const double rho = combaero::density(T, P, X);
  EXPECT_NEAR(h.dP, 0.5 * (0.021 * 0.021 - 0.015 * 0.015) / (rho * A * A), 1e-12 * h.dP);
}

TEST(StationHalfDrop, IsMirrorSymmetric) {
  // #815: a station traversed the other way -- a and b swapped, both flows
  // negated -- is the same physics, so the drop from a to b reverses sign.
  // The momentum flux m u = m^2/(rho A) is positive either way; the m|m|
  // form this replaced gave the SAME sign, a merge in a reversed chain
  // raising the static pressure along the flow.
  const auto X = air();
  const double P = 1.0e5, T = 300.0, A = 2e-3;
  for (const auto& [ma, mb] : {std::pair{0.10, 0.13}, std::pair{0.10, 0.07}, std::pair{-0.04, 0.05}}) {
    const double fwd = combaero::station_half_drop(ma, mb, P, T, X, A, 0.0).dP;
    const double mir = combaero::station_half_drop(-mb, -ma, P, T, X, A, 0.0).dP;
    EXPECT_NEAR(mir, -fwd, 1e-12 * std::abs(fwd));
  }
  // A symmetric split (the side stream leaves both ways): no net drop.
  EXPECT_NEAR(combaero::station_half_drop(-0.05, 0.05, P, T, X, A, 0.0).dP, 0.0, 1e-12);
}

TEST(StationHalfDrop, DerivativesMatchCentralDifferences) {
  const auto X = air();
  const double P = 1.1e5, T = 420.0, A = 2e-3, e = 1e-8;
  for (double kappa : {0.0, 0.75, 1.0}) {
    for (double ma : {-0.05, 0.08}) {
      for (double mb : {-0.03, 0.06}) {
        const auto r = combaero::station_half_drop(ma, mb, P, T, X, A, kappa);
        auto f = [&](double a, double b, double p, double t) {
          return combaero::station_half_drop(a, b, p, t, X, A, kappa).dP;
        };
        const double fa = (f(ma + e, mb, P, T) - f(ma - e, mb, P, T)) / (2 * e);
        const double fb = (f(ma, mb + e, P, T) - f(ma, mb - e, P, T)) / (2 * e);
        const double fP = (f(ma, mb, P + 1.0, T) - f(ma, mb, P - 1.0, T)) / 2.0;
        const double fT = (f(ma, mb, P, T + 1e-3) - f(ma, mb, P, T - 1e-3)) / 2e-3;
        // Absolute floor: round-off of a ~1e2 Pa value over a 1e-8 step,
        // where the analytic slope can be exactly 0 (2 m_b = kappa m_a).
        EXPECT_NEAR(r.d_dm_a, fa, 1e-6 * std::abs(fa) + 1e-4);
        EXPECT_NEAR(r.d_dm_b, fb, 1e-6 * std::abs(fb) + 1e-4);
        EXPECT_NEAR(r.d_dP, fP, 1e-5 * std::abs(fP) + 1e-12);
        EXPECT_NEAR(r.d_dT, fT, 1e-5 * std::abs(fT) + 1e-12);
      }
    }
  }
}

TEST(ChannelEntryDrop, IsTheDynamicHeadPlusTheEntryLoss) {
  const auto X = air();
  const double P = 1.0e5, T = 300.0, A = 1e-3, m = 0.05;
  const double rho = combaero::density(T, P, X);
  const double u = m / (rho * A);
  const auto r = combaero::channel_entry_drop(m, P, T, X, A, 0.5);
  EXPECT_NEAR(r.dP, 1.5 * 0.5 * rho * u * u, 1e-9 * r.dP);
  const double e = 1e-8;
  const double fd = (combaero::channel_entry_drop(m + e, P, T, X, A, 0.5).dP -
                     combaero::channel_entry_drop(m - e, P, T, X, A, 0.5).dP) / (2 * e);
  EXPECT_NEAR(r.d_dm_a, fd, 1e-6 * fd);
}

TEST(ChannelEntryDrop, AReversedFlowDischargesAtStaticContinuity) {
  // Out of the duct into the reservoir: no acceleration to pay, and the
  // dynamic head is lost -- not the mirrored entry (#705).
  const auto X = air();
  const double P = 1.0e5, T = 300.0, A = 1e-3;
  const auto r = combaero::channel_entry_drop(-0.05, P, T, X, A, 0.5);
  EXPECT_EQ(r.dP, 0.0);
  EXPECT_EQ(r.d_dm_a, 0.0);
  EXPECT_EQ(r.d_dP, 0.0);
  EXPECT_EQ(r.d_dT, 0.0);
  // C1 through zero: value and slope both vanish from either side.
  const double e = 1e-9;
  const double slope = combaero::channel_entry_drop(0.05, P, T, X, A, 0.5).d_dm_a;
  EXPECT_NEAR(combaero::channel_entry_drop(e, P, T, X, A, 0.5).d_dm_a, 0.0, 1e-6 * slope);
  EXPECT_EQ(combaero::channel_entry_drop(0.0, P, T, X, A, 0.5).dP, 0.0);
}

// ---- supply-side crossflow ratio (#471) ------------------------------------

TEST(CrossflowRatio, IsTheMeanDuctVelocityOverTheStaticIdealJet) {
  const auto X = air();
  const double ma = 0.30, mb = 0.28, P = 1.05e5, T = 600.0, Pd = 1.0e5, A = 1.5e-3;
  const auto r = combaero::crossflow_velocity_ratio(ma, mb, P, T, X, Pd, A);
  const double rho = combaero::density(T, P, X);
  EXPECT_NEAR(r.U1, 0.5 * (ma + mb) / (rho * A), 1e-12 * r.U1);
  // Low pressure ratio: Bernoulli on the STATIC drive.
  EXPECT_NEAR(r.Vi, std::sqrt(2.0 * (P - Pd) / rho), 2e-2 * r.Vi);
  EXPECT_NEAR(r.U1_over_Vi, r.U1 / r.Vi, 1e-12);
}

TEST(CrossflowRatio, DerivativesMatchCentralDifferences) {
  const auto X = air();
  const double ma = 0.30, mb = 0.28, P = 1.05e5, T = 600.0, Pd = 1.0e5, A = 1.5e-3;
  const auto r = combaero::crossflow_velocity_ratio(ma, mb, P, T, X, Pd, A);
  auto f = [&](double a, double b, double p, double t, double pd) {
    return combaero::crossflow_velocity_ratio(a, b, p, t, X, pd, A).U1_over_Vi;
  };
  const double e = 1e-7, hP = 1.0, hT = 1e-3;
  EXPECT_NEAR(r.d_dm_a, (f(ma + e, mb, P, T, Pd) - f(ma - e, mb, P, T, Pd)) / (2 * e), 1e-5 * std::abs(r.d_dm_a));
  EXPECT_NEAR(r.d_dm_b, (f(ma, mb + e, P, T, Pd) - f(ma, mb - e, P, T, Pd)) / (2 * e), 1e-5 * std::abs(r.d_dm_b));
  const double fP = (f(ma, mb, P + hP, T, Pd) - f(ma, mb, P - hP, T, Pd)) / (2 * hP);
  const double fT = (f(ma, mb, P, T + hT, Pd) - f(ma, mb, P, T - hT, Pd)) / (2 * hT);
  const double fPd = (f(ma, mb, P, T, Pd + hP) - f(ma, mb, P, T, Pd - hP)) / (2 * hP);
  EXPECT_NEAR(r.d_dP, fP, 1e-3 * std::abs(fP));
  EXPECT_NEAR(r.d_dT, fT, 1e-3 * std::abs(fT));
  EXPECT_NEAR(r.d_dP_down, fPd, 1e-3 * std::abs(fPd));
}

TEST(CrossflowRatio, VanishingDriveStaysFinite) {
  const auto X = air();
  const auto r = combaero::crossflow_velocity_ratio(0.3, 0.3, 1.0e5, 300.0, X, 1.0e5, 1e-3);
  EXPECT_TRUE(std::isfinite(r.U1_over_Vi));
  EXPECT_NEAR(r.Vi, combaero::CROSSFLOW_VI_FLOOR, 1e-12);
}

// ---- sonic static state (#637) ------------------------------------------------

TEST(SonicStaticState, TheFluxIsSonicThereAndTheDerivativesAreExact) {
  const auto X = air();
  for (const auto& [T0, G] : std::vector<std::pair<double, double>>{
           {300.0, 250.0}, {1600.0, 80.0}, {900.0, 1.0}}) {
    const auto s = combaero::sonic_static_state(T0, G, X);
    // At P* the chamber's own closure reads Mach 1.
    EXPECT_NEAR(combaero::chamber_static_state(T0, s.P, G, X).M, 1.0, 1e-10);
    const double hT = 1e-4 * T0;
    const double hG = 1e-5 * G;
    const double fdT = (combaero::sonic_static_state(T0 + hT, G, X).P -
                        combaero::sonic_static_state(T0 - hT, G, X).P) / (2.0 * hT);
    const double fdG = (combaero::sonic_static_state(T0, G + hG, X).P -
                        combaero::sonic_static_state(T0, G - hG, X).P) / (2.0 * hG);
    EXPECT_NEAR(s.dP_dTt / fdT, 1.0, 1e-7);
    EXPECT_NEAR(s.dP_dG / fdG, 1.0, 1e-7);
  }
  EXPECT_EQ(combaero::sonic_static_state(300.0, 0.0, X).P, 0.0);
}

// ---- composition partials (#1146) ------------------------------------------
//
// Each kernel's d/dY against central differences in the mass fractions, along
// directions that keep sum Y = 1 (e_k - e_ref): the solver's unknowns. Every
// species is present, so no step crosses Y = 0 (where s() drops a species).

namespace {
std::vector<double> mixed_Y() {
  std::vector<double> X = air();
  X[combaero::species_index_from_name("CH4")] = 0.03;
  X[combaero::species_index_from_name("CO2")] = 0.02;
  X[combaero::species_index_from_name("H2O")] = 0.04;
  const double n = static_cast<double>(X.size());
  for (double& x : X) {
    x = 0.99 * x + 0.01 / n;
  }
  double s = 0.0;
  for (double x : X) {
    s += x;
  }
  for (double& x : X) {
    x /= s;
  }
  return combaero::mole_to_mass(X);
}

std::vector<double> X_of(std::vector<double> Y) {
  double s = 0.0;
  for (double y : Y) {
    s += y;
  }
  for (double& y : Y) {
    y /= s;
  }
  return combaero::mass_to_mole(Y);
}

// max over k of |analytic - FD| along e_k - e_ref, over the largest |FD|.
template <class F>
double dY_error(F value, const std::vector<double>& dY, const std::vector<double>& Y) {
  const std::size_t ref = combaero::species_index_from_name("N2");
  const double h = 1e-7;
  double err = 0.0, scale = 0.0;
  for (std::size_t k = 0; k < Y.size(); ++k) {
    if (k == ref) {
      continue;
    }
    std::vector<double> yp = Y, ym = Y;
    yp[k] += h;
    yp[ref] -= h;
    ym[k] -= h;
    ym[ref] += h;
    const double fd = (value(X_of(yp)) - value(X_of(ym))) / (2.0 * h);
    err = std::max(err, std::abs(dY[k] - dY[ref] - fd));
    scale = std::max(scale, std::abs(fd));
  }
  return err / scale;
}
}  // namespace

TEST(CompositionPartials, StationDropsMatchCentralDifferences) {
  const auto Y = mixed_Y();
  const auto X = X_of(Y);
  const double P = 1.0e5, T = 600.0, A = 2e-4;
  for (double m : {0.012, -0.012}) {
    const auto s = combaero::station_half_drop(0.01, m, P, T, X, A, 0.75);
    ASSERT_EQ(s.d_dY.size(), X.size());
    EXPECT_LT(dY_error([&](const std::vector<double>& x) {
      return combaero::station_half_drop(0.01, m, P, T, x, A, 0.75).dP;
    }, s.d_dY, Y), 1e-7);
    const auto a = combaero::segment_area_change_drop(m, P, T, X, A, 1.5 * A);
    EXPECT_LT(dY_error([&](const std::vector<double>& x) {
      return combaero::segment_area_change_drop(m, P, T, x, A, 1.5 * A).dP;
    }, a.d_dY, Y), 1e-7);
  }
  const auto e = combaero::channel_entry_drop(0.012, P, T, X, A, 0.5);
  EXPECT_LT(dY_error([&](const std::vector<double>& x) {
    return combaero::channel_entry_drop(0.012, P, T, x, A, 0.5).dP;
  }, e.d_dY, Y), 1e-7);
  // Zero drops carry zero partials, one per species.
  for (const auto& z : {combaero::channel_entry_drop(-0.012, P, T, X, A, 0.5),
                        combaero::segment_area_change_drop(0.012, P, T, X, A, A)}) {
    ASSERT_EQ(z.d_dY.size(), X.size());
    for (double d : z.d_dY) {
      EXPECT_EQ(d, 0.0);
    }
  }
}

TEST(CompositionPartials, JetImpulseMatchesCentralDifferencesOnBothBranches) {
  const auto Y = mixed_Y();
  const auto X = X_of(Y);
  const double m = 0.02, Tt = 700.0, P = 1.0e5;
  // Unchoked, just past the sonic point, and deep in the choked branch.
  for (double Pt : {1.2e5, 1.9e5, 3.0e5}) {
    const auto r = jet_impulse(m, Pt, Tt, P, X);
    ASSERT_EQ(r.dJ_dY.size(), X.size());
    EXPECT_LT(dY_error([&](const std::vector<double>& x) {
      return jet_impulse(m, Pt, Tt, P, x).J;
    }, r.dJ_dY, Y), 1e-6) << "Pt " << Pt << " choked " << r.choked;
  }
  // No outflow: zero partials, one per species.
  const auto z = jet_impulse(m, 0.9e5, Tt, P, X);
  ASSERT_EQ(z.dJ_dY.size(), X.size());
  for (double d : z.dJ_dY) {
    EXPECT_EQ(d, 0.0);
  }
}
