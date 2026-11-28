// 1-D supersonic ejector model (critical / subcritical / unchoked jet-pump
// regimes). See include/ejector.h for paper references and the scope note.

#include "ejector.h"

#include "dual_number.h"
#include "flow_drive_floor.h"
#include "input_validation.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <utility>

namespace combaero::solver {

// -----------------------------------------------------------------------------
// Input validation (#769). Every public entry point refuses non-finite input
// and parameters outside their physical domain with std::invalid_argument
// (Python ValueError, which the network solver treats as a failed
// evaluation). Two classes of argument:
//   - parameters (geometry, gamma, r_gas, efficiencies, recovery, smoothing
//     widths, loss coefficients): finite and in their physical domain;
//   - state: in the standalone closures the supply stagnation state must be
//     finite and positive; in the two whole-element residual functions, which
//     take Newton iterates, only finite. A trial point legitimately carries a
//     negative pressure or a reversed flow there (measured on the 'ejector'
//     stress family), and the closures continue through C1 floors (#481).
// The element functions call the unchecked *_raw closures below, so a probe
// outside the closures' physical domain is not refused by them.
// -----------------------------------------------------------------------------
namespace {

[[noreturn]] void refuse(const char* who, const char* what) {
  throw std::invalid_argument(std::string(who) + ": " + what);
}

void require_finite_arg(const char* who, double x, const char* what) {
  if (!validation::is_finite(x)) refuse(who, what);
}

void require_positive_arg(const char* who, double x, const char* what) {
  if (!validation::is_finite_positive(x)) refuse(who, what);
}

void check_gamma(const char* who, double gamma) {
  if (!validation::is_finite(gamma) || !(gamma > 1.0)) {
    refuse(who, "gamma must be finite and > 1");
  }
}

void check_gas(const char* who, double gamma, double r_gas) {
  check_gamma(who, gamma);
  require_positive_arg(who, r_gas, "r_gas must be finite and > 0");
}

// A nozzle efficiency on kinetic energy: 0 < eta <= 1.
void check_efficiency(const char* who, double eta, const char* what) {
  if (!validation::is_finite_in_range(eta, 0.0, 1.0) || !(eta > 0.0)) refuse(who, what);
}

void check_geometry(const char* who, const EjectorGeometry& geom) {
  if (!validation::is_finite(geom.area_ratio_nozzle) || !(geom.area_ratio_nozzle > 1.0)) {
    refuse(who, "area_ratio_nozzle (A_p1/A_t) must be finite and > 1");
  }
  if (!validation::is_finite(geom.area_ratio_mix) ||
      !(geom.area_ratio_mix >= geom.area_ratio_nozzle)) {
    // A_3 = A_p1 (no suction annulus) is the closures' normal-shock limit.
    refuse(who, "area_ratio_mix (A_3/A_t) must be finite and >= area_ratio_nozzle");
  }
}

// Supply stagnation states of the standalone closures.
void check_supplies(const char* who, double p_g, double t_g, double p_e, double t_e) {
  if (!validation::is_finite_positive(p_g) || !validation::is_finite_positive(t_g) ||
      !validation::is_finite_positive(p_e) || !validation::is_finite_positive(t_e)) {
    refuse(who, "p_g, t_g, p_e and t_e must be finite and > 0");
  }
}

void check_nozzle(const char* who, double Pt, double Tt, double area_throat, double area_exit,
                  double gamma, double r_gas, double eta, double eps_frac) {
  require_positive_arg(who, Pt, "Pt must be finite and > 0");
  require_positive_arg(who, Tt, "Tt must be finite and > 0");
  require_positive_arg(who, area_throat, "area_throat must be finite and > 0");
  if (!validation::is_finite(area_exit) || !(area_exit >= area_throat)) {
    refuse(who, "area_exit must be finite and >= area_throat");
  }
  check_gas(who, gamma, r_gas);
  check_efficiency(who, eta, "eta must be finite and in (0, 1]");
  if (!validation::is_finite_nonnegative(eps_frac)) {
    refuse(who, "eps_frac must be finite and >= 0");
  }
}

// A recovery on the lossless mixed stagnation pressure: 0 < recovery <= 1.
// Above 1 the ejector would gain pressure beyond lossless mixing.
void check_recovery(const char* who, double recovery, const char* name) {
  if (!validation::is_finite_in_range(recovery, 0.0, 1.0) || !(recovery > 0.0)) {
    throw std::invalid_argument(std::string(who) + ": " + name +
                                " must be finite and in (0, 1]; above 1 is a pressure gain "
                                "beyond lossless mixing");
  }
}

void check_window(const char* who, double lo, double hi, const char* what) {
  if (!validation::is_finite(lo) || !validation::is_finite(hi) || !(lo < hi)) refuse(who, what);
}

} // namespace

// -----------------------------------------------------------------------------
// The closures, written once against DualN<N> (include/dual_number.h).
//
// Every closure below takes its inputs as duals, gamma and r_gas included, so
// one implementation serves both callers: the public *_and_jacobian functions
// seed (p_g, t_g, p_e, t_e) (or their own inputs) and hold gamma and r_gas
// constant, and the whole element seeds gamma and r_gas too (#1246). gamma and
// R are the secondary port's, and that port's gas moves with the flows when a
// leg reverses or a side stream feeds the port chamber, so the element's rows
// need their partials in both.
//
// Value parity. The closures were written against two dual types: DualN for
// the operating-regime closures and a fixed-seed Dual4 for the critical path,
// whose quotient was a * (1 / b). The critical-path closures keep that
// quotient through ddiv4 below, every expression keeps its operand order, and
// a gamma-only expression the earlier code wrote in doubles, where the
// compiler may fuse a * b + c, keeps that double form for its value
// (revalue). DualN's quotient of two duals is the exact a / b the double one
// was. With the C-D nozzle's exit Mach in its cancellation-free form
// (cd_nozzle_forward) the only exception, the element's residuals and every
// closure value are the earlier code's bit for bit (measured at 3000 random
// element states, #1246).
// -----------------------------------------------------------------------------
namespace {

// The critical path's quotient: a * (1 / b) (see "Value parity" above).
template <int N> DualN<N> ddiv4(const DualN<N>& a, const DualN<N>& b) { return a * (1.0 / b); }
template <int N> DualN<N> ddiv4(double a, const DualN<N>& b) {
  return DualN<N>::constant(a) * (1.0 / b);
}

// a with its value replaced by v: the value of a gamma-only expression the
// earlier code evaluated in doubles, where the compiler may fuse a * b + c
// (and so round it differently from the dual operators), with the dual's
// partials (see "Value parity" above).
template <int N> DualN<N> revalue(DualN<N> a, double v) {
  a.v = v;
  return a;
}

// a^b with both a dual: d(a^b) = b a^(b-1) da + a^b ln(a) db. The value is
// std::pow's, as for a double exponent.
template <int N> DualN<N> dpow_dd(const DualN<N>& a, const DualN<N>& b) {
  DualN<N> r;
  r.v = std::pow(a.v, b.v);
  const double da = b.v * std::pow(a.v, b.v - 1.0);
  const double db = r.v * std::log(a.v);
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * da + b.d[i] * db;
  return r;
}

// dpos_floor with a dual floor eps: a for a >= eps, eps^2 / (2 eps - a)
// below, whose eps-partial is 2 eps (eps - a) / (2 eps - a)^2.
template <int N> DualN<N> dpos_floor_dd(const DualN<N>& a, const DualN<N>& eps) {
  if (a.v >= eps.v) {
    return a;
  }
  DualN<N> r;
  const double den = 2.0 * eps.v - a.v;
  r.v = eps.v * eps.v / den;
  const double da = r.v / den;
  const double de = 2.0 * eps.v * (eps.v - a.v) / (den * den);
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * da + eps.d[i] * de;
  return r;
}

template <int N> DualN<N> area_over_astar_d(const DualN<N>& mach, const DualN<N>& gamma) {
  const DualN<N> gm1 = gamma - 1.0;
  const DualN<N> gp1 = gamma + 1.0;
  return (1.0 / mach) *
         dpow_dd((2.0 / gp1) * (1.0 + 0.5 * gm1 * mach * mach), gp1 / (2.0 * gm1));
}

} // namespace

// The bisection's own A/A*, in doubles: area_over_astar_d is the same
// expression, for its partials (nozzle_exit_mach_d).
static double ejector_area_over_astar_raw(double mach, double gamma) {
  double gm1 = gamma - 1.0;
  double gp1 = gamma + 1.0;
  return (1.0 / mach) *
         std::pow((2.0 / gp1) * (1.0 + 0.5 * gm1 * mach * mach), gp1 / (2.0 * gm1));
}

// Supersonic root of A/A* = area_ratio by bisection on the Mach bracket
// [kSupersonicMachLo, kSupersonicMachHi]. An area ratio outside the bracket's
// image has no root there: the bisection then ends on the bracket end
// (clamp), and out_of_range says which side (#945).
static EjectorNozzleExitMach ejector_nozzle_exit_mach_raw(double area_ratio, double gamma) {
  double lo = ejector::kSupersonicMachLo;
  double hi = ejector::kSupersonicMachHi;
  EjectorNozzleExitMach out{};
  if (area_ratio < ejector_area_over_astar_raw(lo, gamma)) {
    out.out_of_range = -1;
  } else if (area_ratio > ejector_area_over_astar_raw(hi, gamma)) {
    out.out_of_range = 1;
  }
  for (int i = 0; i < ejector::kBisectionIterations; ++i) {
    double mid = 0.5 * (lo + hi);
    if (ejector_area_over_astar_raw(mid, gamma) < area_ratio) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  out.mach = 0.5 * (lo + hi);
  return out;
}

static double ejector_mach_from_area_ratio_supersonic_raw(double area_ratio, double gamma) {
  return ejector_nozzle_exit_mach_raw(area_ratio, gamma).mach;
}

namespace {

// M_p1 as a dual: the bisection's root, constant in every input but gamma.
// In range, A/A*(M_p1, gamma) = area_ratio gives dM_p1/dgamma =
// -(dA/dgamma) / (dA/dM) (implicit function theorem; dA/dM > 0 for M > 1);
// the bisection resolves the root to rounding (200 halvings), so this is the
// derivative of the value used. A clamped M_p1 is the bracket end, a constant.
template <int N>
DualN<N> nozzle_exit_mach_d(double area_ratio, const DualN<N>& gamma) {
  const EjectorNozzleExitMach m = ejector_nozzle_exit_mach_raw(area_ratio, gamma.v);
  DualN<N> r = DualN<N>::constant(m.mach);
  if (m.out_of_range != 0) {
    return r;
  }
  using D2 = DualN<2>;
  const D2 a = area_over_astar_d(D2::seed(m.mach, 0), D2::seed(gamma.v, 1));
  const double k = -a.d[1] / a.d[0];
  for (int i = 0; i < N; ++i) r.d[i] = k * gamma.d[i];
  return r;
}

template <int N>
DualN<N> choked_mass_flow_d(const DualN<N>& Pt, const DualN<N>& Tt, double area_throat,
                            const DualN<N>& gamma, const DualN<N>& r_gas, double eta) {
  const DualN<N> gm1 = gamma - 1.0;
  const DualN<N> gp1 = gamma + 1.0;
  const DualN<N> choke_const = dsqrt(gamma / r_gas * dpow_dd(2.0 / gp1, gp1 / gm1));
  return (Pt * area_throat / dsqrt(Tt)) * choke_const * std::sqrt(eta);
}

} // namespace

static double ejector_choked_mass_flow_raw(double Pt, double Tt, double area_throat,
                                double gamma, double r_gas, double eta) {
  using D1 = DualN<1>;
  return choked_mass_flow_d(D1::constant(Pt), D1::constant(Tt), area_throat,
                            D1::constant(gamma), D1::constant(r_gas), eta)
      .v;
}

static std::tuple<double, double, double>
ejector_choked_mass_flow_and_jacobian_raw(double Pt, double Tt, double area_throat,
                                      double gamma, double r_gas, double eta) {
  double mdot = ejector_choked_mass_flow_raw(Pt, Tt, area_throat, gamma, r_gas, eta);
  // mdot = C * Pt / sqrt(Tt)  =>  d/dPt = mdot/Pt,  d/dTt = -mdot/(2*Tt)
  return {mdot, mdot / Pt, -0.5 * mdot / Tt};
}

namespace {

// Raw forward C-D nozzle flow from stagnation p_up to exit static p_down < p_up:
// the subsonic exit flux smooth-min'd with the sonic-throat cap.
double cd_nozzle_forward_raw(double p_up, double Tt, double p_down, double area_throat,
                             double area_exit, double gamma, double r_gas, double eta,
                             double eps_frac) {
  double cap = ejector_choked_mass_flow_raw(p_up, Tt, area_throat, gamma, r_gas, eta);
  double exit_flux = 0.0;
  double r = p_down / p_up;
  if (r < 1.0) {
    double gm1 = gamma - 1.0;
    double gp1 = gamma + 1.0;
    // r^(-(g-1)/g) - 1 as expm1(-(g-1)/g log1p(x)), x = (p_down - p_up) / p_up:
    // the pow form cancels as r -> 1, a relative error ~ eps / (1 - r) in M^2
    // (1.5e-7 at 1 - r = 1e-9), and so in the flow and its P_py slope.
    const double x = (p_down - p_up) / p_up;
    double mach = std::sqrt(2.0 / gm1 * std::expm1(-gm1 / gamma * std::log1p(x)));
    // Below the critical ratio the throat is choked: hold the exit flux at
    // its sonic maximum (flux'(M = 1) = 0, so C1). The supersonic-exit flux
    // FALLS with M and dipped under the cap at deep expansion (r < 0.16 at
    // A_e/A_t = 1.5), a non-monotone nozzle flow in P_py (#619).
    mach = std::min(mach, 1.0);
    double flux = mach * std::pow(1.0 + 0.5 * gm1 * mach * mach, -gp1 / (2.0 * gm1));
    exit_flux = area_exit * p_up / std::sqrt(Tt) * std::sqrt(gamma / r_gas) * flux * std::sqrt(eta);
  }
  double eps = eps_frac * cap;
  double diff = exit_flux - cap;
  return 0.5 * (exit_flux + cap - std::sqrt(diff * diff + eps * eps));
}

// Forward C-D nozzle flow with the shared zero-drive sub-floor law
// (flow_drive_floor.h; #1259, #1269): below kJetDriveFloor of drive
// p_up - p_down (when p_up > kJetDriveFloor), the flow is the floor's flow
// times the odd cubic h(t) with exponent kJetFloorExponent (1/2). Above the
// floor it is bit-identical to cd_nozzle_forward_raw; at zero drive (r = 1)
// the flow vanishes with finite slope 1.25 * m(floor) / floor.
double cd_nozzle_forward(double p_up, double Tt, double p_down, double area_throat,
                         double area_exit, double gamma, double r_gas, double eta,
                         double eps_frac) {
  const double drive = p_up - p_down;
  if (drive >= kJetDriveFloor || !(p_up > kJetDriveFloor)) {
    return cd_nozzle_forward_raw(p_up, Tt, p_down, area_throat, area_exit, gamma, r_gas, eta,
                                 eps_frac);
  }
  const double m_floor = cd_nozzle_forward_raw(p_up, Tt, p_up - kJetDriveFloor, area_throat,
                                               area_exit, gamma, r_gas, eta, eps_frac);
  const DriveFloorCubic c = drive_floor_cubic(drive / kJetDriveFloor, kJetFloorExponent);
  return m_floor * c.h;
}

} // namespace

static double ejector_cd_nozzle_mass_flow_raw(double Pt, double Tt, double p_static_down,
                                   double area_throat, double area_exit,
                                   double gamma, double r_gas, double eta,
                                   double eps_frac) {
  if (p_static_down <= Pt) {
    return cd_nozzle_forward(Pt, Tt, p_static_down, area_throat, area_exit, gamma, r_gas, eta,
                             eps_frac);
  }
  // Exit static above the supply (#619): the nozzle passes REVERSE flow, the
  // same closure from p_static_down back into Pt. With the sub-floor cubic
  // (#1269), forward(Pt, Pt) = 0 with a finite slope matching the reverse side
  // at r = 1 (when Pt > kJetDriveFloor); for an unphysical reservoir below the
  // floor, retain the smooth-min r = 1 offset for continuity.
  if (Pt > kJetDriveFloor) {
    return -cd_nozzle_forward(p_static_down, Tt, Pt, area_throat, area_exit, gamma, r_gas, eta,
                              eps_frac);
  }
  double cap = ejector_choked_mass_flow_raw(Pt, Tt, area_throat, gamma, r_gas, eta);
  double eps = eps_frac * cap;
  double offset = std::sqrt(cap * cap + eps * eps) - cap;
  return -cd_nozzle_forward(p_static_down, Tt, Pt, area_throat, area_exit, gamma, r_gas, eta,
                            eps_frac) -
         offset;
}

double ejector_q_lambda(double lam, double gamma) {
  if (!validation::is_finite_nonnegative(lam)) {
    refuse("ejector_q_lambda", "lambda must be finite and >= 0");
  }
  check_gamma("ejector_q_lambda", gamma);
  double gm1 = gamma - 1.0;
  double gp1 = gamma + 1.0;
  return std::pow(1.0 - (gm1 / gp1) * lam * lam, 1.0 / gm1) * std::pow(gp1 / 2.0, 1.0 / gm1) * lam;
}

namespace {

// The C-D nozzle (R0's primary, and the jet pump's secondary flux) in duals.
template <int N>
DualN<N> cd_nozzle_d(const DualN<N>& p0, const DualN<N>& t0, const DualN<N>& p_py,
                     double area_throat, double area_exit, const DualN<N>& gamma,
                     const DualN<N>& r_gas, double eta, double eps_frac) {
  using D = DualN<N>;
  const D gm1 = gamma - 1.0;
  const D gp1 = gamma + 1.0;
  const D choke_const = dsqrt(gamma / r_gas * dpow_dd(2.0 / gp1, gp1 / gm1));
  auto cap_of = [&](const D& p_up) {
    return (choke_const * area_throat * std::sqrt(eta)) * p_up / dsqrt(t0);
  };
  // eps = eps_frac * cap depends on (p_up, t0); differentiate it too (FD of
  // the Python value function does), else the Jacobian would be inconsistent.
  auto forward_raw = [&](const D& p_up, const D& p_down) {
    D cap = cap_of(p_up);
    D exit_flux = D::constant(0.0);
    if (p_down.v < p_up.v) {
      // Cancellation-free as r = p_down / p_up -> 1 (see cd_nozzle_forward_raw).
      D x = (p_down - p_up) / p_up;
      D mach_sq = (2.0 / gm1) * dexpm1((-gm1 / gamma) * dlog1p(x));
      D mach = dsqrt(mach_sq);
      if (mach.v > 1.0) mach = D::constant(1.0); // choked: see ejector_cd_nozzle_mass_flow
      D flux = mach * dpow_dd(1.0 + 0.5 * gm1 * mach * mach, -gp1 / (2.0 * gm1));
      exit_flux = (area_exit * dsqrt(gamma / r_gas) * std::sqrt(eta)) * p_up / dsqrt(t0) * flux;
    }
    D eps = eps_frac * cap;
    D diff = exit_flux - cap;
    return 0.5 * (exit_flux + cap - dsqrt(diff * diff + eps * eps));
  };
  auto forward = [&](const D& p_up, const D& p_down) {
    const double drive = p_up.v - p_down.v;
    if (drive >= kJetDriveFloor || !(p_up.v > kJetDriveFloor)) {
      return forward_raw(p_up, p_down);
    }
    // Sub-floor odd cubic (#1259, #1269): p_floor = p_up - kJetDriveFloor
    // carries p_up's dual seeds automatically so every partial (p0, t0, p_py,
    // gamma, r_gas) is exact across the band and at zero drive (r = 1).
    const D p_floor = p_up - kJetDriveFloor;
    const D m_floor = forward_raw(p_up, p_floor);
    const D t = (p_up - p_down) / kJetDriveFloor;
    const DriveFloorCubic c = drive_floor_cubic(t.v, kJetFloorExponent);
    D h = D::constant(c.h);
    for (int i = 0; i < N; ++i) {
      h.d[i] = c.dh * t.d[i];
    }
    return m_floor * h;
  };

  if (p_py.v <= p0.v) {
    return forward(p0, p_py);
  }
  // Reverse flow above the supply (#619, #1269): see ejector_cd_nozzle_mass_flow_raw.
  if (p0.v > kJetDriveFloor) {
    return -forward(p_py, p0);
  }
  D cap = cap_of(p0);
  D eps = eps_frac * cap;
  return -forward(p_py, p0) - (dsqrt(cap * cap + eps * eps) - cap);
}

} // namespace

static EjectorCDNozzleJacobian ejector_cd_nozzle_mass_flow_and_jacobian_raw(
    double Pt_val, double Tt_val, double p_py_val, double area_throat,
    double area_exit, double gamma, double r_gas, double eta, double eps_frac) {
  // Seeds: 0 = p0 (Pt), 1 = t0 (Tt), 2 = p_py.
  using D = DualN<3>;
  const D mdot = cd_nozzle_d(D::seed(Pt_val, 0), D::seed(Tt_val, 1), D::seed(p_py_val, 2),
                             area_throat, area_exit, D::constant(gamma), D::constant(r_gas), eta,
                             eps_frac);
  return {mdot.v, mdot.d[0], mdot.d[1], mdot.d[2]};
}

EjectorCDExitStaticJacobian ejector_cd_nozzle_exit_static_and_jacobian(
    double m_dot, double Pt, double Tt, double area_throat, double area_exit,
    double gamma, double r_gas, double eta, double eps_frac) {
  constexpr const char* who = "ejector_cd_nozzle_exit_static_and_jacobian";
  if (!validation::is_finite_nonnegative(m_dot)) refuse(who, "m_dot must be finite and >= 0");
  check_nozzle(who, Pt, Tt, area_throat, area_exit, gamma, r_gas, eta, eps_frac);
  double gm1 = gamma - 1.0;
  double gp1 = gamma + 1.0;
  double r_crit = std::pow(2.0 / gp1, gamma / gm1);
  // cd_nozzle_mass_flow is monotone decreasing in P_py: bisect for the root.
  const double hi_frac = 1.0 - ejector::kExitStaticTopGap;
  double lo = Pt * r_crit;
  double hi = Pt * hi_frac;
  // An m_dot outside [flow(hi), flow(lo)] has no root in the bracket: the
  // bisection ends on the bracket end it is pushed to (flag and clamp, #945).
  const double m_lo_end = ejector_cd_nozzle_mass_flow_raw(Pt, Tt, lo, area_throat, area_exit,
                                                          gamma, r_gas, eta, eps_frac);
  const double m_hi_end = ejector_cd_nozzle_mass_flow_raw(Pt, Tt, hi, area_throat, area_exit,
                                                          gamma, r_gas, eta, eps_frac);
  int out_of_range = 0;
  if (m_dot > m_lo_end) {
    out_of_range = 1;  // above the choked end: P_py held at P*
  } else if (m_dot < m_hi_end) {
    out_of_range = -1;  // below the near-zero end: P_py held at Pt (1 - gap)
  }
  for (int i = 0; i < ejector::kBisectionIterations; ++i) {
    double mid = 0.5 * (lo + hi);
    double f = ejector_cd_nozzle_mass_flow_raw(Pt, Tt, mid, area_throat, area_exit, gamma, r_gas,
                                           eta, eps_frac);
    if (f > m_dot) {
      lo = mid; // too much flow -> raise P_py
    } else {
      hi = mid;
    }
  }
  double p_py = 0.5 * (lo + hi);
  if (out_of_range != 0) {
    // The clamped function is the bracket end, k Pt (k = r_crit or 1 - gap):
    // flat in m_dot and Tt, so its Jacobian is (0, k, 0), consistent with
    // the value. The implicit-function partials below would describe a root
    // that does not exist.
    const double k = out_of_range > 0 ? r_crit : hi_frac;
    return {p_py, 0.0, k, 0.0, out_of_range};
  }
  // Implicit function theorem: F(P_py) = cd_nozzle(P_py) - m_dot = 0, so
  // dP_py/dx = -(dF/dx)/(dF/dP_py), using the C-D nozzle's own Jacobian.
  EjectorCDNozzleJacobian j = ejector_cd_nozzle_mass_flow_and_jacobian_raw(
      Pt, Tt, p_py, area_throat, area_exit, gamma, r_gas, eta, eps_frac);
  double inv = 1.0 / j.dmdot_dp_py;
  return {p_py, inv, -j.dmdot_dp0 * inv, -j.dmdot_dt0 * inv, 0};
}

namespace {

// The jet pump's mixed-flow stagnation pressure (lossless, recovery applied
// by the caller). gamma is a dual; r_gas cancels (lambda and q are
// dimensionless).
template <int N>
DualN<N> jetpump_discharge_d(DualN<N> p_g, const DualN<N>& t_g, DualN<N> p_e,
                             const DualN<N>& t_e, const DualN<N>& p_py, const DualN<N>& omega,
                             const DualN<N>& gamma, double eta_primary, double eta_secondary) {
  using D = DualN<N>;
  const D gm1 = gamma - 1.0;
  const D gp1 = gamma + 1.0;
  const D k = gm1 / gamma;

  // Both streams expand isentropically to the common mixing static P_py:
  // lambda = sqrt((g+1)/(g-1) * (1 - (P_py/P0)^((g-1)/g))). Subsonic (< 1).
  //
  // Out of domain (#481): a solver probe can put P_py above a supply's Pt
  // (no expansion: lambda^2 < 0), at or below zero, or reverse a flow. The
  // arguments then continue through a C1 positive floor (dpos_floor), so the
  // discharge stays finite and smooth instead of NaN -- which otherwise
  // poisoned the whole blended residual over most of the low-primary range.
  // Exact wherever the closure is defined.
  constexpr double kRatioFloor = 1e-12;
  constexpr double kLambdaSqFloor = 1e-8;
  auto lambda_of = [&](const D& p0) {
    const D r = dpos_floor(p_py / p0, kRatioFloor);
    return dsqrt(dpos_floor((gp1 / gm1) * (1.0 - dpow_dd(r, k)), kLambdaSqFloor));
  };
  // The velocities the nozzles deliver (#588): lambda = sqrt(eta) lambda_is,
  // each stream at the stagnation pressure that velocity implies at P_py,
  // p0 = p_supply pi(lambda_is) / pi(lambda) (= P_py / pi(lambda) in the
  // domain; = p_supply exactly at eta = 1, floors included).
  const D kpi = gamma / gm1;
  auto pi_of = [&](const D& lam) {
    return dpow_dd(dpos_floor(1.0 - (gm1 / gp1) * lam * lam, kRatioFloor), kpi);
  };
  const D lambda1_is = lambda_of(p_g);
  const D lambda2_is = lambda_of(p_e);
  D lambda1 = std::sqrt(eta_primary) * lambda1_is;
  D lambda2 = std::sqrt(eta_secondary) * lambda2_is;
  if (eta_primary != 1.0) {
    p_g = p_g * pi_of(lambda1_is) / pi_of(lambda1);
  }
  if (eta_secondary != 1.0) {
    p_e = p_e * pi_of(lambda2_is) / pi_of(lambda2);
  }
  const D gam = dpos_floor(omega, kRatioFloor);
  D theta21 = t_e / t_g;

  // Kracik & Dvorak mixing (their Eqs. 7-13), identical to the critical path's
  // chain -- only lambda1 is subsonic here instead of supersonic.
  D z1 = lambda1 + 1.0 / lambda1;
  D z2 = lambda2 + 1.0 / lambda2;
  D z3 = (z1 + gam * dsqrt(theta21) * z2) / dsqrt((1.0 + gam) * (1.0 + gam * theta21));
  // Subsonic root (Eq. 12), (z3 - sqrt(z3^2 - 4)) / 2 in the cancellation-free
  // form: a floored lambda makes z3 large, where the difference rounds to 0.
  // z3 >= 2 analytically; the floor only catches rounding at z3 = 2.
  D lambda3 = 2.0 / (z3 + dsqrt(dpos_floor(z3 * z3 - 4.0, kRatioFloor)));

  const D k7 = dpow_dd(gp1 / 2.0, 1.0 / gm1);
  auto q_of = [&](const D& lam) {
    return dpow_dd(dpos_floor(1.0 - (gm1 / gp1) * lam * lam, kRatioFloor), 1.0 / gm1) * k7 * lam;
  };
  D q1 = q_of(lambda1);
  D q2 = q_of(lambda2);
  D q3 = q_of(lambda3);

  return p_g * dsqrt((1.0 + gam) * (1.0 + gam * theta21)) /
         (1.0 + (p_g / p_e) * gam * dsqrt(theta21) * q1 / q2) * q1 / q3;
}

} // namespace

static EjectorJetPumpDischargeJacobian ejector_jetpump_discharge_and_jacobian_raw(
    double p_g_val, double t_g_val, double p_e_val, double t_e_val,
    double p_py_val, double omega_val, double gamma, double recovery_efficiency,
    double eta_primary, double eta_secondary) {
  // Seeds: 0=p_g, 1=t_g, 2=p_e, 3=t_e, 4=p_py, 5=omega.
  using D = DualN<6>;
  const D p03 = jetpump_discharge_d(D::seed(p_g_val, 0), D::seed(t_g_val, 1), D::seed(p_e_val, 2),
                                    D::seed(t_e_val, 3), D::seed(p_py_val, 4),
                                    D::seed(omega_val, 5), D::constant(gamma), eta_primary,
                                    eta_secondary);
  D out = recovery_efficiency * p03;
  return {out.v, out.d[0], out.d[1], out.d[2], out.d[3], out.d[4], out.d[5]};
}

namespace {

// Out of the critical closure's domain -- the secondary's choke pressure at or
// above the primary's nozzle exit pressure (M_py^2 <= 0), or a primary core
// that fills the mixing area (A_sy <= 0) -- the arguments continue through
// dpos_floor: exact wherever the closure is defined. The former hard switch
// M_py = 0 there divided by zero (A_py = inf) and the blended element residual
// went NaN.
constexpr double kMachSqFloor = 1e-8;
constexpr double kAreaFloorFrac = 1e-6;  // of the mixing area ratio

// The critical y-y state (Huang 1999 Eqs. 1-8): the primary's nozzle-exit
// Mach M_p1 (Eq. 2, a bisection in gamma alone), the entrained flow choked at
// y-y (M_sy = 1) at the common pressure P_py = P_sy, the primary core's Mach
// there (Eq. 4 inverted) and the areas (Eqs. 5, 8). One chain for the
// entrainment ratio and both anchor families.
template <int N> struct YYState {
  DualN<N> mach_p1, p_py, mach_py, area_py, area_sy, omega;
};

// fused_exit_factor: Huang's anchor took the nozzle-exit factor
// 1 + (g - 1)/2 M_p1^2 in doubles, the entrainment ratio and Kracik &
// Dvorak's in duals; each keeps its own rounding (revalue).
template <int N>
YYState<N> critical_yy_state(const DualN<N>& p_g, const DualN<N>& t_g, const DualN<N>& p_e,
                             const DualN<N>& t_e, const EjectorGeometry& geom,
                             const DualN<N>& gamma, bool fused_exit_factor = false) {
  using D = DualN<N>;
  const D gm1 = gamma - 1.0;
  const D gp1 = gamma + 1.0;
  const D g_over_gm1 = gamma / gm1;
  YYState<N> s;
  s.mach_p1 = nozzle_exit_mach_d(geom.area_ratio_nozzle, gamma);
  D exit_factor = 1.0 + 0.5 * gm1 * s.mach_p1 * s.mach_p1;
  if (fused_exit_factor) {
    const double gm1v = gm1.v;
    const double mach_p1 = s.mach_p1.v;
    exit_factor = revalue(exit_factor, 1.0 + 0.5 * gm1v * mach_p1 * mach_p1);
  }
  const D p_p1 = ddiv4(p_g, dpow_dd(exit_factor, g_over_gm1));
  s.p_py = ddiv4(p_e, dpow_dd(gp1 / 2.0, g_over_gm1));
  const D core = exit_factor * dpow_dd(ddiv4(s.p_py, p_p1), -gm1 / gamma);
  s.mach_py = dsqrt(dpos_floor(ddiv4(core - 1.0, 0.5 * gm1), kMachSqFloor));
  s.area_py = ejector::phi_p * ddiv4(1.0, s.mach_py) *
              dpow_dd(2.0 / gp1 * (1.0 + 0.5 * gm1 * s.mach_py * s.mach_py), gp1 / (2.0 * gm1));
  s.area_sy = dpos_floor(geom.area_ratio_mix - s.area_py, kAreaFloorFrac * geom.area_ratio_mix);
  s.omega = ddiv4(p_e, p_g) * s.area_sy * dsqrt(ddiv4(t_g, t_e)) *
            std::sqrt(ejector::eta_s / ejector::eta_p);
  return s;
}

// Kracik & Dvorak's mixing closure from the critical y-y state, mixed at the
// critical entrainment ratio (P_c*) or at zero entrainment (the dead-head
// P_b0, their gam = 0): ONE chain, so the two anchors of the subcritical droop
// cannot drift apart (1:1 with _yy_lambda_state + _mixed_flow_stagnation).
// Lossless p03; t03 and lambda3 for the result's value fields.
template <int N> struct MixedStagnation {
  DualN<N> p03, t03, lambda3;
};

template <int N>
MixedStagnation<N> mixed_stagnation_from_yy_d(const DualN<N>& p_g, const DualN<N>& t_g,
                                              const DualN<N>& p_e, const DualN<N>& t_e,
                                              const EjectorGeometry& geom, const DualN<N>& gamma,
                                              const DualN<N>& r_gas, bool dead_head) {
  using D = DualN<N>;
  const D gm1 = gamma - 1.0;
  const D gp1 = gamma + 1.0;
  const YYState<N> s = critical_yy_state(p_g, t_g, p_e, t_e, geom, gamma);

  // Static temperatures and velocities at y-y (Huang Eqs. 9-10, 13-14).
  D t_py = ddiv4(t_g, 1.0 + 0.5 * gm1 * s.mach_py * s.mach_py);
  D t_sy = ddiv4(t_e, gp1 / 2.0);
  D v_py = s.mach_py * dsqrt(gamma * r_gas * t_py);
  D v_sy = dsqrt(gamma * r_gas * t_sy); // M_sy = 1

  D c1_star = dsqrt(2.0 * gamma / gp1 * r_gas * t_g);
  D c2_star = dsqrt(2.0 * gamma / gp1 * r_gas * t_e);
  D lambda1 = ddiv4(v_py, c1_star);
  D lambda2 = ddiv4(v_sy, c2_star); // == 1 identically (M_sy = 1); kept computed for value parity

  D theta21 = ddiv4(t_e, t_g);
  D gam = dead_head ? D::constant(0.0) : s.omega;

  D z1 = lambda1 + ddiv4(1.0, lambda1);
  D z2 = lambda2 + ddiv4(1.0, lambda2);
  D z3 = ddiv4(z1 + gam * dsqrt(theta21) * z2, dsqrt((1.0 + gam) * (1.0 + gam * theta21)));
  // Subsonic root (Eq. 12) in the cancellation-free form 2/(z3 + sqrt(z3^2-4)).
  D lambda3 = ddiv4(2.0, z3 + dsqrt(dpos_floor(z3 * z3 - 4.0, 1e-12)));

  const D k7 = dpow_dd(gp1 / 2.0, 1.0 / gm1);
  auto q_of = [&](const D& lam) {
    D term = dpos_floor(1.0 - (gm1 / gp1) * lam * lam, 1e-12);
    return dpow_dd(term, 1.0 / gm1) * k7 * lam;
  };
  D q1 = q_of(lambda1);
  D q2 = q_of(lambda2); // == 1 identically (lambda2 == 1)
  D q3 = q_of(lambda3);

  // Mixed flow's stagnation pressure and temperature after mixing (Eqs. 7, 11).
  MixedStagnation<N> out;
  out.p03 = ddiv4(ddiv4(p_g * dsqrt((1.0 + gam) * (1.0 + gam * theta21)),
                        1.0 + ddiv4(ddiv4(p_g, p_e) * gam * dsqrt(theta21) * q1, q2)) *
                      q1,
                  q3);
  out.t03 = ddiv4(t_g * (1.0 + gam * theta21), 1.0 + gam);
  out.lambda3 = lambda3;
  return out;
}

// Huang 1999's printed P_c* chain from the critical y-y state (see the
// header): constant-pressure mixing with phi_m, energy, a normal shock in the
// mixed flow C2 through M_m = 1, an isentropic diffuser.
template <int N> struct HuangCritical {
  DualN<N> p_c, t03, mach3_sq;
};

template <int N>
HuangCritical<N> huang_critical_d(const DualN<N>& p_g, const DualN<N>& t_g, const DualN<N>& p_e,
                                  const DualN<N>& t_e, const EjectorGeometry& geom,
                                  const DualN<N>& gamma, const DualN<N>& r_gas,
                                  double phi_m) {
  using D = DualN<N>;
  const D gm1 = gamma - 1.0;
  const D gp1 = gamma + 1.0;
  const YYState<N> s = critical_yy_state(p_g, t_g, p_e, t_e, geom, gamma, true);

  // Static temperatures and velocities at y-y (Huang Eqs. 9-10, 13-14).
  D t_py = ddiv4(t_g, 1.0 + 0.5 * gm1 * s.mach_py * s.mach_py);
  D t_sy = ddiv4(t_e, gp1 / 2.0);
  D v_py = s.mach_py * dsqrt(gamma * r_gas * t_py);
  D v_sy = dsqrt(gamma * r_gas * t_sy);  // M_sy = 1

  // Constant-pressure mixing with phi_m (Eq. 11) and energy (Eq. 12).
  D v_m = ddiv4(phi_m * (v_py + s.omega * v_sy), 1.0 + s.omega);
  const D cp = gamma * r_gas / gm1;
  D t03 = ddiv4(t_g + s.omega * t_e, 1.0 + s.omega);
  D t_m = dpos_floor(t03 - ddiv4(0.5 * v_m * v_m, cp), 1e-6 * t03.v);
  D mach_m_sq = ddiv4(v_m * v_m, gamma * r_gas * t_m);

  // Normal shock in the mixed flow (Eqs. 16-17), C2 through M_m = 1.
  D p3 = s.p_py;
  D mach3_sq = mach_m_sq;
  if (mach_m_sq.v > 1.0) {
    p3 = s.p_py * (1.0 + 2.0 * gamma / gp1 * (mach_m_sq - 1.0));
    mach3_sq = ddiv4(1.0 + 0.5 * gm1 * mach_m_sq, gamma * mach_m_sq - 0.5 * gm1);
  }
  // Isentropic diffuser (Eq. 18).
  HuangCritical<N> out;
  out.p_c = p3 * dpow_dd(1.0 + 0.5 * gm1 * mach3_sq, gamma / gm1);
  out.t03 = t03;
  out.mach3_sq = mach3_sq;
  return out;
}

// The sudden-expansion dead head P_b0 (see the header); t_e does not enter.
template <int N> struct SuddenExpansion {
  DualN<N> p_b0, x;
};

template <int N>
SuddenExpansion<N> sudden_expansion_dead_head_d(const DualN<N>& p_g, const DualN<N>& t_g,
                                                const DualN<N>& p_e, const EjectorGeometry& geom,
                                                const DualN<N>& gamma, const DualN<N>& r_gas) {
  using D = DualN<N>;
  const D gm1 = gamma - 1.0;
  const D gp1 = gamma + 1.0;
  const D k = 0.5 * gm1;
  const double a1 = geom.area_ratio_nozzle;
  const double a3 = geom.area_ratio_mix;

  const D mach_p1 = nozzle_exit_mach_d(a1, gamma);
  // The two gamma-only factors keep the rounding of their double form.
  const double kv = k.v;
  const double gv = gamma.v;
  const double mv = mach_p1.v;
  D p_p1 = ddiv4(p_g, revalue(dpow_dd(1.0 + k * mach_p1 * mach_p1, gamma / gm1),
                              std::pow(1.0 + kv * mv * mv, gv / gm1.v)));
  // Choked primary mass flow per unit throat area (isentropic).
  D m = p_g * dsqrt(ddiv4(gamma, r_gas * t_g)) * dpow_dd(2.0 / gp1, gp1 / (2.0 * gm1));
  D impulse = p_p1 * revalue(a1 * (1.0 + gamma * mach_p1 * mach_p1), a1 * (1.0 + gv * mv * mv)) +
              p_e * (a3 - a1);
  D phi = ddiv4(impulse, m * dsqrt(ddiv4(r_gas * t_g, gamma)));
  // Subsonic root of (k Phi^2 - g^2) X^2 + (Phi^2 - 2g) X - 1 = 0, X = M^2.
  D b = phi * phi - 2.0 * gamma;
  D d = dpos_floor(b * b + 4.0 * (k * phi * phi - gamma * gamma),
                   ejector::dead_head_discriminant_floor);
  D x = ddiv4(2.0, b + dsqrt(d));
  D t3 = ddiv4(t_g, 1.0 + k * x);
  D p3 = ddiv4(m * dsqrt(ddiv4(r_gas * t3, gamma)), a3 * dsqrt(x));
  SuddenExpansion<N> out;
  out.p_b0 = p3 * dpow_dd(1.0 + k * x, gamma / gm1);
  out.x = x;
  return out;
}

// The public critical-path closures seed (p_g, t_g, p_e, t_e) and hold gamma
// and r_gas constant.
using D4 = DualN<4>;
struct Seeds4 {
  D4 p_g, t_g, p_e, t_e;
};
Seeds4 seeds4(double p_g, double t_g, double p_e, double t_e) {
  return {D4::seed(p_g, 0), D4::seed(t_g, 1), D4::seed(p_e, 2), D4::seed(t_e, 3)};
}

} // namespace

static EjectorEntrainmentJacobian ejector_entrainment_ratio_and_jacobian_raw(
    double p_g_val, double t_g_val, double p_e_val, double t_e_val,
    const EjectorGeometry& geom, double gamma) {
  const Seeds4 u = seeds4(p_g_val, t_g_val, p_e_val, t_e_val);
  const YYState<4> s = critical_yy_state(u.p_g, u.t_g, u.p_e, u.t_e, geom, D4::constant(gamma));
  EjectorEntrainmentResult value{
      s.omega.v, s.mach_p1.v, s.mach_py.v, s.area_py.v, s.area_sy.v, s.p_py.v,
  };
  return {value, s.omega.d[0], s.omega.d[1], s.omega.d[2], s.omega.d[3]};
}

EjectorEntrainmentResult ejector_entrainment_ratio(double p_g, double t_g, double p_e,
                                                    double t_e, const EjectorGeometry& geom,
                                                    double gamma) {
  return ejector_entrainment_ratio_and_jacobian(p_g, t_g, p_e, t_e, geom, gamma).value;
}

namespace {

EjectorCriticalPressureJacobian mixed_stagnation_from_yy(
    double p_g_val, double t_g_val, double p_e_val, double t_e_val,
    const EjectorGeometry& geom, double gamma, double r_gas, double recovery_efficiency,
    bool dead_head) {
  const Seeds4 u = seeds4(p_g_val, t_g_val, p_e_val, t_e_val);
  const MixedStagnation<4> m = mixed_stagnation_from_yy_d(
      u.p_g, u.t_g, u.p_e, u.t_e, geom, D4::constant(gamma), D4::constant(r_gas), dead_head);
  const double gm1 = gamma - 1.0;
  const double gp1 = gamma + 1.0;
  const double l3 = m.lambda3.v;
  double mach3_sq = 2.0 * l3 * l3 / (gp1 - gm1 * l3 * l3);
  double mach3 = std::sqrt(mach3_sq);

  EjectorCriticalPressureResult value{
      recovery_efficiency * m.p03.v, m.p03.v, m.t03.v, l3, mach3,
  };
  return {value, recovery_efficiency * m.p03.d[0], recovery_efficiency * m.p03.d[1],
          recovery_efficiency * m.p03.d[2], recovery_efficiency * m.p03.d[3]};
}

} // namespace

static EjectorCriticalPressureJacobian ejector_critical_back_pressure_and_jacobian_raw(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom,
    double gamma, double r_gas, double recovery_efficiency) {
  return mixed_stagnation_from_yy(p_g, t_g, p_e, t_e, geom, gamma, r_gas, recovery_efficiency,
                                  false);
}

static EjectorCriticalPressureJacobian ejector_dead_head_back_pressure_and_jacobian_raw(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom,
    double gamma, double r_gas, double recovery_efficiency) {
  return mixed_stagnation_from_yy(p_g, t_g, p_e, t_e, geom, gamma, r_gas, recovery_efficiency,
                                  true);
}

static double ejector_huang_mixing_coefficient_raw(double area_ratio_mix) {
  if (area_ratio_mix > ejector::phi_m_upper_area_ratio) {
    return ejector::phi_m_large;
  }
  return area_ratio_mix >= ejector::phi_m_lower_area_ratio ? ejector::phi_m_mid
                                                           : ejector::phi_m_small;
}

static EjectorCriticalPressureJacobian ejector_huang_critical_back_pressure_and_jacobian_raw(
    double p_g_val, double t_g_val, double p_e_val, double t_e_val, const EjectorGeometry& geom,
    double gamma, double r_gas, double recovery_efficiency) {
  const double gm1 = gamma - 1.0;
  const double gp1 = gamma + 1.0;
  const Seeds4 u = seeds4(p_g_val, t_g_val, p_e_val, t_e_val);
  const HuangCritical<4> h =
      huang_critical_d(u.p_g, u.t_g, u.p_e, u.t_e, geom, D4::constant(gamma),
                       D4::constant(r_gas), ejector_huang_mixing_coefficient_raw(geom.area_ratio_mix));
  const double mach3_sq = h.mach3_sq.v;
  const double lambda3 = std::sqrt(gp1 * mach3_sq / (2.0 + gm1 * mach3_sq));
  EjectorCriticalPressureResult value{
      recovery_efficiency * h.p_c.v, h.p_c.v, h.t03.v, lambda3, std::sqrt(mach3_sq),
  };
  return {value, recovery_efficiency * h.p_c.d[0], recovery_efficiency * h.p_c.d[1],
          recovery_efficiency * h.p_c.d[2], recovery_efficiency * h.p_c.d[3]};
}

static EjectorCriticalPressureJacobian ejector_sudden_expansion_dead_head_and_jacobian_raw(
    double p_g_val, double t_g_val, double p_e_val, double t_e_val, const EjectorGeometry& geom,
    double gamma, double r_gas, double recovery_efficiency) {
  const double gm1 = gamma - 1.0;
  const double gp1 = gamma + 1.0;
  const Seeds4 u = seeds4(p_g_val, t_g_val, p_e_val, t_e_val);
  const SuddenExpansion<4> e = sudden_expansion_dead_head_d(
      u.p_g, u.t_g, u.p_e, geom, D4::constant(gamma), D4::constant(r_gas));
  const double x = e.x.v;
  const double lambda3 = std::sqrt(gp1 * x / (2.0 + gm1 * x));
  EjectorCriticalPressureResult value{
      recovery_efficiency * e.p_b0.v, e.p_b0.v, t_g_val, lambda3, std::sqrt(x),
  };
  return {value, recovery_efficiency * e.p_b0.d[0], recovery_efficiency * e.p_b0.d[1],
          recovery_efficiency * e.p_b0.d[2], recovery_efficiency * e.p_b0.d[3]};
}

EjectorCriticalPressureResult ejector_critical_back_pressure(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom,
    double gamma, double r_gas, double recovery_efficiency) {
  return ejector_critical_back_pressure_and_jacobian(p_g, t_g, p_e, t_e, geom, gamma, r_gas,
                                                      recovery_efficiency)
      .value;
}

// -----------------------------------------------------------------------------
// Whole-element (f, J): the 4-row operating-regime residual system.
//
// A transliteration of EjectorElement.residuals() (ejector_element.py): the
// closures above evaluated on the element's DualN<11> (the nine unknowns,
// gamma and r_gas), blended by the two smootherstep weights. Kept in 1:1
// correspondence so the two can be cross-validated
// (tests/test_ejector_jacobians.cpp FD-checks this directly).
// -----------------------------------------------------------------------------
namespace {

// _chain: value + sum(partial * grad(input dual)) -- 1:1 with Python _chain.
template <int N>
DualN<N> ejector_chain(double value,
                       std::initializer_list<std::pair<double, DualN<N>>> pairs) {
  DualN<N> r = DualN<N>::constant(value);
  for (const auto& pr : pairs) {
    for (int i = 0; i < N; ++i) r.d[i] += pr.first * pr.second.d[i];
  }
  return r;
}

// C1 clamped smootherstep saturating to EXACTLY 0/1 outside [lo, hi] -- 1:1
// with Python _smootherstep (so s_choke = 1 exactly at the choke threshold and
// critical mode reduces exactly).
template <int N>
DualN<N> ejector_smootherstep(const DualN<N>& t, double lo, double hi) {
  DualN<N> x = (t - lo) / (hi - lo);
  if (x.v <= 0.0) return DualN<N>::constant(0.0);
  if (x.v >= 1.0) return DualN<N>::constant(1.0);
  double val = x.v * x.v * x.v * (x.v * (x.v * 6.0 - 15.0) + 10.0);
  double dval = 30.0 * x.v * x.v * (x.v - 1.0) * (x.v - 1.0);
  return ejector_chain<N>(val, {{dval, x}});
}

// The subcritical droop fades in over (P_b0 - P_c*)/P_c* in [lo, hi].
constexpr double kDroopWindowLo = 1e-3;
constexpr double kDroopWindowHi = 1e-2;

} // namespace

EjectorElementResidualJacobian ejector_element_residuals_and_jacobian(
    double mp_v, double ms_v, double mo_v, double p_g_v, double t_g_v,
    double p_e_v, double t_e_v, double p_out_v, double p_py_v,
    const EjectorGeometry& geom, double area_throat, double area_nozzle_exit,
    double area_secondary, double gamma, double r_gas, double eta_primary,
    double eta_secondary, double recovery_efficiency, double eps_frac,
    double s_choke_lo, double s_choke_hi, double s_sub_lo, double s_sub_hi,
    EjectorAnchorClosure anchor_closure) {
  constexpr const char* who = "ejector_element_residuals_and_jacobian";
  // Newton iterates: finite only (see the validation note at the top).
  for (double x : {mp_v, ms_v, mo_v, p_g_v, t_g_v, p_e_v, t_e_v, p_out_v, p_py_v}) {
    require_finite_arg(who, x, "the flows, pressures and temperatures must be finite");
  }
  check_geometry(who, geom);
  require_positive_arg(who, area_throat, "area_throat must be finite and > 0");
  require_positive_arg(who, area_nozzle_exit, "area_nozzle_exit must be finite and > 0");
  require_positive_arg(who, area_secondary, "area_secondary must be finite and > 0");
  check_gas(who, gamma, r_gas);
  check_efficiency(who, eta_primary, "eta_primary must be finite and in (0, 1]");
  check_efficiency(who, eta_secondary, "eta_secondary must be finite and in (0, 1]");
  check_recovery(who, recovery_efficiency, "recovery_efficiency");
  if (!validation::is_finite_nonnegative(eps_frac)) refuse(who, "eps_frac must be finite and >= 0");
  check_window(who, s_choke_lo, s_choke_hi, "need finite s_choke_lo < s_choke_hi");
  check_window(who, s_sub_lo, s_sub_hi, "need finite s_sub_lo < s_sub_hi");
  using D = DualN<11>;
  const bool huang = anchor_closure == EjectorAnchorClosure::huang1999;

  // Seed the nine Newton unknowns and the gas (see header for the index
  // convention).
  const D mp = D::seed(mp_v, 0);
  const D ms = D::seed(ms_v, 1);
  const D mdot_out = D::seed(mo_v, 2);
  const D p_g = D::seed(p_g_v, 3);
  const D t_g = D::seed(t_g_v, 4);
  const D p_e = D::seed(p_e_v, 5);
  const D t_e = D::seed(t_e_v, 6);
  const D p_out = D::seed(p_out_v, 7);
  const D p_py = D::seed(p_py_v, 8); // owned unknown == mixing-plane static
  const D gam = D::seed(gamma, 9);
  const D rg = D::seed(r_gas, 10);

  // s_choke = smootherstep(mp / cap), cap = choked_mass_flow(Pt_p).
  const D cap = choked_mass_flow_d(p_g, t_g, area_throat, gam, rg, eta_primary);
  D s = ejector_smootherstep<11>(mp / cap, s_choke_lo, s_choke_hi);

  // R0: primary C-D nozzle (choked/unchoked via the closure's own smooth-min).
  const D cd_d =
      cd_nozzle_d(p_g, t_g, p_py, area_throat, area_nozzle_exit, gam, rg, eta_primary, eps_frac);
  D r0 = mp - cd_d;

  // Blend the critical (choked-primary) and jet-pump (unchoked) closures by
  // s = s_choke. The critical closures return NaN for an unchoked primary, and
  // the jet-pump closure is symmetric at the choked end; smootherstep is flat
  // at s in {0, 1} so the inactive branch contributes zero, but 0 * NaN = NaN
  // would still poison the value/Jacobian -- so evaluate each branch ONLY where
  // its weight is nonzero, leaving a hard-zero dual otherwise (C1-consistent
  // with the flat endpoint).
  // ACTUAL ratio for the discharge (R1 pins ms itself). A probe can drive
  // either flow to zero or reverse it: both are floored C1 at 1e-9 of the
  // choked flow, so omega stays finite and positive (#481). The floor moves
  // with the choked flow, and so with every input of it.
  const D m_floor = 1e-9 * cap;
  D omega_eff = dpos_floor_dd(ms, m_floor) / dpos_floor_dd(mp, m_floor);
  const bool crit_active = s.v > 0.0;
  const bool jet_active = s.v < 1.0;

  D omega_crit = D::constant(0.0);
  D p03_crit = D::constant(0.0);
  D s_sub = D::constant(1.0);
  // Subcritical droop (Phase 1A, #619): between P_c* and the dead-head P_b0
  // the choked ejector's entrainment falls linearly to 0 (omega_choked) and the
  // outlet floats; the outlet pin switches in at P_b0 (s_b0), not at P_c*.
  D omega_choked = D::constant(0.0);
  D p03_dead = D::constant(0.0);
  D s_b0 = D::constant(1.0);
  if (crit_active) {
    omega_crit = critical_yy_state(p_g, t_g, p_e, t_e, geom, gam).omega;
    p03_crit = huang ? huang_critical_d(p_g, t_g, p_e, t_e, geom, gam, rg,
                                        ejector_huang_mixing_coefficient_raw(geom.area_ratio_mix))
                           .p_c
                     : mixed_stagnation_from_yy_d(p_g, t_g, p_e, t_e, geom, gam, rg, false).p03;
    s_sub = ejector_smootherstep<11>(p_out / (recovery_efficiency * p03_crit), s_sub_lo, s_sub_hi);

    // The droop needs a window P_b0 > P_c*; it collapses in the M_py < 1
    // corner (an unchoked primary, the jet pump's). d fades the droop out
    // there, back to the plateau-and-pin-at-P_c* form, continuously: P_b0
    // depends on Newton unknowns, so a hard switch would be a jump in F.
    const D p03_b0 = huang ? sudden_expansion_dead_head_d(p_g, t_g, p_e, geom, gam, rg).p_b0
                           : mixed_stagnation_from_yy_d(p_g, t_g, p_e, t_e, geom, gam, rg, true).p03;
    D window = (p03_b0 - p03_crit) / p03_crit;
    D d = ejector_smootherstep<11>(window, kDroopWindowLo, kDroopWindowHi);
    D s_b0_live = D::constant(0.0);
    D omega_sub = omega_crit;
    if (d.v > 0.0) {
      // Tier-1 chord (P_c*, omega_crit) -> (P_b0, 0) in the outlet Pt;
      // window >= kDroopWindowLo keeps the division bounded.
      omega_sub = omega_crit * (recovery_efficiency * p03_b0 - p_out) /
                  (recovery_efficiency * (p03_b0 - p03_crit));
      s_b0_live = ejector_smootherstep<11>(p_out / (recovery_efficiency * p03_b0), s_sub_lo,
                                           s_sub_hi);
    }
    // Plateau below P_c* exactly (s_sub = 0), the chord above it: critical
    // results are unchanged, unlike a smooth-min, which shaves omega_crit.
    omega_choked = omega_crit + d * s_sub * (omega_sub - omega_crit);
    s_b0 = d * s_b0_live + (1.0 - d) * s_sub;
    p03_dead = d * p03_b0 + (1.0 - d) * p03_crit;
  }

  D ms_sec = D::constant(0.0);
  D p03_jet = D::constant(0.0);
  if (jet_active) {
    ms_sec = cd_nozzle_d(p_e, t_e, p_py, area_secondary, area_secondary, gam, rg, eta_secondary,
                         eps_frac);
    p03_jet = jetpump_discharge_d(p_g, t_g, p_e, t_e, p_py, omega_eff, gam, eta_primary,
                                  eta_secondary);
  }

  // R1: blended entrainment  ms = s*omega_crit*mp + (1-s)*ms_sec.
  D ms_model = s * omega_choked * mp + (1.0 - s) * ms_sec;
  D r1 = ms - ms_model;

  // R2: mass conservation.
  D r2 = mdot_out - mp - ms;

  // R3: regime-dependent closure of the owned unknown (see ejector_element.py
  // for the w_pin table). discharge = recovery * blend(P_c*, P03_jet) by
  // s_choke; w_pin = 1 - s*(1 - s_sub) is the weight on the outlet-pin form.
  // Choked: P_py = recovery*P_c* (a diagnostic) while the outlet floats,
  // below P_b0; the outlet pinned at recovery*P_b0 above it, which excludes
  // the spurious double-choked root (Pt_primary below the back pressure is
  // only possible above P_b0, since P_b0 < Pt_primary).
  D discharge = recovery_efficiency *
                (s * ((1.0 - s_b0) * p03_crit + s_b0 * p03_dead) + (1.0 - s) * p03_jet);
  D w_pin = 1.0 - s * (1.0 - s_b0);
  D lhs = w_pin * p_out + (1.0 - w_pin) * p_py;
  D r3 = lhs - discharge;

  EjectorElementResidualJacobian out;
  // The closures above take M_p1 from the same bisection; report its clamp.
  out.nozzle_exit_mach_out_of_range =
      ejector_nozzle_exit_mach_raw(geom.area_ratio_nozzle, gamma).out_of_range;
  const D rows[4] = {r0, r1, r2, r3};
  for (int i = 0; i < 4; ++i) {
    out.residuals[i] = rows[i].v;
    for (int j = 0; j < 9; ++j) out.jacobian[i][j] = rows[i].d[j];
    out.d_dgamma[i] = rows[i].d[9];
    out.d_dr_gas[i] = rows[i].d[10];
  }
  return out;
}

EjectorJetPumpIncompressibleResult ejector_jetpump_incompressible_residuals_and_jacobian(
    double mp, double ms, double mo, double pt_p, double v_p, double pt_s,
    double v_s, double pt_o, double p3, double area_nozzle_exit,
    double area_secondary, double area_mixing, double k_nozzle,
    double k_suction, double k_throat, double k_diffuser, double recovery) {
  constexpr const char* who = "ejector_jetpump_incompressible_residuals_and_jacobian";
  // Newton iterates (v = 1/rho from the port closures included): finite only.
  for (double x : {mp, ms, mo, pt_p, v_p, pt_s, v_s, pt_o, p3}) {
    require_finite_arg(who, x, "the flows, pressures and specific volumes must be finite");
  }
  require_positive_arg(who, area_nozzle_exit, "area_nozzle_exit must be finite and > 0");
  require_positive_arg(who, area_secondary, "area_secondary must be finite and > 0");
  require_positive_arg(who, area_mixing, "area_mixing must be finite and > 0");
  for (double k : {k_nozzle, k_suction, k_throat, k_diffuser}) {
    if (!validation::is_finite_nonnegative(k)) {
      refuse(who, "k_nozzle, k_suction, k_throat and k_diffuser must be finite and >= 0");
    }
  }
  check_recovery(who, recovery, "recovery");
  using D = DualN<9>;
  const D Mp = D::seed(mp, 0), Ms = D::seed(ms, 1), Mo = D::seed(mo, 2);
  const D Ptp = D::seed(pt_p, 3), Vp = D::seed(v_p, 4);
  const D Pts = D::seed(pt_s, 5), Vs = D::seed(v_s, 6);
  const D Pto = D::seed(pt_o, 7), P3 = D::seed(p3, 8);
  const double Ae = area_nozzle_exit, As = area_secondary, A3 = area_mixing;

  // m|m| as m * m * sgn(m): exact derivative 2|m|, C1 through zero.
  auto msq = [](const D &m) { return m * m * (m.v >= 0.0 ? 1.0 : -1.0); };
  const D Mt = Mp + Ms;
  const double sgn_t = (Mt.v >= 0.0) ? 1.0 : -1.0;
  // mt|mt| v_m = |mt| (mp v_p + ms v_s)
  const D flux_t = Mt * sgn_t * (Mp * Vp + Ms * Vs);

  const D R0 = Ptp - P3 - (1.0 + k_nozzle) * 0.5 * Vp * msq(Mp) / (Ae * Ae);
  const D R1 = Pts - P3 - (1.0 + k_suction) * 0.5 * Vs * msq(Ms) / (As * As);
  const D R2 = Mo - Mp - Ms;
  const D q4 = 0.5 * flux_t / (A3 * A3);
  const D p4 = P3 + (Vp * msq(Mp) / Ae + Vs * msq(Ms) / As - flux_t / A3) / A3 - k_throat * q4;
  const D R3 = Pto - recovery * (p4 + (1.0 - k_diffuser) * q4);

  EjectorJetPumpIncompressibleResult out;
  const D rows[4] = {R0, R1, R2, R3};
  for (int r = 0; r < 4; ++r) {
    out.residuals[r] = rows[r].v;
    for (int k = 0; k < 9; ++k) {
      out.jacobian[r][k] = rows[r].d[k];
    }
  }
  out.p4 = p4.v;
  out.q4 = q4.v;
  return out;
}

// -----------------------------------------------------------------------------
// Checked public entry points of the closures the element calls unchecked.
// -----------------------------------------------------------------------------

double ejector_area_over_astar(double mach, double gamma) {
  require_positive_arg("ejector_area_over_astar", mach, "mach must be finite and > 0");
  check_gamma("ejector_area_over_astar", gamma);
  return ejector_area_over_astar_raw(mach, gamma);
}

double ejector_mach_from_area_ratio_supersonic(double area_ratio, double gamma) {
  constexpr const char* who = "ejector_mach_from_area_ratio_supersonic";
  if (!validation::is_finite(area_ratio) || !(area_ratio >= 1.0)) {
    refuse(who, "area_ratio must be finite and >= 1");
  }
  check_gamma(who, gamma);
  return ejector_mach_from_area_ratio_supersonic_raw(area_ratio, gamma);
}

EjectorNozzleExitMach ejector_nozzle_exit_mach(double area_ratio, double gamma) {
  constexpr const char* who = "ejector_nozzle_exit_mach";
  if (!validation::is_finite(area_ratio) || !(area_ratio >= 1.0)) {
    refuse(who, "area_ratio must be finite and >= 1");
  }
  check_gamma(who, gamma);
  return ejector_nozzle_exit_mach_raw(area_ratio, gamma);
}

namespace {
void check_choked(const char* who, double Pt, double Tt, double area_throat, double gamma,
                  double r_gas, double eta) {
  require_positive_arg(who, Pt, "Pt must be finite and > 0");
  require_positive_arg(who, Tt, "Tt must be finite and > 0");
  require_positive_arg(who, area_throat, "area_throat must be finite and > 0");
  check_gas(who, gamma, r_gas);
  check_efficiency(who, eta, "eta must be finite and in (0, 1]");
}

void check_closure(const char* who, double p_g, double t_g, double p_e, double t_e,
                   const EjectorGeometry& geom, double gamma, double r_gas,
                   double recovery_efficiency) {
  check_supplies(who, p_g, t_g, p_e, t_e);
  check_geometry(who, geom);
  check_gas(who, gamma, r_gas);
  check_recovery(who, recovery_efficiency, "recovery_efficiency");
}
} // namespace

double ejector_choked_mass_flow(double Pt, double Tt, double area_throat, double gamma,
                                double r_gas, double eta) {
  check_choked("ejector_choked_mass_flow", Pt, Tt, area_throat, gamma, r_gas, eta);
  return ejector_choked_mass_flow_raw(Pt, Tt, area_throat, gamma, r_gas, eta);
}

std::tuple<double, double, double>
ejector_choked_mass_flow_and_jacobian(double Pt, double Tt, double area_throat, double gamma,
                                      double r_gas, double eta) {
  check_choked("ejector_choked_mass_flow_and_jacobian", Pt, Tt, area_throat, gamma, r_gas, eta);
  return ejector_choked_mass_flow_and_jacobian_raw(Pt, Tt, area_throat, gamma, r_gas, eta);
}

// p_static_down above Pt is reverse flow (#619), so it is only checked > 0.
double ejector_cd_nozzle_mass_flow(double Pt, double Tt, double p_static_down,
                                   double area_throat, double area_exit, double gamma,
                                   double r_gas, double eta, double eps_frac) {
  constexpr const char* who = "ejector_cd_nozzle_mass_flow";
  check_nozzle(who, Pt, Tt, area_throat, area_exit, gamma, r_gas, eta, eps_frac);
  require_positive_arg(who, p_static_down, "p_static_down must be finite and > 0");
  return ejector_cd_nozzle_mass_flow_raw(Pt, Tt, p_static_down, area_throat, area_exit, gamma,
                                         r_gas, eta, eps_frac);
}

EjectorCDNozzleJacobian ejector_cd_nozzle_mass_flow_and_jacobian(
    double Pt, double Tt, double p_py, double area_throat, double area_exit, double gamma,
    double r_gas, double eta, double eps_frac) {
  constexpr const char* who = "ejector_cd_nozzle_mass_flow_and_jacobian";
  check_nozzle(who, Pt, Tt, area_throat, area_exit, gamma, r_gas, eta, eps_frac);
  require_positive_arg(who, p_py, "p_py must be finite and > 0");
  return ejector_cd_nozzle_mass_flow_and_jacobian_raw(Pt, Tt, p_py, area_throat, area_exit,
                                                      gamma, r_gas, eta, eps_frac);
}

EjectorJetPumpDischargeJacobian ejector_jetpump_discharge_and_jacobian(
    double p_g, double t_g, double p_e, double t_e, double p_py, double omega, double gamma,
    double recovery_efficiency, double eta_primary, double eta_secondary) {
  constexpr const char* who = "ejector_jetpump_discharge_and_jacobian";
  check_supplies(who, p_g, t_g, p_e, t_e);
  require_positive_arg(who, p_py, "p_py must be finite and > 0");
  if (!validation::is_finite_nonnegative(omega)) refuse(who, "omega must be finite and >= 0");
  check_gamma(who, gamma);
  check_recovery(who, recovery_efficiency, "recovery_efficiency");
  check_efficiency(who, eta_primary, "eta_primary must be finite and in (0, 1]");
  check_efficiency(who, eta_secondary, "eta_secondary must be finite and in (0, 1]");
  return ejector_jetpump_discharge_and_jacobian_raw(p_g, t_g, p_e, t_e, p_py, omega, gamma,
                                                    recovery_efficiency, eta_primary,
                                                    eta_secondary);
}

EjectorEntrainmentJacobian ejector_entrainment_ratio_and_jacobian(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom, double gamma) {
  constexpr const char* who = "ejector_entrainment_ratio_and_jacobian";
  check_supplies(who, p_g, t_g, p_e, t_e);
  check_geometry(who, geom);
  check_gamma(who, gamma);
  return ejector_entrainment_ratio_and_jacobian_raw(p_g, t_g, p_e, t_e, geom, gamma);
}

EjectorCriticalPressureJacobian ejector_critical_back_pressure_and_jacobian(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom, double gamma,
    double r_gas, double recovery_efficiency) {
  check_closure("ejector_critical_back_pressure_and_jacobian", p_g, t_g, p_e, t_e, geom, gamma,
                r_gas, recovery_efficiency);
  return ejector_critical_back_pressure_and_jacobian_raw(p_g, t_g, p_e, t_e, geom, gamma, r_gas,
                                                         recovery_efficiency);
}

EjectorCriticalPressureJacobian ejector_dead_head_back_pressure_and_jacobian(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom, double gamma,
    double r_gas, double recovery_efficiency) {
  check_closure("ejector_dead_head_back_pressure_and_jacobian", p_g, t_g, p_e, t_e, geom, gamma,
                r_gas, recovery_efficiency);
  return ejector_dead_head_back_pressure_and_jacobian_raw(p_g, t_g, p_e, t_e, geom, gamma, r_gas,
                                                          recovery_efficiency);
}

double ejector_huang_mixing_coefficient(double area_ratio_mix) {
  if (!validation::is_finite(area_ratio_mix) || !(area_ratio_mix > 1.0)) {
    refuse("ejector_huang_mixing_coefficient", "area_ratio_mix must be finite and > 1");
  }
  return ejector_huang_mixing_coefficient_raw(area_ratio_mix);
}

EjectorCriticalPressureJacobian ejector_huang_critical_back_pressure_and_jacobian(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom, double gamma,
    double r_gas, double recovery_efficiency) {
  check_closure("ejector_huang_critical_back_pressure_and_jacobian", p_g, t_g, p_e, t_e, geom,
                gamma, r_gas, recovery_efficiency);
  return ejector_huang_critical_back_pressure_and_jacobian_raw(p_g, t_g, p_e, t_e, geom, gamma,
                                                               r_gas, recovery_efficiency);
}

EjectorCriticalPressureJacobian ejector_sudden_expansion_dead_head_and_jacobian(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom, double gamma,
    double r_gas, double recovery_efficiency) {
  // t_e is unused (the secondary is dead-headed) but still a supply state.
  check_closure("ejector_sudden_expansion_dead_head_and_jacobian", p_g, t_g, p_e, t_e, geom,
                gamma, r_gas, recovery_efficiency);
  return ejector_sudden_expansion_dead_head_and_jacobian_raw(p_g, t_g, p_e, t_e, geom, gamma,
                                                             r_gas, recovery_efficiency);
}

} // namespace combaero::solver
