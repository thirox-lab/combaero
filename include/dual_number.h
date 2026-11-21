#pragma once

// -----------------------------------------------------------------------------
// Forward-mode dual numbers with N seeds, for whole-element analytic (f, J).
//
// The repo's solver rule is that a solver-facing calculation exposes value AND
// derivative from one C++ evaluation, with the derivative obtained by the chain
// rule rather than by finite differences. A dual number carries both: the
// primal value plus one partial per Newton unknown, and every operator applies
// the chain rule as it goes, so a routine written once against DualN<N> yields
// the exact Jacobian rows with no separate derivation to keep in sync.
//
// Extracted from src/ejector.cpp, where it was a translation-unit-local helper
// for the ejector's operating-regime closures, so the junction closure can use
// the same one (issue #271). The ejector's critical-mode path, once on a
// fixed-seed Dual4, is on DualN too since #1246 (gamma and r_gas became seeds).
//
// Seeding: DualN<N>::seed(v, i) is the unknown whose partial index is i;
// DualN<N>::constant(c) is anything the Newton step does not vary.
//
// Branch-on-primal. Several closures choose a branch from a sign, a quadrant or
// a wrap count. Those decisions are taken on the PRIMAL value (`.v`) and the
// derivative is then exact within the chosen branch and one-sided at the
// boundary. That is the right behaviour and is much better than a zero column,
// but it means the Jacobian is a one-sided derivative exactly at a switching
// surface -- document it wherever such a branch is taken.
// -----------------------------------------------------------------------------

#include <array>
#include <cmath>

#include "math_constants.h"

namespace combaero::solver {

template <int N> struct DualN {
  double v = 0.0;
  std::array<double, N> d{};

  static DualN constant(double c) {
    DualN r;
    r.v = c;
    return r;
  }
  static DualN seed(double c, int i) {
    DualN r;
    r.v = c;
    r.d[i] = 1.0;
    return r;
  }
};

// --- arithmetic --------------------------------------------------------------

template <int N> DualN<N> operator+(const DualN<N>& a, const DualN<N>& b) {
  DualN<N> r;
  r.v = a.v + b.v;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] + b.d[i];
  return r;
}
template <int N> DualN<N> operator+(const DualN<N>& a, double c) {
  DualN<N> r = a;
  r.v += c;
  return r;
}
template <int N> DualN<N> operator+(double c, const DualN<N>& a) { return a + c; }

template <int N> DualN<N> operator-(const DualN<N>& a, const DualN<N>& b) {
  DualN<N> r;
  r.v = a.v - b.v;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] - b.d[i];
  return r;
}
template <int N> DualN<N> operator-(const DualN<N>& a, double c) {
  DualN<N> r = a;
  r.v -= c;
  return r;
}
template <int N> DualN<N> operator-(double c, const DualN<N>& a) {
  DualN<N> r;
  r.v = c - a.v;
  for (int i = 0; i < N; ++i) r.d[i] = -a.d[i];
  return r;
}
template <int N> DualN<N> operator-(const DualN<N>& a) { return 0.0 - a; }

template <int N> DualN<N> operator*(const DualN<N>& a, const DualN<N>& b) {
  DualN<N> r;
  r.v = a.v * b.v;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * b.v + a.v * b.d[i];
  return r;
}
template <int N> DualN<N> operator*(const DualN<N>& a, double c) {
  DualN<N> r;
  r.v = a.v * c;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * c;
  return r;
}
template <int N> DualN<N> operator*(double c, const DualN<N>& a) { return a * c; }

template <int N> DualN<N> operator/(const DualN<N>& a, const DualN<N>& b) {
  DualN<N> r;
  r.v = a.v / b.v;
  double inv = 1.0 / b.v;
  for (int i = 0; i < N; ++i) r.d[i] = (a.d[i] - r.v * b.d[i]) * inv;
  return r;
}
template <int N> DualN<N> operator/(const DualN<N>& a, double c) { return a * (1.0 / c); }
template <int N> DualN<N> operator/(double c, const DualN<N>& a) {
  DualN<N> r;
  r.v = c / a.v;
  double coef = -c / (a.v * a.v);
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * coef;
  return r;
}

// --- elementary functions ----------------------------------------------------

template <int N> DualN<N> dsqrt(const DualN<N>& a) {
  DualN<N> r;
  r.v = std::sqrt(a.v);
  double coef = 0.5 / r.v;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * coef;
  return r;
}
template <int N> DualN<N> dpow(const DualN<N>& a, double c) {
  DualN<N> r;
  r.v = std::pow(a.v, c);
  double coef = c * std::pow(a.v, c - 1.0);
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * coef;
  return r;
}

// d/dx exp(u) = exp(u) u'
template <int N> DualN<N> dexp(const DualN<N>& a) {
  DualN<N> r;
  r.v = std::exp(a.v);
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * r.v;
  return r;
}
// d/dx log(u) = u'/u
// C1 positive floor: a itself for a >= eps, eps^2 / (2 eps - a) below. Value
// and slope match at eps; the continuation is positive, monotone, and decays
// only like eps^2/|a|, so it never underflows to 0 (an exponential tail with
// the small eps a domain guard needs does, within a few units of a). Guards a
// closure's sqrt/pow argument that a solver probe can push out of its domain:
// exact inside it, finite and C1 outside (#481), instead of a NaN that
// poisons a whole residual. The same form as the wall htc knee.
template <int N> DualN<N> dpos_floor(const DualN<N>& a, double eps) {
  if (a.v >= eps) {
    return a;
  }
  DualN<N> r;
  const double den = 2.0 * eps - a.v;
  r.v = eps * eps / den;
  const double slope = r.v / den;  // eps^2 / den^2
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * slope;
  return r;
}

template <int N> DualN<N> dlog(const DualN<N>& a) {
  DualN<N> r;
  r.v = std::log(a.v);
  double inv = 1.0 / a.v;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * inv;
  return r;
}
// d/dx log(1 + u) = u'/(1 + u): log near 1 without the cancellation of
// log(1 + u) for |u| ~ machine epsilon.
template <int N> DualN<N> dlog1p(const DualN<N>& a) {
  DualN<N> r;
  r.v = std::log1p(a.v);
  double inv = 1.0 / (1.0 + a.v);
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * inv;
  return r;
}
// d/dx (exp(u) - 1) = exp(u) u'
template <int N> DualN<N> dexpm1(const DualN<N>& a) {
  DualN<N> r;
  r.v = std::expm1(a.v);
  const double e = std::exp(a.v);
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * e;
  return r;
}
// d/dx sin(u) = cos(u) u'
template <int N> DualN<N> dsin(const DualN<N>& a) {
  DualN<N> r;
  r.v = std::sin(a.v);
  double coef = std::cos(a.v);
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * coef;
  return r;
}
// d/dx cos(u) = -sin(u) u'
template <int N> DualN<N> dcos(const DualN<N>& a) {
  DualN<N> r;
  r.v = std::cos(a.v);
  double coef = -std::sin(a.v);
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * coef;
  return r;
}

// d atan2(y, x) = (x dy - y dx) / (x^2 + y^2)
//
// The quadrant is resolved by std::atan2 on the primal values, so the result
// is continuous in the derivative everywhere the function itself is (that is,
// away from the negative-x axis where atan2 jumps by 2*pi -- the jump is a
// constant, so the derivative below stays correct across it).
template <int N> DualN<N> datan2(const DualN<N>& y, const DualN<N>& x) {
  DualN<N> r;
  r.v = std::atan2(y.v, x.v);
  double denom = x.v * x.v + y.v * y.v;
  double coef = 1.0 / denom;
  for (int i = 0; i < N; ++i) r.d[i] = (x.v * y.d[i] - y.v * x.d[i]) * coef;
  return r;
}

// |u|, with the derivative taken on the primal sign. A BRANCH-ON-PRIMAL: at
// u = 0 the one-sided derivative +u' is returned rather than 0, which keeps the
// Newton step informative on the seam instead of dropping the column.
template <int N> DualN<N> dabs(const DualN<N>& a) {
  DualN<N> r;
  r.v = std::abs(a.v);
  double sign = (a.v < 0.0) ? -1.0 : 1.0;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * sign;
  return r;
}

// --- angle wrapping ----------------------------------------------------------
//
// Matlab's wrapToPi / wrapTo2Pi, which the Mynard junction closure uses. Both
// shift the value by an integer multiple of 2*pi, so within a branch the
// derivative is the identity and the partials pass through untouched. Which
// multiple is a BRANCH-ON-PRIMAL, and the derivative is one-sided exactly at a
// wrap boundary.

// To [-pi, pi). NOTE the closed end: at exactly +pi this returns -pi, matching
// the Python reference's `(x + pi) % (2*pi) - pi`. Matlab's wrapToPi is
// (-pi, pi] and returns +pi there. The reference is what the port must match.
template <int N> DualN<N> dwrap_to_pi(const DualN<N>& a) {
  DualN<N> r = a;
  r.v = std::fmod(a.v + M_PI, 2.0 * M_PI);
  if (r.v < 0.0) r.v += 2.0 * M_PI;
  r.v -= M_PI;
  return r;
}

// To [0, 2*pi).
template <int N> DualN<N> dwrap_to_2pi(const DualN<N>& a) {
  DualN<N> r = a;
  r.v = std::fmod(a.v, 2.0 * M_PI);
  if (r.v < 0.0) r.v += 2.0 * M_PI;
  return r;
}

} // namespace combaero::solver
