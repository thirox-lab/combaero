#include "../include/thermo.h"
#include "../include/composition.h"
#include "../include/humidair.h"
#include "../include/thermo_c2_refit_data.h"
#include "../include/thermo_transport_data.h"
#include "../include/transport.h"
#include "../include/input_validation.h"
#include "../include/inner_solve.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace combaero {

double J_per_mol_to_J_per_kg(double value, double molar_mass) {
  return value * 1000.0 /
         molar_mass; // molar_mass is in g/mol, convert to kg/mol
}

// -------------------------------------------------------------
// Species metadata and lookup functions
// -------------------------------------------------------------

std::string species_name(std::size_t idx) {
    if (idx >= species_names.size()) {
        throw std::out_of_range("species_name: index out of bounds");
    }
    return species_names[idx];
}

std::size_t species_index_from_name(const std::string& name) {
    auto it = species_index.find(name);
    if (it == species_index.end()) {
        throw std::invalid_argument("species_index_from_name: species not found: " + name);
    }
    return static_cast<std::size_t>(it->second);
}

std::size_t num_species() {
    return species_names.size();
}

double species_molar_mass(std::size_t idx) {
    if (idx >= molar_masses.size()) {
        throw std::out_of_range("species_molar_mass: index out of bounds");
    }
    return molar_masses[idx];
}

double species_molar_mass_from_name(const std::string& name) {
    return species_molar_mass(species_index_from_name(name));
}


// -------------------------------------------------------------
// NASA 9 Polynomial evaluation
// -------------------------------------------------------------

namespace {
const NASA9_Interval& get_interval(std::size_t species_idx, double T) {
    const auto& nasa = nasa_coeffs[species_idx];
    for (const auto& interval : nasa.intervals) {
        if (T >= interval.T_min && T <= interval.T_max) {
            return interval;
        }
    }
    // Fallback to nearest if out of range
    if (T < nasa.intervals.front().T_min) return nasa.intervals.front();
    return nasa.intervals.back();
}

// --- the C2 re-fit (#549) ----------------------------------------------------

std::atomic<int> g_thermo_model{static_cast<int>(ThermoModel::NASA9_C2)};

// The bridge covering T for this species, or nullptr (NASA9 model, or T
// outside every window).
const C2Bridge* active_bridge(std::size_t species_idx, double T) {
  if (g_thermo_model.load(std::memory_order_relaxed) !=
      static_cast<int>(ThermoModel::NASA9_C2)) {
    return nullptr;
  }
  for (const auto& b : c2_bridges[species_idx]) {
    if (std::abs(T - b.T_join) < b.half_width) {
      return &b;
    }
  }
  return nullptr;
}

double bridge_poly(const C2Bridge& b, double t) {
  double acc = 0.0;
  for (std::size_t k = b.c.size(); k-- > 0;) {
    acc = acc * t + b.c[k];
  }
  return acc;
}

double bridge_poly_dt(const C2Bridge& b, double t) {
  double acc = 0.0;
  for (std::size_t k = b.c.size(); k-- > 1;) {
    acc = acc * t + static_cast<double>(k) * b.c[k];
  }
  return acc;
}

double bridge_poly_dt2(const C2Bridge& b, double t) {
  double acc = 0.0;
  for (std::size_t k = b.c.size(); k-- > 2;) {
    acc = acc * t + static_cast<double>(k * (k - 1)) * b.c[k];
  }
  return acc;
}

// Antiderivative in t of the bridge polynomial.
double bridge_poly_int(const C2Bridge& b, double t) {
  double acc = 0.0;
  for (std::size_t k = b.c.size(); k-- > 0;) {
    acc = acc * t + b.c[k] / static_cast<double>(k + 1);
  }
  return acc * t;
}

// The NASA-9 interval that ends at the bridge's join (its left side).
const NASA9_Interval& left_interval(std::size_t species_idx, double T_join) {
  for (const auto& iv : nasa_coeffs[species_idx].intervals) {
    if (iv.T_max == T_join) {
      return iv;
    }
  }
  return nasa_coeffs[species_idx].intervals.front();
}

// NASA-9 H/R [K] and S/R of one interval (the published closed forms).
double nasa_H_R(const std::array<double, 10>& a, double T) {
  return -a[0] / T + a[1] * std::log(T) + a[2] * T + a[3] * T * T / 2.0 +
         a[4] * T * T * T / 3.0 + a[5] * T * T * T * T / 4.0 +
         a[6] * T * T * T * T * T / 5.0 + a[8];
}

double nasa_S_R(const std::array<double, 10>& a, double T) {
  return -a[0] / (2.0 * T * T) - a[1] / T + a[2] * std::log(T) + a[3] * T +
         a[4] * T * T / 2.0 + a[5] * T * T * T / 3.0 +
         a[6] * T * T * T * T / 4.0 + a[9];
}

// H/R [K] inside a bridge: NASA-9 at the window's left edge plus the
// bridge's integral of cp/R dT, closed form.
double bridge_H_R(std::size_t species_idx, const C2Bridge& b, double T) {
  const double w = b.half_width;
  const double t = (T - b.T_join) / w;
  const auto& L = left_interval(species_idx, b.T_join).coeffs;
  return nasa_H_R(L, b.T_join - w) +
         w * (bridge_poly_int(b, t) - bridge_poly_int(b, -1.0));
}

// S/R inside a bridge: NASA-9 at the left edge plus the integral of cp/(R T)
// dT over [-1, t], by 8-point Gauss-Legendre. The integrand is a polynomial
// over (T_join/w + t), analytic with its pole at t = -20: the rule is exact
// to rounding.
double bridge_S_R(std::size_t species_idx, const C2Bridge& b, double T) {
  static const double xg[4] = {0.1834346424956498, 0.5255324099163290,
                               0.7966664774136267, 0.9602898564975363};
  static const double wg[4] = {0.3626837833783620, 0.3137066458778873,
                               0.2223810344533745, 0.1012285362903763};
  const double w = b.half_width;
  const double t = (T - b.T_join) / w;
  const double c = 0.5 * (t - 1.0), r = 0.5 * (t + 1.0);
  auto f = [&](double tau) { return bridge_poly(b, tau) / (b.T_join + w * tau); };
  double acc = 0.0;
  for (int i = 0; i < 4; ++i) {
    acc += wg[i] * (f(c - r * xg[i]) + f(c + r * xg[i]));
  }
  const auto& L = left_interval(species_idx, b.T_join).coeffs;
  return nasa_S_R(L, b.T_join - w) + w * r * acc;
}

// --- the tails beyond the data (#606) ----------------------------------------

// One end's tail. t = (T - T0) / w >= 0 on its side (w < 0 below the data);
// cp/R = p(t) on t < 1, C2 with NASA-9 at T0 and with p' = p'' = 0 at t = 1,
// and cp_inf = p(1) beyond. Below the data w = -T0 (t = 1 is 0 K) and p is the
// quintic that also takes p(1) to the rigid-rotor value; above it the quartic
// (c5 = 0) that levels off wherever it reaches.
struct Tail {
  double T0 = 0.0;
  double w = 0.0;
  std::array<double, 6> c{};  // p(t), ascending powers
  double H0 = 0.0;            // H/R [K] and S/R of the data at T0
  double S0 = 0.0;
  double cp_inf = 0.0;
  double T_end = 0.0;  // T0 + w, where the fade ends
  double H_end = 0.0;
  double S_end = 0.0;
  // p(t) = (t - t_pole) Q(t) + r with t_pole = -T0 / w (T = 0): S/R in closed
  // form, since cp/(R T) dT = (Q(t) + r / (t - t_pole)) dt.
  std::array<double, 5> Q{};
  double r = 0.0;
};

struct SpeciesTails {
  Tail lo;
  Tail hi;
  double T_lo = 0.0;  // the data's own range
  double T_hi = 0.0;
};

double tail_p(const Tail& tl, double t) {
  double acc = 0.0;
  for (std::size_t k = tl.c.size(); k-- > 0;) {
    acc = acc * t + tl.c[k];
  }
  return acc;
}

double tail_p_dt(const Tail& tl, double t) {
  double acc = 0.0;
  for (std::size_t k = tl.c.size(); k-- > 1;) {
    acc = acc * t + static_cast<double>(k) * tl.c[k];
  }
  return acc;
}

double tail_p_dt2(const Tail& tl, double t) {
  double acc = 0.0;
  for (std::size_t k = tl.c.size(); k-- > 2;) {
    acc = acc * t + static_cast<double>(k * (k - 1)) * tl.c[k];
  }
  return acc;
}

// Integral of p over [0, t].
double tail_p_int(const Tail& tl, double t) {
  double acc = 0.0;
  for (std::size_t k = tl.c.size(); k-- > 0;) {
    acc = acc * t + tl.c[k] / static_cast<double>(k + 1);
  }
  return acc * t;
}

// Integral of Q over [0, t].
double tail_Q_int(const Tail& tl, double t) {
  double acc = 0.0;
  for (std::size_t k = tl.Q.size(); k-- > 0;) {
    acc = acc * t + tl.Q[k] / static_cast<double>(k + 1);
  }
  return acc * t;
}

double tail_H_R(const Tail& tl, double T) {
  const double t = (T - tl.T0) / tl.w;
  if (t >= 1.0) {
    return tl.H_end + tl.cp_inf * (T - tl.T_end);
  }
  return tl.H0 + tl.w * tail_p_int(tl, t);
}

// T > 0 (s_R holds T >= S_R_T_MIN), so the low tail, whose fade ends at 0 K,
// is on t < 1; t rounds to 1 below ~1e-14 T0, where the fade's own branch
// (Q_int(t) + r ln(T/T0), finite at t = 1) still applies, not the end value
// (-inf at T_end = 0 K).
double tail_S_R(const Tail& tl, double T) {
  const double t = (T - tl.T0) / tl.w;
  if (t >= 1.0 && tl.T_end > 0.0) {
    return tl.S_end + tl.cp_inf * std::log(T / tl.T_end);
  }
  return tl.S0 + tail_Q_int(tl, t) + tl.r * std::log(T / tl.T0);
}

// cp_end: the value p(1) the fade must reach (the quintic), or NaN to let the
// quartic level off where it does.
Tail make_tail(const std::array<double, 10>& a, double T0, double w, double cp_end) {
  const double T = T0;
  const double cp = a[0] / (T * T) + a[1] / T + a[2] + a[3] * T + a[4] * T * T +
                    a[5] * T * T * T + a[6] * T * T * T * T;
  const double d1 = -2.0 * a[0] / (T * T * T) - a[1] / (T * T) + a[3] + 2.0 * a[4] * T +
                    3.0 * a[5] * T * T + 4.0 * a[6] * T * T * T;
  const double d2 = 6.0 * a[0] / (T * T * T * T) + 2.0 * a[1] / (T * T * T) + 2.0 * a[4] +
                    6.0 * a[5] * T + 12.0 * a[6] * T * T;
  Tail tl;
  tl.T0 = T0;
  tl.w = w;
  const double c0 = cp;
  const double c1 = w * d1;
  const double c2 = 0.5 * w * w * d2;
  double c3 = 0.0;
  double c4 = 0.0;
  double c5 = 0.0;
  if (std::isnan(cp_end)) {
    // p'(1) = c1 + 2 c2 + 3 c3 + 4 c4 = 0 and p''(1) = 2 c2 + 6 c3 + 12 c4 = 0.
    c4 = 0.5 * (c1 + c2);
    c3 = -(3.0 * c1 + 4.0 * c2) / 3.0;
  } else {
    // And p(1) = cp_end: c3 + c4 + c5 = A, 3 c3 + 4 c4 + 5 c5 = B,
    // 6 c3 + 12 c4 + 20 c5 = C.
    const double A = cp_end - c0 - c1 - c2;
    const double B = -(c1 + 2.0 * c2);
    const double C = -2.0 * c2;
    c5 = 0.5 * (C + 12.0 * A - 6.0 * B);
    c4 = B - 3.0 * A - 2.0 * c5;
    c3 = A - c4 - c5;
  }
  tl.c = {c0, c1, c2, c3, c4, c5};
  tl.cp_inf = tail_p(tl, 1.0);
  tl.H0 = nasa_H_R(a, T0);
  tl.S0 = nasa_S_R(a, T0);
  // Synthetic division of p by (t - t_pole).
  const double tp = -T0 / w;
  tl.Q[4] = c5;
  tl.Q[3] = c4 + tp * tl.Q[4];
  tl.Q[2] = c3 + tp * tl.Q[3];
  tl.Q[1] = c2 + tp * tl.Q[2];
  tl.Q[0] = c1 + tp * tl.Q[1];
  tl.r = c0 + tp * tl.Q[0];
  tl.T_end = T0 + w;
  tl.H_end = tl.H0 + w * tail_p_int(tl, 1.0);
  // At 0 K (the low tail) s has no end value; tail_S_R never needs it there.
  tl.S_end = tl.T_end > 0.0
                 ? tl.S0 + tail_Q_int(tl, 1.0) + tl.r * std::log(tl.T_end / T0)
                 : -std::numeric_limits<double>::infinity();
  return tl;
}

// cp/R at 0 K: see thermo::CP_R_RIGID_* (H2's rotation freezes in the fade).
double cp_R_rigid_rotor(std::size_t species_idx) {
  if (species_names[species_idx] == "H2") {
    return thermo::CP_R_RIGID_ATOM;
  }
  switch (transport_props[species_idx].geometry) {
    case MolecularGeometry::Atom:
      return thermo::CP_R_RIGID_ATOM;
    case MolecularGeometry::Linear:
      return thermo::CP_R_RIGID_LINEAR;
    case MolecularGeometry::Nonlinear:
      return thermo::CP_R_RIGID_NONLINEAR;
  }
  return thermo::CP_R_RIGID_NONLINEAR;
}

const std::vector<SpeciesTails>& species_tails() {
  static const std::vector<SpeciesTails> tails = [] {
    std::vector<SpeciesTails> out;
    out.reserve(nasa_coeffs.size());
    for (std::size_t i = 0; i < nasa_coeffs.size(); ++i) {
      const auto& first = nasa_coeffs[i].intervals.front();
      const auto& last = nasa_coeffs[i].intervals.back();
      SpeciesTails st;
      st.T_lo = first.T_min;
      st.T_hi = last.T_max;
      st.lo = make_tail(first.coeffs, first.T_min, -first.T_min, cp_R_rigid_rotor(i));
      st.hi = make_tail(last.coeffs, last.T_max, thermo::TAIL_FRACTION_HIGH * last.T_max,
                        std::numeric_limits<double>::quiet_NaN());
      out.push_back(st);
    }
    return out;
  }();
  return tails;
}

// The tail covering T for this species, or nullptr inside the data.
const Tail* active_tail(std::size_t species_idx, double T) {
  const SpeciesTails& st = species_tails()[species_idx];
  if (T < st.T_lo) {
    return &st.lo;
  }
  if (T > st.T_hi) {
    return &st.hi;
  }
  return nullptr;
}

// H/R [K] of one species at any T (the tails are linear in T beyond the fade).
double species_H_R(std::size_t species_idx, double T) {
  if (const Tail* tl = active_tail(species_idx, T)) {
    return tail_H_R(*tl, T);
  }
  if (const C2Bridge* b = active_bridge(species_idx, T)) {
    return bridge_H_R(species_idx, *b, T);
  }
  return nasa_H_R(get_interval(species_idx, T).coeffs, T);
}

// Newton on a function increasing in T (h, s, u: their slopes cp, cp/T, cv
// are positive at every T with the tails), safeguarded as rtsafe (Numerical
// Recipes): steps within a factor 4 until the root is bracketed, then
// bisection whenever Newton leaves the bracket or fails to halve the step
// before it -- which breaks the 2-cycle a step in h or s at an interval join
// can trap Newton in. Throws if it does not converge, rather than returning a
// wrong T.
template <class FDF>
double solve_increasing_T(FDF&& f_df, double T_guess, double tol,
                          std::size_t max_iter, const char* name) {
  double T = (std::isfinite(T_guess) && T_guess > 0.0) ? T_guess : 300.0;
  double lo = 0.0;
  double hi = 0.0;
  bool have_lo = false;
  bool have_hi = false;
  double dx_old = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < max_iter; ++i) {
    const auto [f, df] = f_df(T);
    if (f == 0.0) {
      return T;
    }
    if (f > 0.0) {
      hi = T;
      have_hi = true;
    } else {
      lo = T;
      have_lo = true;
    }
    const double T_newton = (df > 0.0 && std::isfinite(df)) ? T - f / df
                                                              : (f > 0.0 ? 0.25 * T : 4.0 * T);
    double T_new = T_newton;
    if (have_lo && have_hi) {
      if (!(T_new >= lo && T_new <= hi) || std::abs(T_new - T) > 0.5 * dx_old) {
        T_new = 0.5 * (lo + hi);
      }
    } else if (have_lo) {
      T_new = std::min(T_new, 4.0 * T);
    } else {
      T_new = std::max(T_new, 0.25 * T);
    }
    // Converged on a small NEWTON step: its result is accurate to far below
    // tol, so the inverse is smooth in its target (a bisection midpoint is
    // only within tol, which finite differences see), and a clamped step
    // shrinking toward 0 K is a target no T reaches. Or on a bracket
    // narrower than tol: the published NASA-9 h and s step slightly at the
    // interval joins (H/R up to 5e-3 K, S/R up to 4e-7 at 6000 K), so a
    // target on a join can have no exact root for Newton to land on.
    // Both take tol or, where it is finer than T resolves, a few ulps of T:
    // above T ~ tol / (4 eps) = 1.1e9 K at tol 1e-6 K a Newton step from one
    // side cannot get below tol, and a probe at such a T threw (#1232).
    const double tol_T = std::max(tol, 4.0 * std::numeric_limits<double>::epsilon() * T);
    if (have_lo && have_hi &&
        hi - lo <= std::max(tol, 4.0 * std::numeric_limits<double>::epsilon() * hi)) {
      return T_new;
    }
    if (T_new == T_newton && std::abs(T_new - T) < tol_T) {
      return T_new;
    }
    dx_old = std::abs(T_new - T);
    T = T_new;
  }
  throw std::runtime_error(std::string(name) + ": no convergence in " +
                           std::to_string(max_iter) + " iterations (last T = " +
                           std::to_string(T) + " K)");
}
}  // namespace

void set_thermo_model(ThermoModel model) {
  g_thermo_model.store(static_cast<int>(model), std::memory_order_relaxed);
}

ThermoModel thermo_model() {
  return static_cast<ThermoModel>(g_thermo_model.load(std::memory_order_relaxed));
}

namespace {
// The breakpoints of one model (cp_breakpoints): every species' data ends and
// its high tail's fade end, and the NASA-9 joins (NASA9) or the C2 bridges'
// window edges (NASA9_C2). Sorted, unique, > 0 K.
std::vector<double> model_breakpoints(bool c2) {
  std::vector<double> out;
  for (std::size_t i = 0; i < nasa_coeffs.size(); ++i) {
    const SpeciesTails& st = species_tails()[i];
    out.push_back(st.T_lo);
    out.push_back(st.T_hi);
    out.push_back(st.hi.T_end);
    if (c2) {
      for (const auto& b : c2_bridges[i]) {
        out.push_back(b.T_join - b.half_width);
        out.push_back(b.T_join + b.half_width);
      }
    } else {
      for (const auto& iv : nasa_coeffs[i].intervals) {
        out.push_back(iv.T_min);
        out.push_back(iv.T_max);
      }
    }
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  out.erase(std::remove_if(out.begin(), out.end(), [](double t) { return !(t > 0.0); }),
            out.end());
  return out;
}
}  // namespace

const std::vector<double>& cp_breakpoints() {
  static const std::vector<double> nasa9 = model_breakpoints(false);
  static const std::vector<double> c2 = model_breakpoints(true);
  return thermo_model() == ThermoModel::NASA9_C2 ? c2 : nasa9;
}

double cp_R(std::size_t species_idx, double T) {
  if (species_idx >= species_names.size()) {
    throw std::out_of_range("cp_R: species index out of bounds");
  }
  if (const Tail* tl = active_tail(species_idx, T)) {
    const double t = (T - tl->T0) / tl->w;
    return t >= 1.0 ? tl->cp_inf : tail_p(*tl, t);
  }
  if (const C2Bridge* b = active_bridge(species_idx, T)) {
    return bridge_poly(*b, (T - b->T_join) / b->half_width);
  }
  const auto& a = get_interval(species_idx, T).coeffs;
  // Cp/R = a1/T^2 + a2/T + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
  return a[0]/(T*T) + a[1]/T + a[2] + a[3]*T + a[4]*T*T + a[5]*T*T*T + a[6]*T*T*T*T;
}

double h_RT(std::size_t species_idx, double T) {
  if (species_idx >= species_names.size()) {
    throw std::out_of_range("h_RT: species index out of bounds");
  }
  T = std::max(T, thermo::T_FLOOR);
  return species_H_R(species_idx, T) / T;
}

double s_R(std::size_t species_idx, double T) {
  if (species_idx >= species_names.size()) {
    throw std::out_of_range("s_R: species index out of bounds");
  }
  // Exact at every T > 0 (the low tail runs to 0 K, s -> -inf as cp ln T);
  // held only below S_R_T_MIN, which no solved T reaches (#935).
  T = std::max(T, thermo::S_R_T_MIN);
  if (const Tail* tl = active_tail(species_idx, T)) {
    return tail_S_R(*tl, T);
  }
  if (const C2Bridge* b = active_bridge(species_idx, T)) {
    return bridge_S_R(species_idx, *b, T);
  }
  const auto& a = get_interval(species_idx, T).coeffs;
  // S/R = -a1/(2*T^2) - a2/T + a3*ln(T) + a4*T + a5*T^2/2 + a6*T^3/3 + a7*T^4/4 + a9
  return nasa_S_R(a, T);
}

// Both terms held at T_FLOOR, as they were: g/RT feeds the equilibrium
// constants, whose callers floor their temperature there.
double g_over_RT(std::size_t species_idx, double T) {
  T = std::max(T, thermo::T_FLOOR);
  return h_RT(species_idx, T) - s_R(species_idx, T);
}

// -------------------------------------------------------------
// Dimensional property evaluation
// -------------------------------------------------------------

double cp_species(std::size_t species_idx, double T) {
  return cp_R(species_idx, T) * thermo::R_GAS;
}

double h_species(std::size_t species_idx, double T) {
  if (species_idx >= species_names.size()) {
    throw std::out_of_range("h_species: species index out of bounds");
  }
  return species_H_R(species_idx, T) * thermo::R_GAS;
}

double s_species(std::size_t species_idx, double T) {
  return s_R(species_idx, T) * thermo::R_GAS;
}

// -------------------------------------------------------------
// Mixture property evaluation
// -------------------------------------------------------------

// The mixture sums skip absent species (X[i] == 0): their term is an exact
// +-0.0, so the result is bit-identical, and a mixture usually carries a few
// of the species (dry air 4 of 15; transport already skipped them). The Fanno
// channel kernel evaluates these ~10^4 times per call.
double cp(double T, const std::vector<double> &X) {
  double cp_mix = 0.0;
  for (std::size_t i = 0; i < X.size(); ++i) {
    if (X[i] != 0.0) {
      cp_mix += X[i] * cp_species(i, T);
    }
  }
  return cp_mix;
}

double h(double T, const std::vector<double> &X) {
  double h_mix = 0.0;
  for (std::size_t i = 0; i < X.size(); ++i) {
    if (X[i] != 0.0) {
      h_mix += X[i] * h_species(i, T);
    }
  }
  return h_mix;
}

double s(double T, const std::vector<double> &X, double P, double P_ref) {
  double s_mix = 0.0;
  for (std::size_t i = 0; i < X.size(); ++i) {
    if (X[i] > 1e-15) {
      s_mix += X[i] * (s_species(i, T) - thermo::R_GAS * std::log(X[i] * P / P_ref));
    }
  }
  return s_mix;
}

double cv(double T, const std::vector<double> &X) {
  return cp(T, X) - thermo::R_GAS;
}

double u(double T, const std::vector<double> &X) {
  return h(T, X) - thermo::R_GAS * T;
}

double density(double T, double P, const std::vector<double> &X) {
  return (P * mwmix(X)) / (thermo::R_GAS * 1000.0 * T);
}

double molar_volume(double T, double P) {
  return (thermo::R_GAS * T) / P;
}

double specific_gas_constant(const std::vector<double> &X) {
  return (thermo::R_GAS * 1000.0) / mwmix(X);
}

double isentropic_expansion_coefficient(double T, const std::vector<double> &X) {
  // cv = cp - R_GAS (cv above): one mixture sum, the same value twice.
  const double c = cp(T, X);
  return c / (c - thermo::R_GAS);
}

double speed_of_sound(double T, const std::vector<double> &X) {
  double gamma = isentropic_expansion_coefficient(T, X);
  double R_spec = specific_gas_constant(X);
  return std::sqrt(gamma * R_spec * T);
}

// -------------------------------------------------------------
// Mass-specific thermodynamic properties
// -------------------------------------------------------------

double cp_mass(double T, const std::vector<double> &X) {
  return cp(T, X) / (mwmix(X) / 1000.0);
}

double cv_mass(double T, const std::vector<double> &X) {
  return cv(T, X) / (mwmix(X) / 1000.0);
}

std::tuple<double, double, double> cp_secant_mass_and_jacobian(double T1, double T2,
                                                               const std::vector<double> &X) {
  const double d = T1 - T2;
  if (std::abs(d) <= 1e-6 * std::max(std::abs(T1), 1.0)) {
    // cp(T1) + cp'(T1) d / 2 + O(d^2): the partials are cp' / 2 each.
    const double half = 0.5 * dcp_dT(0.5 * (T1 + T2), X) / (mwmix(X) / 1000.0);
    return {cp_mass(0.5 * (T1 + T2), X), half, half};
  }
  const double c = (h_mass(T1, X) - h_mass(T2, X)) / d;
  return {c, (cp_mass(T1, X) - c) / d, (c - cp_mass(T2, X)) / d};
}

std::vector<double> cp_secant_mass_dY(double T1, double T2, const std::vector<double> &X) {
  const double d = T1 - T2;
  const bool limit = std::abs(d) <= 1e-6 * std::max(std::abs(T1), 1.0);
  const double c = std::get<0>(cp_secant_mass_and_jacobian(T1, T2, X));
  const double MW = mwmix(X) / 1000.0;  // kg/mol
  std::vector<double> g(X.size(), 0.0);
  for (std::size_t j = 0; j < X.size(); ++j) {
    // The species' own molar secant; the mixture's is sum X_j c_j / MW.
    const double cj = limit ? cp_species(j, 0.5 * (T1 + T2))
                            : (h_species(j, T1) - h_species(j, T2)) / d;
    g[j] = (cj - c * species_molar_mass(j) / 1000.0) / MW;
  }
  return mole_gradient_to_mass(X, g);
}

std::vector<double> specific_gas_constant_dY(const std::vector<double> &X) {
  // R = R_GAS / MW: dR/dX_j = -R W_j / MW.
  const double MW = mwmix(X);
  const double R = specific_gas_constant(X);
  std::vector<double> g(X.size(), 0.0);
  for (std::size_t j = 0; j < X.size(); ++j) {
    g[j] = -R * species_molar_mass(j) / MW;
  }
  return mole_gradient_to_mass(X, g);
}

std::pair<double, std::vector<double>> isentropic_expansion_coefficient_dT_dY(
    double T, const std::vector<double> &X) {
  // gamma = c / (c - R_GAS) with c the molar cp: dgamma/dc = -R_GAS / (c -
  // R_GAS)^2, and c moves with T through dcp_dT and with X_j by cp_j.
  const double c = cp(T, X);
  const double k = -thermo::R_GAS / ((c - thermo::R_GAS) * (c - thermo::R_GAS));
  std::vector<double> g(X.size(), 0.0);
  for (std::size_t j = 0; j < X.size(); ++j) {
    g[j] = k * cp_species(j, T);
  }
  return {k * dcp_dT(T, X), mole_gradient_to_mass(X, g)};
}

double h_mass(double T, const std::vector<double> &X) {
  return h(T, X) / (mwmix(X) / 1000.0);
}

double s_mass(double T, const std::vector<double> &X, double P, double P_ref) {
  return s(T, X, P, P_ref) / (mwmix(X) / 1000.0);
}

double u_mass(double T, const std::vector<double> &X) {
  return u(T, X) / (mwmix(X) / 1000.0);
}

// -------------------------------------------------------------
// Air Properties Bundle
// -------------------------------------------------------------

AirProperties air_properties(double T, double P, double humidity) {
  if (!validation::is_finite_positive(T)) {
    throw std::invalid_argument("air_properties: temperature must be positive");
  }
  if (!validation::is_finite_positive(P)) {
    throw std::invalid_argument("air_properties: pressure must be positive");
  }
  if (!validation::is_finite_nonnegative(humidity) || humidity > 1.0) {
    throw std::invalid_argument(
        "air_properties: humidity must be in range [0, 1]");
  }

  std::vector<double> X;
  if (humidity > 0.0) {
    X = humid_air_composition(T, P, humidity);
  } else {
    X = dry_air();
  }
  return transport_state(T, P, X);
}

// -------------------------------------------------------------
// Thermodynamic State Bundle
// -------------------------------------------------------------

ThermoState thermo_state(double T, double P, const std::vector<double> &X, double P_ref) {
  if (!validation::is_finite_positive(T)) throw std::invalid_argument("thermo_state: temperature must be positive");
  if (!validation::is_finite_positive(P)) throw std::invalid_argument("thermo_state: pressure must be positive");
  if (X.empty()) throw std::invalid_argument("thermo_state: composition must not be empty");
  if (!validation::is_finite_positive(P_ref)) throw std::invalid_argument("thermo_state: reference pressure must be positive");
  ThermoState ts;
  ts.T = T;
  ts.P = P;
  ts.mw = mwmix(X);
  double mw_kg = ts.mw / 1000.0;
  ts.cp_mole = cp(T, X);
  ts.h_mole = h(T, X);
  ts.s_mole = s(T, X, P, P_ref);
  ts.cv_mole = cv(T, X);
  ts.u_mole = u(T, X);
  ts.rho = density(T, P, X);
  ts.gamma = ts.cp_mole / ts.cv_mole;
  ts.a = std::sqrt(ts.gamma * (thermo::R_GAS / mw_kg) * T);
  ts.cp = ts.cp_mole / mw_kg;
  ts.cv = ts.cv_mole / mw_kg;
  ts.h = ts.h_mole / mw_kg;
  ts.s = ts.s_mole / mw_kg;
  ts.u = ts.u_mole / mw_kg;
  return ts;
}

// -------------------------------------------------------------
// Complete State Bundle
// -------------------------------------------------------------

CompleteState complete_state(double T, double P, const std::vector<double> &X, double P_ref) {
  CompleteState cs;
  cs.thermo = thermo_state(T, P, X, P_ref);
  cs.transport = transport_state(T, P, X);
  return cs;
}

// -------------------------------------------------------------
// Derivatives wrt temperature
// -------------------------------------------------------------

double dh_dT(double T, const std::vector<double> &X) { return cp(T, X); }

double ds_dT(double T, const std::vector<double> &X) { return cp(T, X) / T; }

double dcp_R_dT(std::size_t species_idx, double T) {
  if (species_idx >= species_names.size()) {
    throw std::out_of_range("dcp_R_dT: species index out of bounds");
  }
  if (const Tail* tl = active_tail(species_idx, T)) {
    const double t = (T - tl->T0) / tl->w;
    return t >= 1.0 ? 0.0 : tail_p_dt(*tl, t) / tl->w;
  }
  if (const C2Bridge* b = active_bridge(species_idx, T)) {
    return bridge_poly_dt(*b, (T - b->T_join) / b->half_width) / b->half_width;
  }
  const auto& a = get_interval(species_idx, T).coeffs;
  // d(Cp/R)/dT = -2*a1/T^3 - a2/T^2 + a4 + 2*a5*T + 3*a6*T^2 + 4*a7*T^3
  return -2.0*a[0]/(T*T*T) - a[1]/(T*T) + a[3] + 2.0*a[4]*T +
         3.0*a[5]*T*T + 4.0*a[6]*T*T*T;
}

double dcp_dT(double T, const std::vector<double> &X) {
  double dcp_mix = 0.0;
  for (std::size_t i = 0; i < X.size(); ++i) {
    if (X[i] != 0.0) {
      dcp_mix += X[i] * dcp_R_dT(i, T) * thermo::R_GAS;
    }
  }
  return dcp_mix;
}

double d2cp_R_dT2(std::size_t species_idx, double T) {
  if (species_idx >= species_names.size()) {
    throw std::out_of_range("d2cp_R_dT2: species index out of bounds");
  }
  if (const Tail* tl = active_tail(species_idx, T)) {
    const double t = (T - tl->T0) / tl->w;
    return t >= 1.0 ? 0.0 : tail_p_dt2(*tl, t) / (tl->w * tl->w);
  }
  if (const C2Bridge* b = active_bridge(species_idx, T)) {
    return bridge_poly_dt2(*b, (T - b->T_join) / b->half_width) / (b->half_width * b->half_width);
  }
  const auto& a = get_interval(species_idx, T).coeffs;
  // d^2(Cp/R)/dT^2 = 6*a1/T^4 + 2*a2/T^3 + 2*a5 + 6*a6*T + 12*a7*T^2
  const double T2 = T * T;
  return 6.0 * a[0] / (T2 * T2) + 2.0 * a[1] / (T2 * T) + 2.0 * a[4] + 6.0 * a[5] * T +
         12.0 * a[6] * T2;
}

double d2cp_dT2(double T, const std::vector<double> &X) {
  double d2cp_mix = 0.0;
  for (std::size_t i = 0; i < X.size(); ++i) {
    if (X[i] != 0.0) {
      d2cp_mix += X[i] * d2cp_R_dT2(i, T) * thermo::R_GAS;
    }
  }
  return d2cp_mix;
}

double dg_over_RT_dT(double T, const std::vector<double> &X) {
  return -h(T, X) / (thermo::R_GAS * T * T);
}

// -------------------------------------------------------------
// Inverse solvers (T from properties)
// -------------------------------------------------------------

double calc_T_from_h(double h_target, const std::vector<double> &X, double T_guess, double tol, std::size_t max_iter) {
  return solve_increasing_T(
      [&](double T) { return std::make_pair(h(T, X) - h_target, cp(T, X)); }, T_guess,
      tol, max_iter, "calc_T_from_h");
}

double calc_T_from_s(double s_target, double P, const std::vector<double> &X, double T_guess, double tol, std::size_t max_iter) {
  return solve_increasing_T(
      [&](double T) { return std::make_pair(s(T, X, P) - s_target, cp(T, X) / T); },
      T_guess, tol, max_iter, "calc_T_from_s");
}

namespace {

// calc_T_from_cp's scan nodes (thermo.h): CP_INVERSE_SCAN_STEP K apart down
// to CP_INVERSE_LOG_SCAN_T, the same share of T below it; downward they step
// onto CP_INVERSE_T_DATA and CP_INVERSE_LOG_SCAN_T, and both ends clamp.
constexpr double kCpScanRatio = thermo::CP_INVERSE_SCAN_STEP / thermo::CP_INVERSE_LOG_SCAN_T;
// Relative slack in the nearer-node comparison: T0 + 10 - T0 and T0 - (T0 - 10)
// differ in the last bits, and an equal-step tie must go up.
constexpr double kCpScanTieRel = 1.0e-12;
static_assert(thermo::CP_INVERSE_T_MIN < thermo::CP_INVERSE_LOG_SCAN_T &&
                  thermo::CP_INVERSE_LOG_SCAN_T < thermo::CP_INVERSE_T_DATA &&
                  thermo::CP_INVERSE_T_DATA < thermo::CP_INVERSE_T_MAX && kCpScanRatio < 1.0,
              "calc_T_from_cp scan nodes out of order");

double cp_scan_up(double T) {
  if (T < thermo::CP_INVERSE_LOG_SCAN_T) {
    return std::min(T / (1.0 - kCpScanRatio), thermo::CP_INVERSE_LOG_SCAN_T);
  }
  return std::min(T + thermo::CP_INVERSE_SCAN_STEP, thermo::CP_INVERSE_T_MAX);
}

double cp_scan_down(double T) {
  if (T > thermo::CP_INVERSE_T_DATA) {
    return std::max(T - thermo::CP_INVERSE_SCAN_STEP, thermo::CP_INVERSE_T_DATA);
  }
  if (T > thermo::CP_INVERSE_LOG_SCAN_T) {
    return std::max(T - thermo::CP_INVERSE_SCAN_STEP, thermo::CP_INVERSE_LOG_SCAN_T);
  }
  return std::max(T * (1.0 - kCpScanRatio), thermo::CP_INVERSE_T_MIN);
}

// Newton kept inside the bracket [a, b] of cp - cp_target (fa, fb of opposite
// sign; bisection when a step leaves it); converged on a Newton step below
// tol or a bracket narrower than tol.
double polish_cp_root(double cp_target, const std::vector<double> &X, double a, double fa,
                      double b, double fb, double tol, std::size_t max_iter, const char *fn) {
  auto f = [&](double T) { return cp(T, X) - cp_target; };
  double T = a - fa * (b - a) / (fb - fa);
  for (std::size_t i = 0; i < max_iter; ++i) {
    const double fT = f(T);
    if (fT == 0.0) {
      return T;
    }
    if ((fT > 0.0) == (fa > 0.0)) {
      a = T, fa = fT;
    } else {
      b = T;
    }
    if (b - a <= tol) {
      return 0.5 * (a + b);
    }
    const double d = dcp_dT(T, X);
    const double T_newton = (d != 0.0 && std::isfinite(d)) ? T - fT / d : a - 1.0;
    if (T_newton > a && T_newton < b) {
      const double step = std::abs(T_newton - T);
      T = T_newton;
      if (step < tol) {
        return T;
      }
    } else {
      T = 0.5 * (a + b);
    }
  }
  inner_solve::throw_no_convergence(fn, max_iter, "cp residual [J/(mol K)]", f(T));
}

}  // namespace

double calc_T_from_cp(double cp_target, const std::vector<double> &X,
                      double T_guess, double tol, std::size_t max_iter) {
  // cp(T) need not be monotone, so the root is a branch choice: the one
  // nearest T_guess inside [CP_INVERSE_T_MIN, CP_INVERSE_T_MAX] (thermo.h),
  // found by an outward scan and polished on its bracket (#868 R6).
  if (!validation::is_finite_positive(cp_target)) {
    throw std::invalid_argument("calc_T_from_cp: cp_target must be positive and finite");
  }
  if (!validation::is_finite_positive(T_guess)) {
    throw std::invalid_argument("calc_T_from_cp: T_guess must be positive and finite");
  }
  if (!validation::is_finite_positive(tol)) {
    throw std::invalid_argument("calc_T_from_cp: tol must be positive and finite");
  }
  constexpr double T_min = thermo::CP_INVERSE_T_MIN;
  constexpr double T_max = thermo::CP_INVERSE_T_MAX;
  auto f = [&](double T) { return cp(T, X) - cp_target; };
  const double T0 = std::clamp(T_guess, T_min, T_max);
  const double f0 = f(T0);
  if (f0 == 0.0) {
    return T0;
  }
  // Outward scan: each step advances the side whose next node is nearer T0
  // (a tie goes up; within rounding, so that equal steps alternate up, down,
  // as the [200, 6000] K scan did). The first segment whose ends differ in
  // sign brackets the root nearest T0, to within the local step.
  double a = 0.0, b = 0.0, fa = 0.0, fb = 0.0;
  bool bracketed = false;
  double lo = T0, f_lo = f0, hi = T0, f_hi = f0;
  while (lo > T_min || hi < T_max) {
    const double up = hi < T_max ? cp_scan_up(hi) : T_max;
    const double down = lo > T_min ? cp_scan_down(lo) : T_min;
    const bool go_up =
        hi < T_max && (lo <= T_min || up - T0 <= (T0 - down) * (1.0 + kCpScanTieRel));
    if (go_up) {
      const double f_next = f(up);
      if (f_next == 0.0) {
        return up;
      }
      if ((f_next > 0.0) != (f_hi > 0.0)) {
        a = hi, fa = f_hi, b = up, fb = f_next;
        bracketed = true;
        break;
      }
      hi = up, f_hi = f_next;
    } else {
      const double f_next = f(down);
      if (f_next == 0.0) {
        return down;
      }
      if ((f_next > 0.0) != (f_lo > 0.0)) {
        a = down, fa = f_next, b = lo, fb = f_lo;
        bracketed = true;
        break;
      }
      lo = down, f_lo = f_next;
    }
  }
  if (!bracketed) {
    char buf[256];
    std::snprintf(buf, sizeof buf,
                  "calc_T_from_cp: no root in [%g, %g] K (cp_target %.6g J/(mol K), cp "
                  "%.6g at %g K, %.6g at %g K)",
                  T_min, T_max, cp_target, f_lo + cp_target, lo, f_hi + cp_target, hi);
    throw std::runtime_error(buf);
  }
  return polish_cp_root(cp_target, X, a, fa, b, fb, tol, max_iter, "calc_T_from_cp");
}

std::vector<CpRoot> calc_T_from_cp_roots(double cp_target, const std::vector<double> &X,
                                         double T_min, double T_max, double tol,
                                         std::size_t max_iter) {
  // Every root in [T_min, T_max] (#963). cp is split at its extrema into
  // monotone pieces, which the scan nodes subdivide; each sign change of
  // cp - cp_target between neighbouring points is one root.
  if (!validation::is_finite_positive(cp_target)) {
    throw std::invalid_argument("calc_T_from_cp_roots: cp_target must be positive and finite");
  }
  if (!validation::is_finite_positive(tol)) {
    throw std::invalid_argument("calc_T_from_cp_roots: tol must be positive and finite");
  }
  if (!std::isfinite(T_min) || !std::isfinite(T_max) || T_min < thermo::CP_INVERSE_T_MIN ||
      T_max > thermo::CP_INVERSE_T_MAX || !(T_min < T_max)) {
    char buf[160];
    std::snprintf(buf, sizeof buf,
                  "calc_T_from_cp_roots: need %g <= T_min < T_max <= %g K (got %g, %g)",
                  thermo::CP_INVERSE_T_MIN, thermo::CP_INVERSE_T_MAX, T_min, T_max);
    throw std::invalid_argument(buf);
  }
  auto sign = [](double v) { return (v > 0.0) - (v < 0.0); };
  auto dcp = [&](double T) { return dcp_dT(T, X); };
  // Bisect dcp/dT to adjacent doubles between a (sign sa) and b (not sa); the
  // more extreme end of cp: a maximum when cp was rising (sa > 0).
  auto extremum = [&](double a, int sa, double b) {
    for (double m = 0.5 * (a + b); m > a && m < b; m = 0.5 * (a + b)) {
      if (sign(dcp(m)) == sa) {
        a = m;
      } else {
        b = m;
      }
    }
    return (cp(a, X) > cp(b, X)) == (sa > 0) ? a : b;
  };
  // Scan nodes from T_min up, and dcp/dT on them.
  std::vector<double> pts;
  std::vector<double> ds;
  for (double T = T_min;; T = std::min(cp_scan_up(T), T_max)) {
    pts.push_back(T);
    ds.push_back(dcp(T));
    if (T >= T_max) {
      break;
    }
  }
  const std::size_t n = pts.size();
  std::vector<double> ext;
  // Sign changes of dcp/dT between nodes (a node where it is exactly 0, as on
  // argon's constant cp, is passed over): one extremum each.
  std::size_t k_last = n;
  for (std::size_t k = 0; k < n; ++k) {
    const int s = sign(ds[k]);
    if (s == 0) {
      continue;
    }
    if (k_last < n && s != sign(ds[k_last])) {
      ext.push_back(extremum(pts[k_last], sign(ds[k_last]), pts[k]));
    }
    k_last = k;
  }
  // An extremum pair inside one step leaves no sign change on the nodes, but
  // |dcp/dT| dips there: at a node where it is a local minimum (the range
  // ends included) with one sign around it, search the window of the
  // neighbouring nodes for the opposite sign (golden section on sign*dcp/dT)
  // and, when found, bisect on both sides of it.
  constexpr double kInvPhi = 0.6180339887498949;  // (sqrt(5) - 1) / 2
  for (std::size_t k = 0; k < n; ++k) {
    const std::size_t lo = k > 0 ? k - 1 : k;
    const std::size_t hi = k + 1 < n ? k + 1 : k;
    const int s = sign(ds[k]);
    if (s == 0 || sign(ds[lo]) != s || sign(ds[hi]) != s || lo == hi ||
        !(k == lo || std::abs(ds[k]) < std::abs(ds[lo])) ||
        !(k == hi || std::abs(ds[k]) <= std::abs(ds[hi]))) {
      continue;
    }
    auto g = [&](double T) { return s * dcp(T); };
    double a = pts[lo], b = pts[hi];
    double c = b - kInvPhi * (b - a), d = a + kInvPhi * (b - a);
    double gc = g(c), gd = g(d);
    double m = std::numeric_limits<double>::quiet_NaN();
    for (int i = 0; i < thermo::CP_ROOTS_DIP_ITER; ++i) {
      if (gc < 0.0 || gd < 0.0) {
        m = gc < 0.0 ? c : d;
        break;
      }
      if (gc < gd) {
        b = d, d = c, gd = gc;
        c = b - kInvPhi * (b - a);
        gc = g(c);
      } else {
        a = c, c = d, gc = gd;
        d = a + kInvPhi * (b - a);
        gd = g(d);
      }
    }
    if (!std::isnan(m)) {
      ext.push_back(extremum(pts[lo], s, m));
      ext.push_back(extremum(m, -s, pts[hi]));
    }
  }
  std::sort(ext.begin(), ext.end());
  // Merge the extrema in (an extremum on a node replaces it).
  std::vector<bool> is_ext(pts.size(), false);
  for (double e : ext) {
    const auto it = std::lower_bound(pts.begin(), pts.end(), e);
    const auto k = static_cast<std::size_t>(it - pts.begin());
    if (it != pts.end() && *it == e) {
      is_ext[k] = true;
    } else {
      pts.insert(it, e);
      is_ext.insert(is_ext.begin() + static_cast<std::ptrdiff_t>(k), true);
    }
  }
  std::vector<double> fs(pts.size());
  for (std::size_t k = 0; k < pts.size(); ++k) {
    fs[k] = cp(pts[k], X) - cp_target;
    if (is_ext[k] && std::abs(fs[k]) <= thermo::CP_ROOTS_DOUBLE_REL * cp_target) {
      fs[k] = 0.0;
    }
  }
  std::vector<CpRoot> roots;
  for (std::size_t k = 0; k < pts.size(); ++k) {
    if (k > 0 && fs[k - 1] != 0.0 && fs[k] != 0.0 && sign(fs[k - 1]) != sign(fs[k])) {
      const double T = polish_cp_root(cp_target, X, pts[k - 1], fs[k - 1], pts[k], fs[k], tol,
                                      max_iter, "calc_T_from_cp_roots");
      roots.push_back({T, 1.0 / dcp_dT(T, X), false});
    }
    if (fs[k] == 0.0) {
      if (k > 0 && fs[k - 1] == 0.0) {
        std::size_t lo = k - 1, hi = k;
        while (lo > 0 && fs[lo - 1] == 0.0) {
          --lo;
        }
        while (hi + 1 < pts.size() && fs[hi + 1] == 0.0) {
          ++hi;
        }
        char buf[200];
        std::snprintf(buf, sizeof buf,
                      "calc_T_from_cp_roots: cp equals cp_target %.6g J/(mol K) at every scan "
                      "node in [%g, %g] K; the roots are not isolated",
                      cp_target, pts[lo], pts[hi]);
        throw std::domain_error(buf);
      }
      roots.push_back({pts[k],
                       is_ext[k] ? std::numeric_limits<double>::quiet_NaN()
                                 : 1.0 / dcp_dT(pts[k], X),
                       is_ext[k]});
    }
  }
  return roots;
}

double calc_T_from_u(double u_target, const std::vector<double> &X, double T_guess, double tol, std::size_t max_iter) {
  return solve_increasing_T(
      [&](double T) { return std::make_pair(u(T, X) - u_target, cv(T, X)); }, T_guess,
      tol, max_iter, "calc_T_from_u");
}

double calc_T_from_h_mass(double h_mass_target, const std::vector<double> &X, double T_guess, double tol, std::size_t max_iter) {
  return calc_T_from_h(h_mass_target * (mwmix(X) / 1000.0), X, T_guess, tol, max_iter);
}

double calc_T_from_s_mass(double s_mass_target, double P, const std::vector<double> &X, double T_guess, double tol, std::size_t max_iter) {
  return calc_T_from_s(s_mass_target * (mwmix(X) / 1000.0), P, X, T_guess, tol, max_iter);
}

double calc_T_from_u_mass(double u_mass_target, const std::vector<double> &X, double T_guess, double tol, std::size_t max_iter) {
  return calc_T_from_u(u_mass_target * (mwmix(X) / 1000.0), X, T_guess, tol, max_iter);
}

double calc_T_from_sv_mass(double s_mass_target, double v_mass_target, const std::vector<double>& X, double T_guess, double tol, std::size_t max_iter) {
  const double mw_kg = mwmix(X) / 1000.0;
  const double rho = 1.0 / v_mass_target;
  // ds/dT at fixed v = cp/T - R/P dP/dT = cv/T.
  return solve_increasing_T(
      [&](double T) {
        const double P = (rho * thermo::R_GAS * T) / mw_kg;
        return std::make_pair(s_mass(T, X, P) - s_mass_target, cv_mass(T, X) / T);
      },
      T_guess, tol, max_iter, "calc_T_from_sv_mass");
}

double calc_T_from_sh_mass([[maybe_unused]] double s_mass_target, double h_mass_target, const std::vector<double>& X, double T_guess, double tol, std::size_t max_iter) {
  // h = h(T) only for ideal gas. Find T from h first.
  double T = calc_T_from_h_mass(h_mass_target, X, T_guess, tol, max_iter);
  // Then P is determined by s(T, P) = s_target
  return T;
}

// -------------------------------------------------------------
// State-based overloads
// -------------------------------------------------------------

double mwmix(const State &s) { return combaero::mwmix(s.X); }
double cp(const State &s) { return combaero::cp(s.T, s.X); }
double h(const State &s) { return combaero::h(s.T, s.X); }
double s(const State &s) { return combaero::s(s.T, s.X, s.P); }
double cv(const State &s) { return combaero::cv(s.T, s.X); }
double u(const State &s) { return combaero::u(s.T, s.X); }
double density(const State &s) { return combaero::density(s.T, s.P, s.X); }
double specific_gas_constant(const State &s) { return combaero::specific_gas_constant(s.X); }
double isentropic_expansion_coefficient(const State &s) { return combaero::isentropic_expansion_coefficient(s.T, s.X); }
double speed_of_sound(const State &s) { return combaero::speed_of_sound(s.T, s.X); }

} // namespace combaero
