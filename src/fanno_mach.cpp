#include "../include/fanno_mach.h"

#include "../include/composition.h"
#include "../include/dual_number.h"
#include "../include/friction.h"
#include "../include/inner_solve.h"
#include "../include/input_validation.h"
#include "../include/solver_interface.h"
#include "../include/stagnation.h"
#include "../include/thermo.h"
#include "../include/transport.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace combaero {

namespace {

double sound_speed_sq(double T, const std::vector<double>& X) {
  return isentropic_expansion_coefficient(T, X) * specific_gas_constant(X) * T;
}

// d(a^2)/dT analytically: a^2 = gamma R T, gamma = cp/cv with cv = cp - R_u
// (molar), so dgamma/dT = -R_u cp' / cv^2. Exact, so the length integrand
// carries no finite-difference noise and the quadrature can meet its
// tolerance.
double sound_speed_sq_dT(double T, const std::vector<double>& X) {
  const double cp_mol = cp(T, X);
  const double cv_mol = cp_mol - thermo::R_GAS;
  const double gamma = cp_mol / cv_mol;
  const double dgamma = -thermo::R_GAS * dcp_dT(T, X) / (cv_mol * cv_mol);
  return specific_gas_constant(X) * (gamma + T * dgamma);
}

// d^2(a^2)/dT^2 analytically, for the Rayleigh RK4 discrete tangent (#1276).
double sound_speed_sq_d2T(double T, const std::vector<double>& X) {
  const double cp_mol = cp(T, X);
  const double cv_mol = cp_mol - thermo::R_GAS;
  const double dcp = dcp_dT(T, X);
  const double d2cp = d2cp_dT2(T, X);
  const double inv_cv2 = 1.0 / (cv_mol * cv_mol);
  const double dgamma = -thermo::R_GAS * dcp * inv_cv2;
  const double d2gamma = -thermo::R_GAS * (d2cp * cv_mol - 2.0 * dcp * dcp) * (inv_cv2 / cv_mol);
  return specific_gas_constant(X) * (2.0 * dgamma + T * d2gamma);
}

// Static temperature at Mach M from the stagnation temperature:
// h(T) + M^2 a(T)^2 / 2 = h(Tt). F(T) rises with T, F(Tt) >= 0, so a
// safeguarded Newton on a bracket [T_lo, Tt] with F(T_lo) < 0 cannot leave
// the physical branch.
double static_T(double Tt, double M, const std::vector<double>& X, double T_guess = -1.0) {
  const double h0 = h_mass(Tt, X);
  const double R = specific_gas_constant(X);
  auto F = [&](double T) {
    const double a2 = isentropic_expansion_coefficient(T, X) * R * T;
    return h_mass(T, X) + 0.5 * M * M * a2 - h0;
  };
  if (M <= 0.0) {
    return Tt;
  }
  // The bracket must hold the root: F(lo) < 0 <= F(hi). The usual floor,
  // max(20 K, 0.05 Tt), does not wherever the static T lies below it: a
  // stagnation T of 20 K or less (the floor at or above Tt), Tt up to ~24 K
  // near M 1, or M > 1 at any Tt (M 3 at 50 K). There the safeguarded Newton
  // returned the floor itself as the root, with no flag (#1305: a duct cooled
  // 450 -> 9.7 K gave a static T above Tt and a warm-started solve a
  // non-root). The ideal gas holds to 0 K (#626) and F(0+) = h(0) - h(Tt) < 0,
  // so lowering the floor always finds the sign change; where the usual floor
  // already holds it, nothing changes.
  constexpr double kFloorT = 20.0;       // [K]
  constexpr double kFloorFraction = 0.05;  // of Tt; also the lowering factor
  constexpr int kMaxLower = 20;          // 0.05^20 Tt ~ 1e-26 Tt
  double lo = std::max(kFloorT, kFloorFraction * Tt);
  double hi = Tt;
  for (int k = 0; !(F(lo) < 0.0); ++k) {
    if (k == kMaxLower) {
      inner_solve::throw_no_convergence("fanno static_T bracket", kMaxLower, "h residual [J/kg]",
                                        F(lo));
    }
    hi = std::min(hi, lo);
    lo *= kFloorFraction;
  }
  // Ideal-gas guess at constant gamma (or the caller's, from a neighbouring
  // station), then Newton with the bracket. The root does not depend on the
  // start, to the 1e-14 the loop converges to.
  double T = T_guess;
  if (!(T > lo && T < hi)) {
    const double g = isentropic_expansion_coefficient(Tt, X);
    T = Tt / (1.0 + 0.5 * (g - 1.0) * M * M);
  }
  constexpr int kMaxIter = 100;
  for (int it = 0; it < kMaxIter; ++it) {
    const double f = F(T);
    if (f > 0.0) {
      hi = T;
    } else {
      lo = T;
    }
    const double fp = cp_mass(T, X) + 0.5 * M * M * sound_speed_sq_dT(T, X);
    double T_new = (fp > 0.0) ? T - f / fp : 0.5 * (lo + hi);
    if (!(T_new > lo && T_new < hi)) {
      T_new = 0.5 * (lo + hi);
    }
    if (std::abs(T_new - T) <= 1e-14 * T) {
      return T_new;
    }
    T = T_new;
  }
  // A cap exit is a root only if the bracket has closed onto it (#868).
  if (inner_solve::cap_exit_is_root(hi - lo, T)) {
    return T;
  }
  inner_solve::throw_no_convergence("fanno static_T", kMaxIter, "h residual [J/kg]", F(T));
}

double friction(double Re, double e_D, const std::string& model, double f_mult) {
  if (model == "fixed") {
    return f_mult;
  }
  return std::get<0>(solver::friction_and_jacobian(model, Re, e_D).result) * f_mult;
}

// Composite 12-point Gauss-Legendre on kFannoPanels geometric panels. A FIXED
// rule, not an adaptive one: the length is then a smooth function of G and the
// Mach limits, so differences of it (the network Jacobians) carry no noise from
// a refinement pattern that switches as the arguments move. Measured against a
// 64-panel reference over inlet Mach 0.01-0.99: <= 2e-10 at constant f; with a
// Re-dependent f the floor is ~1e-8, set by kinks in mu(T), not by the rule.
constexpr int kFannoPanels = 4;
constexpr std::array<double, 6> kGlX = {0.12523340851146894, 0.3678314989981802,
                                        0.5873179542866175,  0.7699026741943047,
                                        0.9041172563704748,  0.9815606342467192};
constexpr std::array<double, 6> kGlW = {0.24914704581340297, 0.2334925365383549,
                                        0.2031674267230659,  0.16007832854334636,
                                        0.10693932599531876, 0.04717533638651128};

template <typename Fn>
double composite_gl(const Fn& f, double a, double b) {
  // Geometric panel edges: the interval in s = 1/M^2 spans orders of magnitude
  // at low Mach.
  const double ratio = std::pow(b / a, 1.0 / kFannoPanels);
  double total = 0.0;
  double lo = a;
  for (int k = 0; k < kFannoPanels; ++k) {
    const double hi = (k == kFannoPanels - 1) ? b : lo * ratio;
    const double c = 0.5 * (lo + hi);
    const double h = 0.5 * (hi - lo);
    double acc = 0.0;
    for (std::size_t j = 0; j < kGlX.size(); ++j) {
      acc += kGlW[j] * (f(c - h * kGlX[j]) + f(c + h * kGlX[j]));
    }
    total += acc * h;
    lo = hi;
  }
  return total;
}

// Brent's method for a sign change of fn on [a, b]. Without one there is no
// root to converge to and it throws: the bracket searches before each call
// stop at a budget (60 e-folds of flux, for example), and an interval they
// leave unbracketed used to return a point of it as the root (#868).
template <typename Fn>
double brent(const Fn& fn, double a, double b, double fa, double fb, double rel_tol) {
  if (fa != 0.0 && fb != 0.0 && (fa > 0.0) == (fb > 0.0)) {
    char buf[160];
    std::snprintf(buf, sizeof buf,
                  "fanno (Mach) brent: no sign change on [%.6e, %.6e] (f = %.6e, %.6e)", a, b,
                  fa, fb);
    throw std::runtime_error(buf);
  }
  double c = a, fc = fa, d = b - a, e = d;
  for (int it = 0; it < 200; ++it) {
    if ((fb > 0.0) == (fc > 0.0)) {
      c = a;
      fc = fa;
      d = e = b - a;
    }
    if (std::abs(fc) < std::abs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    const double tol = 2.0 * 1e-16 * std::abs(b) + 0.5 * rel_tol * std::abs(b);
    const double m = 0.5 * (c - b);
    if (std::abs(m) <= tol || fb == 0.0) {
      return b;
    }
    if (std::abs(e) >= tol && std::abs(fa) > std::abs(fb)) {
      double p, q, r;
      const double s = fb / fa;
      if (a == c) {
        p = 2.0 * m * s;
        q = 1.0 - s;
      } else {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0) {
        q = -q;
      } else {
        p = -p;
      }
      if (2.0 * p < std::min(3.0 * m * q - std::abs(tol * q), std::abs(e * q))) {
        e = d;
        d = p / q;
      } else {
        d = m;
        e = m;
      }
    } else {
      d = m;
      e = m;
    }
    a = b;
    fa = fb;
    b += (std::abs(d) > tol) ? d : (m > 0.0 ? tol : -tol);
    fb = fn(b);
  }
  return b;
}

// Isentropic static pressure at T from the stagnation state: for an ideal gas
// at fixed composition s(T, P) = s(T, Pref) - R ln(P / Pref), so equal entropy
// gives P in closed form.
double isentropic_P(double Pt, double Tt, double T, const std::vector<double>& X) {
  const double R = specific_gas_constant(X);
  return Pt * std::exp((s_mass(T, X, 101325.0) - s_mass(Tt, X, 101325.0)) / R);
}

double isentropic_flux(double Pt, double Tt, double M, const std::vector<double>& X) {
  if (M <= 0.0) {
    return 0.0;
  }
  const double T = static_T(Tt, M, X);
  const double P = isentropic_P(Pt, Tt, T, X);
  return density(T, P, X) * M * std::sqrt(sound_speed_sq(T, X));
}

// d flux/dM and d flux/dTt (Pt, X fixed), exact. With s = s(T) - R ln P:
//   energy      dT  = (cp(Tt) dTt - M a^2 dM) / (cp(T) + M^2/2 da^2/dT)
//   isentropic  d ln P = cp(T)/(R T) dT - cp(Tt)/(R Tt) dTt
//   d ln rho = d ln P - dT/T,  d ln u = dM/M + (da^2/dT)/(2 a^2) dT.
// Differenced instead, these carried ~3e-8 rounding (flux good to ~3e-15,
// step 1e-7), which the low-drive cancellation in raw_flow amplifies ~1e5.
struct FluxSlopes {
  double dM = 0.0;
  double dTt = 0.0;
};

FluxSlopes isentropic_flux_slopes(double Pt, double Tt, double M, const std::vector<double>& X) {
  const double T = static_T(Tt, M, X);
  const double R = specific_gas_constant(X);
  const double a2 = sound_speed_sq(T, X);
  const double da2 = sound_speed_sq_dT(T, X);
  const double cpT = cp_mass(T, X);
  const double cpTt = cp_mass(Tt, X);
  const double den = cpT + 0.5 * M * M * da2;
  const double flux = isentropic_flux(Pt, Tt, M, X);
  const double dT_dM = -M * a2 / den;
  const double dT_dTt = cpTt / den;
  // d ln flux = (cp/R - 1)/T dT + dM/M + da2/(2 a2) dT - cp(Tt)/(R Tt) dTt
  const double per_dT = (cpT / R - 1.0) / T + 0.5 * da2 / a2;
  FluxSlopes s;
  s.dM = flux * (per_dT * dT_dM + 1.0 / M);
  s.dTt = flux * (per_dT * dT_dTt - cpTt / (R * Tt));
  return s;
}

void check_duct(double L, double D) {
  if (!validation::is_finite_positive(L) || !validation::is_finite_positive(D)) {
    throw std::invalid_argument("fanno (Mach): L and D must be positive and finite");
  }
}

// NaN or inf roughness or f_multiplier used to vanish inside the root
// finders and return a plausible flux (#769). Finite values keep their
// existing meaning (f_multiplier = 0 with "fixed" is frictionless).
void check_friction(double roughness, double f_multiplier) {
  if (!validation::is_finite(roughness) || !validation::is_finite(f_multiplier)) {
    throw std::invalid_argument("fanno (Mach): roughness and f_multiplier must be finite");
  }
}

}  // namespace

FannoMachState fanno_state_at_mach(double G, double Tt, double M,
                                   const std::vector<double>& X) {
  if (!(M > 0.0) || !(G > 0.0) || !(Tt > 0.0)) {
    throw std::invalid_argument("fanno_state_at_mach: G, Tt and M must be positive");
  }
  FannoMachState s;
  s.M = M;
  s.T = static_T(Tt, M, X);
  s.a = std::sqrt(sound_speed_sq(s.T, X));
  s.u = M * s.a;
  s.rho = G / s.u;
  s.P = s.rho * specific_gas_constant(X) * s.T;
  return s;
}

double fanno_dx_dmach(double G, double Tt, double M, const std::vector<double>& X,
                      double D, double roughness, const std::string& friction_model,
                      double f_multiplier) {
  if (!validation::is_finite(D)) {
    throw std::invalid_argument("fanno_dx_dmach: D must be finite");
  }
  check_friction(roughness, f_multiplier);
  const FannoMachState s = fanno_state_at_mach(G, Tt, M, X);
  const double R = specific_gas_constant(X);
  const double a2 = s.a * s.a;
  const double da2_dT = sound_speed_sq_dT(s.T, X);
  // Energy h(T) + M^2 a^2 / 2 = h0, differentiated at fixed h0.
  const double dT_dM = -M * a2 / (cp_mass(s.T, X) + 0.5 * M * M * da2_dT);
  const double du_dM = s.a + M * da2_dT / (2.0 * s.a) * dT_dM;
  const double drho_dM = -G * du_dM / (s.u * s.u);
  const double dP_dM = R * (drho_dM * s.T + s.rho * dT_dM);
  const double Re = G * D / viscosity(s.T, s.P, X);
  const double f = friction(Re, roughness / D, friction_model, f_multiplier);
  return -(dP_dM + G * du_dM) * 2.0 * D / (f * G * s.u);
}

double fanno_length_between(double G, double Tt, double M1, double M2,
                            const std::vector<double>& X, double D, double roughness,
                            const std::string& friction_model, double f_multiplier) {
  if (!(M1 > 0.0) || !(M2 >= M1) || !(M2 <= 1.0)) {
    throw std::invalid_argument("fanno_length_between: need 0 < M1 <= M2 <= 1");
  }
  if (!validation::is_finite(D)) {
    throw std::invalid_argument("fanno_length_between: D must be finite");
  }
  check_friction(roughness, f_multiplier);
  if (M2 == M1) {
    return 0.0;
  }
  // In s = 1/M^2, dx/ds = (M^3 / 2) dx/dM up to sign: bounded at low Mach,
  // where dx/dM itself grows as 1/M^3.
  auto integrand = [&](double s) {
    const double M = 1.0 / std::sqrt(s);
    return 0.5 * M * M * M *
           fanno_dx_dmach(G, Tt, M, X, D, roughness, friction_model, f_multiplier);
  };
  return composite_gl(integrand, 1.0 / (M2 * M2), 1.0 / (M1 * M1));
}

double fanno_sonic_mass_flux(double Pt, double Tt, const std::vector<double>& X) {
  return isentropic_flux(Pt, Tt, 1.0, X);
}

double fanno_inlet_mach(double Pt, double Tt, double G, const std::vector<double>& X) {
  if (G <= 0.0) {
    return 0.0;
  }
  const double G_max = fanno_sonic_mass_flux(Pt, Tt, X);
  if (G >= G_max) {
    return 1.0;
  }
  // Newton on the isentropic flux, which rises monotonically to G_max at
  // M = 1, kept inside the bracket [lo, hi] (bisection when a step leaves it).
  // Starts from the low-Mach limit G ~ G_max M (rho and a near stagnation).
  auto fn = [&](double M) { return isentropic_flux(Pt, Tt, M, X) - G; };
  double lo = 0.0, hi = 1.0;
  double M = std::min(G / G_max, 0.9);
  constexpr int kMaxIter = 60;
  for (int it = 0; it < kMaxIter; ++it) {
    const double f = fn(M);
    if (f > 0.0) {
      hi = M;
    } else {
      lo = M;
    }
    // The exact d flux/dM (#868 R8). The finite-difference stencil it
    // replaces carried ~1e-8 relative rounding, and below M = h it crossed
    // M = 0, read half the slope, and Newton 2-cycled.
    const double fp = isentropic_flux_slopes(Pt, Tt, M, X).dM;
    double M_new = (fp > 0.0) ? M - f / fp : 0.5 * (lo + hi);
    if (!(M_new > lo && M_new < hi)) {
      M_new = 0.5 * (lo + hi);
    }
    if (std::abs(M_new - M) <= 1e-15 * M + 1e-300) {
      return M_new;
    }
    M = M_new;
  }
  // A cap exit is a root only if the bracket has closed onto it (#868).
  if (inner_solve::cap_exit_is_root(hi - lo, M)) {
    return M;
  }
  inner_solve::throw_no_convergence("fanno_inlet_mach", kMaxIter, "relative flux residual",
                                    fn(M) / G);
}

FannoDuctResult fanno_duct(double Pt, double Tt, double G, const std::vector<double>& X,
                           double L, double D, double roughness,
                           const std::string& friction_model, double f_multiplier) {
  check_duct(L, D);
  check_friction(roughness, f_multiplier);
  if (!(G > 0.0)) {
    throw std::invalid_argument("fanno_duct: G must be positive");
  }
  FannoDuctResult out;
  out.G = G;
  const double M_in = fanno_inlet_mach(Pt, Tt, G, X);
  out.inlet = fanno_state_at_mach(G, Tt, M_in, X);
  if (M_in >= 1.0) {
    out.choked = true;
    out.L_star = 0.0;
    out.exit = out.inlet;
    out.Pt_exit = Pt_from_static(out.exit.P, out.exit.T, 1.0, X);
    return out;
  }
  out.L_star = fanno_length_between(G, Tt, M_in, 1.0, X, D, roughness, friction_model,
                                    f_multiplier);
  double M_exit = 1.0;
  if (out.L_star <= L) {
    out.choked = out.L_star < L;
  } else {
    // Exit Mach: the length from M_in reaches L. Newton on the cumulative
    // length, whose slope is dx/dM exactly, inside a bracket.
    auto length_to = [&](double Me) {
      return fanno_length_between(G, Tt, M_in, Me, X, D, roughness, friction_model,
                                  f_multiplier) -
             L;
    };
    M_exit = brent(length_to, M_in, 1.0, -L, out.L_star - L, kFannoMachRelTol);
  }
  out.exit = fanno_state_at_mach(G, Tt, M_exit, X);
  out.Pt_exit = Pt_from_static(out.exit.P, out.exit.T, M_exit, X, 1e-13, 100);
  return out;
}

namespace {

// Physics-derived choked mass-flux scale for bracketing fanno_choked_mass_flux
// and Rayleigh choked_flux (#1337).
//
// In the high-friction (low-inlet-Mach) limit, momentum balance
//   P dP ~ -(f / 2D) G^2 R T dx
// integrated from P(0) = Pt0 to sonic exit P* << Pt0 gives:
//   - Constant / turbulent f: G_turb = Pt0 / sqrt(R T_hi (f_hi L / D + 1)),
//     which also tracks frictionless thermal choking (f = 0 -> Pt0 / sqrt(R T_hi)).
//   - Laminar f = f_mult * (f Re)_0 * mu / (G D):
//     G_lam = Pt0^2 D^2 / (f_mult * (f Re)_0 * mu(T_hi) * R * T_hi * L).
// Taking min(G_turb, G_lam) tracks both turbulent and long-laminar ducts over
// all (L, D, Pt0, Tt) without an arbitrary e-fold cap.
double choked_flux_scale(double Pt0, double Tt0, double Tt_exit, const std::vector<double>& X,
                         double L, double D, double rough, const std::string& model, double fmult,
                         double G_max) {
  const double R = specific_gas_constant(X);
  const double T_hi = std::max(Tt0, Tt_exit);
  double f_hi = 0.0;
  double mu_hi = 0.0;
  if (model == "fixed") {
    f_hi = std::max(fmult, 0.0);
  } else if (fmult > 0.0) {
    mu_hi = viscosity(T_hi, Pt0, X);
    f_hi = friction(G_max * D / mu_hi, rough / D, model, fmult);
  }
  double G_est = Pt0 / std::sqrt(R * T_hi * (f_hi * L / D + 1.0));
  if (model != "fixed" && fmult > 0.0) {
    const double G_lam =
        (Pt0 * Pt0 * D * D) / (fmult * laminar_friction::f_re * mu_hi * R * T_hi * L);
    G_est = std::min(G_est, G_lam);
  }
  return G_est;
}

}  // namespace

double fanno_choked_mass_flux(double Pt, double Tt, const std::vector<double>& X,
                              double L, double D, double roughness,
                              const std::string& friction_model, double f_multiplier) {
  check_duct(L, D);
  check_friction(roughness, f_multiplier);
  const double G_max = fanno_sonic_mass_flux(Pt, Tt, X);
  // L*(G) falls monotonically from infinity (G -> 0) to 0 (G -> G_max).
  auto excess = [&](double logG) {
    const double G = std::exp(logG);
    const double M_in = fanno_inlet_mach(Pt, Tt, G, X);
    if (M_in >= 1.0) {
      return -L;
    }
    return fanno_length_between(G, Tt, M_in, 1.0, X, D, roughness, friction_model,
                                f_multiplier) -
           L;
  };
  double hi = std::log(G_max * (1.0 - 1e-12));
  double f_hi = excess(hi);
  if (f_hi >= 0.0) {
    return G_max;
  }
  // Lowest log(G) at which 1 / M_in^2 ~ (G_max / G)^2 stays finite in double.
  const double log_G_min = std::log(G_max) + 0.5 * std::log(std::numeric_limits<double>::min());
  const double G_est =
      choked_flux_scale(Pt, Tt, Tt, X, L, D, roughness, friction_model, f_multiplier, G_max);
  double lo = std::clamp(std::log(G_est) - 1.0, log_G_min, hi - 1.0);
  double f_lo = excess(lo);
  double step = 2.0;
  while (f_lo <= 0.0 && lo > log_G_min) {
    hi = lo;
    f_hi = f_lo;
    lo = std::max(lo - step, log_G_min);
    f_lo = excess(lo);
    step *= 2.0;
  }
  return std::exp(brent(excess, lo, hi, f_lo, f_hi, 1e-14));
}

namespace {

// Exit pressure per unit mass flux: at fixed (M, Tt) the station state's P is
// G R T / (M a), linear in G, and so is its stagnation pressure (P times the
// isentropic ratio, which depends on T and M only). So for a target pressure
// the flux at any exit Mach is explicit: G = P_target / exit_per_flux(M).
double exit_per_flux(double Me, double Tt, const std::vector<double>& X, bool exit_total) {
  const FannoMachState s = fanno_state_at_mach(1.0, Tt, Me, X);
  if (!exit_total) {
    return s.P;
  }
  const double R = specific_gas_constant(X);
  return s.P * std::exp((s_mass(Tt, X, 101325.0) - s_mass(s.T, X, 101325.0)) / R);
}

// d exit_per_flux/dM and d/dTt, exact (same energy relation as
// isentropic_flux_slopes): ln phi = ln(R T) - ln M - ln(a^2)/2, times the
// isentropic ratio exp((s(Tt) - s(T))/R) for the exit stagnation pressure.
// Differenced, these left ~1e-10 noise, which the low-drive cancellation in
// raw_flow amplified into the drive-floor exponent sigma -- 1e-5 noise in G
// below the floor.
struct PhiSlopes {
  double dM = 0.0;
  double dTt = 0.0;
};

PhiSlopes exit_per_flux_slopes(double Me, double Tt, const std::vector<double>& X,
                               bool exit_total) {
  const double phi = exit_per_flux(Me, Tt, X, exit_total);
  const double T = static_T(Tt, Me, X);
  const double R = specific_gas_constant(X);
  const double a2 = sound_speed_sq(T, X);
  const double da2 = sound_speed_sq_dT(T, X);
  const double cpT = cp_mass(T, X);
  const double cpTt = cp_mass(Tt, X);
  const double den = cpT + 0.5 * Me * Me * da2;
  const double dT_dM = -Me * a2 / den;
  const double dT_dTt = cpTt / den;
  double per_dT = 1.0 / T - 0.5 * da2 / a2;
  double per_dTt = 0.0;
  if (exit_total) {
    per_dT -= cpT / (R * T);
    per_dTt = cpTt / (R * Tt);
  }
  PhiSlopes s;
  s.dM = phi * (per_dT * dT_dM - 1.0 / Me);
  s.dTt = phi * (per_dT * dT_dTt + per_dTt);
  return s;
}

// ---------------------------------------------------------------------------
// Exact tangents of the Fanno pieces at fixed Mach limits (#1311, #1310).
//
// raw_flow's implicit derivatives combine the pieces below (the exit pressure
// per unit flux, the isentropic inlet flux, the length between two Mach
// numbers) in sums that cancel to drive / P at low drive (5e-6 at 1 Pa on
// 2 bar), twice for a composition direction. Differenced, each piece carried
// ~1e-8 relative error, which the cancellation amplified to 1-3 % in dG/dY
// near the 1 Pa floor (#1311), whatever the step: a larger one trades
// round-off for truncation. Forward mode carries round-off only.
//
// Slots: a DualN<N> whose slot i is a state seed (G, Tt0) has comp[i] = 0; a
// composition slot has comp[i] = 1 and takes the X partials along dX. A
// direction v in the mass fractions moves the mole fractions by
//   dX_i = MW (v_i / W_i - X_i sum_k v_k / W_k),
// the step mass_to_mole(Y + t v) takes (it keeps sum X = 1). The mixture's
// molar cp, h, dcp/dT, d2cp/dT2 and MW are linear in X, so each one's
// derivative along dX is the same sum with dX in place of X. For an absent
// species it is the one-sided derivative as the species enters, the limit
// d_ddir's one-sided stencils approximated.
// ---------------------------------------------------------------------------

template <int N>
using DualT = solver::DualN<N>;

template <int N>
struct MixTangent {
  const std::vector<double>* X = nullptr;
  std::vector<double> dX;   // along the composition slot's direction (empty: none)
  double mw = 0.0;          // [g/mol]
  double mw_v = 0.0;        // d MW along dX
  std::array<double, N> comp{};

  bool has_comp() const { return !dX.empty(); }
};

template <int N>
MixTangent<N> mix_tangent(const std::vector<double>& X, const std::vector<double>* v, int slot) {
  MixTangent<N> m;
  m.X = &X;
  m.mw = mwmix(X);
  if (v == nullptr) {
    return m;
  }
  if (v->size() != X.size()) {
    throw std::invalid_argument("fanno_channel_flow: a dY direction has the wrong length");
  }
  double s = 0.0;
  for (std::size_t k = 0; k < X.size(); ++k) {
    s += (*v)[k] / species_molar_mass(k);
  }
  m.dX.assign(X.size(), 0.0);
  for (std::size_t k = 0; k < X.size(); ++k) {
    m.dX[k] = m.mw * ((*v)[k] / species_molar_mass(k) - X[k] * s);
  }
  m.mw_v = mwmix(m.dX);
  m.comp[slot] = 1.0;
  return m;
}

// val with slope d_dT along T and d_dv along the composition slot.
template <int N>
DualT<N> lift_tx(double val, double d_dT, const DualT<N>& T, double d_dv, const MixTangent<N>& m) {
  DualT<N> r;
  r.v = val;
  for (int i = 0; i < N; ++i) {
    r.d[i] = d_dT * T.d[i] + d_dv * m.comp[i];
  }
  return r;
}

template <int N>
DualT<N> gas_constant_t(const MixTangent<N>& m) {
  const double R = specific_gas_constant(*m.X);
  return lift_tx(R, 0.0, DualT<N>::constant(0.0), m.has_comp() ? -R * m.mw_v / m.mw : 0.0, m);
}

// Molar cp [J/(mol K)] and its T-slope, each with its composition partial.
template <int N>
DualT<N> cp_molar_t(const DualT<N>& T, const MixTangent<N>& m) {
  return lift_tx(cp(T.v, *m.X), dcp_dT(T.v, *m.X), T, m.has_comp() ? cp(T.v, m.dX) : 0.0, m);
}

template <int N>
DualT<N> dcp_molar_dT_t(const DualT<N>& T, const MixTangent<N>& m) {
  return lift_tx(dcp_dT(T.v, *m.X), d2cp_dT2(T.v, *m.X), T,
                 m.has_comp() ? dcp_dT(T.v, m.dX) : 0.0, m);
}

// A mass-based property from its molar one: q / MW * 1000.
template <int N>
DualT<N> per_mass_t(const DualT<N>& q_molar, const MixTangent<N>& m) {
  DualT<N> r = q_molar * (1000.0 / m.mw);
  if (m.has_comp()) {
    for (int i = 0; i < N; ++i) {
      r.d[i] -= r.v * (m.mw_v / m.mw) * m.comp[i];
    }
  }
  return r;
}

template <int N>
DualT<N> h_mass_t(const DualT<N>& T, const MixTangent<N>& m) {
  const DualT<N> h_mol =
      lift_tx(h(T.v, *m.X), cp(T.v, *m.X), T, m.has_comp() ? h(T.v, m.dX) : 0.0, m);
  return per_mass_t(h_mol, m);
}

template <int N>
DualT<N> cp_mass_t(const DualT<N>& T, const MixTangent<N>& m) {
  return per_mass_t(cp_molar_t(T, m), m);
}

// a^2 = gamma R T and d(a^2)/dT = R (gamma + T dgamma/dT), gamma = cp/(cp - R_u).
template <int N>
DualT<N> sound_speed_sq_t(const DualT<N>& T, const MixTangent<N>& m) {
  const DualT<N> c = cp_molar_t(T, m);
  return c / (c - thermo::R_GAS) * gas_constant_t(m) * T;
}

template <int N>
DualT<N> sound_speed_sq_dT_t(const DualT<N>& T, const MixTangent<N>& m) {
  const DualT<N> c = cp_molar_t(T, m);
  const DualT<N> cv = c - thermo::R_GAS;
  const DualT<N> dgamma = -thermo::R_GAS * dcp_molar_dT_t(T, m) / (cv * cv);
  return gas_constant_t(m) * (c / cv + T * dgamma);
}

// (s(T1) - s(T2)) / R at fixed composition and pressure: the mixing terms
// cancel, so it is sum_k X_k (s_k(T1) - s_k(T2)) / R_u, linear in X.
template <int N>
DualT<N> entropy_rise_over_R_t(const DualT<N>& T1, const DualT<N>& T2, const MixTangent<N>& m) {
  const std::vector<double>& X = *m.X;
  const double val =
      (s_mass(T1.v, X, 101325.0) - s_mass(T2.v, X, 101325.0)) / specific_gas_constant(X);
  double dv = 0.0;
  if (m.has_comp()) {
    for (std::size_t k = 0; k < X.size(); ++k) {
      if (m.dX[k] != 0.0) {
        dv += m.dX[k] * (s_species(k, T1.v) - s_species(k, T2.v));
      }
    }
    dv /= thermo::R_GAS;
  }
  DualT<N> r;
  r.v = val;
  const double s1 = cp(T1.v, X) / (thermo::R_GAS * T1.v);
  const double s2 = cp(T2.v, X) / (thermo::R_GAS * T2.v);
  for (int i = 0; i < N; ++i) {
    r.d[i] = s1 * T1.d[i] - s2 * T2.d[i] + dv * m.comp[i];
  }
  return r;
}

// static_T with its tangent: h(T) + M^2 a^2(T) / 2 = h(Tt), differentiated
// implicitly. The value is static_T's own.
template <int N>
DualT<N> static_T_t(const DualT<N>& Tt, const DualT<N>& M, const MixTangent<N>& m) {
  const std::vector<double>& X = *m.X;
  const double T_v = static_T(Tt.v, M.v, X);
  const DualT<N> T0 = DualT<N>::constant(T_v);
  // F at the root with T held: its tangent is everything but the T term.
  const DualT<N> F = h_mass_t(T0, m) + 0.5 * M * M * sound_speed_sq_t(T0, m) - h_mass_t(Tt, m);
  const double F_T = cp_mass(T_v, X) + 0.5 * M.v * M.v * sound_speed_sq_dT(T_v, X);
  DualT<N> T;
  T.v = T_v;
  for (int i = 0; i < N; ++i) {
    T.d[i] = -F.d[i] / F_T;
  }
  return T;
}

template <int N>
DualT<N> exit_per_flux_t(const DualT<N>& Me, const DualT<N>& Tt, const MixTangent<N>& m,
                         bool exit_total) {
  const DualT<N> T = static_T_t(Tt, Me, m);
  const DualT<N> u = Me * solver::dsqrt(sound_speed_sq_t(T, m));
  const DualT<N> P = gas_constant_t(m) * T / u;
  if (!exit_total) {
    return P;
  }
  return P * solver::dexp(entropy_rise_over_R_t(Tt, T, m));
}

template <int N>
DualT<N> isentropic_flux_t(double Pt, const DualT<N>& Tt, const DualT<N>& M,
                           const MixTangent<N>& m) {
  const DualT<N> T = static_T_t(Tt, M, m);
  const DualT<N> P = Pt * solver::dexp(entropy_rise_over_R_t(T, Tt, m));
  return P / (gas_constant_t(m) * T) * M * solver::dsqrt(sound_speed_sq_t(T, m));
}

// fanno_dx_dmach at a fixed M, with tangents in G, Tt and composition.
template <int N>
DualT<N> fanno_dx_dmach_t(const DualT<N>& G, const DualT<N>& Tt, double M,
                          const MixTangent<N>& m, double D, double roughness,
                          const std::string& friction_model, double f_multiplier) {
  const DualT<N> Md = DualT<N>::constant(M);
  const DualT<N> R = gas_constant_t(m);
  const DualT<N> T = static_T_t(Tt, Md, m);
  const DualT<N> a2 = sound_speed_sq_t(T, m);
  const DualT<N> a = solver::dsqrt(a2);
  const DualT<N> u = M * a;
  const DualT<N> rho = G / u;
  const DualT<N> P = rho * R * T;
  const DualT<N> da2_dT = sound_speed_sq_dT_t(T, m);
  const DualT<N> dT_dM = -M * a2 / (cp_mass_t(T, m) + 0.5 * M * M * da2_dT);
  const DualT<N> du_dM = a + M * da2_dT / (2.0 * a) * dT_dM;
  const DualT<N> drho_dM = -G * du_dM / (u * u);
  const DualT<N> dP_dM = R * (drho_dM * T + rho * dT_dM);
  DualT<N> f = DualT<N>::constant(f_multiplier);
  if (friction_model != "fixed") {
    const TransportDerivatives tr = transport_and_dT(T.v, P.v, *m.X);
    double dmu_dv = 0.0;
    if (m.has_comp()) {
      const std::vector<double> g = viscosity_dX(T.v, P.v, *m.X);
      for (std::size_t k = 0; k < g.size(); ++k) {
        dmu_dv += g[k] * m.dX[k];
      }
    }
    const DualT<N> mu = lift_tx(tr.mu, tr.dmu_dT, T, dmu_dv, m);
    const DualT<N> Re = G * D / mu;
    const auto [f0, df_dRe] =
        solver::friction_and_jacobian(friction_model, Re.v, roughness / D).result;
    f.v = f0 * f_multiplier;
    for (int i = 0; i < N; ++i) {
      f.d[i] = df_dRe * f_multiplier * Re.d[i];
    }
  }
  return -(dP_dM + G * du_dM) * 2.0 * D / (f * G * u);
}

// fanno_length_between at fixed Mach limits: the same fixed rule applied to
// the integrand's tangent is the exact tangent of the rule's sum.
template <int N>
DualT<N> fanno_length_between_t(const DualT<N>& G, const DualT<N>& Tt, double M1, double M2,
                                const MixTangent<N>& m, double D, double roughness,
                                const std::string& friction_model, double f_multiplier) {
  // The value function's checks, so the tangent refuses what it refuses: a
  // frictionless duct whose inlet Mach passes the exit's by round-off threw
  // there and must not divide by f = 0 here (#1354 review).
  if (!(M1 > 0.0) || !(M2 >= M1) || !(M2 <= 1.0)) {
    throw std::invalid_argument("fanno_length_between: need 0 < M1 <= M2 <= 1");
  }
  if (!validation::is_finite(D)) {
    throw std::invalid_argument("fanno_length_between: D must be finite");
  }
  check_friction(roughness, f_multiplier);
  if (M2 == M1) {
    return DualT<N>::constant(0.0);
  }
  const double a = 1.0 / (M2 * M2);
  const double b = 1.0 / (M1 * M1);
  const double ratio = std::pow(b / a, 1.0 / kFannoPanels);
  DualT<N> total = DualT<N>::constant(0.0);
  double lo = a;
  for (int k = 0; k < kFannoPanels; ++k) {
    const double hi = (k == kFannoPanels - 1) ? b : lo * ratio;
    const double c = 0.5 * (lo + hi);
    const double hw = 0.5 * (hi - lo);
    DualT<N> acc = DualT<N>::constant(0.0);
    for (std::size_t j = 0; j < kGlX.size(); ++j) {
      for (const double s : {c - hw * kGlX[j], c + hw * kGlX[j]}) {
        const double M = 1.0 / std::sqrt(s);
        acc = acc + kGlW[j] * 0.5 * M * M * M *
                        fanno_dx_dmach_t(G, Tt, M, m, D, roughness, friction_model,
                                         f_multiplier);
      }
    }
    total = total + acc * hw;
    lo = hi;
  }
  return total;
}

struct FlowProblem {
  double Pt0, Tt0, P, L, D, rough, fmult;
  const std::vector<double>* X;
  const std::string* model;
  bool exit_total;
  double Me_guess;
  const std::vector<std::vector<double>>* dirs;
  bool with_derivatives;
};

// Length the flux G(Me) = P / exit_per_flux(Me) needs to reach exit Mach Me,
// less L. Falls monotonically with Me: positive for a slow exit (a long duct
// would be needed), negative once G passes the choked flux. A flux the inlet
// cannot supply (sonic inlet) needs no length at all: -L.
double psi(double Me, double P, double Pt0, double Tt0, const FlowProblem& pb) {
  const std::vector<double>& X = *pb.X;
  const double G = P / exit_per_flux(Me, Tt0, X, pb.exit_total);
  if (G >= fanno_sonic_mass_flux(Pt0, Tt0, X)) {
    return -pb.L;
  }
  const double Mi = fanno_inlet_mach(Pt0, Tt0, G, X);
  if (Mi >= Me) {
    return -pb.L;
  }
  return fanno_length_between(G, Tt0, Mi, Me, X, pb.D, pb.rough, *pb.model, pb.fmult) - pb.L;
}

// d g(X)/ds along each direction v in mass-fraction space, Y -> Y + s v (X
// renormalised). An element chains it with the solver's relay: the feeding
// node's dY/dx all lie in the span of the directions it passes (an
// orthonormal basis, usually ONE vector for two mixing streams), so the
// projected gradient sum_i (dg/dv_i) v_i reproduces dg/dY . dY/dx exactly at
// a fraction of a per-species gradient's cost. A unit vector e_k gives the
// per-species partial with the other mass fractions fixed.
//
// Central differences of a smooth function at fixed (G or Me), the largest
// component of the step kFannoDirStep (times step_scale, which a caller near a
// singularity shrinks). Where the central step would take a fraction below
// zero (a species absent or at trace level here but carried by an inflow), a
// one-sided second-order difference, O(h^2) like the central one, at half
// the step where the doubled step is infeasible too. Where neither side is
// feasible (two absent or trace species with opposite-sign components), the
// direction is split into its blocking species and the rest, each of which
// has a feasible stencil, and the parts are summed (dg/dv is linear in v).
//
// A component at or below kFannoDirRoundoff of the largest that would take its
// fraction below zero is dropped from the step first: it contributes nothing
// measurable to the derivative, but kept, its sign alone chose the stencil. A
// -2.2e-20 on an absent species switched a duct's dG/dv to one-sided O(h)
// and moved J by 2.3e-5 relative (#1199).
template <class Fn>
std::vector<double> d_ddir(const std::vector<double>& X,
                           const std::vector<std::vector<double>>& dirs, Fn&& g,
                           double step_scale = 1.0) {
  const std::vector<double> Y = mole_to_mass(X);
  std::vector<double> out(dirs.size(), 0.0);
  // g at X itself, for the one-sided stencils: once per call, if at all.
  bool have_g0 = false;
  double g0 = 0.0;
  auto g_at_x = [&]() {
    if (!have_g0) {
      g0 = g(X);
      have_g0 = true;
    }
    return g0;
  };
  auto at = [&](const std::vector<double>& v, double s) {
    std::vector<double> Ys = Y;
    for (std::size_t k = 0; k < Y.size(); ++k) {
      Ys[k] += s * v[k];
    }
    return g(mass_to_mole(Ys));
  };
  auto feasible = [&](const std::vector<double>& v, double s) {
    for (std::size_t k = 0; k < Y.size(); ++k) {
      if (Y[k] + s * v[k] < 0.0) {
        return false;
      }
    }
    return true;
  };
  std::vector<double> v;
  for (std::size_t i = 0; i < dirs.size(); ++i) {
    v = dirs[i];
    if (v.size() != Y.size()) {
      throw std::invalid_argument("fanno_channel_flow: a dY direction has the wrong length");
    }
    double vmax = 0.0;
    for (double vk : v) {
      vmax = std::max(vmax, std::abs(vk));
    }
    if (!(vmax > 0.0)) {
      continue;
    }
    const double h = step_scale * kFannoDirStep / vmax;
    for (std::size_t k = 0; k < Y.size(); ++k) {
      if (std::abs(v[k]) <= kFannoDirRoundoff * vmax && Y[k] - h * std::abs(v[k]) < 0.0) {
        v[k] = 0.0;
      }
    }
    if (feasible(v, h) && feasible(v, -h)) {
      out[i] = (at(v, h) - at(v, -h)) / (2.0 * h);
      continue;
    }
    // One-sided towards the feasible side s = +-h, second order:
    // (-3 g(0) + 4 g(s) - g(2s)) / (2s). Where 2s is infeasible too (a trace
    // species blocks it), the same stencil at half the step, which needs only
    // s: the feasible set {t : Y + t v >= 0} is an interval holding 0, so s/2
    // is feasible with s. A first-order fallback here put a 1.2e-5 jump into
    // dG/dv as the trace fraction crossed 2h|v_k|.
    bool done = false;
    for (const double s : {h, -h}) {
      if (!feasible(v, s)) {
        continue;
      }
      const double t = feasible(v, 2.0 * s) ? s : 0.5 * s;
      out[i] = (-3.0 * g_at_x() + 4.0 * at(v, t) - at(v, 2.0 * t)) / (2.0 * t);
      done = true;
      break;
    }
    if (done) {
      continue;
    }
    // Neither side feasible: v has opposite-sign components on two species
    // that are absent or below h|v_k| (#1203). It returned 0 here, a J entry
    // of 0 where the true dG/dv was 58-66. dg/dv is linear in v, so split it:
    // v = sum_b v_b e_b + w, with b each species that blocks a step of 2h
    // either way and w the rest. A unit vector e_b is feasible on its
    // positive side. w has no component blocked at h, but is stepped at its
    // own (larger) h, where another component can block it; the recursion
    // splits that too, and each level drops at least one component of w.
    //
    // Only for a state on the simplex. A negative fraction (here only from a
    // direct kernel call; the network clamps Y to [0, 1]) blocks every
    // direction, so no split shrinks and it would recurse forever: such a
    // direction keeps the 0 it had.
    bool on_simplex = true;
    for (const double yk : Y) {
      on_simplex = on_simplex && yk >= 0.0;
    }
    if (!on_simplex) {
      continue;
    }
    std::vector<std::vector<double>> parts;
    std::vector<double> coeff;
    std::vector<double> w = v;
    for (std::size_t k = 0; k < Y.size(); ++k) {
      if (v[k] != 0.0 && Y[k] < 2.0 * h * std::abs(v[k])) {
        std::vector<double> e(Y.size(), 0.0);
        e[k] = 1.0;
        parts.push_back(std::move(e));
        coeff.push_back(v[k]);
        w[k] = 0.0;
      }
    }
    parts.push_back(w);
    coeff.push_back(1.0);
    const std::vector<double> sub = d_ddir(X, parts, g, step_scale);
    double sum = 0.0;
    for (std::size_t p = 0; p < sub.size(); ++p) {
      sum += coeff[p] * sub[p];
    }
    out[i] = sum;
  }
  return out;
}

// The raw (unregularised) solve at drive Pt0 - P > 0.
FannoChannelFlow raw_flow(const FlowProblem& pb) {
  const std::vector<double>& X = *pb.X;
  FannoChannelFlow out;
  // A back pressure at or below zero (a solver probe; vacuum at best) chokes
  // any duct: the choked branch, which P_target -> 0+ already reaches. psi
  // itself would build a non-positive flux there.
  const double choked_test = (pb.P > 0.0) ? psi(1.0, pb.P, pb.Pt0, pb.Tt0, pb) : 0.0;
  if (choked_test >= 0.0) {
    // Choked: even a sonic exit at this back pressure would need more duct
    // than there is, so the exit sits at sonic above P and the flux is the
    // choked one, set by the inlet state alone.
    out.choked = true;
    out.M_exit = 1.0;
    out.G = fanno_choked_mass_flux(pb.Pt0, pb.Tt0, X, pb.L, pb.D, pb.rough, *pb.model, pb.fmult);
    out.M_in = fanno_inlet_mach(pb.Pt0, pb.Tt0, out.G, X);
    if (!pb.with_derivatives) {
      return out;
    }
    out.dG_dP_target = 0.0;
    const double G = out.G;
    const double Mi = out.M_in;
    const FluxSlopes fs = isentropic_flux_slopes(pb.Pt0, pb.Tt0, Mi, X);
    if (!(Mi < 1.0) || !(fs.dM > 0.0)) {
      // The duct is too short to hold any friction at all: the flux is the
      // sonic one of the inlet state, G_max(Pt0, Tt0, Y).
      out.dG_dPt0 = G / pb.Pt0;
      out.dG_dTt0 = isentropic_flux_slopes(pb.Pt0, pb.Tt0, 1.0, X).dTt;
      if (!pb.dirs->empty()) {
        out.dG_ddir.assign(pb.dirs->size(), 0.0);
        for (std::size_t i = 0; i < pb.dirs->size(); ++i) {
          out.dG_ddir[i] = isentropic_flux_t(pb.Pt0, DualT<1>::constant(pb.Tt0),
                                             DualT<1>::constant(1.0),
                                             mix_tangent<1>(X, &(*pb.dirs)[i], 0))
                               .d[0];
        }
      }
      return out;
    }
    // Lam = L*(G; Pt0, Tt0, Y) - L = 0, differentiated implicitly the same way
    // as the unchoked psi below (Leibniz): the inlet limit M_in(G, Pt0, Tt0, Y)
    // of the length to sonic is taken exactly, -dx/dM(M_in) dM_in/d(.), with
    // dM_in/d(.) from the isentropic inlet flux(Pt0, Tt0, M_in, Y) = G, and
    // the length at FIXED limits carries its own exact tangent (in G, Tt0, Y).
    //
    // Differencing L* through a re-solved M_in instead (#1019) broke on a
    // short duct: there the choked flux sits within ~1e-6 of the sonic one
    // (1 - M_in ~ sqrt(f L / D)), so a 1e-6 step in G, Tt0 or Y took the
    // stencil past sonic, where M_in clamps to 1 and L* to 0. At L = 1 um
    // dG/dTt0 came out 17% high and dG/dY up to 65%.
    const std::string& model = *pb.model;
    const double dx_Mi = fanno_dx_dmach(G, pb.Tt0, Mi, X, pb.D, pb.rough, model, pb.fmult);
    // flux is proportional to Pt0 at fixed (Tt0, M), so flux_Pt = G / Pt0.
    const double dMi_dG = 1.0 / fs.dM;
    const double dMi_dPt = -(G / pb.Pt0) / fs.dM;
    const double dMi_dT = -fs.dTt / fs.dM;
    // The length at fixed limits, exact in G and Tt0 (#1311).
    const MixTangent<2> mix_s = mix_tangent<2>(X, nullptr, 0);
    const DualT<2> ell_s =
        fanno_length_between_t(DualT<2>::seed(G, 0), DualT<2>::seed(pb.Tt0, 1), Mi, 1.0, mix_s,
                               pb.D, pb.rough, model, pb.fmult);
    const double Lam_G = ell_s.d[0] - dx_Mi * dMi_dG;
    const double Lam_T = ell_s.d[1] - dx_Mi * dMi_dT;
    const double Lam_Pt = -dx_Mi * dMi_dPt;
    out.dG_dPt0 = -Lam_Pt / Lam_G;
    out.dG_dTt0 = -Lam_T / Lam_G;
    if (pb.dirs->empty()) {
      return out;
    }
    // Along each direction v: Lam_v = ell_v + dx/dM(M_in) flux_v / flux_M.
    const std::vector<std::vector<double>>& dirs = *pb.dirs;
    out.dG_ddir.assign(dirs.size(), 0.0);
    for (std::size_t i = 0; i < dirs.size(); ++i) {
      const MixTangent<1> mix_v = mix_tangent<1>(X, &dirs[i], 0);
      const DualT<1> Tt0_c = DualT<1>::constant(pb.Tt0);
      const double ell_v = fanno_length_between_t(DualT<1>::constant(G), Tt0_c, Mi, 1.0, mix_v,
                                                  pb.D, pb.rough, model, pb.fmult)
                               .d[0];
      const double flux_v =
          isentropic_flux_t(pb.Pt0, Tt0_c, DualT<1>::constant(Mi), mix_v).d[0];
      out.dG_ddir[i] = -(ell_v + dx_Mi * flux_v / fs.dM) / Lam_G;
    }
    return out;
  }
  // Unchoked: the exit Mach is the root of psi, which falls monotonically in
  // Me. A caller's guess (the previous solve) gives a tight bracket, widened
  // geometrically until it holds a sign change; without one the bracket is
  // grown down from Me = 1 by halving.
  //
  // IEEE-754 machine-range guard for the cold halving loop (#1242):
  // fanno_length_between integrates in s = 1 / M^2 up to 1 / M_in^2, and at
  // low Mach M_in ~ Me * (P / Pt0), so 1 / M_in^2 <= DBL_MAX requires
  // Me >= (Pt0 / P) / sqrt(DBL_MAX). Keeping the warm-start loop identical to
  // main preserves bit-identity on every ordinary duct (Me > 2^-20, cold and
  // warm); only the cold halving loop extends below 1e-6 down to kMachMin.
  const double inv_sqrt_max = 1.0 / std::sqrt(std::numeric_limits<double>::max());
  const double kMachMin = (pb.P > 0.0) ? std::max((pb.Pt0 / pb.P) * inv_sqrt_max, inv_sqrt_max)
                                       : inv_sqrt_max;
  auto f = [&](double Me) { return psi(Me, pb.P, pb.Pt0, pb.Tt0, pb); };
  double lo = 0.0, hi = 0.0, f_lo = 0.0, f_hi = 0.0;
  bool bracketed = false;
  if (pb.Me_guess > 0.0 && pb.Me_guess < 1.0) {
    double w = 0.01;
    lo = pb.Me_guess * (1.0 - w);
    hi = std::min(1.0, pb.Me_guess * (1.0 + w));
    f_lo = f(lo);
    f_hi = f(hi);
    while ((f_lo <= 0.0 || f_hi > 0.0) && w < 1.0) {
      w *= 8.0;
      if (f_lo <= 0.0) {
        hi = lo;
        f_hi = f_lo;
        lo = std::max(pb.Me_guess * (1.0 - w), 1e-6);
        f_lo = f(lo);
      } else {
        lo = hi;
        f_lo = f_hi;
        hi = std::min(1.0, pb.Me_guess * (1.0 + w));
        f_hi = f(hi);
      }
    }
    bracketed = f_lo > 0.0 && f_hi <= 0.0;
  }
  // No guess, or a guess the widening could not bracket: it stops at about
  // 3.6x the guess, so a root above that was never enclosed and Brent ran on
  // an interval without a sign change. A short duct warm-started at M 0.01
  // then returned the same flux for every drive (#585).
  if (!bracketed) {
    hi = 1.0;
    f_hi = choked_test;
    lo = 0.5;
    f_lo = f(lo);
  }
  while (f_lo <= 0.0 && lo > kMachMin) {
    hi = lo;
    f_hi = f_lo;
    lo = std::max(lo * 0.5, kMachMin);
    f_lo = f(lo);
  }
  const double Me = brent(f, lo, hi, f_lo, f_hi, 1e-13);
  const double phi = exit_per_flux(Me, pb.Tt0, X, pb.exit_total);
  out.M_exit = Me;
  out.G = pb.P / phi;
  out.M_in = fanno_inlet_mach(pb.Pt0, pb.Tt0, out.G, X);
  if (!pb.with_derivatives) {
    return out;
  }
  // psi(Me; P, Pt0, Tt0) = 0, differentiated implicitly. Its slope in Me stays
  // finite at sonic -- dx/dM -> 0 there, but G(Me) still moves -- so this is
  // well conditioned right up to the choke, where it meets the choked branch.
  //
  // psi(Me, P) = Lam(G, Me) with G = P / phi(Me, Tt0) and
  // Lam = (length from M_in(G, Pt0, Tt0, Y) to Me) - L. Leibniz: the limits
  // are taken exactly -- d Lam/d Me = dx/dM at Me, and every parameter that
  // moves M_in contributes -dx/dM(M_in) dM_in/d(.), with dM_in/d(.) from the
  // isentropic inlet flux(Pt0, Tt0, M_in, Y) = G -- and only the length at
  // FIXED limits carries its own exact tangent (in G, Tt0, Y), with the smooth
  // inlet flux.
  //
  // Differencing the length through a re-solved M_in instead is ill
  // conditioned at low drive: (Me - M_in)/Me ~ drive/P, 7e-6 at 1.5 Pa on
  // 2 bar, which amplifies the inlet solve's ~1e-15 noise to ~4e-10 in the
  // length. A 1e-6 step then left 4e-4 noise in Lam_G -- a 5e-4 jump in G
  // below the drive floor, through sigma -- and several % in dG/dY, where
  // two O(G) terms cancel to 0.4% of G at low Mach.
  const double G = out.G;
  const double Mi = out.M_in;
  const std::string& model = *pb.model;
  const double dx_Mi = fanno_dx_dmach(G, pb.Tt0, Mi, X, pb.D, pb.rough, model, pb.fmult);
  const FluxSlopes fs = isentropic_flux_slopes(pb.Pt0, pb.Tt0, Mi, X);
  const double flux_M = fs.dM;
  const double flux_T = fs.dTt;
  // flux is proportional to Pt0 at fixed (Tt0, M), so flux_Pt = G / Pt0.
  const double dMi_dG = 1.0 / flux_M;
  const double dMi_dPt = -(G / pb.Pt0) / flux_M;
  const double dMi_dT = -flux_T / flux_M;
  // The length at fixed limits, with exact tangents in G and Tt0 (#1311;
  // central differences at 1e-6 until then).
  const MixTangent<2> mix_s = mix_tangent<2>(X, nullptr, 0);
  const DualT<2> ell_s =
      fanno_length_between_t(DualT<2>::seed(G, 0), DualT<2>::seed(pb.Tt0, 1), Mi, Me, mix_s,
                             pb.D, pb.rough, model, pb.fmult);
  const double Lam_G = ell_s.d[0] - dx_Mi * dMi_dG;
  const double Lam_T = ell_s.d[1] - dx_Mi * dMi_dT;
  const double Lam_Pt = -dx_Mi * dMi_dPt;
  const double Lam_Me = fanno_dx_dmach(G, pb.Tt0, Me, X, pb.D, pb.rough, model, pb.fmult);
  // G = P / phi(Me, Tt0)
  const PhiSlopes ps = exit_per_flux_slopes(Me, pb.Tt0, X, pb.exit_total);
  const double phi_M = ps.dM;
  const double phi_T = ps.dTt;
  const double G_P = 1.0 / phi;
  const double G_M = -pb.P * phi_M / (phi * phi);
  const double G_T = -pb.P * phi_T / (phi * phi);
  const double psi_M = Lam_G * G_M + Lam_Me;
  const double psi_P = Lam_G * G_P;
  const double psi_Pt = Lam_Pt;
  const double psi_T = Lam_G * G_T + Lam_T;
  const double dMe_dP = -psi_P / psi_M;
  const double dMe_dPt = -psi_Pt / psi_M;
  const double dMe_dT = -psi_T / psi_M;
  out.dG_dP_target = G_P + G_M * dMe_dP;
  out.dG_dPt0 = G_M * dMe_dPt;
  out.dG_dTt0 = G_T + G_M * dMe_dT;
  if (pb.dirs->empty()) {
    return out;
  }
  // Composition, the same way along each direction v: psi_v = Lam_G G_v +
  // Lam_v at fixed Me, with G_v = -P phi_v / phi^2 and
  // Lam_v = ell_v - dx/dM(M_in) dM_in/dv.
  // Each piece's tangent along v is exact (#1311): d_ddir's differences of
  // them, amplified by the cancellation above, left 1-3 % in dG/dY near the
  // 1 Pa floor.
  const std::vector<std::vector<double>>& dirs = *pb.dirs;
  out.dG_ddir.assign(dirs.size(), 0.0);
  for (std::size_t i = 0; i < dirs.size(); ++i) {
    const MixTangent<1> mix_v = mix_tangent<1>(X, &dirs[i], 0);
    const DualT<1> Me_c = DualT<1>::constant(Me);
    const DualT<1> Tt0_c = DualT<1>::constant(pb.Tt0);
    const double phi_v = exit_per_flux_t(Me_c, Tt0_c, mix_v, pb.exit_total).d[0];
    const double ell_v = fanno_length_between_t(DualT<1>::constant(G), Tt0_c, Mi, Me, mix_v,
                                                pb.D, pb.rough, model, pb.fmult)
                             .d[0];
    const double flux_v =
        isentropic_flux_t(pb.Pt0, Tt0_c, DualT<1>::constant(Mi), mix_v).d[0];
    const double G_v = -pb.P * phi_v / (phi * phi);
    const double psi_v = Lam_G * G_v + ell_v + dx_Mi * flux_v / flux_M;
    out.dG_ddir[i] = G_v - G_M * psi_v / psi_M;
  }
  return out;
}

// Below the drive floor the flux is the floor's scaled by h, but the end
// states were the floor's: M_in and M_exit at G0, up to 80x the flow's own at
// a sub-Pa drive (#1017), which the second-law heat band and the diagnostics
// read. Re-evaluate them at the blended flux: M_in from the isentropic inlet
// at G, M_exit at the floor's exit-to-inlet ratio (at low Mach M scales with
// G at a fixed stagnation state; the ratio moves by ~drive / P). G and its
// partials are untouched, so the residual and its Jacobian are too.
void blend_end_states(FannoChannelFlow& out, const FannoChannelFlow& at, double Pt0, double Tt0,
                      const std::vector<double>& X) {
  out.M_exit_root = at.M_exit;
  if (!(out.G > 0.0) || !(at.M_in > 0.0)) {
    out.M_in = 0.0;
    out.M_exit = 0.0;
    return;
  }
  const double Mi = fanno_inlet_mach(Pt0, Tt0, out.G, X);
  out.M_exit = std::min(1.0, at.M_exit * (Mi / at.M_in));
  out.M_in = Mi;
}

// Low-Mach viscous/inertial duct closure below kFannoFlowDriveFloor (#1216).
// Solves Delta P_inc(G) = drive directly in drive (no Pt0 - P cancellation),
// using friction_f_re_and_jacobian (no laminar re_floor), so G -> Hagen-
// Poiseuille (1 - w(0)) as drive -> 0 across all L/D and pressures.
// d(dp)/dv at fixed G of the low-Mach drop below, dp = F + A (or the
// subsonic-cap form F sqrt(F / (F - 2 A)) where A < 0), along a mass-fraction
// direction v, exact (#1311; d_ddir's difference until then). Composition
// enters only through R, in both prefactors (c_acc and fric_pref are each
// proportional to R), and through each node's mu: with Re = G D / mu,
// d(mu F_Re(Re))/dmu = F_Re - Re F_Re'.
double low_mach_drop_ddir(double G, double c_acc, double fric_pref, int n_nodes,
                          const double* w_gl, const double* Tt_n, const double* mu_n,
                          const std::string& model, double e_D, double D, double Pt0,
                          const std::vector<double>& X, const std::vector<double>& v) {
  const MixTangent<1> m = mix_tangent<1>(X, &v, 0);
  const double R_rel = -m.mw_v / m.mw;  // dR/dv / R
  double I_fric = 0.0;
  double I_fric_v = 0.0;
  for (int i = 0; i < n_nodes; ++i) {
    const std::vector<double> g = viscosity_dX(Tt_n[i], Pt0, X);
    double dmu_v = 0.0;
    for (std::size_t k = 0; k < g.size(); ++k) {
      dmu_v += g[k] * m.dX[k];
    }
    const double Re = G * D / mu_n[i];
    const auto [F_Re, dF_dRe] = solver::friction_f_re_and_jacobian(model, Re, e_D);
    I_fric += w_gl[i] * mu_n[i] * Tt_n[i] * F_Re;
    I_fric_v += w_gl[i] * Tt_n[i] * (F_Re - Re * dF_dRe) * dmu_v;
  }
  const double F = fric_pref * G * I_fric;
  const double F_v = F * R_rel + fric_pref * G * I_fric_v;
  const double A = c_acc * G * G;
  const double A_v = A * R_rel;
  if (A >= 0.0 || !(F > 0.0)) {
    return F_v + A_v;
  }
  const double den = F - 2.0 * A;
  const double sq = std::sqrt(F / den);
  return sq * (F - 3.0 * A) / den * F_v + sq * F / den * A_v;
}

struct LowMachSolution {
  bool valid = false;
  double G = 0.0;
  double dG_ddrive = 0.0;
  double dG_dPt0_fixed_drive = 0.0;
  double dG_dTt0 = 0.0;
  double dG_dTt_exit = 0.0;
  std::vector<double> dG_ddir;
};

LowMachSolution low_mach_duct_flow(double Pt0, double Tt0, double Tt_exit,
                                   const std::vector<double>& X, double drive, bool exit_total,
                                   double L, double D, double rough, const std::string& model,
                                   double fmult, const std::vector<std::vector<double>>& dirs,
                                   bool with_derivatives, bool tied_Tt) {
  LowMachSolution sol;
  if (model == "fixed" || !(fmult > 0.0)) {
    return sol;
  }
  const double e_D = rough / D;
  const auto [f_re_0, df_re_0] = solver::friction_f_re_and_jacobian(model, 0.0, e_D);
  (void)df_re_0;
  if (!(f_re_0 > 0.0)) {
    return sol;  // unknown friction tag
  }

  const double R = specific_gas_constant(X);
  const bool uniform_T = (Tt_exit == Tt0);
  const double c_acc = (R / Pt0) * ((exit_total ? 0.5 : 1.0) * Tt_exit - 0.5 * Tt0);
  const double fric_pref = fmult * L * R / (2.0 * Pt0 * D * D);

  // Precompute station temperatures and viscosities once (1 node at uniform
  // Tt_exit == Tt0 with xi = 0.5 matching the exact integral of xi and 1 - xi
  // over [0, 1]; 12 Gauss-Legendre nodes on [0, 1] otherwise).
  constexpr int kGlNodes = 12;
  const int n_nodes = uniform_T ? 1 : kGlNodes;
  std::array<double, kGlNodes> xi{};
  std::array<double, kGlNodes> w_gl{};
  std::array<double, kGlNodes> Tt_n{};
  std::array<double, kGlNodes> mu_n{};
  std::array<double, kGlNodes> dmu_n{};
  if (uniform_T) {
    xi[0] = 0.5;
    w_gl[0] = 1.0;
    Tt_n[0] = Tt0;
    if (with_derivatives) {
      const TransportDerivatives tr = transport_and_dT(Tt0, Pt0, X);
      mu_n[0] = tr.mu;
      dmu_n[0] = tr.dmu_dT;
    } else {
      mu_n[0] = viscosity(Tt0, Pt0, X);
    }
  } else {
    for (std::size_t j = 0; j < kGlX.size(); ++j) {
      const double x_minus = 0.5 * (1.0 - kGlX[j]);
      const double x_plus = 0.5 * (1.0 + kGlX[j]);
      const double wt = 0.5 * kGlW[j];
      const int i0 = static_cast<int>(2 * j);
      const int i1 = i0 + 1;
      xi[i0] = x_minus;
      w_gl[i0] = wt;
      xi[i1] = x_plus;
      w_gl[i1] = wt;
    }
    for (int i = 0; i < n_nodes; ++i) {
      Tt_n[i] = (1.0 - xi[i]) * Tt0 + xi[i] * Tt_exit;
      if (with_derivatives) {
        const TransportDerivatives tr = transport_and_dT(Tt_n[i], Pt0, X);
        mu_n[i] = tr.mu;
        dmu_n[i] = tr.dmu_dT;
      } else {
        mu_n[i] = viscosity(Tt_n[i], Pt0, X);
      }
    }
  }

  // Combine friction drop F >= 0 and acceleration A = c_acc * G^2.
  // For A >= 0 (adiabatic, heated, or mildly cooled static exit), dp = F + A.
  // For A < 0 (cooling with Rayleigh gain), use the C1 unconditionally-monotone
  // continuation dp = F * sqrt(F / (F - 2 A)): matches F + A in value, first
  // derivatives in F and A, and O(A^2 / F) across A = 0; homogeneous of degree
  // 1 in 1/Pt0; strictly increasing for all G >= 0 since 3 G F' - 2 F > 0
  // across all Re >= 0, so a unique positive root always exists with no fold.
  struct DpEval {
    double dp = 0.0;
    double ddp_dG = 0.0;
    double wF = 1.0;
    double wA = 1.0;
  };
  auto eval_dp = [&](double G) -> DpEval {
    double I_fric = 0.0;
    double I_fric_G = 0.0;
    for (int i = 0; i < n_nodes; ++i) {
      const double Re = G * D / mu_n[i];
      const auto [F_Re, dF_dRe] = solver::friction_f_re_and_jacobian(model, Re, e_D);
      const double coeff = w_gl[i] * mu_n[i] * Tt_n[i];
      I_fric += coeff * F_Re;
      I_fric_G += coeff * (F_Re + Re * dF_dRe);
    }
    const double F = fric_pref * G * I_fric;
    const double F_G = fric_pref * I_fric_G;
    const double A = c_acc * G * G;
    const double A_G = 2.0 * c_acc * G;
    DpEval ev;
    if (A >= 0.0 || !(F > 0.0)) {
      ev.dp = F + A;
      ev.ddp_dG = F_G + A_G;
      ev.wF = 1.0;
      ev.wA = 1.0;
    } else {
      const double den = F - 2.0 * A;
      const double s = std::sqrt(F / den);
      ev.dp = F * s;
      ev.wF = s * (F - 3.0 * A) / den;
      ev.wA = s * F / den;
      ev.ddp_dG = ev.wF * F_G + ev.wA * A_G;
    }
    return ev;
  };

  double G = 0.0;
  DpEval ev_root = eval_dp(0.0);
  if (!(ev_root.ddp_dG > 0.0)) {
    return sol;
  }

  if (drive > 0.0) {
    const double R0 = ev_root.ddp_dG;
    double lo = 0.0;
    double hi = 0.0;
    const double G0 = (c_acc > 0.0)
                          ? (2.0 * drive) / (R0 + std::sqrt(R0 * R0 + 4.0 * c_acc * drive))
                          : drive / R0;
    double g_cur = G0;
    DpEval cur = eval_dp(g_cur);
    if (cur.dp < drive) {
      lo = g_cur;
      hi = g_cur;
      for (int k = 0; k < 100 && cur.dp < drive; ++k) {
        lo = hi;
        hi *= 2.0;
        cur = eval_dp(hi);
      }
      if (cur.dp < drive) {
        return sol;
      }
    } else {
      hi = g_cur;
    }
    G = (lo > 0.0) ? 0.5 * (lo + hi) : g_cur;

    bool converged = false;
    for (int it = 0; it < 80; ++it) {
      const DpEval ev = eval_dp(G);
      const double f = ev.dp - drive;
      if (!(ev.ddp_dG > 0.0)) {
        return sol;
      }
      if (std::abs(f) <= 1e-15 * drive) {
        ev_root = ev;
        converged = true;
        break;
      }
      if (f > 0.0) {
        hi = G;
      } else {
        lo = G;
      }
      double G_new = G - f / ev.ddp_dG;
      if (std::abs(G_new - G) <= 1e-14 * std::abs(G_new) + 1e-300) {
        G = G_new;
        ev_root = eval_dp(G);
        converged = true;
        break;
      }
      if (!(G_new > lo && G_new < hi)) {
        G_new = 0.5 * (lo + hi);
      }
      G = G_new;
    }
    if (!converged || !(ev_root.ddp_dG > 0.0)) {
      return sol;
    }
  }

  const double ddp_dG = ev_root.ddp_dG;
  sol.valid = true;
  sol.G = G;
  sol.dG_ddrive = 1.0 / ddp_dG;
  if (!with_derivatives) {
    return sol;
  }

  // Both F and A (and F * sqrt(F / (F - 2 A))) are homogeneous of degree 1 in
  // 1/Pt0, so d(dp)/dPt0 at fixed G is always -drive / Pt0.
  sol.dG_dPt0_fixed_drive = (G > 0.0) ? (drive / (Pt0 * ddp_dG)) : 0.0;
  if (G > 0.0) {
    const double dA_dTt0 = (-0.5 * R / Pt0) * G * G;
    const double dA_dTt_exit = ((exit_total ? 0.5 : 1.0) * R / Pt0) * G * G;
    double dF_dTt0 = 0.0;
    double dF_dTt_exit = 0.0;
    for (int i = 0; i < n_nodes; ++i) {
      const double Re = G * D / mu_n[i];
      const auto [F_Re, dF_dRe] = solver::friction_f_re_and_jacobian(model, Re, e_D);
      const double dterm_dTt = mu_n[i] * F_Re + Tt_n[i] * (F_Re - Re * dF_dRe) * dmu_n[i];
      const double dnode = fric_pref * G * w_gl[i] * dterm_dTt;
      dF_dTt0 += dnode * (1.0 - xi[i]);
      dF_dTt_exit += dnode * xi[i];
    }
    const double ddp_dTt0 = ev_root.wA * dA_dTt0 + ev_root.wF * dF_dTt0;
    const double ddp_dTt_exit = ev_root.wA * dA_dTt_exit + ev_root.wF * dF_dTt_exit;
    if (tied_Tt) {
      // Combine both ends into Tt0 only for fanno_channel_flow (Tt_exit tied
      // to Tt0); fanno_rayleigh_channel_flow treats Tt0 and Tt_exit as
      // independent arguments even when Tt_exit == Tt0.
      sol.dG_dTt0 = -(ddp_dTt0 + ddp_dTt_exit) / ddp_dG;
      sol.dG_dTt_exit = 0.0;
    } else {
      sol.dG_dTt0 = -ddp_dTt0 / ddp_dG;
      sol.dG_dTt_exit = -ddp_dTt_exit / ddp_dG;
    }
  }

  if (!dirs.empty()) {
    sol.dG_ddir.assign(dirs.size(), 0.0);
    if (G > 0.0) {
      for (std::size_t i = 0; i < dirs.size(); ++i) {
        sol.dG_ddir[i] = -low_mach_drop_ddir(G, c_acc, fric_pref, n_nodes, w_gl.data(),
                                             Tt_n.data(), mu_n.data(), model, e_D, D, Pt0, X,
                                             dirs[i]) /
                         ddp_dG;
      }
    }
  }
  return sol;
}

// Floor log-slope derivatives d(sigma)/dtheta at drive = kFannoFlowDriveFloor
// (where sigma = clamp(-at.dG_dP_target * floor / at.G, 0, kDriveFloorSigmaMax)).
//
// Design rationale (owner-approved FD exception in #1303; follow-up #1317):
//   - We keep the exact drive-derivative sigma rather than a value-based secant
//     or s = 1 blend so dG/dP_target is strictly C1 at drive = 1 Pa with a
//     single compressible solve on residual-only calls (a secant above 1 Pa
//     introduces a slope jump at 1 Pa or re-triggers Pt0 - P cancellation and
//     doubles residual cost; s = 1 breaks C1 on choked and cooled ducts).
//   - Because sigma is built from a first derivative of G, d(sigma)/dtheta is
//     a mixed second partial and the compressible kernel has no second-order
//     path yet. Internal C++ FD is therefore used here as an owner-approved
//     exception until #1317 replaces it with exact nested-dual second-order
//     forward mode on top of #1304's DualN<4> tangents:
//       * Pt0, Tt0, Tt_exit: Schwarz mixed-partial identity
//         d/dtheta (dG/ddrive) = d/ddrive (dG/dtheta) via 4-point O(h^4)
//         central difference in drive around floor (falling back to 2-point
//         central, then 1-sided O(h) near choke/fold boundaries), using
//         warm-started kernel evaluations with dirs empty.
//       * Composition directions: direct d_ddir(X, dirs, sigma(X)) at step
//         kFannoFloorSigmaDirScale * kFannoDirStep (1e-3), avoiding the
//         low-Mach d_ddir cancellation jitter that Schwarz differentiation
//         would amplify.
struct SigmaDerivs {
  double dPt0_fixed = 0.0;
  double dTt0 = 0.0;
  double dTt_exit = 0.0;
  std::vector<double> ddir;
};

template <typename EvalAtDrive, typename EvalAtX>
SigmaDerivs compute_floor_sigma_derivs(const FannoChannelFlow& at, const EvalAtDrive& eval_at_drive,
                                       const EvalAtX& eval_at_X, const std::vector<double>& X,
                                       const std::vector<std::vector<double>>& dirs) {
  SigmaDerivs sd;
  sd.ddir.assign(dirs.size(), 0.0);
  if (!(at.G > 0.0) || at.choked) {
    return sd;
  }
  const double floor = kFannoFlowDriveFloor;
  const double sigma_raw = -at.dG_dP_target * floor / at.G;
  if (!(sigma_raw > 0.0 && sigma_raw < kDriveFloorSigmaMax)) {
    return sd;
  }
  const double hp = kFannoFloorSigmaDriveStep * floor;
  auto valid_nbr = [&](const FannoChannelFlow& r) {
    return !r.choked && r.G > 0.0 && std::abs(r.G - at.G) < 0.25 * at.G;
  };
  const FannoChannelFlow at_p1 = eval_at_drive(floor + hp);
  const FannoChannelFlow at_m1 = eval_at_drive(floor - hp);
  const bool ok_p1 = valid_nbr(at_p1);
  const bool ok_m1 = valid_nbr(at_m1);
  FannoChannelFlow at_p2;
  FannoChannelFlow at_m2;
  bool ok_4pt = false;
  if (ok_p1 && ok_m1) {
    at_p2 = eval_at_drive(floor + 2.0 * hp);
    at_m2 = eval_at_drive(floor - 2.0 * hp);
    ok_4pt = valid_nbr(at_p2) && valid_nbr(at_m2);
  }
  if (ok_p1 || ok_m1) {
    auto dsig = [&](double d0, double dp1, double dm1, double dp2, double dm2) {
      double dg1_dth = 0.0;
      if (ok_4pt) {
        dg1_dth = (8.0 * (dp1 - dm1) - (dp2 - dm2)) / (12.0 * hp);
      } else if (ok_p1 && ok_m1) {
        dg1_dth = (dp1 - dm1) / (2.0 * hp);
      } else if (ok_p1) {
        dg1_dth = (dp1 - d0) / hp;
      } else {
        dg1_dth = (d0 - dm1) / hp;
      }
      return (floor / at.G) * dg1_dth - sigma_raw * (d0 / at.G);
    };
    sd.dPt0_fixed = dsig(at.dG_dPt0 + at.dG_dP_target, at_p1.dG_dPt0 + at_p1.dG_dP_target,
                         at_m1.dG_dPt0 + at_m1.dG_dP_target, at_p2.dG_dPt0 + at_p2.dG_dP_target,
                         at_m2.dG_dPt0 + at_m2.dG_dP_target);
    sd.dTt0 = dsig(at.dG_dTt0, at_p1.dG_dTt0, at_m1.dG_dTt0, at_p2.dG_dTt0, at_m2.dG_dTt0);
    sd.dTt_exit =
        dsig(at.dG_dTt_exit, at_p1.dG_dTt_exit, at_m1.dG_dTt_exit, at_p2.dG_dTt_exit,
             at_m2.dG_dTt_exit);
  }
  if (!dirs.empty()) {
    sd.ddir = d_ddir(
        X, dirs,
        [&](const std::vector<double>& Xk) {
          const FannoChannelFlow r = eval_at_X(Xk);
          if (r.choked || !(r.G > 0.0)) {
            return sigma_raw;
          }
          return std::clamp(-r.dG_dP_target * floor / r.G, 0.0, kDriveFloorSigmaMax);
        },
        kFannoFloorSigmaDirScale);
  }
  return sd;
}

// Evaluate dp_inc(G) and its analytical partials at a given mass flux G (used
// to compute chi = dp_inc(at.G) / kFannoFlowDriveFloor and d(chi)/dtheta for
// the sub-floor choking smoothstep).
struct LowMachDpAtG {
  bool valid = false;
  double dp = 0.0;
  double ddp_dG = 0.0;
  double ddp_dTt0 = 0.0;
  double ddp_dTt_exit = 0.0;
  std::vector<double> ddp_ddir;
};

LowMachDpAtG low_mach_dp_at_G(double Pt0, double Tt0, double Tt_exit, const std::vector<double>& X,
                              double G, bool exit_total, double L, double D, double rough,
                              const std::string& model, double fmult,
                              const std::vector<std::vector<double>>& dirs, bool with_derivatives,
                              bool tied_Tt) {
  LowMachDpAtG out;
  if (model == "fixed" || !(fmult > 0.0) || !(G > 0.0)) {
    return out;
  }
  const double e_D = rough / D;
  const auto [f_re_0, df_re_0] = solver::friction_f_re_and_jacobian(model, 0.0, e_D);
  (void)df_re_0;
  if (!(f_re_0 > 0.0)) {
    return out;
  }
  const double R = specific_gas_constant(X);
  const bool uniform_T = (Tt_exit == Tt0);
  const double c_acc = (R / Pt0) * ((exit_total ? 0.5 : 1.0) * Tt_exit - 0.5 * Tt0);
  const double fric_pref = fmult * L * R / (2.0 * Pt0 * D * D);

  constexpr int kGlNodes = 12;
  const int n_nodes = uniform_T ? 1 : kGlNodes;
  std::array<double, kGlNodes> xi{};
  std::array<double, kGlNodes> w_gl{};
  std::array<double, kGlNodes> Tt_n{};
  std::array<double, kGlNodes> mu_n{};
  std::array<double, kGlNodes> dmu_n{};
  if (uniform_T) {
    xi[0] = 0.5;
    w_gl[0] = 1.0;
    Tt_n[0] = Tt0;
    if (with_derivatives) {
      const TransportDerivatives tr = transport_and_dT(Tt0, Pt0, X);
      mu_n[0] = tr.mu;
      dmu_n[0] = tr.dmu_dT;
    } else {
      mu_n[0] = viscosity(Tt0, Pt0, X);
    }
  } else {
    for (std::size_t j = 0; j < kGlX.size(); ++j) {
      const double x_minus = 0.5 * (1.0 - kGlX[j]);
      const double x_plus = 0.5 * (1.0 + kGlX[j]);
      const double wt = 0.5 * kGlW[j];
      const int i0 = static_cast<int>(2 * j);
      const int i1 = i0 + 1;
      xi[i0] = x_minus;
      w_gl[i0] = wt;
      xi[i1] = x_plus;
      w_gl[i1] = wt;
    }
    for (int i = 0; i < n_nodes; ++i) {
      Tt_n[i] = (1.0 - xi[i]) * Tt0 + xi[i] * Tt_exit;
      if (with_derivatives) {
        const TransportDerivatives tr = transport_and_dT(Tt_n[i], Pt0, X);
        mu_n[i] = tr.mu;
        dmu_n[i] = tr.dmu_dT;
      } else {
        mu_n[i] = viscosity(Tt_n[i], Pt0, X);
      }
    }
  }

  double I_fric = 0.0;
  double I_fric_G = 0.0;
  for (int i = 0; i < n_nodes; ++i) {
    const double Re = G * D / mu_n[i];
    const auto [F_Re, dF_dRe] = solver::friction_f_re_and_jacobian(model, Re, e_D);
    const double coeff = w_gl[i] * mu_n[i] * Tt_n[i];
    I_fric += coeff * F_Re;
    I_fric_G += coeff * (F_Re + Re * dF_dRe);
  }
  const double F = fric_pref * G * I_fric;
  const double F_G = fric_pref * I_fric_G;
  const double A = c_acc * G * G;
  const double A_G = 2.0 * c_acc * G;
  double wF = 1.0;
  double wA = 1.0;
  if (A >= 0.0 || !(F > 0.0)) {
    out.dp = F + A;
    out.ddp_dG = F_G + A_G;
  } else {
    const double den = F - 2.0 * A;
    const double s = std::sqrt(F / den);
    out.dp = F * s;
    wF = s * (F - 3.0 * A) / den;
    wA = s * F / den;
    out.ddp_dG = wF * F_G + wA * A_G;
  }
  out.valid = true;
  if (!with_derivatives) {
    return out;
  }
  const double dA_dTt0 = (-0.5 * R / Pt0) * G * G;
  const double dA_dTt_exit = ((exit_total ? 0.5 : 1.0) * R / Pt0) * G * G;
  double dF_dTt0 = 0.0;
  double dF_dTt_exit = 0.0;
  for (int i = 0; i < n_nodes; ++i) {
    const double Re = G * D / mu_n[i];
    const auto [F_Re, dF_dRe] = solver::friction_f_re_and_jacobian(model, Re, e_D);
    const double dterm_dTt = mu_n[i] * F_Re + Tt_n[i] * (F_Re - Re * dF_dRe) * dmu_n[i];
    const double dnode = fric_pref * G * w_gl[i] * dterm_dTt;
    dF_dTt0 += dnode * (1.0 - xi[i]);
    dF_dTt_exit += dnode * xi[i];
  }
  const double ddp_dTt0 = wA * dA_dTt0 + wF * dF_dTt0;
  const double ddp_dTt_exit = wA * dA_dTt_exit + wF * dF_dTt_exit;
  if (tied_Tt) {
    out.ddp_dTt0 = ddp_dTt0 + ddp_dTt_exit;
    out.ddp_dTt_exit = 0.0;
  } else {
    out.ddp_dTt0 = ddp_dTt0;
    out.ddp_dTt_exit = ddp_dTt_exit;
  }
  if (!dirs.empty()) {
    out.ddp_ddir.assign(dirs.size(), 0.0);
    for (std::size_t i = 0; i < dirs.size(); ++i) {
      out.ddp_ddir[i] = low_mach_drop_ddir(G, c_acc, fric_pref, n_nodes, w_gl.data(),
                                           Tt_n.data(), mu_n.data(), model, e_D, D, Pt0, X,
                                           dirs[i]);
    }
  }
  return out;
}

// Sub-floor blend for both fanno_channel_flow (tied_Tt = true) and
// fanno_rayleigh_channel_flow (tied_Tt = false) below kFannoFlowDriveFloor
// (#1216).
//
// For every registered viscous friction model with f_multiplier > 0, maps the
// normalized low-Mach flux u = low.G / low_1.G in [0, 1] through the
// shape-preserving C1 rational cubic:
//   P(u) = u^3 + (r + s) u^2 v + r u v^2,   v = 1 - u,
//   Q(u) = P(u) + v (v + s u),
//   y(u) = P(u) / Q(u),
//   G_visc = at.G * y(u),
//   r     = low_1.G / at.G > 0,
//   s     = sigma / sigma_low >= 0,
//   sigma = clamp(-at.dG_dP_target * floor / at.G, 0, kDriveFloorSigmaMax),
//   sigma_low = low_1.dG_ddrive * floor / low_1.G > 0.
// Properties (verified symbolically):
//   - y(0) = 0, y'(0) = r  => G_visc(drive) = low.G + O(drive^2) (exact Hagen-
//     Poiseuille as drive -> 0, even when the 1 Pa floor state chokes).
//   - y(1) = 1, y'(1) = s, dy/dr(1) = dy/ds(1) = 0  => G_visc and all five
//     partials match the compressible kernel continuously at drive = floor.
//   - dy/du = (r v^4 + 2(r+s) u v^3 + (r+s^2+s+3) u^2 v^2 + 2(s+1) u^3 v
//             + s u^4) / Q^2 > 0 for all u in [0, 1), r > 0, s >= 0
//     => strictly monotone with 0 <= G_visc <= at.G unconditionally.
//   - dy/dr = u v (Q - P) / Q^2,  dy/ds = u v (u (Q - P) - P) / Q^2.
//   - Retaining ds/dtheta via Schwartz mixed partials (and dh/dsigma on the
//     odd-cubic branch) keeps all five partials exact in the blend interior
//     (0, floor), including cooled exit_total Rayleigh ducts.
//
// When vanishing length (L -> 0) or f_multiplier -> 0 chokes the duct well
// below kFannoFlowDriveFloor, chi = dp_inc(at.G) / floor drops below the exact
// Fanno sonic limit (2 / (gamma + 1))^((gamma + 1) / (gamma - 1)) in
// [0.316, 0.368] (0.3349 for gamma = 1.4). Across the binary octave
// chi in [kFannoLowMachDropLo = 0.125, kFannoLowMachDropHi = 0.25], a C2
// quintic smoothstep w(chi) blends G_visc smoothly into the odd cubic
// G_cub = at.G * h(t, sigma).
FannoChannelFlow apply_sub_floor_blend(const FannoChannelFlow& at, const SigmaDerivs& sig_d,
                                       double Pt0, double Tt0, double Tt_exit,
                                       const std::vector<double>& X, double drive, bool exit_total,
                                       double L, double D, double rough, const std::string& model,
                                       double fmult, const std::vector<std::vector<double>>& dirs,
                                       bool with_derivatives, bool tied_Tt) {
  const double floor = kFannoFlowDriveFloor;
  const double t = drive / floor;
  const double sigma =
      (at.G > 0.0) ? std::clamp(-at.dG_dP_target * floor / at.G, 0.0, kDriveFloorSigmaMax) : 0.5;
  const DriveFloorCubic cubic = drive_floor_cubic(t, sigma);
  const double h = cubic.h;
  const double dh = cubic.dh;
  const double dh_dsigma = 0.5 * t * (t * t - 1.0);
  const bool have_sig_d = with_derivatives && drive > 0.0 && sig_d.ddir.size() == dirs.size();

  if (at.G > 0.0 && model != "fixed" && fmult > 0.0) {
    const LowMachDpAtG dp_at = low_mach_dp_at_G(Pt0, Tt0, Tt_exit, X, at.G, exit_total, L, D, rough,
                                                model, fmult, dirs, with_derivatives, tied_Tt);
    const double chi = dp_at.valid ? (dp_at.dp / floor) : 0.0;
    if (dp_at.valid && chi > kFannoLowMachDropLo) {
      const LowMachSolution low_1 =
          low_mach_duct_flow(Pt0, Tt0, Tt_exit, X, floor, exit_total, L, D, rough, model, fmult,
                             dirs, with_derivatives, tied_Tt);
      if (low_1.valid && low_1.G > 0.0 && low_1.dG_ddrive > 0.0) {
        const LowMachSolution low =
            low_mach_duct_flow(Pt0, Tt0, Tt_exit, X, drive, exit_total, L, D, rough, model, fmult,
                               dirs, with_derivatives, tied_Tt);
        if (low.valid) {
          double w = 1.0;
          double dw_dchi = 0.0;
          if (chi < kFannoLowMachDropHi) {
            const double span = kFannoLowMachDropHi - kFannoLowMachDropLo;
            const double xi = (chi - kFannoLowMachDropLo) / span;
            const double xi2 = xi * xi;
            const double one_m_xi = 1.0 - xi;
            w = xi2 * xi * (10.0 + xi * (-15.0 + 6.0 * xi));
            dw_dchi = 30.0 * xi2 * one_m_xi * one_m_xi / span;
          }

          const double r = low_1.G / at.G;
          const double u = std::clamp(low.G / low_1.G, 0.0, 1.0);
          const double v = 1.0 - u;
          const double sigma_low = low_1.dG_ddrive * floor / low_1.G;
          const double s = sigma / sigma_low;
          const double u2 = u * u;
          const double v2 = v * v;
          const double P = u * (u2 + (r + s) * u * v + r * v2);
          const double Q_m_P = v * (v + s * u);
          const double Q = P + Q_m_P;
          const double inv_Q = 1.0 / Q;
          const double y = P * inv_Q;

          const double G_visc = at.G * y;
          const double G_cub = at.G * h;

          FannoChannelFlow out = at;
          out.choked = false;
          out.G = w * G_visc + (1.0 - w) * G_cub;
          blend_end_states(out, at, Pt0, Tt0, X);
          if (drive == 0.0) {
            out.M_exit_root = 0.0;
          }
          if (!with_derivatives) {
            out.dG_dPt0 = 0.0;
            out.dG_dTt0 = 0.0;
            out.dG_dP_target = 0.0;
            out.dG_dTt_exit = 0.0;
            out.dG_ddir.clear();
            return out;
          }
          const double num_du = r * v2 * v2 + 2.0 * (r + s) * u * v2 * v +
                                (r + s * s + s + 3.0) * u2 * v2 + 2.0 * (s + 1.0) * u2 * u * v +
                                s * u2 * u2;
          const double dy_du = (num_du * inv_Q) * inv_Q;
          const double dy_dr = ((u * v * Q_m_P) * inv_Q) * inv_Q;
          const double dy_ds = ((u * v * (u * Q_m_P - P)) * inv_Q) * inv_Q;
          const double dG_visc_ddrive = at.G * dy_du * (low.dG_ddrive / low_1.G);
          const double dG_cub_ddrive = at.G * dh / floor;
          const double dG_ddrive = w * dG_visc_ddrive + (1.0 - w) * dG_cub_ddrive;

          // Compute d(sigma_low)/dtheta via 2-point central Schwarz mixed
          // partials on low_mach_duct_flow at floor +- 1e-4 Pa (only needed
          // when dy_ds != 0 and sigma > 0; part of the owner-approved #1303
          // internal-FD exception, tracked for exact second-order replacement
          // in #1317).
          SigmaDerivs sig_low_d;
          sig_low_d.ddir.assign(dirs.size(), 0.0);
          if (drive > 0.0 && sigma > 0.0) {
            const double hp_low = 1e-4 * floor;
            const LowMachSolution lp = low_mach_duct_flow(Pt0, Tt0, Tt_exit, X, floor + hp_low,
                                                          exit_total, L, D, rough, model, fmult,
                                                          dirs, true, tied_Tt);
            const LowMachSolution lm = low_mach_duct_flow(Pt0, Tt0, Tt_exit, X, floor - hp_low,
                                                          exit_total, L, D, rough, model, fmult,
                                                          dirs, true, tied_Tt);
            if (lp.valid && lm.valid) {
              auto dsig_l = [&](double d0, double dp_v, double dm_v) {
                const double dg1_dth = (dp_v - dm_v) / (2.0 * hp_low);
                return (floor / low_1.G) * dg1_dth - sigma_low * (d0 / low_1.G);
              };
              sig_low_d.dPt0_fixed = dsig_l(low_1.dG_dPt0_fixed_drive, lp.dG_dPt0_fixed_drive,
                                            lm.dG_dPt0_fixed_drive);
              sig_low_d.dTt0 = dsig_l(low_1.dG_dTt0, lp.dG_dTt0, lm.dG_dTt0);
              sig_low_d.dTt_exit = dsig_l(low_1.dG_dTt_exit, lp.dG_dTt_exit, lm.dG_dTt_exit);
              for (std::size_t i = 0; i < dirs.size(); ++i) {
                sig_low_d.ddir[i] = dsig_l(low_1.dG_ddir[i], lp.dG_ddir[i], lm.dG_ddir[i]);
              }
            }
          }

          auto chain_param = [&](double dG_inc, double dG_inc_1, double dG_comp_1, double dsig_c,
                                 double dsig_l, double ddp_at_dtheta, bool is_Pt0) {
            const double du = (dG_inc - u * dG_inc_1) / low_1.G;
            const double dr = (dG_inc_1 - r * dG_comp_1) / at.G;
            const double ds = (dsig_c - s * dsig_l) / sigma_low;
            const double dG_visc =
                dG_comp_1 * y + at.G * (dy_du * du + dy_dr * dr + dy_ds * ds);
            const double dG_cub = h * dG_comp_1 + at.G * dh_dsigma * dsig_c;
            double dchi = 0.0;
            if (dw_dchi != 0.0) {
              dchi = (dp_at.ddp_dG * dG_comp_1 + ddp_at_dtheta) / floor -
                     (is_Pt0 ? (chi / Pt0) : 0.0);
            }
            return w * dG_visc + (1.0 - w) * dG_cub + dw_dchi * dchi * (G_visc - G_cub);
          };
          out.dG_dP_target = -dG_ddrive;
          const double dG_dPt0_fixed = chain_param(
              low.dG_dPt0_fixed_drive, low_1.dG_dPt0_fixed_drive, at.dG_dPt0 + at.dG_dP_target,
              have_sig_d ? sig_d.dPt0_fixed : 0.0, sig_low_d.dPt0_fixed, 0.0, true);
          out.dG_dPt0 = dG_ddrive + dG_dPt0_fixed;
          out.dG_dTt0 = chain_param(low.dG_dTt0, low_1.dG_dTt0, at.dG_dTt0,
                                    have_sig_d ? sig_d.dTt0 : 0.0, sig_low_d.dTt0, dp_at.ddp_dTt0,
                                    false);
          out.dG_dTt_exit =
              chain_param(low.dG_dTt_exit, low_1.dG_dTt_exit, at.dG_dTt_exit,
                          have_sig_d ? sig_d.dTt_exit : 0.0, sig_low_d.dTt_exit, dp_at.ddp_dTt_exit,
                          false);
          out.dG_ddir.assign(dirs.size(), 0.0);
          for (std::size_t i = 0; i < dirs.size(); ++i) {
            out.dG_ddir[i] = chain_param(low.dG_ddir[i], low_1.dG_ddir[i], at.dG_ddir[i],
                                         have_sig_d ? sig_d.ddir[i] : 0.0, sig_low_d.ddir[i],
                                         dp_at.ddp_ddir[i], false);
          }
          return out;
        }
      }
    }
  }

  // Fallback: odd cubic matching value and slope at the floor (fixed/
  // frictionless ducts or ducts already choked deep below the 1 Pa floor
  // with dp_inc(at.G) / floor <= kFannoLowMachDropLo).
  FannoChannelFlow out = at;
  out.choked = false;
  out.G = at.G * h;
  blend_end_states(out, at, Pt0, Tt0, X);
  if (drive == 0.0) {
    out.M_exit_root = 0.0;  // no root at zero flow: no warm start, as before #1209
  }
  out.dG_dP_target = -at.G * dh / floor;
  out.dG_dPt0 = at.G * dh / floor + h * (at.dG_dPt0 + at.dG_dP_target) +
                (have_sig_d ? at.G * dh_dsigma * sig_d.dPt0_fixed : 0.0);
  out.dG_dTt0 = h * at.dG_dTt0 + (have_sig_d ? at.G * dh_dsigma * sig_d.dTt0 : 0.0);
  out.dG_dTt_exit = h * at.dG_dTt_exit + (have_sig_d ? at.G * dh_dsigma * sig_d.dTt_exit : 0.0);
  for (std::size_t i = 0; i < out.dG_ddir.size(); ++i) {
    out.dG_ddir[i] =
        h * out.dG_ddir[i] + (have_sig_d ? at.G * dh_dsigma * sig_d.ddir[i] : 0.0);
  }
  return out;
}

}  // namespace

FannoChannelFlow fanno_channel_flow(double Pt0, double Tt0, const std::vector<double>& X,
                                    double P_target, bool exit_total, double L, double D,
                                    double roughness, const std::string& friction_model,
                                    double f_multiplier, double M_exit_guess,
                                    const std::vector<std::vector<double>>& dY_directions,
                                    bool with_derivatives) {
  check_duct(L, D);
  check_friction(roughness, f_multiplier);
  FannoChannelFlow out;
  const double drive = Pt0 - P_target;
  static const std::vector<std::vector<double>> kNoDirections;
  const std::vector<std::vector<double>>& dirs = with_derivatives ? dY_directions : kNoDirections;
  // The no-flow branch is for finite states; NaN used to answer G = 0 (#769).
  if (!validation::is_finite(Pt0) || !validation::is_finite(Tt0) ||
      !validation::is_finite(P_target)) {
    throw std::invalid_argument("fanno_channel_flow: Pt0, Tt0 and P_target must be finite");
  }
  // No flow at or below zero drive (the caller reverses the duct for reversed
  // flow). At exactly zero drive the flux's slope is finite (#1209, #1216):
  // only the derivatives need the floor's flow.
  if (!(drive >= 0.0) || !(Pt0 > 0.0) || !(Tt0 > 0.0) || (drive == 0.0 && !with_derivatives)) {
    out.dG_ddir.assign(dirs.size(), 0.0);
    return out;
  }
  FlowProblem pb{Pt0, Tt0, P_target, L, D, roughness, f_multiplier, &X, &friction_model, exit_total,
                 M_exit_guess, &dirs, with_derivatives};
  if (drive >= kFannoFlowDriveFloor) {
    FannoChannelFlow raw = raw_flow(pb);
    raw.M_exit_root = raw.M_exit;
    return raw;
  }
  // Below the floor: evaluate the compressible kernel at kFannoFlowDriveFloor
  // (with derivatives for the C1 handoff slope) and blend smoothly to zero.
  pb.P = Pt0 - kFannoFlowDriveFloor;
  pb.with_derivatives = true;
  const FannoChannelFlow at = raw_flow(pb);
  SigmaDerivs sig_d;
  if (with_derivatives && drive > 0.0) {
    sig_d = compute_floor_sigma_derivs(
        at,
        [&](double d_eval) {
          FlowProblem pb_d = pb;
          pb_d.P = Pt0 - d_eval;
          pb_d.Me_guess = at.M_exit;
          pb_d.dirs = &kNoDirections;
          return raw_flow(pb_d);
        },
        [&](const std::vector<double>& Xk) {
          FlowProblem pb_x = pb;
          pb_x.X = &Xk;
          pb_x.Me_guess = at.M_exit;
          pb_x.dirs = &kNoDirections;
          return raw_flow(pb_x);
        },
        X, dirs);
  }
  return apply_sub_floor_blend(at, sig_d, Pt0, Tt0, Tt0, X, drive, exit_total, L, D, roughness,
                               friction_model, f_multiplier, dirs, with_derivatives, true);
}

// ---------------------------------------------------------------------------
// Fanno-Rayleigh (#551). See fanno_mach.h for the relations.
// ---------------------------------------------------------------------------

namespace {

// Fixed RK4 march on a smooth mesh graded toward the exit (#1337):
//   x_i = L (1 - t_i^p * exp(-alpha * (1 - t_i)^p)),  t_i = 1 - i / N,
//   alpha_back = ln(1 + lam),  lam = f_ref L / D,
//   alpha_fwd  = ln(1 + lam_eff),  lam_eff = lam / (1 + w_in / lam).
// A FIXED step count and C-infinity map, so the march is smooth in all
// arguments and Dual4 differentiates the grid and ODE stages exactly.
//   - Backward sonic march (back_w_inlet): exit is sonic (w_exit = 0), so
//     alpha = ln(1 + lam) scales exit steps as (D / f) / N^p to resolve the
//     O(D / f) sonic acceleration layer on arbitrarily long laminar ducts.
//   - Forward march (march_zeta): smooth rational gate lam / (lam + w_in) on
//     lam. Near choke (w_in ~ gamma lam) alpha_fwd ~ ln(1 + lam / (1 + gamma))
//     resolves the sonic layer; at low Mach (w_in >> lam) lam_eff -> 0 so nodes
//     stay spread across the duct instead of crowding the exit on heated/cooled
//     low-Re ducts.
//
// The state marched is the Fanno-regular coordinate w = 1/M^2 - 1 + ln(M^2) >= 0:
// near M = 1, w ~ (1 - M^2)^2 / 2 = zeta / 2 (regular at sonic); at low Mach
// w ~ 1/M^2, linear in x for Fanno friction across arbitrary f L / D, carrying
// full relative precision at low Mach and removing the theta Riccati blow-up.
constexpr int kRayleighSteps = 32;
constexpr double kRayleighGrading = 3.0;
// Width in log G to which the cold scan locates the largest flux the march
// carries subsonically (rayleigh_raw_flow, #1161): well inside the ~3e-7 by
// which the march's own sonic point falls short of the choke relation's flux.
constexpr double kRayleighSubsonicEndTol = 1e-12;
// A flux counts as carried subsonically while the march's zeta_min exceeds
// this (rayleigh_raw_flow, #1161). On a static-pressure exit the fluxes this
// excludes (within ~1e-10 of the march's sonic flux) span 0.1-8 Pa of drive
// below the choke (1e-6 to 5e-6 of it), now reported choked; on a
// total-pressure exit at most 1e-5 Pa. G there is within 3e-7 of Gc either way.
constexpr double kRayleighZetaFloor = 1e-10;
// A cooled duct's excess can fall through zero and rise again inside one
// ln-2 bracket of the cold scan (#1167; see rayleigh_raw_flow). Its slope at
// a bracket end is a backward difference over this step in log G (0.1% of
// the flux: far above the march's rounding, small against a bracket); a
// minimum is searched to kRayleighDipTol (plus 1.5e-8 relative) in its
// coordinate. The last bracket, which ends at the march's sonic flux, is
// sampled at fluxes 4^-j of that end below it, j = First..Last (6% down to
// 2e-10 of the flux).
constexpr double kRayleighSlopeStep = 1e-3;
constexpr double kRayleighDipTol = 1e-10;
constexpr int kRayleighEndGridFirst = 2;
constexpr int kRayleighEndGridLast = 16;
// A cooled duct's warm bracket widens (by 8x, from 1e-3) only while it spans
// less than this fraction of the guess: so at most to 6.4%.
constexpr double kRayleighWarmTurnSpan = 0.05;
// The march's exit pressure rounds. Over 21 fluxes 1e-9 apart in log G (1e-4
// of the distance to the choke nearer it), its scatter about a quartic, the
// same at a tenth of that spacing (so not truncation), is at most 7e-14 Pt0
// plus 6.2e-14 of the local slope |dE/dlog G| (7749 samples of the rsweep
// and #1225 review generators, 1e-4 Gc to 3e-8 below the choke): 2e-14 Pt0
// typical and up to 6.6e-14 Pt0 below 0.5 Gc and on the fold dips of #1221,
// growing with the slope toward the choke (to 1.7e-9 Pt0 within 1e-6 of it).
// An excess within kRayleighExcessNoise Pt0 of zero has no sign to read, so a
// cooled duct's warm bracket is taken only where both its ends are clear of
// it (#1221). The slope term exceeds the constant only within ~1e-13 in log G
// of a crossing; exact placements there, near the choke, found no miss.
constexpr double kRayleighExcessNoise = 1e-13;
// The Hermite test's derivative quadratic is treated as linear when its t^2
// coefficient is below this fraction of the others (rounding, not a turn).
constexpr double kRayleighHermiteDegenerate = 1e-14;
// The minimiser starts at the cubic's minimum, kept this fraction of the
// window away from its ends.
constexpr double kRayleighDipStartMargin = 0.05;
// Iteration cap of the minimiser (it converges in ~20-40).
constexpr int kRayleighDipMaxIter = 200;
// A cooled bracket whose ends straddle zero is walked in this many equal
// parts of log G, each given the Hermite test (#1190). The rise of a laminar
// fold spans 0.18-0.21 of log G on the measured ducts, against the bracket's
// ln 2. In the kernel, over 23020 random cooled cold calls (seeds 7, 23, 41)
// and the #1190 review's 23k (random, phase-swept and small-diameter ducts),
// the ends' slopes alone (one part) left 14 + 20 straddling folds on a later
// crossing; two parts found the first fall on all of them. Three, four and
// twelve parts found no more (only root rounding differs), at +0.6, +1.1-1.3
// and +5-6 marches per cold call over two.
constexpr int kRayleighStraddleSplit = 2;
// A cooled interior bracket whose ends are both positive and whose Hermite
// test shows no turn is walked in the same equal parts (#1206): a dip of the
// excess below zero narrower than the bracket, near one end of it, leaves the
// ends' values and slopes looking monotone. The gate is a cost filter, not a
// bound: the walk is skipped while |N| (N = drive - F, the needed drive)
// stays below kRayleighBothGate of the drive at both ends. Where N is
// non-negative and rises at most as G^2 (turbulent friction) over a ln-2
// bracket it grows at most 4x, so it cannot reach the drive inside one from
// below a quarter of it. A duct with a second fold after a negative valley
// breaks that (N rises from below zero through the dip): there the walk is
// skipped and the whole bracket's Hermite test is what catches the dip.
// Measured on the dense needed-drive grids of 6000 rsweep ducts (about 70 %
// cooled): of the 27173 (fold, drive) pairs whose dip is narrower than a
// bracket, every bracket that holds the dip has max |N| at its ends >= 0.74
// of the drive. The #1225 review's strongly cooled generator has 101 pairs
// below 0.25 (down to 5e-4 of the drive) in 21 of 800 ducts; their dips are
// wide enough for the whole-bracket test, and gate 0 gives the same results
// there (0 misses over 31794 warm calls). Over 9200 cooled cold and
// 149k warm calls of those ducts (and 6400 more calls placed in their
// fold bands) the walk costs +4.3-4.7 marches per cooled cold call (+8-9
// %, from 51-53); walking every both-positive bracket (gate 0, +19-20) or in
// four parts (+10-11) finds no more first falls.
constexpr double kRayleighBothGate = 0.25;

// Brent's minimiser (golden section with parabolic steps) of fn on (a, b)
// from x0 inside it, stopped as soon as fn reads <= 0: returns that point,
// or NaN if the minimum it converges to (to abs_tol + 1.5e-8 |x|) stays
// positive.
template <typename Fn>
double brent_min_to_zero(const Fn& fn, double a, double b, double x0, double abs_tol) {
  constexpr double kCGold = 0.3819660112501051;  // (3 - sqrt(5)) / 2
  constexpr double kSqrtEps = 1.5e-8;
  double x = (x0 > a && x0 < b) ? x0 : a + kCGold * (b - a);
  double fx = fn(x);
  if (fx <= 0.0) {
    return x;
  }
  double w = x, v = x, fw = fx, fv = fx;
  double d = 0.0, e = 0.0;
  for (int it = 0; it < kRayleighDipMaxIter; ++it) {
    const double xm = 0.5 * (a + b);
    const double tol1 = abs_tol + kSqrtEps * std::abs(x);
    const double tol2 = 2.0 * tol1;
    if (std::abs(x - xm) <= tol2 - 0.5 * (b - a)) {
      break;
    }
    bool golden = true;
    if (std::abs(e) > tol1) {
      const double r = (x - w) * (fx - fv);
      double q = (x - v) * (fx - fw);
      double p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0) {
        p = -p;
      } else {
        q = -q;
      }
      const double e_prev = e;
      e = d;
      if (std::abs(p) < std::abs(0.5 * q * e_prev) && p > q * (a - x) && p < q * (b - x)) {
        d = p / q;
        const double u = x + d;
        if (u - a < tol2 || b - u < tol2) {
          d = (xm >= x) ? tol1 : -tol1;
        }
        golden = false;
      }
    }
    if (golden) {
      e = (x >= xm) ? a - x : b - x;
      d = kCGold * e;
    }
    const double u = (std::abs(d) >= tol1) ? x + d : x + ((d >= 0.0) ? tol1 : -tol1);
    const double fu = fn(u);
    if (fu <= 0.0) {
      return u;
    }
    if (fu <= fx) {
      if (u >= x) {
        a = x;
      } else {
        b = x;
      }
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    } else {
      if (u < x) {
        a = u;
      } else {
        b = u;
      }
      if (fu <= fw || w == x) {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }
  return std::nan("");
}

struct RayleighDuct {
  double G, Tt0, Tt_exit, L, D, rough, fmult;
  const std::vector<double>* X;
  const std::string* model;
  // The march's last static T: the next stage's Newton start (a cost, not a
  // result: see static_T).
  mutable double T_hint = -1.0;
};

// Forward-mode discrete tangent of the fixed 32-step RK4 march (#1276).
// Single DualN<4> instantiation across unchoked and kExit paths (kThroat keeps
// its smooth G* re-solve):
//   seed 0: G,  seed 1: Pt0,  seed 2: Tt0,  seed 3: Tt_exit.
using Dual4 = solver::DualN<4>;

Dual4 lift_dual(double val, double slope, const Dual4& arg) {
  Dual4 r;
  r.v = val;
  for (int i = 0; i < 4; ++i) {
    r.d[i] = slope * arg.d[i];
  }
  return r;
}

// Fanno-regular Mach coordinate w(M^2) = 1/M^2 - 1 + ln(M^2) >= 0 (#1337).
// Near M = 1 (eps = 1 - M^2 -> 0+): w = eps^2/2 + 2 eps^3/3 + 3 eps^4/4 + ...
// (~ zeta / 2), regular at sonic; at low Mach w ~ 1/M^2, linear in x for Fanno
// flow across arbitrary f L / D (removing the Riccati blow-up that caused
// backward theta-RK4 to overshoot below zero on long laminar ducts).
double w_of_mach2(double M2) {
  const double eps = 1.0 - M2;
  if (!(eps > 0.0)) {
    return 0.0;
  }
  if (eps < 1e-4) {
    return eps * eps *
           (0.5 + eps * (2.0 / 3.0 + eps * (0.75 + eps * (0.8 + eps * (5.0 / 6.0)))));
  }
  return eps / M2 + std::log(M2);
}

// Solve d - ln(1 + d) = w for d = 1/M^2 - 1 > 0 (w > 0). Convex and monotone
// for d > 0: Newton converges to machine precision in 2-5 iterations across all
// positive doubles without cancellation as w -> 0+ or w -> infinity.
double d_of_w(double w) {
  double dv = 0.0;
  if (w < 1e-4) {
    const double s = std::sqrt(2.0 * w);
    dv = s * (1.0 + s * (1.0 / 3.0 + s / 36.0));
  } else {
    dv = w + std::log(w + 1.0 + std::log(w + 2.0));
  }
  for (int it = 0; it < 12; ++it) {
    const double fv =
        (dv < 1e-3)
            ? (dv * dv * (0.5 - dv * (1.0 / 3.0 - dv * (0.25 - dv * (0.2 - dv / 6.0)))) - w)
            : (dv - std::log1p(dv) - w);
    const double fp = dv / (1.0 + dv);
    const double step = fv / fp;
    dv -= step;
    if (std::abs(step) <= 1e-15 * dv) {
      break;
    }
  }
  return dv;
}

// M^2 from w; w <= 0 (a probe past sonic) is evaluated at M = 1, so the march
// -- and the choke test -- continue smoothly through it.
double mach2_of_w(double w) {
  if (!(w > 0.0)) {
    return 1.0;
  }
  return 1.0 / (1.0 + d_of_w(w));
}

// M^2(w) with exact derivative d(M^2)/dw = -1 / (d (1 + d)), where d = 1/M^2 - 1.
// For w <= 0 (including the sonic start w = Dual4::constant(0.0) of
// back_w_inlet_dual), returns Dual4::constant(1.0) without 0 * inf = NaN.
Dual4 mach2_of_w_dual(const Dual4& w) {
  if (!(w.v > 0.0)) {
    return Dual4::constant(1.0);
  }
  const double dv = d_of_w(w.v);
  const double inv_M2 = 1.0 + dv;
  return lift_dual(1.0 / inv_M2, -1.0 / (dv * inv_M2), w);
}

// zeta = (1 - M^2)^2 = (d / (1 + d))^2 for w > 0, extended C1 across w = 0 as
// 2 w for w <= 0 (matching zeta(0+) = 0 and dzeta/dw(0+) = 2).
double zeta_of_w(double w) {
  if (!(w > 0.0)) {
    return 2.0 * w;
  }
  const double dv = d_of_w(w);
  const double r = dv / (1.0 + dv);
  return r * r;
}

// Dimensionless friction length lam = f_ref L / D evaluated at exit stagnation
// temperature (ideal-gas viscosity is pressure-independent; Pt0 is passed for
// consistency with choked_flux_scale).
template <typename Scalar>
Scalar rayleigh_friction_length(const RayleighDuct& d, const Scalar& G_val, const Scalar& Tte_val,
                                double Pt0) {
  constexpr bool kIsDual = !std::is_same_v<Scalar, double>;
  if (*d.model == "fixed") {
    const double lam = std::max(d.fmult, 0.0) * d.L / d.D;
    if constexpr (kIsDual) {
      return Dual4::constant(lam);
    } else {
      return lam;
    }
  }
  if (!(d.fmult > 0.0)) {
    if constexpr (kIsDual) {
      return Dual4::constant(0.0);
    } else {
      return 0.0;
    }
  }
  const std::vector<double>& X = *d.X;
  if constexpr (kIsDual) {
    const TransportDerivatives tr = transport_and_dT(Tte_val.v, Pt0, X);
    const Dual4 mu = lift_dual(tr.mu, tr.dmu_dT, Tte_val);
    const Dual4 Re = G_val * d.D / mu;
    const auto [f0, df_dRe] = solver::friction_and_jacobian(*d.model, Re.v, d.rough / d.D).result;
    const Dual4 f_ref = lift_dual(f0 * d.fmult, df_dRe * d.fmult, Re);
    return f_ref * (d.L / d.D);
  } else {
    const double mu = viscosity(Tte_val, Pt0, X);
    const double f_ref = friction(G_val * d.D / mu, d.rough / d.D, *d.model, d.fmult);
    return f_ref * d.L / d.D;
  }
}

// Backward-march stretching exponent (sonic exit, w_exit = 0):
//   alpha = ln(1 + f_ref L / D).
template <typename Scalar>
Scalar rayleigh_grid_alpha_back(const RayleighDuct& d, const Scalar& G_val, const Scalar& Tte_val,
                                double Pt0) {
  constexpr bool kIsDual = !std::is_same_v<Scalar, double>;
  const Scalar lam = rayleigh_friction_length<Scalar>(d, G_val, Tte_val, Pt0);
  if constexpr (kIsDual) {
    return lift_dual(std::log1p(lam.v), 1.0 / (1.0 + lam.v), lam);
  } else {
    return std::log1p(lam);
  }
}

// Forward-march stretching exponent gated by sonic proximity:
//   lam_eff = lam / (1 + w_in / lam),  alpha = ln(1 + lam_eff).
// Matches alpha_back at w_in = 0, stays O(ln(lam)) near choke where w_in ~
// gamma * lam, and vanishes smoothly as lam^2 / w_in at low Mach.
template <typename Scalar>
Scalar rayleigh_grid_alpha_fwd(const RayleighDuct& d, const Scalar& G_val, const Scalar& Tte_val,
                               const Scalar& w_in, double Pt0) {
  constexpr bool kIsDual = !std::is_same_v<Scalar, double>;
  const Scalar lam = rayleigh_friction_length<Scalar>(d, G_val, Tte_val, Pt0);
  if constexpr (kIsDual) {
    if (!(lam.v > 0.0)) {
      return Dual4::constant(0.0);
    }
    const Dual4 lam_eff = lam / (1.0 + w_in / lam);
    return lift_dual(std::log1p(lam_eff.v), 1.0 / (1.0 + lam_eff.v), lam_eff);
  } else {
    if (!(lam > 0.0)) {
      return 0.0;
    }
    const double lam_eff = lam / (1.0 + w_in / lam);
    return std::log1p(lam_eff);
  }
}

template <typename Scalar>
Scalar rayleigh_node_x(int i, double L, const Scalar& alpha) {
  constexpr bool kIsDual = !std::is_same_v<Scalar, double>;
  const double t = 1.0 - static_cast<double>(i) / kRayleighSteps;
  const double s = 1.0 - t;
  const double tp = std::pow(t, kRayleighGrading);
  const double sp = std::pow(s, kRayleighGrading);
  if constexpr (kIsDual) {
    const double e = std::exp(-alpha.v * sp);
    const Dual4 frac = lift_dual(tp * e, -tp * sp * e, alpha);
    return L * (1.0 - frac);
  } else {
    return L * (1.0 - tp * std::exp(-alpha * sp));
  }
}

// Single-source dw/dx = -((1 - M^2) / M^4) d(M^2)/dx templated over Scalar
// (double for primal marches, Dual4 for exact forward-mode discrete tangents).
// Both (1 - M^2) and M^4 cancel identically, leaving
//   dw/dx = -(2 k / (a cv)) [f a / (2 D) + (C / M) dTt/dx],
// whose friction term is -(k / cv) (f / D) -- O(1) across all subsonic Mach.
template <typename Scalar>
Scalar rayleigh_dw_dx_impl(const Scalar& w, const Scalar& Tt, const Scalar& dTt_dx,
                           const Scalar& G_val, const RayleighDuct& d) {
  constexpr bool kIsDual = !std::is_same_v<Scalar, double>;
  const std::vector<double>& X = *d.X;
  auto val_of = [](const Scalar& x) -> double {
    if constexpr (kIsDual) {
      return x.v;
    } else {
      return x;
    }
  };
  auto sqrt_of = [](const Scalar& x) -> Scalar {
    if constexpr (kIsDual) {
      return solver::dsqrt(x);
    } else {
      return std::sqrt(x);
    }
  };

  Scalar M2, inv_M2;
  if constexpr (kIsDual) {
    if (w.v > 0.0) {
      const double dv = d_of_w(w.v);
      const double inv_v = 1.0 + dv;
      inv_M2 = lift_dual(inv_v, inv_v / dv, w);
      M2 = lift_dual(1.0 / inv_v, -1.0 / (dv * inv_v), w);
    } else {
      inv_M2 = Dual4::constant(1.0);
      M2 = Dual4::constant(1.0);
    }
  } else {
    if (w > 0.0) {
      inv_M2 = 1.0 + d_of_w(w);
      M2 = 1.0 / inv_M2;
    } else {
      inv_M2 = 1.0;
      M2 = 1.0;
    }
  }
  const Scalar M = sqrt_of(M2);
  const double T_v = static_T(val_of(Tt), val_of(M), X, d.T_hint);
  d.T_hint = T_v;

  const double R = specific_gas_constant(X);
  const double cpT_v = cp_mass(T_v, X);
  const double cpTt_v = cp_mass(val_of(Tt), X);
  const double a2_v = sound_speed_sq(T_v, X);
  const double da2_v = sound_speed_sq_dT(T_v, X);

  Scalar T, a2, da2, cpT, cpTt;
  if constexpr (kIsDual) {
    const double inv_mw = 1000.0 / mwmix(X);
    const double dcpT_v = dcp_dT(T_v, X) * inv_mw;
    const double dcpTt_v = dcp_dT(Tt.v, X) * inv_mw;
    const double d2a2_v = sound_speed_sq_d2T(T_v, X);
    // Implicit derivative of h(T) + 0.5 * M^2 * a^2(T) - h(Tt) = 0.
    const double k_v = cpT_v + 0.5 * M2.v * da2_v;
    T.v = T_v;
    for (int i = 0; i < 4; ++i) {
      T.d[i] = (cpTt_v * Tt.d[i] - 0.5 * a2_v * M2.d[i]) / k_v;
    }
    a2 = lift_dual(a2_v, da2_v, T);
    da2 = lift_dual(da2_v, d2a2_v, T);
    cpT = lift_dual(cpT_v, dcpT_v, T);
    cpTt = lift_dual(cpTt_v, dcpTt_v, Tt);
  } else {
    T = T_v;
    a2 = a2_v;
    da2 = da2_v;
    cpT = cpT_v;
    cpTt = cpTt_v;
  }

  const Scalar a = sqrt_of(a2);
  const Scalar a_T = 0.5 * da2 / a;
  const Scalar cv = cpT - R;
  const Scalar k = cpT + 0.5 * M2 * da2;
  const Scalar u = M * a;

  Scalar f;
  if constexpr (kIsDual) {
    f = Dual4::constant(d.fmult);
  } else {
    f = d.fmult;
  }
  if (*d.model != "fixed") {
    const double P_v = d.G * R * T_v / val_of(u);
    if constexpr (kIsDual) {
      const TransportDerivatives tr = transport_and_dT(T_v, P_v, X);
      const Dual4 mu = lift_dual(tr.mu, tr.dmu_dT, T);
      const Dual4 Re = G_val * d.D / mu;
      const auto [f0, df_dRe] = solver::friction_and_jacobian(*d.model, Re.v, d.rough / d.D).result;
      f = lift_dual(f0 * d.fmult, df_dRe * d.fmult, Re);
    } else {
      f = friction(G_val * d.D / viscosity(T_v, P_v, X), d.rough / d.D, *d.model, d.fmult);
    }
  }

  const Scalar gamma = cpT / cv;
  const Scalar fric = f * a / (2.0 * d.D);
  const Scalar heat = (cpTt / k) * (a_T + inv_M2 * (R / a - a_T / gamma)) * dTt_dx;
  return -(2.0 * k / (a * cv)) * (fric + heat);
}

// dw/dx at (w, Tt).
double rayleigh_dw_dx(double w, double Tt, double dTt_dx, const RayleighDuct& d) {
  return rayleigh_dw_dx_impl<double>(w, Tt, dTt_dx, d.G, d);
}

struct MarchEnd {
  double M_in = 0.0;
  double w_exit = 0.0;
  double zeta_min = 1.0;    // zeta_of_w(min(w)): <= 0 once the duct chokes
  bool max_inside = false;  // min(w) (max Mach) at a node before the exit
  // zeta at the C1 cubic Hermite minimum around an inside peak: the throat's
  // choke relation. The node minimum itself jumps as the peak moves between
  // nodes, which put 0.6-36% into its derivatives.
  double zeta_peak = 1.0;
};

MarchEnd march_zeta(const RayleighDuct& d, double Pt0) {
  MarchEnd out;
  out.M_in = fanno_inlet_mach(Pt0, d.Tt0, d.G, *d.X);
  double w = w_of_mach2(out.M_in * out.M_in);
  double wmin = w;
  int imin = 0;
  std::array<double, kRayleighSteps + 1> ws{};
  std::array<double, kRayleighSteps + 1> dws{};
  ws[0] = w;
  const double dTt_dx = (d.Tt_exit - d.Tt0) / d.L;
  const double alpha = rayleigh_grid_alpha_fwd<double>(d, d.G, d.Tt_exit, w, Pt0);
  auto x_at = [&](int i) { return rayleigh_node_x<double>(i, d.L, alpha); };
  auto Tt_at = [&](double x) { return d.Tt0 + dTt_dx * x; };
  for (int i = 0; i < kRayleighSteps; ++i) {
    const double x0 = x_at(i);
    const double x1 = x_at(i + 1);
    const double h = x1 - x0;
    const double k1 = rayleigh_dw_dx(w, Tt_at(x0), dTt_dx, d);
    dws[i] = k1;
    const double k2 = rayleigh_dw_dx(w + 0.5 * h * k1, Tt_at(x0 + 0.5 * h), dTt_dx, d);
    const double k3 = rayleigh_dw_dx(w + 0.5 * h * k2, Tt_at(x0 + 0.5 * h), dTt_dx, d);
    const double k4 = rayleigh_dw_dx(w + h * k3, Tt_at(x1), dTt_dx, d);
    w += h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    ws[i + 1] = w;
    if (w < wmin) {
      wmin = w;
      imin = i + 1;
    }
  }
  out.w_exit = w;
  out.zeta_min = zeta_of_w(wmin);
  out.max_inside = imin < kRayleighSteps;
  out.zeta_peak = out.zeta_min;
  if (out.max_inside) {
    // Cubic Hermite minimum on the interval bracketing dw/dx = 0 around imin:
    // fourth-order accurate and C1 across node crossings (by the envelope
    // theorem, matching value and first derivatives when dw/dx = 0 at a node),
    // avoiding 3-node Lagrange stencil jumps when imin shifts.
    dws[kRayleighSteps] = rayleigh_dw_dx(w, d.Tt_exit, dTt_dx, d);
    int k = -1;
    if (imin > 0 && dws[imin - 1] <= 0.0 && dws[imin] >= 0.0) {
      k = imin - 1;
    } else if (dws[imin] <= 0.0 && dws[imin + 1] >= 0.0) {
      k = imin;
    }
    if (k >= 0) {
      const double h = x_at(k + 1) - x_at(k);
      const double y0 = ws[k], y1 = ws[k + 1];
      const double s0 = dws[k], s1 = dws[k + 1];
      const double a = 6.0 * (y0 - y1) + 3.0 * h * (s0 + s1);
      const double b = -6.0 * (y0 - y1) - 4.0 * h * s0 - 2.0 * h * s1;
      const double c = h * s0;
      const double disc = std::max(0.0, b * b - 4.0 * a * c);
      const double r = std::sqrt(disc);
      const double denom = b + r;
      const double t =
          (b >= 0.0) ? ((denom > 0.0) ? (-2.0 * c / denom) : 0.0) : ((-b + r) / (2.0 * a));
      const double tc = std::clamp(t, 0.0, 1.0);
      const double yv = y0 + tc * (h * s0 + tc * ((3.0 * (y1 - y0) - h * (2.0 * s0 + s1)) +
                                                  tc * (2.0 * (y0 - y1) + h * (s0 + s1))));
      out.zeta_peak = zeta_of_w(yv);
    }
  }
  return out;
}

double exit_mach(const MarchEnd& m) { return std::sqrt(mach2_of_w(m.w_exit)); }

// Exit pressure (static, or stagnation if exit_total) of the march at G.
double rayleigh_exit_pressure(const RayleighDuct& d, const MarchEnd& m, bool exit_total) {
  return d.G * exit_per_flux(exit_mach(m), d.Tt_exit, *d.X, exit_total);
}

struct RayleighProblem {
  double Pt0, Tt0, Tt_exit, P, L, D, rough, fmult;
  const std::vector<double>* X;
  const std::string* model;
  bool exit_total;
  double G_guess;
  const std::vector<std::vector<double>>* dirs;
  bool with_derivatives;

  RayleighDuct duct(double G, double Tt0x, double Tt_ex, const std::vector<double>* Xx) const {
    return RayleighDuct{G, Tt0x, Tt_ex, L, D, rough, fmult, Xx, model};
  }
};

// w at the inlet of the trajectory that is sonic exactly at the exit, marched
// BACKWARD from w = 0 (M = 1) at x = L. Going upstream w only rises (dw/dx < 0
// at a reachable sonic exit), smoothly and linearly in x for Fanno friction
// across arbitrary f L / D (#1337), so no stage ever leaves the subsonic
// branch or overshoots M = 0.
double back_w_inlet(const RayleighDuct& d, double Pt0) {
  double w = 0.0;
  const double dTt_dx = (d.Tt_exit - d.Tt0) / d.L;
  const double alpha = rayleigh_grid_alpha_back<double>(d, d.G, d.Tt_exit, Pt0);
  auto x_at = [&](int i) { return rayleigh_node_x<double>(i, d.L, alpha); };
  auto Tt_at = [&](double x) { return d.Tt0 + dTt_dx * x; };
  for (int i = kRayleighSteps; i > 0; --i) {
    const double x0 = x_at(i);
    const double x1 = x_at(i - 1);
    const double h = x1 - x0;
    const double k1 = rayleigh_dw_dx(w, Tt_at(x0), dTt_dx, d);
    const double k2 = rayleigh_dw_dx(w + 0.5 * h * k1, Tt_at(x0 + 0.5 * h), dTt_dx, d);
    const double k3 = rayleigh_dw_dx(w + 0.5 * h * k2, Tt_at(x0 + 0.5 * h), dTt_dx, d);
    const double k4 = rayleigh_dw_dx(w + h * k3, Tt_at(x1), dTt_dx, d);
    w += h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
  }
  return w;
}

double inlet_w(double Pt0, double Tt0, double G, const std::vector<double>& X) {
  const double M = fanno_inlet_mach(Pt0, Tt0, G, X);
  return w_of_mach2(M * M);
}

// w_inlet(G, Pt0, Tt0) with exact implicit derivatives from
// isentropic_flux(Pt0, Tt0, M_in, X) = G.
Dual4 inlet_w_dual(double Pt0, double Tt0, double G, const std::vector<double>& X) {
  const double M = fanno_inlet_mach(Pt0, Tt0, G, X);
  const double M2 = M * M;
  Dual4 out;
  out.v = w_of_mach2(M2);
  if (M <= 0.0 || M >= 1.0) {
    return out;
  }
  const FluxSlopes fs = isentropic_flux_slopes(Pt0, Tt0, M, X);
  // dw/dM = dw/d(M^2) * 2M = -2 (1 - M^2) / M^3
  const double dw_dM = -2.0 * (1.0 - M2) / (M2 * M);
  const double scale = dw_dM / fs.dM;
  out.d[0] = scale;
  out.d[1] = -scale * (G / Pt0);
  out.d[2] = -scale * fs.dTt;
  out.d[3] = 0.0;
  return out;
}

// Exact tangent of rayleigh_dw_dx with respect to (G, Pt0, Tt0, Tt_exit).
Dual4 rayleigh_dw_dx_dual(const Dual4& w, const Dual4& Tt, const Dual4& dTt_dx,
                          const Dual4& G_dual, const RayleighDuct& d) {
  return rayleigh_dw_dx_impl<Dual4>(w, Tt, dTt_dx, G_dual, d);
}

// Forward RK4 march in Dual4 returning exit w with exact partials in
// (G, Pt0, Tt0, Tt_exit) for the unchoked discrete tangent (#1276).
Dual4 march_w_exit_dual(const RayleighDuct& d, double Pt0) {
  const Dual4 G_dual = Dual4::seed(d.G, 0);
  const Dual4 Tt0_dual = Dual4::seed(d.Tt0, 2);
  const Dual4 Tte_dual = Dual4::seed(d.Tt_exit, 3);
  const Dual4 dTt_dx = (Tte_dual - Tt0_dual) / d.L;

  Dual4 w = inlet_w_dual(Pt0, d.Tt0, d.G, *d.X);
  const Dual4 alpha = rayleigh_grid_alpha_fwd<Dual4>(d, G_dual, Tte_dual, w, Pt0);
  auto x_at = [&](int i) { return rayleigh_node_x<Dual4>(i, d.L, alpha); };
  auto Tt_at = [&](const Dual4& x) { return Tt0_dual + dTt_dx * x; };

  for (int i = 0; i < kRayleighSteps; ++i) {
    const Dual4 x0 = x_at(i);
    const Dual4 x1 = x_at(i + 1);
    const Dual4 h = x1 - x0;
    const Dual4 k1 = rayleigh_dw_dx_dual(w, Tt_at(x0), dTt_dx, G_dual, d);
    const Dual4 k2 =
        rayleigh_dw_dx_dual(w + 0.5 * h * k1, Tt_at(x0 + 0.5 * h), dTt_dx, G_dual, d);
    const Dual4 k3 =
        rayleigh_dw_dx_dual(w + 0.5 * h * k2, Tt_at(x0 + 0.5 * h), dTt_dx, G_dual, d);
    const Dual4 k4 = rayleigh_dw_dx_dual(w + h * k3, Tt_at(x1), dTt_dx, G_dual, d);
    w = w + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
  }
  return w;
}

Dual4 back_w_inlet_dual(const RayleighDuct& d, double Pt0) {
  const Dual4 G_dual = Dual4::seed(d.G, 0);
  const Dual4 Tt0_dual = Dual4::seed(d.Tt0, 2);
  const Dual4 Tte_dual = Dual4::seed(d.Tt_exit, 3);
  const Dual4 dTt_dx = (Tte_dual - Tt0_dual) / d.L;
  const Dual4 alpha = rayleigh_grid_alpha_back<Dual4>(d, G_dual, Tte_dual, Pt0);

  Dual4 w = Dual4::constant(0.0);
  auto x_at = [&](int i) { return rayleigh_node_x<Dual4>(i, d.L, alpha); };
  auto Tt_at = [&](const Dual4& x) { return Tt0_dual + dTt_dx * x; };
  for (int i = kRayleighSteps; i > 0; --i) {
    const Dual4 x0 = x_at(i);
    const Dual4 x1 = x_at(i - 1);
    const Dual4 h = x1 - x0;
    const Dual4 k1 = rayleigh_dw_dx_dual(w, Tt_at(x0), dTt_dx, G_dual, d);
    const Dual4 k2 =
        rayleigh_dw_dx_dual(w + 0.5 * h * k1, Tt_at(x0 + 0.5 * h), dTt_dx, G_dual, d);
    const Dual4 k3 =
        rayleigh_dw_dx_dual(w + 0.5 * h * k2, Tt_at(x0 + 0.5 * h), dTt_dx, G_dual, d);
    const Dual4 k4 = rayleigh_dw_dx_dual(w + h * k3, Tt_at(x1), dTt_dx, G_dual, d);
    w = w + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
  }
  return w;
}

Dual4 rayleigh_exit_pressure_dual(const RayleighDuct& d, const Dual4& w_exit, bool exit_total) {
  const Dual4 G_dual = Dual4::seed(d.G, 0);
  const Dual4 Me2 = mach2_of_w_dual(w_exit);
  const Dual4 Me = solver::dsqrt(Me2);
  const double phi_v = exit_per_flux(Me.v, d.Tt_exit, *d.X, exit_total);
  const PhiSlopes ps = exit_per_flux_slopes(Me.v, d.Tt_exit, *d.X, exit_total);
  Dual4 phi;
  phi.v = phi_v;
  for (int i = 0; i < 4; ++i) {
    phi.d[i] = ps.dM * Me.d[i];
  }
  phi.d[3] += ps.dTt;
  return G_dual * phi;
}

// Composition tangents of the march (#1310), on the mixture primitives of
// #1311: the same fixed RK4 grid in w (#1337), one pass per direction v at
// fixed (G, Pt0, Tt0, Tt_exit), so the tangent of the discrete march is
// exact. The grid's stretching moves with composition too: through the
// friction length (mu at Tt_exit) and, forward, the inlet w. They replace
// d_ddir's step-cut differences of the relation, whose cancellation in
// E(G; X) left up to 12 % in dG/dY near the 1 Pa floor.

// M^2 and 1/M^2 from w with their tangents: d(M^2)/dw = -1 / (d (1 + d)) and
// d(1/M^2)/dw = (1 + d) / d, d = 1/M^2 - 1; w <= 0 (a probe past sonic) is
// sonic, constant, as in rayleigh_dw_dx_impl.
template <int N>
std::pair<DualT<N>, DualT<N>> mach2_of_w_t(const DualT<N>& w) {
  if (!(w.v > 0.0)) {
    return {DualT<N>::constant(1.0), DualT<N>::constant(1.0)};
  }
  const double dv = d_of_w(w.v);
  const double inv_v = 1.0 + dv;
  DualT<N> M2 = DualT<N>::constant(1.0 / inv_v);
  DualT<N> inv = DualT<N>::constant(inv_v);
  for (int i = 0; i < N; ++i) {
    M2.d[i] = -w.d[i] / (dv * inv_v);
    inv.d[i] = w.d[i] * inv_v / dv;
  }
  return {M2, inv};
}

template <int N>
DualT<N> rayleigh_dw_dx_t(const DualT<N>& w, const DualT<N>& Tt, double dTt_dx, double G,
                          const RayleighDuct& d, const MixTangent<N>& m) {
  const auto [M2, inv_M2] = mach2_of_w_t(w);
  const DualT<N> M = solver::dsqrt(M2);
  const DualT<N> T = static_T_t(Tt, M, m);
  const DualT<N> R = gas_constant_t(m);
  const DualT<N> cpT = cp_mass_t(T, m);
  const DualT<N> cpTt = cp_mass_t(Tt, m);
  const DualT<N> a2 = sound_speed_sq_t(T, m);
  const DualT<N> da2 = sound_speed_sq_dT_t(T, m);
  const DualT<N> a = solver::dsqrt(a2);
  const DualT<N> a_T = 0.5 * da2 / a;
  const DualT<N> cv = cpT - R;
  const DualT<N> k = cpT + 0.5 * M2 * da2;
  const DualT<N> u = M * a;
  DualT<N> f = DualT<N>::constant(d.fmult);
  if (*d.model != "fixed") {
    const double P_v = G * R.v * T.v / u.v;
    const TransportDerivatives tr = transport_and_dT(T.v, P_v, *m.X);
    double dmu_dv = 0.0;
    if (m.has_comp()) {
      const std::vector<double> g = viscosity_dX(T.v, P_v, *m.X);
      for (std::size_t j = 0; j < g.size(); ++j) {
        dmu_dv += g[j] * m.dX[j];
      }
    }
    const DualT<N> mu = lift_tx(tr.mu, tr.dmu_dT, T, dmu_dv, m);
    const DualT<N> Re = G * d.D / mu;
    const auto [f0, df_dRe] = solver::friction_and_jacobian(*d.model, Re.v, d.rough / d.D).result;
    f.v = f0 * d.fmult;
    for (int i = 0; i < N; ++i) {
      f.d[i] = df_dRe * d.fmult * Re.d[i];
    }
  }
  const DualT<N> gamma = cpT / cv;
  const DualT<N> fric = f * a / (2.0 * d.D);
  const DualT<N> heat = (cpTt / k) * (a_T + inv_M2 * (R / a - a_T / gamma)) * dTt_dx;
  return -(2.0 * k / (a * cv)) * (fric + heat);
}

// rayleigh_friction_length along v: f_ref(Re at Tt_exit) L / D, with mu's
// composition partial; the same Pt0 and Tt_exit as the primal.
DualT<1> rayleigh_friction_length_comp(const RayleighDuct& d, double Pt0, const MixTangent<1>& m) {
  if (*d.model == "fixed") {
    return DualT<1>::constant(std::max(d.fmult, 0.0) * d.L / d.D);
  }
  if (!(d.fmult > 0.0)) {
    return DualT<1>::constant(0.0);
  }
  const std::vector<double>& X = *m.X;
  const double mu_v = viscosity(d.Tt_exit, Pt0, X);
  double dmu_dv = 0.0;
  const std::vector<double> g = viscosity_dX(d.Tt_exit, Pt0, X);
  for (std::size_t j = 0; j < g.size(); ++j) {
    dmu_dv += g[j] * m.dX[j];
  }
  const double Re = d.G * d.D / mu_v;
  const auto [f0, df_dRe] = solver::friction_and_jacobian(*d.model, Re, d.rough / d.D).result;
  DualT<1> lam = DualT<1>::constant(f0 * d.fmult * d.L / d.D);
  lam.d[0] = df_dRe * d.fmult * (-Re * dmu_dv / mu_v) * d.L / d.D;
  return lam;
}

// alpha = ln(1 + lam_eff): lam_eff = lam (backward) or lam / (1 + w_in / lam)
// (forward), as rayleigh_grid_alpha_back / _fwd.
DualT<1> rayleigh_grid_alpha_comp(const DualT<1>& lam, const DualT<1>* w_in) {
  if (w_in != nullptr && !(lam.v > 0.0)) {
    return DualT<1>::constant(0.0);
  }
  const DualT<1> lam_eff = (w_in == nullptr) ? lam : lam / (1.0 + *w_in / lam);
  DualT<1> alpha = DualT<1>::constant(std::log1p(lam_eff.v));
  alpha.d[0] = lam_eff.d[0] / (1.0 + lam_eff.v);
  return alpha;
}

// rayleigh_node_x with alpha's tangent: x = L (1 - t^g exp(-alpha s^g)).
DualT<1> rayleigh_node_x_comp(int i, double L, const DualT<1>& alpha) {
  const double t = 1.0 - static_cast<double>(i) / kRayleighSteps;
  const double s = 1.0 - t;
  const double tp = std::pow(t, kRayleighGrading);
  const double sp = std::pow(s, kRayleighGrading);
  const double e = std::exp(-alpha.v * sp);
  DualT<1> x = DualT<1>::constant(L * (1.0 - tp * e));
  x.d[0] = L * tp * sp * e * alpha.d[0];
  return x;
}

// w_in along v: isentropic_flux(Pt0, Tt0, M_in, X) = G at fixed G, so
// dM_in/dv = -flux_v / flux_M, and dw/dM = -2 (1 - M^2) / M^3.
DualT<1> inlet_w_comp(double Pt0, double Tt0, double G, const MixTangent<1>& m) {
  const std::vector<double>& X = *m.X;
  const double M = fanno_inlet_mach(Pt0, Tt0, G, X);
  const double M2 = M * M;
  DualT<1> out = DualT<1>::constant(w_of_mach2(M2));
  if (M <= 0.0 || M >= 1.0) {
    return out;
  }
  const double flux_v =
      isentropic_flux_t(Pt0, DualT<1>::constant(Tt0), DualT<1>::constant(M), m).d[0];
  const double dw_dM = -2.0 * (1.0 - M2) / (M2 * M);
  out.d[0] = -dw_dM * flux_v / isentropic_flux_slopes(Pt0, Tt0, M, X).dM;
  return out;
}

// The RK4 march along v, forward from the inlet (march_zeta's grid) or
// backward from a sonic exit (back_w_inlet's grid).
DualT<1> march_w_comp(const RayleighDuct& d, double Pt0, const MixTangent<1>& m, bool backward) {
  const double dTt_dx = (d.Tt_exit - d.Tt0) / d.L;
  DualT<1> w = backward ? DualT<1>::constant(0.0) : inlet_w_comp(Pt0, d.Tt0, d.G, m);
  const DualT<1> lam = rayleigh_friction_length_comp(d, Pt0, m);
  const DualT<1> alpha = rayleigh_grid_alpha_comp(lam, backward ? nullptr : &w);
  auto x_at = [&](int i) { return rayleigh_node_x_comp(i, d.L, alpha); };
  auto Tt_at = [&](const DualT<1>& x) { return d.Tt0 + dTt_dx * x; };
  auto step = [&](const DualT<1>& wv, const DualT<1>& x0, const DualT<1>& x1) {
    const DualT<1> h = x1 - x0;
    const DualT<1> k1 = rayleigh_dw_dx_t(wv, Tt_at(x0), dTt_dx, d.G, d, m);
    const DualT<1> k2 =
        rayleigh_dw_dx_t(wv + 0.5 * h * k1, Tt_at(x0 + 0.5 * h), dTt_dx, d.G, d, m);
    const DualT<1> k3 =
        rayleigh_dw_dx_t(wv + 0.5 * h * k2, Tt_at(x0 + 0.5 * h), dTt_dx, d.G, d, m);
    const DualT<1> k4 = rayleigh_dw_dx_t(wv + h * k3, Tt_at(x1), dTt_dx, d.G, d, m);
    return wv + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
  };
  if (backward) {
    for (int i = kRayleighSteps; i > 0; --i) {
      w = step(w, x_at(i), x_at(i - 1));
    }
  } else {
    for (int i = 0; i < kRayleighSteps; ++i) {
      w = step(w, x_at(i), x_at(i + 1));
    }
  }
  return w;
}

// d(relation)/dv at fixed G: the exit pressure G phi(M_exit, Tt_exit)
// unchoked; for an exit choke Psi = w_inlet - w_back(0), as rel_dual.
double rayleigh_relation_ddir(const RayleighDuct& d, double Pt0, bool exit_choked,
                              bool exit_total, const std::vector<double>& v) {
  const MixTangent<1> m = mix_tangent<1>(*d.X, &v, 0);
  if (exit_choked) {
    return inlet_w_comp(Pt0, d.Tt0, d.G, m).d[0] - march_w_comp(d, Pt0, m, true).d[0];
  }
  const DualT<1> w_exit = march_w_comp(d, Pt0, m, false);
  const DualT<1> Me = solver::dsqrt(mach2_of_w_t(w_exit).first);
  return d.G * exit_per_flux_t(Me, DualT<1>::constant(d.Tt_exit), m, exit_total).d[0];
}

// Where the duct chokes, which sets the relation that defines the choked flux:
//   kInlet   G is the inlet's own sonic flux (cooling keeps a sonic inlet
//            subsonic downstream)
//   kExit    Psi = w_inlet(G) - w_back(0) = 0: sonic exactly at L
//   kThroat  zeta_min = 0 inside the duct (cooling after an interior sonic
//            point), from the forward march
enum class ChokeMode : std::uint8_t { kInlet, kExit, kThroat };

struct Choke {
  double G = 0.0;
  ChokeMode mode = ChokeMode::kExit;
};

double choke_relation(ChokeMode mode, double G, double Pt0, const RayleighDuct& d) {
  const std::vector<double>& X = *d.X;
  switch (mode) {
    case ChokeMode::kInlet:
      return G - fanno_sonic_mass_flux(Pt0, d.Tt0, X);
    case ChokeMode::kExit:
      return inlet_w(Pt0, d.Tt0, G, X) - back_w_inlet(d, Pt0);
    case ChokeMode::kThroat:
      return march_zeta(d, Pt0).zeta_peak;
  }
  return 0.0;
}

Choke choked_flux(const RayleighProblem& pb, double Pt0, double Tt0, double Tt_ex,
                  const std::vector<double>& X) {
  const double G_max = fanno_sonic_mass_flux(Pt0, Tt0, X);
  const double hi0 = std::log(G_max * (1.0 - 1e-12));
  auto duct_at = [&](double logG) { return pb.duct(std::exp(logG), Tt0, Tt_ex, &X); };
  if (march_zeta(duct_at(hi0), Pt0).zeta_min >= 0.0) {
    return Choke{G_max, ChokeMode::kInlet};
  }
  // Lowest log(G) at which 1 / M_in^2 ~ (G_max / G)^2 stays finite in double.
  const double log_G_min = std::log(G_max) + 0.5 * std::log(std::numeric_limits<double>::min());
  const double G_est =
      choked_flux_scale(Pt0, Tt0, Tt_ex, X, pb.L, pb.D, pb.rough, *pb.model, pb.fmult, G_max);
  auto solve = [&](ChokeMode mode) {
    auto fn = [&](double logG) { return choke_relation(mode, std::exp(logG), Pt0, duct_at(logG)); };
    double hi = hi0;
    double f_hi = fn(hi);
    double lo = std::clamp(std::log(G_est) - 1.0, log_G_min, hi - 1.0);
    double f_lo = fn(lo);
    double step = 2.0;
    while (f_lo <= 0.0 && lo > log_G_min) {
      hi = lo;
      f_hi = f_lo;
      lo = std::max(lo - step, log_G_min);
      f_lo = fn(lo);
      step *= 2.0;
    }
    return std::exp(brent(fn, lo, hi, f_lo, f_hi, 1e-14));
  };
  // A sonic exit is reachable only where friction and heat still accelerate
  // the flow at M = 1 (dw/dx < 0 there).
  const RayleighDuct d_hi = duct_at(hi0);
  if (rayleigh_dw_dx(0.0, Tt_ex, (Tt_ex - Tt0) / pb.L, d_hi) < 0.0) {
    const double Ge = solve(ChokeMode::kExit);
    // At that flux the forward march must peak AT the exit; a peak inside
    // that already reaches sonic is a throat. (Its value alone cannot tell:
    // near sonic the forward march's own error is ~1e-5 in zeta.)
    const MarchEnd me = march_zeta(pb.duct(Ge, Tt0, Tt_ex, &X), Pt0);
    if (!(me.max_inside && me.zeta_min < 0.0)) {
      return Choke{Ge, ChokeMode::kExit};
    }
  }
  return Choke{solve(ChokeMode::kThroat), ChokeMode::kThroat};
}

FannoChannelFlow rayleigh_raw_flow(const RayleighProblem& pb) {
  const std::vector<double>& X = *pb.X;
  FannoChannelFlow out;
  // H(log G): the exit-pressure excess while the duct carries G subsonically,
  // and -Pt0 (1e-12 - zeta_min) once it chokes, so a sign change brackets
  // the unchoked root or the choke itself.
  // The last few marches, by flux: the subsonic checks and the final state
  // re-read fluxes F has just marched (a pure cache: same G, same march).
  std::array<std::pair<double, MarchEnd>, 4> memo{};
  std::size_t memo_n = 0;
  auto march_at = [&](double G) {
    for (std::size_t i = 0; i < std::min(memo_n, memo.size()); ++i) {
      if (memo[i].first == G) {
        return memo[i].second;
      }
    }
    const MarchEnd m = march_zeta(pb.duct(G, pb.Tt0, pb.Tt_exit, &X), pb.Pt0);
    memo[memo_n % memo.size()] = {G, m};
    ++memo_n;
    return m;
  };
  auto F = [&](double logG) {
    const double G = std::exp(logG);
    const MarchEnd m = march_at(G);
    if (m.zeta_min <= 0.0) {
      return pb.Pt0 * (m.zeta_min - 1e-12);
    }
    return rayleigh_exit_pressure(pb.duct(G, pb.Tt0, pb.Tt_exit, &X), m, pb.exit_total) - pb.P;
  };
  // The largest log G in (lo, hi) the march carries subsonically (zeta_min
  // above kRayleighZetaFloor), given it does at lo and not at hi: Illinois
  // regula falsi on zeta_min, to kRayleighSubsonicEndTol in log G. zeta_min
  // is close to linear in G there (slope ~ -12 to -35 per unit relative
  // flux), so it converges in a few marches; Illinois keeps the far end from
  // pinning a plain regula falsi, and a step outside the bracket bisects.
  // Returns the subsonic end, so the caller's F is the exit-pressure excess.
  auto subsonic_end = [&](double a, double b) {
    double za = march_at(std::exp(a)).zeta_min - kRayleighZetaFloor;
    double zb = march_at(std::exp(b)).zeta_min - kRayleighZetaFloor;
    int side = 0;
    for (int it = 0; it < 100 && b - a > kRayleighSubsonicEndTol; ++it) {
      double c = b - zb * (b - a) / (zb - za);
      if (!(c > a && c < b)) {
        c = 0.5 * (a + b);
      }
      const double zc = march_at(std::exp(c)).zeta_min - kRayleighZetaFloor;
      if (zc > 0.0) {
        a = c;
        za = zc;
        if (side == -1) {
          zb *= 0.5;
        }
        side = -1;
      } else {
        b = c;
        zb = zc;
        if (side == 1) {
          za *= 0.5;
        }
        side = 1;
      }
    }
    return a;
  };
  const double G_max = fanno_sonic_mass_flux(pb.Pt0, pb.Tt0, X);
  bool choked = false;
  ChokeMode mode = ChokeMode::kExit;
  double G = 0.0;
  bool done = false;
  // The fold search of a cooled duct (#1167, #1190, #1191).
  const bool may_turn = pb.Tt_exit < pb.Tt0;
  // The root on the branch continuous from zero flow: the FIRST fall of
  // the excess through zero, scanned up from low flux. The excess need not
  // be monotone: near choke a cooled duct's Rayleigh gain can outgrow its
  // friction (a short duct, entrance-enhanced heat transfer and fully
  // developed f: 0.76 of it on a 1.27-diameter impingement crossflow
  // segment), so it rises again toward the choke. Taking the choke whenever
  // the excess was still positive there, or the root nearest the choke,
  // jumped that duct to its sonic flux at a few Pa of drive. Only with no
  // fall below the choke does the duct choke.
  // A fall and rise of the excess inside one ln-2 bracket (a cooled duct
  // whose needed drive peaks and turns back) leaves both of its ends
  // positive: testing the ends alone read such a duct choked while a warm
  // start found its first crossing (#1167). Only cooling can turn the
  // excess back up (a Rayleigh gain), so on a cooled duct a bracket with
  // both ends positive is also searched: an interior one where the cubic
  // Hermite interpolant of its ends' values and slopes (kRayleighSlopeStep)
  // has an interior minimum, the last one (ending at the march's sonic
  // flux) at fluxes 4^-j of that end below it. Brent's minimiser looks for
  // a point at or below zero near each minimum, and the first fall is then
  // bracketed below that point. A cooled bracket whose ends straddle zero
  // is walked in kRayleighStraddleSplit parts with the same test (#1190),
  // and solved in the first part that holds a fall. An interior cooled
  // bracket with both ends positive and no minimum is walked the same way
  // where the needed drive at an end reaches kRayleighBothGate of the drive
  // (#1206). Heated and adiabatic ducts, and a cooled bracket with no
  // minimum below that gate, scan exactly as before. Not seen: two turns of
  // the excess between the ends of a part whose values and slopes show none.
  // first_fall(top) scans log G up to top, in ln-2 brackets from top - 7, and
  // returns the bracket of the first fall below top (or the root, where a
  // straddle split already solved it). The cold call scans up to the choke;
  // a cooled warm call scans below its bracket (#1191).
  struct FirstFall {
    bool found = false;
    bool past_subsonic = false;  // no fall up to the march's sonic flux
    double lo = 0.0, f_lo = 0.0, hi = 0.0, f_hi = 0.0;
    double root = std::nan("");  // log G of the root, where already solved
  };
  const double log_G_min = std::log(G_max) + 0.5 * std::log(std::numeric_limits<double>::min());
  auto first_fall = [&](double top) {
    double lo = std::max(top - 7.0, log_G_min);
    double f_lo = F(lo);
    double step = 7.0;
    while (f_lo <= 0.0 && lo > log_G_min) {
      lo = std::max(lo - step, log_G_min);
      f_lo = F(lo);
      step *= 2.0;
    }
    if (f_lo <= 0.0) {
      throw std::runtime_error("fanno-rayleigh (first_fall): no sign change");
    }
    // A fall counts only at a flux the march carries subsonically: at the
    // choke itself the march's own error can read as choked, and F encodes
    // that as negative too (it put Brent onto the choke, reported unchoked).
    // A bracket whose upper end the march reads as choked is cut back to the
    // largest flux it still carries subsonically, and the fall is tested
    // there: the march sonics ~3e-7 below the choke relation's flux, so the
    // last bracket's end at the choke always reads choked, and breaking there
    // skipped every root above its lower end, ~0.93 Gc (#1161).
    // With the excess still positive there, no subsonic flux matches the
    // drive, and the duct is choked even where the excess at Gc itself reads
    // negative: a static exit's E(G) is steep enough that the gap between
    // the two (with kRayleighZetaFloor) leaves a band of 0.1-8 Pa of drive
    // below the choke.
    bool found = false;
    double root = std::nan("");  // log G of the root, where already solved
    bool past_subsonic = false;
    // Slopes are dF/dlog G.
    auto slope = [&](double x, double fx) {
      return (fx - F(x - kRayleighSlopeStep)) / kRayleighSlopeStep;
    };
    // A point at or below zero, by minimising F over [a, b] of log G from
    // x0 inside it, or NaN: the first fall then lies between a and it.
    auto dip_in = [&](double a, double b, double x0) {
      return brent_min_to_zero(F, a, b, x0, kRayleighDipTol);
    };
    // Whether a both-positive bracket is walked (kRayleighBothGate, #1206).
    const double drive = pb.Pt0 - pb.P;
    auto walk_both = [&](double fa, double fb) {
      return std::max(std::abs(drive - fa), std::abs(drive - fb)) >= kRayleighBothGate * drive;
    };
    double s_lo = std::nan("");  // slope at lo, NaN not yet read
    double x_prev = std::nan(""), f_prev = 0.0;  // the scan point before lo
    double hi = lo, f_hi = f_lo;
    // The Hermite test of [lo, hi] (#1167): the cubic interpolant of its
    // ends' values and slopes (s_lo read here if not yet, then s_hi, which
    // it returns through its argument). Where the cubic has an interior
    // minimum, F is minimised over [lo, the cubic's next maximum, or hi]: a
    // point at or below zero there, or NaN.
    auto hermite_dip = [&](double& s_hi, bool s_hi_known = false) {
      if (std::isnan(s_lo)) {
        s_lo = slope(lo, f_lo);
      }
      if (!s_hi_known) {
        s_hi = slope(hi, f_hi);
      }
      const double H = hi - lo;
      const double qa = 6.0 * (f_lo - f_hi) + 3.0 * H * (s_lo + s_hi);
      const double qb = -6.0 * (f_lo - f_hi) - 4.0 * H * s_lo - 2.0 * H * s_hi;
      const double qc = H * s_lo;
      double t_min = std::nan(""), t_max = 1.0;
      if (std::abs(qa) > kRayleighHermiteDegenerate * (std::abs(qb) + std::abs(qc))) {
        const double disc = qb * qb - 4.0 * qa * qc;
        if (disc > 0.0) {
          const double r = std::sqrt(disc);
          const double t1 = (-qb - r) / (2.0 * qa), t2 = (-qb + r) / (2.0 * qa);
          for (const double t : {t1, t2}) {
            if (t > 0.0 && t < 1.0) {
              if (2.0 * qa * t + qb > 0.0) {
                t_min = t;
              }
            }
          }
          for (const double t : {t1, t2}) {
            if (t > 0.0 && t < 1.0 && 2.0 * qa * t + qb < 0.0 && t > t_min) {
              t_max = t;
            }
          }
        }
      } else if (qb > 0.0 && -qc / qb > 0.0 && -qc / qb < 1.0) {
        t_min = -qc / qb;
      }
      if (!std::isfinite(t_min)) {
        return std::nan("");
      }
      const double b = lo + t_max * H;
      const double x0 = lo + std::clamp(t_min, kRayleighDipStartMargin * t_max,
                                        (1.0 - kRayleighDipStartMargin) * t_max) *
                                 H;
      return dip_in(lo, b, x0);
    };
    while (hi < top) {
      x_prev = (hi > lo) ? lo : std::nan("");
      f_prev = f_lo;
      lo = hi;
      f_lo = f_hi;
      hi = std::min(top, lo + std::log(2.0));
      const bool cut = !(march_at(std::exp(hi)).zeta_min > kRayleighZetaFloor);
      if (cut) {
        hi = subsonic_end(lo, hi);
      }
      f_hi = F(hi);
      if (f_lo > 0.0 && f_hi <= 0.0) {
        found = true;
        if (may_turn) {
          // A cooled bracket whose ends straddle zero can still hold a
          // fall, a rise and a later fall, and Brent may take the later
          // one (#1190). Its ends' slopes alone miss such a turn wherever
          // the fold sits near an end, so the bracket is walked in
          // kRayleighStraddleSplit equal parts with the same Hermite test
          // on each: the first part that straddles zero, or a point at or
          // below zero at a part's minimum, ends the bracket. On the last
          // bracket (cut == true), a near-choke fold sits 1e-4 to 1e-2 below
          // hi, so also walk the 4^-j end grid if the root lands in the final
          // part. Its first fall is the bracket's (F(lo) > 0), with fewer
          // crossings after it. The whole bracket's root is kept where that
          // part holds it, so a bracket with no turn returns the same root.
          const double whole = brent(F, lo, hi, f_lo, f_hi, 0.0);
          const double end = hi, f_end = f_hi;
          const double part = (end - lo) / kRayleighStraddleSplit;
          const double lo0 = lo;
          bool split_hit = false;
          for (int i = 1; i <= kRayleighStraddleSplit; ++i) {
            if (i < kRayleighStraddleSplit) {
              hi = lo0 + i * part;
              f_hi = F(hi);
            } else {
              hi = end;
              f_hi = f_end;
            }
            double s_hi = 0.0;
            const double x = hermite_dip(s_hi);
            if (std::isfinite(x)) {
              hi = x;
              f_hi = F(x);
              split_hit = true;
              break;
            }
            if (f_hi <= 0.0) {
              split_hit = (i < kRayleighStraddleSplit);
              break;
            }
            lo = hi;
            f_lo = f_hi;
            s_lo = s_hi;
          }
          if (cut && !split_hit) {
            const double G_end = std::exp(end);
            auto x_of_u = [&](double u) { return std::log(G_end - std::exp(u)); };
            auto Fu = [&](double u) { return F(x_of_u(u)); };
            double xa = lo0, fa = F(lo0);
            double xb = lo, fb = f_lo;
            bool have_a = (xb > xa) && fa > 0.0;
            for (int j = kRayleighEndGridFirst; j <= kRayleighEndGridLast + 1; ++j) {
              const double xc =
                  (j > kRayleighEndGridLast) ? end : end + std::log1p(-std::pow(4.0, -j));
              if (!(xc > xb)) {
                continue;
              }
              const double fc = (j > kRayleighEndGridLast) ? f_end : F(xc);
              if (fc <= 0.0) {
                lo = xb;
                f_lo = fb;
                hi = xc;
                f_hi = fc;
                break;
              }
              if (have_a && fb < fa && fb <= fc) {
                const double ua = std::log(G_end - std::exp(xa));
                const double uc = (j > kRayleighEndGridLast)
                                      ? std::log(G_end) - (kRayleighEndGridLast + 1) * std::log(4.0)
                                      : std::log(G_end - std::exp(xc));
                const double ub = std::log(G_end - std::exp(xb));
                const double u = brent_min_to_zero(Fu, uc, ua, ub, kRayleighDipTol);
                if (std::isfinite(u)) {
                  lo = xa;
                  f_lo = fa;
                  hi = x_of_u(u);
                  f_hi = F(hi);
                  break;
                }
              }
              xa = xb;
              fa = fb;
              xb = xc;
              fb = fc;
              have_a = true;
            }
          }
          if (whole >= lo && whole <= hi) {
            root = whole;
          }
        }
        break;
      }
      if (!(may_turn && f_lo > 0.0 && f_hi > 0.0)) {
        s_lo = std::nan("");
        if (cut) {
          past_subsonic = f_hi > 0.0;
          break;
        }
        continue;
      }
      if (cut) {
        // The last bracket ends at the march's sonic flux, where the excess
        // varies on the scale of the distance to it, not of log G: a fold
        // there sits 1e-4 to 1e-2 of the flux below the end (the issue's
        // ducts). Sample at fluxes 4^-j of the end below it, in order, and
        // search each sampled minimum in u = log(G_end - G).
        const double G_end = std::exp(hi);
        auto x_of_u = [&](double u) { return std::log(G_end - std::exp(u)); };
        auto Fu = [&](double u) { return F(x_of_u(u)); };
        double xa = x_prev, fa = f_prev;  // two points back
        double xb = lo, fb = f_lo;        // one point back
        bool have_a = std::isfinite(xa) && fa > 0.0;
        for (int j = kRayleighEndGridFirst; j <= kRayleighEndGridLast + 1 && !found; ++j) {
          const double xc =
              (j > kRayleighEndGridLast) ? hi : hi + std::log1p(-std::pow(4.0, -j));
          if (!(xc > xb)) {
            continue;
          }
          const double fc = (j > kRayleighEndGridLast) ? f_hi : F(xc);
          if (fc <= 0.0) {
            lo = xb;
            f_lo = fb;
            hi = xc;
            f_hi = fc;
            found = true;
            break;
          }
          if (have_a && fb < fa && fb <= fc) {
            // At the end itself (xc = hi) u is -inf: bound the window at the
            // next grid distance, 4^-(Last+1) of the end, instead.
            const double ua = std::log(G_end - std::exp(xa));
            const double uc = (j > kRayleighEndGridLast)
                                  ? std::log(G_end) - (kRayleighEndGridLast + 1) * std::log(4.0)
                                  : std::log(G_end - std::exp(xc));
            const double ub = std::log(G_end - std::exp(xb));
            const double u = brent_min_to_zero(Fu, uc, ua, ub, kRayleighDipTol);
            if (std::isfinite(u)) {
              lo = xa;
              f_lo = fa;
              hi = x_of_u(u);
              f_hi = F(hi);
              found = true;
              break;
            }
          }
          xa = xb;
          fa = fb;
          xb = xc;
          fb = fc;
          have_a = true;
        }
        if (!found) {
          past_subsonic = true;
        }
        break;
      }
      // An interior bracket: the cubic Hermite interpolant of the ends'
      // values and slopes. An interior minimum of it (one turn, or a fall,
      // rise and fall) is searched, over [lo, the cubic's next maximum].
      double s_hi = 0.0;
      const double x = hermite_dip(s_hi);
      if (std::isfinite(x)) {
        hi = x;
        f_hi = F(x);
        found = true;
        break;
      }
      if (walk_both(f_lo, f_hi)) {
        // The cubic of the whole bracket shows no turn, but a dip narrower
        // than it near one end can hide from its ends' slopes (#1206): walk
        // it in kRayleighStraddleSplit equal parts with the same test. A part
        // end at or below zero, or a point at or below zero at a part's
        // minimum, ends the bracket below it; otherwise the scan goes on from
        // the whole bracket's end as before, so a bracket with no dip is
        // scanned to the same root.
        const double lo0 = lo, f_lo0 = f_lo;
        const double end = hi, f_end = f_hi, s_end = s_hi;
        const double part = (end - lo0) / kRayleighStraddleSplit;
        for (int i = 1; i <= kRayleighStraddleSplit && !found; ++i) {
          double s_i = s_end;
          if (i < kRayleighStraddleSplit) {
            hi = lo0 + i * part;
            f_hi = F(hi);
          } else {
            hi = end;
            f_hi = f_end;
          }
          if (f_hi <= 0.0) {
            found = true;
            break;
          }
          const double xi = hermite_dip(s_i, i == kRayleighStraddleSplit);
          if (std::isfinite(xi)) {
            hi = xi;
            f_hi = F(xi);
            found = true;
            break;
          }
          lo = hi;
          f_lo = f_hi;
          s_lo = s_i;
        }
        if (found) {
          break;
        }
        lo = lo0;
        f_lo = f_lo0;
        hi = end;
        f_hi = f_end;
      }
      s_lo = s_hi;
    }
    FirstFall out_ff;
    out_ff.found = found;
    out_ff.past_subsonic = past_subsonic;
    out_ff.lo = lo;
    out_ff.f_lo = f_lo;
    out_ff.hi = hi;
    out_ff.f_hi = f_hi;
    out_ff.root = root;
    return out_ff;
  };
  auto solve_fall = [&](const FirstFall& ff) {
    return std::exp(std::isfinite(ff.root) ? ff.root
                                           : brent(F, ff.lo, ff.hi, ff.f_lo, ff.f_hi, 0.0));
  };
  // A warm start: a tight bracket around the previous flux, widened. A
  // Newton iterate's flux is usually within 0.1% of the last one. On a
  // cooled duct, whose excess can turn, the bracket stops widening at 6.4% of
  // the guess (kRayleighWarmTurnSpan): a wider one could hold a fall, a rise
  // and a later fall, and Brent could take the later one (#1167); the cold
  // scan, which looks for the first, decides those.
  if (pb.G_guess > 0.0 && pb.G_guess < G_max && pb.P > 0.0) {
    const double w_max = may_turn ? kRayleighWarmTurnSpan : 1.0;
    double w = 1e-3;
    double lo = std::log(pb.G_guess * (1.0 - w));
    double hi = std::log(std::min(G_max * (1.0 - 1e-12), pb.G_guess * (1.0 + w)));
    double f_lo = F(lo), f_hi = F(hi);
    while ((f_lo <= 0.0 || f_hi > 0.0) && w < w_max) {
      w *= 8.0;
      if (f_lo <= 0.0) {
        hi = lo;
        f_hi = f_lo;
        lo = std::log(pb.G_guess * std::max(1.0 - w, 1e-6));
        f_lo = F(lo);
      } else {
        lo = hi;
        f_lo = f_hi;
        hi = std::log(std::min(G_max * (1.0 - 1e-12), pb.G_guess * (1.0 + w)));
        f_hi = F(hi);
      }
    }
    // Only a bracket whose upper end is itself unchoked holds a plain root;
    // one that straddles the choke goes the cold way below. On a cooled duct
    // an end whose excess is within the march's noise of zero (both crossings
    // of a fold's dip are roots there) has no sign, and Brent converged on a
    // noise crossing next to it: a warm bracket ending on the dip's rising
    // crossing returned that crossing (#1221). Such a bracket goes the cold
    // way too.
    const double noise = kRayleighExcessNoise * pb.Pt0;
    const bool ends_clear = !may_turn || (f_lo > noise && f_hi < -noise);
    if (f_lo > 0.0 && f_hi <= 0.0 && ends_clear &&
        march_at(std::exp(hi)).zeta_min > kRayleighZetaFloor) {
      // A cooled duct's excess can fall, rise and fall again, and a guess
      // just above a later fall brackets that one (#1191): the warm root then
      // depended on the solver's history. Its root is checked against the
      // cold one by the cold scan itself, run up to the bracket's lower end:
      // a fall below it is the first fall, and the root. Its brackets are
      // placed from lo, not from the choke; both scans walk their brackets
      // in parts (#1190, #1206), so a fold that hides between one scan's
      // bracket ends is not left to the other's, and the check matches the
      // cold root to the root's rounding, which within ~1e-5 of a fold's
      // peak drive is up to 1e-3 (and the derivatives' noise below ~1 Pa).
      // Without a fall below it the warm bracket holds the first fall.
      if (may_turn) {
        const FirstFall ff = first_fall(lo);
        if (ff.found) {
          G = solve_fall(ff);
          done = true;
        }
      }
      if (!done) {
        G = std::exp(brent(F, lo, hi, f_lo, f_hi, 1e-13));
        done = true;
      }
    }
  }
  if (!done) {
    const Choke ck = choked_flux(pb, pb.Pt0, pb.Tt0, pb.Tt_exit, X);
    mode = ck.mode;
    const double Gc = ck.G;
    const MarchEnd mc = march_at(Gc);
    const double Fc =
        (pb.P > 0.0)
            ? rayleigh_exit_pressure(pb.duct(Gc, pb.Tt0, pb.Tt_exit, &X), mc, pb.exit_total) - pb.P
            : 1.0;
    const FirstFall ff = first_fall(std::log(Gc));
    if (ff.found) {
      G = solve_fall(ff);
    } else {
      choked = Fc >= 0.0 || pb.P <= 0.0 || ff.past_subsonic;
      G = Gc;
    }
  }
  const MarchEnd m = march_at(G);
  out.G = G;
  out.choked = choked;
  out.M_in = m.M_in;
  out.M_exit = choked ? 1.0 : exit_mach(m);
  if (!pb.with_derivatives) {
    return out;
  }
  // Exact discrete tangents through the fixed 32-step RK4 march (#1276) via
  // DualN<4> seeded in (G, Pt0, Tt0, Tt_exit), replacing the internal
  // Richardson/central finite differences:
  //   unchoked  E(G; Pt0, Tt0, Tt_exit, Y) - P_target = 0 (one forward dual march)
  //   kInlet    G = fanno_sonic_mass_flux(Pt0, Tt0, Y) (analytic, zero marches)
  //   kExit     Psi = back_theta_inlet - inlet_theta = 0 (one backward dual march)
  //   kThroat   G* re-solve (interior throat relation goes as (G* - G)^(2/3)
  //             and clamps past sonic at RK4 stages, so G* itself is differenced).
  if (choked && mode == ChokeMode::kInlet) {
    out.dG_dPt0 = G / pb.Pt0;
    out.dG_dTt0 = isentropic_flux_slopes(pb.Pt0, pb.Tt0, 1.0, X).dTt;
    out.dG_dTt_exit = 0.0;
    out.dG_dP_target = 0.0;
    if (!pb.dirs->empty()) {
      out.dG_ddir.assign(pb.dirs->size(), 0.0);
      for (std::size_t i = 0; i < pb.dirs->size(); ++i) {
        out.dG_ddir[i] = isentropic_flux_t(pb.Pt0, DualT<1>::constant(pb.Tt0),
                                           DualT<1>::constant(1.0),
                                           mix_tangent<1>(X, &(*pb.dirs)[i], 0))
                             .d[0];
      }
    }
    return out;
  }

  if (choked && mode == ChokeMode::kThroat) {
    // A throat's relation goes as (G* - G)^(2/3) -- infinite slope at the
    // choke, and clamped past it -- so it cannot be differentiated implicitly
    // (1-36% errors). G* itself is smooth: re-solve it. Costly (two choke
    // solves per argument), but only a cooled duct with an interior sonic
    // point takes this branch.
    auto Gs = [&](double Pt0x, double Tt0x, double Ttex, const std::vector<double>& Xx) {
      return choked_flux(pb, Pt0x, Tt0x, Ttex, Xx).G;
    };
    const double ht = 1e-4;
    auto cdt = [&](auto&& g, double x0) {
      const double h = ht * std::abs(x0);
      return (g(x0 + h) - g(x0 - h)) / (2.0 * h);
    };
    out.dG_dPt0 = cdt([&](double v) { return Gs(v, pb.Tt0, pb.Tt_exit, X); }, pb.Pt0);
    out.dG_dTt0 = cdt([&](double v) { return Gs(pb.Pt0, v, pb.Tt_exit, X); }, pb.Tt0);
    out.dG_dTt_exit = cdt([&](double v) { return Gs(pb.Pt0, pb.Tt0, v, X); }, pb.Tt_exit);
    out.dG_dP_target = 0.0;
    if (!pb.dirs->empty()) {
      out.dG_ddir = d_ddir(X, *pb.dirs, [&](const std::vector<double>& Xk) {
        return Gs(pb.Pt0, pb.Tt0, pb.Tt_exit, Xk);
      });
    }
    return out;
  }

  const RayleighDuct d_root = pb.duct(G, pb.Tt0, pb.Tt_exit, &X);
  Dual4 rel_dual;
  if (choked) {
    // ChokeMode::kExit: Psi = inlet_w(Pt0, Tt0, G) - back_w_inlet(G, Tt0, Tt_exit) = 0.
    rel_dual = inlet_w_dual(pb.Pt0, pb.Tt0, G, X) - back_w_inlet_dual(d_root, pb.Pt0);
  } else {
    const Dual4 w_exit = march_w_exit_dual(d_root, pb.Pt0);
    rel_dual = rayleigh_exit_pressure_dual(d_root, w_exit, pb.exit_total);
  }
  const double R_G = rel_dual.d[0];
  out.dG_dPt0 = -rel_dual.d[1] / R_G;
  out.dG_dTt0 = -rel_dual.d[2] / R_G;
  out.dG_dTt_exit = -rel_dual.d[3] / R_G;
  out.dG_dP_target = choked ? 0.0 : 1.0 / R_G;
  if (pb.dirs->empty()) {
    return out;
  }

  // Composition directions: the exact tangent of the relation at fixed G
  // along v (#1310), over the exact discrete tangent R_G.
  out.dG_ddir.assign(pb.dirs->size(), 0.0);
  for (std::size_t i = 0; i < pb.dirs->size(); ++i) {
    const double R_v =
        rayleigh_relation_ddir(d_root, pb.Pt0, choked, pb.exit_total, (*pb.dirs)[i]);
    out.dG_ddir[i] = -R_v / R_G;
  }
  return out;
}

}  // namespace

FannoRayleighMarch fanno_rayleigh_march(double G, double Pt0, double Tt0, double Tt_exit,
                                        const std::vector<double>& X, double L, double D,
                                        double roughness, const std::string& friction_model,
                                        double f_multiplier) {
  check_duct(L, D);
  check_friction(roughness, f_multiplier);
  if (!validation::is_finite_positive(G) || !validation::is_finite_positive(Pt0) ||
      !validation::is_finite_positive(Tt0) || !validation::is_finite_positive(Tt_exit)) {
    throw std::invalid_argument("fanno_rayleigh_march: G, Pt0, Tt0 and Tt_exit must be positive");
  }
  const RayleighDuct d{G, Tt0, Tt_exit, L, D, roughness, f_multiplier, &X, &friction_model};
  const MarchEnd m = march_zeta(d, Pt0);
  FannoRayleighMarch out;
  out.M_in = m.M_in;
  out.zeta_min = m.zeta_min;
  out.zeta_peak = m.zeta_peak;
  out.M_exit = exit_mach(m);
  const FannoMachState s = fanno_state_at_mach(G, Tt_exit, out.M_exit, X);
  out.P_exit = s.P;
  out.T_exit = s.T;
  out.Pt_exit = G * exit_per_flux(out.M_exit, Tt_exit, X, true);
  return out;
}

double fanno_rayleigh_choked_mass_flux(double Pt0, double Tt0, double Tt_exit,
                                       const std::vector<double>& X, double L, double D,
                                       double roughness, const std::string& friction_model,
                                       double f_multiplier) {
  check_duct(L, D);
  check_friction(roughness, f_multiplier);
  if (!validation::is_finite(Pt0) || !validation::is_finite(Tt0) ||
      !validation::is_finite(Tt_exit)) {
    throw std::invalid_argument(
        "fanno_rayleigh_choked_mass_flux: Pt0, Tt0 and Tt_exit must be finite");
  }
  static const std::vector<std::vector<double>> kNone;
  const RayleighProblem pb{Pt0, Tt0, Tt_exit, 0.0, L, D, roughness, f_multiplier, &X,
                           &friction_model, false, -1.0, &kNone, false};
  return choked_flux(pb, Pt0, Tt0, Tt_exit, X).G;
}

FannoChokedAlternative fanno_rayleigh_choked_alternative(
    double Pt0, double Tt0, double Tt_exit, const std::vector<double>& X, double P_target,
    bool exit_total, double L, double D, double roughness, const std::string& friction_model,
    double f_multiplier) {
  check_duct(L, D);
  check_friction(roughness, f_multiplier);
  if (!validation::is_finite_positive(Pt0) || !validation::is_finite_positive(Tt0) ||
      !validation::is_finite_positive(Tt_exit) || !validation::is_finite(P_target)) {
    throw std::invalid_argument(
        "fanno_rayleigh_choked_alternative: Pt0, Tt0 and Tt_exit must be positive, P_target "
        "finite");
  }
  FannoChokedAlternative out;
  if (!(Pt0 > P_target)) {
    // No forward drive: fanno_rayleigh_channel_flow returns no flow, and no
    // forward choked state exists to compare against (nothing evaluated).
    out.exit_margin = std::numeric_limits<double>::quiet_NaN();
    return out;
  }
  static const std::vector<std::vector<double>> kNone;
  const RayleighProblem pb{Pt0, Tt0, Tt_exit, P_target, L, D, roughness, f_multiplier, &X,
                           &friction_model, exit_total, -1.0, &kNone, false};
  // The same exit pressure at the choke as the cold path's choke test in
  // rayleigh_raw_flow (its Fc), so "admissible" is the kernel's own verdict.
  out.G = choked_flux(pb, Pt0, Tt0, Tt_exit, X).G;
  const RayleighDuct d = pb.duct(out.G, Tt0, Tt_exit, &X);
  out.exit_margin = rayleigh_exit_pressure(d, march_zeta(d, Pt0), exit_total) - P_target;
  out.admissible = out.exit_margin >= 0.0;
  return out;
}

FannoChannelFlow fanno_rayleigh_channel_flow(
    double Pt0, double Tt0, double Tt_exit, const std::vector<double>& X, double P_target,
    bool exit_total, double L, double D, double roughness, const std::string& friction_model,
    double f_multiplier, double G_guess, const std::vector<std::vector<double>>& dY_directions,
    bool with_derivatives) {
  check_duct(L, D);
  check_friction(roughness, f_multiplier);
  FannoChannelFlow out;
  const double drive = Pt0 - P_target;
  static const std::vector<std::vector<double>> kNoDirections;
  const std::vector<std::vector<double>>& dirs = with_derivatives ? dY_directions : kNoDirections;
  if (!validation::is_finite(Pt0) || !validation::is_finite(Tt0) ||
      !validation::is_finite(Tt_exit) || !validation::is_finite(P_target)) {
    throw std::invalid_argument(
        "fanno_rayleigh_channel_flow: Pt0, Tt0, Tt_exit and P_target must be finite");
  }
  // As fanno_channel_flow: at exactly zero drive, the cubic's slope (#1209).
  if (!(drive >= 0.0) || !(Pt0 > 0.0) || !(Tt0 > 0.0) || !(Tt_exit > 0.0) ||
      (drive == 0.0 && !with_derivatives)) {
    out.dG_ddir.assign(dirs.size(), 0.0);
    return out;
  }
  RayleighProblem pb{Pt0, Tt0, Tt_exit, P_target, L, D, roughness, f_multiplier, &X,
                     &friction_model, exit_total, G_guess, &dirs, with_derivatives};
  if (drive >= kFannoFlowDriveFloor) {
    FannoChannelFlow raw = rayleigh_raw_flow(pb);
    raw.M_exit_root = raw.M_exit;
    return raw;
  }
  // Below the floor: same C1 low-Mach closure (or odd-cubic fallback) as
  // fanno_channel_flow (#1216).
  pb.P = Pt0 - kFannoFlowDriveFloor;
  pb.with_derivatives = true;
  const FannoChannelFlow at = rayleigh_raw_flow(pb);
  SigmaDerivs sig_d;
  if (with_derivatives && drive > 0.0) {
    sig_d = compute_floor_sigma_derivs(
        at,
        [&](double d_eval) {
          RayleighProblem pb_d = pb;
          pb_d.P = Pt0 - d_eval;
          pb_d.G_guess = at.G;
          pb_d.dirs = &kNoDirections;
          return rayleigh_raw_flow(pb_d);
        },
        [&](const std::vector<double>& Xk) {
          RayleighProblem pb_x = pb;
          pb_x.X = &Xk;
          pb_x.G_guess = at.G;
          pb_x.dirs = &kNoDirections;
          return rayleigh_raw_flow(pb_x);
        },
        X, dirs);
  }
  return apply_sub_floor_blend(at, sig_d, Pt0, Tt0, Tt_exit, X, drive, exit_total, L, D, roughness,
                               friction_model, f_multiplier, dirs, with_derivatives, false);
}

}  // namespace combaero
