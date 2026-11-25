#include "../include/equilibrium.h"
#include "../include/combustion.h"
#include "../include/inner_solve.h"
#include "thermo.h"
#include "composition.h"
#include "thermo_transport_data.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace combaero {

// -------------------------------------------------------------
// Internal utilities
// -------------------------------------------------------------

namespace {

// Internal WGS configuration (species indices)
struct WgsConfig {
  std::size_t i_CO;
  std::size_t i_H2O;
  std::size_t i_CO2;
  std::size_t i_H2;
};

static inline double clamp(double v, double lo, double hi) {
  return std::max(lo, std::min(v, hi));
}

static inline void composition_from_xi(const std::vector<double> &n0, double xi,
                                       const WgsConfig &cfg,
                                       std::vector<double> &n) {
  n = n0;
  n[cfg.i_CO] -= xi;
  n[cfg.i_H2O] -= xi;
  n[cfg.i_CO2] += xi;
  n[cfg.i_H2] += xi;
}

static inline void molefractions(const std::vector<double> &n,
                                 std::vector<double> &X) {
  double nt = 0.0;
  for (double v : n)
    nt += v;

  X.resize(n.size());
  if (nt <= 0.0) {
    double inv = 1.0 / static_cast<double>(n.size());
    for (std::size_t k = 0; k < n.size(); ++k)
      X[k] = inv;
  } else {
    for (std::size_t k = 0; k < n.size(); ++k)
      X[k] = n[k] / nt;
  }
}

// ln Kp of WGS, taken directly (exp then log overflows at a cold T).
static inline double lnKp_WGS(double T, const WgsConfig &cfg) {
  const double d_gRT = (g_over_RT(cfg.i_CO2, T) + g_over_RT(cfg.i_H2, T)) -
                       (g_over_RT(cfg.i_CO, T) + g_over_RT(cfg.i_H2O, T));
  return -d_gRT;
}

// Build WgsConfig from species names
static WgsConfig make_wgs_config() {
  WgsConfig cfg;
  cfg.i_CO = species_index.at("CO");
  cfg.i_H2O = species_index.at("H2O");
  cfg.i_CO2 = species_index.at("CO2");
  cfg.i_H2 = species_index.at("H2");
  return cfg;
}

// -------------------------------------------------------------
// Solve WGS equilibrium at fixed T
// -------------------------------------------------------------

// Extent xi of CO + H2O <-> CO2 + H2 at T from the amounts n0 (#927).
// F(xi) = ln Q - ln Kp with
//   ln Q = ln(n_CO2 + xi) + ln(n_H2 + xi) - ln(n_CO - xi) - ln(n_H2O - xi)
// (no change in moles, so the total cancels) rises from -inf at the lower
// bound (CO2 or H2 used up) to +inf at the upper (CO or H2O used up), with the
// analytic slope
//   dF/dxi = 1/n_CO2 + 1/n_H2 + 1/n_CO + 1/n_H2O.
// Safeguarded Newton (rtsafe) in the bracket, bisecting where a step leaves it
// or does not halve. The former damped Newton took a 1e-6 finite-difference
// step: with CO or H2O a trace the admissible interval is 3e-9..2e-7 wide, the
// stencil read past the bound (a 1e-300 floor made F jump there) and it
// stopped at ln(Q/K) = 12-14: where the shift runs away from the trace, a
// whole shift short (CO 7e-7 instead of 0.015 from CO 1e-8 at 800 K).
static double solve_wgs_extent(const std::vector<double> &n0, double T, const WgsConfig &cfg) {
  const double a_CO2 = n0[cfg.i_CO2];
  const double a_H2 = n0[cfg.i_H2];
  const double a_CO = n0[cfg.i_CO];
  const double a_H2O = n0[cfg.i_H2O];
  const double xi_min = std::max(-a_CO2, -a_H2);
  const double xi_max = std::min(a_CO, a_H2O);
  if (!(xi_max > xi_min)) {
    // Empty interval (a reactant and a product side both absent): nothing can
    // move, and 0 is the answer.
    return clamp(0.0, xi_min, xi_max);
  }
  const double lnKp = lnKp_WGS(T, cfg);
  if (!std::isfinite(lnKp)) {
    inner_solve::throw_no_convergence("solve_wgs_extent", 0, "ln Kp", lnKp);
  }
  const double inf = std::numeric_limits<double>::infinity();
  const double eps = std::numeric_limits<double>::epsilon();
  // F and dF/dxi at xi; noise is F's rounding error (each amount n0 +- xi is
  // rounded to eps (|n0| + |xi|), relative to the amount).
  auto eval = [&](double xi, double &dF, double &noise) {
    const double p_CO2 = a_CO2 + xi;
    const double p_H2 = a_H2 + xi;
    const double r_CO = a_CO - xi;
    const double r_H2O = a_H2O - xi;
    if (!(p_CO2 > 0.0 && p_H2 > 0.0)) return -inf;
    if (!(r_CO > 0.0 && r_H2O > 0.0)) return inf;
    dF = 1.0 / p_CO2 + 1.0 / p_H2 + 1.0 / r_CO + 1.0 / r_H2O;
    const double ax = std::abs(xi);
    noise = eps * ((std::abs(a_CO2) + ax) / p_CO2 + (std::abs(a_H2) + ax) / p_H2 +
                   (std::abs(a_CO) + ax) / r_CO + (std::abs(a_H2O) + ax) / r_H2O +
                   std::abs(lnKp));
    return std::log(p_CO2) + std::log(p_H2) - std::log(r_CO) - std::log(r_H2O) - lnKp;
  };
  double lo = xi_min;
  double hi = xi_max;
  double xi = (lo < 0.0 && 0.0 < hi) ? 0.0 : 0.5 * (lo + hi);
  double dx_old = hi - lo;
  double dx = dx_old;
  double F = 0.0;
  for (int it = 0; it < inner_solve::kWgsExtentMaxIter; ++it) {
    double dF = 0.0;
    double noise = 0.0;
    F = eval(xi, dF, noise);
    if (std::isnan(F)) {
      inner_solve::throw_no_convergence("solve_wgs_extent", it, "ln(Q/K)", F);
    }
    // A root to the tolerance, or to F's rounding (a trace amount at the root
    // carries only eps |n0| / n of relative precision).
    if (std::abs(F) <= EQUILIBRIUM_WGS_FTOL ||
        std::abs(F) <= EQUILIBRIUM_WGS_NOISE_FACTOR * noise) {
      return xi;
    }
    if (F > 0.0) {
      hi = xi;
    } else {
      lo = xi;
    }
    double xi_new = 0.0;
    const bool newton_ok = std::isfinite(F) && dF > 0.0 && std::isfinite(dF);
    const double step = newton_ok ? -F / dF : 0.0;
    if (newton_ok && xi + step > lo && xi + step < hi &&
        std::abs(2.0 * step) <= std::abs(dx_old)) {
      dx_old = dx;
      dx = step;
      xi_new = xi + step;
    } else {
      dx_old = dx;
      dx = 0.5 * (hi - lo);
      xi_new = lo + dx;
    }
    // The bracket cannot be split further: a root to rounding.
    if (!(xi_new > lo && xi_new < hi)) {
      return xi;
    }
    xi = xi_new;
  }
  if (inner_solve::cap_exit_is_root(hi - lo, xi)) {
    return xi;
  }
  inner_solve::throw_no_convergence("solve_wgs_extent", inner_solve::kWgsExtentMaxIter,
                                    "ln(Q/K)", F);
}

// -------------------------------------------------------------
// Adiabatic temperature
// -------------------------------------------------------------

// Root of the energy residual F(T) = H(T, X_eq(T)) - H_target on
// [EQUILIBRIUM_T_MIN, EQUILIBRIUM_T_MAX] (#927). F rises with T (the
// equilibrium heat capacity is positive), so it is bracketed from T_guess --
// the first step is the frozen-composition Newton step, which the equilibrium
// shift makes an overshoot, then twice the secant estimate -- and solved by Brent's method to a
// bracket of EQUILIBRIUM_T_XTOL_REL T. The former damped Newton with a
// finite-difference slope cycled on the S-shaped F of a rich reforming case
// (1500 K burned, the root at 924 K) and stopped at 1390 K after 50
// iterations, reported converged.
//
// Fun provides f(T, conv) (the residual [J/kg] or any consistent unit, conv =
// whether the inner equilibrium converged) and slope_frozen() (dF/dT at the
// composition of its last f). conv_out is false if an inner solve at either
// end of the final bracket did not converge.
template <typename Fun>
static double solve_adiabatic_T(const Fun &Fn, double T_guess, bool &conv_out) {
  const char *fn = "solve_adiabatic_T";
  auto eval = [&](double T, bool &c) {
    const double f = Fn.f(T, c);
    if (!std::isfinite(f)) {
      inner_solve::throw_no_convergence(fn, 0, "energy residual", f);
    }
    return f;
  };
  double a = clamp(T_guess, EQUILIBRIUM_T_MIN, EQUILIBRIUM_T_MAX);
  bool ca = true;
  double fa = eval(a, ca);
  conv_out = ca;
  if (fa == 0.0) return a;
  double step = std::abs(fa) / Fn.slope_frozen();
  if (!(std::isfinite(step) && step > 0.0)) step = EQUILIBRIUM_T_BRACKET_STEP;
  const double dir = fa > 0.0 ? -1.0 : 1.0;
  double b = a;
  double fb = fa;
  bool cb = ca;
  for (int k = 0;; ++k) {
    b = clamp(a + dir * step, EQUILIBRIUM_T_MIN, EQUILIBRIUM_T_MAX);
    fb = eval(b, cb);
    if (fb == 0.0) {
      conv_out = cb;
      return b;
    }
    if ((fa > 0.0) != (fb > 0.0)) break;
    if (b == EQUILIBRIUM_T_MIN || b == EQUILIBRIUM_T_MAX || k >= inner_solve::kAdiabaticMaxIter) {
      // The root at the limit itself, to the T tolerance (no fuel at
      // EQUILIBRIUM_T_MIN: F is rounding there).
      if (std::abs(fb) <= Fn.slope_frozen() * EQUILIBRIUM_T_XTOL_REL * b) {
        conv_out = cb;
        return b;
      }
      inner_solve::throw_no_convergence(fn, k + 1, "energy residual at the T limit", fb);
    }
    // Next step: twice the secant estimate of the remaining distance (an
    // overshoot that brackets the root when the estimate is good); doubling
    // if the secant is not usable. kAdiabaticMaxIter bounds the search.
    const double secant = std::abs(fb) * std::abs(b - a) / std::abs(fa - fb);
    const bool secant_ok = std::isfinite(secant) && std::abs(fb) < std::abs(fa);
    step = secant_ok ? 2.0 * secant : 2.0 * std::abs(b - a);
    if (!(step > 0.0)) step = EQUILIBRIUM_T_BRACKET_STEP;
    a = b;
    fa = fb;
    ca = cb;
  }
  // Brent (1973), as zbrent in Numerical Recipes: b is the best iterate, c
  // keeps the bracket [b, c].
  double c = a;
  double fc = fa;
  bool cc = ca;
  double d = b - a;
  double e = d;
  const double eps = std::numeric_limits<double>::epsilon();
  for (int it = 0; it < inner_solve::kAdiabaticMaxIter; ++it) {
    if ((fb > 0.0) == (fc > 0.0)) {
      c = a;
      fc = fa;
      cc = ca;
      d = b - a;
      e = d;
    }
    if (std::abs(fc) < std::abs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
      std::swap(cb, cc);
      ca = cc;
    }
    const double tol1 = 2.0 * eps * std::abs(b) + 0.5 * EQUILIBRIUM_T_XTOL_REL * std::abs(b);
    const double xm = 0.5 * (c - b);
    if (std::abs(xm) <= tol1 || fb == 0.0) {
      conv_out = cb && cc;
      return b;
    }
    if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
      // Inverse quadratic interpolation (secant when a == c).
      const double s = fb / fa;
      double p = 0.0;
      double q = 0.0;
      if (a == c) {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      } else {
        const double qq = fa / fc;
        const double r = fb / fc;
        p = s * (2.0 * xm * qq * (qq - r) - (b - a) * (r - 1.0));
        q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0) q = -q;
      p = std::abs(p);
      const double min1 = 3.0 * xm * q - std::abs(tol1 * q);
      const double min2 = std::abs(e * q);
      if (2.0 * p < std::min(min1, min2)) {
        e = d;
        d = p / q;
      } else {
        d = xm;
        e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    ca = cb;
    b += std::abs(d) > tol1 ? d : (xm > 0.0 ? tol1 : -tol1);
    fb = eval(b, cb);
  }
  inner_solve::throw_no_convergence(fn, inner_solve::kAdiabaticMaxIter, "bracket [K]", c - b);
}

struct AdiabaticFun {
  const std::vector<double> &n_in;
  double H_mass_target;
  WgsConfig cfg;

  mutable std::vector<double> n;
  mutable std::vector<double> X;
  mutable double T_last = 0.0;

  AdiabaticFun(const std::vector<double> &n_in_, double H_mass_target_,
               const WgsConfig &cfg_)
      : n_in(n_in_), H_mass_target(H_mass_target_), cfg(cfg_), n(n_in_.size()),
        X(n_in_.size()) {}

  double f(double T, bool &conv) const {
    // Inner: solve xi(T) (throws if it does not converge)
    const double xi = solve_wgs_extent(n_in, T, cfg);
    conv = true;
    composition_from_xi(n_in, xi, cfg, n);
    molefractions(n, X);
    T_last = T;
    return h(T, X) / mwmix(X) - H_mass_target;
  }

  double slope_frozen() const { return cp(T_last, X) / mwmix(X); }
};

// -------------------------------------------------------------
// General Steam Reforming + WGS configuration
// -------------------------------------------------------------
// Handles all hydrocarbons CnHm via:
//   CnHm + n*H2O -> n*CO + (n + m/2)*H2
// Plus WGS:
//   CO + H2O <-> CO2 + H2

struct ReformingConfig {
  std::size_t i_H2O;
  std::size_t i_CO;
  std::size_t i_CO2;
  std::size_t i_H2;

  // Hydrocarbon species indices and their C, H counts
  std::vector<std::size_t> hc_indices;
  std::vector<int> hc_C; // number of C atoms
  std::vector<int> hc_H; // number of H atoms
};

static ReformingConfig make_reforming_config() {
  ReformingConfig cfg;
  cfg.i_H2O = species_index.at("H2O");
  cfg.i_CO = species_index.at("CO");
  cfg.i_CO2 = species_index.at("CO2");
  cfg.i_H2 = species_index.at("H2");

  // Find all hydrocarbon species (C > 0, H > 0, O = 0, N = 0)
  for (std::size_t i = 0; i < species_names.size(); ++i) {
    const auto &ms = molecular_structures[i];
    if (ms.C > 0 && ms.H > 0 && ms.O == 0 && ms.N == 0) {
      cfg.hc_indices.push_back(i);
      cfg.hc_C.push_back(static_cast<int>(ms.C));
      cfg.hc_H.push_back(static_cast<int>(ms.H));
    }
  }
  return cfg;
}

} // anonymous namespace

// -------------------------------------------------------------
// State-based WGS equilibrium functions (Public API)
// -------------------------------------------------------------

// WGS equilibrium (isothermal) - equilibrate at input temperature
EquilibriumResult wgs_equilibrium(const State &in) {
  WgsConfig cfg = make_wgs_config();

  // Solve for extent of reaction at fixed T
  const double xi = solve_wgs_extent(in.X, in.T, cfg);

  // Build output composition
  std::vector<double> n_out;
  composition_from_xi(in.X, xi, cfg, n_out);

  std::vector<double> X_out;
  molefractions(n_out, X_out);

  State out;
  out.set_TPX(in.T, in.P, X_out);
  return EquilibriumResult{out, in.T, 0.0, true};
}

// WGS equilibrium (adiabatic) - solve for equilibrium T and composition
EquilibriumResult wgs_equilibrium_adiabatic(const State &in) {
  WgsConfig cfg = make_wgs_config();

  // Target enthalpy (conserved)
  double H_mass_target = h(in.T, in.X) / mwmix(in.X);

  // Solve for adiabatic temperature
  AdiabaticFun F{in.X, H_mass_target, cfg};
  bool conv = true;
  const double T_ad = solve_adiabatic_T(F, in.T, conv);

  // Get equilibrium composition at T_ad
  State at_T_ad;
  at_T_ad.set_TPX(T_ad, in.P, in.X);

  EquilibriumResult inner = wgs_equilibrium(at_T_ad);
  return EquilibriumResult{inner.state, in.T, inner.state.T - in.T, conv};
}

// -------------------------------------------------------------
// General Reforming + WGS equilibrium (handles all hydrocarbons)
// -------------------------------------------------------------

// Equilibrium constant for general steam reforming: CnHm + n*H2O <-> n*CO +
// (n+m/2)*H2 Kp = (y_CO^n * y_H2^(n+m/2)) / (y_CnHm * y_H2O^n) *
// (P/P0)^(delta_n) where delta_n = n + (n+m/2) - 1 - n = m/2 - 1 + n = (m-2)/2
// + n
// (returned as ln Kp: exp then log overflows at a cold T)
static double lnKp_reforming(double T, std::size_t i_hc, int n_C, int n_H,
                             const ReformingConfig &cfg) {
  double g_hc = g_over_RT(i_hc, T);
  double g_H2O = g_over_RT(cfg.i_H2O, T);
  double g_CO = g_over_RT(cfg.i_CO, T);
  double g_H2 = g_over_RT(cfg.i_H2, T);

  // CnHm + n*H2O -> n*CO + (n + m/2)*H2
  double n_H2_prod = n_C + 0.5 * n_H; // n + m/2

  // Delta_G/RT = n*g_CO + (n+m/2)*g_H2 - g_CnHm - n*g_H2O
  double d_gRT = n_C * g_CO + n_H2_prod * g_H2 - g_hc - n_C * g_H2O;
  return -d_gRT;
}

// Solve single hydrocarbon reforming equilibrium at fixed T
// CnHm + n*H2O <-> n*CO + (n+m/2)*H2
static double solve_reforming_isothermal(const std::vector<double> &n0,
                                         double T, double P, std::size_t i_hc,
                                         int n_C, int n_H,
                                         const ReformingConfig &cfg) {
  const double min_species = 1e-12;
  const double m_H2 = n_C + 0.5 * n_H;

  // The extent is bracketed both ways: forward (reforming) until the
  // hydrocarbon or the steam runs out, backward (methanation, the hydrocarbon
  // re-forming from CO + H2) until CO or H2 does. Forward only left a mixture
  // of higher hydrocarbons at 600-800 K 0.2 off in mole fraction: at
  // equilibrium they become methane, which is CH4 reforming run backward.
  const double xi_hi = std::min(n0[i_hc], n0[cfg.i_H2O] / n_C);
  const double xi_lo = -std::min(n0[cfg.i_CO] / n_C, n0[cfg.i_H2] / m_H2);
  if (xi_hi - xi_lo < min_species) {
    return 0.0; // Nothing can move either way
  }

  const double lnKp = lnKp_reforming(T, i_hc, n_C, n_H, cfg);

  // delta_n = products - reactants = (n + (n+m/2)) - (1 + n) = n + m/2 - 1
  double delta_n = n_C + 0.5 * n_H - 1.0;

  const double span = xi_hi - xi_lo;

  // Safeguarded Newton on F(xi) = ln Q - ln Kp over (xi_lo, xi_hi), with the
  // analytic slope (#716)
  //   dF/dxi = nC^2/n_CO + m^2/n_H2 + 1/n_hc + nC^2/n_H2O - dn^2/n_t,  m = nC + nH/2,
  // and bisection where a step leaves the bracket. F rises from -inf (CO or H2
  // used up) to +inf as the hydrocarbon or the steam runs out, so the root is
  // bracketed. The former finite-difference slope used an absolute step of
  // 1e-8 mol: near complete reforming (1e-11 mol left) it stepped past the
  // hydrocarbon and stalled off the root.
  double nt0 = 0.0;
  for (double v : n0)
    nt0 += v;
  auto eval = [&](double xi, double &dF) {
    const double n_hc = n0[i_hc] - xi;
    const double n_H2O = n0[cfg.i_H2O] - n_C * xi;
    const double n_CO = n0[cfg.i_CO] + n_C * xi;
    const double n_H2 = n0[cfg.i_H2] + m_H2 * xi;
    const double nt = nt0 + delta_n * xi;
    // An amount rounded to zero at a bracket end: F's limit there (#927).
    if (!(n_CO > 0.0 && n_H2 > 0.0)) return -std::numeric_limits<double>::infinity();
    if (!(n_hc > 0.0 && n_H2O > 0.0)) return std::numeric_limits<double>::infinity();
    const double lnQ = n_C * std::log(n_CO) + m_H2 * std::log(n_H2) - std::log(n_hc) -
                       n_C * std::log(n_H2O) - delta_n * std::log(nt) +
                       delta_n * std::log(P / EQUILIBRIUM_P_REF);
    dF = n_C * n_C / n_CO + m_H2 * m_H2 / n_H2 + 1.0 / n_hc + n_C * n_C / n_H2O -
         delta_n * delta_n / nt;
    return lnQ - lnKp;
  };
  double lo = xi_lo;
  double hi = xi_hi;
  double xi = 0.5 * (xi_lo + xi_hi);
  double F = 0.0;
  for (int it = 0; it < inner_solve::kReformingExtentMaxIter; ++it) {
    double dF = 0.0;
    F = eval(xi, dF);
    if (std::isnan(F)) {
      // Formerly a silent break to the clamped extent (#868 R9, E2).
      inner_solve::throw_no_convergence("solve_reforming_isothermal", it, "ln(Q/K)", F);
    }
    if (std::abs(F) < 1e-12) {
      return xi;
    }
    if (F > 0.0) {
      hi = xi;
    } else {
      lo = xi;
    }
    double xi_new =
        (std::isfinite(F) && dF > 0.0 && std::isfinite(dF)) ? xi - F / dF : 0.5 * (lo + hi);
    if (!(xi_new > lo && xi_new < hi)) {
      xi_new = 0.5 * (lo + hi);
    }
    if (hi - lo <= 1e-15 * span) {
      return xi_new;
    }
    xi = xi_new;
  }
  if (inner_solve::cap_exit_is_root(hi - lo, xi)) {
    return xi;
  }
  inner_solve::throw_no_convergence("solve_reforming_isothermal",
                                    inner_solve::kReformingExtentMaxIter, "ln(Q/K)", F);
}

// Apply reforming extent to composition
static void apply_reforming(std::vector<double> &n, double xi, std::size_t i_hc,
                            int n_C, int n_H, const ReformingConfig &cfg) {
  n[i_hc] -= xi;
  n[cfg.i_H2O] -= n_C * xi;
  n[cfg.i_CO] += n_C * xi;
  n[cfg.i_H2] += (n_C + 0.5 * n_H) * xi;
}

// One round: reform each hydrocarbon in turn, then WGS, each to its own
// equilibrium at the current composition.
static void reforming_wgs_round(std::vector<double> &n_out, double T, double P,
                                const ReformingConfig &cfg) {
  for (std::size_t k = 0; k < cfg.hc_indices.size(); ++k) {
    std::size_t i_hc = cfg.hc_indices[k];
    int n_C = cfg.hc_C[k];
    int n_H = cfg.hc_H[k];

    double xi = solve_reforming_isothermal(n_out, T, P, i_hc, n_C, n_H, cfg);
    apply_reforming(n_out, xi, i_hc, n_C, n_H, cfg);
  }

  WgsConfig wgs_cfg = make_wgs_config();

  double xi_min = std::max(-n_out[wgs_cfg.i_CO2], -n_out[wgs_cfg.i_H2]);
  double xi_max = std::min(n_out[wgs_cfg.i_CO], n_out[wgs_cfg.i_H2O]);

  if (xi_max > xi_min + 1e-15) {
    const double xi_wgs = solve_wgs_extent(n_out, T, wgs_cfg);

    n_out[wgs_cfg.i_CO] -= xi_wgs;
    n_out[wgs_cfg.i_H2O] -= xi_wgs;
    n_out[wgs_cfg.i_CO2] += xi_wgs;
    n_out[wgs_cfg.i_H2] += xi_wgs;
  }
}

// Largest change of any amount between two rounds, relative to the total.
static double reforming_round_change(const std::vector<double> &before,
                                     const std::vector<double> &after) {
  double total = 0.0, d = 0.0;
  for (std::size_t j = 0; j < after.size(); ++j) {
    total += after[j];
    d = std::max(d, std::abs(after[j] - before[j]));
  }
  return d / std::max(total, 1e-300);
}

// -------------------------------------------------------------
// Joint reforming + WGS equilibrium (#759)
// -------------------------------------------------------------
// One extent per reaction (steam reforming of each hydrocarbon, then WGS),
// solved at once: F_r(xi) = ln Q_r - ln K_r = 0, with
//   ln Q_r = sum_j nu_rj ln n_j - dn_r ln n_t + dn_r ln(P/P0),
//   dF_r/dxi_s = sum_j nu_rj nu_sj / n_j - dn_r dn_s / n_t,
// the Hessian of the Gibbs function in the extents (symmetric, positive
// definite while an inert is present). Rounds of scalar solves only provide
// the start: once the steam is used up, reverse WGS returns it a sliver per
// round and the rounds crawl (rich propane, phi >= 2.5).

struct JointReaction {
  std::vector<std::pair<std::size_t, double>> nu;  // (species, coefficient)
  double dn = 0.0;                                 // change in total moles
  std::size_t key = 0;  // entry of nu whose species no other reaction has
};

// The reactions: reforming of each hydrocarbon (if with_reforming), then WGS.
// to_co2: reform to CO2 (CnHm + 2n H2O = n CO2 + (2n + m/2) H2) and key WGS
// on CO instead of CO2, the same reaction space (#1262). Where CO is the
// scarcer of the two (a lean burn with trace fuel: CO at 4e-12 beside CO2 at
// 4e-4), it is then a key and steps in ln n; as a non-key it was left as the
// difference of the hydrocarbon keys' changes, each quantised to its own ulp,
// a ln Q noise of 1.4e-10 at phi 1e-4, 300 K.
static std::vector<JointReaction> joint_reactions(const ReformingConfig &cfg,
                                                  bool with_reforming, bool to_co2 = false) {
  std::vector<JointReaction> rx;
  if (with_reforming) {
    for (std::size_t k = 0; k < cfg.hc_indices.size(); ++k) {
      const double nC = cfg.hc_C[k];
      const double nH = cfg.hc_H[k];
      if (to_co2) {
        rx.push_back(JointReaction{{{cfg.hc_indices[k], -1.0},
                                    {cfg.i_H2O, -2.0 * nC},
                                    {cfg.i_CO2, nC},
                                    {cfg.i_H2, 2.0 * nC + 0.5 * nH}},
                                   nC + 0.5 * nH - 1.0, 0});
      } else {
        const double m = nC + 0.5 * nH;
        rx.push_back(JointReaction{{{cfg.hc_indices[k], -1.0},
                                    {cfg.i_H2O, -nC},
                                    {cfg.i_CO, nC},
                                    {cfg.i_H2, m}},
                                   m - 1.0, 0});
      }
    }
  }
  rx.push_back(JointReaction{
      {{cfg.i_CO, -1.0}, {cfg.i_H2O, -1.0}, {cfg.i_CO2, 1.0}, {cfg.i_H2, 1.0}}, 0.0,
      to_co2 ? std::size_t{0} : std::size_t{2}});
  return rx;
}

// The same reaction space in a basis for a hydrocarbon-rich mix (#1262).
// With R = CrHs the most abundant hydrocarbon: two reactions of R among the
// oxygen carriers H2O, CO and CO2,
//   R + r H2O  = r CO + (r + s/2) H2      (steam),
//   R + r CO2  = 2r CO + (s/2) H2          (dry),
//   R + 2r H2O = r CO2 + (2r + s/2) H2     (steam to CO2),
// the pair that shares the most abundant carrier (each other carrier the key
// of its own), and each other hydrocarbon's exchange with R,
//   CnHm + x H2 = (n/r) R,  x = (n s/r - m)/2         (key CnHm).
// Every set spans steam reforming of each hydrocarbon and WGS, so the
// equilibrium is the same. With the oxygen nearly gone (a fuel-only or very
// rich burn) every steam-reforming row carries the trace carriers with weight
// 1/n, and the hydrocarbon interconversions they span were left as
// differences of those rows: the scaled Hessian's condition number reached
// 2.5e11 (eigenvalue 3e-11 along CH4 - C2H6), and each extent update
// cancelled 1e-3 mol into 1e-12 mol of steam, a ln Q noise of 1e-9 Newton
// could not get under. Here the interconversions hold no oxygen species, and
// the two scarcer carriers step in ln n as keys, so none of them is left as
// the difference of larger changes. R itself is no key: false (rx unchanged)
// unless it is at least REFORMING_RICH_BASIS_HC_FRACTION of the total.
static bool joint_reactions_rich(const ReformingConfig &cfg, const std::vector<double> &n,
                                 std::vector<JointReaction> &rx) {
  if (cfg.hc_indices.empty()) return false;
  double nt = 0.0;
  for (double v : n) nt += v;
  std::size_t k_ref = 0;
  for (std::size_t k = 1; k < cfg.hc_indices.size(); ++k) {
    if (n[cfg.hc_indices[k]] > n[cfg.hc_indices[k_ref]]) k_ref = k;
  }
  const std::size_t i_ref = cfg.hc_indices[k_ref];
  if (!(n[i_ref] >= REFORMING_RICH_BASIS_HC_FRACTION * nt)) return false;
  const double r = cfg.hc_C[k_ref];
  const double sH = cfg.hc_H[k_ref];
  rx.clear();
  for (std::size_t k = 0; k < cfg.hc_indices.size(); ++k) {
    if (k == k_ref) continue;
    const double nC = cfg.hc_C[k];
    const double x = 0.5 * (nC * sH / r - cfg.hc_H[k]);
    rx.push_back(JointReaction{{{cfg.hc_indices[k], -1.0}, {cfg.i_H2, -x}, {i_ref, nC / r}},
                               nC / r - 1.0 - x, 0});
  }
  // nu entries ordered so that entry 0 is the key.
  const double dn = r + 0.5 * sH - 1.0;
  const JointReaction h2o_co_key_h2o{
      {{cfg.i_H2O, -r}, {i_ref, -1.0}, {cfg.i_CO, r}, {cfg.i_H2, r + 0.5 * sH}}, dn, 0};
  const JointReaction h2o_co_key_co{
      {{cfg.i_CO, r}, {i_ref, -1.0}, {cfg.i_H2O, -r}, {cfg.i_H2, r + 0.5 * sH}}, dn, 0};
  const JointReaction co2_co_key_co2{
      {{cfg.i_CO2, -r}, {i_ref, -1.0}, {cfg.i_CO, 2.0 * r}, {cfg.i_H2, 0.5 * sH}}, dn, 0};
  const JointReaction co2_co_key_co{
      {{cfg.i_CO, 2.0 * r}, {i_ref, -1.0}, {cfg.i_CO2, -r}, {cfg.i_H2, 0.5 * sH}}, dn, 0};
  const JointReaction h2o_co2_key_co2{
      {{cfg.i_CO2, r}, {i_ref, -1.0}, {cfg.i_H2O, -2.0 * r}, {cfg.i_H2, 2.0 * r + 0.5 * sH}},
      dn, 0};
  const JointReaction h2o_co2_key_h2o{
      {{cfg.i_H2O, -2.0 * r}, {i_ref, -1.0}, {cfg.i_CO2, r}, {cfg.i_H2, 2.0 * r + 0.5 * sH}},
      dn, 0};
  const double a_h2o = n[cfg.i_H2O], a_co = n[cfg.i_CO], a_co2 = n[cfg.i_CO2];
  if (a_co >= a_h2o && a_co >= a_co2) {
    rx.push_back(h2o_co_key_h2o);
    rx.push_back(co2_co_key_co2);
  } else if (a_h2o >= a_co2) {
    rx.push_back(h2o_co_key_co);
    rx.push_back(h2o_co2_key_co2);
  } else {
    rx.push_back(co2_co_key_co);
    rx.push_back(h2o_co2_key_h2o);
  }
  // H2 is in every reaction above and the key of none, so its change is the
  // sum of all the extents'. Where it is scarcer than the second hydrocarbon
  // (the pyrolysis of cold heptane at phi 1000, 10 bar: H2 at 4e-12 beside
  // hexane at 0.12), the linear solve cannot place it to better than a ln Q
  // of up to 7e-9 in the oxygen carriers' rows (15 H2 in each), above
  // REFORMING_NEWTON_FLOOR_FTOL, and half the states in 505-525 K were
  // refused (#1262 CI). There, the second most abundant hydrocarbon joins
  // the first as a reference: every other hydrocarbon, and H2 itself, is
  // exchanged with the two (C and H balanced by them), and H2 is substituted
  // out of the oxygen carriers' rows. H2 then keys its own reaction only; the
  // same equilibrium.
  std::size_t k2 = cfg.hc_indices.size();
  for (std::size_t k = 0; k < cfg.hc_indices.size(); ++k) {
    if (k == k_ref) continue;
    if (k2 == cfg.hc_indices.size() || n[cfg.hc_indices[k]] > n[cfg.hc_indices[k2]]) k2 = k;
  }
  if (k2 == cfg.hc_indices.size()) return true;
  const std::size_t i_2 = cfg.hc_indices[k2];
  const double c2 = cfg.hc_C[k2];
  const double h2 = cfg.hc_H[k2];
  const double det = r * h2 - c2 * sH;
  if (!(n[cfg.i_H2] < n[i_2]) || !(std::abs(det) > 0.0)) return true;
  // (c, h) = a (r, sH) + b (c2, h2).
  auto split = [&](double c, double h, double &a, double &b) {
    a = (c * h2 - c2 * h) / det;
    b = (r * h - c * sH) / det;
  };
  double a_h = 0.0, b_h = 0.0;
  split(0.0, 2.0, a_h, b_h);
  std::vector<JointReaction> two;
  for (std::size_t k = 0; k < cfg.hc_indices.size(); ++k) {
    if (k == k_ref || k == k2) continue;
    double a = 0.0, b = 0.0;
    split(cfg.hc_C[k], cfg.hc_H[k], a, b);
    two.push_back(
        JointReaction{{{cfg.hc_indices[k], -1.0}, {i_ref, a}, {i_2, b}}, a + b - 1.0, 0});
  }
  two.push_back(JointReaction{{{cfg.i_H2, -1.0}, {i_ref, a_h}, {i_2, b_h}}, a_h + b_h - 1.0, 0});
  for (std::size_t q = rx.size() - 2; q < rx.size(); ++q) {
    JointReaction o;
    o.key = 0;
    double nu_h2 = 0.0, nu_ref = 0.0;
    for (const auto &p : rx[q].nu) {
      if (p.first == cfg.i_H2) {
        nu_h2 = p.second;
      } else if (p.first == i_ref) {
        nu_ref = p.second;
      } else {
        o.nu.push_back(p);
      }
    }
    o.nu.push_back({i_ref, nu_ref + nu_h2 * a_h});
    o.nu.push_back({i_2, nu_h2 * b_h});
    o.dn = rx[q].dn + nu_h2 * (a_h + b_h - 1.0);
    two.push_back(o);
  }
  rx = two;
  return true;
}

static bool solve_dense(std::vector<double> A, std::vector<double> &B, std::size_t m,
                        std::size_t k);

// The same reaction space in a component basis (#1288; Smith and Missen's
// stoichiometric formulation): the most abundant species that are
// independent in their elements are the components, and every other species
// of the joint solve keys its own formation from them. Each non-key amount
// is n plus the extents' changes, and cancels where it is far below them;
// here the non-keys are the components, the most abundant species there can
// be, so every scarce species is a key and steps in ln n. False if the
// present species do not span the present elements.
static bool joint_reactions_components(const ReformingConfig &cfg,
                                       const std::vector<double> &n,
                                       std::vector<JointReaction> &rx) {
  std::vector<std::size_t> pool = cfg.hc_indices;
  for (std::size_t j : {cfg.i_H2O, cfg.i_CO, cfg.i_CO2, cfg.i_H2}) pool.push_back(j);
  auto formula = [](std::size_t j) {
    const auto &ms = molecular_structures[j];
    return std::vector<double>{static_cast<double>(ms.C), static_cast<double>(ms.H),
                               static_cast<double>(ms.O)};
  };
  double b[3] = {0.0, 0.0, 0.0};
  for (std::size_t j : pool) {
    const auto a = formula(j);
    for (std::size_t e = 0; e < 3; ++e) b[e] += a[e] * n[j];
  }
  std::vector<std::size_t> el;
  for (std::size_t e = 0; e < 3; ++e) {
    if (b[e] > 0.0) el.push_back(e);
  }
  const std::size_t E = el.size();
  if (E == 0) return false;
  auto reduced = [&](std::size_t j) {
    const auto a = formula(j);
    std::vector<double> v(E);
    for (std::size_t k = 0; k < E; ++k) v[k] = a[el[k]];
    return v;
  };
  auto spans_present = [&](std::size_t j) {
    const auto a = formula(j);
    for (std::size_t e = 0; e < 3; ++e) {
      if (a[e] > 0.0 && !(b[e] > 0.0)) return false;
    }
    return true;
  };
  std::vector<std::size_t> order;
  for (std::size_t j : pool) {
    if (n[j] > 0.0) order.push_back(j);
  }
  std::stable_sort(order.begin(), order.end(),
                   [&](std::size_t x, std::size_t y) { return n[x] > n[y]; });
  std::vector<std::size_t> comp;
  std::vector<std::vector<double>> red;
  for (std::size_t j : order) {
    if (comp.size() == E) break;
    std::vector<double> v = reduced(j);
    for (const auto &w : red) {
      std::size_t pk = 0;
      for (std::size_t k = 1; k < E; ++k) {
        if (std::abs(w[k]) > std::abs(w[pk])) pk = k;
      }
      const double f = v[pk] / w[pk];
      for (std::size_t k = 0; k < E; ++k) v[k] -= f * w[k];
    }
    double vmax = 0.0;
    for (double x : v) vmax = std::max(vmax, std::abs(x));
    if (vmax > REFORMING_EP_COMPONENT_RANK_TOL) {
      comp.push_back(j);
      red.push_back(v);
    }
  }
  if (comp.size() < E) return false;
  // C^T nu = a_j: column c of C^T is component c's formula.
  std::vector<double> Ct(E * E);
  for (std::size_t c = 0; c < E; ++c) {
    const auto v = reduced(comp[c]);
    for (std::size_t k = 0; k < E; ++k) Ct[k * E + c] = v[k];
  }
  rx.clear();
  for (std::size_t j : pool) {
    if (std::find(comp.begin(), comp.end(), j) != comp.end() || !spans_present(j)) continue;
    std::vector<double> nu = reduced(j);
    if (!solve_dense(Ct, nu, E, 1)) return false;
    JointReaction r;
    r.nu.push_back({j, -1.0});
    r.dn = -1.0;
    r.key = 0;
    for (std::size_t c = 0; c < E; ++c) {
      // Formulas are small integers: a coefficient within rounding of zero
      // is none (the others are kept as solved, fractions included, so the
      // reaction conserves the elements).
      const double v = nu[c];
      if (std::abs(v) <= REFORMING_EP_COMPONENT_RANK_TOL) continue;
      r.nu.push_back({comp[c], v});
      r.dn += v;
    }
    rx.push_back(r);
  }
  return true;
}

// The smallest present amount that is no reaction's key: the species that
// cancels first as the difference of the extents' changes.
static double joint_min_nonkey(const std::vector<JointReaction> &rx,
                               const std::vector<double> &n) {
  std::vector<bool> key(n.size(), false), in(n.size(), false);
  for (const auto &r : rx) {
    key[r.nu[r.key].first] = true;
    for (const auto &p : r.nu) in[p.first] = true;
  }
  double m = std::numeric_limits<double>::infinity();
  for (std::size_t j = 0; j < n.size(); ++j) {
    if (in[j] && !key[j] && n[j] > 0.0) m = std::min(m, n[j]);
  }
  return m;
}

// ln K_r(T) = -sum_j nu_rj g_j/RT.
static double joint_lnK(const JointReaction &r, double T) {
  double d = 0.0;
  for (const auto &p : r.nu) d += p.second * g_over_RT(p.first, T);
  return -d;
}

// Solve A x = B in place (A m x m, B m x k, row-major) by Gaussian elimination
// with partial pivoting. Returns false if A is singular.
static bool solve_dense(std::vector<double> A, std::vector<double> &B, std::size_t m,
                        std::size_t k) {
  for (std::size_t c = 0; c < m; ++c) {
    std::size_t piv = c;
    for (std::size_t r = c + 1; r < m; ++r) {
      if (std::abs(A[r * m + c]) > std::abs(A[piv * m + c])) piv = r;
    }
    if (!(std::abs(A[piv * m + c]) > 0.0)) return false;
    if (piv != c) {
      for (std::size_t q = 0; q < m; ++q) std::swap(A[c * m + q], A[piv * m + q]);
      for (std::size_t q = 0; q < k; ++q) std::swap(B[c * k + q], B[piv * k + q]);
    }
    for (std::size_t r = c + 1; r < m; ++r) {
      const double f = A[r * m + c] / A[c * m + c];
      if (f == 0.0) continue;
      for (std::size_t q = c; q < m; ++q) A[r * m + q] -= f * A[c * m + q];
      for (std::size_t q = 0; q < k; ++q) B[r * k + q] -= f * B[c * k + q];
    }
  }
  for (std::size_t c = m; c-- > 0;) {
    for (std::size_t q = 0; q < k; ++q) {
      double s = B[c * k + q];
      for (std::size_t r = c + 1; r < m; ++r) s -= A[c * m + r] * B[r * k + q];
      B[c * k + q] = s / A[c * m + c];
    }
  }
  return true;
}

// The reactions whose species are all present: their extent can move both
// ways, so ln Q is finite and the equilibrium is interior.
static std::vector<std::size_t> joint_active(const std::vector<JointReaction> &rx,
                                             const std::vector<double> &n) {
  std::vector<std::size_t> act;
  for (std::size_t r = 0; r < rx.size(); ++r) {
    bool ok = true;
    for (const auto &p : rx[r].nu) ok = ok && n[p.first] > 0.0;
    if (ok) act.push_back(r);
  }
  return act;
}

// Residuals F_a and Jacobian J_ab of the active reactions at n.
static void joint_residual(const std::vector<JointReaction> &rx,
                           const std::vector<std::size_t> &act, const std::vector<double> &n,
                           double P, const std::vector<double> &lnK, std::vector<double> &F,
                           std::vector<double> *J) {
  const std::size_t m = act.size();
  double nt = 0.0;
  for (double v : n) nt += v;
  const double ln_nt = std::log(nt);
  const double ln_P = std::log(P / EQUILIBRIUM_P_REF);
  F.assign(m, 0.0);
  for (std::size_t a = 0; a < m; ++a) {
    const JointReaction &r = rx[act[a]];
    double lnQ = r.dn * (ln_P - ln_nt);
    for (const auto &p : r.nu) lnQ += p.second * std::log(n[p.first]);
    F[a] = lnQ - lnK[act[a]];
  }
  if (J == nullptr) return;
  J->assign(m * m, 0.0);
  for (std::size_t a = 0; a < m; ++a) {
    const JointReaction &ra = rx[act[a]];
    for (std::size_t b = 0; b < m; ++b) {
      const JointReaction &rb = rx[act[b]];
      double s = -ra.dn * rb.dn / nt;
      for (const auto &pa : ra.nu) {
        for (const auto &pb : rb.nu) {
          if (pa.first == pb.first) s += pa.second * pb.second / n[pa.first];
        }
      }
      (*J)[a * m + b] = s;
    }
  }
}

static double max_abs(const std::vector<double> &v) {
  double s = 0.0;
  for (double x : v) s = std::max(s, std::abs(x));
  return s;
}

// Newton on the active extents from n (updated in place). Each step solves
// the diagonally scaled system (the hydrocarbon rows carry 1/n_hc, up to
// 1e30 above the others), is cut back to keep every amount above
// (1 - REFORMING_FRACTION_TO_BOUNDARY) of its current value, and backtracks
// until ||F|| falls. Converged when max |F| <= REFORMING_NEWTON_FTOL, after
// which up to two more steps take it to round-off.
// The scale s = 1/sqrt(J_aa) of a reaction whose key is so small (below
// about 1e-308, subnormal: steam at 4e-315 beside heptane with 1e-20 of air,
// #1287) that its nu_k^2 / n_k term of J_aa overflows to inf. A key is in its
// own reaction only, so the inf is on the diagonal alone, and the finite form
//   s = sqrt(n_k) / sqrt(J'),  J' = nu_k^2 + n_k R,
// R the rest of J_aa, gives the same scale; s / n_k = 1 / sqrt(n_k J') is the
// factor the key's own 1/n_k terms take. The scaled diagonal is 1.
static double subnormal_key_scale(const JointReaction &r, const std::vector<double> &n,
                                  double nt, double *s_over_nk = nullptr) {
  const auto &kp = r.nu[r.key];
  double R = -r.dn * r.dn / nt;
  for (const auto &p : r.nu) {
    if (p.first != kp.first) R += p.second * p.second / n[p.first];
  }
  const double nk = n[kp.first];
  const double Jp = kp.second * kp.second + nk * R;
  // sqrt(n_k) is a normal number (about 1e-159): n_k / J' and n_k J' are
  // not formed, as they would be further below 1e-308, with fewer digits.
  const double rk = std::sqrt(nk), rJ = std::sqrt(Jp);
  if (s_over_nk != nullptr) *s_over_nk = 1.0 / (rk * rJ);
  return rk / rJ;
}

static bool joint_newton(std::vector<double> &n, double P, const std::vector<JointReaction> &rx,
                         const std::vector<std::size_t> &act, const std::vector<double> &lnK,
                         int max_it) {
  const std::size_t m = act.size();
  if (m == 0) return true;
  std::vector<double> F, J, F_try, n_try;
  joint_residual(rx, act, n, P, lnK, F, &J);
  double fnorm = max_abs(F);
  int polish = 0;
  for (int it = 0; it < max_it; ++it) {
    if (!std::isfinite(fnorm)) return false;
    if (fnorm <= REFORMING_NEWTON_FTOL) {
      if (polish >= 2 || fnorm == 0.0) return true;
      ++polish;
    }
    std::vector<double> s(m);
    for (std::size_t a = 0; a < m; ++a) s[a] = 1.0 / std::sqrt(std::max(J[a * m + a], 1e-300));
    // A subnormal key's J_aa is inf (subnormal_key_scale).
    std::vector<bool> inf_diag(m, false);
    for (std::size_t a = 0; a < m; ++a) {
      if (std::isfinite(J[a * m + a])) continue;
      inf_diag[a] = true;
      double nt_a = 0.0;
      for (double v : n) nt_a += v;
      s[a] = subnormal_key_scale(rx[act[a]], n, nt_a);
    }
    std::vector<double> A(m * m), dx(m);
    for (std::size_t a = 0; a < m; ++a) {
      for (std::size_t b = 0; b < m; ++b) A[a * m + b] = s[a] * J[a * m + b] * s[b];
      if (inf_diag[a]) A[a * m + a] = 1.0;
      dx[a] = -s[a] * F[a];
    }
    if (!solve_dense(A, dx, m, 1)) return polish > 0;
    for (std::size_t a = 0; a < m; ++a) dx[a] *= s[a];
    // The step along a curve tangent to the Newton direction: a reaction's
    // key species (its hydrocarbon; CO2 for WGS) moves in ln n, as
    // n exp(dn/n), so it stays positive and a trace amount e^300 off its
    // equilibrium (heptane at 2000 K) closes in a step or two instead of a
    // factor 1/(1 - tau) or 1 + dn/n per step. Where the step raises it, the
    // log step may take it up to REFORMING_TRACE_FRACTION of the total and the
    // linear one beyond (exp would overshoot a species that grows from a trace
    // to a major one), as the trace-species damping of Gordon and McBride
    // (NASA RP-1311, 1994, Eq. 3.2). The extent follows from the key, the
    // other species from the extents.
    std::vector<double> dn(n.size(), 0.0);
    std::vector<bool> is_key(n.size(), false);
    for (std::size_t a = 0; a < m; ++a) {
      is_key[rx[act[a]].nu[rx[act[a]].key].first] = true;
      for (const auto &p : rx[act[a]].nu) dn[p.first] += p.second * dx[a];
    }
    // Fraction to the boundary on the other species (first guess; the line
    // search enforces it on the curved step).
    double alpha = 1.0;
    for (std::size_t j = 0; j < n.size(); ++j) {
      if (dn[j] < 0.0 && !is_key[j]) {
        alpha = std::min(alpha, REFORMING_FRACTION_TO_BOUNDARY * n[j] / -dn[j]);
      }
    }
    double nt = 0.0;
    for (double v : n) nt += v;
    auto curved = [&](double al, std::vector<double> &out) {
      out = n;
      for (std::size_t a = 0; a < m; ++a) {
        const JointReaction &r = rx[act[a]];
        const auto &kp = r.nu[r.key];
        const double nk = n[kp.first];
        const double d = kp.second * al * dx[a];
        if (nk >= REFORMING_TRACE_FRACTION * nt) {
          // A major key (methane in a fuel-only burn) steps linearly, and the
          // extent is taken as it is: from the key's own change it would be
          // quantised to ulp(n_key), and the reaction's trace species (steam
          // at 1e-8 beside methane at 1) then could not move by less than
          // 1e-16 / 1e-8 of themselves -- a residual floor of 1e-9, above
          // REFORMING_NEWTON_FTOL (#1262).
          const double xi = al * dx[a];
          for (const auto &p : r.nu) out[p.first] += p.second * xi;
          continue;
        }
        // The key's change, n_k expm1(d / n_k), and the extent from it
        // directly: from the difference n_new - n_k the extent was quantised
        // to ulp(n_k), and the other species with it. Methane at 2.8e-6 moved
        // trace H2 (7e-12, four per methane) by 2.3e-10 of itself per ulp, a
        // ln Q floor of 1.0e-9 in heptane's row (22 H2), above
        // REFORMING_NEWTON_FLOOR_FTOL: very lean burns below 245 K were
        // refused (#1285).
        double dk = nk * std::expm1(d / nk);
        if (d > 0.0) {
          const double ln_cap = std::log(REFORMING_TRACE_FRACTION * nt);
          dk = std::max(d, nk * std::expm1(std::min(d / nk, ln_cap - std::log(nk))));
        }
        const double xi = dk / kp.second;
        for (const auto &p : r.nu) out[p.first] += p.second * xi;
        // The key itself from exp where it falls steeply: n_k + dk cancels
        // there (its relative error grows as eps e^(-d / n_k)), and a trace
        // e^300 off would close by only e^-37 per step.
        out[kp.first] =
            d / nk < REFORMING_LOG_KEY_EXP_BELOW ? nk * std::exp(d / nk) : nk + dk;
      }
      for (std::size_t j = 0; j < n.size(); ++j) {
        if (!is_key[j] && !(out[j] >= (1.0 - REFORMING_FRACTION_TO_BOUNDARY) * n[j])) {
          return false;
        }
        if (n[j] > 0.0 && !(out[j] > 0.0)) return false;
      }
      return true;
    };
    bool accepted = false;
    for (int ls = 0; ls < 60; ++ls) {
      if (curved(alpha, n_try)) {
        joint_residual(rx, act, n_try, P, lnK, F_try, nullptr);
        const double f_try = max_abs(F_try);
        // Sufficient decrease of ||F||: the Newton direction descends it. A
        // step that does not lower ||F|| at all is no progress, also where
        // alpha is so small that the bound rounds to fnorm itself (it was
        // taken, 200 times over, at alpha 4.5e-13).
        if (std::isfinite(f_try) && f_try <= (1.0 - 1e-4 * alpha) * fnorm && f_try < fnorm) {
          accepted = true;
          break;
        }
      }
      alpha *= 0.5;
    }
    if (!accepted) {
      // Round-off floor: no step lowers ||F||. Converged if each row is
      // within REFORMING_NEWTON_FLOOR_FTOL, the floor a trace species leaves,
      // plus REFORMING_NEWTON_ULP_FLOOR_STEPS ulp steps of its amounts (a
      // subnormal one holds fewer digits; #1287).
      for (std::size_t a = 0; a < m; ++a) {
        double floor_a = REFORMING_NEWTON_FLOOR_FTOL;
        for (const auto &p : rx[act[a]].nu) {
          const double v = n[p.first];
          const double ulp = std::nextafter(v, std::numeric_limits<double>::infinity()) - v;
          floor_a += REFORMING_NEWTON_ULP_FLOOR_STEPS * std::abs(p.second) * ulp / v;
        }
        if (!(std::abs(F[a]) <= floor_a)) return false;
      }
      return true;
    }
    n = n_try;
    joint_residual(rx, act, n, P, lnK, F, &J);
    fnorm = max_abs(F);
  }
  return fnorm <= REFORMING_NEWTON_FTOL;
}

// Element-potential start (Gordon and McBride, NASA RP-1311, 1994, Ch. 2-3)
// for when Newton on the extents does not converge from the rounds: far from
// equilibrium a species that must grow by orders of magnitude (the steam of a
// steam-starved rich burn, 1e-6 -> 5e-2) makes the linearisation in the
// extents useless and the fraction-to-boundary steps tiny. The same
// equilibrium (the species of the active reactions, the others held fixed as
// inerts) is solved for the element potentials pi and ln N, with the
// corrections in ln n_j:
//   d ln n_j = -mu_j + sum_e a_ej pi_e + d ln N,   mu_j = g_j/RT + ln(n_j P/(N P0)),
// damped as RP-1311 Eqs. 3.1-3.2. The result is then written back as extents
// from n (so the elements are conserved exactly) for the Newton polish.
// Returns false if it does not converge or an amount comes back non-positive.
// rebase (optional) is called with the element potentials' amounts before
// they are written back as extents, and may replace rx and act (the caller's
// own, which these reference) by another basis of the same reactions: the
// write-back takes each non-key amount as n plus the extents' changes, so a
// key chosen at the start can leave a carrier the difference of larger
// changes (steam at 9e-15 came back as -7e-15 from 4e-3; #1262).
static bool element_potential_start(
    std::vector<double> &n, double T, double P, const std::vector<JointReaction> &rx,
    const std::vector<std::size_t> &act,
    const std::function<void(const std::vector<double> &)> &rebase = nullptr) {
  std::vector<std::size_t> sp;  // species of the active reactions
  for (std::size_t a : act) {
    for (const auto &p : rx[a].nu) {
      if (std::find(sp.begin(), sp.end(), p.first) == sp.end()) sp.push_back(p.first);
    }
  }
  const std::size_t S = sp.size();
  std::vector<std::vector<double>> aj(S);  // atoms C, H, O
  for (std::size_t q = 0; q < S; ++q) {
    const auto &ms = molecular_structures[sp[q]];
    aj[q] = {static_cast<double>(ms.C), static_cast<double>(ms.H), static_cast<double>(ms.O)};
  }
  double n_inert = 0.0;
  for (std::size_t j = 0; j < n.size(); ++j) {
    if (std::find(sp.begin(), sp.end(), j) == sp.end()) n_inert += n[j];
  }
  std::vector<double> b(3, 0.0);
  for (std::size_t q = 0; q < S; ++q) {
    for (std::size_t e = 0; e < 3; ++e) b[e] += aj[q][e] * n[sp[q]];
  }
  std::vector<std::size_t> el;
  for (std::size_t e = 0; e < 3; ++e) {
    if (b[e] > 0.0) el.push_back(e);
  }
  const std::size_t E = el.size();
  // The element equations in a basis of component species (#1287): with
  // heptane the only major species of a fuel with 1e-20 of air, every element
  // row is dominated by it (C, H and the total moles in a ratio 7 : 16 : 1),
  // and the traces' share of each row, the only thing that sets their
  // distribution, is below the rounding of its sum; their corrections did
  // not converge, the steam of the oxygen row came back 1e10 times its
  // element's total, and hexane, the trace that holds the excess hydrogen,
  // was missed. As Smith and Missen's component (stoichiometric)
  // formulation, the E most abundant species that are independent are the
  // components, C their formula matrix, and each formula is replaced by
  // a''_q = adj(C) a_q -- a linear map, so the same element constraints and
  // the same equilibrium. Formulas are small integers, so each component's
  // a'' is det(C) e_k exactly: a major species contributes exactly nothing to
  // a trace component's row, whose sums (and right-hand side, summed per
  // species) keep their own relative precision. The components are re-chosen
  // as the amounts move.
  std::vector<double> n_init(S);
  for (std::size_t q = 0; q < S; ++q) n_init[q] = n[sp[q]];
  std::vector<std::vector<double>> ap(S, std::vector<double>(E, 0.0));
  std::vector<double> bp(E, 0.0), bscale(E, 0.0);
  std::vector<std::size_t> comp;
  auto choose_components = [&](const std::vector<double> &lnv) -> bool {
    std::vector<std::size_t> order(S);
    for (std::size_t q = 0; q < S; ++q) order[q] = q;
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t x, std::size_t y) { return lnv[x] > lnv[y]; });
    std::vector<std::size_t> c;
    std::vector<std::vector<double>> red;  // the components' reduced formulas
    for (std::size_t q : order) {
      if (c.size() == E) break;
      std::vector<double> v(E);
      for (std::size_t k = 0; k < E; ++k) v[k] = aj[q][el[k]];
      for (const auto &w : red) {
        std::size_t pk = 0;
        for (std::size_t k = 1; k < E; ++k) {
          if (std::abs(w[k]) > std::abs(w[pk])) pk = k;
        }
        const double f = v[pk] / w[pk];
        for (std::size_t k = 0; k < E; ++k) v[k] -= f * w[k];
      }
      double vmax = 0.0;
      for (double x : v) vmax = std::max(vmax, std::abs(x));
      if (vmax > REFORMING_EP_COMPONENT_RANK_TOL) {
        red.push_back(v);
        c.push_back(q);
      }
    }
    if (c.size() < E) return false;
    if (c == comp) return true;
    comp = c;
    // C[i][k]: element el[i] of component k, embedded in a 3 x 3 identity for
    // E < 3 (the E x E block of adj(diag(C, I)) is adj(C)), so one cofactor
    // form serves every E: adj(C)[i][k] = cofactor C_ki.
    double C3[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    for (std::size_t i = 0; i < E; ++i) {
      for (std::size_t k = 0; k < E; ++k) C3[i * 3 + k] = aj[comp[k]][el[i]];
    }
    std::vector<double> adj(E * E, 0.0);
    for (std::size_t i = 0; i < E; ++i) {
      for (std::size_t k = 0; k < E; ++k) {
        const std::size_t r0 = (k + 1) % 3, r1 = (k + 2) % 3;
        const std::size_t c0 = (i + 1) % 3, c1 = (i + 2) % 3;
        adj[i * E + k] = C3[r0 * 3 + c0] * C3[r1 * 3 + c1] - C3[r0 * 3 + c1] * C3[r1 * 3 + c0];
      }
    }
    for (std::size_t q = 0; q < S; ++q) {
      for (std::size_t i = 0; i < E; ++i) {
        double v = 0.0;
        for (std::size_t k = 0; k < E; ++k) v += adj[i * E + k] * aj[q][el[k]];
        ap[q][i] = v;
      }
    }
    for (std::size_t i = 0; i < E; ++i) {
      bp[i] = 0.0;
      bscale[i] = 0.0;
      for (std::size_t q = 0; q < S; ++q) {
        bp[i] += ap[q][i] * n_init[q];
        bscale[i] += std::abs(ap[q][i]) * n_init[q];
      }
    }
    return true;
  };
  std::vector<double> lnn(S), g(S);
  for (std::size_t q = 0; q < S; ++q) {
    lnn[q] = std::log(n[sp[q]]);
    g[q] = g_over_RT(sp[q], T);
  }
  double lnN = 0.0;
  {
    double nt = n_inert;
    for (std::size_t q = 0; q < S; ++q) nt += n[sp[q]];
    lnN = std::log(nt);
  }
  const double ln_P = std::log(P / EQUILIBRIUM_P_REF);
  const double SIZE = 18.420681;  // RP-1311: ln(1e8), a trace species
  bool done = false;
  double last_weighted = std::numeric_limits<double>::infinity();
  for (int it = 0; it < REFORMING_NEWTON_MAX_IT && !done; ++it) {
    if (!choose_components(lnn)) return false;
    const double N = std::exp(lnN);
    std::vector<double> nq(S), mu(S);
    double sum_n = 0.0;
    for (std::size_t q = 0; q < S; ++q) {
      nq[q] = std::exp(lnn[q]);
      mu[q] = g[q] + lnn[q] - lnN + ln_P;
      sum_n += nq[q];
    }
    const std::size_t M = E + 1;
    std::vector<double> A(M * M, 0.0), r(M, 0.0);
    for (std::size_t k = 0; k < E; ++k) {
      double ak_n = 0.0;
      for (std::size_t q = 0; q < S; ++q) {
        const double akq = ap[q][k] * nq[q];
        for (std::size_t i = 0; i < E; ++i) A[k * M + i] += akq * ap[q][i];
        ak_n += akq;
        r[k] += akq * mu[q];
      }
      A[k * M + E] = ak_n;
      A[E * M + k] = ak_n;
      r[k] += bp[k] - ak_n;
    }
    A[E * M + E] = sum_n - N;
    r[E] = N - n_inert - sum_n;
    for (std::size_t q = 0; q < S; ++q) r[E] += nq[q] * mu[q];
    // Symmetric diagonal scaling: a trace component's row is some 1e-20 of a
    // major one's.
    std::vector<double> D(M, 1.0);
    for (std::size_t k = 0; k < E; ++k) {
      if (A[k * M + k] > 0.0) D[k] = 1.0 / std::sqrt(A[k * M + k]);
    }
    D[E] = 1.0 / std::sqrt(N);
    for (std::size_t i = 0; i < M; ++i) {
      for (std::size_t k = 0; k < M; ++k) A[i * M + k] *= D[i] * D[k];
      r[i] *= D[i];
    }
    if (!solve_dense(A, r, M, 1)) return false;
    for (std::size_t i = 0; i < M; ++i) r[i] *= D[i];
    const double dlnN = r[E];
    std::vector<double> dl(S);
    double lam1 = 5.0 * std::abs(dlnN), dmax = 0.0;
    double lam2 = 1.0;
    for (std::size_t q = 0; q < S; ++q) {
      dl[q] = -mu[q] + dlnN;
      for (std::size_t i = 0; i < E; ++i) dl[q] += ap[q][i] * r[i];
      dmax = std::max(dmax, std::abs(dl[q]));
      if (lnn[q] - lnN > -SIZE) {
        lam1 = std::max(lam1, std::abs(dl[q]));
      } else if (dl[q] >= 0.0 && dl[q] - dlnN != 0.0) {
        lam2 = std::min(lam2, std::abs((-(lnn[q] - lnN) - 9.2103404) / (dl[q] - dlnN)));
      }
    }
    const double lam = std::min({1.0, lam1 > 0.0 ? 2.0 / lam1 : 1.0, lam2});
    for (std::size_t q = 0; q < S; ++q) lnn[q] += lam * dl[q];
    lnN += lam * dlnN;
    if (!std::isfinite(lnN)) return false;
    done = lam == 1.0 && std::max(dmax, std::abs(dlnN)) < 1e-9;
    if (lam == 1.0) {
      double w = std::abs(dlnN);
      for (std::size_t q = 0; q < S; ++q) w = std::max(w, nq[q] * std::abs(dl[q]) / N);
      last_weighted = w;
    }
  }
  // Only a start for the Newton polish, which decides convergence: with one
  // element nearly absent (oxygen at 1e-11 in methane) the corrections stall
  // at a 1e-8 noise in ln n above the 1e-9 bound, from a start already that
  // close (#1262).
  // A trace far below the others can keep a correction noise above the bound
  // (1e-6 in ln n, about 1e-6 of a species 1e-20 of the total: lean methane at
  // 192 K, 77 bar; #1285); RP-1311's own test weights each by its amount,
  // n_j |d ln n_j| / N (Eq. 3.5), and below REFORMING_EP_START_CORRECTION
  // that hands over too.
  // The weighted correction is never above the full one (n_j / N <= 1), so
  // it alone decides the handover.
  if (!done && !(last_weighted < REFORMING_EP_START_CORRECTION)) {
    return false;
  }
  // The corrections stop near 1e-10 of each amount, so the amounts balance
  // the elements only to that; the write-back below takes each non-key as n
  // plus the extents' changes, which then carries the whole imbalance, and a
  // trace came back negative (H2 at 4.7e-17 as -5e-17 beside methane at 3e-7
  // in a lean ethane burn at 155 K; #1285). The amounts are first balanced
  // to rounding: n_j exp(sum_e a''_ej lambda_e), lambda by Newton (the
  // element equations alone, in the components' basis, each row to its own
  // sum of |a''| n), which keeps every amount positive. A sum of S terms can
  // leave about S eps, so the loop may end on its step cap instead of
  // REFORMING_EP_BALANCE_RTOL; the element check after Newton decides.
  for (int it = 0; it < REFORMING_EP_BALANCE_IT; ++it) {
    std::vector<double> Jl(E * E, 0.0), rl(E, 0.0);
    bool balanced = true;
    for (std::size_t k = 0; k < E; ++k) {
      double bk = 0.0;
      for (std::size_t q = 0; q < S; ++q) {
        const double akq = ap[q][k] * std::exp(lnn[q]);
        bk += akq;
        for (std::size_t i = 0; i < E; ++i) Jl[k * E + i] += akq * ap[q][i];
      }
      rl[k] = bp[k] - bk;
      if (!(std::abs(rl[k]) <= REFORMING_EP_BALANCE_RTOL * bscale[k])) balanced = false;
    }
    if (balanced) break;
    if (!solve_dense(Jl, rl, E, 1)) return false;
    for (std::size_t q = 0; q < S; ++q) {
      for (std::size_t k = 0; k < E; ++k) lnn[q] += ap[q][k] * rl[k];
    }
  }
  if (rebase) {
    std::vector<double> n_gm = n;
    for (std::size_t q = 0; q < S; ++q) n_gm[sp[q]] = std::exp(lnn[q]);
    rebase(n_gm);
  }
  // Back to extents from the current n: n_key = n_gm,key fixes each extent.
  std::vector<double> n_new = n;
  for (std::size_t a : act) {
    const auto &kp = rx[a].nu[rx[a].key];
    const std::size_t qk = static_cast<std::size_t>(
        std::find(sp.begin(), sp.end(), kp.first) - sp.begin());
    const double xi = (std::exp(lnn[qk]) - n[kp.first]) / kp.second;
    for (const auto &p : rx[a].nu) n_new[p.first] += p.second * xi;
  }
  for (std::size_t a : act) {
    const std::size_t j = rx[a].nu[rx[a].key].first;
    n_new[j] = std::exp(lnn[static_cast<std::size_t>(std::find(sp.begin(), sp.end(), j) -
                                                     sp.begin())]);
  }
  for (std::size_t j : sp) {
    if (!(n_new[j] > 0.0)) return false;
  }
  n = n_new;
  return true;
}

// Absent species of the reactions (#1262). A reaction with an absent species
// is left out of the Newton solve, and a round cannot make one where CO and
// H2 are below 1e-12 (the scalar reforming bracket is empty), so a
// hydrocarbon the equilibrium holds at percent level could stay exactly 0 at
// a state reported converged: ethane at phi 2400, 266 K, 4.2 bar had no
// methane (0.71 at equilibrium) and came out 276 K against 441 K; heptane in
// air stepped 2.55 K between phi 87.45 and 87.475 as isobutane dropped out.
// A reaction with an absent species whose other side is all present has an
// unbounded driving force toward the absent one (ln Q = -inf there): it is
// run that way by REFORMING_SEED_FRACTION of the total (at most half of any
// amount it uses), which seeds the absent species in its own stoichiometry,
// conserving every element exactly. (Seeding each absent species by a
// least-norm combination of the present ones lost conservation where an
// element was a trace: oxygen at 1e-12 left the combination's system with a
// condition number of 1e12, and H/C moved by 2.8e-5.) Returns whether any
// was; runnable (if given) receives whether one could still be.
static bool seed_absent_species(std::vector<double> &n, const std::vector<JointReaction> &rx,
                                bool *runnable = nullptr) {
  double nt = 0.0;
  for (double v : n) nt += v;
  bool seeded = false, can = false;
  for (const auto &r : rx) {
    bool any_absent = false, lhs = true, rhs = true;
    for (const auto &p : r.nu) {
      if (!(n[p.first] > 0.0)) {
        any_absent = true;
        if (p.second < 0.0) lhs = false;
        if (p.second > 0.0) rhs = false;
      }
    }
    if (!any_absent || (!lhs && !rhs)) continue;
    can = true;
    if (runnable != nullptr) continue;
    const double dir = lhs ? 1.0 : -1.0;  // toward the side with the absent species
    double xi = REFORMING_SEED_FRACTION * nt;
    for (const auto &p : r.nu) {
      if (dir * p.second < 0.0) xi = std::min(xi, 0.5 * n[p.first] / std::abs(p.second));
    }
    if (!(xi > 0.0)) continue;
    for (const auto &p : r.nu) n[p.first] = std::max(0.0, n[p.first]) + dir * p.second * xi;
    seeded = true;
  }
  if (runnable != nullptr) *runnable = can;
  return seeded;
}

// Joint reforming + WGS equilibrium at fixed T from n0. Rounds start it
// (REFORMING_START_PASSES, which also make the species of a reaction that can
// move present); Newton then solves the reactions whose species are all
// present. If some reaction was left out, a trial round checks that it cannot
// move; if it can, the round is taken and Newton re-solved. Returns whether
// the residuals met REFORMING_NEWTON_FTOL and no left-out reaction moves.
// active_out (optional) receives the reactions solved by Newton.
// Elements: every step and seed conserves them; a state that does not, to
// round-off of each element's own total, is not a solution (#1262 review,
// #1285 review).
static bool elements_conserved(const std::vector<double> &n0, const std::vector<double> &n) {
  double b0[3] = {0.0, 0.0, 0.0}, b1[3] = {0.0, 0.0, 0.0};
  for (std::size_t j = 0; j < n0.size(); ++j) {
    const auto &ms = molecular_structures[j];
    const double a[3] = {static_cast<double>(ms.C), static_cast<double>(ms.H),
                         static_cast<double>(ms.O)};
    for (std::size_t e = 0; e < 3; ++e) {
      b0[e] += a[e] * n0[j];
      b1[e] += a[e] * n[j];
    }
  }
  for (std::size_t e = 0; e < 3; ++e) {
    if (!(std::abs(b1[e] - b0[e]) <= REFORMING_ELEMENT_RTOL * b0[e])) return false;
  }
  return true;
}

static bool solve_reforming_wgs_joint(const std::vector<double> &n0, double T, double P,
                                      const ReformingConfig &cfg, bool with_reforming,
                                      std::vector<double> &n_out,
                                      std::vector<JointReaction> *rx_out = nullptr,
                                      std::vector<std::size_t> *active_out = nullptr) {
  std::vector<JointReaction> rx = joint_reactions(cfg, with_reforming);
  std::vector<double> lnK;

  n_out = n0;
  std::vector<double> before;
  auto round = [&](std::vector<double> &n) {
    if (with_reforming) {
      reforming_wgs_round(n, T, P, cfg);
    } else {
      const WgsConfig w = make_wgs_config();
      const double lo = std::max(-n[w.i_CO2], -n[w.i_H2]);
      const double hi = std::min(n[w.i_CO], n[w.i_H2O]);
      if (hi > lo + 1e-15) {
        const double xi = solve_wgs_extent(n, T, w);
        n[w.i_CO] -= xi;
        n[w.i_H2O] -= xi;
        n[w.i_CO2] += xi;
        n[w.i_H2] += xi;
      }
    }
  };
  for (int pass = 0; pass < REFORMING_START_PASSES; ++pass) {
    before = n_out;
    round(n_out);
    if (reforming_round_change(before, n_out) <= REFORMING_PASS_RTOL) break;
  }
  // A hydrocarbon not a trace (it needs no log step of its own): the basis
  // in which the interconversions carry no oxygen species, its keys the two
  // scarcer oxygen carriers at n (#1262). Re-chosen whenever the start moves:
  // a carrier that was the scarcest after the rounds can be the most abundant
  // after the element potentials, and a key's own ulp then quantises the
  // others (CO at 5e-18 stepped 1.9e-26 at a time beside steam at 2e-11).
  auto choose_basis = [&](const std::vector<double> &n) {
    rx = joint_reactions(cfg, with_reforming, n[cfg.i_CO] < n[cfg.i_CO2]);
    std::vector<JointReaction> rich;
    if (with_reforming && joint_reactions_rich(cfg, n, rich)) rx = rich;
    if (with_reforming) {
      // Then the component basis where its scarcest non-key is more
      // abundant than this one's (#1288): a non-key far below the extents'
      // changes cancels, and the element potentials' state came back with it
      // negative. H2 fuel with a trace of CH4 at 300 K methanates to CO2 at
      // 2e-30 beside CH4 at 3e-5, CO2 a non-key of the steam-reforming basis
      // (written back as -1e-21); CO fuel with 0.4 % methane at 250 K leaves
      // H2 at 1e-22, a non-key of the rich basis (-4e-19).
      std::vector<JointReaction> comp;
      if (joint_reactions_components(cfg, n, comp) &&
          joint_min_nonkey(comp, n) > joint_min_nonkey(rx, n)) {
        rx = comp;
      }
    }
    lnK.assign(rx.size(), 0.0);
    for (std::size_t r = 0; r < rx.size(); ++r) lnK[r] = joint_lnK(rx[r], T);
  };
  choose_basis(n_out);

  bool conv = false;
  std::vector<std::size_t> act;
  for (int seed_pass = 0; seed_pass <= REFORMING_SEED_PASSES; ++seed_pass) {
  for (int restart = 0; restart <= REFORMING_JOINT_RESTARTS; ++restart) {
    if (restart > 0) choose_basis(n_out);
    act = joint_active(rx, n_out);
    std::vector<double> n_try = n_out;
    // A Newton solve that broke an element's balance is not converged either
    // (#1312): a basis chosen at the start can carry an element's only
    // carriers as components of reactions that move far more. Heavy alkane
    // with trace air and CO at 250-260 K: CO and CO2 were the C and O
    // components while the alkanes disproportionated (heptane 1e-10 -> 0.57),
    // so CO at 5e-7 came out of 0.1-mol cancellations and oxygen's own total
    // moved by 3e-12. The element-potential start lands near equilibrium,
    // where the basis is re-chosen with alkanes as components and oxygen in
    // reactions of its own.
    conv = joint_newton(n_try, P, rx, act, lnK, REFORMING_NEWTON_ROUNDS_IT) &&
           elements_conserved(n0, n_try);
    if (!conv) {
      // Far from equilibrium: restart Newton from the element potentials.
      std::vector<double> n_ep = n_out;
      const bool ep_ok = element_potential_start(
          n_ep, T, P, rx, act, [&](const std::vector<double> &n_gm) {
            choose_basis(n_gm);
            act = joint_active(rx, n_out);
          });
      if (ep_ok) {
        choose_basis(n_ep);
        const std::vector<std::size_t> act_ep = joint_active(rx, n_ep);
        if (joint_newton(n_ep, P, rx, act_ep, lnK, REFORMING_NEWTON_MAX_IT)) {
          n_try = n_ep;
          act = act_ep;
          conv = true;
        }
      }
      if (!conv) {
        // Back to the basis of the state kept.
        choose_basis(n_try);
        act = joint_active(rx, n_try);
      }
    }
    n_out = n_try;
    if (act.size() == rx.size()) break;
    // A left-out reaction: does a round move it?
    std::vector<double> n_round = n_out;
    round(n_round);
    if (reforming_round_change(n_out, n_round) <= REFORMING_JOINT_CHECK_RTOL) break;
    n_out = n_round;
    conv = false;
  }
    // An absent species with a positive affinity: seed it and solve again;
    // one still predicted after the last pass means not converged.
    if (!conv) break;
    if (seed_pass == REFORMING_SEED_PASSES) {
      bool runnable = false;
      seed_absent_species(n_out, rx, &runnable);
      if (runnable) conv = false;  // still driven toward an absent species
      break;
    }
    if (!seed_absent_species(n_out, rx)) break;
    conv = false;
    choose_basis(n_out);
  }
  if (!elements_conserved(n0, n_out)) conv = false;
  if (rx_out != nullptr) *rx_out = rx;
  if (active_out != nullptr) *active_out = act;
  return conv;
}

// Solve general reforming + WGS equilibrium at fixed T. Returns whether it
// converged.
static bool solve_reforming_wgs_isothermal(const std::vector<double> &n0,
                                           double T, double P,
                                           const ReformingConfig &cfg,
                                           std::vector<double> &n_out) {
  return solve_reforming_wgs_joint(n0, T, P, cfg, true, n_out);
}

// Adiabatic solver functor for general reforming + WGS
struct ReformingAdiabaticFun {
  const std::vector<double> &n_in;
  double H_mass_target;
  double P;
  ReformingConfig cfg;

  mutable std::vector<double> n;
  mutable std::vector<double> X;
  mutable double T_last = 0.0;

  ReformingAdiabaticFun(const std::vector<double> &n_in_, double H_mass_target_,
                        double P_, const ReformingConfig &cfg_)
      : n_in(n_in_), H_mass_target(H_mass_target_), P(P_), cfg(cfg_),
        n(n_in_.size()), X(n_in_.size()) {}

  // conv: whether the inner equilibrium at T converged (its flag was dropped
  // before #927).
  double f(double T, bool &conv) const {
    conv = solve_reforming_wgs_isothermal(n_in, T, P, cfg, n);
    molefractions(n, X);
    T_last = T;
    return h(T, X) / mwmix(X) - H_mass_target;
  }

  double slope_frozen() const { return cp(T_last, X) / mwmix(X); }
};

// Public API: General reforming + WGS equilibrium (isothermal)
EquilibriumResult reforming_equilibrium(const State &in) {
  ReformingConfig cfg = make_reforming_config();

  // Check if any hydrocarbon is present
  bool has_hc = false;
  for (std::size_t i_hc : cfg.hc_indices) {
    if (in.X[i_hc] > 1e-10) {
      has_hc = true;
      break;
    }
  }

  if (!has_hc) {
    // No hydrocarbons - fall back to WGS only
    return wgs_equilibrium(in);
  }

  std::vector<double> n_out;
  const bool conv = solve_reforming_wgs_isothermal(in.X, in.T, in.P, cfg, n_out);

  std::vector<double> X_out;
  molefractions(n_out, X_out);

  State out;
  out.set_TPX(in.T, in.P, X_out);
  return EquilibriumResult{out, in.T, 0.0, conv};
}

// Public API: General reforming + WGS equilibrium (adiabatic)
EquilibriumResult reforming_equilibrium_adiabatic(const State &in) {
  ReformingConfig cfg = make_reforming_config();

  // Check if any hydrocarbon is present
  bool has_hc = false;
  for (std::size_t i_hc : cfg.hc_indices) {
    if (in.X[i_hc] > 1e-10) {
      has_hc = true;
      break;
    }
  }

  if (!has_hc) {
    // No hydrocarbons - fall back to WGS only
    return wgs_equilibrium_adiabatic(in);
  }

  // Target enthalpy
  double H_mass_target = h(in.T, in.X) / mwmix(in.X);

  // Solve for adiabatic temperature
  ReformingAdiabaticFun F{in.X, H_mass_target, in.P, cfg};
  bool conv_T = true;
  double T_ad = solve_adiabatic_T(F, in.T, conv_T);

  // conv_T is false when an inner solve at either end of the final bracket
  // failed (#1262 review: five heavy rich states refused for their bracket
  // ends): Brent then also bracketed on that end's unconverged energy
  // residual, and T_ad is off its root by up to 1e-3 K. There, the root is
  // polished by secant steps on converged inner solves near T_ad, and the
  // state is converged if the equilibrium converges at it and the energy
  // balance holds to the bracket's own tolerance.
  if (!conv_T) {
    bool c0 = false;
    double f0 = F.f(T_ad, c0);
    double slope = F.slope_frozen();
    double T_prev = T_ad, f_prev = f0;
    for (int it = 0; it < EQUILIBRIUM_T_POLISH_IT && c0; ++it) {
      const double tol =
          4.0 * F.slope_frozen() * EQUILIBRIUM_T_XTOL_REL * std::max(T_ad, EQUILIBRIUM_T_MIN);
      if (std::abs(f0) <= tol) {
        conv_T = true;
        break;
      }
      const double T_new = T_ad - f0 / slope;
      bool c1 = false;
      const double f1 = F.f(T_new, c1);
      if (!c1 || !std::isfinite(f1)) break;
      T_prev = T_ad;
      f_prev = f0;
      T_ad = T_new;
      f0 = f1;
      c0 = c1;
      if (T_ad != T_prev) slope = (f0 - f_prev) / (T_ad - T_prev);
      if (!(slope > 0.0) || !std::isfinite(slope)) slope = F.slope_frozen();
    }
  }

  // Get equilibrium composition at T_ad
  std::vector<double> n_out;
  const bool conv = solve_reforming_wgs_isothermal(in.X, T_ad, in.P, cfg, n_out) && conv_T;

  std::vector<double> X_out;
  molefractions(n_out, X_out);

  State out;
  out.set_TPX(T_ad, in.P, X_out);
  return EquilibriumResult{out, in.T, T_ad - in.T, conv};
}

// -------------------------------------------------------------
// Combustion + Equilibrium (convenience function)
// -------------------------------------------------------------

EquilibriumResult combustion_equilibrium(const State &in, bool smooth_phi0, double k0) {
  // Step 1: Complete combustion (adiabatic), with the Phi = 1 softmin always on
  // (#1048). The reforming + WGS stage below forms no O2-side dissociation, so
  // without the softmin T keeps complete combustion's kink at Phi = 1 (one-sided
  // slopes +1450 and -960 K per unit Phi for CH4-air). The softmin is about
  // 1/SMOOTHING_K_PHI1 wide, as on the complete path. On its lean side it leaves
  // a trace of fuel, which reforming_equilibrium_adiabatic stops reforming once
  // below 1e-10 (|Phi - 1| ~ 5e-4): a step of about 1e-6 K in T there.
  State burned = complete_combustion(in, smooth_phi0, true, k0, SMOOTHING_K_PHI1);

  // Step 2: Reforming + WGS equilibrium on combustion products
  // For rich mixtures: unburned HC + H2O -> CO + H2
  // For all mixtures: CO + H2O <-> CO2 + H2
  EquilibriumResult result = reforming_equilibrium_adiabatic(burned);
  return EquilibriumResult{result.state, in.T, result.state.T - in.T,
                           result.converged};
}

// -------------------------------------------------------------
// Combustion + equilibrium with analytic derivatives (#716)
// -------------------------------------------------------------

namespace {

// A forward-mode dual with a runtime number of partials: the solver's inputs
// (here the burned moles n0, T, P).
struct DV {
  double v = 0.0;
  std::vector<double> d;
};
DV dv_const(double v, std::size_t N) { return DV{v, std::vector<double>(N, 0.0)}; }
DV operator+(const DV &a, const DV &b) {
  DV r{a.v + b.v, a.d};
  for (std::size_t i = 0; i < r.d.size(); ++i) r.d[i] += b.d[i];
  return r;
}
DV operator*(double c, const DV &a) {
  DV r{c * a.v, a.d};
  for (double &x : r.d) x *= c;
  return r;
}
DV operator*(const DV &a, const DV &b) {
  DV r{a.v * b.v, std::vector<double>(a.d.size())};
  for (std::size_t i = 0; i < r.d.size(); ++i) r.d[i] = a.d[i] * b.v + a.v * b.d[i];
  return r;
}
DV operator/(const DV &a, const DV &b) {
  DV r{a.v / b.v, std::vector<double>(a.d.size())};
  for (std::size_t i = 0; i < r.d.size(); ++i) r.d[i] = (a.d[i] * b.v - a.v * b.d[i]) / (b.v * b.v);
  return r;
}
// f(a) with slope s, lifted.
DV dv_lift(double f, double s, const DV &a) {
  DV r{f, a.d};
  for (double &x : r.d) x *= s;
  return r;
}

}  // namespace

EquilibriumJacobian combustion_equilibrium_and_jacobian(const State &in, bool smooth_phi0,
                                                        double k0) {
  const std::size_t ns = in.X.size();
  const std::size_t iT = ns, iP = ns + 1, N = ns + 2;
  EquilibriumJacobian out;
  const EquilibriumResult res = combustion_equilibrium(in, smooth_phi0, k0);
  out.state = res.state;
  const double T_ad = res.state.T;
  const double P = in.P;

  // Stage 1: complete combustion (phi = 1 softmin on, as the value path).
  const CompleteCombustionJacobian cj =
      complete_combustion_products_and_jacobian(in.X, smooth_phi0, true, k0, SMOOTHING_K_PHI1);
  const std::vector<double> &n0 = cj.X_b;

  // Stage 2: the joint reforming + WGS equilibrium at T_ad (#759), as the
  // value path: reform if the BURNED state has any hydrocarbon, else WGS
  // alone (as reforming_equilibrium_adiabatic).
  const ReformingConfig cfg = make_reforming_config();
  bool has_hc = false;
  for (std::size_t i_hc : cfg.hc_indices) {
    if (n0[i_hc] > 1e-10) has_hc = true;
  }
  std::vector<double> nv;
  std::vector<JointReaction> rx;
  std::vector<std::size_t> act;
  const bool joint_conv = solve_reforming_wgs_joint(n0, T_ad, P, cfg, has_hc, nv, &rx, &act);
  // The derivatives below are the implicit function's at a root: at a state
  // the solve did not converge to they are not those of any F (#1262: ln Q -
  // ln K up to 62 at a cold, very rich burn gave dT/dP 30 % off), so neither
  // the state nor its Jacobian is returned. The network takes the exception
  // as an unphysical probe.
  if (!res.converged || !joint_conv) {
    throw std::runtime_error(
        "combustion_equilibrium_and_jacobian: the reforming + WGS equilibrium did not "
        "converge (T_in = " + std::to_string(in.T) + " K, P = " + std::to_string(P) + " Pa)");
  }
  // A species absent at the solution that the complete-combustion step feeds
  // (the oxygen carriers of a fuel-only mix as air starts to arrive): its
  // reactions are left out, so the derivatives below would keep it at 0,
  // while the solution an instant later holds it in the ratios equilibrium
  // sets among the trace species, linear in its amount (#1251: J 2.260e6
  // against a forward difference of 2.357e6 K per kg/s at exactly zero air).
  // The derivatives are then those at the same state with each such species
  // at a trace of the total (REFORMING_JACOBIAN_TRACE, up by 100 until the
  // solve converges) and equilibrated: the one-sided limit to that order. The
  // value stays the solution's.
  if (act.size() < rx.size()) {
    double nt0 = 0.0;
    for (double v : n0) nt0 += v;
    // Only species of the reactions left out: an absent species no reaction
    // of the solve holds (a hydrocarbon of a lean burn's WGS-only solve) has
    // no extent to differentiate.
    std::vector<bool> in_left_out(ns, false);
    for (std::size_t r = 0; r < rx.size(); ++r) {
      if (std::find(act.begin(), act.end(), r) != act.end()) continue;
      for (const auto &p : rx[r].nu) in_left_out[p.first] = true;
    }
    std::vector<bool> n0v_fed(ns, false);
    bool fed = false;
    for (std::size_t j = 0; j < ns; ++j) {
      if (nv[j] > 0.0 || !in_left_out[j]) continue;
      bool any = false;
      for (std::size_t k = 0; k < ns; ++k) any = any || cj.dXb_dX[j][k] != 0.0;
      if (any) {
        n0v_fed[j] = true;
        fed = true;
      }
    }
    if (fed) {
      // The smallest trace in the ladder at which the solve converges with
      // the fed species' reactions active: 1e-30 converged for methane only
      // (ethane, propane and blends from ~1e-22 up), and the no-trace
      // derivatives kept as a fallback were 34-45 % off for them (#1251
      // review 4). None converging is refused, as an unconverged solve.
      bool got = false;
      for (double trace = REFORMING_JACOBIAN_TRACE; trace <= REFORMING_JACOBIAN_TRACE_MAX;
           trace *= 100.0) {
        std::vector<double> n0v = n0;
        for (std::size_t j = 0; j < ns; ++j) {
          if (n0v_fed[j]) n0v[j] += trace * nt0;
        }
        std::vector<double> nv2;
        std::vector<JointReaction> rx2;
        std::vector<std::size_t> act2;
        if (solve_reforming_wgs_joint(n0v, T_ad, P, cfg, has_hc, nv2, &rx2, &act2) &&
            act2.size() > act.size()) {
          nv = nv2;
          rx = rx2;
          act = act2;
          got = true;
          break;
        }
      }
      if (!got) {
        // The shift alone (no hydrocarbon), left out, and unable to run either
        // way with the fed species seeded at the solution (an absent species
        // on each side): its extent is boxed at zero an instant later (CO =
        // -xi >= 0 and H2 = xi >= 0), so the fed species stays as fed and the
        // derivatives above are exact. A burn that lacks steam or CO2, lean,
        // with only that carrier fed (#1293): pure CO fuel in dry air holds
        // CO2 and none of CO, H2 or steam, and only steam is fed; H2 or NH3
        // in O2/N2 likewise with CO2 fed. Every trace solve converged with the
        // shift still left out, and was refused. With reforming, reactions can
        // combine where none runs alone, so that stays refused.
        bool runnable = has_hc;
        if (!runnable) {
          std::vector<double> nvs = nv;
          for (std::size_t j = 0; j < ns; ++j) {
            if (n0v_fed[j]) nvs[j] += REFORMING_JACOBIAN_TRACE * nt0;
          }
          for (std::size_t r = 0; r < rx.size() && !runnable; ++r) {
            if (std::find(act.begin(), act.end(), r) != act.end()) continue;
            bool fwd = true, bwd = true;
            for (const auto &p : rx[r].nu) {
              if (!(nvs[p.first] > 0.0)) (p.second < 0.0 ? fwd : bwd) = false;
            }
            runnable = fwd || bwd;
          }
        }
        got = !runnable;
      }
      if (!got) {
        throw std::runtime_error(
            "combustion_equilibrium_and_jacobian: no converged trace state for the one-sided "
            "derivative at an absent fed species (T_in = " +
            std::to_string(in.T) + " K, P = " + std::to_string(P) + " Pa)");
      }
    }
  }

  // Implicit function theorem on F(xi; n0, T, P) = 0 at the solution,
  // n = n0 + N xi: J dxi = -(dF/dn0 dn0 + dF/dT dT + dF/dP dP), with
  //   dF_a/dn0_l = nu_al / n_l - dn_a / n_t,
  //   dF_a/dT    = -d ln K_a/dT = -sum_j nu_aj h_j/(R T^2),
  //   dF_a/dP    = dn_a / P.
  // A reaction left out of the solve (a species absent, so it cannot move)
  // keeps its extent: zero derivative.
  const std::size_t m = act.size();
  double nt_sol = 0.0;
  for (double v : nv) nt_sol += v;
  std::vector<double> dxi(m * N, 0.0);  // [a][q]
  if (m > 0) {
    std::vector<double> lnK(rx.size(), 0.0), F, J;
    joint_residual(rx, act, nv, P, lnK, F, &J);  // J only; F unused
    for (std::size_t a = 0; a < m; ++a) {
      const JointReaction &r = rx[act[a]];
      for (std::size_t l = 0; l < ns; ++l) dxi[a * N + l] = r.dn / nt_sol;
      double dlnK_dT = 0.0;
      for (const auto &p : r.nu) {
        dxi[a * N + p.first] -= p.second / nv[p.first];
        dlnK_dT += p.second * h_RT(p.first, T_ad) / T_ad;
      }
      dxi[a * N + iT] = dlnK_dT;      // -dF/dT
      dxi[a * N + iP] = -r.dn / P;    // -dF/dP
    }
    // Scaled symmetrically (s_a = 1/sqrt(J_aa)): a trace species' rows carry
    // 1/n_j, up to 1e30 at the Jacobian's trace state above.
    // A subnormal key (subnormal_key_scale) has J_aa = inf and an inf
    // -nu_k / n_k on its own column of the right-hand side; both are taken
    // scaled, from the finite forms.
    std::vector<double> sc(m);
    for (std::size_t a = 0; a < m; ++a) sc[a] = 1.0 / std::sqrt(std::max(J[a * m + a], 1e-300));
    for (std::size_t a = 0; a < m; ++a) {
      const bool inf_diag = !std::isfinite(J[a * m + a]);
      double s_over_nk = 0.0;
      if (inf_diag) sc[a] = subnormal_key_scale(rx[act[a]], nv, nt_sol, &s_over_nk);
      for (std::size_t q = 0; q < N; ++q) dxi[a * N + q] *= sc[a];
      if (inf_diag) {
        const auto &kp = rx[act[a]].nu[rx[act[a]].key];
        dxi[a * N + kp.first] = sc[a] * rx[act[a]].dn / nt_sol - kp.second * s_over_nk;
      }
    }
    for (std::size_t a = 0; a < m; ++a) {
      for (std::size_t b = 0; b < m; ++b) J[a * m + b] *= sc[a] * sc[b];
      if (!std::isfinite(J[a * m + a])) J[a * m + a] = 1.0;
    }
    solve_dense(J, dxi, m, N);
    for (std::size_t a = 0; a < m; ++a) {
      for (std::size_t q = 0; q < N; ++q) dxi[a * N + q] *= sc[a];
    }
  }
  std::vector<DV> n(ns);
  for (std::size_t j = 0; j < ns; ++j) {
    n[j] = dv_const(nv[j], N);
    n[j].d[j] = 1.0;
  }
  for (std::size_t a = 0; a < m; ++a) {
    for (const auto &p : rx[act[a]].nu) {
      for (std::size_t q = 0; q < N; ++q) n[p.first].d[q] += p.second * dxi[a * N + q];
    }
  }
  const DV T = [&]() {
    DV t = dv_const(T_ad, N);
    t.d[iT] = 1.0;
    return t;
  }();
  // Stage 3: the adiabatic root F(T) = h(T, X_out)/MW_out - H_in = 0.
  DV h_mol = dv_const(0.0, N);
  DV mw = dv_const(0.0, N);
  for (std::size_t j = 0; j < ns; ++j) {
    const DV hj = dv_lift(h_species(j, T_ad), cp_species(j, T_ad), T);
    h_mol = h_mol + n[j] * hj;
    mw = mw + species_molar_mass(j) * n[j];
  }
  const DV H_mass = h_mol / mw;  // n/nt cancels in h/MW
  const double F_T = H_mass.d[iT];
  const double mw_in = mwmix(in.X);
  const double h_in = h(in.T, in.X);
  const double H_in = h_in / mw_in;

  out.dT_dT_in = cp(in.T, in.X) / mw_in / F_T;
  out.dT_dP = -H_mass.d[iP] / F_T;
  out.dT_dX.assign(ns, 0.0);
  for (std::size_t k = 0; k < ns; ++k) {
    double dF = 0.0;
    for (std::size_t j = 0; j < ns; ++j) dF += H_mass.d[j] * cj.dXb_dX[j][k];
    const double dH = (h_species(k, in.T) - H_in * species_molar_mass(k)) / mw_in;
    out.dT_dX[k] = (dH - dF) / F_T;
  }

  // Product mole fractions: n_out(n0, T, P), T = T_ad(X, T_in, P).
  std::vector<double> nvals(ns);
  double ntv = 0.0;
  for (std::size_t j = 0; j < ns; ++j) {
    nvals[j] = n[j].v;
    ntv += n[j].v;
  }
  // Rewrites dn in place: a fresh std::vector<double>(ns) here tripped GCC 14's
  // -Wstringop-overflow at -O3 (a false positive), failing the manylinux wheel.
  auto dX_from_dn = [&](std::vector<double> dn) {
    double dN = 0.0;
    for (double x : dn) dN += x;
    for (std::size_t j = 0; j < dn.size(); ++j) dn[j] = (dn[j] - nvals[j] / ntv * dN) / ntv;
    return dn;
  };
  out.dXout_dX.assign(ns, std::vector<double>(ns, 0.0));
  for (std::size_t k = 0; k < ns; ++k) {
    std::vector<double> dn(ns, 0.0);
    for (std::size_t j = 0; j < ns; ++j) {
      double v = n[j].d[iT] * out.dT_dX[k];
      for (std::size_t l = 0; l < ns; ++l) v += n[j].d[l] * cj.dXb_dX[l][k];
      dn[j] = v;
    }
    const auto dX = dX_from_dn(std::move(dn));
    for (std::size_t j = 0; j < ns; ++j) out.dXout_dX[j][k] = dX[j];
  }
  std::vector<double> dn_T(ns), dn_P(ns);
  for (std::size_t j = 0; j < ns; ++j) {
    dn_T[j] = n[j].d[iT] * out.dT_dT_in;
    dn_P[j] = n[j].d[iT] * out.dT_dP + n[j].d[iP];
  }
  out.dXout_dT_in = dX_from_dn(std::move(dn_T));
  out.dXout_dP = dX_from_dn(std::move(dn_P));
  return out;
}

} // namespace combaero
