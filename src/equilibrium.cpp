#include "../include/equilibrium.h"
#include "../include/combustion.h"
#include "../include/inner_solve.h"
#include "thermo.h"
#include "composition.h"
#include "thermo_transport_data.h"
#include <algorithm>
#include <cmath>
#include <limits>
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
static std::vector<JointReaction> joint_reactions(const ReformingConfig &cfg,
                                                  bool with_reforming) {
  std::vector<JointReaction> rx;
  if (with_reforming) {
    for (std::size_t k = 0; k < cfg.hc_indices.size(); ++k) {
      const double nC = cfg.hc_C[k];
      const double m = nC + 0.5 * cfg.hc_H[k];
      rx.push_back(JointReaction{{{cfg.hc_indices[k], -1.0},
                                  {cfg.i_H2O, -nC},
                                  {cfg.i_CO, nC},
                                  {cfg.i_H2, m}},
                                 m - 1.0, 0});
    }
  }
  rx.push_back(JointReaction{
      {{cfg.i_CO, -1.0}, {cfg.i_H2O, -1.0}, {cfg.i_CO2, 1.0}, {cfg.i_H2, 1.0}}, 0.0, 2});
  return rx;
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
    std::vector<double> A(m * m), dx(m);
    for (std::size_t a = 0; a < m; ++a) {
      for (std::size_t b = 0; b < m; ++b) A[a * m + b] = s[a] * J[a * m + b] * s[b];
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
        double nk_new = nk * std::exp(d / nk);
        if (d > 0.0) {
          const double ln_cap = std::log(REFORMING_TRACE_FRACTION * nt);
          nk_new = std::max(nk + d, std::exp(std::min(std::log(nk) + d / nk, ln_cap)));
          if (nk >= REFORMING_TRACE_FRACTION * nt) nk_new = nk + d;
        }
        const double xi = (nk_new - nk) / kp.second;
        for (const auto &p : r.nu) out[p.first] += p.second * xi;
        out[kp.first] = nk_new;
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
        // Sufficient decrease of ||F||: the Newton direction descends it.
        if (std::isfinite(f_try) && f_try <= (1.0 - 1e-4 * alpha) * fnorm) {
          accepted = true;
          break;
        }
      }
      alpha *= 0.5;
    }
    if (!accepted) {
      // Round-off floor (or no descent): converged if within tolerance.
      return fnorm <= REFORMING_NEWTON_FTOL;
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
static bool element_potential_start(std::vector<double> &n, double T, double P,
                                    const std::vector<JointReaction> &rx,
                                    const std::vector<std::size_t> &act) {
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
  for (int it = 0; it < REFORMING_NEWTON_MAX_IT && !done; ++it) {
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
        const double akq = aj[q][el[k]] * nq[q];
        for (std::size_t i = 0; i < E; ++i) A[k * M + i] += akq * aj[q][el[i]];
        ak_n += akq;
        r[k] += akq * mu[q];
      }
      A[k * M + E] = ak_n;
      A[E * M + k] = ak_n;
      r[k] += b[el[k]] - ak_n;
    }
    A[E * M + E] = sum_n - N;
    r[E] = N - n_inert - sum_n;
    for (std::size_t q = 0; q < S; ++q) r[E] += nq[q] * mu[q];
    if (!solve_dense(A, r, M, 1)) return false;
    const double dlnN = r[E];
    std::vector<double> dl(S);
    double lam1 = 5.0 * std::abs(dlnN), dmax = 0.0;
    double lam2 = 1.0;
    for (std::size_t q = 0; q < S; ++q) {
      dl[q] = -mu[q] + dlnN;
      for (std::size_t i = 0; i < E; ++i) dl[q] += aj[q][el[i]] * r[i];
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
  }
  if (!done) return false;
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

// Joint reforming + WGS equilibrium at fixed T from n0. Rounds start it
// (REFORMING_START_PASSES, which also make the species of a reaction that can
// move present); Newton then solves the reactions whose species are all
// present. If some reaction was left out, a trial round checks that it cannot
// move; if it can, the round is taken and Newton re-solved. Returns whether
// the residuals met REFORMING_NEWTON_FTOL and no left-out reaction moves.
// active_out (optional) receives the reactions solved by Newton.
static bool solve_reforming_wgs_joint(const std::vector<double> &n0, double T, double P,
                                      const ReformingConfig &cfg, bool with_reforming,
                                      std::vector<double> &n_out,
                                      std::vector<JointReaction> *rx_out = nullptr,
                                      std::vector<std::size_t> *active_out = nullptr) {
  const std::vector<JointReaction> rx = joint_reactions(cfg, with_reforming);
  std::vector<double> lnK(rx.size());
  for (std::size_t r = 0; r < rx.size(); ++r) lnK[r] = joint_lnK(rx[r], T);

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

  bool conv = false;
  std::vector<std::size_t> act;
  for (int restart = 0; restart <= REFORMING_JOINT_RESTARTS; ++restart) {
    act = joint_active(rx, n_out);
    std::vector<double> n_try = n_out;
    conv = joint_newton(n_try, P, rx, act, lnK, REFORMING_NEWTON_ROUNDS_IT);
    if (!conv) {
      // Far from equilibrium: restart Newton from the element potentials.
      std::vector<double> n_ep = n_out;
      if (element_potential_start(n_ep, T, P, rx, act) &&
          joint_newton(n_ep, P, rx, act, lnK, REFORMING_NEWTON_MAX_IT)) {
        n_try = n_ep;
        conv = true;
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
  const double T_ad = solve_adiabatic_T(F, in.T, conv_T);

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
  solve_reforming_wgs_joint(n0, T_ad, P, cfg, has_hc, nv, &rx, &act);

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
    solve_dense(J, dxi, m, N);
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
