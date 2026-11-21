// Isentropic stagnation/static conversion utilities.
// Uses the variable-cp thermo backend (thermo.h) throughout.

#include "../include/stagnation.h"
#include "../include/composition.h"
#include "../include/compressible.h"
#include "../include/thermo.h"
#include "../include/transport.h"
#include "../include/input_validation.h"
#include "../include/inner_solve.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace combaero {

// -------------------------------------------------------------
// Kinetic energy / stagnation enthalpy
// -------------------------------------------------------------

double h0_from_static(double h_static_J_per_kg, double v) {
    return h_static_J_per_kg + 0.5 * v * v;
}

double v_from_h0(double h0, double h_static_J_per_kg) {
    double dh = h0 - h_static_J_per_kg;
    return (dh > 0.0) ? std::sqrt(2.0 * dh) : 0.0;
}

// -------------------------------------------------------------
// Mach number
// -------------------------------------------------------------

double mach_number(double v, double T, const std::vector<double>& X) {
    double a = speed_of_sound(T, X);
    if (!validation::is_finite_positive(a)) {
        throw std::invalid_argument("mach_number: speed of sound is non-positive");
    }
    return v / a;
}

// -------------------------------------------------------------
// Stagnation temperature
// -------------------------------------------------------------

// Helper: h_mass [J/kg] at temperature T for composition X
static double h_mass_at(double T, const std::vector<double>& X) {
    double mw_g = mwmix(X);  // g/mol
    return h(T, X) / mw_g * 1000.0;  // J/mol / (g/mol) * 1000 = J/kg
}

// Helper: cp_mass [J/(kg*K)] at temperature T for composition X
static double cp_mass_at(double T, const std::vector<double>& X) {
    double mw_g = mwmix(X);  // g/mol
    return cp(T, X) / mw_g * 1000.0;  // J/(mol*K) / (g/mol) * 1000 = J/(kg*K)
}

double Tt_from_static_v(double T, double v, const std::vector<double>& X,
                         double tol, std::size_t max_iter) {
    if (!validation::is_finite_positive(T)) {
        throw std::invalid_argument("Tt_from_static_v: T must be positive");
    }
    if (!validation::is_finite_nonnegative(v)) {
        throw std::invalid_argument("Tt_from_static_v: v must be non-negative");
    }

    double h_static = h_mass_at(T, X);
    double h0 = h_static + 0.5 * v * v;

    // Newton iteration: find T0 such that h(T0) = h0
    // F(T0) = h(T0) - h0 = 0,  dF/dT0 = cp(T0)
    double T0 = T + 0.5 * v * v / cp_mass_at(T, X);  // ideal-gas initial guess

    // No floor (#935): the root is at or above T > 0 (h0 >= h(T), h rising),
    // and steps of at most half the iterate keep it positive. The 10 K floor
    // returned 10 K for a static T below it at low speed, not a root.
    double F = 0.0;
    for (std::size_t it = 0; it < max_iter; ++it) {
        F = h_mass_at(T0, X) - h0;
        double dF = cp_mass_at(T0, X);

        if (std::abs(dF) < 1e-30) break;

        double dT = -F / dF;
        double max_step = 0.5 * T0;
        if (std::abs(dT) > max_step) {
            dT = (dT > 0.0) ? max_step : -max_step;
        }
        T0 += dT;

        // Absolute below 1 K: h carries the heat of formation, so its
        // rounding resolves T only to ~1e-12 K (a relative test cycled there).
        if (std::abs(dT) < tol * (1.0 + std::abs(T0))) {
            return T0;
        }
    }
    inner_solve::throw_no_convergence("Tt_from_static_v", max_iter, "h residual [J/kg]", F);
}

double Tt_from_static(double T, double M, const std::vector<double>& X,
                       double tol, std::size_t max_iter) {
    if (!validation::is_finite_nonnegative(M)) {
        throw std::invalid_argument("Tt_from_static: M must be non-negative");
    }
    double a = speed_of_sound(T, X);
    double v = M * a;
    return Tt_from_static_v(T, v, X, tol, max_iter);
}

double T_from_stagnation(double Tt, double M, const std::vector<double>& X,
                          double tol, std::size_t max_iter) {
    if (!validation::is_finite_positive(Tt)) {
        throw std::invalid_argument("T_from_stagnation: Tt must be positive");
    }
    if (!validation::is_finite_nonnegative(M)) {
        throw std::invalid_argument("T_from_stagnation: M must be non-negative");
    }
    if (M == 0.0) return Tt;

    double h0 = h_mass_at(Tt, X);

    // Newton iteration: find T such that h(T) + M^2*a(T)^2/2 = h0
    // F(T) = h(T) + M^2*a(T)^2/2 - h0
    // dF/dT = cp(T) + (M^2/2) d(a^2)/dT, with a^2 = gamma R T and
    // d(a^2)/dT = R (gamma + T dgamma/dT), dgamma/dT = -R_u cp'/cv^2 (molar),
    // exact as fanno_mach.cpp's sound_speed_sq_dT. Dropping T dgamma/dT
    // made the convergence linear: it stopped ~1e-8 T short (#868 R3).
    double T = Tt / (1.0 + 0.2 * M * M);  // calorically perfect gas initial guess
    const double R_specific = specific_gas_constant(X);  // J/(kg*K)

    // No floor (#935): F(0) = h(0 K) - h0 < 0 < F(Tt), so the root is
    // positive, and steps of at most half the iterate keep it so. The 10 K
    // floor returned 10 K for a hypersonic static T below it, not a root.
    double F = 0.0;
    for (std::size_t it = 0; it < max_iter; ++it) {
        double a_val = speed_of_sound(T, X);
        double v = M * a_val;
        F = h_mass_at(T, X) + 0.5 * v * v - h0;
        const double cp_mol = cp(T, X);
        const double cv_mol = cp_mol - thermo::R_GAS;
        const double gamma = cp_mol / cv_mol;
        const double dgamma_dT = -thermo::R_GAS * dcp_dT(T, X) / (cv_mol * cv_mol);
        const double da2_dT = R_specific * (gamma + T * dgamma_dT);
        double dF = cp_mass_at(T, X) + 0.5 * M * M * da2_dT;

        if (std::abs(dF) < 1e-30) break;

        double dT = -F / dF;
        double max_step = 0.5 * T;
        if (std::abs(dT) > max_step) {
            dT = (dT > 0.0) ? max_step : -max_step;
        }
        T += dT;

        if (std::abs(dT) < tol * (1.0 + std::abs(T))) {  // absolute below 1 K, as above
            return T;
        }
    }
    inner_solve::throw_no_convergence("T_from_stagnation", max_iter, "h residual [J/kg]", F);
}

// -------------------------------------------------------------
// Stagnation pressure (isentropic)
// -------------------------------------------------------------

double Pt_from_static(double P, double T, double M,
                       const std::vector<double>& X,
                       double tol, std::size_t max_iter) {
    if (!validation::is_finite_positive(P)) {
        throw std::invalid_argument("Pt_from_static: P must be positive");
    }
    if (!validation::is_finite_nonnegative(M)) {
        throw std::invalid_argument("Pt_from_static: M must be non-negative");
    }
    if (M == 0.0) return P;

    constexpr double P_REF = 101325.0;
    double s0 = s(T, X, P, P_REF);
    double T0 = Tt_from_static(T, M, X, tol, max_iter);

    double R_specific = specific_gas_constant(X);  // J/(kg*K)
    double mw_kg = mwmix(X) / 1000.0;  // kg/mol

    double s_T0_Pref = s(T0, X, P_REF, P_REF) / mw_kg;  // J/(kg*K)
    double s0_mass   = s0 / mw_kg;                        // J/(kg*K)
    double ln_ratio = (s_T0_Pref - s0_mass) / R_specific;
    return P_REF * std::exp(ln_ratio);
}

double P_from_stagnation(double Pt, double Tt, double M,
                          const std::vector<double>& X,
                          double tol, std::size_t max_iter) {
    if (!validation::is_finite_positive(Pt)) {
        throw std::invalid_argument("P_from_stagnation: Pt must be positive");
    }
    if (!validation::is_finite_nonnegative(M)) {
        throw std::invalid_argument("P_from_stagnation: M must be non-negative");
    }
    if (M == 0.0) return Pt;

    constexpr double P_REF = 101325.0;
    double s0 = s(Tt, X, Pt, P_REF);
    double T_static = T_from_stagnation(Tt, M, X, tol, max_iter);

    double R_specific = specific_gas_constant(X);  // J/(kg*K)
    double mw_kg = mwmix(X) / 1000.0;  // kg/mol

    double s_T_Pref = s(T_static, X, P_REF, P_REF) / mw_kg;  // J/(kg*K)
    double s0_mass  = s0 / mw_kg;                              // J/(kg*K)
    double ln_ratio = (s_T_Pref - s0_mass) / R_specific;
    return P_REF * std::exp(ln_ratio);
}

// -------------------------------------------------------------
// Adiabatic wall temperature
// -------------------------------------------------------------

double recovery_factor(double Pr, bool turbulent) {
    if (!validation::is_finite_positive(Pr)) {
        throw std::invalid_argument("recovery_factor: Pr must be positive");
    }
    return turbulent ? std::cbrt(Pr) : std::sqrt(Pr);
}

double T_adiabatic_wall(double T_static, double v,
                         double T, double P, const std::vector<double>& X,
                         bool turbulent) {
    // The 10 K floor below turned a NaN T_static into 10 K at v = 0 (#769).
    if (!validation::is_finite(T_static) || !validation::is_finite(v)) {
        throw std::invalid_argument("T_adiabatic_wall: T_static and v must be finite");
    }
    double T_safe = std::max(T_static, 10.0);
    double v_abs = std::abs(v);
    if (v_abs < 1e-12) return T_safe;


    double Pr = prandtl(T, P, X);
    double r = recovery_factor(Pr, turbulent);
    double cp_val = cp_mass_at(T, X);

    return T_static + r * 0.5 * v * v / cp_val;
}

double T_adiabatic_wall_mach(double T_static, double M,
                               double T, double P, const std::vector<double>& X,
                               bool turbulent) {
    if (!validation::is_finite_nonnegative(M)) {
        throw std::invalid_argument("T_adiabatic_wall_mach: M must be non-negative");
    }
    if (M == 0.0) return T_static;

    double a = speed_of_sound(T, X);
    double v = M * a;
    return T_adiabatic_wall(T_static, v, T, P, X, turbulent);
}

std::tuple<double, double> stagnation_from_static(
    double T, double P, double v,
    const std::vector<double>& X,
    double tol, std::size_t max_iter)
{
    if (!validation::is_finite_positive(T))
        throw std::invalid_argument("stagnation_from_static: T must be positive");
    if (!validation::is_finite_positive(P))
        throw std::invalid_argument("stagnation_from_static: P must be positive");
    if (!validation::is_finite_nonnegative(v))
        throw std::invalid_argument("stagnation_from_static: v must be non-negative");
    if (v == 0.0) return {T, P};

    // One Newton loop for Tt.
    double Tt = Tt_from_static_v(T, v, X, tol, max_iter);

    // Pt from isentropic entropy conservation: s(T, P) = s(Tt, Pt).
    // Bypasses Pt_from_static to avoid its internal Tt_from_static call,
    // which would repeat the Newton loop and re-evaluate speed_of_sound(T, X).
    constexpr double P_REF = 101325.0;
    const double mw_kg  = mwmix(X) / 1000.0;
    const double R_spec = specific_gas_constant(X);
    const double s_static = s(T,  X, P,     P_REF) / mw_kg;
    const double s_stag   = s(Tt, X, P_REF, P_REF) / mw_kg;
    double Pt = P_REF * std::exp((s_stag - s_static) / R_spec);

    return {Tt, Pt};
}

namespace {
// Integral of f over [a, b] by 8-point Gauss-Legendre, split at the
// NASA-9 range boundaries (1000 K, 6000 K) where cp has a small jump (#549),
// and at the data's lower end (200 K), where the low tail (#606) joins C2:
// cp''' jumps there, and one rule across it missed the integral by enough
// that its slope in the limit differed from the integrand, which the
// analytic Jacobians use (a momentum chamber at T_static 168 K < 200 K <
// Tt 317 K, 2.2e-6 off its own d/dP; Pt moved by up to 2.2e-6, #1255).
// Used for h(T0) - h(T) and s(T0) - s(T) at the SAME composition: the direct
// differences cancel ~4 significant digits near T0 and left ~1e-9 Pa of noise
// in P0, which finite-difference Jacobian checks with small steps read as
// derivative errors (#579 on Linux).
template <class F>
double gl8_split(F f, double a, double b) {
  static const double xg[4] = {0.1834346424956498, 0.5255324099163290,
                               0.7966664774136267, 0.9602898564975363};
  static const double wg[4] = {0.3626837833783620, 0.3137066458778873,
                               0.2223810344533745, 0.1012285362903763};
  auto seg = [&](double lo, double hi) {
    const double c = 0.5 * (lo + hi), r = 0.5 * (hi - lo);
    double acc = 0.0;
    for (int i = 0; i < 4; ++i) {
      acc += wg[i] * (f(c - r * xg[i]) + f(c + r * xg[i]));
    }
    return acc * r;
  };
  const double sign = b >= a ? 1.0 : -1.0;
  double lo = std::min(a, b), hi = std::max(a, b);
  double total = 0.0;
  for (double bp : {thermo::CP_INVERSE_T_DATA, 1000.0, 6000.0}) {
    if (lo < bp && bp < hi) {
      total += seg(lo, bp);
      lo = bp;
    }
  }
  total += seg(lo, hi);
  return sign * total;
}

}  // namespace

ChamberStaticState chamber_static_state(double Tt, double P, double G,
                                        const std::vector<double> &X) {
  // A non-positive or non-finite Tt (garbage from an upstream state at a
  // solver probe) has no static state: refuse it (#868 R4). The Newton below
  // either ran to its cap on it or, when it happened to converge, returned a
  // negative T. An invalid argument (ValueError, #928), unlike the cap exit;
  // the solver rejects the probe either way.
  if (!validation::is_finite_positive(Tt)) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.6e", Tt);
    throw std::invalid_argument(
        std::string("chamber_static_state: Tt must be finite and > 0 K (got ") + buf + ")");
  }
  const double R = specific_gas_constant(X);
  const double k = (G * R / P) * (G * R / P);  // v^2 = k T^2
  auto cp_of = [&](double t) { return cp_mass(t, X); };
  auto cp_over_T = [&](double t) { return cp_mass(t, X) / t; };
  // f(T) = [h(Tt) - h(T)] - k T^2 / 2, the enthalpy difference integrated.
  double T = Tt;
  for (int it = 0;; ++it) {
    const double f = gl8_split(cp_of, T, Tt) - 0.5 * k * T * T;
    const double fT = -cp_mass(T, X) - k * T;
    const double step = f / fT;
    T -= step;
    if (std::abs(step) <= 1e-13 * Tt) {
      break;
    }
    if (it + 1 == inner_solve::kStaticStateMaxIter) {
      if (inner_solve::cap_exit_is_root(step, Tt)) {
        break;
      }
      inner_solve::throw_no_convergence("chamber_static_state", inner_solve::kStaticStateMaxIter,
                                        "step [K]", step);
    }
  }
  ChamberStaticState out{};
  out.T = T;
  const double cp_T = cp_mass(T, X);
  const double cp_Tt = cp_mass(Tt, X);
  const double fT = -cp_T - k * T;
  // df/dx for x = Tt, P, G; dT/dx = -(df/dx) / fT.
  out.dT_dTt = -cp_Tt / fT;
  out.dT_dP = -(k * T * T / P) / fT;
  out.dT_dG = (G > 0.0) ? -(-k * T * T / G) / fT : 0.0;
  // ln(Pt / P) = (s(Tt) - s(T)) / R = (1/R) integral of cp / T.
  out.lnPtP = gl8_split(cp_over_T, T, Tt) / R;
  const double dln_dTs = -cp_T / (R * T);
  out.dln_dTt = cp_Tt / (R * Tt) + dln_dTs * out.dT_dTt;
  out.dln_dP = dln_dTs * out.dT_dP;
  out.dln_dG = dln_dTs * out.dT_dG;
  const double a = speed_of_sound(T, X);
  out.M = (a > 1e-9) ? G * R * T / (P * a) : 0.0;
  return out;
}

ChamberStaticStateDX chamber_static_state_dX(double Tt, double P, double G,
                                            const std::vector<double> &X) {
  const std::size_t n = X.size();
  ChamberStaticStateDX out;
  out.dT_dX.assign(n, 0.0);
  out.dlnPtP_dX.assign(n, 0.0);
  const ChamberStaticState cs = chamber_static_state(Tt, P, G, X);
  const double T = cs.T;
  const double R = specific_gas_constant(X);
  const double MW = mwmix(X) / 1000.0;  // kg/mol
  const double Ru = R * MW;
  const double cpT = cp_mass(T, X);
  const double v = G * R * T / P;
  const double fT = -cpT - (G * R / P) * (G * R / P) * T;
  double dH = 0.0;
  for (std::size_t j = 0; j < n; ++j) {
    dH += X[j] * (h_species(j, Tt) - h_species(j, T));
  }
  dH /= MW;
  for (std::size_t j = 0; j < n; ++j) {
    const double W = species_molar_mass(j) / 1000.0;
    // h_mass difference and R = R_u / MW both move with X_j.
    const double fX = (h_species(j, Tt) - h_species(j, T) - dH * W) / MW + v * v * W / MW;
    out.dT_dX[j] = -fX / fT;
    out.dlnPtP_dX[j] =
        (s_species(j, Tt) - s_species(j, T)) / Ru - cpT / (R * T) * out.dT_dX[j];
  }
  return out;
}

SonicStaticState sonic_static_state(double Tt, double G,
                                    const std::vector<double> &X) {
  const double R = specific_gas_constant(X);
  const double mw_g = mwmix(X);  // g/mol
  // cp, gamma and their T derivatives on the mass basis.
  auto cp_m = [&](double t) { return cp_mass(t, X); };
  auto gamma_and_d = [&](double t, double &dg) {
    const double c = cp_mass(t, X);
    const double dc = J_per_mol_to_J_per_kg(dcp_dT(t, X), mw_g);
    const double cv = c - R;
    dg = -R * dc / (cv * cv);
    return c / cv;
  };
  // F(T) = [h(Tt) - h(T)] - a(T)^2 / 2, a^2 = gamma R T; F falls in T.
  double dg0 = 0.0;
  const double g0 = gamma_and_d(Tt, dg0);
  double T = 2.0 * Tt / (g0 + 1.0);
  double dFdT = -1.0;
  for (int it = 0;; ++it) {
    double dg = 0.0;
    const double g = gamma_and_d(T, dg);
    const double F = gl8_split(cp_m, T, Tt) - 0.5 * g * R * T;
    dFdT = -cp_mass(T, X) - 0.5 * R * (g + T * dg);
    const double step = F / dFdT;
    T -= step;
    if (std::abs(step) <= 1e-13 * Tt) {
      break;
    }
    if (it + 1 == inner_solve::kStaticStateMaxIter) {
      if (inner_solve::cap_exit_is_root(step, Tt)) {
        break;
      }
      inner_solve::throw_no_convergence("sonic_static_state", inner_solve::kStaticStateMaxIter,
                                        "step [K]", step);
    }
  }
  double dg = 0.0;
  const double g = gamma_and_d(T, dg);
  dFdT = -cp_mass(T, X) - 0.5 * R * (g + T * dg);
  SonicStaticState out{};
  out.T = T;
  out.dT_dTt = cp_mass(Tt, X) / -dFdT;
  const double Gabs = std::abs(G);
  out.P = Gabs * std::sqrt(R * T / g);
  // d ln P* / d T* = 1/(2 T*) - gamma'/(2 gamma).
  const double dlnP_dT = 0.5 / T - 0.5 * dg / g;
  out.dP_dTt = out.P * dlnP_dT * out.dT_dTt;
  out.dP_dG = (Gabs > 0.0) ? out.P / Gabs * (G >= 0.0 ? 1.0 : -1.0) : std::sqrt(R * T / g);
  return out;
}

DuctStaticState duct_static_state(double Tt, double Pt, double G,
                                  const std::vector<double> &X) {
  // An unphysical iterate is evaluated at the floors, finite; the partials
  // in a floored argument are then zero.
  const bool T0_floored = !(Tt > thermo::T_FLOOR);
  const bool P0_floored = !(Pt > thermo::P_FLOOR);
  Tt = T0_floored ? thermo::T_FLOOR : Tt;
  Pt = P0_floored ? thermo::P_FLOOR : Pt;
  const double R = specific_gas_constant(X);
  const double lnP0 = std::log(Pt);
  G = std::abs(G);
  // u = ln(P / Pt), P = Pt e^u: L(u) = lnPtP(Tt, P, G) + u, increasing in u on
  // [ln(P* / Pt), 0]. Relative to Pt, not ln P itself, so that u = 0 is P = Pt
  // to the bit and P <= Pt throughout: exp(log(Pt)) came back one ulp ABOVE
  // Pt at G = 0, a static above its own total, and the in-pipe orifice's
  // drive p1 - p then exceeded the jet's p0 - p at zero flow (#997).
  auto P_of = [&](double u) { return Pt * std::exp(u); };
  auto eval = [&](double u, ChamberStaticState &cs) {
    cs = chamber_static_state(Tt, P_of(u), G, X);
    return cs.lnPtP + u;
  };
  DuctStaticState out{};
  ChamberStaticState cs{};
  const double P_star = critical_pressure_ratio(Tt, Pt, X) * Pt;
  double u_lo = std::log(P_star / Pt);
  double u_hi = 0.0;
  if (eval(u_lo, cs) > 0.0) {
    // Beyond choke: the sonic state, independent of G. T* from the energy
    // balance at M = 1 (sonic_static_state, Newton to 1e-13 Tt) and P* on
    // the isentrope through (Tt, Pt), ln(Pt / P*) = (s(Tt) - s(T*)) / R --
    // the subsonic branch's own entropy integral, so the two meet at G*.
    // P_star above (golden section on G(P), ~1e-8 relative) only decides
    // the branch: as the value its noise read as a derivative error, and
    // its Tt dependence (cp varies) was dropped, which put the in-pipe
    // orifice's Tt column, and the flow columns it reaches, wrong past the
    // approach pipe's sonic limit (#958).
    const SonicStaticState son = sonic_static_state(Tt, G, X);
    const double T_s = son.T;
    const double cp_Tt = cp_mass(Tt, X);
    const double cp_Ts = cp_mass(T_s, X);
    auto cp_over_T = [&](double t) { return cp_mass(t, X) / t; };
    const double P_s = Pt * std::exp(-gl8_split(cp_over_T, T_s, Tt) / R);
    // d ln P* / d Tt = -(cp(Tt) / Tt - cp(T*) / T* dT*/dTt) / R.
    const double dlnP_dTt = -(cp_Tt / Tt - cp_Ts / T_s * son.dT_dTt) / R;
    out.T = T_s;
    out.P = P_s;
    out.rho = P_s / (R * T_s);
    out.M = 1.0;
    out.choked = true;
    out.G_eff = G;
    out.dGeff_dG = 1.0;
    out.dP_dPt = P_s / Pt;
    out.dP_dTt = P_s * dlnP_dTt;
    out.dT_dTt = son.dT_dTt;
    out.dT_dPt = 0.0;
    out.drho_dTt = out.rho * (dlnP_dTt - out.dT_dTt / out.T);
    out.drho_dPt = out.rho / Pt;
    if (T0_floored) {
      out.dT_dTt = out.dP_dTt = out.drho_dTt = 0.0;
    }
    if (P0_floored) {
      out.dT_dPt = out.dP_dPt = out.drho_dPt = 0.0;
    }
    return out;
  }
  double u = u_hi;
  double L = eval(u, cs);
  for (int it = 0; L != 0.0; ++it) {
    if (L > 0.0) {
      u_hi = u;
    } else {
      u_lo = u;
    }
    const double dL = 1.0 + P_of(u) * cs.dln_dP;
    double u_new = (dL > 0.0) ? u - L / dL : 0.5 * (u_lo + u_hi);
    if (!(u_new > u_lo && u_new < u_hi)) {
      u_new = 0.5 * (u_lo + u_hi);
    }
    const double step = u_new - u;
    const bool done = std::abs(step) <= 1e-15 * std::abs(lnP0) + 1e-15;
    u = u_new;
    L = eval(u, cs);
    if (done) {
      break;
    }
    if (it + 1 == inner_solve::kDuctStaticStateMaxIter) {
      // Cap rule (#868 R9): the last step closed to rounding is a root;
      // anything wider is not, and the implicit derivatives below would
      // describe a non-root.
      if (inner_solve::cap_exit_is_root(step, std::max(1.0, std::abs(u + lnP0)))) {
        break;
      }
      inner_solve::throw_no_convergence("duct_static_state", inner_solve::kDuctStaticStateMaxIter,
                                        "ln P step", step);
    }
  }
  const double P = P_of(u);
  const double dL_du = 1.0 + P * cs.dln_dP;
  // In v = ln P = u + ln Pt, L(v) = lnPtP + v - ln Pt: dv/dx = -(dL/dx) /
  // (dL/dv), dL/dT0 = dln_dTt, dL/dP0 = -1/Pt, dL/dG = dln_dG, dL/dv = dL/du.
  const double du_dT0 = -cs.dln_dTt / dL_du;
  const double du_dP0 = (1.0 / Pt) / dL_du;
  const double du_dG = -cs.dln_dG / dL_du;
  out.P = P;
  out.T = cs.T;
  out.rho = P / (R * cs.T);
  out.M = cs.M;
  out.choked = false;
  out.dP_dTt = P * du_dT0;
  out.dP_dPt = P * du_dP0;
  out.dP_dG = P * du_dG;
  out.dT_dTt = cs.dT_dTt + cs.dT_dP * out.dP_dTt;
  out.dT_dPt = cs.dT_dP * out.dP_dPt;
  out.dT_dG = cs.dT_dG + cs.dT_dP * out.dP_dG;
  out.drho_dTt = out.rho * (out.dP_dTt / P - out.dT_dTt / out.T);
  out.drho_dPt = out.rho * (out.dP_dPt / P - out.dT_dPt / out.T);
  out.drho_dG = out.rho * (out.dP_dG / P - out.dT_dG / out.T);
  out.G_eff = G;
  out.dGeff_dG = 1.0;
  if (T0_floored) {
    out.dT_dTt = out.dP_dTt = out.drho_dTt = 0.0;
  }
  if (P0_floored) {
    out.dT_dPt = out.dP_dPt = out.drho_dPt = 0.0;
  }
  return out;
}

DuctStaticState duct_htc_static_state(double Tt, double Pt, double G,
                                      const std::vector<double> &X) {
  G = std::abs(G);
  if (!(Tt > thermo::T_FLOOR) || !(Pt > thermo::P_FLOOR)) {
    return duct_static_state(Tt, Pt, G, X);
  }
  auto sonic_flux = [&](double T, double P) {
    return mass_flux_isentropic(T, P, critical_pressure_ratio(T, P, X) * P, X);
  };
  const double Gs = sonic_flux(Tt, Pt);
  const double w = kDuctHtcSonicMargin;
  const double a = 1.0 - 2.0 * w;
  const double g = G / Gs;
  if (g <= a) {
    return duct_static_state(Tt, Pt, G, X);
  }
  // g_eff = g - (g - a)^2 / (4 w) on [a, 1]: value and slope meet g at a and
  // the cap 1 - w (slope 0) at 1.
  double ge = 1.0 - w;
  double dge = 0.0;
  if (g < 1.0) {
    ge = g - (g - a) * (g - a) / (4.0 * w);
    dge = 1.0 - (g - a) / (2.0 * w);
  }
  DuctStaticState s = duct_static_state(Tt, Pt, ge * Gs, X);
  // G_eff = Gs(Tt, Pt) ge(G / Gs): d/dG = dge, d/dGs = ge - g dge. The sonic
  // flux is proportional to Pt at fixed Tt (ideal gas: P*/Pt depends on Tt
  // only); its Tt slope is differenced (smooth, step 1e-6).
  const double hT = 1e-6 * Tt;
  const double dGs_dT0 = (sonic_flux(Tt + hT, Pt) - sonic_flux(Tt - hT, Pt)) / (2.0 * hT);
  const double dGs_dP0 = Gs / Pt;
  const double dGe_dGs = ge - g * dge;
  s.G_eff = ge * Gs;
  s.dGeff_dG = dge;
  s.dGeff_dTt = dGe_dGs * dGs_dT0;
  s.dGeff_dPt = dGe_dGs * dGs_dP0;
  // Chain the state's G partials through G_eff.
  auto chain = [&](double &d_T0, double &d_P0, double &d_G) {
    d_T0 += d_G * s.dGeff_dTt;
    d_P0 += d_G * s.dGeff_dPt;
    d_G *= s.dGeff_dG;
  };
  chain(s.dT_dTt, s.dT_dPt, s.dT_dG);
  chain(s.dP_dTt, s.dP_dPt, s.dP_dG);
  chain(s.drho_dTt, s.drho_dPt, s.drho_dG);
  return s;
}

DuctStaticStateDY duct_htc_static_state_dY(double Tt, double Pt, double G,
                                           const std::vector<double> &X,
                                           const DuctStaticState &state) {
  const std::size_t n = X.size();
  DuctStaticStateDY out;
  out.dT_dY.assign(n, 0.0);
  out.dP_dY.assign(n, 0.0);
  out.drho_dY.assign(n, 0.0);
  out.dGeff_dY.assign(n, 0.0);
  G = std::abs(G);
  if (!(Tt > thermo::T_FLOOR) || !(Pt > thermo::P_FLOOR)) {
    return out;
  }
  // The caller's solved state. Below the bend it IS the plain state at G
  // (G_eff = G exactly) and nothing is solved again; in the bend, the raw
  // state at G_eff gives its own G partials, before the bend's chain.
  const DuctStaticState &s = state;
  const bool bent = s.G_eff != G;
  const DuctStaticState raw = bent ? duct_static_state(Tt, Pt, s.G_eff, X) : s;
  auto sonic_flux = [&](const std::vector<double> &Xv) {
    return mass_flux_isentropic(Tt, Pt, critical_pressure_ratio(Tt, Pt, Xv) * Pt, Xv);
  };
  const double w = kDuctHtcSonicMargin;
  const double a = 1.0 - 2.0 * w;
  const double T = s.T;
  const double P = s.P;
  const double R = specific_gas_constant(X);
  const double MW = mwmix(X) / 1000.0;  // kg/mol
  const double Ru = R * MW;
  const double cpT = cp_mass(T, X);
  const double v = s.G_eff * R * T / P;
  const double kk = (s.G_eff * R / P) * (s.G_eff * R / P);
  // Mass-specific enthalpy drop, molar species terms over MW.
  // Species enthalpy and entropy drops Tt -> T, once.
  std::vector<double> dh_j(n), ds_j(n);
  double dH = 0.0;
  for (std::size_t j = 0; j < n; ++j) {
    dh_j[j] = h_species(j, Tt) - h_species(j, T);
    ds_j[j] = s_species(j, Tt) - s_species(j, T);
    dH += X[j] * dh_j[j];
  }
  dH /= MW;
  // Jacobian in (T, u = ln P): f_T, f_u; L_T, L_u = 1.
  const double fT = -cpT - kk * T;
  const double fu = kk * T * T;
  const double LT = -cpT / (R * T);
  const double det = fT - fu * LT;
  std::vector<double> dT_dX(n, 0.0), du_dX(n, 0.0), drho_dX(n, 0.0);
  for (std::size_t j = 0; j < n; ++j) {
    const double W = species_molar_mass(j) / 1000.0;
    // dR/dX_j = -R W/MW enters f through -(v^2 / R) dR.
    const double fX = (dh_j[j] - dH * W) / MW + v * v * W / MW;
    const double LX = ds_j[j] / Ru;
    dT_dX[j] = -(fX - fu * LX) / det;
    du_dX[j] = -(fT * LX - LT * fX) / det;
  }
  // The bent flux near choke: G_eff = ge G* with G*(Tt, Pt, X).
  std::vector<double> dGe_dX(n, 0.0);
  if (bent) {
    const double Gs = sonic_flux(X);
    const double g = G / Gs;
    double ge = 1.0 - w;
    double dge = 0.0;
    if (g < 1.0) {
      ge = g - (g - a) * (g - a) / (4.0 * w);
      dge = 1.0 - (g - a) / (2.0 * w);
    }
    const double dGe_dGs = ge - g * dge;
    for (std::size_t j = 0; j < n; ++j) {
      // One-sided off an absent species, as it enters.
      const double hX = 1e-6;
      std::vector<double> Xp = X, Xm = X;
      Xp[j] += hX;
      double dGs = 0.0;
      if (X[j] > hX) {
        Xm[j] -= hX;
        dGs = (sonic_flux(Xp) - sonic_flux(Xm)) / (2.0 * hX);
      } else {
        dGs = (sonic_flux(Xp) - Gs) / hX;
      }
      dGe_dX[j] = dGe_dGs * dGs;
      dT_dX[j] += raw.dT_dG * dGe_dX[j];
      du_dX[j] += raw.dP_dG / P * dGe_dX[j];
    }
  }
  std::vector<double> dP_dX(n, 0.0);
  for (std::size_t j = 0; j < n; ++j) {
    const double W = species_molar_mass(j) / 1000.0;
    dP_dX[j] = P * du_dX[j];
    // rho = P MW / (R_u T)
    drho_dX[j] = s.rho * (du_dX[j] - dT_dX[j] / T + W / MW);
  }
  out.dT_dY = mole_gradient_to_mass(X, dT_dX);
  out.dP_dY = mole_gradient_to_mass(X, dP_dX);
  out.drho_dY = mole_gradient_to_mass(X, drho_dX);
  out.dGeff_dY = mole_gradient_to_mass(X, dGe_dX);
  return out;
}

} // namespace combaero
