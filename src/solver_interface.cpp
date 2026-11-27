#include "../include/solver_interface.h"
#include "../include/tee_junction.h"
#include "../include/composition.h"
#include "../include/transport.h"
#include "../include/combustion.h"
#include "../include/compressible.h"
#include "../include/cooling_correlations.h"
#include "../include/equilibrium.h"
#include "../include/flow_drive_floor.h"
#include "../include/friction.h"
#include "../include/heat_transfer.h"
#include "../include/input_validation.h"
#include "../include/registry.h"
#include "../include/quadrature.h"
#include "../include/stagnation.h"
#include "../include/thermo.h"
#include "../include/math_constants.h"
#include "../include/orifice.h"
#include "../include/orifice_jobson.h"
#include "../include/thermo_transport_data.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <utility>

namespace combaero {
namespace solver {

std::tuple<double, double> mixer_heat_per_unit_flow(double mdot_tot) {
  const double m0 = kMixerHeatMdotFloor;
  const double a = std::abs(mdot_tot);
  const double sgn = mdot_tot >= 0.0 ? 1.0 : -1.0;
  if (a >= m0) {
    return {1.0 / a, -sgn / (a * a)};
  }
  const double u = a / m0;
  return {(4.0 * u * u - 3.0 * u * u * u) / m0, sgn * (8.0 * u - 9.0 * u * u) / (m0 * m0)};
}

namespace {

// T from a mixed specific enthalpy, for the network's mixing nodes. The gas is
// ideal down to 0 K (the thermo's low tail, #626), so above h(0 K) this is the
// exact inverse. Below it -- a solver probe or starting point that withdraws
// more heat than the stream holds -- no T > 0 exists, and the inverse would
// raise and abort the solve (#552 item 3, #607): T continues on h's own line
// through 0 K (cp is constant there, so h and cp at the returned T are the
// continuation's own, C1 with the inverse). The network reports a node at
// T <= 0 as holding no gas state.
double T_from_h_mass_continued(double h_mass, const std::vector<double>& X, double T_guess) {
  const double h_0 = combaero::h_mass(0.0, X);
  if (h_mass > h_0) {
    return combaero::calc_T_from_h_mass(h_mass, X, T_guess);
  }
  return (h_mass - h_0) / combaero::cp_mass(0.0, X);
}

// The streams a mixing node's state is made of (#1232): those arriving (m_i >
// 0), or every stream when none arrives (all m_i <= 0, the #1224 reverse
// mean). The mixer and the combustor kernels weight these alone, and take
// their total as the flow the node's heat is spread over. A stream at m_i = 0
// beside an arriving one is in the mix with weight 0: its partials are those
// of the side it arrives from, so a stream about to flow still moves the mix.
std::vector<bool> streams_in_mix(const std::vector<Stream> &streams) {
  bool any_arriving = false;
  for (const auto &s : streams) {
    any_arriving = any_arriving || s.m_dot > 0.0;
  }
  std::vector<bool> in_mix(streams.size(), true);
  if (any_arriving) {
    for (std::size_t i = 0; i < streams.size(); ++i) {
      in_mix[i] = streams[i].m_dot >= 0.0;
    }
  }
  return in_mix;
}
} // anonymous namespace

// -----------------------------------------------------------------------------
// 1. Incompressible Flow Components
// -----------------------------------------------------------------------------

std::tuple<double, double, double> orifice_mdot_and_jacobian(double dP, double rho,
                                                              double Cd, double area,
                                                              double beta) {
  // Regularized orifice equation: mdot = Cd * E * A * f_reg(dP, rho)
  // where E = 1/sqrt(1 - beta^4) is the velocity-of-approach factor.
  // We use a smooth sqrt to avoid the infinite gradient at dP = 0.
  const double E = velocity_of_approach_factor(beta);

  /// Smooth clamp: rho_eff ~ max(rho, L)
  constexpr double L = 1e-9;
  constexpr double delta = 1e-4;
  double diff = rho - L;
  double sq_dist = std::sqrt(diff * diff + delta);
  double rho_eff = 0.5 * (rho + L + sq_dist);
  double d_rho_eff_d_rho = 0.5 * (1.0 + diff / sq_dist);

  // Pressure regularization
  constexpr double eps2 = 1.0; // Pa^2 (regularisation scale eps = 1 Pa)
  double dP2 = dP * dP;
  double U = dP2 + eps2;
  double U_pow_025 = std::sqrt(std::sqrt(U));
  double U_pow_125 = U * U_pow_025;

  // Mass Flow Calculation
  // mdot = (Cd * E * A * sqrt(2 * rho_eff)) * (dP / U^0.25)
  double geom_const = Cd * E * area;
  double sqrt_2rho = std::sqrt(2.0 * rho_eff);
  double mdot = geom_const * sqrt_2rho * (dP / U_pow_025);

  // Jacobian w.r.t dP
  // d(mdot)/d(dP) = constant * (0.5 * dP^2 + eps^2) / (dP^2 + eps^2)^1.25
  double d_mdot_dp = geom_const * sqrt_2rho * (0.5 * dP2 + eps2) / U_pow_125;

  // Jacobian w.r.t rho
  // d(mdot)/d(rho) = (mdot / (2 * rho_eff)) * d(rho_eff)/d(rho)
  double d_mdot_rho = (mdot / (2.0 * rho_eff)) * d_rho_eff_d_rho;

  return {mdot, d_mdot_dp, d_mdot_rho};
}

std::tuple<double, double> pressure_loss_and_jacobian(double v, double rho,
                                                       double K) {
  // Standard K-factor loss: dP = K * 0.5 * rho * v^2
  double constant = K * 0.5 * rho;
  double dP = constant * v * v;

  // Analytical derivative via chain rule:
  // d(dP)/d(v) = K * rho * v
  double jacobian = K * rho * v;

  return {dP, jacobian};
}

std::tuple<double, double> lossless_pressure_and_jacobian(double P_in,
                                                          double P_out) {
  // A lossless connection preserves total pressure natively.
  // Residual: P_in - P_out = 0.
  // Analytical derivative w.r.t P_out: -1.0
  return {P_in - P_out, -1.0};
}

// -----------------------------------------------------------------------------
// 2. Heat Transfer Components
// -----------------------------------------------------------------------------

CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_dittus_boelter(double Re, double Pr, bool heating) {
  if (Re <= 0.0) {
    return {{0.0, 0.0},
            CorrelationValidity::INVALID,
            "nusselt_dittus_boelter requires Re > 0"};
  }

  double n = heating ? dittus_boelter::coeff_prandtl_heating_exp
                     : dittus_boelter::coeff_prandtl_cooling_exp;

  double constant = dittus_boelter::coeff_outer * std::pow(Pr, n);
  double Nu = constant * std::pow(Re, dittus_boelter::coeff_reynolds_exp);

  // dNu_dRe = coeff_outer * coeff_reynolds_exp * Re^(coeff_reynolds_exp - 1) *
  // Pr^n
  //         = coeff_reynolds_exp * Nu / Re
  double dNu_dRe = dittus_boelter::coeff_reynolds_exp * Nu / Re;

  return {{Nu, dNu_dRe}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_gnielinski(double Re, double Pr, double f) {
  // Exact derivative at fixed f (#446) -- previously a central difference.
  // Below the laminar anchor, including Re <= 0, Gnielinski is the laminar
  // value with zero slope: defined, so returned, but flagged as extrapolated
  // rather than INVALID.
  const auto r = nusselt_gnielinski_with_derivative(Re, Pr, f, 0.0);
  const auto validity = Re >= gnielinski::blend_re_turbulent
                            ? CorrelationValidity::VALID
                            : CorrelationValidity::EXTRAPOLATED;
  return {{r.Nu, r.dNu_dRe}, validity, ""};
}

CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_sieder_tate(double Re, double Pr, double mu_ratio) {
  if (Re <= 0.0) {
    return {{0.0, 0.0},
            CorrelationValidity::INVALID,
            "nusselt_sieder_tate requires Re > 0"};
  }
  double Nu = nusselt_sieder_tate(Re, Pr, mu_ratio);

  // Nu = 0.027 * Re^0.8 * Pr^(1/3) * mu_ratio^0.14
  // dNu/dRe = 0.8 * Nu / Re
  double dNu_dRe = sieder_tate::coeff_reynolds_exp * Nu / Re;
  return {{Nu, dNu_dRe}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_petukhov(double Re, double Pr, double f) {
  if (Re <= 0.0) {
    return {{0.0, 0.0},
            CorrelationValidity::INVALID,
            "nusselt_petukhov requires Re > 0"};
  }
  double Nu = nusselt_petukhov(Re, Pr, f);

  // Central finite difference in Re (holding f fixed), consistent with
  // Gnielinski Jacobian treatment.
  const double eps = std::max(1.0, std::abs(Re) * 1e-6);
  double Nu_plus = nusselt_petukhov(Re + eps, Pr, f);
  double Nu_minus = nusselt_petukhov(Re - eps, Pr, f);
  double dNu_dRe = (Nu_plus - Nu_minus) / (2.0 * eps);
  return {{Nu, dNu_dRe}, CorrelationValidity::VALID, ""};
}

// (f, df/dRe) of the laminar-turbulent blend (friction.h, laminar_friction).
// The laminar branch is Hagen-Poiseuille with its own negligible floor, so
// f Re = 64 and the drop is linear in m down to essentially zero flow.
static std::tuple<double, double> blend_laminar(double Re, double f_turb, double df_turb_dRe) {
  using namespace laminar_friction;
  const double Re_l = std::sqrt(Re * Re + re_floor * re_floor);
  const double dRe_l_dRe = Re / Re_l;
  const double f_lam = f_re / Re_l;
  const double df_lam_dRe = (-f_re / (Re_l * Re_l)) * dRe_l_dRe;
  const double t = std::tanh((Re_l - blend_center) / blend_width);
  const double w = 0.5 * (1.0 + t);
  const double dw_dRe = 0.5 * (1.0 - t * t) / blend_width * dRe_l_dRe;
  const double f = (1.0 - w) * f_lam + w * f_turb;
  const double df = (1.0 - w) * df_lam_dRe + w * df_turb_dRe + dw_dRe * (f_turb - f_lam);
  return {f, df};
}

// (F_Re, dF_Re/dRe) where F_Re(Re) = (1 - w(Re)) * f_re + w(Re) * f_turb(Re) * Re:
// the re_floor = 0 limit of Re * blend_laminar, smooth and cancellation-free
// for all Re >= 0 (#1216).
static std::tuple<double, double> blend_laminar_f_re(double Re, double f_turb,
                                                     double df_turb_dRe) {
  using namespace laminar_friction;
  const double t = std::tanh((Re - blend_center) / blend_width);
  const double w = 0.5 * (1.0 + t);
  const double dw_dRe = 0.5 * (1.0 - t * t) / blend_width;
  const double F_Re = (1.0 - w) * f_re + w * f_turb * Re;
  const double dF_dRe = w * (f_turb + Re * df_turb_dRe) + dw_dRe * (f_turb * Re - f_re);
  return {F_Re, dF_dRe};
}

static std::tuple<double, double> f_turb_and_dRe_haaland(double Re, double e_D) {
  const double re_eps = laminar_friction::turbulent_re_floor;
  const double Re_eff = std::sqrt(Re * Re + re_eps * re_eps);
  const double a = std::pow(e_D / haaland::coeff_roughness, haaland::coeff_exponent);
  const double b = haaland::coeff_reynolds / Re_eff;
  const double arg = a + b;
  const double inv_sqrt_f = haaland::coeff_outer * std::log10(arg);
  const double dU_darg = haaland::coeff_outer / (M_LN10 * arg);
  const double darg_dRe_eff = -haaland::coeff_reynolds / (Re_eff * Re_eff);
  const double dRe_eff_dRe = Re / Re_eff;
  const double darg_dRe = darg_dRe_eff * dRe_eff_dRe;
  const double dU_dRe = dU_darg * darg_dRe;
  const double f_turb = 1.0 / (inv_sqrt_f * inv_sqrt_f);
  const double df_turb_dU = -2.0 / (inv_sqrt_f * inv_sqrt_f * inv_sqrt_f);
  const double df_turb_dRe = df_turb_dU * dU_dRe;
  return {f_turb, df_turb_dRe};
}

static std::tuple<double, double> f_turb_and_dRe_serghides(double Re, double e_D) {
  const double re_eps = laminar_friction::turbulent_re_floor;
  const double Re_eff = std::sqrt(Re * Re + re_eps * re_eps);
  const double r = serghides::coeff_roughness;
  const double cA = serghides::coeff_reynolds_A;
  const double cBC = serghides::coeff_reynolds_BC;
  const double ln10 = M_LN10;

  const double argA = e_D / r + cA / Re_eff;
  const double A = -2.0 * std::log10(argA);
  const double argB = e_D / r + cBC * A / Re_eff;
  const double B = -2.0 * std::log10(argB);
  const double argC = e_D / r + cBC * B / Re_eff;
  const double C = -2.0 * std::log10(argC);

  const double denom = C - 2.0 * B + A;
  const double inv_sqrt_f =
      (std::abs(denom) < 1e-30) ? C : A - (B - A) * (B - A) / denom;

  const double dRe_eff_dRe = Re / Re_eff;
  const double dA_dRe = (2.0 / (ln10 * argA)) * (cA / (Re_eff * Re_eff)) * dRe_eff_dRe;
  const double dargB_dRe = cBC * (dA_dRe * Re_eff - A * dRe_eff_dRe) / (Re_eff * Re_eff);
  const double dB_dRe = (-2.0 / (ln10 * argB)) * dargB_dRe;
  const double dargC_dRe = cBC * (dB_dRe * Re_eff - B * dRe_eff_dRe) / (Re_eff * Re_eff);
  const double dC_dRe = (-2.0 / (ln10 * argC)) * dargC_dRe;

  double d_inv_sqrt_f_dRe;
  if (std::abs(denom) < 1e-30) {
    d_inv_sqrt_f_dRe = dC_dRe;
  } else {
    const double BmA = B - A;
    const double d_BmA_dRe = dB_dRe - dA_dRe;
    const double d_denom_dRe = dC_dRe - 2.0 * dB_dRe + dA_dRe;
    d_inv_sqrt_f_dRe =
        dA_dRe - (2.0 * BmA * d_BmA_dRe * denom - BmA * BmA * d_denom_dRe) /
                     (denom * denom);
  }
  const double f_turb = 1.0 / (inv_sqrt_f * inv_sqrt_f);
  const double df_turb_dRe = (-2.0 / (inv_sqrt_f * inv_sqrt_f * inv_sqrt_f)) * d_inv_sqrt_f_dRe;
  return {f_turb, df_turb_dRe};
}

static std::tuple<double, double> f_turb_and_dRe_colebrook(double Re, double e_D) {
  const double re_eps = laminar_friction::turbulent_re_floor;
  const double Re_eff = std::sqrt(Re * Re + re_eps * re_eps);
  const double f_val = ::friction_colebrook(Re, e_D);
  const double x = 1.0 / std::sqrt(f_val);
  const double ln10 = M_LN10;
  const double arg =
      e_D / serghides::coeff_roughness + serghides::coeff_reynolds_BC * x / Re_eff;
  const double dF_dx = 1.0 + 2.0 * serghides::coeff_reynolds_BC / (ln10 * arg * Re_eff);
  const double dF_dRe_eff =
      (-2.0 / (ln10 * arg)) * (serghides::coeff_reynolds_BC * x / (Re_eff * Re_eff));
  const double dx_dRe = -(dF_dRe_eff * (Re / Re_eff)) / dF_dx;
  const double df_turb_dRe = (-2.0 / (x * x * x)) * dx_dRe;
  return {f_val, df_turb_dRe};
}

static std::tuple<double, double> f_turb_and_dRe_petukhov(double Re) {
  const double re_eps = laminar_friction::turbulent_re_floor;
  const double Re_eff = std::sqrt(Re * Re + re_eps * re_eps);
  const double x = petukhov::coeff_a * std::log(Re_eff) - petukhov::coeff_b;
  const double f_turb = 1.0 / (x * x);
  const double dx_dRe = petukhov::coeff_a * Re / (Re_eff * Re_eff);
  const double df_turb_dRe = (-2.0 / (x * x * x)) * dx_dRe;
  return {f_turb, df_turb_dRe};
}

CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_haaland(double Re, double e_D) {
  const auto [f_turb, df_turb_dRe] = f_turb_and_dRe_haaland(Re, e_D);
  const auto [f, jacobian] = blend_laminar(Re, f_turb, df_turb_dRe);
  return {{f, jacobian}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_serghides(double Re, double e_D) {
  const auto [f_turb, df_turb_dRe] = f_turb_and_dRe_serghides(Re, e_D);
  const auto [f, jacobian] = blend_laminar(Re, f_turb, df_turb_dRe);
  return {{f, jacobian}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_colebrook(double Re, double e_D) {
  const auto [f_turb, df_turb_dRe] = f_turb_and_dRe_colebrook(Re, e_D);
  const auto [f, jacobian] = blend_laminar(Re, f_turb, df_turb_dRe);
  return {{f, jacobian}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_petukhov(double Re) {
  const auto [f_turb, df_turb_dRe] = f_turb_and_dRe_petukhov(Re);
  const auto [f, jacobian] = blend_laminar(Re, f_turb, df_turb_dRe);
  return {{f, jacobian}, CorrelationValidity::VALID, ""};
}

CorrelationResult<std::tuple<double, double>>
friction_and_jacobian(const std::string &tag, double Re, double e_D) {
  // An unknown tag throws (#1213), as Registry::get_friction_model and
  // friction_f_re_at_zero do. Every consumer reads .result, so the INVALID
  // {0, 0} this returned made a typo a frictionless duct that solved.
  return Registry::get().get_friction_model(tag)(Re, e_D);
}

double laminar_blend_f_re_at_zero() {
  using namespace laminar_friction;
  const double w0 = 0.5 * (1.0 + std::tanh(-blend_center / blend_width));
  return f_re * (1.0 - w0);
}

double friction_f_re_at_zero(const std::string &tag) {
  return Registry::get().friction_f_re_at_zero(tag);
}

std::tuple<double, double> friction_f_re_and_jacobian_haaland(double Re, double e_D) {
  const auto [f_turb, df_turb_dRe] = f_turb_and_dRe_haaland(Re, e_D);
  return blend_laminar_f_re(Re, f_turb, df_turb_dRe);
}

std::tuple<double, double> friction_f_re_and_jacobian_serghides(double Re, double e_D) {
  const auto [f_turb, df_turb_dRe] = f_turb_and_dRe_serghides(Re, e_D);
  return blend_laminar_f_re(Re, f_turb, df_turb_dRe);
}

std::tuple<double, double> friction_f_re_and_jacobian_colebrook(double Re, double e_D) {
  const auto [f_turb, df_turb_dRe] = f_turb_and_dRe_colebrook(Re, e_D);
  return blend_laminar_f_re(Re, f_turb, df_turb_dRe);
}

std::tuple<double, double> friction_f_re_and_jacobian_petukhov(double Re) {
  const auto [f_turb, df_turb_dRe] = f_turb_and_dRe_petukhov(Re);
  return blend_laminar_f_re(Re, f_turb, df_turb_dRe);
}

std::tuple<double, double> friction_f_re_and_jacobian(const std::string &tag,
                                                      double Re, double e_D) {
  if (Registry::get().has_friction_model(tag)) {
    return Registry::get().get_friction_f_re_model(tag)(Re, e_D);
  }
  return {0.0, 0.0};
}

// -----------------------------------------------------------------------------
// 4. Thermodynamic & Transport Components
// -----------------------------------------------------------------------------

std::tuple<double, double, double>
density_and_jacobians(double T, double P, const std::vector<double> &X) {
  // Robust clamping with leaky gradient to avoid solver divergence into negative space
  double T_eff = (T > 10.0) ? T : (10.0 + 0.1 * (T - 10.0));
  double P_eff = (P > 1.0) ? P : (1.0 + 0.1 * (P - 1.0));
  double dTeff_dT = (T > 10.0) ? 1.0 : 0.1;
  double dPeff_dP = (P > 1.0) ? 1.0 : 0.1;

  // Calculate base density natively
  double rho = combaero::density(T_eff, P_eff, X);

  // Ideal gas law: rho = P * W / (R * T)
  // d_rho_dT = d(rho)/d(T_eff) * d(T_eff)/dT
  double d_rho_d_T = (-rho / T_eff) * dTeff_dT;

  // d_rho_dP = d(rho)/d(P_eff) * d(P_eff)/dP
  double d_rho_d_P = (rho / P_eff) * dPeff_dP;

  return {rho, d_rho_d_T, d_rho_d_P};
}

std::vector<double> viscosity_mass_gradient(double T, double P,
                                            const std::vector<double> &X) {
  const double T_eff = (T > 10.0) ? T : (10.0 + 0.1 * (T - 10.0));
  return combaero::mole_gradient_to_mass(X, combaero::viscosity_dX(T_eff, P, X));
}

std::tuple<double, double> enthalpy_and_jacobian(double T,
                                                  const std::vector<double> &X) {
  double T_eff = (T > 10.0) ? T : (10.0 + 0.1 * (T - 10.0));
  double dTeff_dT = (T > 10.0) ? 1.0 : 0.1;

  // h(T) = base mass enthalpy
  double h_val = combaero::h_mass(T_eff, X);

  // Analytical derivative of enthalpy with respect to temperature is Cp
  double d_h_d_T = combaero::cp_mass(T_eff, X) * dTeff_dT;

  return {h_val, d_h_d_T};
}

std::tuple<double, double, double>
viscosity_and_jacobians(double T, double P, const std::vector<double> &X) {

  double T_eff = (T > 10.0) ? T : (10.0 + 0.1 * (T - 10.0));
  double P_eff = (P > 1.0) ? P : (1.0 + 0.1 * (P - 1.0));
  double dTeff_dT = (T > 10.0) ? 1.0 : 0.1;
  double dPeff_dP = (P > 1.0) ? 1.0 : 0.1;

  // Analytic d/dT through the C1 collision-integral interpolation and the
  // Wilke rule (#485); it was a +/-1 mK central difference, exact on each
  // linear segment of the old table but discontinuous at every node.
  // Dilute-gas viscosity does not depend on pressure.
  (void)dPeff_dP;
  const combaero::TransportDerivatives tr = combaero::transport_and_dT(T_eff, P_eff, X);
  return {tr.mu, tr.dmu_dT * dTeff_dT, 0.0};
}

std::tuple<double, double>
mach_number_and_jacobian_v(double v, double T, const std::vector<double> &X) {
  double T_eff = (T > 10.0) ? T : (10.0 + 0.1 * (T - 10.0));

  double a = combaero::speed_of_sound(T_eff, X);
  if (a <= 1.0) a = 1.0; // Minimal fallback

  double M = v / a;
  // dM/dv = 1/a
  double dM_dv = 1.0 / a;

  return {M, dM_dv};
}

std::tuple<double, double>
T_adiabatic_wall_and_jacobian_v(double T_static, double v, double T, double P,
                                 const std::vector<double> &X, bool turbulent) {
  double T_s_eff = (T_static > 10.0) ? T_static : (10.0 + 0.1 * (T_static - 10.0));
  double v_eff = (v >= 0.0) ? v : (0.1 * v);
  double T_eff = (T > 10.0) ? T : (10.0 + 0.1 * (T - 10.0));
  double P_eff = (P > 1.0) ? P : (1.0 + 0.1 * (P - 1.0));

  double Pr = combaero::prandtl(T_eff, P_eff, X);
  double r = turbulent ? std::cbrt(Pr) : std::sqrt(Pr);

  double mw_g = combaero::mwmix(X);
  double cp_val = combaero::cp(T_eff, X) / mw_g * 1000.0;

  double T_aw = T_s_eff + r * 0.5 * v_eff * v_eff / cp_val;
  // Use v_eff for derivative to keep it leaky
  double dT_aw_dv = r * v_eff / cp_val;

  return {T_aw, dT_aw_dv};
}

std::tuple<double, double>
Tt_from_static_and_jacobian_M(double T, double M,
                               const std::vector<double> &X) {
  // Enforce non-negative M (Tt_from_static requires M >= 0).
  // Use a smooth floor so the derivative remains informative near M = 0.
  // For solver purposes the useful range is M >= 0; negative values are
  // mapped to zero with a leaky 0.01 gradient.
  double M_eff = (M >= 0.0) ? M : 0.01 * M;
  double dMeff_dM = (M >= 0.0) ? 1.0 : 0.01;

  double T0 = Tt_from_static(T, M_eff, X);

  // Analytical derivative via the definition h(T0, X) = h(T, X) + 0.5 * v^2
  // where v = M_eff * a(T, X).  Differentiating w.r.t. M_eff:
  //   cp_mass(T0) * dT0/dM_eff = M_eff * a(T)^2
  // => dT0/dM_eff = M_eff * a^2 / cp_mass(T0)
  double a = combaero::speed_of_sound(T, X);
  double cp_T0 = combaero::cp_mass(T0, X);
  double dT0_dMeff = (cp_T0 > 1.0) ? (M_eff * a * a / cp_T0) : 0.0;

  double dT0_dM = dT0_dMeff * dMeff_dM;

  return {T0, dT0_dM};
}

std::tuple<double, double>
Pt_from_static_and_jacobian_M(double P, double T, double M,
                               const std::vector<double> &X,
                               double tol, std::size_t max_iter) {
  // Enforce non-negative M with a leaky floor for differentiability.
  double M_eff = (M >= 0.0) ? M : 0.01 * M;
  double dMeff_dM = (M >= 0.0) ? 1.0 : 0.01;

  double P0 = Pt_from_static(P, T, M_eff, X, tol, max_iter);

  // Chain rule: P0 = P * exp((s(T0) - s(T,P)) / R_specific)
  // where s(T0) = s(T0, X, P_REF) and T0 = Tt_from_static(T, M_eff, X).
  //
  // dP0/dM_eff = P0 * d/dM_eff [ln(P0 / P)]
  //            = P0 * d/dT0 [(s(T0) - const) / R_specific] * dT0/dM_eff
  //            = P0 * [cp_mass(T0) / T0 / R_specific] * (dT0/dM_eff)
  //
  // The positive sign comes from: s increases with T (at const P_ref) for ideal
  // gas, so raising T0 increases the entropy difference, increasing P0/P.
  //
  // dT0/dM_eff (analytical, from T0 definition): M_eff * a(T)^2 / cp_mass(T0)
  double a = combaero::speed_of_sound(T, X);
  double T0_val = Tt_from_static(T, M_eff, X);
  double cp_T0 = combaero::cp_mass(T0_val, X);
  double R_spec = combaero::specific_gas_constant(X); // J/(kg*K)

  double dT0_dMeff = (cp_T0 > 1.0) ? (M_eff * a * a / cp_T0) : 0.0;

  // dP0/dT0 = P0 * (+cp_mass(T0) / T0) / R_specific
  double dP0_dT0 = (T0_val > 0.0 && R_spec > 1.0)
                       ? P0 * (cp_T0 / T0_val) / R_spec
                       : 0.0;

  double dP0_dM = dP0_dT0 * dT0_dMeff * dMeff_dM;

  return {P0, dP0_dM};
}

// -----------------------------------------------------------------------------
// 6. Combustion Interfaces
// -----------------------------------------------------------------------------

namespace {

// Complete combustion at constant enthalpy with analytic derivatives (#716):
// h(T_ad, X_b) = h(T_in, X_in) MW_b / MW_in (molar h, mass conserved), so
//   dT_ad/dT_in = cp(T_in, X_in) (MW_b / MW_in) / cp(T_ad, X_b),
//   dT_ad/dX_k  = -[sum_j h_j(T_ad) dXb_jk - d(h_in MW_b / MW_in)/dX_k] / cp_b,
// with dXb/dX from complete_combustion_products_and_jacobian. X is taken
// unnormalised (chain with a dX keeping sum X = 1).
struct CompleteAdiabatic {
  double T_ad = 0.0;
  std::vector<double> X_b;
  double dT_dTin = 0.0;
  std::vector<double> dT_dX;
  std::vector<std::vector<double>> dXb_dX;
};

CompleteAdiabatic complete_adiabatic(double T_in, double P, const std::vector<double> &X_in) {
  // Clamp to a safe physical floor instead of throwing: the solver may probe
  // unphysical intermediates; below the floor T_in has no effect.
  const bool clamped = T_in < thermo::T_FLOOR;
  T_in = std::max(T_in, thermo::T_FLOOR);
  P = std::max(P, kEquilibriumPressureFloor);
  const std::size_t ns = X_in.size();
  CompleteAdiabatic out;
  const CompleteCombustionJacobian pj = combaero::complete_combustion_products_and_jacobian(X_in);
  out.X_b = pj.X_b;
  out.dXb_dX = pj.dXb_dX;
  State in_state;
  in_state.set_TPX(T_in, P, X_in);
  out.T_ad = combaero::complete_combustion(in_state, true, true).T;

  const double mw_in = combaero::mwmix(X_in);
  const double mw_b = combaero::mwmix(out.X_b);
  const double h_in = combaero::h(T_in, X_in);
  const double cp_b = combaero::cp(out.T_ad, out.X_b);
  out.dT_dTin = clamped ? 0.0 : combaero::cp(T_in, X_in) * (mw_b / mw_in) / cp_b;
  std::vector<double> h_ad(ns), W(ns);
  for (std::size_t j = 0; j < ns; ++j) {
    h_ad[j] = combaero::h_species(j, out.T_ad);
    W[j] = combaero::species_molar_mass(j);
  }
  out.dT_dX.assign(ns, 0.0);
  for (std::size_t k = 0; k < ns; ++k) {
    double dhb = 0.0, dmwb = 0.0;
    for (std::size_t j = 0; j < ns; ++j) {
      dhb += h_ad[j] * out.dXb_dX[j][k];
      dmwb += W[j] * out.dXb_dX[j][k];
    }
    const double dtarget = combaero::h_species(k, T_in) * mw_b / mw_in + h_in * dmwb / mw_in -
                           h_in * mw_b * W[k] / (mw_in * mw_in);
    out.dT_dX[k] = -(dhb - dtarget) / cp_b;
  }
  return out;
}

}  // namespace

std::tuple<double, double, std::vector<double>>
adiabatic_T_complete_and_jacobian_T(double T_in, double P,
                                     const std::vector<double> &X_in) {
  const CompleteAdiabatic c = complete_adiabatic(T_in, P, X_in);
  return {c.T_ad, c.dT_dTin, c.X_b};
}

std::tuple<double, double, double, std::vector<double>>
adiabatic_T_equilibrium_and_jacobians(double T_in, double P,
                                       const std::vector<double> &X_in) {
  const bool clamped = T_in < thermo::T_FLOOR;
  T_in = std::max(T_in, thermo::T_FLOOR);
  P = std::max(P, kEquilibriumPressureFloor);
  State in_state;
  in_state.set_TPX(T_in, P, X_in);
  // Analytic (#716): implicit in the equilibrium extents and the adiabatic root.
  const EquilibriumJacobian ej = combaero::combustion_equilibrium_and_jacobian(in_state);
  return {ej.state.T, clamped ? 0.0 : ej.dT_dT_in, ej.dT_dP, ej.state.X};
}

// -----------------------------------------------------------------------------
// 7. Stream-Based Network Solvers
// -----------------------------------------------------------------------------

EntropyGenerationResult entropy_generation(const std::vector<Stream> &inflows,
                                           const std::vector<Stream> &outflows,
                                           double S_heat) {
  if (inflows.empty()) {
    throw std::invalid_argument("entropy_generation: a control volume needs an inflow");
  }
  auto flux = [](const std::vector<Stream> &streams, double &m_sum) {
    double S = 0.0;
    for (const auto &st : streams) {
      if (!(st.m_dot >= 0.0)) {
        throw std::invalid_argument("entropy_generation: a stream's m_dot must be >= 0");
      }
      m_sum += st.m_dot;
      S += st.m_dot * s_mass(st.Tt, mass_to_mole(st.Y), st.Pt);
    }
    return S;
  };
  double m_in = 0.0;
  double m_out = 0.0;
  const double S_in = flux(inflows, m_in);
  const double S_out = flux(outflows, m_out);
  const Stream &ref = *std::max_element(
      inflows.begin(), inflows.end(),
      [](const Stream &a, const Stream &b) { return a.m_dot < b.m_dot; });
  const double sigma = S_out - S_in - S_heat;
  const double scale = specific_gas_constant(mass_to_mole(ref.Y)) * m_in;
  return {sigma, scale > 0.0 ? sigma / scale : 0.0};
}

double heat_entropy_flux(double m_dot, double Tt_in, double Tt_out,
                         const std::vector<double> &Y, double M_in, double M_out) {
  if (!(Tt_in > 0.0 && Tt_out > 0.0)) {
    throw std::invalid_argument("heat_entropy_flux: temperatures must be positive");
  }
  if (!(M_in >= 0.0 && M_out >= 0.0)) {
    throw std::invalid_argument("heat_entropy_flux: Mach numbers must be >= 0");
  }
  if (Tt_in == Tt_out) {
    return 0.0;
  }
  const std::vector<double> X = mass_to_mole(Y);
  // integral cp / T = integral cp / T0 + integral cp (1/T - 1/T0). The first
  // is exactly the constant-pressure entropy rise (the NASA polynomials' own
  // s, no quadrature across their 1000 K break); only the small Mach part is
  // integrated.
  const double at_rest = s_mass(Tt_out, X, 1.0e5) - s_mass(Tt_in, X, 1.0e5);
  if (M_in == 0.0 && M_out == 0.0) {
    return m_dot * at_rest;
  }
  const double span = Tt_out - Tt_in;
  auto mach_part = [&](double T0) {
    const double M = M_in + (M_out - M_in) * (T0 - Tt_in) / span;
    const double T = M > 0.0 ? T_from_stagnation(T0, M, X) : T0;
    return cp_mass(T0, X) * (1.0 / T - 1.0 / T0);
  };
  return m_dot * (at_rest + composite_simpson(mach_part, Tt_in, Tt_out, 64));
}

MixerResult
mixer_from_streams_and_jacobians(const std::vector<Stream> &streams,
                                double Q, double fraction) {
  std::size_t n_streams = streams.size();
  std::size_t n_species = combaero::num_species();

  // The streams that make up the mixed state (#1232): those arriving (m_i >
  // 0), or every stream when none arrives. A stream leaving the node carries
  // the node's gas away; its own state is not part of what the node holds.
  // The signed mean of all streams, sum m_i q_i / sum m_i, left the range of
  // its inputs with streams of both signs and had a pole at zero net flow
  // (Tt -1700 K at -0.02 kg/s net, 1e10 K near zero; a negative mean Y
  // clamped). Over streams of one sign the mean is a convex combination, so
  // Y stays in the simplex and h and Pt within the inputs for any nonzero
  // total: there is no pole, and no band near zero flow where another state
  // stands in (the band keeps only zero partials, below). With no stream arriving (all m_i <= 0) it is the mean of the
  // leaving streams (#1224), as before; the state then jumps from the last
  // arriving stream's own state to that mean as that stream crosses zero.
  const std::vector<bool> in_mix = streams_in_mix(streams);

  // mdot_tot: the total flow of the streams in the mix (all of them when
  // they share a sign). Zero only when no stream flows.
  double mdot_tot = 0.0;
  std::vector<double> Y_mix(n_species, 0.0);
  double H_tot = 0.0;
  double P_total_tot = 0.0;

  std::vector<double> h_stream(n_streams, 0.0);
  std::vector<double> cp_stream(n_streams, 0.0);
  std::vector<std::vector<double>> hk_mass(n_streams,
                                            std::vector<double>(n_species, 0.0));

  for (std::size_t i = 0; i < n_streams; ++i) {
    for (std::size_t k = 0; k < n_species; ++k) {
      hk_mass[i][k] = J_per_mol_to_J_per_kg(combaero::h_species(k, streams[i].Tt),
                                             combaero::species_molar_mass(k));
      h_stream[i] += streams[i].Y[k] * hk_mass[i][k];

      double cpk_mass = J_per_mol_to_J_per_kg(combaero::cp_species(k, streams[i].Tt),
                                               combaero::species_molar_mass(k));
      cp_stream[i] += streams[i].Y[k] * cpk_mass;
    }
    if (!in_mix[i]) {
      continue;
    }
    mdot_tot += streams[i].m_dot;
    for (std::size_t k = 0; k < n_species; ++k) {
      Y_mix[k] += streams[i].m_dot * streams[i].Y[k];
    }
    H_tot += streams[i].m_dot * h_stream[i];
    P_total_tot += streams[i].m_dot * streams[i].Pt;
  }

  // Y, h and Pt are the mean sum_i m_i q_i / m over the streams in the mix
  // at any nonzero m; with no flow (every m_i = 0) stream 0's state stands
  // in. The partials below are that mean's outside the near-zero band |m| <=
  // kMixerNoFlowBand, and 0 inside it: there the mean's own partials, of
  // order (q_i - q) / m, made the first Newton step from a near-zero start go
  // the wrong way (cond 1e18 against 6e6; 7 of 7 such starts stalled that
  // converge with the band, #1238 review). The state is the same mean inside
  // and outside the band, so it is continuous at the band's edge; only its
  // slope is dropped. The heat per unit flow takes |m|. Reversing every
  // stream of one sign gives the same state bit for bit (#1224).
  const bool averaged = mdot_tot != 0.0;
  // A single stream has no difference q_i - q to divide by m: its mixed
  // state is its own (plus the heat, which fades C1 to none at zero flow),
  // so its partials are bounded down to zero flow and kept in the band. A
  // wall-heated duct's own exit state is such a mix; zeroed, it dropped the
  // exit Tt's 1:1 move with the duct's feed at exactly zero flow (#1265).
  const bool one_stream = n_streams == 1;
  const bool with_partials = std::abs(mdot_tot) > kMixerNoFlowBand || one_stream;
  // The weight m_i / m and 1 / m of the partials below. With no flow at all
  // (one stream, m = 0) the stream's own state stands in, weight 1, and
  // every difference q_i - q_ref is exactly 0, so 1 / m is never needed.
  const double inv_m = averaged ? 1.0 / mdot_tot : 0.0;
  if (averaged) {
    for (std::size_t k = 0; k < n_species; ++k) {
      Y_mix[k] /= mdot_tot;
    }
  } else if (n_streams > 0) {
    Y_mix = streams[0].Y;
  }

  // Compute adiabatic mass-averaged enthalpy
  double h_mix_base = averaged ? (H_tot / mdot_tot)
                               : (n_streams > 0 ? h_stream[0] : 0.0);

  // Heat per unit flow: Q / |m|, fading out only near zero flow.
  const auto [phi_Q, dphi_Q] = mixer_heat_per_unit_flow(mdot_tot);

  std::vector<double> normalized_Y_mix = combaero::normalize_fractions(Y_mix);

  // 'fraction' scales the SENSIBLE enthalpy, h - h(T_ref) at the mixture's
  // composition: absolute h includes formation enthalpy, so a fraction of it
  // heated air at ambient and combustion products alike (#481).
  std::vector<double> href_k(n_species, 0.0);
  double h_ref = 0.0;
  for (std::size_t k = 0; k < n_species; ++k) {
    href_k[k] = J_per_mol_to_J_per_kg(combaero::h_species(k, kSensibleEnthalpyRefT),
                                      combaero::species_molar_mass(k));
    h_ref += normalized_Y_mix[k] * href_k[k];
  }

  // Apply energy transfer: delta_h = Q phi + fraction * (h_mix_base - h_ref)
  double delta_h = Q * phi_Q + fraction * (h_mix_base - h_ref);

  double h_mix = h_mix_base + delta_h;
  double Pt_mix = averaged ? (P_total_tot / mdot_tot)
                           : (n_streams > 0 ? streams[0].Pt : 0.0);

  std::vector<double> X_mix = combaero::mass_to_mole(normalized_Y_mix);
  double T_guess = 300.0;
  if (n_streams > 0)
    T_guess = streams[0].Tt;
  if (n_streams > 0 && !in_mix[0]) {
    for (std::size_t i = 0; i < n_streams; ++i) {
      if (in_mix[i] && streams[i].m_dot != 0.0) {
        T_guess = streams[i].Tt;
        break;
      }
    }
  }
  const double Tt_mix = T_from_h_mass_continued(h_mix, X_mix, T_guess);
  const double T_eval = Tt_mix;

  MixerResult res;
  res.Tt_mix = Tt_mix;
  res.Pt_mix = Pt_mix;
  res.Y_mix = normalized_Y_mix;
  res.dTt_mix_d_stream.resize(n_streams);
  res.dPt_mix_d_stream.resize(n_streams);
  res.dY_mix_d_stream.assign(n_species, std::vector<StreamJacobian>(n_streams));

  double cp_mix = 0.0;
  std::vector<double> hk_mix(n_species, 0.0);
  for (std::size_t k = 0; k < n_species; ++k) {
    hk_mix[k] =
        J_per_mol_to_J_per_kg(combaero::h_species(k, T_eval), combaero::species_molar_mass(k));
    double cpk_mix =
        J_per_mol_to_J_per_kg(combaero::cp_species(k, T_eval), combaero::species_molar_mass(k));
    cp_mix += Y_mix[k] * cpk_mix;
  }

  // Store dT/d(delta_h) for reference (not directly used in Jacobians)
  res.dTt_mix_d_delta_h = (cp_mix > 0.0) ? (1.0 / cp_mix) : 0.0;
  res.dTt_mix_dQ = res.dTt_mix_d_delta_h * phi_Q;

  // The flow partials are differences against the mixed state, q_i - q_mix
  // with q_mix = sum_j m_j q_j / m (#1198). Where every flowing stream in the
  // mix (m_j != 0) carries the same q -- one stream, or streams of one gas -- q_mix is
  // that q exactly, but m q / m rounds and left ~1e-16 partials where the
  // exact one is 0. The difference is taken against the common input value
  // there, so it cancels exactly; elsewhere it is against q_mix as before.
  // The mixed state itself (Y_mix, h_mix_base, Pt_mix) is not changed.
  std::vector<double> Y_diff_ref = Y_mix;
  double h_diff_ref = h_mix_base;
  double Pt_diff_ref = Pt_mix;
  {
    std::size_t first = n_streams;
    for (std::size_t i = 0; i < n_streams; ++i) {
      if (streams[i].m_dot != 0.0 && in_mix[i]) {
        first = i;
        break;
      }
    }
    if (first < n_streams) {
      bool same_h = averaged;
      bool same_Pt = averaged;
      std::vector<bool> same_Y(n_species, true);
      for (std::size_t i = first + 1; i < n_streams; ++i) {
        if (streams[i].m_dot == 0.0 || !in_mix[i]) {
          continue;
        }
        same_h = same_h && h_stream[i] == h_stream[first];
        same_Pt = same_Pt && streams[i].Pt == streams[first].Pt;
        for (std::size_t k = 0; k < n_species; ++k) {
          same_Y[k] = same_Y[k] && streams[i].Y[k] == streams[first].Y[k];
        }
      }
      for (std::size_t k = 0; k < n_species; ++k) {
        if (same_Y[k]) {
          Y_diff_ref[k] = streams[first].Y[k];
        }
      }
      if (same_h) {
        h_diff_ref = h_stream[first];
      }
      if (same_Pt) {
        Pt_diff_ref = streams[first].Pt;
      }
    }
  }

  for (std::size_t i = 0; i < n_streams; ++i) {
    StreamJacobian &dT_jac = res.dTt_mix_d_stream[i];
    StreamJacobian &dP_jac = res.dPt_mix_d_stream[i];
    dT_jac.d_Y.assign(n_species, 0.0);
    dP_jac.d_Y.assign(n_species, 0.0);

    // A stream outside the mix (leaving while another arrives) moves nothing.
    if (with_partials && in_mix[i]) {
      const double w = averaged ? streams[i].m_dot * inv_m : 1.0;
      // P_total sensitivities
      dT_jac.d_Pt = 0.0;
      dP_jac.d_Pt = w;
      dP_jac.d_mdot = (streams[i].Pt - Pt_diff_ref) * inv_m;
      dP_jac.d_Tt = 0.0;

      // d(Tt_mix)/d(T_i)
      // Base contribution
      double dT_base_dT = w * cp_stream[i] / cp_mix;

      // Fraction contribution: d(fraction * (h_mix_base - h_ref))/d(T_i)
      // = fraction * (mdot_i / mdot_tot) * cp_i  (h_ref does not move with T)
      double dT_frac_dT = fraction * w * cp_stream[i] / cp_mix;

      dT_jac.d_Tt = dT_base_dT + dT_frac_dT;

      // d(Tt_mix)/d(m_dot_i). h_mix = (1+f) h_base - f h_ref + Q/m_eff, then T
      // at fixed h moves with the composition: -sum_k hk(Tt_mix) dY_k/dm_i.
      // Only the enthalpy terms carry (1+f), not the composition term.
      double h_diff = (h_stream[i] - h_diff_ref);
      double y_sum = 0.0;
      double ref_sum = 0.0;
      for (std::size_t k = 0; k < n_species; ++k) {
        y_sum += hk_mix[k] * (streams[i].Y[k] - Y_diff_ref[k]);
        ref_sum += href_k[k] * (streams[i].Y[k] - Y_diff_ref[k]);
      }
      double dh_mix_dmdot = ((1.0 + fraction) * h_diff - fraction * ref_sum) * inv_m;

      // Q contribution: d(Q phi)/d(mdot_i) = Q phi'
      double dh_Q_dmdot = Q * dphi_Q;

      dT_jac.d_mdot = (dh_mix_dmdot + dh_Q_dmdot - y_sum * inv_m) / cp_mix;

      // d(Tt_mix)/d(Y_i,k): h_base moves by w hk_i, h_ref by w (href_k - h_ref)
      // and the normalised composition by w (e_k - Y), w = m_i / m; the last
      // shifts T at fixed h by -w (hk_mix - h_mix) -- the FINAL h_mix, which
      // the base-only form (h_mix_base) got wrong whenever Q or f was set.
      for (std::size_t k = 0; k < n_species; ++k) {
        dT_jac.d_Y[k] = w *
                        ((1.0 + fraction) * hk_mass[i][k] - fraction * (href_k[k] - h_ref) -
                         hk_mix[k] + h_mix) /
                        cp_mix;
      }

      // Y_mix sensitivities
      for (std::size_t k = 0; k < n_species; ++k) {
        StreamJacobian &dY_jac_k = res.dY_mix_d_stream[k][i];
        dY_jac_k.d_Y.resize(n_species);
        dY_jac_k.d_mdot = (streams[i].Y[k] - Y_diff_ref[k]) * inv_m;
        dY_jac_k.d_Tt = 0.0;
        dY_jac_k.d_Pt = 0.0;
        for (std::size_t j = 0; j < n_species; ++j) {
            double delta = (k == j) ? 1.0 : 0.0;
            dY_jac_k.d_Y[j] = w * (delta - Y_mix[k]);
        }
      }
    } else {
      dT_jac.d_mdot = 0.0;
      dT_jac.d_Tt = 0.0;
      dT_jac.d_Pt = 0.0;
      dP_jac.d_mdot = 0.0;
      dP_jac.d_Tt = 0.0;
      dP_jac.d_Pt = 0.0;
      for (std::size_t k = 0; k < n_species; ++k) {
        StreamJacobian &dY_jac_k = res.dY_mix_d_stream[k][i];
        dY_jac_k.d_mdot = 0.0;
        dY_jac_k.d_Tt = 0.0;
        dY_jac_k.d_Pt = 0.0;
        dY_jac_k.d_Y.assign(n_species, 0.0);
      }
    }
  }

  return res;
}

MixerResult adiabatic_T_complete_and_jacobian_T_from_streams(
    const std::vector<Stream> &streams, double P, double Q, double fraction) {
  // First mix adiabatically to get base state
  MixerResult mix = mixer_from_streams_and_jacobians(streams);

  // The flow the heat is spread over: that of the streams the mixer weights
  // (#1232).
  const std::vector<bool> in_mix = streams_in_mix(streams);
  double mdot_tot = 0.0;
  std::size_t n_species = combaero::num_species();
  for (std::size_t i_s = 0; i_s < streams.size(); ++i_s) {
    if (in_mix[i_s]) {
      mdot_tot += streams[i_s].m_dot;
    }
  }

  // Heat per unit flow: Q / |m|, fading out only near zero flow.
  const auto [phi_Q, dphi_Q] = mixer_heat_per_unit_flow(mdot_tot);

  // Q now; 'fraction' of the PRODUCTS' sensible enthalpy below, once T_ad
  // and the burnt composition are known (#481).
  double delta_h = Q * phi_Q;

  std::vector<double> X_mix = combaero::mass_to_mole(combaero::normalize_fractions(mix.Y_mix));
  const CompleteAdiabatic cad = complete_adiabatic(mix.Tt_mix, P, X_mix);
  double T_ad = cad.T_ad;
  double dT_ad_dT_mix = cad.dT_dTin;
  std::vector<double> X_b = cad.X_b;
  std::vector<double> Y_b = combaero::mole_to_mass(combaero::normalize_fractions(X_b));

  // Post-combustion enthalpy shift: h_out = h_ad + Q/m + f (h_ad - h_ref_b),
  // so dh_out/dT_ad = (1 + f) cp_ad.
  double dT_out_dT_ad = 1.0;
  double cp_out = combaero::cp_mass(T_ad, X_b);
  double T_out = T_ad;
  if (fraction != 0.0) {
    delta_h += fraction * (combaero::h_mass(T_ad, X_b) -
                           combaero::h_mass(kSensibleEnthalpyRefT, X_b));
  }
  if (std::abs(delta_h) > 0.0) {
    double h_ad = combaero::h_mass(T_ad, X_b);
    T_out = T_from_h_mass_continued(h_ad + delta_h, X_b, T_ad);
    cp_out = combaero::cp_mass(T_out, X_b);
    double cp_ad = combaero::cp_mass(T_ad, X_b);
    dT_out_dT_ad = (1.0 + fraction) * cp_ad / cp_out;
  }
  // The heat shift also moves T_out through the burnt composition at fixed
  // h_out: d(h_out) = (1+f) d(h_ad) - f d(h_ref) with h = sum Y_j h_j(T), so
  //   dT_out = [(1+f) cp_ad dT_ad + sum_j w_j dYb_j] / cp_out,
  //   w_j = (1+f) h_j(T_ad) - f h_j(T_ref) - h_j(T_out)   [J/kg].
  // Zero without a shift (T_out = T_ad, f = 0).
  std::vector<double> w_shift(n_species, 0.0);
  if (std::abs(delta_h) > 0.0 && cp_out > 0.0) {
    for (std::size_t j = 0; j < n_species; ++j) {
      const double Mj = combaero::species_molar_mass(j);
      const double h_ad_j = J_per_mol_to_J_per_kg(combaero::h_species(j, T_ad), Mj);
      const double h_ref_j =
          J_per_mol_to_J_per_kg(combaero::h_species(j, kSensibleEnthalpyRefT), Mj);
      const double h_out_j = J_per_mol_to_J_per_kg(combaero::h_species(j, T_out), Mj);
      w_shift[j] = ((1.0 + fraction) * h_ad_j - fraction * h_ref_j - h_out_j) / cp_out;
    }
  }
  const double T_ad_base = T_ad;  // the composition differences are of T_ad
  T_ad = T_out;

  std::size_t n_streams = streams.size();

  (void)T_ad_base;
  // Analytic (#716): X_mix(Y_mix) and Y_b(X_b) chained around dX_b/dX_mix.
  //   dX_i/dY_m = (delta_im - X_i) / (W_m S),  S = sum_l Y_l / W_l;
  //   dY_j/dX_l = (delta_jl W_j - Y_j W_l) / MW_b.
  std::vector<double> dT_ad_dY_mix(n_species, 0.0);
  std::vector<std::vector<double>> dYb_dY_mix(
      n_species, std::vector<double>(n_species, 0.0));
  {
    std::vector<double> W(n_species);
    double S = 0.0;
    const std::vector<double> Ym = combaero::normalize_fractions(mix.Y_mix);
    for (std::size_t m = 0; m < n_species; ++m) {
      W[m] = combaero::species_molar_mass(m);
      S += Ym[m] / W[m];
    }
    const double mw_b = combaero::mwmix(X_b);
    // dXb/dY_m = sum_k dXb/dX_k dX_k/dY_m, then to Y_b.
    for (std::size_t m = 0; m < n_species; ++m) {
      std::vector<double> dX(n_species);
      for (std::size_t i = 0; i < n_species; ++i) {
        dX[i] = (((i == m) ? 1.0 : 0.0) - X_mix[i]) / (W[m] * S);
      }
      double dT = 0.0;
      std::vector<double> dXb(n_species, 0.0);
      for (std::size_t k = 0; k < n_species; ++k) {
        dT += cad.dT_dX[k] * dX[k];
        for (std::size_t l = 0; l < n_species; ++l) {
          dXb[l] += cad.dXb_dX[l][k] * dX[k];
        }
      }
      dT_ad_dY_mix[m] = dT;
      double wdx = 0.0;
      for (std::size_t l = 0; l < n_species; ++l) {
        wdx += W[l] * dXb[l];
      }
      for (std::size_t j = 0; j < n_species; ++j) {
        dYb_dY_mix[j][m] = (W[j] * dXb[j] - Y_b[j] * wdx) / mw_b;
      }
    }
  }

  std::vector<double> dYb_dT_mix(n_species, 0.0);

  // Scale adiabatic derivatives by dT_out/dT_ad for post-combustion shift,
  // plus the burnt-composition term (w_shift above).
  dT_ad_dT_mix *= dT_out_dT_ad;
  for (std::size_t j = 0; j < n_species; ++j) {
    dT_ad_dT_mix += w_shift[j] * dYb_dT_mix[j];
  }
  for (std::size_t k = 0; k < n_species; ++k) {
    dT_ad_dY_mix[k] *= dT_out_dT_ad;
    for (std::size_t j = 0; j < n_species; ++j) {
      dT_ad_dY_mix[k] += w_shift[j] * dYb_dY_mix[j][k];
    }
  }

  MixerResult res;
  res.Tt_mix = T_ad;
  res.Pt_mix = mix.Pt_mix;
  res.Y_mix = Y_b;
  res.dTt_mix_d_delta_h = (cp_out > 0.0) ? (1.0 / cp_out) : 0.0;
  res.dTt_mix_dQ = res.dTt_mix_d_delta_h * phi_Q;
  res.dTt_mix_d_stream.resize(n_streams);
  res.dPt_mix_d_stream = mix.dPt_mix_d_stream;
  res.dY_mix_d_stream.assign(n_species, std::vector<StreamJacobian>(n_streams));

  for (std::size_t i = 0; i < n_streams; ++i) {
    StreamJacobian &dT_stream = res.dTt_mix_d_stream[i];
    dT_stream.d_Y.assign(n_species, 0.0);

    auto compute_chain = [&](double d_prop_T,
                             const std::vector<double> &d_prop_Y) {
      double sum = dT_ad_dT_mix * d_prop_T;
      for (std::size_t k = 0; k < n_species; ++k) {
        sum += dT_ad_dY_mix[k] * d_prop_Y[k];
      }
      return sum;
    };

    std::vector<double> dY_mix_dmdot(n_species), dY_mix_dT(n_species),
        dY_mix_dP_total(n_species);
    for (std::size_t k = 0; k < n_species; ++k) {
      dY_mix_dmdot[k] = mix.dY_mix_d_stream[k][i].d_mdot;
      dY_mix_dT[k] = mix.dY_mix_d_stream[k][i].d_Tt;
      dY_mix_dP_total[k] = mix.dY_mix_d_stream[k][i].d_Pt;
    }

    dT_stream.d_mdot =
        compute_chain(mix.dTt_mix_d_stream[i].d_mdot, dY_mix_dmdot);
    // Q is spread over the throughput: d(Q phi)/dm_i = Q phi', a pure
    // enthalpy shift of the products. Missing, a combustor carrying wall heat
    // had dT/dm = 0 where the difference gave 1658 K per kg/s.
    if (in_mix[i]) {
      dT_stream.d_mdot += res.dTt_mix_d_delta_h * Q * dphi_Q;
    }
    dT_stream.d_Tt = compute_chain(mix.dTt_mix_d_stream[i].d_Tt, dY_mix_dT);
    dT_stream.d_Pt =
        compute_chain(mix.dTt_mix_d_stream[i].d_Pt, dY_mix_dP_total);

    for (std::size_t j = 0; j < n_species; ++j) {
      std::vector<double> dY_j(n_species);
      for (std::size_t k = 0; k < n_species; ++k)
        dY_j[k] = mix.dY_mix_d_stream[k][i].d_Y[j];
      dT_stream.d_Y[j] = compute_chain(mix.dTt_mix_d_stream[i].d_Y[j], dY_j);
    }

    for (std::size_t tgt = 0; tgt < n_species; ++tgt) {
      StreamJacobian &dY_stream = res.dY_mix_d_stream[tgt][i];
      dY_stream.d_Y.assign(n_species, 0.0);

      auto compute_Y_chain = [&](double d_prop_T,
                                 const std::vector<double> &d_prop_Y) {
        double sum = dYb_dT_mix[tgt] * d_prop_T;
        for (std::size_t k = 0; k < n_species; ++k) {
          sum += dYb_dY_mix[tgt][k] * d_prop_Y[k];
        }
        return sum;
      };

      dY_stream.d_mdot =
          compute_Y_chain(mix.dTt_mix_d_stream[i].d_mdot, dY_mix_dmdot);
      dY_stream.d_Tt = compute_Y_chain(mix.dTt_mix_d_stream[i].d_Tt, dY_mix_dT);
      dY_stream.d_Pt =
          compute_Y_chain(mix.dTt_mix_d_stream[i].d_Pt, dY_mix_dP_total);

      for (std::size_t j = 0; j < n_species; ++j) {
        std::vector<double> dY_j(n_species);
        for (std::size_t k = 0; k < n_species; ++k)
          dY_j[k] = mix.dY_mix_d_stream[k][i].d_Y[j];
        dY_stream.d_Y[j] = compute_Y_chain(mix.dTt_mix_d_stream[i].d_Y[j], dY_j);
      }
    }
  }
  return res;
}

namespace {

// The burnt state's partials in the equilibrium's pressure argument P, which
// the stream partials of adiabatic_T_equilibrium_and_jacobians_from_streams
// hold fixed. The combustor kernel chains them through the mixed Pt (#1244).
struct EquilibriumPressurePartials {
  double dTt_dP = 0.0;
  std::vector<double> dY_dP;
};

MixerResult equilibrium_from_streams(const std::vector<Stream> &streams, double P, double Q,
                                     double fraction, EquilibriumPressurePartials *dP_out);

}  // namespace

MixerResult adiabatic_T_equilibrium_and_jacobians_from_streams(
    const std::vector<Stream> &streams, double P, double Q, double fraction) {
  return equilibrium_from_streams(streams, P, Q, fraction, nullptr);
}

namespace {

MixerResult equilibrium_from_streams(const std::vector<Stream> &streams, double P, double Q,
                                     double fraction, EquilibriumPressurePartials *dP_out) {
  // First mix adiabatically to get base state
  MixerResult mix = mixer_from_streams_and_jacobians(streams);

  // The flow the heat is spread over: that of the streams the mixer weights
  // (#1232).
  const std::vector<bool> in_mix = streams_in_mix(streams);
  double mdot_tot = 0.0;
  std::size_t n_species = combaero::num_species();
  for (std::size_t i_s = 0; i_s < streams.size(); ++i_s) {
    if (in_mix[i_s]) {
      mdot_tot += streams[i_s].m_dot;
    }
  }

  // Heat per unit flow: Q / |m|, fading out only near zero flow.
  const auto [phi_Q, dphi_Q] = mixer_heat_per_unit_flow(mdot_tot);

  // Q now; 'fraction' of the PRODUCTS' sensible enthalpy below, once T_ad
  // and the burnt composition are known (#481).
  double delta_h = Q * phi_Q;

  std::vector<double> X_mix = combaero::mass_to_mole(combaero::normalize_fractions(mix.Y_mix));
  const bool T_clamped = mix.Tt_mix < thermo::T_FLOOR;
  const bool P_clamped = P < kEquilibriumPressureFloor;
  State eq_in;
  eq_in.set_TPX(std::max(mix.Tt_mix, thermo::T_FLOOR), std::max(P, kEquilibriumPressureFloor),
                X_mix);
  const EquilibriumJacobian ej = combaero::combustion_equilibrium_and_jacobian(eq_in);
  double T_ad = ej.state.T;
  double dT_ad_dT_mix = T_clamped ? 0.0 : ej.dT_dT_in;
  std::vector<double> X_b = ej.state.X;
  std::vector<double> Y_b = combaero::mole_to_mass(combaero::normalize_fractions(X_b));

  // Post-combustion enthalpy shift: h_out = h_ad + Q/m + f (h_ad - h_ref_b),
  // so dh_out/dT_ad = (1 + f) cp_ad.
  double dT_out_dT_ad = 1.0;
  double cp_out = combaero::cp_mass(T_ad, X_b);
  double T_out = T_ad;
  if (fraction != 0.0) {
    delta_h += fraction * (combaero::h_mass(T_ad, X_b) -
                           combaero::h_mass(kSensibleEnthalpyRefT, X_b));
  }
  if (std::abs(delta_h) > 0.0) {
    double h_ad = combaero::h_mass(T_ad, X_b);
    T_out = T_from_h_mass_continued(h_ad + delta_h, X_b, T_ad);
    cp_out = combaero::cp_mass(T_out, X_b);
    double cp_ad = combaero::cp_mass(T_ad, X_b);
    dT_out_dT_ad = (1.0 + fraction) * cp_ad / cp_out;
  }
  // The heat shift also moves T_out through the burnt composition at fixed
  // h_out: d(h_out) = (1+f) d(h_ad) - f d(h_ref) with h = sum Y_j h_j(T), so
  //   dT_out = [(1+f) cp_ad dT_ad + sum_j w_j dYb_j] / cp_out,
  //   w_j = (1+f) h_j(T_ad) - f h_j(T_ref) - h_j(T_out)   [J/kg].
  // Zero without a shift (T_out = T_ad, f = 0).
  std::vector<double> w_shift(n_species, 0.0);
  if (std::abs(delta_h) > 0.0 && cp_out > 0.0) {
    for (std::size_t j = 0; j < n_species; ++j) {
      const double Mj = combaero::species_molar_mass(j);
      const double h_ad_j = J_per_mol_to_J_per_kg(combaero::h_species(j, T_ad), Mj);
      const double h_ref_j =
          J_per_mol_to_J_per_kg(combaero::h_species(j, kSensibleEnthalpyRefT), Mj);
      const double h_out_j = J_per_mol_to_J_per_kg(combaero::h_species(j, T_out), Mj);
      w_shift[j] = ((1.0 + fraction) * h_ad_j - fraction * h_ref_j - h_out_j) / cp_out;
    }
  }
  const double T_ad_base = T_ad;  // the composition differences are of T_ad
  T_ad = T_out;

  std::size_t n_streams = streams.size();

  (void)T_ad_base;
  // Analytic (#716): X_mix(Y_mix) and Y_b(X_b) chained around the
  // equilibrium's dX_out/dX and dX_out/dT_in.
  //   dX_i/dY_m = (delta_im - X_i) / (W_m S),  S = sum_l Y_l / W_l;
  //   dY_j/dX_l = (delta_jl W_j - Y_j W_l) / MW_b.
  // Per zero-flow stream: T_ad and Y_b along its own mixing direction.
  std::vector<bool> zero_flow_dir(streams.size(), false);
  std::vector<double> zero_flow_dT_ad(streams.size(), 0.0);
  std::vector<std::vector<double>> zero_flow_dYb(streams.size());
  std::vector<double> dT_ad_dY_mix(n_species, 0.0);
  std::vector<std::vector<double>> dYb_dY_mix(
      n_species, std::vector<double>(n_species, 0.0));
  std::vector<double> dYb_dT_mix(n_species, 0.0);
  {
    std::vector<double> W(n_species);
    double S = 0.0;
    const std::vector<double> Ym = combaero::normalize_fractions(mix.Y_mix);
    for (std::size_t m = 0; m < n_species; ++m) {
      W[m] = combaero::species_molar_mass(m);
      S += Ym[m] / W[m];
    }
    const double mw_b = combaero::mwmix(X_b);
    auto to_Y = [&](const std::vector<double> &dXb, std::size_t j) {
      double wdx = 0.0;
      for (std::size_t l = 0; l < n_species; ++l) wdx += W[l] * dXb[l];
      return (W[j] * dXb[j] - Y_b[j] * wdx) / mw_b;
    };
    for (std::size_t m = 0; m < n_species; ++m) {
      std::vector<double> dX(n_species);
      for (std::size_t i = 0; i < n_species; ++i) {
        dX[i] = (((i == m) ? 1.0 : 0.0) - X_mix[i]) / (W[m] * S);
      }
      double dT = 0.0;
      std::vector<double> dXb(n_species, 0.0);
      for (std::size_t k = 0; k < n_species; ++k) {
        dT += ej.dT_dX[k] * dX[k];
        for (std::size_t l = 0; l < n_species; ++l) dXb[l] += ej.dXout_dX[l][k] * dX[k];
      }
      dT_ad_dY_mix[m] = dT;
      for (std::size_t j = 0; j < n_species; ++j) dYb_dY_mix[j][m] = to_Y(dXb, j);
    }
    if (!T_clamped) {
      for (std::size_t j = 0; j < n_species; ++j) dYb_dT_mix[j] = to_Y(ej.dXout_dT_in, j);
    }
    // A stream at exactly zero flow (#1241): where it feeds a species the
    // burned state lacks, the equilibrium has a one-sided slope along its
    // mixing direction that no dT/dX matrix holds (in dry air the slopes
    // along two fuels do not add). Its flow columns are taken from the
    // equilibrium's derivatives along that direction alone; every other
    // column, and every stream with flow, keeps the matrix above.
    for (std::size_t i_s = 0; i_s < streams.size(); ++i_s) {
      if (streams[i_s].m_dot != 0.0) continue;
      std::vector<double> dY(n_species);
      bool any = false;
      for (std::size_t m = 0; m < n_species; ++m) {
        dY[m] = mix.dY_mix_d_stream[m][i_s].d_mdot;
        any = any || dY[m] != 0.0;
      }
      if (!any) continue;
      std::vector<double> dX(n_species, 0.0);
      for (std::size_t m = 0; m < n_species; ++m) {
        if (dY[m] == 0.0) continue;
        for (std::size_t i = 0; i < n_species; ++i) {
          dX[i] += (((i == m) ? 1.0 : 0.0) - X_mix[i]) / (W[m] * S) * dY[m];
        }
      }
      const EquilibriumJacobian ejd = combaero::combustion_equilibrium_and_jacobian(
          eq_in, true, combaero::SMOOTHING_K_PHI0, &dX);
      double dT = 0.0;
      std::vector<double> dXb(n_species, 0.0);
      for (std::size_t k = 0; k < n_species; ++k) {
        dT += ejd.dT_dX[k] * dX[k];
        for (std::size_t l = 0; l < n_species; ++l) dXb[l] += ejd.dXout_dX[l][k] * dX[k];
      }
      std::vector<double> dYb(n_species);
      for (std::size_t j = 0; j < n_species; ++j) dYb[j] = to_Y(dXb, j);
      zero_flow_dir[i_s] = true;
      zero_flow_dT_ad[i_s] = dT;
      zero_flow_dYb[i_s] = dYb;
    }
    if (dP_out != nullptr) {
      // T_out(T_ad(P), Y_b(P)) at the shift, as the T_mix chain below.
      dP_out->dY_dP.assign(n_species, 0.0);
      dP_out->dTt_dP = 0.0;
      if (!P_clamped) {
        for (std::size_t j = 0; j < n_species; ++j) dP_out->dY_dP[j] = to_Y(ej.dXout_dP, j);
        dP_out->dTt_dP = dT_out_dT_ad * ej.dT_dP;
        for (std::size_t j = 0; j < n_species; ++j) {
          dP_out->dTt_dP += w_shift[j] * dP_out->dY_dP[j];
        }
      }
    }
  }

  // Scale adiabatic derivatives by dT_out/dT_ad for post-combustion shift,
  // plus the burnt-composition term (w_shift above).
  dT_ad_dT_mix *= dT_out_dT_ad;
  for (std::size_t j = 0; j < n_species; ++j) {
    dT_ad_dT_mix += w_shift[j] * dYb_dT_mix[j];
  }
  for (std::size_t k = 0; k < n_species; ++k) {
    dT_ad_dY_mix[k] *= dT_out_dT_ad;
    for (std::size_t j = 0; j < n_species; ++j) {
      dT_ad_dY_mix[k] += w_shift[j] * dYb_dY_mix[j][k];
    }
  }

  MixerResult res;
  res.Tt_mix = T_ad;
  res.Pt_mix = mix.Pt_mix;
  res.Y_mix = Y_b;
  res.dTt_mix_d_delta_h = (cp_out > 0.0) ? (1.0 / cp_out) : 0.0;
  res.dTt_mix_dQ = res.dTt_mix_d_delta_h * phi_Q;
  res.dTt_mix_d_stream.resize(n_streams);
  res.dPt_mix_d_stream = mix.dPt_mix_d_stream;
  res.dY_mix_d_stream.assign(n_species, std::vector<StreamJacobian>(n_streams));

  for (std::size_t i = 0; i < n_streams; ++i) {
    StreamJacobian &dT_stream = res.dTt_mix_d_stream[i];
    dT_stream.d_Y.assign(n_species, 0.0);

    auto compute_chain = [&](double d_prop_T,
                             const std::vector<double> &d_prop_Y) {
      double sum = dT_ad_dT_mix * d_prop_T;
      for (std::size_t k = 0; k < n_species; ++k) {
        sum += dT_ad_dY_mix[k] * d_prop_Y[k];
      }
      return sum;
    };

    std::vector<double> dY_mix_dmdot(n_species), dY_mix_dT(n_species),
        dY_mix_dP_total(n_species);
    for (std::size_t k = 0; k < n_species; ++k) {
      dY_mix_dmdot[k] = mix.dY_mix_d_stream[k][i].d_mdot;
      dY_mix_dT[k] = mix.dY_mix_d_stream[k][i].d_Tt;
      dY_mix_dP_total[k] = mix.dY_mix_d_stream[k][i].d_Pt;
    }

    if (zero_flow_dir[i]) {
      double dT_dir = zero_flow_dT_ad[i] * dT_out_dT_ad;
      for (std::size_t j = 0; j < n_species; ++j) dT_dir += w_shift[j] * zero_flow_dYb[i][j];
      dT_stream.d_mdot = dT_ad_dT_mix * mix.dTt_mix_d_stream[i].d_mdot + dT_dir;
    } else {
      dT_stream.d_mdot = compute_chain(mix.dTt_mix_d_stream[i].d_mdot, dY_mix_dmdot);
    }
    // Q is spread over the throughput: d(Q phi)/dm_i = Q phi', a pure
    // enthalpy shift of the products. Missing, a combustor carrying wall heat
    // had dT/dm = 0 where the difference gave 1658 K per kg/s.
    if (in_mix[i]) {
      dT_stream.d_mdot += res.dTt_mix_d_delta_h * Q * dphi_Q;
    }
    dT_stream.d_Tt = compute_chain(mix.dTt_mix_d_stream[i].d_Tt, dY_mix_dT);
    dT_stream.d_Pt =
        compute_chain(mix.dTt_mix_d_stream[i].d_Pt, dY_mix_dP_total);

    for (std::size_t j = 0; j < n_species; ++j) {
      std::vector<double> dY_j(n_species);
      for (std::size_t k = 0; k < n_species; ++k)
        dY_j[k] = mix.dY_mix_d_stream[k][i].d_Y[j];
      dT_stream.d_Y[j] = compute_chain(mix.dTt_mix_d_stream[i].d_Y[j], dY_j);
    }

    for (std::size_t tgt = 0; tgt < n_species; ++tgt) {
      StreamJacobian &dY_stream = res.dY_mix_d_stream[tgt][i];
      dY_stream.d_Y.assign(n_species, 0.0);

      auto compute_Y_chain = [&](double d_prop_T,
                                 const std::vector<double> &d_prop_Y) {
        double sum = dYb_dT_mix[tgt] * d_prop_T;
        for (std::size_t k = 0; k < n_species; ++k) {
          sum += dYb_dY_mix[tgt][k] * d_prop_Y[k];
        }
        return sum;
      };

      dY_stream.d_mdot =
          zero_flow_dir[i]
              ? dYb_dT_mix[tgt] * mix.dTt_mix_d_stream[i].d_mdot + zero_flow_dYb[i][tgt]
              : compute_Y_chain(mix.dTt_mix_d_stream[i].d_mdot, dY_mix_dmdot);
      dY_stream.d_Tt = compute_Y_chain(mix.dTt_mix_d_stream[i].d_Tt, dY_mix_dT);
      dY_stream.d_Pt =
          compute_Y_chain(mix.dTt_mix_d_stream[i].d_Pt, dY_mix_dP_total);

      for (std::size_t j = 0; j < n_species; ++j) {
        std::vector<double> dY_j(n_species);
        for (std::size_t k = 0; k < n_species; ++k)
          dY_j[k] = mix.dY_mix_d_stream[k][i].d_Y[j];
        dY_stream.d_Y[j] = compute_Y_chain(mix.dTt_mix_d_stream[i].d_Y[j], dY_j);
      }
    }
  }
  return res;
}

}  // namespace

MixerResult combustor_residuals_and_jacobians(
    const std::vector<Stream> &streams, double Q, double fraction,
    const PressureLossCorrelation &pressure_loss, bool use_equilibrium) {

  // 1. Get unburned mixture state (T_in, P_in)
  MixerResult mix_in = mixer_from_streams_and_jacobians(streams);
  double T_in = mix_in.Tt_mix;
  // The burn's pressure is that mix's Pt (#1272): the stream partials below
  // (P_out's, and the equilibrium's through dP/d(stream)) are those of
  // mix_in.Pt_mix, so no other P would carry them rightly (dT/dPt_fuel
  // 1.3e-4 against an FD of 0 with P = 1e5 beside a mixed Pt of 2.03e5).
  if (!std::isfinite(mix_in.Pt_mix)) {
    throw std::invalid_argument(
        "combustor_residuals_and_jacobians: pressure must be finite and positive");
  }
  const double P = mix_in.Pt_mix;

  // 2. Get burned state (T_ad, Y_products) including all stream jacobians.
  // P is the streams' mixed total pressure (mix_in.Pt_mix, which the P_out
  // partials below already take it as), so the equilibrium's dependence on P
  // moves with each stream's flow and Pt through mix_in.dPt_mix_d_stream
  // (#1244). Without it, an equilibrium combustor near zero air flow, where
  // that mix's Pt moves 1.7e7 Pa per kg/s, was off by 3e-4 of its row.
  // Complete combustion does not depend on P.
  EquilibriumPressurePartials eq_dP;
  MixerResult res = use_equilibrium
      ? equilibrium_from_streams(streams, P, Q, fraction, &eq_dP)
      : adiabatic_T_complete_and_jacobian_T_from_streams(streams, P, Q, fraction);
  if (use_equilibrium) {
    const std::size_t n_sp = num_species();
    for (std::size_t i = 0; i < streams.size(); ++i) {
      const StreamJacobian &dPm = mix_in.dPt_mix_d_stream[i];
      StreamJacobian &dT = res.dTt_mix_d_stream[i];
      dT.d_mdot += eq_dP.dTt_dP * dPm.d_mdot;
      dT.d_Pt += eq_dP.dTt_dP * dPm.d_Pt;
      dT.d_Tt += eq_dP.dTt_dP * dPm.d_Tt;
      for (std::size_t k = 0; k < n_sp; ++k) dT.d_Y[k] += eq_dP.dTt_dP * dPm.d_Y[k];
      for (std::size_t j = 0; j < n_sp; ++j) {
        StreamJacobian &dY = res.dY_mix_d_stream[j][i];
        const double c = eq_dP.dY_dP[j];
        dY.d_mdot += c * dPm.d_mdot;
        dY.d_Pt += c * dPm.d_Pt;
        dY.d_Tt += c * dPm.d_Tt;
        for (std::size_t k = 0; k < n_sp; ++k) dY.d_Y[k] += c * dPm.d_Y[k];
      }
    }
  }

  double T_ad = res.Tt_mix;
  double theta = (T_in > 0.0) ? (T_ad / T_in - 1.0) : 0.0;

  // 3. Form Context and evaluate Hook
  State reactant_state;
  std::vector<double> X_in = mass_to_mole(mix_in.Y_mix);
  reactant_state.set_TPX(T_in, P, X_in);

  // Compute X_products from Y_products for the context
  std::vector<double> X_products = mass_to_mole(res.Y_mix);

  // The loss sees the flow the combustor burns (#1239): the total of the
  // streams in the mix, the one the mixer above weights and spreads Q over.
  // The signed sum of all streams shrank toward zero, or changed sign, while
  // gas still arrived beside a leaving stream. The network solver points
  // each stream it hands a node by its flow's sign, so there the two totals
  // agree (#1301 is the one known exception); this keeps a direct caller
  // consistent.
  const std::vector<bool> in_mix = streams_in_mix(streams);
  double mdot_total = 0.0;
  for (std::size_t i = 0; i < streams.size(); ++i) {
    if (in_mix[i]) {
      mdot_total += streams[i].m_dot;
    }
  }

  PressureLossContext ctx{
      reactant_state,
      0.0, // phi
      T_ad,
      X_products,
      res.Y_mix, // Y_products
      theta,
      0.0, 0.0, // mdot fuel/air
      mdot_total,
  };

  auto loss_res = pressure_loss(ctx);
  double zeta = std::get<0>(loss_res);
  double dzeta_dtheta = std::get<1>(loss_res);

  // 4. Update Result with P_out
  double P_out = P * (1.0 - zeta);
  res.Pt_mix = P_out;

  // 5. Analytical Chain Rule for dP_out/d_stream
  // dP_out/dx = dP/dx * (1 - zeta) - P * dzeta/dx
  // dzeta/dx = dzeta/dtheta * dtheta/dx
  // dtheta/dx = d(T_ad/T_in - 1)/dx = (1/T_in) * dT_ad/dx - (T_ad/T_in^2) * dT_in/dx

  std::size_t n_streams = streams.size();
  std::size_t n_species = num_species();

  for (std::size_t i = 0; i < n_streams; ++i) {
      auto update_jac = [&](StreamJacobian& dP_out_jac,
                          const StreamJacobian& dT_ad_jac,
                          const StreamJacobian& dT_in_jac,
                          const StreamJacobian& dP_in_jac) {

          auto compute_dtheta = [&](double d_Tad, double d_Tin) {
              return (T_in > 0.0) ? (d_Tad / T_in - (T_ad / (T_in * T_in)) * d_Tin) : 0.0;
          };

          double dtheta_dmdot = compute_dtheta(dT_ad_jac.d_mdot, dT_in_jac.d_mdot);
          double dtheta_dT = compute_dtheta(dT_ad_jac.d_Tt, dT_in_jac.d_Tt);
          double dtheta_dP = compute_dtheta(dT_ad_jac.d_Pt, dT_in_jac.d_Pt);

          dP_out_jac.d_mdot = dP_in_jac.d_mdot * (1.0 - zeta) - P * dzeta_dtheta * dtheta_dmdot;
          dP_out_jac.d_Tt = dP_in_jac.d_Tt * (1.0 - zeta) - P * dzeta_dtheta * dtheta_dT;
          dP_out_jac.d_Pt = dP_in_jac.d_Pt * (1.0 - zeta) - P * dzeta_dtheta * dtheta_dP;

          for (std::size_t k = 0; k < n_species; ++k) {
              double dtheta_dY = compute_dtheta(dT_ad_jac.d_Y[k], dT_in_jac.d_Y[k]);
              dP_out_jac.d_Y[k] = dP_in_jac.d_Y[k] * (1.0 - zeta) - P * dzeta_dtheta * dtheta_dY;
          }
      };

      update_jac(res.dPt_mix_d_stream[i],
                res.dTt_mix_d_stream[i], // dT_ad
                mix_in.dTt_mix_d_stream[i], // dT_in
                mix_in.dPt_mix_d_stream[i]); // dP_in
  }

  return res;
}

OrificeResult orifice_residuals_and_jacobian(double m_dot, double P_total_up,
                                             double P_static_up, double Tt_up,
                                             const std::vector<double> &Y_up,
                                             double P_static_down, double Cd,
                                             double area, double beta,
                                             double P_rho_ref) {
  (void)m_dot;
  // The drive follows approach_from_static: with E applied the approach head
  // is E's, so the drive is the upstream STATIC pressure (ISO's tapping
  // differential); without it, the upstream stagnation pressure.
  const bool from_static = approach_from_static(beta);
  const double dP = (from_static ? P_static_up : P_total_up) - P_static_down;
  // Density reference: the upstream static unless the caller names another.
  const bool own_rho_ref = P_rho_ref > 0.0;
  const double P_rho = own_rho_ref ? P_rho_ref : P_static_up;
  const std::vector<double> X_up = combaero::mass_to_mole(combaero::normalize_fractions(Y_up));
  auto [rho, drho_dT, drho_dP] = density_and_jacobians(Tt_up, P_rho, X_up);

  // Get mass flow and analytical Jacobians from low-level function
  auto [mdot_calc, dmdot_ddP, dmdot_drho] = orifice_mdot_and_jacobian(dP, rho, Cd, area, beta);

  OrificeResult res;
  res.m_dot_calc = mdot_calc;
  res.d_mdot_dP_total_up = from_static ? 0.0 : dmdot_ddP;
  res.d_mdot_dP_static_down = -dmdot_ddP;

  // Chain rule: d(mdot)/d(P_rho) = d(mdot)/d(rho) * d(rho)/d(P)
  const double d_rho_term = dmdot_drho * drho_dP;
  res.d_mdot_dP_static_up = (from_static ? dmdot_ddP : 0.0) + (own_rho_ref ? 0.0 : d_rho_term);
  res.d_mdot_dP_rho_ref = own_rho_ref ? d_rho_term : 0.0;
  res.d_mdot_dTt_up = dmdot_drho * drho_dT;

  // Composition, analytic (#964; a forward difference until then, off a
  // central one by up to 1e-5 of the row): d(mdot)/dY = d(mdot)/d(rho)
  // d(rho)/dY. rho = P MW / (R_u T) with MW = sum X_k W_k, so at fixed (T, P)
  // d(rho)/dX_j = rho W_j / MW, chained to the mass fractions.
  std::vector<double> drho_dX(X_up.size(), 0.0);
  const double mw = combaero::mwmix(X_up);
  for (std::size_t j = 0; j < X_up.size(); ++j) {
    drho_dX[j] = rho * combaero::species_molar_mass(j) / mw;
  }
  const std::vector<double> drho_dY = combaero::mole_gradient_to_mass(X_up, drho_dX);
  res.d_mdot_dY_up.resize(Y_up.size(), 0.0);
  for (std::size_t i = 0; i < Y_up.size(); ++i) {
    res.d_mdot_dY_up[i] = dmdot_drho * drho_dY[i];
  }

  return res;
}

ChannelResult channel_residuals_and_jacobian(double m_dot, double P_total_up,
                                       double P_static_up, double Tt_up,
                                       const std::vector<double> &Y_up,
                                       double P_static_down, double L, double D,
                                       double roughness,
                                       const std::string &friction_model,
                                       double f_multiplier,
                                       bool exit_head_lost) {
  (void)P_total_up;
  (void)P_static_down;
  const std::vector<double> X_up = combaero::mass_to_mole(combaero::normalize_fractions(Y_up));
  auto [rho, drho_dT, drho_dP] = density_and_jacobians(Tt_up, P_static_up, X_up);
  auto [mu, dmu_dT, dmu_dP] = viscosity_and_jacobians(Tt_up, P_static_up, X_up);

  double area = 0.25 * M_PI * D * D;
  double v = m_dot / (rho * area);
  double abs_v = std::abs(v);
  double Re = rho * abs_v * D / mu;

  auto [f, df_dRe] = friction_and_jacobian(friction_model, Re, roughness / D).result;
  f *= f_multiplier;
  df_dRe *= f_multiplier;

  double dP_calc = f * (L / D) * (0.5 * rho * v * abs_v);

  ChannelResult res;
  res.dP_calc = dP_calc;

  // dP = f(Re(|m|)) (L/D) m|m| / (2 rho A^2) is odd in m, so its slope is
  // even: |m| (f + |m| f'(Re) dRe/d|m| / 2) (L/D) / (rho A^2). Signed m here
  // gave reversed flow a slope of the wrong sign (#481).
  double dRe_dmdot = D / (mu * area);
  const double abs_m = std::abs(m_dot);
  if (abs_m > 0.0) {
    res.d_dP_d_mdot =
        (L / D) * (1.0 / (rho * area * area)) * abs_m * (f + 0.5 * abs_m * df_dRe * dRe_dmdot);
  } else {
    // At exactly zero flow |m| f is 0 * inf (f -> 64 / Re), and the formula
    // above read 0 where the drop's slope is laminar: dP = m (f Re) (L/D) /
    // (2 rho A^2 dRe/dm) with f Re finite, so d dP/dm = f_re_0 (L/D) /
    // (2 rho A^2 dRe/dm) = f_re_0 mu L / (2 rho A D^2), 32 mu L / (rho A D^2)
    // for Hagen-Poiseuille. The limit is the one-sided slope from both sides
    // (dP is odd in m). It omits the laminar branch's re_floor, which below
    // Re ~ 1e-6 (|m| ~ 1e-13 kg/s on a 1 cm air duct) turns the drop
    // quadratic: a floor for 64 / Re at Re = 0, not physics, and inside any
    // difference step (#1156; Newton stalled at an exact zero start). dP is
    // unchanged, so J here is the limit from |m| above the floor band, not
    // dP's own (zero) slope inside it. Forming the laminar drop floor-free,
    // m (1 - w) f_re / (dRe/dm), would remove the band but moves F by
    // round-off (19 of 3180 stress F arrays, 1.8e-15): a follow-up.
    res.d_dP_d_mdot = friction_f_re_at_zero(friction_model) * f_multiplier * (L / D) /
                      (2.0 * rho * area * area * dRe_dmdot);
  }

  const double eps_P = 1.0;
  auto [rho_p, drho_dT_p, drho_dP_p] = density_and_jacobians(Tt_up, P_static_up + eps_P, X_up);
  auto [mu_p, dmu_dT_p, dmu_dP_p] = viscosity_and_jacobians(Tt_up, P_static_up + eps_P, X_up);
  auto [rho_m, drho_dT_m, drho_dP_m] = density_and_jacobians(Tt_up, P_static_up - eps_P, X_up);
  auto [mu_m, dmu_dT_m, dmu_dP_m] = viscosity_and_jacobians(Tt_up, P_static_up - eps_P, X_up);
  double v_p = m_dot / (rho_p * area);
  double v_m = m_dot / (rho_m * area);
  double Re_p = rho_p * std::abs(v_p) * D / mu_p;
  double Re_m = rho_m * std::abs(v_m) * D / mu_m;
  auto [f_p, df_dRe_p] = friction_and_jacobian(friction_model, Re_p, roughness / D).result;
  auto [f_m, df_dRe_m] = friction_and_jacobian(friction_model, Re_m, roughness / D).result;
  f_p *= f_multiplier;
  f_m *= f_multiplier;
  double dP_p = f_p * (L / D) * (0.5 * rho_p * v_p * std::abs(v_p));
  double dP_m = f_m * (L / D) * (0.5 * rho_m * v_m * std::abs(v_m));
  res.d_dP_dP_static_up = (dP_p - dP_m) / (2.0 * eps_P);

  const double eps_T = 1e-3;
  auto [rho_T_p, drho_dT_T_p, drho_dP_T_p] = density_and_jacobians(Tt_up + eps_T, P_static_up, X_up);
  auto [mu_T_p, dmu_dT_T_p, dmu_dP_T_p] = viscosity_and_jacobians(Tt_up + eps_T, P_static_up, X_up);
  auto [rho_T_m, drho_dT_T_m, drho_dP_T_m] = density_and_jacobians(Tt_up - eps_T, P_static_up, X_up);
  auto [mu_T_m, dmu_dT_T_m, dmu_dP_T_m] = viscosity_and_jacobians(Tt_up - eps_T, P_static_up, X_up);
  double v_T_p = m_dot / (rho_T_p * area);
  double v_T_m = m_dot / (rho_T_m * area);
  double Re_T_p = rho_T_p * std::abs(v_T_p) * D / mu_T_p;
  double Re_T_m = rho_T_m * std::abs(v_T_m) * D / mu_T_m;
  auto [f_T_p, df_dRe_T_p] = friction_and_jacobian(friction_model, Re_T_p, roughness / D).result;
  auto [f_T_m, df_dRe_T_m] = friction_and_jacobian(friction_model, Re_T_m, roughness / D).result;
  f_T_p *= f_multiplier;
  f_T_m *= f_multiplier;
  double dP_T_p = f_T_p * (L / D) * (0.5 * rho_T_p * v_T_p * std::abs(v_T_p));
  double dP_T_m = f_T_m * (L / D) * (0.5 * rho_T_m * v_T_m * std::abs(v_T_m));
  res.d_dP_dTt_up = (dP_T_p - dP_T_m) / (2.0 * eps_T);

  res.d_dP_dY_up.resize(Y_up.size(), 0.0);
  const double eps_Y = 1e-6;
  for (std::size_t i = 0; i < Y_up.size(); ++i) {
    std::vector<double> Y_plus = Y_up;
    Y_plus[i] += eps_Y;
    auto X_plus = combaero::mass_to_mole(combaero::normalize_fractions(Y_plus));
    double rho_plus = std::get<0>(density_and_jacobians(Tt_up, P_static_up, X_plus));
    double mu_plus = std::get<0>(viscosity_and_jacobians(Tt_up, P_static_up, X_plus));
    double v_plus = m_dot / (rho_plus * area);
    double Re_plus = rho_plus * std::abs(v_plus) * D / mu_plus;
    auto [f_plus, dummy] = friction_and_jacobian(friction_model, Re_plus, roughness / D).result;
    f_plus *= f_multiplier;
    double dP_plus = f_plus * (L / D) * (0.5 * rho_plus * v_plus * std::abs(v_plus));
    res.d_dP_dY_up[i] = (dP_plus - dP_calc) / eps_Y;
    if (exit_head_lost) {
      // q = m|m| / (2 rho A^2) depends on Y only through rho.
      res.d_dP_dY_up[i] += 0.5 * m_dot * abs_m / (area * area) *
                           (1.0 / rho_plus - 1.0 / rho) / eps_Y;
    }
  }

  if (exit_head_lost) {
    // Exit dynamic head, lost at a stagnant sink (#588). Odd in m like the
    // friction drop; its slope |m| / (rho A^2) is even.
    const double q_exit = 0.5 * m_dot * abs_m / (rho * area * area);
    res.dP_calc += q_exit;
    res.d_dP_d_mdot += abs_m / (rho * area * area);
    res.d_dP_dP_static_up += -q_exit / rho * drho_dP;
    res.d_dP_dTt_up += -q_exit / rho * drho_dT;
  }

  return res;
}

using combaero::ChamberStaticState;
using combaero::chamber_static_state;

MomentumChamberResult momentum_chamber_residual_and_jacobian(
    double P, double Pt, double m_dot, double Tt,
    const std::vector<double> &Y, double area, bool incompressible) {
  // Momentum chamber: Pt = Pt(P, T, M), the isentropic stagnation
  // pressure of the chamber's static state. The node carries the MIXED
  // STAGNATION temperature Tt (its inflows mix as (m, Tt, Pt, Y) streams), so
  // the static T is solved here from h(Tt) = h(T) + v^2 / 2 (#576).

  // Soft normalization to handle edge cases (all-zero Y)
  std::vector<double> Y_safe = Y;
  double sum_Y = 0.0;
  for (double y : Y) {
    sum_Y += y;
  }
  // Use softmax-like approach: if sum is too small, use standard air
  const double eps = 1e-10;
  if (sum_Y < eps) {
    // Use standard dry air composition as fallback (N2=0.767, O2=0.233 mass fractions)
    Y_safe = std::vector<double>(Y.size(), 0.0);
    std::size_t n2_idx = combaero::species_index_from_name("N2");
    std::size_t o2_idx = combaero::species_index_from_name("O2");
    Y_safe[n2_idx] = 0.767;
    Y_safe[o2_idx] = 0.233;
  }

  const std::vector<double> X = combaero::mass_to_mole(combaero::normalize_fractions(Y_safe));

  MomentumChamberResult res;

  if (incompressible) {
    // Low-Mach limit (#588): T_static = Tt, rho = P / (R Tt), and
    // Pt = P + q with q = G^2 / (2 rho), G = |m_dot| / A. Even in m_dot like
    // the compressible form; q scales as Tt / P through 1 / rho.
    const double G = (area > 1e-12) ? std::abs(m_dot) / area : 0.0;
    const double sgn = (m_dot >= 0.0) ? 1.0 : -1.0;
    const double R_gas = combaero::specific_gas_constant(X);
    const double rho = P / (R_gas * Tt);
    const double q = 0.5 * G * G / rho;
    res.residual = Pt - P - q;
    res.d_res_dP_total = 1.0;
    res.d_res_dP = -1.0 + q / P;
    res.d_res_dT = -q / Tt;
    res.d_res_dmdot = (area > 1e-12) ? -(G / rho) * sgn / area : 0.0;
    res.T_static = Tt;
    res.dT_static_dT = 1.0;
    res.q_dyn = q;
    res.M = (G / rho) / combaero::speed_of_sound(Tt, X);
    res.rho_static = rho;
    res.drho_dP = rho / P;
    res.drho_dT = -rho / Tt;
    res.drho_dmdot = 0.0;
    // Y moves R alone here: q ~ R and rho ~ 1 / R, with dR/dY analytic
    // (specific_gas_constant_dY). Every species carries its partial, an
    // absent one included: a stream arriving at zero flow brings it, and the
    // relay moves it there (#1245).
    res.d_res_dY.assign(Y.size(), 0.0);
    res.drho_dY.assign(Y.size(), 0.0);
    if (sum_Y >= eps) {
      const std::vector<double> dR = combaero::specific_gas_constant_dY(X);
      for (std::size_t k = 0; k < Y.size(); ++k) {
        res.d_res_dY[k] = -q / R_gas * dR[k];
        res.drho_dY[k] = -rho / R_gas * dR[k];
      }
    }
    return res;
  }

  // Compressible stagnation closure (issue #357): Pt = P * exp((s(Tt) -
  // s(T)) / R), entropy conservation with variable cp, rather than the
  // incompressible Pt = P + 0.5*rho*v^2, which under-predicts the stagnation
  // rise without bound past sonic.
  //
  // Until #576 the node's Tt was used here AS the static T: rho, a and M all
  // came from Tt, so the static state was too hot by Tt / T = 1 + (g-1)/2 M^2.
  // At fixed Pt that capped the mass flux at 92% of choked (peak at M 0.913),
  // and Newton stalled on any chamber asked to pass more. Solving the static
  // T from the energy equation removes the cap.
  //
  // The stagnation rise does not depend on flow direction: the mass flux is
  // taken as |m_dot| / A and the sign is carried back into d(res)/d(m_dot),
  // keeping the residual even in m_dot.
  const double G = (area > 1e-12) ? std::abs(m_dot) / area : 0.0;
  const ChamberStaticState st = chamber_static_state(Tt, P, G, X);
  const double Pt_cl = P * std::exp(st.lnPtP);
  const double sgn = (m_dot >= 0.0) ? 1.0 : -1.0;

  res.residual = Pt - Pt_cl;
  res.d_res_dP_total = 1.0;
  res.T_static = st.T;
  res.dT_static_dP = st.dT_dP;
  res.dT_static_dT = st.dT_dTt;
  res.dT_static_dmdot = (area > 1e-12) ? st.dT_dG * sgn / area : 0.0;
  res.q_dyn = P * std::expm1(st.lnPtP);
  res.M = st.M;
  // Static density, ideal gas: rho = P / (R T), T = T(Tt, P, G).
  const double R_gas = combaero::specific_gas_constant(X);
  const double rho = P / (R_gas * st.T);
  res.rho_static = rho;
  res.drho_dP = rho / P - (rho / st.T) * st.dT_dP;
  res.drho_dT = -(rho / st.T) * st.dT_dTt;
  res.drho_dmdot = (area > 1e-12) ? -(rho / st.T) * st.dT_dG * sgn / area : 0.0;
  res.d_res_dP = -(Pt_cl / P + Pt_cl * st.dln_dP);
  res.d_res_dmdot = (area > 1e-12) ? -Pt_cl * st.dln_dG * sgn / area : 0.0;
  res.d_res_dT = -Pt_cl * st.dln_dTt;

  // Branch selection (#788). At fixed (T0, G) the closure Pt_cl(P) falls to
  // its minimum Pt* at the sonic static pressure P* and rises again below it,
  // so one Pt above Pt* has two static roots: the subsonic one and a
  // supersonic one. No element here accelerates a stream past Mach 1 into a
  // chamber, and Newton landed on the supersonic root ("converged" junctions
  // at Mach 1.46). Below P* the closure is therefore continued by its mirror
  // image about Pt*, g(P) = 2 Pt* - Pt_cl(P): g rises monotonically through
  // P* (C1 there, both slopes zero), so every Pt has exactly ONE root, the
  // subsonic one when Pt >= Pt*. A Pt below Pt* (more flow than the chamber's
  // area carries at that stagnation pressure, a choked duct discharging into
  // it) roots on the continuation: P < P*, M > 1, `bounded`. Such a root is
  // the choked-bounded junction the solver reports as such, or refuses.
  // Pt*'s partials need no dP*/dx: dPt_cl/dP = 0 at P* (the envelope
  // theorem), so dPt*/dx is the closure's own partial at P*.
  res.bounded = st.M > 1.0;
  if (res.bounded) {
    const combaero::SonicStaticState son = combaero::sonic_static_state(Tt, G, X);
    const ChamberStaticState ss = chamber_static_state(Tt, son.P, G, X);
    const double Pt_star = son.P * std::exp(ss.lnPtP);
    res.P_sonic = son.P;
    res.Pt_sonic = Pt_star;
    res.residual = Pt - (2.0 * Pt_star - Pt_cl);
    res.d_res_dP = -res.d_res_dP;
    res.d_res_dmdot =
        (area > 1e-12) ? (Pt_cl * st.dln_dG - 2.0 * Pt_star * ss.dln_dG) * sgn / area : 0.0;
    res.d_res_dT = Pt_cl * st.dln_dTt - 2.0 * Pt_star * ss.dln_dTt;
  }

  // Composition, analytic (#743): Y moves R, h, cp and s together, and the
  // relay reaches this row only through these entries (#517, which used a
  // central difference here). chamber_static_state_dX gives the static T and
  // ln(Pt/P) at fixed Tt, P, G; rho = P MW / (R_u T).
  res.d_res_dY.assign(Y.size(), 0.0);
  res.drho_dY.assign(Y.size(), 0.0);
  res.dT_static_dY.assign(Y.size(), 0.0);
  if (sum_Y >= eps) {
    const ChamberStaticStateDX cx = chamber_static_state_dX(Tt, P, G, X);
    const double MW = combaero::mwmix(X) / 1000.0;
    std::vector<double> g_res(X.size()), g_rho(X.size());
    for (std::size_t k = 0; k < X.size(); ++k) {
      const double W = combaero::species_molar_mass(k) / 1000.0;
      g_res[k] = -Pt_cl * cx.dlnPtP_dX[k];
      g_rho[k] = rho * (W / MW - cx.dT_dX[k] / st.T);
    }
    if (res.bounded) {
      // The continuation 2 Pt* - Pt_cl: Pt* moves with X at fixed P*
      // (envelope, as above).
      const ChamberStaticStateDX cs = chamber_static_state_dX(Tt, res.P_sonic, G, X);
      for (std::size_t k = 0; k < X.size(); ++k) {
        g_res[k] = Pt_cl * cx.dlnPtP_dX[k] - 2.0 * res.Pt_sonic * cs.dlnPtP_dX[k];
      }
    }
    res.d_res_dY = combaero::mole_gradient_to_mass(X, g_res);
    res.drho_dY = combaero::mole_gradient_to_mass(X, g_rho);
    res.dT_static_dY = combaero::mole_gradient_to_mass(X, cx.dT_dX);
    // A species absent at the node keeps its partial: a stream arriving at
    // exactly zero flow is in the mix with weight 0 and the relay moves its
    // species there, as it enters (#1245). Zeroing them gave J a hybrid of
    // the two sides: the arriving side's T and present species, no absent
    // ones (1.3e-3 scaled on a tee's port chamber).
  }

  return res;
}

// -----------------------------------------------------------------------------
// Compressible Flow Components
// -----------------------------------------------------------------------------

namespace {

// d(mdot)/dX_j of the isentropic nozzle at fixed (T0, P0, P_out), the mole
// fractions taken unnormalised (chain with mole_gradient_to_mass) (#964).
// mdot_fwd > 0 is the forward flow, T the outlet temperature at P_out.
//   Outlet on the isentrope through (T0, P0): the mixing entropy is the same
//   at both ends, so sum_k X_k [s_k(T) - s_k(T0)] = R_u ln(P_out / P0) and
//     dT/dX_j = -(s_j(T) - s_j(T0)) T / cp(T)               (molar cp);
//   a j-independent term (the pressure ratio's) drops out of the chain.
//   mdot = A_eff rho v, rho = P_out MW / (R_u T), v = sqrt(2 dH),
//   dH = sum_k X_k [h_k(T0) - h_k(T)] / MW:
//     d ln rho / dX_j = W_j / MW - (dT/dX_j) / T,
//     d dH / dX_j = (dh_j - dH W_j) / MW - cp_mass(T) dT/dX_j.
// Choked, P_out = P* maximises G = rho v over the outlet pressure, so
// dG*/dX = dG/dX at fixed P* (envelope theorem): the same expression at the
// throat, and nozzle_flow's golden-section P* enters only at second order.
std::vector<double> nozzle_mdot_dX(double mdot_fwd, double T0, double T,
                                   const std::vector<double>& X) {
  const std::size_t n = X.size();
  std::vector<double> out(n, 0.0);
  if (!(mdot_fwd > 0.0)) {
    return out;
  }
  const double MW = combaero::mwmix(X);  // g/mol, as species_molar_mass
  const double cp_T = combaero::cp(T, X);  // J/(mol K)
  std::vector<double> dh(n, 0.0);
  double dH = 0.0;  // J/mol of mixture
  for (std::size_t j = 0; j < n; ++j) {
    dh[j] = combaero::h_species(j, T0) - combaero::h_species(j, T);
    dH += X[j] * dh[j];
  }
  if (!(dH > 0.0)) {
    return out;
  }
  // Per mole of mixture throughout: dH / MW is the mass-specific drop and
  // cp_T / MW the mass cp, so the MW factors cancel in d ln dH.
  for (std::size_t j = 0; j < n; ++j) {
    const double W = combaero::species_molar_mass(j);
    const double dT =
        -(combaero::s_species(j, T) - combaero::s_species(j, T0)) * T / cp_T;
    const double dln_rho = W / MW - dT / T;
    const double dln_dH = (dh[j] - dH * W / MW - cp_T * dT) / dH;
    out[j] = mdot_fwd * (dln_rho + 0.5 * dln_dH);
  }
  return out;
}

// The isentropic nozzle's flow and partials in its forward frame (P_back <=
// P0): mdot = A_eff G(P0, T0, P_back), with nozzle_flow's composition
// gradient (nozzle_mdot_dX, unnormalised X) when dm_dX is given.
struct NozzleForward {
  double mdot = 0.0;
  double dP0 = 0.0;
  double dP_back = 0.0;
  double dT0 = 0.0;
};

NozzleForward nozzle_forward(double Tt, double P0_fwd, double P_back_fwd,
                             const std::vector<double>& X, double A_eff,
                             std::vector<double>* dm_dX) {
  auto sol = combaero::nozzle_flow(Tt, P0_fwd, P_back_fwd, A_eff, X);

  // Compute Jacobians analytically using isentropic relations
  // G = rho * v = (P/RT) * sqrt(2*(h0 - h))
  // For isentropic flow: dT/dP = 1/(rho*cp), d_rho/dP = 1/a^2
  // dG/dP = 1/v * (M^2 - 1)

  double T_outlet = sol.outlet.T;
  double v_outlet = sol.v;
  double rho_outlet = sol.outlet.rho();
  double M_outlet = sol.M;
  double cp_outlet = combaero::cp_mass(T_outlet, X);
  double cp0 = combaero::cp_mass(Tt, X);

  // Back-pressure derivative: the exact derivative of the value as written.
  // Unchoked, mdot = A_eff G(Pb) on the isentrope, dG/dPb = (M^2 - 1) / v.
  // Choked, the outlet sits at P* whatever Pb is, so dmdot/dPb = 0. The value
  // is already C1 through P*: G peaks there, so the unchoked slope (M^2 - 1)/v
  // falls to zero as M -> 1 and meets the choked branch's zero. No blend is
  // needed. One was here until #772, a cubic in PR over 0.01 either side of
  // the critical ratio applied to the derivative only: above P* it scaled the
  // slope by a factor that reached 1.5 at the band's edge and jumped back to 1
  // there, so d(mdot)/dPb was up to 46% off a difference of the value (and
  // d(mdot)/dPt, which reuses it, up to 0.3%).
  double d_mdot_dP_back_fwd = 0.0;
  if (!sol.choked && v_outlet > 1e-6) {
    d_mdot_dP_back_fwd = A_eff * (1.0 / v_outlet) * (M_outlet * M_outlet - 1.0);
  }
  const double PR = P_back_fwd / P0_fwd;

  // Upstream pressure derivative: Pt only affects PR and scales G linearly if choked
  // dmdot/dP0 = mdot/Pt - (Pb/Pt) * dmdot/dPb
  double d_mdot_dP0_fwd = (P0_fwd > 1e-6) ? (sol.mdot / P0_fwd - PR * d_mdot_dP_back_fwd) : 0.0;

  // Upstream temperature derivative: mdot ~ 1/sqrt(Tt) roughly
  // More precisely: dG/dT0 = rho * cp0 / v * [(1 - T/Tt) - v^2/(cp*Tt)]
  double d_mdot_dT0_fwd = 0.0;
  if (v_outlet > 1e-6 && Tt > 1e-6) {
    d_mdot_dT0_fwd = A_eff * rho_outlet * cp0 / v_outlet * ((1.0 - T_outlet / Tt) - v_outlet * v_outlet / (cp_outlet * Tt));
  } else if (Tt > 1e-6) {
    // Choked limit: dG/dT0 = -G / (2*Tt)
    d_mdot_dT0_fwd = -sol.mdot / (2.0 * Tt);
  }

  if (dm_dX != nullptr) {
    // nozzle_flow evaluates at its floored stagnation state.
    *dm_dX = nozzle_mdot_dX(sol.mdot, sol.stagnation.T, T_outlet, X);
  }
  return {sol.mdot, d_mdot_dP0_fwd, d_mdot_dP_back_fwd, d_mdot_dT0_fwd};
}

// orifice_compressible_mdot_and_jacobian, with the composition gradient of
// its flow (nozzle_mdot_dX, unnormalised X) when dm_dX is given.
std::tuple<double, double, double, double> orifice_compressible_mdot_impl(
    double Tt, double Pt, double P_back,
    const std::vector<double>& X,
    double Cd, double area, double beta, std::vector<double>* dm_dX) {

  // Compute effective area with velocity-of-approach factor
  const double E = velocity_of_approach_factor(beta);
  double A_eff = Cd * E * area;

  // Handle reverse flow (P_back > Pt)
  bool reverse_flow = (P_back > Pt);
  double T0_fwd = Tt;  // Temperature same in both directions
  double P0_fwd = reverse_flow ? P_back : Pt;
  double P_back_fwd = reverse_flow ? Pt : P_back;

  // Below kJetDriveFloor of drive, the flow is the floor's times the odd
  // cubic h(t), t = drive / floor (flow_drive_floor.h; #1256): the jet's
  // sqrt(drive) has an infinite slope at zero, and nozzle_flow's 1 - PR is a
  // round-off staircase near P0 (J flipped between 0 and -22 on a 1e-10 Pa
  // change at 1 bar). The floor's flow G0 is at P_back = P0 - floor, and the
  // cubic's exponent is the jet's, kJetFloorExponent = 1/2, a constant, so
  //   dm/dPb = -G0 h'(t) / floor,
  //   dm/dP0 =  G0 h'(t) / floor + h (dG0/dP0 + dG0/dPb),
  //   dm/dT0 = h dG0/dT0, dm/dX = h dG0/dX
  // are exact. The flow is continuous at the floor, but its slope steps
  // there by the nozzle's own exponent's distance from 1/2: measured, 4e-7 at
  // 1 bar and 6.1e-5 at 100 bar.
  const double drive = P0_fwd - P_back_fwd;
  const bool sub_floor = drive < kJetDriveFloor && P0_fwd > kJetDriveFloor;
  NozzleForward fwd;
  if (sub_floor) {
    const NozzleForward at =
        nozzle_forward(T0_fwd, P0_fwd, P0_fwd - kJetDriveFloor, X, A_eff, dm_dX);
    const DriveFloorCubic c = drive_floor_cubic(drive / kJetDriveFloor, kJetFloorExponent);
    fwd.mdot = at.mdot * c.h;
    fwd.dP_back = -at.mdot * c.dh / kJetDriveFloor;
    fwd.dP0 = at.mdot * c.dh / kJetDriveFloor + c.h * (at.dP0 + at.dP_back);
    fwd.dT0 = c.h * at.dT0;
    if (dm_dX != nullptr) {
      for (double& v : *dm_dX) {
        v *= c.h;
      }
    }
  } else {
    fwd = nozzle_forward(T0_fwd, P0_fwd, P_back_fwd, X, A_eff, dm_dX);
  }
  double mdot = reverse_flow ? -fwd.mdot : fwd.mdot;

  // Assign results
  double d_mdot_dP0 = fwd.dP0;
  double d_mdot_dP_back = fwd.dP_back;
  double d_mdot_dT0 = fwd.dT0;

  // For reverse flow, swap Jacobian signs appropriately.
  // mdot = -mdot_fwd(P_back, Pt), so all three derivatives are negated,
  // and Pt/P_back roles are swapped.
  if (reverse_flow) {
    std::swap(d_mdot_dP0, d_mdot_dP_back);
    d_mdot_dP0 = -d_mdot_dP0;
    d_mdot_dP_back = -d_mdot_dP_back;
    d_mdot_dT0 = -d_mdot_dT0;
  }

  if (dm_dX != nullptr && reverse_flow) {
    for (double& v : *dm_dX) {
      v = -v;
    }
  }

  return {mdot, d_mdot_dP0, d_mdot_dP_back, d_mdot_dT0};
}

}  // namespace

std::tuple<double, double, double, double> orifice_compressible_mdot_and_jacobian(
    double Tt, double Pt, double P_back,
    const std::vector<double>& X,
    double Cd, double area, double beta) {
  return orifice_compressible_mdot_impl(Tt, Pt, P_back, X, Cd, area, beta, nullptr);
}

OrificeResult orifice_compressible_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double Tt_up,
    const std::vector<double>& Y_up,
    double P_static_down, double Cd, double area, double beta) {

  (void)m_dot;  // Not used in compressible orifice (computed from pressures)
  const std::vector<double> X_up = combaero::mass_to_mole(combaero::normalize_fractions(Y_up));

  // The same rule as the incompressible kernel: with E applied the approach
  // head is E's, so the nozzle expands from the upstream STATIC state, the
  // compressible form of ISO's rho1 at the upstream tap. The static
  // temperature lies on the isentrope through (Tt_up, P_total_up): pairing
  // the static pressure with the stagnation temperature made the flow
  // sqrt(T1/Tt) low behind a node with velocity, 0.3% at pipe Mach 0.18
  // (#645). Where P_static_up = P_total_up (a plenum) T1 = Tt_up exactly.
  const bool from_static = approach_from_static(beta);
  const double P0 = from_static ? P_static_up : P_total_up;
  // Tight: the Jacobian tests check this column to 1e-7.
  const double static_T_tol = 1e-10;
  auto static_T = [&](const std::vector<double>& X) {
    if (!from_static || P_static_up == P_total_up) {
      return Tt_up;
    }
    return combaero::calc_T_from_s(combaero::s(Tt_up, X, P_total_up), P_static_up, X, Tt_up,
                                   static_T_tol);
  };
  const double T1 = static_T(X_up);
  // Along the isentrope, ds = cp/T dT - R/P dP = 0 at both ends:
  //   dT1/dTt = (cp(Tt)/Tt) / (cp(T1)/T1), dT1/dPt = -(R/Pt) / (cp(T1)/T1),
  //   dT1/dP1 = (R/P1) / (cp(T1)/T1).
  double dT1_dTt = 1.0;
  double dT1_dPt = 0.0;
  double dT1_dP1 = 0.0;
  // Applied at P_static_up = P_total_up too: only the value coincides there,
  // and a plenum perturbs its P and Pt unknowns separately.
  if (from_static) {
    const double c1 = combaero::cp(T1, X_up) / T1;
    dT1_dTt = (combaero::cp(Tt_up, X_up) / Tt_up) / c1;
    dT1_dPt = -(combaero::thermo::R_GAS / P_total_up) / c1;
    dT1_dP1 = (combaero::thermo::R_GAS / P_static_up) / c1;
  }

  // Get mass flow and Jacobians from low-level function
  std::vector<double> dm_dX;
  auto [mdot_calc, dmdot_dP0, dmdot_dP_back, dmdot_dT0] =
      orifice_compressible_mdot_impl(T1, P0, P_static_down, X_up, Cd, area, beta, &dm_dX);

  OrificeResult res;
  res.m_dot_calc = mdot_calc;
  res.d_mdot_dP_total_up = from_static ? dmdot_dT0 * dT1_dPt : dmdot_dP0;
  res.d_mdot_dP_static_down = dmdot_dP_back;
  res.d_mdot_dP_static_up = from_static ? dmdot_dP0 + dmdot_dT0 * dT1_dP1 : 0.0;
  res.d_mdot_dP_rho_ref = 0.0;
  res.d_mdot_dTt_up = dmdot_dT0 * dT1_dTt;

  // Composition, analytic (#964; a forward difference until then, off a
  // central one by up to 7e-4 of the row): the nozzle's own gradient at fixed
  // (T1, P0, Pb), plus T1's where it is the approach static temperature. On
  // the isentrope through (Tt_up, P_total_up), as nozzle_mdot_dX,
  //   dT1/dX_j = -(s_j(T1) - s_j(Tt_up)) T1 / cp(T1).
  if (from_static && P_static_up != P_total_up) {
    const double c1 = combaero::cp(T1, X_up);
    for (std::size_t j = 0; j < X_up.size(); ++j) {
      const double dT1 =
          -(combaero::s_species(j, T1) - combaero::s_species(j, Tt_up)) * T1 / c1;
      dm_dX[j] += dmdot_dT0 * dT1;
    }
  }
  res.d_mdot_dY_up = combaero::mole_gradient_to_mass(X_up, dm_dX);

  return res;
}

// -----------------------------------------------------------------------------
// Opt-in compressible orifice models (#854): Jobson (1955), ISO 5167-2 Eq. (5)
// -----------------------------------------------------------------------------

namespace {

// Partial indices: upstream total and approach static pressure, upstream total
// temperature, back pressure, the element's incompressible coefficient, and
// the two numbers composition enters by: the jet exponent n = gamma(Tt) and
// the specific gas constant R (#964).
enum : int { kMdPt = 0, kMdP1 = 1, kMdTt = 2, kMdPb = 3, kMdCd = 4, kMdN = 5, kMdR = 6,
             kMdCount = 7 };
using MD = DualN<kMdCount>;

// gamma(Tt) carrying d/dTt: gamma = cp / (cp - R), so
// dgamma/dT = -R cp'(T) / (cp - R)^2 (molar cp; R the universal constant).
// Seeded in its own direction too: the composition chain (model_result).
MD gamma_of_Tt(double Tt, const std::vector<double>& X) {
  const double c = combaero::cp(Tt, X);
  const double cv = c - combaero::thermo::R_GAS;
  MD n = MD::constant(c / cv);
  n.d[kMdTt] = -combaero::thermo::R_GAS * combaero::dcp_dT(Tt, X) / (cv * cv);
  n.d[kMdN] = 1.0;
  return n;
}

struct ModelEval {
  MD m = MD::constant(0.0);
  MD eps = MD::constant(1.0);
  double Ci = 0.0;
  double C = 0.0;
  bool choked = false;
  bool in_range = true;
};

ModelEval jobson_eval_raw(double Pt, double P1, double Tt, double Pb, double Cd, double area,
                          double beta, double Cv, const std::vector<double>& X) {
  const bool from_static = approach_from_static(beta);
  const double b = from_static ? beta : 0.0;
  const MD pt = MD::seed(Pt, kMdPt);
  const MD tt = MD::seed(Tt, kMdTt);
  const MD pb = MD::seed(Pb, kMdPb);
  const MD cd = MD::seed(Cd, kMdCd);
  // The drive: the approach static with E, the reservoir's total without.
  const MD pd = from_static ? MD::seed(P1, kMdP1) : pt;
  const MD n = gamma_of_Tt(Tt, X);
  const MD R = MD::seed(combaero::specific_gas_constant(X), kMdR);
  const MD ci = orifice_jobson::ci_from_cd(cd, b, Cv);
  const double E = velocity_of_approach_factor(b);

  ModelEval out;
  out.Ci = ci.v;
  out.C = ci.v;
  bool degenerate = false;
  if (pb.v < pd.v) {
    const MD rho0 = pt / (R * tt);
    const auto fl = orifice_jobson::flux(pt, rho0, pd, pb, ci, n, b, Cv);
    out.m = Cv * area * fl.g;
    out.C = fl.C.v;
    out.choked = fl.choked;
    degenerate = fl.degenerate;
    const MD rho1 = rho0 * orifice_jobson::pow_dd(pd / pt, 1.0 / n);
    out.eps = out.m / (cd * (E * area) * dsqrt(2.0 * rho1 * (pd - pb)));
  } else if (pb.v > pd.v) {
    // Reverse band of the raw kernel: a reservoir at the back pressure.
    const MD rho0 = pb / (R * tt);
    const auto fl = orifice_jobson::flux(pb, rho0, pb, pd, ci, n, b, Cv);
    out.m = -(Cv * area) * fl.g;
    out.C = fl.C.v;
    out.choked = fl.choked;
    degenerate = fl.degenerate;
    out.eps = out.m / (-(cd * (E * area)) * dsqrt(2.0 * rho0 * (pb - pd)));
  }
  namespace jb = ::orifice::jobson1955;
  out.in_range = ci.v >= jb::ci_min && ci.v <= jb::ci_max && b <= jb::beta_max && !degenerate;
  return out;
}

// Jobson's flow with the network's zero-drive floor (flow_drive_floor.h;
// #1256): below kJetDriveFloor of drive |pd - pb| (pd the approach static
// with E, else Pt) the flow is the floor's times the odd cubic h(t), t =
// drive / floor, as the nozzle kernel. The floor's flow m0 is evaluated with
// the back pressure (forward) or the drive pressure (reverse) moved to one
// floor from the other, so its partials carry that pressure's onto the one
// it was moved from. The cubic's exponent is the constant kJetFloorExponent,
// so every partial is exact: matching Jobson's own exponent at the floor was
// tried and moved with the approach head q = Pt - P1 at the floor's scale
// (0.5 at q = 0, 0.57 at 1e-4 Pa, 0.81 at 1e-3, beta 0.5), so holding it fixed
// left dm/dPt 22-27% and dm/dP1 7-15% off in the band (#1259 review). The
// flow is continuous at the floor, but its slope steps there by Jobson's
// exponent's distance from 1/2. In a solved state q at the floor is about
// Cd^2 beta^4 / (1 - beta^4) times the drive, 5e-6 to 1e-4 Pa: in solved
// pipe-fed states the step in dm/dPb measured 0.19% (beta 0.5, Cd 0.27), 3.9%
// (beta 0.5, Cd 0.6) and 27% (beta 0.7, Cd 0.6), with no effect on
// convergence in any run. Exact in-band partials are the trade: a matched,
// C1 exponent would need the jet law's second derivatives. At exactly zero
// drive the flow is 0 with the forward cubic's slope.
ModelEval jobson_eval(double Pt, double P1, double Tt, double Pb, double Cd, double area,
                      double beta, double Cv, const std::vector<double>& X) {
  const bool from_static = approach_from_static(beta);
  const double pd = from_static ? P1 : Pt;
  const int kd = from_static ? kMdP1 : kMdPt;
  const bool reverse = Pb > pd;
  const double drive = reverse ? Pb - pd : pd - Pb;
  const double p_up = reverse ? Pb : pd;
  if (!(drive < kJetDriveFloor) || !(p_up > kJetDriveFloor)) {
    return jobson_eval_raw(Pt, P1, Tt, Pb, Cd, area, beta, Cv, X);
  }
  ModelEval at;
  int moved = kMdPb;    // the pressure set one floor from its partner
  int partner = kd;     // the pressure it was set from
  if (reverse) {
    const double p = Pb - kJetDriveFloor;
    at = jobson_eval_raw(from_static ? Pt : p, from_static ? p : P1, Tt, Pb, Cd, area, beta, Cv,
                         X);
    moved = kd;
    partner = kMdPb;
  } else {
    at = jobson_eval_raw(Pt, P1, Tt, pd - kJetDriveFloor, Cd, area, beta, Cv, X);
  }
  // m0 along the floor: its partial by the moved pressure belongs to the
  // partner, which the moved one follows.
  MD m0 = at.m;
  m0.d[partner] += m0.d[moved];
  m0.d[moved] = 0.0;
  // t = drive / floor as a dual number in the two pressures.
  MD t = MD::constant(drive / kJetDriveFloor);
  t.d[kd] = (reverse ? -1.0 : 1.0) / kJetDriveFloor;
  t.d[kMdPb] = (reverse ? 1.0 : -1.0) / kJetDriveFloor;
  const DriveFloorCubic c = drive_floor_cubic(t.v, kJetFloorExponent);
  MD h = MD::constant(c.h);
  for (int i = 0; i < kMdCount; ++i) {
    h.d[i] = c.dh * t.d[i];
  }
  ModelEval out = at;
  out.m = m0 * h;
  return out;
}

ModelEval iso5167_eval(double Pt, double P1, double Tt, double Pb, double Cd, double area,
                       double beta, const std::vector<double>& X) {
  const bool from_static = approach_from_static(beta);
  const double b = from_static ? beta : 0.0;
  const MD pt = MD::seed(Pt, kMdPt);
  const MD tt = MD::seed(Tt, kMdTt);
  const MD pb = MD::seed(Pb, kMdPb);
  const MD pd = from_static ? MD::seed(P1, kMdP1) : pt;
  const MD n = gamma_of_Tt(Tt, X);
  const MD R = MD::seed(combaero::specific_gas_constant(X), kMdR);
  const MD rho0 = pt / (R * tt);
  const MD rho1 = rho0 * orifice_jobson::pow_dd(pd / pt, 1.0 / n);

  const bool forward = pb.v <= pd.v;
  const MD r = forward ? pb / pd : pd / pb;
  const MD rho_up = forward ? rho1 : pb / (R * tt);
  const double b4 = b * b * b * b;
  namespace ie = ::orifice::iso5167_expansibility;
  const double coeff = ie::c0 + ie::c4 * b4 + ie::c8 * b4 * b4;

  ModelEval out;
  out.Ci = Cd;
  out.eps = 1.0 - coeff * (1.0 - orifice_jobson::pow_dd(r, 1.0 / n));
  // ISO's incompressible Eq. (1) at unit C (linear in it), with the network's
  // 1 Pa regularisation of the drive.
  const auto [k1, dk_ddp, dk_drho] = orifice_mdot_and_jacobian(pd.v - pb.v, rho_up.v, 1.0, area, b);
  MD m_inc = MD::constant(Cd * k1);
  for (int i = 0; i < kMdCount; ++i) {
    m_inc.d[i] = Cd * (dk_ddp * (pd.d[i] - pb.d[i]) + dk_drho * rho_up.d[i]);
  }
  m_inc.d[kMdCd] += k1;
  out.m = out.eps * m_inc;
  out.in_range = r.v >= ie::min_pressure_ratio;
  return out;
}

template <typename Eval>
OrificeModelResult model_result(const Eval& eval, const std::vector<double>& Y_up, double Tt) {
  const std::vector<double> X =
      combaero::mass_to_mole(combaero::normalize_fractions(Y_up));
  const ModelEval e = eval(X);
  OrificeModelResult res;
  res.m_dot_calc = e.m.v;
  res.d_mdot_dP_total_up = e.m.d[kMdPt];
  res.d_mdot_dP_static_up = e.m.d[kMdP1];
  res.d_mdot_dP_static_down = e.m.d[kMdPb];
  res.d_mdot_dP_rho_ref = 0.0;
  res.d_mdot_dTt_up = e.m.d[kMdTt];
  res.d_mdot_dCd = e.m.d[kMdCd];
  res.Ci = e.Ci;
  res.contraction = e.C;
  res.expansibility = e.eps.v;
  res.choked = e.choked;
  res.in_range = e.in_range;
  // Composition, analytic (#964; a forward difference until then). It enters
  // both models through n = gamma(Tt) = cp / (cp - R_u) (molar cp = sum_k X_k
  // cp_k(Tt)) and R = R_u / MW (MW = sum_k X_k W_k) only, whose own partials
  // the DualN carries (kMdN, kMdR):
  //   dn/dX_j = -R_u cp_j(Tt) / (cp - R_u)^2,  dR/dX_j = -R W_j / MW,
  // the mole fractions unnormalised, chained to the mass fractions.
  const std::size_t n = X.size();
  std::vector<double> dm_dX(n, 0.0);
  const double c = combaero::cp(Tt, X);
  const double cv = c - combaero::thermo::R_GAS;
  const double R = combaero::specific_gas_constant(X);
  const double MW = combaero::mwmix(X);
  for (std::size_t j = 0; j < n; ++j) {
    const double dn = -combaero::thermo::R_GAS * combaero::cp_species(j, Tt) / (cv * cv);
    const double dR = -R * combaero::species_molar_mass(j) / MW;
    dm_dX[j] = e.m.d[kMdN] * dn + e.m.d[kMdR] * dR;
  }
  res.d_mdot_dY_up = combaero::mole_gradient_to_mass(X, dm_dX);
  return res;
}

}  // namespace

OrificeModelResult orifice_jobson_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double Tt_up,
    const std::vector<double>& Y_up, double P_static_down, double Cd, double area,
    double beta, double Cv) {
  (void)m_dot;
  if (!validation::is_finite_positive(Cd) || !validation::is_finite_positive(Cv) ||
      !validation::is_finite_positive(area)) {
    throw std::invalid_argument(
        "orifice_jobson_residuals_and_jacobian: Cd, Cv and area must be > 0");
  }
  if (!validation::is_finite(beta)) {
    throw std::invalid_argument("orifice_jobson_residuals_and_jacobian: beta must be finite");
  }
  // The kernel's no-flow branch is for finite states; NaN used to land in it.
  if (!validation::is_finite(P_total_up) || !validation::is_finite(P_static_up) ||
      !validation::is_finite(Tt_up) || !validation::is_finite(P_static_down)) {
    throw std::invalid_argument(
        "orifice_jobson_residuals_and_jacobian: pressures and Tt_up must be finite");
  }
  return model_result(
      [&](const std::vector<double>& X) {
        return jobson_eval(P_total_up, P_static_up, Tt_up, P_static_down, Cd, area, beta, Cv,
                           X);
      },
      Y_up, Tt_up);
}

OrificeModelResult orifice_iso5167_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double Tt_up,
    const std::vector<double>& Y_up, double P_static_down, double Cd, double area,
    double beta) {
  (void)m_dot;
  if (!validation::is_finite_positive(Cd) || !validation::is_finite_positive(area)) {
    throw std::invalid_argument("orifice_iso5167_residuals_and_jacobian: Cd and area must be > 0");
  }
  if (!validation::is_finite(beta)) {
    throw std::invalid_argument("orifice_iso5167_residuals_and_jacobian: beta must be finite");
  }
  if (!validation::is_finite(P_total_up) || !validation::is_finite(P_static_up) ||
      !validation::is_finite(Tt_up) || !validation::is_finite(P_static_down)) {
    throw std::invalid_argument(
        "orifice_iso5167_residuals_and_jacobian: pressures and Tt_up must be finite");
  }
  return model_result(
      [&](const std::vector<double>& X) {
        return iso5167_eval(P_total_up, P_static_up, Tt_up, P_static_down, Cd, area, beta, X);
      },
      Y_up, Tt_up);
}

SonicStatic sonic_static_pressure(double Tt, double G, const std::vector<double>& Y) {
  auto X_of = [](const std::vector<double>& Yv) {
    return combaero::mass_to_mole(combaero::normalize_fractions(Yv));
  };
  const auto st = combaero::sonic_static_state(Tt, G, X_of(Y));
  SonicStatic out{};
  out.P = st.P;
  out.T = st.T;
  out.dP_dTt = st.dP_dTt;
  out.dP_dG = st.dP_dG;
  out.dP_dY.assign(Y.size(), 0.0);
  const double eps_Y = 1e-6;
  for (std::size_t i = 0; i < Y.size(); ++i) {
    std::vector<double> Y_plus = Y;
    Y_plus[i] += eps_Y;
    out.dP_dY[i] = (combaero::sonic_static_state(Tt, G, X_of(Y_plus)).P - out.P) / eps_Y;
  }
  return out;
}

namespace {

// d P* / dX_j of duct_static_state's choked branch (#964), the mole fractions
// unnormalised (chain with mole_gradient_to_mass). T* solves the M = 1 energy
// balance F(T, X) = dH(T) - gamma R T / 2 = 0 (mass basis, dH = sum_k X_k
// [h_k(Tt) - h_k(T)] / MW), and ln P* = ln Pt - sum_k X_k [s_k(Tt) - s_k(T*)]
// / R_u on the isentrope through (Tt, Pt). With gamma = cp / (cp - R_u)
// (molar cp) and R = R_u / MW:
//   dF/dX_j = (dh_j - dH W_j) / MW - (T/2) (R dgamma/dX_j + gamma dR/dX_j),
//   dgamma/dX_j = -R_u cp_j(T) / (cp - R_u)^2,  dR/dX_j = -R W_j / MW,
//   dT*/dX_j = -(dF/dX_j) / (dF/dT),
//   dlnP*/dX_j = -(s_j(Tt) - s_j(T*)) / R_u + cp(T*) / (R_u T*) dT*/dX_j.
std::vector<double> choked_duct_dP_dX(double Tt, double T, double P,
                                      const std::vector<double>& X) {
  const std::size_t n = X.size();
  const double Ru = combaero::thermo::R_GAS;
  const double MWg = combaero::mwmix(X);  // g/mol
  const double MW = MWg / 1000.0;         // kg/mol
  const double R = Ru / MW;
  const double c = combaero::cp(T, X);  // molar
  const double cv = c - Ru;
  const double g = c / cv;
  const double dg_dT = -Ru * combaero::dcp_dT(T, X) / (cv * cv);
  const double cp_mass = c / MW;
  const double dF_dT = -cp_mass - 0.5 * R * (g + T * dg_dT);
  std::vector<double> dh(n, 0.0);
  double dH = 0.0;
  for (std::size_t j = 0; j < n; ++j) {
    dh[j] = combaero::h_species(j, Tt) - combaero::h_species(j, T);
    dH += X[j] * dh[j];
  }
  dH /= MW;
  std::vector<double> out(n, 0.0);
  for (std::size_t j = 0; j < n; ++j) {
    const double W = combaero::species_molar_mass(j) / 1000.0;
    const double dg = -Ru * combaero::cp_species(j, T) / (cv * cv);
    const double dR = -R * W / MW;
    const double dF = (dh[j] - dH * W) / MW - 0.5 * T * (R * dg + g * dR);
    const double dT = -dF / dF_dT;
    const double dlnP =
        -(combaero::s_species(j, Tt) - combaero::s_species(j, T)) / Ru + c / (Ru * T) * dT;
    out[j] = P * dlnP;
  }
  return out;
}

}  // namespace

ApproachStatic approach_static_pressure(double Pt, double Tt, double G,
                                        const std::vector<double>& Y,
                                        bool incompressible) {
  const std::vector<double> X = combaero::mass_to_mole(combaero::normalize_fractions(Y));
  ApproachStatic out{};
  // Composition moves R (and, compressible, cp). Analytic (#964; a forward
  // difference until then): dP_dX in unnormalised mole fractions, chained to
  // the mass fractions below.
  std::vector<double> dP_dX(X.size(), 0.0);
  if (incompressible) {
    const double R = combaero::specific_gas_constant(X);
    const double c = 2.0 * G * G * R * Tt;
    const double disc = Pt * Pt - c;
    out.choked = !(disc > 0.0);
    if (out.choked) {
      out.P = 0.5 * Pt;
      out.dP_dPt = 0.5;
      out.dP_dTt = 0.0;
      out.dP_dG = 0.0;
    } else {
      const double r = std::sqrt(disc);
      out.P = 0.5 * (Pt + r);
      out.dP_dPt = 0.5 * (1.0 + Pt / r);
      out.dP_dTt = -0.5 * (G * G * R) / r;
      out.dP_dG = -(G * R * Tt) / r;
      // dP/dR = -G^2 Tt / (2 r), dR/dX_j = -R W_j / MW.
      const double dP_dR = -0.5 * G * G * Tt / r;
      const double mw = combaero::mwmix(X);
      for (std::size_t j = 0; j < X.size(); ++j) {
        dP_dX[j] = dP_dR * (-R * combaero::species_molar_mass(j) / mw);
      }
    }
    out.dP_dY = combaero::mole_gradient_to_mass(X, dP_dX);
    return out;
  }
  const auto st = combaero::duct_static_state(Tt, Pt, G, X);
  out.P = st.P;
  out.dP_dPt = st.dP_dPt;
  out.dP_dTt = st.dP_dTt;
  out.dP_dG = st.dP_dG;
  out.choked = st.choked;
  const bool floored = !(Tt > combaero::thermo::T_FLOOR) || !(Pt > combaero::thermo::P_FLOOR);
  if (floored) {
    // Evaluated at the floors: no composition column, as no Tt or Pt one.
    out.dP_dY.assign(Y.size(), 0.0);
  } else if (st.choked) {
    out.dP_dY = combaero::mole_gradient_to_mass(X, choked_duct_dP_dX(Tt, st.T, st.P, X));
  } else {
    // Subsonic: the implicit-function gradient of the same two equations
    // duct_static_state solves (G_eff = G, so nothing is solved again).
    out.dP_dY = combaero::duct_htc_static_state_dY(Tt, Pt, G, X, st).dP_dY;
  }
  return out;
}

OrificeTapPressure orifice_tap_pressure_and_jacobian(double P_static_up,
                                                     double Pt_up,
                                                     double Pt_down,
                                                     double beta, double C,
                                                     bool compressible) {
  auto [k, dk_dC] = permanent_loss_ratio_and_dC(beta, C);
  const double dw = Pt_up - Pt_down;
  OrificeTapPressure out{};
  out.loss_ratio = k;
  out.p_tap = P_static_up - dw / k;
  out.d_dP_static_up = 1.0;
  out.d_dP_total_up = -1.0 / k;
  out.d_dP_total_down = 1.0 / k;
  out.d_dC = dw * dk_dC / (k * k);
  const double p_floor = ORIFICE_TAP_PRESSURE_FLOOR * P_static_up;
  out.floored = compressible && out.p_tap < p_floor;
  if (out.floored) {
    out.p_tap = p_floor;
    out.d_dP_static_up = ORIFICE_TAP_PRESSURE_FLOOR;
    out.d_dP_total_up = 0.0;
    out.d_dP_total_down = 0.0;
    out.d_dC = 0.0;
  }
  out.in_range = P_static_up > 0.0 &&
                 out.p_tap >= ORIFICE_RECOVERY_MIN_PRESSURE_RATIO * P_static_up;
  return out;
}

// -----------------------------------------------------------------------------
// 5. Tee Junction Elements (Bassett 2001)
// -----------------------------------------------------------------------------

namespace {
// Shared core computation for merging/branching tee wrappers.
// Returns residuals + analytical m_dot Jacobians given pre-computed K values.
struct TeeCore {
    double R_straight;
    double R_branch;
    double dRs_d_mcom;
    double dRs_d_mbranch;
    double dRb_d_mcom;
    double dRb_d_mbranch;
    double q_dyn;
    double q;
};

TeeCore tee_core(double m_dot_com, double m_dot_branch,
                 double dP0_straight, double dP0_branch,
                 double rho, double F_C,
                 double Ks, double dKs_dq,
                 double Kb, double dKb_dq) {
    // Regularised flow ratio: q = m_dot_branch/m_dot_com, safe at m_dot_com=0.
    // Uses q = m_b * m_c / (m_c^2 + eps^2) which gives 0 at m_c=0 and
    // correct sign for both positive and negative m_dot_com.
    constexpr double eps_m = 1e-10; // kg/s
    const double m_com_sq = m_dot_com * m_dot_com;
    const double m_com_sq_safe = m_com_sq + eps_m * eps_m;
    const double q = m_dot_branch * m_dot_com / m_com_sq_safe;

    // Regularised u_com^2 so residual is finite even at zero flow.
    constexpr double eps_v_sq = 1e-12; // (m/s)^2
    const double u_com_sq_safe =
        m_com_sq / (rho * rho * F_C * F_C) + eps_v_sq;
    const double q_dyn = 0.5 * rho * u_com_sq_safe;

    // Residuals
    const double R_straight = dP0_straight - Ks * q_dyn;
    const double R_branch = dP0_branch - Kb * q_dyn;

    // Analytical Jacobians for m_dot (psi/theta fixed; only q varies with flow)
    // dq/dm_dot_com = m_dot_branch * (eps^2 - m_com_sq) / m_com_sq_safe^2
    // dq/dm_dot_branch = m_dot_com / m_com_sq_safe
    const double dq_dmcom =
        m_dot_branch * (eps_m * eps_m - m_com_sq) / (m_com_sq_safe * m_com_sq_safe);
    const double dq_dmbranch = m_dot_com / m_com_sq_safe;

    // dq_dyn/dm_dot_com = m_dot_com / (rho * F_C^2)
    const double dqdyn_dmcom = m_dot_com / (rho * F_C * F_C);

    const double dRs_dmcom =
        -dKs_dq * dq_dmcom * q_dyn - Ks * dqdyn_dmcom;
    const double dRs_dmbranch = -dKs_dq * dq_dmbranch * q_dyn;
    const double dRb_dmcom =
        -dKb_dq * dq_dmcom * q_dyn - Kb * dqdyn_dmcom;
    const double dRb_dmbranch = -dKb_dq * dq_dmbranch * q_dyn;

    return {R_straight, R_branch, dRs_dmcom, dRs_dmbranch,
            dRb_dmcom, dRb_dmbranch, q_dyn, q};
}

// Finite-difference Jacobians for density-dependent variables (T, P_static, Y).
// Evaluates func(rho_perturbed) to get the FD derivative.
void tee_fd_density_jacobians(
    double m_dot_com, double m_dot_branch,
    double dP0_straight, double dP0_branch,
    double P_static, double T,
    const std::vector<double>& X,
    double theta, double psi, double F_C, double blend_k,
    bool is_merging,
    double& dRs_dP, double& dRb_dP,
    double& dRs_dT, double& dRb_dT,
    std::vector<double>& dRs_dY, std::vector<double>& dRb_dY,
    const std::vector<double>& Y) {

    // K depends only on q/psi/theta -- compute once, reuse across all rho evals.
    constexpr double eps_m = 1e-10;
    const double q = m_dot_branch * m_dot_com
                     / (m_dot_com * m_dot_com + eps_m * eps_m);
    double Ks, dKs, Kb, dKb;
    if (is_merging) {
        Ks  = merging_tee_K_straight(q, psi, theta, blend_k);
        dKs = merging_tee_dK_straight_dq(q, psi, theta, blend_k);
        Kb  = merging_tee_K_branch(q, psi, theta, blend_k);
        dKb = merging_tee_dK_branch_dq(q, psi, theta, blend_k);
    } else {
        Ks  = branching_tee_K_straight(q, psi, theta, blend_k);
        dKs = branching_tee_dK_straight_dq(q, psi, theta, blend_k);
        Kb  = branching_tee_K_branch(q, psi, theta, blend_k);
        dKb = branching_tee_dK_branch_dq(q, psi, theta, blend_k);
    }

    auto eval = [&](double rho) {
        auto c = tee_core(m_dot_com, m_dot_branch, dP0_straight, dP0_branch,
                          rho, F_C, Ks, dKs, Kb, dKb);
        return std::make_pair(c.R_straight, c.R_branch);
    };

    // Central FD w.r.t. P_static -- scale with pressure so the step is
    // never negligibly small at high-altitude or high-pressure conditions.
    const double eps_P = std::max(1.0, std::abs(P_static) * 1e-4);
    auto [Rs_pp, Rb_pp] = eval(std::get<0>(density_and_jacobians(T, P_static + eps_P, X)));
    auto [Rs_pm, Rb_pm] = eval(std::get<0>(density_and_jacobians(T, P_static - eps_P, X)));
    dRs_dP = (Rs_pp - Rs_pm) / (2.0 * eps_P);
    dRb_dP = (Rb_pp - Rb_pm) / (2.0 * eps_P);

    // Central FD w.r.t. T
    const double eps_T = 1e-3;
    auto [Rs_tp, Rb_tp] = eval(std::get<0>(density_and_jacobians(T + eps_T, P_static, X)));
    auto [Rs_tm, Rb_tm] = eval(std::get<0>(density_and_jacobians(T - eps_T, P_static, X)));
    dRs_dT = (Rs_tp - Rs_tm) / (2.0 * eps_T);
    dRb_dT = (Rb_tp - Rb_tm) / (2.0 * eps_T);

    // Central FD w.r.t. Y[i] (composition)
    const double eps_Y = 1e-6;
    dRs_dY.assign(Y.size(), 0.0);
    dRb_dY.assign(Y.size(), 0.0);
    for (std::size_t i = 0; i < Y.size(); ++i) {
        std::vector<double> Y_p = Y, Y_m = Y;
        Y_p[i] += eps_Y;
        Y_m[i] -= eps_Y;
        auto X_p = combaero::mass_to_mole(combaero::normalize_fractions(Y_p));
        auto X_m = combaero::mass_to_mole(combaero::normalize_fractions(Y_m));
        auto [Rs_yp, Rb_yp] = eval(std::get<0>(density_and_jacobians(T, P_static, X_p)));
        auto [Rs_ym, Rb_ym] = eval(std::get<0>(density_and_jacobians(T, P_static, X_m)));
        dRs_dY[i] = (Rs_yp - Rs_ym) / (2.0 * eps_Y);
        dRb_dY[i] = (Rb_yp - Rb_ym) / (2.0 * eps_Y);
    }
}
} // anonymous namespace

TeeJunctionResult merging_tee_residuals_and_jacobian(
    double m_dot_com, double m_dot_branch,
    double dP0_straight, double dP0_branch,
    double P_static_com, double T_com,
    const std::vector<double>& Y_com,
    double theta, double psi, double F_C,
    double blend_k) {

    const auto X_com = combaero::mass_to_mole(
        combaero::normalize_fractions(Y_com));
    const auto [rho, drho_dT, drho_dP] =
        density_and_jacobians(T_com, P_static_com, X_com);

    // Flow ratio (regularised) and blended K values
    constexpr double eps_m = 1e-10;
    const double m_com_sq_safe =
        m_dot_com * m_dot_com + eps_m * eps_m;
    const double q = m_dot_branch * m_dot_com / m_com_sq_safe;

    const double Ks = merging_tee_K_straight(q, psi, theta, blend_k);
    const double dKs = merging_tee_dK_straight_dq(q, psi, theta, blend_k);
    const double Kb = merging_tee_K_branch(q, psi, theta, blend_k);
    const double dKb = merging_tee_dK_branch_dq(q, psi, theta, blend_k);

    auto core = tee_core(m_dot_com, m_dot_branch, dP0_straight, dP0_branch,
                         rho, F_C, Ks, dKs, Kb, dKb);

    TeeJunctionResult res;
    res.R_straight = core.R_straight;
    res.R_branch = core.R_branch;
    res.dR_straight_d_mdot_com = core.dRs_d_mcom;
    res.dR_straight_d_mdot_branch = core.dRs_d_mbranch;
    res.dR_branch_d_mdot_com = core.dRb_d_mcom;
    res.dR_branch_d_mdot_branch = core.dRb_d_mbranch;

    tee_fd_density_jacobians(
        m_dot_com, m_dot_branch, dP0_straight, dP0_branch,
        P_static_com, T_com, X_com, theta, psi, F_C, blend_k, true,
        res.dR_straight_dP_static_com, res.dR_branch_dP_static_com,
        res.dR_straight_dT_com, res.dR_branch_dT_com,
        res.dR_straight_dY_com, res.dR_branch_dY_com, Y_com);

    res.K_straight = Ks;
    res.K_branch = Kb;
    res.q = q;
    res.blend_w = blend_weight(q, blend_k);
    res.topology_valid = tee_topology_valid(q);
    auto inp = tee_check_inputs(q, psi, theta);
    res.status = inp.valid() ? CorrelationValidity::VALID
                             : CorrelationValidity::EXTRAPOLATED;

    return res;
}

TeeJunctionResult branching_tee_residuals_and_jacobian(
    double m_dot_com, double m_dot_branch,
    double dP0_straight, double dP0_branch,
    double P_static_com, double T_com,
    const std::vector<double>& Y_com,
    double theta, double psi, double F_C,
    double blend_k) {

    const auto X_com = combaero::mass_to_mole(
        combaero::normalize_fractions(Y_com));
    const auto [rho, drho_dT, drho_dP] =
        density_and_jacobians(T_com, P_static_com, X_com);

    constexpr double eps_m = 1e-10;
    const double m_com_sq_safe =
        m_dot_com * m_dot_com + eps_m * eps_m;
    const double q = m_dot_branch * m_dot_com / m_com_sq_safe;

    const double Ks = branching_tee_K_straight(q, psi, theta, blend_k);
    const double dKs = branching_tee_dK_straight_dq(q, psi, theta, blend_k);
    const double Kb = branching_tee_K_branch(q, psi, theta, blend_k);
    const double dKb = branching_tee_dK_branch_dq(q, psi, theta, blend_k);

    auto core = tee_core(m_dot_com, m_dot_branch, dP0_straight, dP0_branch,
                         rho, F_C, Ks, dKs, Kb, dKb);

    TeeJunctionResult res;
    res.R_straight = core.R_straight;
    res.R_branch = core.R_branch;
    res.dR_straight_d_mdot_com = core.dRs_d_mcom;
    res.dR_straight_d_mdot_branch = core.dRs_d_mbranch;
    res.dR_branch_d_mdot_com = core.dRb_d_mcom;
    res.dR_branch_d_mdot_branch = core.dRb_d_mbranch;

    tee_fd_density_jacobians(
        m_dot_com, m_dot_branch, dP0_straight, dP0_branch,
        P_static_com, T_com, X_com, theta, psi, F_C, blend_k, false,
        res.dR_straight_dP_static_com, res.dR_branch_dP_static_com,
        res.dR_straight_dT_com, res.dR_branch_dT_com,
        res.dR_straight_dY_com, res.dR_branch_dY_com, Y_com);

    res.K_straight = Ks;
    res.K_branch = Kb;
    res.q = q;
    res.blend_w = blend_weight(q, blend_k);
    res.topology_valid = tee_topology_valid(q);
    auto inp = tee_check_inputs(q, psi, theta);
    res.status = inp.valid() ? CorrelationValidity::VALID
                             : CorrelationValidity::EXTRAPOLATED;

    return res;
}

// =============================================================================
// Unified0D compressible tee junction (Section 5-7, sources/junction/)
// =============================================================================

// R enters both tee kernels only as R T (u ~ m R T / (P A), q ~ m^2 R T /
// (P A^2), the merging datum's R T_dat from the suppliers' R T and gamma):
// d/dR at fixed T is (T / R) d/dT (#1231).
static void tee_gas_constant_partials(CompressibleTeeResult &res, const BranchInput &com,
                                      const BranchInput &str, const BranchInput &bra) {
  auto ratio = [](const BranchInput &b) { return b.R_gas > 0.0 ? b.T / b.R_gas : 0.0; };
  res.dR0_dRgas_com = ratio(com) * res.dR0_dT_com;
  res.dR0_dRgas_str = ratio(str) * res.dR0_dT_str;
  res.dR0_dRgas_bra = ratio(bra) * res.dR0_dT_bra;
  res.dR1_dRgas_com = ratio(com) * res.dR1_dT_com;
  res.dR1_dRgas_str = ratio(str) * res.dR1_dT_str;
  res.dR1_dRgas_bra = ratio(bra) * res.dR1_dT_bra;
}

CompressibleTeeResult compressible_branching_tee_rj(
    const BranchInput& com,
    const BranchInput& str,
    const BranchInput& bra)
{
    constexpr double eps_m = 1e-10;

    // Pseudodatum = common branch (single supplier for branching)
    const double m_com = com.m_dot;
    const double m_bra = bra.m_dot;
    const double m_str = str.m_dot;  // = m_com - m_bra, passed by caller
    const double Pc    = com.P_static;
    const double Tc    = com.T;
    const double Rc    = com.R_gas;
    const double Ac    = com.A;

    const double m_com_sq_safe = m_com * m_com + eps_m * eps_m;
    const double u_dat  = m_com * Rc * Tc / (Pc * Ac);
    const double q_ref  = 0.5 * (Pc / (Rc * Tc)) * u_dat * u_dat;
    // A_dat = Ac (derived: m_com / (rho_dat * u_dat) = m_com * Rc * Tc / (Pc * u_dat) = Ac)

    // Flow ratios (regularised, same pattern as tee_core)
    const double x_str = m_str * m_com / m_com_sq_safe * (Ac / str.A);
    const double x_bra = m_bra * m_com / m_com_sq_safe * (Ac / bra.A);

    // Effective inflow angles: phi_j = 0.75 * |theta_j - theta_dat|.
    // (The LaTeX formula 0.75*(pi-(theta_dat-theta_j)) applies the convention that
    // theta_dat is measured as the branch's outward direction; Python passes the
    // "same-axis = 0" convention so the equivalent form is the absolute difference.)
    const double phi_str = 0.75 * std::abs(str.theta - com.theta);
    const double phi_bra = 0.75 * std::abs(bra.theta - com.theta);

    // Both collector legs use one consistent blended turning-loss closure (Finding 5):
    //   R = Pt_com - Pt_leg - K_turn(theta)*q_leg - BETA_EXTRACT*x_leg*K_bc_leg*q_ref
    // where q_leg is the leg-local dynamic head (turning loss, vanishes at theta=0) and
    // K_bc_leg = 1 + x^2 - 2*x*cos(phi) is the incompressible Borda-Carnot kernel
    // (faithful to the FD prototype; phi is constant w.r.t. the unknowns, so its
    // x-derivative is the trivial 2x - 2cos(phi)). The straight leg has theta_str = 0,
    // so K_turn(0) = 0 and its closure collapses to the pure extraction term.
    const double BE = combaero::BETA_EXTRACT;

    // Straight collector (theta_str = 0 -> no turning loss).
    const double cos_phi_str = std::cos(phi_str);
    const double K_bc_str    = 1.0 + x_str * x_str - 2.0 * x_str * cos_phi_str;
    const double dKbc_str_dx = 2.0 * x_str - 2.0 * cos_phi_str;

    const double R_0 = com.Pt - str.Pt - BE * x_str * K_bc_str * q_ref;

    // Branch collector (theta_bra != 0 -> turning loss present).
    const double cos_phi_bra = std::cos(phi_bra);
    const double K_bc        = 1.0 + x_bra * x_bra - 2.0 * x_bra * cos_phi_bra;
    const double dKbc_dx     = 2.0 * x_bra - 2.0 * cos_phi_bra;
    const double K_turn      = combaero::K_turn_div(bra.theta);

    const double Pb    = bra.P_static;
    const double Tb    = bra.T;
    const double Rb    = bra.R_gas;
    const double Ab    = bra.A;
    const double q_bra = 0.5 * Rb * Tb * m_bra * m_bra / (Pb * Ab * Ab);

    const double R_1 = com.Pt - bra.Pt - K_turn * q_bra - BE * x_bra * K_bc * q_ref;

    // Pseudodatum derivatives (q_ref = 0.5 * m_com^2 * Rc * Tc / (Pc * Ac^2)).
    const double dq_dPc    = -q_ref / Pc;
    const double dq_dTc    = q_ref / Tc;
    const double dq_dm_com = m_com * Rc * Tc / (Pc * Ac * Ac);

    // Flow-ratio derivatives (regularised).
    // x_str = m_str * m_com / m_com_sq_safe * psi_str
    // In Python: m_str = m_com - m_branch, so dm_str/dm_com = 1, dm_str/dm_branch = -1.
    // Total dx_str/dm_com includes both the partial w.r.t. m_com and the chain through m_str.
    const double dx_str_dm_com_partial = m_str * (eps_m*eps_m - m_com*m_com) / (m_com_sq_safe*m_com_sq_safe) * (Ac / str.A);
    const double dx_str_dm_str = m_com / m_com_sq_safe * (Ac / str.A);  // dm_str/dm_com = 1
    const double dx_str_dm_com = dx_str_dm_com_partial + dx_str_dm_str;
    const double dx_str_dm_branch = -dx_str_dm_str;  // dm_str/dm_branch = -1
    const double dx_bra_dm_com    = m_bra * (eps_m*eps_m - m_com*m_com) / (m_com_sq_safe*m_com_sq_safe) * (Ac / bra.A);
    const double dx_bra_dm_branch = m_com / m_com_sq_safe * (Ac / bra.A);   // dm_bra/dm_branch = +1

    CompressibleTeeResult res{};
    res.R_0 = R_0;
    res.R_1 = R_1;

    // R_0 Jacobians: R_0 = Pt_com - Pt_str - BETA_EXTRACT*x_str*K_bc_str*q_ref.
    // Pure extraction term (straight leg, no turning loss); phi_str is constant w.r.t.
    // the unknowns. d(x_str*K_bc_str)/dx_str = K_bc_str + x_str*dKbc_str_dx.
    const double d_xKbc_str_dx = K_bc_str + x_str * dKbc_str_dx;
    res.dR0_dPt_com    = 1.0;
    res.dR0_dPt_str    = -1.0;
    res.dR0_dP_com     = -BE * x_str * K_bc_str * dq_dPc;
    res.dR0_dT_com     = -BE * x_str * K_bc_str * dq_dTc;
    res.dR0_dmdot_com  = -BE * (q_ref * d_xKbc_str_dx * dx_str_dm_com + x_str * K_bc_str * dq_dm_com);
    res.dR0_dmdot_branch = -BE * q_ref * d_xKbc_str_dx * dx_str_dm_branch;
    // R_0 has no bra node dependence for branching.
    res.dR0_dP_str = 0.0;
    res.dR0_dT_str = 0.0;
    res.dR0_dP_bra = 0.0;

    // R_1 Jacobians: R_1 = Pt_com - Pt_bra - K_turn*q_bra - BETA_EXTRACT*x_bra*K_bc*q_ref.
    // K_turn and phi_bra (hence K_bc's cos term) are constant w.r.t. the unknowns.
    // Extraction term E = BETA_EXTRACT*x_bra*K_bc*q_ref; d(x_bra*K_bc)/dx_bra = K_bc + x_bra*dKbc_dx.
    const double d_xKbc_dx      = K_bc + x_bra * dKbc_dx;
    const double dq_bra_dm_bra  = Rb * Tb * m_bra / (Pb * Ab * Ab);  // d(q_bra)/d(m_branch)

    res.dR1_dPt_com    = 1.0;
    res.dR1_dPt_bra    = -1.0;
    // Branch-local q_bra (q_bra ~ 1/P_bra, ~ T_bra); turning-loss term only.
    res.dR1_dP_bra     = K_turn * q_bra / Pb;
    res.dR1_dT_bra     = -K_turn * q_bra / Tb;
    // Common datum enters only through q_ref in the extraction term (K_bc is incompressible).
    res.dR1_dP_com     = -BE * x_bra * K_bc * dq_dPc;
    res.dR1_dT_com     = -BE * x_bra * K_bc * dq_dTc;
    // Flows: q_bra depends on m_branch only; x_bra and q_ref depend on m_com; x_bra on m_branch.
    res.dR1_dmdot_com  = -BE * (q_ref * d_xKbc_dx * dx_bra_dm_com + x_bra * K_bc * dq_dm_com);
    res.dR1_dmdot_branch = -K_turn * dq_bra_dm_bra
                           - BE * q_ref * d_xKbc_dx * dx_bra_dm_branch;
    // R_1 has no straight-node dependence for branching.
    res.dR1_dP_str = 0.0;
    res.dR1_dT_str = 0.0;
    res.dR1_dPt_str = 0.0;

    // gamma_eff does not enter the branching closure.
    tee_gas_constant_partials(res, com, str, bra);
    return res;
}

CompressibleTeeResult compressible_merging_tee_rj(
    const BranchInput& com,
    const BranchInput& str,
    const BranchInput& bra)
{
    using combaero::K_dat_j_closed;
    using combaero::cstep_deriv;

    const double eps_m = 1e-12;

    // Supplier branches: str and bra. Collector: com.
    const double m_str = str.m_dot;  // positive (supplier)
    const double m_bra = bra.m_dot;  // positive (supplier)
    const double m_dat = m_str + m_bra;
    const double m_dat_safe = (m_dat < eps_m) ? eps_m : m_dat;

    const double Ps = str.P_static, Ts = str.T, Rs = str.R_gas, gs = str.gamma_eff, As = str.A;
    const double Pb = bra.P_static, Tb = bra.T, Rb = bra.R_gas, gb = bra.gamma_eff, Ab = bra.A;

    // No net supply (a solver probe: every flow zero, or the suppliers' sum
    // reversed): the pseudodatum is stagnant, p0_dat = p_dat and q_ref = 0 --
    // the m_dat -> 0+ limit of the general case. Dividing the suppliers'
    // momentum and enthalpy by the eps floor instead gave a negative datum
    // temperature and NaN (#481).
    if (m_dat < eps_m) {
        CompressibleTeeResult stagnant{};
        stagnant.R_0 = Ps - Pb;
        stagnant.R_1 = Ps - com.Pt;
        stagnant.dR0_dP_str = 1.0;
        stagnant.dR0_dP_bra = -1.0;
        stagnant.dR1_dP_str = 1.0;
        stagnant.dR1_dPt_com = -1.0;
        return stagnant;
    }

    const double cp_s = gs * Rs / (gs - 1.0);
    const double cp_b = gb * Rb / (gb - 1.0);
    const double cp_dat = cp_s;  // reference = straight supplier

    // Per-branch velocities and stagnation enthalpies
    const double u_str = m_str * Rs * Ts / (Ps * As);
    const double u_bra = m_bra * Rb * Tb / (Pb * Ab);
    const double h0_str = cp_s * Ts + 0.5 * u_str * u_str;
    const double h0_bra = cp_b * Tb + 0.5 * u_bra * u_bra;

    // Pseudodatum (mass-weighted average of suppliers)
    const double u_dat   = (m_str * u_str + m_bra * u_bra) / m_dat_safe;
    const double h0_dat  = (m_str * h0_str + m_bra * h0_bra) / m_dat_safe;
    const double T_dat   = (h0_dat - 0.5 * u_dat * u_dat) / cp_dat;
    const double p_dat   = Ps;  // reference = straight supplier static pressure
    const double g_dat   = gs;
    const double rho_dat = p_dat / (Rs * T_dat);
    const double M2_dat  = (u_dat > eps_m) ? u_dat * u_dat / (g_dat * Rs * T_dat) : 0.0;
    const double q_ref   = 0.5 * rho_dat * u_dat * u_dat;
    const double A_dat   = (u_dat > eps_m) ? m_dat_safe / (rho_dat * u_dat) : 1.0;

    // Pseudodatum angle (mass-weighted)
    const double S_S = m_str * std::sin(str.theta) + m_bra * std::sin(bra.theta);
    const double C_S = m_str * std::cos(str.theta) + m_bra * std::cos(bra.theta);
    const double SC_sq = S_S * S_S + C_S * C_S;
    const double theta_dat = std::atan2(S_S, C_S);

    // K_com: loss coefficient from pseudodatum to collector (common) branch.
    // x_com = A_dat / A_com  (m_com = m_dat so flow ratio = 1, area ratio gives x).
    const double x_com_v   = (u_dat > eps_m) ? A_dat / com.A : 0.0;
    const double phi_com_v = 0.75 * std::abs(com.theta - theta_dat);
    // dphi_com / d(theta_dat): sign depends on relative direction
    const double dphi_com_dtheta = (com.theta >= theta_dat) ? -0.75 : 0.75;
    const double K_com_v   = K_dat_j_closed(x_com_v, phi_com_v, M2_dat);

    // Stagnation pressure of pseudodatum: p0_dat = p_dat * Phi_dat^(gamma/(gamma-1))
    const double Phi_dat    = 1.0 + 0.5 * (g_dat - 1.0) * M2_dat;
    const double p0_dat_val = (u_dat > eps_m)
        ? p_dat * std::pow(Phi_dat, g_dat / (g_dat - 1.0)) : p_dat;

    // Residuals (LaTeX Sec 6, merging 2-supplier / 1-collector):
    //   R_sp      = p_str - p_bra              (equal static pressure at merge point)
    //   R_K,com   = p0_dat - p0_com - K * q    (stagnation loss from dat to collector)
    const double R_0 = Ps - Pb;
    const double R_1 = p0_dat_val - com.Pt - K_com_v * q_ref;

    // Complex-step K derivatives for K_com
    auto K_fn = [](auto x, auto phi, auto M2) { return K_dat_j_closed(x, phi, M2); };
    const double dKc_dx   = cstep_deriv(K_fn, 0, x_com_v, phi_com_v, M2_dat);
    const double dKc_dphi = cstep_deriv(K_fn, 1, x_com_v, phi_com_v, M2_dat);
    const double dKc_dM2  = cstep_deriv(K_fn, 2, x_com_v, phi_com_v, M2_dat);

    // Helper: compute d(-K_com * q_ref)/dv from pseudodatum chain.
    // m_j = m_dat, A_j = com.A, dm_j_dv = dm_d (since m_j = m_dat).
    auto dNeg_Kq = [&](double du_d, double dh0_d,
                        double dm_d, double dp_d, double dtheta_d) -> double
    {
        const double dT_d = (dh0_d - u_dat * du_d) / cp_dat;
        const double drho_d = rho_dat * (dp_d / p_dat - dT_d / T_dat);
        const double dM2_d = (u_dat > eps_m)
            ? M2_dat * (2.0 * du_d / u_dat - dT_d / T_dat) : 0.0;
        const double dq_d = (u_dat > eps_m)
            ? q_ref * (drho_d / rho_dat + 2.0 * du_d / u_dat) : 0.0;
        const double dA_d = (u_dat > eps_m)
            ? A_dat * (dm_d / m_dat_safe - drho_d / rho_dat - du_d / u_dat) : 0.0;
        // x_com = A_dat / A_com, so dx_com = dA_d / A_com
        const double dx_j  = dA_d / com.A;
        const double dphi_j = dphi_com_dtheta * dtheta_d;
        const double dK_j  = dKc_dx * dx_j + dKc_dphi * dphi_j + dKc_dM2 * dM2_d;
        return -(dK_j * q_ref + K_com_v * dq_d);
    };

    // Helper: d(p0_dat)/dv via pseudodatum chain.
    auto dp0dat_dv = [&](double du_d, double dh0_d, double dp_d) -> double {
        if (u_dat < eps_m) return dp_d;
        const double dT_d  = (dh0_d - u_dat * du_d) / cp_dat;
        const double dM2_d = M2_dat * (2.0 * du_d / u_dat - dT_d / T_dat);
        return (p0_dat_val / p_dat) * dp_d
             + p0_dat_val * g_dat / 2.0 / Phi_dat * dM2_d;
    };

    CompressibleTeeResult res{};
    res.R_0 = R_0;
    res.R_1 = R_1;

    // R_0 = Ps - Pb: trivial Jacobian (only static pressures of suppliers)
    res.dR0_dPt_str = 0.0; res.dR0_dPt_com = 0.0; res.dR0_dPt_bra = 0.0;
    res.dR0_dP_str  = 1.0; res.dR0_dP_com  = 0.0; res.dR0_dP_bra  = -1.0;
    res.dR0_dT_str  = 0.0; res.dR0_dT_com  = 0.0; res.dR0_dT_bra  = 0.0;
    res.dR0_dmdot_com = 0.0; res.dR0_dmdot_branch = 0.0;

    // R_1 = p0_dat - Pt_com - K_com*q_ref
    // Direct term: dR1/dPt_com = -1
    res.dR1_dPt_com = -1.0;
    res.dR1_dPt_str = 0.0; res.dR1_dPt_bra = 0.0;
    res.dR1_dP_com  = 0.0; res.dR1_dT_com  = 0.0;

    // dR1/dP_str: u_str depends on Ps; p_dat = Ps so dp_d = 1
    {
        const double du_s = -u_str / Ps;
        const double dh0_s = u_str * du_s;
        const double du_d = m_str / m_dat_safe * du_s;
        const double dh0_d = m_str / m_dat_safe * dh0_s;
        res.dR1_dP_str = dp0dat_dv(du_d, dh0_d, 1.0)
                       + dNeg_Kq(du_d, dh0_d, 0.0, 1.0, 0.0);
    }
    // dR1/dT_str
    {
        const double du_s = u_str / Ts;
        const double dh0_s = cp_s + u_str * u_str / Ts;
        const double du_d = m_str / m_dat_safe * du_s;
        const double dh0_d = m_str / m_dat_safe * dh0_s;
        res.dR1_dT_str = dp0dat_dv(du_d, dh0_d, 0.0)
                       + dNeg_Kq(du_d, dh0_d, 0.0, 0.0, 0.0);
    }
    // dR1/dP_bra: p_dat = Ps so dp_d = 0 for bra
    {
        const double du_b = -u_bra / Pb;
        const double dh0_b = u_bra * du_b;
        const double du_d = m_bra / m_dat_safe * du_b;
        const double dh0_d = m_bra / m_dat_safe * dh0_b;
        res.dR1_dP_bra = dp0dat_dv(du_d, dh0_d, 0.0)
                       + dNeg_Kq(du_d, dh0_d, 0.0, 0.0, 0.0);
    }
    // dR1/dT_bra
    {
        const double du_b = u_bra / Tb;
        const double dh0_b = cp_b + u_bra * u_bra / Tb;
        const double du_d = m_bra / m_dat_safe * du_b;
        const double dh0_d = m_bra / m_dat_safe * dh0_b;
        res.dR1_dT_bra = dp0dat_dv(du_d, dh0_d, 0.0)
                       + dNeg_Kq(du_d, dh0_d, 0.0, 0.0, 0.0);
    }
    // dR1/dm_dot_com: m_str = m_com - m_bra, dm_str/dm_com = 1, dm_dat/dm_com = 1
    {
        const double du_d = (2.0 * u_str - u_dat) / m_dat_safe;
        const double dh0_d = (u_str * u_str + h0_str - h0_dat) / m_dat_safe;
        const double dtheta_d = (SC_sq > 1e-30)
            ? (C_S * std::sin(str.theta) - S_S * std::cos(str.theta)) / SC_sq : 0.0;
        res.dR1_dmdot_com = dp0dat_dv(du_d, dh0_d, 0.0)
                          + dNeg_Kq(du_d, dh0_d, 1.0, 0.0, dtheta_d);
    }
    // dR1/dm_dot_branch: dm_str/dm_branch = -1, dm_bra/dm_branch = +1, dm_dat = 0
    {
        const double du_d = (-2.0 * u_str + 2.0 * u_bra) / m_dat_safe;
        const double dh0_d = (-h0_str - u_str*u_str + h0_bra + u_bra*u_bra) / m_dat_safe;
        const double dt_s = (SC_sq > 1e-30)
            ? (C_S * std::sin(str.theta) - S_S * std::cos(str.theta)) / SC_sq : 0.0;
        const double dt_b = (SC_sq > 1e-30)
            ? (C_S * std::sin(bra.theta) - S_S * std::cos(bra.theta)) / SC_sq : 0.0;
        const double dtheta_d = -dt_s + dt_b;
        res.dR1_dmdot_branch = dp0dat_dv(du_d, dh0_d, 0.0)
                             + dNeg_Kq(du_d, dh0_d, 0.0, 0.0, dtheta_d);
    }
    // dR1/dgamma_bra: gamma_bra enters only cp_bra, so h0_bra (#1231).
    {
        const double dcp_b = -Rb / ((gb - 1.0) * (gb - 1.0));
        const double dh0_d = m_bra / m_dat_safe * dcp_b * Tb;
        res.dR1_dgamma_bra = dp0dat_dv(0.0, dh0_d, 0.0)
                           + dNeg_Kq(0.0, dh0_d, 0.0, 0.0, 0.0);
    }
    // dR1/dgamma_str: cp_str (h0_str and the datum's cp, cp_dat = cp_str)
    // and the datum's own gamma (M2_dat, Phi_dat and its exponent).
    {
        const double dcp_s = -Rs / ((gs - 1.0) * (gs - 1.0));
        // T_dat = (h0_dat - u_dat^2 / 2) / cp_dat: the helpers' dT_d is
        // (dh0_d - u_dat du_d) / cp_dat, so cp_dat's own change enters as
        // an equivalent -T_dat dcp_s in dh0_d.
        const double dh0_d = m_str / m_dat_safe * dcp_s * Ts - T_dat * dcp_s;
        double d = dp0dat_dv(0.0, dh0_d, 0.0) + dNeg_Kq(0.0, dh0_d, 0.0, 0.0, 0.0);
        if (u_dat > eps_m) {
            const double gm1 = g_dat - 1.0;
            const double dM2_g = -M2_dat / g_dat;  // M2 = u^2 / (g R T_dat)
            const double dPhi = 0.5 * M2_dat + 0.5 * gm1 * dM2_g;
            d += p0_dat_val * (-std::log(Phi_dat) / (gm1 * gm1) + g_dat / gm1 * dPhi / Phi_dat);
            d -= dKc_dM2 * dM2_g * q_ref;
        }
        res.dR1_dgamma_str = d;
    }
    tee_gas_constant_partials(res, com, str, bra);
    return res;
}

} // namespace solver
} // namespace combaero
