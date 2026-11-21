#include "../include/heat_transfer.h"
#include "../include/cooling_correlations.h" // rib_*, dimple_*, pin_fin_*, impingement_*
#include "../include/composition.h"   // mwmix, mole_gradient_to_mass
#include "../include/correlation_status.h"
#include "../include/friction.h"       // friction_petukhov, friction_colebrook
#include "../include/math_constants.h" // M_PI cross-platform
#include "../include/solver_interface.h" // density_and_jacobians, viscosity_and_jacobians
#include "../include/stagnation.h"     // T_adiabatic_wall
#include "../include/state.h"
#include "../include/thermo.h"    // density, cp_mass, speed_of_sound
#include "../include/transport.h" // viscosity, thermal_conductivity, prandtl
#include "../include/input_validation.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>



namespace combaero {

// -------------------------------------------------------------
// Internal Flow Correlations
// -------------------------------------------------------------

double nusselt_dittus_boelter(double Re, double Pr, bool heating,
                              combaero::CorrelationStatus *status) {
  bool extrapolated = false;
  if (Re < 10000.0) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_dittus_boelter: Re = " + std::to_string(Re) +
                     " is below validated range [10000, inf]."
                     " Extrapolating power-law correlation; check results.");
    }
  }
  if (Pr < 0.6 || Pr > 160.0) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_dittus_boelter: Pr = " + std::to_string(Pr) +
                     " is outside validated range [0.6, 160]."
                     " Extrapolating; check results.");
    }
  }
  if (status) {
    *status = extrapolated ? combaero::CorrelationStatus::Extrapolated
                           : combaero::CorrelationStatus::Valid;
  }
  double n = heating ? dittus_boelter::coeff_prandtl_heating_exp
                     : dittus_boelter::coeff_prandtl_cooling_exp;
  return dittus_boelter::coeff_outer *
         std::pow(Re, dittus_boelter::coeff_reynolds_exp) * std::pow(Pr, n);
}

double nusselt_gnielinski(double Re, double Pr, double f,
                          combaero::CorrelationStatus *status) {
  if (!validation::is_finite_positive(f)) {
    throw std::invalid_argument(
        "nusselt_gnielinski: friction factor must be positive");
  }

  bool extrapolated = false;

  if (Re > 5.0e6) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_gnielinski: Re = " + std::to_string(Re) +
                     " is above validated range [2300, 5e6]."
                     " Extrapolating; check results.");
    }
  }
  if (Pr < 0.5 || Pr > 2000.0) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_gnielinski: Pr = " + std::to_string(Pr) +
                     " is outside validated range [0.5, 2000]."
                     " Extrapolating; check results.");
    }
  }
  if (status) {
    *status = extrapolated ? combaero::CorrelationStatus::Extrapolated
                           : combaero::CorrelationStatus::Valid;
  }

  double f8 = f / 8.0;
  double sqrt_f8 = std::sqrt(f8);
  double Pr_23 = std::pow(Pr, 2.0 / 3.0);
  double denom = 1.0 + gnielinski::coeff_prandtl_factor * sqrt_f8 * (Pr_23 - 1.0);

  if (Re >= gnielinski::blend_re_turbulent) {
    return f8 * (Re - gnielinski::coeff_re_offset) * Pr / denom;
  }

  // Below Re=2300 the raw Gnielinski formula goes negative at Re<1000.
  // Use a C1-smooth cubic Hermite blend between Re=1000 (laminar Nu) and
  // Re=2300 (Gnielinski value + slope), so the Jacobian is continuous.
  if (status && *status == combaero::CorrelationStatus::Valid) {
    *status = combaero::CorrelationStatus::Extrapolated;
  }
  if (!status) {
    combaero::warn("nusselt_gnielinski: Re = " + std::to_string(Re) +
                   " is below validated range [2300, 5e6]."
                   " Using C1 Hermite blend to laminar Nu; check results.");
  }

  // Anchor points for the Hermite cubic
  constexpr double Re0 = gnielinski::blend_re_laminar;   // lower anchor: laminar Nu
  constexpr double Re1 = gnielinski::blend_re_turbulent; // upper anchor: Gnielinski
  const double Nu0 = NU_LAMINAR_CONST_T; // 3.66
  const double dNu0 = 0.0;               // laminar: flat w.r.t. Re
  const double Nu1 = f8 * (Re1 - gnielinski::coeff_re_offset) * Pr / denom;
  const double dNu1 = f8 * Pr / denom; // d(Nu)/d(Re) at Re1

  // Normalised coordinate t in [0, 1]
  double h = Re1 - Re0;
  double t = (Re - Re0) / h;
  t = std::max(0.0, std::min(1.0, t)); // clamp outside [Re0, Re1]

  // Cubic Hermite basis
  double t2 = t * t;
  double t3 = t2 * t;
  double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
  double h10 = t3 - 2.0 * t2 + t;
  double h01 = -2.0 * t3 + 3.0 * t2;
  double h11 = t3 - t2;

  return h00 * Nu0 + h10 * h * dNu0 + h01 * Nu1 + h11 * h * dNu1;
}

double nusselt_gnielinski(double Re, double Pr,
                          combaero::CorrelationStatus *status) {
  return nusselt_gnielinski(Re, Pr, friction_petukhov_clamped(Re), status);
}

NuAndDerivative nusselt_gnielinski_with_derivative(double Re, double Pr,
                                                   double f, double df_dRe) {
  using namespace gnielinski;
  if (Re <= blend_re_laminar) {
    return {NU_LAMINAR_CONST_T, 0.0};
  }
  const double f8 = f / 8.0;
  const double sqrt_f8 = std::sqrt(f8);
  const double pr_term = std::pow(Pr, 2.0 / 3.0) - 1.0;
  const double denom = 1.0 + coeff_prandtl_factor * sqrt_f8 * pr_term;
  // q = f8 / denom carries all the friction dependence; dq/df for the
  // chain rule through a Re-dependent f.
  const double q = f8 / denom;
  const double ddenom_df = coeff_prandtl_factor * pr_term / (16.0 * sqrt_f8);
  const double dq_df = (denom / 8.0 - f8 * ddenom_df) / (denom * denom);

  if (Re >= blend_re_turbulent) {
    const double nu = q * (Re - coeff_re_offset) * Pr;
    const double dnu = q * Pr + (Re - coeff_re_offset) * Pr * dq_df * df_dRe;
    return {nu, dnu};
  }
  // C1 Hermite blend to laminar, as in nusselt_gnielinski; the turbulent
  // anchors depend on f, so a Re-dependent f moves them too.
  const double h = blend_re_turbulent - blend_re_laminar;
  const double nu1 = q * (blend_re_turbulent - coeff_re_offset) * Pr;
  const double dnu1 = q * Pr;
  const double dnu1_df = (blend_re_turbulent - coeff_re_offset) * Pr * dq_df;
  const double ddnu1_df = Pr * dq_df;
  const double t = (Re - blend_re_laminar) / h;
  const double t2 = t * t;
  const double t3 = t2 * t;
  const double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
  const double h01 = -2.0 * t3 + 3.0 * t2;
  const double h11 = t3 - t2;
  const double nu = h00 * NU_LAMINAR_CONST_T + h01 * nu1 + h11 * h * dnu1;
  const double dnu_dt = (6.0 * t2 - 6.0 * t) * NU_LAMINAR_CONST_T +
                        (-6.0 * t2 + 6.0 * t) * nu1 +
                        (3.0 * t2 - 2.0 * t) * h * dnu1;
  const double dnu = dnu_dt / h +
                     (h01 * dnu1_df + h11 * h * ddnu1_df) * df_dRe;
  return {nu, dnu};
}

NuAndDerivative nusselt_gnielinski_smooth_with_derivative(double Re,
                                                          double Pr) {
  return nusselt_gnielinski_with_derivative(
      Re, Pr, friction_petukhov_clamped(Re), friction_petukhov_clamped_dRe(Re));
}

namespace {

// C1 smoothstep weight over [lo, hi] in Re: 0 below, 1 above, zero slope at
// both ends.
struct Weight {
  double w;
  double dw;
};

Weight smoothstep(double Re, double lo, double hi) {
  if (Re <= lo) {
    return {0.0, 0.0};
  }
  if (Re >= hi) {
    return {1.0, 0.0};
  }
  const double s = (Re - lo) / (hi - lo);
  return {s * s * (3.0 - 2.0 * s), 6.0 * s * (1.0 - s) / (hi - lo)};
}

}  // namespace

FrictionAndDerivative friction_turbulent_and_derivative(double Re, double e_D) {
  const FrictionAndDerivative smooth{friction_petukhov_clamped(Re),
                                     friction_petukhov_clamped_dRe(Re)};
  if (!(e_D > 0.0)) {
    return smooth;
  }
  const Weight w = smoothstep(Re, channel_transition::turbulent_start,
                              channel_transition::rough_full);
  if (w.w == 0.0) {
    return smooth;
  }
  const FrictionAndDerivative rough{friction_colebrook(Re, e_D),
                                    friction_colebrook_dRe(Re, e_D)};
  if (w.w == 1.0) {
    return rough;
  }
  return {w.w * rough.f + (1.0 - w.w) * smooth.f,
          w.dw * (rough.f - smooth.f) + w.w * rough.df_dRe +
              (1.0 - w.w) * smooth.df_dRe};
}

FrictionAndDerivative friction_channel_and_derivative(double Re, double e_D) {
  if (!(Re > 0.0)) {
    return {0.0, 0.0};
  }
  const FrictionAndDerivative lam{64.0 / Re, -64.0 / (Re * Re)};
  const Weight w = smoothstep(Re, channel_transition::laminar_end,
                              channel_transition::turbulent_start);
  if (w.w == 0.0) {
    return lam;
  }
  const FrictionAndDerivative turb = friction_turbulent_and_derivative(Re, e_D);
  if (w.w == 1.0) {
    return turb;
  }
  return {w.w * turb.f + (1.0 - w.w) * lam.f,
          w.dw * (turb.f - lam.f) + w.w * turb.df_dRe + (1.0 - w.w) * lam.df_dRe};
}

NuAndDerivative nusselt_channel_gnielinski_and_derivative(
    double Re, double Pr, double f_turb, double df_turb_dRe, double Nu_laminar) {
  const Weight w = smoothstep(Re, channel_transition::laminar_end,
                              channel_transition::turbulent_start);
  if (w.w == 0.0) {
    return {Nu_laminar, 0.0};
  }
  const NuAndDerivative turb =
      nusselt_gnielinski_with_derivative(Re, Pr, f_turb, df_turb_dRe);
  if (w.w == 1.0) {
    return turb;
  }
  return {w.w * turb.Nu + (1.0 - w.w) * Nu_laminar,
          w.dw * (turb.Nu - Nu_laminar) + w.w * turb.dNu_dRe};
}

double nusselt_sieder_tate(double Re, double Pr, double mu_ratio,
                           combaero::CorrelationStatus *status) {
  if (!validation::is_finite_positive(mu_ratio)) {
    throw std::invalid_argument(
        "nusselt_sieder_tate: mu_ratio must be positive");
  }
  bool extrapolated = false;
  if (Re < 10000.0) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_sieder_tate: Re = " + std::to_string(Re) +
                     " is below validated range [10000, inf]."
                     " Extrapolating power-law correlation; check results.");
    }
  }
  if (Pr < 0.7 || Pr > 16700.0) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_sieder_tate: Pr = " + std::to_string(Pr) +
                     " is outside validated range [0.7, 16700]."
                     " Extrapolating; check results.");
    }
  }
  if (status) {
    *status = extrapolated ? combaero::CorrelationStatus::Extrapolated
                           : combaero::CorrelationStatus::Valid;
  }
  return sieder_tate::coeff_outer
         * std::pow(Re, sieder_tate::coeff_reynolds_exp)
         * std::pow(Pr, sieder_tate::coeff_prandtl_exp)
         * std::pow(mu_ratio, sieder_tate::coeff_viscosity_exp);
}

double nusselt_petukhov(double Re, double Pr, double f,
                        combaero::CorrelationStatus *status) {
  if (!validation::is_finite_positive(f)) {
    throw std::invalid_argument(
        "nusselt_petukhov: friction factor must be positive");
  }
  bool extrapolated = false;
  if (Re < 1.0e4 || Re > 5.0e6) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_petukhov: Re = " + std::to_string(Re) +
                     " is outside validated range [1e4, 5e6]."
                     " Extrapolating; check results.");
    }
  }
  if (Pr < 0.5 || Pr > 2000.0) {
    extrapolated = true;
    if (!status) {
      combaero::warn("nusselt_petukhov: Pr = " + std::to_string(Pr) +
                     " is outside validated range [0.5, 2000]."
                     " Extrapolating; check results.");
    }
  }
  if (status) {
    *status = extrapolated ? combaero::CorrelationStatus::Extrapolated
                           : combaero::CorrelationStatus::Valid;
  }
  double f8 = f / 8.0;
  double sqrt_f8 = std::sqrt(f8);
  double Pr_23 = std::pow(Pr, 2.0 / 3.0);
  return f8 * Re * Pr / (petukhov_ht::coeff_offset + petukhov_ht::coeff_prandtl_factor * sqrt_f8 * (Pr_23 - 1.0));
}

double nusselt_petukhov(double Re, double Pr,
                        combaero::CorrelationStatus *status) {
  double f =
      friction_petukhov_clamped(Re); // avoid log(0) in Petukhov
  return nusselt_petukhov(Re, Pr, f, status);
}

// -------------------------------------------------------------
// Helper Functions
// -------------------------------------------------------------

double htc_from_nusselt(double Nu, double k, double L) {
  if (!validation::is_finite_positive(L)) {
    throw std::invalid_argument(
        "htc_from_nusselt: characteristic length must be positive");
  }
  if (!validation::is_finite_positive(k)) {
    throw std::invalid_argument(
        "htc_from_nusselt: thermal conductivity must be positive");
  }
  return Nu * k / L;
}

// -------------------------------------------------------------
// Overall Heat Transfer Coefficient
// -------------------------------------------------------------

double overall_htc(const std::vector<double> &h_values,
                   const std::vector<double> &t_over_k) {
  if (h_values.empty()) {
    throw std::invalid_argument("overall_htc: at least one h value required");
  }

  double R_total = 0.0;

  // Add convective resistances: 1/h
  for (double h : h_values) {
    if (!(h > 0.0)) {
      throw std::invalid_argument("overall_htc: h values must be positive");
    }
    R_total += 1.0 / h;
  }

  // Add conductive resistances: t/k
  for (double tk : t_over_k) {
    if (!validation::is_finite_nonnegative(tk)) {
      throw std::invalid_argument(
          "overall_htc: t/k values must be non-negative");
    }
    R_total += tk;
  }

  return 1.0 / R_total;
}

double overall_htc_wall(double h_inner, double h_outer, double t_wall,
                        double k_wall, double R_fouling) {
  if (!validation::is_finite_positive(k_wall)) {
    throw std::invalid_argument("overall_htc_wall: k_wall must be positive");
  }
  return overall_htc_wall(h_inner, h_outer, std::vector<double>{t_wall / k_wall}, R_fouling);
}

double overall_htc_wall(double h_inner, double h_outer,
                        const std::vector<double> &t_over_k_layers) {
  return overall_htc({h_inner, h_outer}, t_over_k_layers);
}

double overall_htc_wall(double h_inner, double h_outer,
                        const std::vector<double> &t_over_k_layers,
                        double R_fouling) {
  if (!validation::is_finite_nonnegative(R_fouling)) {
    throw std::invalid_argument(
        "overall_htc_wall: R_fouling must be non-negative");
  }
  // Append fouling resistance to layers
  std::vector<double> all_resistances = t_over_k_layers;
  all_resistances.push_back(R_fouling);
  return overall_htc({h_inner, h_outer}, all_resistances);
}

double thermal_resistance(double h, double A) {
  if (!(h > 0.0) || !validation::is_finite_positive(A)) {
    throw std::invalid_argument("thermal_resistance: h and A must be positive");
  }
  return 1.0 / (h * A);
}

double thermal_resistance_wall(double thickness, double k, double A) {
  if (!validation::is_finite_nonnegative(thickness) || !validation::is_finite_positive(k) || !validation::is_finite_positive(A)) {
    throw std::invalid_argument("thermal_resistance_wall: invalid parameters");
  }
  return thickness / (k * A);
}

// -------------------------------------------------------------
// Log Mean Temperature Difference (LMTD)
// -------------------------------------------------------------

double lmtd(double dT1, double dT2) {
  if (!validation::is_finite_positive(dT1) || !validation::is_finite_positive(dT2)) {
    throw std::invalid_argument(
        "lmtd: temperature differences must be positive");
  }

  // If dT1 ~= dT2, use arithmetic mean to avoid 0/0
  double ratio = dT1 / dT2;
  if (std::abs(ratio - 1.0) < 1e-6) {
    return (dT1 + dT2) / 2.0;
  }

  return (dT1 - dT2) / std::log(ratio);
}

double lmtd_counterflow(double T_hot_in, double T_hot_out, double T_cold_in,
                        double T_cold_out) {
  double dT1 = T_hot_in - T_cold_out;
  double dT2 = T_hot_out - T_cold_in;
  return lmtd(dT1, dT2);
}

double lmtd_parallelflow(double T_hot_in, double T_hot_out, double T_cold_in,
                         double T_cold_out) {
  double dT1 = T_hot_in - T_cold_in;
  double dT2 = T_hot_out - T_cold_out;
  return lmtd(dT1, dT2);
}

// -------------------------------------------------------------
// Temperature Profile Through Wall
// -------------------------------------------------------------

std::vector<double>
wall_temperature_profile(double T_hot, double T_cold, double h_hot,
                         double h_cold, const std::vector<double> &t_over_k,
                         double R_fouling, double &q) {

  if (!(h_hot > 0.0) || !(h_cold > 0.0)) {
    throw std::invalid_argument(
        "wall_temperature_profile: HTCs must be positive");
  }

  // Total thermal resistance per unit area [m^2*K/W]
  double R_total = 1.0 / h_hot + 1.0 / h_cold + R_fouling;
  for (double r : t_over_k) {
    if (!validation::is_finite_nonnegative(r)) {
      throw std::invalid_argument(
          "wall_temperature_profile: t/k values must be non-negative");
    }
    R_total += r;
  }

  // Heat flux (constant through all layers)
  q = (T_hot - T_cold) / R_total;

  // Build temperature profile
  // Start from hot side, subtract temperature drops
  std::vector<double> temps;
  temps.reserve(t_over_k.size() + 1);

  // Hot surface temperature (after convective resistance)
  double T = T_hot - q / h_hot;
  temps.push_back(T);

  // Temperature at each layer interface
  for (double r : t_over_k) {
    T -= q * r;
    temps.push_back(T);
  }

  return temps;
}

std::vector<double>
wall_temperature_profile(double T_hot, double T_cold, double h_hot,
                         double h_cold, const std::vector<double> &t_over_k,
                         double R_fouling) {
  double q_unused;
  return wall_temperature_profile(T_hot, T_cold, h_hot, h_cold, t_over_k,
                                  R_fouling, q_unused);
}

// -------------------------------------------------------------
// Heat Flux from Measured Temperature
// -------------------------------------------------------------

double heat_flux_from_T_at_edge(double T_measured, std::size_t edge_idx,
                                double T_hot, double T_cold, double h_hot,
                                double h_cold,
                                const std::vector<double> &t_over_k) {

  std::size_t n_layers = t_over_k.size();
  if (edge_idx > n_layers) {
    throw std::invalid_argument(
        "heat_flux_from_T_at_edge: edge_idx out of range (max = " +
        std::to_string(n_layers) + ")");
  }

  if (!(h_hot > 0.0) || !(h_cold > 0.0)) {
    throw std::invalid_argument(
        "heat_flux_from_T_at_edge: HTCs must be positive");
  }

  // Thermal resistance from hot bulk to edge
  double R_hot_to_edge = 1.0 / h_hot; // convective resistance
  for (std::size_t i = 0; i < edge_idx; ++i) {
    R_hot_to_edge += t_over_k[i];
  }

  // Thermal resistance from edge to cold bulk
  double R_edge_to_cold = 0.0;
  for (std::size_t i = edge_idx; i < n_layers; ++i) {
    R_edge_to_cold += t_over_k[i];
  }
  R_edge_to_cold += 1.0 / h_cold; // convective resistance

  // Heat flux: q = (T_hot - T_measured) / R_hot_to_edge
  //          = (T_measured - T_cold) / R_edge_to_cold
  // Both should give same q in steady state
  // Use the side with larger resistance for better numerical stability
  if (R_hot_to_edge >= R_edge_to_cold) {
    return (T_hot - T_measured) / R_hot_to_edge;
  } else {
    return (T_measured - T_cold) / R_edge_to_cold;
  }
}

double heat_flux_from_T_at_depth(double T_measured, double depth_from_hot,
                                 double T_hot, [[maybe_unused]] double T_cold,
                                 double h_hot, double h_cold,
                                 const std::vector<double> &thicknesses,
                                 const std::vector<double> &conductivities) {

  if (thicknesses.size() != conductivities.size()) {
    throw std::invalid_argument("heat_flux_from_T_at_depth: thicknesses and "
                                "conductivities must have same size");
  }
  if (!validation::is_finite_nonnegative(depth_from_hot)) {
    throw std::invalid_argument(
        "heat_flux_from_T_at_depth: depth must be non-negative");
  }
  if (!(h_hot > 0.0) || !(h_cold > 0.0)) {
    throw std::invalid_argument(
        "heat_flux_from_T_at_depth: HTCs must be positive");
  }

  // Calculate total wall thickness
  double total_thickness = 0.0;
  for (double t : thicknesses) {
    total_thickness += t;
  }

  if (depth_from_hot > total_thickness) {
    throw std::invalid_argument(
        "heat_flux_from_T_at_depth: depth exceeds total wall thickness (" +
        std::to_string(total_thickness) + " m)");
  }

  // Find which layer the depth is in and compute resistance to that point
  double R_hot_to_point = 1.0 / h_hot; // convective resistance
  double cumulative_depth = 0.0;

  for (std::size_t i = 0; i < thicknesses.size(); ++i) {
    double layer_end = cumulative_depth + thicknesses[i];

    if (depth_from_hot <= layer_end) {
      // Point is in this layer
      double depth_in_layer = depth_from_hot - cumulative_depth;
      R_hot_to_point += depth_in_layer / conductivities[i];
      break;
    } else {
      // Point is past this layer
      R_hot_to_point += thicknesses[i] / conductivities[i];
      cumulative_depth = layer_end;
    }
  }

  // q = (T_hot - T_measured) / R_hot_to_point
  return (T_hot - T_measured) / R_hot_to_point;
}

double bulk_T_from_edge_T_and_q(double T_measured, std::size_t edge_idx,
                                double q, double h_hot, double h_cold,
                                const std::vector<double> &t_over_k,
                                const std::string &solve_for) {

  std::size_t n_layers = t_over_k.size();
  if (edge_idx > n_layers) {
    throw std::invalid_argument(
        "bulk_T_from_edge_T_and_q: edge_idx out of range (max = " +
        std::to_string(n_layers) + ")");
  }

  if (!(h_hot > 0.0) || !(h_cold > 0.0)) {
    throw std::invalid_argument(
        "bulk_T_from_edge_T_and_q: HTCs must be positive");
  }

  if (solve_for == "hot") {
    // T_hot = T_measured + q * R_hot_to_edge
    double R_hot_to_edge = 1.0 / h_hot;
    for (std::size_t i = 0; i < edge_idx; ++i) {
      R_hot_to_edge += t_over_k[i];
    }
    return T_measured + q * R_hot_to_edge;
  } else if (solve_for == "cold") {
    // T_cold = T_measured - q * R_edge_to_cold
    double R_edge_to_cold = 0.0;
    for (std::size_t i = edge_idx; i < n_layers; ++i) {
      R_edge_to_cold += t_over_k[i];
    }
    R_edge_to_cold += 1.0 / h_cold;
    return T_measured - q * R_edge_to_cold;
  } else {
    throw std::invalid_argument(
        "bulk_T_from_edge_T_and_q: solve_for must be 'hot' or 'cold'");
  }
}

// -------------------------------------------------------------
// Temperature Sensitivity
// -------------------------------------------------------------

std::pair<double, double> dT_edge_dT_bulk(std::size_t edge_idx, double h_hot,
                                          double h_cold,
                                          const std::vector<double> &t_over_k) {

  std::size_t n_layers = t_over_k.size();
  if (edge_idx > n_layers) {
    throw std::invalid_argument(
        "dT_edge_dT_bulk: edge_idx out of range (max = " +
        std::to_string(n_layers) + ")");
  }

  if (!(h_hot > 0.0) || !(h_cold > 0.0)) {
    throw std::invalid_argument("dT_edge_dT_bulk: HTCs must be positive");
  }

  // Total thermal resistance
  double R_total = 1.0 / h_hot + 1.0 / h_cold;
  for (double r : t_over_k) {
    R_total += r;
  }

  // Resistance from hot bulk to edge
  double R_hot_to_edge = 1.0 / h_hot;
  for (std::size_t i = 0; i < edge_idx; ++i) {
    R_hot_to_edge += t_over_k[i];
  }

  // Resistance from edge to cold bulk
  double R_edge_to_cold = R_total - R_hot_to_edge;

  // Sensitivities
  // T_edge = T_hot * (R_edge_to_cold / R_total) + T_cold * (R_hot_to_edge /
  // R_total)
  double dT_dT_hot = R_edge_to_cold / R_total;
  double dT_dT_cold = R_hot_to_edge / R_total;

  return {dT_dT_hot, dT_dT_cold};
}

double dT_edge_dT_hot(std::size_t edge_idx, double h_hot, double h_cold,
                      const std::vector<double> &t_over_k) {
  return dT_edge_dT_bulk(edge_idx, h_hot, h_cold, t_over_k).first;
}

double dT_edge_dT_cold(std::size_t edge_idx, double h_hot, double h_cold,
                       const std::vector<double> &t_over_k) {
  return dT_edge_dT_bulk(edge_idx, h_hot, h_cold, t_over_k).second;
}

double dT_edge_dq(std::size_t edge_idx, double h_hot,
                  const std::vector<double> &t_over_k) {

  std::size_t n_layers = t_over_k.size();
  if (edge_idx > n_layers) {
    throw std::invalid_argument("dT_edge_dq: edge_idx out of range (max = " +
                                std::to_string(n_layers) + ")");
  }

  if (!(h_hot > 0.0)) {
    throw std::invalid_argument("dT_edge_dq: h_hot must be positive");
  }

  // T_edge = T_hot - q * R_hot_to_edge
  // dT_edge/dq = -R_hot_to_edge
  double R_hot_to_edge = 1.0 / h_hot;
  for (std::size_t i = 0; i < edge_idx; ++i) {
    R_hot_to_edge += t_over_k[i];
  }

  return -R_hot_to_edge;
}

// -------------------------------------------------------------
// NTU-Effectiveness Method
// -------------------------------------------------------------

double effectiveness_counterflow(double NTU, double C_r) {
  if (!(NTU >= 0.0)) {
    throw std::invalid_argument(
        "effectiveness_counterflow: NTU must be non-negative");
  }
  if (!validation::is_finite_nonnegative(C_r) || C_r > 1) {
    throw std::invalid_argument(
        "effectiveness_counterflow: C_r must be in [0, 1]");
  }

  // Special case: C_r = 1 (balanced heat exchanger)
  if (std::abs(C_r - 1.0) < 1e-10) {
    return NTU / (1.0 + NTU);
  }

  // Special case: C_r = 0 (one fluid has infinite capacity, e.g.,
  // condenser/evaporator)
  if (C_r < 1e-10) {
    return 1.0 - std::exp(-NTU);
  }

  // General case
  double exp_term = std::exp(-NTU * (1.0 - C_r));
  return (1.0 - exp_term) / (1.0 - C_r * exp_term);
}

double effectiveness_parallelflow(double NTU, double C_r) {
  if (!(NTU >= 0.0)) {
    throw std::invalid_argument(
        "effectiveness_parallelflow: NTU must be non-negative");
  }
  if (!validation::is_finite_nonnegative(C_r) || C_r > 1) {
    throw std::invalid_argument(
        "effectiveness_parallelflow: C_r must be in [0, 1]");
  }

  // Special case: C_r = 0
  if (C_r < 1e-10) {
    return 1.0 - std::exp(-NTU);
  }

  // General case
  return (1.0 - std::exp(-NTU * (1.0 + C_r))) / (1.0 + C_r);
}

// -------------------------------------------------------------
// State-based convenience functions
// -------------------------------------------------------------

// VDI Heat Atlas transition (vdi_transition, heat_transfer.h) for a
// turbulent-only correlation nu_turb(Re): laminar Nu below 2300, the linear
// interpolation to nu_turb(1e4) up to 1e4, nu_turb above.
template <class Fn>
static double vdi_transition_blend(double Re, bool heating, Fn &&nu_turb) {
  const double Nu_lam = heating ? NU_LAMINAR_CONST_T : NU_LAMINAR_CONST_Q;
  if (Re < vdi_transition::re_laminar) {
    return Nu_lam;
  }
  if (Re >= vdi_transition::re_turbulent) {
    return nu_turb(Re);
  }
  const double g = (Re - vdi_transition::re_laminar) /
                   (vdi_transition::re_turbulent - vdi_transition::re_laminar);
  return (1.0 - g) * Nu_lam + g * nu_turb(vdi_transition::re_turbulent);
}

double nusselt_circular_channel(const State &s, double velocity, double diameter,
                                bool heating, double roughness) {
    if (!validation::is_finite_positive(velocity)) {
        throw std::invalid_argument("nusselt_circular_channel: velocity must be positive");
    }
    if (!validation::is_finite_positive(diameter)) {
        throw std::invalid_argument("nusselt_circular_channel: diameter must be positive");
    }

  // Re = rho * V * D / mu
  double rho = s.rho();
  double mu = s.mu();
  double Re = rho * velocity * diameter / mu;
  double Pr = s.Pr();

  // Laminar below 2300, Gnielinski above 3000, C1 between (#448).
  const double Nu_lam = heating ? NU_LAMINAR_CONST_T : NU_LAMINAR_CONST_Q;
  const auto ft = friction_turbulent_and_derivative(Re, roughness / diameter);
  return nusselt_channel_gnielinski_and_derivative(Re, Pr, ft.f, ft.df_dRe,
                                                    Nu_lam)
      .Nu;
}

double htc_circular_channel(const State &s, double velocity, double diameter,
                            bool heating, double roughness) {
    double Nu = nusselt_circular_channel(s, velocity, diameter, heating, roughness);
    double k = s.k();
    return htc_from_nusselt(Nu, k, diameter);
}

std::tuple<double, double, double>
htc_circular_channel(double T, double P, const std::vector<double> &X,
                     double velocity, double diameter,
                     const std::string &correlation, bool heating,
                     double mu_ratio, double roughness) {
    if (!validation::is_finite_positive(velocity)) {
        throw std::invalid_argument("htc_circular_channel: velocity must be positive");
    }
    if (!validation::is_finite_positive(diameter)) {
        throw std::invalid_argument("htc_circular_channel: diameter must be positive");
    }
    if (!validation::is_finite_positive(T)) {
        throw std::invalid_argument("htc_circular_channel: temperature must be positive");
    }
    if (!validation::is_finite_positive(P)) {
        throw std::invalid_argument("htc_circular_channel: pressure must be positive");
    }

  // Compute thermodynamic and transport properties
  double rho = density(T, P, X);
  double mu = viscosity(T, P, X);
  double k = thermal_conductivity(T, P, X);
  [[maybe_unused]] double cp = cp_mass(T, X);
  double Pr = prandtl(T, P, X);

  // Reynolds number
  double Re = rho * velocity * diameter / mu;

  // Select correlation and compute Nusselt number
  double Nu;

  if (correlation == "gnielinski") {
    // Laminar below 2300, Gnielinski above 3000, C1 between (#448).
    const double Nu_lam = heating ? NU_LAMINAR_CONST_T : NU_LAMINAR_CONST_Q;
    const auto ft = friction_turbulent_and_derivative(Re, roughness / diameter);
    Nu = nusselt_channel_gnielinski_and_derivative(Re, Pr, ft.f, ft.df_dRe,
                                                    Nu_lam)
             .Nu;
  } else if (correlation == "dittus_boelter") {
    // Dittus-Boelter: Re > 10000, 0.6 < Pr < 160
    Nu = vdi_transition_blend(Re, heating, [&](double Re_t) {
      return nusselt_dittus_boelter(Re_t, Pr, heating);
    });
  } else if (correlation == "sieder_tate") {
    // Sieder-Tate: Re > 10000, 0.7 < Pr < 16700
    Nu = vdi_transition_blend(Re, heating, [&](double Re_t) {
      return nusselt_sieder_tate(Re_t, Pr, mu_ratio);
    });
  } else if (correlation == "petukhov") {
    // Petukhov: 1e4 < Re < 5e6, 0.5 < Pr < 2000
    Nu = vdi_transition_blend(Re, heating, [&](double Re_t) {
      return nusselt_petukhov(Re_t, Pr, friction_petukhov(Re_t));
    });
  } else {
    throw std::invalid_argument("htc_circular_channel: unknown correlation '" +
                                correlation + "'. " +
                                "Valid options: 'gnielinski', "
                                "'dittus_boelter', 'sieder_tate', 'petukhov'");
  }

  // Heat transfer coefficient
  double h = htc_from_nusselt(Nu, k, diameter);

  return std::make_tuple(h, Nu, Re);
}

// -------------------------------------------------------------
// Internal helper: build ChannelResult from common fields
// -------------------------------------------------------------

static ChannelResult make_channel_result(double h, double Nu, double Re,
                                         double Pr, double f, double dP,
                                         double M, double T_aw, double T_hot) {
  double q = std::numeric_limits<double>::quiet_NaN();
  if (std::isfinite(T_hot)) {
    q = h * (T_aw - T_hot);
  }
  return ChannelResult{h, Nu, Re, Pr, f, dP, M, T_aw, q};
}


// -------------------------------------------------------------
// channel_smooth
// -------------------------------------------------------------

ChannelResult
channel_smooth(double T, double P, const std::vector<double> &X,
               double velocity, double diameter, double length,
               double T_hot, const std::string &correlation,
                             bool heating, double mu_ratio, double roughness,
                             double Nu_multiplier, double f_multiplier,
                             double flow_area) {
  // Use absolute velocity: Re, Nu, f, and dP are sign-independent.
  // Negative velocity can arise from reversed flow during Newton iteration
  // or from negative mass-flow boundary conditions.
  velocity = std::abs(velocity);
  if (!validation::is_finite_positive(diameter)) {
    throw std::invalid_argument("channel_smooth: diameter must be positive");
  }
  if (!validation::is_finite_nonnegative(length)) {
    throw std::invalid_argument("channel_smooth: length must be non-negative");
  }

  // Fluid properties and their derivatives
  auto [rho, drho_dT, drho_dP] = solver::density_and_jacobians(T, P, X);
  auto [mu, dmu_dT, dmu_dP] = solver::viscosity_and_jacobians(T, P, X);
  double k = thermal_conductivity(T, P, X);
  double Pr = prandtl(T, P, X);

  // Reynolds number
  double Re = (velocity > 0.0) ? rho * velocity * diameter / mu : 0.0;

  // Friction factor: 64/Re below 2300, turbulent above 3000, C1 between;
  // rough walls blend Petukhov into Colebrook over 3000-4000 (#448).
  const double e_D = roughness / diameter;
  const FrictionAndDerivative f_chan = friction_channel_and_derivative(Re, e_D);
  // f_multiplier applies to every friction the correlations see, so the
  // Gnielinski/Petukhov Nu and its derivative use the multiplied value.
  double f = f_chan.f * f_multiplier;

  // Nusselt number
  const double Nu_lam = heating ? NU_LAMINAR_CONST_T : NU_LAMINAR_CONST_Q;
  FrictionAndDerivative f_turb{};
  NuAndDerivative nu_gn{};
  double Nu;
  // Turbulent-only correlations: d(Nu)/dRe and d(Nu)/dPr before Nu_multiplier.
  double dNu_dRe_raw = 0.0;
  double dNu_dPr_raw = 0.0;
  if (correlation == "gnielinski") {
    // Gnielinski sees the TURBULENT friction, never the laminar-blended one.
    f_turb = friction_turbulent_and_derivative(Re, e_D);
    f_turb.f *= f_multiplier;
    f_turb.df_dRe *= f_multiplier;
    nu_gn = nusselt_channel_gnielinski_and_derivative(Re, Pr, f_turb.f,
                                                      f_turb.df_dRe, Nu_lam);
    Nu = nu_gn.Nu;
  } else if (correlation == "dittus_boelter" || correlation == "sieder_tate" ||
             correlation == "petukhov") {
    // Turbulent-only (Re >= 1e4): exact there; laminar below 2300; the VDI
    // Heat Atlas interpolation between (vdi_transition, heat_transfer.h).
    // Each partial is taken at Re_t = max(Re, 1e4), where the correlation is
    // evaluated.
    struct Turb {
      double Nu, dRe, dPr;
    };
    auto turb = [&](double Re_t) -> Turb {
      if (correlation == "dittus_boelter") {
        const double n = heating ? 0.4 : 0.3;
        const double v = nusselt_dittus_boelter(Re_t, Pr, heating);
        return {v, 0.8 * v / Re_t, n * v / Pr};
      }
      if (correlation == "sieder_tate") {
        const double v = nusselt_sieder_tate(Re_t, Pr, mu_ratio);
        return {v, 0.8 * v / Re_t, v / (3.0 * Pr)};
      }
      // Petukhov with the channel friction it sees (f_multiplier included).
      auto f_at = [&](double R) { return friction_channel_and_derivative(R, e_D).f * f_multiplier; };
      const double v = nusselt_petukhov(Re_t, Pr, f_at(Re_t));
      // The stencil straddles 1e4 at the band edge: a status pointer keeps
      // its sub-range point from warning (the value above still does).
      CorrelationStatus quiet{};
      const double eps = std::max(1.0, Re_t * 1e-6);
      const double dRe = (nusselt_petukhov(Re_t + eps, Pr, f_at(Re_t + eps), &quiet) -
                          nusselt_petukhov(Re_t - eps, Pr, f_at(Re_t - eps), &quiet)) /
                         (2.0 * eps);
      const double eps_Pr = std::max(1e-6, Pr * 1e-6);
      const double dPr = (nusselt_petukhov(Re_t, Pr + eps_Pr, f_at(Re_t), &quiet) -
                          nusselt_petukhov(Re_t, Pr - eps_Pr, f_at(Re_t), &quiet)) /
                         (2.0 * eps_Pr);
      return {v, dRe, dPr};
    };
    if (Re < vdi_transition::re_laminar) {
      Nu = Nu_lam;
    } else if (Re < vdi_transition::re_turbulent) {
      const Turb t = turb(vdi_transition::re_turbulent);
      const double span = vdi_transition::re_turbulent - vdi_transition::re_laminar;
      const double g = (Re - vdi_transition::re_laminar) / span;
      Nu = (1.0 - g) * Nu_lam + g * t.Nu;
      dNu_dRe_raw = (t.Nu - Nu_lam) / span;
      dNu_dPr_raw = g * t.dPr;
    } else {
      const Turb t = turb(Re);
      Nu = t.Nu;
      dNu_dRe_raw = t.dRe;
      dNu_dPr_raw = t.dPr;
    }
  } else {
    throw std::invalid_argument("channel_smooth: unknown correlation '" +
                                correlation + "'");
  }

  // Apply Nu multiplier
  Nu *= Nu_multiplier;

  double h = htc_from_nusselt(Nu, k, diameter);
  double dP = (velocity > 0.0)
                  ? f * (length / diameter) * (rho * velocity * velocity / 2.0)
                  : 0.0;

  // Mach number and T_aw -- always continuous
  double a = speed_of_sound(T, X);
  double M = (a > 0.0) ? velocity / a : 0.0;
  const bool turbulent_flow = true;
  double T_aw = T_adiabatic_wall(T, velocity, T, P, X, turbulent_flow);

  // Compute Jacobians
  // Note: For network solver, mdot is the independent variable, not velocity
  // velocity = mdot / (rho * A_cross), where A_cross = pi/4 * D^2
  // We compute derivatives w.r.t. mdot and T

  ChannelResult result = make_channel_result(h, Nu, Re, Pr, f, dP, M, T_aw, T_hot);

  // Only compute Jacobians if flow exists
  if (velocity > 0.0 && Re > 0.0) {
    // The flow area velocity was taken on: the caller's when given (a
    // non-circular channel), else circular (#463).
    double A_cross = (std::isfinite(flow_area) && flow_area > 0.0)
                         ? flow_area
                         : M_PI / 4.0 * diameter * diameter;
    double mdot = rho * velocity * A_cross;

    // Thermal conductivity and Prandtl derivatives via central FD
    const double eps_T_props = 1e-3;
    double k_plus = thermal_conductivity(T + eps_T_props, P, X);
    double k_minus = thermal_conductivity(T - eps_T_props, P, X);
    double dk_dT = (k_plus - k_minus) / (2.0 * eps_T_props);

    double Pr_plus = prandtl(T + eps_T_props, P, X);
    double Pr_minus = prandtl(T - eps_T_props, P, X);
    double dPr_dT = (Pr_plus - Pr_minus) / (2.0 * eps_T_props);

    // Chain rule: dRe/dmdot and dRe/dT
    // Re = rho * v * D / mu = mdot * D / (A_cross * mu)
    // At constant mdot: dRe/dT = -Re * dmu/dT / mu (no rho dependence)
    double dRe_dmdot = diameter / (A_cross * mu);
    double dRe_dT = -Re * dmu_dT / mu;

    // Get Nu and f derivatives w.r.t. Re and Pr (these depend on correlation)
    double dNu_dRe = 0.0;
    double dNu_dPr = 0.0;
    double df_dRe = 0.0;

    // Friction derivative: exact in every regime and band (#448) -- the
    // previous Petukhov formula ignored the Re 3000 clamp and Colebrook used
    // a central difference.
    df_dRe = f_chan.df_dRe * f_multiplier;

    if (correlation == "gnielinski") {
      // Exact total derivative through the laminar/turbulent blend and the
      // friction chain (previously a central difference whose stencil used
      // the unclamped Petukhov friction).
      dNu_dRe = nu_gn.dNu_dRe * Nu_multiplier;

      // dNu/dPr via central FD at fixed Re and turbulent friction
      double eps_Pr = std::max(1e-6, Pr * 1e-6);
      double Nu_Pr_plus = nusselt_channel_gnielinski_and_derivative(
                              Re, Pr + eps_Pr, f_turb.f, f_turb.df_dRe, Nu_lam)
                              .Nu;
      double Nu_Pr_minus = nusselt_channel_gnielinski_and_derivative(
                               Re, Pr - eps_Pr, f_turb.f, f_turb.df_dRe, Nu_lam)
                               .Nu;
      dNu_dPr = (Nu_Pr_plus - Nu_Pr_minus) / (2.0 * eps_Pr) * Nu_multiplier;
    } else {
      // Turbulent-only correlations: partials from the value section (the
      // VDI transition's are exact; laminar's are zero).
      dNu_dRe = dNu_dRe_raw * Nu_multiplier;
      dNu_dPr = dNu_dPr_raw * Nu_multiplier;
    }

    // dh/dmdot = (dh/dNu) * (dNu/dRe) * (dRe/dmdot) = (k/D) * dNu/dRe * dRe/dmdot
    result.dh_dmdot = (k / diameter) * dNu_dRe * dRe_dmdot;

    // dh/dT: Full chain rule h = Nu(Re, Pr) * k(T) / D
    // dh/dT = (1/D) * [k * (dNu/dRe * dRe/dT + dNu/dPr * dPr/dT) + Nu * dk/dT]
    result.dh_dT = (1.0 / diameter) * (k * (dNu_dRe * dRe_dT + dNu_dPr * dPr_dT) + Nu * dk_dT);

    // ddP/dmdot: dP = f * (L/D) * (rho * v^2 / 2)
    // v = mdot / (rho * A), so v^2 = mdot^2 / (rho^2 * A^2)
    // dP = f * (L/D) * mdot^2 / (2 * rho * A^2)
    // d(dP)/dmdot = (df/dRe * dRe/dmdot) * (L/D) * mdot^2/(2*rho*A^2) + f*(L/D)*mdot/(rho*A^2)
    double dP_factor = (length / diameter) / (2.0 * rho * A_cross * A_cross);
    result.ddP_dmdot = df_dRe * dRe_dmdot * dP_factor * mdot * mdot
                     + f * (length / diameter) * mdot / (rho * A_cross * A_cross);
    // mdot = rho * velocity * A_cross at fixed rho, so d(mdot)/d(velocity) = rho * A_cross
    result.ddP_dvelocity = result.ddP_dmdot * rho * A_cross;

    // d(dP)/dT: includes df/dT and drho/dT terms
    // dP = f * (L/D) * mdot^2 / (2*rho*A^2), so d(dP)/dT has df/dT and d(1/rho)/dT = -drho/dT/rho^2
    result.ddP_dT = df_dRe * dRe_dT * dP_factor * mdot * mdot
                  - f * (length / diameter) * velocity * velocity / 2.0 * drho_dT;

    // dq/dmdot and dq/dT: q = h * (T_aw - T_hot)
    // Full product rule: dq/dx = dh/dx * (T_aw - T_hot) + h * dT_aw/dx
    // T_aw derivatives -- always computed so Python relay can use them
    // T_aw = T + r * v^2 / (2*cp), where v = mdot/(rho*A), r = Pr^(1/3)
    double cp_mass_val = cp_mass(T, X);
    double r = std::cbrt(Pr);

    // dT_aw/dmdot = r * v / (cp * rho * A)
    result.dT_aw_dmdot = r * velocity / (cp_mass_val * rho * A_cross);

    // dT_aw/dT at constant mdot via central FD
    const double eps_T_aw = 0.5;
    double rho_plus = density(T + eps_T_aw, P, X);
    double rho_minus = density(T - eps_T_aw, P, X);
    double v_plus = mdot / (rho_plus * A_cross);
    double v_minus = mdot / (rho_minus * A_cross);
    double T_aw_plus = T_adiabatic_wall(T + eps_T_aw, v_plus, T + eps_T_aw, P, X, turbulent_flow);
    double T_aw_minus = T_adiabatic_wall(T - eps_T_aw, v_minus, T - eps_T_aw, P, X, turbulent_flow);
    result.dT_aw_dT = (T_aw_plus - T_aw_minus) / (2.0 * eps_T_aw);
    // At fixed mdot, v ~ 1/rho ~ 1/P (ideal gas): T_aw - T ~ 1/P^2 (#496).
    result.dT_aw_dP = (P > 0.0) ? -2.0 * (T_aw - T) / P : 0.0;

    // Composition (#717), at fixed mdot, T and P. Re = mdot D/(A mu) moves
    // with mu alone; Pr = cp mu/k; h = Nu k/D; T_aw - T = r v^2/(2 cp) with
    // r = Pr^(1/3) and v = mdot/(rho A), rho ~ MW. Mole-fraction partials
    // (cp and MW taken linear in X), then mapped to mass fractions.
    {
      const auto tx = transport_and_dX(T, P, X);
      const std::size_t n = X.size();
      const double MW = mwmix(X) / 1000.0;  // kg/mol
      const double rise = T_aw - T;
      std::vector<double> dh_dX(n, 0.0);
      std::vector<double> dTaw_dX(n, 0.0);
      for (std::size_t j = 0; j < n; ++j) {
        const double W = species_molar_mass(j) / 1000.0;
        const double dcp = (cp_species(j, T) - cp_mass_val * W) / MW;
        const double dRe = -Re * tx.dmu_dX[j] / mu;
        const double dPr = Pr * (dcp / cp_mass_val + tx.dmu_dX[j] / mu - tx.dk_dX[j] / k);
        dh_dX[j] = (k * (dNu_dRe * dRe + dNu_dPr * dPr) + Nu * tx.dk_dX[j]) / diameter;
        dTaw_dX[j] = rise * (dPr / (3.0 * Pr) - 2.0 * W / MW - dcp / cp_mass_val);
      }
      result.dh_dY = mole_gradient_to_mass(X, dh_dX);
      result.dT_aw_dY = mole_gradient_to_mass(X, dTaw_dX);
    }

    if (std::isfinite(T_hot)) {
      double dT_diff = result.T_aw - T_hot;
      result.dq_dmdot = result.dh_dmdot * dT_diff + h * result.dT_aw_dmdot;
      result.dq_dT = result.dh_dT * dT_diff + h * result.dT_aw_dT;
      result.dq_dT_hot = -h;
    }
  }

  return result;
}

// -------------------------------------------------------------
// wall_coupling_and_jacobian
// -------------------------------------------------------------

WallCouplingResult wall_coupling_and_jacobian(
    double h_a, double T_aw_a,
    double h_b, double T_aw_b,
    double t_over_k,
    double A) {
  return wall_coupling_and_jacobian(h_a, T_aw_a, h_b, T_aw_b,
                                    std::vector<double>{t_over_k}, A, 0.0);
}

namespace {
// h above WALL_HTC_KNEE as is; below it knee^2 / (2 knee - h): C1 at the
// knee, monotone, strictly positive and decaying only like 1/|h| (an
// exponential tail underflows to 0 within a few knees and makes R infinite),
// so an iterate whose correlation went to h <= 0 still sees a finite wall
// and a derivative of the h actually used (#481).
std::pair<double, double> floored_htc(double h) {
  if (h >= WALL_HTC_KNEE) {
    return {h, 1.0};
  }
  const double k = WALL_HTC_KNEE;
  const double d = 2.0 * k - h;
  return {k * k / d, k * k / (d * d)};
}
}  // namespace

WallCouplingResult wall_coupling_and_jacobian(
    double h_a, double T_aw_a,
    double h_b, double T_aw_b,
    const std::vector<double> &t_over_k_layers,
    double A,
    double R_fouling) {
  WallCouplingResult result;

  // Total conductive resistance
  double R_wall = R_fouling;
  for (double tk : t_over_k_layers) {
    R_wall += tk;
  }

  const auto [ha, dha] = floored_htc(h_a);
  const auto [hb, dhb] = floored_htc(h_b);
  const double R_a = 1.0 / ha;
  const double R_b = 1.0 / hb;
  double R_total = R_a + R_wall + R_b;
  double U = 1.0 / R_total;

  // Heat transfer rate: Q = U * A * (T_aw_a - T_aw_b)
  double dT = T_aw_a - T_aw_b;
  result.Q = U * A * dT;

  // T_hot on side A (surface temperature)
  result.T_hot = T_aw_a - dT * (R_a / R_total);

  // dQ/dh = A dT dU/dh, dU/dh_a = (R_a / R_total)^2 at the floored h, times
  // the floor's own slope.
  const double factor_a = R_a / R_total;
  const double factor_b = R_b / R_total;
  result.dQ_dT_aw_a = U * A;

  result.dQ_dh_a = (factor_a * factor_a) * A * dT * dha;
  result.dQ_dh_b = (factor_b * factor_b) * A * dT * dhb;

  result.dQ_dT_aw_b = -U * A;
  // dQ/dR_wall = -A dT / R_total^2.
  result.dQ_dR_wall = -U * U * A * dT;

  return result;
}

namespace {
// g(x) = (1 - e^-x) / x and its slope as g', x g' and x^2 g' (x >= 0): the
// last two stay finite as x -> inf (0 and -1) where g' alone underflows.
struct ExchangeG {
  double g;
  double gp;
  double xgp;
  double x2gp;
};

ExchangeG exchange_g(double x) {
  if (!std::isfinite(x)) {
    return {0.0, 0.0, 0.0, -1.0};
  }
  if (x < 1e-4) {
    const double g = 1.0 - x / 2.0 + x * x / 6.0;
    const double gp = -0.5 + x / 3.0 - x * x / 8.0;
    return {g, gp, x * gp, x * x * gp};
  }
  const double em = std::exp(-x);
  const double g = -std::expm1(-x) / x;
  return {g, (em - g) / x, em - g, x * (em - g)};
}

// (UA w)^2 g'(x) with x = UA D, w <= ... any inverse capacity: directly where
// x is small (D may vanish), as (w / D)^2 x^2 g'(x) where it is not (g' may
// underflow).
double scaled_gp(double UA, double w, double D, const ExchangeG& gv, double x) {
  if (x < 1.0) {
    const double t = UA * w;
    return t * t * gv.gp;
  }
  const double r = w / D;
  return r * r * gv.x2gp;
}

// K = Q / dT [W/K] and its partials in UA, C_a, C_b.
struct ExchangeK {
  double K;
  double dK_dUA;
  double dK_dCa;
  double dK_dCb;
};

ExchangeK exchange_conductance(double UA, double C_a, double C_b, bool counterflow) {
  if (!(UA > 0.0)) {
    return {0.0, 1.0, 0.0, 0.0};
  }
  // A stream that does not flow takes no heat: Q -> C dT as C -> 0. With
  // both stopped K(C_a, 0) = K(0, C_b) = 0 for any C, so both partials are
  // 0 (#1245); 0.5 each was the slope along C_a = C_b, not a partial.
  if (C_a <= 0.0 && C_b <= 0.0) {
    return {0.0, 0.0, 0.0, 0.0};
  }
  if (C_a <= 0.0) {
    return {0.0, 0.0, 1.0, 0.0};
  }
  if (C_b <= 0.0) {
    return {0.0, 0.0, 0.0, 1.0};
  }
  const double w_a = std::isinf(C_a) ? 0.0 : 1.0 / C_a;
  const double w_b = std::isinf(C_b) ? 0.0 : 1.0 / C_b;
  if (w_a == 0.0 && w_b == 0.0) {
    return {UA, 1.0, 0.0, 0.0};  // two reservoirs: NTU = 0
  }
  if (!counterflow) {
    const double S = w_a + w_b;
    const double x = UA * S;
    const ExchangeG gv = exchange_g(x);
    const double fa = w_a / S;
    const double fb = w_b / S;
    return {UA * gv.g, gv.g + gv.xgp, -fa * fa * gv.x2gp, -fb * fb * gv.x2gp};
  }
  // Counterflow, symmetric in A and B: evaluated with w_hi >= w_lo (the
  // smaller capacity rate on 'hi') so that x >= 0.
  const bool a_hi = w_a >= w_b;
  const double w_hi = a_hi ? w_a : w_b;
  const double w_lo = a_hi ? w_b : w_a;
  const double D = w_hi - w_lo;
  const double x = UA * D;
  const ExchangeG gv = exchange_g(x);
  const double G = UA * gv.g;
  const double den = 1.0 + G * w_lo;
  const double K = G / den;
  const double dK_dG = 1.0 / (den * den);
  const double dK_dUA = dK_dG * (gv.g + gv.xgp);
  const double dK_dC_hi = -dK_dG * scaled_gp(UA, w_hi, D, gv, x);
  const double Gw = G * w_lo;
  const double dK_dC_lo = dK_dG * scaled_gp(UA, w_lo, D, gv, x) + Gw * Gw * dK_dG;
  return {K, dK_dUA, a_hi ? dK_dC_hi : dK_dC_lo, a_hi ? dK_dC_lo : dK_dC_hi};
}
}  // namespace

WallCouplingResult wall_exchange_and_jacobian(
    double h_a, double T_aw_a, double C_a,
    double h_b, double T_aw_b, double C_b,
    const std::vector<double> &t_over_k_layers,
    double A,
    double R_fouling,
    bool counterflow) {
  // UA and its h derivatives: the linear coupling at a unit difference.
  const WallCouplingResult unit =
      wall_coupling_and_jacobian(h_a, 1.0, h_b, 0.0, t_over_k_layers, A, R_fouling);
  const double UA = unit.Q;
  // The surface temperature as the linear coupling places it at the inlets.
  WallCouplingResult result =
      wall_coupling_and_jacobian(h_a, T_aw_a, h_b, T_aw_b, t_over_k_layers, A, R_fouling);

  const ExchangeK k = exchange_conductance(UA, C_a, C_b, counterflow);
  const double dT = T_aw_a - T_aw_b;
  result.Q = k.K * dT;
  result.dQ_dT_aw_a = k.K;
  result.dQ_dT_aw_b = -k.K;
  result.dQ_dh_a = dT * k.dK_dUA * unit.dQ_dh_a;
  result.dQ_dh_b = dT * k.dK_dUA * unit.dQ_dh_b;
  result.dQ_dC_a = dT * k.dK_dCa;
  result.dQ_dC_b = dT * k.dK_dCb;
  result.dQ_dR_wall = dT * k.dK_dUA * unit.dQ_dR_wall;
  result.UA = UA;
  const double C_min = std::min(C_a, C_b);
  result.effectiveness = (C_min > 0.0 && std::isfinite(C_min)) ? k.K / C_min : 0.0;
  return result;
}

std::array<double, 4> wall_resistance_fixed_point_derivatives(
    double T_hot, double T_cold, double h_hot, double h_cold,
    const std::vector<double> &t_over_k, double R_fouling,
    const std::vector<double> &dr_dTm) {
  const std::size_t n = t_over_k.size();
  if (dr_dTm.size() != n) {
    throw std::invalid_argument(
        "wall_resistance_fixed_point_derivatives: dr_dTm needs one entry per layer");
  }
  std::array<double, 4> out{0.0, 0.0, 0.0, 0.0};
  if (std::all_of(dr_dTm.begin(), dr_dTm.end(), [](double g) { return g == 0.0; })) {
    return out;
  }
  double q = 0.0;
  // Validates the inputs and gives q exactly as the caller's profile.
  (void)wall_temperature_profile(T_hot, T_cold, h_hot, h_cold, t_over_k, R_fouling, q);
  double R_total = 1.0 / h_hot + 1.0 / h_cold + R_fouling;
  for (double r : t_over_k) {
    R_total += r;
  }
  // c_i: the resistance from the hot fluid to layer i's middle.
  std::vector<double> c(n);
  double acc = 1.0 / h_hot;
  for (std::size_t i = 0; i < n; ++i) {
    c[i] = acc + 0.5 * t_over_k[i];
    acc += t_over_k[i];
  }
  // Augmented system [I - G M | G P], M = dTm/dr (n x n), P = dTm/dp (n x 4).
  const std::size_t cols = n + 4;
  std::vector<double> a(n * cols, 0.0);
  const double hh2 = h_hot * h_hot;
  const double hc2 = h_cold * h_cold;
  for (std::size_t i = 0; i < n; ++i) {
    const double g = dr_dTm[i];
    double *row = &a[i * cols];
    for (std::size_t j = 0; j < n; ++j) {
      const double dc = j < i ? 1.0 : (j == i ? 0.5 : 0.0);
      const double M_ij = -q * dc + q * c[i] / R_total;  // dq/dr_j = -q / R_total
      row[j] = (i == j ? 1.0 : 0.0) - g * M_ij;
    }
    const double f = c[i] / R_total;
    row[n + 0] = g * (1.0 - f);                    // dTm_i/dT_hot
    row[n + 1] = g * f;                            // dTm_i/dT_cold
    row[n + 2] = g * (q / hh2) * (1.0 - f);        // dTm_i/dh_hot
    row[n + 3] = g * (-q * f / hc2);               // dTm_i/dh_cold
  }
  // Gaussian elimination with partial pivoting (n is the layer count).
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t piv = k;
    for (std::size_t i = k + 1; i < n; ++i) {
      if (std::abs(a[i * cols + k]) > std::abs(a[piv * cols + k])) {
        piv = i;
      }
    }
    if (!(std::abs(a[piv * cols + k]) > 0.0) || !std::isfinite(a[piv * cols + k])) {
      throw std::runtime_error(
          "wall_resistance_fixed_point_derivatives: the k(T) fixed point is not isolated");
    }
    if (piv != k) {
      for (std::size_t j = 0; j < cols; ++j) {
        std::swap(a[k * cols + j], a[piv * cols + j]);
      }
    }
    for (std::size_t i = k + 1; i < n; ++i) {
      const double m = a[i * cols + k] / a[k * cols + k];
      for (std::size_t j = k; j < cols; ++j) {
        a[i * cols + j] -= m * a[k * cols + j];
      }
    }
  }
  for (std::size_t col = 0; col < 4; ++col) {
    std::vector<double> x(n);
    for (std::size_t ii = n; ii-- > 0;) {
      double v = a[ii * cols + n + col];
      for (std::size_t j = ii + 1; j < n; ++j) {
        v -= a[ii * cols + j] * x[j];
      }
      x[ii] = v / a[ii * cols + ii];
    }
    double sum = 0.0;
    for (double v : x) {
      sum += v;
    }
    out[col] = sum;
  }
  return out;
}

} // namespace combaero
