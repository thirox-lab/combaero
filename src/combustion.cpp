#include "../include/combustion.h"
#include "../include/equilibrium.h"
#include "../include/humidair.h"
#include "../include/thermo.h"
#include "../include/composition.h"
#include "../include/thermo_transport_data.h"
#include "../include/input_validation.h"
#include "../include/math_constants.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <stdexcept>

namespace combaero {

namespace {

double pure_species_enthalpy(std::size_t idx,
                             const double reference_temperature) {
  return h_species(idx, reference_temperature);
}

void validate_fuel_mole_fractions(const std::vector<double> &X_fuel) {
  if (X_fuel.size() != species_names.size()) {
    throw std::invalid_argument(
        "fuel_lhv_*: mixture size does not match number of species");
  }
  for (double x : X_fuel) {
    if (!std::isfinite(x) || x < 0.0) {
      throw std::invalid_argument(
          "fuel_lhv_*: mole fractions must be finite and non-negative");
    }
  }
  const double sum = std::accumulate(X_fuel.begin(), X_fuel.end(), 0.0);
  if (std::abs(sum - 1.0) > 1.0e-5) {
    throw std::invalid_argument("fuel_lhv_*: mole fractions must sum to 1.0");
  }
}

// A NaN passes every `x <= bound` guard (it compares false), so the
// helpers below test finiteness first: a NaN phi or target is refused
// instead of returning NaN, or a valid-looking bracket end.
void require_finite(double value, const char *what) {
  if (!std::isfinite(value)) {
    throw std::invalid_argument(std::string(what) + " must be finite");
  }
}

void require_finite_fractions(const std::vector<double> &X, const char *what) {
  for (double x : X) {
    if (!std::isfinite(x)) {
      throw std::invalid_argument(std::string(what) +
                                  ": mole fractions must be finite");
    }
  }
}

// complete_combustion*: a physical reactant state (T > 0, P > 0, finite).
void require_physical_state(const State &in, const char *what) {
  if (!std::isfinite(in.T) || in.T <= 0.0) {
    throw std::invalid_argument(std::string(what) +
                                ": temperature must be finite and positive");
  }
  if (!std::isfinite(in.P) || in.P <= 0.0) {
    throw std::invalid_argument(std::string(what) +
                                ": pressure must be finite and positive");
  }
}

} // namespace

// Calculate oxygen required for complete combustion [mol O2/mol fuel]
double oxygen_required_per_mol_fuel(std::size_t fuel_index) {
  // Validate fuel index
  if (fuel_index >= molecular_structures.size()) {
    throw std::invalid_argument("Invalid fuel index");
  }

  // Get molecular structure
  const Molecular_Structure &fuel = molecular_structures[fuel_index];

  // For complete combustion:
  // C_x H_y O_z N_w + (x + y/4 - z/2) O2 -> x CO2 + (y/2) H2O + (w/2) N2
  // Oxygen required = x + y/4 - z/2 [mol O2/mol fuel]

  double oxygen_required = static_cast<double>(fuel.C) +
                           static_cast<double>(fuel.H) / 4.0 -
                           static_cast<double>(fuel.O) / 2.0;

  // Ensure non-negative value (for cases where fuel already has excess oxygen)
  return std::max(0.0, oxygen_required);
}

// Calculate oxygen required for complete combustion [kg O2/kg fuel]
double oxygen_required_per_kg_fuel(std::size_t fuel_index) {
  // Validate fuel index
  if (fuel_index >= molar_masses.size()) {
    throw std::invalid_argument("Invalid fuel index");
  }

  // Get molar masses
  double fuel_molar_mass = molar_masses[fuel_index];               // g/mol
  double oxygen_molar_mass = molar_masses[species_index.at("O2")]; // g/mol

  // Calculate molar oxygen requirement
  double molar_oxygen_required =
      oxygen_required_per_mol_fuel(fuel_index); // mol O2/mol fuel

  // Convert to mass basis: (mol O2/mol fuel) * (g O2/mol O2) / (g fuel/mol
  // fuel)
  return molar_oxygen_required * oxygen_molar_mass /
         fuel_molar_mass; // kg O2/kg fuel
}

// Calculate oxygen required for complete combustion of a mixture [mol O2/mol
// mixture]
double oxygen_required_per_mol_mixture(const std::vector<double> &X) {
  require_finite_fractions(X, "oxygen_required_per_mol_mixture");
  // Check if mole fractions sum to approximately 1.0
  double sum = std::accumulate(X.begin(), X.end(), 0.0);
  if (std::abs(sum - 1.0) > 1.0e-5) {
    throw std::invalid_argument("Mole fractions must sum to 1.0");
  }

  double total_oxygen_required = 0.0;

  // Calculate weighted sum of oxygen requirements for each fuel component
  for (size_t i = 0; i < X.size(); ++i) {
    // Skip if mole fraction is zero
    if (X[i] <= 0.0)
      continue;

    // Check if this is a fuel (contains carbon or hydrogen)
    if (molecular_structures[i].C > 0 || molecular_structures[i].H > 0) {
      total_oxygen_required += X[i] * oxygen_required_per_mol_fuel(i);
    }
  }

  return total_oxygen_required; // mol O2/mol mixture
}

// Calculate oxygen required for complete combustion of a mixture [kg O2/kg
// mixture]
double oxygen_required_per_kg_mixture(const std::vector<double> &X) {
  require_finite_fractions(X, "oxygen_required_per_kg_mixture");
  // Check if mole fractions sum to approximately 1.0
  double sum = std::accumulate(X.begin(), X.end(), 0.0);
  if (std::abs(sum - 1.0) > 1.0e-5) {
    throw std::invalid_argument("Mole fractions must sum to 1.0");
  }

  // Calculate mixture molecular weight
  double mixture_mw = mwmix(X); // g/mol

  // Get oxygen molecular weight
  double oxygen_mw = molar_masses[species_index.at("O2")]; // g/mol

  // Calculate molar oxygen requirement for the mixture
  double molar_oxygen_required =
      oxygen_required_per_mol_mixture(X); // mol O2/mol mixture

  // If no oxygen required (no fuel in mixture), return 0
  if (molar_oxygen_required <= 0.0) {
    return 0.0;
  }

  // Convert to mass basis: (mol O2/mol mixture) * (g O2/mol O2) / (g
  // mixture/mol mixture)
  return molar_oxygen_required * oxygen_mw / mixture_mw; // kg O2/kg mixture
}

double equivalence_ratio(const std::vector<double> &X) {
  // Guard for malformed input vector size
  if (X.size() != molecular_structures.size()) {
    throw std::invalid_argument("equivalence_ratio: input vector size (" +
                             std::to_string(X.size()) +
                             ") does not match species table (" +
                             std::to_string(molecular_structures.size()) + ")");
  }
  require_finite_fractions(X, "equivalence_ratio");

  double n_C = 0.0;
  double n_H = 0.0;
  double n_O = 0.0;

  for (size_t i = 0; i < X.size(); ++i) {
    const double X_k = X[i];
    if (X_k <= 1e-15)
      continue;
    const auto &ms = molecular_structures[i];
    n_C += X_k * static_cast<double>(ms.C);
    n_H += X_k * static_cast<double>(ms.H);
    n_O += X_k * static_cast<double>(ms.O);
  }

  const double n_O_stoich = 2.0 * n_C + 0.5 * n_H;

  // Guard: Non-combustible mixture (Pure oxidizer or inert)
  if (n_O_stoich <= 1e-15) {
    return 0.0;
  }

  // Guard: Zero oxygen present (Pure fuel or fuel+inert)
  // Return a large number to indicate extreme rich condition without crashing
  if (n_O <= 1e-15) {
    return 1.0e9;
  }

  return n_O_stoich / n_O;
}

double equivalence_ratio_mass(const std::vector<double> &Y) {
  // Convert mass Y to molar X and call the molar version
  return equivalence_ratio(mass_to_mole(Y));
}

double fuel_lhv_molar(const std::vector<double> &X_fuel,
                      const double reference_temperature) {
  validate_fuel_mole_fractions(X_fuel);

  if (!std::isfinite(reference_temperature) || reference_temperature <= 0.0) {
    throw std::invalid_argument(
        "fuel_lhv_molar: reference_temperature must be finite and positive");
  }

  const std::size_t idx_O2 = species_index.at("O2");
  const std::size_t idx_CO2 = species_index.at("CO2");
  const std::size_t idx_H2O = species_index.at("H2O");
  const std::size_t idx_N2 = species_index.at("N2");

  const double h_O2 = pure_species_enthalpy(idx_O2, reference_temperature);
  const double h_CO2 = pure_species_enthalpy(idx_CO2, reference_temperature);
  const double h_H2O = pure_species_enthalpy(idx_H2O, reference_temperature);
  const double h_N2 = pure_species_enthalpy(idx_N2, reference_temperature);

  double lhv_molar = 0.0;

  for (std::size_t i = 0; i < X_fuel.size(); ++i) {
    const double x_i = X_fuel[i];
    if (x_i <= 0.0) {
      continue;
    }

    const double nu_O2 = oxygen_required_per_mol_fuel(i);
    if (nu_O2 <= 0.0) {
      continue;
    }

    const Molecular_Structure &s = molecular_structures[i];
    const double h_fuel_i = pure_species_enthalpy(i, reference_temperature);

    const double h_react = h_fuel_i + nu_O2 * h_O2;
    const double h_prod = static_cast<double>(s.C) * h_CO2 +
                          0.5 * static_cast<double>(s.H) * h_H2O +
                          0.5 * static_cast<double>(s.N) * h_N2;

    const double delta_h_rxn = h_prod - h_react; // [J/mol species i]
    lhv_molar += x_i * (-delta_h_rxn);
  }

  return std::max(0.0, lhv_molar); // [J/mol fuel]
}

double fuel_lhv_mass(const std::vector<double> &X_fuel,
                     const double reference_temperature) {
  validate_fuel_mole_fractions(X_fuel);
  const double mw_fuel = mwmix(X_fuel); // [g/mol]
  if (!validation::is_finite_positive(mw_fuel)) {
    throw std::invalid_argument(
        "fuel_lhv_mass: fuel molecular weight must be positive");
  }
  return fuel_lhv_molar(X_fuel, reference_temperature) * 1000.0 /
         mw_fuel; // [J/kg fuel]
}

// -------------------------------------------------------------
// Dry air requirements (using standard dry air composition)
// -------------------------------------------------------------

// Calculate dry air required for complete combustion [mol air/mol fuel]
double dryair_required_per_mol_fuel(std::size_t fuel_index) {
  double O2_required = oxygen_required_per_mol_fuel(fuel_index);
  double X_O2_air = dry_air_composition.at("O2");
  return O2_required / X_O2_air;
}

// Calculate dry air required for complete combustion [kg air/kg fuel]
double dryair_required_per_kg_fuel(std::size_t fuel_index) {
  if (fuel_index >= molar_masses.size()) {
    throw std::invalid_argument("Invalid fuel index");
  }

  double fuel_mw = molar_masses[fuel_index]; // g/mol
  std::vector<double> X_air = dry_air();
  double air_mw = mwmix(X_air); // g/mol

  double molar_air_required = dryair_required_per_mol_fuel(fuel_index);
  return molar_air_required * air_mw / fuel_mw;
}

// Calculate dry air required for complete combustion [mol air/mol mixture]
double dryair_required_per_mol_mixture(const std::vector<double> &X) {
  double O2_required = oxygen_required_per_mol_mixture(X);
  double X_O2_air = dry_air_composition.at("O2");
  return O2_required / X_O2_air;
}

// Calculate dry air required for complete combustion [kg air/kg mixture]
double dryair_required_per_kg_mixture(const std::vector<double> &X) {
  double mixture_mw = mwmix(X); // g/mol
  std::vector<double> X_air = dry_air();
  double air_mw = mwmix(X_air); // g/mol

  double molar_air_required = dryair_required_per_mol_mixture(X);
  return molar_air_required * air_mw / mixture_mw;
}

std::vector<double> complete_combustion_to_CO2_H2O(const std::vector<double> &X,
                                                   bool smooth_phi0,
                                                   bool smooth_phi1, double k0,
                                                   double k1) {
  double fuel_burn_fraction = 1.0;
  return complete_combustion_to_CO2_H2O(X, fuel_burn_fraction, smooth_phi0,
                                        smooth_phi1, k0, k1);
}

namespace {

// Zero-fuel smoothing (#1037). With x = k0 a (a >= 0: nu_i X_i >= 0),
//
//   g(x) = ln(1 + e^x) - ln2 sech(x),   n_req = g(x) / k0.
//
// g(0) = 0 and g'(0) = 1/2 (as the softplus ln(1 + e^x) - ln2 it replaces:
// half the fuel burns as phi -> 0), g is C-infinity and increasing, and
// 0 <= g(x) <= x. Unlike that softplus, which tends to x - ln2 and left
// ln2/k0 mol O2 per mol of mixture unburnt at every equivalence ratio,
// g(x) - x = ln(1 + e^-x) - ln2 sech(x) ~ (1 - 2 ln2) e^-x: exact to
// round-off once k0 a > ~40. The transition is about 1/k0 wide in a.
//
// Below x = 1e-3, s = g(x)/x is its series (odd in x after the 1/2, so the
// next omitted term, about 0.06 x^5, is below 1e-16): the direct quotient
// loses digits to the ln2 - ln2 cancellation there.
struct ZeroFuelSmoothing {
  double n_req = 0.0;
  double dnreq_da = 0.0;
  double s = 0.0;      // n_req / a
  double ds_da = 0.0;
};

constexpr double ZERO_FUEL_SERIES_X = 1.0e-3;
constexpr double ZERO_FUEL_C2 = 0.125 + 0.5 * M_LN2;              // g''(0) / 2
constexpr double ZERO_FUEL_C4 = -1.0 / 192.0 - 5.0 * M_LN2 / 24.0;  // g''''(0) / 24

ZeroFuelSmoothing zero_fuel_smoothing(double a, double k0) {
  ZeroFuelSmoothing z;
  const double x = k0 * std::max(a, 0.0);
  if (x < ZERO_FUEL_SERIES_X) {
    // s = 1/2 + C2 x + C4 x^3; at a = 0 these are the one-sided limits.
    z.s = 0.5 + ZERO_FUEL_C2 * x + ZERO_FUEL_C4 * x * x * x;
    const double ds_dx = ZERO_FUEL_C2 + 3.0 * ZERO_FUEL_C4 * x * x;
    z.ds_da = k0 * ds_dx;
    z.n_req = z.s * std::max(a, 0.0);
    z.dnreq_da = z.s + x * ds_dx;
    return z;
  }
  const double em = std::exp(-x);                   // x >= 1e-3: no overflow
  const double em2 = em * em;
  const double sech = 2.0 * em / (1.0 + em2);
  const double th = (1.0 - em2) / (1.0 + em2);
  const double sigma = 1.0 / (1.0 + em);           // d/dx ln(1 + e^x)
  // n_req = a + (g(x) - x) / k0, so the tail is a to round-off.
  z.n_req = std::max(0.0, a + (std::log1p(em) - M_LN2 * sech) / k0);
  z.dnreq_da = sigma + M_LN2 * sech * th;
  z.s = z.n_req / a;
  z.ds_da = (z.dnreq_da * a - z.n_req) / (a * a);
  return z;
}

// Complete combustion to CO2 and H2O on a 1 mol reactant basis, with the
// analytic Jacobian of the product mole fractions in the reactant ones
// (#716). One implementation for the value and its derivative:
//
//   a      = sum_i nu_i X_i, the O2 the fuels need (nu_i > 0 for C/H species);
//   n_req  = a, or (smooth_phi0) g(k0 a) / k0 (zero_fuel_smoothing above),
//            so the heat release is smooth at zero fuel; s = n_req / a;
//   r      = X_O2 / n_req;
//   f      = min(1, r), or (smooth_phi1) its softmin -ln(e^-k1 + e^(-k1 r))/k1
//            at every r: a transition about 1/k1 wide at Phi = 1;
//   fuels burn the fraction f s; O2 left: X_O2 - f n_req (>= 0, as f <= r),
//   which conserves the atoms inside the smoothing band too.
struct CompleteCombustionCore {
  std::vector<double> X_b;
  double burn_fraction = 0.0;
  std::vector<std::vector<double>> dXb_dX;  // [product j][reactant k]
};

CompleteCombustionCore complete_combustion_core(const std::vector<double> &X,
                                                bool smooth_phi0, bool smooth_phi1,
                                                double k0, double k1, bool with_jacobian) {
  const double tol_sum = 1.0e-5;
  const double tol = 1.0e-12;
  require_finite_fractions(X, "complete_combustion");
  const double sum = std::accumulate(X.begin(), X.end(), 0.0);
  if (std::abs(sum - 1.0) > tol_sum) {
    throw std::invalid_argument("Mole fractions must sum to 1.0");
  }
  if (X.size() != molecular_structures.size()) {
    throw std::invalid_argument("Mixture size does not match number of species");
  }
  const std::size_t ns = X.size();
  const std::size_t idx_O2 = species_index.at("O2");
  const std::size_t idx_CO2 = species_index.at("CO2");
  const std::size_t idx_H2O = species_index.at("H2O");
  const std::size_t idx_N2 = species_index.at("N2");

  CompleteCombustionCore out;
  auto unchanged = [&]() {
    out.X_b = X;
    out.burn_fraction = 0.0;
    if (with_jacobian) {
      out.dXb_dX.assign(ns, std::vector<double>(ns, 0.0));
      for (std::size_t j = 0; j < ns; ++j) {
        out.dXb_dX[j][j] = 1.0;
      }
    }
    return out;
  };

  std::vector<double> nu(ns, 0.0);
  double a = 0.0;
  for (std::size_t i = 0; i < ns; ++i) {
    const auto &sp = molecular_structures[i];
    if (sp.C > 0 || sp.H > 0) {
      nu[i] = oxygen_required_per_mol_fuel(i);
    }
    if (X[i] > 0.0) {
      a += nu[i] * X[i];
    }
  }
  const double X_O2 = X[idx_O2];
  // No oxygen, no combustion: the products are the reactants. The value stays
  // that; the Jacobian is the one-sided limit X_O2 -> 0+ below whenever a fuel
  // is present (#1251): a fuel-only mix burns the first oxygen that arrives,
  // so dX_b/dX_O2 is not the identity there. The identity slope gave a
  // combustor fed by fuel alone dT/dm_air = 0 where the difference gave
  // 2.2e6 K per kg/s. Any X_O2 > 0 burns (f <= r keeps O2 left >= 0): a cut
  // at X_O2 <= 1e-12, as before, put a kink of the same size there.
  const bool no_oxygen = !(X_O2 > 0.0);
  if (no_oxygen && (!with_jacobian || !(a > 0.0))) {
    return unchanged();
  }

  // n_req(a), its slope, s = n_req / a and ds/da.
  double n_req = a;
  double dnreq_da = 1.0;
  double s = 1.0;
  double ds_da = 0.0;
  if (smooth_phi0) {
    const ZeroFuelSmoothing z = zero_fuel_smoothing(a, k0);
    n_req = z.n_req;
    dnreq_da = z.dnreq_da;
    s = z.s;
    ds_da = z.ds_da;
  } else if (n_req <= tol) {
    return unchanged();  // no fuel
  }

  // f(r) and its slope.
  const double n_den = std::max(1e-15, n_req);
  const double r = X_O2 / n_den;
  double f = 1.0;
  double df_dr = 0.0;
  if (smooth_phi1) {
    const double x = -k1;
    const double y = -k1 * r;
    const double mx = (x > y) ? x : y;
    f = -(mx + std::log1p(std::exp(-std::abs(x - y)))) / k1;
    // df/dr = 1 / (1 + e^(k1 (r - 1))), evaluated without overflow.
    const double z = k1 * (r - 1.0);
    df_dr = (z > 0.0) ? std::exp(-z) / (1.0 + std::exp(-z)) : 1.0 / (1.0 + std::exp(z));
  } else if (r < 1.0) {
    f = r;
    df_dr = 1.0;
  }
  const double f_eff = f * s;
  // Without fuel nothing burns; f_eff stays the a -> 0+ limit for the slope.
  out.burn_fraction = (a > 0.0) ? f_eff : 0.0;

  // Products (moles, 1 mol reactant basis).
  std::vector<double> n = X;
  double d_CO2 = 0.0, d_H2O = 0.0, d_N2 = 0.0;
  for (std::size_t i = 0; i < ns; ++i) {
    if (nu[i] <= 0.0 || X[i] <= 0.0) {
      continue;
    }
    const auto &sp = molecular_structures[i];
    n[i] = X[i] * (1.0 - f_eff);
    d_CO2 += static_cast<double>(sp.C) * f_eff * X[i];
    d_H2O += 0.5 * static_cast<double>(sp.H) * f_eff * X[i];
    d_N2 += 0.5 * static_cast<double>(sp.N) * f_eff * X[i];
  }
  const double O2_left = X_O2 - f * n_req;
  n[idx_O2] = std::max(0.0, O2_left);
  n[idx_CO2] += d_CO2;
  n[idx_H2O] += d_H2O;
  n[idx_N2] += d_N2;
  double N = 0.0;
  for (double v : n) {
    N += v;
  }
  out.X_b.resize(ns);
  for (std::size_t j = 0; j < ns; ++j) {
    out.X_b[j] = n[j] / N;
  }
  if (!with_jacobian) {
    return out;
  }

  // dn/dX_k. a depends on X_k through nu_k; r through X_O2 and n_req.
  //   df/dX_k = df_dr (delta_k,O2 / n_req - r dnreq_da nu_k / n_req)
  //   df_eff/dX_k = df/dX_k s + f ds_da nu_k
  std::vector<std::vector<double>> dn(ns, std::vector<double>(ns, 0.0));
  for (std::size_t k = 0; k < ns; ++k) {
    const double dnreq = dnreq_da * nu[k];
    const double dr = ((k == idx_O2) ? 1.0 : 0.0) / n_den - (n_req > 1e-15 ? r * dnreq / n_den : 0.0);
    const double df = df_dr * dr;
    const double dfeff = df * s + f * ds_da * nu[k];
    for (std::size_t i = 0; i < ns; ++i) {
      if (nu[i] <= 0.0) {
        continue;
      }
      // n_i = X_i (1 - f_eff), and its products, for every fuel (an absent
      // one enters with the one-sided derivative).
      const double xi = std::max(X[i], 0.0);
      const double dn_i = ((i == k) ? (1.0 - f_eff) : 0.0) - xi * dfeff;
      dn[i][k] += dn_i;
      const auto &sp = molecular_structures[i];
      const double dburn = ((i == k) ? f_eff : 0.0) + xi * dfeff;
      dn[idx_CO2][k] += static_cast<double>(sp.C) * dburn;
      dn[idx_H2O][k] += 0.5 * static_cast<double>(sp.H) * dburn;
      dn[idx_N2][k] += 0.5 * static_cast<double>(sp.N) * dburn;
    }
    // Non-fuel species pass through (O2 handled below).
    if (nu[k] <= 0.0 && k != idx_O2) {
      dn[k][k] += 1.0;
    }
    // O2 left = X_O2 - f n_req.
    if (O2_left > 0.0) {
      dn[idx_O2][k] = ((k == idx_O2) ? 1.0 : 0.0) - df * n_req - f * dnreq;
    } else {
      dn[idx_O2][k] = 0.0;
    }
  }
  // X_b = n / N.
  out.dXb_dX.assign(ns, std::vector<double>(ns, 0.0));
  for (std::size_t k = 0; k < ns; ++k) {
    double dN = 0.0;
    for (std::size_t j = 0; j < ns; ++j) {
      dN += dn[j][k];
    }
    for (std::size_t j = 0; j < ns; ++j) {
      out.dXb_dX[j][k] = (dn[j][k] - out.X_b[j] * dN) / N;
    }
  }
  if (no_oxygen) {
    // The value is the unburnt mix (as without the Jacobian); the slope above
    // is that of the burning branch at r = X_O2 / n_req -> 0.
    out.X_b = X;
    out.burn_fraction = 0.0;
  }
  return out;
}

}  // namespace

std::vector<double> complete_combustion_to_CO2_H2O(const std::vector<double> &X,
                                                   double &fuel_burn_fraction,
                                                   bool smooth_phi0,
                                                   bool smooth_phi1, double k0,
                                                   double k1) {
  CompleteCombustionCore c = complete_combustion_core(X, smooth_phi0, smooth_phi1, k0, k1, false);
  fuel_burn_fraction = c.burn_fraction;
  return c.X_b;
}

CompleteCombustionJacobian complete_combustion_products_and_jacobian(const std::vector<double> &X,
                                                                     bool smooth_phi0,
                                                                     bool smooth_phi1, double k0,
                                                                     double k1) {
  CompleteCombustionCore c = complete_combustion_core(X, smooth_phi0, smooth_phi1, k0, k1, true);
  return {c.X_b, c.dXb_dX};
}

// This function is now redundant as it's merged into the other overload with default arguments,
// but for internal consistency within the .cpp we can keep it if needed or just remove it.
// I'll remove it to avoid ambiguity as the other one now has defaults.

// -------------------------------------------------------------
// Equivalence ratio helpers (mole basis, multi-species fuel)
// -------------------------------------------------------------

// We assume species_names, molar_masses and oxygen_required_per_mol_mixture(X)
// are defined in thermo_transport_data / this translation unit.

// Stoichiometric fuel-to-oxidizer mole ratio (F/O)_st for a fuel mixture
// (mole fractions X_fuel) and an oxidizer mixture X_ox (mole fractions).
//
// Uses:
//   - oxygen_required_per_mol_mixture(X_fuel): moles O2 per mole fuel mixture
//   - O2 mole fraction in oxidizer stream X_ox[O2_index]
static double stoich_f_over_o_mole(const std::vector<double> &X_fuel,
                                   const std::vector<double> &X_ox) {
  if (X_fuel.size() != species_names.size() ||
      X_ox.size() != species_names.size()) {
    throw std::invalid_argument("stoich_f_over_o_mole: size mismatch");
  }

  const double X_O2_ox = X_ox[species_index.at("O2")];
  if (!validation::is_finite_positive(X_O2_ox)) {
    throw std::invalid_argument(
        "stoich_f_over_o_mole: oxidizer has zero O2 mole fraction");
  }

  // moles O2 required per mole of fuel mixture
  const double nu_O2_req = oxygen_required_per_mol_mixture(X_fuel);
  if (!validation::is_finite_positive(nu_O2_req)) {
    throw std::invalid_argument(
        "stoich_f_over_o_mole: oxygen_required_per_mol_mixture returned "
        "non-positive");
  }

  // Let:
  //   nu_O2_req = X_O2_ox * n_ox_st  ->  n_ox_st = nu_O2_req / X_O2_ox
  // then stoich F/O (molar) = n_F / n_Ox = 1 / n_ox_st = X_O2_ox / nu_O2_req
  const double r_st = X_O2_ox / nu_O2_req;
  return r_st;
}

// Given a mixture X_mix that is formed *only* by mixing:
//   - a fuel stream with composition X_fuel
//   - an oxidizer stream with composition X_ox
// we can write:
//   X_mix = (r * X_fuel + 1 * X_ox) / (r + 1)
// where r = n_F / n_Ox is the fuel-to-oxidizer mole ratio.
//
// For any species j with X_fuel[j] != X_ox[j]:
//   X_mix[j] (r + 1) = r X_fuel[j] + X_ox[j]
//   r (X_mix[j] - X_fuel[j]) = X_ox[j] - X_mix[j]
//   r = (X_ox[j] - X_mix[j]) / (X_mix[j] - X_fuel[j])
//
// We pick the first j with |X_fuel[j] - X_ox[j]| large enough.
static double actual_f_over_o_mole(const std::vector<double> &X_mix,
                                   const std::vector<double> &X_fuel,
                                   const std::vector<double> &X_ox) {
  const std::size_t n = species_names.size();
  if (X_mix.size() != n || X_fuel.size() != n || X_ox.size() != n) {
    throw std::invalid_argument("actual_f_over_o_mole: size mismatch");
  }

  const double eps = 1e-12;
  bool found = false;
  double r = 0.0;

  for (std::size_t j = 0; j < n; ++j) {
    const double xf = X_fuel[j];
    const double xo = X_ox[j];
    const double xm = X_mix[j];

    double denom = xm - xf;
    if (std::abs(xf - xo) < eps) {
      continue; // no contrast between fuel and oxidizer for this species
    }
    if (std::abs(denom) < eps) {
      continue; // would be numerically unstable
    }

    double r_j = (xo - xm) / denom;

    // Sanity: r_j must be >= 0 for physical mixing
    if (r_j < -1e-8) {
      continue; // skip clearly unphysical projection
    }

    r = r_j;
    found = true;
    break;
  }

  if (!found) {
    throw std::invalid_argument(
        "actual_f_over_o_mole: could not determine mixing ratio "
        "from X_mix, X_fuel, X_ox (mixture not on mixing line?)");
  }

  return r;
}

double equivalence_ratio_mole(const std::vector<double> &X_mix,
                              const std::vector<double> &X_fuel,
                              const std::vector<double> &X_ox) {
  const double r_act = actual_f_over_o_mole(X_mix, X_fuel, X_ox);
  const double r_st = stoich_f_over_o_mole(X_fuel, X_ox);

  return r_act / r_st; // phi = (F/O)_act / (F/O)_st
}

std::vector<double>
set_equivalence_ratio_mole(double phi, const std::vector<double> &X_fuel,
                           const std::vector<double> &X_ox) {
  const std::size_t n = species_names.size();
  if (X_fuel.size() != n || X_ox.size() != n) {
    throw std::invalid_argument("set_equivalence_ratio_mole: size mismatch");
  }
  require_finite(phi, "set_equivalence_ratio_mole: phi");
  if (!validation::is_finite_positive(phi)) {
    throw std::invalid_argument("set_equivalence_ratio_mole: phi must be > 0");
  }

  const double r_st = stoich_f_over_o_mole(X_fuel, X_ox);
  const double r = phi * r_st; // actual F/O = phi * (F/O)_st

  const double n_F = r;
  const double n_Ox = 1.0;
  const double n_tot = n_F + n_Ox;

  std::vector<double> X_mix(n, 0.0);

  // Fuel-stream contribution
  for (std::size_t k = 0; k < n; ++k) {
    if (X_fuel[k] != 0.0) {
      X_mix[k] += (n_F * X_fuel[k]) / n_tot;
    }
  }

  // Oxidizer-stream contribution
  for (std::size_t k = 0; k < n; ++k) {
    if (X_ox[k] != 0.0) {
      X_mix[k] += (n_Ox * X_ox[k]) / n_tot;
    }
  }

  // Renormalize to guard against small numerical drifts
  double sum = std::accumulate(X_mix.begin(), X_mix.end(), 0.0);
  if (!validation::is_finite_positive(sum)) {
    throw std::invalid_argument(
        "set_equivalence_ratio_mole: resulting mixture has non-positive sum");
  }
  for (double &x : X_mix) {
    x /= sum;
  }

  return X_mix;
}

// -------------------------------------------------------------
// Stoichiometric F/O (mass basis)
// -------------------------------------------------------------

// Stoichiometric fuel-to-oxidizer mass ratio (F/O)_st for a fuel mixture
// (mass fractions Y_fuel) and an oxidizer mixture Y_ox (mass fractions).
//
// Uses:
//   - oxygen_required_per_kg_mixture(X_fuel): kg O2 per kg fuel mixture
//   - O2 mass fraction in oxidizer stream Y_ox[O2_index]
static double stoich_f_over_o_mass(const std::vector<double> &Y_fuel,
                                   const std::vector<double> &Y_ox) {
  if (Y_fuel.size() != species_names.size() ||
      Y_ox.size() != species_names.size()) {
    throw std::invalid_argument("stoich_f_over_o_mass: size mismatch");
  }

  const double Y_O2_ox = Y_ox[species_index.at("O2")];
  if (!validation::is_finite_positive(Y_O2_ox)) {
    throw std::invalid_argument(
        "stoich_f_over_o_mass: oxidizer has zero O2 mass fraction");
  }

  // Convert fuel stream from mass fractions to mole fractions for use in
  // oxygen_required_per_kg_mixture(X_fuel).
  std::vector<double> X_fuel = mass_to_mole(Y_fuel);

  // kg O2 per kg fuel mixture
  const double m_O2_req = oxygen_required_per_kg_mixture(X_fuel);
  if (!validation::is_finite_positive(m_O2_req)) {
    throw std::invalid_argument(
        "stoich_f_over_o_mass: oxygen_required_per_kg_mixture returned "
        "non-positive");
  }

  // Stoichiometric oxidizer mass (kg) per 1 kg fuel:
  //   m_O2_req = Y_O2_ox * m_ox_st  ->  m_ox_st = m_O2_req / Y_O2_ox
  // hence (F/O)_st (mass) = m_fuel / m_ox_st = 1 / (m_O2_req / Y_O2_ox)
  //                        = Y_O2_ox / m_O2_req
  const double r_st = Y_O2_ox / m_O2_req;
  return r_st;
}

// -------------------------------------------------------------
// Actual F/O (mass basis) from a given mixture Y_mix
// -------------------------------------------------------------

// For mass mixing of fuel and oxidizer:
//   Y_mix = (m_F Y_fuel + m_Ox Y_ox) / (m_F + m_Ox),
// let r = m_F / m_Ox. For any species j with Y_fuel[j] != Y_ox[j]:
//
//   Y_mix[j] (r + 1) = r Y_fuel[j] + Y_ox[j]
//   r (Y_mix[j] - Y_fuel[j]) = Y_ox[j] - Y_mix[j]
//   r = (Y_ox[j] - Y_mix[j]) / (Y_mix[j] - Y_fuel[j])
//
// We pick the first j that gives a numerically stable r.
static double actual_f_over_o_mass(const std::vector<double> &Y_mix,
                                   const std::vector<double> &Y_fuel,
                                   const std::vector<double> &Y_ox) {
  const std::size_t n = species_names.size();
  if (Y_mix.size() != n || Y_fuel.size() != n || Y_ox.size() != n) {
    throw std::invalid_argument("actual_f_over_o_mass: size mismatch");
  }

  const double eps = 1e-12;
  bool found = false;
  double r = 0.0;

  for (std::size_t j = 0; j < n; ++j) {
    const double yf = Y_fuel[j];
    const double yo = Y_ox[j];
    const double ym = Y_mix[j];

    if (std::abs(yf - yo) < eps) {
      continue; // no contrast between fuel and oxidizer for this species
    }

    const double denom = ym - yf;
    if (std::abs(denom) < eps) {
      continue; // numerically unstable
    }

    const double r_j = (yo - ym) / denom;

    // Require physically reasonable (non-negative) r_j
    if (r_j < -1e-8) {
      continue;
    }

    r = r_j;
    found = true;
    break;
  }

  if (!found) {
    throw std::invalid_argument(
        "actual_f_over_o_mass: could not determine mixing ratio "
        "from Y_mix, Y_fuel, Y_ox (mixture not on mixing line?)");
  }

  return r;
}

// -------------------------------------------------------------
// Public mass-basis phi functions
// -------------------------------------------------------------

double equivalence_ratio_mass(const std::vector<double> &Y_mix,
                              const std::vector<double> &Y_fuel,
                              const std::vector<double> &Y_ox) {
  const double r_act = actual_f_over_o_mass(Y_mix, Y_fuel, Y_ox);
  const double r_st = stoich_f_over_o_mass(Y_fuel, Y_ox);

  return r_act / r_st; // phi = (F/O)_act / (F/O)_st (mass basis)
}

std::vector<double>
set_equivalence_ratio_mass(double phi, const std::vector<double> &Y_fuel,
                           const std::vector<double> &Y_ox) {
  const std::size_t n = species_names.size();
  if (Y_fuel.size() != n || Y_ox.size() != n) {
    throw std::invalid_argument("set_equivalence_ratio_mass: size mismatch");
  }
  require_finite(phi, "set_equivalence_ratio_mass: phi");
  if (!validation::is_finite_positive(phi)) {
    throw std::invalid_argument("set_equivalence_ratio_mass: phi must be > 0");
  }

  const double r_st = stoich_f_over_o_mass(Y_fuel, Y_ox);
  const double r = phi * r_st; // (F/O)_act = phi * (F/O)_st

  const double m_F = r;
  const double m_Ox = 1.0;
  const double m_tot = m_F + m_Ox;

  std::vector<double> Y_mix(n, 0.0);

  // Fuel-stream contribution
  for (std::size_t k = 0; k < n; ++k) {
    if (Y_fuel[k] != 0.0) {
      Y_mix[k] += (m_F * Y_fuel[k]) / m_tot;
    }
  }

  // Oxidizer-stream contribution
  for (std::size_t k = 0; k < n; ++k) {
    if (Y_ox[k] != 0.0) {
      Y_mix[k] += (m_Ox * Y_ox[k]) / m_tot;
    }
  }

  // Renormalize to guard against any small numerical drift
  double sum = std::accumulate(Y_mix.begin(), Y_mix.end(), 0.0);
  if (!validation::is_finite_positive(sum)) {
    throw std::invalid_argument(
        "set_equivalence_ratio_mass: resulting mixture has non-positive sum");
  }
  for (double &y : Y_mix) {
    y /= sum;
  }

  return Y_mix;
}

double bilger_stoich_mixture_fraction_mass(const std::vector<double> &Y_F,
                                           const std::vector<double> &Y_O) {
  // r_st = (F/O)_st on a MASS basis
  // already implemented as part of equivalence_ratio_mass helpers:
  const double r_st = stoich_f_over_o_mass(Y_F, Y_O);

  if (!validation::is_finite_positive(r_st)) {
    throw std::invalid_argument("bilger_stoich_mixture_fraction_mass: r_st <= 0");
  }

  // For 2-stream fuel/oxidizer mixing we always have:
  //
  //     Z = m_F / (m_F + m_O)
  //     r = m_F / m_O
  //
  //  =>  Z = r / (1 + r)
  //
  // At stoichiometry r = r_st, giving:
  //
  //     Z_st = r_st / (1 + r_st)
  //
  return r_st / (1.0 + r_st);
}

double bilger_Z_from_equivalence_ratio_mass(double phi,
                                            const std::vector<double> &Y_F,
                                            const std::vector<double> &Y_O) {
  require_finite(phi, "bilger_Z_from_equivalence_ratio_mass: phi");
  if (!validation::is_finite_positive(phi)) {
    throw std::invalid_argument(
        "bilger_Z_from_equivalence_ratio_mass: phi must be > 0");
  }

  const double r_st = stoich_f_over_o_mass(Y_F, Y_O);
  if (!validation::is_finite_positive(r_st)) {
    throw std::invalid_argument("bilger_Z_from_equivalence_ratio_mass: r_st <= 0");
  }

  // r = (F/O)_act (mass) = phi * r_st
  const double r = phi * r_st;

  // Two-stream mixing relation on a mass basis:
  //   Z = r / (1 + r)
  return r / (1.0 + r);
}

double equivalence_ratio_from_bilger_Z_mass(double Z,
                                            const std::vector<double> &Y_F,
                                            const std::vector<double> &Y_O) {
  if (!(Z > 0.0 && Z < 1.0)) {  // also refuses a NaN Z
    throw std::invalid_argument(
        "equivalence_ratio_from_bilger_Z_mass: Z must be in (0, 1)");
  }

  const double r_st = stoich_f_over_o_mass(Y_F, Y_O);
  if (!validation::is_finite_positive(r_st)) {
    throw std::invalid_argument("equivalence_ratio_from_bilger_Z_mass: r_st <= 0");
  }

  // Invert Z = r / (1 + r) -> r = Z / (1 - Z)
  const double r = Z / (1.0 - Z);

  // phi = (F/O)_act / (F/O)_st = r / r_st
  return r / r_st;
}

// ============================
// Bilger mixture fraction
// ============================

// Bilger scalar beta(Y) based on element mass fractions of C, H, O.
// Uses:
//   molecular_structures[k].C/H/O  (atom counts per species)
//   molar_masses[k]                (species molar masses, g/mol)
//
// For each element e:
//   Z_e / W_e = sum_k Y_k * n_e,k / W_k
// so Bilger's beta can be written without explicit atomic weights:
//   beta = 2 * sum_k Y_k C_k / W_k
//        + 0.5 * sum_k Y_k H_k / W_k
//        -       sum_k Y_k O_k / W_k
double bilger_beta(const std::vector<double> &Y) {
  if (Y.size() != molar_masses.size() ||
      Y.size() != molecular_structures.size()) {
    throw std::invalid_argument("bilger_beta: mass-fraction vector size "
                                "does not match species data");
  }

  double sum_C = 0.0;
  double sum_H = 0.0;
  double sum_O = 0.0;

  for (std::size_t k = 0; k < Y.size(); ++k) {
    double Yk = Y[k];
    if (Yk == 0.0)
      continue;

    const auto &ms = molecular_structures[k];
    double Wk = molar_masses[k]; // g/mol

    if (!validation::is_finite_positive(Wk)) {
      throw std::runtime_error("bilger_beta: non-positive molar mass");
    }

    // Contributions to Z_e / W_e, see comment above
    sum_C += Yk * static_cast<double>(ms.C) / Wk;
    sum_H += Yk * static_cast<double>(ms.H) / Wk;
    sum_O += Yk * static_cast<double>(ms.O) / Wk;
  }

  // Bilger's scalar (proportional to 2 Z_C/W_C + 0.5 Z_H/W_H - Z_O/W_O)
  return 2.0 * sum_C + 0.5 * sum_H - sum_O;
}

// Bilger mixture fraction Z in [0,1] between two streams F and O,
// defined as a normalized Bilger scalar:
//   Z = (beta - beta_O) / (beta_F - beta_O)
double bilger_mixture_fraction(const std::vector<double> &Y,
                               const std::vector<double> &Y_F,
                               const std::vector<double> &Y_O) {
  if (Y.size() != Y_F.size() || Y.size() != Y_O.size()) {
    throw std::invalid_argument(
        "bilger_mixture_fraction: all mass-fraction vectors "
        "must have the same size");
  }

  const double beta = bilger_beta(Y);
  const double beta_F = bilger_beta(Y_F);
  const double beta_O = bilger_beta(Y_O);

  const double denom = beta_F - beta_O;
  if (std::abs(denom) < 1e-16) {
    throw std::invalid_argument("bilger_mixture_fraction: fuel and oxidizer have "
                             "identical Bilger scalar (denominator ~ 0)");
  }

  double Z = (beta - beta_O) / denom;

  // Clamp to [0,1] to handle small numerical overshoots
  if (Z < 0.0)
    Z = 0.0;
  else if (Z > 1.0)
    Z = 1.0;

  return Z;
}

// Convenience wrapper: inputs are mole fractions X, X_F, X_O.
// These are converted to mass fractions first and then passed to
// the mass-fraction-based Bilger implementation above.
double bilger_mixture_fraction_from_moles(const std::vector<double> &X,
                                          const std::vector<double> &X_F,
                                          const std::vector<double> &X_O) {
  if (X.size() != X_F.size() || X.size() != X_O.size()) {
    throw std::invalid_argument(
        "bilger_mixture_fraction_from_moles: all mole-fraction "
        "vectors must have the same size");
  }

  const auto Y = mole_to_mass(X);
  const auto Y_F = mole_to_mass(X_F);
  const auto Y_O = mole_to_mass(X_O);

  return bilger_mixture_fraction(Y, Y_F, Y_O);
}

// -------------------------------------------------------------
// Stream-based equivalence ratio helpers
// -------------------------------------------------------------

// Stoichiometric mass flow ratio (mdot_fuel / mdot_ox at stoich)
// Delegates to the mass-basis stoich F/O helper already in this file.
static double stoich_fuel_ox_ratio(const std::vector<double> &X_fuel,
                                   const std::vector<double> &X_ox) {
  const std::vector<double> Y_fuel = mole_to_mass(X_fuel);
  const std::vector<double> Y_ox = mole_to_mass(X_ox);
  return stoich_f_over_o_mass(Y_fuel, Y_ox);
}

Stream set_fuel_stream_for_phi(double phi, const Stream &fuel,
                               const Stream &oxidizer) {
  require_finite(phi, "set_fuel_stream_for_phi: phi");
  if (!validation::is_finite_positive(phi)) {
    throw std::invalid_argument(
        "set_fuel_stream_for_phi: phi must be positive");
  }
  if (!validation::is_finite_positive(oxidizer.mdot)) {
    throw std::invalid_argument(
        "set_fuel_stream_for_phi: oxidizer.mdot must be positive");
  }

  // mdot_fuel = phi * FAR_stoich * mdot_ox
  // Works for any fuel/oxidizer combination including blended fuels where
  // all species appear in both streams.
  const double FAR_stoich = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
  const double mdot_fuel = phi * FAR_stoich * oxidizer.mdot;

  Stream result = fuel;
  result.mdot = mdot_fuel;
  return result;
}

// -------------------------------------------------------------
// Inverse solvers for fuel/oxidizer streams (complete combustion only)
// -------------------------------------------------------------

namespace {

// Brent's method: finds x such that f(x) = target on [x_lo, x_hi].
// Combines bisection safety with secant/inverse-quadratic interpolation speed.
// Guaranteed to converge if f is continuous and target is bracketed.
// tol_f: tolerance on |f(x) - target|
// tol_x_rel: relative bracket width tolerance
template <typename Func>
double bisection_solve(Func f, double target, double x_lo, double x_hi,
                       double tol_f, std::size_t max_iter,
                       double tol_x_rel = 1e-10) {
  // Shift to root-finding: g(x) = f(x) - target
  double a = x_lo, b = x_hi;
  double fa = f(a) - target;
  double fb = f(b) - target;

  if (std::abs(fa) < tol_f)
    return a;
  if (std::abs(fb) < tol_f)
    return b;

  if (fa * fb > 0.0) {
    throw std::runtime_error(
        "bisection_solve: target outside achievable range");
  }

  // Ensure |f(b)| <= |f(a)| (b is the best guess)
  if (std::abs(fa) < std::abs(fb)) {
    std::swap(a, b);
    std::swap(fa, fb);
  }

  double c = a, fc = fa;
  bool mflag = true;
  double s = b, fs = fb;
  double d = 0.0;

  for (std::size_t iter = 0; iter < max_iter; ++iter) {
    if (std::abs(fb) < tol_f)
      return b;
    if (std::abs(b - a) < tol_x_rel * std::abs(b))
      return b;

    if (fa != fc && fb != fc) {
      // Inverse quadratic interpolation
      s = a * fb * fc / ((fa - fb) * (fa - fc)) +
          b * fa * fc / ((fb - fa) * (fb - fc)) +
          c * fa * fb / ((fc - fa) * (fc - fb));
    } else {
      // Secant method
      s = b - fb * (b - a) / (fb - fa);
    }

    // Conditions to fall back to bisection
    const double b3 = (3.0 * a + b) / 4.0;
    bool cond1 = !((s > std::min(b3, b) && s < std::max(b3, b)));
    bool cond2 = mflag && std::abs(s - b) >= std::abs(b - c) / 2.0;
    bool cond3 = !mflag && std::abs(s - b) >= std::abs(c - d) / 2.0;
    bool cond4 = mflag && std::abs(b - c) < tol_x_rel * std::abs(b);
    bool cond5 = !mflag && std::abs(c - d) < tol_x_rel * std::abs(b);

    if (cond1 || cond2 || cond3 || cond4 || cond5) {
      s = 0.5 * (a + b);
      mflag = true;
    } else {
      mflag = false;
    }

    fs = f(s) - target;
    d = c;
    c = b;
    fc = fb;

    if (fa * fs < 0.0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }

    if (std::abs(fa) < std::abs(fb)) {
      std::swap(a, b);
      std::swap(fa, fb);
    }
  }

  return b;
}

// Helper: compute Tad for given fuel mdot (oxidizer mdot fixed)
double compute_Tad_fuel(double mdot_fuel, const Stream &fuel,
                        const Stream &oxidizer) {
  Stream fuel_trial = fuel;
  fuel_trial.mdot = mdot_fuel;
  Stream mixed = mix({fuel_trial, oxidizer});
  State burned = complete_combustion(mixed.state);
  return burned.T;
}

// Helper: compute Tad for given oxidizer mdot (fuel mdot fixed)
double compute_Tad_ox(double mdot_ox, const Stream &fuel,
                      const Stream &oxidizer) {
  Stream ox_trial = oxidizer;
  ox_trial.mdot = mdot_ox;
  Stream mixed = mix({fuel, ox_trial});
  State burned = complete_combustion(mixed.state);
  return burned.T;
}

// Helper: compute species mole fraction in burned products
// vary_fuel: true = vary fuel mdot, false = vary oxidizer mdot
// dry_basis: true = remove H2O before computing mole fraction
double compute_species_burned(double mdot_var, const Stream &fuel,
                              const Stream &oxidizer,
                              const std::string &species, bool vary_fuel,
                              bool dry_basis) {
  Stream fuel_trial = fuel;
  Stream ox_trial = oxidizer;
  if (vary_fuel) {
    fuel_trial.mdot = mdot_var;
  } else {
    ox_trial.mdot = mdot_var;
  }
  Stream mixed = mix({fuel_trial, ox_trial});
  std::vector<double> X_burned = complete_combustion_to_CO2_H2O(mixed.state.X, false, false);
  if (dry_basis) {
    X_burned = convert_to_dry_fractions(X_burned);
  }
  return X_burned[species_index.at(species)];
}

} // anonymous namespace

// -------------------------------------------------------------
// Find fuel stream for target property (oxidizer mdot fixed)
// -------------------------------------------------------------

Stream set_fuel_stream_for_Tad(double T_ad_target, const Stream &fuel,
                               const Stream &oxidizer, double tol,
                               std::size_t max_iter, bool lean,
                               double phi_max) {
  require_finite(T_ad_target, "set_fuel_stream_for_Tad: T_ad_target");
  if (!validation::is_finite_positive(oxidizer.mdot)) {
    throw std::invalid_argument(
        "set_fuel_stream_for_Tad: oxidizer.mdot must be positive");
  }
  if (phi_max <= 1.0) {
    throw std::invalid_argument(
        "set_fuel_stream_for_Tad: phi_max must be > 1.0");
  }

  // Tad must be above oxidizer temperature (minimum achievable with zero fuel)
  double T_min = oxidizer.T();
  if (T_ad_target <= T_min) {
    throw std::invalid_argument(
        "set_fuel_stream_for_Tad: T_ad_target must be > oxidizer temperature");
  }
  if (T_ad_target > 5000.0) {
    throw std::invalid_argument(
        "set_fuel_stream_for_Tad: T_ad_target must be <= 5000 K");
  }

  // Estimate stoichiometric fuel mass flow
  double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
  double mdot_stoich = oxidizer.mdot * ratio;

  // Peak Tad is at stoichiometric
  double mdot_peak = mdot_stoich;
  double T_peak = compute_Tad_fuel(mdot_peak, fuel, oxidizer);

  if (T_ad_target > T_peak + tol) {
    throw std::runtime_error("set_fuel_stream_for_Tad: target T_ad exceeds "
                             "maximum achievable (stoichiometric)");
  }

  double mdot_result = 0.0;
  auto f = [&](double mdot) { return compute_Tad_fuel(mdot, fuel, oxidizer); };

  if (lean) {
    // Search on lean side: low fuel to stoichiometric
    // Tad increases with fuel on lean side
    double mdot_lo = mdot_stoich * 0.01;
    double T_lo = compute_Tad_fuel(mdot_lo, fuel, oxidizer);

    if (T_ad_target < T_lo - tol) {
      throw std::runtime_error("set_fuel_stream_for_Tad: target T_ad below "
                               "minimum achievable (lean)");
    }

    mdot_result =
        bisection_solve(f, T_ad_target, mdot_lo, mdot_peak, tol, max_iter);
  } else {
    // Search on rich side: stoichiometric to phi_max
    // Tad decreases with fuel on rich side
    double mdot_hi = mdot_stoich * phi_max;
    double T_hi = compute_Tad_fuel(mdot_hi, fuel, oxidizer);

    if (T_ad_target < T_hi - tol) {
      throw std::runtime_error("set_fuel_stream_for_Tad: target T_ad below "
                               "minimum achievable (rich)");
    }

    mdot_result =
        bisection_solve(f, T_ad_target, mdot_peak, mdot_hi, tol, max_iter);
  }

  Stream result = fuel;
  result.mdot = mdot_result;
  return result;
}

Stream set_fuel_stream_for_O2(double X_O2_target, const Stream &fuel,
                              const Stream &oxidizer, double tol,
                              std::size_t max_iter) {
  require_finite(X_O2_target, "set_fuel_stream_for_O2: X_O2_target");
  if (!validation::is_finite_positive(oxidizer.mdot)) {
    throw std::invalid_argument(
        "set_fuel_stream_for_O2: oxidizer.mdot must be positive");
  }

  double X_O2_ox = oxidizer.X()[species_index.at("O2")];
  if (!validation::is_finite_positive(X_O2_target) || X_O2_target >= X_O2_ox) {
    throw std::invalid_argument(
        "set_fuel_stream_for_O2: X_O2_target must be in (0, X_O2_oxidizer)");
  }

  double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
  double mdot_stoich = oxidizer.mdot * ratio;

  double mdot_lo = mdot_stoich * 0.001;
  double mdot_hi = mdot_stoich * 0.999;

  auto f = [&](double mdot) {
    return compute_species_burned(mdot, fuel, oxidizer, "O2", true, false);
  };
  double mdot_result =
      bisection_solve(f, X_O2_target, mdot_lo, mdot_hi, tol, max_iter);

  Stream result = fuel;
  result.mdot = mdot_result;
  return result;
}

Stream set_fuel_stream_for_O2_dry(double X_O2_dry_target, const Stream &fuel,
                                  const Stream &oxidizer, double tol,
                                  std::size_t max_iter) {
  require_finite(X_O2_dry_target, "set_fuel_stream_for_O2_dry: X_O2_dry_target");
  if (!validation::is_finite_positive(oxidizer.mdot)) {
    throw std::invalid_argument(
        "set_fuel_stream_for_O2_dry: oxidizer.mdot must be positive");
  }

  // Dry O2 upper bound is from dry oxidizer composition
  std::vector<double> X_ox_dry = convert_to_dry_fractions(oxidizer.X());
  double X_O2_ox_dry = X_ox_dry[species_index.at("O2")];
  if (!validation::is_finite_positive(X_O2_dry_target) || X_O2_dry_target >= X_O2_ox_dry) {
    throw std::invalid_argument("set_fuel_stream_for_O2_dry: X_O2_dry_target "
                                "must be in (0, X_O2_dry_oxidizer)");
  }

  double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
  double mdot_stoich = oxidizer.mdot * ratio;

  double mdot_lo = mdot_stoich * 0.001;
  double mdot_hi = mdot_stoich * 0.999;

  auto f = [&](double mdot) {
    return compute_species_burned(mdot, fuel, oxidizer, "O2", true, true);
  };
  double mdot_result =
      bisection_solve(f, X_O2_dry_target, mdot_lo, mdot_hi, tol, max_iter);

  Stream result = fuel;
  result.mdot = mdot_result;
  return result;
}

Stream set_fuel_stream_for_CO2(double X_CO2_target, const Stream &fuel,
                               const Stream &oxidizer, double tol,
                               std::size_t max_iter) {
  require_finite(X_CO2_target, "set_fuel_stream_for_CO2: X_CO2_target");
  if (!validation::is_finite_positive(oxidizer.mdot)) {
    throw std::invalid_argument(
        "set_fuel_stream_for_CO2: oxidizer.mdot must be positive");
  }
  if (!validation::is_finite_positive(X_CO2_target)) {
    throw std::invalid_argument(
        "set_fuel_stream_for_CO2: X_CO2_target must be positive");
  }

  double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
  double mdot_stoich = oxidizer.mdot * ratio;

  double mdot_lo = mdot_stoich * 0.001;
  double mdot_hi = mdot_stoich * 0.999;

  auto f = [&](double mdot) {
    return compute_species_burned(mdot, fuel, oxidizer, "CO2", true, false);
  };
  double mdot_result =
      bisection_solve(f, X_CO2_target, mdot_lo, mdot_hi, tol, max_iter);

  Stream result = fuel;
  result.mdot = mdot_result;
  return result;
}

Stream set_fuel_stream_for_CO2_dry(double X_CO2_dry_target, const Stream &fuel,
                                   const Stream &oxidizer, double tol,
                                   std::size_t max_iter) {
  require_finite(X_CO2_dry_target, "set_fuel_stream_for_CO2_dry: X_CO2_dry_target");
  if (!validation::is_finite_positive(oxidizer.mdot)) {
    throw std::invalid_argument(
        "set_fuel_stream_for_CO2_dry: oxidizer.mdot must be positive");
  }
  if (!validation::is_finite_positive(X_CO2_dry_target)) {
    throw std::invalid_argument(
        "set_fuel_stream_for_CO2_dry: X_CO2_dry_target must be positive");
  }

  double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
  double mdot_stoich = oxidizer.mdot * ratio;

  double mdot_lo = mdot_stoich * 0.001;
  double mdot_hi = mdot_stoich * 0.999;

  auto f = [&](double mdot) {
    return compute_species_burned(mdot, fuel, oxidizer, "CO2", true, true);
  };
  double mdot_result =
      bisection_solve(f, X_CO2_dry_target, mdot_lo, mdot_hi, tol, max_iter);

  Stream result = fuel;
  result.mdot = mdot_result;
  return result;
}

// -------------------------------------------------------------
// Find oxidizer stream for target property (fuel mdot fixed)
// -------------------------------------------------------------

Stream set_oxidizer_stream_for_Tad(double T_ad_target, const Stream &fuel,
                                   const Stream &oxidizer, double tol,
                                   std::size_t max_iter, bool lean,
                                   double phi_max) {
  require_finite(T_ad_target, "set_oxidizer_stream_for_Tad: T_ad_target");
  if (!validation::is_finite_positive(fuel.mdot)) {
    throw std::invalid_argument(
        "set_oxidizer_stream_for_Tad: fuel.mdot must be positive");
  }
  if (phi_max <= 1.0) {
    throw std::invalid_argument(
        "set_oxidizer_stream_for_Tad: phi_max must be > 1.0");
  }

  // Tad must be above fuel temperature (minimum achievable with infinite air)
  double T_min = fuel.T();
  if (T_ad_target <= T_min) {
    throw std::invalid_argument(
        "set_oxidizer_stream_for_Tad: T_ad_target must be > fuel temperature");
  }
  if (T_ad_target > 5000.0) {
    throw std::invalid_argument(
        "set_oxidizer_stream_for_Tad: T_ad_target must be <= 5000 K");
  }

  // Estimate stoichiometric oxidizer mass flow
  double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
  double mdot_stoich = (ratio > 0.0) ? fuel.mdot / ratio : fuel.mdot * 20.0;

  // Peak Tad is at stoichiometric
  double mdot_peak = mdot_stoich;
  double T_peak = compute_Tad_ox(mdot_peak, fuel, oxidizer);

  if (T_ad_target > T_peak + tol) {
    throw std::runtime_error("set_oxidizer_stream_for_Tad: target T_ad exceeds "
                             "maximum achievable (stoichiometric)");
  }

  double mdot_result = 0.0;
  auto f = [&](double mdot) { return compute_Tad_ox(mdot, fuel, oxidizer); };

  if (lean) {
    // Search on lean side: stoichiometric to high oxidizer (1/phi_min factor)
    // Tad decreases with oxidizer on lean side
    double mdot_hi = mdot_stoich * phi_max; // phi_max used as dilution factor
    double T_hi = compute_Tad_ox(mdot_hi, fuel, oxidizer);

    if (T_ad_target < T_hi - tol) {
      throw std::runtime_error("set_oxidizer_stream_for_Tad: target T_ad below "
                               "minimum achievable (lean)");
    }

    mdot_result =
        bisection_solve(f, T_ad_target, mdot_peak, mdot_hi, tol, max_iter);
  } else {
    // Search on rich side: low oxidizer (1/phi_max) to stoichiometric
    // Tad increases with oxidizer on rich side
    double mdot_lo = mdot_stoich / phi_max;
    double T_lo = compute_Tad_ox(mdot_lo, fuel, oxidizer);

    if (T_ad_target < T_lo - tol) {
      throw std::runtime_error("set_oxidizer_stream_for_Tad: target T_ad below "
                               "minimum achievable (rich)");
    }

    mdot_result =
        bisection_solve(f, T_ad_target, mdot_lo, mdot_peak, tol, max_iter);
  }

  Stream result = oxidizer;
  result.mdot = mdot_result;
  return result;
}

Stream set_oxidizer_stream_for_O2(double X_O2_target, const Stream &fuel,
                                  const Stream &oxidizer, double tol,
                                  std::size_t max_iter) {
  require_finite(X_O2_target, "set_oxidizer_stream_for_O2: X_O2_target");
  if (!validation::is_finite_positive(fuel.mdot)) {
    throw std::invalid_argument(
        "set_oxidizer_stream_for_O2: fuel.mdot must be positive");
  }

  double X_O2_ox = oxidizer.X()[species_index.at("O2")];
  if (!validation::is_finite_positive(X_O2_target) || X_O2_target >= X_O2_ox) {
    throw std::invalid_argument("set_oxidizer_stream_for_O2: X_O2_target must "
                                "be in (0, X_O2_oxidizer)");
  }

  double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
  double mdot_stoich = (ratio > 0.0) ? fuel.mdot / ratio : fuel.mdot * 20.0;

  double mdot_lo = mdot_stoich * 1.001; // just above stoich -> low O2
  double mdot_hi = mdot_stoich * 100.0; // very lean -> high O2

  // O2 increases with increasing oxidizer (on lean side)
  auto f = [&](double mdot) {
    return compute_species_burned(mdot, fuel, oxidizer, "O2", false, false);
  };
  double mdot_result =
      bisection_solve(f, X_O2_target, mdot_lo, mdot_hi, tol, max_iter);

  Stream result = oxidizer;
  result.mdot = mdot_result;
  return result;
}

Stream set_oxidizer_stream_for_O2_dry(double X_O2_dry_target,
                                      const Stream &fuel,
                                      const Stream &oxidizer, double tol,
                                      std::size_t max_iter) {
  require_finite(X_O2_dry_target, "set_oxidizer_stream_for_O2_dry: X_O2_dry_target");
  if (!validation::is_finite_positive(fuel.mdot)) {
    throw std::invalid_argument(
        "set_oxidizer_stream_for_O2_dry: fuel.mdot must be positive");
  }

  std::vector<double> X_ox_dry = convert_to_dry_fractions(oxidizer.X());
  double X_O2_ox_dry = X_ox_dry[species_index.at("O2")];
  if (!validation::is_finite_positive(X_O2_dry_target) || X_O2_dry_target >= X_O2_ox_dry) {
    throw std::invalid_argument(
        "set_oxidizer_stream_for_O2_dry: X_O2_dry_target must be in (0, "
        "X_O2_dry_oxidizer)");
  }

  double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
  double mdot_stoich = (ratio > 0.0) ? fuel.mdot / ratio : fuel.mdot * 20.0;

  double mdot_lo = mdot_stoich * 1.001;
  double mdot_hi = mdot_stoich * 100.0;

  auto f = [&](double mdot) {
    return compute_species_burned(mdot, fuel, oxidizer, "O2", false, true);
  };
  double mdot_result =
      bisection_solve(f, X_O2_dry_target, mdot_lo, mdot_hi, tol, max_iter);

  Stream result = oxidizer;
  result.mdot = mdot_result;
  return result;
}

Stream set_oxidizer_stream_for_CO2(double X_CO2_target, const Stream &fuel,
                                   const Stream &oxidizer, double tol,
                                   std::size_t max_iter) {
  require_finite(X_CO2_target, "set_oxidizer_stream_for_CO2: X_CO2_target");
  if (!validation::is_finite_positive(fuel.mdot)) {
    throw std::invalid_argument(
        "set_oxidizer_stream_for_CO2: fuel.mdot must be positive");
  }
  if (!validation::is_finite_positive(X_CO2_target)) {
    throw std::invalid_argument(
        "set_oxidizer_stream_for_CO2: X_CO2_target must be positive");
  }

  double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
  double mdot_stoich = (ratio > 0.0) ? fuel.mdot / ratio : fuel.mdot * 20.0;

  double mdot_lo = mdot_stoich * 1.001; // just above stoich -> high CO2
  double mdot_hi = mdot_stoich * 100.0; // very lean -> low CO2

  // CO2 decreases with increasing oxidizer (dilution)
  auto f = [&](double mdot) {
    return compute_species_burned(mdot, fuel, oxidizer, "CO2", false, false);
  };
  double mdot_result =
      bisection_solve(f, X_CO2_target, mdot_lo, mdot_hi, tol, max_iter);

  Stream result = oxidizer;
  result.mdot = mdot_result;
  return result;
}

Stream set_oxidizer_stream_for_CO2_dry(double X_CO2_dry_target,
                                       const Stream &fuel,
                                       const Stream &oxidizer, double tol,
                                       std::size_t max_iter) {
  require_finite(X_CO2_dry_target, "set_oxidizer_stream_for_CO2_dry: X_CO2_dry_target");
  if (!validation::is_finite_positive(fuel.mdot)) {
    throw std::invalid_argument(
        "set_oxidizer_stream_for_CO2_dry: fuel.mdot must be positive");
  }
  if (!validation::is_finite_positive(X_CO2_dry_target)) {
    throw std::invalid_argument(
        "set_oxidizer_stream_for_CO2_dry: X_CO2_dry_target must be positive");
  }

  double ratio = stoich_fuel_ox_ratio(fuel.X(), oxidizer.X());
  double mdot_stoich = (ratio > 0.0) ? fuel.mdot / ratio : fuel.mdot * 20.0;

  double mdot_lo = mdot_stoich * 1.001;
  double mdot_hi = mdot_stoich * 100.0;

  auto f = [&](double mdot) {
    return compute_species_burned(mdot, fuel, oxidizer, "CO2", false, true);
  };
  double mdot_result =
      bisection_solve(f, X_CO2_dry_target, mdot_lo, mdot_hi, tol, max_iter);

  Stream result = oxidizer;
  result.mdot = mdot_result;
  return result;
}

// -------------------------------------------------------------
// State-based combustion functions
// -------------------------------------------------------------

// Complete combustion (isothermal) - keeps input temperature
State complete_combustion_isothermal(const State &in, bool smooth_phi0,
                                     bool smooth_phi1, double k0, double k1) {
  require_physical_state(in, "complete_combustion_isothermal");
  State out;
  out.T = in.T;
  out.P = in.P;
  out.X = complete_combustion_to_CO2_H2O(in.X, smooth_phi0, smooth_phi1, k0, k1);
  out.set_X(out.X); // Ensure Y is populated
  return out;
}

namespace {

// Adiabatic complete combustion; burn_fraction receives the share of the
// fuel burnt (min(1, 1/phi) unsmoothed).
State complete_combustion_adiabatic(const State &in, bool smooth_phi0,
                                    bool smooth_phi1, double k0, double k1,
                                    double &burn_fraction) {
  require_physical_state(in, "complete_combustion");
  // Get burned composition at constant T first
  std::vector<double> X_burned = complete_combustion_to_CO2_H2O(
      in.X, burn_fraction, smooth_phi0, smooth_phi1, k0, k1);

  // Calculate inlet enthalpy
  double H_in = h(in.T, in.X);

  // Convert H_in (J/mol_reactants) to J/mol_products to strictly conserve
  // energy on a mass basis since chemical reactions do not conserve moles.
  double H_target = H_in * (mwmix(X_burned) / mwmix(in.X));

  // Solve for adiabatic flame temperature: h(T_ad, X_burned) = H_target
  double T_ad = calc_T_from_h(H_target, X_burned, in.T);

  State out;
  out.set_TPX(T_ad, in.P, X_burned);
  return out;
}

}  // namespace

// Complete combustion (adiabatic) - solves for flame temperature
State complete_combustion(const State &in, bool smooth_phi0, bool smooth_phi1,
                         double k0, double k1) {
  double burn_fraction = 0.0;
  return complete_combustion_adiabatic(in, smooth_phi0, smooth_phi1, k0, k1,
                                       burn_fraction);
}

// -------------------------------------------------------------
// Combustion State Dataclass
// -------------------------------------------------------------

// Compute Bilger mixture fraction from mole fraction compositions.
// Delegates to the correct mass-fraction-based implementation.
static double compute_bilger_mixture_fraction(const std::vector<double> &X_mix,
                                              const std::vector<double> &X_fuel,
                                              const std::vector<double> &X_ox) {
  return bilger_mixture_fraction_from_moles(X_mix, X_fuel, X_ox);
}

// Variant 1: From equivalence ratio (phi as input)
CombustionState combustion_state(const std::vector<double> &X_fuel,
                                 const std::vector<double> &X_ox, double phi,
                                 double T_reactants, double P,
                                 const std::string &fuel_name,
                                 CombustionMethod method,
                                 bool smooth_phi0, bool smooth_phi1,
                                 double k0, double k1) {
  CombustionState result;

  // Store phi, fuel name, and method
  result.phi = phi;
  result.fuel_name = fuel_name;
  result.method = method;

  // Compute reactant mixture composition
  std::vector<double> X_reactants =
      set_equivalence_ratio_mole(phi, X_fuel, X_ox);

  // Compute reactant CompleteState
  result.reactants = complete_state(T_reactants, P, X_reactants);

  // Compute products via selected method
  State reactant_state;
  reactant_state.set_TPX(T_reactants, P, X_reactants);
  State product_state;
  if (method == CombustionMethod::Equilibrium) {
    EquilibriumResult eq = combustion_equilibrium(reactant_state);
    product_state = eq.state;
    result.fuel_burn_fraction = eq.converged ? 1.0 : 0.0;
  } else {
    // Rich, complete combustion leaves 1 - 1/phi of the fuel unburnt.
    product_state = complete_combustion_adiabatic(reactant_state, smooth_phi0, smooth_phi1, k0,
                                                  k1, result.fuel_burn_fraction);
  }

  // Compute product CompleteState
  result.products =
      complete_state(product_state.T, product_state.P, product_state.X);

  // Compute mixture fraction
  result.mixture_fraction =
      compute_bilger_mixture_fraction(X_reactants, X_fuel, X_ox);

  return result;
}

namespace {

// The reactant state two streams make: T and P of their adiabatic mix()
// (P the lowest inlet pressure, the products' pressure too) and phi from
// the mass flows, phi = (mdot_fuel / mdot_ox) / FAR_stoich.
struct StreamsInlet {
  double T;
  double P;
  double phi;
};

StreamsInlet streams_inlet(const Stream &fuel_stream, const Stream &ox_stream) {
  if (!std::isfinite(ox_stream.mdot) || ox_stream.mdot <= 0.0) {
    throw std::invalid_argument(
        "combustion_state_from_streams: ox_stream.mdot must be positive");
  }
  if (!std::isfinite(fuel_stream.mdot) || fuel_stream.mdot < 0.0) {
    throw std::invalid_argument(
        "combustion_state_from_streams: fuel_stream.mdot must be >= 0");
  }
  const Stream mixed = mix({fuel_stream, ox_stream});
  const double FAR_stoich =
      stoich_fuel_ox_ratio(fuel_stream.X(), ox_stream.X());
  return {mixed.state.T, mixed.state.P,
          (fuel_stream.mdot / ox_stream.mdot) / FAR_stoich};
}

}  // namespace

// Variant 2: From streams (phi computed from mass flows)
CombustionState combustion_state_from_streams(const Stream &fuel_stream,
                                               const Stream &ox_stream,
                                               const std::string &fuel_name,
                                               CombustionMethod method,
                                               bool smooth_phi0, bool smooth_phi1,
                                               double k0, double k1) {
  const StreamsInlet in = streams_inlet(fuel_stream, ox_stream);
  return combustion_state(fuel_stream.X(), ox_stream.X(), in.phi, in.T, in.P,
                          fuel_name, method, smooth_phi0, smooth_phi1, k0, k1);
}

// -------------------------------------------------------------
// Overloads with user-supplied PressureLossCorrelation
// -------------------------------------------------------------

// Variant 1 with hook: from equivalence ratio
CombustionState combustion_state(const std::vector<double> &X_fuel,
                                 const std::vector<double> &X_ox, double phi,
                                 double T_reactants, double P,
                                 const std::string &fuel_name,
                                 CombustionMethod method,
                                 const PressureLossCorrelation &pressure_loss,
                                 bool smooth_phi0, bool smooth_phi1,
                                 double k0, double k1) {
  // Compute reactant mixture and product state (at inlet P)
  const std::vector<double> X_reactants =
      set_equivalence_ratio_mole(phi, X_fuel, X_ox);
  State reactant_state;
  reactant_state.set_TPX(T_reactants, P, X_reactants);

  State product_state;
  double fuel_burn_fraction = 1.0;
  if (method == CombustionMethod::Equilibrium) {
    EquilibriumResult eq = combustion_equilibrium(reactant_state);
    product_state = eq.state;
    fuel_burn_fraction = eq.converged ? 1.0 : 0.0;
  } else {
    product_state = complete_combustion_adiabatic(reactant_state, smooth_phi0, smooth_phi1, k0,
                                                  k1, fuel_burn_fraction);
  }

  // Build context -- all fields are loop-free
  const double T_ad = product_state.T;
  const double theta = (T_reactants > 0.0) ? (T_ad / T_reactants - 1.0) : 0.0;
  PressureLossContext ctx{
      reactant_state,
      phi,
      T_ad,
      product_state.X,
      product_state.Y,
      theta,
      0.0, // mdot_fuel not available in phi-based call
      0.0, // mdot_air not available in phi-based call
      0.0, // mdot_total not available in phi-based call
  };
  const double delta_P_frac = std::get<0>(pressure_loss(ctx));
  const double P_out = P * (1.0 - delta_P_frac);

  // Assemble result
  CombustionState result;
  result.phi = phi;
  result.fuel_name = fuel_name;
  result.method = method;
  result.fuel_burn_fraction = fuel_burn_fraction;
  result.reactants = complete_state(T_reactants, P, X_reactants);
  result.products = complete_state(product_state.T, P_out, product_state.X);
  result.mixture_fraction =
      compute_bilger_mixture_fraction(X_reactants, X_fuel, X_ox);
  return result;
}

// Variant 2 with hook: from streams
CombustionState combustion_state_from_streams(
    const Stream &fuel_stream, const Stream &ox_stream,
    const std::string &fuel_name, CombustionMethod method,
    const PressureLossCorrelation &pressure_loss, bool smooth_phi0,
    bool smooth_phi1, double k0, double k1) {
  const StreamsInlet in = streams_inlet(fuel_stream, ox_stream);
  const double T_reactants = in.T;
  const double P = in.P;
  const double phi = in.phi;

  // Compute products at inlet P to build the context
  const std::vector<double> X_reactants =
      set_equivalence_ratio_mole(phi, fuel_stream.X(), ox_stream.X());
  State reactant_state;
  reactant_state.set_TPX(T_reactants, P, X_reactants);

  State product_state;
  double fuel_burn_fraction = 1.0;
  if (method == CombustionMethod::Equilibrium) {
    EquilibriumResult eq = combustion_equilibrium(reactant_state);
    product_state = eq.state;
    fuel_burn_fraction = eq.converged ? 1.0 : 0.0;
  } else {
    product_state = complete_combustion_adiabatic(reactant_state, smooth_phi0, smooth_phi1, k0,
                                                  k1, fuel_burn_fraction);
  }

  const double T_ad = product_state.T;
  const double theta = (T_reactants > 0.0) ? (T_ad / T_reactants - 1.0) : 0.0;
  PressureLossContext ctx{reactant_state,
                          phi,
                          T_ad,
                          product_state.X,
                          product_state.Y,
                          theta,
                          fuel_stream.mdot,
                          ox_stream.mdot,
                          fuel_stream.mdot + ox_stream.mdot};
  const double delta_P_frac = std::get<0>(pressure_loss(ctx));
  const double P_out = P * (1.0 - delta_P_frac);

  CombustionState result;
  result.phi = phi;
  result.fuel_name = fuel_name;
  result.method = method;
  result.fuel_burn_fraction = fuel_burn_fraction;
  result.reactants = complete_state(T_reactants, P, X_reactants);
  result.products = complete_state(product_state.T, P_out, product_state.X);
  result.mixture_fraction = compute_bilger_mixture_fraction(
      X_reactants, fuel_stream.X(), ox_stream.X());
  return result;
}
} // namespace combaero
