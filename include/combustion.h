#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "state.h"

namespace combaero {

// -----------------------------------------------------------------
// Global smoothing constants for combustion (see combustion.cpp for rationale).
// k0: the zero-fuel smoothing (smooth_phi0) is about 1/k0 wide in the O2
// demand a [mol O2/mol mixture] and exact once k0 a > ~40 (#1037).
// k1: the softmin at Phi = 1 (smooth_phi1) is about 1/k1 wide in X_O2/n_req.
// -----------------------------------------------------------------
constexpr double SMOOTHING_K_PHI0 = 20000.0;
constexpr double SMOOTHING_K_PHI1 = 20000.0;

// -------------------------------------------------------------
// User-supplied pressure-loss correlation hook
// -------------------------------------------------------------
struct PressureLossContext {
  const State &state_in;                 // inlet T, P, X (static conditions)
  double phi;                            // equivalence ratio [-]
  double T_ad;                           // adiabatic flame temperature [K]
  const std::vector<double> &X_products; // burned gas mole fractions
  const std::vector<double> &Y_products; // burned gas mass fractions
  double theta;      // T_ad/T_in - 1  (dimensionless temperature rise)
  double mdot_fuel;  // fuel mass flow [kg/s]  (0 if phi-based call)
  double mdot_air;   // air mass flow [kg/s]   (0 if phi-based call)
  double mdot_total; // total flow of the arriving streams [kg/s] (#1239)
};

using PressureLossCorrelation =
    std::function<std::tuple<double, double>(const PressureLossContext &)>;

// Combustion calculations - oxygen requirements
double oxygen_required_per_mol_fuel(std::size_t fuel_index);
double oxygen_required_per_kg_fuel(std::size_t fuel_index);
double oxygen_required_per_mol_mixture(const std::vector<double> &X);
double oxygen_required_per_kg_mixture(const std::vector<double> &X);

// Equivalence ratio based on elemental analysis (Phi = O2_needed / O2_available)
double equivalence_ratio(const std::vector<double> &X);
double equivalence_ratio_mass(const std::vector<double> &Y);

// Fuel lower heating value (LHV) from complete combustion to CO2 + H2O(g)
double fuel_lhv_molar(const std::vector<double> &X_fuel,
                      const double reference_temperature = 298.15);
double fuel_lhv_mass(const std::vector<double> &X_fuel,
                     const double reference_temperature = 298.15);

// Combustion calculations - dry air requirements
double dryair_required_per_mol_fuel(std::size_t fuel_index);
double dryair_required_per_kg_fuel(std::size_t fuel_index);
double dryair_required_per_mol_mixture(const std::vector<double> &X);
double dryair_required_per_kg_mixture(const std::vector<double> &X);

// Equivalence ratio (mole basis)
double equivalence_ratio_mole(const std::vector<double> &X_mix,
                               const std::vector<double> &X_fuel,
                               const std::vector<double> &X_ox);

std::vector<double>
set_equivalence_ratio_mole(double phi, const std::vector<double> &X_fuel,
                           const std::vector<double> &X_ox);

// Equivalence ratio (mass basis)
double equivalence_ratio_mass(const std::vector<double> &Y_mix,
                               const std::vector<double> &Y_fuel,
                               const std::vector<double> &Y_ox);

std::vector<double>
set_equivalence_ratio_mass(double phi, const std::vector<double> &Y_fuel,
                           const std::vector<double> &Y_ox);

// Stoichiometric Bilger mixture fraction
double bilger_stoich_mixture_fraction_mass(const std::vector<double> &Y_F,
                                            const std::vector<double> &Y_O);

double equivalence_ratio_from_bilger_Z_mass(double Z,
                                             const std::vector<double> &Y_F,
                                             const std::vector<double> &Y_O);

double bilger_Z_from_equivalence_ratio_mass(double phi,
                                             const std::vector<double> &Y_F,
                                             const std::vector<double> &Y_O);

// Complete combustion to CO2 and H2O.
std::vector<double> complete_combustion_to_CO2_H2O(const std::vector<double> &X,
                                                   bool smooth_phi0 = false,
                                                   bool smooth_phi1 = false,
                                                   double k0 = SMOOTHING_K_PHI0,
                                                   double k1 = SMOOTHING_K_PHI1);

std::vector<double> complete_combustion_to_CO2_H2O(const std::vector<double> &X,
                                                   double &fuel_burn_fraction,
                                                   bool smooth_phi0 = false,
                                                   bool smooth_phi1 = false,
                                                   double k0 = SMOOTHING_K_PHI0,
                                                   double k1 = SMOOTHING_K_PHI1);

// The same products with the analytic Jacobian of their mole fractions in the
// reactant ones, dXb_dX[j][k] = d X_b,j / d X_k (X taken unnormalised: chain it
// with a dX that keeps sum X = 1). One implementation with the value above
// (#716): the softmin at Phi = 1 is applied at every r, and the O2 left is
// X_O2 - f n_req, so the smoothing band conserves atoms and is smooth.
struct CompleteCombustionJacobian {
    std::vector<double> X_b;
    std::vector<std::vector<double>> dXb_dX;  // [product j][reactant k]
};
CompleteCombustionJacobian complete_combustion_products_and_jacobian(
    const std::vector<double> &X, bool smooth_phi0 = true, bool smooth_phi1 = true,
    double k0 = SMOOTHING_K_PHI0, double k1 = SMOOTHING_K_PHI1);

// Mixture fraction (Bilger) utilities
double bilger_beta(const std::vector<double> &Y);
double bilger_mixture_fraction(const std::vector<double> &Y,
                                const std::vector<double> &Y_F,
                                const std::vector<double> &Y_O);

double bilger_mixture_fraction_from_moles(const std::vector<double> &X,
                                           const std::vector<double> &X_F,
                                           const std::vector<double> &X_O);

// Stream-based equivalence ratio helpers
Stream set_fuel_stream_for_phi(double phi, const Stream &fuel,
                                const Stream &oxidizer);

// Inverse solvers for fuel/oxidizer streams
Stream set_fuel_stream_for_Tad(double T_ad_target, const Stream &fuel,
                                const Stream &oxidizer, double tol = 1.0,
                                std::size_t max_iter = 100, bool lean = true,
                                double phi_max = 10.0);

Stream set_fuel_stream_for_O2(double X_O2_target, const Stream &fuel,
                               const Stream &oxidizer, double tol = 1e-6,
                               std::size_t max_iter = 100);

Stream set_fuel_stream_for_O2_dry(double X_O2_dry_target, const Stream &fuel,
                                   const Stream &oxidizer, double tol = 1e-6,
                                   std::size_t max_iter = 100);

Stream set_fuel_stream_for_CO2(double X_CO2_target, const Stream &fuel,
                                const Stream &oxidizer, double tol = 1e-6,
                                std::size_t max_iter = 100);

Stream set_fuel_stream_for_CO2_dry(double X_CO2_dry_target, const Stream &fuel,
                                    const Stream &oxidizer, double tol = 1e-6,
                                    std::size_t max_iter = 100);

Stream set_oxidizer_stream_for_Tad(double T_ad_target, const Stream &fuel,
                                    const Stream &oxidizer, double tol = 1.0,
                                    std::size_t max_iter = 100, bool lean = true,
                                    double phi_max = 10.0);

Stream set_oxidizer_stream_for_O2(double X_O2_target, const Stream &fuel,
                                   const Stream &oxidizer, double tol = 1e-6,
                                   std::size_t max_iter = 100);

Stream set_oxidizer_stream_for_O2_dry(double X_O2_dry_target,
                                       const Stream &fuel,
                                       const Stream &oxidizer, double tol = 1e-6,
                                       std::size_t max_iter = 100);

Stream set_oxidizer_stream_for_CO2(double X_CO2_target, const Stream &fuel,
                                    const Stream &oxidizer, double tol = 1e-6,
                                    std::size_t max_iter = 100);

Stream set_oxidizer_stream_for_CO2_dry(double X_CO2_dry_target,
                                        const Stream &fuel,
                                        const Stream &oxidizer,
                                        double tol = 1e-6,
                                        std::size_t max_iter = 100);

// State-based combustion functions
State complete_combustion(const State &in, bool smooth_phi0 = false, bool smooth_phi1 = false,
                         double k0 = SMOOTHING_K_PHI0, double k1 = SMOOTHING_K_PHI1);
State complete_combustion_isothermal(const State &in, bool smooth_phi0 = false, bool smooth_phi1 = false,
                                    double k0 = SMOOTHING_K_PHI0, double k1 = SMOOTHING_K_PHI1);

// Compute combustion state
CombustionState combustion_state(
    const std::vector<double> &X_fuel, const std::vector<double> &X_ox,
    double phi, double T_reactants, double P, const std::string &fuel_name = "",
    CombustionMethod method = CombustionMethod::Complete,
    bool smooth_phi0 = false, bool smooth_phi1 = false,
    double k0 = SMOOTHING_K_PHI0, double k1 = SMOOTHING_K_PHI1);

CombustionState combustion_state(const std::vector<double> &X_fuel,
                                 const std::vector<double> &X_ox, double phi,
                                 double T_reactants, double P,
                                 const std::string &fuel_name,
                                 CombustionMethod method,
                                 const PressureLossCorrelation &pressure_loss,
                                 bool smooth_phi0 = false, bool smooth_phi1 = false,
                                 double k0 = SMOOTHING_K_PHI0, double k1 = SMOOTHING_K_PHI1);

CombustionState combustion_state_from_streams(
    const Stream &fuel_stream, const Stream &ox_stream,
    const std::string &fuel_name = "",
    CombustionMethod method = CombustionMethod::Complete,
    bool smooth_phi0 = false, bool smooth_phi1 = false,
    double k0 = SMOOTHING_K_PHI0, double k1 = SMOOTHING_K_PHI1);

CombustionState combustion_state_from_streams(
    const Stream &fuel_stream, const Stream &ox_stream,
    const std::string &fuel_name, CombustionMethod method,
    const PressureLossCorrelation &pressure_loss,
    bool smooth_phi0 = false, bool smooth_phi1 = false,
    double k0 = SMOOTHING_K_PHI0, double k1 = SMOOTHING_K_PHI1);

} // namespace combaero
