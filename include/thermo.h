#ifndef THERMO_H
#define THERMO_H

#include "state.h"
#include "composition.h"
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace combaero::thermo {

// Universal gas constant [J/(mol*K)]
constexpr double R_GAS = 8.31446261815324;

// Boltzmann constant [J/K]
constexpr double BOLTZMANN = 1.380649e-23;

// Avogadro's number [1/mol]
constexpr double AVOGADRO = 6.02214076e23;

// Species thermo beyond its NASA-9 data (#606). Past each end of a species'
// data (T_lo = 200 K; T_hi = 6000 or 20000 K) cp/R leaves the end interval's
// polynomial C2 (value, slope and curvature matched) and levels off (zero
// slope and curvature) at the end of a fade; h and s are its exact integrals,
// so cp stays positive and h, s monotone at every T for the solvers and the
// inverses.
//   Below the data the gas is ideal down to 0 K: over the whole [0, T_lo] the
//   quintic fades cp/R to the rigid-rotor value, translation plus fully
//   excited rotation (5/2 for an atom, 7/2 linear, 4 nonlinear), as the
//   vibrations freeze out -- translation alone (5/2) for H2, whose rotation
//   freezes inside the fade (Theta_rot = 85 K; H2O's largest, the next, is
//   40 K). h(0 K) is finite, so a stream's enthalpy has a floor: less is more
//   heat out than it holds, at no T. Not modelled: condensation (a real
//   stream liquefies well above 0 K).
//   Above the data the quartic levels off after a fade of TAIL_FRACTION_HIGH
//   T_hi at whatever constant it reaches; not physics (no dissociation).
constexpr double TAIL_FRACTION_HIGH = 0.1;
// cp/R of an ideal gas with translation and rotation fully excited and the
// vibrations frozen: the low tail's value at 0 K (CP_R_RIGID_ATOM, translation
// alone, also for H2).
constexpr double CP_R_RIGID_ATOM = 2.5;
constexpr double CP_R_RIGID_LINEAR = 3.5;
constexpr double CP_R_RIGID_NONLINEAR = 4.0;
// h_RT (and g_over_RT, the equilibrium's g/RT) holds its value at this T
// below it (1/T); the kernels' stagnation-T inputs are floored here for an
// unphysical solver probe. Not a floor on any solved static T: those go to
// their root, however small (#935).
constexpr double T_FLOOR = 10.0;
// s_R holds its value below this T [K] (ln T), so that a non-positive probe
// T keeps s finite. It lies far below any isentropic root a double can carry
// (T/Tt = (P/Pt)^(R/cp) >= (1e-308)^(1/2.5) ~ 1e-123), so s is exact at every
// T a solve can reach.
constexpr double S_R_T_MIN = 1e-300;
// Below this T [K] the isentropic static-T Newton (solve_T_isentropic) steps
// in ln T, where s is ~ cp ln T and one step nearly lands on the root; above
// it the steps (and so the results) are those of the T Newton (#935). Equal
// to the former floor, so no iterate that stayed above it changes.
constexpr double T_LOG_STEP = 10.0;
// Pressure floor of the duct-state kernels [Pa]: an unphysical solver
// iterate (P <= 0) is evaluated there, finite, instead of returning NaN.
constexpr double P_FLOOR = 1.0;
// calc_T_from_cp searches [CP_INVERSE_T_MIN, CP_INVERSE_T_MAX] K. It scans
// outward from T_guess, each step on the side whose next node is nearer it,
// and polishes the first sign change of cp - target: the root nearest T_guess
// in K, to within the local step. Two roots closer than a step (a cp extremum
// within it of the target) can be missed (#868 R6); calc_T_from_cp_roots,
// which returns every root, finds them (#963).
//   Above the data (6000 K for H2O, NH3, CH4 and the alkanes, the shortest)
//   the high tail levels cp off at a constant that is not physics, so the
//   range ends there. Below the data (200 K) the low tail is the ideal gas
//   down to 0 K (#606), so cp identifies T there too; but the fade has zero
//   slope and curvature at 0 K, so cp - cp(0 K) ~ T^3 and the root's
//   condition number grows as 1/T^3. At CP_INVERSE_T_MIN = 1 K the least
//   sensitive species, CO (cp/R within 6e-10 of 7/2 at 1 K), still fixes T to
//   ~5e-7 relative from cp rounded to the last bit; at 0.1 K only to ~5e-4,
//   and below 1 K every species' cp is within 3e-6 relative of its 0 K value
//   (a mixture can cancel further). The tail's own extrema (near 100-130 K
//   for N2, O2, H2O, CH4, CO, NH3; 6.2 K for CO2) are scanned like the data's.
//   Down to CP_INVERSE_LOG_SCAN_T the nodes are CP_INVERSE_SCAN_STEP K apart;
//   below it, the same step relative to T (a factor 1 - SCAN_STEP /
//   LOG_SCAN_T = 0.9 per node, continuous at the switch), so that each step
//   spans a like share of the T^3 approach to 0 K (44 nodes from 100 K to
//   1 K). The downward scan steps onto CP_INVERSE_T_DATA, the data's lower end
//   where the tail joins, so for a guess in [200, 6000] K and a target with no
//   root below 200 K the nodes, the bracket and the result are those of the
//   former [200, 6000] K scan.
constexpr double CP_INVERSE_T_MIN = 1.0;
constexpr double CP_INVERSE_T_MAX = 6000.0;
constexpr double CP_INVERSE_T_DATA = 200.0;
constexpr double CP_INVERSE_SCAN_STEP = 10.0;
constexpr double CP_INVERSE_LOG_SCAN_T = 100.0;
// calc_T_from_cp_roots reports a cp extremum as a double root when cp there
// is within this share of cp_target: rounding level for cp, and far wider
// than tol in cp, so the two roots a slightly lower (higher) target would
// give are separated by much more than tol when they are reported apart (#963).
constexpr double CP_ROOTS_DOUBLE_REL = 1.0e-13;
// Golden-section steps of calc_T_from_cp_roots' search for an extremum pair
// hidden inside one scan step (a dip of |dcp/dT| between nodes): the window
// of two steps shrinks by 0.618 per step, to ~1e-11 K for 10 K steps.
constexpr int CP_ROOTS_DIP_ITER = 60;

} // namespace combaero::thermo

namespace combaero {

// Conversion factor from J/mol to J/kg
double J_per_mol_to_J_per_kg(double value, double molar_mass);

// Species metadata and lookup functions
std::string species_name(std::size_t idx);
std::size_t species_index_from_name(const std::string& name);
std::size_t num_species();

// Molar mass of a species [g/mol]
double species_molar_mass(std::size_t idx);
double species_molar_mass_from_name(const std::string& name);

// Which species thermo every function below evaluates (#549).
//   NASA9    : the NASA-9 polynomials as published. C0 in cp at
//              the interval joins (1000 K; 6000 K for N2, O2, AR, CO2, CO,
//              H2) but not C1: dcp/dT jumps 1-4 % at 1000 K (H2 20 %) and up
//              to 93 % at 6000 K, so every Jacobian using it steps there.
//   NASA9_C2 : the same polynomials with each join bridged, inside
//              |T - T_join| < 0.05 T_join, by a degree-7 cp/R that is C2 with
//              them at both window edges and whose h and s rises across the
//              window equal NASA-9's. Identical to NASA9 outside the windows
//              (h and s included); inside, cp within 1e-4 at 1000 K (5e-4 at
//              6000 K), h within 1e-5 cp T, s within 3e-5 R. Bridges are
//              generated (include/thermo_c2_refit_data.h). The default.
// Process-wide; thread-safe to read. Switch before solving, not during.
enum class ThermoModel { NASA9, NASA9_C2 };
void set_thermo_model(ThermoModel model);
ThermoModel thermo_model();
// Every T [K] where some species' cp is not smooth under the active model (a
// jump in cp or in one of its first three derivatives): each species' data
// ends (200 K; 6000 or 20000 K), its high tail's fade end, and the NASA-9
// joins (NASA9) or the C2 bridges' window edges, T_join -+ half_width
// (NASA9_C2). Sorted, unique. A Gauss-Legendre rule over T splits there, or
// misses the integral by more than its slope can follow (#1255).
const std::vector<double>& cp_breakpoints();

// NASA polynomial evaluations (dimensionless), with the tails beyond the data
// (thermo::TAIL_FRACTION_LOW/HIGH) under either model.
double cp_R(std::size_t species_idx, double T);
// d(cp/R)/dT of one species [1/K]; 0 on the constant tails.
double dcp_R_dT(std::size_t species_idx, double T);
double h_RT(std::size_t species_idx, double T);
double s_R(std::size_t species_idx, double T);
double g_over_RT(std::size_t species_idx, double T);

// Per-species dimensional properties
double cp_species(std::size_t species_idx, double T);  // [J/(mol*K)]
double h_species(std::size_t species_idx, double T);   // [J/mol]
double s_species(std::size_t species_idx, double T);   // [J/(mol*K)] at P_ref

// Thermodynamic properties
double cp(double T, const std::vector<double>& X);
double h(double T, const std::vector<double>& X);
double s(double T, const std::vector<double>& X, double P, double P_ref = 101325.0);
double cv(double T, const std::vector<double>& X);
double u(double T, const std::vector<double>& X);  // Internal energy [J/mol]
double density(double T, double P, const std::vector<double>& X);
double molar_volume(double T, double P);  // Ideal gas molar volume [m^3/mol]
double specific_gas_constant(const std::vector<double>& X);
double isentropic_expansion_coefficient(double T, const std::vector<double>& X);
double speed_of_sound(double T, const std::vector<double>& X);

// Mass-specific thermodynamic properties
double cp_mass(double T, const std::vector<double>& X);
// Secant cp over [T1, T2], (h(T1) - h(T2)) / (T1 - T2) [J/(kg K)], and its
// partials in T1 and T2; cp at T1 in the limit T2 -> T1. A stream's capacity
// rate over a span it may cross: the variable-property C of an epsilon-NTU
// exchange, which makes C_min (T_hot - T_cold) its exact maximum heat (#627).
std::tuple<double, double, double> cp_secant_mass_and_jacobian(double T1, double T2,
                                                               const std::vector<double>& X);
// The secant cp's gradient in the MASS fractions [J/(kg K)] per species
// (#717): a wall side's capacity rate |m| cp moves with its composition.
std::vector<double> cp_secant_mass_dY(double T1, double T2, const std::vector<double>& X);
// Gradients in the MASS fractions, per species, of R = R_GAS / MW [J/(kg K)]
// and of gamma = cp / (cp - R_GAS) (isentropic_expansion_coefficient), with
// gamma's dT [1/K]: an element taking a node's R and gamma as inputs chains
// them onto the node's composition (#1231, the legacy tee's kernel).
std::vector<double> specific_gas_constant_dY(const std::vector<double>& X);
std::pair<double, std::vector<double>> isentropic_expansion_coefficient_dT_dY(
    double T, const std::vector<double>& X);
double cv_mass(double T, const std::vector<double>& X);
double h_mass(double T, const std::vector<double>& X);
double s_mass(double T, const std::vector<double>& X, double P, double P_ref = 101325.0);
double u_mass(double T, const std::vector<double>& X);

// Compute all air properties at once
AirProperties air_properties(double T, double P, double humidity = 0.0);

// Compute all thermodynamic properties at once
ThermoState thermo_state(double T, double P, const std::vector<double>& X, double P_ref = 101325.0);

// Compute all thermodynamic and transport properties at once
CompleteState complete_state(double T, double P, const std::vector<double>& X, double P_ref = 101325.0);

// Derivatives of thermodynamic properties with respect to temperature, on the
// molar basis (#1279): dh/dT = cp [J/(mol*K)], ds/dT = cp/T [J/(mol*K^2)],
// dcp/dT [J/(mol*K^2)].
double dh_dT(double T, const std::vector<double>& X);
double ds_dT(double T, const std::vector<double>& X);
double dcp_dT(double T, const std::vector<double>& X);
double dg_over_RT_dT(double T, const std::vector<double>& X);

// Inverse calculations (solving for temperature). h, s, u and the sv form are
// strictly increasing in T at every T (the tails keep cp > 0), so each has one
// root: safeguarded Newton (factor-4 steps until bracketed, then bisection
// when a step leaves the bracket). They throw std::runtime_error if they do
// not converge in max_iter, rather than return the last iterate.
double calc_T_from_h(double h_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_s(double s_target, double P, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
// cp(T) is not monotone: the root nearest T_guess in [CP_INVERSE_T_MIN,
// CP_INVERSE_T_MAX] = [1, 6000] K (see there); std::runtime_error when none lies in it,
// std::invalid_argument for a non-finite or non-positive input (#868 R6).
double calc_T_from_cp(double cp_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
// One root of cp(T) = cp_target (calc_T_from_cp_roots).
struct CpRoot {
  double T = 0.0;       // [K]
  double dT_dcp = 0.0;  // 1 / (dcp/dT) at T [K/(J/(mol K))]; NaN for a double root
  bool double_root = false;  // T is a cp extremum and cp there is the target
};
// Every root of cp(T) = cp_target in [T_min, T_max], sorted (#963); the
// opt-in alternative to calc_T_from_cp's nearest-to-guess branch. Empty when
// there is none. Method: dcp/dT is evaluated on calc_T_from_cp's scan nodes
// (from T_min up, clamped at T_max); each sign change is bisected to the
// last bit, a cp extremum. An extremum pair inside one step leaves no sign
// change, but |dcp/dT| dips there: at each node where it is a local minimum
// with one sign around it, a golden-section search of the two neighbouring
// steps looks for the opposite sign (CP_ROOTS_DIP_ITER). Between neighbouring
// nodes and extrema cp is then monotone (barring three or more extrema within
// two steps), so each sign change of cp - target there is one root, polished
// as calc_T_from_cp polishes its bracket. An extremum whose cp is within
// CP_ROOTS_DOUBLE_REL of the target
// is one double root (dT_dcp NaN). std::invalid_argument for a non-finite or
// non-positive cp_target or tol, or a range outside [CP_INVERSE_T_MIN,
// CP_INVERSE_T_MAX] or with T_min >= T_max; std::domain_error when cp equals
// the target on an interval (argon's constant cp: the roots are not isolated).
std::vector<CpRoot> calc_T_from_cp_roots(double cp_target, const std::vector<double>& X,
                                         double T_min = thermo::CP_INVERSE_T_MIN,
                                         double T_max = thermo::CP_INVERSE_T_MAX,
                                         double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_u(double u_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);

// Inverse calculations for mass-specific targets
double calc_T_from_h_mass(double h_mass_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_s_mass(double s_mass_target, double P, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_u_mass(double u_mass_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_sv_mass(double s_mass_target, double v_mass_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_sh_mass(double s_mass_target, double h_mass_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);

// State-based overloads
double mwmix(const State& s);
double cp(const State& s);
double h(const State& s);
double s(const State& s);
double cv(const State& s);
double u(const State& s);
double density(const State& s);
double specific_gas_constant(const State& s);
double isentropic_expansion_coefficient(const State& s);
double speed_of_sound(const State& s);

} // namespace combaero

#endif // THERMO_H
