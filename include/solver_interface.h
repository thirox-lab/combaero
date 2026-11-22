#ifndef SOLVER_INTERFACE_H
#define SOLVER_INTERFACE_H

#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "state.h"
#include "combustion.h"
// -----------------------------------------------------------------------------
// Network Solver Fast-Path Native Interface
// -----------------------------------------------------------------------------
// Functions in this module execute highly optimized `(f, J)` combined
// evaluations required by the global Python `scipy.optimize.root` solver.
// Note: Concepts currently named "Pipe" are scheduled for renaming to "Channel"
// to account for advanced non-circular geometries and aerospace cooling applications.
//
// Every function returns a `std::tuple<double, double>`, corresponding to:
//    (Value, Exact Analytical Jacobian / Derivative)
//
// This architecture natively evaluates chain-rule analytical gradients where
// possible, and tightly bounded finite-differences where strictly necessary,
// entirely shielding the Python runtime from grid discontinuity noise and
// locking the C++/Python return overhead to ~0.2 microseconds.

namespace combaero {
namespace solver {

enum class CorrelationValidity : std::uint8_t { VALID, EXTRAPOLATED, INVALID };

template <typename T> struct CorrelationResult {
  T result;
  CorrelationValidity status;
  std::string message;
};

// -----------------------------------------------------------------------------
// Stream Definitions for Network Solvers
// -----------------------------------------------------------------------------

struct Stream {
  double m_dot;
  double Tt;  // total (stagnation) temperature [K]
  double Pt;  // total pressure [Pa]
  std::vector<double> Y;
};

struct StreamJacobian {
  double d_mdot;
  double d_Tt;
  double d_Pt;
  std::vector<double> d_Y;
};

// The mixed TOTAL state: Tt_mix from the h(Tt) balance, Pt_mix the
// mass-averaged total pressure (#782 names).
struct MixerResult {
  double Tt_mix;
  double Pt_mix;
  std::vector<double> Y_mix;
  double dTt_mix_d_delta_h = 0.0;
  // dTt_mix/dQ [K/W] = dTt_mix_d_delta_h * phi, with phi = delta_h / Q the
  // heat per unit flow (mixer_heat_per_unit_flow: 1/|m| above
  // kMixerHeatMdotFloor, fading to 0 at zero flow below it). A
  // caller relaying a heat input (a wall) needs this, not 1/(cp sum m),
  // which was wrong below the floor and for a zero or negative total (#481).
  double dTt_mix_dQ = 0.0;
  std::vector<StreamJacobian> dTt_mix_d_stream;
  std::vector<StreamJacobian> dPt_mix_d_stream;
  std::vector<std::vector<StreamJacobian>> dY_mix_d_stream;
};

struct ChamberResult {
  std::vector<double> residuals;
  // Local Jacobian: map from eq_idx to {var_name: derivative}
  // var_name is a string like "P", "P_total", "T", "Y[0]"...
  std::vector<std::map<std::string, double>> local_jacobian;
  // Stream Jacobian: map from eq_idx to list of StreamJacobian (one per
  // upstream stream)
  std::vector<std::vector<StreamJacobian>> stream_jacobian;
};

struct OrificeResult {
  double m_dot_calc;
  double d_mdot_dP_total_up;
  double d_mdot_dP_static_down;
  double d_mdot_dP_static_up;
  double d_mdot_dP_rho_ref = 0.0;  // density reference named by P_rho_ref
  double d_mdot_dTt_up;
  std::vector<double> d_mdot_dY_up;
};

struct ChannelResult {
  double dP_calc;
  double d_dP_d_mdot;
  double d_dP_dP_static_up;
  double d_dP_dTt_up;
  std::vector<double> d_dP_dY_up;
};

using PipeResult = ChannelResult;

// -----------------------------------------------------------------------------
// 1. Incompressible Flow Components
// -----------------------------------------------------------------------------

// Calculate Mass Flow and its derivative with respect to differential pressure.
//
// Method: Analytical Chain Rule
// Inputs:
//   dP   : Differential pressure drop across the orifice [Pa]
//   rho  : Fluid density [kg/m^3]
//   Cd   : Discharge coefficient [-]
//   beta : Orifice beta ratio (d/D) [-]
//
// Returns:
//   tuple(mdot, d_mdot_d_dP, d_mdot_d_rho)
//   mdot         : Mass flow rate [kg/s]
//   d_mdot_d_dP  : Jacobian gradient w.r.t. pressure difference [kg/(s*Pa)]
//   d_mdot_d_rho : Jacobian gradient w.r.t. density [kg*m^3/(s*kg)] = [m^3/s]
std::tuple<double, double, double> orifice_mdot_and_jacobian(double dP, double rho,
                                                              double Cd, double area,
                                                              double beta);

// Does the orifice drive from the upstream STATIC pressure? True when the
// velocity-of-approach factor E = 1/sqrt(1 - beta^4) applies (0 < beta < 1):
// E already carries the approach velocity head, so driving from the upstream
// stagnation pressure as well counted it twice -- +6.2% flow at beta = 0.7
// behind a momentum chamber. Without E (beta = 0, a plenum-fed orifice) the
// drive is the upstream stagnation pressure. Behind a plenum node, where
// Pt = P, both rules give the same flow.
inline bool approach_from_static(double beta) { return beta > 0.0 && beta < 1.0; }

// Full orifice evaluation with all derivatives. The drive pressure follows
// approach_from_static. The density is evaluated at P_static_up, or at
// P_rho_ref when that is > 0 (its derivative then goes to d_mdot_dP_rho_ref).
OrificeResult orifice_residuals_and_jacobian(double m_dot, double P_total_up,
                                             double P_static_up, double Tt_up,
                                             const std::vector<double> &Y_up,
                                             double P_static_down, double Cd,
                                             double area, double beta,
                                             double P_rho_ref = -1.0);

// Calculate Pressure Drop and its derivative with respect to velocity.
//
// Method: Analytical Chain Rule
// Inputs:
//   v   : Fluid velocity [m/s]
//   rho : Fluid density [kg/m^3]
//   K   : Forward loss coefficient [-]
//
// Returns:
//   tuple(dP, d_dP_d_v)
//   d_dP_d_v : Jacobian gradient [Pa/(m/s)]
std::tuple<double, double> pressure_loss_and_jacobian(double v, double rho,
                                                       double K);

// Full channel evaluation with all derivatives.
//
// dP_calc is the stagnation-pressure drop: Darcy friction
// f (L/D) m|m| / (2 rho A^2), plus, when exit_head_lost, the exit dynamic head
// m|m| / (2 rho A^2) -- a duct discharging into a stagnant sink (a pressure
// boundary taking outflow at its static pressure) dissipates its exit velocity
// head, exactly as the compressible Fanno form does (#588). rho is the
// upstream static density, the same constant density the friction term uses.
ChannelResult channel_residuals_and_jacobian(double m_dot, double P_total_up,
                                             double P_static_up, double Tt_up,
                                             const std::vector<double> &Y_up,
                                             double P_static_down, double L,
                                             double D, double roughness,
                                             const std::string &friction_model,
                                             double f_multiplier = 1.0,
                                             bool exit_head_lost = false);



// Calculate Lossless Connection Ideal Pressure and Jacobian
//
// Method: Analytical Difference
// Inputs:
//   P_in  : Total pressure entering the element [Pa]
//   P_out : Total pressure exiting the element [Pa]
//
// Returns:
//   tuple(Residual, d_Residual_d_P_out)
//   Residual           : Difference between P_in and P_out [Pa]
//   d_Residual_d_P_out : Jacobian gradient [-]
std::tuple<double, double> lossless_pressure_and_jacobian(double P_in,
                                                          double P_out);

// -----------------------------------------------------------------------------
// 1b. Compressible Flow Components
// -----------------------------------------------------------------------------

// Compressible orifice flow using isentropic nozzle model.
// Applies discharge coefficient and velocity-of-approach factor to nozzle_flow.
//
// Method: analytic Jacobians on the isentrope through the nozzle_flow outlet
// state, the exact derivatives of the returned mdot on both sides of choke.
// Below kJetDriveFloor (1e-3 Pa) of drive |Pt - P_back| the flow is the floor's
// times the odd cubic of flow_drive_floor.h, exponent 1/2 (#1256): a finite slope at zero
// drive, where the jet's is infinite, and no round-off staircase near it.
// Inputs:
//   Tt   : Stagnation temperature [K]
//   Pt   : Stagnation pressure [Pa]
//   P_back : Back pressure [Pa]
//   X    : Mole fractions [-]
//   Cd   : Discharge coefficient [-]
//   area : Geometric throat area [m^2]
//   beta : Orifice beta ratio (d/D) [-]
//
// Returns:
//   tuple(mdot, d_mdot_dPt, d_mdot_dP_back, d_mdot_dTt)
//   mdot           : Mass flow rate [kg/s]
//   d_mdot_dPt     : Jacobian w.r.t. stagnation pressure [kg/(s*Pa)]
//   d_mdot_dP_back : Jacobian w.r.t. back pressure [kg/(s*Pa)] (zero when choked; continuous through P*)
//   d_mdot_dTt     : Jacobian w.r.t. stagnation temperature [kg/(s*K)]
std::tuple<double, double, double, double> orifice_compressible_mdot_and_jacobian(
    double Tt, double Pt, double P_back,
    const std::vector<double>& X,
    double Cd, double area, double beta);

// Full compressible orifice evaluation with all derivatives for network solver.
// Tt_up is the upstream STAGNATION temperature. The nozzle expands from the
// upstream state approach_from_static picks: without E from (Tt_up,
// P_total_up); with E from the upstream static state, P_static_up at the
// temperature on the isentrope through (Tt_up, P_total_up) (#645).
OrificeResult orifice_compressible_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double Tt_up,
    const std::vector<double>& Y_up,
    double P_static_down, double Cd, double area, double beta);

// The opt-in compressible orifice models (#854), with OrificeResult's columns
// plus d(mdot)/d(Cd): Jobson's flow is not linear in Cd, so the network
// cannot take m/Cd for it.
struct OrificeModelResult : OrificeResult {
  double d_mdot_dCd = 0.0;
  double Ci = 0.0;             // incompressible contraction (Jobson); Cd (ISO)
  double contraction = 0.0;    // Jobson's compressible C; 0 for ISO
  double expansibility = 1.0;  // mdot / (Cd E A sqrt(2 rho1 dp)), rho1 static
  bool choked = false;
  bool in_range = true;        // inside the model's stated range (a flag only)
};

// Jobson (1955) with the approach-velocity extension (orifice.h). Inputs as
// orifice_compressible_residuals_and_jacobian; Cd is the element's
// incompressible coefficient in ISO's form, Cv Jobson's velocity coefficient
// (Ci from jobson_ci_from_cd). With E (0 < beta < 1) the approach static
// state is (P_static_up, isentropic from the total); without, a reservoir at
// (P_total_up, Tt_up). The jet expands from the TOTAL state, polytropic with
// n = gamma(Tt_up). Back pressure above the drive: the mirrored flow from
// a reservoir at P_static_down (the raw kernel's reverse band; the network
// evaluates reversed flow from the other node). Analytic (Pt, P1, Tt, Pb,
// Cd) columns by DualN; composition analytic too (#964), through the two
// numbers it enters by, n = gamma(Tt_up) and R, carried as DualN directions.
// in_range: jobson1955::ci_min <= Ci <= ci_max, beta <=
// beta_max and a real root.
OrificeModelResult orifice_jobson_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double Tt_up,
    const std::vector<double>& Y_up, double P_static_down, double Cd, double area,
    double beta, double Cv = 1.0);

// ISO 5167-2:2003 Eq. (1) with Eq. (5)'s expansibility:
//   mdot = Cd E eps A sqrt(2 rho1 (p1 - p2)),
// rho1 the approach static density (isentropic from the total state, kappa
// = gamma(Tt_up)). Eq. (5) is extrapolated below p2/p1 = 0.75 (in_range
// false there; no choke exists in it). The drive is regularised at 1 Pa
// like orifice_residuals_and_jacobian; reversed, eps uses p1/p2.
OrificeModelResult orifice_iso5167_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double Tt_up,
    const std::vector<double>& Y_up, double P_static_down, double Cd, double area,
    double beta);

// An orifice plate between pipes recovers pressure downstream of its vena
// contracta (#645). The network row sees the RECOVERED drop, the stagnation
// difference between the two pipes' nodes, while C is referenced to the
// tapping differential dp = dw / k (k: permanent_loss_ratio_and_dC). This
// maps the recovered drop to the downstream TAP pressure the kernels expand
// to:
//   p_tap = P_static_up - (Pt_up - Pt_down) / k(beta, C).
// Compressible, p_tap is held at ORIFICE_TAP_PRESSURE_FLOOR P_static_up,
// below any critical ratio, so the choked nozzle reads it as choked; the
// incompressible map has no floor. ISO's relation is incompressible:
// p_tap / P_static_up < ORIFICE_RECOVERY_MIN_PRESSURE_RATIO is outside it
// (in_range false), ISO 5167-2's own expansibility limit.
// The static pressure an orifice's approach pipe carries (#651): the flow at
// mass flux G = |m| / A_pipe from the upstream node's stagnation state (Pt,
// Tt). With E applied the orifice drives from this, not from its node's own
// static pressure, so a PlenumNode (P = Pt by definition) and a
// MomentumChamberNode give the same flow, and beta, E and p1 all refer to the
// same pipe. Compressible: duct_static_state (entropy and energy, as the
// momentum chamber's closure). Incompressible: the chamber's low-Mach form
// Pt = P + G^2 R Tt / (2 P), whose subsonic root is
// P = [Pt + sqrt(Pt^2 - 2 G^2 R Tt)] / 2; past its maximum flux it holds
// P = Pt / 2 (choked = true, derivatives of that branch). dP_dY analytic
// (#964): the low-Mach root through R; subsonic, the implicit-function gradient
// of duct_static_state's two equations (duct_htc_static_state_dY); past
// choke, that of its sonic point (M = 1 energy balance and the isentrope).
// The sonic static pressure P* of a chamber at stagnation Tt and mass flux G
// (#637): the floor its static pressure cannot go below while it discharges
// into a reservoir (combaero::sonic_static_state). dP_dY by a forward
// difference (approach_static_pressure's is analytic since #964).
struct SonicStatic {
  double P;
  double T;
  double dP_dTt;
  double dP_dG;
  std::vector<double> dP_dY;
};

SonicStatic sonic_static_pressure(double Tt, double G, const std::vector<double>& Y);

struct ApproachStatic {
  double P;
  double dP_dPt;
  double dP_dTt;
  double dP_dG;
  std::vector<double> dP_dY;
  bool choked;
};

ApproachStatic approach_static_pressure(double Pt, double Tt, double G,
                                        const std::vector<double>& Y,
                                        bool incompressible);

inline constexpr double ORIFICE_TAP_PRESSURE_FLOOR = 0.1;
inline constexpr double ORIFICE_RECOVERY_MIN_PRESSURE_RATIO = 0.75;

struct OrificeTapPressure {
  double p_tap;
  double loss_ratio;        // k = dw / dp
  double d_dP_static_up;
  double d_dP_total_up;    // d/dPt_up
  double d_dP_total_down;  // d/dPt_down
  double d_dC;              // through k(C)
  bool floored;
  bool in_range;
};

OrificeTapPressure orifice_tap_pressure_and_jacobian(double P_static_up,
                                                     double Pt_up,
                                                     double Pt_down,
                                                     double beta, double C,
                                                     bool compressible);

// -----------------------------------------------------------------------------
// 2. Heat Transfer Components
// -----------------------------------------------------------------------------

// Calculate Nusselt Number and its derivative with respect to Reynolds Number.
// Uses the Dittus-Boelter fundamental relation: Nu = 0.023 * Re^0.8 * Pr^n
//
// Method: Analytical Chain Rule
// Inputs:
//   Re      : Reynolds number [-]
//   Pr      : Prandtl number [-]
//   heating : True if fluid is being heated (n=0.4), False if cooled (n=0.3)
//
// Returns:
//   tuple(Nu, d_Nu_d_Re)
//   Nu        : Nusselt number [-]
//   d_Nu_d_Re : Jacobian gradient [-]
CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_dittus_boelter(double Re, double Pr, bool heating = true);

// Calculate Nusselt Number and its derivative with respect to Reynolds Number
// using the Gnielinski correlation.
//
// Method: Finite Central Difference
CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_gnielinski(double Re, double Pr, double f);

// Calculate Nusselt Number and its derivative with respect to Reynolds Number
// using the Sieder-Tate correlation.
//
// Method: Analytical Chain Rule
CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_sieder_tate(double Re, double Pr, double mu_ratio);

// Calculate Nusselt Number and its derivative with respect to Reynolds Number
// using the Petukhov correlation.
//
// Method: Analytical Chain Rule
CorrelationResult<std::tuple<double, double>>
nusselt_and_jacobian_petukhov(double Re, double Pr, double f);

// Calculate Friction Factor and its derivative with respect to Reynolds Number.
// Uses the Haaland explicit pipe friction relation.
//
// Method: Analytical Chain Rule
// Inputs:
//   Re  : Reynolds number [-]
//   e_D : Relative roughness [-]
//
// Returns:
//   tuple(f, d_f_d_Re)
//   f        : Darcy friction factor [-]
//   d_f_d_Re : Jacobian gradient [-]
CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_haaland(double Re, double e_D);

// Calculate Friction Factor and its derivative using the Serghides explicit
// approximation to Colebrook-White (Steffensen acceleration).
//
// Method: Analytical Chain Rule through A, B, C intermediates
// Inputs:
//   Re  : Reynolds number [-]
//   e_D : Relative roughness [-]
//
// Returns:
//   tuple(f, d_f_d_Re)
CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_serghides(double Re, double e_D);

// Calculate Friction Factor and its derivative using the Colebrook-White
// implicit equation (Newton-Raphson converged internally).
//
// Method: Implicit differentiation of converged Colebrook root
// Inputs:
//   Re  : Reynolds number [-]
//   e_D : Relative roughness [-]
//
// Returns:
//   tuple(f, d_f_d_Re)
CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_colebrook(double Re, double e_D);

// Calculate Friction Factor and its derivative using the Petukhov smooth-pipe
// correlation f_turb = (coeff_a * ln(Re_eff) - coeff_b)^(-2), at the turbulent
// fits' regularized Re_eff = sqrt(Re^2 + 64^2), blended with Hagen-Poiseuille
// like the other network kernels (friction.h, laminar_friction; #1210). The
// raw fit has a pole at Re 7.97; friction_petukhov is the fit alone.
//
// Method: Analytical Chain Rule
// Inputs:
//   Re : Reynolds number [-]  (smooth pipe: e_D not used)
//
// Returns:
//   tuple(f, d_f_d_Re)
CorrelationResult<std::tuple<double, double>>
friction_and_jacobian_petukhov(double Re);

// Convenience dispatcher: selects the correct friction_and_jacobian_* overload
// based on a string tag matching FrictionModelLiteral from Python.
//
// tag must be one of: "haaland", "serghides", "colebrook", "petukhov"
//
// Returns:
//   tuple(f, d_f_d_Re)
CorrelationResult<std::tuple<double, double>>
friction_and_jacobian(const std::string &tag, double Re, double e_D);

// lim_{Re -> 0} f Re of the laminar-turbulent blend (laminar_friction in
// friction.h): f_re (1 - w(0)), w(0) = (1 + tanh(-blend_center /
// blend_width)) / 2. The turbulent fits are finite at Re = 0 (their own
// floor), so w f_turb Re -> 0. The limit drops the laminar branch's re_floor,
// which only keeps 64 / Re finite at Re = 0 (#1156).
double laminar_blend_f_re_at_zero();

// lim_{Re -> 0} f Re of friction model `tag` (Registry::friction_f_re_at_zero):
// the slope of a duct's drop at exactly zero flow is f_mult f_re_0 mu L /
// (2 rho A D^2) (#1156). Throws std::invalid_argument for an unknown tag.
double friction_f_re_at_zero(const std::string &tag);


// -----------------------------------------------------------------------------
// 4. Thermodynamic & Transport Components
// -----------------------------------------------------------------------------

// Calculate Mass Density and its analytical derivatives w.r.t. T and P.
// Method: Ideal Gas Law Analytical Chain Rule
// Inputs:
//   T : Temperature [K]
//   P : Pressure [Pa]
//   X : Mole fractions [-]
//
// Returns:
//   tuple(rho, d_rho_d_T, d_rho_d_P)
//   rho       : Density [kg/m^3]
//   d_rho_d_T : Jacobian gradient w.r.t Temperature [(kg/m^3)/K]
//   d_rho_d_P : Jacobian gradient w.r.t Pressure [(kg/m^3)/Pa]
std::tuple<double, double, double>
density_and_jacobians(double T, double P, const std::vector<double> &X);

// Calculate Mass Enthalpy and its analytical derivative w.r.t. T.
// Method: dh/dT is identically the mass-specific heat capacity (cp_mass).
// Inputs:
//   T : Temperature [K]
//   X : Mole fractions [-]
//
// Returns:
//   tuple(h, d_h_d_T)
//   h       : Mass Enthalpy [J/kg]
//   d_h_d_T : Jacobian gradient w.r.t Temperature [J/(kg*K)]
std::tuple<double, double> enthalpy_and_jacobian(double T,
                                                 const std::vector<double> &X);

// Calculate Dynamic Viscosity and its derivatives w.r.t. T and P.
// Method: Central Finite Difference (tightly bracketed)
// Inputs:
//   T : Temperature [K]
//   P : Pressure [Pa]
//   X : Mole fractions [-]
//
// Returns:
//   tuple(mu, d_mu_d_T, d_mu_d_P)
//   mu       : Dynamic viscosity [Pa*s]
//   d_mu_d_T : Jacobian gradient w.r.t Temperature [(Pa*s)/K]
//   d_mu_d_P : Jacobian gradient w.r.t Pressure [(Pa*s)/Pa]
std::tuple<double, double, double>
viscosity_and_jacobians(double T, double P, const std::vector<double> &X);

// d(mu)/dY_k [Pa*s per unit mass fraction] at fixed (T, P), analytic (Wilke's
// rule, combaero::viscosity_dX chained to the mass fractions), at the same
// T >= 10 K leaky floor as viscosity_and_jacobians (#964). X are the mole
// fractions of the state whose mass fractions are differentiated. For the
// Re-dependent orifice Cd: Re ~ 1/mu, so dRe/dY = -Re dmu/dY / mu.
std::vector<double> viscosity_mass_gradient(double T, double P, const std::vector<double> &X);

// -----------------------------------------------------------------------------
// 5. Stagnation State Conversions
// -----------------------------------------------------------------------------

// Calculate Mach number and its derivative wrt Velocity
// Method: Analytical Chain Rule
std::tuple<double, double>
mach_number_and_jacobian_v(double v, double T, const std::vector<double> &X);

// Calculate Adiabatic Wall Temperature and its derivative wrt Velocity:
// T_aw = T_static + r v^2 / (2 cp). T and P: the state the gas properties
// (Pr for the recovery factor r, cp) are read at, normally the same static
// state as T_static.
// Method: Analytical Chain Rule
std::tuple<double, double>
T_adiabatic_wall_and_jacobian_v(double T_static, double v, double T, double P,
                                const std::vector<double> &X,
                                bool turbulent = true);

// Calculate Stagnation Temperature from Static and its derivative wrt Mach
// Method: Central Finite Difference
std::tuple<double, double>
Tt_from_static_and_jacobian_M(double T, double M, const std::vector<double> &X);

// Calculate Stagnation Pressure from Static and its derivative wrt Mach
// Method: Central Finite Difference
std::tuple<double, double>
Pt_from_static_and_jacobian_M(double P, double T, double M,
                              const std::vector<double> &X,
                              double tol = 1e-8, std::size_t max_iter = 50);

// -----------------------------------------------------------------------------
// 6. Combustion Interfaces
// -----------------------------------------------------------------------------

// Calculate Adiabatic Flame Temperature (Complete Combustion) and its
// derivative wrt inlet temperature.
// Method: analytic (#716), dT_ad/dT_in = cp(T_in, X_in) (MW_b/MW_in) / cp(T_ad, X_b);
// smoothed at phi = 0 and phi = 1 (complete_combustion_products_and_jacobian).
// Inputs:
//   T_in : Inlet temperature [K]
//   P    : Pressure [Pa]
//   X_in : Inlet mole fractions [-]
// Returns:
//   tuple(T_ad, dT_ad_dT_in, X_products)
//   T_ad        : Adiabatic flame temperature [K]
//   dT_ad_dT_in : Jacobian gradient [-]
//   X_products  : Mole fractions of products [-]
std::tuple<double, double, std::vector<double>>
adiabatic_T_complete_and_jacobian_T(double T_in, double P,
                                    const std::vector<double> &X_in);

// Calculate Adiabatic Flame Temperature (Equilibrium Combustion) and its
// derivatives wrt inlet temperature and pressure.
// Method: analytic (#716), implicit in the equilibrium extents and the adiabatic root.
// Inputs:
//   T_in : Inlet temperature [K]
//   P    : Pressure [Pa]
//   X_in : Inlet mole fractions [-]
// Returns:
//   tuple(T_ad, dT_ad_dT_in, dT_ad_dP, X_products)
//   T_ad        : Adiabatic flame temperature [K]
//   dT_ad_dT_in : Jacobian gradient wrt Inlet Temperature [-]
//   dT_ad_dP    : Jacobian gradient wrt Pressure [K/Pa]
//   X_products  : Mole fractions of products [-]
std::tuple<double, double, double, std::vector<double>>
adiabatic_T_equilibrium_and_jacobians(double T_in, double P,
                                      const std::vector<double> &X_in);

// -----------------------------------------------------------------------------
// 7. Stream-Based Network Solvers
// -----------------------------------------------------------------------------


// Below this total flow [kg/s] a mixing node takes up only part of the heat Q
// it is given (mixer_heat_per_unit_flow). At or above it Q/m is exact, so the
// heat a node is given is the heat it receives (energy conservation; the
// previous sqrt(m^2 + 1e-6) lost 5e-5 of Q at 0.1 kg/s and 11% at 2 g/s).
// 1 mg/s: below any flow a real network carries, so only a stagnant node (a
// dead end) withholds any; the solver reports it as {node}.Q_withheld.
// Convergence measured the same at 2e-3 and 1e-6 (#481).
constexpr double kMixerHeatMdotFloor = 1e-6;

// Near-zero band of a mixer's total flow [kg/s]: within |m| <= this, the
// mixed state is the same arriving-stream mean as outside it, but its
// partials are 0 (#1232). The mean's own partials grow as 1 / m, and at a
// near-zero start they sent the first Newton step the wrong way. A single
// stream has no such pole and keeps its partials (#1265).
constexpr double kMixerNoFlowBand = 1e-12;

// Relative mismatch combustor_residuals_and_jacobians accepts between its P
// and the streams' mixed Pt (round-off of the same mean computed twice).
constexpr double kCombustorPressureMatchRtol = 1e-12;

// Lowest pressure [Pa] the combustion kernels evaluate a burn at: a solver
// probe below it is burnt at this pressure, with zero pressure partials.
constexpr double kEquilibriumPressureFloor = 100.0;

// The enthalpy rise per unit heat, phi = delta_h / Q [1/(kg/s)], a mixing node
// applies at total flow mdot_tot, and d(phi)/d(mdot_tot). 1/|m| at or above
// kMixerHeatMdotFloor m0; below it phi = (4 u^2 - 3 u^3) / m0, u = |m| / m0,
// so the heat taken up, Q |m| phi = Q (4 u^3 - 3 u^4), fades C1 from all of Q
// at m0 to none at zero flow (phi <= 1.053 / m0, smooth through m = 0). A
// stagnant node therefore holds the temperature of what reaches it, not
// Q / m0 above it: the C1 floor m_eff = (m^2 + m0^2) / (2 m0) this replaced
// put 100 W into a dead end at 1e5 K at the solution, which a [200, 5000] K
// node clamp hid (#626).
std::tuple<double, double> mixer_heat_per_unit_flow(double mdot_tot);

// Reference temperature [K] of the sensible enthalpy h(T) - h(T_ref) that an
// EnergyBoundary 'fraction' scales (at the stream's own composition).
constexpr double kSensibleEnthalpyRefT = 298.15;

// Mix generalized non-reacting incoming streams into a single outgoing state,
// and compute all analytical Jacobians w.r.t upstream mass flows, Temperatures,
// and mass fractions (Y_i).
// Streams are signed (m_i > 0 arrives at the node). Y, h and Pt are the
// mass-weighted mean sum_i m_i q_i / m over the arriving streams, or over
// every stream when none arrives (#1224), and the partials are that mean's; a
// leaving stream beside an arriving one moves nothing (#1232). Over streams
// of one sign the mean stays within its inputs at any nonzero m; with no flow
// at all stream 0's state stands in. Within |m| <= kMixerNoFlowBand the
// partials are 0, except a single stream's (its own state: weight 1).
// Q [W]: absolute heat transfer rate (positive = heating)
// fraction [-]: relative to total enthalpy flow (positive = heating)
MixerResult
mixer_from_streams_and_jacobians(const std::vector<Stream> &streams,
                                double Q = 0.0, double fraction = 0.0);

// Internal entropy generation of a control volume (#545): the second-law
// check the network applies to every element after a solve.
//
//   sigma = sum_out m s(Tt, Pt, Y) - sum_in m s(Tt, Pt, Y) - S_heat
//
// Streams carry STAGNATION states (Tt, Pt, as in the mixer); s(Tt, Pt) =
// s(T, P) by definition of the stagnation state. S_heat
// [W/K] is the entropy the heat brought in, integral dQ / T at the fluid's
// static temperature (heat_entropy_flux), which makes sigma the INTERNAL
// generation: a necessary condition that needs no wall temperature.
//
// sigma_norm = sigma / (R m_in), with R the specific gas constant of the
// largest inflow. For an adiabatic two-port of fixed composition it is
// exactly -ln(Pt_out / Pt_in).
struct EntropyGenerationResult {
  double sigma;       // [W/K]
  double sigma_norm;  // [-]
};

EntropyGenerationResult entropy_generation(const std::vector<Stream> &inflows,
                                           const std::vector<Stream> &outflows,
                                           double S_heat = 0.0);

// Entropy [W/K] that heat added to a stream brings in: integral of dQ / T as
// the stream's total temperature runs from Tt_in to Tt_out, dQ = m_dot
// dh(Tt), at the STATIC temperature T(Tt, M) (T_from_stagnation). The Mach
// number runs linearly in Tt from M_in to M_out; pass M_in == M_out for a
// constant one. At M = 0 it is m_dot [s(Tt_out, P) - s(Tt_in, P)], for any
// P. Heating at Mach M brings in more entropy than at M = 0 -- the static
// temperature is lower -- which is the total-pressure cost of Rayleigh flow,
// and cooling the gain. Variable cp throughout (composite Simpson, n = 64).
double heat_entropy_flux(double m_dot, double Tt_in, double Tt_out,
                         const std::vector<double> &Y, double M_in = 0.0,
                         double M_out = 0.0);

// Combust multiple generalized incoming streams to exactly complete
// equilibrium, outputting the theoretical adiabatic flame temperature and
// product mass fractions. Computes all exact analytical Jacobians w.r.t
// upstream conditions.
// Q [W]: absolute heat transfer rate (positive = heating)
// fraction [-]: of the products' sensible enthalpy at T_ad, h(T_ad) -
// h(kSensibleEnthalpyRefT) (#481) (positive = heating)
MixerResult adiabatic_T_complete_and_jacobian_T_from_streams(
    const std::vector<Stream> &streams, double P, double Q = 0.0, double fraction = 0.0);

// Combust multiple generalized incoming streams to generalized equilibrium,
// outputting the theoretical adiabatic flame temperature and product mass
// fractions. Computes all exact analytical Jacobians w.r.t upstream conditions.
// Q [W]: absolute heat transfer rate (positive = heating)
// fraction [-]: of the products' sensible enthalpy at T_ad (#481)
MixerResult adiabatic_T_equilibrium_and_jacobians_from_streams(
    const std::vector<Stream> &streams, double P, double Q = 0.0, double fraction = 0.0);
// Full combustor evaluation including pressure loss analytical chain rule.
// P is the streams' mixed total pressure (the mixer's Pt_mix), checked to
// kCombustorPressureMatchRtol (std::invalid_argument otherwise); the stream
// partials chain the equilibrium's dT/dP and dY/dP through it (#1244).
MixerResult combustor_residuals_and_jacobians(
    const std::vector<Stream> &streams, double P, double Q, double fraction,
    const PressureLossCorrelation &pressure_loss, bool use_equilibrium = false);

// Momentum chamber residual: Pt = Pt(P, T), the isentropic stagnation
// pressure of the chamber's static state. Tt is the chamber's STAGNATION
// temperature (#576); the static T is solved inside from
// h(Tt) = h(T) + v^2 / 2, v = m_dot / (rho(T, P) * A).
// Returns residual and Jacobian entries for P, Pt (d_res_dP_total), m_dot,
// Tt (the *_dT fields) and Y.
struct MomentumChamberResult {
  double residual;
  double d_res_dP;
  double d_res_dP_total;
  double d_res_dmdot;
  // Sensitivity to the stagnation temperature Tt: through the static T it
  // implies and through s(Tt). Closed form (#576).
  double d_res_dT;
  // Composition sensitivity, per species mass fraction (Y is normalised
  // inside, so any mass-conserving direction of the relay is exact). Y
  // enters through rho, a and s. A species absent at the node keeps its
  // partial: a stream arriving at exactly zero flow brings it (#1245).
  std::vector<double> d_res_dY;
  // The static state the closure solved for (#576): static T [K] and Mach
  // number [-] from (Tt, P, |m_dot| / A). Reporting reads these so the Mach
  // it shows is the one the residual used.
  double T_static = 0.0;
  double M = 0.0;
  // Partials of T_static w.r.t. the closure's inputs P, Tt and m_dot (signed:
  // T_static is even in m_dot), for a wall correlation evaluated at the
  // chamber's static T. The incompressible limit has T_static = Tt.
  double dT_static_dP = 0.0;
  double dT_static_dT = 0.0; // d/dTt
  double dT_static_dmdot = 0.0;
  // ... and w.r.t. the mass fractions (#743), per species; zero in the
  // incompressible limit.
  std::vector<double> dT_static_dY;
  // Static density rho(T_static, P) and its sensitivities to the closure's
  // inputs, for elements that read the chamber as a port (junction kernels):
  // their momentum rows must see the same static state as this row.
  // Isentropic dynamic head Pt - P, computed as P expm1(ln(Pt/P)) rather
  // than as a difference of two pressures (one ULP of Pt would otherwise
  // be its noise floor). Its partials are the residual's, negated, with
  // d/dP shifted by one: dq/dP = -d_res_dP - 1.
  double q_dyn = 0.0;
  double rho_static = 0.0;
  double drho_dP = 0.0;
  double drho_dT = 0.0;     // d/dTt
  double drho_dmdot = 0.0;
  std::vector<double> drho_dY;
  // Branch selection (#788): true when P is below the sonic static pressure
  // P* of the chamber's flux, where the closure is continued by its mirror
  // image about Pt*, 2 Pt* - Pt_cl(P), so that each Pt has one root. Such a
  // root has Pt < Pt*: the chamber cannot carry its flux at that stagnation
  // pressure (a choked duct discharging into it). The static state fields
  // above stay those of the closure at P. P_sonic / Pt_sonic are P* and Pt*
  // [Pa], set only when bounded. Always false in the incompressible limit.
  bool bounded = false;
  double P_sonic = 0.0;
  double Pt_sonic = 0.0;
};

// incompressible (#588): the low-Mach limit of the same closure, for a
// network solved incompressibly. The static state IS the stagnation one
// (T_static = Tt, rho = P / (R Tt)) and Pt = P + (|m_dot| / A)^2 / (2 rho):
// no sonic limit, so no fold for a neighbour to drive the chamber onto. It
// agrees with the compressible closure to O(Ma^2).
MomentumChamberResult momentum_chamber_residual_and_jacobian(
    double P, double Pt, double m_dot, double Tt,
    const std::vector<double> &Y, double area, bool incompressible = false);

// 4. Area change elements: see area_change.h (area_change_flow), called by
// the network element directly (#731).

// -----------------------------------------------------------------------------
// 5. Tee Junction Elements (Bassett 2001)
// -----------------------------------------------------------------------------
// Three-port pressure-loss elements: MergingTee and BranchingTee.
// Both residuals are normalised by 0.5*rho*u_com^2.
//
// Port convention:
//   MergingTee:   MAIN_INLET=A, BRANCH=B (inflow),  MAIN_OUTLET=C (common)
//   BranchingTee: MAIN_OUTLET=A, BRANCH=B (outflow), MAIN_INLET=C (common)
//
// Inputs common to both wrappers:
//   m_dot_com    : total flow at common port C [kg/s] (positive = design direction)
//   m_dot_branch : branch flow at port B [kg/s] (positive = design direction)
//   dP0_straight : p0_upstream - p0_downstream for the straight path [Pa]
//   dP0_branch   : p0_upstream - p0_downstream for the branch path [Pa]
//   P_static_com : static pressure at common port [Pa] (used for density)
//   T_com        : STATIC temperature at common port [K] (pairs with P_static_com)
//   Y_com        : mass fractions at common port [-]
//   theta        : branch angle [rad] (validated range: (0, pi/2])
//   psi          : area ratio F_C/F_B [-] (validated range: >= 0.05)
//   F_C          : main duct cross-section area [m^2]
//   blend_k      : tanh blend sharpness (default 30.0)

struct TeeJunctionResult {
    // Straight path residual and m_dot Jacobians (analytical)
    double R_straight;
    double dR_straight_d_mdot_com;
    double dR_straight_d_mdot_branch;
    // Straight path: density-dependent Jacobians (finite difference)
    double dR_straight_dP_static_com;
    double dR_straight_dT_com;
    std::vector<double> dR_straight_dY_com;

    // Branch path residual and m_dot Jacobians (analytical)
    double R_branch;
    double dR_branch_d_mdot_com;
    double dR_branch_d_mdot_branch;
    // Branch path: density-dependent Jacobians (finite difference)
    double dR_branch_dP_static_com;
    double dR_branch_dT_com;
    std::vector<double> dR_branch_dY_com;

    // Diagnostics
    double K_straight;      // effective straight-path loss coefficient
    double K_branch;        // effective branch-path loss coefficient
    double q;               // flow ratio m_dot_branch / m_dot_com
    double blend_w;         // tanh blend weight (1 = primary topology)
    bool topology_valid;    // true if q in design direction (within 5% margin)
    CorrelationValidity status; // VALID if inputs in validated Bassett range
};

TeeJunctionResult merging_tee_residuals_and_jacobian(
    double m_dot_com, double m_dot_branch,
    double dP0_straight, double dP0_branch,
    double P_static_com, double T_com,
    const std::vector<double>& Y_com,
    double theta, double psi, double F_C,
    double blend_k = 30.0);

TeeJunctionResult branching_tee_residuals_and_jacobian(
    double m_dot_com, double m_dot_branch,
    double dP0_straight, double dP0_branch,
    double P_static_com, double T_com,
    const std::vector<double>& Y_com,
    double theta, double psi, double F_C,
    double blend_k = 30.0);

// -----------------------------------------------------------------------------
// 6. Compressible Unified0D Tee Junction (Mynard & Valen-Sendstad 2015)
// -----------------------------------------------------------------------------
// Replaces the Bassett incompressible K tables with a geometry-derived K to
// O(M^2) that is consistent with mass-flow continuity and stagnation pressure.
// See sources/junction/junction_model.tex for derivation.
//
// Per-branch input struct.  mdot sign: >0 = supplier (inflow to junction).
// P_static and Pt are both passed explicitly (available as state_X.P / Pt in
// the Python solver; C++ does not need to invert the stagnation relation).
struct BranchInput {
    double P_static; // static pressure [Pa]
    double Pt;       // stagnation pressure [Pa]
    double T;        // static temperature [K]
    double m_dot;    // signed mass flow [kg/s]  (>0 = supplier)
    double A;        // cross-section area [m^2]
    double theta;    // centreline angle [rad]
    double gamma_eff;
    double R_gas;    // specific gas constant [J/(kg*K)]
};

// Result: two residuals [Pa] and their Jacobians w.r.t. all coupled unknowns.
// Naming: dR{0|1}_d{Pt|P|T|mdot}_{com|str|bra}
// Pt = stagnation pressure, P = static pressure, T = static temperature.
// Not all fields are nonzero for both tee types -- zero fields are explicitly set.
struct CompressibleTeeResult {
    double R_0;              // first residual [Pa]
    double R_1;              // second residual [Pa]
    // R_0 Jacobians
    double dR0_dPt_com;      // d(R_0)/d(Pt_com)
    double dR0_dP_com;       // d(R_0)/d(P_static_com)
    double dR0_dT_com;       // d(R_0)/d(T_com)
    double dR0_dmdot_com;    // d(R_0)/d(m_dot_com)
    double dR0_dPt_str;      // d(R_0)/d(Pt_str)
    double dR0_dP_str;       // d(R_0)/d(P_static_str)
    double dR0_dT_str;
    double dR0_dPt_bra;
    double dR0_dP_bra;       // d(R_0)/d(P_static_bra)
    double dR0_dT_bra;
    double dR0_dmdot_branch; // d(R_0)/d(m_dot_branch element unknown)
    // R_1 Jacobians
    double dR1_dPt_com;
    double dR1_dP_com;
    double dR1_dT_com;
    double dR1_dmdot_com;
    double dR1_dPt_str;
    double dR1_dP_str;
    double dR1_dT_str;
    double dR1_dPt_bra;
    double dR1_dP_bra;
    double dR1_dT_bra;
    double dR1_dmdot_branch;
    // The ports' gas properties (#1231): d/d(R_gas) and d/d(gamma_eff) of
    // each port, which an element chains onto its nodes' composition. Both
    // kernels take R only in the product R T, so d/dR = (T / R) d/dT.
    double dR0_dRgas_com, dR0_dRgas_str, dR0_dRgas_bra;
    double dR1_dRgas_com, dR1_dRgas_str, dR1_dRgas_bra;
    double dR0_dgamma_com, dR0_dgamma_str, dR0_dgamma_bra;
    double dR1_dgamma_com, dR1_dgamma_str, dR1_dgamma_bra;
};

// Branching tee: common is the single supplier; straight and branch are collectors.
//   R_0 = p0_dat - p0_straight - K_straight * q_ref
//   R_1 = p0_dat - p0_branch   - K_branch   * q_ref
// com.m_dot and bra.m_dot carry the element unknowns m_dot_com, m_dot_branch.
// str.m_dot = com.m_dot - bra.m_dot (derived; pass the current iterate value).
CompressibleTeeResult compressible_branching_tee_rj(
    const BranchInput& com,
    const BranchInput& str,
    const BranchInput& bra);

// Merging tee: straight and branch are suppliers; common is the collector.
//   R_0 = P_static_straight - P_static_branch   (supplier static-pressure equality)
//   R_1 = p0_dat - p0_com - K_com * q_ref
// All three m_dot values are passed (str.m_dot + bra.m_dot == com.m_dot by mass balance).
CompressibleTeeResult compressible_merging_tee_rj(
    const BranchInput& com,
    const BranchInput& str,
    const BranchInput& bra);

} // namespace solver
} // namespace combaero

#endif // SOLVER_INTERFACE_H
