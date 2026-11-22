#pragma once

// -----------------------------------------------------------------------------
// 1-D critical-mode supersonic ejector model.
//
// C++ port of python/combaero/network/_ejector_huang1999.py -- see that
// module's docstring for the full paper-by-paper research narrative and
// validation/ejector/README.md for the accuracy comparison against
// alternative closures. Kept in exact correspondence with the Python
// reference so the two can be cross-validated row-for-row against
// validation/ejector/data/huang1999_reference_data.h (see
// tests/test_ejector_physics.cpp, tests/test_ejector_jacobians.cpp).
//
// References
// ----------
// [Huang 1999] B.J. Huang, J.M. Chang, C.P. Wang, V.A. Petrenko,
//     "A 1-D analysis of ejector performance",
//     International Journal of Refrigeration 22 (1999) 354-364.
//     Local copy: sources/ejector/Huang_1d_analysis_ejector.pdf
//     Source of ejector_entrainment_ratio (Eqs. 1-8) and
//     ejector_choked_mass_flow (Eqs. 1, 7).
// [Kracik & Dvorak 2016] J. Kracik, V. Dvorak,
//     "Development of an Analytical Method for Predicting Flow in a
//     Supersonic Air Ejector", EPJ Web of Conferences 114, 02059 (2016),
//     DOI: 10.1051/epjconf/201611402059.
//     Local copy: sources/ejector/Development_of_an_Analytical_Method_for.pdf
//     Source of the mixing closure ejector_critical_back_pressure
//     implements (their Eqs. 7-13).
//
// Scope. The critical-mode scalar closures (ejector_entrainment_ratio,
// ejector_critical_back_pressure and their Jacobians) model double-choked
// operation only; the additional operating-regime closures below
// (ejector_cd_nozzle_mass_flow, ejector_jetpump_discharge, ...) extend to the
// subcritical-droop and unchoked-primary jet-pump regimes, and
// ejector_element_residuals_and_jacobian composes ALL of them into the
// network element's full 4-row residual system spanning every regime. Design
// + provenance: validation/ejector/OPERATING_REGIMES_DESIGN.md.
//
// Jacobian scope. The critical-mode scalar closures expose derivatives w.r.t.
// the 4 thermodynamic inputs (p_g, t_g, p_e, t_e); geometry and
// recovery_efficiency are element parameters. The operating-regime closures
// additionally carry the P_py / omega / m_dot partials the extended element
// needs, and the whole-element function returns the full Jacobian w.r.t. all
// nine element unknowns AND each row's partials in gamma and r_gas: those are
// the secondary port gas's, which moves with the flows (#1246). The closures
// are written once against DualN (src/ejector.cpp), gamma and r_gas duals too.
//
// The one internal root-find (ejector_mach_from_area_ratio_supersonic, a
// bisection for the primary nozzle-exit Mach number M_p1) depends ONLY on
// area_ratio_nozzle and gamma, so M_p1 has an exact zero derivative w.r.t.
// p_g, t_g, p_e, t_e, and its gamma-partial is the implicit function
// theorem's on A/A*(M_p1, gamma) = area_ratio_nozzle (zero when clamped).
//
// Input validation (#769). Every function below throws std::invalid_argument
// (Python ValueError) for non-finite input and for parameters outside their
// physical domain: gamma > 1, r_gas > 0, nozzle efficiencies in (0, 1],
// recovery in (0, 1] (above 1 is a pressure gain beyond lossless mixing),
// 1 < area_ratio_nozzle <= area_ratio_mix, positive areas, area_exit >=
// area_throat, loss coefficients and eps_frac >= 0, smootherstep windows
// lo < hi. The standalone closures also need finite, positive supply
// stagnation states. The two whole-element residual functions take Newton
// iterates and check those for finiteness only: a trial point can carry a
// negative pressure or a reversed flow, which the closures continue through
// C1 floors (#481).
// -----------------------------------------------------------------------------

#include <array>
#include <tuple>

namespace ejector {
constexpr double eta_p = 0.95; // primary-nozzle isentropic efficiency (Huang 1999, Section 4)
constexpr double eta_s = 0.85; // entrained-flow isentropic efficiency
constexpr double phi_p = 0.88; // primary-core area-loss coefficient at section y-y (Eq. 5)
// Mixing-loss coefficient phi_m on the momentum balance (Huang 1999 Eq. 19),
// stepped in A_3/A_t: phi_m_large above the upper area ratio, phi_m_small
// below the lower one, phi_m_mid between (inclusive).
constexpr double phi_m_large = 0.80;
constexpr double phi_m_mid = 0.82;
constexpr double phi_m_small = 0.84;
constexpr double phi_m_upper_area_ratio = 8.3;
constexpr double phi_m_lower_area_ratio = 6.9;
// C1 floor on the sudden-expansion dead head's discriminant: at or below it the
// mixed flow chokes (M_3 = 1). Same value as _DEAD_HEAD_DISCRIMINANT_FLOOR in
// _ejector_huang1999.py.
constexpr double dead_head_discriminant_floor = 1e-12;
// Bracket of the supersonic area-Mach bisection (nozzle-exit Mach M_p1) and
// the bisections' iteration count (bracket / 2^200: exact to rounding).
constexpr double kSupersonicMachLo = 1.0 + 1e-12;
constexpr double kSupersonicMachHi = 50.0;
constexpr int kBisectionIterations = 200;
// The C-D nozzle exit-static inverse brackets P_py in [P*, Pt (1 - gap)].
constexpr double kExitStaticTopGap = 1e-12;
} // namespace ejector

namespace combaero::solver {

struct EjectorGeometry {
  double area_ratio_nozzle; // A_p1 / A_t
  double area_ratio_mix;    // A_3 / A_t
};

// -----------------------------------------------------------------------------
// Value-only functions (Python 1:1 parity; see tests/test_ejector_physics.cpp)
// -----------------------------------------------------------------------------

// Isentropic area ratio A/A* for a given Mach number.
double ejector_area_over_astar(double mach, double gamma);

// Supersonic root of A/A* = area_ratio by bisection (inverse of Eq. 2).
// A/A* is monotonically increasing for M > 1, so bisection is
// unconditionally convergent.
double ejector_mach_from_area_ratio_supersonic(double area_ratio, double gamma);

// The same root with an out-of-range flag (#945, flag and clamp). The
// bisection's bracket is [ejector::kSupersonicMachLo, kSupersonicMachHi]; an
// area ratio outside its image A/A*(bracket) has no root there and returns the
// bracket end: out_of_range = -1 below A/A*(kSupersonicMachLo), +1 above
// A/A*(kSupersonicMachHi) (e.g. ~7.8e3 at gamma = 5/3), else 0. M_p1 is a
// constant in every Newton unknown, so the clamp's derivative is zero, as the
// element's chain already uses.
struct EjectorNozzleExitMach {
  double mach;      // M_p1, or the bracket end it is clamped to
  int out_of_range; // -1 below, +1 above the bracket's area-ratio range, 0 in it
};
EjectorNozzleExitMach ejector_nozzle_exit_mach(double area_ratio, double gamma);

// Choked (sonic-throat) mass flow rate (Huang 1999 Eqs. 1 and 7).
//
// Method: Analytical Chain Rule (see ejector_choked_mass_flow_and_jacobian).
double ejector_choked_mass_flow(double Pt, double Tt, double area_throat,
                                double gamma, double r_gas, double eta);

// Converging-diverging nozzle mass flow spanning choked and unchoked, the R0
// primary residual of the operating-regime extension (see
// validation/ejector/OPERATING_REGIMES_DESIGN.md sec 3.1). Subsonic exit flux
// (area_exit, exit static p_static_down) smooth-min'd with the sonic-throat cap
// (area_throat); the two are equal at the choke threshold so the min is
// continuous, and the C1 sqrt-smoothing (eps = eps_frac * cap) rounds the
// corner. 1:1 with cd_nozzle_mass_flow in _ejector_huang1999.py.
double ejector_cd_nozzle_mass_flow(double Pt, double Tt, double p_static_down,
                                   double area_throat, double area_exit,
                                   double gamma, double r_gas, double eta,
                                   double eps_frac);

struct EjectorEntrainmentResult {
  double omega;                     // entrainment ratio ms / mp
  double mach_nozzle_exit;          // M_p1
  double mach_hypothetical_throat;  // M_py
  double area_ratio_primary_core;   // A_py / A_t
  double area_ratio_entrained;      // A_sy / A_t
  double p_mixing_pa;               // P_py = P_sy
};

// Critical-mode entrainment ratio omega (Huang 1999 Eqs. 1-8).
//
// Method: Analytical Chain Rule, composed around one precomputed bisection
// root (M_p1) whose derivative w.r.t. p_g/t_g/p_e/t_e is exactly zero --
// see the header comment above.
EjectorEntrainmentResult ejector_entrainment_ratio(double p_g, double t_g,
                                                    double p_e, double t_e,
                                                    const EjectorGeometry& geom,
                                                    double gamma);

// Aerodynamic mass-flow function q(lambda) (Kracik & Dvorak Eq. 9).
double ejector_q_lambda(double lam, double gamma);

struct EjectorCriticalPressureResult {
  double p_c_pa;                  // critical back pressure P_c*
  double p_mixed_stagnation_pa;   // p03, lossless (recovery_efficiency=1) mixed stagnation pressure
  double temp_mixed_stagnation_k; // T03
  double lambda_mixed;            // lambda3, subsonic root
  double mach_mixed;              // Mach number equivalent of lambda3
};

// Critical back pressure P_c* via Kracik & Dvorak's mixing closure
// (their Eqs. 7-13), continuing from the entrainment-ratio state above.
//
// Method: Analytical Chain Rule (same M_p1 treatment as
// ejector_entrainment_ratio).
EjectorCriticalPressureResult ejector_critical_back_pressure(
    double p_g, double t_g, double p_e, double t_e,
    const EjectorGeometry& geom, double gamma, double r_gas,
    double recovery_efficiency);

// -----------------------------------------------------------------------------
// Jacobian-carrying variants: (value, d/dp_g, d/dt_g, d/dp_e, d/dt_e).
// -----------------------------------------------------------------------------

// (m_dot, d_mdot_dPt, d_mdot_dTt), Pt and Tt the primary stagnation state
std::tuple<double, double, double>
ejector_choked_mass_flow_and_jacobian(double Pt, double Tt, double area_throat,
                                      double gamma, double r_gas, double eta);

struct EjectorCDNozzleJacobian {
  double mdot;        // == ejector_cd_nozzle_mass_flow
  double dmdot_dp0;   // d/d(primary stagnation pressure)
  double dmdot_dt0;   // d/d(primary stagnation temperature)
  double dmdot_dp_py; // d/d(exit/mixing static pressure P_py)
};

// C-D nozzle mass flow with its analytic Jacobian w.r.t. the three inputs that
// are Newton-solved unknowns for EjectorElement (Pt, Tt = the primary
// stagnation state, p_static_down = P_py). gamma/r_gas/geometry/eta are held
// constant here (the whole element takes the gamma and r_gas partials too,
// #1246). The dmdot_dp0/dmdot_dt0 fields are d/dPt and d/dTt.
// The secondary entrained flux reuses this with area_throat = area_exit = A_s.
EjectorCDNozzleJacobian ejector_cd_nozzle_mass_flow_and_jacobian(
    double Pt, double Tt, double p_static_down, double area_throat,
    double area_exit, double gamma, double r_gas, double eta, double eps_frac);

struct EjectorCDExitStaticJacobian {
  double p_py;         // exit static that passes m_dot (== cd_nozzle_exit_static)
  double dp_py_dm_dot; // d/d(mass flow)
  double dp_py_dp0;    // d/d(primary stagnation pressure)
  double dp_py_dt0;    // d/d(primary stagnation temperature)
  // Flag and clamp (#945): +1 when m_dot is above the flow at P_py = P* (the
  // choked end; P_py is held at P*), -1 when below the flow at Pt (1 - gap)
  // (P_py held there), 0 in range. A clamped result's partials are those of
  // the clamped function k Pt: (0, k, 0).
  int out_of_range;
};

// Inverse of the C-D nozzle: the exit static P_py that passes m_dot, with its
// implicit-function-theorem Jacobian. The derived mixing pressure of the
// operating-regime element's unchoked branch (OPERATING_REGIMES_DESIGN.md sec
// 6b); 1:1 with cd_nozzle_exit_static in _ejector_huang1999.py. The partials
// are ratios of the C-D nozzle's own Jacobian at the solved P_py:
//   dP_py/dm_dot = 1/(dmdot/dP_py),  dP_py/dPt = -(dmdot/dPt)/(dmdot/dP_py), ...
// An m_dot outside the nozzle's range is clamped and flagged (out_of_range).
EjectorCDExitStaticJacobian ejector_cd_nozzle_exit_static_and_jacobian(
    double m_dot, double Pt, double Tt, double area_throat, double area_exit,
    double gamma, double r_gas, double eta, double eps_frac);

struct EjectorJetPumpDischargeJacobian {
  double p03;      // recovery_efficiency * mixed stagnation pressure
  double dp03_dp_g;
  double dp03_dt_g;
  double dp03_dp_e;
  double dp03_dt_e;
  double dp03_dp_py; // d/d(mixing static pressure P_py)
  double dp03_domega; // d/d(entrainment ratio -- the R1 unknown feeding the mix)
};

// Jet-pump mixed-flow discharge stagnation: the Kracik & Dvorak mixing closure
// (reused from the critical path, ejector_critical_back_pressure) evaluated at
// the SUBSONIC-primary jet-pump state -- both streams expand to the common
// mixing static P_py (lambda1/lambda2 from that expansion), mixed at ratio
// omega. The R3 (outlet-pin) building block for the operating-regime element;
// 1:1 with _mixed_flow_stagnation in _ejector_huang1999.py fed jet-pump
// lambdas. r_gas cancels (lambda/q are dimensionless), so it is not an input.
//
// eta_primary / eta_secondary (#588): the nozzles' kinetic-energy
// efficiencies, the same ones the element's flux rows carry (sqrt(eta) on
// the isentropic flux). Each stream reaches the mixing plane at the velocity
// its lossy nozzle delivers, lambda = sqrt(eta) lambda_is, and at the
// stagnation pressure consistent with that velocity at P_py,
// p0 = P_py / pi(lambda) = p_supply pi(lambda_is) / pi(lambda) -- below the
// supply by exactly the nozzle loss. Mixing isentropic velocities instead
// handed the nozzle losses back: at zero drive the jet pump held a flow with
// zero entropy generation. 1 (the default) is the isentropic form, unchanged.
EjectorJetPumpDischargeJacobian ejector_jetpump_discharge_and_jacobian(
    double p_g, double t_g, double p_e, double t_e, double p_py, double omega,
    double gamma, double recovery_efficiency, double eta_primary = 1.0,
    double eta_secondary = 1.0);

struct EjectorEntrainmentJacobian {
  EjectorEntrainmentResult value;
  double domega_dp_g;
  double domega_dt_g;
  double domega_dp_e;
  double domega_dt_e;
};

EjectorEntrainmentJacobian ejector_entrainment_ratio_and_jacobian(
    double p_g, double t_g, double p_e, double t_e,
    const EjectorGeometry& geom, double gamma);

struct EjectorCriticalPressureJacobian {
  EjectorCriticalPressureResult value;
  double dpc_dp_g;
  double dpc_dt_g;
  double dpc_dp_e;
  double dpc_dt_e;
};

EjectorCriticalPressureJacobian ejector_critical_back_pressure_and_jacobian(
    double p_g, double t_g, double p_e, double t_e,
    const EjectorGeometry& geom, double gamma, double r_gas,
    double recovery_efficiency);

// Dead-head back pressure P_b0 (Phase 1A, OPERATING_REGIMES_DESIGN.md sec 8.1):
// the highest back pressure a choked ejector holds with NO entrainment, the
// upper anchor of the subcritical droop. The SAME mixing closure as P_c*,
// mixed at omega = 0 (one shared chain); 1:1 with dead_head_back_pressure in
// _ejector_huang1999.py. The result struct's p_c_pa field carries P_b0
// (recovered), its p_mixed_stagnation_pa the lossless p03 at omega = 0.
EjectorCriticalPressureJacobian ejector_dead_head_back_pressure_and_jacobian(
    double p_g, double t_g, double p_e, double t_e,
    const EjectorGeometry& geom, double gamma, double r_gas,
    double recovery_efficiency);

// The operating-anchor closures (P_c*, P_b0) the whole element can take; 1:1
// with ANCHOR_CLOSURES in _ejector_huang1999.py. The pair always changes
// together, so the droop window P_b0/P_c* is one closure family's.
enum class EjectorAnchorClosure : int {
  kracik_dvorak = 0, // ejector_critical_back_pressure + ejector_dead_head_back_pressure
  huang1999 = 1,     // ejector_huang_critical_back_pressure + ejector_sudden_expansion_dead_head
};

// phi_m of Huang 1999 Eq. 19 (Chen 2013 Eq. 28) for a mixing area ratio.
double ejector_huang_mixing_coefficient(double area_ratio_mix);

// Critical back pressure P_c* by Huang 1999's own chain, as printed: from the
// critical y-y state, constant-pressure mixing with phi_m (Eq. 11), energy
// (Eq. 12), a normal shock in the mixed flow (Eqs. 16-17), an isentropic
// diffuser (Eq. 18). 1:1 with huang_critical_back_pressure in
// _ejector_huang1999.py, where its fidelity gap and cross-source accuracy are
// recorded. C2 through M_m = 1 (a normal shock's loss is O((M - 1)^3)).
// Result fields: p_c_pa = recovery * P_c*, p_mixed_stagnation_pa = P_c*,
// temp_mixed_stagnation_k = the mixed stagnation T, lambda/mach_mixed = the
// post-shock state.
//
// Method: Analytical Chain Rule (DualN, same M_p1 treatment as
// ejector_entrainment_ratio).
EjectorCriticalPressureJacobian ejector_huang_critical_back_pressure_and_jacobian(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom, double gamma,
    double r_gas, double recovery_efficiency);

// Zero-entrainment back pressure P_b0 as a sudden expansion over A_3 (Keenan &
// Neumann's zero-secondary-flow constant-area mixer): the choked primary
// leaves the nozzle into the mixing tube, the dead-headed suction annulus at
// P_e, mass/momentum/energy over A_3 to a subsonic exit, recovered
// isentropically. The exit Mach is the closed-form subsonic root of a
// quadratic in M^2 (floored C1 at the choked exit). 1:1 with
// sudden_expansion_dead_head in _ejector_huang1999.py. t_e does not enter
// (zero partial). Result fields as ejector_dead_head_back_pressure_and_jacobian.
//
// Method: Analytical Chain Rule (DualN).
EjectorCriticalPressureJacobian ejector_sudden_expansion_dead_head_and_jacobian(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom, double gamma,
    double r_gas, double recovery_efficiency);

// -----------------------------------------------------------------------------
// Whole-element (f, J): the 4-row operating-regime residual system.
//
// Composes the scalar closures above (choked/C-D nozzle mass flow, critical
// entrainment + back pressure, jet-pump discharge) into the four coupled rows
// the network EjectorElement solves, blending the critical/subcritical/unchoked
// regimes by two smootherstep weights, and returns all four residuals with
// their analytic Jacobian w.r.t. the nine Newton unknowns -- forward-mode dual
// (DualN<9>) through the whole blend. This is the C++ home of the assembly
// (matching the whole-element (f, J) practice of MultiPortChamberBase /
// TeeJunctionElement); the Python EjectorElement is a thin relabeling shim.
// 1:1 with EjectorElement.residuals() in ejector_element.py. Full derivation
// and provenance: validation/ejector/OPERATING_REGIMES_DESIGN.md sec 6c/8.3.
//
// Unknown (seed) order -- both `jac` columns and the physical-flow convention:
//   0 mp        primary mass flow (physical, = -port_mdot_primary)
//   1 ms        secondary (entrained) mass flow (physical)
//   2 mdot_out  outlet mass flow (physical)
//   3 p_g       primary stagnation pressure   Pt_primary
//   4 t_g       primary stagnation temperature Tt_primary
//   5 p_e       secondary stagnation pressure  Pt_secondary
//   6 t_e       secondary stagnation temperature Tt_secondary
//   7 p_out     outlet-node stagnation pressure Pt_outlet
//   8 p_py      the owned unknown, repurposed as the mixing-plane static P_py
// and, beside `jacobian`, each row's partials in the gas: d_dgamma and
// d_dr_gas. gamma and r_gas are the secondary port's, and its gas moves with
// the flows (a reversed leg delivers the ejector's own mix there, a side
// stream changes it), so the element chains these onto that node's T0 and
// composition (#1246).
struct EjectorElementResidualJacobian {
  std::array<double, 4> residuals;               // R0..R3
  std::array<std::array<double, 9>, 4> jacobian; // jacobian[row][seed]
  std::array<double, 4> d_dgamma{};              // dR_row / dgamma
  std::array<double, 4> d_dr_gas{};              // dR_row / dr_gas [per J/(kg K)]
  // ejector_nozzle_exit_mach's out_of_range at (area_ratio_nozzle, gamma):
  // nonzero when the nozzle-exit Mach the closures use is a clamped bracket
  // end (#945).
  int nozzle_exit_mach_out_of_range = 0;
};

EjectorElementResidualJacobian ejector_element_residuals_and_jacobian(
    double mp, double ms, double mdot_out, double p_g, double t_g, double p_e,
    double t_e, double p_out, double p_py, const EjectorGeometry& geom,
    double area_throat, double area_nozzle_exit, double area_secondary,
    double gamma, double r_gas, double eta_primary, double eta_secondary,
    double recovery_efficiency, double eps_frac, double s_choke_lo,
    double s_choke_hi, double s_sub_lo, double s_sub_hi,
    EjectorAnchorClosure anchor_closure = EjectorAnchorClosure::huang1999);

// ---------------------------------------------------------------------------
// Incompressible jet pump (#588)
// ---------------------------------------------------------------------------
//
// The low-Mach limit of the ejector, for a network solved incompressibly:
// Sanger's conventional one-dimensional jet-pump theory (NASA TN D-4445,
// 1968, Appendix B, Eqs. 1-13), with each stream at its own constant density
// (specific volume v = 1/rho). Same owned unknown as the compressible element
// (p3 = the mixing-plane static pressure, at the nozzle exit and the suction
// annulus) and the same four rows:
//
//   R0 = Pt_p - p3 - (1 + K_nozzle)  v_p mp|mp| / (2 A_e^2)     nozzle
//   R1 = Pt_s - p3 - (1 + K_suction) v_s ms|ms| / (2 A_s^2)     suction
//   R2 = mo - mp - ms                                           mass
//   R3 = Pt_o - recovery (p4 + (1 - K_diffuser) q4)             outlet
//
// with the constant-area mixing duct A_3 (momentum, mixing complete):
//   p4 = p3 + [v_p mp|mp| / A_e + v_s ms|ms| / A_s - |mt| (mp v_p + ms v_s) / A_3] / A_3
//           - K_throat q4,     q4 = |mt| (mp v_p + ms v_s) / (2 A_3^2)
// mt = mp + ms; the mixed stream's specific volume is volume-additive,
// v_m = (mp v_p + ms v_s) / mt, so mt|mt| v_m = |mt| (mp v_p + ms v_s): no
// division at zero flow. At equal densities this is exactly Sanger's Eq. 11
// (P5 = Pt_o, P2 = Pt_s); the K's are his component losses (default 0) and
// `recovery` multiplies the outlet stagnation pressure as the compressible
// element's does, so matched losses compare like for like. No choking, no
// critical regime: the flow follows the pressure ratio.
//
// Seeds (jacobian columns): 0 mp, 1 ms, 2 mo, 3 Pt_p, 4 v_p, 5 Pt_s, 6 v_s,
// 7 Pt_o, 8 p3. Rows R0, R1, R3 are pressures; R2 is a mass flow.
struct EjectorJetPumpIncompressibleResult {
  std::array<double, 4> residuals;
  std::array<std::array<double, 9>, 4> jacobian;
  double p4 = 0.0;  // mixing-duct exit static [Pa]
  double q4 = 0.0;  // mixing-duct exit dynamic head [Pa]
};

EjectorJetPumpIncompressibleResult ejector_jetpump_incompressible_residuals_and_jacobian(
    double mp, double ms, double mo, double pt_p, double v_p, double pt_s,
    double v_s, double pt_o, double p3, double area_nozzle_exit,
    double area_secondary, double area_mixing, double k_nozzle,
    double k_suction, double k_throat, double k_diffuser, double recovery);

} // namespace combaero::solver
