#ifndef AREA_CHANGE_H
#define AREA_CHANGE_H

#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Area change closures (sudden and conical expansion / contraction), #731.
// -----------------------------------------------------------------------------
//
// Source: I. E. Idel'chik, Handbook of Hydraulic Resistance (AEC-tr-6630,
// 1966), transcribed in validation/area_change/ (signed off 2026-10-25):
//
//   diagram 4-1  sudden expansion: Borda-Carnot (1 - n)^2 for Re >= 3.5e3,
//                zeta(Re, n) tabulated for 10 < Re < 3.5e3, 26/Re for Re < 8;
//   diagram 3-9  sharp contraction: 0.5 (1 - n), Re > 1e4;
//   diagram 3-10 sharp contraction, Re < 1e4 (Karev): zeta(Re, n) tabulated,
//                27/Re for Re < 8;
//   diagram 3-7  conical contraction: zeta' (1 - n), zeta'(l/D, full angle)
//                from diagram 3-6 (friction separate);
//   diagram 5-2  conical diffuser: zeta_exp(full angle, n) = phi (1 - n)^2.
//
// n is always the small area over the large one, in [0, 1]; every zeta is on
// the dynamic head of the SMALL section, and Re = w_small D_h / nu there.
//
// Compressibility. Idelchik's coefficients are incompressible. The element
// carries them to Mach by the control volume they come from:
//
//   expansion   the compressible Borda-Carnot step: the base pressure acts on
//               the whole large face, mass, momentum and total enthalpy are
//               conserved (compressible_borda_carnot). Its loss times the
//               incompressible table's ratio to Borda-Carnot (the low-Re
//               factor of 4-1, or phi of 5-2).
//   contraction a vena contracta of area Cc F_small re-expanding to F_small by
//               the same compressible Borda-Carnot step, with
//               Cc = 1 / (1 + sqrt(zeta_c)), zeta_c the high-Re zeta (3-9 or
//               3-7), so its incompressible limit is that zeta exactly; times
//               3-10's low-Re ratio for a sharp one. It chokes at the vena
//               contracta.
//
// Both reduce to the tables as M -> 0. The gas is calorically perfect across
// the step, at gamma(Tt): within 0.2% of the thermally perfect control volume
// up to M 0.65 at 1500 K (0.03% at 300 K), test_area_change.py.

namespace combaero {

// The low-Re tables' bounds on Re (small section). Below Re 8 both pages give
// a laminar zeta = A / Re; above the last column the high-Re zeta holds.
constexpr double LOW_RE_LAMINAR = 8.0;
constexpr double EXPANSION_LAMINAR_COEFF = 26.0;     // 4-1
constexpr double EXPANSION_RE_TURBULENT = 3500.0;    // 4-1: Borda-Carnot above
constexpr double CONTRACTION_LAMINAR_COEFF = 27.0;   // 3-10
constexpr double CONTRACTION_RE_TURBULENT = 1.0e4;   // 3-10: 0.5 (1 - n) above

// The loss exponent d ln(zeta Re^2) / d ln Re that a low-Re interpolation
// keeps between two table columns at least (see sudden_expansion_zeta). Every
// 3-10 secant is above it already (0.95 at least).
constexpr double LOW_RE_MIN_SECANT_EXPONENT = 0.25;

// Small-section dynamic head [Pa] that sets a branch's drive floor: the floor
// is the incompressible loss there, zeta_inc * AREA_CHANGE_HEAD_FLOOR (about
// 1 Pa for zeta ~ 0.5, as kFannoFlowDriveFloor). Below it the flow is blended
// C1 into an odd function through zero (area_change_flow); the smaller of the
// two directions' floors is the width of the forward/reverse blend
// (area_change_element_flow). A head, not a drive: at a fixed drive a step of
// n 0.99 (zeta 1e-4) would reach M 0.35, where the floor's low-Mach exponent
// no longer holds. 2 Pa is M ~0.005 at 1 bar.
constexpr double AREA_CHANGE_HEAD_FLOOR = 2.0;

struct ZetaAndSlope {
    double zeta = 0.0;
    double dzeta_dRe = 0.0;
};

// Sudden expansion, diagram 4-1, on the small section's dynamic head.
//
//   Re <= 8      26 / Re (laminar: the loss is linear in the flow);
//   Re >= 3500   Borda-Carnot (1 - n)^2;
//   in between   per printed column, a cubic Hermite in ln(zeta Re^2) against
//                ln Re through the printed values, with knot slopes the
//                harmonic mean of the neighbouring secants (Fritsch-Butland),
//                slope 1 at Re 8 (the laminar law) and 2 at Re 3500
//                (Borda-Carnot): zeta and its Re-slope continuous everywhere.
//
// zeta Re^2 is the loss at fixed geometry and fluid, so its slope is the
// solver's: a falling stretch would give one pressure drop two flows. The
// printed n = 0.6 column falls between Re 2000 and 3000 (secant exponent
// -0.26), so a knot is lowered until every secant exponent is at least
// LOW_RE_MIN_SECANT_EXPONENT, working back from the Borda-Carnot anchor:
// zeta(2000) at n 0.6 is 0.41 (printed 0.50) and at n 0.5 0.61 (printed
// 0.65). Every other printed value is matched; the Re 3500 column, printed
// to two decimals, is Borda-Carnot exactly. The smallest local exponent left
// is 0.058 (n >= 0.6 near Re 3100), where the steepest printed drop has to
// turn onto Borda-Carnot's slope 2.
//
// In n, the ratio zeta / (1 - n)^2 is linear between the printed columns
// (0.1 ... 0.6) and the edge column's is held outside them: Borda-Carnot is
// exact at Re >= 3500 for every n, zeta -> 0 as n -> 1, and between columns
// at Re < 8 zeta is within 1% of 26/Re (the page's n-independent law).
ZetaAndSlope sudden_expansion_zeta(double Re, double area_ratio);

// Sudden (sharp) contraction, on the small section's dynamic head:
//
//   Re <= 8      27 / Re (3-10, laminar);
//   Re >= 1e4    0.5 (1 - n) (3-9, Weisbach);
//   in between   3-10 (Karev), interpolated as sudden_expansion_zeta. Every
//                printed value is matched: no secant needs lowering.
//
// The ratio to 0.5 (1 - n) is linear in n between the printed columns and the
// edge column's is held outside them, as for the expansion.
ZetaAndSlope sudden_contraction_zeta(double Re, double area_ratio);

// Conical diffuser, diagram 5-2: zeta_exp(alpha, n), alpha the FULL included
// angle [deg], bilinear in (alpha, n) through the printed table (exact at
// every printed cell). phi = zeta_exp / (1 - n)^2 is held at the n = 0.6 row
// above it, and falls linearly to 0 at alpha = 0 below the 3 deg column.
// Friction is not included (the table's zeta_fr is separate).
double conical_diffuser_zeta(double full_angle_deg, double area_ratio);

// Conical contraction (nozzle), diagram 3-7: zeta' (1 - n), zeta'(l/D, alpha)
// from diagram 3-6 (conical converging bellmouth with end wall), alpha the
// FULL included angle [deg] and l / D the cone's axial length over the small
// diameter. Bilinear in (alpha, l/D) through the printed table (exact at every
// printed cell); l/D is held at the table's ends (0.025, 0.6). Both angle ends
// are the sharp edge, zeta' = 0.5 (3-9). Re > 1e4; friction not included.
double conical_nozzle_zeta(double full_angle_deg, double length_over_D, double area_ratio);

// Compressible Borda-Carnot step from area F1 to F2 = F1 / n, calorically
// perfect gas: with X = 1 + (gamma - 1) M^2 / 2,
//   mass + momentum   M2 sqrt(X2) / (1 + gamma M2^2)
//                       = n M1 sqrt(X1) / (1 + n gamma M1^2),
//   p2 / p1 = (1 + n gamma M1^2) / (1 + gamma M2^2),
//   Pt2 / Pt1 = (p2 / p1) (X2 / X1)^(gamma / (gamma - 1)),
// the subsonic root (it exists for every M1 <= 1). At low Mach
// 1 - Pt2/Pt1 = (1 - n)^2 gamma M1^2 / 2 (1 + O(M1^2)).
struct BordaCarnotStep {
    double M2 = 0.0;
    double pt_ratio = 1.0;   // Pt2 / Pt1
    double dpt_dM1 = 0.0;
    double dpt_dgamma = 0.0;
};
BordaCarnotStep compressible_borda_carnot(double M1, double area_ratio, double gamma);

// -----------------------------------------------------------------------------
// The element's closure
// -----------------------------------------------------------------------------
//
// One branch: flow from the inlet face F_in to the outlet face F_out, fed from
// the stagnation state (Pt0, Tt0, Y) upstream. model is "sharp" or "conical"
// (length [m] sets the cone's half angle, atan((r_large - r_small) / length),
// circular faces). D_h [m] is the small section's hydraulic diameter for Re;
// 0 takes the circular sqrt(4 F_small / pi). incompressible runs the M -> 0
// limit: density Pt0 / (R Tt0), no choking.
//
// Loss at a given flow: the stagnation-pressure drop and where it is.
struct AreaChangeLoss {
    double dPt = 0.0;      // Pt0 - Pt_out [Pa], >= 0
    double zeta = 0.0;     // dPt over the small section's static 0.5 rho w^2
    double mach = 0.0;     // the closure's Mach: small section (expansion) or
                           // vena contracta (contraction); 0 if incompressible
    double Re = 0.0;       // small section, mu at the station's static T
    bool expansion = true;
    bool choked = false;   // the flow was above the sonic flow; dPt is at it
};
AreaChangeLoss area_change_loss(double m_dot_abs, double Pt0, double Tt0,
                                const std::vector<double>& Y, double F_in,
                                double F_out, const std::string& model,
                                double length, double D_h, bool incompressible);

// Flow at a given drive: the mass flow [kg/s] whose loss is Pt0 - Pt_target,
// the form a network residual wants (m - m_calc). It rises monotonically with
// the drive and saturates at the sonic flow once the closure chokes (C1: the
// loss is parametrised by the closure's Mach, whose flow has zero slope at
// M = 1), so it exists for every drive >= 0.
//
// At zero drive the raw law's slope is infinite for a loss quadratic in the
// flow, and a sudden expansion is laminar (finite slope) only below Re 8 --
// under 1e-5 Pa in a few-cm duct, far below what a network resolves; above it
// the flow goes about as drive^(1/2), on which Newton oscillates. So below its
// drive floor (zeta_inc AREA_CHANGE_HEAD_FLOOR) every branch is
// m0 h(drive / floor), m0 the flow at the floor and
// h(t) = t exp((sigma - 1)(t^2 - 1)/2): odd, with a finite positive slope at
// zero, matching value and slope at the floor for the incompressible law's
// exponent sigma there (1/2 for a quadratic loss; from the 4-1
// interpolation's Re-slope for a sudden expansion). sigma is analytic in the
// state, so the partials stay exact; the compressible exponent differs by
// O(M^2) at the floor: a slope step there of at most 1.1e-4, measured
// 2026-10-25 over ducts of 0.2 mm - 0.3 m, n 0.05 - 0.99, 1 - 20 bar.
//
// Derivatives are analytic: implicit in the closure's Mach, with gamma, R and
// mu chained to the mass fractions (viscosity_dX). dm_dY is per species, Y
// taken unnormalised (as normalize_fractions would).
struct AreaChangeFlow {
    double m_dot = 0.0;    // >= 0
    bool choked = false;
    double mach = 0.0;
    double Re = 0.0;
    double zeta = 0.0;
    double dm_dPt0 = 0.0;
    double dm_dTt0 = 0.0;
    double dm_dPt_target = 0.0;
    std::vector<double> dm_dY;
};
AreaChangeFlow area_change_flow(double Pt0, double Tt0, const std::vector<double>& Y,
                                double Pt_target, double F_in, double F_out,
                                const std::string& model, double length, double D_h,
                                bool incompressible, bool with_derivatives = true,
                                bool exit_head_lost = false);

// The network element's flow between its two nodes a (F_a, the from side) and
// b (F_b): m_dot > 0 from a to b. Each direction is area_change_flow fed from
// its own upstream node, odd in the drive Pt_a - Pt_b; the two are blended by
// w = sigmoid(drive / width), width the smaller drive floor of the two,
//   m = w m_ab + (1 - w) m_ba,
// so the flow is C1 through zero drive although the two directions are
// different closures (an expansion one way, a contraction the other). Past
// 30 floors of drive only the upwind branch is evaluated (w is 0 or 1 to
// 1e-13). The diagnostics (choked, mach, Re, zeta) are the upwind branch's.
// exit_head_lost_b / exit_head_lost_a apply static-pressure exit coupling
// (#1128) on the a -> b and b -> a branches respectively.
struct AreaChangeElementFlow {
    double m_dot = 0.0;
    bool choked = false;
    double mach = 0.0;
    double Re = 0.0;
    double zeta = 0.0;
    double dm_dPt_a = 0.0;
    double dm_dTt_a = 0.0;
    double dm_dPt_b = 0.0;
    double dm_dTt_b = 0.0;
    std::vector<double> dm_dY_a;
    std::vector<double> dm_dY_b;
};
AreaChangeElementFlow area_change_element_flow(
    double Pt_a, double Tt_a, const std::vector<double>& Y_a,
    double Pt_b, double Tt_b, const std::vector<double>& Y_b,
    double F_a, double F_b, const std::string& model, double length, double D_h,
    bool incompressible, bool with_derivatives = true,
    bool exit_head_lost_b = false, bool exit_head_lost_a = false);

// -----------------------------------------------------------------------------
// Incompressible kernels at a given density (standalone use)
// -----------------------------------------------------------------------------
//
// dP = zeta_eff m_dot |m_dot| / (2 rho F_small^2), the zeta above, with the
// direction blended by a sigmoid of m_dot / m_scale (forward F0 -> F1).
// Analytic partials in m_dot, rho and mu (mu through Re, sharp expansion only).
struct AreaChangeResult {
    double dP;       // Total pressure drop [Pa], signed
    double dS_dm;    // d(dP)/d(m_dot)  [Pa*s/kg]
    double dS_drho;  // d(dP)/d(rho)    [Pa*m^3/kg]
    double dS_dmu;   // d(dP)/d(mu)     [Pa/(Pa*s)]
};

// D_h: hydraulic diameter of the small section for Re (0 = circular).
AreaChangeResult sharp_area_change(double m_dot, double rho, double mu,
                                   double F0, double F1,
                                   double m_scale = 1e-4,
                                   double D_h = 0.0);

// length: axial length of the cone [m]; circular faces.
AreaChangeResult conical_area_change(double m_dot, double rho, double mu,
                                     double F0, double F1,
                                     double length,
                                     double m_scale = 1e-4);

} // namespace combaero

#endif // AREA_CHANGE_H
