#ifndef FRICTION_H
#define FRICTION_H

#include <string>
#include <unordered_map>

// -------------------------------------------------------------
// Darcy friction factor correlations for turbulent channel flow
// -------------------------------------------------------------
// All functions compute the Darcy friction factor f for use in:
//   Delta_P = f * (L/D) * (rho * v^2/2)
//
// Parameters:
//   Re  : Reynolds number [-] (must be > 0)
//   e_D : relative roughness eps/D [-] (must be >= 0)
//
// Valid for turbulent flow (Re > ~2300). For laminar flow, use f = 64/Re.
//
// References:
// - Colebrook (1939): Journal of the ICE, 11
// - Haaland (1983): ASME J. Fluids Eng., 105(1)
// - Serghides (1984): Chemical Engineering, Nov. 5

// Haaland correlation (1983) - explicit approximation
// 1/sqrt(f) = -1.8 * log10( (eps/D / 3.7)^1.11 + 6.9/Re )
// Accuracy: ~2-3% vs Colebrook-White
namespace haaland {
constexpr double coeff_outer = -1.8;
constexpr double coeff_roughness = 3.7;
constexpr double coeff_exponent = 1.11;
constexpr double coeff_reynolds = 6.9;
} // namespace haaland
double friction_haaland(double Re, double e_D);

// Serghides correlation (1984) - explicit approximation
// Uses Steffensen acceleration on Colebrook iteration.
// 1/sqrt(f) = A - (B-A)^2 / (C - 2B + A)
// Accuracy: <0.3% vs Colebrook-White (very accurate)
namespace serghides {
constexpr double coeff_roughness = 3.7;    // eps/D divisor (same as Colebrook/Haaland)
constexpr double coeff_reynolds_A = 12.0;  // 12/Re term in A
constexpr double coeff_reynolds_BC = 2.51; // 2.51*A/Re term in B, C
} // namespace serghides
double friction_serghides(double Re, double e_D);

// Colebrook-White equation (1939) - implicit, solved iteratively
// 1/sqrt(f) = -2 * log10( eps/D/3.7 + 2.51/(Re*sqrt(f)) )
// The reference standard for turbulent friction factor.
// Uses Haaland as initial guess, Newton-Raphson iteration.
double friction_colebrook(double Re, double e_D, double tol = 1e-10,
                          int max_iter = 20);

// Laminar-turbulent blend of the network friction kernels
// (friction_and_jacobian_{haaland,serghides,colebrook,petukhov}): f = (1 - w) 64/Re +
// w f_turb, w = (1 + tanh((Re - blend_center) / blend_width)) / 2.
// The laminar branch is Hagen-Poiseuille, f Re = 64, and carries its own
// floor, re_floor, far below any flow a solve meets, so the drop stays linear
// in m down to essentially zero flow. The turbulent fits keep the larger
// turbulent_re_floor, which only guards their log at Re = 0 (w ~ 0 there).
// It used to floor the laminar branch too: f saturated at 1 below Re ~ 64, so
// the drop went as m^2 with zero slope at zero flow, and a junction port that
// had to reverse stalled there (the stress ideal_nport family).
namespace laminar_friction {
constexpr double f_re = 64.0;          // Hagen-Poiseuille f Re, circular duct
constexpr double re_floor = 1e-6;      // laminar-branch floor
constexpr double turbulent_re_floor = 64.0;
constexpr double blend_center = 2300.0;
constexpr double blend_width = 400.0;
} // namespace laminar_friction

// Petukhov correlation (1970) - smooth channels only
// f = (coeff_a * ln(Re) - coeff_b)^(-2)
// Valid for 3000 < Re < 5x10^6
// Often used with Gnielinski/Petukhov heat transfer correlations.
// Reference: Petukhov (1970), Advances in Heat Transfer, 6, 503
namespace petukhov {
constexpr double coeff_a = 0.790; // ln(Re) coefficient
constexpr double coeff_b = 1.64;  // offset
// Smooth-pipe heat-transfer paths hold Petukhov friction at clamp_re and
// below (the fit is not meant for laminar Re, and ln(Re) is undefined at 0).
// The hold is a C1 quadratic soft-max of Re over
// [clamp_re - clamp_halfwidth, clamp_re + clamp_halfwidth]: identical to
// max(Re, clamp_re) outside that band, continuous slope through it. The
// hard max() it replaces left a 25% dNu/dRe jump in Gnielinski at Re 3000
// (#446).
constexpr double clamp_re = 3000.0;
constexpr double clamp_halfwidth = 500.0;
} // namespace petukhov
double friction_petukhov(double Re);

// Petukhov friction at the C1-clamped Reynolds number, and its exact
// Re-derivative (0 below clamp_re - clamp_halfwidth).
double friction_petukhov_clamped(double Re);
double friction_petukhov_clamped_dRe(double Re);

// Exact Re-derivative of friction_colebrook, by implicit differentiation of
// the Colebrook equation at its converged root (including the Re smoothing
// friction_colebrook applies internally). Replaces central differences in
// solver Jacobians (#448).
double friction_colebrook_dRe(double Re, double e_D);

// -------------------------------------------------------------
// Channel Roughness Database
// -------------------------------------------------------------
// Absolute roughness eps [m] for common channel materials.
// Used as e_D = eps/D in friction factor correlations.
//
// References:
// - Moody (1944): Transactions of the ASME, 66(8)
// - White (2011): Fluid Mechanics (7th ed.), Table 6.1
// - Munson et al. (2013): Fundamentals of Fluid Mechanics (7th ed.), Table 8.1
// - Crane Co. (2009): Technical Paper No. 410, Table A-24

// Get absolute roughness for a standard channel material [m].
// Material names are case-insensitive.
// Throws std::invalid_argument if material not found.
// Use standard_channel_roughness() to list available materials.
double channel_roughness(const std::string &material);


// Get all standard channel roughness values.
// Returns: map of material name -> roughness [m]
std::unordered_map<std::string, double> standard_channel_roughness();


#endif // FRICTION_H
