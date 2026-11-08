#ifndef COOLING_CORRELATIONS_H
#define COOLING_CORRELATIONS_H

#include "correlation_status.h"

#include <limits>
#include <tuple>
#include <utility>
#include <vector>

namespace combaero::cooling {


// Thermal-hydraulic performance factor (Webb & Eckert 1972)
// eta = (Nu/Nu0) / (f/f0)^(1/3)
// Compares heat transfer gain against pumping power penalty at equal velocity.
// eta > 1: surface is a net improvement; eta < 1: pressure-drop generator.
// Generic: works for ribs, dimples, pin fins, impingement, or any surface.
//
// Parameters:
//   Nu_ratio: Nu_enhanced / Nu_smooth [-]
//   f_ratio : f_enhanced  / f_smooth  [-]
//
// Returns: thermal-hydraulic performance factor eta [-]
//
// Source: Webb, R.L. & Eckert, E.R.G. (1972) Int. J. Heat Mass Transfer 15(8), 1647-1658
double thermal_performance_factor(double Nu_ratio, double f_ratio);


// -------------------------------------------------------------
// Helper Functions
// -------------------------------------------------------------

// Adiabatic wall temperature from cooling effectiveness
// Uses definition: eta = (T_hot - T_aw) / (T_hot - T_coolant)
//
// Parameters:
//   T_hot      : mainstream hot-gas temperature [K]
//   T_coolant  : coolant supply temperature [K]
//   eta        : adiabatic cooling effectiveness [-]
//
// Returns: adiabatic wall temperature [K]
double adiabatic_wall_temperature(double T_hot, double T_coolant, double eta);

// Film/effusion-cooled wall heat flux with single wall layer
// Computes q based on adiabatic wall temperature and thermal resistance network:
//   q = U * (T_aw - T_coolant)
// where U uses both convective sides and wall conduction.
//
// Parameters:
//   T_hot      : mainstream hot-gas driving temperature [K]
//                For M < 0.3: pass static temperature T_static.
//                For M > 0.3 (combustor liner, turbine cooling): pass
//                T_adiabatic_wall() from stagnation.h as T_hot.
//                Using Tt overcorrects (recovery factor r < 1, so T_aw < Tt).
//   T_coolant  : coolant supply temperature [K]
//   h_hot      : hot-side HTC [W/(m^2*K)]
//   h_coolant  : coolant-side HTC [W/(m^2*K)]
//   eta        : adiabatic cooling effectiveness [-]
//   t_wall     : wall thickness [m]
//   k_wall     : wall thermal conductivity [W/(m*K)]
//
// Returns: wall heat flux [W/m^2]
double cooled_wall_heat_flux(double T_hot, double T_coolant,
                             double h_hot, double h_coolant,
                             double eta,
                             double t_wall, double k_wall);


// -------------------------------------------------------------
// Baldauf et al. (2002) - laterally averaged film-cooling effectiveness
// -------------------------------------------------------------
//
// Baldauf, S., Scheurlen, M., Schulz, A. and Wittig, S. (2002). "Correlation
// of Film-Cooling Effectiveness From Thermographic Measurements at Enginelike
// Conditions." ASME J. Turbomachinery 124(4), 686-698. DOI 10.1115/1.1504443
// (conference version ASME GT-2002-30180).
// sources/heat_transfer/film/686_1_Baldauf_Film_2002.pdf (gitignored).
//
// A row of CYLINDRICAL, streamwise-inclined holes on a flat plate. Unlike the
// older exponential-decay correlations it is valid from the ejection point
// rather than only far downstream, and it carries the ADJACENT JET
// INTERACTION -- the lateral hole-spacing effect that drives jet lift-off --
// as a correlated parameter rather than an excluded case.
//
// eta here is the paper's Eq. (1), (T_G - T_AW)/(T_G - T_C), which is exactly
// combaero's own convention, so it feeds adiabatic_wall_temperature()
// unchanged.
//
// STRUCTURE. The correlation collapses every measured curve onto one base
// curve in scaled coordinates, then back-transforms:
//
//   x/D -> xi'      Eq. (40) inverted in closed form
//   xi' -> eta*'    Eq. (36), the turbulence-dependent base curve
//   eta*' -> eta*   Eq. (37), undoing the adjacent-jet downstream branch
//   eta* -> eta_c   Eq. (38), blowing up to this case's peak
//   eta_c -> eta    Eq. (39), backscaling out the geometry normalisation
//
// ANGLE UNITS. The paper tabulates alpha in DEGREES but every trigonometric
// function consumes RADIANS. This was not assumed: seven of the worked
// example's coefficients, spanning cos(a), cos(1.5a), cos(2.3a), cos(2.5a),
// sin(2a) and cos^0.65(a), reproduce the paper's Table 4 on radians and fail
// on degrees by 10-30%. The API therefore takes degrees and converts.
//
// A PUBLISHED INCONSISTENCY, IMPLEMENTED AS PRINTED. Eq. (31) for b_0 does
// not reproduce the paper's own Table 4:
//
//     Eq. (31) as printed -> 0.83612467
//     Table 4             -> 0.61626073        (36% apart)
//
// Every other one of the 17 checkable Table 4 coefficients reproduces to
// better than 2e-6 relative, so this is the paper's, not ours. Eq. (31) was
// confirmed from four independent channels (publisher MathML, the page at
// 600 dpi, the PDF text layer, and a reader's own transcription); no angle
// convention, unit choice or single-token variant reaches Table 4's value,
// and Table 4 is self-consistent with Eq. (32) to the digit. We implement the
// PRINTED equation, because fidelity means implementing what the paper
// states, and `b_0_override` lets a caller test the Table 4 reading without
// the library inventing a second formula -- the paper contains none. Table 4
// gives one number at the Table 3 conditions; reaching it needs
// sin(...) = -0.167 where the printed equation gives +0.470, and no sign,
// unit or angle convention makes that argument negative at alpha = 30 deg.
//
// THE EFFECT IS GOVERNED BY x/D, NOT BY M. An earlier note here said "under
// 5% below M ~ 0.5 and up to 50% at M = 2.5", which misattributes it: b_0
// feeds b_1 through Eq. (32), and b_1 is the DESCENDING-branch gradient, so
// it does nothing near the hole whatever the blowing. Measured by sweeping
// the envelope with b_0_override:
//
//   x/D <= 20            under 1% at every M
//   x/D = 200, M = 2.0   -19%
//   x/D = 400, M = 2.0   -28%
//   worst over the box   -75%  (s/D 5, alpha 90, Tu 0.0035, M 2.5, x/D 400)
//
// Practical consequence: at effusion row spacings (x/D of 6 to 9) the two
// readings differ by under 1%, so no effusion dataset can arbitrate Eq. (31).
// Only far-downstream single-row data could. See
// validation/cooling/extractions/baldauf_2002_film_effectiveness.md.
namespace baldauf2002 {

// Table 2 -- constants of the base curve fit. Valid for all geometry and
// density ratios, which is the whole point of the base-curve collapse.
constexpr double xi_0 = 9.0;
constexpr double eta_0 = 5.8;      // superseded by eta_0T, Eq. (35)
constexpr double a_star = 4.0;     // ascending-branch gradient
constexpr double b_star = 0.7;     // descending-branch gradient (see Eq. 34)
constexpr double c_star = 0.24;    // apex sharpness, fits the apex to eta* = 1

// Eq. (35) rebases eta_0 when the turbulence-dependent slope b*_T moves.
constexpr double eta_0T_c0 = 2.5;
constexpr double eta_0T_c1 = 5.8;

// The paper's stated envelope. Outside it the correlation is extrapolation:
// Tu in particular appears inside exp[2.6 Tu - 0.0012/Tu^2 - 1.76], which
// diverges as Tu -> 0.
constexpr double M_min = 0.2, M_max = 2.5;
constexpr double P_min = 1.2, P_max = 1.8;
constexpr double s_over_D_min = 2.0, s_over_D_max = 5.0;
constexpr double alpha_deg_min = 30.0, alpha_deg_max = 90.0;
constexpr double Tu_min = 0.0035, Tu_max = 0.075;

// The paper's own accuracy: overall RMS deviation of the measurements about
// the correlation, 5% at the apex and 3% on the descending branch.
constexpr double rms_deviation = 0.055;

// Eq. (31)'s b_0 as the worked example's Table 4 reports it, at the Table 3
// conditions. Kept so the discrepancy above is on the record in code, not
// only in prose. NOT used by the correlation.
constexpr double b0_table4_at_table3 = 0.61626073;

} // namespace baldauf2002

// Laterally averaged adiabatic film-cooling effectiveness downstream of one
// row of cylindrical, streamwise-inclined holes.
//
// Parameters:
//   x_over_D  : streamwise distance from the ejection point, in hole diameters
//   M         : blowing rate (rho u)_C / (rho u)_G [-]
//   P         : density ratio rho_C / rho_G [-]
//   alpha_deg : ejection angle to the SURFACE [deg], 30 to 90
//   s_over_D  : lateral hole spacing / diameter [-], 2 to 5
//   Tu        : mainstream turbulence intensity [-], e.g. 0.015 for 1.5%
//
// Returns: laterally averaged eta [-], Eq. (1) convention.
double film_effectiveness_baldauf_2002(double x_over_D, double M, double P,
                                       double alpha_deg, double s_over_D,
                                       double Tu,
                                       CorrelationStatus *status = nullptr,
                                       double b_0_override =
                                           std::numeric_limits<double>::quiet_NaN());

// Solver-facing (f, J): (eta, d eta/dM, d eta/dP).
//
// M and P are the two inputs a network solve varies -- M through the coolant
// mass flow, P through the temperature ratio. Geometry and Tu are fixed per
// element and carry no partials. Analytic via forward-mode dual numbers, not
// finite differences.
std::tuple<double, double, double> film_effectiveness_baldauf_2002_and_derivatives(
    double x_over_D, double M, double P, double alpha_deg, double s_over_D,
    double Tu, CorrelationStatus *status = nullptr,
    double b_0_override = std::numeric_limits<double>::quiet_NaN());

// -------------------------------------------------------------
// Colban, Thole and Bogard (2011) - shaped-hole film effectiveness
// -------------------------------------------------------------
//
// Colban, W.F., Thole, K.A. and Bogard, D. (2011). "A Film-Cooling
// Correlation for Shaped Holes on a Flat-Plate Surface." ASME J.
// Turbomachinery 133(1), 011002. DOI 10.1115/1.4002064.
// sources/heat_transfer/film/2011-Colban-JT-1fygyl9.pdf (gitignored).
//
// One row of FAN-SHAPED or LAIDBACK FAN-SHAPED holes on a flat plate.
// Baldauf above is for cylindrical holes; the two are not interchangeable,
// and the paper's own Fig. 7 is the reason: at M = 2.5 Baldauf's
// cylindrical correlation predicts eta ~ 0.04 where the shaped-hole data
// sit at 0.15-0.5, because a cylindrical jet has lifted off and a shaped
// one has not.
//
// Eq. (19), with xi from Eq. (18):
//
//   eta = 1 / (P/t + C1 M^C2 xi^C3)
//   xi  = (4/pi) (x/D) (P/D) / (M AR)
//
// D is the METERING (inlet) diameter, M is referenced to the metering
// velocity (Eq. 7), AR = A_exit/A_inlet (Eq. 13), t is the hole breakout
// width and P the lateral pitch. At x/D = 0 the correlation returns t/P,
// the adiabatic mixing limit of Eq. (12) -- a boundary condition built into
// the form, not fitted.
//
// eta is the laterally averaged (T_aw - T_jet)/(T_inf - T_jet) of Eqs. (1)
// and (2): the cooled fraction of the driving difference, the same quantity
// as Baldauf's and as adiabatic_wall_temperature() consumes.
//
// FITTED, NOT DERIVED. C1-C3 are a linear regression over 57 data sets from
// four studies (Table 1), all at alpha = 30 deg, DR 1.7-2.0, Tu 0.2-4% and
// external Mach 0.3. Density ratio, surface angle and turbulence are absent
// from Eq. (19) BECAUSE the data never varied them, not because they do not
// matter -- so they are range-checked here even though they do not enter
// the value. See validation/cooling/extractions/colban_2011_shaped_film.md.
namespace colban2011 {

// Table 3, Eq. (19).
constexpr double C1 = 0.1721;
constexpr double C2 = -0.2664;
constexpr double C3 = 0.8749;

// Table 4: the validity limits the paper derives from where its own fit
// keeps R^2 > 0.7. The third is the paper's "coverage intensity", the
// combination it found most limiting -- AR / (M P/D) above 1.17 is where
// adjacent jets interact.
constexpr double M_min = 0.2, M_max = 2.5;
constexpr double t_over_P_min = 0.31, t_over_P_max = 0.65;
constexpr double coverage_intensity_min = 0.17, coverage_intensity_max = 1.17;

// The fitting data's span (Section 5, Table 1). Not stated as limits by the
// paper, but nothing outside them was fitted.
constexpr double AR_min = 2.5, AR_max = 4.7;
constexpr double P_over_D_min = 4.0, P_over_D_max = 8.0;
// Every data set was at 30 deg (Section 3.4) and density ratio 1.7-2.0
// (Section 3.8); neither enters Eq. (19).
constexpr double alpha_deg_fitted = 30.0;
constexpr double density_ratio_min = 1.7, density_ratio_max = 2.0;

} // namespace colban2011

// Laterally averaged adiabatic film-cooling effectiveness downstream of one
// row of shaped (fan or laidback fan) holes.
//
// Parameters:
//   x_over_D      : distance downstream of the hole breakout, in METERING
//                   diameters, >= 0
//   M             : blowing ratio at the metering section [-], 0.2 to 2.5
//   area_ratio    : A_exit / A_inlet [-], >= 1 (2.5 to 4.7 fitted)
//   P_over_D      : lateral pitch / metering diameter [-] (4 to 8 fitted)
//   t_over_P      : coverage ratio, breakout width / pitch [-], 0.31 to 0.65
//   alpha_deg     : surface angle [deg] -- range check only, fitted at 30
//   density_ratio : rho_jet / rho_ext [-] -- range check only, 1.7 to 2.0
//
// Returns: laterally averaged eta [-].
double film_effectiveness_colban_2011(double x_over_D, double M,
                                      double area_ratio, double P_over_D,
                                      double t_over_P, double alpha_deg,
                                      double density_ratio,
                                      CorrelationStatus *status = nullptr);

// Solver-facing (f, J): (eta, d eta/dM). M is the only input a network
// solve varies -- geometry is fixed per element, and the density ratio does
// not enter Eq. (19). Closed form: M enters as M^(C2 - C3) once xi is
// substituted, so the derivative is exact.
std::tuple<double, double> film_effectiveness_colban_2011_and_derivatives(
    double x_over_D, double M, double area_ratio, double P_over_D,
    double t_over_P, double alpha_deg, double density_ratio,
    CorrelationStatus *status = nullptr);

// -------------------------------------------------------------
// Marek and Tacina (1975) - slot film effectiveness against turbulence
// -------------------------------------------------------------
//
// Marek, C.J. and Tacina, R.R. (1975). "Effect of Free-Stream Turbulence on
// Film Cooling." NASA TN D-7958.
// sources/heat_transfer/film/oct23/19750017021.pdf (gitignored).
//
// A two-dimensional SLOT film (a continuous slot across the span, as in a
// combustor liner's cooling ring), not a row of holes: Baldauf and Colban
// above do not apply to it. Juhasz-Marek turbulent mixing, Eq. (4):
//
//   eta = 1 / (1 + C_m x / (M S))
//
// with the mixing coefficient from the free-stream turbulence, Eq. (5):
//
//   C_m = (0.25 +- 0.15) Tu / 100     (Tu in percent, as the report writes it)
//       = (0.25 +- 0.15) Tu           (Tu a fraction, as this kernel takes it)
//
// x is the distance downstream of the slot lip, S the slot height, M =
// rho_s U_s / (rho_H U_H) the blowing ratio and Tu the axial free-stream
// turbulence intensity, here a FRACTION (0.14 = 14%) like every other
// turbulence input in the library (#773; the report gives percent).
// eta = (T_H - T_wall)/(T_H - T_s), the
// cooled fraction of the driving difference as for the hole films. Eq. (4)
// carries no cp ratio; Juhasz and Marek (1971, TN D-6360) wrote one in,
// and Eq. (5) was fitted without it.
//
// FITTED, NOT DERIVED. Eq. (5) is a line through the four per-blockage
// averages of C_m (Tu 7, 14, 23, 35%). The entrainment rate C_m rho_H U_H
// already scales with the hot-stream velocity, so the form has no Re;
// Tu is its only turbulence input. Every test is co-flowing film.
// See validation/cooling/extractions/marek_tacina_1975_d7958.md.
namespace marek_tacina1975 {

// Eq. (5): C_m per unit of turbulence (Tu a fraction), and its stated
// +- band. The report's 0.25 per percent is 0.25 per unit fraction x 100
// / 100: the same line.
constexpr double C_m_per_Tu = 0.25;
constexpr double C_m_per_Tu_band = 0.15;

// The data's span (Table I and p. 9): axial Tu (a fraction, 7 to 35%), blowing ratio
// (m_s/m_H)(H/S) with H = 10.5 cm and S = 0.63 cm, and x/S at the nine
// stations. Not stated as limits by the report; nothing outside was run.
constexpr double Tu_min = 0.07, Tu_max = 0.35;
constexpr double M_min = 0.41, M_max = 1.87;
constexpr double x_over_S_max = 33.8;

// The authors' maximum eta deviation from measurement (p. 6): +-5 K
// hot-gas variation and up to 7 K axial wall conduction.
constexpr double eta_measurement_band = 0.04;

} // namespace marek_tacina1975

// Slot film-cooling effectiveness downstream of a 2-D slot.
//
// Parameters:
//   x_over_S       : distance downstream of the slot lip in slot heights, >= 0
//   M              : blowing ratio rho_s U_s / (rho_H U_H) [-], > 0
//   Tu             : axial free-stream turbulence intensity, a FRACTION in
//                    [0, 1] (0.14 = 14%; 0.07 to 0.35 run). Above 1 throws:
//                    that is a percent passed by mistake.
//   C_m_multiplier : the user's rig-matching factor on Eq. (5)'s C_m [-],
//                    > 0. 1 is the report's line; the band is 0.4 to 1.6.
//                    Never scored.
//
// Returns: eta [-]; 1 at the lip.
double slot_film_effectiveness_marek_tacina_1975(
    double x_over_S, double M, double Tu, double C_m_multiplier = 1.0,
    CorrelationStatus *status = nullptr);

// Solver-facing (f, J): (eta, d eta/dM). M is the only input a network
// solve varies; geometry and Tu are fixed per element. Closed form.
std::tuple<double, double> slot_film_effectiveness_marek_tacina_1975_and_derivatives(
    double x_over_S, double M, double Tu, double C_m_multiplier = 1.0,
    CorrelationStatus *status = nullptr);

// -------------------------------------------------------------
// Multi-row film superposition
// -------------------------------------------------------------
//
// Sellers, J.P. (1963). "Gaseous film cooling with multiple slot injection."
// AIAA Journal 1(9), 2154-2156. As restated in:
//
//   Gao, Z., Qiu, T., Liu, P., Ding, S., Li, Z., Cheng, R. and Yuan, Q.
//   (2025). "A Study on the Film Superposition Method for the Multi-Row Film
//   Cooling of the Turbine Outer Ring." Processes 13, 143, Eqs. (1)-(11).
//   sources/heat_transfer/film/processes-13-00143-v2.pdf (open access).
//
// WHY A CORRECTION IS NEEDED. Sellers assumes the rows are independent. Gao
// measures what that costs: "the Sellers method accumulates prediction
// errors as the number of hole rows increases, leading to an OVERESTIMATION
// of the cooling efficiency." That is worst exactly where effusion lives --
// many closely spaced rows -- so a film module built for a few rows and then
// reused for effusion would be wrong in the regime it is needed most.
//
// Gao's fix is a per-row mainstream temperature correction alpha. It is not
// an invented knob: Eqs. (3) and (4) are an energy balance on the mainstream
// entrained into the boundary layer at each injection, giving
//
//   (T_g - T'_aw)/(T_g - T_aw) = C (m_c/m_g) / (C (m_c/m_g) + 1)
//
// so alpha is the fraction of the film's temperature deficit that survives
// mixing on the way to the next row. alpha = 1 recovers Sellers exactly.
//
// GAO'S COEFFICIENTS, AND WHY THEY ARE STILL NOT A DEFAULT. Eq. (5) gives
//
//   alpha_i = a r / (a r + 1) + b,     r = m_coolant / m_mainstream
//
// and section 4.3.1 DOES print the fitted values: "The empirical
// coefficients in Equation (5), a and b, were determined to be 12 and
// 0.9465, respectively." They are recorded below as gao2025::a_case1 and
// b_case1. (An earlier note here said they were never published; that was
// wrong, and it was wrong because a regex looked for "a = 12" while the
// paper states it in prose.)
//
// They are NOT defaulted, because they cannot reach what a tighter plate
// needs -- and that follows from the constants alone, with no rig
// arithmetic. Since a r / (a r + 1) >= 0 for a > 0 and r >= 0,
//
//     alpha >= b = 0.9465   for EVERY r,
//
// so the published pair can never damp a row by more than 5.35% while
// staying at or below 1 (above 1 the superposition refuses -- alpha > 1
// would create coolant). Murray & Ireland's 5.75 D staggered plate needs a
// per-row alpha of about 0.85 at low blowing and 0.69 near M = 1, both
// BELOW that floor. No choice of r reaches them.
//
// This argument is deliberately independent of how r is scaled: Gao's test
// section dimensions are not stated in the text, so r cannot be put on
// their scale, and any claim that depended on doing so would be
// unverifiable.
//
// The likely reason is structural: Eq. (5) carries NO streamwise-spacing
// term, while Murray shows streamwise spacing to be the dominant driver --
// tripling it to 17.25 D restores plain superposition to under 10% error.
// Gao calibrated on a plate with 10.5 d streamwise spacing; Murray's
// effective spacing is 2.875 D, 3.7x tighter. Two plates with the same
// coolant fraction and different row spacing get the same alpha from this
// form, which cannot be right if spacing dominates.
//
// So alpha stays a caller input defaulting to 1 (plain Sellers). This is
// the tuner slot the validation policy reserves: matching a specific rig is
// the user's job and the harness never scores it. See
// validation/cooling/extractions/gao_alpha_fit_on_murray.md.
namespace film_superposition {

// alpha = 1 is plain Sellers. Anything below it damps the accumulated
// effectiveness, which is the direction Gao's measurements require.
constexpr double alpha_sellers = 1.0;

// Gao's own fitted coefficients for Eq. (5), section 4.3.1. Recorded so the
// published values are in code rather than only in prose, NOT so they can be
// used as defaults -- see the note above for why they do not transfer.
// Calibrated on Case 1 (d = 1.2 mm, spanwise pitch 3.5 d, streamwise spacing
// 10.5 d) at blowing ratios 0.3 and 1.0.
constexpr double gao_a_case1 = 12.0;
constexpr double gao_b_case1 = 0.9465;

} // namespace film_superposition

// Sellers superposition, Gao Eq. (1):
//
//   eta = eta_1 + sum_{i=2..n} eta_i prod_{j<i} (1 - eta_j)
//       = 1 - prod_i (1 - eta_i)
//
// The two forms are algebraically identical; the product form is used
// because it is O(n) and cannot accumulate the pairwise rounding the sum
// form does. Known to OVERESTIMATE as the row count grows.
double film_superposition_sellers(const std::vector<double>& eta_rows);

// (eta, d eta / d eta_i) for the above. d eta/d eta_i = prod_{j != i}
// (1 - eta_j), computed without division so a fully effective row
// (eta_j = 1) does not produce a NaN.
std::pair<double, std::vector<double>> film_superposition_sellers_and_gradient(
    const std::vector<double>& eta_rows);

// Gao Eq. (7), Sellers with the per-row mainstream temperature correction:
//
//   eta = sum_i [ eta_i prod_{j=i..n-1} alpha_j prod_{k=i+1..n} (1 - eta_k) ]
//
// alpha_between_rows has n-1 entries for n rows: alpha_j is the correction
// applied between row j and row j+1. Passing all ones reproduces
// film_superposition_sellers exactly.
double film_superposition_corrected(const std::vector<double>& eta_rows,
                                    const std::vector<double>& alpha_between_rows);

// (eta, d eta / d eta_i) for Eq. (7).
std::pair<double, std::vector<double>> film_superposition_corrected_and_gradient(
    const std::vector<double>& eta_rows,
    const std::vector<double>& alpha_between_rows);

// Gao Eq. (5), the mainstream temperature correction. Both the form and the
// coefficients are published (gao_a_case1, gao_b_case1 above), but a and b
// stay REQUIRED rather than defaulted: the published pair is calibrated on
// one plate and returns alpha > 1 on a tighter one -- see the note above.
//
//   mass_flow_ratio : m_coolant / m_mainstream for this row [-]
double mainstream_temperature_correction(double mass_flow_ratio, double a,
                                         double b);

// Gao Eq. (9): equivalent slot width of a row of holes, s = A_hole / pitch.
// This is what lets one row's measured distribution stand in for another
// spacing, via the X/(M s) scaling of Eq. (8).
double equivalent_slot_width(double hole_area, double pitch);

// Gao Eq. (10): equivalent blowing ratio, M_e = M_0 A_0 / A_e. Normalises a
// hole count and spacing onto the baseline single-row configuration.
double equivalent_blowing_ratio(double M_baseline, double area_baseline,
                                double area_equivalent);

// Panel-averaged adiabatic film effectiveness of an effusion plate: Baldauf
// (2002) per row, Sellers-superposed, averaged over the panel's rows. A
// point at row n sees rows 1..n upstream at distances j * pitch_x, so
//
//   eta_panel = (1/N) sum_{n=1..N} Sellers{ eta_Baldauf(j pitch_x / D), j = 1..n }
//
// This is the alternative closure the Andrews 88-GT-290 scoring reports
// beside the validated two-resistance one (validation/cooling/
// effusion_overall_runner.py): the data there REFUSE it offered alone, and
// the missing physics is a gas-side augmentation. It exists to be selected
// explicitly and reported, not as a default.
//
// `extrapolated` is set when any row is outside Baldauf's envelope (M, P,
// s/D, alpha, Tu); for full-coverage effusion it usually is, through s/D and M.
struct EffusionPanelFilm {
  double eta = 0.0;
  bool extrapolated = false;
};

EffusionPanelFilm effusion_panel_film_effectiveness(int n_rows,
                                                    double pitch_x_over_D,
                                                    double s_over_D, double M,
                                                    double density_ratio,
                                                    double alpha_deg, double Tu);

// Film a downstream wall sees: the laterally averaged adiabatic
// effectiveness of n_rows rows of holes, Sellers-superposed, then averaged
// over the wall's streamwise span [x_start_over_D, x_end_over_D]. Both are
// measured downstream of the LAST row, in hole diameters; row j (0 = last)
// sits a further j * row_pitch_over_D upstream:
//
//   eta(x)  = 1 - prod_j (1 - eta_row(x + j row_pitch))
//   eta_bar = 1/(x_end - x_start) int eta(x) dx
//
// x_start == x_end is the point value. The span integral is fixed composite
// Gauss-Legendre in ln(1 + x/D), where the near-hole decay is steep -- a
// quadrature, not a finite difference, so the derivatives are exact for it.
//
// This is what a ThermalWall's `film` reads (#386): the per-row kernels are
// the published ones, and the superposition and averaging add no fitted
// constant. Sellers is known to overestimate as the row count grows (see
// film_superposition_sellers above).
struct FilmSpanAverage {
  double eta = 0.0;      // eta_bar [-]
  double deta_dM = 0.0;  // d eta_bar / d M
  double deta_dP = 0.0;  // d eta_bar / d density ratio (0 for Colban)
  bool extrapolated = false;
};

// Cylindrical holes: Baldauf (2002) per row. s_over_D is the lateral pitch.
FilmSpanAverage film_span_average_baldauf_2002(
    int n_rows, double row_pitch_over_D, double x_start_over_D,
    double x_end_over_D, double M, double density_ratio, double alpha_deg,
    double s_over_D, double Tu);

// Shaped holes: Colban (2011) per row. D is the METERING diameter and x is
// measured from the breakout, as in film_effectiveness_colban_2011.
FilmSpanAverage film_span_average_colban_2011(
    int n_rows, double row_pitch_over_D, double x_start_over_D,
    double x_end_over_D, double M, double area_ratio, double P_over_D,
    double t_over_P, double alpha_deg, double density_ratio);

// 2-D slot: Marek and Tacina (1975), D-7958 Eqs. (4)-(5). S is the slot
// height and x is measured from the slot lip; one slot, no superposition.
// deta_dP is 0: the density ratio does not enter Eq. (4).
FilmSpanAverage slot_film_span_average_marek_tacina_1975(
    double x_start_over_S, double x_end_over_S, double M, double Tu,
    double C_m_multiplier = 1.0);

// Film from rows in GROUPS, each with its own blowing ratio, density ratio
// and coolant temperature: an effusion liner's panels along one duct (#695).
// Groups are listed upstream first; group g holds n_rows[g] rows, and every
// row (across groups) is row_pitch_over_D from the next. x is measured
// downstream of the LAST row of the LAST group, as above.
//
// With per-row coolant temperatures, superposition is a sequence of layers:
// past row j the adiabatic wall sits at (1 - eta_j) T_before + eta_j T_c,j.
// Unrolled from the mainstream T_ref:
//
//   T_aw(x) = T_ref prod_j (1 - eta_j) + sum_j W_j T_c,j,
//   W_j     = eta_j prod_{k downstream of j} (1 - eta_k),
//
// Sellers exactly when every T_c,j is equal (sum_j W_j = 1 - prod(1 - eta)).
// Averaged over the span:
//
//   T_aw,bar = T_ref (1 - eta_bar) + sum_g weight[g] T_c,g,
//   sum_g weight[g] = eta_bar.
//
// d_weight_dM[g][h] = d weight[g] / d M_h (likewise d_weight_dP for the
// density ratio, deta_dM[h] / deta_dP[h] for eta_bar). Same quadrature as
// film_span_average_*: exact derivatives for it.
struct FilmSpanAverageGroups {
  double eta = 0.0;
  std::vector<double> weight;
  std::vector<double> deta_dM;
  std::vector<double> deta_dP;
  std::vector<std::vector<double>> d_weight_dM;
  std::vector<std::vector<double>> d_weight_dP;
  bool extrapolated = false;
};

// Cylindrical holes, Baldauf (2002) per row; M[g] and density_ratio[g] per
// group, the hole geometry (alpha, s/D, Tu) shared.
FilmSpanAverageGroups film_span_average_baldauf_2002_groups(
    const std::vector<int> &n_rows, double row_pitch_over_D,
    double x_start_over_D, double x_end_over_D, const std::vector<double> &M,
    const std::vector<double> &density_ratio, double alpha_deg,
    double s_over_D, double Tu);

// -------------------------------------------------------------
// Effusion plate INTERNAL heat transfer
// -------------------------------------------------------------
//
// Andrews, Alikhanizadeh, Asere, Hussain, Khoshkbar Azari and Mkpadi (1986),
// "Small Diameter Film Cooling Holes: Wall Convective Heat Transfer",
// ASME 86-GT-225. See
// validation/cooling/extractions/andrews_effusion_internal_h.md.
//
// A coolant hole cools its wall in TWO places, and the paper's central point
// is that the first dominates:
//
//   1. the APPROACH flow over the coolant-side plate surface, converging
//      into the hole (Sparrow's multi-hole correlation, Eq. 15);
//   2. the THROAT, a short tube with a sharp-edged entry (Mills' tabulated
//      data, Eq. 12 with the entry-length factor R_Nu of Eqs. 13 and 14).
//
// Andrews sums them -- "the authors have treated the wall heat transfer as
// the summation of Equations 13 or 14 and 15" -- after putting both on the
// same Nusselt definition, which is Eq. 18's job. Both Nu below are
// therefore referenced to the HOLE DIAMETER and the HOLE INTERNAL AREA
// pi D L, so they add directly.
//
// CONVERTING TO A PLATE-AREA COEFFICIENT. Effusion elements want a
// coefficient per unit plate area, which is what the source's own Fig. 8
// plots. That is
//
//   h_plate = Nu k / D * A_h / A,   A_h = pi D L,  A = X^2 - pi D^2 / 4
//
// and the A_h/A factor is not optional: for Andrews' plate C it is 3.46, so
// omitting it overstates the coefficient by that factor. The conversion is
// left to the caller because it needs no correlation -- only geometry.
namespace andrews1986 {

// Eq. 18's leading constant, 0.881, is Sparrow's multi-hole approach
// correlation Nu_l = 0.881 Re^0.476 Pr^(1/3) rebased from the surface
// dimension l (Eq. 16) onto the hole internal area. Independent of
// pitch-to-diameter ratio, per the paper.
constexpr double sparrow_coefficient = 0.881;
constexpr double sparrow_re_exponent = 0.476;

// Mills' throat data, Eq. 12, is Dittus-Boelter in form with a tabulated
// entry-length multiplier.
constexpr double mills_coefficient = 0.023;
constexpr double mills_re_exponent = 0.8;

// R_Nu is curve-fitted in two branches that meet at L/D = 2. They agree
// there to 1.3e-3, which is how the transcription was checked.
constexpr double mills_branch_L_over_D = 2.0;

} // namespace andrews1986

// Mills' entry-length factor R_Nu for a sharp-edged tube entry, Andrews
// Eqs. (13) and (14). A short hole never reaches a developed profile, so
// R_Nu > 1; it decays to 1 as the hole lengthens.
//
//   L/D <= 2 :  R_Nu = 0.13 z^3 - 0.75 z^2 + 1.04 z + 2.24,       z = L/D
//   L/D >  2 :  R_Nu = 1 - 49.3 w^4 + 58.6 w^3 - 26.5 w^2 + 7.48 w, w = D/L
double mills_entry_length_factor(double L_over_D);

// Sparrow's hole-approach contribution, Andrews Eq. (18): the coolant-side
// plate surface converging into the hole, expressed on the hole's own
// Nusselt definition so it can be summed with the throat term.
//
//   Nu = 0.881 Re^0.476 Pr^(1/3) * X / (pi L)
//
//   Re        : on the HOLE diameter [-]
//   X_over_L  : hole pitch / hole length [-]
double effusion_approach_nusselt(double Re, double Pr, double X_over_L);

// Mills' short-hole throat contribution, Andrews Eq. (12).
//
//   Nu = 0.023 Re^0.8 Pr^(1/3) R_Nu(L/D)
double effusion_throat_nusselt(double Re, double Pr, double L_over_D);

// Andrews Eq. (19): the two summed, both on the hole diameter and hole
// internal area.
//
//   Nu = (0.881 (X/pi L) Re^0.476 + 0.023 Re^0.8 R_Nu) Pr^(1/3)
//
// Multiply by k/D for a hole-area coefficient, and by A_h/A again for the
// plate-area coefficient an effusion element needs.
double effusion_internal_nusselt(double Re, double Pr, double X_over_L,
                                 double L_over_D);

// -------------------------------------------------------------
// Turbulent flat plate -- the unblown reference of full-coverage film data
// -------------------------------------------------------------
//
// Crawford, M.E., Choe, H., Kays, W.M. and Moffat, R.J. (1976). "Full-
// coverage film cooling heat transfer study -- summary of data for normal-
// hole injection and 30 deg slant-hole injection." NASA CR-2648 (Stanford
// HMT-19). sources/heat_transfer/film/ (gitignored).
//
// The Stanford programme (CR-2642, CR-2648, CR-2786, CR-3219) reports every
// full-coverage Stanton number against ONE unblown reference, printed under
// each data table ("Stanton number ratio based on ...", Note 2 of the
// nomenclature): the turbulent flat-plate boundary layer
//
//   St_o Pr^0.4 = 0.0295 Re_x^-0.2 [1 - (xi/(x - x_vo))^0.9]^(-1/9)
//
// with Re_x = (x - x_vo) U/nu from the hydrodynamic virtual origin x_vo and
// xi the unheated starting length. STCR = St/St_o is the measured gas-side
// augmentation of a blown plate, i.e. the h_f/h_0 an effusion element's
// `gas_augmentation` stands for -- so this is the baseline a published
// augmentation ratio is referenced to. Here xi = 0: the thermal boundary
// layer starts with the velocity one (the unheated-starting-length bracket
// is 1).
//
// LOCAL:  Nu_x = St_o Re_x Pr = 0.0295 Re_x^0.8 Pr^0.6
//
// PANEL AVERAGE over x in [x0, x0 + L], integrated in closed form
// (h_x ~ x^-0.2):
//
//   Nu_L = h_bar L / k = (0.0295 / 0.8) Pr^0.6 (Re_x1^0.8 - Re_x0^0.8)
//        = 0.036875 Pr^0.6 Re_L^0.8 [(1 + s)^0.8 - s^0.8],   s = x0 / L
//
// At s = 0 it is the plate-from-the-leading-edge average (0.037 Re_L^0.8 in
// the textbook Pr^(1/3) form).
//
// WHY THIS FORM AND NOT THE TEXTBOOK ONE. The common 0.0296 Re_x^0.8 Pr^(1/3)
// (Colburn) is not held in the source register, and at combustion-gas Pr
// ~0.74 it is 8% above this one -- Pr^0.6 against Pr^(1/3). An augmentation
// ratio only carries over against the reference it was measured on, and
// St/St_o is the one quantity the held full-coverage data give, so the
// reference is theirs. The unheated-starting-length bracket is left out: a
// user rarely knows a separate thermal start, and it only raises h.
//
// ENVELOPE. The Stanford data span Re_x 1.5e5 to 5e6 (CR-2786's summary),
// air (Pr ~0.7), low free-stream turbulence. Below 1.5e5 a real boundary
// layer may be transitional; the formula is still smooth there and is
// evaluated, and `extrapolated` says so -- it never refuses.
namespace crawford1976 {

constexpr double st0_coefficient = 0.0295;
constexpr double st0_re_exponent = -0.2;
constexpr double st0_pr_exponent = -0.4;
constexpr double Re_x_min = 1.5e5;
constexpr double Re_x_max = 5.0e6;

} // namespace crawford1976

// Local turbulent flat-plate Nusselt number on x, and its derivative.
//
//   Re_x : x-Reynolds number from where the boundary layer starts [-], > 0
//   Pr   : gas Prandtl number [-], > 0
//
// Returns (Nu_x, dNu_x/dRe_x) = (0.0295 Re_x^0.8 Pr^0.6, 0.8 Nu_x/Re_x).
// `status` is set to Extrapolated outside Re_x [1.5e5, 5e6].
std::tuple<double, double> flat_plate_turbulent_nusselt_and_derivative(
    double Re_x, double Pr, CorrelationStatus *status = nullptr);

struct FlatPlatePanelNusselt {
    double Nu = 0.0;          // panel-averaged, on the panel length L [-]
    double dNu_dRe_L = 0.0;   // d Nu / d Re_L at fixed x0/L [-]
    double Re_x_start = 0.0;  // Re_x at the panel's leading edge [-]
    double Re_x_end = 0.0;    // Re_x at its trailing edge [-]
    bool extrapolated = false;  // any of [Re_x_start, Re_x_end] outside the envelope
};

// Panel-averaged turbulent flat-plate Nusselt number over [x0, x0 + L].
//
//   Re_L           : Reynolds number on the panel length, rho U L / mu [-], > 0
//   Pr             : gas Prandtl number [-], > 0
//   x_start_over_L : x0 / L, where the boundary layer starts upstream of the
//                    panel's leading edge, in panel lengths [-], >= 0
//
// h_bar = Nu k / L. `extrapolated` is set when the panel's Re_x range leaves
// [1.5e5, 5e6] anywhere (always at x0 = 0, where Re_x starts from zero).
FlatPlatePanelNusselt flat_plate_turbulent_panel_nusselt(double Re_L, double Pr,
                                                         double x_start_over_L);

// Solver-facing (f, J): (Nu_L, dNu_L/dRe_L) at fixed x0/L. Closed form:
// Nu_L ~ Re_L^0.8, so the derivative is 0.8 Nu_L / Re_L exactly.
std::tuple<double, double> flat_plate_turbulent_panel_nusselt_and_derivative(
    double Re_L, double Pr, double x_start_over_L,
    CorrelationStatus *status = nullptr);

} // namespace combaero::cooling

#endif // COOLING_CORRELATIONS_H
