#ifndef ORIFICE_H
#define ORIFICE_H

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <tuple>

// -------------------------------------------------------------
// Orifice Discharge Coefficient (Cd) Correlations
// -------------------------------------------------------------
//
// This module provides Cd correlations for various orifice geometries:
//   - Sharp thin-plate orifices (ISO 5167, Reader-Harris/Gallagher)
//   - Thick-plate orifices (t/d correction per Idelchik)
//   - Rounded-entry orifices (r/d based, Idelchik)
//   - User-defined (tabulated or custom function)
//   - Discharge holes -- bleed, film, effusion (McGreehan & Schotsch 1988)
//
// TWO FAMILIES, TWO SELECTORS. MeteringCdCorrelation covers orifices that sit
// in a pipe and are described by beta = d/D; DischargeCdCorrelation covers
// holes in a wall, described by L/d, r/d and the approach crossflow. See the
// comment on DischargeCdCorrelation for why they are not one enum.
//
// The discharge coefficient Cd relates actual to ideal flow:
//   mdot_actual = Cd * mdot_ideal
//   mdot_ideal  = A * sqrt(2 * rho * dP)
//
// Compressible flow note:
//   For sharp-edged orifices, Cd is weakly dependent on Mach number and
//   remains dominated by geometry and Reynolds number. These incompressible
//   correlations provide adequate Cd values into the choked-flow regime when
//   combined with compressible mass-flow relations (see compressible.h).
//
//   WHO OWNS THE EXPANSION TERM. Two routes exist and they must not be
//   combined:
//     - regime='compressible' solves the isentropic nozzle exactly via
//       combaero::nozzle_flow, choked branch included. Apply NO expansion
//       factor on top of it.
//     - the incompressible form m_dot = Cd A sqrt(2 rho dP) carries no
//       compressibility at all, and is extended to finite pressure ratio by
//       orifice::mcgreehan_schotsch::expansion_factor (Eqs. 4-7).
//   Eq. (5) reproduces nozzle_flow to 0.008%, so using both double-counts.
//
// References:
// - ISO 5167-2:2003 - Orifice plates
// - Reader-Harris & Gallagher (1998) - NEL/ASME correlation
// - Idelchik, I.E. - Handbook of Hydraulic Resistance (3rd ed.)
// - Bohl, W. - Technische Stroemungslehre (declared, NOT implemented)
// - Spink, L.K. - Principles and Practice of Flow Meter Engineering

// -------------------------------------------------------------
// Orifice geometry
// -------------------------------------------------------------

enum class OrificeType {
    SharpThinPlate,   // ISO 5167-type sharp-edged thin plate
    ThickPlate,       // Finite thickness plate (t/d > 0)
    RoundedEntry,     // Rounded inlet edge (r/d > 0)
    Conical,          // Conical (beveled) inlet
    QuarterCircle,    // Quarter-circle (quadrant) edge
    UserDefined       // Custom correlation
};

// Pressure-tap arrangement of a metering orifice (ISO 5167-2:2003, 5.3.2.1).
//
// A MEASUREMENT choice, not physics: it fixes which pressure difference C is
// referenced to, i.e. what a meter at those taps reads. It does not change
// what the plate does to the flow -- the permanent loss is the same plate
// whichever taps report it (#645).
enum class OrificeTap {
    Corner,  // L1 = L2' = 0
    D_D2,    // D and D/2 taps: L1 = 1, L2' = 0.47
    Flange   // L1 = L2' = 25.4 mm / D
};

struct OrificeGeometry {
    double d = 0.0;       // Orifice bore diameter [m]
    double D = 0.0;       // Pipe diameter [m]
    // Taps the metering correlations reference C to. Flange is the default,
    // and was the only arrangement before the tap became selectable.
    OrificeTap tap = OrificeTap::Flange;
    double t = 0.0;       // Plate thickness [m] (for thick plate)
    double r = 0.0;       // Inlet edge radius [m] (for rounded entry)
    double bevel = 0.0;   // Bevel angle [rad] (for conical)

    // Derived quantities
    double beta() const;          // Diameter ratio d/D [-]
    double area() const;          // Orifice area [m^2]
    double t_over_d() const;      // Thickness ratio t/d [-]
    double r_over_d() const;      // Radius ratio r/d [-]

    // Validation
    bool is_valid() const;
};

// -------------------------------------------------------------
// Flow state at orifice
// -------------------------------------------------------------

struct OrificeState {
    double Re_D = 0.0;    // Pipe Reynolds number (based on D) [-]
    double dP = 0.0;      // Differential pressure across orifice [Pa]
    double rho = 0.0;     // Fluid density [kg/m^3]
    double mu = 0.0;      // Dynamic viscosity [Pa*s]

    // Derived quantities
    double Re_d(double beta) const;  // Orifice Reynolds number (based on d)
};

// -------------------------------------------------------------
// Correlation identifiers
// -------------------------------------------------------------

// Cd correlations for a NORMED measurement orifice -- a standardised plate in
// a pipe, where Cd is referenced to the TAPPING differential and every member
// is a function of beta = d/D.
//
// The thick-plate and rounded-entry members that used to live here were
// removed in favour of DischargeCdCorrelation::Idelchik1966*: they computed
// an ISO 5167 Cd and multiplied it by a correction, but a thick-edged or
// rounded orifice is not the normed device, so the ISO base does not apply --
// and Idelchik gives the complete zeta directly, with no ISO base needed.
enum class MeteringCdCorrelation {
    // Sharp thin-plate correlations
    ReaderHarrisGallagher,  // ISO 5167-2:2003 / ASME MFC-3M (most accurate)

    // Special
    Constant,               // Fixed Cd value (for testing/simple cases)
    UserFunction            // User-provided function
};

// Cd families for a DISCHARGE hole -- a bleed, film or effusion hole that
// dumps coolant out of its circuit.
//
// WHY THIS IS A SEPARATE SELECTOR from MeteringCdCorrelation, rather than
// more members on it. The two families do not take the same inputs. A
// metering orifice sits in a pipe, and every correlation above is a function
// of beta = d/D and the pipe Reynolds number. A discharge hole has no pipe to
// form beta with; its Cd is a function of L/d, r/d and the approach
// crossflow, none of which OrificeGeometry/OrificeState can express (there is
// no crossflow term in OrificeState at all). One enum over both would have
// made every caller pass a meaningless D and silently drop the crossflow.
enum class DischargeCdCorrelation {
    // McGreehan & Schotsch (1988), the composite chain of Eqs. (8)-(17).
    // Sharp-to-rounded inlet, finite L/d, inlet-side crossflow.
    McGreehanSchotsch1988,

    // Idelchik (1966), Section IV: a hole in a LARGE WALL (F1 = F2 = inf --
    // plenum to plenum, which is the effusion-plate geometry). zeta is
    // referenced to the hole velocity w0 and DH is the full permanent loss
    // because there is no downstream recovery, so Cd = 1/sqrt(zeta) is exact
    // rather than a convention-dependent conversion.
    //
    // One member per edge type, each mapping to exactly one diagram; the
    // geometry does NOT auto-select. All four share the anchor zeta = 2.85 at
    // zero length and zero radius, so each degrades continuously to sharp.
    Idelchik1966Sharp,      // diagram 4-17, l/Dh <= 0.015
    Idelchik1966Thick,      // diagram 4-18a, deep hole, l/Dh > 0.015
    Idelchik1966Beveled,    // diagram 4-18b
    Idelchik1966Rounded,    // diagram 4-18c

    // Lichtarowicz, Duggins and Markland (1965): a LONG orifice, l/d 2-10,
    // down to Re = 10. This is the low-Reynolds regime an effusion hole
    // actually runs in -- Andrews' own plate C data spans Re 432 to 8700 --
    // where McGreehan-Schotsch is floored at re_min = 1e4 and returns a
    // near-constant.
    Lichtarowicz1965,

    // Fixed Cd. The value belongs to the caller: a measured plate value, or a
    // literature constant such as Florschuetz's 0.79 for a jet plate. Use
    // make_constant_discharge_correlation to set it.
    Constant
};

// -------------------------------------------------------------
// Cd computation - free functions (simple interface)
// -------------------------------------------------------------

// Sharp thin-plate orifice (ISO 5167-2, Reader-Harris/Gallagher), at the
// taps geom.tap. ISO's range is in metering_cd_in_range.
double Cd_sharp_thin_plate(const OrificeGeometry& geom, const OrificeState& state);

// The same Cd and its exact derivative in Re_D (the Re-correction, the
// A-dependent and the tap terms), held below ISO's lower Re_D bound as
// metering_cd_and_dRe describes. For the network orifice's Jacobian: the
// metering Cd moves with the flow through Re_D, a column the row used to miss.
std::tuple<double, double> Cd_sharp_thin_plate_and_dRe(const OrificeGeometry& geom,
                                                       const OrificeState& state);

// (Cd, dCd/dRe_D) of a NORMED metering correlation at the taps geom.tap: the
// one entry point the network element calls, so no metering formula lives in
// Python. ReaderHarrisGallagher only; Constant and UserFunction
// carry no formula of their own and are refused (std::invalid_argument).
//
// ReaderHarrisGallagher: below its lower Re_D bound Re_b (metering_cd_re_min)
// Eq. (4) is not evaluated: its Re terms grow without bound as Re_D falls
// (Cd 1.16 at Re_D 100 on a 4 mm pipe), while a real sharp orifice loses Cd
// in laminar flow. Cd is HELD instead (#873), entered C1 in
// x = ln Re_D over orifice::metering_hold::blend_ln_re = w just below Re_b:
//   Re_D >= Re_b:          the correlation itself (bit-identical to 0.7),
//   Re_b e^-w < Re_D < Re_b: C = C_b + s w (t^2 - 1) / 2,  dC/dRe_D = s t / Re_D,
//   Re_D <= Re_b e^-w:     C = C_b - s w / 2,           dC/dRe_D = 0,
// with t = 1 + ln(Re_D / Re_b) / w, C_b and s = Re_b dC/dRe_D the
// correlation's value and log-slope at Re_b. The slope ramps linearly in
// ln Re_D from s to 0, so Cd is monotone and Cd, dCd/dRe_D are continuous.
// The hold is an out-of-range placeholder; metering_cd_in_range says so.
// The 1.5 cap stays for flange taps in pipes below D ~32 mm.
//
// Stolz (ISO 5167:1980) was removed in 0.9.0 (#876): ISO replaced it with
// Reader-Harris/Gallagher (ISO 5167-1:1991/Amd 1:1998, ISO 5167-2:2003).
std::tuple<double, double> metering_cd_and_dRe(MeteringCdCorrelation id,
                                               const OrificeGeometry& geom,
                                               double Re_D);

// The lower Re_D bound Re_b of Reader-Harris/Gallagher's stated range for
// this geometry: the Re_D below which metering_cd_and_dRe holds Cd (see
// there). ISO 5167-2:2003 5.3.1:
//   corner / D-D/2 taps: 5000 (beta <= 0.56), 16000 beta^2 (beta > 0.56);
//   flange taps:         max(5000, 170 beta^2 D[mm]).
// Refuses (std::invalid_argument) an invalid geometry and every other
// correlation: Constant / UserFunction, whose Cd is the caller's.
double metering_cd_re_min(MeteringCdCorrelation id, const OrificeGeometry& geom);

// Is the plate inside the correlation's stated range? A diagnostic flag,
// never a selector or a refusal.
//   ReaderHarrisGallagher, ISO 5167-2:2003 5.3.1: d >= 12.5 mm,
//     50 <= D <= 1000 mm, 0.1 <= beta <= 0.75, and
//     Re_D >= metering_cd_re_min.
//   Constant / UserFunction: the Cd is the caller's, so this returns true.
bool metering_cd_in_range(MeteringCdCorrelation id, const OrificeGeometry& geom,
                          double Re_D);

// ISO 5167-2:2003 5.4.1 permanent pressure loss of an orifice plate between
// pipes, as a fraction of the tapping differential, and its derivative in C:
//   k = dw / dp = [s - C beta^2] / [s + C beta^2],  s = sqrt(1 - beta^4 (1 - C^2)).
// dw is the loss between the undisturbed pipe upstream and the fully
// recovered pipe downstream (equal diameters). The relation is Bernoulli to
// the vena contracta and a momentum balance after it, so its dp is the
// undisturbed-to-vena-contracta differential: D-D/2 taps measure that, and
// with corner taps' C the loss differs by up to 6% at beta 0.75 (#645).
// Incompressible, like the derivation. Refuses beta outside (0, 1) and C <= 0
// (std::invalid_argument).
std::tuple<double, double> permanent_loss_ratio_and_dC(double beta, double C);

// ISO 5167-1 velocity-of-approach factor E = 1 / sqrt(1 - beta^4). Returns 1
// outside 0 < beta < 1, the network's convention for "no correction": no
// single approach pipe (beta = 0), or a bore at least as wide as the pipe.
// A NaN or inf beta throws std::invalid_argument (#769).
double velocity_of_approach_factor(double beta);


// -------------------------------------------------------------
// Individual correlation implementations
// -------------------------------------------------------------

namespace orifice {

// Fallback Cd values used when a Constant correlation is selected without a
// value. Both are placeholders for a number the caller should supply, not
// correlations, and neither varies with anything.
namespace defaults {
// Sharp thin-plate metering orifice, high Re: the handbook round number.
constexpr double metering_cd = 0.61;
// Plain sharp-edged discharge hole: the round number the discharge-hole
// sources work around. Deliberately NOT spelled as
// mcgreehan_schotsch::cd_reference, which happens to share the value but
// means something else (the reference Cd inside Eq. (16)).
constexpr double discharge_cd = 0.60;
} // namespace defaults

// Tap positions of ISO 5167-2:2003, 5.3.2.1, as L1 = l1/D and L2' = l2'/D.
namespace taps {
constexpr double inch_mm  = 25.4;  // flange taps sit 1 inch from the faces [mm]
constexpr double dd2_L1   = 1.0;   // D and D/2 taps, upstream
constexpr double dd2_L2   = 0.47;  // D and D/2 taps, downstream (ISO's L2')
} // namespace taps

// Reader-Harris/Gallagher (1998) - ISO 5167-2:2003 Eq. (4)
namespace reader_harris {
constexpr double C0        = 0.5961;   // base coefficient
constexpr double C1        = 0.0261;   // beta^2 term
constexpr double C2        = 0.216;    // beta^8 term
constexpr double C3        = 0.000521; // Re-correction coefficient
constexpr double C3_exp    = 0.7;      // Re-correction exponent
constexpr double C4a       = 0.0188;   // second Re term, base
constexpr double C4b       = 0.0063;   // second Re term, A coefficient
constexpr double C4_beta   = 3.5;      // second Re term, beta exponent
constexpr double C4_re     = 0.3;      // second Re term, Re exponent
constexpr double C5a       = 0.043;    // upstream tap term base
constexpr double C5b       = 0.080;    // upstream tap term exp1 coefficient
constexpr double C5c       = 0.123;    // upstream tap term exp2 coefficient
constexpr double C5_exp1   = 10.0;     // exp(-10*L1) coefficient
constexpr double C5_exp2   = 7.0;      // exp(-7*L1) coefficient
constexpr double C5_A_coef = 0.11;     // (1 - 0.11*A) factor
constexpr double C6        = 0.031;    // downstream tap coefficient
constexpr double C6_M_coef = 0.8;      // (M2' - 0.8 M2'^1.1)
constexpr double C6_M_exp  = 1.1;      // M2^1.1 exponent
constexpr double C6_beta   = 1.3;      // beta^1.3 exponent
constexpr double A_coef    = 19000.0;  // A = (19000 beta / Re_D)^0.8
constexpr double A_exp     = 0.8;
constexpr double re_scale  = 1.0e6;    // 10^6 * beta / Re scaling
// Pipes below 71.12 mm (2.8 in) add 0.011 (0.75 - beta) (2.8 - D/25.4).
constexpr double small_pipe_D_mm  = 71.12;
constexpr double small_pipe_coef  = 0.011;
constexpr double small_pipe_beta  = 0.75;
constexpr double small_pipe_D_in  = 2.8;
// Range, ISO 5167-2:2003 5.3.1
constexpr double range_d_min_mm   = 12.5;
constexpr double range_D_min_mm   = 50.0;
constexpr double range_D_max_mm   = 1000.0;
constexpr double range_beta_min   = 0.1;
constexpr double range_beta_max   = 0.75;
constexpr double range_re_min     = 5000.0;
constexpr double range_beta_split = 0.56;     // corner, D-D/2: above, 16000 beta^2
constexpr double range_re_beta2   = 16000.0;
constexpr double range_re_flange  = 170.0;    // flange: Re_D >= 170 beta^2 D[mm]
} // namespace reader_harris

// The C1 hold of the RHG Cd below its lower Re_D bound (#873); the form
// is in metering_cd_and_dRe. blend_ln_re is the width w in ln Re_D over
// which the slope ramps to zero (Re_b e^-0.1 = 0.905 Re_b to Re_b). The held
// value sits s w / 2 above C_b, so w is the trade between that offset (~w)
// and the blend's curvature (~1/w). It is fixed ONCE by the offset: at most
// half of ISO 5167-2:2003 5.3.3.1's stated uncertainty of C at the bound,
// over every in-range RHG geometry. The worst is flange taps, D 50 mm, beta
// 0.75: |s|/C_b = 0.070 against 0.75%, so w <= 0.107; 0.1 is the round value
// under it (offset 0.35% there, <= 0.17% at beta <= 0.6). Never tuned to a
// score.
namespace metering_hold {
constexpr double blend_ln_re = 0.1;
} // namespace metering_hold

// Thickness correction (Idelchik model)
namespace thickness {
constexpr double reattach_coef    = 0.35;   // reattachment benefit coefficient
constexpr double reattach_exp     = 8.0;    // exp(-8*t_d) decay rate
constexpr double blasius_coef     = 0.316;  // Blasius f = 0.316/Re^0.25
constexpr double blasius_exp      = 0.25;
constexpr double friction_coef    = 5.67;   // friction loss calibration factor
constexpr double correction_min   = 0.5;
constexpr double correction_max   = 1.3;
constexpr double thin_plate_limit = 0.02;   // t/d <= this -> no correction
constexpr double re_floor         = 100.0;
} // namespace thickness

// Rounded-entry (Idelchik)
namespace rounded {
constexpr double K_well_rounded   = 0.04;   // K for r/d >= 0.15
constexpr double K_sharp          = 0.5;    // K for r/d = 0 (sharp)
constexpr double radius_threshold = 0.15;   // r/d >= this -> well-rounded
constexpr double re_correction_coef = 0.1;
constexpr double re_correction_ref  = 1.0e5;
constexpr double re_correction_exp  = 0.2;
} // namespace rounded

// -------------------------------------------------------------
// Idelchik (1966) - orifice in a large wall, Section IV
// -------------------------------------------------------------
//
// Idelchik, I.E. "Handbook of Hydraulic Resistance", 1st English edition,
// AEC-tr-6630 (1966). sources/junction/Idelchik.pdf (gitignored, copyrighted) --
// the same copy the junction work digitised diagrams 7-1..7-7 from.
//
// GEOMETRY: F1 = F2 = infinity, a hole in a wall between two large volumes.
// That is the effusion-plate case and it is why these correlations take
// DischargeHoleGeometry, which has no pipe diameter.
//
// zeta = DH / (rho w0^2 / 2), referenced to the HOLE velocity w0, and DH is
// the full permanent loss because nothing recovers downstream. Hence
//   Cd = 1 / sqrt(zeta)
// exactly. This is NOT true of Idelchik's in-a-pipe diagrams (4-13..4-16),
// whose zeta is referenced to the pipe velocity w1 and whose DH is the
// permanent loss rather than ISO 5167's tapping differential.
//
// Every value below was read off the page rendered at 400 dpi with
// `pdftoppm -r 400`, not the PDF text layer, which garbles these tables.
// See validation/cooling/extractions/idelchik_1966_wall_orifice.md.
namespace idelchik {

// Sharp-edged hole, Re >= 1e5 (diagram 4-17). Also the l/Dh = 0 and r/Dh = 0
// anchor of all three edge tables below, which is what makes them continuous
// with the sharp case.
constexpr double zeta_sharp = 2.85;

// Diagram 4-17 is the sharp (thin-walled) hole up to this depth; deeper
// holes are diagram 4-18a's. Only used to flag the range, never to select.
constexpr double sharp_l_over_d_max = 0.015;

// Above this Reynolds number zeta is Re-independent (diagram 4-17 item 1).
constexpr double re_fully_turbulent = 1.0e5;

// Diagram 4-17, low-Re branch: zeta = zeta_phi0(Re) + eps_re(Re).
// Re spans 25 to 1e6 -- three decades below McGreehan-Schotsch's re_min.
constexpr int re_n = 14;
constexpr double re_points[re_n] = {
    2.5e1, 4.0e1, 6.0e1, 1.0e2, 2.0e2, 4.0e2, 1.0e3,
    2.0e3, 4.0e3, 1.0e4, 2.0e4, 1.0e5, 2.0e5, 1.0e6};
constexpr double zeta_phi0[re_n] = {
    1.94, 1.38, 1.14, 0.89, 0.69, 0.54, 0.39,
    0.30, 0.22, 0.15, 0.11, 0.04, 0.01, 0.00};
constexpr double eps_re[re_n] = {
    1.00, 1.05, 1.09, 1.15, 1.23, 1.37, 1.56,
    1.71, 1.88, 2.17, 2.38, 2.56, 2.72, 2.85};

// Diagram 4-18a, thick-walled (deep) hole: zeta = zeta_thick(l/Dh) + lam*l/Dh.
constexpr int thick_n = 12;
constexpr double thick_l_over_d[thick_n] = {
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 4.0};
constexpr double thick_zeta[thick_n] = {
    2.85, 2.72, 2.60, 2.34, 1.95, 1.76, 1.67, 1.62, 1.60, 1.58, 1.55, 1.55};

// Diagram 4-18a low-Re: zeta = zeta_phi0 + k eps_re zeta' + lam l/Dh.
//
// The source PRINTS k = 0.342. That value is 1/2.85 rounded to three figures
// (1/2.85 = 0.350877), and the rounding is the entire 2.5% by which the
// source's own two formulas disagree at high Re: item 2 must reduce to item 1
// as eps_re -> 2.85 and zeta_phi0 -> 0, which forces k = 1/zeta_sharp exactly.
//
// We use the exact value. Using the printed one would leave a 2.5% STEP in
// Cd at the top of the table -- a C0 break in a solver input, and the same
// defect class as the 11% jump this whole correlation replaced. The printed
// figure is kept below so the discrepancy is on the record rather than
// silently corrected.
constexpr double thick_low_re_coef = 1.0 / zeta_sharp;
constexpr double thick_low_re_coef_as_printed = 0.342;

// Diagram 4-18b, beveled edges (bevel angle 40-60 deg per diagram 4-15).
constexpr int beveled_n = 12;
constexpr double beveled_l_over_d[beveled_n] = {
    0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20};
constexpr double beveled_zeta[beveled_n] = {
    2.85, 2.80, 2.70, 2.60, 2.50, 2.41, 2.33, 2.18, 2.08, 1.98, 1.84, 1.80};

// Diagram 4-18c, rounded edges. The leading (0, 2.85) is read from graph c,
// which the tabulated row starts one point after; it is also forced by the
// sharp-edged value, so the curve is continuous at r = 0.
constexpr int rounded_n = 10;
constexpr double rounded_r_over_d[rounded_n] = {
    0.0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.12, 0.16, 0.20};
constexpr double rounded_zeta[rounded_n] = {
    2.85, 2.72, 2.56, 2.40, 2.27, 2.06, 1.88, 1.60, 1.38, 1.37};

// Wall roughness used for the thick-hole friction term lam, which Idelchik
// takes from diagrams 2-2..2-5 as a function of Re and Delta/Dh. Those
// diagrams ARE the Colebrook/Nikuradse family, so friction.h's Haaland
// explicit form reproduces them rather than substituting for them. A drilled
// or laser-cut cooling hole is hydraulically smooth at these Reynolds
// numbers; the term is worth ~2% of zeta at l/Dh = 2.
constexpr double default_roughness_over_d = 0.0;

// Orifice in a straight pipe, F1 = F2 (diagrams 4-13..4-16, Re >= 1e4-1e5).
// zeta1 = DH / (rho w1^2 / 2) is referenced to the PIPE velocity and is the
// permanent loss; n = F0/F1 = beta^2, m = 1 - n.
//   4-13/15/16 sharp, beveled, rounded:  zeta1 n^2 = [m + sqrt(zeta' m)]^2
//   4-14 thick:  zeta1 n^2 = (0.5 + tau sqrt(m)) m + m^2 + lam l/Dh
// At n -> 0 each is the wall diagram's zeta (4-17/4-18), so zeta' and tau are
// taken FROM those tables: zeta' = (sqrt(zeta_wall) - 1)^2 and tau =
// zeta'_wall - 1.5. That reproduces 4-14..4-16's own zeta tables within 3.7%
// for n <= 0.8. Their printed zeta'/tau columns do worse where they disagree
// with the zeta tables beside them: 4-16's zeta' at r/Dh 0.16 (0.06) misses
// its own row by 25%; the wall table's 0.031 is within 3.4%.
// "At lam = 0.02 the values of zeta are completely determined from the
// curves" (4-14): the lam its table is computed at. Geometry only here, so
// the in-pipe ratio carries no Re.
constexpr double pipe_thick_lambda = 0.02;
// 4-14's inlet term and the Borda-Carnot 1 that, with tau, make up the wall
// hole's zeta' at n = 0.
constexpr double pipe_thick_inlet = 0.5;
constexpr double pipe_thick_wall_offset = 1.5;

// Inclined hole, OPT-IN (DischargeHoleGeometry::inclination_deg < 90).
//
// Diagram 3-2 (PDF p. 99, printed p. 93; also Eq. 3-1 on printed p. 82,
// Weisbach's formula): a straight entrance into a conduit mounted flush with
// the wall at an angle delta between the conduit axis and the wall,
// Re = w0 Dh / nu > 1e4:
//   zeta_in = 0.5 + 0.3 cos(delta) + 0.2 cos^2(delta)
// tabulated 0.96, 0.91, 0.81, 0.70, 0.63, 0.56, 0.50 at delta 20..90 deg.
//
// THE COMPOSITION IS THE SOURCE'S, ITS TRANSFER TO AN INCLINED HOLE IS OURS.
// Eq. (4-20) (printed p. 120) builds a wall orifice of arbitrary edge shape
// and thickness as zeta' + 1 + tau + zeta_fr (F1 = F2 = infinity), with
// zeta' the inlet-edge coefficient "determined as zeta from diagrams 3-1 to
// 3-3, and 3-6" -- a range that contains 3-2 -- and, for the thick-walled
// orifice, zeta' = 0.5 (Eq. 4-27; pipe_thick_inlet above). Diagram 4-18a's
// tabulated zeta'(l/Dh) is therefore 0.5 + 1 + tau(l/Dh). Replacing the
// 0.5 by diagram 3-2's zeta_in adds
//   dzeta = zeta_in(delta) - 0.5 = 0.3 cos(delta) + 0.2 cos^2(delta).
// What Idelchik does NOT say: tau (diagram 4-11, the inlet jet's interaction
// with the wall thickness) was measured on normal holes and is kept at its
// normal-hole value here. That, and applying 4-20 to an inclined hole at
// all, is our composition.
constexpr double inclination_cos_coef = 0.3;
constexpr double inclination_cos2_coef = 0.2;
// Diagram 3-2's own Reynolds condition.
constexpr double inclination_re_min = 1.0e4;
// The entrance loss is separable from the exit loss only once the inlet's
// vena contracta has reattached inside the hole: diagram 3-2 is a conduit
// entrance. Diagram 4-18a's zeta' reaches its plateau (1.55, tau ~ 0.05) at
// l/Dh = 2.0 and is still falling steeply below 1.0, so 2.0 is the floor;
// shorter holes are FLAGGED out of range, not refused. l is the length along
// the hole axis (t / sin(delta)).
constexpr double inclination_l_over_d_min = 2.0;

} // namespace idelchik

// -------------------------------------------------------------
// Lichtarowicz, Duggins and Markland (1965) - long orifices
// -------------------------------------------------------------
//
// Lichtarowicz, A., Duggins, R.K. and Markland, E. (1965). "Discharge
// coefficients for incompressible non-cavitating flow through long
// orifices." J. Mech. Engng Sci. 7(2), 210-219.
// sources/orifices/lichtarowicz-et-al-1965-...pdf (gitignored, copyrighted).
//
// WHY IT EARNS A PLACE. It plugs a coverage hole the other two leave:
//
//   correlation          Re range        l/d range
//   McGreehan-Schotsch   >= 1e4          long holes, with crossflow
//   Idelchik 4-18a       25 to 1e6       l/Dh up to 4
//   Lichtarowicz         10 to 2e4       2 to 10        <- this one
//
// A cooling hole sits at low Re: Andrews' effusion plate C spans Re 432 to
// 8718 across its own measured range, where Cd varies by 20%. McGreehan's
// floor makes it blind to that; measured against Lichtarowicz it reads
// +21.7% high at Re = 432.
//
// Equations read off the page rendered at 600 dpi, not the OCR text layer.
namespace lichtarowicz {

// Eq. (7): the ultimate (high-Re) discharge coefficient. Stated to 1.5%.
constexpr double cdu_c0 = 0.827;
constexpr double cdu_c1 = 0.0085;

// For 1.5 <= l/d < 2 the source gives a flat value instead, same accuracy.
constexpr double cdu_short = 0.810;

// Below l/d = 1.5 the source's design recommendation (1) is to AVOID the
// geometry entirely: "the discharge coefficient varies rapidly with l/d
// below this value, and there is the possibility of hysteresis in
// operation." So this is a refusal boundary, not a clamp -- see
// LichtarowiczCorrelation::Cd.
constexpr double l_over_d_min = 1.5;
constexpr double l_over_d_split = 2.0;
constexpr double l_over_d_max = 10.0;

// Eq. (12), the low-Re form:
//   1/Cd = 1/Cdu + (20/Re)(1 + 2.25 l/d)
//          - (0.005 l/d) / (1 + 7.5 (log10(0.00015 Re))^2)
// "fits all but a few points to better than 0.02 in the range of l/d from
// 2 to 10 and of Re from 10 to 2 x 10^4".
constexpr double visc_c0   = 20.0;
constexpr double visc_c1   = 2.25;
constexpr double trans_c0  = 0.005;
constexpr double trans_c1  = 7.5;
constexpr double trans_c2  = 0.00015;

// Eq. (12)'s fitted Re range. Its upper end is NOT the envelope's: design
// recommendation (2) (p. 218) takes Cd "substantially constant for values of
// Re greater than 2 x 10^4", given by Eq. (7), the "ultimate" Cd measured at
// Re = 2e4; Eq. (12) approaches it from below (|Cd12 - Cd7| 1.2e-3 at l/d 2
// to 2.8e-3 at l/d 10 at Re 2e4, about 2e-4 by 1e5). The data the paper
// collects stop at about 1e5 (summary), so re_data_max is the envelope's
// upper edge: discharge_cd_envelope flags Re below re_validated_min or above
// re_data_max.
constexpr double re_validated_min = 10.0;
constexpr double re_validated_max = 2.0e4;
constexpr double re_data_max = 1.0e5;

// Re is held at this floor rather than allowed to reach zero, where the
// 20/Re term diverges and Cd collapses to zero. The source plots Eq. (12)
// down to Re ~ 1 in its Fig. 10, so the floor is below the drawn curve and
// well below the validated range; it exists to keep the solver finite, not
// to express physics.
constexpr double re_floor = 1.0;

} // namespace lichtarowicz

// McGreehan and Schotsch (1988) - composite Cd for a long orifice with
// corner radiusing and inlet crossflow. ASME J. Turbomachinery 110(2),
// 213-217. Equation numbers below are the paper's own.
//
// This is a DIFFERENT configuration from the ISO 5167 family above: a hole
// discharging plenum-to-plenum (beta = 0), not a metering orifice in a pipe
// run. It is the correlation for a gas-turbine jet plate or cooling transfer
// hole. See validation/cooling/extractions/orifice_discharge_coefficient.md
// for the extraction and the checks behind every constant here.
namespace mcgreehan_schotsch {
// Eq. (8), orifice Reynolds baseline: Cd = re_c0 + re_c1/Re
constexpr double re_c0        = 0.5885;
constexpr double re_c1        = 372.0;
// Eq. (9), nozzle Reynolds baseline: Cd = nozzle_c0 - nozzle_c1/sqrt(Re)
constexpr double nozzle_c0    = 0.9981;
constexpr double nozzle_c1    = 4.73;
// Eq. (12), corner-radius effects function f
constexpr double corner_floor = 0.008;
constexpr double corner_coef  = 0.992;
constexpr double corner_exp1  = 5.5;
constexpr double corner_exp2  = 3.5;
// Eq. (14), L/d effects function g
constexpr double length_coef  = 1.3;
constexpr double length_decay = 1.606;
constexpr double length_a     = 0.435;
constexpr double length_b     = 0.021;
// Eq. (17), relative tangential (crossflow) velocity terms
constexpr double c1_exp       = 1.2;
constexpr double c2_coef      = 0.5;
constexpr double c2_exp       = 0.6;
constexpr double c2_cd_exp    = -0.5;
constexpr double c3_coef      = 0.5;
constexpr double c3_exp       = 0.9;
constexpr double rv_cd_exp    = -3.0;
// The paper's reference point: a sharp-edged Cd of 0.60 at Re = 3.2e4, which
// Eq. (8) reproduces to 0.02%. Every /0.6 divisor in Eqs. (1) and (17) is this
// number, so it is defined once here rather than repeated as a literal.
constexpr double cd_reference = 0.6;

// Crossflow regularisation width, in U1/Vi units. See cd_and_derivatives.
//
// Eq. (17) carries Rv^0.6 and Rv^0.9, both with unbounded slope at Rv = 0 --
// and U1/Vi = 0 is the DEFAULT, and the physically correct value for a
// plenum-fed jet plate. The singularity therefore sits exactly where a Newton
// solver spends most of its time, not in a remote corner of the envelope.
//
// The treatment regularises the INPUT, U1/Vi -> sqrt((U1/Vi)^2 + eps^2),
// rather than Rv: U1/Vi is the physical quantity eps should be scaled
// against, and Rv carries a cd_base-dependent factor that would make eps mean
// something different for every geometry.
//
// CHOSEN BY MEASUREMENT, between two walls bisected from opposite directions,
// the same way x_smooth_eps was. Both walls are sourced rather than picked:
//
//   LOWER WALL, eps >= 1.86e-5. The worst |dCd/d(U1/Vi)| anywhere must stay
//   within 10x the PHYSICAL derivative scale, which is 0.663 -- measured on
//   the exact correlation over u >= 0.01, the range Figs. 4-6 actually carry
//   data for. Below the wall the Jacobian entry is stiff enough to dominate
//   a Newton step for no physical reason.
//
//   UPPER WALL, eps <= 4.15e-4. The departure from Eq. (17) exactly must stay
//   within 10% of the scatter the correlation itself sits in. The paper gives
//   NO error statistic anywhere (extraction item 25); the +/-0.02 in Cd comes
//   from digitising the spread of its own Fig. 4 data.
//
//   The cost MUST be measured including u = 0. That is the default, and the
//   physically right value for a plenum-fed jet plate, and it is where the
//   regularisation bites hardest. An earlier sweep that started at u = 0.01
//   reported the cost as 3.3e-7 when it is really 8.6e-4 -- three orders out,
//   and it moved the upper wall by a factor of 21.
//
// | eps   | max |dCd/du| | x physical | cost     | % of bound | verdict |
// |-------|--------------|------------|----------|------------|---------|
// | 0     | unbounded    | --         | 0        | 0%         | stiff   |
// | 1e-6  | 21.5         | 32.4       | 5.5e-5   | 3%         | stiff   |
// | 2e-5  | 6.45         | 9.7        | 3.3e-4   | 16%        | no margin |
// | 3e-5  | 5.47         | 8.2        | 4.2e-4   | 21%        | ok      |  <-- default
// | 1e-4  | 3.35         | 5.1        | 8.6e-4   | 43%        | ok      |
// | 1e-3  | 1.26         | 1.9        | 3.4e-3   | 168%       | infidel |
//
// The default sits 15% across the window, deliberately towards the LOW-
// smoothing end -- below where the Y blend's default sits in its own window
// (29%). Operational experience on this solver is that OVER-smoothing is the
// worse failure: a residual that no longer matches the physics stalls the
// solve ("ghost residuals"), and that bit harder than a stiff Jacobian ever
// did. 3e-5 keeps a factor of 5 of headroom against the fidelity wall while
// still holding a margin against the stiffness one (8.2x of a 10x limit);
// 2e-5 would smooth marginally less but sits at 9.7x, with no room for a
// re-measurement to move it.
//
// There is NO PLATEAU here, and looking for one is how this was first got
// wrong: max |dCd/du| follows a clean eps^-0.4 power law, exactly as the
// Rv^-0.4 singularity predicts, so every eps trades derivative against
// fidelity and the choice has to come from the two walls. An early linear
// sweep appeared to show a plateau at 1.451 -- it had simply never sampled
// u below 0.002, which is where the whole singularity lives.
//
// Honest about what it does: it caps a derivative that physically diverges.
// It removes a mathematical spike, not a physical one.
//
// Passing eps = 0 recovers Eq. (17) exactly, as for the Y smoothing above.
constexpr double rv_smooth_eps = 3.0e-5;
// Eq. (4), orifice adiabatic expansion factor
constexpr double y_orifice_coef = 0.41;
// Eq. (7), the Cd-dependent blend between orifice and nozzle expansion
constexpr double x_blend_cd0   = 0.82;   // below this, Y = Y_o
constexpr double x_blend_slope = 8.333;  // X reaches 1 at Cd = 0.94
// Saturation smoothing widths. The paper's Eq. (7) is a bare linear ramp with
// no clamp at either end, and a HARD clamp would put an exactly-zero
// derivative outside [0.82, 0.94] plus a discontinuous jump of ~0.73 in
// dY/dCd at both knees -- a Newton hazard, and Cd > 0.94 is a real design
// point (r/d >= 0.2 at t/d = 2). These widths saturate smoothly instead.
// Passing eps = 0 recovers the paper's exact hard clamp.
// Cost, measured: max |Y_soft - Y_hard| = 0.0042, i.e. 0.47% of a typical Y.
constexpr double x_smooth_eps  = 0.05;   // in X units
constexpr double s_smooth_eps  = 0.01;   // in pressure-ratio units

// Stated validity floor for Eqs. (8) and (9).
constexpr double re_min       = 1.0e4;
// The end of the data Eq. (14)'s length function was fitted to: Fig. 3 plots
// Lichtarowicz et al.'s long orifices on an L/d axis to 10. The paper states
// no limit, but g keeps falling beyond it (Cd 0.74 at L/d 10, 0.66 at 20,
// 0.49 at 40 for a sharp inlet at high Re), the steepest extrapolation in the
// hole-Cd family, so beyond it is reported (DischargeCdEnvelope::l_over_d).
constexpr double l_over_d_max = 10.0;
// Above this r/d an ASME nozzle Cd is reached and further radiusing buys
// nothing (p.214). Eq. (12) is evaluated as-is; this is the documented knee.
constexpr double r_over_d_nozzle_limit = 0.82;

// Eq. (8). Sharp-edged orifice, plenum-to-plenum. Stated valid Re >= re_min.
double reynolds_baseline(double Re);

// Eq. (9). Nozzle equivalent of Eq. (8). Stated valid Re >= re_min.
double nozzle_baseline(double Re);

// Eq. (12), the f factor. f(0) = 1 (sharp corner, no correction);
// f decays to corner_floor for a fully rounded inlet.
double corner_factor(double r_over_d);

// Eq. (14), the g factor. g(0) = 1 to within 0.05%, so Eq. (13) reduces to
// its input at zero length -- an identity the fitted constants satisfy rather
// than one imposed here.
double length_factor(double L_over_d);

// Eq. (11). Reynolds baseline corrected for inlet corner radius.
double cd_with_corner(double Re, double r_over_d);

// Eqs. (13), (15) and (16). Adds the long-orifice effect.
//
// When r/d > 0 the paper applies its combined-effects correction: a revised
// basic Cd from Eq. (15) with g taken at r/d, and an effective length
// (L/d)' = L/d - r/d from Eq. (16). Eq. (16) is PRINTED as "L/D - r/d"; the
// capital D is an erratum in the paper (D is pipe diameter, undefined for
// plenum-to-plenum flow). See item 17a of the extraction.
double cd_with_corner_and_length(double Re, double r_over_d, double L_over_d);

// Eq. (17). Adds the relative tangential velocity effect.
//
// U1_over_Vi is the ratio of INLET (approach, supply-side) tangential
// velocity to ideal through-flow velocity. It is NOT a discharge-side
// crossflow ratio: feeding it Florschuetz's Gc/Gj -- spent-air crossflow in
// an impingement channel, on the far face of the plate -- is wrong and
// produces a plausible number biased low. A plenum-fed jet plate has
// U1_over_Vi = 0. See decision D8 of the extraction.
//
// Vi must be built from STATIC inlet conditions, not from a total pressure
// that already includes the tangential velocity head (item 22).
//
// Not monotonic: Cd rises above its zero-crossflow value by up to ~5.7% near
// U1_over_Vi ~ 0.09 before falling away. That is the source's own Fig. 4 and
// its data, not an artifact -- see check F of the extraction.
double cd(double Re, double r_over_d, double L_over_d, double U1_over_Vi,
          double eps = rv_smooth_eps);

// Eq. (17) applied to a baseline supplied by the caller, rather than one
// computed from Eqs. (8)-(16).
//
// This is how the source itself uses Eq. (17) in its own validation: Figs. 5
// and 6 plot Rohde's data against curves anchored to "a set baseline point at
// U1/Vi = 0" taken from Rohde's measurements (0.64, 0.73, 0.88), not to the
// chain's prediction for the same geometry -- the chain runs 3-12% higher,
// which is the paper's own remark that Rohde's "basic values are lower".
//
// Use it when a plate's zero-crossflow Cd is KNOWN -- a measured value, or a
// literature one such as Florschuetz's per-configuration Table 1 -- and only
// the crossflow correction is wanted.
double cd_with_crossflow(double cd_base, double U1_over_Vi,
                         double eps = rv_smooth_eps);

// Cd and its derivatives with respect to the two solver unknowns it depends
// on: (Cd, dCd/dRe, dCd/d(U1/Vi)). r/d and L/d are geometry and stay
// constant, so they never enter the Jacobian.
//
// TWO NUMERICAL TREATMENTS, both solver aids and both stated:
//
//  1. The crossflow input is regularised by rv_smooth_eps (above), which
//     bounds the worst dCd/d(U1/Vi) at 5.47 -- 8.2x the physical derivative
//     scale -- where the exact chain is unbounded. This one DOES move the
//     value, worst 4.2e-4 in Cd at U1/Vi = 0, which is 2% of the +/-0.02
//     scatter the correlation's own data sits in.
//
//  2. Below re_min the VALUE stays exactly floored -- Eq. (8) is outside its
//     validity there and diverges -- but the DERIVATIVE is continued from
//     re_min rather than reported as the true zero. A Newton step that
//     wanders below the floor would otherwise see no Re sensitivity at all
//     and stall, and the exact derivative also jumps discontinuously at
//     re_min (a KINK the smoothness scan reports at Re ~ 9908).
//
//     Softening the VALUE was measured and rejected: a soft-max floor buys
//     almost no gradient for real Cd error (eps = 1e4 in Re units recovers
//     only 16% of the live slope while moving Cd by 7.1e-3), and it would
//     distort a region where the correlation does not apply. Continuing the
//     derivative alone changes no reported Cd, and makes dCd/dRe continuous
//     across re_min. The value this returns agrees with cd() to within a
//     few ULP (measured worst 3.9e-16): same expression, different operation
//     order through the dual, so it is not bit-identical.
std::tuple<double, double, double> cd_and_derivatives(double Re, double r_over_d,
                                                      double L_over_d,
                                                      double U1_over_Vi);

// -------------------------------------------------------------
// Adiabatic expansion factor Y, Eqs. (4)-(7)
// -------------------------------------------------------------
//
// The paper's mass flow is Eq. (2):
//     W = Cd * Y * A_a * sqrt(2 g_c P_t1/(R T_t1) (P_t1 - P_s2))
// so Y is what carries compressibility in an otherwise incompressible
// orifice equation.
//
// WHERE THIS BELONGS. Y is for the INCOMPRESSIBLE formulation only. combaero's
// regime='compressible' path solves the isentropic nozzle exactly
// (solver_interface.cpp -> combaero::nozzle_flow, with a choked branch), and
// Eq. (5) reproduces that solve to 0.008% -- verified against it directly.
// Applying Y there would correct for compressibility twice.

// Critical pressure ratio, (2/(g+1))^(g/(g-1)). Below it the isentropic form
// gives DECREASING flow, so Y saturates here rather than following it down.
double critical_pressure_ratio(double gamma);

// Eq. (4). Orifice form: Y_o = 1 - 0.41 (1 - S)/gamma, S = P_s2/P_t1.
double expansion_orifice(double S, double gamma);

// Eq. (5). Nozzle form, the exact isentropic expansion factor. Reproduces
// combaero's own nozzle_flow to 0.008% over 0.6 <= S <= 0.99.
double expansion_nozzle(double S, double gamma);

// Eq. (7). Blend weight X = 8.333 (Cd - 0.82), saturated smoothly into [0, 1].
// eps = 0 gives the paper's exact unsaturated ramp clamped hard; the default
// trades 0.47% in Y for a derivative that is continuous and non-zero. See
// x_smooth_eps.
double expansion_blend_weight(double cd, double eps = x_smooth_eps);

// Eq. (6). Y = (1 - X) Y_o + X Y_n.
//
// An orifice is not a nozzle, and Y_o and Y_n differ by up to 17% at
// S = 0.6; the Cd-dependent blend is what interpolates. A sharp-edged hole
// (Cd < 0.82) gets Y_o; one rounded enough to behave like a nozzle
// (Cd > 0.94) gets Y_n.
double expansion_factor(double cd, double S, double gamma,
                        double eps = x_smooth_eps);
} // namespace mcgreehan_schotsch

// (L1, L2') of a tap arrangement in a pipe of diameter D [m].
std::tuple<double, double> tap_positions(OrificeTap tap, double D);

// Reader-Harris/Gallagher (1998) - ISO 5167-2:2003 Eq. (4), with the
// small-pipe term for D < 71.12 mm. The standard correlation for sharp-edged
// orifices. Held below its lower Re_D bound (metering_cd_and_dRe, #873).
double Cd_ReaderHarrisGallagher(double beta, double Re_D, double D,
                                OrificeTap tap = OrificeTap::Flange);

// McGreehan and Schotsch (1988) - the full chain, Eqs. (8) through (17).
// Convenience wrapper over orifice::mcgreehan_schotsch::cd; see that
// namespace for the per-equation stages and for what U1_over_Vi means.
double Cd_McGreehanSchotsch(double Re, double r_over_d, double L_over_d,
                            double U1_over_Vi,
                            double eps = mcgreehan_schotsch::rv_smooth_eps);

// Loss coefficient K from Cd: K = (1/Cd^2 - 1) * (1 - beta^4)
double K_from_Cd(double Cd, double beta);

// Cd from loss coefficient K
double Cd_from_K(double K, double beta);

// ISO 5167-2:2003 5.3.2.2, Eq. (5): the empirical expansibility of an orifice
// plate, eps = 1 - (c0 + c4 beta^4 + c8 beta^8) [1 - (p2/p1)^(1/kappa)].
// Valid for p2/p1 >= min_pressure_ratio only; its stated uncertainty is
// uncertainty_pct dp / (kappa p1) percent (5.3.3.2).
namespace iso5167_expansibility {
constexpr double c0                 = 0.351;
constexpr double c4                 = 0.256;
constexpr double c8                 = 0.93;
constexpr double min_pressure_ratio = 0.75;
constexpr double uncertainty_pct    = 3.5;
} // namespace iso5167_expansibility

// Jobson (1955), Proc. IMechE 169, 767: the compressible contraction
// coefficient of a sharp-edged orifice from its incompressible one, through
// the force defect on the reservoir wall (Eq. 5), f = 1/Ci - 1/(2 Ci^2).
// The paper evaluates Ci from 0.5 (Borda mouthpiece) to 0.7 (Figs. 3 and 4);
// p. 772 "tentatively" suggests Ci = 0.7 as the theory's upper limit, as it
// fails when the orifice approaches a nozzle (f -> 0.5); the beta
// extension (#854) is checked against ISO 5167-2 Eq. (5) up to its beta 0.75.
// These bounds set a diagnostic flag, never a refusal.
namespace jobson1955 {
constexpr double ci_min   = 0.5;   // Borda mouthpiece, f = 0
constexpr double ci_max   = 0.7;   // the paper's suggested upper limit (p. 772)
constexpr double beta_max = 0.75;  // the beta extension's checked range
} // namespace jobson1955

} // namespace orifice

// -------------------------------------------------------------
// Orifice correlation class (for polymorphic use)
// -------------------------------------------------------------

class OrificeCorrelationBase {
public:
    virtual ~OrificeCorrelationBase() = default;
    virtual double Cd(const OrificeGeometry& geom, const OrificeState& state) const = 0;
    virtual std::string name() const = 0;
};

// Factory function to create correlation objects
// Returns nullptr for UserFunction (use make_user_correlation instead)
std::unique_ptr<OrificeCorrelationBase> make_correlation(MeteringCdCorrelation id);

// Fixed-Cd correlation from an explicit value.
//
// make_correlation(MeteringCdCorrelation::Constant) can only give you the default,
// so this is the way to pin Cd to a chosen number -- a measured plate value,
// or a literature constant such as Florschuetz's 0.79 for a jet plate -- and
// to switch deliberately between a fixed Cd and a computed one.
std::unique_ptr<OrificeCorrelationBase> make_constant_correlation(double Cd);

// User-defined correlation from function
//
// IMPORTANT: User-provided functions should be smooth (C1 continuous) for
// best solver performance. Discontinuous derivatives can cause convergence issues.
//
// LIMITATIONS:
// - Network solver Jacobians assume constant Cd
// - If Cd depends on flow variables (Re, PR), Jacobians are incomplete
// - For production use with variable Cd, analytical Jacobians needed (future work)
using CdFunction = std::function<double(const OrificeGeometry&, const OrificeState&)>;
std::unique_ptr<OrificeCorrelationBase> make_user_correlation(
    CdFunction fn,
    const std::string& name = "UserDefined");

// User-defined correlation from tabulated data (beta, Re_D, Cd)
// Interpolates bilinearly in beta and log(Re_D)
//
// IMPORTANT LIMITATIONS (current implementation):
// 1. Bilinear interpolation is C0 continuous but NOT C1 continuous
//    - Derivative discontinuities at grid points can cause solver issues
//    - For production use, consider using smooth data or dense grids
// 2. Network solver Jacobians assume constant Cd
//    - When Cd = f(Re), Jacobians are incomplete (missing dCd/dRe terms)
//    - May reduce convergence rate in network solvers
// 3. Pressure-ratio dependent Cd creates circular dependency
//    - If Cd = f(PR), iterative solution needed (see solve_orifice_mdot)
//    - Jacobians not available for PR-dependent correlations
//
// BEST PRACTICES:
// - Use sufficiently dense grids (Delta_beta < 0.05, Delta_logRe < 0.2)
// - Ensure smooth underlying data (avoid measurement noise)
// - For critical applications, validate convergence with grid refinement
//
// FUTURE WORK (for production readiness):
// - Implement C1-continuous interpolation (bicubic splines)
// - Add analytical Jacobians with chain rule: d_mdot_ddP including dCd/dRe terms
// - Support for Cd = f(PR) with consistent Jacobians
std::unique_ptr<OrificeCorrelationBase> make_tabulated_correlation(
    const std::vector<double>& beta_values,
    const std::vector<double>& Re_values,
    const std::vector<std::vector<double>>& Cd_table,
    const std::string& name = "Tabulated");

// -------------------------------------------------------------
// Discharge-hole correlation class (for polymorphic use)
// -------------------------------------------------------------

// Geometry of one discharge hole. There is deliberately no pipe diameter:
// nothing here forms beta, and a hole in a wall has no pipe.
struct DischargeHoleGeometry {
    double d = 0.0;   // Hole diameter [m]
    double L = 0.0;   // Hole length ALONG ITS AXIS [m]. Equal to the wall
                      // thickness for a normal hole; t/sin(alpha) for a hole
                      // drilled at angle alpha to the wall.
    double r = 0.0;   // Inlet edge radius [m] (0 = sharp)
    double bevel = 0.0;  // Bevel depth along the axis [m] (0 = not beveled).
                         // Idelchik diagram 4-18b's argument is l/Dh, the
                         // beveled DEPTH over the hole diameter, at a bevel
                         // angle of 40-60 deg; it is not the wall thickness.
    double inclination_deg = 90.0;  // Angle between the hole axis and the
                                    // wall [deg], Idelchik's delta. 90 (the
                                    // default) = a normal hole and NO
                                    // inclination term; below 90 opts in to
                                    // diagram 3-2's inlet loss, composed with
                                    // Idelchik1966Thick only (the others
                                    // refuse it). See
                                    // orifice::idelchik::inclination_*.

    double L_over_d() const;   // Length ratio L/d [-]
    double bevel_over_d() const;  // Bevel depth ratio l/d [-]
    double r_over_d() const;   // Radius ratio r/d [-]
    double area() const;       // Hole area [m^2]
    bool is_inclined() const;  // inclination_deg < 90: the opt-in term applies

    bool is_valid() const;
};

// Idelchik diagram 3-2's inclined-entrance increment over the normal sharp
// entrance, dzeta = 0.3 cos(delta) + 0.2 cos^2(delta) [-] on the hole
// velocity head; exactly 0 at delta = 90 deg. Geometry only (no Re
// dependence). Throws std::invalid_argument outside 0 < delta <= 90.
double discharge_inclination_dzeta(double inclination_deg);

// Flow state at one discharge hole.
struct DischargeHoleState {
    double Re = 0.0;          // Hole Reynolds number, based on d [-]

    // Ratio of INLET (approach, supply-side) tangential velocity to ideal
    // through-flow velocity. A plenum-fed hole has 0. See the note on
    // orifice::mcgreehan_schotsch::cd for why a discharge-side crossflow
    // ratio must not be fed here.
    double U1_over_Vi = 0.0;
};

class DischargeCorrelationBase {
public:
    virtual ~DischargeCorrelationBase() = default;

    virtual double Cd(const DischargeHoleGeometry& hole,
                      const DischargeHoleState& flow) const = 0;

    // Solver-facing (f, J): (Cd, dCd/dRe, dCd/d(U1_over_Vi)). L/d and r/d are
    // geometry and never enter the Jacobian. Analytic, never finite
    // differences -- see the solver rule in CLAUDE.md.
    virtual std::tuple<double, double, double> Cd_and_derivatives(
        const DischargeHoleGeometry& hole,
        const DischargeHoleState& flow) const = 0;

    virtual std::string name() const = 0;
};

// Factory for discharge-hole correlations. Every member is implemented.
std::unique_ptr<DischargeCorrelationBase> make_discharge_correlation(
    DischargeCdCorrelation id);

// Which inputs of a discharge-hole correlation lie outside the range its
// source fitted or tabulated, by variable. INFORMATION to report beside the
// result, never a selector and never a warning: Cd is returned regardless,
// held or continued as each correlation documents, and the caller can do
// nothing better than know which variable is outside. The ranges are the
// sources' own:
//   re        Idelchik1966* below 25 or above 1e6 (table held at its ends);
//             McGreehanSchotsch1988 below re_min (value held); Lichtarowicz1965
//             below 10 (Eq. 12 continued to re_floor) or above re_data_max
//             (1e5, the end of the collected data; Eq. 7 held, flat). 2e4 to
//             1e5 is in range: Eq. 7, design recommendation (2); an inclined
//             Thick hole at or below
//             diagram 3-2's Re 1e4
//   l_over_d  Idelchik1966Sharp above 0.015; Idelchik1966Thick above 4
//             (zeta' held at its table end, only bore friction grows, a gentle
//             fall); McGreehanSchotsch1988 above l_over_d_max (10; Eq. 14's g
//             keeps falling); Lichtarowicz1965 outside 2 to 10 (held at 10);
//             an inclined Thick hole below inclination_l_over_d_min
//   edge      Idelchik1966Beveled / Rounded beyond their bevel or r/d table;
//             an inclined hole on any correlation but Idelchik1966Thick (its
//             Cd refuses the geometry)
// Constant is always inside.
struct DischargeCdEnvelope {
    bool re = false;
    bool l_over_d = false;
    bool edge = false;
    bool in_range() const { return !(re || l_over_d || edge); }
};

DischargeCdEnvelope discharge_cd_envelope(DischargeCdCorrelation id,
                                          const DischargeHoleGeometry& hole,
                                          const DischargeHoleState& flow);

// True when no input lies outside the source's range:
// discharge_cd_envelope(...).in_range(). A FLAG, never a selector.
bool discharge_cd_in_range(DischargeCdCorrelation id,
                           const DischargeHoleGeometry& hole,
                           const DischargeHoleState& flow);

// The hole Reynolds number below which the selected correlation HOLDS its Cd
// value (Re is clamped up to it), so the value's true dCd/dRe there is zero.
// Cd_and_derivatives does not always report that zero: McGreehanSchotsch1988
// and the Idelchik Re tables continue the slope at the floor, a Newton aid
// for a row that evaluates Cd at its own iterate. A caller that differentiates
// a quantity built on the held VALUE -- an implicit-function sensitivity of a
// self-consistent flow (#760, #801) -- must use the zero below this Re.
//   McGreehanSchotsch1988  re_min (1e4)
//   Idelchik1966*          25, the first table knot (Beveled and Rounded do
//                          not read Re at all)
//   Lichtarowicz1965       re_floor (its derivative is already 0 below it)
//   Constant               0 (no Re dependence)
double discharge_cd_re_floor(DischargeCdCorrelation id);

// Idelchik (1966) permanent loss of a discharge hole as an orifice in a
// straight pipe, zeta1 = DH / (rho w1^2 / 2), PIPE-velocity referenced
// (diagrams 4-13..4-16, see orifice::idelchik). The hole's geometry picks the
// diagram:
//   Idelchik1966Sharp                       4-13
//   Idelchik1966Thick, Lichtarowicz1965,
//   McGreehanSchotsch1988                   4-14 at L/d (M-S's edge radius
//                                           is not in the in-pipe shape)
//   Idelchik1966Beveled                     4-15 at bevel/d
//   Idelchik1966Rounded                     4-16 at r/d
// High-Re, incompressible, like the source. Throws std::invalid_argument
// outside 0 < beta < 1 and for Constant (no geometry).
double discharge_pipe_zeta(DischargeCdCorrelation id, const DischargeHoleGeometry& hole,
                           double beta);

// The factor R that carries a plenum-measured discharge Cd into a pipe of
// diameter ratio beta (#645): C_pipe = Cd_wall(Re, hole) R(beta). C_eq(beta)
// solves ISO 5167-2's permanent-loss relation k(C)(1 - beta^4)/(C^2 beta^4)
// = discharge_pipe_zeta, i.e. it is the vena-contracta C that gives the
// source's in-pipe loss through the network's recovery row; R = C_eq(beta) /
// C_eq(0), so R -> 1 as the pipe grows and the Cd stays the plenum one. The
// wall Cd keeps its Re and crossflow dependence; Idelchik gives the in-pipe
// beta shape at high Re only. Geometry only: no state, no Jacobian. Returns
// 1 outside 0 < beta < 1 and for Constant.
double discharge_pipe_cd_ratio(DischargeCdCorrelation id, const DischargeHoleGeometry& hole,
                               double beta);

// Fixed-Cd discharge correlation from an explicit value. This is the way to
// pin a measured or literature Cd; make_discharge_correlation(Constant) can
// only give you orifice::defaults::discharge_cd.
std::unique_ptr<DischargeCorrelationBase> make_constant_discharge_correlation(
    double Cd);

// -------------------------------------------------------------
// Orifice flow calculations (uses incompressible.h internally)
// -------------------------------------------------------------

// Mass flow through orifice given Cd and expansibility factor epsilon
// mdot = Cd * E * epsilon * A * sqrt(2 * rho * dP)
double orifice_mdot(const OrificeGeometry& geom, double Cd, double dP,
                    double rho, double epsilon = 1.0);

// Pressure drop for given mass flow
// dP = (mdot / (Cd * E * epsilon * A))^2 / (2 * rho)
double orifice_dP(const OrificeGeometry& geom, double Cd, double mdot,
                  double rho, double epsilon = 1.0);

// Solve for Cd given measured mdot and dP
// Cd = mdot / (A * sqrt(2 * rho * dP))
double orifice_Cd_from_measurement(const OrificeGeometry& geom,
                                    double mdot, double dP, double rho);

// -------------------------------------------------------------
// Iterative solver for Cd-Re coupling
// -------------------------------------------------------------

// Solve for mass flow rate accounting for Cd-Re_D dependency.
//
// This function iteratively solves the coupled system:
//   mdot = eps * Cd(Re_D) * A * sqrt(2 * rho * Delta_P)
//   Re_D = 4 * mdot / (pi * D * mu)
//
// The iteration continues until mdot converges within the specified tolerance.
//
// Algorithm:
//   1. Start with initial Cd guess (typically 0.61 for sharp orifices)
//   2. Calculate mdot using current Cd and expansibility factor eps
//   3. Update Re_D based on new mdot
//   4. Recalculate Cd using updated Re_D
//   5. Repeat until |mdot_new - mdot_old| < tol * mdot_old
//
// Parameters:
//   geom        : Orifice geometry
//   dP          : Differential pressure [Pa]
//   rho         : Upstream density [kg/m^3]
//   mu          : Dynamic viscosity [Pa*s]
//   P_upstream  : Absolute upstream pressure [Pa] (for compressibility)
//   kappa       : Isentropic exponent cp/cv [-] (0 = incompressible)
//   correlation : Cd correlation to use (default: ReaderHarrisGallagher)
//   tol         : Relative convergence tolerance (default: 1e-6)
//   max_iter    : Maximum iterations (default: 20)
//
// Returns: Converged mass flow rate [kg/s]
//
// Throws: std::runtime_error if iteration fails to converge
//
// Note: For incompressible flow, set kappa = 0 (expansibility eps = 1)
//       For typical applications, convergence occurs in 3-5 iterations
double solve_orifice_mdot(
    const OrificeGeometry& geom,
    double dP,
    double rho,
    double mu,
    double P_upstream = 101325.0,
    double kappa = 0.0,
    MeteringCdCorrelation correlation = MeteringCdCorrelation::ReaderHarrisGallagher,
    double tol = 1e-6,
    int max_iter = 20);

// -------------------------------------------------------------
// Orifice flow result bundle
// -------------------------------------------------------------

// Bundle of orifice flow properties for convenient access
// All properties computed from (geom, dP, T, P, mu, Z) in a single call
struct OrificeFlowResult {
    double mdot;           // Mass flow rate [kg/s]
    double v;              // Velocity through orifice throat [m/s]
    double Re_D;           // Pipe Reynolds number (based on D) [-]
    double Re_d;           // Orifice Reynolds number (based on d) [-]
    double Cd;             // Discharge coefficient [-]
    double epsilon;        // Expansibility factor [-] (1.0 for incompressible)
    double rho_corrected;  // Density corrected for compressibility [kg/m^3]
};

// Compute all orifice flow properties at once with real gas correction
// Parameters:
//   geom        : Orifice geometry
//   dP          : Differential pressure [Pa]
//   T           : Upstream static temperature [K]
//   P           : Absolute upstream static pressure [Pa]
//   mu          : Dynamic viscosity [Pa*s]
//   Z           : Compressibility factor [-] (default: 1.0 = ideal gas)
//   X           : Mole fractions [-] (optional, uses air if empty)
//   kappa       : Isentropic exponent cp/cv [-] (default: 0.0 = incompressible)
//   correlation : Cd correlation to use (default: ReaderHarrisGallagher)
// Returns: OrificeFlowResult struct with all properties
OrificeFlowResult orifice_flow(
    const OrificeGeometry& geom,
    double dP,
    double T,
    double P,
    double mu,
    double Z = 1.0,
    const std::vector<double>& X = {},
    double kappa = 0.0,
    MeteringCdCorrelation correlation = MeteringCdCorrelation::ReaderHarrisGallagher);

// -------------------------------------------------------------
// Utility functions
// -------------------------------------------------------------

// Velocity through orifice from mass flow rate with real gas correction
// v = mdot / (rho_corrected * A) where rho_corrected = rho / Z
double orifice_velocity_from_mdot(double mdot, double rho, double d, double Z = 1.0);

// Orifice area from beta ratio
// A = pi * (D * beta / 2)^2
double orifice_area_from_beta(double D, double beta);

// Beta ratio from diameters
// beta = d / D
double beta_from_diameters(double d, double D);

// Orifice Reynolds number from mass flow rate
// Re_d = 4 * mdot / (pi * d * mu)
double orifice_Re_d_from_mdot(double mdot, double d, double mu);

// -------------------------------------------------------------
// Compressible flow correction
// -------------------------------------------------------------

// Expansibility factor for compressible gas flow through orifices.
//
// The expansibility factor eps accounts for gas expansion as it accelerates
// through the orifice. For incompressible flow, eps = 1.0.
//
// Formula from ISO 5167-2:2003, Section 5.3.2.2:
//   eps = 1 - (0.351 + 0.256*beta^4 + 0.93*beta^8) * [1 - (1 - tau)^(1/kappa)]
//
// where:
//   tau = Delta_P / P_upstream  (pressure ratio)
//   beta = d/D           (diameter ratio)
//   kappa = cp/cv        (isentropic exponent)
//
// Usage in compressible mass flow:
//   mdot = eps * Cd * A * sqrt(2 * rho_upstream * Delta_P)
//
// Valid for:
//   - 0.1 <= beta <= 0.75
//   - tau <= 0.25 (Delta_P/P <= 25%)
//   - Ideal gas behavior
//
// Parameters:
//   beta        : Diameter ratio d/D [-]
//   dP          : Differential pressure [Pa]
//   P_upstream  : Absolute upstream static pressure [Pa]
//   kappa       : Isentropic exponent cp/cv [-]
//
// Returns: Expansibility factor eps [-] (dimensionless, 0 < eps <= 1)
//
// Reference: ISO 5167-2:2003, Measurement of fluid flow by means of
//            pressure differential devices inserted in circular cross-section
//            conduits running full
double expansibility_factor(double beta, double dP, double P_upstream, double kappa);

// -------------------------------------------------------------
// Jobson (1955): compressible contraction of a sharp-edged orifice (#854)
// -------------------------------------------------------------
//
// Jobson's control volume (his Fig. 2): the momentum of the jet balances the
// drive across the orifice area plus a force defect F = f mdot^2/(rho A) on
// the reservoir wall, whose coefficient f is fixed by incompressible flow and
// then carried unchanged to compressible flow. The jet is isentropic from the
// upstream stagnation state (polytropic index n); choked, the jet's first
// throat is sonic and the back pressure acts outside it (Eq. 15), so C keeps
// rising past choke as the jet widens. Mass flow m = Cv C K_N A sqrt(p0 rho0).
//
// APPROACH VELOCITY (#854, not in the paper). In a pipe, the approach flow
// enters the control volume at V1 = mdot beta^2 / (rho1 A) with static p1,
// rho1:  (p1 - p) A + f mdot^2/(rho1 A) = mdot (u - V1),  u from the TOTAL
// state. Its incompressible limit fixes
//   f(Ci, beta) = 1/Ci - 1/(2 Ci^2) - beta^2 + beta^4 / 2,
// Jobson's Eq. 5 at beta = 0. Choked, the same V1 term enters Eq. 15's
// balance, (p1 - p) A - (pc - p) C A + f mdot^2/(rho1 A) = mdot (uc - V1).
// With a velocity coefficient Cv (actual flow = Cv x theoretical) V1 is the
// ACTUAL approach velocity, Cv mdot beta^2 / (rho1 A), the one the pipe's
// static p1 carries, and the incompressible limit then fixes
//   f(Ci, beta, Cv) = 1/Ci - 1/(2 Ci^2) - Cv beta^2 + Cv^2 beta^4 / 2,
// the form above at Cv = 1.

// The force-defect coefficient f(Ci, beta, Cv) above [-].
double jobson_force_defect(double Ci, double beta = 0.0, double Cv = 1.0);

// The incompressible contraction Ci of an orifice whose incompressible flow
// is m = Cd E A sqrt(2 rho1 dp) (ISO's form, E = 1/sqrt(1 - beta^4)), with
// the velocity coefficient Cv as an overall factor (Jobson's Cd = Cv Ci):
// Cv Ci / sqrt(1 - (Cv Ci)^2 beta^4) = Cd E. beta = 0: Ci = Cd / Cv.
double jobson_ci_from_cd(double Cd, double beta = 0.0, double Cv = 1.0);

// The paper's C(r) for a reservoir-fed orifice (beta = 0): Eq. (11) for
// r = p/p0 >= r_c and Eq. (15) below it. Ci in (0, 1], n > 1, 0 < r <= 1.
double jobson_contraction_coefficient(double r, double Ci, double n);

// Expansibility the beta-extended theory implies, comparable with ISO
// 5167-2's eps: mdot / mdot_inc at the same upstream STATIC state (p1, rho1)
// and the same Ci, with mdot_inc = Ci A sqrt(2 rho1 dp) / sqrt(1 - Ci^2 b^4).
// r = p2/p1; the approach velocity's stagnation rise is solved to
// consistency (a bracketed Newton on the jet flux per bore area; the root is
// unique, #957). Subcritical and choked alike, any r in (0, 1]. Where the
// momentum quadratic has no real root (high beta Ci, outside jobson1955's
// range) the closure holds its discriminant at zero, as the network kernel
// does (in_range false there).
double jobson_expansibility(double r, double beta, double Ci, double n);

#endif // ORIFICE_H
