#pragma once

// Pin-fin array heat transfer and friction, from provenanced parameter sets
// (#335).
//
// WHAT A SET IS. A data fit with a validity box, carried in the source's own
// form and basis. No generic pin-fin law is assumed: every printed
// correlation found (Han, Dutta and Ekkad 2012, sec. 4.3; Armstrong and
// Winstanley 1988) is a single-rig fit or collapses one geometry variable at
// most. Outside its box a set still evaluates and flags `extrapolated`; it
// never refuses, because a user supplying their own coefficients knows their
// own hardware.
//
// THE CANONICAL BASIS. Every evaluator takes and returns the same basis, and
// converts a set's native basis exactly, from geometry alone:
//
//   Re_D = mdot D / (mu A_min)          pin diameter, velocity at A_min
//   Nu_D = h D / k                      h on the set's surface (PinNuSurface)
//   f    = dP / (2 rho Vmax^2 N)        per row, Vmax at A_min (Metzger)
//
// so dP = 2 rho Vmax^2 N f. The code #332 removed used rho Vmax^2 / 2, a
// factor of 4 low against this definition (Armstrong and Winstanley 1988,
// nomenclature).
//
// GEOMETRY. One pin per X x S unit cell in both arrangements (pin density
// 1 / (X S D^2)). S is the TRANSVERSE pitch, X the STREAMWISE pitch, H the
// pin height (= channel height), all over D. The minimum flow area is set by
// the smaller of the transverse gap (S - 1) and, for staggered arrays, the
// diagonal gap 2 (sqrt((S/2)^2 + X^2) - 1). Every source geometry here is
// transverse-limited; the diagonal case is carried so a user's dense array is
// not silently mis-scaled.
//
// BELOW THE FITTED Re RANGE the power law continues on the smoothed |Re| and
// the result is flagged. There is no handover baseline: a pin array has no
// laminar analogue of the smooth duct's Gnielinski limit, and inventing one
// would be unsourced.
//
// ROBUSTNESS, as for the rib sets: Re enters as sqrt(Re^2 + eps^2), so values
// are even in Re, derivatives odd and finite through Re = 0; evaluators never
// throw; parameters that can only be mistakes are rejected once by the
// validate_* functions; derivatives are analytic (solver (f, J) rule).
//
// NOT YET CARRIED (declared, so a caller does not assume them):
//   * row-count correction (Metzger et al. 1986 row curve, to be digitised):
//     a geometry with a different N than the set's own is flagged;
//   * channel convergence (Metzger 1986 phi; Brown 1980; Brigham 1984), whose
//     cause is unresolved in the literature;
//   * long pins (Faulkner 1971, an exponential geometry coefficient, a later
//     Nu form);
//   * inline friction from one rig only (Chyu 1990), at one geometry.

#include <string>

#include "rib_correlation.h"

namespace combaero {
namespace cooling {

// Smoothing scale for |Re|, as RATIO_RE_EPS: negligible above Re ~ 100.
constexpr double PIN_FIN_RE_EPS = 1.0;

// A two-segment friction set blends C1 in ln Re over
// [Re_split / PIN_FIN_F_BLEND, Re_split * PIN_FIN_F_BLEND]; each segment is
// exact outside that band.
constexpr double PIN_FIN_F_BLEND = 1.25;

enum class PinArrangement { Staggered, Inline };

// Which surface a set's Nu describes. Array element heat transfer needs Total
// (pin and endwall, area-weighted); Endwall and Pin sets exist because
// sources publish them and they show the split.
enum class PinNuSurface { Total, Endwall, Pin };

enum class PinReBasis {
  // Re_D = mdot D / (mu A_min). Canonical. Metzger, Chyu, Damerow, Lawson.
  DiameterVmax,
  // Re_D' = mdot D' / (mu A'), D' = 4 V / A_t, A' = V / L (VanFossen 1982).
  VanFossenDprime,
};

enum class PinFrictionBasis {
  // f = dP / (2 rho Vmax^2 N). Canonical (Metzger, via Armstrong and
  // Winstanley 1988).
  PerRowVmax,
  // f = dP_T rho / (2 (N - 1) (w / A_min)^2): per contraction-expansion pair,
  // total pressure first to last row (Damerow et al. 1972, Eq. 10).
  PerRowGapVmax,
  // f = 2 dP / (rho Vmax^2 N) = dP / (0.5 rho Vmax^2 N): per row on the
  // dynamic head, 4x the canonical value (Chyu 1990, Eq. 16).
  PerRowHalfHeadVmax,
};

struct PinFinGeometry {
  double S_D = 2.5;  // transverse pitch / D
  double X_D = 2.5;  // streamwise pitch / D
  double H_D = 1.0;  // pin height (channel height) / D
  int N_rows = 10;
  PinArrangement arrangement = PinArrangement::Staggered;
};

// ---- Exact geometry converters (Armstrong and Winstanley 1988, Eqs 8, 10,
// 13, re-derived from the unit cell). All are ratios to the pin diameter or
// to A_min. ----

// Minimum gap per transverse pitch, over D: min(S - 1, diagonal) staggered,
// S - 1 inline.
double pin_fin_min_gap_D(const PinFinGeometry &g);
// A_min / A_frontal = gap / S; equivalently Vchannel / Vmax.
double pin_fin_amin_over_afrontal(const PinFinGeometry &g);
// D' / D = H (4 X S - pi) / (2 X S + pi (H - 1/2)).
double pin_fin_dprime_over_D(const PinFinGeometry &g);
// A' / A_min = (S - pi / (4 X)) / gap.
double pin_fin_aprime_over_amin(const PinFinGeometry &g);
// D_h / D = 4 X H gap / (2 X S + pi (H - 1/2)) (tube-bank D_h = 4 A_min L / A_t).
double pin_fin_dh_over_D(const PinFinGeometry &g);

// Areas per wall per unit cell, over the cell's base (planform) area X S D^2,
// with each pin counted to mid-height (each wall feeds half the pin).
struct PinFinAreaFractions {
  double endwall_exposed = 0.0;  // (X S - pi/4) / (X S)
  double pin = 0.0;              // (pi H / 2) / (X S)
  double pin_over_total = 0.0;   // A_f / A_t, the fin-efficiency weight
};
PinFinAreaFractions pin_fin_area_fractions(const PinFinGeometry &g);

// Pin fin and array efficiency (Armstrong and Winstanley 1988, Eqs 14-16):
//   eta_fin = tanh(m L) / (m L),  m = sqrt(4 h / (k_pin D)),  L = H / 2
//   eta_t   = 1 - (A_f / A_t) (1 - eta_fin)
// D and H in metres, h in W/(m^2 K), k_pin in W/(m K). The derivative is
// analytic and finite at h -> 0 (series), where eta -> 1.
struct PinFinEfficiency {
  double eta_fin = 1.0;
  double eta_t = 1.0;
  double deta_fin_dh = 0.0;
  double deta_t_dh = 0.0;
};
PinFinEfficiency pin_fin_array_efficiency(double h, double k_pin, double D,
                                          double H, double A_f_over_A_t);

// ---- Sets ----

// Nu_native = C * (Re_native)^Re_exp * Pr^Pr_exp * (X/ref)^.. (S/ref)^..
// (H/ref)^.. in the set's native basis, returned in the canonical one.
struct PinFinNuSet {
  std::string name;
  std::string source;
  std::string validity_source;
  RibProvenance provenance = RibProvenance::User;
  PinArrangement arrangement = PinArrangement::Staggered;
  PinNuSurface surface = PinNuSurface::Total;
  PinReBasis re_basis = PinReBasis::DiameterVmax;

  double C = 0.0;
  double Re_exp = 0.0;
  double Pr_exp = 0.0;
  RibTerm term_XD, term_SD, term_HD;

  RibRange valid_Re;  // in the NATIVE basis, as the source states it
  RibRange valid_SD, valid_XD, valid_HD, valid_Nrows;
  double valid_Pr = 0.0;  // 0 means unstated
  StatedAccuracy accuracy_Nu;
};

// f_native = C_i * Re_D^Re_exp_i * (S/ref)^.. (X/ref)^.. (H/ref)^.. on segment
// i; segment 2 applies above Re_split (0 = one segment), blended C1 in ln Re.
// Re is always Re_D (every friction source in hand uses it).
struct PinFinFrictionSet {
  std::string name;
  std::string source;
  std::string validity_source;
  RibProvenance provenance = RibProvenance::User;
  PinArrangement arrangement = PinArrangement::Staggered;
  PinFrictionBasis basis = PinFrictionBasis::PerRowVmax;

  double C1 = 0.0;
  double Re_exp1 = 0.0;
  double C2 = 0.0;
  double Re_exp2 = 0.0;
  double Re_split = 0.0;
  RibTerm term_SD, term_XD, term_HD;

  RibRange valid_Re, valid_SD, valid_XD, valid_HD, valid_Nrows;
  StatedAccuracy accuracy_f;
};

// A data-backed ratio applied to another set's output, e.g. an arrangement
// transfer: Nu_to = Nu_from * C_Nu Re_D^Nu_Re_exp. Applying a ratio measured
// on one rig to a set from another is a transfer the caller chooses; the
// modifier records where the ratio came from.
struct PinFinRatioModifier {
  std::string name;
  std::string source;
  RibProvenance provenance = RibProvenance::User;
  PinArrangement from_arrangement = PinArrangement::Staggered;
  PinArrangement to_arrangement = PinArrangement::Staggered;
  PinNuSurface surface = PinNuSurface::Total;

  double C_Nu = 1.0;
  double Nu_Re_exp = 0.0;
  bool has_f = false;  // false: the source gives no friction ratio
  double C_f = 1.0;
  double f_Re_exp = 0.0;

  RibRange valid_Re, valid_SD, valid_XD, valid_HD;
};

struct PinFinNuResult {
  double Nu = 0.0;         // canonical Nu_D on the set's surface
  double dNu_dRe = 0.0;    // d Nu_D / d Re_D, analytic, odd in Re
  double dNu_dPr = 0.0;    // d Nu_D / d Pr at fixed Re_D, analytic (#1254)
  double Re_native = 0.0;  // the Re the set was evaluated at
  bool extrapolated = false;
};

struct PinFinFrictionResult {
  double f = 0.0;       // canonical, dP / (2 rho Vmax^2 N)
  double df_dRe = 0.0;  // d f / d Re_D
  bool extrapolated = false;
};

struct PinFinModifierResult {
  double ratio_Nu = 1.0;
  double dratio_Nu_dRe = 0.0;
  double ratio_f = 1.0;  // 1 when has_f is false
  double dratio_f_dRe = 0.0;
  bool has_f = false;
  bool extrapolated = false;
};

// Evaluators. Re_D is the canonical Reynolds number; Pr is guarded smoothly.
// None throws.
PinFinNuResult evaluate_pin_fin_nu(const PinFinNuSet &set,
                                   const PinFinGeometry &geom, double Re_D,
                                   double Pr);
PinFinFrictionResult evaluate_pin_fin_friction(const PinFinFrictionSet &set,
                                               const PinFinGeometry &geom,
                                               double Re_D);
PinFinModifierResult evaluate_pin_fin_modifier(const PinFinRatioModifier &mod,
                                               const PinFinGeometry &geom,
                                               double Re_D);

// Reject what can only be a mistake. Hard errors, named by field.
void validate_pin_fin_geometry(const PinFinGeometry &geom);
void validate_pin_fin_nu_set(const PinFinNuSet &set);
void validate_pin_fin_friction_set(const PinFinFrictionSet &set);
void validate_pin_fin_modifier(const PinFinRatioModifier &mod);

// ---- Shipped sets. Equations confirmed from page images; see
// validation/cooling/extractions/pin_fin_sources.md. ----

// Metzger, Shepard and Haley (1986), ASME 86-GT-132, via Armstrong and
// Winstanley (1988) Eq. 2: Nu_D = 0.135 Re_D^0.69 (X/D)^-0.34, pin and
// endwall combined, staggered. Fitted at H/D 1, S/D 2.5, 1.5 <= X/D <= 5,
// 10 rows; the validity box is the review's recommended short-pin limits
// (H/D <= 3, 2 <= S/D <= 4, 1.5 <= X/D <= 5, 1e3 <= Re_D <= 1e5), where it
// found +/-20% against other labs' data (stated). No Pr term: air data.
// The default heat-transfer set.
PinFinNuSet metzger_1986_staggered_nu();

// Metzger, Fan and Shepard (1982), Heat Transfer 1982 vol. 3, pp. 137-142,
// via Armstrong and Winstanley (1988) Eqs 20-21:
//   f = 0.317 Re_D^-0.132  (1e3 < Re_D < 1e4)
//   f = 1.76  Re_D^-0.318  (1e4 < Re_D < 1e5)
// fitted at H/D 1, S/D 2.5, 1.5 <= X/D <= 5 (+/-15%, stated); the review
// extends it to 0.5 <= H/D <= 6, 2 <= S/D <= 4 on Peng's data. The branches
// meet in value at 1e4 but not in slope; the C1 blend supplies the slope.
// The original paper is not accessible here, so fidelity is unchecked; Lawson
// et al. (2011) read about 25% lower at S/D 4. The default friction set.
PinFinFrictionSet metzger_1982_staggered_friction();

// VanFossen (1982), J. Eng. Power 104, 268 (NASA TM-81696), Eq. 16:
// Nu_D' = 0.153 Re_D'^0.685 on D' = 4V/A_t and A' = V/L, pin and endwall
// combined (h assumed equal on both), staggered equilateral arrays,
// H/D 0.5 and 2, 4 rows, 300 < Re_D' < 6e4. Metzger's H/D 1 data, converted
// to D', fall on it (Brigham and VanFossen 1983, Fig. 6).
PinFinNuSet vanfossen_1982_staggered_nu();

// Damerow, Murtaugh and Burggraf (1972), NASA CR-120883, Eq. 18:
// f = 2.06 (X_T/D)^-1.1 Re_D^-0.16, X_T the transverse pitch, defined per
// (N - 1) row pairs on total pressure (converted here to the canonical
// per-row basis). Square-diagonal staggered arrays, X_T/D 4.24 and 7.07,
// H/D 2 and 4 (no height effect), 10 rows. Friction rises above inlet Mach
// 0.36 (not modelled).
PinFinFrictionSet damerow_1972_staggered_friction();

// Chyu, Hsing, Shih and Natarajan (1998), ASME 98-GT-175, as tabulated in
// Han, Dutta and Ekkad (2012) Table 4.7: Nu / Pr^0.4 = a Re^b by arrangement
// and surface. Geometry S/D = X/D = 2.5, H/D = 1 (Lyall 2006, Table 2-1),
// 7 rows; naphthalene mass transfer via the heat/mass analogy. The Re basis
// (D, Vmax) is ASSUMED from Lyall's Re_d labelling and the cross-check with
// Metzger, not read from the paper.
PinFinNuSet chyu_1998_nu(PinArrangement arrangement, PinNuSurface surface);

// The inline-over-staggered Total Nu ratio from the two Chyu (1998) sets,
// exact algebra: (0.068 / 0.320) Re_D^(0.733 - 0.583). Same rig, same
// geometry; no friction ratio (has_f = false).
PinFinRatioModifier chyu_1998_inline_over_staggered();

// Chyu (1990), J. Heat Transfer 112, 926, Table 2: Sh/Sc^0.4 = Nu/Pr^0.4 =
// A Re^B for straight or pin-endwall-FILLET pins, inline or staggered.
// Geometry H/D = 1, S/D = X/D = 2.5, 7 rows, Re 5e3 to 3e4, Re = Umax D/nu
// with Umax = Q/A_min (Eqs 8-9), naphthalene (Sh uncertainty 4%). Measured on
// the PIN surfaces; Table 1 finds endwall/pin Sh 0.89-1.09, so the paper
// treats it as the array value -- the surface is recorded as Pin.
PinFinNuSet chyu_1990_nu(PinArrangement arrangement, bool fillet);

// Chyu (1990) Fig. 6, per-row friction f = 2 dP / (rho Umax^2 N) (Eq. 16),
// digitised (no correlation is printed) and fitted, so provenance is Fitted:
//   inline straight    0.1693           (no Re trend; 3 points, 1.4% rms)
//   inline fillet      0.685 Re^-0.139  (3 points)
//   staggered straight 1.616 Re^-0.187  (8 points, 2.0% rms)
//   staggered fillet   9.387 Re^-0.363  (7 points, 1.9% rms)
// converted to the canonical basis by 1/4. The digitised points are in
// validation/cooling/data/chyu1990/. The inline straight set is the default
// inline friction; it is the ONLY inline friction source in hand. Chyu's
// staggered straight values sit 16-25% below Metzger (1982) at the same
// geometry -- a real disagreement between the two labs, reported, not
// reconciled.
PinFinFrictionSet chyu_1990_friction(PinArrangement arrangement, bool fillet);

// Fillet over straight-pin Nu from Chyu (1990) Table 2, exact algebra; same
// arrangement in and out. Staggered 0.339 Re^0.097 (0.78 / 0.83 / 0.92 at
// Re 5e3 / 1e4 / 3e4, the paper's -25 / -17 / -8%); inline about 0.97-1.
// No friction ratio (use chyu_1990_friction with fillet = true).
PinFinRatioModifier chyu_1990_fillet_over_straight(PinArrangement arrangement);

}  // namespace cooling
}  // namespace combaero
