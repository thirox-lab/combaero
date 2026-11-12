#pragma once

// Parametrised rib-roughness correlations.
//
// Rib correlations in the literature share one algebraic shape: a constant
// times power-law terms in the geometry, with the heat-transfer function
// carrying an extra term in the roughness Reynolds number.
//
//   R = C_R * (e/D / nD)^a * (p/e / nP)^b * (W/H / nW)^c * (alpha/90)^d
//   G = C_G * (same geometry terms) * (e+)^n
//
// Expressing them as data rather than code means a new correlation is a
// parameter set, and a user can supply their own -- which is the expected path
// for real hardware, since no published correlation is precise enough for a
// specific rig.
//
// THE NORMALISER IS PART OF THE DATA, NOT A CONVENTION. (p/e/10)^b and
// (p/e)^b are the same function with different constants: Han's
// R = 3.2 (p/e/10)^0.35 is identically 1.4294 (p/e)^0.35. Applying one
// source's constant under another's convention is wrong by a constant factor
// -- 10^0.35 = 2.24x here -- which never looks like a trend and so survives
// review. Every term therefore carries its own reference value.
//
// R CARRIES NO e+ TERM by construction, so the four-sided friction factor
// this family produces is independent of Reynolds number. That is true of
// both correlations in scope (Han 1988, Rallabandi et al. 2009). A channel
// with smooth walls as well still depends on Re through their f_s -- see
// rib_channel_friction. A correlation whose R varies with e+ cannot be
// expressed here and would need its own implicit solve.

#include <limits>
#include <string>
#include <tuple>

namespace combaero {
namespace cooling {

// How much a parameter set's numbers can be trusted, and on whose authority.
enum class RibProvenance {
  // An equation printed in a source, transcribed and checked against it.
  Extracted,
  // A form WE chose, fitted to data WE digitised from a published figure.
  // The exponents are ours; the source only supplied points.
  Fitted,
  // Supplied by the caller. Carries no claim at all.
  User,
};

// What SHAPE a rib is, beyond its angle. Continuous-angled, V, crossed and
// Lambda ribs at the same alpha are different geometries with different
// secondary-flow structure and measurably different G (Han, Zhang and Lee
// 1991, Table 2: at 60 deg, G spans 1.35 to 3.24 times (e+)^b across shapes),
// so a set binds one shape and a series of another shape is refused rather
// than scored -- aliasing one shape onto another reads as model error.
//
//   Transverse  90 deg, continuous across the wall
//   Parallel    angled; the ribs on the two walls are parallel
//   Crossed     angled; the ribs on the two walls cross each other
//   V           pointing downstream (forward V)
//   Lambda      pointing upstream -- a V with the flow reversed, which is
//               why a V or Lambda set is not `symmetric`
//
// Unspecified binds nothing: it is what a user-built set gets by default,
// and the user then owns the shape question.
enum class RibShape { Unspecified, Transverse, Parallel, Crossed, V, Lambda };

// One power-law term. `reference` is what the variable is divided by before
// the exponent is applied; 1.0 means the raw value.
struct RibTerm {
  double exponent = 0.0;
  double reference = 1.0;
};

// An advisory range. Outside it the evaluator warns; it never refuses, because
// a user supplying their own coefficients knows their own hardware and a band
// belongs to the source's rig rather than to theirs.
//
// hi == 0 means no upper bound (a half-open range, or no range at all when lo
// is 0 too). A single tested value is a CLOSED range with hi == lo: Han, Zhang
// and Lee (1991) tested e/D 0.0625 only, so e/D 0.1 is outside it. Testing
// hi > lo here read every point range as half-open and never flagged a value
// above the point.
struct RibRange {
  double lo = 0.0;
  double hi = 0.0;
  // Geometry is quoted to three significant figures (Taslim AR 3.5: W/H
  // 0.286 against the set's 1/3.5), so a bound is read with this relative
  // slack. Without it a point range flags its own paper's data.
  static constexpr double rel_tol = 5.0e-3;
  bool bounded() const { return hi > 0.0 && hi >= lo; }
};

// WHERE AN ACCURACY FIGURE CAME FROM, and it decides what may be done with
// it. This distinction is not bookkeeping: judging a model against a band
// derived from that same model's error is circular, and it is the defect
// #389 exists to remove. It was removed from the dataset's `uncertainty`
// field in #415 and it was still here, one level down.
//
// Of the three shipped sets, only han_1988_orthogonal's figures are the
// author's own claim. han_park_1988_angled's were measured by this project
// THROUGH evaluate_rib -- through the very code path they would be judging
// -- and rallabandi_2009_high_re's from the printed equation against
// digitised data. Both are useful information and neither is a band.
enum class AccuracyProvenance {
  // No figure available. `value` is NaN and must not be read.
  Unstated,
  // The source states it ("95% of data within 6%"). The ONLY kind that may
  // be used as a band to judge a model against.
  Stated,
  // This project measured it against digitised data. Reportable, never a
  // judging band -- see usable_as_band().
  Measured,
};

struct StatedAccuracy {
  // NaN unless `provenance` says otherwise, so an unstated figure read by
  // accident propagates loudly rather than reading as "perfect agreement".
  // That is what the previous 0.0 sentinel did.
  double value = std::numeric_limits<double>::quiet_NaN();
  AccuracyProvenance provenance = AccuracyProvenance::Unstated;

  // Whether this figure may be used as the band a model is judged against.
  // True only for Stated: a Measured figure is the model's own error, and
  // scoring a model inside its own error answers nothing.
  bool usable_as_band() const {
    return provenance == AccuracyProvenance::Stated;
  }

  static StatedAccuracy stated(double v) {
    return {v, AccuracyProvenance::Stated};
  }
  static StatedAccuracy measured(double v) {
    return {v, AccuracyProvenance::Measured};
  }
  static StatedAccuracy unstated() { return {}; }
};

struct RibCorrelationSet {
  std::string name;
  // Where the coefficients came from.
  std::string source;
  // Where the VALIDITY came from. Usually the same, but not always: a set may
  // borrow a range from a later study that showed the correlation holds
  // further than its author claimed.
  std::string validity_source;
  RibProvenance provenance = RibProvenance::User;

  // False when reversing the flow changes the configuration rather than only
  // its direction -- a forward V-rib becomes an inverted V, which the
  // literature measures as a different geometry. True for 90 deg orthogonal
  // and parallel angled ribs, where reversal is a mirror operation.
  bool symmetric = true;

  // Rib shape this set was fitted to. See RibShape.
  RibShape shape = RibShape::Unspecified;

  // Friction roughness function. Two shapes exist in the sources:
  //
  //   PowerLaw       R = C_R * (e/D)^.. * (p/e)^.. * (W/H)^.. * (alpha/90)^d
  //                  (Han 1988, Rallabandi 2009 -- both have d = 0)
  //
  //   QuadraticAlpha R / [(p/e/10)^0.35 * (W/H)^m] = R_quad_c0 +
  //                      R_quad_c1*(alpha/90) + R_quad_c2*(alpha/90)^2
  //                  where m switches on the rib angle itself:
  //                      m = R_quad_WH_exponent_at_90   if alpha == 90 deg
  //                      m = R_quad_WH_exponent_off_90  otherwise
  //                  (Han and Park 1988, Eq. 4.17). R_eD, R_pe, R_WH and
  //                  R_alpha are IGNORED in this shape -- the quadratic
  //                  fields replace them entirely, they do not compose.
  //
  // The switch on m is a genuine discontinuity the source states, not a
  // numerical artefact: at W/H = 2 it is a 27.5% jump in R exactly at
  // alpha = 90, at W/H = 4 it is 62.5%. No smooth interpolation is given by
  // the source, so none is invented here -- see evaluate_rib's comment for
  // what that means for a solver that traverses this exact angle.
  //
  //   QuadraticAlphaTwoBand  R = R_quad_c0 + R_quad_c1*(alpha/90) +
  //                              R_quad_c2*(alpha/90)^2          if W/H >=
  //                              R_WH_band_boundary ("wide" sub-band,
  //                              1/2 <= W/H < 1)
  //                          R = R_quad_narrow_c0 + R_quad_narrow_c1*(a/90)
  //                              + R_quad_narrow_c2*(a/90)^2     otherwise
  //                              ("narrow" sub-band, 1/4 < W/H < 1/2)
  //                  (Han et al. 1989 narrow-channel extension, Eq. 4.19 /
  //                  Fig. 4.48a). No printed equation for R exists at
  //                  W/H < 1 -- these two quadratics are THIS PROJECT's own
  //                  fit to the figure's drawn correlation line (same
  //                  functional shape as Eq. 4.17's QuadraticAlpha, chosen
  //                  because it reproduces that line to <1% RMS), not an
  //                  extraction. See validation/cooling/extractions/
  //                  han_ribbed.md, item 40 and decision D6. R_eD, R_pe,
  //                  R_WH and R_alpha are IGNORED in this shape too -- there
  //                  is no p/e or W/H power-law term; W/H only selects
  //                  which quadratic applies.
  enum class RAlphaShape { PowerLaw, QuadraticAlpha, QuadraticAlphaTwoBand };
  RAlphaShape R_alpha_shape = RAlphaShape::PowerLaw;
  double C_R = 0.0;
  RibTerm R_eD, R_pe, R_WH, R_alpha;
  double R_quad_c0 = 0.0, R_quad_c1 = 0.0, R_quad_c2 = 0.0;
  double R_quad_WH_exponent_at_90 = 0.0;
  double R_quad_WH_exponent_off_90 = 0.0;
  // Eq. 4.17's own cap: "if W/H > 2, set W/H = 2". 0 means uncapped.
  double R_quad_WH_cap = 0.0;
  // QuadraticAlphaTwoBand only: the second (narrow) sub-band's quadratic,
  // and the W/H value that switches between it and R_quad_c0/c1/c2.
  double R_quad_narrow_c0 = 0.0, R_quad_narrow_c1 = 0.0,
         R_quad_narrow_c2 = 0.0;
  double R_WH_band_boundary = 0.0;

  // Heat-transfer roughness function. G_alpha and G_pe are either fixed
  // constants (Fixed, the default -- Han 1988, Rallabandi 2009) or switch
  // on whether the channel is square (Han and Park 1988, Eq. 4.18):
  //
  //   m = G_shape_alpha_exponent_square, n = G_shape_pe_exponent_square
  //       if W/H == 1 (square)
  //   m = G_shape_alpha_exponent_rect,   n = G_shape_pe_exponent_rect
  //       otherwise (rectangular)
  //
  // G_alpha and G_pe are IGNORED when G_shape_model is SquareVsRectangular.
  // Same caveat as R's switch: a genuine, unsmoothed discontinuity at
  // W/H = 1, up to 27% at alpha = 30, p/e = 20.
  //
  //   NarrowChannelAlphaSwitch (Han et al. 1989, Eq. 4.19, W/H < 1):
  //
  //     G = C * (e+)^n, with
  //       C = G_narrow_C_alpha90    if alpha == 90 deg
  //       C = G_narrow_C_off_axis   otherwise (30 < alpha < 90)
  //       n = G_eplus_exponent (the base 0.35)
  //     and, only below G_narrow_WH_band_boundary (the narrow sub-band,
  //     1/4 < W/H < 1/2):
  //       C *= (W/H)^G_narrow_WH_C_exponent
  //       n *= (W/H)^G_narrow_WH_n_exponent
  //
  //   This is a genuine switch the text states directly (page 378: "C =
  //   2.24 if alpha = 90 deg, and C = 1.80 if 30 deg < alpha < 90 deg"),
  //   not a numerical artefact -- a 24% jump in C at alpha = 90 for any
  //   W/H. Unlike SquareVsRectangular, the switch is on alpha and the
  //   constant itself, not on channel shape and the alpha/p-e exponents;
  //   it needs its own fields because no existing combination expresses an
  //   exponent that is itself a function of W/H. G_alpha, G_pe, G_WH and
  //   G_eD are IGNORED in this shape.
  enum class GShapeModel { Fixed, SquareVsRectangular, NarrowChannelAlphaSwitch };
  GShapeModel G_shape_model = GShapeModel::Fixed;
  double C_G = 0.0;
  RibTerm G_eD, G_pe, G_WH, G_alpha;
  double G_eplus_exponent = 0.0;
  double G_shape_alpha_exponent_square = 0.0, G_shape_alpha_exponent_rect = 0.0;
  double G_shape_pe_exponent_square = 0.0, G_shape_pe_exponent_rect = 0.0;
  // NarrowChannelAlphaSwitch only.
  double G_narrow_C_alpha90 = 0.0, G_narrow_C_off_axis = 0.0;
  double G_narrow_WH_band_boundary = 0.0;
  double G_narrow_WH_C_exponent = 0.0, G_narrow_WH_n_exponent = 0.0;

  // The FOUR-WALL average heat-transfer roughness function, when the source
  // prints its own: G_bar = C_Gbar (e+)^Gbar_eplus_exponent. C_Gbar = 0 means
  // unstated, and RibResult::has_G_bar is then false -- callers fall back to
  // Han's published G_bar = 1.2 G. A source that prints G_bar per
  // configuration (Han, Zhang and Lee 1991: G_bar/G 1.07-1.26 across nine)
  // must not have its own value replaced by that single constant.
  double C_Gbar = 0.0;
  double Gbar_eplus_exponent = 0.0;

  // Advisory validity.
  RibRange valid_Re, valid_eD, valid_pe, valid_WH, valid_alpha, valid_eplus;
  double valid_Pr = 0.0;  // 0 means unstated

  // Accuracy of the correlation, as a fraction (0.06 for "within 6%"),
  // carried WITH where the number came from. See AccuracyProvenance: two of
  // the three shipped sets carry this project's own measurement here, not an
  // author claim, and a band measured from a model must never be used to
  // judge that model.
  StatedAccuracy accuracy_R, accuracy_G;
};

// Han, J.C. (1988), ASME J. Heat Transfer 110, 321, for 90 deg orthogonal ribs
// in two-opposite-wall rectangular channels. Extracted and confirmed; see
// validation/cooling/extractions/han_ribbed.md.
RibCorrelationSet han_1988_orthogonal();

// Rallabandi, A.P., Yang, H. and Han, J.-C. (2009), ASME J. Heat Transfer
// 131(7), 071703, for 45 deg square/sharp-edged ribs in a square channel at
// Reynolds numbers an order of magnitude above han_1988_orthogonal's range --
// "typical of land-based turbines" per the source. A separate regime, not a
// revision: different experiments, and per the source the correlations do
// not agree with Han's in the extended range. Selected explicitly, never
// blended -- see validation/cooling/extractions/han_ribbed_high_re.md.
//
// Sharp-edged ribs only. The source reports round-edged ribs at the same
// conditions instead following Han's correlation, but "coincidentally" --
// its own word, meaning reported as fortuitous rather than physically
// grounded, so it is not encoded as a switch here. A caller with round-edged
// ribs above this Reynolds range chooses the correlation directly, same as
// any other set choice.
RibCorrelationSet rallabandi_2009_high_re();

// Han, J.C. and Park, J.S. (1988), IJHMT 31(1), 183, Eq. 4.17/4.18, for
// broad-aspect-ratio rectangular ducts with angled ribs -- alpha 30-90 deg,
// W/H 1-4. A different paper and configuration from han_1988_orthogonal
// (which is 90 deg only): the two R correlations happen to agree closely
// at alpha=90 (3.44%), but the two G correlations disagree by up to 22%
// there, because Fig. 4.46 and this set describe genuinely different rib
// configurations that are not obliged to agree. See
// validation/cooling/extractions/han_ribbed.md, item 23.
RibCorrelationSet han_park_1988_angled();

// Han, J.C. et al. (1989), Eq. 4.19, via Han, Dutta & Ekkad (2012) 2nd ed.
// Fig. 4.48 -- extends the same rib family to narrow-aspect-ratio channels,
// 1/4 < W/H < 1, alpha 30-90 deg. A different W/H domain from
// han_1988_orthogonal and han_park_1988_angled (both 1-4), not a revision of
// either: W/H = 1 stays with the wide sets by the same figure-driven
// convention that resolved Fig. 4.48's own branch-boundary discontinuity
// (validation/cooling/extractions/han_ribbed.md, item 41).
//
// G is text-extracted (items 36-38): a hard switch on alpha (C = 2.24 at
// 90 deg, 1.80 off-axis) with a further W/H correction below W/H = 1/2.
//
// R has NO PRINTED EQUATION in this source (item 40) -- only a drawn
// correlation line in Fig. 4.48a. The coefficients here are THIS PROJECT's
// own fit to that line (same quadratic-in-alpha shape as Eq. 4.17, <1% RMS
// against the drawn curve), not an extraction, so `provenance` is `Fitted`
// for the set as a whole even though G alone would be `Extracted` -- see
// decision D6 in han_ribbed.md for why a single set-level enum was kept
// rather than splitting it per quantity.
RibCorrelationSet han_1989_narrow_channel();

// Han, J.C., Zhang, Y.M. and Lee, C.P. (1991), ASME J. Heat Transfer 113,
// 590, Table 2: R, G and G_bar for nine rib configurations in ONE rig --
// square channel, e/D 0.0625, P/e 10, in-line ribs, Re 15,000-90,000. Table 2
// verified on four channels (validation/cooling/extractions/
// han_zhang_lee_1991_jht.md). Supported configurations:
//
//   Transverse at 90 deg; Parallel, Crossed, V, Lambda at 45 or 60 deg.
//
// Anything else throws std::invalid_argument -- the paper did not test it,
// and a nearby configuration is not a substitute. The validity range is a
// single point in e/D, P/e and W/H; everything off it is extrapolation.
RibCorrelationSet han_zhang_lee_1991(RibShape shape, double alpha_deg);

// Geometry of the ribbed channel, as the correlation sees it.
struct RibGeometry {
  double e_D = 0.0;        // rib height / hydraulic diameter
  double p_e = 0.0;        // rib pitch / rib height
  double W_H = 1.0;        // channel width / height
  double alpha_deg = 90.0; // rib angle to the flow
};

// The chain, evaluated. Every field is finite for every real input: the
// evaluator never throws from inside a residual, because the solver
// legitimately probes states that are not physical and a throw there kills the
// solve. Bad PARAMETERS are rejected once, by validate_rib_set; bad OPERATING
// POINTS are guarded smoothly.
struct RibResult {
  double R = 0.0;        // friction roughness function
  // Four-sided channel friction factor, FANNING (see rib_channel_friction
  // for a channel's Darcy value with smooth walls as well).
  double f = 0.0;
  double e_plus = 0.0;   // roughness Reynolds number, signed with the flow
  double G = 0.0;        // heat-transfer roughness function
  double St_r = 0.0;     // ribbed-side Stanton number
  // d(St_r)/d(Re). df/dRe of the four-sided f is identically zero for this
  // family -- R carries no e+ term -- so it is not reported.
  double dSt_dRe = 0.0;
  // Four-wall average G, only when the set prints it (C_Gbar > 0); 0 and
  // false otherwise. Evaluated at the same smoothed e+ as G.
  double G_bar = 0.0;
  bool has_G_bar = false;
  bool extrapolated = false;  // outside the set's advisory validity
};

// Evaluate at a Reynolds number. `Re` may be negative (reverse flow) or zero;
// the guards are smooth through both, so the derivative is continuous where
// Newton iterates. `geom.alpha_deg` is NOT guarded the same way: a set's R/G
// terms are polynomials in alpha with no built-in symmetry, so an angle
// outside `valid_alpha` is evaluated raw rather than folded or clamped. See
// validation/cooling/extractions/han_ribbed.md, "Eq. 4.17 outside its stated
// 30-90 degree range" for the numerical characterisation this rests on: the
// raw polynomial stays finite and positive well past the fitted range, but
// its vertex sits at 68 deg (not 90), so it is unbounded and physically
// unmotivated beyond 30-90, not a safe extrapolation. No guard exists yet;
// if one is added, reflection-folding (alpha_eff = 180 - alpha) is the
// bounded, `symmetric`-consistent choice documented there, not the raw value.
RibResult evaluate_rib(const RibCorrelationSet &set, const RibGeometry &geom,
                       double Re);

// Reject a set that cannot be evaluated. Throws with a message naming the
// field. Unlike the validity ranges, these are hard errors: they are mistakes,
// not operating points.
void validate_rib_set(const RibCorrelationSet &set);

// ---------------------------------------------------------------------------
// Channel friction for a duct with some of its walls ribbed (#1106)
// ---------------------------------------------------------------------------
//
// RibResult::f is Han's FOUR-SIDED, FANNING friction factor. Both halves of
// that sentence are definitions in the source:
//
//   FANNING. Han and Park (1988), p. 187, define the measured passage
//   friction as fbar = dP / {4 (dL/D) [G^2 / (2 rho g_c)]}, i.e. Fanning
//   (tau_w = f rho V^2 / 2); their smooth reference f(FD) = 0.046 Re^-0.2
//   (Eq. 8) is the Fanning Blasius form. A Darcy drop needs 4 f.
//
//   FOUR-SIDED. The measured fbar is for two OPPOSITE ribbed walls (width W)
//   and two smooth ones (height H). Han converts it to the all-ribbed channel
//   by area weighting, f = fbar + (H/W)(fbar - f_s) (Eq. 6, "weighted by
//   area"), with f_s = f(FD). R is formed from that f.
//
// So the friction of a channel with a ribbed perimeter fraction phi is the
// same area weighting run forwards:
//
//   f_channel = phi f + (1 - phi) f_s(Re),   f_D = 4 f_channel   (Darcy)
//
// which for two opposite ribbed walls (phi = W/(W+H)) is exactly Han's
// measured fbar again, and for four ribbed walls (phi = 1) is f. One ribbed
// wall (phi = W/(2(W+H))) is the same weighting outside the configuration
// Han measured -- see rib_ribbed_wall_fraction.

// Han's smooth-channel friction f(FD) = 0.046 Re^-0.2, FANNING (Han and Park
// 1988 Eq. 8; Han, Zhang and Lee 1991 Eq. 2): the f_s of the decomposition.
inline constexpr double HAN_SMOOTH_FANNING_COEF = 0.046;
inline constexpr double HAN_SMOOTH_FANNING_RE_EXP = -0.2;
// Darcy f_D = 4 f_F (Fanning): wall shear tau = f_F rho V^2 / 2.
inline constexpr double DARCY_PER_FANNING = 4.0;

// Ribbed fraction of the wetted perimeter of a W x H duct with
// `n_ribbed_walls` ribbed walls, ribs on the walls of width W:
//   n = 2 (two opposite walls, Han's configuration):  W/(W+H)
//   n = 4:                                             1
//   n = 1 (one wall of width W):                       W/(2(W+H))
// Throws std::invalid_argument for any other n (three ribbed walls has no
// unambiguous geometry) or a non-positive W_H. A configuration choice, so it
// is validated, not guarded.
double rib_ribbed_wall_fraction(double W_H, int n_ribbed_walls);

// (f_D, df_D/dRe): the DARCY friction factor of the whole channel,
// f_D = 4 [phi f + (1 - phi) f_s(Re)], for the drop
// dP = f_D (L/D_h) rho V^2 / 2. `Re` may be negative or zero: f_D is even in
// Re and its derivative odd, with f_s evaluated at a smoothly floored |Re| so
// both stay finite through zero flow. f itself carries no Re dependence (R
// has no e+ term), so the derivative is the smooth walls' alone and vanishes
// for n = 4.
std::tuple<double, double> rib_channel_friction(const RibCorrelationSet &set,
                                                const RibGeometry &geom,
                                                double Re, int n_ribbed_walls);

}  // namespace cooling
}  // namespace combaero
