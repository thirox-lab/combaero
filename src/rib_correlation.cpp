#include "rib_correlation.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>

namespace combaero {
namespace cooling {

RibCorrelationSet han_1988_orthogonal() {
  RibCorrelationSet s;
  s.name = "han_1988_orthogonal";
  s.source = "Han, J.C. (1988). ASME J. Heat Transfer 110, 321, via Han, "
             "Dutta & Ekkad (2012) 2nd ed. Eq. 4.15/4.16, Fig. 4.46";
  s.validity_source = s.source;
  s.provenance = RibProvenance::Extracted;
  // 90 deg orthogonal ribs: reversing the flow leaves the geometry unchanged.
  s.symmetric = true;
  s.shape = RibShape::Transverse;

  // R = 3.2 (p/e / 10)^0.35, independent of e+.
  s.C_R = 3.2;
  s.R_pe = {0.35, 10.0};

  // G = 3.7 (e+)^0.28 at Pr ~ 0.7. No geometry terms: Han reports G as a
  // function of e+ alone for this configuration.
  s.C_G = 3.7;
  s.G_eplus_exponent = 0.28;

  s.valid_Re = {10000.0, 60000.0};
  s.valid_eD = {0.047, 0.078};
  s.valid_pe = {10.0, 20.0};
  s.valid_WH = {1.0, 4.0};
  s.valid_alpha = {90.0, 90.0};
  // The correlation is stated for e+ >= 50; no upper bound is given.
  s.valid_eplus = {50.0, 0.0};
  s.valid_Pr = 0.7;

  // Han's own claims, and the only STATED pair in this file -- so the only
  // one the scorecard may judge a model against.
  s.accuracy_R = StatedAccuracy::stated(0.06);  // "95% of data within 6%"
  s.accuracy_G = StatedAccuracy::stated(0.08);  // "95% of data within 8%"
  return s;
}

RibCorrelationSet rallabandi_2009_high_re() {
  RibCorrelationSet s;
  s.name = "rallabandi_2009_high_re";
  s.source = "Rallabandi, A.P., Yang, H. and Han, J.-C. (2009). ASME J. Heat "
             "Transfer 131(7), 071703. Eq. (17)/(18)";
  s.validity_source = s.source;
  s.provenance = RibProvenance::Extracted;
  // 45 deg square/sharp-edged ribs, same family as Han and Park's angled
  // ribs: reversing the flow gives the mirror image, same magnitude in a
  // symmetric duct. Not orthogonal (unlike han_1988_orthogonal), but still
  // symmetric -- see the case A/B/C reverse-flow discussion on #334.
  s.symmetric = true;
  s.shape = RibShape::Parallel;

  // R = 1.13 (e/D)^-0.17 (p/e)^0.38. Raw e/D and p/e, no /10 normaliser --
  // unlike han_1988_orthogonal's R_pe, which divides by 10 because that is
  // how Han's Eq. 4.16 is printed. Two different papers, two different
  // normalisations; carrying han_1988_orthogonal's convention here would be
  // silently wrong by (p/e)^0.38 evaluated at the wrong reference.
  s.C_R = 1.13;
  s.R_eD = {-0.17, 1.0};
  s.R_pe = {0.38, 1.0};

  // G = 1.24 (e/D)^0.014 (p/e)^-0.02 (e+)^0.42. The exponents on e/D and p/e
  // are the resolved reading (extraction item 31): the textbook reprint at
  // Fig. 4.193c prints 0.14 for the e/D exponent, a lost decimal place
  // against the paper's own 0.014, worth 22% in G. Confirmed against the
  // paper directly, not the reprint -- see han_ribbed_high_re.md, D4.
  s.C_G = 1.24;
  s.G_eD = {0.014, 1.0};
  s.G_pe = {-0.02, 1.0};
  s.G_eplus_exponent = 0.42;

  s.valid_Re = {30000.0, 400000.0};
  s.valid_eD = {0.1, 0.18};
  s.valid_pe = {5.0, 10.0};
  // Square channel only (W/H = 1) and 45 deg ribs only: both are how the
  // source's own experiments were run, not modelling choices, so they are
  // recorded as point ranges the same way han_1988_orthogonal pins alpha.
  s.valid_WH = {1.0, 1.0};
  s.valid_alpha = {45.0, 45.0};
  // Derived, not separately stated: the e+ reached at the corners of the
  // Re x e/D x p/e validity box above, through this set's own R -> f -> e+
  // chain. The source's item-8 anchor (e+ = 18,000 at Re = 400K, e/D = 0.18)
  // falls inside this box; it does not fix p/e, which this box does span.
  s.valid_eplus = {542.0, 25340.0};
  s.valid_Pr = 0.7;  // unchanged from Han; the source does not restate it

  // MEASURED, not stated. Rallabandi et al. state no percentage accuracy for
  // Eq. (18) in the extracted text. This is the RMS this project measured,
  // of Eq. (18) itself against 38 digitised points from Fig. 4.193c
  // (validation/cooling/data/han2012/fig4.193c_G_scatter.csv): mean
  // pred/data 0.994, RMS 6.9% -- see han_ribbed_high_re.md, items 27-31.
  // Reportable, but never a band to judge the model against, because it IS
  // the model's error.
  s.accuracy_G = StatedAccuracy::measured(0.069);
  // Unstated: no R data was digitised for this range, so there is nothing
  // to measure it against either. Previously 0.0, which read as "perfect".
  s.accuracy_R = StatedAccuracy::unstated();
  return s;
}

RibCorrelationSet han_park_1988_angled() {
  RibCorrelationSet s;
  s.name = "han_park_1988_angled";
  s.source = "Han, J.C. and Park, J.S. (1988). IJHMT 31(1), 183. Eq. "
             "4.17/4.18, via Han, Dutta & Ekkad (2012) 2nd ed. Fig. 4.47";
  s.validity_source = s.source;
  s.provenance = RibProvenance::Extracted;
  // Parallel angled ribs: reversing the flow gives the mirror image, same
  // magnitude in a symmetric duct, with the secondary-flow direction
  // flipped -- not the reversal case that changes the configuration (a
  // V-shaped rib reversed becomes an inverted V, a different measured
  // geometry). See the case A/B/C reverse-flow discussion on #334.
  s.symmetric = true;
  s.shape = RibShape::Parallel;

  // R/[(P/e/10)^0.35 (W/H)^m] = 12.31 - 27.07(alpha/90) + 17.86(alpha/90)^2
  // m = 0 at alpha = 90 deg (any W/H collapses to 1), m = 0.35 otherwise.
  // The leading constant is 12.31 (Fig. 4.47), not p. 376's printed 12.3 --
  // decision D2 in han_ribbed.md, accepted in favour of the figure.
  s.R_alpha_shape = RibCorrelationSet::RAlphaShape::QuadraticAlpha;
  s.R_quad_c0 = 12.31;
  s.R_quad_c1 = -27.07;
  s.R_quad_c2 = 17.86;
  s.R_pe = {0.35, 10.0};
  s.R_quad_WH_exponent_at_90 = 0.0;
  s.R_quad_WH_exponent_off_90 = 0.35;
  s.R_quad_WH_cap = 2.0;  // Eq. 4.17: "if W/H > 2, set W/H = 2"

  // G = 2.24 (W/H)^0.1 (alpha/90)^m (p/e/10)^n (e+)^0.35. (W/H)^0.1 sits
  // OUTSIDE the square/rectangular switch and applies always -- carried on
  // G_WH, which the switch does not touch. m, n switch on channel shape:
  // square (W/H=1) gets m=0.35, n=0.1; rectangular gets m=n=0, the stated
  // consequence being that alpha and P/e stop mattering for G once the
  // channel is not square.
  s.C_G = 2.24;
  s.G_WH = {0.1, 1.0};
  s.G_shape_model = RibCorrelationSet::GShapeModel::SquareVsRectangular;
  s.G_shape_alpha_exponent_square = 0.35;
  s.G_shape_alpha_exponent_rect = 0.0;
  s.G_shape_pe_exponent_square = 0.1;
  s.G_shape_pe_exponent_rect = 0.0;
  s.G_eplus_exponent = 0.35;

  s.valid_Re = {10000.0, 60000.0};
  s.valid_eD = {0.047, 0.078};
  s.valid_pe = {10.0, 20.0};
  s.valid_WH = {1.0, 4.0};
  s.valid_alpha = {30.0, 90.0};
  s.valid_eplus = {50.0, 0.0};  // e+ >= 50, same floor as han_1988_orthogonal
  s.valid_Pr = 0.7;

  // MEASURED, and the most circular of the three: these are the RMS this
  // project measured scoring THROUGH evaluate_rib against figure 4.47's own
  // digitised data (39 points for R, 115 for G) with a representative
  // geometry (the cloud carries no legend) -- so they are the error of the
  // very code path they would be judging. Reportable, never a band. See
  // validation/cooling/extractions/han_ribbed.md and
  // validation/cooling/data/han2012/fig4.47_*.csv.
  s.accuracy_R = StatedAccuracy::measured(0.105);
  s.accuracy_G = StatedAccuracy::measured(0.088);
  return s;
}

RibCorrelationSet han_1989_narrow_channel() {
  RibCorrelationSet s;
  s.name = "han_1989_narrow_channel";
  s.source = "Han, J.C. et al. (1989), Eq. 4.19 (G), via Han, Dutta & Ekkad "
             "(2012) 2nd ed. Fig. 4.48 (R, this project's own fit -- no "
             "printed R equation exists for W/H < 1)";
  s.validity_source = s.source;
  // R is FITTED by this project to Fig. 4.48a's drawn correlation line (no
  // printed equation exists -- han_ribbed.md item 40); G is text-extracted
  // (items 36-38). One enum per set, so the more cautious value is kept --
  // see decision D6 and the header comment on this factory.
  s.provenance = RibProvenance::Fitted;
  // Parallel angled ribs, same reasoning as han_park_1988_angled: reversal
  // gives the mirror image, same magnitude.
  s.symmetric = true;
  s.shape = RibShape::Parallel;

  // R: no printed equation for W/H < 1 (item 40). Fitted to Fig. 4.48a's own
  // drawn correlation line, one quadratic-in-alpha per W/H sub-band, same
  // shape Eq. 4.17 already uses. RMS against the line itself: 0.99% (wide
  // sub-band, validation/cooling/data/han2012/fig4.48_R_correlation_WH0.5.csv)
  // and 0.62% (narrow sub-band, ..._WH0.25.csv).
  s.R_alpha_shape = RibCorrelationSet::RAlphaShape::QuadraticAlphaTwoBand;
  s.R_quad_c0 = 10.66406939;
  s.R_quad_c1 = -21.43911537;
  s.R_quad_c2 = 15.49713879;
  s.R_quad_narrow_c0 = 8.79873413;
  s.R_quad_narrow_c1 = -14.76786182;
  s.R_quad_narrow_c2 = 12.4724873;
  s.R_WH_band_boundary = 0.5;

  // G = C(e+)^n, Eq. 4.19 (p. 378): C = 2.24 at alpha = 90 deg, 1.80
  // off-axis; below W/H = 1/2 both C and n pick up a (W/H) correction.
  s.G_shape_model = RibCorrelationSet::GShapeModel::NarrowChannelAlphaSwitch;
  // C_G is unused by this shape -- G_narrow_C_alpha90/off_axis replace it --
  // and exempt from validate_rib_set's positivity check, so it is left at
  // its default rather than carrying a dummy value nothing reads.
  s.G_eplus_exponent = 0.35;  // base n, before the narrow-band (W/H)^0.44
  s.G_narrow_C_alpha90 = 2.24;
  s.G_narrow_C_off_axis = 1.80;
  s.G_narrow_WH_band_boundary = 0.5;
  s.G_narrow_WH_C_exponent = -0.76;
  s.G_narrow_WH_n_exponent = 0.44;

  s.valid_Re = {10000.0, 60000.0};
  s.valid_eD = {0.047, 0.078};
  s.valid_pe = {10.0, 20.0};
  // Carried over from han_park_1988_angled's box (same figure family) --
  // the text gives no separate Re/e_D/p_e/Pr range for Eq. 4.19, only the
  // W/H and alpha domain below. Not a source claim; an explicit assumption.
  s.valid_WH = {0.25, 1.0};  // 1/4 < W/H < 1; W/H = 1 stays with the wide
                             // sets (item 41's convention)
  s.valid_alpha = {30.0, 90.0};
  s.valid_eplus = {50.0, 0.0};
  s.valid_Pr = 0.7;

  // MEASURED by this project, never a judging band. R: pooled RMS 5.4%
  // (bias -0.9%/+0.4% per sub-band) against the 16 raw digitised points in
  // fig4.48_R_WH0.25.csv/_WH0.5.csv -- the fitted coefficients above were
  // derived from the drawn LINE, so this is a genuinely separate check
  // against the scatter, same split rallabandi_2009_high_re observes. G:
  // pooled RMS 5.4% (bias +1.9%) against fig4.48_G_WH0.25.csv/_WH0.5.csv,
  // evaluated at alpha = 90 deg uniformly because those series pool
  // multiple rib angles onto one e+-vs-G curve with no per-point angle
  // recorded -- checked, not assumed, that this doesn't hide a 24%-scale
  // defect: a free power-law fit to the pooled scatter ALONE (ignoring the
  // text's alpha switch entirely) already lands within 3.9-5.0% RMS, so the
  // angles are not separable on this axis to begin with.
  s.accuracy_R = StatedAccuracy::measured(0.054);
  s.accuracy_G = StatedAccuracy::measured(0.054);
  return s;
}

namespace {

// Han, Zhang and Lee (1991) JHT 113, 590, Table 2 (p. 595), verified on four
// channels. R = R_a (e+)^0 for every case; G and G_bar are a (e+)^b.
struct HzlRow {
  RibShape shape;
  double alpha_deg;
  const char *tag;
  double R_a, G_a, G_b, Gbar_a, Gbar_b;
};

constexpr HzlRow kHzlTable2[] = {
    {RibShape::Transverse, 90.0, "90", 3.18, 3.97, 0.28, 4.86, 0.28},
    {RibShape::Parallel, 60.0, "60par", 2.05, 1.52, 0.41, 2.41, 0.36},
    {RibShape::Crossed, 60.0, "60crs", 3.18, 3.24, 0.32, 4.57, 0.28},
    {RibShape::V, 60.0, "60vee", 1.72, 1.35, 0.42, 1.76, 0.40},
    {RibShape::Lambda, 60.0, "60lam", 1.42, 1.59, 0.43, 2.12, 0.41},
    {RibShape::Parallel, 45.0, "45par", 3.05, 2.07, 0.36, 3.01, 0.32},
    {RibShape::Crossed, 45.0, "45crs", 4.40, 2.11, 0.37, 2.43, 0.37},
    {RibShape::V, 45.0, "45vee", 2.04, 1.36, 0.43, 1.93, 0.40},
    {RibShape::Lambda, 45.0, "45lam", 1.70, 1.83, 0.41, 2.49, 0.38},
};

}  // namespace

RibCorrelationSet han_zhang_lee_1991(RibShape shape, double alpha_deg) {
  const HzlRow *row = nullptr;
  for (const auto &r : kHzlTable2) {
    if (r.shape == shape && std::abs(r.alpha_deg - alpha_deg) < 1e-9) {
      row = &r;
      break;
    }
  }
  if (row == nullptr) {
    throw std::invalid_argument(
        "han_zhang_lee_1991: no Table 2 row for this shape at alpha = " +
        std::to_string(alpha_deg) +
        " deg. Tested: Transverse at 90; Parallel, Crossed, V, Lambda at 45 "
        "or 60. A nearby configuration is not a substitute.");
  }

  RibCorrelationSet s;
  s.name = std::string("han_zhang_lee_1991_") + row->tag;
  s.source = "Han, J.C., Zhang, Y.M. and Lee, C.P. (1991). ASME J. Heat "
             "Transfer 113, 590. Table 2";
  s.validity_source = s.source;
  s.provenance = RibProvenance::Extracted;
  s.shape = row->shape;
  // A V reversed is a Lambda -- a different measured configuration -- so
  // neither is symmetric under flow reversal. Parallel, crossed and
  // transverse ribs are their own mirror image.
  s.symmetric = !(row->shape == RibShape::V || row->shape == RibShape::Lambda);

  // R = a (e+)^0: independent of e+, the paper's own finding (b = 0 in every
  // row), so it fits the plain power law with no geometry terms.
  s.C_R = row->R_a;
  s.C_G = row->G_a;
  s.G_eplus_exponent = row->G_b;
  s.C_Gbar = row->Gbar_a;
  s.Gbar_eplus_exponent = row->Gbar_b;

  // One rig. Point ranges flag everything off it as extrapolation.
  s.valid_Re = {15000.0, 90000.0};
  s.valid_eD = {0.0625, 0.0625};
  s.valid_pe = {10.0, 10.0};
  s.valid_WH = {1.0, 1.0};
  s.valid_alpha = {row->alpha_deg, row->alpha_deg};
  s.valid_Pr = 0.7;

  // The paper states MEASUREMENT uncertainty (Nu and f under 8%), not an
  // accuracy for its fits, so nothing is claimed here.
  s.accuracy_R = StatedAccuracy::unstated();
  s.accuracy_G = StatedAccuracy::unstated();
  return s;
}

namespace {

void require_positive_reference(const RibTerm &t, const std::string &field) {
  if (!(t.reference > 0.0) || !std::isfinite(t.reference)) {
    throw std::invalid_argument(
        "rib correlation set: " + field +
        ".reference must be finite and positive, got " +
        std::to_string(t.reference) +
        ". A missing reference is indistinguishable from a deliberate 1.0, so "
        "it is rejected rather than defaulted.");
  }
  if (!std::isfinite(t.exponent)) {
    throw std::invalid_argument("rib correlation set: " + field +
                                ".exponent must be finite");
  }
}

}  // namespace

namespace {

// Smooth floor on e+ so the heat-transfer power law is total and its
// derivative continuous through zero. A hard abs() or max() would put a kink
// exactly where Newton iterates; sqrt(x^2 + eps^2) does not. At e+ = 1000 the
// distortion is 1.4e-07.
//
// The floor is only significant below e+ ~ 5, which for a typical rib geometry
// is Re below roughly 700 -- an order of magnitude under the bottom of Han's
// validity (Re = 10,000, e+ >= 50) and a regime where the velocity, the
// pressure drop and the Nusselt number are all negligible anyway. The guard
// exists to keep the solver's arithmetic well behaved while it passes through
// such states, not to predict them.
constexpr double EPLUS_FLOOR = 1.0;

// The bracket (2/f)^(1/2) must stay positive for f to be meaningful. It can
// reach zero for coefficient and geometry combinations a user can supply --
// with C_R = 1.0 it happens at e/D = 0.274, inside plausible input -- so it is
// floored rather than allowed to produce an infinite or negative f.
constexpr double BRACKET_FLOOR = 1e-3;

// The Stanton denominator 1 + (G - R)(f/2)^(1/2) thins as e/D and C_R grow and
// would give a NEGATIVE Stanton number if it crossed zero, which is worse than
// a crash because it looks like a number.
constexpr double ST_DENOM_FLOOR = 1e-3;

// Smooth floor on |Re| in the smooth walls' f_s = 0.046 Re^-0.2, which is
// infinite at zero flow. sqrt(Re^2 + eps^2) keeps f_s and its derivative
// finite and continuous through zero, as EPLUS_FLOOR does for G; the drop
// f_s m|m| still vanishes there. At Re = 1000 the distortion is 1e-7, three
// orders below Han's validity (Re >= 10,000).
constexpr double SMOOTH_RE_FLOOR = 1.0;

double power_term(const RibTerm &t, double value) {
  if (t.exponent == 0.0) {
    return 1.0;
  }
  const double x = value / t.reference;
  if (!(x > 0.0) || !std::isfinite(x)) {
    return 1.0;
  }
  return std::pow(x, t.exponent);
}

double softmin_floor(double x, double floor_value) {
  // Smooth one-sided floor: equals x well above floor_value, approaches
  // floor_value below it, differentiable throughout.
  return 0.5 * (x + std::sqrt(x * x + floor_value * floor_value));
}

bool outside(const RibRange &r, double v) {
  if (!r.bounded()) {
    // A half-open range: lo set, hi unset.
    return r.lo > 0.0 && v < r.lo * (1.0 - RibRange::rel_tol);
  }
  return v < r.lo * (1.0 - RibRange::rel_tol) || v > r.hi * (1.0 + RibRange::rel_tol);
}

// True for an angle within floating-point noise of exactly 90 degrees.
// alpha_deg is a fixed geometry parameter in every current caller, not a
// Newton-iterated state, so an exact-valued switch does not sit in a
// solver's search path today; see evaluate_rib's header comment for what
// changes if that ever stops being true.
bool is_alpha_90(double alpha_deg) {
  return std::abs(alpha_deg - 90.0) < 1e-9;
}

// True for a square channel, within floating-point noise of W/H = 1.
bool is_square_channel(double W_H) { return std::abs(W_H - 1.0) < 1e-9; }

double han_park_R(const RibCorrelationSet &set, double e_D, double W_H,
                  double p_e, double alpha_deg) {
  (void)e_D;  // Eq. 4.17 carries no e/D term.
  const double u = alpha_deg / 90.0;
  const double alpha_poly =
      set.R_quad_c0 + set.R_quad_c1 * u + set.R_quad_c2 * u * u;
  double W_H_capped = W_H;
  if (set.R_quad_WH_cap > 0.0 && W_H_capped > set.R_quad_WH_cap) {
    W_H_capped = set.R_quad_WH_cap;
  }
  const double m = is_alpha_90(alpha_deg) ? set.R_quad_WH_exponent_at_90
                                          : set.R_quad_WH_exponent_off_90;
  const double WH_term =
      (W_H_capped > 0.0 && std::isfinite(W_H_capped) && m != 0.0)
          ? std::pow(W_H_capped, m)
          : 1.0;
  return alpha_poly * power_term(set.R_pe, p_e) * WH_term;
}

// Eq. 4.19's R (han_ribbed.md item 40): no printed equation, so this project
// fitted one quadratic-in-alpha per W/H sub-band to Fig. 4.48a's drawn line.
// Unlike han_park_R there is no p/e or W/H power-law term -- W/H only picks
// which quadratic applies.
double han_narrow_R(const RibCorrelationSet &set, double W_H,
                    double alpha_deg) {
  const double u = alpha_deg / 90.0;
  if (W_H < set.R_WH_band_boundary) {
    return set.R_quad_narrow_c0 + set.R_quad_narrow_c1 * u +
           set.R_quad_narrow_c2 * u * u;
  }
  return set.R_quad_c0 + set.R_quad_c1 * u + set.R_quad_c2 * u * u;
}

// Eq. 4.19's G: a hard switch on alpha == 90 deg (not a continuous function
// of alpha, unlike every other shape in this file), with a further W/H
// correction to both the leading constant and the e+ exponent itself below
// W/H = 1/2. See the GShapeModel::NarrowChannelAlphaSwitch header comment.
double han_narrow_G(const RibCorrelationSet &set, double W_H,
                    double alpha_deg, double e_plus_safe) {
  double C = is_alpha_90(alpha_deg) ? set.G_narrow_C_alpha90
                                    : set.G_narrow_C_off_axis;
  double n = set.G_eplus_exponent;
  if (W_H < set.G_narrow_WH_band_boundary && W_H > 0.0 &&
      std::isfinite(W_H)) {
    C *= std::pow(W_H, set.G_narrow_WH_C_exponent);
    n *= std::pow(W_H, set.G_narrow_WH_n_exponent);
  }
  return C * std::pow(e_plus_safe, n);
}

}  // namespace

RibResult evaluate_rib(const RibCorrelationSet &set, const RibGeometry &geom,
                       double Re) {
  RibResult out;

  const double e_D = geom.e_D;
  const double W_H = geom.W_H;

  if (set.R_alpha_shape ==
      RibCorrelationSet::RAlphaShape::QuadraticAlphaTwoBand) {
    out.R = han_narrow_R(set, W_H, geom.alpha_deg);
  } else if (set.R_alpha_shape ==
             RibCorrelationSet::RAlphaShape::QuadraticAlpha) {
    out.R = han_park_R(set, e_D, W_H, geom.p_e, geom.alpha_deg);
  } else {
    out.R =
        set.C_R * power_term(set.R_eD, e_D) * power_term(set.R_pe, geom.p_e) *
        power_term(set.R_WH, W_H) *
        power_term(set.R_alpha, geom.alpha_deg / 90.0 * set.R_alpha.reference);
  }

  // Invert the wall law for f. The geometry group is
  // (2 e/D) * (2 W/(W + H)); with W_H = W/H it is (2 e/D) * (2 W_H/(W_H + 1)).
  const double geom_group =
      (2.0 * e_D) * (2.0 * W_H / (W_H + 1.0));
  double bracket = 0.0;
  if (geom_group > 0.0 && std::isfinite(geom_group)) {
    bracket = out.R - 2.5 * std::log(geom_group) - 2.5;
  }
  bracket = softmin_floor(bracket, BRACKET_FLOOR);
  out.f = 2.0 / (bracket * bracket);

  // e+ keeps the sign of the flow so a caller can tell direction; the power
  // law below uses the smoothed magnitude.
  out.e_plus = e_D * Re * std::sqrt(out.f / 2.0);
  const double ep_safe =
      std::sqrt(out.e_plus * out.e_plus + EPLUS_FLOOR * EPLUS_FLOOR);

  // The exponent actually applied to ep_safe in out.G -- NOT always
  // set.G_eplus_exponent. NarrowChannelAlphaSwitch's narrow sub-band scales
  // the base exponent by (W/H)^n_exponent, and the Re-derivative below must
  // differentiate THAT exponent, not the unscaled base, or it silently
  // disagrees with central differences whenever the narrow correction is
  // active (W/H and alpha are fixed geometry here, so this is still a
  // constant with respect to Re, same as every other shape).
  double g_eplus_exponent_effective = set.G_eplus_exponent;

  if (set.G_shape_model ==
      RibCorrelationSet::GShapeModel::NarrowChannelAlphaSwitch) {
    if (W_H < set.G_narrow_WH_band_boundary && W_H > 0.0 &&
        std::isfinite(W_H)) {
      g_eplus_exponent_effective *= std::pow(W_H, set.G_narrow_WH_n_exponent);
    }
    out.G = han_narrow_G(set, W_H, geom.alpha_deg, ep_safe);
  } else {
    double G_alpha_term, G_pe_term;
    if (set.G_shape_model ==
        RibCorrelationSet::GShapeModel::SquareVsRectangular) {
      const bool square = is_square_channel(W_H);
      const double m = square ? set.G_shape_alpha_exponent_square
                              : set.G_shape_alpha_exponent_rect;
      const double n = square ? set.G_shape_pe_exponent_square
                              : set.G_shape_pe_exponent_rect;
      G_alpha_term = (m != 0.0) ? std::pow(geom.alpha_deg / 90.0, m) : 1.0;
      G_pe_term = (n != 0.0 && geom.p_e / 10.0 > 0.0)
                      ? std::pow(geom.p_e / 10.0, n)
                      : 1.0;
    } else {
      G_alpha_term = power_term(
          set.G_alpha, geom.alpha_deg / 90.0 * set.G_alpha.reference);
      G_pe_term = power_term(set.G_pe, geom.p_e);
    }

    out.G = set.C_G * power_term(set.G_eD, e_D) * G_pe_term *
            power_term(set.G_WH, W_H) * G_alpha_term *
            std::pow(ep_safe, set.G_eplus_exponent);
  }

  // The source's own four-wall average, when printed. Same smoothed e+ as G.
  if (set.C_Gbar > 0.0) {
    out.G_bar = set.C_Gbar * std::pow(ep_safe, set.Gbar_eplus_exponent);
    out.has_G_bar = true;
  }

  const double root_f2 = std::sqrt(out.f / 2.0);
  const double denom =
      softmin_floor(1.0 + (out.G - out.R) * root_f2, ST_DENOM_FLOOR);
  out.St_r = out.f / (2.0 * denom);

  // dSt/dRe by the chain rule. f is independent of Re for this family, so only
  // G moves: dG/dRe = G * n * (d ep_safe/dRe) / ep_safe, and
  // d ep_safe/dRe = e_plus * (e_D * sqrt(f/2)) / ep_safe.
  const double dep_dRe = e_D * root_f2;
  const double dep_safe_dRe = out.e_plus * dep_dRe / ep_safe;
  const double dG_dRe =
      g_eplus_exponent_effective * out.G / ep_safe * dep_safe_dRe;
  out.dSt_dRe = -out.f * root_f2 * dG_dRe / (2.0 * denom * denom);

  out.extrapolated =
      outside(set.valid_Re, std::abs(Re)) || outside(set.valid_eD, e_D) ||
      outside(set.valid_pe, geom.p_e) || outside(set.valid_WH, W_H) ||
      outside(set.valid_alpha, geom.alpha_deg) ||
      outside(set.valid_eplus, std::abs(out.e_plus));

  return out;
}

namespace {

// An accuracy figure must carry a finite value exactly when it claims to
// have one. The old 0.0-means-unstated sentinel could not express that, and
// a zero read as "perfect agreement" rather than "nothing recorded".
void check_accuracy(const std::string &set_name, const char *field,
                    const StatedAccuracy &a) {
  const bool has_value = std::isfinite(a.value);
  const bool claims_value = a.provenance != AccuracyProvenance::Unstated;
  if (has_value != claims_value) {
    throw std::invalid_argument(
        "rib correlation set '" + set_name + "': " + std::string(field) +
        " must carry a finite value if and only if its provenance is not "
        "Unstated");
  }
  if (claims_value && !(a.value >= 0.0)) {
    throw std::invalid_argument("rib correlation set '" + set_name +
                                "': " + std::string(field) +
                                " must be non-negative");
  }
}

}  // namespace

void validate_rib_set(const RibCorrelationSet &set) {
  if (set.name.empty()) {
    throw std::invalid_argument("rib correlation set: name must not be empty");
  }
  check_accuracy(set.name, "accuracy_R", set.accuracy_R);
  check_accuracy(set.name, "accuracy_G", set.accuracy_G);
  // C_Gbar = 0 means "not printed"; anything else must be a usable fit.
  if (!(set.C_Gbar >= 0.0) || !std::isfinite(set.C_Gbar) ||
      !std::isfinite(set.Gbar_eplus_exponent)) {
    throw std::invalid_argument(
        "rib correlation set '" + set.name +
        "': C_Gbar must be finite and non-negative (0 = not printed), and "
        "Gbar_eplus_exponent finite");
  }
  if (set.source.empty()) {
    throw std::invalid_argument(
        "rib correlation set '" + set.name +
        "': source must not be empty. A coefficient whose origin is not "
        "recorded is what this rebuild exists to remove.");
  }
  if (set.R_alpha_shape == RibCorrelationSet::RAlphaShape::QuadraticAlpha) {
    // C_R plays no part in this shape -- the quadratic fields replace it
    // entirely -- so it is exempt from the positivity check below rather
    // than requiring a caller to set a dummy value nothing reads.
    if (!std::isfinite(set.R_quad_c0) || !std::isfinite(set.R_quad_c1) ||
        !std::isfinite(set.R_quad_c2)) {
      throw std::invalid_argument("rib correlation set '" + set.name +
                                  "': R_quad_c0/c1/c2 must be finite");
    }
    if (!std::isfinite(set.R_quad_WH_exponent_at_90) ||
        !std::isfinite(set.R_quad_WH_exponent_off_90)) {
      throw std::invalid_argument(
          "rib correlation set '" + set.name +
          "': R_quad_WH_exponent_at_90/off_90 must be finite");
    }
    if (!(set.R_quad_WH_cap >= 0.0) || !std::isfinite(set.R_quad_WH_cap)) {
      throw std::invalid_argument(
          "rib correlation set '" + set.name +
          "': R_quad_WH_cap must be finite and non-negative (0 = uncapped)");
    }
  } else if (set.R_alpha_shape ==
             RibCorrelationSet::RAlphaShape::QuadraticAlphaTwoBand) {
    // Same exemption as QuadraticAlpha above: C_R plays no part here either.
    if (!std::isfinite(set.R_quad_c0) || !std::isfinite(set.R_quad_c1) ||
        !std::isfinite(set.R_quad_c2) ||
        !std::isfinite(set.R_quad_narrow_c0) ||
        !std::isfinite(set.R_quad_narrow_c1) ||
        !std::isfinite(set.R_quad_narrow_c2)) {
      throw std::invalid_argument(
          "rib correlation set '" + set.name +
          "': R_quad_c0/c1/c2 and R_quad_narrow_c0/c1/c2 must be finite");
    }
    if (!(set.R_WH_band_boundary > 0.0) ||
        !std::isfinite(set.R_WH_band_boundary)) {
      throw std::invalid_argument(
          "rib correlation set '" + set.name +
          "': R_WH_band_boundary must be finite and positive, got " +
          std::to_string(set.R_WH_band_boundary));
    }
  } else if (!(set.C_R > 0.0) || !std::isfinite(set.C_R)) {
    throw std::invalid_argument("rib correlation set '" + set.name +
                                "': C_R must be finite and positive, got " +
                                std::to_string(set.C_R));
  }
  if (set.G_shape_model !=
      RibCorrelationSet::GShapeModel::NarrowChannelAlphaSwitch) {
    // C_G plays no part in NarrowChannelAlphaSwitch -- G_narrow_C_alpha90/
    // off_axis replace it entirely -- so it is exempt here, same reasoning
    // as C_R's exemption for the quadratic R shapes above.
    if (!(set.C_G > 0.0) || !std::isfinite(set.C_G)) {
      throw std::invalid_argument(
          "rib correlation set '" + set.name +
          "': C_G must be finite and positive, got " +
          std::to_string(set.C_G));
    }
  }
  if (!std::isfinite(set.G_eplus_exponent)) {
    throw std::invalid_argument("rib correlation set '" + set.name +
                                "': G_eplus_exponent must be finite");
  }
  if (set.G_shape_model ==
      RibCorrelationSet::GShapeModel::NarrowChannelAlphaSwitch) {
    if (!std::isfinite(set.G_narrow_C_alpha90) ||
        !std::isfinite(set.G_narrow_C_off_axis)) {
      throw std::invalid_argument(
          "rib correlation set '" + set.name +
          "': G_narrow_C_alpha90/off_axis must be finite");
    }
    if (!(set.G_narrow_WH_band_boundary > 0.0) ||
        !std::isfinite(set.G_narrow_WH_band_boundary)) {
      throw std::invalid_argument(
          "rib correlation set '" + set.name +
          "': G_narrow_WH_band_boundary must be finite and positive, got " +
          std::to_string(set.G_narrow_WH_band_boundary));
    }
    if (!std::isfinite(set.G_narrow_WH_C_exponent) ||
        !std::isfinite(set.G_narrow_WH_n_exponent)) {
      throw std::invalid_argument(
          "rib correlation set '" + set.name +
          "': G_narrow_WH_C_exponent/n_exponent must be finite");
    }
  }
  if (set.G_shape_model == RibCorrelationSet::GShapeModel::SquareVsRectangular) {
    if (!std::isfinite(set.G_shape_alpha_exponent_square) ||
        !std::isfinite(set.G_shape_alpha_exponent_rect) ||
        !std::isfinite(set.G_shape_pe_exponent_square) ||
        !std::isfinite(set.G_shape_pe_exponent_rect)) {
      throw std::invalid_argument(
          "rib correlation set '" + set.name +
          "': G_shape_* exponents must be finite");
    }
  }
  require_positive_reference(set.R_eD, "R_eD");
  require_positive_reference(set.R_pe, "R_pe");
  require_positive_reference(set.R_WH, "R_WH");
  require_positive_reference(set.R_alpha, "R_alpha");
  require_positive_reference(set.G_eD, "G_eD");
  require_positive_reference(set.G_pe, "G_pe");
  require_positive_reference(set.G_WH, "G_WH");
  require_positive_reference(set.G_alpha, "G_alpha");
}

double rib_ribbed_wall_fraction(double W_H, int n_ribbed_walls) {
  if (!(W_H > 0.0) || !std::isfinite(W_H)) {
    throw std::invalid_argument(
        "rib_ribbed_wall_fraction: W_H must be positive and finite");
  }
  // Ribs sit on the walls of width W; the wetted perimeter is 2(W + H).
  const double width_share = W_H / (W_H + 1.0);  // W / (W + H)
  switch (n_ribbed_walls) {
  case 1:
    return 0.5 * width_share;
  case 2:
    return width_share;
  case 4:
    return 1.0;
  default:
    throw std::invalid_argument(
        "n_ribbed_walls must be 1, 2 or 4, got " +
        std::to_string(n_ribbed_walls) +
        ". 2 means two opposite walls, which is the configuration the "
        "correlations were measured on.");
  }
}

std::tuple<double, double> rib_channel_friction(const RibCorrelationSet &set,
                                                const RibGeometry &geom,
                                                double Re, int n_ribbed_walls) {
  const double phi = rib_ribbed_wall_fraction(geom.W_H, n_ribbed_walls);
  const double f_four_sided = evaluate_rib(set, geom, Re).f;

  // f_s at a smoothly floored |Re|; d|Re|_safe/dRe = Re / |Re|_safe.
  const double re_safe = std::sqrt(Re * Re + SMOOTH_RE_FLOOR * SMOOTH_RE_FLOOR);
  const double f_s =
      HAN_SMOOTH_FANNING_COEF * std::pow(re_safe, HAN_SMOOTH_FANNING_RE_EXP);
  const double df_s_dRe =
      HAN_SMOOTH_FANNING_RE_EXP * f_s / re_safe * (Re / re_safe);

  const double smooth_share = 1.0 - phi;
  const double f_fanning = phi * f_four_sided + smooth_share * f_s;
  return {DARCY_PER_FANNING * f_fanning,
          DARCY_PER_FANNING * smooth_share * df_s_dRe};
}

}  // namespace cooling
}  // namespace combaero
