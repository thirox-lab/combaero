#include "pin_fin_correlation.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "input_validation.h"
#include "math_constants.h"  // MSVC compatibility for M_PI

namespace combaero {
namespace cooling {

namespace {

// Smooth one-sided floor for Pr, as in rib_ratio_correlation.cpp.
constexpr double PR_EPS = 1e-4;

double smooth_pr(double Pr) {
  return 0.5 * (Pr + std::sqrt(Pr * Pr + PR_EPS * PR_EPS));
}

// (value / reference)^exponent, 1 for a zero exponent or a bad base;
// geometry terms are fixed per evaluation and need no derivative.
double geometry_term(const RibTerm &t, double value) {
  if (t.exponent == 0.0) {
    return 1.0;
  }
  const double x = value / t.reference;
  if (!(x > 0.0) || !std::isfinite(x)) {
    return 1.0;
  }
  return std::pow(x, t.exponent);
}

struct ValueSlope {
  double v = 0.0;
  double d = 0.0;  // derivative with respect to the smoothed |Re|
};

// C * x^a and its x-derivative; x > 0 always.
ValueSlope power_law(double C, double a, double x) {
  const double v = C * (a == 0.0 ? 1.0 : std::pow(x, a));
  return {v, a * v / x};
}

bool outside(const RibRange &r, double v) {
  if (!r.bounded()) {
    return r.lo > 0.0 && v < r.lo * (1.0 - RibRange::rel_tol);
  }
  return v < r.lo * (1.0 - RibRange::rel_tol) || v > r.hi * (1.0 + RibRange::rel_tol);
}

// C1 smoothstep weight in ln x over [x_lo, x_hi].
ValueSlope blend_weight(double x, double x_lo, double x_hi) {
  if (x >= x_hi) {
    return {1.0, 0.0};
  }
  if (x <= x_lo) {
    return {0.0, 0.0};
  }
  const double span = std::log(x_hi / x_lo);
  const double s = std::log(x / x_lo) / span;
  return {s * s * (3.0 - 2.0 * s), 6.0 * s * (1.0 - s) / (span * x)};
}

bool geometry_ok(const PinFinGeometry &g) {
  return std::isfinite(g.S_D) && std::isfinite(g.X_D) &&
         std::isfinite(g.H_D) && g.S_D > 1.0 && g.X_D > 0.0 && g.H_D > 0.0 &&
         g.N_rows >= 1 && pin_fin_min_gap_D(g) > 0.0;
}

bool geometry_outside(const RibRange &S, const RibRange &X, const RibRange &H,
                      const PinFinGeometry &g) {
  return outside(S, g.S_D) || outside(X, g.X_D) || outside(H, g.H_D);
}

void fail(const std::string &what, const std::string &name,
          const std::string &msg) {
  throw std::invalid_argument(what + " '" + name + "': " + msg);
}

void require_term(const std::string &what, const std::string &name,
                  const RibTerm &t, const char *field) {
  if (!(t.reference > 0.0) || !std::isfinite(t.reference) ||
      !std::isfinite(t.exponent)) {
    fail(what, name,
         std::string(field) +
             " needs a finite positive reference and a finite exponent");
  }
}

void check_accuracy(const std::string &what, const std::string &name,
                    const char *field, const StatedAccuracy &a) {
  const bool has_value = std::isfinite(a.value);
  const bool claims_value = a.provenance != AccuracyProvenance::Unstated;
  if (has_value != claims_value || (claims_value && !(a.value >= 0.0))) {
    fail(what, name,
         std::string(field) +
             " must carry a finite, non-negative value if and only if its "
             "provenance is not Unstated");
  }
}

RibRange range(double lo, double hi) {
  RibRange r;
  r.lo = lo;
  r.hi = hi;
  return r;
}

RibTerm term(double exponent, double reference = 1.0) {
  RibTerm t;
  t.exponent = exponent;
  t.reference = reference;
  return t;
}

}  // namespace

// ---- Geometry ----

double pin_fin_min_gap_D(const PinFinGeometry &g) {
  const double transverse = g.S_D - 1.0;
  if (g.arrangement == PinArrangement::Inline) {
    return transverse;
  }
  const double half = 0.5 * g.S_D;
  const double diagonal = 2.0 * (std::sqrt(half * half + g.X_D * g.X_D) - 1.0);
  return std::min(transverse, diagonal);
}

double pin_fin_amin_over_afrontal(const PinFinGeometry &g) {
  return pin_fin_min_gap_D(g) / g.S_D;
}

double pin_fin_dprime_over_D(const PinFinGeometry &g) {
  const double xs = g.X_D * g.S_D;
  return g.H_D * (4.0 * xs - M_PI) / (2.0 * xs + M_PI * (g.H_D - 0.5));
}

double pin_fin_aprime_over_amin(const PinFinGeometry &g) {
  return (g.S_D - M_PI / (4.0 * g.X_D)) / pin_fin_min_gap_D(g);
}

double pin_fin_dh_over_D(const PinFinGeometry &g) {
  const double xs = g.X_D * g.S_D;
  return 4.0 * g.X_D * g.H_D * pin_fin_min_gap_D(g) /
         (2.0 * xs + M_PI * (g.H_D - 0.5));
}

PinFinAreaFractions pin_fin_area_fractions(const PinFinGeometry &g) {
  PinFinAreaFractions out;
  const double base = g.X_D * g.S_D;
  const double endwall = base - 0.25 * M_PI;
  const double pin = 0.5 * M_PI * g.H_D;
  out.endwall_exposed = endwall / base;
  out.pin = pin / base;
  out.pin_over_total = pin / (endwall + pin);
  return out;
}

PinFinEfficiency pin_fin_array_efficiency(double h, double k_pin, double D,
                                          double H, double A_f_over_A_t) {
  PinFinEfficiency out;
  // z = m L with L = H / 2 and m = sqrt(4 h / (k D)), so z^2 = H^2 h / (k D).
  // d eta_fin / dh = (z sech^2 z - tanh z) H^2 / (2 k D z^3), whose h -> 0
  // limit is -H^2 / (3 k D).
  // NaN used to reach the "no fin model" branch below and answer eta = 1
  // (#769). +inf k_pin (an isothermal pin) and +inf h are limits and stay.
  const auto finite_or_pos_inf = [](double x) {
    return std::isfinite(x) || (std::isinf(x) && x > 0.0);
  };
  if (!finite_or_pos_inf(h) || !finite_or_pos_inf(k_pin) || !validation::is_finite(D) ||
      !validation::is_finite(H) || !validation::is_finite(A_f_over_A_t)) {
    throw std::invalid_argument(
        "pin_fin_array_efficiency: h and k_pin must be finite or +inf; D, H and "
        "A_f_over_A_t finite");
  }
  if (!(k_pin > 0.0) || !(D > 0.0) || !std::isfinite(k_pin) ||
      !(H >= 0.0)) {
    return out;  // no fin model: eta = 1
  }
  const double scale = H * H / (k_pin * D);  // z^2 / h
  const double hp = std::max(h, 0.0);
  const double z = std::sqrt(scale * hp);
  double eta = 1.0;
  double g = -2.0 / 3.0;  // (z sech^2 z - tanh z) / z^3
  if (z < 2e-2) {
    // Series: the closed form cancels catastrophically as z -> 0.
    const double z2 = z * z;
    eta = 1.0 - z2 / 3.0 + 2.0 * z2 * z2 / 15.0 -
          17.0 * z2 * z2 * z2 / 315.0;
    g = -2.0 / 3.0 + 8.0 * z2 / 15.0 - 34.0 * z2 * z2 / 105.0;
  } else {
    const double t = std::tanh(z);
    const double sech2 = 1.0 - t * t;
    eta = t / z;
    g = (z * sech2 - t) / (z * z * z);
  }
  out.eta_fin = eta;
  out.deta_fin_dh = 0.5 * g * scale;
  out.eta_t = 1.0 - A_f_over_A_t * (1.0 - eta);
  out.deta_t_dh = A_f_over_A_t * out.deta_fin_dh;
  return out;
}

// ---- Evaluators ----

PinFinNuResult evaluate_pin_fin_nu(const PinFinNuSet &set,
                                   const PinFinGeometry &geom, double Re_D,
                                   double Pr) {
  PinFinNuResult out;
  const bool geo_ok = geometry_ok(geom);

  // Re_native = s Re_D and Nu_D = Nu_native / r: both exact and linear.
  double s = 1.0;
  double r = 1.0;
  if (set.re_basis == PinReBasis::VanFossenDprime && geo_ok) {
    r = pin_fin_dprime_over_D(geom);
    s = r / pin_fin_aprime_over_amin(geom);
  }

  const double x = std::sqrt(Re_D * Re_D + PIN_FIN_RE_EPS * PIN_FIN_RE_EPS);
  const double dx_dRe = Re_D / x;
  const double xn = s * x;
  const double g = geometry_term(set.term_XD, geom.X_D) *
                   geometry_term(set.term_SD, geom.S_D) *
                   geometry_term(set.term_HD, geom.H_D);
  const double pr = set.Pr_exp == 0.0 ? 1.0
                                      : std::pow(smooth_pr(Pr), set.Pr_exp);
  const ValueSlope N = power_law(set.C * g * pr, set.Re_exp, xn);

  out.Nu = N.v / r;
  out.dNu_dRe = N.d * s / r * dx_dRe;
  // Nu is linear in smooth_pr(Pr)^Pr_exp: dNu/dPr = Nu Pr_exp sp'(Pr) / sp(Pr).
  if (set.Pr_exp != 0.0) {
    const double sp = smooth_pr(Pr);
    const double dsp = 0.5 * (1.0 + Pr / std::sqrt(Pr * Pr + PR_EPS * PR_EPS));
    out.dNu_dPr = sp > 0.0 ? out.Nu * set.Pr_exp * dsp / sp : 0.0;
  }
  out.Re_native = s * Re_D;
  out.extrapolated =
      !geo_ok || geom.arrangement != set.arrangement ||
      outside(set.valid_Re, xn) ||
      geometry_outside(set.valid_SD, set.valid_XD, set.valid_HD, geom) ||
      outside(set.valid_Nrows, static_cast<double>(geom.N_rows));
  return out;
}

PinFinFrictionResult evaluate_pin_fin_friction(const PinFinFrictionSet &set,
                                               const PinFinGeometry &geom,
                                               double Re_D) {
  PinFinFrictionResult out;
  const bool geo_ok = geometry_ok(geom);
  const double x = std::sqrt(Re_D * Re_D + PIN_FIN_RE_EPS * PIN_FIN_RE_EPS);
  const double dx_dRe = Re_D / x;
  const double g = geometry_term(set.term_SD, geom.S_D) *
                   geometry_term(set.term_XD, geom.X_D) *
                   geometry_term(set.term_HD, geom.H_D);

  ValueSlope f = power_law(set.C1 * g, set.Re_exp1, x);
  if (set.Re_split > 0.0) {
    const ValueSlope f2 = power_law(set.C2 * g, set.Re_exp2, x);
    const ValueSlope w = blend_weight(x, set.Re_split / PIN_FIN_F_BLEND,
                                      set.Re_split * PIN_FIN_F_BLEND);
    f = {(1.0 - w.v) * f.v + w.v * f2.v,
         w.d * (f2.v - f.v) + (1.0 - w.v) * f.d + w.v * f2.d};
  }

  // Native -> canonical per-row basis.
  double to_canonical = 1.0;
  if (set.basis == PinFrictionBasis::PerRowGapVmax && geom.N_rows >= 1) {
    to_canonical = static_cast<double>(geom.N_rows - 1) /
                   static_cast<double>(geom.N_rows);
  } else if (set.basis == PinFrictionBasis::PerRowHalfHeadVmax) {
    to_canonical = 0.25;
  }
  out.f = f.v * to_canonical;
  out.df_dRe = f.d * to_canonical * dx_dRe;
  out.extrapolated =
      !geo_ok || geom.arrangement != set.arrangement ||
      outside(set.valid_Re, x) ||
      geometry_outside(set.valid_SD, set.valid_XD, set.valid_HD, geom) ||
      outside(set.valid_Nrows, static_cast<double>(geom.N_rows));
  return out;
}

PinFinModifierResult evaluate_pin_fin_modifier(const PinFinRatioModifier &mod,
                                               const PinFinGeometry &geom,
                                               double Re_D) {
  PinFinModifierResult out;
  const double x = std::sqrt(Re_D * Re_D + PIN_FIN_RE_EPS * PIN_FIN_RE_EPS);
  const double dx_dRe = Re_D / x;
  const ValueSlope rN = power_law(mod.C_Nu, mod.Nu_Re_exp, x);
  out.ratio_Nu = rN.v;
  out.dratio_Nu_dRe = rN.d * dx_dRe;
  out.has_f = mod.has_f;
  if (mod.has_f) {
    const ValueSlope rF = power_law(mod.C_f, mod.f_Re_exp, x);
    out.ratio_f = rF.v;
    out.dratio_f_dRe = rF.d * dx_dRe;
  }
  out.extrapolated =
      !geometry_ok(geom) || geom.arrangement != mod.to_arrangement ||
      outside(mod.valid_Re, x) ||
      geometry_outside(mod.valid_SD, mod.valid_XD, mod.valid_HD, geom);
  return out;
}

// ---- Validation ----

void validate_pin_fin_geometry(const PinFinGeometry &geom) {
  if (!(geom.S_D > 1.0) || !std::isfinite(geom.S_D)) {
    fail("pin-fin geometry", "", "S_D must be finite and > 1 (pins overlap)");
  }
  if (!(geom.X_D > 0.0) || !std::isfinite(geom.X_D)) {
    fail("pin-fin geometry", "", "X_D must be finite and positive");
  }
  if (!(geom.H_D > 0.0) || !std::isfinite(geom.H_D)) {
    fail("pin-fin geometry", "", "H_D must be finite and positive");
  }
  if (geom.N_rows < 1) {
    fail("pin-fin geometry", "", "N_rows must be at least 1");
  }
  if (!(pin_fin_min_gap_D(geom) > 0.0)) {
    fail("pin-fin geometry", "",
         "the diagonal gap is closed: adjacent rows' pins overlap");
  }
}

void validate_pin_fin_nu_set(const PinFinNuSet &set) {
  const std::string w = "pin-fin Nu set";
  if (set.name.empty()) {
    fail(w, set.name, "a set needs a name");
  }
  if (!(set.C > 0.0) || !std::isfinite(set.C)) {
    fail(w, set.name, "C must be finite and positive");
  }
  if (!std::isfinite(set.Re_exp) || !std::isfinite(set.Pr_exp)) {
    fail(w, set.name, "Re_exp and Pr_exp must be finite");
  }
  require_term(w, set.name, set.term_XD, "term_XD");
  require_term(w, set.name, set.term_SD, "term_SD");
  require_term(w, set.name, set.term_HD, "term_HD");
  check_accuracy(w, set.name, "accuracy_Nu", set.accuracy_Nu);
}

void validate_pin_fin_friction_set(const PinFinFrictionSet &set) {
  const std::string w = "pin-fin friction set";
  if (set.name.empty()) {
    fail(w, set.name, "a set needs a name");
  }
  if (!(set.C1 > 0.0) || !std::isfinite(set.C1) ||
      !std::isfinite(set.Re_exp1)) {
    fail(w, set.name, "C1 must be finite and positive, Re_exp1 finite");
  }
  if (!(set.Re_split >= 0.0) || !std::isfinite(set.Re_split)) {
    fail(w, set.name, "Re_split must be finite and >= 0 (0 = one segment)");
  }
  if (set.Re_split > 0.0 &&
      (!(set.C2 > 0.0) || !std::isfinite(set.C2) ||
       !std::isfinite(set.Re_exp2))) {
    fail(w, set.name,
         "a two-segment set needs C2 finite and positive, Re_exp2 finite");
  }
  require_term(w, set.name, set.term_SD, "term_SD");
  require_term(w, set.name, set.term_XD, "term_XD");
  require_term(w, set.name, set.term_HD, "term_HD");
  check_accuracy(w, set.name, "accuracy_f", set.accuracy_f);
}

void validate_pin_fin_modifier(const PinFinRatioModifier &mod) {
  const std::string w = "pin-fin ratio modifier";
  if (mod.name.empty()) {
    fail(w, mod.name, "a modifier needs a name");
  }
  if (!(mod.C_Nu > 0.0) || !std::isfinite(mod.C_Nu) ||
      !std::isfinite(mod.Nu_Re_exp)) {
    fail(w, mod.name, "C_Nu must be finite and positive, Nu_Re_exp finite");
  }
  if (mod.has_f && (!(mod.C_f > 0.0) || !std::isfinite(mod.C_f) ||
                    !std::isfinite(mod.f_Re_exp))) {
    fail(w, mod.name, "C_f must be finite and positive, f_Re_exp finite");
  }
}

// ---- Shipped sets ----

PinFinNuSet metzger_1986_staggered_nu() {
  PinFinNuSet s;
  s.name = "metzger_1986_staggered_nu";
  s.source =
      "Metzger, Shepard and Haley (1986), ASME 86-GT-132, via Armstrong and "
      "Winstanley (1988), J. Turbomach. 110, 94, Eq. 2";
  s.validity_source =
      "Armstrong and Winstanley (1988) recommended short-pin limits; fitted "
      "at H/D 1, S/D 2.5, 1.5 <= X/D <= 5, 10 rows";
  s.provenance = RibProvenance::Extracted;
  s.arrangement = PinArrangement::Staggered;
  s.surface = PinNuSurface::Total;
  s.re_basis = PinReBasis::DiameterVmax;
  s.C = 0.135;
  s.Re_exp = 0.69;
  s.term_XD = term(-0.34);
  s.valid_Re = range(1.0e3, 1.0e5);
  s.valid_SD = range(2.0, 4.0);
  s.valid_XD = range(1.5, 5.0);
  s.valid_HD = range(0.5, 3.0);
  s.valid_Nrows = range(9.5, 10.5);
  s.valid_Pr = 0.7;
  s.accuracy_Nu = StatedAccuracy::stated(0.20);
  return s;
}

PinFinFrictionSet metzger_1982_staggered_friction() {
  PinFinFrictionSet s;
  s.name = "metzger_1982_staggered_friction";
  s.source =
      "Metzger, Fan and Shepard (1982), Heat Transfer 1982 vol. 3, 137-142, "
      "via Armstrong and Winstanley (1988), Eqs 20-21";
  s.validity_source =
      "Armstrong and Winstanley (1988): 0.5 <= H/D <= 6, 2 <= S/D <= 4 "
      "(Peng's data); fitted at H/D 1, S/D 2.5, 1.5 <= X/D <= 5";
  s.provenance = RibProvenance::Extracted;
  s.arrangement = PinArrangement::Staggered;
  s.basis = PinFrictionBasis::PerRowVmax;
  s.C1 = 0.317;
  s.Re_exp1 = -0.132;
  s.C2 = 1.76;
  s.Re_exp2 = -0.318;
  s.Re_split = 1.0e4;
  s.valid_Re = range(1.0e3, 1.0e5);
  s.valid_SD = range(2.0, 4.0);
  s.valid_XD = range(1.5, 5.0);
  s.valid_HD = range(0.5, 6.0);
  s.accuracy_f = StatedAccuracy::stated(0.15);
  return s;
}

PinFinNuSet vanfossen_1982_staggered_nu() {
  PinFinNuSet s;
  s.name = "vanfossen_1982_staggered_nu";
  s.source =
      "VanFossen (1982), J. Eng. Power 104, 268 (NASA TM-81696), Eq. 16";
  s.validity_source =
      "NASA TM-81696 Table I: equilateral arrays, H/D 0.5 and 2, S/D 2 and "
      "4, 4 rows, Re_D' 300 to 60,000";
  s.provenance = RibProvenance::Extracted;
  s.arrangement = PinArrangement::Staggered;
  s.surface = PinNuSurface::Total;
  s.re_basis = PinReBasis::VanFossenDprime;
  s.C = 0.153;
  s.Re_exp = 0.685;
  s.valid_Re = range(300.0, 6.0e4);
  s.valid_SD = range(2.0, 4.0);
  s.valid_XD = range(1.73, 3.47);
  s.valid_HD = range(0.5, 2.0);
  s.valid_Nrows = range(3.5, 4.5);
  s.valid_Pr = 0.7;
  return s;
}

PinFinFrictionSet damerow_1972_staggered_friction() {
  PinFinFrictionSet s;
  s.name = "damerow_1972_staggered_friction";
  s.source =
      "Damerow, Murtaugh and Burggraf (1972), NASA CR-120883, Eqs 10 and 18";
  s.validity_source =
      "CR-120883 Fig. 3: square-diagonal staggered, X_T/D 4.24 and 7.07, "
      "X_L/D 2.12 and 3.54, Z/D 2 to 4, 10 rows; Re 2e3 to 6e4 (Figs 27-28); "
      "f rises above inlet Mach 0.36";
  s.provenance = RibProvenance::Extracted;
  s.arrangement = PinArrangement::Staggered;
  s.basis = PinFrictionBasis::PerRowGapVmax;
  s.C1 = 2.06;
  s.Re_exp1 = -0.16;
  s.term_SD = term(-1.1);
  s.valid_Re = range(2.0e3, 6.0e4);
  s.valid_SD = range(4.24, 7.07);
  s.valid_XD = range(2.12, 3.54);
  s.valid_HD = range(2.0, 4.0);
  s.valid_Nrows = range(9.5, 10.5);
  return s;
}

PinFinNuSet chyu_1998_nu(PinArrangement arrangement, PinNuSurface surface) {
  // Han, Dutta and Ekkad (2012) Table 4.7.
  struct Row {
    PinArrangement arr;
    PinNuSurface surf;
    double a, b;
    const char *tag;
  };
  static const Row rows[] = {
      {PinArrangement::Inline, PinNuSurface::Pin, 0.155, 0.658, "inline_pin"},
      {PinArrangement::Inline, PinNuSurface::Endwall, 0.052, 0.759,
       "inline_endwall"},
      {PinArrangement::Inline, PinNuSurface::Total, 0.068, 0.733,
       "inline_total"},
      {PinArrangement::Staggered, PinNuSurface::Pin, 0.337, 0.585,
       "staggered_pin"},
      {PinArrangement::Staggered, PinNuSurface::Endwall, 0.315, 0.582,
       "staggered_endwall"},
      {PinArrangement::Staggered, PinNuSurface::Total, 0.320, 0.583,
       "staggered_total"},
  };
  for (const Row &r : rows) {
    if (r.arr != arrangement || r.surf != surface) {
      continue;
    }
    PinFinNuSet s;
    s.name = std::string("chyu_1998_") + r.tag;
    s.source =
        "Chyu, Hsing, Shih and Natarajan (1998), ASME 98-GT-175, via Han, "
        "Dutta and Ekkad (2012) Table 4.7";
    s.validity_source =
        "S/D = X/D = 2.5, H/D = 1 (Lyall 2006, Table 2-1); 7 rows; Re span "
        "of Han Fig. 4.130";
    s.provenance = RibProvenance::Extracted;
    s.arrangement = r.arr;
    s.surface = r.surf;
    s.re_basis = PinReBasis::DiameterVmax;
    s.C = r.a;
    s.Re_exp = r.b;
    s.Pr_exp = 0.4;
    s.valid_Re = range(5.0e3, 2.5e4);
    s.valid_SD = range(2.45, 2.55);
    s.valid_XD = range(2.45, 2.55);
    s.valid_HD = range(0.95, 1.05);
    s.valid_Nrows = range(6.5, 7.5);
    s.valid_Pr = 0.7;
    return s;
  }
  throw std::invalid_argument("chyu_1998_nu: unknown arrangement/surface");
}

PinFinRatioModifier chyu_1998_inline_over_staggered() {
  PinFinRatioModifier m;
  m.name = "chyu_1998_inline_over_staggered";
  m.source =
      "Chyu, Hsing, Shih and Natarajan (1998), ASME 98-GT-175: the ratio of "
      "the inline and staggered Total correlations of Han Table 4.7";
  m.provenance = RibProvenance::Extracted;
  m.from_arrangement = PinArrangement::Staggered;
  m.to_arrangement = PinArrangement::Inline;
  m.surface = PinNuSurface::Total;
  m.C_Nu = 0.068 / 0.320;
  m.Nu_Re_exp = 0.733 - 0.583;
  m.has_f = false;
  m.valid_Re = range(5.0e3, 2.5e4);
  m.valid_SD = range(2.45, 2.55);
  m.valid_XD = range(2.45, 2.55);
  m.valid_HD = range(0.95, 1.05);
  return m;
}

PinFinNuSet chyu_1990_nu(PinArrangement arrangement, bool fillet) {
  // Chyu (1990) Table 2.
  double A = 0.0, B = 0.0;
  if (arrangement == PinArrangement::Inline) {
    A = fillet ? 0.403 : 0.463;
    B = fillet ? 0.550 : 0.537;
  } else {
    A = fillet ? 0.234 : 0.690;
    B = fillet ? 0.608 : 0.511;
  }
  PinFinNuSet s;
  s.name = std::string("chyu_1990_") +
           (arrangement == PinArrangement::Inline ? "inline" : "staggered") +
           (fillet ? "_fillet" : "_straight");
  s.source =
      "Chyu (1990), J. Heat Transfer 112, 926, Table 2 (naphthalene, pin "
      "surface; endwall/pin Sh 0.89-1.09 per Table 1)";
  s.validity_source =
      "Chyu (1990): H/D = 1, S/D = X/D = 2.5, 7 rows, Re 5e3 to 3e4";
  s.provenance = RibProvenance::Extracted;
  s.arrangement = arrangement;
  s.surface = PinNuSurface::Pin;
  s.re_basis = PinReBasis::DiameterVmax;
  s.C = A;
  s.Re_exp = B;
  s.Pr_exp = 0.4;
  s.valid_Re = range(5.0e3, 3.0e4);
  s.valid_SD = range(2.45, 2.55);
  s.valid_XD = range(2.45, 2.55);
  s.valid_HD = range(0.95, 1.05);
  s.valid_Nrows = range(6.5, 7.5);
  s.valid_Pr = 0.7;
  return s;
}

PinFinFrictionSet chyu_1990_friction(PinArrangement arrangement, bool fillet) {
  // Fitted to Chyu (1990) Fig. 6 as digitised in
  // validation/cooling/data/chyu1990/; Re ranges are the digitised extremes,
  // rounded outward.
  struct Fit {
    double C, m, re_lo, re_hi, rms;
  };
  Fit fit{};
  if (arrangement == PinArrangement::Inline) {
    fit = fillet ? Fit{0.6851, -0.1388, 1.02e4, 2.21e4, 0.0021}
                 : Fit{0.1693, 0.0, 8.8e3, 2.08e4, 0.0137};
  } else {
    fit = fillet ? Fit{9.3866, -0.3632, 1.04e4, 2.24e4, 0.0188}
                 : Fit{1.6163, -0.1867, 9.2e3, 2.12e4, 0.0202};
  }
  PinFinFrictionSet s;
  s.name = std::string("chyu_1990_") +
           (arrangement == PinArrangement::Inline ? "inline" : "staggered") +
           (fillet ? "_fillet" : "_straight") + "_friction";
  s.source =
      "Chyu (1990), J. Heat Transfer 112, 926, Fig. 6, digitised and fitted "
      "(validation/cooling/data/chyu1990)";
  s.validity_source =
      "Chyu (1990): H/D = 1, S/D = X/D = 2.5, 7 rows; Re span of the "
      "digitised points";
  s.provenance = RibProvenance::Fitted;
  s.arrangement = arrangement;
  s.basis = PinFrictionBasis::PerRowHalfHeadVmax;
  s.C1 = fit.C;
  s.Re_exp1 = fit.m;
  s.valid_Re = range(fit.re_lo, fit.re_hi);
  s.valid_SD = range(2.45, 2.55);
  s.valid_XD = range(2.45, 2.55);
  s.valid_HD = range(0.95, 1.05);
  s.valid_Nrows = range(6.5, 7.5);
  // The fit's own residual against the points it was fitted to: reportable,
  // never a judging band.
  s.accuracy_f = StatedAccuracy::measured(fit.rms);
  return s;
}

PinFinRatioModifier chyu_1990_fillet_over_straight(
    PinArrangement arrangement) {
  const PinFinNuSet straight = chyu_1990_nu(arrangement, false);
  const PinFinNuSet fillet = chyu_1990_nu(arrangement, true);
  PinFinRatioModifier m;
  m.name = std::string("chyu_1990_fillet_over_straight_") +
           (arrangement == PinArrangement::Inline ? "inline" : "staggered");
  m.source =
      "Chyu (1990), J. Heat Transfer 112, 926: the ratio of the fillet and "
      "straight-pin correlations of Table 2";
  m.provenance = RibProvenance::Extracted;
  m.from_arrangement = arrangement;
  m.to_arrangement = arrangement;
  m.surface = PinNuSurface::Pin;
  m.C_Nu = fillet.C / straight.C;
  m.Nu_Re_exp = fillet.Re_exp - straight.Re_exp;
  m.has_f = false;
  m.valid_Re = range(5.0e3, 3.0e4);
  m.valid_SD = range(2.45, 2.55);
  m.valid_XD = range(2.45, 2.55);
  m.valid_HD = range(0.95, 1.05);
  return m;
}

}  // namespace cooling
}  // namespace combaero
