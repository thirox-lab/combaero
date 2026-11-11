#include "rib_ratio_correlation.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "heat_transfer.h"

namespace combaero {
namespace cooling {

namespace {

// Smooth one-sided floor for Pr: equals Pr well above PR_EPS, tends to it
// below, differentiable throughout. Pr is a fluid property and positive in
// practice; the guard only keeps Pr^n finite for a state no fluid reaches.
constexpr double PR_EPS = 1e-4;

double smooth_pr(double Pr) {
  return 0.5 * (Pr + std::sqrt(Pr * Pr + PR_EPS * PR_EPS));
}

// (value / reference)^exponent, 1 for a zero exponent; geometry terms only,
// whose values are fixed per evaluation and need no derivative.
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

// C * geometry * (x / Re_ref)^a, and its x-derivative. x > 0 always.
ValueSlope ratio_law(double C, double geometry, const RibTerm &re_term,
                     double x) {
  const double v =
      C * geometry *
      (re_term.exponent == 0.0
           ? 1.0
           : std::pow(x / re_term.reference, re_term.exponent));
  return {v, re_term.exponent * v / x};
}

ValueSlope baseline_law(const RatioBaseline &b, double x, double pr) {
  const double v =
      b.coeff * std::pow(x, b.re_exponent) * std::pow(pr, b.pr_exponent);
  return {v, b.re_exponent * v / x};
}

bool outside(const RibRange &r, double v) {
  if (!r.bounded()) {
    return r.lo > 0.0 && v < r.lo * (1.0 - RibRange::rel_tol);
  }
  return v < r.lo * (1.0 - RibRange::rel_tol) || v > r.hi * (1.0 + RibRange::rel_tol);
}

// C1 smoothstep weight in ln x over [x_lo, x_hi]: 0 at and below x_lo, 1 at
// and above x_hi, zero slope at both ends.
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

void require_term(const std::string &set_name, const RibTerm &t,
                  const char *field) {
  if (!(t.reference > 0.0) || !std::isfinite(t.reference) ||
      !std::isfinite(t.exponent)) {
    throw std::invalid_argument("rib ratio set '" + set_name + "': " +
                                std::string(field) +
                                " needs a finite positive reference and a "
                                "finite exponent");
  }
}

void require_baseline(const std::string &set_name, const RatioBaseline &b,
                      const char *field) {
  if (!(b.coeff > 0.0) || !std::isfinite(b.coeff) ||
      !std::isfinite(b.re_exponent) || !std::isfinite(b.pr_exponent)) {
    throw std::invalid_argument(
        "rib ratio set '" + set_name + "': " + std::string(field) +
        " must have a finite positive coeff and finite exponents. A ratio "
        "without the baseline it was fitted to is not a number.");
  }
}

void check_accuracy(const std::string &set_name, const char *field,
                    const StatedAccuracy &a) {
  const bool has_value = std::isfinite(a.value);
  const bool claims_value = a.provenance != AccuracyProvenance::Unstated;
  if (has_value != claims_value || (claims_value && !(a.value >= 0.0))) {
    throw std::invalid_argument(
        "rib ratio set '" + set_name + "': " + std::string(field) +
        " must carry a finite, non-negative value if and only if its "
        "provenance is not Unstated");
  }
}

}  // namespace

RibRatioResult evaluate_rib_ratio(const RibRatioSet &set,
                                  const RibGeometry &geom, double Re,
                                  double Pr, const RibRatioOptions &options) {
  RibRatioResult out;

  // Smoothed |Re| and its derivative: even values, odd derivatives.
  const double x = std::sqrt(Re * Re + RATIO_RE_EPS * RATIO_RE_EPS);
  const double dx_dRe = Re / x;
  const double pr = smooth_pr(Pr);

  const double gN = geometry_term(set.Nu_eD, geom.e_D) *
                    geometry_term(set.Nu_pe, geom.p_e) *
                    geometry_term(set.Nu_WH, geom.W_H) *
                    geometry_term(set.Nu_alpha,
                                  geom.alpha_deg / 90.0 * set.Nu_alpha.reference);
  const double gF = geometry_term(set.f_eD, geom.e_D) *
                    geometry_term(set.f_pe, geom.p_e) *
                    geometry_term(set.f_WH, geom.W_H) *
                    geometry_term(set.f_alpha,
                                  geom.alpha_deg / 90.0 * set.f_alpha.reference);

  const double x_hi = set.Re_floor;
  const double x_lo = set.Re_floor / RATIO_BLEND_FACTOR;
  const ValueSlope w = blend_weight(x, x_lo, x_hi);

  // ---- Nu: r * Nu0_source above the floor, k * Nu0_ext below ----
  const ValueSlope rN = ratio_law(set.C_Nu, gN, set.Nu_Re, x);
  const ValueSlope Ns = baseline_law(set.Nu0_source, x, pr);
  const ValueSlope A = {rN.v * Ns.v, rN.d * Ns.v + rN.v * Ns.d};
  out.ratio_Nu = rN.v;

  double nu = A.v;
  double dnu = A.d;
  if (w.v < 1.0) {
    auto ext = [&](double at) -> ValueSlope {
      switch (options.below_floor) {
        case RatioBelowFloor::Gnielinski: {
          const auto g = nusselt_gnielinski_smooth_with_derivative(at, pr);
          return {g.Nu, g.dNu_dRe};
        }
        case RatioBelowFloor::User:
          return baseline_law(options.user_Nu0, at, pr);
        case RatioBelowFloor::SourceBaseline:
        default:
          return baseline_law(set.Nu0_source, at, pr);
      }
    };
    ValueSlope E = ext(x);
    ValueSlope E_floor = ext(x_hi);
    if (!(E_floor.v > 0.0) || !std::isfinite(E_floor.v)) {
      // Unvalidated options (validate_rib_ratio_options rejects them): keep
      // the evaluator total by falling back to the source baseline.
      E = baseline_law(set.Nu0_source, x, pr);
      E_floor = baseline_law(set.Nu0_source, x_hi, pr);
    }
    const ValueSlope rN_floor = ratio_law(set.C_Nu, gN, set.Nu_Re, x_hi);
    const ValueSlope Ns_floor = baseline_law(set.Nu0_source, x_hi, pr);
    // k matches the source form's VALUE at the floor; the blend supplies C1.
    const double k = rN_floor.v * Ns_floor.v / E_floor.v;
    nu = w.v * A.v + (1.0 - w.v) * k * E.v;
    dnu = w.d * (A.v - k * E.v) + w.v * A.d + (1.0 - w.v) * k * E.d;
  }
  out.Nu = nu;
  out.dNu_dRe = dnu * dx_dRe;

  // ---- f: ratio held at its floor value below the range ----
  const double xf_hi = set.Re_floor_f > 0.0 ? set.Re_floor_f : set.Re_floor;
  const ValueSlope wf = blend_weight(x, xf_hi / RATIO_BLEND_FACTOR, xf_hi);
  const ValueSlope rF = ratio_law(set.C_f, gF, set.f_Re, x);
  const ValueSlope rF_floor = ratio_law(set.C_f, gF, set.f_Re, xf_hi);
  const double rf = wf.v * rF.v + (1.0 - wf.v) * rF_floor.v;
  const double drf = wf.d * (rF.v - rF_floor.v) + wf.v * rF.d;
  const ValueSlope F0 = baseline_law(set.f0_source, x, pr);
  out.ratio_f = rf;
  out.f = rf * F0.v;
  out.df_dRe = (drf * F0.v + rf * F0.d) * dx_dRe;

  out.below_floor = x < set.Re_floor;
  out.extrapolated = out.below_floor || outside(set.valid_Re, std::abs(Re)) ||
                     outside(set.valid_eD, geom.e_D) ||
                     outside(set.valid_pe, geom.p_e) ||
                     outside(set.valid_WH, geom.W_H) ||
                     outside(set.valid_alpha, geom.alpha_deg);
  return out;
}

namespace {

// Taslim & Spring (1987): figure 9 markers (C, the Re-normalised ratio) and
// figure 11 means (Fanning f), per configuration, with the tested Re range
// of the Nu data (figures 4/5) and of the friction data (figure 11).
struct TaslimRow {
  double ar, e_D, C_Nu, f_mean, re_nu_lo, re_f_lo, re_hi;
};

// Re bounds are the digitised extremes rounded OUTWARD, so every measured
// point sits inside its own configuration's box.
constexpr TaslimRow kTaslim1987[] = {
    {0.5, 0.125, 3.34316, 0.10146, 28119.0, 10831.0, 108478.0},
    {0.5, 0.250, 3.90148, 0.51645, 21544.0, 20850.0, 61550.0},
    {1.0, 0.083, 3.42519, 0.04680, 26292.0, 21287.0, 102470.0},
    {1.0, 0.167, 3.99269, 0.11610, 27820.0, 21280.0, 103921.0},
    {3.5, 0.053, 3.03740, 0.01498, 36460.0, 52396.0, 209006.0},
    {3.5, 0.107, 3.12189, 0.02090, 33668.0, 51910.0, 192254.0},
    {3.5, 0.161, 3.33653, 0.03162, 37590.0, 52445.0, 217071.0},
};

// Re normaliser for both ratios, the paper's own Re_ref.
constexpr double kTaslimReRef = 1.0e4;
// Declared Fanning friction baseline the constant f is carried on: Han's
// smooth-channel f(FD), the same f_s rib_channel_friction uses.
constexpr double kFanningF0Coeff = HAN_SMOOTH_FANNING_COEF;
constexpr double kFanningF0Exp = HAN_SMOOTH_FANNING_RE_EXP;

}  // namespace

RibRatioSet taslim_spring_1987(double aspect_ratio_taslim, double e_D) {
  const TaslimRow *row = nullptr;
  for (const auto &r : kTaslim1987) {
    if (std::abs(r.ar - aspect_ratio_taslim) < 1e-9 && std::abs(r.e_D - e_D) < 1e-6) {
      row = &r;
      break;
    }
  }
  if (row == nullptr) {
    throw std::invalid_argument(
        "taslim_spring_1987: no two-side configuration at AR " +
        std::to_string(aspect_ratio_taslim) + ", e/D " + std::to_string(e_D) +
        ". Tested: AR 0.5 (e/D 0.125, 0.25), 1.0 (0.083, 0.167), 3.5 (0.053, "
        "0.107, 0.161).");
  }
  char tag[64];
  std::snprintf(tag, sizeof(tag), "taslim_spring_1987_ar%.1f_eD%.3f", row->ar,
                row->e_D);

  RibRatioSet s;
  s.name = tag;
  s.source = "Taslim, M.E. and Spring, S.D. (1987), AIAA-87-2009. Form stated "
             "(Nu_T ~ Re^0.6, D-B normalised); C from Fig. 9 markers, f from "
             "Fig. 11 means (digitised)";
  s.validity_source = s.source;
  s.provenance = RibProvenance::Fitted;
  s.shape = RibShape::Transverse;
  s.symmetric = true;

  s.C_Nu = row->C_Nu;
  s.Nu_Re = {-0.2, kTaslimReRef};
  s.Nu0_source = {0.023, 0.8, 0.4};  // Dittus-Boelter, as Taslim normalises

  // f is Re-independent: f/f0 = C_f (Re/ref)^0.2 cancels f0's Re^-0.2 exactly.
  s.C_f = row->f_mean / (kFanningF0Coeff * std::pow(kTaslimReRef, kFanningF0Exp));
  s.f_Re = {0.2, kTaslimReRef};
  s.f0_source = {kFanningF0Coeff, kFanningF0Exp, 0.0};

  s.Re_floor = row->re_nu_lo;
  s.Re_floor_f = row->re_f_lo;
  s.valid_Re = {std::min(row->re_nu_lo, row->re_f_lo), row->re_hi};
  s.valid_eD = {row->e_D, row->e_D};
  s.valid_pe = {10.0, 10.0};
  s.valid_WH = {1.0 / row->ar, 1.0 / row->ar};
  s.valid_alpha = {90.0, 90.0};
  s.valid_Pr = 0.7;
  // No accuracy is stated; the scorecard reports what the data says.
  s.accuracy_Nu = StatedAccuracy::unstated();
  s.accuracy_f = StatedAccuracy::unstated();
  return s;
}

void validate_rib_ratio_set(const RibRatioSet &set) {
  if (set.name.empty()) {
    throw std::invalid_argument("rib ratio set: name must not be empty");
  }
  if (set.source.empty()) {
    throw std::invalid_argument("rib ratio set '" + set.name +
                                "': source must not be empty");
  }
  if (!(set.C_Nu > 0.0) || !std::isfinite(set.C_Nu) || !(set.C_f > 0.0) ||
      !std::isfinite(set.C_f)) {
    throw std::invalid_argument(
        "rib ratio set '" + set.name +
        "': C_Nu and C_f must be finite and positive -- a non-positive "
        "multiplier makes Nu or f negative");
  }
  if (!(set.Re_floor > 0.0) || !std::isfinite(set.Re_floor)) {
    throw std::invalid_argument("rib ratio set '" + set.name +
                                "': Re_floor must be finite and positive");
  }
  if (!(set.Re_floor_f >= 0.0) || !std::isfinite(set.Re_floor_f)) {
    throw std::invalid_argument(
        "rib ratio set '" + set.name +
        "': Re_floor_f must be finite and non-negative (0 = Re_floor)");
  }
  require_baseline(set.name, set.Nu0_source, "Nu0_source");
  require_baseline(set.name, set.f0_source, "f0_source");
  require_term(set.name, set.Nu_Re, "Nu_Re");
  require_term(set.name, set.Nu_eD, "Nu_eD");
  require_term(set.name, set.Nu_pe, "Nu_pe");
  require_term(set.name, set.Nu_WH, "Nu_WH");
  require_term(set.name, set.Nu_alpha, "Nu_alpha");
  require_term(set.name, set.f_Re, "f_Re");
  require_term(set.name, set.f_eD, "f_eD");
  require_term(set.name, set.f_pe, "f_pe");
  require_term(set.name, set.f_WH, "f_WH");
  require_term(set.name, set.f_alpha, "f_alpha");
  check_accuracy(set.name, "accuracy_Nu", set.accuracy_Nu);
  check_accuracy(set.name, "accuracy_f", set.accuracy_f);
}

void validate_rib_ratio_options(const RibRatioOptions &options) {
  if (options.below_floor == RatioBelowFloor::User) {
    require_baseline("options", options.user_Nu0, "user_Nu0");
  }
}

}  // namespace cooling
}  // namespace combaero
