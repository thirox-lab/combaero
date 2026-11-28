#include "area_change.h"

#include "composition.h"
#include "dual_number.h"
#include "inner_solve.h"
#include "input_validation.h"
#include "math_constants.h"  // MSVC compatibility for M_PI
#include "thermo.h"
#include "transport.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace combaero {

namespace {

using solver::DualN;

constexpr double kEps = 1e-12;

// --- scalar helpers shared by double and DualN -------------------------------

inline double sqrt_(double a) { return std::sqrt(a); }
inline double exp_(double a) { return std::exp(a); }
inline double log_(double a) { return std::log(a); }
inline double log1p_(double a) { return std::log1p(a); }
inline double expm1_(double a) { return std::expm1(a); }

template <int N> DualN<N> sqrt_(const DualN<N>& a) { return solver::dsqrt(a); }
template <int N> DualN<N> exp_(const DualN<N>& a) { return solver::dexp(a); }
template <int N> DualN<N> log_(const DualN<N>& a) {
    DualN<N> r;
    r.v = std::log(a.v);
    const double c = 1.0 / a.v;
    for (int i = 0; i < N; ++i) r.d[i] = c * a.d[i];
    return r;
}
template <int N> DualN<N> log1p_(const DualN<N>& a) {
    DualN<N> r;
    r.v = std::log1p(a.v);
    const double c = 1.0 / (1.0 + a.v);
    for (int i = 0; i < N; ++i) r.d[i] = c * a.d[i];
    return r;
}
template <int N> DualN<N> expm1_(const DualN<N>& a) {
    DualN<N> r;
    r.v = std::expm1(a.v);
    const double c = std::exp(a.v);
    for (int i = 0; i < N; ++i) r.d[i] = c * a.d[i];
    return r;
}

inline double value_of(double a) { return a; }
template <int N> double value_of(const DualN<N>& a) { return a.v; }

// A function of one variable with its slope, lifted onto S: f(a) with
// df = slope da.
inline double lift(double f, double slope, double a) {
    (void)slope;
    (void)a;
    return f;
}
template <int N> DualN<N> lift(double f, double slope, const DualN<N>& a) {
    DualN<N> r;
    r.v = f;
    for (int i = 0; i < N; ++i) r.d[i] = slope * a.d[i];
    return r;
}

// --- Idelchik (1966) diagram 4-1, sudden expansion at low Re -----------------

// Both low-Re tables print 6 area ratios (0.1 ... 0.6) and, with the laminar
// knot at Re 8, 14 Reynolds numbers.
constexpr std::size_t kCols = 6;
constexpr std::size_t kKnots = 14;
constexpr std::array<double, kCols> kTableN = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};

// Diagram 4-1, sudden expansion: as printed, Re 10 ... 3500. The 3500 column
// is Borda-Carnot rounded to two decimals (0.50 for 0.49 at n 0.3); the knot
// takes (1 - n)^2 exactly.
constexpr std::array<double, kKnots> kExpRe = {8.0,    10.0,   15.0,   20.0,  30.0,
                                               40.0,   50.0,   100.0,  200.0, 500.0,
                                               1000.0, 2000.0, 3000.0, 3500.0};
constexpr double kExpZeta[kCols][kKnots - 1] = {
    {3.10, 3.20, 3.00, 2.40, 2.15, 1.95, 1.70, 1.65, 1.70, 2.00, 1.60, 1.00, 0.81},
    {3.10, 3.20, 2.80, 2.20, 1.85, 1.65, 1.40, 1.30, 1.30, 1.60, 1.25, 0.70, 0.64},
    {3.10, 3.10, 2.60, 2.00, 1.60, 1.40, 1.20, 1.10, 1.10, 1.30, 0.95, 0.60, 0.50},
    {3.10, 3.00, 2.40, 1.80, 1.50, 1.30, 1.10, 1.00, 0.85, 1.05, 0.80, 0.40, 0.36},
    {3.10, 2.80, 2.30, 1.65, 1.35, 1.15, 0.90, 0.75, 0.65, 0.90, 0.65, 0.30, 0.25},
    {3.10, 2.70, 2.15, 1.55, 1.25, 1.05, 0.80, 0.60, 0.40, 0.60, 0.50, 0.20, 0.16},
};

// Diagram 3-10 (Karev), sudden contraction: as printed, Re 10 ... 1e4. Of the
// two columns headed 10^4 the last is 3-9's 0.5 (1 - n), the knot here, taken
// exactly; the other (0.50 at n 0.1 against 0.45, equal elsewhere) is not used.
constexpr std::array<double, kKnots> kConRe = {8.0,    10.0,   20.0,   30.0,   40.0,
                                               50.0,   100.0,  200.0,  500.0,  1000.0,
                                               2000.0, 4000.0, 5000.0, 10000.0};
constexpr double kConZeta[kCols][kKnots - 1] = {
    {5.00, 3.20, 2.40, 2.00, 1.80, 1.30, 1.04, 0.82, 0.64, 0.50, 0.80, 0.75, 0.45},
    {5.00, 3.10, 2.30, 1.84, 1.62, 1.20, 0.95, 0.70, 0.50, 0.40, 0.60, 0.60, 0.40},
    {5.00, 2.95, 2.15, 1.70, 1.50, 1.10, 0.85, 0.60, 0.44, 0.30, 0.55, 0.55, 0.35},
    {5.00, 2.80, 2.00, 1.60, 1.40, 1.00, 0.78, 0.50, 0.35, 0.25, 0.45, 0.50, 0.30},
    {5.00, 2.70, 1.80, 1.46, 1.30, 0.90, 0.65, 0.42, 0.30, 0.20, 0.40, 0.42, 0.25},
    {5.00, 2.60, 1.70, 1.35, 1.20, 0.80, 0.56, 0.35, 0.24, 0.15, 0.35, 0.35, 0.20},
};

enum class LowRe { Expansion, Contraction };

const std::array<double, kKnots>& table_re(LowRe t) {
    return t == LowRe::Expansion ? kExpRe : kConRe;
}
double table_laminar(LowRe t) {
    return t == LowRe::Expansion ? EXPANSION_LAMINAR_COEFF : CONTRACTION_LAMINAR_COEFF;
}
// The high-Re zeta the table closes on: Borda-Carnot, or 3-9's 0.5 (1 - n).
double table_high(LowRe t, double n) {
    return t == LowRe::Expansion ? (1.0 - n) * (1.0 - n) : 0.5 * (1.0 - n);
}
double table_printed(LowRe t, std::size_t c, std::size_t i) {
    return t == LowRe::Expansion ? kExpZeta[c][i] : kConZeta[c][i];
}

// One column as a cubic Hermite of y = ln(zeta Re^2) in x = ln Re.
struct HermiteColumn {
    std::array<double, kKnots> x{};
    std::array<double, kKnots> y{};
    std::array<double, kKnots> d{};
};

HermiteColumn build_column(LowRe t, std::size_t c) {
    HermiteColumn h;
    const double n = kTableN[c];
    const auto& re = table_re(t);
    for (std::size_t i = 0; i < kKnots; ++i) {
        double zeta = 0.0;
        if (i == 0) {
            zeta = table_laminar(t) / re[i];
        } else if (i == kKnots - 1) {
            zeta = table_high(t, n);
        } else {
            zeta = table_printed(t, c, i - 1);
        }
        h.x[i] = std::log(re[i]);
        h.y[i] = std::log(zeta * re[i] * re[i]);
    }
    // Every secant exponent at least LOW_RE_MIN_SECANT_EXPONENT, working back
    // from the high-Re anchor: a knot above that line is lowered (4-1 only;
    // every 3-10 secant is above it).
    for (std::size_t i = kKnots - 1; i-- > 0;) {
        const double cap = h.y[i + 1] - LOW_RE_MIN_SECANT_EXPONENT * (h.x[i + 1] - h.x[i]);
        h.y[i] = std::min(h.y[i], cap);
    }
    std::array<double, kKnots - 1> sec{};
    for (std::size_t i = 0; i + 1 < kKnots; ++i) {
        sec[i] = (h.y[i + 1] - h.y[i]) / (h.x[i + 1] - h.x[i]);
    }
    h.d[0] = 1.0;           // the laminar law, zeta Re^2 ~ Re
    h.d[kKnots - 1] = 2.0;  // the high-Re constant zeta, zeta Re^2 ~ Re^2
    for (std::size_t i = 1; i + 1 < kKnots; ++i) {
        const double a = sec[i - 1];
        const double b = sec[i];
        h.d[i] = (a > 0.0 && b > 0.0) ? 2.0 * a * b / (a + b) : 0.0;
    }
    return h;
}

const std::array<HermiteColumn, kCols>& table_columns(LowRe t) {
    auto build = [](LowRe tt) {
        std::array<HermiteColumn, kCols> out;
        for (std::size_t c = 0; c < kCols; ++c) {
            out[c] = build_column(tt, c);
        }
        return out;
    };
    static const std::array<HermiteColumn, kCols> exp_cols = build(LowRe::Expansion);
    static const std::array<HermiteColumn, kCols> con_cols = build(LowRe::Contraction);
    return t == LowRe::Expansion ? exp_cols : con_cols;
}

// A column's ratio to its high-Re zeta, R = zeta / zeta_high(n_c), with dR/dRe
// and d2R/dRe2 (the second for the drive floor's exponent, area_change_flow).
struct Ratio {
    double R = 1.0;
    double dR = 0.0;
    double d2R = 0.0;
};

Ratio column_ratio(LowRe t, double Re, std::size_t c) {
    const double b = table_high(t, kTableN[c]);
    const auto& re = table_re(t);
    if (Re >= re.back()) {
        return {1.0, 0.0, 0.0};
    }
    if (Re <= re.front()) {
        const double R = table_laminar(t) / (Re * b);
        return {R, -R / Re, 2.0 * R / (Re * Re)};
    }
    const HermiteColumn& h = table_columns(t)[c];
    const double x = std::log(Re);
    std::size_t i = 0;
    while (i + 2 < kKnots && x > h.x[i + 1]) {
        ++i;
    }
    const double w = h.x[i + 1] - h.x[i];
    const double u = (x - h.x[i]) / w;
    const double u2 = u * u;
    const double u3 = u2 * u;
    const double y = (2.0 * u3 - 3.0 * u2 + 1.0) * h.y[i] + (u3 - 2.0 * u2 + u) * w * h.d[i] +
                     (-2.0 * u3 + 3.0 * u2) * h.y[i + 1] + (u3 - u2) * w * h.d[i + 1];
    const double dy = ((6.0 * u2 - 6.0 * u) * h.y[i] + (3.0 * u2 - 4.0 * u + 1.0) * w * h.d[i] +
                       (-6.0 * u2 + 6.0 * u) * h.y[i + 1] + (3.0 * u2 - 2.0 * u) * w * h.d[i + 1]) /
                      w;
    const double d2y = ((12.0 * u - 6.0) * h.y[i] + (6.0 * u - 4.0) * w * h.d[i] +
                        (-12.0 * u + 6.0) * h.y[i + 1] + (6.0 * u - 2.0) * w * h.d[i + 1]) /
                       (w * w);
    // R(x) = exp(y - 2x) / b, x = ln Re.
    const double R = std::exp(y - 2.0 * x) / b;
    const double R_x = R * (dy - 2.0);
    const double R_xx = R * ((dy - 2.0) * (dy - 2.0) + d2y);
    return {R, R_x / Re, (R_xx - R_x) / (Re * Re)};
}

// The ratio at area ratio n: linear in n between printed columns, the edge
// column's outside them.
Ratio table_ratio3(LowRe t, double Re, double n) {
    if (n <= kTableN.front()) {
        return column_ratio(t, Re, 0);
    }
    if (n >= kTableN.back()) {
        return column_ratio(t, Re, kCols - 1);
    }
    std::size_t j = 0;
    while (j + 2 < kCols && n > kTableN[j + 1]) {
        ++j;
    }
    const double w = (n - kTableN[j]) / (kTableN[j + 1] - kTableN[j]);
    const Ratio a = column_ratio(t, Re, j);
    const Ratio b = column_ratio(t, Re, j + 1);
    return {(1.0 - w) * a.R + w * b.R, (1.0 - w) * a.dR + w * b.dR,
            (1.0 - w) * a.d2R + w * b.d2R};
}

ZetaAndSlope table_ratio(LowRe t, double Re, double n) {
    const Ratio r = table_ratio3(t, Re, n);
    return {r.R, r.dR};
}

// The incompressible drive's exponent in Re, e = d ln((R + offset) Re^2) / d ln Re,
// and de/dRe: the flow then goes as drive^(1/e). offset is 0 for total-to-total
// and r_exit^2 / zeta_inc when exit_head_lost adds the exit dynamic head (#1128).
std::pair<double, double> table_exponent(LowRe t, double Re, double n,
                                         double offset_over_zeta_inc = 0.0) {
    const Ratio r = table_ratio3(t, Re, n);
    const double c = r.R + offset_over_zeta_inc;
    const double e = 2.0 + Re * r.dR / c;
    const double de = r.dR / c + Re * r.d2R / c - Re * r.dR * r.dR / (c * c);
    return {e, de};
}

void check_area_ratio(double n, const char* who) {
    if (!(n >= 0.0) || n > 1.0) {
        throw std::invalid_argument(std::string(who) + ": area ratio must be in [0, 1]");
    }
}

// --- Idelchik (1966) diagram 5-2, conical diffuser ---------------------------

constexpr std::size_t kConeAngles = 14;
constexpr std::size_t kConeRows = 11;
constexpr std::array<double, kConeAngles> kConeAlpha = {3.0,  6.0,  8.0,  10.0, 12.0,
                                                        14.0, 16.0, 20.0, 24.0, 30.0,
                                                        40.0, 60.0, 90.0, 180.0};
constexpr std::array<double, kConeRows> kConeN = {0.0,  0.05, 0.075, 0.1, 0.15, 0.2,
                                                  0.25, 0.3,  0.4,   0.5, 0.6};
constexpr double kConeZeta[kConeRows][kConeAngles] = {
    {0.03, 0.08, 0.11, 0.15, 0.19, 0.23, 0.27, 0.36, 0.47, 0.65, 0.92, 1.15, 1.10, 1.02},
    {0.03, 0.07, 0.10, 0.14, 0.16, 0.20, 0.24, 0.32, 0.42, 0.58, 0.83, 1.04, 0.99, 0.92},
    {0.03, 0.07, 0.09, 0.13, 0.16, 0.19, 0.23, 0.30, 0.40, 0.55, 0.79, 0.99, 0.95, 0.88},
    {0.03, 0.07, 0.09, 0.12, 0.15, 0.18, 0.22, 0.29, 0.38, 0.52, 0.75, 0.93, 0.89, 0.83},
    {0.02, 0.06, 0.08, 0.11, 0.14, 0.17, 0.20, 0.26, 0.34, 0.46, 0.67, 0.84, 0.79, 0.74},
    {0.02, 0.05, 0.07, 0.10, 0.12, 0.15, 0.17, 0.23, 0.30, 0.41, 0.59, 0.74, 0.70, 0.65},
    {0.02, 0.05, 0.06, 0.08, 0.10, 0.13, 0.15, 0.20, 0.26, 0.35, 0.47, 0.65, 0.62, 0.58},
    {0.02, 0.04, 0.05, 0.07, 0.09, 0.11, 0.13, 0.18, 0.23, 0.31, 0.40, 0.57, 0.54, 0.50},
    {0.01, 0.03, 0.04, 0.06, 0.07, 0.08, 0.10, 0.13, 0.17, 0.23, 0.33, 0.41, 0.39, 0.37},
    {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.09, 0.12, 0.16, 0.23, 0.29, 0.28, 0.26},
    {0.01, 0.01, 0.02, 0.03, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.15, 0.18, 0.17, 0.16},
};

// phi = zeta_exp / (1 - n)^2 at a printed row, linear in alpha.
double cone_phi_row(std::size_t r, double alpha) {
    const double b = (1.0 - kConeN[r]) * (1.0 - kConeN[r]);
    if (alpha <= kConeAlpha.front()) {
        return kConeZeta[r][0] / b * std::max(alpha, 0.0) / kConeAlpha.front();
    }
    if (alpha >= kConeAlpha.back()) {
        return kConeZeta[r][kConeAngles - 1] / b;
    }
    std::size_t i = 0;
    while (i + 2 < kConeAngles && alpha > kConeAlpha[i + 1]) {
        ++i;
    }
    const double w = (alpha - kConeAlpha[i]) / (kConeAlpha[i + 1] - kConeAlpha[i]);
    return ((1.0 - w) * kConeZeta[r][i] + w * kConeZeta[r][i + 1]) / b;
}

double cone_phi(double alpha, double n) {
    if (n >= kConeN.back()) {
        return cone_phi_row(kConeRows - 1, alpha);
    }
    std::size_t j = 0;
    while (j + 2 < kConeRows && n > kConeN[j + 1]) {
        ++j;
    }
    const double w = (n - kConeN[j]) / (kConeN[j + 1] - kConeN[j]);
    return (1.0 - w) * cone_phi_row(j, alpha) + w * cone_phi_row(j + 1, alpha);
}

// --- Idelchik (1966) diagram 3-6, conical converging bellmouth with end wall ---

constexpr std::size_t kNozAngles = 9;
constexpr std::size_t kNozLengths = 6;
constexpr std::array<double, kNozAngles> kNozAlpha = {0.0,  10.0,  20.0,  30.0, 40.0,
                                                      60.0, 100.0, 140.0, 180.0};
constexpr std::array<double, kNozLengths> kNozLD = {0.025, 0.05, 0.075, 0.10, 0.15, 0.60};
constexpr double kNozZeta[kNozLengths][kNozAngles] = {
    {0.50, 0.47, 0.45, 0.43, 0.41, 0.40, 0.42, 0.45, 0.50},
    {0.50, 0.45, 0.41, 0.36, 0.33, 0.30, 0.35, 0.42, 0.50},
    {0.50, 0.42, 0.35, 0.30, 0.26, 0.23, 0.30, 0.40, 0.50},
    {0.50, 0.39, 0.32, 0.25, 0.22, 0.18, 0.27, 0.38, 0.50},
    {0.50, 0.37, 0.27, 0.20, 0.16, 0.15, 0.25, 0.37, 0.50},
    {0.50, 0.27, 0.18, 0.13, 0.11, 0.12, 0.23, 0.36, 0.50},
};

// zeta'(alpha, l/D): bilinear, l/D held at the table's ends.
double nozzle_prime(double alpha, double lD) {
    alpha = std::clamp(alpha, 0.0, 180.0);
    lD = std::clamp(lD, kNozLD.front(), kNozLD.back());
    std::size_t i = 0;
    while (i + 2 < kNozAngles && alpha > kNozAlpha[i + 1]) {
        ++i;
    }
    std::size_t j = 0;
    while (j + 2 < kNozLengths && lD > kNozLD[j + 1]) {
        ++j;
    }
    const double wa = (alpha - kNozAlpha[i]) / (kNozAlpha[i + 1] - kNozAlpha[i]);
    const double wl = (lD - kNozLD[j]) / (kNozLD[j + 1] - kNozLD[j]);
    auto row = [&](std::size_t r) {
        return (1.0 - wa) * kNozZeta[r][i] + wa * kNozZeta[r][i + 1];
    };
    return (1.0 - wl) * row(j) + wl * row(j + 1);
}

// --- compressible Borda-Carnot, on any scalar ---------------------------------

template <class S> struct BordaCarnotT {
    S y;       // M2^2
    S ln_p21;  // ln(p2 / p1)
    S ln_pi;   // ln(Pt2 / Pt1)
};

// The subsonic root as 2 c^2 / (b + sqrt(disc)): no cancellation at low Mach,
// and b = 1 - 2 gamma c^2 > 0 for every M1 <= 1. ln(Pt2/Pt1) from log1p
// terms, each relative-accurate, so 1 - Pt2/Pt1 = -expm1(ln_pi) keeps its
// relative precision as M1 -> 0 (it is O(M1^2)).
template <class S> BordaCarnotT<S> borda_carnot_t(const S& M1, double n, const S& g) {
    const S M1sq = M1 * M1;
    const S a1 = 0.5 * (g - 1.0) * M1sq;
    const S c = n * M1 * sqrt_(1.0 + a1) / (1.0 + n * g * M1sq);
    const S c2 = c * c;
    const S b = 1.0 - 2.0 * g * c2;
    const S disc = 1.0 - 2.0 * (g + 1.0) * c2;
    const S y = 2.0 * c2 / (b + sqrt_(disc));
    const S ln_p21 = log1p_(n * g * M1sq) - log1p_(g * y);
    const S ln_pi = ln_p21 + g / (g - 1.0) * (log1p_(0.5 * (g - 1.0) * y) - log1p_(a1));
    return {y, ln_p21, ln_pi};
}

// --- one branch of the element ------------------------------------------------

struct Branch {
    bool expansion = true;
    bool sharp = true;
    bool low_re = false; // sharp: the low-Re table's ratio (4-1 or 3-10)
    bool exit_head_lost = false; // static-pressure exit coupling (#1128)
    LowRe table = LowRe::Expansion;
    double n = 1.0;      // small / large
    double step = 1.0;   // area ratio of the Borda-Carnot step: n, or Cc
    double Cc = 1.0;     // contraction: vena contracta over the small section
    double phi = 1.0;    // conical expansion: 5-2's ratio to Borda-Carnot
    double zeta_inc = 0.0;  // incompressible high-Re zeta (the low-Re table closes on it)
    double r_exit = 0.0;    // F_small / F_out (n for expansion, 1 for contraction)
    double F_small = 0.0;
    double D = 0.0;      // Re length
    double floor = 0.0;  // drive floor [Pa]: (zeta_inc + exit_head) * AREA_CHANGE_HEAD_FLOOR
};

Branch make_branch(double F_in, double F_out, const std::string& model, double length,
                   double D_h, bool exit_head_lost = false) {
    if (!validation::is_finite_positive(F_in) || !validation::is_finite_positive(F_out)) {
        throw std::invalid_argument("area_change: areas must be positive and finite");
    }
    // D_h <= 0 selects the equivalent-circle diameter below; NaN used to
    // select it too (#769).
    if (!validation::is_finite(D_h)) {
        throw std::invalid_argument("area_change: D_h must be finite");
    }
    if (model != "sharp" && model != "conical") {
        throw std::invalid_argument("area_change: model must be 'sharp' or 'conical'");
    }
    Branch b;
    b.F_small = std::min(F_in, F_out);
    const double F_large = std::max(F_in, F_out);
    b.n = b.F_small / F_large;
    if (!(b.n < 1.0)) {
        throw std::invalid_argument(
            "area_change: equal areas are no area change (use a lossless connection)");
    }
    b.expansion = F_in < F_out;
    b.sharp = model == "sharp";
    b.low_re = b.sharp;
    b.exit_head_lost = exit_head_lost;
    b.r_exit = b.expansion ? b.n : 1.0;
    b.table = b.expansion ? LowRe::Expansion : LowRe::Contraction;
    b.D = (D_h > 0.0) ? D_h : std::sqrt(4.0 * b.F_small / M_PI);
    if (b.expansion) {
        b.step = b.n;
        if (!b.sharp) {
            if (!validation::is_finite_positive(length)) {
                throw std::invalid_argument("area_change: a conical model needs length > 0");
            }
            const double dr = std::sqrt(F_large / M_PI) - std::sqrt(b.F_small / M_PI);
            const double alpha = 2.0 * std::atan2(dr, length) * 180.0 / M_PI;
            b.phi = cone_phi(alpha, b.n);
        }
        b.zeta_inc = b.phi * (1.0 - b.n) * (1.0 - b.n);
    } else {
        if (b.sharp) {
            b.zeta_inc = 0.5 * (1.0 - b.n);
        } else {
            if (!validation::is_finite_positive(length)) {
                throw std::invalid_argument("area_change: a conical model needs length > 0");
            }
            const double D_small = std::sqrt(4.0 * b.F_small / M_PI);
            const double dr = std::sqrt(F_large / M_PI) - std::sqrt(b.F_small / M_PI);
            const double alpha = 2.0 * std::atan2(dr, length) * 180.0 / M_PI;
            b.zeta_inc = nozzle_prime(alpha, length / D_small) * (1.0 - b.n);
        }
        b.Cc = 1.0 / (1.0 + std::sqrt(b.zeta_inc));
        b.step = b.Cc;
    }
    const double zeta_drive = b.zeta_inc + (b.exit_head_lost ? b.r_exit * b.r_exit : 0.0);
    b.floor = zeta_drive * AREA_CHANGE_HEAD_FLOOR;
    return b;
}

// Thermodynamic inputs at the upstream stagnation state, frozen across the
// step, with the partials the derivatives need.
struct Gas {
    std::vector<double> X;
    double g = 1.4;
    double Rg = 287.0;
    double g_T = 0.0;
    std::vector<double> g_X;
    std::vector<double> Rg_X;
    double P = 1e5;  // for the transport calls (unused by them)
};

Gas make_gas(double Tt0, double Pt0, const std::vector<double>& Y, bool with_derivatives) {
    Gas gas;
    gas.X = mass_to_mole(normalize_fractions(Y));
    gas.P = Pt0;
    const double cp_mol = cp(Tt0, gas.X);
    const double cv_mol = cp_mol - thermo::R_GAS;
    gas.g = cp_mol / cv_mol;
    gas.Rg = specific_gas_constant(gas.X);
    if (with_derivatives) {
        const double c = -thermo::R_GAS / (cv_mol * cv_mol);
        gas.g_T = c * dcp_dT(Tt0, gas.X);
        const double MW = mwmix(gas.X);
        gas.g_X.resize(gas.X.size());
        gas.Rg_X.resize(gas.X.size());
        for (std::size_t k = 0; k < gas.X.size(); ++k) {
            gas.g_X[k] = c * cp_species(k, Tt0);
            gas.Rg_X[k] = -gas.Rg * species_molar_mass(k) / MW;
        }
    }
    return gas;
}

// Viscosity lifted onto S: value and dT at T, plus the seed of an additive
// offset (index kMuSeed) that carries d/dX at fixed T.
constexpr int kSeeds = 6;  // the parameter (M or G), Pt0, Tt0, gamma, R, mu offset
using D6 = DualN<kSeeds>;
constexpr int kP = 0, kPt = 1, kTt = 2, kG = 3, kR = 4, kMu = 5;

double mu_at(double T, const Gas& gas) { return viscosity(T, gas.P, gas.X); }
template <int N> DualN<N> mu_at(const DualN<N>& T, const Gas& gas) {
    const TransportDerivatives t = transport_and_dT(T.v, gas.P, gas.X);
    DualN<N> r = lift(t.mu, t.dmu_dT, T);
    if constexpr (N == kSeeds) {
        r.d[kMu] += 1.0;
    }
    return r;
}

template <class S> struct EvalT {
    S loss;         // matched drive: Pt0 - Pt_out (or Pt0 - P_exit if exit_head_lost) [Pa]
    S exit_margin;  // Phi_max(g) * pi_exit - b.step * Phi_1 (exit choke margin, >= 0)
    S G_small;      // mass flux at the small section
    S Re;           // small section (sharp expansion/contraction; else 0)
    double dPt_loss = 0.0; // step's own stagnation loss Pt0 - Pt_exit [Pa]
    double q_small = 0.0;  // static dynamic head at the small section
    double mach = 0.0;     // at the closure's station
};

// Solve subsonic exit Mach Me in [0, 1) from dimensionless isentropic mass
// flux Phi(Me, g) = Me * (1 + (g - 1)/2 Me^2)^(-(g + 1)/(2 (g - 1))) = phi_target.
double solve_exit_mach(double phi_target, double g, double M_guess) {
    if (!(phi_target > 0.0)) {
        return 0.0;
    }
    double lo = 0.0;
    double hi = 1.0;
    double M = std::clamp(M_guess, 1e-15, 1.0 - 1e-12);
    const double exp_pow = -(g + 1.0) / (2.0 * (g - 1.0));
    for (int it = 0; it < inner_solve::kAreaChangeMaxIter; ++it) {
        const double x = 0.5 * (g - 1.0) * M * M;
        const double base = std::exp(exp_pow * std::log1p(x));
        const double f = M * base - phi_target;
        if (f > 0.0) {
            hi = M;
        } else {
            lo = M;
        }
        const double fp = base * (1.0 - M * M) / (1.0 + x);
        double M_new = (fp > 0.0) ? M - f / fp : 0.5 * (lo + hi);
        if (!(M_new > lo && M_new < hi)) {
            M_new = 0.5 * (lo + hi);
        }
        if (std::abs(M_new - M) <= 1e-15 * std::max(1.0, M) || hi - lo <= 1e-15) {
            return M_new;
        }
        M = M_new;
    }
    return M;
}

inline double lift_exit_mach(double Me_v, double phi_e, double g) {
    (void)phi_e;
    (void)g;
    return Me_v;
}

template <int N>
DualN<N> lift_exit_mach(double Me_v, const DualN<N>& phi_e, const DualN<N>& g) {
    const DualN<N> x = 0.5 * (g - 1.0) * (Me_v * Me_v);
    const DualN<N> base = exp_(-(g + 1.0) / (2.0 * (g - 1.0)) * log1p_(x));
    const DualN<N> phi_fixed = Me_v * base;
    const double dPhi_dMe = base.v * (1.0 - Me_v * Me_v) / (1.0 + x.v);
    DualN<N> Me;
    Me.v = Me_v;
    for (int i = 0; i < N; ++i) {
        Me.d[i] = (phi_e.d[i] - phi_fixed.d[i]) / dPhi_dMe;
    }
    return Me;
}

// Compressible: parametrised by the Mach number M at the closure's station
// (the small section for an expansion, the vena contracta for a contraction).
template <class S>
EvalT<S> eval_compressible(const S& M, const S& Pt0, const S& Tt0, const S& g, const S& Rg,
                           const Branch& b, const Gas& gas, bool exit_choked_at_face = false) {
    const S Msq = M * M;
    const S lnX = log1p_(0.5 * (g - 1.0) * Msq);
    const S phi_1 = M * exp_(-(g + 1.0) / (2.0 * (g - 1.0)) * lnX);
    const S Ga = Pt0 * sqrt_(g / (Rg * Tt0)) * phi_1;
    const BordaCarnotT<S> bc = borda_carnot_t(M, b.step, g);
    S dPt_loss = -Pt0 * expm1_(bc.ln_pi);
    S mult = 1.0 + 0.0 * M;
    const S G_small = b.expansion ? Ga : b.Cc * Ga;
    S Re = 0.0 * M;
    const double Pa = value_of(Pt0) * std::exp(-value_of(g) / (value_of(g) - 1.0) * value_of(lnX));
    double q_small = 0.0;
    if (b.low_re) {
        // The low-Re table's ratio at the small section's Re (mu at the
        // closure's station temperature).
        const S Ta = Tt0 / exp_(lnX);
        Re = G_small * b.D / mu_at(Ta, gas);
        const ZetaAndSlope r = table_ratio(b.table, value_of(Re), b.n);
        mult = mult * lift(r.zeta, r.dzeta_dRe, Re);
    }
    if (b.expansion) {
        q_small = 0.5 * value_of(g) * Pa * value_of(Msq);
        if (!b.sharp) {
            mult = mult * b.phi;
        }
    } else {
        // The small section is the step's exit, station 2.
        const double p2 = Pa * std::exp(value_of(bc.ln_p21));
        q_small = 0.5 * value_of(g) * p2 * value_of(bc.y);
    }
    dPt_loss = dPt_loss * mult;
    S loss = dPt_loss;
    S exit_margin = 0.0 * M;
    if (b.exit_head_lost) {
        // Exit stagnation pressure ratio pi_exit = Pt_exit / Pt0 = 1 + mult * expm1(ln_pi).
        // Exit face flux parameter Phi(Me) = b.step * Phi(M) / pi_exit.
        const S pi_exit = 1.0 + expm1_(bc.ln_pi) * mult;
        const S phi_max = exp_(-(g + 1.0) / (2.0 * (g - 1.0)) * log_(0.5 * (g + 1.0)));
        exit_margin = phi_max * pi_exit - b.step * phi_1;
        S lnXe = 0.0 * M;
        if (exit_choked_at_face) {
            lnXe = log1p_(0.5 * (g - 1.0));
        } else {
            const S phi_e = b.step * phi_1 / pi_exit;
            const double Me_v = solve_exit_mach(value_of(phi_e), value_of(g),
                                                std::sqrt(std::max(value_of(bc.y), 0.0)));
            const S Me = lift_exit_mach(Me_v, phi_e, g);
            lnXe = log1p_(0.5 * (g - 1.0) * Me * Me);
        }
        const S q_exit_comp = -Pt0 * pi_exit * expm1_(-g / (g - 1.0) * lnXe);
        loss = dPt_loss + q_exit_comp;
    }
    return {loss, exit_margin, G_small, Re, value_of(dPt_loss), q_small, value_of(M)};
}

// Incompressible (M -> 0): parametrised by the small section's mass flux.
template <class S>
EvalT<S> eval_incompressible(const S& G, const S& Pt0, const S& Tt0, const S& g, const S& Rg,
                             const Branch& b, const Gas& gas) {
    const S rho = Pt0 / (Rg * Tt0);
    const S q = G * G / (2.0 * rho);
    S zeta = b.zeta_inc + 0.0 * G;
    S Re = 0.0 * G;
    if (b.low_re) {
        Re = G * b.D / mu_at(Tt0, gas);
        const ZetaAndSlope r = table_ratio(b.table, value_of(Re), b.n);
        zeta = lift(r.zeta, r.dzeta_dRe, Re) * b.zeta_inc;
    }
    const S dPt_loss = zeta * q;
    const S loss = b.exit_head_lost ? (zeta + b.r_exit * b.r_exit) * q : dPt_loss;
    const double a0 = std::sqrt(value_of(g) * value_of(Rg) * value_of(Tt0));
    const double Ga = value_of(G) / (b.expansion ? 1.0 : b.Cc);
    return {loss, 0.0 * G, G, Re, value_of(dPt_loss), value_of(q), Ga / (value_of(rho) * a0)};
}

template <class S>
EvalT<S> eval_branch(bool incompressible, const S& p, const S& Pt0, const S& Tt0, const S& g,
                     const S& Rg, const Branch& b, const Gas& gas,
                     bool exit_choked_at_face = false) {
    return incompressible ? eval_incompressible(p, Pt0, Tt0, g, Rg, b, gas)
                          : eval_compressible(p, Pt0, Tt0, g, Rg, b, gas, exit_choked_at_face);
}

double exit_margin_at(double M, double Pt0, double Tt0, const Branch& b, const Gas& gas) {
    return eval_compressible<double>(M, Pt0, Tt0, gas.g, gas.Rg, b, gas, true).exit_margin;
}

double solve_exit_choke_mach(double Pt0, double Tt0, const Branch& b, const Gas& gas) {
    double hi = 1.0;
    double lo = 0.5;
    for (int i = 0; i < 200 && exit_margin_at(lo, Pt0, Tt0, b, gas) <= 0.0; ++i) {
        hi = lo;
        lo *= 0.5;
    }
    double u_lo = std::log(lo);
    double u_hi = std::log(hi);
    double u = 0.5 * (u_lo + u_hi);
    for (int it = 0; it < inner_solve::kAreaChangeMaxIter; ++it) {
        const double M = std::exp(u);
        using D1 = DualN<1>;
        const auto e = eval_compressible<D1>(D1::seed(M, 0), D1::constant(Pt0), D1::constant(Tt0),
                                             D1::constant(gas.g), D1::constant(gas.Rg), b, gas,
                                             true);
        const double f = e.exit_margin.v;
        if (f > 0.0) {
            u_lo = u;
        } else {
            u_hi = u;
        }
        const double fp = e.exit_margin.d[0] * M;
        double u_new = (fp < 0.0) ? u - f / fp : 0.5 * (u_lo + u_hi);
        if (!(u_new > u_lo && u_new < u_hi)) {
            u_new = 0.5 * (u_lo + u_hi);
        }
        if (std::abs(u_new - u) <= 1e-15 * std::max(1.0, std::abs(u)) || u_hi - u_lo < 1e-15) {
            return std::exp(u_new);
        }
        u = u_new;
    }
    return std::exp(u);
}

double loss_at(bool inc, double p, double Pt0, double Tt0, const Branch& b, const Gas& gas) {
    return eval_branch<double>(inc, p, Pt0, Tt0, gas.g, gas.Rg, b, gas).loss;
}

// d ln(loss) / d ln(p), for the Newton step on ln p.
double log_slope(bool inc, double p, double Pt0, double Tt0, const Branch& b, const Gas& gas) {
    using D1 = DualN<1>;
    const auto e = eval_branch<D1>(inc, D1::seed(p, 0), D1::constant(Pt0), D1::constant(Tt0),
                                   D1::constant(gas.g), D1::constant(gas.Rg), b, gas);
    return e.loss.d[0] * p / e.loss.v;
}

// The parameter p (Mach, or flux if incompressible) whose loss is the drive,
// on (lo, hi): safeguarded Newton in ln p. The loss rises monotonically in p
// (sudden_expansion_zeta keeps the loss exponent positive), so the root is
// unique.
double solve_param(bool inc, double drive, double hi, double guess, double Pt0, double Tt0,
                   const Branch& b, const Gas& gas) {
    double p_hi = hi;
    double p_lo = std::min(guess, 0.5 * hi);
    for (int i = 0; i < 2000 && loss_at(inc, p_lo, Pt0, Tt0, b, gas) >= drive; ++i) {
        p_hi = p_lo;
        p_lo *= 0.1;
    }
    double u_lo = std::log(p_lo);
    double u_hi = std::log(p_hi);
    double u = std::log(std::clamp(guess, p_lo, p_hi));
    const double ln_drive = std::log(drive);
    for (int it = 0; it < inner_solve::kAreaChangeMaxIter; ++it) {
        const double p = std::exp(u);
        const double f = std::log(loss_at(inc, p, Pt0, Tt0, b, gas)) - ln_drive;
        if (f > 0.0) {
            u_hi = u;
        } else {
            u_lo = u;
        }
        const double fp = log_slope(inc, p, Pt0, Tt0, b, gas);
        double u_new = (fp > 0.0) ? u - f / fp : 0.5 * (u_lo + u_hi);
        if (!(u_new > u_lo && u_new < u_hi)) {
            u_new = 0.5 * (u_lo + u_hi);
        }
        if (std::abs(u_new - u) <= 1e-15 * std::max(1.0, std::abs(u)) || u_hi - u_lo < 1e-15) {
            return std::exp(u_new);
        }
        u = u_new;
    }
    // A cap exit is a root only if the bracket (in ln p) has closed (#868).
    if (inner_solve::cap_exit_is_root(u_hi - u_lo, std::max(1.0, std::abs(u)))) {
        return std::exp(u);
    }
    inner_solve::throw_no_convergence("area_change (solve_param)", inner_solve::kAreaChangeMaxIter,
                                      "bracket width in ln p", u_hi - u_lo);
}

// The flow at a drive >= 0, without a floor. dm/d(drive) and dm/dPt0 at
// fixed drive (the caller converts to Pt_target).
struct RawFlow {
    double m = 0.0;
    bool choked = false;
    double mach = 0.0;
    double Re = 0.0;
    double zeta = 0.0;
    double dm_ddrive = 0.0;
    double dm_dPt0 = 0.0;   // at fixed drive
    double dm_dTt0 = 0.0;   // explicit + through gamma
    double dm_dg = 0.0;
    double dm_dR = 0.0;
    double dm_dmu = 0.0;    // through the viscosity offset
    double T_mu = 0.0;      // where mu was evaluated
    // Re's partials at fixed drive, as the flow's (sharp expansion; else 0).
    double dRe_dPt0 = 0.0;
    double dRe_dTt0 = 0.0;
    double dRe_dg = 0.0;
    double dRe_dR = 0.0;
    double dRe_dmu = 0.0;
};

RawFlow raw_flow(double drive, double Pt0, double Tt0, const Branch& b, const Gas& gas, bool inc,
                 bool with_derivatives) {
    RawFlow out;
    double p = 0.0;
    bool exit_choked_at_face = false;
    const double zeta_est =
        std::max(b.zeta_inc + (b.exit_head_lost ? b.r_exit * b.r_exit : 0.0), 1e-3);
    if (inc) {
        const double rho = Pt0 / (gas.Rg * Tt0);
        const double guess = std::sqrt(2.0 * rho * drive / zeta_est);
        double hi = 2.0 * guess + kEps;
        for (int i = 0; i < 200 && loss_at(true, hi, Pt0, Tt0, b, gas) < drive; ++i) {
            hi *= 4.0;
        }
        p = solve_param(true, drive, hi, guess, Pt0, Tt0, b, gas);
    } else {
        // The sonic limit: closure choke at M = 1, or (at low Re with
        // exit_head_lost) exit-face choke at M_max < 1 where exit_margin = 0.
        double M_max = 1.0;
        if (b.exit_head_lost && exit_margin_at(1.0, Pt0, Tt0, b, gas) < 0.0) {
            M_max = solve_exit_choke_mach(Pt0, Tt0, b, gas);
            exit_choked_at_face = true;
        }
        const double drive_max =
            eval_compressible<double>(M_max, Pt0, Tt0, gas.g, gas.Rg, b, gas, exit_choked_at_face)
                .loss;
        if (drive_max <= drive) {
            p = M_max;
            out.choked = true;
        } else {
            exit_choked_at_face = false;
            const double guess = std::sqrt(2.0 * drive / (gas.g * Pt0 * zeta_est));
            p = solve_param(false, drive, M_max, std::min(guess, 0.5 * M_max), Pt0, Tt0, b, gas);
        }
    }
    const auto e = eval_branch<D6>(inc, D6::seed(p, kP), D6::seed(Pt0, kPt), D6::seed(Tt0, kTt),
                                   D6::seed(gas.g, kG), D6::seed(gas.Rg, kR), b, gas,
                                   out.choked && exit_choked_at_face);
    out.m = b.F_small * e.G_small.v;
    out.mach = e.mach;
    out.Re = e.Re.v;
    out.zeta = (e.q_small > 0.0) ? e.dPt_loss / e.q_small : 0.0;
    out.T_mu = inc ? Tt0 : Tt0 / (1.0 + 0.5 * (gas.g - 1.0) * e.mach * e.mach);
    if (!with_derivatives) {
        return out;
    }
    // loss(p; q) = drive: dp/dq = -L_q / L_p, dp/d(drive) = 1 / L_p.
    // Choked at closure (M = 1): p is pinned (dp = 0).
    // Choked at exit face (exit_margin(M; q) = 0): dp/dq = -margin_q / margin_p.
    const double L_p = e.loss.d[kP];
    auto dp = [&](int q) {
        if (!out.choked) {
            return -e.loss.d[q] / L_p;
        }
        if (exit_choked_at_face) {
            return -e.exit_margin.d[q] / e.exit_margin.d[kP];
        }
        return 0.0;
    };
    auto dm = [&](int q) { return b.F_small * (e.G_small.d[q] + e.G_small.d[kP] * dp(q)); };
    out.dm_ddrive = out.choked ? 0.0 : b.F_small * e.G_small.d[kP] / L_p;
    out.dm_dPt0 = dm(kPt);
    out.dm_dg = dm(kG);
    out.dm_dR = dm(kR);
    out.dm_dmu = dm(kMu);
    out.dm_dTt0 = dm(kTt) + out.dm_dg * gas.g_T;
    auto dRe = [&](int q) { return e.Re.d[q] + e.Re.d[kP] * dp(q); };
    out.dRe_dPt0 = dRe(kPt);
    out.dRe_dg = dRe(kG);
    out.dRe_dR = dRe(kR);
    out.dRe_dmu = dRe(kMu);
    out.dRe_dTt0 = dRe(kTt) + out.dRe_dg * gas.g_T;
    return out;
}

}  // namespace

// --- public helpers ------------------------------------------------------------

ZetaAndSlope sudden_expansion_zeta(double Re, double area_ratio) {
    check_area_ratio(area_ratio, "sudden_expansion_zeta");
    if (!(Re > 0.0)) {
        throw std::invalid_argument("sudden_expansion_zeta: Re must be positive");
    }
    const double b = (1.0 - area_ratio) * (1.0 - area_ratio);
    const ZetaAndSlope r = table_ratio(LowRe::Expansion, Re, area_ratio);
    return {r.zeta * b, r.dzeta_dRe * b};
}

ZetaAndSlope sudden_contraction_zeta(double Re, double area_ratio) {
    check_area_ratio(area_ratio, "sudden_contraction_zeta");
    if (!(Re > 0.0)) {
        throw std::invalid_argument("sudden_contraction_zeta: Re must be positive");
    }
    const double b = 0.5 * (1.0 - area_ratio);
    const ZetaAndSlope r = table_ratio(LowRe::Contraction, Re, area_ratio);
    return {r.zeta * b, r.dzeta_dRe * b};
}

double conical_diffuser_zeta(double full_angle_deg, double area_ratio) {
    check_area_ratio(area_ratio, "conical_diffuser_zeta");
    if (!(full_angle_deg >= 0.0)) {
        throw std::invalid_argument("conical_diffuser_zeta: angle must be >= 0");
    }
    return cone_phi(full_angle_deg, area_ratio) * (1.0 - area_ratio) * (1.0 - area_ratio);
}

double conical_nozzle_zeta(double full_angle_deg, double length_over_D, double area_ratio) {
    check_area_ratio(area_ratio, "conical_nozzle_zeta");
    if (!(full_angle_deg >= 0.0) || !(length_over_D >= 0.0)) {
        throw std::invalid_argument("conical_nozzle_zeta: angle and length must be >= 0");
    }
    return nozzle_prime(full_angle_deg, length_over_D) * (1.0 - area_ratio);
}

BordaCarnotStep compressible_borda_carnot(double M1, double area_ratio, double gamma) {
    check_area_ratio(area_ratio, "compressible_borda_carnot");
    if (!(M1 >= 0.0) || M1 > 1.0) {
        throw std::invalid_argument("compressible_borda_carnot: need 0 <= M1 <= 1");
    }
    if (!(gamma > 1.0)) {
        throw std::invalid_argument("compressible_borda_carnot: need gamma > 1");
    }
    using D2 = DualN<2>;
    const auto bc = borda_carnot_t(D2::seed(M1, 0), area_ratio, D2::seed(gamma, 1));
    const D2 pi = solver::dexp(bc.ln_pi);
    return {std::sqrt(bc.y.v), pi.v, pi.d[0], pi.d[1]};
}

AreaChangeLoss area_change_loss(double m_dot_abs, double Pt0, double Tt0,
                                const std::vector<double>& Y, double F_in, double F_out,
                                const std::string& model, double length, double D_h,
                                bool incompressible) {
    const Branch b = make_branch(F_in, F_out, model, length, D_h);
    AreaChangeLoss out;
    out.expansion = b.expansion;
    // Zero flow answers zero loss; NaN used to as well (#769).
    if (!validation::is_finite(m_dot_abs) || !validation::is_finite(Pt0) ||
        !validation::is_finite(Tt0)) {
        throw std::invalid_argument("area_change_loss: m_dot_abs, Pt0 and Tt0 must be finite");
    }
    if (!(m_dot_abs > 0.0)) {
        return out;
    }
    const Gas gas = make_gas(Tt0, Pt0, Y, false);
    const double G = m_dot_abs / b.F_small;
    double p = G;
    if (!incompressible) {
        // The closure's Mach from its station's flux, isentropic from (Pt0, Tt0).
        const double Ga = b.expansion ? G : G / b.Cc;
        auto flux = [&](double M) {
            return eval_compressible<double>(M, Pt0, Tt0, gas.g, gas.Rg, b, gas).G_small /
                   (b.expansion ? 1.0 : b.Cc);
        };
        if (flux(1.0) <= Ga) {
            p = 1.0;
            out.choked = true;
        } else {
            double lo = 0.0;
            double hi = 1.0;
            double M = std::min(Ga / flux(1.0), 0.9);
            for (int it = 0;; ++it) {
                using D1 = DualN<1>;
                const auto e = eval_compressible<D1>(D1::seed(M, 0), D1::constant(Pt0),
                                                     D1::constant(Tt0), D1::constant(gas.g),
                                                     D1::constant(gas.Rg), b, gas);
                const double scale = b.expansion ? 1.0 : b.Cc;
                const double f = e.G_small.v / scale - Ga;
                if (f > 0.0) {
                    hi = M;
                } else {
                    lo = M;
                }
                const double fp = e.G_small.d[0] / scale;
                double M_new = (fp > 0.0) ? M - f / fp : 0.5 * (lo + hi);
                if (!(M_new > lo && M_new < hi)) {
                    M_new = 0.5 * (lo + hi);
                }
                if (std::abs(M_new - M) <= 1e-15 * M) {
                    M = M_new;
                    break;
                }
                M = M_new;
                if (it + 1 == inner_solve::kAreaChangeMaxIter) {
                    if (inner_solve::cap_exit_is_root(hi - lo, M)) {
                        break;
                    }
                    inner_solve::throw_no_convergence(
                        "area_change_loss (Mach)", inner_solve::kAreaChangeMaxIter,
                        "bracket width", hi - lo);
                }
            }
            p = M;
        }
    }
    const auto e = eval_branch<double>(incompressible, p, Pt0, Tt0, gas.g, gas.Rg, b, gas);
    out.dPt = e.loss;
    out.mach = e.mach;
    out.zeta = (e.q_small > 0.0) ? e.loss / e.q_small : 0.0;
    // Re at the closure station's static T (the stagnation T in the
    // incompressible limit), as raw_flow's T_mu; diagnostic only, the loss
    // above does not read it.
    const double T_mu =
        incompressible ? Tt0 : Tt0 / (1.0 + 0.5 * (gas.g - 1.0) * e.mach * e.mach);
    out.Re = b.low_re ? e.Re : G * b.D / viscosity(T_mu, gas.P, gas.X);
    return out;
}

AreaChangeFlow area_change_flow(double Pt0, double Tt0, const std::vector<double>& Y,
                                double Pt_target, double F_in, double F_out,
                                const std::string& model, double length, double D_h,
                                bool incompressible, bool with_derivatives, bool exit_head_lost) {
    if (!(Pt0 > 0.0) || !(Tt0 > 0.0)) {
        throw std::invalid_argument("area_change_flow: Pt0 and Tt0 must be positive");
    }
    const Branch b = make_branch(F_in, F_out, model, length, D_h, exit_head_lost);
    const Gas gas = make_gas(Tt0, Pt0, Y, with_derivatives);
    // Odd in the drive: a negative drive runs the same branch backwards,
    // m(-d) = -m(d) at the same source state.
    const double drive_signed = Pt0 - Pt_target;
    const double sgn = (drive_signed < 0.0) ? -1.0 : 1.0;
    const double drive = std::abs(drive_signed);

    RawFlow r;
    if (drive >= b.floor) {
        r = raw_flow(drive, Pt0, Tt0, b, gas, incompressible, with_derivatives);
    } else {
        // Below the floor: m = m0 h(t), t = drive / floor, m0 the raw flow at
        // the floor and
        //   h(t) = t exp((sigma - 1)(t^2 - 1) / 2),
        // odd, h(1) = 1, h'(1) = sigma, h'(0) = exp((1 - sigma)/2) > 0, and
        // monotone for every sigma > 0. sigma is the incompressible law's
        // exponent d ln m / d ln(drive) at the floor: 1/2 for a quadratic loss
        // (contraction, cone), 1/e(Re) for a sudden expansion/contraction, e
        // the exponent of (zeta + exit_head) Re^2. It is analytic in the state
        // through Re, so every partial is exact; the compressible law's own
        // exponent differs by O(M^2) at the floor, the size of the slope's step
        // there.
        const RawFlow at = raw_flow(b.floor, Pt0, Tt0, b, gas, incompressible, true);
        double sigma = 0.5;
        double dsig_dRe = 0.0;
        if (b.low_re) {
            const double offset = b.exit_head_lost ? (b.r_exit * b.r_exit) / b.zeta_inc : 0.0;
            const auto [e, de] = table_exponent(b.table, at.Re, b.n, offset);
            sigma = 1.0 / e;
            dsig_dRe = -de / (e * e);
        }
        const double t = drive / b.floor;
        const double E = std::exp(0.5 * (sigma - 1.0) * (t * t - 1.0));
        const double h = t * E;
        const double dh = E * (1.0 + (sigma - 1.0) * t * t);
        const double k = at.m * 0.5 * h * (t * t - 1.0) * dsig_dRe;  // dm / dRe_floor
        r = at;
        r.m = at.m * h;
        r.choked = false;
        r.dm_ddrive = at.m * dh / b.floor;
        r.dm_dPt0 = h * at.dm_dPt0 + k * at.dRe_dPt0;
        r.dm_dTt0 = h * at.dm_dTt0 + k * at.dRe_dTt0;
        r.dm_dg = h * at.dm_dg + k * at.dRe_dg;
        r.dm_dR = h * at.dm_dR + k * at.dRe_dR;
        r.dm_dmu = h * at.dm_dmu + k * at.dRe_dmu;
        const AreaChangeLoss l = area_change_loss(r.m, Pt0, Tt0, Y, F_in, F_out, model, length,
                                                  D_h, incompressible);
        r.mach = l.mach;
        r.Re = l.Re;
        r.zeta = l.zeta;
    }

    AreaChangeFlow out;
    out.m_dot = sgn * r.m;
    out.choked = r.choked;
    out.mach = r.mach;
    out.Re = r.Re;
    out.zeta = r.zeta;
    if (!with_derivatives) {
        return out;
    }
    // (Pt0, drive) -> (Pt0, Pt_target): drive = Pt0 - Pt_target. Odd in the
    // drive: d/d(drive) is even, the rest odd.
    out.dm_dPt_target = -r.dm_ddrive;
    out.dm_dPt0 = sgn * r.dm_dPt0 + r.dm_ddrive;
    out.dm_dTt0 = sgn * r.dm_dTt0;
    // Composition: gamma and R through X, and mu at the closure's
    // temperature (mu_at's offset seed). dX/dY tangent to sum X = 1:
    //   dX_i/dY_k = (delta_ik - X_i) / (W_k S),  S = sum_j Y_j / W_j.
    const std::size_t ns = Y.size();
    std::vector<double> mu_X(ns, 0.0);
    if (r.dm_dmu != 0.0) {
        mu_X = viscosity_dX(r.T_mu, gas.P, gas.X);
    }
    std::vector<double> dm_dX(ns, 0.0);
    double x_dot = 0.0;
    for (std::size_t k = 0; k < ns; ++k) {
        dm_dX[k] = sgn * (r.dm_dg * gas.g_X[k] + r.dm_dR * gas.Rg_X[k] + r.dm_dmu * mu_X[k]);
        x_dot += gas.X[k] * dm_dX[k];
    }
    double S = 0.0;
    for (std::size_t k = 0; k < ns; ++k) {
        S += Y[k] / species_molar_mass(k);
    }
    out.dm_dY.assign(ns, 0.0);
    for (std::size_t k = 0; k < ns; ++k) {
        out.dm_dY[k] = (dm_dX[k] - x_dot) / (species_molar_mass(k) * S);
    }
    return out;
}

AreaChangeElementFlow area_change_element_flow(double Pt_a, double Tt_a,
                                               const std::vector<double>& Y_a, double Pt_b,
                                               double Tt_b, const std::vector<double>& Y_b,
                                               double F_a, double F_b, const std::string& model,
                                               double length, double D_h, bool incompressible,
                                               bool with_derivatives, bool exit_head_lost_b,
                                               bool exit_head_lost_a) {
    // The blend's width: the smaller of the two directions' drive floors.
    const double width =
        std::min(make_branch(F_a, F_b, model, length, D_h, exit_head_lost_b).floor,
                 make_branch(F_b, F_a, model, length, D_h, exit_head_lost_a).floor);
    const double s = (Pt_a - Pt_b) / width;
    double w = 1.0;
    double dw_ds = 0.0;
    if (s < -30.0) {
        w = 0.0;
    } else if (s <= 30.0) {
        w = 1.0 / (1.0 + std::exp(-s));
        dw_ds = w * (1.0 - w);
    }
    AreaChangeElementFlow out;
    out.dm_dY_a.assign(Y_a.size(), 0.0);
    out.dm_dY_b.assign(Y_b.size(), 0.0);
    AreaChangeFlow ab;
    AreaChangeFlow ba;
    if (w > 0.0) {
        ab = area_change_flow(Pt_a, Tt_a, Y_a, Pt_b, F_a, F_b, model, length, D_h, incompressible,
                              with_derivatives, exit_head_lost_b);
    }
    if (w < 1.0) {
        ba = area_change_flow(Pt_b, Tt_b, Y_b, Pt_a, F_b, F_a, model, length, D_h, incompressible,
                              with_derivatives, exit_head_lost_a);
    }
    // In the element's orientation the b -> a flow is negative.
    const double m_ab = ab.m_dot;
    const double m_ba = -ba.m_dot;
    out.m_dot = w * m_ab + (1.0 - w) * m_ba;
    const AreaChangeFlow& up = (w >= 0.5) ? ab : ba;
    out.choked = up.choked;
    out.mach = up.mach;
    out.Re = up.Re;
    out.zeta = up.zeta;
    if (!with_derivatives) {
        return out;
    }
    const double blend = dw_ds / width * (m_ab - m_ba);
    out.dm_dPt_a = blend;
    out.dm_dPt_b = -blend;
    if (w > 0.0) {
        out.dm_dPt_a += w * ab.dm_dPt0;
        out.dm_dPt_b += w * ab.dm_dPt_target;
        out.dm_dTt_a += w * ab.dm_dTt0;
        for (std::size_t k = 0; k < Y_a.size(); ++k) {
            out.dm_dY_a[k] += w * ab.dm_dY[k];
        }
    }
    if (w < 1.0) {
        out.dm_dPt_b -= (1.0 - w) * ba.dm_dPt0;
        out.dm_dPt_a -= (1.0 - w) * ba.dm_dPt_target;
        out.dm_dTt_b -= (1.0 - w) * ba.dm_dTt0;
        for (std::size_t k = 0; k < Y_b.size(); ++k) {
            out.dm_dY_b[k] -= (1.0 - w) * ba.dm_dY[k];
        }
    }
    return out;
}

// --- incompressible standalone kernels ------------------------------------------

namespace {

inline double smooth_abs(double m, double eps = 1e-12) { return std::sqrt(m * m + eps); }

// Logistic sigmoid (w, dw/dx), clamped against overflow.
inline std::pair<double, double> sigmoid(double x) {
    if (x > 30.0) return {1.0, 0.0};
    if (x < -30.0) return {0.0, 0.0};
    const double w = 1.0 / (1.0 + std::exp(-x));
    return {w, w * (1.0 - w)};
}

// zeta_fwd/zeta_rev (and their Re-slopes) blended by direction, then
// dP = zeta m |m| / (2 rho F_small^2).
AreaChangeResult blend_kernel(double m_dot, double rho, double mu, double F0, double F1,
                              double m_scale, double D_small, double zeta_exp, double dzeta_exp_dRe,
                              double zeta_con, double dzeta_con_dRe, double Re) {
    const double F_small = std::min(F0, F1);
    const double m_abs = smooth_abs(m_dot);
    const bool fwd_exp = F0 <= F1;
    const double zeta_fwd = fwd_exp ? zeta_exp : zeta_con;
    const double zeta_rev = fwd_exp ? zeta_con : zeta_exp;
    const auto [w, dw_dx] = sigmoid(m_dot / m_scale);
    const double zeta_eff = w * zeta_fwd + (1.0 - w) * zeta_rev;
    const double K = 1.0 / (2.0 * rho * F_small * F_small + kEps);
    const double dP = K * zeta_eff * m_dot * m_abs;
    const double dm_abs_dm = m_dot / m_abs;
    const double dRe_dm = D_small / (mu * F_small + kEps) * dm_abs_dm;
    const double w_exp = fwd_exp ? w : 1.0 - w;
    // d zeta_eff / dRe at fixed direction weight.
    const double dzeta_dRe = w_exp * dzeta_exp_dRe + (1.0 - w_exp) * dzeta_con_dRe;
    const double dzeta_dm = (zeta_fwd - zeta_rev) * dw_dx / m_scale + dzeta_dRe * dRe_dm;
    const double dS_dm = K * (zeta_eff * (m_abs + m_dot * dm_abs_dm) + dzeta_dm * m_dot * m_abs);
    const double dS_drho = (std::abs(rho) > kEps) ? -dP / rho : 0.0;
    const double dRe_dmu = (std::abs(mu) > kEps) ? -Re / mu : 0.0;
    const double dS_dmu = K * dzeta_dRe * dRe_dmu * m_dot * m_abs;
    return {dP, dS_dm, dS_drho, dS_dmu};
}

}  // namespace

AreaChangeResult sharp_area_change(double m_dot, double rho, double mu, double F0, double F1,
                                   double m_scale, double D_h) {
    if (!(F0 > 0.0) || !(F1 > 0.0)) {
        throw std::invalid_argument("sharp_area_change: areas must be positive");
    }
    const double F_small = std::min(F0, F1);
    const double n = F_small / std::max(std::max(F0, F1), kEps);
    const double D_small = (D_h > 0.0) ? D_h : std::sqrt(std::max(4.0 * F_small / M_PI, 0.0));
    const double Re = smooth_abs(m_dot) * D_small / (mu * F_small + kEps);
    const ZetaAndSlope e = sudden_expansion_zeta(Re, n);
    const ZetaAndSlope c = sudden_contraction_zeta(Re, n);
    return blend_kernel(m_dot, rho, mu, F0, F1, m_scale, D_small, e.zeta, e.dzeta_dRe, c.zeta,
                        c.dzeta_dRe, Re);
}

AreaChangeResult conical_area_change(double m_dot, double rho, double mu, double F0, double F1,
                                     double length, double m_scale) {
    if (!(F0 > 0.0) || !(F1 > 0.0)) {
        throw std::invalid_argument("conical_area_change: areas must be positive");
    }
    const double F_small = std::min(F0, F1);
    const double F_large = std::max(F0, F1);
    const double n = F_small / std::max(F_large, kEps);
    const double dr = std::sqrt(F_large / M_PI) - std::sqrt(F_small / M_PI);
    const double theta = std::atan2(dr, std::max(length, kEps));
    const double zeta_exp = conical_diffuser_zeta(2.0 * theta * 180.0 / M_PI, n);
    const double D_small = std::sqrt(4.0 * F_small / M_PI);
    const double zeta_con =
        conical_nozzle_zeta(2.0 * theta * 180.0 / M_PI, std::max(length, kEps) / D_small, n);
    return blend_kernel(m_dot, rho, mu, F0, F1, m_scale, 1.0, zeta_exp, 0.0, zeta_con, 0.0, 0.0);
}

}  // namespace combaero
