#include "../include/dual_number.h"
#include "../include/correlation_status.h"
#include "../include/friction.h"
#include "../include/humidair.h"
#include "../include/math_constants.h"
#include "../include/orifice.h"
#include "../include/orifice_jobson.h"
#include "../include/input_validation.h"
#include "../include/inner_solve.h"
#include "../include/thermo.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// -------------------------------------------------------------
// OrificeGeometry implementation
// -------------------------------------------------------------

double OrificeGeometry::beta() const {
    if (!combaero::validation::is_finite_positive(D)) {
        throw std::invalid_argument("OrificeGeometry: pipe diameter D must be > 0");
    }
    return d / D;
}

double OrificeGeometry::area() const {
    return M_PI * d * d / 4.0;
}

double OrificeGeometry::t_over_d() const {
    if (d <= 0.0) return 0.0;
    return t / d;
}

double OrificeGeometry::r_over_d() const {
    if (d <= 0.0) return 0.0;
    return r / d;
}

bool OrificeGeometry::is_valid() const {
    if (d <= 0.0 || D <= 0.0) return false;
    if (d >= D) return false;
    if (t < 0.0 || r < 0.0) return false;
    return true;
}

// -------------------------------------------------------------
// OrificeState implementation
// -------------------------------------------------------------

double OrificeState::Re_d(double beta) const {
    // Orifice Reynolds number based on orifice diameter d
    // From continuity: v_orifice = v_pipe / beta^2
    // Re_d = (rho * v_orifice * d) / mu
    //      = (rho * (v_pipe / beta^2) * (D * beta)) / mu
    //      = (rho * v_pipe * D) / mu * (1/beta)
    //      = Re_D / beta
    // However, the conventional definition uses Re_d = Re_D * beta
    // (based on the diameter ratio, not the actual flow velocity)
    return Re_D * beta;
}

// -------------------------------------------------------------
// Namespace: individual correlations
// -------------------------------------------------------------

namespace orifice {

// Clamping helpers for numerical stability near d -> D
inline double clamp_beta(double b) { return std::min(b, 0.98); }
inline double clamp_Cd(double c) { return std::min(c, 1.5); }

std::tuple<double, double> tap_positions(OrificeTap tap, double D) {
    switch (tap) {
        case OrificeTap::Corner:
            return {0.0, 0.0};
        case OrificeTap::D_D2:
            return {taps::dd2_L1, taps::dd2_L2};
        case OrificeTap::Flange: {
            if (!(D > 0.0)) {
                throw std::invalid_argument(
                    "Flange taps need the pipe diameter: L1 = 25.4 mm / D");
            }
            const double L = taps::inch_mm / (D * 1000.0);
            return {L, L};
        }
    }
    throw std::invalid_argument("unknown OrificeTap");
}

namespace {

// Reader-Harris/Gallagher, ISO 5167-2:2003 Eq. (4), and its exact d/dRe_D,
// as printed: no hold and no cap. Only evaluated at Re >= its lower bound.
std::tuple<double, double> rhg_raw_and_dRe(double beta, double Re, double D,
                                           OrificeTap tap) {
    namespace R = reader_harris;

    // beta is held at 0.98 for stability near d -> D; ISO's range ends at
    // 0.75 (metering_cd_in_range says so).
    const double b = clamp_beta(beta);
    const double b2 = b * b;
    const double b4 = b2 * b2;
    const double b8 = b4 * b4;

    const auto [L1, L2] = tap_positions(tap, D);
    const double M2 = 2.0 * L2 / (1.0 - b);

    const double A = std::pow(R::A_coef * b / Re, R::A_exp);
    const double dA = -R::A_exp * A / Re;

    // Re terms: 0.000521 (1e6 b / Re)^0.7 + (0.0188 + 0.0063 A) b^3.5 (1e6 / Re)^0.3
    const double t3 = R::C3 * std::pow(R::re_scale * b / Re, R::C3_exp);
    const double X = std::pow(R::re_scale / Re, R::C4_re);
    const double dX = -R::C4_re * X / Re;
    const double b35 = std::pow(b, R::C4_beta);

    // Upstream tap: (0.043 + 0.080 e^-10L1 - 0.123 e^-7L1)(1 - 0.11 A) b^4/(1 - b^4)
    const double up = (R::C5a + R::C5b * std::exp(-R::C5_exp1 * L1) -
                       R::C5c * std::exp(-R::C5_exp2 * L1)) *
                      b4 / (1.0 - b4);

    double C = R::C0 + R::C1 * b2 - R::C2 * b8 + t3 +
               (R::C4a + R::C4b * A) * b35 * X + up * (1.0 - R::C5_A_coef * A) -
               R::C6 * (M2 - R::C6_M_coef * std::pow(M2, R::C6_M_exp)) *
                   std::pow(b, R::C6_beta);

    // Small pipes, D < 71.12 mm: + 0.011 (0.75 - beta)(2.8 - D/25.4).
    const double D_mm = D * 1000.0;
    if (D_mm > 0.0 && D_mm < R::small_pipe_D_mm) {
        C += R::small_pipe_coef * (R::small_pipe_beta - b) *
             (R::small_pipe_D_in - D_mm / taps::inch_mm);
    }

    const double dC = -R::C3_exp * t3 / Re + b35 * (R::C4b * dA * X + (R::C4a + R::C4b * A) * dX) -
                      up * R::C5_A_coef * dA;
    return {C, dC};
}

// Lower Re_D bound of Reader-Harris/Gallagher, ISO 5167-2:2003 5.3.1.
double rhg_re_min(double beta, double D, OrificeTap tap) {
    namespace R = reader_harris;
    if (tap == OrificeTap::Flange) {
        return std::max(R::range_re_min, R::range_re_flange * beta * beta * D * 1000.0);
    }
    return beta <= R::range_beta_split ? R::range_re_min : R::range_re_beta2 * beta * beta;
}

// The #873 hold: the correlation at Re_D >= Re_b (untouched, so bit-identical
// there), its slope ramped linearly in ln Re_D to 0 over w = blend_ln_re
// below Re_b, and held below that. The form and its derivative are in
// orifice.h at metering_cd_and_dRe. raw(Re) is (C, dC/dRe) as printed.
template <class Raw>
std::tuple<double, double> held_below_bound(const Raw& raw, double Re_D, double Re_b) {
    if (Re_D >= Re_b) {
        return raw(Re_D);
    }
    const auto [C_b, dC_b] = raw(Re_b);
    const double s = dC_b * Re_b;  // dC/d(ln Re_D) at the bound
    const double w = metering_hold::blend_ln_re;
    // Also catches Re_D <= 0 (no flow, reverse-flow iterates). A NaN Re_D
    // falls through to the blend and stays NaN.
    if (Re_D <= Re_b * std::exp(-w)) {
        return {C_b - 0.5 * s * w, 0.0};
    }
    const double t = 1.0 + std::log(Re_D / Re_b) / w;
    return {C_b + 0.5 * s * w * (t * t - 1.0), s * t / Re_D};
}

std::tuple<double, double> rhg_and_dRe(double beta, double Re_D, double D, OrificeTap tap) {
    const auto [C, dC] = held_below_bound(
        [&](double Re) { return rhg_raw_and_dRe(beta, Re, D, tap); }, Re_D,
        rhg_re_min(beta, D, tap));
    // The 1.5 cap (a numerical bound for lossless bypass recovery) stays for
    // geometry far outside ISO's: with the hold it is out of reach for corner
    // and D-D/2 taps (the network's), but flange taps in a pipe below ~32 mm
    // sit many diameters from the plate, where Eq. (4)'s downstream term
    // grows without bound at ANY Re_D (#873 scan).
    if (C >= 1.5) {
        return {1.5, 0.0};
    }
    return {C, dC};
}

} // namespace

double Cd_ReaderHarrisGallagher(double beta, double Re_D, double D, OrificeTap tap) {
    return std::get<0>(rhg_and_dRe(beta, Re_D, D, tap));
}

// Convert between Cd and loss coefficient K
// -------------------------------------------------------------
// McGreehan and Schotsch (1988)
// -------------------------------------------------------------

namespace mcgreehan_schotsch {

// Below its stated floor Eq. (8) does not merely lose accuracy, it diverges:
// 0.5885 + 372/Re reaches Cd = 1.0 at Re = 904 and grows without bound below
// that. So Re is held AT the validity floor rather than merely kept positive.
// The result is constant for Re < re_min, which is a documented "held at the
// edge of validity" extrapolation and not a claim about the laminar regime --
// that regime belongs to Wu, Burton and Schoenau (2002), deliberately not
// implemented (item I5 / decision D7 of the extraction).
double reynolds_baseline(double Re) {
    return re_c0 + re_c1 / std::max(Re, re_min);
}

double nozzle_baseline(double Re) {
    return nozzle_c0 - nozzle_c1 / std::sqrt(std::max(Re, re_min));
}

double corner_factor(double r_over_d) {
    const double rd = std::max(r_over_d, 0.0);
    return corner_floor +
           corner_coef * std::exp(-corner_exp1 * rd - corner_exp2 * rd * rd);
}

double length_factor(double L_over_d) {
    const double ld = std::max(L_over_d, 0.0);
    return (1.0 + length_coef * std::exp(-length_decay * ld * ld)) *
           (length_a + length_b * ld);
}

double cd_with_corner(double Re, double r_over_d) {
    // Eq. (11)
    return 1.0 - corner_factor(r_over_d) * (1.0 - reynolds_baseline(Re));
}

double cd_with_corner_and_length(double Re, double r_over_d, double L_over_d) {
    const double rd = std::max(r_over_d, 0.0);
    double basic   = cd_with_corner(Re, rd);
    double ld      = std::max(L_over_d, 0.0);

    if (rd > 0.0) {
        // Eq. (15): a revised basic Cd, with g evaluated at r/d in place of
        // L/d. The corner radius suppresses inlet separation and so removes
        // part of the long orifice's dynamic-pressure-recovery benefit.
        basic = 1.0 - length_factor(rd) * (1.0 - basic);
        // Eq. (16): the inlet radius is subtracted from the flat length.
        ld = std::max(ld - rd, 0.0);
    }

    // Eq. (13)
    return 1.0 - length_factor(ld) * (1.0 - basic);
}

double cd_with_crossflow(double cd_base, double U1_over_Vi, double eps) {
    double u = std::max(U1_over_Vi, 0.0);

    // eps = 0 recovers Eq. (17) exactly. It is an identity at zero crossflow
    // (C1 = 1, C2 = 0), so returning early keeps that exact and avoids
    // pow(0, fractional) entirely.
    if (eps <= 0.0) {
        if (u == 0.0) return cd_base;
    } else {
        // Regularised: bounds dCd/d(U1/Vi) at 1.451 instead of 1.4e4, for
        // 5.5e-5 on Cd. See rv_smooth_eps for the plateau this sits on.
        u = std::sqrt(u * u + eps * eps);
    }

    const double cd_ratio = cd_base / cd_reference;
    const double Rv       = u * std::pow(cd_ratio, rv_cd_exp);

    const double C1 = std::exp(-std::pow(Rv, c1_exp));
    const double C2 = c2_coef * std::pow(Rv, c2_exp) * std::pow(cd_ratio, c2_cd_exp);
    const double C3 = std::exp(-c3_coef * std::pow(Rv, c3_exp));

    return cd_base * (C1 + C2 * C3);
}

std::tuple<double, double, double> cd_and_derivatives(double Re, double r_over_d,
                                                      double L_over_d,
                                                      double U1_over_Vi) {
    using D = combaero::solver::DualN<2>;  // partial 0 = Re, 1 = U1/Vi

    // Seeding max(Re, re_min) rather than Re IS the Jacobian-only floor
    // continuation. Below the floor the value is the floored one -- correct,
    // Eq. (8) is invalid and divergent there -- while the derivative comes
    // out as the live slope AT re_min instead of the true zero, so a Newton
    // step that wanders below has something to climb back on and dCd/dRe is
    // continuous across re_min. No reported Cd changes.
    const D re = D::seed(std::max(Re, re_min), 0);
    const D rb = re_c0 + re_c1 / re;  // Eq. (8)

    // r/d and L/d are geometry, never solver unknowns, so every factor built
    // from them is an ordinary double and contributes no partials.
    const double rd = std::max(r_over_d, 0.0);
    double ld       = std::max(L_over_d, 0.0);

    D basic = 1.0 - corner_factor(rd) * (1.0 - rb);  // Eq. (11)
    if (rd > 0.0) {
        basic = 1.0 - length_factor(rd) * (1.0 - basic);  // Eq. (15)
        ld    = std::max(ld - rd, 0.0);                   // Eq. (16)
    }
    const D cd_base = 1.0 - length_factor(ld) * (1.0 - basic);  // Eq. (13)

    // Eq. (17). The regularisation also removes the pow(0, fractional) that
    // the scalar path avoids with an early return: u >= rv_smooth_eps > 0
    // always, so Rv never reaches zero and every dpow below is finite.
    const D u_raw = D::seed(std::max(U1_over_Vi, 0.0), 1);
    const D u     = combaero::solver::dsqrt(u_raw * u_raw + rv_smooth_eps * rv_smooth_eps);

    const D cd_ratio = cd_base / cd_reference;
    const D Rv       = u * combaero::solver::dpow(cd_ratio, rv_cd_exp);
    const D C1       = combaero::solver::dexp(-combaero::solver::dpow(Rv, c1_exp));
    const D C2       = c2_coef * combaero::solver::dpow(Rv, c2_exp) * combaero::solver::dpow(cd_ratio, c2_cd_exp);
    const D C3       = combaero::solver::dexp(-c3_coef * combaero::solver::dpow(Rv, c3_exp));
    const D cd_out   = cd_base * (C1 + C2 * C3);

    return {cd_out.v, cd_out.d[0], cd_out.d[1]};
}

double cd(double Re, double r_over_d, double L_over_d, double U1_over_Vi,
          double eps) {
    return cd_with_crossflow(cd_with_corner_and_length(Re, r_over_d, L_over_d),
                             U1_over_Vi, eps);
}

// -------------------------------------------------------------
// Expansion factor Y, Eqs. (4)-(7)
// -------------------------------------------------------------

namespace {

// Smooth min/max. At eps = 0 these reduce EXACTLY to std::min / std::max, so
// the paper's hard clamp is the eps -> 0 limit rather than a separate branch.
// Away from the knee the deviation falls off quadratically; at the knee it is
// eps/2.
double soft_min(double x, double limit, double eps) {
    const double d = x - limit;
    return 0.5 * (x + limit - std::sqrt(d * d + eps * eps));
}

double soft_max(double x, double limit, double eps) {
    const double d = x - limit;
    return 0.5 * (x + limit + std::sqrt(d * d + eps * eps));
}

} // namespace

double critical_pressure_ratio(double gamma) {
    return std::pow(2.0 / (gamma + 1.0), gamma / (gamma - 1.0));
}

double expansion_orifice(double S, double gamma) {
    if (!combaero::validation::is_finite(S) || !combaero::validation::is_finite(gamma)) {
        throw std::invalid_argument("expansion_orifice: S and gamma must be finite");
    }
    // Eq. (4). Written as (1 - S)/gamma, identical to the printed
    // (P_t1 - P_s2)/(gamma P_t1).
    return 1.0 - y_orifice_coef * (1.0 - S) / gamma;
}

double expansion_nozzle(double S, double gamma) {
    // NaN S used to fall through to the 0 branch below (#769).
    if (!combaero::validation::is_finite(S) || !combaero::validation::is_finite(gamma)) {
        throw std::invalid_argument("expansion_nozzle: S and gamma must be finite");
    }
    // Eq. (5). The bracket is 0/0 at S = 1; the limit of
    // (1 - S^((g-1)/g))/(1 - S) there is (g-1)/g, which makes Y_n(1) = 1.
    if (S >= 1.0) return 1.0;
    const double e = (gamma - 1.0) / gamma;
    const double one_minus_S = 1.0 - S;
    if (one_minus_S < 1.0e-9) return 1.0;
    const double ratio = (1.0 - std::pow(S, e)) / one_minus_S;
    const double val = std::pow(S, 2.0 / gamma) * (gamma / (gamma - 1.0)) * ratio;
    return val > 0.0 ? std::sqrt(val) : 0.0;
}

double expansion_blend_weight(double cd, double eps) {
    const double x = x_blend_slope * (cd - x_blend_cd0);
    // Saturate into [0, 1]. The paper prints neither bound: it states the
    // blend applies "for Cd > 0.82" and X reaches 1 at Cd = 0.94, but gives
    // no clamp above. Unclamped, Cd = 0.99 would give X = 1.45 and
    // Y = -0.45 Y_o + 1.45 Y_n, extrapolating past a nozzle. Clamping is a
    // decision the source did not make; see decision D10 of the extraction.
    return soft_max(soft_min(x, 1.0, eps), 0.0, eps);
}

double expansion_factor(double cd, double S, double gamma, double eps) {
    // gamma <= 1 still means "no expansion correction"; NaN meant it too (#769).
    if (!combaero::validation::is_finite(cd) || !combaero::validation::is_finite(S) ||
        !combaero::validation::is_finite(gamma) || !combaero::validation::is_finite(eps)) {
        throw std::invalid_argument("expansion_factor: cd, S, gamma and eps must be finite");
    }
    if (!(gamma > 1.0)) return 1.0;

    // Below the critical pressure ratio the isentropic form turns over and
    // predicts DECREASING flow -- the classic non-physical branch. Saturate S
    // at the choke point, smoothly for the same reason the blend weight is
    // smoothed: a hard clamp here would zero dY/dS discontinuously.
    // Both bounds are SOFT. An earlier revision smoothed the choke point but
    // left a hard std::min(S, 1.0) at the top, which a scan for zero-derivative
    // regions and derivative discontinuities caught: it put a kink of ~0.455 in
    // dY/dS at S = 1 and a floor above it. S > 1 is reverse flow, which a
    // solver iterate can reach, so it gets the same treatment as the rest.
    const double S_eff = soft_max(soft_min(S, 1.0, s_smooth_eps),
                                  critical_pressure_ratio(gamma),
                                  s_smooth_eps);

    const double X = expansion_blend_weight(cd, eps);
    return (1.0 - X) * expansion_orifice(S_eff, gamma)
           + X * expansion_nozzle(S_eff, gamma);
}

} // namespace mcgreehan_schotsch

double Cd_McGreehanSchotsch(double Re, double r_over_d, double L_over_d,
                            double U1_over_Vi, double eps) {
    return mcgreehan_schotsch::cd(Re, r_over_d, L_over_d, U1_over_Vi, eps);
}

double K_from_Cd(double Cd, double beta) {
    if (!combaero::validation::is_finite_positive(Cd) || Cd > 1.5) {
        throw std::invalid_argument("Cd must be in (0, 1.5]");
    }
    const double b = clamp_beta(beta);
    const double beta4 = std::pow(b, 4.0);
    return (1.0 / (Cd * Cd) - 1.0) * (1.0 - beta4);
}

double Cd_from_K(double K, double beta) {
    if (!combaero::validation::is_finite_nonnegative(K)) {
        throw std::invalid_argument("K must be >= 0");
    }
    const double b = clamp_beta(beta);
    const double beta4 = std::pow(b, 4.0);
    return clamp_Cd(1.0 / std::sqrt(1.0 + K / (1.0 - beta4)));
}

} // namespace orifice

// -------------------------------------------------------------
// Main Cd functions (free functions)
// -------------------------------------------------------------

std::tuple<double, double> Cd_sharp_thin_plate_and_dRe(const OrificeGeometry& geom,
                                                       const OrificeState& state) {
    return metering_cd_and_dRe(MeteringCdCorrelation::ReaderHarrisGallagher, geom, state.Re_D);
}

double Cd_sharp_thin_plate(const OrificeGeometry& geom, const OrificeState& state) {
    return std::get<0>(Cd_sharp_thin_plate_and_dRe(geom, state));
}

std::tuple<double, double> metering_cd_and_dRe(MeteringCdCorrelation id,
                                               const OrificeGeometry& geom,
                                               double Re_D) {
    if (!geom.is_valid()) {
        throw std::invalid_argument("Invalid orifice geometry");
    }
    switch (id) {
        case MeteringCdCorrelation::ReaderHarrisGallagher:
            return orifice::rhg_and_dRe(geom.beta(), Re_D, geom.D, geom.tap);
        case MeteringCdCorrelation::Constant:
        case MeteringCdCorrelation::UserFunction:
            break;
    }
    throw std::invalid_argument(
        "metering_cd_and_dRe: only ReaderHarrisGallagher carries a "
        "formula; a Constant or UserFunction Cd belongs to the caller");
}

double metering_cd_re_min(MeteringCdCorrelation id, const OrificeGeometry& geom) {
    if (!geom.is_valid()) {
        throw std::invalid_argument("Invalid orifice geometry");
    }
    switch (id) {
        case MeteringCdCorrelation::ReaderHarrisGallagher:
            return orifice::rhg_re_min(geom.beta(), geom.D, geom.tap);
        case MeteringCdCorrelation::Constant:
        case MeteringCdCorrelation::UserFunction:
            break;
    }
    throw std::invalid_argument(
        "metering_cd_re_min: only ReaderHarrisGallagher carries a Reynolds "
        "range; a Constant or UserFunction Cd belongs to the caller");
}

bool metering_cd_in_range(MeteringCdCorrelation id, const OrificeGeometry& geom,
                          double Re_D) {
    namespace R = orifice::reader_harris;
    if (id != MeteringCdCorrelation::ReaderHarrisGallagher) {
        return true;  // Constant / UserFunction: the caller's Cd
    }
    if (!geom.is_valid()) {
        return false;
    }
    // The bound the hold starts at (#873), defined once.
    if (!(Re_D >= metering_cd_re_min(id, geom))) {
        return false;
    }
    const double d_mm = geom.d * 1000.0;
    const double D_mm = geom.D * 1000.0;
    const double b = geom.beta();
    return !(d_mm < R::range_d_min_mm || D_mm < R::range_D_min_mm || D_mm > R::range_D_max_mm ||
             b < R::range_beta_min || b > R::range_beta_max);
}

std::tuple<double, double> permanent_loss_ratio_and_dC(double beta, double C) {
    if (!(beta > 0.0 && beta < 1.0) || !(C > 0.0)) {
        throw std::invalid_argument(
            "permanent_loss_ratio_and_dC: needs 0 < beta < 1 and C > 0");
    }
    const double b2 = beta * beta;
    const double b4 = b2 * b2;
    const double s = std::sqrt(1.0 - b4 * (1.0 - C * C));
    const double a = C * b2;
    const double k = (s - a) / (s + a);
    // dk/dC = 2 (a s' - s a') / (s + a)^2, with s' = beta^4 C / s, a' = beta^2.
    const double ds = b4 * C / s;
    const double dk = 2.0 * (a * ds - s * b2) / ((s + a) * (s + a));
    return {k, dk};
}

double velocity_of_approach_factor(double beta) {
    // A finite beta outside (0, 1) means "no correction"; NaN used to be
    // read the same way (#769).
    if (!combaero::validation::is_finite(beta)) {
        throw std::invalid_argument("velocity_of_approach_factor: beta must be finite");
    }
    if (!(beta > 0.0 && beta < 1.0)) {
        return 1.0;
    }
    const double b2 = beta * beta;
    return 1.0 / std::sqrt(1.0 - b2 * b2);
}

// -------------------------------------------------------------
// Correlation classes
// -------------------------------------------------------------

namespace {

// An explicit inclination opt-in on a correlation it does not compose with
// is refused, never ignored: only Idelchik1966Thick has the source's
// inlet-edge term (Eq. 4-20's zeta') to put diagram 3-2 into.
void refuse_inclination(const DischargeHoleGeometry& hole, const std::string& who) {
    if (hole.is_inclined()) {
        throw std::invalid_argument(
            who + ": the hole-inclination term (Idelchik diagram 3-2) composes "
            "only with Idelchik1966Thick, whose zeta' carries a sharp normal "
            "inlet it can replace. Use Idelchik1966Thick, or leave "
            "inclination_deg at 90.");
    }
}

class ReaderHarrisGallagherCorrelation : public OrificeCorrelationBase {
public:
    double Cd(const OrificeGeometry& geom, const OrificeState& state) const override {
        return Cd_sharp_thin_plate(geom, state);
    }
    std::string name() const override { return "Reader-Harris/Gallagher (ISO 5167-2)"; }
};

class ConstantCdCorrelation : public OrificeCorrelationBase {
    double Cd_value_;
public:
    explicit ConstantCdCorrelation(double Cd = orifice::defaults::metering_cd)
        : Cd_value_(Cd) {}
    double Cd(const OrificeGeometry&, const OrificeState&) const override {
        return Cd_value_;
    }
    std::string name() const override { return "Constant Cd"; }
};

class UserFunctionCorrelation : public OrificeCorrelationBase {
    CdFunction fn_;
    std::string name_;
public:
    UserFunctionCorrelation(CdFunction fn, std::string name)
        : fn_(std::move(fn)), name_(std::move(name)) {}
    double Cd(const OrificeGeometry& geom, const OrificeState& state) const override {
        return fn_(geom, state);
    }
    std::string name() const override { return name_; }
};

class TabulatedCorrelation : public OrificeCorrelationBase {
    std::vector<double> beta_values_;
    std::vector<double> Re_values_;
    std::vector<std::vector<double>> Cd_table_;
    std::string name_;

public:
    TabulatedCorrelation(std::vector<double> beta_values,
                         std::vector<double> Re_values,
                         std::vector<std::vector<double>> Cd_table,
                         std::string name)
        : beta_values_(std::move(beta_values))
        , Re_values_(std::move(Re_values))
        , Cd_table_(std::move(Cd_table))
        , name_(std::move(name)) {}

    double Cd(const OrificeGeometry& geom, const OrificeState& state) const override {
        return interpolate(geom.beta(), state.Re_D);
    }

    std::string name() const override { return name_; }

private:
    double interpolate(double beta, double Re_D) const {
        // Bilinear interpolation in beta and log(Re_D)
        const double log_Re = std::log10(std::max(Re_D, 1.0));

        // Find beta indices
        auto it_beta = std::lower_bound(beta_values_.begin(), beta_values_.end(), beta);
        std::size_t i_beta = (it_beta == beta_values_.begin()) ? 0 :
                             (it_beta == beta_values_.end()) ? beta_values_.size() - 2 :
                             static_cast<std::size_t>(it_beta - beta_values_.begin() - 1);

        // Find Re indices (in log space)
        std::vector<double> log_Re_values;
        log_Re_values.reserve(Re_values_.size());
        for (double Re : Re_values_) {
            log_Re_values.push_back(std::log10(std::max(Re, 1.0)));
        }
        auto it_Re = std::lower_bound(log_Re_values.begin(), log_Re_values.end(), log_Re);
        std::size_t i_Re = (it_Re == log_Re_values.begin()) ? 0 :
                           (it_Re == log_Re_values.end()) ? log_Re_values.size() - 2 :
                           static_cast<std::size_t>(it_Re - log_Re_values.begin() - 1);

        // Clamp indices
        i_beta = std::min(i_beta, beta_values_.size() - 2);
        i_Re = std::min(i_Re, Re_values_.size() - 2);

        // Interpolation weights
        const double t_beta = (beta - beta_values_[i_beta]) /
                              (beta_values_[i_beta + 1] - beta_values_[i_beta]);
        const double t_Re = (log_Re - log_Re_values[i_Re]) /
                            (log_Re_values[i_Re + 1] - log_Re_values[i_Re]);

        // Clamp weights to [0, 1]
        const double tb = std::max(0.0, std::min(t_beta, 1.0));
        const double tr = std::max(0.0, std::min(t_Re, 1.0));

        // Bilinear interpolation
        const double c00 = Cd_table_[i_beta][i_Re];
        const double c10 = Cd_table_[i_beta + 1][i_Re];
        const double c01 = Cd_table_[i_beta][i_Re + 1];
        const double c11 = Cd_table_[i_beta + 1][i_Re + 1];

        return (1 - tb) * (1 - tr) * c00
             + tb * (1 - tr) * c10
             + (1 - tb) * tr * c01
             + tb * tr * c11;
    }
};

// -------------------------------------------------------------
// Discharge-hole correlations
// -------------------------------------------------------------

// Monotone cubic (Fritsch-Carlson PCHIP) interpolation, returning the value
// AND its exact derivative.
//
// WHY NOT LINEAR. Every Idelchik curve below is a solver input, and linear
// interpolation is only C0: dCd/dRe would jump at each of the 14 table knots.
// That is the hazard class we regularised out of the McGreehan chain, and
// re-introducing 14 of them would be worse than the one we removed.
//
// WHY NOT A NATURAL SPLINE. These curves are monotone and a natural cubic
// overshoots near the flat tails -- zeta_thick is 1.58, 1.55, 1.55 over its
// last three knots, where an unlimited spline dips below the asymptote and
// invents a Cd above the source's. Fritsch-Carlson cannot overshoot, so the
// interpolant stays inside the tabulated envelope by construction.
//
// The derivative is exact, not a difference: the Hermite form is
// differentiated analytically, so it satisfies the (f, J) rule.
//
// Outside the table the value is held at the endpoint and the derivative is
// zero. Every table here ends on a flat tail, so that is the source's own
// behaviour rather than an extrapolation, except at the low-Re end -- see
// IdelchikWallCorrelation::clamped_re.
struct Interpolated {
    double value;
    double derivative;
};

Interpolated pchip(const double* xs, const double* ys, int n, double x) {
    if (n < 2) {
        return {ys[0], 0.0};
    }
    if (x <= xs[0]) {
        return {ys[0], 0.0};
    }
    if (x >= xs[n - 1]) {
        return {ys[n - 1], 0.0};
    }

    // Secant slopes
    std::vector<double> h(n - 1);
    std::vector<double> delta(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        h[i] = xs[i + 1] - xs[i];
        delta[i] = (ys[i + 1] - ys[i]) / h[i];
    }

    // Fritsch-Carlson tangents: zero at every sign change or flat span, and a
    // weighted harmonic mean elsewhere, which is what bounds the overshoot.
    std::vector<double> m(n, 0.0);
    for (int i = 1; i < n - 1; ++i) {
        if (delta[i - 1] * delta[i] > 0.0) {
            const double w1 = 2.0 * h[i] + h[i - 1];
            const double w2 = h[i] + 2.0 * h[i - 1];
            m[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i]);
        }
    }
    // One-sided ends, limited so they cannot introduce a new extremum.
    m[0] = delta[0];
    m[n - 1] = delta[n - 2];

    // Locate the interval
    int k = 0;
    while (k < n - 2 && x > xs[k + 1]) {
        ++k;
    }

    const double s = (x - xs[k]) / h[k];
    const double s2 = s * s;
    const double s3 = s2 * s;

    // Cubic Hermite basis and its derivative in s
    const double h00 = 2.0 * s3 - 3.0 * s2 + 1.0;
    const double h10 = s3 - 2.0 * s2 + s;
    const double h01 = -2.0 * s3 + 3.0 * s2;
    const double h11 = s3 - s2;

    const double d00 = 6.0 * s2 - 6.0 * s;
    const double d10 = 3.0 * s2 - 4.0 * s + 1.0;
    const double d01 = -6.0 * s2 + 6.0 * s;
    const double d11 = 3.0 * s2 - 2.0 * s;

    const double value = h00 * ys[k] + h10 * h[k] * m[k]
                       + h01 * ys[k + 1] + h11 * h[k] * m[k + 1];
    const double dvalue_ds = d00 * ys[k] + d10 * h[k] * m[k]
                           + d01 * ys[k + 1] + d11 * h[k] * m[k + 1];

    return {value, dvalue_ds / h[k]};
}

// Darcy friction factor and dlam/dRe for Idelchik's lam*l/Dh bore-friction
// term (his diagrams 2-2..2-5).
//
// TWO BRANCHES, both exact, blended because a hard switch would be a KINK in
// a solver input:
//   - laminar, Re < 2300: lam = 64/Re, Hagen-Poiseuille, not a correlation;
//   - turbulent, Re > 4000: Haaland's explicit form, using friction.h's own
//     constants (Define Once) rather than a second copy of them.
// Between them a smoothstep in log(Re) carries one to the other. The blend is
// a NUMERICAL DEVICE, not physics: inside 2300 < Re < 4000 neither branch is
// the source's, and the transition regime has no correlation in Idelchik
// either. It moves zeta by at most 0.9% at l/Dh = 2 and nothing at all
// outside that window.
Interpolated darcy_friction(double Re, double e_over_d) {
    const double re = std::max(Re, 1.0);

    const double lam_lam = 64.0 / re;
    const double dlam_lam = -64.0 / (re * re);

    // Haaland: 1/sqrt(f) = -1.8 log10[(e_D/3.7)^1.11 + 6.9/Re]
    const double rough = std::pow(e_over_d / haaland::coeff_roughness,
                                  haaland::coeff_exponent);
    const double X = rough + haaland::coeff_reynolds / re;
    const double A = haaland::coeff_outer * std::log10(X);
    const double lam_turb = 1.0 / (A * A);
    const double dX_dre = -haaland::coeff_reynolds / (re * re);
    const double dA_dre = haaland::coeff_outer * dX_dre / (X * std::log(10.0));
    const double dlam_turb = -2.0 * dA_dre / (A * A * A);

    constexpr double re_lam = 2300.0;
    constexpr double re_turb = 4000.0;
    if (re <= re_lam) {
        return {lam_lam, dlam_lam};
    }
    if (re >= re_turb) {
        return {lam_turb, dlam_turb};
    }

    const double u = (std::log(re) - std::log(re_lam))
                   / (std::log(re_turb) - std::log(re_lam));
    const double w = u * u * (3.0 - 2.0 * u);
    const double dw_du = 6.0 * u * (1.0 - u);
    const double du_dre = 1.0 / (re * (std::log(re_turb) - std::log(re_lam)));

    const double value = (1.0 - w) * lam_lam + w * lam_turb;
    const double deriv = (1.0 - w) * dlam_lam + w * dlam_turb
                       + dw_du * du_dre * (lam_turb - lam_lam);
    return {value, deriv};
}


// Idelchik (1966) wall-orifice correlations, one class per edge type.
//
// zeta is referenced to the hole velocity and carries the full permanent
// loss, so Cd = 1/sqrt(zeta) and dCd/dRe = -0.5 zeta^-1.5 dzeta/dRe.
//
// NONE of these has a crossflow term: Idelchik's geometry is plenum to
// plenum, with no approach velocity to form U1/Vi from. dCd/d(U1_over_Vi) is
// therefore returned as EXACTLY zero -- an honest absence, not a modelling
// simplification. A hole that does see inlet crossflow wants
// McGreehanSchotsch1988, whose Eq. (17) is the term Idelchik lacks.
class IdelchikWallCorrelation : public DischargeCorrelationBase {
public:
    ~IdelchikWallCorrelation() override = default;

    double Cd(const DischargeHoleGeometry& hole,
              const DischargeHoleState& flow) const override {
        return std::get<0>(Cd_and_derivatives(hole, flow));
    }

    std::tuple<double, double, double> Cd_and_derivatives(
        const DischargeHoleGeometry& hole,
        const DischargeHoleState& flow) const override {
        if (hole.is_inclined()) {
            accept_inclination(hole);
        }
        const Interpolated z = zeta(hole, flow);
        if (!combaero::validation::is_finite_positive(z.value)) {
            throw std::invalid_argument(
                "Idelchik wall orifice: non-positive resistance coefficient");
        }
        const double cd = 1.0 / std::sqrt(z.value);
        const double dcd_dre = -0.5 * z.derivative / (z.value * std::sqrt(z.value));
        return {cd, dcd_dre, 0.0};
    }

protected:
    // Whether this edge type composes with the opt-in inclined inlet. Only
    // the thick-walled hole does (see orifice::idelchik::inclination_*),
    // and it applies the term inside its own zeta(); the sharp, beveled and
    // rounded edges refuse: a thin plate has no conduit entrance, and a
    // beveled or rounded inlet already IS the zeta' term diagram 3-2 would
    // replace.
    virtual void accept_inclination(const DischargeHoleGeometry& hole) const {
        refuse_inclination(hole, name());
    }

    // zeta and dzeta/dRe for this edge type.
    virtual Interpolated zeta(const DischargeHoleGeometry& hole,
                              const DischargeHoleState& flow) const = 0;

    // Idelchik's tables stop at Re = 25. Below it the VALUE is held there and
    // the DERIVATIVE is continued, the same treatment and for the same reason
    // as mcgreehan_schotsch::re_min: a Newton step that wanders below the
    // floor must still see a Re sensitivity or it stalls. Held-at-edge is a
    // documented extrapolation, not a claim about creeping flow.
    static double clamped_re(double Re) {
        return std::max(Re, orifice::idelchik::re_points[0]);
    }

    // Diagram 4-17's low-Re pair, interpolated in log(Re) because the table
    // is log-spaced over four and a half decades.
    static Interpolated phi0_at(double Re) {
        static const std::vector<double> lx = log_re_points();
        const Interpolated r = pchip(lx.data(), orifice::idelchik::zeta_phi0,
                                     orifice::idelchik::re_n, std::log(Re));
        return {r.value, r.derivative / Re};
    }

    static Interpolated eps_at(double Re) {
        static const std::vector<double> lx = log_re_points();
        const Interpolated r = pchip(lx.data(), orifice::idelchik::eps_re,
                                     orifice::idelchik::re_n, std::log(Re));
        return {r.value, r.derivative / Re};
    }

    static Interpolated friction_term(const DischargeHoleGeometry& hole, double Re) {
        const Interpolated lam =
            darcy_friction(Re, orifice::idelchik::default_roughness_over_d);
        const double ld = hole.L_over_d();
        return {lam.value * ld, lam.derivative * ld};
    }

private:
    static std::vector<double> log_re_points() {
        std::vector<double> lx(orifice::idelchik::re_n);
        for (int i = 0; i < orifice::idelchik::re_n; ++i) {
            lx[i] = std::log(orifice::idelchik::re_points[i]);
        }
        return lx;
    }
};

// Diagram 4-17: sharp-edged hole, l/Dh <= 0.015.
class IdelchikSharpCorrelation : public IdelchikWallCorrelation {
public:
    std::string name() const override {
        return "Idelchik (1966) sharp-edged hole in a wall, diagram 4-17";
    }

protected:
    Interpolated zeta(const DischargeHoleGeometry&,
                      const DischargeHoleState& flow) const override {
        // Diagram 4-17's table runs the WHOLE range and converges to
        // zeta_sharp = 2.85 at its last knot, Re = 1e6; item 1's
        // "Re >= 1e5: zeta = 2.85" is the coarse statement of the same
        // curve, not a separate branch. Switching to the constant at 1e5
        // would put a 4.5% step in Cd there (2.60 -> 2.85) -- the source
        // says 2.60 at 1e5. So the table runs everywhere and PCHIP holds it
        // at 2.85 above the last knot, which is exact.
        const double re = clamped_re(flow.Re);
        const Interpolated p = phi0_at(re);
        const Interpolated e = eps_at(re);
        return {p.value + e.value, p.derivative + e.derivative};
    }
};

// Diagram 4-18a: thick-walled (deep) hole.
class IdelchikThickCorrelation : public IdelchikWallCorrelation {
public:
    std::string name() const override {
        return "Idelchik (1966) thick-walled hole in a wall, diagram 4-18a";
    }

protected:
    void accept_inclination(const DischargeHoleGeometry&) const override {}

    Interpolated zeta(const DischargeHoleGeometry& hole,
                      const DischargeHoleState& flow) const override {
        const double re = clamped_re(flow.Re);
        // zeta'(l/Dh) is geometry, so it contributes no dzeta/dRe.
        Interpolated zp = pchip(orifice::idelchik::thick_l_over_d,
                                orifice::idelchik::thick_zeta,
                                orifice::idelchik::thick_n,
                                hole.L_over_d());
        // OPT-IN inclined inlet (orifice::idelchik::inclination_*). The
        // table's zeta' is Eq. (4-20)'s 0.5 + 1 + tau with the sharp NORMAL
        // inlet 0.5 of Eq. (4-27); an inclined inlet takes diagram 3-2's
        // value instead, so zeta' gains the excess. It enters zeta' -- not
        // the total -- because item 2 below scales the high-Re zeta' "as
        // for Re >= 1e5" by k eps_re(Re); at the top of the table (k eps_re
        // = 1) the two readings coincide. Geometry: the increment itself has
        // no Re dependence; the Re derivative below carries it through
        // k d(eps_re)/dRe exactly like the rest of zeta'. Skipped entirely
        // for a normal hole, so the default stays bit-identical.
        if (hole.is_inclined()) {
            zp.value += discharge_inclination_dzeta(hole.inclination_deg);
        }
        const Interpolated fr = friction_term(hole, re);

        // One expression over the whole range, for the same reason as the
        // sharp branch. With k = 1/zeta_sharp this reduces EXACTLY to item
        // 1's zeta' + lam l/Dh at the table's last knot, so there is no
        // branch to step across: eps_re(1e6) = 2.85 and zeta_phi0(1e6) = 0
        // give k*2.85*zeta' = zeta'. Above the last knot PCHIP holds both
        // table values, so the high-Re form continues exactly.
        const Interpolated p = phi0_at(re);
        const Interpolated e = eps_at(re);
        const double k = orifice::idelchik::thick_low_re_coef;
        return {p.value + k * e.value * zp.value + fr.value,
                p.derivative + k * e.derivative * zp.value + fr.derivative};
    }
};

// Diagram 4-18b: beveled edges.
class IdelchikBeveledCorrelation : public IdelchikWallCorrelation {
public:
    std::string name() const override {
        return "Idelchik (1966) beveled-edge hole in a wall, diagram 4-18b";
    }

protected:
    Interpolated zeta(const DischargeHoleGeometry& hole,
                      const DischargeHoleState& flow) const override {
        // Diagram 4-18b tabulates zeta directly against the bevel depth and
        // states no Re dependence; the curve is for Re >= 1e5. Below that the
        // source gives no beveled branch, so the value is held rather than
        // borrowed from the sharp one.
        //
        // Only the VALUE of the interpolant: its derivative is dzeta/d(l/Dh),
        // a geometry slope, and zeta()'s derivative is read as dzeta/dRe.
        // Returned as such it gave a Cd constant in Re a dCd/dRe of 0.83,
        // and the self-consistent row (m - m*, #582) a denominator 1 - K
        // dCd/dm of order -1e5: every column but the flow's ~0 (#964).
        (void)flow;
        return {pchip(orifice::idelchik::beveled_l_over_d,
                      orifice::idelchik::beveled_zeta,
                      orifice::idelchik::beveled_n, hole.bevel_over_d()).value,
                0.0};
    }
};

// Diagram 4-18c: rounded edges.
class IdelchikRoundedCorrelation : public IdelchikWallCorrelation {
public:
    std::string name() const override {
        return "Idelchik (1966) rounded-edge hole in a wall, diagram 4-18c";
    }

protected:
    Interpolated zeta(const DischargeHoleGeometry& hole,
                      const DischargeHoleState& flow) const override {
        // No Re dependence: the interpolant's derivative is d/d(r/Dh), a
        // geometry slope, not dzeta/dRe (see the beveled class, #964).
        (void)flow;
        return {pchip(orifice::idelchik::rounded_r_over_d,
                      orifice::idelchik::rounded_zeta,
                      orifice::idelchik::rounded_n, hole.r_over_d()).value,
                0.0};
    }
};

// Lichtarowicz, Duggins and Markland (1965), Eqs. (7) and (12).
//
// The long-orifice, LOW-Reynolds member of the discharge family. Unlike the
// Idelchik classes this is not a table lookup: the source gives closed-form
// expressions, so Cd and dCd/dRe are both analytic with no interpolation.
//
// There is no crossflow term -- the experiments are plenum-fed -- so
// dCd/d(U1_over_Vi) is EXACTLY zero, as for Idelchik.
class LichtarowiczCorrelation : public DischargeCorrelationBase {
public:
    std::string name() const override {
        return "Lichtarowicz, Duggins and Markland (1965) long orifice";
    }

    double Cd(const DischargeHoleGeometry& hole,
              const DischargeHoleState& flow) const override {
        return std::get<0>(Cd_and_derivatives(hole, flow));
    }

    std::tuple<double, double, double> Cd_and_derivatives(
        const DischargeHoleGeometry& hole,
        const DischargeHoleState& flow) const override {
        namespace L = orifice::lichtarowicz;
        refuse_inclination(hole, name());

        const double ld = hole.L_over_d();
        if (ld < L::l_over_d_min) {
            // The source's own design recommendation (1), not our caution:
            // below l/d = 1.5 "the discharge coefficient varies rapidly with
            // l/d ... and there is the possibility of hysteresis". A
            // correlation cannot represent hysteresis, so refusing is the
            // honest answer rather than returning a single-valued Cd for a
            // geometry the source says does not have one.
            throw std::invalid_argument(
                "Lichtarowicz (1965) is not valid below l/d = 1.5, where the "
                "source reports rapid variation and possible hysteresis. Use "
                "Idelchik1966Thick for a short hole, or supply a measured Cd.");
        }

        const double re = std::max(flow.Re, L::re_floor);

        // Eq. (7), with the source's separate flat value for 1.5 <= l/d < 2.
        // l/d beyond 10 is HELD at 10 rather than extrapolated: Eq. (7) is
        // linear and would keep falling without bound.
        const double ld_eff = std::min(ld, L::l_over_d_max);
        const double cdu = (ld_eff < L::l_over_d_split)
                               ? L::cdu_short
                               : L::cdu_c0 - L::cdu_c1 * ld_eff;

        // Eq. (12)
        const double b = L::visc_c0 * (1.0 + L::visc_c1 * ld_eff);
        const double cc = L::trans_c0 * ld_eff;
        const double u = std::log10(L::trans_c2 * re);
        const double den = 1.0 + L::trans_c1 * u * u;

        const double inv_cd = 1.0 / cdu + b / re - cc / den;
        if (!combaero::validation::is_finite_positive(inv_cd)) {
            throw std::invalid_argument(
                "Lichtarowicz (1965): non-positive 1/Cd; check l/d and Re");
        }
        const double cd = 1.0 / inv_cd;

        // d(1/Cd)/dRe = -b/Re^2 + 15 cc u / (den^2 Re ln10), then
        // dCd/dRe = -Cd^2 d(1/Cd)/dRe. Analytic, per the (f, J) rule.
        const double d_inv_dre = -b / (re * re)
                               + (2.0 * L::trans_c1) * cc * u
                                     / (den * den * re * std::log(10.0));
        const double dcd_dre = -cd * cd * d_inv_dre;

        // Below the floor the VALUE is held and the DERIVATIVE reported as
        // zero, because the value genuinely stops changing there. This
        // differs from the Idelchik floor, where the derivative is continued
        // -- there the floor is the edge of a table and the curve is still
        // moving; here it is a guard against 20/Re diverging.
        return {cd, (flow.Re > L::re_floor) ? dcd_dre : 0.0, 0.0};
    }
};

class McGreehanSchotsch1988Correlation : public DischargeCorrelationBase {
public:
    double Cd(const DischargeHoleGeometry& hole,
              const DischargeHoleState& flow) const override {
        refuse_inclination(hole, name());
        warn_below_floor(flow.Re);
        return orifice::mcgreehan_schotsch::cd(flow.Re, hole.r_over_d(),
                                               hole.L_over_d(), flow.U1_over_Vi);
    }

    std::tuple<double, double, double> Cd_and_derivatives(
        const DischargeHoleGeometry& hole,
        const DischargeHoleState& flow) const override {
        refuse_inclination(hole, name());
        warn_below_floor(flow.Re);
        return orifice::mcgreehan_schotsch::cd_and_derivatives(
            flow.Re, hole.r_over_d(), hole.L_over_d(), flow.U1_over_Vi);
    }

    std::string name() const override { return "McGreehan-Schotsch (1988)"; }

private:
    // Below re_min the chain holds Re AT the floor, so Cd stops responding
    // to flow entirely. That is documented behaviour and deliberate -- Eq.
    // (8) diverges below it -- but silently returning a frozen Cd to a
    // caller who does not know is how a plenum-fed hole came to be read
    // +21.7% high at Re = 432 against Lichtarowicz, which IS valid there.
    // Warn once per call rather than refuse: refusing would break existing
    // networks that transit low Re during Newton iteration.
    static void warn_below_floor(double Re) {
        if (Re < orifice::mcgreehan_schotsch::re_min) {
            combaero::warn(
                "McGreehanSchotsch1988: Re = " + std::to_string(Re) +
                " is below the correlation's floor of " +
                std::to_string(orifice::mcgreehan_schotsch::re_min) +
                "; Cd is held at the floor value and no longer responds to "
                "flow. For a long hole at low Re use Lichtarowicz1965 "
                "(valid 10 to 2e4), or Idelchik1966Thick (valid from 25).");
        }
    }
};

class ConstantDischargeCdCorrelation : public DischargeCorrelationBase {
    double Cd_value_;
public:
    explicit ConstantDischargeCdCorrelation(
        double Cd = orifice::defaults::discharge_cd)
        : Cd_value_(Cd) {}

    double Cd(const DischargeHoleGeometry& hole,
              const DischargeHoleState&) const override {
        refuse_inclination(hole, name());
        return Cd_value_;
    }

    // Exactly zero, not a small number: a constant has no Re or crossflow
    // sensitivity, and a solver is entitled to that being exact.
    std::tuple<double, double, double> Cd_and_derivatives(
        const DischargeHoleGeometry& hole,
        const DischargeHoleState&) const override {
        refuse_inclination(hole, name());
        return {Cd_value_, 0.0, 0.0};
    }

    std::string name() const override { return "Constant Cd"; }
};

} // anonymous namespace

std::unique_ptr<OrificeCorrelationBase> make_correlation(MeteringCdCorrelation id) {
    switch (id) {
        case MeteringCdCorrelation::ReaderHarrisGallagher:
            return std::make_unique<ReaderHarrisGallagherCorrelation>();
        case MeteringCdCorrelation::Constant:
            return std::make_unique<ConstantCdCorrelation>();
        case MeteringCdCorrelation::UserFunction:
            return nullptr;  // Use make_user_correlation instead
    }
    return nullptr;
}

std::unique_ptr<OrificeCorrelationBase> make_constant_correlation(double Cd) {
    return std::make_unique<ConstantCdCorrelation>(Cd);
}

std::unique_ptr<OrificeCorrelationBase> make_user_correlation(
    CdFunction fn,
    const std::string& name) {
    return std::make_unique<UserFunctionCorrelation>(std::move(fn), name);
}

std::unique_ptr<OrificeCorrelationBase> make_tabulated_correlation(
    const std::vector<double>& beta_values,
    const std::vector<double>& Re_values,
    const std::vector<std::vector<double>>& Cd_table,
    const std::string& name) {
    return std::make_unique<TabulatedCorrelation>(beta_values, Re_values, Cd_table, name);
}

double DischargeHoleGeometry::L_over_d() const {
    return (d > 0.0) ? L / d : 0.0;
}

double DischargeHoleGeometry::r_over_d() const {
    return (d > 0.0) ? r / d : 0.0;
}

double DischargeHoleGeometry::bevel_over_d() const {
    return (d > 0.0) ? bevel / d : 0.0;
}

double DischargeHoleGeometry::area() const {
    return M_PI * d * d / 4.0;
}

bool DischargeHoleGeometry::is_inclined() const {
    return inclination_deg < 90.0;
}

bool DischargeHoleGeometry::is_valid() const {
    // L = 0 is a legitimate limit (a knife-edge hole), r = 0 is a sharp
    // inlet. Only a non-positive diameter makes the hole meaningless; an
    // inclination is Idelchik's delta, 0 < delta <= 90.
    return d > 0.0 && L >= 0.0 && r >= 0.0 && inclination_deg > 0.0 &&
           inclination_deg <= 90.0;
}

double discharge_inclination_dzeta(double inclination_deg) {
    if (!(inclination_deg > 0.0 && inclination_deg <= 90.0)) {
        throw std::invalid_argument(
            "discharge_inclination_dzeta: the hole-to-wall angle must be in (0, 90] deg");
    }
    if (inclination_deg == 90.0) {
        return 0.0;  // exactly: cos(pi/2) is 6e-17 in floating point
    }
    namespace I = orifice::idelchik;
    const double c = std::cos(inclination_deg * M_PI / 180.0);
    return I::inclination_cos_coef * c + I::inclination_cos2_coef * c * c;
}

std::unique_ptr<DischargeCorrelationBase> make_discharge_correlation(
    DischargeCdCorrelation id) {
    switch (id) {
        case DischargeCdCorrelation::McGreehanSchotsch1988:
            return std::make_unique<McGreehanSchotsch1988Correlation>();
        case DischargeCdCorrelation::Idelchik1966Sharp:
            return std::make_unique<IdelchikSharpCorrelation>();
        case DischargeCdCorrelation::Idelchik1966Thick:
            return std::make_unique<IdelchikThickCorrelation>();
        case DischargeCdCorrelation::Idelchik1966Beveled:
            return std::make_unique<IdelchikBeveledCorrelation>();
        case DischargeCdCorrelation::Idelchik1966Rounded:
            return std::make_unique<IdelchikRoundedCorrelation>();
        case DischargeCdCorrelation::Lichtarowicz1965:
            return std::make_unique<LichtarowiczCorrelation>();
        case DischargeCdCorrelation::Constant:
            return std::make_unique<ConstantDischargeCdCorrelation>();
    }
    return nullptr;
}

std::unique_ptr<DischargeCorrelationBase> make_constant_discharge_correlation(
    double Cd) {
    return std::make_unique<ConstantDischargeCdCorrelation>(Cd);
}

bool discharge_cd_in_range(DischargeCdCorrelation id,
                           const DischargeHoleGeometry& hole,
                           const DischargeHoleState& flow) {
    namespace I = orifice::idelchik;
    namespace L = orifice::lichtarowicz;
    const double re = flow.Re;
    const double ld = hole.L_over_d();
    const bool idelchik_re = re >= I::re_points[0] && re <= I::re_points[I::re_n - 1];
    if (hole.is_inclined()) {
        // Only the thick-walled hole composes with diagram 3-2; inside its
        // own range it also needs 3-2's Re and a hole long enough for the
        // inlet jet to reattach (orifice::idelchik::inclination_*).
        return id == DischargeCdCorrelation::Idelchik1966Thick && idelchik_re &&
               ld <= I::thick_l_over_d[I::thick_n - 1] && re > I::inclination_re_min &&
               ld >= I::inclination_l_over_d_min;
    }
    switch (id) {
        case DischargeCdCorrelation::McGreehanSchotsch1988:
            return re >= orifice::mcgreehan_schotsch::re_min;
        case DischargeCdCorrelation::Idelchik1966Sharp:
            return idelchik_re && ld <= I::sharp_l_over_d_max;
        case DischargeCdCorrelation::Idelchik1966Thick:
            return idelchik_re && ld <= I::thick_l_over_d[I::thick_n - 1];
        case DischargeCdCorrelation::Idelchik1966Beveled:
            return idelchik_re &&
                   hole.bevel_over_d() <= I::beveled_l_over_d[I::beveled_n - 1];
        case DischargeCdCorrelation::Idelchik1966Rounded:
            return idelchik_re &&
                   hole.r_over_d() <= I::rounded_r_over_d[I::rounded_n - 1];
        case DischargeCdCorrelation::Lichtarowicz1965:
            return re >= L::re_validated_min && re <= L::re_validated_max &&
                   ld >= L::l_over_d_split && ld <= L::l_over_d_max;
        case DischargeCdCorrelation::Constant:
            return true;
    }
    return false;
}

double discharge_cd_re_floor(DischargeCdCorrelation id) {
    switch (id) {
        case DischargeCdCorrelation::McGreehanSchotsch1988:
            return orifice::mcgreehan_schotsch::re_min;
        case DischargeCdCorrelation::Idelchik1966Sharp:
        case DischargeCdCorrelation::Idelchik1966Thick:
        case DischargeCdCorrelation::Idelchik1966Beveled:
        case DischargeCdCorrelation::Idelchik1966Rounded:
            return orifice::idelchik::re_points[0];
        case DischargeCdCorrelation::Lichtarowicz1965:
            return orifice::lichtarowicz::re_floor;
        case DischargeCdCorrelation::Constant:
            return 0.0;
    }
    return 0.0;
}

// -------------------------------------------------------------
// Discharge hole as an orifice in a straight pipe (#645)
// -------------------------------------------------------------

namespace {

// zeta1 n^2: the in-pipe loss on the HOLE velocity head (w0 = w1 / n).
// zeta' and tau come from the wall tables (see orifice::idelchik), so n = 0
// returns the wall hole's high-Re zeta exactly.
double pipe_zeta0(DischargeCdCorrelation id, const DischargeHoleGeometry& hole, double n) {
    namespace I = orifice::idelchik;
    if (hole.is_inclined() && id != DischargeCdCorrelation::Idelchik1966Thick) {
        refuse_inclination(hole, "discharge_pipe_zeta");
    }
    const double m = 1.0 - n;
    const auto edge = [m](double zeta_wall) {
        const double s = std::sqrt(zeta_wall) - 1.0;
        const double a = m + s * std::sqrt(m);
        return a * a;
    };
    switch (id) {
        case DischargeCdCorrelation::Idelchik1966Sharp:
            return edge(I::zeta_sharp);
        case DischargeCdCorrelation::Idelchik1966Beveled:
            return edge(pchip(I::beveled_l_over_d, I::beveled_zeta, I::beveled_n,
                              hole.bevel_over_d())
                            .value);
        case DischargeCdCorrelation::Idelchik1966Rounded:
            return edge(pchip(I::rounded_r_over_d, I::rounded_zeta, I::rounded_n,
                              hole.r_over_d())
                            .value);
        case DischargeCdCorrelation::Idelchik1966Thick:
        case DischargeCdCorrelation::Lichtarowicz1965:
        case DischargeCdCorrelation::McGreehanSchotsch1988: {
            const double ld = hole.L_over_d();
            const double tau =
                pchip(I::thick_l_over_d, I::thick_zeta, I::thick_n, ld).value
                - I::pipe_thick_wall_offset;
            // Eq. (4-20)'s inlet zeta' is diagram 3-2's for an inclined
            // Thick hole, exactly as in the wall correlation; tau stays the
            // normal hole's.
            const double inlet =
                hole.is_inclined()
                    ? I::pipe_thick_inlet + discharge_inclination_dzeta(hole.inclination_deg)
                    : I::pipe_thick_inlet;
            return (inlet + tau * std::sqrt(m)) * m + m * m + I::pipe_thick_lambda * ld;
        }
        case DischargeCdCorrelation::Constant:
            break;
    }
    throw std::invalid_argument("discharge_pipe_zeta: a Constant Cd has no hole geometry");
}

// C solving k(C) (1 - beta^4) / C^2 = z0. The left side falls strictly in C
// (k and 1/C^2 both do), so a bracketed Newton converges from any start.
double pipe_equivalent_c(double beta, double z0) {
    const double b4 = std::pow(beta, 4.0);
    const auto g = [beta, b4, z0](double c, double& dg) {
        const auto [k, dk] = permanent_loss_ratio_and_dC(beta, c);
        dg = (1.0 - b4) * (dk / (c * c) - 2.0 * k / (c * c * c));
        return k * (1.0 - b4) / (c * c) - z0;
    };
    double lo = 1e-6;
    double hi = 1e3;
    double c = 1.0 / std::sqrt(z0);
    c = std::min(std::max(c, lo), hi);
    for (int i = 0; i < combaero::inner_solve::kPipeEquivalentCMaxIter; ++i) {
        double dg = 0.0;
        const double f = g(c, dg);
        if (f > 0.0) {
            lo = c;
        } else {
            hi = c;
        }
        double next = (dg < 0.0) ? c - f / dg : 0.5 * (lo + hi);
        if (!(next > lo && next < hi)) {
            next = 0.5 * (lo + hi);
        }
        if (std::abs(next - c) <= 1e-14 * c) {
            return next;
        }
        c = next;
    }
    // A cap exit is a root only if the bracket has closed onto it (#868).
    if (combaero::inner_solve::cap_exit_is_root(hi - lo, c)) {
        return c;
    }
    combaero::inner_solve::throw_no_convergence(
        "pipe_equivalent_c", combaero::inner_solve::kPipeEquivalentCMaxIter, "bracket width",
        hi - lo);
}

} // namespace

double discharge_pipe_zeta(DischargeCdCorrelation id, const DischargeHoleGeometry& hole,
                           double beta) {
    if (!(beta > 0.0 && beta < 1.0)) {
        throw std::invalid_argument("discharge_pipe_zeta: needs 0 < beta < 1");
    }
    const double n = beta * beta;
    return pipe_zeta0(id, hole, n) / (n * n);
}

double discharge_pipe_cd_ratio(DischargeCdCorrelation id, const DischargeHoleGeometry& hole,
                               double beta) {
    if (!(beta > 0.0 && beta < 1.0) || id == DischargeCdCorrelation::Constant) {
        return 1.0;
    }
    const double z0 = pipe_zeta0(id, hole, beta * beta);
    // beta -> 0: k -> 1 and the relation is C^2 = 1 / z0.
    const double c0 = 1.0 / std::sqrt(pipe_zeta0(id, hole, 0.0));
    return pipe_equivalent_c(beta, z0) / c0;
}

// -------------------------------------------------------------
// Flow calculations
// -------------------------------------------------------------

double orifice_mdot(const OrificeGeometry& geom, double Cd, double dP,
                    double rho, double epsilon) {
    if (!combaero::validation::is_finite_nonnegative(dP) || !combaero::validation::is_finite_positive(rho) || !combaero::validation::is_finite_positive(Cd)) {
        throw std::invalid_argument("orifice_mdot: invalid parameters");
    }
    const double beta = geom.beta();
    const double E = velocity_of_approach_factor(beta);
    return Cd * E * epsilon * geom.area() * std::sqrt(2.0 * rho * dP);
}

double orifice_dP(const OrificeGeometry& geom, double Cd, double mdot,
                  double rho, double epsilon) {
    if (!combaero::validation::is_finite_nonnegative(mdot) || !combaero::validation::is_finite_positive(rho) || !combaero::validation::is_finite_positive(Cd)) {
        throw std::invalid_argument("orifice_dP: invalid parameters");
    }
    // From mdot = Cd * E * epsilon * A * sqrt(2 * rho * dP)
    // Solve for dP: dP = (mdot / (Cd * E * epsilon * A))^2 / (2 * rho)
    const double A = geom.area();
    const double beta = geom.beta();
    const double E = velocity_of_approach_factor(beta);
    const double term = mdot / (Cd * E * epsilon * A);
    return term * term / (2.0 * rho);
}

double orifice_Cd_from_measurement(const OrificeGeometry& geom,
                                    double mdot, double dP, double rho) {
    if (!combaero::validation::is_finite_positive(dP) || !combaero::validation::is_finite_positive(rho) || !combaero::validation::is_finite_positive(mdot)) {
        throw std::invalid_argument("orifice_Cd_from_measurement: invalid parameters");
    }
    const double beta = geom.beta();
    const double E = velocity_of_approach_factor(beta);
    double A = geom.area();
    return mdot / (E * A * std::sqrt(2.0 * rho * dP));
}

// -------------------------------------------------------------
// Iterative solver for Cd-Re coupling
// -------------------------------------------------------------

double solve_orifice_mdot(
    const OrificeGeometry& geom,
    double dP,
    double rho,
    double mu,
    double P_upstream,
    double kappa,
    MeteringCdCorrelation correlation,
    double tol,
    int max_iter)
{
    // Input validation
    if (!geom.is_valid()) {
        throw std::invalid_argument("solve_orifice_mdot: invalid geometry");
    }
    if (!combaero::validation::is_finite_nonnegative(dP)) {
        throw std::invalid_argument("solve_orifice_mdot: dP must be non-negative");
    }
    if (!combaero::validation::is_finite_positive(rho)) {
        throw std::invalid_argument("solve_orifice_mdot: rho must be positive");
    }
    if (!combaero::validation::is_finite_positive(mu)) {
        throw std::invalid_argument("solve_orifice_mdot: mu must be positive");
    }
    if (!combaero::validation::is_finite_positive(P_upstream)) {
        throw std::invalid_argument("solve_orifice_mdot: P_upstream must be positive");
    }
    if (!combaero::validation::is_finite_positive(tol)) {
        throw std::invalid_argument("solve_orifice_mdot: tol must be positive");
    }
    if (max_iter < 1) {
        throw std::invalid_argument("solve_orifice_mdot: max_iter must be >= 1");
    }

    // Handle zero pressure drop case
    if (dP == 0.0) {
        return 0.0;
    }

    // Precompute constants
    const double area = geom.area();
    const double beta = geom.beta();
    const double D = geom.D;

    // Dispatch through the same factory the polymorphic API uses, rather than
    // a second switch that only knew the three ISO correlations. That switch
    // was why IdelchikThick/IdelchikRounded/Constant were implemented but
    // unreachable from here (and from Python, which only ever saw three enum
    // members). Built once: the correlation object is stateless in Re.
    std::unique_ptr<OrificeCorrelationBase> corr = make_correlation(correlation);
    if (!corr) {
        // make_correlation returns nullptr only for UserFunction, which cannot
        // be selected by enum -- there is no function to carry with it.
        throw std::invalid_argument(
            "solve_orifice_mdot: UserFunction cannot be selected by enum; "
            "pass the Cd explicitly or use make_user_correlation");
    }

    // Initial guess for Cd (typical value for sharp orifices)
    double Cd = orifice::defaults::metering_cd;

    // Initial guess for mdot (use incompressible formula with initial Cd)
    double mdot = Cd * area * std::sqrt(2.0 * rho * dP);

    // Iteration loop
    for (int iter = 0; iter < max_iter; ++iter) {
        // Calculate expansibility factor if compressible (kappa > 1)
        double epsilon = 1.0;
        if (kappa > 1.0) {
            epsilon = expansibility_factor(beta, dP, P_upstream, kappa);
        }

        // Calculate velocity-of-approach factor (E)
        const double E = velocity_of_approach_factor(beta);

        // Calculate new mass flow rate using current Cd, E, and epsilon
        const double mdot_new = Cd * E * epsilon * area * std::sqrt(2.0 * rho * dP);

        // Check convergence
        const double rel_error = std::abs(mdot_new - mdot) / (mdot + 1e-30);
        if (rel_error < tol) {
            return mdot_new;
        }

        // Update Reynolds number based on new mdot
        // Re_D = (4 * mdot) / (pi * D * mu)
        const double Re_D = (4.0 * mdot_new) / (M_PI * D * mu);

        // Update Cd based on new Reynolds number
        OrificeState iter_state;
        iter_state.Re_D = Re_D;
        iter_state.dP = dP;
        iter_state.rho = rho;
        iter_state.mu = mu;
        Cd = corr->Cd(geom, iter_state);

        // Update mdot for next iteration
        mdot = mdot_new;
    }

    // Failed to converge
    throw std::runtime_error(
        "solve_orifice_mdot: failed to converge after " +
        std::to_string(max_iter) + " iterations");
}

// -------------------------------------------------------------
// Compressible flow correction
// -------------------------------------------------------------

double expansibility_factor(double beta, double dP, double P_upstream, double kappa) {
    // Input validation
    if (!combaero::validation::is_finite_positive(beta) || beta >= 1.0) {
        throw std::invalid_argument("expansibility_factor: beta must be in range (0, 1)");
    }
    if (!combaero::validation::is_finite_positive(P_upstream)) {
        throw std::invalid_argument("expansibility_factor: P_upstream must be positive");
    }
    if (!combaero::validation::is_finite_nonnegative(dP)) {
        throw std::invalid_argument("expansibility_factor: dP must be non-negative");
    }

    if (!combaero::validation::is_finite(kappa)) {
        throw std::invalid_argument("expansibility_factor: kappa must be finite");
    }

    // Incompressible limit: kappa <= 1 or no pressure drop
    if (kappa <= 1.0 || dP <= 0.0) {
        return 1.0;
    }

    // Pressure ratio
    const double tau = dP / P_upstream;

    // Check validity range (ISO 5167-2 recommends tau <= 0.25)
    if (tau > 0.25) {
        // Still compute but this is outside recommended range
        // User should be aware via documentation
    }

    // Compute beta powers
    const double beta2 = beta * beta;
    const double beta4 = beta2 * beta2;
    const double beta8 = beta4 * beta4;

    // ISO 5167-2:2003 expansibility factor formula
    // eps = 1 - (0.351 + 0.256*beta^4 + 0.93*beta^8) * [1 - (1 - tau)^(1/kappa)]
    namespace ie = orifice::iso5167_expansibility;
    const double coeff = ie::c0 + ie::c4 * beta4 + ie::c8 * beta8;
    const double expansion_term = 1.0 - std::pow(1.0 - tau, 1.0 / kappa);

    return 1.0 - coeff * expansion_term;
}

// -------------------------------------------------------------
// Orifice flow result bundle
// -------------------------------------------------------------

OrificeFlowResult orifice_flow(
    const OrificeGeometry& geom,
    double dP,
    double T,
    double P,
    double mu,
    double Z,
    const std::vector<double>& X,
    double kappa,
    MeteringCdCorrelation correlation)
{
    // Input validation
    if (!geom.is_valid()) {
        throw std::invalid_argument("orifice_flow: invalid geometry");
    }
    if (!combaero::validation::is_finite_positive(T)) {
        throw std::invalid_argument("orifice_flow: temperature must be positive");
    }
    if (!combaero::validation::is_finite_positive(P)) {
        throw std::invalid_argument("orifice_flow: pressure must be positive");
    }
    if (!combaero::validation::is_finite_positive(mu)) {
        throw std::invalid_argument("orifice_flow: viscosity must be positive");
    }
    if (!combaero::validation::is_finite_positive(Z)) {
        throw std::invalid_argument("orifice_flow: compressibility factor Z must be positive");
    }
    if (!combaero::validation::is_finite_nonnegative(dP)) {
        throw std::invalid_argument("orifice_flow: differential pressure must be non-negative");
    }

    // Ideal-gas density of the given composition (dry air when X is empty,
    // as documented). It used air's molar mass for every X (#829).
    const double rho_ideal = combaero::density(T, P, X.empty() ? combaero::dry_air() : X);

    // Apply real gas correction
    const double rho_corrected = rho_ideal / Z;

    // Solve for mass flow rate with corrected density
    const double mdot = solve_orifice_mdot(
        geom, dP, rho_corrected, mu, P, kappa, correlation);

    // Calculate expansibility factor
    const double beta = geom.beta();
    const double epsilon = (kappa > 1.0) ?
        expansibility_factor(beta, dP, P, kappa) : 1.0;

    // Calculate velocity through orifice
    const double A = geom.area();
    const double v = mdot / (rho_corrected * A);

    // Calculate Reynolds numbers
    const double D = geom.D;
    const double d = geom.d;
    const double Re_D = (4.0 * mdot) / (M_PI * D * mu);
    const double Re_d = (4.0 * mdot) / (M_PI * d * mu);

    // Get discharge coefficient
    OrificeState state;
    state.Re_D = Re_D;
    state.dP = dP;
    state.rho = rho_corrected;
    state.mu = mu;

    // Third dispatch site, now the same factory as the other two. The
    // `default: auto-select by geometry` arm this replaces is deliberately
    // gone: picking a correlation from r and t behind the caller's back is
    // how a rounded-entry request came back as the Stolz Cd (removed in 0.9.0).
    std::unique_ptr<OrificeCorrelationBase> corr = make_correlation(correlation);
    if (!corr) {
        throw std::invalid_argument(
            "orifice_flow: UserFunction cannot be selected by enum; "
            "use make_user_correlation");
    }
    const double Cd_value = corr->Cd(geom, state);

    // Populate result struct
    OrificeFlowResult result;
    result.mdot = mdot;
    result.v = v;
    result.Re_D = Re_D;
    result.Re_d = Re_d;
    result.Cd = Cd_value;
    result.epsilon = epsilon;
    result.rho_corrected = rho_corrected;

    return result;
}

// -------------------------------------------------------------
// Utility functions
// -------------------------------------------------------------

double orifice_velocity_from_mdot(double mdot, double rho, double d, double Z) {
    if (!combaero::validation::is_finite_nonnegative(mdot)) {
        throw std::invalid_argument("orifice_velocity_from_mdot: mdot must be non-negative");
    }
    if (!combaero::validation::is_finite_positive(rho)) {
        throw std::invalid_argument("orifice_velocity_from_mdot: rho must be positive");
    }
    if (!combaero::validation::is_finite_positive(d)) {
        throw std::invalid_argument("orifice_velocity_from_mdot: d must be positive");
    }
    if (!combaero::validation::is_finite_positive(Z)) {
        throw std::invalid_argument("orifice_velocity_from_mdot: Z must be positive");
    }

    // Apply real gas correction to density
    const double rho_corrected = rho / Z;

    // Calculate area
    const double A = M_PI * d * d / 4.0;

    // v = mdot / (rho_corrected * A)
    return mdot / (rho_corrected * A);
}

double orifice_area_from_beta(double D, double beta) {
    if (!combaero::validation::is_finite_positive(D)) {
        throw std::invalid_argument("orifice_area_from_beta: D must be positive");
    }
    if (!combaero::validation::is_finite_positive(beta) || beta >= 1.0) {
        throw std::invalid_argument("orifice_area_from_beta: beta must be in range (0, 1)");
    }

    // A = pi * (D * beta / 2)^2
    const double d = D * beta;
    return M_PI * d * d / 4.0;
}

double beta_from_diameters(double d, double D) {
    if (!combaero::validation::is_finite_positive(d)) {
        throw std::invalid_argument("beta_from_diameters: d must be positive");
    }
    if (!combaero::validation::is_finite_positive(D)) {
        throw std::invalid_argument("beta_from_diameters: D must be positive");
    }
    if (d >= D) {
        throw std::invalid_argument("beta_from_diameters: d must be < D");
    }

    // beta = d / D
    return d / D;
}

double orifice_Re_d_from_mdot(double mdot, double d, double mu) {
    if (!combaero::validation::is_finite_nonnegative(mdot)) {
        throw std::invalid_argument("orifice_Re_d_from_mdot: mdot must be non-negative");
    }
    if (!combaero::validation::is_finite_positive(d)) {
        throw std::invalid_argument("orifice_Re_d_from_mdot: d must be positive");
    }
    if (!combaero::validation::is_finite_positive(mu)) {
        throw std::invalid_argument("orifice_Re_d_from_mdot: mu must be positive");
    }

    // Re_d = 4 * mdot / (pi * d * mu)
    return (4.0 * mdot) / (M_PI * d * mu);
}

// -------------------------------------------------------------
// Jobson (1955) compressible contraction (#854)
// -------------------------------------------------------------

namespace {
using J1 = combaero::solver::DualN<1>;

void check_jobson_inputs(const char* who, double Ci, double beta, double n) {
    if (!(Ci > 0.0) || !std::isfinite(Ci)) {
        throw std::invalid_argument(std::string(who) + ": Ci must be > 0");
    }
    if (!(beta >= 0.0 && beta < 1.0)) {
        throw std::invalid_argument(std::string(who) + ": beta must be in [0, 1)");
    }
    const double b4 = beta * beta * beta * beta;
    if (!(Ci * Ci * b4 < 1.0)) {
        throw std::invalid_argument(std::string(who) + ": Ci^2 beta^4 must be < 1");
    }
    if (!(n > 1.0) || !std::isfinite(n)) {
        throw std::invalid_argument(std::string(who) + ": n must be > 1");
    }
}
} // namespace

double jobson_force_defect(double Ci, double beta, double Cv) {
    if (!combaero::validation::is_finite_positive(Ci) ||
        !combaero::validation::is_finite_positive(Cv)) {
        throw std::invalid_argument("jobson_force_defect: Ci and Cv must be > 0");
    }
    if (!combaero::validation::is_finite(beta)) {
        throw std::invalid_argument("jobson_force_defect: beta must be finite");
    }
    const double b2 = (beta > 0.0 && beta < 1.0) ? beta * beta : 0.0;
    return 1.0 / Ci - 0.5 / (Ci * Ci) - Cv * b2 + 0.5 * Cv * Cv * b2 * b2;
}

double jobson_ci_from_cd(double Cd, double beta, double Cv) {
    if (!combaero::validation::is_finite_positive(Cd) ||
        !combaero::validation::is_finite_positive(Cv)) {
        throw std::invalid_argument("jobson_ci_from_cd: Cd and Cv must be > 0");
    }
    if (!combaero::validation::is_finite(beta)) {
        throw std::invalid_argument("jobson_ci_from_cd: beta must be finite");
    }
    return combaero::orifice_jobson::ci_from_cd(J1::constant(Cd), beta, Cv).v;
}

double jobson_contraction_coefficient(double r, double Ci, double n) {
    check_jobson_inputs("jobson_contraction_coefficient", Ci, 0.0, n);
    if (!(r > 0.0 && r <= 1.0)) {
        throw std::invalid_argument("jobson_contraction_coefficient: r must be in (0, 1]");
    }
    if (r == 1.0) {
        return Ci;  // the incompressible limit, where the quadratic is 0/0
    }
    // Reservoir at p0 = rho0 = 1: C depends only on r, n and Ci.
    const J1 one = J1::constant(1.0);
    return combaero::orifice_jobson::flux(one, one, one, J1::constant(r), J1::constant(Ci),
                                          J1::constant(n), 0.0)
        .C.v;
}

double jobson_expansibility(double r, double beta, double Ci, double n) {
    check_jobson_inputs("jobson_expansibility", Ci, beta, n);
    if (!(r > 0.0 && r <= 1.0)) {
        throw std::invalid_argument("jobson_expansibility: r must be in (0, 1]");
    }
    if (r == 1.0) {
        return 1.0;
    }
    const double b2 = beta * beta;
    // Upstream static state p1 = rho1 = 1; the approach velocity raises the
    // stagnation state the jet expands from, V1 = g beta^2 / rho1, so g solves
    // h(g) = F(g) - g = 0, F(g) = flux(p0(g), rho0(g)). F falls with g (every
    // point of the #957 sweep), so the root is unique and h' = F' - 1 < -1.
    // At high beta Ci, F' < -1 at the root, where the fixed point g = F(g)
    // cycled instead of converging (#957). Bracketed Newton (rtsafe) with F'
    // by forward mode: h(0) = F(0) > 0, and the upper end doubles from F(0)
    // until h < 0.
    using combaero::inner_solve::kJobsonExpansibilityMaxIter;
    using combaero::inner_solve::kJobsonExpansibilityRelTol;
    auto h_of = [&](double g, double& dh) {
        const J1 V1 = J1::seed(g, 0) * b2;
        const J1 tau = 1.0 + V1 * V1 * (n - 1.0) / (2.0 * n);  // T0/T1
        const J1 p0 = combaero::solver::dpow(tau, n / (n - 1.0));
        const J1 rho0 = combaero::solver::dpow(tau, 1.0 / (n - 1.0));
        const J1 F = combaero::orifice_jobson::flux(p0, rho0, J1::constant(1.0),
                                                    J1::constant(r), J1::constant(Ci),
                                                    J1::constant(n), beta)
                         .g;
        dh = F.d[0] - 1.0;
        return F.v - g;
    };
    double dh = 0.0;
    double lo = 0.0;
    double hi = h_of(0.0, dh);  // F(0) > 0 for r < 1
    double h_hi = h_of(hi, dh);
    for (int k = 0; h_hi > 0.0; ++k) {
        if (k == kJobsonExpansibilityMaxIter || !std::isfinite(h_hi)) {
            combaero::inner_solve::throw_no_convergence(
                "jobson_expansibility", k, "upper bracket residual", h_hi);
        }
        lo = hi;
        hi *= 2.0;
        h_hi = h_of(hi, dh);
    }
    double g = hi;
    double h = h_hi;
    double dx_old = hi - lo;
    double dx = dx_old;
    for (int it = 0; it < kJobsonExpansibilityMaxIter; ++it) {
        if (!std::isfinite(h)) {
            combaero::inner_solve::throw_no_convergence("jobson_expansibility", it,
                                                        "residual F(g) - g", h);
        }
        if (h == 0.0) {
            break;
        }
        if (h > 0.0) {
            lo = g;
        } else {
            hi = g;
        }
        const bool newton_ok = std::isfinite(dh) && dh < 0.0;
        const double step = newton_ok ? -h / dh : 0.0;
        double g_new = 0.0;
        if (newton_ok && g + step > lo && g + step < hi &&
            std::abs(2.0 * step) <= std::abs(dx_old)) {
            dx_old = dx;
            dx = step;
            g_new = g + step;
        } else {
            dx_old = dx;
            dx = 0.5 * (hi - lo);
            g_new = lo + dx;
        }
        // A step at the tolerance, or a bracket that cannot be split further:
        // the root to rounding.
        const bool done = std::abs(g_new - g) <= kJobsonExpansibilityRelTol * std::abs(g_new) ||
                          !(g_new > lo && g_new < hi);
        g = g_new;
        if (done) {
            break;
        }
        h = h_of(g, dh);
        if (it + 1 == kJobsonExpansibilityMaxIter) {
            if (combaero::inner_solve::cap_exit_is_root(hi - lo, g)) {
                break;
            }
            combaero::inner_solve::throw_no_convergence(
                "jobson_expansibility", kJobsonExpansibilityMaxIter, "residual F(g) - g", h);
        }
    }
    const double b4 = b2 * b2;
    const double g_inc = Ci * std::sqrt(2.0 * (1.0 - r)) / std::sqrt(1.0 - Ci * Ci * b4);
    return g / g_inc;
}
