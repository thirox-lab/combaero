#include "crossflow_cd.h"

#include "inner_solve.h"
#include "input_validation.h"
#include "math_constants.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace combaero {

namespace {

// d(Cd)/d(q) carried onto (G_ideal, G_c) through q = G_ideal / G_c.
CrossflowCd from_q(double Cd, double dCd_dq, double q, double G_c) {
  CrossflowCd out;
  out.Cd = Cd;
  out.dCd_dG_ideal = dCd_dq / G_c;
  out.dCd_dG_c = -dCd_dq * q / G_c;
  return out;
}

CrossflowCd capped(double Cd0) {
  CrossflowCd out;
  out.Cd = Cd0;
  out.capped = true;
  return out;
}

struct Table41Plate {
  double xn_d, yn_d, z_d;
  double cd0;       // Table 3.1, up to the break point
  double x_break;   // the printed break point (documentation; the law joins by min)
  double mult;      // multiplier on zeta
  double lin_a;     // linear piece Cd = lin_a x + lin_b (0: none)
  double lin_b;     // printed intercept (documentation; see below)
  double x_lin_end; // where the linear piece hands over to zeta
};

constexpr Table41Plate TABLE_4_1[] = {
    {5.0, 4.0, 1.0, 0.85, 0.63, 1.00, 0.0, 0.0, 0.0},
    {5.0, 4.0, 2.0, 0.85, 0.83, 1.08, 0.0, 0.0, 0.0},
    {5.0, 4.0, 3.0, 0.85, 0.90, 1.11, 0.0, 0.0, 0.0},
    {5.0, 8.0, 1.0, 0.80, 0.54, 1.00, -0.169, 0.893, 1.5},
    {10.0, 4.0, 1.0, 0.76, 0.54, 1.00, -0.128, 0.825, 1.8},
    {10.0, 8.0, 1.0, 0.76, 0.54, 1.00, -0.128, 0.825, 1.8},
};

const Table41Plate* find_plate(double xn_d, double yn_d, double z_d) {
  for (const auto& p : TABLE_4_1) {
    if (std::abs(p.xn_d - xn_d) < 1e-9 && std::abs(p.yn_d - yn_d) < 1e-9 &&
        std::abs(p.z_d - z_d) < 1e-9) {
      return &p;
    }
  }
  return nullptr;
}

}  // namespace

CrossflowCd rogers_hersh_1976_cd(double Cd0, double G_ideal, double G_c, bool suction) {
  if (!(Cd0 > 0.0)) {
    throw std::invalid_argument("rogers_hersh_1976_cd: Cd0 must be positive");
  }
  // G_c <= 0 means no crossflow and G_ideal <= 0 no jet; NaN meant them too.
  if (!validation::is_finite(G_ideal) || !validation::is_finite(G_c)) {
    throw std::invalid_argument("rogers_hersh_1976_cd: G_ideal and G_c must be finite");
  }
  if (!(G_c > 0.0)) {
    return capped(Cd0);
  }
  if (!(G_ideal > 0.0)) {
    return CrossflowCd{};
  }
  const double a = suction ? ROGERS_HERSH_1976_SUCTION_A : ROGERS_HERSH_1976_INJECTION_A;
  const double b = suction ? ROGERS_HERSH_1976_SUCTION_B : ROGERS_HERSH_1976_INJECTION_B;
  // Cd = a (Cd q)^b  =>  Cd = a^(1/(1-b)) q^(b/(1-b)).
  const double p = 1.0 / (1.0 - b);
  const double e = b * p;
  const double q = G_ideal / G_c;
  const double cd = std::pow(a, p) * std::pow(q, e);
  // C1 join to the cap: u = cd / Cd0 is kept below u0 = 1 - w, held at 1
  // above u1 = 1 + w, and the quadratic u - (u - u0)^2 / (4 w) in between.
  const double w = ROGERS_HERSH_1976_CAP_BAND;
  const double u = cd / Cd0;
  if (u >= 1.0 + w) {
    return capped(Cd0);
  }
  const double dcd_dq = e * cd / q;
  if (u <= 1.0 - w) {
    return from_q(cd, dcd_dq, q, G_c);
  }
  const double s = u - (1.0 - w);
  const double f = u - s * s / (4.0 * w);
  const double df_du = 1.0 - s / (2.0 * w);
  return from_q(Cd0 * f, df_du * dcd_dq, q, G_c);
}

bool florschuetz_metzger_su_1982_has_plate(double xn_d, double yn_d, double z_d) {
  return find_plate(xn_d, yn_d, z_d) != nullptr;
}

CrossflowCd florschuetz_metzger_su_1982_cd(double xn_d, double yn_d, double z_d, double G_ideal,
                                           double G_c) {
  const Table41Plate* plate = find_plate(xn_d, yn_d, z_d);
  if (plate == nullptr) {
    throw std::invalid_argument(
        "florschuetz_metzger_su_1982_cd: CR-3630 Table 4.1 covers the B(5,4,1), "
        "B(5,4,2), B(5,4,3), B(5,8,1), B(10,4,1) and B(10,8,1) plates only");
  }
  if (!validation::is_finite(G_ideal) || !validation::is_finite(G_c)) {
    throw std::invalid_argument(
        "florschuetz_metzger_su_1982_cd: G_ideal and G_c must be finite");
  }
  if (!(G_c > 0.0)) {
    return capped(plate->cd0);
  }
  if (!(G_ideal > 0.0)) {
    return CrossflowCd{};
  }
  const double q = G_ideal / G_c;
  const double m = plate->mult * FLORSCHUETZ_METZGER_SU_1982_ZETA_C;
  const double x0 = FLORSCHUETZ_METZGER_SU_1982_ZETA_X0;
  const double n = FLORSCHUETZ_METZGER_SU_1982_ZETA_N;
  const bool has_line = plate->lin_a != 0.0;
  const double A = plate->lin_a;
  // The line through zeta at its end, with the printed slope (its intercept
  // moves 0.0015 / 0.0006 from the printed 0.893 / 0.825), so the pieces
  // meet. Within the printed digits.
  const double end = plate->x_lin_end;
  const double B = has_line ? m * std::pow(end + x0, -n) - A * end : 0.0;
  // The falling part resolved in q, unbounded by cd0: the line's root,
  // or zeta's (Newton; dF/dCd = 1 - n x/(x + x0) > 0, a unique root).
  double cd_f = 0.0;
  double dcdf_dq = 0.0;
  bool on_line = false;
  if (has_line) {
    // Cd = A / (Cd q) + B  =>  q Cd^2 - B q Cd - A = 0, the larger root.
    const double disc = B * B * q * q + 4.0 * A * q;
    if (disc >= 0.0) {
      const double cd = (B * q + std::sqrt(disc)) / (2.0 * q);
      if (cd > 0.0 && 1.0 / (cd * q) <= end) {
        const double F_cd = 1.0 + A / (cd * cd * q);
        const double F_q = A / (cd * q * q);
        cd_f = cd;
        dcdf_dq = -F_q / F_cd;
        on_line = true;
      }
    }
  }
  if (!on_line) {
    double cd = plate->cd0;
    const double cd_max = 4.0 * plate->cd0;
    for (int it = 0;; ++it) {
      const double x = 1.0 / (cd * q);
      const double z = m * std::pow(x + x0, -n);
      const double dz_dx = -n * z / (x + x0);
      const double F_cd = 1.0 + dz_dx / (cd * cd * q);  // dx/dCd = -x/Cd
      double next = cd - (cd - z) / F_cd;
      if (!(next > 0.0)) {
        next = 0.5 * cd;
      }
      if (next > cd_max) {
        next = cd_max;
      }
      const bool done = std::abs(next - cd) <= 1e-15 * plate->cd0;
      const double step = next - cd;
      cd = next;
      if (done) {
        break;
      }
      if (it + 1 == inner_solve::kCrossflowCdMaxIter) {
        if (inner_solve::cap_exit_is_root(step, plate->cd0)) {
          break;
        }
        inner_solve::throw_no_convergence("crossflow_cd (falling branch)",
                                          inner_solve::kCrossflowCdMaxIter, "Cd step", step);
      }
    }
    const double x = 1.0 / (cd * q);
    const double dz_dx = -n * m * std::pow(x + x0, -n) / (x + x0);
    const double F_cd = 1.0 + dz_dx / (cd * cd * q);
    const double F_q = dz_dx / (cd * q * q);  // dx/dq = -x/q
    cd_f = cd;
    dcdf_dq = -F_q / F_cd;
  }
  // The constant joins C1 as for Rogers & Hersh: the falling part below
  // (1 - w) cd0, cd0 above (1 + w) cd0, the quadratic in between.
  const double w = FLORSCHUETZ_METZGER_SU_1982_CAP_BAND;
  const double u = cd_f / plate->cd0;
  if (u >= 1.0 + w) {
    return capped(plate->cd0);
  }
  if (u <= 1.0 - w) {
    return from_q(cd_f, dcdf_dq, q, G_c);
  }
  const double sb = u - (1.0 - w);
  const double f = u - sb * sb / (4.0 * w);
  return from_q(plate->cd0 * f, (1.0 - sb / (2.0 * w)) * dcdf_dq, q, G_c);
}

// ---- Rowbury, Oldfield & Lock (2001): exit crossflow -----------------------

namespace {

const rowbury2001::Row& rowbury_row(Rowbury2001Set set) {
  switch (set) {
    case Rowbury2001Set::CylindricalPressureSurface:
      return rowbury2001::kCylPressure;
    case Rowbury2001Set::CylindricalLeadingEdge:
      return rowbury2001::kCylLeadingEdge;
    case Rowbury2001Set::CylindricalSuctionSurface:
      return rowbury2001::kCylSuction;
    case Rowbury2001Set::FanShapedPressureSurface:
      return rowbury2001::kFanPressure;
    case Rowbury2001Set::FanShapedSuctionSurface:
      return rowbury2001::kFanSuction;
  }
  throw std::invalid_argument("rowbury_2001: unknown coefficient set");
}

const rowbury2001::Envelope& rowbury_envelope(Rowbury2001Set set) {
  switch (set) {
    case Rowbury2001Set::CylindricalPressureSurface:
      return rowbury2001::kEnvCylPressure;
    case Rowbury2001Set::CylindricalLeadingEdge:
      return rowbury2001::kEnvCylLeadingEdge;
    case Rowbury2001Set::CylindricalSuctionSurface:
      return rowbury2001::kEnvCylSuction;
    case Rowbury2001Set::FanShapedPressureSurface:
      return rowbury2001::kEnvFanPressure;
    case Rowbury2001Set::FanShapedSuctionSurface:
      return rowbury2001::kEnvFanSuction;
  }
  throw std::invalid_argument("rowbury_2001: unknown coefficient set");
}

double sin_deg(double deg) { return std::sin(deg * M_PI / 180.0); }

}  // namespace

Rowbury2001Coefficients rowbury_2001_coefficients(Rowbury2001Set set, double theta_deg,
                                                  double alpha_deg) {
  if (!(theta_deg > 0.0 && theta_deg <= 90.0)) {
    throw std::invalid_argument(
        "rowbury_2001: theta (hole axis to the surface) must be in (0, 90] deg");
  }
  if (!(alpha_deg >= 0.0 && alpha_deg <= 90.0)) {
    throw std::invalid_argument(
        "rowbury_2001: alpha (orientation / compound angle) must be in [0, 90] deg");
  }
  const rowbury2001::Row& r = rowbury_row(set);
  // sin(90 deg) is exactly 1 in floating point; sin(0) exactly 0.
  const double st = sin_deg(theta_deg);
  const double sa = sin_deg(alpha_deg);
  Rowbury2001Coefficients c;
  c.A = (r.a1 * st + r.a0) * (1.0 + r.sa * sa);
  c.B = (r.b1 * st + r.b0) * (1.0 + r.sb * sa);
  c.C = r.c;
  if (!(c.A > 0.0 && c.B > 0.0)) {
    throw std::invalid_argument(
        "rowbury_2001: Table 2 gives A = " + std::to_string(c.A) + ", B = " +
        std::to_string(c.B) + " at theta = " + std::to_string(theta_deg) +
        " deg; Eq. 8 needs A > 0 and B > 0 (the hole is shallower than the "
        "set's fitted angles reach).");
  }
  return c;
}

double rowbury_2001_k_floor(Rowbury2001Set set, double theta_deg, double alpha_deg) {
  const Rowbury2001Coefficients c = rowbury_2001_coefficients(set, theta_deg, alpha_deg);
  double k = ROWBURY_2001_K_PERP_MIN;
  if (c.B > 1.0) {
    // (B - 1) A k^-B = M (1 - C)
    const double k_unique =
        std::pow((c.B - 1.0) * c.A / (ROWBURY_2001_UNIQUE_MARGIN * (1.0 - c.C)), 1.0 / c.B);
    k = std::max(k, k_unique);
  }
  return k;
}

bool rowbury_2001_in_range(Rowbury2001Set set, double theta_deg, double alpha_deg) {
  const rowbury2001::Envelope& e = rowbury_envelope(set);
  constexpr double tol = 1e-9;
  return theta_deg >= e.theta_min - tol && theta_deg <= e.theta_max + tol &&
         alpha_deg >= e.alpha_min - tol && alpha_deg <= e.alpha_max + tol;
}

Rowbury2001Cd rowbury_2001_exit_crossflow_cd(double Cd0, double k_perp, Rowbury2001Set set,
                                             double theta_deg, double alpha_deg) {
  if (!(Cd0 > 0.0 && Cd0 <= 1.0)) {
    throw std::invalid_argument("rowbury_2001_exit_crossflow_cd: Cd0 must be in (0, 1]");
  }
  if (!(k_perp >= 0.0) || std::isnan(k_perp)) {
    throw std::invalid_argument("rowbury_2001_exit_crossflow_cd: k_perp must be >= 0");
  }
  const Rowbury2001Coefficients c = rowbury_2001_coefficients(set, theta_deg, alpha_deg);
  const double k_floor = rowbury_2001_k_floor(set, theta_deg, alpha_deg);
  Rowbury2001Cd out;
  out.held = k_perp < k_floor;
  out.k_used = out.held ? k_floor : k_perp;
  // k = +inf (a vanishing mainstream on a flowing hole) gives A k^-B = 0:
  // the -C asymptote, finite.
  const double ak = std::isinf(out.k_used) ? 0.0 : c.A * std::pow(out.k_used, -c.B);
  out.delta = ak - c.C;
  const double x = 1.0 / (Cd0 * Cd0) + out.delta;
  if (!(x > 0.0)) {
    throw std::invalid_argument(
        "rowbury_2001_exit_crossflow_cd: 1/Cd0^2 + delta_out <= 0, no discharge coefficient");
  }
  out.Cd = 1.0 / std::sqrt(x);
  const double cd3 = out.Cd * out.Cd * out.Cd;
  // x = Cd0^-2 + delta: dCd/dCd0 = (-x^-3/2 / 2)(-2 Cd0^-3) = Cd^3 / Cd0^3.
  out.dCd_dCd0 = cd3 / (Cd0 * Cd0 * Cd0);
  // d delta/dk = -B A k^-B / k; dCd/dk = -Cd^3 / 2 d delta/dk. Zero where held.
  if (!out.held && !std::isinf(out.k_used)) {
    out.dCd_dk = 0.5 * cd3 * c.B * ak / out.k_used;
  }
  return out;
}

Rowbury2001KPerp rowbury_2001_k_perp(double G_h, double rho_h, double rho_m_u_m2,
                                     double theta_deg) {
  if (!(G_h >= 0.0) || !std::isfinite(G_h)) {
    throw std::invalid_argument("rowbury_2001_k_perp: G_h must be finite and >= 0");
  }
  validation::require_finite_positive(rho_h, "rowbury_2001_k_perp: rho_h must be > 0");
  validation::require_finite_positive(
      rho_m_u_m2, "rowbury_2001_k_perp: the mainstream momentum flux rho_m u_m^2 must be > 0");
  if (!(theta_deg > 0.0 && theta_deg <= 90.0)) {
    throw std::invalid_argument("rowbury_2001_k_perp: theta must be in (0, 90] deg");
  }
  const double s = sin_deg(theta_deg);
  Rowbury2001KPerp out;
  out.k = G_h * G_h * s / (rho_h * rho_m_u_m2);
  out.dk_dG = 2.0 * G_h * s / (rho_h * rho_m_u_m2);
  out.dk_drho_h = -out.k / rho_h;
  out.dk_dq_m = -out.k / rho_m_u_m2;
  return out;
}

}  // namespace combaero
