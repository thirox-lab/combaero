#pragma once

// -----------------------------------------------------------------------------
// Ideal momentum junction: whole-element (f, J), parameter-free (#517).
//
// A baseline closure, parallel to the shipped Mynard tee (mpce_junction.h). It
// is not fitted to any data and is not meant to beat tuned losses on their own
// data; it is meant to be physically honest and robust where tuned closures are
// out of range.
//
// PHYSICS. One junction static pressure p_jn. Every port i has a fixed outward
// axis n_i (pointing away from the junction) and a signed mass flow m_i,
// positive OUT of the junction, so its velocity vector is
//
//     v_i = u_i n_i,   u_i = m_i / (rho_i A_i)
//
// for inflow and outflow alike. The supplying streams arrive at p_jn (Benson)
// and mix into their mass-weighted mean velocity V = sum_in |m_k| v_k /
// sum_in |m_k|. Each receiving port j then takes its momentum along its own
// axis from V, the misaligned part being lost to the walls:
//
//     supplier k:   P_k - p_jn = 0
//     receiver j:   P_j - p_jn + (m_j / A_j) (u_j - V . n_j) = 0
//
// i.e. (p_jn - P_j) A_j = m_j (u_j - V . n_j), the momentum balance along n_j.
//
//   - Merge (one receiver on the main axis): (p_jn - P_c) A_c = rho_c u_c^2 A_c
//     - sum_k rho_k u_k^2 A_k cos d_k, Blaisdell & Manson's theory (their Eqs
//     12/19; `bassett_cv(..., wall=False)` in validation/junction).
//   - Split (one supplier): V = v_s, and each receiver loses
//     Pt_s - Pt_j = rho/2 |v_j - v_s|^2 >= 0, the vector Borda-Carnot loss --
//     Idelchik's diverging-branch Eq 7-6 with A' = 1, and kappa = 1 (straight
//     run (1 - q)^2 on the straight fraction) for a collinear receiver.
//
// COMPRESSIBLE. The balances are written in static pressure with each port's
// own density, so they are the exact control-volume momentum balance at any
// Mach number (for uniform port profiles). There is no Mach inversion here:
// Pt <-> P is the port node's own closure.
//
// ROLES ARE SMOOTH. Which ports supply and which receive follows the flow
// signs, so the rows would switch form at m_i = 0 (value continuous, slope
// jump m_s cos / (rho A)). Instead every port carries both roles with smooth
// weights
//
//     s(m) = (1 + m / sqrt(m^2 + eps^2)) / 2,  w_out = m s(m),  w_in = w_out - m
//
// with eps = sqrt((IDEAL_ROLE_EPS_REL * |m|_2)^2 + IDEAL_ROLE_EPS_ABS^2), a
// fixed fraction of the junction's throughput. Both weights vanish at m = 0, so
// a dead port neither dilutes V nor departs from p_jn; away from reversal the
// smoothing moves a row by O((eps/m)^2). The rows are C-infinity everywhere,
// including the all-zero state (V's denominator carries IDEAL_ROLE_EPS_ABS),
// so there is no degenerate-state fallback and no direction classification.
//
// LOSS FACTOR (#526). Each receiver's row rearranges exactly to
//
//     Pt_j = Pt_mix,j - L_j,  Pt_mix,j = p_jn + rho_j |V|^2 / 2,
//                             L_j = rho_j / 2 |v_j - V|^2 >= 0,
//
// so at uniform density the junction's mass-weighted Pt dissipation is
//
//     rho/2 sum_sup |m_k| |v_k - V|^2  +  sum_rec m_j L_j
//
// (mixing of the suppliers, then the receivers' vector Borda loss), both sums
// of squares. An optional per-port factor f_j >= 0 scales L_j only:
// Pt_j = Pt_mix,j - f_j L_j, added to the row as (f_j - 1) s_j L_j with s_j
// the smooth receiver weight, so a supplier's row is untouched and f_j = 1 is
// the ideal tee exactly (the term is not formed). The mixing term stays at
// its ideal value, so f_j = 0 is lossless only in a split. L_j >= 0 at any
// density, so no f_j >= 0 makes the mechanical dissipation negative. The
// factor is a rig-matching knob for the user; the parameter-free closure is
// f = 1. The same algebra is Pt = P + rho u^2/2, i.e. the low-Mach reading of
// the loss; the rows themselves stay the exact momentum balance.
//
// FADE AT A DEAD PORT. With f != 1 a receiver carrying vanishing flow still
// has L_j -> rho_j |V|^2 / 2, so its static tends to p_jn + (1 - f) rho|V|^2/2
// while as a supplier it sits at p_jn: a step of (f - 1) supplied heads across
// reversal (measured 4.0 heads at f = 5 over +-1e-4 of the throughput, none at
// f = 1). The role blend alone turned that into a cliff Newton stalled on, or
// landed in and was refused, whenever a branch shut off or reversed. So the
// factor fades to 1 with the port's flow:
//
//     f_eff,j = 1 + (f_j - 1) phi_j,  phi_j = m_j^2 / (m_j^2 + delta^2),
//
// delta = IDEAL_LOSS_FADE_REL of the throughput. A dead port is back at p_jn,
// f_eff lies between 1 and f_j (so >= 0, sign-safe), and f = 1 is untouched.
// A port carrying a share s of the throughput keeps 1 - delta^2/s^2 of its
// factor's effect: 96% at s = 5%, 99% at 10%.
//
// DENSITY. Thermodynamics-free like mpce_junction.h: rho and drho/dP per port
// enter as values, seeded first order on P_i.
//
// JACOBIAN columns: P[0..N-1], outer_mdot[0..N-1], p_jn. outer_mdot is the
// CONNECTING element's flow; port_sign maps it to junction convention.
// -----------------------------------------------------------------------------

#include <cmath>
#include <stdexcept>
#include <vector>

#include "dual_number.h"

namespace combaero::solver {

constexpr int kIdealMaxPorts = 6;
constexpr int kIdealSeeds = 2 * kIdealMaxPorts + 1;

// Width of the supplier/receiver role blend, as a fraction of the junction's
// throughput |m|_2. Not a physics parameter: it only smooths the switch.
constexpr double IDEAL_ROLE_EPS_REL = 1.0e-4;
// Absolute floor [kg/s] so the all-zero state stays smooth.
constexpr double IDEAL_ROLE_EPS_ABS = 1.0e-9;
// Flow share below which a receiver's loss factor fades back to 1 (#526),
// as a fraction of the throughput |m|_2. A smoothing width, not physics.
constexpr double IDEAL_LOSS_FADE_REL = 1.0e-2;

struct IdealJunctionGeometry {
  std::vector<double> area;       // [m^2]
  std::vector<double> axis_rad;   // outward port axis angle [rad]
  std::vector<double> port_sign;  // +1: connecting flow positive OUT of junction
  // Receiver loss factor f_i >= 0 per port (#526); empty = all 1, the
  // parameter-free closure.
  std::vector<double> loss_factor;
};

struct IdealJunctionResult {
  // N port rows [Pa] then the mass row [kg/s].
  std::vector<double> residual;
  // [row][col], columns P[0..N-1], outer_mdot[0..N-1], p_jn.
  std::vector<std::vector<double>> jacobian;
  // Signed port velocity u_i (positive out) [m/s] and the supplied mean
  // velocity vector V [m/s], for diagnostics.
  std::vector<double> velocity;
  double V_x = 0.0;
  double V_y = 0.0;
};

inline IdealJunctionResult ideal_junction_residuals_and_jacobian(
    const std::vector<double>& p_static, const std::vector<double>& rho,
    const std::vector<double>& drho_dp, const std::vector<double>& outer_mdot,
    double p_jn, const IdealJunctionGeometry& geom) {
  const int n = static_cast<int>(p_static.size());
  if (n < 2 || n > kIdealMaxPorts) {
    throw std::invalid_argument(
        "ideal_junction_residuals_and_jacobian: 2..6 ports required");
  }
  const auto sz = static_cast<std::size_t>(n);
  if (rho.size() != sz || drho_dp.size() != sz || outer_mdot.size() != sz ||
      geom.area.size() != sz || geom.axis_rad.size() != sz ||
      geom.port_sign.size() != sz) {
    throw std::invalid_argument(
        "ideal_junction_residuals_and_jacobian: inconsistent port counts");
  }
  if (!geom.loss_factor.empty()) {
    if (geom.loss_factor.size() != sz) {
      throw std::invalid_argument(
          "ideal_junction_residuals_and_jacobian: loss_factor needs one entry "
          "per port (or none)");
    }
    for (double f : geom.loss_factor) {
      if (!(f >= 0.0) || !std::isfinite(f)) {
        throw std::invalid_argument(
            "ideal_junction_residuals_and_jacobian: loss_factor must be finite "
            "and >= 0");
      }
    }
  }

  using D = DualN<kIdealSeeds>;
  std::vector<D> P(sz), M(sz), U(sz), Rho(sz);
  for (int i = 0; i < n; ++i) {
    P[i] = D::seed(p_static[i], i);
    Rho[i] = D::constant(rho[i]);
    Rho[i].d[i] = drho_dp[i];
    M[i] = D::seed(outer_mdot[i], kIdealMaxPorts + i) * geom.port_sign[i];
    U[i] = M[i] / (Rho[i] * geom.area[i]);
  }
  const D Pjn = D::seed(p_jn, 2 * kIdealMaxPorts);

  D m2 = D::constant(0.0);
  for (int i = 0; i < n; ++i) m2 = m2 + M[i] * M[i];
  const D eps2 = m2 * (IDEAL_ROLE_EPS_REL * IDEAL_ROLE_EPS_REL) +
                 IDEAL_ROLE_EPS_ABS * IDEAL_ROLE_EPS_ABS;

  std::vector<D> w_out(sz), w_in(sz), S(sz);
  D w_in_sum = D::constant(IDEAL_ROLE_EPS_ABS);
  D Vx = D::constant(0.0);
  D Vy = D::constant(0.0);
  for (int i = 0; i < n; ++i) {
    // s = smooth step 0 -> 1 across m = 0; s(m) + s(-m) = 1, so
    // w_out - w_in = m and both vanish at m = 0. With r = sqrt(m^2 + eps^2),
    // s = (r + m) / (2 r) and 1 - s = (r - m) / (2 r). The side whose
    // difference cancels (r - m for m > 0, r + m for m < 0; |m| >> eps away
    // from reversal) is formed as eps^2 / (r +- m), the same function to
    // round-off. Formed as m s - m, a receiver's w_in kept only
    // 1e-16 / (eps/m)^2 of its digits (4e-9 .. 2.4e-8 relative at
    // eps/m ~ 1e-4) and, with no supplier, V's numerator and denominator are
    // all such weights: 1e-10 of noise in a 1e5 Pa row, read by a finite
    // difference as a Jacobian miss (#1250).
    const D r = dsqrt(M[i] * M[i] + eps2);
    const D inv_2r = 0.5 / r;
    if (M[i].v >= 0.0) {
      S[i] = (r + M[i]) * inv_2r;
      w_in[i] = -M[i] * (eps2 / (r + M[i])) * inv_2r;
    } else {
      S[i] = (eps2 / (r - M[i])) * inv_2r;
      w_in[i] = -M[i] * (r - M[i]) * inv_2r;
    }
    w_out[i] = M[i] * S[i];
    const double cx = std::cos(geom.axis_rad[i]);
    const double cy = std::sin(geom.axis_rad[i]);
    w_in_sum = w_in_sum + w_in[i];
    Vx = Vx + w_in[i] * U[i] * cx;
    Vy = Vy + w_in[i] * U[i] * cy;
  }
  Vx = Vx / w_in_sum;
  Vy = Vy / w_in_sum;

  IdealJunctionResult out;
  out.residual.assign(sz + 1, 0.0);
  out.jacobian.assign(sz + 1, std::vector<double>((2 * sz) + 1, 0.0));
  auto store = [&](int row, const D& r) {
    out.residual[row] = r.v;
    for (int c = 0; c < n; ++c) {
      out.jacobian[row][c] = r.d[c];
      out.jacobian[row][n + c] = r.d[kIdealMaxPorts + c];
    }
    out.jacobian[row][2 * sz] = r.d[kIdealSeeds - 1];
  };

  for (int j = 0; j < n; ++j) {
    const D V_axial = Vx * std::cos(geom.axis_rad[j]) +
                      Vy * std::sin(geom.axis_rad[j]);
    D row = P[j] - Pjn + w_out[j] / geom.area[j] * (U[j] - V_axial);
    const double f = geom.loss_factor.empty() ? 1.0 : geom.loss_factor[j];
    if (f != 1.0) {
      const D dx = U[j] * std::cos(geom.axis_rad[j]) - Vx;
      const D dy = U[j] * std::sin(geom.axis_rad[j]) - Vy;
      const D L = Rho[j] * 0.5 * (dx * dx + dy * dy);
      const D delta2 = m2 * (IDEAL_LOSS_FADE_REL * IDEAL_LOSS_FADE_REL) +
                       IDEAL_ROLE_EPS_ABS * IDEAL_ROLE_EPS_ABS;
      const D phi = M[j] * M[j] / (M[j] * M[j] + delta2);
      row = row + S[j] * L * phi * (f - 1.0);
    }
    store(j, row);
  }
  D mass = D::constant(0.0);
  for (int i = 0; i < n; ++i) mass = mass + M[i];
  store(n, mass);

  out.velocity.resize(sz);
  for (int i = 0; i < n; ++i) out.velocity[i] = U[i].v;
  out.V_x = Vx.v;
  out.V_y = Vy.v;
  return out;
}

}  // namespace combaero::solver
