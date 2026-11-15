#include "side_stream_momentum.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "composition.h"
#include "compressible.h"
#include "inner_solve.h"
#include "input_validation.h"
#include "solver_interface.h"
#include "stagnation.h"
#include "thermo.h"

namespace combaero {

namespace {
// Floor on the impulse quadratic's discriminant, relative to Pi^2, so a
// choked face stays finite and differentiable (flagged, not hidden).
constexpr double MERGE_DISCRIMINANT_FLOOR = 1e-8;

// d(dP)/dY_k of a drop dP ~ 1/rho at fixed (P, T) (#1146): rho = P MW /
// (R_u T), so d(dP)/dX_j = -dP W_j / MW, mapped onto the mass fractions.
std::vector<double> inverse_density_drop_dY(double dP, const std::vector<double>& X) {
  const double MW = mwmix(X);
  std::vector<double> g(X.size(), 0.0);
  for (std::size_t j = 0; j < X.size(); ++j) {
    g[j] = -dP * species_molar_mass(j) / MW;
  }
  return mole_gradient_to_mass(X, g);
}
}  // namespace

StationHalfDrop station_half_drop(double m_a, double m_b, double P, double T,
                                  const std::vector<double>& X, double area,
                                  double kappa) {
  if (!(area > 0.0) || !(P > 0.0) || !(T > 0.0)) {
    throw std::invalid_argument("station_half_drop: area, P and T must be positive");
  }
  auto [rho, drho_dT, drho_dP] = solver::density_and_jacobians(T, P, X);
  const double k = 0.5 / (rho * area * area);
  // Momentum flux m u = m^2 / (rho A) is positive whichever way the gas
  // moves: m^2, not m|m| (#815; m|m| flipped the drop for reversed flow).
  // The side stream's axial velocity follows the ARRIVING flow (#1178):
  // s = M_in / A_in = max(m_a, 0) - max(-m_b, 0), the mass-weighted mean
  // signed velocity of the streams arriving along the channel, times rho A
  // (M_in = max(m_a,0)^2 - max(-m_b,0)^2 and A_in = max(m_a,0) + max(-m_b,0)
  // share the factor A_in; s = 0 when nothing arrives). Forward s = m_a.
  // C0 across m_a = 0 and m_b = 0. Exactly at a joint the partials are the
  // right-sided ones (the flow at +0) for both flows, so a segment flow that
  // is m_b at one station and m_a at the next reads one side at both.
  const double s = std::max(m_a, 0.0) - std::max(-m_b, 0.0);
  const double ds_da = m_a >= 0.0 ? 1.0 : 0.0;
  const double ds_db = m_b < 0.0 ? 1.0 : 0.0;
  // Written so that forward flow (s = m_a, ds_da = 1, ds_db = 0) evaluates the
  // same operations as the m_a form it replaces, bit for bit.
  const double f = m_b * m_b - m_a * m_a + kappa * s * (m_a - m_b);
  StationHalfDrop out;
  out.dP = f * k;
  out.d_dm_a = (-2.0 * m_a + kappa * ((s + ds_da * m_a) - ds_da * m_b)) * k;
  out.d_dm_b = (2.0 * m_b + kappa * (ds_db * (m_a - m_b) - s)) * k;
  const double d_drho = -out.dP / rho;
  out.d_dP = d_drho * drho_dP;
  out.d_dT = d_drho * drho_dT;
  out.d_dY = inverse_density_drop_dY(out.dP, X);
  return out;
}

StationHalfDrop segment_area_change_drop(double m, double P, double T,
                                         const std::vector<double>& X, double area_from,
                                         double area_to) {
  if (!(area_from > 0.0) || !(area_to > 0.0) || !(P > 0.0) || !(T > 0.0)) {
    throw std::invalid_argument(
        "segment_area_change_drop: areas, P and T must be positive");
  }
  StationHalfDrop out;
  if (area_from == area_to) {
    out.d_dY.assign(X.size(), 0.0);
    return out;
  }
  auto [rho, drho_dT, drho_dP] = solver::density_and_jacobians(T, P, X);
  const double c = 0.5 * (1.0 / (area_to * area_to) - 1.0 / (area_from * area_from)) / rho;
  out.dP = m * m * c;
  out.d_dm_a = 2.0 * m * c;
  const double d_drho = -out.dP / rho;
  out.d_dP = d_drho * drho_dP;
  out.d_dT = d_drho * drho_dT;
  out.d_dY = inverse_density_drop_dY(out.dP, X);
  return out;
}

StationHalfDrop channel_entry_drop(double m, double P, double T,
                                   const std::vector<double>& X, double area,
                                   double K_in) {
  if (!(area > 0.0) || !(P > 0.0) || !(T > 0.0)) {
    throw std::invalid_argument("channel_entry_drop: area, P and T must be positive");
  }
  auto [rho, drho_dT, drho_dP] = solver::density_and_jacobians(T, P, X);
  const double k = 0.5 * (1.0 + K_in) / (rho * area * area);
  // Only an entering flow accelerates out of the reservoir; a reversed one
  // discharges into it at static continuity (see the header).
  const double m_in = std::max(m, 0.0);
  StationHalfDrop out;
  out.dP = m_in * m_in * k;
  out.d_dm_a = 2.0 * m_in * k;  // the entering flow is "m_a"
  const double d_drho = -out.dP / rho;
  out.d_dP = d_drho * drho_dP;
  out.d_dT = d_drho * drho_dT;
  out.d_dY = inverse_density_drop_dY(out.dP, X);
  return out;
}

MergeFaceState chamber_merge_face_state(double m_main, double Tt_main,
                                        const std::vector<double>& X_main,
                                        double m_out, double P, double Tt,
                                        const std::vector<double>& X,
                                        double side_momentum, double area,
                                        bool incompressible) {
  if (!(area > 0.0) || !(P > 0.0) || !(Tt > 0.0) || !(Tt_main > 0.0)) {
    throw std::invalid_argument(
        "chamber_merge_face_state: area, P, Tt and Tt_main must be positive");
  }

  if (incompressible) {
    // Low-Mach limit (#588): rho = P / (R Tt) for both streams, at the
    // chamber pressure. Linear in P_f; every sensitivity is closed form.
    MergeFaceState out;
    const double R = specific_gas_constant(X);
    const double Rm = specific_gas_constant(X_main);
    const double A2 = area * area;
    const double a = m_main * m_main / A2;  // main impulse factor
    const double b = m_out * m_out / A2;    // outlet impulse factor
    const double vm = Rm * Tt_main / P;     // 1 / rho_main
    const double vo = R * Tt / P;           // 1 / rho
    out.P_face = P + b * vo - side_momentum / area - a * vm;
    out.Pt_face = out.P_face + 0.5 * a * vm;
    out.M_face = std::abs(m_main) * vm / area / speed_of_sound(Tt_main, X_main);
    out.choked = false;
    out.dPf_dm_main = -2.0 * m_main * vm / A2;
    out.dPf_dm_out = 2.0 * m_out * vo / A2;
    out.dPf_dJ = -1.0 / area;
    out.dPf_dP = 1.0 - b * vo / P + a * vm / P;
    out.dPf_dT = b * R / P;
    out.dPf_dT_main = -a * Rm / P;
    out.dPtf_dm_main = -m_main * vm / A2;
    out.dPtf_dm_out = out.dPf_dm_out;
    out.dPtf_dJ = out.dPf_dJ;
    out.dPtf_dP = 1.0 - b * vo / P + 0.5 * a * vm / P;
    out.dPtf_dT = out.dPf_dT;
    out.dPtf_dT_main = -0.5 * a * Rm / P;
    // Composition: R = R_u / MW in both specific volumes, dR/dX_j = -R W_j / MW.
    {
      const std::size_t n = X.size();
      const double MW = mwmix(X) / 1000.0;
      const double MWm = mwmix(X_main) / 1000.0;
      std::vector<double> gPf(n), gPfm(n), gPtm(n);
      for (std::size_t j = 0; j < n; ++j) {
        const double W = species_molar_mass(j) / 1000.0;
        gPf[j] = -b * vo * W / MW;
        gPfm[j] = a * vm * W / MWm;
        gPtm[j] = 0.5 * a * vm * W / MWm;
      }
      out.dPf_dY = mole_gradient_to_mass(X, gPf);
      out.dPtf_dY = out.dPf_dY;
      out.dPf_dY_main = mole_gradient_to_mass(X_main, gPfm);
      out.dPtf_dY_main = mole_gradient_to_mass(X_main, gPtm);
    }
    return out;
  }
  const double A = area;
  const double J = side_momentum;

  // Chamber (outlet) side: static T from the chamber's stagnation Tt and the
  // outlet mass flux (#578), the same solve as the chamber's own closure.
  const double R = specific_gas_constant(X);
  const double s_out = m_out >= 0.0 ? 1.0 : -1.0;
  const ChamberStaticState cs = chamber_static_state(Tt, P, std::abs(m_out) / A, X);
  // Outlet momentum flux per area, m|m| / (rho A^2) = m|m| R Ts / (P A^2).
  const double q = m_out * std::abs(m_out) * R * cs.T / (P * A * A);
  const double Pi = P + q - J / A;
  const double dPi_dP = 1.0 + q * (cs.dT_dP / cs.T - 1.0 / P);
  const double dPi_dmout =
      2.0 * std::abs(m_out) * R * cs.T / (P * A * A) + q * (cs.dT_dG * s_out / A) / cs.T;
  const double dPi_dTt = q * cs.dT_dTt / cs.T;
  const double dPi_dJ = -1.0 / A;

  // Main face: static (P_f, T_f) from the impulse and the main stream's
  // stagnation Tt_main (#578; T_main was read as static before):
  //   F1 = P_f + c T_f / P_f - Pi = 0,   c = m_main^2 R_m / A^2
  //   F2 = T_f - Ts(Tt_main, P_f, |m_main| / A) = 0
  // Fixed point from T_f = Tt_main: the quadratic in P_f on its subsonic
  // root, then the static T at that pressure. Contracts on the subsonic
  // branch (the face Mach is well below 1 wherever the quadratic has a root).
  const double Rm = specific_gas_constant(X_main);
  const double c = m_main * m_main * Rm / (A * A);
  const double Gm = std::abs(m_main) / A;
  double T_f = Tt_main;
  double P_f = Pi;
  bool choked = false;
  ChamberStaticState fs{};
  for (int it = 0;; ++it) {
    const double D_raw = Pi * Pi - 4.0 * c * T_f;
    const double D_min = MERGE_DISCRIMINANT_FLOOR * Pi * Pi;
    choked = D_raw < D_min;
    P_f = 0.5 * (Pi + std::sqrt(std::max(D_raw, D_min)));
    fs = chamber_static_state(Tt_main, P_f, Gm, X_main);
    const double dT = fs.T - T_f;
    T_f = fs.T;
    if (std::abs(dT) <= 1e-12 * Tt_main) {
      break;
    }
    if (it + 1 == inner_solve::kMergeFaceMaxIter) {
      if (inner_solve::cap_exit_is_root(dT, Tt_main)) {
        break;
      }
      // A fixed point that has not contracted: (P_f, T_f) is not the face
      // state, and the implicit derivatives below would describe a non-root.
      inner_solve::throw_no_convergence("merge face state (side_stream_momentum)",
                                        inner_solve::kMergeFaceMaxIter, "T step [K]", dT);
    }
  }
  MergeFaceState out;
  out.P_face = P_f;
  out.choked = choked;
  const double a_f = speed_of_sound(T_f, X_main);
  const double u_f = Gm * Rm * T_f / P_f;
  out.M_face = a_f > 0.0 ? u_f / a_f : 0.0;
  // Pt_f = P_f exp((s(Tt_main) - s(T_f)) / R_m), consistent with the main
  // stream's own stagnation state; the entropy difference is the integrated
  // one from the static-state solve (no cancellation, #579).
  const double Pt_f = P_f * std::exp(fs.lnPtP);
  out.Pt_face = Pt_f;

  // Sensitivities of (P_f, T_f) by the implicit function theorem on
  // (F1, F2); with the discriminant floored, P_f follows Pi alone.
  // M = [[F1_Pf, F1_Tf], [F2_Pf, 1]].
  double dPf_dPi, dPf_dc, dPf_dTtm, dPf_dGm;
  if (choked) {
    dPf_dPi = 0.5 * (1.0 + std::sqrt(MERGE_DISCRIMINANT_FLOOR));
    dPf_dc = dPf_dTtm = dPf_dGm = 0.0;
  } else {
    const double F1_Pf = 1.0 - c * T_f / (P_f * P_f);
    const double F1_Tf = c / P_f;
    const double F2_Pf = -fs.dT_dP;
    const double det = F1_Pf - F1_Tf * F2_Pf;
    // d(P_f) = -(F1_x - F1_Tf * F2_x) / det for an input x with partials
    // (F1_x, F2_x); F2 has unit T_f coefficient.
    auto dPf = [&](double F1_x, double F2_x) { return -(F1_x - F1_Tf * F2_x) / det; };
    dPf_dPi = dPf(-1.0, 0.0);
    dPf_dc = dPf(T_f / P_f, 0.0);
    dPf_dTtm = dPf(0.0, -fs.dT_dTt);
    dPf_dGm = dPf(0.0, -fs.dT_dG);
  }
  // T_f = Ts(Tt_main, P_f, Gm): its total derivative per input follows P_f.
  const double dTf_dPf = fs.dT_dP;
  const double cp_f = cp_mass(T_f, X_main);
  const double cp_0 = cp_mass(Tt_main, X_main);
  // dPt_f = (Pt_f / P_f) dP_f + Pt_f dL, dL = -cp_f/(R T_f) dT_f + cp_0/(R Tt_main) dTt_main.
  auto dPt = [&](double dPf_dx, double dTf_dx, double dTtm_dx) {
    return (Pt_f / P_f) * dPf_dx +
           Pt_f * (-cp_f / (Rm * T_f) * dTf_dx + cp_0 / (Rm * Tt_main) * dTtm_dx);
  };

  const double sgn_m = m_main >= 0.0 ? 1.0 : -1.0;
  const double dc_dmm = 2.0 * m_main * Rm / (A * A);
  const double dGm_dmm = sgn_m / A;

  // Inputs that reach the face only through Pi.
  auto through_Pi = [&](double dPi_dx, double& dPf_out, double& dPt_out) {
    dPf_out = dPf_dPi * dPi_dx;
    dPt_out = dPt(dPf_out, dTf_dPf * dPf_out, 0.0);
  };
  through_Pi(dPi_dP, out.dPf_dP, out.dPtf_dP);
  through_Pi(dPi_dmout, out.dPf_dm_out, out.dPtf_dm_out);
  through_Pi(dPi_dJ, out.dPf_dJ, out.dPtf_dJ);
  through_Pi(dPi_dTt, out.dPf_dT, out.dPtf_dT);

  // m_main: through c and through the face mass flux Gm.
  out.dPf_dm_main = dPf_dc * dc_dmm + dPf_dGm * dGm_dmm;
  {
    const double dTf = dTf_dPf * out.dPf_dm_main + fs.dT_dG * dGm_dmm;
    out.dPtf_dm_main = dPt(out.dPf_dm_main, dTf, 0.0);
  }
  // Tt_main: through T_f (F2) and through s(Tt_main) in Pt_f.
  out.dPf_dT_main = dPf_dTtm;
  {
    const double dTf = dTf_dPf * dPf_dTtm + fs.dT_dTt;
    out.dPtf_dT_main = dPt(dPf_dTtm, dTf, 1.0);
  }

  // Composition (#743). Outlet gas X: through Pi (R and the outlet static T
  // in q). Main X: through c = m_main^2 R_m / A^2 (F1), the face static T
  // (F2) and the entropy difference in Pt_f.
  {
    const std::size_t n = X.size();
    const double MW = mwmix(X) / 1000.0;
    const double MWm = mwmix(X_main) / 1000.0;
    const ChamberStaticStateDX cx = chamber_static_state_dX(Tt, P, std::abs(m_out) / A, X);
    const ChamberStaticStateDX fx = chamber_static_state_dX(Tt_main, P_f, Gm, X_main);
    std::vector<double> gPf(n), gPt(n), gPfm(n), gPtm(n);
    for (std::size_t j = 0; j < n; ++j) {
      const double W = species_molar_mass(j) / 1000.0;
      through_Pi(q * (-W / MW + cx.dT_dX[j] / cs.T), gPf[j], gPt[j]);
      const double dc = -c * W / MWm;
      double dPf_x = 0.0;
      if (!choked) {
        const double F1_Pf = 1.0 - c * T_f / (P_f * P_f);
        const double F1_Tf = c / P_f;
        const double det = F1_Pf - F1_Tf * (-fs.dT_dP);
        // F1_x = (T_f / P_f) dc, F2_x = -dTs/dX at fixed P_f.
        dPf_x = -((T_f / P_f) * dc - F1_Tf * (-fx.dT_dX[j])) / det;
      }
      const double dTf = dTf_dPf * dPf_x + fx.dT_dX[j];
      gPfm[j] = dPf_x;
      // ln(P0/P) at the face moves with X directly and through T_f.
      const double dL = fx.dlnPtP_dX[j] - fx.dT_dX[j] * (-cp_f / (Rm * T_f)) +
                        (-cp_f / (Rm * T_f)) * dTf;
      gPtm[j] = (Pt_f / P_f) * dPf_x + Pt_f * dL;
    }
    out.dPf_dY = mole_gradient_to_mass(X, gPf);
    out.dPtf_dY = mole_gradient_to_mass(X, gPt);
    out.dPf_dY_main = mole_gradient_to_mass(X_main, gPfm);
    out.dPtf_dY_main = mole_gradient_to_mass(X_main, gPtm);
  }
  return out;
}

namespace {
// w(Pt, Tt, P) and its partials, analytic (#816). The gas is ideal with
// cp(T): along the isentrope s0(T) - R ln P = s0(Tt) - R ln Pt, so
//   (cp(T)/T) dT = (cp(Tt)/Tt) dTt - R dPt/Pt + R dP/P   (mass-specific).
struct JetW {
  double w = 0.0;
  bool choked = false;
  double dw_dPt = 0.0, dw_dTt = 0.0, dw_dP = 0.0;
  // d/dX_j at fixed (Pt, Tt, P), on the extension of the mixture properties
  // linear in X (h, s0, cp molar) with MW = sum X W (#1146). Empty: no jet.
  std::vector<double> dw_dX;
};

JetW jet_w(double Pt, double Tt, double P, const std::vector<double>& X) {
  JetW out;
  if (!(Pt > 0.0) || !(Tt > 0.0) || !(P > 0.0) || P >= Pt) {
    return out;
  }
  const double R = specific_gas_constant(X);
  const double cp_t = cp_mass(Tt, X);
  // The sonic state ON the isentrope: T* from h(Tt) - h(T*) = a(T*)^2 / 2 (a
  // smooth Newton solve), P* = Pt exp((s0(T*) - s0(Tt)) / R). The golden-
  // section mass-flux maximum (critical_pressure_ratio) locates the same
  // point only to ~1e-8 of Pt on a flat peak, noise the choked w carried.
  const double T_star = T_from_stagnation(Tt, 1.0, X);
  const double P_star =
      Pt * std::exp((s_mass(T_star, X, Pt) - s_mass(Tt, X, Pt)) / R);
  if (P >= P_star) {
    const double M = mach_from_pressure_ratio(Tt, Pt, P, X);
    const double Ts = T_from_stagnation(Tt, M, X);
    out.w = M * speed_of_sound(Ts, X);
    if (!(out.w > 0.0)) {
      return out;
    }
    // w^2 = 2 (h(Tt) - h(Ts)): dw = (cp(Tt) dTt - cp(Ts) dTs) / w, with dTs
    // from the isentrope: dTs = (Ts/cp(Ts)) ((cp(Tt)/Tt) dTt - R dPt/Pt +
    // R dP/P), so cp(Ts) dTs = Ts (...).
    out.dw_dTt = (cp_t - Ts * cp_t / Tt) / out.w;
    out.dw_dPt = Ts * R / Pt / out.w;
    out.dw_dP = -Ts * R / P / out.w;
    // Composition (#1146), molar: the isentrope sum X_j (s0_j(Ts) - s0_j(Tt))
    // = R_u ln(P / Pt) moves Ts by dTs/dX_j = -(s0_j(Ts) - s0_j(Tt)) Ts / c(Ts);
    // w^2 / 2 = (H(Tt) - H(Ts)) / MW with H = sum X_j h_j.
    {
      const std::size_t n = X.size();
      const double MW = mwmix(X);
      const double MWkg = MW / 1000.0;
      const double c_s = cp(Ts, X);
      const double half_w2 = 0.5 * out.w * out.w;
      out.dw_dX.assign(n, 0.0);
      for (std::size_t j = 0; j < n; ++j) {
        const double dTs = -(s_species(j, Ts) - s_species(j, Tt)) * Ts / c_s;
        const double dq = (h_species(j, Tt) - h_species(j, Ts) - c_s * dTs) / MWkg -
                          half_w2 * species_molar_mass(j) / MW;
        out.dw_dX[j] = dq / out.w;
      }
    }
    return out;
  }
  // Sonic: h(Tt) - h(T*) = a(T*)^2 / 2 fixes T*(Tt); P* = Pt exp(-(s0(Tt) -
  // s0(T*))/R) from the isentrope; w = u* + (P* - P) / (rho* u*).
  out.choked = true;
  const double u_star = speed_of_sound(T_star, X);
  const double rho_star = density(T_star, P_star, X);
  out.w = u_star + (P_star - P) / (rho_star * u_star);
  // d(a^2)/dT = R (gamma + T dgamma/dT), gamma = c/(c - R_mol) on molar cp c.
  const double c = cp(T_star, X);
  const double gamma = c / (c - thermo::R_GAS);
  const double dgamma_dT = -thermo::R_GAS * dcp_dT(T_star, X) / ((c - thermo::R_GAS) * (c - thermo::R_GAS));
  const double da2_dT = R * (gamma + T_star * dgamma_dT);
  const double cp_s = cp_mass(T_star, X);
  const double dTs_dTt = cp_t / (cp_s + 0.5 * da2_dT);
  const double du_dTt = da2_dT / (2.0 * u_star) * dTs_dTt;
  // ln P* = ln Pt + (s0(T*) - s0(Tt)) / R.
  const double dlnPs_dTt = (cp_s / T_star * dTs_dTt - cp_t / Tt) / R;
  // w = u* + (P* - P) K, K = R T* / (P* u*) = 1 / (rho* u*).
  const double K = R * T_star / (P_star * u_star);
  const double dPs_dTt = P_star * dlnPs_dTt;
  const double dK_dTt = K * (dTs_dTt / T_star - dlnPs_dTt - du_dTt / u_star);
  out.dw_dP = -K;
  out.dw_dPt = K * P / Pt;  // P* ~ Pt at fixed Tt: d((P* - P) K)/dPt = K P / Pt
  out.dw_dTt = du_dTt + dPs_dTt * K + (P_star - P) * dK_dTt;
  // Composition (#1146), each step at fixed (Pt, Tt, P): T* from
  // E = (H(Tt) - H(T*)) / MW - a^2(T*) / 2 = 0, a^2 = gamma R T,
  // gamma = c / (c - R_u), R = R_u / MW; then ln P*, u*, K as above.
  {
    const std::size_t n = X.size();
    const double MW = mwmix(X);
    const double MWkg = MW / 1000.0;
    const double dH = h(Tt, X) - h(T_star, X);  // [J/mol]
    const double dE_dTs = cp_s + 0.5 * da2_dT;   // -dE/dT*
    out.dw_dX.assign(n, 0.0);
    for (std::size_t j = 0; j < n; ++j) {
      const double Wj = species_molar_mass(j);
      const double cj = cp_species(j, T_star);
      const double dgamma = -thermo::R_GAS * cj / ((c - thermo::R_GAS) * (c - thermo::R_GAS));
      const double da2 = T_star * (dgamma * R - gamma * R * Wj / MW);
      const double dE = (h_species(j, Tt) - h_species(j, T_star)) / MWkg -
                        dH / MWkg * Wj / MW - 0.5 * da2;
      const double dTs = dE / dE_dTs;
      const double dlnPs =
          (s_species(j, T_star) - s_species(j, Tt) + c / T_star * dTs) / thermo::R_GAS;
      const double du = (da2 + da2_dT * dTs) / (2.0 * u_star);
      const double dlnK = -Wj / MW + dTs / T_star - dlnPs - du / u_star;
      out.dw_dX[j] = du + P_star * dlnPs * K + (P_star - P) * K * dlnK;
    }
  }
  return out;
}
}  // namespace

JetImpulse jet_impulse(double m, double Pt, double Tt, double P,
                       const std::vector<double>& X) {
  // jet_w's "no jet" branch is for finite states; NaN used to answer J = 0.
  if (!validation::is_finite(m) || !validation::is_finite(Pt) || !validation::is_finite(Tt) ||
      !validation::is_finite(P)) {
    throw std::invalid_argument("jet_impulse: m, Pt, Tt and P must be finite");
  }
  const JetW base = jet_w(Pt, Tt, P, X);
  JetImpulse out;
  out.J = m * base.w;
  out.dJ_dm = base.w;
  out.choked = base.choked;
  out.dJ_dPt = m * base.dw_dPt;
  out.dJ_dTt = m * base.dw_dTt;
  out.dJ_dP = m * base.dw_dP;
  if (base.dw_dX.empty()) {
    out.dJ_dY.assign(X.size(), 0.0);
  } else {
    out.dJ_dY = mole_gradient_to_mass(X, base.dw_dX);
    for (double& d : out.dJ_dY) {
      d *= m;
    }
  }
  return out;
}

CrossflowRatio crossflow_velocity_ratio(double m_a, double m_b, double P,
                                        double T, const std::vector<double>& X,
                                        double P_down, double area) {
  if (!(area > 0.0) || !(P > 0.0) || !(T > 0.0)) {
    throw std::invalid_argument(
        "crossflow_velocity_ratio: area, P and T must be positive");
  }
  auto [rho, drho_dT, drho_dP] = solver::density_and_jacobians(T, P, X);
  const double m_mean = 0.5 * (m_a + m_b);
  const double sgn = m_mean >= 0.0 ? 1.0 : -1.0;
  const double U1 = std::abs(m_mean) / (rho * area);
  // Unit-flow jet impulse: J = w, the isentropic velocity from (P, T) to P_down.
  const JetImpulse jet = jet_impulse(1.0, P, T, P_down, X);
  const double Vi = std::sqrt(jet.J * jet.J + CROSSFLOW_VI_FLOOR * CROSSFLOW_VI_FLOOR);
  const double dVi_dw = Vi > 0.0 ? jet.J / Vi : 0.0;

  CrossflowRatio out;
  out.U1 = U1;
  out.Vi = Vi;
  out.U1_over_Vi = U1 / Vi;
  const double dU1_dm = 0.5 * sgn / (rho * area);
  out.d_dm_a = dU1_dm / Vi;
  out.d_dm_b = dU1_dm / Vi;
  const double dU1_drho = -U1 / rho;
  const double r_over_Vi = out.U1_over_Vi / Vi;
  out.d_dP = dU1_drho * drho_dP / Vi - r_over_Vi * dVi_dw * jet.dJ_dPt;
  out.d_dT = dU1_drho * drho_dT / Vi - r_over_Vi * dVi_dw * jet.dJ_dTt;
  out.d_dP_down = -r_over_Vi * dVi_dw * jet.dJ_dP;
  return out;
}

}  // namespace combaero
