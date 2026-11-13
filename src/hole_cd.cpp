#include "hole_cd.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "solver_interface.h"

// The hole-Cd chain (see hole_cd.h). The operation order mirrors the Python
// composition it replaced (python/combaero/network/hole_cd.py before this
// port), so that the values are bit-identical: Cd0 from the correlation's
// value function, its partials from Cd_and_derivatives, the exit term's
// chain rule dCd/dx = dCd/dCd0 dCd0/dx (+ dCd/dk dk/dx where dCd/dk != 0).

namespace combaero {

namespace {

// Term 1: Cd0 and its own partials.
void compose_base(const HoleCdModel& model, const DischargeHoleGeometry& hole,
                  const OrificeGeometry& orifice, const HoleCdInputs& in, HoleCdResult& out) {
  switch (model.base) {
    case HoleCdBase::Fixed: {
      const double v = model.fixed_cd < 1.0 ? model.fixed_cd : 1.0;
      out.cd_base = HOLE_CD_FIXED_MIN < v ? v : HOLE_CD_FIXED_MIN;
      return;
    }
    case HoleCdBase::Metering: {
      const auto [cd, dRe_D] = metering_cd_and_dRe(model.metering, orifice, in.Re_D);
      out.cd_base = cd;
      out.dcd_base_dRe_D = dRe_D;
      return;
    }
    case HoleCdBase::Discharge: {
      const std::unique_ptr<DischargeCorrelationBase> c =
          make_discharge_correlation(model.discharge);
      DischargeHoleState flow;
      flow.Re = in.Re;
      flow.U1_over_Vi = in.U1_over_Vi;
      out.cd_base = c->Cd(hole, flow);
      const auto [cd, dRe, dU] = c->Cd_and_derivatives(hole, flow);
      (void)cd;
      out.re_held = in.Re < discharge_cd_re_floor(model.discharge);
      out.dcd_base_dRe = (out.re_held && !model.continue_re_floor) ? 0.0 : dRe;
      out.dcd_base_dU1_over_Vi = dU;
      out.inclination_dzeta = discharge_inclination_dzeta(hole.inclination_deg);
      return;
    }
  }
}

// Term 2: the opt-in exit-crossflow loss on Cd0.
void compose_exit(const HoleCdModel& model, const HoleCdInputs& in, HoleCdResult& out) {
  out.cd_hole = out.cd_base;
  out.dCd_dRe = out.dcd_base_dRe;
  out.dCd_dU1_over_Vi = out.dcd_base_dU1_over_Vi;
  out.dCd_dRe_D = out.dcd_base_dRe_D;
  if (!model.exit_crossflow) {
    return;
  }
  const Rowbury2001Set set = model.exit_set;
  const double th = model.exit_theta_deg;
  const double al = model.exit_alpha_deg;
  out.exit_k_floor = rowbury_2001_k_floor(set, th, al);
  out.exit_floor_above_fit = out.exit_k_floor > rowbury_2001_k_fit_max(set);
  out.exit_in_range = rowbury_2001_in_range(set, th, al);
  if (in.rho_m_u_m2 > 0.0) {
    out.exit_active = true;
    out.k_perp = rowbury_2001_k_perp(in.G_h, in.rho_h, in.rho_m_u_m2, th, al);
    out.exit = rowbury_2001_exit_crossflow_cd(out.cd_base, out.k_perp.k, set, th, al);
    out.cd_hole = out.exit.Cd;
    const double a = out.exit.dCd_dCd0;
    out.dCd_dRe = a * out.dcd_base_dRe;
    out.dCd_dU1_over_Vi = a * out.dcd_base_dU1_over_Vi;
    out.dCd_dRe_D = a * out.dcd_base_dRe_D;
    if (out.exit.dCd_dk != 0.0) {
      const double b = out.exit.dCd_dk;
      out.dCd_dG_h = b * out.k_perp.dk_dG;
      out.dCd_drho_h = b * out.k_perp.dk_drho_h;
      out.dCd_drho_m_u_m2 = b * out.k_perp.dk_dq_m;
    }
    out.exit_k_tier = out.exit.faded                                      ? 1.0
                      : out.exit.k_used > rowbury_2001_k_data_max(set) ? 0.5
                                                                         : 0.0;
  }
  out.exit_ratio = out.cd_base > 0.0 ? out.cd_hole / out.cd_base
                                     : std::numeric_limits<double>::quiet_NaN();
}

// Term 3: the impingement jet-crossflow law, which replaces the hole's Cd
// where it acts.
void compose_jet(const HoleCdModel& model, const HoleCdInputs& in, HoleCdResult& out) {
  out.Cd = out.cd_hole;
  if (model.jet == HoleCdJetCrossflow::Off) {
    return;
  }
  CrossflowCd r;
  if (model.jet == HoleCdJetCrossflow::RogersHersh1976) {
    r = rogers_hersh_1976_cd(out.cd_hole, in.G_ideal, in.G_c, in.suction);
    out.jet_capped = r.capped;
    if (r.capped) {
      return;  // the hole's own Cd and partials hold
    }
  } else {
    r = florschuetz_metzger_su_1982_cd(model.xn_d, model.yn_d, model.z_d, in.G_ideal, in.G_c);
    out.jet_capped = r.capped;
  }
  out.Cd = r.Cd;
  out.jet_acts = true;
  out.dCd_dRe = 0.0;
  out.dCd_dU1_over_Vi = 0.0;
  out.dCd_dRe_D = 0.0;
  out.dCd_dG_h = 0.0;
  out.dCd_drho_h = 0.0;
  out.dCd_drho_m_u_m2 = 0.0;
  // The table's constant, or no drive: Cd ~ G_ideal^1.15 has no slope at
  // G_ideal = 0 that the flow (~ Cd sqrt(dP)) would feel.
  out.jet_partials = !r.capped && in.G_ideal > 0.0;
  if (out.jet_partials) {
    out.dCd_dG_ideal = r.dCd_dG_ideal;
    out.dCd_dG_c = r.dCd_dG_c;
  }
}

}  // namespace

HoleCdResult hole_cd_compose(const HoleCdModel& model, const DischargeHoleGeometry& hole,
                             const OrificeGeometry& orifice, const HoleCdInputs& inputs) {
  HoleCdResult out;
  compose_base(model, hole, orifice, inputs, out);
  compose_exit(model, inputs, out);
  compose_jet(model, inputs, out);
  return out;
}

HoleSupplyDriveResult hole_supply_drive(HoleSupplyDrive drive, double P, double Tt,
                                        const std::vector<double>& Y, double m_dot, double area,
                                        bool incompressible) {
  // The station's dynamic head from the momentum-chamber closure, read as a
  // residual Pt - Pt_cl(P, Tt, m_dot, Y) at Pt = 0: Pt_cl = P + q, and its
  // partials are the residual's, negated.
  const solver::MomentumChamberResult c =
      solver::momentum_chamber_residual_and_jacobian(P, 0.0, m_dot, Tt, Y, area, incompressible);
  HoleSupplyDriveResult out;
  // The closure's own value, -residual at Pt = 0: P + q on the subsonic
  // branch, and past the station's sonic static the closure's bounded
  // continuation 2 Pt* - Pt_cl(P) (#788), so value and partials (the
  // residual's, negated) stay one function on both branches.
  const double pt_station = -c.residual;
  out.q_dyn = pt_station - P;
  if (drive == HoleSupplyDrive::Static) {
    out.P_drive = P;
    out.dP_dP = 1.0;
    out.dP_dY.assign(Y.size(), 0.0);
    return out;
  }
  out.P_drive = pt_station;
  out.dP_dP = -c.d_res_dP;
  out.dP_dTt = -c.d_res_dT;
  out.dP_dm_dot = -c.d_res_dmdot;
  out.dP_dY.resize(c.d_res_dY.size());
  for (std::size_t k = 0; k < c.d_res_dY.size(); ++k) {
    out.dP_dY[k] = -c.d_res_dY[k];
  }
  return out;
}

HoleDriveCdRatio hole_supply_recovery_cd_ratio(double P1, double T1, double Pt1, double Tt1,
                                               double P2, const std::vector<double>& Y, double f,
                                               bool incompressible) {
  if (!(f >= 0.0 && f <= 1.0)) {
    throw std::invalid_argument("hole_supply_recovery_cd_ratio: the recovered fraction f = " +
                                std::to_string(f) + " must lie in [0, 1]");
  }
  HoleDriveCdRatio out;
  out.d_dY.assign(Y.size(), 0.0);
  // The partly recovered drive and its partials: Pf = P1 + f q1.
  const double q1 = Pt1 - P1;
  // f = 1 takes Pt1 itself, so both ends are exact (ratio 1 bit for bit).
  const double Pf = (f == 1.0) ? Pt1 : P1 + f * q1;
  const double dPf_dP1 = 1.0 - f;
  const double dPf_dPt1 = f;
  const double dPf_df = q1;
  if (!(Pf > P2) || !(Pt1 > P2)) {
    return out;
  }
  if (incompressible) {
    // r = sqrt(a / b), a = Pf - P2, b = Pt1 - P2: one density on both sides,
    // so the ratio of the two ideal flows is the ratio of their pressure
    // drops' square roots.
    const double a = Pf - P2;
    const double b = Pt1 - P2;
    const double r = std::sqrt(a / b);
    const double dr_da = 0.5 * r / a;
    out.ratio = r;
    out.d_dP1 = dr_da * dPf_dP1;
    out.d_dPt1 = dr_da * dPf_dPt1 + (-0.5 * r / b);
    out.d_dP2 = 0.5 * r * (1.0 / b - 1.0 / a);
    out.d_df = dr_da * dPf_df;
    return out;
  }
  // Tf = T1 exp(s M), M = ln(Tt1/T1), s = ln(Pf/P1)/L = log1p(f x)/log1p(x),
  // x = q1/P1: the isentrope through the static and stagnation states (see
  // the header). log1p keeps s exact as q1 -> 0; its x-derivative
  // ds/dx = (f/(1 + f x) - s/(1 + x))/L cancels to O(x) there, so below
  // |x| = 1e-5 it is the series f(1-f)(12 - 4x - 16fx + 3x^2 + 6fx^2 +
  // 18f^2x^2)/24 (truncation O(x^3)), and ds/df = x/((1 + f x) L) likewise
  // 1 + x/2 - fx - x^2/12 - fx^2/2 + f^2x^2. At f = 0, s and ds/dx are 0
  // exactly, so the f = 0 ratio and partials are the static-to-total ones
  // bit for bit.
  const double M = (T1 > 0.0 && Tt1 > 0.0) ? std::log(Tt1 / T1) : 0.0;
  const double x = q1 / P1;
  const double L = std::log1p(x);
  const double s = (L != 0.0) ? std::log1p(f * x) / L : f;
  double ds_dx = 0.0, ds_df = 0.0;
  if (std::abs(x) < 1e-5) {
    ds_dx = f * (1.0 - f) *
            (12.0 - 4.0 * x - 16.0 * f * x + 3.0 * x * x + 6.0 * f * x * x +
             18.0 * f * f * x * x) /
            24.0;
    ds_df = 1.0 + 0.5 * x - f * x - x * x / 12.0 - 0.5 * f * x * x + f * f * x * x;
  } else {
    ds_dx = (f / (1.0 + f * x) - s / (1.0 + x)) / L;
    ds_df = x / ((1.0 + f * x) * L);
  }
  // x = Pt1/P1 - 1.
  const double ds_dP1 = ds_dx * (-Pt1 / (P1 * P1));
  const double ds_dPt1 = ds_dx / P1;
  const double Tf = (s == 1.0) ? Tt1 : T1 * std::exp(s * M);
  const double dTf_dT1 = Tf * (1.0 - s) / T1;
  const double dTf_dTt1 = Tf * s / Tt1;
  const double dTf_dP1 = Tf * M * ds_dP1;
  const double dTf_dPt1 = Tf * M * ds_dPt1;
  const double dTf_df = Tf * M * ds_df;
  // Each side's ideal flow from rest at its own (pressure, temperature):
  // beta 0, so the kernel expands from (P_total_up, Tt_up) as given.
  const solver::OrificeResult sr =
      solver::orifice_compressible_residuals_and_jacobian(0.0, Pf, Pf, Tf, Y, P2, 1.0, 1.0, 0.0);
  const solver::OrificeResult t =
      solver::orifice_compressible_residuals_and_jacobian(0.0, Pt1, Pt1, Tt1, Y, P2, 1.0, 1.0, 0.0);
  if (!(t.m_dot_calc > 0.0)) {
    return out;
  }
  const double r = sr.m_dot_calc / t.m_dot_calc;
  const double inv = 1.0 / t.m_dot_calc;
  const double gP = sr.d_mdot_dP_total_up;
  const double gT = sr.d_mdot_dTt_up;
  out.ratio = r;
  out.d_dP1 = (gP * dPf_dP1 + gT * dTf_dP1) * inv;
  out.d_dT1 = gT * dTf_dT1 * inv;
  out.d_dPt1 = (gP * dPf_dPt1 + gT * dTf_dPt1) * inv - r * t.d_mdot_dP_total_up * inv;
  out.d_dTt1 = gT * dTf_dTt1 * inv - r * t.d_mdot_dTt_up * inv;
  out.d_dP2 = (sr.d_mdot_dP_static_down - r * t.d_mdot_dP_static_down) * inv;
  out.d_df = (gP * dPf_df + gT * dTf_df) * inv;
  for (std::size_t k = 0; k < Y.size() && k < sr.d_mdot_dY_up.size(); ++k) {
    out.d_dY[k] = (sr.d_mdot_dY_up[k] - r * t.d_mdot_dY_up[k]) * inv;
  }
  return out;
}

HoleDriveCdRatio hole_static_to_total_cd_ratio(double P1, double T1, double Pt1, double Tt1,
                                               double P2, const std::vector<double>& Y,
                                               bool incompressible) {
  return hole_supply_recovery_cd_ratio(P1, T1, Pt1, Tt1, P2, Y, 0.0, incompressible);
}

}  // namespace combaero
