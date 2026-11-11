#include "hole_cd.h"

#include <limits>
#include <memory>
#include <tuple>

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

}  // namespace combaero
