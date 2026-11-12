#include "../include/combustor_loss.h"

#include "../include/composition.h"
#include "../include/input_validation.h"
#include "../include/thermo.h"
#include "../include/thermo_transport_data.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace combaero {

CombustorLossResult combustor_fraction_loss(double theta, double mdot, double k,
                                            double xi0, double m_ref) {
  // m_ref <= 0 switches the low-flow ramp off; NaN used to as well, and a
  // NaN mdot skipped the ramp (#769).
  if (!validation::is_finite(m_ref) || !validation::is_finite(mdot)) {
    throw std::invalid_argument("combustor_fraction_loss: m_ref and mdot must be finite");
  }
  CombustorLossResult r;
  const double base = k * theta + xi0;
  double g = 1.0;
  double dg_dm = 0.0;
  if (m_ref > 0.0) {
    const double rr = std::abs(mdot) / m_ref;
    if (rr < 1.0) {
      const double s = 1.0 - rr * rr;
      g = 1.0 - s * s * s;
      const double sign = mdot >= 0.0 ? 1.0 : -1.0;
      dg_dm = 6.0 * rr * s * s * sign / m_ref;
    }
  }
  r.xi = base * g;
  r.dxi_dtheta = k * g;
  r.dxi_dmdot = base * dg_dm;
  return r;
}

CombustorLossResult combustor_head_loss(double theta, double mdot, double Tt_in,
                                        double Pt_in,
                                        const std::vector<double> &Y_in,
                                        double k, double zeta0, double area) {
  if (Y_in.size() != molar_masses.size()) {
    throw std::invalid_argument("combustor_head_loss: Y_in size mismatch");
  }
  CombustorLossResult r;
  r.dxi_dY.assign(Y_in.size(), 0.0);
  if (mdot <= 0.0 || area <= 0.0 || Pt_in <= 0.0 || Tt_in <= 0.0) {
    return r;
  }
  // MW = S / U with S = sum(Y), U = sum(Y_k / W_k): mass_to_mole normalises,
  // so this is the MW the density sees [g/mol].
  double S = 0.0;
  double U = 0.0;
  for (std::size_t i = 0; i < Y_in.size(); ++i) {
    S += Y_in[i];
    U += Y_in[i] / molar_masses[i];
  }
  if (S <= 0.0 || U <= 0.0) {
    return r;
  }
  const double MW = S / U;
  // R in J/(mol K) over MW in g/mol: the specific gas constant is R*1000/MW.
  const double q = mdot * mdot * thermo::R_GAS * 1000.0 * Tt_in /
                   (2.0 * area * area * Pt_in * Pt_in * MW);
  const double c = k * theta + zeta0;
  r.xi = c * q;
  r.dxi_dtheta = k * q;
  r.dxi_dmdot = 2.0 * r.xi / mdot;
  r.dxi_dTt_in = r.xi / Tt_in;
  r.dxi_dPt_in = -2.0 * r.xi / Pt_in;
  // dMW/dY_k = 1/U - S/(U^2 W_k); dq/dMW = -q/MW.
  for (std::size_t i = 0; i < Y_in.size(); ++i) {
    const double dMW = 1.0 / U - S / (U * U * molar_masses[i]);
    r.dxi_dY[i] = -r.xi / MW * dMW;
  }
  return r;
}

PressureLossResult pressure_loss_compose(const PressureLossModel &model,
                                         const PressureLossInputs &inputs) {
  const std::vector<double> &x = inputs.x;
  double theta = 0.0;
  if (model.k != 0.0) {
    if (x.size() <= PRESSURE_LOSS_INPUT_THETA) {
      throw std::invalid_argument(
          "pressure_loss_compose: a model with k != 0 reads theta, x is empty");
    }
    theta = x[PRESSURE_LOSS_INPUT_THETA];
  }
  PressureLossResult out;
  out.dxi_dx.assign(x.size(), 0.0);
  // r: the partials by the state (mdot, Tt_in, Pt_in, Y_in); v: the value and
  // its partials by the model inputs x, at the composition the value reads.
  CombustorLossResult r;
  double dxi_dtheta = 0.0;
  if (model.kind == PressureLossKind::Fraction) {
    // A fraction without a design flow does not read the flow.
    const double m = model.m_ref > 0.0 ? inputs.mdot : 0.0;
    r = combustor_fraction_loss(theta, m, model.k, model.base, model.m_ref);
    out.xi = r.xi;
    dxi_dtheta = r.dxi_dtheta;
  } else {
    r = combustor_head_loss(theta, inputs.mdot, inputs.Tt_in, inputs.Pt_in, inputs.Y_in,
                            model.k, model.base, model.area);
    out.xi = r.xi;
    dxi_dtheta = r.dxi_dtheta;
    if (!inputs.X_in.empty()) {
      const CombustorLossResult v =
          combustor_head_loss(theta, inputs.mdot, inputs.Tt_in, inputs.Pt_in,
                              mole_to_mass(inputs.X_in), model.k, model.base, model.area);
      out.xi = v.xi;
      dxi_dtheta = v.dxi_dtheta;
    }
  }
  out.dxi_dmdot = r.dxi_dmdot;
  out.dxi_dTt_in = r.dxi_dTt_in;
  out.dxi_dPt_in = r.dxi_dPt_in;
  out.dxi_dY = std::move(r.dxi_dY);
  if (model.k != 0.0) {
    out.dxi_dx[PRESSURE_LOSS_INPUT_THETA] = dxi_dtheta;
  }
  return out;
}

} // namespace combaero
