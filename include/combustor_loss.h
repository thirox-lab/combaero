#pragma once

// Combustor total-pressure loss correlations, as a FRACTION of the inlet
// total pressure (Pt_out = Pt_in (1 - xi)), with every partial derivative.
//
// The Python correlation classes (combaero.network.pressure_loss) call these,
// and PressureLossElement builds its Jacobian from the partials instead of
// differencing the correlation (#718). A user-supplied Python callable keeps
// the differenced path; that interface is the prototype for loss maps.
//
//   fraction family:  xi = (k theta + xi0) g(|mdot| / m_ref)
//   head family:      xi = (k theta + zeta0) q,  q = 0.5 rho v^2 / Pt_in
//                     v = mdot / (rho area),  rho = Pt_in MW / (R Tt_in)
//                     =>  q = mdot^2 R Tt_in / (2 area^2 Pt_in^2 MW)
//
// The head family is Lefebvre's pressure-loss factor, TOTAL-TOTAL at low
// Mach: dPt/Pt3 = PLF (R/2) (mdot sqrt(Tt3) / (A Pt3))^2, with the
// compressor-delivery (inlet) stagnation state (#782 audit). At a plenum
// inlet Tt = T and Pt = P; a coupled Fanno-Rayleigh form that reduces to it
// at low Mach is #830.
//
// theta = T_burned / T_unburned - 1. g(r) = 1 - (1 - r^2)^3 below r = 1 and 1
// above (C2 at r = 1): a fraction loss quoted at its design flow m_ref fades
// with the dynamic head and vanishes at zero flow (#550). m_ref <= 0 means no
// fade (a constant fraction at every flow).
//
// The composition partial is with respect to the inlet MASS fractions Y_k,
// each taken with the others held fixed (the "{node}.Y[k]" columns the
// network relays): MW = sum(Y) / sum(Y_k / W_k).

#include <cstddef>
#include <vector>

namespace combaero {

struct CombustorLossResult {
  double xi = 0.0;         // fractional total-pressure loss [-]
  double dxi_dtheta = 0.0; // [-]
  double dxi_dmdot = 0.0;  // [1/(kg/s)]
  double dxi_dTt_in = 0.0; // [1/K]
  double dxi_dPt_in = 0.0; // [1/Pa]
  std::vector<double> dxi_dY; // per species [-]; empty for the fraction family
};

// xi = (k theta + xi0) g(|mdot| / m_ref).
CombustorLossResult combustor_fraction_loss(double theta, double mdot, double k,
                                            double xi0, double m_ref);

// xi = (k theta + zeta0) mdot^2 R Tt_in / (2 area^2 Pt_in^2 MW(Y_in)). Zero
// (with zero partials) for mdot <= 0, area <= 0 or Pt_in <= 0.
CombustorLossResult combustor_head_loss(double theta, double mdot, double Tt_in,
                                        double Pt_in,
                                        const std::vector<double> &Y_in,
                                        double k, double zeta0, double area);

// ---------------------------------------------------------------------------
// The composition point of PressureLossElement's loss models (#725): one call
// returns xi and every partial for any built-in model. The element (python/
// combaero/network/components.py) forms the residual Pt_dn - Pt_up (1 - xi),
// gathers the inputs from the network state (theta from the theta source's
// burned/unburned temperatures, the flow-upstream node's Tt, Pt and Y) and
// maps the partials onto unknown names; that needs the topology.
//
// The model's own inputs are a GENERIC vector x, so a model that reads more
// than theta (a loss map over theta and fuel splits, #725) is another kind
// here, not a new field on every caller. Index PRESSURE_LOSS_INPUT_THETA is
// theta for every kind; dxi_dx has the size of x.
//
// Kinds, each the kernel's above (no physics is derived here):
//   Fraction  combustor_fraction_loss: xi = (k theta + base) g(mdot / m_ref)
//   Head      combustor_head_loss: xi = (k theta + base) q at (Tt_in, Pt_in)
// A model with k = 0 does not read theta (x may then be empty), and a
// Fraction without a design flow (m_ref <= 0) does not read mdot.
//
// COMPOSITION. The value's density reads X_in, the mole fractions of the
// state (the network's clamped and normalised composition), and so do the
// partials by x; the partials by the state (mdot, Tt_in, Pt_in, Y_in) are
// taken at the raw mass fractions Y_in, the "{node}.Y[k]" unknowns. MW is
// invariant to the scale of Y, so these are the derivatives of the value
// through that normalisation. With X_in empty everything reads Y_in.

// Index of theta in PressureLossInputs::x, for every kind.
constexpr std::size_t PRESSURE_LOSS_INPUT_THETA = 0;

enum class PressureLossKind { Fraction, Head };

// The resolved model (Python's validation already applied).
struct PressureLossModel {
  PressureLossKind kind = PressureLossKind::Fraction;
  double k = 0.0;     // d(coefficient)/d theta [-]
  double base = 0.0;  // Fraction: xi0 [-]; Head: zeta0 [-]
  double m_ref = 0.0; // Fraction: design flow [kg/s]; <= 0: no fade
  double area = 0.0;  // Head: velocity reference area [m^2]
};

// The state the model reads.
struct PressureLossInputs {
  double mdot = 0.0;         // |m_dot| [kg/s]
  double Tt_in = 0.0;        // inlet stagnation temperature [K]
  double Pt_in = 0.0;        // inlet stagnation pressure [Pa]
  std::vector<double> Y_in;  // inlet mass fractions; the partials are by these
  std::vector<double> X_in;  // mole fractions the value reads; empty: Y_in
  std::vector<double> x;     // the model's inputs; x[PRESSURE_LOSS_INPUT_THETA] = theta
};

struct PressureLossResult {
  double xi = 0.0;            // fractional total-pressure loss [-]
  double dxi_dmdot = 0.0;     // [1/(kg/s)]
  double dxi_dTt_in = 0.0;    // [1/K]
  double dxi_dPt_in = 0.0;    // [1/Pa]
  std::vector<double> dxi_dY; // per species [-]; empty for a Fraction
  std::vector<double> dxi_dx; // per model input; the size of x
};

// One evaluation of the model. Throws std::invalid_argument for a model with
// k != 0 and no theta in x, and what the kernels throw.
PressureLossResult pressure_loss_compose(const PressureLossModel &model,
                                         const PressureLossInputs &inputs);

} // namespace combaero
