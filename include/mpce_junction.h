#pragma once

// -----------------------------------------------------------------------------
// Whole-element (f, J) for the momentum-CV junction (MultiPortChamberElement).
//
// C++ port of python/combaero/network/mpce_element.py::MultiPortChamberElement.residuals
// -- the Mynard physics path only. See that module for the derivation and
// MPCE_CPP_PORT_DESIGN.md for the port's sequence and gates.
//
// SCOPE. This kernel owns the physics residual and its Jacobian. It does NOT
// own the element's guards -- the all-ports-zero fallback, the
// wrong-direction soft barrier, or the degenerate-mask snapping. Those are
// solver policy rather than junction physics, they changed twice recently
// (#302, #303), and the barrier's residual is a continuity row plus a penalty
// whose Jacobian is three lines. They stay in Python, where they are cheap to
// iterate on; the caller reaches this kernel only once a state is a junction.
//
// WHAT THE SEEDING BUYS. The Python assembles its Jacobian from an explicit
// dKQ/dmdot block plus a hand-derived dR/dP column for the COMMON port alone.
// But K is a function of every port's velocity U_j = -mdot_j / (rho_j A_j),
// and rho_j moves with P_j, so every port's static pressure belongs in the
// Jacobian. Measured against central differences at a converged state, the
// two non-common P columns were off by 4.4e-3 and 2.3e-3. Seeding the whole
// element produces them with no derivation to get wrong, so this Jacobian is
// MORE complete than the Python's rather than equal to it. The gate is
// therefore finite differences of the residual, not agreement with the
// Python's analytic Jacobian.
//
// DENSITY. The kernel is thermodynamics-free: it takes rho and drho/dP per
// port as values. A Jacobian is first order, so seeding rho as
// rho_val + (drho/dP)(P - P_val) is exact, and the equation of state stays at
// the call site where the mixture is known.
//
// SEED LAYOUT (10). Chosen so the Jacobian comes out already expressed in the
// solver's own unknowns, with no chain rule left for the shim to apply:
//
//     0..2   P[i]            static pressure at port i
//     3..5   Pt[i]           total pressure at port i
//     6..8   outer_mdot[i]   the CONNECTING element's mass flow, not the
//                            junction-convention one -- port_signs maps them
//     9      Pt_jct          the junction's own unknown
// -----------------------------------------------------------------------------

#include <array>
#include <cmath>

#include "dual_number.h"
#include "math_constants.h"
#include "mynard_junction.h"

namespace combaero::solver {

constexpr int kMpcePorts = kMynardPorts;      // 3
constexpr int kMpceSeeds = 3 * kMpcePorts + 1; // 10
constexpr int kMpceRows = kMpcePorts + 1;      // N port rows + mass

struct MpceResidualJacobian {
  std::array<double, kMpceRows> residual{};
  // [row][seed], in the layout documented above.
  std::array<std::array<double, kMpceSeeds>, kMpceRows> jacobian{};
  // False when the flow pattern is not a junction (no supplier or no
  // collector). The caller owns that case -- see the SCOPE note.
  bool valid = false;
  // Reported so the caller can label ports without re-deriving the split.
  int common_port = -1;
  double k_term_sign = 0.0;
  std::array<double, kMpcePorts> k_per_port{};
};

struct MpceGeometry {
  std::array<double, kMpcePorts> area{};
  std::array<double, kMpcePorts> theta_rad{};
  // +1 when the connecting element's positive direction is OUT of the
  // junction, -1 when it is in. port_mdot[i] = port_sign[i] * outer_mdot[i].
  std::array<double, kMpcePorts> port_sign{};
  double joining_etransfer_alpha = 0.0;
  double eta_scale = 0.0;
  // Ratio of specific heats per port, for the compressible reference head.
  // Follows the same rule as rho and drho_dp above: the equation of state is
  // evaluated at the call site where the mixture is known, and enters here as
  // a value. gamma is a weak function of T, so carrying it as a constant is
  // the same order of approximation as the ideal drho_dp = rho/P already used.
  //
  // Left at zero, the head falls back to the incompressible 0.5*rho*u^2, so a
  // caller that has not been updated keeps its previous behaviour.
  std::array<double, kMpcePorts> gamma{};
};

inline MpceResidualJacobian mpce_residuals_and_jacobian(
    const std::array<double, kMpcePorts>& p_static,
    const std::array<double, kMpcePorts>& p_total,
    const std::array<double, kMpcePorts>& rho,
    const std::array<double, kMpcePorts>& drho_dp,
    const std::array<double, kMpcePorts>& outer_mdot, double pt_jct,
    const MpceGeometry& geom) {
  using D = DualN<kMpceSeeds>;
  MpceResidualJacobian out;

  std::array<D, kMpcePorts> P{};
  std::array<D, kMpcePorts> Pt{};
  std::array<D, kMpcePorts> Rho{};
  std::array<D, kMpcePorts> Mdot{}; // junction convention: positive = OUT
  for (int i = 0; i < kMpcePorts; ++i) {
    P[i] = D::seed(p_static[i], i);
    Pt[i] = D::seed(p_total[i], kMpcePorts + i);
    // rho is a function of P; first order is all a Jacobian carries.
    Rho[i] = D::constant(rho[i]);
    Rho[i].d[i] = drho_dp[i];
    Mdot[i] = D::seed(outer_mdot[i], 2 * kMpcePorts + i) * geom.port_sign[i];
  }
  D Pt_jct = D::seed(pt_jct, 3 * kMpcePorts);

  // Mynard's convention is positive INTO the junction, the opposite of the
  // element's port convention.
  std::array<D, kMpcePorts> U{};
  for (int i = 0; i < kMpcePorts; ++i) {
    U[i] = (0.0 - Mdot[i]) / (Rho[i] * geom.area[i]);
  }

  // The element re-points the AXIAL-BACK port -- the single one whose flow
  // opposes the other N-1 -- along the main duct at pi, whatever the caller
  // declared. Mynard's vessel-direction convention needs it; the declared
  // angles are measured from the main axis instead. Missing this was worth
  // 4293 Pa on the first golden case, so it is not a detail.
  //
  // BRANCH-ON-PRIMAL, and note the 1e-9 dead band: it matches the Python's
  // counts exactly, which is stricter than the closure's own `Q < 0` test.
  // A port inside the band belongs to neither count, so neither arm fires and
  // every angle stays as declared -- the same as the Python.
  std::array<double, kMpcePorts> theta = geom.theta_rad;
  {
    int n_pos = 0;
    int n_neg = 0;
    int arg_max = 0;
    int arg_min = 0;
    for (int i = 0; i < kMpcePorts; ++i) {
      if (U[i].v > 1.0e-9) ++n_pos;
      if (U[i].v < -1.0e-9) ++n_neg;
      if (U[i].v > U[arg_max].v) arg_max = i;
      if (U[i].v < U[arg_min].v) arg_min = i;
    }
    if (n_pos == 1 && n_neg == kMpcePorts - 1) {
      theta[arg_max] = M_PI;
    } else if (n_neg == 1 && n_pos == kMpcePorts - 1) {
      theta[arg_min] = M_PI;
    }
  }

  auto mynard = mynard_junction_loss_coefficient<kMpceSeeds>(
      U, geom.area, theta, geom.joining_etransfer_alpha, geom.eta_scale);
  if (!mynard.valid) return out;

  // Which port normalises, and which way the loss term signs.
  //   separating (one supplier): the common port is the supplier and holds
  //     the HIGHER Pt, so collectors sit below it  -> +1
  //   joining (one collector): the common port is the collector and holds
  //     the LOWER Pt, so suppliers sit above it    -> -1
  const bool separating = (mynard.n_supplier == 1);
  const double k_term_sign = separating ? 1.0 : -1.0;
  int common = -1;
  for (int i = 0; i < kMpcePorts; ++i) {
    const bool is_collector = U[i].v < 0.0;
    if (is_collector != separating) {
      common = i;
      break;
    }
  }

  std::array<D, kMpcePorts> K{};
  for (int j = 0; j < kMpcePorts - 1; ++j) K[mynard.k_port[j]] = mynard.K[j];

  // Common-side reference head.
  //
  // K is tabulated as dPt / q_dyn, so the reference has to be the same notion
  // of dynamic head the port nodes use. While MomentumChamberNode closed
  // incompressibly, 0.5*rho*u^2 was that notion and the two agreed. Once the
  // port node closes on the isentropic stagnation state, the incompressible
  // head drifts from it -- 11% at M = 1.0, 34% at M = 1.4 -- and the two
  // halves of the junction disagree about the same quantity (issue #357).
  //
  // So take the isentropic stagnation rise of the common port's OWN state:
  //   q = P * [(1 + (g-1)/2 * M^2)^(g/(g-1)) - 1]
  // which reduces to 0.5*rho*u^2 as M -> 0, i.e. exactly where the K data was
  // measured, and tracks the real rise above it.
  //
  // Computed from the port's own P, rho and u rather than from Pt - P: the
  // latter is the same number at the solution but couples every port row to
  // the COMMON port's Pt, which breaks the one-Pt-per-row structure the
  // assembled Jacobian relies on (and measured 7x slower on the bc_swap
  // network). Here the dependency set is unchanged from the incompressible
  // head -- P, rho, u of the common port -- so the sparsity is preserved.
  //
  // For an ideal gas a^2 = g*P/rho, hence M^2 = rho*u^2/(g*P); no temperature
  // is needed, which is why gamma alone is threaded through MpceGeometry.
  // u_com^2 rather than |u_com|^2 so no absolute value enters the derivative.
  D u_com = U[common];
  const double g_com = geom.gamma[common];
  D q_dyn_com = Rho[common] * u_com * u_com * 0.5;
  if (g_com > 1.0) {
    D M2 = (Rho[common] * u_com * u_com) / (P[common] * g_com);
    // As P * expm1(g/(g-1) * log1p((g-1)/2 M^2)), not pow(1 + x, e) - 1:
    // the latter cancels to round-off as M -> 0 and returned q = 0 exactly
    // below M^2 ~ 1e-16 (a port near zero flow), while its dual still carried
    // the analytic slope -- with K ~ 1/u_com^2 that left J at 1.7e11 against
    // an FD of 0 (#1252).
    D x = M2 * (0.5 * (g_com - 1.0));
    q_dyn_com = P[common] * dexpm1(dlog1p(x) * (g_com / (g_com - 1.0)));
  }

  for (int i = 0; i < kMpcePorts; ++i) {
    D row = Pt[i] - Pt_jct + K[i] * q_dyn_com * k_term_sign;
    out.residual[i] = row.v;
    out.jacobian[i] = row.d;
    out.k_per_port[i] = K[i].v;
  }
  D mass = Mdot[0] + Mdot[1] + Mdot[2];
  out.residual[kMpcePorts] = mass.v;
  out.jacobian[kMpcePorts] = mass.d;

  out.common_port = common;
  out.k_term_sign = k_term_sign;
  out.valid = true;
  return out;
}

} // namespace combaero::solver
