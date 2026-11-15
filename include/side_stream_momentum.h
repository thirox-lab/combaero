#pragma once

// Momentum of side streams joining a channel, compressible (#465, #471).
//
// A stream that enters a channel brings mass, and only the streamwise part of
// its momentum. The channel flow must make up the rest, which costs static
// pressure. Over one injection station of constant cross-section A, flow m_a
// arriving and m_b leaving, the impulse balance is exact for compressible
// flow:
//
//   P_a A + m_a u_a + J = P_b A + m_b u_b,    u = m / (rho A)
//
// with J the side streams' streamwise momentum (zero for normal injection).
// Density is each station's own (static P, T); nothing here assumes it
// constant. At low Mach it reduces to Florschuetz, Truman and Metzger's
// (1981) P + G^2/rho = const.

#include <vector>

namespace combaero {

// ---------------------------------------------------------------
// One station, half at a time: merge AND bleed (#471)
// ---------------------------------------------------------------
//
// A station is a node where a side stream joins (merge) or leaves (bleed) a
// channel of cross-section A; m_a arrives, m_b leaves along the channel.
// With the side stream's AXIAL velocity written as u_s = kappa * u_a,
// momentum over the station gives the static-pressure DROP
//
//   dP = (m_b^2 - m_a^2 + kappa m_a (m_a - m_b)) / (rho A^2)
//
// at the station node's own density. The momentum flux m u = m^2/(rho A) is
// positive whichever way the gas moves, so the form holds for reversed flow
// and for a station the flow leaves both ways (#815). An m|m| form, used
// until #815, flipped the sign of the drop for reversed flow: a merge in a
// reversed chain RAISED the static pressure along the flow, and a split
// station had no physical root.
//
//   kappa = 0     normal injection (merge): the jets bring no axial momentum.
//                 Florschuetz's P + G^2/rho = const; the #465 impingement term.
//   kappa = 0.75  bleed through the wall: Bassett, Winterbone & Pearson (2001)
//                 separating straight-run coefficient K2/K5, Eq. (15),
//                 K = q^2 - 1.5 q + 0.5 (q = m_b/m_a), theta- and
//                 psi-independent, rated excellent against their data. The
//                 station reproduces it for EVERY q exactly at kappa = 0.75:
//                 (Pt_a - Pt_b)/(rho u_a^2 / 2) = 2(1-q)(kappa - (1+q)/2).
//
// A network segment between two stations carries HALF of each, so every node
// sits at the middle of its station (the centred scheme of #465). The
// segment's own flow then enters both halves at their own densities, which
// keeps the term m^2 (1/rho_from - 1/rho_to)/(2 A^2) that a single-density
// form drops. One state per node makes the scheme consistent to
// O(M^2 dm/m), not exact.
constexpr double STATION_KAPPA_MERGE_NORMAL = 0.0;
constexpr double STATION_KAPPA_BLEED_BASSETT = 0.75;

struct StationHalfDrop {
  double dP = 0.0;      // half the station's static-pressure drop [Pa]
  double d_dm_a = 0.0;  // [Pa/(kg/s)]
  double d_dm_b = 0.0;  // [Pa/(kg/s)]
  double d_dP = 0.0;    // d/d(station static P) [-]
  double d_dT = 0.0;    // d/d(station T) [Pa/K]
  // d/d(mass fractions Y_k) of the gas at fixed (P, T) [Pa per unit mass
  // fraction], one entry per species (#1146). The drop is ~ 1/rho, and
  // rho = P MW / (R_u T), so d/dX_j = -dP W_j / MW, mapped onto the solver's
  // mass fractions by mole_gradient_to_mass. Zero (all entries) where the
  // drop is.
  std::vector<double> d_dY;
};

StationHalfDrop station_half_drop(double m_a, double m_b, double P, double T,
                                  const std::vector<double>& X, double area,
                                  double kappa);

// Entry from a reservoir (Pt) into a duct, static at the entry node:
//   dP = (1 + K_in) max(m, 0)^2 / (2 rho A^2)
// the dynamic head plus the entry loss (K_in = 0.5 sharp, ~0 well rounded).
// Reversed (m < 0) the duct discharges INTO the reservoir: static continuity,
// the dynamic head lost, the segment's own convention at an exit into a
// plenum. So 0, not the mirrored m|m| that charged the outflow (1 + K_in)
// dynamic heads of static drop (#705). C1 through m = 0.
// Area change over a duct segment between two static ends (#812), Bernoulli:
//   dP = P_from - P_to = m^2 / (2 rho) (1 / A_to^2 - 1 / A_from^2)
// even in m (lossless; friction is the segment's own), at the density of the
// gas the segment carries. Zero, with zero partials, for A_from == A_to.
// Derivative in d_dm_a.
StationHalfDrop segment_area_change_drop(double m, double P, double T,
                                         const std::vector<double>& X, double area_from,
                                         double area_to);

StationHalfDrop channel_entry_drop(double m, double P, double T,
                                   const std::vector<double>& X, double area,
                                   double K_in);

// ---------------------------------------------------------------
// Merge chamber: main inlet + side streams -> one outlet (#471)
// ---------------------------------------------------------------
//
// Constant-area control volume of area A, axis = the outlet direction. The
// chamber node carries the outlet state: static P and the STAGNATION T (#578);
// the outlet's static T comes from chamber_static_state, as in the chamber's
// own closure (#576). The main stream arrives along the axis at the MAIN FACE
// with stagnation Tt_main; the face's static (P_f, T_f) solve
//
//   P_f + m_main^2 R T_f / (P_f A^2) = P + m_out^2 / (rho A^2) - J / A
//   h(Tt_main) = h(T_f) + u_f^2 / 2,   u_f = m_main R T_f / (P_f A)
//
// by fixed point (the impulse quadratic on its SUBSONIC root, then the energy
// solve), with analytic sensitivities from the implicit function theorem.
// Pt_face = P_f exp((s(Tt_main) - s(T_f)) / R), the same entropy-based closure
// as the chamber itself (#357). No
// loss coefficient: the stagnation-pressure loss of mixing follows from
// momentum; transverse momentum is reacted by the walls. With no side streams
// (m_out = m_main, J = 0, Tt = Tt_main) the face IS the chamber.
//
// When the impulse cannot be carried subsonically (the main stream would
// have to choke at the face) the discriminant is held at a floor and
// `choked` is set: the result is then a bound, not a solution.
struct MergeFaceState {
  double P_face = 0.0;   // [Pa]
  double Pt_face = 0.0;  // [Pa]
  double M_face = 0.0;   // [-]
  bool choked = false;
  // d/d(m_main, m_out, J, P, Tt, Tt_main), all analytic; the dT/dT_main
  // field suffixes name the stagnation temperatures Tt and Tt_main (implicit function
  // theorem through the static-state solves, #578).
  double dPf_dm_main = 0.0, dPf_dm_out = 0.0, dPf_dJ = 0.0;
  double dPf_dP = 0.0, dPf_dT = 0.0, dPf_dT_main = 0.0;
  double dPtf_dm_main = 0.0, dPtf_dm_out = 0.0, dPtf_dJ = 0.0;
  double dPtf_dP = 0.0, dPtf_dT = 0.0, dPtf_dT_main = 0.0;
  // Gradients in the MASS fractions of the chamber (outlet) gas and of the
  // main stream (#743), one entry per species.
  std::vector<double> dPf_dY, dPf_dY_main, dPtf_dY, dPtf_dY_main;
};

//
// incompressible (#588): the low-Mach limit. Both densities are taken at the
// chamber pressure and the stagnation temperatures (rho = P / (R Tt)), so the
// impulse is linear in P_f and there is no root to choose and nothing to
// choke:
//   P_f  = P + m_out^2 / (rho A^2) - J / A - m_main^2 / (rho_main A^2)
//   Pt_f = P_f + m_main^2 / (2 rho_main A^2)
MergeFaceState chamber_merge_face_state(double m_main, double Tt_main,
                                        const std::vector<double>& X_main,
                                        double m_out, double P, double Tt,
                                        const std::vector<double>& X,
                                        double side_momentum, double area,
                                        bool incompressible = false);

// Streamwise impulse of a jet discharged from stagnation (Pt, Tt) into static
// P, per unit cos(theta): J = m w.
//   unchoked: w = u_is, the isentropic velocity to P (variable cp);
//   choked (P < P*): w = u* + (P* - P) / (rho* u*), the sonic jet's momentum
//                    plus its pressure thrust on the vena-contracta area.
// Cd does not enter: it sets the area, not the jet velocity. w is C1 across
// choking. P >= Pt gives w = 0 (no outflow to carry momentum). The partials
// are analytic on the variable-cp isentrope (#816; they were central
// differences, about 1% off once choked).
struct JetImpulse {
  double J = 0.0;       // [N]
  double dJ_dm = 0.0;   // = w [m/s]
  double dJ_dPt = 0.0, dJ_dTt = 0.0, dJ_dP = 0.0;
  // d/d(mass fractions Y_k) of the jet's gas at fixed (Pt, Tt, P) [N per
  // unit mass fraction], one per species, both branches analytic (#1146).
  std::vector<double> dJ_dY;
  bool choked = false;
};

JetImpulse jet_impulse(double m, double Pt, double Tt, double P,
                       const std::vector<double>& X);

// ---------------------------------------------------------------
// Supply-side crossflow ratio of a duct-fed hole (#471)
// ---------------------------------------------------------------
//
// McGreehan & Schotsch (1988) correlate a hole fed from a duct on U1/Vi:
//   U1 : the duct's tangential (approach) velocity at the hole,
//        (m_a + m_b) / (2 rho A) -- the mean of the flows arriving at and
//        leaving the hole's station, the same mean the bleed station uses;
//   Vi : the IDEAL jet velocity built from the duct's STATIC state (their
//        caution, extraction item 22: "not the total pressure based on the
//        tangential velocity") -- the isentropic velocity from (P, T) to the
//        discharge static P_down, jet_impulse's w with the static state as
//        its stagnation.
// Vi is held above VI_FLOOR (sqrt(Vi^2 + floor^2)) so a vanishing drive gives
// a large, finite ratio rather than a division by zero; the correlation's own
// regularisation handles U1/Vi -> 0.
constexpr double CROSSFLOW_VI_FLOOR = 0.1;  // [m/s]

struct CrossflowRatio {
  double U1_over_Vi = 0.0;
  double U1 = 0.0, Vi = 0.0;  // [m/s]
  double d_dm_a = 0.0, d_dm_b = 0.0;   // [1/(kg/s)]
  double d_dP = 0.0, d_dT = 0.0;       // duct station static P [1/Pa], T [1/K]
  double d_dP_down = 0.0;              // discharge static [1/Pa]
};

CrossflowRatio crossflow_velocity_ratio(double m_a, double m_b, double P,
                                        double T, const std::vector<double>& X,
                                        double P_down, double area);

}  // namespace combaero
