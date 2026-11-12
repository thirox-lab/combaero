#pragma once

#include <vector>

#include "crossflow_cd.h"
#include "orifice.h"

// The discharge coefficient of a hole element, every term composed in one
// place (the hole-Cd family's composition point, #1095 follow-up). The
// network's hole elements (OrificeElement and its plate, film, slot and
// impingement subclasses) resolve their options in Python
// (python/combaero/network/hole_cd.py, _HoleCdChain) and call this once per
// evaluation; the mapping of the partials below onto the network's unknowns
// (the hole Reynolds number onto m_dot and the viscosity's T, P, Y; U1/Vi onto
// the backside duct's flows; the mainstream onto the discharge chamber's
// face) stays there, because it needs the network topology.
//
// The terms compose in this order, each the existing kernel's:
//
//   1. BASE Cd0, one of
//        Fixed      the caller's value, clamped into [HOLE_CD_FIXED_MIN, 1]
//                   (A_eff <= A_geom), no derivatives;
//        Metering   metering_cd_and_dRe (ISO 5167) at the PIPE Reynolds
//                   number Re_D;
//        Discharge  the discharge-hole correlation (orifice.h) at the HOLE
//                   Reynolds number Re and the INLET (supply-side) crossflow
//                   U1/Vi, the INCLINATION term (Idelchik diagram 3-2) inside
//                   the hole geometry (DischargeHoleGeometry::inclination_deg).
//                   Below discharge_cd_re_floor the value is held, so dCd/dRe
//                   is 0 there unless continue_re_floor (a Newton aid for a
//                   search ON Cd, see the floor's comment in orifice.h).
//   2. EXIT crossflow (opt-in, Rowbury et al. 2001), an additive outlet loss
//      on Cd0 at k_perp = rowbury_2001_k_perp(G_h, rho_h, rho_m u_m^2). Off
//      where the mainstream has no momentum (rho_m u_m^2 <= 0): Cd = Cd0.
//   3. JET crossflow law of an impingement row (Rogers & Hersh 1976, or
//      Florschuetz, Metzger & Su 1982), which REPLACES the hole's Cd where it
//      acts. Rogers & Hersh capped: the hole's Cd and its partials hold.
//      Where the law acts, the Cd's partials are the law's alone (dCd/dG_ideal,
//      dCd/dG_c); the hole's Cd enters Rogers & Hersh only as its cap. No
//      drive (G_ideal <= 0) or a capped table: no partials at all.
//
// No physics is derived here: every value and partial is a kernel's, chained.

namespace combaero {

// The floor of a Fixed Cd: a positive effective area.
constexpr double HOLE_CD_FIXED_MIN = 1e-4;

enum class HoleCdBase { Fixed, Metering, Discharge };

enum class HoleCdJetCrossflow { Off, RogersHersh1976, FlorschuetzMetzgerSu1982 };

// The resolved options (Python's "auto" logic and refusals already applied).
struct HoleCdModel {
  HoleCdBase base = HoleCdBase::Fixed;
  double fixed_cd = orifice::defaults::discharge_cd;  // Fixed only
  MeteringCdCorrelation metering = MeteringCdCorrelation::ReaderHarrisGallagher;
  DischargeCdCorrelation discharge = DischargeCdCorrelation::Idelchik1966Thick;
  bool continue_re_floor = false;  // Discharge: keep the slope below the floor

  bool exit_crossflow = false;  // Rowbury et al. (2001) term on
  Rowbury2001Set exit_set = Rowbury2001Set::CylindricalSuctionSurface;
  double exit_theta_deg = 90.0;  // hole inclination to the surface [deg]
  double exit_alpha_deg = 0.0;   // radial orientation [deg]

  HoleCdJetCrossflow jet = HoleCdJetCrossflow::Off;
  double xn_d = 0.0, yn_d = 0.0, z_d = 0.0;  // FlorschuetzMetzgerSu1982 plate
};

// The state the chain reads. Each term reads only its own inputs.
struct HoleCdInputs {
  double Re = 0.0;          // Discharge: hole Reynolds number [-]
  double U1_over_Vi = 0.0;  // Discharge: inlet crossflow ratio [-]
  double Re_D = 0.0;        // Metering: pipe Reynolds number [-]
  double G_h = 0.0;         // exit: the hole's actual exit-plane mass flux [kg/(m^2 s)]
  double rho_h = 0.0;       // exit: coolant exit-plane static density [kg/m^3]
  double rho_m_u_m2 = 0.0;  // exit: mainstream momentum flux [Pa]; <= 0: term off
  double G_ideal = 0.0;     // jet: ideal hole mass velocity [kg/(m^2 s)]
  double G_c = 0.0;         // jet: crossflow mass velocity [kg/(m^2 s)]
  bool suction = false;     // jet: the hole draws from the crossflow
};

struct HoleCdResult {
  double Cd = 0.0;  // the element's Cd, every term composed

  // Partials of Cd by input. Zero for an input the active terms do not read.
  double dCd_dRe = 0.0;
  double dCd_dU1_over_Vi = 0.0;
  double dCd_dRe_D = 0.0;
  double dCd_dG_h = 0.0;
  double dCd_drho_h = 0.0;
  double dCd_drho_m_u_m2 = 0.0;
  double dCd_dG_ideal = 0.0;
  double dCd_dG_c = 0.0;

  // Per-term values.
  double cd_base = 0.0;  // Cd0: term 1 (with the inclination term)
  double cd_hole = 0.0;  // terms 1-2: the hole's own Cd
  double dcd_base_dRe = 0.0;          // term 1 alone (floor applied)
  double dcd_base_dU1_over_Vi = 0.0;  // term 1 alone
  double dcd_base_dRe_D = 0.0;        // term 1 alone
  double inclination_dzeta = 0.0;     // discharge_inclination_dzeta, 0 at 90 deg
  bool re_held = false;               // Discharge: Re below the floor

  bool exit_active = false;  // the exit term evaluated (on, with mainstream)
  double exit_ratio = 0.0;   // cd_hole / cd_base (NaN for cd_base <= 0)
  Rowbury2001KPerp k_perp;   // exit_active only
  Rowbury2001Cd exit;        // exit_active only
  double exit_k_floor = 0.0;       // exit term on: rowbury_2001_k_floor
  double exit_k_tier = 0.0;        // 0 data, 0.5 authors' application, 1 faded
  bool exit_floor_above_fit = false;
  bool exit_in_range = false;      // Table 1 angle span

  bool jet_acts = false;      // the jet law's Cd replaced the hole's
  bool jet_partials = false;  // ... and carries partials (not capped, G_ideal > 0)
  bool jet_capped = false;    // the law's own capped flag
};

// One evaluation of the chain. `hole` is read by the Discharge base,
// `orifice` by the Metering base. Throws what the kernels throw.
HoleCdResult hole_cd_compose(const HoleCdModel& model, const DischargeHoleGeometry& hole,
                             const OrificeGeometry& orifice, const HoleCdInputs& inputs);

// The SUPPLY DRIVE of a duct-fed hole (#1087, supply side): the pressure the
// hole's ideal flow is taken from, the Cd's reference. A duct station holds
// the duct's STATIC pressure P with its gas moving at the mean flow m_dot
// through the duct section `area`.
//   Total   the duct's stagnation pressure at the station, P + q. The rule
//           (owner, 2026-11-12): a hole drives from the upstream total to
//           the downstream static, and its Cd model owns how much of the
//           supply head it recovers (the recovered fraction f, 0 to 1). A
//           static-referenced Cd is converted onto it exactly
//           (hole_static_to_total_cd_ratio).
//   Static  the station static P: the pre-0.12 drive, kept to reproduce
//           earlier results.
// q is the momentum-chamber closure's isentropic dynamic head of the station
// (solver::momentum_chamber_residual_and_jacobian: static T from h(Tt) =
// h(T) + v^2/2 at G = |m_dot| / area, compressible, or its low-Mach limit),
// so the station's total is the one a MomentumChamberNode would hold there;
// past the sonic static, that closure's bounded continuation (#788). It is
// even in m_dot and smooth through m_dot = 0 (a duct-flow reversal).
// Design record: docs/archive/hole_drive_convention.md.
enum class HoleSupplyDrive { Static, Total };

struct HoleSupplyDriveResult {
  double P_drive = 0.0;  // the drive pressure [Pa]
  double q_dyn = 0.0;    // the station's dynamic head Pt - P [Pa], either model
                         // (past the sonic static: the closure's bounded
                         // continuation, #788)
  // Partials of P_drive by the station's static P, stagnation T, the mean
  // duct flow (signed: P_drive is even in m_dot) and its mass fractions.
  double dP_dP = 0.0;
  double dP_dTt = 0.0;
  double dP_dm_dot = 0.0;
  std::vector<double> dP_dY;
};

HoleSupplyDriveResult hole_supply_drive(HoleSupplyDrive drive, double P, double Tt,
                                        const std::vector<double>& Y, double m_dot, double area,
                                        bool incompressible = false);

// The exact conversion of a STATIC-referenced hole Cd onto the TOTAL drive
// (owner's rule, 2026-11-12: a hole drives from the upstream total to the
// downstream static, and its Cd model owns how much supply dynamic head is
// recovered). One measured flow, two ideal flows:
//   Cd_T / Cd_S = m_id(P1, T1 -> P2) / m_id(Pt1, Tt1 -> P2),
// both from the element's own compressible nozzle kernel
// (solver::orifice_compressible_residuals_and_jacobian, Cd 1, beta 0), so
// Cd_T m_id(Pt1, Tt1 -> P2) is Cd_S m_id(P1, T1 -> P2) to round-off: the
// static-referenced model's own flow. (P1, T1) is the supply station's
// STATIC state (McGreehan & Schotsch 1988, p. 215: "static pressures and
// temperature"), (Pt1, Tt1) its stagnation state. Incompressibly the ratio
// is sqrt(1 - X), X = (Pt1 - P1)/(Pt1 - P2); the compressible ratio depends
// on the supply Mach as well (docs/archive/hole_drive_convention.md).
// `incompressible`: the exact ratio of the incompressible kernel (an element
// in the incompressible regime, the solver's incompressible proxy included),
// whose ideal flows are A sqrt(2 rho dp) at one density:
//   Cd_T / Cd_S = sqrt(1 - X),  X = (Pt1 - P1)/(Pt1 - P2),
// i.e. sqrt((P1 - P2)/(Pt1 - P2)); T1, Tt1 and Y do not enter (their
// partials are 0). The compressible ratio tends to it at low supply Mach.
// ratio = 0 (and every partial 0) where P1 <= P2: the static-referenced Cd
// gives no forward flow there, or where m_id(Pt1, Tt1 -> P2) = 0.
struct HoleDriveCdRatio {
  double ratio = 0.0;
  double d_dP1 = 0.0;   // [1/Pa]
  double d_dT1 = 0.0;   // [1/K]
  double d_dPt1 = 0.0;  // [1/Pa]
  double d_dTt1 = 0.0;  // [1/K]
  double d_dP2 = 0.0;   // [1/Pa]
  std::vector<double> d_dY;  // by mass fraction (the same Y on both sides)
};

HoleDriveCdRatio hole_static_to_total_cd_ratio(double P1, double T1, double Pt1, double Tt1,
                                               double P2, const std::vector<double>& Y,
                                               bool incompressible = false);

}  // namespace combaero
