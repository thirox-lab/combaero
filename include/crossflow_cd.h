#pragma once

// Discharge coefficient of a wall hole under a grazing crossflow (#800).
//
// A hole whose exit (or entrance) faces a fast stream passing along the wall
// discharges far less than the same hole into still air: the stream bends
// the jet over and blocks part of the hole. Both laws below are written in
//
//   q = G_ideal / G_c
//
// G_ideal = sqrt(2 rho dP) the ideal (Cd = 1) mass velocity through the hole
// from the supply to the crossflow's static pressure, G_c the crossflow's
// mass velocity. The sources state their laws in the ACTUAL jet velocity
// (Vi = Cd Vi_ideal, Gj = Cd G_ideal), which makes the hole's own flow an
// argument of its Cd: a flow-dependent Cd whose derivative at zero flow is
// infinite. Both are resolved here so that Cd is a function of q alone; the
// hole flow then goes as dP^1.07 under a strong crossflow (a linear
// resistance), smooth through dP = 0, instead of sqrt(dP).
//
// Each law reports where it holds the hole's own (crossflow-free) Cd:
// ``capped``. With no crossflow (G_c = 0) every law is capped.

namespace combaero {

// Rogers & Hersh (1976), NASA CR-2681, Fig. 11 (p. 41): single square-edged
// holes, D 12.7 mm, L/D 1, grazing V_inf 24-85 m/s, Vi/V_inf 0.004-0.9, Cd
// referenced to the grazing STATIC pressure:
//
//   injection (into the grazing stream):  Cd = 0.97 (Vi / V_inf)^0.534
//   suction (out of the grazing stream):  Cd = 0.80 (Vi / V_inf)^0.47
//
// The paper calls injection "outflow" and suction "inflow" (from the cavity's
// side). Resolved in q (equal densities, Vi / V_inf = Gj / G_c):
// Cd = a^(1/(1-b)) q^(b/(1-b)), capped at the hole's own Cd0. The cap is
// joined C1 over +/-ROGERS_HERSH_1976_CAP_BAND of Cd0 (a min() kink there
// stalled Newton on cold CR-3630 arrays): the fit itself below 0.9 Cd0,
// exactly Cd0 above 1.1 Cd0 (every Florschuetz 1981 row sits above 2), at
// most 2.5% under the fit in between.
constexpr double ROGERS_HERSH_1976_CAP_BAND = 0.1;
constexpr double ROGERS_HERSH_1976_INJECTION_A = 0.97;
constexpr double ROGERS_HERSH_1976_INJECTION_B = 0.534;
constexpr double ROGERS_HERSH_1976_SUCTION_A = 0.80;
constexpr double ROGERS_HERSH_1976_SUCTION_B = 0.47;

// Florschuetz, Metzger & Su (1982), NASA CR-3630, Table 4.1 and Eq. 4.2:
// the special Cd tests on the B-size jet plates under crossflow. In
// x = Gc/Gj, Cd is the plate's Table 3.1 value up to a break point, then
//
//   zeta(x) = 1.06 / (x + 0.806)^0.602            (Eq. 4.2)
//
//   B(5,4,1): 0.85 to 0.63, then zeta
//   B(5,4,2): 0.85 to 0.83, then 1.08 zeta
//   B(5,4,3): 0.85 to 0.90, then 1.11 zeta
//   B(5,8,1): 0.80 to 0.54, then -0.169 x + 0.893 to 1.5, then zeta
//   B(10,4,1), B(10,8,1): 0.76 to 0.54, then -0.128 x + 0.825 to 1.8, then
//             zeta (INFERRED by the authors from the (5,x) plates, not tested)
//
// The printed pieces meet only to their printed digits (B(5,8,1)'s line
// gives 0.802 at its 0.54 break against the constant 0.80, and lands 0.0015
// below zeta at 1.5). Resolved implicitly, a jump would leave a band of q
// with no Cd at all, so the law is made continuous within those digits: the
// line keeps its printed slope but passes through zeta at its end (intercept
// 0.8945 for 0.893, 0.8256 for 0.825), and the constant joins the falling
// part C1 over +/-FLORSCHUETZ_METZGER_SU_1982_CAP_BAND of Cd0 (as Rogers &
// Hersh's cap; a min() kink stalled Newton on B(5,4,1) at mc/mj 0.99). Gj = Cd G_ideal makes x depend on Cd:
// x = 1 / (Cd q), solved here (closed form on the line, Newton on zeta) with
// the implicit derivative. The Cd is the table's own, not the caller's: Cd0
// is not an input.
constexpr double FLORSCHUETZ_METZGER_SU_1982_CAP_BAND = 0.1;
constexpr double FLORSCHUETZ_METZGER_SU_1982_ZETA_C = 1.06;
constexpr double FLORSCHUETZ_METZGER_SU_1982_ZETA_X0 = 0.806;
constexpr double FLORSCHUETZ_METZGER_SU_1982_ZETA_N = 0.602;

struct CrossflowCd {
  double Cd = 0.0;
  double dCd_dG_ideal = 0.0;  // [1/(kg/(m^2 s))]
  double dCd_dG_c = 0.0;      // [1/(kg/(m^2 s))]
  bool capped = false;        // the crossflow-free Cd holds
};

// Rogers & Hersh (1976). ``suction``: the hole draws FROM the grazing stream
// (a reversed jet). G_ideal <= 0 gives Cd = 0 unless G_c <= 0 (capped).
CrossflowCd rogers_hersh_1976_cd(double Cd0, double G_ideal, double G_c, bool suction);

// True when Table 4.1 covers the plate (xn/d, yn/d, z/d), to 1e-9.
bool florschuetz_metzger_su_1982_has_plate(double xn_d, double yn_d, double z_d);

// Table 4.1 for a plate it covers; throws std::invalid_argument otherwise.
CrossflowCd florschuetz_metzger_su_1982_cd(double xn_d, double yn_d, double z_d, double G_ideal,
                                           double G_c);

// -------------------------------------------------------------------------
// EXIT (external) crossflow on a film-cooling hole (#1074, #386)
// -------------------------------------------------------------------------
//
// Rowbury, Oldfield & Lock (2001), ASME J. Turbomach. 123, 258-265: an
// additive OUTLET loss on top of the hole's no-crossflow Cd,
//
//   delta_out = 1/Cd^2 - 1/Cd,no^2                                   (Eq. 7)
//   delta_out = A k_perp^-B - C                                      (Eq. 8)
//   k_perp    = rho_h u_h u_h,perp / (rho_m u_m^2)
//
// rho_h u_h the hole's ACTUAL mass flux at its exit plane, u_h,perp = u_h
// sin(theta) the component normal to the surface (theta the hole axis to
// the surface; see "PERPENDICULAR" below), rho_m u_m^2 the mainstream's
// momentum flux. A, B, C are Table 2 (p. 263) per vane region and hole
// shape, functions of theta and the orientation (compound) angle alpha:
//
//   cylindrical, pressure surface (rows 1-6):
//     A = (2.95 sin th - 0.88)(1 + sin al)  B = (0.714 sin th + 0.393)(1 + sin al)  C = 0.10
//   cylindrical, leading edge (rows 7-10):
//     A = 0.079 sin th (1 + sin al)         B = (4.795 sin th - 4.487)(1 + sin al)  C = 0.28
//   cylindrical, suction surface (rows 11-14):
//     A = (0.682 sin th - 0.24)(1 + sin al) B = (0.777 sin th + 0.405)(1 + sin al)  C = 0.10
//   fan-shaped, pressure surface:
//     A = (0.435 sin th - 0.140)(1 + sin al) B = (0.350 sin th + 0.570)(1 + sin al) C = 0.02
//   fan-shaped, suction surface:
//     A = (0.225 sin th - 0.110)(1 + sin al) B = (0.350 sin th + 0.570)(1 + 2.8 sin al) C = 0.05
//   (no fan-shaped leading-edge set: "N/A")
//
// For a fan-shaped hole the authors take the momentum flux at the hole EXIT
// plane (p. 262), i.e. the metering mass flux over the exit/inlet area ratio.
//
// THE TWO ENDS OF Eq. 8 ARE NOT THE PAPER'S RANGE, and are handled here:
//
// * k -> infinity (no mainstream, u_m -> 0): delta -> -C, NOT 0. Eq. 8 has
//   crossover built in (the paper's Fig. 9 Ma_m 0.3 prediction rises above
//   its Ma_m 0 one at PR > 1.5), but its no-crossflow prediction is Cd,no by
//   definition (Fig. 9, "Predict Mm = 0", flat at Lichtarowicz's 0.776). So
//   the limit is discontinuous at u_m = 0 exactly; the caller applies the
//   term only to a mainstream with momentum (rho_m u_m^2 > 0) and returns
//   Cd,no otherwise. No upper k is printed (Fig. 4 spans k_out 0.2-3.3 for
//   row 1 only), so none is imposed.
//
// * k -> 0 (no hole flow, or a mainstream overwhelming it): A k^-B diverges
//   and Cd -> 0, so a solver iterating on the hole's own flow sees a
//   spurious zero-flow root; and with B > 1 (near-normal holes, compound
//   angles) Eq. 7-8 solved on the actual flow have NO root at all below some
//   k (the m = Cd(m) K fixed point turns over). The term is therefore HELD
//   at k_floor = max(ROWBURY_2001_K_PERP_MIN, k_unique) -- value kept, slope
//   zero, the pattern of the discharge correlations' Re floor:
//     ROWBURY_2001_K_PERP_MIN = 0.01: the lowest k_perp the authors evaluate
//       their own correlation at -- Fig. 9's Ma_m 1.2 prediction starts at
//       PR ~1.11, Cd ~0.425, which is k_perp ~0.010 (back-calculated, ours).
//     k_unique (B > 1 only): where (B - 1) A k^-B = M (1 - C), M =
//       ROWBURY_2001_UNIQUE_MARGIN. On the actual flow the fixed point is
//       unique while (B - 1) A k^-B < 1/Cd,no^2 - C; 1/Cd,no^2 >= 1 makes
//       (1 - C) a bound for any Cd,no, and the margin leaves the
//       self-consistent solve a slope d ln Cd / d ln m <= B / (2B - 1) < 1
//       from this term (room for the Re dependence of Cd,no). OUR choice,
//       not the paper's. The hold starts where A k^-B = M (1 - C) / (B - 1):
//       for B near 1 (a normal hole, suction set: B 1.18, k_floor 0.23) only
//       once the term has cut Cd by a third or more; for larger B sooner.
//
// PERPENDICULAR. Fig. 3 measures theta in the plane section and alpha in
// plan view; for alpha = 0 the surface-normal velocity is u_h sin(theta)
// exactly. With a compound angle the true normal component also depends on
// alpha; the authors fold alpha into A, B through (1 + sin alpha) instead,
// so u_h,perp = u_h sin(theta) is used throughout (theta read as the hole
// axis to the surface, the elements' angle_deg).
//
// FITTED ENVELOPE (Table 1, CHTT rows): flagged by rowbury_2001_in_range,
// never refused. Coefficients with A <= 0 or B <= 0 (theta below the set's
// zero crossing) give Eq. 8 no meaning and are refused.
enum class Rowbury2001Set {
  CylindricalPressureSurface,
  CylindricalLeadingEdge,
  CylindricalSuctionSurface,
  FanShapedPressureSurface,
  FanShapedSuctionSurface,
};

namespace rowbury2001 {
// Table 2: A = (a1 sin th + a0)(1 + sa sin al), B = (b1 sin th + b0)(1 + sb sin al), C.
struct Row {
  double a1, a0, sa;
  double b1, b0, sb;
  double c;
};
constexpr Row kCylPressure{2.95, -0.88, 1.0, 0.714, 0.393, 1.0, 0.10};
constexpr Row kCylLeadingEdge{0.079, 0.0, 1.0, 4.795, -4.487, 1.0, 0.28};
constexpr Row kCylSuction{0.682, -0.24, 1.0, 0.777, 0.405, 1.0, 0.10};
constexpr Row kFanPressure{0.435, -0.140, 1.0, 0.350, 0.570, 1.0, 0.02};
constexpr Row kFanSuction{0.225, -0.110, 1.0, 0.350, 0.570, 2.8, 0.05};

// Table 1 angle spans of the rows each set was fitted on [deg]. Fan-shaped
// rows (p. 259: all but rows 5-10 and 11 were also run fan-shaped): pressure
// rows 1-4, suction rows 12-14.
struct Envelope {
  double theta_min, theta_max, alpha_min, alpha_max;
};
constexpr Envelope kEnvCylPressure{20.0, 58.0, 0.0, 7.0};    // rows 1-6
constexpr Envelope kEnvCylLeadingEdge{83.0, 90.0, 60.0, 60.0};  // rows 7-10
constexpr Envelope kEnvCylSuction{30.0, 78.0, 0.0, 6.0};     // rows 11-14
constexpr Envelope kEnvFanPressure{20.0, 30.0, 0.0, 7.0};    // rows 1-4
constexpr Envelope kEnvFanSuction{30.0, 40.0, 0.0, 6.0};     // rows 12-14
}  // namespace rowbury2001

constexpr double ROWBURY_2001_K_PERP_MIN = 0.01;
constexpr double ROWBURY_2001_UNIQUE_MARGIN = 0.5;

struct Rowbury2001Coefficients {
  double A = 0.0;
  double B = 0.0;
  double C = 0.0;
};

// Table 2 at (theta, alpha) [deg]; theta in (0, 90], alpha in [0, 90].
// Throws std::invalid_argument outside those, or where A <= 0 or B <= 0.
Rowbury2001Coefficients rowbury_2001_coefficients(Rowbury2001Set set, double theta_deg,
                                                  double alpha_deg);

// The k_perp below which the term is held (see above).
double rowbury_2001_k_floor(Rowbury2001Set set, double theta_deg, double alpha_deg);

// Inside the set's Table 1 angle span (to 1e-9 deg).
bool rowbury_2001_in_range(Rowbury2001Set set, double theta_deg, double alpha_deg);

struct Rowbury2001Cd {
  double Cd = 0.0;
  double dCd_dCd0 = 0.0;
  double dCd_dk = 0.0;   // 0 where held
  double delta = 0.0;    // delta_out at the k used
  double k_used = 0.0;   // max(k_perp, k_floor)
  bool held = false;     // k_perp < k_floor
};

// Cd with the exit-crossflow loss, from the no-crossflow Cd0 (0 < Cd0 <= 1)
// and k_perp >= 0. Throws where 1/Cd0^2 + delta <= 0 (cannot happen for
// Cd0 <= 1, C < 1).
Rowbury2001Cd rowbury_2001_exit_crossflow_cd(double Cd0, double k_perp, Rowbury2001Set set,
                                             double theta_deg, double alpha_deg);

struct Rowbury2001KPerp {
  double k = 0.0;
  double dk_dG = 0.0;      // [1/(kg/(m^2 s))]
  double dk_drho_h = 0.0;  // [m^3/kg]
  double dk_dq_m = 0.0;    // [1/Pa]
};

// k_perp = G_h^2 sin(theta) / (rho_h rho_m u_m^2), G_h >= 0 the hole's
// actual exit-plane mass flux, rho_h > 0, rho_m u_m^2 > 0 (the caller turns
// the term off at no mainstream momentum).
Rowbury2001KPerp rowbury_2001_k_perp(double G_h, double rho_h, double rho_m_u_m2,
                                     double theta_deg);

}  // namespace combaero
