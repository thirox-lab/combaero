#pragma once

#include <vector>

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
// rho_h u_h the hole's ACTUAL mass flux at its exit plane, rho_h the
// coolant's exit-plane density ("h = coolant at the hole exit plane", p.
// 264), u_h,perp the "component of coolant exit velocity perpendicular to
// the vane surface" (p. 264), rho_m u_m^2 the APPROACHING mainstream's
// momentum flux. A, B, C are Table 2 (p. 263) per vane region and hole
// shape, functions of theta and alpha:
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
// ANGLES (Fig. 3). theta, "conventional angle to surface (inclination), what
// you see on a plane section view": the hole axis's angle to the surface in
// the STREAMWISE-normal plane. alpha, "conventional definition of radial
// orientation, relative to the engine axis (or local wall surface)": drawn
// in the CROSS-STREAM plane, measured from the surface NORMAL (the vertical
// edge of the figure's end box). The plan-view angle (beta) and the "true"
// exit angle (zeta) are derived, not inputs. So the hole axis runs along
// (cot theta, 1, tan alpha) in (streamwise, normal, spanwise) and
//   u_h,perp / u_h = 1 / sqrt(cot^2 theta + 1 + tan^2 alpha)
// (sin theta at alpha = 0). The (1 + sin alpha) factors in Table 2 multiply
// A and B IN ADDITION: the k definition (p. 262 and the nomenclature) is the
// geometric normal component, and nothing in the text says alpha is folded
// into A, B instead. AMBIGUITY, stated: Fig. 3 draws alpha from the normal,
// but "relative to ... local wall surface" could be read as from the wall,
// which would give u_perp / u_h = 0.87 instead of 0.50 for the leading-edge
// rows (theta ~90, alpha 60): k_perp 1.7x larger there. The figure's reading
// is taken. Pressure- and suction-side rows have alpha <= 7 deg, where both
// readings agree to 0.8%.
//
// THE TWO ENDS OF Eq. 8 ARE NOT THE PAPER'S RANGE. Our extrapolations:
//
// * k > ROWBURY_2001_K_PERP_MAX = 3.3, the top of the one k range the paper
//   prints (Fig. 4, row 1: k_out 0.2-3.3; Figs. 5-9 plot Cd against PR, not
//   k). A large k is a jet that dominates the crossflow, where the outlet
//   loss must vanish (Cd -> Cd,no); Eq. 8 instead tends to -C, a fitted
//   offset, and would jump to Cd,no at u_m = 0 exactly (the authors' own
//   "Predict Mm = 0" is Cd,no by definition). Beyond k_max the OFFSET fades:
//     delta = A k^-B - C w(x),  x = ln(k / k_max),
//     w = 1 / (1 + (x / ROWBURY_2001_FADE_WIDTH)^2)
//   w(0) = 1 and w'(0) = 0, so delta and its slope are continuous at k_max;
//   delta -> 0 as k -> infinity (A k^-B decays too), so Cd -> Cd,no and
//   u_m -> 0 is continuous. Low slope: |d delta / d ln k| <= B A k^-B + 0.65
//   C / width. Not monotone where delta(k_max) < 0: the power law is still
//   falling at k_max and the faded offset must come back to 0, so delta dips
//   at most a few hundredths below delta(k_max) before rising; C1 and
//   monotone cannot both hold there, and C1 is what a Newton solve needs.
//   FADE_WIDTH = 1 (one e-fold of k halves the offset) is OUR choice.
//
// * k -> 0 (no hole flow, or a mainstream overwhelming it): A k^-B diverges
//   and Cd -> 0, so a solver iterating on the hole's own flow sees a
//   spurious zero-flow root; and with B > 1 (near-normal holes, compound
//   angles) Eqs. 7-8 solved on the actual flow have NO root at all below some
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
//       not the paper's, and results there depend on it: Hay 1983's 90 deg
//       holes (outside the suction rows' 30-78 deg) move with M. The hold
//       starts where A k^-B = M (1 - C) / (B - 1): for B near 1 (a normal
//       hole, suction set: B 1.18, k_floor 0.23) once the term has cut Cd by
//       a third or more; for larger B sooner. Off Table 1's angles k_floor
//       can exceed k_max (pressure set at theta 90, alpha 60: 2.9), i.e. the
//       whole fitted k range is held: rowbury_2001_k_floor reports it.
//
// Cd > 1. The -C crossover lets Cd exceed Cd,no, and with the leading-edge
// set (C 0.28) inside the fitted k range Cd exceeds 1 for Cd,no above about
// 0.9 (Cd,no 0.98 at theta 90, alpha 60, k 3.3: Cd 1.09). Not capped: the
// source's model allows it and a cap would put a kink in the hole's flow.
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
// rows (p. 259: rows 5-11 are cylindrical only): pressure rows 1-4, suction
// rows 12-14.
struct Envelope {
  double theta_min, theta_max, alpha_min, alpha_max;
};
constexpr Envelope kEnvCylPressure{20.0, 58.0, 0.0, 7.0};       // rows 1-6
constexpr Envelope kEnvCylLeadingEdge{83.0, 90.0, 60.0, 60.0};  // rows 7-10
constexpr Envelope kEnvCylSuction{30.0, 78.0, 0.0, 6.0};        // rows 11-14
constexpr Envelope kEnvFanPressure{20.0, 30.0, 0.0, 7.0};       // rows 1-4
constexpr Envelope kEnvFanSuction{30.0, 40.0, 0.0, 6.0};        // rows 12-14
}  // namespace rowbury2001

constexpr double ROWBURY_2001_K_PERP_MIN = 0.01;
constexpr double ROWBURY_2001_K_PERP_MAX = 3.3;
constexpr double ROWBURY_2001_FADE_WIDTH = 1.0;
constexpr double ROWBURY_2001_UNIQUE_MARGIN = 0.5;

struct Rowbury2001Coefficients {
  double A = 0.0;
  double B = 0.0;
  double C = 0.0;
};

// Table 2 at (theta, alpha) [deg]; theta in (0, 90], alpha in [0, 90).
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
  bool faded = false;    // k_used > ROWBURY_2001_K_PERP_MAX: the offset fades
};

// Cd with the exit-crossflow loss, from the no-crossflow Cd0 > 0 and
// k_perp >= 0 (+inf: Cd0). Throws where 1/Cd0^2 + delta <= 0.
Rowbury2001Cd rowbury_2001_exit_crossflow_cd(double Cd0, double k_perp, Rowbury2001Set set,
                                             double theta_deg, double alpha_deg);

struct Rowbury2001KPerp {
  double k = 0.0;
  double dk_dG = 0.0;      // [1/(kg/(m^2 s))]
  double dk_drho_h = 0.0;  // [m^3/kg]
  double dk_dq_m = 0.0;    // [1/Pa]
};

// k_perp = G_h^2 f / (rho_h rho_m u_m^2), f = u_h,perp / u_h =
// 1 / sqrt(cot^2 theta + 1 + tan^2 alpha); G_h >= 0 the hole's actual
// exit-plane mass flux, rho_h > 0, rho_m u_m^2 > 0 (the caller turns the
// term off at no mainstream momentum), alpha in [0, 90).
Rowbury2001KPerp rowbury_2001_k_perp(double G_h, double rho_h, double rho_m_u_m2,
                                     double theta_deg, double alpha_deg);

// The coolant's exit-plane static density for k_perp: the isentropic
// expansion from the supply (Pt, Tt) to the discharge static P, variable cp,
// s(T, P, X) = s(Tt, Pt, X). The fully expanded state also past choking
// (the sonic exit plane of a choked hole is not resolved). Partials
// analytic; drho_dY over mass fractions (the solver's unknowns).
struct Rowbury2001ExitDensity {
  double rho = 0.0;  // [kg/m^3]
  double T = 0.0;    // [K]
  double drho_dPt = 0.0, drho_dTt = 0.0, drho_dP = 0.0;
  std::vector<double> drho_dY;
};
Rowbury2001ExitDensity rowbury_2001_exit_density(double Pt, double Tt, double P,
                                                 const std::vector<double>& X);

}  // namespace combaero
