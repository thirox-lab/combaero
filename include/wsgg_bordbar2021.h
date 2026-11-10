#pragma once

// -----------------------------------------------------------------------------
// Pressure-dependent weighted-sum-of-gray-gases (WSGG) emissivity of H2O-CO2
// mixtures (#1071): Bordbar, Coelho, Fraga, Franca & Hostikka (2021),
// "Pressure-dependent weighted-sum-of-gray-gases models for heterogeneous
// CO2-H2O mixtures at sub- and super-atmospheric pressure", Int. J. Heat Mass
// Transfer 173, 121207 (CC BY 4.0). Fitted to HITEMP2010 line-by-line total
// emissivities, 300-3000 K, Mr = x_H2O/x_CO2 = 0.1-4, ten total-pressure sets
// 0.1-80 atm. Non-luminous: no soot.
//
// At one pressure set (Eqs. 7-11):
//   eps = sum_i a_i [1 - exp(-kappa_p,i (p_H2O + p_CO2) L)]       (7)
//   a_i = sum_j b_ij (T / 1500 K)^(j-1),  j = 1..5                (9)
//   b_ij = sum_k c_ijk Mr^(4-k),  kappa_p,i = sum_k d_ik Mr^(4-k)  (10, 11)
// four gray gases plus a transparent one, kappa_p in 1/(atm m), the partial
// pressures in atm. The constants are Appendix A's (wsgg_bordbar2021_data.h).
//
// Between the sets (the paper gives no rule; these choices are combaero's):
// - eps, not the coefficients, is interpolated: gray gas i is a different
//   band in each set, so mixing coefficients across sets mixes bands.
// - each bracketing set is evaluated at the ACTUAL partial-pressure path
//   (x_H2O + x_CO2) P L, so the interpolation carries only the pressure
//   broadening the sets differ by; the weight is linear in ln P (the sets are
//   spaced roughly geometrically). d eps/dP is continuous inside an interval
//   and kinks at the nodes.
// - outside 0.1-80 atm the nearest set is used (still at the actual path).
//
// Clamps (the fit's own range; the quartics diverge outside it):
// - T to [300, 3000] K (eps turns negative near 3500-4000 K); d/dT = 0 outside.
// - Mr to [0.1, 4]; x_CO2 = 0 (hydrogen or ammonia products) is Mr = 4 and
//   x_H2O = 0 is Mr = 0.1, both extrapolations of the model; d/dMr = 0 there.
// - eps to [0, 1]: sum_i a_i reaches 1.039 (2 atm), so the optically thick
//   limit can exceed 1.
// -----------------------------------------------------------------------------

#include <array>
#include <vector>

namespace combaero {

// Reference temperature of Eq. 9 [K].
constexpr double WSGG_BORDBAR2021_T_REF = 1500.0;
// Fitted range: temperature [K], mole ratio x_H2O/x_CO2 [-], pressure [atm].
constexpr double WSGG_BORDBAR2021_T_MIN = 300.0;
constexpr double WSGG_BORDBAR2021_T_MAX = 3000.0;
constexpr double WSGG_BORDBAR2021_MR_MIN = 0.1;
constexpr double WSGG_BORDBAR2021_MR_MAX = 4.0;
// Pa per standard atmosphere: the model's pressures are in atm.
constexpr double WSGG_PA_PER_ATM = 101325.0;

// The four gray gases of one pressure set at (T, Mr): weights a_i [-] and
// pressure absorption coefficients kappa_p,i [1/(atm m)] (Eqs. 9-11). For an
// RTE solve; the transparent gas has weight 1 - sum a_i.
struct WsggGrayGases {
  std::array<double, 4> a{};
  std::array<double, 4> kappa{};
};

// p_set_atm must be one of the ten tabulated pressures (0.1, 0.5, 1, 2, 5,
// 10, 20, 30, 50, 80 atm) within 1e-9 relative; T in [300, 3000] K and Mr in
// [0.1, 4]. Unclamped: throws outside the range.
WsggGrayGases wsgg_bordbar2021_gray_gases(double p_set_atm, double T,
                                          double Mr);

// Emissivity and its partial derivatives.
struct WsggEmissivityResult {
  double eps = 0.0;         // [-]
  double deps_dT = 0.0;     // [1/K]
  double deps_dP = 0.0;     // [1/Pa]
  double deps_dx_h2o = 0.0; // [-]
  double deps_dx_co2 = 0.0; // [-]
  double deps_dl_b = 0.0;   // [1/m]
};

// eps of an H2O-CO2 mixture across the pressure sets (see above).
//   T_g   : gas temperature [K], > 0
//   P     : total pressure [Pa], >= 0
//   x_h2o, x_co2 : mole fractions [-], in [0, 1]
//   l_b   : path length or mean beam length [m], >= 0
// Throws on a non-finite or out-of-domain input.
WsggEmissivityResult gas_emissivity_bordbar2021(double T_g, double P,
                                                double x_h2o, double x_co2,
                                                double l_b);

// The same from a mixture's mass fractions Y (one per species), with the
// derivative chained through the mole fractions of H2O and CO2.
struct WsggMixtureEmissivityResult {
  double eps = 0.0;
  double deps_dT = 0.0;            // [1/K]
  double deps_dP = 0.0;            // [1/Pa]
  double deps_dl_b = 0.0;          // [1/m]
  std::vector<double> deps_dY;     // [-], one per species
};

WsggMixtureEmissivityResult
gas_emissivity_bordbar2021_mixture(double T_g, double P,
                                   const std::vector<double> &Y, double l_b);

} // namespace combaero
