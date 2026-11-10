#include "wsgg_bordbar2021.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "thermo_transport_data.h"
#include "wsgg_bordbar2021_data.h"

namespace combaero {

namespace {

namespace wd = wsgg_bordbar2021;

// A tabulated pressure matches a request to this relative tolerance.
constexpr double kSetMatchRelTol = 1e-9;

void require_finite(double v, const char *what) {
  if (!std::isfinite(v)) {
    throw std::invalid_argument(std::string(what) + " must be finite");
  }
}

// Eqs. 9-11 at one set, with d/dT and d/dMr of the weights and kappas.
struct SetGases {
  std::array<double, 4> a{}, da_dT{}, da_dMr{};
  std::array<double, 4> kappa{}, dkappa_dMr{};
};

SetGases set_gases(int s, double T, double Mr) {
  // Mr^(4-k), k = 0..4, and its derivative.
  std::array<double, wd::N_MR_TERMS> mp{};
  std::array<double, wd::N_MR_TERMS> dmp{};
  for (int k = 0; k < wd::N_MR_TERMS; ++k) {
    const int n = wd::N_MR_TERMS - 1 - k;
    mp[k] = std::pow(Mr, n);
    dmp[k] = n > 0 ? n * std::pow(Mr, n - 1) : 0.0;
  }
  const double th = T / WSGG_BORDBAR2021_T_REF;
  std::array<double, wd::N_T_TERMS> tp{};
  std::array<double, wd::N_T_TERMS> dtp{};
  for (int j = 0; j < wd::N_T_TERMS; ++j) {
    tp[j] = std::pow(th, j);
    dtp[j] = j > 0 ? j * std::pow(th, j - 1) / WSGG_BORDBAR2021_T_REF : 0.0;
  }
  SetGases g;
  for (int i = 0; i < wd::N_GRAY_GASES; ++i) {
    double kap = 0.0;
    double dkap = 0.0;
    for (int k = 0; k < wd::N_MR_TERMS; ++k) {
      kap += wd::D[s][i][k] * mp[k];
      dkap += wd::D[s][i][k] * dmp[k];
    }
    double a = 0.0;
    double daT = 0.0;
    double daM = 0.0;
    for (int j = 0; j < wd::N_T_TERMS; ++j) {
      double b = 0.0;
      double db = 0.0;
      for (int k = 0; k < wd::N_MR_TERMS; ++k) {
        b += wd::C[s][i][j][k] * mp[k];
        db += wd::C[s][i][j][k] * dmp[k];
      }
      a += b * tp[j];
      daT += b * dtp[j];
      daM += db * tp[j];
    }
    g.a[i] = a;
    g.da_dT[i] = daT;
    g.da_dMr[i] = daM;
    g.kappa[i] = kap;
    g.dkappa_dMr[i] = dkap;
  }
  return g;
}

// Eq. 7 at one set and a partial-pressure path S = (p_H2O + p_CO2) L
// [atm m]: eps and its derivatives in T, Mr and S.
struct SetEps {
  double eps = 0.0, dT = 0.0, dMr = 0.0, dS = 0.0;
};

SetEps set_eps(int s, double T, double Mr, double S) {
  const SetGases g = set_gases(s, T, Mr);
  SetEps r;
  for (int i = 0; i < wd::N_GRAY_GASES; ++i) {
    const double e = std::exp(-g.kappa[i] * S);
    const double om = 1.0 - e;
    r.eps += g.a[i] * om;
    r.dT += g.da_dT[i] * om;
    r.dMr += g.da_dMr[i] * om + g.a[i] * e * S * g.dkappa_dMr[i];
    r.dS += g.a[i] * g.kappa[i] * e;
  }
  return r;
}

} // namespace

WsggGrayGases wsgg_bordbar2021_gray_gases(double p_set_atm, double T,
                                          double Mr) {
  require_finite(p_set_atm, "p_set_atm");
  require_finite(T, "T");
  require_finite(Mr, "Mr");
  int s = -1;
  for (int n = 0; n < wd::N_PRESSURE_SETS; ++n) {
    if (std::abs(p_set_atm - wd::PRESSURE_SETS_ATM[n]) <=
        kSetMatchRelTol * wd::PRESSURE_SETS_ATM[n]) {
      s = n;
    }
  }
  if (s < 0) {
    throw std::invalid_argument(
        "wsgg_bordbar2021_gray_gases: p_set_atm must be a tabulated pressure "
        "(0.1, 0.5, 1, 2, 5, 10, 20, 30, 50, 80 atm)");
  }
  if (T < WSGG_BORDBAR2021_T_MIN || T > WSGG_BORDBAR2021_T_MAX) {
    throw std::invalid_argument(
        "wsgg_bordbar2021_gray_gases: T must be in [300, 3000] K");
  }
  if (Mr < WSGG_BORDBAR2021_MR_MIN || Mr > WSGG_BORDBAR2021_MR_MAX) {
    throw std::invalid_argument(
        "wsgg_bordbar2021_gray_gases: Mr must be in [0.1, 4]");
  }
  const SetGases g = set_gases(s, T, Mr);
  WsggGrayGases r;
  r.a = g.a;
  r.kappa = g.kappa;
  return r;
}

WsggEmissivityResult gas_emissivity_bordbar2021(double T_g, double P,
                                                double x_h2o, double x_co2,
                                                double l_b) {
  require_finite(T_g, "T_g");
  require_finite(P, "P");
  require_finite(x_h2o, "x_h2o");
  require_finite(x_co2, "x_co2");
  require_finite(l_b, "l_b");
  if (!(T_g > 0.0)) {
    throw std::invalid_argument("gas_emissivity_bordbar2021: T_g must be > 0");
  }
  if (P < 0.0 || l_b < 0.0) {
    throw std::invalid_argument(
        "gas_emissivity_bordbar2021: P and l_b must be >= 0");
  }
  if (x_h2o < 0.0 || x_co2 < 0.0 || x_h2o + x_co2 > 1.0 + 1e-12) {
    throw std::invalid_argument(
        "gas_emissivity_bordbar2021: x_h2o, x_co2 >= 0 with x_h2o + x_co2 <= 1");
  }

  // Clamped temperature and mole ratio (see the header).
  double T = T_g;
  double dT_dTg = 1.0;
  if (T < WSGG_BORDBAR2021_T_MIN || T > WSGG_BORDBAR2021_T_MAX) {
    T = std::clamp(T, WSGG_BORDBAR2021_T_MIN, WSGG_BORDBAR2021_T_MAX);
    dT_dTg = 0.0;
  }
  const double Mr_raw = x_co2 > 0.0
                            ? x_h2o / x_co2
                            : std::numeric_limits<double>::infinity();
  double Mr = Mr_raw;
  double dMr_dxh = 0.0;
  double dMr_dxc = 0.0;
  if (Mr_raw < WSGG_BORDBAR2021_MR_MIN) {
    Mr = WSGG_BORDBAR2021_MR_MIN;
  } else if (Mr_raw > WSGG_BORDBAR2021_MR_MAX) {
    Mr = WSGG_BORDBAR2021_MR_MAX;
  } else {
    dMr_dxh = 1.0 / x_co2;
    dMr_dxc = -x_h2o / (x_co2 * x_co2);
  }

  // Partial-pressure path S = (x_H2O + x_CO2) (P / atm) l_b [atm m].
  const double P_atm = P / WSGG_PA_PER_ATM;
  const double xs = x_h2o + x_co2;
  const double S = xs * P_atm * l_b;

  // Bracketing sets and the ln P weight on the upper one.
  const int n_sets = wd::N_PRESSURE_SETS;
  int lo = 0;
  int hi = 0;
  double w = 0.0;
  double dw_dP = 0.0;
  if (!(P_atm > wd::PRESSURE_SETS_ATM[0])) {
    lo = hi = 0;
  } else if (P_atm >= wd::PRESSURE_SETS_ATM[n_sets - 1]) {
    lo = hi = n_sets - 1;
  } else {
    hi = 1;
    while (wd::PRESSURE_SETS_ATM[hi] < P_atm) {
      ++hi;
    }
    lo = hi - 1;
    const double ln_ratio =
        std::log(wd::PRESSURE_SETS_ATM[hi] / wd::PRESSURE_SETS_ATM[lo]);
    w = std::log(P_atm / wd::PRESSURE_SETS_ATM[lo]) / ln_ratio;
    dw_dP = 1.0 / (P * ln_ratio);
  }
  SetEps e = set_eps(lo, T, Mr, S);
  if (hi != lo && w > 0.0) {
    const SetEps eh = set_eps(hi, T, Mr, S);
    const double d_eps = eh.eps - e.eps;
    e.eps += w * d_eps;
    e.dT += w * (eh.dT - e.dT);
    e.dMr += w * (eh.dMr - e.dMr);
    e.dS += w * (eh.dS - e.dS);
    // The weight's own pressure derivative, added below in Pa.
    dw_dP *= d_eps;
  } else {
    dw_dP = 0.0;
  }

  WsggEmissivityResult r;
  if (e.eps <= 0.0 || e.eps >= 1.0) {
    r.eps = std::clamp(e.eps, 0.0, 1.0);
    return r;
  }
  r.eps = e.eps;
  r.deps_dT = e.dT * dT_dTg;
  const double dS_dx = P_atm * l_b;
  r.deps_dx_h2o = e.dMr * dMr_dxh + e.dS * dS_dx;
  r.deps_dx_co2 = e.dMr * dMr_dxc + e.dS * dS_dx;
  r.deps_dP = e.dS * xs * l_b / WSGG_PA_PER_ATM + dw_dP;
  r.deps_dl_b = e.dS * xs * P_atm;
  return r;
}

WsggMixtureEmissivityResult
gas_emissivity_bordbar2021_mixture(double T_g, double P,
                                   const std::vector<double> &Y, double l_b) {
  const std::size_t n = molar_masses.size();
  if (Y.size() != n) {
    throw std::invalid_argument("gas_emissivity_bordbar2021_mixture: Y needs "
                                "one mass fraction per species");
  }
  const std::size_t ih = static_cast<std::size_t>(species_index.at("H2O"));
  const std::size_t ic = static_cast<std::size_t>(species_index.at("CO2"));
  // x_k = (Y_k / M_k) / N, N = sum_j Y_j / M_j.
  double N = 0.0;
  for (std::size_t k = 0; k < n; ++k) {
    require_finite(Y[k], "Y");
    if (Y[k] < 0.0) {
      throw std::invalid_argument(
          "gas_emissivity_bordbar2021_mixture: Y must be >= 0");
    }
    N += Y[k] / molar_masses[k];
  }
  if (!(N > 0.0)) {
    throw std::invalid_argument(
        "gas_emissivity_bordbar2021_mixture: Y must not be all zero");
  }
  const double xh = Y[ih] / molar_masses[ih] / N;
  const double xc = Y[ic] / molar_masses[ic] / N;
  const WsggEmissivityResult e =
      gas_emissivity_bordbar2021(T_g, P, xh, xc, l_b);
  WsggMixtureEmissivityResult r;
  r.eps = e.eps;
  r.deps_dT = e.deps_dT;
  r.deps_dP = e.deps_dP;
  r.deps_dl_b = e.deps_dl_b;
  // dx_m/dY_j = (delta_mj / M_m - x_m / M_j) / N.
  r.deps_dY.assign(n, 0.0);
  const double g = e.deps_dx_h2o * xh + e.deps_dx_co2 * xc;
  for (std::size_t j = 0; j < n; ++j) {
    r.deps_dY[j] = -g / (molar_masses[j] * N);
  }
  r.deps_dY[ih] += e.deps_dx_h2o / (molar_masses[ih] * N);
  r.deps_dY[ic] += e.deps_dx_co2 / (molar_masses[ic] * N);
  return r;
}

} // namespace combaero
