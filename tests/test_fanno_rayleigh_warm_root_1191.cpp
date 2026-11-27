#include <gtest/gtest.h>

#include "../include/fanno_mach.h"
#include "../include/humidair.h"

#include <initializer_list>
#include <string>

using namespace combaero;

// ---------------------------------------------------------------------------
// A warm start of fanno_rayleigh_channel_flow (G_guess > 0, the element's
// previous flux) brackets the root around the guess. On a cooled duct whose
// needed drive rises again after a fold's valley, the exit-pressure excess
// falls, rises and falls again, and a guess within 6.4 % above the later
// fall bracketed that one: the warm call returned a later crossing while the
// cold scan returned the first fall, so the flux depended on the solver's
// history (#1191). A cooled warm root is now checked by the cold scan run
// up to the warm bracket's lower end.
// ---------------------------------------------------------------------------

namespace {

const std::string kModel = "haaland";

struct Duct {
  double Pt0, Tt0, Tt_exit, L, D;
  bool exit_total;
};

// The issue's duct (total-pressure exit, 6797.70 Pa of drive): the needed
// drive rises to ~7046 Pa near 0.80 Gc, falls to ~6777 Pa near 0.94 Gc and
// rises again toward the choke, so the excess falls through zero at
// 0.666323 Gc (the first fall) and again at 0.957556 Gc.
const Duct kIssue{1.28643e6, 1984.90, 1248.21, 0.15866, 5.8745e-3, true};
constexpr double kIssueDrive = 6797.70;
// The #1190 duct (static exit, 2.145 Pa of drive): first fall at
// 0.0150235 Gc, a later one at 0.021177 Gc inside its laminar fold.
const Duct kLaminar{199111.5, 868.17, 348.98, 0.09662, 9.1945e-3, false};
constexpr double kLaminarDrive = 2.145;

// The issue duct's root is defined to the march's rounding: needed - drive
// changes sign 43 times within -5.3e-12..+6.3e-12 of it in ln G (measured,
// macOS clang), and the warm and cold calls land 0.8e-12 to 1.9e-12 apart.
// 1e-9 leaves 80x the band, far below the 44% jump to the later fall. The
// laminar duct's band is 3.3e-7 (test_fanno_rayleigh_straddle_1190.cpp).
constexpr double kIssueBand = 1e-9;
constexpr double kLaminarBand = 1e-5;

FannoChannelFlow flow(const Duct& d, double drive, double G_guess) {
  return fanno_rayleigh_channel_flow(d.Pt0, d.Tt0, d.Tt_exit, dry_air(), d.Pt0 - drive,
                                     d.exit_total, d.L, d.D, 0.0, kModel, 1.0, G_guess, {},
                                     true);
}

double choked_flux(const Duct& d) {
  return fanno_rayleigh_choked_mass_flux(d.Pt0, d.Tt0, d.Tt_exit, dry_air(), d.L, d.D, 0.0,
                                         kModel);
}

}  // namespace

TEST(FannoRayleighWarmRoot1191, IssueWarmStartsReturnTheFirstFall) {
  const double Gc = choked_flux(kIssue);
  const FannoChannelFlow c = flow(kIssue, kIssueDrive, -1.0);
  ASSERT_FALSE(c.choked);
  EXPECT_NEAR(c.G / Gc, 0.666323, 1e-6);
  // 0.95, 0.99 and 0.999 Gc returned 0.957556 Gc (the later fall) before.
  for (const double s : {0.05, 0.5, 0.9, 0.95, 0.99, 0.999}) {
    const FannoChannelFlow w = flow(kIssue, kIssueDrive, s * Gc);
    EXPECT_FALSE(w.choked) << s;
    EXPECT_NEAR(w.G / c.G, 1.0, kIssueBand) << s;
    // The derivatives follow the returned root: those of the cold call.
    EXPECT_NEAR(w.dG_dP_target / c.dG_dP_target, 1.0, 1e-6) << s;
    EXPECT_NEAR(w.dG_dPt0 / c.dG_dPt0, 1.0, 1e-6) << s;
    EXPECT_NEAR(w.dG_dTt_exit / c.dG_dTt_exit, 1.0, 1e-6) << s;
  }
}

TEST(FannoRayleighWarmRoot1191, LaminarFoldWarmStartAboveTheLaterFall) {
  // A guess at 0.02 Gc, 5.6% below the later fall at 0.021177 Gc, returned
  // that fall before.
  const double Gc = choked_flux(kLaminar);
  const FannoChannelFlow c = flow(kLaminar, kLaminarDrive, -1.0);
  ASSERT_FALSE(c.choked);
  EXPECT_NEAR(c.G / Gc, 0.0150297, 1e-6);
  for (const double s : {0.0145, 0.017, 0.02, 0.0215}) {
    const FannoChannelFlow w = flow(kLaminar, kLaminarDrive, s * Gc);
    EXPECT_FALSE(w.choked) << s;
    EXPECT_NEAR(w.G / c.G, 1.0, kLaminarBand) << s;
  }
}

TEST(FannoRayleighWarmRoot1191, HeatedAndAdiabaticWarmStartsAgreeWithCold) {
  // The check runs on cooled ducts only; heated and adiabatic ducts cannot
  // turn (no Rayleigh gain), so a warm start from anywhere below the inlet's
  // sonic flux returns the cold root. Bit identity with the previous kernel
  // is checked by the stress F,J hashes, not here.
  for (const double Tte : {1984.90, 2400.0}) {
    const Duct d{kIssue.Pt0, kIssue.Tt0, Tte, kIssue.L, kIssue.D, true};
    const double Gc = choked_flux(d);
    const FannoChannelFlow c = flow(d, kIssueDrive, -1.0);
    ASSERT_FALSE(c.choked) << Tte;
    for (const double s : {0.05, 0.5, 0.9, 0.99}) {
      const FannoChannelFlow w = flow(d, kIssueDrive, s * Gc);
      EXPECT_NEAR(w.G / c.G, 1.0, kIssueBand) << Tte << " " << s;
    }
  }
}
