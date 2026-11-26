#include <gtest/gtest.h>

#include "../include/fanno_mach.h"
#include "../include/humidair.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using namespace combaero;

// ---------------------------------------------------------------------------
// A cooled duct's needed drive can rise from zero flow to a peak and fall
// back (the Rayleigh Pt gain), so at one drive it can run subsonic while
// choked flow holds too (#1307). fanno_rayleigh_channel_flow keeps the
// start-up root; fanno_rayleigh_choked_alternative says whether choked flow
// also meets the same boundary values. Case from the issue: Pt0 1e5 Pa,
// Tt0 400 K, L = D = 1 cm, total-pressure exit, air; the band at a 1.5 Pa
// drive is Tt_exit ~389-393 K.
// ---------------------------------------------------------------------------

namespace {

const std::string kModel = "haaland";
constexpr double kPt0 = 1.0e5;
constexpr double kTt0 = 400.0;
constexpr double kL = 0.01;
constexpr double kD = 0.01;

FannoChannelFlow flow(double Te, double drive) {
  return fanno_rayleigh_channel_flow(kPt0, kTt0, Te, dry_air(), kPt0 - drive, true, kL, kD, 0.0,
                                     kModel, 1.0, -1.0, {}, false);
}

FannoChokedAlternative alt(double Te, double drive) {
  return fanno_rayleigh_choked_alternative(kPt0, kTt0, Te, dry_air(), kPt0 - drive, true, kL, kD,
                                           0.0, kModel);
}

}  // namespace

// Inside the band: the start-up root is subsonic, and choked flow at the
// duct's own choked flux also meets the target, an order of magnitude more
// flow.
TEST(FannoRayleighChokedAlternative1307, InsideTheBandChokedFlowAlsoHolds) {
  const FannoChannelFlow sub = flow(391.0, 1.5);
  const FannoChokedAlternative a = alt(391.0, 1.5);
  ASSERT_FALSE(sub.choked);
  EXPECT_TRUE(a.admissible);
  EXPECT_GT(a.exit_margin, 0.0);
  const double Gc = fanno_rayleigh_choked_mass_flux(kPt0, kTt0, 391.0, dry_air(), kL, kD, 0.0,
                                                    kModel);
  EXPECT_NEAR(a.G, Gc, 1e-12 * Gc);
  EXPECT_GT(a.G, 5.0 * sub.G);
}

// Where the channel flow itself chokes, the alternative is that same state:
// its verdict is the cold path's own choke test.
TEST(FannoRayleighChokedAlternative1307, AChokedChannelFlowIsItsOwnAlternative) {
  for (double Te : {380.0, 391.0, 400.0, 420.0}) {
    const FannoChannelFlow choked = flow(Te, 5.0e4);
    const FannoChokedAlternative a = alt(Te, 5.0e4);
    ASSERT_TRUE(choked.choked) << Te;
    EXPECT_TRUE(a.admissible) << Te;
    EXPECT_NEAR(a.G, choked.G, 1e-10 * choked.G) << Te;
  }
}

// Adiabatic, heated, and cooled too little: the needed drive is monotone in
// G up to the choke, so the choked exit pressure is below the target.
TEST(FannoRayleighChokedAlternative1307, NotAdmissibleWhereTheDriveIsMonotone) {
  for (double Te : {400.0, 420.0, 396.0}) {
    for (double drive : {1.0, 5.0}) {
      const FannoChokedAlternative a = alt(Te, drive);
      EXPECT_FALSE(a.admissible) << Te << " " << drive;
      EXPECT_LT(a.exit_margin, 0.0) << Te << " " << drive;
    }
  }
}

// No forward drive: the channel flow is no flow, and there is no forward
// choked state to compare (nothing evaluated).
TEST(FannoRayleighChokedAlternative1307, NoForwardDriveHasNoAlternative) {
  for (double drive : {0.0, -5.0}) {
    const FannoChokedAlternative a = alt(391.0, drive);
    EXPECT_FALSE(a.admissible) << drive;
    EXPECT_EQ(a.G, 0.0) << drive;
    EXPECT_TRUE(std::isnan(a.exit_margin)) << drive;
    EXPECT_EQ(flow(391.0, drive).G, 0.0) << drive;
  }
}

TEST(FannoRayleighChokedAlternative1307, RefusesInvalidInput) {
  const std::vector<double> X = dry_air();
  EXPECT_THROW(fanno_rayleigh_choked_alternative(-1.0, kTt0, 391.0, X, 0.0, true, kL, kD, 0.0,
                                                 kModel),
               std::invalid_argument);
  EXPECT_THROW(fanno_rayleigh_choked_alternative(kPt0, kTt0, 0.0, X, kPt0 - 1.5, true, kL, kD, 0.0,
                                                 kModel),
               std::invalid_argument);
  EXPECT_THROW(fanno_rayleigh_choked_alternative(kPt0, kTt0, 391.0, X, std::nan(""), true, kL, kD,
                                                 0.0, kModel),
               std::invalid_argument);
}
