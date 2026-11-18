#include <gtest/gtest.h>

#include "../include/fanno_mach.h"
#include "../include/humidair.h"

#include <initializer_list>
#include <string>

using namespace combaero;

// ---------------------------------------------------------------------------
// The cold scan of fanno_rayleigh_channel_flow gives each cooled ln-2 bracket
// whose ends are both positive the Hermite test of its ends' values and
// slopes (#1167). A dip of the exit-pressure excess below zero that is
// narrower than the bracket and sits near one end of it leaves those values
// and slopes looking monotone: the scan moved on and returned a later fall
// (#1206). Such a bracket is now walked in parts with the same test wherever
// the needed drive at an end reaches a quarter of the drive. The #1191 warm
// check runs the same scan, so a warm call returns the same first fall.
// ---------------------------------------------------------------------------

namespace {

struct Duct {
  double Pt0, Tt0, Tt_exit, L, D, f_mult;
  bool exit_total;
};

// The issue's duct (dry air, haaland, total-pressure exit): the needed drive
// peaks at ~1536.7 Pa near 0.81 Gc and dips to ~1524.7 Pa near 0.90 Gc. For
// drives in the top 4.5 Pa of that fold the first fall lies in a bracket
// whose ends are both positive; the cold call returned the crossing near
// 0.934-0.940 Gc (1.18-1.22x the first fall).
const Duct kIssue{177138.87, 1137.32, 606.41, 0.105077, 4.18657e-3, 1.0, true};
// The #1219 review duct (f_multiplier 1.04505): the cold call returned
// 0.029521 Gc (1.48x the first fall), while #1191's warm check found the
// first fall: warm and cold disagreed.
const Duct kReview{105787.80, 1036.195, 931.352, 0.086848, 20.196e-3, 1.04505, true};
constexpr double kReviewDrive = 1.1438082;

// The first falls, from a 6000-point grid of the needed drive bisected to
// 1e-15 in ln G (tmp/1206/gt/repro.py). The cold roots land within 3e-11
// (issue) and 2e-9 (review, 1 Pa of drive near a fold's peak) of them; 1e-6
// is far below the 18-48 % jumps to the later falls.
constexpr double kBand = 1e-6;

FannoChannelFlow flow(const Duct& d, double drive, double G_guess) {
  return fanno_rayleigh_channel_flow(d.Pt0, d.Tt0, d.Tt_exit, dry_air(), d.Pt0 - drive,
                                     d.exit_total, d.L, d.D, 0.0, "haaland", d.f_mult, G_guess,
                                     {}, false);
}

double choked_flux(const Duct& d) {
  return fanno_rayleigh_choked_mass_flux(d.Pt0, d.Tt0, d.Tt_exit, dry_air(), d.L, d.D, 0.0,
                                         "haaland", d.f_mult);
}

}  // namespace

TEST(FannoRayleighBothPositive1206, IssueColdCallsReturnTheFirstFall) {
  const double Gc = choked_flux(kIssue);
  const struct {
    double drive, first;
  } cases[] = {{1532.2, 0.767108500}, {1533.5, 0.773092693}, {1535.0, 0.781905364},
               {1536.5, 0.797662498}};
  for (const auto& c : cases) {
    const FannoChannelFlow r = flow(kIssue, c.drive, -1.0);
    EXPECT_FALSE(r.choked) << c.drive;
    EXPECT_NEAR(r.G / (c.first * Gc), 1.0, kBand) << c.drive;
    // Warm starts below, at and above the first fall, and at the later fall.
    for (const double s : {0.5, 0.99, 1.01, 1.3}) {
      const FannoChannelFlow w = flow(kIssue, c.drive, s * r.G);
      EXPECT_FALSE(w.choked) << c.drive << " " << s;
      EXPECT_NEAR(w.G / r.G, 1.0, kBand) << c.drive << " " << s;
    }
  }
}

TEST(FannoRayleighBothPositive1206, ReviewDuctWarmAndColdAgree) {
  const double Gc = choked_flux(kReview);
  const FannoChannelFlow c = flow(kReview, kReviewDrive, -1.0);
  ASSERT_FALSE(c.choked);
  EXPECT_NEAR(c.G / (0.019942574 * Gc), 1.0, kBand);
  // From the cold root's old value (0.029521 Gc) and around the first fall.
  for (const double s : {0.010, 0.0199, 0.0201, 0.029521}) {
    const FannoChannelFlow w = flow(kReview, kReviewDrive, s * Gc);
    EXPECT_FALSE(w.choked) << s;
    EXPECT_NEAR(w.G / c.G, 1.0, kBand) << s;
  }
}
