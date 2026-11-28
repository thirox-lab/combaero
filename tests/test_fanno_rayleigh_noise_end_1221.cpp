#include <gtest/gtest.h>

#include "../include/fanno_mach.h"
#include "../include/humidair.h"

#include <string>

using namespace combaero;

// ---------------------------------------------------------------------------
// A cooled duct's exit-pressure excess can fall through zero, dip below it
// and rise back (a fold). At the dip's rising crossing the excess is a root
// too, and within the march's rounding noise of it its sign is noise. A warm
// bracket that ended there was accepted as a sign change: Brent converged on
// a noise crossing next to that end, the dip's rising crossing, instead of
// the first fall below it (#1221). A warm bracket with an end inside the
// noise now goes the cold way, so the warm call returns the cold root.
// ---------------------------------------------------------------------------

namespace {

struct Duct {
  double Pt0, Tt0, Tt_exit, L, D, roughness, f_mult;
  bool exit_total;
  const char* model;
};

struct Case {
  Duct duct;
  double drive;  // Pt0 - P_target [Pa]
  double first;  // the first fall [kg/(m^2 s)]
  double guess;  // places a warm bracket end on the dip's rising crossing
};

// From the #1206 fold sweep's ground truth (rsweep ducts (seed, k); dense
// needed-drive grid, folds golden-refined, the first fall and rising crossing
// bisected to machine precision against the kernel's own march; no kernel
// bracket logic). Drives 0.99999 of the way up the fold's band, guess 1.001x
// the dip's rising crossing: before #1221, main returned the rising crossing.
const Case kR1008{{105787.80410853265, 1036.195206507452, 931.3523221535137, 0.08684772270785689,
                   0.020196137042498154, 0.0, 1.045051934528771, true, "haaland"},
                  1.1553241394866045,
                  2.9285782769082354,
                  2.9315908706208385 * 1.001};
const Case kR3092{{1926908.6296252967, 1923.1165135646884, 1364.7275630384754, 0.11627250252396268,
                   0.009383313294004455, 0.0, 1.0, true, "haaland"},
                  2.2949718273158304,
                  7.4663566815227771,
                  7.4855788075446235 * 1.001};
// A deeper dip (0.7 mPa, 0.99 of the band; rsweep (4, 823), static exit):
// the guess puts the last widening's upper end (1.064x the guess) on the
// rising crossing. Before #1221, main returned that crossing.
const Case kR4823{{376443.11135940684, 1336.1805388982839, 434.1327334181294, 0.30493706715635316,
                   0.016752827082123504, 1.9024691672720516e-07, 1.0, false, "haaland"},
                  1.1048811150901019,
                  2.7259658442912285,
                  2.8718825003648591 / 1.064};

// The cold roots land within 4e-6 of the references (root rounding near the
// fold's peak, where the excess is flat); the misses were 1e-3 to 5e-2 off.
constexpr double kBand = 1e-5;

FannoChannelFlow flow(const Case& c, double G_guess) {
  const Duct& d = c.duct;
  return fanno_rayleigh_channel_flow(d.Pt0, d.Tt0, d.Tt_exit, dry_air(), d.Pt0 - c.drive,
                                     d.exit_total, d.L, d.D, d.roughness, d.model, d.f_mult,
                                     G_guess, {}, false);
}

}  // namespace

TEST(FannoRayleighNoiseEnd1221, WarmBracketEndOnTheDipsRiseReturnsTheFirstFall) {
  for (const Case* c : {&kR1008, &kR3092, &kR4823}) {
    const FannoChannelFlow cold = flow(*c, -1.0);
    ASSERT_FALSE(cold.choked) << c->drive;
    EXPECT_NEAR(cold.G / c->first, 1.0, kBand) << c->drive;
    const FannoChannelFlow warm = flow(*c, c->guess);
    EXPECT_FALSE(warm.choked) << c->drive;
    EXPECT_NEAR(warm.G / c->first, 1.0, kBand) << c->drive;
    EXPECT_EQ(warm.G, cold.G) << c->drive;
  }
}
