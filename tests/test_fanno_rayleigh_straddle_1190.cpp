#include <gtest/gtest.h>

#include "../include/fanno_mach.h"
#include "../include/humidair.h"

#include <cstddef>
#include <initializer_list>
#include <string>
#include <vector>

using namespace combaero;

// ---------------------------------------------------------------------------
// The cold scan of fanno_rayleigh_channel_flow (no warm start) walks log G up
// in ln-2 brackets and Brent-solves the first bracket whose excess (exit
// pressure minus the target) goes from positive to non-positive. On a cooled
// duct whose needed drive folds, that bracket can hold a fall, a rise and a
// later fall: Brent returned the third crossing, not the first fall from
// zero flow, and a warm start from below it returned another flux (#1190).
// ---------------------------------------------------------------------------

namespace {

const std::string kModel = "haaland";
// Relative tolerance on the issue duct's roots: 30x their measured noise band
// (see IssueDriveReturnsTheFirstFall).
constexpr double kRootBand = 1e-5;

struct Duct {
  double Pt0, Tt0, Tt_exit, L, D;
  bool exit_total;
};

// The issue's duct (static exit): the needed drive peaks at 2.15199 Pa near
// 0.0161 Gc and dips to 2.12681 Pa near 0.0197 Gc, its laminar range.
const Duct kIssue{199111.5, 868.17, 348.98, 0.09662, 9.1945e-3, false};
// Random cooled ducts (seed 41 of the #1190 sweep) whose fold the ends'
// slopes of the straddling bracket do not show: a fold near the choke on a
// total-pressure exit (peak 1536.73 Pa near 0.806 Gc, valley 1524.74 Pa near
// 0.904 Gc), and a laminar fold (peak 2.01274 Pa near 0.00553 Gc).
const Duct kNearChoke{177138.87128520015, 1137.3195556976007, 606.410788967628,
                      0.10507669919469245, 0.004186565305486364, true};
const Duct kLaminar{2330564.792012615, 1122.0160934321832, 972.1390722135684,
                    0.02524235838620738, 0.004121983069797663, true};

FannoChannelFlow flow(const Duct& d, double drive, double G_guess = -1.0) {
  return fanno_rayleigh_channel_flow(d.Pt0, d.Tt0, d.Tt_exit, dry_air(), d.Pt0 - drive,
                                     d.exit_total, d.L, d.D, 0.0, kModel, 1.0, G_guess, {},
                                     false);
}

double choked_flux(const Duct& d) {
  return fanno_rayleigh_choked_mass_flux(d.Pt0, d.Tt0, d.Tt_exit, dry_air(), d.L, d.D, 0.0,
                                         kModel);
}

// The drive the duct needs at flux G: Pt0 - P_exit (Pt_exit for a
// total-pressure exit) of the public march.
double needed(const Duct& d, double G) {
  const FannoRayleighMarch m =
      fanno_rayleigh_march(G, d.Pt0, d.Tt0, d.Tt_exit, dry_air(), d.L, d.D, 0.0, kModel);
  return d.Pt0 - (d.exit_total ? m.Pt_exit : m.P_exit);
}

// Reference: the first flux on an n-point grid over [a, b] Gc whose needed
// drive reaches the target, refined by 60 bisections. The grid is marched
// once per duct and reused across drives.
struct Reference {
  const Duct& d;
  std::vector<double> G, need;
  Reference(const Duct& duct, double a, double b, int n) : d(duct) {
    const double Gc = choked_flux(d);
    for (int i = 0; i <= n; ++i) {
      G.push_back(Gc * (a + (b - a) * i / n));
      need.push_back(needed(d, G.back()));
    }
  }
  double first_fall(double drive) const {
    EXPECT_LT(need.front(), drive);
    std::size_t i = 1;
    while (i + 1 < G.size() && need[i] < drive) {
      ++i;
    }
    double lo = G[i - 1], hi = G[i];
    for (int k = 0; k < 60; ++k) {
      const double m = 0.5 * (lo + hi);
      (needed(d, m) >= drive ? hi : lo) = m;
    }
    return hi;
  }
};

}  // namespace

TEST(FannoRayleighStraddle1190, IssueDriveReturnsTheFirstFall) {
  const double Gc = choked_flux(kIssue);
  const FannoChannelFlow c = flow(kIssue, 2.145);
  EXPECT_FALSE(c.choked);
  // The issue's table: the first fall at 0.0150235 Gc; the third crossing,
  // which the scan returned, at 0.021177 Gc.
  EXPECT_NEAR(c.G / Gc, 0.0150235, 1e-6);
  // The root itself is defined only to the march's rounding at 2 Pa of
  // drive on 2 bar: needed - drive changes sign within up to 3.3e-7 of it in
  // ln G over this duct's drives (measured, macOS clang), and another
  // compiler may pick any of those crossings. kRootBand leaves 30x that,
  // far below the 41% jump to the third crossing.
  EXPECT_NEAR(c.G / Reference(kIssue, 0.005, 0.035, 3000).first_fall(2.145), 1.0, kRootBand);
  // A warm start from below the first fall returns it too.
  EXPECT_NEAR(flow(kIssue, 2.145, 0.0145 * Gc).G / c.G, 1.0, kRootBand);
}

TEST(FannoRayleighStraddle1190, IssueFoldBandReturnsTheFirstFall) {
  // Drives over the band between the fold's valley and its peak: the scan
  // returned the third crossing from 2.1378 Pa up to the peak (81 of 200
  // such drives).
  const double Gc = choked_flux(kIssue);
  const Reference ref(kIssue, 0.005, 0.035, 3000);
  for (int i = 0; i < 24; ++i) {
    const double drive = 2.1270 + (2.1519 - 2.1270) * i / 23.0;
    const FannoChannelFlow c = flow(kIssue, drive);
    EXPECT_FALSE(c.choked) << drive;
    EXPECT_LT(c.G, 0.0162 * Gc) << drive;  // below the peak's flux
    EXPECT_NEAR(c.G / ref.first_fall(drive), 1.0, kRootBand) << drive;
  }
}

TEST(FannoRayleighStraddle1190, FoldsTheEndSlopesMissAreFound) {
  // The cubic Hermite of the straddling bracket's own ends showed no turn on
  // these drives; its two halves, each with the same test, do (one part,
  // the ends alone, returns 1.23x and 1.31x the first fall here). The
  // tolerance follows each root's own noise band, measured as the spread in
  // ln G of needed - drive's sign changes around it: <= 1e-8 on the
  // near-choke duct, so 1e-6; +-1.64e-6 (55 sign changes) on the laminar
  // one, 0.0006 Pa below its fold's peak where the excess is flat, so 1e-4.
  struct Case {
    const Duct* d;
    double drive, a, b;
    int n;
    double tol;
  };
  for (const Case& k : {Case{&kNearChoke, 1526.0, 0.5, 0.93, 430, 1e-6},
                        Case{&kNearChoke, 1528.0, 0.5, 0.93, 430, 1e-6},
                        Case{&kNearChoke, 1530.0, 0.5, 0.93, 430, 1e-6},
                        Case{&kLaminar, 2.01216, 0.004, 0.006, 400, 1e-4}}) {
    const FannoChannelFlow c = flow(*k.d, k.drive);
    EXPECT_FALSE(c.choked) << k.drive;
    EXPECT_NEAR(c.G / Reference(*k.d, k.a, k.b, k.n).first_fall(k.drive), 1.0, k.tol)
        << k.drive;
  }
}
