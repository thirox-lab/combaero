#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "rib_correlation.h"

using combaero::cooling::evaluate_rib;
using combaero::cooling::han_1988_orthogonal;
using combaero::cooling::han_1989_narrow_channel;
using combaero::cooling::han_zhang_lee_1991;
using combaero::cooling::RibShape;
using combaero::cooling::han_park_1988_angled;
using combaero::cooling::rallabandi_2009_high_re;
using combaero::cooling::rib_channel_friction;
using combaero::cooling::rib_ribbed_wall_fraction;
using combaero::cooling::RibCorrelationSet;
using combaero::cooling::RibGeometry;
using combaero::cooling::validate_rib_set;

namespace {

RibGeometry ref_geometry() {
  RibGeometry g;
  g.e_D = 0.047;
  g.p_e = 10.0;
  g.W_H = 1.0;
  g.alpha_deg = 90.0;
  return g;
}

RibGeometry rallabandi_ref_geometry() {
  RibGeometry g;
  g.e_D = 0.1;
  g.p_e = 10.0;
  g.W_H = 1.0;
  g.alpha_deg = 45.0;
  return g;
}

}  // namespace

// The numbers below are the worked values recorded in the CONFIRMED extraction
// at validation/cooling/extractions/han_ribbed.md. They come from the document,
// not from this implementation, so this test can fail.
TEST(RibCorrelationTest, ReproducesTheConfirmedExtraction) {
  const auto set = han_1988_orthogonal();
  const auto g = ref_geometry();

  const auto lo = evaluate_rib(set, g, 10000.0);
  EXPECT_NEAR(lo.R, 3.2000, 1e-4);
  EXPECT_NEAR(lo.f, 0.04576, 1e-5);
  EXPECT_NEAR(lo.e_plus, 71.1, 0.1);
  EXPECT_NEAR(lo.G, 12.21, 0.01);
  EXPECT_NEAR(lo.St_r, 0.00968, 1e-5);

  const auto hi = evaluate_rib(set, g, 60000.0);
  EXPECT_NEAR(hi.e_plus, 426.6, 0.1);
  EXPECT_NEAR(hi.G, 20.16, 0.01);
  EXPECT_NEAR(hi.St_r, 0.00642, 1e-5);
}

// R carries no e+ term, so f cannot depend on Reynolds number. This is a
// structural property of the family, and it is what makes df/d(mdot) zero.
TEST(RibCorrelationTest, FrictionFactorIsIndependentOfReynoldsNumber) {
  const auto set = han_1988_orthogonal();
  const auto g = ref_geometry();
  const double f_ref = evaluate_rib(set, g, 10000.0).f;
  for (double Re : {1.0, 1e3, 1e5, 1e7}) {
    EXPECT_NEAR(evaluate_rib(set, g, Re).f, f_ref, 1e-12);
  }
}

// The solver probes states that are not physical. Every one of these must
// return finite numbers rather than NaN, an infinity, or a complex result.
TEST(RibCorrelationTest, GuardsHoldForStatesTheSolverActuallyProbes) {
  const auto set = han_1988_orthogonal();
  const auto g = ref_geometry();
  for (double Re : {-1e6, -1e4, -1.0, 0.0, 1e-12, 1.0, 1e8}) {
    const auto r = evaluate_rib(set, g, Re);
    EXPECT_TRUE(std::isfinite(r.f)) << "Re = " << Re;
    EXPECT_TRUE(std::isfinite(r.G)) << "Re = " << Re;
    EXPECT_TRUE(std::isfinite(r.St_r)) << "Re = " << Re;
    EXPECT_TRUE(std::isfinite(r.dSt_dRe)) << "Re = " << Re;
    EXPECT_GT(r.f, 0.0) << "Re = " << Re;
    EXPECT_GT(r.St_r, 0.0) << "Re = " << Re;
  }
}

// Reversing the flow must not change the magnitude of a heat-transfer
// quantity: h does not care which way the fluid goes. Only the sign of e+
// records direction.
TEST(RibCorrelationTest, ReverseFlowIsSymmetricInMagnitude) {
  const auto set = han_1988_orthogonal();
  const auto g = ref_geometry();
  const auto fwd = evaluate_rib(set, g, 30000.0);
  const auto rev = evaluate_rib(set, g, -30000.0);

  EXPECT_NEAR(fwd.G, rev.G, 1e-12);
  EXPECT_NEAR(fwd.St_r, rev.St_r, 1e-12);
  EXPECT_NEAR(fwd.e_plus, -rev.e_plus, 1e-9);
  // The derivative of an even function flips sign across zero.
  EXPECT_NEAR(fwd.dSt_dRe, -rev.dSt_dRe, 1e-12);
}

// The analytic derivative must match a central difference THROUGH the guard,
// not only in the far field where the guard is inert. A pin-fin Jacobian that
// was right everywhere except where it mattered is exactly what #328 fixed.
TEST(RibCorrelationTest, StantonDerivativeMatchesCentralDifferences) {
  const auto set = han_1988_orthogonal();
  const auto g = ref_geometry();
  // The last two sit where the e+ floor is active.
  for (double Re : {50000.0, 10000.0, 100.0, 10.0}) {
    const double h = std::max(1e-6, std::abs(Re) * 1e-6);
    const double fd = (evaluate_rib(set, g, Re + h).St_r -
                       evaluate_rib(set, g, Re - h).St_r) /
                      (2.0 * h);
    const double analytic = evaluate_rib(set, g, Re).dSt_dRe;
    const double scale = std::max({std::abs(fd), std::abs(analytic), 1e-14});
    EXPECT_LT(std::abs(analytic - fd) / scale, 1e-5) << "Re = " << Re;
  }
}

// Outside the source's ranges the evaluator reports, it does not refuse: the
// band belongs to the source's rig, not to the caller's hardware.
TEST(RibCorrelationTest, ValidityIsAdvisoryNotEnforced) {
  const auto set = han_1988_orthogonal();
  auto g = ref_geometry();
  EXPECT_FALSE(evaluate_rib(set, g, 30000.0).extrapolated);

  g.e_D = 0.20;  // well past the 0.047-0.078 band
  const auto r = evaluate_rib(set, g, 30000.0);
  EXPECT_TRUE(r.extrapolated);
  EXPECT_TRUE(std::isfinite(r.St_r));
  EXPECT_GT(r.St_r, 0.0);
}

// Bad PARAMETERS are a mistake and are rejected once, loudly -- the opposite
// treatment from a bad operating point.
TEST(RibCorrelationTest, MalformedSetsAreRejected) {
  EXPECT_NO_THROW(validate_rib_set(han_1988_orthogonal()));

  auto no_source = han_1988_orthogonal();
  no_source.source.clear();
  EXPECT_THROW(validate_rib_set(no_source), std::invalid_argument);

  auto zero_reference = han_1988_orthogonal();
  zero_reference.R_pe.reference = 0.0;
  EXPECT_THROW(validate_rib_set(zero_reference), std::invalid_argument);

  auto nan_constant = han_1988_orthogonal();
  nan_constant.C_G = std::nan("");
  EXPECT_THROW(validate_rib_set(nan_constant), std::invalid_argument);

  auto negative = han_1988_orthogonal();
  negative.C_R = -1.0;
  EXPECT_THROW(validate_rib_set(negative), std::invalid_argument);
}

// The numbers below are hand-computed from Eq. (17)/(18) as printed in
// Rallabandi, Yang and Han (2009), independently of this implementation --
// see validation/cooling/extractions/han_ribbed_high_re.md. R carries no e+
// term here either, so it repeats across Re exactly as for han_1988_orthogonal.
TEST(RibCorrelationTest, ReproducesTheConfirmedRallabandiExtraction) {
  const auto set = rallabandi_2009_high_re();
  const auto g = rallabandi_ref_geometry();

  const auto lo = evaluate_rib(set, g, 30000.0);
  EXPECT_NEAR(lo.R, 4.00939, 1e-4);
  EXPECT_NEAR(lo.f, 0.06533, 1e-4);
  EXPECT_NEAR(lo.e_plus, 542.20, 0.1);
  EXPECT_NEAR(lo.G, 16.1350, 0.01);
  EXPECT_NEAR(lo.St_r, 0.010235, 1e-5);

  const auto hi = evaluate_rib(set, g, 400000.0);
  EXPECT_NEAR(hi.e_plus, 7229.37, 0.1);
  EXPECT_NEAR(hi.G, 47.8899, 0.01);
  EXPECT_NEAR(hi.St_r, 0.003658, 1e-5);
}

// han_1988_orthogonal (90 deg) and rallabandi_2009_high_re (45 deg, sharp
// ribs) are different papers a caller must choose between explicitly. Their
// Reynolds-number bands actually OVERLAP (10e3-60e3 vs 30e3-400e3, sharing
// 30e3-60e3) -- Re alone does not distinguish them, which is exactly why
// there is no auto-switching on Re. What is genuinely disjoint is alpha and
// e/D, confirmed here so a future edit cannot silently narrow either band
// into the other's territory without this test noticing.
TEST(RibCorrelationTest, TheTwoSetsDoNotOverlapInAlphaOrEd) {
  const auto han = han_1988_orthogonal();
  const auto rallabandi = rallabandi_2009_high_re();

  EXPECT_EQ(han.valid_alpha.lo, 90.0);
  EXPECT_EQ(han.valid_alpha.hi, 90.0);
  EXPECT_EQ(rallabandi.valid_alpha.lo, 45.0);
  EXPECT_EQ(rallabandi.valid_alpha.hi, 45.0);

  // Rallabandi's study was a square channel only.
  EXPECT_EQ(rallabandi.valid_WH.lo, 1.0);
  EXPECT_EQ(rallabandi.valid_WH.hi, 1.0);

  // e/D bands do not overlap: Han's is small ribs, Rallabandi's is large.
  EXPECT_LT(han.valid_eD.hi, rallabandi.valid_eD.lo);

  // Their Re bands DO overlap -- recorded as a fact, not a defect. See the
  // comment above: this is why the two sets are named and chosen, not
  // switched on Reynolds number.
  EXPECT_GT(han.valid_Re.hi, rallabandi.valid_Re.lo);
}

// Reversing the flow must not change the magnitude, same reasoning as the
// 90 deg set: a 45 deg parallel rib reversed is its mirror image.
TEST(RibCorrelationTest, RallabandiReverseFlowIsSymmetricInMagnitude) {
  const auto set = rallabandi_2009_high_re();
  const auto g = rallabandi_ref_geometry();
  const auto fwd = evaluate_rib(set, g, 100000.0);
  const auto rev = evaluate_rib(set, g, -100000.0);

  EXPECT_NEAR(fwd.G, rev.G, 1e-9);
  EXPECT_NEAR(fwd.St_r, rev.St_r, 1e-9);
  EXPECT_NEAR(fwd.e_plus, -rev.e_plus, 1e-6);
}

// Same guard discipline as the 90 deg set: every state a solver can probe
// must return finite, positive numbers, never NaN or a complex result.
TEST(RibCorrelationTest, RallabandiGuardsHoldForStatesTheSolverActuallyProbes) {
  const auto set = rallabandi_2009_high_re();
  const auto g = rallabandi_ref_geometry();
  for (double Re : {-1e6, -1e4, -1.0, 0.0, 1e-12, 1.0, 1e8}) {
    const auto r = evaluate_rib(set, g, Re);
    EXPECT_TRUE(std::isfinite(r.f)) << "Re = " << Re;
    EXPECT_TRUE(std::isfinite(r.G)) << "Re = " << Re;
    EXPECT_TRUE(std::isfinite(r.St_r)) << "Re = " << Re;
    EXPECT_TRUE(std::isfinite(r.dSt_dRe)) << "Re = " << Re;
    EXPECT_GT(r.f, 0.0) << "Re = " << Re;
    EXPECT_GT(r.St_r, 0.0) << "Re = " << Re;
  }
}

// The analytic Stanton derivative must match a central difference through
// this set's own guards too -- not assumed to inherit correctness from the
// shared evaluate_rib code path just because han_1988_orthogonal's does.
TEST(RibCorrelationTest, RallabandiStantonDerivativeMatchesCentralDifferences) {
  const auto set = rallabandi_2009_high_re();
  const auto g = rallabandi_ref_geometry();
  for (double Re : {200000.0, 50000.0, 100.0, 10.0}) {
    const double h = std::max(1e-6, std::abs(Re) * 1e-6);
    const double fd = (evaluate_rib(set, g, Re + h).St_r -
                       evaluate_rib(set, g, Re - h).St_r) /
                      (2.0 * h);
    const double analytic = evaluate_rib(set, g, Re).dSt_dRe;
    const double scale = std::max({std::abs(fd), std::abs(analytic), 1e-14});
    EXPECT_LT(std::abs(analytic - fd) / scale, 1e-5) << "Re = " << Re;
  }
}

TEST(RibCorrelationTest, RallabandiSetPassesValidation) {
  EXPECT_NO_THROW(validate_rib_set(rallabandi_2009_high_re()));
}


// Hand-computed independently (Python, from Eq. 4.17/4.18 as printed), not
// from this implementation -- see
// validation/cooling/extractions/han_ribbed.md. e+ agrees only after
// applying the documented smooth floor (sqrt(e+^2 + 1)); a raw-e+
// computation differs at the 5th decimal in G by construction, not error.
TEST(RibCorrelationTest, ReproducesTheConfirmedHanParkExtraction) {
  const auto set = han_park_1988_angled();

  RibGeometry square90;
  square90.e_D = 0.0625;
  square90.p_e = 15.0;
  square90.W_H = 1.0;
  square90.alpha_deg = 90.0;
  const auto r1 = evaluate_rib(set, square90, 30000.0);
  EXPECT_NEAR(r1.R, 3.5726760161, 1e-6);
  EXPECT_NEAR(r1.f, 0.0508531225, 1e-6);
  EXPECT_NEAR(r1.e_plus, 298.9820302819, 1e-3);
  EXPECT_NEAR(r1.G, 17.1527507, 1e-4);  // matches only with the e+ floor
  EXPECT_NEAR(r1.St_r, 0.0080325541, 1e-6);

  RibGeometry angled;
  angled.e_D = 0.0625;
  angled.p_e = 15.0;
  angled.W_H = 2.0;
  angled.alpha_deg = 45.0;
  const auto r2 = evaluate_rib(set, angled, 30000.0);
  EXPECT_NEAR(r2.R, 4.7592382829, 1e-6);
  EXPECT_NEAR(r2.f, 0.0440439106, 1e-6);

  // W/H = 3, above Eq. 4.17's own cap: the R quadratic must use W/H = 2, not
  // 3, in its (W/H)^m term. If the cap were dropped this would not match.
  RibGeometry capped;
  capped.e_D = 0.06;
  capped.p_e = 12.0;
  capped.W_H = 3.0;
  capped.alpha_deg = 60.0;
  const auto r3 = evaluate_rib(set, capped, 40000.0);
  EXPECT_NEAR(r3.R, 2.9903078484, 1e-6);
  EXPECT_NEAR(r3.f, 0.0876323069, 1e-6);
}

// Eq. 4.17's own cap ("if W/H > 2, set W/H = 2") must actually engage: R at
// W/H = 3 and W/H = 5 should be identical to R at the capped value 2, not
// merely close, since above the cap the (W/H)^m term stops depending on
// W/H at all.
TEST(RibCorrelationTest, HanParkWHCapActuallyCaps) {
  const auto set = han_park_1988_angled();
  RibGeometry g;
  g.e_D = 0.06;
  g.p_e = 12.0;
  g.alpha_deg = 60.0;  // off 90 deg, so the (W/H)^m term is active (m=0.35)

  g.W_H = 2.0;
  const double R_at_cap = evaluate_rib(set, g, 40000.0).R;
  g.W_H = 3.0;
  const double R_above = evaluate_rib(set, g, 40000.0).R;
  g.W_H = 5.0;
  const double R_further_above = evaluate_rib(set, g, 40000.0).R;

  EXPECT_NEAR(R_above, R_at_cap, 1e-9);
  EXPECT_NEAR(R_further_above, R_at_cap, 1e-9);

  // Below the cap, R must actually vary with W/H -- otherwise the cap logic
  // could be silently always-on and this test would not distinguish it.
  g.W_H = 1.0;
  const double R_below = evaluate_rib(set, g, 40000.0).R;
  EXPECT_GT(std::abs(R_below - R_at_cap), 1e-6);
}

// The R and G switches are genuine discontinuities the source states, not
// numerical artefacts to be smoothed away -- pin the jump sizes so a future
// "helpful" smoothing change cannot land silently.
TEST(RibCorrelationTest, HanParkRAndGSwitchesAreGenuineDiscontinuities) {
  const auto set = han_park_1988_angled();
  RibGeometry g;
  g.e_D = 0.06;
  g.p_e = 12.0;
  g.W_H = 2.0;

  g.alpha_deg = 90.0;
  const double R_at_90 = evaluate_rib(set, g, 40000.0).R;
  g.alpha_deg = 90.0 - 1e-6;
  const double R_just_below_90 = evaluate_rib(set, g, 40000.0).R;
  // A ~28% jump at W/H=2 (2^0.35 - 1), not a rounding-sized difference.
  EXPECT_GT(std::abs(R_just_below_90 / R_at_90 - 1.0), 0.25);

  RibGeometry g2;
  g2.e_D = 0.06;
  g2.p_e = 20.0;
  g2.alpha_deg = 30.0;

  g2.W_H = 1.0;
  const double G_at_square = evaluate_rib(set, g2, 40000.0).G;
  g2.W_H = 1.0 - 1e-6;
  const double G_just_off_square = evaluate_rib(set, g2, 40000.0).G;
  EXPECT_GT(std::abs(G_just_off_square / G_at_square - 1.0), 0.20);
}

// Same reverse-flow symmetry as the other two sets: a parallel angled rib
// reversed is its mirror image.
TEST(RibCorrelationTest, HanParkReverseFlowIsSymmetricInMagnitude) {
  const auto set = han_park_1988_angled();
  RibGeometry g;
  g.e_D = 0.0625;
  g.p_e = 15.0;
  g.W_H = 2.0;
  g.alpha_deg = 45.0;
  const auto fwd = evaluate_rib(set, g, 30000.0);
  const auto rev = evaluate_rib(set, g, -30000.0);

  EXPECT_NEAR(fwd.G, rev.G, 1e-9);
  EXPECT_NEAR(fwd.St_r, rev.St_r, 1e-9);
  EXPECT_NEAR(fwd.e_plus, -rev.e_plus, 1e-6);
}

// Same guard discipline: every state a solver can probe must return finite,
// positive numbers, across the full alpha and W/H range including the two
// switch boundaries themselves.
TEST(RibCorrelationTest, HanParkGuardsHoldForStatesTheSolverActuallyProbes) {
  const auto set = han_park_1988_angled();
  for (double alpha : {30.0, 45.0, 60.0, 89.999, 90.0, 90.001}) {
    for (double W_H : {0.9999, 1.0, 1.0001, 2.0, 4.0, 10.0}) {
      RibGeometry g;
      g.e_D = 0.0625;
      g.p_e = 15.0;
      g.W_H = W_H;
      g.alpha_deg = alpha;
      for (double Re : {-1e6, -1.0, 0.0, 1.0, 1e8}) {
        const auto r = evaluate_rib(set, g, Re);
        EXPECT_TRUE(std::isfinite(r.f))
            << "alpha=" << alpha << " W_H=" << W_H << " Re=" << Re;
        EXPECT_TRUE(std::isfinite(r.G))
            << "alpha=" << alpha << " W_H=" << W_H << " Re=" << Re;
        EXPECT_TRUE(std::isfinite(r.St_r))
            << "alpha=" << alpha << " W_H=" << W_H << " Re=" << Re;
        EXPECT_GT(r.f, 0.0) << "alpha=" << alpha << " W_H=" << W_H;
        EXPECT_GT(r.St_r, 0.0) << "alpha=" << alpha << " W_H=" << W_H;
      }
    }
  }
}

// The analytic Stanton derivative must match a central difference through
// this set's guards too. Perturbing alpha by h risks crossing the alpha=90
// switch itself if alpha is near it, which would corrupt a naive central
// difference through no fault of the Jacobian -- so this sweeps Re (which
// the switch does not depend on) at a fixed, safely-off-switch geometry.
TEST(RibCorrelationTest, HanParkStantonDerivativeMatchesCentralDifferences) {
  const auto set = han_park_1988_angled();
  RibGeometry g;
  g.e_D = 0.0625;
  g.p_e = 15.0;
  g.W_H = 2.0;
  g.alpha_deg = 45.0;
  for (double Re : {50000.0, 20000.0, 100.0, 10.0}) {
    const double h = std::max(1e-6, std::abs(Re) * 1e-6);
    const double fd = (evaluate_rib(set, g, Re + h).St_r -
                       evaluate_rib(set, g, Re - h).St_r) /
                      (2.0 * h);
    const double analytic = evaluate_rib(set, g, Re).dSt_dRe;
    const double scale = std::max({std::abs(fd), std::abs(analytic), 1e-14});
    EXPECT_LT(std::abs(analytic - fd) / scale, 1e-5) << "Re = " << Re;
  }
}

TEST(RibCorrelationTest, HanParkSetPassesValidation) {
  EXPECT_NO_THROW(validate_rib_set(han_park_1988_angled()));

  // C_R plays no role in the quadratic-alpha shape, so it must NOT be
  // rejected for being unset (0.0) the way the plain power-law shape is.
  auto s = han_park_1988_angled();
  EXPECT_EQ(s.C_R, 0.0);
  EXPECT_NO_THROW(validate_rib_set(s));

  auto bad_cap = han_park_1988_angled();
  bad_cap.R_quad_WH_cap = -1.0;
  EXPECT_THROW(validate_rib_set(bad_cap), std::invalid_argument);
}

// Han (1988) (90 deg only) and Han & Park (1988) (angled, this set) both
// claim alpha = 90 deg, and their R correlations were confirmed to agree
// closely there (han_ribbed.md item 23: 3.44%) while their G correlations
// disagree by up to 22% -- different papers, not obliged to agree, and the
// asymmetry between the two functions is itself part of the confirmed
// record. Pinned so a future edit cannot silently "fix" G to match R's
// agreement, which would misrepresent what the sources actually say.
TEST(RibCorrelationTest, HanParkAgreesWithHan1988OnRButNotG) {
  const auto han = han_1988_orthogonal();
  const auto park = han_park_1988_angled();

  RibGeometry g;
  g.e_D = 0.0625;
  g.p_e = 10.0;
  g.W_H = 1.0;
  g.alpha_deg = 90.0;

  const auto r_han = evaluate_rib(han, g, 30000.0);
  const auto r_park = evaluate_rib(park, g, 30000.0);
  const double R_ratio = r_park.R / r_han.R;
  EXPECT_NEAR(R_ratio, 1.0, 0.06);  // within ~6%, close per item 23

  const double G_ratio = r_park.G / r_han.G;
  EXPECT_GT(std::abs(G_ratio - 1.0), 0.05);  // genuinely disagree
}

// Eq. 4.19 (Han et al. 1989), narrow-aspect-ratio channels, W/H < 1.
// R has no printed equation (han_ribbed.md item 40) -- the coefficients here
// are this project's own fit to Fig. 4.48a's drawn line, hand-checked
// against the raw quadratic rather than a stored magic number.
TEST(RibCorrelationTest, ReproducesTheConfirmedHan1989NarrowExtraction) {
  const auto set = han_1989_narrow_channel();

  // Wide sub-band (1/2 <= W/H < 1), alpha = 90 deg: R from the wide
  // quadratic, G = 2.24 (e+)^0.35 exactly (text-extracted, item 37).
  RibGeometry wide;
  wide.e_D = 0.0625;
  wide.p_e = 15.0;
  wide.W_H = 0.5;
  wide.alpha_deg = 90.0;
  const auto r1 = evaluate_rib(set, wide, 30000.0);
  const double u1 = 1.0;  // alpha/90
  const double expected_R1 = 10.66406939 + -21.43911537 * u1 +
                             15.49713879 * u1 * u1;
  EXPECT_NEAR(r1.R, expected_R1, 1e-6);
  EXPECT_NEAR(r1.G, 2.24 * std::pow(r1.e_plus, 0.35), 1e-4);  // e+ floor

  // Narrow sub-band (1/4 < W/H < 1/2), alpha = 45 deg (off-axis, C = 1.80):
  // both C and n carry the (W/H)^-0.76 / (W/H)^0.44 correction.
  RibGeometry narrow;
  narrow.e_D = 0.0625;
  narrow.p_e = 15.0;
  narrow.W_H = 0.25;
  narrow.alpha_deg = 45.0;
  const auto r2 = evaluate_rib(set, narrow, 30000.0);
  const double u2 = 0.5;  // alpha/90
  const double expected_R2 = 8.79873413 + -14.76786182 * u2 +
                             12.4724873 * u2 * u2;
  EXPECT_NEAR(r2.R, expected_R2, 1e-6);
  const double expected_C2 = 1.80 * std::pow(0.25, -0.76);
  const double expected_n2 = 0.35 * std::pow(0.25, 0.44);
  EXPECT_NEAR(r2.G, expected_C2 * std::pow(r2.e_plus, expected_n2),
              1e-4);  // e+ floor
}

// Both switches this set introduces are genuine discontinuities the source
// states, not numerical artefacts -- pin the jump sizes the same way
// HanParkRAndGSwitchesAreGenuineDiscontinuities does for the other set.
TEST(RibCorrelationTest, HanNarrowRAndGSwitchesAreGenuineDiscontinuities) {
  const auto set = han_1989_narrow_channel();

  // R: the two sub-bands are independently fitted quadratics, so there is
  // no reason for them to meet at the boundary -- confirm they do not.
  RibGeometry g;
  g.e_D = 0.0625;
  g.p_e = 15.0;
  g.alpha_deg = 60.0;
  g.W_H = 0.5;  // wide branch (>= boundary)
  const double R_wide_side = evaluate_rib(set, g, 30000.0).R;
  g.W_H = 0.5 - 1e-6;  // narrow branch
  const double R_narrow_side = evaluate_rib(set, g, 30000.0).R;
  EXPECT_GT(std::abs(R_narrow_side / R_wide_side - 1.0), 0.01);

  // G: a ~20% jump in C at alpha == 90 deg (1.80/2.24 - 1 = -19.6%, for any
  // W/H), not a rounding-sized difference.
  RibGeometry g2;
  g2.e_D = 0.06;
  g2.p_e = 12.0;
  g2.W_H = 0.5;
  g2.alpha_deg = 90.0;
  const double G_at_90 = evaluate_rib(set, g2, 40000.0).G;
  g2.alpha_deg = 90.0 - 1e-6;
  const double G_just_below_90 = evaluate_rib(set, g2, 40000.0).G;
  EXPECT_GT(std::abs(G_just_below_90 / G_at_90 - 1.0), 0.15);
}

// Same reverse-flow symmetry as the other three sets: parallel angled ribs
// reversed are the mirror image.
TEST(RibCorrelationTest, HanNarrowReverseFlowIsSymmetricInMagnitude) {
  const auto set = han_1989_narrow_channel();
  RibGeometry g;
  g.e_D = 0.0625;
  g.p_e = 15.0;
  g.W_H = 0.3;
  g.alpha_deg = 45.0;
  const auto fwd = evaluate_rib(set, g, 30000.0);
  const auto rev = evaluate_rib(set, g, -30000.0);

  EXPECT_NEAR(fwd.G, rev.G, 1e-9);
  EXPECT_NEAR(fwd.St_r, rev.St_r, 1e-9);
  EXPECT_NEAR(fwd.e_plus, -rev.e_plus, 1e-6);
}

// Every state a solver can probe must return finite, positive numbers,
// across both switch boundaries.
TEST(RibCorrelationTest, HanNarrowGuardsHoldForStatesTheSolverActuallyProbes) {
  const auto set = han_1989_narrow_channel();
  for (double alpha : {30.0, 45.0, 89.999, 90.0, 90.001}) {
    for (double W_H : {0.25, 0.4999, 0.5, 0.5001, 0.75, 0.9999}) {
      RibGeometry g;
      g.e_D = 0.0625;
      g.p_e = 15.0;
      g.W_H = W_H;
      g.alpha_deg = alpha;
      for (double Re : {-1e6, -1.0, 0.0, 1.0, 1e8}) {
        const auto r = evaluate_rib(set, g, Re);
        EXPECT_TRUE(std::isfinite(r.f))
            << "alpha=" << alpha << " W_H=" << W_H << " Re=" << Re;
        EXPECT_TRUE(std::isfinite(r.G))
            << "alpha=" << alpha << " W_H=" << W_H << " Re=" << Re;
        EXPECT_TRUE(std::isfinite(r.St_r))
            << "alpha=" << alpha << " W_H=" << W_H << " Re=" << Re;
        EXPECT_GT(r.f, 0.0) << "alpha=" << alpha << " W_H=" << W_H;
        EXPECT_GT(r.St_r, 0.0) << "alpha=" << alpha << " W_H=" << W_H;
      }
    }
  }
}

TEST(RibCorrelationTest, HanNarrowStantonDerivativeMatchesCentralDifferences) {
  const auto set = han_1989_narrow_channel();
  RibGeometry g;
  g.e_D = 0.0625;
  g.p_e = 15.0;
  g.W_H = 0.3;
  g.alpha_deg = 45.0;
  for (double Re : {50000.0, 20000.0, 100.0, 10.0}) {
    const double h = std::max(1e-6, std::abs(Re) * 1e-6);
    const double fd = (evaluate_rib(set, g, Re + h).St_r -
                       evaluate_rib(set, g, Re - h).St_r) /
                      (2.0 * h);
    const double analytic = evaluate_rib(set, g, Re).dSt_dRe;
    const double scale = std::max({std::abs(fd), std::abs(analytic), 1e-14});
    EXPECT_LT(std::abs(analytic - fd) / scale, 1e-5) << "Re = " << Re;
  }
}

TEST(RibCorrelationTest, HanNarrowSetPassesValidation) {
  EXPECT_NO_THROW(validate_rib_set(han_1989_narrow_channel()));

  // C_R and C_G play no role in this set's shapes, so they must NOT be
  // rejected for being unset (0.0) the way the plain power-law shapes are.
  auto s = han_1989_narrow_channel();
  EXPECT_EQ(s.C_R, 0.0);
  EXPECT_EQ(s.C_G, 0.0);
  EXPECT_NO_THROW(validate_rib_set(s));

  auto bad_boundary = han_1989_narrow_channel();
  bad_boundary.R_WH_band_boundary = -1.0;
  EXPECT_THROW(validate_rib_set(bad_boundary), std::invalid_argument);

  auto bad_g_boundary = han_1989_narrow_channel();
  bad_g_boundary.G_narrow_WH_band_boundary = 0.0;
  EXPECT_THROW(validate_rib_set(bad_g_boundary), std::invalid_argument);
}

// Adding this set must not perturb any of the other three -- new enum
// values and fields are additive, defaulted to the prior behaviour. Their
// own extraction tests (ReproducesTheConfirmedExtraction,
// ReproducesTheConfirmedRallabandiExtraction,
// ReproducesTheConfirmedHanParkExtraction) running unmodified alongside
// this one, in the same binary, is the proof; this test only pins the
// cheap, explicit version of that claim.
TEST(RibCorrelationTest, HanNarrowChannelDoesNotPerturbTheOtherThreeSets) {
  RibGeometry g;
  g.e_D = 0.0625;
  g.p_e = 15.0;
  g.W_H = 1.0;
  g.alpha_deg = 90.0;

  const auto han = evaluate_rib(han_1988_orthogonal(), g, 30000.0);
  EXPECT_NEAR(han.R, 3.2 * std::pow(15.0 / 10.0, 0.35), 1e-9);

  g.W_H = 2.0;
  g.alpha_deg = 45.0;
  const auto park = evaluate_rib(han_park_1988_angled(), g, 30000.0);
  EXPECT_NEAR(park.R, 4.7592382829, 1e-6);

  g.W_H = 1.0;
  g.alpha_deg = 45.0;
  const auto ra = evaluate_rib(rallabandi_2009_high_re(), g, 100000.0);
  EXPECT_TRUE(std::isfinite(ra.G));
}

// The normaliser is data, not a convention: the same correlation written
// against a raw p/e needs a different constant, and mixing them is wrong by a
// constant factor that never looks like a trend.
TEST(RibCorrelationTest, NormaliserAndConstantAreEquivalentWhenConsistent) {
  auto normalised = han_1988_orthogonal();  // 3.2 * (p/e / 10)^0.35

  auto raw = han_1988_orthogonal();
  raw.C_R = 3.2 * std::pow(10.0, -0.35);  // 1.4294...
  raw.R_pe.reference = 1.0;

  auto g = ref_geometry();
  for (double p_e : {5.0, 10.0, 20.0}) {
    g.p_e = p_e;
    EXPECT_NEAR(evaluate_rib(normalised, g, 30000.0).R,
                evaluate_rib(raw, g, 30000.0).R, 1e-12)
        << "p/e = " << p_e;
  }

  // And the mistake this prevents: the normalised constant used raw.
  auto mixed = han_1988_orthogonal();
  mixed.R_pe.reference = 1.0;  // constant left at 3.2
  g.p_e = 10.0;
  EXPECT_NEAR(evaluate_rib(mixed, g, 30000.0).R /
                  evaluate_rib(normalised, g, 30000.0).R,
              std::pow(10.0, 0.35), 1e-9);
}

// --------------------------------------------------------------------------
// Accuracy provenance (#389). A band measured FROM a model must never be the
// band that model is judged against; the type makes the distinction, and
// validate_rib_set makes it impossible to record half of it.
// --------------------------------------------------------------------------

TEST(AccuracyProvenance, OnlyAnAuthorStatedFigureIsUsableAsABand) {
  // Han's own "95% of data within X%" claims -- the only stated pair among
  // the shipped sets.
  const auto han = combaero::cooling::han_1988_orthogonal();
  EXPECT_EQ(han.accuracy_R.provenance,
            combaero::cooling::AccuracyProvenance::Stated);
  EXPECT_TRUE(han.accuracy_R.usable_as_band());
  EXPECT_TRUE(han.accuracy_G.usable_as_band());
  EXPECT_DOUBLE_EQ(han.accuracy_R.value, 0.06);

  // Measured THROUGH evaluate_rib -- the code path it would be judging.
  const auto hp = combaero::cooling::han_park_1988_angled();
  EXPECT_EQ(hp.accuracy_R.provenance,
            combaero::cooling::AccuracyProvenance::Measured);
  EXPECT_FALSE(hp.accuracy_R.usable_as_band());
  EXPECT_FALSE(hp.accuracy_G.usable_as_band());

  // Unstated carries NaN, not 0.0. The old sentinel read as "perfect".
  const auto ra = combaero::cooling::rallabandi_2009_high_re();
  EXPECT_EQ(ra.accuracy_R.provenance,
            combaero::cooling::AccuracyProvenance::Unstated);
  EXPECT_TRUE(std::isnan(ra.accuracy_R.value));
  EXPECT_FALSE(ra.accuracy_R.usable_as_band());
  // A NaN fails every comparison, so an accidental read cannot pass a
  // tolerance check silently -- which 0.0 would have done.
  EXPECT_FALSE(ra.accuracy_R.value < 0.10);
  EXPECT_FALSE(ra.accuracy_R.value >= 0.10);
  EXPECT_EQ(ra.accuracy_G.provenance,
            combaero::cooling::AccuracyProvenance::Measured);
}

TEST(AccuracyProvenance, ValueAndProvenanceCannotDriftApart) {
  auto s = combaero::cooling::han_1988_orthogonal();
  EXPECT_NO_THROW(combaero::cooling::validate_rib_set(s));

  // A value with no provenance behind it.
  s.accuracy_R = combaero::cooling::StatedAccuracy::unstated();
  s.accuracy_R.value = 0.06;
  EXPECT_THROW(combaero::cooling::validate_rib_set(s), std::invalid_argument);

  // A provenance with no value behind it.
  auto t = combaero::cooling::han_1988_orthogonal();
  t.accuracy_G = combaero::cooling::StatedAccuracy::stated(
      std::numeric_limits<double>::quiet_NaN());
  EXPECT_THROW(combaero::cooling::validate_rib_set(t), std::invalid_argument);

  // A negative band is not a band.
  auto u = combaero::cooling::han_1988_orthogonal();
  u.accuracy_G = combaero::cooling::StatedAccuracy::stated(-0.05);
  EXPECT_THROW(combaero::cooling::validate_rib_set(u), std::invalid_argument);
}

// Han, Zhang and Lee (1991) JHT 113, 590, Table 2 -- every row, against the
// printed coefficients (four-channel verified). G and G_bar are a (e+)^b at
// the chain's own e+; the smooth e+ floor is negligible at these e+.
TEST(RibCorrelationTest, HanZhangLee1991ReproducesTable2) {
  struct Row {
    RibShape shape;
    double alpha, R_a, G_a, G_b, Gbar_a, Gbar_b;
    bool symmetric;
  };
  const Row rows[] = {
      {RibShape::Transverse, 90, 3.18, 3.97, 0.28, 4.86, 0.28, true},
      {RibShape::Parallel, 60, 2.05, 1.52, 0.41, 2.41, 0.36, true},
      {RibShape::Crossed, 60, 3.18, 3.24, 0.32, 4.57, 0.28, true},
      {RibShape::V, 60, 1.72, 1.35, 0.42, 1.76, 0.40, false},
      {RibShape::Lambda, 60, 1.42, 1.59, 0.43, 2.12, 0.41, false},
      {RibShape::Parallel, 45, 3.05, 2.07, 0.36, 3.01, 0.32, true},
      {RibShape::Crossed, 45, 4.40, 2.11, 0.37, 2.43, 0.37, true},
      {RibShape::V, 45, 2.04, 1.36, 0.43, 1.93, 0.40, false},
      {RibShape::Lambda, 45, 1.70, 1.83, 0.41, 2.49, 0.38, false},
  };
  for (const auto &row : rows) {
    const auto s = han_zhang_lee_1991(row.shape, row.alpha);
    EXPECT_NO_THROW(validate_rib_set(s)) << s.name;
    EXPECT_EQ(s.shape, row.shape) << s.name;
    EXPECT_EQ(s.symmetric, row.symmetric) << s.name;
    RibGeometry g;
    g.e_D = 0.0625;
    g.p_e = 10.0;
    g.W_H = 1.0;
    g.alpha_deg = row.alpha;
    for (double Re : {15000.0, 40000.0, 90000.0}) {
      const auto r = evaluate_rib(s, g, Re);
      EXPECT_NEAR(r.R, row.R_a, 1e-12) << s.name;
      EXPECT_NEAR(r.G, row.G_a * std::pow(r.e_plus, row.G_b), 1e-4 * r.G) << s.name;
      ASSERT_TRUE(r.has_G_bar) << s.name;
      EXPECT_NEAR(r.G_bar, row.Gbar_a * std::pow(r.e_plus, row.Gbar_b), 1e-4 * r.G_bar)
          << s.name;
      EXPECT_GT(r.G_bar, r.G) << s.name << ": printed G_bar/G is 1.07-1.26";
      EXPECT_FALSE(r.extrapolated) << s.name << " at its own rig, Re " << Re;
    }
  }
}

// The paper tested nine configurations; a nearby one is not a substitute.
TEST(RibCorrelationTest, HanZhangLee1991RefusesUntestedConfigurations) {
  EXPECT_THROW(han_zhang_lee_1991(RibShape::V, 30.0), std::invalid_argument);
  EXPECT_THROW(han_zhang_lee_1991(RibShape::Transverse, 60.0), std::invalid_argument);
  EXPECT_THROW(han_zhang_lee_1991(RibShape::Parallel, 90.0), std::invalid_argument);
  EXPECT_THROW(han_zhang_lee_1991(RibShape::Unspecified, 90.0), std::invalid_argument);
}

TEST(RibCorrelationTest, HanZhangLee1991GuardsAndDerivative) {
  for (auto shape : {RibShape::Crossed, RibShape::V, RibShape::Lambda}) {
    const auto s = han_zhang_lee_1991(shape, 60.0);
    RibGeometry g;
    g.e_D = 0.0625;
    g.p_e = 10.0;
    g.W_H = 1.0;
    g.alpha_deg = 60.0;
    for (double Re : {-1e6, -1.0, 0.0, 1.0, 1e8}) {
      const auto r = evaluate_rib(s, g, Re);
      EXPECT_TRUE(std::isfinite(r.G) && std::isfinite(r.G_bar) && std::isfinite(r.St_r))
          << s.name << " Re=" << Re;
      EXPECT_GT(r.St_r, 0.0) << s.name << " Re=" << Re;
    }
    for (double Re : {50000.0, 20000.0, 100.0}) {
      const double h = std::max(1e-6, std::abs(Re) * 1e-6);
      const double fd =
          (evaluate_rib(s, g, Re + h).St_r - evaluate_rib(s, g, Re - h).St_r) / (2.0 * h);
      const double an = evaluate_rib(s, g, Re).dSt_dRe;
      EXPECT_LT(std::abs(an - fd) / std::max({std::abs(fd), std::abs(an), 1e-14}), 1e-5)
          << s.name << " Re=" << Re;
    }
  }
}

// Adding shape and printed G_bar must leave the four earlier sets as they
// were: each declares its shape, prints no G_bar, and validates.
TEST(RibCorrelationTest, EarlierSetsDeclareShapeAndPrintNoGbar) {
  const std::pair<RibCorrelationSet, RibShape> sets[] = {
      {han_1988_orthogonal(), RibShape::Transverse},
      {han_park_1988_angled(), RibShape::Parallel},
      {rallabandi_2009_high_re(), RibShape::Parallel},
      {han_1989_narrow_channel(), RibShape::Parallel},
  };
  for (const auto &[s, shape] : sets) {
    EXPECT_EQ(s.shape, shape) << s.name;
    EXPECT_EQ(s.C_Gbar, 0.0) << s.name;
    EXPECT_FALSE(evaluate_rib(s, ref_geometry(), 30000.0).has_G_bar) << s.name;
  }
  EXPECT_EQ(RibCorrelationSet{}.shape, RibShape::Unspecified);

  auto bad = han_zhang_lee_1991(RibShape::V, 60.0);
  bad.C_Gbar = -1.0;
  EXPECT_THROW(validate_rib_set(bad), std::invalid_argument);
}

// #1106: the channel's Darcy friction is Han's area weighting (Han and Park
// 1988 Eq. 6) run forwards, so two ribbed walls give back his measured
// two-wall fbar, x4 for Darcy.
TEST(RibCorrelationTest, ChannelFrictionIsHansAreaWeightingInDarcyForm) {
  const auto s = han_1988_orthogonal();
  for (const double w_h : {0.5, 1.0, 2.0, 4.0}) {
    RibGeometry g = ref_geometry();
    g.W_H = w_h;
    const double re = 30000.0;
    const double f4 = evaluate_rib(s, g, re).f;
    const double f_s = 0.046 * std::pow(re, -0.2);
    const double fbar = (f4 * w_h + f_s) / (w_h + 1.0);
    // Eq. 6 holds for the fbar recovered.
    EXPECT_NEAR(f4, fbar + (fbar - f_s) / w_h, 1e-12);
    const auto [f2, df2] = rib_channel_friction(s, g, re, 2);
    EXPECT_NEAR(f2, 4.0 * fbar, 1e-9 * f2);
    const auto [f_four, df_four] = rib_channel_friction(s, g, re, 4);
    EXPECT_DOUBLE_EQ(f_four, 4.0 * f4);
    EXPECT_EQ(df_four, 0.0);
    const auto [f1, df1] = rib_channel_friction(s, g, re, 1);
    EXPECT_LT(f1, f2);
    EXPECT_LT(df1, 0.0);  // smooth walls' f_s falls with Re
    EXPECT_LT(df1, df2);  // more smooth wall, more Re dependence
  }
}

TEST(RibCorrelationTest, ChannelFrictionDerivativeIsFiniteAndOddThroughZero) {
  const auto s = han_1988_orthogonal();
  for (const double re : {-4.0e4, -1.0, 0.0, 1.0, 2.0e4}) {
    const auto [f, df] = rib_channel_friction(s, ref_geometry(), re, 1);
    const double h = std::max(1e-6 * std::abs(re), 1e-4);
    const double fd = (std::get<0>(rib_channel_friction(s, ref_geometry(), re + h, 1)) -
                       std::get<0>(rib_channel_friction(s, ref_geometry(), re - h, 1))) /
                      (2.0 * h);
    EXPECT_TRUE(std::isfinite(f));
    EXPECT_NEAR(df, fd, 1e-6 * std::max(std::abs(fd), 1e-9));
    EXPECT_DOUBLE_EQ(std::get<1>(rib_channel_friction(s, ref_geometry(), -re, 1)), -df);
  }
}

TEST(RibCorrelationTest, RibbedWallFractionAndRefusals) {
  EXPECT_DOUBLE_EQ(rib_ribbed_wall_fraction(1.0, 2), 0.5);
  EXPECT_DOUBLE_EQ(rib_ribbed_wall_fraction(2.0, 2), 2.0 / 3.0);
  EXPECT_DOUBLE_EQ(rib_ribbed_wall_fraction(2.0, 1), 1.0 / 3.0);
  EXPECT_DOUBLE_EQ(rib_ribbed_wall_fraction(0.25, 4), 1.0);
  EXPECT_THROW(rib_ribbed_wall_fraction(1.0, 3), std::invalid_argument);
  EXPECT_THROW(rib_ribbed_wall_fraction(1.0, 0), std::invalid_argument);
  EXPECT_THROW(rib_ribbed_wall_fraction(0.0, 2), std::invalid_argument);
}
