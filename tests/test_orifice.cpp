#include <gtest/gtest.h>
#include "../include/orifice.h"
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

// Test fixture for orifice tests
class OrificeTest : public ::testing::Test {
protected:
    OrificeGeometry geom;
    OrificeState state;

    void SetUp() override {
        // Standard test case: 50mm orifice in 100mm pipe
        geom.d = 0.050;  // 50 mm
        geom.D = 0.100;  // 100 mm

        // Typical flow conditions
        state.Re_D = 100000.0;
        state.dP = 10000.0;   // 10 kPa
        state.rho = 1.2;      // kg/m^3 (air at ~20degC)
        state.mu = 1.8e-5;    // Pa*s
    }
};

// -------------------------------------------------------------
// Geometry tests
// -------------------------------------------------------------

TEST_F(OrificeTest, GeometryBeta) {
    EXPECT_DOUBLE_EQ(geom.beta(), 0.5);
}

TEST_F(OrificeTest, GeometryArea) {
    const double expected = 3.14159265358979323846 * 0.050 * 0.050 / 4.0;
    EXPECT_NEAR(geom.area(), expected, 1e-10);
}

TEST_F(OrificeTest, GeometryValid) {
    EXPECT_TRUE(geom.is_valid());

    OrificeGeometry invalid;
    invalid.d = 0.0;
    invalid.D = 0.1;
    EXPECT_FALSE(invalid.is_valid());

    invalid.d = 0.15;  // d > D
    invalid.D = 0.1;
    EXPECT_FALSE(invalid.is_valid());
}

// -------------------------------------------------------------
// Sharp thin-plate Cd tests
// -------------------------------------------------------------

TEST_F(OrificeTest, CdSharpThinPlateRange) {
    // Cd should be in reasonable range for sharp-edged orifice
    double Cd = Cd_sharp_thin_plate(geom, state);
    EXPECT_GT(Cd, 0.58);
    EXPECT_LT(Cd, 0.68);
}

TEST_F(OrificeTest, CdSharpThinPlateBetaDependence) {
    // Cd should increase slightly with beta
    OrificeGeometry geom_low, geom_high;
    geom_low.d = 0.020;
    geom_low.D = 0.100;  // beta = 0.2
    geom_high.d = 0.070;
    geom_high.D = 0.100; // beta = 0.7

    double Cd_low = Cd_sharp_thin_plate(geom_low, state);
    double Cd_high = Cd_sharp_thin_plate(geom_high, state);

    // Both should be in valid range
    EXPECT_GT(Cd_low, 0.58);
    EXPECT_LT(Cd_low, 0.68);
    EXPECT_GT(Cd_high, 0.58);
    EXPECT_LT(Cd_high, 0.68);
}

TEST_F(OrificeTest, CdSharpThinPlateReynoldsDependence) {
    // Cd should approach asymptotic value at high Re
    OrificeState state_low, state_high;
    state_low.Re_D = 10000.0;
    state_high.Re_D = 1000000.0;

    double Cd_low = Cd_sharp_thin_plate(geom, state_low);
    double Cd_high = Cd_sharp_thin_plate(geom, state_high);

    // Both should be reasonable
    EXPECT_GT(Cd_low, 0.58);
    EXPECT_GT(Cd_high, 0.58);

    // Difference should be small at high Re
    EXPECT_LT(std::abs(Cd_high - Cd_low), 0.02);
}

// -------------------------------------------------------------
// Idelchik (1966) wall-orifice Cd tests
// -------------------------------------------------------------
//
// Ground truth is Idelchik's own tabulated zeta, read off diagrams 4-17 and
// 4-18 rendered at 400 dpi. Cd = 1/sqrt(zeta) is exact for this geometry, so
// these are equality checks at the knots, not directional ones. The tests
// they replace asserted only that thick > thin and round > sharp, which the
// 11% jump discontinuity at Re = 1e5 passed without complaint.

namespace {

DischargeHoleGeometry wall_hole(double d, double L, double r, double bevel) {
    DischargeHoleGeometry hole;
    hole.d = d;
    hole.L = L;
    hole.r = r;
    hole.bevel = bevel;
    return hole;
}

DischargeHoleState wall_flow(double Re) {
    DischargeHoleState flow;
    flow.Re = Re;
    flow.U1_over_Vi = 0.0;
    return flow;
}

}  // namespace

TEST(IdelchikWallOrifice, SharpEdgedMatchesDiagram417AtHighRe) {
    auto c = make_discharge_correlation(DischargeCdCorrelation::Idelchik1966Sharp);
    // zeta reaches 2.85 at the table's LAST knot, Re = 1e6. Item 1's
    // "Re >= 1e5 -> 2.85" is the coarse statement of the same curve; the
    // table itself says 2.60 at 1e5, and treating the two as separate
    // branches put a 4.5% step there.
    const double expected = 1.0 / std::sqrt(orifice::idelchik::zeta_sharp);
    EXPECT_NEAR(c->Cd(wall_hole(1e-3, 0.0, 0.0, 0.0), wall_flow(1.0e6)), expected, 1e-9);
    EXPECT_NEAR(c->Cd(wall_hole(1e-3, 0.0, 0.0, 0.0), wall_flow(1.0e7)), expected, 1e-9);
}

TEST(IdelchikWallOrifice, LowReBranchReproducesEveryTabulatedPoint) {
    auto c = make_discharge_correlation(DischargeCdCorrelation::Idelchik1966Sharp);
    // Diagram 4-17 item 2: zeta = zeta_phi0(Re) + eps_re(Re), at each of the
    // table's own Re. Interpolation must be exact at its knots.
    for (int i = 0; i < orifice::idelchik::re_n; ++i) {
        const double Re = orifice::idelchik::re_points[i];
        const double zeta = orifice::idelchik::zeta_phi0[i] + orifice::idelchik::eps_re[i];
        const double expected = 1.0 / std::sqrt(zeta);
        EXPECT_NEAR(c->Cd(wall_hole(1e-3, 0.0, 0.0, 0.0), wall_flow(Re)), expected, 1e-10)
            << "Re = " << Re;
    }
}

TEST(IdelchikWallOrifice, RoundedAndBeveledReproduceTheirTables) {
    auto cr = make_discharge_correlation(DischargeCdCorrelation::Idelchik1966Rounded);
    for (int i = 0; i < orifice::idelchik::rounded_n; ++i) {
        const double rd = orifice::idelchik::rounded_r_over_d[i];
        const double expected = 1.0 / std::sqrt(orifice::idelchik::rounded_zeta[i]);
        EXPECT_NEAR(cr->Cd(wall_hole(1e-3, 0.0, rd * 1e-3, 0.0), wall_flow(1.0e6)),
                    expected, 1e-10) << "r/d = " << rd;
    }

    auto cb = make_discharge_correlation(DischargeCdCorrelation::Idelchik1966Beveled);
    for (int i = 0; i < orifice::idelchik::beveled_n; ++i) {
        const double ld = orifice::idelchik::beveled_l_over_d[i];
        const double expected = 1.0 / std::sqrt(orifice::idelchik::beveled_zeta[i]);
        EXPECT_NEAR(cb->Cd(wall_hole(1e-3, 0.0, 0.0, ld * 1e-3), wall_flow(1.0e6)),
                    expected, 1e-10) << "l/d = " << ld;
    }
}

TEST(IdelchikWallOrifice, EveryEdgeTypeDegradesToSharpAtZero) {
    // All three edge tables are anchored on zeta = 2.85, so a hole with no
    // radius, no bevel and no depth must give the sharp-edged answer whichever
    // correlation is asked. This is what makes the selector safe to sweep.
    const double sharp = 1.0 / std::sqrt(orifice::idelchik::zeta_sharp);
    const auto hole = wall_hole(1e-3, 0.0, 0.0, 0.0);
    for (auto id : {DischargeCdCorrelation::Idelchik1966Sharp,
                    DischargeCdCorrelation::Idelchik1966Thick,
                    DischargeCdCorrelation::Idelchik1966Beveled,
                    DischargeCdCorrelation::Idelchik1966Rounded}) {
        auto c = make_discharge_correlation(id);
        EXPECT_NEAR(c->Cd(hole, wall_flow(1.0e6)), sharp, 1e-9) << c->name();
    }
}

TEST(IdelchikWallOrifice, CdIsContinuousAcrossTheWholeReynoldsRange) {
    // The correlation this replaces had an 11.1% JUMP at Re = 1e5 (0.881391
    // -> 0.979326), a C0 break in a solver input. Walk the range finely and
    // assert no step larger than what a smooth curve can produce.
    auto c = make_discharge_correlation(DischargeCdCorrelation::Idelchik1966Sharp);
    const auto hole = wall_hole(1e-3, 0.0, 0.0, 0.0);
    double prev = c->Cd(hole, wall_flow(25.0));
    for (int i = 1; i <= 4000; ++i) {
        const double Re = 25.0 * std::pow(1.0e6 / 25.0, i / 4000.0);
        const double cd = c->Cd(hole, wall_flow(Re));
        EXPECT_LT(std::abs(cd - prev), 2.0e-3) << "step at Re = " << Re;
        prev = cd;
    }
}

TEST(IdelchikWallOrifice, AnalyticDerivativeAgreesWithFiniteDifferences) {
    // The (f, J) rule: dCd/dRe is differentiated out of the Hermite form, so
    // it must match a central difference everywhere, including across table
    // knots where a linear interpolant's derivative would be undefined.
    for (auto id : {DischargeCdCorrelation::Idelchik1966Sharp,
                    DischargeCdCorrelation::Idelchik1966Thick}) {
        auto c = make_discharge_correlation(id);
        const auto hole = wall_hole(1e-3, 2e-3, 0.0, 0.0);
        for (double Re : {5.0e1, 3.0e2, 1.5e3, 7.0e3, 3.0e4, 8.0e4, 5.0e5}) {
            const double h = Re * 1e-6;
            const auto fj = c->Cd_and_derivatives(hole, wall_flow(Re));
            const double fd = (c->Cd(hole, wall_flow(Re + h))
                             - c->Cd(hole, wall_flow(Re - h))) / (2.0 * h);
            EXPECT_NEAR(std::get<1>(fj), fd, std::max(1e-8, std::abs(fd) * 1e-5))
                << c->name() << " at Re = " << Re;
        }
    }
}

TEST(IdelchikWallOrifice, BeveledAndRoundedHaveNoReynoldsSlope) {
    // Diagrams 4-18b and 4-18c tabulate zeta against the edge geometry only,
    // so Cd is constant in Re and dCd/dRe is 0. Their interpolant's slope is
    // d/d(edge/Dh); it was returned as the Re slope (0.83 and 1.53 here),
    // which made the network row's flow columns vanish (#964). Mid-table
    // geometry, where that geometry slope is not zero.
    const struct {
        DischargeCdCorrelation id;
        DischargeHoleGeometry hole;
    } cases[] = {
        {DischargeCdCorrelation::Idelchik1966Beveled, wall_hole(1e-2, 0.0, 0.0, 1e-3)},
        {DischargeCdCorrelation::Idelchik1966Rounded, wall_hole(1e-2, 0.0, 1e-3, 0.0)},
    };
    for (const auto& k : cases) {
        auto c = make_discharge_correlation(k.id);
        for (double Re : {1.0e3, 3.0e4, 2.7e5}) {
            const double h = Re * 1e-6;
            const double fd = (c->Cd(k.hole, wall_flow(Re + h))
                             - c->Cd(k.hole, wall_flow(Re - h))) / (2.0 * h);
            EXPECT_EQ(fd, 0.0) << c->name();
            EXPECT_EQ(std::get<1>(c->Cd_and_derivatives(k.hole, wall_flow(Re))), 0.0)
                << c->name() << " at Re = " << Re;
        }
    }
}

TEST(IdelchikWallOrifice, CrossflowDerivativeIsExactlyZeroNotApproximately) {
    // Idelchik's geometry is plenum to plenum: there is no approach velocity,
    // so the absence of a crossflow term is a property of the source, not a
    // truncation. Exact zero, and a test that says why.
    for (auto id : {DischargeCdCorrelation::Idelchik1966Sharp,
                    DischargeCdCorrelation::Idelchik1966Rounded}) {
        auto c = make_discharge_correlation(id);
        const auto fj = c->Cd_and_derivatives(wall_hole(1e-3, 0.0, 1e-4, 0.0),
                                              wall_flow(5.0e4));
        EXPECT_EQ(std::get<2>(fj), 0.0);
    }
}

TEST(IdelchikWallOrifice, MonotoneInterpolantNeverLeavesTheTabulatedEnvelope) {
    // Fritsch-Carlson cannot overshoot; a natural cubic would, on flat tails
    // such as zeta_thick's 1.58, 1.55, 1.55, and would invent a Cd above the
    // source's. Checked on the beveled and rounded correlations, where zeta
    // IS the interpolant -- no friction term and no Re factor -- so the
    // envelope bound is exact rather than approximate.
    struct Case {
        DischargeCdCorrelation id;
        const double* xs;
        const double* ys;
        int n;
        bool by_radius;
    };
    const Case cases[] = {
        {DischargeCdCorrelation::Idelchik1966Beveled,
         orifice::idelchik::beveled_l_over_d, orifice::idelchik::beveled_zeta,
         orifice::idelchik::beveled_n, false},
        {DischargeCdCorrelation::Idelchik1966Rounded,
         orifice::idelchik::rounded_r_over_d, orifice::idelchik::rounded_zeta,
         orifice::idelchik::rounded_n, true},
    };

    for (const auto& tc : cases) {
        auto c = make_discharge_correlation(tc.id);
        for (int i = 0; i < tc.n - 1; ++i) {
            const double z_min = std::min(tc.ys[i], tc.ys[i + 1]);
            const double z_max = std::max(tc.ys[i], tc.ys[i + 1]);
            for (int k = 0; k <= 20; ++k) {
                const double x = tc.xs[i] + (tc.xs[i + 1] - tc.xs[i]) * k / 20.0;
                const auto hole = tc.by_radius ? wall_hole(1e-3, 0.0, x * 1e-3, 0.0)
                                               : wall_hole(1e-3, 0.0, 0.0, x * 1e-3);
                const double cd = c->Cd(hole, wall_flow(1.0e6));
                const double zeta = 1.0 / (cd * cd);
                EXPECT_GE(zeta, z_min - 1e-9) << c->name() << " at x = " << x;
                EXPECT_LE(zeta, z_max + 1e-9) << c->name() << " at x = " << x;
            }
        }
    }
}

TEST(IdelchikWallOrifice, ThickHoleCdRisesWithDepthThenFrictionTakesOver) {
    // zeta' falls from 2.85 to 1.55 as the hole deepens (vena-contracta
    // reattachment) and then goes flat, while the lam*l/Dh friction term
    // keeps growing. So Cd must rise, peak, and fall -- and the peak is the
    // source's, not a fitted one. Idelchik's table goes flat at l/Dh = 2.0.
    auto c = make_discharge_correlation(DischargeCdCorrelation::Idelchik1966Thick);
    double best_cd = 0.0;
    double best_ld = 0.0;
    for (int i = 0; i <= 200; ++i) {
        const double ld = 4.0 * i / 200.0;
        const double cd = c->Cd(wall_hole(1e-3, ld * 1e-3, 0.0, 0.0), wall_flow(1.0e5));
        if (cd > best_cd) {
            best_cd = cd;
            best_ld = ld;
        }
    }
    EXPECT_GT(best_ld, 1.5);
    EXPECT_LT(best_ld, 3.0);
    // And the far end is below the peak: friction has taken over.
    EXPECT_LT(c->Cd(wall_hole(1e-3, 4e-3, 0.0, 0.0), wall_flow(1.0e5)), best_cd);
}

TEST(IdelchikWallOrifice, LaminarBoreFrictionDominatesADeepHoleAtLowRe) {
    // Falsification found this gap: swapping the laminar branch (64/Re) for
    // the turbulent one changed no test result, so nothing pinned it.
    //
    // It matters. At Re = 25 with l/d = 2, lam = 64/25 = 2.56 and the
    // friction term lam*l/Dh = 5.12 is LARGER than zeta' itself (1.55), so
    // omitting the laminar branch would overstate Cd by more than a factor
    // of two. A deep hole at creeping flow is a Poiseuille pipe, not an
    // orifice.
    auto c = make_discharge_correlation(DischargeCdCorrelation::Idelchik1966Thick);
    const auto deep = wall_hole(1e-3, 2e-3, 0.0, 0.0);   // l/d = 2

    const double cd_creep = c->Cd(deep, wall_flow(25.0));
    // zeta = zeta_phi0(25) + k*eps_re(25)*zeta'(2) + 64/25 * 2
    //      = 1.94 + (1/2.85)*1.00*1.55 + 5.12
    const double zeta_expected = 1.94 + (1.0 / 2.85) * 1.00 * 1.55 + 64.0 / 25.0 * 2.0;
    EXPECT_NEAR(cd_creep, 1.0 / std::sqrt(zeta_expected), 1e-9);

    // Without the laminar branch Haaland would give lam ~ 0.08 here, a
    // friction term of 0.16 instead of 5.12 -- so the Cd must be far below
    // what the turbulent branch alone would predict.
    EXPECT_LT(cd_creep, 0.5 * c->Cd(deep, wall_flow(1.0e5)));

    // A hole with no depth has no bore to rub against, laminar or not.
    const auto flat = wall_hole(1e-3, 0.0, 0.0, 0.0);
    EXPECT_NEAR(c->Cd(flat, wall_flow(25.0)),
                c->Cd(flat, wall_flow(25.0)), 0.0);
    const double zeta_flat = 1.94 + (1.0 / 2.85) * 1.00 * 2.85;
    EXPECT_NEAR(c->Cd(flat, wall_flow(25.0)), 1.0 / std::sqrt(zeta_flat), 1e-9);
}

TEST(IdelchikWallOrifice, FrictionBlendIsContinuousAcrossTheTransitionWindow) {
    // The laminar/turbulent blend over 2300 < Re < 4000 is a numerical
    // device, so it must at least not put a step where it removes one.
    auto c = make_discharge_correlation(DischargeCdCorrelation::Idelchik1966Thick);
    const auto deep = wall_hole(1e-3, 2e-3, 0.0, 0.0);
    double prev = c->Cd(deep, wall_flow(2000.0));
    for (int i = 1; i <= 2000; ++i) {
        const double Re = 2000.0 * std::pow(5000.0 / 2000.0, i / 2000.0);
        const double cd = c->Cd(deep, wall_flow(Re));
        EXPECT_LT(std::abs(cd - prev), 1.0e-3) << "step at Re = " << Re;
        prev = cd;
    }
}

TEST(IdelchikWallOrifice, AgreesWithMcGreehanSchotschAcrossTheSharedRange) {
    // CROSS-SOURCE accuracy, labelled as such per the validation policy:
    // Idelchik (1966, Russian handbook) against McGreehan & Schotsch (1988,
    // gas-turbine cooling paper). These are independent measurements of the
    // same physics and they agree to within 5% over the whole L/d range, and
    // to 0.05% at the sharp-edged baseline. A regression on either side that
    // broke that agreement would be worth knowing about.
    auto c = make_discharge_correlation(DischargeCdCorrelation::Idelchik1966Thick);
    const double Re = 1.0e5;
    for (int i = 0; i < orifice::idelchik::thick_n; ++i) {
        const double ld = orifice::idelchik::thick_l_over_d[i];
        const auto hole = wall_hole(1e-3, ld * 1e-3, 0.0, 0.0);
        const double cd_i = c->Cd(hole, wall_flow(Re));
        const double cd_m = orifice::Cd_McGreehanSchotsch(Re, 0.0, ld, 0.0);
        EXPECT_NEAR(cd_m, cd_i, 0.06 * cd_i) << "l/d = " << ld;
    }
    // The sharp-edged anchor is the tight one: 0.5923 vs 0.5926.
    EXPECT_NEAR(orifice::Cd_McGreehanSchotsch(Re, 0.0, 0.0, 0.0),
                1.0 / std::sqrt(orifice::idelchik::zeta_sharp), 5e-4);
}

TEST(Lichtarowicz1965, ReproducesTheSourcesOwnFigure10) {
    // Fig. 10 plots Eq. (12) at l/d = 2 from Re ~ 1 to 1e6: it climbs from
    // near zero, crosses 0.5 around Re ~ 160, and plateaus on Eq. (7).
    auto c = make_discharge_correlation(DischargeCdCorrelation::Lichtarowicz1965);
    const auto hole = wall_hole(1e-3, 2e-3, 0.0, 0.0);   // l/d = 2

    EXPECT_NEAR(c->Cd(hole, wall_flow(10.0)), 0.0817, 5e-4);
    EXPECT_NEAR(c->Cd(hole, wall_flow(100.0)), 0.4284, 5e-4);
    EXPECT_NEAR(c->Cd(hole, wall_flow(1000.0)), 0.7446, 5e-4);

    // The curve crosses 0.5 between Re = 100 and 300, as Fig. 10 shows.
    EXPECT_LT(c->Cd(hole, wall_flow(100.0)), 0.5);
    EXPECT_GT(c->Cd(hole, wall_flow(300.0)), 0.5);
}

TEST(Lichtarowicz1965, ApproachesEquationSevenAtHighReynolds) {
    // Eq. (12) must reduce to Eq. (7) as the viscous and transition terms
    // die away. Checked across the stated l/d range at the stated upper Re.
    auto c = make_discharge_correlation(DischargeCdCorrelation::Lichtarowicz1965);
    for (double ld : {2.0, 4.0, 6.0, 8.0, 10.0}) {
        const double cdu = orifice::lichtarowicz::cdu_c0
                         - orifice::lichtarowicz::cdu_c1 * ld;
        const auto hole = wall_hole(1e-3, ld * 1e-3, 0.0, 0.0);
        EXPECT_NEAR(c->Cd(hole, wall_flow(1.0e6)), cdu, 2e-3) << "l/d = " << ld;
        // At the validated ceiling it is already within the source's own
        // +/-0.02 of the asymptote.
        EXPECT_NEAR(c->Cd(hole, wall_flow(2.0e4)), cdu, 0.02) << "l/d = " << ld;
    }
    // And the flat branch the source gives for 1.5 <= l/d < 2.
    const auto shortish = wall_hole(1e-3, 1.7e-3, 0.0, 0.0);
    EXPECT_NEAR(c->Cd(shortish, wall_flow(1.0e6)),
                orifice::lichtarowicz::cdu_short, 2e-3);
}

TEST(Lichtarowicz1965, RefusesBelowTheSourcesOwnValidityFloor) {
    // Design recommendation (1): avoid l/d < 1.5, "the discharge coefficient
    // varies rapidly with l/d below this value, and there is the possibility
    // of hysteresis". A single-valued correlation cannot represent
    // hysteresis, so refusing beats returning a number.
    auto c = make_discharge_correlation(DischargeCdCorrelation::Lichtarowicz1965);
    EXPECT_THROW(c->Cd(wall_hole(1e-3, 1.0e-3, 0.0, 0.0), wall_flow(1e4)),
                 std::invalid_argument);
    EXPECT_THROW(c->Cd(wall_hole(1e-3, 1.4e-3, 0.0, 0.0), wall_flow(1e4)),
                 std::invalid_argument);
    // At and above the floor it answers.
    EXPECT_GT(c->Cd(wall_hole(1e-3, 1.5e-3, 0.0, 0.0), wall_flow(1e4)), 0.0);
}

TEST(Lichtarowicz1965, HoldsGeometryAboveTheTabulatedRangeRatherThanExtrapolating) {
    // Eq. (7) is LINEAR in l/d, so extrapolating past 10 walks Cd down
    // without bound and reaches zero near l/d = 97. Held instead.
    auto c = make_discharge_correlation(DischargeCdCorrelation::Lichtarowicz1965);
    const double at10 = c->Cd(wall_hole(1e-3, 10e-3, 0.0, 0.0), wall_flow(1e4));
    for (double ld : {12.0, 20.0, 100.0}) {
        EXPECT_NEAR(c->Cd(wall_hole(1e-3, ld * 1e-3, 0.0, 0.0), wall_flow(1e4)),
                    at10, 1e-12) << "l/d = " << ld;
    }
    EXPECT_GT(at10, 0.7);
}

TEST(Lichtarowicz1965, AnalyticDerivativeAgreesWithFiniteDifferences) {
    auto c = make_discharge_correlation(DischargeCdCorrelation::Lichtarowicz1965);
    for (double ld : {2.0, 5.0, 9.0}) {
        const auto hole = wall_hole(1e-3, ld * 1e-3, 0.0, 0.0);
        for (double Re : {2.0e1, 5.0e1, 5.0e2, 5.0e3, 1.5e4}) {
            const double h = Re * 1e-6;
            const auto fj = c->Cd_and_derivatives(hole, wall_flow(Re));
            const double fd = (c->Cd(hole, wall_flow(Re + h))
                             - c->Cd(hole, wall_flow(Re - h))) / (2.0 * h);
            EXPECT_NEAR(std::get<1>(fj), fd, std::max(1e-12, std::abs(fd) * 1e-5))
                << "l/d = " << ld << " Re = " << Re;
        }
    }
}

TEST(Lichtarowicz1965, HasNoCrossflowTermAndSaysSoExactly) {
    // The experiments are plenum-fed; there is no approach velocity. Exact
    // zero, like Idelchik, not a small number.
    auto c = make_discharge_correlation(DischargeCdCorrelation::Lichtarowicz1965);
    const auto fj = c->Cd_and_derivatives(wall_hole(1e-3, 3e-3, 0.0, 0.0),
                                          wall_flow(5.0e3));
    EXPECT_EQ(std::get<2>(fj), 0.0);
}

TEST(Lichtarowicz1965, IsMonotoneAcrossTheValidatedRange) {
    // Cd must rise with Re throughout the range the source validates,
    // Re = 10 to 2e4: more inertia, less viscous loss.
    auto c = make_discharge_correlation(DischargeCdCorrelation::Lichtarowicz1965);
    const auto hole = wall_hole(1e-3, 4e-3, 0.0, 0.0);
    const double lo = orifice::lichtarowicz::re_validated_min;
    const double hi = orifice::lichtarowicz::re_validated_max;

    double prev = c->Cd(hole, wall_flow(lo));
    for (int i = 1; i <= 2000; ++i) {
        const double Re = lo * std::pow(hi / lo, i / 2000.0);
        const double cd = c->Cd(hole, wall_flow(Re));
        EXPECT_TRUE(std::isfinite(cd));
        EXPECT_GT(cd, prev) << "not monotone at Re = " << Re;
        prev = cd;
    }
    EXPECT_LT(prev, 1.0);
}

TEST(Lichtarowicz1965, StaysFiniteAndBoundedOutsideTheValidatedRange) {
    // BELOW: the 20/Re term would diverge and drive Cd to zero, so Re is
    // held at a floor. The value is held and the derivative reported as
    // zero, because the value genuinely stops changing.
    auto c = make_discharge_correlation(DischargeCdCorrelation::Lichtarowicz1965);
    const auto hole = wall_hole(1e-3, 4e-3, 0.0, 0.0);
    const double at_floor = c->Cd(hole, wall_flow(orifice::lichtarowicz::re_floor));
    EXPECT_TRUE(std::isfinite(at_floor));
    EXPECT_GT(at_floor, 0.0);
    EXPECT_NEAR(c->Cd(hole, wall_flow(1e-9)), at_floor, 1e-12);
    EXPECT_EQ(std::get<1>(c->Cd_and_derivatives(hole, wall_flow(1e-9))), 0.0);

    // ABOVE: Eq. (12) is NOT monotone forever. Its two Re terms pull
    // opposite ways -- 20(1+2.25 l/d)/Re falls without limit, while the
    // log-squared term peaks at Re = 1/0.00015 = 6667 and decays either
    // side. Past Re ~ 8.7e5 the viscous term is spent while the log term is
    // still decaying, so Cd overshoots Eq. (7) and settles back.
    //
    // Recorded rather than smoothed away: the overshoot is 2.2e-4, which is
    // 43x beyond the source's validated ceiling of 2e4 and roughly 100x
    // SMALLER than its own stated accuracy of +/-0.02. It is an artifact of
    // extrapolating the fit, not a defect to correct.
    const double cdu = orifice::lichtarowicz::cdu_c0
                     - orifice::lichtarowicz::cdu_c1 * 4.0;
    double worst = 0.0;
    for (double Re : {1.0e5, 3.0e5, 8.7e5, 1.0e6, 1.0e7, 1.0e9}) {
        const double cd = c->Cd(hole, wall_flow(Re));
        EXPECT_TRUE(std::isfinite(cd));
        worst = std::max(worst, std::abs(cd - cdu));
    }
    EXPECT_LT(worst, 1.0e-3);
    EXPECT_LT(worst, 0.02);   // far inside the source's own accuracy
}

TEST(Lichtarowicz1965, DisagreesWithMcGreehanWhereMcGreehanIsFloored) {
    // CROSS-SOURCE, and the reason this correlation was added. McGreehan's
    // chain floors Re at 1e4, so below that its Cd stops moving; Lichtarowicz
    // is validated from Re = 10 and keeps falling. At Re = 432 -- inside
    // Andrews' own effusion data -- they differ by more than 20%.
    auto cl = make_discharge_correlation(DischargeCdCorrelation::Lichtarowicz1965);
    auto cm = make_discharge_correlation(DischargeCdCorrelation::McGreehanSchotsch1988);
    const auto hole = wall_hole(3.27e-3, 6.3e-3, 0.0, 0.0);   // Andrews plate C

    const double lo_l = cl->Cd(hole, wall_flow(432.0));
    const double lo_m = cm->Cd(hole, wall_flow(432.0));
    EXPECT_GT((lo_m - lo_l) / lo_l, 0.15) << "expected McGreehan high at low Re";

    // And they converge where both are valid.
    const double hi_l = cl->Cd(hole, wall_flow(2.0e4));
    const double hi_m = cm->Cd(hole, wall_flow(2.0e4));
    EXPECT_NEAR(hi_m, hi_l, 0.03 * hi_l);
}


// -------------------------------------------------------------
// Flow calculation tests
// -------------------------------------------------------------

TEST_F(OrificeTest, OrificeMdot) {
    double Cd_val = 0.61;
    double mdot = orifice_mdot(geom, Cd_val, state.dP, state.rho);

    // Expected: Cd * E * A * sqrt(2 * rho * dP)
    const double beta = geom.beta();
    const double E = 1.0 / std::sqrt(1.0 - std::pow(beta, 4.0));
    double expected = Cd_val * E * geom.area() * std::sqrt(2.0 * state.rho * state.dP);
    EXPECT_NEAR(mdot, expected, 1e-10);
}

TEST_F(OrificeTest, OrificeDpRoundTrip) {
    double Cd_val = 0.61;
    double mdot = orifice_mdot(geom, Cd_val, state.dP, state.rho);
    double dP_calc = orifice_dP(geom, Cd_val, mdot, state.rho);
    EXPECT_NEAR(dP_calc, state.dP, 1e-6);
}

TEST_F(OrificeTest, OrificeCdFromMeasurement) {
    double Cd_val = 0.61;
    double mdot = orifice_mdot(geom, Cd_val, state.dP, state.rho);
    double Cd_calc = orifice_Cd_from_measurement(geom, mdot, state.dP, state.rho);
    EXPECT_NEAR(Cd_calc, Cd_val, 1e-10);
}

// -------------------------------------------------------------
// Correlation class tests
// -------------------------------------------------------------

TEST_F(OrificeTest, CorrelationFactory) {
    auto corr = make_correlation(MeteringCdCorrelation::ReaderHarrisGallagher);
    ASSERT_NE(corr, nullptr);
    EXPECT_FALSE(corr->name().empty());

    double Cd_class = corr->Cd(geom, state);
    double Cd_func = Cd_sharp_thin_plate(geom, state);
    EXPECT_DOUBLE_EQ(Cd_class, Cd_func);
}

TEST_F(OrificeTest, UserCorrelation) {
    auto corr = make_user_correlation(
        [](const OrificeGeometry&, const OrificeState&) { return 0.65; },
        "TestConstant");

    ASSERT_NE(corr, nullptr);
    EXPECT_EQ(corr->name(), "TestConstant");
    EXPECT_DOUBLE_EQ(corr->Cd(geom, state), 0.65);
}

// -------------------------------------------------------------
// Namespace function tests
// -------------------------------------------------------------

TEST_F(OrificeTest, KFromCd) {
    double Cd_val = 0.61;
    double beta = geom.beta();
    double K = orifice::K_from_Cd(Cd_val, beta);

    // K should be positive for Cd < 1
    EXPECT_GT(K, 0.0);

    // Round-trip
    double Cd_back = orifice::Cd_from_K(K, beta);
    EXPECT_NEAR(Cd_back, Cd_val, 1e-10);
}

// -------------------------------------------------------------
// Hardening and Numerical Stability tests
// -------------------------------------------------------------

TEST_F(OrificeTest, HardeningHighBeta) {
    // Test that Cd remains finite for beta -> 1.0 (approaching pipe diameter)
    OrificeGeometry geom_extreme;
    geom_extreme.D = 0.100;
    geom_extreme.d = 0.0999999; // beta = 0.9999999

    double Cd_val = Cd_sharp_thin_plate(geom_extreme, state);
    // Should be capped at 1.5, not 500,000
    EXPECT_LE(Cd_val, 1.5);
    EXPECT_GT(Cd_val, 0.5);
}

TEST_F(OrificeTest, KFromCdHighBeta) {
    // Test that K_from_Cd and Cd_from_K are stable and self-consistent near beta = 1.0
    double beta_extreme = 0.99999;
    double Cd_val = 0.95;

    double K = orifice::K_from_Cd(Cd_val, beta_extreme);
    EXPECT_TRUE(std::isfinite(K));
    EXPECT_GE(K, 0.0);

    // Round-trip should be exact even with clamping because both use the same internal clamp
    double Cd_back = orifice::Cd_from_K(K, beta_extreme);
    EXPECT_NEAR(Cd_back, Cd_val, 1e-10);
}

// -------------------------------------------------------------
// Metering Cd held below its lower Re_D bound (#873)
// -------------------------------------------------------------
//
// RHG was evaluated far below its Re_D range, where the Re terms grow
// without bound: a 4 mm, beta 0.5 plate hit the 1.5 cap below Re_D ~ 40. It
// is now held below Re_b = metering_cd_re_min, entered C1 in ln Re_D over
// orifice::metering_hold::blend_ln_re.

namespace {

OrificeGeometry metering_plate(double beta, double D, OrificeTap tap) {
    OrificeGeometry g;
    g.D = D;
    g.d = beta * D;
    g.tap = tap;
    return g;
}

struct MeteringCase {
    MeteringCdCorrelation id;
    OrificeGeometry geom;
};

std::vector<MeteringCase> metering_cases() {
    std::vector<MeteringCase> out;
    for (OrificeTap tap : {OrificeTap::Corner, OrificeTap::D_D2, OrificeTap::Flange}) {
        for (double D : {0.004, 0.1, 1.0}) {
            for (double beta : {0.2, 0.5, 0.6, 0.75}) {
                out.push_back(
                    {MeteringCdCorrelation::ReaderHarrisGallagher, metering_plate(beta, D, tap)});
            }
        }
    }
    return out;
}

}  // namespace

TEST(MeteringCdHold, IsConstantWithZeroSlopeBelowTheBlend) {
    const double w = orifice::metering_hold::blend_ln_re;
    for (const auto& c : metering_cases()) {
        const double Rb = metering_cd_re_min(c.id, c.geom);
        const auto [C_b, dC_b] = metering_cd_and_dRe(c.id, c.geom, Rb);
        const double held = C_b - 0.5 * dC_b * Rb * w;
        EXPECT_GE(held, C_b);  // monotone: the hold never dips below the bound
        for (double Re : {Rb * std::exp(-w) * 0.999, 0.5 * Rb, 1e-3 * Rb, 1.0, 0.0, -3.0}) {
            const auto [C, dC] = metering_cd_and_dRe(c.id, c.geom, Re);
            EXPECT_NEAR(C, held, 1e-13) << "Re_D " << Re;
            EXPECT_EQ(dC, 0.0) << "Re_D " << Re;
        }
    }
}

TEST(MeteringCdHold, IsC1AndItsDerivativeIsExactAcrossTheBlend) {
    const double w = orifice::metering_hold::blend_ln_re;
    for (const auto& c : metering_cases()) {
        const double Rb = metering_cd_re_min(c.id, c.geom);
        const double s = std::abs(std::get<1>(metering_cd_and_dRe(c.id, c.geom, Rb))) * Rb;
        for (double edge : {Rb, Rb * std::exp(-w)}) {
            const auto [Cl, dl] = metering_cd_and_dRe(c.id, c.geom, edge * (1.0 - 1e-11));
            const auto [Ch, dh] = metering_cd_and_dRe(c.id, c.geom, edge * (1.0 + 1e-11));
            EXPECT_NEAR(Cl, Ch, 3e-11 * s + 1e-14);
            EXPECT_NEAR(dl * edge, dh * edge, 1e-8 * s);
        }
        for (double u : {-0.09, -0.05, -0.01, -1e-4, 1e-4, 0.05}) {
            const double Re = Rb * std::exp(u);
            const double h = 1e-6 * Re;
            const double fd = (std::get<0>(metering_cd_and_dRe(c.id, c.geom, Re + h)) -
                               std::get<0>(metering_cd_and_dRe(c.id, c.geom, Re - h))) /
                              (2.0 * h);
            const double an = std::get<1>(metering_cd_and_dRe(c.id, c.geom, Re));
            EXPECT_NEAR(an * Re, fd * Re, 1e-7 * s) << "u " << u;
            if (u < 0.0) {
                // The log-slope ramps linearly from its bound value to 0.
                const double s_b = std::get<1>(metering_cd_and_dRe(c.id, c.geom, Rb)) * Rb;
                EXPECT_NEAR(an * Re, s_b * (1.0 + u / w), 1e-9 * s) << "u " << u;
            }
        }
    }
}

TEST(MeteringCdHold, TheIssueExampleNoLongerClimbsToTheCap) {
    // D 4 mm, beta 0.5, D-D/2: 1.5 at Re_D 1-40, 1.162 at 100, 0.687 at 1000
    // before; 0.6345 everywhere below 0.905 Re_b now.
    const auto g = metering_plate(0.5, 0.004, OrificeTap::D_D2);
    for (double Re : {1.0, 40.0, 100.0, 1000.0}) {
        EXPECT_NEAR(std::get<0>(metering_cd_and_dRe(MeteringCdCorrelation::ReaderHarrisGallagher,
                                                    g, Re)),
                    0.6345, 5e-5);
    }
}

TEST(MeteringCdHold, BoundsAreTheStandardsAndTheFlagSaysBelowThem) {
    using M = MeteringCdCorrelation;
    EXPECT_EQ(metering_cd_re_min(M::ReaderHarrisGallagher, metering_plate(0.5, 0.1, OrificeTap::Corner)), 5000.0);
    EXPECT_NEAR(metering_cd_re_min(M::ReaderHarrisGallagher, metering_plate(0.7, 0.1, OrificeTap::D_D2)), 7840.0, 1e-9);
    EXPECT_NEAR(metering_cd_re_min(M::ReaderHarrisGallagher, metering_plate(0.7, 0.5, OrificeTap::Flange)), 41650.0, 1e-9);
    EXPECT_THROW(metering_cd_re_min(M::Constant, metering_plate(0.5, 0.1, OrificeTap::Corner)),
                 std::invalid_argument);
    const auto g = metering_plate(0.5, 0.1, OrificeTap::D_D2);
    const double Rb = metering_cd_re_min(M::ReaderHarrisGallagher, g);
    EXPECT_TRUE(metering_cd_in_range(M::ReaderHarrisGallagher, g, Rb));
    EXPECT_FALSE(metering_cd_in_range(M::ReaderHarrisGallagher, g, 0.999 * Rb));
    EXPECT_FALSE(metering_cd_in_range(M::ReaderHarrisGallagher, g, 0.0));
    EXPECT_TRUE(metering_cd_in_range(M::Constant, g, 0.0));  // the caller's Cd
}

TEST(IdelchikWallOrifice, DegenerateGeometryStaysFinite) {
    // The beta -> 1 singularity this replaces cannot arise here: a hole in a
    // wall has no pipe to form beta with. What can still be fed in is a
    // radius or depth far outside the table, which must clamp rather than
    // extrapolate into a negative zeta.
    DischargeHoleGeometry hole;
    hole.d = 1e-3;
    hole.r = 1.0;      // r/d = 1000, far past the table's 0.20
    hole.L = 1.0;      // l/d = 1000
    DischargeHoleState flow;
    flow.Re = 1.0e5;

    for (auto id : {DischargeCdCorrelation::Idelchik1966Rounded,
                    DischargeCdCorrelation::Idelchik1966Thick}) {
        auto c = make_discharge_correlation(id);
        const double cd = c->Cd(hole, flow);
        EXPECT_TRUE(std::isfinite(cd)) << c->name();
        EXPECT_GT(cd, 0.0) << c->name();
    }

    // And below the table's Re floor of 25 the value holds rather than
    // running off: Idelchik stops there and so do we.
    auto c = make_discharge_correlation(DischargeCdCorrelation::Idelchik1966Sharp);
    DischargeHoleGeometry sharp;
    sharp.d = 1e-3;
    DischargeHoleState creep;
    creep.Re = 1.0e-3;
    const double cd_floor = c->Cd(sharp, creep);
    creep.Re = orifice::idelchik::re_points[0];
    EXPECT_NEAR(cd_floor, c->Cd(sharp, creep), 1e-12);
}

// -------------------------------------------------------------
// McGreehan and Schotsch (1988)
//
// Every expectation below is anchored OUTSIDE the implementation: a value
// printed elsewhere in the paper, a different author's correlation, a curve
// drawn on the page, a physical identity, or another paper entirely. None is
// a golden value captured from this code. Provenance and the measured
// agreement for each are in
// validation/cooling/extractions/orifice_discharge_coefficient.md.
// -------------------------------------------------------------

namespace ms = orifice::mcgreehan_schotsch;

// The paper states a baseline sharp-edged Cd of 0.60 at Re = 3.2e4 in running
// text on p.213, independently of Eq. (8) on p.214. Eq. (8) must reproduce it.
TEST(McGreehanSchotsch, ReynoldsBaselineReproducesStatedReferencePoint) {
    EXPECT_NEAR(ms::reynolds_baseline(3.2e4), ms::cd_reference, 5e-4);
}

// Eq. (12) against Schoder and Dawson's Eq. (10), dCd/Cd = 3.10 (r/d), which
// the paper cites as "a good approximation for 0 < r/d < 0.1" but does not use.
// A different author's correlation, so this could genuinely have disagreed.
TEST(McGreehanSchotsch, CornerFactorAgreesWithSchoderDawson) {
    for (const double rd : {0.02, 0.05, 0.10}) {
        const double schoder = ms::cd_reference * (1.0 + 3.10 * rd);
        const double got     = 1.0 - ms::corner_factor(rd) * (1.0 - ms::cd_reference);
        EXPECT_NEAR(got, schoder, 0.015 * schoder) << "at r/d = " << rd;
    }
}

// p.214: an ASME nozzle Cd is the limiting value, reached at r/d = 0.82.
// Extrapolating the corner-radius fit that far must land on Eq. (9)'s own
// high-Re nozzle asymptote -- a different equation fitted to different data.
TEST(McGreehanSchotsch, CornerRadiusExtrapolatesToNozzleAsymptote) {
    const double fully_rounded =
        1.0 - ms::corner_factor(ms::r_over_d_nozzle_limit) * (1.0 - ms::re_c0);
    EXPECT_NEAR(fully_rounded, ms::nozzle_c0, 0.005 * ms::nozzle_c0);
}

// Identities the fitted constants must satisfy for the chain to be consistent.
// g(0) = 1 is a property of 1.3 and 0.435 conspiring, not something imposed.
TEST(McGreehanSchotsch, FactorsAreIdentitiesAtZero) {
    EXPECT_NEAR(ms::length_factor(0.0), 1.0, 1e-3);
    EXPECT_DOUBLE_EQ(ms::corner_factor(0.0), 1.0);
}

// Eq. (17) must reduce exactly to its input at zero crossflow.
TEST(McGreehanSchotsch, CrossflowIsIdentityAtZero) {
    const double base = ms::cd_with_corner_and_length(1.0e4, 0.0, 1.0);

    // Eq. (17) is an identity at zero crossflow (C1 = 1, C2 = 0), and eps = 0
    // reproduces that exactly. This asserted only the line below until the
    // regularisation landed; it is kept because the paper's own behaviour is
    // what the correlation must still be able to give.
    EXPECT_DOUBLE_EQ(ms::cd_with_crossflow(base, 0.0, 0.0), base);

    // At the default eps the identity holds to the documented cost instead of
    // exactly -- that is the whole point of the regularisation, and the bound
    // is what makes it acceptable rather than the exactness.
    EXPECT_NE(ms::cd(1.0e4, 0.0, 1.0, 0.0), base);
    EXPECT_NEAR(ms::cd(1.0e4, 0.0, 1.0, 0.0), base, 5.0e-4)
        << "the zero-crossflow departure exceeds what rv_smooth_eps documents";
}

// Eq. (13)/(14) against the curve drawn in Fig. 3, read at ~0.01 in Cd.
// Fig. 3 is drawn for a basic Cd of 0.60, which Eq. (8) gives at Re = 3.2e4.
// L/d = 0.5 is excluded: it sits on the steepest part of the curve, where a
// 0.1 error in reading L/d moves Cd by 0.017 (see check E of the extraction).
TEST(McGreehanSchotsch, LengthEffectMatchesFigure3) {
    const struct { double L_over_d; double Cd_drawn; } points[] = {
        {1.0, 0.77}, {2.0, 0.81}, {3.0, 0.81}, {5.0, 0.78}, {8.0, 0.76}, {10.0, 0.74},
    };
    for (const auto& p : points) {
        const double got = ms::cd_with_corner_and_length(3.2e4, 0.0, p.L_over_d);
        EXPECT_NEAR(got, p.Cd_drawn, 0.02 * p.Cd_drawn) << "at L/d = " << p.L_over_d;
    }
}

// Eq. (17) against the Cd:r,L = 0.6 curve drawn in Figs. 4 and 7. A basic Cd
// of 0.60 is reached at Re = 3.2e4 with r/d = L/d = 0.
TEST(McGreehanSchotsch, CrossflowEffectMatchesFigure4) {
    const struct { double U1_over_Vi; double Cd_drawn; } points[] = {
        {1.0, 0.40}, {2.0, 0.24}, {4.0, 0.12},
    };
    for (const auto& p : points) {
        const double got = ms::cd(3.2e4, 0.0, 0.0, p.U1_over_Vi);
        EXPECT_NEAR(got, p.Cd_drawn, 0.05 * p.Cd_drawn) << "at U1/Vi = " << p.U1_over_Vi;
    }
}

// The rise before the fall is drawn on Fig. 4 and supported by Rohde's own
// result 2 (slanting an orifice into the flow increases Cd). It is deliberate,
// so it gets a test: anyone who "fixes" Eq. (17) into a monotone decay, or
// clamps the factor at 1, breaks this.
TEST(McGreehanSchotsch, CrossflowIsNonMonotonicNearTheOrigin) {
    const double base = ms::cd(3.2e4, 0.0, 0.0, 0.0);

    double peak = base;
    double peak_at = 0.0;
    for (int i = 1; i <= 400; ++i) {
        const double u  = 0.005 * i;
        const double cd = ms::cd(3.2e4, 0.0, 0.0, u);
        if (cd > peak) { peak = cd; peak_at = u; }
    }

    // Fig. 4's drawn peak is ~0.635 at U1/Vi ~ 0.1, against a 0.60 baseline.
    EXPECT_GT(peak, base) << "the rise before the fall has been lost";
    EXPECT_NEAR(peak, 0.635, 0.02 * 0.635);
    EXPECT_NEAR(peak_at, 0.1, 0.05);
    // ...and it really does fall away afterwards.
    EXPECT_LT(ms::cd(3.2e4, 0.0, 0.0, 4.0), 0.5 * base);
}

// The #375 anchor. A bare plenum-fed jet plate -- sharp holes, t/d = 1, no
// supply-side crossflow -- against Florschuetz, Truman and Metzger (1981),
// who recommend C_D = 0.79 and measure 0.73-0.85 (their Table 1). A different
// paper, rig and decade, and the value combaero currently hard-codes.
TEST(McGreehanSchotsch, BareJetPlateAgreesWithFlorschuetzDefault) {
    EXPECT_NEAR(ms::cd(1.0e4, 0.0, 1.0, 0.0), 0.79, 0.01 * 0.79);
}

TEST(McGreehanSchotsch, JetPlateStaysInsideFlorschuetzMeasuredBand) {
    for (const double Re : {5.0e3, 1.0e4, 3.0e4, 7.0e4}) {
        for (const double t_over_d : {1.0, 1.5, 2.0, 3.0}) {
            const double got = ms::cd(Re, 0.0, t_over_d, 0.0);
            EXPECT_GE(got, 0.73) << "Re = " << Re << ", t/d = " << t_over_d;
            EXPECT_LE(got, 0.85) << "Re = " << Re << ", t/d = " << t_over_d;
        }
    }
}

// Eqs. (15)/(16) only engage when a corner radius is present, so the chain has
// a branch at r/d = 0. g(0) = 1.0005 rather than exactly 1, so the branch is
// not perfectly continuous. Bound the step, so that if anyone widens it the
// test says so.
TEST(McGreehanSchotsch, CombinedCorrectionBranchIsEffectivelyContinuous) {
    const double at_zero = ms::cd_with_corner_and_length(3.2e4, 0.0, 2.0);
    const double just_above = ms::cd_with_corner_and_length(3.2e4, 1e-9, 2.0);
    EXPECT_NEAR(just_above, at_zero, 5e-4);
}

// Eqs. (13)+(15)+(16) collapse algebraically to a single product form,
//   Cd = 1 - g(L/d - r/d) * g(r/d) * (1 - Cd:r),
// derived by hand from the two sequential steps the implementation actually
// performs. An independent expression of the same thing: it catches a dropped
// Eq. (16) subtraction, a g() evaluated at the wrong argument, or the printed
// "L/D" being taken literally.
TEST(McGreehanSchotsch, CombinedCorrectionMatchesItsClosedProductForm) {
    const struct { double r_over_d; double L_over_d; } cases[] = {
        {0.1, 1.0}, {0.3, 2.0}, {0.5, 3.0}, {0.05, 0.5},
    };
    for (const auto& c : cases) {
        const double cd_r = ms::cd_with_corner(3.2e4, c.r_over_d);
        const double expected =
            1.0 - ms::length_factor(c.L_over_d - c.r_over_d) *
                      ms::length_factor(c.r_over_d) * (1.0 - cd_r);
        EXPECT_NEAR(ms::cd_with_corner_and_length(3.2e4, c.r_over_d, c.L_over_d),
                    expected, 1e-12)
            << "at r/d = " << c.r_over_d << ", L/d = " << c.L_over_d;
    }
}

// Rounding the inlet must raise Cd monotonically, and at the r/d = 0.82 knee
// the paper says an ASME nozzle is reached -- so a long, well-rounded orifice
// must sit close to 1.
TEST(McGreehanSchotsch, CornerRadiusRaisesCdTowardTheNozzleLimit) {
    double previous = ms::cd_with_corner_and_length(3.2e4, 0.0, 2.0);
    for (const double rd : {0.05, 0.1, 0.2, 0.3, 0.5, 0.82}) {
        const double got = ms::cd_with_corner_and_length(3.2e4, rd, 2.0);
        EXPECT_GT(got, previous) << "not monotone at r/d = " << rd;
        previous = got;
    }
    EXPECT_GT(previous, 0.98);
    EXPECT_LE(previous, 1.0);
}

TEST(McGreehanSchotsch, WrapperMatchesNamespacedChain) {
    EXPECT_DOUBLE_EQ(orifice::Cd_McGreehanSchotsch(2.0e4, 0.1, 1.5, 0.3),
                     ms::cd(2.0e4, 0.1, 1.5, 0.3));
}

TEST(McGreehanSchotsch, StaysFiniteOnSolverExcursions) {
    for (const double Re : {1.0e-3, 1.0, 1.0e3, 1.0e9}) {
        for (const double u : {0.0, 1.0e-12, 50.0}) {
            const double got = ms::cd(Re, 0.0, 1.0, u);
            EXPECT_TRUE(std::isfinite(got)) << "Re = " << Re << ", U1/Vi = " << u;
            EXPECT_GT(got, 0.0);
            EXPECT_LE(got, 1.5);
        }
    }
    // Negative inputs are clamped rather than producing NaN from pow().
    EXPECT_TRUE(std::isfinite(ms::cd(1.0e4, -1.0, -1.0, -1.0)));
}

// A constant Cd must be settable, not just selectable: make_correlation() can
// only hand back the default.
TEST(McGreehanSchotsch, ConstantCorrelationTakesAnExplicitValue) {
    OrificeGeometry g;
    g.d = 0.001;
    g.D = 0.010;
    OrificeState s;
    s.Re_D = 1.0e4;

    auto fixed = make_constant_correlation(0.79);
    ASSERT_NE(fixed, nullptr);
    EXPECT_DOUBLE_EQ(fixed->Cd(g, s), 0.79);
}

// Eq. (8) diverges below its stated floor (it reaches Cd = 1.0 at Re = 904),
// so the chain holds Re at re_min rather than extrapolating into nonsense.
// Constant below the floor is a deliberate, documented choice.
TEST(McGreehanSchotsch, ReynoldsBaselineIsHeldAtItsValidityFloor) {
    // Eq. (8) at its own stated floor: 0.5885 + 372/1e4.
    const double at_floor = ms::reynolds_baseline(ms::re_min);
    EXPECT_NEAR(at_floor, 0.6257, 1e-4);

    for (const double Re : {1.0e-3, 1.0, 9.0e2, 5.0e3}) {
        EXPECT_DOUBLE_EQ(ms::reynolds_baseline(Re), at_floor)
            << "not held at the floor for Re = " << Re;
    }
    // Above the floor it must still respond to Re.
    EXPECT_LT(ms::reynolds_baseline(1.0e5), at_floor);
    EXPECT_NEAR(ms::reynolds_baseline(3.2e4), ms::cd_reference, 5e-4);
}

// Eq. (17) applied to a caller-supplied baseline. This is how the source uses
// it in its own validation, so the chain must decompose exactly this way.
TEST(McGreehanSchotsch, ChainDecomposesIntoBaselineAndCrossflow) {
    for (const double u : {0.0, 0.1, 0.5, 1.0, 3.0}) {
        const double base = ms::cd_with_corner_and_length(3.2e4, 0.1, 2.0);
        EXPECT_DOUBLE_EQ(ms::cd(3.2e4, 0.1, 2.0, u),
                         ms::cd_with_crossflow(base, u))
            << "at U1/Vi = " << u;
    }
}

// Eq. (17) is otherwise validated against a DRAWN curve by
// CrossflowEffectMatchesFigure4, and against Rohde's digitised data by
// python/tests/test_orifice_validation.py. No figure-read test is duplicated
// here for the Fig. 6 baseline curves.

// -------------------------------------------------------------
// Expansion factor Y, Eqs. (4)-(7)
// -------------------------------------------------------------

TEST(McGreehanSchotschY, NozzleFormIsTheIsentropicExpansionFactor) {
    // Eq. (5) IS the isentropic nozzle expansion factor. Checked here against
    // the closed-form isentropic mass flux ratio, derived independently:
    //   Y_n = [mdot_isentropic / mdot_incompressible] at the same dP.
    // (The stronger check, against combaero's own nozzle_flow with real gas
    // properties, lives in python/tests/test_orifice_expansion.py.)
    const double g = 1.4;
    for (const double S : {0.99, 0.95, 0.90, 0.80, 0.70, 0.60}) {
        // Isentropic mass flux, normalised the same way Y is defined.
        const double num = std::sqrt(
            (2.0 * g / (g - 1.0)) *
            (std::pow(S, 2.0 / g) - std::pow(S, (g + 1.0) / g)));
        const double den = std::sqrt(2.0 * (1.0 - S));
        EXPECT_NEAR(ms::expansion_nozzle(S, g), num / den, 1e-10)
            << "at S = " << S;
    }
}

TEST(McGreehanSchotschY, BothFormsTendToOneAtZeroPressureDrop) {
    for (const double g : {1.3, 1.4, 1.67}) {
        EXPECT_NEAR(ms::expansion_orifice(1.0, g), 1.0, 1e-12);
        EXPECT_NEAR(ms::expansion_nozzle(1.0, g), 1.0, 1e-9);
        EXPECT_NEAR(ms::expansion_nozzle(1.0 - 1e-10, g), 1.0, 1e-6);
    }
}

TEST(McGreehanSchotschY, OrificeFormMatchesThePrintedCoefficient) {
    // Eq. (4) printed as 1 - 0.41 (P_t1 - P_s2)/(gamma P_t1). Checked against
    // that literal form rather than the (1-S)/gamma rewrite used internally.
    const double g = 1.4, P_t1 = 4.0e5, P_s2 = 3.0e5;
    const double printed = 1.0 - 0.41 * (P_t1 - P_s2) / (g * P_t1);
    EXPECT_NEAR(ms::expansion_orifice(P_s2 / P_t1, g), printed, 1e-12);
}

TEST(McGreehanSchotschY, BlendWeightReducesToThePaperAtZeroSmoothing) {
    // eps = 0 must recover Eq. (7) clamped hard, exactly.
    for (const double cd : {0.70, 0.80, 0.82, 0.86, 0.90, 0.94, 0.99}) {
        const double x = ms::x_blend_slope * (cd - ms::x_blend_cd0);
        const double hard = std::max(0.0, std::min(1.0, x));
        EXPECT_NEAR(ms::expansion_blend_weight(cd, 0.0), hard, 1e-12)
            << "at Cd = " << cd;
    }
    // X reaches 1 at Cd = 0.94, as the paper's constants imply.
    EXPECT_NEAR(ms::expansion_blend_weight(0.94, 0.0), 1.0, 1e-3);
}

// The reviewer's concern, made a test: a hard clamp puts an exactly-zero
// derivative outside [0.82, 0.94] and a discontinuous jump at both knees.
// The smoothed default must have neither.
TEST(McGreehanSchotschY, SmoothedBlendHasNoDeadZoneAndNoJump) {
    const double g = 1.4, S = 0.7;
    const double h = 1e-6;
    auto dY = [&](double cd, double eps) {
        return (ms::expansion_factor(cd + h, S, g, eps)
                - ms::expansion_factor(cd - h, S, g, eps)) / (2.0 * h);
    };

    // Hard clamp: dead on both sides. This documents what is being avoided.
    EXPECT_DOUBLE_EQ(dY(0.78, 0.0), 0.0);
    EXPECT_DOUBLE_EQ(dY(0.99, 0.0), 0.0);

    // Smoothed: not merely non-zero but with ROOM in it. An earlier version
    // of this assertion tested != 0, which passes at eps = 0.005 where the
    // smallest slope is 1.4e-6 -- numerically a floor. 5e-5 is below the
    // 1.36e-4 the default achieves and above the 2.2e-5 that eps = 0.02 gets,
    // so this assertion is itself a lower wall on eps.
    for (const double cd : {0.70, 0.78, 0.818, 0.822, 0.90, 0.936, 0.944,
                            0.97, 0.99}) {
        EXPECT_GT(std::abs(dY(cd, ms::x_smooth_eps)), 5.0e-5)
            << "derivative is effectively floored at Cd = " << cd;
    }

    // ...and CONTINUOUS across both knees. Tested by how the change scales
    // with the sampling interval rather than against a threshold: a genuine
    // discontinuity holds the same jump however closely you sample, while a
    // continuous derivative's change shrinks with the interval. Measured, the
    // hard clamp sits at 0.73353 for every h; the smoothed one halves each
    // time h halves.
    for (const double knee : {0.82, 0.94}) {
        auto jump = [&](double half, double eps) {
            return std::abs(dY(knee + half, eps) - dY(knee - half, eps));
        };

        // Hard clamp: the jump does not shrink. This is what is being avoided.
        EXPECT_NEAR(jump(0.004, 0.0), jump(0.00025, 0.0), 1e-6);
        EXPECT_NEAR(jump(0.004, 0.0), 0.7335, 1e-3);

        // Smoothed: halving the interval halves the change, to within 5%.
        const double wide = jump(0.002, ms::x_smooth_eps);
        const double half = jump(0.001, ms::x_smooth_eps);
        EXPECT_NEAR(half / wide, 0.5, 0.05)
            << "derivative is not continuous at the knee Cd = " << knee;
    }
}

TEST(McGreehanSchotschY, SmoothingCostIsBounded) {
    // What the smoothing buys must be paid for in a quantified amount of Y.
    double worst = 0.0;
    for (int i = 0; i <= 300; ++i) {
        const double cd = 0.70 + 0.30 * i / 300.0;
        for (const double S : {0.6, 0.7, 0.8, 0.9, 0.99}) {
            worst = std::max(worst,
                             std::abs(ms::expansion_factor(cd, S, 1.4, ms::x_smooth_eps)
                                      - ms::expansion_factor(cd, S, 1.4, 0.0)));
        }
    }
    EXPECT_LT(worst, 0.006) << "smoothing deviation grew: " << worst;
}

// Below the critical pressure ratio the isentropic form turns over and
// predicts decreasing flow. Y must saturate instead.
// The evidence for x_smooth_eps, rather than an assurance that it is sensible.
//
// Two independent properties bound it from opposite directions, and the
// default has to sit between them. Measured walls: continuity is satisfied
// for eps >= 0.0304, fidelity for eps <= 0.0976. The default sits 29% across
// that window. If someone moves the constant, this says which wall they hit.
TEST(McGreehanSchotschY, DefaultSmoothingSitsInsideItsAdmissibleWindow) {
    const double g = 1.4, S = 0.7, h = 1e-6;
    auto dY = [&](double cd, double eps) {
        return (ms::expansion_factor(cd + h, S, g, eps)
                - ms::expansion_factor(cd - h, S, g, eps)) / (2.0 * h);
    };
    // Continuity, as the ratio of the derivative change at two sampling
    // intervals. 0.5 means continuous; 1.0 means a kink.
    auto continuity = [&](double eps) {
        auto jump = [&](double half) {
            return std::abs(dY(0.94 + half, eps) - dY(0.94 - half, eps));
        };
        const double wide = jump(0.002);
        return wide > 0.0 ? jump(0.001) / wide : 1.0;
    };
    // Departure from the paper's exact Eq. (7).
    auto cost = [&](double eps) {
        double worst = 0.0;
        for (int i = 0; i <= 300; ++i) {
            const double cd = 0.70 + 0.30 * i / 300.0;
            for (const double s : {0.6, 0.7, 0.8, 0.9, 0.99}) {
                worst = std::max(worst,
                                 std::abs(ms::expansion_factor(cd, s, g, eps)
                                          - ms::expansion_factor(cd, s, g, 0.0)));
            }
        }
        return worst;
    };

    // The default satisfies both.
    EXPECT_NEAR(continuity(ms::x_smooth_eps), 0.5, 0.05);
    EXPECT_LT(cost(ms::x_smooth_eps), 0.006);

    // Below the window, continuity goes -- the transition becomes narrower
    // than anything a Newton step can resolve.
    EXPECT_GT(std::abs(continuity(0.02) - 0.5), 0.05)
        << "the lower wall has moved; re-measure the window";

    // Above it, fidelity goes.
    EXPECT_GT(cost(0.15), 0.006)
        << "the upper wall has moved; re-measure the window";

    // And the window really is narrow, so the default is not free to drift.
    EXPECT_LT(cost(ms::x_smooth_eps) / 0.006, 0.9);
}

TEST(McGreehanSchotschY, SaturatesAtChoking) {
    const double g = 1.4;
    const double s_crit = ms::critical_pressure_ratio(g);
    EXPECT_NEAR(s_crit, 0.5283, 1e-3);

    auto flux = [&](double S) {
        return ms::expansion_factor(0.99, S, g) * std::sqrt(std::max(1.0 - S, 0.0));
    };
    // The proxy mass flux must not fall away below the choke point.
    const double at_crit = flux(s_crit);
    for (const double S : {0.45, 0.35, 0.20, 0.05}) {
        EXPECT_GT(flux(S), 0.97 * at_crit) << "flow collapses below S* at S = " << S;
    }
}

// Both bounds on the pressure ratio must be soft, not just the choke point.
// A scan for zero-derivative regions and derivative discontinuities caught an
// earlier revision that smoothed S* but left a hard min(S, 1): it put a kink
// of ~0.455 in dY/dS at S = 1 and a floor above it. S > 1 is reverse flow,
// which a solver iterate can reach.
TEST(McGreehanSchotschY, PressureRatioIsSmoothAtBothBounds) {
    const double g = 1.4, cd = 0.90, h = 1e-7;
    auto dS = [&](double S) {
        return (ms::expansion_factor(cd, S + h, g)
                - ms::expansion_factor(cd, S - h, g)) / (2.0 * h);
    };
    auto continuous_at = [&](double S0) {
        auto jump = [&](double half) {
            return std::abs(dS(S0 + half) - dS(S0 - half));
        };
        const double wide = jump(0.004);
        if (wide < 1e-9) return true;   // already flat on both sides
        return jump(0.002) / wide < 0.75;  // shrinks with the interval
    };

    EXPECT_TRUE(continuous_at(1.0)) << "kink at the no-flow bound S = 1";
    EXPECT_TRUE(continuous_at(ms::critical_pressure_ratio(g)))
        << "kink at the choke point";

    // No floor anywhere a solver can wander, including into reverse flow.
    for (const double S : {0.2, 0.4, ms::critical_pressure_ratio(g), 0.7,
                           0.95, 0.999, 1.0, 1.05, 1.2}) {
        EXPECT_NE(dS(S), 0.0) << "derivative floored at S = " << S;
    }
}

TEST(McGreehanSchotschY, BlendSelectsOrificeWhenSharpAndNozzleWhenRounded) {
    const double g = 1.4, S = 0.7;
    // A sharp hole (Cd well under 0.82) is pure orifice.
    EXPECT_NEAR(ms::expansion_factor(0.70, S, g),
                ms::expansion_orifice(S, g), 5e-3);
    // One rounded enough to act like a nozzle is pure nozzle.
    EXPECT_NEAR(ms::expansion_factor(0.995, S, g),
                ms::expansion_nozzle(S, g), 5e-3);
    // They genuinely differ, so the blend is doing work.
    EXPECT_GT(std::abs(ms::expansion_orifice(S, g) - ms::expansion_nozzle(S, g)),
              0.08);
}

TEST(McGreehanSchotschY, StaysFiniteAndBounded) {
    for (const double g : {1.1, 1.4, 1.67}) {
        for (const double S : {1.5, 1.0, 0.5, 0.1, 0.0, -0.1}) {
            for (const double cd : {0.3, 0.8, 1.2}) {
                const double y = ms::expansion_factor(cd, S, g);
                EXPECT_TRUE(std::isfinite(y)) << "g=" << g << " S=" << S;
                EXPECT_GT(y, 0.0);
                EXPECT_LE(y, 1.01);
            }
        }
    }
    EXPECT_DOUBLE_EQ(ms::expansion_factor(0.8, 0.7, 1.0), 1.0);  // degenerate gamma
}

// The crossflow regularisation width, bounded from both sides by measurement
// rather than chosen. Mirrors DefaultSmoothingSitsInsideItsAdmissibleWindow
// for the Y blend, but the hazard is the opposite one: Eq. (17) has an
// UNBOUNDED derivative at U1/Vi = 0, not a dead zone, so the lower wall is a
// cap on |dCd/du| rather than a floor on it.
TEST(McGreehanSchotsch, CrossflowSmoothingSitsInsideItsAdmissibleWindow) {
    const double base = ms::cd(3.2e4, 0.0, 1.0, 0.0);

    // Worst |dCd/du| anywhere. Sampled on a LOG grid: the singularity lives
    // below u = 1e-3 and a linear sweep walks straight past it, which is how
    // an early version of this analysis reported a plateau that is not there.
    auto max_slope = [&](double eps, double lo) {
        const int per_decade = 80;
        const double hi = 10.0;
        const int n = static_cast<int>(std::log10(hi / lo) * per_decade);
        double worst = 0.0;
        for (int i = 0; i <= n; ++i) {
            const double u = lo * std::pow(hi / lo, static_cast<double>(i) / n);
            const double h = u * 1e-4;
            worst = std::max(worst,
                             std::abs((ms::cd_with_crossflow(base, u + h, eps)
                                       - ms::cd_with_crossflow(base, u - h, eps))
                                      / (2.0 * h)));
        }
        return worst;
    };
    // Departure from Eq. (17) exactly, over the range Figs. 4-6 carry data.
    auto cost = [&](double eps) {
        double worst = 0.0;
        for (int i = 0; i <= 240; ++i) {
            const double u = 0.01 * std::pow(1000.0, static_cast<double>(i) / 240.0);
            worst = std::max(worst, std::abs(ms::cd_with_crossflow(base, u, eps)
                                             - ms::cd_with_crossflow(base, u, 0.0)));
        }
        return worst;
    };

    // The physical derivative scale, from the exact correlation over the
    // validated range. Everything below is measured against this, not against
    // a round number.
    const double physical = max_slope(0.0, 0.01);
    EXPECT_NEAR(physical, 0.663, 0.02) << "the physical scale moved; re-measure the window";

    // The data scatter the correlation sits in. The paper states no error
    // statistic at all -- this is digitised from its own Fig. 4.
    const double scatter = 0.02;

    // The default satisfies both walls.
    EXPECT_LT(max_slope(ms::rv_smooth_eps, 1e-9), 10.0 * physical);
    EXPECT_LT(cost(ms::rv_smooth_eps), 0.1 * scatter);

    // Below the window the Jacobian entry goes stiff for no physical reason.
    EXPECT_GT(max_slope(1e-6, 1e-9), 10.0 * physical)
        << "the lower wall has moved; re-measure the window";

    // Above it, the residual stops matching the physics -- the ghost-residual
    // stall that over-smoothing causes, which is the worse failure of the two.
    EXPECT_GT(cost(3e-2), 0.1 * scatter)
        << "the upper wall has moved; re-measure the window";

    // And the default sits towards the LOW-smoothing end deliberately, so it
    // cannot drift upward into the region where ghost residuals appear.
    EXPECT_LT(cost(ms::rv_smooth_eps) / (0.1 * scatter), 0.30);

    // eps = 0 still recovers the paper exactly.
    EXPECT_DOUBLE_EQ(ms::cd_with_crossflow(base, 0.0, 0.0), base);
}

// The (f, J) rule: an analytic derivative is cross-checked against finite
// differences, never shipped on the strength of the derivation alone.
TEST(McGreehanSchotsch, AnalyticDerivativesAgreeWithFiniteDifferences) {
    struct Case { double Re, rd, ld, u; };
    const Case cases[] = {
        {3.2e4, 0.00, 1.0, 0.05}, {3.2e4, 0.00, 1.0, 0.50},
        {1.0e5, 0.10, 2.0, 0.20}, {2.0e4, 0.20, 3.0, 1.00},
        {5.0e4, 0.05, 0.5, 2.00}, {1.5e4, 0.00, 5.0, 0.01},
        {8.0e5, 0.15, 1.5, 0.30},
    };
    for (const auto& c : cases) {
        const auto [cd, dRe, du] = ms::cd_and_derivatives(c.Re, c.rd, c.ld, c.u);

        const double hRe = c.Re * 1e-6;
        const double fdRe = (std::get<0>(ms::cd_and_derivatives(c.Re + hRe, c.rd, c.ld, c.u))
                             - std::get<0>(ms::cd_and_derivatives(c.Re - hRe, c.rd, c.ld, c.u)))
                            / (2.0 * hRe);
        EXPECT_NEAR(dRe, fdRe, 1e-5 * std::max(std::abs(fdRe), 1e-12))
            << "dCd/dRe disagrees with FD at Re = " << c.Re;

        const double hu = c.u * 1e-6;
        const double fdu = (std::get<0>(ms::cd_and_derivatives(c.Re, c.rd, c.ld, c.u + hu))
                            - std::get<0>(ms::cd_and_derivatives(c.Re, c.rd, c.ld, c.u - hu)))
                           / (2.0 * hu);
        EXPECT_NEAR(du, fdu, 1e-5 * std::max(std::abs(fdu), 1e-12))
            << "dCd/d(U1/Vi) disagrees with FD at U1/Vi = " << c.u;

        // The value must be the same correlation, not a reimplementation of
        // it: same expression through the dual, so a few ULP but no more.
        EXPECT_NEAR(cd, ms::cd(c.Re, c.rd, c.ld, c.u), 1e-14);
    }

    // u = 0 is deliberately excluded above: max(U1_over_Vi, 0) makes a
    // central difference one-sided there, so FD measures a forward slope and
    // cannot be compared. The analytic value is zero because the regularised
    // input sqrt(u^2 + eps^2) has zero slope at u = 0 -- a smooth stationary
    // point, not a floor.
    const auto [cd0, dRe0, du0] = ms::cd_and_derivatives(3.2e4, 0.0, 1.0, 0.0);
    EXPECT_DOUBLE_EQ(du0, 0.0);
    EXPECT_LT(dRe0, 0.0) << "Cd must still fall with Re at zero crossflow";
    EXPECT_GT(cd0, 0.0);
}

// The Jacobian-only floor continuation: below re_min the VALUE is held, but
// the derivative is continued from the floor rather than reported as zero.
TEST(McGreehanSchotsch, ReynoldsFloorKeepsValueButContinuesTheDerivative) {
    const auto below = ms::cd_and_derivatives(5.0e3, 0.0, 1.0, 0.2);
    const auto at    = ms::cd_and_derivatives(ms::re_min, 0.0, 1.0, 0.2);

    // Value: floored, exactly as cd() reports it. No Cd changes.
    EXPECT_DOUBLE_EQ(std::get<0>(below), std::get<0>(at));
    EXPECT_NEAR(std::get<0>(below), ms::cd(5.0e3, 0.0, 1.0, 0.2), 1e-14);

    // Derivative: NOT the true zero, but the live slope at the floor, so a
    // Newton step that wanders below has something to climb back on and
    // dCd/dRe is continuous across re_min rather than jumping.
    EXPECT_LT(std::get<1>(below), 0.0) << "the dead zone is back";
    EXPECT_DOUBLE_EQ(std::get<1>(below), std::get<1>(at));
}

// The eps walls were bisected at one operating point (r/d = 0, t/d = 1,
// Re = 3.2e4). They must hold across the geometry range, and there is a
// specific reason to doubt it: Eq. (17)'s Rv carries (cd_base/0.6)^-3, so a
// rounded long hole at Cd ~ 0.95 sees the same U1/Vi as a 4.5x smaller Rv --
// and the regularisation is applied to U1/Vi, not Rv, so its effect in
// Rv-space is geometry-dependent.
//
// Measured, the ratio to the physical scale is nearly geometry-INVARIANT
// (7.8-8.3x across the range) because the Rv stretching moves the regularised
// peak and the physical maximum together. That is why the lower wall is
// robust rather than a coincidence of where it was measured.
TEST(McGreehanSchotsch, CrossflowSmoothingWallsHoldAcrossTheGeometryRange) {
    auto sweep = [](double Re, double rd, double ld, double eps, double lo) {
        const double hi = 10.0;
        const int n = 360;
        double worst = 0.0;
        for (int i = 0; i <= n; ++i) {
            const double u = lo * std::pow(hi / lo, static_cast<double>(i) / n);
            const double h = u * 1e-4;
            worst = std::max(worst, std::abs((ms::cd(Re, rd, ld, u + h, eps)
                                              - ms::cd(Re, rd, ld, u - h, eps))
                                             / (2.0 * h)));
        }
        return worst;
    };
    auto cost = [](double Re, double rd, double ld, double eps) {
        double worst = std::abs(ms::cd(Re, rd, ld, 0.0, eps)
                                - ms::cd(Re, rd, ld, 0.0, 0.0));
        for (int i = 0; i <= 240; ++i) {
            const double u = 1e-9 * std::pow(1e10, static_cast<double>(i) / 240.0);
            worst = std::max(worst, std::abs(ms::cd(Re, rd, ld, u, eps)
                                             - ms::cd(Re, rd, ld, u, 0.0)));
        }
        return worst;
    };

    struct G { double Re, rd, ld; };
    const G geoms[] = {
        {1.0e4, 0.00, 1.0}, {1.0e4, 0.00, 5.0}, {1.0e4, 0.20, 3.0},
        {3.2e4, 0.00, 1.0}, {3.2e4, 0.10, 2.0}, {3.2e4, 0.20, 0.5},
        {1.0e6, 0.00, 1.0}, {1.0e6, 0.20, 3.0},
    };
    for (const auto& g : geoms) {
        const double phys = sweep(g.Re, g.rd, g.ld, 0.0, 0.01);
        const double got  = sweep(g.Re, g.rd, g.ld, ms::rv_smooth_eps, 1e-9);
        EXPECT_LT(got, 10.0 * phys)
            << "stiffness wall broken at Re=" << g.Re << " r/d=" << g.rd
            << " L/d=" << g.ld;
        EXPECT_LT(cost(g.Re, g.rd, g.ld, ms::rv_smooth_eps), 0.1 * 0.02)
            << "fidelity wall broken at Re=" << g.Re << " r/d=" << g.rd
            << " L/d=" << g.ld;
        // The ratio really is near-invariant; if it stops being so, the
        // single-point bisection is no longer a safe way to set eps.
        EXPECT_LT(got / phys, 9.0);
    }
}

TEST(DischargeCdInRange, FlagsEachSourcesOwnEnvelope) {
    using C = DischargeCdCorrelation;
    // Lichtarowicz: l/d 2 to 10 and Re 10 to 2e4 (#465 jet plates sit at
    // t/d 1-3, so the l/d edge is the one that bites).
    EXPECT_TRUE(discharge_cd_in_range(C::Lichtarowicz1965, wall_hole(1e-3, 3e-3, 0.0, 0.0),
                                      wall_flow(1.0e4)));
    EXPECT_FALSE(discharge_cd_in_range(C::Lichtarowicz1965, wall_hole(1e-3, 1.8e-3, 0.0, 0.0),
                                       wall_flow(1.0e4)));
    // Above Eq. (12)'s fitted 2e4 the source takes Cd constant at Eq. (7)
    // (design recommendation (2)): in range up to the end of its collected
    // data, about 1e5; beyond that the held Eq. (7) is an extrapolation.
    EXPECT_TRUE(discharge_cd_in_range(C::Lichtarowicz1965, wall_hole(1e-3, 3e-3, 0.0, 0.0),
                                      wall_flow(3.0e4)));
    EXPECT_TRUE(discharge_cd_in_range(C::Lichtarowicz1965, wall_hole(1e-3, 3e-3, 0.0, 0.0),
                                      wall_flow(1.0e5)));
    EXPECT_FALSE(discharge_cd_in_range(C::Lichtarowicz1965, wall_hole(1e-3, 3e-3, 0.0, 0.0),
                                       wall_flow(1.0e6)));
    EXPECT_FALSE(discharge_cd_in_range(C::Lichtarowicz1965, wall_hole(1e-3, 3e-3, 0.0, 0.0),
                                       wall_flow(5.0)));
    // Idelchik 4-18a: Re 25 to 1e6, l/d up to 4.
    EXPECT_TRUE(discharge_cd_in_range(C::Idelchik1966Thick, wall_hole(1e-3, 1e-3, 0.0, 0.0),
                                      wall_flow(1.0e4)));
    EXPECT_FALSE(discharge_cd_in_range(C::Idelchik1966Thick, wall_hole(1e-3, 5e-3, 0.0, 0.0),
                                       wall_flow(1.0e4)));
    EXPECT_FALSE(discharge_cd_in_range(C::Idelchik1966Thick, wall_hole(1e-3, 1e-3, 0.0, 0.0),
                                       wall_flow(10.0)));
    EXPECT_FALSE(discharge_cd_in_range(C::Idelchik1966Sharp, wall_hole(1e-3, 1e-3, 0.0, 0.0),
                                       wall_flow(1.0e4)));
    // McGreehan-Schotsch: its stated floor Re = 1e4.
    EXPECT_TRUE(discharge_cd_in_range(C::McGreehanSchotsch1988, wall_hole(1e-3, 1e-3, 0.0, 0.0),
                                      wall_flow(2.0e4)));
    EXPECT_FALSE(discharge_cd_in_range(C::McGreehanSchotsch1988,
                                       wall_hole(1e-3, 1e-3, 0.0, 0.0), wall_flow(5.0e3)));
    // ... and beyond L/d 10, the end of the Lichtarowicz data its Eq. (14)
    // was fitted to (Fig. 3), where g keeps falling.
    EXPECT_FALSE(discharge_cd_in_range(C::McGreehanSchotsch1988,
                                       wall_hole(1e-3, 12e-3, 0.0, 0.0), wall_flow(2.0e4)));
    EXPECT_TRUE(discharge_cd_in_range(C::Constant, wall_hole(1e-3, 9e-3, 0.0, 0.0),
                                      wall_flow(1.0)));
}

TEST(DischargeCdEnvelope, NamesTheVariableOutsideTheSourcesRange) {
    using C = DischargeCdCorrelation;
    auto env = [](C id, double L, double Re) {
        return discharge_cd_envelope(id, wall_hole(1e-3, L, 0.0, 0.0), wall_flow(Re));
    };
    // Idelchik 4-18a beyond its l/d 4 table end: l/d only, not Re.
    auto e = env(C::Idelchik1966Thick, 6e-3, 5e4);
    EXPECT_TRUE(e.l_over_d);
    EXPECT_FALSE(e.re);
    EXPECT_FALSE(e.edge);
    EXPECT_FALSE(e.in_range());
    // Below its Re table: Re only.
    e = env(C::Idelchik1966Thick, 2e-3, 10.0);
    EXPECT_TRUE(e.re);
    EXPECT_FALSE(e.l_over_d);
    // McGreehan-Schotsch: Re below its floor (value held), l/d beyond 10.
    e = env(C::McGreehanSchotsch1988, 1e-3, 5e3);
    EXPECT_TRUE(e.re);
    EXPECT_FALSE(e.l_over_d);
    e = env(C::McGreehanSchotsch1988, 40e-3, 5e4);
    EXPECT_FALSE(e.re);
    EXPECT_TRUE(e.l_over_d);
    // Lichtarowicz: l/d beyond 10 (held at 10); Re below 10 or above 1e5;
    // 2e4-1e5 is Eq. (7), in range.
    e = env(C::Lichtarowicz1965, 12e-3, 5e4);
    EXPECT_TRUE(e.l_over_d);
    EXPECT_FALSE(e.re);
    EXPECT_TRUE(env(C::Lichtarowicz1965, 3e-3, 5.0).re);
    e = env(C::Lichtarowicz1965, 3e-3, 2e5);
    EXPECT_TRUE(e.re);
    EXPECT_FALSE(e.l_over_d);
    // Rounded beyond its r/d table: the edge axis.
    e = discharge_cd_envelope(C::Idelchik1966Rounded, wall_hole(1e-3, 0.0, 0.5e-3, 0.0),
                              wall_flow(1e5));
    EXPECT_TRUE(e.edge);
    EXPECT_FALSE(e.re);
    // in_range is the envelope's and nothing else.
    for (C id : {C::Idelchik1966Thick, C::McGreehanSchotsch1988, C::Lichtarowicz1965}) {
        for (double L : {1e-3, 3e-3, 6e-3, 12e-3, 40e-3}) {
            for (double Re : {5.0, 50.0, 5e3, 5e4, 2e6}) {
                const DischargeHoleGeometry h = wall_hole(1e-3, L, 0.0, 0.0);
                EXPECT_EQ(discharge_cd_in_range(id, h, wall_flow(Re)),
                          discharge_cd_envelope(id, h, wall_flow(Re)).in_range());
            }
        }
    }
}

// -------------------------------------------------------------
// ISO 5167-2 permanent pressure loss (#645)
// -------------------------------------------------------------

TEST(PermanentLoss, RatioIsIsosLongFormWithAnExactCDerivative) {
    for (double b : {0.1, 0.3, 0.5, 0.7, 0.75}) {
        for (double C : {0.6, 0.62, 0.8}) {
            auto k_of = [b](double c) {
                const double s = std::sqrt(1.0 - std::pow(b, 4) * (1.0 - c * c));
                return (s - c * b * b) / (s + c * b * b);
            };
            auto [k, dk] = permanent_loss_ratio_and_dC(b, C);
            EXPECT_NEAR(k, k_of(C), 1e-14) << b << " " << C;
            const double h = 1e-6;
            EXPECT_NEAR(dk, (k_of(C + h) - k_of(C - h)) / (2.0 * h), 1e-8) << b << " " << C;
            EXPECT_GT(k, 0.0);
            EXPECT_LT(k, 1.0);
        }
    }
}

TEST(PermanentLoss, RefusesOutsideItsDomain) {
    EXPECT_THROW(permanent_loss_ratio_and_dC(0.0, 0.6), std::invalid_argument);
    EXPECT_THROW(permanent_loss_ratio_and_dC(1.0, 0.6), std::invalid_argument);
    EXPECT_THROW(permanent_loss_ratio_and_dC(0.5, 0.0), std::invalid_argument);
}

TEST(VelocityOfApproach, IsIsoFactorInsideAndOneOutside) {
    const double b = 0.6;
    EXPECT_NEAR(velocity_of_approach_factor(b), 1.0 / std::sqrt(1.0 - b * b * b * b), 1e-15);
    // No approach pipe (beta = 0), or a bore at least as wide as the pipe:
    // the network applies no correction rather than NaN or inf.
    EXPECT_EQ(velocity_of_approach_factor(0.0), 1.0);
    EXPECT_EQ(velocity_of_approach_factor(1.0), 1.0);
    EXPECT_EQ(velocity_of_approach_factor(1.5), 1.0);
    EXPECT_EQ(velocity_of_approach_factor(-0.2), 1.0);
}

namespace {

DischargeHoleGeometry pipe_hole(double L_d, double r_d, double bevel_d) {
    DischargeHoleGeometry h;
    h.d = 0.01;
    h.L = L_d * h.d;
    h.r = r_d * h.d;
    h.bevel = bevel_d * h.d;
    return h;
}

// Worst relative error of discharge_pipe_zeta against a printed row, n <= 0.8.
double worst_vs_row(DischargeCdCorrelation id, const DischargeHoleGeometry& h,
                    const std::vector<double>& n, const std::vector<double>& zeta) {
    double worst = 0.0;
    for (std::size_t i = 0; i < n.size(); ++i) {
        if (n[i] > 0.8) {
            continue;
        }
        const double e = discharge_pipe_zeta(id, h, std::sqrt(n[i])) / zeta[i] - 1.0;
        worst = std::max(worst, std::abs(e));
    }
    return worst;
}

}  // namespace

// Idelchik (1966) diagrams 4-14..4-16, Re >= 1e4-1e5, read at 400 dpi.
// zeta' and tau come from the wall tables (4-18), not the printed columns
// beside these rows; 4.2% is the measured worst over every row below.
TEST(DischargePipe, ZetaReproducesIdelchiksInPipeTables) {
    using C = DischargeCdCorrelation;
    const std::vector<double> n14 = {0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20, 0.25,
                                     0.30, 0.40, 0.50, 0.60, 0.70, 0.80};
    const std::vector<double> n16 = {0.02, 0.04, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60,
                                     0.70, 0.80};
    struct Row {
        C id;
        DischargeHoleGeometry h;
        const std::vector<double>* n;
        std::vector<double> zeta;
    };
    const std::vector<Row> rows = {
        {C::Idelchik1966Thick, pipe_hole(0.0, 0.0, 0.0), &n14,
         {7000, 1670, 730, 400, 245, 96.0, 51.5, 30.0, 18.2, 8.25, 4.00, 2.00, 0.97, 0.42}},
        {C::Idelchik1966Thick, pipe_hole(1.0, 0.0, 0.0), &n14,
         {4260, 1030, 443, 240, 149, 60.0, 31.0, 17.8, 11.1, 5.00, 2.40, 1.20, 0.61, 0.29}},
        {C::Idelchik1966Thick, pipe_hole(3.0, 0.0, 0.0), &n14,
         {3765, 913, 392, 214, 132, 53.5, 27.5, 15.9, 10.0, 4.50, 2.24, 1.17, 0.61, 0.31}},
        {C::Idelchik1966Beveled, pipe_hole(0.0, 0.0, 0.04), &n14,
         {6130, 1480, 636, 345, 214, 86.5, 44.2, 25.6, 15.8, 7.00, 3.36, 1.67, 0.80, 0.34}},
        {C::Idelchik1966Beveled, pipe_hole(0.0, 0.0, 0.16), &n14,
         {4460, 1080, 462, 251, 154, 62.0, 31.6, 18.1, 11.2, 4.80, 2.28, 1.10, 0.50, 0.20}},
        {C::Idelchik1966Rounded, pipe_hole(0.0, 0.04, 0.0), &n16,
         {5510, 1330, 192, 39.0, 14.1, 6.19, 3.00, 1.45, 0.70, 0.29}},
        {C::Idelchik1966Rounded, pipe_hole(0.0, 0.16, 0.0), &n16,
         {3320, 797, 113, 23.0, 7.90, 3.40, 1.60, 0.70, 0.32, 0.12}},
    };
    for (const Row& row : rows) {
        EXPECT_LT(worst_vs_row(row.id, row.h, *row.n, row.zeta), 0.042);
    }
}

TEST(DischargePipe, RatioCarriesTheLossThroughIsosRecoveryRelation) {
    // C_pipe = R c0, c0 = 1/sqrt(zeta at n = 0): ISO's permanent loss with
    // that C is the in-pipe zeta (the network's recovery row then charges it).
    using C = DischargeCdCorrelation;
    for (C id : {C::Idelchik1966Sharp, C::Idelchik1966Thick, C::Idelchik1966Rounded,
                 C::Lichtarowicz1965}) {
        const DischargeHoleGeometry h = pipe_hole(2.0, 0.04, 0.0);
        const double tiny = 1e-4;
        // zeta1 n^2 at n = 1e-8: the hole-velocity loss of the wall limit.
        const double z00 = tiny * tiny * tiny * tiny * discharge_pipe_zeta(id, h, tiny);
        for (double beta : {0.2, 0.5, 0.7, 0.85}) {
            const double c = discharge_pipe_cd_ratio(id, h, beta) / std::sqrt(z00);
            const double k = std::get<0>(permanent_loss_ratio_and_dC(beta, c));
            const double b4 = std::pow(beta, 4.0);
            EXPECT_NEAR(k * (1.0 - b4) / (c * c * b4) / discharge_pipe_zeta(id, h, beta), 1.0,
                        1e-6);
        }
    }
}

TEST(DischargePipe, RatioIsOneWithoutAPipe) {
    using C = DischargeCdCorrelation;
    const DischargeHoleGeometry h = pipe_hole(1.0, 0.0, 0.0);
    EXPECT_EQ(discharge_pipe_cd_ratio(C::Idelchik1966Thick, h, 0.0), 1.0);
    EXPECT_EQ(discharge_pipe_cd_ratio(C::Idelchik1966Thick, h, 1.0), 1.0);
    EXPECT_EQ(discharge_pipe_cd_ratio(C::Constant, h, 0.5), 1.0);
    // Continuous with the plenum Cd as the pipe grows.
    EXPECT_NEAR(discharge_pipe_cd_ratio(C::Idelchik1966Thick, h, 1e-3), 1.0, 1e-6);
    EXPECT_THROW(discharge_pipe_zeta(C::Idelchik1966Thick, h, 0.0), std::invalid_argument);
    EXPECT_THROW(discharge_pipe_zeta(C::Constant, h, 0.5), std::invalid_argument);
}

// -------------------------------------------------------------
// Opt-in hole inclination: Idelchik diagram 3-2 in Eq. (4-20)'s zeta'
// -------------------------------------------------------------

namespace {

DischargeHoleGeometry inclined_hole(double l_over_d, double delta_deg) {
    DischargeHoleGeometry hole = wall_hole(1e-3, l_over_d * 1e-3, 0.0, 0.0);
    hole.inclination_deg = delta_deg;
    return hole;
}

}  // namespace

TEST(HoleInclination, IncrementReproducesDiagram32Table) {
    // Diagram 3-2's printed table (PDF p. 99): zeta_in at delta 20..90 deg.
    const std::vector<std::pair<double, double>> table = {
        {20.0, 0.96}, {30.0, 0.91}, {45.0, 0.81}, {60.0, 0.70},
        {70.0, 0.63}, {80.0, 0.56}, {90.0, 0.50}};
    for (const auto& [delta, zeta_in] : table) {
        EXPECT_NEAR(orifice::idelchik::pipe_thick_inlet + discharge_inclination_dzeta(delta),
                    zeta_in, 0.005)
            << delta;
    }
    EXPECT_EQ(discharge_inclination_dzeta(90.0), 0.0);  // exactly
    EXPECT_NEAR(discharge_inclination_dzeta(30.0), 0.3 * std::sqrt(3.0) / 2.0 + 0.2 * 0.75,
                1e-15);
    EXPECT_THROW(discharge_inclination_dzeta(0.0), std::invalid_argument);
    EXPECT_THROW(discharge_inclination_dzeta(91.0), std::invalid_argument);
}

TEST(HoleInclination, NormalHoleIsBitIdenticalToTheDefault) {
    auto c = make_discharge_correlation(DischargeCdCorrelation::Idelchik1966Thick);
    for (double re : {50.0, 3e3, 2e4, 1e6}) {
        const auto a = c->Cd_and_derivatives(wall_hole(1e-3, 3e-3, 0.0, 0.0), wall_flow(re));
        const auto b = c->Cd_and_derivatives(inclined_hole(3.0, 90.0), wall_flow(re));
        EXPECT_EQ(std::get<0>(a), std::get<0>(b));
        EXPECT_EQ(std::get<1>(a), std::get<1>(b));
    }
}

TEST(HoleInclination, ThirtyDegreesAddsTheIncrementToZetaPrimeAtHighRe) {
    // At the top of the table (k eps_re = 1, zeta_phi0 = 0) the composed
    // zeta is the normal hole's plus dzeta exactly.
    auto c = make_discharge_correlation(DischargeCdCorrelation::Idelchik1966Thick);
    const double re = 1e6;
    const double cd0 = c->Cd(inclined_hole(3.0, 90.0), wall_flow(re));
    const double cd30 = c->Cd(inclined_hole(3.0, 30.0), wall_flow(re));
    const double dz = discharge_inclination_dzeta(30.0);
    EXPECT_NEAR(1.0 / (cd30 * cd30), 1.0 / (cd0 * cd0) + dz, 1e-12);
    EXPECT_LT(cd30, cd0);
}

TEST(HoleInclination, DerivativeContractIsExact) {
    // (Cd, dCd/dRe, dCd/dU) analytic: dCd/dRe against a central difference,
    // and the crossflow slot exactly zero.
    auto c = make_discharge_correlation(DischargeCdCorrelation::Idelchik1966Thick);
    const DischargeHoleGeometry h = inclined_hole(3.0, 30.0);
    for (double re : {300.0, 3e3, 3e4, 3e5}) {
        const auto [cd, dcd, du] = c->Cd_and_derivatives(h, wall_flow(re));
        const double step = 1e-6 * re;
        const double fd =
            (c->Cd(h, wall_flow(re + step)) - c->Cd(h, wall_flow(re - step))) / (2.0 * step);
        EXPECT_NEAR(dcd, fd, 1e-6 * std::abs(cd) / re + 1e-12) << re;
        EXPECT_EQ(du, 0.0);
    }
}

TEST(HoleInclination, OnlyTheThickHoleComposes) {
    using C = DischargeCdCorrelation;
    const DischargeHoleGeometry h = inclined_hole(3.0, 30.0);
    for (C id : {C::Idelchik1966Sharp, C::Idelchik1966Beveled, C::Idelchik1966Rounded,
                 C::Lichtarowicz1965, C::McGreehanSchotsch1988, C::Constant}) {
        auto c = make_discharge_correlation(id);
        EXPECT_THROW(c->Cd(h, wall_flow(2e4)), std::invalid_argument);
        EXPECT_THROW(c->Cd_and_derivatives(h, wall_flow(2e4)), std::invalid_argument);
        EXPECT_FALSE(discharge_cd_in_range(id, h, wall_flow(2e4)));
    }
    EXPECT_THROW(make_constant_discharge_correlation(0.7)->Cd(h, wall_flow(2e4)),
                 std::invalid_argument);
    EXPECT_THROW(discharge_pipe_zeta(C::Lichtarowicz1965, h, 0.5), std::invalid_argument);
}

TEST(HoleInclination, RangeFlagsDiagram32sReAndTheReattachmentFloor) {
    using C = DischargeCdCorrelation;
    EXPECT_TRUE(
        discharge_cd_in_range(C::Idelchik1966Thick, inclined_hole(3.0, 30.0), wall_flow(2e4)));
    // Below diagram 3-2's Re > 1e4 and below l/d 2: flagged, still evaluated.
    EXPECT_FALSE(
        discharge_cd_in_range(C::Idelchik1966Thick, inclined_hole(3.0, 30.0), wall_flow(5e3)));
    EXPECT_FALSE(
        discharge_cd_in_range(C::Idelchik1966Thick, inclined_hole(1.5, 30.0), wall_flow(2e4)));
    // The same short hole, normal, is inside the plain Thick range.
    EXPECT_TRUE(
        discharge_cd_in_range(C::Idelchik1966Thick, inclined_hole(1.5, 90.0), wall_flow(2e4)));
    auto c = make_discharge_correlation(C::Idelchik1966Thick);
    EXPECT_GT(c->Cd(inclined_hole(1.5, 30.0), wall_flow(2e4)), 0.0);
}

TEST(HoleInclination, InPipeLossTakesTheSameInletCoefficient) {
    using C = DischargeCdCorrelation;
    const double beta = 1e-3;  // the wall limit: zeta1 n^2 -> the wall zeta'
    const double n2 = std::pow(beta, 4.0);
    const double z90 =
        n2 * discharge_pipe_zeta(C::Idelchik1966Thick, inclined_hole(3.0, 90.0), beta);
    const double z30 =
        n2 * discharge_pipe_zeta(C::Idelchik1966Thick, inclined_hole(3.0, 30.0), beta);
    EXPECT_NEAR(z30 - z90, discharge_inclination_dzeta(30.0), 1e-3);
}
