// Turbulent flat plate -- NASA CR-2648's unblown reference St_o.
//
// Ground truth: the printed St_o Pr^0.4 = 0.0295 Re_x^-0.2, and the panel
// average as an independent numerical integral of the local value, not this
// code's own closed form.

#include "cooling_correlations.h"

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <tuple>

namespace cc = combaero::cooling;

TEST(FlatPlateTurbulent, LocalIsThePrintedStantonNumber) {
    const double Re = 1.0e6;
    const double Pr = 0.71;
    const double st0 = 0.0295 * std::pow(Re, -0.2) * std::pow(Pr, -0.4);  // CR-2648
    const auto [nu, dnu] = cc::flat_plate_turbulent_nusselt_and_derivative(Re, Pr);
    EXPECT_NEAR(nu / (st0 * Re * Pr), 1.0, 1e-14);
    const double h = 1e-6 * Re;
    const double fd = (std::get<0>(cc::flat_plate_turbulent_nusselt_and_derivative(Re + h, Pr))
                       - std::get<0>(cc::flat_plate_turbulent_nusselt_and_derivative(Re - h, Pr)))
                    / (2.0 * h);
    EXPECT_NEAR(dnu / fd, 1.0, 1e-8);
}

TEST(FlatPlateTurbulent, PanelAverageIsTheIntegralOfTheLocal) {
    // Midpoint rule over x in [x0, x0 + L], h_x = Nu_x k / x, against the
    // closed form. Re_L = 2e5 and x0/L = 0.8 (Re_x 1.6e5 to 3.6e5) put the panel inside the data.
    const double Re_L = 2.0e5;
    const double Pr = 0.7;
    const double s = 0.8;
    const int n = 200000;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        const double xi = s + (i + 0.5) / n;  // x / L
        const double nu_x = std::get<0>(
            cc::flat_plate_turbulent_nusselt_and_derivative(Re_L * xi, Pr));
        sum += nu_x / xi / n;  // (h_x L / k) dx/L
    }
    const auto r = cc::flat_plate_turbulent_panel_nusselt(Re_L, Pr, s);
    EXPECT_NEAR(r.Nu / sum, 1.0, 1e-8);
    EXPECT_DOUBLE_EQ(r.Re_x_start, Re_L * s);
    EXPECT_DOUBLE_EQ(r.Re_x_end, Re_L * (1.0 + s));
    EXPECT_FALSE(r.extrapolated);
}

TEST(FlatPlateTurbulent, FromTheLeadingEdgeIsTheTextbookAverage) {
    // s = 0: (0.0295/0.8) Re_L^0.8 Pr^0.6, i.e. 0.036875 -- the same 0.037
    // the Pr^(1/3) textbook form carries.
    const auto r = cc::flat_plate_turbulent_panel_nusselt(1.0e6, 0.7, 0.0);
    EXPECT_NEAR(r.Nu / (0.036875 * std::pow(1.0e6, 0.8) * std::pow(0.7, 0.6)), 1.0, 1e-14);
    // Re_x starts at 0, below the data: flagged, still evaluated.
    EXPECT_TRUE(r.extrapolated);
}

TEST(FlatPlateTurbulent, DerivativeMatchesFiniteDifference) {
    const double Re_L = 3.0e4;
    const double Pr = 0.72;
    const double s = 1.75;
    const auto [nu, dnu] = cc::flat_plate_turbulent_panel_nusselt_and_derivative(Re_L, Pr, s);
    const double h = 1e-6 * Re_L;
    const double fd =
        (cc::flat_plate_turbulent_panel_nusselt(Re_L + h, Pr, s).Nu
         - cc::flat_plate_turbulent_panel_nusselt(Re_L - h, Pr, s).Nu)
        / (2.0 * h);
    EXPECT_NEAR(dnu / fd, 1.0, 1e-8);
    EXPECT_DOUBLE_EQ(nu, cc::flat_plate_turbulent_panel_nusselt(Re_L, Pr, s).Nu);
}

TEST(FlatPlateTurbulent, EnvelopeIsTheStanfordReSpan) {
    combaero::CorrelationStatus st = combaero::CorrelationStatus::Valid;
    cc::flat_plate_turbulent_nusselt_and_derivative(1.49e5, 0.7, &st);
    EXPECT_EQ(st, combaero::CorrelationStatus::Extrapolated);
    cc::flat_plate_turbulent_nusselt_and_derivative(1.5e5, 0.7, &st);
    EXPECT_EQ(st, combaero::CorrelationStatus::Valid);
    cc::flat_plate_turbulent_nusselt_and_derivative(5.01e6, 0.7, &st);
    EXPECT_EQ(st, combaero::CorrelationStatus::Extrapolated);
    // A panel whose trailing edge leaves the span is flagged too.
    EXPECT_TRUE(cc::flat_plate_turbulent_panel_nusselt(3.0e6, 0.7, 1.0).extrapolated);
    cc::flat_plate_turbulent_panel_nusselt_and_derivative(2.0e5, 0.7, 0.75, &st);
    EXPECT_EQ(st, combaero::CorrelationStatus::Valid);
}

TEST(FlatPlateTurbulent, RefusesNonPhysicalInputs) {
    EXPECT_THROW(cc::flat_plate_turbulent_nusselt_and_derivative(0.0, 0.7),
                 std::invalid_argument);
    EXPECT_THROW(cc::flat_plate_turbulent_nusselt_and_derivative(1e5, 0.0),
                 std::invalid_argument);
    EXPECT_THROW(cc::flat_plate_turbulent_panel_nusselt(1e5, 0.7, -0.1),
                 std::invalid_argument);
    EXPECT_THROW(cc::flat_plate_turbulent_panel_nusselt(-1.0, 0.7, 0.0),
                 std::invalid_argument);
}
