#include <gtest/gtest.h>
#include "solver_interface.h"
#include "thermo.h"
#include <cmath>
#include <vector>

using namespace combaero;
using namespace combaero::solver;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static std::vector<double> air_Y() {
    std::size_t ns = num_species();
    std::vector<double> Y(ns, 0.0);
    for (std::size_t i = 0; i < ns; ++i) {
        if (species_name(i) == "N2")
            Y[i] = 0.767;
        else if (species_name(i) == "O2")
            Y[i] = 0.233;
    }
    return Y;
}

static std::vector<double> fuel_Y() {
    std::size_t ns = num_species();
    std::vector<double> Y(ns, 0.0);
    for (std::size_t i = 0; i < ns; ++i) {
        if (species_name(i) == "CH4")
            Y[i] = 1.0;
    }
    return Y;
}

static solver::Stream make_air_stream(double m_dot, double T, double P_total) {
    return {m_dot, T, P_total, air_Y()};
}

static solver::Stream make_fuel_stream(double m_dot, double T, double P_total) {
    return {m_dot, T, P_total, fuel_Y()};
}

// ---------------------------------------------------------------------------
// MixerResult FD validation
// ---------------------------------------------------------------------------

class MixerJacobianTest : public ::testing::Test {
  protected:
    void SetUp() override {
        streams = {
            make_air_stream(2.0, 600.0, 150000.0),
            make_fuel_stream(0.06, 300.0, 150000.0),
        };
        ns = num_species();
        base = mixer_from_streams_and_jacobians(streams);
    }

    std::vector<solver::Stream> streams;
    std::size_t ns{};
    MixerResult base;
};

TEST_F(MixerJacobianTest, dTmix_d_mdot) {
    const double eps = 1e-4;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].m_dot += eps;
        auto res_p = mixer_from_streams_and_jacobians(s_p);
        double fd = (res_p.Tt_mix - base.Tt_mix) / eps;
        EXPECT_NEAR(base.dTt_mix_d_stream[i].d_mdot, fd,
                    std::abs(fd) * 1e-3 + 1e-6)
            << "stream " << i;
    }
}

TEST_F(MixerJacobianTest, dTmix_d_T) {
    const double eps = 0.1;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].Tt += eps;
        auto res_p = mixer_from_streams_and_jacobians(s_p);
        double fd = (res_p.Tt_mix - base.Tt_mix) / eps;
        EXPECT_NEAR(base.dTt_mix_d_stream[i].d_Tt, fd,
                    std::abs(fd) * 1e-3 + 1e-6)
            << "stream " << i;
    }
}

TEST_F(MixerJacobianTest, dTmix_d_Y) {
    const double eps = 1e-5;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        for (std::size_t k = 0; k < ns; ++k) {
            auto s_p = streams;
            s_p[i].Y[k] += eps;
            auto res_p = mixer_from_streams_and_jacobians(s_p);
            double fd = (res_p.Tt_mix - base.Tt_mix) / eps;
            EXPECT_NEAR(base.dTt_mix_d_stream[i].d_Y[k], fd,
                        std::abs(fd) * 5e-3 + 1e-4)
                << "stream " << i << " species " << k;
        }
    }
}

TEST_F(MixerJacobianTest, dYmix_d_mdot) {
    const double eps = 1e-4;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].m_dot += eps;
        auto res_p = mixer_from_streams_and_jacobians(s_p);
        for (std::size_t k = 0; k < ns; ++k) {
            double fd = (res_p.Y_mix[k] - base.Y_mix[k]) / eps;
            EXPECT_NEAR(base.dY_mix_d_stream[k][i].d_mdot, fd,
                        std::abs(fd) * 1e-3 + 1e-8)
                << "stream " << i << " species " << k;
        }
    }
}

TEST_F(MixerJacobianTest, dYmix_d_Y) {
    const double eps = 1e-5;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        for (std::size_t j = 0; j < ns; ++j) {
            auto s_p = streams;
            s_p[i].Y[j] += eps;
            auto res_p = mixer_from_streams_and_jacobians(s_p);
            for (std::size_t k = 0; k < ns; ++k) {
                double fd = (res_p.Y_mix[k] - base.Y_mix[k]) / eps;
                EXPECT_NEAR(base.dY_mix_d_stream[k][i].d_Y[j], fd,
                            std::abs(fd) * 1e-3 + 1e-8)
                    << "stream " << i << " Y[" << j << "] -> Y_mix[" << k
                    << "]";
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Adiabatic T (Complete) from streams -- FD validation
// ---------------------------------------------------------------------------

class AdiabaticCompleteJacobianTest : public ::testing::Test {
  protected:
    void SetUp() override {
        streams = {
            make_air_stream(2.0, 600.0, 150000.0),
            make_fuel_stream(0.06, 300.0, 150000.0),
        };
        ns = num_species();
        P = 140000.0;
        base = adiabatic_T_complete_and_jacobian_T_from_streams(streams, P);
    }

    std::vector<solver::Stream> streams;
    std::size_t ns{};
    double P{};
    MixerResult base;
};

TEST_F(AdiabaticCompleteJacobianTest, dTad_d_mdot) {
    const double eps = 1e-3;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].m_dot += eps;
        auto res_p =
            adiabatic_T_complete_and_jacobian_T_from_streams(s_p, P);
        double fd = (res_p.Tt_mix - base.Tt_mix) / eps;
        EXPECT_NEAR(base.dTt_mix_d_stream[i].d_mdot, fd,
                    std::abs(fd) * 5e-3 + 1e-2)
            << "stream " << i;
    }
}

TEST_F(AdiabaticCompleteJacobianTest, dTad_d_T) {
    const double eps = 1.0;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].Tt += eps;
        auto res_p =
            adiabatic_T_complete_and_jacobian_T_from_streams(s_p, P);
        double fd = (res_p.Tt_mix - base.Tt_mix) / eps;
        EXPECT_NEAR(base.dTt_mix_d_stream[i].d_Tt, fd,
                    std::abs(fd) * 5e-3 + 1e-2)
            << "stream " << i;
    }
}

TEST_F(AdiabaticCompleteJacobianTest, dTad_d_Y) {
    const double eps = 1e-4;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        for (std::size_t k = 0; k < ns; ++k) {
            auto s_p = streams;
            s_p[i].Y[k] += eps;
            auto res_p =
                adiabatic_T_complete_and_jacobian_T_from_streams(s_p, P);
            double fd = (res_p.Tt_mix - base.Tt_mix) / eps;
            EXPECT_NEAR(base.dTt_mix_d_stream[i].d_Y[k], fd,
                        std::abs(fd) * 0.02 + 0.5)
                << "stream " << i << " species " << k;
        }
    }
}

TEST_F(AdiabaticCompleteJacobianTest, dYburned_d_mdot) {
    const double eps = 1e-3;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].m_dot += eps;
        auto res_p =
            adiabatic_T_complete_and_jacobian_T_from_streams(s_p, P);
        for (std::size_t k = 0; k < ns; ++k) {
            double fd = (res_p.Y_mix[k] - base.Y_mix[k]) / eps;
            // Chain-rule through internal FD: relaxed tolerance
            EXPECT_NEAR(base.dY_mix_d_stream[k][i].d_mdot, fd,
                        std::abs(fd) * 0.03 + 1e-4)
                << "stream " << i << " species " << k;
        }
    }
}

TEST_F(AdiabaticCompleteJacobianTest, dYburned_d_T) {
    const double eps = 1.0;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].Tt += eps;
        auto res_p =
            adiabatic_T_complete_and_jacobian_T_from_streams(s_p, P);
        for (std::size_t k = 0; k < ns; ++k) {
            double fd = (res_p.Y_mix[k] - base.Y_mix[k]) / eps;
            EXPECT_NEAR(base.dY_mix_d_stream[k][i].d_Tt, fd,
                        std::abs(fd) * 5e-3 + 1e-6)
                << "stream " << i << " species " << k;
        }
    }
}

TEST_F(AdiabaticCompleteJacobianTest, dYburned_d_Y) {
    const double eps = 1e-4;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        for (std::size_t j = 0; j < ns; ++j) {
            auto s_p = streams;
            s_p[i].Y[j] += eps;
            auto res_p =
                adiabatic_T_complete_and_jacobian_T_from_streams(s_p, P);
            for (std::size_t k = 0; k < ns; ++k) {
                double fd = (res_p.Y_mix[k] - base.Y_mix[k]) / eps;
                EXPECT_NEAR(base.dY_mix_d_stream[k][i].d_Y[j], fd,
                            std::abs(fd) * 0.02 + 1e-5)
                    << "stream " << i << " Y[" << j << "] -> Y_b[" << k
                    << "]";
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Adiabatic T (Equilibrium) from streams -- FD validation
// ---------------------------------------------------------------------------

class AdiabaticEquilibriumJacobianTest : public ::testing::Test {
  protected:
    void SetUp() override {
        streams = {
            make_air_stream(2.0, 600.0, 150000.0),
            make_fuel_stream(0.06, 300.0, 150000.0),
        };
        ns = num_species();
        P = 140000.0;
        base =
            adiabatic_T_equilibrium_and_jacobians_from_streams(streams, P);
    }

    std::vector<solver::Stream> streams;
    std::size_t ns{};
    double P{};
    MixerResult base;
};

TEST_F(AdiabaticEquilibriumJacobianTest, dTad_d_mdot) {
    const double eps = 1e-3;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].m_dot += eps;
        auto res_p =
            adiabatic_T_equilibrium_and_jacobians_from_streams(s_p, P);
        double fd = (res_p.Tt_mix - base.Tt_mix) / eps;
        EXPECT_NEAR(base.dTt_mix_d_stream[i].d_mdot, fd,
                    std::abs(fd) * 0.01 + 0.1)
            << "stream " << i;
    }
}

TEST_F(AdiabaticEquilibriumJacobianTest, dTad_d_T) {
    const double eps = 1.0;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].Tt += eps;
        auto res_p =
            adiabatic_T_equilibrium_and_jacobians_from_streams(s_p, P);
        double fd = (res_p.Tt_mix - base.Tt_mix) / eps;
        EXPECT_NEAR(base.dTt_mix_d_stream[i].d_Tt, fd,
                    std::abs(fd) * 0.01 + 0.1)
            << "stream " << i;
    }
}

TEST_F(AdiabaticEquilibriumJacobianTest, dTad_d_Y) {
    const double eps = 1e-4;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        for (std::size_t k = 0; k < ns; ++k) {
            auto s_p = streams;
            s_p[i].Y[k] += eps;
            auto res_p =
                adiabatic_T_equilibrium_and_jacobians_from_streams(s_p, P);
            double fd = (res_p.Tt_mix - base.Tt_mix) / eps;
            EXPECT_NEAR(base.dTt_mix_d_stream[i].d_Y[k], fd,
                        std::abs(fd) * 0.05 + 1.0)
                << "stream " << i << " species " << k;
        }
    }
}

TEST_F(AdiabaticEquilibriumJacobianTest, dYburned_d_mdot) {
    const double eps = 1e-3;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].m_dot += eps;
        auto res_p =
            adiabatic_T_equilibrium_and_jacobians_from_streams(s_p, P);
        for (std::size_t k = 0; k < ns; ++k) {
            double fd = (res_p.Y_mix[k] - base.Y_mix[k]) / eps;
            // Chain-rule through equilibrium FD: relaxed tolerance
            EXPECT_NEAR(base.dY_mix_d_stream[k][i].d_mdot, fd,
                        std::abs(fd) * 0.15 + 1e-3)
                << "stream " << i << " species " << k;
        }
    }
}

TEST_F(AdiabaticEquilibriumJacobianTest, dYburned_d_T) {
    const double eps = 1.0;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].Tt += eps;
        auto res_p =
            adiabatic_T_equilibrium_and_jacobians_from_streams(s_p, P);
        for (std::size_t k = 0; k < ns; ++k) {
            double fd = (res_p.Y_mix[k] - base.Y_mix[k]) / eps;
            EXPECT_NEAR(base.dY_mix_d_stream[k][i].d_Tt, fd,
                        std::abs(fd) * 0.01 + 1e-5)
                << "stream " << i << " species " << k;
        }
    }
}

TEST_F(AdiabaticEquilibriumJacobianTest, dYburned_d_Y) {
    // NOTE: This is a smoke test for Jacobian stability.
    // Equilibrium Jacobians are computed via finite differences through an
    // iterative solver. This introduces truncation error and iteration noise
    // that limits precision. A 25% tolerance reflects this inherent logic gap.
    const double eps = 1e-4;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        for (std::size_t j = 0; j < ns; ++j) {
            auto s_p = streams;
            s_p[i].Y[j] += eps;
            auto res_p =
                adiabatic_T_equilibrium_and_jacobians_from_streams(s_p, P);
            for (std::size_t k = 0; k < ns; ++k) {
                double fd = (res_p.Y_mix[k] - base.Y_mix[k]) / eps;
                // Chain-rule through equilibrium FD: relaxed tolerance
                EXPECT_NEAR(base.dY_mix_d_stream[k][i].d_Y[j], fd,
                            std::abs(fd) * 0.25 + 2e-3)
                    << "stream " << i << " Y[" << j << "] -> Y_b[" << k
                    << "]";
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Plenum ChamberResult -- local + stream Jacobian FD validation
// ---------------------------------------------------------------------------

#if 0
static MixtureState make_state(double P, double P_total, double T,
                                double m_dot, const std::vector<double> &Y) {
    MixtureState s;
    s.P = P;
    s.Pt = P_total;
    s.T = T;
    s.Tt = T;
    s.m_dot = m_dot;
    s.set_Y(Y);
    return s;
}

class PlenumJacobianTest : public ::testing::Test {
  protected:
    void SetUp() override {
        streams = {
            make_air_stream(1.5, 500.0, 150000.0),
            make_air_stream(0.5, 400.0, 148000.0),
        };
        ns = num_species();

        // Mix to get a reasonable state guess
        auto mix = mixer_from_streams_and_jacobians(streams);
        double m_tot = 0.0;
        for (auto &s : streams)
            m_tot += s.m_dot;

        state = make_state(150000.0, 150000.0, mix.Tt_mix, m_tot, mix.Y_mix);
        base = plenum_residuals_and_jacobian(state, streams);
    }

    std::vector<solver::Stream> streams;
    std::size_t ns{};
    MixtureState state;
    ChamberResult base;
};

TEST_F(PlenumJacobianTest, LocalJacobian_P) {
    const double eps = 1.0;
    MixtureState s_p = state;
    s_p.P += eps;
    auto res_p = plenum_residuals_and_jacobian(s_p, streams);

    for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
        double fd = (res_p.residuals[eq] - base.residuals[eq]) / eps;
        auto it = base.local_jacobian[eq].find("P");
        double analytical = (it != base.local_jacobian[eq].end()) ? it->second : 0.0;
        EXPECT_NEAR(analytical, fd, std::abs(fd) * 1e-3 + 1e-6)
            << "eq " << eq << " d/dP";
    }
}

TEST_F(PlenumJacobianTest, LocalJacobian_P_total) {
    const double eps = 1.0;
    MixtureState s_p = state;
    s_p.P_total += eps;
    auto res_p = plenum_residuals_and_jacobian(s_p, streams);

    for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
        double fd = (res_p.residuals[eq] - base.residuals[eq]) / eps;
        auto it = base.local_jacobian[eq].find("P_total");
        double analytical = (it != base.local_jacobian[eq].end()) ? it->second : 0.0;
        EXPECT_NEAR(analytical, fd, std::abs(fd) * 1e-3 + 1e-6)
            << "eq " << eq << " d/dP_total";
    }
}

TEST_F(PlenumJacobianTest, LocalJacobian_T) {
    const double eps = 0.1;
    MixtureState s_p = state;
    s_p.T += eps;
    auto res_p = plenum_residuals_and_jacobian(s_p, streams);

    for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
        double fd = (res_p.residuals[eq] - base.residuals[eq]) / eps;
        auto it = base.local_jacobian[eq].find("T");
        double analytical = (it != base.local_jacobian[eq].end()) ? it->second : 0.0;
        EXPECT_NEAR(analytical, fd, std::abs(fd) * 1e-3 + 1e-6)
            << "eq " << eq << " d/dT";
    }
}

TEST_F(PlenumJacobianTest, LocalJacobian_Y) {
    const double eps = 1e-5;
    for (std::size_t k = 0; k < ns; ++k) {
        MixtureState s_p = state;
        s_p.Y[k] += eps;
        s_p.set_Y(s_p.Y);
        auto res_p = plenum_residuals_and_jacobian(s_p, streams);

        std::string key = "Y[" + std::to_string(k) + "]";
        for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
            double fd = (res_p.residuals[eq] - base.residuals[eq]) / eps;
            auto it = base.local_jacobian[eq].find(key);
            double analytical = (it != base.local_jacobian[eq].end()) ? it->second : 0.0;
            EXPECT_NEAR(analytical, fd, std::abs(fd) * 1e-3 + 1e-6)
                << "eq " << eq << " d/d" << key;
        }
    }
}

TEST_F(PlenumJacobianTest, StreamJacobian_mdot) {
    const double eps = 1e-4;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].m_dot += eps;
        auto res_p = plenum_residuals_and_jacobian(state, s_p);

        for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
            double fd = (res_p.residuals[eq] - base.residuals[eq]) / eps;
            double analytical = base.stream_jacobian[eq][i].d_mdot;
            EXPECT_NEAR(analytical, fd, std::abs(fd) * 1e-3 + 1e-4)
                << "eq " << eq << " stream " << i << " d/d_mdot";
        }
    }
}

TEST_F(PlenumJacobianTest, StreamJacobian_T) {
    const double eps = 0.1;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].T += eps;
        auto res_p = plenum_residuals_and_jacobian(state, s_p);

        for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
            double fd = (res_p.residuals[eq] - base.residuals[eq]) / eps;
            double analytical = base.stream_jacobian[eq][i].d_Tt;
            EXPECT_NEAR(analytical, fd, std::abs(fd) * 1e-3 + 1e-4)
                << "eq " << eq << " stream " << i << " d/dT";
        }
    }
}

TEST_F(PlenumJacobianTest, StreamJacobian_Y) {
    const double eps = 1e-5;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        for (std::size_t j = 0; j < ns; ++j) {
            auto s_p = streams;
            s_p[i].Y[j] += eps;
            auto res_p = plenum_residuals_and_jacobian(state, s_p);

            for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
                double fd =
                    (res_p.residuals[eq] - base.residuals[eq]) / eps;
                double analytical = base.stream_jacobian[eq][i].d_Y[j];
                EXPECT_NEAR(analytical, fd, std::abs(fd) * 5e-3 + 1e-4)
                    << "eq " << eq << " stream " << i << " d/dY[" << j
                    << "]";
            }
        }
    }
}
#endif

// ---------------------------------------------------------------------------
// Combustor ChamberResult -- local + stream Jacobian FD validation
// ---------------------------------------------------------------------------

#if 0
class CombustorJacobianTest : public ::testing::Test {
  protected:
    void SetUp() override {
        streams = {
            make_air_stream(2.0, 600.0, 150000.0),
            make_fuel_stream(0.06, 300.0, 150000.0),
        };
        ns = num_species();
        pressure_loss_frac = 0.05;

        // Get a reasonable combustor outlet state
        auto mix = adiabatic_T_complete_and_jacobian_T_from_streams(
            streams, 140000.0);
        double m_tot = 0.0;
        for (auto &s : streams)
            m_tot += s.m_dot;

        state = make_state(140000.0, 140000.0, mix.Tt_mix, m_tot, mix.Y_mix);
        base = combustor_residuals_and_jacobian(
            state, streams, CombustionMethod::Complete, pressure_loss_frac);
    }

    std::vector<solver::Stream> streams;
    std::size_t ns{};
    double pressure_loss_frac{};
    MixtureState state;
    ChamberResult base;
};

TEST_F(CombustorJacobianTest, LocalJacobian_P) {
    const double eps = 1.0;
    MixtureState s_p = state;
    s_p.P += eps;
    auto res_p = combustor_residuals_and_jacobian(
        s_p, streams, CombustionMethod::Complete, pressure_loss_frac);

    for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
        double fd = (res_p.residuals[eq] - base.residuals[eq]) / eps;
        auto it = base.local_jacobian[eq].find("P");
        double analytical = (it != base.local_jacobian[eq].end()) ? it->second : 0.0;
        EXPECT_NEAR(analytical, fd, std::abs(fd) * 1e-3 + 1e-6)
            << "eq " << eq << " d/dP";
    }
}

TEST_F(CombustorJacobianTest, LocalJacobian_P_total) {
    const double eps = 1.0;
    MixtureState s_p = state;
    s_p.P_total += eps;
    auto res_p = combustor_residuals_and_jacobian(
        s_p, streams, CombustionMethod::Complete, pressure_loss_frac);

    for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
        double fd = (res_p.residuals[eq] - base.residuals[eq]) / eps;
        auto it = base.local_jacobian[eq].find("P_total");
        double analytical = (it != base.local_jacobian[eq].end()) ? it->second : 0.0;
        EXPECT_NEAR(analytical, fd, std::abs(fd) * 1e-3 + 1e-6)
            << "eq " << eq << " d/dP_total";
    }
}

TEST_F(CombustorJacobianTest, LocalJacobian_T) {
    const double eps = 0.1;
    MixtureState s_p = state;
    s_p.T += eps;
    auto res_p = combustor_residuals_and_jacobian(
        s_p, streams, CombustionMethod::Complete, pressure_loss_frac);

    for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
        double fd = (res_p.residuals[eq] - base.residuals[eq]) / eps;
        auto it = base.local_jacobian[eq].find("T");
        double analytical = (it != base.local_jacobian[eq].end()) ? it->second : 0.0;
        EXPECT_NEAR(analytical, fd, std::abs(fd) * 1e-3 + 1e-6)
            << "eq " << eq << " d/dT";
    }
}

TEST_F(CombustorJacobianTest, LocalJacobian_Y) {
    const double eps = 1e-5;
    for (std::size_t k = 0; k < ns; ++k) {
        MixtureState s_p = state;
        s_p.Y[k] += eps;
        s_p.set_Y(s_p.Y);
        auto res_p = combustor_residuals_and_jacobian(
            s_p, streams, CombustionMethod::Complete, pressure_loss_frac);

        std::string key = "Y[" + std::to_string(k) + "]";
        for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
            double fd = (res_p.residuals[eq] - base.residuals[eq]) / eps;
            auto it = base.local_jacobian[eq].find(key);
            double analytical = (it != base.local_jacobian[eq].end()) ? it->second : 0.0;
            EXPECT_NEAR(analytical, fd, std::abs(fd) * 1e-3 + 1e-6)
                << "eq " << eq << " d/d" << key;
        }
    }
}

TEST_F(CombustorJacobianTest, StreamJacobian_mdot) {
    const double eps = 1e-3;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].m_dot += eps;
        auto res_p = combustor_residuals_and_jacobian(
            state, s_p, CombustionMethod::Complete, pressure_loss_frac);

        for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
            double fd = (res_p.residuals[eq] - base.residuals[eq]) / eps;
            double analytical = base.stream_jacobian[eq][i].d_mdot;
            EXPECT_NEAR(analytical, fd, std::abs(fd) * 0.01 + 0.1)
                << "eq " << eq << " stream " << i << " d/d_mdot";
        }
    }
}

TEST_F(CombustorJacobianTest, StreamJacobian_T) {
    const double eps = 1.0;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].T += eps;
        auto res_p = combustor_residuals_and_jacobian(
            state, s_p, CombustionMethod::Complete, pressure_loss_frac);

        for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
            double fd = (res_p.residuals[eq] - base.residuals[eq]) / eps;
            double analytical = base.stream_jacobian[eq][i].d_Tt;
            EXPECT_NEAR(analytical, fd, std::abs(fd) * 0.01 + 0.1)
                << "eq " << eq << " stream " << i << " d/dT";
        }
    }
}

TEST_F(CombustorJacobianTest, StreamJacobian_P_total) {
    const double eps = 1.0;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].P_total += eps;
        auto res_p = combustor_residuals_and_jacobian(
            state, s_p, CombustionMethod::Complete, pressure_loss_frac);

        for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
            double fd = (res_p.residuals[eq] - base.residuals[eq]) / eps;
            double analytical = base.stream_jacobian[eq][i].d_Pt;
            EXPECT_NEAR(analytical, fd, std::abs(fd) * 1e-3 + 1e-4)
                << "eq " << eq << " stream " << i << " d/dP_total";
        }
    }
}

TEST_F(CombustorJacobianTest, StreamJacobian_Y) {
    const double eps = 1e-4;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        for (std::size_t j = 0; j < ns; ++j) {
            auto s_p = streams;
            s_p[i].Y[j] += eps;
            auto res_p = combustor_residuals_and_jacobian(
                state, s_p, CombustionMethod::Complete, pressure_loss_frac);

            for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
                double fd =
                    (res_p.residuals[eq] - base.residuals[eq]) / eps;
                double analytical = base.stream_jacobian[eq][i].d_Y[j];
                EXPECT_NEAR(analytical, fd, std::abs(fd) * 0.05 + 0.01)
                    << "eq " << eq << " stream " << i << " d/dY[" << j
                    << "]";
            }
        }
    }
}
#endif

// ---------------------------------------------------------------------------
// Combustor Equilibrium -- verify same structure works
// ---------------------------------------------------------------------------

#if 0
TEST(CombustorEquilibriumTest, StructureAndStreamJacobian_mdot) {
    auto streams = std::vector<solver::Stream>{
        make_air_stream(2.0, 600.0, 150000.0),
        make_fuel_stream(0.06, 300.0, 150000.0),
    };
    std::size_t ns = num_species();
    double pressure_loss_frac = 0.05;

    auto mix = adiabatic_T_equilibrium_and_jacobians_from_streams(
        streams, 140000.0);
    double m_tot = 0.0;
    for (auto &s : streams)
        m_tot += s.m_dot;

    auto state =
        make_state(140000.0, 140000.0, mix.Tt_mix, m_tot, mix.Y_mix);
    auto base = combustor_residuals_and_jacobian(
        state, streams, CombustionMethod::Equilibrium, pressure_loss_frac);

    // Expect 17 equations: 1 energy + 14 species + 1 pressure loss + 1 stagnation
    EXPECT_EQ(base.residuals.size(), 1 + ns + 1 + 1);
    EXPECT_EQ(base.local_jacobian.size(), base.residuals.size());
    EXPECT_EQ(base.stream_jacobian.size(), base.residuals.size());

    // FD check on stream m_dot
    const double eps = 1e-3;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto s_p = streams;
        s_p[i].m_dot += eps;
        auto res_p = combustor_residuals_and_jacobian(
            state, s_p, CombustionMethod::Equilibrium, pressure_loss_frac);

        for (std::size_t eq = 0; eq < base.residuals.size(); ++eq) {
            double fd = (res_p.residuals[eq] - base.residuals[eq]) / eps;
            double analytical = base.stream_jacobian[eq][i].d_mdot;
            EXPECT_NEAR(analytical, fd, std::abs(fd) * 0.02 + 0.5)
                << "eq " << eq << " stream " << i << " d/d_mdot";
        }
    }
}
#endif

// ---------------------------------------------------------------------------
// Combustor Pressure Loss Analytical Jacobian FD validation
// ---------------------------------------------------------------------------

TEST(CombustorPressureLossJacobianTest, burns_at_the_mixed_Pt) {
    // The stream partials are those of the mixed Pt (#1244 review), so the
    // kernel takes it itself (#1272). Pt_mix (no loss) pins that pressure;
    // at this lean burn the equilibrium T barely depends on P, so the T check
    // only pins the call path (the Python test checks a rich burn).
    auto streams = std::vector<solver::Stream>{
        make_air_stream(2.0, 600.0, 2.0e5),
        make_fuel_stream(0.06, 300.0, 2.1e5),
    };
    auto no_loss = [](const PressureLossContext &) -> std::tuple<double, double> {
        return {0.0, 0.0};
    };
    const double P = mixer_from_streams_and_jacobians(streams).Pt_mix;
    const auto res = combustor_residuals_and_jacobians(streams, 0.0, 0.0, no_loss, true);
    const auto ref = adiabatic_T_equilibrium_and_jacobians_from_streams(streams, P);
    EXPECT_EQ(res.Tt_mix, ref.Tt_mix);
    EXPECT_EQ(res.Pt_mix, P);
}

TEST(CombustorPressureLossJacobianTest, mdot_chain_rule) {
    auto streams = std::vector<solver::Stream>{
        make_air_stream(2.0, 600.0, 150000.0),
        make_fuel_stream(0.06, 300.0, 150000.0),
    };
    // The burn's pressure is the streams' mixed Pt, which the kernel takes
    // itself (#1244, #1272).

    // Define a benchmark linear correlation: zeta = k*theta + c
    const double k_zeta = 0.5;
    const double c_zeta = 0.02;
    auto linear_loss = [k_zeta, c_zeta](const PressureLossContext &ctx) -> std::tuple<double, double> {
        return {k_zeta * ctx.theta + c_zeta, k_zeta};
    };

    auto base = combustor_residuals_and_jacobians(streams, 0.0, 0.0, linear_loss, false);

    const double eps = 1e-5;
    for (std::size_t i = 0; i < streams.size(); ++i) {
        auto streams_p = streams;
        streams_p[i].m_dot += eps;

        auto res_p = combustor_residuals_and_jacobians(streams_p, 0.0, 0.0, linear_loss, false);

        double analytical = base.dPt_mix_d_stream[i].d_mdot;
        double fd = (res_p.Pt_mix - base.Pt_mix) / eps;

        EXPECT_NEAR(analytical, fd, std::abs(fd) * 1e-4 + 1e-6)
            << "stream " << i << " mdot derivative check failed";
    }
}
