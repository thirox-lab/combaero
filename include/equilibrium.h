#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include "state.h"
#include <vector>
#include "combustion.h"
#include "thermo.h"

namespace combaero {

// -------------------------------------------------------------
// State-based WGS equilibrium functions
// -------------------------------------------------------------
EquilibriumResult wgs_equilibrium(const State &in);
EquilibriumResult wgs_equilibrium_adiabatic(const State &in);

// -------------------------------------------------------------
// Reforming + WGS equilibrium
// -------------------------------------------------------------
// Higher-fidelity product model for rich hydrocarbon combustion.
// Produces CO, H2, CO2, H2O, N2 (ignoring soot/solid carbon).
//
// The reactions (steam reforming of each hydrocarbon, and WGS) are solved
// SIMULTANEOUSLY (#759): one extent per reaction, Newton on
// F_r = ln Q_r - ln K_r = 0 with the analytic Jacobian; each reaction's own
// species (its hydrocarbon; CO2 for WGS) steps in ln n, the others are kept
// positive by a fraction-to-boundary step. A single round of scalar solves
// leaves reforming short of equilibrium once WGS has moved CO, H2 and H2O --
// at 650-850 K, the pre-reformer range, H2 came out 7-49% low against Cantera
// -- and repeated rounds (#755) crawl once the steam is used up (rich propane
// at phi >= 2.5 stopped 2000 rounds short, 120-320 K hot). Rounds now only
// start the Newton solve: up to REFORMING_START_PASSES of them, stopping early
// once no amount moves by more than REFORMING_PASS_RTOL of the total.
constexpr int REFORMING_START_PASSES = 3;
constexpr double REFORMING_PASS_RTOL = 1e-14;
constexpr int REFORMING_NEWTON_MAX_IT = 200;
// Newton iterations from the rounds' start before it restarts from the
// element potentials (Gordon and McBride, NASA RP-1311): far from equilibrium
// -- the steam of a steam-starved burn must grow 1e-6 -> 5e-2 -- the
// linearisation in the extents is poor and the steps short.
constexpr int REFORMING_NEWTON_ROUNDS_IT = 50;
// Converged when max_r |ln Q_r - ln K_r| is below this (then polished to
// round-off with up to two more steps).
constexpr double REFORMING_NEWTON_FTOL = 1e-10;
// A step keeps at least 1 - this of every amount (other than a reaction's own
// species, which steps in ln n) it reduces.
constexpr double REFORMING_FRACTION_TO_BOUNDARY = 0.99;
// A Newton step may raise a trace species in ln n up to this fraction of the
// total, linearly beyond (RP-1311's trace-species bound, 1e-4).
constexpr double REFORMING_TRACE_FRACTION = 1e-4;
// A reaction with an absent species is left out of the Newton solve; a trial
// round checks it cannot move by more than this of the total, else it is
// taken and the solve repeated, up to REFORMING_JOINT_RESTARTS times.
constexpr double REFORMING_JOINT_CHECK_RTOL = 1e-10;
constexpr int REFORMING_JOINT_RESTARTS = 10;
// The adiabatic temperature (WGS and reforming) is bracketed in
// [EQUILIBRIUM_T_MIN, EQUILIBRIUM_T_MAX] K from its guess (the low end the
// thermo floor: the former 300 K returned 300 K for a colder unreacting inlet) -- the first step the
// frozen-composition Newton step, then EQUILIBRIUM_T_BRACKET_STEP if that is
// not finite, doubling -- and solved by Brent's method to a bracket of
// EQUILIBRIUM_T_XTOL_REL T (#927). No root in the range throws.
constexpr double EQUILIBRIUM_T_MIN = thermo::T_FLOOR;
constexpr double EQUILIBRIUM_T_MAX = 4500.0;
constexpr double EQUILIBRIUM_T_BRACKET_STEP = 100.0;
constexpr double EQUILIBRIUM_T_XTOL_REL = 1e-14;
// The WGS extent is converged when |ln Q - ln Kp| is below EQUILIBRIUM_WGS_FTOL
// or below EQUILIBRIUM_WGS_NOISE_FACTOR times its rounding error (a trace
// amount at the root carries only eps |n0| / n of relative precision).
constexpr double EQUILIBRIUM_WGS_FTOL = 1e-12;
constexpr double EQUILIBRIUM_WGS_NOISE_FACTOR = 8.0;
// Reference pressure of the equilibrium constants [Pa].
constexpr double EQUILIBRIUM_P_REF = 101325.0;
EquilibriumResult reforming_equilibrium(const State &in);
EquilibriumResult reforming_equilibrium_adiabatic(const State &in);

// -------------------------------------------------------------
// High-level combustion interface
// -------------------------------------------------------------

// Complete combustion (with the Phi = 1 softmin, SMOOTHING_K_PHI1, always on:
// #1048), then reforming + WGS equilibrium on the products, adiabatic.
// smooth_phi0 selects the zero-fuel smoothing of the first step.
EquilibriumResult combustion_equilibrium(const State &reactants,
                                         bool smooth_phi0 = true,
                                         double k0 = SMOOTHING_K_PHI0);

// combustion_equilibrium with analytic derivatives (#716): one implicit-function
// solve on the joint extents (reforming per hydrocarbon, WGS; #759) and the
// adiabatic root, chained through complete_combustion_products_and_jacobian.
// X is taken unnormalised (chain with a dX that keeps sum X = 1). A reaction
// with an absent species cannot move and keeps a zero derivative.
struct EquilibriumJacobian {
    State state;                                // as combustion_equilibrium
    double dT_dT_in = 0.0;
    double dT_dP = 0.0;
    std::vector<double> dT_dX;                  // [reactant k]
    std::vector<std::vector<double>> dXout_dX;  // [product j][reactant k]
    std::vector<double> dXout_dT_in;
    std::vector<double> dXout_dP;
};
EquilibriumJacobian combustion_equilibrium_and_jacobian(const State &reactants,
                                                        bool smooth_phi0 = true,
                                                        double k0 = SMOOTHING_K_PHI0);

} // namespace combaero

#endif // EQUILIBRIUM_H
