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
// Partial by design: no NO/NO2 (so no NOx) and nothing dissociates. N2 is
// inert except where NH3 is present: there the decomposition
// 2 NH3 = N2 + 3 H2 is equilibrated with the others (#1314). The reactions
// cover the current species table only; a species added to it is not
// equilibrated (see docs/API_CPP.md, Chemical Equilibrium, and
// python/tests/test_equilibrium_scope.py).
//
// The reactions (steam reforming of each hydrocarbon, WGS, and NH3
// decomposition where NH3 is present) are solved
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
// Where a hydrocarbon is not a trace, Newton runs in another basis of the
// same reactions (#1262): that hydrocarbon's reactions with two of the
// oxygen carriers, keyed on the two scarcer carriers, and every other
// hydrocarbon's exchange with it. With the oxygen nearly gone the
// steam-reforming basis left the hydrocarbon interconversions as differences
// of rows weighted by 1/n_H2O (condition number 2.5e11) and the trace
// carriers as differences of larger changes; 201 of 375 trace-oxygen and
// very rich methane states did not converge, none does now. Likewise WGS is
// keyed on CO where it is the scarcer carbon oxide (a lean burn with trace
// fuel), and a reaction driven toward an absent species (pyrolysis of a fuel
// with no oxygen too, #1251) is run by a trace,
// seeding it, and the solve repeated (REFORMING_SEED_FRACTION). Where a
// species is a non-key of that basis (n plus the extents' changes) and
// scarcer than any species of a component basis -- the most abundant
// independent species as components, each other species keying its own
// formation from them -- the component basis is taken (#1288): H2 or CO fuel
// with a hydrocarbon, cold, left CO2, steam or H2 at 1e-20 to 1e-63 a
// non-key, and it came back negative.
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
// Where no step can lower max |ln Q - ln K| any further (its round-off
// floor), the solve is converged if it is below this (#1262). A non-key trace
// species is the sum of the keys' changes, each rounded to its key's ulp: a
// hydrocarbon key at 1e-6 moves hydrogen at 5e-9 in steps of 1e-22, and
// heptane's 22 H2 in its ln Q left a floor of 1.1e-10 (lean propane, 250 K,
// 20 bar). 1e-9 in ln Q is 1e-9 relative in the trace amount.
constexpr double REFORMING_NEWTON_FLOOR_FTOL = 1e-9;
// Below about 1e-308 an amount is subnormal and holds fewer digits: one ulp
// of steam at 4e-315 is 1.2e-9 of it, 8.5e-9 in the ln Q of a row with 7 of
// it (heptane with 1e-20 of air, 300 K; #1287). Where no step lowers ||F||,
// a row is also accepted within this many such ulp steps of its amounts,
// sum_j |nu_j| ulp(n_j) / n_j, which is about 1e-14 for normal amounts.
constexpr double REFORMING_NEWTON_ULP_FLOOR_STEPS = 2.0;
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
// The joint solve takes the hydrocarbon-rich basis (the most abundant
// hydrocarbon's reactions with the oxygen carriers, every other
// hydrocarbon's exchange with it; the same equilibrium) when that
// hydrocarbon is at least this fraction of the total (#1262). Where H2 is
// scarcer than the second most abundant hydrocarbon (pyrolysis), that one
// joins as a second reference and H2 keys its own exchange with the two.
constexpr double REFORMING_RICH_BASIS_HC_FRACTION = 1e-4;
// The element-potential start is handed to the Newton polish once its full
// corrections in ln n are below this, converged or not (#1262).
constexpr double REFORMING_EP_START_CORRECTION = 1e-6;
// Before its write-back as extents, the element-potential start's amounts are
// balanced to this of each element's total, by up to this many Newton steps
// on the element equations alone (#1285).
constexpr double REFORMING_EP_BALANCE_RTOL = 1e-15;
constexpr int REFORMING_EP_BALANCE_IT = 10;
// A log-stepping key falling by more than this many times its amount in one
// step is taken from exp(d / n_k), not n_k plus its extent, which cancels
// there (#1285).
constexpr double REFORMING_LOG_KEY_EXP_BELOW = -1.0;
// The element-potential start's component species are chosen by amount,
// each independent of those before it: its formula's remainder after
// elimination against theirs exceeds this (formulas are small integers).
constexpr double REFORMING_EP_COMPONENT_RANK_TOL = 1e-9;
// A reaction with an absent species whose other side is present is run
// toward it by this fraction of the total, seeding it, and the solve
// repeated, up to REFORMING_SEED_PASSES times; one still runnable after them
// means not converged (#1262).
constexpr double REFORMING_SEED_FRACTION = 1e-10;
constexpr int REFORMING_SEED_PASSES = 3;
// A solved state's element totals (C, H, O, N) match the input's to this, each
// relative to its own total, or the solve is not converged (#1262). Relative
// to the largest total, it allowed about 1e-5 of hydrogen's own at phi 1e-7
// (#1285 review).
constexpr double REFORMING_ELEMENT_RTOL = 1e-12;
// The Jacobian of a solution with an absent species that the complete-
// combustion step feeds is taken with that species at a trace of the total,
// equilibrated: the one-sided limit (#1251). From this, up by 100 to
// REFORMING_JACOBIAN_TRACE_MAX, the smallest at which the solve converges;
// none means refused.
constexpr double REFORMING_JACOBIAN_TRACE = 1e-20;
constexpr double REFORMING_JACOBIAN_TRACE_MAX = 1e-12;
// From the first converging trace the ladder steps down by 100 until the
// next trace's derivatives agree with it to this, relative to each
// quantity's size (#1241): the trace is then in the one-sided limit's regime
// (its derivatives change linearly in the trace, so a 1e-8 gap over a factor
// 100 leaves the kept ones within about 1e-8 of the limit). Not fitted: a
// limit tolerance well below the 1e-6 the derivative tests ask of J.
constexpr double REFORMING_JACOBIAN_LIMIT_RTOL = 1e-8;
constexpr int REFORMING_JOINT_RESTARTS = 10;
// Which reactions the equilibrium runs (#1241; no threshold): reforming of
// each hydrocarbon, with its backward direction, methanation, where the
// burned state holds a reducing species (H2, CO, a hydrocarbon or NH3) and
// both carbon and hydrogen; NH3 synthesis and decomposition, 2 NH3 = N2 +
// 3 H2, where it holds one and both nitrogen and hydrogen (N then joins the
// element balances, and N2 is not inert). A burned state with no reducing
// species, a complete lean burn or a fuel-free mix, takes the shift alone,
// which cannot move there: bit for bit as before. Up to 0.12.0 both were cut
// at 1e-10 of the burned state (EQUILIBRIUM_HC_MIN_FRACTION and
// EQUILIBRIUM_NH3_MIN_FRACTION, removed): below the cut no CH4 or NH3 formed
// from CO, H2 and N2 (H2/CO at phi 14 and 32.7 bar came out 320 K below
// Cantera, NH3 with CO up to 497 K), and T stepped where a trace crossed it
// (NH3 synthesis switching on: up to 1.64 K at 100 bar).
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
// Secant steps that polish T_ad on converged inner solves where an inner
// solve at a bracket end failed (#1262).
constexpr int EQUILIBRIUM_T_POLISH_IT = 20;
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
// solve on the joint extents (reforming per hydrocarbon, WGS, NH3
// decomposition where NH3 is present; #759, #1314) and the
// adiabatic root, chained through complete_combustion_products_and_jacobian.
// X is taken unnormalised (chain with a dX that keeps sum X = 1). A reaction
// with an absent species cannot move and keeps a zero derivative. Throws
// std::runtime_error where the equilibrium does not converge (#1262): the
// implicit-function derivatives are those of a root only.
//
// Where the burned state lacks a species that an arriving reactant feeds (a
// fuel-free mix as fuel arrives, a fuel-only one as air does), T and X_out
// have one-sided slopes along each direction dX, but in dry air no matrix
// holds them: steam is a trace too, and the slopes along two fuels do not add
// (CH4 and CO at 335 K: 11 %; #1241). dX_dir, if given, is the direction the
// caller will apply the derivatives along (a zero-flow stream's mixing
// direction): the trace states that take the one-sided limit are seeded along
// it, so J dX_dir is that direction's slope. Without it they are seeded as a
// fixed trace of each fed species, and J holds a slope along each reactant
// alone only where those slopes add. Elsewhere dX_dir changes nothing.
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
                                                        double k0 = SMOOTHING_K_PHI0,
                                                        const std::vector<double> *dX_dir = nullptr);

} // namespace combaero

#endif // EQUILIBRIUM_H
