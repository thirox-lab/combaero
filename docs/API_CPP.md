# CombAero C++ API Reference

This document provides the technical reference for the CombAero C++ library. For Python usage, see [API_PYTHON.md](API_PYTHON.md).

## Table of Contents
- [Input validation (input_validation.h)](#input-validation-input_validationh)
- [Thermodynamics (thermo.h)](#thermodynamics-thermoh)
- [Transport Properties (transport.h)](#transport-properties-transporth)
- [Combustion (combustion.h)](#combustion-combustionh)
- [Chemical Equilibrium (equilibrium.h)](#chemical-equilibrium-equilibriumh)
- [Stagnation / Static Conversions (stagnation.h)](#stagnation--static-conversions-stagnationh)
- [Compressible Flow (compressible.h)](#compressible-flow-compressibleh)
- [Incompressible Flow (incompressible.h)](#incompressible-flow-incompressibleh)
- [Area Change (area_change.h)](#area-change-area_changeh)
- [Friction Factor Correlations (friction.h)](#friction-factor-correlations-frictionh)
- [Heat Transfer Correlations (heat_transfer.h)](#heat-transfer-correlations-heat_transferh)
  - [Nusselt Numbers](#nusselt-numbers)
  - [Overall Heat Transfer](#overall-heat-transfer)
  - [Channel Flow Models](#channel-flow-models)
  - [Wall Coupling](#wall-coupling)
  - [Liner Radiation](#liner-radiation-liner_radiationh-793)
  - [Data Structures](#data-structures)
- [Geometry Utilities (geometry.h)](#geometry-utilities-geometryh)
- [Orifice Flow (orifice.h)](#orifice-flow-orificeh)
- [Acoustics (acoustics.h)](#acoustics-acousticsh)
- [Humid Air (humidair.h)](#humid-air-humidairh)
- [Network Solver Interface (solver_interface.h)](#network-solver-interface-solver_interfaceh)
  - [Tee Junction Components](#tee-junction-components-tee_junctionh--solver_interfaceh)

---

## Input validation (input_validation.h)

A domain check written `if (x <= 0) throw` lets NaN through, because every
comparison with NaN is false. Public functions check input with the NaN-safe
predicates below, which NaN and +/-inf fail, and keep the exception type they
already threw (#769):

```cpp
namespace combaero::validation {
bool is_finite(double x);
bool is_finite_positive(double x);                    // finite, x > 0
bool is_finite_nonnegative(double x);                 // finite, x >= 0
bool is_finite_in_range(double x, double lo, double hi);  // finite, lo <= x <= hi

// Throw E(message) when the predicate fails; E defaults to std::invalid_argument.
template <class E = std::invalid_argument> void require_finite(double x, const char* message);
template <class E = std::invalid_argument> void require_finite_positive(double x, const char* message);
template <class E = std::invalid_argument> void require_finite_nonnegative(double x, const char* message);
template <class E = std::invalid_argument>
void require_finite_in_range(double x, double lo, double hi, const char* message);
}
```

Where +inf is a physical limit rather than bad input (NTU, a film coefficient
h, YSZ aging hours) the check refuses NaN only.

A clamp, a root finder or a "no correction" branch can also turn NaN into a
plausible number without any comparison in sight (`std::max(NaN, x)` is `x`
for one argument order). Such arguments refuse non-finite values too, while
finite values keep their existing meaning: `velocity_of_approach_factor` is
still 1 for a finite beta outside (0, 1), `D_h <= 0` still selects the
equivalent-circle diameter in the area-change closures, `m_ref <= 0` still
switches the combustor low-flow ramp off. The same holds for the "neutral
answer" branches (no fin model in `pin_fin_array_efficiency`, Gc/Gj = 0 in
`crossflow_to_jet_ratio_at_x`, Y = 1 for gamma <= 1 in the McGreehan-Schotsch
expansion factors, zero flow in the channel, area-change, jet-impulse,
crossflow-Cd and orifice kernels): finite input keeps them, NaN and inf are
refused. Solver-iterate arguments are refused only when non-finite; their
negative or zero values keep the kernels' no-flow and floor branches.

---

## Thermodynamics (thermo.h)

> [!NOTE]
> **Base Units**: Temperature [K], Pressure [Pa], Mass [kg], Energy [J], Power [W].
> Compositions are `std::vector<double>` of mole fractions (sum to 1.0) unless explicitly marked as `_mass` (mass fractions).

### Composition Utilities (composition.h)
```cpp
// Mass/Mole fraction conversions
double mwmix(const std::vector<double>& X);
std::vector<double> mole_to_mass(const std::vector<double>& X);
std::vector<double> mass_to_mole(const std::vector<double>& Y);
std::vector<double> normalize_fractions(const std::vector<double>& X);  // negatives -> 0, then sum to 1; absent species stay 0
std::vector<double> convert_to_dry_fractions(const std::vector<double>& X);
```

### Species and Mixture Metadata (thermo.h)
```cpp
// Get species metadata
std::string species_name(std::size_t index);
std::size_t species_index_from_name(const std::string& name);
double species_molar_mass(std::size_t index);
double species_molar_mass_from_name(const std::string& name);
std::size_t num_species();
```

### Thermodynamic Properties (Molar Basis)
```cpp
double cp(double T, const std::vector<double>& X);
double cv(double T, const std::vector<double>& X);
double h(double T, const std::vector<double>& X);
double u(double T, const std::vector<double>& X);
double s(double T, const std::vector<double>& X, double P, double P_ref = 101325.0);
```

### Thermodynamic Properties (Mass Basis)
```cpp
double cp_mass(double T, const std::vector<double>& X);
// Secant cp over [T1, T2] and its partials in T1, T2 (cp(T1) as T2 -> T1):
// the variable-property capacity rate of an epsilon-NTU wall exchange (#627).
std::tuple<double, double, double> cp_secant_mass_and_jacobian(
    double T1, double T2, const std::vector<double>& X);
// Its gradient in the MASS fractions, per species: a wall side's capacity
// rate moves with its composition (#717).
std::vector<double> cp_secant_mass_dY(double T1, double T2, const std::vector<double>& X);
double cv_mass(double T, const std::vector<double>& X);
double h_mass(double T, const std::vector<double>& X);
double u_mass(double T, const std::vector<double>& X);
double s_mass(double T, double P, const std::vector<double>& X, double P_ref = 101325.0);
double density(double T, double P, const std::vector<double>& X);
double speed_of_sound(double T, const std::vector<double>& X);
double isentropic_expansion_coefficient(double T, const std::vector<double>& X);
```

### Inverse Solvers
```cpp
double calc_T_from_h(double h_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_s(double s_target, double P, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_cp(double cp_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
struct CpRoot { double T; double dT_dcp; bool double_root; };
std::vector<CpRoot> calc_T_from_cp_roots(double cp_target, const std::vector<double>& X, double T_min = thermo::CP_INVERSE_T_MIN, double T_max = thermo::CP_INVERSE_T_MAX, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_u(double u_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_h_mass(double h_mass_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_s_mass(double s_mass_target, double P, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_u_mass(double u_mass_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
```

h, s, u and the sv form are strictly increasing in T at every T (the thermo
tails below keep cp > 0), so each inverse has one root: safeguarded Newton
(factor-4 steps until bracketed, bisection when a step leaves the bracket,
converged only on a Newton step). They throw `std::runtime_error` when they do
not converge -- a target no positive T reaches -- rather than return the last
iterate (#606). cp is not monotone, so `calc_T_from_cp` returns the root
nearest `T_guess` in [`thermo::CP_INVERSE_T_MIN`, `thermo::CP_INVERSE_T_MAX`]
= [1, 6000] K: an outward scan from the guess (clamped into the range), each
step on the side whose next node is nearer the guess, then bracketed Newton on
the first sign change (#868 R6). Nodes are `thermo::CP_INVERSE_SCAN_STEP` =
10 K apart down to `thermo::CP_INVERSE_LOG_SCAN_T` = 100 K and 10 % of T below
it, and the downward scan steps onto `thermo::CP_INVERSE_T_DATA` = 200 K, so
a guess in [200, 6000] K and a target with no root below 200 K give the result
of the former [200, 6000] K range bit for bit. Above 6000 K the high tail
levels cp off at a constant that is not physics. Below 200 K the low tail is
the ideal gas to 0 K, but cp - cp(0 K) ~ T^3, so the root's condition number
grows as 1/T^3: at 1 K, CO's T is fixed to ~5e-7 relative, at 0.1 K only to
~5e-4 (see the constants' comment in `thermo.h`). No root in the range, or no
convergence in `max_iter`, throws `std::runtime_error`; a non-finite or
non-positive `cp_target`, `T_guess` or `tol` throws `std::invalid_argument`.

`calc_T_from_cp_roots` is the opt-in alternative: every root in
[`T_min`, `T_max`] (within [1, 6000] K), sorted, empty when there is none
(#963). It needs no guess, so it also finds the pairs closer than a scan step
that `calc_T_from_cp` misses (O2's cp(127 K), again at 120.8 K). dcp/dT is
evaluated on the same scan nodes (from `T_min` up); each sign change is
bisected to the last bit, a cp extremum, and where |dcp/dT| has a local
minimum on the nodes without a sign change, a golden-section search
(`thermo::CP_ROOTS_DIP_ITER` steps) looks for an extremum pair inside one step.
Between neighbouring nodes and extrema cp is then monotone, and each sign
change of cp - target is one root, polished as in `calc_T_from_cp`. Each
`CpRoot` carries `dT_dcp` = 1 / (dcp/dT) at the root. An extremum whose cp is
within `thermo::CP_ROOTS_DOUBLE_REL` = 1e-13 of the target is returned once
as a double root (`double_root`, `dT_dcp` NaN). Not covered: three or more
extrema within two steps. A non-finite or non-positive `cp_target` or `tol`,
or a range outside [1, 6000] K or with `T_min >= T_max`, throws
`std::invalid_argument`; cp equal to the target at consecutive scan points
(argon's exactly constant cp up to ~950 K: the roots are not isolated) throws
`std::domain_error`; no convergence in `max_iter` throws `std::runtime_error`.

### Species Thermo Beyond the Data (#606)
Every species' NASA-9 data starts at 200 K and ends at 20 000 K (N2, O2, AR,
CO2, CO, H2) or 6000 K (H2O, NH3, CH4, the alkanes): CEA's own ranges. Past
either end, under both thermo models, cp/R leaves the end interval C2 (value,
slope and curvature matched) and levels off; h and s are its exact integrals.

- **Below the data the gas is ideal down to 0 K** (#626): over the whole
  [0, 200 K] a quintic fades cp/R to the rigid-rotor value at 0 K,
  translation plus fully excited rotation (`thermo::CP_R_RIGID_ATOM` 5/2,
  `_LINEAR` 7/2, `_NONLINEAR` 4), as the vibrations freeze out. H2 goes to
  5/2: its rotation freezes inside the fade (Theta_rot = 85 K). The heavy
  alkanes fall monotonically (n-heptane 15.4 -> 4); a near-rigid species
  carries NASA-9's slope at 200 K up to 0.1 past either end (CH4 dips to 3.95
  at 120 K). h(0 K) is finite, so an enthalpy below it is reached by no T: the
  inverses throw. Not modelled: condensation.
- **Above the data** the quartic levels off after `thermo::TAIL_FRACTION_HIGH`
  (0.1) T_hi at the constant it reaches. Solver plumbing, not physics: no
  dissociation.

s_R is exact at every T > 0 (s -> -inf as cp ln T at 0 K) and holds its value
only below `thermo::S_R_T_MIN` (1e-300 K, below any root a solve can reach), so a
non-positive probe T stays finite. h_RT and g_over_RT hold below
`thermo::T_FLOOR` (10 K) (1/T; the equilibrium floors its T there) (#935).

### Complete State Functions
```cpp
struct AirProperties { /* air properties at once */ };
struct ThermoState { /* thermodynamic properties */ };
struct CompleteState { /* thermodynamic + transport properties */ };

AirProperties air_properties(double T, double P, double humidity = 0.0);
ThermoState thermo_state(double T, double P, const std::vector<double>& X, double P_ref = 101325.0);
CompleteState complete_state(double T, double P, const std::vector<double>& X, double P_ref = 101325.0);
```

### Derivatives
```cpp
double dh_dT(double T, const std::vector<double>& X);
double ds_dT(double T, const std::vector<double>& X);
double dcp_dT(double T, const std::vector<double>& X);
double dg_over_RT_dT(double T, const std::vector<double>& X);
```

> [!NOTE]
> The species polynomials switch range at 1000 K. Values (cp, h, s, k, mu)
> are continuous there to ~1e-9, but slopes are one-sided: dcp/dT jumps
> -0.01% for air and -0.4% for lean CH4 products, and dk/dT -0.01% / -0.1%.
> Pure H2 is the outlier, with a 20% jump in its own dcp/dT. The
> data are left unblended: the jump is far below what Newton needs, and a
> blend would change every species' thermo data. Finite-difference checks
> should not straddle 1000 K.

### State-Based Overloads
```cpp
double mwmix(const State& s);
double cp(const State& s);
double h(const State& s);
double s(const State& s);
double cv(const State& s);
double u(const State& s);
double density(const State& s);
double specific_gas_constant(const State& s);
double isentropic_expansion_coefficient(const State& s);
double speed_of_sound(const State& s);
```

---

## Transport Properties (transport.h)

> [!NOTE]
> All transport properties are evaluated at the local (T, P) state.
> **Units**: Viscosity [Pa·s], Conductivity [W/(m·K)], Diffusivity [m²/s].

```cpp
double viscosity(double T, double P, const std::vector<double>& X);
double thermal_conductivity(double T, double P, const std::vector<double>& X);
double prandtl(double T, double P, const std::vector<double>& X);
double kinematic_viscosity(double T, double P, const std::vector<double>& X);
double thermal_diffusivity(double T, double P, const std::vector<double>& X);
double reynolds_from_state(double rho, double v, double L, double mu);

// Viscosity and conductivity with analytic d/dT (chain through the collision
// integral interpolation, the modified Eucken factors and both mixing rules).
struct TransportDerivatives { double mu, dmu_dT, k, dk_dT; };
TransportDerivatives transport_and_dT(double T, double P, const std::vector<double>& X);

// d mu_mix / dX_k at fixed T, analytic (Wilke), per species, X unnormalised
// (chain with a dX/dY that keeps sum X = 1); one-sided for an absent species.
std::vector<double> viscosity_dX(double T, double P, const std::vector<double>& X);

// mu and k with their mole-fraction gradients at fixed T (#717), same X
// convention; k_mix is the mean of the series and parallel bounds.
struct TransportDerivativesX { double mu, k; std::vector<double> dmu_dX, dk_dX; };
TransportDerivativesX transport_and_dX(double T, double P, const std::vector<double>& X);

// composition.h: a mole-fraction gradient as a mass-fraction gradient at X,
// dg/dY_k = (g_X[k] - sum_i X_i g_X[i]) MW / W_k (the step keeps sum X = 1).
std::vector<double> mole_gradient_to_mass(const std::vector<double>& X,
                                          const std::vector<double>& g_X);

// Omega*(2,2) and dOmega*/dT* (C2 quintic Hermite, see Transport Model below).
struct Omega22Result { double omega, domega_dTstar; };
Omega22Result omega22_and_derivative(double T_star, double delta_star);

// thermo.h: d(cp/R)/dT of one species [1/K].
double dcp_R_dT(std::size_t species_idx, double T);

// thermo.h: which species thermo every property evaluates (#549).
// NASA9: the published NASA-9 polynomials, C0 in cp at the interval
// joins but not C1 (dcp/dT jumps 0.15-100 % at 1000 K, up to 48 % at 6000 K).
// NASA9_C2: each join bridged inside |T - T_join| < 0.05 T_join by a degree-7
// cp/R, C2 with NASA-9 at both window edges, its h and s rises across the
// window equal to NASA-9's. Identical to NASA-9 outside the windows; inside,
// cp within 4.8e-4, h within 1.1e-5 cp T, s within 3.1e-5 R. The bridges are
// generated (thermo_data_generator/generate_c2_refit.py ->
// include/thermo_c2_refit_data.h). NASA9_C2 is the default (#549). Process-wide.
enum class ThermoModel { NASA9, NASA9_C2 };
void set_thermo_model(ThermoModel model);
ThermoModel thermo_model();
```

### Transport Model

Viscosity uses Chapman-Enskog kinetic theory with the Omega*(2,2) collision
integral from `include/collision_integral_data.h` (41 T* x 8 delta*, generated
by `thermo_data_generator/collision_integrals.py`), interpolated by quintic
Hermite in (ln T*, ln Omega*) with TABULATED node slopes and curvatures -- C2
in T, so dmu/dT and dk/dT are analytic and themselves continuously
differentiable (#611; the C1 cubic Hermite of #485 stepped d2mu/dT2 by
0.1-10% at every node) -- and linearly in delta* (constant per species). Mixture viscosity uses the Wilke mixing rule.

- **delta* = 0** (Lennard-Jones): values, slopes and curvatures from quadrature of the
  classical collision integrals; agrees with Kim & Monroe (2014) to <= 1e-4
  from T* = 0.3 to 400, between nodes included. The Monchick-Mason column it
  replaces was off by up to 0.6% at the nodes, and linear interpolation in raw
  T* added up to 2.7% between them (#485).
- **delta* > 0** (Stockmayer, polar species H2O, NH3, CO): Monchick-Mason
  values scaled by the delta* = 0 column's exact/tabulated ratio, slopes
  and curvatures from a natural cubic spline (so the interpolant is that
  spline). No exact polar reference is computed.
- Outside T* 0.1-400: continued linearly in (ln T*, ln Omega*) with the end
  slope; the end curvature is set to zero, so the join is C2.

Thermal conductivity uses the Mason-Monchick modified Eucken formula with
Parker Z_rot temperature correction, matching the Cantera GasTransport model.
Its vibrational heat capacity cv_vib = max(0, cp/R - 5/2 - cv_rot/R) is
clamped by a C2 smooth positive part, exact outside |x| < 0.002
(`CV_VIB_SMOOTH_HALF_WIDTH`, in units of R); the hard clamp stepped dk/dT by
up to 4% (H2 355 K, CH4 178 K, O2 189 K; #611).

Species transport parameters are stored in `transport_props` (thermo_transport_data.h):

```cpp
struct Transport_Props {
    MolecularGeometry geometry;  // Atom, Linear, or Nonlinear
    double well_depth;           // epsilon/k_B [K]
    double diameter;             // sigma [Ang]
    double polarizability;       // alpha [Ang^3]; 0.0 for non-polar
    double dipole_moment;        // mu [Debye]; 0.0 for non-polar
    double z_rot;                // rotational relaxation collision number [-]
};
```

---

## Combustion (combustion.h)

> [!NOTE]
> Complete combustion: every C/H fuel burns to CO2 + H2O(g) (fuel N to N2),
> no dissociation. Rich, all fuels burn the same fraction f = X_O2 / n_req
> (= 1/phi) and the rest passes through unburnt; no CO or H2 is formed (use
> `combustion_equilibrium` for that). `include/combustion.h` is the full list.
>
> Smoothing (`smooth_phi0`, `smooth_phi1`, k0 = k1 = 20000): `smooth_phi1`
> replaces min(1, r) with a softmin about 1/k1 wide at phi = 1;
> `smooth_phi0` replaces the O2 demand a with
> n_req = [ln(1 + e^(k0 a)) - ln2 sech(k0 a)] / k0, a transition about 1/k0
> wide at zero fuel: n_req(0) = 0, half the fuel burns as phi -> 0, and
> n_req = a to round-off once k0 a > ~40, so T_ad matches the unsmoothed
> value there (#1037; the softplus before it left ln2/k0 = 3.5e-5 mol O2 per
> mol unburnt at every phi, -0.3 to -0.6 K on T_ad).

### Stoichiometry & Equivalence Ratio
```cpp
double oxygen_required_per_mol_fuel(std::size_t fuel_index);   // mol O2/mol, x + y/4 - z/2 (>= 0)
double oxygen_required_per_kg_fuel(std::size_t fuel_index);    // kg O2/kg
double oxygen_required_per_mol_mixture(const std::vector<double>& X);
double oxygen_required_per_kg_mixture(const std::vector<double>& X);
double dryair_required_per_mol_fuel(std::size_t fuel_index);   // and _per_kg_fuel, _per_*_mixture

// Elemental phi = (2 n_C + n_H / 2) / n_O over ALL atoms of the mixture
// (CO2/H2O count, so a burnt mixture keeps its phi; dry air reads 0.0019 from
// its CO2). 0 without C/H; 1e9 without O.
double equivalence_ratio(const std::vector<double>& X);
double equivalence_ratio_mass(const std::vector<double>& Y);

// Two-stream phi = (F/O) / (F/O)_st, (F/O)_st from the fuel stream's O2 demand
// and the oxidizer stream's O2 (O2 in the fuel stream, fuel in the oxidizer
// stream are not counted).
double equivalence_ratio_mole(const std::vector<double>& X_mix,
                              const std::vector<double>& X_fuel,
                              const std::vector<double>& X_ox);
double equivalence_ratio_mass(const std::vector<double>& Y_mix,
                              const std::vector<double>& Y_fuel,
                              const std::vector<double>& Y_ox);
std::vector<double> set_equivalence_ratio_mole(double phi, const std::vector<double>& X_fuel,
                                               const std::vector<double>& X_ox);  // phi > 0, finite
std::vector<double> set_equivalence_ratio_mass(double phi, const std::vector<double>& Y_fuel,
                                               const std::vector<double>& Y_ox);

// Bilger mixture fraction (two-stream, mass basis)
double bilger_beta(const std::vector<double>& Y);
double bilger_mixture_fraction(const std::vector<double>& Y, const std::vector<double>& Y_F,
                               const std::vector<double>& Y_O);   // clamped to [0, 1]
double bilger_mixture_fraction_from_moles(const std::vector<double>& X,
                                          const std::vector<double>& X_F,
                                          const std::vector<double>& X_O);
double bilger_stoich_mixture_fraction_mass(const std::vector<double>& Y_F,
                                           const std::vector<double>& Y_O);
double bilger_Z_from_equivalence_ratio_mass(double phi, const std::vector<double>& Y_F,
                                            const std::vector<double>& Y_O);
double equivalence_ratio_from_bilger_Z_mass(double Z, const std::vector<double>& Y_F,
                                            const std::vector<double>& Y_O);

// LHV, H2O as vapour [J/mol fuel mixture], [J/kg fuel mixture]
double fuel_lhv_molar(const std::vector<double>& X_fuel, double reference_temperature = 298.15);
double fuel_lhv_mass(const std::vector<double>& X_fuel, double reference_temperature = 298.15);
```

### Complete Combustion Solvers
```cpp
// T, P finite and > 0, X finite (std::invalid_argument otherwise).
State complete_combustion(const State& in, bool smooth_phi0 = false, bool smooth_phi1 = false,
                          double k0 = SMOOTHING_K_PHI0, double k1 = SMOOTHING_K_PHI1);
State complete_combustion_isothermal(const State& in, bool smooth_phi0 = false,
                                     bool smooth_phi1 = false, double k0 = SMOOTHING_K_PHI0,
                                     double k1 = SMOOTHING_K_PHI1);
std::vector<double> complete_combustion_to_CO2_H2O(const std::vector<double>& X,
                                                   bool smooth_phi0 = false, bool smooth_phi1 = false,
                                                   double k0 = SMOOTHING_K_PHI0,
                                                   double k1 = SMOOTHING_K_PHI1);
// ... and an overload with `double& fuel_burn_fraction` after X.

// The same products with the analytic Jacobian dX_b/dX (#716): one
// implementation with the value. Smoothed, the phi = 1 softmin applies at
// every r and the O2 left is X_O2 - f n_req, so the band conserves atoms.
struct CompleteCombustionJacobian { std::vector<double> X_b; std::vector<std::vector<double>> dXb_dX; };
CompleteCombustionJacobian complete_combustion_products_and_jacobian(
    const std::vector<double>& X, bool smooth_phi0 = true, bool smooth_phi1 = true,
    double k0 = SMOOTHING_K_PHI0, double k1 = SMOOTHING_K_PHI1);
```

### Streams and CombustionState
```cpp
Stream set_fuel_stream_for_phi(double phi, const Stream& fuel, const Stream& oxidizer);
// Brent inverse solvers on complete combustion (targets finite; out-of-range -> invalid_argument):
// set_{fuel,oxidizer}_stream_for_{Tad,O2,O2_dry,CO2,CO2_dry}(target, fuel, oxidizer, tol, max_iter[, lean, phi_max])

// Products at T_ad; fuel_burn_fraction = min(1, 1/phi) for Complete.
// _from_streams: T and P of mix({fuel, ox}) (the lowest inlet P), phi from the mass flows.
CombustionState combustion_state(const std::vector<double>& X_fuel, const std::vector<double>& X_ox,
                                 double phi, double T_reactants, double P,
                                 const std::string& fuel_name = "",
                                 CombustionMethod method = CombustionMethod::Complete,
                                 bool smooth_phi0 = false, bool smooth_phi1 = false,
                                 double k0 = SMOOTHING_K_PHI0, double k1 = SMOOTHING_K_PHI1);
CombustionState combustion_state_from_streams(const Stream& fuel_stream, const Stream& ox_stream,
                                              const std::string& fuel_name = "",
                                              CombustionMethod method = CombustionMethod::Complete,
                                              bool smooth_phi0 = false, bool smooth_phi1 = false,
                                              double k0 = SMOOTHING_K_PHI0,
                                              double k1 = SMOOTHING_K_PHI1);
// Both also take a PressureLossCorrelation after `method`: P_out = P (1 - xi(ctx)).

struct PressureLossContext {
    const State& state_in;                  // inlet T, P, X
    double phi, T_ad;
    const std::vector<double>& X_products;
    const std::vector<double>& Y_products;
    double theta;                           // T_ad / T_in - 1
    double mdot_fuel, mdot_air, mdot_total; // 0 in the phi-based call
};

using PressureLossCorrelation = std::function<std::tuple<double, double>(const PressureLossContext &)>;
```

---

## Combustor Pressure Loss (combustor_loss.h)

The built-in combustor losses as a fraction of the inlet total pressure,
`Pt_out = Pt_in (1 - xi)`, with every partial (#718). The Python classes in
`combaero.network.pressure_loss` call these, and `PressureLossElement` builds
its Jacobian from the partials.

```cpp
struct CombustorLossResult {
    double xi, dxi_dtheta, dxi_dmdot, dxi_dTt_in, dxi_dPt_in;
    std::vector<double> dxi_dY;  // per inlet mass fraction; empty for the fraction family
};
// xi = (k theta + xi0) g(|mdot|/m_ref), g = 1 - (1 - r^2)^3 below r = 1; m_ref <= 0: g = 1
CombustorLossResult combustor_fraction_loss(double theta, double mdot, double k,
                                            double xi0, double m_ref);
// xi = (k theta + zeta0) mdot^2 R Tt_in / (2 area^2 Pt_in^2 MW(Y_in)): Lefebvre's pressure-loss
// factor, total-total at low Mach (the dynamic head at the inlet stagnation state, over Pt_in)
CombustorLossResult combustor_head_loss(double theta, double mdot, double Tt_in, double Pt_in,
                                        const std::vector<double>& Y_in,
                                        double k, double zeta0, double area);
```

---

## Chemical Equilibrium (equilibrium.h)

```cpp
// Combustion equilibrium (adiabatic): complete combustion with the phi = 1
// softmin always on (#1048; T is C1 through stoichiometric), then reforming + WGS
EquilibriumResult combustion_equilibrium(const State& in, bool smooth_phi0 = true,
                                         double k0 = SMOOTHING_K_PHI0);

// With analytic derivatives (#716): one implicit-function solve on the joint
// extents (reforming per hydrocarbon, WGS; #759) and the adiabatic root; a
// reaction with an absent species keeps a zero derivative.
struct EquilibriumJacobian {
    State state; double dT_dT_in, dT_dP; std::vector<double> dT_dX;
    std::vector<std::vector<double>> dXout_dX; std::vector<double> dXout_dT_in, dXout_dP;
};
EquilibriumJacobian combustion_equilibrium_and_jacobian(const State& in, bool smooth_phi0 = true,
                                                        double k0 = SMOOTHING_K_PHI0);

// Shift and reforming reactions
State wgs_equilibrium(const State& in);
State wgs_equilibrium_adiabatic(const State& in);
State reforming_equilibrium(const State& in);
```

---

## Stagnation / Static Conversions (stagnation.h)

```cpp
double Tt_from_static(double T, double M, const std::vector<double>& X);
double Tt_from_static_v(double T, double v, const std::vector<double>& X);
double Pt_from_static(double P, double T, double M, const std::vector<double>& X);
double T_from_stagnation(double Tt, double M, const std::vector<double>& X);
double P_from_stagnation(double Pt, double Tt, double M, const std::vector<double>& X);
// Removed in 0.9.0 (#782): the T0_from_static, T0_from_static_v,
// P0_from_static forwarders (use the Tt/Pt names above), solve_P0_from_mdot
// (compressible.h; use solve_Pt_from_mdot) and T0_/P0_from_static_and_jacobian_M
// (solver_interface.h; use Tt_/Pt_from_static_and_jacobian_M).

double T_adiabatic_wall(double T_static, double v, double T, double P,
                        const std::vector<double>& X, bool turbulent = true);
double recovery_factor(double Pr, bool turbulent = true);

double bulk_velocity(double m_dot, double rho, double area);
double kinetic_energy(double v);

std::tuple<double, double> stagnation_from_static(
    double T, double P, double v,
    const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

// Static state of a duct station from stagnation Tt, static P and mass flux
// G = |m_dot|/A (#576): T from h(Tt) = h(T) + (G R T / P)^2 / 2 (unique root),
// lnPtP = ln(Pt/P) = (s(Tt) - s(T))/R, Mach, and partials in Tt, P and G.
struct ChamberStaticState {
    double T, dT_dTt, dT_dP, dT_dG;
    double lnPtP, dln_dTt, dln_dP, dln_dG;
    double M;
};
ChamberStaticState chamber_static_state(double Tt, double P, double G,
                                        const std::vector<double>& X);
// Its mole-fraction gradients at fixed Tt, P, G (#743): dT_dX, dlnPtP_dX
// (map with mole_gradient_to_mass). Implicit function theorem on the energy
// balance; closed-form species h, s and R.
struct ChamberStaticStateDX { std::vector<double> dT_dX, dlnPtP_dX; };
ChamberStaticStateDX chamber_static_state_dX(double Tt, double P, double G,
                                            const std::vector<double>& X);

// The SONIC static state of flux G at stagnation Tt (#637): T* from
// h(Tt) - h(T*) = a*^2 / 2, P* = G sqrt(R T* / gamma*), with dT*/dTt,
// dP*/dTt and dP*/dG. A chamber discharging into a reservoir cannot sit
// below P*.
struct SonicStaticState { double T, P, dT_dTt, dP_dTt, dP_dG; };
SonicStaticState sonic_static_state(double Tt, double G, const std::vector<double>& X);

// Static state at a duct inlet fed from stagnation (Tt, Pt) at mass flux G
// (#636): chamber_static_state at the static P whose isentropic stagnation
// pressure is Pt (bracketed Newton on the subsonic branch [P*, Pt]). Above
// the choked flux: the sonic state, choked = true, no G dependence. Tt and Pt
// are floored at thermo::T_FLOOR / thermo::P_FLOOR (partials 0 there). The
// network's duct heat transfer is evaluated at this state, not at the feed
// node's (stagnation) T.
struct DuctStaticState {
    double T, P, rho, M;
    double dT_dTt, dT_dPt, dT_dG;
    double dP_dTt, dP_dPt, dP_dG;
    double drho_dTt, drho_dPt, drho_dG;
    bool choked;
    // The mass flux the state was evaluated at, and its partials: G itself
    // (dGeff_dG = 1) from duct_static_state; the bent flux from
    // duct_htc_static_state.
    double G_eff, dGeff_dTt, dGeff_dPt, dGeff_dG;
};
DuctStaticState duct_static_state(double Tt, double Pt, double G,
                                  const std::vector<double>& X);

// duct_htc_static_state's mass-fraction gradients at fixed Tt, Pt, G (#717):
// the implicit energy/entropy balances, plus the near-choke bent flux.
struct DuctStaticStateDY { std::vector<double> dT_dY, dP_dY, drho_dY, dGeff_dY; };
DuctStaticStateDY duct_htc_static_state_dY(double Tt, double Pt, double G,
                                           const std::vector<double>& X,
                                           const DuctStaticState& state);  // the solved htc state
```

---

## Compressible Flow (compressible.h)

> [!IMPORTANT]
> All compressible solvers assume **ideal-gas** behavior with variable properties (gamma and Cp are recalculated at each state).
> **Choking**: Functions return a `choked` flag. If true, the mass flow is limited by the sonic condition at the throat.

### Results Structs
```cpp
struct CompressibleFlowSolution {
    State stagnation;
    State outlet;
    double v;
    double M;
    double mdot;
    bool choked;
};

struct FannoSolution {
    State inlet, outlet;
    double mdot, h0, L, D, f_avg;
    bool choked;
    double L_choke;
    std::vector<FannoStation> profile;
};
```

### Solvers
```cpp
// The isentropic static T goes to its root however small (#935; it was held
// at thermo::T_FLOOR = 10 K, not a root, below P/Pt ~ 1e-5, and before that at
// 200 K, which zeroed the flow from a stagnation state at or below 200 K, #552)
// and has no upper bound (a 6000 K cap read a hotter stream's density at
// 6000 K, #626); likewise the stagnation conversions (stagnation.h). Below
// thermo::T_LOG_STEP the Newton steps in ln T.
CompressibleFlowSolution nozzle_flow(double Tt, double Pt, double P_back,
                                     double A_eff, const std::vector<double>& X,
                                     double tol = 1e-8, std::size_t max_iter = 50);

FannoSolution fanno_channel(double T_in, double P_in, double u_in, double L, double D,
                          double f, const std::vector<double>& X,
                          std::size_t n_steps = 100, bool store_profile = false);

FannoSolution fanno_channel_rough(double T_in, double P_in, double u_in, double L, double D,
                                double roughness, const std::vector<double>& X,
                                const std::string& correlation = "haaland",
                                std::size_t n_steps = 100, bool store_profile = false);
```

### Quasi-1D Nozzle Flow
```cpp
// Area function type: A(x) returning area [m²] at position x [m]
using AreaFunction = std::function<double(double)>;

struct NozzleStation {
    double x, A, P, T, rho, u, M, h;
};

struct NozzleSolution {
    State inlet, outlet;
    double mdot, h0, Tt, Pt;
    bool choked;
    double x_throat, A_throat;
    std::vector<NozzleStation> profile;
};

NozzleSolution nozzle_quasi1d(double Tt, double Pt, double P_exit,
                             const AreaFunction& area_func,
                             double x_start, double x_end,
                             const std::vector<double>& X,
                             std::size_t n_stations = 100);

NozzleSolution nozzle_quasi1d(double Tt, double Pt, double P_exit,
                             const std::vector<std::pair<double, double>>& area_profile,
                             const std::vector<double>& X,
                             std::size_t n_stations = 100);

NozzleSolution nozzle_cd(double Tt, double Pt, double P_exit,
                         double A_inlet, double A_throat, double A_exit,
                         double x_throat, double x_exit,
                         const std::vector<double>& X,
                         std::size_t n_stations = 100);
```

### Inverse Solvers
```cpp
double solve_A_eff_from_mdot(double Tt, double Pt, double P_back, double mdot_target,
                             const std::vector<double>& X,
                             double tol = 1e-8, std::size_t max_iter = 50);

double solve_P_back_from_mdot(double Tt, double Pt, double A_eff, double mdot_target,
                             const std::vector<double>& X,
                             double tol = 1e-8, std::size_t max_iter = 50);

double solve_Pt_from_mdot(double Tt, double P_back, double A_eff, double mdot_target,
                         const std::vector<double>& X,
                         double tol = 1e-8, std::size_t max_iter = 50);
```

### Utility Functions
```cpp
double critical_pressure_ratio(double Tt, double Pt, const std::vector<double>& X,
                               double tol = 1e-8, std::size_t max_iter = 50);

double mach_from_pressure_ratio(double Tt, double Pt, double P,
                               const std::vector<double>& X,
                               double tol = 1e-8, std::size_t max_iter = 50);

double mass_flux_isentropic(double Tt, double Pt, double P,
                             const std::vector<double>& X,
                             double tol = 1e-8, std::size_t max_iter = 50);

double fanno_max_length(double T_in, double P_in, double u_in,
                       double D, double f, const std::vector<double>& X,
                       double tol = 1e-6, std::size_t max_iter = 100);
```

Inner solves that reach their iteration cap without converging throw
`std::runtime_error` naming the function and the last residual, rather than
return the last iterate (#868): the isentropic and Fanno-march static-T
Newtons, the Fanno marches themselves (step budget), `Tt_from_static_v`,
`T_from_stagnation`, `chamber_static_state`, `sonic_static_state`, the merge
face fixed point, `friction_colebrook`, and the bracketed solves (Fanno
`static_T`, `fanno_inlet_mach`, area-change, `pipe_equivalent_c`,
`jobson_expansibility`) unless their bracket or last step has closed to
`inner_solve::kCapAcceptRel` (1e-10, `inner_solve.h`). The static-T solves have
no floor: a root below the former 10 K `T_FLOOR` is returned as such (#935). The Fanno-Mach `brent`
throws when called without a sign change. Audit and the solves left as they
are: `docs/archive/AUDIT_868_INNER_SOLVES.md`.

The inverse solvers refuse non-positive inputs (`std::invalid_argument`);
`solve_P_back_from_mdot` throws `std::runtime_error` above the choked flow and
`solve_A_eff_from_mdot` when `P_back >= Pt`. `fanno_max_length` is
`fanno_length_between(G, Tt, M_in, 1, X, D, 0, "fixed", f)` at the inlet flux
and stagnation temperature (`tol`, `max_iter` unused) and throws for
`u_in <= 0` or `f <= 0`.

### Fanno flow integrated in Mach number (`fanno_mach.h`)

At fixed mass flux G and stagnation temperature Tt the station state is
algebraic in M (energy `h(T) + M^2 a^2/2 = h(Tt)`; `u = M a`, `rho = G/u`,
`P = rho R T`), so momentum leaves one smooth integral,
`dx/dM = -(P' + G u') 2D / (f G u)`, which is 0 at M = 1. L* is exact (no
cutoff Mach, no gradient floor); the integral is a fixed 4-panel x 12-point
Gauss-Legendre rule in `s = 1/M^2`, smooth in its arguments.

```cpp
FannoMachState fanno_state_at_mach(double G, double Tt, double M, const std::vector<double>& X);
double fanno_dx_dmach(double G, double Tt, double M, const std::vector<double>& X, double D,
                      double roughness, const std::string& friction_model,
                      double f_multiplier = 1.0);  // "fixed": constant Darcy f = f_multiplier
double fanno_length_between(double G, double Tt, double M1, double M2, ...);  // [m]
double fanno_inlet_mach(double Pt, double Tt, double G, const std::vector<double>& X);
double fanno_sonic_mass_flux(double Pt, double Tt, const std::vector<double>& X);
FannoDuctResult fanno_duct(double Pt, double Tt, double G, const std::vector<double>& X,
                           double L, double D, double roughness,
                           const std::string& friction_model, double f_multiplier = 1.0);
double fanno_choked_mass_flux(double Pt, double Tt, const std::vector<double>& X, double L,
                              double D, double roughness, const std::string& friction_model,
                              double f_multiplier = 1.0);  // G with L* = L

// Mass flux from (Pt0, Tt0) against P_target -- the duct's exit stagnation
// pressure if exit_total, else its exit static -- with dG/dPt0, dG/dTt0,
// dG/dP_target. Monotone in P_target, saturating exactly at the choked flux
// (choked = true, dG/dP_target = 0); 0 at or below zero drive, and blended C1
// to a finite slope below kFannoFlowDriveFloor = 1 Pa, where M_in and M_exit
// are re-evaluated at the blended G (#1017) and M_exit_root keeps the floor
// root's exit Mach (the warm start for M_exit_guess). Unchoked it solves for
// the exit Mach (P_exit is linear in G at fixed M, so G(M_exit) is explicit)
// and differentiates implicitly, well conditioned up to the choke.
// dY_directions: mass-fraction directions v; dG_ddir returns dG/ds along each
// (Y -> Y + s v, X renormalised), about two length integrals per direction. A
// unit vector e_k gives dG/dY_k at fixed others. ChannelElement passes an
// orthonormal basis of its feeding node's dY/dx -- usually one vector -- and
// projects back onto "{node}.Y[k]" keys. with_derivatives = false computes G,
// choked and the Mach numbers only (a residual-only evaluation).
FannoChannelFlow fanno_channel_flow(double Pt0, double Tt0, const std::vector<double>& X,
                                    double P_target, bool exit_total, double L, double D,
                                    double roughness, const std::string& friction_model,
                                    double f_multiplier = 1.0, double M_exit_guess = -1.0,
                                    const std::vector<std::vector<double>>& dY_directions = {},
                                    bool with_derivatives = true);
```

### Fanno-Rayleigh: friction and heat (`fanno_mach.h`, #551)

The same duct with Tt rising linearly from `Tt0` to `Tt_exit`. Exactly, for a
thermally perfect gas, `P_M + G u_M = -G a cv (1 - M^2) / (k M^2)` with
`k = cp(T) + M^2 a da/dT`, so momentum gives `dM/dx` with a heat term; heat
moves M either way, so the march is in x, of `theta = M^2 - M^4/2`
(`(1 - zeta)/2`, `zeta = (1 - M^2)^2`), regular at sonic and full relative
precision at low Mach. Fixed RK4, 32 steps graded toward the exit.

The choked flux is set by where the duct goes sonic: at its **exit** (a
backward march from M = 1 at L meets the isentropic inlet -- smooth), at an
interior **throat** (cooling decelerates the flow after a sonic point; its
derivatives re-solve the choke), or at its **inlet**. Checked against the
Fanno kernel at zero heat (<= 1e-7 in G) and Rayleigh algebra at zero
friction (<= 1e-8 in exit Mach). About 8 ms per warm call with derivatives,
~3-5x the adiabatic kernel.

```cpp
FannoRayleighMarch fanno_rayleigh_march(double G, double Pt0, double Tt0, double Tt_exit,
                                        const std::vector<double>& X, double L, double D,
                                        double roughness, const std::string& friction_model,
                                        double f_multiplier = 1.0);
// {M_in, M_exit, zeta_min (<= 0: choked), zeta_peak, P_exit, Pt_exit, T_exit}
double fanno_rayleigh_choked_mass_flux(double Pt0, double Tt0, double Tt_exit, ...);
// fanno_channel_flow's contract plus dG_dTt_exit; G_guess warm-starts the root.
FannoChannelFlow fanno_rayleigh_channel_flow(double Pt0, double Tt0, double Tt_exit,
                                             const std::vector<double>& X, double P_target,
                                             bool exit_total, double L, double D,
                                             double roughness, const std::string& friction_model,
                                             double f_multiplier = 1.0, double G_guess = -1.0,
                                             const std::vector<std::vector<double>>& dY_directions = {},
                                             bool with_derivatives = true);
```

The flux is the root on the branch continuous from zero flow, 0 at or below
zero drive. Limit: a cooled duct's Rayleigh gain can outgrow its friction --
not under Reynolds analogy (NTU ~ f L/2D keeps it below half), but in a short
duct with entrance-enhanced heat transfer and fully developed friction (an
impingement crossflow segment, L/D 1.3). At a fixed Tt_exit its exit Pt then
rises with flow and no flow matches a pressure drop; developing-flow friction
is the consistent closure.

### Rocket Nozzle Thrust
```cpp
struct ThrustResult {
    double thrust;
    double specific_impulse;
    double thrust_coefficient;
    double mdot;
    double u_exit;
    double P_exit;
};

ThrustResult nozzle_thrust(const NozzleSolution& sol, double P_amb);

ThrustResult nozzle_thrust(double Tt, double Pt, double P_design, double P_amb,
                           double A_inlet, double A_throat, double A_exit,
                           double x_throat, double x_exit,
                           const std::vector<double>& X,
                           std::size_t n_stations = 100);
```

---

## Incompressible Flow (incompressible.h)

```cpp
double bernoulli_P2(double P1, double v1, double v2, double rho, double dz = 0.0);
double orifice_mdot(double P1, double P2, double A, double Cd, double rho);
double channel_dP(double v, double L, double D, double f, double rho);
```

---

## Area Change (area_change.h)

Sudden and conical expansions and contractions (#731). The closures are
Idelchik's (1966; diagrams 4-1, 3-6/3-7, 3-9, 3-10, 5-2, transcribed in
`validation/area_change/`), every zeta on the small section's dynamic head
and Re there; n is the small area over the large one. Mach enters through the
control volume the closures come from: an expansion is the compressible
Borda-Carnot step times the table's ratio to it, a contraction a vena
contracta (Cc = 1 / (1 + sqrt(zeta_c))) re-expanding by the same step. The
header documents the interpolation, its monotonicity regularisation and the
drive floor; the measurements are in `python/tests/test_area_change.py`.

```cpp
struct ZetaAndSlope { double zeta, dzeta_dRe; };
ZetaAndSlope sudden_expansion_zeta(double Re, double area_ratio);  // 4-1, C1 in Re
ZetaAndSlope sudden_contraction_zeta(double Re, double area_ratio);  // 3-10, 3-9, C1 in Re
double conical_diffuser_zeta(double full_angle_deg, double area_ratio);  // 5-2
double conical_nozzle_zeta(double full_angle_deg, double length_over_D,
                           double area_ratio);                // 3-7 with 3-6

struct BordaCarnotStep { double M2, pt_ratio, dpt_dM1, dpt_dgamma; };
BordaCarnotStep compressible_borda_carnot(double M1, double area_ratio, double gamma);

// Loss at a given flow (diagnostics), fed from the stagnation state upstream.
struct AreaChangeLoss { double dPt, zeta, mach, Re; bool expansion, choked; };
AreaChangeLoss area_change_loss(double m_dot_abs, double Pt0, double Tt0,
                                const std::vector<double>& Y, double F_in, double F_out,
                                const std::string& model, double length, double D_h,
                                bool incompressible);

// Flow at a given drive Pt0 - Pt_target (odd in it), saturating at the sonic
// flow, with analytic partials (dm_dY per species).
struct AreaChangeFlow {
    double m_dot; bool choked; double mach, Re, zeta;
    double dm_dPt0, dm_dTt0, dm_dPt_target;
    std::vector<double> dm_dY;
};
AreaChangeFlow area_change_flow(double Pt0, double Tt0, const std::vector<double>& Y,
                                double Pt_target, double F_in, double F_out,
                                const std::string& model, double length, double D_h,
                                bool incompressible, bool with_derivatives = true);

// The network element's flow from node a (face F_a) to node b: each direction
// fed from its own upstream node, the two blended C1 through zero drive.
struct AreaChangeElementFlow {
    double m_dot; bool choked; double mach, Re, zeta;
    double dm_dPt_a, dm_dTt_a, dm_dPt_b, dm_dTt_b;
    std::vector<double> dm_dY_a, dm_dY_b;
};
AreaChangeElementFlow area_change_element_flow(
    double Pt_a, double Tt_a, const std::vector<double>& Y_a,
    double Pt_b, double Tt_b, const std::vector<double>& Y_b,
    double F_a, double F_b, const std::string& model, double length, double D_h,
    bool incompressible, bool with_derivatives = true);

// Incompressible kernels at a given density (standalone): dP and its partials.
struct AreaChangeResult { double dP, dS_dm, dS_drho, dS_dmu; };
AreaChangeResult sharp_area_change(double m_dot, double rho, double mu, double F0, double F1,
                                   double m_scale = 1e-4, double D_h = 0.0);
AreaChangeResult conical_area_change(double m_dot, double rho, double mu, double F0, double F1,
                                     double length, double m_scale = 1e-4);
```

---

## Friction Factor Correlations (friction.h)

```cpp
double friction_haaland(double Re, double e_D);
double friction_serghides(double Re, double e_D);
double friction_colebrook(double Re, double e_D, double tol = 1e-10, int max_iter = 20);
double friction_petukhov(double Re);
// Petukhov held below Re 3000 by a C1 soft-max over Re 2500-3500 (#446);
// used by every smooth-pipe Gnielinski / Petukhov heat-transfer path.
double friction_petukhov_clamped(double Re);
double friction_petukhov_clamped_dRe(double Re);
```

The functions above are the turbulent fits on their own. The network kernels
`solver::friction_and_jacobian(model, Re, e_D)` ("haaland", "serghides",
"colebrook") blend them with Hagen-Poiseuille, `f = (1 - w) 64/Re + w f_turb`
with `w = (1 + tanh((Re - 2300)/400))/2` (`laminar_friction` in friction.h).
The laminar branch carries its own floor of Re 1e-6, so a channel's drop stays
linear in its flow down to essentially zero flow.

---

## Heat Transfer Correlations (heat_transfer.h)

> [!NOTE]
> **Units**: HTC [W/(m²·K)], q [W/m²].
> **Jacobians**: High-performance solver components return `dh/dmdot` and other derivatives essential for coupled fluid-thermal networks.

### Nusselt Numbers
```cpp
double nusselt_dittus_boelter(double Re, double Pr, bool heating = true);
double nusselt_gnielinski(double Re, double Pr);
double nusselt_gnielinski(double Re, double Pr, double f);
// Value and EXACT Re-derivative (no warnings); df_dRe = slope of the f passed.
NuAndDerivative nusselt_gnielinski_with_derivative(double Re, double Pr,
                                                   double f, double df_dRe = 0.0);
NuAndDerivative nusselt_gnielinski_smooth_with_derivative(double Re, double Pr);

// Channel regime transitions (#448): every switch is a C1 smoothstep, each
// regime EXACT outside its band. channel_smooth, nusselt_circular_channel and
// htc_circular_channel (gnielinski branch) use these, so neither the value nor
// the Jacobian steps at Re 2300 or 4000.
//   laminar -> turbulent  Re 2300-3000   f: 64/Re -> turbulent; Nu: 3.66/4.36 -> Gnielinski
//   smooth  -> rough      Re 3000-4000   f: clamped Petukhov -> Colebrook (e_D > 0)
// Turbulent-only Nu (dittus_boelter, sieder_tate, petukhov; Re >= 1e4) in
// channel_smooth and htc_circular_channel: the VDI Heat Atlas transition
// (vdi_transition) -- linear in Re from laminar Nu at 2300 to the correlation
// at 1e4, exact derivatives; they threw in that band (#481).
FrictionAndDerivative friction_channel_and_derivative(double Re, double e_D);   // Darcy, exact df/dRe
FrictionAndDerivative friction_turbulent_and_derivative(double Re, double e_D);
NuAndDerivative nusselt_channel_gnielinski_and_derivative(
    double Re, double Pr, double f_turb, double df_turb_dRe, double Nu_laminar);
double friction_colebrook_dRe(double Re, double e_D);  // implicit, exact (friction.h)
double nusselt_sieder_tate(double Re, double Pr, double mu_ratio);
```

### Overall Heat Transfer
```cpp
double overall_htc(const std::vector<double>& h_values, const std::vector<double>& t_over_k);
double overall_htc_wall(double h_inner, double h_outer, const std::vector<double>& t_over_k_layers);
```

### Channel Flow Models
```cpp
// Smooth channel
ChannelResult channel_smooth(double T, double P, const std::vector<double>& X,
                              double velocity, double diameter, double length,
                              double T_hot = std::numeric_limits<double>::quiet_NaN(),
                              const std::string& correlation = "gnielinski",
                              bool heating = true, double Nu_multiplier = 1.0,
                              double f_multiplier = 1.0,
                              double flow_area = NaN);  // NaN: pi D^2/4; used by the mdot derivatives (#463)

// Parametrised rib correlations -- rib_correlation.h
//
//   R = C_R * (e/D / nD)^a * (p/e / nP)^b * (W/H / nW)^c * (alpha/90)^d
//   G = C_G * (same geometry terms) * (e+)^n
//
// R carries no e+ term by construction, so f is independent of Reynolds
// number and df/d(mdot) is identically zero for this family.
//
// han_park_1988_angled uses two shapes the plain power law above cannot
// express: R is a quadratic in alpha (Eq. 4.17), and G's alpha/p_e
// exponents switch on channel shape, square vs rectangular (Eq. 4.18).
// Both are genuine discontinuities the source states -- up to 62.5% in R
// at alpha=90, W/H=4; up to 27% in G at W/H=1 -- and neither is smoothed,
// since no smooth transition is stated. See
// validation/cooling/extractions/han_ribbed.md, decision D5, and
// RibCorrelationSet::RAlphaShape / GShapeModel in rib_correlation.h.
//
// han_1989_narrow_channel (Eq. 4.19, W/H < 1) adds two more shapes. R has
// NO PRINTED EQUATION in this source -- only a drawn figure line -- so its
// quadratic-in-alpha coefficients (QuadraticAlphaTwoBand, one per W/H
// sub-band) are this project's OWN FIT to that line, not an extraction;
// provenance is Fitted. G (NarrowChannelAlphaSwitch) IS text-extracted: a
// hard ~20% switch in its leading constant at alpha == 90 deg, plus a
// further W/H correction to both the constant and the e+ exponent itself
// below W/H = 1/2. See han_ribbed.md item 40 and decision D6.
//
// Rib SHAPE is part of a set (#434): continuous-angled, V, crossed and
// Lambda ribs at the same alpha are different geometries, and a set refuses
// a series of another shape rather than aliasing onto it. Unspecified (the
// default for a user-built set) binds nothing.
enum class RibShape { Unspecified, Transverse, Parallel, Crossed, V, Lambda };
// set.shape; set.C_Gbar / set.Gbar_eplus_exponent carry a source's OWN
// printed four-wall G_bar (0 = not printed -> RibResult::has_G_bar false,
// callers fall back to Han's G_bar = 1.2 G).
//
// Named sets, selected explicitly -- there is no auto-switching between
// them, since they cover disjoint regimes rather than one being an update
// of the other:
RibCorrelationSet han_1988_orthogonal();       // 90 deg, Re 10e3-60e3
RibCorrelationSet rallabandi_2009_high_re();   // 45 deg sharp ribs, Re 30e3-400e3
RibCorrelationSet han_park_1988_angled();      // 30-90 deg, W/H 1-4, Re 10e3-60e3
RibCorrelationSet han_1989_narrow_channel();   // 30-90 deg, W/H 1/4-1, Re 10e3-60e3
// Han, Zhang & Lee (1991) Table 2: Transverse at 90; Parallel, Crossed, V,
// Lambda at 45 or 60. One rig (square, e/D 0.0625, P/e 10); anything else
// throws. Prints its own G_bar per configuration.
RibCorrelationSet han_zhang_lee_1991(RibShape shape, double alpha_deg);
void validate_rib_set(const RibCorrelationSet& set);   // throws on a bad set

// Accuracy carries WHERE IT CAME FROM, because only an author's own claim
// may be used as the band a model is judged against. A figure this project
// measured is that model's own error, and scoring a model inside its own
// error answers nothing (#389).
enum class AccuracyProvenance { Unstated, Stated, Measured };
struct StatedAccuracy {
  double value = NaN;                 // NaN unless provenance says otherwise
  AccuracyProvenance provenance = AccuracyProvenance::Unstated;
  bool usable_as_band() const;        // true only for Stated
  static StatedAccuracy stated(double v);
  static StatedAccuracy measured(double v);
  static StatedAccuracy unstated();
};
// set.accuracy_R / set.accuracy_G are StatedAccuracy. Of the four shipped
// sets only han_1988_orthogonal's are Stated (Han's "95% within 6%/8%");
// han_park_1988_angled's and han_1989_narrow_channel's were measured
// through evaluate_rib itself, and rallabandi_2009_high_re's accuracy_R is
// Unstated -- previously 0.0, which read as perfect agreement.
RibResult evaluate_rib(const RibCorrelationSet& set,
                       const RibGeometry& geom, double Re);
// evaluate_rib never throws: Re may be negative or zero and the guards are
// smooth through both, because the solver probes states that are not physical.

// Ratio-form rib sets -- rib_ratio_correlation.h (#444). Nu/Nu0 and f/f0 as
// multipliers on the baseline the SOURCE fitted to, beside the R/G sets, for
// sources that publish ratios (and so that their data is not forced through
// the law-of-the-wall similarity R/G assumes).
//
//   in range:    Nu = r(Re, geom) * Nu0_source        (reproduces the paper)
//   below floor: C1 blend in ln Re over [Re_floor/2, Re_floor] to
//                k * Nu0_ext, k matching the value at Re_floor;
//                Nu0_ext = smooth-pipe Gnielinski by default (laminar 3.66
//                limit at Re -> 0), the source baseline, or a user power law.
//   friction:    ratio held at its floor value below the range.
//
// Never divides by Nu0; Re enters as sqrt(Re^2 + 1) so Nu, f are even and
// their analytic derivatives odd and finite through Re = 0.
struct RatioBaseline { double coeff, re_exponent, pr_exponent; };  // c Re^m Pr^n
enum class RatioBelowFloor { Gnielinski, SourceBaseline, User };
struct RibRatioOptions { RatioBelowFloor below_floor; RatioBaseline user_Nu0; };
RibRatioResult evaluate_rib_ratio(const RibRatioSet& set, const RibGeometry& geom,
                                  double Re, double Pr,
                                  const RibRatioOptions& options = {});
// -> Nu, dNu_dRe, f, df_dRe, ratio_Nu, ratio_f, below_floor, extrapolated
void validate_rib_ratio_set(const RibRatioSet& set);
void validate_rib_ratio_options(const RibRatioOptions& options);
// set.Re_floor_f: friction's own floor (0 = Re_floor), for sources that
// measured f over a different Re range than Nu.
RibRatioSet taslim_spring_1987(double aspect_ratio_taslim, double e_D);
// Taslim & Spring (1987), two ribbed walls, transverse, p/e 10. Seven tested
// configurations: AR 0.5 (e/D 0.125, 0.250), 1.0 (0.083, 0.167), 3.5 (0.053,
// 0.107, 0.161); any other throws. Taslim's AR is height/width: W/H = 1/AR.
// Nu/Nu_DB = C (Re/1e4)^-0.2 (fig. 9), f Fanning passage-average constant
// (fig. 11). Provenance Fitted; accuracy Unstated.

// Enhanced surfaces removed in 0.7.0 for unprovenanced correlations (issue
// #339). Ribbed is re-added on han_1988_orthogonal / rallabandi_2009_high_re
// above, wired through RibbedModel (components.py) into the GUI and network
// solver. Impingement is re-added below, wired through ImpingementModel /
// SingleJetImpingementModel (components.py) into the GUI and network solver
// the same way. Dimpled and pin-fin remain removed, tracked individually
// (#335-336).
```

### Pin-Fin Correlations
```cpp
// pin_fin_correlation.h (#335). Sets carry each source's own form and basis;
// evaluators take and return one canonical basis:
//   Re_D = mdot D / (mu A_min),  Nu_D = h D / k,  f = dP / (2 rho Vmax^2 N)
// Heat transfer and friction are separate sets. Outside a set's box the
// result is flagged `extrapolated`, never refused; evaluators never throw.
enum class PinArrangement { Staggered, Inline };
enum class PinNuSurface { Total, Endwall, Pin };
enum class PinReBasis { DiameterVmax, VanFossenDprime };
enum class PinFrictionBasis { PerRowVmax, PerRowGapVmax };
struct PinFinGeometry { double S_D, X_D, H_D; int N_rows; PinArrangement arrangement; };

// Exact converters (Armstrong & Winstanley 1988, Eqs 8, 10, 13): D'/D, A'/A_min,
// D_h/D, A_min/A_frontal (min of transverse and staggered diagonal gap), and
// per-wall area fractions with pins to mid-height.
double pin_fin_min_gap_D(const PinFinGeometry& g);
double pin_fin_amin_over_afrontal(const PinFinGeometry& g);
double pin_fin_dprime_over_D(const PinFinGeometry& g);
double pin_fin_aprime_over_amin(const PinFinGeometry& g);
double pin_fin_dh_over_D(const PinFinGeometry& g);
PinFinAreaFractions pin_fin_area_fractions(const PinFinGeometry& g);
// eta_fin = tanh(mL)/(mL), L = H/2, m = sqrt(4h/(k D)); eta_t and d/dh analytic.
PinFinEfficiency pin_fin_array_efficiency(double h, double k_pin, double D,
                                          double H, double A_f_over_A_t);

PinFinNuResult evaluate_pin_fin_nu(const PinFinNuSet& set, const PinFinGeometry& g,
                                   double Re_D, double Pr);         // Nu, dNu_dRe
PinFinFrictionResult evaluate_pin_fin_friction(const PinFinFrictionSet& set,
                                               const PinFinGeometry& g,
                                               double Re_D);        // f, df_dRe
PinFinModifierResult evaluate_pin_fin_modifier(const PinFinRatioModifier& mod,
                                               const PinFinGeometry& g, double Re_D);
void validate_pin_fin_geometry(const PinFinGeometry& g);
void validate_pin_fin_nu_set(const PinFinNuSet& set);
void validate_pin_fin_friction_set(const PinFinFrictionSet& set);
void validate_pin_fin_modifier(const PinFinRatioModifier& mod);

PinFinNuSet metzger_1986_staggered_nu();             // default Nu
PinFinFrictionSet metzger_1982_staggered_friction(); // default f, C1 blend at 1e4
PinFinNuSet vanfossen_1982_staggered_nu();           // D' basis
PinFinFrictionSet damerow_1972_staggered_friction(); // per (N-1) rows, wide S/D
PinFinNuSet chyu_1998_nu(PinArrangement a, PinNuSurface s);
PinFinRatioModifier chyu_1998_inline_over_staggered(); // has_f = false
PinFinNuSet chyu_1990_nu(PinArrangement a, bool fillet);          // Table 2, pin surface
PinFinFrictionSet chyu_1990_friction(PinArrangement a, bool fillet); // Fig. 6, Fitted,
// basis PerRowHalfHeadVmax (f = 2 dP/(rho Vmax^2 N), 4x canonical); the only
// inline friction source in hand.
PinFinRatioModifier chyu_1990_fillet_over_straight(PinArrangement a);
```

### Jet Impingement Correlations
```cpp
// impingement_correlation.h -- two independent regimes, matching
// han_impingement.md's own split, wired into ImpingementModel /
// SingleJetImpingementModel (components.py). Converting an array's
// total/mean mass flow into a per-row local Re_j (Florschuetz's own Eq. 7)
// is deliberately NOT implemented here -- ImpingementModel recovers each
// row's own mass flow from its own ConvectiveSurface.area instead, leaving
// row-to-row flow distribution to whoever assembles the chain of elements.

// Single jet (Goldstein, Behbahani and Heppelmann, 1986):
//   Nu_bar = Re^0.76 * (A - |L/D - 7.75|) / (B + C*(R/D)^n)
// n switches with the boundary condition the correlation was fitted under.
enum class ImpingementThermalBC { ConstantHeatFlux, ConstantWallTemperature };
SingleJetImpingementSet goldstein_1986_single_jet();
void validate_single_jet_set(const SingleJetImpingementSet& set);
double single_jet_impingement_nu(const SingleJetImpingementSet& set,
                                 ImpingementThermalBC bc, double Re,
                                 double L_D, double R_D);
// Same formula, plus an analytic d(Nu)/d(Re) -- for a network element's
// wall-coupling Jacobian, per the project's (f, J) rule.
SingleJetImpingementResult single_jet_impingement(
    const SingleJetImpingementSet& set, ImpingementThermalBC bc, double Re,
    double L_D, double R_D);

// Jet array with crossflow (Florschuetz, Truman and Metzger, 1981):
//   Nu = A * Re_j^m * {1 - B*[(z/d)(Gc/Gj)]^n} * Pr^(1/3)
// A, m, B, n are each C*(xn/d)^nx*(yn/d)^ny*(z/d)^nz (Table 4.1). Re_j and
// Gc/Gj are the ROW's own local values, not an array mean.
// JetArrayImpingementResult::dNu_dRe_j is analytic, holding Gc/Gj and
// geometry fixed -- Gc/Gj never depends on Re_j or mdot, only on geometry.
enum class JetHolePattern { Inline, Staggered };
JetArrayCorrelationSet florschuetz_1981_inline();
JetArrayCorrelationSet florschuetz_1981_staggered();
void validate_jet_array_set(const JetArrayCorrelationSet& set);
JetArrayImpingementResult jet_array_impingement_nu(
    const JetArrayCorrelationSet& set, double Re_j, double Gc_Gj, double Pr,
    double xn_d, double yn_d, double z_d);

// Crossflow-to-jet mass flux ratio, Florschuetz's own closed form (Eq. 8),
// depending on (yn/d)(z/d) only -- not xn/d, per the source. C_D is the
// jet-plate discharge coefficient (FLORSCHUETZ_1981_DEFAULT_CD = 0.79
// absent a measured value; combaero has no correlation of its own for a
// jet-plate array yet -- see #375, distinct from orifice.h's pipe-metering
// Cd family, which does not apply to this geometry).
constexpr double FLORSCHUETZ_1981_DEFAULT_CD = 0.79;
double crossflow_to_jet_ratio_at_x(double yn_d, double z_d, double C_D,
                                   double x_over_xn);
double crossflow_to_jet_ratio_at_row(double yn_d, double z_d, double C_D,
                                     int row);  // row 1 is exactly 0
// Incompressible; measured range channel P/P0 >= 0.72 (NASA CR-3630 Table D.1:
// RMS 3.4% there, as the compressible network chain). Below ~0.7: extrapolation (#474).

// One row from NETWORK flows (#465): the generic core of a jet-plate row.
// Gc/Gj = (m_crossflow/m_jet)(pi/4)/((yn/d)(z/d)), Re_j = 4 m_jet/(n pi d mu).
// Which flows count as crossflow is the caller's configuration.
struct JetRowGeometry { double d, xn_d, yn_d, z_d, n_holes; };
struct JetRowHeatTransfer {
    double h, Nu, Re_j, Gc_Gj;
    double dh_dm_jet;        // analytic, includes the Gc/Gj ~ 1/m_jet term
    double dh_dm_crossflow;  // analytic
    bool extrapolated;
    double eta_r;                // impingement_eta_r at this Gc/Gj
    double deta_r_dm_jet, deta_r_dm_crossflow;  // analytic, through Gc/Gj
    bool eta_r_extrapolated;
};
double jet_row_gc_gj_per_flow_ratio(double yn_d, double z_d);
JetRowHeatTransfer jet_row_heat_transfer(const JetArrayCorrelationSet& set,
                                         const JetRowGeometry& geom,
                                         double m_jet, double m_crossflow,
                                         double mu, double k, double Pr);

// Crossflow temperature influence factor (#467), CR-3936 Eq. 2.11:
// eta_r = (T_ref,r - T_j)/(T_m,n - T_j). FITTED (the source publishes no
// correlation) to its Appendix A, 39 inline tests:
//   eta_r = X^n/(1+X^n), X = a (Gc/Gj)(z/d)^r
// RMSE 0.079 in-sample, 0.083 leave-one-geometry-out. 0 at Gc/Gj <= 0 (C1:
// n > 1), -> 1 as the crossflow dominates. extrapolated: xn/d outside 5-10,
// yn/d 4-8, z/d 1-3, or Gc/Gj > 4.44 (rel_tol 5e-3). Staggered arrays and
// the source's eta_r > 1 are outside the form.
constexpr double FLORSCHUETZ_SU_1985_ETA_A = 2.042;
constexpr double FLORSCHUETZ_SU_1985_ETA_R = 1.073;
constexpr double FLORSCHUETZ_SU_1985_ETA_N = 1.368;
struct ImpingementEtaR { double eta, deta_dGc_Gj; bool extrapolated; };
ImpingementEtaR impingement_eta_r(double Gc_Gj, double xn_d, double yn_d,
                                  double z_d);
```

### Side-Stream Momentum (`side_stream_momentum.h`)
```cpp
// Exact (compressible) impulse of streams joining a channel; each station at
// the density of its own static (P, T, X). Derivatives w.r.t. the solver's
// variables, densities through density_and_jacobians.

// One side-stream station, half at a time (#465, #471). m_a arrives, m_b
// leaves along the channel; the side stream's axial velocity is kappa*u_a:
//   dP = (m_b^2 - m_a^2 + kappa m_a (m_a - m_b)) / (rho A^2)   (m^2: any flow direction, #815)
// kappa = 0: normal injection (merge, Florschuetz). kappa = 0.75: bleed,
// Bassett et al. (2001) K2/K5 Eq. 15 exactly. A segment carries half of the
// station at each end, each at its node's density.
constexpr double STATION_KAPPA_MERGE_NORMAL = 0.0;
constexpr double STATION_KAPPA_BLEED_BASSETT = 0.75;
struct StationHalfDrop { double dP, d_dm_a, d_dm_b, d_dP, d_dT; };
StationHalfDrop station_half_drop(double m_a, double m_b, double P, double T,
                                  const std::vector<double>& X, double area, double kappa);
// Area change over a segment between two static ends (#812), Bernoulli, even in m:
//   m^2 / (2 rho) (1/A_to^2 - 1/A_from^2); zero for equal areas.
StationHalfDrop segment_area_change_drop(double m, double P, double T,
                                         const std::vector<double>& X, double area_from,
                                         double area_to);
// Reservoir entry into a duct: (1 + K_in) m|m| / (2 rho A^2).
StationHalfDrop channel_entry_drop(double m, double P, double T,
                                   const std::vector<double>& X, double area, double K_in);

// Supply-side crossflow ratio of a duct-fed hole (McGreehan-Schotsch U1/Vi):
// U1 = (m_a + m_b)/(2 rho A) at the station, Vi = isentropic velocity from
// the duct's STATIC (P, T) to P_down (floored at CROSSFLOW_VI_FLOOR).
struct CrossflowRatio { double U1_over_Vi, U1, Vi, d_dm_a, d_dm_b, d_dP, d_dT, d_dP_down; };
CrossflowRatio crossflow_velocity_ratio(double m_a, double m_b, double P, double T,
                                        const std::vector<double>& X, double P_down,
                                        double area);

// Wall-hole Cd under a grazing crossflow (#800, crossflow_cd.h), in
// q = G_ideal / G_c (ideal hole over crossflow mass velocity), resolved so
// Cd does not read the hole's own flow. Rogers & Hersh (1976) CR-2681:
// 0.97 (Gj/Gc)^0.534 injecting, 0.80 (Gj/Gc)^0.47 suction, joined C1 to Cd0
// over +/-ROGERS_HERSH_1976_CAP_BAND. Florschuetz, Metzger & Su (1982)
// CR-3630 Table 4.1 / Eq. 4.2 for the plates it covers (throws otherwise),
// its pieces made continuous within the printed digits.
constexpr double ROGERS_HERSH_1976_CAP_BAND = 0.1;
constexpr double ROGERS_HERSH_1976_INJECTION_A = 0.97, ROGERS_HERSH_1976_INJECTION_B = 0.534;
constexpr double ROGERS_HERSH_1976_SUCTION_A = 0.80, ROGERS_HERSH_1976_SUCTION_B = 0.47;
constexpr double FLORSCHUETZ_METZGER_SU_1982_CAP_BAND = 0.1;
constexpr double FLORSCHUETZ_METZGER_SU_1982_ZETA_C = 1.06, FLORSCHUETZ_METZGER_SU_1982_ZETA_X0 = 0.806,
                 FLORSCHUETZ_METZGER_SU_1982_ZETA_N = 0.602;
struct CrossflowCd { double Cd, dCd_dG_ideal, dCd_dG_c; bool capped; };
CrossflowCd rogers_hersh_1976_cd(double Cd0, double G_ideal, double G_c, bool suction);
// CR-2681's envelope by variable and branch (INFO, never a selector):
// Vi/V_inf 0.004-0.5 (Fig. 11); L/D 0.12-2.0 injecting (0.6-1.0 above
// Vi/V_inf 0.35, where the L/D independence is unconfirmed), 0.6-1.0 for
// suction (Sect. 5.2.3: the L/D effect is considerable; Fig. 11 is L/D 1);
// V_inf 24-85 m/s. Constants ROGERS_HERSH_1976_*.
struct RogersHersh1976Envelope { bool velocity_ratio, l_over_d, grazing_velocity; bool in_range() const; };
RogersHersh1976Envelope rogers_hersh_1976_envelope(double Vi_over_V_inf, double L_over_D,
                                                   double V_inf, bool suction);
bool florschuetz_metzger_su_1982_has_plate(double xn_d, double yn_d, double z_d);
CrossflowCd florschuetz_metzger_su_1982_cd(double xn_d, double yn_d, double z_d,
                                           double G_ideal, double G_c);

// EXIT (external) crossflow on a film-cooling hole (#1074, crossflow_cd.h):
// Rowbury, Oldfield & Lock (2001) Eqs. 7-8, an additive outlet loss on the
// no-crossflow Cd0, 1/Cd^2 = 1/Cd0^2 + A k_perp^-B - C, k_perp = rho_h u_h
// u_h,perp / (rho_m u_m^2) on the hole's ACTUAL flow, u_h,perp/u_h =
// 1/sqrt(cot^2 theta + 1 + tan^2 alpha) (Fig. 3, alpha from the surface
// normal; Table 1 settles it), rho_h the isentropic exit-plane density.
// A, B, C: Table 2 per set at (theta, alpha). Held below
// rowbury_2001_k_floor (Cd -> 0 and, for B > 1, no root of the hole's flow
// below it). Beyond rowbury_2001_k_fit_max(set) = max(tier 1: the k_perp the
// set's Table 1 rows provably reach, rowbury_2001_k_data_max; tier 2: the
// largest the authors applied it at) -- 4.7 / 4.3 / 3.8 (suction, tier 2,
// Fig. 9; data 0.58) / 0.13 / 0.054, combaero's choice -- the
// offset C fades (w = 1/(1 + (ln(k/k_fit)/FADE_WIDTH)^2), C1), so Cd -> Cd0 as
// k -> inf. Table 1's angle span is a flag, not a refusal. Cd > 1 can occur
// with the leading-edge set at Cd0 near 1 (not capped).
enum class Rowbury2001Set { CylindricalPressureSurface, CylindricalLeadingEdge,
                            CylindricalSuctionSurface, FanShapedPressureSurface,
                            FanShapedSuctionSurface };
constexpr double ROWBURY_2001_K_PERP_MIN = 0.01, ROWBURY_2001_FADE_WIDTH = 1.0,
                 ROWBURY_2001_UNIQUE_MARGIN = 0.5;
double rowbury_2001_k_fit_max(Rowbury2001Set set);   // max(k_data, k_applied)
double rowbury_2001_k_data_max(Rowbury2001Set set);  // tier 1 alone
struct Rowbury2001Coefficients { double A, B, C; };
Rowbury2001Coefficients rowbury_2001_coefficients(Rowbury2001Set set, double theta_deg,
                                                  double alpha_deg);  // throws A, B <= 0
double rowbury_2001_k_floor(Rowbury2001Set set, double theta_deg, double alpha_deg);
bool rowbury_2001_in_range(Rowbury2001Set set, double theta_deg, double alpha_deg);
struct Rowbury2001Cd { double Cd, dCd_dCd0, dCd_dk, delta, k_used; bool held, faded; };
Rowbury2001Cd rowbury_2001_exit_crossflow_cd(double Cd0, double k_perp, Rowbury2001Set set,
                                             double theta_deg, double alpha_deg);
struct Rowbury2001KPerp { double k, dk_dG, dk_drho_h, dk_dq_m; };
Rowbury2001KPerp rowbury_2001_k_perp(double G_h, double rho_h, double rho_m_u_m2,
                                     double theta_deg, double alpha_deg);
// Isentropic expansion (Pt, Tt) -> P, variable cp; analytic partials,
// drho_dY over mass fractions.
struct Rowbury2001ExitDensity { double rho, T, drho_dPt, drho_dTt, drho_dP;
                                std::vector<double> drho_dY; };
Rowbury2001ExitDensity rowbury_2001_exit_density(double Pt, double Tt, double P,
                                                 const std::vector<double>& X);

// The hole-Cd chain in one call (hole_cd.h): the composition point of the
// network's hole elements (python/combaero/network/hole_cd.py resolves the
// options, gathers the inputs and maps the partials onto unknown names).
// Terms, each the kernel's above: 1. base Cd0 -- Fixed (clamped into
// [HOLE_CD_FIXED_MIN, 1]), Metering (metering_cd_and_dRe at Re_D) or
// Discharge (the correlation at the hole Re and the inlet U1/Vi, inclination
// inside DischargeHoleGeometry; dCd/dRe 0 below discharge_cd_re_floor unless
// continue_re_floor); 2. the opt-in exit crossflow (Rowbury) on Cd0, off at
// rho_m_u_m2 <= 0; 3. the impingement jet law, which replaces the hole's Cd
// where it acts (its own partials only; Rogers & Hersh capped: the hole's
// hold). Partials analytic, by input.
constexpr double HOLE_CD_FIXED_MIN = 1e-4;
enum class HoleCdBase { Fixed, Metering, Discharge };
enum class HoleCdJetCrossflow { Off, RogersHersh1976, FlorschuetzMetzgerSu1982 };
struct HoleCdModel { HoleCdBase base; double fixed_cd; MeteringCdCorrelation metering;
                     DischargeCdCorrelation discharge; bool continue_re_floor;
                     bool exit_crossflow; Rowbury2001Set exit_set;
                     double exit_theta_deg, exit_alpha_deg;
                     HoleCdJetCrossflow jet; double xn_d, yn_d, z_d; };
struct HoleCdInputs { double Re, U1_over_Vi, Re_D, G_h, rho_h, rho_m_u_m2, G_ideal, G_c;
                      bool suction; };
struct HoleCdResult {
    double Cd, dCd_dRe, dCd_dU1_over_Vi, dCd_dRe_D, dCd_dG_h, dCd_drho_h,
           dCd_drho_m_u_m2, dCd_dG_ideal, dCd_dG_c;
    double cd_base, cd_hole, dcd_base_dRe, dcd_base_dU1_over_Vi, dcd_base_dRe_D,
           inclination_dzeta;  // per term
    bool re_held, exit_active; double exit_ratio;
    Rowbury2001KPerp k_perp; Rowbury2001Cd exit;
    double exit_k_floor, exit_k_tier;  // tier 0 data, 0.5 authors' use, 1 faded
    bool exit_floor_above_fit, exit_in_range, jet_acts, jet_partials, jet_capped;
};
HoleCdResult hole_cd_compose(const HoleCdModel& model, const DischargeHoleGeometry& hole,
                             const OrificeGeometry& orifice, const HoleCdInputs& inputs);

// Merge chamber (#471): the main-inlet face from the impulse balance
//   P_f + m_main^2 R T_f / (P_f A^2) = P + m_out^2/(rho A^2) - J/A
// with Tt and Tt_main STAGNATION temperatures (#578): the face's static T_f
// solves h(Tt_main) = h(T_f) + u_f^2/2, the outlet's comes from
// chamber_static_state. Subsonic root; Pt_face = P_f exp((s(Tt_main) -
// s(T_f))/R). Sensitivities analytic (the dT/dT_main field suffixes are
// with respect to Tt/Tt_main). `choked` flags a face that cannot carry the
// impulse subsonically.
struct MergeFaceState {
    double P_face, Pt_face, M_face; bool choked;
    double dPf_dm_main, dPf_dm_out, dPf_dJ, dPf_dP, dPf_dT, dPf_dT_main;
    double dPtf_dm_main, dPtf_dm_out, dPtf_dJ, dPtf_dP, dPtf_dT, dPtf_dT_main;
    // Mass-fraction gradients of the chamber gas and the main stream (#743).
    std::vector<double> dPf_dY, dPf_dY_main, dPtf_dY, dPtf_dY_main;
};
// incompressible (#588): the low-Mach limit, densities at the chamber P and
// the stagnation temperatures; linear in P_face, nothing to choke.
MergeFaceState chamber_merge_face_state(
    double m_main, double Tt_main, const std::vector<double>& X_main,
    double m_out, double P, double Tt, const std::vector<double>& X,
    double side_momentum, double area, bool incompressible = false);

// Streamwise jet impulse per unit cos(theta), J = m w: isentropic velocity to
// P (variable cp), or choked: w = u* + (P* - P)/(rho* u*), the sonic point
// on the isentrope at M = 1. Partials analytic (#816). C1 across choking;
// Cd does not enter.
struct JetImpulse { double J, dJ_dm, dJ_dPt, dJ_dTt, dJ_dP; bool choked; };
JetImpulse jet_impulse(double m, double Pt, double Tt, double P,
                       const std::vector<double>& X);
```

```cpp
// Every ejector function throws std::invalid_argument for non-finite input and
// out-of-domain parameters (gamma > 1, r_gas > 0, efficiencies in (0, 1],
// recovery in (0, 1] (above 1: a gain beyond lossless mixing), 1 < A_p1/A_t
// <= A_3/A_t, positive areas, A_exit >= A_throat,
// losses and eps_frac >= 0, blend windows lo < hi); the standalone closures
// also need a finite, positive supply state. The two whole-element residual
// functions check their Newton iterates for finiteness only (#769).
//
// Flag and clamp (#945): the two bisections return their bracket end for a
// target outside the bracket and say so. The nozzle-exit Mach (bracket
// M in [1, 50]) by ejector_nozzle_exit_mach and the element result's
// nozzle_exit_mach_out_of_range; the C-D nozzle exit-static inverse
// (P_py in [P*, Pt (1 - 1e-12)]) by EjectorCDExitStaticJacobian::out_of_range,
// with the clamped function's partials (0, k, 0). -1 below, +1 above, 0 in.
struct EjectorNozzleExitMach { double mach; int out_of_range; };
EjectorNozzleExitMach ejector_nozzle_exit_mach(double area_ratio, double gamma);
//
// Ejector operating anchors (P_c*, P_b0): the end of the critical plateau and
// the far end of the subcritical droop. Two closure families, selected as a
// pair by ejector_element_residuals_and_jacobian's anchor_closure argument:
enum class EjectorAnchorClosure : int { kracik_dvorak = 0, huang1999 = 1 };
// huang1999: P_c* by Huang 1999's printed chain (phi_m mixing, Eq. 19; normal
// shock; isentropic diffuser), C2 through a sonic mixed flow; and P_b0 by the
// sudden expansion over A_3 (zero secondary flow, closed-form subsonic exit).
// Both (value, d/d{p_g, t_g, p_e, t_e}); value.p_c_pa carries the anchor.
double ejector_huang_mixing_coefficient(double area_ratio_mix);
EjectorCriticalPressureJacobian ejector_huang_critical_back_pressure_and_jacobian(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom,
    double gamma, double r_gas, double recovery_efficiency);
EjectorCriticalPressureJacobian ejector_sudden_expansion_dead_head_and_jacobian(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom,
    double gamma, double r_gas, double recovery_efficiency);
// The whole element (4 rows, Jacobian over its nine unknowns) takes the pair
// as its last argument, default huang1999.
EjectorElementResidualJacobian ejector_element_residuals_and_jacobian(
    /* ... flows, port states, geometry, gamma, r_gas, efficiencies, blend
       windows ... */, EjectorAnchorClosure anchor_closure = EjectorAnchorClosure::huang1999);
```

```cpp
// Incompressible jet pump (#588), the ejector's low-Mach limit: Sanger 1968
// (NASA TN D-4445 App. B) with each stream at its own specific volume v = 1/rho.
// Rows R0..R3 (nozzle, suction, mass, outlet) and their Jacobian over
// (mp, ms, mo, Pt_p, v_p, Pt_s, v_s, Pt_o, p3); p4, q4 at the mixing-duct exit.
EjectorJetPumpIncompressibleResult ejector_jetpump_incompressible_residuals_and_jacobian(
    double mp, double ms, double mo, double pt_p, double v_p, double pt_s,
    double v_s, double pt_o, double p3, double area_nozzle_exit,
    double area_secondary, double area_mixing, double k_nozzle,
    double k_suction, double k_throat, double k_diffuser, double recovery);
```

### Wall Coupling
```cpp
WallCouplingResult wall_coupling_and_jacobian(
    double h_a, double T_aw_a, double h_b, double T_aw_b,
    double t_over_k,   // wall_thickness / wall_conductivity [m²·K/W]
    double A = 1.0     // contact area [m²]
);
// Multi-layer overload: t_over_k_layers, A, R_fouling. h is used as given
// at or above WALL_HTC_KNEE = 1e-2 W/(m^2 K); below it the C1, strictly
// positive floor knee^2/(2 knee - h), so a correlation driven to h <= 0 at a
// solver iterate still gives a finite wall. dQ/dh is the slope of the
// floored value.

// What the network uses (#627): the heat two streams exchange, the
// epsilon-NTU integral along the wall, Q = eps(NTU, C_r) C_min (T_aw_a -
// T_aw_b) at the inlets. The linear form above is its NTU -> 0 limit.
// C = |m| cp [W/K]; +inf for a node side (a reservoir). counterflow follows
// from the flow signs and ThermalWall.orientation. Finite at C = 0
// (Q -> C dT), C = inf and C_a = C_b; adds dQ_dC_a, dQ_dC_b, UA and
// effectiveness to the result.
WallCouplingResult wall_exchange_and_jacobian(
    double h_a, double T_aw_a, double C_a,
    double h_b, double T_aw_b, double C_b,
    const std::vector<double>& t_over_k_layers, double A,
    double R_fouling, bool counterflow);
// Both forms also give dQ_dR_wall = dQ/d(t/k) of any one layer (#932).

// A k(T) wall (#932): each layer's t/k at its mean temperature in
// wall_temperature_profile, iterated to a fixed point by the caller. The total
// derivatives of R_wall = sum t/k through that fixed point, by implicit
// differentiation, dr_dTm = d(t/k)_i/dTm_i (0 for a constant k):
std::array<double, 4> wall_resistance_fixed_point_derivatives(
    double T_hot, double T_cold, double h_hot, double h_cold,
    const std::vector<double>& t_over_k, double R_fouling,
    const std::vector<double>& dr_dTm);  // dR_wall/d(T_hot, T_cold, h_hot, h_cold)

// materials.h: k and its analytic dk/dT on every path of
// material_conductivity (data, #882 envelope fade, linear continuation).
ConductivityAndSlope material_conductivity_and_slope(
    const std::string& name, double T, bool extrapolate);  // {k, dk_dT}
```

### Slot Film-Side h (`slot_film_htc.h`, #962)

The heat transfer coefficient under a 2-D slot film, on the adiabatic-wall
temperature. Default: Seban (1960) Fig. 12 as DRAWN, `St_s = 0.135 Re_s^-0.3
(x/s)^-0.3` (curve B; Eq. 13 prints 0.10) and `0.40 Re_s^-0.3 (x/s)^-0.6`
(curve C; printed 0.41), C1-blended over +-0.1 decade at their x/s-37
crossing. Named option: Ballal & Lefebvre (1973) Eq. 6 = Lefebvre & Ballal
Eq. 8.44, `Nu_x = 0.10 Re_s^0.8 (x/s)^0.44`. The constants live in
`namespace seban1960` / `namespace ballal_lefebvre1973`, the validity box
included. The header comment gives the derivation and the misprint
evidence.

```cpp
#include "slot_film_htc.h"

// Span means over [x_start, x_end] / s, the in-box part (x/s <= 300 / 50)
SlotFilmStanton slot_film_stanton_seban_1960(double x_start_over_s,
    double x_end_over_s, double Re_s, double C_near = seban1960::C_near,
    double C_far = seban1960::C_far);   // St, dSt_dRe, fraction_in_box, branch, flags
SlotFilmNusselt slot_film_nusselt_ballal_lefebvre_1973(double x_start_over_s,
    double x_end_over_s, double Re_s);  // Nu_s = h s / k_s, dNu_dRe, ...

// A wall side's h: h_duct + w(M) f_x (h_film - h_duct), slot-fluid properties
// at (T_s, P_s, X_s), G_s = m_s / A_slot, M = G_s / G_gas. Hands back to the
// duct below M_handback (1.0 / 1.3, C1 over SLOT_FILM_HTC_HANDBACK_BAND = 0.1)
// and past x/s max.
SlotFilmHtc slot_film_htc_seban_1960(double x_start, double x_end, double s,
    double G_s, double M, double h_duct, double T_s, double P_s,
    const std::vector<double>& X_s);
SlotFilmHtc slot_film_htc_ballal_lefebvre_1973(/* same */);
// SlotFilmHtc: h, dh_dG_s, dh_dT_s, dh_dM, dh_dh_duct (analytic), h_film,
// Re_s, fraction_in_box, weight_M, branch, x_in_box, Re_in_box, M_in_box
```

### Liner Radiation (`liner_radiation.h`, #793)

Lefebvre & Ballal (2010) Ch. 8: gray, lumped liner radiation. The header's
comment block holds the equations, the book's errata (Eq. 8.9 prints
"2 - exp") and the alternatives that are not shipped (L Eqs. 8.12-8.14;
beam length 3.6 V/A). Each kernel returns its value with analytic partials.

```cpp
// Constants: STEFAN_BOLTZMANN (CODATA 2018), REEVES_EMISSIVITY_COEFF = 290
// [kPa^-1 m^-0.5 K^1.5], LUMINOSITY_HYDROGEN_COEFF = 336,
// BEAM_LENGTH_COEFF_LEFEBVRE = 3.4, BEAM_LENGTH_COEFF_HOWELL = 3.6.

// Eq. 8.11: eps_g = 1 - exp[-290 (P/1000) L (q l_b)^0.5 Tg^-1.5], P in Pa.
GasEmissivityResult gas_emissivity_luminous(double T_g, double P, double q,
                                            double l_b, double L = 1.0);
// {eps, deps_dT, deps_dP, deps_dq, deps_dl_b, deps_dL}

std::tuple<double, double> luminosity_factor_hydrogen(double H_mass_percent);
// Eq. 8.15: (L = 336/H^2, dL/dH)

// Eq. 8.7: R1 = 0.5 sigma (1 + eps_w) eps_g Tg^1.5 (Tg^2.5 - Tw^2.5)
LinerGasRadiationResult liner_gas_radiation(double T_g, double T_w,
                                            double eps_g, double eps_w);
// {R1, dR1_dT_g (fixed eps_g), dR1_dT_w, dR1_deps_g, dR1_deps_w}

// Eq. 8.17 per unit liner area: (R2, dR2/dT_w, dR2/dT_c)
std::tuple<double, double, double> liner_casing_radiation(
    double T_w, double T_c, double eps_w, double eps_c, double area_ratio);
double liner_casing_exchange_factor(double eps_w, double eps_c, double area_ratio);

double mean_beam_length(double volume, double area,
                        double coeff = BEAM_LENGTH_COEFF_LEFEBVRE);  // Eq. 8.10

// q = m(C + H atoms) / m(rest), and dq/dY: the same for reactants and products.
std::tuple<double, std::vector<double>> elemental_fuel_air_ratio(
    const std::vector<double>& Y);

// A wall's local balance with radiation on a face (kind 0 none, 1 gas
// Eq. 8.7 with k = 0.5 sigma (1 + eps_w) eps_g, 2 surface Eq. 8.17 with k the
// exchange factor): Newton on the two surface temperatures, folded into
// effective sides H = h + h_r, T* = (h T_aw + h_r T_src) / H, where h_r is
// the radiation written exactly as a conductance (a divided difference,
// finite at T_s = T_src). Total derivatives by implicit differentiation
// (DualN) in [T_aw_a, h_a, T_aw_b, h_b, R_wall, T_src_a, k_a, T_src_b, k_b].
WallRadiationSides wall_radiation_sides(
    double T_aw_a, double h_a, double T_aw_b, double h_b, double R_wall,
    int kind_a, double T_src_a, double k_a,
    int kind_b, double T_src_b, double k_b);

// The face of a wall whose other side is adiabatic (#970): its own balance
// h (T_aw - T_s) + h_r(T_s) (T_src - T_s) = 0, h_r as above. Kind 0 (or
// k = 0): T_s = T_aw exactly. The left side falls strictly with T_s, so the
// root is unique, between T_aw and T_src (bracketed Newton). Returns T_s,
// h_r, q_rad = h_r (T_src - T_s) [W/m^2] and dT_s/d[T_aw, h, T_src, k]
// (implicit differentiation, DualN).
WallAdiabaticFace wall_adiabatic_face(
    double T_aw, double h, int kind, double T_src, double k);
```

### Pressure-Dependent WSGG Emissivity (`wsgg_bordbar2021.h`, #1071)

Bordbar et al. (2021), Int. J. Heat Mass Transfer 173, 121207 (CC BY 4.0):
non-luminous H2O-CO2, 300-3000 K, Mr = x_H2O/x_CO2 0.1-4, ten pressure sets
0.1-80 atm. The header's comment block holds Eqs. 7-11, the interpolation
between the sets and the clamps; the constants are Appendix A's, in
`wsgg_bordbar2021_data.h`.

```cpp
// Constants: WSGG_BORDBAR2021_T_REF = 1500 K, _T_MIN/_T_MAX = 300/3000 K,
// _MR_MIN/_MR_MAX = 0.1/4, WSGG_PA_PER_ATM = 101325.

// One tabulated set's gray gases (Eqs. 8, 10, 11): a[4] [-], kappa[4] [1/(atm m)].
// Unclamped: throws outside the fitted range or off a tabulated pressure.
WsggGrayGases wsgg_bordbar2021_gray_gases(double p_set_atm, double T, double Mr);

// eps across the sets (linear in ln P, each set at the actual path), clamped.
WsggEmissivityResult gas_emissivity_bordbar2021(double T_g, double P,
                                                double x_h2o, double x_co2,
                                                double l_b);
// {eps, deps_dT, deps_dP, deps_dx_h2o, deps_dx_co2, deps_dl_b}

// The same from mass fractions: {eps, deps_dT, deps_dP, deps_dl_b, deps_dY}.
WsggMixtureEmissivityResult gas_emissivity_bordbar2021_mixture(
    double T_g, double P, const std::vector<double>& Y, double l_b);
```

### Data Structures
```cpp
struct ChannelResult {
    // Primary outputs
    double h;              // Heat transfer coefficient [W/(m²·K)]
    double Nu;             // Nusselt number [-]
    double Re;             // Reynolds number [-]
    double Pr;             // Prandtl number [-]
    double f;              // Friction factor [-]
    double dP;             // Pressure drop [Pa]
    double M;              // Mach number [-]
    double T_aw;           // Adiabatic wall temperature [K]
    double q;              // Heat flux [W/m²] (nan if T_hot not supplied)

    // Jacobians
    double dh_dmdot;       // dh/dmdot [W/(m²·K·kg/s)]
    double dh_dT;          // dh/dT [W/(m²·K²)]
    double ddP_dmdot;      // d(dP)/dmdot [Pa·s/kg]
    double ddP_dvelocity;  // d(dP)/d(velocity) [Pa·s/m] - chain from this, not
                           // ddP_dmdot, which uses the correlation's own flow area
    double ddP_dT;         // d(dP)/dT [Pa/K]
    double dT_aw_dmdot;    // dT_aw/dmdot [K·s/kg]
    double dT_aw_dT;       // dT_aw/dT [-] (approx 1 at low Mach)
    double dT_aw_dP;       // dT_aw/dP [K/Pa] at fixed mdot, T (#496)
    double dq_dmdot;       // dq/dmdot [W·s/(m²·kg)]
    double dq_dT;          // dq/dT [W/(m²·K)]
    double dq_dT_hot;     // dq/dT_hot [W/(m²·K)]
    // Mass-fraction gradients at fixed mdot, T, P (#717); channel_smooth
    // fills them, other correlations leave them empty (read as zero).
    std::vector<double> dh_dY, dT_aw_dY;
};

struct WallCouplingResult {
    double Q;              // Heat transfer rate [W]
    double dQ_dh_a;        // ∂Q/∂h_a [m²·K/W]
    double dQ_dh_b;        // ∂Q/∂h_b [m²·K/W]
    double dQ_dT_aw_a;     // ∂Q/∂T_aw_a [W/K]
    double dQ_dT_aw_b;     // ∂Q/∂T_aw_b [W/K]
};
```

---

## Geometry Utilities (geometry.h)

```cpp
double hydraulic_diameter(double A, double P_wetted);         // 4 A / P
double hydraulic_diameter_rect(double a, double b);            // 2ab / (a + b)
double hydraulic_diameter_annulus(double D_outer, double D_inner);  // Do - Di
double circular_area(double D);
double annular_area(double D_outer, double D_inner);           // pi (Do^2 - Di^2) / 4
double cylinder_volume(double D, double L);
// Every function and struct method throws std::invalid_argument unless
// D_outer > D_inner >= 0 and every size it uses is finite and positive (NaN
// and inf included, #769); the struct methods check the members they read.

struct Tube { double L, D; area(); volume(); perimeter(); };
struct Annulus {
    double L, D_inner, D_outer;   // diameters [m]
    int n_azimuthal_max = 5;      // acoustic solvers
    D_mean(); gap();              // gap = (Do - Di) / 2, half the hydraulic diameter
    area(); volume(); circumference();  // circumference = pi D_mean
    radius_inner(); radius_outer();
};
struct CanAnnularFlowGeometry {   // primary zone + transition + annulus
    double L, D_inner, D_outer, L_primary, D_primary, L_transition;
    // Annulus helpers plus area_primary, volume_primary, volume_transition,
    // volume_total, length_total
};

double residence_time(double V, double Q);
double residence_time_mdot(double V, double mdot, double rho);
// overloads for Tube, Annulus, CanAnnularFlowGeometry; space_velocity(Q, V)
```

`channel_roughness` lives in friction.h.

---

## Orifice Flow (orifice.h)

### Geometry and State
```cpp
enum class OrificeType {
    SharpThinPlate, ThickPlate, RoundedEntry, Conical, QuarterCircle, UserDefined
};

// Pressure taps (ISO 5167-2:2003 5.3.2.1). A measurement choice: it fixes
// what C is referenced to, not what the plate does to the flow.
enum class OrificeTap { Corner, D_D2, Flange };

struct OrificeGeometry {
    double d, D;
    OrificeTap tap = OrificeTap::Flange;
    double t, r, bevel;

    double beta() const;          // Diameter ratio d/D [-]
    double area() const;          // Orifice area [m²]
    double t_over_d() const;      // Thickness ratio t/d [-]
    double r_over_d() const;      // Radius ratio r/d [-]
    bool is_valid() const;
};

struct OrificeState {
    double Re_D, dP, rho, mu;

    double Re_d(double beta) const;  // Orifice Reynolds number (based on d)
};
```

### Cd Correlations
```cpp
// NORMED measurement orifice: a standardised plate in a pipe, Cd referenced
// to the TAPPING differential, every member a function of beta = d/D.
enum class MeteringCdCorrelation {
    ReaderHarrisGallagher,          // ISO 5167-2:2003
    Constant, UserFunction
};
// Stolz (ISO 5167:1980) and Cd_Stolz were removed in 0.9.0 (#876): ISO
// replaced it with Reader-Harris/Gallagher (ISO 5167-1:1991/Amd 1:1998,
// ISO 5167-2:2003). Use ReaderHarrisGallagher / Cd_ReaderHarrisGallagher.

// Reader-Harris/Gallagher at geom.tap, with the small-pipe term below 71.12 mm.
double Cd_sharp_thin_plate(const OrificeGeometry& geom, const OrificeState& state);
// (Cd, dCd/dRe_D), analytic; held below ISO's lower Re_D bound (below).
// The network orifice's flow column uses it (dRe_D/dm = Re_D/|m|).
std::tuple<double, double> Cd_sharp_thin_plate_and_dRe(const OrificeGeometry& geom,
                                                       const OrificeState& state);
// (Cd, dCd/dRe_D) of ReaderHarrisGallagher at geom.tap; the network
// element's only metering entry point. Constant/UserFunction throw. RHG below
// Re_b = metering_cd_re_min is HELD (#873, 0.8.0), entered C1 in ln Re_D over
// orifice::metering_hold::blend_ln_re = 0.1 just below Re_b (the slope ramps
// linearly to 0); at Re_D >= Re_b nothing changed. The form and its
// derivative are in include/orifice.h. A 1.5 cap stays for flange taps below
// D ~32 mm.
std::tuple<double, double> metering_cd_and_dRe(MeteringCdCorrelation id,
                                               const OrificeGeometry& geom, double Re_D);
// RHG's lower Re_D bound, ISO 5167-2:2003 5.3.1 (5000, 16000 b^2 above
// beta 0.56, flange max(5000, 170 b^2 D[mm])). Throws std::invalid_argument
// for every other correlation (Constant/UserFunction carry no formula).
double metering_cd_re_min(MeteringCdCorrelation id, const OrificeGeometry& geom);
// ISO 5167-2:2003 5.3.1 range (Reader-Harris/Gallagher); Constant and
// UserFunction return true.
bool metering_cd_in_range(MeteringCdCorrelation id, const OrificeGeometry& geom,
                          double Re_D);
// ISO 5167-2:2003 5.4.1: permanent loss / tapping differential of a plate
// between equal pipes, k = [s - C b^2] / [s + C b^2], s = sqrt(1 - b^4 (1 - C^2)),
// and dk/dC. Derived for vena-contracta (D-D/2) taps. Throws outside
// 0 < beta < 1 or C <= 0.
std::tuple<double, double> permanent_loss_ratio_and_dC(double beta, double C);
// ISO 5167-1 velocity-of-approach factor E = 1 / sqrt(1 - beta^4); 1 outside
// 0 < beta < 1 (no approach pipe, or bore >= pipe).
double velocity_of_approach_factor(double beta);

// Jobson (1955), Proc. IMechE 169, 767: compressible contraction of a
// sharp-edged orifice from its incompressible Ci (#854). Force defect with the
// approach-velocity extension, f = 1/Ci - 1/(2 Ci^2) - Cv b^2 + Cv^2 b^4 / 2
// (Eq. 5 at beta = 0). Constants in orifice::jobson1955 (Ci 0.5-0.7, beta
// <= 0.75: a diagnostic range) and orifice::iso5167_expansibility (Eq. 5).
double jobson_force_defect(double Ci, double beta = 0.0, double Cv = 1.0);
// Ci of an orifice whose incompressible flow is Cd E A sqrt(2 rho1 dp):
// Cv Ci / sqrt(1 - (Cv Ci)^2 beta^4) = Cd E.
double jobson_ci_from_cd(double Cd, double beta = 0.0, double Cv = 1.0);
// The paper's C(r), reservoir-fed: Eq. (11) for r >= r_c, Eq. (15) below.
double jobson_contraction_coefficient(double r, double Ci, double n);
// mdot / mdot_inc at the same upstream static state, r = p2/p1, the beta
// extension solved to consistency (bracketed Newton, unique root, any r in
// (0, 1], #957): comparable with ISO 5167-2 Eq. (5).
double jobson_expansibility(double r, double beta, double Ci, double n);

// Fixed Cd from an explicit value. make_correlation(MeteringCdCorrelation::Constant)
// can only return the default, so this is how a chosen constant is pinned.
std::unique_ptr<OrificeCorrelationBase> make_constant_correlation(double Cd);
```

### Discharge holes: a separate selector

A hole in a wall that dumps coolant out of its circuit is not a metering
orifice: there is no pipe to form `beta` with, and `Cd` depends on `L/d`,
`r/d` and the approach crossflow instead. It gets its own selector and its own
geometry/state pair.

```cpp
enum class DischargeCdCorrelation {
    McGreehanSchotsch1988,   // cooling hole with inlet crossflow
    Idelchik1966Sharp,       // diagram 4-17, Re 25..1e6
    Idelchik1966Thick,       // diagram 4-18a, deep hole
    Idelchik1966Beveled,     // diagram 4-18b
    Idelchik1966Rounded,     // diagram 4-18c
    Lichtarowicz1965,        // long orifice, l/d 2-10, Re >= 10 (Eq. 12 fitted to 2e4)
    Constant
};

struct DischargeHoleGeometry {
    double d, L, r, bevel;   // [m]; no pipe diameter, by design
    double inclination_deg = 90.0;  // hole axis to wall [deg]; < 90 opts in
                                    // to Idelchik diagram 3-2 (Thick only)
    double L_over_d() const;
    double r_over_d() const;
    double bevel_over_d() const;
    double area() const;
    bool is_inclined() const;  // inclination_deg < 90
    bool is_valid() const;
};

// Idelchik diagram 3-2's inclined-inlet increment over the normal sharp
// inlet: 0.3 cos(delta) + 0.2 cos^2(delta), exactly 0 at 90 deg. Geometry
// only; throws outside 0 < delta <= 90.
double discharge_inclination_dzeta(double inclination_deg);

struct DischargeHoleState {
    double Re;            // based on the hole diameter
    double U1_over_Vi;    // INLET tangential velocity ratio; 0 for a plenum
};

class DischargeCorrelationBase {
public:
    virtual double Cd(const DischargeHoleGeometry&, const DischargeHoleState&) const = 0;
    // (Cd, dCd/dRe, dCd/d(U1_over_Vi)), analytic
    virtual std::tuple<double, double, double> Cd_and_derivatives(
        const DischargeHoleGeometry&, const DischargeHoleState&) const = 0;
    virtual std::string name() const = 0;
};

std::unique_ptr<DischargeCorrelationBase> make_discharge_correlation(DischargeCdCorrelation id);
std::unique_ptr<DischargeCorrelationBase> make_constant_discharge_correlation(double Cd);

// Which input lies outside the source's own range, by variable: INFO to
// report beside the result, never a selector and never a warning (Cd is
// still returned, held or continued). re / l_over_d / edge (bevel or r/d
// table; an inclined hole on a correlation that refuses it). Thick beyond l/d
// 4, McGreehan-Schotsch beyond l/d 10 (mcgreehan_schotsch::l_over_d_max) or
// below Re 1e4, Lichtarowicz below Re 10, above 1e5 (the end of its data;
// 2e4-1e5 is Eq. 7, in range) or outside l/d 2-10.
struct DischargeCdEnvelope {
    bool re = false, l_over_d = false, edge = false;
    bool in_range() const;
};
DischargeCdEnvelope discharge_cd_envelope(DischargeCdCorrelation id,
                                          const DischargeHoleGeometry&, const DischargeHoleState&);
// discharge_cd_envelope(...).in_range()
bool discharge_cd_in_range(DischargeCdCorrelation id,
                           const DischargeHoleGeometry&, const DischargeHoleState&);

// Re below which the Cd VALUE is held (Re clamped up to it): its true dCd/dRe
// is 0 there. Cd_and_derivatives continues the slope below it for
// McGreehan-Schotsch and Idelchik (a Newton aid); a sensitivity built on the
// held value takes the zero. 0 for Constant.
double discharge_cd_re_floor(DischargeCdCorrelation id);

// The hole as an orifice in a straight pipe (Idelchik 4-13..4-16): permanent
// loss on the PIPE velocity head, and the factor R with C_pipe = Cd_wall R
// that gives it through ISO 5167-2's recovery relation (R -> 1 as beta -> 0).
// Geometry only; R = 1 outside 0 < beta < 1 and for Constant.
double discharge_pipe_zeta(DischargeCdCorrelation id, const DischargeHoleGeometry&, double beta);
double discharge_pipe_cd_ratio(DischargeCdCorrelation id, const DischargeHoleGeometry&,
                               double beta);
```

For the Idelchik members `zeta` is referenced to the hole velocity and carries
the full permanent loss, so `Cd = 1/sqrt(zeta)` exactly.

**Inclined holes, opt-in.** By default a hole's angle enters only through
its length `L = t/sin(angle)`. `inclination_deg < 90` adds diagram 3-2's
inclined-inlet loss to `Idelchik1966Thick`'s inlet coefficient zeta' (Eq.
4-20 builds the wall hole as zeta' + 1 + tau + zeta_fr with zeta' from
diagrams 3-1..3-3; 4-18a's table has the normal sharp inlet's 0.5 in it), so
zeta' gains `discharge_inclination_dzeta`; tau keeps its normal-hole value
(our composition). It inherits zeta''s low-Re scaling (4-18a item 2), so
`dCd/dRe` stays analytic. Every other correlation, and `discharge_pipe_zeta`
on them, throws for an inclined hole; `discharge_cd_in_range` also needs
diagram 3-2's Re > 1e4 and `l/d >= orifice::idelchik::inclination_l_over_d_min`
(2, where 4-18a's zeta' has plateaued: the inlet jet has reattached). The
in-pipe loss takes the same inlet coefficient. Constants:
`orifice::idelchik::inclination_*`. See
`validation/cooling/extractions/idelchik_1966_wall_orifice.md` for the tables,
the cross-source agreement with McGreehan-Schotsch, and why the in-a-pipe
diagrams are a different conversion.

### Plenum-to-plenum holes: McGreehan and Schotsch (1988)

A different configuration from the ISO 5167 family above -- a cooling transfer
hole or jet-plate hole, not a metering orifice in a pipe run. ASME J.
Turbomachinery 110(2), 213-217, Eqs. (8)-(17). Each stage is exposed
separately so a check binds to one equation rather than the whole chain.

```cpp
namespace orifice::mcgreehan_schotsch {
double reynolds_baseline(double Re);                 // Eq. (8)
double nozzle_baseline(double Re);                   // Eq. (9)
double corner_factor(double r_over_d);               // Eq. (12), f
double length_factor(double L_over_d);               // Eq. (14), g
double cd_with_corner(double Re, double r_over_d);   // Eq. (11)
double cd_with_corner_and_length(double Re, double r_over_d,
                                 double L_over_d);   // Eqs. (13), (15), (16)
double cd(double Re, double r_over_d, double L_over_d, double U1_over_Vi,
          double eps = rv_smooth_eps);               // Eq. (17)
double cd_with_crossflow(double cd_base, double U1_over_Vi,
                         double eps = rv_smooth_eps);  // Eq. (17) alone

// Solver-facing form: (Cd, dCd/dRe, dCd/d(U1/Vi)). r/d and L/d are geometry
// and carry no partials. Two numerical treatments, both stated in the header:
// the crossflow input is regularised by rv_smooth_eps (Eq. (17) is unbounded
// in slope at U1/Vi = 0, which is the default), and below re_min the value
// stays floored while the derivative is continued from the floor.
std::tuple<double, double, double> cd_and_derivatives(double Re,
                                                      double r_over_d,
                                                      double L_over_d,
                                                      double U1_over_Vi);

// Adiabatic expansion factor, for the INCOMPRESSIBLE form only.
// regime='compressible' already solves the isentropic nozzle exactly via
// nozzle_flow; Eq. (5) reproduces that to 0.008%, so using both double-counts.
double critical_pressure_ratio(double gamma);
double expansion_orifice(double S, double gamma);                 // Eq. (4)
double expansion_nozzle(double S, double gamma);                  // Eq. (5)
double expansion_blend_weight(double cd, double eps = x_smooth_eps);  // Eq. (7)
double expansion_factor(double cd, double S, double gamma,
                        double eps = x_smooth_eps);               // Eq. (6)
}

// The full chain
double orifice::Cd_McGreehanSchotsch(double Re, double r_over_d,
                                     double L_over_d, double U1_over_Vi);
```

`U1_over_Vi` is the INLET (approach, supply-side) tangential velocity ratio,
0 for a plenum-fed plate. It is not a discharge-side crossflow ratio; see the
header and `validation/cooling/extractions/orifice_discharge_coefficient.md`.
`Cd` is not monotonic in it, and is held constant below `Re = 1e4`, its stated
validity floor -- both deliberate, both from the source.

### Film cooling: Baldauf et al. (2002)

Laterally averaged adiabatic film-cooling effectiveness downstream of one row
of **cylindrical**, streamwise-inclined holes. Valid from the ejection point
rather than only far downstream, and it carries the **adjacent jet
interaction** -- the lateral hole-spacing effect driving jet lift-off -- as a
correlated parameter rather than an excluded case.

```cpp
// eta = (T_G - T_AW)/(T_G - T_C), the same convention as
// adiabatic_wall_temperature(), so the two compose directly.
double film_effectiveness_baldauf_2002(
    double x_over_D, double M, double P, double alpha_deg, double s_over_D,
    double Tu, CorrelationStatus *status = nullptr,
    double b_0_override = std::numeric_limits<double>::quiet_NaN());

// Solver-facing (f, J): (eta, d eta/dM, d eta/dP), analytic via dual numbers.
std::tuple<double, double, double>
film_effectiveness_baldauf_2002_and_derivatives(
    double x_over_D, double M, double P, double alpha_deg, double s_over_D,
    double Tu, CorrelationStatus *status = nullptr,
    double b_0_override = std::numeric_limits<double>::quiet_NaN());
```

`b_0_override` replaces Eq. (31)'s `b_0`; NaN (the default) uses Eq. (31) as
printed. It exists so the paper's Table 4 reading can be *tested* -- the
library offers no second formula because the paper contains none. **The
discrepancy is governed by `x/D`, not by `M`**: `b_0` reaches the model only
through `b_1 = b_0 / (1 + M^-3)` (Eq. 32), the descending-branch gradient, so
it is under 1% at every `M` for `x/D <= 20` and reaches -28% by `x/D = 400`
(worst over the envelope, -75%). At effusion row spacings no effusion dataset
can arbitrate it.

Stated envelope (`baldauf2002::`): `M` 0.2-2.5, `P` 1.2-1.8, `s/D` 2-5,
`alpha` 30-90 deg, `Tu` 0.0035-0.075. The paper's own RMS deviation is 5.5%.

**Reported, not enforced.** Outside the envelope both functions still answer,
finitely and with an unchanged value -- a network solve transits odd states
during Newton iteration, and refusing there would break convergence rather
than protect anyone. They now say so, following
`nusselt_dittus_boelter`'s contract: pass a `CorrelationStatus*` to take the
flag silently, or leave it null to get a warning per out-of-range input
through the global warning handler.

`alpha_deg` is the ejection angle to the **surface** and is converted
internally; the paper's trigonometry is in radians. That is not an
assumption -- seven of its Table 4 coefficients reproduce on radians and fail
on degrees by 10-30%.

**Eq. (31) is implemented as printed and contradicts the paper's own
Table 4** by 36% in `b_0`. Under 5% effect below `M ~ 0.5`, up to 50% at
`M = 2.5`. See
`validation/cooling/extractions/baldauf_2002_film_effectiveness.md`.

### Film cooling: Colban, Thole and Bogard (2011), shaped holes

Laterally averaged adiabatic effectiveness downstream of one row of
**shaped** (fan or laidback fan) holes, Eq. (19):

```cpp
// eta = 1 / (P/t + C1 M^C2 xi^C3),  xi = (4/pi)(x/D)(P/D)/(M AR)
double film_effectiveness_colban_2011(
    double x_over_D, double M, double area_ratio, double P_over_D,
    double t_over_P, double alpha_deg, double density_ratio,
    CorrelationStatus *status = nullptr);

// Solver-facing (f, J): (eta, d eta/dM), closed form.
std::tuple<double, double> film_effectiveness_colban_2011_and_derivatives(
    double x_over_D, double M, double area_ratio, double P_over_D,
    double t_over_P, double alpha_deg, double density_ratio,
    CorrelationStatus *status = nullptr);
```

`D` is the metering diameter, `M` the metering-section blowing ratio,
`area_ratio = A_exit/A_inlet`, `t_over_P` breakout width over pitch. At
`x/D = 0` it returns `t/P`, Eq. (12)'s mixing limit, built into the form.
Coefficients and limits live in `colban2011::` (Table 3, Table 4, and the
fitted spans of Table 1). `alpha_deg` and `density_ratio` are range-checked
only: the fit never varied them (30 deg, DR 1.7-2.0), which is why they are
absent from the equation, not because they do not matter. Same
reported-not-enforced contract as Baldauf. See
`validation/cooling/extractions/colban_2011_shaped_film.md`.

### Film cooling: Marek and Tacina (1975), 2-D slot

Effectiveness downstream of a continuous **slot** (a liner's cooling ring),
Juhasz-Marek mixing with the mixing coefficient from the free-stream
turbulence, NASA TN D-7958 Eqs. (4) and (5):

```cpp
// eta = 1 / (1 + C_m x/(M S)),  C_m = 0.25 Tu (x C_m_multiplier), Tu a fraction
double slot_film_effectiveness_marek_tacina_1975(
    double x_over_S, double M, double Tu, double C_m_multiplier = 1.0,
    CorrelationStatus *status = nullptr);

// Solver-facing (f, J): (eta, d eta/dM), closed form.
std::tuple<double, double> slot_film_effectiveness_marek_tacina_1975_and_derivatives(
    double x_over_S, double M, double Tu, double C_m_multiplier = 1.0,
    CorrelationStatus *status = nullptr);
```

`S` is the slot height, `M = rho_s U_s/(rho_H U_H)`, `Tu` the axial
free-stream turbulence as a FRACTION (0.14 = 14%; the report writes Eq. 5
as 0.25 Tu/100 with Tu in percent, the same line; above 1 throws, #773).
`C_m_multiplier` is the user's
rig-matching factor (Eq. 5's band is 0.4-1.6 of the line); the harness never
scores it. Eq. (4) has no cp ratio, unlike Juhasz-Marek (1971). The run's
span (Tu 0.07-0.35, M 0.41-1.87, x/S up to 33.8) and the authors' +-0.04 eta
measurement band live in `marek_tacina1975::`. See
`validation/cooling/extractions/marek_tacina_1975_d7958.md`.

### Effusion plate internal heat transfer

Andrews 86-GT-225. A coolant hole cools its wall in two places and the
paper's point is that the approach dominates at the Reynolds numbers
effusion runs at.

```cpp
// Mills' entry-length factor, Eqs. (13) and (14). Two curve fits meeting
// at L/D = 2; decays to 1 as the hole lengthens.
double mills_entry_length_factor(double L_over_D);

// Sparrow's hole approach, Eq. (18): 0.881 Re^0.476 Pr^(1/3) X/(pi L)
double effusion_approach_nusselt(double Re, double Pr, double X_over_L);

// Mills' short-hole throat, Eq. (12): 0.023 Re^0.8 Pr^(1/3) R_Nu
double effusion_throat_nusselt(double Re, double Pr, double L_over_D);

// Eq. (19): the two summed.
double effusion_internal_nusselt(double Re, double Pr, double X_over_L,
                                 double L_over_D);
```

`Re` is on the **hole diameter**, and every `Nu` is referenced to the
**hole internal area** `pi D L` -- which is what makes the two summable.
For the coefficient per unit **plate** area that an effusion element needs:

```
h_plate = Nu * k / D * A_h / A,   A_h = pi D L,  A = X^2 - pi D^2 / 4
```

That factor is 3.46 for Andrews' plate C, so it is not optional. The
conversion is left to the caller because it needs only geometry.

### Turbulent flat plate (`cooling_correlations.h`)

The unblown reference of the Stanford full-coverage data: NASA CR-2648's
`St_o Pr^0.4 = 0.0295 Re_x^-0.2` (constants in `cooling::crawford1976`).

```cpp
// Local: (Nu_x, dNu_x/dRe_x), Nu_x = 0.0295 Re_x^0.8 Pr^0.6.
std::tuple<double, double> flat_plate_turbulent_nusselt_and_derivative(
    double Re_x, double Pr, CorrelationStatus* status = nullptr);

// Averaged over x in [x0, x0 + L], on L:
//   Nu_L = (0.0295/0.8) Pr^0.6 Re_L^0.8 [(1 + s)^0.8 - s^0.8],  s = x0/L
struct FlatPlatePanelNusselt { double Nu, dNu_dRe_L, Re_x_start, Re_x_end; bool extrapolated; };
FlatPlatePanelNusselt flat_plate_turbulent_panel_nusselt(double Re_L, double Pr,
                                                         double x_start_over_L);
std::tuple<double, double> flat_plate_turbulent_panel_nusselt_and_derivative(
    double Re_L, double Pr, double x_start_over_L, CorrelationStatus* status = nullptr);
```

Outside the data's Re_x 1.5e5-5e6 the status (or `extrapolated`) says so and
the value is still returned; no warning is printed. `EffusionPlateElement`'s
`gas_baseline="flat_plate"` is the network use.

Scored against Andrews' own Fig. 8 at **-9.9% bias, 11.6% RMSE** (-13.5% /
14.7% while air's viscosity was 11.8% low, #485), labelled
**accuracy**: the correlations are the 1986 paper and the data the 1988
one, and same lab plus different study is cross-source. See
`validation/cooling/extractions/andrews_effusion_internal_h.md`.

### Multi-row film superposition

```cpp
// Sellers, Gao Eq. (1):  eta = 1 - prod_i (1 - eta_i)
double film_superposition_sellers(const std::vector<double>& eta_rows);

// Gao Eq. (7): Sellers with a per-row mainstream temperature correction.
// alpha_between_rows has n-1 entries; all ones reproduces Sellers exactly.
double film_superposition_corrected(const std::vector<double>& eta_rows,
                                    const std::vector<double>& alpha_between_rows);

// (eta, d eta / d eta_i) for each of the two above. Built from explicit
// partial products rather than dividing the total, so a fully effective row
// (eta_j = 1) gives 0, not NaN.
std::pair<double, std::vector<double>> film_superposition_sellers_and_gradient(
    const std::vector<double>& eta_rows);
std::pair<double, std::vector<double>> film_superposition_corrected_and_gradient(
    const std::vector<double>& eta_rows,
    const std::vector<double>& alpha_between_rows);

// Panel-averaged effusion film (#471): Baldauf per row, Sellers-superposed,
// averaged over the rows -- (1/N) sum_n Sellers{eta(j pitch_x/D), j=1..n}.
// The alternative closure Andrews 88-GT-290 refuses offered alone; select it
// explicitly. extrapolated flags Baldauf's envelope.
struct EffusionPanelFilm { double eta; bool extrapolated; };
EffusionPanelFilm effusion_panel_film_effectiveness(int n_rows, double pitch_x_over_D,
                                                    double s_over_D, double M,
                                                    double density_ratio,
                                                    double alpha_deg, double Tu);

// Film a DOWNSTREAM wall sees (#386): n_rows rows Sellers-superposed, averaged
// over [x_start, x_end] measured from the last row (in D); x_start == x_end is
// the point value. Fixed composite Gauss-Legendre in ln(1 + x/D); exact
// derivatives. Read by ThermalWall(film=...).
struct FilmSpanAverage { double eta, deta_dM, deta_dP; bool extrapolated; };
FilmSpanAverage film_span_average_baldauf_2002(int n_rows, double row_pitch_over_D,
                                               double x_start_over_D, double x_end_over_D,
                                               double M, double density_ratio,
                                               double alpha_deg, double s_over_D, double Tu);
FilmSpanAverage film_span_average_colban_2011(int n_rows, double row_pitch_over_D,
                                              double x_start_over_D, double x_end_over_D,
                                              double M, double area_ratio, double P_over_D,
                                              double t_over_P, double alpha_deg,
                                              double density_ratio);  // deta_dP = 0
// One 2-D slot, Marek & Tacina (1975); x from the slot lip in slot heights.
FilmSpanAverage slot_film_span_average_marek_tacina_1975(double x_start_over_S,
                                                         double x_end_over_S, double M,
                                                         double Tu,
                                                         double C_m_multiplier = 1.0);  // deta_dP = 0

// Rows in GROUPS (an effusion liner's panels, #695), upstream first, each at
// its own M, density ratio and coolant temperature, layered row by row:
// T_aw,bar = T_ref (1 - eta) + sum_g weight[g] T_c,g, sum_g weight[g] = eta
// (Sellers when the T_c,g are equal). d_weight_dM[g][h] = d weight[g]/d M_h.
struct FilmSpanAverageGroups {
  double eta;
  std::vector<double> weight, deta_dM, deta_dP;
  std::vector<std::vector<double>> d_weight_dM, d_weight_dP;
  bool extrapolated;
};
FilmSpanAverageGroups film_span_average_baldauf_2002_groups(
    const std::vector<int>& n_rows, double row_pitch_over_D, double x_start_over_D,
    double x_end_over_D, const std::vector<double>& M,
    const std::vector<double>& density_ratio, double alpha_deg, double s_over_D,
    double Tu);

// Gao Eq. (5): the published FORM of alpha. a and b are REQUIRED -- the
// published values 12 and 0.9465 (gao_a_case1/b_case1) do not transfer:
// alpha >= b for every r, so they cannot reach the 0.69-0.85 a tighter
// plate needs.
double mainstream_temperature_correction(double mass_flow_ratio, double a, double b);

double equivalent_slot_width(double hole_area, double pitch);           // Eq. (9)
double equivalent_blowing_ratio(double M0, double A0, double Ae);        // Eq. (10)
```

Both superposition forms also have `..._and_gradient` variants returning
`(eta, d eta/d eta_i)`, built from partial products so a fully effective row
does not divide by zero.

**Sellers overestimates, and worsens with row count** -- Gao's measurement,
and the reason effusion cannot simply reuse a few-row film model. `alpha` is
the correction: an energy balance on mainstream entrained at each injection
(Gao Eqs. 3-4), bounded in [0,1], equal to 1 for uncorrected Sellers. See
`validation/cooling/extractions/film_superposition.md`.

---

## Acoustics (acoustics.h)

```cpp
std::vector<AcousticMode> tube_axial_modes(const Tube& tube, double c,
                                           BoundaryCondition upstream,
                                           BoundaryCondition downstream,
                                           int n_max = 5);

double helmholtz_frequency(double V, double A_neck, double L_neck, double c,
                           double end_correction = 0.85);

double acoustic_impedance(double rho, double c);
double sound_pressure_level(double p_rms, double p_ref = 20e-6);
```

### Can-annular and annular-duct modes (acoustics.h, can_annular_solvers.h)

```cpp
struct CanAnnularGeometry { int n_cans; double length_can, area_can, radius_plenum, area_plenum; };
struct BlochMode { int m_azimuthal; double frequency; int n_cans; std::string symmetry_type() const; };
CanAnnularGeometry to_acoustic_geometry(const CanAnnularFlowGeometry&, int n_cans,
                                        double L_can, double D_can);

// Bloch-Floquet: per m = 0..N/2, roots of the junction volume-flow balance
// (can + two plenum sectors). Default method: magnitude minimization + Muller.
std::vector<BlochMode> can_annular_eigenmodes(const CanAnnularGeometry&, double c_can,
    double c_plenum, double rho_can, double rho_plenum, double f_max = 2000,
    BoundaryCondition bc_can_top = BoundaryCondition::Closed);
enum class CanAnnularSolverMethod { MagnitudeMinimization, ArgumentPrinciple };
std::vector<BlochMode> can_annular_eigenmodes_with_method(..., CanAnnularSolverMethod);  // C++ only
std::vector<BlochMode> solve_magnitude_minimization(...);   // C++ only
std::vector<BlochMode> solve_argument_principle(...);       // C++ only

// Uniform thin annular duct, rigid at x = 0, bc_ends at x = L (Open = closed-open).
// Rigid ends include the axially uniform n_axial = 0 mode for m >= 1. f_max exclusive.
struct AnnularMode { int m_azimuthal, n_axial; double frequency; };
std::vector<AnnularMode> annular_duct_eigenmodes(const Annulus&, double c, double rho,
    double f_max, BoundaryCondition bc_ends = BoundaryCondition::Closed);
std::vector<AnnularMode> annular_duct_modes_analytical(const Annulus&, double c,
    double f_max, BoundaryCondition bc_ends = BoundaryCondition::Closed);
using AnnularDuctGeometry [[deprecated]] = Annulus;  // since v0.2.0, no removal version set
```

---

## Humid Air (humidair.h)

```cpp
double humidity_ratio(double T, double P, double RH);
double dewpoint(double T, double P, double RH);  // RH in (0, 1]; RH = 1 -> T; RH > 1, RH <= 0 or non-finite throw invalid_argument
double humid_air_density(double T, double P, double RH);

struct HumidAir {
    State state;      // public member
    double rh = 0.0;  // RH used at the last set_*; not updated if state.T/P change
    HumidAir& set_TP_RH(double T, double P, double RH);
    HumidAir& set_TP_dewpoint(double T, double P, double Tdp);
    HumidAir& set_TP_omega(double T, double P, double w);
    double RH() const;
    double dewpoint() const;
    double wet_bulb() const;
    double humidity_ratio() const;
    double h_mass() const;   // ASHRAE psychrometric, J/kg dry air
};
```

Saturation pressure is Huang (2018), over ice at T <= 273.15 K (RH below
0 C is relative to ice), accepted on 163.15-383.15 K. Mixing is ideal-gas
with no enhancement factor: humidity ratio is ~0.5% low against ASHRAE
RP-1485 at 1 atm, ~2% at 5 bar and ~13% at 30 bar on cold saturated air;
density and `humid_air_enthalpy` are ideal-gas as well.
`humid_air_enthalpy_nasa9` is on the NASA-9 datum (elements at 298.15 K),
not the psychrometric one.

---

## Network Solver Interface (solver_interface.h)

Fast-path native interface for network solvers with combined residual and Jacobian evaluations.

### Result Types
```cpp
enum class CorrelationValidity : std::uint8_t { VALID, EXTRAPOLATED, INVALID };

template <typename T> struct CorrelationResult {
    T result;
    CorrelationValidity status;
    std::string message;
};

struct Stream {
    double m_dot;
    double Tt;  // total temperature [K] (T before #782)
    double Pt;  // total pressure [Pa] (P_total before #782)
    std::vector<double> Y;
};

struct StreamJacobian {
    double d_mdot;
    double d_Tt;
    double d_Pt;
    std::vector<double> d_Y;
};

struct OrificeResult {
    double m_dot_calc;
    double d_mdot_dP_total_up;
    double d_mdot_dP_static_down;
    double d_mdot_dTt_up;
    std::vector<double> d_mdot_dY_up;
};

struct ChannelResult {
    double dP_calc;
    double d_dP_d_mdot;
    double d_dP_dP_static_up;
    double d_dP_dTt_up;
    std::vector<double> d_dP_dY_up;
};

// The mixed TOTAL state (#782 names; T_mix/P_total_mix before 0.8): Tt_mix
// from the h(Tt) balance, Pt_mix the mass-averaged total pressure. The
// combustor wrappers return the same struct.
struct MixerResult {
    double Tt_mix;
    double Pt_mix;
    std::vector<double> Y_mix;
    double dTt_mix_d_delta_h;
    double dTt_mix_dQ;  // = dTt_mix_d_delta_h * phi (mixer_heat_per_unit_flow)
    std::vector<StreamJacobian> dTt_mix_d_stream;
    std::vector<StreamJacobian> dPt_mix_d_stream;
    std::vector<std::vector<StreamJacobian>> dY_mix_d_stream;
};
```

### Incompressible Flow Components
```cpp
// With a velocity of approach (0 < beta < 1, approach_from_static) the drive
// is P_static_up - P_static_down: E = 1/sqrt(1 - beta^4) already carries the
// approach head. Without one it is P_total_up - P_static_down. Both kernels
// follow the rule.
inline bool approach_from_static(double beta);

// Orifice flow with discharge coefficient. Density at P_static_up, or at
// P_rho_ref when > 0 (its derivative is then d_mdot_dP_rho_ref). Composition
// columns analytic through rho = P MW / (R_u T) (#964).
OrificeResult orifice_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double Tt_up,
    const std::vector<double>& Y_up, double P_static_down, double Cd,
    double area, double beta = 0.0, double P_rho_ref = -1.0);

// d(mu)/dY_k [Pa*s] at fixed (T, P), analytic (Wilke's rule chained to the
// mass fractions; X the state's mole fractions): the composition column of a
// Reynolds-dependent orifice Cd, Re ~ 1 / mu (#964).
std::vector<double> viscosity_mass_gradient(double T, double P, const std::vector<double>& X);

// Channel flow with Darcy friction; exit_head_lost adds the exit dynamic
// head m|m| / (2 rho A^2), dissipated at a stagnant sink (#588)
ChannelResult channel_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double Tt_up,
    const std::vector<double>& Y_up, double P_static_down, double L, double D,
    double roughness, const std::string& friction_model,
    double f_multiplier = 1.0, bool exit_head_lost = false);

// Flow restriction (K-factor)
ChannelResult restriction_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double Tt_up,
    const std::vector<double>& Y_up, double P_static_down, double K);
```

### Momentum Chamber Closure
```cpp
// Pt = Pt_cl(P, Tt, |m_dot| / A), the isentropic stagnation pressure of the
// static state solved from h(Tt) = h(T) + v^2 / 2 (#357, #576); the low-Mach
// Pt = P + G^2 / (2 rho) when incompressible (#588). Residual Pt - Pt_cl with
// analytic partials in P, Pt, m_dot, Tt and Y, plus the static state
// (T_static, M, rho_static, q_dyn) and its partials.
//
// One root per Pt (#788): below the sonic static pressure P* of the flux the
// compressible closure is continued by its mirror image about its minimum,
// 2 Pt* - Pt_cl(P), monotone in P and C1 through P*. `bounded` is true there
// (Pt < Pt* at the root: the area cannot carry the flux at that Pt), with
// P_sonic = P*, Pt_sonic = Pt* [Pa] (zero otherwise); the static-state
// fields stay the closure's at P.
MomentumChamberResult momentum_chamber_residual_and_jacobian(
    double P, double Pt, double m_dot, double Tt,
    const std::vector<double>& Y, double area, bool incompressible = false);
```

### Compressible Flow Components
```cpp
// Compressible orifice flow using isentropic nozzle model
std::tuple<double, double, double, double> orifice_compressible_mdot_and_jacobian(
    double Tt, double Pt, double P_back, const std::vector<double>& X,
    double Cd, double area, double beta);

// Full compressible orifice evaluation for network solver. Tt_up is the
// upstream stagnation temperature. The nozzle expands from the state
// approach_from_static picks: (Tt_up, P_total_up) without E; with E the
// upstream static state, P_static_up at the isentrope's temperature through
// (Tt_up, P_total_up). Composition columns analytic (#964): the nozzle's
// isentrope at fixed (T0, P0, Pb) -- choked, at fixed P* by the envelope
// theorem -- plus the approach static temperature's on its isentrope.
OrificeResult orifice_compressible_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double Tt_up,
    const std::vector<double>& Y_up, double P_static_down, double Cd,
    double area, double beta);

// Opt-in compressible orifice models (#854). OrificeResult's columns plus
// d_mdot_dCd (Jobson is not linear in Cd), Ci, contraction (Jobson's C),
// expansibility (mdot / Cd E A sqrt(2 rho1 dp)), choked and in_range (a flag).
// Analytic (Pt, P1, Tt, Pb, Cd) columns by DualN, and composition through
// n = gamma(Tt_up) and R as two more DualN directions (#964). Jet polytropic,
// n = gamma(Tt_up).
struct OrificeModelResult : OrificeResult { /* d_mdot_dCd, Ci, contraction,
    expansibility, choked, in_range */ };
// Jobson (1955) with the approach-velocity extension, sub- and supercritical
// (Eq. 15 with the same V1 term), C1 through r_c. include/orifice_jobson.h
// holds the one templated body the library functions above share.
OrificeModelResult orifice_jobson_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double Tt_up,
    const std::vector<double>& Y_up, double P_static_down, double Cd, double area,
    double beta, double Cv = 1.0);
// ISO 5167-2:2003 Eq. (1) with Eq. (5); in_range false below p2/p1 = 0.75,
// where Eq. (5) is extrapolated (it has no choke).
OrificeModelResult orifice_iso5167_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double Tt_up,
    const std::vector<double>& Y_up, double P_static_down, double Cd, double area,
    double beta);

// The sonic static pressure P* = G sqrt(R T* / gamma*) of a chamber at
// stagnation Tt and mass flux G (#637), T* from h(Tt) - h(T*) = a*^2 / 2: the
// floor of a junction port's chamber discharging into a pressure boundary.
SonicStatic sonic_static_pressure(double Tt, double G, const std::vector<double>& Y);

// The static pressure of an orifice's (or incompressible channel's) approach
// pipe (#651): the flow at mass flux G from stagnation (Pt, Tt), with its
// derivatives. Compressible: duct_static_state. Incompressible: the
// momentum chamber's low-Mach closure Pt = P + G^2 R Tt / (2 P), root
// [Pt + sqrt(Pt^2 - 2 G^2 R Tt)] / 2; past its maximum flux P = Pt / 2, choked.
// dP_dY analytic (#964), subsonic and past the sonic limit.
ApproachStatic approach_static_pressure(double Pt, double Tt, double G,
                                        const std::vector<double>& Y,
                                        bool incompressible);

// An orifice between pipes (#645): the tap pressure the kernels expand to,
// from the recovered (stagnation) drop between the two pipes' nodes,
//   p_tap = P_static_up - (Pt_up - Pt_down) / k(beta, C),
// with its derivatives (d_dP_total_up/down are d/dPt_up and d/dPt_down; d_dC
// through k). Compressible, floored at
// ORIFICE_TAP_PRESSURE_FLOOR * P_static_up (0.1, below any critical ratio);
// in_range is p_tap / P_static_up >= ORIFICE_RECOVERY_MIN_PRESSURE_RATIO (0.75).
OrificeTapPressure orifice_tap_pressure_and_jacobian(double P_static_up,
                                                     double Pt_up,
                                                     double Pt_down,
                                                     double beta, double C,
                                                     bool compressible);

// The compressible channel's network residual is m - m_calc with m_calc from
// fanno_channel_flow (fanno_mach.h): mass flow a duct passes from its feeding
// node's stagnation state against the other node's pressure, saturating at
// the choked flux, with derivatives. See "Fanno flow integrated in Mach number".
```

### Tee Junction Components (tee_junction.h + solver_interface.h)
```cpp
// Constants
constexpr double TEE_Q_LO     = 0.0;        // lower bound of validated flow ratio
constexpr double TEE_Q_HI     = 1.0;        // upper bound of validated flow ratio
constexpr double TEE_PSI_MIN  = 0.05;       // minimum area ratio psi = A_branch / A_com
constexpr double TEE_THETA_MAX = M_PI / 2;  // maximum branch angle

// Validity check (non-throwing)
struct TeeInputStatus {
    bool q_in_range;
    bool psi_valid;
    bool theta_valid;
    bool valid() const;
    CorrelationStatus status() const;
};
TeeInputStatus tee_check_inputs(double q, double psi, double theta);

// Raw K-coefficient functions (Bassett 2001, Table 2 + Eq 33/34 angle corrections).
// Pure math, no input guards. K2 and K5 are psi/theta-independent per Eq 15.
double K1(double q, double psi, double theta);
double K2(double q);                                   // = q^2 - 1.5*q + 0.5
double K3(double q, double psi, double theta);
double K4(double q, double psi, double theta);
double K5(double q);                                   // same as K2; separating type 3
double K6(double q, double psi, double theta);
double K7(double q, double psi, double theta);
double K8(double q, double psi, double theta);
double K9(double q, double psi, double theta);
double K10(double q, double psi, double theta);
double K11(double q, double psi, double theta);
double K12(double q, double psi, double theta);
// Analytic q-derivatives for each.
double dK1_dq(double q, double psi, double theta);
double dK3_dq(double q, double psi, double theta);
double dK4_dq(double q, double psi, double theta);
double dK5_dq(double q);
double dK6_dq(double q, double psi, double theta);
double dK7_dq(double q, double psi, double theta);
double dK8_dq(double q, double psi, double theta);
double dK9_dq(double q, double psi, double theta);
double dK10_dq(double q, double psi, double theta);
double dK11_dq(double q, double psi, double theta);
double dK12_dq(double q, double psi, double theta);

// Smooth blend helpers
double soft_lower(double x, double lo);              // smooth max(x, lo)
double blend_weight(double r, double k = 30.0);      // 0.5*(1+tanh(k*r))
double blend_weight_deriv(double r, double k = 30.0); // d(blend_weight)/dr

// Blended K functions (safe for any real inputs, never NaN)
double merging_tee_K_straight(double q, double psi, double theta, double blend_k = 30.0);
double merging_tee_K_branch(double q, double psi, double theta, double blend_k = 30.0);
double branching_tee_K_straight(double q, double psi, double theta, double blend_k = 30.0);
double branching_tee_K_branch(double q, double psi, double theta, double blend_k = 30.0);

// Helper utilities
double tee_flow_ratio(double m_dot_branch, double m_dot_com);
bool   tee_topology_valid(double q, double epsilon = 0.05);

// Solver result struct (declared in solver_interface.h)
struct TeeJunctionResult {
    double R_straight, R_branch;
    double dR_straight_d_mdot_com, dR_straight_d_mdot_branch;
    double dR_straight_dP_static_com, dR_straight_dT_com;
    std::vector<double> dR_straight_dY_com;
    double dR_branch_d_mdot_com, dR_branch_d_mdot_branch;
    double dR_branch_dP_static_com, dR_branch_dT_com;
    std::vector<double> dR_branch_dY_com;
    double K_straight, K_branch, q, blend_w;
    bool topology_valid;
    CorrelationValidity status;
};

// Solver interface
TeeJunctionResult merging_tee_residuals_and_jacobian(
    double m_dot_com, double m_dot_branch,
    double dP0_straight, double dP0_branch,
    double P_static_com, double T_com, const std::vector<double>& Y_com,
    double theta, double psi, double F_C, double blend_k = 30.0);

TeeJunctionResult branching_tee_residuals_and_jacobian(
    double m_dot_com, double m_dot_branch,
    double dP0_straight, double dP0_branch,
    double P_static_com, double T_com, const std::vector<double>& Y_com,
    double theta, double psi, double F_C, double blend_k = 30.0);
```

Port convention: MAIN_INLET=A, BRANCH=B, MAIN_OUTLET=C; `m_dot_com` = common-port
mass flow, `m_dot_branch` = branch mass flow, `F_C` = cross-sectional area [m^2] at
the common port, `psi` = A_branch / A_com, `theta` = branch angle [rad].

### Momentum-CV Junction: whole-element (f, J)

Backs `MultiPortChamberElement`. Unlike the tee functions above, this returns the whole
element's residual vector and its full Jacobian from one seeded evaluation --
forward-mode dual numbers over every unknown, so there is no separate
derivation to keep in sync (`include/mpce_junction.h`).

```cpp
struct MpceGeometry {
    std::array<double, 3> area, theta_rad, port_sign;
    double joining_etransfer_alpha, eta_scale;
};

struct MpceResidualJacobian {
    std::array<double, 4> residual;                    // 3 port rows + mass
    std::array<std::array<double, 10>, 4> jacobian;    // [row][seed]
    bool valid;                                        // false: not a junction
    int common_port;
    double k_term_sign;
    std::array<double, 3> k_per_port;
};

MpceResidualJacobian mpce_residuals_and_jacobian(
    const std::array<double, 3>& p_static, const std::array<double, 3>& p_total,
    const std::array<double, 3>& rho,      const std::array<double, 3>& drho_dp,
    const std::array<double, 3>& outer_mdot, double pt_jct,
    const MpceGeometry& geom);
```

Seed order, matching the Jacobian's columns:

| columns | unknown |
|---|---|
| 0-2 | `p_static[i]` |
| 3-5 | `p_total[i]` |
| 6-8 | `outer_mdot[i]` |
| 9 | `pt_jct` |

`outer_mdot` is the CONNECTING element's mass flow, not the junction-convention
one; `geom.port_sign` maps between them, so the Jacobian arrives already
expressed in the solver's own unknowns.

**Physics only.** The kernel does not own the degenerate-state guards or the
wrong-direction soft barrier -- those are solver policy and live in
`MultiPortChamberElement`. It reports `valid = false` for a flow pattern that is not a
junction in any regime rather than guessing.

**Thermodynamics stays at the call site**: pass `rho` and `drho/dP` per port.
A Jacobian is first order, so seeding `rho + (drho/dP)(P - P_val)` is exact.
For combaero's mixtures `drho/dP = rho/P` holds to 1e-11.

Residual sign convention: `R_straight = dP0_straight + K_straight * q_dyn`,
`R_branch = dP0_branch + K_branch * q_dyn`.

### Ideal Momentum Junction: whole-element (f, J)

A parameter-free junction closure (#517), parallel to the Mynard closure above, not
a replacement: a baseline for configurations no tuned closure covers. Nothing in it
is fitted (`include/ideal_junction.h`).

```cpp
struct IdealJunctionGeometry {
    std::vector<double> area, axis_rad, port_sign;  // axis_rad: OUTWARD port axis
    std::vector<double> loss_factor;  // receiver loss factor f_i >= 0; empty = all 1
};

struct IdealJunctionResult {
    std::vector<double> residual;               // N port rows [Pa] + mass [kg/s]
    std::vector<std::vector<double>> jacobian;  // [row][col]
    std::vector<double> velocity;               // u_i, + out of the junction
    double V_x, V_y;                            // supplied mean velocity
};

IdealJunctionResult ideal_junction_residuals_and_jacobian(
    const std::vector<double>& p_static, const std::vector<double>& rho,
    const std::vector<double>& drho_dp, const std::vector<double>& outer_mdot,
    double p_jn, const IdealJunctionGeometry& geom);  // 2..6 ports
```

Rows, with `v_i = u_i n_i` and `V` the supplied mass-weighted mean velocity vector:

| port | row |
|---|---|
| supplier k | `P_k - p_jn = 0` |
| receiver j | `P_j - p_jn + (m_j / A_j)(u_j - V . n_j) = 0` |
| mass | `sum m_i = 0` |

- **Merge** (one receiver): Blaisdell & Manson's theory -- common static pressure,
  momentum projected on the outlet axis, no wall term.
- **Split** (one supplier): every receiver loses `rho/2 |v_j - v_s|^2 >= 0`, the vector
  Borda-Carnot loss -- Idelchik's diverging-branch Eq 7-6 with `A' = 1`; a collinear
  equal-area receiver gives `(1 - q)^2` on its own flow fraction.
- **Compressible:** written in static pressure with each port's own density, so the
  momentum balance is exact at any Mach number for uniform port profiles; there is no
  Mach inversion in the kernel.
- **Roles are smooth.** Every port carries both roles with weights `w_out = m s(m)`,
  `w_in = w_out - m`, `s` a smooth step of width `IDEAL_ROLE_EPS_REL = 1e-4` of the
  throughput. Rows are C-infinity everywhere, the all-zero state included, so there is no
  degenerate-state fallback; the smoothing moves a row by `O((eps/m)^2)` (measured
  <= 5e-6 of the collector head over the B&M merge grid; worst at a supplier carrying 5%
  of the flow at psi = 16).

Columns: `p_static[0..N-1]`, `outer_mdot[0..N-1]`, `p_jn`. Density is passed as values
(`rho`, `drho/dP`), as for the Mynard kernel.

**Loss factor (#526).** Each receiver row rearranges exactly to
`Pt_j = Pt_mix,j - L_j`, with `Pt_mix,j = p_jn + rho_j |V|^2/2` the mixed supplied
stream and `L_j = rho_j/2 |v_j - V|^2 >= 0`. At uniform density the junction's
mass-weighted Pt dissipation is the suppliers' mixing term
`rho/2 sum |m_k| |v_k - V|^2` plus `sum m_j L_j`, both sums of squares.

- `loss_factor` scales `L_j` per port: `Pt_j = Pt_mix,j - f_j L_j`, added to the row
  as `(f_j - 1) s_j L_j`. A supplying port's row is untouched; a port keeps its factor
  when it reverses, and the factor acts only while the port receives.
- `f = 1` (or an empty vector) is the parameter-free closure, bit for bit.
- The factor fades to 1 with the port's flow, `f_eff = 1 + (f - 1) m^2/(m^2 + delta^2)`,
  with `delta = IDEAL_LOSS_FADE_REL = 1e-2` of the throughput.
  - Without it, a receiver carrying vanishing flow sat at `p_jn + (1 - f) rho|V|^2/2`
    while as a supplier it sits at `p_jn`: a step of `(f - 1)` supplied heads across
    reversal, which stalled or misled the solve whenever a branch shut off.
  - A port with share `s` keeps `1 - delta^2/s^2` of the factor's effect: 96% at 5%,
    99% at 10%.
- The mixing term is not scaled, so `f = 0` is lossless only in a split.
- Any `f >= 0` keeps the mechanical dissipation non-negative. Negative or non-finite
  factors are refused.
- `f` is the user's rig-matching knob; the scored closure is `f = 1`.

DOF accounting per junction: +1 unknown (`P_jct`), +N+1 residuals = net +N
equations. Combined with port-MCN mass-row skip in the solver (the junction's
sum-mass residual replaces the N port-MCN continuity rows), the system stays
square.

### Mixing and Combustion Components
```cpp
// Stream mixing with optional heat transfer. T is the exact inverse of the
// mixed enthalpy down to 0 K (the gas is ideal there); below h(0 K), a heat
// withdrawal past what the streams hold, T continues on h's line through
// 0 K (C1) instead of raising (#552), and the network reports a node at
// T <= 0 in __T_clamped__. The combustor forms do the same for their
// post-combustion heat shift.
MixerResult mixer_from_streams_and_jacobians(
    const std::vector<Stream>& streams, double Q = 0.0, double fraction = 0.0);

// Complete combustion of mixed streams, then a heat shift: Q [W] and
// 'fraction' of the PRODUCTS' sensible enthalpy at T_ad (#481)
MixerResult adiabatic_T_complete_and_jacobian_T_from_streams(
    const std::vector<Stream>& streams, double P, double Q = 0.0, double fraction = 0.0);

// The same with reforming + water-gas-shift equilibrium of the products
MixerResult adiabatic_T_equilibrium_and_jacobians_from_streams(
    const std::vector<Stream>& streams, double P, double Q = 0.0, double fraction = 0.0);

// Second law of a control volume (#545). Streams at their stagnation state
// (T = total temperature); S_heat [W/K] the entropy the heat brought in.
// sigma = sum_out m s - sum_in m s - S_heat; sigma_norm = sigma / (R m_in),
// exactly -ln(Pt_out / Pt_in) for an adiabatic two-port.
EntropyGenerationResult entropy_generation(
    const std::vector<Stream>& inflows, const std::vector<Stream>& outflows,
    double S_heat = 0.0);

// integral dQ / T [W/K] as a stream's total temperature runs Tt_in -> Tt_out,
// at the STATIC temperature (Mach linear in Tt from M_in to M_out). At M = 0
// exactly m [s(Tt_out, P) - s(Tt_in, P)]; the Mach part by composite Simpson.
double heat_entropy_flux(double m_dot, double Tt_in, double Tt_out,
                         const std::vector<double>& Y, double M_in = 0.0,
                         double M_out = 0.0);

// All three add Q as Q / m_tot exactly (energy-conserving) for total flows of
// at least kMixerHeatMdotFloor = 1e-6 kg/s; below it Q phi(m) with
// (phi, dphi/dm) = mixer_heat_per_unit_flow(m): the heat taken up,
// Q (4 u^3 - 3 u^4) with u = |m| / m0, fades C1 to none at zero flow, so a
// stagnant node holds the temperature of what reaches it.
std::tuple<double, double> mixer_heat_per_unit_flow(double mdot_tot);

// Combustion of mixed streams with a pressure-loss correlation
MixerResult combustor_residuals_and_jacobians(
    const std::vector<Stream>& streams, double P, double Q, double fraction,
    const PressureLossCorrelation& pressure_loss, bool use_equilibrium = false);
```
