# CombAero Python API Reference

This document provides the high-level reference for the `combaero` Python package. For the C++ technical reference, see [API_CPP.md](API_CPP.md).

## Table of Contents
- [Core Thermodynamics](#core-thermodynamics)
- [Combustion & Equilibrium](#combustion--equilibrium)
- [Flow Regimes (Symmetric API)](#flow-regimes-symmetric-api)
- [Advanced Flow Functions](#advanced-flow-functions)
- [Comprehensive Orifice Functions](#comprehensive-orifice-functions)
- [Heat Transfer](#heat-transfer)
  - [Correlations](#correlations)
  - [Network Heat Transfer](#network-heat-transfer)
- [Acoustics](#acoustics)
- [Network Solver](#network-solver)
- [NetworkRunner (GUI-JSON Programmatic Driver)](#networkrunner-gui-json-programmatic-driver)
- [Geometry & Materials](#geometry--materials)
- [Advanced Thermodynamics](#advanced-thermodynamics)
- [Psychrometrics (Humid Air)](#psychrometrics-humid-air)

---

## Core Thermodynamics

### State Object
The `State` object is the central data structure, managing T [K], P [Pa], and composition.

```python
import combaero as cb

state = cb.State()
state.TPX = (300.0, 101325.0, "O2:0.21, N2:0.79")
# OR setting mass fractions directly
state.Y = y_vector

print(f"Mole fractions: {state.X}")
print(f"Mass fractions: {state.Y}")

print(f"Density: {state.rho} kg/m³")
print(f"Enthalpy: {state.h} J/mol")
print(f"Viscosity: {state.mu} Pa·s")
```

### MixtureState (Network Data)
The `MixtureState` struct is used for network node properties and carries both static and stagnation conditions.
**Units**: P [Pa], T [K], m_dot [kg/s].

```python
# Constructor: MixtureState(P, Pt, T, Tt, m_dot, Y); P, T static, Pt, Tt total
node_state = cb.MixtureState(1e5, 1.1e5, 300, 300, 1.5, y_vec)

print(node_state.P)       # Static [Pa]
print(node_state.m_dot)   # Mass flow [kg/s]
print(node_state.density()) # Static density [kg/m³]
```

### Species Module
The `combaero.species` module provides high-level helpers for composition vectors.

```python
import combaero as cb

# Get standard compositions (returns numpy array)
X_air = cb.species.dry_air()        # Mole fractions
Y_air = cb.mole_to_mass(X_air)   # Mass fractions (standard for network)

# Humid air at T [K], P [Pa], RH [0-1]
Y_humid = cb.species.humid_air_mass(300, 101325, 0.5)

# Pure species or from mapping
Y_ch4 = cb.species.pure_species("CH4")
Y_mix = cb.species.from_mapping({"O2": 0.20, "N2": 78, "AR": 0.02})

# Explicit conversion
Y = cb.species.to_mass(X_air)
X = cb.species.to_mole(Y_air)
```

### Species Lookup (Legacy)
```python
# Low-level metadata from core
num = cb.num_species()
name = cb.species_name(0)
idx = cb.species_index_from_name("H2O")
```

### Species Thermo Model (#549)
```python
cb.thermo_model()               # "nasa9_c2" (default): each interval join bridged by a C2 cp
cb.set_thermo_model("nasa9")     # the published NASA-9 polynomials
cb.set_thermo_model("nasa9_c2")
```
NASA-9 is C0 in cp at its interval joins (1000 K; 6000 K for N2, O2, AR, CO2,
CO, H2) but not C1: dcp/dT jumps there, 0.15-100 % at 1000 K and up to 48 % at
6000 K, and every solver Jacobian that uses it steps with it. `"nasa9_c2"`
bridges each join inside |T - T_join| < 0.05 T_join with a degree-7 cp/R that
is C2 with NASA-9 at both window edges and whose h and s rises across the window
equal NASA-9's: everything is identical to NASA-9 outside the windows (formation
enthalpies, absolute entropies, 298.15 K reference), and inside cp moves by at
most 4.8e-4 (1e-4 at 1000 K), h by 1.1e-5 cp T, s by 3.1e-5 R. The setting is
process-wide; switch before building and solving a network, not during. The
re-fit is the default (#549): with it the full suite passes and converges in
the same evaluations, so `"nasa9"` stays for comparison against the published
data.

### Beyond the Data (#606)
The data runs from 200 K to 20 000 K (N2, O2, AR, CO2, CO, H2) or 6000 K (H2O,
NH3, CH4, the alkanes), CEA's own ranges. Past either end, under both models,
cp/R fades C2 from the end interval to a constant: over 0.25 T_lo (200 -> 150
K) below and 0.1 T_hi above (6000 -> 6600 K, 20 000 -> 22 000 K), h and s its
exact integrals. cp stays positive and h, s, u monotone at every T, for the
solvers and the inverses -- not physics (no dissociation, no condensation).
Before, the end polynomials ran on: cp/R -25 for CH4 at 50 K, -241 for H2O at
20 000 K. `calc_T_from_h/_s/_u` (and their `_mass` and `sv` forms) raise
`RuntimeError` for a target no positive T reaches instead of returning the last
iterate.

### Inner Solves Report a Cap Exit (#868)
A kernel's inner iterative solve that reaches its iteration cap without
meeting its tolerance raises `RuntimeError` (naming the function and its last
residual) instead of returning the last iterate as if it were a root:
`Tt_from_static` / `Tt_from_static_v`,
`T_from_stagnation`, `Pt_from_static`, `P_from_stagnation`,
`stagnation_from_static`, the isentropic expansion behind
`mach_from_pressure_ratio`, `mass_flux_isentropic` and `nozzle_flow`,
`friction_colebrook`, `fanno_channel` / `fanno_channel_rough`, the static-state
and merge-face kernels the network elements call, and the Fanno/Rayleigh
choke searches (`fanno_rayleigh_channel_flow` with no sign change to bracket).
With the defaults these exits are not reached on physical inputs: what raises
is a caller's own small `max_iter`, an unreachable target, or a solver probe at
an unphysical state, which the network solver already treats as a rejected
step. The static temperature has no floor (#935): an
expansion whose isentropic T lies below the former 10 K `T_FLOOR` (P/Pt below
~1e-5 from 300 K) returns that root, not 10 K, and so do `T_from_stagnation`
at hypersonic M and `Tt_from_static` from a static T below 10 K; the mass flux
and Mach built on it follow.
`fanno_inlet_mach` below ~1e-10 of the sonic flux now returns the root (it was
up to 73 % high: the slope stencil crossed M = 0); its Newton now takes the
exact flux slope. `T_from_stagnation` (and `P_from_stagnation`) keeps the
T dgamma/dT term in its slope and converges to rounding (results move by up to
~1e-9 relative); `State.set_VH_mass` is the enthalpy inverse (u + R T = h, as
`set_HP_mass`). `calc_T_from_cp` returns the root nearest `T_guess` in
[1, 6000] K. Above 6000 K the tails level cp off at a constant, so cp no longer
identifies T; below the 200 K NASA-9 data the ideal-gas tail applies (#606),
where cp - cp(0 K) ~ T^3, so a root near 1 K is ill-conditioned (~5e-7
relative for CO, ~5e-4 at 0.1 K, hence the 1 K end). cp(T) need not be
monotone (N2 has a minimum near 214 K and a tail maximum near 127 K), so the
guess picks the branch: the range is scanned outward from the guess (clamped
into it), each step on the side whose next node is nearer it, 10 K apart and
10 % of T below 100 K, and the first sign change polished. No root in the
range raises `RuntimeError` ("no root in [1, 6000] K"); a non-finite or
non-positive `cp_target`, `T_guess` or `tol`
raises `ValueError`. To get every root instead,
`calc_T_from_cp_roots(cp_target, X, T_min=1.0, T_max=6000.0)` returns a sorted
list of `CpRoot` (`T` [K], `dT_dcp` = 1/(dcp/dT) [K/(J/(mol K))],
`double_root`), empty when there is none (#963). It splits cp at its extrema
(including a pair inside one scan step), so it finds the close pairs the
nearest-root scan misses; a target equal to cp at an extremum (to 1e-13) is
one root with `double_root=True` and `dT_dcp` NaN. An invalid input or range,
and argon's exactly constant cp (roots not isolated), raise `ValueError`. The audit of every inner solve, with what was left
unchanged: `docs/archive/AUDIT_868_INNER_SOLVES.md`.

### Compression and Expansion (#893)
The outlet total temperature of a compressor or turbine at a total-to-total
pressure ratio, with exactly one efficiency in (0, 1]. The composition is frozen.

```python
r = cb.compression(T1, X, PR, eta_poly=0.9)   # PR = P2/P1 >= 1
r = cb.compression(T1, X, PR, eta_is=0.86)
t = cb.expansion(T1, X, PR, eta_poly=0.9)     # PR = P1/P2 >= 1
r.T2, r.dh        # outlet Tt [K]; h_mass(T2) - h_mass(T1) [J/kg], > 0 compressing
r.eta_poly, r.eta_is, r.T2s   # both efficiencies of the path; the isentropic T2
r.dT2_dT1, r.dT2_dPR, r.dT2_deta, r.dT2_dY    # the Jacobian of T2 (and ddh_*)
```

- **Polytropic is exact** for combaero's variable cp, with no path
  integration: S(T2) − S(T1) = (R/η_p) ln PR for compression and
  −η_p R ln PR for expansion, with S = Σ X_k s_k at the reference pressure.
- **Isentropic** uses the end states: h2 = h1 + (h2s − h1)/η_s for
  compression and h1 − η_s (h1 − h2s) for expansion.
- **The pressure level drops out.** At a fixed PR, T2 does not depend on
  the inlet pressure.
- **Refused:** PR < 1, an efficiency outside (0, 1] (above 1 the path would
  destroy entropy), and neither or both efficiencies given.

Checked against:
- Cantera on the same NASA-9 data. Isentropic, through Cantera's SP/HP setters:
  3e-8 on humid air at PR 2-40 and on combustion products at PR 2-20.
  Polytropic, against the path ODE on Cantera's cp: 1.6e-9.
- The constant-gamma closed forms on argon (cp = 5/2 R).
- The limit of n isentropic stages, which converges at first order.


---

## Combustion & Equilibrium

### Complete Combustion
Fuel C/H to CO2 + H2O(g), no dissociation; rich, every fuel burns 1/phi and
the rest stays unburnt (no CO/H2). Model and smoothing: the Combustion note in
[API_CPP.md](API_CPP.md). Invalid input (NaN, T or P <= 0, non-finite phi or
target) raises `ValueError`.
```python
X = cb.set_equivalence_ratio_mole(0.8, X_fuel, X_air)   # phi > 0
# ADIABATIC flame temperature and products. The Python defaults smooth both
# phi = 0 and phi = 1 (smooth_phi0=smooth_phi1=True; the C++ default is off).
burned = cb.complete_combustion(T=300.0, X=X, P=101325.0)
print(f"Adiabatic T: {burned.T} K")
exact = cb.complete_combustion(300.0, X, 101325.0, smooth_phi0=False, smooth_phi1=False)

# Isothermal products, the bare product vector, and with dX_b/dX
products = cb.complete_combustion_isothermal(300.0, X, 101325.0)
X_b = cb.complete_combustion_to_CO2_H2O(X)               # smooth=False (phi = 0 only)
X_b, f = cb.complete_combustion_to_CO2_H2O_with_fraction(X)
jac = cb.complete_combustion_products_and_jacobian(X)    # .X_b, .dXb_dX[j][k]

phi = cb.equivalence_ratio(X)                            # elemental, all atoms
phi = cb.equivalence_ratio_mole(X, X_fuel, X_air)        # two-stream
lhv = cb.fuel_lhv_mass(X_fuel)                           # J/kg, H2O vapour
st = cb.combustion_state(X_fuel, X_air, phi=1.2, T_reactants=300.0, P=1e5)
st.fuel_burn_fraction                                    # 1/1.2 for Complete
```
Also: `set_equivalence_ratio_mass`, `bilger_stoich_mixture_fraction_mass`,
`bilger_Z_from_equivalence_ratio_mass`, `equivalence_ratio_from_bilger_Z_mass`,
`oxygen_required_per_*`, `combustion_state_from_streams`,
`set_fuel_stream_for_phi` and the `set_{fuel,oxidizer}_stream_for_{Tad,O2,O2_dry,CO2,CO2_dry}`
inverse solvers.

### Chemical Equilibrium
```python
# T, P constant equilibrium
eq_state = cb.combustion_equilibrium(T=300.0, X=X, P=101325.0, smooth=True)

# Water-gas shift
wgs = cb.wgs_equilibrium(state)
```

---

## Flow Regimes (Symmetric API)

CombAero provides symmetric submodules for incompressible and compressible flow.

```python
from combaero import incompressible as flow   # OR: from combaero import compressible as flow

sol = flow.channel_flow(T=400, P=2e5, X=air, u=10.0, L=2.0, D=0.05, f=0.02)
print(f"Pressure drop: {sol.dP} Pa")
print(f"Outlet Mach: {sol.M}") # nan for incompressible
```

### FlowSolution Fields
- `mdot`: Mass flow rate [kg/s]
- `v`: bulk velocity [m/s]
- `dP`: Pressure drop [Pa]
- `M`: Mach number
- `choked`: Boolean flag

---

## Advanced Flow Functions

### Inverse Solvers
```python
# Solve for unknown given mass flow rate (invert nozzle_flow)
A_eff = cb.solve_A_eff_from_mdot(Tt=300, Pt=2e5, P_back=1e5, mdot_target=0.1, X=air)
# Subsonic root; RuntimeError above the choked flow (0.0467 kg/s here)
P_back = cb.solve_P_back_from_mdot(Tt=300, Pt=2e5, A_eff=1e-4, mdot_target=0.04, X=air)
Pt = cb.solve_Pt_from_mdot(Tt=300, P_back=1e5, A_eff=1e-4, mdot_target=0.1, X=air)
```

Non-positive temperatures, pressures, areas or target flows raise `ValueError`.
`fanno_max_length` integrates the constant-`f` Fanno flow in Mach number up to
`M = 1` (the `fanno_length_between` quadrature); it raises `ValueError` for
`u_in <= 0` or `f <= 0`, where `L*` is infinite.

### Utility Functions
```python
# Critical pressure ratio and Mach number
P_crit_ratio = cb.critical_pressure_ratio(Tt=300, Pt=2e5, X=air)  # P*/Pt
Mach = cb.mach_from_pressure_ratio(Tt=300, Pt=2e5, P=1.5e5, X=air)

# Mass flux and flow analysis
G = cb.mass_flux_isentropic(Tt=300, Pt=2e5, P=1.5e5, X=air)
L_max = cb.fanno_max_length(T_in=300, P_in=2e5, u_in=50, D=0.05, f=0.02, X=air)
```

A duct's inlet static state from the stagnation state feeding it and its mass
flux (#636), with partials in Tt, Pt and G; the sonic state above the choked
flux (`choked`). Network ducts evaluate their heat transfer here, so their
T_aw lies between the static and the stagnation temperature:

```python
s = cb.duct_static_state(Tt=1600, Pt=2e5, G=150.0, X=air)
s.T, s.P, s.rho, s.M, s.choked, s.dT_dTt, s.dP_dG
```

Its composition gradients (#717), one entry per species in the MASS
fractions, at fixed Tt, Pt and G, for the wall relay behind a combustor:

```python
s = cb.duct_htc_static_state(Tt=1600, Pt=2e5, G=150.0, X=air)
d = cb.duct_htc_static_state_dY(Tt=1600, Pt=2e5, G=150.0, X=air, state=s)
d.dT_dY, d.dP_dY, d.drho_dY, d.dGeff_dY
```

### Fanno flow in Mach number
Integrated in M rather than x: the state is algebraic in M at fixed mass flux,
so the duct length is a smooth integral that vanishes at sonic -- exact L*,
no cutoff. Fed isentropically from a stagnation state:

```python
G_ch = cb.fanno_choked_mass_flux(Pt=2e5, Tt=300, X=air, L=1.0, D=0.02,
                                 roughness=1e-5, friction_model="haaland")  # kg/(m^2 s)
r = cb.fanno_duct(2e5, 300, 0.8 * G_ch, air, 1.0, 0.02, 1e-5, "haaland")
r.inlet.M, r.exit.M, r.exit.P, r.Pt_exit, r.L_star, r.choked
cb.fanno_length_between(G=300, Tt=300, M1=0.3, M2=1.0, X=air, D=0.02,
                        roughness=0.0, friction_model="fixed", f_multiplier=0.02)
```
Agrees with the x-march (`fanno_channel_rough`, 20000 steps) to 1e-9 at
constant f and to ~1e-8 with a Re-dependent f, the floor set by kinks in
mu(T).

### Rocket Nozzle Functions
```python
# Converging-diverging nozzle analysis
sol = cb.nozzle_cd(Tt=300, Pt=2e5, P_exit=1e5,
                   A_inlet=1e-3, A_throat=5e-4, A_exit=2e-3,
                   x_throat=0.1, x_exit=0.2, X=air)
sol.Tt, sol.Pt, sol.mdot, sol.choked   # stagnation state the nozzle expands from

# Thrust calculations
thrust = cb.nozzle_thrust_cd(Tt=300, Pt=2e5, P_design=1e5, P_amb=101325,
                             A_inlet=1e-3, A_throat=5e-4, A_exit=2e-3,
                             x_throat=0.1, x_exit=0.2, X=air)
print(f"Thrust: {thrust.thrust} N")
print(f"Specific impulse: {thrust.specific_impulse} s")
```

---

## Comprehensive Orifice Functions

### Discharge Coefficient Correlations
```python
# Normed metering orifice (ISO 5167 family, a function of beta = d/D)
geom.tap = cb.OrificeTap.Corner  # Corner, D_D2 or Flange (default)
Cd = cb.Cd_sharp_thin_plate(geom, state)  # Reader-Harris/Gallagher at geom.tap
Cd, dCd_dRe = cb.Cd_sharp_thin_plate_and_dRe(geom, state)  # analytic Re_D derivative
Cd, dCd_dRe = cb.metering_cd_and_dRe(cb.MeteringCdCorrelation.ReaderHarrisGallagher, geom, Re_D)
# MeteringCdCorrelation.Stolz (ISO 5167:1980) was removed in 0.9.0 (#876).
ok = cb.metering_cd_in_range(cb.MeteringCdCorrelation.ReaderHarrisGallagher, geom, Re_D)
# RHG's lower Re_D bound (ISO 5167-2 5.3.1). Below it the RHG Cd is HELD (#873,
# 0.8.0), entered C1 in ln Re_D over 0.1 just below the bound; metering_cd_in_range
# is False there. It used to climb to a 1.5 cap.
Re_b = cb.metering_cd_re_min(cb.MeteringCdCorrelation.ReaderHarrisGallagher, geom)
# ISO 5167-2 permanent loss / tapping differential between pipes, and dk/dC
k, dk_dC = cb.permanent_loss_ratio_and_dC(beta, Cd)
E = cb.velocity_of_approach_factor(beta)  # 1/sqrt(1-beta^4); 1 outside 0<beta<1

# Jobson (1955) compressible contraction of a sharp-edged orifice (#854)
f = cb.jobson_force_defect(Ci, beta=0.0, Cv=1.0)  # 1/Ci - 1/(2 Ci^2) - Cv b^2 + Cv^2 b^4/2
Ci = cb.jobson_ci_from_cd(Cd, beta, Cv=1.0)  # Cv Ci / sqrt(1 - (Cv Ci)^2 b^4) = Cd E
C = cb.jobson_contraction_coefficient(r, Ci, n=1.4)  # Eqs. 11 / 15, r = p/p0, reservoir-fed
eps = cb.jobson_expansibility(r, beta, Ci, 1.4)  # vs ISO 5167-2 Eq. (5), r = p2/p1
# (approach-velocity consistency by bracketed Newton: any r in (0, 1], #957)

# Discharge hole in a wall -- no pipe, no beta. Separate selector, separate
# geometry. Idelchik (1966) diagrams 4-17/4-18, valid Re 25 to 1e6.
hole = cb.DischargeHoleGeometry(d=1e-3, L=2e-3, r=0.0)
flow = cb.DischargeHoleState(Re=5e4)
Cd = cb.discharge_cd(cb.DischargeCdCorrelation.Idelchik1966Thick, hole, flow)

# Solver-facing form: (Cd, dCd/dRe, dCd/d(U1_over_Vi)), analytic throughout
Cd, dCd_dRe, dCd_dU = cb.discharge_cd_and_derivatives(
    cb.DischargeCdCorrelation.Idelchik1966Thick, hole, flow
)

# Inside the source's own Re and l/d range? A flag, never a selector.
cb.discharge_cd_in_range(cb.DischargeCdCorrelation.Lichtarowicz1965, hole, flow)

# Re below which the Cd VALUE is held (true dCd/dRe = 0 there; the solver
# form above continues the slope for McGreehan-Schotsch and Idelchik)
re_floor = cb.discharge_cd_re_floor(cb.DischargeCdCorrelation.McGreehanSchotsch1988)  # 1e4

# An inclined hole, OPT-IN: Idelchik diagram 3-2's
# inlet loss on the thick-walled hole's zeta', dzeta = 0.3 cos + 0.2 cos^2 of
# the hole-to-wall angle (0 at 90). Idelchik1966Thick only; the others raise.
# Flagged out of range below Re 1e4 or l/d 2.
tilted = cb.DischargeHoleGeometry(d=1e-3, L=4e-3, inclination_deg=30.0)
Cd30 = cb.discharge_cd(cb.DischargeCdCorrelation.Idelchik1966Thick, tilted, flow)
dzeta = cb.discharge_inclination_dzeta(30.0)  # 0.41

# The same hole as an orifice in a straight pipe (Idelchik 4-13..4-16):
# permanent loss on the pipe velocity head, and the Cd factor the network
# applies between pipes (C_pipe = Cd_wall * R, R -> 1 as beta -> 0)
zeta1 = cb.discharge_pipe_zeta(cb.DischargeCdCorrelation.Idelchik1966Thick, hole, 0.5)
R = cb.discharge_pipe_cd_ratio(cb.DischargeCdCorrelation.Idelchik1966Thick, hole, 0.5)
```

There is deliberately no auto-selection from geometry: it is how a
rounded-entry request used to come back as Stolz.

#### Plenum-to-plenum holes: McGreehan and Schotsch (1988)

The correlations above are ISO 5167 pipe metering -- a single hole IN a pipe
run. For a cooling transfer hole or a jet-plate hole discharging
plenum-to-plenum, use the composite chain of McGreehan and Schotsch (1988),
ASME J. Turbomachinery 110(2), 213-217, which covers Reynolds number, inlet
corner radius, orifice length and inlet crossflow:

```python
# A bare, sharp-edged jet plate, plate thickness = hole diameter
cb.mcgreehan_schotsch_1988_cd(Re=1e4, r_over_d=0.0, L_over_d=1.0)   # 0.785

# A radiused, long cooling hole
cb.mcgreehan_schotsch_1988_cd(Re=3e4, r_over_d=0.1, L_over_d=3.0)   # 0.889

# Fed from a duct rather than a plenum: supply-side crossflow
cb.mcgreehan_schotsch_1988_cd(Re=3e4, r_over_d=0.0, L_over_d=1.0,
                              U1_over_Vi=0.5)                       # 0.770
```

`U1_over_Vi` is the **inlet (approach, supply-side)** tangential velocity over
the ideal through-flow velocity, and defaults to 0 for a plenum-fed plate. It
is *not* a discharge-side crossflow ratio -- do not pass Florschuetz's
`Gc/Gj`, which is spent air in the impingement channel on the far face of the
plate. `V_i` must be built from static inlet conditions, not from a total
pressure that already contains the tangential velocity head.

Two behaviours that look like bugs and are not:

- **Not monotonic in `U1_over_Vi`.** `Cd` rises up to ~5.7% above its
  zero-crossflow value near `U1_over_Vi ~ 0.09` before falling away. This is
  drawn in the source's Fig. 4 and supported by Rohde (NASA TN D-5467), whose
  result 2 is that slanting an orifice into the flow increases `Cd`.
- **Constant below `Re = 1e4`.** That is the correlation's stated floor, and
  below it the underlying equation diverges (it reaches `Cd = 1.0` at
  `Re = 904`), so the value is held at the floor rather than extrapolated.

When the plate's zero-crossflow `Cd` is already known -- measured, or a
literature value such as Florschuetz's per-configuration Table 1 -- apply only
the crossflow correction:

```python
cb.mcgreehan_schotsch_1988_crossflow_cd(cd_base=0.79, U1_over_Vi=0.5)  # 0.794
cb.mcgreehan_schotsch_1988_crossflow_cd(cd_base=0.79, U1_over_Vi=2.0)  # 0.541
```

Note the first: at a 0.79 baseline, `U1/Vi = 0.5` still sits inside the rise,
because a higher baseline pushes the peak to larger `U1/Vi` (`R_v` carries a
`(Cd/0.6)^-3` factor). The decay is well established by `U1/Vi = 2`.

This is how the source uses Eq. (17) in its own validation. Its Figs. 5 and 6
anchor to "a set baseline point at `U1/Vi = 0`" taken from Rohde's
measurements (0.64, 0.73, 0.88), not to the chain's prediction for the same
geometry -- the chain runs 3-12% higher, which is the paper's own remark that
Rohde's "basic values are lower".

**Measured agreement.** Against Rohde's data as replotted in Fig. 6
(`t/d = 0.51`, baseline 0.64, 11 points), Eq. (17) scores bias `+5.95%`,
RMS `7.54%` -- reading high, and increasingly so with crossflow (`+3%` below
`U1/Vi = 0.4`, `+16%` at 1.41). Scored by
`validation/cooling/orifice_runner.py`.

#### Film cooling: Baldauf et al. (2002)

Laterally averaged adiabatic effectiveness downstream of one row of
cylindrical, streamwise-inclined holes -- valid from the ejection point, and
carrying the adjacent jet interaction rather than excluding it.

```python
eta = cb.film_effectiveness_baldauf_2002(
    x_over_D=20.0, M=2.0, P=1.2, alpha_deg=30.0, s_over_D=3.0, Tu=0.015
)

# Solver form: (eta, deta/dM, deta/dP), analytic
eta, deta_dM, deta_dP = cb.film_effectiveness_baldauf_2002_and_derivatives(
    20.0, 2.0, 1.2, 30.0, 3.0, 0.015
)
```

`eta` uses the same convention as `adiabatic_wall_temperature`, so they
compose. `deta/dM` changes sign through the lift-off peak, which is physical
and is what a solver needs to push blowing the right way.

Envelope: `M` 0.2-2.5, `P` 1.2-1.8, `s/D` 2-5, `alpha` 30-90 deg,
`Tu` 0.0035-0.075; the paper's RMS deviation is 5.5%. Eq. (31) is
implemented as printed and disagrees with the paper's own worked example --
see `validation/cooling/extractions/baldauf_2002_film_effectiveness.md`.

`b_0=<float>` replaces Eq. (31)'s `b_0` (default `None` = as printed), so the
paper's Table 4 reading can be tested. Its effect is governed by `x/D`, not
`M`: under 1% for `x/D <= 20` at any blowing, -28% by `x/D = 400`.

Outside the envelope the value is unchanged -- it still answers, finitely --
but each out-of-range input now raises a warning through the global handler,
one per parameter:

```python
with cb.suppress_warnings():          # or set_warning_handler to collect them
    eta = cb.film_effectiveness_baldauf_2002(20.0, 3.0, 1.0, 30.0, 7.37, 0.05)
# unsuppressed, that call reports M, P and s_over_D as extrapolated
```

#### Film cooling: Colban, Thole and Bogard (2011), shaped holes

Laterally averaged adiabatic effectiveness downstream of one row of
**shaped** (fan or laidback fan) holes. Baldauf above is for cylindrical
holes; they are not interchangeable -- at `M = 2.5` a cylindrical jet has
lifted off and a shaped one has not, and the paper's Fig. 7 has Baldauf
near 0.04 where the shaped-hole data run 0.15-0.5.

```python
eta = cb.film_effectiveness_colban_2011(
    x_over_D=10.0, M=1.5, area_ratio=3.5, P_over_D=6.0, t_over_P=0.49
)  # alpha_deg=30.0, density_ratio=1.7 by default

# Solver form: (eta, deta/dM), closed form
eta, deta_dM = cb.film_effectiveness_colban_2011_and_derivatives(
    10.0, 1.5, 3.5, 6.0, 0.49
)
```

`D` is the **metering** (inlet) diameter and `M` is referenced to the
metering velocity; `area_ratio = A_exit/A_inlet`, `t_over_P` is breakout
width over pitch. `eta = t/P` at the breakout (Eq. 12). `alpha_deg` and
`density_ratio` do not enter Eq. (19) -- every data set was at 30 deg and
DR 1.7-2.0 -- so they only drive the range check.

Limits (Table 4, warned not enforced): `M` 0.2-2.5, `t/P` 0.31-0.65,
`AR/(M P/D)` 0.17-1.17; fitted spans `AR` 2.5-4.7, `P/D` 4-8. On its own
Figs. 6-7 data (Gritsch configuration B, which it was fitted on) it scores
3.6% RMS at `M = 0.5` and 12.1% at `M = 2.5`. See
`validation/cooling/extractions/colban_2011_shaped_film.md`.

#### Film cooling: Marek and Tacina (1975), 2-D slot

Effectiveness downstream of a continuous **slot**, as in a combustor
liner's cooling ring; Baldauf and Colban are for rows of holes and do not
apply. NASA TN D-7958: `eta = 1/(1 + C_m x/(M S))` (Eq. 4) with
`C_m = 0.25 Tu/100` (Eq. 5, Tu in percent), so the free-stream turbulence
sets the mixing. `Tu` is passed as a FRACTION (0.14 = 14%), like every
turbulence input in combaero (#773); above 1 raises.

```python
eta = cb.slot_film_effectiveness_marek_tacina_1975(x_over_S=15.0, M=1.0, Tu=0.14)

# Solver form: (eta, deta/dM), closed form
eta, deta_dM = cb.slot_film_effectiveness_marek_tacina_1975_and_derivatives(15.0, 1.0, 0.14)

# Matching a rig is the user's job: scale C_m (Eq. 5's band is 0.4-1.6)
eta = cb.slot_film_effectiveness_marek_tacina_1975(15.0, 1.0, 0.14, C_m_multiplier=1.3)
```

`M = rho_s U_s/(rho_H U_H)`; `x/S` is measured from the slot lip, where
`eta = 1`. Warned, not enforced, outside the run: `Tu` 0.07-0.35, `M` 0.41-1.87,
`x/S` up to 33.8. Every test was co-flowing film.

Scored on the report's own Table I (140 cooled points, the 72% plate's
16.20 cm station not printed): with the authors' per-blockage C_m the kernel
puts 136 of 140 points inside the +-30% the report states (p. 8); with
C_m from Eq. (5) at the measured Tu the eta RMS is 0.070, 57 of 140 inside
the authors' +-0.04 measurement band. See
`validation/cooling/extractions/marek_tacina_1975_d7958.md`.

#### Film-side h near a slot: Seban (1960) (#962)

The heat transfer coefficient UNDER a slot film, on the adiabatic-wall
temperature of the same flow (`q = h (T_aw - T_w)`, T_aw from the film
effectiveness). Near the slot the boundary layer is the slot fluid's, so h
follows the slot flow and the distance from the lip, not the duct's flow.

```python
# Seban 1960 Fig. 12, St_s = h/(rho_s u_s cp_s), span mean over [x_start, x_end]/s
st = cb.slot_film_stanton_seban_1960(0.0, 25.0, Re_s=5000.0)
st.St, st.dSt_dRe, st.branch, st.fraction_in_box, st.x_in_box, st.Re_in_box

# The wall side's h, dimensional, with the hand-back to the duct's h
r = cb.slot_film_htc_seban_1960(x_start=0.0, x_end=0.05, s=2e-3, G_s=60.0, M=2.0,
        h_duct=150.0, T_s=600.0, P_s=1e6, X_s=cb.species.dry_air())
r.h, r.dh_dG_s, r.dh_dT_s, r.dh_dM, r.dh_dh_duct, r.h_film, r.Re_s, r.weight_M

# The book's m > 1.3 wall jet (named option): h s / k_s, and dimensional
n = cb.slot_film_nusselt_ballal_lefebvre_1973(18.0, 18.0, Re_s=1.81e4)
r = cb.slot_film_htc_ballal_lefebvre_1973(0.0261, 0.0261, 1.45e-3, 485.7, 2.54,
        0.0, 880.0, 3.04e6, cb.species.dry_air())
```

**What it represents.** Seban's two curves on slot-fluid properties
(`Re_s = rho_s u_s s / mu_s`):

- near the slot, curve B: `St_s = 0.135 Re_s^-0.3 (x/s)^-0.3`;
- downstream, curve C: `St_s = 0.40 Re_s^-0.3 (x/s)^-0.6`;
- joined by a C1 blend in `ln(x/s)` over +-0.1 decade around their crossing
  (x/s 37.4), within 1.5 % of the lower curve there.

**0.135, not the printed 0.10.** Seban's Eq. 13 prints 0.10 and Eq. 14 0.41.
The 0.10 contradicts his own curve B, his Fig. 9 and his data. The drawn
B and C cross at x/s 37, where his text switches; the printed pair would
cross at 110. Fig. 9 back-solves to 0.140-0.144, and the digitised Fig. 12
points fit 0.134. So the constants are the curves as DRAWN:
`cb.SEBAN1960_C_NEAR` 0.135 and `cb.SEBAN1960_C_FAR` 0.40. The printed
values are `..._PRINTED`, and passing them reproduces the print.
Lefebvre & Ballal's Eq. 8.41, `Nu_x = 0.069 [Re_s x/s]^0.7`, is the printed
0.10 times Pr 0.69, so it inherits the misprint. The Cranfield SME4
report's 0.090 is Seban's curve and data times 0.69. See
`validation/cooling/extractions/seban_1960_slot_film_h.md`.

**Validity box** (Seban's heat-transfer runs; not stated as limits):

- x/s about 4-300;
- mass-velocity ratio m 1.06-8.91;
- Re_s about 1.4e3-8e3;
- air into air at a small temperature ratio.

Where the form stops, h hands back to the duct's own h:

- **Past the far end** (x/s 300; for the wall jet x/s 50, where Ballal &
  Lefebvre switch to Colburn), the span is mixed by length:
  `h = f_x h_film + (1 - f_x) h_duct`.
- **Below the hand-back m**: m 1.0 for Seban, who found plain flat-plate h
  below m = 1 (his Fig. 8; his one identified run near 1, m 1.10, already
  sits 25 % below curve B); m 1.3 for the wall jet. The step is a C1
  smoothstep over 0.1 in m. That band is numerical, not physics.
- **Elsewhere outside the box** (x/s < 4, Re_s, m above 8.91) the form is
  applied and flagged (`x_in_box`, `Re_in_box`, `M_in_box`).

`h_film` is the span MEAN of h(x) over the wall's span: the curve is
integrated, not a point value. `x_start == x_end` is the point value; the
point at the lip is refused.

**The wall-jet option** (`slot_film_*_ballal_lefebvre_1973`) is the book's
Eq. 8.43/8.44 (SME4 Eq. 6), `Nu_x = 0.10 Re_s^0.8 (x/s)^0.44` with k of the
slot fluid. It reproduces the book's 1926 W/m^2K from the book's own k_a and
Re_s. It sits 20-67 % ABOVE the Cranfield rig's own data at low Re (SME4
Fig. 2), and only Pai & Whitelaw lie on it. Seban's curves were fitted for m
up to 8.9 and his Fig. 12 is the high-injection data, so Seban is the
recommended form for all m > 1. Where the film is cooler than the wall
(radiation-heated liners), a high h is not conservative.

**Scored** (`python/tests/test_slot_film_htc_validation_962.py`):

- *Fidelity*, on Seban's own Fig. 12 points (38, of which 5 are clean;
  overplotted, so this is an upper bound): bias -1.0 %, RMS 11 %. The
  printed 0.10 is 25 % low on the same points.
- *Accuracy, cross-source*, on SME4 Fig. 1 (Hartnett et al. and Pai &
  Whitelaw, 61 points, m < 1.3): curve B at Pr 0.69 gives bias -2.5 % and
  RMS 18 %. These data lie 32 % above 0.069.

Matching a specific rig is the user's job.

#### Film a downstream wall sees (#386)

```python
f = cb.film_span_average_colban_2011(n_rows=2, row_pitch_over_D=10.0,
        x_start_over_D=0.0, x_end_over_D=50.0, M=1.0,
        area_ratio=3.5, P_over_D=6.0, t_over_P=0.5)
f = cb.film_span_average_baldauf_2002(2, 10.0, 0.0, 50.0, M=1.0,
        density_ratio=1.8, alpha_deg=30.0, s_over_D=3.0, Tu=0.05)
f = cb.slot_film_span_average_marek_tacina_1975(0.0, 30.0, M=1.0, Tu=0.14)  # x in slot heights
f.eta, f.deta_dM, f.deta_dP, f.extrapolated
```

Rows Sellers-superposed, then averaged over `[x_start, x_end]` measured
from the LAST row (in D); `x_start == x_end` is the point value. The span
integral is fixed composite Gauss-Legendre in `ln(1 + x/D)`, so the
derivatives are exact for it (`deta_dP` is 0 for Colban). This is what
`ThermalWall(film=...)` reads.

```python
g = cb.film_span_average_baldauf_2002_groups([2, 3, 2], 12.0, 0.0, 40.0,
        M=[0.8, 1.3, 2.0], density_ratio=[1.4, 1.6, 1.9],
        alpha_deg=30.0, s_over_D=7.0, Tu=0.05)
g.eta, g.weight, g.deta_dM, g.d_weight_dM    # per group (upstream first)
```

Rows in groups (an effusion liner's panels, #695), each at its own M,
density ratio and coolant temperature, layered row by row: past row j the
wall sits at `(1 - eta_j) T_before + eta_j T_c,j`, so
`T_aw = T_ref (1 - eta) + sum_g weight[g] T_c,g` with `sum(weight) = eta`
(Sellers when every T_c is the same). One group is
`film_span_average_baldauf_2002` exactly.

#### Effusion plate internal heat transfer

```python
nu = cb.effusion_internal_nusselt(Re=3000.0, Pr=0.727,
                                  X_over_L=15.24 / 6.3, L_over_D=6.3 / 3.27)

# the two terms separately, if you want to see which dominates
nu_a = cb.effusion_approach_nusselt(3000.0, 0.727, 15.24 / 6.3)
nu_t = cb.effusion_throat_nusselt(3000.0, 0.727, 6.3 / 3.27)
r    = cb.mills_entry_length_factor(6.3 / 3.27)     # entry-length factor
```

`Re` is on the hole diameter and `Nu` on the **hole internal area**. For a
plate-area coefficient multiply by `k/D` and then by `A_h/A` with
`A_h = pi D L` and `A = X**2 - pi D**2/4` -- a factor of 3.46 for Andrews'
plate C, so not optional.

The approach term leads at low `Re` (exponent 0.476) and the throat
overtakes it (0.8); the crossover falls inside a typical effusion plate's
own operating range, which is why both are needed.

#### Turbulent flat plate (the unblown reference of full-coverage data)

```python
nu_x, dnu_dre = cb.flat_plate_turbulent_nusselt_and_derivative(Re_x=4e5, Pr=0.71)
fp = cb.flat_plate_turbulent_panel_nusselt(Re_L=3e4, Pr=0.71, x_start_over_L=1.2)
fp.Nu, fp.dNu_dRe_L, fp.Re_x_start, fp.Re_x_end, fp.extrapolated
nu, dnu = cb.flat_plate_turbulent_panel_nusselt_and_derivative(3e4, 0.71, 1.2)
h_bar = fp.Nu * k / L
```

NASA CR-2648's St_o, printed under every Stanford full-coverage table as the
reference of St/St_o: `St_o Pr^0.4 = 0.0295 Re_x^-0.2`, so
`Nu_x = 0.0295 Re_x^0.8 Pr^0.6`. The panel form integrates it over
`[x0, x0 + L]`: `Nu_L = (0.0295/0.8) Pr^0.6 Re_L^0.8 [(1+s)^0.8 - s^0.8]`,
`s = x0/L`, no unheated starting length. `extrapolated` flags a Re_x range
leaving the data's 1.5e5-5e6 (air, low Tu); never refused.

##### On the element

`EffusionPlateElement.internal_heat_transfer(state_in)` applies the above to
the panel's own geometry and appears in `diagnostics()` after a solve:

```python
diag = result["__element_diag__"]["panel"]
diag["h_plate_area"]   # W/m^2 K on the approach area, Andrews' definition
diag["h_hole_area"]    # W/m^2 K on the hole internal area
diag["approach_to_hole_area_ratio"]  # approach / hole internal, 3.46 for Andrews' plate C
diag["Nu_approach"], diag["Nu_throat"]
```

**Both coefficients are reported because they differ by
`approach_to_hole_area_ratio`** (named `area_ratio` before #773), and
using one where the other is meant is a factor-of-three error. For a panel
energy balance: `Q = h_plate_area * area_approach_total * dT`.

The external film is **not** included here, and `overall_effectiveness()`
below explains why that is a measured decision rather than a gap.

##### The plate's own wall in a network (#471)

The plate IS the wall the coolant passes through, so it owns it: no
`ThermalWall` connects to it, and its diagnostics carry the wall solution
(`wall_heat_transfer`). The node it DISCHARGES into decides the gas side:

| Discharge into | Gas side |
|---|---|
| `MomentumChamberNode` (flow) | by default (`gas_baseline="chamber"`) the chamber's own surface correlation is the unblown `h_gas_unblown` (so the chamber's surface model IS the baseline); with `gas_baseline="flat_plate"` it is a turbulent boundary layer restarted `gas_start_length` upstream of the panel (see below); which baseline that is, and where, is reported beside it as `h_gas_unblown_model` (e.g. `SmoothModel(gnielinski)`, fully developed) and `h_gas_unblown_Re` / `_Pr` / `_Nu` (#476); its velocity gives `blowing_ratio`; `gas_augmentation` (the caller's, 1.0) scales h. `T_gas` is the approaching MAIN stream when the chamber has one -- its own state is the mixed outlet, diluted by the effused coolant |
| a plenum (state, no flow) | the imposed `gas_heat_flux` [W/m^2], 0 by default (adiabatic) |

**No coolant flow (#998)** (below `cb.MIXER_HEAT_MDOT_FLOOR`, the network's
no-through-flow floor): `Re_hole` is 0 and the coolant-side `Nu`/`h` are NaN.
A gas-side chamber below the same floor has no gas flow.

- **Hot gas, no coolant** (a chamber with gas flow): the adiabatic limit,
  q = 0 and both faces at the gas side's `T_adiabatic_wall`, no film.
- **No coolant and no gas-side flow** (a plenum, or a chamber with none): no
  cooling path. The wall temperatures are NaN, `no_cooling_path` is 1, and
  the solver refuses the root (`__success__` False, INCONSISTENT) with
  `EffusionPlateElement 'id': no cooling path: wall T undetermined (heats
  until failure)` in `__message__`, the plates listed in
  `__wall_undetermined__`. Radiation, backside convection without hole
  through-flow and conduction to the mounts are not modelled; the comment at
  `wall_heat_transfer` says how each could be, and #898 is the path for a
  locally starved plate.
- **Temperature limit, always reported.** `T_wall_limit` [K] is the hot end
  of the plate material's #882 envelope (its solidus) when `material=` names
  one from `cb.list_materials()`, else `cb.GENERIC_METAL_SOLIDUS` (1533.15
  K). `wall_above_limit` is 1 when a face is hotter, 0 when not, NaN when the
  wall temperature is undefined; the solve lists such plates in
  `__wall_above_limit__` (reported only). `material` sets the limit only: the
  conduction stays `wall_conductivity`, a constant k. `EffusionLiner` takes
  `material` for its panels and reports `no_cooling_path`, `T_wall_limit` and
  `wall_above_limit` over them.

```python
panel = EffusionPlateElement("eff", "coolant", "liner", hole_diameter=3.27e-3,
                             wall_thickness=6.3e-3, pitch=15.24e-3,
                             panel_length=0.152, panel_width=0.152,
                             wall_conductivity=20.0,      # one layer, t/k
                             gas_augmentation=1.0,        # the caller's
                             gas_film="none",             # or "baldauf_sellers"
                             gas_heat_flux=0.0,           # plenum discharge only
                             material=None,               # T_wall_limit only (#998)
                             gas_baseline="chamber",      # or "flat_plate"
                             gas_start_length=None,       # [m], flat_plate only
                             cd_inclination="none",       # or "idelchik_3_2"
                             radiation=None)              # or GasRadiation(beam_length=...)
diag["T_wall_hot"], diag["T_wall_cold"], diag["q_wall"], diag["eta_overall"]
```

- **Energy-neutral for the network.** The wall heat leaves the gas and
  returns with the effusing coolant into the same node, so it is an output,
  not a source.
- **Film off by default.** `gas_film="baldauf_sellers"` adds C++'s
  `effusion_panel_film_effectiveness` (Baldauf per row, Sellers over the
  panel's rows; needs `panel_length`), flagged `film_extrapolated` outside
  Baldauf's envelope. Scored on Andrews 88-GT-290 the data refuse it offered
  alone; the missing physics is the gas-side augmentation.
- **Plenum-fed.** McGreehan-Schotsch's U1/Vi stays 0 and the coolant side is
  Andrews' still plenum; `coolant_crossflow_ignored` flags a supply node a
  channel runs through. A channel-fed liner segment is a separate element.
- With `wall_conductivity` large the wall reduces exactly to
  `overall_effectiveness` below.
- **Gas-side baseline, opt-in flat plate.** `gas_baseline="flat_plate"` with
  `gas_start_length` [m] >= 0 (where the hot boundary layer restarts -- a
  bleed, a starter-film lip, a liner step -- to the panel's leading edge)
  replaces the chamber's coefficient by `cb.flat_plate_turbulent_panel_nusselt`
  over x in [x0, x0 + panel_length] at the chamber's static state and
  velocity; `T_gas`, `U_gas` and the jet ratios stay the chamber's.
  `h_gas_unblown_model` names it, `h_gas_unblown_x_start` / `_x_end` [m] and
  `_Re_x_start` / `_Re_x_end` give the range, and
  `h_gas_unblown_extrapolated` is 1 where Re_x leaves the Stanford data's
  1.5e5-5e6 (evaluated anyway). Needs `panel_length` and a chamber discharge
  (refused on a plenum); the chamber forbids `gas_start_length`. This is the
  baseline published full-coverage St/St_o ratios use, so pick by which one
  an augmentation you carry over was quoted against. On Behrendt 2008's
  reacting rig (`validation/cooling/effusion_network_runner.py`, accuracy,
  cross-source) it moves the mean over-prediction from +22% to +13% at
  `gas_augmentation` 1.0; nothing here recommends an augmentation.
- **Gas radiation on the hot face, opt-in (#1071).** `radiation=GasRadiation(
  beam_length=...)` (the `ThermalWall` face spec; a `CasingRadiation` is
  refused) adds Lefebvre & Ballal Eq. 8.7 with Eq. 8.11's emissivity, solved
  with the wall in C++'s `wall_radiation_sides`. The gas radiates at `T_gas`
  (the film lowers the convective driving temperature only), the chamber's
  static P and that stream's composition, unless the spec fixes them.
  Reported: `q_rad_hot` [W/m^2] (included in `q_wall`) and `eps_g`; 0 without
  coolant (the adiabatic plate sits at `T_gas`), NaN on a plenum discharge,
  whose imposed `gas_heat_flux` stays the whole gas side. Off, the wall solve
  is unchanged. On Behrendt 2008 (accuracy, cross-source, rig beam length and
  a 0.9 wall emissivity) it moves the flat-plate mean over-prediction from
  +12.5% to +8.8%.
- **Hole-angle Cd term, opt-in.** By default the angle enters the Cd only
  through the drilled length `t/sin(angle_deg)`. `cd_inclination="idelchik_3_2"`
  (`EffusionPlateElement` and `FilmCoolingElement`, `correlation="IdelchikThick"`
  only -- any other correlation, `"fixed"` included, raises) adds Idelchik
  diagram 3-2's inclined-inlet loss to the thick hole's inlet coefficient,
  reported as `Cd_inclination_dzeta` (0.0 when off); outside Re > 1e4 and
  l/d >= 2 `Cd_in_range` reads 0. Cross-source accuracy on Hay 1983 and
  Gritsch 1998 (`validation/cooling/hole_cd_inclination_runner.py`): see the
  runner; off is the default and bit-identical.

##### Duct-fed: the effusion liner (#471)

A liner's coolant usually runs along a backside duct and bleeds through the
wall, so the holes see a SUPPLY-side crossflow. `EffusionLiner` builds that
configuration from N stations:

```python
from combaero.network import EffusionLiner

liner = EffusionLiner("ln", n_segments=4, length=0.152, width=0.152, duct_height=0.03,
                      hole_diameter=3.27e-3, wall_thickness=6.3e-3,
                      pitch_x=15.24e-3, pitch_y=15.24e-3,
                      radiation=None)   # or GasRadiation(...): every panel's hot face
liner.add_to(net, coolant_in="cin", coolant_out="cout", gas="liner_chamber")
liner.summarize(res["__element_diag__"])   # bleed, U1/Vi, Cd, P_backside, T_wall per station
```

- **Stations are Bassett bleeds** (`CrossflowSegmentElement`, kappa 0.75):
  the backside static pressure rises along the duct as it slows, so the
  downstream panels see more drive.
- **Each panel is duct-fed**: `EffusionPlateElement(crossflow_segments=(prev,
  next), crossflow_area=A)` takes McGreehan-Schotsch's supply-side
  `U1/Vi`, with `U1` the duct's mean velocity at the station and `Vi` the
  isentropic jet velocity from the station's STATIC state (C++
  `crossflow_velocity_ratio`). Only `McGreehanSchotsch` has a crossflow term,
  so a duct-fed panel refuses the others.
- **What sets U1/Vi** is mainly the duct's own pressure drop against the
  holes' drive, roughly `sqrt(dp_duct/dp_hole)`, not the duct size: a
  pressure-driven duct carries more flow when made bigger. Flags:
  `crossflow_cd_beyond_8pct` (> 0.2) and `crossflow_cd_degraded` (> 0.35),
  the Rohde-scored bounds.
- **Flagged, provisional:** the coolant-side heat transfer is still Andrews'
  plenum-fed correlation (`coolant_ht_plenum_assumed`), and the stations
  bleed through holes far smaller than the duct, beyond Bassett's measured
  psi = 1-3 (`station_psi_beyond_bassett`).
- **A film provider (#695).** `add_to` registers the panels' film on the
  network under the liner's id, so a wall downstream on the gas path can
  name it: `ThermalWall(..., film="ln")`. Each panel is a group of
  `round(panel_length / pitch_x)` Baldauf rows at its own M, DR and coolant
  temperature (`film_span_average_baldauf_2002_groups`).
  `film_direction="with_coolant"` (default) puts station 1 upstream on the gas
  side; `"against_coolant"` reverses the panels. The undiluted gas (#694)
  excludes every panel's stream when all discharge into one node.

##### Overall effectiveness, as an output

```python
out = panel.overall_effectiveness(state_in, h_gas_unblown=250.0,
                                  T_gas=1600.0, U_gas=80.0,
                                  gas_augmentation=1.0, eta_film=0.0)
out["eta_overall"]       # (Tg - Tw)/(Tg - Tc)
out["T_wall"]            # K
out["T_adiabatic_wall"]  # K, what eta_film sets
out["resistance_ratio"]  # h_internal / h_gas
out["velocity_ratio"]    # jet / mainstream -- reported, never applied
out["blowing_ratio"], out["momentum_flux_ratio"]
```

```
h_gas = h_gas_unblown * gas_augmentation
eta   = (h_i + h_gas * eta_film) / (h_i + h_gas)
```

An overall-effectiveness **correlation** must never be the closure here
-- it already contains the internal convection this computes, and the two
would double-count. So `eta` is an output and the gas side is the
caller's input.

**A film is TWO numbers, and the signature says so.** The standard form
is `q = h_f (T_aw - T_w)` with `T_aw = Tg - eta_f (Tg - Tc)`, so a film
owes both a driving temperature (`eta_film`) and a conductance
(`gas_augmentation = h_f/h_0`). `film_effectiveness_baldauf_2002` gives
only the first -- an adiabatic wall passes no heat, so it measures no
coefficient. **Passing `eta_film` while leaving `gas_augmentation` at 1.0
over-predicts**; both default to their no-film values so that omitting
the pair is consistent rather than half-right.

**Mind what `h_gas_unblown` is measured against.** An augmentation ratio
is meaningless without its baseline. For Andrews 88-GT-290's rig the
required `gas_augmentation` (plate B, with Baldauf's film) is 2.7-4.7
against a fully developed Dittus-Boelter and 1.5-2.7 against the same duct
with a thermal-entry correction -- a factor of 1.76 from that choice alone.
Both plates, both baselines, with and without the film:
`validation/cooling/extractions/andrews_effusion_overall_eta.md` (#476). Published
film-cooling ratios are usually referenced to a flat-plate turbulent
boundary layer at the same x, a third baseline again.

**`gas_augmentation` and `internal_Nu_multiplier` are yours, never
fitted here**, matching `Nu_multiplier` on `ConvectiveSurface`. At their
defaults the closure scores +4.0% on Andrews' effusion plate C and +24.6%
on plate B; setting `gas_augmentation = 2.12` closes plate B at G = 0.6
exactly. Note that `internal_Nu_multiplier = 0.471` also closes it -- by
halving the coolant side, which is the wrong direction, since that
correlation runs 10.4% LOW against Andrews' Fig. 8. Reaching the target
is not the same as being the right dial.

**How good `eta` is depends on a regime this cannot predict.** Andrews'
two plates need gas-side coefficients differing by 1.7x, and that ratio
survives any choice of film model. Pass `U_gas` and read the reported
ratios; if `velocity_ratio` exceeds about 1, treat `eta` as an upper
bound unless `gas_augmentation` already accounts for the jets. No
threshold is applied, because the data supports none -- see
`validation/cooling/extractions/andrews_effusion_overall_eta.md`.

#### Multi-row film superposition

```python
per_row = [cb.film_effectiveness_baldauf_2002(x, M, P, 30.0, 3.0, 0.015)
           for x in row_distances]

eta = cb.film_superposition_sellers(per_row)                  # Gao Eq. (1)
eta = cb.film_superposition_corrected(per_row, alphas)        # Gao Eq. (7)
```

`alphas` has one fewer entry than `per_row`: `alpha_j` is the fraction of the
film's temperature deficit surviving between row `j` and `j+1`. All ones is
plain Sellers.

Sellers **overestimates, and worsens as rows accumulate** -- which is why
`alpha` exists and why effusion cannot reuse a few-row film model unchanged.
`cb.mainstream_temperature_correction(r, a, b)` gives Gao's published form
for it. `a` and `b` are required rather than defaulted even though the paper *does*
print them (12 and 0.9465, section 4.3.1). Since `a r/(a r + 1) >= 0`, that
pair gives `alpha >= 0.9465` for **every** `r`, so it cannot reach the
0.69-0.85 per-row damping a tighter-pitched plate needs -- an argument that
holds however `r` is scaled. Eq. (5) carries no streamwise-spacing term,
and that is the variable Murray shows to dominate.

##### Compressibility: the expansion factor

The paper's mass flow is `W = Cd * Y * A * sqrt(2 rho_t1 (P_t1 - P_s2))`. `Y`
carries compressibility in an otherwise incompressible orifice equation:

```python
S = P_s2 / P_t1
cb.mcgreehan_schotsch_1988_expansion_factor(cd=0.80, S=0.7, gamma=1.4)  # 0.912
```

**Only for the incompressible formulation.** `regime='compressible'` already
solves the isentropic nozzle exactly (`nozzle_flow`, choked branch included),
and Eq. (5) reproduces that solve to 0.008% -- applying `Y` on top of it
corrects for compressibility twice, worth over 8% at `S = 0.7`.

The paper carries two forms because an orifice is not a nozzle, and blends
between them on `Cd`:

```python
cb.mcgreehan_schotsch_1988_expansion_orifice(S=0.7, gamma=1.4)  # 0.912, Eq. (4)
cb.mcgreehan_schotsch_1988_expansion_nozzle(S=0.7, gamma=1.4)   # 0.824, Eq. (5)
```

A sharp hole (`Cd < 0.82`) gets the orifice form, one rounded enough to behave
like a nozzle (`Cd > 0.94`) gets the nozzle form. Two departures from the
printed equations, both deliberate and both recorded as decisions D11/D12 in
the extraction:

- **The blend weight saturates smoothly.** The paper prints `X = 8.333(Cd -
  0.82)` with no upper bound, which reaches 1.45 for a well-rounded long hole
  and extrapolates past a nozzle. A hard clamp would give an exactly-zero
  `dY/dCd` outside the blend and a discontinuous jump of 0.734 at each knee.
  Pass `eps=0` for the paper's exact hard clamp; the default costs 0.47% in
  `Y` and keeps the derivative continuous and non-zero.
- **`Y` saturates at the critical pressure ratio.** Below `S*` the isentropic
  form predicts *decreasing* flow. The paper has no choked branch.

Provenance, the checks behind every constant, and the two errata found in
the paper are in
`validation/cooling/extractions/orifice_discharge_coefficient.md`.

### Geometry and Flow Analysis
```python
# Area calculations
area = cb.orifice_area(d=0.01)
area = cb.orifice_area_from_beta(beta=0.5, D=0.02)

# Cd-K conversions
K = cb.orifice_K_from_Cd(Cd=0.65)
Cd = cb.orifice_Cd_from_K(K=0.5)

# Flow state analysis
state = cb.orifice_flow_state(P1=2e5, P2=1e5, rho=1.2, mu=1.8e-5, d=0.01, D=0.02)
Re_d = cb.orifice_Re_d_from_mdot(mdot=0.1, d=0.01, mu=1.8e-5)
```

### Advanced Flow Functions
```python
# Thermodynamic orifice flow
sol = cb.orifice_flow_thermo(T=300, P=2e5, X=air, m_dot=0.1, area=1e-4, Cd=0.65)

# Impedance with flow effects
# (normalised Z/(rho c) referred to the panel velocity; inertia uses t + 0.85 d)
Z = cb.orifice_impedance_with_flow(freq=500.0, u_bias=2.0, u_grazing=10.0, d_orifice=1e-3,
                                   l_orifice=1e-3, porosity=0.05, Cd=0.7, rho=1.2, c=340.0)

```

### Pressure and Flow Calculations
```python
# Pressure drop
dP = cb.orifice_dP(mdot=0.1, area=1e-4, Cd=0.65, rho=1.2)
dP = cb.orifice_dP_Cd(mdot=0.1, area=1e-4, Cd=0.65, rho=1.2)

# Mass flow
mdot = cb.orifice_mdot(P1=2e5, P2=1e5, area=1e-4, Cd=0.65, rho=1.2)
mdot = cb.orifice_mdot_Cd(P1=2e5, P2=1e5, area=1e-4, Cd=0.65, rho=1.2)

# Velocity and heat transfer
v = cb.orifice_velocity(mdot=0.1, area=1e-4, rho=1.2)
v = cb.orifice_velocity_from_mdot(mdot=0.1, area=1e-4, rho=1.2)
Q = cb.orifice_Q(mdot=0.1, h1=3e5, h2=2.8e5)
```

---

## Heat Transfer

### Correlations
```python
Nu = cb.nusselt_gnielinski(Re=1e5, Pr=0.7)
h = cb.htc_from_nusselt(Nu, k=0.026, L=0.05)

# Multi-layer wall temperature profile
temps, q = cb.wall_temperature_profile(T_hot=1200, T_cold=300, h_hot=200, h_cold=20, t_over_k=[0.01/50, 0.05/0.5])
```

#### Liner radiation kernels (#793)

Lefebvre & Ballal (2010), *Gas Turbine Combustion* 3rd ed., Ch. 8. Each
returns its value with analytic partials (the solver's (f, J) rule).

```python
# Luminous gas emissivity, Eq. 8.11 (Reeves Eq. 8.9 with L = 1). P in Pa.
e = cb.gas_emissivity_luminous(T_g=2280.0, P=3.04e6, q=0.0588, l_b=0.0792, L=1.7)
e.eps, e.deps_dT, e.deps_dP, e.deps_dq, e.deps_dl_b, e.deps_dL   # 0.609, ...

L, dL_dH = cb.luminosity_factor_hydrogen(13.9)   # Eq. 8.15, 336 / H^2 (H mass %)
l_b = cb.mean_beam_length(V, A)                   # Eq. 8.10, 3.4 V/A
l_b = cb.mean_beam_length(V, A, cb.BEAM_LENGTH_COEFF_HOWELL)   # 3.6 V/A

# Gas to liner, Eq. 8.7, partials at fixed eps_g
r = cb.liner_gas_radiation(T_g=2280.0, T_w=1640.0, eps_g=e.eps, eps_w=0.7)
r.R1, r.dR1_dT_g, r.dR1_dT_w, r.dR1_deps_g, r.dR1_deps_w

# Liner to casing, Eq. 8.17 (view factor 1), per unit liner area
R2, dR2_dTw, dR2_dTc = cb.liner_casing_radiation(T_w=1603.0, T_c=880.0,
                                                  eps_w=0.7, eps_c=0.4, area_ratio=0.7)

# Fuel/air ratio from the elements (C + H over the rest), and dq/dY
q, dq_dY = cb.elemental_fuel_air_ratio(Y)
```

`cb.STEFAN_BOLTZMANN`, `cb.REEVES_EMISSIVITY_COEFF` (290),
`cb.LUMINOSITY_HYDROGEN_COEFF` (336), `cb.BEAM_LENGTH_COEFF_LEFEBVRE` (3.4)
and `cb.BEAM_LENGTH_COEFF_HOWELL` (3.6) are the constants. The luminosity
alternatives, Eqs. 8.12-8.14 in the fuel C/H ratio (Lefebvre & Herbert
7.53 (C/H - 5.5)^0.84, its later 3 (C/H - 5.2)^0.75, Kretschmer & Odgers
0.0691 (C/H - 1.82)^2.71), are not shipped: they disagree with each other
and with Eq. 8.15, so pass `L` directly to use one. On a ThermalWall, see
"Radiation (liner)" below.

### Network Heat Transfer

#### ConvectiveSurface
Defines convective heat transfer properties on network elements. Supports different channel models.

```python
import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, PinFinModel, RibbedModel, SmoothModel

# Smooth channel with Gnielinski correlation (default)
surface = ConvectiveSurface(
    area=np.pi * 0.04 * 1.0,  # Surface area [m²]
    model=SmoothModel(correlation="gnielinski")
)

# Ribbed walls (see "Ribbed Channels" below)
ribbed = ConvectiveSurface(
    area=2.5,
    model=RibbedModel(e_D=0.06, p_e=10.0, alpha_deg=90.0, W_H=1.0, n_ribbed_walls=2),
)

# Staggered pin-fin array (see "Pin-Fin Channels" below); area is the
# endwall's BASE (planform) area
pins = ConvectiveSurface(
    area=0.01,
    model=PinFinModel(pin_diameter=0.005, S_D=2.5, X_D=2.5, H_D=1.0, N_rows=10),
)
```

#### ThermalWall
Thermal coupling between two network elements through a shared wall.

**What it represents.** One shared wall *patch* with two sides: element A's
fluid on one face, element B's on the other, the heat crossing the metal
straight through (1-D, normal to the wall, through its layers).

- **Area:** the shared part of each element's perimeter, not its whole wetted
  perimeter -- `contact_area`, or `area_a` / `area_b` when the faces differ
  (#887, below).
- **Several walls per element:** one per face it shares (a rectangular duct
  has four faces; a round one may share part of its perimeter). A face with no
  wall is adiabatic.
- **h:** each side uses its element's own h (the channel correlation's
  perimeter-mean, with its `Nu_multiplier`), applied to this patch; there is
  no face-by-face h.
- **Not modelled: conduction along the metal.** Walls of elements in series
  (a duct and the next one downstream) exchange heat only through the fluid.
  That axial path conducts on the metal's cross-section (thickness x width)
  over the distance between segments, `G ~ k t w / dx` -- usually
  negligible, but not for thick metal, short segments or a steep gradient
  (#898).

```python
from combaero.network import ThermalWall, WallLayer

# Couple hot and cold channels through a wall
wall = ThermalWall(
    id="coupling_wall",
    element_a="hot_channel",      # First element ID
    element_b="cold_channel",     # Second element ID
    layers=[WallLayer(thickness=0.002, conductivity=25.0)],
    contact_area=None,            # Optional: override surface area
    orientation="aligned",        # both from_nodes at the same end of the wall
)

# Add to network
network.add_wall(wall)

# A wall whose two faces differ (#887): a hole in a thick wall, a finned face,
# a curved wall. Each side's convective resistance is on its own face; the
# layers conduct on wall_area's area: "log_mean" (default; the area runs
# linearly through the stack and each layer takes the log mean of its faces,
# exact for a thick cylinder, layer by layer), "slab_a" or "slab_b".
# 1/UA = 1/(h_a A_a) + sum t/(k A_m) + R_fouling/A_a + 1/(h_b A_b).
# Both areas or neither, and not with contact_area.
tube = ThermalWall(
    id="tube_wall", element_a="gas", element_b="coolant",
    layers=[WallLayer(thickness=0.015, conductivity=20.0)],
    area_a=2 * 3.14159 * 0.025 * 0.3, area_b=2 * 3.14159 * 0.010 * 0.3,
)
```

**The heat is the two streams' exchange (#627)**: the epsilon-NTU integral
along the wall, `Q = eps(NTU, C_r) C_min (T_aw,a - T_aw,b)` at each stream's
inlet, with `NTU = UA / C_min` and `C = |m| cp` (cp the secant over the inlet
span, `cb.cp_secant_mass_and_jacobian`, so `C_min dT` is the exact maximum). Behind a combustor the side's composition moves with the unknowns too:
`ChannelResult.dh_dY` and `dT_aw_dY` (smooth channels) and
`cb.cp_secant_mass_dY` carry h, T_aw and C to the relay (#717).
It never exceeds what the weaker stream can take, and a stream that stops
takes none. A node side is a well-mixed reservoir (`C = inf`). The linear
`U A dT` it replaced is its small-NTU limit; it moved up to 2.3x the heat at
NTU ~ 1 and drove a stalled stream below 0 K.

- **The arrangement is not declared, it follows.** `orientation` is geometry:
  `"aligned"` (the default) puts both elements' `from_node` at the same end of
  the wall, `"opposed"` at opposite ends. With the flow signs it decides
  parallel or counter flow, so a flow that reverses flips the arrangement --
  continuously, since it passes zero flow, where both forms meet.
- **Model limit:** one lumped exchanger at constant (secant) cp. Over a wide
  span at high NTU the parallel outlets can cross by a little (20 K in a
  1300 K span at NTU 5.4); segment the wall for more.
- **Reported:** `{wall}.UA` [W/K], `{wall}.NTU`, `{wall}.effectiveness`
  (`Q / (C_min dT)`), `{wall}.counterflow` (1.0 or 0.0).

Each side's heat goes into that element's `to_node` (or the node itself when
the wall names a node). A **boundary** cannot take it -- its state is fixed --
so heat a wall puts there leaves the network with the stream crossing that
boundary, and the result reports it:

```python
res["coupling_wall.Q"]              # W, side a -> side b
res["coupling_wall.Q_to_boundary"]  # W, the part that landed on boundaries (0 if none)
res["outlet.Q_wall_out"]            # W, wall heat leaving through boundary "outlet"
```

With it the network's energy balance closes:
`sum_out m h(T_upstream) + sum_boundaries Q_wall_out = sum_in m h(T_in)`
(`python/tests/test_energy_conservation_cooling.py` checks it to 1e-9 for
every cooling configuration).

- **Reversed flow (#481).** Streams follow the SIGN of an element's flow: a
  2-port element with `m_dot < 0` delivers into its `from_node`, at its
  `to_node`'s state, and a wall puts that side's heat into the node the flow
  goes into. Junction and tee elements keep their declared directions, which
  their own models enforce.
- **k(T) layers** (a `material` from `cb.list_materials()`) take their
  conductivity at the wall temperatures of the same evaluation, iterated to
  a fixed point -- not from the previous one.
- **h <= 0** from a correlation past its range is floored smoothly at
  `cb.WALL_HTC_KNEE` (1e-2 W/(m^2 K)); the wall never refuses it.
- With `thermal_coupling_enabled = False` the walls report nothing.

Every converged solve in the test suite is also checked for this balance
(`python/tests/conftest.py`, `_closure_check.py`).

##### A side with no stream: adiabatic or heat loss (#970)

Side B may have no stream: `element_b=None` with `side_b`. The layers run
from side A's stream to the free side; the free side is B only.

```python
from combaero.network import HeatLoss

# An insulated backside: no heat crosses; the wall sits at side A's T_aw
# (with radiation_a=GasRadiation, at the face's radiation balance).
net.add_wall(ThermalWall("w1", "hot_duct", None, side_b="adiabatic",
                         layers=[WallLayer(0.0016, 16.0)], contact_area=0.0012))

# A lossy backside: a fixed T_ext behind h [W/(m^2 K)]
net.add_wall(ThermalWall("w2", "hot_duct", None,
                         side_b=HeatLoss(T_ext=300.0, h=10.0),
                         layers=[WallLayer(0.0016, 16.0)], contact_area=0.0012))
```

- **Adiabatic:** Q = 0 and the network solves exactly as without the wall
  (the duct is not wall-heated). Reported: `{wall}.Q` = 0, `T_hot`,
  `T_interface` (the whole stack at the face temperature), `T_aw_a`, `h_a`,
  the film and radiation keys of side A. The face temperature with
  radiation is `_solver_tools.wall_adiabatic_face` (C++):
  `h (T_aw - T_s) + h_r (T_src - T_s) = 0`.
- **HeatLoss(T_ext, h):** side A's stream exchanges with the reservoir at
  `T_ext` through `1/(h_a A_a) + sum t/(k A_m) + R_fouling/A_a + 1/(h A_b)`:
  the epsilon-NTU integral against a reservoir,
  `Q = C_a (1 - exp(-UA/C_a)) (T_aw_a - T_ext)` (C++
  `wall_exchange_and_jacobian` with `C_b = inf`), bounded and none into a
  stream that does not flow; analytic Jacobian. The heat leaves the network:
  `{wall}.Q_loss` [W] (> 0 out), booked by the ledger (`__ledger__["Q_loss"]`)
  and in the wall's second-law check as `Q / T_ext` into the reservoir.
  `h` is the user's input (radiation from that face linearised into it).
- **Not with:** `radiation_b`, a film on side B, and for an adiabatic wall
  `radiation_a=CasingRadiation` (its heat would leave to a casing no stream
  carries; use HeatLoss). Per-side areas, k(T) layers, `GasRadiation` on
  side A and a film on side A work as on a two-stream wall.
- Not in the GUI yet: exporting such a wall to the GUI refuses it.

##### Radiation (liner) (#793)

Off by default. A face can carry radiation: `radiation_a` / `radiation_b`.

```python
from combaero.network import CasingRadiation, GasRadiation

net.add_wall(ThermalWall(
    id="liner", element_a="primary_zone", element_b="annulus",
    layers=[WallLayer(0.0012, 26.0)],
    # gas -> liner, Lefebvre Eq. 8.7 with eps_g from Eq. 8.11
    radiation_a=GasRadiation(beam_length=0.08, emissivity_wall=0.7,
                             luminosity=cb.luminosity_factor_hydrogen(13.9)[0]),
    # liner -> casing, Eq. 8.17
    radiation_b=CasingRadiation(area_ratio=0.7, emissivity_wall=0.7,
                                emissivity_casing=0.4),
))
```

**What it represents.** The Lefebvre-Herbert liner method as Lefebvre &
Ballal print it (Ch. 8):

- `GasRadiation`: `R1 = 0.5 sigma (1 + eps_w) eps_g Tg^1.5 (Tg^2.5 - Tw^2.5)`
  into the face, with `eps_g = 1 - exp[-290 P L (q l_b)^0.5 Tg^-1.5]` (P in
  kPa). `Tg` is the side's stream Tt, `P` its static pressure and `q` its
  fuel/air ratio from its elements (`cb.elemental_fuel_air_ratio`), unless
  `T_gas`, `pressure` or `fuel_air_ratio` is given.
- `CasingRadiation`: `R2 = sigma eps_w eps_c / (eps_c + eps_w (1 - eps_c)
  A_w/A_c) (Tw^4 - Tc^4)` out of the face, to a casing at the side's stream
  Tt (the book's T3) unless `T_casing` is given. The heat goes into that
  side's stream: in steady state a casing that loses nothing outward hands
  what it absorbs to the annulus air.
- Each face's radiation adds to its convection in the local wall balance,
  `h (T_aw - T_s) + R1 = conduction = h (T_s - T_aw) + R2`, solved for the
  two surface temperatures at the streams' inlet states. It is then folded
  into each side's conductance and driving temperature, `H = h + h_r`,
  `T* = (h T_aw + h_r T_src) / H` (`h_r` is the radiation written exactly as
  a conductance), and the streams' epsilon-NTU exchange runs on those.
  Between two reservoirs `Q / A` is the book's balance exactly; with finite
  capacity rates the linearised radiation integrates along the wall like the
  convection.
- k(T) layers and the radiation share one fixed point; the Jacobian is
  analytic through both (the gas Tt, P and composition included).
- **Reported:** `{wall}.q_rad_a` / `q_rad_b` [W/m^2 of that face],
  `Q_rad_a` / `Q_rad_b` [W], `T_surface_a` / `T_surface_b` [K], and
  `eps_g_a` / `eps_g_b` for a gas face. `{wall}.UA` includes the radiative
  conductances.

**What it does not represent.**

- **Gray and lumped.** One mean gas temperature, one beam length, one
  emissivity per surface. No spectral bands, no zoning (Hottel), no
  wall-to-wall view factors inside the liner, no axial radiation between
  segments. `alpha_g / eps_g = (Tg/Tw)^1.5` (Eq. 8.6) and the wall's
  `0.5 (1 + eps_w)` (Eq. 8.5) are the book's approximations.
- **A kerosine-product law.** Eq. 8.9 is Reeves' (1956) fit for an
  industrial gas turbine burning kerosine, written in `q`, with no stated
  range. It is not a composition-based emissivity and does not carry over to
  natural gas, hydrogen or ammonia products (no CO2, or much more H2O per
  `q`). A pure-hydrogen flame has `q` from its H alone.
- **No soot model.** Luminosity enters only through `L`. The default `L = 1`
  is the non-luminous form; `cb.luminosity_factor_hydrogen(H)` gives
  Lefebvre's engine fit (Eq. 8.15). The C/H forms (Eqs. 8.12-8.14) disagree
  with it and with each other.
- **Your geometry and your emissivities.** `beam_length` (3.4 V/A of the
  gas volume the face sees, `cb.mean_beam_length`; 0.6-0.9 D_L for tubular
  liners, 1.0/1.2 D_L for annular inner/outer liners per the book; the
  textbook value is 3.6 V/A), `area_ratio`, `emissivity_wall` (0.7 Nimonic,
  0.8 stainless steel; a TBC surface differs) and `emissivity_casing` are
  inputs the user owns. The defaults are the book's example values, not a
  recommendation.

Fidelity (printed quantities): the book's three chained examples -- uncooled
1640/1603 K, film-cooled at x/s 18 1283/1265 K, 0.5 mm TBC 1304/1236/1220 K
-- are reproduced from their printed inputs through `ThermalWall`
(`python/tests/test_liner_radiation_793.py`; the uncooled case comes out
1.6 K above print, the book having rounded `eps_g` to 0.61 before solving).
As a network with combaero's own convection (Gnielinski, h 523 against the
book's 565 W/m^2K) the uncooled liner is 1633/1597 K: a method comparison,
not fidelity. The film and TBC cases are networks too, with a slot's
film-side h (#962, `SlotFilmElement(film_htc=...)`; see the slot film
below). They are a method comparison, not fidelity
(`python/tests/test_slot_film_htc_validation_962.py`):

| Case | Book (printed) | Seban, the default | Eq. 8.44, named option |
|---|---|---|---|
| Film-cooled, x/s 18 | 1283/1265 K | 1285/1267 K | 1269/1252 K |
| With a 0.5 mm TBC | 1304/1236/1220 K | 1308/1240/1224 K | 1288/1224/1209 K |

- The book uses the wall jet (Eq. 8.44) at m = 2.54, with h = 1926.
- Seban's h on combaero's slot air is 1632, which puts the wall about 10 K
  hotter at the book's T_aw.
- Eq. 8.44 on combaero's properties gives 2068 (k_a 0.060 against the
  book's 0.0553), about 4 K cooler.
- combaero's Marek-Tacina eta (0.801 against the book's 0.789) lowers T_aw
  by 18 K.

##### A film on one side (#386)

Film holes are a `FilmCoolingElement`: an orifice of `n_rows * holes_per_row`
holes from the coolant supply into the gas, and a film PROVIDER. It owns no
wall and is not a wall side -- the surface it cools lies downstream on the
gas path. The wall on that gas duct names it:

```python
from combaero.network import FilmCoolingElement

net.add_element(FilmCoolingElement(
    "film", "coolant", "gas", hole_diameter=1e-3,   # METERING D for shaped
    holes_per_row=30, lateral_pitch=6e-3, wall_thickness=2e-3,
    n_rows=2, row_pitch=0.01, angle_deg=30.0,
    hole_type="shaped", area_ratio=3.5, breakout_width=3e-3,  # Colban
    # hole_type="cylindrical", turbulence_intensity=0.05      # Baldauf
    # cd_inclination="idelchik_3_2"  # opt-in hole-angle Cd term (IdelchikThick)
))
net.add_wall(ThermalWall(
    id="liner", element_a="hot_duct", element_b="backside",
    layers=[WallLayer(0.002, 25.0)],
    film="film", film_side="a",          # the gas duct downstream of the holes
    film_start=0.0, film_end=0.05,       # [m] downstream of the LAST row
))
```

- **A film slot is a provider too** (`SlotFilmElement`, #714): a continuous
  slot of `slot_height` across `span`, as in a combustor liner's cooling
  ring. Film: Marek & Tacina (1975), `eta = 1/(1 + C_m x/(M S))` with
  `C_m = 0.25 Tu/100`; the wall's span is measured from the slot lip.
  `turbulence_intensity` is the axial free-stream Tu as a FRACTION, like
  every element's (#773), default 0.14 = 14%
  (`SLOT_FILM_DEFAULT_TU`: D-7958's 52% blockage case, where Eq. 5 fits that
  report best, inside D-6360's combustor 10-20%) -- a working example; a
  liner's Tu comes from what is upstream of it. `C_m_multiplier` is the
  user's rig-matching factor, never scored. Flow: an orifice of the slot's
  area, `lip_thickness` long, its discharge correlation fed the hydraulic
  diameter `2SW/(S+W)` (give a measured Cd with `correlation="fixed"`).

  ```python
  net.add_element(SlotFilmElement("ring", "annulus", "gas",
      slot_height=2e-3, span=0.25, lip_thickness=1e-3, turbulence_intensity=0.14))
  ```
- **A slot's film-side h** (opt-in, #962): `SlotFilmElement(...,
  film_htc="seban1960")` gives the filmed wall side the film's own h, the
  mean of Seban's curve over the wall's span (`film_end > 0`), in place of
  the duct's. The film's h comes from the slot flow and the slot fluid's
  properties at its supply plenum's state. It hands back to the duct's h
  past x/s 300 and below m = 1 (kernel section above).
  `"ballal_lefebvre1973"` is the book's m > 1.3 wall jet (Eq. 8.44), a named
  option. `"none"` (the default) keeps the duct's h. The duct's
  `Nu_multiplier` does not act on the film's h. The Jacobian reaches the slot
  flow, its coolant temperature (through mu, cp, k) and the gas flow
  (through m). The slot fluid's composition is not chained, as the density
  ratio's is not. **Reported:** `{wall}.film_h` (the side's h),
  `film_h_duct`, `film_h_film` (the film's span mean over its in-box
  part), `film_h_fraction_in_box`, `film_h_weight_M`, `film_h_Re_s`,
  `film_h_branch` (`cb.SLOT_FILM_BRANCH_*`: 0 curve B, 1 the B/C blend, 2
  curve C, 3 wall jet, -1 none), `film_h_in_box` and its parts
  `film_h_x_in_box` / `film_h_Re_in_box` / `film_h_M_in_box`.

  ```python
  net.add_element(SlotFilmElement("ring", "annulus", "gas",
      slot_height=2e-3, span=0.25, lip_thickness=1e-3, film_htc="seban1960"))
  ```
- **An effusion panel is a provider too** (`film="<panel id>"`, needs
  `panel_length`): its rows, Baldauf each `pitch_x` apart and Sellers over
  them, read past its LAST row. Naming it selects that film whatever the
  panel's own `gas_film` says. A full-coverage panel usually runs outside
  Baldauf's s/D and M envelope (`film_extrapolated`), and at high M Baldauf
  has the jets lifted off, so the film it carries downstream is small.
- **An effusion liner is a provider too** (`film="<liner id>"`, #695): its
  panels, one group each, see the liner's section above. Each panel's flow
  and station temperature reach the Jacobian; `{wall}.T_coolant_film` is the
  film's weighted coolant temperature and `{wall}.film_blowing_ratio` /
  `film_density_ratio` the flow-weighted means.
- **By default only the driving temperature changes.** The covered side's
  T_aw becomes `T_aw - eta_bar (T_aw - T_coolant)`; its h, stream capacity
  and heat booking stay the duct's. A slot can also carry its own h
  (below). `eta_bar` is `cb.film_span_average_*` over the
  wall's span: the hole type's kernel per row, Sellers over the rows.
- **Blowing ratio** is the metering mass flux over the duct's gas mass flux;
  **density ratio** is `(W_c/W_g)(T_g/T_c)` at a common static pressure.
  No forward coolant flow, no film.
- **Jacobian:** T_aw,f reaches the provider's flow, its supply node's
  temperature and the undiluted gas temperature (below), all relayed
  analytically.
- **Reported:** `{wall}.eta_film`, `{wall}.film_blowing_ratio`,
  `{wall}.film_density_ratio`, `{wall}.film_extrapolated` (the kernel's
  envelope), `{wall}.T_coolant_film`, `{wall}.T_aw_unfilmed` (the T_aw the
  film is read against), `{wall}.film_gas_undiluted` (1 when that is the
  undiluted gas), and on every wall `{wall}.T_aw_a` / `{wall}.T_aw_b`.
- **Not scored as a coupled wall.** Each kernel is checked on its own
  author's data; no source we hold measures a heated wall downstream of a
  film row.
- **The gas is the gas before the coolant mixed in (#694).** Baldauf and
  Colban define `eta = (T_inf - T_aw)/(T_inf - T_c)` with T_inf the
  mainstream before injection. When the covered side is fed by the node the
  holes discharge into (a plenum or momentum chamber), T_inf is that node's
  other inflows mixed without the holes' stream (its energy boundaries
  included, wall heat into it not): the side's T_aw shifts by
  `T_inf - T_mix` and DR takes T_inf. Otherwise (a combustor, a boundary, a
  side further downstream) the side's own state stands in, flagged by
  `film_gas_undiluted = 0`. Mixing-out further downstream in a confined
  duct is outside both correlations: where `eta` falls below the coolant's
  mixed-in share, a filmed wall reads HOTTER than the same wall unfilmed
  (which assumes the coolant fully mixed) -- the conservative side.

#### Element Integration
Elements with convective surfaces support heat transfer calculations.

```python
from combaero.network import ChannelElement

# Channel with convective heat transfer
channel = ChannelElement(
    id="heated_channel",
    from_node="inlet",
    to_node="outlet",
    diameter=0.04,
    length=1.0,
    roughness=0.0,
    surface=surface  # ConvectiveSurface instance
)

# Get heat transfer coefficient and temperature
h, T = channel.htc_and_T(state)  # Returns (htc [W/(m²·K)], T [K])
```

#### Thermal Coupling Toggle
Enable/disable thermal coupling globally for debugging or fast solves.

```python
network = FlowNetwork()
network.thermal_coupling_enabled = True   # Enable wall heat transfer
network.thermal_coupling_enabled = False  # Disable (faster, no heat transfer)
```

---

## Acoustics

```python
# Acoustic properties bundle
props = cb.acoustic_properties(f=1000, rho=1.2, c=340, p_rms=1.0)
print(f"SPL: {props.spl} dB")

# Cavity resonators
f_helm = cb.helmholtz_frequency(V=0.001, A_neck=1e-4, L_neck=0.01, c=340)

# Duct modes
tube = cb.Tube(L=1.0, D=0.1)
modes = cb.tube_axial_modes(tube, c=340, upstream=cb.BoundaryCondition.Closed,
                            downstream=cb.BoundaryCondition.Open)
```

---

## Programmatic Network Design

For automated design and agents, `combaero.network` provides a pure-Python interface to build and solve networks without the GUI.

```python
from combaero.network import FlowNetwork, NetworkSolver, OrificeElement, PressureBoundary, PlenumNode
import combaero.species as species

# 1. Initialize Network
net = FlowNetwork()

# 2. Define Boundaries & Nodes
net.add_node(PressureBoundary("inlet", Pt=5e5, Tt=400, Y=species.dry_air_mass()))
plenum = PlenumNode("plenum")       # its P (= Pt) is solved; T comes from the inflows
plenum.initial_guess = {"plenum.P": 4e5, "plenum.Pt": 4e5}   # optional seed
net.add_node(plenum)
net.add_node(PressureBoundary("outlet", Pt=1e5, Tt=300))

# 3. Add Elements
net.add_element(OrificeElement("feed", "inlet", "plenum", diameter=0.01, Cd=0.65, correlation="fixed"))
net.add_element(OrificeElement("exhaust", "plenum", "outlet", diameter=0.015, Cd=0.62, correlation="fixed"))

# 4. Solve
solver = NetworkSolver(net)
results = solver.solve(method="hybr")

# 5. Access Results: a flat dict keyed "<id>.<quantity>"
print(f"Plenum Pressure: {results['plenum.P']} Pa")
print(f"System Mass Flow: {results['feed.m_dot']} kg/s")
```

> [!TIP]
> **Agent Hint**: Use `net.to_dict()` to get a JSON-serializable representation of the network, which can be saved or sent to the GUI.

---

## Network Solver

### What a solve reports

`NetworkSolver.solve` returns the solved unknowns plus a set of `__dunder__`
bookkeeping keys. Five of them say what happened:

| key | type | meaning |
|---|---|---|
| `__success__` | `bool` | converged **and** consistent. Unchanged; keep using it for "can I trust this result". |
| `__converged__` | `bool` | the root finder reached the residual tolerance, whatever the consistency checks then said |
| `__consistent__` | `bool \| None` | the elements' own physical-consistency verdict. **`None` means not checked** -- either the solve never converged, or no element in this network has a verifier. It is not a pass. |
| `__inconsistent_elements__` | `list[str]` | ids of the elements that rejected the solution |
| `__P_nonpositive__` | `list[str]` | nodes whose pressure (`P` or `Pt`) is at or below zero at the root: not a gas state, though an element's continuation can close its row there (a Fanno kernel below the sonic static). Reported INCONSISTENT, checked before any diagnostics evaluate the state (#810). |
| `__T_clamped__` | `list[str]` | nodes with flow whose temperature is at or below 0 K at the solution: more heat out than the stream holds as a gas (the thermo is an ideal gas down to 0 K). Reported INCONSISTENT (`__success__` False). The name is kept from the [200, 5000] K node clamp this replaced (#626). |
| `__wall_undetermined__` | `list[str]` | effusion plates with no cooling path at the root (#998): no coolant flow through the holes and no gas flow on the other side, so the wall temperature is undetermined (`no_cooling_path` 1, wall keys NaN). Reported INCONSISTENT, the plate named in `__message__`. |
| `__wall_above_limit__` | `list[str]` | effusion plates with a face above `T_wall_limit`, their material's envelope hot end (#998; `wall_above_limit` 1). Reported only. |
| `__fraction_unanchored__` | `list[str]` | fraction losses without `m_ref` on their zero-flow blend in a driven flow component, not a dead end (#1004; `fraction_unanchored` 1): the loss law is not met at that root. Reported only, never a demotion. Checked once on the outermost solve's result, after every retry has restored the boundaries; on a failed solve the losses flagged at the best iterate are named in `__message__`. |
| `__heat_faded__` | `list[str]` | heated nodes at a root inside the heat fade below `MIXER_HEAT_MDOT_FLOOR`, taking up part of their heat: neither flowing nor stagnant, a root the fade makes. Reported INCONSISTENT. |
| `__supersonic_chambers__` | `list[str]` | compressible momentum chambers at a root past the sonic limit of their area, REFUSED: the chamber's Pt is below the sonic Pt* of its mass flux at its area, so no state of that area carries the flow (since #788 the closure has one root per Pt, so this is the only way past the sonic point; the supersonic root Newton once landed on is gone). Not a choked-bounded junction (below): the reason is in `__message__` (`'e1' (LosslessConnectionElement) feeds it, not a choked channel`, `row(s) [...] read its static pressure`, a channel wider than the chamber). Usually the element downstream is no restriction for that flow, or the chamber's area is too small for it. Reported INCONSISTENT. |
| `__choked_bounded_junctions__` | `list[str]` | momentum chambers ACCEPTED past their sonic limit (#788): a choked compressible channel of no larger area discharges into the chamber, which is its under-expanded exit -- the "choked bounded" exit a pressure boundary takes (#636) -- and no row but the chamber's own reads its static pressure. The root is the one a `PlenumNode` in its place gives; the expansion loss Pt* - Pt sits in the channel's control volume (its sigma). The chamber reports as that plenum: `P` = `Pt`, `mach` 0, `T` = `Tt`, with `choked_bounded` 1, `P_sonic` / `Pt_sonic` (the sonic state of its flux at its area) and `expansion_loss` = Pt* - Pt; the channel reports its own choked exit (`mach_out` 1). `__x_solution__` keeps the chamber's P unknown, on the closure's continuation. |
| `__outcome__` | `SolveOutcome` | why it ended, in one machine-readable value |
| `__worst_residuals__` | `list[dict]` | the rows carrying the residual, largest first |
| `__second_law__` | `dict[str, dict]` | per element: `verdict` (`"pass"`, `"fail"`, `"not_checked"` with a `reason`), `sigma` [W/K], `sigma_norm` [-] and its band `sigma_norm_low`/`_high`; chain members carry their `group`. A control volume with a stream on the incompressible closure widens the band by that closure's own resolution, `sigma_norm_resolution` = sum of m x_i (x_max + x_i/2) over the inflow, x = (Pt - P)/P per stream (order M^4; absent on compressible ports): each port carries its energy at its own constant density, so a split lossless in its own model sits on the bound. A converged root with any element at `"fail"` (`sigma_norm_high < -SIGMA_TOL = -1e-8`: it creates flow work on its most favourable heat path) is reported INCONSISTENT (`__success__` False, #598), with the element named in `__message__`. |
| `__second_law_failures__` | `list[str]` | ids whose verdict is `"fail"` |
| `__second_law_walls__` | `dict[str, dict]` | per ThermalWall, the wall as its own control volume (#551): `verdict`, `sigma_wall` [W/K] = the entropy the heat carries into side B's stream minus what it took from side A's, each integrated along that stream's own temperature path; `sigma_wall_norm` (over the mean \|S\| of the two sides) and its band `_low`/`_high`. Negative when heat runs from the stream that is colder on average. A necessary condition only (end temperatures can cross while the path means do not). **Report only**: it never changes `__success__`. |
| `__second_law_wall_failures__` | `list[str]` | wall ids whose verdict is `"fail"` |
| `__incompressible_validity__` | `dict[str, dict]` | per element on an incompressible closure: `verdict` (`"ok"`, `"warn"`, `"not_checked"` with a `reason`), its highest internal `mach` and the `limit`. **Report only** (#588). |
| `__incompressible_warnings__` | `list[str]` | ids whose verdict is `"warn"` |
| `__ledger__` | `dict` | mass and energy over the network's control volume: `streams` (one row per crossing, `m_dot_in` > 0 entering), `m_dot_in/out`, `H_in/out`, `Q_given`, `Q_wall_out`, `Q_withheld`, `Q_loss` (heat to a wall's HeatLoss side, #970), `mass_closure(_rel)`, `energy_closure(_rel)`. Empty when the root finder did not converge. A converged root with either relative closure above `combaero.network.second_law.LEDGER_CLOSURE_TOL = 1e-6` is reported INCONSISTENT (`__success__` False, #552); the suite's accepted solves close to 4.2e-8 (mass) and 1.1e-7 (energy). |
| `__caught_exceptions__` | `list[dict]` | every exception the solve caught instead of raising, at most 50: `phase` (`start`: an attempt's starting point; `probe` / `probe_jacobian`: a trial point the root finder rejected; `diagnostics`: a node's or element's diagnostics on a FAILED solve's final state, which then reports none for it -- on a successful solve that error raises (#810); a `wall_free_seed:` prefix for the wall-free seed solve), `where` (the node or element being evaluated, or the non-finite Jacobian rows and columns), `type`, `message`, `frames` (innermost three). Empty on a clean solve. |
| `__caught_exception_count__` | `int` | all caught exceptions, including those beyond the 50 kept |

**Newton-homotopy last retry (#642), opt-in.** `solver.homotopy_retry = True` (the GUI backend, `NetworkRunner` and the stress harness set it) gives a cold solve that did not converge one more try:
- **The homotopy:** H(x, λ) = F(x) − (1 − λ) F(x0), from the cold guess x0, which solves it at λ = 0. λ rises adaptively, each step warm-started from the last, and the λ = 1 step is a plain solve with every check.
- **Budget and limits:** at most `solver.homotopy_budget` seconds within the caller's timeout (default `NEWTON_HOMOTOPY_BUDGET` = 5 s; `None` gives it all the time the timeout leaves, which the GUI uses so the user's max solve time is the only limit). It does not run for a root the checks rejected (`__converged__` true), and it stops when no step from the guess succeeds, the case of a singular Jacobian at x0.
- **What it gives:** 3 of the 23 stress random draws that failed every other start now solve within their 20 s.
- **Why it is off by default:** every non-converging solve pays its budget. A junction scorecard with expected failures went from 8 s to 50 s.

**Boundary-condition ramp retry (#642), opt-in.** `solver.bc_ramp_retry = True` (the GUI backend, `NetworkRunner` and the stress harness set it) runs after the homotopy, in the time the timeout leaves:
- **The ramp:** at λ every `PressureBoundary` sits at Pm + λ (Pt − Pm), with Pm the mean of their Pt, and every `MassFlowBoundary` flow and user `EnergyBoundary` heat at λ times its own. The gentle end (`BC_RAMP_LAM0` = 0.05, halved down to `BC_RAMP_LAM_MIN` until it solves cold) is solved with the usual retries. λ then rises, each solve warm-started from the last root, and the λ = 1 solve is a plain solve with every check. The boundaries are restored whatever happens.
- **Starts cold:** a cold start drops the Fanno kernels' warm-start guesses (`_last_G`, `_last_M_exit`) left by the failed attempts. One kept a draw's gentle end from solving.
- **Budget:** with the ramp on, an uncapped homotopy (`homotopy_budget = None`) takes half of the time left instead of all of it, and a capped one at most `NEWTON_HOMOTOPY_BUDGET_BEFORE_RAMP` = 3 s (#1018). Its stress rescues (random 60, 91, 136) take under 2 s; the 5 s it spent on the draws it did not rescue was the time the ramp then lacked (random 119). Not 2 s: each homotopy step's timeout is what the cap leaves, and a solve runs its LM fallback only with more than 2 s left.
- **What it gives:** it rescues 3 of the 21 stress random draws the homotopy left unsolved (31, 34, 35), in 13-46 s. In the 20 s harness that is random 34, plus legacy_tee 15. Draws that still fail now use more of their budget. Since the compressible seed (#971), random 34 converges cold; random 35 is the draw the ramp still rescues.

**Resistance-network cold start, on by default.** `solver.resistance_seed = True` replaces the cold start's node pressures and element flows with the root of the network as quadratic resistances:
- **The network:** each element is dPt = K m|m|, with K its incompressible start-guess estimate (`_seed_resistance`) at the reference density. Pressure boundaries are fixed, a `MassFlowBoundary` injects its m_dot (a sink draws it), and every other node conserves mass.
- **A channel's K:** friction at f = 0.02, plus one exit head only where its residual loses it: the declared exit is a pressure boundary that takes the flow at its static pressure (#360). Between plenums a channel carries friction alone. Charging the head on every segment started the ten 1.22 cm segments of the D-7156 film rig 14 kPa high against a 13 Pa root, and the cold solve stalled (#971). A `CrossflowSegmentElement` always charges the head: its static stations lose it by construction. Before the topology is resolved the exit is unknown and the head is charged as before. A network built from the estimate keeps it: the stress splitter sizes its supply from it.
- **Compressible elements:** the incompressible K misses choking. A short duct between plenums passes the choked flux on a small Pt difference. So once the quadratic root is found, any compressible channel or orifice whose seed flux reaches `SEED_COMPRESSIBLE_FLUX_GATE` (0.3) of the sonic flux at its upstream Pt has its K re-fitted. The new K is the drop over the flow the residual's own C++ kernel passes at the seed's pressures (`_seed_flow`: Fanno for a channel, the orifice's compressible model). The root is then re-solved, for at most `SEED_COMPRESSIBLE_PASSES` (4) passes, stopping once ln K moves less than `SEED_COMPRESSIBLE_RTOL`. Below the gate the start is unchanged. The kernels run at the reference temperature: the flow component's mean boundary Tt (#899, #990).
- **Why it always has a start:** the network is monotone, so it has exactly one solution, found by damped fixed-point iteration on the conductances (`SEED_MAX_ITER`, `SEED_RELAX`, `SEED_MDOT_FLOOR`; a lossless connection gets `SEED_LOSSLESS_FRACTION` of the smallest K).
- **The gate:** it applies when every element offers an estimate (orifice, channel, pressure loss with a shipped correlation, area change, lossless connection) and every node is a plenum, momentum chamber, combustor, wall or boundary. Any other network (junctions, ejectors, plates, ...) keeps the propagated start unchanged.
- **Why:** the propagated start sets each node's pressure from its neighbours one element at a time. In a large network of resistances decades apart, it leaves a duct carrying MPa at the guessed flow. 48 of 160 stress random draws stalled there even made incompressible and adiabatic.
- **What it gives** (stress harness, 20 s, `--jobs 4`):
  - `random`: 67 -> 113 of 160 converged, in about half the total time;
  - `hotcold` / `hotcold_mild`: about 4x faster;
  - `impingement`: 4 of 80 draws, slower from it, cross the 20 s budget (all converge at 60 s).

  It is a start guess only; it never enters a residual.

**Starting-point safety net.** The root finder always rejected a trial point whose evaluation raised (a penalty residual); the starting point was the one evaluation it did not cover, so a property refusing a physical-but-extreme state there aborted `solve()`, retries included. Now a `RuntimeError`, `ValueError` or `ArithmeticError` at an attempt's starting point ends that attempt with outcome `ERROR` (a `RuntimeWarning` names where), and the retries try their other starts. A starting point whose residual is non-finite counts as one that could not be evaluated (`FloatingPointError`, naming the rows' owners); it used to be reported converged, since a NaN norm passes no comparison. Code bugs (`TypeError`, `AttributeError`, `KeyError`, `IndexError`) always raise. Strict mode raises as before: `solver.raise_internal_errors = True`, or the environment variable `COMBAERO_SOLVER_STRICT=1` when the attribute is left `None` (the test suite runs strict, `python/tests/conftest.py`). For a stress or discovery run, either go strict to stop at the first finding, or leave it lenient and collect `__caught_exceptions__` across the solves: it also lists the rejected probes, which strict mode does not change.

**Second law per element (#545).** Each element's control volume must not
destroy entropy: `sigma = sum_out m s - sum_in m s - integral dQ/T >= 0`, with
the streams at their stagnation states (`cb.s_mass(Tt, X, Pt)`, ideal mixing
included) and the heat taken at the fluid's static temperature
(`cb.heat_entropy_flux`), which makes `sigma` the internal generation. For an
adiabatic two-port `sigma_norm = sigma / (R m)` is exactly `-ln(Pt_out/Pt_in)`:
total pressure does not rise. Heating at Mach `M` costs total pressure and
cooling may raise it (Rayleigh); the heat term allows exactly that. Inside an
element only its end Mach numbers are known, so the heat is evaluated along a
linear Mach path and bracketed by constant-Mach paths at either end; a
`"fail"` holds even for the most favourable one (`sigma_norm_high <
-SIGMA_TOL`, `SIGMA_TOL = 1e-8`). A chain of crossflow segments is one control
volume (a segment's momentum spans half of each station). The legacy
`TeeJunctionElement` is not checked. Each element's `sigma_norm` also appears
in its diagnostics.
A failing verdict demotes the solve to INCONSISTENT (#598). An element whose
own closure explains the failure says why through `second_law_hint(sol) -> str`,
appended to `__message__`: `IdealMomentumTeeElement` names a `loss_factor`
below ~M^2/4 at its port Mach, where its low-Mach loss form creates flow work.

**Incompressible validity (#588).** Each element running an incompressible
closure is checked against one Mach limit, `solver.incompressible_mach_limit`
(default `combaero.network.compressibility.INCOMPRESSIBLE_MACH_LIMIT = 0.2`),
at its own highest internal Mach: an orifice's jet (`mach_throat`), a
channel's higher end (`mach_in`, `mach_out`); an element declares which in
`validity_mach_keys`. A momentum chamber on its incompressible closure is
judged at its own `mach` and listed by its node id; the GUI attributes an
inserted port chamber to its junction. The incompressible closures are the low-Mach limits of
the compressible ones, drifting as Ma^2: measured at the same flow, a channel's
drop is 4.0% off at Ma 0.2 (13% at 0.4), an orifice's 2.9% at a jet Ma of 0.2
(5% at 0.26). The limit is a rule of thumb, not a guarantee; the solution and
`__success__` are unchanged.

```python
solver = NetworkSolver(net)
solver.incompressible_mach_limit = 0.3   # default 0.2
r = solver.solve()
for eid in r["__incompressible_warnings__"]:
    v = r["__incompressible_validity__"][eid]
    print(eid, v["mach"], ">", v["limit"])
```

`__success__` is False in two quite different situations, and reading it alone
cannot tell them apart: Newton never got there, or Newton got there and a
junction rejected the root as unphysical. The second reports a *small*
`__final_norm__` next to `success=False`, which looks like a contradiction
until `__converged__` and `__consistent__` are read separately.

```python
from combaero.network import SolveOutcome

result = solver.solve()
if result["__outcome__"] == SolveOutcome.INCONSISTENT:
    print("converged, then rejected by", result["__inconsistent_elements__"])
elif not result["__success__"]:
    worst = result["__worst_residuals__"][0]
    print(f"{result['__outcome__']}: {worst['name']} carries {worst['residual']:.3e}")
```

`SolveOutcome` is a `StrEnum`, so it compares and serialises as a plain
string: `CONVERGED`, `INCONSISTENT`, `NO_PROGRESS`, `RESIDUAL_TOO_LARGE`,
`TIMEOUT`, `ERROR`, `NOT_CONVERGED`, `NO_UNKNOWNS`.

> [!TIP]
> Prefer `__outcome__` over matching on `__message__`. Part of that text comes
> from SciPy and can be reworded without notice; `__outcome__` is classified
> once, by the code that knows the answer.

Note there is no `__residual_norm__`; the norm is `__final_norm__`.


### Basic Usage
```python
from combaero.network import FlowNetwork, NetworkSolver, OrificeElement

graph = FlowNetwork()
graph.add_element(OrificeElement("ori1", "nodeA", "nodeB", diameter=0.011284))  # Reader-Harris/Gallagher Cd

solver = NetworkSolver(graph)
# init_strategy options: "default", "analytical_pt_prop", "homotopy",
# "continuation", "outletref_warmstart" ("incompressible_warmstart" was
# removed in 0.8.0 and raises ValueError).
# "default" auto-upgrades to "analytical_pt_prop" for networks that
# contain a MultiPortChamberBase; pass another strategy or an
# explicit x0 to opt out.
# "outletref_warmstart" solves the outlet-referenced incompressible
# proxy (densities at the downstream static) and warm-starts the
# compressible solve from it directly -- the auto-retry's seed as a
# primary strategy, for networks known to stall the cold path.
# "default" also starts there (falling back to the cold path, auto-retry
# and all, if it fails) when an element declares prefers_warm_start;
# ImpingementCrossflowElement does: impingement stress family 37/40 ->
# 40/40 converged, 167 s -> 33 s, same roots (#706).
# auto_retry (default True): failed cold solves on compressible
# MPCE networks are retried once from an outlet-referenced
# incompressible warm start (densities at the downstream static);
# on MPCE networks incompressible throughout, from a compressible
# warm start (same root at low Mach, #588; reported as
# "compressible_warmstart").
# The primary attempt gets 40% of `timeout`, the retry the rest.
# A cold solve of a network with walls that still fails is then
# retried from the same network solved WITHOUT its walls -- the
# message says "Converged from the wall-free flow solution ...".
# A network whose prefers_warm_start start failed retries from the
# inlet-referenced seed, not from the outlet-referenced one it just tried
# (#1018).
# An iterate the model refuses (an exception, NaN/inf) is answered with a
# residual ten times the last physical one, so the trust region rejects
# the step and shrinks; a failed solve lists the rejected probes in its
# message. Each evaluation repeats the state propagation until states read
# ahead of its order (wall back-edges, recirculation) settle, so the
# residual is a function of x alone.
# Stall detection ends doomed phases early: when hybr's best |F|
# plateaus far from tolerance it hands over to the LM fallback, and
# when the LM fallback plateaus too the attempt fails fast so the
# auto-retry gets the budget. Far from tolerance a phase also hands
# over when its recent rate projects tolerance beyond 3x the time the
# attempt has left (60 s without a timeout). Cold attempts only -- seeded solves
# (explicit x0 or a warm-start seed) are the last rung and keep
# their full budget. stall_handover=False (default True, #864) turns
# every handover off for the call, its retries and nested solves
# included: each cold rung then runs to its own share of the budget.
# The GUI's "Stall handover" checkbox sets it; the results table keeps
# it with the solver options (solver:stall_handover).
results = solver.solve(method="hybr", init_strategy="homotopy")
```

### Network Nodes
```python
from combaero.network import (
    PressureBoundary, MassFlowBoundary, FuelBoundary, WallNode, PlenumNode,
    CombustorNode, MomentumChamberNode, EnergyBoundary
)

# Boundary conditions
inlet = PressureBoundary("inlet", Pt=2e5, Tt=300, Y=air)
outlet = PressureBoundary("outlet", Pt=1e5, Tt=300, Y=air)
mass_in = MassFlowBoundary("mass_in", m_dot=0.1, Tt=400, Y=air)

# Fuel controller (#892): the fuel flow that puts combustor "burner" at a
# target exit Tt. The flow is a solver unknown; two roots exist (lean, rich).
# branch="lean" is the default, "rich" for RQL or reformers. A root on the other
# side of stoichiometric, a target not above the oxidizer's Tt, or a target
# above the peak is refused with the reason. Reports fuel.m_dot and fuel.phi.
# `combustor` names the controlled node, like a probe: the fuel may reach it
# through any elements (fuel line, premix plenum), as long as all of it does.
# Air through a pressure boundary: a cold solve first holds the fuel at its
# start, solves the flows, re-sizes the fuel on that air, then the controller.
fuel = FuelBoundary("fuel", combustor="burner", Tt_target=1500.0, branch="lean", Tt=300, Y=ch4)

# Closed-end boundary (zero mass flow)
# Typical use: symmetry plane of a ring manifold, or any dead-end branch.
wall = WallNode("wall")

# Internal nodes: their P, T and Y are solved (or derived), not given
junction = PlenumNode("junction")
combustor = CombustorNode("combustor", method="complete")  # or "equilibrium"
momentum = MomentumChamberNode("momentum", area=0.01, regime="compressible")
# A chamber narrower than an adjoining pipe (#1023) is a contraction or
# expansion with no loss: FlowNetwork.validate warns (UserWarning naming the
# chamber, the pipe and both areas; momentum.pipe_area_warnings(net) returns the
# messages, and the GUI's resolve preview shows them), suggesting an explicit
# AreaChangeElement. Never for the DEFAULT_CHAMBER_AREA placeholder, a merge
# chamber's area inherited from an orifice main inlet, or a pipe at the GUI's
# default diameter; equal is to NEIGHBOUR_AREA_RTOL on the equivalent diameter.
combustor.add_energy_boundary(EnergyBoundary("energy", Q=50000))  # heat addition [W]
combustor.add_energy_boundary(EnergyBoundary("loss", fraction=-0.05))  # -5% of the SENSIBLE enthalpy
junction.add_energy_boundary(EnergyBoundary("heater", Tt_target=650.0))  # set outlet Tt [K]: Q solved
```

`EnergyBoundary(id, Q=0.0, fraction=0.0, Tt_target=None)` attaches to a
`PlenumNode`, `MomentumChamberNode` or `CombustorNode`. Sign: > 0 heats the
fluid. Three modes:

| mode | set | the solver |
|---|---|---|
| heat flow | `Q` [W] | derives the outlet temperature |
| fraction | `fraction` [-] | of the sensible enthalpy (below) |
| outlet temperature (#792) | `Tt_target` [K], the node's TOTAL temperature | solves the heat flow |

- **`Tt_target`** makes the heat an unknown, `{id}.Q` [W], closed by one row
  `Tt_node - Tt_target = 0`, scaled as a temperature (row kind `"T"`, see
  `NetworkSolver._build_residual_scales`). Its Jacobian column is C++'s
  `dTt_mix/dQ`, relayed to every row that reads the node's temperature. It is
  refused (`ValueError`) together with a nonzero `Q` or `fraction` on the same
  boundary, as a second `Tt_target` on one node, and on any other node type.
  A fixed-`Q` boundary beside it on the same node is fine: the target's `Q`
  is what the rest leaves. On a `CombustorNode` it is the heat taken out (or
  put in) AFTER the reaction to reach that exit Tt; it does not move the fuel
  flow. A node with no through-flow (below `cb.MIXER_HEAT_MDOT_FLOOR`) has no
  heat that sets its temperature: the solve does not succeed, its message
  names the boundary, `__heat_undetermined__` lists it and its `Q` is
  reported as NaN.
- **Reported heat**, in every mode: `{id}.Q` [W] per boundary (the set `Q`,
  the solved `Q`, or its share of `{node}.Q_fraction`) and `{node}.Q`, their
  sum. The second-law check and the balance ledger (#552) use the solved heat.
  `{node}.Q_withheld` and `{node}.Q_fraction` [W] are always reported on a
  node with an energy boundary, 0.0 when nothing is withheld or no fraction
  heat applied.
- **A result's keys follow the network, never the operating point.** A key
  whose value is undefined at a point (a per-unit-flow quantity at zero flow,
  a merge coefficient on a split, a choking anchor of an unchoked ejector) is
  present and NaN; the GUI's `/solve` writes it as `null`.
  `python/tests/test_result_keys_invariant.py` solves every palette network at
  its own point, at zero, reversed and choked flow and requires the same keys.

- **`fraction`** scales the inflow's sensible enthalpy, `h(T) - h(298.15 K)`
  at its own composition (`cb.SENSIBLE_ENTHALPY_REF_T`); on a combustor, the
  products'. It is not a fraction of absolute enthalpy, which includes
  formation enthalpy (#481). The heat it applied is reported as
  `{node}.Q_fraction` [W].
  It must be above -1, per boundary and summed over a node (#1024): at -1
  it removes all the sensible enthalpy.
- **`Q`** is spread over the node's real inflow. Only a near-stagnant node
  (below `cb.MIXER_HEAT_MDOT_FLOOR` = 1 mg/s) takes up less, fading C1 to none
  at zero flow (`cb.mixer_heat_per_unit_flow(m)` returns the enthalpy rise per
  unit heat and its derivative), so a dead end holds the temperature of what
  reaches it; what it withholds is reported as `{node}.Q_withheld` [W], as is
  wall heat aimed at a `WallNode`.
- **A `MassFlowBoundary` between elements** is an injection: its `m_dot`
  joins the node at its own `Tt` and `Y`.
- **Where the network ends** is a property of the node class, read by graph
  validation, the solver and the GUI export (#545):

  | attribute | True on | meaning |
  |---|---|---|
  | `is_boundary` | `PressureBoundary`, `MassFlowBoundary`, `WallNode` | the network may end here: one connected element, or flow one way only |
  | `is_flow_boundary` | `PressureBoundary`, `MassFlowBoundary` | mass and enthalpy cross the control volume here, at the state the node imposes |
  | `state_at_rest` | `PlenumNode`, `CombustorNode`, `WallNode`, `PressureBoundary`, `MassFlowBoundary` | the node's own definition holds Pt = P, so a `ChannelElement` reports its duct's own state at this end (#888) |

  A `WallNode` is a boundary nothing crosses. An injecting `MassFlowBoundary`
  is still a flow boundary.
- **A discrete loss (`PressureLossElement`) follows the flow** (#550).
  - "In" is the upstream node of the flow as it runs, and the correlation sees `|m_dot|`.
  - Through zero flow the forward and reverse forms blend over `FLOW_REVERSAL_M_SCALE` (1e-4 kg/s).
- **Orifices and incompressible channels read their upstream side along the flow** (#554), so an element's declared direction is only a sign convention:
  - **Orifices:** the kernel is evaluated from the node the flow comes from: its stagnation state, its density, and the velocity of approach of the pipe on its side (`beta`, or `beta_reverse` from the to-node's pipe).
    - The direction comes from the pressures, blended through equal Pt over `ORIFICE_REVERSAL_P_SCALE` (0.1 Pa).
    - An orifice into a dead end (no flow boundary beyond its to-node) is never reversed: every root has zero flow there.
    - **A zero bore is a closed restriction** (#1024). `OrificeElement(diameter=0)` or `area=0`, `EffectiveAreaConnectionElement(effective_area=0)` and `DiameterDischargeCoefficientConnectionElement(diameter=0)` are walls (`element.closed` is True). The row is `m_dot = 0` with the constant Jacobian 1 (row kind "mdot"), and no kernel is evaluated. The diagnostics are `closed = 1`, the pressures either side, and zero flow, area and velocity. The resolve preview reports `closed = 1`, and the start guess seeds the flow at zero. A part of the network that only closed elements join to the rest solves when it holds a `PressureBoundary` of its own. A part they cut off from every `PressureBoundary` is refused at `FlowNetwork.validate()`, naming the closed elements, because nothing sets the pressure of the gas trapped there. The #1004 dead-end test and an orifice's own dead-end test do not search past a closed element. A negative or non-finite bore, area or effective area is refused at construction.
    - **With a pipe upstream (E applied) the orifice drives from that pipe's static state** (#651), the node's (Pt, Tt) at the pipe's mass flux, not from the node's own static pressure: a `PlenumNode` (P = Pt) and a `MomentumChamberNode` then give the same flow. Reported as `P_approach`. The row then solves its own consistency, m* = m_calc(m*, x).
    - **Between pipes the orifice charges its permanent loss** (#645). With a pipe on the side the flow comes from (beta) and a wider one on the side it goes to, a metering or `fixed` orifice's drop is the nodes' stagnation difference, ISO 5167-2's k(beta, C) times the tapping differential, and the kernel expands to the implied tap pressure. Into a node with no pipe beyond it, the jet's head is lost and the drop is the tapping differential. The row's Cd is at D-D/2 taps; `OrificeElement(tap=...)` only picks the `dP_tap` it reports.
    - **The bore against its pipe** (#788, 0.8.0). The pipes that set beta (the single upstream channel forward, the single to-side channel in reverse) are compared with the bore to `ORIFICE_BORE_PIPE_RTOL` (1e-6 relative). Equal on either side is **no orifice**: zero loss, the row is `Pt_in - Pt_out = 0` with Jacobian +1/-1 (row kind "p"), `beta = 1`, and the resolve preview reports `pass_through = 1`; this is how an orifice is switched off without editing the network. A bore larger than either pipe is refused (`ValueError`). Until 0.8.0 a bore at or above the pipe silently became a wall hole (beta = 0, and only in that direction); that fallback is gone. Write a lumped Cd*A restriction as `OrificeElement(..., in_pipe=False)` (GUI: untick "In a pipe"): no pipe diameter, beta = 0 both ways, no recovery. Its metering Cd (Reader-Harris/Gallagher) is then evaluated as for a plenum-fed orifice, at D = 1 m; the old fallback used D = 10 d.
    - **The bore against an adjacent momentum chamber** (#1016). The same rule with areas: the bore's own geometric area (for `EffusionPlateElement`, `FilmCoolingElement`, `SlotFilmElement` and `ImpingementPlateElement` the total hole area; never Cd A) against the area of a `MomentumChamberNode` the element opens into or draws from, compared on the equivalent diameter to `ORIFICE_BORE_PIPE_RTOL`. Larger is refused at `resolve_topology` (`ValueError` naming the element, the chamber and both areas), whether or not `in_pipe`; equal is no orifice (the pass-through row above; the resolve preview's `pass_through` source names the chamber); smaller is unchanged, even with Cd A close to the chamber area. Only a chamber area that states the geometry is compared: a user's area or Dh, or one inherited from a channel or another area-bearing neighbour. Not compared: the `DEFAULT_CHAMBER_AREA` placeholder of a chamber nothing sizes, a merge chamber's area inherited from an orifice as its `main_inlet`, and, with `in_pipe=False`, an area inherited from a channel (that pipe's end; the bore against the pipe is the rule above, which `in_pipe=False` waives). `EffectiveAreaConnectionElement` (a Cd A, no bore) and `DiameterDischargeCoefficientConnectionElement` are not checked.
  - **Incompressible channels:** the friction and exit-head density is the upstream node's along the flow, blended through zero flow over `FLOW_REVERSAL_M_SCALE`. It is taken at the duct's own inlet static state, the node's (Pt, T) at the duct's flux (#651), not at the node's static pressure.
  - **Cold start:** on networks driven by pressure boundaries alone, the cold-start pressure guess of an interior node is clamped to the boundaries' range.
  - A fraction loss (`ConstantFractionLoss`, `LinearThetaFractionLoss`) takes an optional design flow `m_ref` [kg/s]: at and above it the loss is exactly the fraction; below it the loss scales with the dynamic head (`1 - (1 - r^2)^3`, `r = |m_dot|/m_ref`, C2) and vanishes at zero flow.
  - Set `m_ref` wherever the flow may run low or reverse. A constant fraction holds a pressure difference at zero flow, and a reversing flow cannot be carried through it.
  - **A fraction loss without `m_ref` on its zero-flow blend is flagged** (#1004). Where the network gives it a smaller drop than `xi Pt`, its only root is the blend at `|m_dot|` below `FLOW_REVERSAL_M_SCALE`, and the solve converges there. `fraction_unanchored` (every `PressureLossElement`) is 1 there when its flow component is driven and the loss is not a dead end. Driven is read from the inputs, never an iterate: pressure boundaries whose Pt differ by more than 1e-9 of the largest, a `MassFlowBoundary` feeding a nonzero `m_dot`, a `FuelBoundary` whose `Tt_target` is above the coldest other boundary Tt, or a compressor or turbine at a pressure ratio other than 1. A dead end is topology: one side of the loss reaches no pressure boundary, fuel controller or nonzero mass-flow boundary without it. 0 otherwise (any head loss, any loss with `m_ref`, zero drive, a dead end); NaN where its flow is not finite. `PressureLossElement.unanchored_fraction` says whether the loss is one that can be. The GUI's Inspector offers "use the solved flow as m_ref": at `m_dot = m_ref` the loss is exactly `xi`, so the solved point is unchanged.
  - **The four built-in loss classes evaluate in C++** (`cb.combustor_fraction_loss`, `cb.combustor_head_loss`, #718). `corr.evaluate(ctx, Y_in)` returns a `CombustorLossResult` (`xi` and its partials in theta, m_dot, Tt_in, Pt_in and each inlet `Y_k`), and `PressureLossElement` builds its Jacobian from it. Calling the class still returns `(xi, dxi/dtheta)`.
  - **The context is the inlet STAGNATION state** (#782 audit): a `PressureLossContext` carries `Tt_in`/`Pt_in`, and the head loss is Lefebvre's total-total low-Mach pressure-loss factor, `xi = zeta 0.5 rho v^2 / Pt_in` with `rho` at (`Tt_in`, `Pt_in`). Until 0.8 the context's `P_in` was the inlet's STATIC pressure, so the loss came out high at a moving inlet (+8% at M 0.3); at a plenum inlet nothing changes. The 0.7 spellings `ctx.T_in`/`ctx.P_in`, `combustor_head_loss(T_in=, P_in=)` and `CombustorLossResult.dxi_dT_in`/`dxi_dP_in` were removed in 0.9.0; a hand-built context without `Tt_in`/`Pt_in` raises `AttributeError`. A coupled Fanno-Rayleigh / impulse form whose low-Mach limit is this one is #830.
  - **A user-supplied callable** returning `(xi, dxi/dtheta)` still works. Its other dependences are central differences; that interface is the prototype for loss maps (#725).
- **A `MomentumChamberNode` carries its stagnation temperature `Tt`** (#576). Its inflows mix as `(m, Tt, Pt, Y)` streams. The closure solves the static T from `h(T0) = h(T) + v^2/2`, `v = m_dot / (rho(T, P) A)`, and closes `Pt = P exp((s(T0) - s(T)) / R)`.
  - **One root per Pt** (#788). At fixed `(T0, |m|/A)` that `Pt(P)` falls to its minimum `Pt*` at the sonic static pressure `P*` and rises again below it, so one Pt had a subsonic and a supersonic static root, and Newton could land on the supersonic one. Below `P*` the closure is continued by its mirror image about `Pt*`, `2 Pt* - Pt(P)`: monotone in P through the sonic point (C1 there), so the root is subsonic whenever `Pt >= Pt*`. A `Pt < Pt*` roots on the continuation (the result's `bounded`, with `P_sonic` / `Pt_sonic`): accepted as a choked-bounded junction, or refused (`__choked_bounded_junctions__`, `__supersonic_chambers__`). The partials of `Pt*` need no `dP*/dx` (the closure's own slope is zero at `P*`). The incompressible closure has no sonic point and is unchanged.
  - The node's `mach` and its reported properties (`rho`, `a`, `velocity`, `Re`, plus its static `T`) come from that static state.
  - Since 0.8 every node's result key `T`, and every element's `T_in`/`T_out`, is the STATIC temperature, like `P`; the total is `Tt`, `Tt_in`/`Tt_out` (#782, see [Static T in network results](#static-t-in-network-results-08)).
  - `MomentumChamberNode.static_density(state)` returns `(rho, {name: d rho / d name})`. Junction elements (`MultiPortChamberBase.port_densities`) take their port densities from it, so a port's momentum rows and its node's Pt describe the same flow.
  - The chamber now passes up to the true choked flux. Before, T0 read as static capped it at 92 % of choked.
  - The merge-chamber face state (`cb.chamber_merge_face_state`) takes both temperatures as stagnation and solves the face's and the outlet's static states; its sensitivities are analytic (#578).
  - The chamber's heat-transfer correlation, and an effusion plate's gas-side density, use the chamber's static state (#578).
  - Still read as static: an incompressible `ChannelElement` whose upstream node is a chamber (#578). The legacy `TeeJunctionElement` now hands its kernel the static T of a chamber port at the arm's flow (#782); so do `CrossflowSegmentElement`'s station and entry drops and a duct-fed `EffusionPlateElement`'s `U1/Vi`, at the station's mean flow.
- **A `MomentumChamberNode` has a `regime`** (#588), default `DEFAULT_REGIME` like the elements; the GUI gives every chamber (drawn, or inserted at a junction port, an edge or a merge) the network's global regime.
  - `"compressible"`: the closure above.
  - `"incompressible"`: its low-Mach limit. The static state is the stagnation one (`T_static = T0`, `rho = P / (R T0)`) and `Pt = P + (m/A)^2 / (2 rho)`: no sonic limit, so an incompressible network keeps none hidden in its chambers. The merge face takes both densities at the chamber pressure, `P_face = P + m_out^2/(rho A^2) - J/A - m_main^2/(rho_main A^2)`, linear with nothing to choke.
  - Both agree with the compressible forms to O(Ma^2): 0.25% of the head at Ma 0.2, the gap quartering as Ma halves.
  - A junction's reference head (`MultiPortChamberElement`'s constant K) follows its port chamber's regime.
  - A script that validates the compressible closure must ask for it: `MomentumChamberNode(..., regime="compressible")`.

#### Areas against their neighbours (0.11)

One rule, applied wherever an element's area meets a pipe or chamber. Equal means within the 1e-6 band on the equivalent diameter.

| the area is ... | then | where |
|---|---|---|
| equal to its neighbour | a pass-through: zero loss, the row `Pt_in - Pt_out` | orifice bore vs pipe (#788) or chamber (#1016); `AreaChangeElement` F0 == F1 (#1023) |
| larger than a geometry allows | refused, naming both areas | bore above its pipe or chamber; explicit F0/F1, F_C/F_branch against the attached pipe; ejector `mixing_area` above its outlet pipe (#1023) |
| a contraction with no loss | a `UserWarning` from `validate()` | a `MomentumChamberNode` narrower than an adjoining pipe (#1023) |
| zero on a restriction | closed: no flow, `closed = 1` | orifice, effective-area and Cd D connections (#1024) |

Buildable but unusual geometry stays accepted, for example a tee branch wider than its main. A zero area on a duct or node is refused (below).

#### Geometry the network refuses (#1024)

Each refusal names the element or node, the input and the limit, and raises `ValueError` at `FlowNetwork.validate()` unless noted. Every check is a method `_geometry_error()` on its class, which returns the message or "". The GUI's resolve preview reports that message as the item's `refusal` before Solve.

| what | refused |
|---|---|
| `ChannelElement` | a diameter (given or inherited) that is not positive and finite; a given `Dh` that is not positive and finite |
| `MomentumChamberNode`, `CombustorNode` | a given area or `Dh` that is not positive and finite; a given `Dh` above `2 sqrt(A/pi)` |
| `EffusionPlateElement` | `pitch_y <= d`; `pitch_x <= d / sin(angle_deg)`, the hole's footprint on the surface (rows taken inline: the panel gives no arrangement) |
| `FilmCoolingElement` | `lateral_pitch <= d`; with `n_rows > 1`, `row_pitch <= d / sin(angle_deg)` (rows taken inline) |
| `ImpingementPlateElement` | `pi/4 / (xn_d yn_d) >= 1`; `yn_d <= 1`; `xn_d <= 1` inline; staggered (`pattern="staggered"`, the next row half a spanwise pitch over) `min(sqrt(xn_d^2 + (yn_d/2)^2), 2 xn_d) <= 1` |
| `CompressorElement`, `TurbineElement` | a non-finite `pressure_ratio` (at construction) |
| `CombustorNode` | a `method` other than `"complete"` / `"equilibrium"` (at construction and at `validate()`) |
| `EnergyBoundary` | `fraction <= -1`, per boundary (at construction) and summed over a node's boundaries |

- **Holes must leave material between them.** Equal pitches are refused too, because touching holes leave no plate. If a valid staggered film or effusion design is refused, it gets its own issue.
- **The isoperimetric bound** (chambers and combustors). `Dh = 4A/P`, and no closed section of area `A` has a shorter perimeter than the circle, so `Dh <= 2 sqrt(A/pi)`. The tolerance is `ISOPERIMETRIC_DH_RTOL` (1e-3 relative), which absorbs an area typed to three significant figures. A `ChannelElement` is not bounded: a channel that stands for a share of a wider duct wets only its part of the wall, and keeps the duct's Dh, as the C3X vane example's hole segments do.
- **A fraction of -1 removes all the sensible enthalpy.** The node then lands on the 298.15 K reference whatever it is fed. On a combustor the report divided by zero. Below -1 the node goes below the reference temperature, and from a hot inflow below 0 K. To set an outlet temperature, use `Tt_target`.
- **To shut a path**, use a closed restriction (above), not a zero-area duct.

#### Merge chamber: main inlet plus side streams (#471)

A `MomentumChamberNode` with `main_inlet` declared accepts further inflows as
SIDE STREAMS, still with exactly one outlet. The chamber's axis is its outlet
direction; the main stream arrives along it; each side stream brings
`m u_jet cos(theta)` of axial momentum and discharges at the chamber's static
pressure. A constant-area control volume closes it:

    P_face A + m_main u_main + sum(m_s u_jet cos theta) = P A + m_out u_out

```python
liner = MomentumChamberNode("liner", main_inlet="duct")   # area inherited from "duct"
panel = EffusionPlateElement("eff", "coolant", "liner", hole_diameter=0.8e-3,
                             wall_thickness=2e-3, pitch=6e-3, panel_area=0.01,
                             angle_deg=30.0)               # jets at 30 deg to the gas
```

- **The node's `(P, Pt)` is the chamber (outlet) state.** The main-inlet
  element sees the MAIN-FACE state the impulse balance gives
  (`cb.chamber_merge_face_state`): exact for compressible flow, the face
  static pressure on the subsonic root and its Pt from the same
  entropy-based closure as the chamber. With no side streams the face is the
  chamber, so declaring a single-inlet chamber's main inlet changes nothing.
- **Jets are compressible too.** A side stream's impulse is `cb.jet_impulse`:
  the isentropic velocity from its supply (Pt, Tt) to the chamber pressure,
  or, choked, the sonic momentum plus its pressure thrust.
- **No loss coefficient.** The stagnation-pressure loss of mixing follows
  from momentum; transverse momentum is reacted by the walls and its kinetic
  energy dissipated. Normal injection costs the main stream its mixing loss;
  inclined jets push it, ejector-like.
- **Directions live only inside the junction.** Edges carry scalars and a
  channel is 1D along its own axis, so no angle is ever passed on.
- **One area.** The main inlet's port area is the chamber area, inherited
  from the main-inlet element; an explicit area that disagrees is refused.
- **The injector owns its jet.** Orifice-type elements carry
  `injection_angle_deg` (90 by default, no axial momentum);
  `EffusionPlateElement` uses its hole inclination `angle_deg`. The jet
  velocity is the element's own, from its supply's stagnation state. A
  crossflow-segment station reads the same declaration (#813).
- **Refused:** two inflows without `main_inlet`, a `main_inlet` that does not
  flow in, and any split (more than one outlet needs a loss closure: use a
  tee).
- In the GUI, inflows on the chamber's "s" handle are side streams; the one on
  its flow handle is the main inlet.

### Network Elements

**Writing or subclassing an element.** Every element declares the units of its residual rows for the solver's row scaling: `residual_scale_kind = "p"` (Pa) or `"mdot"` (kg/s), a class attribute or a property; a junction (`MultiPortChamberBase`) declares `row_scale_kinds()`, one per row. There is no default. `FlowNetwork.validate()` raises `TypeError` unless the declaration sits on the class that defines `residuals` or a subclass of it, so a subclass that rewrites `residuals` re-declares its units. A wrong kind weights the row 1e5-1e9 wrong in the norm the root finder minimises (RESIDUAL-ROW-SCALING in `NetworkSolver._build_residual_scales`). A strict solve (see the safety net above) also checks every declaration against the Jacobian at the converged solution and raises `RowKindMismatchError` (a `TypeError`, in `combaero.network.solver`) on a misfiled row; `solver.check_row_kinds(x)` returns the per-row report (`owner`, `declared`, `measured`, `verdict` = agree / DISAGREE / unclear / untestable) for any `x`, best a converged one.

```python
from combaero.network import (
    CompressorElement, TurbineElement,
    OrificeElement, ChannelElement, EffectiveAreaConnectionElement,
    LosslessConnectionElement, PortDischargeElement, DiameterDischargeCoefficientConnectionElement,
    TeeJunctionElement, VortexElement, EffusionPlateElement,
)
from combaero.network.mpce_element import ConstantKTeeElement, MultiPortChamberElement

# Flow elements
# Cd is the "fixed" correlation's value only (#872, 0.8.0): with any other
# correlation (default "ReaderHarrisGallagher") the correlation computes it and a
# Cd given too raises ValueError; "fixed" without a Cd raises as well.
# correlation="Stolz" was removed in 0.9.0 (#876) and raises ValueError.
orifice = OrificeElement("orifice", "node1", "node2", diameter=0.011284, regime="compressible")
measured = OrificeElement("measured", "node1", "node2", Cd=0.65, diameter=0.011284, correlation="fixed")
# The bore against the pipe it sits in (#788, 0.8.0): equal (to
# ORIFICE_BORE_PIPE_RTOL = 1e-6 relative) is NO orifice -- zero loss, the row is
# Pt_in - Pt_out = 0 and the flow passes unchanged -- and larger is refused. A
# lumped Cd*A restriction is an orifice with no pipe diameter: in_pipe=False
# (beta = 0 both ways, no recovery, whatever channels adjoin it).
lumped = OrificeElement("lumped", "node1", "node2", Cd=0.8, diameter=0.05, correlation="fixed", in_pipe=False)
# The same against a MomentumChamberNode it opens into or draws from (#1016):
# the bore's area (a plate's total hole area) larger than the chamber's is
# refused, equal is no orifice; a defaulted chamber area is not compared.
# How the compressible regime expands the flow (#854); the incompressible regime
# ignores it. compressible_model=None (default, 0.8.0, #861) selects ONCE at
# build, never per iterate: "jobson1955" for a Cd referenced to incompressible
# flow (ReaderHarrisGallagher, IdelchikThick/Beveled/Rounded, Lichtarowicz,
# McGreehanSchotsch) on a sharp-edged bore (edge_radius = bevel_depth = 0) at
# beta <= cb.JOBSON_1955_BETA_MAX (0.75) in both directions; "nozzle" for
# "fixed", a rounded/beveled edge or a wider bore. The cooling-plate and
# connection subclasses keep the nozzle. An explicit choice always wins.
# "nozzle": isentropic nozzle at Cd E A, flat past choke. "jobson1955" (any
# correlation): Jobson's sharp-edged jet, which widens
# as p2/p1 falls, through choke and beyond, with the approach-velocity extension
# in a pipe; Ci from the correlation's Cd = Cv Ci (Cv default 1, Jobson used
# 0.982; Cv != 1 needs compressible_model="jobson1955" set explicitly).
# "iso5167" (ReaderHarrisGallagher only): ISO 5167-2:2003
# Eq. (1) with Eq. (5)'s expansibility. Diagnostics: "compressible_model", "Ci",
# "expansibility", "model_in_range" (Jobson: 0.5 <= Ci <= 0.7, beta <= 0.75;
# ISO: p2/p1 >= 0.75; a flag, never a refusal), and for Jobson "contraction_C", "Cv".
# Kernels (internal solver tools): orifice_jobson_residuals_and_jacobian(...,
# beta, Cv) and orifice_iso5167_residuals_and_jacobian(..., beta) return an
# OrificeModelResult: OrificeResult's columns plus d_mdot_dCd, Ci, contraction,
# expansibility, choked, in_range.
sharp = OrificeElement("sharp", "node1", "node2", diameter=0.01, compressible_model="jobson1955", Cv=0.982)
channel = ChannelElement("channel", "node2", "node3", length=2.0, diameter=0.05, roughness=1e-4, regime="compressible")
# Compressible channel (#481): residual m - m_calc, m_calc = A * fanno_channel_flow(
# feeding node's Pt/Tt, other node's Pt). The duct is fed isentropically from the
# node the flow comes FROM (reversed flow: the to_node), discharges with its exit
# head lost into a PressureBoundary (coupling "static"/default) or recovered into
# anything else, and its flow SATURATES at the choked flow -- a choked duct
# converges on its choked flow, never past it. Diagnostics: "choked",
# "M_in_duct", "M_exit_duct". Where the feeding node mixes gases, its Jacobian
# carries d/dY along the directions that composition moves in (#489).
# On a ThermalWall (#551) the duct carries its heat in the momentum
# (fanno_rayleigh_channel_flow, Tt linear from the feeding node's to the duct's
# exit Tt); diagnostics "rayleigh_in_momentum" = 1. Into a plenum or chamber fed
# by this duct alone, with no wall of its own, the exit Tt is that node's. Into
# any other node (one mixing other streams, one with its own wall, a boundary)
# the solver forms the duct's own exit state from its stream and its own wall
# heat, through the same mixer and relay as a node (#672). Only where no flow
# enters the receiving node does the duct fall back to adiabatic Fanno flow,
# reporting the omitted first-order Rayleigh change ("rayleigh_omitted_frac",
# flag "rayleigh_omitted").
# At an end on a node at rest (#888; `state_at_rest`, Pt = P by the node's own
# definition: PlenumNode, CombustorNode, PressureBoundary, MassFlowBoundary),
# "mach_in" / "mach_out" and the inlet reference block ("velocity", "Re",
# "rho_ref", "mu_ref", ...) are the duct's own state there, not the node's (a
# reservoir's P = Pt, T = Tt; an outflow boundary's static P with the arriving
# Tt): compressible, the Fanno station at that end's Mach (so they equal
# "M_in_duct" / "M_exit_duct"); incompressible, the residual's density state,
# Tt at the approach static P where the node holds a total pressure. An end on
# a moving node (MomentumChamberNode, Pt > P) reports its node's static state;
# ribbed / pin-fin surfaces and CrossflowSegmentElement keep the node states.
# "P_in" ... "Tt_out" stay the nodes'. Reporting only; the heat-transfer
# correlation was already at the duct's inlet static state (#636).

# Effusion (multi-perforated) wall panel. Geometry is given the way a plate
# is designed -- pitch, hole diameter, wall thickness, inclination -- and the
# hole count follows from the panel area.
panel = EffusionPlateElement(
    "panel", "coolant", "gas",
    hole_diameter=0.6e-3, wall_thickness=1.0e-3,
    pitch=3.0e-3, panel_area=0.02 * 0.02, angle_deg=30.0,
)
# One panel is one coolant pressure and one gas pressure, so it cannot show
# coolant migration WITHIN itself. For a wall along a feed channel, hang a
# panel off each channel segment -- the node mass balance then gives
# m_channel_in = m_channel_out + m_effusion, i.e. coolant flow as f(x). The
# resolution is your choice of segment count.

# Connection elements
effective_area = EffectiveAreaConnectionElement("ea", "node3", "node4", diameter=0.015958)
lossless = LosslessConnectionElement("lossless", "node4", "node5")
# validate() refuses a closed loop of lossless connections only -- two in
# parallel between the same nodes, or a ring (#988): Pt is equal all round
# it, so its circulation is not determined. In the GUI a plain (untyped)
# edge is a lossless connection. The same refusal covers any closed loop made
# only of elements whose Pt relation does not depend on their own flow:
# lossless connections, free vortices, a splitter's reference outlet,
# fixed-ratio compressors/turbines, and fraction losses (ConstantFractionLoss,
# LinearThetaFractionLoss) WITHOUT m_ref. Round such a loop the Pt ratios
# either fail to close (no solution) or close exactly (the split is not
# determined). A flow-dependent loss on one path (orifice, channel, head
# loss), or an m_ref on a fraction loss, fixes it.
# Every PressureBoundary fixes its Pt, so validate() treats all of them as one
# node (#1003): a path of these elements from one boundary to another is such
# a loop too and is refused, naming the two boundaries. In the GUI this is a
# plenum with plain edges to two outlets: equal outlet Pt leaves the split
# free (an outlet can even run backwards), unequal leaves no solution. Put a
# flow-dependent loss on each outlet path, or merge the outlets into one
# boundary. These elements couple the boundary's Pt, never its static P (only
# a channel, crossflow segment or port discharge loses an exit head there),
# so a single such link into one boundary is fine.
# A momentum chamber discharging into a PressureBoundary (#637): static
# coupling (the chamber's static P is the boundary's, its head lost) and,
# choked, under-expanded at the sonic static pressure of its flux. Inflow and
# a coupling="total" boundary keep Pt equal. The GUI inserts it where a
# junction port is wired straight to a boundary.
discharge = PortDischargeElement("discharge", "chamber", "outlet_boundary")
area_cd = DiameterDischargeCoefficientConnectionElement("area_cd", "node5", "node6", diameter=0.011284, Cd=0.7)

# Compressor and turbine at a fixed total-to-total pressure ratio (#893):
# Pt_out = PR Pt_in (compressor) or Pt_in / PR (turbine); the network sets the
# flow. The stream leaves at cb.compression / cb.expansion's outlet Tt, with
# exactly one of eta_poly and eta_is. Reports power [W] > 0 (absorbed /
# delivered), dh [J/kg] signed into the flow, Tt_out, eta_poly, eta_is. The
# ledger and the second law count the work as enthalpy with no entropy. A
# reversed-flow root (surge) is refused; between two pressure boundaries
# nothing sets the flow, so that is refused too.
comp = CompressorElement("comp", "inlet_plenum", "burner", pressure_ratio=8.0, eta_poly=0.88)
turb = TurbineElement("turb", "burner", "turbine_exit", pressure_ratio=4.0, eta_is=0.9)

# Flow splitter (#903): one inlet node to N outlet nodes in given mass-flow
# fractions, each outlet a mass-flow controller. A configuration, not an
# element: add_to adds one SplitterOutletElement per outlet ("{id}__{outlet}"),
# a two-port from the inlet node. Each carries one unknown, its m_dot, and one
# row: a metered outlet m_i - f_i sum_j m_j = 0 (linear, exact Jacobian 1 - f_i
# and -f_i), which holds no pressure, so its drop Pt_in - Pt_out is free and
# reported as the valve the rig would need; the REFERENCE outlet Pt_out = Pt_in,
# its valve open (only N - 1 fractions are independent: with a valve on every
# outlet the total flow, or the inlet pressure at a fixed supply, would have
# no equation). The reference must be the path needing the most pressure (a
# metered outlet needing more would need a pump, a root the second-law check
# refuses), and the SOLVE FINDS IT (auto_reference=True, default): it starts
# from reference=, else the remainder, else first_guess (the largest seed
# resistance K f^2 among the paths); after a root refused for pumps at
# metered outlets (or a stalled iterate showing them) it moves the reference
# to the outlet whose node needs the most Pt there and solves again from that
# root, inlet raised to the new reference's Pt. An outer loop, never a switch
# inside Newton; at most N - 1 swaps per splitter, never back to one tried.
# At a fixed supply flow one swap is exact. Reported: __reference_swaps__ =
# {splitter: {"reference", "start", "swapped_from": [...]}},
# __reference_swap_solves__, __reference_swap_time__ [s], each outlet's
# diagnostics "reference_swaps", and a sentence in __message__. The swap stays
# on the network (the elements' reference flags). With no feasible reference
# the first refusal is returned (flags restored, the swaps tried listed, the
# outlet to make the reference named); auto_reference=False keeps the given
# reference and the refusal. Cost (stress splitter family, seeds 0-99, started
# from the LEAST demanding path): 98/100 solve, 1 swap each, 2 solves, wall
# time x3.4 (median) against a right start. Each outlet carries the inlet's Tt
# and Y; Pt follows its own path. One outlet may be the remainder (None, 1 -
# the others); fractions must lie in (0, 1) and sum to 1 within
# SPLITTER_FRACTION_SUM_TOL (1e-9), never renormalised.
# Diagnostics per outlet: fraction, fraction_achieved, m_dot_splitter (the
# throughput), is_reference, valve_dPt [Pa], valve_needs_pump (dPt below
# -SIGMA_TOL Pt_in), and through a positive drop valve_CdA [m^2] (an
# isentropic jet from the inlet's (Tt, Pt) to the outlet's static P, held at
# the sonic pressure: valve_choked), with areas= also valve_Cd = CdA / area and
# valve_zeta = 1 / Cd^2 (the K over the dynamic head at that area).
from combaero.network import FlowSplitter
split = FlowSplitter(
    "split", inlet="coolant_supply",
    outlets={"hole1": 0.2, "hole2": 0.3, "hole3": None},   # hole3: the remainder
    reference="hole1",                    # the start; default: the remainder
    auto_reference=True,                  # move it if infeasible (default)
    areas={"hole2": 2e-4, "hole3": 2e-4},  # optional, for valve_Cd / valve_zeta
)
split.add_to(net)                          # the nodes must be in net already
# split.summarize(result["__element_diag__"]) -> {"m_dot", "reference" (in
# use), "reference_start", "needs_pump", "outlets": {node: {m_dot, fraction,
# valve_dPt, valve_CdA, ...}}}

# Three-port tee junction (Unified0D compressible model, Mynard & Valen-Sendstad 2015)
import math
# Merging tee: two inlets (straight_node, branch_node) -> one outlet (common_node)
tee_merge = TeeJunctionElement(
    "tee", common_node="outlet", straight_node="inlet_s", branch_node="inlet_b",
    theta=math.pi / 2.0, F_C=0.01, psi=1.0, tee_type="merging",
)
# Branching tee: one inlet (common_node) -> two outlets (straight_node, branch_node)
tee_branch = TeeJunctionElement(
    "tee", common_node="inlet", straight_node="outlet_s", branch_node="outlet_b",
    theta=math.pi / 2.0, F_C=0.01, psi=1.0, tee_type="branching",
)
# Solved unknowns: tee.m_dot_com (total), tee.m_dot_branch (branch arm)
# Straight flow is implicit: m_dot_straight = m_dot_com - m_dot_branch
# A given F_C (common and straight arm) or F_branch that contradicts the single
# pipe at its arm is refused at resolve_topology (#1023, as AreaChangeElement's
# faces). A branch wider than the main (psi < 1) is accepted.

# Momentum-CV junction (PDF spec, supersedes K-closure for n>3 manifolds and
# high-Mach / ejector behaviour; see sources/junction/momentum cv implementation guide.pdf).
# N port-MCNs feed a single junction, whose closure carries every port's loss --
# connect the ports with plain elements (channels, LosslessConnectionElement).
mc_com = MomentumChamberNode("mc_com", area=0.01)
mc_str = MomentumChamberNode("mc_str", area=0.01)
mc_bra = MomentumChamberNode("mc_bra", area=0.008)

# MultiPortChamberBase is the ABSTRACT BASE since 0.6.0 -- it owns the port
# machinery and no longer carries a residual. Use MultiPortChamberElement (the Mynard
# closure) or ConstantKTeeElement (fixed handbook K).
jct = MultiPortChamberElement(
    "jct",
    inlet_nodes=["mc_com"],
    outlet_nodes=["mc_str", "mc_bra"],
    inlet_angles_deg=[0.0],
    outlet_angles_deg=[0.0, 90.0],   # geometric branch angles per outlet port
    flow_direction="branch",         # 1 supplier + N-1 collectors
    # port_areas (length = N_in + N_out) inherited from connecting channels if not given
)
# Do NOT add a turning-loss element on a port: the Mynard closure already
# carries the lateral turning loss, and the two together double-count it
# (Bassett K6 at 90 deg, psi = 1, q = 0.5: 0.867 alone, which is Bassett's
# value; 1.249 with the removed BorderCarnotLossElement; #272).
# Solved unknowns: jct.Pt_jct -- for this closure the junction TOTAL pressure
# (jct.P_jct before #773; the old initial_guess key was removed in 0.8.0).
# Result keys (#773): P_jct is the
# junction STATIC pressure (Pt_jct minus the common port's dynamic head,
# i.e. the common port's static), Pt_jct the TOTAL; every junction element
# reports both. Per-port mass flows
# live on the connecting channels / loss elements; the junction reads them via
# the graph. The N+1 residuals = N port total-pressure relations
# (Pt_i - Pt_jct +/- K_i * q_dyn_com) + global sum-of-port-mdots = 0, computed
# in C++ (see mpce_residuals_and_jacobian in API_CPP.md).

# Constant-K junction (the "simplest model" tier): fixed handbook loss
# coefficients instead of the Mynard closure. K is referenced to the
# common-leg dynamic head (Idelchik convention) and does not vary with the
# flow split. K_ports keys are port indices (inlets first, then outlets);
# the common port (single inlet for 'branch', single outlet for 'merge')
# is ignored. Exposed in the GUI as the tee's "Constant K (handbook)"
# junction model.
from combaero.network.mpce_element import ConstantKTeeElement
jct_k = ConstantKTeeElement(
    "jct",
    inlet_nodes=["mc_com"],
    outlet_nodes=["mc_str", "mc_bra"],
    flow_direction="branch",
    K_ports={1: 0.4, 2: 1.0},   # K_straight, K_branch
)
# A constant K belongs to one split: the branch K grows like (q psi)^2, so a
# value for an equal-area tee is several times too low for a narrow branch.
# bassett_design_K gives Bassett 2001's analytical pair (dividing K5/K6,
# joining K11/K12 angle-corrected) at a DESIGN branch share; the GUI uses it
# for any K left blank (design_branch_share, default 0.5).
K_str, K_bra = ConstantKTeeElement.bassett_design_K(
    "branch", q_branch=0.5, psi=4.0, theta_deg=90.0   # -> (0.0, 3.47)
)
# For a MERGING tee Bassett's design pair can be net-negative at shallow
# angles (psi = 1, theta <= 60 deg; at 30 deg and q = 0.5 both K are -0.04,
# against +0.05 measured): the solve then converges to a junction that
# raises total pressure and the second-law check rejects it (#779). Set
# K_ports explicitly there, or use MultiPortChamberElement.
# design_net_loss is the pair's mass-weighted loss at its design share,
# (1 - q) K_straight + q K_branch per common-leg head; below zero the pair
# creates total pressure (one negative K is physical, a negative sum is not).
# The GUI refuses such a pair at build time with this cause.
ConstantKTeeElement.design_net_loss(0.5, *ConstantKTeeElement.bassett_design_K(
    "merge", q_branch=0.5, psi=1.0, theta_deg=30.0))   # -> -0.04

# Ideal momentum junction (#517): parameter-free, a baseline beside the tuned
# closures rather than a replacement. Suppliers sit at one junction STATIC
# pressure p_jn (the solved unknown and result key P_jct; static, not Pt).
# Every joining closure shares this one-static-pressure property (Mynard and
# Bassett too), so with equal supply pressures and equal upstream legs a
# merging tee splits exactly 50/50 at any angle -- by construction, and
# within the scatter of Blaisdell & Manson's measured statics. Each receiver
# obeys
# P_j - p_jn + (m_j/A_j)(u_j - V.n_j) = 0, V the supplied mass-weighted
# velocity vector. Merge = Blaisdell & Manson's theory; split = vector
# Borda-Carnot loss on every receiver (Idelchik Eq 7-6, A' = 1). No
# flow_direction: roles follow the flow smoothly, 2 to 6 ports. Angles as for
# the other tees (outward axis = angle for outlets, pi - angle for inlets).
# Any port may converge reversed (the solver feeds a reversed junction port by
# the sign of its flow); on an all-pressure network that can make more than
# one root admissible. Range: subsonic supplying ports -- a supplier that would
# choke (Wang 2014, a = 2.44, high q and M3) is driven supersonic and the solve
# stalls at Mach 1; not refused or flagged. Kernel: ideal_junction_residuals_and_jacobian
# (API_CPP.md).
from combaero.network import IdealMomentumTeeElement
jct_ideal = IdealMomentumTeeElement(
    "jct",
    inlet_nodes=["mc_com"],
    outlet_nodes=["mc_str", "mc_bra"],
    inlet_angles_deg=[0.0],
    outlet_angles_deg=[0.0, 90.0],
)
# Diagnostics: P_jct (= p_jn, static) and Pt_jct (P_jct plus the common
# port's dynamic head; p_jn was a duplicate of P_jct until #773);
# port_i_K referenced to the converged common port's head;
# for a 3-port tee the geometric flows m_dot_com/straight/branch (signed in the
# wired direction) and n_reversed_ports, plus K_straight/K_branch (dividing) and
# K11/K12 (joining), each with a value only when the converged regime is that
# one and NaN otherwise (every key is always present). Exposed in the GUI as the
# tee's "Ideal momentum (parameter-free)" junction model.
#
# Every multi-port junction also reports port_i_outer_area (the connected
# element's own area), port_i_connected_area_ratio (that over the port area;
# port_i_area_ratio before #773) and
# port_area_mismatch (ports off by more than PORT_AREA_MISMATCH_TOL = 5%): an
# area step there is not modelled, and this closure, coupled on static
# pressure, then creates total pressure across it (#536).
#
# loss_factor (#526): scales each RECEIVER's loss against the mixed supplied
# stream, Pt_j = Pt_mix - f rho/2 |v_j - V|^2; the suppliers' mixing loss is
# not scaled. One value for all ports or one per port (inlets, then outlets).
# Default 1.0 = the parameter-free closure, bit for bit; any f >= 0 keeps the
# mechanical dissipation >= 0. In a split (no mixing) f scales the whole loss,
# and below f ~ M^2/4 the energy check refuses it (the loss is written at low
# Mach, the port Pt is isentropic). A rig-matching knob for the user; the
# harness scores only f = 1.
jct_tuned = IdealMomentumTeeElement(
    "jct", inlet_nodes=["mc_com"], outlet_nodes=["mc_str", "mc_bra"],
    inlet_angles_deg=[0.0], outlet_angles_deg=[0.0, 90.0], loss_factor=0.8,
)

# Free-vortex swirl redistribution from r_in to r_out (#661); keyword-only
vortex = VortexElement(
    "vortex", "node_in", "node_out",
    r_in=0.05,        # inlet radius [m] (0 = on-axis: no angular momentum)
    r_out=0.10,       # outlet radius [m]; r_out < r_in is inward flow
    shaft_speed_rpm=3000,  # shaft speed [rpm]; sets the default inlet swirl
                      # (omega_rpm before #773; removed in 0.8.0)
    v_theta_in=None,  # inlet swirl [m/s]; None = solid body, omega * r_in
)
# Residual: Pt_out - Pt_in = 0 (no work: Pt and Tt conserved).
# r * V_theta is conserved, so V_theta_out = V_theta_in * r_in / r_out; the
# swirl's dynamic head turns into static pressure outward (and back inward).
# Diagnostics: m_dot, shaft_speed_rpm, v_theta_in/out [m/s], P_swirl_in/out [Pa]
# (static at each radius, isentropic from the stagnation state),
# dP_swirl_static [Pa], M_theta_in/out, Pt_in, Pt_out.
# r_c and n (the pre-#661 Vatistas core) were removed in 0.8.0; the
# cb.vatistas_* functions below are unchanged.

# Supersonic ejector across all operating regimes -- critical (double-choked),
# subcritical droop, and unchoked-primary subsonic jet pump -- in one C1 residual
# system (Huang 1999 entrainment ratio; anchors by anchor_closure, Huang 1999 by
# default or Kracik & Dvorak 2016; see validation/ejector/README.md and OPERATING_REGIMES_DESIGN.md).
from combaero.network.ejector_element import EjectorElement
ejector = EjectorElement(
    "ej1", primary_node="mc_primary", secondary_node="mc_secondary", outlet_node="mc_outlet",
    throat_area=3.14e-5,      # primary nozzle throat area A_t [m^2]
    nozzle_exit_area=1.0e-4,  # primary nozzle exit area A_p1 [m^2], must exceed throat_area
    mixing_area=8.0e-4,       # mixing section A_3 [m^2], > nozzle_exit_area; None = inherit
                              # from the outlet pipe (else DEFAULT_EJECTOR_MIXING_AREA);
                              # a given one larger than the single outlet pipe is
                              # refused (#1023), a wider pipe (a diffuser exit) is not
    recovery_efficiency=1.0,  # in (0, 1]: multiplies the lossless mixed stagnation P; not fitted
    regime="compressible",    # default DEFAULT_REGIME ("compressible" since #585), see below
    anchor_closure="huang1999",  # the (P_c*, P_b0) pair (default): Huang's printed P_c*
                              # chain + the sudden-expansion dead head, 2-3x lower curve
                              # error cross-source; or "kracik_dvorak", the lossless pair
                              # (the only one up to 0.6.0; validation/ejector/README.md)
)
# On a failed solve whose outlet boundary sits above the dead head P_b0, the
# solver's message says so: a real ejector backflows there, outside the model.
# regime="incompressible" (#588): the ejector's low-Mach limit, Sanger's 1-D jet
# pump (NASA TN D-4445, 1968, App. B) with each stream at its own constant density
# -- no choking, no critical regime; the flow follows the pressure ratio. Same
# owned unknown (P_py, the mixing-plane static) and four rows: nozzle and suction
# Bernoulli with losses, mass, and the constant-area mixing-duct momentum to an
# outlet Pt = recovery_efficiency * (p4 + (1 - k_diffuser) q4). Loss coefficients
# k_nozzle, k_suction (default: the compressible closure's own nozzle efficiencies
# as K = 1/eta - 1, 0.053 / 0.176, so both regimes agree at low Mach), k_throat,
# k_diffuser (default 0). It reproduces Sanger's Eq. 11 to machine precision; with
# his calibrated K's his theory scores 9.3% on his own Table II (D-4759). The
# compressible ejector reduces to it as Ma^2 (gap 0.1% at a jet Ma of 0.03, equal
# and unequal stream temperatures). Its diagnostics report entrainment_ratio, Sanger's head
# ratio N = (Pt_o - Pt_s)/(Pt_p - Pt_o), efficiency M N, and mach_jet, which the
# incompressible validity report (#592) judges: a gas ejector of any real drive
# runs its jet past Ma 0.2 and is flagged.
# Solved unknowns: the three port mass flows + ej1.P_jct, the single owned scalar,
# repurposed as the mixing-plane static pressure P_py. Four residual rows: primary
# C-D nozzle (choked or unchoked), blended entrainment, mass conservation, and a
# regime-dependent outlet/discharge closure. The regime is picked smoothly by two
# smootherstep weights (s_choke on mp/mdot_choked, s_sub on outlet.Pt/P_c*): in
# choked mode the outlet floats free below P_b0 and P_jct is the diagnostic
# recovered P_c*; in the jet-pump regime (and choked beyond P_b0) the outlet is
# pinned to the ejector's discharge and P_py floats. Working fluid: combaero's real-fluid
# EOS (combustion/air species only, no halocarbon refrigerants); gamma/R evaluated
# live at the entrained-flow choking plane. The whole 4-row (f, J) is assembled in
# C++ (full analytic Jacobian, no FD fallback). diagnostics() reports the actual
# entrainment_ratio = ms/mp (Huang's omega; the internal API keeps that name),
# a critical_mode flag, nozzle_exit_mach_out_of_range (#945: +1 or -1 when
# area_ratio_nozzle lies outside A/A* of the nozzle-exit Mach bisection's
# bracket [1, 50], so M_p1 is that bracket end, clamped; 0 otherwise), and (in
# critical mode) P_c_star and entrainment_ratio_critical. Pressure keys
# follow P = static, Pt = total (#773): P_jct is the mixing-plane static P_py,
# Pt_out the outlet total;
# there is no Pt_jct, as the plane carries two streams at one static
# pressure and different velocities, hence no single total.
# verify_solution_consistent(sol) returns True unconditionally now that every
# regime is modeled.
# Warm start: a compressible ejector seeds its own cold start at every
# resolve_topology (ejector.seed_warm_start(graph); the solver calls it): each
# port traced to its boundary, port chamber P/Pt, P_jct (= P_c*, or just below
# the back pressure for a mass-driven primary that cannot choke) and the port
# chains' flows from the closed forms, only where you gave no initial_guess.
# It replaces its own seeds on a re-solve and seeds nothing incompressible.
# Choked SUBCRITICAL operation (#619): between P_c* and the dead-head P_b0 (the
# same mixing closure at zero entrainment) the entrainment droops linearly from
# omega_crit to 0 (Tier-1, design doc sec 8.1) and the outlet floats; the outlet
# pin switches in at P_b0. diagnostics() reports P_b0 next to P_c_star.
# Beyond P_b0 (backflow) is outside the model.
```

### Area Change

`AreaChangeElement(id, from_node, to_node, F0=None, F1=None, model_type="sharp" | "conical", length=None, D_h=0.0, regime="compressible")`
is a sudden or conical expansion or contraction from face `F0` to `F1` (#731).

- **Closures: Idelchik (1966),** transcribed in `validation/area_change/` and scored on every printed cell (`python/tests/test_area_change.py`):
  - sudden expansion, diagram 4-1: `cb.sudden_expansion_zeta(Re, n)` returns `(zeta, dzeta_dRe)`. It is 26/Re below Re 8, a C1 interpolation of the printed table to Re 3500, and Borda-Carnot `(1 - n)^2` above. Three of the 78 printed cells move: two to keep the loss rising with the flow (n 0.6 and 0.5 at Re 2000), and one rounded Borda-Carnot value;
  - sharp contraction, diagrams 3-10 and 3-9: `cb.sudden_contraction_zeta(Re, n)`: 27/Re below Re 8, Karev's table (every printed cell) to Re 1e4, `0.5 (1 - n)` above;
  - conical diffuser, diagram 5-2: `cb.conical_diffuser_zeta(full_angle_deg, n)`, exact at every printed cell. Friction is not included;
  - conical contraction, diagram 3-7: `cb.conical_nozzle_zeta(full_angle_deg, l/D, n)` = zeta' (1 - n), zeta' from diagram 3-6 (exact at every printed cell). Friction is not included.
- **Mach.** The tables are incompressible; the element carries them by their control volume.
  - An expansion is the compressible Borda-Carnot step (`cb.compressible_borda_carnot`) times the table's ratio to it. At n 0.8 and M 0.62 that is 1.62x Borda-Carnot.
  - A contraction is a vena contracta re-expanding by the same step. It chokes, and the element's flow then stays at the sonic flow.
  - `regime="incompressible"` is the M -> 0 limit.
  - The gas is calorically perfect across the step, at gamma(Tt): within 0.2% of the thermally perfect control volume to M 0.65 at 1500 K.
- **Residual: mass flow,** `m - m_calc(Pt_from, Pt_to, ...)`, as the compressible channel's and orifice's.
  - The flow rises monotonically with the drive and saturates at the sonic flow.
  - Each direction is fed from its own upstream node (reversed flow reads the `to_node`). The two are blended C1 through zero drive.
  - Below a branch's drive floor (zeta times 2 Pa of head) the flow is a C1 odd blend with a finite slope at zero.
  - Every partial is analytic, composition included (`cb.area_change_element_flow`).
- **Equal faces are no area change** (#1023, until 0.10.0 refused): F0 == F1 to `NEIGHBOUR_AREA_RTOL` (1e-6, #788's bore-pipe band) on the equivalent diameter `sqrt(F0 / F1)` is a pass-through, as a bore equal to its pipe is. Zero loss, the row is `Pt_in - Pt_out = 0` with Jacobian +1/-1 (row kind "p"), and the resolve preview reports `pass_through = 1`; diagnostics report `dPt_loss = zeta = 0` and the duct's own Mach as `mach_loss_ref`. Two in parallel are a flow-independent loop (refused at validate, as two lossless connections are).
- **An explicit face against its pipe** (#1023): a given `F0` (`F1`) is compared with the single `ChannelElement` at the `from_node` (`to_node`), and one that differs by more than `NEIGHBOUR_AREA_RTOL` on the equivalent diameter is refused at `resolve_topology` (`ValueError` naming the element, the pipe and both areas). Equal is consistent; an inherited face is the pipe's by construction. Not compared: a node that joins more than one pipe, a pipe at the GUI's `DEFAULT_CHANNEL_DIAMETER`, and a face whose node is a `MomentumChamberNode` or `CombustorNode` with its own (user) area or Dh -- the face meets that chamber, and the chamber against the pipe is the chamber's warning (see `MomentumChamberNode`).
- **Diagnostics** (upwind node's state):
  - `dPt_loss`: the closure's stagnation loss (formerly `dP_static`);
  - `zeta`: on the small section's static head;
  - `mach_loss_ref`: the small section's Mach for an expansion, the vena contracta's for a contraction;
  - `Re_small`, `expansion` and `choked`.

### Flow-Area Inference

Components that need a reference flow area infer one from the topology when you
do not supply it. Inference searches the neighbourhood of the attached nodes:
channels first (so networks that resolved before resolve identically), then any
other area-bearing element (`AreaChangeElement`, `TeeJunctionElement`,
`PressureLossElement`), then adjacent nodes whose
area you set explicitly (`MomentumChamberNode`, `CombustorNode`). An area that
was itself auto-sized is never used as a source -- inheriting a placeholder
would launder one guess into another.

What happens when nothing can be inferred depends on what the area means:

| Component | Parameter | Unresolvable |
|---|---|---|
| `AreaChangeElement` | `F0`, `F1` | **raises** |
| `ChannelElement` | `diameter` | **raises**; with `default_diameter=` set (the GUI passes `DEFAULT_CHANNEL_DIAMETER` = 0.1 m, #570), that default, reported as `area_source = "default"` |
| `TeeJunctionElement` | `F_C` | **raises** |
| `MomentumChamberNode` | `area` | nominal default |
| `CombustorNode` | `area` | nominal default |
| `PressureLossElement` | `area` | nominal default, warns if a head-loss correlation reads it |

The first group raises because the area *is* the geometry being modelled -- a
default there invents a contraction or expansion the network does not contain,
and the solver either stalls on it or converges to a confidently wrong answer.
The second group defaults because the area is only a dynamic-head reference and
a network can legitimately never read it.

**An inferred flow area never becomes a wetted area** (#1022). A
`MomentumChamberNode`, `CombustorNode` or `PressureLossElement` keeps the
`surface=ConvectiveSurface(area=...)` you gave it, whatever its flow area
resolves to; with no surface given, `surface.area` stays 0 (no convective
surface), and a `ThermalWall` on it is refused by name. Give the wetted area
explicitly to wall such a component. The GUI, which has no wetted-area field
on these components, takes it from the wall drawn to them (that side's area,
else the wall's area).

```python
# Unresolvable -> a named error saying which parameter clears it
AreaChangeElement("expansion", "plenum_a", "plenum_b")
# ValueError: AreaChangeElement 'expansion': cannot determine upstream area F0
# -- no neighbouring channel, area-bearing element or node with a known area
# was found at 'plenum_a'. Set upstream area F0 explicitly ...

# Either fix clears it: set the value, or give a neighbour a known area
AreaChangeElement("expansion", "plenum_a", "plenum_b", F0=0.05, F1=0.08)
MomentumChamberNode("chamber", area=0.15)  # neighbour the area change can read
```

### Tuned Constants in the Junction Closure

`MultiPortChamberElement`'s Mynard closure carries three tuned constants, each named
and documented at its definition in `combaero.network._mynard2010`:

| constant | value | what it is | switch |
|---|---|---|---|
| `MYNARD_ETA_A0`, `MYNARD_ETA_A1` | 0.8, -0.2 | Mynard 2015 Eq 36 energy-transfer factor, CFD-fitted | `eta_scale` |
| `FLOW_RATIO_DAMPING` | 0.02 | regulariser from the Matlab reference, not physics | -- |
| `MultiPortChamberElement.DEFAULT_JOINING_ETRANSFER_ALPHA` | 0.2 | combaero's joining-side correction, fitted by `validation/junction/calibrate_etransfer.py` | pass `joining_etransfer_alpha=0.0` |

Under this repo's policy a tuned constant earns its place only by improving
agreement on the digitised validation data; the on/off tables live on issue
#271. `eta_scale` exists for that measurement:

```python
MultiPortChamberElement(..., eta_scale=0.0)   # default: energy-transfer term off
MultiPortChamberElement(..., eta_scale=1.0)   # Mynard's Eq 36 term on

# Mynard's published model -- the dividing-streamline recovery is a combaero
# addition (on by default) and has to be switched off too:
junction_loss_coefficient(U, A, theta, joining_etransfer_alpha=0.0,
                          eta_scale=1.0, dividing_streamline_recovery=0.0)
```

That configuration is held to Mynard 2015's own figures (Figs 4, 6-11) by
`python/tests/test_junction_mynard_fidelity.py`; see
`validation/junction/README.md`, "Fidelity".

Changing any of these values is a retune and needs a before/after table
before it lands, not a silent edit -- `python/tests/test_junction_tuned_constants.py`
pins the documented values.

### The Junction Soft-Barrier Weight

Distinct from the tuned constants above: this is a **numerical** parameter, not
physics.

When `strict=False` and a port flows against its declared direction,
`MultiPortChamberElement` replaces the physics with a continuity residual plus a
one-sided quadratic penalty, `alpha * max(0, -e_i * mdot_i)^2`. The penalty
shares its row with the continuity relation, so it balances against a pressure
error rather than driving the offending flow to zero, and has a fixed point at
`slack* = sqrt(dP / alpha)`. A solve that reaches it parks there.

`alpha` is `soft_penalty_alpha`, default `DEFAULT_SOFT_PENALTY_ALPHA = 1e11`
Pa/(kg/s)^2. It carries units, so it places the fixed point at a sensible
fraction of the flow for one network size only. The solver used to derive it
from the network's scales; that hand-off was removed once the junction seed
stopped starting ports in the wrong basin and it no longer changed any outcome
from 1e-6 to 1e4 times the reference size (#272). Set it explicitly for a
network far from ~0.1 kg/s that parks a port in the barrier:

```python
element.soft_penalty_alpha = 5.0e7
```

Raising `alpha` shrinks the fixed point and never destabilises the solve --
the response is monotone and saturates -- so when in doubt, larger.

### Combustion Integration
```python
from combaero.network.combustion import (
    combustion_from_streams, combustion_from_phi, mix_streams, stoichiometric_products
)

# Create combustion from streams
# cb.MassStream(m_dot, Tt, Pt, Y): total temperature and pressure
fuel_stream = cb.MassStream(0.01, 300.0, 2e5, fuel_y)
oxidizer_stream = cb.MassStream(0.2, 300.0, 2e5, air_y)

result = combustion_from_streams([fuel_stream, oxidizer_stream], phi=1.0)
print(f"Products: {result.products.Y}")
print(f"Adiabatic T: {result.T_ad} K")

# Or from equivalence ratio
result = combustion_from_phi(phi=0.8, T=300, P=2e5, fuel="CH4", oxidizer="air")
```

### MixtureState for Network Data
```python
import combaero as cb
from combaero import MixtureState

air_Y = list(cb.species.dry_air_mass())

# Constructor: MixtureState(P, Pt, T, Tt, m_dot, Y); P, T static, Pt, Tt total
node_state = MixtureState(P=1e5, Pt=1.1e5, T=292.0, Tt=300.0, m_dot=1.5, Y=air_Y)

print(f"Static pressure: {node_state.P} Pa")
print(f"Total pressure: {node_state.Pt} Pa")
print(f"Mass flow: {node_state.m_dot} kg/s")
print(f"Static density: {node_state.density()} kg/m^3")
```

The network's own node state is `combaero.network.NetworkMixtureState`, which
derives the static `T` from `Tt`; build it with keywords
(`NetworkMixtureState(P=..., Pt=..., Tt=..., m_dot=..., Y=...)`).

### Advanced (f, J) Interface
High-accuracy analytical Jacobians are available for performance-critical solver loops.
```python
# Exact residuals and derivatives wrt (P_tot, P_static, Tt, Y); internal
# (combaero._solver_tools), used by the network elements. Tt_up is the
# upstream TOTAL temperature (T_up before #782); res.d_mdot_dTt_up,
# res.d_dP_dTt_up its columns.
from combaero import _solver_tools as st
res_ori = st.orifice_residuals_and_jacobian(m_dot, P_tot, P_stat, Tt, Y, P_down, Cd, area)
res_chan = st.channel_residuals_and_jacobian(m_dot, P_tot, P_stat, Tt, Y, P_down, L, D, roughness, correlation,
                                             f_multiplier=1.0, exit_head_lost=False)  # True adds the exit head (#588)
# Mix and burn streams (cb.MassStream), Q [W] and 'fraction' after combustion
res_comb = st.combustor_residuals_and_jacobians(streams, P, Q=0.0, fraction=0.0,
                                                pressure_loss=loss, use_equilibrium=False)
```

The combustor returns a `MixerResult`: the products' total temperature
`Tt_mix`, `Y_mix` and total pressure `Pt_mix`, with their derivatives per
stream (`dTt_mix_d_stream`, `dPt_mix_d_stream`, `dY_mix_d_stream`). The
pre-0.8 spellings `T_mix`, `P_total_mix`, `dT_mix_d_stream`,
`dP_total_mix_d_stream`, `dT_mix_d_delta_h` and `dT_mix_dQ` were removed in
0.9.0.

### Tee Junction Interface

Two solver-facing APIs are available:

**Unified0D compressible model** (current, used by `TeeJunctionElement`): stagnation-pressure
residuals with mass-flow continuity, consistent with compressible duct elements. Supports
arbitrary per-branch gamma and R_gas.

**Legacy Bassett 2001 model**: empirical incompressible K tables retained for direct
low-level access; not used by the network solver.

#### Unified0D compressible interface

`BranchInput` holds per-branch thermodynamic state. `CompressibleTeeResult` holds the
two residuals and the full Jacobian.

```python
import combaero._core as _core
import math

# Build per-branch state (P_static [Pa], Pt [Pa], T [K], m_dot [kg/s],
#                         A [m^2], theta [rad], gamma_eff [-], R_gas [J/kg/K])
com = _core.BranchInput(P_static=1.95e5, Pt=2.0e5, T=400.0, m_dot=0.5,
                         A=0.01, theta=0.0, gamma_eff=1.4, R_gas=287.0)
str_b = _core.BranchInput(P_static=1.88e5, Pt=1.95e5, T=400.0, m_dot=0.3,
                           A=0.01, theta=0.0, gamma_eff=1.4, R_gas=287.0)
bra = _core.BranchInput(P_static=1.82e5, Pt=1.90e5, T=400.0, m_dot=0.2,
                         A=0.01, theta=math.pi/2, gamma_eff=1.4, R_gas=287.0)

# Branching: common=supplier, straight+branch=collectors
res = _core.compressible_branching_tee_rj(com=com, str=str_b, bra=bra)
# res.R_0, res.R_1           -- residuals [Pa]
# res.dR0_dPt_com, res.dR0_dPt_str, ...  -- full Jacobian fields

# Merging: straight+branch=suppliers, common=collector
res = _core.compressible_merging_tee_rj(com=com, str=str_b, bra=bra)
```

#### Legacy Bassett 2001 interface

```python
import combaero._core as _core
import math

Y = [0.0] * 15
Y[0] = 0.767  # N2
Y[1] = 0.233  # O2

# Check input validity (non-throwing)
status = _core.tee_check_inputs(q=0.5, psi=1.0, theta=math.pi/2)
# status.valid, status.q_in_range, status.psi_valid, status.theta_valid

# Raw K-coefficient functions (Bassett 2001 Table 2 + Eq 33/34 angle corrections).
# K2 and K5 take only q (psi/theta-independent per Eq 15); the rest take (q, psi, theta).
# Provided for diagnostics and validation against measured data.
k1  = _core.tee_K1(q, psi, theta)
k2  = _core.tee_K2(q)
k3  = _core.tee_K3(q, psi, theta)
k4  = _core.tee_K4(q, psi, theta)
k5  = _core.tee_K5(q)               # straight arm, separating tee
k6  = _core.tee_K6(q, psi, theta)   # branch arm, separating tee
k7  = _core.tee_K7(q, psi, theta)
k8  = _core.tee_K8(q, psi, theta)
k9  = _core.tee_K9(q, psi, theta)
k10 = _core.tee_K10(q, psi, theta)
k11 = _core.tee_K11(q, psi, theta)  # straight arm, joining tee
k12 = _core.tee_K12(q, psi, theta)  # branch arm, joining tee

# Blended K functions (always finite, smooth at q=0 and across topology reversal)
ks = _core.merging_tee_K_straight(q, psi, theta)
kb = _core.merging_tee_K_branch(q, psi, theta)

# Legacy solver residuals and full Jacobians
res = _core.merging_tee_residuals_and_jacobian(
    m_dot_com=0.5, m_dot_branch=0.2,
    dP0_straight=150.0, dP0_branch=200.0,
    P_static_com=2e5, T_com=600.0, Y_com=Y,
    theta=math.pi/2, psi=1.0, F_C=0.01,
)
# res.R_straight, res.R_branch
# res.dR_straight_d_mdot_com, res.dR_straight_d_mdot_branch
# res.dR_straight_dP_static_com, res.dR_straight_dT_com
# res.topology_valid, res.status  (CorrelationValidity.VALID or EXTRAPOLATED)

res = _core.branching_tee_residuals_and_jacobian(
    m_dot_com=0.5, m_dot_branch=0.15,
    dP0_straight=100.0, dP0_branch=180.0,
    P_static_com=1.5e5, T_com=500.0, Y_com=Y,
    theta=math.pi/3, psi=1.2, F_C=0.008,
)
```

---

## Vatistas n-Vortex Model

Reference: Vatistas, Kozel, Mih (1991), *Exp. Fluids* 11, 73–76.

Shape parameter `n >= 1` (default 2.0). `n=2` gives the best fit to most
experimental swirl data.  Closed-form antiderivatives for the pressure integral
are used for `n=1` and `n=2`; composite Simpson quadrature for general `n`.
All Jacobians are analytical.

### Free functions

```python
import combaero as cb

# Normalised shape functions (no units, r_bar = r / r_c)
v0   = cb.vatistas_v0_bar(r_bar, n)           # tangential velocity shape
dv0  = cb.vatistas_dv0_bar_drbar(r_bar, n)    # d(V0_bar)/d(r_bar)
vr   = cb.vatistas_vr_bar(r_bar, n)           # radial velocity shape (< 0)
dvr  = cb.vatistas_dvr_bar_drbar(r_bar, n)    # d(Vr_bar)/d(r_bar)
I    = cb.vatistas_pressure_integral(r_bar, n) # integral_0^r V0^2/r' dr'
dI   = cb.vatistas_d_pressure_integral_drbar(r_bar, n)  # V0_bar^2/r_bar (FTC)

# Dimensional functions (Gamma [m^2/s], r_c [m], rho [kg/m^3])
Vt   = cb.vatistas_v_theta(r, Gamma, r_c, n)  # [m/s]
Vt, dVt_dr, dVt_dG, dVt_drc = cb.vatistas_v_theta_and_jacobians(r, Gamma, r_c, n)

dP   = cb.vatistas_delta_p(r, rho, Gamma, r_c, n)  # P(r)-P(0) [Pa]
dP, ddP_dr, ddP_dG, ddP_drc = cb.vatistas_delta_p_and_jacobians(r, rho, Gamma, r_c, n)
```

### `VatistasVortex` class

```python
import numpy as np
import combaero as cb

vortex = cb.VatistasVortex(Gamma=0.5, r_c=0.05, n=2.0)

r = np.linspace(0, 0.2, 200)
V  = vortex.V_theta(r)          # [m/s], vectorised
dV = vortex.dV_theta_dr(r)      # [1/s], vectorised

Vmax = vortex.V_theta_max()     # peak tangential velocity at r = r_c

dP   = vortex.delta_P(r=0.1, rho=1.2)    # [Pa], scalar
dPbar = vortex.delta_P_bar(r=0.1)        # normalised [0, 1)
```

---

## NetworkRunner (GUI-JSON Programmatic Driver)

`NetworkRunner` loads a network saved from the GUI and runs it from Python
scripts or Jupyter notebooks — no web server required.
It configures the solver exactly as the GUI does (`gui/backend/solver_setup.py`):
the Newton-homotopy last retry on, bounded only by the solve's timeout
(default 30 s, `DEFAULT_SOLVE_TIMEOUT`), and the drawing's incompressible Mach
limit. A saved drawing solves the same from a script as on the canvas.

A worked notebook, [`python/examples/orifice_flow_from_dp.ipynb`](../python/examples/orifice_flow_from_dp.ipynb), loads the shipped ISO 5167 orifice example and sweeps a table of supply pressures and pressure differences into mass flows.

### Loading a network

```python
from gui.backend.runner import NetworkRunner

# From a file downloaded from the GUI
runner = NetworkRunner.from_file("my_network.json")

# From an already-loaded dict (normalises GUI camelCase keys automatically)
import json
with open("my_network.json") as f:
    runner = NetworkRunner.from_dict(json.load(f))
```

Node labels are read from the GUI's label field.  If a node has no label the
node ID is used as the fallback key for overrides and result lookup.

### Single solve with boundary-condition overrides

```python
# Override format: "<label>.<attribute>": value
result = runner.solve({
    "air_inlet.m_dot": 1.2,   # kg/s
    "air_inlet.Tt": 650.0,    # K
    "outlet.Pt": 200_000.0,   # Pa
})

print(result.success, result.final_norm)

# Scalar result by label.quantity (resolves label → node ID automatically)
T_combustor = result.get("combustor.Tt")         # K
phi          = result.get("combustor.phi")
m_dot_out    = result.get("hot_channel.m_dot")  # kg/s

# Full thermodynamic state dict for a node
state = result.node_state("combustor")
# {"T": 1738.0, "P": 195000.0, "Pt": 200000.0, "Tt": 1760.0, ...}

# Full-detail DataFrame (matches GUI CSV export, unit-annotated columns).
# Every row carries the solve's status -- success, converged, consistent,
# outcome, final_norm -- so filter on it before trusting the numbers: a
# timed-out or second-law-rejected solve still produces a full frame.
# A pressure boundary's P is its Pt (a reservoir: static = total).
df = result.to_dataframe()
df = df[df["success"]]
```

The GUI's **Export CSV** writes this frame for the displayed result (#783):
it solves again only when no solve of the canvas is held, and the header
`X-Combaero-Export-Source` (`displayed` / `solved`) says which. Its response
headers `X-Combaero-Success` / `X-Combaero-Outcome` repeat the status (the
GUI warns when the exported solve failed), and any text cell that starts
with `= + - @`, a tab or a carriage return gets a leading `'` so a spreadsheet
does not run it as a formula. Unit labels are ASCII (`kg/m^3`, `Pa*s`).

Every quantity column carries its unit, `<key> [<unit>]`, dimensionless ones
included (`entrainment_ratio [-]`, `port_0_P [Pa]`, `rows_T_hot [K]`); only
the meta columns (`id`, `label`, `type`, `kind`, the solve status, ...) and the
text-valued keys (`area_source`, `second_law`, `tap`, ...) are bare. The units
come from `gui/backend/units.py` (`UNIT_MAP`, `unit_of`). Until 0.8.0 about 130
keys had no unit there and kept a bare column (#774). The results table below
keeps quantity names unit-free and carries the units in `df.attrs["units"]`.

Since 0.8.0 the frame also carries the case's inputs (#784), after every
result column, so an export reads back into the full results table on its
own:

- `input:<key> [<unit>]` on each entity's row: its node or edge data, nested
  fields dotted (`input:L [m]`, `input:composition.source`,
  `input:layers.0.thickness [m]`), plus `input:type`. Every numeric input
  has a unit (`INPUT_UNIT_MAP` in `gui/backend/units.py`); text is bare.
- `initial_guess:<key>` on the entity's row.
- one row of `type` `network` (id `_network`): the network-wide settings as
  `input:` columns and the solver options (and warm start) as `solver:<key>`.
- a row of `type` `edge` for an edge with inputs and no result row; a flow
  edge's inputs go on its auto-linked element row.

Rows the builder made (junctions, port chambers, auto links, an array's rows
and segments) carry a derived `kind` and `label` and a `parent` column naming
the drawn node or edge they belong to (#783). List-valued cells (`rows_*`,
`T_interface`) are lists in the frame and strict JSON arrays in the CSV
(`null` for a non-finite value): `pd.read_csv(path, converters={col:
json.loads})` reads them.

The result rows and their `<key> [<unit>]` columns are unchanged; a reader
that iterates over all rows now also meets the `network` row (filter on
`type`). The full format is in [GUI_TECHNICAL.md](GUI_TECHNICAL.md#csv-export).

### Results table: one row per case

`sweep()`, `load_results()` and `append_case()` (`gui.backend.results`, #784)
share one shape: one row per case, columns a 3-level MultiIndex
`(kind, entity, quantity)` with kind `inputs` or `outputs`.

```
kind          inputs                        outputs
entity        inlet          chan           _solve   chan    outlet
quantity      m_dot   Tt     L     D   ...  success  m_dot   T      ...
0             0.5     300.0  1.0   0.1      True     0.5     300.0
1             1.0     300.0  1.0   0.1      True     1.0     300.0
```

- **Rows**: labelled `0, 1, 2...`; `sweep()` uses the index of `params`
  (set it to names for readable rows), `append_case(..., name=...)` names one.
- **`df["outputs"]`** is the plain results table, `(entity, quantity)` columns:
  `df["outputs"]["chan", "m_dot"]`. Its `_solve` entity holds the status
  (`success`, `converged`, `consistent`, `outcome`, `final_norm`).
- **`df["inputs"]`** is what defined the case, as the network file gives it
  (a field left at its default is absent): boundary conditions, geometry,
  element and mixture settings (`composition.source`, ...), each entity's
  `type`, and the network-wide physics settings under `_network`
  (`global_regime`, `Nu_multiplier`, ambient conditions, ...).
- **Not inputs**: initial guesses, the warm start and the solver options
  (`method`, `init_strategy`, `timeout`, `incompressible_mach_limit`) do not
  change a converged answer. They are in `df.attrs["initial_guesses"]`,
  `df.attrs["warm_start"]` (the case a sweep row started from) and
  `df.attrs["solver_options"]`, each keyed by case label.
- **Units**: `df.attrs["units"][("outputs", "chan", "m_dot")] == "kg/s"`
  (ASCII, as in the export). A table re-read from its own `df.to_csv()` gets
  them back from the quantity names, but not the other attrs: keep the long
  export for full persistence.
- **Entities** are named by label, falling back to the id. Two entities with
  the same name are refused with a message naming both ids. `_solve` and
  `_network` are reserved.

```python
import pandas as pd
from gui.backend.results import append_case, load_results, varied_inputs

params = pd.DataFrame({"fuel_inlet.m_dot": [0.020, 0.025, 0.030, 0.035]})
df = runner.sweep(params)                      # every quantity of every entity
df = runner.sweep(params, metrics=["combustor.Tt", "combustor.phi"])  # only these outputs
T_exit = df["outputs"]["combustor", "Tt"]
ok = df["outputs"]["_solve", "success"]

varied_inputs(df)                              # the inputs that differ between rows

# One more case: inputs are read from the network the result was solved on
df = append_case(df, runner.solve({"fuel_inlet.m_dot": 0.04}), name="rich")

# A GUI CSV export reads back in full: inputs, outputs, units, initial
# guesses and solver options
one = load_results("combaero_results.csv")

df.to_csv("study.csv")                         # and back: load_results("study.csv")
```

| Function | Returns | Description |
|---|---|---|
| `load_results(source)` | `DataFrame` | A GUI CSV export, a table saved with `to_csv()`, a `to_dataframe()` frame, a `NetworkResult`, a table of this shape, or a 0.7 full-detail sweep frame, as the table |
| `append_case(df, source, name=None)` | `DataFrame` | `df` with a case appended (a `NetworkResult` reads its full inputs itself); refuses a name already used |
| `varied_inputs(df)` | `DataFrame` | The `inputs` columns whose values differ between rows |
| `NetworkRunner.sweep(params, metrics=None, ...)` | `DataFrame` | One case per row of `params`; a row whose solve raised keeps its inputs, `_solve.outcome == "exception"` |

Until 0.8.0 `sweep()` returned a flat frame: with `metrics`, the override
columns, `success`, `final_norm` and one column per metric string; without,
`to_dataframe()` stacked per case with a `_sweep_index` column.
`load_results()` still reads the latter.

### Injecting labels for unlabelled networks

GUI exports may omit node labels.  Inject them before constructing the runner:

```python
import json

with open("network.json") as f:
    schema = json.load(f)

label_map = {
    "node_1778698958490": "air_inlet",
    "node_1778698961265": "fuel_inlet",
    "node_1778699014636": "combustor",
}
for node in schema["nodes"]:
    if node["id"] in label_map:
        node["data"]["label"] = label_map[node["id"]]

runner = NetworkRunner.from_dict(schema)
```

### From a script to the GUI (#902)

`network_to_gui` turns a `FlowNetwork` built in Python into a GUI network
file. The GUI then opens it, and it solves to the same results.

```python
import json
from gui.backend.export import network_to_gui, UnsupportedForGUIError

doc = network_to_gui(
    net,
    layout="rows",                     # or {id: (x, y)}, or a callable id -> (x, y)
    labels={"burner": "Burner"},       # optional card labels
)
with open("my_network.json", "w") as f:
    json.dump(doc, f, indent=2)
```

- **Supported (phase 1):**
  - nodes: `PressureBoundary`, `MassFlowBoundary`, `PlenumNode`,
    `MomentumChamberNode` (side streams on its side inlet), `CombustorNode`;
  - elements: `ChannelElement`, `OrificeElement`, `PressureLossElement` (all
    four correlations), `LosslessConnectionElement`;
  - `ThermalWall`;
  - smooth convective surfaces only.
- **Compositions** are written as exact mass fractions.
- **Anything without a GUI form** raises `UnsupportedForGUIError`, which lists
  every object and why; nothing is dropped. Examples:
  - a ribbed or pin-fin surface;
  - a channel with an explicit convective area (the GUI derives it);
  - a head-loss area different from the element's area (the GUI has one);
  - a static-coupled `PressureBoundary`;
  - an element type not listed above.
- `layout="rows"` places each node in a column by its flow depth.
- Handles follow each node type (a plenum's `flow-left` / `flow-right`, a
  lossless connection's id-less handles, a merge chamber's `side-target`), so
  the canvas draws every edge (#901).

### NetworkResult attributes

| Attribute | Type | Description |
|---|---|---|
| `success` | `bool` | True when solver converged |
| `message` | `str` | Human-readable solver status |
| `final_norm` | `float \| None` | Residual norm at convergence |

| Method | Returns | Description |
|---|---|---|
| `get(key)` | `float` | Scalar by `<id>.<qty>` or `<label>.<qty>` |
| `node_state(label)` | `dict` | Full state dict for a node |
| `to_dataframe()` | `DataFrame` | Unit-annotated full-detail export |
| `swap_boundary(label, new_type)` | `NetworkRunner` | Retype a boundary node (mass ↔ pressure) and return a new runner pre-seeded with the full solved state |

### Boundary condition swap

```python
result = runner.solve({"air_inlet.m_dot": 1.0})

# Re-run the same point with a pressure BC (Pt auto-populated from solved state)
p_runner = result.swap_boundary("air_inlet", "pressure_boundary")
result2 = p_runner.solve()  # starts from solved state, converges immediately

# Sweep inlet pressure around the design point
import pandas as pd, numpy as np
design_Pt = result.node_state("air_inlet")["Pt"]
params = pd.DataFrame({"air_inlet.Pt": design_Pt * np.linspace(0.90, 1.10, 11)})
sweep_df = p_runner.sweep(params, metrics=["combustor.Tt"])
```

Only `mass_boundary` ↔ `pressure_boundary` swaps are supported.  The returned
runner carries the full solved state (pressures, flows, junction states) as a
warm start so `solve()` requires no additional `init_strategy`.

---

## Geometry & Materials

### Geometry Classes
```python
# Basic geometry
tube = cb.Tube(L=1.0, D=0.1)
annulus = cb.Annulus(L=0.5, D_inner=0.1, D_outer=0.2)   # DIAMETERS [m]
annulus.area(), annulus.gap(), annulus.D_mean(), annulus.circumference()
# Methods raise ValueError unless D_outer > D_inner >= 0 and the sizes they
# read are finite and positive (a swapped annulus used to give a negative area)

# Can-annular flow geometry (primary zone, transition, annular section)
can_geo = cb.CanAnnularFlowGeometry(
    L=0.3, D_inner=0.4, D_outer=0.6, L_primary=0.2, D_primary=0.1, L_transition=0.05
)
```

### Material Properties
```python
# List available materials
materials = cb.list_materials()

# Thermal conductivity [W/(m·K)]
k_ss = cb.k_stainless_steel_316  # Gen3 CSP measured, 533-1093 K
k_310 = cb.k_stainless_steel_310  # Gen3 CSP measured, 533-1073 K
k_304 = cb.k_stainless_steel_304  # CINDAS Report 45 Table 2 (Ho & Chu 1977), 200-1672 K (solidus)
k_inconel = cb.k_inconel718
k_haynes = cb.k_haynes230     # HAYNES 230 brochure alpha*cp*rho, 298-1273 K
k_hx = cb.k_hastelloy_x       # HASTELLOY X brochure + data sheet, 300-1273 K
k_188 = cb.k_haynes188        # HAYNES 188 brochure table, 300-1273 K
k_263 = cb.k_haynes263        # HAYNES 263 brochure table, 300-1273 K
k_25 = cb.k_haynes25          # HAYNES 25 (L-605) brochure table, 298-1273 K
k_625 = cb.k_inconel625       # Special Metals INCONEL 625 Table 3, 116-1255 K (Gen3 CSP 14-25% higher, #881)
k_75 = cb.k_nimonic75         # Virgamet N75 table (RT = Special Metals), 293-1273 K
k_cu = cb.k_copper            # high-purity copper, NSRDS-NBS 8 Table 2b, 100-1200 K
k_al2o3 = cb.k_alumina        # dense Al2O3, NSRDS-NBS 8 Table 14, 100-2000 K
k_1_4923 = cb.k_x22crmov12_1  # X22CrMoV12-1 data sheet, 20 and 650 C, 293-923 K
k_bc = cb.k_bondcoat_nicocraly  # NiCoCrAlY bond coat (Zhu & Miller 1999), 848-1373 K
k_tbc = cb.k_tbc_ysz
# Out-of-range or NaN T raises ValueError. get_material_conductivity(name, T)
# answers any finite T > 0 (solver robustness) and raises for NaN or inf T.

# Outside its valid range a k_* function raises; extrapolate=True opts into a
# linear continuation of the end slope instead (unsourced). The database lookup
# (and WallLayer) follows the k(T) envelope by default (#882): the end slope
# ramped linearly to 0 over 100 K, then flat, from 200 K (or T_min, if lower) to
# the solidus (1533.15 K where none is printed; T_max for alumina and ysz),
# held beyond. Not sourced: material_in_range says which part is data, and a
# wall reports "<wall>.material_in_range" (0 when any layer is outside).
k = cb.k_hastelloy_x(1400.0, extrapolate=True)
k = cb.get_material_conductivity("hastelloy_x", 1400.0)                    # envelope: 29.70
k = cb.get_material_conductivity("hastelloy_x", 1400.0, extrapolate=True)  # end slope: 31.24
cb.material_in_range("hastelloy_x", 1400.0)                                 # False
# k and its analytic dk/dT [W/(m*K^2)] on the same paths (#932): the table
# segment's slope in range, the envelope's fade outside, 0 where it holds.
k, dk_dT = cb.get_material_conductivity_and_slope("ss310", 900.0)

# The record: sourced range, source, and the solidus [K] its source prints
# (NaN where none does), the hot end of the #882 envelope.
info = cb.get_material_info("haynes230")  # T_min, T_max, T_solidus = 1574.15, solidus_source, ...
# The hot end for a metal whose source prints no solidus [K]: the lowest
# sourced solidus in the database. An effusion plate's limit without a material.
cb.GENERIC_METAL_SOLIDUS  # 1533.15
```

### Advanced Geometry Functions
```python
# Annulus: pi (Do^2 - Di^2) / 4, Dh = Do - Di (diameters, D_outer > D_inner >= 0)
area_ann = cb.annular_area(D_outer=0.2, D_inner=0.1)
dh_ann = cb.hydraulic_diameter_annulus(D_outer=0.2, D_inner=0.1)

# Residence times
t_res = cb.residence_time(V=0.01, Q=0.1)
t_res_ann = cb.residence_time_annulus(annulus, Q=0.1)
t_res_can = cb.residence_time_can_annular(can_geo, Q=0.1)
```

### Can-annular and annular-duct acoustics
```python
# N cans on an annular plenum, Bloch-Floquet: modes per azimuthal order m = 0..N/2
geom = cb.to_acoustic_geometry(can_geo, n_cans=24, L_can=0.5, D_can=0.1)
modes = cb.can_annular_eigenmodes(geom, c_can=500, c_plenum=550, rho_can=1.2,
                                  rho_plenum=1.0, f_max=2000,
                                  bc_can_top=cb.BoundaryCondition.Closed)

# Uniform thin annular duct, rigid at x = 0; bc_ends is the far end (Open =
# closed-open quarter-wave modes). Closed ends include n_axial = 0 azimuthal modes.
modes = cb.annular_duct_eigenmodes(annulus, c=500, rho=1.2, f_max=1000)
modes = cb.annular_duct_modes_analytical(annulus, c=500, f_max=1000)
# cb.AnnularDuctGeometry is a deprecated alias of cb.Annulus
```

### Rib Correlations (parametrised)

Rib correlations are **data, not code**. A parameter set carries the
coefficients, their normalisers, the validity band, the stated accuracy, and
where all of it came from -- so a built-in set and a set you tuned on your own
rig are structurally distinguishable.

```python
import combaero as cb

s = cb.han_1988_orthogonal()          # Han (1988), 90 deg orthogonal ribs
g = cb.RibGeometry(e_D=0.047, p_e=10.0, W_H=1.0, alpha_deg=90.0)

r = cb.evaluate_rib(s, g, Re=10_000)
r.R, r.f, r.e_plus, r.G, r.St_r       # 3.2000, 0.04576, 71.1, 12.21, 0.00968
r.extrapolated                        # outside the set's advisory validity
```

**Two named sets, chosen explicitly.** They cover disjoint Reynolds-number
regimes reported by different papers, not one superseding the other -- so
there is no auto-switching between them:

```python
s90 = cb.han_1988_orthogonal()        # 90 deg, Re 10,000-60,000
s45 = cb.rallabandi_2009_high_re()    # 45 deg sharp ribs, Re 30,000-400,000

g = cb.RibGeometry(e_D=0.14, p_e=7.5, W_H=1.0, alpha_deg=45.0)
r = cb.evaluate_rib(s45, g, Re=100_000)
r.G                                    # 35.40 -- see han_ribbed_high_re.md
```

`rallabandi_2009_high_re` is 45 deg and sharp-edged ribs only: the source
reports round-edged ribs at the same conditions instead following Han's
correlation, but calls the agreement "coincidental" rather than physically
grounded, so it is not wired as an automatic edge-profile switch.

**A third set covers angled ribs generally**, `han_park_1988_angled` (Eq.
4.17/4.18, `alpha` 30-90 deg, `W/H` 1-4):

```python
angled = cb.han_park_1988_angled()
g = cb.RibGeometry(e_D=0.06, p_e=15.0, W_H=2.0, alpha_deg=60.0)
r = cb.evaluate_rib(angled, g, Re=30_000)
r.R, r.G                              # 3.2332, 18.4408
```

Its `R` and `G` carry genuine, unsmoothed discontinuities the source
states rather than a numerical artefact: `R`'s `(W/H)^m` term switches at
`alpha == 90 deg` (up to 62.5% at `W/H = 4`), and `G`'s `alpha`/`p_e`
exponents switch on whether the channel is square (`W/H == 1`, up to 27%
at low `alpha`). Neither the printed equations nor the extraction record a
smooth transition between the branches, so none is invented -- a caller
whose solver traverses either boundary exactly needs its own guard.

`RAlphaShape` and `GShapeModel` (both importable from `combaero`) are what
make this representable as data: setting `R_alpha_shape` to
`QuadraticAlpha` and populating `R_quad_c0/c1/c2` plus the two
`R_quad_WH_exponent_*` fields lets a user-supplied set describe its own
angle-dependent `R`, the same way `RibTerm`'s normaliser lets one describe
its own geometry dependence.

**A fourth set extends the same family to narrow channels**,
`han_1989_narrow_channel` (Eq. 4.19, `alpha` 30-90 deg, `W/H` 1/4-1):

```python
narrow = cb.han_1989_narrow_channel()
g = cb.RibGeometry(e_D=0.0625, p_e=15.0, W_H=0.3, alpha_deg=45.0)
r = cb.evaluate_rib(narrow, g, Re=30_000)
r.R, r.G
```

Its source has **no printed equation for `R`** below `W/H = 1` -- only a
drawn figure line -- so `R_alpha_shape = QuadraticAlphaTwoBand`'s
coefficients (`R_quad_c0/c1/c2` for `1/2 <= W/H < 1`,
`R_quad_narrow_c0/c1/c2` for `1/4 < W/H < 1/2`, switched by
`R_WH_band_boundary`) are this project's own fit to that line, and the
set's `provenance` is `Fitted` rather than `Extracted` -- see
`validation/cooling/extractions/han_ribbed.md`, item 40. `G`
(`GShapeModel.NarrowChannelAlphaSwitch`) *is* text-extracted: a genuine
~20% jump in its leading constant at `alpha == 90 deg`
(`G_narrow_C_alpha90` vs `G_narrow_C_off_axis`), plus a further `(W/H)`
correction to both the constant and the `e+` exponent itself below
`G_narrow_WH_band_boundary`.

**Rib shape, and nine configurations from one paper.** A set declares the
`RibShape` it was fitted to (`Transverse`, `Parallel`, `Crossed`, `V`,
`Lambda`; `Unspecified` binds nothing) because V, crossed and Lambda ribs at
the same angle are different geometries. `han_zhang_lee_1991` returns Han,
Zhang & Lee (1991)'s Table 2 row for a tested configuration -- one rig,
square channel, `e/D` 0.0625, `P/e` 10 -- and raises for anything else:

```python
vee = cb.han_zhang_lee_1991(cb.RibShape.V, 60.0)
g = cb.RibGeometry(e_D=0.0625, p_e=10.0, W_H=1.0, alpha_deg=60.0)
r = cb.evaluate_rib(vee, g, Re=30_000)
r.G, r.G_bar, r.has_G_bar           # its own printed G_bar, not 1.2 G
vee.symmetric                       # False: a V reversed is a Lambda
```

**Ratio-form sets** (`RibRatioSet`, #444) carry `Nu/Nu0` and `f/f0` as
multipliers on the smooth-duct baseline the source normalised by --
`Nu0_source` is required, because a ratio without its reference is not a
number. In the fitted range they reproduce the paper; below `Re_floor` they
hand over smoothly to Gnielinski by default, so Nu keeps a laminar floor as
Re -> 0 rather than collapsing with Dittus-Boelter:

```python
s = cb.RibRatioSet()
s.name, s.source = "my_rig", "measured 2026-10"
s.C_Nu, s.Nu_Re = 2.5, cb.RibTerm(-0.2, 1.0e4)       # Nu/Nu0 = 2.5 (Re/1e4)^-0.2
s.C_f, s.f_Re = 8.0, cb.RibTerm(0.1, 1.0e4)           # f/f0  = 8 (Re/1e4)^0.1
s.Nu0_source = cb.RatioBaseline(0.023, 0.8, 0.4)     # Dittus-Boelter
s.f0_source = cb.RatioBaseline(0.046, -0.2)
s.Re_floor = 1.0e4
cb.validate_rib_ratio_set(s)
r = cb.evaluate_rib_ratio(s, g, Re=30_000, Pr=0.7)
r.Nu, r.dNu_dRe, r.f, r.df_dRe, r.below_floor
```

Shipped: `cb.taslim_spring_1987(aspect_ratio_taslim, e_D)`, one set per
two-side configuration Taslim & Spring (1987) tested -- AR 0.5 (e/D 0.125,
0.250), 1.0 (0.083, 0.167), 3.5 (0.053, 0.107, 0.161), transverse ribs,
p/e 10; anything else raises. Their AR is height/width, so `W/H = 1/AR`.
Friction is a constant Fanning `f` over its own measured range
(`Re_floor_f`, which may differ from Nu's `Re_floor`):

```python
t = cb.taslim_spring_1987(1.0, 0.083)
g1 = cb.RibGeometry(e_D=0.083, p_e=10.0, W_H=1.0, alpha_deg=90.0)
cb.evaluate_rib_ratio(t, g1, Re=50_000, Pr=0.7).Nu
```

`RibRatioOptions.below_floor` picks the handover (`Gnielinski`,
`SourceBaseline`, or `User` with your own `user_Nu0`). It only acts below the
fitted range: inside it, Nu is the paper's whatever you choose -- matching a
specific rig across the range is what `Nu_multiplier` is for.

**Supply your own.** Real hardware needs it -- no published correlation is
precise enough for a specific rig:

```python
mine = cb.han_1988_orthogonal()
mine.name = "rig_3"
mine.source = "measured 2026-09"
mine.provenance = cb.RibProvenance.User   # not Extracted -- the claim differs
mine.C_G = 4.1
cb.validate_rib_set(mine)                 # rejects malformed sets, loudly
```

Three things worth knowing:

- **The normaliser is data.** `RibTerm(exponent, reference)` divides by
  `reference` before applying the exponent. `3.2 (p/e/10)^0.35` and
  `1.4294 (p/e)^0.35` are the same function; using one constant under the
  other convention is wrong by `10^0.35 = 2.24x`, uniformly, which never looks
  like a trend.
- **Validity is advisory.** `r.extrapolated` reports; nothing refuses. A band
  belongs to the source's rig, not to yours.
- **`St_r` is the ribbed side.** Combining it with the smooth walls is the
  caller's job, and depends on how many walls are ribbed.

Bad *parameters* raise from `validate_rib_set`. Bad *operating points* never
raise: reverse flow, zero flow and extreme values are guarded smoothly, because
a solver probes states that are not physical and a throw inside a residual
kills the solve.

### Ribbed Channels

```python
from combaero.network import ChannelElement, ConvectiveSurface, RibbedModel

surface = ConvectiveSurface(
    area=0.1,
    model=RibbedModel(
        e_D=0.06, p_e=10.0, alpha_deg=90.0,
        W_H=1.0,                  # aspect ratio: a Dh does not determine it
        n_ribbed_walls=2,         # 2 = two OPPOSITE walls, as Han measured
    ),
)
```

`n_ribbed_walls` accepts 1, 2 or 4. Three is rejected -- it has no unambiguous
geometry.

**Two asymmetries, both from the source rather than convenience:**

- **Friction needs no wall weighting.** The correlation's `f` is already the
  four-sided channel value, so the element uses it directly. Nothing multiplies
  pipe friction.
- **Heat transfer does.** The correlation gives the ribbed side; the smooth
  walls come from the base correlation and the channel average is the
  area-weighted combination. The result exposes `h_ribbed` and `h_smooth`
  separately, because the average hides a modelling choice worth seeing.

**A documented gap.** Plain smooth walls give `h_s/h_r` around 0.42 where Han's
own channel average implies 0.70 -- ribs enhance the adjacent smooth wall by
10-50% too, which no correlation here covers. The channel average is therefore
about **20% below** Han's measurement for the two-ribbed-wall square case.

```python
model = RibbedModel(..., smooth_wall_Nu_multiplier=1.67)   # reproduces Han
```

That knob exists because the ribbed side already has a better one: change `C_G`
on the parameter set, which records what you changed and why. The smooth walls
come from Gnielinski and have no set of their own.

### Pin-Fin Channels

`PinFinModel` puts a pin-fin array on a `ChannelElement`'s `ConvectiveSurface`
(#335). It uses the pin-fin sets below: heat transfer from `nu_set`
(default `metzger_1986_staggered_nu()`), the drop from `f_set` (default
`metzger_1982_staggered_friction()`), optionally transferred by a
`modifier`.

```python
from combaero.network import ChannelElement, ConvectiveSurface, PinFinModel

pins = PinFinModel(pin_diameter=0.005, S_D=2.5, X_D=2.5, H_D=1.0, N_rows=10,
                   k_pin=20.0)     # W/(m K), for the fin efficiency
ch = ChannelElement("te", "A", "B", length=0.125, diameter=0.009,
                    surface=ConvectiveSurface(area=0.01, model=pins))
```

- **The element's flow area is the unobstructed channel.** Vmax follows
  exactly from the array geometry (`A_min/A_frontal`), and
  `Re_D = rho Vmax D / mu`.
- **`ConvectiveSurface.area` is the endwall's base (planform) area.** The
  returned `h` is referenced to it and includes the fin efficiency of
  pins fed from each wall to mid-height:
  `h_eff = h (A_endwall_exposed + eta_fin A_pin) / A_base`. `h_array` and
  `eta_fin` on the result show both halves.
- **The drop is the correlation's own per-row form,**
  `dP = 2 rho Vmax^2 N f(Re_D)`, with an analytic Jacobian in mass flow,
  upstream T and P (through rho and mu). The code removed in 0.7.0 used
  `rho V^2/2` here, a factor of 4 low.
- **Inline arrays** default to `chyu_1990_friction(Inline)`, the only
  inline friction source in hand (one geometry: H/D 1, S/D = X/D = 2.5).
  Staggered friction is never substituted: a staggered `f_set` on an inline
  array raises. Heat transfer defaults to Chyu's inline set, or use
  `modifier=cb.chyu_1998_inline_over_staggered()` to transfer a staggered
  set.
- **Known gaps:**
  - `T_aw` is borrowed from the smooth correlation, as the impingement
    paths do.
  - `dh_dT` is the smooth result's relative property sensitivity scaled to
    `h`; the pin correlations expose no temperature derivative of their
    own.

### Impingement-Cooled Walls (jet-plate arrays)

A jet plate fed from a plenum, impinging on a target wall, the spent air
leaving through the gap as crossflow. **This is one configuration**, and
`ImpingementArray` builds it as one object (#465):

```python
from combaero.network import ImpingementArray, WallLayer

arr = ImpingementArray("ia", n_rows=10, d_jet=0.00254, xn_d=5.0, yn_d=4.0, z_d=2.0,
                       span=0.122, plate_thickness=0.00254)
arr.add_to(net, supply="supply", exit="exit")            # rows ia__p{i}, ia__x{i}, ia__c{i}
arr.add_wall(net, "w", hot_element="hot", layers=[WallLayer(0.001, 20.0)],
             radiation_hot=None,       # GasRadiation / CasingRadiation, hot face
             radiation_target=None)    # e.g. CasingRadiation toward the jet plate
res = NetworkSolver(net).solve()
arr.summarize(res["__element_diag__"])   # m_dot, dP, Re_j/Gc_Gj ranges, rows_Nu, ...
```

`add_wall` puts one `ThermalWall` per row on that row's target footprint,
each with the same `radiation_hot` / `radiation_target` mapped to its
`radiation_a` / `radiation_b` (#1071); the
hot side is evaluated once per row at its own state, so a hot side that changes
along the array needs its own segments. Heating the crossflow lowers its
density, which raises the jets' momentum cost and the supply non-uniformity
beyond Florschuetz's isothermal model (1.88 vs 1.61 on a 10-row (5,4,2) array
against 1200 K gas) -- the direction is physical; the isothermal source cannot
score it.

**Initial (bypass) crossflow (#467).** A stream entering the gap's upstream
end, ahead of row 1, is declared by passing its node:

```python
arr = ImpingementArray("ia", ..., bypass_length=None, bypass_entry_K=0.0)
arr.add_to(net, supply="supply", exit="exit", bypass="byp")   # adds ia__x0: byp -> ia__c1
```

- The lead-in `ia__x0` is an `ImpingementCrossflowElement(bypass_inlet=True)`,
  `bypass_length` long (None: `xn/2`, where Florschuetz, Metzger and Su's
  array begins). Row 1's `Gc/Gj` is then the bypass over its jets.
- **What feeds it decides the physics.** Another array's last gap segment
  (arrays in series, `a.add_to(..., exit="joint")`,
  `b.add_to(..., bypass="joint")`): the stream continues with its
  momentum, and the gaps must match. Anything else: the bypass is
  accelerated from the node's total pressure, so a feed duct's dynamic head
  arrives with it, through a plenum or a momentum chamber alike;
  `bypass_entry_K` (0.5) is for an abrupt entry from a large plenum.
- **Flow split checked against the authors' own data**, NASA CR-3630
  Table D.1: 42 tests, 11 geometries, mc/mj 0 to 1
  (`python/tests/test_impingement_bypass_validation.py`). Worst row per
  test in Gj/Gj_mean: z/d 2-3 median 0.021, max 0.100; z/d 1 up to 0.27,
  where the authors' own constant-Cd model misses by up to 1.2 (the hole Cd
  falls with crossflow there, and wall shear matters). P/P0 mean error
  0.007.
- **Heat transfer:** Florschuetz et al. (1981) at the local Gc/Gj, which is
  CR-3630's own usage. A bypass at another temperature moves each row's
  reference temperature (Florschuetz and Su 1985, NASA CR-3936):
  `T_aw = T_j + eta_r phi_0 (T_0 - T_j)`, with `phi_0` the bypass share of
  the crossflow approaching the row and `eta_r` from `cb.impingement_eta_r`,
  a closure fitted to CR-3936's inline data (RMSE 0.083 in eta_r leaving
  one geometry out). Only the bypass's share: the wall's own heating of the
  spent air is inside Florschuetz's h already. The summary reports
  `m_dot_bypass`, `bypass_to_jet_flow`, `bypass_Tt_minus_supply` and the
  rows' `T_aw`, `eta_r` and `bypass_fraction`, and sets
  `heat_transfer_outside_source` when a bypass more than 1 K off the
  supply meets a staggered array or a row outside the inline data.
- **Starved rows are reported unresolved (#729).** A strong bypass leaves
  the upstream rows with almost no jet flow. Their plate drop is then a few
  pascals against a gap dynamic head of kilopascals, below what the
  crossflow chain resolves: their flows come out alternating in sign from
  row to row. `rows_unresolved` (1.0 per row) and `n_rows_unresolved` mark
  a row whose plate total-pressure drop is below 1% of the gap's dynamic
  head (`UNRESOLVED_ROW_DPT_FRACTION`). Read its flow, `Nu` and `eta_r` as
  unresolved: the surface there is cooled by the gap's channel flow, which
  the jet correlation does not model. On Table D.2 the flag marks exactly
  the 8 rows the jet Nu misses by -94%, and no other row.

Spent air leaving through the target (#468) is a different configuration,
with its own element to come. The building blocks below are public for
hand-built networks, one pair per spanwise row:

```python
from combaero.network import ImpingementCrossflowElement, ImpingementPlateElement

d, span = 0.00254, 0.122
for i in range(1, 11):
    net.add_element(ImpingementPlateElement(
        f"p{i}", "supply", f"c{i}", d_jet=d, xn_d=5.0, yn_d=4.0, z_d=2.0,
        span=span, plate_thickness=d, row=i,   # row: diagnostics only
    ))
    net.add_element(ImpingementCrossflowElement(
        f"x{i}", f"c{i}", f"c{i+1}" if i < 10 else "exit",
        length=5.0 * d, height=2.0 * d, span=span,
    ))
```

```
supply plenum ---+--------------+--------------+
                 |              |              |
            [Plate 1]      [Plate 2]      [Plate 3]
                 |              |              |
crossflow:      c1 --[Cross]-- c2 --[Cross]-- c3 --[Cross]-- exit
```

**`ImpingementPlateElement`** is an orifice (`n_holes = round(span/(yn d))`
holes in parallel) whose flow is the row's jet flow, plus Florschuetz,
Truman and Metzger's (1981) heat transfer on the target footprint
`n_holes xn yn d^2`. Its **Gc/Gj comes from the network**, not the
uniform-supply closed form:

    Gc/Gj = (m_c / m_j) (pi/4) / ((yn/d)(z/d))

where `m_j` is the plate's own flow and `m_c` every other inflow to its
`to_node`, i.e. the crossflow approaching the row. A non-uniform supply, a row
of different geometry or a bleed therefore shows up in the heat transfer.
`h` depends on a neighbour's flow, so the element returns `dh_dsources` and
the solver relays it into the Jacobian (the `network_flow_inputs()` opt-in).
`T_aw` is the supply plenum temperature, Florschuetz's own reference, moved
by a bypass at another temperature through CR-3936's eta_r (see
`ImpingementArray` above): the element reads the bypass node's temperature
(`network_temperature_inputs()`) and the lead-in's flow, and returns
`dT_aw_dsources` / `dT_aw_dT_nodes` for the solver to relay.
The physics is C++'s `cb.jet_row_heat_transfer(set, cb.JetRowGeometry(...),
m_jet, m_crossflow, mu, k, Pr)`, which returns `h` and its analytic
derivatives with respect to both flows; the element only decides which flows
are crossflow.
- Hole Cd: default `'fixed'` at Florschuetz's 0.79 (measured 0.73-0.85), what an unset `Cd` takes; a `Cd` given with a correlation raises (#872);
  alternatives `'IdelchikThick'`, `'Lichtarowicz'` (long hole, l/d 2-10)
  and `'McGreehanSchotsch'` (supply-side `U1/Vi = 0`, never Gc/Gj).
  `Cd_in_range` reports whether the plate is inside the chosen source's range.
- **Crossflow on the holes (#800), `cd_crossflow`** (also on `ImpingementArray`
  and in the GUI). The gap's crossflow bends the jets over and lowers their Cd,
  steeply at large Gc/Gj. In q = G_ideal / G_c (ideal hole over crossflow mass
  velocity, G_ideal = sqrt(2 rho dP)):
  - `'rogers_hersh_1976'` (default): Rogers & Hersh (1976, NASA CR-2681) Fig. 11,
    0.97 (Gj/Gc)^0.534 injecting, 0.80 (Gj/Gc)^0.47 for a reversed jet,
    resolved explicitly in q and capped (C1 over +/-10%) at the hole's own
    Cd. Below Gc/Gj ~1.5 at Cd 0.79 the hole's own Cd holds: Florschuetz's
    arrays (Gc/Gj <= 0.77) are bit-identical. Cross-source for any plate.
  - `'florschuetz_metzger_su_1982'`: CR-3630 Table 4.1 / Eq. 4.2, the
    authors' own law for B(5,4,1-3), B(5,8,1) and (inferred) B(10,x,1);
    refused for other plates. The law behind Table D.1's Gj, so the one to
    score D.1 like for like. Its own Cd replaces the hole's.
  - `'none'`: the hole's own Cd at any crossflow (Eq. 8 / Eq. 4.10's model).

  Resolved in q, the row flow goes as dP^1.07 under a strong crossflow (a
  linear resistance) instead of sqrt(dP), which removed the alternating-sign
  stagnant rows of CR-3630 B(5,4,1) (#799). Diagnostics add `Cd` (the one
  used) beside `Cd_own`. A cold solve that fails is retried once from the
  constant-Cd root (`NetworkSolver`, "[from the constant-Cd solution]"). With
  a `'florschuetz_metzger_su_1982'` plate that stage runs first (#1021): the
  law's root has a narrow basin the cold seed misses.
- Diagnostics: `n_holes`, `Re_j`, `Gc_Gj`, `Nu`, `St`, `htc`, `T_aw`, `eta_r`,
  `bypass_fraction`, `eta_r_extrapolated`, `surface_extrapolated`, and
  `Gc_Gj_closed_form` (Eq. 8) when `row` is given.
- `cb.impingement_eta_r(Gc_Gj, xn_d, yn_d, z_d)` returns `eta`,
  `deta_dGc_Gj` and `extrapolated`; the constants are
  `cb.FLORSCHUETZ_SU_1985_ETA_A/R/N`.
- Not modelled: the temperature sensitivity of the correlation's properties
  (`dh_dT = 0`) and the self-heated crossflow's own temperature beyond what
  Florschuetz's h holds.

**`ImpingementCrossflowElement`** is a `ChannelElement` (area `height *
span`) plus the momentum the jets cost: they enter with no streamwise
momentum, so each merge drops the static pressure by
`(m_b^2 - m_a^2)/(rho A^2)` (momentum flux m^2/(rho A), any flow direction,
#815). The drop is split half-and-half between
the segments either side, so each crossflow node holds the static pressure at
its row centre. **Without this term every row sees the same pressure
difference and the supply stays uniform.** That puts Gc/Gj at twice Eq. 8 on
Florschuetz's strongest-crossflow geometry. The term is C++'s
`cb.station_half_drop` with `kappa = 0`, each station at the density of its
own node (see `CrossflowSegmentElement` below).
`bypass_inlet=True, entry_K=0.0` makes it the lead-in of an initial
crossflow (#467): a station with nothing injected when another gap's segment
arrives at its `from_node` (the stream continues; the cross-sections must
match), otherwise an entry from that node's total pressure.

#### Crossflow segments: merge and bleed stations (#471)

`ImpingementCrossflowElement` is one case of `CrossflowSegmentElement`, a
channel segment whose end nodes are side-stream STATIONS. Over a station, with
the side stream's axial velocity `kappa * u_arriving`:

    dP_static = (m_b^2 - m_a^2 + kappa m_a (m_a - m_b)) / (rho A^2)

| `kappa` | station | source |
|---|---|---|
| `cb.STATION_KAPPA_MERGE_NORMAL` = 0 | jets merging normally (impingement) | Florschuetz's P + G^2/rho = const |
| `cb.STATION_KAPPA_BLEED_BASSETT` = 0.75 | bleed through the wall (effusion holes) | Bassett et al. (2001) K2/K5, Eq. 15, reproduced for every split |

**Joining jets bring their own axial momentum (#813).** Every element that
discharges INTO a station declares it, the hook a merge chamber reads:
`injection_momentum` = m w cos(theta), w the jet's velocity from its supply's
stagnation state to the station's static pressure, theta its
`injection_angle_deg` (an `EffusionPlateElement`'s `angle_deg`). It enters the
station as `- J / A`, half on each segment. Normal jets (90 deg, the default,
and every impingement row) bring none, so `kappa = 0` stations are unchanged;
`kappa` stays the receiving duct's own coefficient (Bassett's bleed). The
second-law ledger counts a hole whose jet brings momentum as part of the
chain's control volume, entering at its supply's Pt.

```python
from combaero.network import CrossflowSegmentElement

seg = CrossflowSegmentElement("s1", "b1", "b2", length=0.02, area=4e-3, Dh=0.01,
                              from_kappa=cb.STATION_KAPPA_BLEED_BASSETT, prev_seg="s0",
                              to_kappa=cb.STATION_KAPPA_BLEED_BASSETT, next_seg="s2")
first = CrossflowSegmentElement("s0", "plenum", "b1", length=0.01, area=4e-3,
                                entry_K=0.5,          # reservoir entry, sharp
                                to_kappa=cb.STATION_KAPPA_BLEED_BASSETT, next_seg="s1")
```

- **Taper (#812).** `area` is the from-end's section; `area_to` (default: the
  same) the to-end's. Each end's station half-drop, joining jets and
  static-to-total closure use its own section, friction the mean; Pt carries
  across bar friction (an isentropic area change). Incompressible, a from-end
  station adds Bernoulli's `cb.segment_area_change_drop` (m^2/(2 rho)
  (1/A_to^2 - 1/A_from^2), even in m). The incompressible regime stops
  holding where the dynamic head is a sizeable fraction of P (a 2:1 taper
  reversed at M ~0.55 has no root there; compressible solves it).
- **Centred.** A segment carries half of the station at each end, each at its
  node's own density, so every node holds its station's mid static pressure.
  Station nodes must be `PlenumNode`s (Pt = P).
- **Ends are explicit.** `prev_seg`/`next_seg` name the neighbour segments;
  `entry_K` makes `from_node` a reservoir; a segment into a plenum with no
  `to_kappa` dumps its dynamic head (static continuity). A reversed flow
  through an `entry_K` end is such an exit too: the entry term applies to
  inflow only.
- **Compressible: the Fanno duct runs between the duct's own end states
  (#802).** A station node holds a STATIC pressure, so it is not the duct's
  stagnation pressure. The duct's inlet static is `P_from` less half the
  from-station, and its outlet static is `P_to` plus half the to-station (at
  an exit, the node's P). Each end's Pt is the stagnation pressure of that
  static state at the segment's flux, from the momentum chamber closure. An
  end with no station is a terminal node. Where it feeds the segment, its own
  Pt is used: a reservoir entry also adds only the loss `K_in/(1 + K_in)` of
  `channel_entry_drop` to the outlet. Where the segment discharges into it, a
  plenum or boundary takes static continuity, and a momentum chamber takes Pt
  continuity. An end static below the sonic static P* of the flux is choked:
  its Pt is continued C1 below Pt*, never on the supersonic branch. Before #802
  the node's static P was fed as the duct's P0, a dynamic head low, so exit
  segments choked early: CR-3630 D.2 B(5,4,1) at mc/mj 1.01 runs its exit
  segment at 81% of its true choked flow, not at 99.7% of the false one. The
  diagnostics `choked`, `M_in_duct` and `M_exit_duct` report this duct; they
  used to report every exit segment choked.
- **Bleed scored on Bassett's own measured K5** (Fig. 7c, 43 points): MAE
  0.039 in K overall, 0.021 at q >= 0.6 (the few-percent bleed of an effusion
  station). Bassett measured psi = 1-3; a station bleeding through holes far
  smaller than the duct is an extrapolation in psi, which K5 is independent
  of by derivation.

**One configuration only.** These elements cover spent air leaving down the
crossflow channel. A hand-wired network describing a different
configuration is refused at set-up, naming it:
- a stream feeding the chain that is neither a jet row nor a segment
  declared `bypass_inlet=True` (#467);
- spent air leaving through the target, impingement-effusion or film (#468);
- two plates merging at one node.

**Measured against the source** (`python/tests/test_impingement_plate_validation.py`):
- **Flow model:** over the 27 Fig. 6 geometries x 10 rows, the network
  reproduces Florschuetz's own 1D flow model (Eq. 8, and the cosh jet
  distribution behind it) to a max 3.7% in Gc/Gj (bias -0.5%) and 3.8% in
  Gj/Gj_mean.
- **Fig. 6:** the chain's Nu/Nu1 scores +4.0% bias and 6.9% MAE over 242
  points. The correlation at the paper's own abscissae scores +3.7% / 6.8%.
- **Row 1 vs Fig. 5:** absolute Nu1 matches as in #461.

### Impingement Channels

```python
from combaero.network import ConvectiveSurface, ImpingementModel, SingleJetImpingementModel

surface = ConvectiveSurface(
    area=0.01,   # this ROW's own target-plate footprint
    model=ImpingementModel(
        d_jet=0.002, xn_d=8.0, yn_d=6.0, z_d=2.0,
        row=3,       # 1-indexed, counting from upstream
    ),
)
```

For a jet plate fed from a plenum use the elements above: they carry the
row's own jet flow and read Gc/Gj from the network. `ImpingementModel` on a
`ChannelElement` remains for a row whose flow you impose yourself.

**One element models one spanwise row**, not a whole array -- there is no
single "channel Nu" for a jet array the way there is for a smooth or ribbed
duct, since downstream rows see progressively more crossflow than upstream
ones. Model a real array by chaining `n_rows` separate elements, each with
its own `row`; row 1 sees zero crossflow by definition (`Gc/Gj = 0`,
Florschuetz's own `Nu1`).

**The element's mass flow is this row's jet flow; `area` counts the
holes.** `area` is the row's target-plate footprint, so the row has
`n = area / (xn_d d_jet * yn_d d_jet)` holes, each carrying `m_dot / n`, and
`Re_j = 4 (m_dot / n) / (pi d_jet mu)` -- Florschuetz's jet mass velocity on
the hole area. Until #460 the jet flow was taken as `rho v * area`, which
made Re_j scale with the target area: driven with the paper's own flow, the
element then missed Florschuetz's Fig. 5 by +268% on average; it now matches
it to 5.2%, like the correlation. Whether every row gets the same total flow (a uniform-supply
approximation) or a row-dependent one (Florschuetz's own Eq. 7, deliberately
not implemented -- see `impingement_correlation.h`'s module comment) is a
choice for whoever assembles the chain, not something this element decides.

**Default convective area.** A `ChannelElement` whose surface was given no
`area` gets one that follows the model (#462): an impingement row gets its
target footprint matched to the channel's crossflow cross-section,
`A_flow * (xn/d) / (z/d)`, so the hole count follows from the channel; a single
jet gets Goldstein's averaging disc, `pi (R_D d_jet)^2`; every other surface
keeps the wetted wall `pi D L`. The channel diagnostics report the surface's
own `Re_surface` and `surface_extrapolated`.

**Pressure drop is not modelled here.** Model the jet plate's own orifice
loss with a proper `OrificeElement` upstream, using a real
discharge-coefficient correlation -- duplicating it inside this heat-transfer
model would risk double-counting it. `f`/`dP` reported by this element are
the crossflow's own plain smooth-duct values (Han's book gives no
impingement-specific friction correlation, unlike ribs' own `R`/`f`
relationship), borrowed to compute a `T_aw` Eq. 4.2/4.3's own recovery-factor
model would otherwise have to supply.

**A single free jet** (Goldstein, Behbahani and Heppelmann 1986) is a
separate model, since there is no array and no crossflow:

```python
surface = ConvectiveSurface(
    area=0.001,   # this jet's own target patch
    model=SingleJetImpingementModel(d_jet=0.003, L_D=7.75, R_D=5.0),
)
```

`R_D` is a real modelling choice, not a detail: the correlation's `Nu` is
the AVERAGE over the disc of radius `R = R_D d_jet` (Han's Eq. 4.1), so the
convective area should be that disc.

Both models expose the same wall-coupling derivatives (`dh_dmdot`, `dh_dT`,
`dT_aw_dmdot`, `dT_aw_dT`) and an `extrapolated` flag (array only -- a single
jet has no stated validity range to be outside of), matching `RibbedModel`'s
contract.

### Jet Impingement Correlations (parametrised)

Two independent regimes, matching
`validation/cooling/extractions/han_impingement.md`'s own split, wired into
the network elements above (`ImpingementModel`, `SingleJetImpingementModel`).

**Single jet** (Goldstein, Behbahani and Heppelmann, 1986): one free round
jet on a flat plate, no crossflow, no array.

```python
import combaero as cb

s = cb.goldstein_1986_single_jet()
nu_q = cb.single_jet_impingement_nu(
    s, cb.ImpingementThermalBC.ConstantHeatFlux, Re=25_000, L_D=7.75, R_D=5.0)
nu_t = cb.single_jet_impingement_nu(
    s, cb.ImpingementThermalBC.ConstantWallTemperature, Re=25_000, L_D=7.75, R_D=5.0)
nu_q, nu_t                            # 59.93, 55.71 -- Han's own check point (60, 56)
```

`L_D = 7.75` is the optimum spacing -- structural, not a separately fitted
fact: `A` itself is the numerator's value there, for any `Re` or `R_D`.

**Jet array with crossflow** (Florschuetz, Truman and Metzger, 1981): a
staggered or inline array fed from a common plenum, where every row's Nu is
degraded by the crossflow accumulated from every row upstream of it.

```python
inline = cb.florschuetz_1981_inline()

result = cb.jet_array_impingement_nu(
    inline, Re_j=10_000, Gc_Gj=0.3, Pr=0.7, xn_d=10.0, yn_d=6.0, z_d=2.0)
result.Nu, result.extrapolated        # 27.93, False
```

`Re_j` and `Gc_Gj` are that ROW's own local values, not an array mean --
Florschuetz correlates "the individual spanwise row jet Reynolds number".
`Gc_Gj` at a row has its own closed form (the paper's own Eq. 8), needing
only geometry and a jet-plate discharge coefficient:

```python
cb.crossflow_to_jet_ratio_at_row(yn_d=8.0, z_d=2.0, C_D=cb.FLORSCHUETZ_1981_DEFAULT_CD, row=1)   # 0.0 -- Nu1's own definition
cb.crossflow_to_jet_ratio_at_row(yn_d=8.0, z_d=2.0, C_D=cb.FLORSCHUETZ_1981_DEFAULT_CD, row=10)  # 0.404
```

It depends on `(yn/d)(z/d)` only, not `xn/d` -- the source states the flow
distribution is independent of streamwise hole spacing and hole pattern.

**Pressure range (#474).** The model is incompressible. Its measured range, NASA CR-3630 Table D.1, reaches a channel pressure of 0.72 of the jet plenum's. Over that range it matches the measured distributions (rows 2-10, 11 geometries) with bias -0.0% and RMS 3.4%, the same as the compressible network chain (-0.6%, 3.4%). Below about 0.7 the two part by 4-11% and no data decide between them, so treat it as extrapolation there. The heat-transfer correlation itself held within 16% up to jet Mach ~1 on the data the paper excluded for compressibility.
`FLORSCHUETZ_1981_DEFAULT_CD` (0.79) is the paper's own recommended default
absent a measured value, and remains the default. To predict `C_D` from the
plate's own geometry instead, pass `mcgreehan_schotsch_1988_cd` (see above);
for a plenum-fed plate leave its `U1_over_Vi` at 0, since Florschuetz's
`Gc/Gj` is discharge-side crossflow and must not be fed to it:

```python
C_D = cb.mcgreehan_schotsch_1988_cd(Re=1e4, r_over_d=0.0, L_over_d=t_over_d)
cb.crossflow_to_jet_ratio_at_row(yn_d=8.0, z_d=2.0, C_D=C_D, row=10)
```

At `t/d = 1` the two agree to 0.6%, so this changes little; at `t/d = 0.5` the
correlation gives 0.688 and the predicted `Gc/Gj` rises by up to 12% at
downstream rows (about 1% in `Nu`). The ISO 5167 `Cd_sharp_thin_plate` family
in `orifice.h` remains inapplicable here -- it models a hole in a pipe run.

**Two named sets**, chosen explicitly by hole pattern -- their coefficients
differ, not just their validity:

```python
inl = cb.florschuetz_1981_inline()      # xn/d 5-15
stg = cb.florschuetz_1981_staggered()   # xn/d 5-10, genuinely tighter
```

Bad *parameters* raise from `validate_single_jet_set`/`validate_jet_array_set`.
Bad *operating points* never raise: reverse flow and out-of-range crossflow
ratios are guarded smoothly, matching the ribs' policy above.

Leading-edge/curved-surface impingement (Section 4.1.4) and the simpler,
looser forms (Eq. 4.6 Kercher-Tabakoff -- graphical, not closed-form anyway;
Eq. 4.7/4.8 -- Florschuetz's own less-tight alternate) are out of scope,
deferred per the extraction's I3.

### Pin-Fin Correlations (parametrised)

Pin-fin arrays (#335) as provenanced sets in each source's own form, with
exact geometry-only conversion to one canonical basis: `Re_D` on the pin
diameter and the velocity at the minimum flow area, `Nu_D = h D / k`, and
`f = dP / (2 rho Vmax^2 N)` (per row, so `dP = 2 rho Vmax^2 N f`). Heat
transfer and friction are separate sets, so any Nu set can pair with any
friction set.

```python
import combaero as cb

g = cb.PinFinGeometry(S_D=2.5, X_D=2.5, H_D=1.0, N_rows=10)   # staggered
nu = cb.evaluate_pin_fin_nu(cb.metzger_1986_staggered_nu(), g, Re_D=2.0e4, Pr=0.7)
f = cb.evaluate_pin_fin_friction(cb.metzger_1982_staggered_friction(), g, Re_D=2.0e4)
nu.Nu, nu.dNu_dRe, nu.extrapolated     # 91.78, 0.00317, False
f.f, f.extrapolated                    # 0.0755, False

# VanFossen works in its own D' basis; the result comes back canonical.
cb.evaluate_pin_fin_nu(cb.vanfossen_1982_staggered_nu(), g, 2.0e4, 0.7).Nu   # 90.91
```

The second call is flagged `extrapolated`: VanFossen's arrays had 4 rows, not
10. A set outside its validity box still evaluates; it never refuses.

| set | what | box |
|---|---|---|
| `metzger_1986_staggered_nu()` (default) | `0.135 Re_D^0.69 (X/D)^-0.34`, pin + endwall | H/D <= 3, 2 <= S/D <= 4, 1.5 <= X/D <= 5, Re 1e3-1e5, 10 rows |
| `metzger_1982_staggered_friction()` (default) | `0.317 Re^-0.132` / `1.76 Re^-0.318`, C1 blend at 1e4 | 0.5 <= H/D <= 6, 2 <= S/D <= 4 |
| `vanfossen_1982_staggered_nu()` | `0.153 Re_D'^0.685` on D' = 4V/A_t | H/D 0.5-2, 4 rows |
| `damerow_1972_staggered_friction()` | `2.06 (S/D)^-1.1 Re^-0.16`, per (N-1) rows | S/D 4.24-7.07 |
| `chyu_1998_nu(arrangement, surface)` | Han Table 4.7, staggered or inline, pin / endwall / total | S/D = X/D = 2.5, H/D = 1 |
| `chyu_1998_inline_over_staggered()` | inline/staggered Nu ratio, no friction ratio | Chyu's geometry |
| `chyu_1990_nu(arrangement, fillet=False)` | Chyu 1990 Table 2, pin surface, straight or fillet pins | S/D = X/D = 2.5, H/D = 1, Re 5e3-3e4 |
| `chyu_1990_friction(arrangement, fillet=False)` | Chyu 1990 Fig. 6, digitised and **Fitted**; inline straight = 0.1693/4, Re-independent | same, Re about 9e3-2.2e4 |
| `chyu_1990_fillet_over_straight(arrangement)` | fillet/straight Nu ratio (staggered 0.78-0.92), no friction ratio | same |

**Inline friction** comes from one rig at one geometry: Chyu (1990) Fig. 6,
digitised (points in `validation/cooling/data/chyu1990/`). Chyu's staggered
straight-pin friction sits 16-25% below Metzger's at the same geometry, a
disagreement between the labs that is reported, not reconciled.

**Fin efficiency** of conducting pins (each wall feeds the pin to mid-height):

```python
frac = cb.pin_fin_area_fractions(g)                 # per wall, over base area
e = cb.pin_fin_array_efficiency(h=2000.0, k_pin=20.0, D=1e-3, H=1e-3,
                                A_f_over_A_t=frac.pin_over_total)
e.eta_fin, e.eta_t, e.deta_t_dh                     # 0.968, 0.993, analytic
```

**Supply your own** by copying the default and saying so:

```python
mine = cb.metzger_1986_staggered_nu()
mine.name, mine.source = "rig_7", "measured 2026-10"
mine.provenance = cb.RibProvenance.User
mine.C = 0.150
cb.validate_pin_fin_nu_set(mine)
```

To use a set in a network, see **Pin-Fin Channels** (`PinFinModel`).

Not yet carried, and declared rather than assumed: row-count correction (a
geometry with a different `N_rows` is flagged), channel convergence, long pins
(H/D > 3), and fillets beyond Chyu's single geometry.

### Enhanced Cooling Surfaces

Removed in 0.7.0. The pin-fin, dimple, rib and impingement correlations could
not be traced to their cited sources -- the rib friction multiplier was 4-5x
below the only rib datum in the repository. Provenanced rib, jet
impingement and pin-fin correlations are re-added above (issues #334, #337,
#335); see issue #339 for the rebuild.

`channel_smooth` and the base convective correlations (Gnielinski,
Dittus-Boelter, Sieder-Tate, Petukhov) are unaffected, as are the user-set
`Nu_multiplier` and `f_multiplier` knobs on `ConvectiveSurface`.

---

## Advanced Thermodynamics

### Complete State Functions
```python
# Air properties bundle
air = cb.air_properties(T=300, P=101325, humidity=0.0)
print(f"Density: {air.density} kg/m³")
print(f"Viscosity: {air.viscosity} Pa·s")
print(f"Thermal conductivity: {air.thermal_conductivity} W/(m·K)")

# Complete thermodynamic state
thermo = cb.thermo_state(T=300, P=101325, X=air)
complete = cb.complete_state(T=300, P=101325, X=air)
```

### Derivatives
```python
# Temperature derivatives
dh_dT = cb.dh_dT(T=300, X=air)
ds_dT = cb.ds_dT(T=300, X=air)
dcp_dT = cb.dcp_dT(T=300, X=air)
dg_dT = cb.dg_over_RT_dT(T=300, X=air)
```

### Advanced Inverse Solvers
```python
# Molar-basis inverses (targets in J/mol or J/(mol*K))
T_from_h   = cb.calc_T_from_h(h_target=3e5, X=air)            # J/mol -> K
T_from_s   = cb.calc_T_from_s(s_target=7000, P=101325, X=air) # J/(mol*K), Pa -> K
T_from_cp  = cb.calc_T_from_cp(cp_target=29.1, X=air)         # J/(mol*K) -> K (29 is below air's cp)
roots_cp   = cb.calc_T_from_cp_roots(cp_target=29.1, X=air)  # every root: [CpRoot(T=...), ...]
T_from_u   = cb.calc_T_from_u(u_target=2.2e5, X=air)          # J/mol -> K

# Mass-specific inverses (targets in J/kg or J/(kg*K))
T_from_h_m = cb.calc_T_from_h_mass(h_mass_target=3e5, X=air)
T_from_s_m = cb.calc_T_from_s_mass(s_mass_target=7000, P=101325, X=air)
T_from_u_m = cb.calc_T_from_u_mass(u_mass_target=2.2e5, X=air)

# Flash calculations — solve T from two independent properties
T_from_sv  = cb.calc_T_from_sv_mass(s_mass_target=7000, v_mass_target=0.8, X=air)
T_from_sh  = cb.calc_T_from_sh_mass(s_mass_target=7000, h_mass_target=3e5, X=air)
```

### Flow Analysis
```python
# Bulk velocity and kinetic energy
v = cb.bulk_velocity(m_dot=0.5, rho=1.2, area=1e-3)   # m/s
KE = cb.kinetic_energy(v=v)                             # J/kg

# Dimensionless numbers
# From full mixture state (recomputes rho, mu internally):
Re = cb.reynolds(T=300, P=101325, X=air, V=v, L=0.05)
# From pre-computed state properties (use when CompleteState is available):
Re = cb.reynolds_from_state(rho=1.2, v=v, L=0.05, mu=1.8e-5)

# Stagnation conditions from static state and velocity
Tt, Pt = cb.stagnation_from_static(T=300, P=101325, v=v, X=air)  # (K, Pa)

# One at a time, from the Mach number (variable cp), and back
Tt = cb.Tt_from_static(T=300, M=0.5, X=air)
Pt = cb.Pt_from_static(P=101325, T=300, M=0.5, X=air)
Tt = cb.Tt_from_static_v(T=300, v=v, X=air)
T = cb.T_from_stagnation(Tt=Tt, M=0.5, X=air)
P = cb.P_from_stagnation(Pt=Pt, Tt=Tt, M=0.5, X=air)
```

### Stagnation names: Tt/Pt (#782)

`P`/`T` are static and `Pt`/`Tt` total (stagnation) everywhere in the API:
the pressure rule of #773 holds for temperature too. The spellings below were
deprecated in 0.8.0 and **removed in 0.9.0**; an old keyword now raises
`TypeError`, an old function name or field `AttributeError`.

| Removed in 0.9.0 | Use |
|---|---|
| `T0=`, `P0=` keywords (`nozzle_flow`, `nozzle_cd`, `nozzle_thrust_cd`, `solve_A_eff_from_mdot`, `solve_P_back_from_mdot`, `critical_pressure_ratio`, `mach_from_pressure_ratio`, `mass_flux_isentropic`, `duct_static_state`, `duct_htc_static_state(_dY)`, `T_from_stagnation`, `P_from_stagnation`); `compressible.nozzle_flow(T=, P=)` | `Tt=`, `Pt=` |
| `T0_from_static`, `T0_from_static_v`, `P0_from_static`, `solve_P0_from_mdot` | `Tt_from_static`, `Tt_from_static_v`, `Pt_from_static`, `solve_Pt_from_mdot` |
| `heat_entropy_flux(T0_in=, T0_out=)`; `combustor_head_loss(T_in=, P_in=)`; `chamber_merge_face_state(T_main=, T=)` | `Tt_in=`, `Tt_out=`; `Tt_in=`, `Pt_in=`; `Tt_main=`, `Tt=` |
| fields `MixtureState.T_total`/`.P_total` (and constructor keywords), `MassStream.T`/`.P_total`, `StreamJacobian.d_T`/`.d_P_total`, `NozzleSolution.T0`/`.P0`, `DuctStaticState.d*_dT0`/`d*_dP0`, `CombustorLossResult.dxi_dT_in`/`dxi_dP_in`, `MixerResult.T_mix`, `.P_total_mix`, `.dT_mix_*`, `.dP_total_mix_d_stream` | the same with `Tt`/`Pt` (`MixtureState.Tt`, `MassStream.Pt`, `dT_dTt`, `dxi_dTt_in`, `Tt_mix`, ...) |
| internal solver tools: old keywords of `momentum_chamber_residual_and_jacobian`, `orifice_tap_pressure_and_jacobian`, the ejector kernels (`p0=`, `t0=`), `sonic_static_pressure`, `approach_static_pressure`, `orifice_compressible_mdot_and_jacobian`, `T_up=` of the orifice/channel residuals; fields `SonicStatic.dP_dT0`, `OrificeResult.d_mdot_dT_up`, `ChannelSolverResult.d_dP_dT_up`; `T0_/P0_from_static_and_jacobian_M` | `Pt=`, `Tt=`, `Pt_up=`, `Tt_up=`; `dP_dTt`, `d_mdot_dTt_up`, `d_dP_dTt_up`; `Tt_/Pt_from_static_and_jacobian_M` |

Kept: `Pt0`/`Tt0` of the Fanno functions (`fanno_channel_flow`,
`fanno_rayleigh_*`, `area_change_*`) -- the `t` already says total and the
`0` is the inlet station, paired with `Tt_exit` -- and `T_from_stagnation` /
`P_from_stagnation`, whose names already say what they take. Network result
keys and `NetworkMixtureState` are the next section.

### Static T in network results (0.8)

`P` and `T` are static, `Pt` and `Tt` total, in network results too (#782):

| Key | Until 0.7 | Since 0.8 |
|---|---|---|
| `"<node>.T"` | stagnation T | **static** T; the stagnation T is `"<node>.Tt"` |
| element `T_in`, `T_out` (`result["__element_diag__"][id]`) | stagnation T of the end nodes | **static** T of the end nodes, as `P_in`/`P_out`; `Tt_in`/`Tt_out` unchanged |
| `"<chamber>.T_static"` | a momentum chamber's static T | `"<chamber>.T"` (the alias was removed in 0.9.0) |
| junction `port_<i>_T` (`__element_diag__`) | the port's stagnation T | **static** T of the port; the stagnation T is the new `port_<i>_Tt` |
| a crossflow station's `.T`, `.Pt`, `.mach` | Tt, the static P, 0 | the station's moving state: static T and Mach at the station's mean flow, its true Pt (the node's `.P` is the static P at the middle of the station) |
| an area `CombustorNode`'s `.T` | Tt | the static T at G = m / area (its `.P` = `.Pt`, a large slow volume) |
| `NetworkMixtureState.T` | a field, filled with `Tt` | a computed property, the static T |

The static T is regime-aware. A plenum (`P = Pt`) and an incompressible
node report `T = Tt`, the low-Mach limit their closures solve. A moving
compressible node's `T` lies on the isentrope through (`Tt`, `Pt`) at `P`;
a momentum chamber reports the static T its closure solved, which is that
point.

The same key name changed meaning, so code that read `"<node>.T"` as the
total temperature keeps running with a different number. Read `Tt` for the
total temperature. 0.8 warned once (`FutureWarning`) on the first such read;
the warning was removed in 0.9.0 and `solve()` returns a plain `dict`.

`NetworkMixtureState(P, Pt, *, Tt, m_dot=0.0, Y, incompressible=False)`:
since 0.9.0 `Tt`, `m_dot` and `Y` are keyword-only and `T` is no longer a
constructor argument or settable (it is the computed static T). The 0.8
signature `(P, Pt, T, Tt, m_dot, Y)` accepted `T` with a `DeprecationWarning`;
a positional call written for it now raises `TypeError` instead of shifting
its arguments.

No residual reads `NetworkMixtureState.T` (#825), so solutions are
unchanged. The residuals already form the static state wherever a node
moves (the chamber closure, junction port densities, crossflow stations;
#576, #782). The reported diagnostics did not: an element's inlet and
outlet properties (`rho`, `a`, Mach, transport) were evaluated at
(`Tt`, `P`), up to 11% off in T at a moving compressible node. They take
the static state since #824. Three residual inputs still read `Tt` at a
moving node, as closure choices: the legacy tee's `gamma_eff`, the
ejector's choke-plane gamma (its sources take inlet conditions) and a
loss correlation's `T_in`.

A residual's Jacobian names a node's stagnation-T column `"<node>.Tt"`, the
quantity it is (#782). A custom element that still keys it `"<node>.T"`, the
0.7 name (relayed with a `DeprecationWarning` in 0.8), raises `KeyError`
naming `"<node>.Tt"` since 0.9.0.

### The Stolz metering Cd (#876)

The Stolz correlation (ISO 5167:1980) was deprecated in 0.8.0 and **removed in
0.9.0**. ISO replaced it with Reader-Harris/Gallagher in ISO 5167-1:1991/Amd
1:1998 (and ISO 5167-2:2003, 2022). Reader-Harris/Gallagher has a sourced range
of use (`metering_cd_in_range`) and fits the EEC/API database better over
Stolz's own range (Reader-Harris & Sattary 1996, Table 10: s 0.245 % against
0.390 %).

| Removed in 0.9.0 | Use | Now |
|---|---|---|
| `OrificeElement(correlation="Stolz")` | `correlation="ReaderHarrisGallagher"` (the default) | `ValueError` naming the replacement |
| `cb.MeteringCdCorrelation.Stolz`, C++ `Cd_Stolz` | `cb.MeteringCdCorrelation.ReaderHarrisGallagher`, `Cd_ReaderHarrisGallagher` | `AttributeError` |
| GUI orifice "Stolz (ISO 5167:1980, deprecated)" and `SolveResponse.deprecations` | Reader-Harris/Gallagher | a saved file with `"correlation": "Stolz"` is refused with a message naming `ReaderHarrisGallagher` |

### Al 6061 conductivity (#761)

`cb.k_aluminum_6061` and the database material `"al6061"` cited no source;
deprecated in 0.8.0, removed in 0.9.0. Give the wall layer its conductivity
(`material="custom"`) instead; `get_material_conductivity("al6061", T)` raises
`ValueError`.

---

## Psychrometrics (Humid Air)

```python
air = cb.HumidAir()
air.set_TP_RH(300.0, 101325.0, 0.5)

print(f"Humidity ratio: {air.humidity_ratio} kg/kg")
print(f"Dewpoint: {air.dewpoint} K")
print(f"Underlying state P: {air.state.P} Pa")
```

`cb.dewpoint(T, P, RH)` returns the dewpoint [K] for RH in (0, 1]; at RH = 1
it is the dry-bulb `T` itself. A supersaturated RH > 1 raises `ValueError`
(the vapour-pressure step refuses it; the dewpoint solve, clamped to T, would
otherwise return `T`), as do RH <= 0 and a non-finite RH.
