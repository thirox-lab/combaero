# combaero Units Reference

This document defines the SI-based unit system used throughout the combaero library.
All functions use consistent units to avoid conversion errors.

**Auto-generated from `include/units_data.h` - do not edit manually.**

---

## Base SI Units

| Quantity           | Unit   | Symbol |
|--------------------|--------|--------|
| Temperature        | Kelvin | K      |
| Pressure           | Pascal | Pa     |
| Mass               | kilogram | kg   |
| Length             | meter  | m      |
| Time               | second | s      |
| Amount of substance| mole   | mol    |

---

## Derived Units

| Quantity             | Unit        | Symbol     |
|----------------------|-------------|------------|
| Energy               | Joule       | J = kg*m^2/s^2 |
| Power                | Watt        | W = J/s    |
| Force                | Newton      | N = kg*m/s^2 |
| Dynamic viscosity    | Pascal-second | Pa*s     |
| Kinematic viscosity  | -           | m^2/s      |
| Thermal conductivity | -           | W/(m*K)    |
| Specific heat        | -           | J/(mol*K) or J/(kg*K) |
| Enthalpy             | -           | J/mol or J/kg |
| Entropy              | -           | J/(mol*K) or J/(kg*K) |
| Density              | -           | kg/m^3     |
| Velocity             | -           | m/s        |
| Area                 | -           | m^2        |
| Volume               | -           | m^3        |
| Mass flow rate       | -           | kg/s       |
| Volumetric flow rate | -           | m^3/s      |

---

## Physical Constants

| Constant             | Value              | Unit       |
|----------------------|--------------------|------------|
| Universal gas constant R | 8.31446261815324 | J/(mol*K)  |
| Boltzmann constant   | 1.380649e-23       | J/K        |
| Avogadro's number    | 6.02214076e23      | 1/mol      |
| Standard gravity g0  | 9.80665            | m/s^2      |

---

## Function Units by Module

### thermo.h

#### Species Data

| Function                       | Input Units | Output Unit |
|--------------------------------|-------------|-------------|
| `species_name`                 | index: -    | str         |
| `species_index_from_name`      | name: str   | -           |
| `species_molar_mass`           | index: -    | g/mol       |
| `species_molar_mass_from_name` | name: str   | g/mol       |
| `num_species`                  | -           | -           |

#### Dimensionless NASA Polynomials

| Function    | Input Units | Output Unit |
|-------------|-------------|-------------|
| `cp_R`      | T: K        | - (Cp/R)    |
| `h_RT`      | T: K        | - (H/(R*T)) |
| `s_R`       | T: K        | - (S/R)     |
| `g_over_RT` | T: K        | - (G/(R*T)) |

#### Per-Species Properties

| Function     | Input Units | Output Unit |
|--------------|-------------|-------------|
| `cp_species` | T: K        | J/(mol*K)   |
| `h_species`  | T: K        | J/mol       |
| `s_species`  | T: K        | J/(mol*K)   |

#### Mixture Properties (Molar Basis)

| Function           | Input Units                        | Output Unit               |
|--------------------|------------------------------------|---------------------------|
| `cp`               | T: K, X: mol/mol                   | J/(mol*K)                 |
| `cv`               | T: K, X: mol/mol                   | J/(mol*K)                 |
| `h`                | T: K, X: mol/mol                   | J/mol                     |
| `u`                | T: K, X: mol/mol                   | J/mol                     |
| `s`                | T: K, X: mol/mol, P: Pa, P_ref: Pa | J/(mol*K)                 |
| `dh_dT`            | T: K, X: mol/mol                   | J/(mol*K)                 |
| `ds_dT`            | T: K, X: mol/mol                   | J/(mol*K^2)               |
| `dcp_dT`           | T: K, X: mol/mol                   | J/(mol*K^2)               |
| `dcp_R_dT`         | species_idx: -, T: K               | 1/K                       |
| `set_thermo_model` | name: -                            | -                         |
| `thermo_model`     | -                                  | - ('nasa9' or 'nasa9_c2') |

#### Mixture Properties (Mass/Other Basis)

| Function                           | Input Units                                                                                                                                            | Output Unit                                            |
|------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------|
| `density`                          | T: K, P: Pa, X: mol/mol                                                                                                                                | kg/m^3                                                 |
| `molar_volume`                     | T: K, P: Pa                                                                                                                                            | m^3/mol                                                |
| `specific_gas_constant`            | X: mol/mol                                                                                                                                             | J/(kg*K)                                               |
| `isentropic_expansion_coefficient` | T: K, X: mol/mol                                                                                                                                       | - (gamma)                                              |
| `speed_of_sound`                   | T: K, X: mol/mol                                                                                                                                       | m/s                                                    |
| `cp_mass`                          | T: K, X: mol/mol                                                                                                                                       | J/(kg*K)                                               |
| `cp_secant_mass_and_jacobian`      | T1: K, T2: K, X: mol/mol                                                                                                                               | cp: J/(kg*K), dcp_dT1: J/(kg*K^2), dcp_dT2: J/(kg*K^2) |
| `cp_secant_mass_dY`                | T1: K, T2: K, X: mol/mol                                                                                                                               | J/(kg*K) (per species)                                 |
| `cv_mass`                          | T: K, X: mol/mol                                                                                                                                       | J/(kg*K)                                               |
| `h_mass`                           | T: K, X: mol/mol                                                                                                                                       | J/kg                                                   |
| `s_mass`                           | T: K, X: mol/mol, P: Pa, P_ref: Pa                                                                                                                     | J/(kg*K)                                               |
| `u_mass`                           | T: K, X: mol/mol                                                                                                                                       | J/kg                                                   |
| `air_properties`                   | T: K, P: Pa, humidity: - (0-1)                                                                                                                         | AirProperties struct                                   |
| `thermo_state`                     | T: K, P: Pa, X: mol/mol, P_ref: Pa (default 101325)                                                                                                    | ThermoState struct                                     |
| `transport_state`                  | T: K, P: Pa, X: mol/mol                                                                                                                                | TransportState struct                                  |
| `complete_state`                   | T: K, P: Pa, X: mol/mol, P_ref: Pa (default 101325)                                                                                                    | CompleteState struct                                   |
| `combustion_state`                 | X_fuel: mol/mol, X_ox: mol/mol, phi: -, T_reactants: K, P: Pa, fuel_name: str (default ''), method, smooth_phi0: bool, smooth_phi1: bool, k0: -, k1: - | CombustionState struct                                 |
| `combustion_state_from_streams`    | fuel_stream: Stream, ox_stream: Stream, fuel_name: str (default ''), method, smooth_phi0: bool, smooth_phi1: bool, k0: -, k1: -                        | CombustionState struct                                 |

#### Inverse Solvers

| Function               | Input Units                                                               | Output Unit   |
|------------------------|---------------------------------------------------------------------------|---------------|
| `calc_T_from_h`        | h_target: J/mol, X: mol/mol                                               | K             |
| `calc_T_from_s`        | s_target: J/(mol*K), P: Pa, X: mol/mol                                    | K             |
| `calc_T_from_cp`       | cp_target: J/(mol*K), X: mol/mol                                          | K             |
| `calc_T_from_cp_roots` | cp_target: J/(mol*K), X: mol/mol, T_min: K, T_max: K, tol: K, max_iter: - | CpRoot        |
| `CpRoot::T`            | K                                                                         | K             |
| `CpRoot::dT_dcp`       | K/(J/(mol*K))                                                             | K/(J/(mol*K)) |
| `CpRoot::double_root`  | -                                                                         | -             |
| `calc_T_from_u`        | u_target: J/mol, X: mol/mol                                               | K             |
| `calc_T_from_h_mass`   | h_mass_target: J/kg, X: mol/mol                                           | K             |
| `calc_T_from_s_mass`   | s_mass_target: J/(kg*K), P: Pa, X: mol/mol                                | K             |
| `calc_T_from_u_mass`   | u_mass_target: J/kg, X: mol/mol                                           | K             |
| `calc_T_from_sv_mass`  | s_mass_target: J/(kg*K), v_mass_target: m^3/kg, X: mol/mol                | K             |
| `calc_T_from_sh_mass`  | s_mass_target: J/(kg*K), h_mass_target: J/kg, X: mol/mol                  | K             |

---

### transport.h

#### Transport Properties

| Function                       | Input Units                           | Output Unit                   |
|--------------------------------|---------------------------------------|-------------------------------|
| `viscosity`                    | T: K, P: Pa, X: mol/mol               | Pa*s                          |
| `viscosity_dX`                 | T: K, P: Pa, X: mol/mol               | Pa*s (per species)            |
| `viscosity_mass_gradient`      | T: K, P: Pa, X: mol/mol               | Pa*s (per unit mass fraction) |
| `transport_and_dT`             | T: K, P: Pa, X: mol/mol               | TransportDerivatives          |
| `TransportDerivatives::mu`     | -                                     | Pa*s                          |
| `TransportDerivatives::dmu_dT` | -                                     | Pa*s/K                        |
| `TransportDerivatives::k`      | -                                     | W/(m*K)                       |
| `TransportDerivatives::dk_dT`  | -                                     | W/(m*K^2)                     |
| `omega22_and_derivative`       | T_star: -, delta_star: -              | Omega22Result                 |
| `Omega22Result::omega`         | -                                     | -                             |
| `Omega22Result::domega_dTstar` | -                                     | -                             |
| `thermal_conductivity`         | T: K, P: Pa, X: mol/mol               | W/(m*K)                       |
| `prandtl`                      | T: K, P: Pa, X: mol/mol               | - (Pr)                        |
| `kinematic_viscosity`          | T: K, P: Pa, X: mol/mol               | m^2/s                         |
| `thermal_diffusivity`          | T: K, P: Pa, X: mol/mol               | m^2/s                         |
| `reynolds`                     | T: K, P: Pa, X: mol/mol, V: m/s, L: m | - (Re)                        |
| `reynolds_from_state`          | rho: kg/m^3, v: m/s, L: m, mu: Pa*s   | - (Re)                        |
| `peclet`                       | T: K, P: Pa, X: mol/mol, V: m/s, L: m | - (Pe)                        |

---

### compressible.h

#### Nozzle Flow

| Function                   | Input Units                                                  | Output Unit              |
|----------------------------|--------------------------------------------------------------|--------------------------|
| `nozzle_flow`              | Tt: K, Pt: Pa, P_back: Pa, A_eff: m^2, X: mol/mol            | CompressibleFlowSolution |
| `nozzle_quasi1d`           | Tt: K, Pt: Pa, P_exit: Pa, x: m, A: m^2, X: mol/mol          | NozzleSolution           |
| `nozzle_cd`                | Tt: K, Pt: Pa, P_exit: Pa, A: m^2, x: m, X: mol/mol          | NozzleSolution           |
| `solve_A_eff_from_mdot`    | Tt: K, Pt: Pa, P_back: Pa, mdot_target: kg/s, X: mol/mol     | m^2                      |
| `solve_P_back_from_mdot`   | Tt: K, Pt: Pa, A_eff: m^2, mdot_target: kg/s, X: mol/mol     | Pa                       |
| `solve_Pt_from_mdot`       | Tt: K, P_back: Pa, A_eff: m^2, mdot_target: kg/s, X: mol/mol | Pa                       |
| `critical_pressure_ratio`  | Tt: K, Pt: Pa, X: mol/mol                                    | - (P*/Pt)                |
| `mach_from_pressure_ratio` | Tt: K, Pt: Pa, P: Pa, X: mol/mol                             | - (M)                    |
| `mass_flux_isentropic`     | Tt: K, Pt: Pa, P: Pa, X: mol/mol                             | kg/(m^2*s)               |

#### Fanno Flow

| Function                          | Input Units                                                                                                                                                                                                | Output Unit        |
|-----------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------|
| `fanno_channel`                   | T_in: K, P_in: Pa, u_in: m/s, L: m, D: m, f: -, X: mol/mol                                                                                                                                                 | FannoSolution      |
| `fanno_channel_rough`             | T_in: K, P_in: Pa, u_in: m/s, L: m, D: m, roughness: m, X: mol/mol, correlation: str, f_multiplier: - (default 1.0)                                                                                        | FannoSolution      |
| `fanno_max_length`                | T_in: K, P_in: Pa, u_in: m/s, D: m, f: -, X: mol/mol                                                                                                                                                       | m                  |
| `fanno_state_at_mach`             | G: kg/(m^2*s), Tt: K, M: -, X: mol/mol                                                                                                                                                                     | FannoMachState     |
| `fanno_dx_dmach`                  | G: kg/(m^2*s), Tt: K, M: -, X: mol/mol, D: m, roughness: m, friction_model: str, f_multiplier: -                                                                                                           | m                  |
| `fanno_length_between`            | G: kg/(m^2*s), Tt: K, M1: -, M2: -, X: mol/mol, D: m, roughness: m, friction_model: str, f_multiplier: -                                                                                                   | m                  |
| `fanno_inlet_mach`                | Pt: Pa, Tt: K, G: kg/(m^2*s), X: mol/mol                                                                                                                                                                   | -                  |
| `fanno_sonic_mass_flux`           | Pt: Pa, Tt: K, X: mol/mol                                                                                                                                                                                  | kg/(m^2*s)         |
| `fanno_duct`                      | Pt: Pa, Tt: K, G: kg/(m^2*s), X: mol/mol, L: m, D: m, roughness: m, friction_model: str, f_multiplier: -                                                                                                   | FannoDuctResult    |
| `fanno_channel_flow`              | Pt0: Pa, Tt0: K, X: mol/mol, P_target: Pa, exit_total: bool, L: m, D: m, roughness: m, friction_model: str, f_multiplier: -, M_exit_guess: -, dY_directions: kg/kg, with_derivatives: bool                 | FannoChannelFlow   |
| `fanno_choked_mass_flux`          | Pt: Pa, Tt: K, X: mol/mol, L: m, D: m, roughness: m, friction_model: str, f_multiplier: -                                                                                                                  | kg/(m^2*s)         |
| `fanno_rayleigh_march`            | G: kg/(m^2*s), Pt0: Pa, Tt0: K, Tt_exit: K, X: mol/mol, L: m, D: m, roughness: m, friction_model: str, f_multiplier: -                                                                                     | FannoRayleighMarch |
| `fanno_rayleigh_choked_mass_flux` | Pt0: Pa, Tt0: K, Tt_exit: K, X: mol/mol, L: m, D: m, roughness: m, friction_model: str, f_multiplier: -                                                                                                    | kg/(m^2*s)         |
| `fanno_rayleigh_channel_flow`     | Pt0: Pa, Tt0: K, Tt_exit: K, X: mol/mol, P_target: Pa, exit_total: bool, L: m, D: m, roughness: m, friction_model: str, f_multiplier: -, G_guess: kg/(m^2*s), dY_directions: kg/kg, with_derivatives: bool | FannoChannelFlow   |

#### Thrust

| Function        | Input Units               | Output Unit  |
|-----------------|---------------------------|--------------|
| `nozzle_thrust` | NozzleSolution, P_amb: Pa | ThrustResult |

---

### incompressible.h

#### Thermo-Aware High-Level Functions

| Function                | Input Units                                                                 | Output Unit                |
|-------------------------|-----------------------------------------------------------------------------|----------------------------|
| `channel_flow`          | T: K, P: Pa, X: mol/mol, u: m/s, L: m, D: m, f: -                           | IncompressibleFlowSolution |
| `channel_flow_rough`    | T: K, P: Pa, X: mol/mol, u: m/s, L: m, D: m, roughness: m, correlation: str | IncompressibleFlowSolution |
| `orifice_flow_thermo`   | T: K, P: Pa, X: mol/mol, P_back: Pa, A: m^2, Cd: -                          | IncompressibleFlowSolution |
| `channel_pressure_drop` | T: K, P: Pa, X: mol/mol, v: m/s, D: m, L: m, roughness: m, correlation: str | tuple(Pa, -, -)            |

#### Bernoulli & Orifice

| Function                      | Input Units                                            | Output Unit |
|-------------------------------|--------------------------------------------------------|-------------|
| `bernoulli_P2`                | P1: Pa, v1: m/s, v2: m/s, rho: kg/m^3, dz: m, g: m/s^2 | Pa          |
| `bernoulli_v2`                | P1: Pa, P2: Pa, v1: m/s, rho: kg/m^3, dz: m, g: m/s^2  | m/s         |
| `orifice_mdot`                | P1: Pa, P2: Pa, A: m^2, Cd: -, rho: kg/m^3             | kg/s        |
| `orifice_Q`                   | P1: Pa, P2: Pa, A: m^2, Cd: -, rho: kg/m^3             | m^3/s       |
| `orifice_velocity`            | P1: Pa, P2: Pa, rho: kg/m^3                            | m/s         |
| `orifice_area`                | mdot: kg/s, P1: Pa, P2: Pa, Cd: -, rho: kg/m^3         | m^2         |
| `orifice_dP`                  | mdot: kg/s, A: m^2, Cd: -, rho: kg/m^3                 | Pa          |
| `expansibility_factor`        | beta: -, dP: Pa, P_upstream: Pa, kappa: -              | -           |
| `pressure_loss`               | v: m/s, rho: kg/m^3, zeta: -                           | Pa          |
| `velocity_from_pressure_loss` | dP: Pa, rho: kg/m^3, zeta: -                           | m/s         |
| `zeta_from_Cd`                | Cd: -                                                  | -           |
| `Cd_from_zeta`                | zeta: -                                                | -           |

---

### friction.h

#### Friction Factor Correlations

| Function             | Input Units   | Output Unit |
|----------------------|---------------|-------------|
| `friction_haaland`   | Re: -, e_D: - | - (f)       |
| `friction_serghides` | Re: -, e_D: - | - (f)       |
| `friction_colebrook` | Re: -, e_D: - | - (f)       |
| `friction_petukhov`  | Re: -         | - (f)       |

---

### orifice.h

#### Discharge Coefficients

| Function                                     | Input Units                                            | Output Unit             |
|----------------------------------------------|--------------------------------------------------------|-------------------------|
| `Cd_sharp_thin_plate`                        | OrificeGeometry, OrificeState                          | - (Cd)                  |
| `Cd_sharp_thin_plate_and_dRe`                | OrificeGeometry, OrificeState                          | - (Cd), - (dCd/dRe_D)   |
| `metering_cd_and_dRe`                        | correlation: -, geom: OrificeGeometry, Re_D: -         | - (Cd), - (dCd/dRe_D)   |
| `metering_cd_in_range`                       | correlation: -, geom: OrificeGeometry, Re_D: -         | - (bool)                |
| `metering_cd_re_min`                         | correlation: - (RHG only), geom: OrificeGeometry       | - (Re_D)                |
| `permanent_loss_ratio_and_dC`                | beta: -, C: -                                          | tuple(-, -)             |
| `velocity_of_approach_factor`                | beta: -                                                | -                       |
| `jobson_force_defect`                        | Ci: -, beta: -, Cv: -                                  | - (f)                   |
| `jobson_ci_from_cd`                          | Cd: -, beta: -, Cv: -                                  | - (Ci)                  |
| `jobson_contraction_coefficient`             | r: -, Ci: -, n: -                                      | - (C)                   |
| `jobson_expansibility`                       | r: -, beta: -, Ci: -, n: -                             | - (eps)                 |
| `JOBSON_1955_BETA_MAX`                       | -                                                      | -                       |
| `Cd`                                         | OrificeGeometry, OrificeState                          | - (Cd)                  |
| `mcgreehan_schotsch_1988_cd`                 | Re: -, r_over_d: -, L_over_d: -, U1_over_Vi: -, eps: - | - (Cd)                  |
| `mcgreehan_schotsch_1988_cd_and_derivatives` | Re: -, r_over_d: -, L_over_d: -, U1_over_Vi: -         | - (Cd, dCd/dRe, dCd/du) |
| `mcgreehan_schotsch_1988_crossflow_cd`       | cd_base: - (Cd), U1_over_Vi: -, eps: -                 | - (Cd)                  |
| `mcgreehan_schotsch_1988_expansion_orifice`  | S: -, gamma: -                                         | - (Y)                   |
| `mcgreehan_schotsch_1988_expansion_nozzle`   | S: -, gamma: -                                         | - (Y)                   |
| `mcgreehan_schotsch_1988_expansion_factor`   | cd: - (Cd), S: -, gamma: -, eps: -                     | - (Y)                   |

---

### humidair.h

#### Humid Air Properties

| Function                          | Input Units              | Output Unit  |
|-----------------------------------|--------------------------|--------------|
| `saturation_vapor_pressure`       | T: K                     | Pa           |
| `vapor_pressure`                  | T: K, RH: - (0-1)        | Pa           |
| `humidity_ratio`                  | T: K, P: Pa, RH: - (0-1) | kg/kg        |
| `water_vapor_mole_fraction`       | T: K, P: Pa, RH: - (0-1) | mol/mol      |
| `humid_air_composition`           | T: K, P: Pa, RH: - (0-1) | mol/mol      |
| `dewpoint`                        | T: K, P: Pa, RH: - (0-1) | K            |
| `relative_humidity_from_dewpoint` | T: K, Tdp: K, P: Pa      | - (0-1)      |
| `wet_bulb_temperature`            | T: K, P: Pa, RH: - (0-1) | K            |
| `humid_air_enthalpy`              | T: K, P: Pa, RH: - (0-1) | J/kg dry air |
| `humid_air_density`               | T: K, P: Pa, RH: - (0-1) | kg/m^3       |
| `HumidAir::set_TP_RH`             | T: K, P: Pa, RH: - (0-1) | HumidAir&    |
| `HumidAir::set_TP_dewpoint`       | T: K, P: Pa, Tdp: K      | HumidAir&    |
| `HumidAir::set_TP_omega`          | T: K, P: Pa, w: kg/kg    | HumidAir&    |
| `HumidAir::RH`                    | -                        | - (0-1)      |
| `HumidAir::dewpoint`              | -                        | K            |
| `HumidAir::wet_bulb`              | -                        | K            |
| `HumidAir::humidity_ratio`        | -                        | kg/kg        |
| `HumidAir::h_mass`                | -                        | J/kg dry air |
| `HumidAir::TP_RH`                 | -                        | (K, Pa, -)   |
| `HumidAir::TP_dewpoint`           | -                        | (K, Pa, K)   |

---

### combustion.h

#### Stoichiometry

| Function                                    | Input Units                                                    | Output Unit                                                                                                    |
|---------------------------------------------|----------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------|
| `oxygen_required_per_mol_fuel`              | -                                                              | mol O2/mol fuel                                                                                                |
| `oxygen_required_per_kg_fuel`               | -                                                              | kg O2/kg fuel                                                                                                  |
| `oxygen_required_per_mol_mixture`           | X: mol/mol                                                     | mol O2/mol mix                                                                                                 |
| `complete_combustion_products_and_jacobian` | X: mol/mol, smooth_phi0: bool, smooth_phi1: bool, k0: -, k1: - | CompleteCombustionJacobian (X_b: mol/mol, dXb_dX: -)                                                           |
| `combustion_equilibrium_and_jacobian`       | T: K, X: mol/mol, P: Pa, smooth: bool, k0: -                   | EquilibriumJacobian (state, dT_dT_in: -, dT_dP: K/Pa, dT_dX: K, dXout_dX: -, dXout_dT_in: 1/K, dXout_dP: 1/Pa) |
| `oxygen_required_per_kg_mixture`            | X: mol/mol                                                     | kg O2/kg mix                                                                                                   |
| `fuel_lhv_molar`                            | X_fuel: mol/mol, reference_temperature: K                      | J/mol fuel                                                                                                     |
| `fuel_lhv_mass`                             | X_fuel: mol/mol, reference_temperature: K                      | J/kg fuel                                                                                                      |
| `dryair_required_per_mol_fuel`              | -                                                              | mol air/mol fuel                                                                                               |
| `dryair_required_per_kg_fuel`               | -                                                              | kg air/kg fuel                                                                                                 |
| `dryair_required_per_mol_mixture`           | X: mol/mol                                                     | mol air/mol mix                                                                                                |
| `dryair_required_per_kg_mixture`            | X: mol/mol                                                     | kg air/kg mix                                                                                                  |

#### Equivalence Ratio

| Function                     | Input Units        | Output Unit |
|------------------------------|--------------------|-------------|
| `equivalence_ratio_mole`     | X: mol/mol         | - (phi)     |
| `set_equivalence_ratio_mole` | phi: -, X: mol/mol | mol/mol     |
| `equivalence_ratio_mass`     | Y: kg/kg           | - (phi)     |
| `set_equivalence_ratio_mass` | phi: -, Y: kg/kg   | kg/kg       |

#### Mixture Fraction (Bilger)

| Function                               | Input Units      | Output Unit |
|----------------------------------------|------------------|-------------|
| `bilger_beta`                          | Y: kg/kg         | -           |
| `bilger_mixture_fraction`              | Y: kg/kg         | - (Z)       |
| `bilger_stoich_mixture_fraction_mass`  | Y: kg/kg         | - (Z_st)    |
| `equivalence_ratio_from_bilger_Z_mass` | Z: -, Y: kg/kg   | - (phi)     |
| `bilger_Z_from_equivalence_ratio_mass` | phi: -, Y: kg/kg | - (Z)       |

#### Complete Combustion

| Function                         | Input Units                                                                                  | Output Unit |
|----------------------------------|----------------------------------------------------------------------------------------------|-------------|
| `complete_combustion_to_CO2_H2O` | X: mol/mol, smooth: bool (smooth_phi0, default False)                                        | mol/mol     |
| `complete_combustion`            | T: K, X: mol/mol, P: Pa (default 101325), smooth_phi0: bool, smooth_phi1: bool, k0: -, k1: - | State       |
| `complete_combustion_isothermal` | T: K, X: mol/mol, P: Pa (default 101325), smooth_phi0: bool, smooth_phi1: bool, k0: -, k1: - | State       |

#### Stream Solvers

| Function                          | Input Units               | Output Unit |
|-----------------------------------|---------------------------|-------------|
| `set_fuel_stream_for_phi`         | phi: -                    | Stream      |
| `set_fuel_stream_for_Tad`         | T_ad_target: K            | Stream      |
| `set_fuel_stream_for_O2`          | X_O2_target: mol/mol      | Stream      |
| `set_fuel_stream_for_O2_dry`      | X_O2_dry_target: mol/mol  | Stream      |
| `set_fuel_stream_for_CO2`         | X_CO2_target: mol/mol     | Stream      |
| `set_fuel_stream_for_CO2_dry`     | X_CO2_dry_target: mol/mol | Stream      |
| `set_oxidizer_stream_for_Tad`     | T_ad_target: K            | Stream      |
| `set_oxidizer_stream_for_O2`      | X_O2_target: mol/mol      | Stream      |
| `set_oxidizer_stream_for_O2_dry`  | X_O2_dry_target: mol/mol  | Stream      |
| `set_oxidizer_stream_for_CO2`     | X_CO2_target: mol/mol     | Stream      |
| `set_oxidizer_stream_for_CO2_dry` | X_CO2_dry_target: mol/mol | Stream      |

---

### equilibrium.h

#### Chemical Equilibrium

| Function                          | Input Units                         | Output Unit |
|-----------------------------------|-------------------------------------|-------------|
| `wgs_equilibrium`                 | State (T: K, P: Pa, X: mol/mol)     | State       |
| `wgs_equilibrium_adiabatic`       | State                               | State       |
| `reforming_equilibrium`           | State                               | State       |
| `reforming_equilibrium_adiabatic` | State                               | State       |
| `combustion_equilibrium`          | State, smooth: bool (default False) | State       |

---

### state.h

#### State Properties

| Function       | Input Units | Output Unit |
|----------------|-------------|-------------|
| `State::T`     | -           | K           |
| `State::P`     | -           | Pa          |
| `State::X`     | -           | mol/mol     |
| `State::Y`     | -           | kg/kg       |
| `State::set_X` | X: mol/mol  | -           |
| `State::set_Y` | Y: kg/kg    | -           |
| `State::mw`    | -           | g/mol       |
| `State::cp`    | -           | J/(kg*K)    |
| `State::cv`    | -           | J/(kg*K)    |
| `State::h`     | -           | J/kg        |
| `State::u`     | -           | J/kg        |
| `State::s`     | -           | J/(kg*K)    |
| `State::rho`   | -           | kg/m^3      |
| `State::R`     | -           | J/(kg*K)    |
| `State::gamma` | -           | - (gamma)   |
| `State::a`     | -           | m/s         |
| `State::mu`    | -           | Pa*s        |
| `State::k`     | -           | W/(m*K)     |
| `State::nu`    | -           | m^2/s       |
| `State::Pr`    | -           | - (Pr)      |
| `State::alpha` | -           | m^2/s       |

#### Stream Properties

| Function           | Input Units | Output Unit |
|--------------------|-------------|-------------|
| `Stream::mdot`     | -           | kg/s        |
| `Stream::T`        | -           | K           |
| `Stream::P`        | -           | Pa          |
| `Stream::X`        | -           | mol/mol     |
| `Stream::set_mdot` | mdot: kg/s  | -           |
| `Stream::set_T`    | T: K        | -           |
| `Stream::set_P`    | P: Pa       | -           |
| `Stream::set_X`    | X: mol/mol  | -           |

---

### stagnation.h

#### Stagnation/Static Conversions

| Function                 | Input Units                                                                                   | Output Unit                                         |
|--------------------------|-----------------------------------------------------------------------------------------------|-----------------------------------------------------|
| `kinetic_energy`         | v: m/s                                                                                        | J/kg                                                |
| `bulk_velocity`          | m_dot: kg/s, rho: kg/m^3, area: m^2                                                           | m/s                                                 |
| `h0_from_static`         | h_static: J/kg, v: m/s                                                                        | J/kg                                                |
| `v_from_h0`              | h0: J/kg, h_static: J/kg                                                                      | m/s                                                 |
| `mach_number`            | v: m/s, T: K, X: mol/mol                                                                      | - (M)                                               |
| `Tt_from_static`         | T: K, M: -, X: mol/mol                                                                        | K                                                   |
| `Tt_from_static_v`       | T: K, v: m/s, X: mol/mol                                                                      | K                                                   |
| `T_from_stagnation`      | Tt: K, M: -, X: mol/mol                                                                       | K                                                   |
| `Pt_from_static`         | P: Pa, T: K, M: -, X: mol/mol                                                                 | Pa                                                  |
| `P_from_stagnation`      | Pt: Pa, Tt: K, M: -, X: mol/mol                                                               | Pa                                                  |
| `recovery_factor`        | Pr: -                                                                                         | - (r)                                               |
| `T_adiabatic_wall`       | T_static: K, v: m/s, T: K, P: Pa, X: mol/mol                                                  | K                                                   |
| `T_adiabatic_wall_mach`  | T_static: K, M: -, T: K, P: Pa, X: mol/mol                                                    | K                                                   |
| `stagnation_from_static` | T: K, P: Pa, v: m/s, X: mol/mol                                                               | tuple(K, Pa)                                        |
| `entropy_generation`     | inflows: MassStream (m_dot: kg/s, Tt: K, Pt: Pa, Y: kg/kg), outflows: MassStream, S_heat: W/K | EntropyGenerationResult (sigma: W/K, sigma_norm: -) |
| `heat_entropy_flux`      | m_dot: kg/s, Tt_in: K, Tt_out: K, Y: kg/kg, M_in: -, M_out: -                                 | W/K                                                 |

### composition.h - Composition Utilities

| Function                   | Input Units             | Output Unit |
|----------------------------|-------------------------|-------------|
| `mwmix`                    | X: mol/mol              | g/mol       |
| `molar_mass`               | X: mol/mol              | g/mol       |
| `mole_to_mass`             | X: mol/mol              | kg/kg       |
| `mass_to_mole`             | Y: kg/kg                | mol/mol     |
| `normalize_fractions`      | X: mol/mol              | mol/mol     |
| `convert_to_dry_fractions` | X: mol/mol              | mol/mol     |
| `equivalence_ratio`        | mole_fractions: mol/mol | [-]         |

### materials.h - Material Thermal Conductivity

| Function                              | Input Units                                                | Output Unit                                                                                      |
|---------------------------------------|------------------------------------------------------------|--------------------------------------------------------------------------------------------------|
| `k_inconel718`                        | T: K, extrapolate: bool (default False)                    | W/(m*K)                                                                                          |
| `k_haynes230`                         | T: K, extrapolate: bool (default False)                    | W/(m*K)                                                                                          |
| `k_hastelloy_x`                       | T: K, extrapolate: bool (default False)                    | W/(m*K)                                                                                          |
| `k_haynes188`                         | T: K, extrapolate: bool (default False)                    | W/(m*K)                                                                                          |
| `k_haynes263`                         | T: K, extrapolate: bool (default False)                    | W/(m*K)                                                                                          |
| `k_haynes25`                          | T: K, extrapolate: bool (default False)                    | W/(m*K)                                                                                          |
| `k_inconel625`                        | T: K, extrapolate: bool (default False)                    | W/(m*K)                                                                                          |
| `k_nimonic75`                         | T: K, extrapolate: bool (default False)                    | W/(m*K)                                                                                          |
| `k_copper`                            | T: K, extrapolate: bool (default False)                    | W/(m*K)                                                                                          |
| `k_alumina`                           | T: K, extrapolate: bool (default False)                    | W/(m*K)                                                                                          |
| `k_x22crmov12_1`                      | T: K, extrapolate: bool (default False)                    | W/(m*K)                                                                                          |
| `k_bondcoat_nicocraly`                | T: K, extrapolate: bool (default False)                    | W/(m*K)                                                                                          |
| `k_stainless_steel_304`               | T: K, extrapolate: bool (default False)                    | W/(m*K)                                                                                          |
| `k_stainless_steel_310`               | T: K, extrapolate: bool (default False)                    | W/(m*K)                                                                                          |
| `k_stainless_steel_316`               | T: K, extrapolate: bool (default False)                    | W/(m*K)                                                                                          |
| `k_tbc_ysz`                           | T: K, hours: h (default 0), is_ebpvd: bool (default False) | W/(m*K)                                                                                          |
| `list_materials`                      | -                                                          | list[str]                                                                                        |
| `get_material_conductivity`           | name: str, T: K, extrapolate: bool (default False)         | W/(m*K)                                                                                          |
| `get_material_conductivity_and_slope` | name: str, T: K, extrapolate: bool (default False)         | tuple: k [W/(m*K)], dk_dT [W/(m*K^2)]                                                            |
| `get_material_info`                   | name: str                                                  | dict: T_min, T_max, T_solidus, T_envelope_min, T_envelope_max [K]; source, notes, solidus_source |
| `material_in_range`                   | name: str, T: K                                            | bool                                                                                             |
| `GENERIC_METAL_SOLIDUS`               | -                                                          | K                                                                                                |

### cooling_correlations.h - Advanced Cooling Correlations

| Function                     | Input Units             | Output Unit |
|------------------------------|-------------------------|-------------|
| `thermal_performance_factor` | Nu_ratio: -, f_ratio: - | -           |

### rib_correlation.h - Parametrised rib correlations

| Function                                           | Input Units                                                                                 | Output Unit       |
|----------------------------------------------------|---------------------------------------------------------------------------------------------|-------------------|
| `han_1988_orthogonal`                              | -                                                                                           | RibCorrelationSet |
| `rallabandi_2009_high_re`                          | -                                                                                           | RibCorrelationSet |
| `han_park_1988_angled`                             | -                                                                                           | RibCorrelationSet |
| `han_1989_narrow_channel`                          | -                                                                                           | RibCorrelationSet |
| `han_zhang_lee_1991`                               | shape: RibShape, alpha_deg: deg                                                             | RibCorrelationSet |
| `validate_rib_set`                                 | correlation_set: RibCorrelationSet                                                          | None              |
| `evaluate_rib`                                     | correlation_set: RibCorrelationSet, geometry: RibGeometry, Re: -                            | RibResult         |
| `RibGeometry::e_D`                                 | -                                                                                           | -                 |
| `RibGeometry::p_e`                                 | -                                                                                           | -                 |
| `RibGeometry::W_H`                                 | -                                                                                           | -                 |
| `RibGeometry::alpha_deg`                           | -                                                                                           | deg               |
| `RibTerm::exponent`                                | -                                                                                           | -                 |
| `RibTerm::reference`                               | -                                                                                           | -                 |
| `RibResult::R`                                     | -                                                                                           | -                 |
| `RibResult::f`                                     | -                                                                                           | -                 |
| `RibResult::e_plus`                                | -                                                                                           | -                 |
| `RibResult::G`                                     | -                                                                                           | -                 |
| `RibResult::St_r`                                  | -                                                                                           | -                 |
| `RibResult::dSt_dRe`                               | -                                                                                           | -                 |
| `RibResult::G_bar`                                 | -                                                                                           | -                 |
| `RibResult::has_G_bar`                             | -                                                                                           | -                 |
| `RibCorrelationSet::name`                          | -                                                                                           | -                 |
| `RibCorrelationSet::source`                        | -                                                                                           | -                 |
| `RibCorrelationSet::validity_source`               | -                                                                                           | -                 |
| `RibCorrelationSet::provenance`                    | -                                                                                           | -                 |
| `RibCorrelationSet::symmetric`                     | -                                                                                           | -                 |
| `RibCorrelationSet::shape`                         | -                                                                                           | -                 |
| `RibCorrelationSet::R_alpha_shape`                 | -                                                                                           | -                 |
| `RibCorrelationSet::C_R`                           | -                                                                                           | -                 |
| `RibCorrelationSet::R_eD`                          | -                                                                                           | -                 |
| `RibCorrelationSet::R_pe`                          | -                                                                                           | -                 |
| `RibCorrelationSet::R_WH`                          | -                                                                                           | -                 |
| `RibCorrelationSet::R_alpha`                       | -                                                                                           | -                 |
| `RibCorrelationSet::R_quad_c0`                     | -                                                                                           | -                 |
| `RibCorrelationSet::R_quad_c1`                     | -                                                                                           | -                 |
| `RibCorrelationSet::R_quad_c2`                     | -                                                                                           | -                 |
| `RibCorrelationSet::R_quad_WH_exponent_at_90`      | -                                                                                           | -                 |
| `RibCorrelationSet::R_quad_WH_exponent_off_90`     | -                                                                                           | -                 |
| `RibCorrelationSet::R_quad_WH_cap`                 | -                                                                                           | -                 |
| `RibCorrelationSet::R_quad_narrow_c0`              | -                                                                                           | -                 |
| `RibCorrelationSet::R_quad_narrow_c1`              | -                                                                                           | -                 |
| `RibCorrelationSet::R_quad_narrow_c2`              | -                                                                                           | -                 |
| `RibCorrelationSet::R_WH_band_boundary`            | -                                                                                           | -                 |
| `RibCorrelationSet::G_shape_model`                 | -                                                                                           | -                 |
| `RibCorrelationSet::C_G`                           | -                                                                                           | -                 |
| `RibCorrelationSet::G_eD`                          | -                                                                                           | -                 |
| `RibCorrelationSet::G_pe`                          | -                                                                                           | -                 |
| `RibCorrelationSet::G_WH`                          | -                                                                                           | -                 |
| `RibCorrelationSet::G_alpha`                       | -                                                                                           | -                 |
| `RibCorrelationSet::G_eplus_exponent`              | -                                                                                           | -                 |
| `RibCorrelationSet::G_shape_alpha_exponent_square` | -                                                                                           | -                 |
| `RibCorrelationSet::G_shape_alpha_exponent_rect`   | -                                                                                           | -                 |
| `RibCorrelationSet::G_shape_pe_exponent_square`    | -                                                                                           | -                 |
| `RibCorrelationSet::G_shape_pe_exponent_rect`      | -                                                                                           | -                 |
| `RibCorrelationSet::G_narrow_C_alpha90`            | -                                                                                           | -                 |
| `RibCorrelationSet::G_narrow_C_off_axis`           | -                                                                                           | -                 |
| `RibCorrelationSet::G_narrow_WH_band_boundary`     | -                                                                                           | -                 |
| `RibCorrelationSet::G_narrow_WH_C_exponent`        | -                                                                                           | -                 |
| `RibCorrelationSet::G_narrow_WH_n_exponent`        | -                                                                                           | -                 |
| `RibCorrelationSet::C_Gbar`                        | -                                                                                           | -                 |
| `RibCorrelationSet::Gbar_eplus_exponent`           | -                                                                                           | -                 |
| `RibCorrelationSet::valid_Re`                      | -                                                                                           | -                 |
| `RibCorrelationSet::valid_eD`                      | -                                                                                           | -                 |
| `RibCorrelationSet::valid_pe`                      | -                                                                                           | -                 |
| `RibCorrelationSet::valid_WH`                      | -                                                                                           | -                 |
| `RibCorrelationSet::valid_alpha`                   | -                                                                                           | deg               |
| `RibCorrelationSet::valid_eplus`                   | -                                                                                           | -                 |
| `RibCorrelationSet::valid_Pr`                      | -                                                                                           | -                 |
| `RibCorrelationSet::accuracy_R`                    | -                                                                                           | -                 |
| `RibCorrelationSet::accuracy_G`                    | -                                                                                           | -                 |
| `StatedAccuracy::value`                            | -                                                                                           | -                 |
| `StatedAccuracy::provenance`                       | -                                                                                           | -                 |
| `StatedAccuracy::usable_as_band`                   | -                                                                                           | -                 |
| `StatedAccuracy::stated`                           | -                                                                                           | -                 |
| `StatedAccuracy::measured`                         | -                                                                                           | -                 |
| `StatedAccuracy::unstated`                         | -                                                                                           | -                 |
| `AccuracyProvenance::Unstated`                     | -                                                                                           | -                 |
| `AccuracyProvenance::Stated`                       | -                                                                                           | -                 |
| `AccuracyProvenance::Measured`                     | -                                                                                           | -                 |
| `AccuracyProvenance::name`                         | -                                                                                           | -                 |
| `AccuracyProvenance::value`                        | -                                                                                           | -                 |
| `RibRange::lo`                                     | -                                                                                           | -                 |
| `RibRange::hi`                                     | -                                                                                           | -                 |
| `RibResult::extrapolated`                          | -                                                                                           | -                 |
| `RibProvenance::Extracted`                         | -                                                                                           | -                 |
| `RibProvenance::Fitted`                            | -                                                                                           | -                 |
| `RibProvenance::User`                              | -                                                                                           | -                 |
| `RibShape::Unspecified`                            | -                                                                                           | -                 |
| `RibShape::Transverse`                             | -                                                                                           | -                 |
| `RibShape::Parallel`                               | -                                                                                           | -                 |
| `RibShape::Crossed`                                | -                                                                                           | -                 |
| `RibShape::V`                                      | -                                                                                           | -                 |
| `RibShape::Lambda`                                 | -                                                                                           | -                 |
| `RibShape::name`                                   | -                                                                                           | -                 |
| `RibShape::value`                                  | -                                                                                           | -                 |
| `evaluate_rib_ratio`                               | correlation_set: RibRatioSet, geometry: RibGeometry, Re: -, Pr: -, options: RibRatioOptions | RibRatioResult    |
| `validate_rib_ratio_set`                           | correlation_set: RibRatioSet                                                                | None              |
| `validate_rib_ratio_options`                       | options: RibRatioOptions                                                                    | None              |
| `taslim_spring_1987`                               | aspect_ratio_taslim: -, e_D: -                                                              | RibRatioSet       |
| `RatioBaseline::coeff`                             | -                                                                                           | -                 |
| `RatioBaseline::re_exponent`                       | -                                                                                           | -                 |
| `RatioBaseline::pr_exponent`                       | -                                                                                           | -                 |
| `RatioBelowFloor::Gnielinski`                      | -                                                                                           | -                 |
| `RatioBelowFloor::SourceBaseline`                  | -                                                                                           | -                 |
| `RatioBelowFloor::User`                            | -                                                                                           | -                 |
| `RatioBelowFloor::name`                            | -                                                                                           | -                 |
| `RatioBelowFloor::value`                           | -                                                                                           | -                 |
| `RibRatioOptions::below_floor`                     | -                                                                                           | -                 |
| `RibRatioOptions::user_Nu0`                        | -                                                                                           | -                 |
| `RibRatioResult::Nu`                               | -                                                                                           | -                 |
| `RibRatioResult::dNu_dRe`                          | -                                                                                           | -                 |
| `RibRatioResult::f`                                | -                                                                                           | -                 |
| `RibRatioResult::df_dRe`                           | -                                                                                           | -                 |
| `RibRatioResult::ratio_Nu`                         | -                                                                                           | -                 |
| `RibRatioResult::ratio_f`                          | -                                                                                           | -                 |
| `RibRatioResult::below_floor`                      | -                                                                                           | -                 |
| `RibRatioResult::extrapolated`                     | -                                                                                           | -                 |
| `RibRatioSet::name`                                | -                                                                                           | -                 |
| `RibRatioSet::source`                              | -                                                                                           | -                 |
| `RibRatioSet::validity_source`                     | -                                                                                           | -                 |
| `RibRatioSet::provenance`                          | -                                                                                           | -                 |
| `RibRatioSet::shape`                               | -                                                                                           | -                 |
| `RibRatioSet::symmetric`                           | -                                                                                           | -                 |
| `RibRatioSet::C_Nu`                                | -                                                                                           | -                 |
| `RibRatioSet::Nu_Re`                               | -                                                                                           | -                 |
| `RibRatioSet::Nu_eD`                               | -                                                                                           | -                 |
| `RibRatioSet::Nu_pe`                               | -                                                                                           | -                 |
| `RibRatioSet::Nu_WH`                               | -                                                                                           | -                 |
| `RibRatioSet::Nu_alpha`                            | -                                                                                           | -                 |
| `RibRatioSet::C_f`                                 | -                                                                                           | -                 |
| `RibRatioSet::f_Re`                                | -                                                                                           | -                 |
| `RibRatioSet::f_eD`                                | -                                                                                           | -                 |
| `RibRatioSet::f_pe`                                | -                                                                                           | -                 |
| `RibRatioSet::f_WH`                                | -                                                                                           | -                 |
| `RibRatioSet::f_alpha`                             | -                                                                                           | -                 |
| `RibRatioSet::Nu0_source`                          | -                                                                                           | -                 |
| `RibRatioSet::f0_source`                           | -                                                                                           | -                 |
| `RibRatioSet::Re_floor`                            | -                                                                                           | -                 |
| `RibRatioSet::Re_floor_f`                          | -                                                                                           | -                 |
| `RibRatioSet::valid_Re`                            | -                                                                                           | -                 |
| `RibRatioSet::valid_eD`                            | -                                                                                           | -                 |
| `RibRatioSet::valid_pe`                            | -                                                                                           | -                 |
| `RibRatioSet::valid_WH`                            | -                                                                                           | -                 |
| `RibRatioSet::valid_alpha`                         | -                                                                                           | -                 |
| `RibRatioSet::valid_Pr`                            | -                                                                                           | -                 |
| `RibRatioSet::accuracy_Nu`                         | -                                                                                           | -                 |
| `RibRatioSet::accuracy_f`                          | -                                                                                           | -                 |
| `RibProvenance::name`                              | -                                                                                           | -                 |
| `RibProvenance::value`                             | -                                                                                           | -                 |
| `RAlphaShape::PowerLaw`                            | -                                                                                           | -                 |
| `RAlphaShape::QuadraticAlpha`                      | -                                                                                           | -                 |
| `RAlphaShape::QuadraticAlphaTwoBand`               | -                                                                                           | -                 |
| `RAlphaShape::name`                                | -                                                                                           | -                 |
| `RAlphaShape::value`                               | -                                                                                           | -                 |
| `GShapeModel::Fixed`                               | -                                                                                           | -                 |
| `GShapeModel::SquareVsRectangular`                 | -                                                                                           | -                 |
| `GShapeModel::NarrowChannelAlphaSwitch`            | -                                                                                           | -                 |
| `GShapeModel::name`                                | -                                                                                           | -                 |
| `GShapeModel::value`                               | -                                                                                           | -                 |

### impingement_correlation.h - Jet impingement correlations

| Function                                        | Input Units                                                                                                                    | Output Unit                     |
|-------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------|---------------------------------|
| `goldstein_1986_single_jet`                     | -                                                                                                                              | SingleJetImpingementSet         |
| `validate_single_jet_set`                       | correlation_set: SingleJetImpingementSet                                                                                       | None                            |
| `single_jet_impingement_nu`                     | correlation_set: SingleJetImpingementSet, bc: ImpingementThermalBC, Re: -, L_D: -, R_D: -                                      | -                               |
| `single_jet_impingement`                        | correlation_set: SingleJetImpingementSet, bc: ImpingementThermalBC, Re: -, L_D: -, R_D: -                                      | SingleJetImpingementResult      |
| `florschuetz_1981_inline`                       | -                                                                                                                              | JetArrayCorrelationSet          |
| `florschuetz_1981_staggered`                    | -                                                                                                                              | JetArrayCorrelationSet          |
| `validate_jet_array_set`                        | correlation_set: JetArrayCorrelationSet                                                                                        | None                            |
| `jet_array_impingement_nu`                      | correlation_set: JetArrayCorrelationSet, Re_j: -, Gc_Gj: -, Pr: -, xn_d: -, yn_d: -, z_d: -                                    | JetArrayImpingementResult       |
| `crossflow_to_jet_ratio_at_x`                   | yn_d: -, z_d: -, C_D: -, x_over_xn: -                                                                                          | -                               |
| `crossflow_to_jet_ratio_at_row`                 | yn_d: -, z_d: -, C_D: -, row: -                                                                                                | -                               |
| `jet_row_gc_gj_per_flow_ratio`                  | yn_d: -, z_d: -                                                                                                                | -                               |
| `jet_row_heat_transfer`                         | correlation_set: JetArrayCorrelationSet, geometry: JetRowGeometry, m_jet: kg/s, m_crossflow: kg/s, mu: Pa*s, k: W/(m*K), Pr: - | JetRowHeatTransfer              |
| `JetRowGeometry::d`                             | m                                                                                                                              | m                               |
| `JetRowGeometry::xn_d`                          | -                                                                                                                              | -                               |
| `JetRowGeometry::yn_d`                          | -                                                                                                                              | -                               |
| `JetRowGeometry::z_d`                           | -                                                                                                                              | -                               |
| `JetRowGeometry::n_holes`                       | -                                                                                                                              | -                               |
| `JetRowHeatTransfer::h`                         | W/(m^2*K)                                                                                                                      | W/(m^2*K)                       |
| `JetRowHeatTransfer::Nu`                        | -                                                                                                                              | -                               |
| `JetRowHeatTransfer::Re_j`                      | -                                                                                                                              | -                               |
| `JetRowHeatTransfer::Gc_Gj`                     | -                                                                                                                              | -                               |
| `JetRowHeatTransfer::dh_dm_jet`                 | W/(m^2*K)/(kg/s)                                                                                                               | W/(m^2*K)/(kg/s)                |
| `JetRowHeatTransfer::dh_dm_crossflow`           | W/(m^2*K)/(kg/s)                                                                                                               | W/(m^2*K)/(kg/s)                |
| `JetRowHeatTransfer::extrapolated`              | -                                                                                                                              | -                               |
| `JetRowHeatTransfer::eta_r`                     | -                                                                                                                              | -                               |
| `JetRowHeatTransfer::deta_r_dm_jet`             | 1/(kg/s)                                                                                                                       | 1/(kg/s)                        |
| `JetRowHeatTransfer::deta_r_dm_crossflow`       | 1/(kg/s)                                                                                                                       | 1/(kg/s)                        |
| `JetRowHeatTransfer::eta_r_extrapolated`        | -                                                                                                                              | -                               |
| `combustor_fraction_loss`                       | theta: -, mdot: kg/s, k: -, xi0: -, m_ref: kg/s                                                                                | CombustorLossResult             |
| `combustor_head_loss`                           | theta: -, mdot: kg/s, Tt_in: K, Pt_in: Pa, Y_in: kg/kg, k: -, zeta0: -, area: m^2                                              | CombustorLossResult             |
| `CombustorLossResult::xi`                       | -                                                                                                                              | -                               |
| `CombustorLossResult::dxi_dtheta`               | -                                                                                                                              | -                               |
| `CombustorLossResult::dxi_dmdot`                | 1/(kg/s)                                                                                                                       | 1/(kg/s)                        |
| `CombustorLossResult::dxi_dTt_in`               | 1/K                                                                                                                            | 1/K                             |
| `CombustorLossResult::dxi_dPt_in`               | 1/Pa                                                                                                                           | 1/Pa                            |
| `CombustorLossResult::dxi_dY`                   | -                                                                                                                              | -                               |
| `impingement_eta_r`                             | Gc_Gj: -, xn_d: -, yn_d: -, z_d: -                                                                                             | ImpingementEtaR                 |
| `ImpingementEtaR::eta`                          | -                                                                                                                              | -                               |
| `ImpingementEtaR::deta_dGc_Gj`                  | -                                                                                                                              | -                               |
| `ImpingementEtaR::extrapolated`                 | -                                                                                                                              | -                               |
| `station_half_drop`                             | m_a: kg/s, m_b: kg/s, P: Pa, T: K, X: -, area: m^2, kappa: -                                                                   | StationHalfDrop                 |
| `channel_entry_drop`                            | m: kg/s, P: Pa, T: K, X: -, area: m^2, K_in: -                                                                                 | StationHalfDrop                 |
| `segment_area_change_drop`                      | m: kg/s, P: Pa, T: K, X: -, area_from: m^2, area_to: m^2                                                                       | StationHalfDrop                 |
| `StationHalfDrop::dP`                           | Pa                                                                                                                             | Pa                              |
| `StationHalfDrop::d_dm_a`                       | Pa/(kg/s)                                                                                                                      | Pa/(kg/s)                       |
| `StationHalfDrop::d_dm_b`                       | Pa/(kg/s)                                                                                                                      | Pa/(kg/s)                       |
| `StationHalfDrop::d_dP`                         | -                                                                                                                              | -                               |
| `StationHalfDrop::d_dT`                         | Pa/K                                                                                                                           | Pa/K                            |
| `CROSSFLOW_VI_FLOOR`                            | m/s                                                                                                                            | m/s                             |
| `crossflow_velocity_ratio`                      | m_a: kg/s, m_b: kg/s, P: Pa, T: K, X: -, P_down: Pa, area: m^2                                                                 | CrossflowRatio                  |
| `CrossflowRatio::U1_over_Vi`                    | -                                                                                                                              | -                               |
| `CrossflowRatio::U1`                            | m/s                                                                                                                            | m/s                             |
| `CrossflowRatio::Vi`                            | m/s                                                                                                                            | m/s                             |
| `CrossflowRatio::d_dm_a`                        | 1/(kg/s)                                                                                                                       | 1/(kg/s)                        |
| `CrossflowRatio::d_dm_b`                        | 1/(kg/s)                                                                                                                       | 1/(kg/s)                        |
| `CrossflowRatio::d_dP`                          | 1/Pa                                                                                                                           | 1/Pa                            |
| `CrossflowRatio::d_dT`                          | 1/K                                                                                                                            | 1/K                             |
| `CrossflowRatio::d_dP_down`                     | 1/Pa                                                                                                                           | 1/Pa                            |
| `rogers_hersh_1976_cd`                          | Cd0: -, G_ideal: kg/(m^2*s), G_c: kg/(m^2*s), suction: -                                                                       | CrossflowCd                     |
| `florschuetz_metzger_su_1982_cd`                | xn_d: -, yn_d: -, z_d: -, G_ideal: kg/(m^2*s), G_c: kg/(m^2*s)                                                                 | CrossflowCd                     |
| `florschuetz_metzger_su_1982_has_plate`         | xn_d: -, yn_d: -, z_d: -                                                                                                       | -                               |
| `CrossflowCd::Cd`                               | -                                                                                                                              | -                               |
| `CrossflowCd::dCd_dG_ideal`                     | 1/(kg/(m^2*s))                                                                                                                 | 1/(kg/(m^2*s))                  |
| `CrossflowCd::dCd_dG_c`                         | 1/(kg/(m^2*s))                                                                                                                 | 1/(kg/(m^2*s))                  |
| `CrossflowCd::capped`                           | -                                                                                                                              | -                               |
| `rowbury_2001_coefficients`                     | set: -, theta_deg: deg, alpha_deg: deg                                                                                         | Rowbury2001Coefficients         |
| `Rowbury2001Coefficients::A`                    | -                                                                                                                              | -                               |
| `Rowbury2001Coefficients::B`                    | -                                                                                                                              | -                               |
| `Rowbury2001Coefficients::C`                    | -                                                                                                                              | -                               |
| `rowbury_2001_k_floor`                          | set: -, theta_deg: deg, alpha_deg: deg                                                                                         | - (k_perp)                      |
| `rowbury_2001_in_range`                         | set: -, theta_deg: deg, alpha_deg: deg                                                                                         | - (bool)                        |
| `rowbury_2001_exit_crossflow_cd`                | Cd0: -, k_perp: -, set: -, theta_deg: deg, alpha_deg: deg                                                                      | Rowbury2001Cd                   |
| `Rowbury2001Cd::Cd`                             | -                                                                                                                              | -                               |
| `Rowbury2001Cd::dCd_dCd0`                       | -                                                                                                                              | -                               |
| `Rowbury2001Cd::dCd_dk`                         | -                                                                                                                              | -                               |
| `Rowbury2001Cd::delta`                          | -                                                                                                                              | -                               |
| `Rowbury2001Cd::k_used`                         | -                                                                                                                              | -                               |
| `Rowbury2001Cd::held`                           | -                                                                                                                              | - (bool)                        |
| `rowbury_2001_k_perp`                           | G_h: kg/(m^2*s), rho_h: kg/m^3, rho_m_u_m2: Pa, theta_deg: deg, alpha_deg: deg                                                 | Rowbury2001KPerp                |
| `rowbury_2001_exit_density`                     | Pt: Pa, Tt: K, P: Pa, X: mol/mol                                                                                               | Rowbury2001ExitDensity          |
| `Rowbury2001ExitDensity::rho`                   | kg/m^3                                                                                                                         | kg/m^3                          |
| `Rowbury2001ExitDensity::T`                     | K                                                                                                                              | K                               |
| `Rowbury2001ExitDensity::drho_dPt`              | kg/(m^3*Pa)                                                                                                                    | kg/(m^3*Pa)                     |
| `Rowbury2001ExitDensity::drho_dTt`              | kg/(m^3*K)                                                                                                                     | kg/(m^3*K)                      |
| `Rowbury2001ExitDensity::drho_dP`               | kg/(m^3*Pa)                                                                                                                    | kg/(m^3*Pa)                     |
| `Rowbury2001ExitDensity::drho_dY`               | kg/m^3                                                                                                                         | kg/m^3                          |
| `Rowbury2001Cd::faded`                          | -                                                                                                                              | - (bool)                        |
| `ROWBURY_2001_K_PERP_MAX`                       | -                                                                                                                              | -                               |
| `ROWBURY_2001_FADE_WIDTH`                       | -                                                                                                                              | - (ln k)                        |
| `Rowbury2001KPerp::k`                           | -                                                                                                                              | -                               |
| `Rowbury2001KPerp::dk_dG`                       | 1/(kg/(m^2*s))                                                                                                                 | 1/(kg/(m^2*s))                  |
| `Rowbury2001KPerp::dk_drho_h`                   | m^3/kg                                                                                                                         | m^3/kg                          |
| `Rowbury2001KPerp::dk_dq_m`                     | 1/Pa                                                                                                                           | 1/Pa                            |
| `Rowbury2001Set::CylindricalPressureSurface`    | -                                                                                                                              | -                               |
| `Rowbury2001Set::CylindricalLeadingEdge`        | -                                                                                                                              | -                               |
| `Rowbury2001Set::CylindricalSuctionSurface`     | -                                                                                                                              | -                               |
| `Rowbury2001Set::FanShapedPressureSurface`      | -                                                                                                                              | -                               |
| `Rowbury2001Set::FanShapedSuctionSurface`       | -                                                                                                                              | -                               |
| `Rowbury2001Set::name`                          | -                                                                                                                              | -                               |
| `Rowbury2001Set::value`                         | -                                                                                                                              | -                               |
| `ROWBURY_2001_K_PERP_MIN`                       | -                                                                                                                              | -                               |
| `ROWBURY_2001_UNIQUE_MARGIN`                    | -                                                                                                                              | -                               |
| `ROGERS_HERSH_1976_CAP_BAND`                    | -                                                                                                                              | -                               |
| `FLORSCHUETZ_METZGER_SU_1982_CAP_BAND`          | -                                                                                                                              | -                               |
| `ROGERS_HERSH_1976_INJECTION_A`                 | -                                                                                                                              | -                               |
| `ROGERS_HERSH_1976_INJECTION_B`                 | -                                                                                                                              | -                               |
| `ROGERS_HERSH_1976_SUCTION_A`                   | -                                                                                                                              | -                               |
| `ROGERS_HERSH_1976_SUCTION_B`                   | -                                                                                                                              | -                               |
| `STATION_KAPPA_MERGE_NORMAL`                    | -                                                                                                                              | -                               |
| `STATION_KAPPA_BLEED_BASSETT`                   | -                                                                                                                              | -                               |
| `WALL_HTC_KNEE`                                 | -                                                                                                                              | W/(m^2*K)                       |
| `MIXER_HEAT_MDOT_FLOOR`                         | -                                                                                                                              | kg/s                            |
| `mixer_heat_per_unit_flow`                      | mdot_tot: kg/s                                                                                                                 | phi: s/kg, dphi_dmdot: s^2/kg^2 |
| `SENSIBLE_ENTHALPY_REF_T`                       | -                                                                                                                              | K                               |
| `duct_static_state`                             | Tt: K, Pt: Pa, G: kg/(m^2*s), X: mol/mol                                                                                       | DuctStaticState                 |
| `DuctStaticState::T`                            | K                                                                                                                              | K                               |
| `DuctStaticState::P`                            | Pa                                                                                                                             | Pa                              |
| `DuctStaticState::rho`                          | kg/m^3                                                                                                                         | kg/m^3                          |
| `DuctStaticState::M`                            | -                                                                                                                              | -                               |
| `DuctStaticState::dT_dTt`                       | -                                                                                                                              | -                               |
| `DuctStaticState::dT_dPt`                       | K/Pa                                                                                                                           | K/Pa                            |
| `DuctStaticState::dT_dG`                        | K/(kg/(m^2*s))                                                                                                                 | K/(kg/(m^2*s))                  |
| `DuctStaticState::dP_dTt`                       | Pa/K                                                                                                                           | Pa/K                            |
| `DuctStaticState::dP_dPt`                       | -                                                                                                                              | -                               |
| `DuctStaticState::dP_dG`                        | Pa/(kg/(m^2*s))                                                                                                                | Pa/(kg/(m^2*s))                 |
| `DuctStaticState::drho_dTt`                     | kg/(m^3*K)                                                                                                                     | kg/(m^3*K)                      |
| `DuctStaticState::drho_dPt`                     | kg/(m^3*Pa)                                                                                                                    | kg/(m^3*Pa)                     |
| `DuctStaticState::drho_dG`                      | s/m                                                                                                                            | s/m                             |
| `DuctStaticState::choked`                       | -                                                                                                                              | -                               |
| `DuctStaticState::G_eff`                        | kg/(m^2*s)                                                                                                                     | kg/(m^2*s)                      |
| `DuctStaticState::dGeff_dTt`                    | kg/(m^2*s*K)                                                                                                                   | kg/(m^2*s*K)                    |
| `DuctStaticState::dGeff_dPt`                    | kg/(m^2*s*Pa)                                                                                                                  | kg/(m^2*s*Pa)                   |
| `DuctStaticState::dGeff_dG`                     | -                                                                                                                              | -                               |
| `duct_htc_static_state`                         | Tt: K, Pt: Pa, G: kg/(m^2*s), X: mol/mol                                                                                       | DuctStaticState                 |
| `DUCT_HTC_SONIC_MARGIN`                         | -                                                                                                                              | -                               |
| `duct_htc_static_state_dY`                      | Tt: K, Pt: Pa, G: kg/(m^2*s), X: mol/mol, state: DuctStaticState                                                               | DuctStaticStateDY               |
| `DuctStaticStateDY::dT_dY`                      | K (per species)                                                                                                                | K (per species)                 |
| `DuctStaticStateDY::dP_dY`                      | Pa (per species)                                                                                                               | Pa (per species)                |
| `DuctStaticStateDY::drho_dY`                    | kg/m^3 (per species)                                                                                                           | kg/m^3 (per species)            |
| `DuctStaticStateDY::dGeff_dY`                   | kg/(m^2*s) (per species)                                                                                                       | kg/(m^2*s) (per species)        |
| `chamber_merge_face_state`                      | m_main: kg/s, Tt_main: K, X_main: -, m_out: kg/s, P: Pa, Tt: K, X: -, side_momentum: N, area: m^2, incompressible: -           | MergeFaceState                  |
| `MergeFaceState::P_face`                        | Pa                                                                                                                             | Pa                              |
| `MergeFaceState::Pt_face`                       | Pa                                                                                                                             | Pa                              |
| `MergeFaceState::M_face`                        | -                                                                                                                              | -                               |
| `MergeFaceState::choked`                        | -                                                                                                                              | -                               |
| `MergeFaceState::dPf_dm_main`                   | Pa/(kg/s)                                                                                                                      | Pa/(kg/s)                       |
| `MergeFaceState::dPf_dm_out`                    | Pa/(kg/s)                                                                                                                      | Pa/(kg/s)                       |
| `MergeFaceState::dPf_dJ`                        | 1/m^2                                                                                                                          | 1/m^2                           |
| `MergeFaceState::dPf_dP`                        | -                                                                                                                              | -                               |
| `MergeFaceState::dPf_dT`                        | Pa/K                                                                                                                           | Pa/K                            |
| `MergeFaceState::dPf_dT_main`                   | Pa/K                                                                                                                           | Pa/K                            |
| `MergeFaceState::dPtf_dm_main`                  | Pa/(kg/s)                                                                                                                      | Pa/(kg/s)                       |
| `MergeFaceState::dPtf_dm_out`                   | Pa/(kg/s)                                                                                                                      | Pa/(kg/s)                       |
| `MergeFaceState::dPtf_dJ`                       | 1/m^2                                                                                                                          | 1/m^2                           |
| `MergeFaceState::dPtf_dP`                       | -                                                                                                                              | -                               |
| `MergeFaceState::dPtf_dT`                       | Pa/K                                                                                                                           | Pa/K                            |
| `MergeFaceState::dPtf_dT_main`                  | Pa/K                                                                                                                           | Pa/K                            |
| `MergeFaceState::dPf_dY`                        | Pa (per species)                                                                                                               | Pa (per species)                |
| `MergeFaceState::dPf_dY_main`                   | Pa (per species)                                                                                                               | Pa (per species)                |
| `MergeFaceState::dPtf_dY`                       | Pa (per species)                                                                                                               | Pa (per species)                |
| `MergeFaceState::dPtf_dY_main`                  | Pa (per species)                                                                                                               | Pa (per species)                |
| `jet_impulse`                                   | m: kg/s, Pt: Pa, Tt: K, P: Pa, X: -                                                                                            | JetImpulse                      |
| `JetImpulse::J`                                 | N                                                                                                                              | N                               |
| `JetImpulse::dJ_dm`                             | m/s                                                                                                                            | m/s                             |
| `JetImpulse::dJ_dPt`                            | N/Pa                                                                                                                           | N/Pa                            |
| `JetImpulse::dJ_dTt`                            | N/K                                                                                                                            | N/K                             |
| `JetImpulse::dJ_dP`                             | N/Pa                                                                                                                           | N/Pa                            |
| `JetImpulse::choked`                            | -                                                                                                                              | -                               |
| `FLORSCHUETZ_1981_DEFAULT_CD`                   | -                                                                                                                              | -                               |
| `FLORSCHUETZ_SU_1985_ETA_A`                     | -                                                                                                                              | -                               |
| `FLORSCHUETZ_SU_1985_ETA_R`                     | -                                                                                                                              | -                               |
| `FLORSCHUETZ_SU_1985_ETA_N`                     | -                                                                                                                              | -                               |
| `SingleJetImpingementSet::name`                 | -                                                                                                                              | -                               |
| `SingleJetImpingementSet::source`               | -                                                                                                                              | -                               |
| `SingleJetImpingementSet::A`                    | -                                                                                                                              | -                               |
| `SingleJetImpingementSet::B`                    | -                                                                                                                              | -                               |
| `SingleJetImpingementSet::C`                    | -                                                                                                                              | -                               |
| `SingleJetImpingementSet::Re_exponent`          | -                                                                                                                              | -                               |
| `SingleJetImpingementSet::n_const_heat_flux`    | -                                                                                                                              | -                               |
| `SingleJetImpingementSet::n_const_wall_temp`    | -                                                                                                                              | -                               |
| `JetArrayGeometricFit::C`                       | -                                                                                                                              | -                               |
| `JetArrayGeometricFit::nx`                      | -                                                                                                                              | -                               |
| `JetArrayGeometricFit::ny`                      | -                                                                                                                              | -                               |
| `JetArrayGeometricFit::nz`                      | -                                                                                                                              | -                               |
| `JetArrayCorrelationSet::name`                  | -                                                                                                                              | -                               |
| `JetArrayCorrelationSet::source`                | -                                                                                                                              | -                               |
| `JetArrayCorrelationSet::validity_source`       | -                                                                                                                              | -                               |
| `JetArrayCorrelationSet::pattern`               | -                                                                                                                              | -                               |
| `JetArrayCorrelationSet::A_fit`                 | -                                                                                                                              | -                               |
| `JetArrayCorrelationSet::m_fit`                 | -                                                                                                                              | -                               |
| `JetArrayCorrelationSet::B_fit`                 | -                                                                                                                              | -                               |
| `JetArrayCorrelationSet::n_fit`                 | -                                                                                                                              | -                               |
| `JetArrayCorrelationSet::valid_Re_j`            | -                                                                                                                              | -                               |
| `JetArrayCorrelationSet::valid_Gc_Gj`           | -                                                                                                                              | -                               |
| `JetArrayCorrelationSet::valid_xn_d`            | -                                                                                                                              | -                               |
| `JetArrayCorrelationSet::valid_yn_d`            | -                                                                                                                              | -                               |
| `JetArrayCorrelationSet::valid_z_d`             | -                                                                                                                              | -                               |
| `JetArrayCorrelationSet::valid_aspect_ratio`    | -                                                                                                                              | -                               |
| `JetArrayCorrelationSet::standard_error`        | -                                                                                                                              | -                               |
| `JetArrayImpingementResult::Nu`                 | -                                                                                                                              | -                               |
| `JetArrayImpingementResult::dNu_dRe_j`          | -                                                                                                                              | -                               |
| `JetArrayImpingementResult::dNu_dGc_Gj`         | -                                                                                                                              | -                               |
| `JetArrayImpingementResult::extrapolated`       | -                                                                                                                              | -                               |
| `SingleJetImpingementResult::Nu`                | -                                                                                                                              | -                               |
| `SingleJetImpingementResult::dNu_dRe`           | -                                                                                                                              | -                               |
| `ImpingementRange::lo`                          | -                                                                                                                              | -                               |
| `ImpingementRange::hi`                          | -                                                                                                                              | -                               |
| `ImpingementThermalBC::ConstantHeatFlux`        | -                                                                                                                              | -                               |
| `ImpingementThermalBC::ConstantWallTemperature` | -                                                                                                                              | -                               |
| `ImpingementThermalBC::name`                    | -                                                                                                                              | -                               |
| `ImpingementThermalBC::value`                   | -                                                                                                                              | -                               |
| `JetHolePattern::Inline`                        | -                                                                                                                              | -                               |
| `JetHolePattern::Staggered`                     | -                                                                                                                              | -                               |
| `JetHolePattern::name`                          | -                                                                                                                              | -                               |
| `JetHolePattern::value`                         | -                                                                                                                              | -                               |

### acoustics.h - Acoustic Properties

| Function                    | Input Units                                                                      | Output Unit               |
|-----------------------------|----------------------------------------------------------------------------------|---------------------------|
| `acoustic_properties`       | f: Hz, rho: kg/m^3, c: m/s, p_rms: Pa (default 20e-6), p_ref: Pa (default 20e-6) | AcousticProperties struct |
| `wavelength`                | f: Hz, c: m/s                                                                    | m                         |
| `frequency_from_wavelength` | lambda: m, c: m/s                                                                | Hz                        |
| `acoustic_impedance`        | rho: kg/m^3, c: m/s                                                              | Pa*s/m                    |
| `sound_pressure_level`      | p_rms: Pa, p_ref: Pa (default 20e-6)                                             | dB                        |
| `particle_velocity`         | p: Pa, rho: kg/m^3, c: m/s                                                       | m/s                       |

### orifice.h - Orifice Flow Utilities

| Function                     | Input Units                                                              | Output Unit              |
|------------------------------|--------------------------------------------------------------------------|--------------------------|
| `orifice_flow`               | geom: OrificeGeometry, dP: Pa, T: K, P: Pa, mu: Pa*s, Z: - (default 1.0) | OrificeFlowResult struct |
| `orifice_velocity_from_mdot` | mdot: kg/s, rho: kg/m^3, d: m, Z: - (default 1.0)                        | m/s                      |
| `orifice_area_from_beta`     | D: m, beta: -                                                            | m^2                      |
| `beta_from_diameters`        | d: m, D: m                                                               | -                        |
| `orifice_Re_d_from_mdot`     | mdot: kg/s, d: m, mu: Pa*s                                               | - (Reynolds number)      |

### solver_interface.h - Compressible Flow Elements

| Function                                      | Input Units                                                                                                              | Output Unit                                 |
|-----------------------------------------------|--------------------------------------------------------------------------------------------------------------------------|---------------------------------------------|
| `orifice_compressible_mdot_and_jacobian`      | Tt: K, Pt: Pa, P_back: Pa, X: mol/mol, Cd: -, area: m^2, beta: -                                                         | tuple(kg/s, kg/(s*Pa), kg/(s*Pa), kg/(s*K)) |
| `orifice_compressible_residuals_and_jacobian` | m_dot: kg/s, P_total_up: Pa, P_static_up: Pa, Tt_up: K, Y_up: kg/kg, P_static_down: Pa, Cd: -, area: m^2, beta: -        | OrificeResult                               |
| `orifice_jobson_residuals_and_jacobian`       | m_dot: kg/s, P_total_up: Pa, P_static_up: Pa, Tt_up: K, Y_up: kg/kg, P_static_down: Pa, Cd: -, area: m^2, beta: -, Cv: - | OrificeModelResult                          |
| `orifice_iso5167_residuals_and_jacobian`      | m_dot: kg/s, P_total_up: Pa, P_static_up: Pa, Tt_up: K, Y_up: kg/kg, P_static_down: Pa, Cd: -, area: m^2, beta: -        | OrificeModelResult                          |
| `sonic_static_pressure`                       | Tt: K, G: kg/(m^2*s), Y: kg/kg                                                                                           | SonicStatic                                 |
| `approach_static_pressure`                    | Pt: Pa, Tt: K, G: kg/(m^2*s), Y: kg/kg, incompressible: bool                                                             | ApproachStatic                              |
| `orifice_tap_pressure_and_jacobian`           | P_static_up: Pa, P_total_up: Pa, P_total_down: Pa, beta: -, C: -, compressible: bool                                     | OrificeTapPressure                          |

### geometry.h - Geometric Utilities

| Function                     | Input Units            | Output Unit  |
|------------------------------|------------------------|--------------|
| `circular_area`              | D: m                   | m^2          |
| `annular_area`               | D_outer: m, D_inner: m | m^2          |
| `cylinder_volume`            | D: m, L: m             | m^3          |
| `channel_roughness`          | material: str          | m            |
| `standard_channel_roughness` | -                      | dict[str, m] |
| `hydraulic_diameter`         | A: m^2, P_wetted: m    | m            |
| `hydraulic_diameter_rect`    | a: m, b: m             | m            |
| `hydraulic_diameter_annulus` | D_outer: m, D_inner: m | m            |

### area_change.h - Area Change Elements (Idelchik 1966, #731)

| Function                     | Input Units                                                                                                                                         | Output Unit                                                                                                                                                                                |
|------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `sudden_expansion_zeta`      | Re: -, area_ratio: -                                                                                                                                | tuple (zeta: -, dzeta_dRe: -)                                                                                                                                                              |
| `sudden_contraction_zeta`    | Re: -, area_ratio: -                                                                                                                                | tuple (zeta: -, dzeta_dRe: -)                                                                                                                                                              |
| `conical_diffuser_zeta`      | full_angle_deg: deg, area_ratio: -                                                                                                                  | -                                                                                                                                                                                          |
| `conical_nozzle_zeta`        | full_angle_deg: deg, length_over_D: -, area_ratio: -                                                                                                | -                                                                                                                                                                                          |
| `compressible_borda_carnot`  | M1: -, area_ratio: -, gamma: -                                                                                                                      | BordaCarnotStep (M2: -, pt_ratio: -, dpt_dM1: -, dpt_dgamma: -)                                                                                                                            |
| `area_change_loss`           | m_dot_abs: kg/s, Pt0: Pa, Tt0: K, Y: -, F_in: m^2, F_out: m^2, model: -, length: m, D_h: m, incompressible: bool                                    | AreaChangeLoss (dPt: Pa, zeta: -, mach: -, Re: -, expansion: bool, choked: bool)                                                                                                           |
| `area_change_flow`           | Pt0: Pa, Tt0: K, Y: -, Pt_target: Pa, F_in: m^2, F_out: m^2, model: -, length: m, D_h: m, incompressible: bool, with_derivatives: bool              | AreaChangeFlow (m_dot: kg/s, choked: bool, mach: -, Re: -, zeta: -, dm_dPt0: kg/(s*Pa), dm_dTt0: kg/(s*K), dm_dPt_target: kg/(s*Pa), dm_dY: kg/s)                                          |
| `area_change_element_flow`   | Pt_a: Pa, Tt_a: K, Y_a: -, Pt_b: Pa, Tt_b: K, Y_b: -, F_a: m^2, F_b: m^2, model: -, length: m, D_h: m, incompressible: bool, with_derivatives: bool | AreaChangeElementFlow (m_dot: kg/s, choked: bool, mach: -, Re: -, zeta: -, dm_dPt_a: kg/(s*Pa), dm_dTt_a: kg/(s*K), dm_dPt_b: kg/(s*Pa), dm_dTt_b: kg/(s*K), dm_dY_a: kg/s, dm_dY_b: kg/s) |
| `compression`                | T1: K, X: mol/mol, PR: -, eta_poly: -, eta_is: -                                                                                                    | TurbomachinePath (T2: K, dh: J/kg, eta_poly: -, eta_is: -, T2s: K, dT2_dT1: -, dT2_dPR: K, dT2_deta: K, dT2_dY: K, ddh_dT1: J/(kg*K), ddh_dPR: J/kg, ddh_deta: J/kg, ddh_dY: J/kg)         |
| `expansion`                  | T1: K, X: mol/mol, PR: -, eta_poly: -, eta_is: -                                                                                                    | TurbomachinePath (T2: K, dh: J/kg, eta_poly: -, eta_is: -, T2s: K, dT2_dT1: -, dT2_dPR: K, dT2_deta: K, dT2_dY: K, ddh_dT1: J/(kg*K), ddh_dPR: J/kg, ddh_deta: J/kg, ddh_dY: J/kg)         |
| `TurbomachinePath::T2`       | -                                                                                                                                                   | K                                                                                                                                                                                          |
| `TurbomachinePath::dh`       | -                                                                                                                                                   | J/kg                                                                                                                                                                                       |
| `TurbomachinePath::eta_poly` | -                                                                                                                                                   | -                                                                                                                                                                                          |
| `TurbomachinePath::eta_is`   | -                                                                                                                                                   | -                                                                                                                                                                                          |
| `TurbomachinePath::T2s`      | -                                                                                                                                                   | K                                                                                                                                                                                          |
| `TurbomachinePath::dT2_dT1`  | -                                                                                                                                                   | -                                                                                                                                                                                          |
| `TurbomachinePath::dT2_dPR`  | -                                                                                                                                                   | K                                                                                                                                                                                          |
| `TurbomachinePath::dT2_deta` | -                                                                                                                                                   | K                                                                                                                                                                                          |
| `TurbomachinePath::dT2_dY`   | -                                                                                                                                                   | K                                                                                                                                                                                          |
| `TurbomachinePath::ddh_dT1`  | -                                                                                                                                                   | J/(kg*K)                                                                                                                                                                                   |
| `TurbomachinePath::ddh_dPR`  | -                                                                                                                                                   | J/kg                                                                                                                                                                                       |
| `TurbomachinePath::ddh_deta` | -                                                                                                                                                   | J/kg                                                                                                                                                                                       |
| `TurbomachinePath::ddh_dY`   | -                                                                                                                                                   | J/kg                                                                                                                                                                                       |
| `sharp_area_change`          | m_dot: kg/s, rho: kg/m^3, mu: Pa*s, F0: m^2, F1: m^2, m_scale: kg/s, D_h: m                                                                         | AreaChangeResult (dP: Pa, dS_dm: Pa*s/kg, dS_drho: Pa*m^3/kg, dS_dmu: Pa/(Pa*s))                                                                                                           |
| `conical_area_change`        | m_dot: kg/s, rho: kg/m^3, mu: Pa*s, F0: m^2, F1: m^2, length: m, m_scale: kg/s                                                                      | AreaChangeResult (dP: Pa, dS_dm: Pa*s/kg, dS_drho: Pa*m^3/kg, dS_dmu: Pa/(Pa*s))                                                                                                           |

### tee_junction.h + solver_interface.h - Tee Junction Elements

| Function                               | Input Units                                                                                                                                                              | Output Unit                                                                                                                                                                                             |
|----------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `merging_tee_residuals_and_jacobian`   | m_dot_com: kg/s, m_dot_branch: kg/s, dP0_straight: Pa, dP0_branch: Pa, P_static_com: Pa, T_com: K, Y_com: kg/kg, theta: rad, psi: -, F_C: m^2, blend_k: - (default 30.0) | TeeJunctionResult (R_straight: Pa, R_branch: Pa, dR/dm_dot: Pa*s/kg, dR/dP: -, dR/dT: Pa/K, dR/dY: Pa, K_straight: -, K_branch: -, q: -, blend_w: -, topology_valid: bool, status: CorrelationValidity) |
| `branching_tee_residuals_and_jacobian` | m_dot_com: kg/s, m_dot_branch: kg/s, dP0_straight: Pa, dP0_branch: Pa, P_static_com: Pa, T_com: K, Y_com: kg/kg, theta: rad, psi: -, F_C: m^2, blend_k: - (default 30.0) | TeeJunctionResult (R_straight: Pa, R_branch: Pa, dR/dm_dot: Pa*s/kg, dR/dP: -, dR/dT: Pa/K, dR/dY: Pa, K_straight: -, K_branch: -, q: -, blend_w: -, topology_valid: bool, status: CorrelationValidity) |
| `tee_K1`                               | q: -, psi: -, theta: rad                                                                                                                                                 | - (K)                                                                                                                                                                                                   |
| `tee_K2`                               | q: -                                                                                                                                                                     | - (K)                                                                                                                                                                                                   |
| `tee_K3`                               | q: -, psi: -, theta: rad                                                                                                                                                 | - (K)                                                                                                                                                                                                   |
| `tee_K4`                               | q: -, psi: -, theta: rad                                                                                                                                                 | - (K)                                                                                                                                                                                                   |
| `tee_K5`                               | q: -                                                                                                                                                                     | - (K)                                                                                                                                                                                                   |
| `tee_K6`                               | q: -, psi: -, theta: rad                                                                                                                                                 | - (K)                                                                                                                                                                                                   |
| `tee_K7`                               | q: -, psi: -, theta: rad                                                                                                                                                 | - (K)                                                                                                                                                                                                   |
| `tee_K8`                               | q: -, psi: -, theta: rad                                                                                                                                                 | - (K)                                                                                                                                                                                                   |
| `tee_K9`                               | q: -, psi: -, theta: rad                                                                                                                                                 | - (K)                                                                                                                                                                                                   |
| `tee_K10`                              | q: -, psi: -, theta: rad                                                                                                                                                 | - (K)                                                                                                                                                                                                   |
| `tee_K11`                              | q: -, psi: -, theta: rad                                                                                                                                                 | - (K)                                                                                                                                                                                                   |
| `tee_K12`                              | q: -, psi: -, theta: rad                                                                                                                                                 | - (K)                                                                                                                                                                                                   |

### mpce_junction.h - Momentum-CV junction, whole-element (f, J)

| Function                      | Input Units                                                                                                                                                                                                                                                                                                                       | Output Unit                                                                                                                                                                                                                              |
|-------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `mpce_residuals_and_jacobian` | p_static: list[Pa], p_total: list[Pa], rho: list[kg/m^3], drho_dp: list[kg/(m^3*Pa)], outer_mdot: list[kg/s] (the CONNECTING element's flow; geom.port_sign maps it to junction convention), pt_jct: Pa, geom: MpceGeometry (area: list[m^2], theta_rad: list[rad], port_sign: list[-], joining_etransfer_alpha: -, eta_scale: -) | MpceResidualJacobian (residual: list[Pa] with a final kg/s mass row, jacobian: list[list[d(row)/d(seed)]] over seeds [p_static x3, p_total x3, outer_mdot x3, pt_jct], valid: bool, common_port: -, k_term_sign: -, k_per_port: list[-]) |

### ideal_junction.h - Ideal momentum junction, whole-element (f, J) (#517)

| Function                                | Input Units                                                                                                                                                                                                                                                                                                                                                          | Output Unit                                                                                                                                                                                                 |
|-----------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `ideal_junction_residuals_and_jacobian` | p_static: list[Pa], rho: list[kg/m^3], drho_dp: list[kg/(m^3*Pa)], outer_mdot: list[kg/s] (the CONNECTING element's flow; geom.port_sign maps it to junction convention, + out), p_jn: Pa, geom: IdealJunctionGeometry (area: list[m^2], axis_rad: list[rad] outward port axis, port_sign: list[-], loss_factor: list[-] receiver loss factor f >= 0, empty = all 1) | IdealJunctionResult (residual: list[Pa] with a final kg/s mass row, jacobian: list[list[d(row)/d(col)]] over columns [p_static x N, outer_mdot x N, p_jn], velocity: list[m/s] (+ out), V_x: m/s, V_y: m/s) |

### vatistas.h - Vatistas n-Vortex Model

| Function                             | Input Units                                   | Output Unit                             |
|--------------------------------------|-----------------------------------------------|-----------------------------------------|
| `vatistas_v0_bar`                    | r_bar: -, n: -                                | - (V0_bar)                              |
| `vatistas_dv0_bar_drbar`             | r_bar: -, n: -                                | -                                       |
| `vatistas_vr_bar`                    | r_bar: -, n: -                                | - (Vr_bar)                              |
| `vatistas_dvr_bar_drbar`             | r_bar: -, n: -                                | -                                       |
| `vatistas_pressure_integral`         | r_bar: -, n: -                                | - (integral)                            |
| `vatistas_d_pressure_integral_drbar` | r_bar: -, n: -                                | -                                       |
| `vatistas_v_theta`                   | r: m, Gamma: m^2/s, r_c: m, n: -              | m/s                                     |
| `vatistas_v_theta_and_jacobians`     | r: m, Gamma: m^2/s, r_c: m, n: -              | tuple(m/s, 1/s, (m/s)/(m^2/s), (m/s)/m) |
| `vatistas_delta_p`                   | r: m, rho: kg/m^3, Gamma: m^2/s, r_c: m, n: - | Pa                                      |
| `vatistas_delta_p_and_jacobians`     | r: m, rho: kg/m^3, Gamma: m^2/s, r_c: m, n: - | tuple(Pa, Pa/m, Pa/(m^2/s), Pa/m)       |
| `VatistasVortex::V_theta`            | r: m                                          | m/s                                     |
| `VatistasVortex::V_theta_max`        | -                                             | m/s                                     |
| `VatistasVortex::dV_theta_dr`        | r: m                                          | 1/s                                     |
| `VatistasVortex::delta_P`            | r: m, rho: kg/m^3                             | Pa                                      |
| `VatistasVortex::delta_P_bar`        | r: m                                          | -                                       |

### heat_transfer.h - Heat Transfer Correlations

| Function                                                    | Input Units                                                                                                                                         | Output Unit                                                       |
|-------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------|
| `nusselt_dittus_boelter`                                    | Re: -, Pr: -, heating: bool                                                                                                                         | - (Nu)                                                            |
| `nusselt_gnielinski`                                        | Re: -, Pr: -, f: - (optional)                                                                                                                       | - (Nu)                                                            |
| `nusselt_sieder_tate`                                       | Re: -, Pr: -, mu_ratio: -                                                                                                                           | - (Nu)                                                            |
| `nusselt_petukhov`                                          | Re: -, Pr: -, f: - (optional)                                                                                                                       | - (Nu)                                                            |
| `htc_from_nusselt`                                          | Nu: -, k: W/(m*K), L: m                                                                                                                             | W/(m^2*K)                                                         |
| `nusselt_circular_channel`                                  | State, V: m/s, D: m                                                                                                                                 | - (Nu)                                                            |
| `htc_circular_channel`                                      | State, V: m/s, D: m                                                                                                                                 | W/(m^2*K)                                                         |
| `htc_circular_channel`                                      | T: K, P: Pa, X: mol/mol, velocity: m/s, diameter: m, correlation: str, heating: bool, mu_ratio: -, roughness: m                                     | tuple(W/(m^2*K), -, -)                                            |
| `lmtd`                                                      | dT1: K, dT2: K                                                                                                                                      | K                                                                 |
| `lmtd_counterflow`                                          | T_hot_in: K, T_hot_out: K, T_cold_in: K, T_cold_out: K                                                                                              | K                                                                 |
| `lmtd_parallelflow`                                         | T_hot_in: K, T_hot_out: K, T_cold_in: K, T_cold_out: K                                                                                              | K                                                                 |
| `overall_htc`                                               | h_values: W/(m^2*K), t_over_k: m^2*K/W                                                                                                              | W/(m^2*K)                                                         |
| `overall_htc_wall`                                          | h_inner, h_outer: W/(m^2*K), t_over_k_layers: m^2*K/W                                                                                               | W/(m^2*K)                                                         |
| `overall_htc_wall`                                          | h_inner, h_outer: W/(m^2*K), t_wall: m, k_wall: W/(m*K)                                                                                             | W/(m^2*K)                                                         |
| `thermal_resistance`                                        | h: W/(m^2*K), A: m^2                                                                                                                                | K/W                                                               |
| `thermal_resistance_wall`                                   | t: m, k: W/(m*K), A: m^2                                                                                                                            | K/W                                                               |
| `heat_rate`                                                 | U: W/(m^2*K), A: m^2, dT: K                                                                                                                         | W                                                                 |
| `heat_flux`                                                 | U: W/(m^2*K), dT: K                                                                                                                                 | W/m^2                                                             |
| `heat_transfer_area`                                        | Q: W, U: W/(m^2*K), dT: K                                                                                                                           | m^2                                                               |
| `heat_transfer_dT`                                          | Q: W, U: W/(m^2*K), A: m^2                                                                                                                          | K                                                                 |
| `wall_temperature_profile`                                  | T_hot: K, T_cold: K, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W                                                                                    | K (vector)                                                        |
| `wall_resistance_fixed_point_derivatives`                   | T_hot: K, T_cold: K, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W, R_fouling: m^2*K/W, dr_dTm: m^2/W                                                 | dR_wall/d(T_hot, T_cold) [m^2/W], /d(h_hot, h_cold) [m^4*K^2/W^2] |
| `WallCouplingResult::dQ_dR_wall`                            | -                                                                                                                                                   | W^2/(m^2*K)                                                       |
| `STEFAN_BOLTZMANN`                                          | -                                                                                                                                                   | W/(m^2*K^4)                                                       |
| `REEVES_EMISSIVITY_COEFF`                                   | -                                                                                                                                                   | kPa^-1*m^-0.5*K^1.5                                               |
| `LUMINOSITY_HYDROGEN_COEFF`                                 | -                                                                                                                                                   | - (L = coeff / H[%]^2)                                            |
| `BEAM_LENGTH_COEFF_LEFEBVRE`                                | -                                                                                                                                                   | -                                                                 |
| `BEAM_LENGTH_COEFF_HOWELL`                                  | -                                                                                                                                                   | -                                                                 |
| `gas_emissivity_luminous`                                   | T_g: K, P: Pa, q: -, l_b: m, L: -                                                                                                                   | GasEmissivityResult                                               |
| `GasEmissivityResult`                                       | -                                                                                                                                                   | eps: -, partials                                                  |
| `GasEmissivityResult::eps`                                  | -                                                                                                                                                   | -                                                                 |
| `GasEmissivityResult::deps_dT`                              | -                                                                                                                                                   | 1/K                                                               |
| `GasEmissivityResult::deps_dP`                              | -                                                                                                                                                   | 1/Pa                                                              |
| `GasEmissivityResult::deps_dq`                              | -                                                                                                                                                   | -                                                                 |
| `GasEmissivityResult::deps_dl_b`                            | -                                                                                                                                                   | 1/m                                                               |
| `GasEmissivityResult::deps_dL`                              | -                                                                                                                                                   | -                                                                 |
| `luminosity_factor_hydrogen`                                | H_mass_percent: %                                                                                                                                   | L: -, dL_dH: 1/%                                                  |
| `liner_gas_radiation`                                       | T_g: K, T_w: K, eps_g: -, eps_w: -                                                                                                                  | LinerGasRadiationResult                                           |
| `LinerGasRadiationResult`                                   | -                                                                                                                                                   | R1: W/m^2, partials                                               |
| `LinerGasRadiationResult::R1`                               | -                                                                                                                                                   | W/m^2                                                             |
| `LinerGasRadiationResult::dR1_dT_g`                         | -                                                                                                                                                   | W/(m^2*K)                                                         |
| `LinerGasRadiationResult::dR1_dT_w`                         | -                                                                                                                                                   | W/(m^2*K)                                                         |
| `LinerGasRadiationResult::dR1_deps_g`                       | -                                                                                                                                                   | W/m^2                                                             |
| `LinerGasRadiationResult::dR1_deps_w`                       | -                                                                                                                                                   | W/m^2                                                             |
| `liner_casing_radiation`                                    | T_w: K, T_c: K, eps_w: -, eps_c: -, area_ratio: -                                                                                                   | R2: W/m^2, dR2_dT_w: W/(m^2*K), dR2_dT_c: W/(m^2*K)               |
| `liner_casing_exchange_factor`                              | eps_w: -, eps_c: -, area_ratio: -                                                                                                                   | W/(m^2*K^4)                                                       |
| `mean_beam_length`                                          | volume: m^3, area: m^2, coeff: -                                                                                                                    | m                                                                 |
| `elemental_fuel_air_ratio`                                  | Y: kg/kg                                                                                                                                            | q: -, dq_dY: -                                                    |
| `wall_radiation_sides`                                      | T_aw_a, T_aw_b, T_src_a, T_src_b: K, h_a, h_b: W/(m^2*K), R_wall: m^2*K/W, k_a, k_b: W/(m^2*K^4)                                                    | WallRadiationSides                                                |
| `wall_adiabatic_face`                                       | T_aw, T_src: K, h: W/(m^2*K), kind: -, k: W/(m^2*K^4)                                                                                               | WallAdiabaticFace                                                 |
| `ntu`                                                       | U: W/(m^2*K), A: m^2, C_min: W/K                                                                                                                    | -                                                                 |
| `capacity_ratio`                                            | C_min: W/K, C_max: W/K                                                                                                                              | -                                                                 |
| `effectiveness_counterflow`                                 | NTU: -, C_r: -                                                                                                                                      | -                                                                 |
| `effectiveness_parallelflow`                                | NTU: -, C_r: -                                                                                                                                      | -                                                                 |
| `heat_rate_from_effectiveness`                              | epsilon: -, C_min: W/K, T_hot_in: K, T_cold_in: K                                                                                                   | W                                                                 |
| `adiabatic_wall_temperature`                                | T_hot: K, T_coolant: K, eta: -                                                                                                                      | K                                                                 |
| `film_superposition_sellers`                                | eta_rows: -                                                                                                                                         | - (eta)                                                           |
| `effusion_panel_film_effectiveness`                         | n_rows: -, pitch_x_over_D: -, s_over_D: -, M: -, density_ratio: -, alpha_deg: deg, Tu: -                                                            | EffusionPanelFilm                                                 |
| `EffusionPanelFilm::eta`                                    | -                                                                                                                                                   | -                                                                 |
| `EffusionPanelFilm::extrapolated`                           | -                                                                                                                                                   | -                                                                 |
| `film_span_average_baldauf_2002`                            | n_rows: -, row_pitch_over_D: -, x_start_over_D: -, x_end_over_D: -, M: -, density_ratio: -, alpha_deg: deg, s_over_D: -, Tu: -                      | FilmSpanAverage                                                   |
| `film_span_average_colban_2011`                             | n_rows: -, row_pitch_over_D: -, x_start_over_D: -, x_end_over_D: -, M: -, area_ratio: -, P_over_D: -, t_over_P: -, alpha_deg: deg, density_ratio: - | FilmSpanAverage                                                   |
| `SEBAN1960_C_NEAR`                                          | -                                                                                                                                                   | - (curve B as drawn)                                              |
| `SEBAN1960_C_NEAR_PRINTED`                                  | -                                                                                                                                                   | - (Eq. 13 as printed)                                             |
| `SEBAN1960_C_FAR`                                           | -                                                                                                                                                   | - (curve C as drawn)                                              |
| `SEBAN1960_C_FAR_PRINTED`                                   | -                                                                                                                                                   | - (Eq. 14 as printed)                                             |
| `SLOT_FILM_HTC_HANDBACK_BAND`                               | -                                                                                                                                                   | - (mass-velocity ratio)                                           |
| `SLOT_FILM_BRANCH_NEAR`                                     | -                                                                                                                                                   | - (code)                                                          |
| `SLOT_FILM_BRANCH_BLEND`                                    | -                                                                                                                                                   | - (code)                                                          |
| `SLOT_FILM_BRANCH_FAR`                                      | -                                                                                                                                                   | - (code)                                                          |
| `SLOT_FILM_BRANCH_WALL_JET`                                 | -                                                                                                                                                   | - (code)                                                          |
| `SLOT_FILM_BRANCH_NONE`                                     | -                                                                                                                                                   | - (code)                                                          |
| `slot_film_stanton_seban_1960`                              | x_start_over_s: -, x_end_over_s: -, Re_s: -, C_near: -, C_far: -                                                                                    | SlotFilmStanton                                                   |
| `SlotFilmStanton`                                           | -                                                                                                                                                   | St: -, partials, box state                                        |
| `SlotFilmStanton::St`                                       | -                                                                                                                                                   | -                                                                 |
| `SlotFilmStanton::dSt_dRe`                                  | -                                                                                                                                                   | -                                                                 |
| `SlotFilmStanton::fraction_in_box`                          | -                                                                                                                                                   | -                                                                 |
| `SlotFilmStanton::branch`                                   | -                                                                                                                                                   | - (code)                                                          |
| `SlotFilmStanton::x_in_box`                                 | -                                                                                                                                                   | -                                                                 |
| `SlotFilmStanton::Re_in_box`                                | -                                                                                                                                                   | -                                                                 |
| `slot_film_nusselt_ballal_lefebvre_1973`                    | x_start_over_s: -, x_end_over_s: -, Re_s: -                                                                                                         | SlotFilmNusselt                                                   |
| `SlotFilmNusselt`                                           | -                                                                                                                                                   | Nu_s: -, partials, box state                                      |
| `SlotFilmNusselt::Nu_s`                                     | -                                                                                                                                                   | -                                                                 |
| `SlotFilmNusselt::dNu_dRe`                                  | -                                                                                                                                                   | -                                                                 |
| `SlotFilmNusselt::fraction_in_box`                          | -                                                                                                                                                   | -                                                                 |
| `SlotFilmNusselt::branch`                                   | -                                                                                                                                                   | - (code)                                                          |
| `SlotFilmNusselt::x_in_box`                                 | -                                                                                                                                                   | -                                                                 |
| `slot_film_htc_seban_1960`                                  | x_start, x_end, s: m, G_s: kg/(m^2*s), M: -, h_duct: W/(m^2*K), T_s: K, P_s: Pa, X_s: mol/mol                                                       | SlotFilmHtc                                                       |
| `slot_film_htc_ballal_lefebvre_1973`                        | x_start, x_end, s: m, G_s: kg/(m^2*s), M: -, h_duct: W/(m^2*K), T_s: K, P_s: Pa, X_s: mol/mol                                                       | SlotFilmHtc                                                       |
| `SlotFilmHtc`                                               | -                                                                                                                                                   | h: W/(m^2*K), partials, box state                                 |
| `SlotFilmHtc::h`                                            | -                                                                                                                                                   | W/(m^2*K)                                                         |
| `SlotFilmHtc::dh_dG_s`                                      | -                                                                                                                                                   | W*s/(kg*K)                                                        |
| `SlotFilmHtc::dh_dT_s`                                      | -                                                                                                                                                   | W/(m^2*K^2)                                                       |
| `SlotFilmHtc::dh_dM`                                        | -                                                                                                                                                   | W/(m^2*K)                                                         |
| `SlotFilmHtc::dh_dh_duct`                                   | -                                                                                                                                                   | -                                                                 |
| `SlotFilmHtc::h_film`                                       | -                                                                                                                                                   | W/(m^2*K)                                                         |
| `SlotFilmHtc::Re_s`                                         | -                                                                                                                                                   | -                                                                 |
| `SlotFilmHtc::fraction_in_box`                              | -                                                                                                                                                   | -                                                                 |
| `SlotFilmHtc::weight_M`                                     | -                                                                                                                                                   | -                                                                 |
| `SlotFilmHtc::branch`                                       | -                                                                                                                                                   | - (code)                                                          |
| `SlotFilmHtc::x_in_box`                                     | -                                                                                                                                                   | -                                                                 |
| `SlotFilmHtc::Re_in_box`                                    | -                                                                                                                                                   | -                                                                 |
| `SlotFilmHtc::M_in_box`                                     | -                                                                                                                                                   | -                                                                 |
| `slot_film_span_average_marek_tacina_1975`                  | x_start_over_S: -, x_end_over_S: -, M: -, Tu: -, C_m_multiplier: -                                                                                  | FilmSpanAverage                                                   |
| `film_span_average_baldauf_2002_groups`                     | n_rows: -, row_pitch_over_D: -, x_start_over_D: -, x_end_over_D: -, M: -, density_ratio: -, alpha_deg: deg, s_over_D: -, Tu: -                      | FilmSpanAverageGroups                                             |
| `FilmSpanAverageGroups::eta`                                | -                                                                                                                                                   | -                                                                 |
| `FilmSpanAverageGroups::weight`                             | -                                                                                                                                                   | -                                                                 |
| `FilmSpanAverageGroups::deta_dM`                            | -                                                                                                                                                   | -                                                                 |
| `FilmSpanAverageGroups::deta_dP`                            | -                                                                                                                                                   | -                                                                 |
| `FilmSpanAverageGroups::d_weight_dM`                        | -                                                                                                                                                   | -                                                                 |
| `FilmSpanAverageGroups::d_weight_dP`                        | -                                                                                                                                                   | -                                                                 |
| `FilmSpanAverageGroups::extrapolated`                       | -                                                                                                                                                   | -                                                                 |
| `FilmSpanAverage::eta`                                      | -                                                                                                                                                   | -                                                                 |
| `FilmSpanAverage::deta_dM`                                  | -                                                                                                                                                   | -                                                                 |
| `FilmSpanAverage::deta_dP`                                  | -                                                                                                                                                   | -                                                                 |
| `FilmSpanAverage::extrapolated`                             | -                                                                                                                                                   | -                                                                 |
| `film_superposition_sellers_and_gradient`                   | eta_rows: -                                                                                                                                         | - (eta), - (d eta/d eta_i)                                        |
| `film_superposition_corrected`                              | eta_rows: -, alpha_between_rows: -                                                                                                                  | - (eta)                                                           |
| `film_superposition_corrected_and_gradient`                 | eta_rows: -, alpha_between_rows: -                                                                                                                  | - (eta), - (d eta/d eta_i)                                        |
| `mills_entry_length_factor`                                 | L_over_D: -                                                                                                                                         | -                                                                 |
| `effusion_approach_nusselt`                                 | Re: -, Pr: -, X_over_L: -                                                                                                                           | -                                                                 |
| `effusion_throat_nusselt`                                   | Re: -, Pr: -, L_over_D: -                                                                                                                           | -                                                                 |
| `effusion_internal_nusselt`                                 | Re: -, Pr: -, X_over_L: -, L_over_D: -                                                                                                              | -                                                                 |
| `flat_plate_turbulent_nusselt_and_derivative`               | Re_x: -, Pr: -                                                                                                                                      | - (Nu_x), - (dNu_x/dRe_x)                                         |
| `flat_plate_turbulent_panel_nusselt`                        | Re_L: -, Pr: -, x_start_over_L: -                                                                                                                   | FlatPlatePanelNusselt                                             |
| `FlatPlatePanelNusselt::Nu`                                 | -                                                                                                                                                   | -                                                                 |
| `FlatPlatePanelNusselt::dNu_dRe_L`                          | -                                                                                                                                                   | -                                                                 |
| `FlatPlatePanelNusselt::Re_x_start`                         | -                                                                                                                                                   | -                                                                 |
| `FlatPlatePanelNusselt::Re_x_end`                           | -                                                                                                                                                   | -                                                                 |
| `FlatPlatePanelNusselt::extrapolated`                       | -                                                                                                                                                   | -                                                                 |
| `flat_plate_turbulent_panel_nusselt_and_derivative`         | Re_L: -, Pr: -, x_start_over_L: -                                                                                                                   | - (Nu_L), - (dNu_L/dRe_L)                                         |
| `mainstream_temperature_correction`                         | mass_flow_ratio: -, a: -, b: -                                                                                                                      | - (alpha)                                                         |
| `equivalent_slot_width`                                     | hole_area: m^2, pitch: m                                                                                                                            | m                                                                 |
| `equivalent_blowing_ratio`                                  | M_baseline: -, area_baseline: m^2, area_equivalent: m^2                                                                                             | - (M)                                                             |
| `film_effectiveness_baldauf_2002`                           | x_over_D: -, M: -, P: -, alpha_deg: deg, s_over_D: -, Tu: -                                                                                         | - (eta)                                                           |
| `film_effectiveness_baldauf_2002_and_derivatives`           | x_over_D: -, M: -, P: -, alpha_deg: deg, s_over_D: -, Tu: -                                                                                         | - (eta), - (deta/dM), - (deta/dP)                                 |
| `film_effectiveness_colban_2011`                            | x_over_D: -, M: -, area_ratio: -, P_over_D: -, t_over_P: -, alpha_deg: deg, density_ratio: -                                                        | - (eta)                                                           |
| `film_effectiveness_colban_2011_and_derivatives`            | x_over_D: -, M: -, area_ratio: -, P_over_D: -, t_over_P: -, alpha_deg: deg, density_ratio: -                                                        | - (eta), - (deta/dM)                                              |
| `slot_film_effectiveness_marek_tacina_1975`                 | x_over_S: -, M: -, Tu: -, C_m_multiplier: -                                                                                                         | - (eta)                                                           |
| `slot_film_effectiveness_marek_tacina_1975_and_derivatives` | x_over_S: -, M: -, Tu: -, C_m_multiplier: -                                                                                                         | - (eta), - (deta/dM)                                              |
| `cooled_wall_heat_flux`                                     | T_hot: K, T_coolant: K, h_hot: W/(m^2*K), h_coolant: W/(m^2*K), eta: -, t_wall: m, k_wall: W/(m*K)                                                  | W/m^2                                                             |
| `heat_flux_from_T_at_edge`                                  | T_measured: K, edge_idx, T_hot: K, T_cold: K, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W                                                           | W/m^2                                                             |
| `heat_flux_from_T_at_depth`                                 | T_measured: K, depth: m, T_hot: K, T_cold: K, h_hot, h_cold: W/(m^2*K), thicknesses: m, conductivities: W/(m*K)                                     | W/m^2                                                             |
| `bulk_T_from_edge_T_and_q`                                  | T_measured: K, edge_idx, q: W/m^2, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W, solve_for: str                                                      | K                                                                 |
| `dT_edge_dT_hot`                                            | edge_idx, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W                                                                                               | - (dT_edge/dT_hot)                                                |
| `dT_edge_dT_cold`                                           | edge_idx, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W                                                                                               | - (dT_edge/dT_cold)                                               |
| `dT_edge_dT_bulk`                                           | edge_idx, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W                                                                                               | - (dT/dT_hot, dT/dT_cold)                                         |
| `dT_edge_dq`                                                | edge_idx, h_hot: W/(m^2*K), t_over_k: m^2*K/W                                                                                                       | K*m^2/W (dT_edge/dq)                                              |

### acoustics.h - Acoustic Mode Analysis

| Function                             | Input Units                                                                                                                                                                                          | Output Unit               |
|--------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------|
| `tube_axial_modes`                   | Tube (L: m, D: m), c: m/s, upstream: BC, downstream: BC, n_max                                                                                                                                       | Hz (vector)               |
| `annulus_axial_modes`                | Annulus (L: m, D_inner: m, D_outer: m), c: m/s, upstream: BC, downstream: BC, n_max                                                                                                                  | Hz (vector)               |
| `annulus_azimuthal_modes`            | Annulus, c: m/s, m_max                                                                                                                                                                               | Hz (vector)               |
| `annulus_modes`                      | Annulus, c: m/s, upstream: BC, downstream: BC, n_max, m_max                                                                                                                                          | Hz (vector)               |
| `modes_in_range`                     | modes, f_min: Hz, f_max: Hz                                                                                                                                                                          | Hz (vector)               |
| `closest_mode`                       | modes, f_target: Hz                                                                                                                                                                                  | AcousticMode              |
| `min_mode_separation`                | modes                                                                                                                                                                                                | Hz                        |
| `axial_mode_upstream`                | f0: Hz, M: -                                                                                                                                                                                         | Hz                        |
| `axial_mode_downstream`              | f0: Hz, M: -                                                                                                                                                                                         | Hz                        |
| `axial_mode_split`                   | f0: Hz, M: -                                                                                                                                                                                         | (Hz, Hz)                  |
| `helmholtz_frequency`                | V: m^3, A_neck: m^2, L_neck: m, c: m/s, end_correction: -                                                                                                                                            | Hz                        |
| `strouhal`                           | f: Hz, L: m, u: m/s                                                                                                                                                                                  | - (St)                    |
| `frequency_from_strouhal`            | St: -, L: m, u: m/s                                                                                                                                                                                  | Hz                        |
| `quarter_wave_frequency`             | L: m, c: m/s                                                                                                                                                                                         | Hz                        |
| `orifice_impedance_with_flow`        | freq: Hz, u_bias: m/s, u_grazing: m/s, d_orifice: m, l_orifice: m, porosity: -, Cd: -, rho: kg/m^3, c: m/s                                                                                           | complex [-]               |
| `quarter_wave_resonator_tmm`         | freq: Hz, L_tube: m, A_duct: m^2, A_tube: m^2, d_orifice: m, l_orifice: m, porosity: -, Cd: -, u_bias: m/s, u_grazing: m/s, rho: kg/m^3, c: m/s                                                      | TransferMatrix            |
| `absorption_from_impedance_norm`     | z_norm: complex [-]                                                                                                                                                                                  | - (alpha)                 |
| `liner_sdof_impedance_norm`          | freq: Hz, orifice: LinerOrificeGeometry, cavity: LinerCavity, flow: LinerFlowState, medium: AcousticMedium                                                                                           | complex [-]               |
| `liner_sdof_absorption`              | freq: Hz, orifice: LinerOrificeGeometry, cavity: LinerCavity, flow: LinerFlowState, medium: AcousticMedium                                                                                           | - (alpha)                 |
| `sweep_liner_sdof_absorption`        | freqs: Hz (vector), orifice: LinerOrificeGeometry, cavity: LinerCavity, flow: LinerFlowState, medium: AcousticMedium                                                                                 | - (alpha vector)          |
| `liner_2dof_serial_impedance_norm`   | freq: Hz, face_orifice: LinerOrificeGeometry, septum_orifice: LinerOrificeGeometry, depth_1: m, depth_2: m, face_flow: LinerFlowState, septum_flow: LinerFlowState, medium: AcousticMedium           | complex [-]               |
| `liner_2dof_serial_absorption`       | freq: Hz, face_orifice: LinerOrificeGeometry, septum_orifice: LinerOrificeGeometry, depth_1: m, depth_2: m, face_flow: LinerFlowState, septum_flow: LinerFlowState, medium: AcousticMedium           | - (alpha)                 |
| `sweep_liner_2dof_serial_absorption` | freqs: Hz (vector), face_orifice: LinerOrificeGeometry, septum_orifice: LinerOrificeGeometry, depth_1: m, depth_2: m, face_flow: LinerFlowState, septum_flow: LinerFlowState, medium: AcousticMedium | - (alpha vector)          |
| `is_whistling_risk`                  | freq: Hz, u_bias: m/s, d_orifice: m                                                                                                                                                                  | bool                      |
| `annular_duct_modes_analytical`      | geom: Annulus, c: m/s, f_max: Hz, bc_ends: BoundaryCondition                                                                                                                                         | AnnularMode vector        |
| `to_acoustic_geometry`               | flow_geom: CanAnnularFlowGeometry, n_cans: -, L_can: m, D_can: m                                                                                                                                     | CanAnnularGeometry struct |
| `half_wave_frequency`                | L: m, c: m/s                                                                                                                                                                                         | Hz                        |
| `stokes_layer`                       | nu: m^2/s, f: Hz                                                                                                                                                                                     | m                         |
| `thermal_layer`                      | alpha: m^2/s, f: Hz                                                                                                                                                                                  | m                         |
| `effective_viscothermal_layer`       | delta_nu: m, delta_kappa: m, gamma: -                                                                                                                                                                | m                         |
| `helmholtz_Q`                        | V: m^3, A_neck: m^2, L_neck: m, nu: m^2/s, alpha: m^2/s, gamma: -, f: Hz                                                                                                                             | - (Q)                     |
| `tube_Q`                             | L: m, D: m, nu: m^2/s, alpha: m^2/s, gamma: -, f: Hz                                                                                                                                                 | - (Q)                     |
| `damping_ratio`                      | Q: -                                                                                                                                                                                                 | - (zeta)                  |
| `bandwidth`                          | f0: Hz, Q: -                                                                                                                                                                                         | Hz                        |

### geometry.h - Residence Time

| Function                          | Input Units                                           | Output Unit |
|-----------------------------------|-------------------------------------------------------|-------------|
| `residence_time`                  | V: m^3, Q: m^3/s                                      | s           |
| `residence_time_tube`             | tube: Tube, Q: m^3/s                                  | s           |
| `residence_time_annulus`          | annulus: Annulus, Q: m^3/s                            | s           |
| `residence_time_can_annular`      | geom: CanAnnularFlowGeometry, Q: m^3/s                | s           |
| `residence_time_mdot`             | V: m^3, mdot: kg/s, rho: kg/m^3                       | s           |
| `residence_time_mdot_can_annular` | geom: CanAnnularFlowGeometry, mdot: kg/s, rho: kg/m^3 | s           |
| `space_velocity`                  | Q: m^3/s, V: m^3                                      | 1/s         |

### state.h - Joint Getters / Setters (C++ & Python tuples)

| Function             | Input Units             | Output Unit             |
|----------------------|-------------------------|-------------------------|
| `State::TPX`         | -                       | tuple(K, Pa, mol/mol)   |
| `State::TPY`         | -                       | tuple(K, Pa, kg/kg)     |
| `State::TP`          | -                       | tuple(K, Pa)            |
| `State::set_TP`      | T: K, P: Pa             | -                       |
| `State::set_TPX`     | T: K, P: Pa, X: mol/mol | -                       |
| `State::set_TPY`     | T: K, P: Pa, Y: kg/kg   | -                       |
| `State::DP_mass`     | -                       | tuple(kg/m^3, Pa)       |
| `State::set_DP_mass` | rho: kg/m^3, P: Pa      | -                       |
| `State::HP_mass`     | -                       | tuple(J/kg, Pa)         |
| `State::set_HP_mass` | h: J/kg, P: Pa          | -                       |
| `State::SP_mass`     | -                       | tuple(J/(kg*K), Pa)     |
| `State::set_SP_mass` | s: J/(kg*K), P: Pa      | -                       |
| `State::UP_mass`     | -                       | tuple(J/kg, Pa)         |
| `State::set_UP_mass` | u: J/kg, Pa             | -                       |
| `State::UV_mass`     | -                       | tuple(J/kg, m^3/kg)     |
| `State::set_UV_mass` | u: J/kg, v: m^3/kg      | -                       |
| `State::SV_mass`     | -                       | tuple(J/(kg*K), m^3/kg) |
| `State::set_SV_mass` | s: J/(kg*K), v: m^3/kg  | -                       |
| `State::PV_mass`     | -                       | tuple(Pa, m^3/kg)       |
| `State::set_PV_mass` | P: Pa, v: m^3/kg        | -                       |
| `State::VH_mass`     | -                       | tuple(m^3/kg, J/kg)     |
| `State::set_VH_mass` | v: m^3/kg, h: J/kg      | -                       |
| `State::SH_mass`     | -                       | tuple(J/(kg*K), J/kg)   |
| `State::set_SH_mass` | s: J/(kg*K), h: J/kg    | -                       |

### heat_transfer.h - ChannelResult fields

| Function                       | Input Units | Output Unit             |
|--------------------------------|-------------|-------------------------|
| `ChannelResult::h`             | -           | W/(m^2*K)               |
| `ChannelResult::Nu`            | -           | - (Nu)                  |
| `ChannelResult::Re`            | -           | - (Re)                  |
| `ChannelResult::Pr`            | -           | - (Pr)                  |
| `ChannelResult::f`             | -           | - (friction/loss coeff) |
| `ChannelResult::dP`            | -           | Pa                      |
| `ChannelResult::M`             | -           | - (M)                   |
| `ChannelResult::T_aw`          | -           | K                       |
| `ChannelResult::q`             | -           | W/m^2                   |
| `ChannelResult::dh_dmdot`      | -           | W/(m^2*K*kg/s)          |
| `ChannelResult::dh_dT`         | -           | W/(m^2*K^2)             |
| `ChannelResult::ddP_dmdot`     | -           | Pa*s/kg                 |
| `ChannelResult::ddP_dvelocity` | -           | Pa*s/m                  |
| `ChannelResult::ddP_dT`        | -           | Pa/K                    |
| `ChannelResult::dT_aw_dmdot`   | -           | K*s/kg                  |
| `ChannelResult::dT_aw_dT`      | -           | - (dimensionless)       |
| `ChannelResult::dT_aw_dP`      | -           | K/Pa                    |
| `ChannelResult::dq_dmdot`      | -           | W*s/(m^2*kg)            |
| `ChannelResult::dq_dT`         | -           | W/(m^2*K)               |
| `ChannelResult::dq_dT_hot`     | -           | W/(m^2*K)               |
| `ChannelResult::dh_dY`         | -           | W/(m^2*K) (per species) |
| `ChannelResult::dT_aw_dY`      | -           | K (per species)         |

### heat_transfer.h - Channel flow functions (HTC + pressure drop)

| Function         | Input Units                                                                                                                                                                             | Output Unit   |
|------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|
| `channel_smooth` | T: K, P: Pa, X: mol/mol, velocity: m/s, diameter: m, length: m, T_hot: K, correlation: str, heating: bool, mu_ratio: -, roughness: m, Nu_multiplier: -, f_multiplier: -, flow_area: m^2 | ChannelResult |

### pin_fin_correlation.h - Pin-fin arrays (#335)

| Function                                | Input Units                                                            | Output Unit          |
|-----------------------------------------|------------------------------------------------------------------------|----------------------|
| `PIN_FIN_RE_EPS`                        | -                                                                      | -                    |
| `PIN_FIN_F_BLEND`                       | -                                                                      | -                    |
| `pin_fin_min_gap_D`                     | geometry: PinFinGeometry                                               | -                    |
| `pin_fin_amin_over_afrontal`            | geometry: PinFinGeometry                                               | -                    |
| `pin_fin_dprime_over_D`                 | geometry: PinFinGeometry                                               | -                    |
| `pin_fin_aprime_over_amin`              | geometry: PinFinGeometry                                               | -                    |
| `pin_fin_dh_over_D`                     | geometry: PinFinGeometry                                               | -                    |
| `pin_fin_area_fractions`                | geometry: PinFinGeometry                                               | PinFinAreaFractions  |
| `pin_fin_array_efficiency`              | h: W/(m^2*K), k_pin: W/(m*K), D: m, H: m, A_f_over_A_t: -              | PinFinEfficiency     |
| `evaluate_pin_fin_nu`                   | correlation_set: PinFinNuSet, geometry: PinFinGeometry, Re_D: -, Pr: - | PinFinNuResult       |
| `evaluate_pin_fin_friction`             | correlation_set: PinFinFrictionSet, geometry: PinFinGeometry, Re_D: -  | PinFinFrictionResult |
| `evaluate_pin_fin_modifier`             | modifier: PinFinRatioModifier, geometry: PinFinGeometry, Re_D: -       | PinFinModifierResult |
| `validate_pin_fin_geometry`             | geometry: PinFinGeometry                                               | None                 |
| `validate_pin_fin_nu_set`               | correlation_set: PinFinNuSet                                           | None                 |
| `validate_pin_fin_friction_set`         | correlation_set: PinFinFrictionSet                                     | None                 |
| `validate_pin_fin_modifier`             | modifier: PinFinRatioModifier                                          | None                 |
| `metzger_1986_staggered_nu`             | -                                                                      | PinFinNuSet          |
| `metzger_1982_staggered_friction`       | -                                                                      | PinFinFrictionSet    |
| `vanfossen_1982_staggered_nu`           | -                                                                      | PinFinNuSet          |
| `damerow_1972_staggered_friction`       | -                                                                      | PinFinFrictionSet    |
| `chyu_1998_nu`                          | arrangement: PinArrangement, surface: PinNuSurface                     | PinFinNuSet          |
| `chyu_1998_inline_over_staggered`       | -                                                                      | PinFinRatioModifier  |
| `chyu_1990_nu`                          | arrangement: PinArrangement, fillet: bool                              | PinFinNuSet          |
| `chyu_1990_friction`                    | arrangement: PinArrangement, fillet: bool                              | PinFinFrictionSet    |
| `chyu_1990_fillet_over_straight`        | arrangement: PinArrangement                                            | PinFinRatioModifier  |
| `PinFinGeometry::S_D`                   | -                                                                      | -                    |
| `PinFinGeometry::X_D`                   | -                                                                      | -                    |
| `PinFinGeometry::H_D`                   | -                                                                      | -                    |
| `PinFinGeometry::N_rows`                | -                                                                      | -                    |
| `PinFinGeometry::arrangement`           | -                                                                      | -                    |
| `PinFinAreaFractions::endwall_exposed`  | -                                                                      | -                    |
| `PinFinAreaFractions::pin`              | -                                                                      | -                    |
| `PinFinAreaFractions::pin_over_total`   | -                                                                      | -                    |
| `PinFinEfficiency::eta_fin`             | -                                                                      | -                    |
| `PinFinEfficiency::eta_t`               | -                                                                      | -                    |
| `PinFinEfficiency::deta_fin_dh`         | -                                                                      | m^2*K/W              |
| `PinFinEfficiency::deta_t_dh`           | -                                                                      | m^2*K/W              |
| `PinFinNuSet::name`                     | -                                                                      | -                    |
| `PinFinNuSet::source`                   | -                                                                      | -                    |
| `PinFinNuSet::validity_source`          | -                                                                      | -                    |
| `PinFinNuSet::provenance`               | -                                                                      | -                    |
| `PinFinNuSet::arrangement`              | -                                                                      | -                    |
| `PinFinNuSet::surface`                  | -                                                                      | -                    |
| `PinFinNuSet::re_basis`                 | -                                                                      | -                    |
| `PinFinNuSet::C`                        | -                                                                      | -                    |
| `PinFinNuSet::Re_exp`                   | -                                                                      | -                    |
| `PinFinNuSet::Pr_exp`                   | -                                                                      | -                    |
| `PinFinNuSet::term_XD`                  | -                                                                      | -                    |
| `PinFinNuSet::term_SD`                  | -                                                                      | -                    |
| `PinFinNuSet::term_HD`                  | -                                                                      | -                    |
| `PinFinNuSet::valid_Re`                 | -                                                                      | -                    |
| `PinFinNuSet::valid_SD`                 | -                                                                      | -                    |
| `PinFinNuSet::valid_XD`                 | -                                                                      | -                    |
| `PinFinNuSet::valid_HD`                 | -                                                                      | -                    |
| `PinFinNuSet::valid_Nrows`              | -                                                                      | -                    |
| `PinFinNuSet::valid_Pr`                 | -                                                                      | -                    |
| `PinFinNuSet::accuracy_Nu`              | -                                                                      | -                    |
| `PinFinFrictionSet::name`               | -                                                                      | -                    |
| `PinFinFrictionSet::source`             | -                                                                      | -                    |
| `PinFinFrictionSet::validity_source`    | -                                                                      | -                    |
| `PinFinFrictionSet::provenance`         | -                                                                      | -                    |
| `PinFinFrictionSet::arrangement`        | -                                                                      | -                    |
| `PinFinFrictionSet::basis`              | -                                                                      | -                    |
| `PinFinFrictionSet::C1`                 | -                                                                      | -                    |
| `PinFinFrictionSet::Re_exp1`            | -                                                                      | -                    |
| `PinFinFrictionSet::C2`                 | -                                                                      | -                    |
| `PinFinFrictionSet::Re_exp2`            | -                                                                      | -                    |
| `PinFinFrictionSet::Re_split`           | -                                                                      | -                    |
| `PinFinFrictionSet::term_SD`            | -                                                                      | -                    |
| `PinFinFrictionSet::term_XD`            | -                                                                      | -                    |
| `PinFinFrictionSet::term_HD`            | -                                                                      | -                    |
| `PinFinFrictionSet::valid_Re`           | -                                                                      | -                    |
| `PinFinFrictionSet::valid_SD`           | -                                                                      | -                    |
| `PinFinFrictionSet::valid_XD`           | -                                                                      | -                    |
| `PinFinFrictionSet::valid_HD`           | -                                                                      | -                    |
| `PinFinFrictionSet::valid_Nrows`        | -                                                                      | -                    |
| `PinFinFrictionSet::accuracy_f`         | -                                                                      | -                    |
| `PinFinRatioModifier::name`             | -                                                                      | -                    |
| `PinFinRatioModifier::source`           | -                                                                      | -                    |
| `PinFinRatioModifier::provenance`       | -                                                                      | -                    |
| `PinFinRatioModifier::from_arrangement` | -                                                                      | -                    |
| `PinFinRatioModifier::to_arrangement`   | -                                                                      | -                    |
| `PinFinRatioModifier::surface`          | -                                                                      | -                    |
| `PinFinRatioModifier::C_Nu`             | -                                                                      | -                    |
| `PinFinRatioModifier::Nu_Re_exp`        | -                                                                      | -                    |
| `PinFinRatioModifier::has_f`            | -                                                                      | -                    |
| `PinFinRatioModifier::C_f`              | -                                                                      | -                    |
| `PinFinRatioModifier::f_Re_exp`         | -                                                                      | -                    |
| `PinFinRatioModifier::valid_Re`         | -                                                                      | -                    |
| `PinFinRatioModifier::valid_SD`         | -                                                                      | -                    |
| `PinFinRatioModifier::valid_XD`         | -                                                                      | -                    |
| `PinFinRatioModifier::valid_HD`         | -                                                                      | -                    |
| `PinFinNuResult::Nu`                    | -                                                                      | -                    |
| `PinFinNuResult::dNu_dRe`               | -                                                                      | -                    |
| `PinFinNuResult::Re_native`             | -                                                                      | -                    |
| `PinFinNuResult::extrapolated`          | -                                                                      | -                    |
| `PinFinFrictionResult::f`               | -                                                                      | -                    |
| `PinFinFrictionResult::df_dRe`          | -                                                                      | -                    |
| `PinFinFrictionResult::extrapolated`    | -                                                                      | -                    |
| `PinFinModifierResult::ratio_Nu`        | -                                                                      | -                    |
| `PinFinModifierResult::dratio_Nu_dRe`   | -                                                                      | -                    |
| `PinFinModifierResult::ratio_f`         | -                                                                      | -                    |
| `PinFinModifierResult::dratio_f_dRe`    | -                                                                      | -                    |
| `PinFinModifierResult::has_f`           | -                                                                      | -                    |
| `PinFinModifierResult::extrapolated`    | -                                                                      | -                    |
| `PinArrangement::Staggered`             | -                                                                      | -                    |
| `PinArrangement::Inline`                | -                                                                      | -                    |
| `PinArrangement::name`                  | -                                                                      | -                    |
| `PinArrangement::value`                 | -                                                                      | -                    |
| `PinNuSurface::Total`                   | -                                                                      | -                    |
| `PinNuSurface::Endwall`                 | -                                                                      | -                    |
| `PinNuSurface::Pin`                     | -                                                                      | -                    |
| `PinNuSurface::name`                    | -                                                                      | -                    |
| `PinNuSurface::value`                   | -                                                                      | -                    |
| `PinReBasis::DiameterVmax`              | -                                                                      | -                    |
| `PinReBasis::VanFossenDprime`           | -                                                                      | -                    |
| `PinReBasis::name`                      | -                                                                      | -                    |
| `PinReBasis::value`                     | -                                                                      | -                    |
| `PinFrictionBasis::PerRowVmax`          | -                                                                      | -                    |
| `PinFrictionBasis::PerRowGapVmax`       | -                                                                      | -                    |
| `PinFrictionBasis::PerRowHalfHeadVmax`  | -                                                                      | -                    |
| `PinFrictionBasis::name`                | -                                                                      | -                    |
| `PinFrictionBasis::value`               | -                                                                      | -                    |

### correlation_status.h - Extrapolation validity utilities

| Function                                        | Input Units                                                                                                                        | Output Unit                           |
|-------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------|
| `mixture_h`                                     | -                                                                                                                                  | -                                     |
| `adiabatic_T_wgs`                               | -                                                                                                                                  | -                                     |
| `normalize_fractions`                           | -                                                                                                                                  | -                                     |
| `convert_to_dry_fractions`                      | -                                                                                                                                  | -                                     |
| `solve_orifice_mdot`                            | dP: Pa, rho: kg/m^3, mu: Pa*s, P_upstream: Pa                                                                                      | kg/s                                  |
| `dry_air`                                       | -                                                                                                                                  | X: mol/mol                            |
| `humid_air_enthalpy_nasa9`                      | T: K, P: Pa, RH: - (0-1)                                                                                                           | J/kg dry air                          |
| `complete_combustion_to_CO2_H2O_with_fraction`  | X: mol/mol, smooth: bool                                                                                                           | (X_b: mol/mol, fuel burn fraction: -) |
| `HumidAir::rh`                                  | -                                                                                                                                  | 0-1                                   |
| `HumidAir::state`                               | -                                                                                                                                  | -                                     |
| `AirProperties::P`                              | -                                                                                                                                  | Pa                                    |
| `AirProperties::Pr`                             | -                                                                                                                                  | -                                     |
| `AirProperties::T`                              | -                                                                                                                                  | K                                     |
| `AirProperties::a`                              | -                                                                                                                                  | m/s                                   |
| `AirProperties::alpha`                          | -                                                                                                                                  | m^2/s                                 |
| `AirProperties::cp`                             | -                                                                                                                                  | J/(kg*K)                              |
| `AirProperties::cv`                             | -                                                                                                                                  | J/(kg*K)                              |
| `AirProperties::gamma`                          | -                                                                                                                                  | -                                     |
| `AirProperties::k`                              | -                                                                                                                                  | W/(m*K)                               |
| `AirProperties::mu`                             | -                                                                                                                                  | Pa*s                                  |
| `AirProperties::nu`                             | -                                                                                                                                  | m^2/s                                 |
| `AirProperties::rho`                            | -                                                                                                                                  | kg/m^3                                |
| `ThermoState::P`                                | -                                                                                                                                  | Pa                                    |
| `ThermoState::T`                                | -                                                                                                                                  | K                                     |
| `ThermoState::a`                                | -                                                                                                                                  | m/s                                   |
| `ThermoState::cp`                               | -                                                                                                                                  | J/(kg*K)                              |
| `ThermoState::cp_mole`                          | -                                                                                                                                  | J/(mol*K)                             |
| `ThermoState::cv`                               | -                                                                                                                                  | J/(kg*K)                              |
| `ThermoState::cv_mole`                          | -                                                                                                                                  | J/(mol*K)                             |
| `ThermoState::gamma`                            | -                                                                                                                                  | -                                     |
| `ThermoState::h`                                | -                                                                                                                                  | J/kg                                  |
| `ThermoState::h_mole`                           | -                                                                                                                                  | J/mol                                 |
| `ThermoState::mw`                               | -                                                                                                                                  | g/mol                                 |
| `ThermoState::rho`                              | -                                                                                                                                  | kg/m^3                                |
| `ThermoState::s`                                | -                                                                                                                                  | J/(kg*K)                              |
| `ThermoState::s_mole`                           | -                                                                                                                                  | J/(mol*K)                             |
| `ThermoState::u`                                | -                                                                                                                                  | J/kg                                  |
| `ThermoState::u_mole`                           | -                                                                                                                                  | J/mol                                 |
| `TransportState::P`                             | -                                                                                                                                  | Pa                                    |
| `TransportState::Pr`                            | -                                                                                                                                  | -                                     |
| `TransportState::T`                             | -                                                                                                                                  | K                                     |
| `TransportState::a`                             | -                                                                                                                                  | m/s                                   |
| `TransportState::alpha`                         | -                                                                                                                                  | m^2/s                                 |
| `TransportState::cp`                            | -                                                                                                                                  | J/(kg*K)                              |
| `TransportState::cv`                            | -                                                                                                                                  | J/(kg*K)                              |
| `TransportState::gamma`                         | -                                                                                                                                  | -                                     |
| `TransportState::k`                             | -                                                                                                                                  | W/(m*K)                               |
| `TransportState::mu`                            | -                                                                                                                                  | Pa*s                                  |
| `TransportState::nu`                            | -                                                                                                                                  | m^2/s                                 |
| `TransportState::rho`                           | -                                                                                                                                  | kg/m^3                                |
| `CompleteState::thermo`                         | -                                                                                                                                  | -                                     |
| `CompleteState::transport`                      | -                                                                                                                                  | -                                     |
| `EquilibriumResult::T_in`                       | -                                                                                                                                  | K                                     |
| `EquilibriumResult::converged`                  | -                                                                                                                                  | -                                     |
| `EquilibriumResult::delta_T`                    | -                                                                                                                                  | K                                     |
| `EquilibriumResult::state`                      | -                                                                                                                                  | -                                     |
| `CombustionState::fuel_burn_fraction`           | -                                                                                                                                  | 0-1                                   |
| `CombustionState::fuel_name`                    | -                                                                                                                                  | -                                     |
| `CombustionState::method`                       | -                                                                                                                                  | -                                     |
| `CombustionState::mixture_fraction`             | -                                                                                                                                  | -                                     |
| `CombustionState::phi`                          | -                                                                                                                                  | -                                     |
| `CombustionState::products`                     | -                                                                                                                                  | -                                     |
| `CombustionState::reactants`                    | -                                                                                                                                  | -                                     |
| `CompressibleFlowSolution::M`                   | -                                                                                                                                  | -                                     |
| `CompressibleFlowSolution::choked`              | -                                                                                                                                  | -                                     |
| `CompressibleFlowSolution::mdot`                | -                                                                                                                                  | kg/s                                  |
| `CompressibleFlowSolution::outlet`              | -                                                                                                                                  | -                                     |
| `CompressibleFlowSolution::stagnation`          | -                                                                                                                                  | -                                     |
| `CompressibleFlowSolution::v`                   | -                                                                                                                                  | m/s                                   |
| `NozzleStation::A`                              | -                                                                                                                                  | m^2                                   |
| `NozzleStation::M`                              | -                                                                                                                                  | -                                     |
| `NozzleStation::P`                              | -                                                                                                                                  | Pa                                    |
| `NozzleStation::T`                              | -                                                                                                                                  | K                                     |
| `NozzleStation::h`                              | -                                                                                                                                  | J/kg                                  |
| `NozzleStation::rho`                            | -                                                                                                                                  | kg/m^3                                |
| `NozzleStation::u`                              | -                                                                                                                                  | m/s                                   |
| `NozzleStation::x`                              | -                                                                                                                                  | m                                     |
| `NozzleSolution::A_throat`                      | -                                                                                                                                  | m^2                                   |
| `NozzleSolution::Pt`                            | -                                                                                                                                  | Pa                                    |
| `NozzleSolution::Tt`                            | -                                                                                                                                  | K                                     |
| `NozzleSolution::choked`                        | -                                                                                                                                  | -                                     |
| `NozzleSolution::h0`                            | -                                                                                                                                  | J/kg                                  |
| `NozzleSolution::inlet`                         | -                                                                                                                                  | -                                     |
| `NozzleSolution::mdot`                          | -                                                                                                                                  | kg/s                                  |
| `NozzleSolution::outlet`                        | -                                                                                                                                  | -                                     |
| `NozzleSolution::profile`                       | -                                                                                                                                  | -                                     |
| `NozzleSolution::x_throat`                      | -                                                                                                                                  | m                                     |
| `FannoStation::M`                               | -                                                                                                                                  | -                                     |
| `FannoStation::P`                               | -                                                                                                                                  | Pa                                    |
| `FannoStation::Re`                              | -                                                                                                                                  | -                                     |
| `FannoStation::T`                               | -                                                                                                                                  | K                                     |
| `FannoStation::f`                               | -                                                                                                                                  | -                                     |
| `FannoStation::h`                               | -                                                                                                                                  | J/kg                                  |
| `FannoStation::rho`                             | -                                                                                                                                  | kg/m^3                                |
| `FannoStation::s`                               | -                                                                                                                                  | J/(kg*K)                              |
| `FannoStation::u`                               | -                                                                                                                                  | m/s                                   |
| `FannoStation::x`                               | -                                                                                                                                  | m                                     |
| `FannoSolution::D`                              | -                                                                                                                                  | m                                     |
| `FannoSolution::L`                              | -                                                                                                                                  | m                                     |
| `FannoSolution::L_choke`                        | -                                                                                                                                  | m                                     |
| `FannoSolution::Re_in`                          | -                                                                                                                                  | -                                     |
| `FannoSolution::choked`                         | -                                                                                                                                  | -                                     |
| `FannoSolution::f`                              | -                                                                                                                                  | -                                     |
| `FannoSolution::f_avg`                          | -                                                                                                                                  | -                                     |
| `FannoSolution::h0`                             | -                                                                                                                                  | J/kg                                  |
| `FannoSolution::inlet`                          | -                                                                                                                                  | -                                     |
| `FannoSolution::mdot`                           | -                                                                                                                                  | kg/s                                  |
| `FannoSolution::outlet`                         | -                                                                                                                                  | -                                     |
| `FannoSolution::profile`                        | -                                                                                                                                  | -                                     |
| `FannoMachState::M`                             | -                                                                                                                                  | -                                     |
| `FannoMachState::T`                             | -                                                                                                                                  | K                                     |
| `FannoMachState::P`                             | -                                                                                                                                  | Pa                                    |
| `FannoMachState::rho`                           | -                                                                                                                                  | kg/m^3                                |
| `FannoMachState::u`                             | -                                                                                                                                  | m/s                                   |
| `FannoMachState::a`                             | -                                                                                                                                  | m/s                                   |
| `FannoDuctResult::choked`                       | -                                                                                                                                  | -                                     |
| `FannoDuctResult::G`                            | -                                                                                                                                  | kg/(m^2*s)                            |
| `FannoDuctResult::L_star`                       | -                                                                                                                                  | m                                     |
| `FannoDuctResult::inlet`                        | -                                                                                                                                  | FannoMachState                        |
| `FannoDuctResult::exit`                         | -                                                                                                                                  | FannoMachState                        |
| `FannoDuctResult::Pt_exit`                      | -                                                                                                                                  | Pa                                    |
| `FANNO_MACH_REL_TOL`                            | -                                                                                                                                  | -                                     |
| `FannoChannelFlow::G`                           | -                                                                                                                                  | kg/(m^2*s)                            |
| `FannoChannelFlow::choked`                      | -                                                                                                                                  | -                                     |
| `FannoChannelFlow::M_in`                        | -                                                                                                                                  | -                                     |
| `FannoChannelFlow::M_exit`                      | -                                                                                                                                  | -                                     |
| `FannoChannelFlow::M_exit_root`                 | -                                                                                                                                  | -                                     |
| `FannoChannelFlow::dG_dPt0`                     | -                                                                                                                                  | kg/(m^2*s*Pa)                         |
| `FannoChannelFlow::dG_dTt0`                     | -                                                                                                                                  | kg/(m^2*s*K)                          |
| `FannoChannelFlow::dG_dP_target`                | -                                                                                                                                  | kg/(m^2*s*Pa)                         |
| `FannoChannelFlow::dG_dTt_exit`                 | -                                                                                                                                  | kg/(m^2*s*K)                          |
| `FannoChannelFlow::dG_ddir`                     | -                                                                                                                                  | kg/(m^2*s)                            |
| `FannoRayleighMarch::M_in`                      | -                                                                                                                                  | -                                     |
| `FannoRayleighMarch::M_exit`                    | -                                                                                                                                  | -                                     |
| `FannoRayleighMarch::zeta_min`                  | -                                                                                                                                  | -                                     |
| `FannoRayleighMarch::zeta_peak`                 | -                                                                                                                                  | -                                     |
| `FannoRayleighMarch::P_exit`                    | -                                                                                                                                  | Pa                                    |
| `FannoRayleighMarch::Pt_exit`                   | -                                                                                                                                  | Pa                                    |
| `FannoRayleighMarch::T_exit`                    | -                                                                                                                                  | K                                     |
| `FANNO_FLOW_DRIVE_FLOOR`                        | -                                                                                                                                  | Pa                                    |
| `ThrustResult::P_exit`                          | -                                                                                                                                  | Pa                                    |
| `ThrustResult::mdot`                            | -                                                                                                                                  | kg/s                                  |
| `ThrustResult::specific_impulse`                | -                                                                                                                                  | s                                     |
| `ThrustResult::thrust`                          | -                                                                                                                                  | N                                     |
| `ThrustResult::thrust_coefficient`              | -                                                                                                                                  | -                                     |
| `ThrustResult::u_exit`                          | -                                                                                                                                  | m/s                                   |
| `nozzle_thrust_cd`                              | Tt: K, Pt: Pa, P_design: Pa, P_amb: Pa, A_inlet: m^2, A_throat: m^2, A_exit: m^2, x_throat: m, x_exit: m, X: mol/mol               | ThrustResult                          |
| `TransferMatrix::T11`                           | -                                                                                                                                  | -                                     |
| `TransferMatrix::T12`                           | -                                                                                                                                  | -                                     |
| `TransferMatrix::T21`                           | -                                                                                                                                  | -                                     |
| `TransferMatrix::T22`                           | -                                                                                                                                  | -                                     |
| `CanAnnularGeometry::area_can`                  | -                                                                                                                                  | m^2                                   |
| `CanAnnularGeometry::area_plenum`               | -                                                                                                                                  | m^2                                   |
| `CanAnnularGeometry::length_can`                | -                                                                                                                                  | m                                     |
| `CanAnnularGeometry::n_cans`                    | -                                                                                                                                  | -                                     |
| `CanAnnularGeometry::radius_plenum`             | -                                                                                                                                  | m                                     |
| `BlochMode::frequency`                          | -                                                                                                                                  | Hz                                    |
| `BlochMode::m_azimuthal`                        | -                                                                                                                                  | -                                     |
| `BlochMode::n_cans`                             | -                                                                                                                                  | -                                     |
| `BlochMode::symmetry_type`                      | -                                                                                                                                  | -                                     |
| `can_annular_eigenmodes`                        | geom: CanAnnularGeometry, c_can: m/s, c_plenum: m/s, rho_can: kg/m^3, rho_plenum: kg/m^3, f_max: Hz, bc_can_top: BoundaryCondition | BlochMode vector (frequency: Hz)      |
| `AnnularDuctGeometry::D_inner`                  | -                                                                                                                                  | m                                     |
| `AnnularDuctGeometry::D_mean`                   | -                                                                                                                                  | m                                     |
| `AnnularDuctGeometry::D_outer`                  | -                                                                                                                                  | m                                     |
| `AnnularDuctGeometry::L`                        | -                                                                                                                                  | m                                     |
| `AnnularDuctGeometry::area`                     | -                                                                                                                                  | m^2                                   |
| `AnnularDuctGeometry::circumference`            | -                                                                                                                                  | m                                     |
| `AnnularDuctGeometry::gap`                      | -                                                                                                                                  | m                                     |
| `AnnularDuctGeometry::length`                   | -                                                                                                                                  | m                                     |
| `AnnularDuctGeometry::n_azimuthal_max`          | -                                                                                                                                  | -                                     |
| `AnnularDuctGeometry::radius_inner`             | -                                                                                                                                  | m                                     |
| `AnnularDuctGeometry::radius_outer`             | -                                                                                                                                  | m                                     |
| `AnnularDuctGeometry::volume`                   | -                                                                                                                                  | m^3                                   |
| `AnnularMode::frequency`                        | -                                                                                                                                  | Hz                                    |
| `AnnularMode::m_azimuthal`                      | -                                                                                                                                  | -                                     |
| `AnnularMode::mode_type`                        | -                                                                                                                                  | -                                     |
| `AnnularMode::n_axial`                          | -                                                                                                                                  | -                                     |
| `annular_duct_eigenmodes`                       | geom: Annulus, c: m/s, rho: kg/m^3, f_max: Hz, bc_ends: BoundaryCondition                                                          | AnnularMode vector (frequency: Hz)    |
| `AcousticMedium::c`                             | -                                                                                                                                  | -                                     |
| `AcousticMedium::rho`                           | -                                                                                                                                  | -                                     |
| `LinerOrificeGeometry::Cd`                      | -                                                                                                                                  | -                                     |
| `LinerOrificeGeometry::d_orifice`               | -                                                                                                                                  | -                                     |
| `LinerOrificeGeometry::l_orifice`               | -                                                                                                                                  | -                                     |
| `LinerOrificeGeometry::porosity`                | -                                                                                                                                  | -                                     |
| `LinerFlowState::u_bias`                        | -                                                                                                                                  | -                                     |
| `LinerFlowState::u_grazing`                     | -                                                                                                                                  | -                                     |
| `LinerCavity::depth`                            | -                                                                                                                                  | -                                     |
| `overall_htc_wall_multilayer`                   | -                                                                                                                                  | -                                     |
| `Tube::D`                                       | -                                                                                                                                  | m                                     |
| `Tube::L`                                       | -                                                                                                                                  | m                                     |
| `Tube::area`                                    | -                                                                                                                                  | m^2                                   |
| `Tube::perimeter`                               | -                                                                                                                                  | m                                     |
| `Tube::volume`                                  | -                                                                                                                                  | m^3                                   |
| `Annulus::D_inner`                              | -                                                                                                                                  | m                                     |
| `Annulus::D_mean`                               | -                                                                                                                                  | m                                     |
| `Annulus::D_outer`                              | -                                                                                                                                  | m                                     |
| `Annulus::L`                                    | -                                                                                                                                  | m                                     |
| `Annulus::area`                                 | -                                                                                                                                  | m^2                                   |
| `Annulus::circumference`                        | -                                                                                                                                  | m                                     |
| `Annulus::gap`                                  | -                                                                                                                                  | m                                     |
| `Annulus::length`                               | -                                                                                                                                  | m                                     |
| `Annulus::n_azimuthal_max`                      | -                                                                                                                                  | -                                     |
| `Annulus::radius_inner`                         | -                                                                                                                                  | m                                     |
| `Annulus::radius_outer`                         | -                                                                                                                                  | m                                     |
| `Annulus::volume`                               | -                                                                                                                                  | m^3                                   |
| `CanAnnularFlowGeometry::D_inner`               | -                                                                                                                                  | m                                     |
| `CanAnnularFlowGeometry::D_mean`                | -                                                                                                                                  | m                                     |
| `CanAnnularFlowGeometry::D_outer`               | -                                                                                                                                  | m                                     |
| `CanAnnularFlowGeometry::D_primary`             | -                                                                                                                                  | m                                     |
| `CanAnnularFlowGeometry::L`                     | -                                                                                                                                  | m                                     |
| `CanAnnularFlowGeometry::L_primary`             | -                                                                                                                                  | m                                     |
| `CanAnnularFlowGeometry::L_transition`          | -                                                                                                                                  | m                                     |
| `CanAnnularFlowGeometry::area`                  | -                                                                                                                                  | m^2                                   |
| `CanAnnularFlowGeometry::area_primary`          | -                                                                                                                                  | m^2                                   |
| `CanAnnularFlowGeometry::circumference`         | -                                                                                                                                  | m                                     |
| `CanAnnularFlowGeometry::gap`                   | -                                                                                                                                  | m                                     |
| `CanAnnularFlowGeometry::length_total`          | -                                                                                                                                  | m                                     |
| `CanAnnularFlowGeometry::volume`                | -                                                                                                                                  | m^3                                   |
| `CanAnnularFlowGeometry::volume_primary`        | -                                                                                                                                  | m^3                                   |
| `CanAnnularFlowGeometry::volume_total`          | -                                                                                                                                  | m^3                                   |
| `CanAnnularFlowGeometry::volume_transition`     | -                                                                                                                                  | m^3                                   |
| `AcousticMode::frequency`                       | -                                                                                                                                  | Hz                                    |
| `AcousticMode::label`                           | -                                                                                                                                  | -                                     |
| `AcousticMode::n_axial`                         | -                                                                                                                                  | -                                     |
| `AcousticMode::n_azimuthal`                     | -                                                                                                                                  | -                                     |
| `AcousticProperties::frequency`                 | -                                                                                                                                  | Hz                                    |
| `AcousticProperties::impedance`                 | -                                                                                                                                  | Pa*s/m                                |
| `AcousticProperties::particle_velocity`         | -                                                                                                                                  | m/s                                   |
| `AcousticProperties::spl`                       | -                                                                                                                                  | dB                                    |
| `AcousticProperties::wavelength`                | -                                                                                                                                  | m                                     |
| `channel_dP`                                    | -                                                                                                                                  | Pa                                    |
| `channel_dP_mdot`                               | -                                                                                                                                  | Pa                                    |
| `channel_velocity`                              | -                                                                                                                                  | m/s                                   |
| `channel_mdot`                                  | -                                                                                                                                  | kg/s                                  |
| `channel_area`                                  | -                                                                                                                                  | m^2                                   |
| `channel_volume`                                | -                                                                                                                                  | m^3                                   |
| `channel_roughness`                             | -                                                                                                                                  | m                                     |
| `pressure_drop_channel`                         | -                                                                                                                                  | Pa                                    |
| `nusselt_channel`                               | -                                                                                                                                  | -                                     |
| `htc_channel`                                   | -                                                                                                                                  | W/(m^2*K)                             |
| `dynamic_pressure`                              | -                                                                                                                                  | Pa                                    |
| `velocity_from_q`                               | -                                                                                                                                  | m/s                                   |
| `OrificeGeometry::D`                            | -                                                                                                                                  | m                                     |
| `OrificeGeometry::area`                         | -                                                                                                                                  | m^2                                   |
| `OrificeGeometry::beta`                         | -                                                                                                                                  | -                                     |
| `OrificeGeometry::bevel`                        | -                                                                                                                                  | rad                                   |
| `OrificeGeometry::d`                            | -                                                                                                                                  | m                                     |
| `OrificeGeometry::is_valid`                     | -                                                                                                                                  | -                                     |
| `OrificeGeometry::r`                            | -                                                                                                                                  | m                                     |
| `OrificeGeometry::r_over_d`                     | -                                                                                                                                  | -                                     |
| `OrificeGeometry::t`                            | -                                                                                                                                  | m                                     |
| `OrificeGeometry::t_over_d`                     | -                                                                                                                                  | -                                     |
| `OrificeGeometry::tap`                          | -                                                                                                                                  | -                                     |
| `OrificeTap::Corner`                            | -                                                                                                                                  | -                                     |
| `OrificeTap::D_D2`                              | -                                                                                                                                  | -                                     |
| `OrificeTap::Flange`                            | -                                                                                                                                  | -                                     |
| `OrificeTap::name`                              | -                                                                                                                                  | -                                     |
| `OrificeTap::value`                             | -                                                                                                                                  | -                                     |
| `OrificeState::Re_D`                            | -                                                                                                                                  | -                                     |
| `OrificeState::Re_d`                            | -                                                                                                                                  | -                                     |
| `OrificeState::dP`                              | -                                                                                                                                  | Pa                                    |
| `OrificeState::mu`                              | -                                                                                                                                  | Pa*s                                  |
| `OrificeState::rho`                             | -                                                                                                                                  | kg/m^3                                |
| `MeteringCdCorrelation::Constant`               | -                                                                                                                                  | -                                     |
| `MeteringCdCorrelation::ReaderHarrisGallagher`  | -                                                                                                                                  | -                                     |
| `MeteringCdCorrelation::name`                   | -                                                                                                                                  | -                                     |
| `MeteringCdCorrelation::value`                  | -                                                                                                                                  | -                                     |
| `MeteringCdCorrelation::UserFunction`           | -                                                                                                                                  | -                                     |
| `DischargeCdCorrelation::McGreehanSchotsch1988` | -                                                                                                                                  | -                                     |
| `DischargeCdCorrelation::Idelchik1966Sharp`     | -                                                                                                                                  | -                                     |
| `DischargeCdCorrelation::Idelchik1966Thick`     | -                                                                                                                                  | -                                     |
| `DischargeCdCorrelation::Idelchik1966Beveled`   | -                                                                                                                                  | -                                     |
| `DischargeCdCorrelation::Idelchik1966Rounded`   | -                                                                                                                                  | -                                     |
| `DischargeCdCorrelation::Lichtarowicz1965`      | -                                                                                                                                  | -                                     |
| `DischargeCdCorrelation::Constant`              | -                                                                                                                                  | -                                     |
| `DischargeCdCorrelation::name`                  | -                                                                                                                                  | -                                     |
| `DischargeCdCorrelation::value`                 | -                                                                                                                                  | -                                     |
| `DischargeHoleGeometry::d`                      | m                                                                                                                                  | -                                     |
| `DischargeHoleGeometry::L`                      | m                                                                                                                                  | -                                     |
| `DischargeHoleGeometry::r`                      | m                                                                                                                                  | -                                     |
| `DischargeHoleGeometry::bevel`                  | m                                                                                                                                  | -                                     |
| `DischargeHoleGeometry::inclination_deg`        | deg                                                                                                                                | -                                     |
| `DischargeHoleGeometry::L_over_d`               | -                                                                                                                                  | - (L/d)                               |
| `DischargeHoleGeometry::r_over_d`               | -                                                                                                                                  | - (r/d)                               |
| `DischargeHoleGeometry::bevel_over_d`           | -                                                                                                                                  | - (l/d)                               |
| `DischargeHoleGeometry::area`                   | -                                                                                                                                  | m^2                                   |
| `DischargeHoleGeometry::is_valid`               | -                                                                                                                                  | -                                     |
| `DischargeHoleGeometry::is_inclined`            | -                                                                                                                                  | - (bool)                              |
| `DischargeHoleState::Re`                        | -                                                                                                                                  | -                                     |
| `DischargeHoleState::U1_over_Vi`                | -                                                                                                                                  | -                                     |
| `discharge_cd`                                  | correlation: -, hole: m, flow: -, Cd_constant: -                                                                                   | - (Cd)                                |
| `discharge_cd_in_range`                         | correlation: -, hole: m, flow: -                                                                                                   | - (bool)                              |
| `discharge_cd_re_floor`                         | correlation: -                                                                                                                     | - (Re)                                |
| `discharge_inclination_dzeta`                   | inclination_deg: deg                                                                                                               | - (dzeta)                             |
| `discharge_pipe_zeta`                           | correlation: -, hole: m, beta: -                                                                                                   | -                                     |
| `discharge_pipe_cd_ratio`                       | correlation: -, hole: m, beta: -                                                                                                   | -                                     |
| `discharge_cd_and_derivatives`                  | correlation: -, hole: m, flow: -, Cd_constant: -                                                                                   | - (Cd), - (dCd/dRe), - (dCd/d(U1/Vi)) |
| `orifice_mdot_Cd`                               | -                                                                                                                                  | kg/s                                  |
| `orifice_dP_Cd`                                 | -                                                                                                                                  | Pa                                    |
| `orifice_Cd_from_measurement`                   | -                                                                                                                                  | kg/s                                  |
| `orifice_K_from_Cd`                             | -                                                                                                                                  | -                                     |
| `orifice_Cd_from_K`                             | -                                                                                                                                  | -                                     |
| `OrificeFlowResult::Cd`                         | -                                                                                                                                  | -                                     |
| `OrificeFlowResult::Re_D`                       | -                                                                                                                                  | -                                     |
| `OrificeFlowResult::Re_d`                       | -                                                                                                                                  | -                                     |
| `OrificeFlowResult::epsilon`                    | -                                                                                                                                  | -                                     |
| `OrificeFlowResult::mdot`                       | -                                                                                                                                  | kg/s                                  |
| `OrificeFlowResult::rho_corrected`              | -                                                                                                                                  | kg/m^3                                |
| `OrificeFlowResult::v`                          | -                                                                                                                                  | m/s                                   |
| `SonicStatic::P`                                | -                                                                                                                                  | Pa                                    |
| `SonicStatic::T`                                | -                                                                                                                                  | K                                     |
| `SonicStatic::dP_dTt`                           | -                                                                                                                                  | Pa/K                                  |
| `SonicStatic::dP_dG`                            | -                                                                                                                                  | Pa/(kg/(m^2*s))                       |
| `SonicStatic::dP_dY`                            | -                                                                                                                                  | Pa                                    |
| `ApproachStatic::P`                             | -                                                                                                                                  | Pa                                    |
| `ApproachStatic::dP_dPt`                        | -                                                                                                                                  | -                                     |
| `ApproachStatic::dP_dTt`                        | -                                                                                                                                  | Pa/K                                  |
| `ApproachStatic::dP_dG`                         | -                                                                                                                                  | Pa/(kg/(m^2*s))                       |
| `ApproachStatic::dP_dY`                         | -                                                                                                                                  | Pa                                    |
| `ApproachStatic::choked`                        | -                                                                                                                                  | - (bool)                              |
| `OrificeTapPressure::p_tap`                     | -                                                                                                                                  | Pa                                    |
| `OrificeTapPressure::loss_ratio`                | -                                                                                                                                  | -                                     |
| `OrificeTapPressure::d_dP_static_up`            | -                                                                                                                                  | -                                     |
| `OrificeTapPressure::d_dP_total_up`             | -                                                                                                                                  | -                                     |
| `OrificeTapPressure::d_dP_total_down`           | -                                                                                                                                  | -                                     |
| `OrificeTapPressure::d_dC`                      | -                                                                                                                                  | Pa                                    |
| `OrificeTapPressure::floored`                   | -                                                                                                                                  | - (bool)                              |
| `OrificeTapPressure::in_range`                  | -                                                                                                                                  | - (bool)                              |
| `orifice_flow_state`                            | geom: OrificeGeometry, dP: Pa, T: K, P: Pa, X: mol/mol, Z: -, kappa: -, correlation: -                                             | OrificeFlowResult struct              |
| `FlowSolution::Cd`                              | -                                                                                                                                  | -                                     |
| `FlowSolution::L_choke`                         | -                                                                                                                                  | m                                     |
| `FlowSolution::M`                               | -                                                                                                                                  | -                                     |
| `FlowSolution::P_out`                           | -                                                                                                                                  | Pa                                    |
| `FlowSolution::Re`                              | -                                                                                                                                  | -                                     |
| `FlowSolution::T_out`                           | -                                                                                                                                  | K                                     |
| `FlowSolution::choked`                          | -                                                                                                                                  | -                                     |
| `FlowSolution::dP`                              | -                                                                                                                                  | Pa                                    |
| `FlowSolution::f`                               | -                                                                                                                                  | -                                     |
| `FlowSolution::h0`                              | -                                                                                                                                  | J/kg                                  |
| `FlowSolution::rho`                             | -                                                                                                                                  | kg/m^3                                |
| `UnitInfo::input`                               | -                                                                                                                                  | -                                     |
| `UnitInfo::output`                              | -                                                                                                                                  | -                                     |
| `is_well_behaved`                               | v: any, lo: any, hi: any                                                                                                           | bool                                  |
| `set_warning_handler`                           | handler: callable(str)                                                                                                             | -                                     |
| `get_warning_handler`                           | -                                                                                                                                  | callable(str)                         |
| `suppress_warnings`                             | -                                                                                                                                  | context manager                       |

---

## Dimensionless Quantities

| Quantity                | Symbol | Unit     | Definition                    |
|-------------------------|--------|----------|-------------------------------|
| Mach number             | M      | -        | v / a                         |
| Reynolds number         | Re     | -        | rho*V*L / mu                  |
| Prandtl number          | Pr     | -        | mu*Cp / k                     |
| Peclet number           | Pe     | -        | V*L / alpha                   |
| Isentropic exponent     | gamma  | -        | Cp / Cv                       |
| Equivalence ratio       | phi    | -        | (F/A) / (F/A)_stoich          |
| Mixture fraction        | Z      | -        | Bilger definition             |
| Discharge coefficient   | Cd     | -        | mdot_actual / mdot_ideal      |
| Friction factor (Darcy) | f      | -        | dP / (L/D * rho*v^2/2)        |
| Pressure ratio          | -      | -        | P / P0                        |
| Diameter ratio          | beta   | -        | d / D                         |

---

## Common Reference Values

| Quantity                | Value          | Unit   | Notes                    |
|-------------------------|----------------|--------|--------------------------|
| Standard pressure       | 101325         | Pa     | 1 atm                    |
| Standard temperature    | 298.15         | K      | 25 C                     |
| Standard gravity        | 9.80665        | m/s^2  | Used for Isp             |
| Sea level air density   | ~1.225         | kg/m^3 | At 15 C, 101325 Pa       |

---

## Summary: Key Unit Conventions

1. **Temperature**: Always Kelvin (K)
2. **Pressure**: Always Pascal (Pa)
3. **Molar mass**: g/mol (historical convention; convert to kg/mol for mass-basis calculations)
4. **Thermodynamic properties**: Molar basis (J/mol, J/(mol*K)) in thermo functions
5. **Compressible flow**: Mass basis (J/kg, J/(kg*K)) in flow solutions
6. **Fractions**: mole fractions X (mol/mol), mass fractions Y (kg/kg)
7. **Relative humidity**: Fraction (0-1), not percentage
8. **Angles**: Radians (rad)
