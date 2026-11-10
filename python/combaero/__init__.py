"""Public Python API for the combaero package.

This module re-exports selected functions from the compiled _core extension.
It is written to work both when running from the source tree and when
combaero is installed as a wheel.
"""

from importlib import import_module as _import_module

try:  # Python 3.8+
    from importlib.metadata import PackageNotFoundError as _PackageNotFoundError
    from importlib.metadata import version as _pkg_version
except ModuleNotFoundError:  # pragma: no cover - very old Python
    from importlib_metadata import PackageNotFoundError as _PackageNotFoundError
    from importlib_metadata import version as _pkg_version  # type: ignore[import-not-found]


def _load_version() -> str:
    """Return the installed distribution version, or a safe fallback.

    When running from a source tree without an installed wheel, the
    distribution metadata may not be available; in that case we expose a
    placeholder string instead of raising.
    """

    try:
        return _pkg_version("combaero")
    except _PackageNotFoundError:
        return "0.0.0+local"


__version__: str = _load_version()


try:
    # Preferred: local extension in the same package (installed wheel or
    # in-tree build where _core was successfully built next to this file).
    from ._core import (
        AccuracyProvenance,
        RibCorrelationSet,
        StatedAccuracy,
        RibGeometry,
        GShapeModel,
        RAlphaShape,
        RatioBaseline,
        RatioBelowFloor,
        RibProvenance,
        RibRatioOptions,
        RibRatioResult,
        RibRatioSet,
        RibShape,
        RibRange,
        RibResult,
        RibTerm,
        evaluate_rib,
        evaluate_rib_ratio,
        han_1988_orthogonal,
        han_1989_narrow_channel,
        han_zhang_lee_1991,
        han_park_1988_angled,
        rallabandi_2009_high_re,
        validate_rib_ratio_options,
        validate_rib_ratio_set,
        validate_rib_set,
        taslim_spring_1987,
        PIN_FIN_RE_EPS,
        PIN_FIN_F_BLEND,
        PinArrangement,
        PinNuSurface,
        PinReBasis,
        PinFrictionBasis,
        PinFinGeometry,
        PinFinAreaFractions,
        PinFinEfficiency,
        PinFinNuSet,
        PinFinFrictionSet,
        PinFinRatioModifier,
        PinFinNuResult,
        PinFinFrictionResult,
        PinFinModifierResult,
        pin_fin_min_gap_D,
        pin_fin_amin_over_afrontal,
        pin_fin_dprime_over_D,
        pin_fin_aprime_over_amin,
        pin_fin_dh_over_D,
        pin_fin_area_fractions,
        pin_fin_array_efficiency,
        evaluate_pin_fin_nu,
        evaluate_pin_fin_friction,
        evaluate_pin_fin_modifier,
        validate_pin_fin_geometry,
        validate_pin_fin_nu_set,
        validate_pin_fin_friction_set,
        validate_pin_fin_modifier,
        metzger_1986_staggered_nu,
        metzger_1982_staggered_friction,
        vanfossen_1982_staggered_nu,
        damerow_1972_staggered_friction,
        chyu_1998_nu,
        chyu_1998_inline_over_staggered,
        chyu_1990_nu,
        chyu_1990_friction,
        chyu_1990_fillet_over_straight,
        FLORSCHUETZ_1981_DEFAULT_CD,
        FLORSCHUETZ_SU_1985_ETA_A,
        FLORSCHUETZ_SU_1985_ETA_N,
        FLORSCHUETZ_SU_1985_ETA_R,
        CombustorLossResult,
        combustor_fraction_loss,
        combustor_head_loss,
        ImpingementEtaR,
        ImpingementRange,
        ImpingementThermalBC,
        JetArrayCorrelationSet,
        JetArrayGeometricFit,
        JetArrayImpingementResult,
        JetHolePattern,
        JetRowGeometry,
        JetRowHeatTransfer,
        StationHalfDrop,
        CrossflowRatio,
        JetImpulse,
        MergeFaceState,
        SingleJetImpingementResult,
        SingleJetImpingementSet,
        crossflow_to_jet_ratio_at_row,
        crossflow_to_jet_ratio_at_x,
        florschuetz_1981_inline,
        florschuetz_1981_staggered,
        goldstein_1986_single_jet,
        impingement_eta_r,
        jet_array_impingement_nu,
        jet_row_gc_gj_per_flow_ratio,
        jet_row_heat_transfer,
        station_half_drop,
        crossflow_velocity_ratio,
        FannoMachState,
        FannoDuctResult,
        FANNO_MACH_REL_TOL,
        fanno_state_at_mach,
        fanno_dx_dmach,
        fanno_length_between,
        fanno_inlet_mach,
        fanno_sonic_mass_flux,
        fanno_duct,
        fanno_choked_mass_flux,
        TransportDerivatives,
        Omega22Result,
        transport_and_dT,
        omega22_and_derivative,
        dcp_R_dT,
        set_thermo_model,
        thermo_model,
        FannoChannelFlow,
        FANNO_FLOW_DRIVE_FLOOR,
        fanno_channel_flow,
        FannoRayleighMarch,
        fanno_rayleigh_march,
        fanno_rayleigh_choked_mass_flux,
        fanno_rayleigh_channel_flow,
        CROSSFLOW_VI_FLOOR,
        channel_entry_drop,
        segment_area_change_drop,
        ROGERS_HERSH_1976_CAP_BAND,
        FLORSCHUETZ_METZGER_SU_1982_CAP_BAND,
        ROGERS_HERSH_1976_INJECTION_A,
        ROGERS_HERSH_1976_INJECTION_B,
        ROGERS_HERSH_1976_SUCTION_A,
        ROGERS_HERSH_1976_SUCTION_B,
        CrossflowCd,
        rogers_hersh_1976_cd,
        florschuetz_metzger_su_1982_has_plate,
        florschuetz_metzger_su_1982_cd,
        Rowbury2001Set,
        ROWBURY_2001_K_PERP_MIN,
        ROWBURY_2001_UNIQUE_MARGIN,
        Rowbury2001Coefficients,
        Rowbury2001Cd,
        Rowbury2001KPerp,
        rowbury_2001_coefficients,
        rowbury_2001_k_floor,
        rowbury_2001_in_range,
        rowbury_2001_exit_crossflow_cd,
        rowbury_2001_k_perp,
        rowbury_2001_k_fit_max,
        rowbury_2001_k_data_max,
        ROWBURY_2001_FADE_WIDTH,
        Rowbury2001ExitDensity,
        rowbury_2001_exit_density,
        STATION_KAPPA_MERGE_NORMAL,
        STATION_KAPPA_BLEED_BASSETT,
        WALL_HTC_KNEE,
        STEFAN_BOLTZMANN,
        REEVES_EMISSIVITY_COEFF,
        LUMINOSITY_HYDROGEN_COEFF,
        BEAM_LENGTH_COEFF_LEFEBVRE,
        BEAM_LENGTH_COEFF_HOWELL,
        GasEmissivityResult,
        LinerGasRadiationResult,
        gas_emissivity_luminous,
        luminosity_factor_hydrogen,
        liner_gas_radiation,
        liner_casing_radiation,
        liner_casing_exchange_factor,
        mean_beam_length,
        elemental_fuel_air_ratio,
        SEBAN1960_C_NEAR,
        SEBAN1960_C_NEAR_PRINTED,
        SEBAN1960_C_FAR,
        SEBAN1960_C_FAR_PRINTED,
        SLOT_FILM_HTC_HANDBACK_BAND,
        SLOT_FILM_BRANCH_NEAR,
        SLOT_FILM_BRANCH_BLEND,
        SLOT_FILM_BRANCH_FAR,
        SLOT_FILM_BRANCH_WALL_JET,
        SLOT_FILM_BRANCH_NONE,
        SlotFilmStanton,
        SlotFilmNusselt,
        SlotFilmHtc,
        slot_film_stanton_seban_1960,
        slot_film_nusselt_ballal_lefebvre_1973,
        slot_film_htc_seban_1960,
        slot_film_htc_ballal_lefebvre_1973,
        MIXER_HEAT_MDOT_FLOOR,
        mixer_heat_per_unit_flow,
        SENSIBLE_ENTHALPY_REF_T,
        chamber_merge_face_state,
        duct_static_state,
        duct_htc_static_state,
        duct_htc_static_state_dY,
        DUCT_HTC_SONIC_MARGIN,
        jet_impulse,
        single_jet_impingement,
        single_jet_impingement_nu,
        validate_jet_array_set,
        validate_single_jet_set,
        AcousticMedium,
        AcousticMode,
        AcousticProperties,
        AirProperties,
        AnnularDuctGeometry,
        AnnularMode,
        BlochMode,
        BoundaryCondition,
        CanAnnularGeometry,
        Cd_from_zeta,
        Cd_sharp_thin_plate,
        Cd_sharp_thin_plate_and_dRe,
        MeteringCdCorrelation,
        OrificeTap,
        metering_cd_and_dRe,
        metering_cd_in_range,
        metering_cd_re_min,
        permanent_loss_ratio_and_dC,
        velocity_of_approach_factor,
        jobson_force_defect,
        jobson_ci_from_cd,
        jobson_contraction_coefficient,
        jobson_expansibility,
        JOBSON_1955_BETA_MAX,
        DischargeCdCorrelation,
        DischargeHoleGeometry,
        DischargeHoleState,
        CorrelationValidity,
        CorrelationResult,
        CombustionMethod,
        CombustionState,
        CompleteState,
        CompressibleFlowSolution,
        EquilibriumResult,
        FannoSolution,
        FannoStation,
        IncompressibleFlowSolution,
        LinerCavity,
        LinerFlowState,
        LinerOrificeGeometry,
        NozzleSolution,
        NozzleStation,
        OrificeFlowResult,
        OrificeGeometry,
        OrificeState,
        OrificeResult,
        OrificeTapPressure,
        ChannelSolverResult,
        MixtureState,
        MomentumChamberResult,
        AreaChangeResult,
        AreaChangeLoss,
        AreaChangeFlow,
        AreaChangeElementFlow,
        TurbomachinePath,
        BordaCarnotStep,
        State,
        Stream,
        Registry,
        ThermoState,
        ThrustResult,
        TransferMatrix,
        TransportState,
        UnitInfo,
        absorption_from_impedance_norm,
        acoustic_impedance,
        acoustic_properties,
        adiabatic_T_wgs,
        adiabatic_wall_temperature,
        air_properties,
        all_units,
        annular_area,
        annular_duct_eigenmodes,
        annular_duct_modes_analytical,
        annulus_axial_modes,
        annulus_azimuthal_modes,
        annulus_modes,
        axial_mode_downstream,
        axial_mode_split,
        axial_mode_upstream,
        bandwidth,
        bernoulli_P2,
        bernoulli_v2,
        beta_from_diameters,
        bilger_stoich_mixture_fraction_mass,
        bilger_Z_from_equivalence_ratio_mass,
        bulk_T_from_edge_T_and_q,
        calc_T_from_cp,
        calc_T_from_cp_roots,
        CpRoot,
        calc_T_from_h,
        calc_T_from_h_mass,
        calc_T_from_s,
        calc_T_from_s_mass,
        calc_T_from_u,
        calc_T_from_u_mass,
        calc_T_from_sv_mass,
        calc_T_from_sh_mass,
        dh_dT,
        ds_dT,
        dcp_dT,
        can_annular_eigenmodes,
        capacity_ratio,
        channel_smooth,
        ChannelResult,
        WallCouplingResult,
        closest_mode,
        combustion_equilibrium,
        combustion_state,
        combustion_state_from_streams,
        common_name,
        complete_combustion,
        complete_combustion_isothermal,
        complete_combustion_to_CO2_H2O,
        complete_combustion_to_CO2_H2O_with_fraction,
        complete_state,
        convert_to_dry_fractions,
        cooled_wall_heat_flux,
        equivalent_blowing_ratio,
        effusion_approach_nusselt,
        effusion_internal_nusselt,
        effusion_throat_nusselt,
        flat_plate_turbulent_nusselt_and_derivative,
        flat_plate_turbulent_panel_nusselt,
        flat_plate_turbulent_panel_nusselt_and_derivative,
        FlatPlatePanelNusselt,
        mills_entry_length_factor,
        equivalent_slot_width,
        film_superposition_corrected,
        film_superposition_corrected_and_gradient,
        film_superposition_sellers,
        effusion_panel_film_effectiveness,
        EffusionPanelFilm,
        film_span_average_baldauf_2002,
        film_span_average_baldauf_2002_groups,
        film_span_average_colban_2011,
        slot_film_span_average_marek_tacina_1975,
        FilmSpanAverage,
        FilmSpanAverageGroups,
        film_superposition_sellers_and_gradient,
        mainstream_temperature_correction,
        film_effectiveness_baldauf_2002,
        film_effectiveness_baldauf_2002_and_derivatives,
        film_effectiveness_colban_2011,
        film_effectiveness_colban_2011_and_derivatives,
        slot_film_effectiveness_marek_tacina_1975,
        slot_film_effectiveness_marek_tacina_1975_and_derivatives,
        cp,
        cp_mass,
        cp_secant_mass_and_jacobian,
        cp_secant_mass_dY,
        critical_pressure_ratio,
        cv,
        cv_mass,
        damping_ratio,
        channel_flow,
        channel_flow_rough,
        fanno_channel,
        fanno_channel_rough,
        density,
        dewpoint,
        dry_air,
        dT_edge_dq,
        dT_edge_dT_bulk,
        dT_edge_dT_cold,
        dT_edge_dT_hot,
        dynamic_pressure,
        effective_viscothermal_layer,
        effectiveness_counterflow,
        effectiveness_parallelflow,
        equivalence_ratio,
        equivalence_ratio_from_bilger_Z_mass,
        equivalence_ratio_mass,
        equivalence_ratio_mole,
        expansibility_factor,
        fanno_max_length,
        # Nozzle flow is already channel-like
        formula,
        formula_to_name,
        frequency_from_strouhal,
        frequency_from_wavelength,
        friction_colebrook,
        friction_haaland,
        friction_petukhov,
        friction_serghides,
        fuel_lhv_mass,
        fuel_lhv_molar,
        get_units,
        h,
        h_mass,
        half_wave_frequency,
        has_units,
        heat_flux,
        heat_flux_from_T_at_depth,
        heat_flux_from_T_at_edge,
        heat_rate,
        heat_rate_from_effectiveness,
        heat_transfer_area,
        heat_transfer_dT,
        helmholtz_frequency,
        helmholtz_Q,
        htc_from_nusselt,
        humid_air_composition,
        humid_air_density,
        humid_air_enthalpy,
        humid_air_enthalpy_nasa9,
        humidity_ratio,
        HumidAir,
        hydraulic_diameter,
        hydraulic_diameter_annulus,
        hydraulic_diameter_rect,
        input_units,
        is_whistling_risk,
        isentropic_expansion_coefficient,
        k_bondcoat_nicocraly,
        k_hastelloy_x,
        k_haynes188,
        k_haynes230,
        k_haynes263,
        k_haynes25,
        k_inconel625,
        k_nimonic75,
        k_copper,
        k_alumina,
        k_x22crmov12_1,
        k_inconel718,
        k_stainless_steel_304,
        k_stainless_steel_310,
        k_stainless_steel_316,
        k_tbc_ysz,
        kinematic_viscosity,
        kinetic_energy,
        liner_2dof_serial_absorption,
        liner_2dof_serial_impedance_norm,
        liner_sdof_absorption,
        liner_sdof_impedance_norm,
        list_functions_with_units,
        list_materials,
        get_material_conductivity,
        get_material_conductivity_and_slope,
        get_material_info,
        material_in_range,
        GENERIC_METAL_SOLIDUS,
        lmtd,
        lmtd_counterflow,
        lmtd_parallelflow,
        mach_from_pressure_ratio,
        mach_number,
        mass_flux_isentropic,
        MassStream,
        EntropyGenerationResult,
        entropy_generation,
        heat_entropy_flux,
        mass_to_mole,
        discharge_cd,
        discharge_cd_and_derivatives,
        discharge_cd_in_range,
        discharge_cd_re_floor,
        discharge_inclination_dzeta,
        discharge_pipe_cd_ratio,
        discharge_pipe_zeta,
        mcgreehan_schotsch_1988_cd,
        mcgreehan_schotsch_1988_cd_and_derivatives,
        mcgreehan_schotsch_1988_crossflow_cd,
        mcgreehan_schotsch_1988_expansion_factor,
        mcgreehan_schotsch_1988_expansion_nozzle,
        mcgreehan_schotsch_1988_expansion_orifice,
        min_mode_separation,
        mix,
        mixture_h,
        modes_in_range,
        molar_volume,
        mole_to_mass,
        mwmix,
        name_to_formula,
        normalize_fractions,
        nozzle_cd,
        nozzle_flow,
        nozzle_thrust,
        nozzle_thrust_cd,
        ntu,
        num_species,
        species_name,
        species_index_from_name,
        species_molar_mass,
        species_molar_mass_from_name,
        nusselt_dittus_boelter,
        nusselt_gnielinski,
        nusselt_petukhov,
        nusselt_sieder_tate,
        orifice_area,
        orifice_area_from_beta,
        orifice_Cd_from_K,
        orifice_Cd_from_measurement,
        orifice_dP,
        orifice_dP_Cd,
        orifice_flow,
        orifice_flow_state,
        orifice_flow_thermo,
        orifice_impedance_with_flow,
        orifice_K_from_Cd,
        orifice_mdot,
        orifice_mdot_Cd,
        orifice_Q,
        orifice_Re_d_from_mdot,
        orifice_velocity,
        orifice_velocity_from_mdot,
        output_units,
        overall_htc,
        overall_htc_wall,
        overall_htc_wall_multilayer,
        oxygen_required_per_kg_fuel,
        oxygen_required_per_kg_mixture,
        oxygen_required_per_mol_fuel,
        oxygen_required_per_mol_mixture,
        particle_velocity,
        Pt_from_static,
        P_from_stagnation,
        peclet,
        prandtl,
        pressure_loss,
        quarter_wave_frequency,
        quarter_wave_resonator_tmm,
        recovery_factor,
        reforming_equilibrium,
        reforming_equilibrium_adiabatic,
        relative_humidity_from_dewpoint,
        residence_time,
        residence_time_annulus,
        residence_time_can_annular,
        residence_time_mdot,
        residence_time_mdot_can_annular,
        residence_time_tube,
        reynolds,
        thermal_performance_factor,
        s,
        s_mass,
        saturation_vapor_pressure,
        get_warning_handler,
        is_well_behaved,
        set_equivalence_ratio_mass,
        set_warning_handler,
        set_equivalence_ratio_mole,
        set_fuel_stream_for_CO2,
        set_fuel_stream_for_CO2_dry,
        set_fuel_stream_for_O2,
        set_fuel_stream_for_O2_dry,
        set_fuel_stream_for_phi,
        set_fuel_stream_for_Tad,
        set_oxidizer_stream_for_CO2,
        set_oxidizer_stream_for_CO2_dry,
        set_oxidizer_stream_for_O2,
        set_oxidizer_stream_for_O2_dry,
        set_oxidizer_stream_for_Tad,
        solve_A_eff_from_mdot,
        solve_orifice_mdot,
        solve_Pt_from_mdot,
        solve_P_back_from_mdot,
        sound_pressure_level,
        space_velocity,
        specific_gas_constant,
        speed_of_sound,
        stagnation_from_static,
        bulk_velocity,
        stokes_layer,
        strouhal,
        sweep_liner_2dof_serial_absorption,
        sweep_liner_sdof_absorption,
        Tt_from_static,
        Tt_from_static_v,
        T_adiabatic_wall,
        T_adiabatic_wall_mach,
        T_from_stagnation,
        thermal_conductivity,
        thermal_diffusivity,
        thermal_layer,
        thermal_resistance,
        thermal_resistance_wall,
        thermo_state,
        to_acoustic_geometry,
        transport_state,
        tube_axial_modes,
        tube_Q,
        u,
        u_mass,
        vapor_pressure,
        velocity_from_pressure_loss,
        velocity_from_q,
        viscosity,
        wall_temperature_profile,
        water_vapor_mole_fraction,
        wavelength,
        wet_bulb_temperature,
        wgs_equilibrium,
        wgs_equilibrium_adiabatic,
        zeta_from_Cd,
        geometry,
        h0_from_static,
        v_from_h0,
        sharp_area_change,
        conical_area_change,
        sudden_expansion_zeta,
        sudden_contraction_zeta,
        conical_diffuser_zeta,
        conical_nozzle_zeta,
        compressible_borda_carnot,
        area_change_loss,
        area_change_flow,
        area_change_element_flow,
        compression,
        expansion,
        complete_combustion_products_and_jacobian,
        combustion_equilibrium_and_jacobian,
        CompleteCombustionJacobian,
        EquilibriumJacobian,
        # Channel/circular functions
        htc_circular_channel,
        nusselt_circular_channel,
        channel_dP,
        channel_dP_mdot,
        channel_velocity,
        channel_mdot,
        channel_pressure_drop,
        pressure_drop_channel,
        circular_area,
        cylinder_volume,
        channel_roughness,
        standard_channel_roughness,
        # Vatistas n-vortex model
        vatistas_v0_bar,
        vatistas_dv0_bar_drbar,
        vatistas_vr_bar,
        vatistas_dvr_bar_drbar,
        vatistas_pressure_integral,
        vatistas_d_pressure_integral_drbar,
        vatistas_v_theta,
        vatistas_v_theta_and_jacobians,
        vatistas_delta_p,
        vatistas_delta_p_and_jacobians,
    )

    molar_mass = mwmix
    Tube = geometry.Tube
    Annulus = geometry.Annulus
    CanAnnularFlowGeometry = geometry.CanAnnularFlowGeometry
except (ModuleNotFoundError, ImportError) as e:
    print(f"DEBUG: Import failed with {type(e).__name__}: {e}")
    # Fallback: attempt to import from an installed combaero package that
    # already has _core available, then re-export the symbols.
    _core = _import_module("combaero._core")
    mixture_h = _core.mixture_h
    adiabatic_T_wgs = _core.adiabatic_T_wgs
    cp = _core.cp
    cp_mass = _core.cp_mass
    h = _core.h
    h_mass = _core.h_mass
    s = _core.s
    s_mass = _core.s_mass
    cv = _core.cv
    cv_mass = _core.cv_mass
    density = _core.density
    specific_gas_constant = _core.specific_gas_constant
    isentropic_expansion_coefficient = _core.isentropic_expansion_coefficient
    speed_of_sound = _core.speed_of_sound
    molar_mass = _core.mwmix
    molar_volume = _core.molar_volume
    viscosity = _core.viscosity
    thermal_conductivity = _core.thermal_conductivity
    prandtl = _core.prandtl
    thermal_diffusivity = _core.thermal_diffusivity
    reynolds = _core.reynolds
    peclet = _core.peclet
    mole_to_mass = _core.mole_to_mass
    mass_to_mole = _core.mass_to_mole
    normalize_fractions = _core.normalize_fractions
    convert_to_dry_fractions = _core.convert_to_dry_fractions
    equivalence_ratio_mole = _core.equivalence_ratio_mole
    equivalence_ratio = _core.equivalence_ratio
    equivalence_ratio_mole = _core.equivalence_ratio_mole
    equivalence_ratio_mass = _core.equivalence_ratio_mass
    set_equivalence_ratio_mass = _core.set_equivalence_ratio_mass
    bilger_stoich_mixture_fraction_mass = _core.bilger_stoich_mixture_fraction_mass
    bilger_Z_from_equivalence_ratio_mass = _core.bilger_Z_from_equivalence_ratio_mass
    equivalence_ratio_from_bilger_Z_mass = _core.equivalence_ratio_from_bilger_Z_mass
    expansibility_factor = _core.expansibility_factor
    solve_orifice_mdot = _core.solve_orifice_mdot
    dry_air = _core.dry_air
    humid_air_composition = _core.humid_air_composition
    humid_air_density = _core.humid_air_density
    humid_air_enthalpy = _core.humid_air_enthalpy
    humid_air_enthalpy_nasa9 = _core.humid_air_enthalpy_nasa9
    humidity_ratio = _core.humidity_ratio
    HumidAir = _core.HumidAir
    saturation_vapor_pressure = _core.saturation_vapor_pressure
    vapor_pressure = _core.vapor_pressure
    water_vapor_mole_fraction = _core.water_vapor_mole_fraction
    wet_bulb_temperature = _core.wet_bulb_temperature
    dewpoint = _core.dewpoint
    relative_humidity_from_dewpoint = _core.relative_humidity_from_dewpoint
    formula_to_name = _core.formula_to_name
    name_to_formula = _core.name_to_formula
    common_name = _core.common_name
    formula = _core.formula
    calc_T_from_h = _core.calc_T_from_h
    calc_T_from_h_mass = _core.calc_T_from_h_mass
    calc_T_from_s = _core.calc_T_from_s
    calc_T_from_s_mass = _core.calc_T_from_s_mass
    calc_T_from_cp = _core.calc_T_from_cp
    calc_T_from_cp_roots = _core.calc_T_from_cp_roots
    CpRoot = _core.CpRoot
    calc_T_from_u = _core.calc_T_from_u
    calc_T_from_u_mass = _core.calc_T_from_u_mass
    calc_T_from_sv_mass = _core.calc_T_from_sv_mass
    calc_T_from_sh_mass = _core.calc_T_from_sh_mass
    dh_dT = _core.dh_dT
    ds_dT = _core.ds_dT
    dcp_dT = _core.dcp_dT
    num_species = _core.num_species
    species_name = _core.species_name
    species_index_from_name = _core.species_index_from_name
    species_molar_mass = _core.species_molar_mass
    species_molar_mass_from_name = _core.species_molar_mass_from_name
    oxygen_required_per_mol_fuel = _core.oxygen_required_per_mol_fuel
    oxygen_required_per_kg_fuel = _core.oxygen_required_per_kg_fuel
    oxygen_required_per_mol_mixture = _core.oxygen_required_per_mol_mixture
    oxygen_required_per_kg_mixture = _core.oxygen_required_per_kg_mixture
    fuel_lhv_molar = _core.fuel_lhv_molar
    fuel_lhv_mass = _core.fuel_lhv_mass
    complete_combustion_to_CO2_H2O = _core.complete_combustion_to_CO2_H2O
    complete_combustion_to_CO2_H2O_with_fraction = (
        _core.complete_combustion_to_CO2_H2O_with_fraction
    )
    # State-based API
    MixtureState = _core.MixtureState
    State = _core.State
    Stream = _core.Stream
    Registry = _core.Registry
    mix = _core.mix
    # Inverse solvers - find fuel stream
    set_fuel_stream_for_Tad = _core.set_fuel_stream_for_Tad
    set_fuel_stream_for_O2 = _core.set_fuel_stream_for_O2
    set_fuel_stream_for_O2_dry = _core.set_fuel_stream_for_O2_dry
    set_fuel_stream_for_CO2 = _core.set_fuel_stream_for_CO2
    set_fuel_stream_for_CO2_dry = _core.set_fuel_stream_for_CO2_dry
    # Inverse solvers - find oxidizer stream
    set_oxidizer_stream_for_Tad = _core.set_oxidizer_stream_for_Tad
    set_oxidizer_stream_for_O2 = _core.set_oxidizer_stream_for_O2
    set_oxidizer_stream_for_O2_dry = _core.set_oxidizer_stream_for_O2_dry
    set_oxidizer_stream_for_CO2 = _core.set_oxidizer_stream_for_CO2
    set_oxidizer_stream_for_CO2_dry = _core.set_oxidizer_stream_for_CO2_dry
    CombustionMethod = _core.CombustionMethod
    CombustionState = _core.CombustionState
    complete_combustion = _core.complete_combustion
    complete_combustion_isothermal = _core.complete_combustion_isothermal
    wgs_equilibrium = _core.wgs_equilibrium
    wgs_equilibrium_adiabatic = _core.wgs_equilibrium_adiabatic
    reforming_equilibrium = _core.reforming_equilibrium
    reforming_equilibrium_adiabatic = _core.reforming_equilibrium_adiabatic
    combustion_equilibrium = _core.combustion_equilibrium
    sharp_area_change = _core.sharp_area_change
    conical_area_change = _core.conical_area_change
    AreaChangeResult = _core.AreaChangeResult
    AreaChangeLoss = _core.AreaChangeLoss
    AreaChangeFlow = _core.AreaChangeFlow
    BordaCarnotStep = _core.BordaCarnotStep
    sudden_expansion_zeta = _core.sudden_expansion_zeta
    sudden_contraction_zeta = _core.sudden_contraction_zeta
    conical_diffuser_zeta = _core.conical_diffuser_zeta
    conical_nozzle_zeta = _core.conical_nozzle_zeta
    compressible_borda_carnot = _core.compressible_borda_carnot
    area_change_loss = _core.area_change_loss
    area_change_flow = _core.area_change_flow
    area_change_element_flow = _core.area_change_element_flow
    compression = _core.compression
    expansion = _core.expansion
    complete_combustion_products_and_jacobian = _core.complete_combustion_products_and_jacobian
    combustion_equilibrium_and_jacobian = _core.combustion_equilibrium_and_jacobian
    CompleteCombustionJacobian = _core.CompleteCombustionJacobian
    EquilibriumJacobian = _core.EquilibriumJacobian
    AreaChangeElementFlow = _core.AreaChangeElementFlow
    TurbomachinePath = _core.TurbomachinePath
    # Compressible flow
    CompressibleFlowSolution = _core.CompressibleFlowSolution
    NozzleStation = _core.NozzleStation
    NozzleSolution = _core.NozzleSolution
    FannoStation = _core.FannoStation
    FannoSolution = _core.FannoSolution
    IncompressibleFlowSolution = _core.IncompressibleFlowSolution
    nozzle_flow = _core.nozzle_flow
    solve_A_eff_from_mdot = _core.solve_A_eff_from_mdot
    solve_P_back_from_mdot = _core.solve_P_back_from_mdot
    solve_Pt_from_mdot = _core.solve_Pt_from_mdot
    critical_pressure_ratio = _core.critical_pressure_ratio
    mach_from_pressure_ratio = _core.mach_from_pressure_ratio
    mach_number = _core.mach_number
    mass_flux_isentropic = _core.mass_flux_isentropic
    nozzle_cd = _core.nozzle_cd
    fanno_max_length = _core.fanno_max_length
    # Thrust
    ThrustResult = _core.ThrustResult
    nozzle_thrust = _core.nozzle_thrust
    nozzle_thrust_cd = _core.nozzle_thrust_cd
    # Friction
    friction_haaland = _core.friction_haaland
    friction_serghides = _core.friction_serghides
    friction_colebrook = _core.friction_colebrook
    friction_petukhov = _core.friction_petukhov
    # Vatistas n-vortex model
    vatistas_v0_bar = _core.vatistas_v0_bar
    vatistas_dv0_bar_drbar = _core.vatistas_dv0_bar_drbar
    vatistas_vr_bar = _core.vatistas_vr_bar
    vatistas_dvr_bar_drbar = _core.vatistas_dvr_bar_drbar
    vatistas_pressure_integral = _core.vatistas_pressure_integral
    vatistas_d_pressure_integral_drbar = _core.vatistas_d_pressure_integral_drbar
    vatistas_v_theta = _core.vatistas_v_theta
    vatistas_v_theta_and_jacobians = _core.vatistas_v_theta_and_jacobians
    vatistas_delta_p = _core.vatistas_delta_p
    vatistas_delta_p_and_jacobians = _core.vatistas_delta_p_and_jacobians
    # Heat transfer
    nusselt_dittus_boelter = _core.nusselt_dittus_boelter
    nusselt_gnielinski = _core.nusselt_gnielinski
    nusselt_sieder_tate = _core.nusselt_sieder_tate
    nusselt_petukhov = _core.nusselt_petukhov
    htc_from_nusselt = _core.htc_from_nusselt
    overall_htc = _core.overall_htc
    overall_htc_wall = _core.overall_htc_wall
    overall_htc_wall_multilayer = _core.overall_htc_wall_multilayer
    thermal_resistance = _core.thermal_resistance
    thermal_resistance_wall = _core.thermal_resistance_wall
    lmtd = _core.lmtd
    lmtd_counterflow = _core.lmtd_counterflow
    lmtd_parallelflow = _core.lmtd_parallelflow
    heat_rate = _core.heat_rate
    heat_flux = _core.heat_flux
    heat_transfer_area = _core.heat_transfer_area
    heat_transfer_dT = _core.heat_transfer_dT
    wall_temperature_profile = _core.wall_temperature_profile
    STEFAN_BOLTZMANN = _core.STEFAN_BOLTZMANN
    REEVES_EMISSIVITY_COEFF = _core.REEVES_EMISSIVITY_COEFF
    LUMINOSITY_HYDROGEN_COEFF = _core.LUMINOSITY_HYDROGEN_COEFF
    BEAM_LENGTH_COEFF_LEFEBVRE = _core.BEAM_LENGTH_COEFF_LEFEBVRE
    BEAM_LENGTH_COEFF_HOWELL = _core.BEAM_LENGTH_COEFF_HOWELL
    GasEmissivityResult = _core.GasEmissivityResult
    LinerGasRadiationResult = _core.LinerGasRadiationResult
    gas_emissivity_luminous = _core.gas_emissivity_luminous
    luminosity_factor_hydrogen = _core.luminosity_factor_hydrogen
    liner_gas_radiation = _core.liner_gas_radiation
    liner_casing_radiation = _core.liner_casing_radiation
    liner_casing_exchange_factor = _core.liner_casing_exchange_factor
    mean_beam_length = _core.mean_beam_length
    elemental_fuel_air_ratio = _core.elemental_fuel_air_ratio
    SEBAN1960_C_NEAR = _core.SEBAN1960_C_NEAR
    SEBAN1960_C_NEAR_PRINTED = _core.SEBAN1960_C_NEAR_PRINTED
    SEBAN1960_C_FAR = _core.SEBAN1960_C_FAR
    SEBAN1960_C_FAR_PRINTED = _core.SEBAN1960_C_FAR_PRINTED
    SLOT_FILM_HTC_HANDBACK_BAND = _core.SLOT_FILM_HTC_HANDBACK_BAND
    SLOT_FILM_BRANCH_NEAR = _core.SLOT_FILM_BRANCH_NEAR
    SLOT_FILM_BRANCH_BLEND = _core.SLOT_FILM_BRANCH_BLEND
    SLOT_FILM_BRANCH_FAR = _core.SLOT_FILM_BRANCH_FAR
    SLOT_FILM_BRANCH_WALL_JET = _core.SLOT_FILM_BRANCH_WALL_JET
    SLOT_FILM_BRANCH_NONE = _core.SLOT_FILM_BRANCH_NONE
    SlotFilmStanton = _core.SlotFilmStanton
    SlotFilmNusselt = _core.SlotFilmNusselt
    SlotFilmHtc = _core.SlotFilmHtc
    slot_film_stanton_seban_1960 = _core.slot_film_stanton_seban_1960
    slot_film_nusselt_ballal_lefebvre_1973 = _core.slot_film_nusselt_ballal_lefebvre_1973
    slot_film_htc_seban_1960 = _core.slot_film_htc_seban_1960
    slot_film_htc_ballal_lefebvre_1973 = _core.slot_film_htc_ballal_lefebvre_1973

    ntu = _core.ntu
    capacity_ratio = _core.capacity_ratio
    effectiveness_counterflow = _core.effectiveness_counterflow
    effectiveness_parallelflow = _core.effectiveness_parallelflow
    heat_rate_from_effectiveness = _core.heat_rate_from_effectiveness
    heat_flux_from_T_at_edge = _core.heat_flux_from_T_at_edge
    heat_flux_from_T_at_depth = _core.heat_flux_from_T_at_depth
    bulk_T_from_edge_T_and_q = _core.bulk_T_from_edge_T_and_q
    dT_edge_dT_hot = _core.dT_edge_dT_hot
    dT_edge_dT_cold = _core.dT_edge_dT_cold
    dT_edge_dT_bulk = _core.dT_edge_dT_bulk
    dT_edge_dq = _core.dT_edge_dq
    # Acoustics
    Tube = _core.geometry.Tube
    Annulus = _core.geometry.Annulus
    CanAnnularFlowGeometry = _core.geometry.CanAnnularFlowGeometry
    BoundaryCondition = _core.BoundaryCondition
    AcousticMode = _core.AcousticMode
    to_acoustic_geometry = _core.to_acoustic_geometry
    tube_axial_modes = _core.tube_axial_modes
    annulus_axial_modes = _core.annulus_axial_modes
    annulus_azimuthal_modes = _core.annulus_azimuthal_modes
    annulus_modes = _core.annulus_modes
    modes_in_range = _core.modes_in_range
    closest_mode = _core.closest_mode
    min_mode_separation = _core.min_mode_separation
    # Mean flow correction
    axial_mode_upstream = _core.axial_mode_upstream
    axial_mode_downstream = _core.axial_mode_downstream
    axial_mode_split = _core.axial_mode_split
    # Helmholtz resonator
    helmholtz_frequency = _core.helmholtz_frequency
    # Strouhal
    strouhal = _core.strouhal
    frequency_from_strouhal = _core.frequency_from_strouhal
    # Convenience
    quarter_wave_frequency = _core.quarter_wave_frequency
    half_wave_frequency = _core.half_wave_frequency
    # Viscothermal
    stokes_layer = _core.stokes_layer
    thermal_layer = _core.thermal_layer
    effective_viscothermal_layer = _core.effective_viscothermal_layer
    # Quality factor
    helmholtz_Q = _core.helmholtz_Q
    tube_Q = _core.tube_Q
    damping_ratio = _core.damping_ratio
    bandwidth = _core.bandwidth
    # Acoustic properties
    AcousticProperties = _core.AcousticProperties
    acoustic_properties = _core.acoustic_properties
    wavelength = _core.wavelength
    frequency_from_wavelength = _core.frequency_from_wavelength
    acoustic_impedance = _core.acoustic_impedance
    sound_pressure_level = _core.sound_pressure_level
    particle_velocity = _core.particle_velocity
    # Residence time
    residence_time = _core.residence_time
    residence_time_tube = _core.residence_time_tube
    residence_time_annulus = _core.residence_time_annulus
    residence_time_can_annular = _core.residence_time_can_annular
    residence_time_mdot = _core.residence_time_mdot
    residence_time_mdot_can_annular = _core.residence_time_mdot_can_annular
    space_velocity = _core.space_velocity
    # Incompressible flow
    bernoulli_P2 = _core.bernoulli_P2
    bernoulli_v2 = _core.bernoulli_v2
    orifice_mdot = _core.orifice_mdot
    orifice_Q = _core.orifice_Q
    orifice_velocity = _core.orifice_velocity
    orifice_area = _core.orifice_area
    orifice_dP = _core.orifice_dP
    channel_dP = _core.channel_dP
    channel_dP_mdot = _core.channel_dP_mdot
    channel_velocity = _core.channel_velocity
    channel_mdot = _core.channel_mdot
    channel_pressure_drop = _core.channel_pressure_drop
    pressure_drop_channel = _core.channel_pressure_drop
    dynamic_pressure = _core.dynamic_pressure
    velocity_from_q = _core.velocity_from_q
    hydraulic_diameter = _core.hydraulic_diameter
    hydraulic_diameter_rect = _core.hydraulic_diameter_rect
    hydraulic_diameter_annulus = _core.hydraulic_diameter_annulus
    channel_area = _core.circular_area
    channel_volume = _core.cylinder_volume
    channel_roughness = _core.channel_roughness
    standard_channel_roughness = _core.standard_channel_roughness
    annular_area = _core.annular_area
    # Orifice Cd correlations
    OrificeGeometry = _core.OrificeGeometry
    OrificeState = _core.OrificeState
    MeteringCdCorrelation = _core.MeteringCdCorrelation
    OrificeTap = _core.OrificeTap
    DischargeCdCorrelation = _core.DischargeCdCorrelation
    DischargeHoleGeometry = _core.DischargeHoleGeometry
    DischargeHoleState = _core.DischargeHoleState
    CorrelationValidity = _core.CorrelationValidity
    CorrelationResult = _core.CorrelationResult
    Cd_sharp_thin_plate = _core.Cd_sharp_thin_plate
    Cd_sharp_thin_plate_and_dRe = _core.Cd_sharp_thin_plate_and_dRe
    metering_cd_and_dRe = _core.metering_cd_and_dRe
    metering_cd_in_range = _core.metering_cd_in_range
    metering_cd_re_min = _core.metering_cd_re_min
    permanent_loss_ratio_and_dC = _core.permanent_loss_ratio_and_dC
    velocity_of_approach_factor = _core.velocity_of_approach_factor
    jobson_force_defect = _core.jobson_force_defect
    jobson_ci_from_cd = _core.jobson_ci_from_cd
    jobson_contraction_coefficient = _core.jobson_contraction_coefficient
    jobson_expansibility = _core.jobson_expansibility
    JOBSON_1955_BETA_MAX = _core.JOBSON_1955_BETA_MAX
    discharge_cd = _core.discharge_cd
    discharge_cd_and_derivatives = _core.discharge_cd_and_derivatives
    discharge_cd_in_range = _core.discharge_cd_in_range
    discharge_cd_re_floor = _core.discharge_cd_re_floor
    discharge_inclination_dzeta = _core.discharge_inclination_dzeta
    discharge_pipe_cd_ratio = _core.discharge_pipe_cd_ratio
    discharge_pipe_zeta = _core.discharge_pipe_zeta
    mcgreehan_schotsch_1988_cd = _core.mcgreehan_schotsch_1988_cd
    mcgreehan_schotsch_1988_cd_and_derivatives = _core.mcgreehan_schotsch_1988_cd_and_derivatives
    mcgreehan_schotsch_1988_crossflow_cd = _core.mcgreehan_schotsch_1988_crossflow_cd
    mcgreehan_schotsch_1988_expansion_factor = _core.mcgreehan_schotsch_1988_expansion_factor
    mcgreehan_schotsch_1988_expansion_nozzle = _core.mcgreehan_schotsch_1988_expansion_nozzle
    mcgreehan_schotsch_1988_expansion_orifice = _core.mcgreehan_schotsch_1988_expansion_orifice
    orifice_mdot_Cd = _core.orifice_mdot_Cd
    orifice_dP_Cd = _core.orifice_dP_Cd
    orifice_Cd_from_measurement = _core.orifice_Cd_from_measurement
    orifice_K_from_Cd = _core.orifice_K_from_Cd
    orifice_Cd_from_K = _core.orifice_Cd_from_K
    # Orifice flow utilities
    OrificeFlowResult = _core.OrificeFlowResult
    orifice_flow = _core.orifice_flow
    orifice_flow_state = _core.orifice_flow_state
    orifice_velocity_from_mdot = _core.orifice_velocity_from_mdot
    orifice_area_from_beta = _core.orifice_area_from_beta
    beta_from_diameters = _core.beta_from_diameters
    orifice_Re_d_from_mdot = _core.orifice_Re_d_from_mdot
    # Units query API
    UnitInfo = _core.UnitInfo
    get_units = _core.get_units
    input_units = _core.input_units
    output_units = _core.output_units
    has_units = _core.has_units
    list_functions_with_units = _core.list_functions_with_units
    all_units = _core.all_units
    # Channel flow (HTC + pressure drop)
    ChannelResult = _core.ChannelResult
    channel_smooth = _core.channel_smooth
    geometry = _core.geometry

    htc_circular_channel = _core.htc_circular_channel
    nusselt_circular_channel = _core.nusselt_circular_channel
    sharp_area_change = _core.sharp_area_change
    conical_area_change = _core.conical_area_change
    AreaChangeResult = _core.AreaChangeResult
    AreaChangeLoss = _core.AreaChangeLoss
    AreaChangeFlow = _core.AreaChangeFlow
    BordaCarnotStep = _core.BordaCarnotStep
    sudden_expansion_zeta = _core.sudden_expansion_zeta
    sudden_contraction_zeta = _core.sudden_contraction_zeta
    conical_diffuser_zeta = _core.conical_diffuser_zeta
    conical_nozzle_zeta = _core.conical_nozzle_zeta
    compressible_borda_carnot = _core.compressible_borda_carnot
    area_change_loss = _core.area_change_loss
    area_change_flow = _core.area_change_flow
    area_change_element_flow = _core.area_change_element_flow
    compression = _core.compression
    expansion = _core.expansion
    complete_combustion_products_and_jacobian = _core.complete_combustion_products_and_jacobian
    combustion_equilibrium_and_jacobian = _core.combustion_equilibrium_and_jacobian
    CompleteCombustionJacobian = _core.CompleteCombustionJacobian
    EquilibriumJacobian = _core.EquilibriumJacobian
    AreaChangeElementFlow = _core.AreaChangeElementFlow
    TurbomachinePath = _core.TurbomachinePath

# Terminology aliases - always available
try:
    htc_channel = htc_circular_channel
    nusselt_channel = nusselt_circular_channel
    pressure_drop_channel = channel_pressure_drop
    channel_area = circular_area
    channel_volume = cylinder_volume
except NameError:
    # Fallback for complex loading scenarios or missing core extension
    pass

# Submodule imports - always available regardless of _core load path.
from . import compressible, heat_transfer, incompressible, network, species, _solver_tools
from ._flow_solution import FlowSolution
from .vortex import VatistasVortex


# ---------------------------------------------------------------------------
# State: wrap _core.State to accept Cantera-style composition strings
# e.g. set_TPX(300, 101325, "CH4:1, O2:2, N2:7.52")
# ---------------------------------------------------------------------------
import numpy as _np


def _coerce_comp(comp: "_np.ndarray | str") -> "_np.ndarray":
    if isinstance(comp, str):
        return species._parse_comp_str(comp)
    return comp


_StateBase = State


class State(_StateBase):
    def set_X(self, X: "_np.ndarray | str") -> "State":
        return super().set_X(_coerce_comp(X))

    def set_Y(self, Y: "_np.ndarray | str") -> "State":
        return super().set_Y(_coerce_comp(Y))

    def set_TPX(self, T: float, P: float, X: "_np.ndarray | str") -> "State":
        return super().set_TPX(T, P, _coerce_comp(X))

    def set_TPY(self, T: float, P: float, Y: "_np.ndarray | str") -> "State":
        return super().set_TPY(T, P, _coerce_comp(Y))

    @property
    def X(self) -> "_np.ndarray":
        return super().X

    @X.setter
    def X(self, value: "_np.ndarray | str") -> None:
        super().set_X(_coerce_comp(value))

    @property
    def Y(self) -> "_np.ndarray":
        return super().Y

    @Y.setter
    def Y(self, value: "_np.ndarray | str") -> None:
        super().set_Y(_coerce_comp(value))

    @property
    def TPX(self) -> tuple:
        return super().TPX

    @TPX.setter
    def TPX(self, value: tuple) -> None:
        T, P, X = value
        super().set_TPX(T, P, _coerce_comp(X))

    @property
    def TPY(self) -> tuple:
        return super().TPY

    @TPY.setter
    def TPY(self, value: tuple) -> None:
        T, P, Y = value
        super().set_TPY(T, P, _coerce_comp(Y))


# ---------------------------------------------------------------------------
# suppress_warnings context manager
# ---------------------------------------------------------------------------
from collections.abc import Generator as _Generator
from contextlib import contextmanager as _contextmanager


@_contextmanager
def suppress_warnings() -> _Generator[None, None, None]:
    """Context manager that silences all combaero correlation warnings.

    Useful for batch solver runs where out-of-range extrapolation warnings
    would flood the output.  The original handler is restored on exit.

    Example::

        with combaero.suppress_warnings():
            result = solver.run()
    """
    try:
        original = get_warning_handler()
    except Exception:
        original = None
    try:
        set_warning_handler(lambda msg: None)
        yield
    finally:
        if original is not None:
            set_warning_handler(original)
        else:
            set_warning_handler(None)


__all__ = [
    "AccuracyProvenance",
    "RibCorrelationSet",
    "StatedAccuracy",
    "RibGeometry",
    "GShapeModel",
    "RAlphaShape",
    "RatioBaseline",
    "RatioBelowFloor",
    "RibProvenance",
    "RibRatioOptions",
    "RibRatioResult",
    "RibRatioSet",
    "RibShape",
    "RibRange",
    "RibResult",
    "RibTerm",
    "evaluate_rib",
    "evaluate_rib_ratio",
    "han_1988_orthogonal",
    "han_1989_narrow_channel",
    "han_zhang_lee_1991",
    "han_park_1988_angled",
    "rallabandi_2009_high_re",
    "validate_rib_ratio_options",
    "validate_rib_ratio_set",
    "validate_rib_set",
    "taslim_spring_1987",
    "PIN_FIN_RE_EPS",
    "PIN_FIN_F_BLEND",
    "PinArrangement",
    "PinNuSurface",
    "PinReBasis",
    "PinFrictionBasis",
    "PinFinGeometry",
    "PinFinAreaFractions",
    "PinFinEfficiency",
    "PinFinNuSet",
    "PinFinFrictionSet",
    "PinFinRatioModifier",
    "PinFinNuResult",
    "PinFinFrictionResult",
    "PinFinModifierResult",
    "pin_fin_min_gap_D",
    "pin_fin_amin_over_afrontal",
    "pin_fin_dprime_over_D",
    "pin_fin_aprime_over_amin",
    "pin_fin_dh_over_D",
    "pin_fin_area_fractions",
    "pin_fin_array_efficiency",
    "evaluate_pin_fin_nu",
    "evaluate_pin_fin_friction",
    "evaluate_pin_fin_modifier",
    "validate_pin_fin_geometry",
    "validate_pin_fin_nu_set",
    "validate_pin_fin_friction_set",
    "validate_pin_fin_modifier",
    "metzger_1986_staggered_nu",
    "metzger_1982_staggered_friction",
    "vanfossen_1982_staggered_nu",
    "damerow_1972_staggered_friction",
    "chyu_1998_nu",
    "chyu_1998_inline_over_staggered",
    "chyu_1990_nu",
    "chyu_1990_friction",
    "chyu_1990_fillet_over_straight",
    "FLORSCHUETZ_1981_DEFAULT_CD",
    "FLORSCHUETZ_SU_1985_ETA_A",
    "FLORSCHUETZ_SU_1985_ETA_N",
    "FLORSCHUETZ_SU_1985_ETA_R",
    "CombustorLossResult",
    "combustor_fraction_loss",
    "combustor_head_loss",
    "ImpingementEtaR",
    "ImpingementRange",
    "ImpingementThermalBC",
    "JetArrayCorrelationSet",
    "JetArrayGeometricFit",
    "JetArrayImpingementResult",
    "JetHolePattern",
    "JetRowGeometry",
    "JetRowHeatTransfer",
    "StationHalfDrop",
    "CrossflowRatio",
    "JetImpulse",
    "MergeFaceState",
    "SingleJetImpingementResult",
    "SingleJetImpingementSet",
    "crossflow_to_jet_ratio_at_row",
    "crossflow_to_jet_ratio_at_x",
    "florschuetz_1981_inline",
    "florschuetz_1981_staggered",
    "goldstein_1986_single_jet",
    "impingement_eta_r",
    "jet_array_impingement_nu",
    "jet_row_gc_gj_per_flow_ratio",
    "jet_row_heat_transfer",
    "station_half_drop",
    "crossflow_velocity_ratio",
    "FannoMachState",
    "FannoDuctResult",
    "FANNO_MACH_REL_TOL",
    "fanno_state_at_mach",
    "fanno_dx_dmach",
    "fanno_length_between",
    "fanno_inlet_mach",
    "fanno_sonic_mass_flux",
    "fanno_duct",
    "fanno_choked_mass_flux",
    "TransportDerivatives",
    "Omega22Result",
    "transport_and_dT",
    "omega22_and_derivative",
    "dcp_R_dT",
    "set_thermo_model",
    "thermo_model",
    "FannoChannelFlow",
    "FANNO_FLOW_DRIVE_FLOOR",
    "fanno_channel_flow",
    "FannoRayleighMarch",
    "fanno_rayleigh_march",
    "fanno_rayleigh_choked_mass_flux",
    "fanno_rayleigh_channel_flow",
    "CROSSFLOW_VI_FLOOR",
    "channel_entry_drop",
    "segment_area_change_drop",
    "ROGERS_HERSH_1976_CAP_BAND",
    "FLORSCHUETZ_METZGER_SU_1982_CAP_BAND",
    "ROGERS_HERSH_1976_INJECTION_A",
    "ROGERS_HERSH_1976_INJECTION_B",
    "ROGERS_HERSH_1976_SUCTION_A",
    "ROGERS_HERSH_1976_SUCTION_B",
    "CrossflowCd",
    "rogers_hersh_1976_cd",
    "florschuetz_metzger_su_1982_has_plate",
    "florschuetz_metzger_su_1982_cd",
    "Rowbury2001Set",
    "ROWBURY_2001_K_PERP_MIN",
    "ROWBURY_2001_UNIQUE_MARGIN",
    "Rowbury2001Coefficients",
    "Rowbury2001Cd",
    "Rowbury2001KPerp",
    "rowbury_2001_coefficients",
    "rowbury_2001_k_floor",
    "rowbury_2001_in_range",
    "rowbury_2001_exit_crossflow_cd",
    "rowbury_2001_k_perp",
    "rowbury_2001_k_fit_max",
    "rowbury_2001_k_data_max",
    "ROWBURY_2001_FADE_WIDTH",
    "Rowbury2001ExitDensity",
    "rowbury_2001_exit_density",
    "STATION_KAPPA_MERGE_NORMAL",
    "STATION_KAPPA_BLEED_BASSETT",
    "WALL_HTC_KNEE",
    "STEFAN_BOLTZMANN",
    "REEVES_EMISSIVITY_COEFF",
    "LUMINOSITY_HYDROGEN_COEFF",
    "BEAM_LENGTH_COEFF_LEFEBVRE",
    "BEAM_LENGTH_COEFF_HOWELL",
    "GasEmissivityResult",
    "LinerGasRadiationResult",
    "gas_emissivity_luminous",
    "luminosity_factor_hydrogen",
    "liner_gas_radiation",
    "liner_casing_radiation",
    "liner_casing_exchange_factor",
    "mean_beam_length",
    "elemental_fuel_air_ratio",
    "SEBAN1960_C_NEAR",
    "SEBAN1960_C_NEAR_PRINTED",
    "SEBAN1960_C_FAR",
    "SEBAN1960_C_FAR_PRINTED",
    "SLOT_FILM_HTC_HANDBACK_BAND",
    "SLOT_FILM_BRANCH_NEAR",
    "SLOT_FILM_BRANCH_BLEND",
    "SLOT_FILM_BRANCH_FAR",
    "SLOT_FILM_BRANCH_WALL_JET",
    "SLOT_FILM_BRANCH_NONE",
    "SlotFilmStanton",
    "SlotFilmNusselt",
    "SlotFilmHtc",
    "slot_film_stanton_seban_1960",
    "slot_film_nusselt_ballal_lefebvre_1973",
    "slot_film_htc_seban_1960",
    "slot_film_htc_ballal_lefebvre_1973",
    "MIXER_HEAT_MDOT_FLOOR",
    "mixer_heat_per_unit_flow",
    "SENSIBLE_ENTHALPY_REF_T",
    "chamber_merge_face_state",
    "duct_static_state",
    "duct_htc_static_state",
    "duct_htc_static_state_dY",
    "DUCT_HTC_SONIC_MARGIN",
    "jet_impulse",
    "single_jet_impingement",
    "single_jet_impingement_nu",
    "validate_jet_array_set",
    "validate_single_jet_set",
    "FlowSolution",
    "compressible",
    "geometry",
    "heat_transfer",
    "incompressible",
    "network",
    "species",
    "_solver_tools",
    # Core Constants
    "cp",
    "cp_mass",
    "cp_secant_mass_and_jacobian",
    "cp_secant_mass_dY",
    "h",
    "h_mass",
    "s",
    "s_mass",
    "u_mass",
    "cv",
    "cv_mass",
    "u",
    "u_mass",
    "density",
    "specific_gas_constant",
    "isentropic_expansion_coefficient",
    "speed_of_sound",
    "molar_volume",
    "viscosity",
    "thermal_conductivity",
    "kinematic_viscosity",
    "prandtl",
    "thermal_diffusivity",
    "reynolds",
    "peclet",
    "mole_to_mass",
    "mass_to_mole",
    "normalize_fractions",
    "convert_to_dry_fractions",
    "equivalence_ratio_mole",
    "equivalence_ratio",
    "equivalence_ratio_mole",
    "equivalence_ratio_mass",
    "set_equivalence_ratio_mass",
    "equivalence_ratio_from_bilger_Z_mass",
    "bilger_stoich_mixture_fraction_mass",
    "bilger_Z_from_equivalence_ratio_mass",
    "expansibility_factor",
    "solve_orifice_mdot",
    "dry_air",
    "humid_air_composition",
    "humid_air_density",
    "humid_air_enthalpy",
    "humid_air_enthalpy_nasa9",
    "humidity_ratio",
    "saturation_vapor_pressure",
    "vapor_pressure",
    "water_vapor_mole_fraction",
    "wet_bulb_temperature",
    "dewpoint",
    "relative_humidity_from_dewpoint",
    "formula_to_name",
    "name_to_formula",
    "common_name",
    "formula",
    "calc_T_from_h",
    "calc_T_from_h_mass",
    "calc_T_from_s",
    "calc_T_from_s_mass",
    "calc_T_from_cp",
    "calc_T_from_cp_roots",
    "CpRoot",
    "calc_T_from_u",
    "calc_T_from_u_mass",
    "calc_T_from_sv_mass",
    "calc_T_from_sh_mass",
    "oxygen_required_per_mol_fuel",
    "oxygen_required_per_kg_fuel",
    "oxygen_required_per_mol_mixture",
    "oxygen_required_per_kg_mixture",
    "fuel_lhv_molar",
    "fuel_lhv_mass",
    "complete_combustion_to_CO2_H2O",
    "complete_combustion_to_CO2_H2O_with_fraction",
    # State-based API
    "State",
    "MixtureState",
    "Stream",
    "HumidAir",
    "AirProperties",
    "air_properties",
    "ThermoState",
    "thermo_state",
    "TransportState",
    "transport_state",
    "CompleteState",
    "complete_state",
    "EquilibriumResult",
    "CombustionMethod",
    "CombustionState",
    "combustion_state",
    "combustion_state_from_streams",
    "mix",
    # Materials database
    "k_inconel718",
    "k_bondcoat_nicocraly",
    "k_hastelloy_x",
    "k_haynes188",
    "k_haynes230",
    "k_haynes263",
    "k_haynes25",
    "k_inconel625",
    "k_nimonic75",
    "k_copper",
    "k_alumina",
    "k_x22crmov12_1",
    "k_stainless_steel_304",
    "k_stainless_steel_310",
    "k_stainless_steel_316",
    "k_tbc_ysz",
    "list_materials",
    # Advanced cooling correlations
    "thermal_performance_factor",
    "adiabatic_wall_temperature",
    "cooled_wall_heat_flux",
    "equivalent_blowing_ratio",
    "effusion_approach_nusselt",
    "effusion_internal_nusselt",
    "effusion_throat_nusselt",
    "flat_plate_turbulent_nusselt_and_derivative",
    "flat_plate_turbulent_panel_nusselt",
    "flat_plate_turbulent_panel_nusselt_and_derivative",
    "FlatPlatePanelNusselt",
    "mills_entry_length_factor",
    "equivalent_slot_width",
    "film_superposition_corrected",
    "film_superposition_corrected_and_gradient",
    "film_superposition_sellers",
    "effusion_panel_film_effectiveness",
    "EffusionPanelFilm",
    "film_span_average_baldauf_2002",
    "film_span_average_baldauf_2002_groups",
    "film_span_average_colban_2011",
    "slot_film_span_average_marek_tacina_1975",
    "FilmSpanAverage",
    "FilmSpanAverageGroups",
    "film_superposition_sellers_and_gradient",
    "mainstream_temperature_correction",
    "film_effectiveness_baldauf_2002",
    "film_effectiveness_baldauf_2002_and_derivatives",
    "film_effectiveness_colban_2011",
    "film_effectiveness_colban_2011_and_derivatives",
    "slot_film_effectiveness_marek_tacina_1975",
    "slot_film_effectiveness_marek_tacina_1975_and_derivatives",
    # Channel flow (HTC + pressure drop)
    "ChannelResult",
    "channel_smooth",
    "heat_transfer",
    # Inverse solvers - find fuel stream
    "set_fuel_stream_for_Tad",
    "set_fuel_stream_for_phi",
    "set_fuel_stream_for_O2",
    "set_fuel_stream_for_O2_dry",
    "set_fuel_stream_for_CO2",
    "set_fuel_stream_for_CO2_dry",
    # Inverse solvers - find oxidizer stream
    "set_oxidizer_stream_for_Tad",
    "set_oxidizer_stream_for_O2",
    "set_oxidizer_stream_for_O2_dry",
    "set_oxidizer_stream_for_CO2",
    "set_oxidizer_stream_for_CO2_dry",
    "complete_combustion",
    "complete_combustion_isothermal",
    "wgs_equilibrium",
    "wgs_equilibrium_adiabatic",
    "reforming_equilibrium",
    "reforming_equilibrium_adiabatic",
    "combustion_equilibrium",
    # Compressible flow
    "CompressibleFlowSolution",
    "NozzleStation",
    "NozzleSolution",
    "FannoStation",
    "FannoSolution",
    "nozzle_flow",
    "solve_A_eff_from_mdot",
    "solve_P_back_from_mdot",
    "solve_Pt_from_mdot",
    "critical_pressure_ratio",
    "mach_from_pressure_ratio",
    "mass_flux_isentropic",
    "nozzle_cd",
    "fanno_max_length",
    # Thrust
    "ThrustResult",
    "nozzle_thrust",
    "nozzle_thrust_cd",
    # Friction
    "friction_haaland",
    "friction_serghides",
    "friction_colebrook",
    "friction_petukhov",
    "is_whistling_risk",
    "num_species",
    "species_name",
    "species_index_from_name",
    "species_molar_mass",
    "species_molar_mass_from_name",
    # Heat transfer
    "nusselt_dittus_boelter",
    "nusselt_gnielinski",
    "nusselt_sieder_tate",
    "nusselt_petukhov",
    "nusselt_circular_channel",
    "nusselt_channel",
    "htc_from_nusselt",
    "htc_circular_channel",
    "htc_channel",
    "orifice_impedance_with_flow",
    "quarter_wave_resonator_tmm",
    "TransferMatrix",
    "CanAnnularGeometry",
    "BlochMode",
    "can_annular_eigenmodes",
    "AnnularDuctGeometry",
    "AnnularMode",
    "annular_duct_eigenmodes",
    "annular_duct_modes_analytical",
    "AcousticMedium",
    "LinerOrificeGeometry",
    "LinerFlowState",
    "LinerCavity",
    "absorption_from_impedance_norm",
    "liner_sdof_impedance_norm",
    "liner_sdof_absorption",
    "sweep_liner_sdof_absorption",
    "liner_2dof_serial_impedance_norm",
    "liner_2dof_serial_absorption",
    "sweep_liner_2dof_serial_absorption",
    "overall_htc",
    "overall_htc_wall",
    "overall_htc_wall_multilayer",
    "thermal_resistance",
    "thermal_resistance_wall",
    "lmtd",
    "lmtd_counterflow",
    "lmtd_parallelflow",
    "heat_rate",
    "heat_flux",
    "heat_transfer_area",
    "heat_transfer_dT",
    "wall_temperature_profile",
    "ntu",
    "capacity_ratio",
    "effectiveness_counterflow",
    "effectiveness_parallelflow",
    "heat_rate_from_effectiveness",
    "heat_flux_from_T_at_edge",
    "heat_flux_from_T_at_depth",
    "bulk_T_from_edge_T_and_q",
    "dT_edge_dT_hot",
    "dT_edge_dT_cold",
    "dT_edge_dT_bulk",
    "dT_edge_dq",
    # Acoustics
    "Tube",
    "Annulus",
    "CanAnnularFlowGeometry",
    "BoundaryCondition",
    "AcousticMode",
    "to_acoustic_geometry",
    "tube_axial_modes",
    "annulus_axial_modes",
    "annulus_azimuthal_modes",
    "annulus_modes",
    "modes_in_range",
    "closest_mode",
    "min_mode_separation",
    # Mean flow correction
    "axial_mode_upstream",
    "axial_mode_downstream",
    "axial_mode_split",
    # Helmholtz resonator
    "helmholtz_frequency",
    # Strouhal
    "strouhal",
    "frequency_from_strouhal",
    # Convenience
    "quarter_wave_frequency",
    "half_wave_frequency",
    # Viscothermal
    "stokes_layer",
    "thermal_layer",
    "effective_viscothermal_layer",
    # Quality factor
    "helmholtz_Q",
    "tube_Q",
    "damping_ratio",
    "bandwidth",
    # Acoustic properties
    "AcousticProperties",
    "acoustic_properties",
    "wavelength",
    "frequency_from_wavelength",
    "acoustic_impedance",
    "sound_pressure_level",
    "particle_velocity",
    # Residence time
    "residence_time",
    "residence_time_tube",
    "residence_time_annulus",
    "residence_time_can_annular",
    "residence_time_mdot",
    "residence_time_mdot_can_annular",
    "space_velocity",
    # Incompressible flow
    "bernoulli_P2",
    "bernoulli_v2",
    "orifice_mdot",
    "orifice_Q",
    "orifice_velocity",
    "orifice_area",
    "orifice_dP",
    "channel_dP",
    "channel_dP_mdot",
    "channel_velocity",
    "channel_mdot",
    "channel_pressure_drop",
    "pressure_drop_channel",
    "channel_area",
    "channel_volume",
    "channel_roughness",
    "standard_channel_roughness",
    "circular_area",
    "cylinder_volume",
    "dynamic_pressure",
    # Orifice Cd correlations
    "OrificeGeometry",
    "OrificeState",
    "MeteringCdCorrelation",
    "OrificeTap",
    "DischargeCdCorrelation",
    "DischargeHoleGeometry",
    "DischargeHoleState",
    "CorrelationValidity",
    "CorrelationResult",
    "Cd_sharp_thin_plate",
    "Cd_sharp_thin_plate_and_dRe",
    "metering_cd_and_dRe",
    "metering_cd_in_range",
    "metering_cd_re_min",
    "permanent_loss_ratio_and_dC",
    "velocity_of_approach_factor",
    "jobson_force_defect",
    "jobson_ci_from_cd",
    "jobson_contraction_coefficient",
    "jobson_expansibility",
    "JOBSON_1955_BETA_MAX",
    "discharge_cd",
    "discharge_cd_and_derivatives",
    "discharge_cd_in_range",
    "discharge_cd_re_floor",
    "discharge_inclination_dzeta",
    "discharge_pipe_cd_ratio",
    "discharge_pipe_zeta",
    "mcgreehan_schotsch_1988_cd",
    "mcgreehan_schotsch_1988_cd_and_derivatives",
    "mcgreehan_schotsch_1988_crossflow_cd",
    "mcgreehan_schotsch_1988_expansion_factor",
    "mcgreehan_schotsch_1988_expansion_nozzle",
    "mcgreehan_schotsch_1988_expansion_orifice",
    "orifice_mdot_Cd",
    "orifice_dP_Cd",
    "orifice_Cd_from_measurement",
    "orifice_K_from_Cd",
    "orifice_Cd_from_K",
    # Orifice flow utilities
    "OrificeFlowResult",
    "orifice_flow",
    "orifice_flow_state",
    "orifice_velocity_from_mdot",
    "orifice_area_from_beta",
    "beta_from_diameters",
    "orifice_Re_d_from_mdot",
    # Units query API
    "UnitInfo",
    "get_units",
    "input_units",
    "output_units",
    "has_units",
    "list_functions_with_units",
    "all_units",
    "sharp_area_change",
    "conical_area_change",
    "AreaChangeResult",
    "AreaChangeLoss",
    "AreaChangeFlow",
    "BordaCarnotStep",
    "sudden_expansion_zeta",
    "sudden_contraction_zeta",
    "conical_diffuser_zeta",
    "conical_nozzle_zeta",
    "compressible_borda_carnot",
    "area_change_loss",
    "area_change_flow",
    "area_change_element_flow",
    "compression",
    "expansion",
    "complete_combustion_products_and_jacobian",
    "combustion_equilibrium_and_jacobian",
    "CompleteCombustionJacobian",
    "EquilibriumJacobian",
    "AreaChangeElementFlow",
    "TurbomachinePath",
    # Vatistas n-vortex model
    "vatistas_v0_bar",
    "vatistas_dv0_bar_drbar",
    "vatistas_vr_bar",
    "vatistas_dvr_bar_drbar",
    "vatistas_pressure_integral",
    "vatistas_d_pressure_integral_drbar",
    "vatistas_v_theta",
    "vatistas_v_theta_and_jacobians",
    "vatistas_delta_p",
    "vatistas_delta_p_and_jacobians",
    "VatistasVortex",
    "__version__",
    # -- reachable on the module but previously undeclared ----------------
    #
    # These were exported all along; `__all__` is hand-maintained and had
    # simply fallen behind, which also hid them from test_units_sync.
    # Grouped here rather than merged above so the omission stays legible.
    "Cd_from_zeta",
    "ChannelSolverResult",
    "IncompressibleFlowSolution",
    "MassStream",
    "EntropyGenerationResult",
    "entropy_generation",
    "heat_entropy_flux",
    "MomentumChamberResult",
    "OrificeResult",
    "OrificeTapPressure",
    "Pt_from_static",
    "P_from_stagnation",
    "Registry",
    "Tt_from_static",
    "Tt_from_static_v",
    "T_adiabatic_wall",
    "T_adiabatic_wall_mach",
    "T_from_stagnation",
    "WallCouplingResult",
    "adiabatic_T_wgs",
    "annular_area",
    "channel_flow",
    "channel_flow_rough",
    "dcp_dT",
    "dh_dT",
    "ds_dT",
    "fanno_channel",
    "fanno_channel_rough",
    "get_material_conductivity",
    "get_material_conductivity_and_slope",
    "get_material_info",
    "material_in_range",
    "GENERIC_METAL_SOLIDUS",
    "get_warning_handler",
    "h0_from_static",
    "hydraulic_diameter",
    "hydraulic_diameter_annulus",
    "hydraulic_diameter_rect",
    "is_well_behaved",
    "kinetic_energy",
    "bulk_velocity",
    "stagnation_from_static",
    "mach_number",
    "mixture_h",
    "molar_mass",
    "mwmix",
    "orifice_flow_thermo",
    "pressure_loss",
    "recovery_factor",
    "set_equivalence_ratio_mole",
    "set_warning_handler",
    "suppress_warnings",
    "v_from_h0",
    "velocity_from_pressure_loss",
    "velocity_from_q",
    "vortex",
    "zeta_from_Cd",
]
