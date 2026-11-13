"""Unit catalogue for result-export column headers.

The SI unit of every result key, for the CSV produced by ``/export``
(``P [Pa]``, ``T [K]`` ...) and, through ``unit_of`` and
``scripts/gui_result_keys.py``, for the GUI: the frontend's
``QUANTITY_CATALOGUE`` display units are checked against these
(``quantities.test.ts``).

Keep this module dependency-free; it is imported by ``main.py`` during the
CSV-streaming response.
"""

import re

# ---------------------------------------------------------------------------
# Scalar quantity -> unit mapping
# ---------------------------------------------------------------------------

UNIT_MAP: dict[str, str] = {
    # --- Thermodynamic state ---
    "T": "K",
    "T_total": "K",
    "T_in": "K",
    "T_out": "K",
    "Tt_in": "K",
    "Tt_out": "K",
    "T_aw": "K",
    "T_hot": "K",
    "T_interface": "K",
    "P": "Pa",
    "P_total": "Pa",
    "P_in": "Pa",
    "P_out": "Pa",
    "Pt_in": "Pa",
    "Pt_out": "Pa",
    "dP": "Pa",
    "dPt": "Pa",
    "pr_total": "-",
    "pr_static": "-",
    # --- Flow ---
    "m_dot": "kg/s",
    "velocity": "m/s",
    "v": "m/s",
    "a": "m/s",
    "mach": "-",
    "mach_in": "-",
    "mach_out": "-",
    "mach_throat": "-",
    "mach_loss_ref": "-",
    "dPt_loss": "Pa",  # AreaChangeElement: the closure's stagnation loss (#731)
    "exit_head_recovered": "Pa",  # head a node at rest recovers from a Pt-coupled exit
    "choked": "-",
    # --- Node stagnation ---
    "Tt": "K",
    "Pt": "Pa",
    # --- Transport ---
    "mu": "Pa\u00b7s",
    "k": "W/m/K",
    "nu": "m\u00b2/s",
    "Pr": "-",
    "Re": "-",
    "Dh": "m",
    # --- Thermo per-mass (node/element state) ---
    "rho": "kg/m\u00b3",
    "h": "J/kg",
    "s": "J/kg/K",
    "u": "J/kg",
    "cp": "J/kg/K",
    "cv": "J/kg/K",
    "gamma": "-",
    "mw": "g/mol",  # State::mw (include/units_data.h)
    # --- Thermo reference-state block (Phase 2, element inlet state) ---
    "rho_ref": "kg/m\u00b3",
    "h_ref": "J/kg",
    "s_ref": "J/kg/K",
    "u_ref": "J/kg",
    "cp_ref": "J/kg/K",
    "cv_ref": "J/kg/K",
    "gamma_ref": "-",
    "a_ref": "m/s",
    "mw_ref": "g/mol",
    "mu_ref": "Pa\u00b7s",
    "k_ref": "W/m/K",
    "nu_ref": "m\u00b2/s",
    "Pr_ref": "-",
    # --- Heat transfer ---
    "Nu": "-",
    "htc": "W/m\u00b2/K",
    "h_a": "W/m\u00b2/K",
    "h_b": "W/m\u00b2/K",
    "Q": "W",
    "Q_to_boundary": "W",  # wall heat that left with a stream into a boundary
    "Q_wall_out": "W",
    "Q_fraction": "W",  # heat an EnergyBoundary fraction applied
    "Q_withheld": "W",
    "Tt_target": "K",  # heat source (#792): the set outlet TOTAL temperature
    # --- Tee junction ---
    "m_dot_com": "kg/s",
    "m_dot_straight": "kg/s",
    "m_dot_branch": "kg/s",
    "mass_flow_ratio": "kg/kg",
    "K_straight": "-",
    "K_branch": "-",
    "n_reversed_ports": "-",  # ideal momentum tee: ports run against the wiring
    # --- Element diagnostics ---
    "xi": "-",
    "theta": "-",
    "fraction_unanchored": "-",  # a fraction loss without m_ref on its zero-flow blend (#1004)
    "closed": "-",  # 1 for a restriction of zero area: a wall, no flow (#1024)
    "zeta": "-",
    "area_ratio": "-",  # Phase 2 name
    "is_correlation": "-",
    # --- The hole-Cd group (hole_cd.HOLE_CD_DIAGNOSTIC_KEYS): the same keys
    # on every hole element; a term that is off reads applied 0. The
    # *_model / Cd_correlation keys are text (TEXT_KEYS).
    "Cd_base": "-",  # the base Cd, before any exit-crossflow term
    "Cd_in_range": "-",  # 1 inside the discharge correlation's source range
    # 1 per variable outside the hole Cd's source range (info, not a warning)
    "Cd_extrapolated_Re": "-",
    "Cd_extrapolated_l_d": "-",
    "Cd_extrapolated_edge": "-",
    "Cd_inclination_applied": "-",
    "Cd_inclination_dzeta": "-",  # inclined-inlet increment on zeta' (0 when off or at 90 deg)
    "Cd_inlet_crossflow_applied": "-",
    "Cd_inlet_crossflow_U1_over_Vi": "-",  # a duct-fed hole's U1/Vi, 0 plenum-fed
    "Cd_inlet_crossflow_U1_over_Vi_max": "-",  # effusion liner: max over panels
    "Cd_inlet_crossflow_beyond_8pct": "-",
    "Cd_inlet_crossflow_degraded": "-",
    "Cd_inlet_crossflow_ignored": "-",
    "Cd_supply_drive_applied": "-",  # 1 where a duct-fed hole drives from its duct total (#1087)
    "Cd_supply_drive_head_fraction": "-",  # duct head over the full drive (info), 0 not duct-fed
    "Cd_supply_recovery_q1": "Pa",  # a duct-fed hole's supply head Pt1 - P1, 0 not duct-fed
    "Cd_supply_recovery_fraction": "-",  # recovered fraction f applied (lean recovery: #1124)
    "Cd_supply_recovery_ratio": "-",  # Cd_T / Cd_base from P1 + f q1 onto the total, 1 when none
    "Cd_exit_crossflow_applied": "-",
    "Cd_exit_crossflow_ratio": "-",  # Cd / Cd_base, 1 when off
    "Cd_exit_crossflow_rho_u2": "Pa",  # mainstream momentum flux rho_m u_m^2
    "Cd_exit_crossflow_k_perp": "-",
    "Cd_exit_crossflow_k_used": "-",
    "Cd_exit_crossflow_k_floor": "-",
    "Cd_exit_crossflow_faded": "-",
    "Cd_exit_crossflow_floor_above_fit": "-",
    "Cd_exit_crossflow_delta": "-",
    "Cd_exit_crossflow_held": "-",
    "Cd_exit_crossflow_in_range": "-",
    "Cd_exit_crossflow_extrapolated_velocity_ratio": "-",  # Rogers & Hersh Vi/V_inf
    "Cd_exit_crossflow_extrapolated_l_d": "-",
    "Cd_exit_crossflow_extrapolated_grazing_velocity": "-",
    # Deprecated since 0.12.0, removed in 0.13.0 (hole_cd.HOLE_CD_DIAGNOSTIC_ALIASES):
    # still emitted beside the group; nothing in the GUI reads them.
    "Cd_own": "-",
    "U1_over_Vi": "-",
    "U1_over_Vi_max": "-",
    "coolant_crossflow_ignored": "-",
    "crossflow_cd_beyond_8pct": "-",
    "crossflow_cd_degraded": "-",
    # --- Impingement plate / crossflow (#465) ---
    "n_holes": "-",
    "Re_j": "-",
    "Gc_Gj": "-",
    "Gc_Gj_closed_form": "-",
    "St": "-",
    "surface_extrapolated": "-",
    "material_in_range": "-",
    "dP_momentum": "Pa",
    # --- Effusion plate wall (#471) ---
    "h_internal_plate_area": "W/m\u00b2/K",
    "h_gas_unblown": "W/m\u00b2/K",
    "h_gas": "W/m\u00b2/K",
    "gas_augmentation": "-",
    "T_gas": "K",
    "eta_film": "-",
    # Film-side h near a slot (#962), on the wall the film covers
    "film_h": "W/m\u00b2/K",
    "film_h_duct": "W/m\u00b2/K",
    "film_h_film": "W/m\u00b2/K",
    "film_h_fraction_in_box": "-",
    "film_h_weight_M": "-",
    "film_h_Re_s": "-",
    "film_h_branch": "-",
    "film_h_in_box": "-",
    "film_h_x_in_box": "-",
    "film_h_Re_in_box": "-",
    "film_h_M_in_box": "-",
    "eta_overall": "-",
    "T_adiabatic_wall": "K",
    "T_wall_hot": "K",
    "T_wall_cold": "K",
    "q_wall": "W/m\u00b2",
    "Q_wall": "W",
    "blowing_ratio": "-",
    "velocity_ratio": "-",
    "density_ratio": "-",
    "U_gas": "m/s",
    "gas_heat_flux_imposed": "W/m\u00b2",
    # --- Effusion liner (#471) ---
    "m_dot_out": "kg/s",
    "m_dot_bleed": "kg/s",
    "bleed_fraction": "-",
    "n_segments": "-",
    "bleed_nonuniformity": "-",
    "station_psi_min": "-",
    "T_wall_hot_max": "K",
    "n_rows": "-",
    "Re_j_min": "-",
    "Re_j_max": "-",
    "Gc_Gj_max": "-",
    "Nu_mean": "-",
    "htc_mean": "W/m\u00b2/K",
    "jet_flow_nonuniformity": "-",
    # --- Impingement array bypass (#467) ---
    "m_dot_bypass": "kg/s",
    "bypass_to_jet_flow": "-",
    "bypass_Tt_minus_supply": "K",
    "heat_transfer_outside_source": "-",
    "bypass_reversed": "-",
    "bypass_fraction": "-",
    "eta_r": "-",  # CR-3936 recovery-temperature effectiveness (#713)
    "eta_r_extrapolated": "-",
    "f": "-",
    "Cd": "-",
    "phi": "-",
    # --- Fuel controller (#892); its Tt_target shares the heat source's ---
    "Tt_without_fuel": "K",
    "rich_branch": "-",
    # --- Compressor / turbine (#893) ---
    "power": "W",
    "dh": "J/kg",
    "Tt_out_isentropic": "K",
    "eta_poly": "-",
    "eta_is": "-",
    "pressure_ratio": "-",
    "ln_m_dot": "-",
    # --- Flow splitter (#903): per outlet, and the node's own ---
    "fraction": "-",
    "fraction_achieved": "-",
    "m_dot_splitter": "kg/s",
    "is_reference": "-",
    "valve_dPt": "Pa",
    "valve_CdA": "m\u00b2",
    "valve_Cd": "-",
    "valve_zeta": "-",
    "valve_choked": "-",
    "valve_needs_pump": "-",
    "n_outlets": "-",
    "reference_outlet": "-",
    "reference_start": "-",
    "reference_swapped": "-",
    "reference_swaps": "-",
    "n_needs_pump": "-",
    "valve_dPt_max": "Pa",
    "fraction_error_max": "-",
    # --- Per-row lists of the arrays (``rows_<key>``, unit of <key>) ---
    "rows_P_backside": "Pa",  # effusion liner: each panel's P_in
    "rows_unresolved": "-",  # impingement array: 1.0 marks a starved row
    # --- Geometry ---
    "area": "m\u00b2",
    "diameter": "m",
    "length": "m",
    "roughness": "m",
    # Until 0.8.0 the keys below had a unit in the GUI only and an unlabelled
    # export column (#774). Flags (0/1) and counts are "-".
    # --- Solver views ---
    "sigma_norm": "-",
    "sigma_norm_high": "-",
    "sigma_norm_low": "-",
    "incompressible_mach": "-",
    "incompressible_mach_limit": "-",
    # Channel / Fanno-Rayleigh duct (#551, #636)
    "M_in_duct": "-",
    "M_exit_duct": "-",
    "f_surface": "-",
    "Re_surface": "-",
    "heat_area_ratio": "-",
    "heat_area_exceeds_wetted": "-",
    "rayleigh_in_momentum": "-",
    "rayleigh_omitted": "-",
    "rayleigh_omitted_frac": "-",
    # Area change (#731)
    "F0": "m\u00b2",
    "F1": "m\u00b2",
    "Re_small": "-",
    "expansion": "-",
    # Orifice approach and ISO 5167 taps (#649/#645/#651)
    "beta": "-",
    "beta_flow": "-",
    "E": "-",
    "Cd_pipe_ratio": "-",
    "P_approach": "Pa",
    "dP_tap": "Pa",
    "dP_permanent": "Pa",
    "loss_ratio": "-",
    "pressure_recovery": "-",
    "recovery_in_range": "-",
    # Opt-in compressible orifice models (#854)
    "Ci": "-",
    "contraction_C": "-",
    "expansibility": "-",
    "model_in_range": "-",
    "Cv": "-",
    # Discharge holes, film and slot (#386, #714)
    "G_coolant": "kg/m\u00b2/s",
    "G_jet": "kg/m\u00b2/s",
    "G_slot": "kg/m\u00b2/s",
    "dP_drive": "Pa",
    "L_over_d": "-",
    "hole_length": "m",
    "hydraulic_diameter": "m",
    "slot_height": "m",
    # A fraction on every element (#773); the slot film reports it.
    "turbulence_intensity": "-",
    "is_ingesting": "-",
    # Effusion plate (#471)
    "Nu_approach": "-",
    "Nu_internal": "-",
    "Nu_throat": "-",
    "Re_hole": "-",
    "h_hole_area": "W/m\u00b2/K",
    "h_plate_area": "W/m\u00b2/K",
    "h_gas_unblown_Re": "-",
    "h_gas_unblown_Pr": "-",
    "h_gas_unblown_Nu": "-",
    # The flat-plate gas-side baseline: where its panel average was taken
    "h_gas_unblown_x_start": "m",
    "h_gas_unblown_x_end": "m",
    "h_gas_unblown_Re_x_start": "-",
    "h_gas_unblown_Re_x_end": "-",
    "h_gas_unblown_extrapolated": "-",
    # Gas radiation on the plate's hot face (#1071)
    "q_rad_hot": "W/m\u00b2",
    "eps_g": "-",
    "T_gas_from_main_inlet": "-",
    "film_extrapolated": "-",
    "gas_side_chamber": "-",
    "area_approach_total": "m\u00b2",
    "area_hole_total": "m\u00b2",
    "approach_to_hole_area_ratio": "-",  # approach area / hole area (#773)
    "coolant_ht_plenum_assumed": "-",
    "hole_count_exact": "-",
    "internal_Nu_multiplier": "-",
    "pitch_over_d": "-",
    "porosity": "-",
    "wall_conductivity": "W/m/K",
    "T_wall_limit": "K",  # the plate material's envelope hot end (#998)
    "wall_above_limit": "-",
    "no_cooling_path": "-",
    # Effusion liner (#471)
    "station_psi_beyond_bassett": "-",
    # Impingement array (#465)
    "n_rows_unresolved": "-",
    # Tee and chamber junctions (#517, #536)
    "K11": "-",
    "K12": "-",
    "P_jct": "Pa",
    "Pt_jct": "Pa",
    "n_ports": "-",
    "port_area_mismatch": "-",
    # Ejector (#270, #619)
    "m_dot_primary": "kg/s",
    "m_dot_secondary": "kg/s",
    "m_dot_outlet": "kg/s",
    "entrainment_ratio": "-",  # m_dot_secondary / m_dot_primary (#773)
    "entrainment_ratio_critical": "-",
    "critical_mode": "-",
    "nozzle_exit_mach_out_of_range": "-",  # M_p1 clamped: +1 above, -1 below (#945)
    "s_choke": "-",
    "r_gas": "J/kg/K",
    "P_c_star": "Pa",
    "P_b0": "Pa",
    # Incompressible ejector: Sanger's jet pump (#588)
    "regime_incompressible": "-",
    "mach_jet": "-",
    "head_ratio_N": "-",
    "efficiency": "-",
    # Vortex (#661)
    "shaft_speed_rpm": "rpm",
    "v_theta_in": "m/s",
    "v_theta_out": "m/s",
    "M_theta_in": "-",
    "M_theta_out": "-",
    "P_swirl_in": "Pa",
    "P_swirl_out": "Pa",
    "dP_swirl_static": "Pa",
    # Thermal wall exchange (#627)
    "T_aw_a": "K",
    "T_aw_b": "K",
    "UA": "W/K",
    "NTU": "-",
    "effectiveness": "-",
    "counterflow": "-",
    # Thermal wall radiation (#793), per face: flux on that face, rate,
    # surface temperature, gas emissivity
    "q_rad_a": "W/m\u00b2",
    "q_rad_b": "W/m\u00b2",
    "Q_rad_a": "W",
    "Q_rad_b": "W",
    "T_surface_a": "K",
    "T_surface_b": "K",
    "eps_g_a": "-",
    "eps_g_b": "-",
}

# Per-port keys of the junctions, ``port_<i>_<quantity>`` (i = the element's
# port index) and the tee's hidden port steps, ``step_<port>_<quantity>``.
_PORT_UNITS: dict[str, str] = {
    "P": "Pa",
    "T": "K",  # static (#782, 0.8)
    "Tt": "K",
    "m_dot": "kg/s",
    "area": "m\u00b2",
    "outer_area": "m\u00b2",
    "connected_area_ratio": "-",  # connected element area / port area (#773)
    "sign": "-",
    "K": "-",
    "flow_fraction": "-",
}
_STEP_UNITS: dict[str, str] = {"zeta": "-", "dPt": "Pa", "connected_area_ratio": "-"}
_PORT_RE = re.compile(r"^port_\d+_(\w+)$")
_STEP_RE = re.compile(r"^step_(?:common|straight|branch)_(\w+)$")
# The arrays' per-row lists (impingement array, effusion liner): ``rows_<key>``
# carries <key>'s unit.
_ROWS_RE = re.compile(r"^rows_(\w+)$")


def unit_of(key: str) -> str | None:
    """The SI unit of a result key (UNIT_MAP, then the ``port_<i>_*``,
    ``step_<port>_*`` and ``rows_*`` patterns), or None when none is known.
    The export labels every column this returns a unit for; the GUI shows the
    same."""
    if key in UNIT_MAP:
        return UNIT_MAP[key]
    m = _PORT_RE.match(key)
    if m and m.group(1) in _PORT_UNITS:
        return _PORT_UNITS[m.group(1)]
    m = _STEP_RE.match(key)
    if m and m.group(1) in _STEP_UNITS:
        return _STEP_UNITS[m.group(1)]
    m = _ROWS_RE.match(key)
    if m and m.group(1) in UNIT_MAP:
        return UNIT_MAP[m.group(1)]
    return None


# Result keys whose value is text (a source, a model name, a verdict), not a
# quantity: their export column carries no unit. Listed so that a test can
# tell them from a numeric key nobody gave a unit (test_gui_export_units.py).
TEXT_KEYS: frozenset[str] = frozenset(
    {
        "area_source",  # where an element's area came from (#773)
        # The hole-Cd group's model names (hole_cd.HOLE_CD_DIAGNOSTIC_TEXT_KEYS).
        "Cd_correlation",
        "Cd_inclination_model",
        "Cd_inlet_crossflow_model",
        "Cd_supply_drive_model",
        "Cd_exit_crossflow_model",
        "compressible_model",  # the orifice expansion model (#854)
        "F0_source",
        "F1_source",
        "h_gas_unblown_model",
        "incompressible_validity",
        "second_law",  # the #545 verdict and its context
        "second_law_group",
        "second_law_reason",
        "tap",  # ISO 5167 tap arrangement
    }
)
_PORT_TEXT = frozenset({"area_source"})


def is_text_key(key: str) -> bool:
    """True for a result key whose value is text and has no unit."""
    if key in TEXT_KEYS:
        return True
    m = _PORT_RE.match(key)
    return bool(m and m.group(1) in _PORT_TEXT)


# Composition keys render as ``Y[<label>] [-]``. ``<label>`` may be a species
# name (e.g. ``N2``, ``H2O``, ``CH4``) or a numeric index for forward
# compatibility with older exports.
_COMPOSITION_RE = re.compile(r"^[YX]\[[A-Za-z0-9_+\-()]+\]$")


# ---------------------------------------------------------------------------
# Input (network data) key -> unit mapping (#784)
# ---------------------------------------------------------------------------

# The unit of every numeric field a node, edge or the solver settings can
# carry (gui/backend/schemas.py), by field name. A nested field
# (``composition.ambient_T``, ``surface.d_jet``, ``layers.0.thickness``) takes
# its last segment's unit. Kept apart from UNIT_MAP because some names mean
# something else as an input: a discrete loss's ``k`` is a dimensionless slope,
# not a conductivity. test_results_table.py fails on a numeric field missing
# here.
INPUT_UNIT_MAP: dict[str, str] = {
    # Boundary conditions
    "m_dot": "kg/s",
    "Pt": "Pa",
    "Tt": "K",
    "Q": "W",
    "Tt_target": "K",
    "fraction": "-",
    "m_dot_guess": "kg/s",  # fuel boundary (#892)
    "area_a": "m^2",  # thermal wall faces (#887)
    "area_b": "m^2",
    "pressure_ratio": "-",  # compressor / turbine (#893)
    "efficiency": "-",
    # Composition and ambient
    "relative_humidity": "-",
    "ambient_T": "K",
    "ambient_P": "Pa",
    "ambient_RH": "-",
    # Geometry: lengths
    "L": "m",
    "D": "m",
    "Dh": "m",
    "D_h": "m",
    "diameter": "m",
    "length": "m",
    "width": "m",
    "span": "m",
    "roughness": "m",
    "thickness": "m",
    "wall_thickness": "m",
    "plate_thickness": "m",
    "lip_thickness": "m",
    "edge_radius": "m",
    "bevel_depth": "m",
    "duct_height": "m",
    "slot_height": "m",
    "hole_diameter": "m",
    "d_jet": "m",
    "pin_diameter": "m",
    "channel_height": "m",
    "pitch_x": "m",
    "pitch_y": "m",
    "lateral_pitch": "m",
    "row_pitch": "m",
    "breakout_width": "m",
    "panel_length": "m",
    "panel_width": "m",
    "gas_start_length": "m",
    "bypass_length": "m",
    "r_in": "m",
    "r_out": "m",
    "film_start": "m",
    "film_end": "m",
    # Geometry: areas
    "area": "m^2",
    "F0": "m^2",
    "F1": "m^2",
    "F_C": "m^2",
    "F_branch": "m^2",
    "throat_area": "m^2",
    "nozzle_exit_area": "m^2",
    "mixing_area": "m^2",
    # Angles
    "angle_deg": "deg",
    "exit_crossflow_orientation_deg": "deg",  # Rowbury 2001 compound angle alpha (#1074)
    "alpha_deg": "deg",
    "theta_deg": "deg",
    # Speeds
    "shaft_speed_rpm": "rpm",
    "v_theta_in": "m/s",
    # Materials and heat transfer
    "conductivity": "W/m/K",
    "wall_conductivity": "W/m/K",
    "k_pin": "W/m/K",
    "R_fouling": "m^2*K/W",
    "gas_heat_flux": "W/m^2",
    # Radiation on a wall face or a plate's hot face (#1071)
    "beam_length": "m",
    "emissivity_wall": "-",
    "emissivity_casing": "-",
    "luminosity": "-",
    "T_gas": "K",
    "fuel_air_ratio": "-",
    "pressure": "Pa",
    "T_casing": "K",
    # Flow references
    "m_ref": "kg/s",
    # Dimensionless: coefficients, ratios, multipliers, counts and flags
    "Cd": "-",
    "supply_recovery": "-",  # HoleCd: recovered fraction f of a duct-fed hole's supply head
    "Cv": "-",  # orifice: Jobson velocity coefficient
    "C_D": "-",
    "xi": "-",
    "xi0": "-",
    "k": "-",  # discrete loss: slope of xi or zeta in theta
    "zeta": "-",
    "zeta0": "-",
    "entry_K": "-",
    "bypass_entry_K": "-",
    "K_straight": "-",
    "K_branch": "-",
    "k_nozzle": "-",
    "k_suction": "-",
    "k_throat": "-",
    "k_diffuser": "-",
    "loss_factor": "-",
    "psi": "-",
    "design_branch_share": "-",
    "joining_etransfer_alpha": "-",
    "recovery_efficiency": "-",
    "area_ratio": "-",
    "turbulence_intensity": "-",
    "gas_augmentation": "-",
    "Nu_multiplier": "-",
    "f_multiplier": "-",
    "internal_Nu_multiplier": "-",
    "smooth_wall_Nu_multiplier": "-",
    "C_m_multiplier": "-",
    "e_D": "-",
    "p_e": "-",
    "W_H": "-",
    "xn_d": "-",
    "yn_d": "-",
    "z_d": "-",
    "L_D": "-",
    "R_D": "-",
    "S_D": "-",
    "X_D": "-",
    "n_rows": "-",
    "n_segments": "-",
    "holes_per_row": "-",
    "n_ribbed_walls": "-",
    "N_rows": "-",
    "row": "-",
    "is_staggered": "-",
    "in_pipe": "-",
    # Flow splitter (#903): fractions[i], the reference outlet's index, areas[i]
    "fractions": "-",
    "reference": "-",
    "auto_reference": "-",
    "areas": "m^2",
}


def input_unit(key: str) -> str | None:
    """The ASCII unit of an input key (a network data field, dotted when
    nested), or None for a text field or an unknown one. Custom mixture
    fractions (``composition.custom_fractions.<species>``) are ``-``."""
    if ".custom_fractions." in f".{key}":
        return "-"
    return INPUT_UNIT_MAP.get(key.rsplit(".", 1)[-1])


def column_unit(column: str) -> str | None:
    """The ASCII unit an export column for ``column`` is labelled with
    (``kg/m^3``, ``Pa*s``, ``-`` for a composition), or None for a meta
    column, a text key or a key with no known unit."""
    if column in _META_COLUMNS:
        return None
    unit = unit_of(column)
    if unit is not None:
        return unit.translate(_ASCII_UNITS)
    if _COMPOSITION_RE.match(column):
        return "-"
    return None


def label_with_unit(column: str) -> str:
    """Return ``'<column> [<unit>]'`` when a unit is known, else ``column``.

    Meta columns such as ``id``, ``type``, ``label``, ``kind``, the solve
    status and ``is_boundary`` are passed through unchanged. Units are written
    in ASCII (``kg/m^3``, ``Pa*s``).
    """
    unit = column_unit(column)
    return column if unit is None else f"{column} [{unit}]"


# Export headers are ASCII: the CSV is UTF-8 without a byte-order mark, and
# Excel read "kg/m\u00b3" as mojibake. The GUI keeps the typeset units
# (UNIT_MAP is also the source of resultKeys.generated.json).
_ASCII_UNITS = str.maketrans({"\u00b2": "^2", "\u00b3": "^3", "\u00b7": "*"})

_META_COLUMNS: frozenset[str] = frozenset(
    {
        "id",
        "label",
        "parent",
        "type",
        "kind",
        "success",
        "converged",
        "consistent",
        "outcome",
        "final_norm",
        "is_boundary",
        "combustion_method",
        "ref_location",
    }
)
