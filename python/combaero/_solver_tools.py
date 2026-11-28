"""Internal solver tools re-exported from _core.

These are exact Jacobian helpers and internal residual calculators used by the network solver.
They are not part of the public API and may change without notice.
"""

from . import _core

# Re-export all internal solver tools from _core
Pt_from_static_and_jacobian_M = _core.Pt_from_static_and_jacobian_M
Tt_from_static_and_jacobian_M = _core.Tt_from_static_and_jacobian_M
T_adiabatic_wall_and_jacobian_v = _core.T_adiabatic_wall_and_jacobian_v
adiabatic_T_complete_and_jacobian_T = _core.adiabatic_T_complete_and_jacobian_T
adiabatic_T_complete_and_jacobian_T_from_streams = (
    _core.adiabatic_T_complete_and_jacobian_T_from_streams
)
adiabatic_T_equilibrium_and_jacobians = _core.adiabatic_T_equilibrium_and_jacobians
adiabatic_T_equilibrium_and_jacobians_from_streams = (
    _core.adiabatic_T_equilibrium_and_jacobians_from_streams
)
channel_residuals_and_jacobian = _core.channel_residuals_and_jacobian
combustor_residuals_and_jacobians = _core.combustor_residuals_and_jacobians
density_and_jacobians = _core.density_and_jacobians
ejector_cd_nozzle_mass_flow_and_jacobian = _core.ejector_cd_nozzle_mass_flow_and_jacobian
ejector_cd_nozzle_exit_static_and_jacobian = _core.ejector_cd_nozzle_exit_static_and_jacobian
ejector_choked_mass_flow_and_jacobian = _core.ejector_choked_mass_flow_and_jacobian
ejector_critical_back_pressure_and_jacobian = _core.ejector_critical_back_pressure_and_jacobian
ejector_dead_head_back_pressure_and_jacobian = _core.ejector_dead_head_back_pressure_and_jacobian
ejector_element_residuals_and_jacobian = _core.ejector_element_residuals_and_jacobian
ejector_huang_critical_back_pressure_and_jacobian = (
    _core.ejector_huang_critical_back_pressure_and_jacobian
)
ejector_huang_mixing_coefficient = _core.ejector_huang_mixing_coefficient
ejector_nozzle_exit_mach = _core.ejector_nozzle_exit_mach
ejector_sudden_expansion_dead_head_and_jacobian = (
    _core.ejector_sudden_expansion_dead_head_and_jacobian
)
ejector_jetpump_incompressible_residuals_and_jacobian = (
    _core.ejector_jetpump_incompressible_residuals_and_jacobian
)
ejector_entrainment_ratio_and_jacobian = _core.ejector_entrainment_ratio_and_jacobian
ejector_jetpump_discharge_and_jacobian = _core.ejector_jetpump_discharge_and_jacobian
enthalpy_and_jacobian = _core.enthalpy_and_jacobian
friction_and_jacobian = _core.friction_and_jacobian
friction_and_jacobian_colebrook = _core.friction_and_jacobian_colebrook
friction_and_jacobian_haaland = _core.friction_and_jacobian_haaland
friction_and_jacobian_petukhov = _core.friction_and_jacobian_petukhov
friction_and_jacobian_serghides = _core.friction_and_jacobian_serghides
lossless_pressure_and_jacobian = _core.lossless_pressure_and_jacobian
mach_number_and_jacobian_v = _core.mach_number_and_jacobian_v
mixer_from_streams_and_jacobians = _core.mixer_from_streams_and_jacobians
momentum_chamber_residual_and_jacobian = _core.momentum_chamber_residual_and_jacobian
nusselt_and_jacobian_dittus_boelter = _core.nusselt_and_jacobian_dittus_boelter
nusselt_and_jacobian_gnielinski = _core.nusselt_and_jacobian_gnielinski
nusselt_and_jacobian_petukhov = _core.nusselt_and_jacobian_petukhov
nusselt_and_jacobian_sieder_tate = _core.nusselt_and_jacobian_sieder_tate
orifice_compressible_mdot_and_jacobian = _core.orifice_compressible_mdot_and_jacobian
orifice_compressible_residuals_and_jacobian = _core.orifice_compressible_residuals_and_jacobian
orifice_iso5167_residuals_and_jacobian = _core.orifice_iso5167_residuals_and_jacobian
orifice_jobson_residuals_and_jacobian = _core.orifice_jobson_residuals_and_jacobian
orifice_mdot_and_jacobian = _core.orifice_mdot_and_jacobian
orifice_residuals_and_jacobian = _core.orifice_residuals_and_jacobian
orifice_tap_pressure_and_jacobian = _core.orifice_tap_pressure_and_jacobian
approach_static_pressure = _core.approach_static_pressure
discharge_total_pressure = _core.discharge_total_pressure
sonic_static_pressure = _core.sonic_static_pressure
pressure_loss_and_jacobian = _core.pressure_loss_and_jacobian
viscosity_and_jacobians = _core.viscosity_and_jacobians
viscosity_mass_gradient = _core._viscosity_mass_gradient
wall_coupling_and_jacobian = _core.wall_coupling_and_jacobian
wall_coupling_and_jacobian_multilayer = _core.wall_coupling_and_jacobian_multilayer
wall_exchange_and_jacobian = _core.wall_exchange_and_jacobian
wall_resistance_fixed_point_derivatives = _core.wall_resistance_fixed_point_derivatives
wall_radiation_sides = _core.wall_radiation_sides
WallRadiationSides = _core.WallRadiationSides
wall_adiabatic_face = _core.wall_adiabatic_face
WallAdiabaticFace = _core.WallAdiabaticFace

__all__ = [
    "Pt_from_static_and_jacobian_M",
    "Tt_from_static_and_jacobian_M",
    "T_adiabatic_wall_and_jacobian_v",
    "adiabatic_T_complete_and_jacobian_T",
    "adiabatic_T_complete_and_jacobian_T_from_streams",
    "adiabatic_T_equilibrium_and_jacobians",
    "adiabatic_T_equilibrium_and_jacobians_from_streams",
    "channel_residuals_and_jacobian",
    "combustor_residuals_and_jacobians",
    "density_and_jacobians",
    "ejector_cd_nozzle_mass_flow_and_jacobian",
    "ejector_cd_nozzle_exit_static_and_jacobian",
    "ejector_choked_mass_flow_and_jacobian",
    "ejector_critical_back_pressure_and_jacobian",
    "ejector_dead_head_back_pressure_and_jacobian",
    "ejector_element_residuals_and_jacobian",
    "ejector_jetpump_incompressible_residuals_and_jacobian",
    "ejector_entrainment_ratio_and_jacobian",
    "ejector_jetpump_discharge_and_jacobian",
    "ejector_nozzle_exit_mach",
    "enthalpy_and_jacobian",
    "friction_and_jacobian",
    "friction_and_jacobian_colebrook",
    "friction_and_jacobian_haaland",
    "friction_and_jacobian_petukhov",
    "friction_and_jacobian_serghides",
    "lossless_pressure_and_jacobian",
    "mach_number_and_jacobian_v",
    "mixer_from_streams_and_jacobians",
    "momentum_chamber_residual_and_jacobian",
    "nusselt_and_jacobian_dittus_boelter",
    "nusselt_and_jacobian_gnielinski",
    "nusselt_and_jacobian_petukhov",
    "nusselt_and_jacobian_sieder_tate",
    "orifice_compressible_mdot_and_jacobian",
    "orifice_compressible_residuals_and_jacobian",
    "orifice_iso5167_residuals_and_jacobian",
    "orifice_jobson_residuals_and_jacobian",
    "orifice_mdot_and_jacobian",
    "orifice_residuals_and_jacobian",
    "orifice_tap_pressure_and_jacobian",
    "approach_static_pressure",
    "discharge_total_pressure",
    "sonic_static_pressure",
    "pressure_loss_and_jacobian",
    "viscosity_and_jacobians",
    "viscosity_mass_gradient",
    "wall_coupling_and_jacobian",
    "wall_coupling_and_jacobian_multilayer",
    "wall_exchange_and_jacobian",
    "wall_resistance_fixed_point_derivatives",
    "wall_radiation_sides",
    "WallRadiationSides",
    "wall_adiabatic_face",
    "WallAdiabaticFace",
]
