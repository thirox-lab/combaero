import math
from typing import Annotated, Literal, get_args

from pydantic import (
    AliasChoices,
    BaseModel,
    BeforeValidator,
    ConfigDict,
    Field,
    field_validator,
    model_validator,
)

from combaero.network.components import DEFAULT_REGIME

#: Outlets a GUI flow splitter may have (#903); its ports are out0 ... out31.
SPLITTER_MAX_OUTLETS = 32


class SmoothModelData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["smooth"] = "smooth"


class RibbedModelData(BaseModel):
    """Rib-roughened walls, on a provenanced correlation set.

    Restored in 0.7.0 after the previous rib correlation was removed for
    failing provenance review -- its friction multiplier was 4-5x below the
    only rib datum in the repository. See issue #334.
    """

    model_config = ConfigDict(extra="ignore")
    type: Literal["ribbed"] = "ribbed"
    e_D: float = 0.06
    p_e: float = 10.0
    alpha_deg: float = 90.0
    # Aspect ratio. Carried here rather than on the channel because a
    # hydraulic diameter does not determine one -- a 2:1 duct and a square
    # duct can share a Dh.
    W_H: float = 1.0
    # 1, 2 or 4. Two means two OPPOSITE walls, the configuration Han measured.
    n_ribbed_walls: int = 2
    # A user knob on the smooth walls' contribution. Plain smooth walls give
    # h_s/h_r ~ 0.42 where Han's own channel average implies 0.70, because ribs
    # enhance the adjacent smooth wall too. About 1.67 reproduces Han.
    smooth_wall_Nu_multiplier: float = 1.0


class ImpingementModelData(BaseModel):
    """Jet array impingement cooling, one row (Florschuetz, Truman and
    Metzger 1981), on a real crossflow term.

    The channel is the row's crossflow channel and its mass flow is the row's
    jet flow. The convective area defaults to the row's target footprint
    matched to the channel cross-section, A_flow * (xn/d) / (z/d), so the
    hole count follows from the channel (#462).

    Restored in 0.7.0 after the previous impingement correlation was removed
    for citing this same paper while taking no crossflow input at all,
    making it structurally unable to be that correlation. See issue #337.

    One element models one spanwise row. A full array is a chain of these,
    each with its own ``row`` -- there is no single "channel Nu" for a whole
    array the way there is for a smooth or ribbed duct, since downstream
    rows see progressively more crossflow than upstream ones.
    """

    model_config = ConfigDict(extra="ignore")
    type: Literal["impingement"] = "impingement"
    d_jet: float = 0.002
    xn_d: float = 8.0
    yn_d: float = 6.0
    z_d: float = 2.0
    # 1-indexed, counting from upstream. Row 1 sees zero crossflow by
    # definition (Gc/Gj = 0, Florschuetz's own Nu1).
    row: int = 1
    # Jet-plate discharge coefficient. 0.79 is the source's own recommended
    # default absent a measured value; combaero has no discharge-coefficient
    # correlation of its own for a jet-plate array yet (issue #375).
    C_D: float = 0.79


class SingleJetImpingementModelData(BaseModel):
    """A single free round jet (Goldstein, Behbahani and Heppelmann 1986).

    Unlike ImpingementModelData, there is no array and no crossflow. Nu is
    the AVERAGE over the disc of radius R_D * d_jet around the stagnation
    point (Han Eq. 4.1, "an average heat-transfer coefficient correlation"),
    and the channel's convective area defaults to exactly that disc (#462).
    The channel's mass flow is the jet's flow.
    """

    model_config = ConfigDict(extra="ignore")
    type: Literal["single_jet_impingement"] = "single_jet_impingement"
    bc: Literal["constant_heat_flux", "constant_wall_temperature"] = "constant_heat_flux"
    d_jet: float = 0.003
    # Jet-to-target-plate spacing / d_jet. 7.75 is the correlation's own
    # optimum spacing.
    L_D: float = 7.75
    # Radius of the averaging disc / d_jet. 5.0 matches the source's own
    # closed-form check point (Re 25,000, L/D 7.75: Nu 60 / 56).
    R_D: float = 5.0


class PinFinModelData(BaseModel):
    """A pin-fin array on provenanced correlation sets (#335).

    Restored after the previous pin-fin correlation was removed for failing
    provenance review: among other things its pressure drop was a factor of 4
    low against Metzger's own friction definition.

    Field names are the pre-0.7.0 ones, so a network saved then loads with its
    values intact; ``channel_height`` / ``pin_diameter`` is the pin H/D. The
    defaults are a working example inside the Metzger sets' validity box.
    Inline arrays use Chyu (1990)'s inline friction, the only inline source
    in hand, measured at one geometry (H/D 1, S/D = X/D = 2.5).
    """

    model_config = ConfigDict(extra="ignore")
    type: Literal["pin_fin"] = "pin_fin"
    pin_diameter: float = 0.005
    channel_height: float = 0.005
    S_D: float = 2.5
    X_D: float = 2.5
    N_rows: int = 10
    is_staggered: bool = True
    # Pin conductivity for the fin efficiency [W/(m K)]; 20 is a nickel
    # superalloy at metal temperature.
    k_pin: float = 20.0


# Dimpled surfaces were removed in 0.7.0 -- their correlation could not be
# traced to its cited source -- and remain deferred (issue #339). Saved
# networks carrying one are rejected with a message naming the reason; see
# graph_builder. Ribbed, impingement and pin-fin were restored on provenanced
# replacements (#334, #337, #335).
_SurfaceUnion = (
    SmoothModelData
    | RibbedModelData
    | ImpingementModelData
    | SingleJetImpingementModelData
    | PinFinModelData
)

#: Surface types that existed and were removed, and why (see map_surface_model).
REMOVED_SURFACE_TYPES = frozenset({"dimpled"})


def removed_surface_message(surface_type: object) -> str:
    return (
        f"Surface type {surface_type!r} was removed in 0.7.0: its correlation "
        "could not be traced to its cited source. See issue #339. Use a "
        "smooth surface."
    )


def _refuse_unknown_surface(v: object) -> object:
    # Without this a removed or unknown type reaches the union, which refuses
    # it as one literal mismatch per member -- five errors, none saying the
    # type was removed or what to pick.
    if isinstance(v, dict) and "type" in v:
        known = {get_args(m.model_fields["type"].annotation)[0] for m in get_args(_SurfaceUnion)}
        t = v["type"]
        if t in REMOVED_SURFACE_TYPES:
            raise ValueError(removed_surface_message(t))
        if t not in known:
            raise ValueError(f"unknown surface type {t!r}; expected one of {sorted(known)}")
    return v


SurfaceModelData = Annotated[_SurfaceUnion, BeforeValidator(_refuse_unknown_surface)]


# --- Node Data Definitions ---


class PlenumData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class CompositionData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    mode: Literal["mole", "mass"] = "mole"
    source: Literal["dry_air", "humid_air", "fuel", "custom"] = "humid_air"
    custom_fractions: dict[str, float] | None = None
    relative_humidity: float | None = None  # None = use global or 0.6
    ambient_T: float | None = None  # None = use global or 288.15
    ambient_P: float | None = None  # None = use global or 101325.0


class MassBoundaryData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    m_dot: float = 1.0
    Tt: float = 300.0
    composition: CompositionData = Field(default_factory=CompositionData)
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class FuelBoundaryData(BaseModel):
    """A fuel controller (#892): the fuel flow that puts ``combustor_id``'s
    exit at ``Tt_target``, on the ``branch`` asked for."""

    model_config = ConfigDict(extra="ignore")
    combustor_id: str | None = None
    Tt_target: float = 1500.0
    branch: Literal["lean", "rich"] = "lean"
    Tt: float = 300.0
    composition: CompositionData = Field(
        default_factory=lambda: CompositionData(source="fuel", mode="mole")
    )
    m_dot_guess: float | None = None
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class PressureBoundaryData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    Pt: float = 101325.0
    Tt: float = 300.0
    composition: CompositionData = Field(default_factory=CompositionData)
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class WallData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    label: str | None = None
    initial_guess: dict[str, float] = Field(default_factory=dict)


class HeatSourceData(BaseModel):
    """A heat source/sink wired by one edge to a plenum, momentum chamber or
    combustor (#792): an EnergyBoundary on that node. ``mode`` picks which
    value is set; the others are kept (a mode switch loses nothing) but not
    used. Defaults are a working example: 1 kW into the node."""

    model_config = ConfigDict(extra="ignore")
    mode: Literal["Q", "fraction", "Tt_target"] = "Q"
    Q: float = 1000.0  # [W], > 0 heats the fluid
    fraction: float = 0.01  # [-] of the inflow's sensible enthalpy
    Tt_target: float = 600.0  # [K], the node's TOTAL temperature
    label: str | None = None

    @field_validator("Tt_target")
    @classmethod
    def _positive_target(cls, v: float) -> float:
        if not (math.isfinite(v) and v > 0.0):
            raise ValueError("Tt_target must be a positive total temperature [K]")
        return v


class ConstantFractionLossData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["constant_fraction"] = "constant_fraction"
    xi: float = 0.03


class LinearThetaFractionLossData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["linear_theta_fraction"] = "linear_theta_fraction"
    k: float = 0.5
    xi0: float = 0.02


class ConstantHeadLossData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["constant_head"] = "constant_head"
    zeta: float = 5.0


class LinearThetaHeadLossData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["linear_theta_head"] = "linear_theta_head"
    k: float = 1.0
    zeta0: float = 3.0


PressureLossData = (
    ConstantFractionLossData
    | LinearThetaFractionLossData
    | ConstantHeadLossData
    | LinearThetaHeadLossData
)


class CombustorData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    method: Literal["complete", "equilibrium"] = "complete"
    area: float | None = None  # None = derive from Dh (pi/4 * Dh^2)
    Dh: float | None = None
    surface: SurfaceModelData = Field(default_factory=SmoothModelData)
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class MomentumChamberData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    area: float | None = None  # None = derive from Dh (pi/4 * Dh^2)
    Dh: float | None = None  # None = inherit from upstream channel (Dh = D)
    surface: SurfaceModelData = Field(default_factory=SmoothModelData)
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


# --- Element Data Definitions ---


class ChannelData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    L: float = 1.0
    D: float | None = None  # None = inherit from upstream element geometry
    Dh: float | None = None  # None = circular (Dh = D); override for non-circular ducts
    roughness: float = 1e-5
    friction_model: Literal["haaland", "serghides", "colebrook", "petukhov"] = "haaland"
    surface: SurfaceModelData = Field(default_factory=SmoothModelData)
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0
    regime: Literal["default", "incompressible", "compressible"] = "default"
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class OrificeData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    diameter: float | None = None  # None = the default bore (a bore never inherits)
    # Read for correlation "fixed" only (graph_builder._fixed_cd, #872).
    Cd: float = 0.6
    # Normed metering (ISO 5167, function of beta) and discharge-hole
    # (Idelchik / McGreehan-Schotsch, function of L/d and r/d) families. They
    # are not interchangeable; there is deliberately no "Auto" member.
    correlation: Literal[
        "ReaderHarrisGallagher",
        "IdelchikThick",
        "IdelchikBeveled",
        "IdelchikRounded",
        "McGreehanSchotsch",
        "Lichtarowicz",
        "fixed",
    ] = "ReaderHarrisGallagher"
    plate_thickness: float = 0.0  # t [m], the hole length for a wall orifice
    edge_radius: float = 0.0  # r [m]
    bevel_depth: float = 0.0  # bevel depth along the hole axis [m]
    # The pressure taps whose differential is REPORTED (dP_tap). Measurement
    # only: the solved flow never depends on it (#649/#645).
    tap: Literal["Corner", "Flange", "D_D2"] = "Flange"
    # False: no pipe diameter -- a lumped Cd*A restriction (a hole in a wall,
    # beta = 0) whatever channels adjoin it (#788). True: the adjoining pipe
    # sets beta; a bore equal to it is no orifice, a larger one is refused.
    in_pipe: bool = True
    regime: Literal["default", "incompressible", "compressible"] = "default"
    # How the compressible regime expands the flow (#854). "default" lets the
    # element select at build (#861, 0.8.0): jobson1955 for an
    # incompressible-referenced Cd on a sharp-edged bore at beta <= 0.75,
    # else the nozzle. A file saved without the key gets that selection; an
    # explicit "nozzle", "jobson1955" or "iso5167" (RHG only) is kept.
    compressible_model: Literal["default", "nozzle", "jobson1955", "iso5167"] = "default"
    # Jobson's velocity coefficient, Cd = Cv Ci; read for jobson1955 only.
    Cv: float = Field(default=1.0, gt=0.0, le=1.0)
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None

    @model_validator(mode="before")
    @classmethod
    def refuse_removed_correlations(cls, data: dict) -> dict:
        # A file saved before the removal names a correlation that no longer
        # exists. Say which one and what to pick, rather than a bare Literal
        # mismatch; silently swapping in another Cd would change the result.
        if isinstance(data, dict) and data.get("correlation") == "Miller":
            raise ValueError(
                "The 'Miller' orifice correlation was removed: no published "
                "source matched its formula. Pick 'ReaderHarrisGallagher' for "
                "this orifice (Stolz, which 'Miller' had computed, was removed "
                "in 0.9.0, #876)."
            )
        if isinstance(data, dict) and data.get("correlation") == "Stolz":
            raise ValueError(
                "The 'Stolz' orifice correlation (ISO 5167:1980) was removed in "
                "0.9.0 (#876): ISO replaced it with Reader-Harris/Gallagher. Pick "
                "'ReaderHarrisGallagher' for this orifice."
            )
        return data

    @model_validator(mode="before")
    @classmethod
    def migrate_area_to_diameter(cls, data: dict) -> dict:
        if isinstance(data, dict) and "area" in data and "diameter" not in data:
            import math

            area_val = data.pop("area")
            # 0 is a closed restriction (#1024), as a 0 diameter is.
            if area_val is not None and area_val >= 0:
                data["diameter"] = math.sqrt(4.0 * area_val / math.pi)
        return data


class EffusionLinerData(BaseModel):
    """Effusion liner (#471): a backside coolant duct whose wall bleeds
    through effusion holes into the gas. Expanded into N stations: a Bassett
    bleed plus a duct-fed McGreehan-Schotsch panel each.

    Defaults: Andrews et al.'s plate C holes and pitch (3.27 mm, 6.3 mm wall,
    15.24 mm), a hole t/d inside Rohde's crossflow-tested 0.51-4, over a 30 mm
    backside duct. The crossflow ratio U1/Vi is set mainly by the duct's own
    pressure drop against the holes' drive (roughly sqrt(dp_duct/dp_hole)),
    not by the duct size: a pressure-driven duct carries more flow when made
    bigger. With 200 Pa along the duct and ~3.8 kPa across the holes it peaks
    at 0.27, inside the Rohde-scored band.
    """

    model_config = ConfigDict(extra="ignore")
    n_segments: int = Field(default=4, ge=1)
    length: float = Field(default=0.152, gt=0.0)
    width: float = Field(default=0.152, gt=0.0)
    duct_height: float = Field(default=0.03, gt=0.0)
    hole_diameter: float = Field(default=3.27e-3, gt=0.0)
    wall_thickness: float = Field(default=6.3e-3, gt=0.0)
    pitch_x: float = Field(default=15.24e-3, gt=0.0)
    pitch_y: float = Field(default=15.24e-3, gt=0.0)
    angle_deg: float = Field(default=90.0, gt=0.0, le=90.0)
    entry_K: float = Field(default=0.5, ge=0.0)
    roughness: float = Field(default=0.0, ge=0.0)
    wall_conductivity: float = Field(default=20.0, gt=0.0)
    gas_film: Literal["none", "baldauf_sellers"] = "none"
    gas_augmentation: float = Field(default=1.0, gt=0.0)
    turbulence_intensity: float = Field(default=0.05, gt=0.0)
    gas_heat_flux: float = 0.0
    # Which way the hot gas runs past the panels, for the film the liner lays
    # on a wall downstream (#695).
    film_direction: Literal["with_coolant", "against_coolant"] = "with_coolant"
    # The plates' alloy, for their wall temperature limit only (#998); None
    # is the generic metal cap.
    material: str | None = None
    label: str | None = None


class EffusionPlateData(BaseModel):
    """Effusion / film plate, plenum-fed (#471): it owns its wall.

    The node it discharges into decides the gas side: a momentum chamber
    (flow) gives the gas coefficient, temperature and blowing ratio; a plenum
    (state, no flow) takes the imposed gas_heat_flux. Defaults are Andrews et
    al.'s plate C (86-GT-225 / 88-GT-290), the geometry the closure is scored
    on: 3.27 mm holes, 6.3 mm wall, 15.24 mm square pitch, 152 mm panel.
    """

    model_config = ConfigDict(extra="ignore")
    hole_diameter: float = Field(default=3.27e-3, gt=0.0)
    wall_thickness: float = Field(default=6.3e-3, gt=0.0)
    pitch_x: float = Field(default=15.24e-3, gt=0.0)
    pitch_y: float = Field(default=15.24e-3, gt=0.0)
    panel_length: float = Field(default=0.152, gt=0.0)
    panel_width: float = Field(default=0.152, gt=0.0)
    angle_deg: float = Field(default=90.0, gt=0.0, le=90.0)
    correlation: Literal["IdelchikThick", "Lichtarowicz", "McGreehanSchotsch", "fixed"] = (
        "IdelchikThick"
    )
    Cd: float = 0.6
    internal_Nu_multiplier: float = Field(default=1.0, gt=0.0)
    wall_conductivity: float = Field(default=20.0, gt=0.0)
    gas_film: Literal["none", "baldauf_sellers"] = "none"
    gas_augmentation: float = Field(default=1.0, gt=0.0)
    # What gas_augmentation multiplies: the chamber's own h (default) or a
    # turbulent flat plate restarted gas_start_length [m] upstream of the
    # panel (a bleed, a starter-film lip, a step). The length is read only
    # for the flat plate, so switching back keeps it without refusing.
    gas_baseline: Literal["chamber", "flat_plate"] = "chamber"
    gas_start_length: float | None = Field(default=None, ge=0.0)
    # Opt-in hole-angle Cd term (Idelchik diagram 3-2 on IdelchikThick's
    # inlet coefficient); read only with IdelchikThick, like Cd with 'fixed'.
    cd_inclination: Literal["none", "idelchik_3_2"] = "none"
    turbulence_intensity: float = Field(default=0.05, gt=0.0)
    gas_heat_flux: float = 0.0
    # The plate's alloy, for its wall temperature limit only (#998): the hot
    # end of its k(T) envelope. None is the generic metal cap.
    material: str | None = None
    label: str | None = None


class SlotFilmData(BaseModel):
    """A 2-D film slot (#714): coolant from the upstream node into the gas
    through a continuous slot, and a film PROVIDER. It owns no wall: a
    thermal wall on the gas duct downstream names it (``film``), the span
    measured from the slot lip.

    Defaults are Marek & Tacina's D-7958 rig, on which the film is scored:
    a 0.63 cm slot across the 37 cm duct, Tu 0.14 (14%, the case Eq. 5
    fits best). ``turbulence_intensity`` is a FRACTION, as on every other
    element (#773). The 1 mm lip is not from the report; it only sets the orifice's
    length.
    """

    model_config = ConfigDict(extra="ignore")
    slot_height: float = Field(default=6.3e-3, gt=0.0)
    span: float = Field(default=0.37, gt=0.0)
    lip_thickness: float = Field(default=1.0e-3, gt=0.0)
    turbulence_intensity: float = Field(default=0.14, ge=0.0, le=1.0)
    C_m_multiplier: float = Field(default=1.0, gt=0.0)
    correlation: Literal["IdelchikThick", "Lichtarowicz", "McGreehanSchotsch", "fixed"] = (
        "IdelchikThick"
    )
    Cd: float = Field(default=0.6, gt=0.0, le=1.0)
    # Film-side h on the wall the film covers (#962): "none" keeps the duct's
    # own h; "seban1960" (recommended) or the book's wall jet.
    film_htc: Literal["none", "seban1960", "ballal_lefebvre1973"] = "none"
    label: str | None = None


class FilmCoolingData(BaseModel):
    """Film-cooling holes (#386): coolant from the upstream node into the gas,
    and a film PROVIDER. It owns no wall: a thermal wall on the gas duct
    downstream names it (``film``) and that side's driving temperature
    becomes T_aw - eta (T_aw - T_coolant).

    Defaults are Colban et al.'s Eq. (19) at Gritsch et al.'s configuration
    B, the geometry it is scored on (#686): laidback fan-shaped, AR 3.5,
    P/D 6, t/P 0.49, 30 deg, Tu 4%. The 1 mm metering diameter only sets the
    scale; the correlation is in P/D, t/P and x/D.
    """

    model_config = ConfigDict(extra="ignore")
    hole_type: Literal["shaped", "cylindrical"] = "shaped"
    hole_diameter: float = Field(default=1.0e-3, gt=0.0)
    holes_per_row: int = Field(default=20, ge=1)
    lateral_pitch: float = Field(default=6.0e-3, gt=0.0)
    n_rows: int = Field(default=1, ge=1)
    row_pitch: float | None = Field(default=None, gt=0.0)
    wall_thickness: float = Field(default=3.0e-3, gt=0.0)
    angle_deg: float = Field(default=30.0, gt=0.0, le=90.0)
    area_ratio: float = Field(default=3.5, ge=1.0)
    breakout_width: float = Field(default=2.94e-3, gt=0.0)
    turbulence_intensity: float = Field(default=0.04, gt=0.0)
    correlation: Literal["IdelchikThick", "Lichtarowicz", "McGreehanSchotsch", "fixed"] = (
        "IdelchikThick"
    )
    Cd: float = 0.6
    # Opt-in hole-angle Cd term, read only with IdelchikThick (see
    # EffusionPlateData).
    cd_inclination: Literal["none", "idelchik_3_2"] = "none"
    label: str | None = None


class ImpingementArrayData(BaseModel):
    """Impingement array, channel exit (#465): ONE configuration.

    A jet plate fed from the upstream node, impinging on a target, the spent
    air leaving down the gap to the downstream node (Florschuetz, Truman and
    Metzger 1981). The builder expands it into the row chain. The defaults are
    a Florschuetz test geometry (0.1 in holes, (5, 4, 2), 10 rows) and sit
    inside the correlation's range at a few kPa across the array.
    """

    model_config = ConfigDict(extra="ignore")
    n_rows: int = Field(default=10, ge=1)
    d_jet: float = Field(default=0.00254, gt=0.0)
    xn_d: float = Field(default=5.0, gt=0.0)
    yn_d: float = Field(default=4.0, gt=0.0)
    z_d: float = Field(default=2.0, gt=0.0)
    span: float = Field(default=0.122, gt=0.0)
    plate_thickness: float = Field(default=0.00254, gt=0.0)
    pattern: Literal["inline", "staggered"] = "inline"
    # Discharge-hole family only: a jet plate has no pipe, so no beta.
    correlation: Literal["fixed", "IdelchikThick", "Lichtarowicz", "McGreehanSchotsch"] = "fixed"
    Cd: float = 0.79  # Florschuetz's own default, cb.FLORSCHUETZ_1981_DEFAULT_CD
    # The holes' Cd under the gap's crossflow (#800), ImpingementPlateElement's.
    cd_crossflow: Literal["rogers_hersh_1976", "florschuetz_metzger_su_1982", "none"] = (
        "rogers_hersh_1976"
    )
    Nu_multiplier: float = Field(default=1.0, gt=0.0)
    roughness: float = Field(default=0.0, ge=0.0)
    # Initial (bypass) crossflow (#467), used only when the bypass handle is
    # wired: the lead-in to row 1 (None = xn/2, where CR-3630's array begins)
    # and the entry loss on the gap's dynamic head (0 along a duct).
    bypass_length: float | None = Field(default=None, gt=0.0)
    bypass_entry_K: float = Field(default=0.0, ge=0.0)
    label: str | None = None


class AreaChangeData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    model_type: Literal["sharp", "conical"] = "sharp"
    F0: float | None = None  # None = inherit from upstream channel (A = pi/4*D^2)
    F1: float | None = None  # None = inherit from downstream channel
    length: float | None = None
    D_h: float = 0.0
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class LosslessConnectionData(BaseModel):
    pass


class TurbomachineData(BaseModel):
    """A compressor or turbine at a fixed total-to-total pressure ratio (#893):
    Pt_out / Pt_in = pressure_ratio (compressor) or its inverse (turbine)."""

    model_config = ConfigDict(extra="ignore")
    pressure_ratio: float = Field(4.0, ge=1.0, allow_inf_nan=False)  # finite (#1024)
    efficiency_type: Literal["polytropic", "isentropic"] = "polytropic"
    efficiency: float = Field(0.88, gt=0.0, le=1.0)
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class FlowSplitterData(BaseModel):
    """A flow splitter (#903): one inlet ("port-in-target") to N outlets
    ("port-out<i>-source"), outlet i at fractions[i] of the throughput.

    At most one fraction is None: the remainder, 1 - the others. The rest
    are in (0, 1) and, without a remainder, sum to 1 within
    SPLITTER_FRACTION_SUM_TOL; FlowSplitter refuses anything else, never
    renormalising. ``reference`` is the index of the outlet whose valve is
    open (Pt_out = Pt_in) to start from: the remainder when None, else the
    first guess. With ``auto_reference`` (default) the solve moves it to the
    outlet needing the most pressure when the start is infeasible. ``areas[i]`` is outlet
    i's valve reference area [m^2] for the reported Cd and zeta (None: the
    effective area CdA only).
    """

    model_config = ConfigDict(extra="ignore")
    fractions: list[float | None] = Field(
        default_factory=lambda: [0.5, None], min_length=2, max_length=SPLITTER_MAX_OUTLETS
    )
    reference: int | None = Field(default=None, ge=0)
    areas: list[float | None] | None = None
    auto_reference: bool = True
    label: str | None = None

    @model_validator(mode="after")
    def _indices_fit(self) -> "FlowSplitterData":
        n = len(self.fractions)
        if self.reference is not None and self.reference >= n:
            raise ValueError(f"reference {self.reference} is not an outlet (0 to {n - 1})")
        if self.areas is not None and len(self.areas) != n:
            raise ValueError(f"areas has {len(self.areas)} entries for {n} outlets")
        return self


class MPCETeeData(BaseModel):
    """Mynard-based multi-port chamber tee element with constrained topology.

    3-port junction using Mynard's Unified0D pseudosupplier residual. The
    ``flow_direction`` field declares whether this junction is a merge
    (joining flow, 2 inlets + 1 outlet) or a branch (separating flow,
    1 inlet + 2 outlets); the runtime element raises if the observed
    flow pattern disagrees.
    """

    model_config = ConfigDict(extra="ignore")
    label: str | None = None
    flow_direction: Literal["merge", "branch"] = "branch"
    # Junction closure model. "mynard" (default) is the full Unified0D
    # pseudosupplier physics; "constant_k" is the simplest-model tier:
    # fixed handbook loss coefficients referenced to the common-leg
    # dynamic head (Idelchik convention), K independent of the split;
    # "ideal_momentum" is the parameter-free baseline (#517): one junction
    # static pressure, momentum balance on every receiving port, no tuned
    # loss. Its flow_direction sets the port wiring only -- roles follow the
    # flow.
    junction_model: Literal["mynard", "constant_k", "ideal_momentum"] = "mynard"
    # constant_k only: arm loss coefficients [-], referenced to the common-leg
    # head. None = Bassett 2001's analytical pair at design_branch_share, the
    # resolved areas and theta (ConstantKTeeElement.bassett_design_K): a
    # constant K belongs to one split, and the branch's grows like (q psi)^2.
    # The builder refuses a pair, typed or default, whose mass-weighted loss at
    # design_branch_share is negative (ConstantKTeeElement.design_net_loss, #779).
    K_straight: float | None = None
    K_branch: float | None = None
    design_branch_share: float = Field(default=0.5, gt=0.0, lt=1.0)
    # ideal_momentum only (#526): receiver loss factor f. 1.0 is the
    # parameter-free closure. The floor is the robustness sweep's: down to
    # 0.25 convergence inside Mach 0.3 stays >= 98.4% (99.5% at 1.0), at 0 it
    # is 64%, and a split below f ~ M^2/4 is refused by the energy check. No
    # ceiling: f > 1 is only more loss than the ideal.
    loss_factor: float = Field(default=1.0, ge=0.25)
    theta_deg: float = 90.0  # branch angle [deg], converted to radians in graph_builder
    F_C: float | None = None  # None = inherit from common/straight-arm channel
    F_branch: float | None = None  # None = inherit from branch-arm channel
    psi: float = 1.0  # fallback ratio used only when F_branch is None and not inherited
    # Joining-side etransfer correction scale (combaero extension to faithful
    # Mynard). Default 0.2 calibrated against Bassett K11/K12 analytical +
    # Idelchik 1966 tabulated values. Set to 0.0 to disable the correction
    # (faithful Mynard); set to a custom value for re-calibrated workflows.
    joining_etransfer_alpha: float | None = None
    initial_guess: dict[str, float] = Field(default_factory=dict)

    @field_validator("theta_deg")
    @classmethod
    def _branch_angle_in_range(cls, v: float) -> float:
        # The rule the frontend has enforced since #141. The closures read the
        # angle only through cos/sin, so an unbounded value is not refused
        # there: 900 silently becomes a 180 deg tee, and 0 lays the branch on
        # the straight run. A negative angle is the mirror image.
        if not 0.0 < abs(v) <= 90.0:
            raise ValueError(f"tee branch angle must satisfy 0 < |theta| <= 90 deg, got {v}")
        return v


class EjectorData(BaseModel):
    """Critical-mode supersonic ejector (3-port: primary + secondary inlets,
    one outlet). Wraps ``combaero.network.ejector_element.EjectorElement``.

    Geometry is three absolute areas (the ratios are computed by the
    element). Validation (throat > 0, nozzle_exit > throat, mixing >
    nozzle_exit, 0 < recovery_efficiency <= 1) is enforced by the element
    constructor at build time; the frontend mirrors it for early feedback.
    The recovery bound is also declared on the field, so a request is refused
    by name before the network is built.
    """

    model_config = ConfigDict(extra="ignore")
    label: str | None = None
    throat_area: float = 3.14e-5  # primary nozzle throat area A_t [m^2]
    nozzle_exit_area: float = 1.0e-4  # primary nozzle exit area A_p1 [m^2], > throat
    # Constant-area mixing section A_3 [m^2], > nozzle_exit. None = inherit from
    # the outlet pipe, else EjectorElement.DEFAULT_EJECTOR_MIXING_AREA (#559).
    mixing_area: float | None = None
    # Multiplies the lossless mixed stagnation pressure to give the critical
    # back pressure P_c*. 1.0 = no artificial loss (not fitted to data); above
    # 1 is a pressure gain beyond lossless mixing (#769).
    recovery_efficiency: float = Field(default=1.0, gt=0.0, le=1.0)
    # The (P_c*, P_b0) pair ending the plateau and the droop (compressible):
    # Kracik & Dvorak's lossless mixing, or Huang 1999's printed chain with the
    # sudden-expansion dead head (validation/ejector/README.md).
    anchor_closure: Literal["kracik_dvorak", "huang1999"] = "huang1999"
    # Incompressible regime only (#588): Sanger's loss coefficients. None on the
    # nozzle / suction = the compressible nozzle efficiencies as K = 1/eta - 1.
    k_nozzle: float | None = Field(default=None, ge=0.0)
    k_suction: float | None = Field(default=None, ge=0.0)
    k_throat: float | None = Field(default=None, ge=0.0)
    k_diffuser: float | None = Field(default=None, ge=0.0)
    initial_guess: dict[str, float] = Field(default_factory=dict)

    @field_validator("throat_area", "nozzle_exit_area", mode="before")
    @classmethod
    def _unset_area_is_the_default(cls, v: object, info: object) -> object:
        # The Inspector shows an unset throat / nozzle exit as its default
        # (#559); a cleared field arrives as null and must mean the same here.
        if v is None:
            return cls.model_fields[info.field_name].default  # type: ignore[attr-defined]
        return v


class VortexData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    label: str | None = None
    # A free vortex (#661): Pt and Tt conserved, r * V_theta carried from
    # r_in to r_out. Files saved before #661 also carry r_c and n (the old
    # Vatistas core); extra="ignore" drops them, and the GUI notes the change.
    r_out: float = 0.10  # outlet radius [m]
    r_in: float = 0.05  # inlet radius [m]; 0 = on-axis, no angular momentum
    # Shaft speed [rpm]; None = use the global solver setting. Named omega_rpm
    # before #773 (omega is 2 pi f, not a speed in rpm): files saved with the
    # old name still load.
    shaft_speed_rpm: float | None = Field(
        default=None, validation_alias=AliasChoices("shaft_speed_rpm", "omega_rpm")
    )
    v_theta_in: float | None = None  # inlet swirl [m/s]; None = solid body, omega * r_in
    initial_guess: dict[str, float] = Field(default_factory=dict)


class DiscreteLossData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    correlation_type: Literal[
        "constant_fraction", "constant_head", "linear_theta_fraction", "linear_theta_head"
    ] = "constant_fraction"
    xi: float = 0.03
    k: float = 0.001
    xi0: float = 0.02
    #: Design mass flow [kg/s] of a fraction loss (#550): below it the loss
    #: scales with the dynamic head and vanishes at zero flow, so the flow can
    #: run low or reverse. None keeps a constant fraction at every flow.
    m_ref: float | None = None
    zeta: float = 1.0
    zeta0: float = 1.0
    area: float | None = None
    theta_source: str | None = None
    surface: SurfaceModelData = Field(default_factory=SmoothModelData)
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0
    label: str | None = None
    initial_guess: dict[str, float] = Field(default_factory=dict)


class WallLayerData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    thickness: float = 0.003
    conductivity: float = 20.0
    material: str = "generic"


class ThermalWallData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["thermal"] = "thermal"
    # Support both single-layer (legacy) and multi-layer
    thickness: float | None = None
    conductivity: float | None = None
    layers: list[WallLayerData] | None = None
    area: float | None = 0.05
    # Per-side areas (#887): the wall's two faces [m^2], both or neither;
    # set, they replace `area`, and `wall_area` is the layers' conduction area.
    area_a: float | None = Field(default=None, gt=0.0)
    area_b: float | None = Field(default=None, gt=0.0)
    wall_area: Literal["log_mean", "slab_a", "slab_b"] = "log_mean"
    R_fouling: float = 0.0
    # How the two elements lie along the wall, from_node -> to_node: with the
    # flow signs it decides parallel or counter flow (ThermalWall, #627).
    orientation: Literal["aligned", "opposed"] = "aligned"
    # A film covering one side (#386): the id of a film_cooling, slot_film or
    # effusion_plate node, the side it covers, and the wall's span [m]
    # downstream of that provider's last row.
    film: str | None = None
    film_side: Literal["a", "b"] = "a"
    film_start: float = Field(default=0.0, ge=0.0)
    film_end: float | None = Field(default=None, ge=0.0)

    @model_validator(mode="before")
    @classmethod
    def migrate_legacy_fields(cls, data: dict) -> dict:
        if isinstance(data, dict) and ("layers" not in data or data["layers"] is None):
            # If we have legacy fields but no layers, migrate them
            t = data.get("thickness", 0.003)
            k = data.get("conductivity", 20.0)
            # Ensure they are not None
            t = t if t is not None else 0.003
            k = k if k is not None else 20.0
            data["layers"] = [{"thickness": t, "conductivity": k}]
        return data


# Max solve time [s] when a request names none, for the GUI and NetworkRunner
# alike. The frontend store's solverSettings.timeout mirrors it; kept equal by
# python/tests/test_solve_timeout_default_sync.py. None in a request still
# means no limit.
DEFAULT_SOLVE_TIMEOUT = 30.0


class SolverSettings(BaseModel):
    model_config = ConfigDict(extra="ignore")
    # The library-wide default (combaero.network.components.DEFAULT_REGIME),
    # which the frontend store mirrors; test_regime_default_sync.py keeps the
    # three equal. A payload without solver settings -- headless runs, tests,
    # older clients -- gets the same physics as a freshly drawn GUI network.
    global_regime: Literal["incompressible", "compressible"] = DEFAULT_REGIME
    init_strategy: Literal[
        "default",
        "analytical_pt_prop",
        "outletref_warmstart",
        "homotopy",
        "continuation",
    ] = "default"

    @field_validator("init_strategy", mode="before")
    @classmethod
    def _migrate_init_strategy(cls, v: object) -> object:
        # The inlet-referenced incompressible warm start was removed from the
        # library in 0.8.0. A file saved with it gets its outlet-referenced
        # successor, as the frontend does on load.
        return "outletref_warmstart" if v == "incompressible_warmstart" else v

    method: Literal[
        "hybr",
        "lm",
        "broyden1",
        "broyden2",
        "anderson",
        "linearmixing",
        "diagbroyden",
        "excitingmixing",
        "krylov",
        "df-sane",
    ] = "hybr"
    timeout: float | None = DEFAULT_SOLVE_TIMEOUT
    # NetworkSolver.solve(stall_handover=...) (#864): False lets every cold
    # attempt grind to its share of the budget instead of handing over on a
    # stall. A file saved without it gets the default.
    stall_handover: bool = True
    # Global shaft speed [rpm]; None = disabled. Was omega_rpm before #773;
    # the old name still loads.
    shaft_speed_rpm: float | None = Field(
        default=None, validation_alias=AliasChoices("shaft_speed_rpm", "omega_rpm")
    )
    Nu_multiplier: float | None = (
        None  # global heat-transfer scale; stacks multiplicatively with per-element
    )
    f_multiplier: float | None = (
        None  # global friction scale; stacks multiplicatively with per-element
    )
    # Mach above which an element on an incompressible closure is reported
    # outside its validity (#588, report only). None = the library default,
    # combaero.network.compressibility.INCOMPRESSIBLE_MACH_LIMIT (0.2).
    incompressible_mach_limit: float | None = Field(default=None, gt=0.0)
    ambient_T: float | None = None  # ISO 2314: 288.15 K; overrides per-node value when humid_air
    ambient_P: float | None = None  # ISO 2314: 101325 Pa; overrides per-node value when humid_air
    ambient_RH: float | None = (
        None  # ISO 2314: 0.6; overrides per-node relative_humidity when humid_air
    )


# --- React Flow Wrapper Schemas ---


class NodePosition(BaseModel):
    model_config = ConfigDict(extra="ignore")
    x: float
    y: float


class ReactFlowNode(BaseModel):
    id: str
    type: str
    position: NodePosition
    data: dict


class ReactFlowEdge(BaseModel):
    id: str
    source: str
    target: str
    sourceHandle: str | None = None
    targetHandle: str | None = None
    type: str | None = None
    data: dict | None = None


class NetworkGraphSchema(BaseModel):
    nodes: list[ReactFlowNode]
    edges: list[ReactFlowEdge]
    solver_settings: SolverSettings = Field(default_factory=SolverSettings)

    @model_validator(mode="before")
    @classmethod
    def _normalise_solver_settings_key(cls, data: object) -> object:
        if isinstance(data, dict) and "solverSettings" in data and "solver_settings" not in data:
            data = dict(data)
            data["solver_settings"] = data.pop("solverSettings")
        return data


# --- Results Schemas ---


class StateResult(BaseModel):
    model_config = ConfigDict(extra="allow")
    T: float
    P: float
    Pt: float | None = None
    Tt: float | None = None
    rho: float | None = None
    h: float | None = None
    s: float | None = None
    mach: float | None = None
    Y: list[float]
    X: list[float] | None = None


class NodeResult(BaseModel):
    model_config = ConfigDict(extra="allow")
    state: StateResult
    success: bool = True


class ElementResult(BaseModel):
    model_config = ConfigDict(extra="allow")
    m_dot: float = 0.0
    success: bool = True


class ConvergencePoint(BaseModel):
    eval: int
    t: float
    norm: float


class WorstResidual(BaseModel):
    name: str
    residual: float


class SecondLawFlag(BaseModel):
    """An element whose control volume destroys entropy (#545), on the canvas."""

    element_id: str
    owner_id: str
    sigma_norm: float
    sigma_norm_high: float


class IncompressibleFlag(BaseModel):
    """An element on an incompressible closure above the Mach limit (#588)."""

    element_id: str
    owner_id: str
    mach: float
    limit: float


class SolveResponse(BaseModel):
    success: bool
    #: The two questions ``success`` conflates (NetworkSolver's
    #: ``__converged__``/``__consistent__``/``__outcome__``): a refused root
    #: (outcome "inconsistent") converged, a failed solve did not.
    converged: bool = False
    consistent: bool | None = None
    outcome: str = ""
    message: str = ""
    final_norm: float | None = None
    node_results: dict[str, NodeResult] = {}
    element_results: dict[str, ElementResult] = {}
    edge_results: dict[str, dict] = {}
    convergence_history: list[ConvergencePoint] = []
    worst_residuals: list[WorstResidual] = []
    solver_settings_used: dict = {}
    lm_started_at_eval: int | None = None
    #: Elements that fail the second-law check (#545), each
    #: attributed to the canvas node that owns it.
    second_law_failures: list[SecondLawFlag] = []
    #: Report only (#588): elements on an incompressible closure running above
    #: the solver's Mach limit, each attributed to its canvas node.
    incompressible_warnings: list[IncompressibleFlag] = []
    #: Mass and energy over the network's control volume (empty if the root
    #: finder did not converge).
    ledger: dict = {}
