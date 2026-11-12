"""The hole discharge-coefficient chain: one composition point for the Cd of
every hole element (OrificeElement and its plate, film, slot and impingement
subclasses).

A hole's Cd composes, in this order:

  1. the BASE Cd: the 'fixed' value, a normed metering correlation (ISO
     5167), or a discharge-hole correlation (Idelchik, McGreehan-Schotsch,
     Lichtarowicz) at the hole Reynolds number, with
       - the INCLINATION term (Idelchik diagram 3-2, CD_INCLINATION_MODELS)
         inside the hole geometry the correlation sees, and
       - the INLET (supply-side) crossflow U1/Vi of a duct-fed hole, read by
         McGreehan-Schotsch only;
  2. the EXIT crossflow term (Rowbury et al. 2001, CD_EXIT_CROSSFLOW_MODELS),
     an additive outlet loss on the base Cd;
  3. the JET crossflow law of an impingement row (_CD_CROSSFLOW_LAWS:
     Rogers & Hersh 1976 or Florschuetz, Metzger & Su 1982), which replaces
     the hole's Cd where it acts.

HoleCd is the family's public spec: the model choices (correlation, fixed
Cd, and each term's model) travel as one frozen object, ``hole_cd=HoleCd(...)``
on every hole element and builder, like GasRadiation for radiation. The
element keeps the geometry (bore, angle, hole type, the exit term's region and
orientation); a field left at None/"auto" is that element's default.

_HoleCdChain holds the resolved options, gathers the chain's inputs from the
network state, and maps the partials onto unknown names; the composition
itself -- every term's value and the chain rule between them -- is one C++
call, hole_cd_compose (include/hole_cd.h), over the kernels of orifice.h and
crossflow_cd.h. Every hole element holds one (``_hole_cd``) and delegates to
it. The element supplies the geometry and topology the chain reads (its bore,
hole count and Reynolds number, its nodes, its duct or crossflow sources);
that topology is why the inputs and the mapping stay here.
"""

import dataclasses
import math
import warnings
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import combaero as cb
from combaero import _solver_tools

if TYPE_CHECKING:
    from .components import NetworkMixtureState, OrificeElement
    from .graph import FlowNetwork


# The discharge-hole family: a hole in a wall, no pipe, no beta. 'fixed' and
# the ISO 5167 metering correlations are deliberately not in this map.
# The normed metering family (ISO 5167): Cd and dCd/dRe_D from C++'s
# metering_cd_and_dRe, so no metering formula lives here. Miller was removed:
# no published source matched the formula it carried.
_METERING_SELECTORS = {
    "ReaderHarrisGallagher": cb.MeteringCdCorrelation.ReaderHarrisGallagher,
}

_DISCHARGE_SELECTORS = {
    "IdelchikThick": cb.DischargeCdCorrelation.Idelchik1966Thick,
    "IdelchikBeveled": cb.DischargeCdCorrelation.Idelchik1966Beveled,
    "IdelchikRounded": cb.DischargeCdCorrelation.Idelchik1966Rounded,
    "McGreehanSchotsch": cb.DischargeCdCorrelation.McGreehanSchotsch1988,
    "Lichtarowicz": cb.DischargeCdCorrelation.Lichtarowicz1965,
}

# The inclined-hole inlet loss of a cooling hole's Cd.
# "idelchik_3_2": Idelchik diagram 3-2's inclined entrance put into Eq.
# (4-20)'s inlet coefficient of the thick-walled hole, see
# orifice::idelchik::inclination_* in include/orifice.h. It composes with
# 'IdelchikThick' only; any other correlation refuses it.
# "none": the angle enters the Cd only through the drilled length t/sin.
# "auto" (the default since the 2026-11-10 hole-Cd audit): "idelchik_3_2"
# for an INCLINED (< 90 deg) CYLINDRICAL hole on 'IdelchikThick', "none"
# otherwise -- a normal hole has no inclined inlet (the term is exactly zero
# there), so its element reads "none". ACCURACY, CROSS-SOURCE (Gritsch 1998,
# Hay 1983; neither is Idelchik's): the term moves every 30 deg cylindrical
# cell toward the data, and it worsens shaped holes (Gritsch fan -2 -> -11 %,
# laidback 0 -> -10 % bias without crossflow) because no shape term offsets
# it, so a shaped hole stays off until one exists.
CD_INCLINATION_MODELS = ("auto", "none", "idelchik_3_2")

# The EXIT (external) crossflow loss of a cooling hole's Cd, OPT-IN (#1074;
# default "none": the discharge node's gas velocity never reaches the Cd).
# "rowbury_2001": Rowbury, Oldfield & Lock (2001) Eqs. 7-8, an additive
# outlet loss on top of the hole's no-crossflow Cd (the correlation's, with
# any inclination term, or a 'fixed' one), 1/Cd^2 = 1/Cd0^2 + A k^-B - C,
# k the perpendicular outlet momentum flux ratio on the hole's actual flow;
# see rowbury_2001_* in include/crossflow_cd.h. The mainstream is the
# discharge MomentumChamberNode's gas, the same stream the plate's blowing
# ratio reads.
CD_EXIT_CROSSFLOW_MODELS = ("none", "rowbury_2001")
#: The Table 2 region; a flat wall takes the suction surface's, as the
#: authors do for Gritsch et al. (1998) in their Fig. 9.
CD_EXIT_CROSSFLOW_REGIONS = ("suction_surface", "pressure_surface", "leading_edge")

#: ImpingementPlateElement's crossflow Cd laws (#800).
_CD_CROSSFLOW_LAWS = ("rogers_hersh_1976", "florschuetz_metzger_su_1982", "none")

#: HoleCd.correlation: 'fixed', the normed metering family (ISO 5167) and the
#: discharge-hole family. Which of them an element takes is the element's
#: rule (the cooling plates refuse the metering family at validate).
HOLE_CD_CORRELATIONS = ("fixed", *_METERING_SELECTORS, *_DISCHARGE_SELECTORS)
#: HoleCd.inlet_crossflow: the supply-side crossflow U1/Vi McGreehan-Schotsch
#: reads. "auto": read from the backside duct of a duct-fed
#: EffusionPlateElement (its ``crossflow_segments``), 0 for every other hole.
HOLE_CD_INLET_CROSSFLOW_MODELS = ("auto", "none", "mcgreehan_schotsch")
#: HoleCd.exit_crossflow: the crossflow past the hole's exit. Rowbury et al.
#: (2001) on a film or effusion hole (opt-in, "auto" = "none"); the jet laws
#: of an impingement row ("auto" = "rogers_hersh_1976", #800).
HOLE_CD_EXIT_CROSSFLOW_MODELS = (
    "auto",
    "none",
    "rowbury_2001",
    "rogers_hersh_1976",
    "florschuetz_metzger_su_1982",
)
#: HoleCd.supply_drive: the pressure a hole's ideal flow (the Cd's reference)
#: is taken from (#1087, supply side; C++ hole_supply_drive). It differs from
#: the supply node's own Pt only on a DUCT-FED hole, whose supply node is a
#: duct station holding the duct's static pressure:
#:   "total"   the duct's stagnation pressure at the station, static plus the
#:             dynamic head of its mean flow: the reference of the plenum-fed
#:             Cd data (coolant total to exit static) and the hole-Cd audit's
#:             option (a) (2026-11-10). On any other hole, the supply node's
#:             own Pt, which is its total already.
#:   "static"  the station static: McGreehan & Schotsch's (1988) own
#:             convention, which their Eq. 17 inlet-crossflow Cd is defined
#:             on. Duct-fed holes only.
#:   "auto"    "static" on a duct-fed hole, "total" (its node's Pt) on every
#:             other. The duct-fed default is open, not settled (owner,
#:             2026-11-12): cross-source, the total drive helps co-aligned
#:             30 deg holes (Gritsch 1998, one lab: bias -33 % -> -8 %) and
#:             worsens normal ones (pooled RMS 17.5 % -> 29.6 %, including
#:             McGreehan & Schotsch's own Fig. 6, 14.5 % -> 48 %). Route to a
#:             decided default: Gritsch et al. 2001 (inclination x
#:             orientation x supply crossflow, wanted), and a check whether
#:             the published inlet-crossflow correlations convert to one
#:             common drive convention. A fraction of the head in between is
#:             not offered (no source gives it); Cd_supply_drive_head_fraction
#:             reports how much of the drive the choice decides.
HOLE_CD_SUPPLY_DRIVES = ("auto", "static", "total")
#: Deprecation of the per-element hole-Cd keywords (HoleCd replaces them).
HOLE_CD_KWARGS_DEPRECATED_IN = "0.12.0"
HOLE_CD_KWARGS_REMOVED_IN = "0.13.0"

#: The hole-Cd diagnostic group (_HoleCdChain.diagnostics): the same keys on
#: every hole element (OrificeElement and its subclasses), present whether or
#: not a term is on. The four ``*_model`` / ``Cd_correlation`` keys are text.
HOLE_CD_DIAGNOSTIC_KEYS = (
    "Cd",
    "Cd_correlation",
    "Cd_base",
    "Cd_in_range",
    "Cd_extrapolated_Re",
    "Cd_extrapolated_l_d",
    "Cd_extrapolated_edge",
    "Cd_inclination_model",
    "Cd_inclination_applied",
    "Cd_inclination_dzeta",
    "Cd_inlet_crossflow_model",
    "Cd_inlet_crossflow_applied",
    "Cd_inlet_crossflow_U1_over_Vi",
    "Cd_inlet_crossflow_beyond_8pct",
    "Cd_inlet_crossflow_degraded",
    "Cd_inlet_crossflow_ignored",
    "Cd_supply_drive_model",
    "Cd_supply_drive_applied",
    "Cd_supply_drive_head_fraction",
    "Cd_exit_crossflow_model",
    "Cd_exit_crossflow_applied",
    "Cd_exit_crossflow_ratio",
    "Cd_exit_crossflow_rho_u2",
    "Cd_exit_crossflow_k_perp",
    "Cd_exit_crossflow_k_used",
    "Cd_exit_crossflow_k_floor",
    "Cd_exit_crossflow_delta",
    "Cd_exit_crossflow_held",
    "Cd_exit_crossflow_faded",
    "Cd_exit_crossflow_floor_above_fit",
    "Cd_exit_crossflow_in_range",
    "Cd_exit_crossflow_extrapolated_velocity_ratio",
    "Cd_exit_crossflow_extrapolated_l_d",
    "Cd_exit_crossflow_extrapolated_grazing_velocity",
)
#: The group's text keys (the rest are numbers).
HOLE_CD_DIAGNOSTIC_TEXT_KEYS = (
    "Cd_correlation",
    "Cd_inclination_model",
    "Cd_inlet_crossflow_model",
    "Cd_supply_drive_model",
    "Cd_exit_crossflow_model",
)
#: The pre-0.12 names that shipped in 0.11.1, each its group key. Emitted
#: beside the group, where the element emitted them before, for one release.
HOLE_CD_DIAGNOSTIC_ALIASES = {
    "Cd_own": "Cd_base",
    "U1_over_Vi": "Cd_inlet_crossflow_U1_over_Vi",
    "crossflow_cd_beyond_8pct": "Cd_inlet_crossflow_beyond_8pct",
    "crossflow_cd_degraded": "Cd_inlet_crossflow_degraded",
    "coolant_crossflow_ignored": "Cd_inlet_crossflow_ignored",
}
#: The EffusionLiner summary's pre-0.12 names, each its new name.
LINER_CD_DIAGNOSTIC_ALIASES = {
    "U1_over_Vi_max": "Cd_inlet_crossflow_U1_over_Vi_max",
    "rows_U1_over_Vi": "rows_Cd_inlet_crossflow_U1_over_Vi",
    "crossflow_cd_beyond_8pct": "Cd_inlet_crossflow_beyond_8pct",
    "crossflow_cd_degraded": "Cd_inlet_crossflow_degraded",
}
HOLE_CD_DIAGNOSTIC_ALIASES_DEPRECATED_IN = "0.12.0"
HOLE_CD_DIAGNOSTIC_ALIASES_REMOVED_IN = "0.13.0"


@dataclass(frozen=True)
class HoleCd:
    """The discharge coefficient of a hole element, as one spec: the model
    choices only. The geometry (bore, hole count, angle, hole type, the exit
    term's region and orientation) stays on the element.

    Taken as ``hole_cd=`` by OrificeElement, EffusionPlateElement,
    FilmCoolingElement, SlotFilmElement, ImpingementPlateElement and the
    EffusionLiner and ImpingementArray builders. A field left at None or
    "auto" is THAT ELEMENT's default, which differs per element (see each
    field); every default is the element's behaviour before the spec.

    Attributes
    ----------
    correlation : str | None
        One of HOLE_CD_CORRELATIONS. None: the element's default --
        'ReaderHarrisGallagher' on OrificeElement (or its own ``correlation``
        keyword), 'IdelchikThick' on the effusion, film and slot elements,
        'McGreehanSchotsch' on an EffusionLiner's duct-fed panels, 'fixed'
        on an impingement row.
    Cd : float | None
        The 'fixed' correlation's value and nothing else (#872): given with
        any other correlation it raises. None with 'fixed': the element's
        documented default (Florschuetz's 0.79 on an impingement row) or a
        refusal where there is none.
    inclination : str
        The inclined-hole inlet term, one of CD_INCLINATION_MODELS. "auto":
        Idelchik diagram 3-2 on an inclined cylindrical IdelchikThick film or
        effusion hole, off otherwise (and on every element without an
        inclined hole). "idelchik_3_2" forces it (IdelchikThick, film or
        effusion only); "none" switches it off.
    inlet_crossflow : str
        The supply-side crossflow, one of HOLE_CD_INLET_CROSSFLOW_MODELS.
        "auto": McGreehan-Schotsch's U1/Vi from the backside duct of a
        duct-fed EffusionPlateElement, 0 (plenum-fed) everywhere else.
        "mcgreehan_schotsch" requires that duct and that correlation. "none"
        is no supply-crossflow term: U1/Vi is held at 0, so a duct-fed panel
        is taken as plenum-fed by the caller's choice, and the duct-fed
        panel's refusal of a correlation without the crossflow term (it
        would drop the duct silently) is lifted.
    exit_crossflow : str
        The crossflow past the hole exit, one of
        HOLE_CD_EXIT_CROSSFLOW_MODELS. "auto": 'rogers_hersh_1976' on an
        impingement row (#800), 'none' on every other hole. 'rowbury_2001'
        (film and effusion holes discharging into a flowing chamber, #1074);
        'rogers_hersh_1976' and 'florschuetz_metzger_su_1982' (impingement
        rows only); 'none' everywhere.
    supply_drive : str
        The pressure the hole's ideal flow is taken from, one of
        HOLE_CD_SUPPLY_DRIVES. It matters on a duct-fed EffusionPlateElement
        only, whose supply node is a duct station holding the duct's static:
        "total" drives from the duct's stagnation pressure there (static
        plus the dynamic head of its mean flow), "static" ("auto") from the
        station static (McGreehan & Schotsch's convention). Every other hole
        drives from its supply node's Pt, its total already ("auto" =
        "total"; "static" is refused there).
    """

    correlation: str | None = None
    Cd: float | None = None
    inclination: str = "auto"
    inlet_crossflow: str = "auto"
    exit_crossflow: str = "auto"
    supply_drive: str = "auto"

    def __post_init__(self) -> None:
        if self.correlation is not None and self.correlation not in HOLE_CD_CORRELATIONS:
            raise ValueError(
                f"HoleCd: correlation {self.correlation!r} is not one of "
                f"{list(HOLE_CD_CORRELATIONS)}."
            )
        if self.Cd is not None:
            if not math.isfinite(self.Cd):
                raise ValueError(f"HoleCd: Cd={self.Cd!r} must be a finite number.")
            if self.correlation is not None and self.correlation != "fixed":
                raise ValueError(
                    f"HoleCd: Cd={self.Cd!r} with correlation={self.correlation!r}. The "
                    "correlation computes its own Cd and a given one would be ignored. "
                    "For a fixed discharge coefficient pass correlation='fixed'; to use "
                    "the correlation, drop Cd."
                )
        for name, allowed in (
            ("inclination", CD_INCLINATION_MODELS),
            ("inlet_crossflow", HOLE_CD_INLET_CROSSFLOW_MODELS),
            ("exit_crossflow", HOLE_CD_EXIT_CROSSFLOW_MODELS),
            ("supply_drive", HOLE_CD_SUPPLY_DRIVES),
        ):
            value = getattr(self, name)
            if value not in allowed:
                raise ValueError(f"HoleCd: {name} {value!r} is not one of {list(allowed)}.")


#: Each deprecated per-element keyword's HoleCd field.
_LEGACY_HOLE_CD_FIELDS = {
    "correlation": "correlation",
    "Cd": "Cd",
    "cd_crossflow": "exit_crossflow",
}


def _hole_cd_spec(
    name: str, hole_cd: "HoleCd | None", legacy: dict[str, Any], stacklevel: int
) -> HoleCd:
    """The element's HoleCd: ``hole_cd``, or one built from the deprecated
    per-element keywords (``legacy``, None = not given) with one
    DeprecationWarning naming them. Both together raise."""
    if hole_cd is not None and not isinstance(hole_cd, HoleCd):
        raise ValueError(f"{name}: hole_cd is a HoleCd, not {type(hole_cd).__name__}")
    given = {k: v for k, v in legacy.items() if v is not None}
    if not given:
        return hole_cd if hole_cd is not None else HoleCd()
    kws = ", ".join(given)
    if hole_cd is not None:
        raise ValueError(
            f"{name}: hole_cd=HoleCd(...) and the deprecated keyword(s) {kws} together. "
            "Put every hole-Cd choice in the HoleCd."
        )
    fields = {_LEGACY_HOLE_CD_FIELDS[k]: v for k, v in given.items()}
    suggest = ", ".join(f"{f}={v!r}" for f, v in fields.items())
    warnings.warn(
        f"Hole-Cd keyword(s) {kws} of {name} are deprecated since "
        f"{HOLE_CD_KWARGS_DEPRECATED_IN} and will be removed in {HOLE_CD_KWARGS_REMOVED_IN}: "
        f"pass hole_cd=HoleCd({suggest}) instead.",
        DeprecationWarning,
        stacklevel=stacklevel,
    )
    try:
        return HoleCd(**fields)
    except ValueError as exc:
        raise ValueError(f"{name}: {exc}") from exc


def _with_correlation(spec: HoleCd, default: str, name: str) -> HoleCd:
    """``spec`` with its correlation resolved to the element's ``default``
    where it leaves it at None; a refusal that default causes (a Cd with a
    correlation other than 'fixed') names the element and says so."""
    if spec.correlation is not None:
        return spec
    try:
        return dataclasses.replace(spec, correlation=default)
    except ValueError as exc:
        raise ValueError(
            f"{name}: {exc} (correlation {default!r} is this element's default, taken "
            "because the HoleCd leaves correlation at None)"
        ) from exc


def _refuse_choice(name: str, field_name: str, value: str, allowed: tuple[str, ...]) -> None:
    """Refuse a HoleCd choice this element has no such term for."""
    if value not in allowed:
        raise ValueError(
            f"{name}: HoleCd({field_name}={value!r}) does not apply to this element; "
            f"it takes {list(allowed)}."
        )


def _resolve_inlet_crossflow(name: str, spec: HoleCd, duct_fed: bool) -> bool:
    """Whether the hole reads a supply-duct U1/Vi (_HoleCdChain.supply_duct)."""
    model = spec.inlet_crossflow
    if model == "auto":
        return duct_fed
    if model == "none":
        return False
    if not duct_fed:
        raise ValueError(
            f"{name}: HoleCd(inlet_crossflow='mcgreehan_schotsch') reads the supply duct's "
            "crossflow U1/Vi, and this hole is plenum-fed (U1/Vi = 0); a duct-fed "
            "EffusionPlateElement names its crossflow_segments."
        )
    if spec.correlation != "McGreehanSchotsch":
        raise ValueError(
            f"{name}: HoleCd(inlet_crossflow='mcgreehan_schotsch') is McGreehan-Schotsch's "
            f"own crossflow term; it does not compose with correlation={spec.correlation!r}."
        )
    return True


def _resolve_supply_drive(name: str, spec: HoleCd, duct_fed: bool) -> str:
    """The supply drive a hole is built with ("auto": "static" duct-fed, see
    HOLE_CD_SUPPLY_DRIVES; "total", its node's Pt, otherwise), "static"
    refused on a hole that is not duct-fed: its supply node carries its own
    total pressure, so there is no station static to drive from."""
    drive = spec.supply_drive
    if drive == "auto":
        drive = "static" if duct_fed else "total"
    if drive == "static" and not duct_fed:
        raise ValueError(
            f"{name}: HoleCd(supply_drive='static') is a duct-fed hole's drive from its "
            "duct station's static pressure; this hole is not duct-fed (its supply node "
            "carries its own total pressure). A duct-fed EffusionPlateElement names its "
            "crossflow_segments."
        )
    return drive


@dataclass
class _DrivePair:
    """The two pressures a hole's flow runs between (#1087), with their
    derivatives by unknown name.

    ``p_up`` is the supply side's: the drive the forward flow is taken from
    (HoleCd.supply_drive), and the pressure a reversed jet discharges
    against. ``p_dn`` is the discharge side's: the pressure the forward jet
    discharges against, and the one a reversed jet draws from. The flow's
    direction is the sign of ``p_up - p_dn``, so the forward and reverse
    forms vanish together and the flow passes zero continuously.

    ``supply`` is True where p_up is not the supply node's Pt (a duct-fed
    hole's reconstructed duct total, or a static drive from a station that
    carries its own Pt), ``exit_static`` where p_dn is the discharge node's
    STATIC pressure rather than its Pt (a cooling hole into a flowing
    chamber). ``q_supply`` is the duct's dynamic head at the station (0 when
    not duct-fed)."""

    p_up: float
    dp_up: dict[str, float]
    p_dn: float
    dp_dn: dict[str, float]
    supply: bool
    exit_static: bool
    q_supply: float = 0.0


# Unknown-name key of d(Cd)/d(own m_dot) THROUGH THE DISCHARGE CHAMBER'S
# through-flow (the hole's flow is one of the chamber's inflows). The
# self-consistent solves vary the hole flow at fixed chamber state, so this
# part is not the inner slope; the residuals fold it into the m_dot column
# as an ordinary dependency (_fold_chamber_own_flow).
_DCD_DM_VIA_CHAMBER = "__dCd_dm_via_chamber__"


def _own_flow_key(name: str, element_id: str) -> str:
    """Fold _DCD_DM_VIA_CHAMBER into the element's own m_dot column."""
    return f"{element_id}.m_dot" if name == _DCD_DM_VIA_CHAMBER else name


def _resolve_cd_inclination(
    name: str, cd_inclination: str, correlation: str, hole_type: str, angle_deg: float
) -> str:
    """The inclination model a cooling hole is built with: "auto" resolved
    from the correlation, the hole type and the angle, an explicit choice
    kept (and refused by _cd_inclination_deg where it cannot compose)."""
    if cd_inclination not in CD_INCLINATION_MODELS:
        raise ValueError(
            f"{name}: HoleCd inclination {cd_inclination!r} is not one of "
            f"{list(CD_INCLINATION_MODELS)}."
        )
    if cd_inclination != "auto":
        return cd_inclination
    if correlation == "IdelchikThick" and hole_type == "cylindrical" and angle_deg < 90.0:
        return "idelchik_3_2"
    return "none"


def _cd_inclination_deg(
    name: str, cd_inclination: str, correlation: str, angle_deg: float
) -> float:
    """The hole-to-wall angle the discharge correlation sees: 90 (no term)
    for "none", the hole angle for "idelchik_3_2", which needs
    'IdelchikThick'. ``cd_inclination`` is already resolved."""
    if cd_inclination == "none":
        return 90.0
    if correlation != "IdelchikThick":
        raise ValueError(
            f"{name}: HoleCd(inclination={cd_inclination!r}) is Idelchik diagram 3-2's "
            "inclined inlet, composed into the thick-walled hole's inlet coefficient "
            f"(Eq. 4-20); it does not compose with correlation={correlation!r}. Use "
            "correlation='IdelchikThick', or inclination='none' (or 'auto')."
        )
    return float(angle_deg)


@dataclass(frozen=True)
class _ExitCrossflowSpec:
    """The resolved opt-in: Table 2 set, angles [deg], and the exit/inlet
    area ratio the momentum flux is taken at (1, or a fan's)."""

    coefficient_set: Any
    theta_deg: float
    alpha_deg: float
    exit_area_ratio: float


def _exit_crossflow_spec(
    name: str,
    model: str,
    region: str,
    alpha_deg: float,
    theta_deg: float,
    fan_shaped: bool,
    exit_area_ratio: float,
) -> _ExitCrossflowSpec | None:
    """None (term off) unless ``model`` opts in; refuses a region Table 2
    has no set for, and angles where its A or B is not positive."""
    if model not in CD_EXIT_CROSSFLOW_MODELS:
        raise ValueError(
            f"{name}: HoleCd exit_crossflow {model!r} is not one of "
            f"{list(CD_EXIT_CROSSFLOW_MODELS)}."
        )
    if region not in CD_EXIT_CROSSFLOW_REGIONS:
        raise ValueError(
            f"{name}: exit_crossflow_region {region!r} is not one of "
            f"{list(CD_EXIT_CROSSFLOW_REGIONS)}."
        )
    if not (math.isfinite(alpha_deg) and 0.0 <= alpha_deg < 90.0):
        raise ValueError(
            f"{name}: exit_crossflow_orientation_deg = {alpha_deg!r} (alpha, the radial "
            "orientation from the surface normal, Rowbury et al. Fig. 3) must be in [0, 90)."
        )
    if model == "none":
        return None
    S = cb.Rowbury2001Set
    if fan_shaped:
        if region == "leading_edge":
            raise ValueError(
                f"{name}: Rowbury et al. (2001) Table 2 has no fan-shaped leading-edge set "
                "(N/A): use 'suction_surface' or 'pressure_surface'."
            )
        cset = (
            S.FanShapedSuctionSurface if region == "suction_surface" else S.FanShapedPressureSurface
        )
    else:
        cset = {
            "suction_surface": S.CylindricalSuctionSurface,
            "pressure_surface": S.CylindricalPressureSurface,
            "leading_edge": S.CylindricalLeadingEdge,
        }[region]
    try:
        cb.rowbury_2001_coefficients(cset, theta_deg, alpha_deg)
    except ValueError as exc:
        raise ValueError(f"{name}: HoleCd(exit_crossflow='rowbury_2001'): {exc}") from exc
    return _ExitCrossflowSpec(cset, float(theta_deg), float(alpha_deg), float(exit_area_ratio))


def _resolve_jet_row(name: str, spec: HoleCd, xn_d: float, yn_d: float, z_d: float) -> str:
    """An impingement row's jet crossflow law from its HoleCd ("auto" =
    'rogers_hersh_1976', #800), refused on a plate CR-3630 Table 4.1 does
    not cover. The holes are normal to the plate (no inclination term) and
    plenum-fed (no supply crossflow)."""
    _refuse_choice(name, "inclination", spec.inclination, ("auto", "none"))
    _refuse_choice(name, "exit_crossflow", spec.exit_crossflow, ("auto",) + _CD_CROSSFLOW_LAWS)
    _resolve_inlet_crossflow(name, spec, False)
    _resolve_supply_drive(name, spec, False)
    law = "rogers_hersh_1976" if spec.exit_crossflow == "auto" else spec.exit_crossflow
    _check_jet_crossflow_law(name, law, xn_d, yn_d, z_d)
    return law


def _check_jet_crossflow_law(name: str, law: str, xn_d: float, yn_d: float, z_d: float) -> None:
    """Refuse an unknown impingement crossflow law, and CR-3630 Table 4.1's
    on a plate the authors did not test."""
    if law not in _CD_CROSSFLOW_LAWS:
        raise ValueError(
            f"{name}: HoleCd(exit_crossflow={law!r}) is not one of {sorted(_CD_CROSSFLOW_LAWS)}"
        )
    if law == "florschuetz_metzger_su_1982" and not (
        cb.florschuetz_metzger_su_1982_has_plate(xn_d, yn_d, z_d)
    ):
        raise ValueError(
            f"{name}: HoleCd(exit_crossflow='florschuetz_metzger_su_1982') is "
            f"CR-3630 Table 4.1, which covers the plates the authors tested only: "
            f"B(5,4,1-3), B(5,8,1), B(10,4,1), B(10,8,1); not B({xn_d:g},{yn_d:g},"
            f"{z_d:g}). Use 'rogers_hersh_1976' (generic) or 'none'."
        )


# The geometry a base that does not read it is handed (C++ hole_cd_compose
# reads `hole` for the discharge family, `orifice` for the metering one).
_NO_HOLE = cb.DischargeHoleGeometry()
_NO_ORIFICE = cb.OrificeGeometry()


@dataclass
class _JetContext:
    """What the jet law's partials need to reach the unknowns."""

    up: Any
    up_id: str
    down_id: str
    dP: float
    rho: float
    drho_dT: float
    drho_dP: float
    G_i: float
    m_c: float


@dataclass
class _Context:
    """One evaluation's inputs' own derivatives by name (see _evaluate)."""

    base: str = "fixed"  # "fixed", "metering" or "discharge"
    Re: float = 0.0  # the Reynolds number the base read (hole, or pipe Re_D)
    re_d_moves: bool = False  # metering: Re_D follows the flow (pipe, |m| > 0)
    mu_fallback: bool = False  # metering: mu took the 1.8e-5 fallback
    du: dict[str, float] = field(default_factory=dict)  # d(U1/Vi)/d(name)
    q: float = 0.0  # mainstream momentum flux
    dq: dict[str, float] = field(default_factory=dict)
    exit_density: Any = None
    # The supply total the exit density expands from, by name: the supply
    # node's Pt, or a duct-fed hole's drive (#1087).
    dpt_up: dict[str, float] | None = None
    a_exit: float = 0.0
    jet: _JetContext | None = None


class _HoleCdChain:
    """The Cd of one hole element: its resolved options and every term.

    Options (resolved at the element's construction, see the module
    docstring for the order the terms compose in):
      inclination_deg   the hole-to-wall angle the discharge correlation sees
                        (90: no inclination term).
      inclination_model the inclination model configured, "auto" resolved
                        (an explicit "idelchik_3_2" at 90 deg stays itself).
      exit_xf           the exit-crossflow term's resolved spec, None: off.
      supply_duct       True on a hole whose supply crossflow U1/Vi is read
                        from a backside duct (a duct-fed EffusionPlateElement,
                        its ``crossflow_segments`` set) with the term on.
      duct_fed          True on a duct-fed hole whatever the term
                        (``inlet_crossflow="none"`` reads no U1/Vi, but the
                        diagnostics still report the duct's).
      jet_crossflow     the impingement crossflow law, "none" for every other
                        hole. ``jet_crossflow_suspended`` switches it off for
                        the solver's constant-Cd stage.
      inlet_crossflow_ignored  set by the element's topology: a supply
                        crossflow the Cd does not read (Cd_inlet_crossflow_ignored).

    Topology (set by the element's resolve_topology): the exit term's
    discharge node and approaching mainstream.
    """

    def __init__(self, owner: "OrificeElement") -> None:
        self._el = owner
        self.inclination_deg = 90.0
        self.inclination_model = "none"
        self.exit_xf: _ExitCrossflowSpec | None = None
        self.supply_duct = False
        self.duct_fed = False
        self.jet_crossflow = "none"
        self.jet_crossflow_suspended = False
        self.inlet_crossflow_ignored = False
        self.exit_node: Any = None
        self.exit_mainstream: tuple[str, ...] = ()
        self.exit_main_src: str | None = None
        # The drive pair (#1087, resolve_drive): the supply drive ("auto"
        # resolved), whether the supply station holds the duct's static (a
        # PlenumNode station, Pt = P) and whether the hole discharges against
        # its chamber's static pressure.
        self.supply_drive = "total"
        self.station_holds_static = False
        self.exit_static = False

    # --- options ------------------------------------------------------------

    def configure_cooling_hole(
        self,
        name: str,
        spec: HoleCd,
        *,
        hole_type: str,
        angle_deg: float,
        exit_crossflow_region: str,
        exit_crossflow_orientation_deg: float,
        exit_area_ratio: float,
        duct_fed: bool = False,
    ) -> tuple[str, str]:
        """Resolve an inclined cooling hole's options from its HoleCd (the
        film and effusion elements, ``spec.correlation`` already resolved):
        the inclination term ("auto" from the correlation, hole type and
        angle), the exit-crossflow term ("auto" = "none"; Rowbury et al.
        2001 opt-in) and the supply crossflow ("auto" = ``duct_fed``).
        Returns the inclination and exit-crossflow models in use."""
        correlation = str(spec.correlation)
        _refuse_choice(
            name, "exit_crossflow", spec.exit_crossflow, ("auto",) + CD_EXIT_CROSSFLOW_MODELS
        )
        model = _resolve_cd_inclination(name, spec.inclination, correlation, hole_type, angle_deg)
        self.inclination_deg = _cd_inclination_deg(name, model, correlation, angle_deg)
        self.inclination_model = model
        exit_model = "none" if spec.exit_crossflow == "auto" else spec.exit_crossflow
        self.exit_xf = _exit_crossflow_spec(
            name,
            exit_model,
            exit_crossflow_region,
            float(exit_crossflow_orientation_deg),
            angle_deg,
            hole_type == "shaped",
            exit_area_ratio,
        )
        self.duct_fed = bool(duct_fed)
        self.supply_duct = _resolve_inlet_crossflow(name, spec, duct_fed)
        self.supply_drive = _resolve_supply_drive(name, spec, duct_fed)
        return model, exit_model

    def configure_plain(self, name: str, spec: HoleCd) -> None:
        """A hole with no inclination, exit or supply crossflow term (a
        plain orifice, a film slot): the spec's terms must be "auto" or
        "none"."""
        _refuse_choice(name, "inclination", spec.inclination, ("auto", "none"))
        _refuse_choice(name, "exit_crossflow", spec.exit_crossflow, ("auto", "none"))
        self.supply_duct = _resolve_inlet_crossflow(name, spec, False)
        self.supply_drive = _resolve_supply_drive(name, spec, False)

    # --- the drive pair (#1087) ---------------------------------------------

    def resolve_drive(self, graph: "FlowNetwork", exit_static: bool) -> None:
        """Which pressures the hole's flow runs between, from the topology.

        Supply side: a duct-fed hole's station is a PlenumNode holding the
        duct's static (CrossflowSegmentElement's rule), so its total is
        rebuilt from the station's static state and mean flow; a station
        that is a MomentumChamberNode carries its own Pt. Exit side
        (``exit_static``, the film, effusion and slot elements): a hole
        discharging into a MomentumChamberNode takes its direction from the
        chamber's STATIC pressure, the pressure its jet discharges against,
        and a reversed jet draws from that static (#1087: none of the
        chamber's dynamic head is taken as recovered, as no source gives the
        fraction). A plenum's P is its Pt, so nothing changes there."""
        from .components import MomentumChamberNode, PlenumNode

        el = self._el
        src = graph.nodes.get(el.from_node)
        self.station_holds_static = self.duct_fed and isinstance(src, PlenumNode)
        self.exit_static = bool(exit_static) and isinstance(
            graph.nodes.get(el.to_node), MomentumChamberNode
        )

    def _station_drive(
        self,
        state_in: "NetworkMixtureState",
        flows: dict[str, float] | None,
        drive: "cb.HoleSupplyDrive",
    ) -> tuple[Any, dict[str, float]]:
        """C++ hole_supply_drive at a duct-fed hole's station (its static P,
        T0 and composition, the mean of the arriving and leaving segment
        flows, as U1/Vi reads them) and P_drive's derivatives by name."""
        el = self._el
        flows = flows or {}
        prev_id, next_id = el.crossflow_segments
        m_a = flows.get(prev_id, 0.0) if prev_id else 0.0
        m_b = flows.get(next_id, 0.0) if next_id else 0.0
        r = cb.hole_supply_drive(
            drive,
            float(state_in.P),
            float(state_in.Tt),
            list(state_in.Y),
            0.5 * (m_a + m_b),
            float(el.crossflow_area),
            el.regime == "incompressible",
        )
        d: dict[str, float] = {}
        coeff_m = 0.5 * float(r.dP_dm_dot)
        if coeff_m != 0.0:
            for elem, _eid in el._xf_elems:
                names = elem.unknowns()
                for local, c in elem.flow_jac_at_node(
                    el.from_node, list(range(len(names)))
                ).items():
                    d[names[local]] = d.get(names[local], 0.0) + coeff_m * c
        for key, v in (
            (f"{el.from_node}.P", float(r.dP_dP)),
            (f"{el.from_node}.Tt", float(r.dP_dTt)),
        ):
            if v != 0.0:
                d[key] = d.get(key, 0.0) + v
        for k, v in enumerate(r.dP_dY):
            if v != 0.0:
                key = f"{el.from_node}.Y[{k}]"
                d[key] = d.get(key, 0.0) + float(v)
        return r, d

    def drive_pair(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None,
    ) -> "_DrivePair | None":
        """The pressures the hole's flow runs between (_DrivePair), or None
        where they are the supply node's Pt and the discharge node's Pt --
        every hole that is not duct-fed with a moved drive and does not
        discharge into a flowing chamber, which then evaluates exactly as
        before."""
        el = self._el
        supply = False
        p_up, dp_up = float(state_in.Pt), {f"{el.from_node}.Pt": 1.0}
        q = 0.0
        if self.duct_fed and self.station_holds_static:
            if self.supply_drive == "total":
                r, dp_up = self._station_drive(state_in, flows, cb.HoleSupplyDrive.Total)
                p_up, q, supply = float(r.P_drive), float(r.q_dyn), True
            else:
                # The static drive reads the station as the hole always did
                # (its Pt, = P at a PlenumNode root); q for the diagnostics.
                q = self.station_head(state_in, flows)
        elif self.duct_fed:
            # A station carrying its own Pt (a MomentumChamberNode).
            q = float(state_in.Pt) - float(state_in.P)
            if self.supply_drive == "static":
                p_up, dp_up, supply = float(state_in.P), {f"{el.from_node}.P": 1.0}, True
        if not (supply or self.exit_static):
            return None
        if self.exit_static:
            p_dn, dp_dn = float(state_out.P), {f"{el.to_node}.P": 1.0}
        else:
            p_dn, dp_dn = float(state_out.Pt), {f"{el.to_node}.Pt": 1.0}
        return _DrivePair(p_up, dp_up, p_dn, dp_dn, supply, self.exit_static, q)

    def station_head(
        self, state_in: "NetworkMixtureState", flows: dict[str, float] | None
    ) -> float:
        """The duct's dynamic head at a duct-fed hole's station [Pa] (0 when
        not duct-fed): the momentum-chamber closure at its mean flow, or the
        station's own Pt - P where it carries one."""
        if not self.duct_fed:
            return 0.0
        if not self.station_holds_static:
            return float(state_in.Pt) - float(state_in.P)
        r, _ = self._station_drive(state_in, flows, cb.HoleSupplyDrive.Static)
        return float(r.q_dyn)

    def drive_pressures(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None,
    ) -> tuple[float, float]:
        """(p_up, p_dn) of drive_pair as numbers: the supply node's and the
        discharge node's Pt where it is None."""
        pair = self.drive_pair(state_in, state_out, flows)
        if pair is None:
            return float(state_in.Pt), float(state_out.Pt)
        return pair.p_up, pair.p_dn

    # --- the composed Cd ----------------------------------------------------

    def cd(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None = None,
    ) -> float:
        """The element's Cd: the hole's (hole_cd), under the impingement
        crossflow law where one acts."""
        return float(self._evaluate(state_in, state_out, flows, jet=True)[0].Cd)

    def dcd(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None = None,
        continue_re_floor: bool = False,
    ) -> dict[str, float]:
        """d(cd)/d(unknown) by name: the crossflow law's own partials where
        it acts (its Cd replaces the hole's), the hole's otherwise. The
        own-flow column through the discharge chamber rides under
        _DCD_DM_VIA_CHAMBER (see there)."""
        r, ctx = self._evaluate(
            state_in, state_out, flows, jet=True, continue_re_floor=continue_re_floor
        )
        return self._dcd_by_name(r, ctx, state_in)

    # --- one evaluation of the chain (C++ hole_cd_compose) --------------------

    def _evaluate(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None,
        *,
        jet: bool,
        continue_re_floor: bool = False,
    ) -> tuple[Any, "_Context"]:
        """The chain's inputs from the network state, composed by C++'s
        hole_cd_compose (include/hole_cd.h). Returns its HoleCdResult and the
        context the partials' mapping onto unknown names needs (each input's
        own derivatives by name). ``jet=False``: the hole's own Cd, no
        jet-crossflow law."""
        el = self._el
        geom = el._orifice_geom
        model = cb.HoleCdModel()
        inp = cb.HoleCdInputs()
        ctx = _Context()
        hole = _NO_HOLE
        orifice = _NO_ORIFICE
        if not el.use_correlation or geom is None:
            # A correlation with no geometry (before resolve_topology, or
            # bore == pipe, no orifice) has no Cd: the seed's stands in (#872),
            # clamped <= 1 (A_eff <= A_geom) in C++.
            model.base = cb.HoleCdBase.Fixed
            model.fixed_cd = el._seed_Cd
        else:
            metering = _METERING_SELECTORS.get(el.correlation)
            if metering is not None:
                # The ISO 5167 family reads the PIPE Reynolds number Re_D.
                model.base = cb.HoleCdBase.Metering
                model.metering = metering
                ctx.base = "metering"
                orifice = geom
                D_up = geom.D
                if D_up > 0.0:
                    mu = cb.transport_state(state_in.Tt, state_in.P, state_in.X).mu
                    ctx.mu_fallback = not mu > 0.0
                    mu = 1.8e-5 if ctx.mu_fallback else mu
                    m_dot_abs = abs(state_in.m_dot)
                    if m_dot_abs > 1e-12:
                        inp.Re_D = (4.0 * m_dot_abs) / (math.pi * D_up * mu)
                        ctx.re_d_moves = True
                        ctx.Re = inp.Re_D
                    else:
                        inp.Re_D = 1e4
                else:
                    inp.Re_D = 1e5  # typical reference when no upstream channel
            elif el.correlation in _DISCHARGE_SELECTORS:
                # The discharge-hole family reads the HOLE Reynolds number
                # (Re_D was wrong by 1-8% in a pipe, and froze a plenum-fed
                # hole at the 1e5 fallback), and McGreehan-Schotsch the
                # supply-side crossflow U1/Vi (0 for a plenum-fed hole).
                model.base = cb.HoleCdBase.Discharge
                model.discharge = _DISCHARGE_SELECTORS[el.correlation]
                model.continue_re_floor = continue_re_floor
                hole = self.discharge_hole()
                inp.U1_over_Vi, ctx.du = self.supply_crossflow(state_in, state_out, flows)
                inp.Re = ctx.Re = el._hole_reynolds(state_in)
                ctx.base = "discharge"
            else:
                raise ValueError(
                    f"OrificeElement: unknown correlation {el.correlation!r}. "
                    "The 'Auto' arm was removed: it picked a correlation from the "
                    "geometry behind the caller's back, which is how a "
                    "rounded-entry request came back as Stolz. Name one of "
                    "'fixed', 'ReaderHarrisGallagher', "
                    "'IdelchikThick', 'IdelchikBeveled', 'IdelchikRounded', "
                    "'McGreehanSchotsch', 'Lichtarowicz'."
                )
        spec = self.exit_xf
        if spec is not None and geom is not None:
            model.exit_crossflow = True
            model.exit_set = spec.coefficient_set
            model.exit_theta_deg = spec.theta_deg
            model.exit_alpha_deg = spec.alpha_deg
            self._exit_inputs(spec, state_in, state_out, inp, ctx, flows)
        if jet and self.jet_crossflow != "none" and not self.jet_crossflow_suspended:
            self._jet_inputs(state_in, state_out, flows, model, inp, ctx)
        return cb.hole_cd_compose(model, hole, orifice, inp), ctx

    def _exit_inputs(
        self,
        spec: _ExitCrossflowSpec,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        inp: Any,
        ctx: "_Context",
        flows: dict[str, float] | None = None,
    ) -> None:
        """The exit-crossflow term's inputs (#1074, rowbury_2001_* in
        crossflow_cd.h): k_perp = G_h^2 f / (rho_h rho_m u_m^2), f = u_perp /
        u_h (Fig. 3):
          G_h      the hole's ACTUAL mass flux, |m| / (A_holes AR), AR the
                   exit/inlet area ratio of a fan-shaped hole (the authors
                   take the exit-plane flux, p. 262), 1 otherwise;
          rho_h    the coolant's exit-plane static density, the isentropic
                   expansion from the supply (Pt, Tt) to the discharge
                   static P (rowbury_2001_exit_density; the supply total is
                   a duct-fed hole's drive, #1087): with it the
                   authors' own Fig. 9 predictions are reproduced (the
                   supply-Tt density left Ma_m 1.2 0.008 high);
          rho_m u_m^2  the APPROACHING mainstream's momentum flux, G_m^2 /
                   rho_m: G_m the chamber's main-inlet flow over its area,
                   rho_m the main stream's density at the chamber's main
                   face (mainstream_momentum_flux).
        No mainstream momentum: the term is off, Cd = Cd0. Eq. 8's own limit
        there is -C; the offset fades beyond the fitted k range, so Cd -> Cd0
        continuously as u_m -> 0 (include/crossflow_cd.h). Reverse
        (ingesting) flow takes the same Cd at |m|: the source covers
        injection only."""
        q, ctx.dq = self.mainstream_momentum_flux(state_out)
        ctx.q = q
        if not q > 0.0:
            return
        pt_up = state_in.Pt
        pair = self.drive_pair(state_in, state_out, flows) if self.duct_fed else None
        if pair is not None and pair.supply:
            pt_up, ctx.dpt_up = pair.p_up, pair.dp_up
        ex = cb.rowbury_2001_exit_density(pt_up, state_in.Tt, state_out.P, state_in.X)
        ctx.exit_density = ex
        ctx.a_exit = self._el.area * spec.exit_area_ratio
        inp.G_h = abs(float(state_in.m_dot)) / ctx.a_exit
        inp.rho_h = float(ex.rho)
        inp.rho_m_u_m2 = q

    def _jet_inputs(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None,
        model: Any,
        inp: Any,
        ctx: "_Context",
    ) -> None:
        """The impingement crossflow law's inputs (#800). G_ideal = sqrt(2 rho
        dP) along the flow: supply -> channel injecting, channel -> supply
        for a reversed jet (CR-2681's suction branch), rho the upstream
        static density. G_c = |m_c| / channel_area, m_c every crossflow
        source's flow into ``to_node``."""
        el = self._el
        forward = state_in.Pt >= state_out.P
        up, down = (state_in, state_out) if forward else (state_out, state_in)
        dP = float(up.Pt) - float(down.P)
        rho, drho_dT, drho_dP = _solver_tools.density_and_jacobians(
            float(up.Tt), float(up.P), list(up.X)
        )
        G_i = math.sqrt(2.0 * rho * dP) if dP > 0.0 else 0.0
        m_c = float(sum((flows or {}).get(eid, 0.0) for eid in el._crossflow_sources))
        if self.jet_crossflow == "rogers_hersh_1976":
            model.jet = cb.HoleCdJetCrossflow.RogersHersh1976
        else:
            model.jet = cb.HoleCdJetCrossflow.FlorschuetzMetzgerSu1982
            model.xn_d, model.yn_d, model.z_d = el.xn_d, el.yn_d, el.z_d
        inp.G_ideal = G_i
        inp.G_c = abs(m_c) / el.channel_area
        inp.suction = not forward
        ctx.jet = _JetContext(
            up=up,
            up_id=el.from_node if forward else el.to_node,
            down_id=el.to_node if forward else el.from_node,
            dP=dP,
            rho=rho,
            drho_dT=float(drho_dT),
            drho_dP=float(drho_dP),
            G_i=G_i,
            m_c=m_c,
        )

    # --- the partials by unknown name -----------------------------------------

    def _dcd_by_name(
        self, r: Any, ctx: "_Context", state_in: "NetworkMixtureState"
    ) -> dict[str, float]:
        """C++'s dCd/d(input) mapped onto the network's unknowns: each input's
        derivatives by name, chained. The jet law's own partials where its Cd
        replaces the hole's (none at all where it is capped or has no drive),
        the hole's otherwise: the discharge-hole family through the hole
        Reynolds number and U1/Vi, the ISO metering family through Re_D,
        'fixed' none; the exit-crossflow term on top.

        BELOW THE CORRELATION'S Re FLOOR the Cd value is held, so its true
        dCd/dRe is 0 -- in the m_dot column and in the viscosity columns
        alike. C++ continues the slope there as a Newton aid; that is only
        wanted by a search ON Cd (``continue_re_floor=True``, the
        self-consistent inner solve). Every sensitivity of a flow built on
        the held value takes the zero: with the continued slope, the
        implicit-function denominator 1 - K dCd/dm of a McGreehan-Schotsch
        panel below Re = 1e4 scaled its whole row by ~0.987 (#760, #801).
        hole_cd_compose applies it."""
        if r.jet_acts:
            return self._jet_dcd_by_name(r, ctx.jet) if r.jet_partials else {}
        el = self._el
        out: dict[str, float] = {}
        m = float(state_in.m_dot)
        if ctx.base == "discharge":
            # Through the hole Reynolds number (linear in m_dot, so dRe/dm =
            # Re/m; C++ zeroes dCd/dRe below the correlation's floor unless
            # continue_re_floor, see discharge_cd_re_floor) and, for a
            # duct-fed hole, through U1/Vi.
            Re = ctx.Re
            dCd_dRe = float(r.dCd_dRe)
            if abs(m) > 1e-12:
                out[f"{el.id}.m_dot"] = dCd_dRe * Re / abs(m) * (1.0 if m > 0 else -1.0)
                if el._orifice_geom.d > 0.0:
                    for name, d in self._re_viscosity_dnames(state_in, Re).items():
                        out[name] = out.get(name, 0.0) + dCd_dRe * d
            dCd_dU = float(r.dCd_dU1_over_Vi)
            for name, d in ctx.du.items():
                out[name] = out.get(name, 0.0) + dCd_dU * d
        elif ctx.base == "metering" and ctx.re_d_moves:
            # Through the pipe Reynolds number Re_D = 4 |m| / (pi D mu),
            # linear in the flow. The row's flow column missed it: 1% off on
            # a 2 cm orifice at 1600 K. Held fixed where Re_D is (no upstream
            # pipe, |m| <= 1e-12).
            Re_D = ctx.Re
            dCd_dRe = float(r.dCd_dRe_D)
            out[f"{el.id}.m_dot"] = dCd_dRe * Re_D / abs(m) * (1.0 if m > 0 else -1.0)
            if not ctx.mu_fallback:
                for name, d in self._re_viscosity_dnames(state_in, Re_D).items():
                    out[name] = out.get(name, 0.0) + dCd_dRe * d
        if r.exit_active and r.exit.dCd_dk != 0.0:
            self._exit_dcd_by_name(r, ctx, m, out)
        return out

    def follows_own_flow(self) -> bool:
        """True when Cd moves with this element's own flow through a hole
        Reynolds number with an analytic derivative: the discharge-hole
        family on a correlation. The ISO metering family also reads Re_D but
        carries no dCd/dRe (its documented gap), so it keeps the iterate form.
        """
        el = self._el
        if self.exit_xf is not None and el._orifice_geom is not None:
            # The exit-crossflow loss reads the hole's own flow through
            # k_perp, whatever the Cd0 (a 'fixed' one included).
            return True
        return (
            el.use_correlation
            and el._orifice_geom is not None
            and _DISCHARGE_SELECTORS.get(el.correlation) is not None
        )

    # --- 1. the base Cd -----------------------------------------------------

    def discharge_hole(self) -> "cb.DischargeHoleGeometry":
        """One real hole, as the discharge-hole correlations see it."""
        el = self._el
        hole = cb.DischargeHoleGeometry(
            d=el._orifice_geom.d,
            L=el.plate_thickness,
            r=el.edge_radius,
        )
        hole.bevel = el.bevel_depth
        if self.inclination_deg < 90.0:
            hole.inclination_deg = self.inclination_deg
        return hole

    def _re_viscosity_dnames(self, state_in: "NetworkMixtureState", Re: float) -> dict[str, float]:
        """d(Re)/d(upstream T, P) through mu: Re ~ |m| / mu, so
        dRe/dx = -Re dmu/dx / mu. The Cd's temperature column, missing until
        a heated plenum feeding a Re-dependent orifice showed its row 1e-3 off
        in the column that moves the plenum's T. And through mu's composition
        (Wilke's rule, analytic in C++): the mass-fraction columns, missing
        until #964 -- with a Reader-Harris/Gallagher Cd behind a combustor
        the row's fuel-flow column was 7e-4 off. Empty where Re used the mu
        fallback."""
        el = self._el
        X = list(state_in.X)
        mu, dmu_dT, dmu_dP = _solver_tools.viscosity_and_jacobians(state_in.Tt, state_in.P, X)
        if not mu > 0.0:
            return {}
        out = {f"{el.from_node}.Tt": -Re * float(dmu_dT) / mu}
        if dmu_dP != 0.0:
            out[f"{el.from_node}.P"] = -Re * float(dmu_dP) / mu
        dmu_dY = _solver_tools.viscosity_mass_gradient(state_in.Tt, state_in.P, X)
        for k, d in enumerate(dmu_dY):
            if d != 0.0:
                out[f"{el.from_node}.Y[{k}]"] = -Re * float(d) / mu
        return out

    def supply_crossflow(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None,
        measure: bool = False,
    ) -> tuple[float, dict[str, float]]:
        """McGreehan-Schotsch's supply-side U1/Vi and d(U1/Vi)/d(name).

        0 for a plenum-fed hole. A duct-fed one (``supply_duct``: an
        EffusionPlateElement with crossflow segments) reads the backside
        duct's mean velocity at its supply station over the static-referenced
        ideal jet velocity. ``measure=True`` computes the duct's ratio on a
        duct-fed hole whose term is off (``inlet_crossflow="none"``), for the
        diagnostics only.
        """
        if not (self.supply_duct or (measure and self.duct_fed)):
            return 0.0, {}
        from .components import _moving_static_T

        el = self._el
        flows = flows or {}
        prev_id, next_id = el.crossflow_segments
        m_a = flows.get(prev_id, 0.0) if prev_id else 0.0
        m_b = flows.get(next_id, 0.0) if next_id else 0.0
        # The duct's STATIC state at the station (#782): its node carries T0,
        # and the station's gas moves at the same mean flow U1 is built from.
        st = _moving_static_T(
            state_in,
            0.5 * (m_a + m_b),
            el.crossflow_area,
            el.regime == "incompressible",
        )
        r = cb.crossflow_velocity_ratio(
            m_a, m_b, state_in.P, st.T, state_in.X, state_out.P, el.crossflow_area
        )
        d_dm_T = 0.5 * r.d_dT * st.dm
        d: dict[str, float] = {}
        for elem, eid in el._xf_elems:
            coeff_m = (r.d_dm_a if eid == prev_id else r.d_dm_b) + d_dm_T
            names = elem.unknowns()
            for local, c in elem.flow_jac_at_node(el.from_node, list(range(len(names)))).items():
                d[names[local]] = d.get(names[local], 0.0) + coeff_m * c
        for key, v in (
            (f"{el.from_node}.P", r.d_dP + r.d_dT * st.dP),
            (f"{el.from_node}.Tt", r.d_dT * st.dT0),
            (f"{el.to_node}.P", r.d_dP_down),
        ):
            d[key] = d.get(key, 0.0) + v
        for k, dT_dY in st.dY.items():
            key = f"{el.from_node}.Y[{k}]"
            d[key] = d.get(key, 0.0) + r.d_dT * dT_dY
        return float(r.U1_over_Vi), d

    def pipe_cd_ratio(self, beta: float) -> float:
        """C_pipe / Cd_wall of a discharge-hole Cd carried into an approach
        pipe of this beta (#645); 1 for any other correlation."""
        el = self._el
        selector = _DISCHARGE_SELECTORS.get(el.correlation)
        if not (el.use_correlation and selector is not None):
            return 1.0
        return float(cb.discharge_pipe_cd_ratio(selector, self.discharge_hole(), beta))

    def validate(self) -> None:
        """Lichtarowicz refuses short holes (l/d < 1.5, the source's own
        hysteresis warning). Ask it once here, so the refusal names this
        element at set-up rather than surfacing mid-solve. The limit lives
        in C++ only."""
        el = self._el
        if el.correlation == "Lichtarowicz" and el._orifice_geom is not None:
            try:
                cb.discharge_cd(
                    _DISCHARGE_SELECTORS[el.correlation],
                    self.discharge_hole(),
                    cb.DischargeHoleState(Re=1.0e4),
                )
            except ValueError as exc:
                raise ValueError(
                    f"{type(el).__name__} {el.id!r}: plate_thickness/d = "
                    f"{el.plate_thickness / el._orifice_geom.d:.3g}. {exc}"
                ) from exc

    # --- 2. the exit-crossflow term (Rowbury et al. 2001) ---------------------

    def _exit_dcd_by_name(self, r: Any, ctx: "_Context", m: float, out: dict[str, float]) -> None:
        """Add the exit-crossflow term's dCd/dk dk/dx to ``out``: C++'s
        dCd/dG_h, dCd/drho_h and dCd/d(rho_m u_m^2) through each input's
        derivatives by name (G_h = |m| / a_exit, rho_h's isentropic exit
        state, the mainstream's momentum flux at the chamber's main face)."""
        el = self._el
        if m != 0.0:
            key = f"{el.id}.m_dot"
            out[key] = (
                out.get(key, 0.0) + float(r.dCd_dG_h) * (1.0 if m > 0.0 else -1.0) / ctx.a_exit
            )
        ex = ctx.exit_density
        c_rho = float(r.dCd_drho_h)
        if ctx.dpt_up is None:
            terms = [(f"{el.from_node}.Pt", float(ex.drho_dPt))]
        else:
            terms = [(k, float(ex.drho_dPt) * v) for k, v in ctx.dpt_up.items()]
        terms += [(f"{el.from_node}.Tt", ex.drho_dTt), (f"{el.to_node}.P", ex.drho_dP)]
        for key, d in terms:
            out[key] = out.get(key, 0.0) + c_rho * float(d)
        for k, d in enumerate(ex.drho_dY):
            if d != 0.0:
                key = f"{el.from_node}.Y[{k}]"
                out[key] = out.get(key, 0.0) + c_rho * float(d)
        c_q = float(r.dCd_drho_m_u_m2)
        for name, d in ctx.dq.items():
            out[name] = out.get(name, 0.0) + c_q * d

    def resolve_exit_mainstream(self, graph: "FlowNetwork") -> None:
        """The discharge node and the elements whose flow into it is the
        APPROACHING mainstream: its main inlet, or without one every inflow
        that is not a cooling hole (the holes' own discharge is not the
        crossflow they blow into)."""
        node = graph.nodes.get(self._el.to_node)
        self.exit_node = node
        self.exit_main_src = None
        main = getattr(node, "main_inlet", None)
        if main is not None:
            # A chamber with several inflows must declare its main inlet, so
            # for a hole discharging into a flowing chamber this is the case.
            self.exit_mainstream = (main,)
            elem = graph.elements.get(main)
            self.exit_main_src = getattr(elem, "from_node", None)
            return
        # Without one the chamber's only inflow is this hole: no mainstream
        # (refused at validate).
        self.exit_mainstream = ()

    def mainstream_momentum_flux(
        self, state_out: "NetworkMixtureState"
    ) -> tuple[float, dict[str, float]]:
        """rho_m u_m^2 of the approaching mainstream and its derivatives by
        name, or (0, {}) where the discharge node has none.

        G_m = (the main inlet's flow into the chamber) / A, rho_m the main
        stream's density at the chamber's MAIN FACE: the face's static state
        (MomentumChamberNode.main_face_state, C++ chamber_merge_face_state)
        on the main stream's isentrope from (Pt_face, Tt_main, Y_main) to
        P_face. Not the chamber's own static density, whose Tt and Y are the
        mixed outlet's, holes included. d/d(own m_dot) through the face
        (the outlet flow and the side streams' impulse include this hole's)
        is keyed _DCD_DM_VIA_CHAMBER (see there)."""
        from .components import MomentumChamberNode

        node = self.exit_node
        if not isinstance(node, MomentumChamberNode) or not node.area > 0.0:
            return 0.0, {}
        main = getattr(node, "main_inlet", None)
        main_st, dm = None, {}
        for eid, st, names in getattr(node, "_inflows", []):
            if eid == main:
                main_st, dm = st, names
        if main_st is None or float(main_st.m_dot) == 0.0:
            return 0.0, {}
        m_m = float(main_st.m_dot)
        P_f, Pt_f, dPf, dPtf = node.main_face_state(state_out)
        ex = cb.rowbury_2001_exit_density(Pt_f, float(main_st.Tt), P_f, main_st.X)
        rho = float(ex.rho)
        drho: dict[str, float] = {}
        for src, coef in ((dPf, float(ex.drho_dP)), (dPtf, float(ex.drho_dPt))):
            for var, d in src.items():
                drho[var] = drho.get(var, 0.0) + coef * d
        src_node = self.exit_main_src
        if src_node is not None:
            key = f"{src_node}.Tt"
            drho[key] = drho.get(key, 0.0) + float(ex.drho_dTt)
            for k, d in enumerate(ex.drho_dY):
                if d != 0.0:
                    key = f"{src_node}.Y[{k}]"
                    drho[key] = drho.get(key, 0.0) + float(d)
        A = node.area
        q = (m_m / A) ** 2 / rho
        own = f"{self._el.id}.m_dot"
        dq: dict[str, float] = {}
        dq_dm = 2.0 * m_m / (A * A * rho)
        for var, coeff in dm.items():
            if coeff != 0.0:
                key = _DCD_DM_VIA_CHAMBER if var == own else var
                dq[key] = dq.get(key, 0.0) + dq_dm * coeff
        for var, d in drho.items():
            key = _DCD_DM_VIA_CHAMBER if var == own else var
            dq[key] = dq.get(key, 0.0) - q / rho * d
        return q, dq

    def validate_exit_crossflow(self) -> None:
        """The exit-crossflow opt-in needs a mainstream: a flowing discharge
        chamber. Refused, not ignored, on anything else (a plenum has no gas
        velocity, and the term would silently do nothing)."""
        from .components import MomentumChamberNode

        if self.exit_xf is None:
            return
        el = self._el
        node = self.exit_node
        if not (isinstance(node, MomentumChamberNode) and node.area > 0.0):
            raise ValueError(
                f"{type(el).__name__} {el.id!r}: HoleCd(exit_crossflow="
                f"{getattr(el, 'cd_exit_crossflow', '')!r}) reads the mainstream's momentum "
                f"flux from the discharge node, which must be a MomentumChamberNode with an "
                f"area; {el.to_node!r} has no gas velocity."
            )
        if not self.exit_mainstream:
            raise ValueError(
                f"{type(el).__name__} {el.id!r}: HoleCd(exit_crossflow=...) needs an approaching "
                f"mainstream, and {el.to_node!r} is fed only by cooling holes (declare "
                "its main_inlet, or feed it a gas stream)."
            )

    # --- 3. the impingement jet-crossflow law (#800) --------------------------

    def _jet_dcd_by_name(self, r: Any, jc: "_JetContext") -> dict[str, float]:
        """The jet law's dCd/dG_ideal and dCd/dG_c by name: G_ideal through
        the drive dP and the upstream density (T, P and, #964, composition),
        G_c through every crossflow source's flow into ``to_node``."""
        el = self._el
        out: dict[str, float] = {}

        def add(name: str, v: float) -> None:
            if v != 0.0:
                out[name] = out.get(name, 0.0) + v

        rho, dP, G_i, up, up_id, down_id = jc.rho, jc.dP, jc.G_i, jc.up, jc.up_id, jc.down_id
        dGi_ddP = rho / G_i
        dGi_drho = dP / G_i
        dcd_dGi = float(r.dCd_dG_ideal)
        add(f"{up_id}.Pt", dcd_dGi * dGi_ddP)
        add(f"{down_id}.P", -dcd_dGi * dGi_ddP)
        add(f"{up_id}.Tt", dcd_dGi * dGi_drho * jc.drho_dT)
        add(f"{up_id}.P", dcd_dGi * dGi_drho * jc.drho_dP)
        # rho = P MW / (R_u T) through the composition (#964): at fixed (T, P)
        # d(rho)/dY_k = rho (1 - MW / W_k) on the normalised mass fractions
        # the kernels differentiate (the same along any sum-preserving change).
        W = float(cb.mwmix(up.X))
        c_rho = dcd_dGi * dGi_drho * rho
        for k in range(len(up.Y)):
            add(f"{up_id}.Y[{k}]", c_rho * (1.0 - W / float(cb.species_molar_mass(k))))
        m_c = jc.m_c
        if m_c != 0.0:
            dcd_dmc = float(r.dCd_dG_c) * (1.0 if m_c > 0.0 else -1.0) / el.channel_area
            for eid in el._crossflow_sources:
                src = el._source_elements[eid]
                names = src.unknowns()
                for local, coeff in src.flow_jac_at_node(
                    el.to_node, list(range(len(names)))
                ).items():
                    add(names[local], dcd_dmc * coeff)
        return out

    # --- range flags and diagnostics ------------------------------------------

    def in_range(self, state_in: "NetworkMixtureState") -> bool:
        """Is the Cd inside its source's stated range?

        A flag for diagnostics, never a selector. True for 'fixed' (the value
        is the caller's). The normed metering family asks C++'s
        metering_cd_in_range at the pipe Reynolds number (ISO 5167-2:2003
        5.3.1 for Reader-Harris/Gallagher, below whose lower Re_D bound the
        Cd is held (#873) and this reads False). A plenum-fed metering plate reads False: it has no pipe.
        """
        el = self._el
        metering = _METERING_SELECTORS.get(el.correlation)
        if metering is not None and el._orifice_geom is not None:
            D_up = el._orifice_geom.D
            mu = cb.transport_state(state_in.Tt, state_in.P, state_in.X).mu
            m = abs(float(state_in.m_dot))
            Re_D = 4.0 * m / (math.pi * D_up * mu) if (mu > 0.0 and D_up > 0.0) else 0.0
            in_pipe = bool(el.upstream_diameter)
            return in_pipe and bool(cb.metering_cd_in_range(metering, el._orifice_geom, Re_D))
        selector = _DISCHARGE_SELECTORS.get(el.correlation)
        if selector is None or el._orifice_geom is None:
            return True
        return bool(
            cb.discharge_cd_in_range(
                selector,
                self.discharge_hole(),
                cb.DischargeHoleState(Re=el._hole_reynolds(state_in)),
            )
        )

    def envelope(self, state_in: "NetworkMixtureState") -> dict[str, float]:
        """Which input of a discharge-hole Cd lies outside its source's range,
        by variable (C++ discharge_cd_envelope): 1.0 per variable that is.

        INFORMATION stated in the results, never a warning: outside, the Cd is
        held or continued as the correlation documents, and the user can do
        nothing better than know which variable it is. ``Cd_extrapolated_l_d``
        on a thick hole beyond l/d 4 is a gentle friction-only continuation;
        on McGreehan-Schotsch beyond l/d 10 its length function keeps falling
        steeply. 0.0 for a metering or 'fixed' Cd (``Cd_in_range`` covers the
        metering family's ISO range).
        """
        el = self._el
        out = {"Cd_extrapolated_Re": 0.0, "Cd_extrapolated_l_d": 0.0, "Cd_extrapolated_edge": 0.0}
        selector = _DISCHARGE_SELECTORS.get(el.correlation)
        if selector is None or el._orifice_geom is None:
            return out
        env = cb.discharge_cd_envelope(
            selector,
            self.discharge_hole(),
            cb.DischargeHoleState(Re=el._hole_reynolds(state_in)),
        )
        out["Cd_extrapolated_Re"] = float(env.re)
        out["Cd_extrapolated_l_d"] = float(env.l_over_d)
        out["Cd_extrapolated_edge"] = float(env.edge)
        return out

    def inclination_dzeta(self) -> float:
        """The inclined-inlet increment on zeta' (0.0 when off or at 90 deg)."""
        return float(cb.discharge_inclination_dzeta(self.inclination_deg))

    def diagnostics(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None,
    ) -> dict[str, float | str]:
        """The hole-Cd diagnostic group: every key in HOLE_CD_DIAGNOSTIC_KEYS,
        on every hole element, whether or not a term is on (a term that is
        off reads its model "none", applied 0 and neutral values). One
        evaluation of the chain; the group's ``Cd`` is the element's.

        The terms, in the order they compose (see the module docstring):
          base        Cd_correlation, Cd_base (term 1, before any exit or
                      jet crossflow term), the range flags;
          inclination Cd_inclination_model / _applied / _dzeta;
          inlet       Cd_inlet_crossflow_model / _applied / _U1_over_Vi (a
                      duct-fed hole's U1/Vi, reported where the term is off;
                      0 plenum-fed) / _beyond_8pct / _degraded (the
                      Rohde-scored bounds, 0 where the term is off) /
                      _ignored (a supply crossflow the Cd does not read);
          exit        Cd_exit_crossflow_model / _applied / _ratio (Cd /
                      Cd_base, 1 when off), the Rowbury et al. (2001) state
                      (NaN k values where it is not evaluated) and the Rogers
                      & Hersh envelope flags (0 where the law does not act).
        """
        el = self._el
        r, ctx = self._evaluate(state_in, state_out, flows, jet=True)
        nan = math.nan
        # 2. inlet (supply) crossflow, McGreehan-Schotsch.
        u1_vi, beyond, degraded = 0.0, 0.0, 0.0
        if self.duct_fed:
            u1_vi, _ = self.supply_crossflow(state_in, state_out, flows, measure=True)
            on = self.supply_duct
            beyond = float(on and u1_vi > 0.2)
            degraded = float(on and u1_vi > 0.35)
        # 3. exit crossflow: Rowbury et al. (2001), or an impingement row's
        # jet law (HoleCd.exit_crossflow names both).
        rowbury = self.exit_xf is not None and el._orifice_geom is not None
        exit_on = bool(r.exit_active)
        if self.exit_xf is not None:
            exit_model = "rowbury_2001"
            exit_applied = exit_on
        else:
            exit_model = self.jet_crossflow
            exit_applied = bool(r.jet_acts)
        if rowbury:
            ratio = float(r.exit_ratio)
        elif exit_model != "none":
            ratio = float(r.Cd) / float(r.cd_base) if r.cd_base > 0.0 else nan
        else:
            ratio = 1.0
        out: dict[str, float | str] = {
            "Cd": float(r.Cd),
            "Cd_correlation": str(el.correlation),
            "Cd_base": float(r.cd_base),
            "Cd_in_range": float(self.in_range(state_in)),
            **self.envelope(state_in),
            "Cd_inclination_model": self.inclination_model,
            "Cd_inclination_applied": float(self.inclination_deg < 90.0),
            "Cd_inclination_dzeta": self.inclination_dzeta(),
            "Cd_inlet_crossflow_model": "mcgreehan_schotsch" if self.supply_duct else "none",
            "Cd_inlet_crossflow_applied": float(self.supply_duct),
            "Cd_inlet_crossflow_U1_over_Vi": float(u1_vi),
            "Cd_inlet_crossflow_beyond_8pct": beyond,
            "Cd_inlet_crossflow_degraded": degraded,
            "Cd_inlet_crossflow_ignored": float(self.inlet_crossflow_ignored),
            **self.supply_drive_diagnostics(state_in, state_out, flows),
            "Cd_exit_crossflow_model": exit_model,
            "Cd_exit_crossflow_applied": float(exit_applied),
            "Cd_exit_crossflow_ratio": ratio,
            "Cd_exit_crossflow_rho_u2": float(ctx.q) if rowbury else 0.0,
            # The actual k_perp, and the one the term was evaluated at (the
            # floor where held).
            "Cd_exit_crossflow_k_perp": float(r.k_perp.k) if exit_on else nan,
            "Cd_exit_crossflow_k_used": float(r.exit.k_used) if exit_on else nan,
            "Cd_exit_crossflow_k_floor": float(r.exit_k_floor) if rowbury else nan,
            "Cd_exit_crossflow_delta": float(r.exit.delta) if exit_on else 0.0,
            "Cd_exit_crossflow_held": float(r.exit.held) if exit_on else 0.0,
            # Where k_used sits against the set's k tiers (info): 0 inside
            # what its Table 1 rows provably reach (tier 1, data), 0.5 beyond
            # that but inside the authors' own application of the set (tier
            # 2, Fig. 9 for the suction set; Eq. 8 unfaded), 1 beyond both:
            # the offset fades.
            "Cd_exit_crossflow_faded": float(r.exit_k_tier),
            # Off Table 1's angles the floor can lie above the fitted k
            # range: the whole fitted range is then held (info).
            "Cd_exit_crossflow_floor_above_fit": float(r.exit_floor_above_fit),
            "Cd_exit_crossflow_in_range": float(r.exit_in_range) if rowbury else 1.0,
            **self.jet_crossflow_envelope(state_in, state_out, flows, float(r.cd_hole)),
        }
        return out

    def supply_drive_diagnostics(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None,
    ) -> dict[str, float | str]:
        """The supply-drive keys of the group (#1087):
        Cd_supply_drive_model          the drive in use ("auto" resolved);
        Cd_supply_drive_applied        1 where it is not the supply node's
                                       own Pt (a duct-fed hole's rebuilt
                                       duct total, or a static drive from a
                                       station that carries its own Pt);
        Cd_supply_drive_head_fraction  INFORMATION: the share of the full
                                       drive (duct total less the discharge
                                       pressure) that is the duct's dynamic
                                       head, q / (Pt_duct - p_dn): how much
                                       of the drive the unsourced recovered
                                       fraction decides (0 "static" and 1
                                       "total" bound it). 0 when not
                                       duct-fed; NaN where the duct total
                                       does not exceed the discharge
                                       pressure."""
        pair = self.drive_pair(state_in, state_out, flows)
        fraction = 0.0
        if self.duct_fed:
            q = pair.q_supply if pair is not None else self.station_head(state_in, flows)
            p_dn = pair.p_dn if pair is not None else float(state_out.Pt)
            pt_duct = float(state_in.P) + q if self.station_holds_static else float(state_in.Pt)
            fraction = q / (pt_duct - p_dn) if pt_duct > p_dn else math.nan
        return {
            "Cd_supply_drive_model": self.supply_drive,
            "Cd_supply_drive_applied": float(pair is not None and pair.supply),
            "Cd_supply_drive_head_fraction": float(fraction),
        }

    def legacy_diagnostics(self, group: dict[str, float | str], keys: tuple[str, ...]) -> dict:
        """The deprecated pre-0.12 names among ``keys`` (each a key of
        HOLE_CD_DIAGNOSTIC_ALIASES), valued from ``group``: emitted beside
        the group for one release, where the element emitted them before."""
        return {k: group[HOLE_CD_DIAGNOSTIC_ALIASES[k]] for k in keys}

    def jet_crossflow_envelope(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None,
        base_cd: float,
    ) -> dict[str, float]:
        """Where Rogers & Hersh's law ACTS (not capped), which of its inputs
        lie outside CR-2681's envelope, by variable (C++
        rogers_hersh_1976_envelope): the velocity ratio Vi/V_inf = Cd q (the
        law's equal-density reading), the hole's L/D against the branch's own
        band (a reversed jet is the suction branch, whose L/D effect CR-2681
        calls considerable), and the grazing velocity G_c / rho of the
        channel. INFORMATION stated in the results, never a
        warning. 0.0 everywhere for another law, or where the hole's own Cd
        holds."""
        out = {
            "Cd_exit_crossflow_extrapolated_velocity_ratio": 0.0,
            "Cd_exit_crossflow_extrapolated_l_d": 0.0,
            "Cd_exit_crossflow_extrapolated_grazing_velocity": 0.0,
        }
        if self.jet_crossflow != "rogers_hersh_1976" or self.jet_crossflow_suspended:
            return out
        el = self._el
        forward = state_in.Pt >= state_out.P
        up, down = (state_in, state_out) if forward else (state_out, state_in)
        dP = float(up.Pt) - float(down.P)
        rho = float(cb.density(float(up.Tt), float(up.P), list(up.X)))
        G_i = math.sqrt(2.0 * rho * dP) if dP > 0.0 else 0.0
        m_c = float(sum((flows or {}).get(eid, 0.0) for eid in el._crossflow_sources))
        G_c = abs(m_c) / el.channel_area
        r = cb.rogers_hersh_1976_cd(base_cd, G_i, G_c, not forward)
        if r.capped or not (G_i > 0.0 and G_c > 0.0):
            return out
        rho_c = float(cb.density(float(state_out.T), float(state_out.P), list(state_out.X)))
        env = cb.rogers_hersh_1976_envelope(
            float(r.Cd) * G_i / G_c, el.plate_thickness / el.d_jet, G_c / rho_c, not forward
        )
        out["Cd_exit_crossflow_extrapolated_velocity_ratio"] = float(env.velocity_ratio)
        out["Cd_exit_crossflow_extrapolated_l_d"] = float(env.l_over_d)
        out["Cd_exit_crossflow_extrapolated_grazing_velocity"] = float(env.grazing_velocity)
        return out
