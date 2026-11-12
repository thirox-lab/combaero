import copy
import dataclasses
import math
import warnings
from abc import ABC, abstractmethod
from collections.abc import Callable, Iterator, MutableMapping
from dataclasses import dataclass, field
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any, Literal, NamedTuple, get_args

from scipy.optimize import brentq as _brentq

import combaero as cb
from combaero import _solver_tools

from .hole_cd import (
    _METERING_SELECTORS,
    HoleCd,
    _hole_cd_spec,
    _HoleCdChain,
    _own_flow_key,
    _resolve_jet_row,
    _with_correlation,
)
from .pressure_loss import ConstantFractionLoss, LinearThetaFractionLoss, PressureLossContext

# Physics configuration types for network introspectability
CompressibilityLiteral = Literal["incompressible", "compressible"]
FrictionModelLiteral = Literal["haaland", "colebrook", "serghides", "petukhov"]
HeatTransferModelLiteral = Literal[
    "none", "gnielinski", "dittus_boelter", "sieder_tate", "petukhov"
]

if TYPE_CHECKING:
    from .graph import FlowNetwork

CombustionMethodLiteral = Literal["complete", "equilibrium"]


# ============================================================================
# Convective Heat Transfer Models and Surface (Phase 3)
# ============================================================================


def _safe_rho(rho_raw: float, rho_min: float = 0.01) -> tuple[float, float]:
    """Smooth lower bound on density for numerical robustness.

    Returns (rho_safe, d_rho_safe / d_rho_raw) for Jacobian chain rule.

    Uses softplus: rho_safe = rho_min + rho_min * ln(1 + exp((rho_raw - rho_min) / rho_min))
    - When rho_raw >> rho_min: rho_safe ~= rho_raw  (transparent)
    - When rho_raw -> 0:        rho_safe -> rho_min * (1 + ln2)  (bounded)
    - When rho_raw -> -inf:     rho_safe -> rho_min  (floor)
    - Derivative is always in (0, 1]: solver always sees a gradient.
    """
    z = (rho_raw - rho_min) / rho_min
    if z > 20.0:  # overflow guard
        return rho_raw, 1.0
    if z < -20.0:  # underflow guard
        return rho_min, 0.0
    exp_z = math.exp(z)
    rho_safe = rho_min + rho_min * math.log1p(exp_z)
    d_safe = exp_z / (1.0 + exp_z)  # sigmoid - always in (0, 1)
    return rho_safe, d_safe


@dataclass
class SmoothModel:
    """Parameters for channel_smooth."""

    correlation: str = "gnielinski"  # "gnielinski" | "dittus_boelter" | "sieder_tate" | "petukhov"
    mu_ratio: float = 1.0  # mu_bulk / mu_wall (Sieder-Tate)
    roughness: float = 0.0  # absolute roughness [m]


@dataclass
class RibbedModel:
    """Rib-roughened walls, evaluated from a provenanced parameter set.

    The correlation gives the RIBBED-SIDE heat transfer. Combining it with the
    smooth walls is this element's job, because how many walls are ribbed is a
    design choice rather than a property of the correlation -- the channel is
    the internal wall by definition, and the other side of it is a different
    channel.

    Parameters
    ----------
    correlation_set : object
        A ``combaero.RibCorrelationSet``. Defaults to Han (1988) for 90 deg
        orthogonal ribs. Supply your own to match measured data: the set
        records its own source and provenance class, so a tuned coefficient
        cannot pass for a published one.
    e_D, p_e, alpha_deg : float
        Rib height / hydraulic diameter, pitch / height, and angle to the flow.
    W_H : float
        Channel width / height. It lives here rather than on the element
        because a channel is defined by its hydraulic diameter, which does not
        determine an aspect ratio -- a 2:1 duct and a square one can share a
        Dh. The correlation needs the ratio for both the wall-law geometry
        group and the ribbed/smooth area split, so the rib model carries it.
    n_ribbed_walls : int
        How many of the four walls carry ribs. 2 means two OPPOSITE walls,
        which is the configuration Han measured.
    smooth_wall_Nu_multiplier : float
        A user knob on the SMOOTH walls' contribution to the channel average.
        It encodes nothing and defaults to 1.0.

        It exists because the ribbed side already has a better knob -- change
        ``C_G`` on the parameter set, which records what you did -- while the
        smooth walls come from the base Gnielinski correlation and have none.

        Its practical use is the documented gap below. A plain smooth wall
        gives ``h_s/h_r`` around 0.42, where Han's own channel average implies
        0.70, because ribs enhance the adjacent smooth wall by 10-50% as well.
        Setting this to about **1.67** reproduces Han's measured average for
        the two-ribbed-wall square-channel case.
    """

    correlation_set: object = None
    e_D: float = 0.0
    p_e: float = 0.0
    alpha_deg: float = 90.0
    W_H: float = 1.0
    n_ribbed_walls: int = 2
    smooth_wall_Nu_multiplier: float = 1.0


@dataclass
class ImpingementModel:
    """Jet array impingement cooling for ONE spanwise row, evaluated from a
    provenanced parameter set (Florschuetz, Truman and Metzger 1981).

    A real array's rows see progressively more crossflow as spent air from
    upstream rows accumulates -- there is no single "channel Nu" for a whole
    array the way there is for a smooth or ribbed duct. Model the array by
    chaining ``n_rows`` separate elements, each with its own ``ImpingementModel``
    and ``row``, the same way a real duct is built from segments elsewhere in
    this codebase. This element does not do that chaining or the row-to-row
    mass-flow bookkeeping for you: see ``row`` below for what it assumes
    instead.

    Parameters
    ----------
    correlation_set : object
        A ``combaero.JetArrayCorrelationSet``. Defaults to
        ``florschuetz_1981_inline()``. Supply ``florschuetz_1981_staggered()``
        for a staggered hole pattern, or your own tuned set.
    d_jet : float
        Jet hole diameter [m]. This row's characteristic length for both
        ``Re_j`` and ``h = Nu * k / d_jet``.
    xn_d, yn_d, z_d : float
        Streamwise hole spacing, spanwise hole spacing, and channel height
        (jet-plate-to-target-plate gap), each normalised by ``d_jet``.
    row : int
        This row's position, 1-indexed counting from upstream. Row 1 sees
        zero crossflow by definition (``Gc/Gj = 0``, matching Florschuetz's
        own ``Nu1``); ``crossflow_to_jet_ratio_at_row`` computes the rest
        from geometry alone, so no other row's state is needed here.
    C_D : float
        Jet-plate discharge coefficient, feeding the ``Gc/Gj`` closed form.
        Defaults to the source's own recommendation absent a measured value
        (``FLORSCHUETZ_1981_DEFAULT_CD = 0.79``) -- combaero has no
        discharge-coefficient correlation of its own for a jet-plate array
        yet, see issue #375.

    What this element assumes about mass flow. The ELEMENT's mass flow is
    this row's jet flow, and ``self.area`` (on the enclosing
    ``ConvectiveSurface``) is THIS ROW's target-plate footprint, which counts
    its holes: ``n = area / (xn_d d_jet * yn_d d_jet)``, so each hole carries
    ``m_dot / n`` and ``Re_j = 4 (m_dot/n) / (pi d_jet mu)`` -- Florschuetz's
    jet mass velocity on the hole area. (Until #460 the flow was taken as
    ``rho v * area``, which made it scale with the target area.) Whether every
    row gets
    the same total flow (a uniform-supply approximation) or a row-dependent
    one (Florschuetz's own Eq. 7, deliberately not implemented -- see
    ``impingement_correlation.h``'s module comment) is up to how the caller
    assembles the chain of elements, not this model.

    What this element does NOT model. The jet plate's own orifice pressure
    loss: model that with a proper ``OrificeElement`` upstream, using a real
    discharge-coefficient correlation -- duplicating it here with a
    hand-rolled formula would risk double-counting it. Friction/pressure
    drop reported here (``f``, ``dP``) is the CROSSFLOW's own plain
    smooth-duct value: Han's book gives no friction correlation for
    impingement the way it does for ribs (there is no ``R``/``f``
    relationship in Eq. 4.9), so this is the best available proxy for the
    spent-air flow along the channel, not a claim that impingement leaves
    friction unchanged.
    """

    correlation_set: object = None
    d_jet: float = 0.0
    xn_d: float = 0.0
    yn_d: float = 0.0
    z_d: float = 0.0
    row: int = 1
    C_D: float = cb.FLORSCHUETZ_1981_DEFAULT_CD
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0


@dataclass
class SingleJetImpingementModel:
    """A single free round jet impinging on a flat plate (Goldstein,
    Behbahani and Heppelmann 1986), at a representative radial position.

    Unlike ``ImpingementModel``, there is no array and no crossflow -- this
    is one jet, one target patch. The correlation's Nu is the AVERAGE over the
    disc of radius ``R`` around the stagnation point (Han: "an average
    heat-transfer coefficient correlation"; Eq. 4.1's ``Nu_bar``), so ``R_D``
    sets the patch it describes. A ``ChannelElement`` left without a
    convective area gets exactly that disc, ``pi (R_D d_jet)^2`` (#462).

    Parameters
    ----------
    correlation_set : object
        A ``combaero.SingleJetImpingementSet``. Defaults to
        ``goldstein_1986_single_jet()``.
    bc : object
        A ``combaero.ImpingementThermalBC``. The correlation's ``(R/D)``
        exponent switches on which surface boundary condition it was fitted
        under; this does not change what the element computes for, only
        which of the source's two curves it reads from.
    d_jet : float
        Jet diameter [m].
    L_D : float
        Jet-to-target-plate spacing / ``d_jet``. Defaults to 7.75, the
        source's own optimum spacing -- a defensible default, not a claim
        about any particular rig.
    R_D : float
        Radius of the averaging disc / ``d_jet``: Nu is Goldstein's average
        over the disc of radius ``R`` around the stagnation point. No
        default: the source's closed-form check point uses 5.0, but that is
        a worked example, not a universal choice.

    Mass flow. The ELEMENT's mass flow is the jet's flow, through its one
    hole: ``Re = 4 m_dot / (pi d_jet mu)``, Goldstein's nozzle Reynolds
    number. ``self.area`` is the target patch the reported h applies to; it
    does not set the flow. (Until #460 it did, as ``rho v * area``.)

    Pressure drop is not modelled here for the same reason as
    ``ImpingementModel``: model the nozzle's own loss with a proper
    ``OrificeElement`` upstream, not a formula duplicated inside this
    heat-transfer model. ``f``/``dP`` reported here are the target-side
    channel's own plain smooth-duct values, borrowed for a T_aw computation
    Han's Eq. 4.2/4.3 recovery-factor model would otherwise have to supply
    -- see ``han_impingement.md`` items 6-7, not yet wired in.
    """

    correlation_set: object = None
    bc: object = cb.ImpingementThermalBC.ConstantHeatFlux
    d_jet: float = 0.0
    L_D: float = 7.75
    R_D: float = 0.0
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0


@dataclass
class PinFinModel:
    """A pin-fin array, evaluated from provenanced parameter sets (#335).

    Heat transfer and friction are separate sets, each in its source's own
    form, evaluated in one canonical basis: ``Re_D`` on the pin diameter and
    the velocity at the minimum flow area, ``Nu_D = h D / k``, and
    ``f = dP / (2 rho Vmax^2 N)``. See ``pin_fin_correlation.h``.

    Parameters
    ----------
    nu_set : object
        A ``combaero.PinFinNuSet`` on the Total surface. Defaults to
        ``metzger_1986_staggered_nu()`` for staggered arrays and
        ``chyu_1998_nu(Inline, Total)`` for inline ones.
    f_set : object
        A ``combaero.PinFinFrictionSet``. Defaults to
        ``metzger_1982_staggered_friction()`` for staggered arrays and
        ``chyu_1990_friction(Inline)`` for inline ones -- the only inline
        friction source in hand, fitted to its digitised Fig. 6 at one
        geometry (H/D 1, S/D = X/D = 2.5). Staggered friction is never
        substituted for inline.
    modifier : object
        An optional ``combaero.PinFinRatioModifier`` applied to ``nu_set``
        (and to ``f_set`` when it carries a friction ratio), e.g.
        ``chyu_1998_inline_over_staggered()`` to transfer a staggered set to
        an inline array. Its ``from_arrangement`` must match ``nu_set``'s and
        its ``to_arrangement`` this array's.
    pin_diameter : float
        Pin diameter D [m].
    S_D, X_D, H_D : float
        Transverse pitch, streamwise pitch and pin height (= channel height),
        each over D.
    N_rows : int
        Pin rows along the flow. Sets the pressure drop and is checked against
        the sets' own row counts (no row correction is applied yet).
    arrangement : object
        ``combaero.PinArrangement``.
    k_pin : float
        Pin thermal conductivity [W/(m K)] for the fin efficiency. 20 is a
        nickel superalloy at turbine metal temperatures; ``math.inf`` treats
        the pins as isothermal.

    Area and velocity conventions. ``ConvectiveSurface.area`` is the BASE
    (planform) area of the endwall this surface couples to; the effective
    coefficient returned is referenced to it,
    ``h_eff = h (A_endwall_exposed + eta_fin A_pin_half) / A_base``, with each
    wall feeding its pins to mid-height. The element's flow area is the
    UNOBSTRUCTED channel cross-section; Vmax follows from it exactly through
    the array geometry.

    The defaults are a working example inside the Metzger sets' box
    (5 mm pins, H/D 1, S/D = X/D = 2.5, 10 rows, Inconel-class pins), not a
    generic "typical" pin fin.
    """

    nu_set: object = None
    f_set: object = None
    modifier: object = None
    pin_diameter: float = 0.005
    S_D: float = 2.5
    X_D: float = 2.5
    H_D: float = 1.0
    N_rows: int = 10
    arrangement: object = cb.PinArrangement.Staggered
    k_pin: float = 20.0

    def __post_init__(self) -> None:
        if not (self.pin_diameter > 0.0):
            raise ValueError("PinFinModel: pin_diameter must be positive [m]")
        cb.validate_pin_fin_geometry(self.geometry())
        f_arr = self.resolved_f_set().arrangement
        if f_arr != self.arrangement:
            raise ValueError(
                f"PinFinModel: f_set {self.resolved_f_set().name!r} is "
                f"{f_arr.name}, but this array is {self.arrangement.name}. "
                "Staggered and inline friction differ by about 1.6x (Chyu "
                "1990, Fig. 6) and are not substituted for each other."
            )
        if self.modifier is not None:
            nu_arr = self.resolved_nu_set().arrangement
            if self.modifier.from_arrangement != nu_arr:
                raise ValueError(
                    f"PinFinModel: modifier {self.modifier.name!r} transfers from "
                    f"{self.modifier.from_arrangement.name}, but nu_set is "
                    f"{nu_arr.name}."
                )
            if self.modifier.to_arrangement != self.arrangement:
                raise ValueError(
                    f"PinFinModel: modifier {self.modifier.name!r} transfers to "
                    f"{self.modifier.to_arrangement.name}, but this array is "
                    f"{self.arrangement.name}."
                )

    def geometry(self):
        return cb.PinFinGeometry(
            S_D=self.S_D,
            X_D=self.X_D,
            H_D=self.H_D,
            N_rows=self.N_rows,
            arrangement=self.arrangement,
        )

    def resolved_nu_set(self):
        if self.nu_set is not None:
            return self.nu_set
        if self.modifier is not None:
            # A transfer modifier starts from its own source arrangement.
            if self.modifier.from_arrangement == cb.PinArrangement.Staggered:
                return cb.metzger_1986_staggered_nu()
            return cb.chyu_1998_nu(self.modifier.from_arrangement, cb.PinNuSurface.Total)
        if self.arrangement == cb.PinArrangement.Inline:
            return cb.chyu_1998_nu(cb.PinArrangement.Inline, cb.PinNuSurface.Total)
        return cb.metzger_1986_staggered_nu()

    def resolved_f_set(self):
        if self.f_set is not None:
            return self.f_set
        if self.arrangement == cb.PinArrangement.Inline:
            return cb.chyu_1990_friction(cb.PinArrangement.Inline)
        return cb.metzger_1982_staggered_friction()


@dataclass
class _ImpingementChannelResult:
    """What either impingement model returns.

    Deliberately NOT a ``ChannelResult``: same rationale as
    ``_RibbedChannelResult`` -- this path computes its own derivatives and
    (for the array case) an ``extrapolated`` flag the C++ correlation
    already carries.
    """

    h: float
    Nu: float
    Re: float
    Pr: float
    f: float
    dP: float
    T_aw: float
    extrapolated: bool = False
    dh_dmdot: float = 0.0
    dh_dT: float = 0.0
    dT_aw_dmdot: float = 0.0
    dT_aw_dT: float = 0.0
    dT_aw_dP: float = 0.0


@dataclass
class _PinFinChannelResult:
    """What a pin-fin array returns.

    ``h`` is the effective coefficient on the BASE area (fin efficiency
    included); ``h_array`` is the correlation's own coefficient on the total
    pin + endwall area, before the fin model, so both halves of that
    modelling step are visible.
    """

    h: float
    h_array: float
    eta_fin: float
    Nu: float
    Re: float
    Pr: float
    f: float
    dP: float
    T_aw: float
    extrapolated: bool = False
    dh_dmdot: float = 0.0
    dh_dT: float = 0.0
    dT_aw_dmdot: float = 0.0
    dT_aw_dT: float = 0.0
    dT_aw_dP: float = 0.0


@dataclass
class _RibbedChannelResult:
    """What a ribbed channel returns.

    Deliberately NOT a ``ChannelResult``: that type is the C++ correlation's
    output and carries derivative fields this path computes differently. It
    also exposes ``h_ribbed`` and ``h_smooth`` separately, because the channel
    average hides a 20% modelling choice and a user should be able to see both
    halves of it.
    """

    h: float
    h_ribbed: float
    h_smooth: float
    Nu: float
    Re: float
    Pr: float
    f: float
    dP: float
    T_aw: float
    e_plus: float
    extrapolated: bool
    # Derivatives the solver's wall-coupling path reads. They are not
    # decoration: without them a ribbed channel joined by a ThermalWall raises
    # AttributeError mid-solve, which is how the coupled combustor example
    # found this after the unit tests missed it -- none of them coupled a wall.
    dh_dmdot: float = 0.0
    dh_dT: float = 0.0
    dT_aw_dmdot: float = 0.0
    dT_aw_dT: float = 0.0
    dT_aw_dP: float = 0.0


ChannelModel = (
    SmoothModel | RibbedModel | ImpingementModel | SingleJetImpingementModel | PinFinModel
)


def _stanton(res: Any) -> float:
    """St = Nu / (Re Pr) of a heat-transfer result, at the correlation's own Re
    and Pr (the jet Re for impingement, the pin Re for pin fins), so it is the
    Stanton number of the h it reports (#895). 0 without a defined Re Pr."""
    if res is None:
        return 0.0
    re = float(getattr(res, "Re", 0.0) or 0.0)
    pr = float(getattr(res, "Pr", 0.0) or 0.0)
    nu = float(getattr(res, "Nu", 0.0) or 0.0)
    return nu / (re * pr) if re > 0.0 and pr > 0.0 and math.isfinite(nu) else 0.0


def _max_over_min(flows: list[float]) -> float:
    """max / min of a set of row flows, the nonuniformity an array reports.
    inf when a row is starved while another flows (unbounded, but defined);
    NaN when none flows or any runs backwards (no ratio of the two means
    anything) -- a key with no value is NaN, never absent."""
    lo, hi = min(flows), max(flows)
    if lo < 0.0 or hi <= 0.0:
        return math.nan
    return float(hi / lo) if lo > 0.0 else math.inf


@dataclass
class ConvectiveSurface:
    """Convective surface description attached to a NetworkElement.

    Attributes
    ----------
    area : float
        Wetted convective area [m^2]. Default 0.0 (disabled).
    model : ChannelModel
        Model-specific subclass holding geometry parameters.
    heating : bool | None
        True if fluid is being heated, False if cooled, None = auto-detect
        from sign of (T_hot - T_fluid).
    Nu_multiplier : float
        Empirical correction factor on Nusselt number (default 1.0).
    f_multiplier : float
        Empirical correction factor on friction factor (default 1.0).
    """

    area: float = 0.0  # A_conv [m^2] - 0 disables
    model: ChannelModel = field(default_factory=SmoothModel)
    heating: bool | None = None  # None = auto-detect
    Nu_multiplier: float = 1.0  # empirical correction on Nu
    f_multiplier: float = 1.0  # empirical correction on f
    # The calling element's flow area for the current evaluation (see
    # htc_and_T); not a parameter of the surface.
    _flow_area: float | None = field(default=None, init=False, repr=False, compare=False)

    def _channel_mdot(self, rho: float, velocity: float, diameter: float) -> float:
        """The solver's mass flow for this surface: through the CHANNEL
        cross-section -- the element's own flow area when it passed one
        (``flow_area`` on ``htc_and_T``), else ``pi Dh^2 / 4`` as
        ``channel_smooth`` assumes. A non-circular channel given a separate
        ``Dh`` needs the real area: ``pi Dh^2/4`` is then off by
        ``(Dh/diameter)^2`` (#460).

        Every ``dh_dmdot`` this class returns is a derivative with respect to
        the element's own ``m_dot`` -- that is what the solver's wall coupling
        multiplies it by -- and ``ChannelElement.htc_and_T`` derives the
        velocity from that mass flow and the channel area, the same convention
        ``channel_smooth`` uses. The convective ``area`` is a heat-transfer
        area, not a flow area. Using it here scaled the ribbed and impingement
        derivatives by ``A_cross / A_surface`` (#456).
        """
        area = self._flow_area if self._flow_area else math.pi / 4.0 * diameter * diameter
        return rho * abs(velocity) * area

    def _flow_area_arg(self) -> float:
        """The element's flow area for channel_smooth's mass-flow derivatives,
        NaN when it passed none (channel_smooth then assumes pi D^2/4, #463)."""
        return self._flow_area if getattr(self, "_flow_area", None) else math.nan

    def _ribbed_result(self, T, P, X, velocity, diameter, length, T_hot, heating):
        """Ribbed-channel heat transfer and pressure drop.

        Two things are asymmetric here and both follow from the source rather
        than from convenience.

        FRICTION NEEDS NO WALL WEIGHTING. The ``f`` in the roughness function's
        definition is already the equivalent four-sided channel friction
        factor, so the correlation returns a channel-level value directly. The
        element does not multiply pipe friction by anything -- correlations own
        their ``f``, which is what #331 established after a round-trip through
        a restated friction factor moved a drop by -24.7%.

        HEAT TRANSFER DOES. The correlation gives the ribbed side; the smooth
        walls come from the base correlation, and the channel average is the
        area-weighted combination.

        A DOCUMENTED GAP. Using the plain smooth correlation for the smooth
        walls gives h_s/h_r around 0.42, where Han's own reported channel
        average implies 0.70 -- because ribs enhance the adjacent smooth wall
        by 10-50% as well, which no correlation here covers. The channel
        average is therefore about 20% below Han's measurement for the
        two-ribbed-wall square case. That is a knowable, quantified
        under-prediction rather than an invented constant, and
        ``smooth_wall_Nu_multiplier`` is how a user closes it: about 1.67
        reproduces Han.
        """
        model = self.model
        rib_set = model.correlation_set or cb.han_1988_orthogonal()

        rho, _ = _safe_rho(cb.density(T, P, X))
        cs = cb.complete_state(T, P, X)
        mu = cs.transport.mu
        Re = rho * velocity * diameter / mu if mu > 0.0 else 0.0
        Pr = cs.transport.Pr

        geom = cb.RibGeometry(
            e_D=model.e_D,
            p_e=model.p_e,
            W_H=model.W_H,
            alpha_deg=model.alpha_deg,
        )
        rib = cb.evaluate_rib(rib_set, geom, Re)

        # Ribbed side, from the correlation's Stanton number.
        k = cs.transport.k
        cp = cs.thermo.cp
        h_ribbed = rib.St_r * rho * abs(velocity) * cp

        # Smooth walls, from the base correlation, with the user's knob.
        smooth = cb.channel_smooth(
            T,
            P,
            X,
            velocity,
            diameter,
            length,
            T_hot=T_hot,
            heating=heating,
            Nu_multiplier=model.smooth_wall_Nu_multiplier,
            f_multiplier=1.0,
            flow_area=self._flow_area_arg(),
        )

        frac_ribbed, frac_smooth = self._ribbed_wall_fractions()
        h_avg = frac_ribbed * h_ribbed + frac_smooth * smooth.h

        # The correlation's f is already the four-sided channel value.
        dP = rib.f * (length / diameter) * 0.5 * rho * velocity * abs(velocity)

        # Wall-coupling derivatives. The ribbed side is h_r = St(Re) rho|v| cp,
        # and BOTH Re and |v| are proportional to the mass flow, so
        #     dh_r/dmdot = (Re dSt/dRe + St) rho |v| cp / mdot.
        # The St term dominates (h_r grows roughly as mdot^0.8); it was missing,
        # and the mass flow was taken on the convective area (#456). The
        # smooth side brings its own derivative from the base correlation;
        # both are area-weighted exactly as h is.
        mdot = self._channel_mdot(rho, velocity, diameter)
        dh_ribbed_dmdot = (
            (Re * rib.dSt_dRe + rib.St_r) * rho * abs(velocity) * cp / mdot if mdot else 0.0
        )
        dh_dmdot = (
            frac_ribbed * dh_ribbed_dmdot + frac_smooth * smooth.dh_dmdot
        ) * self.Nu_multiplier
        # Temperature sensitivity of the ribbed side is not exposed by the
        # correlation, so only the smooth side contributes. Recorded rather
        # than approximated from a stand-in, the same gap the array path had.
        dh_dT = frac_smooth * smooth.dh_dT * self.Nu_multiplier

        return _RibbedChannelResult(
            dh_dmdot=dh_dmdot,
            dh_dT=dh_dT,
            dT_aw_dmdot=smooth.dT_aw_dmdot,
            dT_aw_dT=smooth.dT_aw_dT,
            dT_aw_dP=smooth.dT_aw_dP,
            h=h_avg * self.Nu_multiplier,
            h_ribbed=h_ribbed,
            h_smooth=smooth.h,
            Nu=h_avg * diameter / k if k > 0.0 else 0.0,
            Re=Re,
            Pr=Pr,
            f=rib.f * self.f_multiplier,
            dP=dP * self.f_multiplier,
            T_aw=smooth.T_aw,
            e_plus=rib.e_plus,
            extrapolated=rib.extrapolated,
        )

    def _ribbed_wall_fractions(self) -> tuple[float, float]:
        """Area fractions of the ribbed and smooth walls, for a W x H duct.

        Han's configuration is two OPPOSITE walls, which for a duct of width W
        and height H are the two of width W. The ribbed fraction is then
        W/(W+H) -- the same area weighting the four-sided friction conversion
        uses, so friction and heat transfer are treated alike.

        One and four ribbed walls follow from the same picture. Anything else
        is rejected rather than interpolated: "three ribbed walls" has no
        unambiguous geometry.
        """
        model = self.model
        W_H = model.W_H
        n = model.n_ribbed_walls
        if n == 4:
            return 1.0, 0.0
        if n == 2:
            ribbed = W_H / (W_H + 1.0)
            return ribbed, 1.0 - ribbed
        if n == 1:
            ribbed = W_H / (2.0 * (W_H + 1.0))
            return ribbed, 1.0 - ribbed
        raise ValueError(
            f"n_ribbed_walls must be 1, 2 or 4, got {n}. 2 means two opposite "
            "walls, which is the configuration the correlations were measured "
            "on."
        )

    def _impingement_result(self, T, P, X, velocity, diameter, length, T_hot, heating):
        """Jet array impingement, one row (Florschuetz, Truman and Metzger 1981).

        See ``ImpingementModel``'s docstring for the mass-flow and
        pressure-drop conventions this follows.
        """
        model = self.model
        jet_set = model.correlation_set or cb.florschuetz_1981_inline()

        rho, _ = _safe_rho(cb.density(T, P, X))
        cs = cb.complete_state(T, P, X)
        mu = cs.transport.mu
        k = cs.transport.k
        Pr = cs.transport.Pr

        # The element's mass flow IS this row's jet flow (Florschuetz: G_j is
        # the jet mass velocity on the hole area). The convective area only
        # counts the holes. It used to supply the flow as well, as
        # rho v * area -- wrong by A_surface/A_cross (#460).
        mdot_total = self._channel_mdot(rho, velocity, diameter)
        hole_footprint = model.xn_d * model.d_jet * model.yn_d * model.d_jet
        n_holes = self.area / hole_footprint if hole_footprint > 0.0 else 0.0
        mdot_per_hole = mdot_total / n_holes if n_holes > 0.0 else 0.0
        hole_area = math.pi / 4.0 * model.d_jet**2
        v_jet = mdot_per_hole / (rho * hole_area) if rho * hole_area > 0.0 else 0.0
        Re_j = rho * v_jet * model.d_jet / mu if mu > 0.0 else 0.0

        Gc_Gj = cb.crossflow_to_jet_ratio_at_row(model.yn_d, model.z_d, model.C_D, model.row)
        jet = cb.jet_array_impingement_nu(
            jet_set, Re_j, Gc_Gj, Pr, model.xn_d, model.yn_d, model.z_d
        )
        h = jet.Nu * k / model.d_jet if model.d_jet > 0.0 else 0.0

        # T_aw, f and dP are not covered by Eq. 4.9 at all -- Han's book gives
        # no impingement-specific friction/recovery model, unlike ribs' own
        # R/f relationship. Borrowed wholesale from the plain smooth
        # correlation applied to the crossflow's own bulk state, per
        # ImpingementModel's docstring.
        smooth = cb.channel_smooth(
            T,
            P,
            X,
            velocity,
            diameter,
            length,
            T_hot=T_hot,
            heating=heating,
            Nu_multiplier=1.0,
            f_multiplier=1.0,
            flow_area=self._flow_area_arg(),
        )

        # Re_j is proportional to the element's mass flow; the per-hole split
        # above is this row's own bookkeeping, not the solver's m_dot (#456).
        mdot = self._channel_mdot(rho, velocity, diameter)
        dRe_j_dmdot = abs(Re_j / mdot) if mdot else 0.0
        dh_dmdot = jet.dNu_dRe_j * dRe_j_dmdot * (k / model.d_jet if model.d_jet > 0.0 else 0.0)
        dh_dmdot *= self.Nu_multiplier
        # Temperature sensitivity of the correlation itself is not exposed;
        # only the borrowed smooth-side contributes, the same documented gap
        # ribs accepted for their own T-sensitivity.
        dh_dT = smooth.dh_dT * self.Nu_multiplier

        return _ImpingementChannelResult(
            h=h * self.Nu_multiplier,
            Nu=jet.Nu,
            Re=Re_j,
            Pr=Pr,
            f=smooth.f * self.f_multiplier,
            dP=smooth.dP * self.f_multiplier,
            T_aw=smooth.T_aw,
            extrapolated=jet.extrapolated,
            dh_dmdot=dh_dmdot,
            dh_dT=dh_dT,
            dT_aw_dmdot=smooth.dT_aw_dmdot,
            dT_aw_dT=smooth.dT_aw_dT,
            dT_aw_dP=smooth.dT_aw_dP,
        )

    def _single_jet_impingement_result(self, T, P, X, velocity, diameter, length, T_hot, heating):
        """A single free jet (Goldstein, Behbahani and Heppelmann 1986).

        See ``SingleJetImpingementModel``'s docstring for the mass-flow and
        pressure-drop conventions this follows.
        """
        model = self.model
        jet_set = model.correlation_set or cb.goldstein_1986_single_jet()

        rho, _ = _safe_rho(cb.density(T, P, X))
        cs = cb.complete_state(T, P, X)
        mu = cs.transport.mu
        k = cs.transport.k

        # The element's mass flow IS the jet's flow (Goldstein: nozzle Re on
        # the jet's own exit area), not rho v * the target-patch area (#460).
        mdot_total = self._channel_mdot(rho, velocity, diameter)
        hole_area = math.pi / 4.0 * model.d_jet**2
        v_jet = mdot_total / (rho * hole_area) if rho * hole_area > 0.0 else 0.0
        Re = rho * v_jet * model.d_jet / mu if mu > 0.0 else 0.0

        jet = cb.single_jet_impingement(jet_set, model.bc, Re, model.L_D, model.R_D)
        h = jet.Nu * k / model.d_jet if model.d_jet > 0.0 else 0.0

        smooth = cb.channel_smooth(
            T,
            P,
            X,
            velocity,
            diameter,
            length,
            T_hot=T_hot,
            heating=heating,
            Nu_multiplier=1.0,
            f_multiplier=1.0,
            flow_area=self._flow_area_arg(),
        )

        mdot = self._channel_mdot(rho, velocity, diameter)
        dRe_dmdot = abs(Re / mdot) if mdot else 0.0
        dh_dmdot = jet.dNu_dRe * dRe_dmdot * (k / model.d_jet if model.d_jet > 0.0 else 0.0)
        dh_dmdot *= self.Nu_multiplier
        dh_dT = smooth.dh_dT * self.Nu_multiplier

        return _ImpingementChannelResult(
            h=h * self.Nu_multiplier,
            Nu=jet.Nu,
            Re=Re,
            Pr=cs.transport.Pr,
            f=smooth.f * self.f_multiplier,
            dP=smooth.dP * self.f_multiplier,
            T_aw=smooth.T_aw,
            dh_dmdot=dh_dmdot,
            dh_dT=dh_dT,
            dT_aw_dmdot=smooth.dT_aw_dmdot,
            dT_aw_dT=smooth.dT_aw_dT,
            dT_aw_dP=smooth.dT_aw_dP,
        )

    def _pin_fin_terms(self, Re: float, Pr: float):
        """Nu and f (with modifier) at canonical Re_D, and their Re-slopes."""
        model = self.model
        geom = model.geometry()
        nu = cb.evaluate_pin_fin_nu(model.resolved_nu_set(), geom, Re, Pr)
        fr = cb.evaluate_pin_fin_friction(model.resolved_f_set(), geom, Re)
        Nu, dNu = nu.Nu, nu.dNu_dRe
        f, df = fr.f, fr.df_dRe
        extrapolated = nu.extrapolated or fr.extrapolated
        if model.modifier is not None:
            m = cb.evaluate_pin_fin_modifier(model.modifier, geom, Re)
            Nu, dNu = Nu * m.ratio_Nu, dNu * m.ratio_Nu + Nu * m.dratio_Nu_dRe
            if m.has_f:
                f, df = f * m.ratio_f, df * m.ratio_f + f * m.dratio_f_dRe
            extrapolated = extrapolated or m.extrapolated
        return Nu, dNu, f, df, extrapolated

    def _pin_fin_result(self, T, P, X, velocity, diameter, length, T_hot, heating):
        """Pin-fin array: Metzger et al. (and others) via pin_fin_correlation.h.

        ``velocity`` is the bulk velocity through the UNOBSTRUCTED channel
        cross-section; Vmax = velocity / (A_min/A_frontal), exactly. The
        correlation's h lives on the total pin + endwall area; the fin model
        turns it into ``h`` on this surface's BASE area:

            h_eff = h (f_endwall + eta_fin(h) f_pin)

        with the per-wall fractions of ``pin_fin_area_fractions``.

        ``dh_dmdot`` uses the element's mass flow through the channel
        cross-section, ``rho v pi Dh^2 / 4`` -- the same convention as
        ``channel_smooth``. T_aw and its derivatives are borrowed from the
        smooth correlation, as the impingement paths do; ``dh_dT`` is the
        smooth result's relative property sensitivity scaled to ``h_eff``, a
        documented proxy: the pin correlations expose no temperature
        derivative of their own.
        """
        model = self.model
        D = model.pin_diameter
        geom = model.geometry()

        rho, _ = _safe_rho(cb.density(T, P, X))
        cs = cb.complete_state(T, P, X)
        mu = cs.transport.mu
        k = cs.transport.k
        Pr = cs.transport.Pr

        amin = cb.pin_fin_amin_over_afrontal(geom)
        v_max = velocity / amin
        Re = rho * v_max * D / mu if mu > 0.0 else 0.0

        Nu, dNu, f, _, extrapolated = self._pin_fin_terms(Re, Pr)
        h_array = Nu * k / D
        dh_array_dRe = dNu * k / D

        frac = cb.pin_fin_area_fractions(geom)
        eff = cb.pin_fin_array_efficiency(
            h_array, model.k_pin, D, model.H_D * D, frac.pin_over_total
        )
        h_eff = h_array * (frac.endwall_exposed + eff.eta_fin * frac.pin)
        dh_eff_dh = frac.endwall_exposed + frac.pin * (eff.eta_fin + h_array * eff.deta_fin_dh)

        smooth = cb.channel_smooth(
            T,
            P,
            X,
            velocity,
            diameter,
            length,
            T_hot=T_hot,
            heating=heating,
            Nu_multiplier=1.0,
            f_multiplier=1.0,
            flow_area=self._flow_area_arg(),
        )

        mdot = self._channel_mdot(rho, velocity, diameter)
        dRe_dmdot = abs(Re / mdot) if mdot else 0.0
        dh_dmdot = dh_eff_dh * dh_array_dRe * dRe_dmdot * self.Nu_multiplier
        dh_dT = smooth.dh_dT * (h_eff / smooth.h) * self.Nu_multiplier if smooth.h > 0.0 else 0.0

        dP = 2.0 * rho * v_max * abs(v_max) * model.N_rows * f

        return _PinFinChannelResult(
            h=h_eff * self.Nu_multiplier,
            h_array=h_array,
            eta_fin=eff.eta_fin,
            Nu=Nu,
            Re=Re,
            Pr=Pr,
            f=f * self.f_multiplier,
            dP=dP * self.f_multiplier,
            T_aw=smooth.T_aw,
            extrapolated=extrapolated,
            dh_dmdot=dh_dmdot,
            dh_dT=dh_dT,
            dT_aw_dmdot=smooth.dT_aw_dmdot,
            dT_aw_dT=smooth.dT_aw_dT,
            dT_aw_dP=smooth.dT_aw_dP,
        )

    def htc_and_T(
        self,
        T: float,
        P: float,
        X: list[float],
        velocity: float,
        diameter: float,
        length: float,
        T_hot: float = math.nan,
        flow_area: float | None = None,
    ):
        """Compute heat transfer coefficient and adiabatic wall temperature.

        Parameters
        ----------
        T : float
            Bulk static temperature [K].
        P : float
            Bulk static pressure [Pa].
        X : list[float]
            Mole fractions [-].
        velocity : float
            Bulk flow velocity [m/s].
        diameter : float
            Hydraulic diameter [m].
        length : float
            Channel length [m].
        T_hot : float, optional
            Wall temperature [K]. Used for auto-detection of heating/cooling.
        flow_area : float, optional
            The element's flow cross-section [m^2], from which it computed
            ``velocity``. Lets the surface recover the element's mass flow
            exactly; without it ``pi diameter^2 / 4`` is assumed.

        Returns
        -------
        ChannelResult | None
            Full ChannelResult with h, T_aw, and Jacobians (dh_dmdot, dh_dT, etc.),
            or None if area=0. Access convective area via ``self.area``.
        """

        if self.area == 0.0 or abs(velocity) < 1e-12:
            return None
        self._flow_area = flow_area

        # Auto-detect heating direction from T_hot - T sign
        if self.heating is not None:
            heating = self.heating
        elif math.isfinite(T_hot):
            heating = T_hot >= T
        else:
            heating = True  # default when T_hot unknown

        # Dispatch to appropriate C++ channel_* function based on model type
        if isinstance(self.model, SmoothModel):
            result = cb.channel_smooth(
                T,
                P,
                X,
                velocity,
                diameter,
                length,
                T_hot=T_hot,
                correlation=self.model.correlation,
                heating=heating,
                mu_ratio=self.model.mu_ratio,
                roughness=self.model.roughness,
                Nu_multiplier=self.Nu_multiplier,
                f_multiplier=self.f_multiplier,
                flow_area=self._flow_area_arg(),
            )
        elif isinstance(self.model, RibbedModel):
            result = self._ribbed_result(T, P, X, velocity, diameter, length, T_hot, heating)
        elif isinstance(self.model, ImpingementModel):
            result = self._impingement_result(T, P, X, velocity, diameter, length, T_hot, heating)
        elif isinstance(self.model, SingleJetImpingementModel):
            result = self._single_jet_impingement_result(
                T, P, X, velocity, diameter, length, T_hot, heating
            )
        elif isinstance(self.model, PinFinModel):
            result = self._pin_fin_result(T, P, X, velocity, diameter, length, T_hot, heating)
        else:
            raise TypeError(
                f"Unsupported channel model {type(self.model).__name__}. "
                "Enhanced-surface correlations were removed in 0.7.0 pending "
                "provenanced replacements; see issue #339."
            )

        return result


class _DuctHtc:
    """A surface result evaluated at a duct's own static state, its partials
    chained onto what the solver differentiates: the feed node's stagnation
    T and Pt, and the element's flow (#636). Every other attribute is the
    surface result's."""

    __slots__ = (
        "_r",
        "dh_dT",
        "dh_dmdot",
        "dh_dPt",
        "dT_aw_dT",
        "dT_aw_dmdot",
        "dT_aw_dP",
        "dT_aw_dPt",
        "_ctx",
        "_dY",
    )

    def __init__(
        self,
        r: Any,
        s: "cb.DuctStaticState",
        area: float,
        ctx: tuple[float, float, float, list[float], float, "cb.DuctStaticState"] | None = None,
    ) -> None:
        self._r = r
        # (T0, P0, G, X) the static state was taken at: the composition
        # partials need its gradient too (#717), computed on first use.
        self._ctx = ctx
        self._dY: tuple[list[float], list[float]] | None = None
        # The surface's partials are at fixed m, in the static (T, P) it was
        # given; T_s and P_s move with T0, P0 and G = |m| / A.
        dh_dT, dTaw_dT = float(r.dh_dT), float(r.dT_aw_dT)
        dTaw_dP = float(getattr(r, "dT_aw_dP", 0.0))
        # The surface's own flow partials are per unit of the flow it was
        # handed, area * G_eff, which moves with the true flow (dGeff_dG) and,
        # near choke, with T0 and P0 through the sonic flux (#551).
        dh_m, dTaw_m = float(r.dh_dmdot), float(r.dT_aw_dmdot)
        self.dh_dT = dh_dT * s.dT_dTt + dh_m * area * s.dGeff_dTt
        self.dh_dmdot = dh_m * s.dGeff_dG + dh_dT * s.dT_dG / area
        self.dh_dPt = dh_dT * s.dT_dPt + dh_m * area * s.dGeff_dPt
        self.dT_aw_dT = dTaw_dT * s.dT_dTt + dTaw_dP * s.dP_dTt + dTaw_m * area * s.dGeff_dTt
        self.dT_aw_dmdot = dTaw_m * s.dGeff_dG + (dTaw_dT * s.dT_dG + dTaw_dP * s.dP_dG) / area
        self.dT_aw_dPt = dTaw_dT * s.dT_dPt + dTaw_dP * s.dP_dPt + dTaw_m * area * s.dGeff_dPt
        # The static P is the duct's, a function of the feed's Pt: the feed
        # node's own static P does not enter.
        self.dT_aw_dP = 0.0

    def _composition(self) -> tuple[list[float], list[float]]:
        """(dh/dY, dT_aw/dY) of the feed's mass fractions at fixed T0, P0 and
        flow (#717): the surface's own partials at fixed static state, plus
        the static T, P and the bent flux moving with the composition. Empty
        when the surface provides none (read as zero)."""
        if self._dY is None:
            r = self._r
            dh_Y = list(getattr(r, "dh_dY", None) or [])
            dTaw_Y = list(getattr(r, "dT_aw_dY", None) or [])
            if not dh_Y or self._ctx is None:
                self._dY = ([], [])
            else:
                T0, P0, G, X, area, st = self._ctx
                sd = cb.duct_htc_static_state_dY(T0, P0, G, X, st)
                dh_T, dTaw_T = float(r.dh_dT), float(r.dT_aw_dT)
                dTaw_P = float(getattr(r, "dT_aw_dP", 0.0))
                dh_m, dTaw_m = float(r.dh_dmdot), float(r.dT_aw_dmdot)
                self._dY = (
                    [
                        a + dh_T * t + dh_m * area * g
                        for a, t, g in zip(dh_Y, sd.dT_dY, sd.dGeff_dY, strict=True)
                    ],
                    [
                        a + dTaw_T * t + dTaw_P * q + dTaw_m * area * g
                        for a, t, q, g in zip(dTaw_Y, sd.dT_dY, sd.dP_dY, sd.dGeff_dY, strict=True)
                    ],
                )
        return self._dY

    @property
    def dh_dY(self) -> list[float]:
        return self._composition()[0]

    @property
    def dT_aw_dY(self) -> list[float]:
        return self._composition()[1]

    def __getattr__(self, name: str) -> Any:
        return getattr(self._r, name)

    def __dir__(self) -> list[str]:
        # The pass-through fields too (dP, q and their partials, unchained:
        # no solver path reads them from a heat-transfer result).
        return sorted(set(dir(self._r)) | set(self.__slots__))


def _duct_htc(
    surface: "ConvectiveSurface",
    state: "NetworkMixtureState",
    area: float,
    diameter: float,
    length: float,
    T_hot: float = math.nan,
) -> "_DuctHtc | None":
    """The surface correlation at a duct's inlet static state (#636).

    ``state`` is the feed node's: the function reads its ``Tt`` and ``Pt``
    (a plenum's gas is at rest, so there Pt = P; a chamber carries its
    stagnation ``Tt``, #578). The duct's static state comes from them and the
    mass flux by the isentropic acceleration (C++ ``duct_static_state``).
    Handing the node's Tt to the correlation as static put T_aw above the
    stagnation temperature, by
    about (1 + r)(Tt - T): 2.6 K at Mach 0.1, 22 K at 0.3, 51 K at 0.5 for a
    1600 K stream.
    """
    G = abs(float(state.m_dot)) / area
    # The flux bent C1 just below the sonic one (#551): a duct choked at its
    # inlet sits at the sonic flux, where the isentropic state has an
    # infinite slope in G and the wall coupling could not converge. Below
    # 99.8% of the sonic flux this is the plain duct state.
    s = cb.duct_htc_static_state(float(state.Tt), float(state.Pt), G, state.X)
    # Past choke (a solver iterate beyond the duct's capacity) the velocity
    # is the bent flux's too: G / rho at the iterate's G sent T_aw, and the
    # heat, without bound (53 caught mixer failures in the stress sweep, 43%
    # of one network's duct evaluations past choke).
    velocity = s.G_eff / s.rho
    r = surface.htc_and_T(
        T=s.T,
        P=s.P,
        X=state.X,
        velocity=velocity,
        diameter=diameter,
        length=length,
        T_hot=T_hot,
        flow_area=area,
    )
    ctx = (float(state.Tt), float(state.Pt), G, list(state.X), area, s)
    return None if r is None else _DuctHtc(r, s, area, ctx)


class _ChamberHtc:
    """A surface result evaluated at a momentum chamber's static state, its
    partials chained onto the node's own T0 and static P and its throughput:
    the closure's T_static moves with all three. Every other attribute is
    the surface result's."""

    __slots__ = (
        "_r",
        "dh_dT",
        "dh_dmdot",
        "dh_dP",
        "dT_aw_dT",
        "dT_aw_dmdot",
        "dT_aw_dP",
        "dh_dY",
        "dT_aw_dY",
    )

    def __init__(self, r: Any, c: Any) -> None:
        self._r = r
        dh_dT, dTaw_dT = float(r.dh_dT), float(r.dT_aw_dT)
        self.dh_dT = dh_dT * c.dT_static_dT
        self.dh_dmdot = float(r.dh_dmdot) + dh_dT * c.dT_static_dmdot
        self.dh_dP = dh_dT * c.dT_static_dP
        self.dT_aw_dT = dTaw_dT * c.dT_static_dT
        self.dT_aw_dmdot = float(r.dT_aw_dmdot) + dTaw_dT * c.dT_static_dmdot
        self.dT_aw_dP = float(getattr(r, "dT_aw_dP", 0.0)) + dTaw_dT * c.dT_static_dP
        # Composition (#743): the surface's own, plus the static T moving
        # with it. Empty when the surface provides none.
        dh_Y = list(getattr(r, "dh_dY", None) or [])
        dTaw_Y = list(getattr(r, "dT_aw_dY", None) or [])
        dTs_Y = list(getattr(c, "dT_static_dY", None) or [0.0] * len(dh_Y))
        self.dh_dY = [a + dh_dT * t for a, t in zip(dh_Y, dTs_Y, strict=False)]
        self.dT_aw_dY = [a + dTaw_dT * t for a, t in zip(dTaw_Y, dTs_Y, strict=False)]

    def __getattr__(self, name: str) -> Any:
        return getattr(self._r, name)

    def __dir__(self) -> list[str]:
        return sorted(set(dir(self._r)) | set(self.__slots__))


def _summed_inflow_jac(
    upstream_states: list["NetworkMixtureState"],
) -> tuple[list[str], dict[str, float]]:
    """(upstream element ids, d(total inflow)/d(unknown name)) of a node.

    Each upstream state carries d(its flow into the node)/d(unknown) as
    ``_m_dot_jac_names``; a multi-source element (tee) hands the node one
    state per source, so its entries are counted once per element.
    """
    ids: list[str] = []
    jac: dict[str, float] = {}
    seen: set[str] = set()
    for s in upstream_states:
        if hasattr(s, "_element_id"):
            ids.append(s._element_id)
            eid = s._element_id
            if eid not in seen and hasattr(s, "_m_dot_jac_names"):
                seen.add(eid)
                for var, coeff in s._m_dot_jac_names.items():
                    jac[var] = jac.get(var, 0.0) + coeff
    return ids, jac


# ============================================================================
# WallConnection (Phase 4/5)
# ============================================================================


@dataclass
class WallLayer:
    """Description of a single physical layer in a multi-layer ThermalWall.

    Attributes
    ----------
    thickness : float
        Thickness of the layer [m].
    conductivity : float
        Thermal conductivity of the material [W/(m*K)].
    material : str
        Material name (e.g., 'inconel718', 'haynes230', or 'custom').
    extrapolate : bool
        Outside the material's sourced range: False (default) follows the k(T)
        envelope (#882), True continues the end slope linearly. Neither is
        sourced; ``in_range`` says whether the last evaluation was.
    in_range : bool
        Whether the layer's last conductivity came from its material's sourced
        range (always True for a 'custom' layer). Set by update_conductivity.
    """

    thickness: float
    conductivity: float
    material: str = "custom"
    extrapolate: bool = False
    in_range: bool = field(default=True, init=False, compare=False)

    def __post_init__(self) -> None:
        # An unknown name used to be swallowed at every evaluation, so a typo
        # silently kept the layer's manual conductivity.
        if self.material.lower() not in ("generic", "custom") and (
            self.material not in cb.list_materials()
        ):
            raise ValueError(
                f"WallLayer: unknown material {self.material!r}; expected one of "
                f"{cb.list_materials()} or 'custom'."
            )

    def update_conductivity(self, T: float) -> None:
        """Update conductivity based on temperature if a database material is selected."""

        # 'generic' is legacy, 'custom' is new name for manual input
        if self.material.lower() not in ("generic", "custom"):
            # By default the binding follows the material's k(T) envelope (#882),
            # bounded and C1, so a known material does not raise at any finite
            # T > 0. extrapolate continues the end slope and raises if that
            # reaches k <= 0.
            self.conductivity = cb.get_material_conductivity(
                self.material, T, extrapolate=self.extrapolate
            )
            self.in_range = bool(cb.material_in_range(self.material, T))

    @property
    def r_val(self) -> float:
        """Thermal resistance per unit area [m^2*K/W]."""
        # Defensive k-clamping to prevent division by zero/near-zero
        k_safe = max(1e-3, self.conductivity)
        return self.thickness / k_safe


class _WallCoupling:
    """A WallCouplingResult read through per-side areas (#887) and k(T)
    layers (#932): the kernel's fields as they are, plus ``kT`` -- the part of
    dQ/d(T_aw_a, T_aw_b, h_a, h_b) that moves through the layers' k(T) -- and
    dQ/dh_b in side B's own h (the kernel saw h_b A_b / A_a)."""

    _FIELDS = (
        "Q",
        "T_hot",
        "dQ_dC_a",
        "dQ_dC_b",
        "UA",
        "effectiveness",
        "dQ_dR_wall",
    )

    def __init__(
        self,
        res: Any,
        dQ_dh_b_factor: float = 1.0,
        kT: tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0),
    ) -> None:
        for name in self._FIELDS:
            setattr(self, name, getattr(res, name))
        self.dQ_dT_aw_a = res.dQ_dT_aw_a + kT[0]
        self.dQ_dT_aw_b = res.dQ_dT_aw_b + kT[1]
        self.dQ_dh_a = res.dQ_dh_a + kT[2]
        self.dQ_dh_b = (res.dQ_dh_b + kT[3]) * dQ_dh_b_factor


@dataclass
class GasRadiation:
    """Gas-to-wall radiation on one ThermalWall face (#793): Lefebvre & Ballal
    (2010) Eq. 8.7, ``R1 = 0.5 sigma (1 + eps_w) eps_g Tg^1.5 (Tg^2.5 -
    Tw^2.5)``, with the luminous emissivity of Eq. 8.11, ``eps_g = 1 -
    exp[-290 P L (q l_b)^0.5 Tg^-1.5]`` (P in kPa).

    Gray and lumped: one gas temperature, one beam length, one emissivity per
    surface. Eq. 8.11 is Reeves' kerosine-product fit times a luminosity
    factor; it is not a composition-based emissivity.

    Attributes
    ----------
    beam_length : float
        Mean beam length l_b [m] of the gas volume the face sees,
        ``cb.mean_beam_length(V, A)`` (3.4 V/A, Eq. 8.10). The user's geometry.
    emissivity_wall : float
        The face's emissivity eps_w (Lefebvre: about 0.7 Nimonic, 0.8
        stainless steel at liner temperatures).
    luminosity : float
        Luminosity factor L (1: non-luminous, Eq. 8.9). For a liquid fuel,
        ``cb.luminosity_factor_hydrogen(H)[0]`` = 336/H^2 (Eq. 8.15).
    T_gas : float | None
        Mean gas temperature [K]; None (default): the side's stream Tt.
    fuel_air_ratio : float | None
        q [-]; None (default): the side's composition, by its elements
        (``cb.elemental_fuel_air_ratio``).
    pressure : float | None
        Gas pressure [Pa]; None (default): the side's static pressure.
    emissivity_model : str
        The gas emissivity eps_g (#1071). ``"lefebvre"`` (default): Eq. 8.11
        above. ``"bordbar2021"``: Bordbar et al. (2021) pressure-dependent
        WSGG, ``cb.gas_emissivity_bordbar2021_mixture``, from the side's
        H2O and CO2 mole fractions, non-luminous (no soot), 300-3000 K,
        x_H2O/x_CO2 0.1-4, 0.1-80 atm, each clamped at its end. It takes no
        ``fuel_air_ratio`` and no ``luminosity`` (both raise when set).
    """

    beam_length: float
    emissivity_wall: float = 0.7
    luminosity: float = 1.0
    T_gas: float | None = None
    fuel_air_ratio: float | None = None
    pressure: float | None = None
    emissivity_model: str = "lefebvre"

    EMISSIVITY_MODELS = ("lefebvre", "bordbar2021")

    def __post_init__(self) -> None:
        if self.emissivity_model not in self.EMISSIVITY_MODELS:
            raise ValueError(
                f"GasRadiation: emissivity_model must be one of {self.EMISSIVITY_MODELS}, "
                f"not {self.emissivity_model!r}"
            )
        if self.emissivity_model == "bordbar2021":
            if self.fuel_air_ratio is not None:
                raise ValueError(
                    "GasRadiation: emissivity_model 'bordbar2021' takes the H2O and CO2 of the "
                    "side's composition; fuel_air_ratio does not apply"
                )
            if self.luminosity != 1.0:
                raise ValueError(
                    "GasRadiation: emissivity_model 'bordbar2021' is non-luminous (no soot); "
                    "luminosity does not apply"
                )
        if not (math.isfinite(self.beam_length) and self.beam_length > 0.0):
            raise ValueError(f"GasRadiation: beam_length must be > 0, not {self.beam_length}")
        if not 0.0 < self.emissivity_wall <= 1.0:
            raise ValueError("GasRadiation: emissivity_wall must be in (0, 1]")
        if not (math.isfinite(self.luminosity) and self.luminosity >= 0.0):
            raise ValueError("GasRadiation: luminosity must be >= 0")
        if self.T_gas is not None and not self.T_gas > 0.0:
            raise ValueError("GasRadiation: T_gas must be > 0")
        if self.fuel_air_ratio is not None and not self.fuel_air_ratio >= 0.0:
            raise ValueError("GasRadiation: fuel_air_ratio must be >= 0")
        if self.pressure is not None and not self.pressure >= 0.0:
            raise ValueError("GasRadiation: pressure must be >= 0")


@dataclass
class CasingRadiation:
    """Wall-to-casing radiation on one ThermalWall face (#793): Lefebvre &
    Ballal (2010) Eq. 8.17, a gray two-surface enclosure with view factor 1,
    ``R2 = sigma eps_w eps_c / (eps_c + eps_w (1 - eps_c) A_w/A_c) (Tw^4 -
    Tc^4)`` per unit wall area.

    The heat goes into the face's stream: in steady state a casing that loses
    nothing outward hands what it absorbs to the annulus air.

    Attributes
    ----------
    area_ratio : float
        A_w / A_c, wall over casing area (tubular: the diameter ratio;
        tubo-annular about 0.8). The user's geometry.
    emissivity_wall : float
        The face's emissivity eps_w (0.7 Nimonic, 0.8 stainless steel).
    emissivity_casing : float
        eps_c (Lefebvre's example: 0.4, aluminium casing).
    T_casing : float | None
        Casing temperature [K]; None (default): the side's stream Tt, as the
        book takes T3.
    """

    area_ratio: float
    emissivity_wall: float = 0.7
    emissivity_casing: float = 0.4
    T_casing: float | None = None

    def __post_init__(self) -> None:
        if not (math.isfinite(self.area_ratio) and self.area_ratio > 0.0):
            raise ValueError(f"CasingRadiation: area_ratio must be > 0, not {self.area_ratio}")
        for name in ("emissivity_wall", "emissivity_casing"):
            if not 0.0 < getattr(self, name) <= 1.0:
                raise ValueError(f"CasingRadiation: {name} must be in (0, 1]")
        if self.T_casing is not None and not self.T_casing >= 0.0:
            raise ValueError("CasingRadiation: T_casing must be >= 0")


@dataclass
class _RadiationFace:
    """One face's radiation, ready for cb's wall_radiation_sides: kind, the
    source temperature, the coefficient k (in side A's area units), and
    their partials in the side's stream (Tt, P, Y)."""

    kind: int = 0
    T_src: float = 0.0
    k: float = 0.0
    dTsrc_dTt: float = 0.0
    dk_dTt: float = 0.0
    dk_dP: float = 0.0
    dk_dY: Any = None
    eps_g: float = math.nan


def _radiation_face(
    spec: "GasRadiation | CasingRadiation | None",
    gas: "tuple[float, float, Any] | None",
    scale: float,
) -> _RadiationFace:
    """A face's radiation spec at its stream state ``gas`` = (Tt, P, Y);
    ``scale`` = this face's area over side A's (per-side areas, #887)."""
    if spec is None:
        return _RadiationFace()
    needs_gas = (
        (
            spec.T_gas is None
            or spec.fuel_air_ratio is None
            or spec.pressure is None
            or spec.emissivity_model == "bordbar2021"
        )
        if isinstance(spec, GasRadiation)
        else spec.T_casing is None
    )
    if needs_gas and gas is None:
        raise ValueError("ThermalWall radiation: the face's stream state is needed")
    if isinstance(spec, GasRadiation):
        Tt, P, Y = gas if gas is not None else (math.nan, math.nan, None)
        T_g = float(spec.T_gas) if spec.T_gas is not None else float(Tt)
        P_g = float(spec.pressure) if spec.pressure is not None else float(P)
        c = 0.5 * cb.STEFAN_BOLTZMANN * (1.0 + spec.emissivity_wall) * scale
        from_stream = spec.T_gas is None
        if spec.emissivity_model == "bordbar2021":
            w = cb.gas_emissivity_bordbar2021_mixture(T_g, P_g, list(Y), spec.beam_length)
            eps, deps_dT, deps_dP = w.eps, w.deps_dT, w.deps_dP
            dk_dY = [c * float(v) for v in w.deps_dY]
        else:
            if spec.fuel_air_ratio is not None:
                q, dq_dY = float(spec.fuel_air_ratio), None
            else:
                q, dq_dY = cb.elemental_fuel_air_ratio(list(Y))
            e = cb.gas_emissivity_luminous(T_g, P_g, q, spec.beam_length, spec.luminosity)
            eps, deps_dT, deps_dP = e.eps, e.deps_dT, e.deps_dP
            dk_dY = None if dq_dY is None else [c * e.deps_dq * float(v) for v in dq_dY]
        return _RadiationFace(
            kind=1,
            T_src=T_g,
            k=c * eps,
            dTsrc_dTt=1.0 if from_stream else 0.0,
            dk_dTt=c * deps_dT if from_stream else 0.0,
            dk_dP=c * deps_dP if spec.pressure is None else 0.0,
            dk_dY=dk_dY,
            eps_g=float(eps),
        )
    T_c = float(spec.T_casing) if spec.T_casing is not None else float(gas[0])
    F = cb.liner_casing_exchange_factor(
        spec.emissivity_wall, spec.emissivity_casing, spec.area_ratio
    )
    return _RadiationFace(
        kind=2,
        T_src=T_c,
        k=F * scale,
        dTsrc_dTt=1.0 if spec.T_casing is None else 0.0,
    )


@dataclass(frozen=True)
class HeatLoss:
    """A ThermalWall side with no stream: a fixed external temperature behind
    a heat transfer coefficient (#970), ``ThermalWall(element_b=None,
    side_b=HeatLoss(T_ext=..., h=...))``.

    The wall's one stream exchanges with ``T_ext`` through the series
    resistance ``1/(h_a A_a) + sum t/(k A_m) + R_fouling/A_a + 1/(h A_b)``,
    the epsilon-NTU integral against a reservoir: bounded by C_a (T_aw_a -
    T_ext), none into a stream that does not flow. The heat leaves the network
    (``{wall}.Q_loss`` [W], > 0 out) and the ledger books it so.

    Attributes
    ----------
    T_ext : float
        External temperature [K], > 0.
    h : float
        External heat transfer coefficient [W/(m^2*K)], > 0 (radiation
        linearised into it, if any: an input the user owns).
    """

    T_ext: float
    h: float

    def __post_init__(self) -> None:
        if not (math.isfinite(float(self.T_ext)) and float(self.T_ext) > 0.0):
            raise ValueError(
                f"HeatLoss: T_ext must be a positive temperature [K], not {self.T_ext}"
            )
        if not (math.isfinite(float(self.h)) and float(self.h) > 0.0):
            raise ValueError(
                f"HeatLoss: h must be > 0 [W/(m^2*K)], not {self.h}; an insulated side "
                "with no loss is side_b='adiabatic'"
            )


class _AdiabaticWallResult:
    """The coupling of a wall with an adiabatic side (#970): no heat, the
    stream side's face at its own balance (``cb wall_adiabatic_face``)."""

    def __init__(self, T_s: float, radiation: dict[str, float]) -> None:
        self.Q = 0.0
        self.T_hot = T_s
        self.UA = 0.0
        self.effectiveness = 0.0
        self.radiation = radiation


class _RadiativeWallCoupling:
    """A wall's coupling with radiation on a face (#793): the kernel's result
    on the effective sides, with dQ/d(T_aw, h) back in the convective inputs
    and the radiation's own partials in each side's stream: ``rad_dQ_dTt_a``,
    ``rad_dQ_dP_a``, ``rad_dQ_dY_a`` (None when the composition does not
    enter), and side B's. ``radiation`` holds the reported values."""

    def __init__(self, **kw: Any) -> None:
        self.__dict__.update(kw)


@dataclass
class ThermalWall:
    """Thermal coupling between two elements through a shared multi-layer wall.

    What it represents. One shared wall PATCH with two sides: element A's
    fluid wets one face, element B's the other, and the heat crosses the metal
    straight through, normal to the wall (1-D, through the layers).

    - Its area is the shared part of each element's perimeter, not the whole
      wetted perimeter: ``contact_area``, or per-side ``area_a`` / ``area_b``
      when the two faces differ (#887).
    - An element may carry several walls, one per face it shares (the four
      faces of a rectangular duct, or part of a round one's perimeter). A
      face with no wall is adiabatic.
    - Each side's h is its element's own (the channel correlation's
      perimeter-mean h, its ``Nu_multiplier`` included), applied to this
      patch: there is no face-by-face h.
    - Along the wall the exchange is the two streams' epsilon-NTU integral
      (below); parallel or counter flow follows from ``orientation`` and the
      flow signs.
    - No heat flows ALONG the metal: walls of elements in series (a duct and
      the next one downstream) exchange heat only through the fluid. Their
      axial path would conduct on the metal's cross-section, thickness x width,
      over the distance between segments -- usually negligible, but not for
      thick metal, short segments or a steep gradient (#898).

    Attributes
    ----------
    id : str
        Unique identifier for this wall connection.
    element_a : str
        Element ID for side A.
    element_b : str | None
        Element ID for side B; None for a wall with no stream on side B,
        which ``side_b`` then describes (#970).
    side_b : "adiabatic" | HeatLoss | None
        Side B without a stream (#970), with ``element_b=None``; None (the
        default) for a stream on each side. ``"adiabatic"``: no heat crosses
        the wall, which sits at side A's own balance -- its T_aw, or with
        ``radiation_a=GasRadiation`` the face's radiation balance
        (``wall_adiabatic_face``) -- and changes nothing in the network.
        ``HeatLoss(T_ext, h)``: side A's stream exchanges with a fixed
        ``T_ext`` behind ``h`` through the layers (epsilon-NTU against a
        reservoir), and that heat leaves the network (``{wall}.Q_loss``).
        The free side is B only: the layers run from the stream to it. Not
        with ``radiation_b`` or a film on side B; an adiabatic wall also not
        with ``radiation_a=CasingRadiation``, whose heat would leave to a
        casing no stream carries (use HeatLoss).
    layers : list[WallLayer]
        Stack of wall layers ordered from side A to side B.
    contact_area : float | None
        Override effective area [m^2]. If None, uses min(A_conv_a, A_conv_b).
    area_a, area_b : float | None
        Per-side areas [m^2] (#887): the wall's face on side A and on side B,
        for a wall whose two faces differ -- a hole in a thick wall, a finned
        face, a curved wall. Both or neither; not with ``contact_area``. Each
        side's convective resistance is formed on its own face, the layers'
        conduction on ``wall_area``'s area, and the fouling on side A's:
        ``1/UA = 1/(h_a A_a) + sum t/(k A_m) + R_fouling/A_a + 1/(h_b A_b)``.
        Neither (the default): one area for all, as before.
    wall_area : {"log_mean", "slab_a", "slab_b"}
        The conduction area of each layer with per-side areas. "log_mean"
        (the default): the area runs linearly through the stack from A_a to
        A_b (as a cylinder's does with its radius) and each layer conducts on
        the log mean of its two faces, which is exact for a thick cylinder,
        ``Q = 2 pi k L dT / ln(r_o / r_i)``, layer by layer. "slab_a" /
        "slab_b": a plain slab on side A's or side B's area.
    R_fouling : float
        Additional fouling resistance [m^2*K/W] (applied to cold side).
    orientation : {"aligned", "opposed"}
        How the two elements lie along the wall, as declared (from_node ->
        to_node): "aligned" puts both from_nodes at the same end. Geometry, not
        an arrangement: with the flow signs it decides whether the streams run
        parallel or counter along the wall, so the arrangement follows a flow
        that reverses (#627). Irrelevant for a node side, which is well mixed.
    film : str | None
        Id of a film provider (a FilmCoolingElement, or an
        EffusionPlateElement with panel_length) whose film covers this
        wall (#386). The ``film_side`` side's adiabatic wall temperature
        becomes ``T_aw - eta_bar (T_aw - T_coolant)``; its h, stream and heat
        booking stay its own. That side must be the gas duct downstream of
        the holes (an element with a flow area) or a MomentumChamberNode.
    film_side : {"a", "b"}
        The side the film covers.
    film_start, film_end : float
        The wall's span along the gas path [m], measured downstream of the
        provider's LAST row; eta_bar is the film averaged over it.
        ``film_end`` defaults to ``film_start`` (the point value).
    radiation_a, radiation_b : GasRadiation | CasingRadiation | None
        Radiation on side A's / side B's face (#793), default None (none):
        ``GasRadiation`` from that side's gas to the wall (Lefebvre Eq. 8.7)
        or ``CasingRadiation`` from the wall to a casing (Eq. 8.17). Each
        adds to its face's convection in the local wall balance, solved for
        the two surface temperatures; the exchange then runs on each side's
        combined conductance and driving temperature. All of it is booked to
        the two streams: side A's gives it up, side B's takes it.

    The heat is the two streams' exchange, the epsilon-NTU integral along the
    wall at each stream's inlet temperature (``wall_exchange_and_jacobian``):
    bounded by C_min (T_hot,in - T_cold,in), and none into a stream that does
    not flow. The linear U A dT it replaced is its small-NTU limit.
    """

    id: str
    element_a: str
    element_b: str | None
    layers: list[WallLayer] = field(default_factory=list)
    contact_area: float | None = None
    area_a: float | None = None
    area_b: float | None = None
    wall_area: Literal["log_mean", "slab_a", "slab_b"] = "log_mean"
    R_fouling: float = 0.0
    orientation: Literal["aligned", "opposed"] = "aligned"
    film: str | None = None
    film_side: Literal["a", "b"] = "a"
    film_start: float = 0.0
    film_end: float | None = None
    radiation_a: GasRadiation | CasingRadiation | None = None
    radiation_b: GasRadiation | CasingRadiation | None = None
    side_b: Literal["adiabatic"] | HeatLoss | None = None

    # Internal cache for temperature-dependent property iteration
    _last_profile: list[float] | None = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        if (self.area_a is None) != (self.area_b is None):
            raise ValueError(
                f"ThermalWall '{self.id}': give both area_a and area_b (the wall's two "
                "faces), or neither"
            )
        if self.area_a is not None:
            if self.contact_area is not None:
                raise ValueError(
                    f"ThermalWall '{self.id}': contact_area and per-side areas (area_a, "
                    "area_b) are two descriptions of the same wall; give one"
                )
            for name, a in (("area_a", self.area_a), ("area_b", self.area_b)):
                if not (math.isfinite(float(a)) and float(a) > 0.0):
                    raise ValueError(f"ThermalWall '{self.id}': {name} must be > 0, not {a}")
        if self.wall_area not in ("log_mean", "slab_a", "slab_b"):
            raise ValueError(
                f"ThermalWall '{self.id}': wall_area must be 'log_mean', 'slab_a' or "
                f"'slab_b', not {self.wall_area!r}"
            )
        if self.orientation not in ("aligned", "opposed"):
            raise ValueError(
                f"ThermalWall '{self.id}': orientation must be 'aligned' or 'opposed', "
                f"not {self.orientation!r}"
            )
        if self.film is not None:
            if self.film_side not in ("a", "b"):
                raise ValueError(f"ThermalWall '{self.id}': film_side must be 'a' or 'b'")
            if self.film_end is None:
                self.film_end = self.film_start
            if not 0.0 <= self.film_start <= self.film_end:
                raise ValueError(
                    f"ThermalWall '{self.id}': the film span needs "
                    "0 <= film_start <= film_end [m, downstream of the last row]"
                )
        for name in ("radiation_a", "radiation_b"):
            spec = getattr(self, name)
            if spec is not None and not isinstance(spec, (GasRadiation, CasingRadiation)):
                raise ValueError(
                    f"ThermalWall '{self.id}': {name} must be a GasRadiation, a "
                    f"CasingRadiation or None, not {type(spec).__name__}"
                )
        self._check_side_b()

    def _check_side_b(self) -> None:
        """A side B with no stream (#970): what it is and what it excludes."""
        if self.element_b is not None:
            if self.side_b is not None:
                raise ValueError(
                    f"ThermalWall '{self.id}': side_b={self.side_b!r} describes a side B "
                    f"with no stream, but element_b is {self.element_b!r}; give "
                    "element_b=None, or drop side_b"
                )
            return
        if not (self.side_b == "adiabatic" or isinstance(self.side_b, HeatLoss)):
            raise ValueError(
                f"ThermalWall '{self.id}': element_b=None needs side_b='adiabatic' or "
                f"side_b=HeatLoss(T_ext=..., h=...), not {self.side_b!r}"
            )
        if self.radiation_b is not None:
            raise ValueError(
                f"ThermalWall '{self.id}': side B has no stream, so no radiation_b; "
                "linearise any radiation from that face into HeatLoss's h"
            )
        if self.film is not None and self.film_side != "a":
            raise ValueError(
                f"ThermalWall '{self.id}': side B has no stream for a film to cover; "
                "film_side must be 'a'"
            )
        if self.side_b == "adiabatic" and isinstance(self.radiation_a, CasingRadiation):
            raise ValueError(
                f"ThermalWall '{self.id}': an adiabatic wall with CasingRadiation on "
                "side A would lose heat to a casing that no stream carries; use "
                "side_b=HeatLoss(...) for a lossy wall, or a casing stream as side B"
            )

    @property
    def one_sided(self) -> bool:
        """True when side B has no stream (``side_b``, #970)."""
        return self.element_b is None

    @property
    def adiabatic(self) -> bool:
        """True for an adiabatic side B: the wall passes no heat (#970)."""
        return self.element_b is None and self.side_b == "adiabatic"

    @property
    def heat_loss(self) -> HeatLoss | None:
        """Side B's HeatLoss, or None (#970)."""
        return self.side_b if isinstance(self.side_b, HeatLoss) else None

    @property
    def endpoints(self) -> tuple[str, ...]:
        """The wall's endpoints that carry a stream: (A, B), or (A,) with a
        side B that has none (#970)."""
        if self.element_b is None:
            return (self.element_a,)
        return (self.element_a, self.element_b)

    @property
    def heated_endpoints(self) -> tuple[str, ...]:
        """The endpoints whose stream takes or gives this wall's heat: none
        for an adiabatic wall (#970)."""
        return () if self.adiabatic else self.endpoints

    def adiabatic_coupling(
        self,
        h_a: float,
        T_aw_a: float,
        A_conv_a: float,
        gas_a: "tuple[float, float, Any] | None" = None,
    ) -> _AdiabaticWallResult:
        """An adiabatic wall's state (#970): no heat, and the wall at side A's
        face balance, ``h (T_aw - T_s) + h_r (T_src - T_s) = 0`` (C++
        ``wall_adiabatic_face``); without radiation T_s = T_aw_a. No heat
        crosses the layers, so the whole stack sits at T_s."""
        face = _radiation_face(self.radiation_a, gas_a, 1.0)
        h_p = max(float(h_a), cb.WALL_HTC_KNEE)
        r = _solver_tools.wall_adiabatic_face(float(T_aw_a), h_p, face.kind, face.T_src, face.k)
        T_s = float(r.T_s)
        for layer in self.layers:
            if layer.material.lower() not in ("generic", "custom"):
                layer.update_conductivity(T_s)
                layer.in_range = bool(cb.material_in_range(layer.material, T_s))
        self._last_profile = [T_s] * (len(self.layers) + 1)
        radiation: dict[str, float] = {}
        if face.kind:
            if self.area_a is not None:
                A = float(self.area_a)
            elif self.contact_area is not None:
                A = float(self.contact_area)
            else:
                A = float(A_conv_a)
            radiation["q_rad_a"] = float(r.q_rad)
            radiation["Q_rad_a"] = float(r.q_rad) * A
            radiation["T_surface_a"] = T_s
            if face.kind == 1:
                radiation["eps_g_a"] = face.eps_g
        return _AdiabaticWallResult(T_s, radiation)

    def heat_loss_coupling(
        self,
        h_a: float,
        T_aw_a: float,
        A_conv_a: float,
        C_a: float = math.inf,
        gas_a: "tuple[float, float, Any] | None" = None,
    ) -> cb.WallCouplingResult:
        """A HeatLoss wall's exchange (#970): ``compute_coupling`` with side B
        the external reservoir (h, T_ext, C = inf), so Q > 0 leaves side A's
        stream for T_ext, bounded by C_a (T_aw_a - T_ext)."""
        hl = self.heat_loss
        if hl is None:
            raise ValueError(f"ThermalWall '{self.id}': no HeatLoss side")
        return self.compute_coupling(
            h_a,
            T_aw_a,
            A_conv_a,
            float(hl.h),
            float(hl.T_ext),
            math.inf,
            C_a,
            math.inf,
            False,
            gas_a=gas_a,
            gas_b=None,
        )

    @property
    def film_side_id(self) -> str | None:
        """The wall endpoint the film covers, or None without a film."""
        if self.film is None:
            return None
        return self.element_a if self.film_side == "a" else self.element_b

    def _layer_area_ratios(self) -> list[float]:
        """A_a / A_m per layer (#887): side A's area over the area the layer
        conducts on. With the area linear in depth from A_a to A_b, a layer
        between faces A_i and A_o conducts on their log mean (exact for a
        cylinder); a zero-thickness stack or equal faces give the slab."""
        A_a, A_b = float(self.area_a), float(self.area_b)
        if self.wall_area == "slab_a" or A_a == A_b:
            return [1.0] * len(self.layers)
        if self.wall_area == "slab_b":
            return [A_a / A_b] * len(self.layers)
        total = sum(float(L.thickness) for L in self.layers)
        out: list[float] = []
        depth = 0.0
        for L in self.layers:
            if total <= 0.0:
                out.append(1.0)
                continue
            A_i = A_a + (A_b - A_a) * depth / total
            depth += float(L.thickness)
            A_o = A_a + (A_b - A_a) * depth / total
            A_m = A_i if A_o == A_i else (A_o - A_i) / math.log(A_o / A_i)
            out.append(A_a / A_m)
        return out

    def is_counterflow(self, sign_a: float, sign_b: float) -> bool:
        """Whether streams of these flow signs (along each element's declared
        from -> to) run opposite ways along the wall."""
        return (sign_a * sign_b < 0.0) != (self.orientation == "opposed")

    def compute_coupling(
        self,
        h_a: float,
        T_aw_a: float,
        A_conv_a: float,
        h_b: float,
        T_aw_b: float,
        A_conv_b: float,
        C_a: float = math.inf,
        C_b: float = math.inf,
        counterflow: bool = False,
        gas_a: "tuple[float, float, Any] | None" = None,
        gas_b: "tuple[float, float, Any] | None" = None,
    ) -> cb.WallCouplingResult:
        """Compute heat transfer rate and wall temperature using analytical Jacobians.

        Parameters
        ----------
        h_a, h_b : float
            Convective HTCs on sides A and B [W/(m^2*K)].
        T_aw_a, T_aw_b : float
            Adiabatic wall temperatures on sides A and B [K], at each stream's
            inlet.
        A_conv_a, A_conv_b : float
            Convective areas on both sides [m^2].
        C_a, C_b : float
            Capacity rates |m| cp [W/K]; inf (the default) for a node side, a
            reservoir. Two reservoirs give the linear U A dT.
        counterflow : bool
            The streams run opposite ways along the wall (``is_counterflow``).
        gas_a, gas_b : tuple | None
            Each side's stream state (Tt [K], P [Pa], Y), read by that side's
            radiation (``radiation_a`` / ``radiation_b``) only.

        Returns
        -------
        WallCouplingResult
            Object containing Q [W], T_hot [K], and analytical Jacobians.
        """
        # Per-side areas (#887): everything per unit of side A's area. Side B's
        # convective resistance 1/(h_b A_b) is 1/(h_b' A_a) with h_b' = h_b
        # A_b/A_a, and a layer's t/(k A_m) is (t/k)(A_a/A_m) / A_a, so the one-
        # area kernel gives the exact Q, wall temperatures and Jacobian
        # (dQ/dh_b = dQ/dh_b' A_b/A_a). Without them, the one area as before.
        per_side = self.area_a is not None
        if per_side:
            A_eff = float(self.area_a)
            r_b = float(self.area_b) / A_eff
            h_b = float(h_b) * r_b
            layer_scale = self._layer_area_ratios()
        else:
            A_eff = self.contact_area if self.contact_area is not None else min(A_conv_a, A_conv_b)
            layer_scale = [1.0] * len(self.layers)

        # The profile only places the layer temperatures; it refuses h <= 0,
        # which a correlation past its range can return at a solver iterate,
        # so it sees the coupling's own floor (WALL_HTC_KNEE).
        h_pa = max(float(h_a), cb.WALL_HTC_KNEE)
        h_pb = max(float(h_b), cb.WALL_HTC_KNEE)

        if self.radiation_a is not None or self.radiation_b is not None:
            return self._radiative_coupling(
                h_a,
                T_aw_a,
                h_b,
                T_aw_b,
                h_pa,
                h_pb,
                C_a,
                C_b,
                counterflow,
                A_eff,
                layer_scale,
                r_b if per_side else 1.0,
                gas_a,
                gas_b,
            )

        def profile(t_over_k: list[float]) -> list[float]:
            prof, _q = cb.wall_temperature_profile(
                T_aw_a, T_aw_b, h_pa, h_pb, t_over_k, self.R_fouling
            )
            return [float(tp) for tp in prof]

        def resistances() -> list[float]:
            return [L.r_val * f for L, f in zip(self.layers, layer_scale, strict=True)]

        t_over_k_layers = resistances()
        prof = profile(t_over_k_layers)
        kT: tuple[float, float, float, float] | None = None
        if any(L.material.lower() not in ("generic", "custom") for L in self.layers):
            # k(T) at THIS call's wall temperatures, iterated to its fixed
            # point, so Q depends on the inputs alone. It used to be lagged
            # from the previous call -- any trial point the solver probed --
            # which made the residual depend on evaluation history (#481).
            for _ in range(30):
                for i, layer in enumerate(self.layers):
                    layer.update_conductivity(0.5 * (prof[i] + prof[i + 1]))
                new = resistances()
                change = max(
                    abs(a - b) / max(abs(b), 1e-300)
                    for a, b in zip(new, t_over_k_layers, strict=True)
                )
                t_over_k_layers = new
                prof = profile(t_over_k_layers)
                if change < 1e-12:
                    break
            # Whether each layer's k came from its material's sourced range at
            # its mean temperature in the profile just placed (#882).
            for i, layer in enumerate(self.layers):
                if layer.material.lower() not in ("generic", "custom"):
                    layer.in_range = bool(
                        cb.material_in_range(layer.material, 0.5 * (prof[i] + prof[i + 1]))
                    )
            # The Jacobian through k(T) (#932): Q moves with each layer's
            # t/k, and t/k with the layer's mean temperature, itself a
            # function of the inputs and of every t/k. The C++ side
            # differentiates the fixed point; here only each layer's
            # d(t/k)/dT = -t (A_a/A_m) (dk/dT) / k^2 from the material's
            # analytic slope (envelope included), 0 for a constant k.
            dr_dTm: list[float] = []
            for i, (layer, f) in enumerate(zip(self.layers, layer_scale, strict=True)):
                if layer.material.lower() in ("generic", "custom"):
                    dr_dTm.append(0.0)
                    continue
                k, dk = cb.get_material_conductivity_and_slope(
                    layer.material, 0.5 * (prof[i] + prof[i + 1]), extrapolate=layer.extrapolate
                )
                dr_dTm.append(-float(layer.thickness) * f * dk / (k * k))
            dR = _solver_tools.wall_resistance_fixed_point_derivatives(
                T_aw_a, T_aw_b, h_pa, h_pb, t_over_k_layers, self.R_fouling, dr_dTm
            )
            # The profile saw h at the WALL_HTC_KNEE floor: flat below it.
            kT = (
                float(dR[0]),
                float(dR[1]),
                float(dR[2]) if float(h_a) >= cb.WALL_HTC_KNEE else 0.0,
                float(dR[3]) if float(h_b) >= cb.WALL_HTC_KNEE else 0.0,
            )

        res = _solver_tools.wall_exchange_and_jacobian(
            h_a,
            T_aw_a,
            C_a,
            h_b,
            T_aw_b,
            C_b,
            t_over_k_layers,
            A_eff,
            self.R_fouling,
            counterflow,
        )
        self._last_profile = prof
        if kT is not None:
            g = res.dQ_dR_wall
            kT = (g * kT[0], g * kT[1], g * kT[2], g * kT[3])
        if per_side or kT is not None:
            return _WallCoupling(
                res, dQ_dh_b_factor=r_b if per_side else 1.0, kT=kT or (0.0, 0.0, 0.0, 0.0)
            )
        return res

    def _radiative_coupling(
        self,
        h_a: float,
        T_aw_a: float,
        h_b: float,
        T_aw_b: float,
        h_pa: float,
        h_pb: float,
        C_a: float,
        C_b: float,
        counterflow: bool,
        A_eff: float,
        layer_scale: list[float],
        r_b: float,
        gas_a: "tuple[float, float, Any] | None",
        gas_b: "tuple[float, float, Any] | None",
    ) -> _RadiativeWallCoupling:
        """compute_coupling with radiation on a face (#793), per unit of side
        A's area (h_b, and side B's radiation coefficient, scaled by r_b).

        The local balance at the streams' inlet states gives the two surface
        temperatures and, per side, the effective conductance H = h + h_r and
        driving temperature T* (cb wall_radiation_sides); the exchange is
        wall_exchange_and_jacobian on those, so between two reservoirs Q / A
        is the face's convection plus R1 (or R2) exactly. With k(T) layers
        the conductive resistance R_wall and the radiation share one fixed
        point; they couple only through the scalar R_wall, so the total
        derivative is the radiation's plus a scalar implicit term."""
        face_a = _radiation_face(self.radiation_a, gas_a, 1.0)
        face_b = _radiation_face(self.radiation_b, gas_b, r_b)

        def resistances() -> list[float]:
            return [L.r_val * f for L, f in zip(self.layers, layer_scale, strict=True)]

        def solve(R_wall: float) -> Any:
            return _solver_tools.wall_radiation_sides(
                T_aw_a,
                h_pa,
                T_aw_b,
                h_pb,
                R_wall,
                face_a.kind,
                face_a.T_src,
                face_a.k,
                face_b.kind,
                face_b.T_src,
                face_b.k,
            )

        def profile(rad: Any, t_over_k: list[float]) -> list[float]:
            prof, _q = cb.wall_temperature_profile(
                rad.T_a, rad.T_b, rad.H_a, rad.H_b, t_over_k, self.R_fouling
            )
            return [float(tp) for tp in prof]

        t_over_k_layers = resistances()
        rad = solve(sum(t_over_k_layers) + self.R_fouling)
        prof = profile(rad, t_over_k_layers)
        kT_dR: list[float] | None = None
        if any(L.material.lower() not in ("generic", "custom") for L in self.layers):
            # k(T) at this call's wall temperatures, iterated jointly with the
            # radiation's surface temperatures to their common fixed point.
            for _ in range(60):
                for i, layer in enumerate(self.layers):
                    layer.update_conductivity(0.5 * (prof[i] + prof[i + 1]))
                new = resistances()
                change = max(
                    abs(a - b) / max(abs(b), 1e-300)
                    for a, b in zip(new, t_over_k_layers, strict=True)
                )
                t_over_k_layers = new
                rad = solve(sum(t_over_k_layers) + self.R_fouling)
                prof = profile(rad, t_over_k_layers)
                if change < 1e-12:
                    break
            for i, layer in enumerate(self.layers):
                if layer.material.lower() not in ("generic", "custom"):
                    layer.in_range = bool(
                        cb.material_in_range(layer.material, 0.5 * (prof[i] + prof[i + 1]))
                    )
            dr_dTm: list[float] = []
            for i, (layer, f) in enumerate(zip(self.layers, layer_scale, strict=True)):
                if layer.material.lower() in ("generic", "custom"):
                    dr_dTm.append(0.0)
                    continue
                k, dk = cb.get_material_conductivity_and_slope(
                    layer.material, 0.5 * (prof[i] + prof[i + 1]), extrapolate=layer.extrapolate
                )
                dr_dTm.append(-float(layer.thickness) * f * dk / (k * k))
            # dR_wall/d(T_hot, T_cold, h_hot, h_cold) at fixed effective sides.
            kT_dR = list(
                _solver_tools.wall_resistance_fixed_point_derivatives(
                    rad.T_a, rad.T_b, rad.H_a, rad.H_b, t_over_k_layers, self.R_fouling, dr_dTm
                )
            )

        res = _solver_tools.wall_exchange_and_jacobian(
            rad.H_a,
            rad.T_a,
            C_a,
            rad.H_b,
            rad.T_b,
            C_b,
            t_over_k_layers,
            A_eff,
            self.R_fouling,
            counterflow,
        )
        self._last_profile = prof

        # Chain: Q(s, R_wall), s = (H_a, T*_a, H_b, T*_b) = S(p, R_wall),
        # p = [T_aw_a, h_a, T_aw_b, h_b, R_wall, T_src_a, k_a, T_src_b, k_b].
        S = [list(row) for row in rad.d]
        dQ_ds = [res.dQ_dh_a, res.dQ_dT_aw_a, res.dQ_dh_b, res.dQ_dT_aw_b]
        n_p = len(S[0])
        dR_dp = [0.0] * n_p
        if kT_dR is not None:
            # R_wall = K(s) with s = S(p, R_wall): dR/dp = K.S_p / (1 - K.S_R).
            K = [kT_dR[2], kT_dR[0], kT_dR[3], kT_dR[1]]
            KS_R = sum(K[i] * S[i][4] for i in range(4))
            for j in range(n_p):
                if j != 4:
                    dR_dp[j] = sum(K[i] * S[i][j] for i in range(4)) / (1.0 - KS_R)
        dQ_dp = [0.0] * n_p
        for j in range(n_p):
            if j == 4:
                continue
            ds = [S[i][j] + S[i][4] * dR_dp[j] for i in range(4)]
            dQ_dp[j] = sum(dQ_ds[i] * ds[i] for i in range(4)) + res.dQ_dR_wall * dR_dp[j]

        def stream_partials(face: _RadiationFace, j_src: int, j_k: int) -> tuple:
            d_src, d_k = dQ_dp[j_src], dQ_dp[j_k]
            dY = None
            if face.dk_dY is not None:
                dY = [d_k * v for v in face.dk_dY]
            return d_src * face.dTsrc_dTt + d_k * face.dk_dTt, d_k * face.dk_dP, dY

        dTt_a, dP_a, dY_a = stream_partials(face_a, 5, 6)
        dTt_b, dP_b, dY_b = stream_partials(face_b, 7, 8)
        # The radiative fluxes: per unit of each face's own area, and in W.
        radiation: dict[str, float] = {}
        if face_a.kind:
            radiation["q_rad_a"] = float(rad.q_rad_a)
            radiation["Q_rad_a"] = float(rad.q_rad_a) * A_eff
            radiation["T_surface_a"] = float(rad.T_s_a)
            if face_a.kind == 1:
                radiation["eps_g_a"] = face_a.eps_g
        if face_b.kind:
            radiation["q_rad_b"] = float(rad.q_rad_b) / r_b
            radiation["Q_rad_b"] = float(rad.q_rad_b) * A_eff
            radiation["T_surface_b"] = float(rad.T_s_b)
            if face_b.kind == 1:
                radiation["eps_g_b"] = face_b.eps_g
        return _RadiativeWallCoupling(
            Q=res.Q,
            T_hot=res.T_hot,
            dQ_dC_a=res.dQ_dC_a,
            dQ_dC_b=res.dQ_dC_b,
            UA=res.UA,
            effectiveness=res.effectiveness,
            dQ_dR_wall=res.dQ_dR_wall,
            dQ_dT_aw_a=dQ_dp[0],
            # The profile saw h at the WALL_HTC_KNEE floor: flat below it.
            dQ_dh_a=dQ_dp[1] if float(h_a) >= cb.WALL_HTC_KNEE else 0.0,
            dQ_dT_aw_b=dQ_dp[2],
            dQ_dh_b=(dQ_dp[3] if float(h_b) >= cb.WALL_HTC_KNEE else 0.0) * r_b,
            rad_dQ_dTt_a=dTt_a,
            rad_dQ_dP_a=dP_a,
            rad_dQ_dY_a=dY_a,
            rad_dQ_dTt_b=dTt_b,
            rad_dQ_dP_b=dP_b,
            rad_dQ_dY_b=dY_b,
            radiation=radiation,
        )


class EnergyBoundary:
    """Energy source/sink that attaches to mixing nodes.

    Convention: Q > 0 / fraction > 0 means heat **into** the fluid (heating),
                Q < 0 / fraction < 0 means heat **out of** the fluid (cooling).

    Two additive modes:
      - Q [W]: absolute heat transfer rate, spread over the node's inflow
      - fraction [-]: of the inflow's SENSIBLE enthalpy, h(T) - h(298.15 K)
        at its own composition (cb.SENSIBLE_ENTHALPY_REF_T); on a combustor,
        of the products' at T_ad. Not of absolute enthalpy, which includes
        formation enthalpy (#481). fraction=-0.05 removes 5% of it.

    Both can be set simultaneously. The heat a fraction applied is reported
    as {node}.Q_fraction [W]. C++ handles the conversion to delta_h and all
    Jacobian corrections.

    Third mode, a set outlet temperature (#792):
      - Tt_target [K]: the node's TOTAL temperature after the heat. The heat
        flow becomes a solver unknown, ``{id}.Q`` [W], closed by one residual
        Tt_node - Tt_target = 0 (a temperature row). Exclusive: Q and
        fraction must stay 0, and a node takes at most one such boundary. On a
        CombustorNode it is the heat taken out (or put in) AFTER the reaction
        to reach that exit Tt; it does not change the fuel flow.

    In every mode the solve reports the heat a boundary applied as ``{id}.Q``
    [W] -- its Q, its solved Q, or its share of ``{node}.Q_fraction`` -- and
    their sum per node as ``{node}.Q``. ``heat`` is the Q part alone (what
    enters the mixer as Q; the fraction enters as a fraction).
    """

    def __init__(
        self,
        id: str,
        Q: float = 0.0,
        fraction: float = 0.0,
        Tt_target: float | None = None,
    ) -> None:
        self.id = id
        self.Q = Q  # [W]
        self.fraction = fraction  # [-]
        self.Tt_target = Tt_target  # [K], total; None = Q/fraction mode
        # Target mode: the solver writes its unknown's current value here.
        self.Q_solved = 0.0  # [W]
        self.check()

    @property
    def is_target(self) -> bool:
        """True when the heat flow is an unknown closed by Tt_target."""
        return self.Tt_target is not None

    @property
    def unknown_name(self) -> str:
        """The solver unknown of the target mode, ``{id}.Q`` [W]."""
        return f"{self.id}.Q"

    @property
    def heat(self) -> float:
        """The heat flow this boundary applies [W] (> 0 heats the fluid): Q,
        or in target mode the solver's current value of its unknown."""
        return float(self.Q_solved) if self.is_target else float(self.Q)

    def check(self) -> None:
        """Refuse an over-specified or meaningless target, or a fraction that
        removes all the sensible enthalpy or more (ValueError)."""
        _check_sensible_fraction(f"EnergyBoundary '{self.id}'", float(self.fraction))
        if self.Tt_target is None:
            return
        T = float(self.Tt_target)
        if not (math.isfinite(T) and T > 0.0):
            raise ValueError(
                f"EnergyBoundary '{self.id}': Tt_target must be a positive, finite "
                f"total temperature [K], got {self.Tt_target!r}."
            )
        if self.Q != 0.0 or self.fraction != 0.0:
            raise ValueError(
                f"EnergyBoundary '{self.id}': Tt_target is set together with "
                f"Q={self.Q!r} / fraction={self.fraction!r}. With a target outlet "
                "temperature the heat flow is solved for, so Q and fraction must "
                "be 0 (over-specified); use a second boundary for a fixed heat."
            )


#: Node types that can honour an EnergyBoundary(Tt_target=...): the ones
#: whose derived Tt is a mixer (or combustor) output the heat shifts.
def _energy_target_node_types() -> tuple[type, ...]:
    return (PlenumNode, MomentumChamberNode, CombustorNode)


def target_energy_boundaries(node: Any) -> list[EnergyBoundary]:
    """The node's EnergyBoundaries in Tt_target mode (at most one, validated)."""
    return [eb for eb in getattr(node, "energy_boundaries", None) or [] if eb.is_target]


def validate_energy_boundaries(node: Any) -> None:
    """Network-build checks of a node's EnergyBoundaries (#792): each one's own
    consistency, at most one Tt_target per node, and a Tt_target only on a
    node type whose energy balance can honour it."""
    ebs = getattr(node, "energy_boundaries", None) or []
    for eb in ebs:
        eb.check()
    if len(ebs) > 1:
        # The node's mixer takes their sum (#1024).
        _check_sensible_fraction(
            f"Node '{node.id}' ({', '.join(repr(eb.id) for eb in ebs)})",
            sum(float(eb.fraction) for eb in ebs),
        )
    targets = [eb for eb in ebs if eb.is_target]
    if not targets:
        return
    if not isinstance(node, _energy_target_node_types()):
        raise ValueError(
            f"EnergyBoundary '{targets[0].id}' sets Tt_target on {type(node).__name__} "
            f"'{node.id}', which cannot honour it: only a PlenumNode, "
            "MomentumChamberNode or CombustorNode derives its Tt from an energy "
            "balance the heat enters."
        )
    if len(targets) > 1:
        raise ValueError(
            f"Node '{node.id}' has {len(targets)} EnergyBoundaries with Tt_target "
            f"({', '.join(repr(eb.id) for eb in targets)}): one node has one outlet "
            "temperature, so at most one boundary may set it."
        )


def _check_sensible_fraction(what: str, fraction: float) -> None:
    """Refuse a sensible-enthalpy fraction at or below -1 (#1024).

    -1 removes all the sensible enthalpy: the node lands on the 298.15 K
    reference whatever it is fed (a set temperature, not a heat removal), and
    a combustor's report divided by 1 + fraction = 0. Below -1 it removes
    more than there is: below the reference temperature, and from a hot
    enough inflow below 0 K, where the solve failed on a raw thermo error.
    A set outlet temperature is EnergyBoundary(Tt_target=...)."""
    if not (math.isfinite(fraction) and fraction > -1.0):
        raise ValueError(
            f"{what}: fraction {fraction!r} removes {-fraction:.0%} of the inflow's "
            "sensible enthalpy (h - h(298.15 K)); at -1 or below that is all of it or "
            "more, a fixed 298.15 K or a temperature below the reference, not a heat "
            "removal. Use a fraction above -1, or Tt_target=... to set the outlet "
            "temperature."
        )


def _attach_energy_boundary(node: Any, eb: EnergyBoundary) -> None:
    """Shared add_energy_boundary: refuse a second Tt_target at once."""
    if eb.is_target and target_energy_boundaries(node):
        raise ValueError(
            f"Node '{node.id}' already has an EnergyBoundary with Tt_target "
            f"('{target_energy_boundaries(node)[0].id}'); '{eb.id}' would set its "
            "outlet temperature a second time."
        )
    node.energy_boundaries.append(eb)


def _energy_target_rows(
    node: Any, state: "NetworkMixtureState", res: list[float], jac: dict[int, dict[str, float]]
) -> None:
    """Append the Tt_target rows (#792): Tt_node - Tt_target = 0 [K].

    The node's Tt is derived, not an unknown: the ``{node}.T`` key is relayed
    by the solver through the propagation chain, whose ``T`` is the mixer's
    (combustor's) output TOTAL temperature -- and with it the column of the
    boundary's own ``{eb}.Q`` unknown, dT_mix/dQ from C++.
    """
    for eb in target_energy_boundaries(node):
        jac[len(res)] = {f"{node.id}.Tt": 1.0}
        res.append(float(state.Tt) - float(eb.Tt_target))


def _applied_heat(node: Any) -> float:
    """Sum of the heat [W] the node's EnergyBoundaries apply (fraction aside)."""
    return sum(eb.heat for eb in getattr(node, "energy_boundaries", None) or [])


@dataclass(init=False)
class NetworkMixtureState:
    """A node's gas state as the residuals see it.

    ``Tt`` is the node's stagnation temperature, the T the solver carries
    (the derived mixing T of an interior node, a boundary's ``Tt``). The
    solver's Jacobian keys ``f"{node}.Tt"`` are derivatives with respect to
    this ``Tt`` (``f"{node}.T"`` until 0.7; #782).

    ``T`` is the STATIC temperature, computed (#782): ``Tt`` at an
    incompressible node (the low-Mach limit its closure solves) and where
    ``P >= Pt`` (a plenum); otherwise the T on the isentrope through
    (``Tt``, ``Pt``) at the static ``P``, which at a converged momentum
    chamber is its closure's static T. No residual reads it
    (``test_network_state_static_t_unread`` keeps it so); results and
    diagnostics do.

    The property methods are evaluated at (``Tt``, static ``P``); where a
    static T is meant (density, cp, transport of a moving gas) that is the
    stagnation value (#824).

    ``Tt``, ``m_dot`` and ``Y`` are keyword-only: the 0.7 signature
    ``(P, Pt, T, Tt, m_dot, Y)`` took a ``T`` that was always ``Tt``, and a
    positional call written for it now fails instead of shifting arguments.
    """

    P: float
    Pt: float
    Tt: float
    m_dot: float
    Y: list[float]
    incompressible: bool

    def __init__(
        self,
        P: float,
        Pt: float,
        *,
        Tt: float,
        m_dot: float = 0.0,
        Y: list[float],
        incompressible: bool = False,
    ) -> None:
        self.P = P
        self.Pt = Pt
        self.Tt = Tt
        self.m_dot = m_dot
        self.Y = Y
        self.incompressible = incompressible

    @property
    def T(self) -> float:
        """Static temperature [K] (see the class docstring)."""
        # A node whose static state is not on its own (Tt, Pt) isentrope: a
        # crossflow station or an area combustor holds Pt = P while its gas
        # moves. The solver sets this for result reporting only (#782 audit).
        override = self.__dict__.get("_T_static")
        if override is not None:
            return float(override)
        Tt, P, Pt = float(self.Tt), float(self.P), float(self.Pt)
        if self.incompressible or not Pt > P or not Tt > 0.0 or not P > 0.0:
            return Tt
        X = self.X
        return float(cb.calc_T_from_s(cb.s(Tt, X, Pt), P, X, Tt, 1e-10))

    @property
    def X(self) -> list[float]:
        """Cached mole fractions converted from mass fractions with defensive clamping."""
        try:
            return self._X  # type: ignore[return-value]
        except AttributeError:
            # Defensive guard against unphysical intermediate solver states:
            # clamp to [0, 1]. Not to a 1e-12 floor, which put a trace of
            # every species into every state and made every property call
            # evaluate all 15 (see normalize_fractions).
            Y_safe = [max(0.0, min(1.0, float(yi))) for yi in self.Y]
            y_sum = sum(Y_safe)
            if y_sum > 0:
                Y_safe = [yi / y_sum for yi in Y_safe]
            x = list(cb.mass_to_mole(Y_safe))
            object.__setattr__(self, "_X", x)
            return x

    def density(self) -> float:
        """Density [kg/m^3] at (Tt, static P)."""
        return cb.density(self.Tt, self.P, self.X)

    def enthalpy(self) -> float:
        """Specific enthalpy [J/kg] at Tt (equals total_enthalpy)."""
        return cb.h_mass(self.Tt, self.X)

    def total_enthalpy(self) -> float:
        """Total (stagnation) specific enthalpy [J/kg]."""
        return cb.h_mass(self.Tt, self.X)

    def cp(self) -> float:
        """Specific heat capacity at constant pressure [J/(kg*K)] at Tt."""
        return cb.cp_mass(self.Tt, self.X)

    def speed_of_sound(self) -> float:
        """Speed of sound [m/s] at Tt."""
        return cb.speed_of_sound(self.Tt, self.X)

    def gamma(self) -> float:
        """Ratio of specific heats (cp/cv) [-] at Tt."""
        cp = self.cp()
        cv = cb.cv_mass(self.Tt, self.X)
        return cp / cv if cv > 0 else 1.4


class _MovingStaticT(NamedTuple):
    """Static T of the gas at a node moving at ``|m| / A`` and its partials:
    d/d(node T0) (the solver's ``{node}.Tt``), d/d(static P), d/d(m) (signed)
    and d/d(Y[k]) for the species that move it."""

    T: float
    dT0: float
    dP: float
    dm: float
    dY: dict[int, float]


def _moving_static_T(
    state: NetworkMixtureState, m_dot: float, area: float, incompressible: bool = False
) -> _MovingStaticT:
    """The static T a kernel documented as STATIC must be given (#782).

    A node carries its stagnation temperature. Where a kernel pairs the node's
    static P with a velocity m / (rho A) of its own, the static T is the
    momentum chamber closure's, h(T0) = h(T) + (G R T / P)^2 / 2 with
    G = |m| / A (#576) -- the same state a MomentumChamberNode solves for
    itself. Incompressible, or with no area, it is the low-Mach limit
    T = T0.
    """
    T0 = float(state.Tt)
    if incompressible or not area > 0.0:
        return _MovingStaticT(T0, 1.0, 0.0, 0.0, {})
    c = _solver_tools.momentum_chamber_residual_and_jacobian(
        state.P, state.Pt, m_dot, T0, state.Y, area, False
    )
    dY = {k: float(d) for k, d in enumerate(c.dT_static_dY) if d != 0.0}
    return _MovingStaticT(
        float(c.T_static),
        float(c.dT_static_dT),
        float(c.dT_static_dP),
        float(c.dT_static_dmdot),
        dY,
    )


class _MovingTotalP(NamedTuple):
    """Stagnation pressure of gas at static P moving at ``|m| / A`` and its
    partials: d/d(static P), d/d(T0), d/d(m) (signed) and d/d(Y[k])."""

    Pt: float
    dP: float
    dT0: float
    dm: float
    dY: dict[int, float]


def _moving_total_P(
    P: float, T0: float, Y: list[float], m_dot: float, area: float
) -> _MovingTotalP:
    """The stagnation pressure of a duct's static state (#802).

    The momentum chamber closure (#357, #576), as ``_moving_static_T``: the
    static T from h(T0) = h(T) + (G R T / P)^2 / 2, G = |m| / A, then
    Pt = P exp((s(T0) - s(T)) / R). C++ returns the closure as a residual
    P_total - P0 with its analytic partials; P0 is read with P_total = 0 and
    the partials are the residual's, negated.
    """
    c = _solver_tools.momentum_chamber_residual_and_jacobian(P, 0.0, m_dot, T0, Y, area, False)
    return _MovingTotalP(
        float(P) + float(c.q_dyn),
        -float(c.d_res_dP),
        -float(c.d_res_dT),
        -float(c.d_res_dmdot),
        {k: -float(d) for k, d in enumerate(c.d_res_dY) if d != 0.0},
    )


def _station_node_partials(r: "_StaticChainedDrop", node: str) -> list[tuple[str, float]]:
    """A station term's partials in its node's P, T (T0) and Y, keyed by name."""
    out = [(f"{node}.P", r.d_dP), (f"{node}.Tt", r.d_dT)]
    out.extend((f"{node}.Y[{k}]", d) for k, d in r.d_dY.items())
    return out


class _StaticChainedDrop:
    """A ``StationHalfDrop`` evaluated at a node's moving static T, its partials
    chained onto what that T depends on (#782): the node's T0, its static P,
    its composition, and the flow the static state was taken at -- the mean
    of the station's arriving (a) and leaving (b) flows when ``two_sided``,
    else the arriving flow alone."""

    __slots__ = ("dP", "d_dm_a", "d_dm_b", "d_dP", "d_dT", "d_dY", "extra")

    def __init__(self, r: Any, st: _MovingStaticT, two_sided: bool) -> None:
        #: Partials by unknown name beyond the station's own flows and node
        #: state: the joining jets' axial momentum (#813).
        self.extra: dict[str, float] = {}
        dm = float(r.d_dT) * st.dm
        self.dP = float(r.dP)
        self.d_dm_a = float(r.d_dm_a) + (0.5 * dm if two_sided else dm)
        self.d_dm_b = float(r.d_dm_b) + (0.5 * dm if two_sided else 0.0)
        self.d_dP = float(r.d_dP) + float(r.d_dT) * st.dP
        self.d_dT = float(r.d_dT) * st.dT0
        self.d_dY = {k: float(r.d_dT) * d for k, d in st.dY.items()}


class _ReservoirExitInflow:
    """Incompressible acceleration out of a segment's exit node at rest when
    the flow reverses: ``channel_entry_drop`` evaluated at -m (K_in = 0), as
    a momentum term of this segment's signed m. The static P at the duct end
    is the node's P less the head, so the term is a negative drop."""

    __slots__ = ("dP", "d_dm_a", "d_dm_b", "d_dP", "d_dT", "d_dY", "extra")

    def __init__(self, r: Any) -> None:
        self.extra: dict[str, float] = {}
        self.dP = -float(r.dP)
        self.d_dm_a = float(r.d_dm_a)  # d(-q(-m))/dm = q'(-m)
        self.d_dm_b = 0.0
        self.d_dP = -float(r.d_dP)
        self.d_dT = -float(r.d_dT)
        self.d_dY: dict[int, float] = {}


#: Flow regime every element defaults to when none is given. It is the GUI's
#: default too (gui/backend/schemas.py SolverSettings.global_regime and the
#: frontend store), and test_regime_default_sync.py keeps the three equal: a
#: script and a drawing of the same network must solve the same physics.
#: Compressible since #585: incompressible flow is its low-Mach limit, so the
#: default no longer has to be chosen per network.
DEFAULT_REGIME: Literal["incompressible", "compressible"] = "compressible"

#: Nominal cross-section for a chamber-like node whose area could not be
#: inferred and which nothing in the network measures against. It is a
#: bookkeeping placeholder, not a geometry: components for which the area IS
#: the physics (AreaChangeElement, ChannelElement, TeeJunctionElement) refuse
#: to default and raise instead. Any consumer of a defaulted area reports
#: ``area_source = "default"`` in its diagnostics (issue #262).
DEFAULT_CHAMBER_AREA = 0.1

#: Diameter [m] a ChannelElement falls back to when it opts in
#: (``default_diameter``) and neither its own diameter nor any neighbour's
#: area is known (#570). The GUI opts in for every channel whose D is unset;
#: the Python API does not, so a script keeps #262's refusal. Reported as
#: ``area_source = "default"``.
DEFAULT_CHANNEL_DIAMETER = 0.1

#: Bore diameter [m] an orifice falls back to when the user set neither its
#: diameter nor its area. Unlike a chamber area, a bore IS the physics -- it
#: should be set explicitly -- but an orifice is a restriction, not a channel
#: area a neighbour can measure against, so there is nothing to infer; the
#: default keeps the network solvable and is reported as ``area_source =
#: "default"`` (issue #559).
DEFAULT_ORIFICE_BORE = 0.08

#: Relative tolerance on "orifice bore == pipe diameter" (#788). A bore that
#: matches the pipe it sits in to within this is NO orifice: zero loss, the
#: flow passes unchanged (beta = 1). A bore larger than the pipe by more than
#: this is refused. Wide enough to absorb a diameter round-tripped through an
#: area or a unit conversion; at beta = 1 - 1e-6 an in-pipe orifice's
#: permanent loss is nil anyway, so nothing physical hides inside it. The
#: same tolerance compares the geometric bore (a plate's total hole area)
#: with a momentum chamber it opens into or draws from (#1016), on the
#: equivalent diameter sqrt(A_bore / A_chamber).
ORIFICE_BORE_PIPE_RTOL = 1e-6

#: The Cd an orifice-family start guess assumes when a correlation computes
#: the Cd (#872): the value the unset Cd defaulted to before 0.8.0, so every
#: start guess is unchanged. Never a result.
_ORIFICE_SEED_CD = 0.6

#: Relative tolerance on "an element's area equals its neighbour's" (#1023):
#: #788's bore-pipe band, compared on the equivalent diameter
#: sqrt(A / A_neighbour). Equal within it is a pass-through (an area change
#: with F0 == F1) or a consistent face (an explicit F0, F1, F_C or F_branch
#: against the pipe it is attached to); outside it an explicit face is
#: refused, as is an ejector mixing area above its outlet pipe, and a
#: momentum chamber below an adjoining pipe is warned about.
NEIGHBOUR_AREA_RTOL = ORIFICE_BORE_PIPE_RTOL

#: Relative tolerance of the isoperimetric bound on a hydraulic diameter
#: (#1024). Dh = 4 A / P, and no section of area A has a shorter perimeter
#: than the circle, so Dh <= 2 sqrt(A / pi), the circle's diameter. A Dh above
#: it by more than this is refused: no section has it. The band absorbs an
#: area typed to three significant figures (a 20 mm circle as 3.14e-4 m^2 is
#: 2.6e-4 over); an impossible Dh is off by tens of percent. Checked on
#: MomentumChamberNode and CombustorNode, whose section is closed.
ISOPERIMETRIC_DH_RTOL = 1e-3


def _positive_finite(value: float | None) -> bool:
    """True for a number that is finite and above zero."""
    return value is not None and math.isfinite(value) and value > 0.0


def _section_error(
    what: str, area: float | None, area_given: bool, Dh: float | None, Dh_given: bool
) -> str:
    """Why a node's flow section cannot be built (#1024), or "".

    ``what`` names the component. A given area must be positive and finite
    (zero is no section; leave it unset to inherit one), a given Dh the same,
    and with a positive area a given Dh must satisfy the isoperimetric bound
    Dh <= 2 sqrt(A / pi) (ISOPERIMETRIC_DH_RTOL)."""
    if area_given and not _positive_finite(area):
        return (
            f"{what}: area = {area!r} m^2. A flow section needs a positive, finite area; "
            "leave it unset to inherit one from a neighbour."
        )
    if Dh_given and not _positive_finite(Dh):
        return f"{what}: hydraulic diameter Dh = {Dh!r} m must be positive and finite."
    if Dh_given and _positive_finite(area):
        D_max = 2.0 * math.sqrt(float(area) / math.pi)
        if float(Dh) > D_max * (1.0 + ISOPERIMETRIC_DH_RTOL):
            return (
                f"{what}: hydraulic diameter Dh = {float(Dh):.6g} m is larger than any "
                f"section of area A = {float(area):.6g} m^2 allows, 2 sqrt(A/pi) = "
                f"{D_max:.6g} m (Dh = 4A/P, and the circle has the shortest perimeter P "
                "for its area). Lower Dh or enlarge the area."
            )
    return ""


def _hole_spacing_error(
    what: str,
    d: float,
    angle_deg: float,
    spanwise: tuple[str, float],
    streamwise: tuple[str, float] | None,
) -> str:
    """Why rows of holes of diameter ``d`` overlap at these pitches (#1024),
    or "": INLINE rows, the arrangement assumed where none is given.

    Holes must leave material between them: neighbours along a row are the
    spanwise pitch apart and must clear d; consecutive rows are the
    streamwise pitch apart and must clear the hole's footprint on the
    surface, d / sin(angle) (an inclined hole's ellipse is that long).
    Equal is refused: touching holes leave no plate. ``streamwise`` None:
    one row, nothing to clear."""
    name, pitch = spanwise
    if not pitch > d:
        return (
            f"{what}: spanwise {name} = {pitch:.6g} m is not larger than the hole "
            f"diameter d = {d:.6g} m, so neighbouring holes in a row touch or overlap "
            f"and leave no material between them. {name} must exceed d."
        )
    if streamwise is None:
        return ""
    name, pitch = streamwise
    footprint = d / math.sin(math.radians(angle_deg))
    if not pitch > footprint:
        return (
            f"{what}: streamwise {name} = {pitch:.6g} m is not larger than the holes' "
            f"footprint on the surface, d / sin(angle) = {footprint:.6g} m (d = {d:.6g} "
            f"m at {angle_deg:g} deg), so consecutive rows touch or overlap and leave "
            f"no material between them (rows taken inline). {name} must exceed "
            f"{footprint:.6g} m."
        )
    return ""


#: Mass flow [kg/s] over which a directional loss blends its forward and
#: reverse forms, logistic in m_dot / scale (C++ ``sigmoid`` in
#: area_change.cpp). Outside a few scales either form holds alone.
FLOW_REVERSAL_M_SCALE = 1e-4


#: Pressure [Pa] over which an orifice blends its forward and reverse forms,
#: logistic in (Pt_from - Pt_to) / scale (#554). Outside ~30 scales -- 3 Pa --
#: either form holds alone, so an orifice with any real drop is unchanged.
ORIFICE_REVERSAL_P_SCALE = 0.1

#: Relative width of the smooth max(P_boundary, P*) a port discharge holds
#: its chamber at (#637): sqrt(dP^2 + (w P_b)^2). The value moves by at most
#: w P_b / 2 (5 Pa at 1 bar), and only where the exit is exactly at choke.
PORT_DISCHARGE_CHOKE_BLEND = 1e-4

#: Omitted Rayleigh total-pressure change, relative to the inlet Pt, above
#: which a wall-heated compressible channel is flagged (#551): the momentum
#: kernel is adiabatic Fanno flow, so heating's cost is not in the solve.
RAYLEIGH_OMISSION_FLAG = 1e-3

#: Relay key prefix of a wall-heated duct's own exit state (#672): its stream
#: and its own wall heat, mixed alone where the node it feeds also takes other
#: streams or heat. No "." in it: the Jacobian splits keys at the first dot.
RAYLEIGH_EXIT_PREFIX = "__rayleigh_exit__"

#: Ratio of a smooth-wall channel's heat-transfer area to its friction-wetted
#: area 4 A L / Dh above which the channel is flagged (#551): heat then flows
#: over surface its friction never sees, which breaks Reynolds analogy -- in a
#: compressible duct a cooling gain can then outgrow the friction, making the
#: duct a pump at fixed exit Tt. A stress-harness crossflow segment carried
#: 3.9x (the plate's width on a round duct) and could not be solved.
HEAT_AREA_EXCESS_FLAG = 1.1


_ROW_KINDS = ("p", "mdot")


def _defining_class(cls: type, attr: str) -> type | None:
    return next((k for k in cls.__mro__ if attr in k.__dict__), None)


def declared_row_kinds(element: "NetworkElement") -> list[str]:
    """The element's residual row kinds, one per row (RESIDUAL-ROW-SCALING).

    Raises TypeError unless the units are declared where the residual is:
    ``residual_scale_kind`` (``row_scale_kinds()`` for a junction) must be
    defined on the class that defines ``residuals`` or on a subclass of it.
    A subclass that rewrites its residual therefore cannot inherit its
    parent's units unseen; it re-declares them, if only to confirm them.
    """
    cls = type(element)
    junction = isinstance(element, MultiPortChamberBase)
    attr = "row_scale_kinds" if junction else "residual_scale_kind"
    res_owner = _defining_class(cls, "residuals")
    kind_owner = _defining_class(cls, attr)
    if (
        kind_owner is None
        or kind_owner is NetworkElement
        or not (res_owner is not None and issubclass(kind_owner, res_owner))
    ):
        raise TypeError(
            f"{cls.__name__} ({element.id!r}): declare the units of its residual rows "
            f"with {attr} on {res_owner.__name__ if res_owner else cls.__name__} or a "
            f"subclass (found on {kind_owner.__name__ if kind_owner else 'none'}); "
            "see RESIDUAL-ROW-SCALING in NetworkSolver._build_residual_scales."
        )
    kinds = (
        list(element.row_scale_kinds())
        if junction
        else [element.residual_scale_kind] * element.n_equations()
    )
    bad = [k for k in kinds if k not in _ROW_KINDS]
    if bad:
        raise TypeError(
            f"{cls.__name__} ({element.id!r}): row kinds must be 'p' or 'mdot', got {bad}"
        )
    return kinds


def _pressure_blend(dpt: float) -> tuple[float, float]:
    """(w, dw/d(dpt)): the forward weight of an orifice, 1 for Pt_from >> Pt_to."""
    x = dpt / ORIFICE_REVERSAL_P_SCALE
    if x > 30.0:
        return 1.0, 0.0
    if x < -30.0:
        return 0.0, 0.0
    w = 1.0 / (1.0 + math.exp(-x))
    return w, w * (1.0 - w) / ORIFICE_REVERSAL_P_SCALE


def _flow_blend(m_dot: float, scale: float = FLOW_REVERSAL_M_SCALE) -> tuple[float, float]:
    """(w, dw/dm): the forward weight of a directional loss, 1 for m >> scale."""
    x = m_dot / scale
    if x > 30.0:
        return 1.0, 0.0
    if x < -30.0:
        return 0.0, 0.0
    w = 1.0 / (1.0 + math.exp(-x))
    return w, w * (1.0 - w) / scale


# ---------------------------------------------------------------------------
# Flow-area inference
# ---------------------------------------------------------------------------
# Elements and nodes that need a reference flow area infer one from the
# topology when the user did not supply it.  Historically every such site
# searched for a neighbouring ChannelElement only and, failing that, fell back
# to a hard-coded constant (0.01 / 0.02 / 0.1 m^2).  Those constants trace to
# "something had to be there", not to any geometry: with a non-channel
# neighbour -- a momentum chamber, a combustor, an ejector outlet -- the
# fallback fabricated an area change that never existed, and the resulting
# residual is unsatisfiable at the real mass flow, so the solver stalls with a
# message that blames the solver rather than the geometry (issue #262).
#
# _infer_flow_area consults every area-bearing neighbour instead.  Channels
# keep the highest priority so that networks which resolved before resolve to
# the same area now.


def _known_node_area(node: "NetworkNode | None") -> float | None:
    """Flow area of a node, but only when it is genuinely known.

    A node whose area was auto-sized rather than user-supplied
    (``_auto_area``) carries a placeholder, not a measurement; inferring from
    it would launder one guess into another.
    """
    if node is None or getattr(node, "_auto_area", True):
        return None
    area = getattr(node, "area", None)
    return float(area) if area is not None and area > 0.0 else None


def _element_area_facing(elem: "NetworkElement", node_id: str) -> float | None:
    """Flow area the element presents to ``node_id``, or None if it has none.

    Orifices are deliberately excluded: a bore is a restriction, not the
    channel area a neighbour should inherit.
    """
    if isinstance(elem, ChannelElement):
        if elem.diameter is not None:
            return math.pi * (elem.diameter / 2.0) ** 2
        return None
    if isinstance(elem, AreaChangeElement):
        # The element changes area across itself, so the face matters.
        if elem.to_node == node_id:
            return float(elem.F1) if elem.F1 is not None else None
        if elem.from_node == node_id:
            return float(elem.F0) if elem.F0 is not None else None
        return None
    if isinstance(elem, TeeJunctionElement):
        if elem.F_C is None:
            return None
        # Branch arm carries F_C/psi; straight and common arms carry F_C.
        if node_id == elem.branch_node:
            return float(elem.F_C) / elem.psi
        return float(elem.F_C)
    if isinstance(elem, PressureLossElement):
        return float(elem.area) if elem.area is not None else None
    # Elements defined outside this module (e.g. the ejector's outlet duct)
    # declare their own faces.
    facing = getattr(elem, "area_facing", None)
    if callable(facing):
        return facing(node_id)
    return None


def _infer_flow_area(
    graph: "FlowNetwork",
    node_ids: "list[str]",
    exclude: "NetworkElement | None" = None,
) -> "tuple[float | None, str]":
    """Infer a flow area from the neighbourhood of ``node_ids``.

    Returns ``(area, provenance)``; ``area`` is None when nothing in the
    neighbourhood knows one, and ``provenance`` names the source so callers
    can report where a resolved area came from.

    Channels are searched first across all of ``node_ids`` so that a network
    which previously resolved from a channel still resolves to that channel's
    area regardless of what else is attached.
    """
    neighbours: list[tuple[str, NetworkElement]] = []
    for node_id in node_ids:
        for elem in graph.get_upstream_elements(node_id) + graph.get_downstream_elements(node_id):
            if elem is exclude:
                continue
            neighbours.append((node_id, elem))

    for node_id, elem in neighbours:
        if isinstance(elem, ChannelElement):
            area = _element_area_facing(elem, node_id)
            if area is not None:
                return area, f"channel '{elem.id}'"

    for node_id, elem in neighbours:
        area = _element_area_facing(elem, node_id)
        if area is not None:
            return area, f"{type(elem).__name__} '{elem.id}'"

    for node_id in node_ids:
        area = _known_node_area(graph.nodes.get(node_id))
        if area is not None:
            return area, f"node '{node_id}'"

    for node_id, elem in neighbours:
        for other in dict.fromkeys([*elem.all_source_nodes(), *elem.all_sink_nodes()]):
            if other == node_id:
                continue
            area = _known_node_area(graph.nodes.get(other))
            if area is not None:
                return area, f"node '{other}'"

    return None, ""


def _unresolved_area_error(
    component: "NetworkNode | NetworkElement",
    parameter: str,
    node_ids: "list[str]",
) -> ValueError:
    """The error raised when a load-bearing flow area cannot be inferred.

    Names the component, the parameter that clears the error, and where the
    search looked -- a defaulted area here fabricates geometry, so failing
    loudly is the only honest option (issue #262).
    """
    where = ", ".join(f"'{n}'" for n in node_ids)
    return ValueError(
        f"{type(component).__name__} '{component.id}': cannot determine "
        f"{parameter} -- no neighbouring channel, area-bearing element or "
        f"node with a known area was found at {where}. Set {parameter} "
        f"explicitly on this component (or give a neighbour a known area). "
        f"A default would fabricate a geometry that does not exist and the "
        f"solver would stall on it."
    )


def _attached_pipe(
    graph: "FlowNetwork", node_id: str, exclude: "NetworkElement"
) -> "tuple[ChannelElement, float] | None":
    """The pipe an element's face at ``node_id`` is attached to, with its area
    (#1023), or None.

    The single ChannelElement at the node other than ``exclude``, in either
    direction, whose diameter states the geometry: never the GUI's
    DEFAULT_CHANNEL_DIAMETER placeholder. With no channel, or more than one
    (a node that joins pipes of their own sizes), the face has no one pipe
    to agree with.
    """
    pipes = [
        e
        for e in graph.get_upstream_elements(node_id) + graph.get_downstream_elements(node_id)
        if e is not exclude and isinstance(e, ChannelElement)
    ]
    if len(pipes) != 1:
        return None
    pipe = pipes[0]
    if (
        not (pipe.diameter and pipe.diameter > 0.0)
        or getattr(pipe, "_area_source", "") == "default"
    ):
        return None
    return pipe, math.pi * (pipe.diameter / 2.0) ** 2


def _chamber_states_own_area(node: "NetworkNode | None") -> bool:
    """A momentum chamber or combustor at a face whose area is the user's own
    (an area or a Dh): the face then meets that chamber, not the pipe beyond
    it, and a chamber against its pipe is MomentumChamberNode's warning
    (#1023)."""
    return isinstance(node, (MomentumChamberNode, CombustorNode)) and node._area_source in (
        "user",
        "derived from Dh",
    )


def _check_face_against_pipe(
    element: "NetworkElement",
    graph: "FlowNetwork",
    name: str,
    area: float,
    node_ids: "list[str]",
) -> None:
    """Refuse an explicit face area ``name`` of ``element`` that contradicts
    the pipe attached at one of ``node_ids`` (#1023): their equivalent
    diameters differ by more than NEIGHBOUR_AREA_RTOL. Equal is consistent.
    A node whose chamber states its own area is skipped
    (_chamber_states_own_area)."""
    for node_id in dict.fromkeys(node_ids):
        if _chamber_states_own_area(graph.nodes.get(node_id)):
            continue
        hit = _attached_pipe(graph, node_id, element)
        if hit is None:
            continue
        pipe, A_pipe = hit
        if abs(math.sqrt(area / A_pipe) - 1.0) <= NEIGHBOUR_AREA_RTOL:
            continue
        raise ValueError(
            f"{type(element).__name__} '{element.id}': {name} = {area:.6g} m^2 contradicts "
            f"the channel '{pipe.id}' it is attached to at '{node_id}' "
            f"(D = {pipe.diameter:.6g} m, A = {A_pipe:.6g} m^2). The face is that pipe's "
            f"own: leave {name} unset to inherit it, or make the two equal. A change of "
            "area between them is an AreaChangeElement of its own."
        )


class NetworkNode(ABC):
    #: True when this node supplies a meaningful temperature-rise ratio theta
    #: (T_burned / T_unburned - 1).  Read by ``PressureLossElement`` with
    #: theta-aware correlations.  Overridden on combustion-capable nodes
    #: (``CombustorNode``).  Nodes advertising ``has_theta = True`` must
    #: expose ``_T_unburned`` and ``_T_burned`` float attributes after
    #: ``compute_derived_state`` has been called.
    has_theta: bool = False

    #: The network ends here, so the node may be a dead end: one connected
    #: element, or flow in one direction only. The single definition graph
    #: validation, the solver and the GUI read (#545).
    is_boundary: bool = False

    #: Mass and enthalpy cross the network's control volume here, at the state
    #: the node imposes. Every flow boundary is a boundary; a WallNode is a
    #: boundary nothing crosses. A MassFlowBoundary between two elements is
    #: still a flow boundary: it injects its stream into the network.
    is_flow_boundary: bool = False

    #: The node's own definition holds its gas at rest, Pt = P: a plenum, a
    #: combustor, a wall node, and the boundaries (a reservoir at an inflow,
    #: a static exit pressure at an outflow). Its state is not the state in a
    #: duct meeting it, so a ChannelElement reports its own duct state at
    #: such an end (#888). A momentum chamber carries its moving static state
    #: (Pt > P) and keeps False.
    state_at_rest: bool = False

    def __init__(self, id: str):
        self.id = id

    @property
    def has_convective_surface(self) -> bool:
        """True if this node has a convective heat transfer surface."""
        return False

    @abstractmethod
    def unknowns(self) -> list[str]:
        """Names of unknowns this node contributes to the solver."""
        pass

    @abstractmethod
    def residuals(
        self, state: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        """Returns (residuals, local_jacobian)."""
        pass

    def compute_derived_state(
        self, upstream_states: list[NetworkMixtureState]
    ) -> tuple[float, list[float], Any]:
        """
        Computes (Tt, Y, Jac) for nodes based on upstream conditions.
        Jac contains derivatives d(Tt)/d(stream_i) and d(Y)/d(stream_i).
        By default, it just takes the first upstream state (pass-through).
        """
        if not upstream_states:
            return 300.0, list(cb.mole_to_mass(cb.species.dry_air())), None

        up = upstream_states[0]
        return up.Tt, up.Y, None

    @abstractmethod
    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """Called automatically by FlowNetwork to resolve neighbors."""
        pass

    def mach(self, state: NetworkMixtureState) -> float:
        """Returns the local Mach number at this node. Default 0.0 (stagnation)."""
        return 0.0

    def diagnostics(self, state: NetworkMixtureState) -> dict[str, float]:
        """Compute generalized node-level diagnostics (e.g., thermo properties)."""
        if state.P <= 0:
            return {}

        try:
            cs = cb.complete_state(state.T, state.P, state.X)
            return {
                "Tt": float(state.Tt),
                "Pt": float(state.Pt),
                "h": cs.thermo.h,
                "s": cs.thermo.s,
                "u": cs.thermo.u,
                "rho": cs.thermo.rho,
                "gamma": cs.thermo.gamma,
                "a": cs.thermo.a,
                "cp": cs.thermo.cp,
                "cv": cs.thermo.cv,
                "mw": cs.thermo.mw,
                "mu": cs.transport.mu,
                "k": cs.transport.k,
                "Pr": cs.transport.Pr,
                "nu": cs.transport.nu,
            }
        except Exception:
            ts = cb.thermo_state(state.T, state.P, state.X)
            return {
                "Tt": float(state.Tt),
                "Pt": float(state.Pt),
                "h": ts.h,
                "s": ts.s,
                "u": ts.u,
                "rho": ts.rho,
                "gamma": ts.gamma,
                "a": ts.a,
                "cp": ts.cp,
                "cv": ts.cv,
                "mw": ts.mw,
            }


# ---------------------------------------------------------------------------
# Shared diagnostic helpers
# ---------------------------------------------------------------------------


def _element_pressure_block(
    state_in: "NetworkMixtureState",
    state_out: "NetworkMixtureState",
    mach_in: float,
    mach_out: float,
) -> dict[str, float]:
    """Return the universal 14-field pressure/Mach block for any element."""
    p_ratio_total = state_out.Pt / state_in.Pt if state_in.Pt > 0 else 0.0
    p_ratio_static = state_out.P / state_in.P if state_in.P > 0 else 0.0
    return {
        "P_in": float(state_in.P),
        "P_out": float(state_out.P),
        "T_in": float(state_in.T),
        "T_out": float(state_out.T),
        "Pt_in": float(state_in.Pt),
        "Pt_out": float(state_out.Pt),
        "Tt_in": float(state_in.Tt),
        "Tt_out": float(state_out.Tt),
        "mach_in": float(mach_in),
        "mach_out": float(mach_out),
        "dPt": float(state_in.Pt - state_out.Pt),
        "dP": float(state_in.P - state_out.P),
        "pr_total": float(p_ratio_total),
        "pr_static": float(p_ratio_static),
    }


def _element_reference_block(
    cs: Any,
    *,
    m_dot: float,
    area: float,
    Dh: float | None = None,
    location: str = "inlet",
) -> dict[str, float | str]:
    """Return the thermo/transport reference-state block evaluated at the inlet.

    All keys use the ``_ref`` suffix so it is unambiguous which state they
    refer to.  ``velocity`` and ``Re`` are always included; Re is 0 when
    geometry is insufficient to compute it.
    """
    rho = cs.thermo.rho
    v = (abs(m_dot)) / ((rho) * (area)) if area > 0 and rho > 0 else 0.0
    mu = cs.transport.mu
    re = (rho) * (v) * (Dh) / (mu) if (Dh and Dh > 0 and v > 0 and mu > 0) else 0.0
    return {
        "ref_location": location,
        "velocity": v,
        "Re": re,
        "rho_ref": float(rho),
        "h_ref": float(cs.thermo.h),
        "s_ref": float(cs.thermo.s),
        "u_ref": float(cs.thermo.u),
        "cp_ref": float(cs.thermo.cp),
        "cv_ref": float(cs.thermo.cv),
        "gamma_ref": float(cs.thermo.gamma),
        "a_ref": float(cs.thermo.a),
        "mw_ref": float(cs.thermo.mw),
        "mu_ref": float(mu),
        "k_ref": float(cs.transport.k),
        "nu_ref": float(cs.transport.nu),
        "Pr_ref": float(cs.transport.Pr),
    }


class LazyJacobian(MutableMapping):
    """An element's local Jacobian, computed on first access (#489).

    The solver's residual-only evaluations never read an element's Jacobian;
    it is assembled only when the root finder asks for J (hybr: at the start
    and on restarts). An element whose derivatives are costly -- the
    compressible channel's implicit Fanno derivatives -- returns this from
    ``residuals`` instead of a dict, and pays for them only when they are
    read. To any reader it is the ``{eq_idx: {unknown: d}}`` dict.
    """

    def __init__(self, compute: Callable[[], dict[int, dict[str, float]]]) -> None:
        self._compute: Callable[[], dict[int, dict[str, float]]] | None = compute
        self._jac: dict[int, dict[str, float]] = {}

    def _resolved(self) -> dict[int, dict[str, float]]:
        if self._compute is not None:
            compute, self._compute = self._compute, None
            self._jac = compute()
        return self._jac

    def __getitem__(self, key: int) -> dict[str, float]:
        return self._resolved()[key]

    def __setitem__(self, key: int, value: dict[str, float]) -> None:
        self._resolved()[key] = value

    def __delitem__(self, key: int) -> None:
        del self._resolved()[key]

    def __iter__(self) -> Iterator[int]:
        return iter(self._resolved())

    def __len__(self) -> int:
        return len(self._resolved())


class NetworkElement(ABC):
    #: Diagnostics keys holding the element's highest internal Mach numbers,
    #: for the incompressible validity report (#588). Empty: the element has
    #: no incompressible closure, or does not report one.
    validity_mach_keys: tuple[str, ...] = ()

    #: Units of the residual rows, "p" (Pa) or "mdot" (kg/s), for the solver's
    #: row scaling (RESIDUAL-ROW-SCALING, see NetworkSolver._build_residual_scales).
    #: REQUIRED, with no default: declare it on the class that defines
    #: ``residuals`` (a class attribute, or a property when the units follow
    #: the configuration). Junctions declare ``row_scale_kinds()`` instead.
    #: ``declared_row_kinds`` enforces both at network validation: a silent
    #: "mdot" default misfiled pressure rows 1e5-1e9 too heavily three times
    #: (MPCE impulse rows, the compressible channel #481, four Pa elements #614).
    residual_scale_kind: str | None = None

    #: Start a compressible ``"default"`` solve from the outlet-referenced
    #: incompressible proxy instead of cold, falling back to the cold path if
    #: that fails (#706). Declared only where it was measured to pay: a chain
    #: of compressible momentum stations converges slowly from the cold
    #: start, which conditions those coupled rows poorly.
    prefers_warm_start: bool = False

    def __init__(self, id: str, from_node: str, to_node: str):
        self.id = id
        self.from_node = from_node
        self.to_node = to_node

    @property
    def has_convective_surface(self) -> bool:
        """True if this element has a convective heat transfer surface."""
        return False

    @abstractmethod
    def unknowns(self) -> list[str]:
        """Names of unknowns this element contributes to the solver."""
        pass

    @abstractmethod
    def residuals(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        """
        Returns (residuals, local_jacobian).
        local_jacobian is a dict mapping equation_index to a dict of {unknown_name: partial_derivative},
        or a LazyJacobian that computes it when first read.
        """
        pass

    @abstractmethod
    def n_equations(self) -> int:
        """Number of residual equations contributed."""
        pass

    @abstractmethod
    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """Called automatically by FlowNetwork to resolve neighbors."""
        pass

    def all_source_nodes(self) -> list[str]:
        """All node IDs this element draws flow FROM."""
        return [self.from_node]

    def all_sink_nodes(self) -> list[str]:
        """All node IDs this element delivers flow TO."""
        return [self.to_node]

    def flow_at_node(self, node_id: str, x: Any, indices: list[int]) -> float:
        """Mass flow this element contributes at the given connected node."""
        return float(x[indices[0]])

    def flow_jac_at_node(self, node_id: str, indices: list[int]) -> dict[int, float]:
        """d(mass_flow_at_node)/d(solver_unknown): {global_index: coefficient}."""
        return {indices[0]: 1.0}

    def _seed_resistance(self, rho: float, P_ref: float, m_ref: float) -> float | None:
        """START GUESS ONLY: K in dPt ~ K m|m| [Pa s^2/kg^2], an incompressible
        estimate at the reference density, for the solver's resistance-network
        cold start (``NetworkSolver._resistance_network_guess``). Never used in
        a residual. None (the default): no estimate, and a network holding the
        element keeps the propagated start."""
        return None

    def _seed_flow_applies(self) -> bool:
        """START GUESS ONLY: whether ``_seed_flow`` re-fits this element's
        seed resistance to its own compressible kernel (#971)."""
        return False

    def _seed_flow(
        self, Pt_up: float, Pt_dn: float, Tt: float, X: list[float], Y: list[float]
    ) -> float | None:
        """START GUESS ONLY: |m_dot| [kg/s] the element passes from Pt_up to
        Pt_dn at Tt (mole fractions X, mass fractions Y); None for no
        estimate."""
        return None

    def validate(self) -> None:
        """Perform element-specific validation checks.
        Should raise ValueError with a clear message if validation fails.
        """
        return None

    def htc_and_T(self, state: NetworkMixtureState):
        """Compute heat transfer coefficient and adiabatic wall temperature."""
        return None

    def diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float]:
        """Compute diagnostic properties for this element (e.g. Mach, P_ratio)."""
        return {}


class PlenumNode(NetworkNode):
    """
    A stagnation volume where v ~ 0, so Pt = P_static.
    For Phase 1: pressure-only node.
    For Phase 2+: automatically handles mixing when multiple upstream connections exist.
    """

    state_at_rest: bool = True  # Pt = P by definition

    def __init__(self, id: str):
        super().__init__(id)
        self.upstream_elements = []
        self.energy_boundaries: list[EnergyBoundary] = []

    def add_energy_boundary(self, eb: EnergyBoundary) -> None:
        """Attach an energy source/sink to this plenum."""
        _attach_energy_boundary(self, eb)

    def unknowns(self) -> list[str]:
        # Pure Pressure-Flow: Temperature and Composition are derived forward, not unknowns.
        # A Tt_target EnergyBoundary adds its heat flow (#792).
        return [f"{self.id}.P", f"{self.id}.Pt"] + [
            eb.unknown_name for eb in target_energy_boundaries(self)
        ]

    def compute_derived_state(
        self, upstream_states: list[NetworkMixtureState]
    ) -> tuple[float, list[float], Any]:
        """Derived T and Y for a plenum (simple mixing + energy boundaries)."""

        # Each inflow's supply state: a crossflow station reads the axial
        # momentum its joining jets bring from it (#813), as a merge chamber
        # does.
        self._inflows = [(getattr(s, "_element_id", None), s) for s in upstream_states]

        if not upstream_states:
            return 300.0, list(cb.mole_to_mass(cb.species.dry_air())), None

        streams = [cb.MassStream(s.m_dot, s.Tt, s.Pt, s.Y) for s in upstream_states]
        Q_total = _applied_heat(self)
        fraction_total = sum(eb.fraction for eb in self.energy_boundaries)

        mix_res = _solver_tools.mixer_from_streams_and_jacobians(
            streams, Q=Q_total, fraction=fraction_total
        )

        return mix_res.Tt_mix, mix_res.Y_mix, mix_res

    def residuals(
        self, state: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Base residual: Pt = P for plenum
        res = [state.Pt - state.P]
        jac = {0: {f"{self.id}.P": -1.0, f"{self.id}.Pt": 1.0}}
        _energy_target_rows(self, state, res, jac)
        return res, jac

    def diagnostics(self, state: NetworkMixtureState) -> dict[str, float]:
        # Plenums: v ~ 0, so stagnation = static
        ts = cb.thermo_state(state.T, state.P, state.X)
        return {
            "Tt": float(state.Tt),
            "Pt": float(state.P),
            "h": ts.h,
            "s": ts.s,
            "u": ts.u,
            "rho": ts.rho,
            "gamma": ts.gamma,
            "a": ts.a,
            "cp": ts.cp,
            "cv": ts.cv,
            "mw": ts.mw,
        }

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Store upstream elements for mixing calculations
        self.upstream_elements = graph.get_upstream_elements(self.id)


class MomentumChamberNode(NetworkNode):
    """
    Similar to a plenum, but preserves dynamic pressure.
    Enforces Pt = P_static + 0.5 * rho * v^2 instead of Pt = P_static.
    Supports directional flow vectors (port angles) for momentum conservation.
    Automatically handles mixing when multiple upstream connections exist.
    """

    # to_dict serialises the given geometry, not what resolve_topology
    # inherited into the same attributes (#1099).
    _init_kwarg_attrs = {
        "area": lambda n: None if n._auto_area else n.area,
        "Dh": lambda n: n.Dh if n._Dh_source == "user" else None,
    }

    def __init__(
        self,
        id: str,
        area: float | None = None,
        port_angles_deg: dict[str, float] | None = None,
        length: float | None = None,
        surface: ConvectiveSurface | None = None,
        t_hot: float | None = None,
        Dh: float | None = None,
        main_inlet: str | None = None,
        regime: CompressibilityLiteral = DEFAULT_REGIME,
    ):
        super().__init__(id)
        # Which closure (#588): the compressible stagnation closure (#357,
        # #576) or its low-Mach limit, Pt = P + (m/A)^2 / (2 rho) at
        # T_static = T0 -- the network's regime, so that an incompressible
        # network has no sonic limit hidden in its chambers.
        self.regime = regime
        # MERGE CHAMBER (#471). With main_inlet declared the chamber accepts
        # further inflows as SIDE STREAMS, still one outlet. Its axis is the
        # outlet direction; the main stream arrives along it, each side stream
        # brings m u_jet cos(theta) of axial momentum (the injecting element's
        # injection_momentum) and discharges at the chamber static pressure.
        # The node's (P, Pt) is the chamber/outlet state; the main-inlet
        # element sees the main-face state from chamber_merge_face_state.
        # Splits stay refused: one momentum equation cannot set two outlet
        # pressures.
        self.main_inlet = main_inlet
        self._auto_area = area is None
        self.area = area if area is not None else 0.0
        #: Where self.area came from once resolved (issue #262).
        self._area_source = "user" if area is not None else ""
        self.port_angles_deg: dict[str, float] = port_angles_deg or {}
        self.length = length
        self.Dh = Dh
        #: Where self.Dh came from; set here (not in resolve) so a second
        #: resolve pass does not relabel an inherited Dh as user-set.
        self._Dh_source = "user" if Dh is not None else ""
        self.surface = surface or ConvectiveSurface()
        self.t_hot = t_hot
        self.upstream_elements = []
        self.energy_boundaries: list[EnergyBoundary] = []
        #: Set by the solver on a converged root where this chamber is a
        #: choked-bounded junction (#788): the reported static state, else
        #: None. Reporting only; no residual reads it.
        self._choked_bounded: dict[str, float] | None = None

    @property
    def has_convective_surface(self) -> bool:
        return self.surface.area > 0.0

    def add_energy_boundary(self, eb: EnergyBoundary) -> None:
        """Attach an energy source/sink to this momentum chamber."""
        _attach_energy_boundary(self, eb)

    def set_port_angle(self, element_id: str, angle_deg: float) -> None:
        """Set the angle of the flow for a specific connected element."""
        self.port_angles_deg[element_id] = angle_deg

    def get_port_angle(self, element_id: str) -> float:
        """Get the angle of the flow for a specific connected element (defaults to 0.0)."""
        return self.port_angles_deg.get(element_id, 0.0)

    def unknowns(self) -> list[str]:
        # Pure Pressure-Flow: Temperature and Composition are derived forward, not unknowns.
        # A Tt_target EnergyBoundary adds its heat flow (#792).
        return [f"{self.id}.P", f"{self.id}.Pt"] + [
            eb.unknown_name for eb in target_energy_boundaries(self)
        ]

    def compute_derived_state(
        self, upstream_states: list[NetworkMixtureState]
    ) -> tuple[float, list[float], Any]:
        """Derived T and Y for a momentum chamber (simple mixing + energy boundaries)."""

        # Store total mass flow for use in residuals
        self._total_m_dot = sum(s.m_dot for s in upstream_states) if upstream_states else 1.0

        # Named m_dot Jacobian: {var_name: d(m_total)/d(var)}.
        self._upstream_element_ids, self._upstream_m_dot_jac = _summed_inflow_jac(upstream_states)

        # Merge chamber: keep each inflow's supply state for the main-face
        # offset, which needs the chamber's own state (main_face_offset).
        self._inflows = [
            (getattr(s, "_element_id", None), s, dict(getattr(s, "_m_dot_jac_names", {})))
            for s in upstream_states
        ]

        if not upstream_states:
            return 300.0, list(cb.mole_to_mass(cb.species.dry_air())), None

        streams = [cb.MassStream(s.m_dot, s.Tt, s.Pt, s.Y) for s in upstream_states]
        Q_total = _applied_heat(self)
        fraction_total = sum(eb.fraction for eb in self.energy_boundaries)

        mix_res = _solver_tools.mixer_from_streams_and_jacobians(
            streams, Q=Q_total, fraction=fraction_total
        )

        return mix_res.Tt_mix, mix_res.Y_mix, mix_res

    def _closure(self, state: NetworkMixtureState, m_dot: float | None = None) -> Any:
        """The chamber closure at ``state`` in this chamber's regime (#588)."""
        m = float(getattr(self, "_total_m_dot", 0.0)) if m_dot is None else m_dot
        return _solver_tools.momentum_chamber_residual_and_jacobian(
            state.P,
            state.Pt,
            m,
            state.Tt,
            state.Y,
            self.area,
            self.regime == "incompressible",
        )

    def htc_and_T(self, state: NetworkMixtureState):
        """Compute heat transfer coefficient for the momentum chamber."""
        if self.surface.area == 0.0:
            return None
        return self._surface_htc_and_T(state, self.surface)

    def _gas_side_htc_and_T(self, state: NetworkMixtureState):
        """The chamber's coefficient as another component's gas side.

        An effusion plate on this chamber owns the wall and its area; it reads
        only h (per unit area) and T_aw from the chamber's surface model, so
        the chamber needs no wetted area of its own. Its ``surface.area`` is 0
        unless given (#1022), and the surface's zero-area gate would then hide
        the coefficient; a unit-area copy of the surface evaluates it.
        """
        surface = self.surface
        if surface.area == 0.0:
            surface = copy.copy(surface)
            surface.area = 1.0
        return self._surface_htc_and_T(state, surface)

    def _surface_htc_and_T(self, state: NetworkMixtureState, surface: ConvectiveSurface):
        T_hot = self.t_hot if self.t_hot is not None else math.nan
        m_dot_total = getattr(self, "_total_m_dot", 0.0)
        # The correlation takes the BULK STATIC state; state.Tt is the
        # chamber's stagnation T, so use the static state its closure solved
        # (#578). T_aw then sits between T_static and Tt by the recovery
        # factor, as it should.
        T_static = float(state.Tt)
        c = None
        if self.area > 0:
            c = self._closure(state, m_dot_total)
            T_static = float(c.T_static)
            rho, _ = _safe_rho(float(c.rho_static))
            u = m_dot_total / (rho * self.area)
        else:
            u = 1.0

        diameter = self.Dh if self.Dh is not None else math.sqrt(4.0 * self.area / math.pi)
        length = self.length if self.length is not None else diameter

        r = surface.htc_and_T(
            T=T_static,
            P=state.P,
            X=state.X,
            velocity=u,
            diameter=diameter,
            length=length,
            T_hot=T_hot,
            flow_area=self.area,
        )
        # The partials are at the static T; T_static moves with the node's T0,
        # its P and its throughput.
        return None if r is None or c is None else _ChamberHtc(r, c)

    def residuals(
        self, state: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:

        # Momentum chamber: Pt = P_static + 0.5 * rho * v^2
        # Use total mass flow computed during compute_derived_state
        m_dot_total = getattr(self, "_total_m_dot", 0.0)

        # Use C++ function for residual and analytical Jacobian
        result = self._closure(state, m_dot_total)

        res = [result.residual]
        jac: dict[int, dict[str, float]] = {
            0: {
                f"{self.id}.P": result.d_res_dP,
                f"{self.id}.Pt": result.d_res_dP_total,
                # T is not an unknown at this node; the solver relays the
                # entry through the propagation chain (see the "." branch in
                # _residuals_and_jacobian). The compressible closure depends on
                # T through T0, a(T) and s(T), not just through rho, so the
                # term is no longer small enough to omit.
                f"{self.id}.Tt": result.d_res_dT,
            }
        }
        # Composition, relayed the same way: a node mixing gases of
        # different molar mass moves rho, a and s with every split (#517).
        for k, d in enumerate(result.d_res_dY):
            if d != 0.0:
                jac[0][f"{self.id}.Y[{k}]"] = d

        # Add Jacobian entries for upstream element mass flows using the
        # named Jacobian built in compute_derived_state. This correctly handles
        # multi-unknown elements (e.g. TeeJunctionElement with m_dot_com /
        # m_dot_branch) by mapping d(residual)/d(m_total) through the chain
        # d(m_total)/d(each_unknown) provided by flow_jac_at_node.
        for var_name, coeff in getattr(self, "_upstream_m_dot_jac", {}).items():
            if coeff != 0.0:
                jac[0][var_name] = jac[0].get(var_name, 0.0) + result.d_res_dmdot * coeff

        _energy_target_rows(self, state, res, jac)
        return res, jac

    def has_side_streams(self) -> bool:
        """A merge: a declared main inlet and at least one OTHER inflow. The
        main itself need not be flowing in -- a reversed main keeps its face
        (#493)."""
        if self.main_inlet is None:
            return False
        return any(eid != self.main_inlet for eid, _, _ in getattr(self, "_inflows", []))

    def main_face_state(
        self, state: NetworkMixtureState, m_main: float | None = None
    ) -> tuple[float, float, dict[str, float], dict[str, float]]:
        """(P_face, Pt_face, dP_face/d(name), dPt_face/d(name)) for the main inlet.

        Exact (compressible) impulse over the constant-area chamber, C++'s
        ``chamber_merge_face_state``: the face static pressure on the subsonic
        root, the face Pt from the same stagnation closure as the chamber.
        Derivatives are keyed by unknown name; ``"<node>.Tt"`` keys are relayed
        by the solver.

        ``m_main`` is the main element's signed flow into the chamber, for a
        main that is not among the inflows -- REVERSED (#493). The impulse
        holds for either sign (the x-momentum flux through the face is
        m_main^2 / (rho A) both ways), so the face is kept: the fluid at the
        face is then the chamber's, and the outlet carries the side streams
        less the main's outflow. It used to switch off, a 245 Pa step in the
        main's exit pressure as its flow crossed zero.
        """
        elems = {e.id: e for e in self.upstream_elements}
        main_state, main_names = None, {}
        S, dS = 0.0, {}
        for eid, st, names in getattr(self, "_inflows", []):
            if eid == self.main_inlet:
                main_state, main_names = st, names
                continue
            hook = getattr(elems.get(eid), "injection_momentum", None)
            if hook is None:
                continue
            J, dJ = hook(st, state)
            S += J
            for k, v in dJ.items():
                dS[k] = dS.get(k, 0.0) + v
        out_names = dict(getattr(self, "_upstream_m_dot_jac", {}))
        m_out = float(getattr(self, "_total_m_dot", 0.0))
        main_T_key = None
        if main_state is not None:
            m_main = float(main_state.m_dot)
            T_main = float(main_state.Tt)  # stagnation: the face solves its static T (#578)
            X_main = main_state.X
        else:
            # Reversed (or idle) main: chamber fluid leaves through the face.
            m_main = float(m_main) if m_main is not None else 0.0
            T_main = float(state.Tt)
            X_main = state.X
            main_names = {f"{self.main_inlet}.m_dot": 1.0}
            m_out += m_main
            out_names[f"{self.main_inlet}.m_dot"] = (
                out_names.get(f"{self.main_inlet}.m_dot", 0.0) + 1.0
            )
            main_T_key = f"{self.id}.Tt"
        r = cb.chamber_merge_face_state(
            m_main,
            T_main,
            X_main,
            m_out,
            state.P,
            state.Tt,
            state.X,
            S,
            self.area,
            self.regime == "incompressible",
        )
        self._face_choked = bool(r.choked)
        main_src = elems[self.main_inlet].from_node if self.main_inlet in elems else None

        # Composition (#743): the chamber gas's and the main stream's.
        main_Y_node = self.id if main_T_key is not None else main_src

        def chain(d_mm, d_mo, d_J, d_P, d_T, d_Tm, d_Y=(), d_Ym=()) -> dict[str, float]:
            out: dict[str, float] = {}

            def add(k: str, v: float) -> None:
                out[k] = out.get(k, 0.0) + v

            for k, v in main_names.items():
                add(k, d_mm * v)
            for k, v in out_names.items():
                add(k, d_mo * v)
            for k, v in dS.items():
                add(k, d_J * v)
            add(f"{self.id}.P", d_P)
            add(f"{self.id}.Tt", d_T)
            if main_T_key is not None:
                add(main_T_key, d_Tm)
            elif main_src is not None:
                add(f"{main_src}.Tt", d_Tm)
            for k, v in enumerate(d_Y):
                if v != 0.0:
                    add(f"{self.id}.Y[{k}]", v)
            if main_Y_node is not None:
                for k, v in enumerate(d_Ym):
                    if v != 0.0:
                        add(f"{main_Y_node}.Y[{k}]", v)
            return out

        return (
            float(r.P_face),
            float(r.Pt_face),
            chain(
                r.dPf_dm_main,
                r.dPf_dm_out,
                r.dPf_dJ,
                r.dPf_dP,
                r.dPf_dT,
                r.dPf_dT_main,
                r.dPf_dY,
                r.dPf_dY_main,
            ),
            chain(
                r.dPtf_dm_main,
                r.dPtf_dm_out,
                r.dPtf_dJ,
                r.dPtf_dP,
                r.dPtf_dT,
                r.dPtf_dT_main,
                r.dPtf_dY,
                r.dPtf_dY_main,
            ),
        )

    def mach(self, state: NetworkMixtureState) -> float:
        """Computes Mach number using internal total mass flow and area."""

        m_dot_total = getattr(self, "_total_m_dot", 0.0)
        if self.area <= 0 or m_dot_total <= 1e-12:
            return 0.0
        if getattr(self, "_choked_bounded", None) is not None:
            # A choked-bounded junction (#788) is reported as the plenum its
            # root is: the choked jet's excess lost in it, the gas at rest.
            return 0.0

        # The closure solves the static state from the chamber's stagnation
        # state.Tt, and its Mach is the one the residual used (#576).
        return float(self._closure(state, m_dot_total).M)

    def static_density(self, state: NetworkMixtureState) -> tuple[float, dict[str, float]]:
        """Static density of the chamber and d(rho)/d(name), for elements that
        read this chamber as a junction PORT (#576).

        The closure reads the chamber's stagnation ``state.Tt`` and solves
        the static T from (Tt, P, |m_dot| / A); rho evaluated at Tt would
        disagree with it. A junction's
        momentum rows must use the same static state, or the port's Pt and
        the junction's impulse describe two different flows (Wang 2014 MAE
        0.13 -> 0.19, rising with Mach, with only the chamber side fixed).

        Keys: ``{id}.P`` (an unknown), ``{id}.Tt`` / ``{id}.Y[k]`` (relayed by
        the solver) and the upstream flow unknowns that set this chamber's
        through-flow (the chain the chamber's own residual uses). With no area
        the chamber has no velocity: rho(T0, P) with the ideal-gas columns.
        """
        if not self.area > 0.0:
            rho = float(state.density())
            d0: dict[str, float] = {f"{self.id}.P": rho / state.P, f"{self.id}.Tt": -rho / state.Tt}
            W = float(cb.mwmix(state.X))
            for k, y in enumerate(state.Y):
                if y > 0.0:
                    d0[f"{self.id}.Y[{k}]"] = -rho * W / float(cb.species_molar_mass(k))
            return rho, d0
        m_dot_total = getattr(self, "_total_m_dot", 0.0)
        r = self._closure(state, m_dot_total)
        d: dict[str, float] = {f"{self.id}.P": r.drho_dP, f"{self.id}.Tt": r.drho_dT}
        for k, v in enumerate(r.drho_dY):
            if v != 0.0:
                d[f"{self.id}.Y[{k}]"] = v
        if r.drho_dmdot != 0.0:
            for var, coeff in getattr(self, "_upstream_m_dot_jac", {}).items():
                if coeff != 0.0:
                    d[var] = d.get(var, 0.0) + r.drho_dmdot * coeff
        return float(r.rho_static), d

    def diagnostics(self, state: NetworkMixtureState) -> dict[str, float | str]:
        m_dot_total = getattr(self, "_total_m_dot", 0.0)
        # Properties at the STATIC state the closure solved from the
        # chamber's stagnation state.Tt (#576).
        T_static = float(state.Tt)
        bounded = getattr(self, "_choked_bounded", None)
        if bounded is None and self.area > 0:
            T_static = float(self._closure(state, m_dot_total).T_static)
        # A choked-bounded junction (#788) is reported as the plenum its root
        # is (the solver reports P = Pt): the gas at rest, T = Tt.
        cs = cb.complete_state(T_static, state.P, state.X)

        rho = cs.thermo.rho
        u = (
            (abs(m_dot_total)) / ((rho) * (self.area))
            if rho > 0 and self.area > 0 and bounded is None
            else 0.0
        )
        diameter = self.Dh if self.Dh is not None else math.sqrt(4.0 * self.area / math.pi)
        re = (rho) * (u) * (diameter) / (cs.transport.mu) if u > 0 and cs.transport.mu > 0 else 0.0
        mach = u / cs.thermo.a if cs.thermo.a > 0 else 0.0
        Tt, Pt = (state.Tt, state.Pt)

        # Not evaluated (no convective surface, no flow, or a choked-bounded
        # junction): NaN, which the GUI shows as no value. A 0 read as zero
        # heat transfer, and T_aw = Tt as a real recovery temperature.
        nu_val = math.nan
        st_val = math.nan
        htc_val = math.nan
        t_aw_val = math.nan
        # A choked-bounded junction has no surface state of its area: its
        # surface correlation is not evaluated (and no wall may read it, or
        # the junction is refused, #788).
        if self.has_convective_surface and bounded is None:
            h_res = self.htc_and_T(state)
            if h_res is not None:
                nu_val = h_res.Nu
                htc_val = h_res.h
                t_aw_val = h_res.T_aw
                st_val = _stanton(h_res)
        extra: dict[str, float] = {}
        if bounded is not None:
            extra = {
                "choked_bounded": 1.0,
                "P_sonic": float(bounded["P_sonic"]),
                "Pt_sonic": float(bounded["Pt_sonic"]),
                "expansion_loss": float(bounded["Pt_sonic"]) - float(Pt),
            }

        return {
            **extra,
            "Tt": float(Tt),
            "Pt": float(Pt),
            "T_static": T_static,
            "h": cs.thermo.h,
            "s": cs.thermo.s,
            "u": cs.thermo.u,
            "rho": rho,
            "gamma": cs.thermo.gamma,
            "a": cs.thermo.a,
            "cp": cs.thermo.cp,
            "cv": cs.thermo.cv,
            "mw": cs.thermo.mw,
            "mu": cs.transport.mu,
            "k": cs.transport.k,
            "Pr": cs.transport.Pr,
            "nu": cs.transport.nu,
            "Re": re,
            "Dh": diameter,
            "velocity": u,
            "mach": mach,
            "Nu": nu_val,
            "St": st_val,
            "htc": htc_val,
            "T_aw": t_aw_val,
            "area": float(self.area or 0.0),
            "area_source": self._area_source,
        }

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        self.upstream_elements = graph.get_upstream_elements(self.id)
        if self.main_inlet is not None:
            # Merge chamber: one area, and the main inlet's port area IS it.
            # Inherit it from the main-inlet element; an explicit area that
            # disagrees is refused in validate(), never silently reconciled.
            main = graph.elements.get(self.main_inlet)
            main_area = getattr(main, "area", None) if main is not None else None
            self._main_inlet_area = main_area if main_area else None
            if self._auto_area and self._main_inlet_area:
                self.area = self._main_inlet_area
                self.Dh = getattr(main, "Dh", None) or 2.0 * math.sqrt(self.area / math.pi)
                self._area_source = f"inherited from {self.main_inlet}"
                self._Dh_source = self._area_source
                return
        # Inherit hydraulic diameter from upstream channel when not user-specified
        if self.Dh is None:
            for elem in self.upstream_elements:
                if isinstance(elem, ChannelElement) and elem.diameter is not None:
                    self.Dh = elem.diameter
                    self._Dh_source = f"channel '{elem.id}'"
                    break
        # Derive cross-section area from Dh when area was not user-specified
        if self._auto_area and self.Dh is not None:
            self.area = math.pi * (self.Dh / 2.0) ** 2
            self._area_source = "derived from Dh" if self._Dh_source == "user" else self._Dh_source
        elif self._auto_area:
            # No channel to inherit from: widen the search to any area-bearing
            # neighbour before falling back (issue #262).
            area, source = _infer_flow_area(graph, [self.id])
            if area is not None:
                self.area = area
                self.Dh = 2.0 * math.sqrt(area / math.pi)
                self._Dh_source = source
                self._area_source = source
            else:
                # Unlike AreaChangeElement, this area is a dynamic-head
                # reference rather than the geometry being modelled, and a
                # network can legitimately never read it -- so a nominal
                # default stands. It is recorded, not silent: anything that
                # consumes it reports the provenance in diagnostics.
                self.area = DEFAULT_CHAMBER_AREA
                self._area_source = "default"

    def _geometry_error(self) -> str:
        """Why this chamber's section cannot be built (#1024), or "": a given
        area or Dh not positive, or a Dh above 2 sqrt(A/pi)."""
        return _section_error(
            f"MomentumChamberNode {self.id!r}",
            self.area,
            not self._auto_area,
            self.Dh,
            self._Dh_source == "user",
        )

    def validate(self) -> None:
        if error := self._geometry_error():
            raise ValueError(error)
        if self.main_inlet is None:
            return
        ids = [e.id for e in self.upstream_elements]
        if self.main_inlet not in ids:
            raise ValueError(
                f"MomentumChamberNode {self.id!r}: main_inlet {self.main_inlet!r} "
                f"does not flow into it (inflows: {sorted(ids)})."
            )
        main_area = getattr(self, "_main_inlet_area", None)
        if not self._auto_area and main_area and abs(self.area / main_area - 1.0) > 1e-9:
            raise ValueError(
                f"MomentumChamberNode {self.id!r}: area {self.area:.6g} m^2 differs "
                f"from its main inlet {self.main_inlet!r} ({main_area:.6g} m^2). "
                "A merge chamber is a constant-area control volume; leave the "
                "area unset to inherit it."
            )

    def _area_is_a_statement(self, graph: "FlowNetwork") -> bool:
        """Whether the resolved area states the geometry (#1023, #1016's
        rule): a user's area or Dh, or one inherited from a channel or another
        area-bearing neighbour. Not the DEFAULT_CHAMBER_AREA placeholder of a
        chamber nothing sizes, nor a merge chamber's area inherited from an
        orifice main inlet (a bore, never an area a neighbour measures
        against)."""
        source = self._area_source
        if source in ("", "default") or not (self.area and self.area > 0.0):
            return False
        if self._auto_area and self.main_inlet is not None:
            main = graph.elements.get(self.main_inlet)
            if isinstance(main, OrificeElement) and source == f"inherited from {self.main_inlet}":
                return False
        return True

    def pipe_area_warnings(self, graph: "FlowNetwork") -> list[str]:
        """One message per adjoining pipe wider than this chamber (#1023).

        The chamber's closure carries the flow through its own area with no
        loss, so a pipe wider than it is a contraction into the chamber (or an
        expansion out of it) that costs nothing. Someone may mean that, so it
        is a warning, not a refusal; FlowNetwork.validate issues it and the
        GUI's resolve preview shows it. Compared on the equivalent diameter to
        NEIGHBOUR_AREA_RTOL, only where the area is a statement
        (_area_is_a_statement) and against a pipe whose diameter is not the
        GUI's default."""
        if not self._area_is_a_statement(graph):
            return []
        provenance = "" if self._area_source == "user" else f", {self._area_source}"
        out: list[str] = []
        seen: set[str] = set()
        for e in graph.get_upstream_elements(self.id) + graph.get_downstream_elements(self.id):
            if e.id in seen or not isinstance(e, ChannelElement):
                continue
            seen.add(e.id)
            if not (e.diameter and e.diameter > 0.0) or getattr(e, "_area_source", "") == "default":
                continue
            A_pipe = math.pi * (e.diameter / 2.0) ** 2
            if math.sqrt(self.area / A_pipe) >= 1.0 - NEIGHBOUR_AREA_RTOL:
                continue
            into = e.to_node == self.id
            out.append(
                f"MomentumChamberNode '{self.id}': area {self.area:.6g} m^2{provenance} is "
                f"smaller than the adjoining pipe '{e.id}' (D = {e.diameter:.6g} m, "
                f"A = {A_pipe:.6g} m^2): the flow "
                + (
                    "contracts from the pipe into the chamber"
                    if into
                    else "expands from the chamber into the pipe"
                )
                + " with no loss. If that area change is real, model it explicitly "
                "with an AreaChangeElement between them; otherwise make the chamber "
                "at least the pipe's area."
            )
        return out


class PressureBoundary(NetworkNode):
    """
    Supplies fixed stagnation pressure and temperature. Contributes no unknowns and no residuals.
    Essential as a reference pressure and absolute mass sink/source for the network.
    """

    is_boundary: bool = True
    is_flow_boundary: bool = True
    state_at_rest: bool = True

    def __init__(
        self,
        id: str,
        Pt: float = 101325.0,
        Tt: float = 300.0,
        Y: list[float] | None = None,
        coupling: str = "auto",
    ) -> None:
        super().__init__(id)
        self.Pt = Pt
        self.Tt = Tt
        self.Y = Y
        #: How the supplied pressure couples to a duct meeting this boundary.
        #:
        #: ``"total"``   -- it is the stagnation pressure at the connection.
        #:                  Correct for an INFLOW: a reservoir supplying the
        #:                  network, where Pt is what the reservoir holds.
        #: ``"static"``  -- it is the static pressure at the duct face, and the
        #:                  duct's exit dynamic head is dissipated in the
        #:                  expansion. The standard model for a bare duct
        #:                  discharging into a plenum or to atmosphere.
        #: ``"auto"``    -- infer from flow direction: inflow takes total,
        #:                  outflow takes static.
        #:
        #: ``auto`` is well defined wherever a boundary has one role. It cannot
        #: decide for a boundary whose flow REVERSES during a solve, or that
        #: serves both roles -- inflow wants total, outflow wants static, and
        #: there is no single right answer. That case is what the explicit
        #: setting is for; whoever builds such a network knows whether the exit
        #: is a plain opening or a diffuser, and the solver does not.
        #:
        #: Pinning stagnation pressure at an outflow caps the mass flux at the
        #: sonic value FOR THAT PRESSURE, a constraint the physical problem
        #: never imposed. Measured on the combustor of #351: 1.34x an
        #: impossible ceiling under total coupling, M = 0.752 and comfortable
        #: under static. See issue #360.
        self.coupling: str = coupling

    def exit_head_lost(self, is_outflow: bool) -> bool:
        """Whether a duct meeting this boundary loses its exit dynamic head."""
        if self.coupling == "static":
            return True
        if self.coupling == "total":
            return False
        return bool(is_outflow)

    def unknowns(self) -> list[str]:
        return []

    def residuals(
        self, state: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        return [], {}

    def compute_derived_state(
        self, upstream_states: list[NetworkMixtureState]
    ) -> tuple[float, list[float], Any]:
        # Boundary nodes define their own state

        Y = self.Y if self.Y is not None else list(cb.mole_to_mass(cb.species.dry_air()))
        return self.Tt, Y, None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass


class MassFlowBoundary(NetworkNode):
    """
    Supplies fixed mass flow and stagnation temperature. Its pressure floats to satisfy the flow equations.
    Cannot exclusively define a network, a pressure reference is also required.
    """

    is_boundary: bool = True
    is_flow_boundary: bool = True
    state_at_rest: bool = True

    def __init__(
        self,
        id: str,
        m_dot: float = 0.1,
        Tt: float = 300.0,
        Y: list[float] | None = None,
    ) -> None:
        super().__init__(id)
        self.m_dot = m_dot
        self.Tt = Tt
        self.Y = Y

    def unknowns(self) -> list[str]:
        # Pressure floats to whatever is required to push the defined mass flow
        return [f"{self.id}.P", f"{self.id}.Pt"]

    def residuals(
        self, state: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Stagnation assumptions: Pt = P_static
        res = [state.Pt - state.P]
        jac = {0: {f"{self.id}.P": -1.0, f"{self.id}.Pt": 1.0}}
        return res, jac

    def compute_derived_state(
        self, upstream_states: list[NetworkMixtureState]
    ) -> tuple[float, list[float], Any]:
        """
        Computes (Tt, Y, Jac) for a MassFlowBoundary.
        - If upstream_states exists, it acts as a SINK (Outlet): inherits state via mixing.
        - If upstream_states is empty, it acts as a SOURCE (Inlet): uses user settings.
        """

        if upstream_states:
            # SINK behavior: Automatically mix upstream streams using C++ core logic
            # This ensures energy conservation and correct species transport at outlets.
            streams = [cb.MassStream(s.m_dot, s.Tt, s.Pt, s.Y) for s in upstream_states]
            mix_res = _solver_tools.mixer_from_streams_and_jacobians(streams)
            return mix_res.Tt_mix, mix_res.Y_mix, mix_res

        # SOURCE behavior: Use developer/user-defined boundary constants
        Y = self.Y if self.Y is not None else list(cb.mole_to_mass(cb.species.dry_air()))
        return self.Tt, Y, None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass


FuelBranchLiteral = Literal["lean", "rich"]


class FuelBoundary(MassFlowBoundary):
    """A fuel supply whose mass flow puts a combustor's exit at a target total
    temperature (#892): a fuel controller.

    The fuel mass flow is a solver unknown (as ``{id}.ln_m_dot``, so it stays
    positive), closed by one
    temperature row, ``Tt(combustor) - Tt_target = 0`` [K]. Its column comes
    through the solver's relay of the combustor's derived Tt, the
    combustor's own d(Tt)/d(inflow) from C++, so the Jacobian is analytic.

    With reaction the exit Tt rises with the equivalence ratio to a peak near
    stoichiometric, then falls: a reachable target has two roots, and
    ``branch`` chooses one -- "lean" (the default) or "rich" (rich-quench-lean
    staging, reformers). The start is on that branch, and a solution found on
    the other side of stoichiometric is refused, never returned. A target
    below the unfuelled mixture's Tt or above the peak has no root and is
    refused with the reason.

    ``combustor`` names the combustor it controls, like a probe names what it
    reads: the fuel may reach it by any route -- an injector, a fuel line, a
    premix plenum the air also enters -- as long as all of it does. A route
    that lets part of the fuel reach another node without passing the
    combustor, and a second controller or a ``Tt_target`` heat boundary on the
    same combustor, are refused.
    """

    def __init__(
        self,
        id: str,
        combustor: str,
        Tt_target: float,
        branch: FuelBranchLiteral = "lean",
        Tt: float = 300.0,
        Y: list[float] | None = None,
        m_dot_guess: float | None = None,
    ) -> None:
        if branch not in ("lean", "rich"):
            raise ValueError(
                f"FuelBoundary '{id}': branch must be 'lean' or 'rich', not {branch!r}"
            )
        if not (math.isfinite(Tt_target) and Tt_target > 0.0):
            raise ValueError(f"FuelBoundary '{id}': Tt_target must be a positive temperature [K]")
        super().__init__(id, m_dot=0.0 if m_dot_guess is None else float(m_dot_guess), Tt=Tt, Y=Y)
        self.combustor = combustor
        self.Tt_target = float(Tt_target)
        self.branch: FuelBranchLiteral = branch
        self.m_dot_guess = m_dot_guess
        #: The combustor's current Tt, written by the solver before residuals.
        self._target_Tt: float | None = None
        #: Why the start could not be put on the requested branch, if it could not.
        self._guess_note = ""
        #: The oxidizer at the start state [(m_dot, Tt, Y)], for diagnosing a failed solve.
        self._start_oxidizer: list[tuple[float, float, list[float]]] = []
        #: ln(m_dot) the flow is held at instead of the target, while the
        #: solver's fuel-fixed stage solves the flows first (None: the target).
        self._frozen_ln_m: float | None = None

    @property
    def m_dot_unknown(self) -> str:
        """The unknown is ln(m_dot) [ln(kg/s)]: the fuel flow stays positive
        through every Newton step (the solver has no bounds), and a step is
        a ratio, which suits flows that span decades."""
        return f"{self.id}.ln_m_dot"

    def unknowns(self) -> list[str]:
        return [*super().unknowns(), self.m_dot_unknown]

    def residuals(
        self, state: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        res, jac = super().residuals(state)
        if self._frozen_ln_m is not None:
            # The fuel-fixed stage: the flow held, the row kept in kelvin.
            jac[len(res)] = {self.m_dot_unknown: self.Tt_target}
            res.append(
                (math.log(max(float(self.m_dot), 1e-300)) - self._frozen_ln_m) * self.Tt_target
            )
            return res, jac
        Tt = self._target_Tt if self._target_Tt is not None else self.Tt_target
        jac[len(res)] = {f"{self.combustor}.Tt": 1.0}
        res.append(float(Tt) - self.Tt_target)
        return res, jac

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Something to burn (#892): a composition with no combustible species
        # has no root at any target; say so before a solve wanders.
        Y = self.Y if self.Y is not None else list(cb.mole_to_mass(cb.species.dry_air()))
        X = list(cb.mass_to_mole(Y))
        if float(cb.oxygen_required_per_kg_mixture(X)) <= 0.0:
            names = [n for n, x in zip(cb.species.names, X, strict=True) if x > 1e-12]
            what = (
                "an oxidizer, not a fuel"
                if "O2" in names and len(names) == 1
                else "nothing combustible"
            )
            raise ValueError(
                f"FuelBoundary '{self.id}': its composition ({', '.join(names)}) is {what}: "
                "no fuel flow can raise the combustor's Tt. Give it a fuel (e.g. CH4, H2)."
            )
        node = graph.nodes.get(self.combustor)
        if not isinstance(node, CombustorNode):
            raise ValueError(
                f"FuelBoundary '{self.id}': combustor {self.combustor!r} is not a CombustorNode "
                "in this network"
            )
        # One temperature target per combustor: a second FuelBoundary would
        # leave the split between them free (a "converged" answer with an
        # arbitrary pilot/main ratio), a Tt_target heat boundary a second
        # equation for the same Tt.
        others = [
            n.id
            for n in graph.nodes.values()
            if isinstance(n, FuelBoundary) and n is not self and n.combustor == self.combustor
        ]
        if others:
            raise ValueError(
                f"FuelBoundary '{self.id}': combustor {self.combustor!r} already has a fuel "
                f"controller ({', '.join(others)}). One target Tt has one fuel flow: feed the "
                "other fuels as MassFlowBoundary (a set flow, or a set share of this one)."
            )
        if target_energy_boundaries(node):
            raise ValueError(
                f"FuelBoundary '{self.id}': combustor {self.combustor!r} also has a Tt_target "
                "heat boundary: both would set its Tt. Keep one."
            )
        # Any route, as long as ALL of the fuel reaches the combustor: then the
        # target pins the fuel flow, and the solver can tell the combustor's
        # oxidizer from its inflow less the fuel. A path that leaves for
        # anywhere else (a bypass, a second burner) is refused.
        stack, seen, reaches = [self.id], {self.id}, False
        while stack:
            nid = stack.pop()
            downstream = graph.get_downstream_elements(nid)
            if nid != self.id and not downstream:
                raise ValueError(
                    f"FuelBoundary '{self.id}': part of its fuel can reach {nid!r} without "
                    f"passing {self.combustor!r}; all of it must reach the combustor it controls."
                )
            for e in downstream:
                if e.to_node == self.combustor:
                    reaches = True
                    continue
                if e.to_node not in seen:
                    seen.add(e.to_node)
                    stack.append(e.to_node)
        if not reaches:
            raise ValueError(
                f"FuelBoundary '{self.id}': no flow path from it to {self.combustor!r}"
            )

    def _fuel_path_elements(self, graph: "FlowNetwork") -> set[str]:
        """The elements that carry fuel only: from this boundary through
        nodes no other flow enters, up to the combustor or the first node
        where the fuel meets other streams (a premix plenum's air)."""
        out: set[str] = set()
        stack = [self.id]
        while stack:
            nid = stack.pop()
            for e in graph.get_downstream_elements(nid):
                out.add(e.id)
                to = e.to_node
                if to != self.combustor and len(graph.get_upstream_elements(to)) == 1:
                    stack.append(to)
        return out


class WallNode(NetworkNode):
    """Closed-end boundary that imposes zero mass flow.

    The Pt = P residual enforces zero velocity at the dead end.
    The solver's automatic mass-balance equation (applied to every non-PressureBoundary
    node) forces m_dot = 0 on all connected elements without any additional logic here.
    """

    is_boundary: bool = True
    state_at_rest: bool = True

    def __init__(self, id: str) -> None:
        super().__init__(id)

    def unknowns(self) -> list[str]:
        return [f"{self.id}.P", f"{self.id}.Pt"]

    def residuals(
        self, state: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        res = [state.Pt - state.P]
        jac = {0: {f"{self.id}.P": -1.0, f"{self.id}.Pt": 1.0}}
        return res, jac

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass


def _zero_loss(ctx: object) -> tuple[float, float]:
    return 0.0, 0.0


class CombustorNode(NetworkNode):
    """
    Combustion chamber adding the heat of reaction to the flow it mixes.
    Evaluates adiabatic combustion using the C++ core backend based on
    the selected CombustionMethodLiteral.
    Automatically handles mixing when multiple upstream connections exist.

    ``method``:
      - ``"complete"`` (default): complete combustion to CO2 + H2O.
      - ``"equilibrium"``: PARTIAL equilibrium -- complete combustion, then
        reforming + water-gas shift among CO, H2, CO2, H2O and the fuel. No
        dissociation (no O, OH, H, NO), by design: identical to "complete"
        when lean, cooler when rich. A combustor rarely reaches full
        equilibrium in its residence time, least of all for NO (whose
        endothermic equilibrium formation would remove heat the real flame
        keeps). See
        ``CombustionMethod`` in include/state.h for the decision record.

    Pressure loss is **not** a property of the combustor node.  Attach a
    :class:`PressureLossElement` on an adjacent edge (upstream diffuser or
    downstream liner exit) to apply combustor pressure loss; theta-aware
    correlations (``LinearThetaFractionLoss`` / ``LinearThetaHeadLoss``)
    automatically pick up ``_T_unburned`` / ``_T_burned`` from this node.
    """

    state_at_rest: bool = True  # Pt = P by definition

    # Supplies theta = T_burned / T_unburned - 1 to adjacent PressureLossElement.
    has_theta: bool = True

    # to_dict serialises the given geometry, not what resolve_topology
    # inherited into the same attributes (#1099).
    _init_kwarg_attrs = {
        "area": lambda n: None if n._auto_area else n.area,
        "Dh": lambda n: n.Dh if n._Dh_source == "user" else None,
    }

    def __init__(
        self,
        id: str,
        method: CombustionMethodLiteral = "complete",
        area: float | None = None,
        Dh: float | None = None,
        surface: ConvectiveSurface | None = None,
        t_hot: float | None = None,
    ):
        super().__init__(id)
        # Only "equilibrium" was tested, so any other string (a typo) ran as
        # "complete" without a word (#1024).
        _check_combustion_method(id, method)
        self.method = method
        self._auto_area = area is None
        self.area = area if area is not None else 0.0
        #: Where self.area came from once resolved (issue #262).
        self._area_source = "user" if area is not None else ""
        self.Dh = Dh
        #: Where self.Dh came from; set here (not in resolve) so a second
        #: resolve pass does not relabel an inherited Dh as user-set.
        self._Dh_source = "user" if Dh is not None else ""
        self.surface = surface or ConvectiveSurface()
        self.t_hot = t_hot
        self.upstream_elements = []
        self.energy_boundaries: list[EnergyBoundary] = []
        # Populated by compute_derived_state so PressureLossElement can read theta.
        self._T_unburned: float = 300.0
        self._T_burned: float = 300.0

    @property
    def has_convective_surface(self) -> bool:
        return self.surface.area > 0.0

    def add_energy_boundary(self, eb: EnergyBoundary) -> None:
        """Attach an energy source/sink to this combustor (post-combustion)."""
        _attach_energy_boundary(self, eb)

    def unknowns(self) -> list[str]:
        # Pure Pressure-Flow: Temperature and Composition are derived forward, not unknowns.
        # A Tt_target EnergyBoundary adds its heat flow (#792).
        return [f"{self.id}.P", f"{self.id}.Pt"] + [
            eb.unknown_name for eb in target_energy_boundaries(self)
        ]

    def compute_derived_state(
        self, upstream_states: list[NetworkMixtureState]
    ) -> tuple[float, list[float], Any]:
        """Derived T and Y for a combustor (Reaction + Mixing + energy boundaries)."""

        # Store total mass flow and unburned temperature for use in diagnostics
        self._total_m_dot = sum(s.m_dot for s in upstream_states) if upstream_states else 1.0
        # T_unburned is the ENTHALPY-mixed reactant temperature (#718), the
        # same mix the C++ combustor context defines theta with; a mass-weighted
        # Tt read ~4% off in theta for air at 700 K with CH4 at 300 K.
        self._T_unburned = 300.0
        self._unburned_streams = []
        if upstream_states and self._total_m_dot > 0:
            mix_u = _solver_tools.mixer_from_streams_and_jacobians(
                [cb.MassStream(s.m_dot, s.Tt, s.Pt, s.Y) for s in upstream_states]
            )
            self._T_unburned = float(mix_u.Tt_mix)
            # The streams behind T_unburned, for its Jacobian in a theta-sourced
            # PressureLossElement: ({unknown: dm/d(unknown)}, source node,
            # dT_u/dm, dT_u/dTt, dT_u/dY). T_u is not a relayed node property
            # (#481 C1).
            self._unburned_streams = [
                (
                    dict(getattr(s, "_m_dot_jac_names", {})),
                    getattr(s, "_src_node", None),
                    float(d.d_mdot),
                    float(d.d_Tt),
                    list(d.d_Y),
                )
                for s, d in zip(upstream_states, mix_u.dTt_mix_d_stream, strict=True)
            ]

        # Upstream element ids, and d(throughput)/d(unknown) for the wall
        # relay of the node's own htc, which reads the throughput.
        self._upstream_element_ids, self._upstream_m_dot_jac = _summed_inflow_jac(upstream_states)

        if not upstream_states:
            # Default fallback
            return 300.0, list(cb.mole_to_mass(cb.species.dry_air())), None

        streams = [cb.MassStream(s.m_dot, s.Tt, s.Pt, s.Y) for s in upstream_states]
        if self._total_m_dot > 0.0:
            P_ref = sum(s.m_dot * s.Pt for s in upstream_states) / self._total_m_dot
        else:
            P_ref = upstream_states[0].Pt

        Q_total = _applied_heat(self)
        fraction_total = sum(eb.fraction for eb in self.energy_boundaries)

        # Combustor only: mixing + combustion, NO pressure loss (that's on the edge).
        try:
            mix_res = _solver_tools.combustor_residuals_and_jacobians(
                streams,
                P_ref,
                Q=Q_total,
                fraction=fraction_total,
                pressure_loss=_zero_loss,
                use_equilibrium=(self.method == "equilibrium"),
            )
        except Exception as exc:
            # Unphysical intermediate state during Newton iteration (e.g. extreme phi).
            # Re-raise so the solver penalty path can guide the step back to physics.
            raise RuntimeError(
                f"CombustorNode '{self.id}': combustion call failed "
                f"(T_unburned={self._T_unburned:.1f} K)"
            ) from exc
        self._last_mix_res = mix_res
        # Expose burned temperature for adjacent PressureLossElement (theta source).
        self._T_burned = float(mix_res.Tt_mix)
        return mix_res.Tt_mix, mix_res.Y_mix, mix_res

    def residuals(
        self, state: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Stagnation constraint: Pt = P (combustor is a large-area, low-velocity volume).
        # The dynamic pressure 0.5*rho*v^2 is negligible compared to P at combustor scales,
        # so P ~= Pt is an excellent approximation (same as PlenumNode).
        # Pressure loss is applied by an adjacent PressureLossElement, not here.
        res = [state.Pt - state.P]
        jac = {0: {f"{self.id}.P": -1.0, f"{self.id}.Pt": 1.0}}
        _energy_target_rows(self, state, res, jac)
        return res, jac

    def static_T(self, state: NetworkMixtureState) -> float:
        """The static T of the gas moving through ``area`` (#782 audit): the
        node holds Pt = P (a large, slow volume), so its own isentrope gives
        Tt; with an area, the duct closure at G = |m| / area does. Tt with no
        area or no flow. Reporting only."""
        m = abs(float(getattr(self, "_total_m_dot", 0.0)))
        if not (self.area > 0.0 and m > 0.0 and state.Tt > 0.0 and state.Pt > 0.0):
            return float(state.Tt)
        s = cb.duct_htc_static_state(float(state.Tt), float(state.Pt), m / self.area, state.X)
        return float(s.T)

    def diagnostics(self, state: NetworkMixtureState) -> dict[str, float | str]:
        if state.P <= 0 or state.Tt <= 0:
            return {}
        cs = cb.complete_state(state.T, state.P, state.X)

        m_dot_total = getattr(self, "_total_m_dot", 0.0)
        rho = cs.thermo.rho
        u = (abs(m_dot_total)) / ((rho) * (self.area)) if rho > 0 and self.area > 0 else 0.0
        diameter = self.Dh if self.Dh is not None else math.sqrt(4.0 * self.area / math.pi)
        re = (rho) * (u) * (diameter) / (cs.transport.mu) if u > 0 and cs.transport.mu > 0 else 0.0
        mach = u / cs.thermo.a if cs.thermo.a > 0 else 0.0
        Tt, Pt = (state.Tt, state.Pt)

        # Not evaluated (no convective surface, no flow, or a choked-bounded
        # junction): NaN, which the GUI shows as no value. A 0 read as zero
        # heat transfer, and T_aw = Tt as a real recovery temperature.
        nu_val = math.nan
        st_val = math.nan
        htc_val = math.nan
        t_aw_val = math.nan
        if self.has_convective_surface:
            h_res = self.htc_and_T(state)
            if h_res is not None:
                nu_val = h_res.Nu
                htc_val = h_res.h
                t_aw_val = h_res.T_aw
                st_val = _stanton(h_res)

        try:
            phi = float(cb.equivalence_ratio(state.X))
        except Exception as e:
            warnings.warn(
                f"CombustorNode {self.id!r}: equivalence ratio not evaluated ({e})",
                RuntimeWarning,
                stacklevel=2,
            )
            phi = 0.0

        T_u = getattr(self, "_T_unburned", state.Tt)
        theta = (state.Tt / T_u) - 1.0 if T_u > 0 else 0.0

        return {
            "Tt": float(Tt),
            "Pt": float(Pt),
            "m_dot": m_dot_total,
            "h": cs.thermo.h,
            "s": cs.thermo.s,
            "u": cs.thermo.u,
            "rho": rho,
            "gamma": cs.thermo.gamma,
            "a": cs.thermo.a,
            "cp": cs.thermo.cp,
            "cv": cs.thermo.cv,
            "mw": cs.thermo.mw,
            "mu": cs.transport.mu,
            "k": cs.transport.k,
            "Pr": cs.transport.Pr,
            "nu": cs.transport.nu,
            "Re": re,
            "Dh": diameter,
            "velocity": u,
            "mach": mach,
            "Nu": nu_val,
            "St": st_val,
            "htc": htc_val,
            "T_aw": t_aw_val,
            "phi": phi,
            "theta": theta,
            "area": float(self.area or 0.0),
            "area_source": self._area_source,
        }

    def htc_and_T(self, state: NetworkMixtureState):
        """Compute heat transfer coefficient for the combustor wall."""
        if self.surface.area == 0.0:
            return None

        T_hot = self.t_hot if self.t_hot is not None else math.nan
        diameter = self.Dh if self.Dh is not None else 0.01
        length = 0.1  # Combustors are modelled as a single zone
        if not self.area > 0:
            return self.surface.htc_and_T(
                T=state.Tt,
                P=state.P,
                X=state.X,
                velocity=1.0,
                diameter=diameter,
                length=length,
                T_hot=T_hot,
                flow_area=self.area,
            )
        # The htc reads the node's Tt and Pt (its closure is Pt = P, so this
        # Pt is the node's static P); the gas sweeping the wall at m / (rho A)
        # is at the static state that flux implies, as in a plenum-fed duct.
        # Handing the stagnation T to
        # the correlation as static put T_aw above Tt: +4.7 K at Mach 0.15,
        # +46 K at Mach 0.49 for a 1430 K stream.
        s = copy.copy(state)
        s.m_dot = float(getattr(self, "_total_m_dot", 0.0))
        return _duct_htc(self.surface, s, self.area, diameter, length, T_hot)

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        self.upstream_elements = graph.get_upstream_elements(self.id)
        # Inherit hydraulic diameter from upstream channel when not user-specified
        if self.Dh is None:
            for elem in self.upstream_elements:
                if isinstance(elem, ChannelElement) and elem.diameter is not None:
                    self.Dh = elem.diameter
                    self._Dh_source = f"channel '{elem.id}'"
                    break
        # Derive cross-section area from Dh when area was not user-specified
        if self._auto_area and self.Dh is not None:
            self.area = math.pi * (self.Dh / 2.0) ** 2
            self._area_source = "derived from Dh" if self._Dh_source == "user" else self._Dh_source
        elif self._auto_area:
            # No channel to inherit from: widen the search to any area-bearing
            # neighbour before falling back (issue #262).
            area, source = _infer_flow_area(graph, [self.id])
            if area is not None:
                self.area = area
                self.Dh = 2.0 * math.sqrt(area / math.pi)
                self._Dh_source = source
                self._area_source = source
            else:
                # Unlike AreaChangeElement, this area is a dynamic-head
                # reference rather than the geometry being modelled, and a
                # network can legitimately never read it -- so a nominal
                # default stands. It is recorded, not silent: anything that
                # consumes it reports the provenance in diagnostics.
                self.area = DEFAULT_CHAMBER_AREA
                self._area_source = "default"

    def _geometry_error(self) -> str:
        """Why this combustor's section cannot be built (#1024), or "": a
        given area or Dh not positive, or a Dh above 2 sqrt(A/pi)."""
        return _section_error(
            f"CombustorNode {self.id!r}",
            self.area,
            not self._auto_area,
            self.Dh,
            self._Dh_source == "user",
        )

    def validate(self) -> None:
        _check_combustion_method(self.id, self.method)
        if error := self._geometry_error():
            raise ValueError(error)


def _check_combustion_method(node_id: str, method: object) -> None:
    """Refuse a combustion method outside CombustionMethodLiteral (#1024)."""
    allowed = get_args(CombustionMethodLiteral)
    if method not in allowed:
        raise ValueError(
            f"CombustorNode {node_id!r}: method {method!r} is not one of {list(allowed)}."
        )


# The hole Cd -- its selector tables, options (CD_INCLINATION_MODELS,
# CD_EXIT_CROSSFLOW_MODELS, _CD_CROSSFLOW_LAWS), every term and its
# diagnostics -- lives in hole_cd.py (_HoleCdChain).

#: Why Stolz was removed in 0.9.0 (#876), for the refusal of a 'Stolz' Cd.
STOLZ_REMOVAL_REASON = (
    "the Stolz Cd (ISO 5167:1980) was removed in 0.9.0 (#876): use "
    "'ReaderHarrisGallagher'. ISO replaced Stolz with Reader-Harris/Gallagher in "
    "ISO 5167-1:1991/Amd 1:1998 (and ISO 5167-2:2003), which has a sourced range "
    "of use and fits the EEC/API database better over Stolz's own range "
    "(Reader-Harris & Sattary 1996, Table 10: s 0.245 % vs 0.390 %)."
)


# Below this Cd an orifice kernel is evaluated at Cd = 1 and scaled: the
# kernels are linear in Cd and refuse A_eff <= 0 (OrificeElement._kernel_from).
_KERNEL_CD_MIN = 1e-12

#: OrificeElement.compressible_model (#854): the isentropic nozzle, Jobson
#: (1955) with the approach-velocity extension, or ISO 5167-2:2003's Eq. (5)
#: expansibility. Without an explicit choice the element selects one at
#: build (#861): see _JOBSON_DEFAULT_CORRELATIONS.
ORIFICE_COMPRESSIBLE_MODELS = ("nozzle", "jobson1955", "iso5167")

#: The Cd correlations whose Cd is referenced to INCOMPRESSIBLE flow, so it
#: carries no expansion of its own (#861 audit): ISO 5167-2's C (expansion is
#: its separate eps), Idelchik's diagram 4-18 (M <= 0.3, compressibility a
#: separate factor), Lichtarowicz (liquid rigs) and McGreehan-Schotsch's
#: baseline chain (its Y_o removes the expansion). For these the default
#: compressible model is jobson1955 on a sharp-edged bore at beta <=
#: JOBSON_1955_BETA_MAX, chosen once at build and never per iterate. Every
#: other case keeps the nozzle: 'fixed' (the user's measurement basis is
#: unknown) and McGreehan-Schotsch's
#: crossflow term (Rohde-fitted, nozzle-referenced), which only the duct-fed
#: EffusionPlateElement feeds; a plain orifice is plenum-fed, U1/Vi = 0.
_JOBSON_DEFAULT_CORRELATIONS = frozenset(
    {
        "ReaderHarrisGallagher",
        "IdelchikThick",
        "IdelchikBeveled",
        "IdelchikRounded",
        "Lichtarowicz",
        "McGreehanSchotsch",
    }
)

# Jobson's flow is not linear in Cd (Ci enters the force defect as 1/Ci^2):
# below this Cd it is evaluated here and scaled linearly to zero, the same
# trace-Cd treatment as _KERNEL_CD_MIN with a floor the closure stays
# well-conditioned at.
_JOBSON_CD_MIN = 1e-6

# OrificeElement's inner flow solves (#861): a stop at the iteration cap is a
# root when |a - m_calc(a)| is within these [-, kg/s]; otherwise the root is
# solved on a bracket (_root_on_bracket).
_OWN_FLOW_RTOL = 1e-9
_OWN_FLOW_ATOL = 1e-15


def _root_on_bracket(
    h: Callable[[float], float],
    a_end: float,
    lo: float | None,
    hi: float | None,
    probes: Callable[[], list[float]],
) -> float:
    """A root of h(a) = a - m_calc(a) after an inner solve stopped at its
    iteration cap at a_end (#861).

    a_end is returned unchanged when it already is a root to _OWN_FLOW_RTOL /
    _OWN_FLOW_ATOL, so a solve that met its root keeps its bits. Otherwise a
    bracket is formed from the signs seen on the way (lo: h < 0, hi: h > 0),
    completed from ``probes()`` where one side is missing, and the root is
    solved on it by brentq. With no bracket to be had, a_end is returned.
    """
    h_end = h(a_end)
    if abs(h_end) <= _OWN_FLOW_RTOL * max(abs(a_end), abs(a_end - h_end)) + _OWN_FLOW_ATOL:
        return a_end
    if lo is None or hi is None:
        for p in probes():
            hp = h(p)
            if hp == 0.0:
                return p
            if hp < 0.0 and lo is None:
                lo = p
            elif hp > 0.0 and hi is None:
                hi = p
            if lo is not None and hi is not None:
                break
        else:
            return a_end
    assert lo is not None and hi is not None
    return float(_brentq(h, lo, hi, xtol=1e-13 * max(abs(lo), abs(hi), 1e-30)))


# The taps the network's metering Cd is evaluated at. ISO 5167-2's
# permanent-loss relation (permanent_loss_ratio_and_dC) is Bernoulli to the
# vena contracta plus a momentum balance after it, so its tapping
# differential is the undisturbed-to-vena-contracta one, which D-D/2 taps
# measure. With corner taps' C the same relation books up to 6% more loss at
# beta 0.75 (tmp/probe_tap_invariance.py, #645). The tap is measurement, so
# the element's own `tap` only picks the differential it reports.
_NETWORK_TAP = cb.OrificeTap.D_D2

_TAPS = {
    "Corner": cb.OrificeTap.Corner,
    "D_D2": cb.OrificeTap.D_D2,
    "Flange": cb.OrificeTap.Flange,
}

# Partials key carrying d(m_calc)/d(Cd) through _kernel_from and
# _flow_and_partials. Between pipes m_calc is no longer linear in Cd: the
# permanent-loss ratio k(beta, C) moves the tap pressure too (#645).
_DM_DCD = "__dm_dCd__"


def _is_closed_restriction(elem: object) -> bool:
    """True for a restriction of zero area (OrificeElement.closed, #1024): a
    wall in the flow path, through which every root carries no flow."""
    return isinstance(elem, OrificeElement) and elem.closed


class OrificeElement(NetworkElement):
    """
    Orifice flow element with incompressible or compressible formulation.

    - regime='incompressible': m_dot = Cd * A * sqrt(2 * rho * dP)
    - regime='compressible': Uses isentropic nozzle flow, choked below P* (C1 through it)

    The discharge coefficient Cd is determined by the 'correlation' parameter:
      - 'fixed': Uses the user-supplied Cd value directly; Cd is required.
      - 'ReaderHarrisGallagher': Sharp thin-plate (ISO 5167-2:2003 / RHG).
      - 'IdelchikThick': Deep hole in a wall, Idelchik diagram 4-18a
        (requires plate_thickness). Valid Re 25 to 1e6.
      - 'IdelchikBeveled': Beveled-edge hole, diagram 4-18b
        (requires bevel_depth).
      - 'IdelchikRounded': Rounded-edge hole, diagram 4-18c
        (requires edge_radius).
      - 'McGreehanSchotsch': Long cooling hole (1988). Plenum-fed here, so
        its inlet-crossflow term sees U1/Vi = 0; only EffusionPlateElement's
        duct-fed panels feed it a crossflow.
      - 'Lichtarowicz': LONG orifice, Lichtarowicz, Duggins and Markland
        (1965), l/d 2-10 and Re 10 to 2e4 (requires plate_thickness).
        Refused below l/d = 1.5, where the source reports hysteresis.

    The first two are NORMED metering correlations: Cd is referenced to the
    tapping differential and is a function of beta = d/D. Both are evaluated
    at D-D/2 taps, the ones ISO's permanent-loss relation is derived for: the
    tap choice is a measurement convention and must not move the solved flow
    (#645), so ``tap`` only selects the differential reported as ``dP_tap``.

    Between pipes (a pipe on the side the flow comes from sets beta, and a
    pipe on the side it goes to recovers the pressure) a metering or 'fixed'
    orifice charges ISO 5167-2's PERMANENT loss, not the tapping
    differential: the row's drop is the stagnation difference between the
    two nodes, mapped to the tap pressure through k(beta, C)
    (orifice_tap_pressure_and_jacobian). Discharging into a node with no
    pipe beyond it, the jet's head is lost there and the drop is the tapping
    differential. The discharge-hole family keeps that form in a pipe too,
    until its in-pipe diagrams land. The last five are
    DISCHARGE correlations for a hole in a wall, where zeta is referenced to
    the hole velocity and there is no beta. They are not interchangeable, and
    there is deliberately no 'Auto' arm choosing between them from geometry.

    ``Cd`` is the 'fixed' correlation's value and nothing else (#872, 0.8.0):
    with any other correlation the Cd is the correlation's, and passing
    ``Cd`` too raises ValueError rather than being silently ignored, as it
    was before 0.8.0. 'fixed' without a ``Cd`` raises as well.

    ``hole_cd=HoleCd(...)`` (hole_cd.py) is the hole-Cd family's spec, the
    form every cooling and impingement element takes. Here it is an
    alternative to ``correlation``/``Cd``, which stay this element's own
    keywords (the metering use, not deprecated); a None field of the spec
    takes them (None ``correlation``: 'ReaderHarrisGallagher'), and a spec
    and a keyword that disagree raise. A plain orifice has no inclination,
    exit- or supply-crossflow term: those fields must stay "auto"/"none".

    The bore against the pipe it sits in (#788). The pipe diameter comes from
    the single channel on the side the flow comes from (``in_pipe=True``, the
    default); per direction, beta = d / D:
      - d < D: an orifice in a pipe, as above.
      - d == D (to ORIFICE_BORE_PIPE_RTOL, on either side): NO orifice. Zero
        loss, the flow passes unchanged: the row is Pt_in - Pt_out = 0, like
        LosslessConnectionElement. How an orifice is switched off without
        editing the network.
      - d > D: refused; there is no such geometry.
      - ``in_pipe=False``: no pipe diameter, whatever channels adjoin it.
        beta = 0 both ways, no pressure recovery: a hole in a wall, i.e. a
        lumped Cd*A restriction. Write a lumped restriction this way; a bore
        at or above the pipe is no longer read as one.

    The bore against a momentum chamber on its flow path (#1016), the same
    rule with areas (_bore_against_chambers): the bore's own area, never Cd
    A, larger than the chamber it opens into or draws from is refused,
    equal is no orifice, smaller is unchanged. Only an area that states the
    geometry is compared (_chamber_area_statement): never the default
    placeholder. The plate subclasses check their total hole area.

    ``compressible_model`` (#854) picks how the compressible regime expands
    the flow; the incompressible regime ignores it (each reduces to Cd E A
    sqrt(2 rho dP) there). An explicit choice always wins. Left at None it
    is selected ONCE at build from the correlation and the geometry (#861,
    0.8.0), never switched per iterate: 'jobson1955' for a Cd referenced to
    incompressible flow (ReaderHarrisGallagher, the Idelchik family,
    Lichtarowicz, McGreehanSchotsch) on a sharp-edged bore (edge_radius = 0,
    bevel_depth = 0) at beta <= JOBSON_1955_BETA_MAX in both flow
    directions; 'nozzle' otherwise ('fixed', a rounded or beveled
    edge, a wider bore). A Ci above Jobson's 0.7 is flagged by
    model_in_range, not refused. Before resolve_topology the attribute holds
    the selection for a plenum-fed bore (beta = 0).
      - 'nozzle': the isentropic nozzle at A_eff = Cd E A, choked at P*.
        Its expansion is a nozzle's, 3-9% below a sharp plate's at 10-20%
        pressure drop, and flat past choke.
      - 'jobson1955': Jobson (1955), the contraction of a sharp-edged jet
        widening as the pressure ratio falls, through choke and beyond, with
        the approach-velocity extension of #854 in a pipe. Its only input is
        the incompressible contraction Ci from the correlation's Cd, Cd = Cv
        Ci (``Cv``, the velocity coefficient, default 1; a Cv other than 1
        needs compressible_model='jobson1955' set explicitly). Physics,
        nothing fitted; its stated range (0.5 <= Ci <= 0.7, beta <= 0.75)
        is flagged in diagnostics (model_in_range), not enforced.
      - 'iso5167': ISO 5167-2:2003 Eq. (1) with Eq. (5)'s empirical
        expansibility, for a Reader-Harris/Gallagher metering plate only.
        Valid for p2/p1 >= 0.75 (flagged below it, where Eq. (5) is
        extrapolated and has no choke).
    """

    # The jet: its vena-contracta Mach from the pressure ratio, not the bore's.
    validity_mach_keys = ("mach_throat",)

    # Whether an unset compressible_model is selected from the correlation
    # and geometry (#861). The cooling and connection subclasses keep the
    # nozzle: their Cd laws are outside the #861 audit.
    _selects_compressible_model = True

    # The Cd a 'fixed' correlation takes when none is given (#872). None:
    # 'fixed' requires a Cd. A subclass with a documented default (the
    # impingement plate's Florschuetz 0.79) sets it.
    _fixed_Cd_default: float | None = None
    # The Cd the resistance start guess assumes when the Cd is a
    # correlation's (#872). A start guess only: the solved flow does not
    # depend on it.
    _seed_Cd_default: float = _ORIFICE_SEED_CD

    # to_dict serialises the given bore and compressible model, not the ones
    # resolve_topology defaulted or selected into the same attributes (#1099).
    _init_kwarg_attrs = {
        "diameter": lambda e: e.diameter if e._area_source == "user" else None,
        "area": lambda e: e.area if e._area_source == "user" else None,
        "compressible_model": lambda e: (
            e.compressible_model if e._compressible_model_explicit else None
        ),
    }

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        Cd: float | None = None,
        diameter: float | None = None,
        regime: Literal["incompressible", "compressible"] = DEFAULT_REGIME,
        correlation: str | None = None,
        plate_thickness: float = 0.0,
        edge_radius: float = 0.0,
        bevel_depth: float = 0.0,
        area: float | None = None,
        injection_angle_deg: float = 90.0,
        tap: str = "Flange",
        in_pipe: bool = True,
        compressible_model: str | None = None,
        Cv: float = 1.0,
        hole_cd: HoleCd | None = None,
    ):
        super().__init__(id, from_node, to_node)
        # The Cd choice: this element's own ``correlation``/``Cd`` keywords
        # (the metering interface, not deprecated here; None:
        # 'ReaderHarrisGallagher'), or a HoleCd (hole_cd.py), whose None
        # fields take those keywords. A spec and a keyword that disagree
        # raise. The cooling and impingement subclasses resolve their own
        # HoleCd and pass its correlation and Cd here.
        if hole_cd is not None:
            if not isinstance(hole_cd, HoleCd):
                raise ValueError(
                    f"{type(self).__name__} {id!r}: hole_cd is a HoleCd, not "
                    f"{type(hole_cd).__name__}"
                )
            for what, kw, val in (
                ("correlation", correlation, hole_cd.correlation),
                ("Cd", Cd, hole_cd.Cd),
            ):
                if kw is not None and val is not None and kw != val:
                    raise ValueError(
                        f"{type(self).__name__} {id!r}: {what}={kw!r} and "
                        f"hole_cd=HoleCd({what}={val!r}) disagree; give it once."
                    )
            if hole_cd.correlation is not None:
                correlation = hole_cd.correlation
            if hole_cd.Cd is not None:
                Cd = hole_cd.Cd
        if correlation is None:
            correlation = "ReaderHarrisGallagher"
        # False: no pipe diameter, a lumped Cd*A restriction (#788).
        self.in_pipe = bool(in_pipe)
        # How the compressible regime expands the flow (#854); ignored when
        # incompressible, where every model reduces to m = Cd E A
        # sqrt(2 rho dP). None: selected at build (#861), see
        # _default_compressible_model.
        if compressible_model is not None and compressible_model not in ORIFICE_COMPRESSIBLE_MODELS:
            raise ValueError(
                f"OrificeElement {id!r}: compressible_model {compressible_model!r} is not "
                f"one of {list(ORIFICE_COMPRESSIBLE_MODELS)}."
            )
        if compressible_model == "iso5167" and correlation != "ReaderHarrisGallagher":
            raise ValueError(
                f"OrificeElement {id!r}: compressible_model 'iso5167' is ISO 5167-2:2003's "
                "Eq. (5) expansibility, defined with its own discharge coefficient: use "
                f"correlation='ReaderHarrisGallagher', not {correlation!r}."
            )
        if not (0.0 < Cv <= 1.0):
            raise ValueError(f"OrificeElement {id!r}: Cv = {Cv!r} must be in (0, 1].")
        if Cv != 1.0 and compressible_model != "jobson1955":
            raise ValueError(
                f"OrificeElement {id!r}: Cv is Jobson's velocity coefficient and applies "
                "to compressible_model='jobson1955' only: set it explicitly."
            )
        # An explicit model always wins; None is resolved below and again in
        # resolve_topology once beta is known (#861).
        self._compressible_model_explicit = compressible_model is not None
        self.compressible_model = compressible_model or "nozzle"
        self.Cv = float(Cv)
        # bore == pipe (#788) or bore area == an adjacent momentum chamber's
        # (#1016): no orifice, set in resolve_topology, with the reason.
        self._pass_through = False
        self._pass_through_reason = ""
        if tap not in _TAPS:
            raise ValueError(f"OrificeElement {id!r}: tap {tap!r} is not one of {sorted(_TAPS)}.")
        # The taps whose differential diagnostics() reports as dP_tap. Never
        # the solved flow: that is evaluated at _NETWORK_TAP (#645).
        self.tap = tap
        # Angle of the discharged jet to the axis of a merge chamber it feeds
        # as a side stream (#471). 90 = normal injection, no axial momentum.
        self.injection_angle_deg = float(injection_angle_deg)
        # The hole's Cd: every term, its options and diagnostics (hole_cd.py).
        # Its options are all off here; the cooling-hole and impingement
        # subclasses configure theirs.
        self._hole_cd = _HoleCdChain(self)

        # A bore of zero is a CLOSED restriction (#1024); a negative or
        # non-finite one is no geometry.
        for what, size in (("diameter", diameter), ("area", area)):
            if size is not None and not (math.isfinite(size) and size >= 0.0):
                raise ValueError(
                    f"{type(self).__name__} {id!r}: {what} = {size!r} must be finite and "
                    "non-negative (0 closes the restriction: no flow)."
                )
        if diameter is not None:
            self.diameter: float | None = diameter
            self.area: float | None = math.pi * (diameter / 2.0) ** 2
        elif area is not None:
            self.area = area
            self.diameter = math.sqrt(4.0 * area / math.pi)
        else:
            # Both None: fall back to DEFAULT_ORIFICE_BORE in resolve_topology
            # (a bore cannot be inferred from a neighbour).
            self.diameter = None
            self.area = None
        # Provenance of the bore: "user" when a diameter or area was given, else
        # set to "default" in resolve_topology (issue #559).
        self._area_source = "user" if (diameter is not None or area is not None) else ""

        if correlation == "Miller":
            raise ValueError(
                f"OrificeElement {id!r}: the 'Miller' correlation was removed. "
                "No published source matched its formula (and the network's "
                "'Miller' had been computing Stolz). Use 'ReaderHarrisGallagher' "
                "(Stolz was removed in 0.9.0, #876)."
            )
        if correlation == "Stolz":
            raise ValueError(f"OrificeElement {id!r}: {STOLZ_REMOVAL_REASON}")
        name = type(self).__name__
        if correlation == "fixed":
            if Cd is None:
                Cd = type(self)._fixed_Cd_default
            if Cd is None:
                raise ValueError(
                    f"{name} {id!r}: correlation='fixed' needs a Cd. Pass Cd=..., or "
                    "pick a correlation that computes it."
                )
            if not math.isfinite(Cd):
                raise ValueError(f"{name} {id!r}: Cd={Cd!r} must be a finite number.")
        elif Cd is not None:
            raise ValueError(
                f"{name} {id!r}: Cd={Cd!r} with correlation={correlation!r}. The "
                f"correlation computes its own Cd and a given one would be ignored. "
                "For a fixed discharge coefficient pass correlation='fixed'; to use "
                "the correlation, drop Cd."
            )
        # The 'fixed' correlation's Cd; None when a correlation computes it.
        self.Cd: float | None = Cd
        # The resolved spec: correlation and Cd as used, the terms as given
        # (a subclass replaces it with its own). Read by FlowNetwork.to_dict.
        try:
            self.hole_cd = dataclasses.replace(
                hole_cd if hole_cd is not None else HoleCd(), correlation=correlation, Cd=Cd
            )
        except ValueError as exc:
            raise ValueError(f"{name} {id!r}: {exc}") from exc
        self._hole_cd.configure_plain(f"{type(self).__name__} {id!r}", self.hole_cd)
        self.regime = regime
        self.correlation = correlation
        self.use_correlation = correlation != "fixed"
        self.plate_thickness = plate_thickness
        self.edge_radius = edge_radius
        # Bevel DEPTH along the hole axis, Idelchik diagram 4-18b's l/Dh at a
        # bevel angle of 40-60 deg. Not the wall thickness.
        self.bevel_depth = bevel_depth
        self.upstream_diameter: float | None = None
        self._upstream_source = ""
        self.downstream_diameter: float | None = None
        # Velocity-of-approach ratio; recomputed in resolve_topology. Initialized
        # here so diagnostics() is safe before resolve runs.
        self.beta = 0.0
        # Pressure recovery into a downstream pipe, per direction (#645); set
        # in resolve_topology. Subclasses that skip it never recover.
        self._recovery_fwd = False
        self._recovery_rev = False
        # A discharge-hole Cd carried into the approach pipe, per direction
        # (#645): C_pipe = Cd_wall R. 1 without a pipe; set in resolve_topology.
        self._pipe_cd_ratio_fwd = 1.0
        self._pipe_cd_ratio_rev = 1.0
        # OrificeGeometry built in resolve_topology
        self._orifice_geom: object | None = None
        # The build-time selection for a plenum-fed bore (beta = 0);
        # resolve_topology re-selects with the pipes' beta.
        if not self._compressible_model_explicit:
            self.compressible_model = self._default_compressible_model(0.0)

    def _default_compressible_model(self, beta: float) -> str:
        """The compressible model an unset ``compressible_model`` resolves to
        (#861): 'jobson1955' for an incompressible-referenced Cd correlation
        (_JOBSON_DEFAULT_CORRELATIONS) on a sharp-edged bore at beta <=
        JOBSON_1955_BETA_MAX, else 'nozzle'. A function of construction data
        and beta only, so it is fixed for a built network."""
        sharp = self.edge_radius == 0.0 and self.bevel_depth == 0.0
        if (
            type(self)._selects_compressible_model
            and self.correlation in _JOBSON_DEFAULT_CORRELATIONS
            and sharp
            and beta <= cb.JOBSON_1955_BETA_MAX
        ):
            return "jobson1955"
        return "nozzle"

    @property
    def _seed_Cd(self) -> float:
        """The Cd a start guess assumes: the 'fixed' Cd, else
        ``_seed_Cd_default`` standing in for the correlation's (#872)."""
        return self.Cd if self.Cd is not None else type(self)._seed_Cd_default

    def _seed_resistance(self, rho: float, P_ref: float, m_ref: float) -> float | None:
        """Start guess: m = Cd A sqrt(2 rho dP); a correlation's Cd is
        approximated by ``_seed_Cd_default`` (#872). No orifice (bore ==
        pipe): no loss."""
        if self._pass_through:
            return 0.0
        if not self.area:
            return None
        return 1.0 / (2.0 * rho * (self._seed_Cd * self.area) ** 2)

    def _seed_flow_applies(self) -> bool:
        """A compressible orifice (#971): its kernel chokes, which the
        incompressible K does not."""
        return self.regime == "compressible" and not self._pass_through and bool(self.area)

    def _seed_flow(
        self, Pt_up: float, Pt_dn: float, Tt: float, X: list[float], Y: list[float]
    ) -> float | None:
        """The flow the residual's own compressible kernel (``_kernel_from``)
        passes plenum to plenum at the seed's Cd: no approach velocity, no
        recovering pipe."""
        model = getattr(self, "compressible_model", "nozzle")
        cd = self._seed_Cd
        try:
            if model == "jobson1955":
                r = _solver_tools.orifice_jobson_residuals_and_jacobian(
                    0.0,
                    Pt_up,
                    Pt_up,
                    Tt,
                    Y,
                    Pt_dn,
                    max(cd, _JOBSON_CD_MIN),
                    self.area,
                    0.0,
                    self.Cv,
                )
            elif model == "iso5167":
                r = _solver_tools.orifice_iso5167_residuals_and_jacobian(
                    0.0, Pt_up, Pt_up, Tt, Y, Pt_dn, cd, self.area, 0.0
                )
            else:
                r = _solver_tools.orifice_compressible_residuals_and_jacobian(
                    0.0, Pt_up, Pt_up, Tt, Y, Pt_dn, cd, self.area, 0.0
                )
        except (ValueError, RuntimeError):
            return None
        return abs(float(r.m_dot_calc))

    @property
    def closed(self) -> bool:
        """A bore of zero area: the restriction is CLOSED (#1024), a wall.

        Its row is m = 0, exact and linear (constant Jacobian 1 in its own
        flow); no kernel is evaluated, and the solved network carries no flow
        through it in either direction. A part of the network that only
        closed elements join to the rest keeps its own pressure boundaries;
        one they cut off from every PressureBoundary is refused at
        FlowNetwork.validate(), since nothing sets its pressure."""
        return self.area is not None and self.area == 0.0

    def _resolve_closed(self) -> None:
        """Topology of a closed restriction (#1024): no pipe, no beta, no
        recovery, no Cd geometry. Nothing flows, so nothing is discovered."""
        self.beta = 0.0
        self.beta_reverse = 0.0
        self._orifice_geom = None
        self._pass_through = False
        self._recovery_fwd = False
        self._recovery_rev = False
        self._pipe_cd_ratio_fwd = 1.0
        self._pipe_cd_ratio_rev = 1.0
        self._reverse_pipe_diameter = None
        self._reverse_pipe_source = ""

    @property
    def residual_scale_kind(self) -> str:
        """RESIDUAL-ROW-SCALING: m - m_calc, kg/s (every orifice-family
        closure); a total-pressure equality, Pa, where the bore equals the
        pipe and there is no orifice (#788). Must change with the residual."""
        return "p" if self._pass_through else "mdot"

    def _chamber_area_statement(
        self, graph: "FlowNetwork", node_id: str
    ) -> "tuple[float, str] | None":
        """The area of the momentum chamber at ``node_id`` when it is a
        statement about the geometry, with its provenance; None otherwise
        (#1016).

        A user's area, a user's Dh, or one inherited from a channel or another
        area-bearing neighbour is a statement. The DEFAULT_CHAMBER_AREA
        placeholder is not: a chamber nothing sizes says nothing about the
        bore. Nor is an area inherited from an orifice (a merge chamber's
        main inlet): a bore is never the area a neighbour measures against,
        and against itself it would read as equal. And a chamber that took
        its area from a channel is that pipe's end: the bore against the pipe
        is #788's rule, whose waiver is in_pipe=False (a lumped Cd*A
        restriction, no pipe diameter), so such an area is not checked then.
        """
        node = graph.nodes.get(node_id)
        if not isinstance(node, MomentumChamberNode):
            return None
        source = node._area_source
        if source in ("", "default") or not (node.area and node.area > 0.0):
            return None
        if not self.in_pipe and source.startswith("channel '"):
            return None
        if node._auto_area and node.main_inlet is not None:
            main = graph.elements.get(node.main_inlet)
            if isinstance(main, OrificeElement) and source == f"inherited from {node.main_inlet}":
                return None
        return float(node.area), source

    def _bore_against_chambers(self, graph: "FlowNetwork") -> str:
        """The geometric bore against a momentum chamber on either side of the
        flow path (#1016), #788's rule with areas: the bore's own area (a
        plate's total hole area, never Cd A) larger than the chamber's is
        refused, equal (to ORIFICE_BORE_PIPE_RTOL on the equivalent diameter,
        sqrt(A_bore / A_chamber)) is no orifice. Returns the
        pass-through's reason, or "" when the bore is smaller (or no chamber's
        area is a statement, see _chamber_area_statement)."""
        bore = float(self.area or 0.0)
        if not bore > 0.0:
            return ""
        what = "bore" if type(self) is OrificeElement else "total hole area"
        reason = ""
        for node_id in dict.fromkeys((self.from_node, self.to_node)):
            stated = self._chamber_area_statement(graph, node_id)
            if stated is None:
                continue
            A, source = stated
            # Compared as equivalent diameters, sqrt(A), so the band is
            # #788's exactly where the chamber is a pipe's end.
            d_ratio = math.sqrt(bore / A)
            if d_ratio > 1.0 + ORIFICE_BORE_PIPE_RTOL:
                side = "opens into" if node_id == self.to_node else "draws from"
                provenance = "" if source == "user" else f", {source}"
                raise ValueError(
                    f"{type(self).__name__} {self.id!r}: {what} {bore:.6g} m^2 is larger "
                    f"than the momentum chamber it {side} ({node_id!r}, A = {A:.6g} m^2"
                    f"{provenance}). There is no such geometry: the jet would be wider "
                    "than the chamber. A bore area equal to the chamber's is no "
                    "orifice (zero loss); otherwise make the chamber at least the "
                    "bore's area or the bore smaller."
                )
            if abs(d_ratio - 1.0) <= ORIFICE_BORE_PIPE_RTOL:
                reason = f"{what} = chamber '{node_id}' area: no orifice"
        return reason

    def _resolve_plate_against_chambers(self, graph: "FlowNetwork") -> None:
        """The plates' (effusion, film, slot, impingement) check against an
        adjacent momentum chamber (#1016): their resolve_topology skips the
        pipe discovery, so this is the whole of #788's rule for them."""
        self._pass_through = False
        self._pass_through_reason = ""
        reason = self._bore_against_chambers(graph)
        if reason:
            self._set_pass_through(reason)

    def _set_pass_through(self, reason: str) -> None:
        """No orifice (#788, #1016): beta = 1, no recovery, no pipe Cd."""
        self._pass_through = True
        self._pass_through_reason = reason
        self.beta = 1.0
        self._recovery_fwd = False
        self._recovery_rev = False
        self._pipe_cd_ratio_fwd = 1.0
        self._pipe_cd_ratio_rev = 1.0

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        self._to_node_obj = graph.nodes.get(self.to_node)
        self._to_side_closed = self._side_is_closed(graph, self.to_node)
        if self.closed:
            self._resolve_closed()
            return
        # Evaluate both sides of the node for geometry discovery
        upstream_elements = graph.get_upstream_elements(self.from_node)
        downstream_elements = graph.get_downstream_elements(self.to_node)

        # Ensure we don't pick ourselves if it's a tight chain
        upstream_channels = [
            e for e in upstream_elements if isinstance(e, ChannelElement) and e.id != self.id
        ]
        downstream_channels = [
            e for e in downstream_elements if isinstance(e, ChannelElement) and e.id != self.id
        ]

        if len(upstream_channels) == 1:
            self.upstream_diameter = upstream_channels[0].diameter
            # Provenance of the pipe diameter that sets beta (#559).
            self._upstream_source = f"channel '{upstream_channels[0].id}'"

        if len(downstream_channels) == 1:
            self.downstream_diameter = downstream_channels[0].diameter

        # Bore must be explicitly set by the user; fall back to the default when
        # unspecified (an orifice bore cannot be inferred from a neighbour -- it
        # is a restriction, not the channel area). upstream_diameter is kept only
        # for the velocity-of-approach Beta correction.
        if self.diameter is None:
            self.diameter = DEFAULT_ORIFICE_BORE
            self.area = math.pi * (self.diameter / 2.0) ** 2
            self._area_source = "default"

        # Pre-compute beta for velocity-of-approach factor

        d_bore = math.sqrt(4.0 * self.area / math.pi)
        self.beta = 0.0
        self._orifice_geom = None
        self._pass_through = False
        self._pass_through_reason = ""

        # Reversed flow approaches from the to-node's side (#554): the single
        # pipe there, declared either way, sets the velocity of approach.
        self.beta_reverse = 0.0
        to_pipes = [
            e
            for e in graph.get_upstream_elements(self.to_node)
            + graph.get_downstream_elements(self.to_node)
            if isinstance(e, ChannelElement) and e.id != self.id
        ]
        from_pipes = [
            e
            for e in graph.get_upstream_elements(self.from_node)
            + graph.get_downstream_elements(self.from_node)
            if isinstance(e, ChannelElement) and e.id != self.id
        ]
        # The reverse approach pipe's diameter: its static state is the
        # orifice's p1 in reverse (#651).
        self._reverse_pipe_diameter = None
        self._reverse_pipe_source = ""

        if not self.in_pipe:
            # No pipe diameter (#788): a lumped Cd*A restriction, a hole in a
            # wall. Whatever channels adjoin it set no beta and recover
            # nothing, in either direction.
            self.upstream_diameter = None
            self._upstream_source = ""
            to_pipes, from_pipes = [], []

        # The pipe that sets beta, per direction: the single upstream channel
        # forward, the single to-side channel in reverse. The bore against it
        # (#788): equal is no orifice, larger is refused. There is no silent
        # fallback to a wall hole any more; that is in_pipe=False, said.
        pipes = []
        if self.upstream_diameter and self.upstream_diameter > 0:
            pipes.append((self.upstream_diameter, self._upstream_source))
        if len(to_pipes) == 1 and to_pipes[0].diameter and to_pipes[0].diameter > 0:
            pipes.append((to_pipes[0].diameter, f"channel '{to_pipes[0].id}'"))
        pass_through = ""
        for D, source in pipes:
            if d_bore > D * (1.0 + ORIFICE_BORE_PIPE_RTOL):
                raise ValueError(
                    f"OrificeElement {self.id!r}: bore {d_bore:.6g} m is larger than "
                    f"the pipe it sits in ({source}, D = {D:.6g} m). There is no such "
                    "geometry. A bore equal to the pipe is no orifice (zero loss); for "
                    "a lumped Cd*A restriction set in_pipe=False (no pipe diameter, "
                    "beta = 0)."
                )
            if abs(d_bore - D) <= ORIFICE_BORE_PIPE_RTOL * D:
                pass_through = "bore = pipe D: no orifice"
        # The bore's area against a momentum chamber on its flow path (#1016),
        # in_pipe or not: a hole in a wall cannot be wider than the chamber
        # behind it either.
        pass_through = self._bore_against_chambers(graph) or pass_through
        if pass_through:
            # beta = 1: no restriction, and no orifice physics to set up.
            self._set_pass_through(pass_through)
            return

        if len(to_pipes) == 1 and to_pipes[0].diameter:
            self.beta_reverse = d_bore / to_pipes[0].diameter
            self._reverse_pipe_diameter = to_pipes[0].diameter
            self._reverse_pipe_source = f"channel '{to_pipes[0].id}'"

        if self.upstream_diameter and self.upstream_diameter > 0:
            self.beta = d_bore / self.upstream_diameter

        # The compressible model, selected once here from the correlation and
        # geometry when the user set none (#861): Jobson's range must hold for
        # the beta of either flow direction, since one model serves both.
        if not self._compressible_model_explicit:
            self.compressible_model = self._default_compressible_model(
                max(self.beta, self.beta_reverse)
            )

        # Pressure recovery (#645), per flow direction: the approach pipe sets
        # beta, and a single wider pipe on the far side recovers the jet's
        # pressure, so the row charges ISO's permanent loss. A plain node
        # beyond (no pipe) takes the jet's head as lost. Every correlation:
        # a discharge-hole Cd enters through its in-pipe C below.
        def _recovers(pipes: list[ChannelElement]) -> bool:
            return len(pipes) == 1 and bool(pipes[0].diameter) and d_bore < pipes[0].diameter

        self._recovery_fwd = self.beta > 0.0 and _recovers(to_pipes)
        self._recovery_rev = self.beta_reverse > 0.0 and _recovers(from_pipes)

        if self.use_correlation:
            # Build geometry descriptor for Cd correlation.
            # D=0 when no upstream channel known (plenum connection) -> beta=0 -> RHG extrapolates to ~0.597.
            D_up = (
                self.upstream_diameter
                if self.upstream_diameter and self.upstream_diameter > 0
                else 1.0
            )
            if D_up <= d_bore:
                D_up = d_bore * 10.0
            geom = cb.OrificeGeometry()
            geom.d = d_bore
            geom.D = D_up
            geom.t = self.plate_thickness
            geom.r = self.edge_radius
            geom.tap = _NETWORK_TAP
            self._orifice_geom = geom

        # A discharge-hole Cd is plenum-measured (Idelchik's wall, M-S,
        # Lichtarowicz): in a pipe it becomes the vena-contracta C that gives
        # Idelchik's in-pipe loss (4-13..4-16) through the recovery row,
        # C_pipe = Cd_wall R(beta), R -> 1 as the pipe grows (#645). It used
        # to take E on the plenum Cd and lose the whole jet head: 16-20% low
        # flow at beta 0.5 between pipes.
        self._pipe_cd_ratio_fwd = self._hole_cd.pipe_cd_ratio(self.beta)
        self._pipe_cd_ratio_rev = self._hole_cd.pipe_cd_ratio(self.beta_reverse)

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def _hole_reynolds(self, state_in: "NetworkMixtureState") -> float:
        """Reynolds number of the flow through ONE hole, based on its diameter.

        The discharge-hole correlations (Idelchik, McGreehan-Schotsch) are
        functions of this, not of the pipe Reynolds number. For a plain
        orifice there is one hole carrying all the flow; EffusionPlateElement
        overrides this to divide by the hole count.
        """
        d = self._orifice_geom.d if self._orifice_geom is not None else 0.0
        if d <= 0.0:
            return 1.0e5
        mu = cb.transport_state(state_in.Tt, state_in.P, state_in.X).mu
        if mu <= 0.0:
            mu = 1.8e-5
        m_dot_hole = abs(state_in.m_dot) / self._hole_count()
        if m_dot_hole <= 1e-12:
            # No flow yet: a mid-range seed, not the top of the table, so the
            # first Newton step starts on a responsive part of the curve.
            return 1.0e4
        return (4.0 * m_dot_hole) / (math.pi * d * mu)

    def _hole_count(self) -> float:
        """Holes carrying the element's mass flow in parallel. One, here."""
        return 1.0

    def _effective_Cd(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None = None,
    ) -> float:
        """The element's Cd, every term composed (_HoleCdChain.cd)."""
        return self._hole_cd.cd(state_in, state_out, flows)

    def injection_momentum(
        self, state_in: "NetworkMixtureState", state_chamber: "NetworkMixtureState"
    ) -> tuple[float, dict[str, float]]:
        """Streamwise momentum this jet brings into a merge chamber (#471).

        J = m w cos(theta), with w from C++'s ``jet_impulse``: the isentropic
        velocity from the supply's (Pt, Tt) to the chamber static pressure,
        or, choked, the sonic momentum plus its pressure thrust. Cd sets the
        jet's area, not its velocity, so it does not enter. Returns (J,
        dJ/d(name)) keyed by unknown name. Normal injection (90 deg) brings
        none.
        """
        cos_t = math.cos(math.radians(self.injection_angle_deg))
        if abs(cos_t) < 1e-12:
            return 0.0, {}
        jet = cb.jet_impulse(
            float(state_in.m_dot), state_in.Pt, state_in.Tt, state_chamber.P, state_in.X
        )
        return jet.J * cos_t, {
            f"{self.id}.m_dot": jet.dJ_dm * cos_t,
            f"{self.from_node}.Pt": jet.dJ_dPt * cos_t,
            f"{self.from_node}.Tt": jet.dJ_dTt * cos_t,
            f"{self.to_node}.P": jet.dJ_dP * cos_t,
        }

    def validate(self) -> None:
        self._hole_cd.validate()

    def _dCd_dnames(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None = None,
        continue_re_floor: bool = False,
    ) -> dict[str, float]:
        """d(Cd)/d(unknown) by name, every term (_HoleCdChain.dcd). The
        own-flow column through the discharge chamber rides under
        hole_cd._DCD_DM_VIA_CHAMBER."""
        return self._hole_cd.dcd(state_in, state_out, flows, continue_re_floor)

    def residuals(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None = None,
    ) -> tuple[list[float], dict[int, dict[str, float]]]:

        if self.closed:
            # Closed (#1024): a wall, m = 0. Exact and linear, so the Jacobian
            # is the constant 1 in its own flow; no kernel sees A = 0.
            return [float(state_in.m_dot)], {0: {f"{self.id}.m_dot": 1.0}}

        if self._pass_through:
            # Bore == pipe (#788): no orifice. Zero loss, the flow passes
            # unchanged; its rate is the network's. Exact and linear, so
            # the Jacobian is the constant +-1 of a lossless connection.
            return [state_in.Pt - state_out.Pt], {
                0: {f"{self.from_node}.Pt": 1.0, f"{self.to_node}.Pt": -1.0}
            }

        m_dot = state_in.m_dot
        effective_cd = self._effective_Cd(state_in, state_out, flows)

        # A discharge-hole Cd rises with the hole Reynolds number, which is
        # linear in the flow. Evaluated at the solver's iterate, the row
        # m - Cd(Re(m)) K is non-monotone in m -- slope +1 below the table's
        # Re floor, NEGATIVE while Cd climbs (IdelchikThick: 0.37 at Re 25 to
        # 0.86 at Re 1000), then rising to the root -- and |F| has a local
        # minimum at the floor that traps Newton (#582). So the element solves
        # its own consistency, m* = Cd(Re(m*)) K, and the row is m - m*(x):
        # linear in m, same root, same physics.
        # With a velocity of approach the drive p1 = p(Pt, G = |m| / A_pipe)
        # reads the flow too (#651), and at the iterate's m the forward kernel
        # ran BACKWARDS wherever 0 < Pt - p2 < q_pipe: the row's slope in Pt
        # flipped and Newton walked away from the root. Solved out the same
        # way, m* = m_calc(m*, x) has the sign of the stagnation drive.
        if self._approach_reads_flow():
            return self._flow_consistent_residuals(state_in, state_out, flows, effective_cd)
        if self._cd_follows_own_flow():
            if self._flow_nonlinear_in_cd():
                # _self_consistent_residuals takes m = Cd K (#582); Jobson's
                # flow is not linear in Cd, so solve m* = m_calc(m*) on the
                # kernel's own d(m)/d(Cd) instead (#854).
                return self._flow_consistent_residuals(state_in, state_out, flows, effective_cd)
            return self._self_consistent_residuals(state_in, state_out, flows, effective_cd)

        m_calc, partials = self._flow_and_partials(state_in, state_out, effective_cd)
        dm_dcd = partials.pop(_DM_DCD, 0.0)
        res = [m_dot - m_calc]

        # Assemble Jacobian with respect to all node and element unknowns
        jac: dict[int, dict[str, float]] = {0: {f"{self.id}.m_dot": 1.0}}
        for name, d in partials.items():
            jac[0][name] = jac[0].get(name, 0.0) - d

        # The Cd moves with the flow (Re) and, duct-fed, with the crossflow
        # (U1/Vi): d(m_calc)/d(x) += d(m_calc)/d(Cd) dCd/dx. That is m_calc/Cd
        # where m_calc is linear in Cd; between pipes k(beta, C) adds to it.
        if effective_cd > 0.0:
            for name, dcd in self._dCd_dnames(state_in, state_out, flows).items():
                name = _own_flow_key(name, self.id)
                jac[0][name] = jac[0].get(name, 0.0) - dm_dcd * dcd

        return res, jac

    def _flow_nonlinear_in_cd(self) -> bool:
        """True when m_calc is not linear in Cd even without a recovering
        pipe: Jobson's compressible closure (#854)."""
        return (
            self.regime == "compressible"
            and getattr(self, "compressible_model", "nozzle") == "jobson1955"
        )

    def _cd_follows_own_flow(self) -> bool:
        """True when Cd moves with this element's own flow through a hole
        Reynolds number with an analytic derivative
        (_HoleCdChain.follows_own_flow)."""
        return self._hole_cd.follows_own_flow()

    def _side_is_closed(self, graph: "FlowNetwork", start: str) -> bool:
        """True when the network beyond ``start`` (not crossing back through
        this element) holds no flow boundary: a dead end, closed by walls.

        Every root then carries zero flow through this element, so it never
        has to read that side as a source (#554). Reading it anyway made a
        dead-end orifice depend on a stagnant node's temperature, which runs
        like Q/m to the clamp as the flow vanishes, and Newton stalled.
        """
        seen = {start}
        todo = [start]
        while todo:
            nid = todo.pop()
            node = graph.nodes.get(nid)
            if node is None or getattr(node, "is_flow_boundary", False):
                return False
            for elem in graph.get_upstream_elements(nid) + graph.get_downstream_elements(nid):
                # A closed restriction is a wall (#1024): nothing flows past it.
                if elem is self or _is_closed_restriction(elem):
                    continue
                for other in dict.fromkeys([*elem.all_source_nodes(), *elem.all_sink_nodes()]):
                    if other not in seen:
                        seen.add(other)
                        todo.append(other)
        return True

    def _kernel_from(
        self,
        src: "NetworkMixtureState",
        dst: "NetworkMixtureState",
        src_id: str,
        dst_id: str,
        cd: float,
    ) -> tuple[float, dict[str, float]]:
        """C++ (m_calc, d m_calc / d name) for flow from ``src`` into ``dst``.

        The kernels take the upstream node's stagnation state and the
        downstream static pressure. The incompressible density reference is
        the upstream static by default; the solver's compressible-seed proxy
        switches it to the downstream static (_incompressible_p_ref =
        "outlet"), which tracks the compressible solution far better on
        blow-down networks (6-7x closer junction pressures,
        tmp/outlet_ref_seed_experiment.py).

        Into a recovering pipe (#645) the kernels expand to the TAP pressure
        that the nodes' stagnation difference implies through ISO's
        permanent-loss ratio, not to the downstream node's static. m_calc is
        then no longer linear in Cd; d(m_calc)/d(Cd) rides along under
        _DM_DCD.
        """
        # The velocity of approach is the upstream pipe's along the flow.
        forward = src_id == self.from_node
        beta = getattr(self, "beta", 0.0) if forward else getattr(self, "beta_reverse", 0.0)
        recovers = getattr(self, "_recovery_fwd" if forward else "_recovery_rev", False)
        compressible = self.regime == "compressible"
        # A discharge-hole Cd in its approach pipe (#645): a constant factor,
        # so d(m)/d(Cd_wall) = R d(m)/d(C_pipe).
        ratio = getattr(self, "_pipe_cd_ratio_fwd" if forward else "_pipe_cd_ratio_rev", 1.0)
        cd = cd * ratio
        # Both kernels are exactly linear in Cd (A_eff = Cd E A). A crossflow
        # Cd underflows to exactly 0 where the hole's pressure drop vanishes
        # (U1/Vi -> inf in McGreehan-Schotsch's exp(-Rv^n)), and the kernels
        # refuse A_eff <= 0: evaluate at Cd = 1 there and scale, which is
        # exact, so a trace Cd gives a trace flow instead of an exception.
        model = getattr(self, "compressible_model", "nozzle") if compressible else "nozzle"
        if model == "jobson1955":
            # Not linear in Cd: a trace Cd is scaled linearly from the floor.
            cd_k = cd if cd >= _JOBSON_CD_MIN else _JOBSON_CD_MIN
        else:
            cd_k = cd if cd >= _KERNEL_CD_MIN else 1.0
        scale = cd / cd_k
        # With a velocity of approach the kernels drive from the upstream
        # STATIC pressure (approach_from_static), since E carries the approach
        # head. That static state is the approach PIPE's (#651): from the
        # node's stagnation state at the pipe's mass flux, so a PlenumNode
        # (P = Pt) and a MomentumChamberNode give the same flow.
        p1, dp1 = self._approach_static(src, src_id, forward, beta)
        p_down = dst.P
        tap = None
        if recovers:
            tap = _solver_tools.orifice_tap_pressure_and_jacobian(
                p1, src.Pt, dst.Pt, beta, max(cd, _KERNEL_CD_MIN), compressible
            )
            p_down = float(tap.p_tap)
        if model == "jobson1955":
            r = _solver_tools.orifice_jobson_residuals_and_jacobian(
                src.m_dot, src.Pt, p1, src.Tt, src.Y, p_down, cd_k, self.area, beta, self.Cv
            )
        elif model == "iso5167":
            r = _solver_tools.orifice_iso5167_residuals_and_jacobian(
                src.m_dot, src.Pt, p1, src.Tt, src.Y, p_down, cd_k, self.area, beta
            )
        elif compressible:
            r = _solver_tools.orifice_compressible_residuals_and_jacobian(
                src.m_dot, src.Pt, p1, src.Tt, src.Y, p_down, cd_k, self.area, beta
            )
        else:
            _p_ref_outlet = getattr(self, "_incompressible_p_ref", "inlet") == "outlet"
            r = _solver_tools.orifice_residuals_and_jacobian(
                src.m_dot,
                src.Pt,
                p1,
                src.Tt,
                src.Y,
                p_down,
                cd_k,
                self.area,
                beta=beta,
                # Density at the downstream static: its sensitivity moves there.
                P_rho_ref=dst.P if _p_ref_outlet else -1.0,
            )
        g = float(r.d_mdot_dP_static_down)
        d: dict[str, float] = {
            f"{src_id}.Pt": float(r.d_mdot_dP_total_up),
            f"{src_id}.Tt": float(r.d_mdot_dTt_up),
            f"{dst_id}.P": float(r.d_mdot_dP_rho_ref) + (0.0 if tap else g),
        }
        # d(m_calc)/d(p1): the kernel's own column and, between pipes, the tap
        # pressure's through g. Spread over p1's dependencies below.
        dm_dp1 = float(r.d_mdot_dP_static_up)
        # The nozzle and ISO kernels are linear in Cd; Jobson's is not, and
        # returns its own column (exact from the floor up, the linear scale
        # below it).
        if model == "jobson1955" and cd >= _JOBSON_CD_MIN:
            dm_dcd = float(r.d_mdot_dCd)
        else:
            dm_dcd = float(r.m_dot_calc) / cd_k
        if tap is not None:
            dm_dp1 += g * float(tap.d_dP_static_up)
            # p_tap(p1, Pt_up, Pt_down, C): its other columns through g.
            for name, dp in (
                (f"{src_id}.Pt", tap.d_dP_total_up),
                (f"{dst_id}.Pt", tap.d_dP_total_down),
            ):
                if dp != 0.0:
                    d[name] = d.get(name, 0.0) + g * float(dp)
            dm_dcd += scale * g * float(tap.d_dC)
        # Outlet-referenced and plenum-fed, the flow does not read the source
        # static at all: leave the column out rather than carry a zero.
        if dm_dp1 != 0.0:
            for name, dp in dp1.items():
                d[name] = d.get(name, 0.0) + dm_dp1 * dp
        if tap is not None and d[f"{dst_id}.P"] == 0.0:
            # The tap pressure reads the far node's Pt, not its static.
            del d[f"{dst_id}.P"]
        for i, val in enumerate(r.d_mdot_dY_up):
            if val != 0.0:
                key = f"{src_id}.Y[{i}]"
                d[key] = d.get(key, 0.0) + float(val)
        if scale != 1.0:
            d = {k: v * scale for k, v in d.items()}
        d[_DM_DCD] = dm_dcd * ratio
        return float(r.m_dot_calc) * scale, d

    def _approach_static(
        self,
        src: "NetworkMixtureState",
        src_id: str,
        forward: bool,
        beta: float,
    ) -> tuple[float, dict[str, float]]:
        """The upstream static pressure the kernels drive from, and its
        derivatives by name (#651).

        With E (0 < beta < 1) it is the approach pipe's: C++'s
        approach_static_pressure from the node's (Pt, Tt) at G = |m| / A_pipe,
        so the node's own static pressure no longer matters. A PlenumNode
        reports P = Pt by definition, which gave the orifice p1 and rho1 a
        dynamic head too high (1.5% flow at beta 0.7, pipe Mach 0.18).
        Without E the node's static is returned (the kernels then read Pt).
        """
        D = self.upstream_diameter if forward else getattr(self, "_reverse_pipe_diameter", None)
        if not (0.0 < beta < 1.0) or not D:
            return float(src.P), {f"{src_id}.P": 1.0}
        A = math.pi * D * D / 4.0
        m = float(src.m_dot)
        a = _solver_tools.approach_static_pressure(
            float(src.Pt), float(src.Tt), abs(m) / A, src.Y, self.regime != "compressible"
        )
        sgn = 1.0 if m >= 0.0 else -1.0
        out = {
            f"{src_id}.Pt": float(a.dP_dPt),
            f"{src_id}.Tt": float(a.dP_dTt),
            f"{self.id}.m_dot": float(a.dP_dG) * sgn / A,
        }
        for i, val in enumerate(a.dP_dY):
            if val != 0.0:
                out[f"{src_id}.Y[{i}]"] = float(val)
        return float(a.P), out

    def _flow_and_partials(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        cd: float,
    ) -> tuple[float, dict[str, float]]:
        """Signed m_calc (positive from -> to) and its partials, FOLLOWING THE FLOW.

        The kernels read the upstream node: its stagnation state drives, its
        density (incompressible) or (Pt, Tt) (compressible) sets the flow,
        against the downstream static. A reversed orifice used to keep reading
        its declared inlet -- the wrong side's density and drive -- and in the
        compressible regime the kernel swapped a static for a stagnation
        pressure, so it did not converge at all (#554). Forward and reverse
        forms are blended through equal Pt over ORIFICE_REVERSAL_P_SCALE:
        the direction comes from the pressures, so m_calc does not depend on
        m (a flow-based choice jumped m_calc between the two forms as m
        crossed zero at fixed pressures, and hybr stalled on it). Outside
        ~3 Pa either form holds alone, so a forward orifice is evaluated
        exactly as before. An orifice into a dead end (no flow boundary beyond its
        to-node, see _side_is_closed) is never reversed: every root has zero
        flow there.
        """
        # Into a dead end (no flow boundary beyond, a WallNode included) every
        # root has zero flow here: the forward form alone is exact there.
        if getattr(self, "_to_side_closed", False):
            return self._kernel_from(state_in, state_out, self.from_node, self.to_node, cd)
        # A node without a physical state (no positive P and T) is no source.
        if not (float(state_out.P) > 0.0 and float(state_out.Tt) > 0.0):
            return self._kernel_from(state_in, state_out, self.from_node, self.to_node, cd)
        w, dw = _pressure_blend(float(state_in.Pt) - float(state_out.Pt))
        if w == 1.0:
            return self._kernel_from(state_in, state_out, self.from_node, self.to_node, cd)
        m_r, d_r = self._kernel_from(state_out, state_in, self.to_node, self.from_node, cd)
        m_r, d_r = -m_r, {k: -v for k, v in d_r.items()}
        if w == 0.0:
            return m_r, d_r
        m_f, d_f = self._kernel_from(state_in, state_out, self.from_node, self.to_node, cd)
        out: dict[str, float] = {}
        # dict.fromkeys, not a set union: a set of str iterates in hash order,
        # randomised per process, and the key order sets how assemble() sums
        # coinciding terms -- a 1-ulp Jacobian change that moved the path.
        for k in dict.fromkeys([*d_f, *d_r]):
            out[k] = w * d_f.get(k, 0.0) + (1.0 - w) * d_r.get(k, 0.0)
        g = (m_f - m_r) * dw
        out[f"{self.from_node}.Pt"] = out.get(f"{self.from_node}.Pt", 0.0) + g
        out[f"{self.to_node}.Pt"] = out.get(f"{self.to_node}.Pt", 0.0) - g
        return w * m_f + (1.0 - w) * m_r, out

    def _self_consistent_residuals(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None,
        cd_iterate: float,
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        """Row m - m*(x), with m* = Cd(Re(m*)) K(x) solved here (#582).

        K = m_calc / Cd from C++ (m_calc is linear in Cd). Newton on the
        magnitude a = |m*|: g(a) = a - Cd(a) |K|, g' = 1 - |K| dCd/da. At the
        root |K| dCd/da = dlnCd/dlnRe, ~0.2 at its steepest, so g' > 0 and
        the solve is well posed. Sensitivities by the implicit function:
          dm*/dx = (dm_calc/dx|_Cd + K dCd/dx) / (1 - K dCd/dm).
        """
        m_dot = float(state_in.m_dot)
        cd0 = max(cd_iterate, 1e-4)
        K = self._flow_and_partials(state_in, state_out, cd0)[0] / cd0
        sgn = 1.0 if K >= 0.0 else -1.0
        aK = abs(K)
        probe = copy.copy(state_in)
        m_key = f"{self.id}.m_dot"

        def cd_and_slope(a: float, continued: bool = True) -> tuple[float, float]:
            probe.m_dot = a
            cd = self._effective_Cd(probe, state_out, flows)
            if not a > 0:
                return cd, 0.0
            d = self._dCd_dnames(probe, state_out, flows, continue_re_floor=continued)
            return cd, d.get(m_key, 0.0)

        a = aK * cd0
        lo: float | None = None  # an a with g(a) < 0
        hi: float | None = None  # an a with g(a) > 0
        done = False
        for _ in range(50):
            cd, slope = cd_and_slope(a)
            g = a - cd * aK
            if g < 0.0:
                lo = a
            elif g > 0.0:
                hi = a
            gp = 1.0 - aK * slope
            # Newton; off the consistent root (|K| dCd/da < 0.2 there) g' can
            # vanish, where the contracting fixed-point step is used instead.
            step_to = cd * aK if gp <= 0.1 else a - g / gp
            step_to = max(step_to, 0.0)
            done = abs(step_to - a) <= 1e-13 * max(a, 1e-30)
            a = step_to
            if done:
                break
        if not done:
            # The fixed-point step contracts at 1 - g' per iteration: where
            # |K| dCd/da sits just under 1 (a long hole at low Re, steep Cd(Re))
            # 50 steps end up to 1.7% short of the root (#861). g(0) < 0 and,
            # with Cd <= 1, g(|K|) >= 0, so the root is bracketed.
            a = _root_on_bracket(
                lambda t: t - cd_and_slope(t)[0] * aK,
                a,
                lo,
                hi,
                lambda: [0.0, aK, 2.0 * aK, 4.0 * aK],
            )
        # The search above may lean on C++'s slope continued below the Re
        # floor; the sensitivities below are of the HELD value (#760).
        cd_star, cd_m = cd_and_slope(a, continued=False)
        probe.m_dot = sgn * a
        m_star, partial = self._flow_and_partials(probe, state_out, cd_star)
        # Linear in Cd here: the discharge-hole family never recovers (#645),
        # so K above is exact and d(m)/d(Cd) = K.
        partial.pop(_DM_DCD, None)
        denom = 1.0 - aK * cd_m
        # Cd's other dependencies (a duct-fed hole's U1/Vi), at m*.
        # The exit-crossflow term's chamber through-flow reads the ITERATE's
        # own m_dot, not m*: an ordinary column here (_DCD_DM_VIA_CHAMBER).
        for name, dcd in self._dCd_dnames(probe, state_out, flows).items():
            if name != m_key:
                name = _own_flow_key(name, self.id)
                partial[name] = partial.get(name, 0.0) + K * dcd

        row: dict[str, float] = {m_key: 1.0}
        for name, d in partial.items():
            row[name] = row.get(name, 0.0) - d / denom
        return [m_dot - m_star], {0: row}

    def _approach_reads_flow(self) -> bool:
        """Does m_calc read this element's own flow through the approach
        pipe's static state (#651)? With E and a pipe, in either direction."""
        fwd = 0.0 < self.beta < 1.0 and bool(self.upstream_diameter)
        rev = 0.0 < getattr(self, "beta_reverse", 0.0) < 1.0 and bool(
            getattr(self, "_reverse_pipe_diameter", None)
        )
        return fwd or rev

    def _flow_consistent_residuals(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None,
        cd_iterate: float,
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        """Row m - m*(x), with m* = m_calc(m*, x) solved here (#651).

        m_calc reads the element's own flow through the approach pipe's static
        state, p1 = p(Pt, |m| / A_pipe), and, for the discharge-hole family,
        through Cd(Re(m)) as in _self_consistent_residuals (#582). Newton on
        h(a) = a - m_calc(a): the approach makes d(m_calc)/da <= 0, so h' >= 1
        from it; where a rising Cd pulls h' below 0.1 the contracting
        fixed-point step is taken instead. A metering Cd stays at the
        iterate (its documented form), so its m column is carried through.
        Sensitivities by the implicit function:
          dm*/dx = d(m_calc)/dx|_m / (1 - d(m_calc)/dm).
        """
        m_key = f"{self.id}.m_dot"
        own_cd = self._cd_follows_own_flow()
        p_in, p_out = copy.copy(state_in), copy.copy(state_out)

        def evaluate(
            a: float,
            continued: bool = True,
        ) -> tuple[float, float, dict[str, float], float, dict[str, float]]:
            p_in.m_dot = a
            p_out.m_dot = a
            cd = self._effective_Cd(p_in, p_out, flows) if own_cd else cd_iterate
            m, part = self._flow_and_partials(p_in, p_out, cd)
            dm_dcd = part.pop(_DM_DCD, 0.0)
            dcd = (
                self._dCd_dnames(p_in, p_out, flows, continue_re_floor=continued)
                if (own_cd and cd > 0.0)
                else {}
            )
            slope = part.get(m_key, 0.0) + dm_dcd * dcd.get(m_key, 0.0)
            return m, slope, part, dm_dcd, dcd

        a = evaluate(float(state_in.m_dot))[0]
        # The signs of h(a) = a - m_calc(a) seen on the way bracket the root.
        lo: float | None = None  # an a with h(a) < 0
        hi: float | None = None  # an a with h(a) > 0
        done = False
        for _ in range(60):
            m, slope, _, _, _ = evaluate(a)
            h = a - m
            if h < 0.0:
                lo = a
            elif h > 0.0:
                hi = a
            hp = 1.0 - slope
            a_new = m if hp <= 0.1 else a - h / hp
            done = abs(a_new - a) <= 1e-13 * max(abs(a_new), 1e-30)
            a = a_new
            if done:
                break
        if not done:
            # m_calc(a) can carry a kink: the approach pipe's choke (flat
            # beyond it), Jobson's held discriminant (beta past its range).
            # Newton across one 2-cycles to the cap, and the a it stopped at
            # is no root: the row was off by the cycle's width (6% on a
            # stress bleed station; a GUI net's orifice passed 33 kg/s where
            # its consistent flow was 38, #861) and its Jacobian matched
            # neither side. With a Cd rising steeply in Re the fixed-point
            # step also stops short of the root, approaching from one side.
            # Solve on a bracket instead, only when the stop is no root, so
            # every solve that met its root before is untouched.
            def h_of(t: float) -> float:
                return t - evaluate(t)[0]

            def probes() -> list[float]:
                m0 = evaluate(0.0)[0]  # the flow at zero approach velocity
                return [0.0, m0, 2.0 * m0, 4.0 * m0]

            a = _root_on_bracket(h_of, a, lo, hi, probes)
        a = self._drive_signed_root(evaluate, a, float(state_in.Pt) - float(state_out.Pt))
        # Sensitivities of the HELD Cd below its Re floor (#760).
        m_star, slope, part, dm_dcd, dcd = evaluate(a, continued=False)
        hp = 1.0 - slope
        part.pop(m_key, None)
        if own_cd:
            cd_cols: dict[str, float] = {}
            for k, v in dcd.items():
                if k != m_key:
                    k = _own_flow_key(k, self.id)
                    cd_cols[k] = cd_cols.get(k, 0.0) + v
        else:
            cd_cols = self._dCd_dnames(state_in, state_out, flows) if cd_iterate > 0.0 else {}
        row: dict[str, float] = {m_key: 1.0}
        for name, d in part.items():
            row[name] = row.get(name, 0.0) - d / hp
        for name, dc in cd_cols.items():
            row[name] = row.get(name, 0.0) - dm_dcd * dc / hp
        return [float(state_in.m_dot) - m_star], {0: row}

    @staticmethod
    def _drive_signed_root(
        evaluate: Callable[[float], tuple[float, ...]],
        a: float,
        dPt: float,
    ) -> float:
        """The consistent flow on the half-line of its drive (#875).

        m_calc reads |m| through the approach pipe, so h(a) = a - m_calc(a)
        has a second, spurious root against the drive: an iterate above the
        flux where the approach static falls to the downstream pressure
        gives m_calc < 0, the solve follows it below zero, and the approach
        pipe's choke (flat m_calc there) holds it at a reversed flow tens of
        times the true one, where the row reads no pressure at all. The
        physical root has the sign of m0 = m_calc(0), the flow at zero
        approach velocity. Tested only where a disagrees with the
        stagnation difference ``dPt`` (cheap), m0 decides; a root of m0's
        sign keeps its bits. Otherwise h is solved on [0, ...] probing out
        from m0, where h(0) = -m0 and the approach only lowers m_calc."""
        if a != 0.0 and (a > 0.0) == (dPt > 0.0):
            return a
        m0 = float(evaluate(0.0)[0])
        if m0 == 0.0 or a == 0.0 or (a > 0.0) == (m0 > 0.0):
            return a

        def h(t: float) -> float:
            return t - float(evaluate(t)[0])

        end = 0.0
        for k in (1.0, 2.0, 4.0, 8.0, 16.0):
            t = k * m0
            ht = h(t)
            if ht == 0.0:
                return t
            if (ht > 0.0) == (m0 > 0.0):
                lo, hi = (end, t) if m0 > 0.0 else (t, end)
                return float(_brentq(h, lo, hi, xtol=1e-13 * max(abs(t), 1e-30)))
            end = t
        return a

    def n_equations(self) -> int:
        return 1

    def diagnostics(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
    ) -> dict[str, float | str]:
        if state_in.P <= 0 or state_in.Tt <= 0:
            return {}
        if self.closed:
            # Closed (#1024): no flow, no jet, no Cd. The pressures either side
            # are reported; the drop across it is held by the wall.
            return {
                "m_dot": float(state_in.m_dot),
                "closed": 1.0,
                **_element_pressure_block(state_in, state_out, mach_in=0.0, mach_out=0.0),
                "area": 0.0,
                "velocity": 0.0,
                "mach_throat": 0.0,
                "diameter": 0.0,
                "beta": 0.0,
                "area_source": self._area_source,
            }
        cs = cb.complete_state(state_in.T, state_in.P, state_in.X)

        # Throat Mach: isentropic pressure ratio from the UPSTREAM total to the
        # downstream static -- upstream along the flow, which may run against
        # the declared direction (#554).
        gamma = cs.thermo.gamma
        pr_crit = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
        up, dn = (state_out, state_in) if state_out.Pt > state_in.Pt else (state_in, state_out)
        p_ratio_throat = dn.P / up.Pt if up.Pt > 0 else 1.0
        if p_ratio_throat >= 1.0:
            mach_throat = 0.0
        elif p_ratio_throat <= pr_crit:
            mach_throat = 1.0
        else:
            mach_throat = float(cb.mach_from_pressure_ratio(up.Tt, up.Pt, dn.P, up.X))

        # Physical velocity / Mach at the geometric bore (not Cd-adjusted)
        Dh = math.sqrt(4.0 * self.area / math.pi) if self.area > 0 else 0.0
        ref = _element_reference_block(
            cs, m_dot=state_in.m_dot, area=self.area, Dh=Dh, location="inlet"
        )
        v_in = ref["velocity"]
        mach_in = v_in / cs.thermo.a if cs.thermo.a > 0 else 0.0

        # Outlet Mach (use outlet state for accuracy)
        try:
            cs_out = cb.complete_state(state_out.T, state_out.P, state_out.X)
            v_out = (
                (abs(state_in.m_dot)) / ((cs_out.thermo.rho) * (self.area))
                if self.area > 0 and cs_out.thermo.rho > 0
                else 0.0
            )
            mach_out = v_out / cs_out.thermo.a if cs_out.thermo.a > 0 else 0.0
        except Exception:
            mach_out = 0.0

        cd = float(self._effective_Cd(state_in, state_out, flows))
        return {
            "m_dot": float(state_in.m_dot),
            **_element_pressure_block(state_in, state_out, mach_in=mach_in, mach_out=mach_out),
            **ref,
            "mach_throat": float(mach_throat),
            "Cd": cd,
            "is_correlation": float(self.use_correlation),
            "Cd_in_range": float(self._hole_cd.in_range(state_in)),
            **self._hole_cd.envelope(state_in),
            "diameter": float(self.diameter or 0.0),
            "beta": float(self.beta),
            "area_source": self._area_source,
            **self._tap_diagnostics(state_in, state_out, cd),
            **self._model_diagnostics(state_in, state_out, cd),
        }

    def _model_diagnostics(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        cd: float,
    ) -> dict[str, float | str]:
        """The compressible model's state along the flow (#854).

        Empty for the nozzle and in the incompressible regime, so a nozzle
        orifice reports exactly what it did before #854. Ci is the incompressible
        contraction the model used (Jobson) or C (ISO); contraction_C is
        Jobson's compressible C; expansibility is m over the incompressible
        Cd E A sqrt(2 rho1 dP) at the approach static state; model_in_range
        is 1 inside the model's stated range (Jobson: 0.5 <= Ci <= 0.7, beta
        <= 0.75; ISO: p2/p1 >= 0.75). A flag, never a refusal.
        """
        model = getattr(self, "compressible_model", "nozzle")
        if model == "nozzle" or self.regime != "compressible" or self._pass_through:
            return {}
        forward = not (state_out.Pt > state_in.Pt)
        up, dn = (state_in, state_out) if forward else (state_out, state_in)
        up_id = self.from_node if forward else self.to_node
        beta = self.beta if forward else getattr(self, "beta_reverse", 0.0)
        recovers = self._recovery_fwd if forward else self._recovery_rev
        ratio = getattr(self, "_pipe_cd_ratio_fwd" if forward else "_pipe_cd_ratio_rev", 1.0)
        c = max(cd * ratio, _JOBSON_CD_MIN)
        p1 = self._approach_static(up, up_id, forward, beta)[0]
        p_down = float(dn.P)
        if recovers:
            p_down = float(
                _solver_tools.orifice_tap_pressure_and_jacobian(
                    p1, float(up.Pt), float(dn.Pt), beta, c, True
                ).p_tap
            )
        args = (float(up.m_dot), float(up.Pt), p1, float(up.Tt), up.Y, p_down, c, self.area, beta)
        if model == "jobson1955":
            r = _solver_tools.orifice_jobson_residuals_and_jacobian(*args, self.Cv)
        else:
            r = _solver_tools.orifice_iso5167_residuals_and_jacobian(*args)
        out: dict[str, float | str] = {
            "compressible_model": model,
            "Ci": float(r.Ci),
            "expansibility": float(r.expansibility),
            "model_in_range": float(r.in_range),
        }
        if model == "jobson1955":
            out["contraction_C"] = float(r.contraction)
            out["Cv"] = self.Cv
        return out

    def _tap_diagnostics(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        cd: float,
    ) -> dict[str, float | str]:
        """What a meter at ``self.tap`` reads, and the loss the row charged.

        beta_flow and E are the velocity of approach along the flow, which
        reverses with it (#554). dP_tap is the tapping differential at the
        element's ``tap``. The row
        works at _NETWORK_TAP, so for a metering correlation the reported one
        is scaled by (C_network / C_tap)^2 at the same flow. Between pipes
        (#645) dP_permanent is the stagnation drop the row charged and
        loss_ratio = dP_permanent / dP_tap at the network taps;
        recovery_in_range is False where p_tap/p1 < 0.75, outside ISO's
        incompressible relation.
        """
        forward = not (state_out.Pt > state_in.Pt)
        up, dn = (state_in, state_out) if forward else (state_out, state_in)
        beta = self.beta if forward else getattr(self, "beta_reverse", 0.0)
        recovers = self._recovery_fwd if forward else self._recovery_rev
        ratio = getattr(self, "_pipe_cd_ratio_fwd" if forward else "_pipe_cd_ratio_rev", 1.0)
        out: dict[str, float | str] = {
            "tap": self.tap,
            "pressure_recovery": float(recovers),
            # The approach along the flow: beta (above) is the declared one.
            "beta_flow": float(beta),
            "E": float(cb.velocity_of_approach_factor(beta)),
            # C_pipe / Cd: a discharge-hole Cd carried into its pipe (#645).
            "Cd_pipe_ratio": float(ratio),
        }
        cd = cd * ratio
        up_id = self.from_node if forward else self.to_node
        p1 = self._approach_static(up, up_id, forward, beta)[0]
        if 0.0 < beta < 1.0:
            # The approach pipe's static pressure the row drives from (#651).
            out["P_approach"] = p1
        p_up = p1 if 0.0 < beta < 1.0 else float(up.Pt)
        if recovers and cd > 0.0:
            t = _solver_tools.orifice_tap_pressure_and_jacobian(
                p1, float(up.Pt), float(dn.Pt), beta, cd, self.regime == "compressible"
            )
            dp_net = p_up - float(t.p_tap)
            out["dP_permanent"] = float(up.Pt) - float(dn.Pt)
            out["loss_ratio"] = float(t.loss_ratio)
            out["recovery_in_range"] = float(t.in_range)
        else:
            dp_net = p_up - float(dn.P)
        metering = _METERING_SELECTORS.get(self.correlation)
        if (
            metering is not None
            and self._orifice_geom is not None
            and _TAPS[self.tap] != _NETWORK_TAP
        ):
            g = self._orifice_geom
            at_tap = cb.OrificeGeometry()
            at_tap.d, at_tap.D, at_tap.t, at_tap.r = g.d, g.D, g.t, g.r
            at_tap.tap = _TAPS[self.tap]
            mu = cb.transport_state(state_in.T, state_in.P, state_in.X).mu
            m = abs(float(state_in.m_dot))
            Re_D = 4.0 * m / (math.pi * g.D * mu) if (mu > 0.0 and m > 1e-12) else 1e4
            c_net = float(cb.metering_cd_and_dRe(metering, g, Re_D)[0])
            c_tap = float(cb.metering_cd_and_dRe(metering, at_tap, Re_D)[0])
            dp_net *= (c_net / c_tap) ** 2
        out["dP_tap"] = dp_net
        return out


def _effusion_wall(
    h_i: float, T_c: float, h_gas: float, T_aw_gas: float, t_over_k: float
) -> tuple[float, float, float]:
    """(T_wall_hot, T_wall_cold, q) of a plate between gas and coolant.

    Per unit area, one conduction layer, through C++'s wall coupling. The one
    wall core both EffusionPlateElement.overall_effectiveness (no conduction)
    and wall_heat_transfer use.
    """
    w = _solver_tools.wall_coupling_and_jacobian_multilayer(
        h_gas, T_aw_gas, h_i, T_c, [t_over_k], 1.0, 0.0
    )
    q = float(w.Q)
    return float(w.T_hot), T_c + q / h_i, q


#: The numeric keys EffusionPlateElement.wall_heat_transfer reports on a
#: gas-side chamber (besides gas_side_chamber and the coolant-side ones).
_EFFUSION_CHAMBER_SIDE_KEYS = (
    "h_gas_unblown",
    "h_gas_unblown_Re",
    "h_gas_unblown_Pr",
    "h_gas_unblown_Nu",
    "gas_augmentation",
    "h_gas",
    "T_gas",
    "T_gas_from_main_inlet",
    "eta_film",
    "film_extrapolated",
    "T_adiabatic_wall",
    "q_wall",
    "Q_wall",
    "T_wall_hot",
    "T_wall_cold",
    "eta_overall",
    "U_gas",
    "blowing_ratio",
    "density_ratio",
    "velocity_ratio",
)

#: The numeric keys a flat-plate gas side (`gas_baseline="flat_plate"`) adds
#: to _EFFUSION_CHAMBER_SIDE_KEYS: where the restarted boundary layer's panel
#: average was taken, and whether its Re_x left the correlation's data span.
_EFFUSION_FLAT_PLATE_KEYS = (
    "h_gas_unblown_x_start",
    "h_gas_unblown_x_end",
    "h_gas_unblown_Re_x_start",
    "h_gas_unblown_Re_x_end",
    "h_gas_unblown_extrapolated",
)

#: What `h_gas_unblown_model` names for the flat-plate baseline.
_FLAT_PLATE_MODEL = "flat_plate_turbulent(CR-2648 St_o, panel average)"


def _chamber_gas_side(node: Any, state: "NetworkMixtureState") -> dict[str, Any] | None:
    """What a discharge node offers a wall as its gas side, or None.

    Only a MomentumChamberNode has flow, hence a coefficient: its own surface
    correlation (the user's choice on that node) is the UNBLOWN h, and its
    velocity the gas velocity. The gas TEMPERATURE is the stream the wall
    sees: a merge chamber's state is the MIXED outlet, already diluted by
    the side streams -- with liner effusion at 20-30% of the flow that is a
    large, unconservative drop -- so when the chamber has a main inlet the
    approaching main stream's temperature is used instead. A plenum (state,
    no flow) returns None.
    """
    if not isinstance(node, MomentumChamberNode) or node.area <= 0.0:
        return None
    # No gas flow below the network's no-through-flow floor (#792, #998): a
    # zero-flow root leaves ~1e-13 kg/s in the chamber, whose coefficient
    # would make a stagnant chamber an adiabatic gas side.
    if not abs(float(getattr(node, "_total_m_dot", 0.0))) >= cb.MIXER_HEAT_MDOT_FLOOR:
        return None
    gas = node._gas_side_htc_and_T(state)
    if gas is None:
        return None
    # The chamber's static density (#578), from its closure on state.Tt.
    rho_g, _ = _safe_rho(node.static_density(state)[0])
    T_gas, from_main, Y_gas = float(gas.T_aw), False, list(state.Y)
    if node.main_inlet is not None and node.has_side_streams():
        for eid, st, _names in getattr(node, "_inflows", []):
            if eid == node.main_inlet:
                T_gas, from_main, Y_gas = float(st.Tt), True, list(st.Y)
    model = node.surface.model
    corr = getattr(model, "correlation", None)
    return {
        "h_unblown": float(gas.h),
        # What the augmentation multiplies (#476): the chamber surface's own
        # model, at the chamber's Re and Pr. The smooth-duct correlations are
        # fully developed (channel_smooth takes no entry-length correction).
        "baseline_model": f"{type(model).__name__}({corr})" if corr else type(model).__name__,
        "Re_unblown": float(gas.Re),
        "Pr_unblown": float(gas.Pr),
        "Nu_unblown": float(gas.Nu),
        "T_gas": T_gas,
        "T_from_main_inlet": from_main,
        "Y_gas": Y_gas,
        "U_gas": abs(getattr(node, "_total_m_dot", 0.0)) / (rho_g * node.area),
        "rho_gas": rho_g,
    }


class EffusionPlateElement(OrificeElement):
    """A multi-perforated (effusion) wall panel: N holes discharging in parallel.

    Coolant enters from `from_node` and leaves through the wall into
    `to_node`, so mass leaves the coolant circuit -- which is what makes this
    an element rather than a `ConvectiveSurface` model.

    GEOMETRY IS GIVEN THE WAY A PLATE IS DESIGNED, not as an area: pitch,
    hole diameter, wall thickness and inclination angle. The hole count
    follows from the panel area and the pitch, and is rounded to a whole
    number of holes -- `hole_count_exact` and `pitch_actual` record what the
    rounding did, because the rounded count is what the element flows.

        n_holes  = round(panel_area / (pitch_x * pitch_y))
        L_hole   = wall_thickness / sin(angle)     <- the drilled length
        A_total  = n_holes * pi * d^2 / 4
        porosity = A_total / panel_area

    HOMOGENISATION, AND WHEN IT BREAKS. One element is one panel with one
    coolant pressure and one gas pressure, so it cannot represent coolant
    migration WITHIN itself. van de Noort and Ireland (2022) show this is not
    a small effect: with uniform inlet AND outlet pressure, holes at the end
    of an array can pass ~75% of what a central hole passes, and under a
    spanwise pressure gradient one row took ~10% of the total while others
    took ~20% each.

    The answer is to use more panels, not a cleverer one. Because the solver
    balances mass over every element at a node, a panel hung off each segment
    of a coolant channel reproduces their flow network directly:

        coolant:  P0 --[Channel]-- n1 --[Channel]-- n2 --[Channel]-- ...
                                    |               |
                            [EffusionPanel]  [EffusionPanel]
                                    |               |
        gas:                       g1              g2

    At each node, `m_channel_in = m_channel_out + m_effusion`, which is the
    coolant mass flow falling along the wall -- the effusion CHANNEL case,
    with no channel-specific element needed. Resolution is the caller's
    choice of segment count. Compare `dP_drive` across panels to see whether
    one panel is smearing a gradient that deserves several.

    INGESTION. If the gas pressure exceeds the coolant pressure the panel
    ingests hot gas. That is a real failure mode (van de Noort's CMF > 0.5),
    and a homogenised panel would otherwise average it into a healthy net
    outflow, so `diagnostics()` reports `is_ingesting` rather than staying
    silent about it.

    DISCHARGE COEFFICIENT. Default `IdelchikThick`: diagram 4-18a is a
    thick-walled hole in a large wall between two plena, which is exactly a
    plenum-fed effusion plate, and it is valid down to Re = 25. For a panel
    fed by a channel rather than a plenum the approach flow is a CROSSFLOW,
    which is McGreehan-Schotsch's `U1/Vi` term -- use `'McGreehanSchotsch'`
    there. Note that neither carries a velocity-of-approach `beta` factor:
    a plenum has no approach velocity to correct for. van de Noort applies
    one built on a half-pitch square inlet area, which is worth knowing about
    but is a different convention; at their own pitch it is a 0.3% effect.

    THE HOLE ANGLE enters the Cd through the drilled length t/sin(angle) and,
    with IdelchikThick, through Idelchik diagram 3-2's inclined-inlet loss,
    dzeta = 0.3 cos(angle) + 0.2 cos^2(angle), added to diagram 4-18a's
    inlet coefficient (Eq. 4-20); reported as ``Cd_inclination_dzeta``
    (exactly 0 at 90 deg). ``HoleCd(inclination='auto')`` (default)
    applies it to an inclined hole with IdelchikThick (a 90 deg hole
    resolves to 'none') and leaves it off with any other correlation;
    ``'idelchik_3_2'`` forces it (any other correlation refuses it) and
    ``'none'`` switches it off. Outside its range (Re > 1e4, l/d >= 2) the
    Cd is flagged by ``Cd_in_range``, not refused.

    THE CD's MODEL CHOICES are one spec, ``hole_cd=HoleCd(...)``: the
    correlation (default IdelchikThick), a fixed Cd, the inclination term,
    the supply crossflow ("auto": McGreehan-Schotsch's U1/Vi on a duct-fed
    panel) and the exit crossflow ("auto": none; 'rowbury_2001' opts in).
    The keywords ``correlation`` and ``Cd`` are deprecated since 0.12.0
    (removed in 0.13.0); ``inlet_crossflow='none'`` holds U1/Vi at 0 on a
    duct-fed panel (no term) and lifts the refusal of a correlation without
    the crossflow term. ``exit_crossflow_region`` and ``exit_crossflow_orientation_deg`` are
    the hole's geometry and stay here.

    THERMAL, INTERNAL SIDE. `internal_heat_transfer()` gives the coolant-side
    coefficient from Andrews 86-GT-225: the hole APPROACH flow over the
    coolant-side surface plus the THROAT, summed. The approach term dominates
    at effusion Reynolds numbers, which is the paper's own headline and the
    reason a throat-only treatment recovers only a fraction of the measured
    coefficient.

    THE PLATE OWNS ITS WALL (#471). `wall_heat_transfer()` -- reported in
    the diagnostics -- solves the wall with the coolant side above and a gas
    side the DISCHARGE NODE decides: a MomentumChamberNode (flow) supplies
    the gas temperature and velocity, and by default its own surface
    correlation as the unblown coefficient; a plenum (state, no flow) takes
    the imposed `gas_heat_flux`. No ThermalWall connects to this element.

    THE UNBLOWN BASELINE IS SELECTABLE (`gas_baseline`). 'chamber' (default)
    is the chamber's own correlation -- for a smooth duct, fully developed.
    'flat_plate' is a turbulent boundary layer restarted `gas_start_length`
    upstream of the panel (a bleed, a starter-film lip, a liner step),
    averaged over `panel_length`: CR-2648's St_o, the unblown reference the
    Stanford full-coverage St/St_o data -- the measured augmentation -- are
    published against. On two rigs a duct correlation was the wrong
    reference: Andrews 88-GT-290 (plate at x/D_h 1.5, entry-corrected h_0
    1.76x Dittus-Boelter) and Behrendt 2008 (boundary layer bled 21 d before
    the sample; the duct's Dh 319 mm gives about half the flat plate's h).
    Neither makes a `gas_augmentation` value a default; which baseline a
    literature augmentation is quoted against decides which one to pick.

    AN OVERALL EFFECTIVENESS CORRELATION IS NEVER THE CLOSURE. Andrews
    (88-GT-290) states internal and film cooling are NOT additive -- "the
    film cooling reduces the mean gas temperature adjacent to the wall...
    which in turn reduces the heat flux removed by the internal wall
    cooling" -- and an overall correlation already contains the internal
    convection computed here. The overall effectiveness is an OUTPUT. An
    adiabatic film (Baldauf + Sellers, `gas_film='baldauf_sellers'`) is
    selectable, off by default: scored on Andrews, the data refuse it offered
    alone; the missing physics is the gas-side augmentation
    (`gas_augmentation`, the caller's).

    PLENUM-FED. This 2-port plate leaves McGreehan-Schotsch's U1/Vi at 0 and
    uses Andrews' still-plenum coolant side (as does a duct-fed panel with
    HoleCd(inlet_crossflow="none"), which still reports the duct's U1/Vi);
    `coolant_crossflow_ignored`
    flags a supply node a channel runs through.

    GAS RADIATION (#1071), opt-in. `radiation=GasRadiation(beam_length=...)`
    adds Lefebvre & Ballal's Eq. 8.7 (emissivity Eq. 8.11) on the hot face,
    from the gas the convection sees (T_gas, the chamber's static P, that
    stream's composition, unless the spec fixes them), through the same C++
    kernel as a ThermalWall face (`wall_radiation_sides`). Off (default), the
    wall solve is today's. Reported: `q_rad_hot`, `eps_g`.

    WALL TEMPERATURE LIMIT (#998). `material` names the plate's alloy from
    the database (`cb.list_materials()`) for ONE purpose: its #882 envelope's
    hot end (the solidus) is `T_wall_limit`, and `wall_above_limit` flags a
    face hotter than that. Without a material the limit is
    `cb.GENERIC_METAL_SOLIDUS`, the hot end #882 gives a metal with no printed
    solidus. The conduction stays `wall_conductivity`, a constant k.

    NO COOLING PATH (#998). With no coolant flow and no gas-side flow (a
    plenum, or a chamber with none) nothing takes heat from the plate, so its
    temperature is undetermined: `no_cooling_path` is 1 and the solver
    refuses the root, naming the plate.
    """

    _selects_compressible_model = False  # keeps the nozzle (#861)
    # The deprecated hole-Cd keywords (HoleCd replaces them): not serialized.
    _deprecated_init_kwargs = ("correlation", "Cd")

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        hole_diameter: float,
        wall_thickness: float,
        pitch: float | None = None,
        panel_area: float | None = None,
        pitch_x: float | None = None,
        pitch_y: float | None = None,
        panel_length: float | None = None,
        panel_width: float | None = None,
        angle_deg: float = 90.0,
        correlation: str | None = None,
        Cd: float | None = None,
        edge_radius: float = 0.0,
        internal_Nu_multiplier: float = 1.0,
        wall_conductivity: float = 20.0,
        gas_film: Literal["none", "baldauf_sellers"] = "none",
        gas_augmentation: float = 1.0,
        turbulence_intensity: float = 0.05,
        gas_heat_flux: float = 0.0,
        crossflow_segments: tuple[str | None, str | None] | None = None,
        crossflow_area: float | None = None,
        material: str | None = None,
        gas_baseline: Literal["chamber", "flat_plate"] = "chamber",
        gas_start_length: float | None = None,
        radiation: GasRadiation | None = None,
        exit_crossflow_region: Literal[
            "suction_surface", "pressure_surface", "leading_edge"
        ] = "suction_surface",
        exit_crossflow_orientation_deg: float = 0.0,
        hole_cd: HoleCd | None = None,
    ) -> None:
        # The hole Cd's model choices (hole_cd.py): a HoleCd, or the
        # deprecated per-element keywords; the correlation defaults to
        # IdelchikThick.
        spec = _with_correlation(
            _hole_cd_spec(
                f"EffusionPlateElement {id!r}",
                hole_cd,
                {"correlation": correlation, "Cd": Cd},
                stacklevel=3,
            ),
            "IdelchikThick",
            f"EffusionPlateElement {id!r}",
        )
        if hole_diameter <= 0.0:
            raise ValueError("EffusionPlateElement: hole_diameter must be positive")
        # Gas radiation onto the hot face (#1071), opt-in: the same spec and
        # kernel as a ThermalWall face (#793). A casing faces no effusion
        # plate's hot side, so only GasRadiation.
        if radiation is not None and not isinstance(radiation, GasRadiation):
            raise ValueError(
                "EffusionPlateElement: radiation is a GasRadiation (the hot face "
                f"sees the gas), not {type(radiation).__name__}"
            )
        self.radiation = radiation
        # WHICH UNBLOWN COEFFICIENT `gas_augmentation` MULTIPLIES (#476). The
        # chamber's own surface correlation (default), or a turbulent boundary
        # layer restarted `gas_start_length` [m] upstream of the panel's
        # leading edge -- a bleed, a starter-film lip, a liner step -- and
        # averaged over the panel. See wall_heat_transfer.
        if gas_baseline not in ("chamber", "flat_plate"):
            raise ValueError("EffusionPlateElement: gas_baseline is 'chamber' or 'flat_plate'")
        if gas_baseline == "flat_plate":
            if gas_start_length is None or not gas_start_length >= 0.0:
                raise ValueError(
                    "EffusionPlateElement: gas_baseline='flat_plate' needs "
                    "gas_start_length >= 0, the distance [m] from where the hot-side "
                    "boundary layer starts (a bleed, a starter-film lip, a step) to "
                    "the panel's leading edge"
                )
            if panel_length is None:
                raise ValueError(
                    "EffusionPlateElement: gas_baseline='flat_plate' needs "
                    "panel_length, the streamwise extent the coefficient is averaged over"
                )
        elif gas_start_length is not None:
            raise ValueError(
                "EffusionPlateElement: gas_start_length is the flat-plate baseline's "
                "input; the chamber baseline takes the chamber's own correlation "
                "(set gas_baseline='flat_plate', or drop gas_start_length)"
            )
        self.gas_baseline = gas_baseline
        self.gas_start_length = None if gas_start_length is None else float(gas_start_length)
        # The plate's alloy, for its temperature limit only (#998): the hot
        # end of its #882 envelope, or the generic metal cap without one.
        if material is not None and material not in cb.list_materials():
            raise ValueError(
                f"EffusionPlateElement: unknown material {material!r}; expected one "
                f"of {cb.list_materials()}, or None for the generic metal cap."
            )
        self.material = material
        self.T_wall_limit = (
            float(cb.get_material_info(material)["T_envelope_max"])
            if material is not None
            else float(cb.GENERIC_METAL_SOLIDUS)
        )
        if crossflow_segments is not None and not (crossflow_area and crossflow_area > 0.0):
            raise ValueError(
                "EffusionPlateElement: a crossflow-fed panel needs the backside "
                "duct's crossflow_area"
            )
        # DUCT-FED (#471): the panel's supply node is a station of a backside
        # duct; crossflow_segments = (arriving, leaving) segment ids there. The
        # duct's mean velocity at the station, over the static-referenced ideal
        # jet velocity, is McGreehan-Schotsch's U1/Vi.
        self.crossflow_segments = crossflow_segments
        self.crossflow_area = crossflow_area
        self._xf_elems: list[tuple[NetworkElement, str]] = []
        if wall_conductivity <= 0.0:
            raise ValueError("EffusionPlateElement: wall_conductivity must be positive")
        if gas_augmentation <= 0.0:
            raise ValueError("EffusionPlateElement: gas_augmentation must be positive")
        if gas_film not in ("none", "baldauf_sellers"):
            raise ValueError("EffusionPlateElement: gas_film is 'none' or 'baldauf_sellers'")
        if gas_film == "baldauf_sellers" and panel_length is None:
            raise ValueError(
                "EffusionPlateElement: gas_film='baldauf_sellers' needs panel_length "
                "to count the hole rows the film builds over"
            )
        # The plate's own wall (#471): conduction through it, and the gas
        # side the discharge node decides -- see wall_heat_transfer().
        self.wall_conductivity = float(wall_conductivity)
        self.gas_film = gas_film
        self.gas_augmentation = float(gas_augmentation)
        self.turbulence_intensity = float(turbulence_intensity)
        self.gas_heat_flux = float(gas_heat_flux)
        self.panel_length = panel_length
        if internal_Nu_multiplier <= 0.0:
            raise ValueError("EffusionPlateElement: internal_Nu_multiplier must be positive")
        if wall_thickness <= 0.0:
            raise ValueError("EffusionPlateElement: wall_thickness must be positive")
        if not 0.0 < angle_deg <= 90.0:
            raise ValueError(
                "EffusionPlateElement: angle_deg must be in (0, 90]; it is the "
                "hole inclination to the wall PLANE, 90 being a normal hole"
            )

        # The coolant-side tuner, matching `Nu_multiplier` on ConvectiveSurface.
        # Andrews 86-GT-225's correlation runs -10.4% against his own Fig. 8,
        # so a user matching their own plate needs this; the harness never
        # sets, recommends or scores it. See docs/VALIDATION_POLICY.md.
        self.internal_Nu_multiplier = float(internal_Nu_multiplier)

        px = pitch_x if pitch_x is not None else pitch
        py = pitch_y if pitch_y is not None else pitch
        if px is None or py is None:
            raise ValueError(
                "EffusionPlateElement: give pitch (square array) or both pitch_x and pitch_y"
            )
        if px <= 0.0 or py <= 0.0:
            raise ValueError("EffusionPlateElement: pitch must be positive")

        if panel_area is None:
            if panel_length is None or panel_width is None:
                raise ValueError(
                    "EffusionPlateElement: give panel_area, or both panel_length and panel_width"
                )
            panel_area = panel_length * panel_width
        if panel_area <= 0.0:
            raise ValueError("EffusionPlateElement: panel_area must be positive")

        cell_area = px * py
        n_exact = panel_area / cell_area
        n_holes = int(round(n_exact))
        if n_holes < 1:
            raise ValueError(
                f"EffusionPlateElement: the panel holds {n_exact:.3g} holes at "
                f"this pitch, which rounds to none. Enlarge the panel or "
                f"reduce the pitch."
            )

        self.hole_diameter = hole_diameter
        self.wall_thickness = wall_thickness
        self.pitch_x = px
        self.pitch_y = py
        self.panel_area = panel_area
        self.angle_deg = angle_deg
        self.n_holes = n_holes
        # What the rounding cost: the caller asked for a pitch, and the whole
        # number of holes implies a slightly different one.
        self.hole_count_exact = n_exact
        self.pitch_actual = math.sqrt(panel_area / n_holes)

        # Drilled length along the hole axis. An inclined hole is longer than
        # the wall is thick, which is the whole reason effusion holes are
        # inclined: more internal surface for the same wall.
        alpha = math.radians(angle_deg)
        self.hole_length = wall_thickness / math.sin(alpha)

        hole_area = math.pi * hole_diameter * hole_diameter / 4.0
        self.total_hole_area = n_holes * hole_area
        self.porosity = self.total_hole_area / panel_area

        super().__init__(
            id,
            from_node,
            to_node,
            Cd=spec.Cd,
            area=self.total_hole_area,
            correlation=spec.correlation,
            plate_thickness=self.hole_length,
            edge_radius=edge_radius,
        )
        # OrificeElement derived an equivalent single-bore diameter from the
        # TOTAL area. Keep it for the flow equation, but the correlations must
        # see one real hole -- see _hole_diameter_for_correlation.
        self.equivalent_bore = self.diameter
        # Discharged into a merge chamber, the holes' inclination to the wall
        # is the jets' angle to the gas axis (holes pointing downstream).
        self.injection_angle_deg = float(angle_deg)
        # The hole Cd's options (hole_cd.py): the inclined-inlet term
        # (CD_INCLINATION_MODELS, "auto" resolved: an effusion hole is
        # cylindrical), the opt-in exit-crossflow term
        # (CD_EXIT_CROSSFLOW_MODELS), and a duct-fed panel's supply crossflow
        # ("auto": read where the panel is duct-fed). ``cd_inclination`` and
        # ``cd_exit_crossflow`` hold the models in use, "auto" resolved.
        self.exit_crossflow_region = exit_crossflow_region
        self.exit_crossflow_orientation_deg = float(exit_crossflow_orientation_deg)
        self.cd_inclination, self.cd_exit_crossflow = self._hole_cd.configure_cooling_hole(
            f"EffusionPlateElement {id!r}",
            spec,
            hole_type="cylindrical",
            angle_deg=angle_deg,
            exit_crossflow_region=exit_crossflow_region,
            exit_crossflow_orientation_deg=self.exit_crossflow_orientation_deg,
            exit_area_ratio=1.0,
            duct_fed=crossflow_segments is not None,
        )
        self.hole_cd = dataclasses.replace(spec, Cd=self.Cd)

    def _hole_count(self) -> float:
        return float(self.n_holes)

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """No upstream-diameter discovery: a wall panel has no pipe, hence no
        beta. The geometry the correlations need is the single hole, which is
        known at construction. The discharge node decides the gas side."""
        self._orifice_geom = cb.OrificeGeometry()
        self._orifice_geom.d = self.hole_diameter
        self._orifice_geom.D = 0.0
        self._orifice_geom.t = self.hole_length
        self._orifice_geom.r = self.edge_radius
        self.beta = 0.0
        self._discharge_node = graph.nodes.get(self.to_node)
        self._hole_cd.resolve_exit_mainstream(graph)
        # A supply node a channel runs THROUGH feeds the holes with a
        # crossflow. This 2-port plate is plenum-fed: McGreehan-Schotsch's
        # U1/Vi stays 0 and Andrews' coolant side assumes a still plenum, so
        # say so rather than model it (the 3-port liner segment will).
        through = (
            [
                e
                for e in graph.get_upstream_elements(self.from_node)
                + graph.get_downstream_elements(self.from_node)
                if isinstance(e, ChannelElement)
            ]
            if self.from_node in graph.nodes
            else []
        )
        # And a duct-fed panel whose supply-crossflow term is switched off
        # (HoleCd(inlet_crossflow="none")): the caller's choice, flagged.
        self._coolant_crossflow_ignored = (
            len(through) >= 2 and self.crossflow_segments is None
        ) or (self.crossflow_segments is not None and not self._hole_cd.supply_duct)
        self._xf_elems = []
        if self.crossflow_segments is not None:
            for eid in self.crossflow_segments:
                if eid is not None:
                    self._xf_elems.append((graph.elements[eid], eid))
        self._resolve_plate_against_chambers(graph)

    def network_flow_inputs(self) -> list[tuple[str, str]]:
        """The backside segments' flows at this panel's supply station."""
        return [(eid, self.from_node) for _, eid in self._xf_elems]

    def validate(self) -> None:
        super().validate()
        if self.correlation in _METERING_SELECTORS:
            raise ValueError(
                f"EffusionPlateElement {self.id!r}: {self.correlation!r} is a "
                "NORMED metering correlation for a standardised plate in a "
                "pipe, referenced to the tapping differential and a function "
                "of beta = d/D. An effusion panel has no pipe and no beta. "
                "Use 'IdelchikThick' (plenum-fed), 'McGreehanSchotsch' "
                "(channel-fed, with approach crossflow), or 'fixed'."
            )
        if error := self._geometry_error():
            raise ValueError(error)
        node = getattr(self, "_discharge_node", None)
        if self.gas_baseline == "flat_plate" and not (
            isinstance(node, MomentumChamberNode) and node.area > 0.0
        ):
            raise ValueError(
                f"EffusionPlateElement {self.id!r}: gas_baseline='flat_plate' needs a "
                "flowing gas side, a MomentumChamberNode with an area, to discharge "
                f"into; {self.to_node!r} has no gas velocity (a plenum takes the "
                "imposed gas_heat_flux, and no baseline applies)."
            )
        self._hole_cd.validate_exit_crossflow()
        if self.porosity >= 1.0:
            raise ValueError(
                f"EffusionPlateElement {self.id!r}: porosity is "
                f"{self.porosity:.3f}; the holes do not fit in the panel."
            )
        if (
            self.crossflow_segments is not None
            and self.correlation != "McGreehanSchotsch"
            and self.hole_cd.inlet_crossflow != "none"
        ):
            raise ValueError(
                f"EffusionPlateElement {self.id!r}: a duct-fed panel needs the "
                "supply-side crossflow term, which only 'McGreehanSchotsch' has "
                f"({self.correlation!r} would silently treat the duct as a plenum)."
            )

    def _geometry_error(self) -> str:
        """Why the holes cannot be drilled at this pitch (#1024), or "".

        Holes must leave material between them. The panel's arrangement is
        not given, so its rows are taken inline: the spanwise pitch must
        exceed d, the streamwise one the hole's footprint on the surface,
        d / sin(angle). Equal is refused too: touching holes leave no plate."""
        return _hole_spacing_error(
            f"EffusionPlateElement {self.id!r}",
            self.hole_diameter,
            self.angle_deg,
            ("pitch_y", self.pitch_y),
            ("pitch_x", self.pitch_x),
        )

    def internal_heat_transfer(self, state_in: NetworkMixtureState) -> dict[str, float]:
        """Coolant-side heat transfer for this panel, Andrews 86-GT-225.

        The hole approach and the throat, summed (Eq. 19). Returns the two
        Nusselt numbers, the Reynolds number they share, and the coefficient
        expressed BOTH ways -- because which area an `h` belongs to is
        exactly the kind of thing that goes wrong silently:

          `h_hole_area`  : on the hole internal surface, pi d L_hole
          `h_plate_area` : on the approach area per hole, cell - hole, which
                           is Andrews' own definition and what his Fig. 8
                           plots

        The two differ by `approach_to_hole_area_ratio` -- 3.46 for Andrews' plate C -- so
        using one where the other is meant is a factor-of-three error, not a
        refinement. For a panel energy balance use
        `Q = h_plate_area * area_approach_total * dT`.

        Properties are taken at the COOLANT inlet state. The coolant heats up
        through the hole, so a large temperature rise makes this a first
        approximation; Andrews' own rig was near-ambient.

        NON-SQUARE ARRAYS. Andrews' plates are square and Eq. (18) carries
        the pitch explicitly, so a rectangular array is represented by its
        equivalent square pitch, `pitch_actual`. Fine while the aspect ratio
        is near one; a strongly rectangular array is outside what the
        correlation was fitted on.
        """
        ts = cb.transport_state(state_in.Tt, state_in.P, state_in.X)
        mu = ts.mu if ts.mu > 0.0 else 1.8e-5
        k = ts.k
        pr = ts.Pr

        # No coolant flow below the network's no-through-flow floor (#792,
        # #998): a zero-flow root ends with |m_dot| ~ 1e-31, not 0, and the
        # coefficient there (Re ~ 1e-25) put the cold face at 5.7e11 K.
        if not abs(state_in.m_dot) >= cb.MIXER_HEAT_MDOT_FLOOR:
            return {}
        m_hole = abs(state_in.m_dot) / self.n_holes
        Re = 4.0 * m_hole / (math.pi * self.hole_diameter * mu)

        X = self.pitch_actual
        nu_approach = cb.effusion_approach_nusselt(Re, pr, X / self.hole_length)
        nu_throat = cb.effusion_throat_nusselt(Re, pr, self.hole_length / self.hole_diameter)

        area_hole = math.pi * self.hole_diameter * self.hole_length
        area_approach = self.panel_area / self.n_holes - math.pi * (self.hole_diameter**2) / 4.0

        # `internal_Nu_multiplier` is the user's rig-matching knob and is 1.0
        # unless they set it. Applied to the SUM, so the split between the
        # approach and throat terms stays the correlation's.
        nu_internal = (nu_approach + nu_throat) * self.internal_Nu_multiplier
        h_hole = nu_internal * k / self.hole_diameter
        return {
            "Re_hole": float(Re),
            "Nu_approach": float(nu_approach),
            "Nu_throat": float(nu_throat),
            "Nu_internal": float(nu_internal),
            "internal_Nu_multiplier": float(self.internal_Nu_multiplier),
            "h_hole_area": float(h_hole),
            "h_plate_area": float(h_hole * area_hole / area_approach),
            "approach_to_hole_area_ratio": float(area_approach / area_hole),
            "area_approach_total": float(area_approach * self.n_holes),
            "area_hole_total": float(area_hole * self.n_holes),
        }

    def overall_effectiveness(
        self,
        state_in: NetworkMixtureState,
        h_gas_unblown: float,
        T_gas: float,
        U_gas: float | None = None,
        gas_augmentation: float = 1.0,
        eta_film: float = 0.0,
    ) -> dict[str, float]:
        """Overall cooling effectiveness of this panel, as an OUTPUT.

            h_gas = h_gas_unblown * gas_augmentation
            eta   = (h_i + h_gas eta_film) / (h_i + h_gas)
            T_wall = T_gas - eta (T_gas - T_coolant)

        The gas side is the caller's; `eta` and the wall temperature are
        the outputs. An overall-effectiveness CORRELATION must never be the
        closure here -- it would already contain the internal convection
        this computes, and the two would double-count.

        THE FILM IS TWO NUMBERS, NOT ONE, and the signature says so. The
        standard film-cooling form is `q = h_f (T_aw - T_w)` with
        `T_aw = T_gas - eta_film (T_gas - T_c)`, so a film owes you both a
        driving temperature (`eta_film`) AND a conductance
        (`gas_augmentation = h_f/h_0`). An adiabatic effectiveness
        correlation such as `film_effectiveness_baldauf_2002` supplies only
        the first: an adiabatic wall passes no heat, so the experiment
        fixes `T_aw` and measures no coefficient. Supplying `eta_film`
        while leaving `gas_augmentation` at 1.0 therefore OVER-PREDICTS,
        and both default to their no-film values so that omitting the pair
        is consistent rather than half-right.

        `gas_augmentation` IS THE CALLER'S, NEVER FITTED HERE. It matches
        `Nu_multiplier` on `ConvectiveSurface`: the library states what the
        correlations give, and matching a specific rig is the user's job.
        See docs/VALIDATION_POLICY.md.

        MIND WHAT `h_gas_unblown` IS MEASURED AGAINST. An augmentation
        ratio is meaningless without its baseline, and mixing baselines is
        a factor-level error, not a refinement. For Andrews 88-GT-290's rig
        the required `gas_augmentation` (plate B, with Baldauf's film) is
        2.7-4.7 against a fully developed Dittus-Boelter and 1.5-2.7
        against the same duct with a thermal-entry correction -- a factor
        of 1.76 from that choice alone (both plates, both baselines:
        validation/cooling/extractions/andrews_effusion_overall_eta.md). Published film-cooling ratios are usually referenced to a
        flat-plate turbulent boundary layer at the same x, which is a third
        baseline again -- the element's `gas_baseline='flat_plate'` in a
        network. On Behrendt 2008's reacting rig it halves the chamber
        baseline's over-prediction (+21% to +12%, extractions/
        behrendt_2008_effusion_eta_tot.md) and does not close it.

        WHAT SCORING IT AGAINST ANDREWS SHOWED. At `gas_augmentation = 1.0`
        and `eta_film = 0.0` the closure scores +3.1% on his effusion plate
        C and +23.7% on plate B -- and plate C's agreement is a
        CANCELLATION, a real film raising `eta` against its augmentation
        lowering it, not evidence that either is absent. The two plates
        need gas-side coefficients differing by 1.7x, and that ratio
        survives any choice of film model (scaling Baldauf's `eta_film`
        from 0 to 1.25x moves it only 1.96 to 1.60). So the missing physics
        is a `gas_augmentation` correlation for full-coverage effusion,
        which needs a HEATED-wall measurement; a better film correlation
        cannot supply it. See
        `validation/cooling/extractions/andrews_effusion_overall_eta.md`.

        Pass `U_gas` to get the jet ratios reported beside the result.
        They are REPORTED AND NEVER APPLIED -- no threshold is offered,
        because the requirement collapses on none of the velocity ratio,
        the blowing ratio or the momentum flux ratio. If `velocity_ratio`
        exceeds about 1 -- Andrews' plate B ejects at 53 m/s into a
        26.8 m/s crossflow -- treat `eta` as an upper bound unless
        `gas_augmentation` already accounts for it.

        Returns `{}` when there is no coolant flow or no gas side. With
        `U_gas` omitted the jet ratios are absent rather than guessed.
        """
        internal = self.internal_heat_transfer(state_in)
        if not internal or h_gas_unblown <= 0.0 or gas_augmentation <= 0.0:
            return {}
        if not 0.0 <= eta_film < 1.0:
            raise ValueError(
                "EffusionPlateElement: eta_film must be in [0, 1); it is an "
                "adiabatic effectiveness, not an overall one"
            )

        h_i = internal["h_plate_area"]
        h_gas = h_gas_unblown * gas_augmentation
        T_c = state_in.Tt
        # The two-resistance closure eta = (h_i + h_gas eta_film)/(h_i + h_gas)
        # is the wall core below with no conduction resistance.
        T_aw = T_gas - eta_film * (T_gas - T_c)
        t_wall, _, _ = _effusion_wall(h_i, T_c, h_gas, T_aw, 0.0)
        eta = (T_gas - t_wall) / (T_gas - T_c) if T_gas != T_c else 0.0
        out = {
            "eta_overall": float(eta),
            "T_wall": float(t_wall),
            "T_adiabatic_wall": float(T_gas - eta_film * (T_gas - T_c)),
            "h_internal_plate_area": float(h_i),
            "internal_Nu_multiplier": float(self.internal_Nu_multiplier),
            "h_gas_unblown": float(h_gas_unblown),
            "gas_augmentation": float(gas_augmentation),
            "h_gas": float(h_gas),
            "eta_film": float(eta_film),
            "resistance_ratio": float(h_i / h_gas),
            "q_flux": float(h_gas * (T_gas - eta_film * (T_gas - T_c) - t_wall)),
        }
        if U_gas is None or U_gas <= 0.0:
            return out

        rho_c = cb.density(T_c, state_in.P, state_in.X)
        area_hole = math.pi * (self.hole_diameter**2) / 4.0
        mass_flux_c = abs(state_in.m_dot) / self.n_holes / area_hole
        u_jet = mass_flux_c / rho_c
        rho_g = rho_c * T_c / T_gas  # same static pressure, ideal gas
        blowing = mass_flux_c / (rho_g * U_gas)
        out.update(
            {
                "u_jet": float(u_jet),
                "velocity_ratio": float(u_jet / U_gas),
                "blowing_ratio": float(blowing),
                "momentum_flux_ratio": float(blowing * blowing / (rho_c / rho_g)),
                "density_ratio": float(rho_c / rho_g),
            }
        )
        return out

    def film_provider_error(self) -> str | None:
        """Why this panel cannot cover a downstream wall, or None (#386)."""
        if self.panel_length is None:
            return (
                "an effusion panel's film needs panel_length, to count the hole rows it builds over"
            )
        return None

    def film_span_average(
        self, x_start: float, x_end: float, M: float, density_ratio: float
    ) -> Any:
        """The panel's film over [x_start, x_end] [m] downstream of its last row.

        A FILM PROVIDER like FilmCoolingElement (#386): Baldauf (2002) per row,
        the rows ``pitch_x`` apart, Sellers over them -- the same pieces as the
        panel's own ``gas_film='baldauf_sellers'``, read past the panel instead
        of averaged over it. Naming the panel as a wall's ``film`` selects it,
        whatever ``gas_film`` says for the panel's own wall. Sellers
        overestimates as rows accumulate, and a full-coverage panel is usually
        outside Baldauf's s/D and M envelope: ``extrapolated`` says so.
        """
        error = self.film_provider_error()
        if error is not None:
            raise ValueError(f"EffusionPlateElement {self.id!r}: {error}")
        D = self.hole_diameter
        return cb.film_span_average_baldauf_2002(
            max(1, round(self.panel_length / self.pitch_x)),
            self.pitch_x / D,
            x_start / D,
            x_end / D,
            M,
            density_ratio,
            self.angle_deg,
            self.pitch_y / D,
            self.turbulence_intensity,
        )

    def wall_heat_transfer(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float]:
        """The plate's own wall, with the gas side its discharge node decides.

        The plate IS the wall the coolant flows through, so it owns it: no
        ThermalWall connects to it. Its heat goes from the gas into the
        effusing coolant, which discharges back into the same node, so the
        network's energy balance is unchanged; the wall temperatures and the
        heat flux are outputs.

        COOLANT SIDE: Andrews 86-GT-225 on the plate area (`internal_heat_
        transfer`), at the supply temperature.

        GAS SIDE, by discharge node:

        * a MomentumChamberNode has flow. It gives the gas temperature and
          its velocity the blowing ratio. The UNBLOWN coefficient is
          `gas_baseline`'s: the chamber's own surface correlation
          ('chamber', default), or ('flat_plate') CR-2648's turbulent
          flat-plate St_o averaged over x in [gas_start_length,
          gas_start_length + panel_length] from where the boundary layer
          restarts, at the chamber's static state and velocity
          (`h_gas_unblown_x_start`/`_x_end`/`_Re_x_start`/`_Re_x_end`
          report the range; `h_gas_unblown_extrapolated` is 1 where its Re_x
          leaves the Stanford data's 1.5e5-5e6, evaluated anyway).
          `gas_augmentation` (the caller's, 1.0 unless set;
          never fitted -- see docs/VALIDATION_POLICY.md) scales the
          coefficient; Andrews 88-GT-290 shows that augmentation, not a
          film correlation, is the missing physics. `gas_film` stays 'none'
          by default for the same reason: Baldauf + Sellers
          ('baldauf_sellers') is selectable and reported with its envelope
          flag, but the data refuse it offered alone.
        * anything else (a plenum) has a state but no flow, so no
          coefficient: the gas side is the imposed `gas_heat_flux` [W/m^2],
          0 by default (adiabatic).

        RADIATION (#1071, `radiation`, opt-in, a gas-side chamber only): Eq.
        8.7 from the gas onto the hot face, at the gas the convection sees --
        T_gas, the chamber's static P and that stream's composition, each
        unless the GasRadiation fixes it. The film lowers the convective
        driving temperature only; the gas radiates at T_gas, as on a
        ThermalWall face. `q_rad_hot` [W/m^2] and `eps_g` are reported (NaN
        where no gas side is defined, 0 without coolant: the adiabatic wall
        then sits at T_gas). On a plenum the imposed `gas_heat_flux` stays
        the whole gas side; radiation is not added to it (both keys NaN).

        WALL: one layer, `wall_thickness / wall_conductivity`, through C++'s
        wall_coupling_and_jacobian; with radiation, C++'s
        wall_radiation_sides (the surface temperatures solved with the T^4
        term, no outer iteration).

        WITHOUT COOLANT FLOW the keys stay. On a gas-side chamber the plate is
        then adiabatic (no coolant takes heat): q = 0, both faces at the gas
        side's adiabatic wall temperature, no film. On a plenum, or a chamber
        with no gas flow, the wall temperatures are undefined (NaN): nothing
        sets them, and `no_cooling_path` is 1 (the solver refuses the root).

        LIMIT (#998): `T_wall_limit` is the plate material's envelope hot end
        (`cb.GENERIC_METAL_SOLIDUS` without a material); `wall_above_limit`
        is 1 when the hotter face is above it, 0 when not, NaN when the wall
        temperature is undefined.
        """
        internal = self.internal_heat_transfer(state_in)
        no_coolant = not internal
        h_i = math.nan if no_coolant else internal["h_plate_area"]
        T_c = float(state_in.Tt)
        A = self.panel_area
        t_over_k = self.wall_thickness / self.wall_conductivity
        out: dict[str, float] = {
            "h_internal_plate_area": float(h_i),
            "wall_conductivity": float(self.wall_conductivity),
            "coolant_crossflow_ignored": float(getattr(self, "_coolant_crossflow_ignored", False)),
            # Duct-fed: the Cd sees the crossflow, but the coolant-side heat
            # transfer is still Andrews' still-plenum correlation (provisional
            # until a crossflow-supply source; #785).
            "coolant_ht_plenum_assumed": float(self.crossflow_segments is not None),
        }
        node = getattr(self, "_discharge_node", None)
        gas = _chamber_gas_side(node, state_out)
        if gas is not None and self.gas_baseline == "flat_plate":
            gas = self._flat_plate_gas_side(node, state_out, gas)
        # NO COOLING PATH (#998): no coolant through the holes and no gas
        # flow on the other side. Nothing takes heat from the plate, so its
        # temperature is undetermined: an imposed flux q != 0 has no steady
        # state, and q = 0 sets no T. A real plate here heats until it fails,
        # and the solver refuses the root, naming the plate. Not modelled,
        # deliberately (no use case: an uncooled plate in a combustor does not
        # survive). Where to start if one comes:
        # * RADIATION: an emissivity and a surroundings T on each face. The
        #   smallest term that always gives a steady T, since it grows as
        #   T^4 (a ThermalWall face already carries GasRadiation, #793).
        # * BACKSIDE CONVECTION on the coolant face without hole through-flow:
        #   a user h, or a sourced plenum-circulation or natural-convection
        #   correlation (register first: scripts/sources.py find).
        # * CONDUCTION to the mounts, or laterally within the plate. A plate
        #   starved over part of its area, or parallel patches of one plate
        #   with some near zero flow, need the lateral path to the cooled
        #   region: that is #898's wall-to-wall conduction (coarse 2-D walls
        #   from 1-D ThermalWalls), applied to effusion plates, rather than a
        #   plate-specific no-flow model.
        out["no_cooling_path"] = float(no_coolant and gas is None)
        if gas is None and isinstance(node, MomentumChamberNode) and node.area > 0.0:
            # A gas-side chamber with no gas flow: no coefficient, so nothing
            # on the gas side is defined. Its keys stay (the topology, not the
            # flow, decides which side this is), NaN.
            out["gas_side_chamber"] = 1.0
            for key in _EFFUSION_CHAMBER_SIDE_KEYS:
                out[key] = math.nan
            if self.gas_baseline == "flat_plate":
                for key in _EFFUSION_FLAT_PLATE_KEYS:
                    out[key] = math.nan
            if self.radiation is not None:
                out["q_rad_hot"] = out["eps_g"] = math.nan
            return self._with_wall_limit(out)
        if gas is None:
            q = self.gas_heat_flux
            t_cold = math.nan if no_coolant else T_c + q / h_i
            out.update(
                {
                    "gas_side_chamber": 0.0,
                    "gas_heat_flux_imposed": float(q),
                    "q_wall": float(q),
                    "Q_wall": float(q * A),
                    "T_wall_cold": float(t_cold),
                    "T_wall_hot": float(t_cold + q * t_over_k),
                }
            )
            if self.radiation is not None:
                out["q_rad_hot"] = out["eps_g"] = math.nan
            return self._with_wall_limit(out)

        rho_c, _ = _safe_rho(cb.density(T_c, state_out.P, state_in.X))
        U_g, rho_g, T_g = gas["U_gas"], gas["rho_gas"], gas["T_gas"]
        hole_area = math.pi * self.hole_diameter**2 / 4.0
        G_jet = abs(state_in.m_dot) / (self.n_holes * hole_area)
        blowing = G_jet / (rho_g * U_g) if U_g > 0.0 else math.inf
        eta_film, film_extrapolated = 0.0, False
        if self.gas_film == "baldauf_sellers" and math.isfinite(blowing) and not no_coolant:
            film = cb.effusion_panel_film_effectiveness(
                max(1, round(self.panel_length / self.pitch_x)),
                self.pitch_x / self.hole_diameter,
                self.pitch_y / self.hole_diameter,
                blowing,
                rho_c / rho_g,
                self.angle_deg,
                self.turbulence_intensity,
            )
            eta_film, film_extrapolated = float(film.eta), bool(film.extrapolated)
        T_aw_g = T_g - eta_film * (T_g - T_c)
        h_g = gas["h_unblown"] * self.gas_augmentation
        if no_coolant:
            # HOT GAS, NO COOLANT (#998): the adiabatic limit is the honest
            # answer -- nothing takes the heat, so both faces sit at the gas
            # side's T_aw. In a combustor that is above any alloy's solidus:
            # wall_above_limit flags it (_with_wall_limit). The same three
            # terms as above (radiation, backside convection, conduction)
            # would pull it below T_aw; none is modelled.
            t_hot, t_cold, q = T_aw_g, T_aw_g, 0.0  # adiabatic: nothing takes the heat
        elif self.radiation is None:
            t_hot, t_cold, q = _effusion_wall(h_i, T_c, h_g, T_aw_g, t_over_k)
        if self.radiation is not None:
            face = _radiation_face(self.radiation, (T_g, float(state_out.P), gas["Y_gas"]), 1.0)
            out["eps_g"] = face.eps_g
            if no_coolant:
                # No film without coolant, so T_aw_g is T_gas: the gas
                # radiates onto a face at its own temperature, q_rad 0.
                out["q_rad_hot"] = 0.0
            else:
                w = _solver_tools.wall_radiation_sides(
                    T_aw_g, h_g, T_c, h_i, t_over_k, face.kind, face.T_src, face.k, 0, 0.0, 0.0
                )
                t_hot, t_cold, q = float(w.T_s_a), float(w.T_s_b), float(w.q)
                out["q_rad_hot"] = float(w.q_rad_a)
        out.update(
            {
                "gas_side_chamber": 1.0,
                "h_gas_unblown": float(gas["h_unblown"]),
                "h_gas_unblown_model": gas["baseline_model"],
                "h_gas_unblown_Re": gas["Re_unblown"],
                "h_gas_unblown_Pr": gas["Pr_unblown"],
                "h_gas_unblown_Nu": gas["Nu_unblown"],
                "gas_augmentation": float(self.gas_augmentation),
                "h_gas": float(h_g),
                "T_gas": float(T_g),
                "T_gas_from_main_inlet": float(gas["T_from_main_inlet"]),
                "eta_film": eta_film,
                "film_extrapolated": float(film_extrapolated),
                "T_adiabatic_wall": float(T_aw_g),
                "q_wall": float(q),
                "Q_wall": float(q * A),
                "T_wall_hot": float(t_hot),
                "T_wall_cold": float(t_cold),
                "eta_overall": float((T_g - t_hot) / (T_g - T_c)) if T_g != T_c else 0.0,
                "U_gas": float(U_g),
                "blowing_ratio": float(blowing),
                "density_ratio": float(rho_c / rho_g),
                "velocity_ratio": float(G_jet / rho_c / U_g) if U_g > 0.0 else math.inf,
            }
        )
        if self.gas_baseline == "flat_plate":
            out.update({k: gas[k] for k in _EFFUSION_FLAT_PLATE_KEYS})
        return self._with_wall_limit(out)

    def _flat_plate_gas_side(
        self, node: "MomentumChamberNode", state: NetworkMixtureState, gas: dict[str, Any]
    ) -> dict[str, Any]:
        """The chamber's gas side with the unblown coefficient of a turbulent
        flat plate restarted `gas_start_length` upstream of the panel.

        The chamber still supplies the stream -- T_gas, U_gas, rho -- and only
        the coefficient's correlation changes: CR-2648's St_o averaged over
        x in [x0, x0 + panel_length] (C++ `flat_plate_turbulent_panel_
        nusselt`), properties at the chamber's static state like its own
        surface correlation (#578)."""
        L = float(self.panel_length)
        x0 = float(self.gas_start_length)
        T_static = float(node._closure(state, node._total_m_dot).T_static)
        ts = cb.transport_state(T_static, state.P, state.X)
        re_L = gas["rho_gas"] * gas["U_gas"] * L / ts.mu
        fp = cb.flat_plate_turbulent_panel_nusselt(re_L, ts.Pr, x0 / L)
        return {
            **gas,
            "h_unblown": float(fp.Nu * ts.k / L),
            "baseline_model": _FLAT_PLATE_MODEL,
            "Re_unblown": float(re_L),
            "Pr_unblown": float(ts.Pr),
            "Nu_unblown": float(fp.Nu),
            "h_gas_unblown_x_start": x0,
            "h_gas_unblown_x_end": x0 + L,
            "h_gas_unblown_Re_x_start": float(fp.Re_x_start),
            "h_gas_unblown_Re_x_end": float(fp.Re_x_end),
            "h_gas_unblown_extrapolated": float(fp.extrapolated),
        }

    def _with_wall_limit(self, out: dict[str, float]) -> dict[str, float]:
        """The wall solution plus its temperature limit (#998): the plate
        material's envelope hot end, and whether the hotter face is above it
        (NaN when the wall temperature is undefined). Always present."""
        faces = (float(out["T_wall_hot"]), float(out["T_wall_cold"]))
        out["T_wall_limit"] = self.T_wall_limit
        # Not max(): it drops a NaN in its second argument.
        out["wall_above_limit"] = (
            float(any(t > self.T_wall_limit for t in faces))
            if all(math.isfinite(t) for t in faces)
            else math.nan
        )
        return out

    def diagnostics(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
    ) -> dict[str, float]:
        base = super().diagnostics(state_in, state_out, flows)
        m_dot = abs(state_in.m_dot)
        dP_drive = state_in.Pt - state_out.P
        base.update(self._hole_cd.supply_crossflow_diagnostics(state_in, state_out, flows))
        out = {
            "n_holes": float(self.n_holes),
            "hole_count_exact": float(self.hole_count_exact),
            "porosity": float(self.porosity),
            "hole_length": float(self.hole_length),
            "L_over_d": float(self.hole_length / self.hole_diameter),
            "pitch_over_d": float(self.pitch_actual / self.hole_diameter),
            # Andrews et al. (1988) correlate effusion cooling on the coolant
            # mass flow per unit PLATE area, not per hole. Their measurements
            # span G = 0.1 to 1.6 kg/s/m^2.
            "G_coolant": float(m_dot / self.panel_area),
            "dP_drive": float(dP_drive),
            # A homogenised panel would otherwise average an ingesting hole
            # into a healthy net outflow. Say it instead.
            "is_ingesting": float(dP_drive <= 0.0),
            **self._hole_cd.inclination_diagnostics(),
            **self._hole_cd.exit_crossflow_diagnostics(state_in, state_out, flows),
        }
        # Coolant-side heat transfer, when there is flow to carry it, and the
        # plate's own wall with the gas side its discharge node decides.
        out.update(self.internal_heat_transfer(state_in) or self._internal_without_flow())
        out.update(self.wall_heat_transfer(state_in, state_out))
        out.update(base)
        return out

    def _internal_without_flow(self) -> dict[str, float]:
        """internal_heat_transfer's keys without coolant flow: Re 0, the
        geometry as ever, the Nusselt numbers and h undefined (NaN) -- the
        keys stay, whatever the operating point."""
        area_hole = math.pi * self.hole_diameter * self.hole_length
        area_approach = self.panel_area / self.n_holes - math.pi * (self.hole_diameter**2) / 4.0
        return {
            "Re_hole": 0.0,
            "Nu_approach": math.nan,
            "Nu_throat": math.nan,
            "Nu_internal": math.nan,
            "internal_Nu_multiplier": float(self.internal_Nu_multiplier),
            "h_hole_area": math.nan,
            "h_plate_area": math.nan,
            "approach_to_hole_area_ratio": float(area_approach / area_hole),
            "area_approach_total": float(area_approach * self.n_holes),
            "area_hole_total": float(area_hole * self.n_holes),
        }


class FilmCoolingElement(OrificeElement):
    """Rows of film-cooling holes: coolant from `from_node` into the gas at `to_node`.

    The holes are an orifice (the flow) and a FILM PROVIDER (the film). The
    element owns no wall and is not a wall side: the surface the film cools
    lies downstream along the gas path, so a ThermalWall on the gas duct
    there names it with ``film=`` and its gas-side adiabatic wall
    temperature becomes

        T_aw,f = T_aw,gas - eta_bar (T_aw,gas - T_coolant)

    with eta_bar the film over the wall's span (#386). That wall's h, stream
    and heat booking stay the duct's.

    FILM, by hole type:

    * ``"cylindrical"``: Baldauf (2002) per row (`s/D` from
      ``lateral_pitch``, `Tu` from ``turbulence_intensity``).
    * ``"shaped"``: Colban (2011) Eq. 19 per row. ``hole_diameter`` is the
      METERING diameter; ``area_ratio`` is A_exit/A_inlet and
      ``breakout_width`` the hole's exit width t.

    Rows superpose by Sellers, ``row_pitch`` apart; a wall's span is measured
    from the LAST row. The blowing ratio is the metering mass flux over the
    gas duct's, the density ratio is taken at a common static pressure.
    Neither kernel was fitted on a heated wall downstream of a film row: the
    coupled wall temperature is a composition of published pieces, not a
    scored prediction.

    FLOW: an orifice of n_rows * holes_per_row holes in parallel, each of
    the drilled length ``wall_thickness / sin(angle_deg)``, plenum-fed.
    ``HoleCd(inclination='auto')`` (default) adds Idelchik diagram 3-2's
    inclined-inlet loss to the IdelchikThick Cd of a CYLINDRICAL hole (see
    EffusionPlateElement) and leaves a shaped hole without it: the plain-hole
    Cd on the metering area has no shape term, and the inclination term
    alone moves the shaped holes away from Gritsch 1998's data (ACCURACY,
    CROSS-SOURCE: fan -2 -> -11 %, laidback 0 -> -10 % bias without
    crossflow).
    ``'idelchik_3_2'`` forces it, ``'none'`` switches it off; reported as
    ``Cd_inclination_dzeta``. The Cd's choices are ``hole_cd=HoleCd(...)``
    as on EffusionPlateElement (the per-element keywords are deprecated
    since 0.12.0, removed in 0.13.0); the holes are plenum-fed.
    """

    _selects_compressible_model = False  # keeps the nozzle (#861)
    # The deprecated hole-Cd keywords (HoleCd replaces them): not serialized.
    _deprecated_init_kwargs = ("correlation", "Cd")

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        hole_diameter: float,
        holes_per_row: int,
        lateral_pitch: float,
        wall_thickness: float,
        n_rows: int = 1,
        row_pitch: float | None = None,
        angle_deg: float = 30.0,
        hole_type: Literal["cylindrical", "shaped"] = "cylindrical",
        area_ratio: float | None = None,
        breakout_width: float | None = None,
        turbulence_intensity: float = 0.05,
        correlation: str | None = None,
        Cd: float | None = None,
        edge_radius: float = 0.0,
        exit_crossflow_region: Literal[
            "suction_surface", "pressure_surface", "leading_edge"
        ] = "suction_surface",
        exit_crossflow_orientation_deg: float = 0.0,
        hole_cd: HoleCd | None = None,
    ) -> None:
        name = f"FilmCoolingElement {id!r}"
        # The hole Cd's model choices (hole_cd.py): a HoleCd, or the
        # deprecated per-element keywords; IdelchikThick by default.
        spec = _with_correlation(
            _hole_cd_spec(
                name,
                hole_cd,
                {"correlation": correlation, "Cd": Cd},
                stacklevel=3,
            ),
            "IdelchikThick",
            name,
        )
        if hole_diameter <= 0.0 or lateral_pitch <= 0.0 or wall_thickness <= 0.0:
            raise ValueError(
                f"{name}: hole_diameter, lateral_pitch and wall_thickness must be positive"
            )
        if holes_per_row < 1 or n_rows < 1:
            raise ValueError(f"{name}: holes_per_row and n_rows must be at least 1")
        if n_rows > 1 and not (row_pitch and row_pitch > 0.0):
            raise ValueError(f"{name}: more than one row needs a positive row_pitch")
        if not 0.0 < angle_deg <= 90.0:
            raise ValueError(f"{name}: angle_deg is the hole angle to the SURFACE, in (0, 90]")
        if hole_type not in ("cylindrical", "shaped"):
            raise ValueError(f"{name}: hole_type is 'cylindrical' or 'shaped'")
        if hole_type == "shaped":
            if area_ratio is None or breakout_width is None:
                raise ValueError(
                    f"{name}: shaped holes need area_ratio (A_exit/A_inlet) and "
                    "breakout_width (Colban's t)"
                )
            if area_ratio < 1.0 or not 0.0 < breakout_width <= lateral_pitch:
                raise ValueError(
                    f"{name}: area_ratio must be >= 1 and breakout_width in (0, lateral_pitch]"
                )
        if turbulence_intensity <= 0.0:
            raise ValueError(f"{name}: turbulence_intensity must be positive")

        self.hole_diameter = float(hole_diameter)
        self.holes_per_row = int(holes_per_row)
        self.n_rows = int(n_rows)
        self.n_holes = self.n_rows * self.holes_per_row
        self.lateral_pitch = float(lateral_pitch)
        self.row_pitch = float(row_pitch) if row_pitch else 0.0
        self.wall_thickness = float(wall_thickness)
        self.angle_deg = float(angle_deg)
        self.hole_type = hole_type
        self.area_ratio = area_ratio
        self.breakout_width = breakout_width
        self.turbulence_intensity = float(turbulence_intensity)
        self.hole_length = self.wall_thickness / math.sin(math.radians(angle_deg))
        hole_area = math.pi * self.hole_diameter**2 / 4.0
        self.total_hole_area = self.n_holes * hole_area

        super().__init__(
            id,
            from_node,
            to_node,
            Cd=spec.Cd,
            area=self.total_hole_area,
            correlation=spec.correlation,
            plate_thickness=self.hole_length,
            edge_radius=edge_radius,
        )
        self.equivalent_bore = self.diameter
        # Into a merge chamber the jets enter at the hole angle to the gas axis.
        self.injection_angle_deg = float(angle_deg)
        # The hole Cd's options (hole_cd.py): the inclined-inlet term
        # (CD_INCLINATION_MODELS, "auto" resolved from the correlation and
        # the hole type) and the opt-in exit-crossflow term
        # (CD_EXIT_CROSSFLOW_MODELS); a shaped hole takes Table 2's
        # fan-shaped sets and the exit-plane flux. Plenum-fed.
        self.exit_crossflow_region = exit_crossflow_region
        self.exit_crossflow_orientation_deg = float(exit_crossflow_orientation_deg)
        self.cd_inclination, self.cd_exit_crossflow = self._hole_cd.configure_cooling_hole(
            name,
            spec,
            hole_type=hole_type,
            angle_deg=angle_deg,
            exit_crossflow_region=exit_crossflow_region,
            exit_crossflow_orientation_deg=self.exit_crossflow_orientation_deg,
            exit_area_ratio=float(area_ratio) if hole_type == "shaped" else 1.0,
        )
        self.hole_cd = dataclasses.replace(spec, Cd=self.Cd)

    def _hole_count(self) -> float:
        return float(self.n_holes)

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """A wall's holes: no pipe, hence no beta; the correlations see one hole."""
        self._orifice_geom = cb.OrificeGeometry()
        self._orifice_geom.d = self.hole_diameter
        self._orifice_geom.D = 0.0
        self._orifice_geom.t = self.hole_length
        self._orifice_geom.r = self.edge_radius
        self.beta = 0.0
        self._hole_cd.resolve_exit_mainstream(graph)
        self._resolve_plate_against_chambers(graph)

    def validate(self) -> None:
        super().validate()
        self._hole_cd.validate_exit_crossflow()
        if self.correlation in _METERING_SELECTORS:
            raise ValueError(
                f"FilmCoolingElement {self.id!r}: {self.correlation!r} is a NORMED "
                "metering correlation for a plate in a pipe; film holes have no "
                "pipe and no beta. Use 'IdelchikThick', 'McGreehanSchotsch' or 'fixed'."
            )
        if error := self._geometry_error():
            raise ValueError(error)

    def _geometry_error(self) -> str:
        """Why the holes cannot be drilled at these pitches (#1024), or "".

        Holes must leave material between them: lateral_pitch > d, and with
        more than one row row_pitch > d / sin(angle), the footprint on the
        surface. The rows' stagger is not given, so they are taken inline (a
        staggered design this refuses gets its own issue). d is the metering
        diameter of a shaped hole."""
        return _hole_spacing_error(
            f"FilmCoolingElement {self.id!r}",
            self.hole_diameter,
            self.angle_deg,
            ("lateral_pitch", self.lateral_pitch),
            ("row_pitch", self.row_pitch) if self.n_rows > 1 else None,
        )

    def film_span_average(
        self, x_start: float, x_end: float, M: float, density_ratio: float
    ) -> Any:
        """The film over [x_start, x_end] [m] downstream of the last row.

        Returns C++'s FilmSpanAverage: eta_bar, d eta_bar/dM, d eta_bar/d
        density ratio, and the kernel's envelope flag.
        """
        D = self.hole_diameter
        pitch = self.row_pitch / D
        if self.hole_type == "shaped":
            return cb.film_span_average_colban_2011(
                self.n_rows,
                pitch,
                x_start / D,
                x_end / D,
                M,
                float(self.area_ratio),
                self.lateral_pitch / D,
                float(self.breakout_width) / self.lateral_pitch,
                self.angle_deg,
                density_ratio,
            )
        return cb.film_span_average_baldauf_2002(
            self.n_rows,
            pitch,
            x_start / D,
            x_end / D,
            M,
            density_ratio,
            self.angle_deg,
            self.lateral_pitch / D,
            self.turbulence_intensity,
        )

    def diagnostics(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
    ) -> dict[str, float]:
        base = super().diagnostics(state_in, state_out, flows)
        dP_drive = state_in.Pt - state_out.P
        base.update(
            {
                "n_holes": float(self.n_holes),
                "hole_length": float(self.hole_length),
                "L_over_d": float(self.hole_length / self.hole_diameter),
                "G_jet": float(abs(state_in.m_dot) / self.total_hole_area),
                "dP_drive": float(dP_drive),
                "is_ingesting": float(dP_drive <= 0.0),
                **self._hole_cd.inclination_diagnostics(),
                **self._hole_cd.exit_crossflow_diagnostics(state_in, state_out, flows),
            }
        )
        return base


# The default free-stream turbulence of a slot film, a FRACTION (#714, #773):
# 0.14 = 14%, D-7958's 52% blockage plate (axial Tu 14%), where Eq. (5)'s
# line fits that report's own Table I best (bias +0.017 in eta), inside the
# 10-20% Juhasz and Marek
# (1971, TN D-6360) measured in a combustor. A working example, not a
# universal value: a liner's Tu comes from what is upstream of it.
SLOT_FILM_DEFAULT_TU = 0.14


class SlotFilmElement(OrificeElement):
    """A 2-D film slot (a combustor liner's cooling ring): coolant from
    `from_node` into the gas at `to_node` through a continuous slot of height
    ``slot_height`` across ``span``.

    FILM: Marek and Tacina (1975), NASA TN D-7958, eta = 1/(1 + C_m x/(M S))
    with C_m = 0.25 Tu/100 (Eqs. 4-5, Tu in percent), ``turbulence_intensity``
    the axial free-stream Tu as a FRACTION like every other element's (0.14 =
    14%, as the kernel takes it too) and ``C_m_multiplier`` the user's rig-matching
    factor (never scored). As for FilmCoolingElement, the element owns no
    wall: a ThermalWall on the gas duct downstream names it with ``film=``,
    and that wall's span is measured from the slot lip. The blowing ratio is
    the slot's mass flux over the gas duct's. Every D-7958 test was
    co-flowing film.

    FLOW: an orifice of the slot's area, ``lip_thickness`` long. The
    discharge correlations are for round holes; the slot enters as its
    hydraulic diameter 2 S W/(S + W), and as the number of such holes that
    carries its area, so the hole Reynolds number is the slot's own. A
    measured slot Cd is better given as ``hole_cd=HoleCd(correlation='fixed',
    Cd=...)`` (``correlation``/``Cd`` keywords: deprecated since 0.12.0,
    removed in 0.13.0). A slot has no inclination, exit- or supply-crossflow
    term.

    FILM-SIDE h (#962), opt-in with ``film_htc``. By default ("none") the
    film only lowers the filmed side's driving temperature to T_aw and the
    side keeps its duct's own h. Near a slot the boundary layer is the slot
    fluid's, with its own, higher h, so with ``film_htc`` the filmed wall
    side uses the film's h, the mean over the wall's span [film_start,
    film_end] of h(x) (an element mean, as every ThermalWall side's h), on
    T_aw from the film effectiveness:

    - "seban1960" (recommended): Seban (1960) Fig. 12 on slot-fluid
      properties, St_s = h/(G_s cp_s) = 0.135 Re_s^-0.3 (x/s)^-0.3 near the
      slot (curve B as drawn; Eq. 13 prints 0.10, which contradicts his own
      curve and data) and 0.40 Re_s^-0.3 (x/s)^-0.6 beyond (curve C),
      joined by a C1 blend at their crossing (x/s 37). His data: x/s about
      4-300, m 1.06-8.91, Re_s about 1.4e3-8e3, air into air.
    - "ballal_lefebvre1973": the m > 1.3 wall-jet Nu_x = 0.10 Re_s^0.8
      (x/s)^0.44 (Cranfield SME4 Eq. 6, Lefebvre & Ballal Eq. 8.44), the
      book's own method. It sits 20-67 % above the Cranfield rig's own data
      at low Re (SME4 Fig. 2), so it is a named option for reproducing the
      book, not a recommendation.

    Past the form's far end (x/s 300; 50 for the wall jet) and below its
    hand-back mass-velocity ratio (m 1.0, Seban's flat-plate finding; 1.3
    for the wall jet; a C1 step over 0.1 in m) the side hands back to the
    duct's own h. The wall reports which branch and box state applied
    (``film_h_*``). Matching a specific rig is the user's job; the duct's
    ``Nu_multiplier`` does not act on the film's h.
    """

    _selects_compressible_model = False  # keeps the nozzle (#861)
    # The deprecated hole-Cd keywords (HoleCd replaces them): not serialized.
    _deprecated_init_kwargs = ("correlation", "Cd")

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        slot_height: float,
        span: float,
        lip_thickness: float,
        turbulence_intensity: float = SLOT_FILM_DEFAULT_TU,
        C_m_multiplier: float = 1.0,
        correlation: str | None = None,
        Cd: float | None = None,
        edge_radius: float = 0.0,
        film_htc: Literal["none", "seban1960", "ballal_lefebvre1973"] = "none",
        hole_cd: HoleCd | None = None,
    ) -> None:
        name = f"SlotFilmElement {id!r}"
        # The slot Cd's model choices (hole_cd.py): a HoleCd, or the
        # deprecated per-element keywords; IdelchikThick by default. A slot
        # has no inclination, exit or supply crossflow term.
        spec = _with_correlation(
            _hole_cd_spec(name, hole_cd, {"correlation": correlation, "Cd": Cd}, stacklevel=3),
            "IdelchikThick",
            name,
        )
        if film_htc not in SLOT_FILM_HTC_MODELS:
            raise ValueError(
                f"{name}: film_htc must be one of {sorted(SLOT_FILM_HTC_MODELS)}, not {film_htc!r}"
            )
        self.film_htc = film_htc
        if slot_height <= 0.0 or span <= 0.0 or lip_thickness <= 0.0:
            raise ValueError(f"{name}: slot_height, span and lip_thickness must be positive")
        if slot_height > span:
            raise ValueError(f"{name}: a slot is wider than it is high (span >= slot_height)")
        if not 0.0 <= turbulence_intensity <= 1.0:
            raise ValueError(
                f"{name}: turbulence_intensity is the axial Tu as a FRACTION in [0, 1] "
                f"(0.14 = 14%); got {turbulence_intensity}"
            )
        if C_m_multiplier <= 0.0:
            raise ValueError(f"{name}: C_m_multiplier must be positive")
        self.slot_height = float(slot_height)
        self.span = float(span)
        self.lip_thickness = float(lip_thickness)
        self.turbulence_intensity = float(turbulence_intensity)
        self.C_m_multiplier = float(C_m_multiplier)
        self.total_hole_area = self.slot_height * self.span
        self.hydraulic_diameter = (
            2.0 * self.slot_height * self.span / (self.slot_height + self.span)
        )
        super().__init__(
            id,
            from_node,
            to_node,
            Cd=spec.Cd,
            area=self.total_hole_area,
            correlation=spec.correlation,
            plate_thickness=self.lip_thickness,
            edge_radius=edge_radius,
        )
        self._hole_cd.configure_plain(name, spec)
        self.hole_cd = dataclasses.replace(spec, Cd=self.Cd)
        self.equivalent_bore = self.hydraulic_diameter
        # Into a merge chamber the film enters along the gas axis.
        self.injection_angle_deg = 0.0

    def _hole_count(self) -> float:
        # Holes of the hydraulic diameter that carry the slot's area: the
        # hole Reynolds number 4 m/(pi Dh mu) per hole is the slot's m Dh/(A mu).
        return self.total_hole_area / (math.pi * self.hydraulic_diameter**2 / 4.0)

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """A slot in a wall: no pipe, hence no beta; the correlations see the
        hydraulic diameter."""
        self._orifice_geom = cb.OrificeGeometry()
        self._orifice_geom.d = self.hydraulic_diameter
        self._orifice_geom.D = 0.0
        self._orifice_geom.t = self.lip_thickness
        self._orifice_geom.r = self.edge_radius
        self.beta = 0.0
        self._resolve_plate_against_chambers(graph)

    def validate(self) -> None:
        super().validate()
        if self.correlation in _METERING_SELECTORS:
            raise ValueError(
                f"SlotFilmElement {self.id!r}: {self.correlation!r} is a NORMED "
                "metering correlation for a plate in a pipe; a film slot has no "
                "pipe and no beta. Use 'IdelchikThick', 'McGreehanSchotsch' or 'fixed'."
            )

    def film_span_average(
        self, x_start: float, x_end: float, M: float, density_ratio: float
    ) -> Any:
        """The film over [x_start, x_end] [m] downstream of the slot lip:
        C++'s FilmSpanAverage (eta_bar, d eta_bar/dM, 0 for the density ratio,
        which Eq. (4) does not contain, and the kernel's envelope flag)."""
        S = self.slot_height
        return cb.slot_film_span_average_marek_tacina_1975(
            x_start / S,
            x_end / S,
            M,
            self.turbulence_intensity,
            self.C_m_multiplier,
        )

    def film_side_htc(
        self,
        x_start: float,
        x_end: float,
        G_s: float,
        M: float,
        h_duct: float,
        T_s: float,
        P_s: float,
        X_s: list[float],
    ) -> Any:
        """The filmed side's h (#962) over [x_start, x_end] [m] downstream of
        the slot lip, or None with ``film_htc="none"``: C++'s SlotFilmHtc
        (h, dh/dG_s, dh/dT_s, dh/dM, dh/dh_duct and the box state). G_s is
        the slot's mass flux, M = G_s / G_gas, T_s, P_s, X_s the slot fluid's
        state (its supply plenum's)."""
        if self.film_htc == "none":
            return None
        kernel = (
            cb.slot_film_htc_seban_1960
            if self.film_htc == "seban1960"
            else cb.slot_film_htc_ballal_lefebvre_1973
        )
        return kernel(x_start, x_end, self.slot_height, G_s, M, h_duct, T_s, P_s, X_s)

    def diagnostics(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
    ) -> dict[str, float]:
        base = super().diagnostics(state_in, state_out, flows)
        dP_drive = state_in.Pt - state_out.P
        base.update(
            {
                "slot_height": float(self.slot_height),
                "hydraulic_diameter": float(self.hydraulic_diameter),
                "G_slot": float(abs(state_in.m_dot) / self.total_hole_area),
                "turbulence_intensity": float(self.turbulence_intensity),
                "dP_drive": float(dP_drive),
                "is_ingesting": float(dP_drive <= 0.0),
            }
        )
        return base


# SlotFilmElement.film_htc choices (#962): the duct's own h, Seban 1960
# (recommended), Ballal & Lefebvre 1973 wall jet (the book's method).
SLOT_FILM_HTC_MODELS = frozenset({"none", "seban1960", "ballal_lefebvre1973"})


class FilmSource(NamedTuple):
    """One coolant stream of a film provider (#386, #695): the element that
    carries it, the node that supplies it, and its metering hole area [m^2]."""

    element_id: str
    supply_node: str
    hole_area: float


@dataclass
class _FilmGroups:
    """A single-source provider's FilmSpanAverage in the grouped layout of
    cb.FilmSpanAverageGroups: one group whose weight is eta."""

    eta: float
    weight: list[float]
    deta_dM: list[float]
    deta_dP: list[float]
    d_weight_dM: list[list[float]]
    d_weight_dP: list[list[float]]
    extrapolated: bool


def film_sources(provider: Any) -> list[FilmSource]:
    """A film provider's coolant streams: its own ``film_sources()``, or the
    single element it is (its flow, supply node and hole area)."""
    own = getattr(provider, "film_sources", None)
    if own is not None:
        return list(own())
    return [FilmSource(provider.id, provider.from_node, float(provider.total_hole_area))]


def film_discharge_node(provider: Any) -> str | None:
    """The one node every source of the provider discharges into, or None."""
    own = getattr(provider, "film_discharge_node", None)
    if own is not None:
        return own()
    return provider.to_node


def film_groups(
    provider: Any, x_start: float, x_end: float, M: list[float], density_ratio: list[float]
) -> Any:
    """The provider's film over the span, one group per source: its own
    ``film_span_average_groups``, or its ``film_span_average`` as one group."""
    own = getattr(provider, "film_span_average_groups", None)
    if own is not None:
        return own(x_start, x_end, M, density_ratio)
    f = provider.film_span_average(x_start, x_end, M[0], density_ratio[0])
    eta, dM, dP = float(f.eta), float(f.deta_dM), float(f.deta_dP)
    return _FilmGroups(eta, [eta], [dM], [dP], [[dM]], [[dP]], bool(f.extrapolated))


class _FilmedHtc:
    """A wall side's surface result under a film (#386): T_aw replaced by

        T_aw,f = T_ref (1 - eta) + sum_g w_g T_c,g,   sum_g w_g = eta,

    one term per coolant source g (an effusion liner's panels, #695; a single
    source gives T_ref - eta (T_ref - T_c)). eta and w_g depend on every
    source's M_g = G_jet,g / G_gas and DR_g = (W_c,g / W_g)(T_gs / T_c,g),
    T_gs the gas's STATIC T at G_gas (the duct closure at (T_g, Pt, G_gas),
    #782 audit): it moves with T_g, the side's Pt and its flow
    (``dlnTs`` = d ln T_gs / d(T_g, Pt, m_gas)).

    T_ref is the side's own T_aw, shifted by ``T_shift`` = T_inf - T_mix to
    the undiluted gas when the solver has it (#694): eta is defined against
    the gas before the coolant mixed in. Partials are chained onto the side's
    own T, flow, P and Pt; each source's flow and coolant temperature come
    back as the lists ``dT_aw_dfilm_mdot`` and ``dT_aw_dT_coolant``, and the
    undiluted gas temperature as ``dT_aw_dT_inf``, for the solver to relay.
    h and the rest are the side's.
    """

    __slots__ = (
        "_r",
        "T_aw",
        "dT_aw_dT",
        "dT_aw_dmdot",
        "dT_aw_dP",
        "dT_aw_dPt",
        "dT_aw_dfilm_mdot",
        "dT_aw_dT_coolant",
        "dT_aw_dT_inf",
        "film",
        "_k",
    )

    def __init__(
        self,
        r: Any,
        groups: Any,
        M: list[float],
        DR: list[float],
        T_c: list[float],
        m_c: list[float],
        m_gas: float,
        dlnTs: tuple[float, float, float],
        film: dict[str, float],
        T_shift: float | None = None,
    ) -> None:
        self._r = r
        # With the undiluted gas, T_ref no longer moves with the side's own
        # (mixed) T: its d/dT_mix of 1 is taken off and goes to T_inf.
        undiluted = T_shift is not None
        T_ref = float(r.T_aw) + (T_shift if undiluted else 0.0)
        eta = float(groups.eta)
        w = [float(v) for v in groups.weight]
        n = len(w)
        k = 1.0 - eta
        self.T_aw = T_ref * k + sum(w[g] * T_c[g] for g in range(n))
        # dT_aw / d(M_h), d(DR_h) at fixed T_ref and T_c.
        a_M = [
            -T_ref * float(groups.deta_dM[h])
            + sum(T_c[g] * float(groups.d_weight_dM[g][h]) for g in range(n))
            for h in range(n)
        ]
        a_P = [
            -T_ref * float(groups.deta_dP[h])
            + sum(T_c[g] * float(groups.d_weight_dP[g][h]) for g in range(n))
            for h in range(n)
        ]
        # The gas temperature in DR: the side's own T, or T_inf when undiluted,
        # through its static T at G_gas, which also moves with Pt and the flow.
        a_DR = sum(a_P[h] * DR[h] for h in range(n))
        dlnTs_dT, dlnTs_dPt, dlnTs_dm = dlnTs
        d_gas_T = a_DR * dlnTs_dT
        self.dT_aw_dT = k * (float(r.dT_aw_dT) - float(undiluted)) + (0.0 if undiluted else d_gas_T)
        self.dT_aw_dT_inf = (k + d_gas_T) if undiluted else 0.0
        self.dT_aw_dmdot = (
            k * float(r.dT_aw_dmdot)
            - sum(a_M[h] * M[h] / m_gas for h in range(n))
            + a_DR * dlnTs_dm
        )
        self.dT_aw_dP = k * float(getattr(r, "dT_aw_dP", 0.0))
        self.dT_aw_dPt = k * float(getattr(r, "dT_aw_dPt", 0.0)) + a_DR * dlnTs_dPt
        self.dT_aw_dfilm_mdot = [a_M[h] * M[h] / m_c[h] for h in range(n)]
        self.dT_aw_dT_coolant = [w[h] - a_P[h] * DR[h] / T_c[h] for h in range(n)]
        self.film = film
        self._k = k

    @property
    def dT_aw_dY(self) -> list[float]:
        """The side's T_aw composition gradient through T_ref (#717). The
        density ratio's dependence on the gas molar mass is not chained."""
        return [self._k * float(v) for v in (getattr(self._r, "dT_aw_dY", None) or [])]

    def __getattr__(self, name: str) -> Any:
        return getattr(self._r, name)

    def __dir__(self) -> list[str]:
        return sorted(set(dir(self._r)) | set(self.__slots__))


class _FilmSideHtc:
    """A filmed wall side with the film's own h (#962): the slot kernel's

        h = h_duct + w(M) f_x (h_film - h_duct),

    h_film from the slot flow (G_s = m_c / A_slot) and the slot fluid's T,
    M = G_s / G_gas. The duct's h partials scale by dh/dh_duct = 1 - w f_x;
    the gas flow also enters through M (dM/dm_gas = -M / m_gas). Each
    source's flow and coolant temperature come back as ``dh_dfilm_mdot`` and
    ``dh_dT_coolant`` for the solver to relay, as the film's T_aw does. The
    slot fluid's composition is not chained (as the density ratio's is not).
    T_aw and the rest are the filmed side's."""

    __slots__ = (
        "_r",
        "_kd",
        "h",
        "dh_dT",
        "dh_dmdot",
        "dh_dP",
        "dh_dPt",
        "dh_dfilm_mdot",
        "dh_dT_coolant",
    )

    def __init__(
        self, r: Any, res: Any, M: float, m_c: float, m_gas: float, hole_area: float
    ) -> None:
        self._r = r
        kd = float(res.dh_dh_duct)
        self._kd = kd
        self.h = float(res.h)
        self.dh_dT = kd * float(r.dh_dT)
        self.dh_dmdot = kd * float(r.dh_dmdot) - float(res.dh_dM) * M / m_gas
        self.dh_dP = kd * float(getattr(r, "dh_dP", 0.0))
        self.dh_dPt = kd * float(getattr(r, "dh_dPt", 0.0))
        self.dh_dfilm_mdot = [float(res.dh_dG_s) / hole_area + float(res.dh_dM) * M / m_c]
        self.dh_dT_coolant = [float(res.dh_dT_s)]

    @property
    def dh_dY(self) -> list[float]:
        """The duct h's composition gradient, scaled by dh/dh_duct."""
        return [self._kd * float(v) for v in (getattr(self._r, "dh_dY", None) or [])]

    @property
    def dh_dsources(self) -> dict[str, float] | None:
        own = getattr(self._r, "dh_dsources", None)
        if not own:
            return own
        return {k: self._kd * float(v) for k, v in own.items()}

    def __getattr__(self, name: str) -> Any:
        return getattr(self._r, name)

    def __dir__(self) -> list[str]:
        return sorted(set(dir(self._r)) | set(self.__slots__))


def film_on_wall_side(
    provider: Any,
    side: Any,
    r: Any,
    gas_state: "NetworkMixtureState",
    coolant: list[tuple["NetworkMixtureState", float]],
    m_gas: float,
    x_start: float,
    x_end: float,
    gas_inf: tuple[float, list[float]] | None = None,
) -> Any:
    """A wall side's surface result ``r`` under ``provider``'s film (#386).

    ``side`` is the wall's gas-side object: an element with a flow ``area``
    (a duct) or a MomentumChamberNode. ``m_gas`` is the flow its htc reads
    (the duct's own, the chamber's throughput). ``coolant`` is, per source of
    ``film_sources(provider)``, its supply node's state and its signed flow.
    A source without forward flow means no film: r comes back unchanged.

    ``gas_inf`` is (Tt, Y) of the gas the holes discharge into before their
    coolant mixed in (#694), when the side is fed by that node. eta is
    defined against it (Baldauf 2002, Colban 2011: eta = (T_inf - T_aw) /
    (T_inf - T_c)), so T_aw shifts by T_inf - T_mix and DR takes the
    undiluted gas. Without it the side's own, mixed state stands in; that
    reads the wall cooler by (1 - eta) times the mixed-in coolant's share of
    T_inf - T_c. Mixing-out far downstream in a confined duct is outside
    both correlations either way.
    """
    area = float(getattr(side, "area", 0.0) or 0.0)
    mg = abs(float(m_gas))
    m_c = [float(m) for _, m in coolant]
    if r is None or mg <= 0.0 or area <= 0.0 or not m_c or min(m_c) <= 0.0:
        return r
    sources = film_sources(provider)
    G_gas = mg / area
    M = [m / src.hole_area / G_gas for m, src in zip(m_c, sources, strict=True)]
    T_c = [float(st.Tt) for st, _ in coolant]
    if gas_inf is None:
        T_g = float(gas_state.Tt)
        X_g = list(gas_state.X)
        T_shift = None
    else:
        T_g = float(gas_inf[0])
        X_g = list(cb.mass_to_mole(list(gas_inf[1])))
        T_shift = T_g - float(gas_state.Tt)
    W_g = float(cb.mwmix(X_g))
    # DR is rho_c / rho_g where the jets meet the gas: the gas moves at G_gas,
    # so rho_g is at its STATIC T, the duct closure's at (T_g, Pt, G_gas)
    # (#782 audit; the effusion path takes the chamber's static density). The
    # coolant comes from a plenum at rest: its Tt.
    s_g = cb.duct_htc_static_state(T_g, float(gas_state.Pt), G_gas, X_g)
    T_gs = float(s_g.T)
    dlnTs = (
        float(s_g.dT_dTt) / T_gs,
        float(s_g.dT_dPt) / T_gs,
        float(s_g.dT_dG) / (T_gs * area),
    )
    DR = [
        float(cb.mwmix(st.X)) / W_g * (T_gs / Tc) for (st, _), Tc in zip(coolant, T_c, strict=True)
    ]
    f = film_groups(provider, x_start, x_end, M, DR)
    eta = float(f.eta)
    m_tot = sum(m_c)
    film = {
        "eta_film": eta,
        # Flow-weighted over the sources (one source: its own).
        "film_blowing_ratio": sum(m * v for m, v in zip(m_c, M, strict=True)) / m_tot,
        "film_density_ratio": sum(m * v for m, v in zip(m_c, DR, strict=True)) / m_tot,
        "film_extrapolated": float(bool(f.extrapolated)),
        # The coolant temperature the film carries: sum w_g T_c,g / eta.
        "T_coolant_film": (
            sum(float(wg) * Tc for wg, Tc in zip(f.weight, T_c, strict=True)) / eta
            if eta > 0.0
            else sum(m * Tc for m, Tc in zip(m_c, T_c, strict=True)) / m_tot
        ),
        "T_aw_unfilmed": float(r.T_aw) + (T_shift or 0.0),
        "film_gas_undiluted": float(T_shift is not None),
    }
    filmed = _FilmedHtc(r, f, M, DR, T_c, m_c, mg, dlnTs, film, T_shift)
    side_htc = getattr(provider, "film_side_htc", None)
    if side_htc is None or len(sources) != 1:
        return filmed
    st_c = coolant[0][0]
    G_s = m_c[0] / sources[0].hole_area
    res = side_htc(
        x_start, x_end, G_s, M[0], float(r.h), float(st_c.Tt), float(st_c.P), list(st_c.X)
    )
    if res is None:
        return filmed
    # Which branch and box state applied (#962): the side's h with the film,
    # the duct's own, the film's span mean over its in-box part, its share of
    # the span, the hand-back weight in m, Re_s and the in-box flags.
    film.update(
        {
            "film_h": float(res.h),
            "film_h_duct": float(r.h),
            "film_h_film": float(res.h_film),
            "film_h_fraction_in_box": float(res.fraction_in_box),
            "film_h_weight_M": float(res.weight_M),
            "film_h_Re_s": float(res.Re_s),
            "film_h_branch": float(res.branch),
            "film_h_in_box": float(bool(res.x_in_box and res.Re_in_box and res.M_in_box)),
            "film_h_x_in_box": float(bool(res.x_in_box)),
            "film_h_Re_in_box": float(bool(res.Re_in_box)),
            "film_h_M_in_box": float(bool(res.M_in_box)),
        }
    )
    return _FilmSideHtc(filmed, res, M[0], m_c[0], mg, sources[0].hole_area)


@dataclass
class _ImpingementPlateResult:
    """What ``ImpingementPlateElement.htc_and_T`` returns to the wall code.

    ``dh_dsources`` is the part no other element has: the derivative of ``h``
    with respect to each crossflow source's mass flow at this plate's
    ``to_node``, keyed by element id. The solver relays it through the
    element's ``network_flow_inputs()`` (#465).
    """

    h: float
    Nu: float
    Re: float
    Pr: float
    T_aw: float
    Gc_Gj: float
    extrapolated: bool = False
    dh_dmdot: float = 0.0
    dh_dT: float = 0.0
    dT_aw_dmdot: float = 0.0
    dT_aw_dT: float = 1.0
    dh_dsources: dict[str, float] = field(default_factory=dict)
    # T_aw is the CR-3936 reference temperature (#467): it moves with the
    # neighbours' flows (the bypass share, Gc/Gj) and the bypass temperature.
    dT_aw_dsources: dict[str, float] = field(default_factory=dict)
    dT_aw_dT_nodes: dict[str, float] = field(default_factory=dict)
    eta_r: float = 0.0
    bypass_fraction: float = 0.0
    eta_r_extrapolated: bool = False


class ImpingementPlateElement(OrificeElement):
    """One spanwise row of an impingement jet plate, with its target wall (#465).

    The element's mass flow IS the row's jet flow: it leaves the supply plenum
    (``from_node``) through ``n_holes`` holes and joins the crossflow channel
    between jet plate and target (``to_node``). So the same element carries
    the orifice flow AND the Florschuetz, Truman and Metzger (1981) heat
    transfer on the target, which is what ``ImpingementModel`` on a
    ``ChannelElement`` cannot do: there the element's flow is the channel's.

    A full array is a chain, one plate per row:

        supply plenum ---+---------------+---------------+
                         |               |               |
                    [Plate 1]       [Plate 2]       [Plate 3]
                         |               |               |
        crossflow:      c1 --[Channel]-- c2 --[Channel]-- c3 --[Channel]-- exit

    CROSSFLOW FROM THE NETWORK. Florschuetz's crossflow-to-jet mass velocity
    ratio is

        Gc/Gj = (m_c / m_j) (pi/4) / ((yn/d)(z/d))

    with ``m_j`` this row's jet flow and ``m_c`` the crossflow APPROACHING the
    row: every inflow to ``to_node`` except this element's own. In the chain
    above that is the channel from the previous crossflow node, which carries
    rows 1..i-1 -- Florschuetz's Eq. 8 evaluated at x - xn/2. The ratio is
    read from the solved flows, not from the uniform-supply closed form, so a
    non-uniform supply, a row of different geometry, or a bleed shows up in
    the heat transfer. ``h`` therefore depends on a NEIGHBOUR's mass flow; the
    element exposes that as ``dh_dsources`` and the solver relays it into the
    Jacobian.

    An INITIAL crossflow (flow entering upstream of row 1, #467) enters
    through an ``ImpingementCrossflowElement(bypass_inlet=True)`` and counts
    in ``m_c`` like any other. The flow distribution is checked against
    Florschuetz, Metzger and Su (1982, NASA CR-3630) Table D.1. The heat
    transfer stays Florschuetz et al. (1981), who had none, at the local
    Gc/Gj -- CR-3630's own usage. A bypass at ANOTHER temperature moves the
    target's reference temperature (Florschuetz and Su 1985, NASA CR-3936,
    Eq. 2.11; see "THE TARGET" below). ``Gc/Gj`` beyond 0.8 is
    flagged by the set's own range check.

    THE TARGET. ``surface.area`` is the target footprint of the holes,
    ``n_holes * xn * yn``. ``T_aw`` is the reference temperature for h:

        T_aw = T_j + eta_r(Gc/Gj, z/d) phi_0 (T_0 - T_j)

    T_j the supply (Florschuetz's own reference, plenum-fed so static and
    total coincide), T_0 the bypass node's temperature, phi_0 the bypass
    share of the crossflow approaching this row (its lead-in flow over m_c),
    eta_r ``cb.impingement_eta_r``: a closure FITTED to CR-3936 Appendix A
    (the source publishes none). Only the bypass's share: the wall's own
    heating of the spent air is already inside Florschuetz's h, whose tests
    had no initial crossflow, so counting it again would double it. With no
    bypass, or one at T_j, T_aw is T_j exactly. ``eta_r_extrapolated`` in
    the diagnostics flags geometry or Gc/Gj outside the inline data; a
    staggered array is outside it too (the summary says so). The heat
    leaves the wall into ``to_node``: the spent air carries it downstream.

    DISCHARGE COEFFICIENT. The hole's own Cd: default ``'fixed'`` at
    Florschuetz's 0.79 (measured 0.73-0.85 across their plates), the value
    an unset ``Cd`` takes (``hole_cd=HoleCd(correlation=..., Cd=...)``). ``Cd``
    is the 'fixed' value only: passed with a
    correlation it raises (#872). Selectable alternatives, never chosen
    automatically: ``'IdelchikThick'``,
    ``'Lichtarowicz'`` (long hole, l/d 2-10) and ``'McGreehanSchotsch'``
    (whose ``U1/Vi`` is the SUPPLY-side approach velocity, zero for a
    plenum). ``Cd_in_range`` in the diagnostics says whether the plate is
    inside the chosen correlation's source range.

    CROSSFLOW ON THE HOLE (#800). The crossflow passing the jets' exits bends
    them over and lowers the Cd, steeply at large Gc/Gj: Rogers & Hersh
    (1976, NASA CR-2681) measured Cd ~ 0.08 at Vi/V_inf 0.008 (a hole
    discharging into a grazing stream 125 times faster), and Florschuetz,
    Metzger & Su (1982, NASA CR-3630 Table 4.1) on these very plates.
    ``hole_cd=HoleCd(exit_crossflow=...)`` selects the law ("auto" is the
    default below; ``cd_crossflow=`` is deprecated since 0.12.0, removed in
    0.13.0, and the ``cd_crossflow`` attribute reads the law in use), in q =
    G_ideal / G_c (C++
    ``rogers_hersh_1976_cd`` / ``florschuetz_metzger_su_1982_cd``):

    * ``'rogers_hersh_1976'`` (DEFAULT): CR-2681's single-hole fits,
      0.97 (Gj/Gc)^0.534 injecting, 0.80 (Gj/Gc)^0.47 for a reversed jet,
      capped at the hole's own Cd. Generic (any hole); below Gc/Gj ~1.5 at
      Cd 0.79 the cap holds, so Florschuetz's arrays (Gc/Gj <= 0.77) are
      unchanged. Cross-source for any impingement plate: ACCURACY.
    * ``'florschuetz_metzger_su_1982'``: CR-3630 Table 4.1, the authors'
      own law for the plates they tested (B(5,4,z), B(5,8,1); B(10,x,1)
      inferred by them). Its own Cd replaces the hole's. Refused for any
      other plate. The law the authors reduced Table D.1 with, so the one
      to score D.1 like for like.
    * ``'none'``: the hole's own Cd at any crossflow (the constant-Cd
      model of Florschuetz's Eq. 8 and CR-3630's Eq. 4.10).

    Resolved explicitly in the ideal jet velocity, the law makes the row's
    flow go as dP^1.07 under a strong crossflow (a linear resistance) instead
    of sqrt(dP): the stagnant rows of CR-3630 B(5,4,1) no longer alternate
    in sign (#799). G_ideal = sqrt(2 rho dP) from the supply's (jet's) static
    density and the supply-to-channel drop; G_c the approaching crossflow
    over the channel section ``n_holes yn z d^2``.

    NOT MODELLED. The temperature sensitivity of the correlation's
    properties (``dh_dT`` is 0, the same documented gap the channel
    impingement path has), the self-heated crossflow's temperature beyond
    what Florschuetz's plenum-referenced h holds (only a bypass's share
    moves T_aw), and the jet plate's own heat pick-up.

    Parameters
    ----------
    d_jet : float
        Hole diameter [m].
    xn_d, yn_d, z_d : float
        Streamwise pitch, spanwise pitch and plate-to-target gap over d_jet.
    span : float
        Plate span [m]; ``n_holes = round(span / (yn_d d_jet))``.
    plate_thickness : float
        Jet plate thickness [m]: the hole length the Cd correlations read.
    pattern : {'inline', 'staggered'}
        Selects Florschuetz's inline or staggered set when
        ``correlation_set`` is None.
    row : int or None
        Only for diagnostics: when given, ``Gc_Gj_closed_form`` reports Eq. 8
        at this row beside the network value.
    Nu_multiplier : float
        The user's rig-matching knob. 1.0 unless they set it.
    """

    @property
    def residual_scale_kind(self) -> str:
        """RESIDUAL-ROW-SCALING: its residual is the orifice's, m - m_calc,
        kg/s; Pt_in - Pt_out, Pa, where the holes' total area equals an
        adjacent momentum chamber's and there is no plate (#1016)."""
        return "p" if self._pass_through else "mdot"

    # htc_and_T's dh_dmdot comes from C++ in the signed jet flow (#481).
    _htc_dmdot_is_signed = True

    _selects_compressible_model = False  # keeps the nozzle (#861)
    # The deprecated hole-Cd keywords (HoleCd replaces them): not serialized.
    _deprecated_init_kwargs = ("correlation", "Cd", "cd_crossflow")
    # Florschuetz's measured 0.79: the Cd of an unset Cd with correlation
    # 'fixed' (#872), and the start guess's Cd under a correlation, as it
    # was before 0.8.0.
    _fixed_Cd_default = cb.FLORSCHUETZ_1981_DEFAULT_CD
    _seed_Cd_default = cb.FLORSCHUETZ_1981_DEFAULT_CD

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        d_jet: float,
        xn_d: float,
        yn_d: float,
        z_d: float,
        span: float,
        plate_thickness: float,
        pattern: Literal["inline", "staggered"] = "inline",
        correlation: str | None = None,
        Cd: float | None = None,
        correlation_set: object | None = None,
        row: int | None = None,
        Nu_multiplier: float = 1.0,
        edge_radius: float = 0.0,
        cd_crossflow: Literal["rogers_hersh_1976", "florschuetz_metzger_su_1982", "none"]
        | None = None,
        hole_cd: HoleCd | None = None,
    ) -> None:
        # The hole Cd's model choices (hole_cd.py): a HoleCd, or the
        # deprecated per-element keywords; 'fixed' (0.79) by default, the
        # jet crossflow law ("auto") Rogers & Hersh.
        spec = _with_correlation(
            _hole_cd_spec(
                f"ImpingementPlateElement {id!r}",
                hole_cd,
                {"correlation": correlation, "Cd": Cd, "cd_crossflow": cd_crossflow},
                stacklevel=3,
            ),
            "fixed",
            f"ImpingementPlateElement {id!r}",
        )
        # Refused first, as before: the law and the plate it needs.
        jet_law = _resolve_jet_row(f"ImpingementPlateElement {id!r}", spec, xn_d, yn_d, z_d)
        for name, val in (
            ("d_jet", d_jet),
            ("xn_d", xn_d),
            ("yn_d", yn_d),
            ("z_d", z_d),
            ("span", span),
            ("plate_thickness", plate_thickness),
            ("Nu_multiplier", Nu_multiplier),
        ):
            if val <= 0.0:
                raise ValueError(f"ImpingementPlateElement: {name} must be positive")
        if pattern not in ("inline", "staggered"):
            raise ValueError("ImpingementPlateElement: pattern is 'inline' or 'staggered'")

        n_exact = span / (yn_d * d_jet)
        n_holes = int(round(n_exact))
        if n_holes < 1:
            raise ValueError(
                f"ImpingementPlateElement: the span holds {n_exact:.3g} holes at this "
                "spanwise pitch, which rounds to none."
            )

        self.d_jet = d_jet
        self.xn_d = xn_d
        self.yn_d = yn_d
        self.z_d = z_d
        self.span = span
        self.pattern = pattern
        self.row = row
        self.Nu_multiplier = float(Nu_multiplier)
        self.n_holes = n_holes
        # The crossflow channel's section, as Florschuetz's Gc/Gj reads it:
        # (pi/4) / ((yn/d)(z/d)) is the jet row's area over it.
        self.channel_area = n_holes * yn_d * z_d * d_jet * d_jet
        self.hole_count_exact = n_exact
        if correlation_set is None:
            correlation_set = (
                cb.florschuetz_1981_inline()
                if pattern == "inline"
                else cb.florschuetz_1981_staggered()
            )
        self.correlation_set = correlation_set

        super().__init__(
            id,
            from_node,
            to_node,
            Cd=spec.Cd,
            area=n_holes * math.pi * d_jet * d_jet / 4.0,
            correlation=spec.correlation,
            plate_thickness=plate_thickness,
            edge_radius=edge_radius,
        )
        # The crossflow law on the hole's Cd (#800), the chain's last term.
        self._hole_cd.jet_crossflow = jet_law
        self.hole_cd = dataclasses.replace(spec, Cd=self.Cd)
        # The target footprint of the holes this element flows.
        self.surface = ConvectiveSurface(area=n_holes * xn_d * yn_d * d_jet * d_jet)
        self._row_geometry = cb.JetRowGeometry(
            d=d_jet, xn_d=xn_d, yn_d=yn_d, z_d=z_d, n_holes=float(n_holes)
        )
        self._crossflow_sources: list[str] = []
        # The bypass lead-in upstream of this row, and the node it starts at.
        self._bypass_segment: str | None = None
        self._bypass_node: str | None = None
        self._bypass_outlet: str | None = None

    @property
    def cd_crossflow(self) -> str:
        """The crossflow law on the hole's Cd (_CD_CROSSFLOW_LAWS), held by
        the hole-Cd chain."""
        return self._hole_cd.jet_crossflow

    @property
    def _cd_crossflow_suspended(self) -> bool:
        """Set by the solver's constant-Cd stage only (NetworkSolver
        ._solve_constant_cd_staged): the law off, the hole's own Cd."""
        return self._hole_cd.jet_crossflow_suspended

    @_cd_crossflow_suspended.setter
    def _cd_crossflow_suspended(self, value: bool) -> None:
        self._hole_cd.jet_crossflow_suspended = bool(value)

    @property
    def has_convective_surface(self) -> bool:
        return True

    def _hole_count(self) -> float:
        return float(self.n_holes)

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """A plate has no pipe, so no beta; the correlations see one hole.
        The crossflow sources are every other inflow to ``to_node``."""
        self._orifice_geom = cb.OrificeGeometry()
        self._orifice_geom.d = self.d_jet
        self._orifice_geom.D = 0.0
        self._orifice_geom.t = self.plate_thickness
        self._orifice_geom.r = self.edge_radius
        self.beta = 0.0
        self._crossflow_sources = [
            e.id for e in graph.get_upstream_elements(self.to_node) if e.id != self.id
        ]
        self._source_elements = {eid: graph.elements[eid] for eid in self._crossflow_sources}
        self._topology_error = self._configuration_error(graph)
        seg = self._find_bypass(graph)
        self._bypass_segment = seg.id if seg is not None else None
        self._bypass_node = seg.from_node if seg is not None else None
        self._bypass_outlet = seg.to_node if seg is not None else None
        self._resolve_plate_against_chambers(graph)

    def _find_bypass(self, graph: "FlowNetwork") -> "ImpingementCrossflowElement | None":
        """The bypass_inlet segment up the crossflow chain, if any."""
        node, seen = self.to_node, set()
        while node not in seen:
            seen.add(node)
            segments = [
                e
                for e in graph.get_upstream_elements(node)
                if isinstance(e, ImpingementCrossflowElement)
            ]
            if not segments:
                return None
            if segments[0].bypass_inlet:
                return segments[0]
            node = segments[0].from_node
        return None

    def _configuration_error(self, graph: "FlowNetwork") -> str | None:
        """Why this wiring is not the configuration the model covers, or None.

        Florschuetz et al. (1981) is ONE configuration: spent air leaving down
        the channel, the crossflow generated by the rows themselves. A
        hand-wired network can describe others that look alike but are not:

        * spent air leaving through the target or between the jets (#468):
          the row's node must drain through one ImpingementCrossflowElement;
        * an external, bypass crossflow (#467): every stream feeding the
          crossflow chain upstream must be another row's jets, unless it
          enters through a segment declared ``bypass_inlet`` -- the initial
          crossflow of Florschuetz, Metzger and Su (1982), stated, not
          inferred from what happens to be wired;
        * two plates merging at one node would read each other's jets as
          crossflow.
        """
        out = graph.get_downstream_elements(self.to_node)
        if len(out) != 1 or not isinstance(out[0], ImpingementCrossflowElement):
            return (
                f"its node {self.to_node!r} must drain through exactly one "
                "ImpingementCrossflowElement. Spent air leaving through the "
                "target (impingement-effusion/film) or between the jets is a "
                "different configuration with no model yet (#468)."
            )
        node, seen, via = self.to_node, set(), None
        while node not in seen:
            seen.add(node)
            inflow = [e for e in graph.get_upstream_elements(node) if e.id != self.id]
            plates = [e for e in inflow if isinstance(e, ImpingementPlateElement)]
            segments = [e for e in inflow if isinstance(e, ImpingementCrossflowElement)]
            if len(plates) > (0 if node == self.to_node else 1):
                return (
                    f"node {node!r} has more than one plate merging; each row "
                    "needs its own crossflow node."
                )
            if via is not None and via.bypass_inlet and plates:
                return (
                    f"the bypass inlet {via.id!r} starts at {node!r}, where a "
                    "jet row also discharges; the bypass enters the gap's "
                    "upstream end, ahead of row 1."
                )
            other = [e for e in inflow if e not in plates and e not in segments]
            if other and not (via is not None and via.bypass_inlet):
                return (
                    f"crossflow node {node!r} is fed by {other[0].id!r}, which is "
                    "not a jet row. An initial (bypass) crossflow enters through "
                    "a segment declared bypass_inlet=True (#467)."
                )
            if not segments:
                break
            via = segments[0]
            node = via.from_node
        return None

    def validate(self) -> None:
        super().validate()
        if self.correlation in _METERING_SELECTORS:
            raise ValueError(
                f"ImpingementPlateElement {self.id!r}: {self.correlation!r} is a "
                "metering correlation for a plate in a pipe; a jet plate has no "
                "pipe and no beta. Use 'fixed', 'IdelchikThick', 'Lichtarowicz' "
                "or 'McGreehanSchotsch'."
            )
        if getattr(self, "_topology_error", None):
            raise ValueError(
                f"ImpingementPlateElement {self.id!r}: covers jet rows whose spent "
                "air leaves down the crossflow channel (Florschuetz et al. "
                f"1981) only; {self._topology_error}"
            )
        if error := self._geometry_error():
            raise ValueError(error)

    def _geometry_error(self) -> str:
        """Why the jet holes cannot be drilled at these pitches (#1024), or "".

        The holes' total area is below the plate's, pi/4 / (xn/d yn/d) < 1,
        and holes leave material between them: yn/d > 1 along a row; between
        rows the nearest neighbour, xn/d > 1 inline, and staggered (the
        offset row half a spanwise pitch over, Florschuetz's pattern) both
        sqrt((xn/d)^2 + (yn/(2d))^2) > 1 and 2 xn/d > 1 (the same column two
        rows on). Equal is refused: touching holes leave no plate."""
        name = f"ImpingementPlateElement {self.id!r}"
        xn, yn, d = self.xn_d, self.yn_d, self.d_jet
        open_area = math.pi / 4.0 / (xn * yn)
        if not open_area < 1.0:
            return (
                f"{name}: the jet holes' area is {open_area:.3g} times the plate's, "
                f"pi/4 / (xn_d yn_d) at xn_d = {xn:g}, yn_d = {yn:g}: the holes do not "
                "fit in the plate. It must be below 1."
            )
        if self.pattern == "inline":
            return _hole_spacing_error(
                name, d, 90.0, ("yn_d * d_jet", yn * d), ("xn_d * d_jet", xn * d)
            )
        if error := _hole_spacing_error(name, d, 90.0, ("yn_d * d_jet", yn * d), None):
            return error
        nearest = min(math.hypot(xn, 0.5 * yn), 2.0 * xn)
        if not nearest > 1.0:
            return (
                f"{name}: staggered at xn_d = {xn:g}, yn_d = {yn:g}, the nearest holes "
                f"of neighbouring rows are {nearest * d:.6g} m apart (min(sqrt(xn^2 + "
                f"(yn/2)^2), 2 xn) = {nearest:.4g} d), not more than the hole diameter "
                f"d = {d:.6g} m: they touch or overlap and leave no material between them."
            )
        return ""

    def network_flow_inputs(self) -> list[tuple[str, str]]:
        """(element id, node id) of every flow ``htc_and_T`` reads besides its
        own: the crossflow sources, each read at ``to_node``, and the bypass
        lead-in at its outlet (the bypass share of the crossflow, #467)."""
        out = [(eid, self.to_node) for eid in self._crossflow_sources]
        if self._bypass_segment is not None and self._bypass_segment not in self._crossflow_sources:
            out.append((self._bypass_segment, self._bypass_outlet))
        return out

    def network_temperature_inputs(self) -> list[str]:
        """Nodes whose total temperature ``htc_and_T`` reads: the bypass's."""
        return [self._bypass_node] if self._bypass_node is not None else []

    def residuals(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        """The orifice flow; ``flows`` carries the approaching crossflow the
        hole's Cd reads (the crossflow law, see the class docstring)."""
        return super().residuals(state_in, state_out, flows)

    def htc_and_T(
        self,
        state: NetworkMixtureState,
        flows: dict[str, float] | None = None,
        temperatures: dict[str, float] | None = None,
    ) -> _ImpingementPlateResult:
        """Target-side h from this row's jets and the network's crossflow,
        and the reference temperature with any bypass.

        ``state`` is the supply (``from_node``) state with ``m_dot`` set to
        this element's flow; ``flows`` maps each crossflow source to its flow
        into ``to_node`` (and the bypass lead-in to its outlet flow);
        ``temperatures`` the bypass node to its total temperature. Without
        ``flows`` the row sees no crossflow. The physics and its flow
        derivatives are C++'s ``jet_row_heat_transfer``.
        """
        flows = flows or {}
        cs = cb.complete_state(state.Tt, state.P, state.X)
        tr = cs.transport
        m_c = float(sum(flows.get(eid, 0.0) for eid in self._crossflow_sources))
        row = cb.jet_row_heat_transfer(
            self.correlation_set, self._row_geometry, float(state.m_dot), m_c, tr.mu, tr.k, tr.Pr
        )
        mult = self.Nu_multiplier
        T_j = float(state.Tt)
        # CR-3936's reference temperature, on the bypass's share only (see
        # the class docstring): T_aw = T_j + eta_r phi_0 (T_0 - T_j).
        T_0 = T_j
        if self._bypass_node is not None and temperatures:
            T_0 = float(temperatures.get(self._bypass_node, T_j))
        m_b = float(flows.get(self._bypass_segment, 0.0)) if self._bypass_segment else 0.0
        dT0 = T_0 - T_j
        has_share = m_b > 0.0 and m_c > 0.0
        phi = m_b / m_c if has_share else 0.0
        eta = float(row.eta_r)
        dT_aw_dsources: dict[str, float] = {}
        if has_share and dT0 != 0.0:
            # d/dm_s of eta(G) m_b / m_c over every crossflow source s, and
            # d/dm_b on the bypass lead-in (one and the same at row 1).
            d_src = (row.deta_r_dm_crossflow * phi - eta * m_b / (m_c * m_c)) * dT0
            dT_aw_dsources = dict.fromkeys(self._crossflow_sources, d_src)
            dT_aw_dsources[self._bypass_segment] = (
                dT_aw_dsources.get(self._bypass_segment, 0.0) + eta / m_c * dT0
            )
        return _ImpingementPlateResult(
            h=row.h * mult,
            Nu=row.Nu,
            Re=row.Re_j,
            Pr=tr.Pr,
            T_aw=T_j + eta * phi * dT0,
            Gc_Gj=row.Gc_Gj,
            extrapolated=bool(row.extrapolated),
            dh_dmdot=row.dh_dm_jet * mult,
            dh_dsources=dict.fromkeys(self._crossflow_sources, row.dh_dm_crossflow * mult),
            dT_aw_dmdot=row.deta_r_dm_jet * phi * dT0,
            dT_aw_dT=1.0 - eta * phi,
            dT_aw_dsources=dT_aw_dsources,
            dT_aw_dT_nodes=({self._bypass_node: eta * phi} if has_share else {}),
            eta_r=eta,
            bypass_fraction=phi,
            eta_r_extrapolated=bool(row.eta_r_extrapolated),
        )

    def diagnostics(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
        temperatures: dict[str, float] | None = None,
    ) -> dict[str, float]:
        base = super().diagnostics(state_in, state_out, flows)
        if not base:
            return base
        r = self.htc_and_T(state_in, flows=flows, temperatures=temperatures)
        own_cd = self._hole_cd.hole_cd(state_in, state_out, flows)
        out = {
            "n_holes": float(self.n_holes),
            "Cd_own": float(own_cd),
            **self._hole_cd.jet_crossflow_envelope(state_in, state_out, flows, own_cd),
            "Re_j": float(r.Re),
            "Nu": float(r.Nu),
            "St": _stanton(r),
            "htc": float(r.h),
            "T_aw": float(r.T_aw),
            "Gc_Gj": float(r.Gc_Gj),
            "surface_extrapolated": float(r.extrapolated),
            "eta_r": float(r.eta_r),
            "bypass_fraction": float(r.bypass_fraction),
            "eta_r_extrapolated": float(r.eta_r_extrapolated),
        }
        if self.row is not None:
            # Florschuetz's uniform-supply Eq. 8 at this row, for comparison:
            # their constant-Cd model, at the hole's own Cd.
            out["Gc_Gj_closed_form"] = float(
                cb.crossflow_to_jet_ratio_at_row(self.yn_d, self.z_d, own_cd, self.row)
            )
        return {**base, **out}


class EffectiveAreaConnectionElement(OrificeElement):
    """
    An orifice element with user-specified effective area (Cd * A product).

    This element calculates orifice flow using the incompressible flow equation:
        m_dot = A_eff * sqrt(2 * rho * dP)

    where A_eff is the effective area (the product Cd * A).

    **Incompressible Formulation**: The solver uses `cb.orifice_mdot_and_jacobian()`,
    which applies the incompressible Bernoulli equation with regularization for numerical
    stability. Density is evaluated at upstream conditions. For compressible flows with
    significant Mach number effects, use a nozzle element instead (future feature).

    **Effective Area Interpretation**: The effective area represents the combined effect
    of discharge coefficient and geometric area. For example:
    - A_eff = 0.01 m^2 could represent Cd=1.0 * A=0.01 m^2
    - Or equivalently: Cd=0.8 * A=0.0125 m^2
    - The product (Cd * A) is what matters for flow calculation

    Unlike LosslessConnectionElement, this element produces pressure drop proportional to flow rate.
    Use this when you need a simple area-based flow restriction without detailed geometry modeling.

    The effective area is user-specified and already accounts for all geometric effects
    (entrance/exit losses, contraction, vena contracta, etc.), so upstream/downstream
    geometry discovery is not performed.
    """

    _selects_compressible_model = False  # keeps the nozzle (#861)
    _init_kwarg_attrs = {"effective_area": "area"}  # stored as self.area (#1099)

    def __init__(self, id: str, from_node: str, to_node: str, effective_area: float) -> None:
        """
        Initialize with effective area (Cd * A product).

        Args:
            id: Element identifier
            from_node: Upstream node ID
            to_node: Downstream node ID
            effective_area: Effective area for flow calculation (m^2).
                           This is the product Cd * A, where Cd is the discharge coefficient
                           and A is the geometric area. Stored internally as self.area.

        Example:
            >>> conn = EffectiveAreaConnectionElement("conn1", "inlet", "outlet", 0.01)
            >>> conn.area  # 0.01 m^2 (effective area = Cd * A)
            >>> conn.Cd    # 1.0 (normalized discharge coefficient)

        Note:
            The effective area already includes all loss coefficients. For example, if you
            have a sharp-edged orifice with geometric area 0.0125 m^2 and Cd=0.8, you would
            specify effective_area=0.01 m^2 (0.8 * 0.0125).
        """

        if not (math.isfinite(effective_area) and effective_area >= 0.0):
            raise ValueError(
                f"EffectiveAreaConnectionElement {id!r}: effective_area = {effective_area!r} "
                "must be finite and non-negative (0 closes the restriction: no flow)."
            )
        diameter = math.sqrt(4.0 * effective_area / math.pi)
        super().__init__(id, from_node, to_node, Cd=1.0, diameter=diameter, correlation="fixed")

        # Force self.area to be precisely the effective area
        # (to remove any floating point math.sqrt / pi precision issues)
        self.area = effective_area

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """
        Skip upstream/downstream geometry discovery.

        Unlike OrificeElement, we do not discover upstream/downstream diameters because
        the effective area is user-specified and already accounts for all geometric effects.
        No Beta ratio correction is needed or applied.
        """
        # Intentionally empty - effective area is pre-computed by user
        pass

    def n_equations(self) -> int:
        """Return the number of equations (1: mass flow balance)."""
        return 1


class DiameterDischargeCoefficientConnectionElement(OrificeElement):
    """
    An orifice element with user-specified physical diameter and discharge coefficient or loss coefficient.

    This element calculates orifice flow using the incompressible flow equation:
        m_dot = A * Cd * sqrt(2 * rho * dP)

    where A is the computed physical area and Cd is the discharge coefficient.

    **Incompressible Formulation**: Uses the Bernoulli equation with density evaluated
    at upstream conditions. Valid for low Mach number flows (M < 0.3 typically).

    **Parameters**: User provides either Cd (discharge coefficient) or zeta (loss coefficient).
    The other parameter is calculated automatically using the relationship:
        zeta = 1/Cd^2 - 1  or  Cd = 1/sqrt(zeta + 1)

    **Effective Area**: The effective area is calculated as A_eff = A * Cd.

    Unlike LosslessConnectionElement, this element produces pressure drop proportional to flow rate.
    Use this when you have a physical diameter measurement and want to specify loss characteristics
    via discharge coefficient or loss coefficient.
    """

    _selects_compressible_model = False  # keeps the nozzle (#861)

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        diameter: float,
        Cd: float | None = None,
        zeta: float | None = None,
    ) -> None:
        """
        Initialize with physical diameter and either Cd or zeta.

        Args:
            id: Element identifier
            from_node: Upstream node ID
            to_node: Downstream node ID
            diameter: Physical geometric diameter (m)
            Cd: Discharge coefficient (optional, 0 < Cd <= 1)
            zeta: Loss coefficient (optional, zeta >= 0)

        Note:
            Exactly one of Cd or zeta must be provided. If both are provided, zeta is ignored.
            If Cd is provided, zeta is calculated as zeta = 1/Cd^2 - 1.
            If zeta is provided, Cd is calculated as Cd = 1/sqrt(zeta + 1).

        Example:
            >>> # Using Cd
            >>> conn = DiameterDischargeCoefficientConnectionElement("conn1", "inlet", "outlet", 0.1, Cd=0.8)
            >>> conn.diameter  # 0.1 m (physical diameter)
            >>> conn.Cd    # 0.8 (discharge coefficient)

            >>> # Using zeta
            >>> conn = DiameterDischargeCoefficientConnectionElement("conn1", "inlet", "outlet", 0.1, zeta=0.5625)
            >>> conn.diameter  # 0.1 m (physical diameter)
            >>> conn.Cd    # 0.8 (calculated from zeta)

        Raises:
            ValueError: If neither Cd nor zeta is provided, or if invalid values are given.
        """
        if Cd is not None and zeta is not None:
            # If both are provided, use Cd and ignore zeta
            zeta = None
        elif Cd is None and zeta is None:
            raise ValueError("Either Cd or zeta must be provided")

        if Cd is not None:
            if not (0 < Cd <= 1):
                raise ValueError("Cd must be in range (0, 1]")
        else:  # zeta is provided
            if zeta < 0:
                raise ValueError("zeta must be >= 0")
            # Calculate Cd from zeta: Cd = 1/sqrt(zeta + 1)
            Cd = 1.0 / (zeta + 1.0) ** 0.5

        # Pass physical diameter and Cd directly to parent
        super().__init__(id, from_node, to_node, Cd=Cd, diameter=diameter, correlation="fixed")

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """
        Skip upstream/downstream geometry discovery.

        The physical area and discharge coefficient are user-specified and already account
        for all geometric effects. No Beta ratio correction is needed or applied.
        """
        # Intentionally empty - area and Cd are pre-computed by user
        pass

    def n_equations(self) -> int:
        """Return the number of equations (1: mass flow balance)."""
        return 1


class LosslessConnectionElement(NetworkElement):
    """
    An ideal connection with no friction or momentum loss.
    Inherently preserves total pressure: Pt_in = Pt_out.
    Does not require geometric parameters.
    """

    #: RESIDUAL-ROW-SCALING: a total-pressure equality, Pa.
    residual_scale_kind = "p"

    def __init__(self, id: str, from_node: str, to_node: str):
        super().__init__(id, from_node, to_node)
        self.area: float | None = None
        self.diameter: float | None = None

    def _seed_resistance(self, rho: float, P_ref: float, m_ref: float) -> float | None:
        """Start guess: no loss; zero, which the seed raises to a small
        fraction of the network's smallest resistance."""
        return 0.0

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass  # Zero-loss elements do not need upstream geometry

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def residuals(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Residual: Pt_in - Pt_out = 0
        res = [state_in.Pt - state_out.Pt]
        jac = {0: {f"{self.from_node}.Pt": 1.0, f"{self.to_node}.Pt": -1.0}}
        return res, jac

    def get_spatial_profile(self, n_steps: int = 100):
        """
        Returns a mock spatial profile array representing the zero-loss state.
        This provides structural symmetry with realistic pipes (e.g. Fanno/Rough)
        for GUI plotters.
        """

        # Simplified placeholder for test completion validation
        # Normally would utilize the true state nodes
        profile = []
        for _ in range(n_steps):
            st = cb.IncompressibleStation()
            st.x = 0.0
            st.P = 100000.0
            st.T = 300.0
            st.rho = 1.2
            st.v = 1.0
            st.M = 0.0
            st.h = 300000.0
            profile.append(st)
        return profile

    def diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float | str]:
        if state_in.P <= 0 or state_in.Tt <= 0:
            return {}

        cs_in = cb.complete_state(state_in.T, state_in.P, state_in.X)
        ref = _element_reference_block(
            cs_in, m_dot=state_in.m_dot, area=self.area or 0.0, Dh=self.diameter, location="inlet"
        )
        v_in = ref["velocity"]
        mach_in = v_in / cs_in.thermo.a if cs_in.thermo.a > 0 and v_in > 0 else 0.0

        return {
            "m_dot": float(state_in.m_dot),
            **_element_pressure_block(state_in, state_out, mach_in=mach_in, mach_out=mach_in),
            **ref,
        }

    def n_equations(self) -> int:
        return 1


#: |sum of a flow splitter's fractions - 1| above which the set is refused
#: (#903), never renormalised. A tolerance for the round-off of typed decimal
#: fractions (0.1 + 0.2 + 0.7 is 1 + 2.2e-16), not for a mis-typed set.
SPLITTER_FRACTION_SUM_TOL = 1e-9


class SplitterOutletElement(NetworkElement):
    """One outlet of a flow splitter (#903): a mass-flow controller.

    A splitter takes one inlet node's flow to N outlet nodes in given
    fractions f_i (summing to 1). ``FlowSplitter`` builds it as N of these
    two-port elements, all from the inlet node, so the solver, the state
    propagation, the second-law check and the ledger see plain two-ports.
    Each outlet carries one unknown, its flow ``m_dot``, and one row:

    * a **metered** outlet: ``m_i - f_i sum_j m_j = 0`` over the splitter's
      own outlets, linear, with the exact Jacobian ``1 - f_i`` on its own flow
      and ``-f_i`` on each sibling's. The row holds no pressure, so the
      outlet's total-pressure drop ``Pt_in - Pt_out`` is free: it is whatever
      the path downstream leaves to the valve, and it is reported as the
      metering the valve would need (``valve_dPt``, ``valve_CdA`` and, with
      an ``area``, ``valve_Cd`` and ``valve_zeta``).
    * the **reference** outlet: ``Pt_in - Pt_out = 0``, its valve fully open.
      Of N fractions summing to 1 only N - 1 are independent; with a valve on
      every outlet the total flow (or, fed at a fixed flow, the inlet
      pressure) would have no equation. The reference's own fraction follows
      from the others. It must be the most resistive path: a metered outlet
      that needs MORE total pressure than the reference delivers needs a pump
      (``valve_dPt < 0``), and the second-law check refuses that root. With
      ``auto_reference`` the solve then moves the reference
      (``flow_splitter.solve_with_reference_swaps``); without, or when no
      reference works, the refusal names the outlet to make the reference.

    Each outlet carries the inlet's Tt and composition (a two-port moves its
    upstream state); Pt follows each outlet's own path.

    Args:
        id: Element id (``FlowSplitter`` uses ``{splitter}__{outlet node}``).
        from_node: The splitter's inlet node, shared by all its outlets.
        to_node: This outlet's node.
        fraction: f_i in (0, 1), of the splitter's throughput.
        outlets: Ids of all the splitter's outlet elements, this one included.
        reference: True for the one unmetered outlet (``Pt_out = Pt_in``).
        area: Valve reference area [m^2] for ``valve_Cd`` and ``valve_zeta``;
            None reports the effective area ``valve_CdA`` only.
        splitter: The splitter's id, for messages and reports.
        auto_reference: Let ``NetworkSolver.solve`` move the reference to the
            outlet that needs the most pressure when a solve is refused for a
            pump at a metered outlet (``flow_splitter.solve_with_reference_swaps``).
            The same for every outlet of one splitter.
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        fraction: float,
        outlets: list[str],
        reference: bool = False,
        area: float | None = None,
        splitter: str | None = None,
        auto_reference: bool = True,
    ):
        super().__init__(id, from_node, to_node)
        self.auto_reference = bool(auto_reference)
        self.fraction = float(fraction)
        self.outlets = list(outlets)
        self.reference = bool(reference)
        self.area = area
        self.splitter = splitter if splitter is not None else id
        #: {outlet element id: its to_node}, resolved in resolve_topology.
        self._outlet_nodes: dict[str, str] = {}
        self._reference_id: str = ""

    @property
    def residual_scale_kind(self) -> str:
        """RESIDUAL-ROW-SCALING: a metered outlet's row is a mass-flow
        balance (kg/s); the reference's a total-pressure equality (Pa)."""
        return "p" if self.reference else "mdot"

    def _seed_resistance(self, rho: float, P_ref: float, m_ref: float) -> float | None:
        """Start guess: no loss, as a lossless connection. The fractions are
        linear rows, which the first Newton step meets exactly."""
        return 0.0

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def n_equations(self) -> int:
        return 1

    def _siblings(self) -> list[str]:
        return [eid for eid in self.outlets if eid != self.id]

    def network_flow_inputs(self) -> list[tuple[str, str]]:
        """The sibling outlets' flows, at the shared inlet node."""
        return [(eid, self.from_node) for eid in self._siblings()]

    def validate(self) -> None:
        where = f"SplitterOutletElement '{self.id}' (splitter '{self.splitter}')"
        if not (math.isfinite(self.fraction) and 0.0 < self.fraction < 1.0):
            raise ValueError(
                f"{where}: fraction must be in (0, 1), got {self.fraction}. An outlet "
                "that carries nothing is no outlet: remove it."
            )
        if len(set(self.outlets)) != len(self.outlets) or len(self.outlets) < 2:
            raise ValueError(f"{where}: needs at least two distinct outlets, got {self.outlets}.")
        if self.id not in self.outlets:
            raise ValueError(f"{where}: its own id is not among the splitter's outlets.")
        if self.area is not None and not (math.isfinite(self.area) and self.area > 0.0):
            raise ValueError(f"{where}: area must be positive, got {self.area}.")

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """Check the splitter as a whole: its outlets exist, share this inlet
        and outlet list, hold exactly one reference, and their fractions sum
        to 1 (within SPLITTER_FRACTION_SUM_TOL, never renormalised)."""
        where = f"Flow splitter '{self.splitter}'"
        members: list[SplitterOutletElement] = []
        for eid in self.outlets:
            e = graph.elements.get(eid)
            if not isinstance(e, SplitterOutletElement):
                raise ValueError(
                    f"{where}: outlet element '{eid}' is not in the network as a "
                    "SplitterOutletElement."
                )
            if (
                e.from_node != self.from_node
                or e.outlets != self.outlets
                or e.auto_reference != self.auto_reference
            ):
                raise ValueError(
                    f"{where}: outlet '{eid}' does not share the inlet '{self.from_node}', "
                    "the outlet list and the auto_reference setting of the others."
                )
            members.append(e)
        refs = [e.id for e in members if e.reference]
        if len(refs) != 1:
            raise ValueError(
                f"{where}: needs exactly one reference outlet (its valve open, "
                f"Pt_out = Pt_in), found {len(refs)}: {refs}."
            )
        total = math.fsum(e.fraction for e in members)
        if not abs(total - 1.0) <= SPLITTER_FRACTION_SUM_TOL:
            raise ValueError(
                f"{where}: the fractions sum to {total!r}, not 1 (tolerance "
                f"{SPLITTER_FRACTION_SUM_TOL:g}); they are not renormalised. "
                "Declare one outlet the remainder to take what the others leave."
            )
        self._outlet_nodes = {e.id: e.to_node for e in members}
        self._reference_id = refs[0]

    def _throughput(self, m_self: float, flows: dict[str, float] | None) -> tuple[float, list[str]]:
        """(sum of the splitter's outlet flows, the sibling ids)."""
        sibs = self._siblings()
        flows = flows or {}
        missing = [eid for eid in sibs if eid not in flows]
        if missing:
            raise ValueError(
                f"SplitterOutletElement '{self.id}': the flows of its sibling outlets "
                f"{missing} were not passed (flows=); the solver passes them."
            )
        return m_self + math.fsum(float(flows[eid]) for eid in sibs), sibs

    def residuals(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        if self.reference:
            # The open valve: no loss, Pt_in - Pt_out = 0.
            res = [state_in.Pt - state_out.Pt]
            jac = {0: {f"{self.from_node}.Pt": 1.0, f"{self.to_node}.Pt": -1.0}}
            return res, jac
        m_self = float(state_in.m_dot)
        M, sibs = self._throughput(m_self, flows)
        f = self.fraction
        row: dict[str, float] = {f"{self.id}.m_dot": 1.0 - f}
        for eid in sibs:
            row[f"{eid}.m_dot"] = -f
        return [m_self - f * M], {0: row}

    def valve(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState, m_dot: float
    ) -> dict[str, float]:
        """The metering this outlet's valve would need at the solved state.

        ``valve_dPt`` = Pt_in - Pt_out [Pa]. With flow through a drop,
        ``valve_CdA`` [m^2] is the effective area of an isentropic jet from
        the inlet's stagnation state (Tt_in, Pt_in) to the outlet's static
        pressure, held at the sonic pressure below it (``valve_choked``), the
        flux from C++ ``mass_flux_isentropic``. With an ``area``, ``valve_Cd``
        = CdA / area and ``valve_zeta`` = 1 / Cd^2 (``cb.zeta_from_Cd``), which
        at low Mach is dPt over the dynamic head at that area.
        """
        from .second_law import SIGMA_TOL

        dPt = float(state_in.Pt - state_out.Pt)
        # A rise in Pt is a pump where the second-law check sees one: an
        # adiabatic two-port's sigma_norm is -ln(Pt_out / Pt_in) ~ dPt / Pt_in.
        pump = not self.reference and state_in.Pt > 0.0 and dPt < -SIGMA_TOL * state_in.Pt
        out: dict[str, float] = {"valve_dPt": dPt, "valve_needs_pump": 1.0 if pump else 0.0}
        if self.reference or not (dPt > 0.0 and m_dot > 0.0):
            return out
        Tt, Pt, X = float(state_in.Tt), float(state_in.Pt), state_in.X
        P_star = float(cb.critical_pressure_ratio(Tt, Pt, X)) * Pt
        P_exit = float(state_out.P)
        choked = P_exit <= P_star
        G = float(cb.mass_flux_isentropic(Tt, Pt, max(P_exit, P_star), X))
        if not G > 0.0:
            return out
        CdA = m_dot / G
        out["valve_CdA"] = CdA
        out["valve_choked"] = 1.0 if choked else 0.0
        if self.area is not None:
            Cd = CdA / float(self.area)
            out["valve_Cd"] = Cd
            out["valve_zeta"] = float(cb.zeta_from_Cd(Cd))
        return out

    def diagnostics(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
    ) -> dict[str, float | str]:
        if state_in.P <= 0 or state_in.Tt <= 0:
            return {}
        m = float(state_in.m_dot)
        cs_in = cb.complete_state(state_in.T, state_in.P, state_in.X)
        ref = _element_reference_block(
            cs_in, m_dot=m, area=self.area or 0.0, Dh=None, location="inlet"
        )
        out: dict[str, float | str] = {
            "m_dot": m,
            **_element_pressure_block(state_in, state_out, mach_in=0.0, mach_out=0.0),
            **ref,
            "fraction": self.fraction,
            "is_reference": 1.0 if self.reference else 0.0,
            **self.valve(state_in, state_out, m),
        }
        if flows is not None:
            M, _ = self._throughput(m, flows)
            out["m_dot_splitter"] = M
            if M != 0.0:
                out["fraction_achieved"] = m / M
        return out

    def second_law_hint(self, sol: dict[str, Any]) -> str:
        """Why the second-law check refused this outlet: it needs its total
        pressure raised to carry its fraction (a pump, not a valve). Names the
        outlet that should be the reference: the one needing the most."""
        if self.reference:
            return ""
        Pt_in = sol.get(f"{self.from_node}.Pt")
        Pt_out = sol.get(f"{self.to_node}.Pt")
        if Pt_in is None or Pt_out is None or not float(Pt_out) > float(Pt_in):
            return ""
        need = {
            eid: float(sol.get(f"{nid}.Pt", float("nan")))
            for eid, nid in self._outlet_nodes.items()
        }
        worst = max(need, key=lambda k: need[k] if math.isfinite(need[k]) else -math.inf)
        ref_node = self._outlet_nodes.get(self._reference_id, self._reference_id)
        return (
            f"Flow splitter '{self.splitter}': outlet '{self.to_node}' needs its total "
            f"pressure raised by {float(Pt_out) - float(Pt_in):.4g} Pa to carry fraction "
            f"{self.fraction:g}, which no valve does. The reference outlet "
            f"'{ref_node}' is not the most resistive path; make "
            f"'{self._outlet_nodes.get(worst, worst)}' the reference."
        )


class TurbomachineElement(NetworkElement):
    """A compressor or turbine at a fixed total-to-total pressure ratio (#893).

    The pressure ratio is a parameter and the network sets the flow: one
    unknown, ``m_dot``, and one row, ``Pt_out - f Pt_in = 0`` (f = PR for a
    compressor, 1/PR for a turbine). The stream leaves at the outlet total
    temperature of ``cb.compression`` / ``cb.expansion`` with exactly one of
    ``eta_poly`` and ``eta_is``: the solver delivers it to the downstream node
    at that Tt (``outlet_stagnation``), with the Jacobian of Tt_out in Tt_in
    and the composition. At a fixed PR the outlet Tt does not depend on the
    pressure level, so it has no pressure column. The composition is frozen.

    Reported: ``power`` [W] > 0 (absorbed by a compressor, delivered by a
    turbine), ``dh`` [J/kg] signed into the flow, ``Tt_out``, both efficiencies
    and ``pressure_ratio``. A root with the flow reversed (surge, a turbine
    driven backwards) is refused, never returned. Between two pressure
    boundaries nothing sets the flow: refused.
    """

    #: RESIDUAL-ROW-SCALING: a total-pressure relation, Pa.
    residual_scale_kind = "p"
    #: True for a turbine (expansion, PR = P_in / P_out).
    expansion = False

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        pressure_ratio: float,
        eta_poly: float | None = None,
        eta_is: float | None = None,
    ) -> None:
        super().__init__(id, from_node, to_node)
        kind = type(self).__name__
        if (eta_poly is None) == (eta_is is None):
            raise ValueError(f"{kind} '{id}': give exactly one of eta_poly and eta_is")
        eta = eta_poly if eta_poly is not None else eta_is
        if not (0.0 < float(eta) <= 1.0):
            raise ValueError(
                f"{kind} '{id}': efficiency {eta} is outside (0, 1]: above 1 the "
                "machine would destroy entropy"
            )
        if not math.isfinite(float(pressure_ratio)):
            # An infinite ratio passed the >= 1 test and fixed no pressure (#1024).
            raise ValueError(
                f"{kind} '{id}': pressure_ratio {pressure_ratio} must be a finite number "
                ">= 1 (it is Pt_out/Pt_in for a compressor, Pt_in/Pt_out for a turbine)"
            )
        if not (float(pressure_ratio) >= 1.0):
            raise ValueError(
                f"{kind} '{id}': pressure_ratio {pressure_ratio} < 1 (it is Pt_out/Pt_in "
                "for a compressor, Pt_in/Pt_out for a turbine)"
            )
        self.pressure_ratio = float(pressure_ratio)
        self.eta_poly = None if eta_poly is None else float(eta_poly)
        self.eta_is = None if eta_is is None else float(eta_is)
        self.area: float | None = None

    @property
    def pt_factor(self) -> float:
        """Pt_out / Pt_in."""
        return 1.0 / self.pressure_ratio if self.expansion else self.pressure_ratio

    def path(self, Tt_in: float, Y: list[float]) -> Any:
        """The C++ path (cb.compression / cb.expansion) from Tt_in, Y (mass)."""
        fn = cb.expansion if self.expansion else cb.compression
        return fn(
            float(Tt_in),
            list(cb.mass_to_mole(list(Y))),
            self.pressure_ratio,
            eta_poly=self.eta_poly,
            eta_is=self.eta_is,
        )

    def outlet_stagnation(
        self, state_in: NetworkMixtureState
    ) -> tuple[float, float, list[float], float]:
        """(Tt_out, dTt_out/dTt_in, dTt_out/dY, Pt_out/Pt_in) of the stream
        this element delivers: the solver mixes it into the outlet node."""
        p = self.path(state_in.Tt, state_in.Y)
        return float(p.T2), float(p.dT2_dT1), list(p.dT2_dY), self.pt_factor

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        if isinstance(graph.nodes.get(self.from_node), PressureBoundary) and isinstance(
            graph.nodes.get(self.to_node), PressureBoundary
        ):
            raise ValueError(
                f"{type(self).__name__} '{self.id}' runs between two pressure boundaries: "
                "it fixes their pressure ratio, so nothing sets its flow. Put a lossy "
                "element (orifice, channel) or a mass-flow boundary in its path."
            )

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def residuals(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        f = self.pt_factor
        res = [state_out.Pt - f * state_in.Pt]
        jac = {0: {f"{self.to_node}.Pt": 1.0, f"{self.from_node}.Pt": -f}}
        return res, jac

    def diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float | str]:
        if state_in.Tt <= 0.0:
            return {}
        p = self.path(state_in.Tt, state_in.Y)
        m = float(state_in.m_dot)
        return {
            "m_dot": m,
            "pressure_ratio": self.pressure_ratio,
            "Tt_in": float(state_in.Tt),
            "Tt_out": float(p.T2),
            "Tt_out_isentropic": float(p.T2s),
            "dh": float(p.dh),
            "power": abs(m * float(p.dh)),
            "eta_poly": float(p.eta_poly),
            "eta_is": float(p.eta_is),
            "Pt_in": float(state_in.Pt),
            "Pt_out": float(state_out.Pt),
        }

    def n_equations(self) -> int:
        return 1


class CompressorElement(TurbomachineElement):
    """A compressor (#893): Pt_out = PR Pt_in, Tt_out from cb.compression.
    See TurbomachineElement."""

    expansion = False


class TurbineElement(TurbomachineElement):
    """A turbine (#893): Pt_out = Pt_in / PR, Tt_out from cb.expansion.
    See TurbomachineElement."""

    expansion = True


class PortDischargeElement(NetworkElement):
    """A junction port's momentum chamber discharging into a PressureBoundary.

    The link the GUI inserts where a junction port (ejector, tee) is wired
    straight to a boundary (#637). A LosslessConnectionElement held the
    chamber's TOTAL pressure at the boundary's, i.e. made the reservoir
    recover the jet's whole head: no root once the chamber is sonic, which is
    exactly a critical-mode ejector's outlet. What the jet does instead (#360's
    exit-head rule for a channel):

    - **out of the chamber, subsonic:** its STATIC pressure is the boundary's
      and the dynamic head is lost;
    - **out of the chamber, choked:** it leaves under-expanded, at the sonic
      static pressure P* of its flux (C++ sonic_static_pressure), above the
      boundary's: P_chamber = max(P_boundary, P*), smoothed over
      PORT_DISCHARGE_CHOKE_BLEND;
    - **into the chamber:** the boundary is a stagnation source, Pt equal.

    A boundary with coupling="total" keeps the Pt equality both ways. The two
    directions blend in the link's flow over FLOW_REVERSAL_M_SCALE, where they
    differ only by the chamber's dynamic head, which vanishes with the flow.
    An incompressible chamber has no sonic limit, so no floor.
    """

    #: RESIDUAL-ROW-SCALING: a pressure equality, Pa.
    residual_scale_kind = "p"

    def __init__(self, id: str, from_node: str, to_node: str):
        super().__init__(id, from_node, to_node)
        self.area: float | None = None
        self.diameter: float | None = None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        self._chamber = graph.nodes.get(self.from_node)
        self._boundary = graph.nodes.get(self.to_node)

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def _outflow_row(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> tuple[float, dict[str, float], float]:
        """(residual, Jacobian, P*) of the discharge out of the chamber."""
        c, b = self.from_node, self.to_node
        boundary = getattr(self, "_boundary", None)
        static = boundary is not None and getattr(boundary, "exit_head_lost", None) is not None
        if not (static and boundary.exit_head_lost(True)):
            return float(state_in.Pt - state_out.Pt), {f"{c}.Pt": 1.0, f"{b}.Pt": -1.0}, 0.0
        P_b = float(state_out.Pt)
        chamber = getattr(self, "_chamber", None)
        A = float(getattr(chamber, "area", 0.0) or 0.0)
        compressible = getattr(chamber, "regime", "compressible") == "compressible"
        if not (compressible and A > 0.0):
            return float(state_in.P) - P_b, {f"{c}.P": 1.0, f"{b}.Pt": -1.0}, 0.0
        m = float(state_in.m_dot)
        s = _solver_tools.sonic_static_pressure(float(state_in.Tt), abs(m) / A, state_in.Y)
        P_star = float(s.P)
        eps = PORT_DISCHARGE_CHOKE_BLEND * P_b
        d = P_star - P_b
        root = math.sqrt(d * d + eps * eps)
        floor = 0.5 * (P_b + P_star + root)
        w_star = 0.5 * (1.0 + d / root)  # d floor / d P*
        jac = {
            f"{c}.P": 1.0,
            f"{b}.Pt": -(1.0 - w_star),
            f"{self.id}.m_dot": -w_star * float(s.dP_dG) * (1.0 if m >= 0.0 else -1.0) / A,
            # T and Y are relayed through the chamber's propagation chain,
            # as the chamber's own closure row does.
            f"{c}.Tt": -w_star * float(s.dP_dTt),
        }
        for k, dy in enumerate(s.dP_dY):
            if dy != 0.0:
                jac[f"{c}.Y[{k}]"] = -w_star * float(dy)
        return float(state_in.P) - floor, jac, P_star

    def residuals(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        c, b = self.from_node, self.to_node
        r_in, j_in = float(state_in.Pt - state_out.Pt), {f"{c}.Pt": 1.0, f"{b}.Pt": -1.0}
        # The branch follows the FLOW, this row's own unknown. A blend in Pt
        # (the orifice's, whose flow is computed from the pressures) admits a
        # spurious root here: inflow from the boundary with the chamber's
        # static pinned to it and its Pt above it.
        m_name = f"{self.id}.m_dot"
        w, dw = _flow_blend(float(state_in.m_dot))
        if w == 0.0:
            return [r_in], {0: j_in}
        r_out, j_out, _ = self._outflow_row(state_in, state_out)
        if w == 1.0:
            return [r_out], {0: j_out}
        jac: dict[str, float] = {}
        for k in dict.fromkeys([*j_in, *j_out]):  # deterministic order (see _flow_and_partials)
            jac[k] = w * j_out.get(k, 0.0) + (1.0 - w) * j_in.get(k, 0.0)
        jac[m_name] = jac.get(m_name, 0.0) + (r_out - r_in) * dw
        return [w * r_out + (1.0 - w) * r_in], {0: jac}

    def diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float | str]:
        if state_in.P <= 0 or state_in.Tt <= 0:
            return {}
        _, _, P_star = self._outflow_row(state_in, state_out)
        chamber = getattr(self, "_chamber", None)
        mach = 0.0
        if chamber is not None and float(getattr(chamber, "area", 0.0) or 0.0) > 0.0:
            # The chamber's own Mach: the face this link discharges.
            mach = float(chamber._closure(state_in, float(state_in.m_dot)).M)
        out: dict[str, float | str] = {
            "m_dot": float(state_in.m_dot),
            **_element_pressure_block(state_in, state_out, mach_in=mach, mach_out=mach),
        }
        if P_star > 0.0:
            # The sonic static pressure of the chamber's flux, and whether the
            # discharge is choked (under-expanded into the boundary).
            out["P_sonic"] = P_star
            out["exit_choked"] = float(P_star > float(state_out.Pt))
        return out

    def n_equations(self) -> int:
        return 1


# ---------------------------------------------------------------------------
# Discrete pressure loss element: combustor-theta aware or plain cold-flow.
# ---------------------------------------------------------------------------


def _correlation_is_linear_theta(correlation: Callable) -> bool:
    """True if the correlation's value depends on ``ctx.theta``.

    Detected by invoking the correlation twice with different theta values on
    a throwaway context and checking whether the returned xi differs, or by
    checking the ``dxi_dtheta`` return.  The callable is also considered
    theta-aware if its returned ``dxi_dtheta`` is nonzero.
    """
    try:
        ctx_zero = PressureLossContext(
            theta=0.0,
            Tt_in=300.0,
            Pt_in=101325.0,
            X_in=[1.0],
            T_ad=300.0,
            phi=0.0,
            mdot_total=1.0,
            mdot_fuel=0.0,
            mdot_air=0.0,
            X_products=[],
            Y_products=[],
        )
        xi0, dxi_dtheta = correlation(ctx_zero)
    except Exception:
        return False
    if dxi_dtheta != 0.0:
        return True
    ctx_one = SimpleNamespace(**vars(ctx_zero))
    ctx_one.theta = 1.0
    try:
        xi1, _ = correlation(ctx_one)
    except Exception:
        return False
    return abs(xi1 - xi0) > 1e-12


class PressureLossElement(NetworkElement):
    """
    Discrete pressure-loss edge.  Enforces ``Pt_out = Pt_in * (1 - xi)``
    where ``xi`` is supplied by a user correlation callable.

    The loss follows the FLOW: "in" is the upstream node of the flow as it
    runs, and the correlation sees ``|m_dot|`` and that node's state. A
    reversed element used to apply its loss in the declared direction, which
    for a constant fraction raised total pressure by ``1 / (1 - xi)`` (#550),
    and a head loss returned no loss at all. Through zero flow the forward
    and reverse forms are blended (``FLOW_REVERSAL_M_SCALE``), as the area
    change does.

    Four standard correlations are provided in :mod:`combaero.network.pressure_loss`:

    - :class:`ConstantFractionLoss` -- fixed fractional total-pressure drop.
    - :class:`ConstantHeadLoss` -- Euler loss coefficient (dynamic head).
    - :class:`LinearThetaFractionLoss` -- ``xi = k * theta + xi0``.
    - :class:`LinearThetaHeadLoss` -- ``xi = (k * theta + zeta0) * q / P_in``.

    Theta-aware correlations need a reference node with ``has_theta = True``
    (e.g. :class:`CombustorNode`).  Discovery precedence:

    1. ``theta_source=<node_id>`` constructor kwarg (explicit override).
    2. ``to_node`` endpoint if it has ``has_theta = True``.
    3. ``from_node`` endpoint if it has ``has_theta = True``.
    4. Cold-flow fallback: ``theta = 0``.  If the correlation is linear-theta
       a :class:`UserWarning` is emitted at resolution time.

    Placing the element between two ``has_theta = True`` nodes without an
    explicit ``theta_source`` raises ``ValueError`` during
    :meth:`FlowNetwork.validate`.
    """

    #: RESIDUAL-ROW-SCALING: the row is Pt_dn - Pt_up (1 - xi), in Pa. It was
    #: filed as a mass-flow row by default, weighting it ref_p/ref_mdot
    #: (1e5-1e9) too heavily; LM then traded every pressure for it (#614).
    residual_scale_kind = "p"

    # to_dict serialises the given area, not the inherited one (#1099).
    _init_kwarg_attrs = {"area": lambda e: e.area if e._area_source == "user" else None}

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        correlation: Callable[[Any], tuple[float, float]],
        theta_source: str | None = None,
        area: float | None = None,
        surface: ConvectiveSurface | None = None,
    ):
        super().__init__(id, from_node, to_node)
        self.correlation = correlation
        self.theta_source = theta_source
        self.area = area
        #: Where self.area came from once resolved (issue #262).
        self._area_source = "user" if area is not None else ""
        self.surface = surface or ConvectiveSurface()
        # Resolved during resolve_topology.  Empty string sentinel disallowed
        # so callers can cleanly test ``if self._theta_source_resolved is not None``.
        self._theta_source_resolved: str | None = None
        self._both_endpoints_theta: bool = False

    @property
    def has_convective_surface(self) -> bool:
        return self.surface.area > 0.0

    def htc_and_T(self, state: NetworkMixtureState):
        """Convective HTC on the discrete loss element.

        Treats the element as a short duct of length ``L = Dh`` (L/D = 1),
        where ``Dh = sqrt(4*A/pi)`` is derived from the flow area. Returns
        ``None`` if either the convective surface or the flow area is missing.
        """

        if self.surface.area == 0.0 or not self.area or self.area <= 0.0:
            return None

        diameter = math.sqrt(4.0 * self.area / math.pi)
        length = diameter  # Short duct approximation: L/D = 1

        # At the duct's own static state, not the feed node's (#636).
        return _duct_htc(self.surface, state, self.area, diameter, length)

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def _seed_resistance(self, rho: float, P_ref: float, m_ref: float) -> float | None:
        """Start guess from the shipped correlations: a fraction loss xi P at
        its design flow (or the reference flow), a head loss zeta q at its
        area. Any other callable: no estimate."""
        c = self.correlation
        xi = getattr(c, "xi", getattr(c, "xi0", None))
        if xi is not None:
            m = getattr(c, "m_ref", None) or m_ref
            return float(xi) * P_ref / (m * m)
        zeta = getattr(c, "zeta", getattr(c, "zeta0", None))
        area = getattr(c, "area", None)
        if zeta is not None and area:
            return float(zeta) / (2.0 * rho * area * area)
        return None

    def _seed_fraction(self) -> tuple[float, float | None] | None:
        """(xi, m_ref) of a shipped fraction correlation, for the cold start's
        own law dPt = xi Pt_up (``NetworkSolver._resistance_network_guess``,
        #989): a flow-independent ratio, faded as the dynamic head below a
        design flow m_ref. None for a head loss, any other callable, and a
        linear-theta loss that reads a theta source: its xi = k theta + xi0
        waits on the combustion the seed does not know (a burner's exit duct
        seeded at xi0 started 190 times its flow and the cold solve stalled), so
        it keeps the reference-flow K (``_seed_resistance``)."""
        c = self.correlation
        xi = getattr(c, "xi", getattr(c, "xi0", None))
        if xi is None:
            return None
        if getattr(c, "k", 0.0) and self._theta_source_resolved is not None:
            return None
        m_ref = getattr(c, "m_ref", None)
        return float(xi), (float(m_ref) if m_ref else None)

    def n_equations(self) -> int:
        return 1

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Stash a graph reference so residuals/diagnostics can read the theta
        # source node's burned/unburned temperatures without a re-pass.
        self._graph_ref = graph

        # 0. Area inheritance: when no user-specified area, inherit from the
        #    nearest upstream ChannelElement diameter.  This makes discrete loss
        #    work after any element, not just combustor/plenum nodes.
        if self.area is None:
            for elem in graph.get_upstream_elements(self.from_node):
                if isinstance(elem, ChannelElement) and elem.diameter is not None:
                    self.area = math.pi / 4.0 * elem.diameter**2
                    self._area_source = f"channel '{elem.id}'"
                    break
            if self.area is None:
                # Widen to any area-bearing neighbour (a combustor or momentum
                # chamber upstream is the common case) before defaulting.
                area, source = _infer_flow_area(graph, [self.from_node, self.to_node], exclude=self)
                if area is not None:
                    self.area = area
                    self._area_source = source
            if self.area is None:
                # A head-loss correlation reads this area as its velocity
                # reference, so a wrong value rescales the loss rather than
                # failing. Warn only when a correlation will actually consume
                # it -- for a fraction-based correlation nothing reads it and
                # a nominal value is harmless (issue #262).
                self.area = DEFAULT_CHAMBER_AREA
                self._area_source = "default"
                if hasattr(self.correlation, "area"):
                    warnings.warn(
                        f"PressureLossElement '{self.id}': no flow area could be "
                        f"inferred from the topology, so the nominal "
                        f"{DEFAULT_CHAMBER_AREA} m^2 default is used as the "
                        f"velocity reference for "
                        f"{type(self.correlation).__name__}. The head loss scales "
                        f"as 1/area^2 -- set 'area' explicitly to make it "
                        f"meaningful.",
                        UserWarning,
                        stacklevel=2,
                    )
            # Propagate updated area into head-loss correlations. The convective
            # surface is left as given: a wetted area is not a flow area, and
            # an unset surface (0) means no convective surface (#1022).
            if hasattr(self.correlation, "area"):
                self.correlation.area = self.area

        # 1. Explicit override wins.
        if self.theta_source is not None:
            src = graph.nodes.get(self.theta_source)
            if src is None:
                raise ValueError(
                    f"PressureLossElement '{self.id}': theta_source "
                    f"'{self.theta_source}' is not in the network."
                )
            if not getattr(src, "has_theta", False):
                raise ValueError(
                    f"PressureLossElement '{self.id}': theta_source "
                    f"'{self.theta_source}' is a '{type(src).__name__}' which "
                    f"does not provide theta (has_theta=False)."
                )
            self._theta_source_resolved = self.theta_source
            return

        to_has_theta = getattr(graph.nodes[self.to_node], "has_theta", False)
        from_has_theta = getattr(graph.nodes[self.from_node], "has_theta", False)

        if to_has_theta and from_has_theta:
            # Ambiguous: raised at validate() for a cleaner user-facing error.
            self._both_endpoints_theta = True
            return
        if to_has_theta:
            self._theta_source_resolved = self.to_node
            return
        if from_has_theta:
            self._theta_source_resolved = self.from_node
            return

        # 4. Cold-flow fallback.  Warn if the correlation depends on theta.
        if _correlation_is_linear_theta(self.correlation):
            warnings.warn(
                f"PressureLossElement '{self.id}' uses a linear-theta "
                f"correlation ({type(self.correlation).__name__}) but neither "
                f"endpoint ('{self.from_node}' upstream, '{self.to_node}' "
                f"downstream) has has_theta=True. The k-term will have no "
                f"effect; xi falls back to cold-flow (theta=0).",
                UserWarning,
                stacklevel=2,
            )

    def _theta_temperatures(
        self, state_in: NetworkMixtureState, graph: "FlowNetwork | None"
    ) -> tuple[float, float, float]:
        """(theta, T_in, T_ad) the loss reads: the theta source's burned and
        unburned temperatures, T_in the unburned one; without a theta source
        theta = 0 and T_in = T_ad = the inlet's Tt."""
        if self._theta_source_resolved is not None and graph is not None:
            src = graph.nodes[self._theta_source_resolved]
            T_unburned = float(getattr(src, "_T_unburned", state_in.Tt))
            T_burned = float(getattr(src, "_T_burned", state_in.Tt))
            theta = (T_burned / T_unburned - 1.0) if T_unburned > 0 else 0.0
            return theta, T_unburned, T_burned
        T = float(state_in.Tt)
        return 0.0, T, T

    def _build_ctx(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        graph: "FlowNetwork | None",
    ) -> SimpleNamespace:
        """Build a duck-typed PressureLossContext from the current states (for
        a user callable; a built-in model reads ``_composed``)."""
        theta, T_ctx, T_ad = self._theta_temperatures(state_in, graph)
        # The inlet STAGNATION state: Lefebvre's loss factor is total-total
        # (#782 audit); until 0.8 P_in was the static P.
        return PressureLossContext(
            theta=theta,
            Tt_in=T_ctx,
            Pt_in=float(state_in.Pt),
            X_in=list(state_in.X),
            T_ad=T_ad,
            phi=0.0,
            mdot_total=float(state_in.m_dot),
            mdot_fuel=0.0,
            mdot_air=0.0,
            X_products=[],
            Y_products=[],
        )

    def _composed(
        self,
        model: "cb.PressureLossModel",
        state_in: NetworkMixtureState,
        theta: float,
        T_in: float,
        m_abs: float,
    ) -> "cb.PressureLossResult":
        """A built-in model in C++ (``cb.pressure_loss_compose``, #725) at the
        flow-upstream node's state: the value's density reads the state's
        composition X, the partials are by its mass-fraction unknowns Y (a
        fraction model reads neither). Positional: (mdot, Tt_in, Pt_in, Y_in,
        X_in, x), theta the first model input."""
        inp = cb.PressureLossInputs(
            m_abs, T_in, float(state_in.Pt), list(state_in.Y), list(state_in.X), [theta]
        )
        return cb.pressure_loss_compose(model, inp)

    def _directed(
        self,
        up: NetworkMixtureState,
        dn: NetworkMixtureState,
        up_id: str,
        dn_id: str,
        m_dot: float,
    ) -> tuple[float, Callable[[], dict[str, float]]]:
        """Residual and lazy Jacobian row of the loss with the flow up -> dn."""
        graph = getattr(self, "_graph_ref", None)
        model_of = getattr(self.correlation, "pressure_loss_model", None)
        theta, T_in, T_ad = self._theta_temperatures(up, graph)
        ctx: SimpleNamespace | None = None
        partials: Any = None
        try:
            if model_of is not None:
                partials = self._composed(model_of(), up, theta, T_in, abs(float(m_dot)))
                xi, dxi_dtheta = partials.xi, partials.dxi_dtheta
            else:
                ctx = self._build_ctx(up, dn, graph)
                ctx.mdot_total = abs(float(m_dot))
                xi, dxi_dtheta = self.correlation(ctx)
        except Exception as exc:
            raise RuntimeError(
                f"PressureLossElement '{self.id}' correlation raised: {exc}"
            ) from exc
        # Pt_dn - Pt_up * (1 - xi) = 0, in the direction the flow runs.
        res = dn.Pt - up.Pt * (1.0 - xi)
        m_sign = 1.0 if m_dot >= 0.0 else -1.0

        def row() -> dict[str, float]:
            return self._jacobian(
                up, graph, T_in, T_ad, xi, dxi_dtheta, up_id, dn_id, m_sign, partials, ctx
            )[0]

        return res, row

    def residuals(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        m = float(state_in.m_dot)
        w, dw_dm = _flow_blend(m)
        parts: list[tuple[float, float, Callable[[], dict[str, float]]]] = []
        if w > 0.0:
            r_f, j_f = self._directed(state_in, state_out, self.from_node, self.to_node, m)
            parts.append((w, r_f, j_f))
        if w < 1.0:
            # Written in the forward orientation (Pt_to - Pt_from - xi Pt_to):
            # blended as Pt_from - Pt_to (1 - xi), the Pt difference cancels at
            # w = 1/2 and the row loses its pressures.
            r_r, j_r = self._directed(state_out, state_in, self.to_node, self.from_node, m)

            def j_rev(j_r: Callable[[], dict[str, float]] = j_r) -> dict[str, float]:
                return {k: -v for k, v in j_r().items()}

            parts.append((1.0 - w, -r_r, j_rev))
        res = [sum(wt * r for wt, r, _ in parts)]

        # The Jacobian is analytic for a built-in model and differenced for a
        # user callable (one or two calls per variable and species), so it is
        # computed only when the solver reads it (LazyJacobian, #489).
        def jacobian() -> dict[int, dict[str, float]]:
            out: dict[str, float] = {}
            for wt, _, row in parts:
                for k, v in row().items():
                    out[k] = out.get(k, 0.0) + wt * v
            if len(parts) == 2 and dw_dm:
                key = f"{self.id}.m_dot"
                out[key] = out.get(key, 0.0) + dw_dm * (parts[0][1] - parts[1][1])
            return {0: out}

        return res, LazyJacobian(jacobian)

    def _jacobian(
        self,
        state_in: NetworkMixtureState,
        graph: "FlowNetwork | None",
        T_in: float,
        T_ad: float,
        xi: float,
        dxi_dtheta: float,
        up_id: str,
        dn_id: str,
        m_sign: float,
        partials: Any,
        ctx: SimpleNamespace | None,
    ) -> dict[int, dict[str, float]]:
        """d(Pt_dn - Pt_up (1 - xi))/dx with ``state_in`` the flow-upstream
        node's state; the correlation saw ``|m_dot|``, so its m column carries
        the flow's sign. ``partials`` is a built-in model's composed result
        (analytic); without it ``ctx`` is the user callable's context, and its
        other dependences are differenced."""
        jac: dict[int, dict[str, float]] = {
            0: {
                f"{up_id}.Pt": -(1.0 - xi),
                f"{dn_id}.Pt": 1.0,
            }
        }

        # ----- Analytical theta Jacobian (linear-theta correlations) -----
        # d_res/d_T_burned = -Pt_in * dxi_dtheta / Tt_in.
        if self._theta_source_resolved is not None and dxi_dtheta != 0.0 and T_in > 0:
            jac[0][f"{self._theta_source_resolved}.Tt"] = state_in.Pt * dxi_dtheta / T_in

        # ----- The other dependences: m_dot, Pt_in, Tt_in, composition. A
        # built-in model returns them analytically from C++
        # (pressure_loss_compose, #725); a user callable (the loss-map
        # prototype) is differenced.
        def _dxi_dvar(attr: str, val: float, step_scale: float) -> float:
            if partials is not None:
                return {
                    "mdot_total": partials.dxi_dmdot,
                    "Pt_in": partials.dxi_dPt_in,
                    "Tt_in": partials.dxi_dTt_in,
                }[attr]
            assert ctx is not None
            eps = step_scale if val == 0.0 else max(abs(val) * 1e-6, step_scale)
            saved = getattr(ctx, attr)
            try:
                setattr(ctx, attr, saved + eps)
                xi_p, _ = self.correlation(ctx)
                setattr(ctx, attr, saved - eps)
                xi_m, _ = self.correlation(ctx)
            finally:
                setattr(ctx, attr, saved)
            return (xi_p - xi_m) / (2.0 * eps)

        # Only bother if xi actually varies with the velocity/density variables
        # (skip for constant fraction to avoid unneeded correlation calls).
        m_abs = ctx.mdot_total if ctx is not None else 0.0
        dxi_dmdot = _dxi_dvar("mdot_total", m_abs, 1e-6)
        if dxi_dmdot != 0.0:
            jac[0][f"{self.id}.m_dot"] = state_in.Pt * dxi_dmdot * m_sign

        dxi_dPt_in = _dxi_dvar("Pt_in", float(state_in.Pt), 1.0)
        if dxi_dPt_in != 0.0:
            # Pt_in = state_in.Pt: into the same column as the -(1 - xi) term.
            key = f"{up_id}.Pt"
            jac[0][key] = jac[0].get(key, 0.0) + state_in.Pt * dxi_dPt_in

        if self._theta_source_resolved is None:
            # When no theta source, Tt_in = state_in.Tt (a derived state).
            # Emit a relay-routable Tt Jacobian on the upstream node.
            dxi_dTt_in = _dxi_dvar("Tt_in", T_in, 0.01)
            if dxi_dTt_in != 0.0:
                jac[0][f"{up_id}.Tt"] = state_in.Pt * dxi_dTt_in
        elif graph is not None and T_in > 0:
            self._add_unburned_jacobian(
                jac[0], graph, T_in, T_ad, state_in.Pt, dxi_dtheta, _dxi_dvar
            )

        # Composition of the inlet state (a head loss's density): per species,
        # other mass fractions fixed, relayed as "{node}.Y[k]". A combustor's
        # burned composition moves with its fuel flow; this column was
        # missing (#481 C1).
        Y_in = list(state_in.Y)
        if partials is not None:
            for k, d in enumerate(partials.dxi_dY):
                if d != 0.0 and Y_in[k] > 0.0:
                    jac[0][f"{up_id}.Y[{k}]"] = state_in.Pt * d
            return jac
        assert ctx is not None
        X_saved = ctx.X_in
        try:
            for k, y in enumerate(Y_in):
                if not y > 0.0:
                    continue
                # Central where y can take the step down, forward above a
                # trace species: a central step shrunk to 0.5 y underflowed
                # to zero at a denormal y (5e-324) and divided by it.
                h = 1e-6
                signs = (1.0, -1.0) if y >= 2.0 * h else (1.0, 0.0)
                xi_pm = []
                for sgn in signs:
                    Yk = list(Y_in)
                    Yk[k] += sgn * h
                    ctx.X_in = list(cb.mass_to_mole(Yk))
                    xi_pm.append(self.correlation(ctx)[0])
                d = (xi_pm[0] - xi_pm[1]) / ((signs[0] - signs[1]) * h)
                if d != 0.0:
                    jac[0][f"{up_id}.Y[{k}]"] = state_in.Pt * d
        finally:
            ctx.X_in = X_saved

        return jac

    def _add_unburned_jacobian(
        self,
        row: dict[str, float],
        graph: "FlowNetwork",
        T_u: float,
        T_b: float,
        Pt_in: float,
        dxi_dtheta: float,
        dxi_dvar: Callable[[str, float, float], float],
    ) -> None:
        """d(res)/d(unknowns) through the source's UNBURNED temperature (#481 C1).

        T_u enters twice: theta = T_b/T_u - 1 (d theta/dT_u = -T_b/T_u^2) and
        as the correlation's reference temperature Tt_in (a head loss's
        density). T_u is the enthalpy mix of the source's inflows (#718), so
        dT_u = sum(dT_u/dm_i dm_i + dT_u/dTt_i dTt_i + dT_u/dY_i . dY_i), the
        partials from the C++ mixer: the flows' unknowns directly, each Tt_i
        and Y_i through its node's relayed "{src}.Tt" and "{src}.Y[k]". It was
        missing: d/d(fuel flow) read 0 against 95 by differences.
        """
        src = graph.nodes[self._theta_source_resolved]
        streams = getattr(src, "_unburned_streams", None)
        if not streams:
            return
        dres_dTu = Pt_in * (dxi_dtheta * (-T_b / (T_u * T_u)) + dxi_dvar("Tt_in", T_u, 0.01))
        if dres_dTu == 0.0:
            return
        for names, src_node, dTu_dm, dTu_dTt, dTu_dY in streams:
            for var, coeff in names.items():
                row[var] = row.get(var, 0.0) + dres_dTu * dTu_dm * coeff
            if src_node is not None:
                key = f"{src_node}.Tt"
                row[key] = row.get(key, 0.0) + dres_dTu * dTu_dTt
                for k, d in enumerate(dTu_dY):
                    if d != 0.0:
                        key = f"{src_node}.Y[{k}]"
                        row[key] = row.get(key, 0.0) + dres_dTu * d

    @property
    def unanchored_fraction(self) -> bool:
        """A fraction loss with no design flow (#1004): ``ConstantFractionLoss``
        or ``LinearThetaFractionLoss``, ``m_ref`` None and a nonzero loss. It
        holds dPt = xi Pt at every flow but zero, so where the network gives
        it a smaller drop its only root is the forward/reverse blend at
        |m_dot| < FLOW_REVERSAL_M_SCALE, where the loss law is not met. The
        solver's ``fraction_unanchored`` key flags it there."""
        c = self.correlation
        if isinstance(c, ConstantFractionLoss):
            return c.m_ref is None and c.xi != 0.0
        if isinstance(c, LinearThetaFractionLoss):
            return c.m_ref is None and (c.k != 0.0 or c.xi0 != 0.0)
        return False

    def diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float | str]:
        if state_in.P <= 0 or state_in.Tt <= 0 or state_out.P <= 0 or state_out.Tt <= 0:
            return {}

        graph = getattr(self, "_graph_ref", None)
        # xi of the flow as it runs: from its upstream node, at |m_dot|.
        forward = state_in.m_dot >= 0.0
        up, dn = (state_in, state_out) if forward else (state_out, state_in)
        theta, T_in, _ = self._theta_temperatures(up, graph)
        model_of = getattr(self.correlation, "pressure_loss_model", None)
        try:
            if model_of is not None:
                xi = self._composed(model_of(), up, theta, T_in, abs(float(state_in.m_dot))).xi
            else:
                ctx = self._build_ctx(up, dn, graph)
                ctx.mdot_total = abs(float(state_in.m_dot))
                xi, _ = self.correlation(ctx)
        except Exception:
            xi = 0.0

        area = getattr(self, "area", None) or getattr(self.correlation, "area", None)
        Dh = math.sqrt(4.0 * area / math.pi) if area and area > 0 else None

        cs_in = cb.complete_state(state_in.T, state_in.P, state_in.X)
        ref = _element_reference_block(
            cs_in, m_dot=state_in.m_dot, area=area or 0.0, Dh=Dh, location="inlet"
        )
        v_in = ref["velocity"]
        mach_in = v_in / cs_in.thermo.a if cs_in.thermo.a > 0 and v_in > 0 else 0.0

        cs_out = cb.complete_state(state_out.T, state_out.P, state_out.X)
        rho_out = cs_out.thermo.rho
        v_out = (
            (abs(state_out.m_dot)) / ((rho_out) * (area))
            if area and area > 0 and rho_out > 0
            else 0.0
        )
        mach_out = v_out / cs_out.thermo.a if cs_out.thermo.a > 0 and v_out > 0 else 0.0

        diag: dict[str, float | str] = {
            "m_dot": float(state_in.m_dot),
            **_element_pressure_block(state_in, state_out, mach_in=mach_in, mach_out=mach_out),
            **ref,
            "xi": float(xi),
            "theta": float(theta),
            "area": float(area or 0.0),
            "area_source": self._area_source,
            # 1 where an unanchored fraction loss sits on its zero-flow blend
            # in a driven component, not a dead end (#1004). The element alone
            # cannot tell (it needs the network's inputs and topology): NaN until the
            # solver sets it; 0 for any other loss, which is anchored.
            "fraction_unanchored": math.nan if self.unanchored_fraction else 0.0,
        }

        if self.has_convective_surface:
            h_res = self.htc_and_T(state_in)
            if h_res is not None:
                diag["Nu"] = float(h_res.Nu)
                diag["St"] = _stanton(h_res)
                diag["htc"] = float(h_res.h)
                diag["T_aw"] = float(h_res.T_aw)
                diag["f"] = float(h_res.f)
            else:
                # No correlation evaluated (no flow): the keys stay, undefined.
                for key in ("Nu", "St", "htc", "T_aw", "f"):
                    diag[key] = math.nan

        return diag


# START GUESS ONLY: the turbulent friction factor of a duct's seed
# resistance (ChannelElement._seed_resistance); the residual computes its own.
SEED_FRICTION_FACTOR = 0.02


class _SurfaceFeed(NamedTuple):
    """The duct's static (T, P, rho) a ribbed or pin-fin drop reads at a
    reservoir end (#926), each with its partials keyed by unknown name."""

    T: float
    P: float
    rho: float
    drho: dict[str, float]
    dT: dict[str, float]
    dP: dict[str, float]


class ChannelElement(NetworkElement):
    """Duct with frictional pressure drop and optional convective surface.

    Couples STAGNATION pressures: Pt_up - Pt_down = dP_friction. The
    static-versus-stagnation distinction lives in the node constitutive
    relation, not here (see the note in ``residuals``).

    Pressure drop
        ``regime="incompressible"`` uses Darcy-Weisbach with a friction factor
        selected by ``friction_model`` from Haaland, Colebrook, Serghides or
        Petukhov; ``regime="compressible"`` uses Fanno flow. Both are computed
        in C++ -- ``friction.h`` owns the correlations and
        ``solver_interface.h`` the (f, J) entry points. Nothing about the
        friction factor is restated in this file.

    Convective surfaces
        A ``ConvectiveSurface`` adds heat transfer. Ribbed and pin-fin surfaces
        also OWN the drop: their correlations return the channel's friction
        factor (ribs: the four-sided Darcy value; pins: Metzger's per-row
        ``dP / (2 rho Vmax^2 N)``), and the element uses it directly rather
        than as a multiplier on pipe friction -- see ``_ribbed_residuals`` and
        ``_pin_fin_residuals``. Impingement surfaces leave the drop to the
        smooth path (model the jet plate with an ``OrificeElement``).
        Provenance lives with each correlation set: ``rib_correlation.h``,
        ``pin_fin_correlation.h``, ``impingement_correlation.h``.

    Known gaps
        The array correlations expose no dP sensitivity to static pressure or
        composition, so those Jacobian columns are absent rather than
        approximated; ``test_channel_element_jacobian_fd.py`` pins this as a
        strict xfail carrying the magnitude.
    """

    # The duct at both ends; the higher one decides.
    validity_mach_keys = ("mach_in", "mach_out")
    #: Whether diagnostics read the duct's own state at an end on a node at
    #: rest (#888, ``_rest_end_states``), from the node-to-node Fanno solution.
    #: A crossflow segment solves its duct between its own station states
    #: (#802), so it keeps the node states there.
    _rest_end_state_from_duct: bool = True

    # to_dict serialises the given diameter and Dh, not the inherited ones
    # (#1099).
    _init_kwarg_attrs = {
        "diameter": lambda e: e.diameter if e._area_source == "user" else None,
        "Dh": lambda e: e.Dh if e._Dh_source == "user" else None,
    }

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        length: float,
        diameter: float | None = None,
        Dh: float | None = None,
        roughness: float = 1e-5,
        regime: CompressibilityLiteral = DEFAULT_REGIME,
        friction_model: FrictionModelLiteral = "haaland",
        htc_model: HeatTransferModelLiteral = "none",
        t_hot: float | None = None,
        surface: ConvectiveSurface | None = None,
        default_diameter: float | None = None,
    ):
        super().__init__(id, from_node, to_node)
        self.length = length
        self.diameter: float | None = diameter
        #: Opt-in fallback [m] for a diameter that is neither set nor
        #: inferable (#570). None (the default) keeps the refusal of #262: the
        #: diameter sets friction and Mach, so a guessed value changes the
        #: flow. The GUI opts in with DEFAULT_CHANNEL_DIAMETER so that a
        #: freshly dropped chain of inheriting channels solves; the result is
        #: reported as ``area_source = "default"``, never as a user value.
        self.default_diameter = default_diameter
        # Dh: user-specified hydraulic diameter (for non-circular ducts).
        # None until resolved: defaults to diameter (circular assumption).
        self.Dh: float | None = Dh if Dh is not None else diameter
        self.roughness = roughness
        self.area: float | None = math.pi * (diameter / 2) ** 2 if diameter is not None else None
        # Provenance of the diameter/area: "user" when given, else the
        # _infer_flow_area source string once resolved (issue #559).
        self._area_source = "user" if diameter is not None else ""
        #: Dh provenance (#559): "user" for a non-circular override, else it is
        #: the circular Dh = D (whatever D's own source).
        self._Dh_source = "user" if Dh is not None else "derived from D (circular)"

        self.regime = regime
        self.friction_model = friction_model
        self.htc_model = htc_model
        self.t_hot = t_hot
        self.surface = surface or ConvectiveSurface()

    @property
    def has_convective_surface(self) -> bool:
        return True

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def _surface_rest_feed(
        self, state_in: NetworkMixtureState, m_dot: float, area: float
    ) -> "_SurfaceFeed | None":
        """The duct's static state a ribbed or pin-fin drop reads where the
        channel leaves a node at rest (#926), else None (the node state).

        Those residuals read the from node's (Tt, P). On a node whose gas is
        at rest (``state_at_rest``: a plenum, a combustor, a PressureBoundary
        or MassFlowBoundary) that is the reservoir, P = Pt: its density is
        rho(Tt, Pt), about M^2 / 2 above the duct's static density (13% at
        Mach 0.5, 46% at 0.9), and the drop came out low by as much. Here the
        surface reads the duct's own inlet static state instead: the
        isentropic acceleration from the node's (Tt, Pt) to the duct's mass
        flux |m| / A (C++ ``duct_htc_static_state``, the state its heat
        transfer already reads, #636). Partials in m_dot, the node's Tt and
        Pt, and its mass fractions are the state's own analytic ones. An end
        where a boundary takes the flow at its static pressure (an outflow
        losing its head) keeps the node state, as ``_rest_end_states`` does;
        a moving node (a momentum chamber) carries its own static state and
        keeps it too. A crossflow segment solves between its own stations.
        """
        node = getattr(self, "_upstream_node", None)
        if not self._rest_end_state_from_duct or area <= 0.0:
            return None
        if not bool(getattr(node, "state_at_rest", False)):
            return None
        if m_dot < 0.0 and hasattr(node, "exit_head_lost") and bool(node.exit_head_lost(True)):
            return None
        G = abs(float(m_dot)) / area
        X = list(state_in.X)
        s = cb.duct_htc_static_state(float(state_in.Tt), float(state_in.Pt), G, X)
        sd = cb.duct_htc_static_state_dY(float(state_in.Tt), float(state_in.Pt), G, X, s)
        dG_dm = (1.0 if m_dot >= 0.0 else -1.0) / area
        m_key = f"{self.id}.m_dot"
        t_key, p_key = f"{self.from_node}.Tt", f"{self.from_node}.Pt"

        def chain(d_dG: float, d_dTt: float, d_dPt: float, d_dY: Any) -> dict[str, float]:
            out = {m_key: d_dG * dG_dm, t_key: d_dTt, p_key: d_dPt}
            for i, v in enumerate(d_dY):
                if v != 0.0:
                    out[f"{self.from_node}.Y[{i}]"] = float(v)
            return out

        return _SurfaceFeed(
            T=float(s.T),
            P=float(s.P),
            rho=float(s.rho),
            drho=chain(s.drho_dG, s.drho_dTt, s.drho_dPt, sd.drho_dY),
            dT=chain(s.dT_dG, s.dT_dTt, s.dT_dPt, sd.dT_dY),
            dP=chain(s.dP_dG, s.dP_dTt, s.dP_dPt, sd.dP_dY),
        )

    def _ribbed_residuals(self, state_in: NetworkMixtureState, state_out: NetworkMixtureState):
        """Ribbed channel: the correlation owns the friction factor.

        No multiplier on pipe friction. The `f` the roughness function is
        defined against is already the equivalent four-sided channel value, so
        the correlation returns what the channel needs and the element uses it
        directly. Routing it through a locally restated smooth friction factor
        is what #331 removed, after the round-trip moved a drop by -24.7%.

        The Jacobian is simpler than it looks. `R` carries no `e+` term in this
        correlation family, so **f does not depend on Reynolds number** and
        therefore not on mass flow: `df/d(mdot) = 0`. The drop is
        `f (L/D) rho v |v| / 2` with `v = mdot / (rho A)`, so

            dP          = f (L/D) mdot |mdot| / (2 rho A^2)
            d(dP)/dmdot = f (L/D) |mdot| / (rho A^2)

        which is even in `mdot` and so does not flip sign at zero -- the drop
        itself carries the sign. That is the odd-quantity case: magnitude from
        the correlation, direction from the flow.

        `dP` also depends on upstream temperature and static pressure through
        `rho(T, P)` -- `d(dP)/d(rho) = -dP/rho`, chained through
        `density_and_jacobians`' analytic `d(rho)/dT`, `d(rho)/dP` and
        `_safe_rho`'s own floor derivative. This was a real, measured gap
        (issue #378): the analytic `A.T` column was exactly 0.0 against a
        finite-difference value of roughly -486 for a representative case,
        because `_safe_rho(state_in.density())` computed `rho` without ever
        differentiating it. Composition (`Y`) sensitivity through molecular
        weight remains an accepted gap -- no correlation-level Jacobian for
        it exists yet, and no other channel-friction path in this file
        exposes one either.

        Where the channel leaves a node at rest (a plenum or a boundary
        reservoir), rho and mu are the duct's own static state there, not the
        node's (#926, ``_surface_rest_feed``): rho then moves with m_dot and
        the node's Tt and Pt (and its Y, through the static state).
        """
        m_dot = state_in.m_dot
        model = self.surface.model
        f_mult = self.surface.f_multiplier

        area = self.area or 0.0
        feed = self._surface_rest_feed(state_in, m_dot, area)
        if feed is None:
            rho_raw, drho_raw_dT, drho_raw_dP = _solver_tools.density_and_jacobians(
                state_in.Tt, state_in.P, state_in.X
            )
            rho, drho_draw = _safe_rho(rho_raw)
            drho = {
                f"{self.from_node}.Tt": drho_draw * drho_raw_dT,
                f"{self.from_node}.P": drho_draw * drho_raw_dP,
            }
        else:
            # The duct's own static state at a reservoir end (#926).
            rho, drho_draw = _safe_rho(feed.rho)
            drho = {k: drho_draw * v for k, v in feed.drho.items()}
        if area <= 0.0:
            return [state_in.Pt - state_out.Pt], {
                0: {
                    f"{self.from_node}.Pt": 1.0,
                    f"{self.to_node}.Pt": -1.0,
                }
            }

        dh = self.Dh or self.diameter or 1.0
        if feed is None:
            mu = cb.complete_state(state_in.Tt, state_in.P, state_in.X).transport.mu
        else:
            mu = cb.complete_state(feed.T, feed.P, state_in.X).transport.mu
        Re = abs(m_dot) * dh / (area * mu) if mu > 0.0 else 0.0
        Re = math.copysign(Re, m_dot)

        rib = cb.evaluate_rib(
            model.correlation_set or cb.han_1988_orthogonal(),
            cb.RibGeometry(
                e_D=model.e_D,
                p_e=model.p_e,
                W_H=model.W_H,
                alpha_deg=model.alpha_deg,
            ),
            Re,
        )
        f = rib.f * f_mult

        coeff = f * (self.length / dh) / (2.0 * rho * area * area)
        dP = coeff * m_dot * abs(m_dot)
        d_dP_d_mdot = 2.0 * coeff * abs(m_dot)
        q_exit, dq_dmdot = self._exit_head(m_dot, rho, area)
        dP += q_exit
        d_dP_d_mdot += dq_dmdot
        # dP = K * mdot|mdot| / rho for K independent of rho (the exit head
        # too), so d(dP)/d(rho) = -dP/rho.
        d_dP_drho = -dP / rho if rho > 0.0 else 0.0

        res = [state_in.Pt - state_out.Pt - dP]
        row: dict[str, float] = {
            f"{self.id}.m_dot": -d_dP_d_mdot,
            f"{self.from_node}.Pt": 1.0,
            f"{self.to_node}.Pt": -1.0,
        }
        for k, v in drho.items():
            d = d_dP_drho * v
            row[k] = row[k] - d if k in row else -d
        return res, {0: row}

    def _pin_fin_residuals(self, state_in: NetworkMixtureState, state_out: NetworkMixtureState):
        """Pin-fin array: the correlation owns the drop, per row.

            dP = 2 rho Vmax^2 N f(Re_D),   Vmax = mdot / (rho A_min)
               = 2 N f mdot |mdot| / (rho A_min^2)
            Re_D = mdot D / (mu A_min)                          (signed)

        with ``A_min = A_channel * (A_min/A_frontal)`` from the array geometry
        and ``f`` the canonical per-row factor of ``pin_fin_correlation.h``.
        Unlike ribs, pin friction DOES depend on Re, so

            d(dP)/dmdot = 2N/(rho A_min^2) (2 |mdot| f + mdot |mdot| f' D/(mu A_min))

        which is even in mdot, as the ribbed path's. Upstream T and P enter
        through rho (``-dP/rho``) and through mu in Re
        (``dP/df * f' * (-Re/mu) * dmu``), both from analytic or tabulated
        property derivatives (``density_and_jacobians``,
        ``viscosity_and_jacobians``). Composition sensitivity is the same
        accepted gap as every other channel-friction path here.

        Where the channel leaves a node at rest (a plenum or a boundary
        reservoir), rho and mu are the duct's own static state there, not the
        node's (#926, ``_surface_rest_feed``), chained onto m_dot and the
        node's Tt, Pt and Y (mu's direct composition term stays the gap).
        """
        m_dot = state_in.m_dot
        model = self.surface.model
        f_mult = self.surface.f_multiplier

        area = self.area or 0.0
        feed = self._surface_rest_feed(state_in, m_dot, area)
        if feed is None:
            rho_raw, drho_raw_dT, drho_raw_dP = _solver_tools.density_and_jacobians(
                state_in.Tt, state_in.P, state_in.X
            )
            rho, drho_draw = _safe_rho(rho_raw)
            drho = {
                f"{self.from_node}.Tt": drho_draw * drho_raw_dT,
                f"{self.from_node}.P": drho_draw * drho_raw_dP,
            }
        else:
            # The duct's own static state at a reservoir end (#926).
            rho, drho_draw = _safe_rho(feed.rho)
            drho = {k: drho_draw * v for k, v in feed.drho.items()}
        if area <= 0.0:
            return [state_in.Pt - state_out.Pt], {
                0: {
                    f"{self.from_node}.Pt": 1.0,
                    f"{self.to_node}.Pt": -1.0,
                }
            }

        if feed is None:
            mu, dmu_dT, dmu_dP = _solver_tools.viscosity_and_jacobians(
                state_in.Tt, state_in.P, state_in.X
            )
            dmu = {f"{self.from_node}.Tt": dmu_dT, f"{self.from_node}.P": dmu_dP}
        else:
            mu, dmu_dT, dmu_dP = _solver_tools.viscosity_and_jacobians(feed.T, feed.P, state_in.X)
            # mu at the duct's static (T, P), both functions of the feed's
            # Tt, Pt, the flux and the composition. mu's own composition
            # dependence at fixed (T, P) is the accepted gap below.
            dmu = {
                k: dmu_dT * feed.dT.get(k, 0.0) + dmu_dP * feed.dP.get(k, 0.0)
                for k in dict.fromkeys([*feed.dT, *feed.dP])
            }
        D = model.pin_diameter
        a_min = area * cb.pin_fin_amin_over_afrontal(model.geometry())
        Re = m_dot * D / (a_min * mu) if mu > 0.0 else 0.0

        _, _, f, df_dRe, _ = self.surface._pin_fin_terms(Re, 0.7)
        f *= f_mult
        df_dRe *= f_mult

        coeff = 2.0 * model.N_rows / (rho * a_min * a_min)
        dP = coeff * f * m_dot * abs(m_dot)
        dRe_dmdot = D / (a_min * mu) if mu > 0.0 else 0.0
        d_dP_d_mdot = coeff * (2.0 * abs(m_dot) * f + m_dot * abs(m_dot) * df_dRe * dRe_dmdot)
        # The exit head at the channel's own area, not the pins' A_min.
        q_exit, dq_dmdot = self._exit_head(m_dot, rho, area)
        dP += q_exit
        d_dP_d_mdot += dq_dmdot

        # Through rho (dP ~ 1/rho at fixed f, the exit head too) and through
        # mu in Re (the pin friction only).
        d_dP_drho = -dP / rho if rho > 0.0 else 0.0
        d_dP_dmu = coeff * m_dot * abs(m_dot) * df_dRe * (-Re / mu) if mu > 0.0 else 0.0

        res = [state_in.Pt - state_out.Pt - dP]
        row: dict[str, float] = {
            f"{self.id}.m_dot": -d_dP_d_mdot,
            f"{self.from_node}.Pt": 1.0,
            f"{self.to_node}.Pt": -1.0,
        }
        for k in dict.fromkeys([*drho, *dmu]):
            d = d_dP_drho * drho.get(k, 0.0) + d_dP_dmu * dmu.get(k, 0.0)
            row[k] = row[k] - d if k in row else -d
        return res, {0: row}

    @property
    def residual_scale_kind(self) -> str:
        """Units of this element's residual row (RESIDUAL-ROW-SCALING, see
        NetworkSolver._build_residual_scales): "mdot" for the compressible
        mass-flow form m - m_calc (#481), "p" for the pressure-drop forms
        (incompressible, ribbed, pin-fin). Must change with the residual."""
        if self.regime == "compressible" and not (
            self.surface and isinstance(self.surface.model, (RibbedModel, PinFinModel))
        ):
            return "mdot"
        return "p"

    def _exit_head(self, m_dot: float, rho: float, area: float) -> tuple[float, float]:
        """(q, dq/dm_dot) of the exit dynamic head this duct loses, else zeros.

        q = m|m| / (2 rho A^2): lost where the flow discharges into a pressure
        boundary taking it at its static pressure (#360), as the compressible
        Fanno form and the plain incompressible path do (#588). Odd in m, its
        slope even; it scales as 1/rho, so a caller's d(dP)/d(rho) = -dP/rho
        stays exact with q included.
        """
        if area <= 0.0 or rho <= 0.0 or not self._exit_head_lost(m_dot):
            return 0.0, 0.0
        inv = 1.0 / (rho * area * area)
        return 0.5 * m_dot * abs(m_dot) * inv, abs(m_dot) * inv

    def _exit_head_lost(self, m_dot: float) -> bool:
        """Whether this channel's exit dynamic head is lost where it discharges.

        The exit is the node the flow goes INTO: to_node forward, from_node
        when reversed. Only a PressureBoundary declares a coupling; anything
        else (a plenum, a momentum chamber, a junction port) carries its own
        constitutive relation and receives the stagnation pressure. See #360.
        """
        attr = "_downstream_node" if m_dot >= 0.0 else "_upstream_node"
        node = getattr(self, attr, None)
        if node is None or not hasattr(node, "exit_head_lost"):
            return False
        return bool(node.exit_head_lost(True))

    def _compressible_flow(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        f_mult: float,
        derivatives: bool = True,
    ) -> tuple[float, dict[str, float], Any]:
        """Fanno mass flow m_calc and d(m_calc)/d(node unknowns) (#481).

        The duct is fed isentropically from the stagnation state of the node
        the flow comes FROM and discharges against the other node's
        stagnation pressure, matched to its exit static pressure (head lost)
        or exit stagnation pressure (recovered). The flux rises monotonically
        as that pressure falls and saturates at the choked flux: it exists for
        every state, so there is no infeasible region and no barrier.

        derivatives=False returns m_calc only (empty derivs): the kernel then
        skips its implicit derivatives. With derivatives, d/dY of the feeding
        node is included where its composition moves with the solve: the
        solver sets _composition_directions -- per such node, an orthonormal
        basis of the directions it moves in -- before it reads a Jacobian.
        The kernel differentiates along that basis and the gradient is
        projected back onto per-species keys, exact for every relay column.
        Standalone (unset): per species present, other fractions fixed.
        """
        A = self.area
        forward = state_in.Pt >= state_out.Pt
        src, dst = (state_in, state_out) if forward else (state_out, state_in)
        src_id, dst_id = (
            (self.from_node, self.to_node) if forward else (self.to_node, self.from_node)
        )
        lost = self._exit_head_lost(1.0 if forward else -1.0)
        known = getattr(self, "_composition_directions", None)
        if not derivatives:
            dirs: list[Any] = []
        elif known is None:
            n_sp = len(src.Y)
            dirs = [[float(j == k) for j in range(n_sp)] for k, y in enumerate(src.Y) if y > 0.0]
        else:
            dirs = known.get(src_id, [])
        # Wall heat in the momentum (#551, #672): Tt rises linearly from the
        # feeding node's to this duct's exit Tt.
        exit_state = (
            None
            if getattr(self, "_rayleigh_suspended", False)
            else self._rayleigh_exit_state(dst_id, dst, src_id, src)
        )
        heated = exit_state is not None
        if heated:
            exit_Tt, exit_keys = exit_state
            flow = cb.fanno_rayleigh_channel_flow(
                src.Pt,
                src.Tt,
                exit_Tt,
                src.X,
                dst.Pt,
                not lost,
                self.length,
                self.diameter,
                self.roughness,
                self.friction_model,
                f_mult,
                getattr(self, "_last_G", -1.0),
                [list(v) for v in dirs],
                derivatives,
            )
            # Warm start only: the flow does not depend on it.
            if not flow.choked:
                self._last_G = float(flow.G)
        else:
            flow = cb.fanno_channel_flow(
                src.Pt,
                src.Tt,
                src.X,
                dst.Pt,
                not lost,
                self.length,
                self.diameter,
                self.roughness,
                self.friction_model,
                f_mult,
                getattr(self, "_last_M_exit", -1.0),
                [list(v) for v in dirs],
                derivatives,
            )
            # Warm start only: the flow does not depend on it. The root's exit
            # Mach, which below the drive floor is the floor's, not the
            # reported M_exit at the blended flux (#1017).
            if not flow.choked:
                self._last_M_exit = float(flow.M_exit_root)
        sign = 1.0 if forward else -1.0
        if not derivatives:
            return sign * A * flow.G, {}, flow
        derivs = {
            f"{src_id}.Pt": sign * A * flow.dG_dPt0,
            f"{src_id}.Tt": sign * A * flow.dG_dTt0,
            f"{dst_id}.Pt": sign * A * flow.dG_dP_target,
        }
        if heated:
            for key, coeff in exit_keys:
                derivs[key] = derivs.get(key, 0.0) + sign * A * flow.dG_dTt_exit * coeff
        # Composition of the feeding node: the directional derivatives
        # projected back onto species, g = sum_i (dG/dv_i) v_i, chained by the
        # solver with dY/dx from its mixing relay.
        if dirs:
            g = [0.0] * len(src.Y)
            for d, v in zip(flow.dG_ddir, dirs, strict=True):
                for k, vk in enumerate(v):
                    g[k] += d * float(vk)
            for k, gk in enumerate(g):
                if gk != 0.0:
                    derivs[f"{src_id}.Y[{k}]"] = sign * A * gk
        return sign * A * flow.G, derivs, flow

    def choked_discharge_into(
        self, nid: str, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> str:
        """The empty string when this channel discharges CHOKED into node ``nid`` -- its
        Fanno flow at the choked flux, independent of the pressure it
        discharges against -- else why not (#788, a choked-bounded junction).
        Only the plain compressible channel: a crossflow segment solves
        between stations, and a ribbed or pin-fin channel's drop is a
        correlation with no choke."""
        if type(self) is not ChannelElement:
            return f"'{self.id}' ({type(self).__name__}) feeds it, not a plain channel"
        if self.regime != "compressible":
            return f"'{self.id}' is incompressible"
        if self.surface and isinstance(self.surface.model, (RibbedModel, PinFinModel)):
            return f"'{self.id}' carries a ribbed or pin-fin drop, not a Fanno flow"
        f_mult = self.surface.f_multiplier if self.surface else 1.0
        _, _, flow = self._compressible_flow(state_in, state_out, f_mult, derivatives=False)
        into = self.to_node if state_in.Pt >= state_out.Pt else self.from_node
        if into != nid:
            return f"'{self.id}' does not discharge into it"
        if not bool(flow.choked):
            return f"'{self.id}' is not choked"
        return ""

    def residuals(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> list[float]:

        m_dot = state_in.m_dot

        # User-set empirical correction on f. Not a correlation: it encodes
        # nothing and defaults to 1.0, so it is inert unless a caller reaches
        # for it. Correlation-derived multipliers were removed in 0.7.0; this
        # one is a tuning knob and stays. See issue #339.
        f_mult = self.surface.f_multiplier if self.surface else 1.0

        if self.surface and isinstance(self.surface.model, RibbedModel):
            return self._ribbed_residuals(state_in, state_out)
        if self.surface and isinstance(self.surface.model, PinFinModel):
            return self._pin_fin_residuals(state_in, state_out)

        if self.regime == "compressible":
            # Mass-flow form, like the compressible orifice: m - m_calc, with
            # m_calc the Fanno flow the end states drive (#481). The former
            # drop-given-m form had no physical drop past choke and patched
            # one in with a barrier, which converged choked ducts 5-20% above
            # their choked flow.
            # The residual needs the flow only; its implicit derivatives are
            # computed if and when the Jacobian is read (LazyJacobian, #489).
            m_calc, _, _ = self._compressible_flow(state_in, state_out, f_mult, derivatives=False)

            def jacobian() -> dict[int, dict[str, float]]:
                _, derivs, _ = self._compressible_flow(state_in, state_out, f_mult)
                jac_c: dict[str, float] = {f"{self.id}.m_dot": 1.0}
                for name, d in derivs.items():
                    jac_c[name] = jac_c.get(name, 0.0) - d
                return {0: jac_c}

            return [m_dot - m_calc], LazyJacobian(jacobian)
        else:
            # Use incompressible Darcy-Weisbach formulation. Density
            # reference: upstream static by default, downstream static when
            # the solver's compressible-seed proxy sets
            # _incompressible_p_ref = "outlet" (see OrificeElement).
            # The exit dynamic head follows the same sink coupling as the
            # compressible form (#360): lost where the flow discharges into a
            # pressure boundary taking it at its static pressure. Without it
            # the two regimes disagreed by O(1), not O(Ma^2), at Ma -> 0 (#588).
            return self._incompressible_residuals(state_in, state_out, f_mult)

    def _incompressible_drop(
        self,
        src: NetworkMixtureState,
        dst: NetworkMixtureState,
        src_id: str,
        dst_id: str,
        m_dot: float,
        f_mult: float,
    ) -> tuple[float, float, dict[str, float]]:
        """(dP, d dP/dm, {name: d dP/d name}) with ``src`` as the upstream node.

        Darcy-Weisbach at the upstream static density (the downstream static
        under the solver's compressible-seed proxy, _incompressible_p_ref =
        "outlet"), plus the exit head where the flow discharges into a static
        sink (#360, #588). dP is odd in m_dot.
        """
        p_ref_outlet = getattr(self, "_incompressible_p_ref", "inlet") == "outlet"
        if p_ref_outlet:
            p_rho, dp_rho, dp_rho_dm = dst.P, {f"{dst_id}.P": 1.0}, 0.0
        else:
            # The duct's own inlet static state (#651): the node's (Pt, T) at
            # the duct's mass flux. A PlenumNode reports P = Pt by definition,
            # so its static read the duct's density a dynamic head high (3% in
            # friction at Mach 0.18). A node with the duct's area already
            # carries this state.
            A = math.pi * self.diameter * self.diameter / 4.0
            a = _solver_tools.approach_static_pressure(
                float(src.Pt), float(src.Tt), abs(m_dot) / A, src.Y, True
            )
            p_rho = float(a.P)
            dp_rho = {f"{src_id}.Pt": float(a.dP_dPt), f"{src_id}.Tt": float(a.dP_dTt)}
            for i, val in enumerate(a.dP_dY):
                if val != 0.0:
                    dp_rho[f"{src_id}.Y[{i}]"] = float(val)
            dp_rho_dm = float(a.dP_dG) * (1.0 if m_dot >= 0.0 else -1.0) / A
        r = _solver_tools.channel_residuals_and_jacobian(
            m_dot,
            src.Pt,
            p_rho,
            src.Tt,
            src.Y,
            dst.P,
            self.length,
            self.diameter,
            self.roughness,
            self.friction_model,
            f_mult,
            self._exit_head_lost(m_dot),
        )
        g = float(r.d_dP_dP_static_up)
        d: dict[str, float] = {f"{src_id}.Tt": float(r.d_dP_dTt_up)}
        for i, val in enumerate(r.d_dP_dY_up):
            if val != 0.0:
                d[f"{src_id}.Y[{i}]"] = float(val)
        for name, dp in dp_rho.items():
            d[name] = d.get(name, 0.0) + g * dp
        return float(r.dP_calc), float(r.d_dP_d_mdot) + g * dp_rho_dm, d

    def _incompressible_residuals(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        f_mult: float,
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        """Pt_from - Pt_to = dP(m), the drop FOLLOWING THE FLOW (#554).

        The density (and composition, temperature) is the upstream node's
        along the flow, which for a reversed channel is the declared outlet;
        reading the declared inlet regardless made a reversed duct pass the
        wrong flow. Forward and reverse forms blend through zero flow over
        FLOW_REVERSAL_M_SCALE, as the area change and the pressure loss do.

        A channel connects stagnation pressures: Pt_up - Pt_down = dP_friction
        (+ q_exit at a static sink). The static-vs-stagnation distinction lives
        in the node closures (PlenumNode: Pt = P; MomentumChamberNode: Pt = P +
        0.5 rho v^2); referencing the downstream static instead would let an
        inline chamber leak its dynamic head as a free pressure gain.
        """
        m_dot = state_in.m_dot
        w, dw = _flow_blend(m_dot)
        if w == 1.0:
            dP, dP_dm, dd = self._incompressible_drop(
                state_in, state_out, self.from_node, self.to_node, m_dot, f_mult
            )
        elif w == 0.0:
            dP, dP_dm, dd = self._incompressible_drop(
                state_out, state_in, self.to_node, self.from_node, m_dot, f_mult
            )
        else:
            f_dP, f_dm, f_d = self._incompressible_drop(
                state_in, state_out, self.from_node, self.to_node, m_dot, f_mult
            )
            r_dP, r_dm, r_d = self._incompressible_drop(
                state_out, state_in, self.to_node, self.from_node, m_dot, f_mult
            )
            dP = w * f_dP + (1.0 - w) * r_dP
            dP_dm = w * f_dm + (1.0 - w) * r_dm + dw * (f_dP - r_dP)
            dd = {
                k: w * f_d.get(k, 0.0) + (1.0 - w) * r_d.get(k, 0.0)
                for k in dict.fromkeys([*f_d, *r_d])
            }
        res = [state_in.Pt - state_out.Pt - dP]
        row: dict[str, float] = {
            f"{self.id}.m_dot": -dP_dm,
            f"{self.from_node}.Pt": 1.0,
            f"{self.to_node}.Pt": -1.0,
        }
        for k, v in dd.items():
            row[k] = row.get(k, 0.0) - v
        return res, {0: row}

    def _rest_end_states(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flow: Any,
    ) -> list[tuple[float, float, list[float]] | None]:
        """The duct's own static (T, P, X) at each end on a node at rest
        (#888): [from end, to end], None where the end keeps its node state.

        A node whose own definition holds its gas at rest (``state_at_rest``,
        Pt = P: a plenum, a combustor, a PressureBoundary or MassFlowBoundary)
        reports a state that is not the duct's: a reservoir (P = Pt, T = Tt)
        where the duct leaves it, a static pressure with the arriving Tt where
        it enters. Read as the duct's end state it put mach_in at the
        reservoir density (0.58 where the duct carries 0.89) and mach_out
        above 1 in an unchoked duct. At such an end the diagnostics read the
        state the duct's own residual solves instead:

        * compressible: the Fanno station at that end's Mach (M_in at the
          kernel's source, M_exit at its exit) for the kernel's flux and the
          duct's stagnation temperature there (its exit Tt when the wall heat
          is in the momentum), so the reported Mach is the kernel's;
        * incompressible: the residual's own density state (#651), (Tt, the
          approach static P at the flux) where the node holds a total
          pressure; where a boundary takes the flow at its static pressure
          (an outflow losing its head) its state is the duct face's.

        * ribbed or pin-fin surface: the from end's state its own drop reads
          (``_surface_rest_feed``, #926); the to end keeps the node state.

        An end on a moving node (a momentum chamber, Pt > P) carries its own
        static state and is left alone. Reporting only: no residual reads
        this.
        """
        ends: list[tuple[float, float, list[float]] | None] = [None, None]
        nodes = (getattr(self, "_upstream_node", None), getattr(self, "_downstream_node", None))
        # A choked-bounded junction (#788) is reported as a plenum (at rest):
        # that end reads the duct's own Fanno exit, as at a node at rest.
        at_rest = [
            bool(getattr(n, "state_at_rest", False))
            or getattr(n, "_choked_bounded", None) is not None
            for n in nodes
        ]
        m = float(state_in.m_dot)
        area = self.area or 0.0
        if not self._rest_end_state_from_duct or not any(at_rest):
            return ends
        if area <= 0.0 or m == 0.0:
            return ends
        if self.surface is not None and isinstance(self.surface.model, (RibbedModel, PinFinModel)):
            # The state its drop reads at the from end (#926); the to end,
            # which its residual does not read, keeps the node state.
            feed = self._surface_rest_feed(state_in, m, area)
            if feed is not None:
                ends[0] = (feed.T, feed.P, list(state_in.X))
            return ends
        if self.regime == "compressible":
            if flow is None or not float(flow.G) > 0.0:
                return ends
            # The kernel's source is the higher-Pt end (_compressible_flow).
            forward = state_in.Pt >= state_out.Pt
            src, dst = (state_in, state_out) if forward else (state_out, state_in)
            src_id, dst_id = (
                (self.from_node, self.to_node) if forward else (self.to_node, self.from_node)
            )
            exit_Tt = float(src.Tt)
            if not getattr(self, "_rayleigh_suspended", False):
                rx = self._rayleigh_exit_state(dst_id, dst, src_id, src)
                if rx is not None:
                    exit_Tt = float(rx[0])
            X = list(src.X)
            for i in (0, 1):
                if not at_rest[i]:
                    continue
                at_src = (i == 0) == forward
                M, Tt = (
                    (float(flow.M_in), float(src.Tt)) if at_src else (float(flow.M_exit), exit_Tt)
                )
                if not 0.0 < M <= 1.0:
                    continue
                st = cb.fanno_state_at_mach(float(flow.G), Tt, M, X)
                ends[i] = (float(st.T), float(st.P), X)
            return ends
        G = abs(m) / area
        for i, st in enumerate((state_in, state_out)):
            if not at_rest[i]:
                continue
            outflow = (i == 1) == (m >= 0.0)
            node = nodes[i]
            if hasattr(node, "exit_head_lost") and bool(node.exit_head_lost(outflow)):
                continue
            a = _solver_tools.approach_static_pressure(float(st.Pt), float(st.Tt), G, st.Y, True)
            ends[i] = (float(st.Tt), float(a.P), list(st.X))
        return ends

    def diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float | str]:
        # The duct's own Fanno solution (compressible): its end states at a
        # node at rest (#888) and the M_in/M_exit_duct keys below.
        flow = None
        if self.regime == "compressible":
            f_mult = self.surface.f_multiplier if self.surface else 1.0
            _, _, flow = self._compressible_flow(state_in, state_out, f_mult, derivatives=False)
        end_in, end_out = self._rest_end_states(state_in, state_out, flow)
        cs_in = (
            cb.complete_state(*end_in)
            if end_in is not None
            else cb.complete_state(state_in.T, state_in.P, state_in.X)
        )
        ref = _element_reference_block(
            cs_in,
            m_dot=state_in.m_dot,
            area=self.area,
            Dh=self.Dh or self.diameter,
            location="inlet",
        )
        v_in = ref["velocity"]
        re_in = ref["Re"]
        mach_in = v_in / cs_in.thermo.a if cs_in.thermo.a > 0 and v_in > 0 else 0.0

        cs_out = (
            cb.complete_state(*end_out)
            if end_out is not None
            else cb.complete_state(state_out.T, state_out.P, state_out.X)
        )
        rho_out = cs_out.thermo.rho
        v_out = (
            (abs(state_out.m_dot)) / ((rho_out) * (self.area))
            if self.area > 0 and rho_out > 0
            else 0.0
        )
        mach_out = v_out / cs_out.thermo.a if cs_out.thermo.a > 0 and v_out > 0 else 0.0

        # Already at the duct's own inlet static state (#636), not the node's.
        h_res = self.htc_and_T(state_in)
        surface_diag: dict[str, float] = {}
        enhanced = self.surface is not None and isinstance(
            self.surface.model, (RibbedModel, PinFinModel)
        )
        if h_res is not None:
            Nu, htc, T_aw = h_res.Nu, h_res.h, h_res.T_aw
            # The surface correlation's own Re (jet Re_j for impingement, pin
            # Re_D for pin fins) and whether it is outside its validity box.
            surface_diag["Re_surface"] = float(h_res.Re)
            surface_diag["surface_extrapolated"] = float(
                bool(getattr(h_res, "extrapolated", False))
            )
            # The heat-transfer correlation's own f, at its own Re (the duct's
            # static state, #636). It drives the pressure drop only on an
            # enhanced (ribbed, pin-fin) surface, whose residual uses it.
            surface_diag["f_surface"] = float(h_res.f)
        else:
            Nu, htc, T_aw = 0.0, 0.0, state_in.Tt
            # No correlation evaluated (no flow): the keys stay, undefined.
            surface_diag = {
                "Re_surface": math.nan,
                "surface_extrapolated": math.nan,
                "f_surface": math.nan,
            }
        if enhanced and h_res is not None:
            f = h_res.f
        else:
            # Report the friction factor the residual actually uses: C++'s
            # laminar-turbulent blend at the inlet Re, Poiseuille 64/Re with
            # its own re_floor, times the user's f multiplier (#634). Every
            # channel carries a convective surface, so this used to report
            # the heat-transfer correlation's f at ITS Re (16% apart at Re
            # 70), and its old floor sqrt(Re^2 + 64^2) gave f = 1 at Re = 0
            # where the row uses 64/Re.
            dh_diag = self.Dh or self.diameter or 1.0
            e_D = self.roughness / dh_diag if dh_diag > 0 else 0.0
            f_mult = self.surface.f_multiplier if self.surface else 1.0
            f = f_mult * cb._core.friction_and_jacobian(self.friction_model, re_in, e_D).result[0]

        duct: dict[str, float] = {}
        if flow is not None:
            # The duct's own Fanno solution: inlet and exit Mach inside the
            # duct (mach_in / mach_out at an end on a node at rest), and
            # whether it is choked -- its flow then independent of the back
            # pressure.
            duct = {
                "choked": float(flow.choked),
                "M_in_duct": float(flow.M_in),
                "M_exit_duct": float(flow.M_exit),
            }
            if getattr(self, "_wall_heated", False):
                into = self.to_node if float(state_in.m_dot) >= 0.0 else self.from_node
                into_state, from_id, from_state = (
                    (state_out, self.from_node, state_in)
                    if into == self.to_node
                    else (state_in, self.to_node, state_out)
                )
                if self._rayleigh_exit_state(into, into_state, from_id, from_state) is not None:
                    # The heat is in the momentum: nothing omitted (#551, #672).
                    duct["rayleigh_in_momentum"] = 1.0
                else:
                    duct.update(self._rayleigh_omission(state_in, state_out, flow, ref))
        # Heat-transfer area beyond the friction-wetted area (#551). Ribbed and
        # pin-fin surfaces carry their own friction, so only a smooth wall.
        if self.surface is not None and isinstance(self.surface.model, SmoothModel):
            dh_w = self.Dh or self.diameter
            if dh_w and self.area and self.length:
                wetted = 4.0 * self.area * self.length / dh_w
                ratio = float(self.surface.area) / wetted
                duct["heat_area_ratio"] = ratio
                duct["heat_area_exceeds_wetted"] = float(ratio > HEAT_AREA_EXCESS_FLAG)
        return {
            "m_dot": float(state_in.m_dot),
            **_element_pressure_block(state_in, state_out, mach_in=mach_in, mach_out=mach_out),
            **ref,
            **duct,
            "Dh": float(self.Dh or self.diameter or 0.0),
            "diameter": float(self.diameter or 0.0),
            "area_source": self._area_source,
            "Nu": float(Nu),
            "St": _stanton(h_res),
            "htc": float(htc),
            "T_aw": float(T_aw),
            "f": float(f),
            **surface_diag,
        }

    def _rayleigh_exit_state(
        self,
        dst_id: str,
        dst: NetworkMixtureState,
        src_id: str,
        src: NetworkMixtureState,
    ) -> tuple[float, list[tuple[str, float]]] | None:
        """This wall-heated duct's exit Tt into ``dst_id``, with the relay keys
        its derivatives chain through (key, coefficient), or ``None``
        (adiabatic Fanno).

        Into a plenum or chamber fed by this duct alone, with no wall of its
        own, that node's Tt IS the exit Tt (#551). Into any other node (one
        mixing other streams, one with its own wall, a boundary) the solver
        derives the exit state from this duct's stream and its own wall heat
        alone, through the same mixer and relay (#672): ``_rayleigh_exit`` =
        (node the flow enters, exit Tt, node it comes from), keyed
        ``RAYLEIGH_EXIT_PREFIX + id``. The solver sets it where the flow goes
        in by its SIGN; the kernel's exit is the lower-Pt end. Where the two
        disagree (an iterate on its way through reversal) the duct's own
        temperature rise is applied to the kernel's source, Tt_src + (Tt_exit
        - Tt_from): without it the duct dropped to adiabatic there, a jump
        in its residual that cost a stress draw its cold solve.
        """
        if dst_id in getattr(self, "_rayleigh_ends", ()):
            return (float(dst.Tt), [(f"{dst_id}.Tt", 1.0)]) if float(dst.Tt) > 0.0 else None
        if dst_id not in getattr(self, "_rayleigh_mix_ends", ()):
            return None
        rx = getattr(self, "_rayleigh_exit", None)
        if rx is None or not rx[1] > 0.0:
            return None
        key = f"{RAYLEIGH_EXIT_PREFIX}{self.id}.Tt"
        if rx[0] == dst_id:
            return float(rx[1]), [(key, 1.0)]
        # Flow sign and Pt disagree: rx[0] is this kernel's source, rx[2] its exit.
        if rx[0] != src_id or rx[2] != dst_id:
            return None
        exit_Tt = float(src.Tt) + float(rx[1]) - float(dst.Tt)
        if not exit_Tt > 0.0:
            return None
        return exit_Tt, [(key, 1.0), (f"{dst_id}.Tt", -1.0), (f"{src_id}.Tt", 1.0)]

    def _rayleigh_omission(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flow: Any,
        ref: dict[str, float],
    ) -> dict[str, float]:
        """The total-pressure cost of the wall heat this duct's momentum leaves
        out (#551), as a flag: the compressible kernel is Fanno flow at the
        INLET Tt, and the heat only raises the downstream node's Tt. To first
        order in the heating, Rayleigh flow loses
        dPt / Pt = (gamma M^2 / 2) ln(Tt_exit / Tt_in), here with M^2 the mean
        of the duct's inlet and exit values; cooling gains it. O(M^2): it
        matters at high Mach with strong heating (0.06 of Pt at M 0.25-0.73 and
        Tt ratio 2.3-4), not in a low-Mach heat exchanger."""
        m = float(state_in.m_dot)
        up, dn = (state_in, state_out) if m >= 0.0 else (state_out, state_in)
        if not (float(up.Tt) > 0.0 and float(dn.Tt) > 0.0):
            return {}
        gamma = float(ref.get("gamma_ref", 1.4))
        m2 = 0.5 * (float(flow.M_in) ** 2 + float(flow.M_exit) ** 2)
        frac = 0.5 * gamma * m2 * math.log(float(dn.Tt) / float(up.Tt))
        return {
            "rayleigh_omitted_frac": frac,
            "rayleigh_dPt_omitted": frac * float(up.Pt),
            # Flagged above RAYLEIGH_OMISSION_FLAG of the inlet Pt.
            "rayleigh_omitted": float(abs(frac) > RAYLEIGH_OMISSION_FLAG),
        }

    def get_spatial_profile(
        self,
        state_in: NetworkMixtureState,
        n_steps: int = 100,
    ) -> list:
        """
        Compute and return the spatial flow profile array along the channel length.
        Requires the solved inlet NetworkMixtureState.
        """

        rho = cb.density(state_in.Tt, state_in.P, state_in.X)
        u = cb.channel_velocity(state_in.m_dot, self.diameter, rho)

        if self.regime == "incompressible":
            res = cb.channel_flow_rough(
                state_in.Tt,
                state_in.P,
                state_in.X,
                u,
                self.length,
                self.diameter,
                self.roughness,
                self.friction_model,
                n_steps,
                True,
            )
            return res.profile

        elif self.regime == "compressible":
            # The x-march from the duct's inlet STATIC state, which the Mach
            # solution gives (isentropic from the feeding node's stagnation
            # state). It called fanno_channel with its arguments out of order
            # and raised for every compressible channel (#481).
            G = abs(state_in.m_dot) / self.area
            if G <= 0.0:
                return []
            duct = cb.fanno_duct(
                state_in.Pt,
                state_in.Tt,
                G,
                state_in.X,
                self.length,
                self.diameter,
                self.roughness,
                self.friction_model,
                self.surface.f_multiplier if self.surface else 1.0,
            )
            res = cb.fanno_channel_rough(
                duct.inlet.T,
                duct.inlet.P,
                duct.inlet.u,
                self.length,
                self.diameter,
                self.roughness,
                state_in.X,
                self.friction_model,
                self.surface.f_multiplier if self.surface else 1.0,
                n_steps,
                True,
            )
            return res.profile

        return []

    def htc_and_T(self, state: NetworkMixtureState):
        """Compute heat transfer coefficient and adiabatic wall temperature.

        Uses the element's ConvectiveSurface to compute HTC and adiabatic wall temperature.

        Parameters
        ----------
        state : NetworkMixtureState
            Flow state at the element inlet.

        Returns
        -------
        ChannelResult | None
            Full ChannelResult with h, T_aw, and Jacobians (dh_dmdot, dh_dT, etc.),
            or None if surface.area = 0. Access convective area via ``self.surface.area``.
        """
        if self.surface.area == 0.0:
            return None

        # Use nan for T_hot if not specified (matches C++ default)
        T_hot = self.t_hot if self.t_hot is not None else math.nan

        # At the duct's own static state, not the feed node's (#636).
        return _duct_htc(
            self.surface, state, self.area, self.Dh or self.diameter, self.length, T_hot
        )

    def _seed_resistance(self, rho: float, P_ref: float, m_ref: float) -> float | None:
        """Start guess: a turbulent friction factor, plus the exit head where
        the residual loses it -- the declared exit a static-coupled pressure
        boundary (``_exit_head_lost``, #360). Charged on every segment, a
        chain of N short ducts between plenums started N heads high: the
        D-7156 film rig's ten 1.22 cm segments 14 kPa against a 13 Pa root,
        and its cold solve stalled (#971)."""
        dh = self.Dh or self.diameter
        if not dh or not self.area or not self.length:
            return None
        f = SEED_FRICTION_FACTOR * getattr(self.surface, "f_multiplier", 1.0)
        head = 1.0 if self._seed_exit_head() else 0.0
        return (f * self.length / dh + head) / (2.0 * rho * self.area**2)

    def _seed_exit_head(self) -> bool:
        """Whether the start guess charges one exit head: where the residual
        loses it at the declared exit. Before topology is resolved the exit
        is unknown, and the head is charged, as it always was: an estimate
        made while a network is built (the stress splitter sizes its supply
        from it) keeps its value."""
        if not hasattr(self, "_downstream_node"):
            return True
        return self._exit_head_lost(1.0)

    def _seed_flow_applies(self) -> bool:
        """The compressible mass-flow form (#481): its Fanno kernel saturates
        at the choked flux, which the incompressible K does not."""
        return self.residual_scale_kind == "mdot"

    def _seed_flow(
        self, Pt_up: float, Pt_dn: float, Tt: float, X: list[float], Y: list[float]
    ) -> float | None:
        """The flow the residual's own adiabatic Fanno kernel passes in the
        declared direction (``_compressible_flow``), no derivatives (#971)."""
        f_mult = self.surface.f_multiplier if self.surface else 1.0
        try:
            flow = cb.fanno_channel_flow(
                Pt_up,
                Tt,
                X,
                Pt_dn,
                not self._exit_head_lost(1.0),
                self.length,
                self.diameter,
                self.roughness,
                self.friction_model,
                f_mult,
                -1.0,
                [],
                False,
            )
        except (ValueError, RuntimeError):
            return None
        return self.area * float(flow.G)

    def n_equations(self) -> int:
        return 1

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Cached for the exit-coupling lookup in residuals (#360); done before
        # the early return so it is set even when the diameter is explicit.
        self._downstream_node = graph.nodes.get(self.to_node)
        self._upstream_node = graph.nodes.get(self.from_node)
        # A side of a ThermalWall: its heat is deposited on the node the flow
        # goes into, so the end nodes' Tt ratio is this duct's heating (#551).
        # An adiabatic wall (#970) heats nothing.
        walls = list(getattr(graph, "walls", {}).values())
        self._wall_heated = any(self.id in w.heated_endpoints for w in walls)
        # The ends whose Tt IS this duct's exit Tt when the flow goes into
        # them: a plenum or chamber with this duct and one other connection,
        # and no wall of its own -- it mixes nothing else and receives only
        # this duct's heat (#551). At any other end (a boundary, a node mixing
        # several streams or heat) the solver forms the duct's own exit state
        # (_rayleigh_exit_state, #672). Either way the Fanno-Rayleigh kernel
        # carries the heat in the momentum.
        self._rayleigh_ends: set[str] = set()
        if self._wall_heated:
            on_walls = {e for w in walls for e in w.heated_endpoints}
            for nid in (self.from_node, self.to_node):
                node = graph.nodes.get(nid)
                if not isinstance(node, (PlenumNode, MomentumChamberNode)) or nid in on_walls:
                    continue
                n_conn = len(graph.get_upstream_elements(nid)) + len(
                    graph.get_downstream_elements(nid)
                )
                if n_conn == 2:
                    self._rayleigh_ends.add(nid)
            # Every other end: the solver derives this duct's own exit state
            # there (#672).
            self._rayleigh_mix_ends = {self.from_node, self.to_node} - self._rayleigh_ends
        else:
            self._rayleigh_mix_ends = set()
        self._rayleigh_exit = None
        if self.diameter is not None:
            self._fill_convective_area()
            return
        # Inherit diameter from the nearest geometry source. The channel's
        # area sets its friction and Mach behaviour, so a defaulted diameter
        # silently rewrites the flow; unresolvable is an error (issue #262).
        area, source = _infer_flow_area(graph, [self.from_node, self.to_node], exclude=self)
        if area is not None:
            self.diameter = math.sqrt(4.0 * area / math.pi)
            self._area_source = source
        if self.diameter is None:
            if getattr(self, "_topo_pass", 0) >= 1:
                if self.default_diameter is None:
                    raise _unresolved_area_error(self, "diameter", [self.from_node, self.to_node])
                self.diameter = float(self.default_diameter)
                self._area_source = "default"
            else:
                # Defer: an upstream tee may not have resolved F_C yet this pass.
                self._topo_pass = 1
                return
        self.area = math.pi * (self.diameter / 2) ** 2
        if self.Dh is None:
            self.Dh = self.diameter
        # Update convective surface area (was 0 if diameter was deferred)
        self._fill_convective_area()

    def default_convective_area(self) -> float:
        """The convective area a surface gets when none was set, by model.

        * Impingement array (one row): the row's target footprint matched to
          this channel's crossflow cross-section. The channel is z wide in
          the jet direction, so its span is ``A_flow / z`` and the footprint
          ``xn * span = A_flow * (xn/d) / (z/d)``; the hole count then
          follows from the channel itself (#462).
        * Single jet: the disc of radius ``R`` that Goldstein et al. (1986)
          average Nu over (Han Eq. 4.1 is an AVERAGE Nu out to R/D),
          ``pi (R_D d_jet)^2``.
        * Anything else: the wetted wall, ``pi D L``.
        """
        model = self.surface.model if self.surface else None
        if isinstance(model, ImpingementModel) and self.area and model.z_d > 0.0:
            return self.area * model.xn_d / model.z_d
        if isinstance(model, SingleJetImpingementModel):
            return math.pi * (model.R_D * model.d_jet) ** 2
        return math.pi * (self.diameter or 0.0) * self.length

    def _fill_convective_area(self) -> None:
        if self.surface and self.surface.area == 0.0:
            self.surface.area = self.default_convective_area()

    def _geometry_error(self) -> str:
        """Why this duct's section cannot be built (#1024), or "".

        The diameter (the area's, given or inherited) must be positive and
        finite: a duct of no area is no duct, and the Fanno kernel refused it
        mid-solve. A zero-area restriction is a closed orifice instead. A
        given Dh must be positive and finite.

        NOT the isoperimetric bound the nodes check (Dh <= 2 sqrt(A/pi)): it
        holds for a closed wetted perimeter only. A channel standing for a
        share of a wider duct, its area a fraction of the duct's and only its
        part of the wall wetted (a sector), keeps the duct's Dh = 4A/P above
        it, as the C3X vane example's hole segments do (Dh = D, area f A)."""
        name = f"{type(self).__name__} {self.id!r}"
        if self.diameter is not None and not _positive_finite(self.diameter):
            return (
                f"{name}: diameter D = {self.diameter!r} m. A duct needs a positive, "
                "finite diameter (area); to shut a path, close it with a zero-bore "
                "OrificeElement instead."
            )
        if self._Dh_source == "user" and not _positive_finite(self.Dh):
            return f"{name}: hydraulic diameter Dh = {self.Dh!r} m must be positive and finite."
        return ""

    def validate(self) -> None:
        if error := self._geometry_error():
            raise ValueError(error)


class CrossflowSegmentElement(ChannelElement):
    """A channel segment between side-stream STATIONS: merges and bleeds (#465, #471).

    A ``ChannelElement`` (friction, stagnation-pressure coupling) plus the
    momentum of side streams joining or leaving at its end nodes. Over one
    station, ``m_a`` arriving and ``m_b`` leaving along the channel of area
    A, with the side stream's AXIAL velocity ``kappa * u_a``, momentum gives
    the static-pressure drop

        dP = (m_b^2 - m_a^2 + kappa m_a (m_a - m_b)) / (rho A^2)

    * ``kappa = 0``: normal injection (impingement jets). Florschuetz's
      P + G^2/rho = const.
    * ``kappa = 0.75`` (``cb.STATION_KAPPA_BLEED_BASSETT``): bleed through the
      wall. Bassett, Winterbone & Pearson's (2001) separating straight-run
      coefficient K2/K5 (Eq. 15), theta- and psi-independent, reproduced for
      every split.

    JETS that join at a station bring the axial momentum their element
    declares (#813): ``injection_momentum``, m w cos(theta), the hook a
    merge chamber reads. Over the station, (P_a - P_b) A gains -J. Normal
    jets (90 deg, every impingement row) bring none.

    CENTRED: a segment carries half of the station at each end, each at its
    node's own STATIC density (C++ ``station_half_drop``), so every node holds
    the static pressure at the middle of its station. The node carries T0;
    the density is taken at the static T of the station's mean flow, the
    momentum chamber closure's (#782), the same for both halves. Ends are
    explicit, wired by whoever builds the chain:

    * ``from_kappa`` with ``prev_seg``: a station at ``from_node``; the flow
      arriving there is ``prev_seg``'s (0 if None: a closed end).
    * ``to_kappa`` with ``next_seg``: a station at ``to_node``; the flow
      leaving it is ``next_seg``'s.
    * ``entry_K``: ``from_node`` is a reservoir; the segment carries the
      entry acceleration ``(1 + K_in) m^2/(2 rho A^2)`` at the duct (to-node)
      static density instead of a from-station. Reversed, the duct discharges into
      the reservoir at static continuity, the same as an exit: no term.
    * an exit into a plenum (no to-station) loses the dynamic head: static
      continuity, as any Pt-coupled channel into a plenum. Reversed, that node
      at rest is a reservoir: the gas accelerates out of it (one dynamic head,
      no loss; compressible, the duct end takes its Pt).

    Station nodes must be PlenumNodes (Pt = P): the node pressure is the
    channel's static pressure, which is what a side stream sees.

    TAPER (#812): ``area`` is the from-end's section, ``area_to`` (None: the
    same) the to-end's -- a tapering annulus or a growing liner. Each end's
    station half-drop, joining jets and static-to-total closure read its own
    section; friction runs on the mean. Compressible, the ends' closures carry
    the area change (Pt continuous bar friction); incompressible, a from-end
    station adds Bernoulli's m^2/(2 rho)(1/A_to^2 - 1/A_from^2) (C++
    ``segment_area_change_drop``; a reservoir entry already accelerates to
    the to-end's section). Two segments meeting at a station with different
    sections are a sudden area change at the station, each half at its own.

    FRICTION follows ``ChannelElement``: Darcy-Weisbach on the
    equivalent-area diameter (the C++ channel's circular convention, #463).

    COMPRESSIBLE, the friction is the Fanno flow of the duct between its own
    end states (#802): the station half-drops move the end STATIC pressures
    (inlet P_from - half, outlet P_to + half, or the exit node's P), and the
    kernel runs between those statics' stagnation pressures at the segment's
    flux (``_duct_ends``). The incompressible balance as Ma -> 0 (#585).
    Feeding a station's static P to the kernel as P0, as before #802, put the
    duct's choke a dynamic head too early.
    """

    @property
    def residual_scale_kind(self) -> str:
        """RESIDUAL-ROW-SCALING: the segment branches its residual on the
        channel's regime rule (m - m_calc compressible, a pressure drop
        otherwise), so its units are the channel's."""
        return ChannelElement.residual_scale_kind.fget(self)

    # Its duct runs between its own station states (#802): see ChannelElement.
    _rest_end_state_from_duct = False

    # The channel's area is re-derived from a diameter, which need not give
    # the given area back bit for bit, and a constant segment sets area_to
    # to it (#1099).
    _init_kwarg_attrs = {"area": "_given_area", "area_to": "_given_area_to"}

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        length: float,
        area: float,
        Dh: float | None = None,
        roughness: float = 0.0,
        regime: CompressibilityLiteral = DEFAULT_REGIME,
        friction_model: FrictionModelLiteral = "haaland",
        from_kappa: float | None = None,
        to_kappa: float | None = None,
        prev_seg: str | None = None,
        next_seg: str | None = None,
        entry_K: float | None = None,
        area_to: float | None = None,
    ) -> None:
        if length <= 0.0 or area <= 0.0:
            raise ValueError(f"{type(self).__name__}: length and area must be positive")
        if area_to is not None and not area_to > 0.0:
            raise ValueError(f"{type(self).__name__}: area_to must be positive")
        if entry_K is not None and from_kappa is not None:
            raise ValueError(
                f"{type(self).__name__} {id!r}: from_node is either a reservoir "
                "(entry_K) or a station (from_kappa), not both"
            )
        super().__init__(
            id,
            from_node,
            to_node,
            length=length,
            # A taper's friction runs on the mean section (#812).
            diameter=math.sqrt(
                4.0 * (area if area_to is None else 0.5 * (area + area_to)) / math.pi
            ),
            Dh=Dh,
            roughness=roughness,
            regime=regime,
            friction_model=friction_model,
        )
        # Each end's own section (#812): the stations' half-drops, the jets'
        # momentum and the ends' static-to-total closures read it. A constant
        # segment keeps the one area the channel derived, bit for bit.
        if area_to is None:
            self.area_from = self.area_to = self.area
        else:
            self.area_from = float(area)
            self.area_to = float(area_to)
        self.from_kappa = from_kappa
        self.to_kappa = to_kappa
        self.prev_seg = prev_seg
        self.next_seg = next_seg
        self.entry_K = entry_K
        self._given_area = area
        self._given_area_to = area_to
        self._prev: NetworkElement | None = None
        self._next: NetworkElement | None = None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        super().resolve_topology(graph)
        self._prev = graph.elements.get(self.prev_seg) if self.prev_seg else None
        self._next = graph.elements.get(self.next_seg) if self.next_seg else None
        self._elements = graph.elements

    def station_report(
        self, state: NetworkMixtureState, m_prev: float, m_self: float
    ) -> tuple[str, float, float, float] | None:
        """The moving state of this segment's FROM station, for reporting
        (#782 audit): (node id, static T, Pt, Mach). A station node holds the
        static P at the middle of its station and carries T0; its gas moves at
        the station's mean flow, the one its half-drops use. None without a
        from-station."""
        if self._prev is None or self.from_kappa is None:
            return None
        A = float(self.area_from)
        m = 0.5 * (float(m_prev) + float(m_self))
        incompressible = self.regime == "incompressible"
        st = _moving_static_T(state, m, A, incompressible)
        T, P, X = float(st.T), float(state.P), state.X
        if not (A > 0.0 and T > 0.0 and P > 0.0):
            return None
        rho = float(cb.density(T, P, X))
        u = abs(m) / (rho * A)
        a = float(cb.speed_of_sound(T, X))
        M = u / a if a > 0.0 else 0.0
        Pt = P + 0.5 * rho * u * u if incompressible else float(cb.Pt_from_static(P, T, M, X))
        return self.from_node, T, Pt, M

    def _station_jets(
        self, node: Any, station: NetworkMixtureState
    ) -> tuple[float, dict[str, float]]:
        """Axial momentum the side streams joining at a station bring (#813).

        Each joining element declares it (``injection_momentum``, the same
        hook a merge chamber reads): J = m w cos(theta), the jet's velocity
        from its supply to the station's static pressure. 90 deg (impingement
        jets, normal holes) brings none, so such a station is untouched. A
        segment of the chain is not a side stream.
        """
        J, dJ = 0.0, {}
        elements = getattr(self, "_elements", {})
        for eid, st in getattr(node, "_inflows", ()):
            el = elements.get(eid)
            if el is None or isinstance(el, CrossflowSegmentElement):
                continue
            hook = getattr(el, "injection_momentum", None)
            if hook is None:
                continue
            j, dj = hook(st, station)
            if j == 0.0 and not dj:
                continue
            J += j
            for k, v in dj.items():
                dJ[k] = dJ.get(k, 0.0) + v
        return J, dJ

    def _with_jets(self, term: "_StaticChainedDrop", node: Any, station: Any, area: float) -> None:
        """Fold the station's joining-jet momentum into this segment's half.

        Over the station, (P_a - P_b) A = m_b u_b - m_a u_a - J: the jets'
        axial momentum J lowers the static drop. Centred, each segment
        carries half of it.
        """
        J, dJ = self._station_jets(node, station)
        if J == 0.0 and not dJ:
            return
        w = -0.5 / area
        term.dP += w * J
        for k, v in dJ.items():
            term.extra[k] = term.extra.get(k, 0.0) + w * v

    def network_flow_inputs(self) -> list[tuple[str, str]]:
        """The neighbour segments' flows, each read at the station it shares."""
        out = []
        if self._prev is not None:
            out.append((self._prev.id, self.from_node))
        if self._next is not None:
            out.append((self._next.id, self.to_node))
        return out

    def _momentum_terms(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None,
    ) -> list[tuple["_StaticChainedDrop", str, str | None, str, NetworkElement | None, str]]:
        """Each momentum term: (C++ result, node it sits at, which derivative is
        this segment's own flow ('a' or 'b' or None), neighbour's node, the
        neighbour element, kind: 'entry', 'from' or 'to')."""
        flows = flows or {}
        m = float(state_in.m_dot)
        incompressible = self.regime == "incompressible"
        to_station = self.to_kappa is not None and self._next is not None
        # Each kernel forms rho from the node's static P and a STATIC T (#782).
        # The node carries T0; its gas moves at the station's mean flow (the
        # one arriving and the one leaving), so both segments sharing a station
        # see the same static state. A duct node with no station downstream
        # moves at this segment's own flow.
        terms = []
        if self.entry_K is not None:
            m_b = flows.get(self._next.id, m) if to_station else m
            st = _moving_static_T(state_out, 0.5 * (m + m_b), self.area_to, incompressible)
            r = cb.channel_entry_drop(m, state_out.P, st.T, state_out.X, self.area_to, self.entry_K)
            nb = self._next if to_station else None
            terms.append(
                (
                    _StaticChainedDrop(r, st, nb is not None),
                    self.to_node,
                    "a",
                    self.to_node,
                    nb,
                    "entry",
                )
            )
        if self.from_kappa is not None:
            m_a = flows.get(self._prev.id, 0.0) if self._prev is not None else 0.0
            st = _moving_static_T(state_in, 0.5 * (m_a + m), self.area_from, incompressible)
            r = cb.station_half_drop(
                m_a, m, state_in.P, st.T, state_in.X, self.area_from, self.from_kappa
            )
            t_from = _StaticChainedDrop(r, st, True)
            self._with_jets(t_from, getattr(self, "_upstream_node", None), state_in, self.area_from)
            terms.append(
                (
                    t_from,
                    self.from_node,
                    "b",
                    self.from_node,
                    self._prev,
                    "from",
                )
            )
        if to_station:
            m_b = flows.get(self._next.id, m)
            st = _moving_static_T(state_out, 0.5 * (m + m_b), self.area_to, incompressible)
            r = cb.station_half_drop(
                m, m_b, state_out.P, st.T, state_out.X, self.area_to, self.to_kappa
            )
            t_to = _StaticChainedDrop(r, st, True)
            self._with_jets(t_to, getattr(self, "_downstream_node", None), state_out, self.area_to)
            terms.append((t_to, self.to_node, "a", self.to_node, self._next, "to"))
        elif (
            incompressible and m < 0.0 and (self.from_kappa is not None or self.entry_K is not None)
        ):
            # Reversed out of an exit node at rest (a plenum, a boundary): the
            # node is a reservoir and the gas accelerates out of it, one
            # dynamic head at the duct's static density, no entry loss. The
            # compressible duct ends take that node's own Pt (``_duct_ends``);
            # without this term the incompressible reversed flow ran ~5x the
            # compressible one.
            node = getattr(self, "_downstream_node", None)
            if node is None or isinstance(node, PlenumNode) or node.is_boundary:
                st = _moving_static_T(state_out, m, self.area_to, True)
                r = cb.channel_entry_drop(-m, state_out.P, st.T, state_out.X, self.area_to, 0.0)
                terms.append(
                    (_ReservoirExitInflow(r), self.to_node, "a", self.to_node, None, "exit_in")
                )
        if incompressible and self.from_kappa is not None and self.area_from != self.area_to:
            # Bernoulli between the two end statics (#812): a station's static
            # at the from-end, the to-end's static (station or exit). A
            # reservoir entry already accelerates to the to-end's section.
            # Compressible, the ends' closures carry it (``_duct_ends``).
            up_id, up = (self.from_node, state_in) if m >= 0.0 else (self.to_node, state_out)
            st = _moving_static_T(up, m, self.area_from, True)
            r = cb.segment_area_change_drop(m, up.P, st.T, up.X, self.area_from, self.area_to)
            terms.append((_StaticChainedDrop(r, st, False), up_id, "a", up_id, None, "area"))
        return terms

    def _exit_head_lost(self, m_dot: float) -> bool:
        """Never on top of the segment's own convention.

        Station nodes hold the channel's STATIC pressure, so an exit with no
        to-station already discharges at static continuity: the dynamic head
        is lost by construction. ChannelElement's exit head at a pressure
        boundary would charge it a second time (#588).
        """
        return False

    def _seed_exit_head(self) -> bool:
        """Always: static continuity at the stations loses the head by
        construction (``_exit_head_lost``), so the start guess charges it."""
        return True

    def _seed_flow_applies(self) -> bool:
        """No: the stations hold static pressures, not the Pt the plain
        channel's kernel reads."""
        return False

    def _momentum_drop_value(self, state_in, state_out, flows) -> float:
        return float(sum(t[0].dP for t in self._momentum_terms(state_in, state_out, flows)))

    def residuals(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        if self.residual_scale_kind == "mdot":
            return self._compressible_station_residuals(state_in, state_out, flows)
        res, jac = super().residuals(state_in, state_out)
        row = jac[0]
        own = f"{self.id}.m_dot"
        for r, node, own_side, nb_node, nb, _ in self._momentum_terms(state_in, state_out, flows):
            res[0] -= r.dP
            d_own = r.d_dm_a if own_side == "a" else r.d_dm_b
            row[own] = row.get(own, 0.0) - d_own
            if nb is not None:
                d_nb = r.d_dm_b if own_side == "a" else r.d_dm_a
                names = nb.unknowns()
                for local, coeff in nb.flow_jac_at_node(nb_node, list(range(len(names)))).items():
                    row[names[local]] = row.get(names[local], 0.0) - d_nb * coeff
            for key, d in _station_node_partials(r, node):
                row[key] = row.get(key, 0.0) - d
            for key, d in r.extra.items():
                row[key] = row.get(key, 0.0) - d
        return res, jac

    def _compressible_station_residuals(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None,
    ) -> tuple[list[float], "LazyJacobian"]:
        """m - m_calc: the Fanno flow between the duct's own end states.

        The compressible channel row is a MASS FLOW (#481); the station terms
        are pressures. Subtracting them from that row, as the incompressible
        form does, mixed kg/s with Pa and solved a meaningless equation: a
        1 um bleed segment held 1 kPa at zero flow (#585). The station terms
        instead move the duct's end STATIC pressures, and the Fanno kernel
        runs between the stagnation pressures of those static states
        (``_duct_ends``, #802). As Ma -> 0 the row is the incompressible
        Pt_in - Pt_out = dP_friction + dP_momentum.
        """
        f_mult = self.surface.f_multiplier if self.surface else 1.0
        terms = self._momentum_terms(state_in, state_out, flows)
        src, dst, d_in, d_out = self._duct_ends(state_in, state_out, terms)
        m_calc, _, _ = self._compressible_flow(src, dst, f_mult, derivatives=False)

        def jacobian() -> dict[int, dict[str, float]]:
            _, derivs, _ = self._compressible_flow(src, dst, f_mult)
            # The kernel's two end Pt are the duct ends', chained below.
            k_in = derivs.pop(f"{self.from_node}.Pt", 0.0)
            k_out = derivs.pop(f"{self.to_node}.Pt", 0.0)
            own = f"{self.id}.m_dot"
            row: dict[str, float] = {own: 1.0}
            for name, d in derivs.items():
                row[name] = row.get(name, 0.0) - d
            for k, partials in ((k_in, d_in), (k_out, d_out)):
                for name, d in partials.items():
                    row[name] = row.get(name, 0.0) - k * d
            return {0: row}

        return [state_in.m_dot - m_calc], LazyJacobian(jacobian)

    def _term_partials(self, term: tuple) -> dict[str, float]:
        """A momentum term's partials, keyed by unknown: this segment's flow,
        the neighbour segment's (through its flow at the shared station), and
        the station node's P, T0 and Y."""
        r, node, own_side, nb_node, nb, _ = term
        out: dict[str, float] = {f"{self.id}.m_dot": r.d_dm_a if own_side == "a" else r.d_dm_b}
        if nb is not None:
            d_nb = r.d_dm_b if own_side == "a" else r.d_dm_a
            names = nb.unknowns()
            for local, coeff in nb.flow_jac_at_node(nb_node, list(range(len(names)))).items():
                out[names[local]] = out.get(names[local], 0.0) + d_nb * coeff
        for key, d in _station_node_partials(r, node):
            out[key] = out.get(key, 0.0) + d
        for key, d in r.extra.items():
            out[key] = out.get(key, 0.0) + d
        return out

    def _duct_ends(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        terms: list[tuple],
    ) -> tuple[NetworkMixtureState, NetworkMixtureState, dict[str, float], dict[str, float]]:
        """The Fanno kernel's two ends: stagnation pressures of the duct's own
        end states, and their partials by unknown (#802).

        A station node holds the STATIC pressure at the middle of its station
        (Pt = P, a plenum). Half the station lies between it and the duct, so
        the duct's end statics are

            inlet   P_from - dP_half(from station)
            outlet  P_to + dP_half(to station)   (the node's P at an exit:
                                                  static continuity)

        and each end's stagnation pressure is that static state's, moving at
        this segment's flux |m| / A: the momentum chamber closure at the gas's
        T0 and composition (``_moving_total_P``). Before #802 the kernel was
        fed the node's P as the duct's P0, a dynamic head low: CR-3630 D.2
        B(5,4,1)'s exit segment ran at 99.7% of that false choke, 81% of the
        true one.

        The gas is the one this segment carries: the from-node's (T0, Y) for
        m >= 0, the to-node's for m < 0. Both ends use it, so at m = 0 both
        stagnation pressures are the statics and the row stays continuous
        through reversal. At the root the kernel's exit static is the outlet
        static exactly (adiabatic; a wall-heated segment's closure omits its
        own heat, second order in the dynamic head).

        An end with no station (a reservoir entry, an exit) is a TERMINAL
        node. Where it feeds the segment it is a reservoir: its own Pt. Where
        the segment discharges into it, a node at rest (a plenum, a boundary:
        Pt = P) takes the closure of its P, static continuity with the head
        lost; a moving node (a momentum chamber) takes its Pt, as a
        ChannelElement into a chamber. At rest the two meet at m = 0 (C1).
        A reservoir entry's dynamic head is now the outlet closure's, so only
        the loss K_in/(1 + K_in) of ``channel_entry_drop`` is added to the
        outlet. A segment with neither stations nor an entry is a
        ChannelElement: Pt to Pt, unchanged.

        An end static below the sonic static P* of the flux is a choked end:
        its Pt is continued below Pt*, never on the supersonic branch (see
        ``closure``).
        """
        m = float(state_in.m_dot)
        kinds = {t[5]: t for t in terms}
        if self.from_kappa is None and self.entry_K is None and "to" not in kinds:
            return (
                state_in,
                state_out,
                {f"{self.from_node}.Pt": 1.0},
                {f"{self.to_node}.Pt": 1.0},
            )
        gas_id, gas = (self.from_node, state_in) if m >= 0.0 else (self.to_node, state_out)
        own = f"{self.id}.m_dot"

        def closure(e: float, de: dict[str, float], A: float) -> tuple[float, dict[str, float]]:
            # A: this end's own section (#812).
            if m == 0.0:
                return e, dict(de)  # at rest: Pt = P (the closure is even in m)

            def chained(c: _MovingTotalP, w: float, d: dict[str, float]) -> None:
                # w * (c's own partials in m, T0 and Y) into d.
                d[own] = d.get(own, 0.0) + w * c.dm
                d[f"{gas_id}.Tt"] = d.get(f"{gas_id}.Tt", 0.0) + w * c.dT0
                for k, v in c.dY.items():
                    d[f"{gas_id}.Y[{k}]"] = d.get(f"{gas_id}.Y[{k}]", 0.0) + w * v

            # Subsonic branch only. At fixed flux, Pt(P) has its minimum at the
            # sonic static P* and rises again below it (the supersonic branch,
            # #637): there a falling end static would RAISE its Pt, reverse the
            # drive and admit spurious roots. Below P* the end is choked and its
            # Pt is continued as Pt* (1 - u^2), u = (P* - e) / P*: C1 at P*
            # (dPt/dP = 0 on both sides) and monotone, so the end keeps a slope
            # in its static. Holding Pt* instead added a plateau of its own to
            # the kernel's choke plateau: a free exit pressure ran off to
            # -3e7 Pa on it. Only an iterate sees this: at a root the kernel
            # saturates at its choked flow for any target below Pt*.
            T0 = float(gas.Tt)
            s = _solver_tools.sonic_static_pressure(T0, abs(m) / A, gas.Y)
            Ps = float(s.P)
            if e >= Ps:
                c = _moving_total_P(e, T0, gas.Y, m, A)
                d = {name: c.dP * v for name, v in de.items()}
                chained(c, 1.0, d)
                return c.Pt, d
            cs = _moving_total_P(Ps, T0, gas.Y, m, A)
            u = (Ps - e) / Ps
            w = 1.0 - u * u
            d = {name: 2.0 * cs.Pt * u / Ps * v for name, v in de.items()}
            # Through P*: u's own dependence plus Pt*'s (whose dP is ~0 there).
            k_ps = -2.0 * cs.Pt * u * e / (Ps * Ps) + w * cs.dP
            sgn = 1.0 if m >= 0.0 else -1.0
            d[own] = d.get(own, 0.0) + k_ps * float(s.dP_dG) * sgn / A
            d[f"{gas_id}.Tt"] = d.get(f"{gas_id}.Tt", 0.0) + k_ps * float(s.dP_dTt)
            for k, v in enumerate(s.dP_dY):
                if v != 0.0:
                    d[f"{gas_id}.Y[{k}]"] = d.get(f"{gas_id}.Y[{k}]", 0.0) + k_ps * float(v)
            chained(cs, w, d)
            return cs.Pt * w, d

        def terminal(
            node_id: str, state: NetworkMixtureState, node: Any, feeds: bool, A: float
        ) -> tuple[float, dict[str, float]]:
            at_rest = node is None or isinstance(node, PlenumNode) or node.is_boundary
            if feeds or not at_rest:
                return float(state.Pt), {f"{node_id}.Pt": 1.0}
            return closure(float(state.P), {f"{node_id}.P": 1.0}, A)

        if self.from_kappa is not None:
            t = kinds["from"]
            de = {f"{self.from_node}.P": 1.0}
            for name, v in self._term_partials(t).items():
                de[name] = de.get(name, 0.0) - v
            Pt_in, d_in = closure(float(state_in.P) - t[0].dP, de, self.area_from)
        else:
            Pt_in, d_in = terminal(
                self.from_node,
                state_in,
                getattr(self, "_upstream_node", None),
                m >= 0.0,
                self.area_from,
            )

        if "to" in kinds:
            t = kinds["to"]
            de = {f"{self.to_node}.P": 1.0}
            for name, v in self._term_partials(t).items():
                de[name] = de.get(name, 0.0) + v
            Pt_out, d_out = closure(float(state_out.P) + t[0].dP, de, self.area_to)
        else:
            Pt_out, d_out = terminal(
                self.to_node,
                state_out,
                getattr(self, "_downstream_node", None),
                m < 0.0,
                self.area_to,
            )
        if "entry" in kinds:
            t = kinds["entry"]
            share = self.entry_K / (1.0 + self.entry_K)
            Pt_out += share * t[0].dP
            for name, v in self._term_partials(t).items():
                d_out[name] = d_out.get(name, 0.0) + share * v

        src = copy.copy(state_in)
        src.Pt = Pt_in
        dst = copy.copy(state_out)
        dst.Pt = Pt_out
        return src, dst, d_in, d_out

    def htc_and_T(self, state: NetworkMixtureState, flows: dict[str, float] | None = None):
        return super().htc_and_T(state)

    def diagnostics(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
    ) -> dict[str, float | str]:
        out = super().diagnostics(state_in, state_out)
        terms = self._momentum_terms(state_in, state_out, flows)
        out["dP_momentum"] = float(sum(t[0].dP for t in terms))
        if self.regime == "compressible":
            # The duct the row solves, between its own end states (#802), not
            # ChannelElement's node-to-node one.
            src, dst, _, _ = self._duct_ends(state_in, state_out, terms)
            f_mult = self.surface.f_multiplier if self.surface else 1.0
            _, _, flow = self._compressible_flow(src, dst, f_mult, derivatives=False)
            out["choked"] = float(flow.choked)
            out["M_in_duct"] = float(flow.M_in)
            out["M_exit_duct"] = float(flow.M_exit)
        return out


class ImpingementCrossflowElement(CrossflowSegmentElement):
    """The crossflow channel between a jet plate and its target, one row pitch (#465).

    A ``CrossflowSegmentElement`` with normal-injection stations
    (``kappa = 0``): the jets that merge at each row arrive with NO
    streamwise momentum, so the crossflow must accelerate them -- the
    momentum equation behind Florschuetz, Truman and Metzger's (1981) flow
    distribution (Eqs. 7-8), P + G_c^2/rho = const. Without it every row sees
    the same pressure difference and the supply stays uniform; with it the
    downstream rows draw more. Area ``height * span``.

    Chain it with plenum crossflow nodes, one per row:

        c1 --[Crossflow]-- c2 --[Crossflow]-- c3 --[Crossflow]-- exit
         |                  |                  |
      [Plate 1]          [Plate 2]          [Plate 3]

    The neighbours are found from the chain itself: the crossflow segment
    arriving at ``from_node`` (none at row 1: a closed end) and the one
    leaving ``to_node`` (none into the exit, where the dynamic head is lost).
    Centred, each node holds its row's mid-merge static pressure: putting the
    whole merge downstream instead over-fed the downstream rows (Gc/Gj at row
    10 of Florschuetz's strongest-crossflow geometry 13% below Eq. 8, against
    a few percent centred).

    BYPASS INLET (#467). ``bypass_inlet=True`` makes this the lead-in from
    the gap's upstream end to row 1, carrying an initial crossflow. Its
    ``from_node`` has no jets, and what arrives there decides the physics:

    * another impingement gap's last segment (arrays in series): the stream
      continues with its momentum, a station with nothing injected. The two
      gaps must have the same cross-section;
    * anything else: the segment accelerates the bypass from ``from_node``'s
      TOTAL pressure, ``entry_K`` the entry loss. A feed channel keeps its Pt
      into a plenum or a momentum chamber alike, so its dynamic head arrives
      with it: 0 for a bypass arriving along a duct, 0.5 for an abrupt
      sharp-edged entry from a large plenum.

    Florschuetz, Metzger and Su (1982, NASA CR-3630) delivered the initial
    crossflow fully developed, with its momentum, after 16 to 95 hydraulic
    diameters: the first case, or a chamber fed by a long enough duct.
    """

    # A 10-row (5,4,2) array cold: 24 residual evaluations, 1.4 s; from the
    # outlet-referenced proxy: 16, 0.5 s. With an initial crossflow (#467)
    # 97-145 cold against 16-17 warm (#706).
    prefers_warm_start = True

    # The parent's ``entry_K`` is resolved from the topology (None unless a
    # bypass inlet fed from a reservoir); the argument is ``bypass_entry_K``
    # (#1099).
    _init_kwarg_attrs = {"entry_K": "bypass_entry_K"}

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        length: float,
        height: float,
        span: float,
        roughness: float = 0.0,
        regime: CompressibilityLiteral = DEFAULT_REGIME,
        friction_model: FrictionModelLiteral = "haaland",
        bypass_inlet: bool = False,
        entry_K: float = 0.0,
    ) -> None:
        if length <= 0.0 or height <= 0.0 or span <= 0.0:
            raise ValueError(
                "ImpingementCrossflowElement: length, height and span must be positive"
            )
        if entry_K < 0.0:
            raise ValueError("ImpingementCrossflowElement: entry_K must be non-negative")
        super().__init__(
            id,
            from_node,
            to_node,
            length=length,
            area=height * span,
            Dh=2.0 * height * span / (height + span),
            roughness=roughness,
            regime=regime,
            friction_model=friction_model,
            from_kappa=None if bypass_inlet else cb.STATION_KAPPA_MERGE_NORMAL,
            to_kappa=cb.STATION_KAPPA_MERGE_NORMAL,
            entry_K=entry_K if bypass_inlet else None,
        )
        self.height = height
        self.span = span
        self.bypass_inlet = bool(bypass_inlet)
        self.bypass_entry_K = float(entry_K)
        self._inlet_error: str | None = None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        prev = [
            e
            for e in graph.get_upstream_elements(self.from_node)
            if isinstance(e, ImpingementCrossflowElement)
        ]
        nxt = [
            e
            for e in graph.get_downstream_elements(self.to_node)
            if isinstance(e, ImpingementCrossflowElement)
        ]
        self.prev_seg = prev[0].id if prev else None
        self.next_seg = nxt[0].id if nxt else None
        self._inlet_error = None
        if self.bypass_inlet:
            # Continuing an upstream gap: a station with nothing injected, so
            # the arriving stream keeps its momentum. Otherwise a reservoir.
            continues = bool(prev)
            self.from_kappa = cb.STATION_KAPPA_MERGE_NORMAL if continues else None
            self.entry_K = None if continues else self.bypass_entry_K
            if continues and abs(prev[0].area / self.area - 1.0) > 1e-9:
                self._inlet_error = (
                    f"it continues {prev[0].id!r}, whose gap cross-section differs "
                    f"({prev[0].area:.4g} vs {self.area:.4g} m^2). Gaps in series "
                    "must match; join different gaps through an area change."
                )
        super().resolve_topology(graph)

    def validate(self) -> None:
        super().validate()
        if self._inlet_error:
            raise ValueError(f"ImpingementCrossflowElement {self.id!r}: {self._inlet_error}")


class AreaChangeElement(NetworkElement):
    """
    Area change between two ducts: a sharp or conical expansion or contraction.

    The closures are Idelchik's (1966; diagrams 4-1, 3-9, 5-2), carried to Mach
    by their control volume, in C++ (``area_change.h``, #731): an expansion is
    the compressible Borda-Carnot step times the table's ratio to it (the low-Re
    factor of 4-1, the cone's phi of 5-2); a contraction is a vena contracta
    re-expanding by the same step, which chokes. ``regime="incompressible"``
    runs their M -> 0 limit.

    The residual is in mass-flow form, m - m_calc(Pt_from, Pt_to, ...), as the
    compressible channel and orifice are (#481): the flow the drive supports
    rises monotonically with it and saturates at the sonic flow, so it exists for
    every state. Each direction is fed from its own upstream node, and the two
    are blended C1 through zero drive (``area_change_element_flow``).

    Equal faces (#1023), F0 == F1 to NEIGHBOUR_AREA_RTOL on the equivalent
    diameter, are no area change: a pass-through, as a bore equal to its pipe
    is no orifice (#788). Zero loss, the row Pt_in - Pt_out with the constant
    Jacobian of a lossless connection. An explicit F0 or F1 that contradicts
    the pipe attached at its face is refused (_check_face_against_pipe).
    """

    #: Incompressible validity (#588): the closure's own Mach -- the small
    #: section's for an expansion, the vena contracta's for a contraction.
    validity_mach_keys = ("mach_loss_ref",)

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        F0: float | None = None,
        F1: float | None = None,
        model_type: Literal["sharp", "conical"] = "sharp",
        length: float | None = None,
        D_h: float | None = 0.0,
        regime: Literal["compressible", "incompressible"] = "compressible",
    ):
        super().__init__(id, from_node, to_node)
        self.F0: float | None = F0
        self.F1: float | None = F1
        self.model_type = model_type
        self.length = length if length is not None else 0.0
        self.D_h = D_h
        self.regime = regime
        # Where a resolved area came from, for diagnostics (issue #262).
        self._F0_source = "user" if F0 is not None else ""
        self._F1_source = "user" if F1 is not None else ""

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    @property
    def _pass_through(self) -> bool:
        """Equal faces, F0 == F1 to NEIGHBOUR_AREA_RTOL on the equivalent
        diameter: no area change (#1023). Follows F0/F1 wherever they were
        set or resolved."""
        f0, f1 = self.F0, self.F1
        if not (f0 and f1 and f0 > 0.0 and f1 > 0.0):
            return False
        return abs(math.sqrt(f0 / f1) - 1.0) <= NEIGHBOUR_AREA_RTOL

    @property
    def residual_scale_kind(self) -> str:
        """RESIDUAL-ROW-SCALING: m_dot - m_calc, in kg/s (#731; it was a Pa
        row); a total-pressure equality, Pa, where the faces are equal and
        there is no area change (#1023). Must change with the residual."""
        return "p" if self._pass_through else "mdot"

    def _seed_resistance(self, rho: float, P_ref: float, m_ref: float) -> float | None:
        """Start guess: a Borda-Carnot head at the smaller area; none for
        equal faces (no area change, #1023)."""
        if not self.F0 or not self.F1:
            return None
        if self._pass_through:
            return 0.0
        a, b = sorted((self.F0, self.F1))
        return (1.0 - a / b) ** 2 / (2.0 * rho * a * a)

    def n_equations(self) -> int:
        return 1

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # An explicit face against the pipe at it (#1023): before inference,
        # so an inherited face (equal by construction) is never compared.
        if self._F0_source == "user" and self.F0:
            _check_face_against_pipe(self, graph, "F0", float(self.F0), [self.from_node])
        if self._F1_source == "user" and self.F1:
            _check_face_against_pipe(self, graph, "F1", float(self.F1), [self.to_node])
        # F0/F1 ARE the physics here: their ratio is the area change the
        # element exists to model. A defaulted value invents a contraction or
        # expansion the network does not contain, so an unresolvable area is
        # an error rather than a guess (issue #262).
        if self.F0 is None:
            self.F0, self._F0_source = _infer_flow_area(graph, [self.from_node], exclude=self)
            if self.F0 is None:
                raise _unresolved_area_error(self, "upstream area F0", [self.from_node])
        if self.F1 is None:
            self.F1, self._F1_source = _infer_flow_area(graph, [self.to_node], exclude=self)
            if self.F1 is None:
                raise _unresolved_area_error(self, "downstream area F1", [self.to_node])

    def validate(self) -> None:
        if (self.F0 or 0.0) <= 0:
            raise ValueError(
                f"AreaChangeElement '{self.id}' has invalid Upstream Area F0={self.F0}. Must be > 0."
            )
        if (self.F1 or 0.0) <= 0:
            raise ValueError(
                f"AreaChangeElement '{self.id}' has invalid Downstream Area F1={self.F1}. Must be > 0."
            )
        # F0 == F1 is no area change: a pass-through (#1023), no longer refused.
        if self.D_h is not None and self.D_h < 0:
            raise ValueError(
                f"AreaChangeElement '{self.id}' has invalid Hydraulic Diameter D_h={self.D_h}. Must be >= 0."
            )
        if self.model_type == "conical" and not self.length > 0.0:
            raise ValueError(
                f"AreaChangeElement '{self.id}' is conical: it needs length > 0 "
                f"(got {self.length})."
            )

    def _flow(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        derivatives: bool,
    ) -> Any:
        return cb.area_change_element_flow(
            float(state_in.Pt),
            float(state_in.Tt),
            state_in.Y,
            float(state_out.Pt),
            float(state_out.Tt),
            state_out.Y,
            float(self.F0),
            float(self.F1),
            self.model_type,
            float(self.length),
            float(self.D_h or 0.0),
            self.regime == "incompressible",
            derivatives,
        )

    def residuals(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> tuple[list[float], LazyJacobian]:
        if self._pass_through:
            # Equal faces (#1023): no area change. Zero loss, the flow passes
            # unchanged; its rate is the network's. Exact and linear, so the
            # Jacobian is the constant +-1 of a lossless connection.
            return [state_in.Pt - state_out.Pt], LazyJacobian(
                lambda: {0: {f"{self.from_node}.Pt": 1.0, f"{self.to_node}.Pt": -1.0}}
            )
        # The loss couples stagnation to stagnation, as ChannelElement does
        # (#175, #541): for a plenum or boundary downstream Pt == P.
        m_calc = float(self._flow(state_in, state_out, derivatives=False).m_dot)

        def jacobian() -> dict[int, dict[str, float]]:
            f = self._flow(state_in, state_out, derivatives=True)
            a, b = self.from_node, self.to_node
            jac: dict[str, float] = {
                f"{self.id}.m_dot": 1.0,
                f"{a}.Pt": -float(f.dm_dPt_a),
                f"{b}.Pt": -float(f.dm_dPt_b),
                f"{a}.Tt": -float(f.dm_dTt_a),
                f"{b}.Tt": -float(f.dm_dTt_b),
            }
            for node, grad in ((a, f.dm_dY_a), (b, f.dm_dY_b)):
                for i, val in enumerate(grad):
                    if val != 0.0:
                        jac[f"{node}.Y[{i}]"] = -float(val)
            return {0: jac}

        return [state_in.m_dot - m_calc], LazyJacobian(jacobian)

    def diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float | str]:
        f0 = float(self.F0)
        f1 = float(self.F1)
        f_small = min(f0, f1)
        m_dot = float(state_in.m_dot)
        # The upwind node feeds the closure, and the face it enters by.
        up, F_in, F_out = (state_in, f0, f1) if m_dot >= 0.0 else (state_out, f1, f0)
        if self._pass_through:
            return self._pass_through_diagnostics(state_in, state_out, f_small)
        try:
            cs_in = cb.complete_state(state_in.T, state_in.P, state_in.X)
            loss = cb.area_change_loss(
                abs(m_dot),
                float(up.Pt),
                float(up.Tt),
                up.Y,
                F_in,
                F_out,
                self.model_type,
                float(self.length),
                float(self.D_h or 0.0),
                self.regime == "incompressible",
            )
        except Exception:
            return {"dPt_loss": 0.0, "area_ratio": f1 / f0}
        a_in = cs_in.thermo.a
        rho_in = cs_in.thermo.rho
        v_in = abs(m_dot) / (rho_in * f0) if rho_in > 0 else 0.0
        mach_in = v_in / a_in if a_in > 0 else 0.0
        try:
            cs_out = cb.complete_state(state_out.T, state_out.P, state_out.X)
            v_out = abs(m_dot) / (cs_out.thermo.rho * f1) if cs_out.thermo.rho > 0 else 0.0
            mach_out = v_out / cs_out.thermo.a if cs_out.thermo.a > 0 else 0.0
        except Exception:
            mach_out = 0.0
        # The section the loss is referenced to, and its Re: the smaller one.
        Dh_small = self.D_h or math.sqrt(4.0 * f_small / math.pi)
        ref = _element_reference_block(
            cs_in, m_dot=m_dot, area=f_small, Dh=Dh_small, location="small_section"
        )
        return {
            "m_dot": m_dot,
            **_element_pressure_block(state_in, state_out, mach_in=mach_in, mach_out=mach_out),
            **ref,
            # The closure's stagnation-pressure loss at the solved flow (>= 0),
            # on the upwind node's state; zeta is it over the small section's
            # static dynamic head, as Idelchik defines it.
            "dPt_loss": float(loss.dPt),
            "area_ratio": f1 / f0,
            "zeta": float(loss.zeta),
            "mach_loss_ref": float(loss.mach),
            "Re_small": float(loss.Re),
            "expansion": float(loss.expansion),
            "choked": float(loss.choked),
            "F0": f0,
            "F1": f1,
            "F0_source": self._F0_source,
            "F1_source": self._F1_source,
        }

    def _pass_through_diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState, area: float
    ) -> dict[str, float | str]:
        """Equal faces (#1023): no closure runs. Zero loss and zeta, no
        choke; the Mach and Re the validity check and the report read are the
        duct's own, at the face area, on the inlet state."""
        m_dot = float(state_in.m_dot)
        f0, f1 = float(self.F0), float(self.F1)
        try:
            cs_in = cb.complete_state(state_in.T, state_in.P, state_in.X)
        except Exception:
            return {"dPt_loss": 0.0, "area_ratio": f1 / f0}
        Dh = self.D_h or math.sqrt(4.0 * area / math.pi)
        ref = _element_reference_block(
            cs_in, m_dot=m_dot, area=area, Dh=Dh, location="small_section"
        )
        a = cs_in.thermo.a
        mach = ref["velocity"] / a if a > 0 else 0.0
        return {
            "m_dot": m_dot,
            **_element_pressure_block(state_in, state_out, mach_in=mach, mach_out=mach),
            **ref,
            "dPt_loss": 0.0,
            "area_ratio": f1 / f0,
            "zeta": 0.0,
            "mach_loss_ref": float(mach),
            "Re_small": float(ref.get("Re", 0.0)),
            "expansion": 0.0,
            "choked": 0.0,
            "F0": f0,
            "F1": f1,
            "F0_source": self._F0_source,
            "F1_source": self._F1_source,
        }


class TeeJunctionElement(NetworkElement):
    """
    Three-port tee junction element using Bassett 2001 pressure-loss coefficients.

    Superseded by ``MultiPortChamberBase`` / ``MultiPortChamberElement`` (momentum-CV
    junction, #177). Retained as the M -> 0 regression baseline for the
    junction validation suite (see ``validation/junction/models/
    tee_junction_element_network.py``) and for the existing Python tests
    under ``python/tests/test_tee_network.py``. Not exposed through the GUI:
    new user-drawn junctions use ``mpce_tee`` instead.

    Merging (tee_type="merging"): straight_node and branch_node are inlets;
    common_node is the outlet.

    Branching (tee_type="branching"): common_node is the inlet; straight_node
    and branch_node are outlets.

    Unknowns: m_dot_com (total common flow) and m_dot_branch (branch flow).
    Straight flow is implicit: m_dot_straight = m_dot_com - m_dot_branch.

    psi = F_C / F_branch (common-arm area / lateral-branch area, default 1.0).
    Bassett 2001 derives K5/K11 (straight-arm coefficients) under the assumption
    F_straight = F_C; the straight-arm area is therefore not a free parameter.
    """

    #: RESIDUAL-ROW-SCALING: total-pressure balances, Pa (both rows).
    residual_scale_kind = "p"

    # to_dict serialises the given areas, not the inherited ones (#1099).
    _init_kwarg_attrs = {
        "F_C": lambda e: e.F_C if e._area_source == "user" else None,
        "F_branch": lambda e: e._F_branch if e._F_branch_user else None,
    }

    def __init__(
        self,
        id: str,
        common_node: str,
        straight_node: str,
        branch_node: str,
        theta: float,
        F_C: float | None = None,
        F_branch: float | None = None,
        psi: float = 1.0,  # fallback ratio when F_branch is None and not inherited
        tee_type: Literal["merging", "branching"] = "merging",
        blend_k: float = 30.0,
    ):
        from_node = straight_node if tee_type == "merging" else common_node
        to_node = common_node if tee_type == "merging" else straight_node
        super().__init__(id, from_node, to_node)
        self.common_node = common_node
        self.straight_node = straight_node
        self.branch_node = branch_node
        self.theta = theta
        self.F_C: float | None = F_C
        # Provenance of F_C: "user" when given, else set during resolve_topology
        # (channel / branch channel / _infer_flow_area source) (issue #559).
        self._area_source = "user" if F_C is not None else ""
        self._F_branch: float | None = (
            F_branch  # resolved branch area; psi derived from F_C/_F_branch
        )
        #: Whether F_branch was given (#1023): only a given one is compared
        #: with the branch pipe.
        self._F_branch_user = F_branch is not None
        self.psi = psi
        self.tee_type = tee_type
        self.blend_k = blend_k

    def all_source_nodes(self) -> list[str]:
        if self.tee_type == "merging":
            return [self.straight_node, self.branch_node]
        return [self.common_node]

    def all_sink_nodes(self) -> list[str]:
        if self.tee_type == "merging":
            return [self.common_node]
        return [self.straight_node, self.branch_node]

    def flow_at_node(self, node_id: str, x: Any, indices: list[int]) -> float:
        m_com = float(x[indices[0]])
        m_branch = float(x[indices[1]])
        if node_id == self.common_node:
            return m_com
        if node_id == self.branch_node:
            return m_branch
        return m_com - m_branch  # straight_node

    def flow_jac_at_node(self, node_id: str, indices: list[int]) -> dict[int, float]:
        if node_id == self.common_node:
            return {indices[0]: 1.0}
        if node_id == self.branch_node:
            return {indices[1]: 1.0}
        return {indices[0]: 1.0, indices[1]: -1.0}  # straight_node

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot_com", f"{self.id}.m_dot_branch"]

    def n_equations(self) -> int:
        return 2

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Step 0 (#1023): an explicit F_C or F_branch against the pipe
        # attached at its arm. F_C is both the common and the straight arm's
        # area (Bassett: F_straight = F_C).
        if self._area_source == "user" and self.F_C:
            _check_face_against_pipe(
                self, graph, "F_C", float(self.F_C), [self.common_node, self.straight_node]
            )
        if self._F_branch_user and self._F_branch:
            _check_face_against_pipe(
                self, graph, "F_branch", float(self._F_branch), [self.branch_node]
            )
        # Step 1: resolve F_C from common or straight arm channels.
        # Branch arm is only used as last resort (back-calculate F_C = F_branch * psi).
        if self.F_C is None:
            found_fc = False
            for node_id in (self.common_node, self.straight_node):
                for e in graph.get_upstream_elements(node_id) + graph.get_downstream_elements(
                    node_id
                ):
                    if e is self:
                        continue
                    if isinstance(e, ChannelElement) and e.diameter is not None:
                        self.F_C = math.pi * (e.diameter / 2.0) ** 2
                        self._area_source = f"channel '{e.id}'"
                        found_fc = True
                        break
                if found_fc:
                    break
            if not found_fc:
                # Last resort: use branch arm channel to back-calculate F_C = F_branch * psi
                for e in graph.get_upstream_elements(
                    self.branch_node
                ) + graph.get_downstream_elements(self.branch_node):
                    if e is self:
                        continue
                    if isinstance(e, ChannelElement) and e.diameter is not None:
                        f_b = math.pi * (e.diameter / 2.0) ** 2
                        if self._F_branch is None:
                            self._F_branch = f_b
                        self.F_C = f_b * self.psi
                        self._area_source = f"branch channel '{e.id}'"
                        found_fc = True
                        break
            if not found_fc:
                # Widen to any area-bearing neighbour before giving up: the
                # common/straight arms may be fed by a chamber or an area
                # change rather than a channel (issue #262).
                area, source = _infer_flow_area(
                    graph, [self.common_node, self.straight_node], exclude=self
                )
                if area is not None:
                    self.F_C = area
                    self._area_source = source
                    found_fc = True
            if not found_fc:
                # F_C sets every loss reference in the junction closure, so a
                # default would silently rescale all of them.
                if getattr(self, "_topo_pass", 0) >= 1:
                    raise _unresolved_area_error(
                        self, "common-arm area F_C", [self.common_node, self.straight_node]
                    )
                self._topo_pass = 1

        # Step 2: resolve F_branch independently from the branch arm channel.
        if self._F_branch is None:
            for e in graph.get_upstream_elements(self.branch_node) + graph.get_downstream_elements(
                self.branch_node
            ):
                if e is self:
                    continue
                if isinstance(e, ChannelElement) and e.diameter is not None:
                    self._F_branch = math.pi * (e.diameter / 2.0) ** 2
                    break

        # Step 3: recompute psi from resolved areas.
        if self.F_C is not None and self._F_branch is not None and self._F_branch > 0:
            self.psi = self.F_C / self._F_branch

        # The port nodes: a MomentumChamberNode port's gas moves, so the
        # kernel gets its static T (#782).
        self._port_nodes = {
            nid: graph.nodes.get(nid)
            for nid in (self.common_node, self.straight_node, self.branch_node)
        }

        # Step 4: propagate duct area to auto-area MCN nodes at each arm.
        # MCN nodes at straight/branch ports cannot see an upstream channel in
        # their own resolve_topology (the upstream is this tee element), so they
        # fall back to the 0.1 m^2 sentinel.  The tee already knows F_C and
        # F_branch, so we inject the correct geometry here.
        if self.F_C is not None and self.F_C > 0:
            D_main = math.sqrt(4.0 * self.F_C / math.pi)
            for nid in (self.common_node, self.straight_node):
                _node = graph.nodes.get(nid)
                if isinstance(_node, MomentumChamberNode) and _node._auto_area and _node.Dh is None:
                    _node.Dh = D_main
                    _node.area = self.F_C
        if self._F_branch is not None and self._F_branch > 0:
            D_branch = math.sqrt(4.0 * self._F_branch / math.pi)
            _bra_node = graph.nodes.get(self.branch_node)
            if (
                isinstance(_bra_node, MomentumChamberNode)
                and _bra_node._auto_area
                and _bra_node.Dh is None
            ):
                _bra_node.Dh = D_branch
                _bra_node.area = self._F_branch

    def validate(self) -> None:
        if self.tee_type not in ("merging", "branching"):
            raise ValueError(
                f"TeeJunctionElement '{self.id}': tee_type must be 'merging' or 'branching'."
            )
        if not (abs(self.theta) <= math.pi / 2.0):
            raise ValueError(
                f"TeeJunctionElement '{self.id}': |theta| must be <= pi/2, got {self.theta}."
            )
        if (self.F_C or 0.0) <= 0.0:
            raise ValueError(f"TeeJunctionElement '{self.id}': F_C must be > 0, got {self.F_C}.")
        if self.psi <= 0.0:
            raise ValueError(f"TeeJunctionElement '{self.id}': psi must be > 0, got {self.psi}.")

    def _make_branch_input(
        self,
        state: NetworkMixtureState,
        m_dot: float,
        A: float,
        theta: float,
        node_id: str,
    ) -> tuple["cb._core.BranchInput", _MovingStaticT]:
        """The kernel's port input and the static T it carries (#782).

        ``BranchInput.T`` is a STATIC temperature: the kernel forms
        rho = P / (R T) and u = m / (rho A) from it beside the static P. The
        node carries T0. At a MomentumChamberNode port the gas moves, and its
        static T is the chamber closure's at this arm's flow and area (the
        port's own state when its area is the arm's). A plenum or boundary
        port holds gas at rest: static = total.
        """
        node = getattr(self, "_port_nodes", {}).get(node_id)
        if isinstance(node, MomentumChamberNode):
            st = _moving_static_T(state, m_dot, A, node.regime == "incompressible")
        else:
            st = _MovingStaticT(float(state.Tt), 1.0, 0.0, 0.0, {})
        bi = cb._core.BranchInput()
        bi.P_static = state.P
        bi.Pt = state.Pt
        bi.T = st.T
        bi.m_dot = m_dot
        bi.A = A
        bi.theta = theta
        bi.gamma_eff = state.gamma()
        bi.R_gas = float(cb.specific_gas_constant(state.X))
        return bi, st

    def residuals(
        self,
        state_com: NetworkMixtureState,
        state_straight: NetworkMixtureState,
        state_branch: NetworkMixtureState,
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        m_com = state_com.m_dot
        m_branch = state_branch.m_dot
        A_branch = self._F_branch if self._F_branch is not None else self.F_C / max(self.psi, 1e-6)

        if self.tee_type == "merging":
            # Suppliers: straight (m_straight > 0) and branch (m_branch > 0)
            # Collector: common (m_com > 0 is the outlet)
            # BranchInput sign: >0 = supplier. Straight and branch supply; common collects.
            ports = self._branch_inputs(
                state_com, state_straight, state_branch, -1.0, m_com, m_branch, A_branch
            )
            res = cb._core.compressible_merging_tee_rj(*(p[0] for p in ports))
            # R_0 = P_str - P_bra  (R_sp: equal static pressure at merge point)
            # R_1 = p0_dat - Pt_com - K_com * q_ref  (R_K,com: loss to collector)
            jac: dict[int, dict[str, float]] = {
                0: {
                    f"{self.straight_node}.P": res.dR0_dP_str,
                    f"{self.branch_node}.P": res.dR0_dP_bra,
                },
                1: {
                    f"{self.common_node}.Pt": res.dR1_dPt_com,
                    f"{self.straight_node}.P": res.dR1_dP_str,
                    f"{self.straight_node}.Tt": res.dR1_dT_str,
                    f"{self.branch_node}.P": res.dR1_dP_bra,
                    f"{self.branch_node}.Tt": res.dR1_dT_bra,
                    f"{self.id}.m_dot_com": res.dR1_dmdot_com,
                    f"{self.id}.m_dot_branch": res.dR1_dmdot_branch,
                },
            }
        else:
            # Supplier: common. Collectors: straight and branch.
            ports = self._branch_inputs(
                state_com, state_straight, state_branch, 1.0, m_com, m_branch, A_branch
            )
            res = cb._core.compressible_branching_tee_rj(*(p[0] for p in ports))
            jac = {
                0: {
                    f"{self.id}.m_dot_com": res.dR0_dmdot_com,
                    f"{self.id}.m_dot_branch": res.dR0_dmdot_branch,
                    f"{self.common_node}.Pt": res.dR0_dPt_com,
                    f"{self.common_node}.P": res.dR0_dP_com,
                    f"{self.common_node}.Tt": res.dR0_dT_com,
                    f"{self.straight_node}.Pt": res.dR0_dPt_str,
                },
                1: {
                    f"{self.id}.m_dot_com": res.dR1_dmdot_com,
                    f"{self.id}.m_dot_branch": res.dR1_dmdot_branch,
                    f"{self.common_node}.Pt": res.dR1_dPt_com,
                    f"{self.common_node}.P": res.dR1_dP_com,
                    f"{self.common_node}.Tt": res.dR1_dT_com,
                    f"{self.branch_node}.Pt": res.dR1_dPt_bra,
                    f"{self.branch_node}.P": res.dR1_dP_bra,
                    f"{self.branch_node}.Tt": res.dR1_dT_bra,
                },
            }

        self._chain_static_T(res, jac, ports)
        return [res.R_0, res.R_1], jac

    def _branch_inputs(
        self,
        state_com: NetworkMixtureState,
        state_straight: NetworkMixtureState,
        state_branch: NetworkMixtureState,
        com_sign: float,
        m_com: float,
        m_branch: float,
        A_branch: float,
    ) -> list[tuple[Any, _MovingStaticT, str, str, dict[str, float]]]:
        """(BranchInput, static T, node, kernel suffix, d(arm m_dot)/d(unknown))
        for the common, straight and branch arms. ``com_sign`` is the
        BranchInput sign of the common arm: +1 supplier (branching), -1
        collector (merging)."""
        u_com, u_bra = f"{self.id}.m_dot_com", f"{self.id}.m_dot_branch"
        arms = (
            (
                state_com,
                com_sign * m_com,
                self.F_C,
                0.0,
                self.common_node,
                "com",
                {u_com: com_sign},
            ),
            (
                state_straight,
                m_com - m_branch,
                self.F_C,
                0.0,
                self.straight_node,
                "str",
                {u_com: 1.0, u_bra: -1.0},
            ),
            (
                state_branch,
                m_branch,
                A_branch,
                abs(self.theta),
                self.branch_node,
                "bra",
                {u_bra: 1.0},
            ),
        )
        out = []
        for state, m, A, theta, node, suffix, dm in arms:
            bi, st = self._make_branch_input(state, m, A, theta, node)
            out.append((bi, st, node, suffix, dm))
        return out

    @staticmethod
    def _chain_static_T(
        res: Any,
        jac: dict[int, dict[str, float]],
        ports: list[tuple[Any, _MovingStaticT, str, str, dict[str, float]]],
    ) -> None:
        """The kernel's d/d(port static T), chained onto what that T depends
        on (#782): the node's T0 (the solver's ``{node}.Tt``), its static P,
        its composition and the arm's flow."""
        for row in (0, 1):
            r = jac[row]
            for _, st, node, suffix, dm in ports:
                dR_dT = float(getattr(res, f"dR{row}_dT_{suffix}"))
                r.pop(f"{node}.Tt", None)
                if dR_dT == 0.0:
                    continue
                r[f"{node}.Tt"] = dR_dT * st.dT0
                if st.dP != 0.0:
                    r[f"{node}.P"] = r.get(f"{node}.P", 0.0) + dR_dT * st.dP
                for name, coeff in dm.items():
                    r[name] = r.get(name, 0.0) + dR_dT * st.dm * coeff
                for k, d in st.dY.items():
                    key = f"{node}.Y[{k}]"
                    r[key] = r.get(key, 0.0) + dR_dT * d

    def diagnostics(
        self,
        state_com: NetworkMixtureState,
        state_straight: NetworkMixtureState,
        state_branch: NetworkMixtureState,
    ) -> dict[str, float | str]:
        m_com = state_com.m_dot
        m_branch = state_branch.m_dot
        if self.tee_type == "merging":
            res = cb._core.merging_tee_residuals_and_jacobian(
                m_dot_com=m_com,
                m_dot_branch=m_branch,
                dP0_straight=state_straight.Pt - state_com.P,
                dP0_branch=state_branch.Pt - state_com.P,
                P_static_com=state_com.P,
                T_com=state_com.Tt,
                Y_com=state_com.Y,
                theta=abs(self.theta),
                psi=self.psi,
                F_C=self.F_C,
                blend_k=self.blend_k,
            )
        else:
            res = cb._core.branching_tee_residuals_and_jacobian(
                m_dot_com=m_com,
                m_dot_branch=m_branch,
                dP0_straight=state_com.Pt - state_straight.P,
                dP0_branch=state_com.Pt - state_branch.P,
                P_static_com=state_com.P,
                T_com=state_com.Tt,
                Y_com=state_com.Y,
                theta=abs(self.theta),
                psi=self.psi,
                F_C=self.F_C,
                blend_k=self.blend_k,
            )
        is_extrapolated = res.status != cb._core.CorrelationValidity.VALID
        return {
            "m_dot_com": float(m_com),
            "m_dot_straight": float(m_com - m_branch),
            "m_dot_branch": float(m_branch),
            "mass_flow_ratio": float(res.q),
            "K_straight": float(res.K_straight),
            "K_branch": float(res.K_branch),
            "correlation_extrapolated": 1.0 if is_extrapolated else 0.0,
            "F_C": float(self.F_C or 0.0),
            "area_source": self._area_source,
        }


# ---------------------------------------------------------------------------
# Momentum-CV junction machinery (momentum cv implementation guide.pdf). The
# residuals live in the subclasses (mpce_element.py), whose closures carry the
# port losses, so no separate turning-loss element belongs on a port (#272).
# ---------------------------------------------------------------------------


class MultiPortChamberBase(NetworkElement):
    """
    Momentum-CV junction element (N >= 2 ports).

    ABSTRACT SINCE 0.6.0. This class owns the topology and port machinery --
    port ordering, the inlet/outlet sign map, area and angle resolution, the
    connecting-element wiring checks -- and no longer carries a residual of
    its own. Its impulse-function model was deprecated in 0.5.0 and removed
    here; ``MultiPortChamberElement`` supersedes it (issue #271).

    Subclasses own one scalar unknown ``{id}.P_jct`` and emit N per-port
    relations plus a global mass residual ``sum_i mdot_port_i = 0``, with a
    per-port orientation derived from the topology (positive = flow out of
    the junction). The per-port relation, and with it the port losses, is
    the subclass's closure -- Mynard's for ``MultiPortChamberElement``, fixed
    K for ``ConstantKTeeElement``. Do not add a turning-loss element on a
    port as well: that double-counts the turning loss (#272).

    Topology contract (see PDF Section 2 + addendum task #10):

    - Each port must be a :class:`MomentumChamberNode` (carries P, Pt at the
      port face). The MCN is marked ``_is_junction_port = True`` during
      :meth:`resolve_topology` so the solver skips its mass-balance row -- the
      junction's R_mass takes over conservation at the port (otherwise the
      MCN's mass row degenerates to 0 = 0 once the junction declares the
      same flow with opposite sign).
    - Exactly one element (channel, loss element, etc.) must connect to each
      port-MCN on the outside. Its m_dot is the port's mass-flow source; the
      junction reads it via the graph at residual-eval time, and reports it
      back through ``flow_at_node`` so the state propagation gives every
      port-MCN its real face flow (collector ports would otherwise see zero
      and their Pt = P + 0.5*rho*v^2 closure would degenerate to Pt = P).
    - Sign convention per port is fixed by the inlet/outlet declaration in
      the constructor. ``inlet_nodes`` carry sign ``-1`` (canonical flow INTO
      junction); ``outlet_nodes`` carry sign ``+1`` (canonical flow OUT). The
      connecting outer element's ``from_node``/``to_node`` is then checked
      against the declared direction at :meth:`resolve_topology` and a clear
      error is raised on mismatch. Runtime flow may still reverse (mdot < 0
      at any port); only the canonical direction is fixed for graph topology.

    Port areas are inherited from the connecting :class:`ChannelElement`'s
    diameter unless explicitly given.

    See ``sources/junction/momentum cv implementation guide.pdf`` Sections 2-3
    and ``sources/junction/junction_model_v3_addendum.md`` Finding 6 for the
    motivation (the K-closure's structural K_run ~ 0.30 ceiling and its
    inability to represent ejector / merge-split direction natively).
    """

    #: |connected element area / port area - 1| above which a port is reported
    #: as an unmodelled area step (``port_area_mismatch``).
    PORT_AREA_MISMATCH_TOL = 0.05

    # Opt-in for the solver's junction split seed (`_junction_split_guess`).
    # True here because a chamber junction's port flows follow from the port
    # pressure differences: through the element's own closure where it offers
    # `closure_split_seed` (#272), else a Bernoulli share. Subclasses whose
    # port flows are set by their own internal physics must turn it off -- see
    # `EjectorElement`.
    seeds_ports_by_pressure_split: bool = True

    # to_dict serialises the constructor's per-side angles and the given
    # port areas, not the merged list or the inherited areas (#1099).
    _init_kwarg_attrs = {
        "inlet_angles_deg": lambda e: e.port_angles_deg[: len(e.inlet_nodes)],
        "outlet_angles_deg": lambda e: e.port_angles_deg[len(e.inlet_nodes) :],
        "port_areas": lambda e: [
            a if src == "user" else None
            for a, src in zip(e.port_areas, e._port_area_sources, strict=True)
        ],
    }

    def __init__(
        self,
        id: str,
        inlet_nodes: list[str],
        outlet_nodes: list[str],
        inlet_angles_deg: list[float] | None = None,
        outlet_angles_deg: list[float] | None = None,
        port_areas: list[float] | None = None,
    ):
        """Build a momentum-CV junction.

        Args:
            id: Element id.
            inlet_nodes: Port-MCN ids where flow CANONICALLY enters the
                junction (junction is graph-downstream of each). Flow can
                reverse at runtime; this just sets the canonical topology
                direction for the graph validator.
            outlet_nodes: Port-MCN ids where flow CANONICALLY leaves the
                junction.
            inlet_angles_deg: Geometric branch angles for the inlet ports
                (length must match inlet_nodes). Default 0 each.
            outlet_angles_deg: Same, for outlet ports. Default 0 each.
            port_areas: Per-port cross-section areas [m^2], ordered as
                ``inlet_nodes + outlet_nodes``. If None, inherited from
                connecting channels at resolve time.
        """
        if not inlet_nodes:
            raise ValueError(f"MultiPortChamberBase '{id}': need >= 1 inlet port.")
        if not outlet_nodes:
            raise ValueError(f"MultiPortChamberBase '{id}': need >= 1 outlet port.")
        n_in = len(inlet_nodes)
        n_out = len(outlet_nodes)

        # NetworkElement base class requires from_node/to_node; use the first
        # inlet and the first outlet as placeholders. The all_source_nodes /
        # all_sink_nodes overrides below carry the real topology.
        super().__init__(id, inlet_nodes[0], outlet_nodes[0])

        self.inlet_nodes = list(inlet_nodes)
        self.outlet_nodes = list(outlet_nodes)
        self.port_nodes = self.inlet_nodes + self.outlet_nodes  # ordered: inlets first
        self.N = len(self.port_nodes)

        # Canonical sign convention (PDF Section 2.2: positive = out of junction):
        #   inlet ports: flow into junction at canonical orientation -> -1
        #   outlet ports: flow out of junction at canonical orientation -> +1
        self._port_signs: list[float] = [-1.0] * n_in + [+1.0] * n_out

        inlet_angles = list(inlet_angles_deg) if inlet_angles_deg is not None else [0.0] * n_in
        outlet_angles = list(outlet_angles_deg) if outlet_angles_deg is not None else [0.0] * n_out
        if len(inlet_angles) != n_in:
            raise ValueError(
                f"MultiPortChamberBase '{id}': inlet_angles_deg has "
                f"{len(inlet_angles)} entries, need {n_in}."
            )
        if len(outlet_angles) != n_out:
            raise ValueError(
                f"MultiPortChamberBase '{id}': outlet_angles_deg has "
                f"{len(outlet_angles)} entries, need {n_out}."
            )
        self.port_angles_deg = inlet_angles + outlet_angles

        self.port_areas: list[float | None] = (
            list(port_areas) if port_areas is not None else [None] * self.N
        )
        if len(self.port_areas) != self.N:
            raise ValueError(
                f"MultiPortChamberBase '{id}': port_areas has "
                f"{len(self.port_areas)} entries, need {self.N}."
            )

        # Resolved during resolve_topology:
        #   self._port_element_ids[i] = id of the element connected at port i
        #   self._outer_areas[i] = that element's own flow area, when it has one
        self._port_element_ids: list[str] = [""] * self.N
        self._outer_areas: list[float | None] = [None] * self.N
        #: Per-port area provenance: "user" (set at construction), "inherited
        #: from <id>" (set in resolve_topology), or "" until resolved (#559).
        self._port_area_sources: list[str] = [
            "user" if a is not None else "" for a in self.port_areas
        ]

    #: The owned unknown's name, after the result-key convention (#773):
    #: ``P_jct`` where it is the junction STATIC pressure (ideal momentum tee,
    #: ejector mixing plane), ``Pt_jct`` where it is the TOTAL (Mynard,
    #: constant K). Before #773 every subclass called it ``P_jct``.
    JUNCTION_UNKNOWN: str = "P_jct"

    @property
    def junction_unknown(self) -> str:
        """The owned unknown's full name, ``{id}.P_jct`` or ``{id}.Pt_jct``."""
        return f"{self.id}.{self.JUNCTION_UNKNOWN}"

    def unknowns(self) -> list[str]:
        return [self.junction_unknown]

    def n_equations(self) -> int:
        return self.N + 1  # N impulse + 1 sum-mass

    def row_scale_kinds(self) -> list[str]:
        """Per-row scale kind ("p" or "mdot") for the solver's row-scaling
        vector (RESIDUAL-ROW-SCALING, see NetworkSolver._build_residual_scales),
        in the same order as `residuals()`'s returned list.

        Default matches this class's own impulse-CV rows: N pressure-
        magnitude rows (P_i + rho_i u_i^2 - P_jct) followed by one mass-
        magnitude row (sum of port mdots). Subclasses whose residual rows
        have a DIFFERENT type pattern (e.g. a row expressed as a mass-flow
        residual rather than a pressure residual) must override this --
        the solver applies the wrong scale silently otherwise, which has
        previously caused severe convergence stalls (see solver.py's
        scaling block for the historical "doomed primary" stall class this
        guards against).
        """
        return ["p"] * self.N + ["mdot"]

    def all_source_nodes(self) -> list[str]:
        # Inlet ports = nodes the junction draws flow FROM (canonical orientation).
        return list(self.inlet_nodes)

    def all_sink_nodes(self) -> list[str]:
        # Outlet ports = nodes the junction delivers flow TO.
        return list(self.outlet_nodes)

    def _port_throughflow_terms(self, node_id: str) -> list[tuple[str, float]]:
        """Outer-element terms whose weighted sum is the junction throughflow
        at ``node_id`` (positive = canonical direction through the port).

        Supplier (inlet) ports and collector ports of a single-supplier
        junction carry exactly their own outer element's flow. A collector
        of a multi-supplier junction (e.g. merge tee outlet) receives the
        sum of the supplier feeds -- consistent with how the solver's state
        propagation composes the collector's total flow from one stream per
        supplier port.
        """
        if node_id not in self.port_nodes:
            return []
        i = self.port_nodes.index(node_id)
        if self._port_signs[i] < 0 or len(self.inlet_nodes) == 1:
            return [(self._port_element_ids[i], 1.0)]
        return [(self._port_element_ids[j], 1.0) for j in range(self.N) if self._port_signs[j] < 0]

    def flow_at_node(self, node_id: str, x: Any, indices: list[int]) -> float:
        # Port-MCN mass rows are skipped by the solver (_is_junction_port);
        # this is consumed by the state propagation instead, so collector-port
        # MCNs see their real face flow in the Pt = P + 0.5*rho*v^2 closure
        # and mix upstream streams with true mass weights. `indices` are this
        # element's own unknowns (P_jct); the port flows live on the outer
        # elements, resolved via the _port_outer_mdot_idx map the solver
        # stashes when it builds the unknown vector.
        idx_map = getattr(self, "_port_outer_mdot_idx", None)
        if not idx_map:
            return 0.0
        total = 0.0
        for outer_id, coeff in self._port_throughflow_terms(node_id):
            global_idx = idx_map.get(outer_id)
            if global_idx is not None:
                total += coeff * float(x[global_idx])
        return total

    def flow_jac_at_node(self, node_id: str, indices: list[int]) -> dict[int, float]:
        idx_map = getattr(self, "_port_outer_mdot_idx", None)
        if not idx_map:
            return {}
        jac: dict[int, float] = {}
        for outer_id, coeff in self._port_throughflow_terms(node_id):
            global_idx = idx_map.get(outer_id)
            if global_idx is not None:
                jac[global_idx] = jac.get(global_idx, 0.0) + coeff
        return jac

    def port_densities(self, states: list[Any]) -> list[tuple[float, dict[str, float]]]:
        """Per port: (static density, d(rho)/d(name)) for the momentum rows.

        From each port chamber's own closure (``MomentumChamberNode.
        static_density``), so the junction and the port node agree on the
        port's static state (#576). A state stand-in without a temperature or
        composition (unit tests) keeps its ``density()`` with d(rho)/dP =
        rho/P; one with them but no resolved port node gets the ideal-gas
        columns at its own T.
        """
        nodes = getattr(self, "_port_node_objs", None) or []
        out: list[tuple[float, dict[str, float]]] = []
        for i, st in enumerate(states):
            port = self.port_nodes[i]
            T = float(getattr(st, "Tt", 0.0) or 0.0)
            X = getattr(st, "X", None)
            P = float(st.P)
            if i < len(nodes) and T > 0.0 and X is not None and P > 0.0:
                out.append(nodes[i].static_density(st))
                continue
            rho = float(st.density())
            d: dict[str, float] = {f"{port}.P": rho / P if P > 0.0 else 0.0}
            if T > 0.0 and X is not None:
                d[f"{port}.Tt"] = -rho / T
                W = float(cb.mwmix(X))
                for k, y in enumerate(st.Y):
                    if y > 0.0:
                        d[f"{port}.Y[{k}]"] = -rho * W / float(cb.species_molar_mass(k))
            out.append((rho, d))
        return out

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # 1. For each port: locate the (unique) connecting element on the
        #    OUTSIDE of the port-MCN and determine the orientation sign.
        self._port_node_objs: list[MomentumChamberNode] = []
        for i, port_id in enumerate(self.port_nodes):
            port_node = graph.nodes.get(port_id)
            if port_node is None:
                raise ValueError(
                    f"MultiPortChamberBase '{self.id}': port node "
                    f"'{port_id}' is not in the network."
                )
            if not isinstance(port_node, MomentumChamberNode):
                raise ValueError(
                    f"MultiPortChamberBase '{self.id}': port '{port_id}' "
                    f"must be a MomentumChamberNode, got "
                    f"'{type(port_node).__name__}'."
                )
            self._port_node_objs.append(port_node)
            # Mark the port-MCN so the solver skips its mass row.
            port_node._is_junction_port = True
            port_node._junction_id = self.id

            outside_elems = [
                e
                for e in (
                    graph.get_upstream_elements(port_id) + graph.get_downstream_elements(port_id)
                )
                if e is not self
            ]
            if len(outside_elems) != 1:
                raise ValueError(
                    f"MultiPortChamberBase '{self.id}': port '{port_id}' "
                    f"must have exactly one outside element (channel, loss "
                    f"element, etc.), found {len(outside_elems)}."
                )
            outer = outside_elems[0]
            self._port_element_ids[i] = outer.id

            # Validate the outer element's orientation matches the declared
            # inlet/outlet for this port. The sign convention was fixed at
            # __init__ from the inlet/outlet split:
            #   inlet port (sign = -1):  outer.to_node should be the port-MCN
            #     (outer feeds flow INTO the port, junction takes it from there).
            #   outlet port (sign = +1): outer.from_node should be the port-MCN
            #     (outer takes flow OUT of port, junction sends flow to there).
            expected_sign = self._port_signs[i]
            if outer.from_node == port_id:
                topo_sign = +1.0  # outer takes flow OUT of port = +1 outflow
            elif outer.to_node == port_id:
                topo_sign = -1.0  # outer feeds flow INTO port = -1 outflow
            else:
                raise ValueError(
                    f"MultiPortChamberBase '{self.id}': connecting element "
                    f"'{outer.id}' at port '{port_id}' has no direct from/to "
                    f"link to that port (multi-port outer elements are not "
                    f"supported)."
                )
            if topo_sign != expected_sign:
                role = "inlet" if expected_sign < 0 else "outlet"
                raise ValueError(
                    f"MultiPortChamberBase '{self.id}': port '{port_id}' is "
                    f"declared as an {role} (sign={expected_sign:+.0f}) but the "
                    f"connecting element '{outer.id}' is wired in the opposite "
                    f"direction (sign={topo_sign:+.0f}). For an inlet port, the "
                    f"outer element should feed INTO the port-MCN (outer.to_node "
                    f"== port). For an outlet port, the outer element should take "
                    f"flow OUT of the port-MCN (outer.from_node == port)."
                )

            # 2. Inherit port area from connecting channel if not set.
            outer_area: float | None = None
            if isinstance(outer, ChannelElement) and outer.diameter is not None:
                outer_area = math.pi * (outer.diameter / 2.0) ** 2
            self._outer_areas[i] = outer_area
            # Only the unset ports inherit; this branch is idempotent across
            # resolve passes (after it runs once port_areas[i] is no longer None),
            # and user-set ports keep the "user" source from construction.
            if self.port_areas[i] is None and outer_area is not None:
                self.port_areas[i] = outer_area
                self._port_area_sources[i] = f"inherited from {outer.id}"
            if self.port_areas[i] is None:
                raise ValueError(
                    f"MultiPortChamberBase '{self.id}': could not infer "
                    f"area at port '{port_id}'. Provide port_areas explicitly "
                    f"or attach a ChannelElement with diameter set."
                )

            # 3. Give auto-sized port-MCNs a real flow area. Nodes resolve
            #    before elements (FlowNetwork.resolve_all_topology), and a
            #    collector port has no upstream channel to inherit Dh from,
            #    so it would otherwise keep MomentumChamberNode's 0.1 m^2
            #    fallback and its Pt = P + 0.5*rho*v^2 closure would see a
            #    near-zero face velocity. Mirrors the MCN auto path (Dh ->
            #    area; the wetted surface is never written, #1022) exactly
            #    so a later node-resolve pass is idempotent and collector
            #    ports match what inlet ports already get from
            #    upstream-channel Dh inheritance.
            if port_node._auto_area and port_node.Dh is None:
                area_i = float(self.port_areas[i])
                port_node.Dh = 2.0 * math.sqrt(area_i / math.pi)
                port_node.area = area_i

    def diagnostics(
        self,
        states: list[NetworkMixtureState],
        P_jct: float,
        port_mdots: list[float] | None = None,
    ) -> dict[str, float]:
        """Per-port bookkeeping: pressures, temperatures, areas, signs, flows.

        MACHINERY, not model. It reports the port map and the states handed
        in and computes no physics, which is why it survives the 0.6.0
        removal of this class's residual. ``MultiPortChamberElement.diagnostics``
        calls it through ``super()`` and adds its closure quantities on top.

        ``P_jct`` here is the owned unknown as handed in. The result-key
        convention (#773) is ``P_jct`` = junction STATIC and ``Pt_jct`` =
        junction TOTAL pressure; a subclass whose unknown is a total pressure
        (Mynard, constant K) overwrites ``P_jct`` with the static one.
        ``port_{i}_connected_area_ratio`` is the connected element's area over
        the port's own area (1 = no hidden step).
        """
        diag: dict[str, float | str] = {"P_jct": float(P_jct), "n_ports": float(self.N)}
        port_area_sources = getattr(self, "_port_area_sources", None) or [""] * self.N
        for i in range(self.N):
            diag[f"port_{i}_P"] = float(states[i].P)
            diag[f"port_{i}_T"] = float(states[i].T)  # static (#782, 0.8)
            diag[f"port_{i}_Tt"] = float(states[i].Tt)
            diag[f"port_{i}_area"] = float(self.port_areas[i] or 0.0)
            diag[f"port_{i}_area_source"] = port_area_sources[i]
            diag[f"port_{i}_sign"] = float(self._port_signs[i])
            if port_mdots is not None and i < len(port_mdots):
                diag[f"port_{i}_m_dot"] = float(port_mdots[i])
        # A port whose area differs from the element connected to it leaves an
        # area step that no element models (#536). Reported, never corrected.
        n_mismatch = 0
        outer_areas = (
            getattr(self, "_outer_areas", None) or [None] * self.N
        )  # unresolved: nothing to compare
        for i in range(self.N):
            outer_area = outer_areas[i]
            area = self.port_areas[i]
            if outer_area is None or not area:
                continue
            ratio = float(outer_area) / float(area)
            diag[f"port_{i}_outer_area"] = float(outer_area)
            diag[f"port_{i}_connected_area_ratio"] = ratio
            if abs(ratio - 1.0) > self.PORT_AREA_MISMATCH_TOL:
                n_mismatch += 1
        diag["port_area_mismatch"] = float(n_mismatch)
        return diag


# ---------------------------------------------------------------------------
# Vortex element: free-vortex swirl redistribution (#661)
# ---------------------------------------------------------------------------


def _swirl_static_pressure(Pt: float, Tt: float, v: float, X: list[float]) -> tuple[float, float]:
    """(static P [Pa], swirl Mach [-]) of a stream at stagnation (Pt, Tt)
    moving at speed v: h = h_t - v^2/2, then isentropic from the stagnation
    state. Variable cp throughout, from C++."""
    if v == 0.0:
        return Pt, 0.0
    h = cb.h_mass(Tt, X) - 0.5 * v * v
    T = cb.calc_T_from_h_mass(h, X, Tt)
    M = abs(v) / cb.speed_of_sound(T, X)
    return float(cb.P_from_stagnation(Pt, Tt, M, X)), float(M)


class VortexElement(NetworkElement):
    """
    Free-vortex swirl redistribution between an inlet radius r_in and an
    outlet radius r_out (#661).

    No work is done on the gas: total pressure and total temperature are
    conserved, and the angular momentum r * V_theta is carried from r_in to
    r_out, so V_theta_out = V_theta_in * r_in / r_out. What changes is the
    static pressure: swirl kinetic energy turns into static pressure going
    outward and back again going inward. That static state is reported in the
    diagnostics; it is not a network unknown.

    The inlet swirl V_theta_in [m/s] is the user's; it defaults to solid-body
    rotation at the shaft speed ``shaft_speed_rpm`` (``omega_rpm`` before
    #773; the old name was removed in 0.8.0): V_theta_in = 2 pi N / 60 * r_in.

    Until #661 this element added the static radial-equilibrium rise of a
    Vatistas n-vortex to the TOTAL pressure at constant Tt: neither a free
    vortex (no Pt change) nor rotor work (Tt rises with Pt), and it destroyed
    entropy for every omega > 0. Its core radius ``r_c`` and shape ``n`` were
    ignored from #661 and removed in 0.8.0; the Vatistas functions themselves
    (``cb.vatistas_*``) are unchanged.
    """

    #: RESIDUAL-ROW-SCALING: the row is Pt_out - Pt_in, in Pa. It was
    #: filed as a mass-flow row by default, weighting it ref_p/ref_mdot
    #: (1e5-1e9) too heavily; LM then traded every pressure for it (#614).
    residual_scale_kind = "p"

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        *,
        r_out: float = 0.10,
        r_in: float = 0.05,
        shaft_speed_rpm: float = 0.0,
        v_theta_in: float | None = None,
    ):
        super().__init__(id, from_node, to_node)
        if r_out <= 0.0:
            raise ValueError(f"Outlet radius r_out must be > 0; got {r_out}")
        if r_in < 0.0:
            raise ValueError(f"Inlet radius r_in must be >= 0; got {r_in}")
        if r_in == 0.0 and v_theta_in:
            raise ValueError(
                "Swirl at r_in = 0 carries no angular momentum: give r_in > 0 "
                f"for v_theta_in = {v_theta_in} m/s"
            )
        self.r_out = float(r_out)
        self.r_in = float(r_in)
        self.shaft_speed_rpm = float(shaft_speed_rpm)
        self.v_theta_in = None if v_theta_in is None else float(v_theta_in)

    def swirl(self) -> tuple[float, float]:
        """(V_theta at r_in, V_theta at r_out) [m/s]; r * V_theta is conserved."""
        if self.v_theta_in is not None:
            v_in = self.v_theta_in
        else:
            v_in = self.shaft_speed_rpm * 2.0 * math.pi / 60.0 * self.r_in
        return v_in, v_in * self.r_in / self.r_out

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def n_equations(self) -> int:
        return 1

    def residuals(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        res = [state_out.Pt - state_in.Pt]
        jac: dict[int, dict[str, float]] = {
            0: {f"{self.to_node}.Pt": 1.0, f"{self.from_node}.Pt": -1.0}
        }
        return res, jac

    def diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float]:
        out: dict[str, float] = {
            "m_dot": float(state_in.m_dot),
            "shaft_speed_rpm": float(self.shaft_speed_rpm),
            **_element_pressure_block(state_in, state_out, mach_in=0.0, mach_out=0.0),
        }
        v_in, v_out = self.swirl()
        out["v_theta_in"] = float(v_in)
        out["v_theta_out"] = float(v_out)
        # The stagnation state is the same at both radii; read it upstream.
        up = state_in if float(state_in.m_dot) >= 0.0 else state_out
        if not (up.Pt > 0.0 and up.Tt > 0.0):
            return out
        try:
            p_in, m_in = _swirl_static_pressure(up.Pt, up.Tt, v_in, up.X)
            p_out, m_out = _swirl_static_pressure(up.Pt, up.Tt, v_out, up.X)
        except (RuntimeError, ValueError):
            # The swirl's kinetic energy exceeds the stream's enthalpy.
            return out
        out.update(
            {
                "P_swirl_in": p_in,
                "P_swirl_out": p_out,
                "dP_swirl_static": p_out - p_in,
                "M_theta_in": m_in,
                "M_theta_out": m_out,
            }
        )
        return out
