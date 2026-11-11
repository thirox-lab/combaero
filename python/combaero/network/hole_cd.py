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

_HoleCdChain holds the resolved options and evaluates the value, its
derivatives by unknown name, the range and envelope flags and each term's
diagnostics. Every hole element holds one (``_hole_cd``) and delegates to it;
the physics is C++'s (crossflow_cd.h, orifice.h), this is the glue. The
element supplies the geometry and topology the chain reads (its bore, hole
count and Reynolds number, its nodes, its duct or crossflow sources).
"""

import math
from dataclasses import dataclass
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
            f"{name}: cd_inclination {cd_inclination!r} is not one of "
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
            f"{name}: cd_inclination={cd_inclination!r} is Idelchik diagram 3-2's inclined "
            "inlet, composed into the thick-walled hole's inlet coefficient (Eq. 4-20); "
            f"it does not compose with correlation={correlation!r}. Use "
            "correlation='IdelchikThick', or cd_inclination='none' (or 'auto')."
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
            f"{name}: cd_exit_crossflow {model!r} is not one of {list(CD_EXIT_CROSSFLOW_MODELS)}."
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
        raise ValueError(f"{name}: cd_exit_crossflow='rowbury_2001': {exc}") from exc
    return _ExitCrossflowSpec(cset, float(theta_deg), float(alpha_deg), float(exit_area_ratio))


def _check_jet_crossflow_law(law: str, xn_d: float, yn_d: float, z_d: float) -> None:
    """Refuse an unknown impingement crossflow law, and CR-3630 Table 4.1's
    on a plate the authors did not test."""
    if law not in _CD_CROSSFLOW_LAWS:
        raise ValueError(
            f"ImpingementPlateElement: cd_crossflow {law!r} is not one of "
            f"{sorted(_CD_CROSSFLOW_LAWS)}"
        )
    if law == "florschuetz_metzger_su_1982" and not (
        cb.florschuetz_metzger_su_1982_has_plate(xn_d, yn_d, z_d)
    ):
        raise ValueError(
            "ImpingementPlateElement: cd_crossflow 'florschuetz_metzger_su_1982' is "
            f"CR-3630 Table 4.1, which covers the plates the authors tested only: "
            f"B(5,4,1-3), B(5,8,1), B(10,4,1), B(10,8,1); not B({xn_d:g},{yn_d:g},"
            f"{z_d:g}). Use 'rogers_hersh_1976' (generic) or 'none'."
        )


class _HoleCdChain:
    """The Cd of one hole element: its resolved options and every term.

    Options (resolved at the element's construction, see the module
    docstring for the order the terms compose in):
      inclination_deg   the hole-to-wall angle the discharge correlation sees
                        (90: no inclination term).
      exit_xf           the exit-crossflow term's resolved spec, None: off.
      supply_duct       True on a hole whose supply crossflow U1/Vi is read
                        from a backside duct (a duct-fed EffusionPlateElement,
                        its ``crossflow_segments`` set); 0 otherwise.
      jet_crossflow     the impingement crossflow law, "none" for every other
                        hole. ``jet_crossflow_suspended`` switches it off for
                        the solver's constant-Cd stage.

    Topology (set by the element's resolve_topology): the exit term's
    discharge node and approaching mainstream.
    """

    def __init__(self, owner: "OrificeElement") -> None:
        self._el = owner
        self.inclination_deg = 90.0
        self.exit_xf: _ExitCrossflowSpec | None = None
        self.supply_duct = False
        self.jet_crossflow = "none"
        self.jet_crossflow_suspended = False
        self.exit_node: Any = None
        self.exit_mainstream: tuple[str, ...] = ()
        self.exit_main_src: str | None = None

    # --- options ------------------------------------------------------------

    def configure_cooling_hole(
        self,
        name: str,
        *,
        correlation: str,
        hole_type: str,
        angle_deg: float,
        cd_inclination: str,
        cd_exit_crossflow: str,
        exit_crossflow_region: str,
        exit_crossflow_orientation_deg: float,
        exit_area_ratio: float,
    ) -> str:
        """Resolve an inclined cooling hole's inclination and exit-crossflow
        options (the film and effusion elements); returns the inclination
        model in use, "auto" resolved."""
        model = _resolve_cd_inclination(name, cd_inclination, correlation, hole_type, angle_deg)
        self.inclination_deg = _cd_inclination_deg(name, model, correlation, angle_deg)
        self.exit_xf = _exit_crossflow_spec(
            name,
            cd_exit_crossflow,
            exit_crossflow_region,
            float(exit_crossflow_orientation_deg),
            angle_deg,
            hole_type == "shaped",
            exit_area_ratio,
        )
        return model

    # --- the composed Cd ----------------------------------------------------

    def cd(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None = None,
    ) -> float:
        """The element's Cd: the hole's (hole_cd), under the impingement
        crossflow law where one acts."""
        base = self.hole_cd(state_in, state_out, flows)
        if self.jet_crossflow == "none":
            return base
        return self.jet_crossflow_cd(state_in, state_out, flows, base)[0]

    def dcd(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None = None,
        continue_re_floor: bool = False,
    ) -> dict[str, float]:
        """d(cd)/d(unknown) by name: the crossflow law's own partials where
        it acts (its Cd replaces the hole's), the hole's otherwise."""
        if self.jet_crossflow != "none" and not self.jet_crossflow_suspended:
            base = self.hole_cd(state_in, state_out, flows)
            _, d = self.jet_crossflow_cd(state_in, state_out, flows, base)
            if d is not None:
                return d
        return self.hole_dcd(state_in, state_out, flows, continue_re_floor)

    def hole_cd(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None = None,
    ) -> float:
        """The hole's own Cd: from correlation or user-supplied fixed value,
        with the opt-in exit-crossflow loss on top (exit_crossflow)."""
        cd0 = self.base_cd(state_in, state_out, flows)
        if self.exit_xf is None or self._el._orifice_geom is None:
            return cd0
        return self.exit_crossflow(cd0, state_in, state_out)[0]

    def hole_dcd(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None = None,
        continue_re_floor: bool = False,
    ) -> dict[str, float]:
        """d(hole_cd)/d(unknown) with the opt-in exit-crossflow term: the
        chain rule dCd/dCd0 dCd0/dx + dCd/dk dk/dx (exit_crossflow). The
        own-flow column through the discharge chamber rides under
        _DCD_DM_VIA_CHAMBER (see there)."""
        d0 = self.base_dcd(state_in, state_out, flows, continue_re_floor)
        if self.exit_xf is None or self._el._orifice_geom is None:
            return d0
        cd0 = self.base_cd(state_in, state_out, flows)
        _, r, dk = self.exit_crossflow(cd0, state_in, state_out)
        if r is None:
            return d0
        out = {name: float(r.dCd_dCd0) * v for name, v in d0.items()}
        if r.dCd_dk != 0.0:
            for name, v in dk.items():
                out[name] = out.get(name, 0.0) + float(r.dCd_dk) * v
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

    def base_cd(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None = None,
    ) -> float:
        """The Cd without the exit-crossflow term: the correlation's (with any
        opt-in inclination term) or the user's fixed value."""
        el = self._el
        if not el.use_correlation or el._orifice_geom is None:
            # Clamp manual Cd strictly <= 1.0 to strictly preserve vena contracta physics (A_eff <= A_geom).
            # A correlation with no geometry (before resolve_topology, or
            # bore == pipe, no orifice) has no Cd: the seed's stands in (#872).
            return max(1e-4, min(el._seed_Cd, 1.0))

        # Estimate Re_D from current m_dot and upstream viscosity.
        D_up = el._orifice_geom.D
        if D_up > 0.0:
            ts = cb.transport_state(state_in.Tt, state_in.P, state_in.X)
            mu = ts.mu if ts.mu > 0.0 else 1.8e-5
            m_dot_abs = abs(state_in.m_dot)
            Re_D = (4.0 * m_dot_abs) / (math.pi * D_up * mu) if m_dot_abs > 1e-12 else 1e4
        else:
            Re_D = 1e5  # typical reference when no upstream channel

        flow_state = cb.OrificeState()
        flow_state.Re_D = Re_D
        flow_state.dP = max(state_in.Pt - state_out.P, 1.0)
        flow_state.rho = cb.density(state_in.Tt, state_in.P, state_in.X)
        flow_state.mu = cb.transport_state(state_in.Tt, state_in.P, state_in.X).mu

        # The ISO 5167 correlations below are functions of the PIPE Reynolds
        # number Re_D, but the discharge-hole family (Idelchik,
        # McGreehan-Schotsch) is a function of the HOLE Reynolds number.
        # Feeding Re_D to those was wrong by 1-8% in a pipe, and by up to 38%
        # for a plenum-fed hole where D_up = 0 freezes Re_D at the 1e5
        # fallback and the correlation stops responding to flow entirely.
        Re_hole = el._hole_reynolds(state_in)

        # Correlation selection logic
        metering = _METERING_SELECTORS.get(el.correlation)
        if metering is not None:
            return float(cb.metering_cd_and_dRe(metering, el._orifice_geom, flow_state.Re_D)[0])
        elif el.correlation in _DISCHARGE_SELECTORS:
            # The discharge-hole family: a hole in a wall, not a plate in a
            # pipe, so DischargeHoleGeometry and no pipe diameter -- the
            # reason the old ThickPlate/RoundedEntry arms were wrong: they
            # multiplied an ISO 5167 metering Cd by a correction factor. The
            # supply-side crossflow U1/Vi (0 for a plenum-fed hole) only moves
            # McGreehan-Schotsch; the others have no crossflow term.
            u1_vi, _ = self.supply_crossflow(state_in, state_out, flows)
            return float(
                cb.discharge_cd(
                    _DISCHARGE_SELECTORS[el.correlation],
                    self.discharge_hole(),
                    cb.DischargeHoleState(Re=Re_hole, U1_over_Vi=u1_vi),
                )
            )
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

    def base_dcd(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None = None,
        continue_re_floor: bool = False,
    ) -> dict[str, float]:
        """d(Cd)/d(unknown) of the Cd WITHOUT the exit-crossflow term, keyed
        by name, analytic.

        Only the discharge-hole family carries it (C++'s
        discharge_cd_and_derivatives): through the hole Reynolds number
        (linear in m_dot, so dRe/dm = Re/m) and, for a duct-fed hole, through
        U1/Vi. 'fixed' has none, and the normed metering correlations keep
        their documented gap.

        BELOW THE CORRELATION'S Re FLOOR the Cd value is held, so its true
        dCd/dRe is 0 -- in the m_dot column and in the viscosity columns
        alike. C++ continues the slope there as a Newton aid; that is only
        wanted by a search ON Cd (``continue_re_floor=True``, the
        self-consistent inner solve). Every sensitivity of a flow built on
        the held value takes the zero: with the continued slope, the
        implicit-function denominator 1 - K dCd/dm of a McGreehan-Schotsch
        panel below Re = 1e4 scaled its whole row by ~0.987 (#760, #801).
        """
        el = self._el
        if el._orifice_geom is None or not el.use_correlation:
            return {}
        selector = _DISCHARGE_SELECTORS.get(el.correlation)
        if selector is None:
            return self._metering_dcd(state_in)
        u1_vi, du = self.supply_crossflow(state_in, state_out, flows)
        Re = el._hole_reynolds(state_in)
        _, dCd_dRe, dCd_dU = cb.discharge_cd_and_derivatives(
            selector, self.discharge_hole(), cb.DischargeHoleState(Re=Re, U1_over_Vi=u1_vi)
        )
        if not continue_re_floor and Re < cb.discharge_cd_re_floor(selector):
            dCd_dRe = 0.0
        out: dict[str, float] = {}
        m = float(state_in.m_dot)
        if abs(m) > 1e-12:
            out[f"{el.id}.m_dot"] = float(dCd_dRe) * Re / abs(m) * (1.0 if m > 0 else -1.0)
            if el._orifice_geom.d > 0.0:
                for name, d in self._re_viscosity_dnames(state_in, Re).items():
                    out[name] = out.get(name, 0.0) + float(dCd_dRe) * d
        for name, d in du.items():
            out[name] = out.get(name, 0.0) + float(dCd_dU) * d
        return out

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

    def _metering_dcd(self, state_in: "NetworkMixtureState") -> dict[str, float]:
        """d(Cd)/d(own m_dot) of the ISO metering family, through the pipe
        Reynolds number Re_D = 4 |m| / (pi D mu), linear in the flow. The
        row's flow column missed it: 1% off on a 2 cm orifice at 1600 K. Held
        fixed where base_cd holds Re_D (no upstream pipe, |m| <= 1e-12).
        """
        el = self._el
        D_up = el._orifice_geom.D
        m = float(state_in.m_dot)
        if not (D_up > 0.0) or abs(m) <= 1e-12:
            return {}
        mu = cb.transport_state(state_in.Tt, state_in.P, state_in.X).mu
        mu_fallback = not mu > 0.0
        mu = 1.8e-5 if mu_fallback else mu
        Re = 4.0 * abs(m) / (math.pi * D_up * mu)
        metering = _METERING_SELECTORS.get(el.correlation)
        if metering is None:
            return {}
        dCd_dRe = float(cb.metering_cd_and_dRe(metering, el._orifice_geom, Re)[1])
        out = {f"{el.id}.m_dot": dCd_dRe * Re / abs(m) * (1.0 if m > 0 else -1.0)}
        if not mu_fallback:
            for name, d in self._re_viscosity_dnames(state_in, Re).items():
                out[name] = out.get(name, 0.0) + dCd_dRe * d
        return out

    def supply_crossflow(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None,
    ) -> tuple[float, dict[str, float]]:
        """McGreehan-Schotsch's supply-side U1/Vi and d(U1/Vi)/d(name).

        0 for a plenum-fed hole. A duct-fed one (``supply_duct``: an
        EffusionPlateElement with crossflow segments) reads the backside
        duct's mean velocity at its supply station over the static-referenced
        ideal jet velocity.
        """
        if not self.supply_duct:
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

    def exit_crossflow(
        self,
        cd0: float,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
    ) -> tuple[float, Any, dict[str, float]]:
        """(Cd, C++'s Rowbury2001Cd or None, dk_perp/d(name)) of the opt-in
        exit-crossflow term (#1074, rowbury_2001_* in crossflow_cd.h).

        k_perp = G_h^2 f / (rho_h rho_m u_m^2), f = u_perp / u_h (Fig. 3):
          G_h      the hole's ACTUAL mass flux, |m| / (A_holes AR), AR the
                   exit/inlet area ratio of a fan-shaped hole (the authors
                   take the exit-plane flux, p. 262), 1 otherwise;
          rho_h    the coolant's exit-plane static density, the isentropic
                   expansion from the supply (Pt, Tt) to the discharge
                   static P (rowbury_2001_exit_density): with it the
                   authors' own Fig. 9 predictions are reproduced (the
                   supply-Tt density left Ma_m 1.2 0.008 high);
          rho_m u_m^2  the APPROACHING mainstream's momentum flux, G_m^2 /
                   rho_m: G_m the chamber's main-inlet flow over its area,
                   rho_m the main stream's density at the chamber's main
                   face (mainstream_momentum_flux).
        No mainstream momentum: Cd0. Eq. 8's own limit there is -C; the
        offset fades beyond the fitted k range, so Cd -> Cd0 continuously as
        u_m -> 0 (include/crossflow_cd.h). Reverse (ingesting) flow takes
        the same Cd at |m|: the source covers injection only.
        """
        el = self._el
        spec = self.exit_xf
        q, dq = self.mainstream_momentum_flux(state_out)
        if spec is None or not q > 0.0:
            return cd0, None, {}
        ex = cb.rowbury_2001_exit_density(state_in.Pt, state_in.Tt, state_out.P, state_in.X)
        rho_h = float(ex.rho)
        m = float(state_in.m_dot)
        a_exit = el.area * spec.exit_area_ratio
        kp = cb.rowbury_2001_k_perp(abs(m) / a_exit, rho_h, q, spec.theta_deg, spec.alpha_deg)
        r = cb.rowbury_2001_exit_crossflow_cd(
            cd0, kp.k, spec.coefficient_set, spec.theta_deg, spec.alpha_deg
        )
        dk: dict[str, float] = {}
        if m != 0.0:
            dk[f"{el.id}.m_dot"] = kp.dk_dG * (1.0 if m > 0.0 else -1.0) / a_exit
        dk_drho = float(kp.dk_drho_h)
        for key, d in (
            (f"{el.from_node}.Pt", ex.drho_dPt),
            (f"{el.from_node}.Tt", ex.drho_dTt),
            (f"{el.to_node}.P", ex.drho_dP),
        ):
            dk[key] = dk.get(key, 0.0) + dk_drho * float(d)
        for k, d in enumerate(ex.drho_dY):
            if d != 0.0:
                key = f"{el.from_node}.Y[{k}]"
                dk[key] = dk.get(key, 0.0) + dk_drho * float(d)
        for name, d in dq.items():
            dk[name] = dk.get(name, 0.0) + float(kp.dk_dq_m) * d
        return float(r.Cd), r, dk

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
                f"{type(el).__name__} {el.id!r}: cd_exit_crossflow="
                f"{getattr(el, 'cd_exit_crossflow', '')!r} reads the mainstream's momentum "
                f"flux from the discharge node, which must be a MomentumChamberNode with an "
                f"area; {el.to_node!r} has no gas velocity."
            )
        if not self.exit_mainstream:
            raise ValueError(
                f"{type(el).__name__} {el.id!r}: cd_exit_crossflow needs an approaching "
                f"mainstream, and {el.to_node!r} is fed only by cooling holes (declare "
                "its main_inlet, or feed it a gas stream)."
            )

    # --- 3. the impingement jet-crossflow law (#800) --------------------------

    def jet_crossflow_cd(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None,
        base_cd: float,
    ) -> tuple[float, dict[str, float] | None]:
        """(Cd, dCd/d(name)) under the crossflow law; None for the partials
        where the law holds the hole's own Cd (the base partials then apply).

        G_ideal = sqrt(2 rho dP) along the flow: supply -> channel injecting,
        channel -> supply for a reversed jet (CR-2681's suction branch), rho
        the upstream static density. G_c = |m_c| / channel_area, m_c every
        crossflow source's flow into ``to_node``.
        """
        if self.jet_crossflow == "none" or self.jet_crossflow_suspended:
            return base_cd, None
        el = self._el
        forward = state_in.Pt >= state_out.P
        up, down = (state_in, state_out) if forward else (state_out, state_in)
        up_id, down_id = (el.from_node, el.to_node) if forward else (el.to_node, el.from_node)
        dP = float(up.Pt) - float(down.P)
        rho, drho_dT, drho_dP = _solver_tools.density_and_jacobians(
            float(up.Tt), float(up.P), list(up.X)
        )
        G_i = math.sqrt(2.0 * rho * dP) if dP > 0.0 else 0.0
        flows = flows or {}
        m_c = float(sum(flows.get(eid, 0.0) for eid in el._crossflow_sources))
        G_c = abs(m_c) / el.channel_area
        if self.jet_crossflow == "rogers_hersh_1976":
            r = cb.rogers_hersh_1976_cd(base_cd, G_i, G_c, not forward)
        else:
            r = cb.florschuetz_metzger_su_1982_cd(el.xn_d, el.yn_d, el.z_d, G_i, G_c)
        if r.capped and self.jet_crossflow == "rogers_hersh_1976":
            return base_cd, None
        out: dict[str, float] = {}
        if r.capped or not G_i > 0.0:
            # The table's constant, or no drive: Cd ~ G_i^1.15 has no slope
            # at G_i = 0 that the flow (~ Cd sqrt(dP)) would feel.
            return float(r.Cd), out

        def add(name: str, v: float) -> None:
            if v != 0.0:
                out[name] = out.get(name, 0.0) + v

        dGi_ddP = rho / G_i
        dGi_drho = dP / G_i
        dcd_dGi = float(r.dCd_dG_ideal)
        add(f"{up_id}.Pt", dcd_dGi * dGi_ddP)
        add(f"{down_id}.P", -dcd_dGi * dGi_ddP)
        add(f"{up_id}.Tt", dcd_dGi * dGi_drho * float(drho_dT))
        add(f"{up_id}.P", dcd_dGi * dGi_drho * float(drho_dP))
        # rho = P MW / (R_u T) through the composition (#964): at fixed (T, P)
        # d(rho)/dY_k = rho (1 - MW / W_k) on the normalised mass fractions
        # the kernels differentiate (the same along any sum-preserving change).
        W = float(cb.mwmix(up.X))
        c_rho = dcd_dGi * dGi_drho * rho
        for k in range(len(up.Y)):
            add(f"{up_id}.Y[{k}]", c_rho * (1.0 - W / float(cb.species_molar_mass(k))))
        if m_c != 0.0:
            dcd_dmc = float(r.dCd_dG_c) * (1.0 if m_c > 0.0 else -1.0) / el.channel_area
            for eid in el._crossflow_sources:
                src = el._source_elements[eid]
                names = src.unknowns()
                for local, coeff in src.flow_jac_at_node(
                    el.to_node, list(range(len(names)))
                ).items():
                    add(names[local], dcd_dmc * coeff)
        return float(r.Cd), out

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

    def inclination_diagnostics(self) -> dict[str, float]:
        """The inclined-inlet increment on zeta' (0.0 when off or at 90 deg)."""
        return {"Cd_inclination_dzeta": float(cb.discharge_inclination_dzeta(self.inclination_deg))}

    def supply_crossflow_diagnostics(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None,
    ) -> dict[str, float]:
        """A duct-fed hole's U1/Vi and the Rohde-scored bounds of
        McGreehan-Schotsch's crossflow term (orifice_discharge_coefficient.md):
        about +8% while U1/Vi < 0.2, degrading badly beyond 0.35, worst for
        sharp edges. Empty for a plenum-fed hole."""
        if not self.supply_duct:
            return {}
        u1_vi, _ = self.supply_crossflow(state_in, state_out, flows)
        return {
            "U1_over_Vi": float(u1_vi),
            "crossflow_cd_beyond_8pct": float(u1_vi > 0.2),
            "crossflow_cd_degraded": float(u1_vi > 0.35),
        }

    def exit_crossflow_diagnostics(
        self, state_in: "NetworkMixtureState", state_out: "NetworkMixtureState", flows: Any
    ) -> dict[str, float]:
        """The opt-in exit-crossflow term's state; empty when off, so a
        default element reports exactly what it did before."""
        el = self._el
        spec = self.exit_xf
        if spec is None or el._orifice_geom is None:
            return {}
        cd0 = self.base_cd(state_in, state_out, flows)
        cd, r, _ = self.exit_crossflow(cd0, state_in, state_out)
        q, _ = self.mainstream_momentum_flux(state_out)
        cset, th, al = spec.coefficient_set, spec.theta_deg, spec.alpha_deg
        k_floor = float(cb.rowbury_2001_k_floor(cset, th, al))
        k = math.nan
        if r is not None:
            ex = cb.rowbury_2001_exit_density(state_in.Pt, state_in.Tt, state_out.P, state_in.X)
            g_h = abs(float(state_in.m_dot)) / (el.area * spec.exit_area_ratio)
            k = float(cb.rowbury_2001_k_perp(g_h, float(ex.rho), q, th, al).k)
        return {
            "Cd_no_exit_crossflow": float(cd0),
            "Cd_exit_crossflow_ratio": float(cd / cd0) if cd0 > 0.0 else math.nan,
            "exit_crossflow_rho_u2": float(q),
            # The actual k_perp, and the one the term was evaluated at (the
            # floor where held).
            "exit_crossflow_k_perp": k,
            "exit_crossflow_k_used": float(r.k_used) if r is not None else math.nan,
            "exit_crossflow_k_floor": k_floor,
            "exit_crossflow_delta": float(r.delta) if r is not None else 0.0,
            "exit_crossflow_held": float(r.held) if r is not None else 0.0,
            # Where k_used sits against the set's k tiers (info): 0 inside
            # what its Table 1 rows provably reach (tier 1, data), 0.5 beyond
            # that but inside the authors' own application of the set (tier
            # 2, Fig. 9 for the suction set; Eq. 8 unfaded), 1 beyond both:
            # the offset fades.
            "exit_crossflow_faded": (
                0.0
                if r is None
                else 1.0
                if r.faded
                else 0.5
                if r.k_used > cb.rowbury_2001_k_data_max(cset)
                else 0.0
            ),
            # Off Table 1's angles the floor can lie above the fitted k
            # range: the whole fitted range is then held (info).
            "exit_crossflow_floor_above_fit": float(k_floor > cb.rowbury_2001_k_fit_max(cset)),
            "exit_crossflow_in_range": float(cb.rowbury_2001_in_range(cset, th, al)),
        }

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
            "Cd_crossflow_extrapolated_velocity_ratio": 0.0,
            "Cd_crossflow_extrapolated_l_d": 0.0,
            "Cd_crossflow_extrapolated_grazing_velocity": 0.0,
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
        out["Cd_crossflow_extrapolated_velocity_ratio"] = float(env.velocity_ratio)
        out["Cd_crossflow_extrapolated_l_d"] = float(env.l_over_d)
        out["Cd_crossflow_extrapolated_grazing_velocity"] = float(env.grazing_velocity)
        return out
