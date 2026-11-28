"""1-D critical-mode supersonic ejector model: Huang's entrainment ratio plus
a Kracik & Dvorak mixing closure for the critical back pressure.

Canonical implementation for `EjectorElement` (in this package) and for the
validation/reference-model work in `validation/ejector/` -- lives here rather
than in that dev-only tree so it ships inside the `combaero` wheel (the same
packaging lesson as `_mynard2010.py`: PyPI installs cannot import
`validation.*`). `validation/ejector/models/huang1999.py` re-exports this
module unchanged; that is where the full paper-by-paper accuracy comparison
and research narrative live (see its docstring and
`validation/ejector/README.md`'s Accuracy section) -- this docstring stays
focused on what a caller of this module needs to know.

References
----------
[Huang 1999] B.J. Huang, J.M. Chang, C.P. Wang, V.A. Petrenko,
    "A 1-D analysis of ejector performance",
    International Journal of Refrigeration 22 (1999) 354-364.
    Local copy: sources/ejector/Huang_1d_analysis_ejector.pdf
    Source of `entrainment_ratio` (Eqs. 1-8) and `choked_mass_flow`
    (Eqs. 1, 7). NOT the source of `critical_back_pressure`'s equations --
    see below; Huang's own P_c* chain (Eqs. 9-18) was tested against all 31
    digitized Table 3/4 rows and found systematically low (mean 25%, max
    35%) against each row's T_c*-derived target, for reasons confirmed NOT
    to be gamma or a transcription error (see validation/ejector/README.md).

[Kracik & Dvorak 2016] J. Kracik, V. Dvorak,
    "Development of an Analytical Method for Predicting Flow in a Supersonic
    Air Ejector", EPJ Web of Conferences 114, 02059 (2016),
    DOI: 10.1051/epjconf/201611402059.
    Local copy: sources/ejector/Development_of_an_Analytical_Method_for.pdf
    Source of the mixing closure `critical_back_pressure` implements (their
    Eqs. 7-13, a lambda/q(lambda)/z(lambda) formalism after Christianovic,
    Kiselev and Abramovich): mass, momentum and energy are solved together
    directly from the y-y state to the fully-mixed, subsonic state, with NO
    separate loss coefficient on momentum and NO explicit shock -- conjugate
    supersonic/subsonic roots of the same conservation laws play the shock's
    role, the same way a Fanno/Rayleigh line has two roots. Tested against
    the same 31 rows: mean 6.2%, max 13.2% -- the strongest of four
    candidates evaluated (see validation/ejector/README.md for the full
    comparison, including why Akbarnejad & Ziabasharhagh 2025's momentum-
    dropping alternative was evaluated and not adopted).

Scope
-----
Critical (double-choked) operation is the core: both the primary flow (choked
at the nozzle throat) and the entrained flow (choked at the hypothetical throat
y-y) are sonic, so the critical entrainment ratio omega_crit is fixed and
independent of back pressure (`entrainment_ratio`, `critical_back_pressure`).

The subcritical droop is now also modeled as closed forms (Phase 1A of the
operating-regime extension, ``validation/ejector/OPERATING_REGIMES_DESIGN.md``):
`dead_head_back_pressure` is the zero-entrainment upper anchor P_b0 (the SAME
Kracik & Dvorak closure at omega = 0 -- no new constant), and
`subcritical_entrainment_ratio` / `blended_entrainment_ratio` droop omega from
omega_crit at P_c* to 0 at P_b0 with a C1 smooth-min through the critical point.

The unchoked-primary regime is modeled by `jet_pump_entrainment_ratio` (Phase
1B): a Keenan constant-area mixing solve that reuses the SAME Kracik & Dvorak
mixing closure with a subsonic (rather than supersonic) primary core, so the
jet-pump and ejector branches meet continuously at the primary-choke threshold.
Validated to reduce to the lossless incompressible jet-pump closed form (Sanger,
NASA TN D-4445) as Mach -> 0.

A second anchor pair is selectable (``ANCHOR_CLOSURES``, the element's
``anchor_closure``): `huang_critical_back_pressure` (Huang's own printed chain,
Eqs. 9-18) with `sudden_expansion_dead_head` (Keenan & Neumann's
zero-secondary-flow mixer). Cross-source it scores 2-3x closer on measured air
rigs; as printed, the chain does not reproduce Huang's or Chen's own numbers
(validation/ejector/README.md).

Still out of scope of THIS closed form: back-flow (P_b beyond the forward-flow
window, omega < 0) -- see the design note.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# Paper-calibrated loss coefficients (Huang 1999, Section 4).
ETA_P = 0.95  # primary-nozzle isentropic efficiency
ETA_S = 0.85  # entrained-flow isentropic efficiency
PHI_P = 0.88  # primary-core area-loss coefficient at section y-y (Eq. 5)


@dataclass(frozen=True)
class EjectorGeometry:
    """Ejector area ratios (the only geometry omega depends on).

    Parameters
    ----------
    area_ratio_nozzle : float
        Nozzle exit / throat area, A_p1 / A_t.
    area_ratio_mix : float
        Constant-area section / nozzle throat, A_3 / A_t.
    """

    area_ratio_nozzle: float  # A_p1 / A_t
    area_ratio_mix: float  # A_3 / A_t

    @classmethod
    def from_count(
        cls,
        *,
        count: int,
        area_throat_single: float,
        area_exit_single: float,
        area_mix: float,
    ) -> EjectorGeometry:
        """Geometry for `count` identical primary nozzles feeding one shared
        mixing chamber (a peripheral multi-nozzle ejector, e.g. Kracik &
        Dvorak 2016, rather than Huang's single centered nozzle).

        Parameters
        ----------
        count : int
            Number of identical, independently-choking primary nozzles.
        area_throat_single, area_exit_single : float
            Throat and exit area of ONE nozzle (not the aggregate).
        area_mix : float
            The shared constant-area mixing-chamber cross section A_3 -- one
            physical duct, not per-nozzle.

        Notes
        -----
        Relies on Kracik & Dvorak's assumption that each nozzle can be
        treated independently: the aggregate throat/exit areas are just
        `count` times the single-nozzle values. area_ratio_nozzle
        (A_p1/A_t) is then invariant to count -- every nozzle contributes
        proportionally to both. area_ratio_mix (A_3/A_t) is NOT: the
        aggregate throat area grows with count while the shared mixing
        chamber does not, so the same physical A_3 gives a smaller ratio
        as count increases. `count=1` reduces exactly to the direct
        single-nozzle constructor.
        """
        area_throat_total = count * area_throat_single
        return cls(
            area_ratio_nozzle=area_exit_single / area_throat_single,
            area_ratio_mix=area_mix / area_throat_total,
        )


@dataclass(frozen=True)
class EjectorResult:
    """Outputs of a critical-mode evaluation."""

    omega: float  # entrainment ratio ms / mp
    mach_nozzle_exit: float  # M_p1
    mach_hypothetical_throat: float  # M_py
    area_ratio_primary_core: float  # A_py / A_t
    area_ratio_entrained: float  # A_sy / A_t
    p_mixing_pa: float  # P_py = P_sy, held constant through the mixing zone


@dataclass(frozen=True)
class CriticalPressureResult:
    """Outputs of the critical-back-pressure chain (Kracik & Dvorak's mixing
    closure, their Eqs. 7-13)."""

    p_c_pa: float  # critical back pressure P_c* = recovery_efficiency * p_mixed_stagnation_pa
    p_mixed_stagnation_pa: (
        float  # p03, the lossless (recovery_efficiency=1) mixed stagnation pressure
    )
    temp_mixed_stagnation_k: float  # T03, mixed-flow stagnation temperature
    lambda_mixed: float  # lambda3, dimensionless mixed-flow velocity (subsonic root)
    mach_mixed: float  # M3, Mach number equivalent of lambda3
    recovery_efficiency: float  # coefficient actually applied (echoed back)


@dataclass(frozen=True)
class DeadHeadResult:
    """Outputs of the zero-entrainment (dead-head) back-pressure evaluation --
    the same Kracik & Dvorak mixing closure as `CriticalPressureResult`, but
    taken at omega = 0 (their gam = 0)."""

    p_b0_pa: float  # dead-head back pressure P_b0 = recovery_efficiency * p03(omega=0)
    p_b0_stagnation_pa: float  # p03 at omega = 0, before recovery_efficiency
    lambda_mixed: float  # lambda3 = subsonic conjugate of lambda1 (z3 = z1)
    mach_mixed: float  # M3, Mach number equivalent of lambda3
    recovery_efficiency: float  # coefficient actually applied (echoed back)


@dataclass(frozen=True)
class JetPumpResult:
    """Outputs of the unchoked-primary subsonic jet-pump evaluation."""

    omega: float  # entrainment ratio ms/mp at the solved operating point
    p_mix_pa: float  # solved common mixing-inlet static pressure p1
    lambda_primary: float  # lambda1: primary core dimensionless velocity (< 1, subsonic)
    lambda_secondary: float  # lambda2: entrained-stream dimensionless velocity
    mach_mixed: float  # M3, Mach number of the mixed flow
    boundary: str  # "" (interior root), "back_pressure" (P_b >= dead-head), or
    #                "primary_choke" (P_b <= choke limit -> blends to the ejector)
    recovery_efficiency: float  # coefficient actually applied (echoed back)


def _area_over_astar(mach: float, gamma: float) -> float:
    """Isentropic area ratio A/A* for a given Mach number."""
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0
    return (1.0 / mach) * ((2.0 / gp1) * (1.0 + 0.5 * gm1 * mach * mach)) ** (gp1 / (2.0 * gm1))


def _mach_from_area_ratio_supersonic(
    area_ratio: float, gamma: float, *, iterations: int = 200
) -> float:
    """Supersonic root of A/A* = area_ratio by bisection.

    A/A* is monotonically increasing for M > 1, so bisection is unconditionally
    convergent. Used for the nozzle-exit Mach M_p1 (Eq. 2).
    """
    if area_ratio < 1.0:
        raise ValueError("area ratio A/A* must be >= 1 for a supersonic root")
    lo, hi = 1.0 + 1e-12, 50.0
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        if _area_over_astar(mid, gamma) < area_ratio:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def choked_mass_flow(
    p0: float,
    t0: float,
    area_throat: float,
    gamma: float,
    r_gas: float,
    eta: float = ETA_P,
) -> float:
    """Choked (sonic-throat) mass flow rate (Huang 1999 Eqs. 1 and 7).

    Same closed form serves the primary nozzle throat (`eta=ETA_P`) and the
    entrained flow's own choking at the hypothetical throat y-y
    (`eta=ETA_S`, `area_throat` -> A_sy); `entrainment_ratio` uses the ratio
    form of this equation internally (Eq. 1 / Eq. 7, where the choked-flow
    constant cancels) rather than calling this function twice, since it
    never needs an absolute area. Unlike `entrainment_ratio`, this DOES need
    R, since the result is an absolute mass flow rather than a ratio.

    Parameters
    ----------
    p0, t0 : float
        Stagnation pressure [Pa] and temperature [K] upstream of the throat.
    area_throat : float
        Absolute throat area [m^2].
    gamma, r_gas : float
        Ratio of specific heats [-] and specific gas constant [J/(kg K)].
    eta : float, default ETA_P
        Isentropic efficiency coefficient.
    """
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0
    choke_const = math.sqrt(gamma / r_gas * (2.0 / gp1) ** (gp1 / gm1))
    return (p0 * area_throat / math.sqrt(t0)) * choke_const * math.sqrt(eta)


def cd_nozzle_mass_flow(
    p0: float,
    t0: float,
    p_static_down: float,
    area_throat: float,
    area_exit: float,
    gamma: float,
    r_gas: float,
    eta: float = ETA_P,
    eps_frac: float = 1.0e-3,
) -> float:
    """Converging-diverging nozzle mass flow spanning choked and unchoked.

    The R0 primary residual for the operating-regime extension needs one C1
    relation for both a choked (supersonic) and an unchoked (subsonic venturi)
    primary. When unchoked, the flow is subsonic throughout and the exit
    (area_exit, static `p_static_down`) sets the mass flow; when the exit
    pressure is low enough that the throat reaches sonic, the flow chokes and
    the mass flow saturates at the sonic-throat value (`choked_mass_flow`,
    area_throat). The two branches are equal at the choke threshold (the throat
    passes the same mass either way), so their minimum is continuous; a C1
    smooth-min (sqrt-rounded, the MultiPortChamberElement idiom) rounds the corner for the
    Newton Jacobian. See OPERATING_REGIMES_DESIGN.md sec 3.1.

    Parameters mirror `choked_mass_flow` plus `p_static_down` (the exit/mixing
    static pressure P_py) and `area_exit` (nozzle exit area A_p1).

    Above the supply (``p_static_down > p0``) the nozzle passes REVERSE flow:
    the same closure from ``p_static_down`` back into ``p0``, shifted by the
    smooth-min's r = 1 value so the flow is continuous there (#619). The former
    zero exit flux was flat in P_py, which let zero flow at any P_py above both
    supplies close the ejector rows -- a spurious root at P_py ~ 1e13 Pa.
    """
    # Shared zero-drive sub-floor law (include/flow_drive_floor.h; #1259, #1269):
    # below kJetDriveFloor (1e-3 Pa) of drive, scale the floor's flow by the odd
    # cubic h(t) = 1.25 t - 0.25 t^3 (exponent kJetFloorExponent = 1/2).
    jet_floor = 1.0e-3

    def _forward_raw(p_up: float, p_down: float) -> float:
        cap = choked_mass_flow(p_up, t0, area_throat, gamma, r_gas, eta)
        r = p_down / p_up
        if r >= 1.0:
            exit_flux = 0.0
        else:
            gm1 = gamma - 1.0
            gp1 = gamma + 1.0
            x = (p_down - p_up) / p_up
            mach = math.sqrt(2.0 / gm1 * math.expm1((-gm1 / gamma) * math.log1p(x)))
            # Below the critical ratio the throat is choked: hold the exit flux at
            # its sonic maximum (C1, flux'(M = 1) = 0). The supersonic-exit flux
            # falls with M and dipped under the cap at deep expansion (#619).
            mach = min(mach, 1.0)
            flux = mach * (1.0 + 0.5 * gm1 * mach * mach) ** (-gp1 / (2.0 * gm1))
            exit_flux = (
                area_exit * p_up / math.sqrt(t0) * math.sqrt(gamma / r_gas) * flux * math.sqrt(eta)
            )
        eps = eps_frac * cap
        return 0.5 * (exit_flux + cap - math.sqrt((exit_flux - cap) ** 2 + eps * eps))

    def _forward(p_up: float, p_down: float) -> float:
        drive = p_up - p_down
        if drive >= jet_floor or not (p_up > jet_floor):
            return _forward_raw(p_up, p_down)
        m_floor = _forward_raw(p_up, p_up - jet_floor)
        t = drive / jet_floor
        h = 1.25 * t - 0.25 * t * t * t
        return m_floor * h

    if p_static_down <= p0:
        return _forward(p0, p_static_down)
    if p0 > jet_floor:
        return -_forward(p_static_down, p0)
    cap0 = choked_mass_flow(p0, t0, area_throat, gamma, r_gas, eta)
    offset = math.sqrt(cap0 * cap0 + (eps_frac * cap0) ** 2) - cap0
    return -_forward(p_static_down, p0) - offset


def cd_nozzle_exit_static(
    m_dot: float,
    p0: float,
    t0: float,
    area_throat: float,
    area_exit: float,
    gamma: float,
    r_gas: float,
    eta: float = ETA_P,
    eps_frac: float = 1.0e-3,
    iterations: int = 200,
) -> float:
    """Invert `cd_nozzle_mass_flow` for the exit static pressure P_py that
    passes `m_dot`.

    The derived mixing pressure of the operating-regime element's unchoked
    branch: given a forced primary flow and supply, the C-D nozzle's exit
    static is the P_py the element mixes at (OPERATING_REGIMES_DESIGN.md sec
    6b). `cd_nozzle_mass_flow` is monotone DECREASING in P_py (higher exit
    pressure -> less expansion -> less flow: cap at P_py -> p0*r_crit, zero at
    P_py -> p0), so bisection on ``cd_nozzle_mass_flow(P_py) - m_dot`` in
    ``(p0*r_crit, p0)`` is unconditionally convergent for any 0 < m_dot < cap.

    Value only (the analytic implicit-function derivatives live in the C++
    ``ejector_cd_nozzle_exit_static_and_jacobian``, as ratios of the C-D
    nozzle's own Jacobian: dP_py/dm_dot = 1/(dmdot/dP_py), etc.).
    """
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0
    r_crit = (2.0 / gp1) ** (gamma / gm1)
    lo = p0 * r_crit  # cd_nozzle == cap here (max flow)
    hi = p0 * (1.0 - 1e-12)  # cd_nozzle -> 0 here (no flow)
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        if (
            cd_nozzle_mass_flow(p0, t0, mid, area_throat, area_exit, gamma, r_gas, eta, eps_frac)
            > m_dot
        ):
            lo = mid  # too much flow -> raise P_py
        else:
            hi = mid
    return 0.5 * (lo + hi)


def entrainment_ratio(
    p_g: float,
    t_g: float,
    p_e: float,
    t_e: float,
    geom: EjectorGeometry,
    gamma: float,
) -> EjectorResult:
    """Critical-mode entrainment ratio omega (Eqs. 1-8).

    Parameters
    ----------
    p_g, t_g : float
        Primary stagnation pressure [Pa] and temperature [K] at the nozzle inlet.
    p_e, t_e : float
        Entrained-flow stagnation pressure [Pa] and temperature [K] at the
        suction port.
    geom : EjectorGeometry
        Area ratios A_p1/A_t and A_3/A_t.
    gamma : float
        Ratio of specific heats of the working fluid (constant, ideal-gas).

    Notes
    -----
    Derivation of the closed form:

      * Primary flow chokes at the throat; the supersonic nozzle-exit Mach M_p1
        follows from the area ratio A_p1/A_t (Eq. 2), and P_p1 from the
        isentropic relation (Eq. 3).
      * The entrained flow chokes at the hypothetical throat y-y (M_sy = 1), so
        its static pressure there is P_sy = P_e / ((g+1)/2)^(g/(g-1)) (Eq. 6).
      * The two streams share a common pressure at y-y (assumption 6):
        P_py = P_sy. The primary core reaches y-y isentropically from the nozzle
        exit, which fixes M_py (Eq. 4).
      * The primary-core area there is A_py/A_t = phi_p * (A/A*)(M_py); this
        follows from Eq. 5 because (A/A*)(M_p1) = A_p1/A_t by construction.
      * The entrained-flow area fills the remainder: A_sy/A_t = A_3/A_t
        - A_py/A_t (Eq. 8).
      * Forming omega = ms/mp from Eqs. 1 and 7, the choked mass-flow constant
        cancels, leaving:
            omega = (P_e/P_g)(A_sy/A_t) sqrt(T_g/T_e) sqrt(eta_s/eta_p).
    """
    if p_g <= 0.0 or p_e <= 0.0 or t_g <= 0.0 or t_e <= 0.0:
        raise ValueError(
            f"entrainment_ratio: pressures and temperatures must be positive, "
            f"got p_g={p_g}, t_g={t_g}, p_e={p_e}, t_e={t_e}. A solver-driven "
            f"Newton iterate has wandered into an unphysical intermediate "
            f"state; the caller's solver is expected to catch this and back "
            f"off (see NetworkSolver's penalty-residual recovery)."
        )
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0
    g_over_gm1 = gamma / gm1

    # Primary flow: nozzle-exit Mach and static pressure (Eqs. 2, 3).
    mach_p1 = _mach_from_area_ratio_supersonic(geom.area_ratio_nozzle, gamma)
    p_p1 = p_g / (1.0 + 0.5 * gm1 * mach_p1 * mach_p1) ** g_over_gm1

    # Entrained flow chokes at y-y (M_sy = 1); pressure matching P_py = P_sy.
    p_sy = p_e / ((gp1 / 2.0) ** g_over_gm1)
    p_py = p_sy

    # Invert Eq. 4 for the primary-core Mach at y-y:
    #   1 + (g-1)/2 M_py^2 = (1 + (g-1)/2 M_p1^2) (P_py/P_p1)^(-(g-1)/g)
    core = (1.0 + 0.5 * gm1 * mach_p1 * mach_p1) * (p_py / p_p1) ** (-gm1 / gamma)
    mach_py = math.sqrt(max((core - 1.0) / (0.5 * gm1), 0.0))

    # Areas at y-y (Eqs. 5, 8).
    area_py = PHI_P * _area_over_astar(mach_py, gamma)  # A_py / A_t
    area_sy = geom.area_ratio_mix - area_py  # A_sy / A_t

    omega = (p_e / p_g) * area_sy * math.sqrt(t_g / t_e) * math.sqrt(ETA_S / ETA_P)
    return EjectorResult(
        omega=omega,
        mach_nozzle_exit=mach_p1,
        mach_hypothetical_throat=mach_py,
        area_ratio_primary_core=area_py,
        area_ratio_entrained=area_sy,
        p_mixing_pa=p_py,
    )


def _q_lambda(lam: float, gamma: float) -> float:
    """Aerodynamic mass-flow function q(lambda) (Kracik & Dvorak Eq. 9).

    q(1) = 1 by construction (sonic reference); q -> 0 as lambda approaches
    the vacuum limit sqrt((g+1)/(g-1)).
    """
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0
    return (1.0 - (gm1 / gp1) * lam * lam) ** (1.0 / gm1) * (gp1 / 2.0) ** (1.0 / gm1) * lam


@dataclass(frozen=True)
class _YYState:
    """The y-y mixing-plane state both back-pressure closures start from.

    lambda1/lambda2 are the primary-core and entrained-stream dimensionless
    velocities at y-y, each referenced to its own stream's critical velocity;
    theta21 = T_e/T_g; omega is the critical-mode entrainment ratio. mach_py
    (primary core Mach at y-y) is a geometry + pressure-matching quantity,
    independent of the realized entrainment -- which is why the SAME state
    serves both P_c* (mixing at gam = omega) and P_b0 (mixing at gam = 0)."""

    lambda1: float
    lambda2: float
    theta21: float
    mach_py: float
    omega: float


def _yy_lambda_state(
    p_g: float,
    t_g: float,
    p_e: float,
    t_e: float,
    geom: EjectorGeometry,
    gamma: float,
    r_gas: float,
) -> _YYState:
    """Build the y-y dimensionless-velocity state (Huang Eqs. 9-10, 13-14).

    Shared by `critical_back_pressure` and `dead_head_back_pressure` so the
    y-y construction is single-sourced; the only thing the two differ in is
    the mixing ratio gam handed to `_mixed_flow_stagnation`."""
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0

    entr = entrainment_ratio(p_g, t_g, p_e, t_e, geom, gamma)
    mach_py = entr.mach_hypothetical_throat

    # Static temperatures and velocities at y-y; entrained flow chokes there
    # (M_sy = 1).
    t_py = t_g / (1.0 + 0.5 * gm1 * mach_py * mach_py)
    t_sy = t_e / (gp1 / 2.0)
    v_py = mach_py * math.sqrt(gamma * r_gas * t_py)
    v_sy = math.sqrt(gamma * r_gas * t_sy)  # M_sy = 1

    # Dimensionless velocities, each referenced to its OWN stream's critical
    # velocity anchored at its own (unchanged, isentropic) stagnation temp.
    c1_star = math.sqrt(2.0 * gamma / gp1 * r_gas * t_g)
    c2_star = math.sqrt(2.0 * gamma / gp1 * r_gas * t_e)
    return _YYState(
        lambda1=v_py / c1_star,
        lambda2=v_sy / c2_star,
        theta21=t_e / t_g,  # ratio of stagnation temperatures (Eq. 8)
        mach_py=mach_py,
        omega=entr.omega,
    )


def _mixed_flow_stagnation(
    p_g: float,
    t_g: float,
    p_e: float,
    state: _YYState,
    gam: float,
    gamma: float,
) -> tuple[float, float, float, float]:
    """Kracik & Dvorak mixed-flow stagnation state (their Eqs. 7-13).

    Returns (p03, t03, lambda3, mach3). `gam` is the mixing ratio: at
    gam = omega this is the critical-mode mix (P_c*); at gam = 0 the entrained
    terms vanish, z3 -> z1, lambda3 -> the subsonic conjugate of lambda1, and
    p03 -> p_g * q(lambda1)/q(lambda3) -- the dead-head pressure P_b0."""
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0
    lambda1 = state.lambda1
    lambda2 = state.lambda2
    theta21 = state.theta21

    z1 = lambda1 + 1.0 / lambda1
    z2 = lambda2 + 1.0 / lambda2
    z3 = (z1 + gam * math.sqrt(theta21) * z2) / math.sqrt((1.0 + gam) * (1.0 + gam * theta21))
    lambda3 = (z3 - math.sqrt(z3 * z3 - 4.0)) / 2.0  # subsonic root (Eq. 12)

    q1 = _q_lambda(lambda1, gamma)
    q2 = _q_lambda(lambda2, gamma)
    q3 = _q_lambda(lambda3, gamma)

    # Mixed flow's stagnation pressure and temperature after mixing (Eqs. 7, 11).
    p03 = (
        p_g
        * math.sqrt((1.0 + gam) * (1.0 + gam * theta21))
        / (1.0 + (p_g / p_e) * gam * math.sqrt(theta21) * q1 / q2)
        * q1
        / q3
    )
    t03 = t_g * (1.0 + gam * theta21) / (1.0 + gam)

    # Mach number equivalent of lambda3 (standard lambda<->Mach conversion).
    mach3_sq = 2.0 * lambda3 * lambda3 / (gp1 - gm1 * lambda3 * lambda3)
    mach3 = math.sqrt(mach3_sq)
    return p03, t03, lambda3, mach3


def critical_back_pressure(
    p_g: float,
    t_g: float,
    p_e: float,
    t_e: float,
    geom: EjectorGeometry,
    gamma: float,
    r_gas: float,
    recovery_efficiency: float = 1.0,
) -> CriticalPressureResult:
    """Critical back pressure P_c* via Kracik & Dvorak's mixing closure.

    Continues from the critical-mode mixing state already found by
    `entrainment_ratio` (Huang's Eqs. 1-8): mass, momentum and energy are
    solved together directly from the y-y state to the fully-mixed, subsonic
    state (Kracik & Dvorak's Eqs. 7-13), with no separate loss coefficient on
    momentum and no explicit shock -- the conjugate subsonic root of the same
    conservation laws plays the shock's role.

    Parameters
    ----------
    p_g, t_g, p_e, t_e, geom, gamma : see `entrainment_ratio`.
    r_gas : float
        Specific gas constant of the working fluid [J/(kg K)]. Unlike omega,
        the mixed-flow velocity and temperature are absolute quantities, so
        the mass-flow constant no longer cancels R out.
    recovery_efficiency : float, default 1.0
        Multiplies the lossless mixed stagnation pressure p03 to give P_c*.
        1.0 (the default) applies no artificial loss -- p03 is already the
        fully-recovered stagnation pressure of the mixed flow, so no
        separate diffuser step is needed. Exposed for a caller to tune
        against their own data; not fitted to any reference dataset here.

    Notes
    -----
      * Static temperatures and velocities at y-y follow exactly as in
        `entrainment_ratio` (Eqs. 9-10, 13-14 of Huang 1999).
      * lambda = c/c* is dimensionless velocity referenced to EACH STREAM'S
        OWN critical/sonic velocity c* = sqrt(2g/(g+1) R T0), anchored to
        that stream's own constant stagnation temperature (T_g or T_e) --
        unlike Mach number, lambda needs no local static temperature.
      * z(lambda) = lambda + 1/lambda (Eq. 13) combines the two streams'
        momentum into the mixed flow's z(lambda3) (Eq. 12); inverting
        z(lambda)=z3 has two roots (lambda3, 1/lambda3) -- the smaller
        (subsonic) root is taken, exactly as a normal shock's downstream
        state is the subsonic conjugate of its upstream state.
      * p03 (Eq. 7) is the mixed flow's stagnation pressure after mixing --
        already P_c* at recovery_efficiency=1, since a stagnation pressure
        is by definition the pressure after isentropic recovery to zero
        velocity.
    """
    # Critical mode: the mixing ratio is the critical-mode entrainment ratio
    # (Kracik & Dvorak's Gamma = ms/mp = omega).
    state = _yy_lambda_state(p_g, t_g, p_e, t_e, geom, gamma, r_gas)
    p03, t03, lambda3, mach3 = _mixed_flow_stagnation(p_g, t_g, p_e, state, state.omega, gamma)
    p_c = recovery_efficiency * p03

    return CriticalPressureResult(
        p_c_pa=p_c,
        p_mixed_stagnation_pa=p03,
        temp_mixed_stagnation_k=t03,
        lambda_mixed=lambda3,
        mach_mixed=mach3,
        recovery_efficiency=recovery_efficiency,
    )


def dead_head_back_pressure(
    p_g: float,
    t_g: float,
    p_e: float,
    t_e: float,
    geom: EjectorGeometry,
    gamma: float,
    r_gas: float,
    recovery_efficiency: float = 1.0,
) -> DeadHeadResult:
    """Zero-entrainment (dead-head) back pressure P_b0.

    P_b0 is the highest back pressure the ejector can hold with NO secondary
    flow -- the upper anchor of the subcritical droop (omega -> 0), the mirror
    of `critical_back_pressure`'s lower anchor (omega = omega_crit). It is NOT
    a new correlation or fitted constant: it is Kracik & Dvorak's SAME mixing
    closure evaluated at omega = 0 (their gam = 0), reusing the identical y-y
    state. At gam = 0 the entrained terms drop out, z3 collapses to z1, so
    lambda3 is the subsonic conjugate of lambda1 and
    ``p03 = p_g * q(lambda1)/q(lambda3)`` -- see OPERATING_REGIMES_DESIGN.md
    sec 2.2. Parameters are exactly those of `critical_back_pressure`.

    Valid (P_b0 > P_c*) only where the primary core is still supersonic at y-y
    (M_py > 1, a healthy motive ejector); the degenerate M_py < 1 corner
    coincides with an unchoked primary and is handled by the jet-pump regime,
    not this closed form (OPERATING_REGIMES_DESIGN.md sec 2.3).
    """
    state = _yy_lambda_state(p_g, t_g, p_e, t_e, geom, gamma, r_gas)
    p03, _t03, lambda3, mach3 = _mixed_flow_stagnation(p_g, t_g, p_e, state, 0.0, gamma)
    p_b0 = recovery_efficiency * p03

    return DeadHeadResult(
        p_b0_pa=p_b0,
        p_b0_stagnation_pa=p03,
        lambda_mixed=lambda3,
        mach_mixed=mach3,
        recovery_efficiency=recovery_efficiency,
    )


#: The two operating-anchor closures an EjectorElement can take (P_c*, P_b0):
#:   "kracik_dvorak": `critical_back_pressure` + `dead_head_back_pressure`
#:       (lossless mixing; 14-63% high on P_c* and 62-274% on P_b0
#:       cross-source, validation/ejector/README.md);
#:   "huang1999": `huang_critical_back_pressure` + `sudden_expansion_dead_head`.
ANCHOR_CLOSURES = ("kracik_dvorak", "huang1999")


def huang_mixing_coefficient(area_ratio_mix: float) -> float:
    """phi_m, Huang 1999 Eq. 19: the mixing-loss coefficient on the momentum
    balance, stepped in A_3/A_t (calibrated by Huang on R141b, A_3/A_t
    6.4-10.9; Chen 2013 Eq. 28 uses the same)."""
    if area_ratio_mix > 8.3:
        return 0.80
    return 0.82 if area_ratio_mix >= 6.9 else 0.84


def _mixed_stagnation_temperature(t_g: float, t_e: float, omega: float) -> float:
    return (t_g + omega * t_e) / (1.0 + omega)


def huang_critical_back_pressure(
    p_g: float,
    t_g: float,
    p_e: float,
    t_e: float,
    geom: EjectorGeometry,
    gamma: float,
    r_gas: float,
    recovery_efficiency: float = 1.0,
) -> CriticalPressureResult:
    """Critical back pressure P_c* by Huang 1999's own chain, as printed.

    From the critical y-y state (`entrainment_ratio`): constant-pressure mixing
    at P_py with the loss coefficient phi_m (Eq. 11, `huang_mixing_coefficient`),
    energy (Eq. 12), a normal shock in the mixed flow (Eqs. 16-17) and an
    isentropic diffuser (Eq. 18).

    Fidelity and accuracy are recorded in validation/ejector/README.md: as
    printed, this chain lands 16-29% below Huang's own computed P_c* and 8-18%
    below Chen 2013's, for reasons not found (fourteen hypotheses tested). On
    rigs that are not its source it is the best P_c* scored (-19 to +17%
    cross-source, against +14 to +63% for `critical_back_pressure`).

    The shock branch is C2 through M_m = 1: the stagnation-pressure loss of a
    normal shock is O((M - 1)^3), so the supersonic and subsonic forms meet
    with equal value, slope and curvature.
    """
    g = gamma
    gm1 = g - 1.0
    entr = entrainment_ratio(p_g, t_g, p_e, t_e, geom, g)
    om, mach_py, p_py = entr.omega, entr.mach_hypothetical_throat, entr.p_mixing_pa
    t_py = t_g / (1.0 + 0.5 * gm1 * mach_py * mach_py)
    t_sy = t_e * 2.0 / (g + 1.0)
    v_py = mach_py * math.sqrt(g * r_gas * t_py)
    v_sy = math.sqrt(g * r_gas * t_sy)  # M_sy = 1
    v_m = huang_mixing_coefficient(geom.area_ratio_mix) * (v_py + om * v_sy) / (1.0 + om)
    cp = g * r_gas / gm1
    t03 = _mixed_stagnation_temperature(t_g, t_e, om)
    t_m = t03 - 0.5 * v_m * v_m / cp
    mach_m_sq = v_m * v_m / (g * r_gas * t_m)
    if mach_m_sq > 1.0:
        p3 = p_py * (1.0 + 2.0 * g / (g + 1.0) * (mach_m_sq - 1.0))
        mach3_sq = (1.0 + 0.5 * gm1 * mach_m_sq) / (g * mach_m_sq - 0.5 * gm1)
    else:
        p3, mach3_sq = p_py, mach_m_sq
    p_c = p3 * (1.0 + 0.5 * gm1 * mach3_sq) ** (g / gm1)
    lambda3 = math.sqrt((g + 1.0) * mach3_sq / (2.0 + gm1 * mach3_sq))
    return CriticalPressureResult(
        p_c_pa=recovery_efficiency * p_c,
        p_mixed_stagnation_pa=p_c,
        temp_mixed_stagnation_k=t03,
        lambda_mixed=lambda3,
        mach_mixed=math.sqrt(mach3_sq),
        recovery_efficiency=recovery_efficiency,
    )


#: C1 floor on the dead-head discriminant: at or below it the mixed flow chokes
#: (M_3 = 1); see `sudden_expansion_dead_head`.
_DEAD_HEAD_DISCRIMINANT_FLOOR = 1e-12


def sudden_expansion_dead_head(
    p_g: float,
    t_g: float,
    p_e: float,
    t_e: float,
    geom: EjectorGeometry,
    gamma: float,
    r_gas: float,
    recovery_efficiency: float = 1.0,
) -> DeadHeadResult:
    """Zero-entrainment back pressure P_b0 as a sudden expansion over A_3.

    Keenan & Neumann's zero-secondary-flow constant-area mixer: the choked
    primary leaves the nozzle (M_p1, P_p1 on A_p1) into the mixing tube, the
    dead-headed suction annulus A_3 - A_p1 sits at the secondary pressure P_e,
    and mass, momentum and energy over A_3 give the subsonic exit, recovered
    isentropically. No fitted constant, and A_3/A_t enters (the Kracik &
    Dvorak dead head, `dead_head_back_pressure`, is independent of it).
    Cross-source it lands 3-39% above the measured or CFD breakdown, and 3-5%
    on the rig added after its form was fixed (validation/ejector/README.md).

    The exit Mach follows in closed form. With the impulse per mass flow
    I/m = sqrt(R T_g / gamma) Phi and k = (gamma - 1)/2,

        Phi = (1 + gamma M^2) / (M sqrt(1 + k M^2)),

    squared, is a quadratic in X = M^2:
    (k Phi^2 - gamma^2) X^2 + (Phi^2 - 2 gamma) X - 1 = 0, whose subsonic root
    is X = 2 / (B + sqrt(D)), B = Phi^2 - 2 gamma, D = B^2 + 4 (k Phi^2 -
    gamma^2). D = 0 is the choked exit (X = 1); below it no subsonic exit
    exists, and D is floored C1 so the exit chokes continuously. ``t_e`` does
    not enter: the secondary is dead-headed.
    """
    del t_e
    g = gamma
    gm1 = g - 1.0
    k = 0.5 * gm1
    a1, a3 = geom.area_ratio_nozzle, geom.area_ratio_mix
    m_p1 = _mach_from_area_ratio_supersonic(a1, g)
    p_p1 = p_g / (1.0 + k * m_p1 * m_p1) ** (g / gm1)
    # Choked primary mass flow per unit throat area (isentropic).
    m = p_g * math.sqrt(g / (r_gas * t_g)) * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * gm1))
    impulse = p_p1 * a1 * (1.0 + g * m_p1 * m_p1) + p_e * (a3 - a1)
    phi = impulse / (m * math.sqrt(r_gas * t_g / g))
    b = phi * phi - 2.0 * g
    d = b * b + 4.0 * (k * phi * phi - g * g)
    eps = _DEAD_HEAD_DISCRIMINANT_FLOOR
    if d < eps:
        d = eps * eps / (2.0 * eps - d)
    x = 2.0 / (b + math.sqrt(d))
    t3 = t_g / (1.0 + k * x)
    p3 = m * math.sqrt(r_gas * t3 / g) / (a3 * math.sqrt(x))
    p_b0 = p3 * (1.0 + k * x) ** (g / gm1)
    return DeadHeadResult(
        p_b0_pa=recovery_efficiency * p_b0,
        p_b0_stagnation_pa=p_b0,
        lambda_mixed=math.sqrt((g + 1.0) * x / (2.0 + gm1 * x)),
        mach_mixed=math.sqrt(x),
        recovery_efficiency=recovery_efficiency,
    )


def subcritical_entrainment_ratio(
    omega_crit: float,
    p_c_pa: float,
    p_b0_pa: float,
    p_b_pa: float,
) -> float:
    """Tier-1 (linear) subcritical entrainment droop omega_sub(P_b).

    Between the two physical anchors the model already produces --
    (P_c*, omega_crit) at the critical point and (P_b0, 0) at dead-head -- the
    realized entrainment droops as the back pressure rises. Tier 1 is the
    straight chord between them:

        omega_sub(P_b) = omega_crit * (P_b0 - P_b) / (P_b0 - P_c*)

    Exact at both anchors, monotone non-increasing in P_b, and introduces no
    new constant. The near-linear subcritical droop is well documented
    empirically (Henry et al. 2007 HEFAT2007 Fig. 5, air; ESDU 86030); a
    Tier-2 mixing-curve inversion replaces this chord only if validation shows
    it is too coarse (OPERATING_REGIMES_DESIGN.md sec 3.2).

    Returns > omega_crit for P_b < P_c* (the plateau side -- the smooth-min in
    `blended_entrainment_ratio` clips it back to omega_crit) and < 0 for
    P_b > P_b0 (backflow, out of scope). Requires P_b0 > P_c*.
    """
    if p_b0_pa <= p_c_pa:
        raise ValueError(
            f"subcritical_entrainment_ratio: need P_b0 > P_c*, got "
            f"P_b0={p_b0_pa}, P_c*={p_c_pa}. A collapsed droop window (P_b0 <= "
            f"P_c*) is the degenerate M_py < 1 / unchoked-primary corner, which "
            f"the jet-pump regime handles, not this closed form."
        )
    return omega_crit * (p_b0_pa - p_b_pa) / (p_b0_pa - p_c_pa)


def blended_entrainment_ratio(
    omega_crit: float,
    p_c_pa: float,
    p_b0_pa: float,
    p_b_pa: float,
    eps_frac: float = 1.0e-3,
) -> float:
    """Realized entrainment omega_eff across the critical/subcritical boundary.

    The realized entrainment is the smaller of the choke-limited plateau value
    omega_crit and the back-pressure-limited droop omega_sub(P_b), joined by a
    C1 smooth-min so a Newton solver sees a continuous derivative through the
    critical point P_c* (the same sqrt-smoothing idiom as MultiPortChamberElement's
    soft barrier):

        omega_eff = 0.5 * (a + b - sqrt((a - b)^2 + eps^2)),  a = omega_crit,
                    b = omega_sub(P_b),  eps = eps_frac * omega_crit

    Below P_c* (b > a) it returns omega_crit (flat plateau); above P_c* (b < a)
    it returns omega_sub, drooping to 0 at P_b0. The rounding is bounded:
    |omega_eff - min(a, b)| <= eps/2, worst exactly at the corner a = b.
    """
    omega_sub = subcritical_entrainment_ratio(omega_crit, p_c_pa, p_b0_pa, p_b_pa)
    eps = eps_frac * omega_crit
    return 0.5 * (omega_crit + omega_sub - math.sqrt((omega_crit - omega_sub) ** 2 + eps * eps))


def _lambda_from_expansion(p_stag: float, p_static: float, gamma: float) -> float:
    """Dimensionless velocity lambda = V/c* for isentropic expansion of a
    stream from stagnation `p_stag` to static `p_static`.

    lambda^2 = ((g+1)/(g-1)) * (1 - (p_static/p_stag)^((g-1)/g)); r_gas- and
    temperature-independent (it cancels between V and the critical velocity
    c* = sqrt(2 g/(g+1) R T0)). lambda < 1 for a subsonic (unchoked) stream."""
    gm1, gp1 = gamma - 1.0, gamma + 1.0
    return math.sqrt(max((gp1 / gm1) * (1.0 - (p_static / p_stag) ** (gm1 / gamma)), 0.0))


def jet_pump_entrainment_ratio(
    p_g: float,
    t_g: float,
    p_e: float,
    t_e: float,
    geom: EjectorGeometry,
    gamma: float,
    p_b: float,
    recovery_efficiency: float = 1.0,
) -> JetPumpResult:
    """Unchoked-primary subsonic jet-pump entrainment ratio omega.

    When the primary is UNCHOKED (its supply pressure is too low to choke the
    nozzle throat) there is no supersonic core and the device is a subsonic jet
    pump, not an ejector. The double-choked construction (`entrainment_ratio`)
    does not apply: its geometric area-filling gives an unbounded, unphysical
    omega (the M_py < 1 degeneracy, see OPERATING_REGIMES_DESIGN.md sec 2.3).

    This closure is the Keenan-Neumann-Lustwerk constant-area mixing solve:
    both streams expand isentropically from their stagnation states to a common
    mixing-inlet static pressure p1 (the primary through A_p1, the entrained
    flow through A_3 - A_p1), so omega = ms/mp = f(p1); the mass+momentum+energy
    mixing to the fully-mixed state is EXACTLY Kracik & Dvorak's closure
    (`_mixed_flow_stagnation`, verified identical to a direct momentum balance),
    reused here with a subsonic lambda1 instead of supersonic -- so the jet-pump
    and ejector branches share one mixing model and meet continuously at the
    primary-choke threshold. The physical p1 is the one whose mixed flow, after
    diffuser recovery, discharges at the back pressure p_b:
    recovery_efficiency * p03(p1) = p_b, solved for p1.

    Parameters
    ----------
    p_g, t_g, p_e, t_e, geom, gamma : see `entrainment_ratio`.
    p_b : float
        Back (discharge) pressure [Pa] the mixed flow must recover to.
    recovery_efficiency : float, default 1.0
        Diffuser recovery multiplier on the mixed stagnation pressure, same
        convention as `critical_back_pressure`. 1.0 = ideal (no fitted loss),
        per the calibration policy. NOTE: audited against Sanger (NASA TN
        D-4445) water-jet-pump data, the lossless form overpredicts the head
        ratio by ~1.4x (a real jet pump delivers ~70% of ideal head), so
        recovery_efficiency ~= 0.7 is the data-anchored value a caller would
        set for quantitative accuracy -- left unfitted here (default 1.0).

    Boundaries (set on `JetPumpResult.boundary`):
      * "back_pressure": p_b is at or above the zero-entrainment discharge, so
        no forward entrainment is possible (omega clamped to 0; backflow is out
        of scope).
      * "primary_choke": p_b is at or below the discharge at which the primary
        chokes; the operating point has left the jet-pump regime for the
        ejector regime, and omega is returned at that choke boundary (the value
        the Phase 2 blend hands over to the critical/subcritical closure).
    """
    if not (p_g > 0.0 and p_e > 0.0 and t_g > 0.0 and t_e > 0.0 and p_b > 0.0):
        raise ValueError(
            f"jet_pump_entrainment_ratio: pressures and temperatures must be "
            f"positive, got p_g={p_g}, t_g={t_g}, p_e={p_e}, t_e={t_e}, p_b={p_b}."
        )
    gm1, gp1 = gamma - 1.0, gamma + 1.0
    k = gm1 / gamma
    b = geom.area_ratio_nozzle / geom.area_ratio_mix  # A_p1 / A_3
    if not (0.0 < b < 1.0):
        raise ValueError(
            f"jet_pump_entrainment_ratio: need 0 < A_p1/A_3 < 1, got {b} from "
            f"area_ratio_nozzle={geom.area_ratio_nozzle}, "
            f"area_ratio_mix={geom.area_ratio_mix}."
        )
    theta21 = t_e / t_g
    area_ratio_ann = (1.0 - b) / b  # A_s / A_p1, A_s = A_3 - A_p1

    def omega_of(p1: float) -> float:
        # ms/mp: ratio of compressible mass fluxes (p1 and R_gas cancel) times
        # the annulus/nozzle area ratio.
        r_g = p1 / p_g
        r_e = p1 / p_e
        flux_p = r_g ** (-k) * math.sqrt(max(1.0 - r_g**k, 0.0)) / math.sqrt(t_g)
        flux_s = r_e ** (-k) * math.sqrt(max(1.0 - r_e**k, 0.0)) / math.sqrt(t_e)
        return (flux_s / flux_p) * area_ratio_ann if flux_p > 0.0 else 0.0

    def evaluate(p1: float) -> tuple[float, float, float, float, float]:
        lam1 = _lambda_from_expansion(p_g, p1, gamma)
        lam2 = _lambda_from_expansion(p_e, p1, gamma)
        omega = omega_of(p1)
        state = _YYState(lambda1=lam1, lambda2=lam2, theta21=theta21, mach_py=0.0, omega=omega)
        p03, _t03, _lam3, mach3 = _mixed_flow_stagnation(p_g, t_g, p_e, state, omega, gamma)
        return omega, p03, lam1, lam2, mach3

    # p1 window: below min(p_g, p_e) the entrained flow can flow forward; above
    # p_g * (2/(g+1))^(g/(g-1)) the PRIMARY (higher stagnation) is unchoked.
    r_crit = (2.0 / gp1) ** (gamma / gm1)
    p_hi = p_e * (1.0 - 1e-9)  # secondary just able to move
    p_floor = p_g * r_crit  # primary choke boundary
    if p_floor >= p_hi:
        raise ValueError(
            "jet_pump_entrainment_ratio: no unchoked-primary window "
            f"(primary chokes at {p_floor:.1f} Pa >= secondary inflow limit {p_hi:.1f} Pa); "
            "the primary supply is strong enough that this is an ejector, not a jet pump."
        )

    def residual(p1: float) -> float:
        return recovery_efficiency * evaluate(p1)[1] - p_b

    # p03 (hence the residual) is monotone in p1, so the two ends bracket the
    # feasible back-pressure range. Outside it, clamp to the named boundary.
    res_hi = residual(p_hi)
    res_lo = residual(p_floor)
    if res_hi <= 0.0:  # p_b >= zero-entrainment discharge: no forward pumping
        omega, _p03, lam1, lam2, mach3 = evaluate(p_hi)
        return JetPumpResult(0.0, p_hi, lam1, lam2, mach3, "back_pressure", recovery_efficiency)
    if res_lo >= 0.0:  # p_b <= choke-limit discharge: hand over to the ejector
        omega, _p03, lam1, lam2, mach3 = evaluate(p_floor)
        return JetPumpResult(
            omega, p_floor, lam1, lam2, mach3, "primary_choke", recovery_efficiency
        )

    # residual increases with p1 (higher mixing pressure -> higher discharge),
    # with residual(p_hi) > 0 > residual(p_floor); root where it crosses zero.
    lo, hi = p_floor, p_hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if residual(mid) > 0.0:
            hi = mid  # root is at lower p1
        else:
            lo = mid  # root is at higher p1
    p1 = 0.5 * (lo + hi)
    omega, _p03, lam1, lam2, mach3 = evaluate(p1)
    return JetPumpResult(omega, p1, lam1, lam2, mach3, "", recovery_efficiency)
