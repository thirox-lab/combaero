import contextlib
import math
import os
import time
import traceback
import warnings
from collections.abc import Callable
from enum import StrEnum
from typing import Any, Literal

import numpy as np
import scipy.sparse as sp
from scipy.optimize import root

import combaero as cb
from combaero import _solver_tools

from .components import (
    FLOW_REVERSAL_M_SCALE,
    RAYLEIGH_EXIT_PREFIX,
    ChannelElement,
    CombustorNode,
    CrossflowSegmentElement,
    EnergyBoundary,
    FuelBoundary,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    MultiPortChamberBase,
    NetworkMixtureState,
    NetworkNode,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
    TurbomachineElement,
    WallNode,
    _is_closed_restriction,
    film_discharge_node,
    film_on_wall_side,
    film_sources,
    target_energy_boundaries,
)
from .graph import FlowNetwork

# Relay key of the gas a film's holes discharge into, before their coolant
# mixes in (#694): FILM_INF_PREFIX + provider id.
FILM_INF_PREFIX = "__film_inf__"


def _station_ends(diag: dict[str, Any], element: Any, station_pt: dict[str, tuple]) -> None:
    """A crossflow segment ending at a station reports the station's true Pt
    (#782 audit): the node holds Pt = P, its gas moves along the chain. dPt
    and the total pressure ratio follow. Only the chain's own segments: a
    plate's or hole's jets discharge into the station's STATIC pressure and
    do not recover the crossflow's head, so their drop stays the drop to P."""
    if getattr(element, "station_report", None) is None:
        return
    ends = (
        ("Pt_in", getattr(element, "from_node", None)),
        ("Pt_out", getattr(element, "to_node", None)),
    )
    changed = False
    for key, nid in ends:
        if nid in station_pt and key in diag:
            diag[key] = float(station_pt[nid][0])
            changed = True
    if changed and "Pt_in" in diag and "Pt_out" in diag:
        if "dPt" in diag:
            diag["dPt"] = float(diag["Pt_in"] - diag["Pt_out"])
        if "pr_total" in diag:
            diag["pr_total"] = float(diag["Pt_out"] / diag["Pt_in"]) if diag["Pt_in"] > 0 else 0.0


class SolverTimeoutError(Exception):
    """Raised when the solver exceeds the specified wall-clock timeout."""


class SolveOutcome(StrEnum):
    """Why a solve ended, in one machine-readable value.

    A ``StrEnum``, so it compares and serialises as a plain string and an
    existing ``== "converged"`` keeps working.

    ``NetworkSolver`` used to report only ``__success__`` and a free-text
    ``__message__``, so every consumer that needed to know WHY a solve failed
    matched on substrings of that message -- and the messages come partly from
    SciPy, which is free to reword them. This does that classification once, at
    the point that knows the answer.
    """

    #: A root was found, and it passed every element's consistency check (or
    #: no element had one).
    CONVERGED = "converged"
    #: A root was found and then rejected by an element's own physical
    #: consistency check. ``__inconsistent_elements__`` names them.
    INCONSISTENT = "inconsistent"
    #: The minimiser stopped improving before reaching tolerance.
    NO_PROGRESS = "no_progress"
    #: The method reported success at a non-zero local minimum of ||F||.
    RESIDUAL_TOO_LARGE = "residual_too_large"
    #: The wall-clock budget ran out.
    TIMEOUT = "timeout"
    #: An exception was raised while evaluating residuals.
    ERROR = "error"
    #: Did not converge, with no more specific reason available.
    NOT_CONVERGED = "not_converged"
    #: Nothing to solve: the network is empty or fully constrained.
    NO_UNKNOWNS = "no_unknowns"


#: Environment switch for the starting-point safety net: any value but "" or
#: "0" makes ``NetworkSolver.raise_internal_errors`` default to True, so a
#: numerical exception at the starting point raises instead of returning an
#: ERROR result. The test suite sets it (python/tests/conftest.py) so the net
#: can never hide a bug there; a stress harness can set it to stop at the first
#: finding, or leave it off and read ``__caught_exceptions__``.
SOLVER_STRICT_ENV = "COMBAERO_SOLVER_STRICT"

#: Exceptions a physical-but-extreme state can raise from the property and
#: closure layer (C++ std::runtime_error -> RuntimeError, std::invalid_argument
#: -> ValueError). TypeError, AttributeError, KeyError, IndexError are code
#: bugs, never states: they always raise.
_STATE_EXCEPTIONS = (RuntimeError, ValueError, ArithmeticError)

# A FuelBoundary's combustor with less oxidizer than this [kg/s] has nothing
# to burn with (#892). With every flow vanishing, the Tt row is met by an
# undetermined state, so such a "root" is refused, never returned.
_FUEL_MIN_OXIDIZER_FLOW = 1e-9

# A fuel controller's target proven unreachable before the solve (#1114) must
# lie at least this far [K] outside the window its fixed oxidizer allows; one
# at the window's ends is left to the solve. The smoothing at zero fuel moves
# the combustor's Tt by under 1e-5 K (stress fuels, 300-900 K).
_FUEL_PROOF_MARGIN_K = 0.1

# Two PressureBoundary total pressures closer than this, relative to the
# larger, drive no flow (#1004, the fraction-loss check's "driven" rule). The
# values are the user's inputs, so equal ones are equal floats; the tolerance
# only absorbs a unit round trip in the GUI (bar or psi to Pa: ~1e-16
# relative) and is far below any drive a network means: 1e-9 of 1 bar is
# 1e-4 Pa.
_FRACTION_DRIVE_PT_RTOL = 1e-9

# A node's Q_withheld below this fraction of the heat given is the rounding
# of m * (1/m), not heat withheld: reported as 0.0. The real fade below the
# 1 mg/s uptake floor passes it within 4e-7 (relative) below the floor
# (6e-12 at m = 0.999999 mg/s), so only heat at that level is dropped.
_WITHHELD_RTOL = 1e-12

#: Default wall-clock budget [s] (``NetworkSolver.homotopy_budget``) of the
#: Newton-homotopy last retry (#642), within the
#: caller's timeout. Inside the solver its stress rescues took 1.2-4.2 s; every
#: solve that does not converge pays it -- unbounded, a junction scorecard with
#: expected failures went from 8 s to 106 s (10 s: 87 s).
NEWTON_HOMOTOPY_BUDGET = 5.0
#: Its cap when the boundary-condition ramp follows (#1018). Its stress
#: rescues take under 2 s (random 60, 91, 136: 0.5-1.8 s for the whole
#: solve); the 5 s it spent on draws it did not rescue was what the ramp then
#: lacked to finish random 119 within the harness's 20 s. Not 2 s: a step's
#: timeout is the cap's remainder, and a solve runs its LM fallback only with
#: more than 2 s left, which random 91's last step needs.
NEWTON_HOMOTOPY_BUDGET_BEFORE_RAMP = 3.0

#: Boundary-condition ramp retry (#642, ``NetworkSolver.bc_ramp_retry``): the
#: gentle end it first solves cold, halved down to BC_RAMP_LAM_MIN until one
#: does, and its solve cap. On the stress random draws the homotopy left
#: unsolved it rescued 3 of 21 in 6-12 solves (15-34 s).
BC_RAMP_LAM0 = 0.05
BC_RAMP_LAM_MIN = 1.0e-3
BC_RAMP_MAX_SOLVES = 40

#: Element attributes that only warm-start a kernel's own inner solve (the
#: Fanno and Fanno-Rayleigh march). A retry that starts cold drops them.
_KERNEL_WARM_STARTS = ("_last_G", "_last_M_exit")

#: Resistance-network cold start (``NetworkSolver.resistance_seed``): the
#: damped fixed-point iteration's cap and relaxation, the flow floor [kg/s]
#: that keeps a zero-flow element's conductance finite, and the fraction of
#: the network's smallest resistance a lossless connection gets. On 160
#: stress random draws it took the cold path from 67 to 110 converged.
SEED_MAX_ITER = 200
SEED_RELAX = 0.5
SEED_MDOT_FLOOR = 1.0e-4
SEED_LOSSLESS_FRACTION = 1.0e-2
#: Passes that re-fit a compressible element's seed K to the flow its own
#: kernel passes at the seed's pressures (#971), and the change in ln K under
#: which they stop.
SEED_COMPRESSIBLE_PASSES = 4
SEED_COMPRESSIBLE_RTOL = 0.05
#: Seed flux, as a fraction of the sonic flux at the upstream Pt, from which
#: an element's seed K is re-fitted to its compressible kernel (about M 0.2).
SEED_COMPRESSIBLE_FLUX_GATE = 0.3
#: Times a fraction loss may change state in the seed's Picard loop (ratio or
#: resistance; #989) before it keeps the one it is in.
SEED_FRACTION_MAX_SWITCHES = 8
_SEED_NODE_TYPES = (
    PlenumNode,
    MomentumChamberNode,
    CombustorNode,
    PressureBoundary,
    MassFlowBoundary,
    WallNode,
)

#: At most this many exceptions keep their full record per solve; the rest are
#: counted (``__caught_exception_count__``).
_MAX_CAUGHT_RECORDS = 50

#: RESIDUAL-ROW-SCALING check (``NetworkSolver.check_row_kinds``): a row's
#: units are those of the unknown family carrying an EXACT unit coefficient,
#: | |J_ij| - 1 | <= this. A residual written X - f(...) carries +-1 on its
#: own unknown (a Pa row on a pressure, a kg/s row on a flow). The first
#: rule, any coefficient in [0.5, 2], took a coincidental gain for one: a
#: 1 um compressible crossflow segment's kg/s row has dm/dPt ~ 1.7 kg/(s Pa)
#: and no flow coefficient near 1, and read as Pa. Exact, on 349k suite rows
#: and 5.5k stress rows (tmp/band_sum.py): 0 false DISAGREE at any tolerance
#: from 1e-1 to 1e-9, and more rows decided than the band (unclear 2.5k ->
#: 0.5k on the suite).
_ROW_KIND_UNIT_TOL = 1e-3


class RowKindMismatchError(TypeError):
    """A residual row's declared units disagree with its Jacobian
    (RESIDUAL-ROW-SCALING). A code bug, raised by strict solves only."""


def _classify_solver_message(message: str) -> SolveOutcome:
    """The one place SciPy's free-text failure message is interpreted."""
    text = (message or "").lower()
    if "not making good progress" in text:
        return SolveOutcome.NO_PROGRESS
    return SolveOutcome.NOT_CONVERGED

    pass


class SolverStallError(SolverTimeoutError):
    """Raised when a solve phase's best residual plateaus far from the root.

    Subclasses SolverTimeoutError so every existing handler treats it as
    "this phase is over, keep the best iterate". In the hybr phase the
    LM fallback then runs with all the remaining budget instead of
    waiting for hybr to burn its fixed share; in the LM fallback phase
    the primary attempt fails fast so the outlet-referenced warm-start
    auto-retry gets the budget instead of a doomed LM grind.
    """

    pass


# Stall-handover tuning. Validated offline against recorded convergence
# histories (tmp/stall_detector_sim.py) and live (tmp/stall_handover_
# verify.py, 2026-07-05): with a 5 s window and a 10% relative-improvement
# floor the detector fired at 15.5 s on a 16 kg/s 20 bar 5-tee network
# where hybr plateaus at |F| ~ 1e4 and LM-from-cold converges in ~6 s
# (103.5 s -> 23.0 s total), and never fired on networks where hybr
# converges steadily (one-tee 0.4, 5-tee 0.18). The far-from-tolerance
# gate (100x tol) exists because slow-but-steady grinding NEAR the root
# is not a stall: the outlet-ref seeded retry on the 5-tee at 0.18 kg/s
# crawls at |F| ~ 1e-2 for a while and then converges -- cutting it over
# to LM (which cannot finish that network) turned a success into a
# failure during live verification.
_STALL_WINDOW_S = 5.0
_STALL_MIN_REL_IMPROVEMENT = 0.10
_STALL_FAR_FACTOR = 100.0
# Below the far gate a crawl is not a stall, but a DEAD plateau is: the best
# |F| not improved at all over this many windows. The 5-row hot impingement
# array (#605) parks LM at |F| = 0.0718, just under the gate, with not one
# better point in 175 s; it was only caught when its plateau happened to sit
# at 0.138, above it. Faster thermo moved it below (#626).
_STALL_DEAD_WINDOWS = 3
# ... and a CRAWL is: under this relative improvement across those windows.
# With the duct heat transfer at the duct's static state (#636) the same array
# no longer parks dead but creeps, 0.5% every 3 s at |F| ~ 0.075 (1.7% per
# 15 s), for 200 s of LM before the warm-start retry solved it in 0.3 s. LM
# converging near a root is quadratic, orders of magnitude faster than this.
_STALL_NEAR_MIN_REL = 0.05
# Far from the root a crawl that beats the 10% floor is not a stall either,
# but it can still never arrive: LM on the jobson1955 jet_chain stress net
# improves just over 10% per 5 s at |F| ~ 2 and times out at 20 s, while the
# warm-start retry solves it in 6 s (#864). So far from tol, a phase also
# hands over when the time to tol projected from its log-rate over the window
# exceeds _STALL_PROJ_FACTOR x the attempt's remaining budget (or a fixed
# horizon without a timeout). Measured on 14 stress families x 80 seeds under
# both orifice models plus the audit set: no success lost, one gained, 2-3%
# less time; a 2.5 s window lost one, so the window stays at 5 s.
_STALL_PROJ_FACTOR = 3.0
_STALL_PROJ_HORIZON_S = 60.0
# A heated node's share of its heat taken up (cb.mixer_heat_per_unit_flow) is
# 1 at or above MIXER_HEAT_MDOT_FLOOR and 0 when stagnant. A root strictly
# between, by more than this, sits in the fade below the floor: one the fade
# made (in a driven net, the hotter the node the less it flows), not the
# network's. A dead end's root takes up < 1e-12 of its heat (#626).
_HEAT_FADE_TOL = 1e-6
# A compressible momentum chamber's closure Pt = P0(P, T0, |m|/A) had two
# static roots for one (Pt, m): the subsonic one and a supersonic one below the
# sonic P*. No element here accelerates a stream past Mach 1 (a constant-area
# channel chokes at 1; an orifice discharges at the chamber's static P), so a
# chamber on the supersonic root was an artifact of the root finder, not a
# flow. Since #788 the closure is continued past P* so that each Pt has one
# root; a chamber past its sonic point (Mach above 1 by more than this) has a
# Pt below the sonic Pt* of its flux: a choked-bounded junction, or refused.
_SUPERSONIC_CHAMBER_TOL = 1e-6


def _stall_detected(
    samples: list[tuple[float, float]],
    window: float = _STALL_WINDOW_S,
    min_rel: float = _STALL_MIN_REL_IMPROVEMENT,
    tol: float = 1e-3,
    far_factor: float = _STALL_FAR_FACTOR,
) -> bool:
    """Return True when the best-|F| trace has plateaued far from the root.

    Phase-agnostic: used for both the hybr phase and the LM fallback
    phase (the LM phase feeds phase-relative times so the grace period
    applies from the start of the phase, not the solve).

    ``samples`` is the per-eval history of ``(t_elapsed_s, best_norm)``
    with best_norm non-increasing. Stall = the newest best norm is still
    far above tolerance (``> far_factor * tol``) and improved by less
    than ``min_rel`` relative to the best norm from ``window`` seconds
    ago; or, nearer the root, it has improved by less than
    ``_STALL_NEAR_MIN_REL`` over ``_STALL_DEAD_WINDOWS`` windows (a dead
    plateau or a crawl).
    Histories shorter than the window never stall (grace period).
    """
    if not samples:
        return False
    t_now, best_now = samples[-1]
    if t_now < window:
        return False
    if best_now <= far_factor * tol:
        # Near the root: a dead plateau or a crawl over several windows.
        t_dead = t_now - _STALL_DEAD_WINDOWS * window
        if t_dead < 0.0:
            return False
        for tv, bv in reversed(samples):
            if tv <= t_dead:
                return bv <= 0.0 or (bv - best_now) / bv < _STALL_NEAR_MIN_REL
        return False
    best_then = None
    for tv, bv in reversed(samples):
        if tv <= t_now - window:
            best_then = bv
            break
    if best_then is None or best_then <= 0.0:
        return False
    return (best_then - best_now) / best_then < min_rel


def _crawl_too_slow(
    samples: list[tuple[float, float]],
    remaining: float | None,
    window: float = _STALL_WINDOW_S,
    tol: float = 1e-3,
    far_factor: float = _STALL_FAR_FACTOR,
) -> bool:
    """Return True when a far-from-root crawl cannot reach tol in time.

    ``samples`` as for ``_stall_detected``. The projected time to ``tol`` is
    ``ln(best_now / tol) / r`` with ``r = ln(best_then / best_now) / window``
    the log-rate over the last window; it fires when that exceeds
    ``_STALL_PROJ_FACTOR`` x ``remaining`` seconds of the attempt
    (``_STALL_PROJ_HORIZON_S`` when the attempt has no timeout). Near the root
    (``<= far_factor * tol``) it never fires: convergence there turns
    quadratic and a log-rate projection says nothing (#864).
    """
    if not samples:
        return False
    t_now, best_now = samples[-1]
    if t_now < window or best_now <= far_factor * tol:
        return False
    best_then = None
    for tv, bv in reversed(samples):
        if tv <= t_now - window:
            best_then = bv
            break
    if best_then is None or best_then <= 0.0:
        return False
    rate = math.log(best_then / best_now) / window
    budget = _STALL_PROJ_HORIZON_S if remaining is None else max(remaining, 0.0)
    if rate <= 0.0:
        return True
    return math.log(best_now / tol) / rate > _STALL_PROJ_FACTOR * budget


class NetworkSolver:
    """
    Orchestrates the numerical solution of a fluid flow network using scipy.optimize.root.
    Translates the graph topology into a flat array of unknowns and residuals.
    Nodes define their Pressure (P) and Total Pressure (Pt) as unknowns,
    while Temperature (T) and Composition (Y) are forward-propagated using
    topological ordering and automatic differentiation (chain-rule relay).
    """

    SUPPORTED_METHODS = [
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
    ]

    @staticmethod
    def _is_pressure_unknown(name: str) -> bool:
        return name.endswith(".P") or name.endswith(".Pt")

    @staticmethod
    def _is_mdot_unknown(name: str) -> bool:
        return (
            name.endswith(".m_dot") or name.endswith(".m_dot_com") or name.endswith(".m_dot_branch")
        )

    @staticmethod
    def _to_incompressible_regime(regime: str) -> str | None:
        if regime == "compressible":
            return "incompressible"
        return None

    def _compressible_element_overrides(self) -> dict[str, str]:
        """Map compressible element IDs to their incompressible fallback regimes."""
        overrides: dict[str, str] = {}
        for elem_id, element in self.network.elements.items():
            regime = getattr(element, "regime", None)
            if not isinstance(regime, str):
                continue
            fallback = self._to_incompressible_regime(regime)
            if fallback is not None:
                overrides[elem_id] = fallback
        return overrides

    def _reference_scales(self, x0_use: np.ndarray) -> tuple[float, float]:
        pressure_vals = [
            abs(float(v))
            for n, v in zip(self.unknown_names, x0_use, strict=False)
            if self._is_pressure_unknown(n)
        ]
        mdot_vals = [
            abs(float(v))
            for n, v in zip(self.unknown_names, x0_use, strict=False)
            if self._is_mdot_unknown(n)
        ]

        ref_p_abs = float(np.median(pressure_vals)) if pressure_vals else 1.0e5
        ref_mdot = float(np.median(mdot_vals)) if mdot_vals else 1.0

        # Pressure-residual scale: use boundary DP (not absolute P)
        # to match actual residual magnitudes in pressure equations
        if len(pressure_vals) >= 2:
            p_spread = max(pressure_vals) - min(pressure_vals)
            ref_p_residual = max(p_spread, ref_p_abs * 1e-4)
        else:
            ref_p_residual = ref_p_abs * 0.01

        ref_p_residual = max(ref_p_residual, 1.0)
        ref_mdot = max(ref_mdot, 1e-6)
        return ref_p_residual, ref_mdot

    def _build_residual_scales(self, x0_use: np.ndarray) -> np.ndarray:
        """Per-row residual scales D_f, one per residual row, by the row's UNITS.

        RESIDUAL-ROW-SCALING -- the canonical description; every place that
        declares a row's units or compares Jacobians points here.

        What the solver iterates on:
            F_scaled[i] = F_real[i] / D_f[i]            (rows,    see inv_D_f)
            x_real[j]   = x_scaled[j] * D_x[j]          (columns, D_x = |x0|)
            J_scaled    = diag(1/D_f) @ J_real @ diag(D_x)
        hybr/LM accept or reject steps on ||F_scaled||, so D_f decides how
        much each equation counts.

        D_f[i] is chosen by what row i MEASURES, not by which object owns it:
            "p"    -> ref_p    [Pa]   (boundary pressure spread)
            "mdot" -> ref_mdot [kg/s] (median flow unknown)
            "T"    -> ref_T    [K]    (_reference_T_scale: an
                                       EnergyBoundary Tt_target row, #792)
        ref_p / ref_mdot is ~1e5, so a row filed under the wrong kind is
        weighted ~1e5 too heavily or too lightly and Newton stalls. This has
        happened repeatedly:
          * MPCE impulse rows filed as "mdot" (2026-07): 187 s stalls.
          * the compressible ChannelElement switched to a mass-flow residual
            (m - m_calc) while still filed as "p" (#481): 697 evaluations
            and a failed solve where 15 sufficed.

        Who declares what (in residual-row order):
          * nodes: their own rows "p", then the mass-conservation row "mdot"
            (skipped for junction ports);
          * junctions (MultiPortChamberBase): row_scale_kinds(), one per row;
          * every other element: ``residual_scale_kind``, a class attribute,
            or a property when the units follow the configuration
            (ChannelElement: "mdot" when compressible, "p" otherwise).
        There is NO default. ``declared_row_kinds`` (components.py) raises at
        network validation unless the kind is declared on the class that
        defines ``residuals`` or a subclass of it, so a subclass that rewrites
        its residual must re-declare its units. Until then a silent "mdot"
        default filed PressureLoss, AreaChange, BorderCarnotLoss (removed in 0.8.0) and Vortex
        rows (all Pa) 1e5-1e9 too heavily (#614).
        Changing an element's residual form means changing its declared kind
        in the same commit. A strict solve (the test suite) then checks every
        declaration against the Jacobian at the converged solution
        (``check_row_kinds``) and raises ``RowKindMismatchError`` on a
        misfiled row: with #616 reverted, 51 suite tests fail naming the
        element, where the suite had passed with the bug in place.

        Testing Jacobians against finite differences has the same trap: a row
        mixes kg/s and Pa columns, so compare J_ij |x_j| (scaled variables),
        never raw entries normalised by the row's largest -- that hid a 100%
        error in a Pa column next to d/dm = 1 (#481).
        """
        ref_p, ref_mdot = self._reference_scales(x0_use)
        _row_scale_kind_map = {"p": ref_p, "mdot": ref_mdot, "T": self._reference_T_scale()}
        return np.asarray(
            [_row_scale_kind_map[kind] for _owner, kind in self._row_kinds(x0_use)], dtype=float
        )

    def _reference_T_scale(self) -> float:
        """Row scale [K] of a temperature row (an EnergyBoundary Tt_target,
        #792): the spread of the boundary and target total temperatures, at
        least 1% of their largest and at least 1 K -- the temperature analogue
        of ref_p (boundary pressure spread with a relative floor)."""
        temps = [
            float(n.Tt)
            for n in self.network.nodes.values()
            if isinstance(n, PressureBoundary | MassFlowBoundary)
            and math.isfinite(float(getattr(n, "Tt", math.nan)))
        ]
        temps += [
            float(eb.Tt_target)
            for n in self.network.nodes.values()
            for eb in target_energy_boundaries(n)
        ]
        if not temps:
            return 1.0
        return max(max(temps) - min(temps), 1e-2 * max(abs(t) for t in temps), 1.0)

    def _row_kinds(self, x: np.ndarray) -> list[tuple[str, str]]:
        """(owner, declared kind) per residual row, in residual-row order.

        RESIDUAL-ROW-SCALING: the one list both the row scales and the
        strict-mode check (``check_row_kinds``) read, so the check tests
        exactly the units the solver weights by.
        """
        rows: list[tuple[str, str]] = []

        for node_id, node in self.network.nodes.items():
            if isinstance(node, PressureBoundary):
                continue

            state = self._get_node_state(node, x)
            node_res, _ = node.residuals(state)
            # The node's own rows are pressure rows, except the trailing
            # Tt_target rows of its EnergyBoundaries (temperature, #792) and a
            # FuelBoundary's combustor-Tt row (temperature, #892).
            n_T = len(target_energy_boundaries(node)) + (1 if isinstance(node, FuelBoundary) else 0)
            rows.extend((f"node '{node_id}'", "p") for _ in node_res[: len(node_res) - n_T])
            rows.extend((f"node '{node_id}' Tt_target", "T") for _ in range(n_T))

            # Node mass conservation equation appended in _residuals_and_jacobian,
            # unless the node is a junction port (then the MultiPortChamberBase's
            # sum-mass residual covers conservation; mass row is skipped).
            if not getattr(node, "_is_junction_port", False):
                rows.append((f"node '{node_id}' mass", "mdot"))

        from .components import declared_row_kinds

        for element in self.network.elements.values():
            # MultiPortChamberBase's rows are a MIX of pressure- and
            # mass-flow-magnitude residuals, and which row is which depends
            # on the subclass's own physics (row_scale_kinds()). The base
            # class's own impulse rows, (P_i + rho_i u_i^2) - P_jct, are all
            # PRESSURE-magnitude, followed by one sum-mass row -- but a
            # subclass with different residual physics (e.g. EjectorElement)
            # can have the opposite pattern. Classifying all N+1 rows as
            # ref_mdot (the pre-2026-07-07 behavior) over-weighted the
            # impulse rows by ref_p/ref_mdot (~1e5) in the scaled system
            # hybr/LM optimize -- the root cause of the "doomed primary"
            # stall class on MPCE networks (three-port fixture: cold hybr
            # 187 s stall -> 2.6 s clean convergence with correct scales;
            # 5-tee slow pockets and the 16 kg/s net likewise solve cold
            # with no stall handover). Getting a subclass's row types wrong
            # silently reintroduces the same class of stall, hence
            # row_scale_kinds() being a real per-subclass override point
            # rather than a hard-coded N-then-1 pattern.
            # RESIDUAL-ROW-SCALING: every element declares its rows' units
            # (row_scale_kinds() for a junction, residual_scale_kind
            # otherwise); there is no default to fall back on (see this
            # method's docstring and declared_row_kinds).
            label = f"{type(element).__name__} '{element.id}'"
            rows.extend((label, k) for k in declared_row_kinds(element))

        return rows

    def check_row_kinds(self, x: np.ndarray) -> list[dict[str, Any]]:
        """Each residual row's declared units against its Jacobian at ``x``.

        RESIDUAL-ROW-SCALING. A pressure row (Pa) written X - f(...) has a
        coefficient of exactly +-1 on some pressure unknown; a mass-flow row
        (kg/s) on some flow unknown. Per row, ``measured`` is the family with
        a coefficient within ``_ROW_KIND_UNIT_TOL`` (1e-3) of unit magnitude,
        when exactly one has one (a gain that only happens to be near 1, such
        as dm/dPt of a very short segment, is not exact).
        ``verdict``:
          * ``agree`` / ``DISAGREE``: measured equals / differs from declared;
          * ``unclear``: both or neither family has such a coefficient;
          * ``untestable``: the row has no pressure or no flow unknown.
        Only DISAGREE is evidence of a misfiled row; the strict solve
        (``raise_internal_errors`` / COMBAERO_SOLVER_STRICT) raises
        ``RowKindMismatchError`` on it at every converged solution. Evaluate
        at a converged ``x``: far from a root a coefficient can be off by
        orders of magnitude. Pressure unknowns are named ``.P``, ``.Pt`` or
        ``.P_*``; flow unknowns ``.m_dot*``.
        """
        x = np.asarray(x, dtype=float)
        _, J = self._residuals_and_jacobian(x)
        J = J.tocsr()
        rows = self._row_kinds(x)
        if len(rows) != J.shape[0]:
            raise RowKindMismatchError(
                f"{len(rows)} declared row kinds for {J.shape[0]} residual rows"
            )
        last = [name.rsplit(".", 1)[-1] for name in self.unknown_names]
        is_p = [s in ("P", "Pt") or s.startswith("P_") for s in last]
        is_m = [s.startswith("m_dot") for s in last]
        report: list[dict[str, Any]] = []
        for i, (owner, declared) in enumerate(rows):
            lo, hi = J.indptr[i], J.indptr[i + 1]
            has_p = has_m = p_unit = m_unit = False
            for j, v in zip(J.indices[lo:hi], J.data[lo:hi], strict=True):
                unit = abs(abs(v) - 1.0) <= _ROW_KIND_UNIT_TOL
                if is_p[j]:
                    has_p, p_unit = True, p_unit or unit
                elif is_m[j]:
                    has_m, m_unit = True, m_unit or unit
            measured = "p" if p_unit and not m_unit else "mdot" if m_unit and not p_unit else None
            if declared == "T":
                # A temperature row (#792) is written on a DERIVED Tt, never on
                # an unknown of its own, so no unit coefficient identifies it.
                verdict = "untestable"
            elif measured is not None:
                verdict = "agree" if measured == declared else "DISAGREE"
            else:
                verdict = "unclear" if has_p and has_m else "untestable"
            report.append(
                {
                    "row": i,
                    "owner": owner,
                    "declared": declared,
                    "measured": measured,
                    "verdict": verdict,
                }
            )
        return report

    def _assert_row_kinds(self, sol: dict[str, Any]) -> None:
        """Strict mode: raise if a converged solution shows a misfiled row."""
        if not (self._strict() and sol.get("__converged__") and "__x_solution__" in sol):
            return
        x = np.asarray(sol["__x_solution__"], dtype=float)
        if len(x) != len(self.unknown_names):
            return
        bad = [r for r in self.check_row_kinds(x) if r["verdict"] == "DISAGREE"]
        if bad:
            lines = "; ".join(
                f"row {r['row']} of {r['owner']}: declared {r['declared']!r}, "
                f"Jacobian says {r['measured']!r}"
                for r in bad[:10]
            )
            raise RowKindMismatchError(
                f"RESIDUAL-ROW-SCALING: {len(bad)} row(s) filed under the wrong units "
                f"({lines}). Fix the element's residual_scale_kind / row_scale_kinds."
            )

    def __init__(self, network: FlowNetwork):
        self.network = network
        # The starting-point safety net: None follows SOLVER_STRICT_ENV.
        self.raise_internal_errors: bool | None = None
        # Every exception the solve caught (start, rejected probes), with
        # where it was raised; returned as __caught_exceptions__.
        self._caught_exceptions: list[dict[str, Any]] = []
        self._caught_count = 0
        self._eval_where = ""
        self._eval_phase = "start"
        self.unknown_names: list[str] = []
        # Mapping from node/element ID to index in the flattened x array
        self._unknown_indices: dict[str, list[int]] = {}
        # EnergyBoundary(Tt_target) heat unknowns per node (#792):
        # {node id: [(boundary, global index)]}, set by _build_x0.
        self._target_q_columns: dict[str, list[tuple[EnergyBoundary, int]]] = {}
        # Mapping from unknown name to global index
        self._name_to_index: dict[str, int] = {}
        # Canonical default composition (Mass Fractions)
        self._default_Y: list[float] = list(cb.species.dry_air_mass())
        self._n_species: int = len(self._default_Y)
        self._topological_order: list[str] = []
        self._derived_states: dict[str, tuple[float, list[float], Any]] = {}
        # The gas a two-port element delivers into a PressureBoundary, keyed
        # (boundary id, element id): (Tt, Y), from the element's own stream and
        # wall heat. Reporting only: a PressureBoundary's state is its set
        # Tt/Y, which no residual reads for an outflow, but an element's exit
        # diagnostics (Tt_out, T_out, mach_out) did. Recorded on the final
        # propagation only (_record_boundary_arrivals).
        self._boundary_arrivals: dict[tuple[str, str], tuple[float, list[float]]] = {}
        self._record_boundary_arrivals: bool = False
        # The undiluted gas each film provider discharges into (#694), keyed
        # FILM_INF_PREFIX + provider id: (T, Y), with the previous pass's for
        # a back-edge read, as for the nodes.
        self._film_inf_states: dict[str, tuple[float, list[float]]] = {}
        self._prev_film_inf_states: dict[str, tuple[float, list[float]]] = {}
        # Nodes whose T(h) was at or below 0 K in the last propagation: less
        # enthalpy than the stream holds at 0 K, no gas state.
        self._T_no_gas: set[str] = set()
        # Optional initial-guess overrides keyed by unknown name; consulted
        # by _build_x0 after user-provided initial_guess dicts and before
        # the topological p_guess / mdot_guess propagators. Populated by
        # init_strategy branches (e.g. analytical_pt_prop) and cleared
        # once the solve returns.
        self._init_overrides: dict[str, float] = {}
        # Mach above which an element on an incompressible closure is
        # reported outside its validity (#588; report only).
        from .compressibility import INCOMPRESSIBLE_MACH_LIMIT

        self.incompressible_mach_limit: float = INCOMPRESSIBLE_MACH_LIMIT
        # Opt-in last retry of a cold solve that did not converge: a Newton
        # homotopy from the cold guess (#642). Off by default -- every
        # non-converging solve pays up to its budget -- and turned on by the
        # GUI backend and the stress harness. homotopy_budget [s] caps it within
        # the caller's timeout; None gives it all the time the timeout leaves
        # (the GUI: the user's max solve time is the only limit), or
        # NEWTON_HOMOTOPY_BUDGET when there is no timeout either.
        self.homotopy_retry: bool = False
        self.homotopy_budget: float | None = NEWTON_HOMOTOPY_BUDGET
        # Opt-in retry after the homotopy: continuation in the boundary
        # conditions (#642). Every pressure boundary moves from the mean of
        # their Pt to its own, every flow and heat boundary from zero to its
        # own, each step warm-started from the last. It runs in the time the
        # timeout leaves; the GUI, NetworkRunner and the stress harness turn
        # it on.
        self.bc_ramp_retry: bool = False
        # Cold start from a quadratic-resistance network solve (dPt = K m|m|,
        # mass conserved), where every element offers a seed resistance
        # (``_seed_resistance``); other networks keep the propagated start.
        self.resistance_seed: bool = True
        # solve(stall_handover=False) turns this off for the duration of the
        # call, nested solves and the wall-free seed solver included (#864).
        self._stall_handover: bool = True

    def _flow_components(self) -> list[tuple[list[str], list[str]]]:
        """(node ids, element ids) of each connected flow component, in the
        network's insertion order (#899). Only elements connect: a thermal
        wall, a film provider or a node a FuelBoundary names couples no flow.
        Insertion order, never a set's: a set of ids iterates by hash(str),
        which Python randomises per process (see _propagate_pressure_guess)."""
        net = self.network
        seen: set[str] = set()
        comps: list[tuple[list[str], list[str]]] = []
        for nid in net.nodes:
            if nid in seen:
                continue
            members = net._reachable_from(nid)
            seen |= members
            nodes = [n for n in net.nodes if n in members]
            elems = [
                eid
                for eid, e in net.elements.items()
                if any(n in members for n in (*e.all_source_nodes(), *e.all_sink_nodes()))
            ]
            comps.append((nodes, elems))
        return comps

    @contextlib.contextmanager
    def _component_view(self, node_ids: list[str], elem_ids: list[str]) -> Any:
        """``self.network`` as one flow component while the cold-start
        estimates are made (#899): the same node and element objects in a
        topology of their own, no walls. Topology is not re-resolved."""
        full = self.network
        sub = FlowNetwork()
        for nid in node_ids:
            sub.add_node(full.nodes[nid])
        for eid in elem_ids:
            sub.add_element(full.elements[eid])
        sub.thermal_coupling_enabled = full.thermal_coupling_enabled
        sub.film_providers = full.film_providers
        self.network = sub
        try:
            yield
        finally:
            self.network = full

    def _infer_reference_state(self) -> dict[str, Any]:
        """Derive sensible default values for unknowns from boundary nodes.

        Uses the **minimum** ``PressureBoundary`` pressure so that
        interior nodes start near the downstream end of the network,
        consistent with the previous hard-coded 1 atm default for
        low-pressure networks while giving a physically meaningful
        starting point for high-pressure ones.  Temperature and
        composition are averaged
        across all boundary nodes that carry them.  Mass-flow defaults
        to the total prescribed mass flow divided by the number of flow
        elements.

        Returns a dict with keys ``'P'``, ``'T'``, ``'Y'``, ``'m_dot'``.
        """
        pressures: list[float] = []
        temperatures: list[float] = []
        compositions: list[list[float]] = []
        total_mdot = 0.0

        for node in self.network.nodes.values():
            if isinstance(node, PressureBoundary):
                pressures.append(getattr(node, "Pt", 101325.0))
                T = getattr(node, "Tt", None)
                if T is not None:
                    temperatures.append(T)
                Y = getattr(node, "Y", None)
                if Y is not None:
                    compositions.append(list(Y))

            elif isinstance(node, MassFlowBoundary):
                T = getattr(node, "Tt", None)
                if T is not None:
                    temperatures.append(T)
                Y = getattr(node, "Y", None)
                if Y is not None:
                    compositions.append(list(Y))
                # A FuelBoundary's flow is solved, not given (#892): its
                # m_dot here is whatever the last start left in it, and
                # counted it set the reference from that (#899): the C3X
                # burner's cold solve 21 s, without it 15 s.
                if not isinstance(node, FuelBoundary):
                    total_mdot += getattr(node, "m_dot", 0.0)

        ref_P = float(min(pressures)) if pressures else 101325.0
        ref_T = float(np.mean(temperatures)) if temperatures else 300.0
        if compositions:
            ref_Y = [float(v) for v in np.mean(compositions, axis=0)]
        else:
            ref_Y = list(self._default_Y)

        n_elems = max(len(self.network.elements), 1)
        if total_mdot > 0:
            ref_mdot = total_mdot / n_elems
        else:
            ref_mdot = self._bernoulli_reference_mdot(pressures, ref_P, ref_T, ref_Y)

        return {"P": ref_P, "T": ref_T, "Y": ref_Y, "m_dot": ref_mdot}

    def _bernoulli_reference_mdot(
        self,
        pressures: list[float],
        ref_P: float,
        ref_T: float,
        ref_Y: list[float],
    ) -> float:
        """Reference mass flow for a network driven only by pressures.

        With no ``MassFlowBoundary`` anywhere, nothing tells the solver the
        scale of the flow, and this used to return a hard 0.1 kg/s. The level
        is then set entirely by the imposed pressure differences, which is
        exactly what a Bernoulli estimate reads off:

            m = A sqrt(2 rho dP)

        the same estimate ``_propagate_analytical_pt_prop`` already makes for a
        ``ChannelElement``; it was simply never made for the reference level.

        Measured on the random-boundary-condition sweep
        (``validation/junction/random_robustness.py``), over junctions driven
        by three pressure boundaries: the implied level spans 2e-3 to 36 kg/s,
        the fixed 0.1 was off by more than a decade in 27 of 60 draws, and
        replacing it with this estimate takes convergence from 80% to 93% with
        stalls falling from 18% to 3%.

        Falls back to the old constant when the network offers no area or no
        pressure spread to work from.
        """
        if len(pressures) < 2:
            return 0.1
        # Per ELEMENT, not the whole network's spread. The same step
        # `_propagate_pressure_guess` uses for its BFS: in a chain the total
        # drop is shared out, and charging all of it to one element
        # overestimates the flow by roughly the square root of the chain
        # length. Measured: the global spread broke a bypass scenario that
        # converges either side of this change.
        n_elems = max(len(self.network.elements), 1)
        dp = float(max(pressures) - min(pressures)) / n_elems
        if dp <= 0.0:
            return 0.1
        areas = [
            float(a)
            for a in (
                getattr(obj, "area", None)
                for obj in (*self.network.nodes.values(), *self.network.elements.values())
            )
            if a is not None and float(a) > 0.0
        ]
        if not areas:
            return 0.1
        try:
            rho = float(cb.density(ref_T, ref_P, list(cb.mass_to_mole(ref_Y))))
        except Exception:  # noqa: BLE001 -- a reference estimate must not raise
            return 0.1
        if rho <= 0.0:
            return 0.1
        # The median area, so one unusually small or large port does not set
        # the scale for the whole network.
        area = float(np.median(areas))
        return max(area * math.sqrt(2.0 * rho * dp), 1e-9)

    def _component_drive(self, node_ids: list[str], elem_ids: list[str]) -> str:
        """Why a flow component is driven, from its INPUTS only (#1004), or ""
        when nothing drives it. Never reads an iterate.

        - (a) its PressureBoundary Pt differ by more than
          _FRACTION_DRIVE_PT_RTOL of the largest;
        - (b) a MassFlowBoundary (not a fuel controller) feeds a nonzero m_dot;
        - (c) a FuelBoundary's target asks for fuel: Tt_target above the
          coldest other boundary Tt in the component (the fuel flow itself is
          an unknown, so its target is read; which boundary feeds the
          combustor is not known before the solve, so the coldest);
        - (d) a compressor or turbine at a pressure ratio other than 1 sets a
          Pt ratio whatever its flow.

        Nothing else raises Pt from its own inputs: a free vortex conserves
        it (#661), junctions, ejectors and tees redistribute what the
        boundaries supply, a splitter valve or metered outlet only divides a
        flow, heat release lowers Pt and there is no gravity term."""
        net = self.network
        pts: list[float] = []
        temps: list[float] = []
        fuels: list[FuelBoundary] = []
        for nid in node_ids:
            n = net.nodes[nid]
            if isinstance(n, PressureBoundary):
                pts.append(float(n.Pt))
                if getattr(n, "Tt", None) is not None:
                    temps.append(float(n.Tt))
            elif isinstance(n, FuelBoundary):
                fuels.append(n)
            elif isinstance(n, MassFlowBoundary):
                if getattr(n, "Tt", None) is not None:
                    temps.append(float(n.Tt))
                if float(n.m_dot) != 0.0:
                    return f"MassFlowBoundary '{nid}' feeds {float(n.m_dot):g} kg/s"
        if len(pts) >= 2 and max(pts) - min(pts) > _FRACTION_DRIVE_PT_RTOL * max(map(abs, pts)):
            return f"pressure boundaries at {min(pts):g} to {max(pts):g} Pa"
        for f in fuels:
            if temps and float(f.Tt_target) > min(temps):
                return f"FuelBoundary '{f.id}' targets Tt {float(f.Tt_target):g} K"
        for eid in elem_ids:
            e = net.elements[eid]
            if isinstance(e, TurbomachineElement) and float(e.pressure_ratio) != 1.0:
                return f"{type(e).__name__} '{eid}' at pressure ratio {float(e.pressure_ratio):g}"
        return ""

    def _dead_end(self, eid: str) -> bool:
        """True when the element's flow is zero by mass balance, from the
        graph alone (#1004): the part of the network on one side of it,
        without it, reaches no node that passes flow -- no PressureBoundary,
        no fuel controller, no MassFlowBoundary feeding a nonzero m_dot. A
        closed restriction (#1024) passes none: the search stops at it."""
        net = self.network
        el = net.elements[eid]
        touching: dict[str, list[str]] = {}
        for oid, other in net.elements.items():
            if oid == eid or _is_closed_restriction(other):
                continue
            for n in (*other.all_source_nodes(), *other.all_sink_nodes()):
                touching.setdefault(n, []).append(oid)

        def opens(nid: str) -> bool:
            n = net.nodes[nid]
            if isinstance(n, (PressureBoundary, FuelBoundary)):
                return True
            return isinstance(n, MassFlowBoundary) and float(n.m_dot) != 0.0

        for start in (el.from_node, el.to_node):
            seen: set[str] = set()
            queue = [start]
            found = False
            while queue and not found:
                nid = queue.pop()
                if nid in seen:
                    continue
                seen.add(nid)
                if opens(nid):
                    found = True
                    break
                for oid in touching.get(nid, []):
                    other = net.elements[oid]
                    queue.extend(
                        n
                        for n in (*other.all_source_nodes(), *other.all_sink_nodes())
                        if n not in seen
                    )
            if not found:
                return True
        return False

    def _fraction_unanchored(self, sol: dict[str, Any]) -> dict[str, dict[str, Any]]:
        """Per unanchored fraction loss (``PressureLossElement.unanchored_fraction``,
        #1004): its |m_dot| in ``sol``, its component's drive
        (``_component_drive``), whether it is a dead end (``_dead_end``) and
        the flag: 1 where |m_dot| is below FLOW_REVERSAL_M_SCALE, on the blend
        of its forward and reverse forms, in a driven component and not in a
        dead end; 0 otherwise; NaN where its flow is not finite. Only the
        flow is read from the result; the rest is inputs and topology."""
        cand = [
            eid
            for eid, e in self.network.elements.items()
            if isinstance(e, PressureLossElement) and e.unanchored_fraction
        ]
        if not cand:
            return {}
        out: dict[str, dict[str, Any]] = {}
        for node_ids, elem_ids in self._flow_components():
            mine = [eid for eid in cand if eid in elem_ids]
            if not mine:
                continue
            drive = self._component_drive(node_ids, elem_ids)
            for eid in mine:
                v = sol.get(f"{eid}.m_dot")
                try:
                    m = abs(float(v)) if v is not None else math.nan
                except (TypeError, ValueError):
                    m = math.nan
                dead = self._dead_end(eid)
                if not math.isfinite(m):
                    flag = math.nan
                else:
                    flag = 1.0 if (m < FLOW_REVERSAL_M_SCALE and drive and not dead) else 0.0
                out[eid] = {"m_dot": m, "drive": drive, "dead_end": dead, "flag": flag}
        return out

    def _mark_fraction_unanchored(self, sol: dict[str, Any]) -> dict[str, Any]:
        """The #1004 check on a top-level solve's result, once, after every
        retry, homotopy and boundary ramp has restored the user's boundaries:
        sets ``fraction_unanchored`` on each PressureLossElement's result and
        ``__fraction_unanchored__``, and on a failed solve names the flagged
        losses (at the best iterate) in ``__message__``. Reported only: it
        never demotes a solve."""
        try:
            fu = self._fraction_unanchored(sol)
        except Exception:  # noqa: BLE001 -- a diagnostic never breaks a solve
            fu = {}
        diags = sol.get("__element_diag__") or {}
        flagged: list[str] = []
        for eid, info in fu.items():
            diag = diags.get(eid)
            if diag is not None and "fraction_unanchored" in diag:
                diag["fraction_unanchored"] = info["flag"]
                sol[f"{eid}.fraction_unanchored"] = info["flag"]
            if info["flag"] == 1.0:
                flagged.append(eid)
        flagged.sort()
        sol["__fraction_unanchored__"] = flagged
        if flagged and not sol.get("__success__", False):
            sol["__message__"] = (
                f"{sol.get('__message__', '')} PressureLossElement "
                + ", ".join(f"'{eid}'" for eid in flagged)
                + ": a fraction loss without m_ref at ~zero flow (its forward/"
                "reverse blend) at this point in a driven network: it holds dPt = "
                "xi Pt at every flow but zero, so a network that gives it a smaller "
                "drop has no root but zero flow through it. Set its m_ref (design "
                "flow), or use a head loss."
            ).strip()
        return sol

    def _propagate_pressure_guess(self, ref: dict[str, Any]) -> dict[str, float]:
        """Estimate per-node pressures by walking the graph from boundaries.

        Starting from every ``PressureBoundary``, pressures are propagated
        both downstream (subtracting dP) and upstream (adding dP) per
        element via BFS.  ``MassFlowBoundary`` nodes are *not* seeded -
        their pressure is discovered by reverse propagation from known
        outlets, ensuring the pressure gradient points in the correct
        flow direction.
        """
        p_guess: dict[str, float] = {}

        # Seed only PressureBoundary nodes (known pressures).
        # MassFlowBoundary nodes have unknown pressure - let BFS
        # discover them by propagating upstream from known outlets.
        for nid, node in self.network.nodes.items():
            if isinstance(node, PressureBoundary):
                p_guess[nid] = getattr(node, "Pt", ref["P"])

        # Estimated per-element dP: spread across the elements that HAVE a
        # drop, or 0.1 % of ref_P. A LosslessConnectionElement imposes Pt
        # continuity, so it carries none (#272): applying a hop's dP across one
        # put a junction port above its own supply (a straight outlet 50 Pa
        # over the inlet in Bassett's Fig 7b network), and the junction seed
        # read flows in the wrong direction from it.
        seed_pressures = list(p_guess.values())
        lossy = [
            e
            for e in self.network.elements.values()
            if not isinstance(e, (LosslessConnectionElement, TurbomachineElement))
        ]
        if len(seed_pressures) >= 2:
            dp_est = (max(seed_pressures) - min(seed_pressures)) / max(len(lossy), 1)
        else:
            dp_est = ref["P"] * 0.001

        def hop(elem: Any) -> float:
            return 0.0 if isinstance(elem, LosslessConnectionElement) else dp_est

        # BFS propagation from known nodes.
        # Seed the queue from the dict, NOT from `visited`: iterating a set of
        # node-ID strings follows hash(str), which Python randomises per
        # process, and every BFS hop applies dp_est -- so the propagated x0
        # pressures, and with them the Newton basin, depended on
        # PYTHONHASHSEED. Measured on the junction scorecard (issue #271):
        # mfb_two_pb converged 618 or 629 of 691 depending on the seed, and
        # Bassett's own Fig 7b case landed on a different root in 5 of 10
        # identical processes. Dicts are insertion-ordered, so this is
        # deterministic.
        visited = set(p_guess.keys())
        queue = list(p_guess.keys())
        while queue:
            nid = queue.pop(0)
            for elem in self.network.get_downstream_elements(nid):
                to_id = elem.to_node
                if to_id not in visited:
                    # A compressor or turbine sets its pressure ratio (#893).
                    p_guess[to_id] = (
                        p_guess[nid] * elem.pt_factor
                        if isinstance(elem, TurbomachineElement)
                        else p_guess[nid] - hop(elem)
                    )
                    visited.add(to_id)
                    queue.append(to_id)
            for elem in self.network.get_upstream_elements(nid):
                from_id = elem.from_node
                if from_id not in visited:
                    p_guess[from_id] = (
                        p_guess[nid] / elem.pt_factor
                        if isinstance(elem, TurbomachineElement)
                        else p_guess[nid] + hop(elem)
                    )
                    visited.add(from_id)
                    queue.append(from_id)

        # The walk follows the DECLARED directions, so an element declared
        # against its flow walks a node out of the physical range: a plenum
        # between 1 and 2 bar boundaries, wired 1 bar -> p -> 2 bar, was
        # guessed at 0.5 bar, below both, where every element feeds it and
        # Newton stalls (#554). Driven by pressure boundaries alone, a passive
        # interior node lies within their range, so the guess is clamped
        # there. Not when a mass-flow boundary or an ejector can push a node
        # outside it.
        if len(seed_pressures) >= 2 and not any(
            isinstance(n, MassFlowBoundary) for n in self.network.nodes.values()
        ):
            from .ejector_element import EjectorElement

            # Nor when an ejector or a compressor (#893) can raise a node above them.
            if not any(
                isinstance(e, (EjectorElement, TurbomachineElement))
                for e in self.network.elements.values()
            ):
                lo, hi = min(seed_pressures), max(seed_pressures)
                for nid in p_guess:
                    if not isinstance(self.network.nodes[nid], PressureBoundary):
                        p_guess[nid] = min(max(p_guess[nid], lo), hi)

        # WallNode dead-ends carry no flow, so the correct pressure equals the
        # upstream node pressure (no friction drop).  Override the BFS estimate.
        for nid, node in self.network.nodes.items():
            if isinstance(node, WallNode):
                for elem in self.network.get_upstream_elements(nid):
                    p_guess[nid] = p_guess.get(elem.from_node, p_guess.get(nid, ref["P"]))
                    break

        return p_guess

    def _propagate_mdot_guess(self, ref: dict[str, Any]) -> dict[str, float]:
        """Estimate per-element mass flow by propagating boundary flows.

        Mass flow is seeded from ``MassFlowBoundary`` nodes and propagated
        in topological order. Merge nodes sum incoming flows; split nodes
        divide the available flow equally among downstream elements.

        For compressible elements, the propagated value is capped using
        an isentropic estimate at a conservative target Mach number
        (``M=0.3``), based on the element area and reference state.
        This keeps startup states safely subsonic while preserving
        topology-derived split behavior.

        Returns a dict mapping element ID to estimated ``m_dot`` [kg/s].
        """
        mach_target = 0.3
        mdot_guess: dict[str, float] = {}
        node_flow: dict[str, float] = {
            nid: float(getattr(node, "m_dot", 0.0))
            for nid, node in self.network.nodes.items()
            if isinstance(node, MassFlowBoundary)
        }

        ref_X = list(cb.mass_to_mole(ref["Y"]))
        gamma_ref = float(cb.isentropic_expansion_coefficient(ref["T"], ref_X))
        rho_ref = float(cb.density(ref["T"], ref["P"], ref_X))
        a_ref = float(cb.speed_of_sound(ref["T"], ref_X))

        term = 1.0 + 0.5 * (gamma_ref - 1.0) * mach_target**2
        exponent = -(gamma_ref + 1.0) / (2.0 * (gamma_ref - 1.0))
        mdot_ratio = mach_target * term**exponent

        for nid in self._compute_topological_order():
            node = self.network.nodes[nid]

            if isinstance(node, MassFlowBoundary):
                available = node_flow[nid]
            else:
                available = node_flow.get(nid, ref["m_dot"])

            down_elems = self.network.get_downstream_elements(nid)
            if not down_elems:
                continue

            share = available / len(down_elems)
            # Only substitute a nonzero reference guess for elements whose
            # source is unknown.  If the share is zero because an upstream
            # MassFlowBoundary is genuinely OFF (m_dot=0), keep zero.
            if abs(share) < 1e-12 and not isinstance(node, MassFlowBoundary):
                share = ref["m_dot"]

            from_mass_boundary = isinstance(node, MassFlowBoundary)
            for elem in down_elems:
                elem_share = share
                regime = getattr(elem, "regime", None)

                # Dead-end channels (to_node is a WallNode) carry zero flow.
                # Start Newton at m_dot=0 so the mass-balance constraint is
                # trivially satisfied from the first iteration.
                # So does a closed restriction (#1024), whatever lies beyond.
                if isinstance(
                    self.network.nodes.get(elem.to_node), WallNode
                ) or _is_closed_restriction(elem):
                    elem_share = 0.0

                # Only apply the Mach cap when the flow is estimated (pressure-
                # boundary seeded).  When a MassFlowBoundary seeds the element
                # directly, the mass flow is a known constraint and capping it
                # produces an initial guess that is too far from the solution,
                # causing hybr to converge to a spurious local minimum.
                elif regime == "compressible" and not from_mass_boundary:
                    area = getattr(elem, "area", None)
                    if area is not None and area > 0.0:
                        mdot_cap = rho_ref * a_ref * area * mdot_ratio
                        elem_share = min(abs(elem_share), mdot_cap)
                        if share < 0:
                            elem_share = -elem_share

                # Accumulate rather than overwrite: multi-source elements (e.g.
                # merging tee) are visited once per source node and each visit
                # contributes its share, so the total equals the sum of all inflows.
                mdot_guess[elem.id] = mdot_guess.get(elem.id, 0.0) + elem_share
                node_flow[elem.to_node] = node_flow.get(elem.to_node, 0.0) + elem_share

        return mdot_guess

    def _propagate_analytical_pt_prop(self, ref: dict[str, Any]) -> dict[str, float]:
        """Compute topology-aware initial guesses for MPCE compressible cold-start.

        For each ``ChannelElement`` with a signed pressure gradient across
        it, seed the m_dot unknown via a Bernoulli estimate
        ``m_dot = sign(dP) * A * sqrt(2 * rho_ref * |dP|)``.

        For each ``MultiPortChamberBase`` (and subclasses), seed the
        junction unknown (``P_jct`` or ``Pt_jct``) from the propagated Pt at the junction's
        "common" port (the arm on the single-inlet side for a branch or
        the single-outlet side for a merge). The default
        ``ref["P"] = min(boundary Pt)`` is almost never near the true
        junction static, which is what motivated this strategy (see
        LHS-32 audit, ``tmp/mpce_cold_start_audit.py``).

        Returns a dict mapping unknown name to seed value. Consumed by
        ``_build_x0`` via ``self._init_overrides``.
        """
        from .components import MultiPortChamberBase

        p_guess = self._propagate_pressure_guess(ref)
        if not p_guess:
            return {}
        # Reference density from Pt/RT at the propagated mean pressure and
        # ambient reference temperature. R_air is a stand-in for the local
        # gas; MPCE audit shows this level of approximation is enough to
        # get Newton into the correct basin.
        R_air = 287.0
        Tt = float(ref.get("T", 300.0))
        Pt_ref = float(np.mean(list(p_guess.values())))
        rho_ref = max(Pt_ref / (R_air * Tt), 1e-3)

        overrides: dict[str, float] = {}
        for elem_id, elem in self.network.elements.items():
            if isinstance(elem, ChannelElement):
                pt_from = p_guess.get(elem.from_node)
                pt_to = p_guess.get(elem.to_node)
                if pt_from is None or pt_to is None:
                    continue
                dp = pt_from - pt_to
                if dp == 0.0:
                    continue
                diameter = getattr(elem, "diameter", None)
                if diameter is None or diameter <= 0.0:
                    continue
                area = math.pi * (float(diameter) / 2.0) ** 2
                mdot_mag = area * math.sqrt(2.0 * rho_ref * abs(dp))
                overrides[f"{elem_id}.m_dot"] = math.copysign(mdot_mag, dp)
            elif isinstance(elem, MultiPortChamberBase):
                inlets = list(getattr(elem, "inlet_nodes", []))
                outlets = list(getattr(elem, "outlet_nodes", []))
                if len(inlets) == 1:
                    common = inlets[0]
                elif len(outlets) == 1:
                    common = outlets[0]
                else:
                    continue
                pt = p_guess.get(common)
                if pt is None:
                    continue
                overrides[elem.junction_unknown] = pt
                if getattr(elem, "seeds_ports_by_pressure_split", False):
                    overrides.update(
                        self._junction_split_guess(elem, common, p_guess, rho_ref, ref)
                    )
                # A closure whose P_jct is a STATIC pressure (#517) moves the
                # seed off the common port's Pt by that port's head.
                seed_P_jct = getattr(elem, "seed_P_jct", None)
                common_elem = elem._port_element_ids[elem.port_nodes.index(common)]
                m_common = overrides.get(f"{common_elem}.m_dot")
                if seed_P_jct is not None and m_common is not None:
                    overrides[elem.junction_unknown] = seed_P_jct(
                        pt, float(m_common), rho_ref, elem.port_nodes.index(common)
                    )
        return overrides

    def _junction_split_guess(
        self,
        elem: Any,
        common: str,
        p_guess: dict[str, float],
        rho_ref: float,
        ref: dict[str, Any],
    ) -> dict[str, float]:
        """Seed a junction's port flows so x0 conserves mass and reflects the
        imposed pressures.

        The Bernoulli seed above only reaches `ChannelElement`. Junctions in
        both the validation harness and `gui/backend/graph_builder.py` are
        wired with `LosslessConnectionElement`, which got nothing, so every
        port fell back to the flat reference flow: on a three-port junction
        that is 0.1 kg/s in against 0.2 kg/s out, with a 1:1 split whatever the
        boundary pressures say (issue #271, defect 11).

        What this changes is only the SPLIT and the CONTINUITY. The total is
        taken from whatever the topological propagator already put on the
        common port's element, so no new flow-scale heuristic is introduced
        and a MassFlowBoundary upstream still sets the level. The non-common
        ports are then given a Bernoulli share of that total,
        ``A_i * sqrt(2 rho dPt_i)``, normalised to sum to it -- and every port
        is seeded in its own canonical direction, since the sign mapping
        ``port_mdots[i] = port_signs[i] * outer_mdot[i]`` makes a positive
        outer flow the correct direction at that port.

        Ports whose imposed pressure difference is unknown or zero fall back
        to an area-weighted share, which is the incompressible answer when
        nothing distinguishes them.
        """
        port_nodes = list(getattr(elem, "port_nodes", []))
        port_elems = list(getattr(elem, "_port_element_ids", []))
        if len(port_nodes) != len(port_elems) or common not in port_nodes:
            return {}
        areas = list(getattr(elem, "port_areas", [])) or [1.0] * len(port_nodes)
        if len(areas) != len(port_nodes):
            return {}

        common_i = port_nodes.index(common)
        pt_com = p_guess[common]
        mdot_guess = self._propagate_mdot_guess(ref)
        total = abs(float(mdot_guess.get(port_elems[common_i], ref["m_dot"])))
        if total <= 0.0:
            return {}

        # The element's OWN closure, inverted on the guessed port pressures,
        # where it offers one (#272): a junction's port drops are K_i(split)
        # times the common head, so their ratio fixes the split and either one
        # the level. The Bernoulli share below assumes each port's drop drives
        # its own flow, which a dividing straight leg contradicts -- its total
        # pressure RISES for a lateral fraction under 1/2 (K5 < 0) -- and it
        # seeded q ~ 0.25 whatever the target, on the far side of that sign
        # change, from where Newton collapsed onto a dead lateral. The level is
        # taken from the closure only on a pressure-driven network: an imposed
        # flow sets it otherwise.
        closure_seed = getattr(elem, "closure_split_seed", None)
        if closure_seed is not None:
            dpt = [None if p_guess.get(n) is None else p_guess[n] - pt_com for n in port_nodes]
            # The common flow is KNOWN when the common port's own element ends
            # at a MassFlowBoundary (source or sink): the closure then solves
            # for the split at that flow, and the common port's pressure --
            # only a guess there -- drops out.
            m_imposed = None
            common_elem = self.network.elements.get(port_elems[common_i])
            if common_elem is not None:
                for end in (common_elem.from_node, common_elem.to_node):
                    node = self.network.nodes.get(end)
                    if isinstance(node, MassFlowBoundary):
                        m_imposed = abs(float(node.m_dot))
            prefer = None
            if m_imposed is not None:
                generic = self._bernoulli_shares(
                    port_nodes, areas, common_i, pt_com, p_guess, rho_ref
                )
                prefer = generic.get(next(i for i in range(len(port_nodes)) if i != common_i))
            seed = closure_seed(common_i, dpt, rho_ref, m_common=m_imposed, prefer=prefer)
            if seed is not None:
                shares, m_common = seed
                total = m_common
                out: dict[str, float] = {}
                if port_elems[common_i]:
                    out[f"{port_elems[common_i]}.m_dot"] = total
                for i, share in shares.items():
                    if port_elems[i]:
                        out[f"{port_elems[i]}.m_dot"] = share * total
                return out

        out = {}
        if port_elems[common_i]:
            out[f"{port_elems[common_i]}.m_dot"] = total
        for i, share in self._bernoulli_shares(
            port_nodes, areas, common_i, pt_com, p_guess, rho_ref
        ).items():
            if port_elems[i]:
                out[f"{port_elems[i]}.m_dot"] = share * total
        return out

    @staticmethod
    def _bernoulli_shares(
        port_nodes: list[str],
        areas: list[float],
        common_i: int,
        pt_com: float,
        p_guess: dict[str, float],
        rho_ref: float,
    ) -> dict[int, float]:
        """Each non-common port's share of the common flow, A_i sqrt(2 rho
        dPt_i) normalised; area-weighted where no pressure difference is known."""
        others = [i for i in range(len(port_nodes)) if i != common_i]
        shares = []
        for i in others:
            pt_i = p_guess.get(port_nodes[i])
            dpt = 0.0 if pt_i is None else abs(pt_com - pt_i)
            area = abs(float(areas[i])) or 1.0
            shares.append(area * math.sqrt(2.0 * rho_ref * dpt) if dpt > 0.0 else 0.0)
        if sum(shares) <= 0.0:
            shares = [abs(float(areas[i])) or 1.0 for i in others]
        total = sum(shares)
        return {i: s / total for i, s in zip(others, shares, strict=True)}

    def _resistance_network_guess(
        self, ref: dict[str, Any]
    ) -> tuple[dict[str, float], dict[str, float]] | None:
        """Node pressures and element flows of the network as quadratic
        resistances, for the cold start.

        Each element is dPt = K m|m| with K its incompressible
        ``_seed_resistance`` at the reference density; pressure boundaries are
        fixed, mass-flow boundaries inject (a sink draws) their m_dot, every
        other node conserves mass. A monotone network, so it has one solution;
        damped Picard on the conductances G = 1 / (K |m|) finds it. The
        propagated start sets a node's pressure from its neighbours one element
        at a time, and in a large network of resistances decades apart it
        leaves a duct carrying MPa at the guessed flow -- where 48 of 160
        stress random draws stalled even incompressible and adiabatic.

        None when an element offers no resistance, a node is of a type the
        mass balance here does not cover, or there is no pressure boundary.
        """
        net = self.network
        if not all(isinstance(n, _SEED_NODE_TYPES) for n in net.nodes.values()):
            return None
        fixed = {
            nid: float(n.Pt) for nid, n in net.nodes.items() if isinstance(n, PressureBoundary)
        }
        if not fixed:
            return None
        walls = {nid for nid, n in net.nodes.items() if isinstance(n, WallNode)}
        # A closed restriction (#1024) is a wall in the path: no resistance to
        # fit, no flow; it is seeded at zero below.
        shut = [e.id for e in net.elements.values() if _is_closed_restriction(e)]
        els = [
            e
            for e in net.elements.values()
            if e.from_node not in walls and e.to_node not in walls and e.id not in shut
        ]
        # A compressor or turbine (#893) is a pressure ratio, not a resistance:
        # with one end at a fixed pressure it fixes the other (in a chain, in
        # turn), and it leaves the linear network; its flow is the balance at
        # the end it fixed. With neither end fixed the seed does not apply.
        machines = [e for e in els if isinstance(e, TurbomachineElement)]
        boundary = set(fixed)
        changed = True
        while changed:
            changed = False
            for e in machines:
                if e.from_node in fixed and e.to_node not in fixed:
                    fixed[e.to_node] = fixed[e.from_node] * e.pt_factor
                    changed = True
                elif e.to_node in fixed and e.from_node not in fixed:
                    fixed[e.from_node] = fixed[e.to_node] / e.pt_factor
                    changed = True
        if any(e.from_node not in fixed or e.to_node not in fixed for e in machines):
            return None
        els = [e for e in els if not isinstance(e, TurbomachineElement)]
        if not els:
            return None
        P_ref = float(np.mean(list(fixed.values())))
        X_ref = list(cb.mass_to_mole(ref["Y"]))
        # The flow component's own boundary Tt (#990): D-7156 run 19's Mach
        # 0.7 coolant passage (297 K) seeded at a wall-joined mean of 446 K
        # started its pressures 19% high and took the warm-start retry.
        T_seed = float(ref["T"])
        rho = float(cb.density(T_seed, P_ref, X_ref))
        m_ref = max(float(ref["m_dot"]), SEED_MDOT_FLOOR)
        K_list = []
        for e in els:
            k = e._seed_resistance(rho, P_ref, m_ref)
            if k is None or not math.isfinite(k) or k < 0.0:
                return None
            K_list.append(k)
        K = np.array(K_list)
        if not np.any(K > 0.0):
            return None
        # A fraction loss is a pressure ratio, Pt_dn = Pt_up (1 - xi), at any
        # flow (#989): it fixes the drop and leaves the flow to the rest of the
        # network. Its K at the reference flow, xi P_ref / m_ref^2, has no
        # scale -- stress random 3832186739 started its 117 kg/s through one
        # at 0.56 kg/s. The Picard loop below carries it as that ratio once
        # its resistance would take more than xi Pt_up; below that it stays a
        # resistance (with a design flow m_ref, the dynamic head the loss fades
        # with). Between two fixed ends the flow is not its to set: its K.
        frac: dict[int, tuple[float, float | None]] = {}
        for k, e in enumerate(els):
            law = e._seed_fraction() if isinstance(e, PressureLossElement) else None
            if (
                law is not None
                and law[0] > 0.0
                and not (e.from_node in fixed and e.to_node in fixed)
            ):
                frac[k] = law
        sized = [k for k in range(len(els)) if k not in frac and K[k] > 0.0]
        K[K == 0.0] = SEED_LOSSLESS_FRACTION * (K[sized] if sized else K[K > 0.0]).min()

        free = [nid for nid in net.nodes if nid not in fixed and nid not in walls]
        if not free and not machines:
            return None
        idx = {nid: i for i, nid in enumerate(free)}
        inj = np.zeros(len(free))
        for nid in free:
            node = net.nodes[nid]
            if isinstance(node, MassFlowBoundary):
                sink = net.get_upstream_elements(nid) and not net.get_downstream_elements(nid)
                inj[idx[nid]] += -float(node.m_dot) if sink else float(node.m_dot)
        a_ix = [idx.get(e.from_node) for e in els]
        b_ix = [idx.get(e.to_node) for e in els]
        a_fix = np.array([fixed.get(e.from_node, 0.0) for e in els])
        b_fix = np.array([fixed.get(e.to_node, 0.0) for e in els])

        # A fraction loss's state in the loop: a resistance (2) while the drop
        # it takes stays under xi Pt_up, else a forward (+1) or reverse (-1)
        # ratio, a row of its own with its flow as the unknown. With a design
        # flow m_ref the resistance is the fade 1 - (1 - r^2)^3, K = xi Pt_up
        # (3 - 3 r^2 + r^4) / m_ref^2 at r = |m| / m_ref < 1 (``_m_ref_arg``).
        # Without one it is the reference-flow K: its residual holds a drop
        # inside +-xi Pt_up at zero flow (blending the two directions), where
        # a seed at no flow left random 10 to the warm-start retry (11 s
        # against 0.1 s), and where a lossless loop closes on it.
        mode = dict.fromkeys(frac, 2)
        switches = dict.fromkeys(frac, 0)
        pu0 = np.array(
            [
                max(a_fix[k] if a_ix[k] is None else P_ref, b_fix[k] if b_ix[k] is None else P_ref)
                for k in range(len(els))
            ]
        )

        def picard(
            K: np.ndarray, m: np.ndarray
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
            """(p, m, pa, pb): the quadratic-resistance root for K, from m."""
            p = np.full(len(free), P_ref)
            if not free:
                # Every pressure fixed (by boundaries and machines): each
                # element's flow follows from its own drop, dPt = K m |m|.
                dp = a_fix - b_fix
                return p, np.sign(dp) * np.sqrt(np.abs(dp) / K), a_fix, b_fix
            pa, pb = a_fix, b_fix
            pu = pu0
            for _ in range(SEED_MAX_ITER):
                K_it = K
                if frac:
                    K_it = K.copy()
                    for k, (xi, mr) in frac.items():
                        if mode[k] == 2 and mr is not None:
                            r2 = min(abs(m[k]) / mr, 1.0) ** 2
                            K_it[k] = xi * pu[k] * (3.0 - 3.0 * r2 + r2 * r2) / (mr * mr)
                G = 1.0 / (K_it * (np.abs(m) + SEED_MDOT_FLOOR))
                ratio = [k for k in frac if mode[k] in (1, -1)]
                for k in frac:
                    if mode[k] != 2:
                        G[k] = 0.0
                n_free = len(free)
                L = np.zeros((n_free + len(ratio), n_free + len(ratio)))
                b = np.zeros(n_free + len(ratio))
                b[:n_free] = inj
                for k, g in enumerate(G):
                    i, j = a_ix[k], b_ix[k]
                    if i is not None:
                        L[i, i] += g
                        if j is not None:
                            L[i, j] -= g
                        else:
                            b[i] += g * b_fix[k]
                    if j is not None:
                        L[j, j] += g
                        if i is not None:
                            L[j, i] -= g
                        else:
                            b[j] += g * a_fix[k]
                # A ratio's flow leaves its from node and enters its to node;
                # its row is Pt_dn - (1 - xi) Pt_up = 0 in the flow's direction.
                tied = np.zeros(n_free, dtype=bool)
                for c, k in enumerate(ratio):
                    col = n_free + c
                    i, j = a_ix[k], b_ix[k]
                    for node, sign in ((i, 1.0), (j, -1.0)):
                        if node is not None:
                            L[node, col] += sign
                            tied[node] = True
                    keep = 1.0 - frac[k][0]
                    up, dn = ((i, a_fix[k]), (j, b_fix[k]))[:: mode[k]]
                    for (node, fix), coef in ((dn, 1.0), (up, -keep)):
                        if node is not None:
                            L[col, node] += coef
                        else:
                            b[col] -= coef * fix
                # A free node no element reaches (or one only through walls)
                # keeps the reference pressure.
                diag = np.diag(L)[:n_free].copy()
                lone = np.flatnonzero((diag == 0.0) & ~tied)
                L[lone, lone] = 1.0
                b[lone] = P_ref
                try:
                    z = np.linalg.solve(L, b)
                except np.linalg.LinAlgError:
                    z = None
                if z is None or not np.all(np.isfinite(z)):
                    if frac:
                        # Ratios the rest cannot close (two in parallel, say):
                        # the fraction losses keep their reference-flow K.
                        frac.clear()
                        return picard(K, np.full(len(els), m_ref))
                    return None
                p = z[:n_free]
                pa = np.array([p[i] if i is not None else a_fix[k] for k, i in enumerate(a_ix)])
                pb = np.array([p[j] if j is not None else b_fix[k] for k, j in enumerate(b_ix)])
                m_new = G * (pa - pb)
                for c, k in enumerate(ratio):
                    m_new[k] = z[n_free + c]
                if not np.all(np.isfinite(m_new)):
                    return None
                pu = np.maximum(pa, pb)
                switched = False
                for k, (xi, mr) in frac.items():
                    was = mode[k]
                    if switches[k] > SEED_FRACTION_MAX_SWITCHES:
                        continue  # cycling between states: it keeps this one
                    if was == 2:
                        # The resistance takes more than the ratio would.
                        if abs(pa[k] - pb[k]) > xi * pu[k]:
                            mode[k] = 1 if pa[k] > pb[k] else -1
                    elif was * m_new[k] < 0.0:
                        mode[k] = 2
                    elif mr is not None:
                        if abs(m_new[k]) < mr:
                            mode[k] = 2
                    elif K[k] * m_new[k] * m_new[k] < xi * pu[k]:
                        mode[k] = 2  # it would take less at this flow
                    if mode[k] != was:
                        switches[k] += 1
                        switched = True
                done = np.max(np.abs(m_new - m)) <= 1e-6 * (np.max(np.abs(m_new)) + SEED_MDOT_FLOOR)
                m = SEED_RELAX * m + (1.0 - SEED_RELAX) * m_new
                if done and not switched:
                    m = m_new
                    break
            return p, m, pa, pb

        root = picard(K, np.full(len(els), m_ref))
        if root is None:
            return None
        p, m, pa, pb = root
        # Compressible elements (#971): the incompressible K misses what their
        # residual's own kernel does past low Mach -- a short duct between
        # plenums passes the choked flux on a small Pt difference. An element
        # whose seed flux reaches SEED_COMPRESSIBLE_FLUX_GATE of the sonic flux
        # at its upstream Pt joins the re-fit (and stays in it); each pass sets
        # its K to the drop over the flow its kernel passes at the seed's
        # pressures (``_seed_flow``), then re-solves. Below the gate the
        # incompressible K stands, so a low-Mach network keeps its start.
        kernels = [k for k, e in enumerate(els) if e._seed_flow_applies()]
        Y_ref = [float(v) for v in ref["Y"]]
        active: set[int] = set()
        for _ in range(SEED_COMPRESSIBLE_PASSES if kernels else 0):
            for k in kernels:
                if k in active:
                    continue
                hi = float(max(pa[k], pb[k]))
                area = float(getattr(els[k], "area", 0.0) or 0.0)
                if hi > 0.0 and area > 0.0:
                    sonic = area * float(cb.fanno_sonic_mass_flux(hi, T_seed, X_ref))
                    if abs(m[k]) >= SEED_COMPRESSIBLE_FLUX_GATE * sonic:
                        active.add(k)
            if not active:
                break
            K_new = K.copy()
            for k in active:
                hi, lo = max(pa[k], pb[k]), min(pa[k], pb[k])
                if not hi - lo > 0.0:
                    continue
                m_c = els[k]._seed_flow(float(hi), float(lo), T_seed, X_ref, Y_ref)
                if m_c is not None and m_c > 0.0 and math.isfinite(m_c):
                    K_new[k] = (hi - lo) / (m_c * m_c)
            moved = np.max(np.abs(np.log(K_new / K)))
            K = K_new
            root = picard(K, m)
            if root is None:
                return None
            p, m, pa, pb = root
            if moved <= SEED_COMPRESSIBLE_RTOL:
                break
        # An imposed flow past what the resistances pass at the drawn
        # pressures can put a node below zero; the start stays positive.
        p_floor = 0.5 * min(fixed.values())
        p_map = {nid: max(float(p[idx[nid]]), p_floor) for nid in free}
        m_map = {e.id: float(mi) for e, mi in zip(els, m, strict=True)}
        m_map.update(dict.fromkeys(shut, 0.0))
        for nid, P in fixed.items():
            if nid not in boundary:
                p_map[nid] = P  # set by a machine's pressure ratio
        for e in machines:
            # Its flow: the balance at an end it fixed that is no boundary.
            end = e.to_node if e.to_node not in boundary else e.from_node
            sign = 1.0 if end == e.to_node else -1.0
            out = sum(m_map.get(d.id, 0.0) for d in net.get_downstream_elements(end) if d is not e)
            into = sum(m_map.get(u.id, 0.0) for u in net.get_upstream_elements(end) if u is not e)
            m_map[e.id] = sign * (out - into)
        return p_map, m_map

    def _fuel_out(self, fnode: "FuelBoundary", sol_dict: dict[str, Any]) -> float:
        """The fuel flowing out of the boundary [kg/s]: its unknown at a root,
        and what actually left it where the mass rows do not close."""
        return sum(
            max(float(sol_dict.get(f"{e.id}.m_dot", 0.0)), 0.0)
            for e in self.network.get_downstream_elements(fnode.id)
        )

    def _fuel_inflows(
        self, fnode: "FuelBoundary", x: np.ndarray, sol_dict: dict[str, Any]
    ) -> tuple[list[tuple[float, float, list[float]]], list[float]]:
        """The oxidizer the fuel's combustor burns it with, as one stream
        [(m_dot, Tt, Y)] (empty when none flows), and the fuel's Y.

        From the combustor's whole forward inflow less the injected fuel's mass
        and enthalpy, so it holds for any route the fuel takes to the
        combustor -- a direct injector, a fuel line, a premix plenum the air
        also enters -- as long as all of the fuel reaches it, which
        FuelBoundary.resolve_topology checks."""
        Y_f = list(fnode.Y if fnode.Y is not None else cb.mole_to_mass(cb.species.dry_air()))
        n = len(Y_f)
        m_u, H_u, mY_u = 0.0, 0.0, [0.0] * n
        for e in self.network.elements.values():
            if e.to_node != fnode.combustor:
                continue
            m = float(sol_dict.get(f"{e.id}.m_dot", 0.0))
            if m <= 0.0:
                continue
            st = self._report_state(self.network.nodes[e.from_node], x)
            Y = list(st.Y)
            m_u += m
            H_u += m * float(cb.h_mass(float(st.Tt), list(cb.mass_to_mole(Y))))
            mY_u = [a_ + m * y for a_, y in zip(mY_u, Y, strict=True)]
        # The fuel that flows out of the boundary: its unknown at a root, and
        # what actually left it at a point where the mass rows do not close.
        m_f = self._fuel_out(fnode, sol_dict)
        m_ox = m_u - m_f
        if m_ox <= max(1e-12, 1e-9 * m_u):
            return [], Y_f
        Y_ox = [max(a_ - m_f * y, 0.0) / m_ox for a_, y in zip(mY_u, Y_f, strict=True)]
        total = sum(Y_ox)
        Y_ox = [y / total for y in Y_ox]
        h_f = float(cb.h_mass(float(fnode.Tt), list(cb.mass_to_mole(Y_f))))
        X_ox = list(cb.mass_to_mole(Y_ox))
        T_ox = float(cb.calc_T_from_h_mass((H_u - m_f * h_f) / m_ox, X_ox, T_guess=float(fnode.Tt)))
        return [(m_ox, T_ox, Y_ox)], Y_f

    @staticmethod
    def _as_stream(m: float, T: float, Y: list[float]) -> Any:
        st = cb.Stream()
        st.set_T(T)
        st.set_P(101325.0)
        st.set_X(list(cb.mass_to_mole(Y)))
        st.set_mdot(m)
        return st

    def _fuel_peak_T(
        self,
        fnode: "FuelBoundary",
        x: np.ndarray,
        sol_dict: dict[str, Any],
        streams: list[tuple[float, float, list[float]]] | None = None,
    ) -> float:
        """The combustor's Tt with exactly the stoichiometric injected fuel:
        the peak of the reachable window (complete combustion). ``streams``
        replaces the oxidizer at the point."""
        at_point, Y_f = self._fuel_inflows(fnode, x, sol_dict)
        return self._fuel_stoich_T(fnode, at_point if streams is None else streams, Y_f)

    def _fuel_stoich_T(
        self,
        fnode: "FuelBoundary",
        streams: list[tuple[float, float, list[float]]],
        Y_f: list[float],
    ) -> float:
        """``_fuel_peak_T`` for given oxidizer ``streams`` and fuel ``Y_f``."""
        if not streams:
            return float("nan")
        try:
            ox = cb.mix([self._as_stream(m, T, Y) for m, T, Y in streams])
            Y_ox = list(cb.mole_to_mass(list(ox.X)))
            o2_spare = Y_ox[list(cb.species.names).index("O2")] - float(
                cb.oxygen_required_per_kg_mixture(list(ox.X))
            )
            need_f = float(cb.oxygen_required_per_kg_mixture(list(cb.mass_to_mole(Y_f))))
            if o2_spare <= 0.0 or need_f <= 0.0:
                return float("nan")
            m_st = float(ox.mdot) * o2_spare / need_f
            mixed = cb.mix([ox, self._as_stream(m_st, float(fnode.Tt), Y_f)])
            return float(cb.complete_combustion(mixed.T, list(mixed.X), mixed.P).T)
        except (ValueError, RuntimeError):
            return float("nan")

    def _fuel_reactant_view(
        self, fnode: "FuelBoundary", x: np.ndarray, sol_dict: dict[str, Any]
    ) -> tuple[float, float]:
        """(phi, Tt with no injected fuel) at the combustor, from the unburnt
        reactants (#892).

        phi = O2 required by the reactants' combustible species / O2 they
        carry (oxygen_required_per_kg_mixture over Y_O2). The products'
        elemental phi would count the H of water vapour in humid air as fuel
        (0.026 at 60% RH). The Tt with no injected fuel is the oxidizer
        inflows burnt alone -- above their mixing T when they carry fuel."""
        streams, Y_f = self._fuel_inflows(fnode, x, sol_dict)
        if not streams:
            return float("nan"), float("nan")
        try:
            names = list(cb.species.names)
            i_o2 = names.index("O2")
            m_f = self._fuel_out(fnode, sol_dict)
            m_tot = sum(m for m, _, _ in streams) + m_f
            Y_u = [
                (sum(m * Y[k] for m, _, Y in streams) + m_f * Y_f[k]) / m_tot
                for k in range(len(names))
            ]
            o2_have = Y_u[i_o2]
            o2_need = float(cb.oxygen_required_per_kg_mixture(list(cb.mass_to_mole(Y_u))))
            phi = o2_need / o2_have if o2_have > 0.0 else float("inf")
            T_without = self._fuel_T_without(streams)
        except (ValueError, RuntimeError):
            return float("nan"), float("nan")
        return phi, T_without

    def _fuel_T_without(self, streams: list[tuple[float, float, list[float]]]) -> float:
        """The oxidizer streams burnt alone: their mixing T, or above it when
        they carry fuel of their own (complete combustion)."""
        ox = cb.mix([self._as_stream(m, T, Y) for m, T, Y in streams])
        if float(cb.oxygen_required_per_kg_mixture(list(ox.X))) > 0.0:
            return float(cb.complete_combustion(ox.T, list(ox.X), ox.P).T)
        return float(ox.T)

    def _fuel_window_reason(
        self,
        fnode: "FuelBoundary",
        streams: list[tuple[float, float, list[float]]],
        Y_f: list[float],
        T_without: float,
        where: str,
    ) -> str:
        """Why a fuel controller's target is outside its reachable window
        with this oxidizer (#892): not above the Tt with no injected fuel, or
        above the stoichiometric peak; "" when it is inside."""
        fid = fnode.id
        if math.isfinite(T_without) and fnode.Tt_target <= T_without:
            T_mix = float(cb.mix([self._as_stream(m, T, Y) for m, T, Y in streams]).T)
            carries = " (the oxidizer carries fuel of its own)" if T_without > T_mix + 1.0 else ""
            return (
                f"FuelBoundary '{fid}': the target {fnode.Tt_target:.1f} K is not "
                f"above what the combustor reaches with no injected fuel, "
                f"{T_without:.1f} K{carries}{where}: no reacting root (it would need fuel as "
                "a coolant, combustion efficiency ~ 0)."
            )
        T_peak = self._fuel_stoich_T(fnode, streams, Y_f)
        if math.isfinite(T_peak) and fnode.Tt_target > T_peak:
            return (
                f"FuelBoundary '{fid}': the target {fnode.Tt_target:.1f} K "
                f"is above the peak, {T_peak:.1f} K at stoichiometric: the hottest "
                f"this fuel burns in this oxidizer (complete combustion){where}."
            )
        return ""

    def _fuel_guess(
        self, node: "FuelBoundary", mdot_guess: dict[str, float], ref: dict[str, Any]
    ) -> float:
        """A FuelBoundary's start on its branch (#892): the fuel that takes the
        combustor's other inflows (their guessed flows, at their sources' Tt and
        composition) to the target on the chosen side of stoichiometric, by
        complete combustion. An explicit ``m_dot_guess`` wins."""
        comb = node.combustor
        fuel_path = node._fuel_path_elements(self.network)
        inflow = [
            e for e in self.network.elements.values() if e.to_node == comb and e.id not in fuel_path
        ]
        m_air_guessed = sum(abs(mdot_guess.get(e.id, ref["m_dot"])) for e in inflow)
        m_air = m_air_guessed or ref["m_dot"]
        sources = self._fuel_oxidizer_sources(node)
        # Known before the solve: no source can supply oxidizer -- none
        # upstream, or only fixed flows of zero. Every flow could then vanish
        # and the Tt row be met by an undetermined state (a "root" Linux CI
        # returned as a success), so it is refused here.
        if all(
            isinstance(n, MassFlowBoundary)
            and not isinstance(n, FuelBoundary)
            and abs(float(n.m_dot)) <= _FUEL_MIN_OXIDIZER_FLOW
            for n in sources
        ):
            raise ValueError(
                f"FuelBoundary '{node.id}': the combustor {comb!r} gets no oxidizer flow "
                "(no source upstream of it supplies any): the fuel has nothing to burn with."
            )
        temps = [float(n.Tt) for n in sources if getattr(n, "Tt", None) is not None]
        T_ox = float(np.mean(temps)) if temps else float(ref["T"])
        Y_ox = next((n.Y for n in sources if getattr(n, "Y", None) is not None), ref["Y"])
        Y_f = (
            node.Y
            if node.Y is not None
            else list(cb.species.to_mass(cb.species.pure_species("CH4")))
        )
        # Kept for the diagnosis of a failed solve: at a point where the mass
        # rows do not close (a premix plenum), the at-point oxidizer is not
        # the network's.
        node._start_oxidizer = (
            [(m_air_guessed, T_ox, list(Y_ox))]
            if sources and m_air_guessed > _FUEL_MIN_OXIDIZER_FLOW
            else []
        )
        if node.m_dot_guess is not None:
            return float(node.m_dot_guess)
        try:
            m = self._fuel_flow_for(node, m_air, T_ox, list(Y_ox), list(Y_f))
            node._guess_note = ""
            # The unknown is ln(m_dot): a start at ~0 leaves Newton no slope.
            return max(m, 1e-9 * m_air)
        except (ValueError, RuntimeError) as exc:
            # No root on that branch at the guessed oxidizer state: start lean
            # and small; the solve then reports whether the target is
            # reachable, with this reason.
            method = getattr(self.network.nodes.get(comb), "method", "complete")
            node._guess_note = (
                f"No {node.branch} root at the start ({exc}; oxidizer {T_ox:.1f} K, "
                f"{method} combustion)."
            )
            return 0.01 * m_air

    def _fuel_oxidizer_sources(self, node: "FuelBoundary") -> list[NetworkNode]:
        """The oxidizer's sources: every flow boundary upstream of the
        combustor that is not this fuel (a premix plenum's air, or direct)."""
        seen: set[str] = set()
        stack = [node.combustor]
        sources: list[NetworkNode] = []
        while stack:
            nid = stack.pop()
            for e in self.network.get_upstream_elements(nid):
                up = e.from_node
                if up in seen or up == node.id:
                    continue
                seen.add(up)
                upn = self.network.nodes[up]
                if isinstance(upn, (PressureBoundary, MassFlowBoundary)):
                    sources.append(upn)
                else:
                    stack.append(up)
        return sources

    def _fuel_flow_for(
        self,
        node: "FuelBoundary",
        m_air: float,
        T_ox: float,
        Y_ox: list[float],
        Y_f: list[float],
    ) -> float:
        """The fuel flow that burns ``m_air`` of oxidizer to the target on the
        node's branch by the combustor's own model (mix, then
        complete_combustion with its smooth phi blending, or
        combustion_equilibrium for an equilibrium combustor): bisection, exact
        to 1e-6 K; lean on [0, stoichiometric], rich on [stoichiometric, 10 x].
        It finds a target a fraction of a kelvin above the no-fuel Tt, which
        a coarse search returns as no fuel.

        Rich of stoichiometric the two models part far: at phi 2.2, 1500 K by
        equilibrium is 2020 K complete, so a complete-combustion start for an
        equilibrium combustor began at about twice the root's fuel flow and
        reached the rich root only through the LM fallback, or not at all."""
        comb = self.network.nodes.get(node.combustor)
        equilibrium = getattr(comb, "method", "complete") == "equilibrium"
        ox = self._as_stream(m_air, T_ox, Y_ox)
        X_ox = list(ox.X)
        own = float(cb.oxygen_required_per_kg_mixture(X_ox))
        o2 = Y_ox[list(cb.species.names).index("O2")] - own
        need = float(cb.oxygen_required_per_kg_mixture(list(cb.mass_to_mole(Y_f))))
        if o2 <= 0.0 or need <= 0.0:
            raise ValueError("no oxygen to spare, or nothing to burn")
        m_st = m_air * o2 / need

        def burnt(m: float) -> float:
            if m <= 0.0:
                if own > 0.0:
                    return float(cb.complete_combustion(ox.T, X_ox, ox.P).T)
                return float(ox.T)
            mixed = cb.mix([ox, self._as_stream(m, float(node.Tt), Y_f)])
            if equilibrium:
                return float(cb.combustion_equilibrium(mixed.T, list(mixed.X), mixed.P).state.T)
            return float(cb.complete_combustion(mixed.T, list(mixed.X), mixed.P).T)

        lo, hi = (0.0, m_st) if node.branch == "lean" else (m_st, 10.0 * m_st)
        f_lo = burnt(lo) - node.Tt_target
        f_hi = burnt(hi) - node.Tt_target
        if f_lo * f_hi > 0.0:
            raise ValueError(
                f"target outside [{burnt(lo):.1f}, {burnt(hi):.1f}] K on the {node.branch} side"
            )
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            f_mid = burnt(mid) - node.Tt_target
            if abs(f_mid) < 1e-6 or hi - lo <= 1e-15 * m_st:
                return mid
            if (f_mid > 0.0) == (f_lo > 0.0):
                lo, f_lo = mid, f_mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def _fuel_fixed_oxidizer(
        self, node: "FuelBoundary"
    ) -> list[tuple[float, float, list[float]]] | None:
        """The oxidizer a fuel controller's combustor burns, as fixed-flow
        streams [(m_dot, Tt, Y)], when the inputs alone fix it (#1114); None
        otherwise.

        Fixed means: everything upstream of the combustor is an in-tree of
        plain OrificeElement / LosslessConnectionElement edges (no wall on
        any) through PlenumNodes with no energy boundary, whose leaves are
        this fuel and fixed-flow MassFlowBoundaries (m_dot above
        _FUEL_MIN_OXIDIZER_FLOW, one edge each); and the combustor, with no
        energy boundary, has exactly one other edge, its outflow. Mass
        balance from the leaves then fixes every flow in the tree, all
        towards the combustor whatever the pressures, and those edges and
        nodes conserve total enthalpy and composition. So at ANY root the
        combustor burns the fuel with exactly these streams, and with the
        single outflow nothing else enters it. A PressureBoundary leaf is
        not fixed: its flow can reverse and carry other gas in."""
        net = self.network
        comb = net.nodes.get(node.combustor)
        if type(comb) is not CombustorNode or comb.energy_boundaries:
            return None
        walled: set[str] = set()
        for w in net.walls.values():
            walled.update(e for e in (w.element_a, w.element_b) if e)
        touching: dict[str, list[Any]] = {}
        for e in net.elements.values():
            for nid in {*e.all_source_nodes(), *e.all_sink_nodes()}:
                touching.setdefault(nid, []).append(e)

        def plain(e: Any) -> bool:
            return (
                type(e) in (OrificeElement, LosslessConnectionElement)
                and e.id not in walled
                and not _is_closed_restriction(e)
            )

        outflow = [e for e in touching.get(comb.id, []) if e.to_node != comb.id]
        if len(outflow) != 1 or comb.id in outflow[0].all_sink_nodes():
            return None
        stack = [e for e in touching.get(comb.id, []) if e.to_node == comb.id]
        seen: set[str] = set()
        fuel_seen = False
        streams: list[tuple[float, float, list[float]]] = []
        while stack:
            e = stack.pop()
            up = e.from_node
            if not plain(e) or up in seen or up == comb.id:
                return None
            seen.add(up)
            n = net.nodes[up]
            edges = touching.get(up, [])
            if n is node or type(n) is MassFlowBoundary:
                if len(edges) != 1 or edges[0] is not e:
                    return None
                if n is node:
                    fuel_seen = True
                    continue
                if n.Tt is None or n.Y is None or float(n.m_dot) <= _FUEL_MIN_OXIDIZER_FLOW:
                    return None
                streams.append((float(n.m_dot), float(n.Tt), list(n.Y)))
            elif type(n) is PlenumNode and not getattr(n, "energy_boundaries", None):
                for other in edges:
                    if other is e:
                        continue
                    if other.to_node != up:
                        return None
                    stack.append(other)
            else:
                return None
        return streams if fuel_seen and streams else None

    def _fuel_proven_unreachable(self, node: "FuelBoundary") -> tuple[str, float]:
        """Why no root of the solve can meet this fuel controller's target on
        its branch, proven from the inputs alone (#1114), or "" when that
        cannot be proven (the solve then decides, as before); and the Tt with
        no injected fuel.

        Only with the oxidizer fixed (``_fuel_fixed_oxidizer``). The
        combustor's Tt is then a function T(r) of the fuel to oxidizer ratio
        alone (ideal gas: mixing and complete combustion do not depend on P),
        and, complete combustion being exothermic and its heat release per kg
        of fuel above the fuel's own sensible heat up to the flame:

        - lean, T(r) rises from T(0), the Tt with no injected fuel, to
          stoichiometric, so a target below T(0) has no lean root (the
          partial equilibrium is complete combustion when lean);
        - T(r) for every r is at most T_hi, complete combustion WITHOUT the
          phi = 0 / phi = 1 smoothing at stoichiometric (the smoothing burns
          less, rich fuel only dilutes), so a target above T_hi has no root;
        - rich, complete combustion is at least the reactants' mixing T, so
          a target below min(oxidizer Tt, fuel Tt) has no rich root.

        Not proven, and left to the solve: a rich target below T(0) above
        that minimum (fuel as a coolant has a root at large phi: phi ~ 25 for
        798 K with air at 800 K and C3H8 at 300 K), and a rich target for an
        equilibrium combustor, whose reforming is P-dependent and at very rich
        phi runs above complete combustion. A root on the other branch is
        refused after the solve anyway. A margin of _FUEL_PROOF_MARGIN_K keeps
        targets at the window's ends to the solve."""
        streams = self._fuel_fixed_oxidizer(node)
        if streams is None:
            return "", math.nan
        try:
            Y_f = list(node.Y if node.Y is not None else cb.mole_to_mass(cb.species.dry_air()))
            ox = cb.mix([self._as_stream(m, T, Y) for m, T, Y in streams])
            X_ox = list(ox.X)
            Y_ox = list(cb.mole_to_mass(X_ox))
            o2 = Y_ox[list(cb.species.names).index("O2")] - float(
                cb.oxygen_required_per_kg_mixture(X_ox)
            )
            need = float(cb.oxygen_required_per_kg_mixture(list(cb.mass_to_mole(Y_f))))
            if o2 <= 0.0 or need <= 0.0:
                return "", math.nan
            m_st = float(ox.mdot) * o2 / need
            mixed = cb.mix([ox, self._as_stream(m_st, float(node.Tt), Y_f)])
            T_hi = float(
                cb.complete_combustion(
                    mixed.T, list(mixed.X), mixed.P, smooth_phi0=False, smooth_phi1=False
                ).T
            )
            T_without = self._fuel_T_without(streams)
        except (ValueError, RuntimeError):
            return "", math.nan
        complete = getattr(self.network.nodes[node.combustor], "method", "complete") == "complete"
        target = float(node.Tt_target)
        if node.branch == "lean":
            proven = (
                target > T_hi + _FUEL_PROOF_MARGIN_K or target < T_without - _FUEL_PROOF_MARGIN_K
            )
        else:
            proven = complete and (
                target > T_hi + _FUEL_PROOF_MARGIN_K
                or target < min(float(ox.T), float(node.Tt)) - _FUEL_PROOF_MARGIN_K
            )
        if not proven:
            return "", math.nan
        return self._fuel_window_reason(node, streams, Y_f, T_without, ""), T_without

    def _fuel_refused_before_solve(self) -> dict[str, Any] | None:
        """A failed result, without a Newton search, when a fuel controller's
        target is proven unreachable from the inputs (#1114,
        ``_fuel_proven_unreachable``); None otherwise. The reason is the one
        a failed solve reports."""
        fuels = [n for n in self.network.nodes.values() if isinstance(n, FuelBoundary)]
        if not fuels:
            return None
        self.network.resolve_all_topology()
        proofs = {f.id: self._fuel_proven_unreachable(f) for f in fuels}
        reasons = [why for why, _ in proofs.values() if why]
        if not reasons:
            return None
        self.network.validate()
        self._caught_exceptions = []
        self._caught_count = 0
        self._init_overrides = {}
        x0 = self._build_x0()
        sol: dict[str, Any] = {
            "__success__": False,
            "__message__": (
                "Refused before the solve (the combustor's inflows are fixed by the "
                f"inputs): {' '.join(reasons)}"
            ),
            "__iterations__": 0,
            "__converged__": False,
            "__consistent__": None,
            "__inconsistent_elements__": [],
            "__P_nonpositive__": [],
            "__T_clamped__": [],
            "__heat_faded__": [],
            "__supersonic_chambers__": [],
            "__choked_bounded_junctions__": [],
            "__outcome__": SolveOutcome.NOT_CONVERGED,
            "__final_norm__": float("inf"),
            "__worst_residuals__": [],
            "__unknown_names__": list(self.unknown_names),
            "__x_solution__": list(x0),
        }
        for f in fuels:
            sol[f"{f.id}.m_dot"] = float("nan")
            sol[f"{f.id}.phi"] = float("nan")
            sol[f"{f.id}.Tt_target"] = float(f.Tt_target)
            sol[f"{f.id}.Tt_without_fuel"] = proofs[f.id][1]
            sol[f"{f.id}.rich_branch"] = 1.0 if f.branch == "rich" else 0.0
        return self._with_caught(sol)

    def _guess_start(
        self,
    ) -> tuple[
        dict[str, Any], dict[str, float], dict[str, float], dict[str, float], dict[str, float]
    ]:
        """The cold start's estimates for the network ``self.network`` holds:
        (reference state, node pressures, element mass flows, FuelBoundary
        starts, tee branch flows). ``_build_x0`` runs it once per flow
        component (#899)."""
        ref = self._infer_reference_state()
        p_guess = self._propagate_pressure_guess(ref)
        mdot_guess = self._propagate_mdot_guess(ref)
        if self.resistance_seed:
            seeded = self._resistance_network_guess(ref)
            if seeded is not None:
                p_guess.update(seeded[0])
                mdot_guess.update(seeded[1])

        # A FuelBoundary's start (#892), from the other guessed inflows of its
        # combustor; the element that carries the fuel starts at that flow too,
        # or Newton begins with a mass balance decades out.
        fuel_start: dict[str, float] = {}
        for nid, node in self.network.nodes.items():
            if isinstance(node, FuelBoundary):
                g = self._fuel_guess(node, mdot_guess, ref)
                fuel_start[nid] = g
                node.m_dot = g
                for eid in node._fuel_path_elements(self.network):
                    mdot_guess[eid] = g

        # Bernoulli-based initial guess for TeeJunctionElement unknowns.
        # The generic mdot_guess propagator seeds from ref["m_dot"]=0.1 (arbitrary
        # when no MassFlowBoundary exists), but tee flows scale as A*sqrt(2*rho*dP)
        # which can be O(1) kg/s for typical pressure differences.  A factor-of-10
        # error in the starting mass flow prevents Newton from converging.
        from .components import TeeJunctionElement as _TeeJE

        _ref_X_tee = list(cb.mass_to_mole(ref["Y"]))
        _rho_tee = float(cb.density(ref["T"], ref["P"], _ref_X_tee))
        _a_tee = float(cb.speed_of_sound(ref["T"], _ref_X_tee))
        _gamma_tee = float(cb.isentropic_expansion_coefficient(ref["T"], _ref_X_tee))
        _mach_cap_target = 0.3
        _term_tee = 1.0 + 0.5 * (_gamma_tee - 1.0) * _mach_cap_target**2
        _exp_tee = -(_gamma_tee + 1.0) / (2.0 * (_gamma_tee - 1.0))
        _mdot_ratio_tee = _mach_cap_target * _term_tee**_exp_tee

        def _node_pt(n: str) -> float:
            node = self.network.nodes[n]
            if isinstance(node, PressureBoundary):
                return float(node.Pt)
            return p_guess.get(n, ref["P"])

        _tee_branch_guess: dict[str, float] = {}
        for _eid, _elem in self.network.elements.items():
            if not isinstance(_elem, _TeeJE):
                continue
            _src = _elem.all_source_nodes()
            _snk = _elem.all_sink_nodes()
            _pt_src = max(_node_pt(n) for n in _src)
            _pt_snk = min(_node_pt(n) for n in _snk)
            _dP_total = max(_pt_src - _pt_snk, 1.0)
            # F_C may be None before topology is resolved (e.g. when _build_x0
            # is called directly before NetworkSolver.solve resolves topology).
            # Skip the Bernoulli estimate and fall back to the propagated guess.
            if not _elem.F_C:
                continue
            # Mach cap: the tee duct cannot carry more than M=0.3 at reference
            # conditions.  Without this cap the Bernoulli estimate can overestimate
            # the tee flow by 3x when downstream channel losses are included in the
            # boundary-to-boundary dP, causing Newton to start at a wildly
            # inconsistent mass balance.
            _mdot_cap_tee = _rho_tee * _a_tee * _elem.F_C * _mdot_ratio_tee
            _bernoulli_com = min(
                _elem.F_C * math.sqrt(2.0 * _rho_tee * _dP_total),
                _mdot_cap_tee,
            )
            # Prefer a propagated mdot (e.g. from MassFlowBoundary) over the
            # Bernoulli estimate when the former is larger -- the propagated value
            # is exact when a mass-flow BC seeds the network, while the Bernoulli
            # estimate can be a factor of ~6 too small when only one pressure
            # boundary exists (dP from p_guess is only 0.1% of ref_P).
            mdot_guess[_eid] = max(_bernoulli_com, mdot_guess.get(_eid, 0.0))
            # Branch guess: Bernoulli estimate for the branch arm pressure drop.
            _pt_br = _node_pt(_elem.branch_node)
            if _elem.tee_type == "merging":
                _dP_br = max(_pt_br - _pt_snk, 1.0)
            else:
                _dP_br = max(_pt_src - _pt_br, 1.0)
            _bernoulli_br = min(
                _elem.F_C * math.sqrt(2.0 * _rho_tee * _dP_br),
                mdot_guess[_eid] * 0.75,
            )
            _tee_branch_guess[_eid] = max(_bernoulli_br, mdot_guess[_eid] * 0.25)
        return ref, p_guess, mdot_guess, fuel_start, _tee_branch_guess

    def _build_x0(self) -> np.ndarray:
        """
        Constructs the initial guess vector by gathering all unknowns.

        Pressure unknowns are initialised via topological propagation
        from boundary nodes with an estimated per-element drop.
        Temperature and composition are averaged across boundaries.
        Mass-flow is propagated from ``MassFlowBoundary`` nodes through
        the topology, with conservative capping for compressible elements.

        Returns a flat numpy array.
        """
        self.unknown_names = []
        self._unknown_indices = {}
        x0_list = []

        # --- Infer sensible defaults from boundary nodes ----------------
        # Per flow component (#899): a part no element connects to the rest
        # -- a coolant supply -- gets its own reference state. One global
        # reference mixed them: a disconnected 358 K coolant duct beside 19
        # near-sonic 808 K gas paths pulled the reference to 583 K, every gas
        # path started 18% high (density), and the solve landed nine chambers
        # on the supersonic root. A network of one component takes the
        # global path unchanged.
        comps = self._flow_components()
        ref_of: dict[str, dict[str, Any]] = {}
        if len(comps) <= 1:
            ref, p_guess, mdot_guess, fuel_start, _tee_branch_guess = self._guess_start()
        else:
            p_guess, mdot_guess, fuel_start, _tee_branch_guess = {}, {}, {}, {}
            for node_ids, elem_ids in comps:
                with self._component_view(node_ids, elem_ids):
                    r, pg, mg, fs, tb = self._guess_start()
                p_guess.update(pg)
                mdot_guess.update(mg)
                fuel_start.update(fs)
                _tee_branch_guess.update(tb)
                for i in (*node_ids, *elem_ids):
                    ref_of[i] = r
            ref = ref_of[comps[0][0][0]]

        # EnergyBoundary(Tt_target) heat unknowns (#792): {name: (node, eb)}
        target_q = {
            eb.unknown_name: (nid, eb)
            for nid, n in self.network.nodes.items()
            for eb in target_energy_boundaries(n)
        }

        # Iterate through interior nodes and gather unknowns
        for node_id, node in self.network.nodes.items():
            if not isinstance(node, PressureBoundary):
                unknowns = node.unknowns()
                if unknowns:
                    R = ref_of.get(node_id, ref)
                    start_idx = len(x0_list)
                    for unk in unknowns:
                        self.unknown_names.append(unk)
                        guess_dict = getattr(node, "initial_guess", {})
                        if unk in guess_dict:
                            x0_list.append(guess_dict[unk])
                        elif unk in self._init_overrides:
                            x0_list.append(self._init_overrides[unk])
                        elif unk.endswith(".P") or unk.endswith(".Pt"):
                            x0_list.append(p_guess.get(node_id, R["P"]))
                        elif unk.endswith(".T") or unk.endswith(".Tt"):
                            x0_list.append(R["T"])
                        elif ".Y[" in unk:
                            idx = int(unk.split(".Y[")[1].replace("]", ""))
                            x0_list.append(R["Y"][idx])
                        elif unk in target_q:
                            # Refined below, from the propagated start.
                            x0_list.append(0.0)
                        elif isinstance(node, FuelBoundary) and unk == node.m_dot_unknown:
                            x0_list.append(math.log(max(fuel_start[node_id], 1e-12)))
                        else:
                            x0_list.append(R["m_dot"])
                    self._unknown_indices[node_id] = list(range(start_idx, len(x0_list)))

        # Iterate through elements and gather unknowns (m_dot)
        for elem_id, element in self.network.elements.items():
            unknowns = element.unknowns()
            if unknowns:
                R = ref_of.get(elem_id, ref)
                start_idx = len(x0_list)
                for unk in unknowns:
                    self.unknown_names.append(unk)
                    guess_dict = getattr(element, "initial_guess", {})
                    if unk in guess_dict:
                        x0_list.append(guess_dict[unk])
                    elif unk in self._init_overrides:
                        x0_list.append(self._init_overrides[unk])
                    elif unk.endswith(".m_dot") or unk.endswith(".m_dot_com"):
                        x0_list.append(mdot_guess.get(elem_id, R["m_dot"]))
                    elif unk.endswith(".m_dot_branch"):
                        x0_list.append(
                            _tee_branch_guess.get(
                                elem_id, mdot_guess.get(elem_id, R["m_dot"]) * 0.5
                            )
                        )
                    elif unk.endswith(".P_jct") or unk.endswith(".Pt_jct"):
                        # MultiPortChamberBase: P_jct should be close to the
                        # ambient port static pressure. Use the max boundary Pt
                        # as a stand-in (P_jct >= Pt at every port since the
                        # impulse equation gives P_jct = Pt + q with q >= 0).
                        x0_list.append(R["P"])
                    else:
                        x0_list.append(1.0)
                self._unknown_indices[elem_id] = list(range(start_idx, len(x0_list)))

        # Build name to index mapping
        self._name_to_index = {name: i for i, name in enumerate(self.unknown_names)}

        # MultiPortChamberBase ports carry their throughflow on the outer
        # connecting elements' m_dot unknowns. Stash the resolved global
        # indices on each junction so flow_at_node / flow_jac_at_node can
        # report real port flows during state propagation: collector-port
        # MCNs need them for the Pt = P + 0.5*rho*v^2 closure and for
        # mass-weighted mixing of merging streams.
        from .components import MultiPortChamberBase as _MPCElem

        for element in self.network.elements.values():
            if isinstance(element, _MPCElem):
                element._port_outer_mdot_idx = {
                    outer_id: self._name_to_index[f"{outer_id}.m_dot"]
                    for outer_id in element._port_element_ids
                    if f"{outer_id}.m_dot" in self._name_to_index
                }

        # Pre-compute topological order for forward state propagation
        self._topological_order = self._compute_topological_order()

        # Set up flow-dependent wall EnergyBoundaries for thermal coupling
        # The wall EBs live on the nodes, so drop any a previous setup left:
        # each would keep the heat of that attempt's last iterate (a retry,
        # or a second solver on the same network, used to add a second EB
        # beside the stale one -- the heat counted twice).
        for node in self.network.nodes.values():
            ebs = getattr(node, "energy_boundaries", None)
            if ebs:
                ebs[:] = [eb for eb in ebs if eb.id != f"_wall_{node.id}"]
        self._wall_ebs: dict[str, EnergyBoundary] = {}
        if self.network.thermal_coupling_enabled and self.network.walls:
            self._setup_wall_energy_boundaries()

        x0 = np.array(x0_list, dtype=float)
        self._target_q_columns = {
            nid: [(eb, self._name_to_index[name])] for name, (nid, eb) in target_q.items()
        }
        self._seed_chambers_past_sonic(x0)
        if target_q:
            self._seed_target_heat(x0, target_q)
        return x0

    def _seed_chambers_past_sonic(self, x0: np.ndarray) -> None:
        """Start a compressible momentum chamber whose start lies past its
        sonic point on its closure's root (#788).

        The closure is monotone in P with ONE root per Pt, but its slope is
        zero at the sonic static pressure P*. A start with P = Pt that the
        flux puts below P* (a choked channel into a junction whose
        downstream fixes a low Pt: a bore equal to the pipe) made the first
        Newton step cross that fold and overshoot to a negative P, and the
        cold solve stalled. Only such a start moves; a chamber with a user
        initial guess for P, or a start on the subsonic side, keeps its x0.
        """
        chambers = [
            (nid, node)
            for nid, node in self.network.nodes.items()
            if isinstance(node, MomentumChamberNode)
            and node.regime != "incompressible"
            and node.area
            and f"{nid}.P" in self._name_to_index
            and f"{nid}.P" not in getattr(node, "initial_guess", {})
        ]
        if not chambers:
            return
        # The probe propagation must leave no trace: a network with a
        # back-edge reads the previous pass's states (lagged reads), so the
        # first residual would otherwise depend on this probe.
        _kept = ("_derived_states", "_prev_derived_states", "_film_inf_states")
        _kept += ("_prev_film_inf_states", "_rev_junctions")
        saved = {k: self.__dict__[k] for k in _kept if k in self.__dict__}
        try:
            self._propagate_states(x0, with_relay=False)
            self._seed_chamber_roots(x0, chambers)
        except Exception:  # noqa: BLE001 -- a start estimate never breaks a solve
            pass
        finally:
            for k in _kept:
                if k in saved:
                    self.__dict__[k] = saved[k]
                else:
                    self.__dict__.pop(k, None)

    def _seed_chamber_roots(
        self, x0: np.ndarray, chambers: list[tuple[str, MomentumChamberNode]]
    ) -> None:
        """``_seed_chambers_past_sonic`` on the propagated start."""
        from scipy.optimize import brentq

        for nid, node in chambers:
            try:
                state = self._get_node_state(node, x0)
                c = node._closure(state)
                Pt = float(state.Pt)
                if not c.bounded or not (math.isfinite(Pt) and Pt > 0.0):
                    continue
                P_star = float(c.P_sonic)

                def r(
                    lnP: float,
                    state: NetworkMixtureState = state,
                    node: MomentumChamberNode = node,
                ) -> float:
                    state.P = math.exp(lnP)
                    return float(node._closure(state).residual)

                lo, hi = math.log(1e-4 * P_star), math.log(max(Pt, P_star))
                if not (r(lo) > 0.0 > r(hi)):
                    continue
                x0[self._name_to_index[f"{nid}.P"]] = math.exp(brentq(r, lo, hi, xtol=1e-10))
            except (ValueError, RuntimeError):
                continue

    def _seed_target_heat(
        self, x0: np.ndarray, target_q: dict[str, tuple[str, EnergyBoundary]]
    ) -> None:
        """Start each Tt_target heat unknown (#792) at the heat that takes the
        node's propagated start temperature (with no target heat) to its
        target at the start flow: Q0 = m (h(Tt_target) - h(T0)). It also sets
        the column scale (D_x = |x0|), so it is floored at the heat of 1 K of
        that flow. A user initial_guess for the unknown is kept."""
        try:
            self._propagate_states(x0, with_relay=False)
        except Exception:  # noqa: BLE001 -- a start estimate never breaks a solve
            return
        for name, (nid, eb) in target_q.items():
            if name in getattr(self.network.nodes[nid], "initial_guess", {}):
                continue
            derived = self._derived_states.get(nid)
            if derived is None:
                continue
            T0, Y0, _ = derived
            m = abs(sum(m for m, _, _ in self._node_inflows(nid, x0)))
            m = max(m, 1e-3)
            try:
                X = list(cb.mass_to_mole(list(Y0)))
                cp = float(cb.cp_mass(float(T0), X))
                q = m * (float(cb.h_mass(float(eb.Tt_target), X)) - float(cb.h_mass(float(T0), X)))
            except (ValueError, RuntimeError):
                continue
            floor = m * cp * 1.0
            if not math.isfinite(q):
                continue
            x0[self._name_to_index[name]] = q if abs(q) >= floor else math.copysign(floor, q)

    def _compute_topological_order(self) -> list[str]:
        """Simple topological sort starting from boundary nodes."""
        order = []
        visited = set()
        queue = []
        # Boundaries are starting points
        for nid, node in self.network.nodes.items():
            if node.is_flow_boundary:
                queue.append(nid)
                visited.add(nid)

        while queue:
            nid = queue.pop(0)
            order.append(nid)
            for elem in self.network.get_downstream_elements(nid):
                for to_id in elem.all_sink_nodes():
                    if to_id not in visited:
                        up_elems = self.network.get_upstream_elements(to_id)
                        # Ready when all source nodes of all upstream elements are visited
                        if all(src in visited for e in up_elems for src in e.all_source_nodes()):
                            visited.add(to_id)
                            queue.append(to_id)

        # Cleanup: append any unvisited (cycles or islands)
        for nid in self.network.nodes:
            if nid not in visited:
                order.append(nid)
        return order

    def _setup_wall_energy_boundaries(self) -> None:
        """Create flow-dependent EnergyBoundary objects for wall-coupled nodes.

        For each node that has upstream elements participating in a WallConnection,
        a dedicated EnergyBoundary is created and attached to the node.  Its Q is
        updated every Newton iteration inside ``_propagate_states``.
        """

        def _endpoint_to_node_ids(endpoint_id: str) -> list[str] | None:
            # Either end of an element can receive the heat: the one its
            # flow goes into, which reverses with the flow (#481).
            if endpoint_id in self.network.elements:
                e = self.network.elements[endpoint_id]
                return [e.to_node, e.from_node]
            if endpoint_id in self.network.nodes:
                return [endpoint_id]
            return None

        # Build a mapping: node_id -> set of wall IDs that affect it
        affected_nodes: dict[str, set[str]] = {}
        for wall in self.network.walls.values():
            # A wall with a side B of no stream (#970) heats side A's alone;
            # an adiabatic one heats nothing.
            ends = [_endpoint_to_node_ids(e) for e in wall.heated_endpoints]

            if not ends or any(n is None for n in ends):
                continue

            for n in (n for ns in ends for n in ns):
                affected_nodes.setdefault(n, set()).add(wall.id)

        for nid in affected_nodes:
            node = self.network.nodes[nid]
            if not hasattr(node, "energy_boundaries"):
                continue
            eb = EnergyBoundary(id=f"_wall_{nid}", Q=0.0)
            self._wall_ebs[nid] = eb
            node.add_energy_boundary(eb)

    @staticmethod
    def _wall_gas(
        wall: Any, state_a: NetworkMixtureState, state_b: NetworkMixtureState
    ) -> dict[str, tuple[float, float, list[float]]]:
        """The streams' (Tt, P, Y) a wall's radiation reads (#793); nothing
        for a wall without radiation, whose coupling call stays as it was."""
        if wall.radiation_a is None and wall.radiation_b is None:
            return {}
        return {
            "gas_a": (float(state_a.Tt), float(state_a.P), list(state_a.Y)),
            "gas_b": (float(state_b.Tt), float(state_b.P), list(state_b.Y)),
        }

    def _wall_capacity(
        self, eid: str, state: NetworkMixtureState, T_other: float
    ) -> tuple[float, float, float, float, float, Callable[[], np.ndarray] | None]:
        """(C, dC/dm, dC/dT_own, dC/dT_other, flow sign, dC/dY_own) of one wall
        side; dC/dY_own is a callable, evaluated only when the relay needs it.

        C = |m| cp [W/K] of the stream feeding it, cp the secant over the span
        from its own inlet T to the other side's (cb.cp_secant_mass_and_jacobian):
        the variable-property capacity, so C_min dT is the exact maximum heat
        and the weaker stream at most reaches the other's inlet T. cp at its
        own inlet let a 1600 K stream's high cp overstate what it gives up and
        the parallel outlets cross by 43 K. dC/dm w.r.t. the element's own flow
        unknown; the T partials w.r.t. the two feed nodes' T. A node side is a
        well-mixed reservoir: C = inf (#627)."""
        if eid in self.network.nodes:
            return math.inf, 0.0, 0.0, 0.0, 1.0, None
        m = float(state.m_dot)
        X = list(cb.mass_to_mole(list(state.Y)))
        cp, dcp_own, dcp_other = cb.cp_secant_mass_and_jacobian(float(state.Tt), float(T_other), X)
        sign = -1.0 if m < 0.0 else 1.0
        a = abs(m)
        # Its composition moves cp too (#717): on demand, for the relay.
        Tt = float(state.Tt)

        def dC_dY() -> np.ndarray:
            return np.asarray(cb.cp_secant_mass_dY(Tt, float(T_other), X)) * a

        return a * cp, sign * cp, a * dcp_own, a * dcp_other, sign, dC_dY

    @staticmethod
    def _wall_side_dQ_dY(
        wall_result: Any, ch: Any, dC_dY: Callable[[], np.ndarray] | None, side_a: bool
    ) -> np.ndarray | None:
        """dQ/dY of one wall side's feed composition (#717): through its h and
        T_aw (a surface that gives no composition partials contributes none)
        and its capacity rate. None when nothing moves."""
        dQ_dh = wall_result.dQ_dh_a if side_a else wall_result.dQ_dh_b
        dQ_dTaw = wall_result.dQ_dT_aw_a if side_a else wall_result.dQ_dT_aw_b
        dQ_dC = wall_result.dQ_dC_a if side_a else wall_result.dQ_dC_b
        out = None
        for coeff, get in (
            (dQ_dh, lambda: getattr(ch, "dh_dY", None)),
            (dQ_dTaw, lambda: getattr(ch, "dT_aw_dY", None)),
            (dQ_dC, dC_dY),
        ):
            if get is None or coeff == 0.0:
                continue
            grad = get()
            if grad is None or len(grad) == 0:
                continue
            term = coeff * np.asarray(grad, dtype=float)
            out = term if out is None else out + term
        return out

    def _evaluate_walls_for_node(
        self,
        nid: str,
        up_elems: list,
        x: np.ndarray,
    ) -> list[tuple]:
        """Evaluate all wall couplings affecting *nid* and return contributions.

        Returns a list of ``(wall, Q_inject, wall_result, elem_a_id, elem_b_id,
        is_side_a, ch_result_a, ch_result_b)`` tuples.  Also updates the wall
        ``EnergyBoundary`` Q on the node so that ``compute_derived_state`` picks
        it up.
        """
        wall_contributions: list[tuple] = []
        processed_walls: set[str] = set()

        for wall in self.network.walls.values():
            if wall.id in processed_walls:
                continue
            if wall.one_sided:
                # Side B has no stream (#970): adiabatic passes no heat; a
                # HeatLoss side takes side A's heat out of the network.
                if wall.heat_loss is not None:
                    wc = self._evaluate_heat_loss_wall(wall, nid, x)
                    if wc is not None:
                        wall_contributions.append(wc)
                        processed_walls.add(wall.id)
                continue

            # Endpoints may be elements or nodes
            elem_a_id = wall.element_a
            elem_b_id = wall.element_b

            endpoint_a_is_node = elem_a_id in self.network.nodes
            endpoint_b_is_node = elem_b_id in self.network.nodes

            obj_a = self.network.elements.get(elem_a_id) or self.network.nodes.get(elem_a_id)
            obj_b = self.network.elements.get(elem_b_id) or self.network.nodes.get(elem_b_id)

            if obj_a is None or obj_b is None:
                continue

            # The node each side's heat goes into, and the node feeding it,
            # by the sign of that side's flow (#481).
            feed_a, target_node_a = self._flow_ends(obj_a, x)
            feed_b, target_node_b = self._flow_ends(obj_b, x)

            # This wall does not affect the current node
            if nid not in (target_node_a, target_node_b):
                continue

            # Skip if both sides feed/are the same node (net Q = 0)
            if target_node_a == target_node_b:
                continue

            # Determine which side feeds the current node
            is_side_a = target_node_a == nid

            # Get states for both sides (previous-iteration fallback for back-edges)
            state_a = self._get_node_state_with_prev(self.network.nodes[feed_a], x)
            state_b = self._get_node_state_with_prev(self.network.nodes[feed_b], x)

            # Set mass flow rates from solver vector for elements
            if not endpoint_a_is_node:
                m_indices_a = self._unknown_indices.get(elem_a_id, [])
                if m_indices_a:
                    state_a.m_dot = x[m_indices_a[0]]
            if not endpoint_b_is_node:
                m_indices_b = self._unknown_indices.get(elem_b_id, [])
                if m_indices_b:
                    state_b.m_dot = x[m_indices_b[0]]

            # Call htc_and_T() on both objects, a film applied to its side (#386)
            ch_result_a = self._wall_side_htc(wall, obj_a, state_a, x)
            ch_result_b = self._wall_side_htc(wall, obj_b, state_b, x)

            # Only proceed if both elements have heat transfer surfaces
            if ch_result_a is not None and ch_result_b is not None:
                # Extract scalars from ChannelResult
                h_a = ch_result_a.h
                T_aw_a = ch_result_a.T_aw
                h_b = ch_result_b.h
                T_aw_b = ch_result_b.T_aw

                # The streams' exchange (epsilon-NTU): their capacity rates and
                # whether they run counter along the wall, from the flow signs
                # and the declared orientation (#627).
                C_a, dCa_dm, dCa_dTa, dCa_dTb, sign_a, dCa_dYa = self._wall_capacity(
                    elem_a_id, state_a, float(state_b.Tt)
                )
                C_b, dCb_dm, dCb_dTb, dCb_dTa, sign_b, dCb_dYb = self._wall_capacity(
                    elem_b_id, state_b, float(state_a.Tt)
                )
                wall_result = wall.compute_coupling(
                    h_a,
                    T_aw_a,
                    obj_a.surface.area,
                    h_b,
                    T_aw_b,
                    obj_b.surface.area,
                    C_a,
                    C_b,
                    wall.is_counterflow(sign_a, sign_b),
                    **self._wall_gas(wall, state_a, state_b),
                )

                # Convention: Q > 0 means heat flows A->B
                Q_inject = -wall_result.Q if is_side_a else wall_result.Q

                # Store ChannelResults for Jacobian relay
                wall_contributions.append(
                    (
                        wall,
                        Q_inject,
                        wall_result,
                        elem_a_id,
                        elem_b_id,
                        is_side_a,
                        ch_result_a,
                        ch_result_b,
                        (dCa_dm, dCa_dTa, dCa_dTb, dCb_dm, dCb_dTb, dCb_dTa, dCa_dYa, dCb_dYb),
                    )
                )
                processed_walls.add(wall.id)

        # Set accumulated wall Q on the node's wall EnergyBoundary
        if nid in self._wall_ebs:
            self._wall_ebs[nid].Q = sum(c[1] for c in wall_contributions)

        return wall_contributions

    def _one_sided_wall_inputs(
        self, wall: Any, x: np.ndarray, lagged: bool = True, nid: str | None = None
    ) -> tuple[Any, str, str, NetworkMixtureState, Any] | None:
        """(side A object, its feed node, its target node, its state, its
        surface result) of a wall whose side B has no stream (#970); None
        without a surface, or when side A's heat does not go into ``nid``
        (if given)."""
        obj_a = self.network.elements.get(wall.element_a) or self.network.nodes.get(wall.element_a)
        if obj_a is None:
            return None
        feed_a, target_a = self._flow_ends(obj_a, x)
        if nid is not None and nid != target_a:
            return None
        node = self.network.nodes[feed_a]
        state_a = (
            self._get_node_state_with_prev(node, x) if lagged else self._get_node_state(node, x)
        )
        if wall.element_a not in self.network.nodes:
            m_indices_a = self._unknown_indices.get(wall.element_a, [])
            if m_indices_a:
                state_a.m_dot = x[m_indices_a[0]]
        ch_a = self._wall_side_htc(wall, obj_a, state_a, x, lagged=lagged)
        if ch_a is None:
            return None
        return obj_a, feed_a, target_a, state_a, ch_a

    def _evaluate_heat_loss_wall(self, wall: Any, nid: str, x: np.ndarray) -> tuple | None:
        """A HeatLoss wall's contribution to node ``nid`` (#970), in the
        layout of ``_evaluate_walls_for_node``: side A's stream gives Q to
        the external reservoir (C++ ``wall_exchange_and_jacobian`` with
        C_b = inf), so the node takes -Q. Side B carries no unknowns: its
        entries are None and its capacity partials zero."""
        ins = self._one_sided_wall_inputs(wall, x, nid=nid)
        if ins is None:
            return None
        obj_a, _feed_a, _target_a, state_a, ch_a = ins
        hl = wall.heat_loss
        C_a, dCa_dm, dCa_dTa, _dCa_dText, _sign_a, dCa_dYa = self._wall_capacity(
            wall.element_a, state_a, float(hl.T_ext)
        )
        gas_a = (
            (float(state_a.Tt), float(state_a.P), list(state_a.Y))
            if wall.radiation_a is not None
            else None
        )
        wall_result = wall.heat_loss_coupling(
            ch_a.h, ch_a.T_aw, obj_a.surface.area, C_a, gas_a=gas_a
        )
        return (
            wall,
            -wall_result.Q,
            wall_result,
            wall.element_a,
            None,
            True,
            ch_a,
            None,
            (dCa_dm, dCa_dTa, 0.0, 0.0, 0.0, 0.0, dCa_dYa, None),
        )

    def _report_one_sided_wall(
        self, wid: str, wall: Any, x: np.ndarray, sol_dict: dict[str, Any]
    ) -> None:
        """Results of a wall whose side B has no stream (#970), the same
        evaluation as the residuals'. Adiabatic: Q = 0 and the wall at side
        A's face balance. HeatLoss: Q > 0 leaves side A's stream for T_ext,
        reported again as ``{wid}.Q_loss``, the heat leaving the network."""
        ins = self._one_sided_wall_inputs(wall, x, lagged=False)
        if ins is None:
            return
        obj_a, _feed_a, tgt_a, state_a, ch_a = ins
        gas_a = (
            (float(state_a.Tt), float(state_a.P), list(state_a.Y))
            if wall.radiation_a is not None
            else None
        )
        hl = wall.heat_loss
        if hl is None:
            res: Any = wall.adiabatic_coupling(ch_a.h, ch_a.T_aw, obj_a.surface.area, gas_a)
            C_a = math.inf
        else:
            C_a = self._wall_capacity(wall.element_a, state_a, float(hl.T_ext))[0]
            res = wall.heat_loss_coupling(ch_a.h, ch_a.T_aw, obj_a.surface.area, C_a, gas_a=gas_a)
        Q = float(res.Q)
        sol_dict[f"{wid}.Q"] = Q
        sol_dict[f"{wid}.Q_loss"] = Q
        sol_dict[f"{wid}.T_hot"] = float(res.T_hot)
        sol_dict[f"{wid}.UA"] = float(res.UA)
        sol_dict[f"{wid}.NTU"] = float(res.UA) / C_a if 0.0 < C_a < math.inf else 0.0
        sol_dict[f"{wid}.effectiveness"] = float(res.effectiveness)
        sol_dict[f"{wid}.counterflow"] = 0.0
        # Side A's heat on a node that cannot take it, as for two streams.
        q_to_boundary = 0.0
        if Q != 0.0 and tgt_a not in self._wall_ebs:
            if self.network.nodes[tgt_a].is_flow_boundary:
                q_to_boundary = -Q
                key = f"{tgt_a}.Q_wall_out"
            else:
                key = f"{tgt_a}.Q_withheld"
            sol_dict[key] = sol_dict.get(key, 0.0) - Q
        sol_dict[f"{wid}.Q_to_boundary"] = q_to_boundary
        sol_dict[f"{wid}.h_a"] = float(ch_a.h)
        sol_dict[f"{wid}.T_aw_a"] = float(ch_a.T_aw)
        if hl is not None:
            sol_dict[f"{wid}.h_b"] = float(hl.h)
            sol_dict[f"{wid}.T_aw_b"] = float(hl.T_ext)
        for key, val in (getattr(ch_a, "film", None) or {}).items():
            sol_dict[f"{wid}.{key}"] = float(val)
        for key, val in (getattr(res, "radiation", None) or {}).items():
            sol_dict[f"{wid}.{key}"] = float(val)
        sol_dict[f"{wid}.T_interface"] = list(wall._last_profile or [])
        sol_dict[f"{wid}.material_in_range"] = float(all(layer.in_range for layer in wall.layers))

    # Propagation passes per evaluation when a state is read ahead of the
    # order (wall back-edges, recirculation): enough for the contraction seen
    # on wall-coupled networks (x1000 per pass), bounded for safety.
    _MAX_PROPAGATION_PASSES = 50

    def _two_port_mdot_indices(self) -> dict[str, int]:
        """element id -> index of its own m_dot, for plain 2-port elements
        (whose flow may reverse); multi-port elements keep their declared
        directions, which their own models enforce."""
        cache = getattr(self, "_two_port_cache", None)
        if cache is not None and cache[0] is self.unknown_names:
            return cache[1]
        out: dict[str, int] = {}
        for eid, e in self.network.elements.items():
            ind = self._unknown_indices.get(eid, [])
            if (
                len(ind) == 1
                and self.unknown_names[ind[0]] == f"{eid}.m_dot"
                and e.all_source_nodes() == [e.from_node]
                and e.all_sink_nodes() == [e.to_node]
            ):
                out[eid] = ind[0]
        self._two_port_cache = (self.unknown_names, out)
        return out

    def _reversed_elements(self, x: np.ndarray) -> set[str]:
        return {eid for eid, i in self._two_port_mdot_indices().items() if x[i] < 0.0}

    def _reversed_junctions(self, x: np.ndarray) -> dict[str, tuple[list[float], list[int]]]:
        """Junctions with at least one port running against its declared
        direction: id -> (port flows, positive OUT of the junction; the
        global index of each port's outer m_dot).

        A junction feeds its port nodes in the declared direction otherwise
        (flow_at_node). At a reversed port that counts the stream leaving the
        node as a negative inflow, which cancels the real inflow from the
        reversed outer element: the node saw no flow and its Pt fell to P
        (0 Pa against a 711 Pa head, #517). Junctions listed here are fed by
        the sign of their port flows instead -- see _junction_streams. A
        junction with several inlets AND several outlets is always listed
        (flow_at_node credits each outlet with the whole throughput); other
        fully forward junctions are not, so their propagation is unchanged.
        """
        from .components import MultiPortChamberBase

        out: dict[str, tuple[list[float], list[int]]] = {}
        for eid, e in self.network.elements.items():
            if not isinstance(e, MultiPortChamberBase):
                continue
            idx_map = getattr(e, "_port_outer_mdot_idx", None) or {}
            idxs = [idx_map.get(oid) for oid in e._port_element_ids]
            if any(i is None for i in idxs):
                continue
            outer = [float(x[i]) for i in idxs]
            # Several inlets AND several outlets are always fed by the port
            # flows: flow_at_node gives each sink every source's own flow,
            # right for a merge (one outlet) but crediting EACH outlet with
            # the whole throughput here -- an outlet chamber's velocity head
            # at the total flow, and a 10 kPa jump in its closure as its own
            # port flow crossed zero (the stress ideal_nport family).
            several_each = len(e.inlet_nodes) >= 2 and len(e.outlet_nodes) >= 2
            if all(m >= 0.0 for m in outer) and not several_each:
                continue
            flows = [float(s) * m for s, m in zip(e._port_signs, outer, strict=True)]
            out[eid] = (flows, [int(i) for i in idxs])
        return out

    def _declared_direction_hints(self, x: np.ndarray) -> str:
        """Why a failed solve may have no root: a junction with a declared
        flow_direction (MultiPortChamberElement, ConstantKTeeElement) whose
        best iterate holds a port at zero or against that direction. Such a
        junction does not reverse a port, so an operating point that needs it
        has no root; the stress tee3 draws that stalled all solve with
        IdealMomentumTeeElement, with that port reversed.
        """
        from .components import MultiPortChamberBase

        hints = []
        for eid, e in self.network.elements.items():
            direction = getattr(e, "flow_direction", None)
            if not isinstance(e, MultiPortChamberBase) or direction not in ("branch", "merge"):
                continue
            idx_map = getattr(e, "_port_outer_mdot_idx", None) or {}
            idxs = [idx_map.get(oid) for oid in e._port_element_ids]
            if any(i is None for i in idxs):
                continue
            # Positive out of the junction, as in _reversed_junctions.
            flows = [float(sg) * float(x[i]) for sg, i in zip(e._port_signs, idxs, strict=True)]
            scale = max(abs(f) for f in flows)
            if not scale > 0.0:
                continue
            inlets = set(e.inlet_nodes)
            bad = [
                node
                for node, f in zip(e.port_nodes, flows, strict=True)
                if ((f > -1e-3 * scale) if node in inlets else (f < 1e-3 * scale))
            ]
            if bad:
                hints.append(
                    f"{type(e).__name__} '{eid}' (flow_direction='{direction}') holds "
                    f"port(s) {bad} at zero or against its declared direction at the "
                    "best iterate: a declared-direction junction does not reverse a "
                    "port, so this operating point may have no root for it. Check the "
                    "boundary pressures, change flow_direction, or use "
                    "IdealMomentumTeeElement, which reverses ports."
                )
        # Elements that can tell from their boundaries alone that the
        # operating point is outside their model (an ejector's outlet above
        # its dead head).
        for e in self.network.elements.values():
            hint = getattr(e, "no_root_hint", None)
            if callable(hint):
                text = hint(self.network)
                if text:
                    hints.append(text)
        return " ".join(hints)

    def _junction_streams(
        self, elem: Any, nid: str, flows: list[float], idxs: list[int]
    ) -> list[tuple[str, float, dict[int, float]]]:
        """Streams a reversed-port junction delivers INTO port node nid, by
        the sign of the port flows: (supplying port node, mass, d mass/d x).

        A port that supplies the junction (flow into it) gets nothing from
        it. A receiving port j gets its own flow m_j, drawn from the
        suppliers k (flow m_k < 0):

        - one supplier: one stream of m_j from it -- what a dividing junction
          does in its declared direction;
        - one receiver: one stream of |m_k| from each supplier -- what a
          joining junction does in its declared direction;
        - several of each: m_j |m_k| / sum|m| from each supplier, so every
          receiver carries the supplied mixture at its own flow.

        The derivatives are in the outer m_dot unknowns (m_i = sign_i x_i).
        """
        i = elem.port_nodes.index(nid)
        m_i = flows[i]
        if not m_i > 0.0:
            return []
        signs = elem._port_signs
        sup = [k for k, m in enumerate(flows) if m < 0.0]
        rec = [j for j, m in enumerate(flows) if m > 0.0]
        if not sup:
            return []
        if len(sup) == 1:
            k = sup[0]
            return [(elem.port_nodes[k], m_i, {idxs[i]: float(signs[i])})]
        if len(rec) == 1:
            return [(elem.port_nodes[k], -flows[k], {idxs[k]: -float(signs[k])}) for k in sup]
        total = -sum(flows[k] for k in sup)
        out = []
        for k in sup:
            a_k = -flows[k]
            jac: dict[int, float] = {idxs[i]: float(signs[i]) * a_k / total}
            for k2 in sup:
                # d(m_i a_k / S)/d a_k2, then d a_k2 / d x_k2 = -sign_k2.
                d_a = m_i * ((total - a_k) if k2 == k else -a_k) / (total * total)
                jac[idxs[k2]] = jac.get(idxs[k2], 0.0) - float(signs[k2]) * d_a
            out.append((elem.port_nodes[k], m_i * a_k / total, jac))
        return out

    def _flow_ends(self, obj: Any, x: np.ndarray) -> tuple[str, str]:
        """(node feeding the flow, node it delivers into) of a wall endpoint:
        a node is both; an element by the sign of its flow (#481)."""
        if obj.id in self.network.nodes:
            return obj.id, obj.id
        i = self._two_port_mdot_indices().get(obj.id)
        if i is not None and x[i] < 0.0:
            return obj.to_node, obj.from_node
        return obj.from_node, obj.to_node

    def _signed_order(
        self,
        reversed_ids: set[str],
        reversed_junctions: dict[str, tuple[list[float], list[int]]] | None = None,
    ) -> list[str]:
        """Propagation order with each reversed element pointing the way its
        flow goes; nodes left in a cycle follow in network order. A junction
        with a reversed port points from the ports that supply it to the
        ports it supplies."""
        nodes = list(self.network.nodes)
        preds: dict[str, set[str]] = {n: set() for n in nodes}
        succ: dict[str, set[str]] = {n: set() for n in nodes}
        reversed_junctions = reversed_junctions or {}
        for e in self.network.elements.values():
            if e.id in reversed_junctions:
                flows = reversed_junctions[e.id][0]
                srcs = [n for n, m in zip(e.port_nodes, flows, strict=True) if m < 0.0]
                sinks = [n for n, m in zip(e.port_nodes, flows, strict=True) if m > 0.0]
            elif e.id in reversed_ids:
                srcs, sinks = [e.to_node], [e.from_node]
            else:
                srcs, sinks = e.all_source_nodes(), e.all_sink_nodes()
            for a in srcs:
                for b in sinks:
                    if a != b:
                        preds[b].add(a)
                        succ[a].add(b)
        order: list[str] = []
        done: set[str] = set()
        queue = [n for n in nodes if not preds[n]]
        while queue:
            n = queue.pop(0)
            if n in done:
                continue
            done.add(n)
            order.append(n)
            for m in sorted(succ[n], key=nodes.index):
                if m not in done and preds[m] <= done:
                    queue.append(m)
        order += [n for n in nodes if n not in done]
        return order

    def _propagate_states(self, x: np.ndarray, with_relay: bool = True) -> dict:
        """Derive every node's T and Y (and the relay) at x.

        with_relay=False skips the sensitivity relay (dT/dx, dY/dx), which
        only the Jacobian uses and which dominates a residual evaluation on
        mixing and wall networks; the states are the same either way.

        Streams follow the SIGN of the flow: a reversed element feeds its
        from_node from its to_node (#481). When a pass reads a state ahead of
        its order (a wall back-edge, a recirculation loop) the pass is
        repeated from its own result until the temperatures settle, so the
        residual depends on x alone and not on the previous evaluation.

        The relay iterates with them: a node read ahead of its order takes
        that node's sensitivities from the previous pass, and the passes go on
        until those settle too. Each pass used to start the relay from empty,
        so a back-edge dropped its whole chain from the Jacobian -- 19% of
        d(row)/d(cooling flow) in a two-wall series network (#481).
        """
        # A FuelBoundary's flow is an unknown (#892): every reader of its m_dot
        # (propagation, mass rows, ledger) sees the current iterate.
        for node in self.network.nodes.values():
            if isinstance(node, FuelBoundary):
                node.m_dot = math.exp(float(x[self._name_to_index[node.m_dot_unknown]]))
        self._prev_derived_states = getattr(self, "_derived_states", {})
        self._prev_film_inf_states = self._film_inf_states
        reversed_ids = self._reversed_elements(x)
        self._rev_junctions = self._reversed_junctions(x)
        order = (
            self._signed_order(reversed_ids, self._rev_junctions)
            if reversed_ids or self._rev_junctions
            else self._topological_order
        )
        relay: dict = {}
        for _ in range(self._MAX_PROPAGATION_PASSES):
            self._lagged_read = False
            self._relay_borrowed: set[str] = set()
            self._derived_states = {}
            self._film_inf_states = {}
            self._T_no_gas: set[str] = set()
            prev_relay = relay
            relay = self._propagate_pass(x, order, reversed_ids, with_relay, prev_relay)
            if not self._lagged_read:
                break
            change = 0.0
            for n, (T, _, _) in self._derived_states.items():
                prev = self._prev_derived_states.get(n)
                T_prev = float(prev[0]) if prev is not None else 0.0
                change = max(change, abs(float(T) - T_prev) / max(abs(float(T)), 1.0))
            for k, (T, _) in self._film_inf_states.items():
                prev_f = self._prev_film_inf_states.get(k)
                T_prev = float(prev_f[0]) if prev_f is not None else 0.0
                change = max(change, abs(float(T) - T_prev) / max(abs(float(T)), 1.0))
            self._prev_derived_states = self._derived_states
            self._prev_film_inf_states = self._film_inf_states
            # The relay is at its fixed point when every sensitivity a pass
            # BORROWED from the previous one (read ahead of order, or a node
            # through its own walls) equals what this pass computed for it.
            # Junction networks read ahead nearly every evaluation through
            # chains that carry nothing: they settle in the first pass.
            if change < 1e-13 and (
                not with_relay or self._relay_settled(relay, prev_relay, self._relay_borrowed)
            ):
                break
        return relay

    @staticmethod
    def _relay_settled(relay: dict, prev: dict, borrowed: set[str]) -> bool:
        """Each BORROWED node's relayed T and Pt sensitivities equal what the
        pass used for it (the previous pass's; missing = 0) to 1e-9 -- the
        Jacobian's own accuracy. The composition sensitivities iterate with
        the same contraction and are not compared."""
        for nid in borrowed:
            cols = relay.get(nid, {})
            pcols = prev.get(nid, {})
            for idx in cols.keys() | pcols.keys():
                pkg, ppkg = cols.get(idx, {}), pcols.get(idx, {})
                for key in ("T", "Pt", "Pt_mix"):
                    a, b = pkg.get(key, 0.0), ppkg.get(key, 0.0)
                    if abs(a - b) > 1e-9 * max(abs(a), abs(b)):
                        return False
        return True

    def _relay_mix(
        self,
        out: dict,
        nid: str,
        stream_info: list[tuple],
        mix_res: Any,
        wall_contributions: list[tuple],
        x: np.ndarray,
        relay: dict,
        up: Any,
    ) -> None:
        """Chain one mixer's sensitivities into ``out`` (relay[node] layout):
        its streams (by flow into ``nid`` and through each stream's source
        node) and the wall heat it takes. Used for every node, and for a
        heated duct's exit state (its own stream and its own wall heat, #672).
        """
        n_species = self._n_species
        # Iterate over stream_info (one entry per upstream stream, including
        # multiple entries for multi-port elements like TeeJunctionElement).
        for i, (s_up, elem, src_nid, coeffs) in enumerate(stream_info):
            if elem is None:
                continue  # a fixed injection: nothing to relay
            # A turbomachine's stream arrives at Tt_out(Tt_in, Y), f Pt_in (#893).
            tm = getattr(s_up, "_turbomachine", None)
            t_jac = mix_res.dTt_mix_d_stream[i]
            pt_jac = mix_res.dPt_mix_d_stream[i]
            # dY_mix_d_stream is indexed by [species][stream]
            y_jacs = [mix_res.dY_mix_d_stream[k][i] for k in range(n_species)]

            # 1. Direct dependency on upstream mass flow (if it's a solver unknown)
            m_indices = self._unknown_indices.get(elem.id)
            if m_indices:
                # flow_jac_at_node maps each unknown index to its coefficient
                # in d(stream_i.m_dot)/d(x). For 2-port: {idx: 1.0}. For tee:
                # may involve 2 unknowns (m_dot_com and m_dot_branch).
                # The stream's flow is the one INTO this node for a
                # single-source element (flow_at_node(nid), as set in
                # stream_info above), so its derivative must be taken at
                # nid too. Taking it at src_nid gave a branching tee's
                # straight sink d/d(m_com) only and dropped the -1 on
                # m_branch.
                jac_node = nid if len(elem.all_source_nodes()) == 1 else src_nid
                flow_jac = (
                    coeffs if coeffs is not None else elem.flow_jac_at_node(jac_node, m_indices)
                )
                for idx, coeff in flow_jac.items():
                    node_relay = out.setdefault(
                        idx, {"T": 0.0, "Y": np.zeros(n_species), "Pt_mix": 0.0}
                    )
                    node_relay["T"] += t_jac.d_mdot * coeff
                    node_relay["Pt_mix"] += pt_jac.d_mdot * coeff
                    for k in range(n_species):
                        node_relay["Y"][k] += y_jacs[k].d_mdot * coeff

            # 2. Recursive dependency on upstream P/T/Y via src_nid relay
            if src_nid in relay:
                for idx, sens_src in up(src_nid).items():
                    node_relay = out.setdefault(
                        idx, {"T": 0.0, "Y": np.zeros(n_species), "Pt_mix": 0.0}
                    )
                    sens_up = sens_src
                    if tm is not None:
                        dT_dT, dT_dY, f_pt = tm
                        sens_up = {
                            "T": dT_dT * sens_src["T"] + float(np.dot(dT_dY, sens_src["Y"])),
                            "Y": sens_src["Y"],
                            "Pt": f_pt * sens_src.get("Pt", 0.0),
                        }
                    # dT/dx = sum_i( dT/dT_in_i * dT_in_i/dx + dT/dP_in_i * dP_in_i/dx + dT/dY_in_i * dY_in_i/dx )
                    node_relay["T"] += t_jac.d_Tt * sens_up["T"]
                    # Note: stagnation Pt doesn't usually have T/Y dependency, but we relay it if present
                    node_relay["T"] += t_jac.d_Pt * sens_up.get("Pt", 0.0)
                    node_relay["T"] += np.dot(t_jac.d_Y, sens_up["Y"])

                    # dP_total_mix/dx
                    node_relay["Pt_mix"] += pt_jac.d_Tt * sens_up["T"]
                    node_relay["Pt_mix"] += pt_jac.d_Pt * sens_up.get("Pt", 0.0)
                    node_relay["Pt_mix"] += np.dot(pt_jac.d_Y, sens_up["Y"])

                    # dY_k/dx
                    for k in range(n_species):
                        node_relay["Y"][k] += y_jacs[k].d_Tt * sens_up["T"]
                        node_relay["Y"][k] += y_jacs[k].d_Pt * sens_up.get("Pt", 0.0)
                        node_relay["Y"][k] += np.dot(y_jacs[k].d_Y, sens_up["Y"])

        # Wall coupling relay extension: dT_node/dx += dT_mix/dQ * dQ/dx
        # Full chain rule with channel Jacobians:
        #   dQ/dmdot = dQ/dh_a * dh_a/dmdot + dQ/dT_aw_a * dT_aw_a/dmdot + (side B terms)
        #   dQ/dT    = dQ/dh_a * dh_a/dT    + dQ/dT_aw_a * dT_aw_a/dT    + (side B terms)
        if wall_contributions:
            # Wall Q is a heat rate [W]: dT/dQ = (1/cp) / m_eff, m_eff
            # the flow the mixer spreads Q over (dTt_mix_dQ, exact also
            # below kMixerHeatMdotFloor and for a non-positive total).
            dTt_mix_dQ = mix_res.dTt_mix_dQ
            for (
                _wall,
                _Qi,
                wall_result,
                ea_id,
                eb_id,
                is_side_a,
                ch_a,
                ch_b,
                (dCa_dm, dCa_dTa, dCa_dTb, dCb_dm, dCb_dTb, dCb_dTa, dCa_dYa, dCb_dYb),
            ) in wall_contributions:
                sign = -1.0 if is_side_a else 1.0

                # Side A contributions
                obj_a = self.network.elements.get(ea_id) or self.network.nodes.get(ea_id)
                from_a = self._flow_ends(obj_a, x)[0]

                # Temperature coupling: dT_node/dT_upstream_a
                if from_a in relay:
                    # dQ/dT_aw_a * dT_aw_a/dT + dQ/dh_a * dh_a/dT
                    factor_T_a = (
                        dTt_mix_dQ
                        * sign
                        * (
                            wall_result.dQ_dT_aw_a * ch_a.dT_aw_dT
                            + wall_result.dQ_dh_a * ch_a.dh_dT
                            + wall_result.dQ_dC_a * dCa_dTa
                            + wall_result.dQ_dC_b * dCb_dTa
                        )
                    )
                    # And its composition (#717): h, T_aw and C move with Y.
                    # The radiation's gas temperature, from the stream (#793)
                    rad_T_a = getattr(wall_result, "rad_dQ_dTt_a", 0.0)
                    if rad_T_a:
                        factor_T_a += dTt_mix_dQ * sign * rad_T_a
                    ups_a = up(from_a)
                    factor_Y_a = self._wall_side_dQ_dY(wall_result, ch_a, dCa_dYa, side_a=True)
                    # and its fuel/air ratio, from the composition (#793)
                    rad_Y_a = getattr(wall_result, "rad_dQ_dY_a", None)
                    if rad_Y_a is not None:
                        arr = np.asarray(rad_Y_a, dtype=float)
                        factor_Y_a = arr if factor_Y_a is None else factor_Y_a + arr
                    for idx, sens_up in ups_a.items():
                        nr = out.setdefault(idx, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0})
                        nr["T"] += factor_T_a * sens_up["T"]
                        if factor_Y_a is not None:
                            nr["T"] += dTt_mix_dQ * sign * float(np.dot(factor_Y_a, sens_up["Y"]))

                # Mass-flow coupling: dT_node/dmdot_a (direct)
                for idx_mdot_a, dm_a in self._htc_flow_columns(ea_id):
                    # dQ/dh_a * dh_a/dmdot + dQ/dT_aw_a * dT_aw_a/dmdot
                    factor_mdot_a = (
                        dTt_mix_dQ
                        * sign
                        * self._htc_mdot_sign(obj_a, x)
                        * (
                            wall_result.dQ_dh_a * ch_a.dh_dmdot
                            + wall_result.dQ_dT_aw_a * ch_a.dT_aw_dmdot
                        )
                    )
                    # The capacity rate |m| cp (signed in the unknown)
                    factor_mdot_a += dTt_mix_dQ * sign * wall_result.dQ_dC_a * dCa_dm
                    nr = out.setdefault(idx_mdot_a, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0})
                    nr["T"] += factor_mdot_a * dm_a

                # Neighbour flows the side-A htc reads (#465)
                self._relay_flow_inputs(out, obj_a, ch_a, dTt_mix_dQ * sign * wall_result.dQ_dh_a)
                # The film's provider flow and coolant temperature (#386)
                self._relay_film(
                    out, _wall, ch_a, dTt_mix_dQ * sign * wall_result.dQ_dT_aw_a, relay, up
                )
                # and the film's own h on that side (#962)
                self._relay_film_h(
                    out, _wall, ch_a, dTt_mix_dQ * sign * wall_result.dQ_dh_a, relay, up
                )
                # A reference temperature on neighbours' flows and T (#467)
                self._relay_reference_inputs(
                    out, obj_a, ch_a, dTt_mix_dQ * sign * wall_result.dQ_dT_aw_a, relay, up
                )

                # Pressure coupling (#496): the side's T_aw moves with its feed
                # node's static P at fixed flow (v ~ 1/rho); h does not.
                # A momentum chamber's h does too, through its static T.
                idx_P_a = self._name_to_index.get(f"{from_a}.P")
                dQ_dP_a = wall_result.dQ_dT_aw_a * getattr(
                    ch_a, "dT_aw_dP", 0.0
                ) + wall_result.dQ_dh_a * getattr(ch_a, "dh_dP", 0.0)
                # The radiation's gas pressure (#793)
                rad_P_a = getattr(wall_result, "rad_dQ_dP_a", 0.0)
                if rad_P_a:
                    dQ_dP_a += rad_P_a
                if idx_P_a is not None and dQ_dP_a != 0.0:
                    nr = out.setdefault(idx_P_a, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0})
                    nr["T"] += dTt_mix_dQ * sign * dQ_dP_a
                # A duct's T_aw moves with its feed's Pt instead (#636).
                idx_Pt_a = self._name_to_index.get(f"{from_a}.Pt")
                # and so does its h, through the duct's static T.
                dQ_dPt_a = wall_result.dQ_dT_aw_a * getattr(
                    ch_a, "dT_aw_dPt", 0.0
                ) + wall_result.dQ_dh_a * getattr(ch_a, "dh_dPt", 0.0)
                if idx_Pt_a is not None and dQ_dPt_a != 0.0:
                    nr = out.setdefault(idx_Pt_a, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0})
                    nr["T"] += dTt_mix_dQ * sign * dQ_dPt_a

                # Side B contributions: none from a side with no stream (#970)
                if eb_id is None:
                    continue
                obj_b = self.network.elements.get(eb_id) or self.network.nodes.get(eb_id)
                from_b = self._flow_ends(obj_b, x)[0]

                # Temperature coupling: dT_node/dT_upstream_b
                if from_b in relay:
                    factor_T_b = (
                        dTt_mix_dQ
                        * sign
                        * (
                            wall_result.dQ_dT_aw_b * ch_b.dT_aw_dT
                            + wall_result.dQ_dh_b * ch_b.dh_dT
                            + wall_result.dQ_dC_b * dCb_dTb
                            + wall_result.dQ_dC_a * dCa_dTb
                        )
                    )
                    # The radiation's gas temperature, from the stream (#793)
                    rad_T_b = getattr(wall_result, "rad_dQ_dTt_b", 0.0)
                    if rad_T_b:
                        factor_T_b += dTt_mix_dQ * sign * rad_T_b
                    ups_b = up(from_b)
                    factor_Y_b = self._wall_side_dQ_dY(wall_result, ch_b, dCb_dYb, side_a=False)
                    # and its fuel/air ratio, from the composition (#793)
                    rad_Y_b = getattr(wall_result, "rad_dQ_dY_b", None)
                    if rad_Y_b is not None:
                        arr = np.asarray(rad_Y_b, dtype=float)
                        factor_Y_b = arr if factor_Y_b is None else factor_Y_b + arr
                    for idx, sens_up in ups_b.items():
                        nr = out.setdefault(idx, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0})
                        nr["T"] += factor_T_b * sens_up["T"]
                        if factor_Y_b is not None:
                            nr["T"] += dTt_mix_dQ * sign * float(np.dot(factor_Y_b, sens_up["Y"]))

                # Mass-flow coupling: dT_node/dmdot_b (direct)
                for idx_mdot_b, dm_b in self._htc_flow_columns(eb_id):
                    factor_mdot_b = (
                        dTt_mix_dQ
                        * sign
                        * self._htc_mdot_sign(obj_b, x)
                        * (
                            wall_result.dQ_dh_b * ch_b.dh_dmdot
                            + wall_result.dQ_dT_aw_b * ch_b.dT_aw_dmdot
                        )
                    )
                    # The capacity rate |m| cp (signed in the unknown)
                    factor_mdot_b += dTt_mix_dQ * sign * wall_result.dQ_dC_b * dCb_dm
                    nr = out.setdefault(idx_mdot_b, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0})
                    nr["T"] += factor_mdot_b * dm_b

                self._relay_flow_inputs(out, obj_b, ch_b, dTt_mix_dQ * sign * wall_result.dQ_dh_b)
                self._relay_film(
                    out, _wall, ch_b, dTt_mix_dQ * sign * wall_result.dQ_dT_aw_b, relay, up
                )
                # and the film's own h on that side (#962)
                self._relay_film_h(
                    out, _wall, ch_b, dTt_mix_dQ * sign * wall_result.dQ_dh_b, relay, up
                )
                # A reference temperature on neighbours' flows and T (#467)
                self._relay_reference_inputs(
                    out, obj_b, ch_b, dTt_mix_dQ * sign * wall_result.dQ_dT_aw_b, relay, up
                )

                # Pressure coupling (#496): the side's T_aw moves with its feed
                # node's static P at fixed flow (v ~ 1/rho); h does not.
                # A momentum chamber's h does too, through its static T.
                idx_P_b = self._name_to_index.get(f"{from_b}.P")
                dQ_dP_b = wall_result.dQ_dT_aw_b * getattr(
                    ch_b, "dT_aw_dP", 0.0
                ) + wall_result.dQ_dh_b * getattr(ch_b, "dh_dP", 0.0)
                # The radiation's gas pressure (#793)
                rad_P_b = getattr(wall_result, "rad_dQ_dP_b", 0.0)
                if rad_P_b:
                    dQ_dP_b += rad_P_b
                if idx_P_b is not None and dQ_dP_b != 0.0:
                    nr = out.setdefault(idx_P_b, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0})
                    nr["T"] += dTt_mix_dQ * sign * dQ_dP_b
                # A duct's T_aw moves with its feed's Pt instead (#636).
                idx_Pt_b = self._name_to_index.get(f"{from_b}.Pt")
                # and so does its h, through the duct's static T.
                dQ_dPt_b = wall_result.dQ_dT_aw_b * getattr(
                    ch_b, "dT_aw_dPt", 0.0
                ) + wall_result.dQ_dh_b * getattr(ch_b, "dh_dPt", 0.0)
                if idx_Pt_b is not None and dQ_dPt_b != 0.0:
                    nr = out.setdefault(idx_Pt_b, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0})
                    nr["T"] += dTt_mix_dQ * sign * dQ_dPt_b

    def _propagate_pass(
        self,
        x: np.ndarray,
        order: list[str],
        reversed_ids: set[str],
        with_relay: bool = True,
        prev_relay: dict | None = None,
    ) -> dict:
        """One pass of the state propagation in ``order`` (see _propagate_states)."""
        # relay[node_id][global_unknown_index] = {'T': dT/dx, 'Y': [dY/dx], 'Pt_mix': dP_total_mix/dx}
        relay = {nid: {} for nid in self.network.nodes}
        # A node read ahead of its order (a back-edge) gives the previous
        # pass's sensitivities, as it gives its previous state -- and so does
        # the node being processed, read through its own walls (a self-loop:
        # its relay is still being built). The passes iterate T' = a + b T'
        # to its fixed point.
        done: set[str] = set()
        current: list[str] = [""]
        prev_relay = prev_relay or {}

        def up(n: str) -> dict:
            if n in done and n != current[0]:
                return relay[n]
            self._relay_borrowed.add(n)
            return prev_relay.get(n, {})

        has_walls = self.network.thermal_coupling_enabled and bool(self.network.walls)
        film_discharges = self._film_discharges() if has_walls else {}
        # Each pass derives the heated ducts' exit states afresh (#672): one
        # that no node takes this pass (no flow into it) reads as adiabatic.
        for el in self.network.elements.values():
            if getattr(el, "_rayleigh_mix_ends", None):
                el._rayleigh_exit = None

        for nid in order:
            self._eval_where = f"node '{nid}'"
            node = self.network.nodes[nid]
            done.add(nid)
            current[0] = nid

            # Elements delivering INTO nid by the sign of their flow (#481). A
            # junction with a reversed port delivers by the sign of its port
            # flows instead, below (#517).
            rev_jct = getattr(self, "_rev_junctions", {})
            up_elems = [
                e
                for e in self.network.get_upstream_elements(nid)
                if e.id not in reversed_ids and e.id not in rev_jct
            ]
            rev_in = [e for e in self.network.get_downstream_elements(nid) if e.id in reversed_ids]
            jct_in = [
                e
                for e in self.network.get_upstream_elements(nid)
                + self.network.get_downstream_elements(nid)
                if e.id in rev_jct
            ]

            # Boundaries with NO inflow are pure sources
            if node.is_flow_boundary and not up_elems and not rev_in and not jct_in:
                T, Y, _ = node.compute_derived_state([])
                self._derived_states[nid] = (T, Y, None)
                if not with_relay:
                    continue
                # Seed relay for boundary's own unknowns
                node_indices = self._unknown_indices.get(nid, [])
                unknown_names = self.unknown_names
                for global_idx in node_indices:
                    name = unknown_names[global_idx]
                    if name.endswith(".Pt"):
                        relay[nid].setdefault(
                            global_idx, {"T": 0.0, "Y": np.zeros(self._n_species), "Pt": 0.0}
                        )["Pt"] = 1.0
                continue

            # stream_info: (state, elem, src_nid) - one entry per upstream stream.
            # Multi-port elements (e.g. TeeJunctionElement) contribute multiple entries.
            stream_info: list[tuple] = []
            for elem in up_elems:
                indices = self._unknown_indices.get(elem.id, [])
                # Compute named m_dot Jacobian for this element once (flow into nid).
                # Stored on each upstream state so MomentumChamberNode can build the
                # correct d_res/dm_dot Jacobian even when the element has named unknowns
                # like tee.m_dot_com / tee.m_dot_branch instead of a plain elem.m_dot.
                elem_jac_names: dict[str, float] = {}
                if indices:
                    for col, coeff in elem.flow_jac_at_node(nid, indices).items():
                        if col < len(self.unknown_names):
                            elem_jac_names[self.unknown_names[col]] = coeff
                _single_source = len(elem.all_source_nodes()) == 1
                for src_nid in elem.all_source_nodes():
                    state_up = self._get_node_state_with_prev(self.network.nodes[src_nid], x)
                    if indices:
                        # For single-source elements (channels, branching tees) use the
                        # flow INTO the current sink node (nid) so that each sink's
                        # _total_m_dot reflects its actual partial flow, not the
                        # supplier's total.  For multi-source elements (merging tees) each
                        # source contributes its own partial flow, so keep flow_at_node
                        # called with src_nid.
                        if _single_source:
                            state_up.m_dot = elem.flow_at_node(nid, x, indices)
                        else:
                            state_up.m_dot = elem.flow_at_node(src_nid, x, indices)
                        state_up._element_id = elem.id
                        state_up._m_dot_jac_names = elem_jac_names
                    state_up._src_node = src_nid
                    if isinstance(elem, TurbomachineElement):
                        # A compressor or turbine delivers its stream at its
                        # outlet stagnation state (#893); _relay_mix chains
                        # Tt_out through Tt_in and Y, and Pt_out through Pt_in.
                        Tt_out, dT_dT, dT_dY, f_pt = elem.outlet_stagnation(state_up)
                        state_up.Tt = Tt_out
                        state_up.Pt = f_pt * float(state_up.Pt)
                        state_up._turbomachine = (dT_dT, np.asarray(dT_dY, dtype=float), f_pt)
                    stream_info.append((state_up, elem, src_nid, None))
            # A reversed element delivers -m_dot from its to_node (#481).
            for elem in rev_in:
                i_m = self._two_port_mdot_indices()[elem.id]
                state_up = self._get_node_state_with_prev(self.network.nodes[elem.to_node], x)
                state_up.m_dot = -float(x[i_m])
                state_up._element_id = elem.id
                state_up._m_dot_jac_names = {f"{elem.id}.m_dot": -1.0}
                state_up._src_node = elem.to_node
                stream_info.append((state_up, elem, elem.to_node, {i_m: -1.0}))
            # A junction with a reversed port, by the sign of its port flows.
            for elem in jct_in:
                flows, idxs = rev_jct[elem.id]
                streams = self._junction_streams(elem, nid, flows, idxs)
                # The node takes d(total inflow)/dx once per element (it
                # de-duplicates by element id), so every stream carries the
                # element's TOTAL, as flow_jac_at_node does for a merge.
                total_jac: dict[str, float] = {}
                for _, _, jac in streams:
                    for col, coeff in jac.items():
                        name = self.unknown_names[col]
                        total_jac[name] = total_jac.get(name, 0.0) + coeff
                for src_nid, m, jac in streams:
                    state_up = self._get_node_state_with_prev(self.network.nodes[src_nid], x)
                    state_up.m_dot = m
                    state_up._element_id = elem.id
                    state_up._m_dot_jac_names = total_jac
                    state_up._src_node = src_nid
                    stream_info.append((state_up, elem, src_nid, jac))
            # A MassFlowBoundary between elements injects its own stream: the
            # mass row already counts its m_dot, so the mixing must too, at
            # its own Tt and Y (#481). Fixed, so it has no relay column.
            if (
                isinstance(node, MassFlowBoundary)
                and self.network.get_upstream_elements(nid)
                and self.network.get_downstream_elements(nid)
            ):
                own = self._get_node_state(node, x)
                inj = NetworkMixtureState(
                    P=own.P,
                    Pt=own.Pt,
                    Tt=float(node.Tt),
                    m_dot=float(node.m_dot),
                    Y=node.Y if node.Y is not None else self._default_Y,
                )
                stream_info.append((inj, None, nid, {}))
            # The gas a film's holes discharge into, before their coolant
            # mixes in (#694): the T_inf the film's eta is defined against.
            for pid in film_discharges.get(nid, ()):
                if self._film_unmix(nid, pid, stream_info, x, relay, up, with_relay):
                    done.add(f"{FILM_INF_PREFIX}{pid}")

            no_inflow = not stream_info and not node.is_flow_boundary
            up_states = [s for s, _, _, _ in stream_info]

            # Wall Coupling: evaluate walls BEFORE compute_derived_state
            # so that Q flows through the EnergyBoundary abstraction.
            wall_contributions: list[tuple] = []
            if has_walls:
                wall_contributions = self._evaluate_walls_for_node(nid, up_elems, x)

            # A Tt_target boundary's heat is its unknown's value at x (#792).
            target_cols = self._target_q_columns.get(nid, ())
            for eb, i_q in target_cols:
                eb.Q_solved = float(x[i_q])

            T, Y, mix_res = node.compute_derived_state(up_states)
            if no_inflow:
                # A node whose flows all leave it carries no flow: continuous
                # with a vanishing last inflow. Nodes used to set 1 kg/s here,
                # a 4.3 kPa jump in a momentum chamber's dynamic head as its
                # feed reversed (#481). Its T/Y stay the nodes' fallback: a
                # continuous fallback (the sources it would be fed from,
                # weighted 1/|m|) was measured and REJECTED -- it doubled the
                # cold-start time of hot impingement arrays and failed one,
                # with no gain elsewhere.
                node._total_m_dot = 0.0
                node._inflows = []
                node._upstream_m_dot_jac = {}
                node._upstream_element_ids = []
            # No clamp: the node carries its T(h) as is. The gas is ideal down
            # to 0 K (the thermo's low tail); at or below it T is the mixer's
            # continuation on h's line, derivatives consistent with it: the
            # stream holds less enthalpy than at 0 K, no gas state, which the
            # post-solve check reports. A [200, 5000] K clamp here used to
            # flatten T (zero sensitivity), discard the energy, and flag a node
            # fed by a clamped one or not by the last bit.
            T = float(T)
            if not T > 0.0:
                self._T_no_gas.add(nid)
            self._derived_states[nid] = (T, Y, mix_res)

            if self._record_boundary_arrivals and isinstance(node, PressureBoundary):
                for s_up, elem, _src, _coeffs in stream_info:
                    if elem is None or not hasattr(elem, "to_node") or not s_up.m_dot > 0.0:
                        continue
                    own = [
                        wc for wc in wall_contributions if (wc[3] if wc[5] else wc[4]) == elem.id
                    ]
                    arrive = _solver_tools.mixer_from_streams_and_jacobians(
                        [cb.MassStream(s_up.m_dot, s_up.Tt, s_up.Pt, s_up.Y)],
                        Q=sum(float(wc[1]) for wc in own),
                        fraction=0.0,
                    )
                    self._boundary_arrivals[(nid, elem.id)] = (
                        float(arrive.Tt_mix),
                        list(arrive.Y_mix),
                    )

            # A wall-heated duct into a node whose T is not its exit Tt (the
            # node mixes other streams or heat, or is a boundary): its own exit
            # state, from its stream and its own wall heat alone, through the
            # same mixer and relay as a node (#672). The duct reads it for its
            # Fanno-Rayleigh momentum.
            for s_up, elem, src_nid, coeffs in stream_info:
                if elem is None or nid not in getattr(elem, "_rayleigh_mix_ends", ()):
                    continue
                own = [wc for wc in wall_contributions if (wc[3] if wc[5] else wc[4]) == elem.id]
                exit_mix = _solver_tools.mixer_from_streams_and_jacobians(
                    [cb.MassStream(s_up.m_dot, s_up.Tt, s_up.Pt, s_up.Y)],
                    Q=sum(float(wc[1]) for wc in own),
                    fraction=0.0,
                )
                elem._rayleigh_exit = (nid, float(exit_mix.Tt_mix), src_nid)
                if with_relay:
                    key = f"{RAYLEIGH_EXIT_PREFIX}{elem.id}"
                    relay[key] = {}
                    self._relay_mix(
                        relay[key],
                        nid,
                        [(s_up, elem, src_nid, coeffs)],
                        exit_mix,
                        own,
                        x,
                        relay,
                        up,
                    )

            # Store wall coupling debug info on the node
            if wall_contributions:
                dTt_mix_dQ = mix_res.dTt_mix_dQ if mix_res else 0.0
                node._wall_couplings = []
                for (
                    wall,
                    Q_inject,
                    wall_result,
                    ea_id,
                    eb_id,
                    is_side_a,
                    ch_a,
                    ch_b,
                    _caps,
                ) in wall_contributions:
                    # dQ_inject/dQ_wall = -1 if is_side_a, +1 otherwise
                    dT_dQ_signed = -dTt_mix_dQ if is_side_a else dTt_mix_dQ
                    node._wall_couplings.append(
                        {
                            "wall": wall,
                            "wall_result": wall_result,
                            "Q_inject": Q_inject,
                            "is_side_a": is_side_a,
                            "elem_a_id": ea_id,
                            "elem_b_id": eb_id,
                            "dTt_mix_dQ": dT_dQ_signed,
                            "ch_result_a": ch_a,
                            "ch_result_b": ch_b,
                        }
                    )

            if not with_relay:
                continue

            # Sensitivity Relay (Chain-Rule)
            if mix_res:
                self._relay_mix(
                    relay[nid], nid, stream_info, mix_res, wall_contributions, x, relay, up
                )
                # The Tt_target heat unknown's own column (#792): dT_mix/dQ
                # from C++ (the heat per unit flow over cp, fading to 0 at no
                # flow). The composition does not move with Q: the mixer's
                # Y_mix is Q-free and a combustor shifts its products' T at
                # frozen composition.
                for _eb, i_q in target_cols:
                    pkg = relay[nid].setdefault(
                        i_q, {"T": 0.0, "Y": np.zeros(self._n_species), "Pt_mix": 0.0}
                    )
                    pkg["T"] += float(mix_res.dTt_mix_dQ)

            # 3. Add own unknowns (Pt) to the relay
            node_unks = self._unknown_indices.get(nid, [])
            for unk_idx in node_unks:
                unk_name = self.unknown_names[unk_idx]
                if unk_name.endswith(".Pt"):
                    node_relay = relay[nid].setdefault(
                        unk_idx, {"T": 0.0, "Y": np.zeros(self._n_species), "Pt": 0.0}
                    )
                    node_relay["Pt"] = 1.0

        return relay

    @staticmethod
    def _main_face(element: Any, node_out: Any, state_out: NetworkMixtureState) -> Any:
        """Hand ``state_out`` the main-inlet face of a merge chamber (#471).

        The chamber's (P, Pt) is its outlet state; the declared main inlet
        arrives at the face state the impulse balance gives. Returns the face
        (for the Jacobian chain) or None when it does not apply.
        """
        if getattr(node_out, "main_inlet", None) != element.id:
            return None
        if not node_out.has_side_streams():
            return None
        # The element's own signed flow into the chamber: a reversed main
        # keeps its face (#493).
        face = node_out.main_face_state(state_out, m_main=float(state_out.m_dot))
        state_out.P = face[0]
        state_out.Pt = face[1]
        return face

    @staticmethod
    def _chain_main_face(elem_jac: dict, node_id: str, face: Any) -> None:
        """Replace the element's d/d(chamber P, Pt) by the face's chain.

        The element was handed the face (P_face, Pt_face) in place of the
        chamber's (P, Pt), so its entries for those names are really
        d(res)/d(P_face) and d(res)/d(Pt_face); the face depends on the
        chamber state and the flows through jac_P / jac_Pt.
        """
        _, _, jac_P, jac_Pt = face
        for row in elem_jac.values():
            c_P = row.pop(f"{node_id}.P", 0.0)
            c_Pt = row.pop(f"{node_id}.Pt", 0.0)
            for coeff, offs in ((c_P, jac_P), (c_Pt, jac_Pt)):
                if coeff == 0.0:
                    continue
                for k, v in offs.items():
                    row[k] = row.get(k, 0.0) + coeff * v

    def _network_flows(self, obj: Any, x: np.ndarray) -> dict[str, float] | None:
        """Flows an element reads besides its own, or None if it reads none.

        Opt-in: an element that defines ``network_flow_inputs()`` -- a list of
        (element id, node id) pairs -- gets each element's flow at that node,
        taken from the unknown vector, as ``flows=`` on ``residuals``,
        ``htc_and_T`` and ``diagnostics``. ImpingementPlateElement reads the
        crossflow approaching its row; ImpingementCrossflowElement reads the
        jets merging at its outlet (#465). Residual Jacobians name the
        neighbours' unknowns directly; an htc's neighbour sensitivity comes
        back as ``dh_dsources`` and is relayed by ``_relay_flow_inputs``.
        """
        inputs = getattr(obj, "network_flow_inputs", None)
        if inputs is None:
            return None
        flows: dict[str, float] = {}
        for eid, node_id in inputs():
            indices = self._unknown_indices.get(eid, [])
            if indices:
                flows[eid] = float(self.network.elements[eid].flow_at_node(node_id, x, indices))
        return flows

    def _node_inflows(self, nid: str, x: np.ndarray) -> list[tuple[float, float, list[float]]]:
        """(m_dot, Tt, Y) of every stream INTO nid, as the propagation mixes
        them: by flow sign for 2-port elements, a MassFlowBoundary's own
        injection included."""
        two = self._two_port_mdot_indices()
        out: list[tuple[float, float, list[float]]] = []

        def state(n: str) -> NetworkMixtureState:
            return self._get_node_state(self.network.nodes[n], x)

        rev_jct = self._reversed_junctions(x)
        for e in self.network.get_upstream_elements(nid) + self.network.get_downstream_elements(
            nid
        ):
            if e.id in rev_jct:
                for src, m, _ in self._junction_streams(e, nid, *rev_jct[e.id]):
                    st = state(src)
                    out.append((m, float(st.Tt), list(st.Y)))
        for e in self.network.get_upstream_elements(nid):
            if e.id in rev_jct:
                continue
            if e.id in two:
                m = float(x[two[e.id]])
                if m > 0.0:
                    st = state(e.from_node)
                    out.append((m, float(st.Tt), list(st.Y)))
                continue
            ind = self._unknown_indices.get(e.id, [])
            single = len(e.all_source_nodes()) == 1
            for src in e.all_source_nodes():
                try:
                    m = float(e.flow_at_node(nid if single else src, x, ind))
                except (IndexError, TypeError):
                    m = 0.0
                st = state(src)
                out.append((m, float(st.Tt), list(st.Y)))
        for e in self.network.get_downstream_elements(nid):
            if e.id in two and x[two[e.id]] < 0.0:
                st = state(e.to_node)
                out.append((-float(x[two[e.id]]), float(st.Tt), list(st.Y)))
        node = self.network.nodes[nid]
        if (
            isinstance(node, MassFlowBoundary)
            and self.network.get_upstream_elements(nid)
            and self.network.get_downstream_elements(nid)
        ):
            Y = node.Y if node.Y is not None else self._default_Y
            out.append((float(node.m_dot), float(node.Tt), list(Y)))
        return out

    def _htc_mdot_sign(self, obj: Any, x: np.ndarray) -> float:
        """Sign that turns an htc's own-flow derivative into one in the
        signed unknown. Surface correlations take |m_dot| (the element hands
        them abs(m_dot)), so for a reversed element d/dm = -d/d|m| (#481);
        an element whose derivative is already signed says so with
        ``_htc_dmdot_is_signed``."""
        if getattr(obj, "_htc_dmdot_is_signed", False):
            return 1.0
        i = self._two_port_mdot_indices().get(getattr(obj, "id", None))
        return -1.0 if i is not None and x[i] < 0.0 else 1.0

    def _htc_flow_columns(self, obj_id: str) -> list[tuple[int, float]]:
        """(unknown index, d(flow the htc reads)/d(unknown)) of a wall side.

        An element's htc reads its own flow, its first unknown. A node's reads
        its throughput, the sum of its inflows: its first unknown is its own
        P, which took the node's dh_dmdot until a walled combustor or chamber
        node showed a Jacobian off by its whole flow column.
        """
        if obj_id in self.network.nodes:
            jac = getattr(self.network.nodes[obj_id], "_upstream_m_dot_jac", None) or {}
            return [
                (self._name_to_index[name], float(c))
                for name, c in jac.items()
                if name in self._name_to_index
            ]
        ind = self._unknown_indices.get(obj_id, [])
        return [(ind[0], 1.0)] if ind else []

    def _htc_and_T(
        self, obj: Any, state: NetworkMixtureState, x: np.ndarray, lagged: bool = True
    ) -> Any:
        flows = self._network_flows(obj, x)
        temps = self._network_temperatures(obj, x, lagged)
        if flows is None:
            return obj.htc_and_T(state)
        if temps is None:
            return obj.htc_and_T(state, flows=flows)
        return obj.htc_and_T(state, flows=flows, temperatures=temps)

    def _network_temperatures(
        self, obj: Any, x: np.ndarray, lagged: bool = True
    ) -> dict[str, float] | None:
        """Node temperatures an element's htc reads, or None if it reads none.

        Opt-in, like ``network_flow_inputs``: an element that defines
        ``network_temperature_inputs()`` -- node ids -- gets each node's total
        temperature as ``temperatures=``. ImpingementPlateElement reads its
        bypass node's (#467). ``lagged`` reads a node not yet derived this
        pass from the previous one, as the wall reads do; the T_aw
        sensitivity comes back as ``dT_aw_dT_nodes``, relayed by
        ``_relay_reference_inputs``.
        """
        inputs = getattr(obj, "network_temperature_inputs", None)
        if inputs is None:
            return None
        out: dict[str, float] = {}
        for nid in inputs():
            node = self.network.nodes[nid]
            st = (
                self._get_node_state_with_prev(node, x) if lagged else self._get_node_state(node, x)
            )
            out[nid] = float(st.Tt)
        return out

    def _film_discharges(self) -> dict[str, list[str]]:
        """{node: film providers discharging into it} for the walls' films
        whose node only mixes (a plenum or momentum chamber, #694). A
        combustor's or a boundary's state is not its inflows mixed, so no
        undiluted gas is taken there: the film reads the side's own state. A
        provider whose sources discharge into different nodes (a liner with a
        gas node per station, #695) has none either."""
        out: dict[str, list[str]] = {}
        for wall in self.network.walls.values():
            if wall.film is None:
                continue
            to_node = film_discharge_node(self.network.film_provider(wall.film))
            node = self.network.nodes.get(to_node) if to_node is not None else None
            if isinstance(node, PlenumNode | MomentumChamberNode):
                pids = out.setdefault(to_node, [])
                if wall.film not in pids:
                    pids.append(wall.film)
        return out

    def _film_unmix(
        self,
        nid: str,
        pid: str,
        stream_info: list[tuple],
        x: np.ndarray,
        relay: dict,
        up: Any,
        with_relay: bool,
    ) -> bool:
        """Mix ``nid``'s inflows except provider ``pid``'s (every one of its
        sources) into the undiluted gas its film sees (#694), with the node's
        energy boundaries and the same mixer and relay as the node; wall heat
        into the node is not included. False (nothing stored) when the holes
        do not deliver into nid now or nothing else does."""
        own = {src.element_id for src in film_sources(self.network.film_provider(pid))}
        others = [st for st in stream_info if st[1] is None or st[1].id not in own]
        if len(others) == len(stream_info) or not others:
            return False
        if sum(float(st[0].m_dot) for st in others) <= 0.0:
            return False
        node = self.network.nodes[nid]
        ebs = getattr(node, "energy_boundaries", [])
        mix = _solver_tools.mixer_from_streams_and_jacobians(
            [cb.MassStream(st.m_dot, st.Tt, st.Pt, st.Y) for st, _, _, _ in others],
            Q=sum(eb.Q for eb in ebs),
            fraction=sum(eb.fraction for eb in ebs),
        )
        key = f"{FILM_INF_PREFIX}{pid}"
        self._film_inf_states[key] = (float(mix.Tt_mix), list(mix.Y_mix))
        if with_relay:
            relay[key] = {}
            self._relay_mix(relay[key], nid, others, mix, [], x, relay, up)
        return True

    def _film_inf(self, provider: Any, obj: Any, x: np.ndarray, lagged: bool) -> Any:
        """(T, Y) of the undiluted gas the film on ``obj`` sees, or None when
        ``obj`` is not fed by the node the holes discharge into (#694)."""
        to_node = film_discharge_node(provider)
        if to_node is None or self._flow_ends(obj, x)[0] != to_node:
            return None
        key = f"{FILM_INF_PREFIX}{provider.id}"
        if key in self._film_inf_states:
            return self._film_inf_states[key]
        if lagged and key in self._prev_film_inf_states:
            self._lagged_read = True
            return self._prev_film_inf_states[key]
        return None

    def _wall_side_htc(
        self,
        wall: Any,
        obj: Any,
        state: NetworkMixtureState,
        x: np.ndarray,
        lagged: bool = True,
    ) -> Any:
        """A wall side's surface result, under the wall's film if it covers
        this side (#386). Every place that reads a wall side goes through
        here, so the residual, the relay and the reported values agree.

        The film reads each source's flow and its supply node's temperature;
        ``lagged`` reads those nodes as the propagation does (a back-edge gives
        the previous pass's state).
        """
        r = self._htc_and_T(obj, state, x, lagged) if obj.has_convective_surface else None
        if r is None or wall.film is None or wall.film_side_id != obj.id:
            return r
        provider = self.network.film_provider(wall.film)
        coolant = []
        for src in film_sources(provider):
            ind = self._unknown_indices.get(src.element_id, [])
            if not ind:
                return r
            node = self.network.nodes[src.supply_node]
            c_state = (
                self._get_node_state_with_prev(node, x) if lagged else self._get_node_state(node, x)
            )
            coolant.append((c_state, float(x[ind[0]])))
        m_gas = (
            float(getattr(obj, "_total_m_dot", 0.0))
            if obj.id in self.network.nodes
            else float(state.m_dot)
        )
        return film_on_wall_side(
            provider,
            obj,
            r,
            state,
            coolant,
            m_gas,
            wall.film_start,
            wall.film_end,
            gas_inf=self._film_inf(provider, obj, x, lagged),
        )

    def _relay_film(
        self, out: dict, wall: Any, ch: Any, factor: float, relay: dict, up: Any
    ) -> None:
        """Relay a filmed side's T_aw onto each source's flow and its supply
        node's temperature (#386, #695). ``factor`` is dT_node/dT_aw of that
        side."""
        d_ms = getattr(ch, "dT_aw_dfilm_mdot", None)
        if d_ms is None or factor == 0.0:
            return
        provider = self.network.film_provider(wall.film)
        n_species = self._n_species
        for src, d_m, d_T in zip(film_sources(provider), d_ms, ch.dT_aw_dT_coolant, strict=True):
            ind = self._unknown_indices.get(src.element_id, [])
            if ind and d_m != 0.0:
                nr = out.setdefault(ind[0], {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0})
                nr["T"] += factor * d_m
            if src.supply_node in relay and d_T != 0.0:
                for idx, sens_up in up(src.supply_node).items():
                    nr = out.setdefault(idx, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0})
                    nr["T"] += factor * d_T * sens_up["T"]
        # The undiluted gas temperature (#694), through its own mix's relay.
        d_inf = float(getattr(ch, "dT_aw_dT_inf", 0.0))
        if d_inf != 0.0:
            for idx, sens_up in up(f"{FILM_INF_PREFIX}{provider.id}").items():
                nr = out.setdefault(idx, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0})
                nr["T"] += factor * d_inf * sens_up["T"]

    def _relay_film_h(
        self, out: dict, wall: Any, ch: Any, factor: float, relay: dict, up: Any
    ) -> None:
        """Relay a filmed side's film-side h (#962) onto each source's flow
        and its supply node's temperature, as _relay_film does its T_aw.
        ``factor`` is dT_node/dh of that side."""
        d_ms = getattr(ch, "dh_dfilm_mdot", None)
        if d_ms is None or factor == 0.0:
            return
        provider = self.network.film_provider(wall.film)
        n_species = self._n_species
        for src, d_m, d_T in zip(film_sources(provider), d_ms, ch.dh_dT_coolant, strict=True):
            ind = self._unknown_indices.get(src.element_id, [])
            if ind and d_m != 0.0:
                nr = out.setdefault(ind[0], {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0})
                nr["T"] += factor * d_m
            if src.supply_node in relay and d_T != 0.0:
                for idx, sens_up in up(src.supply_node).items():
                    nr = out.setdefault(idx, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0})
                    nr["T"] += factor * d_T * sens_up["T"]

    def _relay_reference_inputs(
        self, out: dict, obj: Any, ch: Any, factor: float, relay: dict, up: Any
    ) -> None:
        """Relay a side's T_aw onto the neighbour flows (``dT_aw_dsources``)
        and node temperatures (``dT_aw_dT_nodes``) it reads (#467): an
        impingement row's reference temperature moves with its bypass share
        and the bypass temperature. ``factor`` is dT_node/dT_aw of that side.
        """
        if factor == 0.0:
            return
        n_species = self._n_species
        d_src = getattr(ch, "dT_aw_dsources", None)
        inputs = getattr(obj, "network_flow_inputs", None)
        if d_src and inputs is not None:
            for eid, node_id in inputs():
                indices = self._unknown_indices.get(eid, [])
                if not indices or eid not in d_src:
                    continue
                elem = self.network.elements[eid]
                for col, coeff in elem.flow_jac_at_node(node_id, indices).items():
                    nr = out.setdefault(col, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0})
                    nr["T"] += factor * d_src[eid] * coeff
        for nid, d_T in (getattr(ch, "dT_aw_dT_nodes", None) or {}).items():
            if nid not in relay or d_T == 0.0:
                continue
            for idx, sens_up in up(nid).items():
                nr = out.setdefault(idx, {"T": 0.0, "Y": np.zeros(n_species), "Pt": 0.0})
                nr["T"] += factor * d_T * sens_up["T"]

    def _relay_flow_inputs(self, node_relay: dict, obj: Any, ch: Any, dT_dh: float) -> None:
        """Relay ``dh_dsources`` onto the neighbours' flow unknowns.

        The wall's own-flow path carries ``dh_dmdot`` to the element's first
        unknown only; an htc that reads a neighbour's flow would otherwise
        lose that column of the Jacobian.
        """
        dh_ds = getattr(ch, "dh_dsources", None)
        inputs = getattr(obj, "network_flow_inputs", None)
        if not dh_ds or inputs is None:
            return
        for eid, node_id in inputs():
            indices = self._unknown_indices.get(eid, [])
            if not indices or eid not in dh_ds:
                continue
            elem = self.network.elements[eid]
            for col, coeff in elem.flow_jac_at_node(node_id, indices).items():
                nr = node_relay.setdefault(
                    col, {"T": 0.0, "Y": np.zeros(self._n_species), "Pt": 0.0}
                )
                nr["T"] += dT_dh * dh_ds[eid] * coeff

    def _get_node_state_with_prev(self, node: NetworkNode, x: np.ndarray) -> NetworkMixtureState:
        """The node's state, its T/Y from the PREVIOUS pass if this pass has
        not derived it yet (a back-edge: a wall or stream read ahead of the
        propagation order). Such a read is flagged, and _propagate_states
        repeats the pass until the lagged states settle, so the residual is a
        function of x alone (#481). Pressures come from x with the same
        floors as ``_get_node_state``.
        """
        state = self._get_node_state(node, x)
        if node.is_flow_boundary or node.id in self._derived_states:
            return state
        self._lagged_read = True
        if node.id in self._prev_derived_states:
            T_prev, Y_prev, _ = self._prev_derived_states[node.id]
            state.Tt = float(T_prev)
            state.Y = Y_prev
        return state

    def _report_moving_nodes(
        self, x: np.ndarray, sol: dict[str, Any]
    ) -> dict[str, tuple[float, float]]:
        """Static T of the nodes that hold Pt = P while their gas moves (#782
        audit): a crossflow station (the segment's mean-flow state, as its
        half-drops form it) and a combustor with an area. Kept for
        ``_report_state``; returns the stations' (Pt, Mach). Reporting only."""
        self._report_T: dict[str, float] = {}
        station_pt: dict[str, tuple[float, float]] = {}
        for eid, el in self.network.elements.items():
            report = getattr(el, "station_report", None)
            prev = getattr(el, "_prev", None)
            if report is None or prev is None:
                continue
            node = self.network.nodes[el.from_node]
            try:
                r = report(
                    self._get_node_state(node, x),
                    float(sol.get(f"{prev.id}.m_dot", 0.0)),
                    float(sol.get(f"{eid}.m_dot", 0.0)),
                )
            except (ValueError, RuntimeError):
                r = None
            if r is not None:
                nid, T_s, Pt_s, M_s = r
                self._report_T[nid] = T_s
                station_pt[nid] = (Pt_s, M_s)
        for nid, node in self.network.nodes.items():
            static_T = getattr(node, "static_T", None)
            if static_T is None or nid in self._report_T:
                continue
            with contextlib.suppress(ValueError, RuntimeError):
                self._report_T[nid] = float(static_T(self._get_node_state(node, x)))
        return station_pt

    def _report_state(self, node: NetworkNode, x: np.ndarray) -> NetworkMixtureState:
        """``_get_node_state`` with the reported static T of a moving node
        that holds Pt = P (``_report_moving_nodes``), and the reported static
        P of a choked-bounded junction (#788, ``_choked_bounded_junction``).
        For diagnostics only: no residual reads either."""
        state = self._get_node_state(node, x)
        bounded = getattr(node, "_choked_bounded", None)
        if bounded is not None:
            state.P = float(bounded["P"])
        T_s = getattr(self, "_report_T", {}).get(node.id)
        if T_s is not None:
            state.__dict__["_T_static"] = T_s
        return state

    def _choked_bounded_junction(
        self, nid: str, x: np.ndarray
    ) -> tuple[dict[str, float] | None, str]:
        """(report, "") when momentum chamber ``nid``, at a root past its
        sonic limit, is a CHOKED-BOUNDED junction (#788); (None, why not)
        otherwise.

        The chamber's Pt is below the sonic Pt* of its mass flux at its area:
        no state of that area carries the flow at that Pt. It is the
        under-expanded exit of a choked duct when
          * exactly one element discharges into it, a plain compressible
            channel (the Fanno flow), choked (``choked_discharge_into``), of
            no larger area than the chamber: its sonic exit has Pt* at least
            the chamber's, and the jet loses the excess as it expands, as at
            a pressure boundary taking a choked duct ("choked bounded",
            #636); and
          * no row but the chamber's own reads its static pressure (its
            Jacobian column): the static P is then the closure's alone, the
            network's other rows are those of a plenum in its place, and the
            root is the plenum's.
        Mass and Tt pass the chamber unchanged; the loss Pt* - Pt sits in the
        channel's control volume (sigma = m R ln(Pt_in / Pt) > 0).

        The report states the chamber as the plenum its root is: P = Pt,
        the gas at rest, with the sonic state of its flux (P*, Pt*) and the
        expansion loss Pt* - Pt beside it. The channel reports its own
        choked exit (mach_out 1).
        """
        node = self.network.nodes[nid]
        two = self._two_port_mdot_indices()
        feeds = []
        for e in self.network.get_upstream_elements(nid) + self.network.get_downstream_elements(
            nid
        ):
            idx = two.get(e.id)
            if idx is None:
                return None, f"'{e.id}' ({type(e).__name__}) is not a two-port element"
            m = float(x[idx])
            if (e.to_node == nid and m > 0.0) or (e.from_node == nid and m < 0.0):
                feeds.append(e)
        if len(feeds) != 1:
            return None, f"{len(feeds)} elements discharge into it, not one choked channel"
        feed = feeds[0]
        if not isinstance(feed, ChannelElement):
            return None, f"'{feed.id}' ({type(feed).__name__}) feeds it, not a choked channel"
        st_in = self._get_node_state(self.network.nodes[feed.from_node], x)
        st_out = self._get_node_state(self.network.nodes[feed.to_node], x)
        st_in.m_dot = st_out.m_dot = float(x[two[feed.id]])
        why = feed.choked_discharge_into(nid, st_in, st_out)
        if why:
            return None, why
        if not float(feed.area or 0.0) <= node.area * (1.0 + 1e-9):
            return None, (
                f"'{feed.id}' ({float(feed.area or 0.0):.4g} m^2) is wider than the "
                f"chamber ({node.area:.4g} m^2), which would choke first"
            )
        _, J = self._residuals_and_jacobian(x)
        col = self.unknown_names.index(f"{nid}.P")
        column = J.tocsc()[:, col] if sp.issparse(J) else sp.csc_matrix(np.asarray(J)[:, [col]])
        own = set(self._unknown_indices.get(nid, []))
        readers = sorted(
            {self.unknown_names[i] for i, v in zip(column.indices, column.data, strict=True) if v}
            - {self.unknown_names[i] for i in own}
        )
        if readers:
            return None, f"row(s) {readers} read its static pressure"
        state = self._get_node_state(node, x)
        c = node._closure(state)
        if not c.bounded or not c.Pt_sonic > 0.0 or not state.Pt > 0.0:
            return None, "no sonic state at its flux"
        return {
            "P": float(state.Pt),
            "P_sonic": float(c.P_sonic),
            "Pt_sonic": float(c.Pt_sonic),
        }, ""

    def _get_node_state(self, node: NetworkNode, x: np.ndarray) -> NetworkMixtureState:
        """Constructs a NetworkMixtureState for a given node based on the current solver vector x."""
        default_Y = self._default_Y

        # Determine boundaries
        if isinstance(node, PressureBoundary):
            guess = getattr(node, "initial_guess", {})
            P_tot = guess.get(f"{node.id}.Pt", getattr(node, "Pt", 101325.0))
            P_stat = guess.get(f"{node.id}.P", P_tot)
            T_tot = guess.get(f"{node.id}.Tt", getattr(node, "Tt", 300.0))
            Y = getattr(node, "Y", None)
            if Y is None:
                Y = default_Y
            return NetworkMixtureState(P=float(P_stat), Pt=float(P_tot), Tt=float(T_tot), Y=Y)

        if isinstance(node, MassFlowBoundary):
            guess = getattr(node, "initial_guess", {})
            m = getattr(node, "m_dot", 0.0)
            T_tot = guess.get(f"{node.id}.Tt", getattr(node, "Tt", 300.0))
            Y = getattr(node, "Y", None)

            # If this boundary has upstream connections (sink), use derived mixed state
            if node.id in self._derived_states:
                T_tot_derived, Y_derived, _ = self._derived_states[node.id]
                T_tot = T_tot_derived
                Y = Y_derived

            if Y is None:
                Y = default_Y

            p_val = guess.get(f"{node.id}.P", 101325.0)
            ptot_val = guess.get(f"{node.id}.Pt", 101325.0)

            indices = self._unknown_indices.get(node.id, [])
            unknowns = node.unknowns()
            for i, unk in zip(indices, unknowns, strict=False):
                if unk.endswith(".P"):
                    p_val = x[i]
                elif unk.endswith(".Pt"):
                    ptot_val = x[i]

            return NetworkMixtureState(
                P=float(p_val), Pt=float(ptot_val), Tt=float(T_tot), m_dot=float(m), Y=Y
            )

        # For interior nodes, unpack from x
        indices = self._unknown_indices.get(node.id, [])
        unknowns = node.unknowns()

        # Default state
        state_dict = {
            "P": 101325.0,
            "Pt": 101325.0,
            "Tt": 300.0,
            "m_dot": 0.0,
            "Y": default_Y,
        }

        # Override with derived state if available
        if node.id in self._derived_states:
            T, Y, _ = self._derived_states[node.id]
            state_dict["Tt"] = T
            state_dict["Y"] = Y

        for i, unk in zip(indices, unknowns, strict=False):
            var_name = unk.split(".")[-1]
            if var_name in state_dict:
                state_dict[var_name] = x[i]
            elif var_name.startswith("Y[") and var_name.endswith("]"):
                idx_str = var_name[2:-1]
                idx = int(idx_str)
                # default_Y is shared, so we must copy it if we intend to mutate it.
                if state_dict["Y"] is default_Y:
                    state_dict["Y"] = list(default_Y)
                state_dict["Y"][idx] = x[i]

        # P1-4: Clamp fractions to [0, 1] before creating state to avoid NaNs in thermo conversion
        Y_raw = state_dict["Y"]
        Y_safe = [max(0.0, min(1.0, float(yi))) for yi in Y_raw]
        y_sum = sum(Y_safe)
        if y_sum > 0:
            Y_safe = [yi / y_sum for yi in Y_safe]

        # P5: Floor on static pressure to prevent NaN in thermo calls at near-zero P.
        P_safe = max(float(state_dict["P"]), 1000.0)

        return NetworkMixtureState(
            P=P_safe,
            Pt=float(state_dict["Pt"]),
            Tt=float(state_dict["Tt"]),
            m_dot=float(state_dict["m_dot"]),
            Y=Y_safe,
            incompressible=getattr(node, "regime", None) == "incompressible",
        )

    def _residuals_and_jacobian(
        self, x: np.ndarray, *, compute_jacobian: bool = True
    ) -> tuple[np.ndarray, Any]:
        """
        Evaluates the global residual vector and its sparse Jacobian.

        LAZY JACOBIAN (#489): the evaluation collects every element's and
        node's derivative dicts (they come with the residuals anyway); only
        the sparse ASSEMBLY -- name lookup, relay chaining, coo->csr -- is
        deferred. With compute_jacobian=False it is kept for the point just
        evaluated, and _jacobian_at_last(x) assembles it if the root finder
        then asks for J at that same point, so asking costs no second
        evaluation.
        """
        # Pure Pressure-Flow: Propagate states forward first. The relay
        # (state sensitivities) only for an immediate assembly; a deferred one
        # rebuilds it at the same x if J is asked for.
        self._eval_where = ""
        relay = self._propagate_states(x, with_relay=compute_jacobian)

        res = []
        # (first row, {eq_idx: {unknown: d}}, chains "Pt" relays)
        terms: list[tuple[int, dict, bool]] = []
        # (row, element, node id, unknown indices, sign) for the mass rows
        mass_terms: list[tuple[int, Any, str, list[int], float]] = []
        # (row, column, coefficient) for a FuelBoundary's own m_dot unknown (#892)
        fuel_terms: list[tuple[int, int, float]] = []

        # 1. Node Residuals
        for node_id, node in self.network.nodes.items():
            self._eval_where = f"node '{node_id}'"
            if isinstance(node, PressureBoundary):
                continue

            start_res_idx = len(res)

            # Node thermodynamic/stagnation constraints
            state = self._get_node_state(node, x)
            if isinstance(node, FuelBoundary):
                # The combustor's derived Tt, for the fuel controller's row (#892).
                node._target_Tt = float(
                    self._get_node_state(self.network.nodes[node.combustor], x).Tt
                )
            node_res, node_jac = node.residuals(state)
            res.extend(node_res)
            terms.append((start_res_idx, node_jac, False))

            # Mass Conservation: Sum(m_dot_in) - Sum(m_dot_out) = 0
            #
            # Skip for port-MCNs of a MultiPortChamberBase: the junction's
            # sum-mdot residual replaces the per-port mass row. Without this
            # skip the row degenerates to 0=0 (both the channel and the junction
            # report the same flow at the node with opposite signs, producing a
            # zero row in the Jacobian).
            if getattr(node, "_is_junction_port", False):
                continue
            mass_res_idx = len(res)
            upstream_elems = self.network.get_upstream_elements(node_id)
            downstream_elems = self.network.get_downstream_elements(node_id)

            m_dot_in = 0.0
            m_dot_out = 0.0

            if isinstance(node, MassFlowBoundary):
                # If connected as a sink (has upstream, no downstream), subtract from balance.
                # Otherwise (source or injection) add to balance.
                sink = bool(upstream_elems and not downstream_elems)
                if sink:
                    m_dot_out += getattr(node, "m_dot", 0.0)
                else:
                    m_dot_in += getattr(node, "m_dot", 0.0)
                if isinstance(node, FuelBoundary):
                    # d(m_dot)/d(ln m_dot) = m_dot
                    fuel_terms.append(
                        (
                            mass_res_idx,
                            self._name_to_index[node.m_dot_unknown],
                            (-1.0 if sink else 1.0) * float(node.m_dot),
                        )
                    )

            for elem in upstream_elems:
                indices = self._unknown_indices.get(elem.id)
                if indices:
                    m_dot_in += elem.flow_at_node(node_id, x, indices)
                    mass_terms.append((mass_res_idx, elem, node_id, indices, 1.0))
                else:
                    from_node = self.network.nodes[elem.from_node]
                    m_dot_in += getattr(from_node, "m_dot", 0.0)
                    if isinstance(from_node, FuelBoundary):
                        fuel_terms.append(
                            (
                                mass_res_idx,
                                self._name_to_index[from_node.m_dot_unknown],
                                float(from_node.m_dot),
                            )
                        )

            for elem in downstream_elems:
                indices = self._unknown_indices.get(elem.id)
                if indices:
                    m_dot_out += elem.flow_at_node(node_id, x, indices)
                    mass_terms.append((mass_res_idx, elem, node_id, indices, -1.0))
                else:
                    to_node = self.network.nodes[elem.to_node]
                    m_dot_out += getattr(to_node, "m_dot", 0.0)

            res.append(m_dot_in - m_dot_out)

        # 2. Element Residuals
        from .components import MultiPortChamberBase, TeeJunctionElement

        for elem_id, element in self.network.elements.items():
            self._eval_where = f"element '{elem_id}'"
            start_res_idx = len(res)

            m_indices = self._unknown_indices.get(elem_id)

            if isinstance(element, TeeJunctionElement):
                nodes = self.network.nodes
                state_com = self._get_node_state(nodes[element.common_node], x)
                state_straight = self._get_node_state(nodes[element.straight_node], x)
                state_branch = self._get_node_state(nodes[element.branch_node], x)
                if m_indices:
                    state_com.m_dot = float(x[m_indices[0]])
                    state_branch.m_dot = float(x[m_indices[1]])
                    state_straight.m_dot = float(x[m_indices[0]]) - float(x[m_indices[1]])
                elem_res, elem_jac = element.residuals(state_com, state_straight, state_branch)
            elif isinstance(element, MultiPortChamberBase):
                # N-port momentum-CV junction: one P_jct unknown plus N+1 residuals.
                # Port mdots are sourced from each connecting element's m_dot
                # unknown, sign-mapped to junction convention (positive = out).
                nodes = self.network.nodes
                port_states = [self._get_node_state(nodes[pid], x) for pid in element.port_nodes]
                P_jct_val = float(x[m_indices[0]]) if m_indices else 0.0
                port_mdots: list[float] = []
                for i, outer_id in enumerate(element._port_element_ids):
                    outer_indices = self._unknown_indices.get(outer_id, [])
                    outer_mdot = float(x[outer_indices[0]]) if outer_indices else 0.0
                    port_mdots.append(element._port_signs[i] * outer_mdot)
                elem_res, elem_jac = element.residuals(port_states, P_jct_val, port_mdots)
            else:
                node_in = self.network.nodes[element.from_node]
                node_out = self.network.nodes[element.to_node]

                state_in = self._get_node_state(node_in, x)
                state_out = self._get_node_state(node_out, x)

                # Unpack element-specific unknowns (like m_dot) into the state
                if m_indices:
                    state_in.m_dot = x[m_indices[0]]
                    state_out.m_dot = x[m_indices[0]]

                face = self._main_face(element, node_out, state_out)

                flows = self._network_flows(element, x)
                if flows is None:
                    elem_res, elem_jac = element.residuals(state_in, state_out)
                else:
                    elem_res, elem_jac = element.residuals(state_in, state_out, flows=flows)
                if face is not None:
                    self._chain_main_face(elem_jac, node_out.id, face)
            res.extend(elem_res)
            terms.append((start_res_idx, elem_jac, True))

        n_res = len(res)

        def assemble() -> Any:
            nonlocal relay
            if not compute_jacobian:
                # Same x, same states: only the sensitivities are new.
                relay = self._propagate_states(x)
            # Before any element Jacobian is read (LazyJacobians compute on
            # first read): for each node whose composition moves with x, the
            # directions it moves in -- where a costly d/dY is worth computing,
            # and the only directions it is needed along.
            directions = self._composition_directions(relay, x)
            for element in self.network.elements.values():
                element._composition_directions = directions
            rows: list[int] = []
            cols: list[int] = []
            data: list[float] = []
            name_to_index = self._name_to_index
            for start, jac_dict, chain_pt in terms:
                for eq_idx, var_derivs in jac_dict.items():
                    row = start + eq_idx
                    for unk_name, deriv in var_derivs.items():
                        if unk_name in name_to_index:
                            rows.append(row)
                            cols.append(name_to_index[unk_name])
                            data.append(deriv)
                        elif "." in unk_name:
                            # Relay for derived properties (T, Y, Pt_mix; and
                            # Pt for element residuals)
                            parts = unk_name.split(".")
                            nid, prop = parts[0], parts[1]
                            if prop == "T":
                                # The 0.7 name of the stagnation-T column
                                # (#782); relayed with a warning in 0.8.
                                raise KeyError(
                                    f"Jacobian key {unk_name!r}: the stagnation-T "
                                    f"column is {nid + '.Tt'!r} since 0.8; the "
                                    "'.T' spelling was removed in 0.9.0"
                                )
                            if nid in relay:
                                for unk_idx, sens_pkg in relay[nid].items():
                                    if prop == "Tt" and "T" in sens_pkg:
                                        rows.append(row)
                                        cols.append(unk_idx)
                                        data.append(deriv * sens_pkg["T"])
                                    elif prop.startswith("Y") and "Y" in sens_pkg:
                                        # Handle species composition relay
                                        if "[" in prop:
                                            s_idx = int(prop.split("[")[1].replace("]", ""))
                                            rows.append(row)
                                            cols.append(unk_idx)
                                            data.append(deriv * sens_pkg["Y"][s_idx])
                                    elif prop == "Pt_mix" and "Pt_mix" in sens_pkg:
                                        rows.append(row)
                                        cols.append(unk_idx)
                                        data.append(deriv * sens_pkg["Pt_mix"])
                                    elif chain_pt and prop == "Pt" and "Pt" in sens_pkg:
                                        rows.append(row)
                                        cols.append(unk_idx)
                                        data.append(deriv * sens_pkg["Pt"])
            for row, elem, node_id, indices, sign in mass_terms:
                for col, coeff in elem.flow_jac_at_node(node_id, indices).items():
                    rows.append(row)
                    cols.append(col)
                    data.append(sign * coeff)
            for row, col, coeff in fuel_terms:
                rows.append(row)
                cols.append(col)
                data.append(coeff)
            return sp.coo_matrix((data, (rows, cols)), shape=(n_res, len(x))).tocsr()

        if compute_jacobian:
            self._deferred_jacobian = None
            return self._offset(np.array(res, dtype=float)), assemble()
        self._deferred_jacobian = (x.copy(), assemble)
        return self._offset(np.array(res, dtype=float)), None

    def _offset(self, res: np.ndarray) -> np.ndarray:
        """The residual less the Newton homotopy's offset (1 - lam) F(x0), if one
        is active (#642); the Jacobian is the same."""
        off = getattr(self, "_homotopy_offset", None)
        return res if off is None else res - off

    @staticmethod
    def _composition_directions(relay: dict[str, Any], x: np.ndarray) -> dict[str, list[Any]]:
        """Per node whose derived composition moves with the unknowns, an
        orthonormal basis of the directions it moves in.

        A node's dY/dx_j, over all unknowns j, span at most (inflow
        compositions - 1) directions: ONE for a plenum mixing two gases. An
        element differentiating along that basis and projecting back,
        g = sum_i (dG/dv_i) v_i, has g . dY/dx_j = dG/dY . dY/dx_j exactly for
        every j, at one directional derivative instead of one per species.

        A node fed by one gas has dY/dx = 0 up to rounding (~1e-17/m); the
        columns are weighted by |x_j| (the change in Y over a 100% change of
        x_j) and directions below 1e-10 of that, or 1e-9 of the largest, are
        dropped.
        """
        directions: dict[str, list[Any]] = {}
        for nid, cols in relay.items():
            vecs = [
                pkg["Y"] * max(abs(x[idx]), 1e-3)
                for idx, pkg in cols.items()
                if pkg.get("Y") is not None
            ]
            if not vecs:
                continue
            M = np.column_stack(vecs)
            if float(np.max(np.abs(M))) <= 1e-10:
                continue
            U, S, _ = np.linalg.svd(M, full_matrices=False)
            keep = max(1e-10, 1e-9 * S[0]) < S
            if np.any(keep):
                directions[nid] = [U[:, i].copy() for i in np.flatnonzero(keep)]
        return directions

    def _jacobian_at_last(self, x: np.ndarray) -> Any:
        """The deferred Jacobian of the last residual-only evaluation, if it
        was at exactly x; else None (the caller evaluates afresh)."""
        pending = getattr(self, "_deferred_jacobian", None)
        if pending is None or not np.array_equal(pending[0], x):
            return None
        return pending[1]()

    def _residuals(self, x: np.ndarray) -> np.ndarray:
        """
        Evaluates the global residual vector given the current unknown vector x.
        """
        res, _ = self._residuals_and_jacobian(x, compute_jacobian=False)
        return res

    def _strict(self) -> bool:
        """Whether the starting-point safety net re-raises (see SOLVER_STRICT_ENV)."""
        if self.raise_internal_errors is not None:
            return bool(self.raise_internal_errors)
        return os.environ.get(SOLVER_STRICT_ENV, "") not in ("", "0")

    def _record_exception(self, exc: BaseException, phase: str | None = None) -> None:
        """Keep a caught exception for ``__caught_exceptions__``: the phase
        (``start``: the attempt's first point; ``probe``: a trial point the
        root finder then rejected), the node or element being evaluated, the
        type, the message and the innermost frames."""
        self._caught_count += 1
        if len(self._caught_exceptions) >= _MAX_CAUGHT_RECORDS:
            return
        frames = traceback.extract_tb(exc.__traceback__)[-3:]
        self._caught_exceptions.append(
            {
                "phase": phase or self._eval_phase,
                "where": self._eval_where,
                "type": type(exc).__name__,
                "message": str(exc)[:500],
                "frames": [f"{f.filename.rsplit('/', 1)[-1]}:{f.lineno} {f.name}" for f in frames],
            }
        )

    def _record_diagnostics_exception(self, exc: BaseException, where: str) -> None:
        """Keep a diagnostics exception on a FAILED solve's final state
        (phase ``diagnostics``, ``where`` the node or element): that state
        reports no diagnostics for it instead of raising out of solve (#810)."""
        saved = self._eval_where
        self._eval_where = where
        try:
            self._record_exception(exc, phase="diagnostics")
        finally:
            self._eval_where = saved

    def _failed_start_result(self, x0_use: np.ndarray, exc: BaseException) -> dict[str, Any]:
        """The result of an attempt whose starting point could not be
        evaluated: not converged, outcome ERROR, nothing derived. The retries
        in ``solve`` then try their other starts."""
        where = f" at {self._eval_where}" if self._eval_where else ""
        return {
            "__success__": False,
            "__message__": (
                f"The starting point could not be evaluated{where}: "
                f"{type(exc).__name__}: {str(exc)[:300]}"
            ),
            "__iterations__": 0,
            "__converged__": False,
            "__consistent__": None,
            "__inconsistent_elements__": [],
            "__P_nonpositive__": [],
            "__T_clamped__": [],
            "__heat_faded__": [],
            "__supersonic_chambers__": [],
            "__choked_bounded_junctions__": [],
            "__outcome__": SolveOutcome.ERROR,
            "__final_norm__": float("inf"),
            "__worst_residuals__": [],
            "__unknown_names__": list(self.unknown_names),
            "__x_solution__": list(x0_use),
        }

    def _with_caught(self, sol: dict[str, Any]) -> dict[str, Any]:
        """``sol`` with this solve's caught-exception record attached; in
        strict mode, after checking its row kinds (``check_row_kinds``)."""
        self._assert_row_kinds(sol)
        sol["__caught_exceptions__"] = list(self._caught_exceptions)
        sol["__caught_exception_count__"] = self._caught_count
        return sol

    def solve(
        self,
        method: str = "hybr",
        timeout: float | None = None,
        options: dict[str, Any] | None = None,
        use_jac: bool = True,
        x0: np.ndarray | None = None,
        init_strategy: Literal[
            "default",
            "outletref_warmstart",
            "homotopy",
            "analytical_pt_prop",
        ] = "default",
        lambda_steps: list[float] | None = None,
        auto_retry: bool = True,
        stall_handover: bool = True,
    ) -> dict[str, float]:
        """
        Executes the root finding algorithm to solve the network.
        Returns a dictionary mapping unknown names to their solved values.

        For cold solves (``x0=None``, ``init_strategy`` 'default' or
        'analytical_pt_prop') on COMPRESSIBLE networks, a failed solve is
        automatically retried once from an outlet-referenced incompressible warm
        start (disable via ``auto_retry=False``): the network is
        re-solved in the incompressible regime with element densities
        evaluated at the DOWNSTREAM static pressure, and that solution
        seeds the compressible retry. On blow-down networks this proxy
        tracks the compressible solution 6-7x closer than the classic
        inlet-referenced incompressible solve and rescues the cold
        path's non-monotone 'slow pockets' (5-tee GUI network at
        0.32 kg/s: cold needs 510 evals/161 s, the seeded retry 96
        evals/30 s; the seed even converges cases where the
        inlet-referenced seed fails -- 2026-07-05 outlet-ref seed
        experiment, tmp/outlet_ref_seed_experiment.py). When a
        ``timeout`` is given and the network holds a
        ``MultiPortChamberBase``, the primary attempt gets 40% of the
        budget, the proxy at most 25%, and the seeded retry the remainder
        (at least 30%). Compressible tee networks are best given 90-120 s
        total. Any other network keeps the whole budget for its cold
        attempt and retries only in the time left, the proxy taking at most
        half of it. Until #585 the retry was limited to junction networks; a
        duct-fed effusion liner, whose cold compressible solve stalls, is
        seeded to its root in one pass from the incompressible one.

        The mirror case: a network INCOMPRESSIBLE throughout (every
        two-regime element incompressible) whose cold solve fails is retried
        from a compressible warm start -- the same network solved with those
        elements compressible, whose root it shares at low Mach (#588). Its
        mass-flow residual forms converge from cold starts where the
        pressure-drop forms stall (a jet pump behind an incompressible pipe:
        zero entrainment from every cold seed, the physical root in a few
        steps from the compressible one). Same budget rules; the strategy is
        reported as ``compressible_warmstart``.

        Pass ``init_strategy='outletref_warmstart'`` to run the proxy
        seed directly as the primary initialization (no cold attempt
        first) -- useful on networks known to stall the cold path.
        Certified audit (2026-07): 16/16 same-root on branch
        topologies, 11/16 on merge -- prefer 'default' /
        'analytical_pt_prop' for merge networks.

        A cold solve of a network with ``ThermalWall`` s that still fails
        is retried once more from the same network solved WITHOUT its walls
        (``thermal_coupling_enabled = False``): the walls move a lot of heat
        into small flows, and a cold Newton step that reverses one of them
        leaves a node heated with almost nothing to carry the heat. The
        wall-less flows are a start on the right side of that (#471:
        impingement array under a hot duct, 6/9 cold vs 9/9 seeded).

        A cold solve of a network with channels carrying their wall heat in
        the momentum (Fanno-Rayleigh, #551) runs in two stages: those channels
        adiabatic first (this whole chain, 60% of the budget), then the exact
        problem from that root; the plain cold chain runs only if that fails.
        Cold, the guess's node temperatures put extreme heating or cooling into
        their first iterates (an impingement crossflow chain: 23 s cold, 4.7 s
        staged).

        A cold solve of a network with a jet plate under CR-3630's Table 4.1
        Cd (``cd_crossflow='florschuetz_metzger_su_1982'``) also runs in two
        stages: every crossflow-Cd plate at its own constant Cd first (60% of
        the budget), then the exact problem from that root; the plain cold
        chain only if that fails (#1021: B(5,4,1) at mc/mj 0.99, 25 s through
        the retries, staged a few s). Under the default Rogers-Hersh law the
        constant-Cd stage stays a retry after a failed cold chain.

        A numerical exception while evaluating an attempt's STARTING point
        (a property or closure refusing a physical-but-extreme state) does not
        abort the solve: that attempt ends with outcome ERROR and the retries
        above try their other starts. Exceptions at trial points the root
        finder then rejects were already caught. Every caught exception is
        returned in ``__caught_exceptions__`` (phase, the node or element
        being evaluated, type, message, innermost frames; at most 50, all
        counted in ``__caught_exception_count__``). Set
        ``raise_internal_errors = True`` on the solver, or the environment
        variable ``COMBAERO_SOLVER_STRICT=1``, to make the starting point raise
        as it used to (the test suite runs strict). Code bugs (TypeError,
        AttributeError, KeyError, IndexError) always raise. A strict solve
        also checks each converged solution's residual rows against their
        declared units and raises ``RowKindMismatchError`` on a misfiled one
        (``check_row_kinds``; RESIDUAL-ROW-SCALING).

        ``stall_handover`` (default True): a COLD attempt whose best |F|
        stalls (``_stall_detected``) hands the rest of its budget to the next
        rung -- hybr to the LM fallback, the LM fallback to the warm-start
        retry. False turns every stall handover off for this call, its
        retries and nested solves included: each cold rung then runs to its
        own share of the budget (#864).

        The fraction-loss check (#1004, ``fraction_unanchored``,
        ``__fraction_unanchored__``) runs once on the result of the outermost
        call, after every retry has restored the boundaries.
        """
        depth = getattr(self, "_solve_depth", 0)
        self._solve_depth = depth + 1
        try:
            # A fuel controller's target proven unreachable (#1114): refused
            # without the Newton search and its retries (up to 8 s).
            refused = self._fuel_refused_before_solve() if depth == 0 else None
            sol = refused or self._solve_entry(
                method=method,
                timeout=timeout,
                options=options,
                use_jac=use_jac,
                x0=x0,
                init_strategy=init_strategy,
                lambda_steps=lambda_steps,
                auto_retry=auto_retry,
                stall_handover=stall_handover,
            )
        finally:
            self._solve_depth = depth
        if depth == 0:
            sol = self._mark_fraction_unanchored(sol)
        return sol

    def _solve_entry(
        self,
        method: str = "hybr",
        timeout: float | None = None,
        options: dict[str, Any] | None = None,
        use_jac: bool = True,
        x0: np.ndarray | None = None,
        init_strategy: Literal[
            "default",
            "outletref_warmstart",
            "homotopy",
            "analytical_pt_prop",
        ] = "default",
        lambda_steps: list[float] | None = None,
        auto_retry: bool = True,
        stall_handover: bool = True,
    ) -> dict[str, float]:
        """``solve`` without the outermost-call checks; recursive calls go
        through ``solve``."""
        if not getattr(self, "_in_reference_swap", False):
            # Flow splitters whose reference the solve finds (#903): an outer
            # loop around this whole solve, see solve_with_reference_swaps.
            from .flow_splitter import auto_reference_groups, solve_with_reference_swaps

            if auto_reference_groups(self.network):
                self._in_reference_swap = True
                try:
                    return solve_with_reference_swaps(
                        self,
                        {
                            "method": method,
                            "timeout": timeout,
                            "options": options,
                            "use_jac": use_jac,
                            "x0": x0,
                            "init_strategy": init_strategy,
                            "lambda_steps": lambda_steps,
                            "auto_retry": auto_retry,
                            "stall_handover": stall_handover,
                        },
                    )
                finally:
                    self._in_reference_swap = False
        if not stall_handover and self._stall_handover:
            self._stall_handover = False
            try:
                return self.solve(
                    method=method,
                    timeout=timeout,
                    options=options,
                    use_jac=use_jac,
                    x0=x0,
                    init_strategy=init_strategy,
                    lambda_steps=lambda_steps,
                    auto_retry=auto_retry,
                )
            finally:
                self._stall_handover = True
        if (
            (self.homotopy_retry or self.bc_ramp_retry)
            and x0 is None
            and auto_retry
            and not getattr(self, "_in_homotopy", False)
            and not getattr(self, "_in_bc_ramp", False)
        ):
            _h0 = time.time()
            self._in_homotopy = True
            try:
                sol = self.solve(
                    method=method,
                    timeout=timeout,
                    options=options,
                    use_jac=use_jac,
                    init_strategy=init_strategy,
                    lambda_steps=lambda_steps,
                )
            finally:
                self._in_homotopy = False
            # Only a solve that did not converge: a root the checks rejected
            # (an unphysical junction root, an open ledger) is a different
            # failure, and a homotopy from the same guess tends to find it again.
            if sol.get("__success__", False) or sol.get("__converged__", False):
                return sol
            if self.homotopy_retry:
                # With the ramp to follow, an uncapped homotopy takes half of
                # the time left, not all of it.
                cap = self.homotopy_budget
                if cap is None and self.bc_ramp_retry and timeout is not None:
                    self.homotopy_budget = 0.5 * max(timeout - (time.time() - _h0), 0.0)
                elif cap is not None and self.bc_ramp_retry:
                    self.homotopy_budget = min(cap, NEWTON_HOMOTOPY_BUDGET_BEFORE_RAMP)
                try:
                    sol = self._newton_homotopy_retry(sol, method, timeout, options, use_jac, _h0)
                finally:
                    self.homotopy_budget = cap
                if sol.get("__success__", False):
                    return sol
            if self.bc_ramp_retry:
                return self._bc_ramp_retry(sol, method, timeout, options, use_jac, _h0)
            return sol
        kw: dict[str, Any] = {
            "method": method,
            "timeout": timeout,
            "options": options,
            "use_jac": use_jac,
            "init_strategy": init_strategy,
            "lambda_steps": lambda_steps,
        }
        if x0 is None and auto_retry and not getattr(self, "_in_rayleigh_stage", False):
            heated = self._rayleigh_heated_elements()
            if heated:
                return self._solve_rayleigh_staged(heated, kw)
        staged_first: list[Any] = []
        if x0 is None and auto_retry and not getattr(self, "_in_fuel_fixed_stage", False):
            # A fuel controller whose oxidizer comes through a pressure
            # boundary (#892): its start is sized on a guessed air flow, which
            # can be decades off, so the fuel-fixed stage runs FIRST and the
            # plain cold chain only if it fails -- after a failed cold attempt
            # the C3X rig network (GUI budget 30 s) had 0.2 s left.
            self.network.resolve_all_topology()
            guessed = [
                n
                for n in self.network.nodes.values()
                if isinstance(n, FuelBoundary)
                and n.m_dot_guess is None
                and any(isinstance(b, PressureBoundary) for b in self._fuel_oxidizer_sources(n))
            ]
            if guessed:
                staged_first = guessed
                _t0 = time.time()
                staged = self._solve_fuel_fixed_staged(
                    guessed, {"__message__": ""}, method, timeout, options, use_jac, _t0
                )
                if staged is not None:
                    return self._with_caught(staged)
                if timeout is not None:
                    kw["timeout"] = max(timeout - (time.time() - _t0), 1.0)
                    timeout = kw["timeout"]
        cd_first = False
        if x0 is None and auto_retry and not getattr(self, "_in_constant_cd_stage", False):
            # A plate under CR-3630's Table 4.1 Cd (#1021): its stagnant rows'
            # flow goes as dP^1.07, with no slope at dP = 0 to pull a backwards
            # cold seed round, and the root's basin is narrow (B(5,4,1) at
            # mc/mj 0.99: 1 of 3 starts at the root +/-1% converge). The
            # constant-Cd root is a start inside it, so that stage runs FIRST
            # and the plain cold chain only if it fails. Rogers-Hersh solves
            # cold and keeps its path.
            plates = self._crossflow_cd_plates()
            if any(el.cd_crossflow == "florschuetz_metzger_su_1982" for el in plates):
                cd_first = True
                _t0 = time.time()
                staged = self._solve_constant_cd_staged(
                    plates,
                    {"__message__": ""},
                    method,
                    None if timeout is None else 0.6 * timeout,
                    options,
                    use_jac,
                    _t0,
                )
                if staged is not None:
                    return self._with_caught(staged)
                if timeout is not None:
                    kw["timeout"] = max(timeout - (time.time() - _t0), 1.0)
                    timeout = kw["timeout"]
        self._caught_exceptions = []
        self._caught_count = 0
        _t0 = time.time()
        sol = self._solve_flow_retries(
            method=method,
            timeout=timeout,
            options=options,
            use_jac=use_jac,
            x0=x0,
            init_strategy=init_strategy,
            lambda_steps=lambda_steps,
            auto_retry=auto_retry,
        )
        if sol.get("__success__", False) or not auto_retry or x0 is not None:
            return self._with_caught(sol)
        if not getattr(self, "_in_fuel_fixed_stage", False):
            fuels = [
                n
                for n in self.network.nodes.values()
                if isinstance(n, FuelBoundary) and n not in staged_first
            ]
            if fuels:
                staged = self._solve_fuel_fixed_staged(
                    fuels, sol, method, timeout, options, use_jac, _t0
                )
                if staged is not None:
                    return self._with_caught(staged)
        if not cd_first and not getattr(self, "_in_constant_cd_stage", False):
            plates = self._crossflow_cd_plates()
            if plates:
                staged = self._solve_constant_cd_staged(
                    plates, sol, method, timeout, options, use_jac, _t0
                )
                if staged is not None:
                    return self._with_caught(staged)
        # Zero-flow fallback: a network whose only solution is zero flow -- equal
        # boundary pressures, nothing imposing a flow -- stalls because dP/dm -> 0
        # at m=0, so the default non-zero flow seed cannot descend to it. Seeding
        # the flows at zero lands on it directly. Gated to networks with no
        # pressure driving potential so it adds no attempt to ordinary failed
        # solves, and accepted only if it genuinely converges, so it never
        # invents a zero-flow answer for a network that actually carries flow.
        if self._no_pressure_driving_potential():
            zsol = self._zero_flow_retry(
                method=method,
                timeout=timeout,
                t0=_t0,
                options=options,
                use_jac=use_jac,
                primary=sol,
            )
            if zsol is not None:
                return self._with_caught(zsol)
        if not (self.network.thermal_coupling_enabled and self.network.walls):
            return self._with_caught(sol)
        remaining = None
        if timeout is not None:
            remaining = timeout - (time.time() - _t0)
            if remaining <= 0.0:
                return self._with_caught(sol)
        seed = self._wall_free_seed(method=method, timeout=remaining, use_jac=use_jac)
        if seed is None:
            return self._with_caught(sol)
        _primary_diag = getattr(self, "_diagnostic_data", None)
        if timeout is not None:
            remaining = max(timeout - (time.time() - _t0), 0.0)
        retry = self._solve_impl(
            method=method, timeout=remaining, options=options, use_jac=use_jac, x0=seed
        )
        primary_msg = str(sol.get("__message__", "")).strip()
        if retry.get("__success__", False):
            retry["__message__"] = (
                "Converged from the wall-free flow solution after the cold "
                f"attempt failed ({primary_msg[:120]})."
            )
            return self._with_caught(retry)
        r_norm, p_norm = retry.get("__final_norm__"), sol.get("__final_norm__")
        if r_norm is not None and p_norm is not None and float(r_norm) < float(p_norm):
            return self._with_caught(retry)
        if _primary_diag is not None:
            self._diagnostic_data = _primary_diag
        sol["__message__"] = f"{primary_msg} [wall-free warm-start retry also failed]"
        return self._with_caught(sol)

    def _rayleigh_heated_elements(self) -> list[Any]:
        """Channels whose cold solve is staged (#551): those that carry their
        wall heat into a node of their own. A duct whose exit state the solver
        derives (#672) is not: staging gives the adiabatic stage 60% of the
        budget, and a stress draw whose cold solve needed 17 s then lost it.
        None with thermal coupling off: no wall moves heat (the wall-free seed
        of a staged solve runs with it off)."""
        self.network.resolve_all_topology()
        if not self.network.thermal_coupling_enabled:
            return []
        return [el for el in self.network.elements.values() if getattr(el, "_rayleigh_ends", None)]

    def _crossflow_cd_plates(self) -> list[Any]:
        """Jet plates whose hole Cd follows the crossflow (#800)."""
        from .components import ImpingementPlateElement

        return [
            el
            for el in self.network.elements.values()
            if isinstance(el, ImpingementPlateElement) and el.cd_crossflow != "none"
        ]

    def _solve_fuel_fixed_staged(
        self,
        fuels: list["FuelBoundary"],
        failed: dict[str, Any],
        method: str,
        timeout: float | None,
        options: dict[str, Any] | None,
        use_jac: bool,
        t0: float,
    ) -> dict[str, Any] | None:
        """Retry of a failed cold solve with fuel controllers (#892): first
        with each fuel flow held at its start (the full retry chain), then the
        fuel re-sized on the oxidizer that root carries, and the controllers
        solved from there. None if either stage fails.

        A fuel start is sized on the guessed oxidizer flow, which can be
        decades off: the C3X rig's burner carries 27 kg/s against the generic
        0.01 kg/s reference, so the controller started 2700x lean and the cold
        solve stalled at |F| 5e4, where the same network with its fuel fixed
        solves.
        """

        def left() -> float | None:
            return None if timeout is None else max(timeout - (time.time() - t0), 1.0)

        if timeout is not None and timeout - (time.time() - t0) <= 1.0:
            return None
        x_start = self._build_x0()
        idx = {n: i for i, n in enumerate(self.unknown_names)}
        self._in_fuel_fixed_stage = True
        was = [f._frozen_ln_m for f in fuels]
        try:
            for f in fuels:
                f._frozen_ln_m = float(x_start[idx[f.m_dot_unknown]])
            seed = self.solve(method=method, timeout=left(), options=options, use_jac=use_jac)
        finally:
            for f, v in zip(fuels, was, strict=True):
                f._frozen_ln_m = v
            self._in_fuel_fixed_stage = False
        if not seed.get("__success__", False) or seed.get("__x_solution__") is None:
            return None
        x = np.asarray(seed["__x_solution__"], dtype=float).copy()
        idx = {n: i for i, n in enumerate(self.unknown_names)}
        # The node states at the seed, for the oxidizer each fuel burns with.
        self._residuals_and_jacobian(x, compute_jacobian=False)
        for f in fuels:
            streams, Y_f = self._fuel_inflows(f, x, seed)
            if not streams:
                return None
            m_ox, T_ox, Y_ox = streams[0]
            try:
                m = max(self._fuel_flow_for(f, m_ox, T_ox, Y_ox, Y_f), 1e-9 * m_ox)
            except (ValueError, RuntimeError):
                return None  # no root on its branch: the primary's reason stands
            x[idx[f.m_dot_unknown]] = math.log(m)
            for eid in f._fuel_path_elements(self.network):
                k = idx.get(f"{eid}.m_dot")
                if k is not None:
                    x[k] = m
        sol = self.solve(method=method, timeout=left(), options=options, use_jac=use_jac, x0=x)
        if not sol.get("__success__", False):
            return None
        primary_msg = str(failed.get("__message__", "")).strip()
        why = f"; the primary attempt failed with: {primary_msg}" if primary_msg else ""
        sol["__message__"] = (
            f"{str(sol.get('__message__', '')).strip()} [from the fuel-fixed solution: "
            "the flows solved with each fuel flow held at its start, then the fuel "
            f"controllers{why}]"
        )
        return sol

    def _solve_constant_cd_staged(
        self,
        plates: list[Any],
        failed: dict[str, Any],
        method: str,
        timeout: float | None,
        options: dict[str, Any] | None,
        use_jac: bool,
        t0: float,
    ) -> dict[str, Any] | None:
        """Retry of a failed cold solve with crossflow-Cd jet plates (#800):
        first with every plate at its own constant Cd (the full retry chain),
        then the exact problem from that root. None if either stage fails.

        A retry, not a stage of every solve: a network that converges cold
        keeps its path. CR-3630 D.2 B(5,4,1) at mc/mj 1.01 with a 20 K
        heated bypass stalled cold at |F| 1.8e-2 (Jacobian exact there to
        3e-7); from the constant-Cd root, whose flows differ only at the
        crossflow-dominated rows, it converged in 1.1 s.

        The exception is a network with a plate under Table 4.1
        (``'florschuetz_metzger_su_1982'``, #1021): ``solve`` runs this stage
        first there, with 60% of the budget, and the cold chain only if it
        fails. ``failed`` is then empty.
        """

        def left() -> float | None:
            return None if timeout is None else max(timeout - (time.time() - t0), 1.0)

        if timeout is not None and timeout - (time.time() - t0) <= 1.0:
            return None
        self._in_constant_cd_stage = True
        was = [el._cd_crossflow_suspended for el in plates]
        try:
            for el in plates:
                el._cd_crossflow_suspended = True
            seed = self.solve(method=method, timeout=left(), options=options, use_jac=use_jac)
        finally:
            for el, flag in zip(plates, was, strict=True):
                el._cd_crossflow_suspended = flag
            self._in_constant_cd_stage = False
        if not seed.get("__success__", False) or seed.get("__x_solution__") is None:
            return None
        sol = self.solve(
            method=method,
            timeout=left(),
            options=options,
            use_jac=use_jac,
            x0=np.asarray(seed["__x_solution__"]),
        )
        if not sol.get("__success__", False):
            return None
        primary_msg = str(failed.get("__message__", "")).strip()
        why = (
            f"the primary attempt failed with: {primary_msg}"
            if primary_msg
            else "run first for the Table 4.1 Cd"
        )
        sol["__message__"] = (
            f"{str(sol.get('__message__', '')).strip()} [from the constant-Cd solution: "
            f"jet plates solved with their own Cd, then with the crossflow Cd; {why}]"
        )
        return sol

    def _solve_rayleigh_staged(self, heated: list[Any], kw: dict[str, Any]) -> dict[str, Any]:
        """Cold solve of a network with Fanno-Rayleigh channels, in two stages
        (#551): first with their heat left out of the momentum (adiabatic
        Fanno, the full retry chain), then the exact problem from that root.

        Cold, a guess's node temperatures put extreme heating or cooling into
        every heated duct's first iterates, and an impingement crossflow chain
        stalled for 22 s; from the adiabatic root, whose flows differ by the
        O(M^2) Rayleigh change only, it converged in 0.7 s. If the exact stage
        fails, the plain cold chain runs in the time left.
        """
        t0 = time.time()
        timeout = kw.get("timeout")
        self._in_rayleigh_stage = True
        # Restore each duct's own flag, never clear it: a solve nested inside
        # stage 1 (the wall-free seed solves the same elements) used to switch
        # the heat back on under the outer stage (#672).
        was = [getattr(el, "_rayleigh_suspended", False) for el in heated]
        try:
            for el in heated:
                el._rayleigh_suspended = True
            seed = self.solve(**{**kw, "timeout": None if timeout is None else 0.6 * timeout})
        finally:
            for el, flag in zip(heated, was, strict=True):
                el._rayleigh_suspended = flag
            self._in_rayleigh_stage = False

        def left() -> float | None:
            return None if timeout is None else max(timeout - (time.time() - t0), 1.0)

        if seed.get("__success__", False) and seed.get("__x_solution__") is not None:
            sol = self.solve(**{**kw, "timeout": left()}, x0=np.asarray(seed["__x_solution__"]))
            if sol.get("__success__", False):
                sol["__message__"] = (
                    f"{str(sol.get('__message__', '')).strip()} [from the adiabatic-momentum "
                    "solution: heated channels solved without, then with, their heat in "
                    "the momentum]"
                )
                return sol
        self._in_rayleigh_stage = True
        try:
            return self.solve(**{**kw, "timeout": left()})
        finally:
            self._in_rayleigh_stage = False

    def _bc_ramp_retry(
        self,
        failed: dict[str, Any],
        method: str,
        timeout: float | None,
        options: dict[str, Any] | None,
        use_jac: bool,
        t0: float,
    ) -> dict[str, Any]:
        """Last retry of a failed cold solve (#642): continuation in the
        boundary conditions. At lam every PressureBoundary sits at
        Pm + lam (Pt - Pm), Pm the mean of their Pt, and every
        MassFlowBoundary m_dot and user EnergyBoundary Q / fraction at lam
        times its own. The gentle end (BC_RAMP_LAM0, halved down to
        BC_RAMP_LAM_MIN until it solves) is solved cold with the usual
        retries; lam then rises, each solve warm-started from the last root,
        the step doubling after a solve (at most 0.5) and halving after a
        failure. The lam = 1 solve is a plain warm solve with every physical
        check. The boundaries are restored whatever happens.

        Complementary to the Newton homotopy: that one starts from the cold
        guess, a state with flow; this one from a network with almost no
        drive, which is easy for some networks and is itself the hard
        near-zero-drive regime for others (15 of the 21 random draws it did
        not rescue have no gentle end that solves).
        """
        pbs = [n for n in self.network.nodes.values() if isinstance(n, PressureBoundary)]
        # A FuelBoundary's flow is an unknown (#892), not a drive to ramp.
        mfbs = [
            n
            for n in self.network.nodes.values()
            if isinstance(n, MassFlowBoundary) and not isinstance(n, FuelBoundary)
        ]
        wall_ebs = {id(eb) for eb in getattr(self, "_wall_ebs", {}).values()}
        ebs = [
            eb
            for n in self.network.nodes.values()
            for eb in (getattr(n, "energy_boundaries", None) or [])
            if id(eb) not in wall_ebs
        ]
        if not pbs and not mfbs:
            return failed
        pt0 = [float(n.Pt) for n in pbs]
        m0 = [float(n.m_dot) for n in mfbs]
        q0 = [(float(eb.Q), float(eb.fraction)) for eb in ebs]
        pm = float(np.mean(pt0)) if pt0 else 0.0

        def set_lam(lam: float) -> None:
            for n, pt in zip(pbs, pt0, strict=True):
                n.Pt = pt if lam >= 1.0 else pm + lam * (pt - pm)
            for n, m in zip(mfbs, m0, strict=True):
                n.m_dot = m if lam >= 1.0 else lam * m
            for eb, (q, f) in zip(ebs, q0, strict=True):
                eb.Q, eb.fraction = (q, f) if lam >= 1.0 else (lam * q, lam * f)

        def left() -> float | None:
            return None if timeout is None else timeout - (time.time() - t0)

        saved = {k: getattr(self, k) for k in ("_last_x", "_diagnostic_data") if hasattr(self, k)}
        x: np.ndarray | None = None
        lam_ok, lam, step, solves = 0.0, BC_RAMP_LAM0, BC_RAMP_LAM0, 0
        self._in_bc_ramp = True
        try:
            while solves < BC_RAMP_MAX_SOLVES:
                remaining = left()
                if remaining is not None and remaining < 0.5:
                    break
                set_lam(lam)
                if x is None:
                    # A cold start inherits nothing from the failed attempts:
                    # the Fanno kernels' warm-start guesses left by them kept
                    # a stress draw's gentle end from solving (random 34,
                    # rescued once they are dropped).
                    for el in self.network.elements.values():
                        for k in _KERNEL_WARM_STARTS:
                            if k in vars(el):
                                delattr(el, k)
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        r = (
                            self.solve(
                                method=method, timeout=remaining, options=options, use_jac=use_jac
                            )
                            if x is None
                            else self.solve(
                                method=method,
                                timeout=remaining,
                                options=options,
                                use_jac=use_jac,
                                x0=x,
                            )
                        )
                except _STATE_EXCEPTIONS:
                    r = {}
                solves += 1
                if r.get("__success__", False) and r.get("__x_solution__") is not None:
                    x, lam_ok = np.asarray(r["__x_solution__"], dtype=float), lam
                    if lam >= 1.0:
                        r["__message__"] = (
                            f"{str(r.get('__message__', '')).strip()} [boundary-condition "
                            f"ramp in {solves} solves after the cold solve and its retries "
                            f"failed: {str(failed.get('__message__', '')).strip()[:120]}]"
                        )
                        return r
                    step = min(2.0 * step, 0.5)
                    lam = min(1.0, lam + step)
                elif x is None:
                    lam *= 0.5  # an even gentler cold start
                    if lam < BC_RAMP_LAM_MIN:
                        break
                else:
                    step *= 0.5
                    if step < BC_RAMP_LAM_MIN:
                        break
                    lam = lam_ok + step
        finally:
            set_lam(1.0)
            self._in_bc_ramp = False
        for k, v in saved.items():
            setattr(self, k, v)
        return failed

    def _newton_homotopy_retry(
        self,
        failed: dict[str, Any],
        method: str,
        timeout: float | None,
        options: dict[str, Any] | None,
        use_jac: bool,
        t0: float,
    ) -> dict[str, Any]:
        """Last retry of a failed cold solve (#642): a Newton homotopy from the
        cold guess x0, H(x, lam) = F(x) - (1 - lam) F(x0).

        At lam = 0 the guess itself is the root -- a state WITH flow, unlike a
        ramp of the boundary pressures, which starts at zero drive and had no
        solution to start from in 13 of 23 failed random draws. lam rises
        adaptively (doubling after a solved step, halving after a failed one,
        down to 1/128), each step warm-started from the last; the lam = 1 step
        is a plain solve with every physical check. It rescued 9 of those 23
        draws (the pressure ramp 7). It cannot help where the Jacobian at x0
        is itself singular -- the solve then stalls at lam ~ 0, and it stops
        once no step from the guess has succeeded down to 1/16. It runs for at
        most ``homotopy_budget`` seconds within the caller's time left (default
        NEWTON_HOMOTOPY_BUDGET; None: all of the time left).
        """
        left0 = None if timeout is None else timeout - (time.time() - t0)
        if left0 is not None and left0 < 1.0:
            return failed
        cap = self.homotopy_budget
        if left0 is None:
            budget = NEWTON_HOMOTOPY_BUDGET if cap is None else cap
        else:
            budget = left0 if cap is None else min(left0, cap)
        h_t0 = time.time()
        try:
            x0 = self._build_x0()
            F0 = np.asarray(self._residuals(x0), dtype=float).copy()
        except _STATE_EXCEPTIONS:
            return failed
        if not np.all(np.isfinite(F0)):
            return failed
        x, lam_ok, step = x0, 0.0, 0.25
        # A failed homotopy leaves the solver as the cold solve left it: its
        # last iterate and diagnostics describe the attempt that is returned.
        saved = {k: getattr(self, k) for k in ("_last_x", "_diagnostic_data") if hasattr(self, k)}
        self._in_homotopy = True
        try:
            while step >= 1.0 / 128:
                left = budget - (time.time() - h_t0)
                # No step at all from the guess: its Jacobian is singular
                # there, which no homotopy can fix.
                if left < 0.5 or (lam_ok == 0.0 and step < 1.0 / 16):
                    break
                lam = min(1.0, lam_ok + step)
                self._homotopy_offset = None if lam >= 1.0 else (1.0 - lam) * F0
                try:
                    r = self._solve_flow_retries(
                        method=method,
                        timeout=left,
                        options=options,
                        use_jac=use_jac,
                        x0=x,
                        init_strategy="default",
                        lambda_steps=None,
                        auto_retry=False,
                    )
                except _STATE_EXCEPTIONS:
                    r = {}
                ok = bool(r.get("__success__") if lam >= 1.0 else r.get("__converged__"))
                if ok and r.get("__x_solution__") is not None:
                    x, lam_ok = np.asarray(r["__x_solution__"], dtype=float), lam
                    if lam >= 1.0:
                        r["__message__"] = (
                            f"{str(r.get('__message__', '')).strip()} [Newton homotopy from "
                            "the cold guess after the cold solve and its retries failed: "
                            f"{str(failed.get('__message__', '')).strip()[:120]}]"
                        )
                        return self._with_caught(r)
                    step = min(2.0 * step, 0.5)
                else:
                    step *= 0.5
        finally:
            self._homotopy_offset = None
            self._in_homotopy = False
        for k, v in saved.items():
            setattr(self, k, v)
        return failed

    def _wall_free_seed(
        self, method: str, timeout: float | None, use_jac: bool
    ) -> np.ndarray | None:
        """This network's solution with its walls switched off, ordered as
        this solver's unknowns; ``None`` when it fails or the unknowns differ."""
        net = self.network
        net.thermal_coupling_enabled = False
        try:
            seed_solver = NetworkSolver(net)
            seed_solver.raise_internal_errors = self.raise_internal_errors
            seed_solver._stall_handover = self._stall_handover
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r = seed_solver.solve(method=method, timeout=timeout, use_jac=use_jac)
        finally:
            net.thermal_coupling_enabled = True
        for rec in r.get("__caught_exceptions__", []):
            if len(self._caught_exceptions) < _MAX_CAUGHT_RECORDS:
                self._caught_exceptions.append({**rec, "phase": f"wall_free_seed:{rec['phase']}"})
        self._caught_count += int(r.get("__caught_exception_count__", 0))
        if not r.get("__success__", False):
            return None
        by_name = dict(zip(seed_solver.unknown_names, r["__x_solution__"], strict=True))
        if set(by_name) != set(self.unknown_names):
            return None
        return np.array([by_name[k] for k in self.unknown_names], dtype=float)

    def _no_pressure_driving_potential(self) -> bool:
        """True when the boundaries cannot drive flow by a pressure difference:
        every ``PressureBoundary`` is at the same total pressure and no
        ``MassFlowBoundary`` imposes a flow. A heat source can still drive flow,
        so the zero-flow retry this gates stays self-validating (accepted only
        on convergence); this only keeps that retry off ordinary failed solves
        that plainly have a driving potential."""
        pressures: list[float] = []
        for node in self.network.nodes.values():
            if isinstance(node, MassFlowBoundary) and abs(float(node.m_dot)) > 0.0:
                return False
            if isinstance(node, PressureBoundary):
                pressures.append(float(node.Pt))
        if not pressures:
            return False
        lo, hi = min(pressures), max(pressures)
        return (hi - lo) <= 1e-6 * max(abs(hi), abs(lo), 1.0)

    def _zero_flow_seed(self) -> np.ndarray | None:
        """The default initial guess with every flow unknown set to zero.

        ``None`` if the network has no flow unknowns (nothing to seed)."""
        x0 = np.array(self._build_x0(), dtype=float)
        touched = False
        for i, name in enumerate(self.unknown_names):
            if name.endswith(".m_dot"):
                x0[i] = 0.0
                touched = True
        return x0 if touched else None

    def _zero_flow_retry(
        self,
        method: str,
        timeout: float | None,
        t0: float,
        options: dict[str, Any] | None,
        use_jac: bool,
        primary: dict[str, float],
    ) -> dict[str, float] | None:
        """Retry a failed cold solve from an all-zero-flow seed.

        Returns the retry only if it genuinely converges (so a network that
        actually carries flow, where zero flow leaves a non-zero residual, is
        never reported as zero flow); otherwise ``None`` and the primary's
        diagnostics are restored. Cheap: on the zero-flow solution it lands in a
        step or two, and on any other network it fails fast."""
        seed = self._zero_flow_seed()
        if seed is None:
            return None
        remaining = None
        if timeout is not None:
            # A couple of seconds is plenty to confirm |F|=0 at the seed; grant
            # that floor even if the primary used the budget (the backend's hard
            # deadline has its own buffer on top).
            remaining = max(timeout - (time.time() - t0), 2.0)
        primary_diag = getattr(self, "_diagnostic_data", None)
        retry = self._solve_impl(
            method=method, timeout=remaining, options=options, use_jac=use_jac, x0=seed
        )
        if retry.get("__success__", False):
            primary_msg = str(primary.get("__message__", "")).strip()
            retry["__message__"] = (
                "Converged to zero flow after the initial attempt failed "
                f"({primary_msg[:120]}): the network has no net driving potential."
            )
            return retry
        if primary_diag is not None:
            self._diagnostic_data = primary_diag
        return None

    def _solve_flow_retries(
        self,
        method: str,
        timeout: float | None,
        options: dict[str, Any] | None,
        use_jac: bool,
        x0: np.ndarray | None,
        init_strategy: str,
        lambda_steps: list[float] | None,
        auto_retry: bool,
    ) -> dict[str, float]:
        """The cold solve with its warm-start auto-retry (see solve).

        A network with compressible elements retries from the incompressible
        proxy; one that is incompressible throughout retries from the
        compressible proxy (#588) -- the same root at low Mach, reached from
        the other side.
        """
        has_compressible = bool(self._compressible_element_overrides())
        all_incompressible = not has_compressible and bool(self._incompressible_element_ids())
        retry_applicable = (
            auto_retry
            and x0 is None
            and init_strategy in ("default", "analytical_pt_prop")
            and (has_compressible or all_incompressible)
        )
        if not retry_applicable:
            return self._solve_impl(
                method=method,
                timeout=timeout,
                options=options,
                use_jac=use_jac,
                x0=x0,
                init_strategy=init_strategy,
                lambda_steps=lambda_steps,
            )

        _t0 = time.time()
        warm_seeded = False
        # Warm-start proxy roots the second-law check refused, for the
        # failure message (_proxy_refusal_note).
        self._proxy_refusals = []
        if (
            init_strategy == "default"
            and has_compressible
            and any(getattr(e, "prefers_warm_start", False) for e in self.network.elements.values())
        ):
            # Elements that converge slowly from the cold start declare it
            # (#706): try the outlet-referenced warm start first, and fall
            # back to the cold path below, with its own auto-retry, only if
            # it fails. Its warnings are held back the same way.
            warm_timeout = None if timeout is None else 0.5 * timeout
            with warnings.catch_warnings(record=True) as warm_warnings:
                warnings.simplefilter("always")
                warm = self._solve_impl(
                    method=method,
                    timeout=warm_timeout,
                    options=options,
                    use_jac=use_jac,
                    x0=None,
                    init_strategy="outletref_warmstart",
                    lambda_steps=lambda_steps,
                )
            if warm.get("__success__", False):
                for w in warm_warnings:
                    warnings.warn_explicit(w.message, w.category, w.filename, w.lineno)
                return warm
            # Its proxy converged: the outlet-referenced retry below would
            # solve from that same seed again (#1018: 7 s of a 25 s solve).
            warm_seeded = bool(getattr(self, "_outletref_seeded", False))
            if timeout is not None:
                timeout = max(timeout - (time.time() - _t0), 0.0)
            _t0 = time.time()
        # Junction (MultiPortChamberBase) networks reserve budget for the
        # retry up front: their cold path has slow pockets the seeded retry
        # rescues (see solve). Any other network keeps its whole budget for
        # the cold attempt, which is what converges it nearly always, and
        # retries only in the time left (#585).
        split_budget = any(
            isinstance(e, MultiPortChamberBase) for e in self.network.elements.values()
        )
        primary_timeout = 0.4 * timeout if timeout is not None and split_budget else timeout
        # The cold attempt's warnings are held back: its "did not converge"
        # is not the outcome when the seeded retry runs (a coupled combustor
        # liner warned on every operating point the retry then converged).
        with warnings.catch_warnings(record=True) as primary_warnings:
            warnings.simplefilter("always")
            sol = self._solve_impl(
                method=method,
                timeout=primary_timeout,
                options=options,
                use_jac=use_jac,
                x0=None,
                init_strategy=init_strategy,
                lambda_steps=lambda_steps,
            )

        def replay(retrying: bool) -> dict[str, float]:
            for w in primary_warnings:
                if retrying and "did not converge" in str(w.message):
                    continue
                warnings.warn_explicit(w.message, w.category, w.filename, w.lineno)
            return sol

        if sol.get("__success__", False):
            return replay(retrying=False)
        left = None if timeout is None else timeout - (time.time() - _t0)
        if not split_budget and left is not None and left <= 0.0:
            return replay(retrying=False)

        primary_norm = sol.get("__final_norm__")
        primary_msg = str(sol.get("__message__", "")).strip()
        # A primary that converged but was refused (outcome INCONSISTENT: a
        # second-law root, no cooling path, ...) still gets the retry -- it can
        # reach another root -- but its message already says "Converged, but
        # ...": the retry's notes must not call it a failed solve.
        primary_converged = bool(sol.get("__converged__", False))
        primary_what = "converged to a refused root" if primary_converged else "solve failed"
        warnings.warn(
            f"Cold '{init_strategy}' {primary_what} (|F|={primary_norm:.3e}); "
            + (
                "auto-retrying from a compressible warm start"
                if all_incompressible
                else "auto-retrying from an outlet-referenced incompressible warm start"
            )
            + ". Pass auto_retry=False to disable.",
            RuntimeWarning,
            stacklevel=3,
        )
        # Heavy networks (high Re, thermal surfaces) can need tens of
        # seconds even for the incompressible proxy (observed 35 s on a
        # 16 kg/s 20 bar 5-tee net), so give it up to 60 s / 30% of budget.
        proxy_timeout = 60.0 if timeout is None else max(min(0.3 * timeout, 60.0), 10.0)
        if not split_budget and left is not None:
            proxy_timeout = min(proxy_timeout, 0.5 * left)
        _primary_diag = getattr(self, "_diagnostic_data", None)
        if all_incompressible:
            seed_x0 = self._compressible_seed(use_jac=use_jac, timeout=proxy_timeout)
            if seed_x0 is None:
                if _primary_diag is not None:
                    self._diagnostic_data = _primary_diag
                if not primary_converged:
                    sol["__message__"] = (
                        f"{primary_msg} [compressible warm-start auto-retry "
                        "skipped: the compressible proxy solve did not converge]"
                    )
                return replay(retrying=False)
            _seed_kind = "compressible"
        elif warm_seeded:
            seed_x0 = None  # tried first, above: on to the inlet-referenced seed
        else:
            seed_x0 = self._outlet_ref_incompressible_seed(use_jac=use_jac, timeout=proxy_timeout)
            _seed_kind = "outlet-referenced"
        if seed_x0 is None:
            # The outlet-referenced system is stiffer and can time out where
            # the classic inlet-referenced proxy solves in a handful of
            # evaluations; its (coarser) seed still rescues heavy MFB-driven
            # networks, so fall back to it before giving up.
            seed_x0 = self._outlet_ref_incompressible_seed(
                use_jac=use_jac, timeout=proxy_timeout, p_ref="inlet"
            )
            _seed_kind = (
                "inlet-referenced (the outlet-referenced warm start ran first)"
                if warm_seeded
                else "inlet-referenced (outlet-ref proxy did not converge)"
            )
        if seed_x0 is None:
            # Restore the primary attempt's diagnostics -- the failed proxy
            # solves overwrote them, and we are returning the primary result.
            if _primary_diag is not None:
                self._diagnostic_data = _primary_diag
            # Skipped only matters to a solve that did not converge; on a
            # refused root it read as "did not converge" after "Converged,
            # but ..." (the no-cooling-path refusal, #998).
            if not primary_converged:
                sol["__message__"] = f"{primary_msg} [incompressible warm-start auto-retry " + (
                    "skipped: the outlet-referenced warm start ran first and "
                    "the inlet-referenced proxy solve did not converge]"
                    if warm_seeded
                    else "skipped: neither the outlet- nor the inlet-referenced "
                    "proxy solve converged]"
                )
            sol["__message__"] = self._with_proxy_refusal_note(sol.get("__message__", ""))
            return replay(retrying=False)

        remaining = None
        if timeout is not None:
            # The proxy took at most half of what was left, so the seeded
            # retry keeps at least the other half.
            floor = 0.3 * timeout if split_budget else 0.0
            remaining = max(timeout - (time.time() - _t0), floor)
        replay(retrying=True)
        retry_sol = self._solve_impl(
            method=method,
            timeout=remaining,
            options=options,
            use_jac=use_jac,
            x0=seed_x0,
            init_strategy="default",
        )
        # Stamp what actually ran so the GUI badge and headless tooling
        # report the effective strategy.
        _diag = getattr(self, "_diagnostic_data", None)
        if isinstance(_diag, dict):
            _ssu = _diag.get("solver_settings_used")
            if isinstance(_ssu, dict):
                _ssu["init_strategy"] = (
                    "compressible_warmstart"
                    if _seed_kind == "compressible"
                    else "outletref_warmstart"
                    if _seed_kind.startswith("outlet")
                    else "inletref_warmstart"
                )
                _ssu["auto_retry"] = True
        if retry_sol.get("__success__", False):
            proxy = (
                "compressible" if _seed_kind == "compressible" else f"{_seed_kind} incompressible"
            )
            retry_sol["__message__"] = (
                f"Converged via {proxy} warm-start "
                f"auto-retry after the '{init_strategy}' attempt "
                f"{'reached a refused root' if primary_converged else 'failed'} "
                f"({primary_msg[:120]})."
            )
            return retry_sol
        primary_end = (
            "converged to a refused root"
            if primary_converged
            else f"failed at |F|={primary_norm:.3e}"
        )
        retry_sol["__message__"] = (
            f"{str(retry_sol.get('__message__', '')).strip()} "
            f"[{_seed_kind} warm-start auto-retry; the primary "
            f"'{init_strategy}' attempt {primary_end} "
            f"with: {primary_msg[:120]}]"
        )
        retry_sol["__message__"] = self._with_proxy_refusal_note(retry_sol["__message__"])
        # Neither attempt converged, so hand back whichever got closer. The
        # retry used to be returned unconditionally, which could replace a
        # near-miss with something far worse and still warn that the best
        # iterate was being returned.
        retry_norm = retry_sol.get("__final_norm__")
        if (
            primary_norm is not None
            and retry_norm is not None
            and float(retry_norm) > float(primary_norm)
        ):
            if _primary_diag is not None:
                self._diagnostic_data = _primary_diag
            sol["__message__"] = (
                f"{primary_msg} [{_seed_kind} warm-start auto-retry also "
                f"failed, at |F|={float(retry_norm):.3e}, so the closer "
                f"primary result is returned]"
                if not primary_converged
                else f"{primary_msg} [{_seed_kind} warm-start auto-retry found "
                f"no admissible root either (|F|={float(retry_norm):.3e})]"
            )
            sol["__message__"] = self._with_proxy_refusal_note(sol["__message__"])
            return sol
        return retry_sol

    def _record_proxy_refusal(self, kind: str, proxy_sol: dict[str, Any]) -> None:
        """Keep a warm-start proxy root the second-law check refused (the
        root converged, its __second_law_failures__ name the elements) for
        the solve's failure message: a refused proxy is otherwise reported
        as one that "did not converge"."""
        bad = list(proxy_sol.get("__second_law_failures__", []) or [])
        if not bad or bool(proxy_sol.get("__success__", False)):
            return
        verdicts = proxy_sol.get("__second_law__", {}) or {}
        worst = {
            eid: float((verdicts.get(eid) or {}).get("sigma_norm_high", math.nan)) for eid in bad
        }
        refusals = getattr(self, "_proxy_refusals", None)
        if refusals is None:
            refusals = self._proxy_refusals = []
        refusals.append((kind, worst))

    def _with_proxy_refusal_note(self, message: str) -> str:
        """``message`` with the warm-start proxy roots the second law refused
        (each element and its sigma_norm_high). A duct-fed hole among them
        points at #1125: booked at the duct total, a station whose duct flow
        reverses manufactures total pressure."""
        refusals = getattr(self, "_proxy_refusals", None) or []
        if not refusals:
            return message
        parts = []
        duct_fed = False
        for kind, worst in refusals:
            els = ", ".join(f"{eid} (sigma_norm {s:+.2e})" for eid, s in worst.items())
            parts.append(
                f"the {kind} incompressible proxy root was refused by the second law at {els}"
            )
            for eid in worst:
                el = self.network.elements.get(eid)
                chain = getattr(el, "_hole_cd", None)
                duct_fed = (
                    duct_fed
                    or bool(getattr(chain, "duct_fed", False))
                    or isinstance(el, CrossflowSegmentElement)
                )
        note = "; ".join(parts)
        if duct_fed:
            note += (
                ": a backside duct whose flow reverses at a station, where the "
                "station's reversed-flow model manufactures total pressure once its "
                "duct-fed holes are booked at the duct total (known limitation, "
                "issue #1125)"
            )
        return f"{str(message).strip()} [{note}]"

    def _incompressible_element_ids(self) -> list[str]:
        """Elements running the incompressible form of a two-regime closure."""
        return [
            eid
            for eid, e in self.network.elements.items()
            if getattr(e, "regime", None) == "incompressible"
        ]

    def _compressible_seed(self, use_jac: bool = True, timeout: float = 30.0) -> np.ndarray | None:
        """Solve the network with its incompressible elements compressible.

        The mirror of the incompressible proxy, for a network that is
        incompressible throughout: at low Mach the two regimes share their
        root to O(Ma^2) (#588), so the compressible solution is a near-root
        of the incompressible system -- and its mass-flow residual forms
        converge from cold starts where the pressure-drop forms stall (a jet
        pump behind an incompressible pipe sat at zero entrainment from every
        cold seed, and converged in a few steps from the compressible root).
        Returns the solution vector, or ``None`` when the proxy fails.
        """
        ids = self._incompressible_element_ids()
        if not ids:
            return None
        seed: np.ndarray | None = None
        # Each element may ask for proxy-only settings (the ejector takes its
        # wider anchor window); restored with the regime below.
        saved: list[tuple[Any, str, Any]] = []
        try:
            for eid in ids:
                el = self.network.elements[eid]
                overrides = getattr(el, "compressible_proxy_overrides", None)
                for attr, value in (overrides() if callable(overrides) else {}).items():
                    saved.append((el, attr, getattr(el, attr)))
                    setattr(el, attr, value)
                el.regime = "compressible"
            self.network.resolve_all_topology()
            # Started as a plain compressible solve would be (cold-start
            # auto-upgrades included): that is the solve known to converge.
            proxy_sol = self._solve_impl(
                method="hybr",
                timeout=timeout,
                use_jac=use_jac,
                x0=None,
                init_strategy="default",
            )
            if bool(proxy_sol.get("__success__", False)):
                candidate = getattr(self, "_last_x", None)
                if candidate is not None:
                    seed = np.asarray(candidate, dtype=float)
        except Exception as e:
            warnings.warn(
                f"Compressible proxy solve raised; skipping warm-start auto-retry: {e}",
                RuntimeWarning,
                stacklevel=2,
            )
            seed = None
        finally:
            for eid in ids:
                self.network.elements[eid].regime = "incompressible"
            for el, attr, value in saved:
                setattr(el, attr, value)
            self.network.resolve_all_topology()
        return seed

    def _outlet_ref_incompressible_seed(
        self,
        use_jac: bool = True,
        timeout: float = 30.0,
        p_ref: str = "outlet",
    ) -> np.ndarray | None:
        """Solve the incompressible proxy; return its solution vector.

        Temporarily flips every compressible element to its incompressible
        fallback regime and sets ``_incompressible_p_ref = p_ref`` on all
        channel/orifice elements. ``p_ref='outlet'`` (default) evaluates
        densities at the downstream static -- the seed that tracks the
        compressible solution best; ``p_ref='inlet'`` is the classic,
        easier-to-solve proxy used as a fallback. Returns the proxy
        solution vector for use as a compressible warm start, or ``None``
        when the proxy fails.
        """
        overrides = self._compressible_element_overrides()
        if not overrides:
            return None
        original_regimes = {
            elem_id: getattr(self.network.elements[elem_id], "regime", None)
            for elem_id in overrides
        }
        ref_elems = [
            e
            for e in self.network.elements.values()
            if isinstance(e, (ChannelElement, OrificeElement))
        ]
        seed: np.ndarray | None = None
        try:
            for elem_id, fallback_regime in overrides.items():
                self.network.elements[elem_id].regime = fallback_regime
            for e in ref_elems:
                e._incompressible_p_ref = p_ref
            self.network.resolve_all_topology()
            self._in_warmstart_proxy = True
            # No maxfev cap: the outlet-referenced system is stiffer than
            # the classic inlet-referenced one (couples pressures through
            # the downstream density) and can need a few hundred cheap
            # incompressible evaluations. A tight timeout would push the
            # solve into the slower LM fallback phase and starve it.
            proxy_sol = self._solve_impl(
                method="hybr",
                timeout=timeout,
                use_jac=use_jac,
                x0=None,
                init_strategy="default",
            )
            if bool(proxy_sol.get("__success__", False)):
                candidate = getattr(self, "_last_x", None)
                if candidate is not None:
                    seed = np.asarray(candidate, dtype=float)
            else:
                self._record_proxy_refusal(f"{p_ref}-referenced", proxy_sol)
        except Exception as e:
            warnings.warn(
                "Outlet-referenced incompressible proxy solve raised; "
                f"skipping warm-start auto-retry: {e}",
                RuntimeWarning,
                stacklevel=2,
            )
            seed = None
        finally:
            self._in_warmstart_proxy = False
            for elem_id, regime in original_regimes.items():
                self.network.elements[elem_id].regime = regime
            for e in ref_elems:
                e._incompressible_p_ref = "inlet"
            self.network.resolve_all_topology()
        return seed

    def _solve_impl(
        self,
        method: str = "hybr",
        timeout: float | None = None,
        options: dict[str, Any] | None = None,
        use_jac: bool = True,
        x0: np.ndarray | None = None,
        init_strategy: Literal[
            "default",
            "outletref_warmstart",
            "homotopy",
            "analytical_pt_prop",
        ] = "default",
        lambda_steps: list[float] | None = None,
    ) -> dict[str, float]:
        """
        Single solve attempt (no auto-retry). See ``solve`` for the
        public entry point and argument semantics.

        Args:
            method: The scipy.optimize.root method to use.
            timeout: Maximum wall-clock time in seconds before stopping
                and returning the best iterate.
            options: Dictionary of solver options (e.g., maxiter, maxfev,
                xtol).
            use_jac: Whether to supply the analytical Jacobian to the
                solver (default True).
            x0: Optional initial guess vector.  When ``None`` (default),
                an initial guess is built automatically.  Pass a
                previous solution's ``_last_x`` for warm-starting
                parameter sweeps.
            init_strategy: Initialization strategy. ``default`` uses
                direct x0 construction; for networks containing a
                ``MultiPortChamberBase`` AND at least two
                ``PressureBoundary`` nodes it auto-upgrades to
                ``analytical_pt_prop`` (32/32 vs 28/32 certified-root
                convergence on the 2026-07 inverse-design audit,
                ``tmp/mpce_audit_v2_runner.py``) -- pass another
                strategy or an explicit ``x0`` to opt out. The >= 2 PB
                requirement exists because the strategy's Pt
                propagation needs a boundary pressure gradient: on
                MFB-driven networks it seeds reversed channel m_dots
                from path-length artifacts.
                ``analytical_pt_prop`` seeds topology-aware initial
                guesses (channel Bernoulli m_dot + MPCE junction P_jct
                at the common port). ``outletref_warmstart`` solves the outlet-referenced
                incompressible proxy (densities at the DOWNSTREAM
                static; the auto-retry's seed) and warm-starts from it
                directly -- 16/16 same-root on branch topologies,
                11/16 on merge in the certified audit.
            lambda_steps: Optional list of mass flow scaling factors for
                ``init_strategy='homotopy'``. Defaults to 10 steps
                from 0.1 to 1.0.
        """
        if method not in self.SUPPORTED_METHODS:
            raise ValueError(
                f"Method '{method}' is not supported. Supported methods are: {', '.join(self.SUPPORTED_METHODS)}"
            )

        if init_strategy not in (
            "default",
            "outletref_warmstart",
            "homotopy",
            "continuation",
            "analytical_pt_prop",
        ):
            raise ValueError(
                f"init_strategy must be one of: 'default', 'outletref_warmstart', "
                f"'homotopy', 'continuation', 'analytical_pt_prop'; got {init_strategy!r}."
            )

        options = {} if options is None else dict(options)

        # Normalise caller option names across methods.
        # scipy.root uses different names for the iteration limit:
        # - hybr: maxfev
        # - lm and quasi-Newton/fixed-point methods: maxiter
        # Translate before applying defaults to avoid stale duplicates.
        if method != "hybr" and "maxfev" in options:
            options.setdefault("maxiter", options.pop("maxfev"))

        # Ensure topology is resolved before solving
        self.network.resolve_all_topology()
        self.network.validate()

        # MPCE networks auto-upgrade the default initialization to the
        # analytical Pt-propagation seeds: on the certified inverse-design
        # audit (tmp/mpce_audit_v2_*, 2026-07-04) analytical_pt_prop
        # converged 32/32 to the certified root vs 28/32 for the plain
        # cold start, and was the fastest strategy overall. Explicitly
        # passing any other init_strategy (or an x0) opts out.
        # A warm-start proxy solve passes "default" and keeps the plain
        # cold start (flag below).
        # The upgrade additionally requires >= 2 PressureBoundary nodes:
        # the strategy's Pt propagation needs a boundary pressure
        # gradient to be meaningful. On MFB-driven networks with a
        # single PB, propagated per-node Pt differences are path-length
        # artifacts and the Bernoulli sign logic seeds REVERSED channel
        # m_dots, parking the first iterate inside the junctions'
        # soft-barrier region (observed on a real 5-tee GUI network:
        # |F(x0)| = 1.15e6 of pure penalty, hybr immediately stuck,
        # while the plain cold start converges in 16 evals).
        if (
            x0 is None
            and init_strategy == "default"
            and not getattr(self, "_in_warmstart_proxy", False)
        ):
            from .components import MultiPortChamberBase as _MPCElem

            _n_pb = sum(isinstance(n, PressureBoundary) for n in self.network.nodes.values())
            if _n_pb >= 2 and any(isinstance(e, _MPCElem) for e in self.network.elements.values()):
                init_strategy = "analytical_pt_prop"

        # analytical_pt_prop seeds channel m_dot + MPCE P_jct through the
        # regular _build_x0 path via self._init_overrides. _build_x0 reads
        # the dict once and clears it, so no cross-solve leakage.
        if x0 is None and init_strategy == "analytical_pt_prop":
            comps = self._flow_components()
            if len(comps) <= 1:
                _ref = self._infer_reference_state()
                self._init_overrides = self._propagate_analytical_pt_prop(_ref)
            else:
                # Per flow component, as _build_x0's reference state (#899).
                self._init_overrides = {}
                for node_ids, elem_ids in comps:
                    with self._component_view(node_ids, elem_ids):
                        self._init_overrides.update(
                            self._propagate_analytical_pt_prop(self._infer_reference_state())
                        )
        else:
            self._init_overrides = {}

        # Build initial guess to populate unknown_names early
        x0_auto = self._build_x0()
        self._init_overrides = {}
        if not self.unknown_names:
            return {
                "__success__": True,
                "__message__": "Network is empty or fully constrained (no unknowns).",
                "__iterations__": 0,
                "__converged__": True,
                "__consistent__": None,
                "__inconsistent_elements__": [],
                "__outcome__": SolveOutcome.NO_UNKNOWNS,
                "__worst_residuals__": [],
            }

        # Set default iteration limit based on problem size (matches hybr's own default
        # of 200*(n+1)) now that we know n.  Applied after the early-return check so
        # the network can't be empty when we read len(unknown_names).
        _n_unk = len(self.unknown_names)
        _default_iters = max(500, 200 * (_n_unk + 1))
        if method == "hybr":
            options.setdefault("maxfev", _default_iters)
            # Reduce initial trust-region step bound (hybr default is 100).
            # Smaller steps prevent large oscillating Newton steps near
            # ill-conditioned Jacobians (e.g., compressible flow near choking).
            options.setdefault("factor", 10)
        else:
            options.setdefault("maxiter", _default_iters)

        warmstart_x0: np.ndarray | None = None
        if x0 is None and init_strategy == "outletref_warmstart":
            # Explicit selection of the auto-retry's seed: solve the
            # outlet-referenced incompressible proxy and warm-start the
            # compressible solve from it. Same budget formula as the
            # retry path. When the proxy fails, warn and proceed cold.
            _proxy_timeout = 60.0 if timeout is None else max(min(0.3 * timeout, 60.0), 10.0)
            warmstart_x0 = self._outlet_ref_incompressible_seed(
                use_jac=use_jac, timeout=_proxy_timeout
            )
            self._outletref_seeded = warmstart_x0 is not None
            if warmstart_x0 is None:
                warnings.warn(
                    "Outlet-referenced incompressible proxy did not converge. "
                    "Proceeding with default initialization.",
                    RuntimeWarning,
                    stacklevel=2,
                )

        if x0 is None and init_strategy == "homotopy":
            # Identify all MassFlowBoundary nodes and their target m_dot
            targets = {
                nid: node.m_dot
                for nid, node in self.network.nodes.items()
                if isinstance(node, MassFlowBoundary) and not isinstance(node, FuelBoundary)
            }

            if targets:
                # Sequential solve loop with increasing mass flow lambda.
                # Use more steps for better robustness near choking.
                if lambda_steps is None:
                    lambda_steps = [0.2, 0.4, 0.6, 0.8, 0.9, 1.0]

                current_x0 = None
                last_sol = {"__success__": False}
                # The whole ladder shares the caller's timeout as a
                # deadline. Previously each rung got the full budget, so
                # a slow ladder could run to ~6x the requested timeout --
                # blowing straight through the GUI's hard (soft + 10 s)
                # request timeout.
                _ladder_t0 = time.time()

                for lam in lambda_steps:
                    # Update boundaries
                    for nid, target_mdot in targets.items():
                        self.network.nodes[nid].m_dot = lam * target_mdot

                    # Solve intermediate step.
                    # Limit maxfev/maxiter for intermediate steps to stay fast.
                    # This is decremental to convergence: trade speed versus robustness near choking.
                    step_options = dict(options)
                    # options already has maxfev/_default_iters set above; honour
                    # it.  The setdefault here is only reached if the caller
                    # explicitly zeroed out maxfev/maxiter, which should not happen.

                    rung_timeout = (
                        None if timeout is None else max(timeout - (time.time() - _ladder_t0), 1.0)
                    )
                    last_sol = self._solve_impl(
                        method=method,
                        timeout=rung_timeout,
                        options=step_options,
                        use_jac=use_jac,
                        x0=current_x0,
                        init_strategy="default",
                    )

                    # Update x0 for next step anyway (best guess propagation)
                    current_x0 = self._last_x

                    if not last_sol["__success__"]:
                        # If an intermediate step fails badly, we might still
                        # want to continue if progress was made, but for now
                        # we stop to avoid compounding errors.
                        # Note: we still use the failed result as x0 if it's the best we have.
                        break

                # Restore original m_dot (lam=1.0 loop already did this if it reached it)
                for nid, target_mdot in targets.items():
                    self.network.nodes[nid].m_dot = target_mdot

                # If homotopy reached lam=1.0 and succeeded, return that solution
                if last_sol["__success__"]:
                    return last_sol

                # If it failed at lam=1.0 but current_x0 is much better than cold start,
                # we will naturally use it as warmstart_x0 below.
                warmstart_x0 = current_x0
        elif x0 is None and init_strategy == "continuation":
            raise ValueError(
                "continuation init: no prior converged solution provided via 'x0'. "
                "The caller is responsible for persisting and passing the previous 'x' vector."
            )

        # --- Initial guess: user-supplied (warm-start) or auto-built ----
        if x0 is not None:
            if len(x0) != len(x0_auto):
                raise ValueError(
                    f"Warm-start x0 length ({len(x0)}) does not match "
                    f"number of unknowns ({len(x0_auto)})."
                )
            x0_use = np.asarray(x0, dtype=float)
        elif warmstart_x0 is not None:
            x0_use = warmstart_x0 if len(warmstart_x0) == len(x0_auto) else x0_auto
        else:
            x0_use = x0_auto

        # Seed forward-propagated thermo/composition states once before the
        # first residual/scaling pass by ensuring the first residual evaluation
        # triggers the propagation internally.

        # Dimension check. The starting point is the one evaluation the root
        # finder's penalty does not cover: an exception here used to abort the
        # whole solve, retries included (#552: the mixer's T(h) inverse
        # refusing a heat withdrawal larger than the stream holds).
        self._eval_phase = "start"
        try:
            res0 = self._residuals(x0_use)
            if not np.all(np.isfinite(res0)):
                # A NaN start raised nothing and was then reported converged:
                # |F| = NaN fails every "<" in the best-iterate tracking and
                # the "> tol" success guard alike. Seen on an ejector seeded
                # from its incompressible proxy, whose P_jct (-4e4 Pa) means
                # something else in the compressible form (stress sweep,
                # 2026-10-17).
                bad = np.flatnonzero(~np.isfinite(res0))
                try:
                    kinds = self._row_kinds(x0_use)
                    owners = sorted({kinds[i][0] for i in bad if i < len(kinds)})
                except Exception:  # noqa: BLE001 -- naming only
                    owners = []
                self._eval_where = f"rows of {', '.join(owners[:4])}" if owners else ""
                raise FloatingPointError(f"non-finite residual in {len(bad)} row(s) at the start")
        except _STATE_EXCEPTIONS as exc:
            self._record_exception(exc)
            if self._strict():
                raise
            warnings.warn(
                f"NetworkSolver: the starting point could not be evaluated"
                f"{' at ' + self._eval_where if self._eval_where else ''} "
                f"({type(exc).__name__}: {str(exc)[:200]}); this attempt ends "
                "with outcome ERROR.",
                RuntimeWarning,
                stacklevel=2,
            )
            return self._failed_start_result(x0_use, exc)
        finally:
            self._eval_phase = "probe"
        if len(x0_use) != len(res0):
            raise ValueError(f"System not square: {len(x0_use)} unknowns vs {len(res0)} equations.")

        # --- Variable & residual scaling (RESIDUAL-ROW-SCALING) ----------
        # Unknowns:  x_real = x_scaled * D_x            (D_x = |x0|)
        # Residuals: F_scaled = F_real / D_f = F_real * inv_D_f
        # D_f is per row, by the row's units: see _build_residual_scales,
        # which is the canonical description. (This comment used to say
        # F_real * D_f, the inverse of what the code does.)
        D_x = np.abs(x0_use).copy()
        D_x[D_x < 1e-12] = 1.0  # avoid division by zero for near-zero guesses
        D_f = self._build_residual_scales(x0_use)
        inv_D_x = 1.0 / D_x
        inv_D_f = 1.0 / D_f

        x0_scaled = x0_use * inv_D_x  # should be ~1.0

        # For compressible networks hybr can oscillate near choking.  Reserve
        # 35% of the wall-clock budget for an LM fallback that has better
        # global convergence via (J^TJ + lI) regularisation.
        # Incompressible inner solves (regime overrides applied) have no
        # compressible elements and keep the full budget.
        # Tee junctions couple two supplier stagnation pressures whose
        # common mode leaves hybr's Jacobian near-singular; LM's regularised
        # step resolves it, so route tee networks through the same fallback.
        # Same near-singularity affects MultiPortChamberBase networks (the
        # P_jct unknown couples N port-pressure rows that share a common-mode
        # direction at low-Mach).
        from .components import MultiPortChamberBase as _MPCElem
        from .components import TeeJunctionElement as _TeeJE

        _has_tee = any(isinstance(e, _TeeJE) for e in self.network.elements.values())
        _has_mpce = any(isinstance(e, _MPCElem) for e in self.network.elements.values())
        # Merge chambers (#471) are junctions too: side-stream momentum feeds
        # back on the main inflow, and from a cold guess hybr can step into
        # reversed main flow where the merge is not a merge. LM from x0
        # converges (2026-10-07, 30-deg effusion into a liner chamber).
        _has_merge = any(
            getattr(n, "main_inlet", None) is not None for n in self.network.nodes.values()
        )
        _has_compressible = method == "hybr" and (
            bool(self._compressible_element_overrides()) or _has_tee or _has_mpce or _has_merge
        )
        _hybr_budget = timeout * 0.65 if (_has_compressible and timeout is not None) else timeout
        # Mutable so the LM fallback block can extend the effective limit.
        _timeout_eff: list[float | None] = [_hybr_budget]

        # Convergence history: record on every ||F|| improvement and at most
        # once per 0.1 s to show oscillation without excessive overhead.
        # Cap at 1 000 entries (~50 KB JSON) regardless of problem size.
        _history: list[dict] = []
        _last_record_t: list[float] = [0.0]
        _RECORD_INTERVAL = 0.1
        _eval_count: list[int] = [0]
        _lm_started_at_eval: list[int | None] = [None]
        # Per-eval (t, best|F|) traces for the stall detector; checked every
        # _STALL_CHECK_EVERY evals. The hybr trace uses solve-relative
        # times; the LM-fallback trace uses phase-relative times so the
        # detector's grace period applies from the start of the LM phase.
        _stall_samples: list[tuple[float, float]] = []
        _lm_stall_samples: list[tuple[float, float]] = []
        _STALL_CHECK_EVERY = 10
        _stall_fired: list[bool] = [False]
        _lm_stall_fired: list[bool] = [False]
        _lm_start_t: list[float] = [0.0]
        # Stall detection (both phases) only makes sense for COLD
        # attempts, where a later rung (LM fallback / warm-start retry)
        # exists to receive the freed budget. Seeded solves (explicit
        # x0 or a warmstart proxy seed) are the LAST rung and must keep
        # their full grind: on the certified audit the seeded LM phase
        # can plateau far from tolerance for > 5 s and still dig out to
        # the root (4/32 outletref cases), and the hybr handover firing
        # mid-seeded-solve turned two 20 s successes into failures
        # (audit cases 13/19) because LM cannot always finish from the
        # same seed.
        # solve(stall_handover=False) turns both handovers off (#864).
        _cold_attempt = x0 is None and warmstart_x0 is None and self._stall_handover

        def _stall_reason(samples: list[tuple[float, float]], t_solve: float) -> str:
            # The plateau/crawl detector first, then the projected-time rule;
            # the remaining budget is the attempt's, from solve-relative time.
            if _stall_detected(samples, tol=_RESIDUAL_TOL):
                return (
                    f"improved less than {_STALL_MIN_REL_IMPROVEMENT:.0%} over the "
                    f"last {_STALL_WINDOW_S:.0f} s"
                )
            remaining = None if timeout is None else timeout - t_solve
            if _crawl_too_slow(samples, remaining, tol=_RESIDUAL_TOL):
                return (
                    f"projects past {_STALL_PROJ_FACTOR:.0f}x the remaining "
                    "budget to reach tolerance"
                )
            return ""

        # State for timeout and best iterate tracking (in real space)
        start_time = time.perf_counter()
        best_x = x0_use.copy()
        best_res_norm = np.linalg.norm(res0)
        last_exception: Exception | None = None
        # Last PHYSICAL evaluation (scaled F, scaled J): what a rejected
        # probe is answered with (see the except branch below).
        last_ok: list[Any] = [None, None]

        # LAZY JACOBIAN (#489). The residual callback returns F only and the
        # Jacobian is a separate callable, assembled only when the root
        # finder asks for it. MINPACK hybrj asks at the start and on
        # restarts, running Broyden updates in between: 74 of 401
        # evaluations on a 36-case mixing study. Returning (F, J) from one
        # callback (scipy's MemoizeJac) assembled J on every evaluation.
        def _jac_for(m: str) -> Any:
            if not use_jac or m not in ("hybr", "lm"):
                return False
            # Both MINPACK drivers get the trace-flushed Jacobian: see
            # flushed_jacobian_wrapper (#1020 for lmder, #1027 for hybrj).
            return flushed_jacobian_wrapper

        n_rejected = [0]

        # The last (x_scaled, F_scaled) answered: scipy evaluates F at x0
        # once to check its shape and the root finder again to start.
        last_call: list[Any] = [None, None]

        def residuals_wrapper(x_scaled: np.ndarray) -> Any:
            if last_call[0] is not None and np.array_equal(last_call[0], x_scaled):
                return last_call[1]
            f = _residuals_scaled(x_scaled)
            last_call[0], last_call[1] = x_scaled.copy(), f
            return f

        def _residuals_scaled(x_scaled: np.ndarray) -> Any:
            nonlocal best_x, best_res_norm, last_exception

            # Check timeout
            if _timeout_eff[0] is not None and (time.perf_counter() - start_time) > _timeout_eff[0]:
                raise SolverTimeoutError(f"Solver timed out after {_timeout_eff[0]:.1f} seconds.")

            # Un-scale to real space
            x_real = x_scaled * D_x

            try:
                # Residuals only: the Jacobian is jacobian_wrapper's job.
                res, _ = self._residuals_and_jacobian(x_real, compute_jacobian=False)
                if not np.all(np.isfinite(res)):
                    raise FloatingPointError("non-finite residual at a solver iterate")

                # Track best iterate (in real space, unscaled residual)
                res_norm = np.linalg.norm(res)
                _improved = res_norm < best_res_norm
                if _improved:
                    best_res_norm = res_norm
                    best_x = x_real.copy()

                # Record convergence history point
                _t = time.perf_counter() - start_time
                if (_improved or _t - _last_record_t[0] >= _RECORD_INTERVAL) and len(
                    _history
                ) < 1000:
                    _history.append({"eval": _eval_count[0], "t": round(_t, 3), "norm": res_norm})
                    _last_record_t[0] = _t
                _eval_count[0] += 1

                # Stall-triggered handover to the LM fallback: when hybr's
                # best |F| plateaus there is no point burning the rest of
                # its budget share -- LM (which is what actually cracks
                # these networks) should get the time instead. Active only
                # while a compressible/tee LM fallback exists and has not
                # started yet, and only on cold attempts -- see
                # _cold_attempt above.
                if (
                    method == "hybr"
                    and _has_compressible
                    and _cold_attempt
                    and _lm_started_at_eval[0] is None
                ):
                    _stall_samples.append((_t, float(best_res_norm)))
                    if _eval_count[0] % _STALL_CHECK_EVERY == 0 and (
                        _why := _stall_reason(_stall_samples, _t)
                    ):
                        _stall_fired[0] = True
                        raise SolverStallError(
                            f"hybr stalled: best |F|={best_res_norm:.3e} "
                            f"{_why}; handing over to the LM fallback."
                        )

                # Same detector for the LM fallback phase (phase-relative
                # times): when LM also plateaus far from tolerance the
                # primary attempt is doomed -- fail fast so the caller's
                # outlet-referenced warm-start auto-retry gets the budget
                # instead of an unbounded LM grind (observed 180 s flat at
                # |F| ~ 4e3 on nets the retry then solves in seconds).
                # Cold attempts only -- see _cold_attempt above.
                elif _has_compressible and _cold_attempt and _lm_started_at_eval[0] is not None:
                    _lm_stall_samples.append((_t - _lm_start_t[0], float(best_res_norm)))
                    if _eval_count[0] % _STALL_CHECK_EVERY == 0 and (
                        _why := _stall_reason(_lm_stall_samples, _t)
                    ):
                        _lm_stall_fired[0] = True
                        raise SolverStallError(
                            f"LM fallback stalled: best |F|={best_res_norm:.3e} "
                            f"{_why}; abandoning this attempt."
                        )

                # Scale residuals
                res_scaled = res * inv_D_f
                last_ok[0] = res_scaled
                return res_scaled
            except SolverTimeoutError:
                raise
            except Exception as e:
                last_exception = e
                self._record_exception(e, "probe")
                # Return a large penalty residual so the trust-region
                # solver can shrink its step instead of terminating.
                # Warn once so developers can trace unexpected failures.
                if not getattr(self, "_penalty_warning_issued", False):
                    warnings.warn(
                        "NetworkSolver hit an unphysical intermediate state; "
                        "returning penalty residual to allow step reduction. "
                        "Enable debug logging to see the original exception.",
                        RuntimeWarning,
                        stacklevel=2,
                    )
                    self._penalty_warning_issued = True
                n_rejected[0] += 1
                # A penalty the trust region must REJECT: the last physical
                # residual's direction at ten times its norm (at least 10 in
                # scaled units); jacobian_wrapper answers a Jacobian request
                # there with the last physical Jacobian. The step then
                # reads as a clear increase and the region shrinks. The former
                # F = x - x_best (J = I) was SMALL near x_best, so a step into
                # the unphysical region read as an improvement and was
                # accepted (#481).
                f_ref = last_ok[0]
                n_f = len(x_scaled)
                if f_ref is None:
                    f_ref = np.ones(n_f)
                ref_norm = float(np.linalg.norm(f_ref))
                direction = f_ref / ref_norm if ref_norm > 0.0 else np.ones(n_f) / np.sqrt(n_f)
                target = 10.0 * max(ref_norm, 1.0)
                return direction * target

        def jacobian_wrapper(x_scaled: np.ndarray) -> np.ndarray:
            """Scaled Jacobian at x_scaled, assembled on request only."""
            if _timeout_eff[0] is not None and (time.perf_counter() - start_time) > _timeout_eff[0]:
                raise SolverTimeoutError(f"Solver timed out after {_timeout_eff[0]:.1f} seconds.")
            try:
                x_real = x_scaled * D_x
                # Usually the point F was just evaluated at: assemble its
                # deferred Jacobian instead of evaluating again.
                jac = self._jacobian_at_last(x_real)
                if jac is not None and np.array_equal(last_call[0], x_scaled):
                    f_scaled = last_call[1]
                else:
                    res, jac = self._residuals_and_jacobian(x_real)
                    f_scaled = res * inv_D_f
                if not np.all(np.isfinite(jac.data)) or not np.all(np.isfinite(f_scaled)):
                    # Name the rows and columns: the last element evaluated
                    # is not where a non-finite entry came from.
                    coo = jac.tocoo()
                    bad = ~np.isfinite(coo.data)
                    rows = sorted(
                        set(coo.row[bad].tolist())
                        | set(np.flatnonzero(~np.isfinite(f_scaled)).tolist())
                    )
                    cols = sorted(set(coo.col[bad].tolist()))
                    names = self.unknown_names
                    self._eval_where = (
                        f"Jacobian rows {[names[i] for i in rows[:4]]} "
                        f"cols {[names[j] for j in cols[:4]]}"
                    )
                    raise FloatingPointError("non-finite Jacobian at a solver iterate")
                # J_scaled = diag(inv_D_f) @ J @ diag(D_x)
                jac_scaled = (inv_D_f[:, None] * jac.toarray()) * D_x[None, :]
                # The root finder asks for J at the point it has ACCEPTED, so
                # that point's residual is the reference a later rejected
                # probe is answered from (the penalty above).
                last_ok[0], last_ok[1] = f_scaled, jac_scaled
                return jac_scaled
            except SolverTimeoutError:
                raise
            except Exception as e:
                # A point the model refuses: the last physical Jacobian, as
                # for the rejected-probe residual above.
                self._record_exception(e, "probe_jacobian")
                return last_ok[1] if last_ok[1] is not None else np.eye(len(x_scaled))

        def flushed_jacobian_wrapper(x_scaled: np.ndarray) -> np.ndarray:
            """The scaled Jacobian for MINPACK, trace entries flushed to 0.

            lmder detects a singular Jacobian by an EXACT zero on R's diagonal
            (lmpar's nsing test) and otherwise divides by that pivot. A column
            that is singular but for an entry at the underflow level passes
            the test with a pivot of ~1e-308; Q^T F / R_jj overflows to inf,
            the step to NaN, and the trust radius it updates (delta = pnorm /
            0.5) becomes NaN for good: every later iterate is NaN (#1020).
            hybrj's dogleg has the same guard: it replaces only an exactly
            zero R_jj (by eps times the column's largest R entry), so the
            Gauss-Newton step Q^T F / 7.6e-309 overflowed, alpha = delta /
            inf = 0 times inf gave a NaN step, delta went NaN, and the last 8
            of hybr's 56 evaluations on stress effusion 105 were NaN (#1027).
            The pivot was in the Jacobian hybrj had just requested, not in a
            Broyden-updated R. An entry below eps times its column's largest
            is below the column's own Householder rounding (QR is backward
            stable column by column), so zeroing it changes nothing QR can
            resolve and lets lmpar take its rank-deficient branch and dogleg
            its zero-pivot guard. The case: an effusion plate whose
            McGreehan-Schotsch crossflow Cd = exp(-Rv^n) underflows at a far
            iterate, so the only row that sets a chamber's pressure level
            carries derivatives of 1e-308.
            """
            jac = jacobian_wrapper(x_scaled)
            col_max = np.max(np.abs(jac), axis=0)
            trace = np.abs(jac) < np.finfo(float).eps * col_max[None, :]
            if np.any(trace & (jac != 0.0)):
                jac = np.where(trace, 0.0, jac)
            return jac

        # Solve with timing
        _RESIDUAL_TOL = 1e-3
        _EXACT_ROOT_TOL = 1e-12
        solve_start_time = time.time()
        outcome: SolveOutcome = SolveOutcome.NOT_CONVERGED
        try:
            sol = root(
                residuals_wrapper, x0_scaled, method=method, options=options, jac=_jac_for(method)
            )
            solve_end_time = time.time()
            self._wall_time = solve_end_time - solve_start_time
            final_x = best_x
            success = sol.success
            message = sol.message
            outcome = SolveOutcome.CONVERGED if success else _classify_solver_message(message)
            # Efficiently use tracked norm instead of re-evaluating
            final_norm = float(best_res_norm)
            # MINPACK's progress test misfires at an EXACT root: once |F| is 0
            # (or within rounding of it) its actual/predicted reduction ratio
            # is 0/0 and hybr reports "not making good progress" from the
            # root itself. A residual at machine level IS the solution (#582:
            # a row linear in its own unknown reaches |F| = 1e-19 in one step).
            if not success and final_norm <= _EXACT_ROOT_TOL:
                success = True
                outcome = SolveOutcome.CONVERGED
                message = (
                    f"Converged (|F|={final_norm:.1e}; {method} stopped on its "
                    f"progress test at the root)."
                )
            # Guard against methods (e.g. lm) that report success at a
            # non-zero local minimum of ||F||^2. Written NaN-safe: a NaN
            # norm is not a converged one.
            if success and not final_norm <= _RESIDUAL_TOL:
                success = False
                outcome = SolveOutcome.RESIDUAL_TOO_LARGE
                message = (
                    f"Solver reported success but |F|={final_norm:.3e} "
                    f"exceeds residual tolerance ({_RESIDUAL_TOL})."
                )
            if not success and last_exception is not None and _eval_count[0] == 0:
                # Not one physical evaluation since the start: that is a
                # failing model, not a search that strayed.
                raise last_exception
            if not success and last_exception is not None:
                # Rejected probes are part of the search, not the outcome:
                # say so, keep the classification of why the solve stopped.
                message = (
                    f"{message} [{n_rejected[0]} unphysical probe(s) rejected; "
                    f"last: {type(last_exception).__name__}: {str(last_exception)[:160]}]"
                )
        except SolverTimeoutError as e:
            final_x = best_x
            success = False
            outcome = SolveOutcome.TIMEOUT
            message = str(e)
            final_norm = float(best_res_norm)
        except Exception as e:
            # Fallback for unexpected errors during residuals evaluation
            final_x = best_x
            success = False
            outcome = SolveOutcome.ERROR
            message = f"Unexpected error during residual evaluation: {e}"
            final_norm = float(best_res_norm)

        # Automatic fallback: if the chosen method failed and isn't hybr,
        # retry with hybr (Powell's method) which is more robust for square
        # nonlinear systems.
        if not success and method != "hybr":
            try:
                sol2 = root(residuals_wrapper, x0_scaled, method="hybr", jac=_jac_for("hybr"))
                fallback_x = best_x
                fallback_norm = float(best_res_norm)
                if sol2.success and fallback_norm < _RESIDUAL_TOL:
                    final_x = fallback_x
                    success = True
                    outcome = SolveOutcome.CONVERGED
                    message = f"Converged after fallback to hybr (original {method} failed)."
                    final_norm = fallback_norm
            except Exception:
                pass  # keep original failure

        # LM fallback for hybr oscillation on compressible networks.
        # Levenberg-Marquardt's (J^TJ + lI) regularisation prevents the
        # ill-conditioned Newton step that causes hybr to bounce near a
        # choked solution without making progress.  Start from the best
        # iterate hybr found; restore the full timeout so LM gets the
        # remaining 35% of the budget. Also for ANY network where hybr
        # reported success at a non-root (its xtol test on a Broyden-updated
        # Jacobian): an 8-row impingement array with the epsilon-NTU wall heat
        # stopped there at |F| = 1.7e-3, where LM from its best iterate reaches
        # 3e-11 (#627). Only a failed solve gets the extra phase.
        if not success and (_has_compressible or outcome == SolveOutcome.RESIDUAL_TOO_LARGE):
            _elapsed = time.perf_counter() - start_time
            _remaining = (timeout - _elapsed) if timeout is not None else None
            if _remaining is None or _remaining > 2.0:
                _lm_started_at_eval[0] = _eval_count[0]
                _lm_start_t[0] = time.perf_counter() - start_time
                _timeout_eff[0] = timeout
                # For MPCE networks the impulse + sum-mass residual structure
                # leaves hybr near-singular at low Mach; its "best iterate" is
                # often a worse LM starting point than the original x0. Restart
                # LM from x0 in that case to avoid getting trapped in hybr's
                # poor basin. For tee / compressible networks -- and merge
                # chambers, which this used to restart too -- the warm start
                # from hybr is better: on the GUI-default effusion liner hybr
                # stopped at |F| ~ 1e-2 with an exact Jacobian and a full
                # Newton step still cutting |F| by a third, and the cold LM
                # restart threw that away (#582).
                _lm_x0_scaled = x0_scaled if _has_mpce else best_x * inv_D_x
                try:
                    root(
                        residuals_wrapper,
                        _lm_x0_scaled,
                        method="lm",
                        jac=_jac_for("lm"),
                        options={"maxiter": _default_iters},
                    )
                    _lm_norm = float(best_res_norm)
                    if _lm_norm < _RESIDUAL_TOL:
                        final_x = best_x
                        success = True
                        outcome = SolveOutcome.CONVERGED
                        _hybr_end = "stalled at" if _stall_fired[0] else "|F|="
                        message = (
                            f"Converged with LM fallback "
                            f"(hybr {_hybr_end} {final_norm:.3e}, "
                            f"LM |F|={_lm_norm:.3e})."
                        )
                        final_norm = _lm_norm
                except SolverTimeoutError:
                    pass
                except Exception:
                    pass
                if not success and _lm_stall_fired[0]:
                    message = (
                        f"{message} LM fallback also stalled at |F|={float(best_res_norm):.3e}."
                    )

        # Honour the promise made in the non-convergence warning: return the
        # BEST iterate, not the one that happened to be current when the
        # primary phase ended.
        #
        # `final_x` and `final_norm` are captured immediately after the primary
        # root() call. Every later phase -- the hybr fallback, and above all
        # the LM fallback -- keeps evaluating through the same wrapper, so
        # `best_x` and `best_res_norm` go on improving; but those phases only
        # re-point `final_x` when they REACH the convergence tolerance. An
        # improvement that falls short of it was being thrown away. Measured
        # over 38 non-converged junction solves before this guard: 30 returned
        # a state worse than the best they had evaluated, by a median of 5.8x
        # and up to 2e5x, and in every case inspected the better point was
        # found after the LM fallback started.
        #
        # Placed before the consistency verification below so the junction
        # checks judge the state that is actually returned.
        if float(best_res_norm) < final_norm:
            final_x = best_x
            final_norm = float(best_res_norm)

        # Post-solve physical-consistency verification for junctions.
        # Two known ways a junction net converges (|F| ~ 1e-10) onto an
        # unphysical exact root: (a) MultiPortChamberElement strict=False soft-barrier
        # mode, where the one-sided penalty alpha * mdot^2 on a slightly
        # reversed port cancels that row's Pt-continuity mismatch
        # (exposed by constant-K warm starts on the certified audit,
        # 2026-07-04: a merge feed reversed to -0.0035 kg/s while the
        # strict residual raises at the same state); (b) v1 MPCE mirror
        # roots -- the impulse rows are even in the flows, so the
        # sign-flipped image of any root is also exact, and the image of
        # a physical root collects flow from low-Pt ports into a
        # higher-Pt port (a passive junction manufacturing flow work).
        # Each element's verify_solution_consistent encodes its own
        # criterion; demote failures to honest non-convergence.
        #
        # The verdict is REPORTED as well as applied. "Not checked" is not the
        # same as "fine": a network whose elements carry no verifier gets
        # ``__consistent__ = None``, never True, so a caller cannot read a
        # silent absence of checking as a clean bill of health.
        # The root finder's own verdict, before the consistency check can
        # demote it. ``__success__`` folds the two together and must keep
        # doing so; this is the half that says whether Newton got there.
        converged = bool(success)
        _bad_junctions: list[str] = []
        _checkable = [
            element
            for element in self.network.elements.values()
            if hasattr(element, "verify_solution_consistent")
        ]
        consistent: bool | None = None
        if success:
            _sol_names = dict(zip(self.unknown_names, final_x, strict=False))
            _bad_junctions = [
                eid
                for eid, element in self.network.elements.items()
                if hasattr(element, "verify_solution_consistent")
                and not element.verify_solution_consistent(_sol_names)
            ]
            if _checkable:
                consistent = not _bad_junctions
            if _bad_junctions:
                success = False
                outcome = SolveOutcome.INCONSISTENT
                message = (
                    "Converged to an unphysical artifact root at "
                    f"junction(s) {_bad_junctions}: the solution failed "
                    "the junction's physical-consistency verification "
                    "(wrong-direction port flow under a soft barrier, an "
                    "energy-inconsistent mirror root, or a state where the "
                    "junction's own loss closure is not net dissipative). "
                    "Treating as non-converged; "
                )
                # Suggest another start only where one is left to try: a root
                # the topology-aware seed itself converged to is the model's
                # answer for this geometry, and re-solving lands on it again
                # (#536: a tee port narrower than its channel, reproduced from
                # every start).
                if x0 is not None or init_strategy != "analytical_pt_prop":
                    message += (
                        "try init_strategy='analytical_pt_prop' or a different initial guess."
                    )
                else:
                    message += (
                        "this root was reached from the topology-aware start, so "
                        "a different start rarely changes it. Check the junction's "
                        "geometry first: a port area that differs from the "
                        "connected element's is the usual cause."
                    )

        # A node at or below zero pressure is not a gas state, whatever the
        # residuals say: an element's continuation can close its row there (a
        # Fanno kernel's C1 continuation below the sonic static: a crossflow
        # array on its false choke "converged" with the free exit at -113 kPa,
        # #810). Not admissible. Checked here, before anything evaluates the
        # state: the element diagnostics below raise on a converged state they
        # cannot evaluate, and that raise used to escape solve() before the
        # second-law check (NaN sigma) could demote the root. NaN-safe.
        _P_nonpositive: list[str] = []
        if success:
            _P_nonpositive = sorted(
                {
                    name.rsplit(".", 1)[0]
                    for name, val in zip(self.unknown_names, final_x, strict=False)
                    if name.rsplit(".", 1)[-1] in ("P", "Pt")
                    and name.rsplit(".", 1)[0] in self.network.nodes
                    and not float(val) > 0.0
                }
            )
            if _P_nonpositive:
                success = False
                consistent = False
                outcome = SolveOutcome.INCONSISTENT
                message = (
                    f"Converged with node pressure(s) at or below zero at "
                    f"{_P_nonpositive}: not a gas state. Usually an element "
                    "between them and the supply cannot pass the imposed flow "
                    "(choked or too restrictive); check the flow against it."
                )

        # A node at T(h) <= 0 K had more heat taken out than its flow holds at
        # 0 K (the gas is ideal down to it): the residuals can close on the
        # continuation while the state is not a gas. Not admissible, and
        # reported, not silently accepted (#481 C4; ``__T_clamped__`` keeps its
        # name). A stagnant node carries no enthalpy anywhere -- its heat is
        # reported as Q_withheld -- so only a node with inflow counts. Until
        # this, a [200, 5000] K clamp discarded the energy instead; a node
        # heated past 5000 K is now simply hot (the thermo's tails carry it),
        # energy conserved.
        _T_clamped: list[str] = []
        if success:
            self._propagate_states(final_x, with_relay=False)
            _T_clamped = [
                nid
                for nid in self._derived_states
                if nid in self._T_no_gas
                and sum(m for m, _, _ in self._node_inflows(nid, final_x))
                >= cb.MIXER_HEAT_MDOT_FLOOR
            ]
            if _T_clamped:
                success = False
                consistent = False
                outcome = SolveOutcome.INCONSISTENT
                message = (
                    f"Converged with node temperature(s) at or below 0 K at "
                    f"{_T_clamped}: more heat out than the stream holds as a gas. "
                    "Check heat inputs against the flow they heat."
                )

        # A heated node at a root inside the heat fade below
        # MIXER_HEAT_MDOT_FLOOR: neither flowing (all its heat taken up) nor
        # stagnant (none). In a driven net heating a stream past ~1e5 K makes
        # it flow less, until the fade, which takes up less heat the less it
        # flows, offers a second, spurious root (100 kW into 1 mg/s of air
        # under a 10 kPa drive: 5.8e7 K, 10 kW withheld).
        _heat_faded: list[str] = []
        if success:
            for nid, node in self.network.nodes.items():
                q_given = sum(
                    float(eb.heat) for eb in getattr(node, "energy_boundaries", None) or []
                )
                if q_given == 0.0:
                    continue
                m_in = sum(m for m, _, _ in self._node_inflows(nid, final_x))
                phi, _ = cb.mixer_heat_per_unit_flow(m_in)
                taken = m_in * phi
                if _HEAT_FADE_TOL < taken < 1.0 - _HEAT_FADE_TOL:
                    _heat_faded.append(nid)
            if _heat_faded:
                success = False
                consistent = False
                outcome = SolveOutcome.INCONSISTENT
                message = (
                    f"Converged with heated node(s) {_heat_faded} taking up only part of "
                    f"their heat: their flow is below MIXER_HEAT_MDOT_FLOOR "
                    f"({cb.MIXER_HEAT_MDOT_FLOOR:g} kg/s) but not stagnant, a root the "
                    "fade of the heat makes, not the network's. Check heat inputs "
                    "against the flow they heat."
                )

        # An EnergyBoundary(Tt_target) on a node with (next to) no through-flow
        # (#792): Q enters as Q phi(m) / cp, so below MIXER_HEAT_MDOT_FLOOR the
        # target row hardly sees Q and at zero flow not at all -- any Q "meets"
        # the target, or none does. Such a Q is not an answer: a converged solve
        # is demoted, a failed one says why, and neither reports the iterate's
        # Q (it is NaN below).
        _target_stagnant: list[str] = []
        for nid, cols in self._target_q_columns.items():
            try:
                m_in = sum(m for m, _, _ in self._node_inflows(nid, final_x))
            except Exception:  # noqa: BLE001 -- a diagnostic never breaks a solve
                continue
            if not m_in >= cb.MIXER_HEAT_MDOT_FLOOR:
                _target_stagnant.extend(eb.id for eb, _ in cols)
        if _target_stagnant:
            note = (
                f"EnergyBoundary {_target_stagnant} sets Tt_target on a node with no "
                f"through-flow (below MIXER_HEAT_MDOT_FLOOR, {cb.MIXER_HEAT_MDOT_FLOOR:g} "
                "kg/s): the heat that sets an outlet temperature is undetermined "
                "without flow, so no Q is reported. Drive flow through the node, or "
                "set Q instead."
            )
            if success:
                success = False
                consistent = False
                outcome = SolveOutcome.INCONSISTENT
                message = f"Converged, but {note[0].lower()}{note[1:]}"
            else:
                message = f"{message} {note}"

        # A momentum chamber on the supersonic root of its closure (the
        # states are those of the propagation above). Seen on a channel into a
        # junction before a wide orifice: the junction at Mach 1.46, the
        # channel "leaving" at Mach 1.67 unchoked, the orifice's static
        # pressure rising across it -- reported as converged and consistent.
        #
        # Since #788 the closure has ONE root per Pt (continued past the sonic
        # point, C++ momentum_chamber_residual_and_jacobian), so a chamber is
        # only ever there when its Pt is below the sonic Pt* of its flow: its
        # area cannot carry the flow at that stagnation pressure. That is a
        # CHOKED-BOUNDED junction, the under-expanded exit a pressure boundary
        # takes as "choked bounded" (#636), when a choked duct of no larger
        # area discharges into it and nothing reads its static pressure; it is
        # accepted and reported as such. Anything else is refused.
        _supersonic: list[str] = []
        _supersonic_notes: list[str] = []
        _choked_bounded: list[str] = []
        for node in self.network.nodes.values():
            if isinstance(node, MomentumChamberNode):
                node._choked_bounded = None
        if success:
            for nid, node in self.network.nodes.items():
                if (
                    not isinstance(node, MomentumChamberNode)
                    or node.regime == "incompressible"
                    or not node.area
                    or nid in self._T_no_gas
                ):
                    continue
                M = float(node.mach(self._get_node_state(node, final_x)))
                if M > 1.0 + _SUPERSONIC_CHAMBER_TOL:
                    bounded, why = self._choked_bounded_junction(nid, final_x)
                    if bounded is not None:
                        _choked_bounded.append(nid)
                        node._choked_bounded = bounded
                        continue
                    _supersonic.append(nid)
                    _supersonic_notes.append(f"'{nid}' ({why})")
            if _supersonic:
                success = False
                consistent = False
                outcome = SolveOutcome.INCONSISTENT
                message = (
                    f"Converged with momentum chamber(s) {', '.join(_supersonic_notes)} "
                    "below the sonic limit of its area: its total pressure is less than "
                    "the least any flow of its mass flux has at its area (the sonic "
                    "Pt*), so the chamber has no subsonic state, and none is a choked "
                    "duct's under-expanded exit (a choked-bounded junction: one choked "
                    "compressible channel of no larger area feeding the chamber, and no "
                    "element reading its static pressure). Usually the element downstream "
                    "is no restriction for that flow -- check its geometry -- or the "
                    "chamber's area is smaller than the flow it is given needs."
                )

        if not success:
            try:
                _hint = self._declared_direction_hints(final_x)
            except Exception:  # noqa: BLE001 -- a diagnostic never breaks a solve
                _hint = ""
            if _hint:
                message = f"{message} {_hint}"
            warnings.warn(
                f"NetworkSolver did not converge: {message} "
                f"(|F|={final_norm:.3e}). "
                "Returning best iterate.",
                stacklevel=2,
            )

        # Store solution for warm-start reuse
        self._last_x = final_x.copy()

        # Decompose final residual vector into per-unknown contributions.
        # One extra _residuals() call; same cost as a single Newton step.
        # Guard with try-except: mocked or patched _residuals_and_jacobian in
        # tests may raise, and diagnostic collection must never break the caller.
        try:
            _res_final = self._residuals(final_x)
            _breakdown = {
                name: float(abs(r)) for name, r in zip(self.unknown_names, _res_final, strict=False)
            }
            _worst = sorted(_breakdown.items(), key=lambda kv: kv[1], reverse=True)[:10]
        except Exception:
            _breakdown = {}
            _worst = []

        # Store diagnostic data
        self._diagnostic_data = {
            "x0": x0_use.copy(),
            "solution": final_x.copy(),
            "converged": success,
            "final_norm": final_norm,
            "message": message,
            "unknown_names": self.unknown_names.copy(),
            "solver_method": method,
            "solver_options": options.copy() if options else {},
            "wall_time": getattr(self, "_wall_time", None),
            "convergence_history": _history,
            "residual_breakdown": _breakdown,
            "worst_residuals": [{"name": n, "residual": v} for n, v in _worst],
            "lm_started_at_eval": _lm_started_at_eval[0],
            "stall_handover": _stall_fired[0],
            "lm_stall": _lm_stall_fired[0],
            "solver_settings_used": {
                "method": method,
                "init_strategy": init_strategy,
                "timeout": timeout,
                "stall_handover": self._stall_handover,
                **options,
            },
        }

        sol_dict = dict(zip(self.unknown_names, final_x, strict=False))
        # A choked-bounded junction reports as the plenum its root is (P = Pt);
        # its P unknown is the closure's continuation past the sonic point
        # (#788), which no other row reads. __x_solution__ keeps the unknown.
        for nid in _choked_bounded:
            sol_dict[f"{nid}.P"] = float(self.network.nodes[nid]._choked_bounded["P"])

        # Ensure _derived_states matches final_x (best iterate), not the
        # last trial point the solver happened to evaluate.  This is cheap:
        # forward T/Y propagation only, no residuals or Jacobian.
        self._boundary_arrivals = {}
        self._record_boundary_arrivals = True
        try:
            self._propagate_states(final_x)
        finally:
            self._record_boundary_arrivals = False
        from combaero import mass_to_mole

        station_pt = self._report_moving_nodes(final_x, sol_dict)
        for nid, (T, Y, _) in self._derived_states.items():
            # The node's stagnation T; ".T" is its static T (#782), set below once
            # the state is evaluable (a node with no gas state keeps T = Tt).
            sol_dict[f"{nid}.Tt"] = T
            sol_dict[f"{nid}.T"] = T
            for i, yi in enumerate(Y):
                sol_dict[f"{nid}.Y[{i}]"] = float(yi)

            # Store mole fractions as well (cheap, prevents redundant conversions)
            X = mass_to_mole(np.array(Y))
            for i, xi in enumerate(X):
                sol_dict[f"{nid}.X[{i}]"] = float(xi)

            # Node-specific diagnostics (e.g. mach, transport/thermo properties).
            # A node at or below 0 K has no thermodynamic state to report:
            # flagged, not evaluated. Any other node a failed solve leaves
            # unevaluable reports nothing, as the elements below do; on a
            # successful solve that is a bug and raises.
            node = self.network.nodes[nid]
            if nid in self._T_no_gas:
                sol_dict[f"{nid}.no_gas_state"] = 1.0
                continue
            try:
                state = self._report_state(node, final_x)
                sol_dict[f"{nid}.mach"] = node.mach(state)
                diag = dict(node.diagnostics(state)) if hasattr(node, "diagnostics") else {}
                # The static T: a momentum chamber's is the one its closure
                # solved (its diagnostics' "T_static", reported as "<node>.T"
                # only; the "<chamber>.T_static" alias went in 0.9.0); any
                # other node's is the state's (Tt at a plenum or
                # incompressible node).
                T_static = diag.pop("T_static", None)
                for key, val in diag.items():
                    sol_dict[f"{nid}.{key}"] = val
                sol_dict[f"{nid}.T"] = float(state.T if T_static is None else T_static)
                if nid in station_pt:
                    sol_dict[f"{nid}.Pt"], sol_dict[f"{nid}.mach"] = station_pt[nid]
            except (ValueError, RuntimeError) as exc:
                if success:
                    raise
                self._record_diagnostics_exception(exc, f"node '{nid}'")

        # Element-specific diagnostics (e.g. throat Mach, P-ratio)
        from .components import (
            MultiPortChamberBase as _MultiPortChamberBase,
        )
        from .components import (
            TeeJunctionElement as _TeeJunctionElement,
        )

        sol_dict["__element_diag__"] = {}
        for eid, element in self.network.elements.items():
            try:
                m_indices = self._unknown_indices.get(eid, [])
                if isinstance(element, _TeeJunctionElement):
                    nodes = self.network.nodes
                    state_com = self._report_state(nodes[element.common_node], final_x)
                    state_straight = self._report_state(nodes[element.straight_node], final_x)
                    state_branch = self._report_state(nodes[element.branch_node], final_x)
                    if m_indices:
                        state_com.m_dot = float(final_x[m_indices[0]])
                        state_branch.m_dot = float(final_x[m_indices[1]])
                        state_straight.m_dot = float(final_x[m_indices[0]]) - float(
                            final_x[m_indices[1]]
                        )
                    diag = element.diagnostics(state_com, state_straight, state_branch)
                elif isinstance(element, _MultiPortChamberBase):
                    nodes = self.network.nodes
                    port_states = [
                        self._report_state(nodes[pid], final_x) for pid in element.port_nodes
                    ]
                    P_jct_val = float(final_x[m_indices[0]]) if m_indices else 0.0
                    # Per-port mdots in junction convention (positive = out of
                    # junction), same sign-mapping the residual call uses.
                    port_mdots: list[float] = []
                    for i, outer_id in enumerate(element._port_element_ids):
                        outer_indices = self._unknown_indices.get(outer_id, [])
                        outer_mdot = float(final_x[outer_indices[0]]) if outer_indices else 0.0
                        port_mdots.append(element._port_signs[i] * outer_mdot)
                    diag = element.diagnostics(port_states, P_jct_val, port_mdots)
                else:
                    state_in = self._report_state(self.network.nodes[element.from_node], final_x)
                    state_out = self._report_state(self.network.nodes[element.to_node], final_x)
                    # Inject solved m_dot so velocity-dependent diagnostics (e.g. Mach) are correct
                    if m_indices:
                        m_solved = float(final_x[m_indices[0]])
                        state_in.m_dot = m_solved
                        state_out.m_dot = m_solved
                    # An end at a PressureBoundary that RECEIVES this element's
                    # flow reports the gas arriving, not the boundary's set Tt/Y.
                    for end_id, end_state in (
                        (element.to_node, state_out),
                        (element.from_node, state_in),
                    ):
                        arrival = self._boundary_arrivals.get((end_id, eid))
                        if arrival is not None:
                            end_state.Tt = arrival[0]
                            end_state.Y = arrival[1]
                    self._main_face(element, self.network.nodes[element.to_node], state_out)
                    flows = self._network_flows(element, final_x)
                    temps = self._network_temperatures(element, final_x, lagged=False)
                    if flows is None:
                        diag = element.diagnostics(state_in, state_out)
                    elif temps is None:
                        diag = element.diagnostics(state_in, state_out, flows=flows)
                    else:
                        diag = element.diagnostics(
                            state_in, state_out, flows=flows, temperatures=temps
                        )
            except (ValueError, RuntimeError) as exc:
                # A failed solve can end on a state no element can evaluate
                # (a negative pressure from a diverging iterate). Raising here
                # threw away the auto-retry that would have recovered it (#536:
                # a tee with port steps, recovered by the warm-start retry), so
                # a FAILED solve reports no diagnostics for that element
                # instead, and records why in __caught_exceptions__ (#810),
                # strict mode or not: a failed state is expected to be
                # unevaluable, as a rejected probe is. A converged state that
                # cannot be evaluated is a bug and still raises; the
                # admissibility checks above demote the non-physical roots
                # (#810: a negative node pressure) before this runs.
                if success:
                    raise
                self._record_diagnostics_exception(exc, f"element '{eid}'")
                diag = {}
            _station_ends(diag, element, station_pt)
            sol_dict["__element_diag__"][eid] = diag
            for key, val in diag.items():
                sol_dict[f"{eid}.{key}"] = val

        # Wall-specific thermodynamics. With coupling off the walls applied no
        # heat, so they report none: a Q computed here would be heat the
        # solution never saw (#481).

        walls = self.network.walls if self.network.thermal_coupling_enabled else {}
        for wid, wall in walls.items():
            if wall.one_sided:
                self._report_one_sided_wall(wid, wall, final_x, sol_dict)
                continue
            obj_a = self.network.elements.get(wall.element_a) or self.network.nodes.get(
                wall.element_a
            )
            obj_b = self.network.elements.get(wall.element_b) or self.network.nodes.get(
                wall.element_b
            )

            feed_a, tgt_a = self._flow_ends(obj_a, final_x)
            feed_b, tgt_b = self._flow_ends(obj_b, final_x)
            state_a = self._get_node_state(self.network.nodes[feed_a], final_x)
            state_b = self._get_node_state(self.network.nodes[feed_b], final_x)

            if wall.element_a in self.network.elements:
                m_indices_a = self._unknown_indices.get(wall.element_a, [])
                if m_indices_a:
                    state_a.m_dot = final_x[m_indices_a[0]]
            if wall.element_b in self.network.elements:
                m_indices_b = self._unknown_indices.get(wall.element_b, [])
                if m_indices_b:
                    state_b.m_dot = final_x[m_indices_b[0]]

            ch_result_a = self._wall_side_htc(wall, obj_a, state_a, final_x, lagged=False)
            ch_result_b = self._wall_side_htc(wall, obj_b, state_b, final_x, lagged=False)

            if ch_result_a and ch_result_b:
                # The same exchange as in the residual evaluation (#627)
                C_a, _, _, _, sign_a, _ = self._wall_capacity(
                    wall.element_a, state_a, float(state_b.Tt)
                )
                C_b, _, _, _, sign_b, _ = self._wall_capacity(
                    wall.element_b, state_b, float(state_a.Tt)
                )
                counterflow = wall.is_counterflow(sign_a, sign_b)
                wall_res = wall.compute_coupling(
                    ch_result_a.h,
                    ch_result_a.T_aw,
                    obj_a.surface.area,
                    ch_result_b.h,
                    ch_result_b.T_aw,
                    obj_b.surface.area,
                    C_a,
                    C_b,
                    counterflow,
                    **self._wall_gas(wall, state_a, state_b),
                )

                sol_dict[f"{wid}.Q"] = float(wall_res.Q)
                sol_dict[f"{wid}.T_hot"] = float(wall_res.T_hot)
                # The exchange: UA, NTU on the smaller capacity rate, and the
                # effectiveness Q / (C_min dT) (0 between two reservoirs).
                C_min = min(C_a, C_b)
                sol_dict[f"{wid}.UA"] = float(wall_res.UA)
                sol_dict[f"{wid}.NTU"] = (
                    float(wall_res.UA) / C_min if 0.0 < C_min < math.inf else 0.0
                )
                sol_dict[f"{wid}.effectiveness"] = float(wall_res.effectiveness)
                sol_dict[f"{wid}.counterflow"] = float(counterflow)

                # Heat injected into a node that cannot take it -- a boundary,
                # whose state is fixed -- leaves the network with the stream
                # crossing that boundary: the network is the control volume
                # and the far side is out of scope. Report it, so the balance
                # closes: the stream leaves with m h(T_upstream) + Q_wall_out.
                # Q > 0 flows A -> B, so A's target receives -Q and B's +Q.
                q_to_boundary = 0.0
                if tgt_a != tgt_b:
                    for tgt, q_in in ((tgt_a, -float(wall_res.Q)), (tgt_b, float(wall_res.Q))):
                        if tgt in self._wall_ebs:
                            continue
                        if self.network.nodes[tgt].is_flow_boundary:
                            q_to_boundary += q_in
                            key = f"{tgt}.Q_wall_out"
                        else:
                            # A node that takes no heat (a WallNode): no
                            # stream carries it anywhere, so it is withheld,
                            # not 'out' (#481).
                            key = f"{tgt}.Q_withheld"
                        sol_dict[key] = sol_dict.get(key, 0.0) + q_in
                sol_dict[f"{wid}.Q_to_boundary"] = q_to_boundary
                sol_dict[f"{wid}.h_a"] = float(ch_result_a.h)
                sol_dict[f"{wid}.h_b"] = float(ch_result_b.h)
                sol_dict[f"{wid}.T_aw_a"] = float(ch_result_a.T_aw)
                sol_dict[f"{wid}.T_aw_b"] = float(ch_result_b.T_aw)
                # The film on its side (#386): eta_bar, M, DR, envelope flag.
                for ch in (ch_result_a, ch_result_b):
                    for key, val in (getattr(ch, "film", None) or {}).items():
                        sol_dict[f"{wid}.{key}"] = float(val)
                # The radiation on its faces (#793): q_rad [W/m^2 of that
                # face], Q_rad [W], the face's surface T and eps_g.
                for key, val in (getattr(wall_res, "radiation", None) or {}).items():
                    sol_dict[f"{wid}.{key}"] = float(val)

                # Temperature profile [K], side A to side B: the one
                # compute_coupling just placed for these inputs (k(T) at its
                # fixed point, h floored as the coupling floors it).
                sol_dict[f"{wid}.T_interface"] = list(wall._last_profile or [])
                # 1 when every layer's k came from its material's sourced range,
                # 0 when any followed the k(T) envelope (#882).
                sol_dict[f"{wid}.material_in_range"] = float(
                    all(layer.in_range for layer in wall.layers)
                )

        # Heat an energy boundary gives a node that the node does not take up:
        # below a 1 mg/s floor the heat taken up fades to none at zero flow
        # (cb.mixer_heat_per_unit_flow), so only a stagnant node withholds any
        # (#481). And the heat a 'fraction' boundary applied.
        for nid, node in self.network.nodes.items():
            ebs = getattr(node, "energy_boundaries", None) or []
            if not ebs:
                continue
            for eb in ebs:
                if eb.is_target and eb.unknown_name in self._name_to_index:
                    eb.Q_solved = float(final_x[self._name_to_index[eb.unknown_name]])
            # An undetermined Tt_target heat (no through-flow) applies none.
            q_given = sum(float(eb.heat) for eb in ebs if eb.id not in _target_stagnant)
            frac = sum(float(eb.fraction) for eb in ebs)
            inflows = self._node_inflows(nid, final_x)
            m_in = sum(m for m, _, _ in inflows)
            # Always reported, 0.0 when nothing is withheld: which keys a
            # result has follows the network, never a solved value. Above the
            # floor phi = 1/m_in, and m_in * (1/m_in) can round to 1 - 1 ulp
            # (m_in = 7.3), so that remainder is snapped to zero rather than
            # reported as heat withheld.
            withheld = 0.0
            if q_given != 0.0:
                phi, _ = cb.mixer_heat_per_unit_flow(m_in)
                withheld = q_given * (1.0 - m_in * phi)
                if abs(withheld) <= _WITHHELD_RTOL * abs(q_given):
                    withheld = 0.0
            sol_dict[f"{nid}.Q_withheld"] = sol_dict.get(f"{nid}.Q_withheld", 0.0) + withheld
            if frac != 0.0 and not isinstance(node, CombustorNode) and m_in > 0.0:
                Y_mix = np.zeros(self._n_species)
                H = 0.0
                for m, Tt, Y in inflows:
                    X = cb.mass_to_mole(list(Y))
                    H += m * float(cb.h_mass(Tt, X))
                    Y_mix += m * np.asarray(Y, dtype=float)
                X_mix = cb.mass_to_mole(list(Y_mix / m_in))
                h_ref = float(cb.h_mass(cb.SENSIBLE_ENTHALPY_REF_T, X_mix))
                sol_dict[f"{nid}.Q_fraction"] = frac * (H - m_in * h_ref)
            elif frac != 0.0 and isinstance(node, CombustorNode) and m_in > 0.0:
                # On a combustor the fraction is of the PRODUCTS' sensible
                # enthalpy at T_ad (#481): h_out - h_ref = (1+f)(h_ad - h_ref)
                # + Q phi, so the heat it applied is f m (h_out - h_ref - Q
                # phi)/(1+f) from the solved products. Unreported, the ledger
                # saw an unexplained imbalance and demoted the root.
                st = self._get_node_state(node, final_x)
                X_b = list(st.X)
                h_out = float(cb.h_mass(float(st.Tt), X_b))
                h_ref = float(cb.h_mass(cb.SENSIBLE_ENTHALPY_REF_T, X_b))
                phi, _ = cb.mixer_heat_per_unit_flow(m_in)
                sol_dict[f"{nid}.Q_fraction"] = (
                    frac * m_in * (h_out - h_ref - q_given * phi) / (1.0 + frac)
                )
            # No fraction, or no flow for it to be a fraction of: none applied.
            sol_dict.setdefault(f"{nid}.Q_fraction", 0.0)
            # The heat each user boundary applied, in every mode (#792): the
            # declared Q, a Tt_target's SOLVED Q (its unknown at final_x, NaN
            # when no through-flow determines it), plus its share of the
            # fraction heat; summed per node as {node}.Q. > 0 heats the fluid.
            # A wall's heat is its wall's {wid}.Q.
            q_frac = float(sol_dict.get(f"{nid}.Q_fraction", 0.0))
            user_ebs = [eb for eb in ebs if eb.id != f"_wall_{nid}"]
            for eb in user_ebs:
                share = q_frac * float(eb.fraction) / frac if frac != 0.0 else 0.0
                undetermined = eb.id in _target_stagnant
                sol_dict[f"{eb.id}.Q"] = math.nan if undetermined else float(eb.heat) + share
            if user_ebs:
                sol_dict[f"{nid}.Q"] = float(sum(sol_dict[f"{eb.id}.Q"] for eb in user_ebs))

        # Second law per element and the network's control-volume ledger
        # (#545). Both demote, below. Evaluated where the root finder
        # converged; a diverged iterate has no balance worth reporting.
        verdicts: dict[str, dict[str, Any]] = {}
        wall_verdicts: dict[str, dict[str, Any]] = {}
        ledger: dict[str, Any] = {}
        if converged:
            from .second_law import SecondLawReport

            report = SecondLawReport(self, final_x, sol_dict)
            verdicts = report.element_verdicts()
            wall_verdicts = report.wall_verdicts()
            try:
                ledger = report.ledger()
            except (ValueError, RuntimeError) as exc:
                # a report, never a reason for the solve to raise
                ledger = {"error": str(exc)}
        for eid, verdict in verdicts.items():
            if "sigma_norm" in verdict:
                sol_dict[f"{eid}.sigma_norm"] = verdict["sigma_norm"]
                sol_dict["__element_diag__"].setdefault(eid, {})["sigma_norm"] = verdict[
                    "sigma_norm"
                ]
        sol_dict["__second_law__"] = verdicts
        sol_dict["__second_law_failures__"] = sorted(
            eid for eid, v in verdicts.items() if v["verdict"] == "fail"
        )
        # Each ThermalWall as its own control volume (#551), kept apart from the
        # per-element verdicts the GUI attributes to canvas elements.
        sol_dict["__second_law_walls__"] = wall_verdicts
        sol_dict["__second_law_wall_failures__"] = sorted(
            wid for wid, v in wall_verdicts.items() if v["verdict"] == "fail"
        )
        sol_dict["__ledger__"] = ledger

        # A root whose control volume does not close in mass or energy is not
        # a solution of the network, whatever |F| says (#552). The bound sits
        # 9x above the worst closure of any accepted solve in the suite.
        if success and ledger and "error" not in ledger:
            from .second_law import LEDGER_CLOSURE_TOL

            _rel = {kind: ledger.get(f"{kind}_closure_rel") for kind in ("mass", "energy")}
            _open = {
                kind: float(v)
                for kind, v in _rel.items()
                if v is not None and not abs(float(v)) <= LEDGER_CLOSURE_TOL  # NaN is open
            }
            if _open:
                success = False
                consistent = False
                outcome = SolveOutcome.INCONSISTENT
                message = (
                    "Converged to a root whose network balance does not close: "
                    + ", ".join(f"{k} {v:+.2e} relative" for k, v in _open.items())
                    + f" (bound {LEDGER_CLOSURE_TOL:.0e}). See __ledger__."
                )
                warnings.warn(f"NetworkSolver: {message}", stacklevel=2)

        # A root on which some element destroys entropy -- manufactures flow
        # work on its most favourable heat path -- is not a physical answer
        # (#598). Demoted like an open ledger. On the stress sweep at the
        # switch, the roots this newly refuses were an ejector entraining more
        # power than its primary jet carries and a lossless (loss_factor = 0)
        # ideal-tee split at Mach, both outside their models. The walls' own
        # check stays a report: it was not part of that audit.
        # Not inside the adiabatic first stage of a heated-channel solve
        # (_solve_rayleigh_staged): its root leaves the heat out of the
        # momentum on purpose -- the very sigma < 0 of #551 -- and is only the
        # seed of the exact stage, which is judged in full. Refusing it sent
        # that stage through the whole retry chain (a grid example 3.5 -> 14 s).
        _seed_stage = any(
            getattr(el, "_rayleigh_suspended", False) for el in self.network.elements.values()
        )
        _sigma_bad = list(sol_dict["__second_law_failures__"])
        if success and _sigma_bad and not _seed_stage:
            from .second_law import SIGMA_TOL

            success = False
            consistent = False
            outcome = SolveOutcome.INCONSISTENT
            _worst_sigma = {
                eid: float(verdicts[eid].get("sigma_norm_high", float("nan"))) for eid in _sigma_bad
            }
            message = (
                "Converged to a root that violates the second law at "
                + ", ".join(f"{eid} (sigma_norm_high {s:+.2e})" for eid, s in _worst_sigma.items())
                + f": the element creates flow work (bound -{SIGMA_TOL:.0e}). "
                "See __second_law__."
            )
            # An element that knows why its own closure creates flow work at
            # this root names the cause (an ideal tee below f ~ M^2/4).
            for eid in _sigma_bad:
                hint = getattr(self.network.elements.get(eid), "second_law_hint", None)
                text = hint(sol_dict) if callable(hint) else ""
                if text:
                    message = f"{message} {text}"
            warnings.warn(f"NetworkSolver: {message}", stacklevel=2)

        # Compressors and turbines (#893): a root with the flow reversed is a
        # surging compressor or a turbine driven backwards, which the fixed
        # pressure ratio cannot describe: refused, never returned.
        for eid, el in self.network.elements.items():
            if not isinstance(el, TurbomachineElement):
                continue
            m_tm = float(sol_dict.get(f"{eid}.m_dot", 0.0))
            if converged and m_tm < -FLOW_REVERSAL_M_SCALE:
                # The cause first: the ledger and second law, which price the
                # machine forward, may already have demoted this root.
                kind = "compressor" if not el.expansion else "turbine"
                why = (
                    f"{type(el).__name__} '{eid}': the flow through it is reversed "
                    f"({m_tm:.4g} kg/s); a {kind} at a fixed pressure ratio "
                    f"{el.pressure_ratio:g} has no reversed operating point "
                    "(surge, or driven backwards). Check what drives the flow."
                )
                message = why if success else f"{why} {message}"
                success = False
                consistent = False
                outcome = SolveOutcome.INCONSISTENT
                warnings.warn(f"NetworkSolver: {why}", stacklevel=2)

        # Fuel controllers (#892): report the fuel flow and the equivalence ratio
        # the combustor reached, and refuse a root on the other side of
        # stoichiometric -- a target has a lean and a rich root, and the one
        # asked for is the answer, never the other.
        for fid, fnode in self.network.nodes.items():
            if not isinstance(fnode, FuelBoundary):
                continue
            phi, T_without = self._fuel_reactant_view(fnode, final_x, sol_dict)
            sol_dict[f"{fid}.m_dot"] = float(fnode.m_dot)
            sol_dict[f"{fid}.phi"] = phi
            sol_dict[f"{fid}.Tt_target"] = float(fnode.Tt_target)
            sol_dict[f"{fid}.Tt_without_fuel"] = T_without
            sol_dict[f"{fid}.rich_branch"] = 1.0 if fnode.branch == "rich" else 0.0
            other = (fnode.branch == "lean" and phi > 1.0 + 1e-6) or (
                fnode.branch == "rich" and phi < 1.0 - 1e-6
            )
            if success and not other:
                streams, _ = self._fuel_inflows(fnode, final_x, sol_dict)
                if sum(m for m, _, _ in streams) <= _FUEL_MIN_OXIDIZER_FLOW:
                    success = False
                    consistent = False
                    outcome = SolveOutcome.INCONSISTENT
                    message = (
                        f"FuelBoundary '{fid}': the combustor {fnode.combustor!r} gets no "
                        "oxidizer flow at this root: the fuel has nothing to burn with, and "
                        "the target is met only by vanishing flows."
                    )
                    warnings.warn(f"NetworkSolver: {message}", stacklevel=2)
            if success and other:
                success = False
                consistent = False
                outcome = SolveOutcome.INCONSISTENT
                message = (
                    f"FuelBoundary '{fid}': the root found is {'rich' if phi > 1.0 else 'lean'} "
                    f"(phi {phi:.3f}), not the {fnode.branch} one asked for. Give "
                    "m_dot_guess on the requested side, or check that the target is "
                    "reachable on it."
                )
                note = getattr(fnode, "_guess_note", "")
                if note:
                    message = f"{message} {note}"
                warnings.warn(f"NetworkSolver: {message}", stacklevel=2)
            if not success:
                streams, Y_f = self._fuel_inflows(fnode, final_x, sol_dict)
                m_ox = sum(m for m, _, _ in streams)
                where = ""
                start_ox = getattr(fnode, "_start_oxidizer", [])
                if m_ox <= max(_FUEL_MIN_OXIDIZER_FLOW, 1e-6 * float(fnode.m_dot)) and start_ox:
                    # No oxidizer at this point, but the network has one: the
                    # point is far from closing its mass rows. Diagnose with
                    # the oxidizer at the start state.
                    streams, where = start_ox, " (oxidizer at the start state)"
                    m_ox = sum(m for m, _, _ in streams)
                    try:
                        T_without = self._fuel_T_without(streams)
                    except (ValueError, RuntimeError):
                        T_without = float("nan")
                if m_ox <= max(_FUEL_MIN_OXIDIZER_FLOW, 1e-6 * float(fnode.m_dot)):
                    message = (
                        f"{message} FuelBoundary '{fid}': the combustor {fnode.combustor!r} gets no "
                        "oxidizer flow at this point: the fuel has nothing to burn with."
                    )
                else:
                    reason = self._fuel_window_reason(fnode, streams, Y_f, T_without, where)
                    if reason:
                        message = f"{message} {reason}"

        # An effusion plate with no cooling path (#998): no coolant through its
        # holes and no gas flow on its other side, so nothing takes heat from
        # it and its wall temperature is undetermined (NaN in its keys). The
        # flow root is fine; the plate's answer is not one. Refused like a
        # Tt_target with no through-flow (#792): a converged solve is demoted,
        # a failed one says why. What would model it: the comment at
        # EffusionPlateElement.wall_heat_transfer.
        _wall_undetermined = sorted(
            eid
            for eid, diag in (sol_dict.get("__element_diag__") or {}).items()
            if float(diag.get("no_cooling_path", 0.0)) == 1.0
        )
        if _wall_undetermined:
            note = (
                "EffusionPlateElement "
                + ", ".join(f"'{eid}'" for eid in _wall_undetermined)
                + ": no cooling path: wall T undetermined (heats until failure). No "
                "coolant flows through its holes and its gas side has no flow (a "
                "plenum, or a chamber with none), so nothing takes heat from the "
                "plate; radiation, backside convection and conduction to its mounts "
                "are not modelled. Drive coolant through it, or discharge it into a "
                "flowing chamber."
            )
            if success:
                success = False
                consistent = False
                outcome = SolveOutcome.INCONSISTENT
                message = f"Converged, but {note}"
                warnings.warn(f"NetworkSolver: {message}", stacklevel=2)
            else:
                message = f"{message} {note}"
        # Reported only (#998): plates with a face above their material's
        # limit -- in hot gas with no coolant, the adiabatic T_aw.
        _wall_above_limit = sorted(
            eid
            for eid, diag in (sol_dict.get("__element_diag__") or {}).items()
            if float(diag.get("wall_above_limit", 0.0)) == 1.0
        )

        # Incompressible validity (#588): each element on an incompressible
        # closure against one Mach limit. REPORTED only.
        validity: dict[str, dict[str, Any]] = {}
        if converged:
            from .compressibility import validity_verdicts

            validity = validity_verdicts(
                self.network,
                sol_dict.get("__element_diag__") or {},
                float(self.incompressible_mach_limit),
                {nid: sol_dict.get(f"{nid}.mach") for nid in self.network.nodes},
            )
        sol_dict["__incompressible_validity__"] = validity
        sol_dict["__incompressible_warnings__"] = sorted(
            eid for eid, v in validity.items() if v["verdict"] == "warn"
        )

        sol_dict["__complete_states__"] = self.extract_complete_states(sol_dict)
        sol_dict["__success__"] = success
        sol_dict["__message__"] = message
        sol_dict["__final_norm__"] = final_norm
        # Why the solve ended, split into the two questions ``__success__``
        # conflates: did the root finder get there, and is what it found
        # physically admissible. ``__success__`` keeps its meaning -- both.
        sol_dict["__converged__"] = bool(converged)
        sol_dict["__consistent__"] = consistent
        sol_dict["__inconsistent_elements__"] = list(_bad_junctions)
        sol_dict["__P_nonpositive__"] = list(_P_nonpositive)
        sol_dict["__T_clamped__"] = list(_T_clamped)
        sol_dict["__heat_faded__"] = list(_heat_faded)
        sol_dict["__heat_undetermined__"] = list(_target_stagnant)
        sol_dict["__wall_undetermined__"] = _wall_undetermined
        sol_dict["__wall_above_limit__"] = _wall_above_limit
        sol_dict["__supersonic_chambers__"] = list(_supersonic)
        sol_dict["__choked_bounded_junctions__"] = list(_choked_bounded)
        sol_dict["__outcome__"] = outcome
        # The rows carrying the residual, largest first. Already computed for
        # the diagnostic payload; a failed solve is far easier to read with
        # them than with a single norm.
        sol_dict["__worst_residuals__"] = [{"name": n, "residual": v} for n, v in _worst]
        sol_dict["__unknown_names__"] = list(self.unknown_names)
        sol_dict["__x_solution__"] = list(final_x)
        # Convergence history: per-eval (eval_idx, t_elapsed_s, residual_norm)
        # triples accumulated by the residuals wrapper. Throttled at 0.1s
        # intervals (or whenever the iterate improves), capped at 1000 points.
        # Used by GUI diagnostics to plot the |F| vs iteration trace -- helps
        # diagnose slow-converging networks ("Newton overshoots 5 steps then
        # damps" vs "monotone but many iterations" tell different stories).
        sol_dict["__convergence_history__"] = list(_history)
        return sol_dict

    def extract_complete_states(self, result: dict[str, float] | None = None) -> dict[str, Any]:
        """
        Extract CompleteState for all network nodes.

        Leverages the existing C++ CompleteState structure for maximum accuracy
        and performance. By default, it uses the solver's internal cached derived
        states populated during solve().

        Args:
            result: Optional dictionary containing node T, P, X and Y. If omitted,
                uses internal solver state from the last solve() call.

        Returns:
            Dictionary mapping node IDs to CompleteState objects.
        """
        import combaero as cb
        from combaero import mass_to_mole

        complete_states = {}

        if result is None:
            # OPTIMIZED: Use cached derived states from _propagate_states
            # Note: This is primarily for manual diagnostic exploration after solve()
            if not self._derived_states:
                return {}

            for nid, (T, Y, _) in self._derived_states.items():
                # For P, we rely on node.P or total P if defined in the network
                node = self.network.nodes[nid]
                P = getattr(node, "P", getattr(node, "Pt", 101325.0))
                X = mass_to_mole(np.array(Y))
                complete_states[nid] = cb.complete_state(T, P, X)
            return complete_states

        # DICT PATH: Consumption from solve() result (as used in the GUI)
        # Try to get composition from first node that has data to use as default
        default_X = None
        for node_id in self.network.nodes:
            # Prefer pre-computed X if available
            if f"{node_id}.X[0]" in result:
                X_list = []
                idx = 0
                while f"{node_id}.X[{idx}]" in result:
                    X_list.append(result[f"{node_id}.X[{idx}]"])
                    idx += 1
                if X_list:
                    default_X = np.array(X_list)
                    break
            elif f"{node_id}.Y[0]" in result:
                Y_list = []
                idx = 0
                while f"{node_id}.Y[{idx}]" in result:
                    Y_list.append(result[f"{node_id}.Y[{idx}]"])
                    idx += 1
                if Y_list:
                    default_X = mass_to_mole(np.array(Y_list))
                    break

        if default_X is None:
            default_X = cb.humid_air_composition(350.0, 101325.0, 0.0)

        for node_id in self.network.nodes:
            try:
                T = result.get(f"{node_id}.T", result.get(f"{node_id}.Tt", 300.0))
                P = result.get(f"{node_id}.P", result.get(f"{node_id}.Pt", 101325.0))

                # Use pre-computed X if available, otherwise fallback to Y-parsing
                X = None
                if f"{node_id}.X[0]" in result:
                    X_list = []
                    idx = 0
                    while f"{node_id}.X[{idx}]" in result:
                        X_list.append(result[f"{node_id}.X[{idx}]"])
                        idx += 1
                    if X_list:
                        X = np.array(X_list)

                if X is None and f"{node_id}.Y[0]" in result:
                    Y_list = []
                    idx = 0
                    while f"{node_id}.Y[{idx}]" in result:
                        Y_list.append(result[f"{node_id}.Y[{idx}]"])
                        idx += 1
                    if Y_list:
                        X = mass_to_mole(np.array(Y_list))

                X = X if X is not None else default_X

                if not (math.isnan(T) or math.isnan(P)):
                    complete_states[node_id] = cb.complete_state(T, P, X)
                else:
                    complete_states[node_id] = None
            except Exception:
                complete_states[node_id] = None

        return complete_states

    def get_diagnostic_report(self) -> dict[str, Any]:
        """
        Generate comprehensive diagnostic report comparing initial guess to final solution.

        Returns:
            Dictionary containing detailed analysis of solver behavior including
            x0 vs solution comparison, convergence metrics, and per-variable analysis.
        """
        if not hasattr(self, "_diagnostic_data"):
            raise RuntimeError("No diagnostic data available. Run solve() first.")

        data = self._diagnostic_data
        x0 = data["x0"]
        solution = data["solution"]
        unknown_names = data["unknown_names"]

        # Calculate differences and statistics
        differences = solution - x0
        abs_differences = np.abs(differences)
        rel_differences = np.abs(differences / (x0 + 1e-10))  # Avoid division by zero

        # Per-variable analysis
        variable_analysis = {}
        for i, name in enumerate(unknown_names):
            variable_analysis[name] = {
                "x0": float(x0[i]),
                "solution": float(solution[i]),
                "difference": float(differences[i]),
                "abs_difference": float(abs_differences[i]),
                "rel_difference": float(rel_differences[i]),
            }

        # Overall statistics
        overall_stats = {
            "max_abs_difference": float(np.max(abs_differences)),
            "mean_abs_difference": float(np.mean(abs_differences)),
            "median_abs_difference": float(np.median(abs_differences)),
            "max_rel_difference": float(np.max(rel_differences)),
            "mean_rel_difference": float(np.mean(rel_differences)),
            "median_rel_difference": float(np.median(rel_differences)),
            "std_difference": float(np.std(differences)),
        }

        # Convergence analysis
        convergence_analysis = {
            "converged": data["converged"],
            "final_norm": data["final_norm"],
            "message": data["message"],
            "wall_time": data["wall_time"],
            "solver_method": data["solver_method"],
            "solver_options": data["solver_options"],
            "num_unknowns": len(unknown_names),
        }

        # Identify most changed variables
        most_changed_indices = np.argsort(abs_differences)[-5:]  # Top 5
        most_changed = [
            {
                "name": unknown_names[i],
                "abs_difference": float(abs_differences[i]),
                "rel_difference": float(rel_differences[i]),
                "x0": float(x0[i]),
                "solution": float(solution[i]),
            }
            for i in reversed(most_changed_indices)
        ]

        return {
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "variable_analysis": variable_analysis,
            "overall_statistics": overall_stats,
            "convergence_analysis": convergence_analysis,
            "most_changed_variables": most_changed,
            "unknown_names": unknown_names,
        }
