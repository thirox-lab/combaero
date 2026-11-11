"""
Pressure-loss correlation classes for CombustorNode.

Two families of correlations are provided:

**Fraction family** -- loss expressed as a fraction of inlet total pressure:
    dPt = xi * Pt_in        =>   Pt_out = Pt_in * (1 - xi)
  Classes: ConstantFractionLoss, LinearThetaFractionLoss

**Head family** -- loss expressed as a multiple of the inlet dynamic head:
    dPt = zeta * 0.5 * rho_in * v_in^2, rho_in and v_in at the inlet
    STAGNATION state (Tt_in, Pt_in): Lefebvre's pressure-loss factor, a
    total-total low-Mach form (#782 audit; a coupled Fanno-Rayleigh form is #830)
  Classes: ConstantHeadLoss, LinearThetaHeadLoss  (Phase 1b, needs mdot_total in ctx)

All callables return ``(loss_value, d_loss_value_d_theta)`` where
``loss_value`` is the *fractional* total-pressure loss (dimensionless)
consumed by ``Pt_out = Pt_in * (1 - loss_value)``.

WHERE THE LOSS IS EVALUATED. The four built-in classes are C++ models: each
resolves to a ``cb.PressureLossModel`` (``pressure_loss_model()``) that one C++
call, ``cb.pressure_loss_compose``, evaluates with every partial (m_dot, Tt_in,
Pt_in, each inlet mass fraction, and each generic model input, theta first;
#725). ``PressureLossElement`` calls it directly, and ``evaluate(ctx, Y_in)``
returns the same ``PressureLossResult``. The C++
``combustor_residuals_and_jacobians`` also accepts a loss callable, but
``CombustorNode`` passes it a zero loss, so that chain does not run in a
network.

A USER CALLABLE STAYS PYTHON (#718). Any object returning ``(xi, dxi/dtheta)``
from a context still works. Its other dependences are central differences,
which are what let it be arbitrary Python. That interface is the prototype for
user-supplied loss maps (e.g. a scipy interpolant over theta and fuel splits;
the input design for the splits is open, #725).

Theta = T_ad / Tt_in - 1  (dimensionless temperature-rise ratio).

THE CONTEXT (``PressureLossContext``) carries the inlet's stagnation state as
``Tt_in``/``Pt_in``. Until 0.8 these were ``T_in``/``P_in`` (``P_in`` the
inlet's STATIC pressure); the old names were removed in 0.9.0.
"""

from dataclasses import dataclass
from types import SimpleNamespace
from typing import Protocol, runtime_checkable

from .._core import (
    PressureLossInputs,
    PressureLossKind,
    PressureLossModel,
    PressureLossResult,
    mole_to_mass,
    pressure_loss_compose,
)

# PressureLossModel(kind, k, base, m_ref, area), positional: it is built on
# every evaluation, and pybind11's keyword parsing doubles its cost.
_FRACTION = PressureLossKind.Fraction
_HEAD = PressureLossKind.Head


@runtime_checkable
class PressureLossCallable(Protocol):
    """Protocol for pressure-loss correlation objects.

    The context object passed at solve time is a ``PressureLossContext``.
    """

    def __call__(self, ctx: object) -> tuple[float, float]: ...


class PressureLossContext(SimpleNamespace):
    """What a loss correlation is evaluated at: ``theta``, ``mdot_total``,
    the inlet stagnation state ``Tt_in``/``Pt_in``, ``X_in``, ``T_ad``, ...

    ``T_in``/``P_in`` were the 0.7 names (``P_in`` was the inlet's static
    pressure), removed in 0.9.0.
    """


def _ctx_inlet(ctx: object) -> tuple[float, float]:
    """(Tt_in, Pt_in) of a context. A hand-built 0.7 context with only
    ``T_in``/``P_in`` is refused with a message naming the new fields."""
    if not (hasattr(ctx, "Tt_in") and hasattr(ctx, "Pt_in")):
        raise AttributeError(
            "a pressure-loss context needs Tt_in and Pt_in (the inlet stagnation "
            "state); T_in/P_in were removed in 0.9.0 (#782)"
        )
    return float(ctx.Tt_in), float(ctx.Pt_in)  # type: ignore[attr-defined]


# ---------------------------------------------------------------------------
# Fraction family:  dPt/Pt_in = xi   (total-pressure ratio convention)
# ---------------------------------------------------------------------------


def _m_ref_arg(m_ref: float | None) -> float:
    """The design flow for C++ (0 = no fade), validated.

    A fraction loss is a head loss quoted at its design point: below it the
    drop scales with the dynamic head (``1 - (1 - r^2)^3``, C2 at r = 1), and
    at zero flow it vanishes. Without this a constant fraction holds a
    pressure difference at zero flow, and a reversing flow must jump through
    it -- the residual steps by xi Pt where the flow changes sign, which the
    solver cannot cross (#550).
    """
    if m_ref is None:
        return 0.0
    if m_ref <= 0.0:
        raise ValueError(f"m_ref must be positive, got {m_ref}")
    return float(m_ref)


def _evaluate_fraction(
    model: PressureLossModel, ctx: object, reads_theta: bool
) -> PressureLossResult:
    """A fraction model at a context: theta (when the class reads it) and,
    below a design flow only, m_dot; a context may omit the rest."""
    x = [float(ctx.theta)] if reads_theta else []  # type: ignore[attr-defined]
    m = float(ctx.mdot_total) if model.m_ref > 0.0 else 0.0  # type: ignore[attr-defined]
    return pressure_loss_compose(model, PressureLossInputs(mdot=m, x=x))


def _evaluate_head(
    model: PressureLossModel, ctx: object, Y_in: list[float] | None
) -> PressureLossResult:
    """A head model at a context. Without ``Y_in`` the composition is the
    context's ``X_in``, converted (the density the value reads either way)."""
    Tt_in, Pt_in = _ctx_inlet(ctx)
    Y = list(Y_in) if Y_in is not None else list(mole_to_mass(list(ctx.X_in)))  # type: ignore[attr-defined]
    inp = PressureLossInputs(
        mdot=float(ctx.mdot_total),  # type: ignore[attr-defined]
        Tt_in=Tt_in,
        Pt_in=Pt_in,
        Y_in=Y,
        x=[float(ctx.theta)],  # type: ignore[attr-defined]
    )
    return pressure_loss_compose(model, inp)


@dataclass
class ConstantFractionLoss:
    """Constant fractional total-pressure loss.

    Formula:  dPt/Pt_in = xi   =>   Pt_out = Pt_in * (1 - xi)

    Parameters
    ----------
    xi : float
        Fractional total-pressure loss [-].  E.g. 0.03 = 3% of Pt_in.
    m_ref : float | None
        Design mass flow [kg/s]. Below it the loss scales with the dynamic
        head and vanishes at zero flow (``_m_ref_arg``); at and above
        it the loss is exactly ``xi``. Set it wherever the flow may run low
        or reverse. ``None`` keeps a constant fraction at every flow.

    Notes
    -----
    ``d(dP/P)/dTheta = 0`` -- loss is independent of temperature rise,
    so the pressure Jacobian contribution is zero.
    """

    xi: float = 0.03
    m_ref: float | None = None

    def pressure_loss_model(self) -> PressureLossModel:
        """The C++ model (``cb.pressure_loss_compose``)."""
        return PressureLossModel(_FRACTION, 0.0, self.xi, _m_ref_arg(self.m_ref), 0.0)

    def evaluate(self, ctx: object, Y_in: list[float] | None = None) -> PressureLossResult:
        """xi and its partials, from C++."""
        return _evaluate_fraction(self.pressure_loss_model(), ctx, reads_theta=False)

    def __call__(self, ctx: object) -> tuple[float, float]:
        r = self.evaluate(ctx)
        return r.xi, r.dxi_dtheta


@dataclass
class LinearThetaFractionLoss:
    """Total-pressure loss fraction linear in temperature-rise ratio Theta.

    Formula:  dPt/Pt_in = k * Theta + xi0
              Theta   = T_ad / Tt_in - 1

    Parameters
    ----------
    k : float
        Sensitivity of fractional loss to temperature rise [-/-].
        Positive k means more heat release = higher pressure loss.
    xi0 : float
        Cold-flow base loss fraction [-] (value at Theta = 0).
    m_ref : float | None
        Design mass flow [kg/s], as for :class:`ConstantFractionLoss`.

    Notes
    -----
    Returns the analytical derivative ``d(dP/P)/dTheta = k``, which the
    ``PressureLossElement`` Jacobian uses directly for its theta column.
    """

    k: float = 0.5
    xi0: float = 0.02
    m_ref: float | None = None

    def pressure_loss_model(self) -> PressureLossModel:
        """The C++ model (``cb.pressure_loss_compose``)."""
        return PressureLossModel(_FRACTION, self.k, self.xi0, _m_ref_arg(self.m_ref), 0.0)

    def evaluate(self, ctx: object, Y_in: list[float] | None = None) -> PressureLossResult:
        """xi and its partials, from C++."""
        return _evaluate_fraction(self.pressure_loss_model(), ctx, reads_theta=True)

    def __call__(self, ctx: object) -> tuple[float, float]:
        r = self.evaluate(ctx)
        return r.xi, r.dxi_dtheta


# ---------------------------------------------------------------------------
# Head family:  dP = zeta * 0.5 * rho_in * v_in^2  (Euler loss coefficient)
#
# v_in = mdot_total / (rho_in * area), rho_in the ideal-gas density at
# (Tt_in, Pt_in, inlet composition): Lefebvre's total-total form. The return
# value is still a *fraction* of Pt_in, xi = zeta * 0.5 rho_in v_in^2 / Pt_in,
# so Pt_out = Pt_in (1 - xi) holds for both families. The partials in m_dot,
# Tt_in, Pt_in and composition come from C++ with the value
# (``cb.pressure_loss_compose`` over ``cb.combustor_head_loss``).
# ---------------------------------------------------------------------------


@dataclass
class ConstantHeadLoss:
    """Constant Euler head-loss coefficient.

    Formula:  dP = zeta * 0.5 * rho_in * v_in^2
              v_in = mdot_total / (rho_in * area)

    Parameters
    ----------
    zeta : float
        Euler loss coefficient [-].
    area : float
        Combustor inlet cross-sectional area [m^2].

    Notes
    -----
    Internally converts to a pressure fraction xi = dPt/Pt_in so both
    families share Pt_out = Pt_in (1 - xi).
    """

    zeta: float = 5.0
    area: float = 0.1

    def pressure_loss_model(self) -> PressureLossModel:
        """The C++ model (``cb.pressure_loss_compose``)."""
        return PressureLossModel(_HEAD, 0.0, self.zeta, 0.0, self.area)

    def evaluate(self, ctx: object, Y_in: list[float] | None = None) -> PressureLossResult:
        """xi and its partials, from C++."""
        return _evaluate_head(self.pressure_loss_model(), ctx, Y_in)

    def __call__(self, ctx: object) -> tuple[float, float]:
        r = self.evaluate(ctx)
        return r.xi, r.dxi_dtheta


@dataclass
class LinearThetaHeadLoss:
    """Euler head-loss coefficient linear in temperature-rise ratio Theta.

    Formula:  dP = (k * Theta + zeta0) * 0.5 * rho_in * v_in^2
              Theta = T_ad / Tt_in - 1
              v_in  = mdot_total / (rho_in * area)

    Parameters
    ----------
    k : float
        Sensitivity of loss coefficient to temperature rise [-/-].
    zeta0 : float
        Cold-flow Euler loss coefficient [-] (value at Theta = 0).
    area : float
        Combustor inlet cross-sectional area [m^2].

    Notes
    -----
    ``d(xi)/dTheta = k * q_frac``, q_frac = 0.5*rho_in*v_in^2 / Pt_in; the
    other partials come with it from ``evaluate``.
    """

    k: float = 1.0
    zeta0: float = 3.0
    area: float = 0.1

    def pressure_loss_model(self) -> PressureLossModel:
        """The C++ model (``cb.pressure_loss_compose``)."""
        return PressureLossModel(_HEAD, self.k, self.zeta0, 0.0, self.area)

    def evaluate(self, ctx: object, Y_in: list[float] | None = None) -> PressureLossResult:
        """xi and its partials, from C++."""
        return _evaluate_head(self.pressure_loss_model(), ctx, Y_in)

    def __call__(self, ctx: object) -> tuple[float, float]:
        r = self.evaluate(ctx)
        return r.xi, r.dxi_dtheta
