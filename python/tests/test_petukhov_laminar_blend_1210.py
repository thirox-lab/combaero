"""friction_model="petukhov" in a network is blended with Hagen-Poiseuille (#1210).

The raw fit f = (0.790 ln Re - 1.64)^-2 has a pole at Re 7.97: a channel's drop
spiked ~700x and turned non-monotone at low flow, a zero-flow seed stalled, and
the Fanno-Rayleigh march's needed drive jumped by up to 400x near 5e-4 of the
choked flux. The network kernel now blends it like haaland / serghides /
colebrook, with the laminar branch's own floor (Re 1e-6), so there is no
artificial floor on the drop either.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import combaero as cb
from combaero import _solver_tools
from combaero.network import (
    ChannelElement,
    CrossflowSegmentElement,
    FlowNetwork,
    NetworkSolver,
    PressureBoundary,
)

AIR = list(cb.species.dry_air())
Y_AIR = list(cb.mole_to_mass(AIR))


def _drop(m: float) -> tuple[float, float]:
    r = _solver_tools.channel_residuals_and_jacobian(
        m, 1e5, 1e5, 300.0, Y_AIR, 0.99e5, 1.0, 0.01, 1e-5, "petukhov", 1.0, False
    )
    return float(r.dP_calc), float(r.d_dP_d_mdot)


def test_the_low_flow_drop_is_linear_down_to_the_laminar_floor() -> None:
    """Hagen-Poiseuille from m = 1e-9 (Re ~ 7e-3) through the old pole (m ~ 1.15e-6)
    to Re ~ 700: no saturating floor, no spike."""
    X = cb.mass_to_mole(Y_AIR)
    rho = cb.density(300.0, 1e5, X)
    mu = cb.transport_state(300.0, 1e5, X).mu
    A = math.pi * 1e-4 / 4.0
    slope_hp = 32.0 * mu * 1.0 / (rho * A * 1e-4)
    for m in (1e-9, 1e-7, 1.15e-6, 1.2e-6, 3e-6, 1e-4):
        dP, dP_dm = _drop(m)
        assert dP == pytest.approx(slope_hp * m, rel=1e-4), m
        assert dP_dm == pytest.approx(slope_hp, rel=1e-4), m


def test_the_drop_rises_with_flow() -> None:
    prev = -1.0
    for m in np.geomspace(1e-12, 1e-1, 1500):
        dP, dP_dm = _drop(float(m))
        assert dP > prev, m
        assert dP_dm > 0.0, m
        prev = dP


@pytest.mark.parametrize("kind", ["channel", "segment"])
@pytest.mark.parametrize("dp", [1.0, 100.0])
def test_a_zero_flow_seed_converges(kind: str, dp: float) -> None:
    """Stalled at m ~ 7e-9 / 7e-7 after #1208 (J = 0 at the fit's f Re -> 0)."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=1e5 + dp, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("b", Pt=1e5, Tt=300.0, Y=Y_AIR))
    if kind == "channel":
        g.add_element(
            ChannelElement(
                "c",
                "a",
                "b",
                length=1.0,
                diameter=0.01,
                regime="incompressible",
                friction_model="petukhov",
            )
        )
    else:
        g.add_element(
            CrossflowSegmentElement(
                "c",
                "a",
                "b",
                length=1.0,
                area=math.pi * 1e-4 / 4.0,
                regime="incompressible",
                friction_model="petukhov",
            )
        )
    s = NetworkSolver(g)
    s._build_x0()
    r = s.solve(timeout=20, x0=s._zero_flow_seed())
    cold = NetworkSolver(g).solve(timeout=20)
    assert r["__success__"]
    assert cold["__success__"]
    assert r["c.m_dot"] == pytest.approx(cold["c.m_dot"], rel=1e-8)


def test_the_fanno_rayleigh_needed_drive_is_monotone_at_low_flux() -> None:
    """The #1204 review's duct: local Re 3.5 .. 11.5 crossed the pole at 5e-4 Gc."""
    Pt0, Tt0, Tte = 111759.76, 1390.05, 428.78
    L, D, rough, fm = 12.80e-3, 3.922e-3, 4.6e-6, 1.26
    Gc = cb._core.fanno_rayleigh_choked_mass_flux(Pt0, Tt0, Tte, AIR, L, D, rough, "petukhov", fm)
    prev = -1.0
    for frac in np.geomspace(1e-4, 1e-2, 120):
        r = cb._core.fanno_rayleigh_march(
            frac * Gc, Pt0, Tt0, Tte, AIR, L, D, rough, "petukhov", fm
        )
        drive = Pt0 - r.Pt_exit
        assert drive > prev, frac
        prev = drive
