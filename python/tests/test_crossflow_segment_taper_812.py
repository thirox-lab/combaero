"""A crossflow segment whose ends have their own sections (#812).

``CrossflowSegmentElement(area=A_from, area_to=A_to)``: a tapered annulus or
a growing liner. Each end's station half-drop, its joining jets' momentum and
its static-to-total closure read that end's section; friction runs on the mean
section. Compressible, the ends' closures carry the area change (Pt continuous
bar friction: isentropic); incompressible, where the node pressures are
statics, the Bernoulli term m^2/(2 rho)(1/A_to^2 - 1/A_from^2) does
(C++ ``segment_area_change_drop``). A constant segment is unchanged bit for
bit.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    CrossflowSegmentElement,
    FlowNetwork,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
X_AIR = list(cb.species.dry_air())


def _net(A1: float, A2: float | None, regime: str, L: float = 0.05, dp: float = 1e3) -> FlowNetwork:
    """Reservoir --s0 (A1)--> b1 --s1 (A1 -> A2)--> b2 --s2 (A2)--> exit, no
    side streams: two stations at which nothing joins."""
    A2_ = A1 if A2 is None else A2
    g = FlowNetwork()
    for n, pt in (("res", 1e5 + dp), ("out", 1e5)):
        b = PressureBoundary(n)
        b.Pt, b.Tt, b.Y = pt, 300.0, Y_AIR
        g.add_node(b)
    g.add_node(PlenumNode("b1"))
    g.add_node(PlenumNode("b2"))
    kw = {"length": L, "regime": regime}
    g.add_element(
        CrossflowSegmentElement(
            "s0", "res", "b1", area=A1, entry_K=0.0, to_kappa=0.0, next_seg="s1", **kw
        )
    )
    taper = {} if A2 is None else {"area_to": A2}
    g.add_element(
        CrossflowSegmentElement(
            "s1", "b1", "b2", area=A1, from_kappa=0.0, prev_seg="s0", to_kappa=0.0,
            next_seg="s2", **taper, **kw,
        )
    )  # fmt: skip
    g.add_element(
        CrossflowSegmentElement("s2", "b2", "out", area=A2_, from_kappa=0.0, prev_seg="s1", **kw)
    )
    return g


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    return s, r


@pytest.mark.parametrize(("A1", "A2"), [(2e-3, 1e-3), (1e-3, 2e-3)])
def test_a_frictionless_taper_is_bernoulli(A1: float, A2: float) -> None:
    """Incompressible, segments 0.1 mm long (friction negligible), 100 Pa
    drive: the chain passes Pt_res - P_out = q(A2), so m = A2 sqrt(2 rho dp),
    whatever the taper upstream. A constant-section s1 would leave b2 at
    b1's static and the flow at A1's."""
    _, r = _solve(_net(A1, A2, "incompressible", L=1e-4, dp=1e2))
    rho = cb.density(300.0, 1e5, X_AIR)
    assert float(r["s1.m_dot"]) == pytest.approx(A2 * math.sqrt(2.0 * rho * 1e2), rel=3e-3)


@pytest.mark.parametrize(("A1", "A2"), [(2e-3, 1e-3), (1e-3, 2e-3)])
def test_the_regimes_agree_at_low_mach(A1: float, A2: float) -> None:
    """Compressible (end closures at each end's section) and incompressible
    (the Bernoulli term) are built independently; at a throat Mach of
    0.04-0.07 they agree to O(M^2)."""
    _, rc = _solve(_net(A1, A2, "compressible", dp=1e2))
    _, ri = _solve(_net(A1, A2, "incompressible", dp=1e2))
    assert float(rc["s1.m_dot"]) == pytest.approx(float(ri["s1.m_dot"]), rel=2e-3)


def test_a_constant_segment_is_unchanged() -> None:
    """``area_to`` unset, or equal to ``area``: the same solution bit for bit."""
    for regime in ("incompressible", "compressible"):
        _, a = _solve(_net(2e-3, None, regime))
        _, b = _solve(_net(2e-3, 2e-3, regime))
        for k in ("s1.m_dot", "b1.P", "b2.P"):
            assert float(a[k]) == float(b[k])


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.parametrize("reverse", [False, True])
def test_the_taper_jacobian_matches_finite_differences(regime: str, reverse: bool) -> None:
    """Forward and reversed (the exit above the reservoir), a contraction.
    Reversed incompressible at 1 kPa: at 3 kPa the reversed jet into the 2:1
    narrow end runs at M ~0.55, where the incompressible regime's node
    densities differ by the dynamic head itself (17%) and the acceleration
    (outlet density) and the recovery (station density) no longer cancel:
    no root. The compressible regime solves it."""
    dp = (-1e3 if regime == "incompressible" else -3e3) if reverse else 3e3
    g = _net(2e-3, 1e-3, regime, dp=dp)
    s, r = _solve(g)
    assert (float(r["s1.m_dot"]) < 0.0) == reverse
    x0 = np.array(r["__x_solution__"], dtype=float)
    for x in (x0, x0 * (1.0 + 0.002 * np.sin(np.arange(len(x0))))):
        J = s._residuals_and_jacobian(x)[1].toarray()
        fd = fd_jacobian(s, x, rel=1e-7, floor=1e-12)
        cols = np.maximum(np.abs(x), 1e-12)[None, :]
        scale = np.maximum(np.max(np.abs(fd * cols), axis=1, keepdims=True), 1e-300)
        assert float(np.max(np.abs((J - fd) * cols) / scale)) < 1e-5


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_a_taper_passes_the_second_law(regime: str) -> None:
    """Total pressure carries across the taper bar friction: no segment
    creates it. A constant-section segment standing in for the taper (the
    #812 prototype's finding) read sigma_norm about -0.06."""
    for A1, A2 in ((2e-3, 1e-3), (1e-3, 2e-3)):
        _, r = _solve(_net(A1, A2, regime, dp=5e3))
        verdicts = {k: v["verdict"] for k, v in r["__second_law__"].items()}
        assert set(verdicts.values()) == {"pass"}, verdicts


def test_the_area_change_kernel() -> None:
    """C++ segment_area_change_drop: Bernoulli, even in m, zero at equal
    areas, analytic partials."""
    r = cb.segment_area_change_drop(0.1, 1e5, 300.0, X_AIR, 2e-3, 1e-3)
    rho = cb.density(300.0, 1e5, X_AIR)
    assert r.dP == pytest.approx(0.01 / (2 * rho) * (1 / 1e-6 - 1 / 4e-6), rel=1e-12)
    assert cb.segment_area_change_drop(-0.1, 1e5, 300.0, X_AIR, 2e-3, 1e-3).dP == r.dP
    assert cb.segment_area_change_drop(0.1, 1e5, 300.0, X_AIR, 2e-3, 2e-3).dP == 0.0
    h = 1e-6
    fd_m = (
        cb.segment_area_change_drop(0.1 + h, 1e5, 300.0, X_AIR, 2e-3, 1e-3).dP
        - cb.segment_area_change_drop(0.1 - h, 1e5, 300.0, X_AIR, 2e-3, 1e-3).dP
    ) / (2 * h)
    fd_P = (
        cb.segment_area_change_drop(0.1, 1e5 + 1.0, 300.0, X_AIR, 2e-3, 1e-3).dP
        - cb.segment_area_change_drop(0.1, 1e5 - 1.0, 300.0, X_AIR, 2e-3, 1e-3).dP
    ) / 2.0
    fd_T = (
        cb.segment_area_change_drop(0.1, 1e5, 300.0 + 1e-3, X_AIR, 2e-3, 1e-3).dP
        - cb.segment_area_change_drop(0.1, 1e5, 300.0 - 1e-3, X_AIR, 2e-3, 1e-3).dP
    ) / 2e-3
    assert r.d_dm_a == pytest.approx(fd_m, rel=1e-8)
    assert r.d_dP == pytest.approx(fd_P, rel=1e-6)
    assert r.d_dT == pytest.approx(fd_T, rel=1e-6)
