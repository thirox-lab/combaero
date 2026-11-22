"""Ideal momentum junction kernel (#517) against the variant already scored.

The merge side of the ideal junction is Blaisdell & Manson's theory, which the
alpha record scores as ``BM theory`` (``joining_variants.bassett_cv(...,
wall=False)``). Those numbers are the issue's pre-registered baseline, so the
kernel must reproduce that variant, not just a formula re-derived beside it.
The C++ tests (tests/test_ideal_junction.cpp) hold the closed forms, energy,
and the Jacobian; this file holds the binding to the scored code.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from combaero import _core
from validation.junction import joining_variants as jv

RHO = 1.2
# Role-smoothing bias, relative to the collector's dynamic head (see the C++
# tests for the measured size).
SMOOTH_TOL = 5e-5


def _kernel_K(U: np.ndarray, A: np.ndarray, th: np.ndarray) -> dict[int, float]:
    """K_j = (Pt_j - Pt_c) / (rho u_c^2 / 2) from the kernel, incompressible.

    joining_variants: axes point away from the junction (the kernel's
    axis_rad as is) and U is positive INTO it, so m = -rho U A.
    """
    m = [-RHO * float(u) * float(a) for u, a in zip(U, A, strict=True)]
    geom = _core.IdealJunctionGeometry()
    geom.area = [float(a) for a in A]
    geom.axis_rad = [float(t) for t in th]
    geom.port_sign = [1.0] * len(A)
    n = len(A)
    r = _core.ideal_junction_residuals_and_jacobian([0.0] * n, [RHO] * n, [0.0] * n, m, 0.0, geom)
    P = [-x for x in r.residual[:n]]
    Pt = [p + 0.5 * RHO * float(u) ** 2 for p, u in zip(P, U, strict=True)]
    c = int(np.flatnonzero(U < 0)[0])
    head = 0.5 * RHO * float(U[c]) ** 2
    return {j: (Pt[j] - Pt[c]) / head for j in range(n) if j != c}


@pytest.mark.parametrize("psi", [1.0, 2.12, 4.0, 7.11, 16.0])
@pytest.mark.parametrize("deflection_deg", [15, 30, 45, 60, 75, 90, 105, 120, 135, 150])
def test_merge_is_the_scored_bm_theory_variant(psi: float, deflection_deg: float) -> None:
    for q in np.linspace(0.05, 0.95, 10):
        U, A, th = jv.tee(float(q), psi, deflection_deg)
        ref = jv.MODELS["BM theory"](U, A, th)
        got = _kernel_K(U, A, th)
        for j, K in ref.items():
            assert got[j] == pytest.approx(K, abs=SMOOTH_TOL), (q, j)


def test_split_branch_is_idelchik_7_6_with_A_prime_one() -> None:
    """Diverging branch: zeta = 1 + (w_b/w_c)^2 - 2 (w_b/w_c) cos a, referenced
    to the arriving head (Idelchik Section VII para 9, Eq 7-6, K'_b = 0 for a
    constant-section main, A' = 1)."""
    for psi in (1.0, 3.0):
        for a_deg in (30.0, 45.0, 60.0, 90.0):
            for qb in (0.1, 0.4, 0.7):
                a = math.radians(a_deg)
                U = np.array([1.0, -(1.0 - qb), -qb * psi])  # + into the junction
                A = np.array([1.0, 1.0, 1.0 / psi])
                th = np.array([math.pi, 0.0, a])
                m = [-RHO * float(u) * float(x) for u, x in zip(U, A, strict=True)]
                geom = _core.IdealJunctionGeometry()
                geom.area, geom.axis_rad, geom.port_sign = list(A), list(th), [1.0] * 3
                r = _core.ideal_junction_residuals_and_jacobian(
                    [0.0] * 3, [RHO] * 3, [0.0] * 3, m, 0.0, geom
                )
                P = [-x for x in r.residual[:3]]
                Pt = [p + 0.5 * RHO * float(u) ** 2 for p, u in zip(P, U, strict=True)]
                w = qb * psi  # w_b / w_c
                zeta = 1.0 + w * w - 2.0 * w * math.cos(a)
                assert (Pt[0] - Pt[2]) / (0.5 * RHO) == pytest.approx(zeta, abs=SMOOTH_TOL)


def test_binding_reports_velocities_and_column_layout() -> None:
    geom = _core.IdealJunctionGeometry()
    geom.area = [1e-3, 1e-3, 5e-4]
    geom.axis_rad = [math.pi, 0.0, math.pi / 2]
    geom.port_sign = [-1.0, 1.0, 1.0]  # inlet element flows INTO the junction
    r = _core.ideal_junction_residuals_and_jacobian(
        [1e5] * 3, [1.2] * 3, [1.2e-5] * 3, [0.01, 0.006, 0.004], 1e5, geom
    )
    assert len(r.residual) == 4
    assert len(r.jacobian) == 4 and all(len(row) == 7 for row in r.jacobian)
    assert r.velocity[0] < 0.0  # junction convention: inflow is negative
    assert r.residual[3] == pytest.approx(0.0, abs=1e-15)  # -0.01 + 0.006 + 0.004
    assert r.jacobian[3][3] == pytest.approx(-1.0)  # mass row vs outer_mdot[0]
    assert r.jacobian[0][6] == pytest.approx(-1.0)  # port row vs p_jn


def _split_Pt_drop(f: list[float] | None) -> tuple[float, float]:
    """(Pt_supply - Pt_straight, Pt_supply - Pt_branch) on a 90-deg split, psi 2."""
    U = np.array([1.0, -0.6, -0.4 * 2.0])  # + into the junction
    A = np.array([1.0, 1.0, 0.5])
    m = [-RHO * float(u) * float(a) for u, a in zip(U, A, strict=True)]
    geom = _core.IdealJunctionGeometry()
    geom.area, geom.axis_rad, geom.port_sign = list(A), [math.pi, 0.0, math.pi / 2], [1.0] * 3
    if f is not None:
        geom.loss_factor = f
    r = _core.ideal_junction_residuals_and_jacobian([0.0] * 3, [RHO] * 3, [0.0] * 3, m, 0.0, geom)
    Pt = [-x + 0.5 * RHO * float(u) ** 2 for x, u in zip(r.residual[:3], U, strict=True)]
    return Pt[0] - Pt[1], Pt[0] - Pt[2]


def test_loss_factor_scales_each_receivers_loss_through_the_binding() -> None:
    """#526: empty and all-ones are the closure; a factor scales its own port."""
    base = _split_Pt_drop(None)
    assert _split_Pt_drop([1.0, 1.0, 1.0]) == base
    scaled = _split_Pt_drop([0.0, 0.5, 1.5])
    # Each factor fades to 1 with its port's flow share (IDEAL_LOSS_FADE_REL
    # = 1e-2 of |m|_2), so the effective one is 1 + (f - 1) phi.
    m = np.array([-RHO * 1.0, RHO * 0.6, RHO * 0.8 * 0.5])  # rho U A, + out
    d2 = 1e-4 * float(np.sum(m * m))
    phi = m * m / (m * m + d2)
    assert scaled[0] == pytest.approx((1.0 - 0.5 * phi[1]) * base[0], rel=1e-6)
    assert scaled[1] == pytest.approx((1.0 + 0.5 * phi[2]) * base[1], rel=1e-6)


@pytest.mark.parametrize("f", [[1.0, -0.1, 1.0], [1.0, 1.0]])
def test_loss_factor_refuses_negative_or_misshaped(f: list[float]) -> None:
    with pytest.raises(ValueError, match="loss_factor"):
        _split_Pt_drop(f)


# #1250: the sweep point "jct pattern 0001 ratio 3" (Case ideal_nport 2/1/0),
# the kernel's own inputs there. Every port flows OUT (the inlet arm reversed),
# so no port supplies and V is a ratio of receiver weights w_in ~ -eps^2/(4 m),
# 1e-8 of the flow. Formed as m s - m they kept 4e-9 .. 2.4e-8 of relative
# error, and the rows carried 1.2e-5 .. 6.8e-5 Pa of noise on 1.8e5 Pa: a
# central difference at 1e-5 |m| read 7.7e-6 .. 8.2e-5 (scaled) off J.
_N1250 = {
    "p": [1.9806129994e05, 1.8888301935e05, 1.8817445859e05, 1.8499773995e05],
    "rho": [2.3000283783, 2.1934436703, 2.1852153599, 2.1483250485],
    "outer": [-0.23296145806, 0.23296145806, 0.23296145806, 0.69888437417],
    "sign": [-1.0, 1.0, 1.0, 1.0],
    "angle_deg": [90.0, 45.0, 135.0, 135.0],
    "area": 0.001963495408493621,
    "p_jn": 1.8064321995e05,
}


def _kernel_1250(outer: list[float]) -> object:
    c = _N1250
    geom = _core.IdealJunctionGeometry()
    geom.area = [c["area"]] * 4
    geom.axis_rad = [
        math.radians(a) if s > 0.0 else math.pi - math.radians(a)
        for a, s in zip(c["angle_deg"], c["sign"], strict=True)
    ]
    geom.port_sign = list(c["sign"])
    return _core.ideal_junction_residuals_and_jacobian(
        c["p"], c["rho"], [0.0] * 4, outer, c["p_jn"], geom
    )


@pytest.mark.parametrize("col", range(4))
def test_no_supplier_rows_are_smooth_and_match_fd_1250(col: int) -> None:
    """#1250: with no supplier the rows are smooth to round-off in every flow,
    so a small-step difference meets J."""
    outer = list(_N1250["outer"])
    m = abs(outer[col])
    J = np.array(_kernel_1250(outer).jacobian)[:4, 4 + col]

    def F(d: float) -> np.ndarray:
        o = list(outer)
        o[col] += d
        return np.array(_kernel_1250(o).residual[:4])

    scale = max(abs(p) for p in _N1250["p"])
    # Noise: departure from a quadratic over +-2e-5 |m| (was 8e-10 of the
    # row scale; round-off is ~1e-16).
    ds = np.linspace(-2e-5, 2e-5, 41) * m
    Fs = np.array([F(d) for d in ds])
    for i in range(4):
        fit = np.polyval(np.polyfit(ds, Fs[:, i], 2), ds)
        assert np.max(np.abs(Fs[:, i] - fit)) < 1e-14 * scale, i
    h = 1e-5 * m
    d1 = (F(h) - F(-h)) / (2.0 * h)
    d2 = (F(2.0 * h) - F(-2.0 * h)) / (4.0 * h)
    fd = (4.0 * d1 - d2) / 3.0
    assert np.max(np.abs(J - fd)) * m / scale < 1e-8
