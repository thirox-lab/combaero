"""A merge chamber whose main inlet reverses keeps its face (#493).

The chamber hands its declared main inlet the face state of its impulse
balance. That balance holds for either sign of the main flow (the x-momentum
flux through the face is m_main^2 / (rho A) both ways), but the face used to
switch off once the main was no longer an inflow: at m_main -> 0+ the face sits
below the chamber pressure (the side streams carry momentum out), so the main's
exit pressure stepped 245 Pa and its residual 3.5 kPa as its flow crossed zero.
"""

from __future__ import annotations

import numpy as np
import pytest
import test_effusion_wall as tew
from _fd_jacobian import fd_jacobian


@pytest.fixture(scope="module")
def liner() -> tuple:
    s, r = tew._liner()
    return s, np.array(r["__x_solution__"])


def _duct_row_and_face(s, x: np.ndarray, m: float) -> tuple[float, tuple | None]:
    j = s.unknown_names.index("duct.m_dot")
    xp = x.copy()
    xp[j] = m
    s._residuals_and_jacobian(xp, compute_jacobian=False)
    el = s.network.elements["duct"]
    chamber = s.network.nodes["ch"]
    st_in = s._get_node_state(s.network.nodes["gas_in"], xp)
    st_out = s._get_node_state(chamber, xp)
    st_in.m_dot = st_out.m_dot = m
    face = s._main_face(el, chamber, st_out)
    res, _ = el.residuals(st_in, st_out)
    return res[0], face


def test_the_face_is_continuous_through_a_reversing_main(liner: tuple) -> None:
    s, x = liner
    m_star = x[s.unknown_names.index("duct.m_dot")]
    r_pos, f_pos = _duct_row_and_face(s, x, 1e-9 * m_star)
    r_neg, f_neg = _duct_row_and_face(s, x, -1e-9 * m_star)
    assert f_neg is not None  # it switched off here
    assert f_neg[0] == pytest.approx(f_pos[0], abs=1e-3)
    assert r_neg == pytest.approx(r_pos, abs=1e-3)
    # Further reversed: still a face, moving smoothly away.
    _, f_far = _duct_row_and_face(s, x, -0.1 * m_star)
    assert f_far is not None and f_far[0] < f_pos[0]


@pytest.mark.parametrize("frac", [1.0, -0.2, -1e-3])
def test_the_jacobian_holds_with_the_main_reversed(liner: tuple, frac: float) -> None:
    s, x = liner
    j = s.unknown_names.index("duct.m_dot")
    xp = x.copy()
    xp[j] = frac * x[j]
    _, J = s._residuals_and_jacobian(xp)
    # Flows step 1e-5, pressures and temperatures 1e-7: a compressible
    # channel's residual repeats only to ~1e-12 (its Fanno iteration is
    # warm-started), which a 1e-6 flow step differences into 3e-6 of noise.
    rel = np.array([1e-5 if n.endswith(".m_dot") else 1e-7 for n in s.unknown_names])
    fd = fd_jacobian(s, xp, rel=rel)
    # RESIDUAL-ROW-SCALING: compare J_ij |x_j| (see _build_residual_scales).
    cols = np.maximum(np.abs(xp), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    assert float(np.max(np.abs(Js - fds) / scale)) < 1e-6
