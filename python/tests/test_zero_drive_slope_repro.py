"""Zero-drive slopes, fixed for #1193, #1256 and #1209.

A row whose drive was exactly zero got no slope in its drive's pressures:
an orifice at dp = 0 (#1256, the trigger of #1193's chaotic Brayton start)
and a compressible duct at zero Fanno drive (#1209). Below kJetDriveFloor
(1e-3 Pa, flow_drive_floor.h) the jet kernels are now the floor's flow times a
C1 odd cubic with a finite slope at zero, and the Fanno kernels return their
sub-floor cubic's slope at exactly zero drive.
"""

from __future__ import annotations

import types
import warnings

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    CrossflowSegmentElement,
    FlowNetwork,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)

_Y_AIR = list(cb.species.dry_air_mass())


def _fresh(net: FlowNetwork) -> NetworkSolver:
    s = NetworkSolver(net)
    net.resolve_all_topology()
    s._build_x0()
    return s


# ---------------------------------------------------------------------------
# #1193 / #1256: the Brayton cold solve was chaotic in its start. Its start puts
# c_in and t_out at exactly the boundaries' Pt, where both orifice rows'
# pressure partials were 0 (|dp| < ~1.5e-10 Pa) or a round-off staircase
# (to ~2e-9 Pa), so J at x0 jumped between 0 and -22 on a 1e-10 Pa change.
# An FD Jacobian at the J requests within 1e-3 Pa of that point alone
# converged 40/40 seeds (analytic 35/40) and 14/14 offsets (analytic 10/14).
# ---------------------------------------------------------------------------


def _brayton_solver(offset: float, unknown: str = "t_out.Pt") -> NetworkSolver:
    from test_turbomachine_elements_893 import _brayton

    s = NetworkSolver(_brayton())
    orig = s._build_x0

    def bx(self: NetworkSolver) -> np.ndarray:
        x = np.array(orig(), dtype=float)
        x[list(self.unknown_names).index(unknown)] += offset
        return x

    s._build_x0 = types.MethodType(bx, s)
    return s


@pytest.mark.parametrize(
    "offset",
    [-3e-10, pytest.param(3e-10, marks=pytest.mark.slow)],  # 3.5 s: its false root
)
def test_1193_brayton_converges_from_a_sub_nanopascal_start_offset(offset: float) -> None:
    """t_out.Pt moved 3e-10 Pa off the exit's Pt: the cold solve ended on a
    false root (|F| 51.8) before #1256 on macOS arm64; an FD Jacobian at the
    zero-drop requests alone converged it (40/40 seeds at 1e-14 relative,
    against 35/40). The sub-floor cubic gives the rows that slope."""
    s = _brayton_solver(offset)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=60)
    assert r["__success__"], r["__message__"]


def test_1193_orifice_rows_keep_a_pressure_slope_at_the_zero_drop_start() -> None:
    """At x0 every orifice row's partials on its node's P and Pt were exactly
    0 (true one-sided slope: unbounded, ~K / sqrt(dp)). The sub-floor cubic
    gives them its finite slope; the row is not blind to its drive there."""
    from test_turbomachine_elements_893 import _brayton

    s = _fresh(_brayton())
    x = np.array(s._build_x0(), dtype=float)
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()
    names = list(s.unknown_names)
    kinds = s._row_kinds(x)
    for el, node in (("intake", "c_in"), ("nozzle", "t_out")):
        i = next(k for k, (o, _) in enumerate(kinds) if o == f"OrificeElement '{el}'")
        cols = [names.index(f"{node}.P"), names.index(f"{node}.Pt")]
        assert np.any(J[i, cols] != 0.0), f"{el}: J row blind to {node}.P/Pt at dp = 0"


# ---------------------------------------------------------------------------
# #1209 / #1216: a compressible duct at exactly zero drive had dG/dPt = 0;
# with the low-Mach viscous closure (#1216), both forward (into the plenum,
# exit head kept) and reverse (into the static-coupled boundary) recover the
# exact Hagen-Poiseuille slope 1.54119e-5 kg/(s Pa) at zero drive (before
# #1216, the static-coupled 1 Pa cubic was 0.82% low, giving a mean of
# 1.53487e-5).
# ---------------------------------------------------------------------------


def _zero_drive_net(kind: str) -> FlowNetwork:
    def pb(n: str, Pt: float) -> PressureBoundary:
        b = PressureBoundary(n)
        b.Pt, b.Tt, b.Y = Pt, 300.0, _Y_AIR
        return b

    g = FlowNetwork()
    g.add_node(pb("a", 1.01e5))
    g.add_node(pb("b", 1.0e5))
    g.add_node(PlenumNode("p"))
    if kind == "channel":
        g.add_element(
            ChannelElement("d1", "a", "p", length=1.0, diameter=0.01, regime="compressible")
        )
    else:
        g.add_element(
            CrossflowSegmentElement(
                "d1", "a", "p", length=1.0, area=np.pi * 0.005**2, regime="compressible"
            )
        )
    g.add_element(ChannelElement("d2", "p", "b", length=1.0, diameter=0.01, regime="compressible"))
    return g


@pytest.mark.parametrize("kind", ["channel", "crossflow_segment"])
def test_1209_duct_flow_row_has_its_sub_floor_slope_at_zero_drive(kind: str) -> None:
    s = _fresh(_zero_drive_net(kind))
    names = list(s.unknown_names)
    x = np.array(s._build_x0(), dtype=float)
    for j, n in enumerate(names):
        if n.endswith("m_dot"):
            x[j] = 0.0
    x[names.index("p.P")] = x[names.index("p.Pt")] = 1.01e5  # drive of d1 exactly 0
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()
    i = next(k for k, (o, _) in enumerate(s._row_kinds(x)) if o.endswith("'d1'"))
    j = names.index("p.Pt")
    h = 1e-3
    xp, xm = x.copy(), x.copy()
    xp[j] += h
    xm[j] -= h
    fd = (s._residuals(xp)[i] - s._residuals(xm)[i]) / (2.0 * h)
    assert fd == pytest.approx(1.5412e-5, rel=1e-3)  # Hagen-Poiseuille slope (#1216)
    assert J[i, j] == pytest.approx(fd, rel=1e-4)
