"""J at a boundary-fed flow of exactly zero (#1245).

A stream at exactly zero flow beside an arriving one is in its node's mix
with weight 0 (streams_in_mix, #1232): its partials are those of the side it
arrives from. J there must equal that one-sided difference. It equalled
neither side, in two places:

* the momentum chamber closure zeroed d(res)/dY_k (and drho_dY, dT_static_dY)
  for a species absent at the node, though the zero-flow stream brings it and
  the relay moves it: J took the arriving side's T and present species and
  none of the absent ones (1.3e-3 scaled on a tee's port chamber, 2.1e-2 on
  the ejector's outlet chamber);
* a ConvectiveSurface returned no result below 1e-12 m/s, so a wall with a
  duct at rest on one side dropped its whole coupling: Q is 0 either way, but
  its flow partial read 0 against a forward slope of about -0.21 (4.0e-2
  scaled on the walled combustor's exit orifice row);
* with that gate gone, both sides of a wall at exactly zero reached the
  exchange's both-stopped branch, which returned dK/dC = 0.5 per side (the
  slope along C_a = C_b) where K(C, 0) = 0 for every C: both partials are 0.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from _jacobian_sweep import Case, Report, _dense, _start, evaluate, point

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FlowNetwork,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)


@pytest.mark.parametrize(
    ("case", "desc"),
    [
        # Tee port chambers: the absent species' closure partials.
        (Case("legacy_tee_ports", 2, 2, 1), "lc.m_dot = 0"),
        # A coolant duct at rest under a wall.
        (Case("walled_combustor", 2, 2, 1), "cold.m_dot = 0"),
    ],
)
def test_exact_zero_boundary_flow_matches_the_arriving_side(case: Case, desc: str) -> None:
    x0, pt = point(case, desc)
    assert pt.kind == "exact_zero"
    rep = Report()
    evaluate(case, pt, x0, rep)
    assert rep.compared == 1
    assert not rep.mismatches(), "\n".join(f.message for f in rep.mismatches())


def test_surface_at_rest_has_its_laminar_coefficient() -> None:
    """A duct surface at zero velocity gives the laminar limit, the value its
    correlation approaches from either side, not no result."""
    surf = ConvectiveSurface(area=0.05)
    X = list(cb.mass_to_mole(list(cb.species.dry_air_mass())))
    r0 = surf.htc_and_T(T=300.0, P=2e5, X=X, velocity=0.0, diameter=0.03, length=0.5)
    r1 = surf.htc_and_T(T=300.0, P=2e5, X=X, velocity=1e-9, diameter=0.03, length=0.5)
    assert r0 is not None and math.isfinite(r0.h) and r0.h > 0.0
    assert r0.h == pytest.approx(r1.h, rel=1e-9)


def _both_ducts_of_a_wall_at_rest() -> tuple[NetworkSolver, np.ndarray, list[int]]:
    """A hot and a cold duct coupled by a wall, each into a chamber an
    unwalled bypass keeps flowing; both walled ducts at exactly zero."""
    Y = list(cb.species.dry_air_mass())
    D, L = 0.02, 0.3
    g = FlowNetwork()
    for side, T in (("h", 900.0), ("c", 300.0)):
        g.add_node(PressureBoundary(f"{side}i", Pt=3e5, Tt=T, Y=Y))
        g.add_node(PressureBoundary(f"{side}o", Pt=2.5e5, Tt=300.0, Y=Y))
        g.add_node(MomentumChamberNode(f"{side}m", area=math.pi * D * D / 4))
        g.add_element(
            ChannelElement(
                f"{side}ch", f"{side}i", f"{side}m", length=L, diameter=D,
                regime="incompressible", surface=ConvectiveSurface(area=math.pi * D * L),
            )
        )  # fmt: skip
        g.add_element(
            ChannelElement(
                f"{side}b", f"{side}i", f"{side}m", length=L, diameter=D, regime="incompressible"
            )
        )
        g.add_element(
            ChannelElement(
                f"{side}x", f"{side}m", f"{side}o", length=0.05, diameter=D,
                regime="incompressible",
            )
        )  # fmt: skip
    g.add_wall(
        ThermalWall(id="w", element_a="hch", element_b="cch", layers=[WallLayer(0.002, 25.0)])
    )
    s = NetworkSolver(g)
    g.resolve_all_topology()
    s._build_x0()
    x = _start(s).copy()
    cols = [s.unknown_names.index("hch.m_dot"), s.unknown_names.index("cch.m_dot")]
    x[cols] = 0.0
    return s, x, cols


def test_both_ducts_of_a_wall_at_rest_match_the_forward_difference() -> None:
    """Each stopped duct's column against the forward difference, Richardson
    extrapolated, row-scaled; before the fix rows read 8.4e5 where both
    one-sided differences are 0."""
    s, x, cols = _both_ducts_of_a_wall_at_rest()
    F, J = s._residuals_and_jacobian(x)
    F, J = np.asarray(F), _dense(J)
    scale = np.maximum(np.abs(F), np.max(np.abs(J) * np.abs(x)[None, :], axis=1))
    for j in cols:
        h = 1e-8
        d = []
        for step in (h, 2.0 * h):
            xp = x.copy()
            xp[j] = step
            d.append((np.asarray(s._residuals_and_jacobian(xp)[0]) - F) / step)
        fwd = 2.0 * d[0] - d[1]
        err = np.abs(J[:, j] - fwd) * 1e-3 / scale  # a column scale of 1 g/s
        assert err.max() < 1e-6, (s.unknown_names[j], int(err.argmax()), err.max())
