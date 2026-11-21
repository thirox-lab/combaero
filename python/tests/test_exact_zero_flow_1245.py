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
  scaled on the walled combustor's exit orifice row).
"""

from __future__ import annotations

import math

import pytest
from _jacobian_sweep import Case, Report, evaluate, point

import combaero as cb
from combaero.network import ConvectiveSurface


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
