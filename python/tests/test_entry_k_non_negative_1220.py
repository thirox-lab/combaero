"""A crossflow segment refuses a negative or non-finite entry_K (#1220).

A negative K is an entry that adds total pressure. That breaks the second law,
so it is invalid input, not an extrapolation. Before #1220 only
ImpingementCrossflowElement checked the sign. The base CrossflowSegmentElement,
and the EffusionLiner builder that feeds it, passed any value through: K = -0.5
raised the flow by 35 %, and K = -1 divided by zero in K / (1 + K).
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network import FlowNetwork
from combaero.network import components as cbn
from combaero.network.effusion_liner import EffusionLiner

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _segment(entry_K: float | None) -> cbn.CrossflowSegmentElement:
    return cbn.CrossflowSegmentElement("s", "a", "b", length=0.1, area=1.0e-4, entry_K=entry_K)


@pytest.mark.parametrize("K", [-0.5, -1.0, -1.0e-12, math.nan, math.inf, -math.inf])
def test_segment_refuses(K: float) -> None:
    with pytest.raises(ValueError, match=r"'s': entry_K must be finite and non-negative"):
        _segment(K)


@pytest.mark.parametrize("K", [None, 0.0, 0.5, 1.0e3])
def test_segment_accepts(K: float | None) -> None:
    assert _segment(K).entry_K == K


def test_effusion_liner_builder_refuses() -> None:
    g = FlowNetwork()
    g.add_node(cbn.PressureBoundary("cin", Pt=1.5e5, Tt=600.0, Y=Y_AIR))
    g.add_node(cbn.PressureBoundary("cout", Pt=1.1e5, Tt=600.0, Y=Y_AIR))
    g.add_node(cbn.PressureBoundary("gas", Pt=1.0e5, Tt=1500.0, Y=Y_AIR))
    liner = EffusionLiner(
        "ln", n_segments=3, length=0.152, width=0.152, duct_height=0.01, hole_diameter=3.27e-3,
        wall_thickness=6.3e-3, pitch_x=15.24e-3, pitch_y=15.24e-3, angle_deg=30.0, entry_K=-0.5,
    )  # fmt: skip
    with pytest.raises(ValueError, match="entry_K must be finite and non-negative"):
        liner.add_to(g, "cin", "cout", "gas")
