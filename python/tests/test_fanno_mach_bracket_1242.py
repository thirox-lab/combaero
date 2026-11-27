"""Regression tests for #1242: adiabatic Fanno Mach bracket on near-zero starts.

On main (e427b904), adiabatic `raw_flow` (src/fanno_mach.cpp) stopped halving
the lower exit Mach bound at 1e-6 (lo = 2^-20 = 9.537e-7 <= 1e-6), throwing
`RuntimeError: fanno (Mach) brent: no sign change on [9.536743e-07, 1.907349e-06]`
whenever the unchoked exit Mach at the 1 Pa drive floor fell below 1e-6 -- both
during intermediate high-Tt probes of `net_deadend` (`sc=1e-11`, 11 caught
exceptions on main) and on normal-temperature narrow laminar ducts.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    EnergyBoundary,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

AIR_X = list(cb.species.dry_air())
AIR_Y = list(cb.mole_to_mass(cb.species.dry_air()))
CH4_Y = [0.0] * len(AIR_Y)
CH4_Y[cb.species.species_index("CH4")] = 1.0


def _net_deadend() -> FlowNetwork:
    """Heated dead-end branch network from #1242 (review1238/scen.py)."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=1.2e5, Tt=500.0, Y=AIR_Y))
    g.add_node(PressureBoundary("b", Pt=1.0e5, Tt=300.0, Y=AIR_Y))
    g.add_node(PlenumNode("M"))
    g.add_node(PlenumNode("D"))
    g.add_element(OrificeElement("e1", "a", "M", Cd=0.7, diameter=0.01, correlation="fixed"))
    g.add_element(OrificeElement("e2", "M", "b", Cd=0.7, diameter=0.01, correlation="fixed"))
    g.add_element(ChannelElement("e3", "M", "D", length=0.5, diameter=0.02))
    g.add_element(OrificeElement("e5", "D", "M", Cd=0.7, diameter=0.006, correlation="fixed"))
    g.nodes["D"].add_energy_boundary(EnergyBoundary("q", Q=1.0e3))
    g.add_node(PressureBoundary("f", Pt=1.15e5, Tt=300.0, Y=CH4_Y))
    g.add_element(OrificeElement("e4", "f", "M", Cd=0.7, diameter=0.003, correlation="fixed"))
    return g


@pytest.mark.parametrize(
    ("Pt0", "Tt0", "drive", "exit_total", "L", "D", "rough", "model", "fmult"),
    [
        # Exact failing intermediate calls captured on main during net_deadend (sc=1e-11).
        (1.110425e5, 5.953859e5, 1.0, True, 0.5, 0.02, 1e-5, "haaland", 1.0),
        (1.157294e5, 3.076902e6, 1.0, True, 0.5, 0.02, 1e-5, "haaland", 1.0),
        (1.157294e5, 8.25e6, 1.0, True, 0.5, 0.02, 1e-5, "haaland", 1.0),
        # Normal-temperature narrow laminar duct (Hagen-Poiseuille Me ~ 1.17e-7 < 1e-6).
        (1.0e5, 1000.0, 1.0, True, 0.1, 1.0e-4, 0.0, "haaland", 1.0),
        (1.0e5, 1000.0, 1.0, False, 0.1, 1.0e-4, 0.0, "haaland", 1.0),
    ],
)
def test_fanno_channel_flow_brackets_below_1e_6_mach(
    Pt0: float,
    Tt0: float,
    drive: float,
    exit_total: bool,
    L: float,
    D: float,
    rough: float,
    model: str,
    fmult: float,
) -> None:
    r = cb.fanno_channel_flow(Pt0, Tt0, AIR_X, Pt0 - drive, exit_total, L, D, rough, model, fmult)
    assert not r.choked
    assert r.G > 0.0
    assert 0.0 < r.M_exit < 1.0e-6
    assert math.isfinite(r.dG_dP_target) and r.dG_dP_target < 0.0
    assert math.isfinite(r.dG_dPt0) and r.dG_dPt0 > 0.0
    L_check = cb.fanno_length_between(r.G, Tt0, r.M_in, r.M_exit, AIR_X, D, rough, model, fmult)
    assert L_check == pytest.approx(L, rel=1e-8)


@pytest.mark.parametrize("sc", [1e-14, 1e-11])
def test_deadend_near_zero_start_raises_no_fanno_brent_exception(sc: float) -> None:
    """#1242 reproducer: starting mass flows at 1e-11 or 1e-14 kg/s must not
    trigger any `fanno (Mach) brent: no sign change` exception inside the
    solver's residual evaluations (main caught 11 at sc=1e-11)."""
    g = _net_deadend()
    g.resolve_all_topology()
    s = NetworkSolver(g)
    x0 = np.array(s._build_x0(), dtype=float)
    for i, name in enumerate(s.unknown_names):
        if "m_dot" in name:
            x0[i] = sc
    sol = s.solve(x0=x0)
    assert sol.get("__caught_exception_count__", 0) == 0
    if sc == 1e-14:
        assert sol["__converged__"] is True
