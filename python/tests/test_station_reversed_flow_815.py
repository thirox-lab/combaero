"""A crossflow station holds for flow in either direction (#815).

Over a station, momentum gives P_a - P_b = (m_b^2 - m_a^2 + kappa (m_a - m_b)
s) / (rho A^2), s the arriving streams' mean (m_a forward; #1178 for the
other signs): the momentum flux m u = m^2/(rho A) is positive
whichever way the gas moves. ``station_half_drop`` used m|m|, which is the same
for forward flow and flips the drop's sign for reversed flow: a merge in a
reversed chain RAISED the static pressure along the flow (a pump), and a
station the flow leaves both ways (a side stream splitting) had no physical
root -- incompressible did not converge, compressible landed on a root the
second law rejected.
"""

from __future__ import annotations

import pytest

import combaero as cb
from combaero.network import (
    CrossflowSegmentElement,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _boundaries(g: FlowNetwork, pts: dict[str, float]) -> None:
    for n, pt in pts.items():
        b = PressureBoundary(n)
        b.Pt, b.Tt, b.Y = pt, 300.0, Y_AIR
        g.add_node(b)


def _chain(mirror: bool, regime: str) -> FlowNetwork:
    """res -> b1 -> b2 -> out with a normal jet joining at b1 and at b2. The
    mirror builds the same chain from ``out`` to ``res``: every segment flow
    is negative, the physics identical."""
    g = FlowNetwork()
    _boundaries(g, {"res": 1.02e5, "out": 1.0e5, "sup": 1.05e5})
    g.add_node(PlenumNode("b1"))
    g.add_node(PlenumNode("b2"))
    kw = {"length": 0.05, "area": 2e-3, "regime": regime}
    if not mirror:
        g.add_element(
            CrossflowSegmentElement(
                "s0", "res", "b1", entry_K=0.0, to_kappa=0.0, next_seg="s1", **kw
            )
        )
        g.add_element(
            CrossflowSegmentElement(
                "s1", "b1", "b2", from_kappa=0.0, prev_seg="s0", to_kappa=0.0,
                next_seg="s2", **kw,
            )
        )  # fmt: skip
        g.add_element(
            CrossflowSegmentElement("s2", "b2", "out", from_kappa=0.0, prev_seg="s1", **kw)
        )
    else:
        # The reservoir is the to-end the reversed flow comes from: a node at
        # rest feeding the segment, its own Pt (#809), the entry_K = 0 above.
        g.add_element(CrossflowSegmentElement("s2", "out", "b2", to_kappa=0.0, next_seg="s1", **kw))
        g.add_element(
            CrossflowSegmentElement(
                "s1", "b2", "b1", from_kappa=0.0, prev_seg="s2", to_kappa=0.0,
                next_seg="s0", **kw,
            )
        )  # fmt: skip
        g.add_element(
            CrossflowSegmentElement("s0", "b1", "res", from_kappa=0.0, prev_seg="s1", **kw)
        )
    for i, b in ((1, "b1"), (2, "b2")):
        g.add_element(
            OrificeElement(
                f"h{i}", "sup", b, Cd=0.7, diameter=0.012, correlation="fixed", regime=regime
            )
        )
    return g


def _solve(g: FlowNetwork) -> dict:
    r = NetworkSolver(g).solve()
    assert r["__success__"], r.get("__message__")
    return r


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_a_reversed_chain_is_its_mirror(regime: str) -> None:
    """Compressible: identical to round-off. Incompressible: to 0.5%, the
    reservoir end's density reference (the entry's at the station, the
    reversed inflow's at the reservoir). Under m|m| the mirror ran 1.7x the
    flow with the stations below the exit pressure."""
    fwd = _solve(_chain(False, regime))
    mir = _solve(_chain(True, regime))
    tol = 1e-9 if regime == "compressible" else 6e-3
    for i in range(3):
        assert -float(mir[f"s{i}.m_dot"]) == pytest.approx(float(fwd[f"s{i}.m_dot"]), rel=tol)
    for n in ("b1", "b2"):
        assert float(mir[f"{n}.P"]) == pytest.approx(float(fwd[f"{n}.P"]), rel=tol * 1e-2)


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.parametrize("p_sup", [1.1e5, 1.3e5])
def test_a_split_station_has_a_physical_root(regime: str, p_sup: float) -> None:
    """#815: a normal jet into the middle of a chain whose ends sit at one
    pressure splits both ways, symmetrically, and the root passes the
    second law."""
    g = FlowNetwork()
    _boundaries(g, {"res": 1.0e5, "sup": p_sup, "out": 1.0e5})
    g.add_node(PlenumNode("st"))
    kw = {"length": 0.02, "area": 2e-3, "regime": regime}
    g.add_element(
        CrossflowSegmentElement("s0", "res", "st", entry_K=0.0, to_kappa=0.0, next_seg="s1", **kw)
    )
    g.add_element(CrossflowSegmentElement("s1", "st", "out", from_kappa=0.0, prev_seg="s0", **kw))
    g.add_element(
        OrificeElement("h", "sup", "st", Cd=0.8, diameter=0.02, regime=regime, correlation="fixed")
    )
    r = _solve(g)
    m0, m1 = float(r["s0.m_dot"]), float(r["s1.m_dot"])
    assert m0 < 0.0 < m1
    assert -m0 == pytest.approx(m1, rel=0.02)
    assert {v["verdict"] for v in r["__second_law__"].values()} == {"pass"}
