"""A crossflow chain's exit into a CombustorNode is an exit at rest (#1136).

A CombustorNode declares ``state_at_rest`` (Pt = P), as a plenum and every
boundary do. The chain asked "a PlenumNode or a boundary" instead, so a
combustor exit was MOVING to the compressible terminal (the duct's exit total
met the combustor's Pt = P: the head recovered) and at rest to the
incompressible row (static continuity forward; reversed, without the
reservoir acceleration a plenum exit carries). The regimes differed by -21 %
forward and +179 % reversed at 2 kPa. Now every terminal reader asks the
node's own definition (``_terminal_at_rest``): a combustor exit is a plenum
exit in both regimes.

The probe: a reservoir entry, one bleed station, the exit segment into the
node, and a channel from the node to the exit.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian, scaled_error

import combaero as cb
from combaero.network import (
    ChannelElement,
    CombustorNode,
    CrossflowSegmentElement,
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkMixtureState,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    WallNode,
)

Y = list(cb.mole_to_mass(cb.species.dry_air()))
A_DUCT = 2.0e-4
BLEED = cb.STATION_KAPPA_BLEED_BASSETT
P0 = 1.0e5


def _pb(name: str, Pt: float) -> PressureBoundary:
    return PressureBoundary(name, Pt=Pt, Tt=300.0, Y=Y)


def _chain(kind: str, regime: str, p_sup: float, p_exit: float) -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(_pb("sup", p_sup))
    g.add_node(_pb("exit", p_exit))
    g.add_node(_pb("sink", P0))
    g.add_node(PlenumNode("s1"))
    g.add_node(CombustorNode("nd") if kind == "combustor" else PlenumNode("nd"))
    g.add_element(
        CrossflowSegmentElement(
            "s0x",
            "sup",
            "s1",
            length=0.05,
            area=A_DUCT,
            entry_K=0.5,
            to_kappa=BLEED,
            next_seg="s1x",
            regime=regime,
        )
    )
    g.add_element(
        CrossflowSegmentElement(
            "s1x",
            "s1",
            "nd",
            length=0.05,
            area=A_DUCT,
            from_kappa=BLEED,
            prev_seg="s0x",
            regime=regime,
        )
    )
    g.add_element(OrificeElement("bl", "s1", "sink", diameter=0.006, Cd=0.6, correlation="fixed"))
    D = math.sqrt(4.0 * A_DUCT / math.pi)
    g.add_element(ChannelElement("tail", "nd", "exit", length=0.05, diameter=D, regime=regime))
    return g


def _ends(dp: float, reversed_: bool) -> tuple[float, float]:
    return (P0, P0 + dp) if reversed_ else (P0 + dp, P0)


def _solve(g: FlowNetwork) -> dict:
    r = NetworkSolver(g).solve()
    assert r["__success__"], r.get("__message__")
    return r


CASES = [(dp, rev) for dp in (2000.0, 8000.0) for rev in (False, True)]


@pytest.mark.parametrize(("dp", "reversed_"), CASES)
@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_combustor_exit_is_a_plenum_exit(dp: float, reversed_: bool, regime: str) -> None:
    """Both regimes, both directions: the combustor exit solves the plenum
    exit's flow. On main compressible forward ran 26.7 % above it (the head
    recovered) and incompressible reversed 2.8x (no reservoir acceleration)."""
    ps, pe = _ends(dp, reversed_)
    mc = float(_solve(_chain("combustor", regime, ps, pe))["s1x.m_dot"])
    mp = float(_solve(_chain("plenum", regime, ps, pe))["s1x.m_dot"])
    assert (mc < 0.0) == reversed_
    assert abs(mc - mp) <= 1e-6 * abs(mp), (mc, mp)


@pytest.mark.parametrize(("dp", "reversed_"), CASES)
def test_combustor_exit_regime_gap_is_the_plenum_gap(dp: float, reversed_: bool) -> None:
    """The regime gap was -21 % forward and +120 to +179 % reversed; it is now
    the plenum exit's (below 4 % at 8 kPa reversed, exit Mach ~0.3)."""
    ps, pe = _ends(dp, reversed_)
    inc = float(_solve(_chain("combustor", "incompressible", ps, pe))["s1x.m_dot"])
    com = float(_solve(_chain("combustor", "compressible", ps, pe))["s1x.m_dot"])
    assert abs((inc - com) / com) < 0.04, (inc, com)


@pytest.mark.parametrize(
    ("node", "at_rest"),
    [
        (PlenumNode("n"), True),
        (CombustorNode("n"), True),
        (PressureBoundary("n", Pt=P0, Tt=300.0), True),
        (MassFlowBoundary("n", m_dot=0.01, Tt=300.0), True),
        (WallNode("n"), True),
        (MomentumChamberNode("n", area=A_DUCT), False),
        (None, True),
    ],
)
def test_terminal_at_rest_is_the_nodes_own_definition(node: object, at_rest: bool) -> None:
    assert CrossflowSegmentElement._terminal_at_rest(node) is at_rest
    if node is not None:
        assert at_rest is bool(node.state_at_rest)


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.parametrize("m", [0.01, -0.01])
def test_combustor_exit_terms_and_terminal(regime: str, m: float) -> None:
    """The momentum terms and the compressible terminal of a combustor exit
    are a plenum exit's: incompressible reversed carries the reservoir
    acceleration (``exit_in``) and no exit head; compressible forward the
    exit end is the closure of the node's static P, not its Pt."""
    for kind in ("combustor", "plenum"):
        g = _chain(kind, regime, P0 + 2000.0, P0)
        for el in g.elements.values():
            el.resolve_topology(g)
        seg = g.elements["s1x"]
        st_in = NetworkMixtureState(1.01e5, 1.01e5, Tt=300.0, m_dot=m, Y=Y)
        st_out = NetworkMixtureState(1.0e5, 1.0e5, Tt=300.0, m_dot=m, Y=Y)
        terms = seg._momentum_terms(st_in, st_out, {"s0x": m})
        kinds = [t[5] for t in terms]
        assert ("exit_in" in kinds) == (regime == "incompressible" and m < 0.0), (kind, kinds)
        assert "exit_head" not in kinds
        _, dst, _, d_out = seg._duct_ends(st_in, st_out, terms)
        if m > 0.0:
            # static continuity: the exit end's Pt is the closure of P, above P
            assert "nd.P" in d_out and "nd.Pt" not in d_out, (kind, d_out)
            assert dst.Pt > st_out.P
        else:
            assert d_out == {"nd.Pt": 1.0}


def test_seed_charges_the_head_into_a_combustor() -> None:
    g = _chain("combustor", "incompressible", P0 + 2000.0, P0)
    for el in g.elements.values():
        el.resolve_topology(g)
    assert g.elements["s1x"]._seed_exit_head() is True


def test_combustor_station_is_refused() -> None:
    """A station node holds the channel's static pressure and its joining
    streams: a CombustorNode station is refused at validate, as before."""
    g = _chain("plenum", "compressible", P0 + 2000.0, P0)
    g.nodes.pop("s1")
    g.add_node(CombustorNode("s1"))
    with pytest.raises(ValueError, match="Station nodes must be PlenumNodes"):
        NetworkSolver(g).solve()


def _fresh_jacobian_error(g: FlowNetwork, x: np.ndarray) -> float:
    # A fresh solver: the residual depends on solver state after a solve (#1142).
    s = NetworkSolver(g)
    s._build_x0()
    J = s._residuals_and_jacobian(x)[1]
    fd = fd_jacobian(s, x, floor=1e-10)
    return float(np.max(scaled_error(J, fd, x)))


@pytest.mark.parametrize("point", ["root", "perturbed", "zero_flow"])
@pytest.mark.parametrize("reversed_", [False, True])
@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_combustor_exit_jacobian(regime: str, reversed_: bool, point: str) -> None:
    ps, pe = _ends(2000.0, reversed_)
    g = _chain("combustor", regime, ps, pe)
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    x = np.array(r["__x_solution__"], dtype=float)
    if point == "perturbed":
        x = x * (1.0 + 1e-3 * np.random.default_rng(7).uniform(-1.0, 1.0, x.size))
    elif point == "zero_flow":
        x[s.unknown_names.index("s1x.m_dot")] = 0.0
    assert _fresh_jacobian_error(_chain("combustor", regime, ps, pe), x) < 1e-5


@pytest.mark.parametrize(("dp", "reversed_"), CASES)
@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_combustor_exit_second_law(dp: float, reversed_: bool, regime: str) -> None:
    ps, pe = _ends(dp, reversed_)
    r = _solve(_chain("combustor", regime, ps, pe))
    assert not r["__second_law_failures__"], r["__second_law__"]
