"""A crossflow chain's exit into a MomentumChamberNode carries its head.

A chain's exit segment (no to-station) into a moving chamber is a terminal
node: compressible, the duct's exit TOTAL meets the chamber's Pt
(``_duct_ends``, #802). Incompressible, the station static met the chamber's
Pt instead -- the duct's head lost, the chamber's invented -- so the regimes
differed O(q), 4.7 % in flow at exit Mach 0.07 against 0.8 % on a plenum
exit. Now the incompressible row carries the exit head
m^2 / (2 rho A_to^2) (``_chamber_exit_head``): the Ma -> 0 limit of the
compressible one, in both directions. The plenum and boundary exits, which
lose the head by convention, are untouched.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian, scaled_error

import combaero as cb
from combaero.network import (
    ChannelElement,
    CrossflowSegmentElement,
    FlowNetwork,
    ImpingementArray,
    MomentumChamberNode,
    NetworkMixtureState,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

D = 0.00254
Y = list(cb.mole_to_mass(cb.species.dry_air()))
SPAN = 0.122
Z_D = 2.0
A_GAP = Z_D * D * SPAN
A_DUCT = 2.0e-4
BLEED = cb.STATION_KAPPA_BLEED_BASSETT


def _pb(name: str, Pt: float) -> PressureBoundary:
    return PressureBoundary(name, Pt=Pt, Tt=300.0, Y=Y)


def _impingement(kind: str, regime: str, dp: float) -> FlowNetwork:
    """A 5-row impingement array whose crossflow exits into ``kind`` ('chamber'
    or 'plenum'), then an orifice to the exit."""
    g = FlowNetwork()
    g.add_node(_pb("supply", 1.0e5 + dp))
    g.add_node(_pb("exit", 1.0e5))
    arr = ImpingementArray(
        "ia", n_rows=5, d_jet=D, xn_d=5.0, yn_d=4.0, z_d=Z_D, span=SPAN, plate_thickness=D
    )
    node = (
        MomentumChamberNode("ch", area=A_GAP, regime=regime)
        if kind == "chamber"
        else PlenumNode("ch")
    )
    g.add_node(node)
    arr.add_to(g, "supply", "ch")
    d_out = 0.9 * math.sqrt(4.0 * A_GAP / math.pi)
    g.add_element(OrificeElement("out", "ch", "exit", Cd=0.95, correlation="fixed", diameter=d_out))
    for e in g.elements.values():
        if hasattr(e, "regime"):
            e.regime = regime
    return g


def _bleed_chain(regime: str, p_sup: float, p_exit: float, entry_only: bool = False) -> FlowNetwork:
    """A reservoir entry, one bleed station, and an exit segment into a
    chamber that a channel drains. ``entry_only``: one entry segment straight
    into the chamber (no station)."""
    g = FlowNetwork()
    g.add_node(_pb("sup", p_sup))
    g.add_node(_pb("exit", p_exit))
    g.add_node(MomentumChamberNode("mc", area=A_DUCT, regime=regime))
    D_duct = math.sqrt(4.0 * A_DUCT / math.pi)
    if entry_only:
        g.add_element(
            CrossflowSegmentElement(
                "s1x", "sup", "mc", length=0.05, area=A_DUCT, entry_K=0.5, regime=regime
            )
        )
    else:
        g.add_node(_pb("sink", 1.0e5))
        g.add_node(PlenumNode("s1"))
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
                "mc",
                length=0.05,
                area=A_DUCT,
                from_kappa=BLEED,
                prev_seg="s0x",
                regime=regime,
            )
        )
        g.add_element(
            OrificeElement("bl", "s1", "sink", diameter=0.006, Cd=0.6, correlation="fixed")
        )
    g.add_element(ChannelElement("tail", "mc", "exit", length=0.05, diameter=D_duct, regime=regime))
    return g


def _flow(g: FlowNetwork, eid: str) -> float:
    r = NetworkSolver(g).solve()
    assert r["__success__"], r.get("__message__")
    return float(r[f"{eid}.m_dot"])


def _gap(build, eid: str) -> float:
    inc = _flow(build("incompressible"), eid)
    com = _flow(build("compressible"), eid)
    return (inc - com) / com


@pytest.mark.parametrize("dp", [1508.0, 6031.0])  # exit Mach ~0.05, ~0.10
def test_chamber_exit_regimes_agree_like_a_plenum_exit(dp: float) -> None:
    """The regime gap on a chamber exit is the plenum exit's, O(Ma^2): on main
    it was 4.7 % at Mach 0.07, the lost head."""
    ch = _gap(lambda r: _impingement("chamber", r, dp), "ia__x5")
    pl = _gap(lambda r: _impingement("plenum", r, dp), "ia__x5")
    assert abs(ch) < 1.35 * abs(pl) + 5e-4, (ch, pl)


def test_chamber_exit_gap_scales_as_mach_squared() -> None:
    """Doubling the exit Mach quadruples the gap (O(Ma^2)); an O(q) gap,
    relative to an O(q) drive, does not fall with Mach."""
    lo = _gap(lambda r: _impingement("chamber", r, 1508.0), "ia__x5")
    hi = _gap(lambda r: _impingement("chamber", r, 6031.0), "ia__x5")
    assert 3.0 < hi / lo < 5.0, (lo, hi)
    assert abs(lo) < 6e-3


@pytest.mark.parametrize(
    ("p_sup", "p_exit", "entry_only"),
    [
        (1.02e5, 1.0e5, False),
        (1.0e5, 1.02e5, False),  # reversed: the chamber is a moving reservoir
        (1.02e5, 1.0e5, True),
        (1.0e5, 1.02e5, True),
    ],
)
def test_bleed_chain_into_chamber_regimes_agree(
    p_sup: float, p_exit: float, entry_only: bool
) -> None:
    """Forward and reversed, through a station or a reservoir entry: the duct
    end's total is the chamber's Pt in both regimes."""
    g = _gap(lambda r: _bleed_chain(r, p_sup, p_exit, entry_only), "s1x")
    if p_exit > p_sup:
        assert _flow(_bleed_chain("incompressible", p_sup, p_exit, entry_only), "s1x") < 0.0
    assert abs(g) < 5e-3, g


def _check_jacobian(g: FlowNetwork, perturb: float) -> None:
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    x = np.array(r["__x_solution__"], dtype=float)
    if perturb:
        rng = np.random.default_rng(7)
        x = x * (1.0 + perturb * rng.uniform(-1.0, 1.0, x.size))
    J = s._residuals_and_jacobian(x)[1]
    fd = fd_jacobian(s, x, floor=1e-10)
    assert float(np.max(scaled_error(J, fd, x))) < 1e-5


@pytest.mark.parametrize("perturb", [0.0, 1e-3])
@pytest.mark.parametrize(
    ("p_sup", "p_exit", "entry_only"),
    [
        (1.02e5, 1.0e5, False),
        (1.0e5, 1.02e5, False),
        (1.02e5, 1.0e5, True),
        (1.0e5, 1.02e5, True),
    ],
)
def test_chamber_exit_head_jacobian(
    p_sup: float, p_exit: float, entry_only: bool, perturb: float
) -> None:
    _check_jacobian(_bleed_chain("incompressible", p_sup, p_exit, entry_only), perturb)


@pytest.mark.parametrize("perturb", [0.0, 1e-3])
def test_impingement_chamber_exit_jacobian(perturb: float) -> None:
    _check_jacobian(_impingement("chamber", "incompressible", 2955.0), perturb)


@pytest.mark.parametrize("kind", ["plenum", "boundary", "chamber"])
@pytest.mark.parametrize("m", [0.01, -0.01])
def test_only_a_chamber_exit_carries_the_head(kind: str, m: float) -> None:
    """A plenum or boundary exit keeps static continuity: no exit-head term,
    either direction. A chamber exit carries it, incompressible only."""
    g = FlowNetwork()
    g.add_node(_pb("sup", 1.02e5))
    g.add_node(PlenumNode("s1"))
    if kind == "plenum":
        g.add_node(PlenumNode("out"))
    elif kind == "boundary":
        g.add_node(_pb("out", 1.0e5))
    else:
        g.add_node(MomentumChamberNode("out", area=A_DUCT))
    g.add_element(ChannelElement("feed", "sup", "s1", length=0.05, diameter=0.016))
    for regime in ("incompressible", "compressible"):
        seg = CrossflowSegmentElement(
            f"s_{regime}", "s1", "out", length=0.05, area=A_DUCT, from_kappa=BLEED, regime=regime
        )
        g.add_element(seg)
        seg.resolve_topology(g)
        st = NetworkMixtureState(1.01e5, 1.01e5, Tt=300.0, m_dot=m, Y=Y)
        kinds = [t[5] for t in seg._momentum_terms(st, st, {})]
        expected = kind == "chamber" and regime == "incompressible"
        assert ("exit_head" in kinds) == expected, kinds


def _hot_merge(regime: str, p_sup: float, p_exit: float, T_hot: float = 900.0) -> FlowNetwork:
    """One entry segment (no station) into a merge chamber that also takes a
    hot side stream: the chamber's gas (388 K mixed) is not the segment's."""
    g = FlowNetwork()
    g.add_node(_pb("sup", p_sup))
    g.add_node(PressureBoundary("hot", Pt=p_sup, Tt=T_hot, Y=Y))
    g.add_node(_pb("exit", p_exit))
    g.add_node(MomentumChamberNode("mc", area=A_DUCT, regime=regime, main_inlet="s1x"))
    D_duct = math.sqrt(4.0 * A_DUCT / math.pi)
    g.add_element(
        CrossflowSegmentElement(
            "s1x", "sup", "mc", length=0.05, area=A_DUCT, entry_K=0.5, regime=regime
        )
    )
    g.add_element(OrificeElement("h", "hot", "mc", diameter=0.008, Cd=0.6, correlation="fixed"))
    g.add_element(
        ChannelElement("tail", "mc", "exit", length=0.05, diameter=1.3 * D_duct, regime=regime)
    )
    return g


@pytest.mark.parametrize("dp", [500.0, 2000.0, 8000.0])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_hot_merge_chamber_behind_an_entry_regimes_agree(dp: float) -> None:
    """Behind a reservoir entry the exit head is taken at the entry term's
    density, so (1 + K) q - q cancels to the entry loss K q the compressible
    duct ends carry. Both take the carried gas's, the supply's (#1137): the
    gap is then the one-gas chain's O(Ma^2) gap (side stream at the supply's
    300 K), +0.059 / +0.145 % against +0.067 / +0.181 % at 2 / 8 kPa.

    With only the exit head on the carried gas it was -7.9 %. With both on
    the chamber's mixed gas (#1139) the gap drifted with the side stream's
    T (+0.18, +0.09, +0.02, -0.11 % at 300, 600, 900, 1500 K, 8 kPa): its
    0.02 % at 900 K was a crossing, not agreement."""
    g = _gap(lambda r: _hot_merge(r, 1.0e5 + dp, 1.0e5), "s1x")
    same = _gap(lambda r: _hot_merge(r, 1.0e5 + dp, 1.0e5, T_hot=300.0), "s1x")
    assert abs(g) < 2e-3, g
    assert abs(g - same) < 0.25 * abs(same) + 1e-5, (g, same)


@pytest.mark.parametrize("perturb", [0.0, 1e-3])
@pytest.mark.parametrize(("p_sup", "p_exit"), [(1.02e5, 1.0e5), (1.0e5, 1.02e5)])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_hot_merge_chamber_exit_jacobian(p_sup: float, p_exit: float, perturb: float) -> None:
    _check_jacobian(_hot_merge("incompressible", p_sup, p_exit), perturb)


@pytest.mark.parametrize("kind", ["plenum", "boundary", "chamber", "station"])
def test_seed_charges_the_head_only_where_the_exit_loses_it(kind: str) -> None:
    """The start guess charges one exit head wherever the exit loses it, not
    into a MomentumChamberNode, which takes it with the exit total."""
    g = _bleed_chain("incompressible", 1.02e5, 1.0e5)
    if kind != "chamber":
        g = FlowNetwork()
        g.add_node(_pb("sup", 1.02e5))
        g.add_node(PlenumNode("s1"))
        g.add_node(_pb("out", 1.0e5) if kind == "boundary" else PlenumNode("out"))
        g.add_element(
            CrossflowSegmentElement(
                "s1x",
                "s1",
                "out",
                length=0.05,
                area=A_DUCT,
                from_kappa=BLEED,
                to_kappa=BLEED if kind == "station" else None,
                next_seg="s2x" if kind == "station" else None,
            )
        )
        if kind == "station":
            g.add_element(
                CrossflowSegmentElement(
                    "s2x", "out", "sup", length=0.05, area=A_DUCT, from_kappa=BLEED, prev_seg="s1x"
                )
            )
    for el in g.elements.values():
        el.resolve_topology(g)
    seg = g.elements["s1x"]
    assert seg._seed_exit_head() == (kind != "chamber")


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_the_seed_change_moves_no_root(regime: str, monkeypatch: pytest.MonkeyPatch) -> None:
    """Seed only: the old start guess (head always charged) and the new one
    reach the same root to solver tolerance."""
    builds = [
        lambda: _bleed_chain(regime, 1.02e5, 1.0e5),
        lambda: _bleed_chain(regime, 1.0e5, 1.02e5),
        lambda: _hot_merge(regime, 1.02e5, 1.0e5),
    ]
    new = CrossflowSegmentElement._seed_exit_head
    for build in builds:
        xs = []
        for seed in (lambda self: True, new):
            monkeypatch.setattr(CrossflowSegmentElement, "_seed_exit_head", seed)
            r = NetworkSolver(build()).solve()
            assert r["__success__"], r.get("__message__")
            xs.append(np.array(r["__x_solution__"], dtype=float))
        rel = np.abs(xs[0] - xs[1]) / np.maximum(np.abs(xs[0]), 1e-12)
        assert float(np.max(rel)) < 1e-5
