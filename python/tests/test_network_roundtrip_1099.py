"""FlowNetwork.to_dict -> from_dict reproduces the network bit for bit (#1099).

to_dict reads each constructor argument back from the attribute of the same
name. Where resolve_topology rewrites that attribute (an inherited area, a
selected compressible model, the crossflow gap's resolved entry_K) the
resolved value was serialised as if the user had given it, or, for the
impingement gap, None was passed back into a float argument and the load
failed. Each class now names the attribute that holds the argument
(``_init_kwarg_attrs``).

The check: the reloaded network's residuals and Jacobian at the original
network's x0 are bit-identical, over every stress family and the impingement
configurations of the issue.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    EnergyBoundary,
    FlowNetwork,
    ImpingementArray,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
from networks import FAMILIES  # noqa: E402

D = 0.00254
SPAN = 12 * 4.0 * D
Y_AIR = cb.mole_to_mass(cb.species.dry_air())


def _pb(name: str, Pt: float, Tt: float = 300.0) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y_AIR
    return b


def _array(id: str = "ia", n_rows: int = 5, **kw) -> ImpingementArray:
    return ImpingementArray(
        id, n_rows=n_rows, d_jet=D, xn_d=5.0, yn_d=4.0, z_d=2.0, span=SPAN, plate_thickness=D, **kw
    )


def _fj(solver: NetworkSolver, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    F, J = solver._residuals_and_jacobian(x)
    Jd = J.toarray() if hasattr(J, "toarray") else np.asarray(J)
    return np.asarray(F, dtype=float), np.asarray(Jd, dtype=float)


def _assert_round_trip(net: FlowNetwork, resolved_first: bool = False) -> FlowNetwork:
    """Reload ``net`` from its dict; F and J at its x0 must be bit-identical."""
    if resolved_first:
        net.resolve_all_topology()
    net2 = FlowNetwork.from_dict(net.to_dict())
    net.resolve_all_topology()
    net2.resolve_all_topology()
    s1, s2 = NetworkSolver(net), NetworkSolver(net2)
    x0 = np.array(s1._build_x0(), dtype=float)
    s2._build_x0()
    assert s2.unknown_names == s1.unknown_names
    rng = np.random.default_rng(0)
    for x in (x0, x0 * (1.0 + 0.02 * rng.standard_normal(x0.size))):
        (F1, J1), (F2, J2) = _fj(s1, x), _fj(s2, x)
        assert np.array_equal(F1, F2, equal_nan=True)
        assert np.array_equal(J1, J2, equal_nan=True)
    return net2


# --- the issue's configurations ----------------------------------------------


def _plain() -> FlowNetwork:
    net = FlowNetwork()
    net.add_node(_pb("supply", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    _array().add_to(net, "supply", "exit")
    return net


def _bypass(entry_K: float = 0.5) -> FlowNetwork:
    net = FlowNetwork()
    net.add_node(_pb("supply", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    net.add_node(MassFlowBoundary("feed_in", m_dot=0.002, Tt=300.0, Y=Y_AIR))
    net.add_node(MomentumChamberNode("byp"))
    net.add_element(ChannelElement("feed", "feed_in", "byp", length=0.2, diameter=0.02))
    _array(bypass_entry_K=entry_K).add_to(net, "supply", "exit", bypass="byp")
    return net


def _series() -> FlowNetwork:
    net = FlowNetwork()
    net.add_node(_pb("supply", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    net.add_node(PlenumNode("joint"))
    _array("a", n_rows=4).add_to(net, "supply", "joint")
    _array("b", n_rows=4, bypass_entry_K=0.5).add_to(net, "supply", "exit", bypass="joint")
    return net


@pytest.mark.parametrize("resolved_first", [False, True])
@pytest.mark.parametrize("build", [_plain, _bypass, _series], ids=["plain", "bypass", "series"])
def test_impingement_round_trip(build, resolved_first: bool) -> None:
    _assert_round_trip(build(), resolved_first)


@pytest.mark.parametrize("resolved_first", [False, True])
def test_gap_serialises_the_entry_K_argument(resolved_first: bool) -> None:
    """The gap's argument, not the parent's resolved entry_K: None for a jet
    segment and for a bypass inlet continuing an upstream gap."""
    net = _series()
    if resolved_first:
        net.resolve_all_topology()
    kw = net.to_dict()["elements"]
    assert kw["a__x1"]["kwargs"]["entry_K"] == 0.0
    assert kw["b__x0"]["kwargs"]["bypass_inlet"] is True
    assert kw["b__x0"]["kwargs"]["entry_K"] == 0.5
    net2 = _assert_round_trip(net)
    lead = net2.elements["b__x0"]
    # It continues gap a, so the parent's entry_K resolves to None again.
    assert lead.bypass_entry_K == 0.5 and lead.entry_K is None and lead.prev_seg == "a__x4"


def test_inherited_geometry_is_serialised_as_not_given() -> None:
    """An area resolve_topology inherited stays inherited after a reload: it
    is not frozen into a user value (nor a 0.0 placeholder before resolve)."""
    net = FlowNetwork()
    net.add_node(_pb("in", 2.0e5))
    net.add_node(_pb("out", 1.0e5))
    net.add_node(MomentumChamberNode("mc"))
    net.add_element(ChannelElement("c1", "in", "mc", length=0.5, diameter=0.03))
    net.add_element(ChannelElement("c2", "mc", "out", length=0.5))
    for resolved in (False, True):
        if resolved:
            net.resolve_all_topology()
        d = net.to_dict()
        assert d["nodes"]["mc"]["kwargs"]["area"] is None
        assert d["nodes"]["mc"]["kwargs"]["Dh"] is None
        assert d["elements"]["c2"]["kwargs"]["diameter"] is None
        assert d["elements"]["c2"]["kwargs"]["Dh"] is None
        assert d["elements"]["c1"]["kwargs"]["diameter"] == 0.03
    net2 = _assert_round_trip(net)
    assert net2.nodes["mc"]._area_source == net.nodes["mc"]._area_source != "user"


def test_energy_boundaries_are_serialised() -> None:
    """Attached after construction, so not among the node's arguments."""
    net = FlowNetwork()
    net.add_node(_pb("in", 2.0e5))
    net.add_node(_pb("out", 1.0e5))
    net.add_node(PlenumNode("p"))
    net.nodes["p"].add_energy_boundary(EnergyBoundary("q", Q=-2000.0))
    net.add_element(OrificeElement("o1", "in", "p", Cd=0.7, diameter=0.01, correlation="fixed"))
    net.add_element(OrificeElement("o2", "p", "out", Cd=0.7, diameter=0.01, correlation="fixed"))
    net2 = _assert_round_trip(net)
    (eb,) = net2.nodes["p"].energy_boundaries
    assert (eb.id, eb.Q, eb.fraction, eb.Tt_target) == ("q", -2000.0, 0.0, None)


# --- every stress family ------------------------------------------------------


def _mutated_bore(net: FlowNetwork) -> bool:
    """The random family shrinks some bores after construction, setting area
    and diameter apart from pi d^2 / 4 (_fit_orifices_into_their_chambers):
    a state no constructor argument reproduces bit for bit."""
    return any(
        type(e) is OrificeElement
        and e.diameter is not None
        and e.area != math.pi * (e.diameter / 2.0) ** 2
        for e in net.elements.values()
    )


@pytest.mark.parametrize("seed", range(6))
@pytest.mark.parametrize("family", sorted(FAMILIES))
def test_stress_family_round_trip(family: str, seed: int) -> None:
    net, _ = FAMILIES[family](seed)
    if _mutated_bore(net):
        pytest.skip("bore set after construction, not a constructor state")
    try:
        net.resolve_all_topology()
        NetworkSolver(net)._build_x0()
    except Exception as e:  # noqa: BLE001
        pytest.skip(f"draw does not build: {type(e).__name__}")
    _assert_round_trip(net)


@pytest.mark.xfail(
    strict=True,
    reason="#1099 follow-up: a node's initial_guess is solver state set after "
    "construction and is not serialised; F/J agree, x0 does not",
)
def test_initial_guess_is_serialised() -> None:
    net, _ = FAMILIES["ejector"](0)
    net.resolve_all_topology()
    net2 = FlowNetwork.from_dict(net.to_dict())
    net2.resolve_all_topology()
    x1 = np.array(NetworkSolver(net)._build_x0(), dtype=float)
    x2 = np.array(NetworkSolver(net2)._build_x0(), dtype=float)
    assert np.array_equal(x1, x2)
