"""``T`` is the static temperature in results and on the state (#782, 0.8).

``<node>.T`` and an element's ``T_in``/``T_out`` held the stagnation T
until 0.7; they are static now, as ``P``/``P_in`` already were, and the
stagnation value is ``<node>.Tt``/``Tt_in``. The static T is regime-aware:
an incompressible node's is its ``Tt`` (the low-Mach limit its closure
solves), a compressible one's lies on the isentrope through (Tt, Pt) at P.

0.8 warned once on the first read of such a key and kept a deprecated
``<chamber>.T_static`` alias; both went in 0.9.0.
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkMixtureState,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)
from combaero.network import solver as solver_module
from combaero.network.components import MultiPortChamberBase

Y_AIR = list(cb.species.dry_air_mass())
X_AIR = list(cb.species.dry_air())


def _chamber_network(regime: str) -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("inlet", m_dot=0.5, Tt=600.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("chamber", area=3e-3, regime=regime))
    g.add_node(PressureBoundary("outlet", Pt=101325.0))
    for eid, a, b in (("orf1", "inlet", "chamber"), ("orf2", "chamber", "outlet")):
        g.add_element(
            OrificeElement(
                eid, a, b, Cd=0.6, diameter=0.04, correlation="fixed", regime="incompressible"
            )
        )
    return g


def test_a_moving_compressible_node_reports_a_static_T_below_its_Tt() -> None:
    sol = NetworkSolver(_chamber_network("compressible")).solve()
    assert sol["__success__"]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        T = sol["chamber.T"]
    Tt, P, Pt = sol["chamber.Tt"], sol["chamber.P"], sol["chamber.Pt"]
    assert Tt - 1.0 > T
    # On the isentrope through (Tt, Pt): the solved closure's static state.
    assert cb.s(T, X_AIR, P) == pytest.approx(cb.s(Tt, X_AIR, Pt), rel=1e-8)
    # The element leaving the chamber reports the chamber's static and total T.
    diag = sol["__element_diag__"]["orf2"]
    assert diag["T_in"] == pytest.approx(T, rel=1e-12)  # isentrope vs closure: round-off
    assert diag["Tt_in"] == Tt


def test_the_result_is_a_plain_dict_without_the_old_T_static_key() -> None:
    """0.9.0: no read-time warning and no "<chamber>.T_static" alias."""
    sol = NetworkSolver(_chamber_network("compressible")).solve()
    assert type(sol) is dict
    assert type(sol["__element_diag__"]["orf2"]) is dict
    assert "chamber.T_static" not in sol
    assert "chamber.T" in sol
    for name in ("_STATIC_T_WARNED", "_StaticTResult", "_warn_static_t", "_warn_jacobian_T_key"):
        assert not hasattr(solver_module, name)


def test_an_incompressible_node_and_a_plenum_report_T_equal_to_Tt() -> None:
    sol = NetworkSolver(_chamber_network("incompressible")).solve()
    assert sol["__success__"]
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        assert sol["chamber.T"] == sol["chamber.Tt"]
        assert sol["inlet.T"] == sol["inlet.Tt"]
    assert sol["chamber.P"] < sol["chamber.Pt"]  # it moves; the regime decides


def test_the_state_computes_its_static_T() -> None:
    plenum = NetworkMixtureState(P=2.0e5, Pt=2.0e5, Tt=600.0, Y=Y_AIR)
    assert plenum.T == 600.0
    moving = NetworkMixtureState(P=1.8e5, Pt=2.0e5, Tt=600.0, Y=Y_AIR)
    assert moving.T < 600.0
    assert cb.s(moving.T, X_AIR, 1.8e5) == pytest.approx(cb.s(600.0, X_AIR, 2.0e5), rel=1e-10)
    low_mach = NetworkMixtureState(P=1.8e5, Pt=2.0e5, Tt=600.0, Y=Y_AIR, incompressible=True)
    assert low_mach.T == 600.0


def test_the_0_7_constructor_and_setter_are_refused() -> None:
    """0.9.0: the 0.7 signature (P, Pt, T, Tt, m_dot, Y) and the T setter
    are gone. Tt, m_dot and Y are keyword-only, so a positional 0.7 call
    fails instead of shifting its arguments."""
    with pytest.raises(TypeError):
        NetworkMixtureState(2.0e5, 2.0e5, 600.0, 600.0, 0.1, Y_AIR)  # type: ignore[misc]
    with pytest.raises(TypeError):
        NetworkMixtureState(P=2.0e5, Pt=2.0e5, T=650.0, Y=Y_AIR)  # type: ignore[call-arg]
    state = NetworkMixtureState(P=2.0e5, Pt=2.0e5, Tt=600.0, m_dot=0.1, Y=Y_AIR)
    with pytest.raises(AttributeError):
        state.T = 610.0  # type: ignore[misc]
    assert state.Tt == 600.0
    with pytest.raises(TypeError, match="Tt"):
        NetworkMixtureState(P=2.0e5, Pt=2.0e5, Y=Y_AIR)  # type: ignore[call-arg]


def test_a_custom_element_keyed_by_the_0_7_name_is_refused(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A residual Jacobian keyed "<node>.T" (the 0.7 name of the stagnation-T
    column) was relayed as "<node>.Tt" with a warning in 0.8; since 0.9.0 it
    raises, naming the new key, rather than dropping the column silently."""
    # Two streams of different Tt mix in the chamber, so its Tt moves with
    # the flow split and the element leaving it has a live Tt column.
    g = FlowNetwork()
    g.add_node(PressureBoundary("hot", Pt=2.0e5, Tt=600.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cold", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("chamber", area=3e-3, regime="compressible"))
    g.add_node(PressureBoundary("outlet", Pt=1.0e5))
    for eid, a, b, d in (
        ("h", "hot", "chamber", 0.03),
        ("c", "cold", "chamber", 0.02),
        ("o", "chamber", "outlet", 0.04),
    ):
        g.add_element(OrificeElement(eid, a, b, Cd=0.6, diameter=d, correlation="fixed"))
    solver = NetworkSolver(g)
    x0 = np.array(solver._build_x0(), dtype=float)
    solver._residuals_and_jacobian(x0)  # the 0.9 names assemble

    def old_names(jac: dict) -> dict:
        return {
            row: {(k[:-1] if k.endswith(".Tt") else k): v for k, v in cols.items()}
            for row, cols in jac.items()
        }

    renamed = 0
    for el in g.elements.values():
        original = el.residuals

        def residuals(*args: Any, _original: Any = original, **kwargs: Any) -> Any:
            nonlocal renamed
            res, jac = _original(*args, **kwargs)
            renamed += sum(k.endswith(".Tt") for cols in jac.values() for k in cols)
            return res, old_names(jac)

        monkeypatch.setattr(el, "residuals", residuals)
    with pytest.raises(KeyError, match=r"(?s)Jacobian key '.*\.T'.*\.Tt"):
        solver._residuals_and_jacobian(x0)
    assert renamed > 0  # the network has a T column to rename


def _tee3(seed: int) -> tuple[Any, dict, str]:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
    from networks import FAMILIES

    net, _ = FAMILIES["tee3"](seed)
    sol = NetworkSolver(net).solve(timeout=30.0)
    assert sol["__converged__"]
    eid = next(k for k, e in net.elements.items() if isinstance(e, MultiPortChamberBase))
    return net, sol, eid


def test_a_junction_reports_each_port_static_and_total_T() -> None:
    """#824: a junction's port_<i>_T is the port's static T (it was the
    total), port_<i>_Tt the total -- the port nodes' own .T and .Tt."""
    net, sol, eid = _tee3(2)  # compressible: the ports move
    d = sol["__element_diag__"][eid]
    for i, port in enumerate(net.elements[eid].port_nodes):
        assert d[f"port_{i}_Tt"] == sol[f"{port}.Tt"]
        # The isentrope through (Tt, Pt) and the chamber closure meet at the
        # converged Pt residual, not to round-off.
        assert d[f"port_{i}_T"] == pytest.approx(sol[f"{port}.T"], rel=1e-10)
        assert d[f"port_{i}_T"] < d[f"port_{i}_Tt"]


def test_an_incompressible_junction_port_T_is_its_Tt() -> None:
    net, sol, eid = _tee3(5)
    d = sol["__element_diag__"][eid]
    for i in range(net.elements[eid].N):
        assert dict.__getitem__(d, f"port_{i}_T") == dict.__getitem__(d, f"port_{i}_Tt")


def test_an_element_reports_its_outlet_flow_at_the_static_state() -> None:
    """#824: the outlet Mach an element reports is formed at its outlet
    node's static state, not at (Tt, P)."""
    sol = NetworkSolver(_chamber_network("compressible")).solve()
    orf = sol["__element_diag__"]["orf1"]  # its outlet is the moving chamber
    T, P = sol["chamber.T"], sol["chamber.P"]
    area = math.pi * 0.04**2 / 4.0
    v = 0.5 / (cb.density(T, P, X_AIR) * area)
    assert orf["mach_out"] == pytest.approx(v / cb.speed_of_sound(T, X_AIR), rel=1e-9)


def _stress(family: str, seed: int) -> tuple[Any, dict]:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
    from networks import FAMILIES

    net, _ = FAMILIES[family](seed)
    sol = NetworkSolver(net).solve(timeout=30.0)
    assert sol["__converged__"]
    return net, sol


def test_a_crossflow_station_reports_its_moving_state() -> None:
    """A crossflow station holds Pt = P (its node is the static P at the
    middle of the station) while its gas moves at the station's mean flow
    (#782 audit). It reports that flow's static T, its true Pt and Mach; the
    segment leaving it reports the same Pt, a plate discharging into it the
    drop to the station's STATIC P its jets see."""
    from combaero.network.components import CrossflowSegmentElement, ImpingementPlateElement

    net, sol = _stress("impingement", 1)
    diag = sol["__element_diag__"]
    stations = {
        e.from_node: eid
        for eid, e in net.elements.items()
        if isinstance(e, CrossflowSegmentElement) and e.prev_seg and e.from_kappa is not None
    }
    assert stations
    for nid, seg in stations.items():
        T, Tt, P, Pt = (sol[f"{nid}.{k}"] for k in ("T", "Tt", "P", "Pt"))
        assert Tt > T and Pt > P and sol[f"{nid}.mach"] > 0.05
        assert cb.s(T, X_AIR, P) == pytest.approx(cb.s(Tt, X_AIR, Pt), rel=1e-9)
        assert diag[seg]["Pt_in"] == Pt
        assert diag[seg]["T_in"] == T
    plates = [e for e in net.elements.values() if isinstance(e, ImpingementPlateElement)]
    for pl in plates:
        if pl.to_node in stations:
            assert diag[pl.id]["Pt_out"] == sol[f"{pl.to_node}.P"]


def test_a_closed_end_station_reports_its_moving_state() -> None:
    """Row 1 of an impingement gap is a station with nothing arriving (a
    closed end): its gas moves at half the leaving segment's flow, the mean
    its half-drop uses. It reports that moving state like every other
    station; it reported Tt, its P and Mach 0 before (B1)."""
    from combaero.network.components import CrossflowSegmentElement

    net, sol = _stress("impingement", 1)
    diag = sol["__element_diag__"]
    closed = [
        e
        for e in net.elements.values()
        if isinstance(e, CrossflowSegmentElement) and e.from_kappa is not None and not e.prev_seg
    ]
    assert closed
    for seg in closed:
        nid = seg.from_node
        T, Tt, P, Pt, M = (sol[f"{nid}.{k}"] for k in ("T", "Tt", "P", "Pt", "mach"))
        m = 0.5 * sol[f"{seg.id}.m_dot"]
        rho = cb.density(T, P, X_AIR)
        assert pytest.approx(abs(m) / (rho * seg.area_from) / cb.speed_of_sound(T, X_AIR)) == M
        assert Tt > T and Pt > P and M > 0.0
        assert cb.s(T, X_AIR, P) == pytest.approx(cb.s(Tt, X_AIR, Pt), rel=1e-9)
        assert diag[seg.id]["Pt_in"] == Pt
        assert diag[seg.id]["T_in"] == T


def test_a_to_only_station_reports_its_moving_state() -> None:
    """A segment's to-station whose next segment does not take it as a
    from-station (B1): the station of the arriving and the next segment,
    reported moving at their mean flow, the state the arriving segment's
    half-drop and a duct-fed panel there read. It reported Tt, its P and 0
    before."""
    from combaero.network import CrossflowSegmentElement, EffusionPlateElement, HoleCd

    kappa = cb.STATION_KAPPA_BLEED_BASSETT
    g = FlowNetwork()
    for n, p in (("in", 1.3e5), ("out", 1.2e5), ("gas", 1.0e5)):
        g.add_node(PressureBoundary(n, Pt=p, Tt=290.0, Y=Y_AIR))
    from combaero.network import PlenumNode

    g.add_node(PlenumNode("n1"))
    g.add_element(
        CrossflowSegmentElement(
            "s0", "in", "n1", length=0.03, area=3e-4, entry_K=0.0, to_kappa=kappa, next_seg="s1"
        )
    )
    g.add_element(CrossflowSegmentElement("s1", "n1", "out", length=0.03, area=3e-4))
    g.add_element(
        EffusionPlateElement(
            "h", "n1", "gas", hole_diameter=1e-3, wall_thickness=3e-3, pitch_x=0.06,
            pitch_y=0.06, panel_length=0.06, panel_width=0.06, angle_deg=30.0,
            hole_cd=HoleCd(correlation="McGreehanSchotsch"), crossflow_segments=("s0", "s1"),
        )
    )  # fmt: skip
    sol = NetworkSolver(g).solve(timeout=60.0)
    assert sol["__success__"]
    T, Tt, P, Pt, M = (sol[f"n1.{k}"] for k in ("T", "Tt", "P", "Pt", "mach"))
    m = 0.5 * (sol["s0.m_dot"] + sol["s1.m_dot"])
    rho = cb.density(T, P, X_AIR)
    assert pytest.approx(abs(m) / (rho * 3e-4) / cb.speed_of_sound(T, X_AIR), rel=1e-9) == M
    assert M > 0.1 and Tt > T and Pt > P
    assert cb.s(T, X_AIR, P) == pytest.approx(cb.s(Tt, X_AIR, Pt), rel=1e-9)
    assert sol["__element_diag__"]["s0"]["Pt_out"] == Pt


def test_an_area_combustor_reports_its_static_T() -> None:
    """A combustor with an area holds Pt = P but its gas moves through the
    area: it reports the duct closure's static T at G = m / area (#782
    audit), and its Mach at that state."""
    from combaero.network import ChannelElement, CombustorNode

    Yf = [0.0] * len(Y_AIR)
    Yf[cb.species_index_from_name("CH4")] = 1.0
    g = FlowNetwork()
    g.add_node(PressureBoundary("ai", Pt=1.6e5, Tt=700.0, Y=Y_AIR))
    g.add_node(PressureBoundary("fi", Pt=1.6e5, Tt=300.0, Y=Yf))
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(CombustorNode("hn", area=1e-3))
    g.add_element(ChannelElement("oa", "ai", "hn", length=0.3, diameter=0.02))
    g.add_element(ChannelElement("of", "fi", "hn", length=0.3, diameter=0.003))
    g.add_element(ChannelElement("hout", "hn", "ho", length=0.4, diameter=0.02))
    sol = NetworkSolver(g).solve()
    assert sol["__converged__"]
    Tt, Pt = sol["hn.Tt"], sol["hn.Pt"]
    m = sol["oa.m_dot"] + sol["of.m_dot"]
    X = [sol[f"hn.X[{k}]"] for k in range(len(Y_AIR))]
    s = cb.duct_htc_static_state(Tt, Pt, m / 1e-3, X)
    assert sol["hn.T"] == pytest.approx(float(s.T), rel=1e-12)
    assert sol["hn.T"] < Tt
