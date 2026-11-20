"""Second-law check per element and the control-volume ledger (#545).

Every converged solve reports:

* ``__second_law__``: per element, the internal entropy generation of its
  control volume (``sigma_norm``; for an adiabatic two-port exactly
  ``-ln(Pt_out/Pt_in)``) and a verdict -- pass, fail, or not_checked with the
  reason.
* ``__ledger__``: mass and energy over the network's control volume, from
  the flow boundaries, the heat given, and the heat leaving or withheld.

Both demote the solve when they fail: the ledger since #552, an element's
second law since #598.

The networks are the suite's own, one per element kind; the broken element
is #541's area change, which coupled the downstream static pressure and so
created the dynamic head of a downstream chamber as total pressure.
"""

from __future__ import annotations

import dataclasses
import io
import math
import warnings
from contextlib import redirect_stdout

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    AreaChangeElement,
    ChannelElement,
    ConvectiveSurface,
    EnergyBoundary,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)
from combaero.network.components import LazyJacobian
from combaero.network.second_law import SIGMA_TOL

with redirect_stdout(io.StringIO()):
    import test_area_change as area_cases
    import test_combustor_network as combustor_cases
    import test_effusion_liner as liner_cases
    import test_impingement_array as array_cases
    import test_merge_chamber as merge_cases
    import test_tee_network as tee_cases

    from validation.junction import random_robustness as rr

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _solve(net: FlowNetwork, **kw) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return NetworkSolver(net).solve(**kw)


# --- networks, one per element kind -----------------------------------------


def _heated_pair(Q_sign: float = 1.0) -> FlowNetwork:
    """Two channels coupled by a wall; the hot one cools, the cold one heats."""
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("hot_in", m_dot=0.1, Tt=800.0, Y=Y_AIR))
    g.add_node(MassFlowBoundary("cold_in", m_dot=0.05, Tt=400.0, Y=Y_AIR))
    for n in ("hot_out", "cold_out"):
        g.add_node(PressureBoundary(n, Pt=2e5, Tt=300.0, Y=Y_AIR))
    for n in ("hot_p", "cold_p"):
        g.add_node(PlenumNode(n))
    area = math.pi * 0.04
    for name in ("hot", "cold"):
        g.add_element(
            ChannelElement(
                name,
                f"{name}_in",
                f"{name}_p",
                length=1.0,
                diameter=0.04,
                roughness=0.0,
                surface=ConvectiveSurface(area=area),
            )
        )
        g.add_element(
            OrificeElement(
                f"{name}_o",
                f"{name}_p",
                f"{name}_out",
                Cd=0.8,
                diameter=0.0357,
                correlation="fixed",
            )
        )
    g.add_wall(
        ThermalWall(id="w", element_a="hot", element_b="cold", layers=[WallLayer(0.002, 25.0)])
    )
    return g


def _junction() -> FlowNetwork:
    import random

    return rr.build(rr.sample(random.Random(3)))


NETWORKS = {
    "orifices_and_heat": lambda: _heated_pair(),
    "area_change_step": lambda: area_cases._step_into_chamber(1e-3, 4e-4).network,
    "impingement_array": lambda: array_cases._net(array_cases._array(4), hot=True),
    "effusion_liner": lambda: liner_cases._net(liner_cases._liner(), chamber=True),
    "merge_chamber": lambda: merge_cases._liner(),
    "combustor": lambda: combustor_cases._make_simple_combustor_net(0.05)[0],
    "mpce_junction": _junction,
    "legacy_tee": lambda: tee_cases._branching_net(),
}


@pytest.mark.parametrize("name", NETWORKS)
def test_every_element_of_a_working_network_passes(name: str) -> None:
    """No element of a physical solution manufactures flow work, and every
    element gets a verdict -- never a silent gap."""
    net = NETWORKS[name]()
    r = _solve(net)
    assert r["__converged__"], r["__message__"]
    verdicts = r["__second_law__"]
    assert set(verdicts) == set(net.elements)
    assert r["__second_law_failures__"] == []
    for eid, v in verdicts.items():
        if v["verdict"] == "not_checked":
            assert v["reason"], eid
        else:
            assert v["verdict"] == "pass", (eid, v)
            assert v["sigma_norm_low"] <= v["sigma_norm"] <= v["sigma_norm_high"]


@pytest.mark.parametrize("name", NETWORKS)
def test_the_ledger_closes_over_the_network(name: str) -> None:
    net = NETWORKS[name]()
    led = _solve(net)["__ledger__"]
    assert abs(led["mass_closure_rel"]) < 1e-6, led
    assert abs(led["energy_closure_rel"]) < 1e-6, led
    assert led["m_dot_in"] > 0.0 and led["H_in"] != 0.0


def test_a_two_port_reports_its_total_pressure_ratio() -> None:
    r = _solve(area_cases._step_into_chamber(1e-3, 4e-4).network)
    v = r["__second_law__"]["step"]
    assert v["sigma_norm"] == pytest.approx(
        -math.log(r["down.Pt"] / r["up.Pt"]), rel=1e-9, abs=1e-12
    )
    assert r["step.sigma_norm"] == v["sigma_norm"]


def test_a_chain_is_one_control_volume_ending_at_its_boundary() -> None:
    """A segment's momentum spans half of each station, so the chain is the
    control volume; the last segment's exit station is the outlet boundary,
    and the stream into it leaves the chain."""
    net = array_cases._net(array_cases._array(4))
    v = _solve(net)["__second_law__"]
    chain = sorted(e for e in net.elements if e.startswith("ia__x"))
    groups = {v[e]["group"] for e in chain}
    assert groups == {"+".join(chain)}
    assert v[chain[0]]["verdict"] == "pass"


def test_the_legacy_tee_is_checked_as_a_junction() -> None:
    """The legacy tee is a three-port control volume like the #517 junctions
    (#1231): a leg into it at its port node's state, a leg out of it at the
    mix it delivers. A forward branching tee dissipates."""
    v = _solve(tee_cases._branching_net())["__second_law__"]["tee"]
    assert v["verdict"] == "pass"
    assert "reason" not in v
    assert v["sigma_norm"] > 0.0
    assert v["sigma_norm_low"] <= v["sigma_norm"] <= v["sigma_norm_high"]


# --- heat ---------------------------------------------------------------------


def test_heat_counts_where_it_matters() -> None:
    """The heated channel's stream gains entropy from its heat; without the
    heat term that entropy would read as generated in the cold channel and
    destroyed in the hot one."""
    r = _solve(_heated_pair())
    assert r["w.Q"] > 0.0
    v = r["__second_law__"]
    assert v["hot"]["verdict"] == v["cold"]["verdict"] == "pass"


def test_without_the_heat_term_a_heated_channel_would_fail(monkeypatch) -> None:
    """Falsification in the suite: drop the heat's entropy and the cooled
    channel reads as destroying entropy."""
    monkeypatch.setattr(cb, "heat_entropy_flux", lambda *a, **k: 0.0)
    v = _solve(_heated_pair())["__second_law__"]
    assert v["hot"]["verdict"] == "fail"


# --- the broken element -------------------------------------------------------


class _Pre541AreaChange(AreaChangeElement):
    """#541 as it was: the loss coupled to the downstream STATIC pressure, so
    a downstream chamber's dynamic head appeared as free total pressure. In
    the mass-flow form (#731): the flow is driven to the downstream static
    pressure instead of its stagnation pressure."""

    residual_scale_kind = "mdot"  # still the mass-flow row, kg/s

    def residuals(self, state_in, state_out):
        res, jac = super().residuals(state_in, dataclasses.replace(state_out, Pt=state_out.P))
        to_pt, to_p = f"{self.to_node}.Pt", f"{self.to_node}.P"

        def remapped() -> dict[int, dict[str, float]]:
            row = dict(jac[0])
            row[to_p] = row.get(to_p, 0.0) + row.pop(to_pt, 0.0)
            return {0: row}

        return res, LazyJacobian(remapped)


def test_the_541_area_change_is_caught_at_that_element_only() -> None:
    solver = area_cases._step_into_chamber(1e-3, 4e-4)
    net = solver.network
    old = net.elements["step"]
    net.elements["step"] = _Pre541AreaChange(
        "step", old.from_node, old.to_node, F0=old.F0, F1=old.F1, model_type="sharp", length=0.05
    )
    net.resolve_all_topology()
    r = _solve(net)
    assert r["__converged__"]
    assert r["down.Pt"] > r["up.Pt"]  # the defect: total pressure rose across a loss
    assert r["__second_law_failures__"] == ["step"]
    assert r["__second_law__"]["step"]["sigma_norm_high"] < -SIGMA_TOL
    # the root is refused, and the message names the element (#598)
    assert r["__success__"] is False
    assert r["__outcome__"] == "inconsistent"
    assert "second law" in r["__message__"] and "step" in r["__message__"]


# --- the ledger ---------------------------------------------------------------


def test_the_ledger_counts_an_injection_at_its_own_state() -> None:
    """A MassFlowBoundary between two elements injects its stream: that stream,
    at its own Tt and Y, is what crosses the control volume there."""
    Y_inj = list(np.array(Y_AIR) * 0.8)
    Y_inj[cb.species.indices["CH4"]] += 0.2
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=2e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("out", Pt=1e5, Tt=300.0, Y=Y_AIR))
    g.add_node(MassFlowBoundary("inj", m_dot=0.05, Tt=800.0, Y=Y_inj))
    g.add_element(OrificeElement("e1", "a", "inj", Cd=0.8, diameter=0.02, correlation="fixed"))
    g.add_element(OrificeElement("e2", "inj", "out", Cd=0.8, diameter=0.03, correlation="fixed"))
    led = _solve(g)["__ledger__"]
    inj = [row for row in led["streams"] if row["via"] == "inj"]
    assert inj == [
        {
            "node": "inj",
            "via": "inj",
            "m_dot_in": 0.05,
            "H_in": pytest.approx(0.05 * cb.h_mass(800.0, cb.mass_to_mole(Y_inj))),
        }
    ]
    assert abs(led["energy_closure_rel"]) < 1e-9


def test_the_ledger_closes_on_a_plenum_heated_past_the_old_clamp() -> None:
    """100 kW into ~0.5 g/s: 1.2e5 K. A [200, 5000] K clamp threw ~190 kW of
    200 away here and the ledger showed the gap; without it the energy is
    conserved (#626)."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=1.1e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("b", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    p = PlenumNode("p")
    p.add_energy_boundary(EnergyBoundary(id="heat", Q=1e5))
    g.add_node(p)
    g.add_element(OrificeElement("o1", "a", "p", Cd=0.8, diameter=0.01, correlation="fixed"))
    g.add_element(OrificeElement("o2", "p", "b", Cd=0.8, diameter=0.01, correlation="fixed"))
    r = _solve(g, auto_retry=False)
    assert r["__success__"]
    assert r["p.Tt"] > 1e5
    assert abs(r["__ledger__"]["energy_closure_rel"]) < 1e-12


# --- the report never takes the solve down ---------------------------------


def test_a_control_volume_that_cannot_be_evaluated_is_not_checked(monkeypatch) -> None:
    """Stress-tested: a wall taking more heat than a stream holds made the
    heat term raise, and the converged solve raised with it."""

    def refuse(*a, **k):
        raise ValueError("heat_entropy_flux: temperatures must be positive")

    monkeypatch.setattr(cb, "heat_entropy_flux", refuse)
    r = _solve(_heated_pair())
    assert r["__converged__"]
    v = r["__second_law__"]["hot"]
    assert v["verdict"] == "not_checked"
    assert "cannot be evaluated" in v["reason"]
    assert r["__second_law__"]["hot_o"]["verdict"] == "pass"  # the rest still checked


def test_heat_beyond_the_streams_enthalpy_is_refused_by_name() -> None:
    """calc_T_from_h_mass raises for an enthalpy no positive T reaches (#606);
    the check refuses it with the heat and the stream named."""
    from combaero.network.second_law import SecondLawReport, _Stream

    report = SecondLawReport.__new__(SecondLawReport)
    with pytest.raises(ValueError, match="more than the stream's enthalpy"):
        report._heated(_Stream(0.005, 1600.0, 2e5, Y_AIR), -18000.0)


# --- each ThermalWall as its own control volume (#551) -----------------------


def _hot_cold_wall(
    D, L, area_factor, Tt_hot, Tt_cold, pr, layer, cold_backwards: bool = False
) -> FlowNetwork:
    """A hot and a cold compressible channel coupled by a conductive wall --
    the #551 stress topology (tmp/stress/gen.py hot_cold_*), parameters fixed.
    ``cold_backwards`` declares the cold chain co -> cm -> ci, so its flows
    are negative (#627)."""
    from combaero.network import MomentumChamberNode

    g = FlowNetwork()
    P = 3e5
    g.add_node(PressureBoundary("hi", Pt=P, Tt=Tt_hot, Y=Y_AIR))
    g.add_node(PressureBoundary("ci", Pt=P, Tt=Tt_cold, Y=Y_AIR))
    g.add_node(PressureBoundary("ho", Pt=P * pr, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("co", Pt=P * pr, Tt=300.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("hm", area=math.pi * D * D / 4))
    g.add_node(MomentumChamberNode("cm", area=math.pi * D * D / 4))
    area = math.pi * D * L * area_factor
    for s in ("h", "c"):
        ends = (f"{s}m", f"{s}i") if (s == "c" and cold_backwards) else (f"{s}i", f"{s}m")
        g.add_element(
            ChannelElement(
                f"{s}ch",
                *ends,
                length=L,
                diameter=D,
                roughness=0.0,
                regime="compressible",
                surface=ConvectiveSurface(area=area),
            )
        )
        x_ends = (f"{s}o", f"{s}m") if (s == "c" and cold_backwards) else (f"{s}m", f"{s}o")
        g.add_element(
            ChannelElement(
                f"{s}x",
                *x_ends,
                length=0.02,
                diameter=D,
                roughness=0.0,
                regime="compressible",
            )
        )
    g.add_wall(ThermalWall(id="w", element_a="hch", element_b="cch", layers=[WallLayer(*layer)]))
    return g


def test_the_wall_no_longer_cools_a_hot_stream_below_the_cold_inlet():
    """#551 seed 23. The linear wall U A dT_in (the small-NTU limit) moved
    more heat than the streams hold at NTU ~ 1: it cooled the 1600 K stream
    to 298 K against a 400 K cold inlet, and the wall's entropy generation was
    negative on every Mach path. The epsilon-NTU exchange (#627) bounds it: in
    parallel flow both outlets stay between the inlets, and the wall passes."""
    r = _solve(_hot_cold_wall(0.006, 0.2, 5.0, 1600.0, 400.0, 0.7, (0.0005, 200.0)), timeout=20.0)
    assert r["__success__"]
    assert r["w.counterflow"] == 0.0 and 0.9 < r["w.NTU"] < 1.5
    assert 400.0 < r["cm.Tt"] < r["hm.Tt"] < 1600.0
    w = r["__second_law_walls__"]["w"]
    assert w["verdict"] == "pass"
    assert w["sigma_wall"] > 0.0
    assert r["__second_law_wall_failures__"] == []


def test_wall_check_passes_a_physical_heat_exchange():
    """A modest wall (1400 K / 300 K streams, 2 mm at k = 10): heat runs hot to
    cold and the wall generates entropy."""
    r = _solve(_hot_cold_wall(0.012, 1.0, 1.0, 1400.0, 300.0, 0.75, (0.002, 10.0)), timeout=20.0)
    assert r["__success__"]
    w = r["__second_law_walls__"]["w"]
    assert w["verdict"] == "pass"
    assert w["sigma_wall"] > 0.0
    assert r["__second_law_wall_failures__"] == []
