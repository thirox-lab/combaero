"""Jacobian coverage sweep: the analytic network Jacobian against finite
differences where the hand-written tests do not go.

Every recent Jacobian bug sat in a region no test visited: a legacy tee's
composition partials with a reversed leg (#1233) or a side stream into its
common node (#1235), an equilibrium combustor at a zero-flow fuel stream
(#1241), the mixer's near-zero band (#1238). The stress families
(``scripts/stress/networks.py``) are mostly one gas at fixed temperatures,
where composition and temperature partials are exact zeros. This sweep takes
their networks, and a few builders of its own for the classes they never draw
(``EXTRA``), and

* randomises every boundary's gas (``GASES``) and Tt (``TEMPS``) and every
  combustor's method, per variant (variant 0 is the family's own draw);
* injects side streams (a reservoir of another gas and temperature, through
  an orifice) into interior nodes, junction ports first;
* evaluates F and J at many points (``points``): the start vector, seeded
  perturbations of it, every sign pattern of each element's legs at two
  magnitude ratios (so a derived leg, a tee's straight run, reverses too; a
  MultiPortChamberBase junction's legs are its arms), one element's flows
  near zero, one flow at exactly zero, and every flow at a node near zero on
  both sides of the mixer's no-flow band (``kMixerNoFlowBand``).

Each point gets a fresh, resolved solver: a fresh NetworkSolver on an
unresolved network is a different model, and a reused one keeps warm starts.

Comparison
----------
Columns are compared in scaled variables (``scaled_error``: per row, max_j
|J_ij - fd_ij| s_j over max_j |fd_ij| s_j), with one change: a flow column's
scale s_j is the larger of |x_j| and its start value, so a near-zero flow's
column still counts (at 1e-9 kg/s, J_ij |x_j| would hide any error in it).

A cheap screen (``fd_jacobian(adaptive=False)``, flow steps ``_flow_step``)
flags rows above ``TOL``; ``_confirm`` re-differences them: pressure columns
on ``fd_jacobian``'s drive-aware ladder (retried from a 100x larger step when
its estimate does not resolve the entry), every other column on a Richardson
ladder over six decades of step (``_ladder``, the window ``_pick`` chooses);
a near-zero flow one-sided away from zero within |x| of the point, an
exactly-zero one on both sides from 1e-6 of its scale down, each row taking
the side it agrees with.
A row is a **mismatch** only when the FD's own estimate (its window's spread,
and never below the round-off ``ROUNDOFF |F_i| / h``) is below a tenth of
the miss and within ``FD_RESOLVED`` of its value; otherwise it is
**unresolved**, counted per element and never silently passed.

Exclusion rules
---------------
Each is a point or entry without a derivative to check (none is a
tolerance); every excluded point is counted per rule in ``Report.excluded``
(``zero_between`` points are never built).

``exact_zero``     a leg exactly 0 at a point not built for it: a derived
                   leg (a tee's straight run at m_com = m_branch) is a kink,
                   the stream's direction and the mixer's arriving set switch
                   there. Near-zero points use distinct magnitudes per
                   unknown so a derived leg stays off zero.
``zero_between``   an exactly-zero flow between two interior nodes: at zero
                   the stream counts as arriving (m >= 0) in the mix of both
                   ends, so J is by construction a hybrid of the two sides.
                   An exactly-zero flow with one boundary end has one
                   interior mix, and J must match the one-sided FD on the
                   side it agrees with (``exact_zero`` points, #1241).
``in_band``        a node whose arriving flow is inside kMixerNoFlowBand has
                   its mixer partials set to 0 by design (#1232), while its
                   state is still the flow-weighted mean: FD there measures
                   the design choice. Evaluated (F, J finite), not compared.
``nonfinite``      F or J not finite, or an exception, at the point or on a
                   stencil: counted per family and exception.
A ladder level is used only while its whole stencil keeps every leg's sign
(``_signature``), so no difference crosses a sign change, own or derived;
pressure columns keep a duct's stencil on its side of the Fanno drive floor.

Mismatches open on main are listed in ``OPEN``, each with a reproducer; the
slow test tier fails on any other. Run as a script for the coverage table and
the mismatch list over every seed and variant::

    uv run python python/tests/_jacobian_sweep.py --out tmp/jacsweep.json
"""

from __future__ import annotations

import argparse
import contextlib
import json
import math
import random
import sys
import time
import warnings
import zlib
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parents[1] / "scripts" / "stress"))

from _fd_jacobian import (  # noqa: E402
    FDJacobian,
    _evaluate,
    _pick,
    fd_jacobian,
    is_pressure,
    scaled_error,
    worst,
)
from networks import FAMILIES  # noqa: E402

import combaero as cb  # noqa: E402
import combaero.network as cbn  # noqa: E402
from combaero.network import (  # noqa: E402
    ChannelElement,
    CombustorNode,
    CompressorElement,
    FilmCoolingElement,
    FlowNetwork,
    FuelBoundary,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PortDischargeElement,
    PressureBoundary,
    SlotFilmElement,
    TeeJunctionElement,
)
from combaero.network.ejector_element import EjectorElement  # noqa: E402

#: The mixer's near-zero band [kg/s] (kMixerNoFlowBand, solver_interface.h).
MIXER_NO_FLOW_BAND = 1e-12
#: A scaled row error above this is a miss (the FD tests' usual bound).
TOL = 1e-6
#: Near-zero magnitudes [kg/s]: inside the band (summed over at most 9
#: streams of 1-4x it stays inside) and three decades outside it, both far
#: from the heat floor kMixerHeatMdotFloor (1e-6, a C1 joint).
NEAR_ZERO_IN_BAND = 1e-13
NEAR_ZERO_OUT_BAND = 1e-9
#: One element's flows near zero: outside the band, and a small real flow.
ELEMENT_NEAR_ZERO = (1e-9, 3e-5)
#: Magnitude ratios on an element's last leg with its sign patterns: a tee's
#: straight run is m_com - m_branch, reversed only when the branch outweighs
#: the common leg.
RATIOS = (1.0, 3.0)
#: Sign patterns per element at most (a junction with many ports).
MAX_PATTERNS = 8
#: Points per case at most: the start and its perturbations, then the
#: mandatory points (every element's first reversal, every boundary-fed flow
#: at zero at a mixer, every node near zero outside the band), then a seeded
#: sample of the rest. Below the mandatory count, a sample of those.
MAX_POINTS = 40

_Y_AIR = list(cb.species.dry_air_mass())


def _species_mix(**mass: float) -> list[float]:
    Y = [0.0] * len(_Y_AIR)
    for k, v in mass.items():
        Y[cb.species.indices[k]] = v
    s = sum(Y)
    return [y / s for y in Y]


def _with_ch4(frac: float) -> list[float]:
    i = cb.species.indices["CH4"]
    return [(1.0 - frac) * a + (frac if k == i else 0.0) for k, a in enumerate(_Y_AIR)]


#: Boundary gases: air, a lean fuel-air premix, an N2-rich mix, vitiated
#: (burnt) air and a rich fuel-air mix.
GASES: dict[str, list[float]] = {
    "air": _Y_AIR,
    "air_ch4": _with_ch4(0.04),
    "n2rich": _species_mix(N2=0.95, O2=0.05),
    "vitiated": _species_mix(N2=0.72, O2=0.12, CO2=0.08, H2O=0.07, AR=0.01),
    "rich_ch4": _with_ch4(0.15),
}
TEMPS = (260.0, 300.0, 450.0, 700.0, 1100.0)


# ---------------------------------------------------------------------------
# Builders for what the stress families never draw
# ---------------------------------------------------------------------------


def _pb(g: FlowNetwork, nid: str, Pt: float, Tt: float = 300.0, Y: list | None = None) -> str:
    g.add_node(PressureBoundary(nid, Pt=Pt, Tt=Tt, Y=Y or _Y_AIR))
    return nid


def legacy_tee_ports(seed: int) -> tuple[FlowNetwork, list[str]]:
    """A Bassett tee with a chamber on every leg, so each leg delivers into
    an interior node whose gas the tee's rows read (the #1233/#1235 pattern;
    the stress family wires its branching tee straight to boundaries)."""
    r = random.Random(seed)
    g = FlowNetwork()
    tee_type = r.choice(["branching", "merging"])
    P = r.choice([1.5e5, 3e5])
    F_C = r.choice([2e-3, 1e-2])
    psi = r.choice([1.0, 2.0, 4.0])
    s = r.choice([0.01, 0.04])
    areas = {"pc": F_C, "ps": F_C, "pb": F_C / psi}
    for nid, a in areas.items():
        g.add_node(MomentumChamberNode(nid, area=a))
    if tee_type == "branching":
        _pb(g, "c", P, 400.0)
        _pb(g, "s", P * (1 - 0.3 * s))
        _pb(g, "b", P * (1 - s))
        links = (("lc", "c", "pc"), ("ls", "ps", "s"), ("lb", "pb", "b"))
    else:
        _pb(g, "s", P, 400.0)
        _pb(g, "b", P * (1 - 0.3 * s), 800.0)
        _pb(g, "c", P * (1 - s))
        links = (("ls", "s", "ps"), ("lb", "b", "pb"), ("lc", "pc", "c"))
    for eid, a, b in links:
        port = b if b.startswith("p") else a
        d = math.sqrt(4.0 * areas[port] / math.pi) * 0.6
        if r.random() < 0.5:
            g.add_element(LosslessConnectionElement(eid, a, b))
        else:
            g.add_element(OrificeElement(eid, a, b, Cd=0.8, diameter=d, correlation="fixed"))
    g.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="pc",
            straight_node="ps",
            branch_node="pb",
            theta=math.radians(r.choice([45.0, 90.0])),
            F_C=F_C,
            psi=psi,
            tee_type=tee_type,
        )
    )
    return g, ["legacy_tee_ports", tee_type]


def film(seed: int) -> tuple[FlowNetwork, list[str]]:
    """Film holes or a film slot from a coolant reservoir into a chamber the
    hot gas passes (test_hole_supply_drive's chamber network)."""
    r = random.Random(seed)
    g = FlowNetwork()
    kind = r.choice(["holes", "slot"])
    P = 5e5
    g.add_node(MassFlowBoundary("gas_in", m_dot=r.choice([0.5, 2.0]), Tt=1500.0, Y=_Y_AIR))
    _pb(g, "exit", P, 1500.0)
    _pb(g, "cool", P * (1.0 + r.choice([0.01, 0.03])), 600.0)
    g.add_node(MomentumChamberNode("ch", area=0.01, main_inlet="feed"))
    g.add_element(LosslessConnectionElement("feed", "gas_in", "ch"))
    g.add_element(LosslessConnectionElement("tail", "ch", "exit"))
    if kind == "holes":
        g.add_element(
            FilmCoolingElement(
                "h", "cool", "ch", hole_diameter=1e-3, holes_per_row=10,
                lateral_pitch=3e-3, wall_thickness=3e-3, angle_deg=30.0,
            )
        )  # fmt: skip
    else:
        g.add_element(
            SlotFilmElement("h", "cool", "ch", slot_height=1e-3, span=0.05, lip_thickness=1e-3)
        )
    return g, ["film", kind]


def machine(seed: int) -> tuple[FlowNetwork, list[str]]:
    """A compressor between two orifices (the turbine is drawn by
    fuel_controller)."""
    r = random.Random(seed)
    g = FlowNetwork()
    PR = r.choice([1.2, 1.6])
    _pb(g, "in", 1.0e5, r.choice([260.0, 300.0]))
    _pb(g, "out", 1.0e5 * PR * 0.95)
    g.add_node(PlenumNode("a"))
    g.add_node(PlenumNode("b"))
    g.add_element(OrificeElement("oi", "in", "a", Cd=0.8, diameter=0.03, correlation="fixed"))
    g.add_element(CompressorElement("comp", "a", "b", pressure_ratio=PR, eta_is=0.85))
    g.add_element(OrificeElement("oo", "b", "out", Cd=0.8, diameter=0.03, correlation="fixed"))
    return g, ["machine", "compressor"]


def port_discharge(seed: int) -> tuple[FlowNetwork, list[str]]:
    """A chamber discharging into a boundary (#637), behind a pipe."""
    r = random.Random(seed)
    g = FlowNetwork()
    _pb(g, "hi", 5e5, 300.0)
    g.add_node(MomentumChamberNode("c", area=5e-4))
    _pb(g, "lo", r.choice([4.0e5, 1.0e5]))
    D = math.sqrt(4.0 * 5e-4 / math.pi)
    g.add_element(ChannelElement("l", "hi", "c", length=2.0, diameter=D))
    g.add_element(PortDischargeElement("d", "c", "lo"))
    return g, ["port_discharge"]


def combustor(seed: int) -> tuple[FlowNetwork, list[str]]:
    """Air and a fuel stream into a combustor, complete or equilibrium; the
    fuel stream's flow is an unknown (an orifice), so the points take it to
    zero and reverse it: a fuel-free mix (#1241)."""
    r = random.Random(seed)
    g = FlowNetwork()
    method = r.choice(["equilibrium", "complete"])
    P = r.choice([2e5, 5e5])
    _pb(g, "air", P, r.choice([300.0, 650.0]))
    _pb(
        g,
        "fuel",
        P * 1.05,
        300.0,
        GASES["rich_ch4"] if r.random() < 0.3 else list(cb.species.pure_species("CH4")),
    )
    _pb(g, "out", P * 0.9)
    g.add_node(CombustorNode("burner", method=method))
    g.add_element(
        OrificeElement(
            "a_in", "air", "burner", Cd=0.8, diameter=0.02, correlation="fixed", in_pipe=False
        )
    )
    g.add_element(
        OrificeElement(
            "inj", "fuel", "burner", Cd=0.7, diameter=0.002, correlation="fixed", in_pipe=False
        )
    )
    g.add_element(
        OrificeElement(
            "exit", "burner", "out", Cd=0.8, diameter=0.025, correlation="fixed", in_pipe=False
        )
    )
    return g, ["combustor", f"method:{method}"]


#: The sweep's own builders, beside the stress FAMILIES.
EXTRA: dict[str, Callable[[int], tuple[FlowNetwork, list[str]]]] = {
    "legacy_tee_ports": legacy_tee_ports,
    "film": film,
    "machine": machine,
    "port_discharge": port_discharge,
    "combustor": combustor,
}
BUILDERS = {**FAMILIES, **EXTRA}


# ---------------------------------------------------------------------------
# Cases: a builder's draw, its gases and side streams
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Case:
    """One network: a builder's seed, a variant (0: as drawn; otherwise
    boundary gases, temperatures and combustor methods redrawn) and how many
    side streams to inject."""

    family: str
    seed: int
    variant: int = 0
    side_streams: int = 0

    def label(self) -> str:
        return f"{self.family}:{self.seed}/v{self.variant}/s{self.side_streams}"

    def key(self) -> int:
        """A seed from the label (str hashes differ per process)."""
        return zlib.crc32(self.label().encode())


def _interior(net: FlowNetwork) -> list[str]:
    return [
        nid
        for nid, n in net.nodes.items()
        if isinstance(n, PlenumNode | MomentumChamberNode | CombustorNode)
    ]


def _port_nodes(net: FlowNetwork) -> list[str]:
    """Interior nodes an element with more than two nodes (a junction)
    touches."""
    interior = set(_interior(net))
    out: list[str] = []
    for e in net.elements.values():
        nodes = e.all_source_nodes() + e.all_sink_nodes()
        if len(nodes) > 2:
            out += [n for n in nodes if n in interior and n not in out]
    return out


def build(case: Case) -> tuple[FlowNetwork, list[str]]:
    """The case's network, variant applied, unresolved."""
    net, tags = BUILDERS[case.family](case.seed)
    tags = list(tags)
    if case.variant == 0 and case.side_streams == 0:
        return net, tags
    r = random.Random(case.key())
    if case.variant:
        for n in net.nodes.values():
            if isinstance(n, FuelBoundary):
                continue  # its Tt_target was drawn for its own gas
            if isinstance(n, PressureBoundary | MassFlowBoundary):
                n.Tt = r.choice(TEMPS)
                n.Y = list(GASES[r.choice(list(GASES))])
            elif isinstance(n, CombustorNode):
                n.method = r.choice(["complete", "equilibrium"])
        tags.append("gases")
    # A MultiPortChamberBase port takes exactly one outside element (the
    # network refuses a second), so side streams go to the other nodes.
    mpcb = {
        n
        for e in net.elements.values()
        if isinstance(e, cbn.MultiPortChamberBase)
        for n in e.all_source_nodes() + e.all_sink_nodes()
    }
    ports = [n for n in _port_nodes(net) if n not in mpcb]
    targets = ports + [n for n in _interior(net) if n not in ports and n not in mpcb]
    p_hi = max(
        (float(n.Pt) for n in net.nodes.values() if isinstance(n, PressureBoundary)),
        default=2e5,
    )
    for k, nid in enumerate(targets[: case.side_streams]):
        src = f"side{k}"
        _pb(net, src, p_hi * 1.05, r.choice(TEMPS), list(GASES[r.choice(list(GASES))]))
        net.add_element(
            OrificeElement(
                f"ss{k}", src, nid, Cd=0.7, diameter=0.003, correlation="fixed", in_pipe=False
            )
        )
        tags.append(f"side_stream:{nid}")
    return net, tags


def resolved_solver(case: Case) -> NetworkSolver:
    """A fresh solver on a freshly built, resolved network."""
    net, _ = build(case)
    s = NetworkSolver(net)
    net.resolve_all_topology()
    s._build_x0()
    return s


# ---------------------------------------------------------------------------
# Points
# ---------------------------------------------------------------------------


@dataclass
class Point:
    kind: str
    desc: str
    x: np.ndarray
    #: A flow unknown set to exactly zero: compared one-sided (``_confirm``).
    zero_col: int | None = None


def _is_flow(name: str) -> bool:
    return name.rsplit(".", 1)[-1].startswith("m_dot")


def _own_flow_cols(s: NetworkSolver) -> dict[str, list[int]]:
    """Each element's own signed flow unknowns (``.m_dot*``; a FuelBoundary's
    ``ln_m_dot`` is a logarithm, not a signed flow)."""
    names = s.unknown_names
    out: dict[str, list[int]] = {}
    for eid in s.network.elements:
        cols = [j for j in s._unknown_indices.get(eid, []) if _is_flow(names[j])]
        if cols:
            out[eid] = cols
    return out


def _flow_cols(s: NetworkSolver) -> dict[str, list[int]]:
    """Each element's legs as flow unknowns: its own, or for a junction with
    none of its own (the MultiPortChamberBase family: its port flows are its
    arms'), the first arm's flow at each port."""
    own = _own_flow_cols(s)
    out = dict(own)
    net = s.network
    for eid, e in net.elements.items():
        nodes = e.all_source_nodes() + e.all_sink_nodes()
        if eid in own or len(nodes) <= 2:
            continue
        arms: list[int] = []
        for nid in nodes:
            for a in net.get_upstream_elements(nid) + net.get_downstream_elements(nid):
                if a.id != eid and a.id in own:
                    arms.append(own[a.id][0])
                    break
        if arms:
            out[eid] = arms
    return out


def _legs(e: Any, x: np.ndarray, cols: list[int]) -> list[float]:
    """An element's signed leg flows: its flow unknowns, and a legacy tee's
    straight run (m_com - m_branch)."""
    legs = [float(x[j]) for j in cols]
    if isinstance(e, TeeJunctionElement) and len(cols) == 2:
        legs.append(legs[0] - legs[1])
    return legs


def _node_flow_cols(s: NetworkSolver) -> dict[str, list[int]]:
    """Per interior node, the own flow unknowns of the elements touching it."""
    fc = _own_flow_cols(s)
    out: dict[str, list[int]] = {}
    for eid, e in s.network.elements.items():
        if eid not in fc:
            continue
        for nid in e.all_source_nodes() + e.all_sink_nodes():
            out.setdefault(nid, [])
            out[nid] += [j for j in fc[eid] if j not in out[nid]]
    interior = set(_interior(s.network))
    return {n: c for n, c in out.items() if n in interior}


def _start(s: NetworkSolver) -> np.ndarray:
    """The start vector with resting flows moved off zero (``exact_zero``):
    the reference every point perturbs, and the flow columns' scale."""
    x0 = np.asarray(s._build_x0(), dtype=float)
    flows = sorted({j for c in _flow_cols(s).values() for j in c})
    nz = [abs(x0[j]) for j in flows if x0[j] != 0.0]
    ref = float(np.median(nz)) if nz else 0.1
    for j in flows:
        if x0[j] == 0.0:
            x0[j] = 0.1 * ref
    return x0


def points(s: NetworkSolver, seed: int, max_points: int | None = MAX_POINTS) -> list[Point]:
    """The sweep's points for one resolved solver (module docstring), the
    optional ones sampled down to ``max_points`` with ``seed``."""
    r = random.Random(seed)
    base = _start(s)
    names = s.unknown_names
    net = s.network
    fc = _flow_cols(s)
    flows = sorted({j for c in fc.values() for j in c})
    pts = [Point("x0", "start vector", base.copy())]
    for k in range(2):
        x = base.copy()
        for j in range(x.size):
            if is_pressure(names[j]):
                x[j] *= 1.0 + 0.01 * r.uniform(-1.0, 1.0)
            elif j in flows:
                x[j] *= math.exp(0.5 * r.uniform(-1.0, 1.0))
        pts.append(Point("perturb", f"perturbation {k}", x))
    mandatory: list[Point] = []
    optional: list[Point] = []
    for eid, cols in fc.items():
        n = len(cols)
        pats = list(range(1, 2**n))
        if len(pats) > MAX_PATTERNS:
            pats = pats[:n] + r.sample(pats[n:], MAX_PATTERNS - n)
        first = True
        for p in pats:
            for ratio in RATIOS if n > 1 else (1.0,):
                x = base.copy()
                for q, j in enumerate(cols):
                    sgn = -1.0 if (p >> q) & 1 else 1.0
                    x[j] = sgn * base[j] * (ratio if q == n - 1 else 1.0)
                pt = Point("signs", f"{eid} pattern {p:0{n}b} ratio {ratio:g}", x)
                (mandatory if first else optional).append(pt)
                first = False
        for mag in ELEMENT_NEAR_ZERO:
            for sgn in (1.0, -1.0):
                x = base.copy()
                for q, j in enumerate(cols):
                    # Distinct magnitudes keep a derived leg (a tee's
                    # straight run) off zero.
                    x[j] = sgn * mag * (1.0 + 0.37 * q)
                optional.append(Point("element_near_zero", f"{eid} flows {sgn * mag:g}", x))
    node_cols = _node_flow_cols(s)
    mixing = {n for n, c in node_cols.items() if len(c) >= 2}
    interior = set(_interior(net))
    for eid, cols in _own_flow_cols(s).items():
        ends = net.elements[eid].all_source_nodes() + net.elements[eid].all_sink_nodes()
        if len(ends) != 2 or sum(n in interior for n in ends) != 1:
            continue  # zero_between (or a junction's own legs)
        for j in cols:
            x = base.copy()
            x[j] = 0.0
            pt = Point("exact_zero", f"{names[j]} = 0", x, zero_col=j)
            (mandatory if any(n in mixing for n in ends) else optional).append(pt)
    for nid, cols in node_cols.items():
        if len(cols) > 9:
            continue
        for mag in (NEAR_ZERO_OUT_BAND, NEAR_ZERO_IN_BAND):
            x = base.copy()
            for q, j in enumerate(cols):
                x[j] = r.choice((1.0, -1.0)) * mag * (1.0 + 0.37 * q)
            pt = Point("node_near_zero", f"{nid} flows +-{mag:g}", x)
            (mandatory if mag == NEAR_ZERO_OUT_BAND else optional).append(pt)
    if max_points is not None and len(pts) + len(mandatory) + len(optional) > max_points:
        # The mandatory points first; below their count (a quick tier) a
        # seeded sample of them.
        room = max(max_points - len(pts), 0)
        if len(mandatory) > room:
            mandatory, optional = r.sample(mandatory, room), []
        else:
            optional = r.sample(optional, min(room - len(mandatory), len(optional)))
    return pts + mandatory + optional


# ---------------------------------------------------------------------------
# Finite differences
# ---------------------------------------------------------------------------

#: A non-pressure column's confirming ladder halves from LADDER_TOP of its
#: scale (a quarter of |x| at most, central) to LADDER_BOTTOM of it: a noisy
#: residual (an equilibrium solve) resolves only at the large steps, one that
#: curves on a small scale only at the small ones, and ``_pick`` finds the
#: window either way. A near-zero flow (``_flow_step`` at its |x| / 4 cap)
#: takes a one-sided ladder away from zero instead, its stencil within |x| of
#: the point: it can curve on |x| itself (a junction's coefficients on its
#: port-flow ratios, a duct's near-zero regime), and a longer step measures a
#: secant that can agree with itself.
LADDER_TOP, LADDER_BOTTOM = 1e-3, 1e-9
MAX_LEVELS = 32
#: A near-zero flow's one-sided ladder: from |x| / 2 down to this times |x|.
NEAR_ZERO_BOTTOM = 1e-4
#: An exactly-zero flow's one-sided ladders: from ZERO_TOP to ZERO_BOTTOM of
#: its column's scale.
ZERO_TOP, ZERO_BOTTOM = 1e-6, 1e-12
#: The FD's estimate must be within this fraction of its value (and below a
#: tenth of the miss) for a miss to count: a window outside the asymptotic
#: range can agree with itself to a fraction of a value far from J.
FD_RESOLVED = 1e-2
#: One residual evaluation's round-off relative to |F_i|, with margin: no
#: difference resolves below ROUNDOFF |F_i| / h (a row of fixed boundary
#: pressures reads exact zeros at every level).
ROUNDOFF = 1e-13


def column_scale(s: NetworkSolver, x: np.ndarray, x0: np.ndarray) -> np.ndarray:
    """|x_j|, a flow column's at least its start value's: at a near-zero
    flow J_ij |x_j| would hide any error in the column, and the start value
    is the flow change a Newton step makes there."""
    sc = np.abs(x).copy()
    for j, n in enumerate(s.unknown_names):
        if _is_flow(n):
            sc[j] = max(sc[j], abs(x0[j]))
    return np.maximum(sc, 1e-12)


def _flow_step(xj: float, scj: float) -> float:
    """A flow column's screen step: 1e-6 of its scale, at most a quarter of
    its distance to zero."""
    return min(1e-6 * scj, 0.25 * abs(xj)) if xj != 0.0 else 1e-6 * scj


def _steps(s: NetworkSolver, x: np.ndarray, x0: np.ndarray) -> tuple[np.ndarray, float]:
    """Per-column ``rel`` and the ``floor`` for fd_jacobian: 1e-6 |x| as
    usual (1e-9 floor), a flow column ``_flow_step``."""
    sc = column_scale(s, x, x0)
    h = np.maximum(1e-6 * np.abs(x), 1e-9)
    for j, n in enumerate(s.unknown_names):
        if _is_flow(n):
            h[j] = _flow_step(x[j], sc[j])
    floor = float(min(h.min(), 1e-9))
    rel = np.where(x != 0.0, h / np.where(x != 0.0, np.abs(x), 1.0), 1.0)
    return rel, floor


def _signature(s: NetworkSolver, x: np.ndarray) -> tuple[bool, ...]:
    """The sign of every leg (``_legs``): where one flips, a stream changes
    direction, a mixer its arriving set, a junction its pattern."""
    fc = _flow_cols(s)
    return tuple(
        leg > 0.0 for eid, e in s.network.elements.items() for leg in _legs(e, x, fc.get(eid, []))
    )


def _ladder(
    s: NetworkSolver, x: np.ndarray, j: int, h0: float, levels: int, side: float = 0.0
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Column j on a ladder halving from ``h0``, each level Richardson-
    extrapolated: (R per level, base steps, allowed). ``side`` 0: central;
    +-1: one-sided on that side (stencil x, x + side 2h). A level is allowed
    only while its whole stencil keeps every leg's sign (``_signature``):
    one across a derived leg's zero differences across a kink, and can agree
    with its neighbours there."""
    x_ref = x.copy()
    if x_ref[j] == 0.0:  # an exact zero: its own leg takes the side's sign
        x_ref[j] = (side or 1.0) * 1e-300
    sig = _signature(s, x_ref)
    F0 = _evaluate(s, x)[0] if side else None
    hs, D, ok = [], [], []
    for k in range(levels):
        h = h0 / 2.0**k
        if side:
            xp, far = x.copy(), x.copy()
            xp[j] += side * h
            far[j] += side * 2.0 * h
            D.append((_evaluate(s, xp)[0] - F0) / (xp[j] - x[j]))
            ok.append(_signature(s, far) == sig)
        else:
            xp, xm, fp, fm = x.copy(), x.copy(), x.copy(), x.copy()
            xp[j] += h
            xm[j] -= h
            fp[j] += 2.0 * h
            fm[j] -= 2.0 * h
            D.append((_evaluate(s, xp)[0] - _evaluate(s, xm)[0]) / (xp[j] - xm[j]))
            ok.append(_signature(s, fp) == sig and _signature(s, fm) == sig)
        hs.append(h)
    Da = np.array(D)
    # Central: (4 D(h) - D(2h)) / 3 cancels the h^2 term; one-sided:
    # 2 D(h) - D(2h) the h term. A level's R uses D(h) and D(2h): both allowed.
    R = (4.0 * Da[1:] - Da[:-1]) / 3.0 if not side else 2.0 * Da[1:] - Da[:-1]
    okA = np.array(ok)
    return R, np.array(hs[1:]), okA[1:] & okA[:-1]


def _raw_est(fd: FDJacobian, x: np.ndarray, F0: np.ndarray) -> np.ndarray:
    """fd_jacobian's estimates back in J's units, floored at round-off."""
    C = np.asarray(fd)
    xs = np.maximum(np.abs(x), 1e-12)
    row = np.maximum(np.max(np.abs(C) * xs[None, :], axis=1), 1e-300)
    E = np.asarray(fd.est, dtype=float) * row[:, None] / xs[None, :]
    steps = np.maximum(np.abs(np.asarray(fd.steps, dtype=float)), 1e-300)
    return np.maximum(np.nan_to_num(E, nan=np.inf), ROUNDOFF * F0[:, None] / steps)


@dataclass
class Confirmed:
    """``_confirm``'s result: per checked row (verdict, column, scaled error,
    scaled FD estimate), and the FD as an FDJacobian whose ``why`` names
    each entry's step."""

    rows: dict[int, tuple[str, int, float, float]]
    fd: FDJacobian
    scale: np.ndarray


def _confirm(
    case: Case,
    x: np.ndarray,
    x0: np.ndarray,
    rows: list[int],
    J: np.ndarray,
    zero_col: int | None = None,
) -> Confirmed:
    """Re-difference ``rows`` (module docstring, "Comparison") on a fresh
    resolved solver. ``zero_col``: that flow is exactly zero; each row takes
    the one-sided difference (either side) it agrees with."""
    s = resolved_solver(case)
    names = s.unknown_names
    rel, floor = _steps(s, x, x0)
    F0 = np.abs(_evaluate(s, x)[0])
    fd = fd_jacobian(s, x, rel, floor)
    C = np.asarray(fd).copy()
    E = _raw_est(fd, x, F0)
    why = np.array(fd.bounds, dtype=object)
    steps = np.asarray(fd.steps, dtype=float).copy()
    sc = column_scale(s, x, x0)
    p_cols = [j for j, n in enumerate(names) if is_pressure(n)]
    if any(E[i, j] > FD_RESOLVED * abs(C[i, j]) for i in rows for j in p_cols):
        fd_b = fd_jacobian(resolved_solver(case), x, rel * 100.0, floor)
        E_b = _raw_est(fd_b, x, F0)
        for i in rows:
            for j in p_cols:
                if E_b[i, j] < E[i, j]:
                    C[i, j], E[i, j] = fd_b[i, j], E_b[i, j]
                    why[i, j], steps[i, j] = f"x100 {fd_b.bounds[i, j]}", fd_b.steps[i, j]
    for j, n in enumerate(names):
        if is_pressure(n):
            continue
        base = max(abs(x[j]) * rel[j], floor)
        if j == zero_col:
            # No scale but the column's: from ZERO_TOP of it down, so the
            # levels sit where the stream is a trace (a large one-sided step
            # from zero measures a secant across whatever the flow does
            # there, and can agree with itself).
            h0 = ZERO_TOP * sc[j]
            n_lv = int(math.ceil(math.log2(ZERO_TOP / ZERO_BOTTOM))) + 1
            ladders = {sd: _ladder(s, x, j, h0, n_lv, sd) for sd in (1.0, -1.0)}
        elif _is_flow(n) and 4.0 * base > 0.5 * abs(x[j]):
            # Near zero: one-sided away from zero, the stencil within |x| of
            # the point (the residual can curve on |x| itself).
            sd = math.copysign(1.0, x[j])
            h0 = 0.5 * abs(x[j])
            n_lv = int(math.ceil(math.log2(0.5 / NEAR_ZERO_BOTTOM))) + 1
            ladders = {sd: _ladder(s, x, j, h0, n_lv, sd)}
        else:
            h0 = min(LADDER_TOP * sc[j], 0.25 * abs(x[j])) if _is_flow(n) else LADDER_TOP * sc[j]
            n_lv = int(math.ceil(math.log2(max(h0 / (LADDER_BOTTOM * sc[j]), 4.0)))) + 1
            ladders = {0.0: _ladder(s, x, j, h0, min(n_lv, MAX_LEVELS))}
        for i in rows:
            best = None
            for sd, (R, hs, ok) in ladders.items():
                if not ok.any():
                    continue
                v, k, e, _ = _pick(R[:, i], hs, ok)
                e = max(float(np.nan_to_num(e, nan=np.inf)), ROUNDOFF * F0[i] / hs[k])
                if best is None or abs(J[i, j] - v) < abs(J[i, j] - best[0]):
                    how = "central" if sd == 0.0 else f"one-sided {sd:+.0f}"
                    best = (v, e, f"{how} ladder", hs[k])
            if best is None:  # every level crosses a leg's zero
                C[i, j], E[i, j], why[i, j], steps[i, j] = J[i, j], np.inf, "no level", np.nan
            else:
                C[i, j], E[i, j], why[i, j], steps[i, j] = best
    out: dict[int, tuple[str, int, float, float]] = {}
    for i in rows:
        row = max(float(np.max(np.abs(C[i] * sc))), 1e-300)
        d = np.abs((J[i] - C[i]) * sc) / row
        j = int(np.argmax(d))
        est = float(E[i, j]) * sc[j] / row
        if d[j] <= TOL:
            verdict = "ok"
        elif est <= 0.1 * d[j] and E[i, j] <= FD_RESOLVED * abs(C[i, j]):
            verdict = "mismatch"
        else:
            verdict = "unresolved"
        out[i] = (verdict, j, float(d[j]), est)
    # E in scaled units, as fd_jacobian reports ``est``.
    rows_sc = np.maximum(np.max(np.abs(C) * sc[None, :], axis=1), 1e-300)
    view = C.view(FDJacobian)
    view.steps, view.bounds, view.names = steps, why, list(names)
    view.est = E * sc[None, :] / rows_sc[:, None]
    return Confirmed(out, view, sc)


# ---------------------------------------------------------------------------
# Coverage and findings
# ---------------------------------------------------------------------------


@dataclass
class Coverage:
    """Compared points that reached a class: in all, with a leg reversed
    against the start, with differing compositions and with differing
    temperatures (an element: across the nodes it touches; a node: over the
    streams arriving at it)."""

    reached: int = 0
    reversed: int = 0
    composition: int = 0
    temperature: int = 0

    def add(self, rev: bool, comp: bool, temp: bool) -> None:
        self.reached += 1
        self.reversed += rev
        self.composition += comp
        self.temperature += temp


@dataclass
class Finding:
    case: str
    point: str
    kind: str
    verdict: str
    element: str
    row_kind: str
    column: str
    error: float
    fd_estimate: float
    #: ``worst``'s message for the point's mismatching rows.
    message: str = ""
    #: The network's classes, and "equilibrium" when a combustor runs it.
    tags: list[str] = field(default_factory=list)


@dataclass
class Report:
    coverage: dict[str, Coverage] = field(default_factory=dict)
    findings: list[Finding] = field(default_factory=list)
    excluded: dict[str, int] = field(default_factory=dict)
    compared: int = 0
    build_refused: list[str] = field(default_factory=list)
    nonfinite: dict[str, int] = field(default_factory=dict)
    #: Wall time per family [s].
    seconds: dict[str, float] = field(default_factory=dict)

    def cov(self, cls: str) -> Coverage:
        return self.coverage.setdefault(cls, Coverage())

    def exclude(self, rule: str) -> None:
        self.excluded[rule] = self.excluded.get(rule, 0) + 1

    def mismatches(self) -> list[Finding]:
        return [f for f in self.findings if f.verdict == "mismatch"]


def _node_state(s: NetworkSolver, nid: str) -> tuple[float, list[float]] | None:
    n = s.network.nodes[nid]
    if isinstance(n, PressureBoundary | MassFlowBoundary):
        Y = n.Y if n.Y is not None else s._default_Y
        return float(n.Tt), list(Y)
    st = s._derived_states.get(nid)
    return None if st is None else (float(st[0]), list(st[1]))


def _differ(states: list[tuple[float, list[float]]]) -> tuple[bool, bool]:
    """(composition differs, temperature differs) over the states."""
    if len(states) < 2:
        return False, False
    Ys = np.array([y for _, y in states])
    Ts = np.array([t for t, _ in states])
    return bool(np.max(Ys.max(0) - Ys.min(0)) > 1e-6), bool(Ts.max() - Ts.min() > 1.0)


def _record_coverage(rep: Report, s: NetworkSolver, x: np.ndarray, x0: np.ndarray) -> None:
    s._propagate_states(x)
    net = s.network
    fc = _flow_cols(s)
    rev_el: dict[str, bool] = {}
    for eid, e in net.elements.items():
        cols = fc.get(eid, [])
        legs, legs0 = _legs(e, x, cols), _legs(e, x0, cols)
        rev = any(a * b < 0.0 for a, b in zip(legs, legs0, strict=True)) or any(
            x[j] < 0.0 for j in cols
        )
        rev_el[eid] = rev
        ends = e.all_source_nodes() + e.all_sink_nodes()
        comp, temp = _differ([st for n in ends if (st := _node_state(s, n)) is not None])
        rep.cov(type(e).__name__).add(rev, comp, temp)
    for nid, n in net.nodes.items():
        adj = net.get_upstream_elements(nid) + net.get_downstream_elements(nid)
        rev = any(rev_el.get(e.id, False) for e in adj)
        try:
            inflows = s._node_inflows(nid, x)
        except Exception:  # noqa: BLE001 -- coverage bookkeeping only
            inflows = []
        comp, temp = _differ([(t, y) for m, t, y in inflows if m > 0.0])
        rep.cov(type(n).__name__).add(rev, comp, temp)
    for w in net.walls.values():
        rep.cov(type(w).__name__).add(False, False, True)


def _in_band(s: NetworkSolver, x: np.ndarray) -> bool:
    """Whether an interior node's arriving flow (all of its flow when none
    arrives, as the mixer weighs it) is inside the mixer band."""
    for nid in _interior(s.network):
        try:
            flows = [m for m, _, _ in s._node_inflows(nid, x)]
        except Exception:  # noqa: BLE001
            continue
        m = sum(f for f in flows if f > 0.0) or sum(abs(f) for f in flows)
        if flows and abs(m) <= MIXER_NO_FLOW_BAND:
            return True
    return False


def _refused(rep: Report, case: Case, exc: Exception, where: str) -> None:
    rep.exclude("nonfinite")
    key = f"{case.family} ({where}): {type(exc).__name__}: {str(exc)[:80]}"
    rep.nonfinite[key] = rep.nonfinite.get(key, 0) + 1


def _findings(
    rep: Report, case: Case, pt: Point, s: NetworkSolver, J: np.ndarray, c: Confirmed
) -> None:
    """The point's misses, the worst per row owner and verdict; mismatches
    carry ``worst``'s message (row, column, scaled error, step)."""
    kinds = s._row_kinds(pt.x)
    labels = [f"{o} [{k}]" for o, k in kinds]
    names = s.unknown_names
    bad = [i for i, v in c.rows.items() if v[0] == "mismatch"]
    message = ""
    if bad:
        Jm = np.asarray(c.fd).copy()
        Jm[bad] = J[bad]
        message = worst(Jm, c.fd, c.scale, labels)[1]
    net = s.network
    tags = sorted({type(e).__name__ for e in net.elements.values()})
    if any(getattr(n, "method", None) == "equilibrium" for n in net.nodes.values()):
        tags.append("equilibrium")
    by_owner: dict[tuple[str, str], Finding] = {}
    for i, (verdict, j, e, est) in c.rows.items():
        if verdict == "ok":
            continue
        if verdict == "unresolved":
            rep.exclude("unresolved")
        owner, rk = kinds[i]
        f = Finding(case.label(), pt.desc, pt.kind, verdict, owner, rk, names[j], e, est, tags=tags)
        if verdict == "mismatch":
            f.message = message
        key = (owner, verdict)
        if key not in by_owner or by_owner[key].error < e:
            by_owner[key] = f
    rep.findings += list(by_owner.values())


def evaluate(
    case: Case,
    pt: Point,
    x0: np.ndarray,
    rep: Report,
    *,
    coverage_only: bool = False,
    perturb_J: Callable[[NetworkSolver, np.ndarray], None] | None = None,
) -> None:
    """One point: the exclusion rules, coverage, the FD screen and the
    confirmation. ``perturb_J`` edits the analytic J in place before the
    comparison: the tests' check that the sweep flags a wrong entry."""
    with _quiet():
        try:
            s = resolved_solver(case)
            F, J = s._residuals_and_jacobian(pt.x)
            J, F = _dense(J), np.asarray(F, dtype=float)
        except Exception as exc:  # noqa: BLE001 -- counted, a finding of its own kind
            _refused(rep, case, exc, "point")
            return
        if not (np.all(np.isfinite(F)) and np.all(np.isfinite(J))):
            _refused(rep, case, FloatingPointError("non-finite F or J"), "point")
            return
        fc = _flow_cols(s)
        x_chk = pt.x.copy()
        if pt.zero_col is not None:
            x_chk[pt.zero_col] = 1e-300
        if any(
            leg == 0.0
            for eid, e in s.network.elements.items()
            for leg in _legs(e, x_chk, fc.get(eid, []))
        ):
            rep.exclude("exact_zero")
            return
        if _in_band(s, pt.x):
            rep.exclude("in_band")
            return
        _record_coverage(rep, s, pt.x, x0)
        rep.compared += 1
        if coverage_only:
            return
        if perturb_J is not None:
            perturb_J(s, J)
        try:
            if pt.zero_col is not None:
                rows = list(range(F.size))
            else:
                rel, floor = _steps(s, pt.x, x0)
                fd = fd_jacobian(resolved_solver(case), pt.x, rel, floor, adaptive=False)
                sc = column_scale(s, pt.x, x0)
                err = scaled_error(J * sc[None, :], np.asarray(fd) * sc[None, :], np.ones_like(sc))
                rows = [int(i) for i in np.flatnonzero(err > TOL)]
            if not rows:
                return
            c = _confirm(case, pt.x, x0, rows, J, pt.zero_col)
        except Exception as exc:  # noqa: BLE001 -- a stencil the model refuses
            _refused(rep, case, exc, "stencil")
            return
        _findings(rep, case, pt, s, J, c)


def _dense(J: Any) -> np.ndarray:
    return J.toarray() if hasattr(J, "toarray") else np.asarray(J)


@contextlib.contextmanager
def _quiet() -> Iterator[None]:
    with warnings.catch_warnings(), np.errstate(all="ignore"):
        warnings.simplefilter("ignore")
        yield


def case_points(case: Case, max_points: int | None = MAX_POINTS) -> tuple[np.ndarray, list[Point]]:
    """The case's start vector (``_start``) and its points."""
    with _quiet():
        s = resolved_solver(case)
        return _start(s), points(s, case.key(), max_points)


def point(case: Case, desc: str) -> tuple[np.ndarray, Point]:
    """One named point of a case (its ``Point.desc``), for a reproducer."""
    x0, pts = case_points(case, None)
    return x0, next(p for p in pts if p.desc == desc)


def run_case(
    case: Case, rep: Report, *, max_points: int | None = MAX_POINTS, coverage_only: bool = False
) -> None:
    try:
        x0, pts = case_points(case, max_points)
    except Exception as exc:  # noqa: BLE001 -- a variant the network refuses
        rep.build_refused.append(f"{case.label()}: {type(exc).__name__}: {str(exc)[:120]}")
        return
    t0 = time.perf_counter()
    for pt in pts:
        evaluate(case, pt, x0, rep, coverage_only=coverage_only)
    rep.seconds[case.family] = rep.seconds.get(case.family, 0.0) + time.perf_counter() - t0


# ---------------------------------------------------------------------------
# Case sets
# ---------------------------------------------------------------------------

#: Networks above this many unknowns are left out: FD at many points of a
#: large network costs minutes, and the small draws reach the same classes.
MAX_UNKNOWNS = 24
#: (variant, side streams) per seed: as drawn, gases only, gases and one
#: side stream, gases and two.
VARIANTS = ((0, 0), (1, 0), (2, 1), (3, 2))
#: Seeds every builder is swept at, and the extra draws that reach a class
#: or pattern those miss: random 8 a WallNode, 9 and 20 a CombustorNode
#: (20 with a WallNode); fuel_controller 7 a ThermalWall, 8 a TurbineElement.
SEEDS = range(3)
#: The slow test tier: the first seed (and EXTRA_SEEDS) at MAX_POINTS_CI
#: points a case; the script runs SEEDS at MAX_POINTS.
EXTRA_SEEDS: dict[str, tuple[int, ...]] = {"random": (8, 9, 20), "fuel_controller": (7, 8)}
#: The slow test tier (CI time): the first seed and the extra draws but
#: fuel_controller 7 (a wall-coupled equilibrium burner, 22 ms a residual),
#: gases redrawn without and with two side streams, MAX_POINTS_CI points a
#: case. The script runs SEEDS, EXTRA_SEEDS and VARIANTS at MAX_POINTS.
SEEDS_CI = (0,)
EXTRA_SEEDS_CI: dict[str, tuple[int, ...]] = {
    "random": (2, 8, 9),  # 2: an AreaChangeElement, 8: a WallNode, 9: a combustor
    "fuel_controller": (8,),
    "film": (1,),  # film holes (seed 0 is a slot)
    "tee3": (1,),  # a MultiPortChamberElement (seed 0 is a ConstantKTee)
    "series": (1,),  # a VortexElement
}
VARIANTS_CI = ((1, 0), (3, 2))
MAX_POINTS_CI = 8


def full_cases(
    seeds: tuple[int, ...] | range = SEEDS,
    extra: dict[str, tuple[int, ...]] = EXTRA_SEEDS,
    variants: tuple[tuple[int, int], ...] = VARIANTS,
    families: list[str] | None = None,
    *,
    small_only: bool = True,
) -> list[Case]:
    """Every builder (or ``families``) at ``seeds`` and ``extra``, each in
    ``variants``; small networks only (MAX_UNKNOWNS) unless ``small_only``
    is False (a test that skips the large ones itself)."""
    out = []
    for fam in families or list(BUILDERS):
        for seed in dict.fromkeys((*seeds, *extra.get(fam, ()))):
            for v, ss in variants:
                c = Case(fam, seed, v, ss)
                if not small_only or _small(c):
                    out.append(c)
    return out


def ci_cases(families: list[str] | None = None, *, small_only: bool = True) -> list[Case]:
    """The slow test tier's cases."""
    return full_cases(SEEDS_CI, EXTRA_SEEDS_CI, VARIANTS_CI, families, small_only=small_only)


def _small(case: Case) -> bool:
    try:
        with _quiet():
            return len(resolved_solver(case).unknown_names) <= MAX_UNKNOWNS
    except Exception:  # noqa: BLE001
        return True  # run_case records the refusal


def sweep(
    cases: list[Case], *, max_points: int | None = MAX_POINTS, coverage_only: bool = False
) -> Report:
    rep = Report()
    for c in cases:
        run_case(c, rep, max_points=max_points, coverage_only=coverage_only)
    return rep


# ---------------------------------------------------------------------------
# Open mismatches on main
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Open:
    """A mismatch pattern open on main: the issue (a placeholder until
    filed), what it is, which findings it covers, and one reproducer (a case
    and a point desc) that the fast tier runs as a strict xfail."""

    issue: str
    what: str
    covers: Callable[[Finding], bool]
    case: Case
    point: str
    #: The reproducer takes over a second (a Fanno network): the slow tier.
    slow: bool = False


def _fam(f: Finding) -> str:
    return f.case.split(":", 1)[0]


def _cls(f: Finding) -> str:
    return f.element.split(" ", 1)[0]


#: Each pattern is measured in the PR that added the sweep; the sweep fails
#: on any mismatch none of them covers.
OPEN: tuple[Open, ...] = (
    Open(
        "#1241",
        "equilibrium combustor, a fuel stream at exactly zero flow (a fuel-free mix): "
        "6.8e-4 .. 8e-1 on the fuel column",
        lambda f: (
            "equilibrium" in f.tags and f.kind == "exact_zero" and f.column.startswith("inj.")
        ),
        Case("combustor", 1, 0, 0),
        "inj.m_dot = 0",
    ),
    Open(
        "#ISSUE-EQ-DEGENERATE-MIX",
        "equilibrium combustor at a near-zero oxidizer flow or with a stream reversed "
        "(complete combustion is clean at the same points): 1.3e-5 .. 7.5e-1",
        lambda f: "equilibrium" in f.tags and f.kind != "exact_zero",
        Case("combustor", 1, 0, 0),
        "a_in flows 3e-05",
    ),
    Open(
        "#ISSUE-EXACT-ZERO",
        "J at an exactly-zero flow with one boundary end matches neither one-sided "
        "difference (tee ports, ConstantKTee, ejector outlet chamber, crossflow segments, "
        "a wall-coupled channel): 1e-6 .. 5e-2",
        lambda f: f.kind == "exact_zero",
        Case("legacy_tee_ports", 2, 2, 1),
        "lc.m_dot = 0",
    ),
    Open(
        "#ISSUE-EJECTOR-REVERSED",
        "EjectorElement rows (and its chambers') with a leg reversed, one gas included: "
        "1e-6 .. 1.8e-1",
        lambda f: _fam(f) == "ejector",
        Case("ejector", 0, 0, 0),
        "cs pattern 1 ratio 1",
    ),
    Open(
        "#ISSUE-MPCE-NEAR-ZERO-ARM",
        "MultiPortChamberElement momentum row's port-pressure column with one arm near "
        "zero flow: 1.5e-1",
        lambda f: _cls(f) == "MultiPortChamberElement",
        Case("tee3", 2, 0, 0),
        "arm0 flows 3e-05",
        slow=True,
    ),
    Open(
        "#ISSUE-IMPINGEMENT-PLATE-GASES",
        "ImpingementPlateElement flow row's chamber-pressure column once the boundary "
        "gases differ, at the start vector already: 7.7e-3",
        lambda f: _cls(f) == "ImpingementPlateElement",
        Case("impingement_bypass", 0, 1, 0),
        "start vector",
        slow=True,
    ),
    Open(
        "#ISSUE-WALL-CHANNEL-PT",
        "a wall-coupled high-Mach channel's flow row on its outlet chamber's Pt, one gas "
        "included (hotcold families): 1e-6 .. 4e-6",
        lambda f: _fam(f) in ("hotcold", "hotcold_mild") and f.kind != "exact_zero",
        Case("hotcold_mild", 0, 1, 0),
        "start vector",
        slow=True,
    ),
    Open(
        "#ISSUE-IDEAL-TEE-REVERSED",
        "IdealMomentumTeeElement with an outlet arm reversed: 1.2e-6",
        lambda f: _cls(f) == "IdealMomentumTeeElement",
        Case("ideal_nport", 2, 1, 0),
        "jct pattern 0001 ratio 3",
        slow=True,
    ),
)


def open_issue(f: Finding) -> str | None:
    """The open pattern a mismatch belongs to, or None (a new one)."""
    return next((o.issue for o in OPEN if o.covers(f)), None)


# ---------------------------------------------------------------------------
# Classes the coverage gate expects
# ---------------------------------------------------------------------------


def network_classes() -> dict[str, type]:
    """Every element and node class exported from combaero.network (and the
    ejector, which lives in its own module)."""
    out: dict[str, type] = {"EjectorElement": EjectorElement}
    for name in cbn.__all__:
        obj = getattr(cbn, name)
        if isinstance(obj, type) and issubclass(obj, cbn.NetworkElement | cbn.NetworkNode):
            out[name] = obj
    return out


def table(rep: Report) -> str:
    """The coverage table (Markdown)."""
    lines = [
        "| class | reached | reversed | differing composition | differing temperature |",
        "|---|---:|---:|---:|---:|",
    ]
    known = network_classes()
    for cls in sorted(known) + sorted(k for k in rep.coverage if k not in known):
        c = rep.coverage.get(cls, Coverage())
        lines.append(f"| {cls} | {c.reached} | {c.reversed} | {c.composition} | {c.temperature} |")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--families", default="all", help="comma list, or 'all'")
    p.add_argument("--max-points", type=int, default=MAX_POINTS)
    p.add_argument("--seeds", default=f"{SEEDS.start}:{SEEDS.stop}", help="lo:hi")
    p.add_argument("--coverage-only", action="store_true")
    p.add_argument("--out", type=Path, help="JSON report")
    a = p.parse_args(argv)
    lo, hi = (int(v) for v in a.seeds.split(":"))
    cases = full_cases(range(lo, hi))
    if a.families != "all":
        cases = [c for c in cases if c.family in a.families.split(",")]
    t0 = time.time()
    rep = sweep(cases, max_points=a.max_points, coverage_only=a.coverage_only)
    dt = time.time() - t0
    print(table(rep))
    print(
        f"\n{len(cases)} cases, {rep.compared} points compared in {dt:.0f} s; excluded {rep.excluded}"
    )
    for r in rep.build_refused:
        print(f"  refused: {r}")
    for k, v in sorted(rep.nonfinite.items()):
        print(f"  nonfinite {v}x {k}")
    for f in sorted(rep.findings, key=lambda f: (f.verdict, -f.error)):
        print(
            f"{f.verdict:10s} {f.error:.2e} {f.element} [{f.row_kind}] col {f.column} | {f.case} | "
            f"{f.point} | est {f.fd_estimate:.1e}"
        )
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(
            json.dumps(
                {
                    "coverage": {k: vars(v) for k, v in rep.coverage.items()},
                    "findings": [vars(f) for f in rep.findings],
                    "excluded": rep.excluded,
                    "compared": rep.compared,
                    "build_refused": rep.build_refused,
                    "nonfinite": rep.nonfinite,
                    "cases": len(cases),
                    "family_seconds": rep.seconds,
                    "seconds": dt,
                },
                indent=1,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
