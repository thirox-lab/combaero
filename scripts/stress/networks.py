"""Stress-network generators: random networks built cold, over and past the
ranges the test suite exercises.

The test suite builds each element near its validated operating points. Cold
random networks find what it cannot: spurious roots reported as converged,
NaN states, misfiled residual rows, regimes with no root (#614, #619-#624).

Every generator takes a seed and returns ``(network, tags)``; FAMILIES maps a
family name to its generator. Some draws have no root by construction (a
legacy tee with a large straight-run drop, an ejector beyond its dead-head
back pressure); a failure is a finding to triage, not automatically a bug.

  random        layered networks of orifices, channels, area changes, losses,
                lossless links, plena, momentum chambers and combustors, with
                back edges, dead ends, injections, reversing branches, heat
                and thermal walls
  hotcold       hot and cold high-Mach channels coupled by a conductive wall
  hotcold_mild  the same with a wall that moves only modest heat
  tee3          a 3-port junction (MPCE, constant-K or ideal), merge or branch
  ideal_nport   an ideal momentum junction with 4-6 ports
  ejector       an ejector between three channels, compressible or jet pump
  legacy_tee    the Bassett TeeJunctionElement
  impingement   an impingement array, optionally under a wall-coupled hot duct
  impingement_bypass  an impingement array under a pressure-fed bypass crossflow
  effusion      an effusion plate, a channel-fed chain of plates, or a liner
  pin_fin       a pin-fin channel, optionally wall-coupled to a hot channel
  bleed_station a crossflow bleed or injection station
  series        vortex, connection, loss and orifice elements in series
  jet_chain     a crossflow chain fed by inclined jets (#813), from a reservoir
                or a closed end
  splitter      a flow splitter (#903), 2-6 outlets at drawn fractions, fed
                from a reservoir through an orifice or at a fixed flow
  fuel_controller  a FuelBoundary on a combustor (#892): lean, rich, near-peak and
                unreachable targets, fuel routes, a wall or turbine downstream;
                each draw tags its expected outcome (expect:...)
"""

from __future__ import annotations

import math
import random
from collections.abc import Callable

import combaero as cb
from combaero.network import (
    AreaChangeElement,
    ChannelElement,
    CombustorNode,
    ConstantKTeeElement,
    ConvectiveSurface,
    CrossflowSegmentElement,
    DiameterDischargeCoefficientConnectionElement,
    EffectiveAreaConnectionElement,
    EffusionLiner,
    EffusionPlateElement,
    EnergyBoundary,
    FlowNetwork,
    FlowSplitter,
    FuelBoundary,
    HoleCd,
    IdealMomentumTeeElement,
    ImpingementArray,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    MultiPortChamberElement,
    OrificeElement,
    PinFinModel,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
    SmoothModel,
    TeeJunctionElement,
    ThermalWall,
    TurbineElement,
    VortexElement,
    WallLayer,
    WallNode,
)
from combaero.network.components import ORIFICE_BORE_PIPE_RTOL
from combaero.network.ejector_element import EjectorElement
from combaero.network.graph import _flow_independent_dpt
from combaero.network.pressure_loss import ConstantFractionLoss, ConstantHeadLoss

Builder = Callable[[int], tuple[FlowNetwork, list[str]]]

Y_AIR = list(cb.species.dry_air_mass())
Y_CH4 = list(cb.species.pure_species("CH4"))
Y_MIX = [0.97 * a + 0.03 * b for a, b in zip(Y_AIR, Y_CH4, strict=True)]
_Y_N2 = [0.0] * len(Y_AIR)
_Y_N2[cb.species.indices["N2"]] = 0.95
_Y_N2[cb.species.indices["O2"]] = 0.05
COMPOSITIONS = {
    "air": Y_AIR,
    "air_ch4": [0.95 * a + 0.05 * b for a, b in zip(Y_AIR, Y_CH4, strict=True)],
    "n2rich": _Y_N2,
}


def _pb(g: FlowNetwork, nid: str, Pt: float, Tt: float = 300.0, Y: list | None = None) -> str:
    g.add_node(PressureBoundary(nid, Pt=Pt, Tt=Tt, Y=Y or Y_AIR))
    return nid


# ---------------------------------------------------------------------------
# random: layered networks of the basic elements
# ---------------------------------------------------------------------------


class _Layered:
    """Inlets, 1-5 layers of 1-3 interior nodes, outlets; every node gets at
    least one inflow and one outflow, then the extras. Node kinds are
    decided last (a momentum chamber needs exactly one in and one out)."""

    def __init__(self, seed: int) -> None:
        self.r = random.Random(seed)
        self.net = FlowNetwork()
        self.n = 0
        self.tags: list[str] = []
        self.surfaced: list[str] = []
        self.pending: dict[str, str] = {}
        self.specs: list[tuple[str, str, str, str | None]] = []
        self.heat: tuple[str, float] | None = None

    def nid(self, prefix: str) -> str:
        self.n += 1
        return f"{prefix}{self.n}"

    def node(self, kind: str) -> str:
        nid = self.nid("n")
        self.pending[nid] = kind
        return nid

    def element(self, a: str, b: str, kind: str | None = None) -> None:
        self.specs.append((self.nid("e"), a, b, kind))

    def _make_node(self, nid: str, kind: str) -> None:
        if kind == "plenum":
            node = PlenumNode(nid)
        elif kind == "mcn":
            node = MomentumChamberNode(nid, area=self.r.choice([1e-3, 4e-3, 2e-2]))
        else:
            node = CombustorNode(nid, method="complete")
        self.net.add_node(node)

    def _make_element(self, eid: str, a: str, b: str, kind: str | None) -> None:
        r = self.r
        kind = kind or r.choice(
            ["orifice", "orifice", "channel", "channel", "channel_hm", "area", "loss", "lossless"]
        )
        na, nb = self.net.nodes[a], self.net.nodes[b]
        both_mcn = isinstance(na, MomentumChamberNode) and isinstance(nb, MomentumChamberNode)
        if kind == "area" and not both_mcn:
            kind = "orifice"
        if kind == "lossless" and isinstance(na, PlenumNode) and isinstance(nb, PlenumNode):
            kind = "orifice"
        if kind == "orifice":
            d = r.choice([0.002, 0.01, 0.03, 0.08, 0.3])
            e = OrificeElement(
                eid,
                a,
                b,
                Cd=r.uniform(0.55, 0.9),
                diameter=d,
                correlation="fixed",
                regime=r.choice(["incompressible", "compressible"]),
            )
        elif kind in ("channel", "channel_hm"):
            hm = kind == "channel_hm"
            D = r.choice([0.004, 0.006]) if hm else r.choice([0.01, 0.03, 0.1])
            L = r.choice([0.02, 0.3, 2.0, 8.0])
            surf = ConvectiveSurface(area=math.pi * D * L) if r.random() < 0.5 else None
            e = ChannelElement(
                eid,
                a,
                b,
                length=L,
                diameter=D,
                roughness=r.choice([0.0, 1e-5, 5e-4]),
                regime="compressible" if (hm or r.random() < 0.5) else "incompressible",
                surface=surf,
            )
            if surf is not None:
                self.surfaced.append(eid)
            if hm:
                self.tags.append("high_mach_channel")
        elif kind == "area":
            model = r.choice(["sharp", "conical"])
            if na.area == nb.area:
                # Equal faces are no area change (#731): a lossless joint,
                # after the same draw. Since #1023 the element would be the
                # same pass-through; the joint stays, so these seeds keep
                # the network they drew.
                e = LosslessConnectionElement(eid, a, b)
            else:
                e = AreaChangeElement(
                    eid, a, b, F0=na.area, F1=nb.area, model_type=model, length=0.05
                )
                self.tags.append("area_change")
        elif kind == "loss":
            # Between two combustors the network refuses an undeclared theta
            # source, whatever the correlation; name the upstream one. No
            # draw, so no seed changes.
            both_theta = getattr(na, "has_theta", False) and getattr(nb, "has_theta", False)
            e = PressureLossElement(
                eid,
                a,
                b,
                correlation=ConstantFractionLoss(xi=r.choice([0.01, 0.05, 0.2])),
                theta_source=a if both_theta else None,
            )
        else:
            e = LosslessConnectionElement(eid, a, b)
        if _flow_independent_dpt(e) is not None and _closes_a_flow_independent_loop(self.net, a, b):
            # A loop of elements whose Pt relation ignores their flow (lossless
            # joints, fraction losses without m_ref) has no root or an
            # undetermined split, and the network refuses it (#988 and its
            # follow-up): a head loss instead, which scales with the flow at
            # every flow (a fraction loss with m_ref is flat above it). No
            # draw, so no other seed changes (of 0:200, random 5, 14, 19, 34,
            # 45, 47, 48, 81, 88, 95, 99, 113, 156 and 192 change). The same
            # for a path of them between two pressure boundaries (#1003): of
            # 0:200, random 21, 37, 39, 43, 48, 49, 52, 58, 60, 75, 78, 89,
            # 121, 125, 130, 131, 139, 148, 151, 178, 187 and 196 change.
            both_theta = getattr(na, "has_theta", False) and getattr(nb, "has_theta", False)
            e = PressureLossElement(
                eid,
                a,
                b,
                correlation=ConstantHeadLoss(zeta=_LOOP_BREAK_ZETA, area=_LOOP_BREAK_AREA),
                area=_LOOP_BREAK_AREA,
                theta_source=a if both_theta else None,
            )
        self.net.add_element(e)

    def build(self) -> FlowNetwork:
        r, net = self.r, self.net
        # Draw order is part of a seed's identity: keep it, or every recorded
        # finding stops reproducing.
        P_hi = r.choice([1.3e5, 2e5, 5e5, 2e6])
        n_in, n_out = r.randint(1, 3), r.randint(1, 3)
        n_layers, width = r.randint(1, 5), r.randint(1, 3)
        inlets = []
        for _ in range(n_in):
            Y = COMPOSITIONS[r.choice(list(COMPOSITIONS))]
            Tt = r.choice([250.0, 300.0, 600.0, 1200.0])
            nid = self.nid("in")
            if r.random() < 0.35:
                net.add_node(MassFlowBoundary(nid, m_dot=r.choice([0.005, 0.05, 0.5]), Tt=Tt, Y=Y))
                self.tags.append("mfb_inlet")
            else:
                net.add_node(PressureBoundary(nid, Pt=P_hi * r.uniform(0.9, 1.1), Tt=Tt, Y=Y))
            inlets.append(nid)
        layers = [inlets]
        for _ in range(n_layers):
            layer = []
            for _ in range(width):
                kind = r.choices(["plenum", "mcn", "comb"], [0.5, 0.4, 0.1])[0]
                layer.append(self.node(kind))
                if kind == "comb":
                    self.tags.append("combustor")
                    f = self.nid("fuel")
                    net.add_node(
                        MassFlowBoundary(f, m_dot=r.choice([0.001, 0.01]), Tt=300.0, Y=Y_CH4)
                    )
                    self.element(f, layer[-1], "orifice")
            layers.append(layer)
        outlets = []
        for _ in range(n_out):
            nid = self.nid("out")
            net.add_node(PressureBoundary(nid, Pt=1e5 * r.uniform(0.95, 1.05), Tt=300.0, Y=Y_AIR))
            outlets.append(nid)
        layers.append(outlets)
        for up, dn in zip(layers, layers[1:], strict=False):
            for b in dn:
                self.element(r.choice(up), b)
            for a in up:
                if not any(sp[1] == a for sp in self.specs) or r.random() < 0.3:
                    self.element(a, r.choice(dn))
        mids = [n for layer in layers[1:-1] for n in layer]
        if len(layers) > 3 and r.random() < 0.3:
            a, b = r.choice(layers[-2]), r.choice(layers[1])
            if a != b:
                self.element(a, b, "orifice")
                self.tags.append("back_edge")
        if mids and r.random() < 0.25:
            w = self.nid("wall")
            net.add_node(WallNode(w))
            self.element(r.choice(mids), w, r.choice(["orifice", "channel"]))
            self.tags.append("wall_dead_end")
        if mids and r.random() < 0.3:
            m = self.nid("inj")
            net.add_node(
                MassFlowBoundary(
                    m,
                    m_dot=r.choice([0.002, 0.02, 0.2]),
                    Tt=r.choice([300.0, 900.0]),
                    Y=r.choice(list(COMPOSITIONS.values())),
                )
            )
            a = r.choice(mids)
            self.element(a, m, "orifice")
            self.element(m, r.choice([n for n in mids + outlets if n != a]), "orifice")
            self.tags.append("injection")
        if mids and r.random() < 0.25:
            rv = self.nid("rev")
            net.add_node(PressureBoundary(rv, Pt=P_hi * 1.3, Tt=r.choice([300.0, 700.0]), Y=Y_AIR))
            self.element(r.choice(mids), rv, "orifice")
            self.tags.append("reverse_branch")
        if mids and r.random() < 0.4:
            node = r.choice(mids)
            Q = r.choice([5e3, 5e4, -2e3, -2e4])
            self.heat = (node, Q)
            self.tags.append("heat+" if Q > 0 else "heat-")
        for nid, kind in self.pending.items():
            n_in = sum(sp[2] == nid for sp in self.specs)
            n_out = sum(sp[1] == nid for sp in self.specs)
            if kind == "mcn" and not (n_in == 1 and n_out == 1):
                kind = "plenum"
            self._make_node(nid, kind)
            if kind == "mcn":
                self.tags.append("mcn")
        for sp in self.specs:
            self._make_element(*sp)
        if self.heat:
            p, Q = self.heat
            net.nodes[p].add_energy_boundary(EnergyBoundary(self.nid("q"), Q=Q))
        if len(self.surfaced) >= 2 and r.random() < 0.7:
            a, b = r.sample(self.surfaced, 2)
            net.add_wall(
                ThermalWall(
                    id=self.nid("w"),
                    element_a=a,
                    element_b=b,
                    layers=[WallLayer(0.002, r.choice([1.0, 25.0, 200.0]))],
                )
            )
            self.tags.append("thermal_wall")
        _fit_orifices_into_their_chambers(net)
        _lump_orifices_at_or_above_their_pipe(net)
        return net


#: The head loss that replaces an element closing a flow-independent loop:
#: zeta 10 keeps it below M ~0.2 at a 20 % drop; the area is a mid-size duct.
_LOOP_BREAK_ZETA = 10.0
_LOOP_BREAK_AREA = 4e-3


def _closes_a_flow_independent_loop(net: FlowNetwork, a: str, b: str) -> bool:
    """Whether a and b are already joined by elements whose Pt relation does
    not depend on their flow alone (``_flow_independent_dpt``: lossless
    joints, fraction losses without m_ref; #988). Every PressureBoundary
    fixes its Pt, so all of them count as one node, as in
    ``FlowNetwork._flow_independent_loop`` (#1003): two boundaries are
    always joined."""

    def key(n: str) -> str:
        return "<pressure boundaries>" if isinstance(net.nodes.get(n), PressureBoundary) else n

    a, b = key(a), key(b)
    adj: dict[str, list[str]] = {}
    for e in net.elements.values():
        if _flow_independent_dpt(e) is not None:
            adj.setdefault(key(e.from_node), []).append(key(e.to_node))
            adj.setdefault(key(e.to_node), []).append(key(e.from_node))
    seen, stack = {a}, [a]
    while stack:
        n = stack.pop()
        if n == b:
            return True
        for nb in adj.get(n, []):
            if nb not in seen:
                seen.add(nb)
                stack.append(nb)
    return False


#: Cd A over the chamber area a bore drawn at or above its momentum chamber is
#: shrunk to (#1016): the ratio the issue's experiment used, so the bore's
#: own area is 0.56-0.91 of the chamber's at the family's Cd 0.55-0.9.
_BORE_CHAMBER_CDA_RATIO = 0.5


def _fit_orifices_into_their_chambers(net: FlowNetwork) -> None:
    """Fit the draws' bores into the momentum chambers they adjoin (#1016).
    The family draws bore and chamber area independently, so a bore can be
    wider than the chamber it opens into or draws from: a geometry the
    network now refuses (equal is no orifice). Shrink such a bore to Cd A =
    _BORE_CHAMBER_CDA_RATIO times the smaller adjoining chamber area, with no
    random draw, so every other seed draws the network it drew before. Runs
    before _lump_orifices_at_or_above_their_pipe, which then sees the fitted
    bore."""
    for e in net.elements.values():
        if type(e) is not OrificeElement or e.diameter is None:
            continue
        areas = [
            n.area
            for n in (net.nodes[e.from_node], net.nodes[e.to_node])
            if isinstance(n, MomentumChamberNode) and n.area
        ]
        if not areas or math.sqrt(e.area / min(areas)) < 1.0 - ORIFICE_BORE_PIPE_RTOL:
            continue
        e.area = _BORE_CHAMBER_CDA_RATIO * min(areas) / e.Cd
        e.diameter = math.sqrt(4.0 * e.area / math.pi)


def _lump_orifices_at_or_above_their_pipe(net: FlowNetwork) -> None:
    """Declare the draws' lumped restrictions (#788). The family draws an
    orifice as a lumped Cd*A restriction, so its bore may meet or exceed a
    neighbouring channel's diameter. Until 0.8.0 the element silently read
    that as a wall hole (beta = 0); now a bore equal to its pipe is no orifice
    and a larger one is refused. Mark those draws in_pipe=False -- the same
    wall hole, said -- with no random draw, so every seed draws the network
    it drew before. The pipes are the ones that set beta: the single
    upstream channel forward, the single to-side channel in reverse. Rows
    are bit-identical except where a draw sat between a pipe narrower than
    its bore and a wider one: that was a wall hole one way and an in-pipe
    orifice the other, and is a wall hole both ways now (17 of seeds
    0-299)."""
    for e in net.elements.values():
        if type(e) is not OrificeElement or e.diameter is None:
            continue
        up = [c for c in net.get_upstream_elements(e.from_node) if isinstance(c, ChannelElement)]
        to = [
            c
            for c in net.get_upstream_elements(e.to_node) + net.get_downstream_elements(e.to_node)
            if isinstance(c, ChannelElement)
        ]
        pipes = [g[0].diameter for g in (up, to) if len(g) == 1 and g[0].diameter]
        if any(e.diameter >= D * (1.0 - ORIFICE_BORE_PIPE_RTOL) for D in pipes):
            e.in_pipe = False


def layered_random(seed: int) -> tuple[FlowNetwork, list[str]]:
    gen = _Layered(seed)
    net = gen.build()
    return net, gen.tags


def _hot_cold(
    r: random.Random,
    D: float,
    L: float,
    P: float,
    Th: float,
    Tc: float,
    pr: float,
    area: float,
    wall: WallLayer,
    roughness: tuple[float, float],
) -> FlowNetwork:
    g = FlowNetwork()
    _pb(g, "hi", P, Th)
    _pb(g, "ci", P, Tc)
    _pb(g, "ho", P * pr)
    _pb(g, "co", P * pr)
    g.add_node(MomentumChamberNode("hm", area=math.pi * D * D / 4))
    g.add_node(MomentumChamberNode("cm", area=math.pi * D * D / 4))
    for s, rough in zip(("h", "c"), roughness, strict=True):
        g.add_element(
            ChannelElement(
                f"{s}ch",
                f"{s}i",
                f"{s}m",
                length=L,
                diameter=D,
                roughness=rough,
                regime="compressible",
                surface=ConvectiveSurface(area=area),
            )
        )
        g.add_element(
            ChannelElement(
                f"{s}x",
                f"{s}m",
                f"{s}o",
                length=0.02,
                diameter=D,
                roughness=0.0,
                regime="compressible",
            )
        )
    g.add_wall(ThermalWall(id="w", element_a="hch", element_b="cch", layers=[wall]))
    return g


def hot_cold(seed: int) -> tuple[FlowNetwork, list[str]]:
    """Strong heating on one side, strong cooling on the other, Mach 0.3-0.7."""
    r = random.Random(seed)
    D, L = r.choice([0.004, 0.006, 0.008]), r.choice([0.2, 0.5, 1.0])
    g = _hot_cold(
        r,
        D,
        L,
        r.choice([3e5, 1e6]),
        r.choice([1200.0, 1600.0]),
        r.choice([300.0, 400.0]),
        r.choice([0.55, 0.7, 0.85]),
        math.pi * D * L * r.choice([1.0, 5.0]),
        WallLayer(0.0005, 200.0),
        (0.0, 0.0),
    )
    return g, ["hotcold"]


def hot_cold_mild(seed: int) -> tuple[FlowNetwork, list[str]]:
    """Heated and cooled channels at Mach 0.2-0.6 through a weak wall."""
    r = random.Random(1000 + seed)
    D, L = r.choice([0.008, 0.012, 0.02]), r.choice([0.3, 1.0])
    P, Th, pr = r.choice([3e5, 1e6]), r.choice([900.0, 1400.0]), r.choice([0.6, 0.75, 0.9])
    rough = (r.choice([0.0, 2e-5]), r.choice([0.0, 2e-5]))
    wall = WallLayer(r.choice([0.002, 0.005]), r.choice([1.0, 10.0]))
    g = _hot_cold(r, D, L, P, Th, 300.0, pr, math.pi * D * L, wall, rough)
    return g, ["hotcold_mild"]


# ---------------------------------------------------------------------------
# element families
# ---------------------------------------------------------------------------


def _arm(r: random.Random, eid: str, a: str, b: str, D: float, regime: str) -> tuple:
    """A random outer element for a junction port."""
    kind = r.choice(["lossless", "channel", "channel", "orifice"])
    if kind == "lossless":
        return LosslessConnectionElement(eid, a, b), kind
    if kind == "channel":
        return ChannelElement(
            eid,
            a,
            b,
            length=r.choice([0.05, 0.5, 5.0]),
            diameter=D,
            roughness=r.choice([0.0, 1e-5]),
            regime=regime,
        ), kind
    return OrificeElement(
        eid,
        a,
        b,
        Cd=r.uniform(0.6, 0.9),
        diameter=D * r.uniform(0.5, 0.95),
        correlation="fixed",
        regime=regime,
    ), kind


def tee3(seed: int) -> tuple[FlowNetwork, list[str]]:
    r = random.Random(seed)
    g = FlowNetwork()
    model = r.choice(["mpce", "constk", "ideal"])
    direction = r.choice(["branch", "merge"])
    regime = r.choice(["compressible", "incompressible"])
    D = r.choice([0.01, 0.03, 0.08])
    psi = math.exp(r.uniform(math.log(0.5), math.log(6.0)))
    A = math.pi * D * D / 4
    A_bra, D_bra = A / psi, D / math.sqrt(psi)
    P = r.choice([1.5e5, 3e5, 1e6])
    Tt = r.choice([250.0, 300.0, 600.0, 900.0])
    dP = P * r.choice([0.005, 0.02, 0.08, 0.2])
    theta = r.choice([30.0, 60.0, 90.0, 135.0])
    tags = [f"tee:{model}", direction, regime]
    for pid, a in (("pc", A), ("ps", A), ("pb", A_bra)):
        g.add_node(MomentumChamberNode(pid, area=a, regime=regime))
    mfb = r.random() < 0.4
    if direction == "branch":
        if mfb:
            g.add_node(
                MassFlowBoundary(
                    "bc", m_dot=r.choice([0.01, 0.1, 1.0]) * A / 7.85e-4, Tt=Tt, Y=Y_AIR
                )
            )
            tags.append("mfb")
        else:
            _pb(g, "bc", P, Tt)
        _pb(g, "bs", P - dP * r.uniform(0.2, 1.0))
        _pb(g, "bb", P - dP * r.uniform(0.2, 1.2))
        wiring = (("bc", "pc", D), ("ps", "bs", D), ("pb", "bb", D_bra))
        ins, outs, ang_in, ang_out, areas = ["pc"], ["ps", "pb"], [0.0], [0.0, theta], [A, A, A_bra]
    else:
        Tb = r.choice([300.0, Tt, 1200.0])
        if mfb:
            g.add_node(MassFlowBoundary("bs", m_dot=0.05 * A / 7.85e-4, Tt=Tt, Y=Y_AIR))
            g.add_node(MassFlowBoundary("bb", m_dot=0.03 * A_bra / 7.85e-4, Tt=Tb, Y=Y_MIX))
            tags.append("mfb")
        else:
            _pb(g, "bs", P, Tt)
            _pb(g, "bb", P * r.uniform(0.98, 1.05), Tb, Y_MIX)
        _pb(g, "bc", P - dP)
        wiring = (("bs", "ps", D), ("bb", "pb", D_bra), ("pc", "bc", D))
        ins, outs, ang_in, ang_out, areas = ["ps", "pb"], ["pc"], [0.0, theta], [0.0], [A, A_bra, A]
    for i, (a, b, d) in enumerate(wiring):
        el, kind = _arm(r, f"arm{i}", a, b, d, regime)
        g.add_element(el)
        tags.append(f"arm:{kind}")
    common = {
        "id": "jct",
        "inlet_nodes": ins,
        "outlet_nodes": outs,
        "inlet_angles_deg": ang_in,
        "outlet_angles_deg": ang_out,
        "port_areas": areas,
    }
    if model == "mpce":
        j = MultiPortChamberElement(**common, flow_direction=direction, strict=False)
    elif model == "constk":
        K = (
            {1: r.uniform(0.0, 1.5), 2: r.uniform(0.2, 2.0)}
            if direction == "branch"
            else {0: r.uniform(0.0, 1.5), 1: r.uniform(0.2, 2.0)}
        )
        j = ConstantKTeeElement(**common, flow_direction=direction, K_ports=K)
    else:
        j = IdealMomentumTeeElement(**common, loss_factor=r.choice([0.0, 0.5, 1.0, 1.5]))
    g.add_element(j)
    return g, tags


def ideal_nport(seed: int) -> tuple[FlowNetwork, list[str]]:
    r = random.Random(seed)
    g = FlowNetwork()
    n_in = r.randint(1, 3)
    n_out = r.randint(4 - n_in if n_in < 3 else 1, 6 - n_in)
    regime = r.choice(["compressible", "incompressible"])
    D = r.choice([0.02, 0.05])
    A = math.pi * D * D / 4
    P = r.choice([2e5, 5e5])
    ins, outs = [], []
    for i in range(n_in):
        g.add_node(MomentumChamberNode(f"pi{i}", area=A))
        _pb(g, f"bi{i}", P * r.uniform(0.99, 1.02), r.choice([300.0, 700.0]))
        g.add_element(
            ChannelElement(f"ci{i}", f"bi{i}", f"pi{i}", length=0.5, diameter=D, regime=regime)
        )
        ins.append(f"pi{i}")
    for i in range(n_out):
        g.add_node(MomentumChamberNode(f"po{i}", area=A))
        _pb(g, f"bo{i}", P * r.uniform(0.9, 0.99))
        g.add_element(
            ChannelElement(f"co{i}", f"po{i}", f"bo{i}", length=0.5, diameter=D, regime=regime)
        )
        outs.append(f"po{i}")
    g.add_element(
        IdealMomentumTeeElement(
            id="jct",
            inlet_nodes=ins,
            outlet_nodes=outs,
            inlet_angles_deg=[r.choice([0.0, 45.0, 90.0]) for _ in ins],
            outlet_angles_deg=[r.choice([0.0, 45.0, 90.0, 135.0]) for _ in outs],
            loss_factor=r.choice([0.0, 1.0]),
        )
    )
    return g, ["ideal_nport", f"ports:{n_in}+{n_out}", regime]


def ejector(seed: int) -> tuple[FlowNetwork, list[str]]:
    """Port chambers carry their boundary pressure as initial_guess; no
    other seeding, so the element's own warm start is what is exercised."""
    r = random.Random(seed)
    g = FlowNetwork()
    regime = r.choice(["compressible", "incompressible"])
    D = 0.05
    A = math.pi * D * D / 4
    At = r.choice([1.5e-5, 3.14e-5, 6e-5])
    Ae = At * r.uniform(1.2, 4.0)
    # A mixing duct wider than its outlet pipe is refused (#1023): such a
    # draw is fitted to the pipe, with no random draw, so every other seed
    # draws what it drew before (of 0:100, ejector 1, 11, 68 and 93 change).
    Am = min(Ae * r.uniform(3.0, 12.0), A)
    if regime == "compressible":
        Pp = r.choice([3e5, 7e5, 1.5e6])
        Ps = Pp * r.uniform(0.02, 0.2)
        Po = Ps * r.uniform(1.1, 2.5)
    else:
        Pp = r.choice([2e5, 4e5])
        Ps = 1e5
        # Inside the jet pump's operating range (#619): below its shut-off
        # head, Sanger's rows at zero suction flow with no losses,
        # dPt_o = 2 (Pt_p - Pt_s) b (1 - b/2), b = A_e / A_3. Above it the
        # suction reverses, which the model was not derived for.
        b = Ae / Am
        Po = Ps + 2.0 * (Pp - Ps) * b * (1.0 - 0.5 * b) * r.uniform(0.05, 0.8)
    for nid, p, tt in (("bp", Pp, r.choice([300.0, 440.0])), ("bs", Ps, 280.0), ("bo", Po, 300.0)):
        _pb(g, nid, p, tt)
    for nid, p in (("mp", Pp), ("ms", Ps), ("mo", Po)):
        mc = MomentumChamberNode(nid, area=A, regime=regime)
        mc.initial_guess = {f"{nid}.P": p, f"{nid}.Pt": p}
        g.add_node(mc)
    for eid, a, b in (("cp", "bp", "mp"), ("cs", "bs", "ms"), ("co", "mo", "bo")):
        g.add_element(
            ChannelElement(
                eid, a, b, length=r.choice([0.05, 0.2, 1.0]), diameter=D, regime="incompressible"
            )
        )
    g.add_element(
        EjectorElement(
            id="ej",
            primary_node="mp",
            secondary_node="ms",
            outlet_node="mo",
            throat_area=At,
            nozzle_exit_area=Ae,
            mixing_area=Am,
            regime=regime,
        )
    )
    return g, ["ejector", regime]


def legacy_tee(seed: int) -> tuple[FlowNetwork, list[str]]:
    """The Bassett tee; its straight run caps at K ~ 0.3, so a large imposed
    straight drop has no root."""
    r = random.Random(seed)
    g = FlowNetwork()
    P = r.choice([1.5e5, 2.1e5, 6e5])
    F_C = r.choice([2e-3, 1e-2, 3e-2])
    psi = r.choice([1.0, 2.0, 4.0])
    theta = math.radians(r.choice([30.0, 60.0, 90.0]))
    s = r.choice([0.002, 0.01, 0.04])
    Tt = r.choice([300.0, 400.0, 800.0])
    if r.random() < 0.5:
        _pb(g, "c", P, Tt)
        _pb(g, "s", P * (1 - 0.05 * s))
        _pb(g, "b", P * (1 - s))
        g.add_element(
            TeeJunctionElement(
                id="tee",
                common_node="c",
                straight_node="s",
                branch_node="b",
                theta=theta,
                F_C=F_C,
                psi=psi,
                tee_type="branching",
            )
        )
        return g, ["legacy_tee", "branching"]
    _pb(g, "is", P, Tt)
    _pb(g, "ib", P * (1 - 0.5 * s), r.choice([300.0, 900.0]))
    _pb(g, "o", P * (1 - s))
    g.add_node(MomentumChamberNode("ps", area=F_C))
    g.add_node(MomentumChamberNode("pb", area=F_C / psi))
    g.add_node(PlenumNode("pc"))
    g.add_element(LosslessConnectionElement("ls", "is", "ps"))
    g.add_element(LosslessConnectionElement("lb", "ib", "pb"))
    g.add_element(LosslessConnectionElement("lc", "pc", "o"))
    g.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="pc",
            straight_node="ps",
            branch_node="pb",
            theta=theta,
            F_C=F_C,
            psi=psi,
            tee_type="merging",
        )
    )
    return g, ["legacy_tee", "merging"]


def impingement(seed: int) -> tuple[FlowNetwork, list[str]]:
    r = random.Random(seed)
    g = FlowNetwork()
    D = r.choice([0.00127, 0.00254, 0.00508])
    xn, yn, zd = r.choice([(5.0, 4.0, 1.0), (5.0, 4.0, 2.0), (10.0, 8.0, 2.0), (15.0, 4.0, 3.0)])
    n = r.randint(1, 10)
    span = 12 * yn * D
    P = r.choice([1.0e5, 3e5, 1e6])
    drive = r.choice([0.005, 0.03, 0.15, 0.5])
    _pb(g, "sup", P * (1 + drive), r.choice([300.0, 600.0]))
    _pb(g, "ex", P)
    ImpingementArray(
        "ia",
        n_rows=n,
        d_jet=D,
        xn_d=xn,
        yn_d=yn,
        z_d=zd,
        span=span,
        plate_thickness=D,
        pattern=r.choice(["inline", "staggered"]),
        hole_cd=HoleCd(correlation=r.choice(["fixed", "IdelchikThick", "McGreehanSchotsch"])),
    ).add_to(g, "sup", "ex")
    tags = ["impingement", f"rows:{n}"]
    if r.random() < 0.5:
        _pb(g, "hi", P * 1.02, r.choice([900.0, 1500.0]))
        _pb(g, "ho", P)
        prev = "hi"
        for i in range(n):
            nxt = f"h{i}" if i < n - 1 else "ho"
            if nxt != "ho":
                g.add_node(PlenumNode(nxt))
            g.add_element(
                ChannelElement(
                    f"g{i}",
                    prev,
                    nxt,
                    length=xn * D,
                    diameter=0.02,
                    roughness=0.0,
                    # The duct's own wetted area: the plate's width here
                    # (xn D span, 3.9x) gave it heat over surface its friction
                    # never sees, a cooling-driven pump (#551).
                    surface=ConvectiveSurface(area=math.pi * 0.02 * xn * D),
                )
            )
            g.add_wall(
                ThermalWall(
                    id=f"w{i}",
                    element_a=f"g{i}",
                    element_b=f"ia__p{i + 1}",
                    layers=[WallLayer(thickness=0.001, conductivity=r.choice([1.0, 20.0, 200.0]))],
                )
            )
            prev = nxt
        tags.append("walls")
    return g, tags



def impingement_bypass(seed: int) -> tuple[FlowNetwork, list[str]]:
    """An impingement array under an initial (bypass) crossflow (#467, #800).

    The bypass is drawn from a feed at a pressure between the exit's and
    above the supply's, so the upstream rows run from strong jets to
    stagnant or reversed ones (Gc/Gj up to ~100): the regime where the
    holes' crossflow Cd acts and where a constant Cd let the stagnant rows
    alternate in sign (#799). The impingement family has no bypass and never
    leaves the Cd's cap.
    """
    r = random.Random(seed)
    g = FlowNetwork()
    D = r.choice([0.00127, 0.00254])
    xn, yn, zd = r.choice([(5.0, 4.0, 1.0), (5.0, 4.0, 2.0), (5.0, 8.0, 1.0), (10.0, 8.0, 1.0)])
    n = r.randint(3, 10)
    P = r.choice([1.0e5, 3e5])
    drive = r.choice([0.005, 0.03, 0.15])
    feed = r.choice([0.3, 0.8, 1.0, 1.1])
    _pb(g, "sup", P * (1 + drive), 300.0)
    _pb(g, "ex", P)
    _pb(g, "byp", P * (1 + feed * drive), r.choice([300.0, 320.0]))
    ImpingementArray(
        "ia",
        n_rows=n,
        d_jet=D,
        xn_d=xn,
        yn_d=yn,
        z_d=zd,
        span=12 * yn * D,
        plate_thickness=D,
        hole_cd=HoleCd(correlation=r.choice(["fixed", "IdelchikThick"])),
        bypass_entry_K=r.choice([0.0, 0.5]),
    ).add_to(g, "sup", "ex", bypass="byp")
    return g, ["impingement_bypass", f"rows:{n}", f"feed:{feed}"]

def effusion(seed: int) -> tuple[FlowNetwork, list[str]]:
    r = random.Random(seed)
    g = FlowNetwork()
    kind = r.choice(["plate", "chain", "liner"])
    P = r.choice([1.0e5, 5e5, 2e6])
    drive = r.choice([0.01, 0.04, 0.2])
    Tg = r.choice([300.0, 1500.0, 2000.0])
    geo = r.choice(
        [
            {"hole_diameter": 2.16e-3, "wall_thickness": 6.3e-3, "pitch": 15.24e-3},
            {"hole_diameter": 0.6e-3, "wall_thickness": 1e-3, "pitch": 3e-3},
            {"hole_diameter": 3.27e-3, "wall_thickness": 6.3e-3, "pitch": 15.24e-3},
        ]
    )
    ang = r.choice([30.0, 60.0, 90.0])
    if kind == "plate":
        _pb(g, "c", P * (1 + drive), r.choice([300.0, 700.0]))
        _pb(g, "gas", P, Tg)
        g.add_element(
            EffusionPlateElement(
                "panel", "c", "gas", panel_area=r.choice([4e-4, 0.01, 0.3]), angle_deg=ang, **geo
            )
        )
    elif kind == "chain":
        _pb(g, "feed", P * (1 + drive), 600.0)
        _pb(g, "dump", P * (1 + 0.1 * drive))
        prev = "feed"
        for i in range(r.randint(2, 6)):
            g.add_node(PlenumNode(f"n{i}"))
            g.add_element(ChannelElement(f"s{i}", prev, f"n{i}", length=0.02, diameter=6e-3))
            _pb(g, f"g{i}", P, Tg)
            g.add_element(
                EffusionPlateElement(
                    f"p{i}", f"n{i}", f"g{i}", panel_area=4e-4, angle_deg=ang, **geo
                )
            )
            prev = f"n{i}"
        g.add_element(ChannelElement("sx", prev, "dump", length=0.02, diameter=6e-3))
    else:
        _pb(g, "cin", P * (1 + drive), 600.0)
        _pb(g, "cout", P * (1 + 0.8 * drive))
        _pb(g, "gas", P, Tg)
        EffusionLiner(
            "ln",
            n_segments=r.randint(1, 6),
            length=0.152,
            width=0.152,
            duct_height=r.choice([0.01, 0.03]),
            hole_diameter=3.27e-3,
            wall_thickness=6.3e-3,
            pitch_x=15.24e-3,
            pitch_y=15.24e-3,
        ).add_to(g, "cin", "cout", "gas")
    return g, ["effusion", kind]


def pin_fin(seed: int) -> tuple[FlowNetwork, list[str]]:
    r = random.Random(seed)
    g = FlowNetwork()
    P = r.choice([2e5, 1.5e6])
    drive = r.choice([0.002, 0.013, 0.1])
    _pb(g, "a", P * (1 + drive), r.choice([300.0, 700.0]))
    _pb(g, "b", P)
    model = PinFinModel(
        arrangement=r.choice([cb.PinArrangement.Staggered, cb.PinArrangement.Inline])
    )
    g.add_element(
        ChannelElement(
            "ch",
            "a",
            "b",
            length=r.choice([0.02, 0.05, 0.2]),
            diameter=r.choice([0.005, 0.009, 0.02]),
            regime=r.choice(["compressible", "incompressible"]),
            surface=ConvectiveSurface(area=0.01, model=model),
        )
    )
    tags = ["pin_fin"]
    if r.random() < 0.6:
        _pb(g, "hi", P * (1 + drive), r.choice([1000.0, 1600.0]))
        _pb(g, "ho", P)
        g.add_element(
            ChannelElement(
                "hot",
                "hi",
                "ho",
                length=0.05,
                diameter=0.02,
                surface=ConvectiveSurface(area=0.01, model=SmoothModel(correlation="gnielinski")),
            )
        )
        g.add_wall(
            ThermalWall(
                id="w",
                element_a="hot",
                element_b="ch",
                layers=[WallLayer(thickness=0.002, conductivity=r.choice([2.0, 25.0]))],
            )
        )
        tags.append("wall")
    return g, tags


def bleed_station(seed: int) -> tuple[FlowNetwork, list[str]]:
    r = random.Random(seed)
    g = FlowNetwork()
    regime = r.choice(["compressible", "incompressible"])
    P = r.choice([1e5, 5e5])
    drive = r.choice([100.0, 1e3, 1e4])
    bleed = r.random() < 0.6
    kappa = cb.STATION_KAPPA_BLEED_BASSETT if bleed else cb.STATION_KAPPA_MERGE_NORMAL
    L = r.choice([1e-6, 0.05, 0.5])
    _pb(g, "in", P + drive)
    _pb(g, "out", P)
    _pb(g, "lat", P - drive if bleed else P + 2 * drive, r.choice([300.0, 600.0]))
    g.add_node(PlenumNode("b"))
    g.add_element(
        CrossflowSegmentElement(
            "s0",
            "in",
            "b",
            length=L,
            area=1e-3,
            entry_K=0.0,
            to_kappa=kappa,
            next_seg="s1",
            regime=regime,
        )
    )
    g.add_element(
        CrossflowSegmentElement(
            "s1", "b", "out", length=L, area=1e-3, from_kappa=kappa, prev_seg="s0", regime=regime
        )
    )
    a, b = ("b", "lat") if bleed else ("lat", "b")
    g.add_element(
        OrificeElement(
            "lat_o",
            a,
            b,
            diameter=r.choice([0.005, 0.01, 0.03]),
            Cd=0.6,
            correlation="fixed",
            regime=regime,
        )
    )
    return g, ["bleed_station", "bleed" if bleed else "inject", regime]



def jet_chain(seed: int) -> tuple[FlowNetwork, list[str]]:
    """A crossflow chain fed by inclined jets (#813): each station's hole
    declares m w cos(theta), which the station takes into its momentum.

    1-10 stations, each with its own angle (90 = normal, no momentum); the
    chain's upstream end is a reservoir entry or a closed end (the jets alone
    drive it, an impingement gap's topology). Exercises the station's jet
    term, its Jacobian and the ledger's grouping of momentum-bearing holes.
    """
    r = random.Random(seed)
    g = FlowNetwork()
    regime = r.choice(["compressible", "incompressible"])
    n = r.randint(1, 10)
    P = r.choice([1e5, 5e5, 2e6])
    drive = r.choice([0.005, 0.03, 0.15, 0.5])
    A = r.choice([2e-4, 2e-3, 2e-2])
    L = r.choice([0.005, 0.02, 0.1])
    d_hole = r.choice([0.002, 0.005, 0.01]) * math.sqrt(A / 2e-3)
    closed = r.random() < 0.4
    _pb(g, "sup", P * (1.0 + drive), r.choice([300.0, 600.0]))
    _pb(g, "ex", P)
    if not closed:
        _pb(g, "res", P * (1.0 + drive * r.choice([0.0, 0.3, 1.0])), r.choice([300.0, 900.0]))
    for i in range(1, n + 1):
        g.add_node(PlenumNode(f"b{i}"))
    if not closed:
        g.add_element(
            CrossflowSegmentElement(
                "s0", "res", "b1", length=L, area=A, entry_K=r.choice([0.0, 0.5]),
                to_kappa=0.0, next_seg="s1", regime=regime,
            )
        )  # fmt: skip
    angles = []
    for i in range(1, n + 1):
        last = i == n
        g.add_element(
            CrossflowSegmentElement(
                f"s{i}", f"b{i}", "ex" if last else f"b{i + 1}", length=L, area=A,
                from_kappa=0.0, prev_seg=None if (closed and i == 1) else f"s{i - 1}",
                to_kappa=None if last else 0.0, next_seg=None if last else f"s{i + 1}",
                regime=regime,
            )
        )  # fmt: skip
        theta = r.choice([90.0, 60.0, 30.0, 15.0])
        angles.append(theta)
        g.add_element(
            OrificeElement(
                f"h{i}", "sup", f"b{i}", diameter=d_hole, Cd=0.7, correlation="fixed",
                injection_angle_deg=theta, regime=regime,
            )
        )  # fmt: skip
    tags = ["jet_chain", f"stations:{n}", regime, "closed" if closed else "reservoir"]
    tags.append("angled" if any(a < 90.0 for a in angles) else "normal")
    return g, tags

def series(seed: int) -> tuple[FlowNetwork, list[str]]:
    r = random.Random(seed)
    g = FlowNetwork()
    P = r.choice([1.2e5, 3e5, 1e6])
    _pb(g, "a", P, r.choice([300.0, 800.0]))
    _pb(g, "b", P * r.uniform(0.5, 0.97))
    nodes = ["a"] + [f"n{i}" for i in range(r.randint(2, 5) - 1)] + ["b"]
    for nid in nodes[1:-1]:
        g.add_node(PlenumNode(nid))
    kinds = []
    for i, (x, y) in enumerate(zip(nodes, nodes[1:], strict=False)):
        k = r.choice(["vortex", "ea", "dcd", "dcd_zeta", "loss", "orifice"])
        if i == 0 and k in ("vortex", "loss"):
            k = "orifice"  # something must set the flow
        if k == "vortex":
            el = VortexElement(
                f"e{i}",
                x,
                y,
                shaft_speed_rpm=r.choice([500.0, 3000.0, 15000.0]),
                r_out=r.choice([0.05, 0.1, 0.3]),
                r_in=r.choice([0.0, 0.03]),
            )
            # The retired Vatistas n (#661), still drawn so every later draw
            # in the seed keeps its value.
            r.choice([1.0, 2.0, 4.0])
        elif k == "ea":
            el = EffectiveAreaConnectionElement(f"e{i}", x, y, r.choice([1e-5, 1e-3, 2e-2]))
        elif k == "dcd":
            el = DiameterDischargeCoefficientConnectionElement(
                f"e{i}", x, y, r.choice([0.003, 0.01, 0.05]), Cd=r.uniform(0.6, 0.95)
            )
        elif k == "dcd_zeta":
            el = DiameterDischargeCoefficientConnectionElement(
                f"e{i}", x, y, r.choice([0.003, 0.01, 0.05]), zeta=r.uniform(0.0, 5.0)
            )
        elif k == "loss":
            el = PressureLossElement(
                f"e{i}", x, y, correlation=ConstantFractionLoss(xi=r.choice([0.01, 0.1]))
            )
        else:
            el = OrificeElement(
                f"e{i}",
                x,
                y,
                Cd=0.8,
                diameter=r.choice([0.005, 0.02]),
                correlation="fixed",
                regime=r.choice(["compressible", "incompressible"]),
            )
        g.add_element(el)
        kinds.append(k)
    return g, ["series", *sorted(set(kinds))]


def splitter(seed: int) -> tuple[FlowNetwork, list[str]]:
    """A flow splitter (#903): 2-6 outlets at drawn fractions (sometimes one
    the remainder), each through an orifice or a channel to the sink. The
    reference is the path the incompressible seed estimate says needs the
    most pressure at its fraction; a draw where the estimate is wrong is
    refused by the second-law check (a finding to triage). Fed from a
    reservoir through an orifice, or at a fixed flow sized for a reference
    drop of 5-50% of the sink pressure."""
    r = random.Random(seed)
    g = FlowNetwork()
    regime = r.choice(["compressible", "incompressible"])
    supply = r.choice(["pressure", "mass"])
    n = r.randint(2, 6)
    T = r.choice([300.0, 700.0])
    P_out = r.choice([1e5, 3e5])
    _pb(g, "out", P_out)
    w = [r.uniform(0.2, 1.0) for _ in range(n)]
    fr = [x / sum(w) for x in w]
    outs = [f"o{i}" for i in range(n)]
    arms = []
    for i, o in enumerate(outs):
        g.add_node(PlenumNode(o))
        D = r.choice([0.005, 0.01, 0.02])
        el, kind = _arm(r, f"a{i}", o, "out", D, regime)
        if kind == "lossless":  # something must take the drop
            el = OrificeElement(
                f"a{i}", o, "out", Cd=0.8, diameter=D, correlation="fixed", regime=regime
            )
        if isinstance(el, OrificeElement):
            el.in_pipe = False  # a lumped restriction, whatever channels share the sink
        g.add_element(el)
        arms.append(el)
    rho = float(cb.density(T, P_out, cb.species.dry_air()))
    K = [float(e._seed_resistance(rho, P_out, 0.01) or 0.0) for e in arms]
    ref = max(range(n), key=lambda i: K[i] * fr[i] ** 2)
    remainder = r.random() < 0.5
    fractions: dict[str, float | None] = dict(zip(outs, fr, strict=True))
    if remainder:
        fractions[outs[ref]] = None
    if supply == "pressure":
        _pb(g, "res", P_out * r.uniform(1.05, 1.6), T)
        g.add_node(PlenumNode("in"))
        g.add_element(
            OrificeElement(
                "feed", "res", "in", Cd=0.8, diameter=0.02, correlation="fixed", regime=regime
            )
        )
    else:
        dP = P_out * r.uniform(0.05, 0.5)
        m_ref = math.sqrt(dP / K[ref]) if K[ref] > 0.0 else 0.01
        g.add_node(MassFlowBoundary("in", m_dot=m_ref / fr[ref], Tt=T, Y=Y_AIR))
    FlowSplitter("sp", "in", fractions, reference=outs[ref]).add_to(g)
    tags = ["splitter", supply, regime, f"outlets:{n}"]
    return g, tags + (["remainder"] if remainder else [])


# ---------------------------------------------------------------------------
# fuel_controller: a FuelBoundary on a combustor (#892, #1049)
# ---------------------------------------------------------------------------

#: The cases, one per seed in turn (seed % len), so every range of
#: len(FUEL_CASES) consecutive seeds draws each once.
FUEL_CASES = (
    "lean",  # a target between the unfuelled Tt and the peak, lean root
    "rich",  # the same on the rich side (rich-quench-lean, reformers)
    "near_peak",  # a few K below the peak, where the two roots meet
    "below",  # below the unfuelled Tt: no root
    "above",  # above the peak: no root
    "rich_no_root",  # rich asked, only a lean root within phi 10
    "leak",  # part of the fuel bypasses the combustor
    "coupled_wall",  # the burnt gas heats a coolant channel through a wall
    "coupled_turbine",  # the burnt gas drives a turbine at a fixed pressure ratio
)
#: The fuel's routes to the combustor, in turn per len(FUEL_CASES) seeds.
FUEL_ROUTES = ("direct", "injector", "premix")
#: The expected-outcome tags run.py checks: converged on the branch, or
#: refused (an exception out of solve(), or no success) with the message part.
EXPECT_CONVERGED = "expect:converged:"
EXPECT_REFUSED = "expect:refused:"
#: A refusal before any solve (resolve_topology): the harness test asserts it.
EXPECT_REFUSED_AT_BUILD = "expect:refused_at_build:"

_FUEL_SPECIES = ("CH4", "CH4", "H2", "C3H8")


def _stream(m: float, T: float, P: float, Y: list[float]) -> cb.Stream:
    st = cb.Stream()
    st.set_T(T)
    st.set_P(P)
    st.set_X(list(cb.mass_to_mole(Y)))
    st.set_mdot(m)
    return st


def _burnt_at_phi(
    phi: float, T_air: float, T_fuel: float, P: float, Y_f: list[float], method: str
) -> float:
    """Tt of 1 kg/s of dry air with fuel at equivalence ratio phi, burnt by the
    combustor's model (mix, then complete or partial-equilibrium combustion),
    as the solver's fuel start computes it."""
    o2 = Y_AIR[cb.species.indices["O2"]]
    m_st = o2 / float(cb.oxygen_required_per_kg_mixture(list(cb.mass_to_mole(Y_f))))
    if phi <= 0.0:
        return T_air
    mixed = cb.mix([_stream(1.0, T_air, P, Y_AIR), _stream(phi * m_st, T_fuel, P, Y_f)])
    if method == "equilibrium":
        return float(cb.combustion_equilibrium(mixed.T, list(mixed.X), mixed.P).state.T)
    return float(cb.complete_combustion(mixed.T, list(mixed.X), mixed.P).T)


def fuel_controller(seed: int) -> tuple[FlowNetwork, list[str]]:
    """A FuelBoundary controlling a combustor's exit Tt (#892, #1049).

    Air (a reservoir through an orifice, or a fixed flow) and a fuel
    controller meet in a combustor, which exits to an outlet through an
    orifice. The case (FUEL_CASES, by seed) sets the target from the
    combustor's own model at the draw's air and fuel: a phi on the branch, a
    few K under the peak, or past either end of the reachable window. The
    route (FUEL_ROUTES): the fuel straight in, through an injector orifice, or
    into a premix plenum the air also enters. The coupled cases put a
    wall-cooled channel or a fixed-PR turbine downstream.

    Every draw declares its expected outcome as a tag (EXPECT_*): converged
    on the asked branch, or refused with a part of the reason.
    """
    case = FUEL_CASES[seed % len(FUEL_CASES)]
    route = FUEL_ROUTES[(seed // len(FUEL_CASES)) % len(FUEL_ROUTES)]
    r = random.Random(seed)
    g = FlowNetwork()
    method = r.choice(["complete", "equilibrium"])
    fuel = r.choice(_FUEL_SPECIES)
    Y_f = list(cb.species.to_mass(cb.species.pure_species(fuel)))
    T_air = r.choice([300.0, 450.0, 650.0, 800.0])
    T_fuel = r.choice([300.0, 400.0])
    supply = r.choice(["pressure", "mass"])
    P = r.choice([1.5e5, 5e5, 2e6])
    if case == "coupled_turbine":
        P = max(P, 5e5)
    m_air = r.choice([0.05, 0.5, 5.0]) * P / 5e5
    branch = r.choice(["lean", "rich"]) if case in ("near_peak", "below", "above") else "lean"
    if case in ("rich", "rich_no_root") or (case.startswith("coupled") and r.random() < 0.3):
        branch = "rich"

    def burnt(phi: float) -> float:
        return _burnt_at_phi(phi, T_air, T_fuel, P, Y_f, method)

    T_peak = max(burnt(1.0 + 0.005 * k) for k in range(-10, 11))
    T_rich_end = burnt(10.0)  # the solver's rich bracket ends at phi 10
    tags = [
        "fuel_controller", f"case:{case}", f"route:{route}", f"method:{method}",
        f"fuel:{fuel}", f"branch:{branch}", f"supply:{supply}",
    ]  # fmt: skip
    if case in ("lean", "rich", "coupled_wall", "coupled_turbine", "leak"):
        phi = r.uniform(0.25, 0.85) if branch == "lean" else r.uniform(1.2, 2.5)
        target = burnt(phi)
        expect = EXPECT_CONVERGED + branch
    elif case == "near_peak":
        target = T_peak - r.choice([1.0, 3.0, 10.0, 30.0])
        expect = EXPECT_CONVERGED + branch
    elif case == "below":
        target = T_air - r.choice([2.0, 20.0, 150.0])
        expect = EXPECT_REFUSED + "with no injected fuel"
    elif case == "above":
        target = T_peak + r.choice([3.0, 20.0, 300.0])
        expect = EXPECT_REFUSED + "above the peak"
    else:  # rich_no_root
        lo, hi = T_air, T_rich_end
        if hi - lo > 40.0:
            target = lo + r.uniform(0.2, 0.8) * (hi - lo)
            expect = EXPECT_REFUSED + "not the rich one asked for"
        else:  # the rich branch reaches below the oxidizer: a rich root
            target = burnt(r.uniform(1.2, 2.5))
            expect = EXPECT_CONVERGED + "rich"
            tags.append("rich_root_exists")
    if case == "leak":
        expect = EXPECT_REFUSED_AT_BUILD + "without passing"
    tags.append(expect)

    if supply == "pressure":
        _pb(g, "air", P, T_air)
    else:
        g.add_node(MassFlowBoundary("air", m_dot=m_air, Tt=T_air, Y=Y_AIR))
    g.add_node(
        FuelBoundary("fuel", combustor="burner", Tt_target=target, branch=branch, Tt=T_fuel, Y=Y_f)
    )
    g.add_node(CombustorNode("burner", method=method))
    # Air orifice sized for ~ m_air at a 4 % drop (rho ~ P / (R T)).
    rho = P / (287.0 * T_air)
    d_air = math.sqrt(4.0 * m_air / (0.8 * math.pi * math.sqrt(2.0 * rho * 0.04 * P)))
    into = "burner"
    if route == "premix":
        g.add_node(PlenumNode("premix"))
        into = "premix"
    if supply == "pressure":
        g.add_element(
            OrificeElement(
                "a_in", "air", into, Cd=0.8, diameter=d_air, correlation="fixed", in_pipe=False
            )
        )
    else:
        g.add_element(LosslessConnectionElement("a_in", "air", into))
    if case == "leak":
        g.add_node(PlenumNode("fsplit"))
        g.add_element(LosslessConnectionElement("inj", "fuel", "fsplit"))
        g.add_element(
            OrificeElement(
                "f_main", "fsplit", into, Cd=0.7, diameter=0.004, correlation="fixed", in_pipe=False
            )
        )
    elif route == "direct":
        g.add_element(LosslessConnectionElement("inj", "fuel", into))
    else:
        g.add_element(
            OrificeElement(
                "inj",
                "fuel",
                into,
                Cd=0.7,
                diameter=r.choice([0.002, 0.005, 0.01]) * math.sqrt(m_air / 0.5),
                correlation="fixed",
                in_pipe=False,
            )
        )
    if route == "premix":
        g.add_element(
            OrificeElement(
                "holder", "premix", "burner", Cd=0.8, diameter=1.5 * d_air,
                correlation="fixed", in_pipe=False,
            )
        )  # fmt: skip

    # Downstream: an exit orifice sized for a ~10 % drop at the burnt density.
    T_out_guess = max(target, T_air)
    d_exit = 1.2 * d_air * (T_out_guess / T_air) ** 0.25
    P_out = P * r.uniform(0.75, 0.9)
    if case == "coupled_wall":
        D, L = 2.0 * d_exit, r.choice([0.1, 0.5])
        _pb(g, "out", P_out)
        g.add_element(
            ChannelElement(
                "hot", "burner", "out", length=L, diameter=D, regime="compressible",
                surface=ConvectiveSurface(area=math.pi * D * L),
            )
        )  # fmt: skip
        _pb(g, "ci", P * 1.05, 300.0)
        _pb(g, "co", P)
        g.add_element(
            ChannelElement(
                "cold", "ci", "co", length=L, diameter=D, regime="compressible",
                surface=ConvectiveSurface(area=math.pi * D * L),
            )
        )  # fmt: skip
        g.add_wall(
            ThermalWall(
                id="w",
                element_a="hot",
                element_b="cold",
                layers=[WallLayer(0.002, r.choice([1.0, 25.0, 200.0]))],
            )
        )
    elif case == "coupled_turbine":
        PR = r.choice([2.0, 3.0, 4.0])
        g.add_node(PlenumNode("t_out"))
        _pb(g, "out", min(P_out / PR, 0.8 * P / PR))
        g.add_element(TurbineElement("turb", "burner", "t_out", PR, eta_poly=0.9))
        g.add_element(
            OrificeElement(
                "nozzle", "t_out", "out", Cd=0.8, diameter=2.0 * d_exit,
                correlation="fixed", in_pipe=False,
            )
        )  # fmt: skip
    else:
        _pb(g, "out", P_out)
        g.add_element(
            OrificeElement(
                "exit", "burner", "out", Cd=0.8, diameter=d_exit, correlation="fixed",
                in_pipe=False,
            )
        )  # fmt: skip
    if case == "leak":
        # The bypass: part of the fuel reaches the outlet unburnt.
        g.add_element(
            OrificeElement(
                "f_leak", "fsplit", "out", Cd=0.7, diameter=0.002, correlation="fixed",
                in_pipe=False,
            )
        )  # fmt: skip
    return g, tags


FAMILIES: dict[str, Builder] = {
    "random": layered_random,
    "hotcold": hot_cold,
    "hotcold_mild": hot_cold_mild,
    "tee3": tee3,
    "ideal_nport": ideal_nport,
    "ejector": ejector,
    "legacy_tee": legacy_tee,
    "impingement": impingement,
    "impingement_bypass": impingement_bypass,
    "effusion": effusion,
    "pin_fin": pin_fin,
    "bleed_station": bleed_station,
    "series": series,
    "jet_chain": jet_chain,
    "splitter": splitter,
    "fuel_controller": fuel_controller,
}
