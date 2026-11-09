"""A ``FlowNetwork`` as a GUI network file (#902).

``network_to_gui(net)`` is the inverse of ``graph_builder.build_network_from_schema``
for the node and element classes listed in ``SUPPORTED``: a network built in
a script opens in the GUI, and ``NetworkRunner.from_dict(network_to_gui(net))``
solves to the same results. Anything without a GUI form is refused with the
list of what and why -- never dropped.

Phase 1 covers the boundaries (the fuel controller too), plenum, momentum chamber, combustor, channel,
orifice, pressure loss (all four correlations), lossless connection, the
compressor and turbine, the flow splitter (#903: one node over its outlet
elements) and the thermal wall, with smooth convective surfaces.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any

import combaero as cb
from combaero.network import (
    CasingRadiation,
    ChannelElement,
    CombustorNode,
    FlowNetwork,
    FuelBoundary,
    GasRadiation,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
    SplitterOutletElement,
    ThermalWall,
    TurbomachineElement,
)
from combaero.network.components import OrificeElement, SmoothModel
from combaero.network.pressure_loss import (
    ConstantFractionLoss,
    ConstantHeadLoss,
    LinearThetaFractionLoss,
    LinearThetaHeadLoss,
)

FILE_VERSION = "1.2"

#: The handle a flow edge uses at each end, per GUI node type: (outlet, inlet).
#: Kept equal to gui/frontend/src/utils/nodeHandles.generated.json by
#: python/tests/test_gui_export.py (the frontend sources do not ship in the
#: wheel, so the export cannot read the table at run time).
FLOW_HANDLES: dict[str, tuple[str | None, str | None]] = {
    "pressure_boundary": ("flow-source", "flow-target"),
    "mass_boundary": ("flow-source", "flow-target"),
    # A fuel controller (#892) only supplies: no inlet.
    "fuel_boundary": ("flow-source", None),
    "plenum": ("flow-right", "flow-left"),
    "momentum_chamber": ("flow-source", "flow-target"),
    "combustor": ("flow-source", "flow-target"),
    "channel": ("flow-source", "flow-target"),
    "orifice": ("flow-source", "flow-target"),
    "discrete_loss": ("flow-source", "flow-target"),
    "lossless_connection": (None, None),
    "compressor": ("flow-source", "flow-target"),
    "turbine": ("flow-source", "flow-target"),
}
#: A flow splitter's (#903) ports: one inlet, outlet i at SPLITTER_OUT.format(i).
SPLITTER_IN = "port-in-target"
SPLITTER_OUT = "port-out{}-source"
#: A merge chamber's side streams arrive at this handle (graph_builder).
SIDE_INLET = "side-target"

SUPPORTED = (
    "PressureBoundary",
    "MassFlowBoundary",
    "FuelBoundary",
    "PlenumNode",
    "MomentumChamberNode",
    "CombustorNode",
    "ChannelElement",
    "OrificeElement",
    "PressureLossElement",
    "LosslessConnectionElement",
    "CompressorElement",
    "TurbineElement",
    "SplitterOutletElement",
    "ThermalWall",
)


class UnsupportedForGUIError(ValueError):
    """The network holds objects the GUI cannot represent; ``problems`` lists them."""

    def __init__(self, problems: list[str]) -> None:
        super().__init__(
            "network_to_gui: no GUI form for:\n  - "
            + "\n  - ".join(problems)
            + f"\nSupported: {', '.join(SUPPORTED)} (smooth surfaces)."
        )
        self.problems = problems


def _composition(Y: list[float] | None) -> dict[str, Any]:
    """Exact: the mass fractions by species name ("custom", mass mode)."""
    if Y is None:
        return {"source": "dry_air", "mode": "mole", "relative_humidity": 0}
    fractions = {n: float(y) for n, y in zip(cb.species.names, Y, strict=True) if y > 0.0}
    return {"source": "custom", "mode": "mass", "custom_fractions": fractions}


def _surface(obj: Any, owner: str, problems: list[str], *, own_area: bool) -> dict[str, Any]:
    s = getattr(obj, "surface", None)
    if s is None:
        return {}
    if not isinstance(s.model, SmoothModel):
        problems.append(f"{owner}: {type(s.model).__name__} surface (phase 1 exports smooth only)")
    elif s.model.correlation != "gnielinski" or s.model.roughness or s.model.mu_ratio != 1.0:
        problems.append(f"{owner}: SmoothModel with non-default settings (the GUI has none)")
    if own_area and s.area:
        problems.append(
            f"{owner}: an explicit convective area ({s.area:g} m^2); the GUI derives it "
            "from the geometry"
        )
    return {"Nu_multiplier": float(s.Nu_multiplier), "f_multiplier": float(s.f_multiplier)}


def _regime(obj: Any) -> str:
    r = getattr(obj, "regime", None)
    return r if r in ("incompressible", "compressible") else "default"


def _loss(elem: PressureLossElement, owner: str, problems: list[str]) -> dict[str, Any]:
    c = elem.correlation
    if isinstance(c, ConstantFractionLoss):
        d: dict[str, Any] = {"correlation_type": "constant_fraction", "xi": float(c.xi)}
        if c.m_ref is not None:
            d["m_ref"] = float(c.m_ref)
    elif isinstance(c, LinearThetaFractionLoss):
        d = {"correlation_type": "linear_theta_fraction", "k": float(c.k), "xi0": float(c.xi0)}
        if c.m_ref is not None:
            d["m_ref"] = float(c.m_ref)
    elif isinstance(c, ConstantHeadLoss):
        d = {"correlation_type": "constant_head", "zeta": float(c.zeta), "area": float(c.area)}
    elif isinstance(c, LinearThetaHeadLoss):
        d = {
            "correlation_type": "linear_theta_head",
            "k": float(c.k),
            "zeta0": float(c.zeta0),
            "area": float(c.area),
        }
    else:
        problems.append(f"{owner}: loss correlation {type(c).__name__}")
        return {}
    if elem.theta_source:
        d["theta_source"] = elem.theta_source
    return d


def _node(nid: str, node: Any, problems: list[str]) -> tuple[str, dict[str, Any]] | None:
    if isinstance(node, PressureBoundary):
        if getattr(node, "coupling", "auto") != "auto":
            problems.append(
                f"{nid}: PressureBoundary coupling={node.coupling!r} (the GUI uses auto)"
            )
        return "pressure_boundary", {
            "Pt": float(node.Pt),
            "Tt": float(node.Tt),
            "composition": _composition(node.Y),
        }
    # Before MassFlowBoundary, its base: exported as one it would lose the
    # controller and fix the fuel flow at whatever it last was.
    if isinstance(node, FuelBoundary):
        d = {
            "combustor_id": node.combustor,
            "Tt_target": float(node.Tt_target),
            "branch": node.branch,
            "Tt": float(node.Tt),
            "composition": _composition(node.Y),
        }
        if node.m_dot_guess is not None:
            d["m_dot_guess"] = float(node.m_dot_guess)
        return "fuel_boundary", d
    if isinstance(node, MassFlowBoundary):
        return "mass_boundary", {
            "m_dot": float(node.m_dot),
            "Tt": float(node.Tt),
            "composition": _composition(node.Y),
        }
    if isinstance(node, MomentumChamberNode):
        d: dict[str, Any] = _surface(node, nid, problems, own_area=False)
        if not node._auto_area:
            d["area"] = float(node.area)
        if node.Dh is not None:
            d["Dh"] = float(node.Dh)
        if node.length is not None:
            problems.append(f"{nid}: MomentumChamberNode length (no GUI field)")
        return "momentum_chamber", d
    if isinstance(node, CombustorNode):
        d = {"method": node.method, **_surface(node, nid, problems, own_area=False)}
        if not node._auto_area:
            d["area"] = float(node.area)
        if getattr(node, "Dh", None) is not None:
            d["Dh"] = float(node.Dh)
        if getattr(node, "energy_boundaries", None):
            problems.append(
                f"{nid}: energy boundaries on a combustor (export them as heat sources)"
            )
        return "combustor", d
    if type(node) is PlenumNode:
        if getattr(node, "energy_boundaries", None):
            problems.append(f"{nid}: energy boundaries on a plenum (export them as heat sources)")
        return "plenum", {}
    problems.append(f"{nid}: node type {type(node).__name__}")
    return None


def _element(eid: str, el: Any, problems: list[str]) -> tuple[str, dict[str, Any]] | None:
    if isinstance(el, ChannelElement):
        if el.diameter is None:
            problems.append(f"{eid}: a channel without a diameter (export it after resolving)")
            return None
        d: dict[str, Any] = {
            "L": float(el.length),
            "D": float(el.diameter),
            "roughness": float(el.roughness),
            "friction_model": el.friction_model,
            "regime": _regime(el),
            **_surface(el, eid, problems, own_area=True),
        }
        if el.Dh is not None and abs(el.Dh - el.diameter) > 1e-15 * max(el.Dh, 1.0):
            d["Dh"] = float(el.Dh)
        return "channel", d
    if isinstance(el, OrificeElement):
        d = {
            "correlation": el.correlation,
            "plate_thickness": float(el.plate_thickness),
            "edge_radius": float(el.edge_radius),
            "bevel_depth": float(el.bevel_depth),
            "tap": el.tap,
            "in_pipe": bool(el.in_pipe),
            "regime": _regime(el),
            "compressible_model": el.compressible_model
            if el._compressible_model_explicit
            else "default",
            "Cv": float(el.Cv),
        }
        if el.diameter is not None:
            d["diameter"] = float(el.diameter)
        if el.correlation == "fixed":
            d["Cd"] = float(el.Cd)
        if el.injection_angle_deg != 90.0:
            problems.append(f"{eid}: orifice injection_angle_deg != 90 (no GUI field)")
        return "orifice", d
    if isinstance(el, PressureLossElement):
        d = _loss(el, eid, problems)
        # The GUI's loss element has ONE area: a head-loss correlation's and the
        # element's flow area are the same field there.
        c_area = getattr(el.correlation, "area", None)
        if c_area is not None and el.area != c_area:
            problems.append(
                f"{eid}: a head-loss area ({c_area:g} m^2) different from the element's "
                f"({el.area}); the GUI has one area -- give the element area={c_area:g}"
            )
        if el.area is not None:
            d["area"] = float(el.area)
        d.update(_surface(el, eid, problems, own_area=False))
        return "discrete_loss", d
    if isinstance(el, LosslessConnectionElement):
        return "lossless_connection", {}
    if isinstance(el, TurbomachineElement):
        polytropic = el.eta_poly is not None
        return ("turbine" if el.expansion else "compressor"), {
            "pressure_ratio": float(el.pressure_ratio),
            "efficiency_type": "polytropic" if polytropic else "isentropic",
            "efficiency": float(el.eta_poly if polytropic else el.eta_is),
        }
    problems.append(f"{eid}: element type {type(el).__name__}")
    return None


def _splitters(
    net: FlowNetwork, problems: list[str]
) -> dict[str, tuple[str, list[str], dict[str, Any]]]:
    """The flow splitters (#903) as GUI nodes: {splitter id: (inlet node, the
    outlet nodes in port order, the node's data)}. Each is one canvas node
    over its SplitterOutletElements, in the order of their ``outlets``."""
    out: dict[str, tuple[str, list[str], dict[str, Any]]] = {}
    for el in net.elements.values():
        if not isinstance(el, SplitterOutletElement) or el.splitter in out:
            continue
        sid = el.splitter
        if sid in net.nodes or sid in net.elements:
            problems.append(f"{sid}: a flow splitter whose id is also a node or element id")
            continue
        members = [net.elements.get(eid) for eid in el.outlets]
        if not all(isinstance(m, SplitterOutletElement) for m in members):
            problems.append(f"{sid}: a flow splitter whose outlets are not all in the network")
            continue
        refs = [i for i, m in enumerate(members) if m.reference]
        areas = [None if m.area is None else float(m.area) for m in members]
        data: dict[str, Any] = {
            "fractions": [float(m.fraction) for m in members],
            "reference": refs[0] if len(refs) == 1 else None,
            "auto_reference": bool(el.auto_reference),
        }
        if any(a is not None for a in areas):
            data["areas"] = areas
        out[sid] = (el.from_node, [m.to_node for m in members], data)
    return out


def _radiation(spec: GasRadiation | CasingRadiation | None) -> dict[str, Any] | None:
    """A face's radiation (#1071) as the GUI's WallRadiationData; None: off."""
    if isinstance(spec, GasRadiation):
        return {
            "kind": "gas",
            "beam_length": float(spec.beam_length),
            "emissivity_wall": float(spec.emissivity_wall),
            "luminosity": float(spec.luminosity),
            "T_gas": spec.T_gas,
            "fuel_air_ratio": spec.fuel_air_ratio,
            "pressure": spec.pressure,
        }
    if isinstance(spec, CasingRadiation):
        return {
            "kind": "casing",
            "area_ratio": float(spec.area_ratio),
            "emissivity_wall": float(spec.emissivity_wall),
            "emissivity_casing": float(spec.emissivity_casing),
            "T_casing": spec.T_casing,
        }
    return None


def _wall(w: ThermalWall, problems: list[str]) -> dict[str, Any]:
    if w.element_b is None:
        problems.append(f"{w.id}: a wall with no stream on side B (side_b, #970; no GUI field)")
    if w.contact_area is None and w.area_a is None:
        problems.append(
            f"{w.id}: a wall without contact_area (the GUI stores the area; give it explicitly)"
        )
    layers = [
        {
            "thickness": float(L.thickness),
            "conductivity": float(L.conductivity),
            "material": L.material,
        }
        for L in w.layers
    ]
    if any(L.extrapolate for L in w.layers):
        problems.append(f"{w.id}: a layer with extrapolate=True (no GUI field)")
    d: dict[str, Any] = {
        "type": "thermal",
        "layers": layers,
        "area": float(w.contact_area) if w.contact_area is not None else None,
        "R_fouling": float(w.R_fouling),
        **(
            {"area_a": float(w.area_a), "area_b": float(w.area_b), "wall_area": w.wall_area}
            if w.area_a is not None
            else {}
        ),
        "orientation": w.orientation,
    }
    if layers:
        d["thickness"] = layers[0]["thickness"]
        d["conductivity"] = layers[0]["conductivity"]
    if w.film is not None:
        d.update(
            film=w.film,
            film_side=w.film_side,
            film_start=float(w.film_start),
            film_end=float(w.film_end) if w.film_end is not None else None,
        )
    for side in ("a", "b"):
        rad = _radiation(getattr(w, f"radiation_{side}"))
        if rad is not None:
            d[f"radiation_{side}"] = rad
    return d


def _rows_layout(
    entities: list[str], edges: list[tuple[str, str]]
) -> dict[str, tuple[float, float]]:
    """Columns by longest path from a source (left to right), rows in order."""
    succ: dict[str, list[str]] = {e: [] for e in entities}
    indeg = dict.fromkeys(entities, 0)
    for a, b in edges:
        succ[a].append(b)
        indeg[b] += 1
    depth = dict.fromkeys(entities, 0)
    queue = [e for e in entities if indeg[e] == 0]
    seen = 0
    while queue:
        a = queue.pop(0)
        seen += 1
        for b in succ[a]:
            depth[b] = max(depth[b], depth[a] + 1)
            indeg[b] -= 1
            if indeg[b] == 0:
                queue.append(b)
    columns: dict[int, list[str]] = {}
    for e in entities:  # a cycle leaves its members at their reached depth
        columns.setdefault(depth[e], []).append(e)
    return {e: (220.0 * c, 140.0 * r) for c, col in columns.items() for r, e in enumerate(col)}


def network_to_gui(
    net: FlowNetwork,
    layout: str | Mapping[str, tuple[float, float]] | Callable[[str], tuple[float, float]] = "rows",
    labels: Mapping[str, str] | None = None,
    example: Mapping[str, Any] | None = None,
    solver_settings: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """The GUI network file of ``net`` (a dict; ``json.dumps`` it to save).

    layout: ``"rows"`` (columns by flow depth), a mapping id -> (x, y), or a
    callable id -> (x, y). labels: optional card labels by id. example: an
    ``example`` block (worked examples, #793). Raises ``UnsupportedForGUIError``
    listing every object without a GUI form.
    """
    labels = labels or {}
    problems: list[str] = []
    kinds: dict[str, str] = {}
    data: dict[str, dict[str, Any]] = {}
    for nid, node in net.nodes.items():
        r = _node(nid, node, problems)
        if r:
            kinds[nid], data[nid] = r
    splitters = _splitters(net, problems)
    for sid, (_, _, sdata) in splitters.items():
        kinds[sid], data[sid] = "flow_splitter", sdata
    for eid, el in net.elements.items():
        if isinstance(el, SplitterOutletElement):
            continue  # its splitter's node (above)
        r = _element(eid, el, problems)
        if r:
            kinds[eid], data[eid] = r
    walls = {wid: _wall(w, problems) for wid, w in net.walls.items()}
    if problems:
        raise UnsupportedForGUIError(problems)

    def side_inlet(eid: str, to: str) -> bool:
        node = net.nodes[to]
        return (
            isinstance(node, MomentumChamberNode)
            and node.main_inlet is not None
            and node.main_inlet != eid
        )

    flow: list[tuple[str, str, str | None, str | None]] = []
    for sid, (inlet, outs, _) in splitters.items():
        flow.append((inlet, sid, FLOW_HANDLES[kinds[inlet]][0], SPLITTER_IN))
        for i, b in enumerate(outs):
            th = SIDE_INLET if side_inlet(f"{sid}__{b}", b) else FLOW_HANDLES[kinds[b]][1]
            flow.append((sid, b, SPLITTER_OUT.format(i), th))
    for eid, el in net.elements.items():
        if isinstance(el, SplitterOutletElement):
            continue
        a, b = el.from_node, el.to_node
        flow.append((a, eid, FLOW_HANDLES[kinds[a]][0], FLOW_HANDLES[kinds[eid]][1]))
        th = SIDE_INLET if side_inlet(eid, b) else FLOW_HANDLES[kinds[b]][1]
        flow.append((eid, b, FLOW_HANDLES[kinds[eid]][0], th))

    order = (
        list(net.nodes)
        + list(splitters)
        + [e for e, el in net.elements.items() if not isinstance(el, SplitterOutletElement)]
    )
    if layout == "rows":
        pos = _rows_layout(order, [(a, b) for a, b, _, _ in flow])
    elif callable(layout):
        pos = {e: tuple(layout(e)) for e in order}
    else:
        missing = [e for e in order if e not in layout]
        if missing:
            raise ValueError(f"network_to_gui: no position for {missing}")
        pos = {e: tuple(layout[e]) for e in order}

    nodes = [
        {
            "id": e,
            "type": kinds[e],
            "position": {"x": float(pos[e][0]), "y": float(pos[e][1])},
            "data": {"id": e, **({"label": labels[e]} if e in labels else {}), **data[e]},
        }
        for e in order
    ]
    edges: list[dict[str, Any]] = [
        {
            "id": f"e_{a}_{b}",
            "source": a,
            "sourceHandle": sh,
            "target": b,
            "targetHandle": th,
            "type": "default",
            "data": {},
        }
        for a, b, sh, th in flow
    ]
    for wid, w in net.walls.items():
        edges.append(
            {
                "id": wid,
                "source": w.element_a,
                "sourceHandle": "thermal-source",
                "target": w.element_b,
                "targetHandle": "thermal-target",
                "type": "thermal",
                "data": walls[wid],
            }
        )
    doc: dict[str, Any] = {"version": FILE_VERSION}
    if example is not None:
        doc["example"] = dict(example)
    doc["nodes"] = nodes
    doc["edges"] = edges
    if solver_settings is not None:
        doc["solverSettings"] = dict(solver_settings)
    return doc


__all__ = ["FLOW_HANDLES", "SUPPORTED", "UnsupportedForGUIError", "network_to_gui"]
