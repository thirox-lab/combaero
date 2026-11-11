from __future__ import annotations

import inspect
import warnings
from collections import deque
from typing import TYPE_CHECKING, Any, Self

if TYPE_CHECKING:
    from .components import NetworkElement, NetworkNode, ThermalWall


class FlowNetwork:
    """
    Manages the topology of the flow network graph, storing nodes and elements.
    Provides methods for elements to discover their upstream/downstream neighbors.
    """

    def __init__(self) -> None:
        self.nodes: dict[str, NetworkNode] = {}
        self.elements: dict[str, NetworkElement] = {}

        # Thermal coupling walls dictionary
        self.walls: dict[str, ThermalWall] = {}
        # Film providers that are not one element (#695): an EffusionLiner
        # registers its panels' film here under its own id.
        self.film_providers: dict[str, Any] = {}
        self.thermal_coupling_enabled: bool = True

        # Directed graph mapping: node_id -> list of element_ids
        self._upstream_of_node: dict[str, list[str]] = {}
        self._downstream_of_node: dict[str, list[str]] = {}

    def add_node(self, node: NetworkNode) -> None:
        """Register a new node in the network."""
        if node.id in self.nodes:
            raise ValueError(f"Node '{node.id}' already exists in topology.")
        self.nodes[node.id] = node
        self._upstream_of_node[node.id] = []
        self._downstream_of_node[node.id] = []

    def add_element(self, element: NetworkElement) -> None:
        """
        Register a new flow element. Connects the element to its source and sink
        nodes, ensuring they exist in the network.
        """
        if element.id in self.elements:
            raise ValueError(f"Element '{element.id}' already exists in topology.")

        all_nodes = element.all_source_nodes() + element.all_sink_nodes()
        for nid in all_nodes:
            if nid not in self.nodes:
                raise ValueError(f"Unknown node '{nid}' for element '{element.id}'.")

        if len(all_nodes) != len(set(all_nodes)):
            raise ValueError(
                f"Element '{element.id}' has duplicate nodes among its source/sink nodes."
            )

        self.elements[element.id] = element

        for src in element.all_source_nodes():
            self._downstream_of_node[src].append(element.id)

        for snk in element.all_sink_nodes():
            self._upstream_of_node[snk].append(element.id)

    def add_wall(self, wall: ThermalWall) -> None:
        """Register a thermal coupling wall between two elements.

        Parameters
        ----------
        wall : ThermalWall
            Wall connection defining thermal coupling between two elements.

        Raises
        ------
        ValueError
            If wall ID already exists or if element_a/element_b don't exist.
        """
        if wall.id in self.walls:
            raise ValueError(f"Wall '{wall.id}' already exists in network.")

        if wall.element_a not in self.elements and wall.element_a not in self.nodes:
            raise ValueError(f"Unknown element_a/node_a '{wall.element_a}' for wall '{wall.id}'.")

        # element_b None: side B has no stream (side_b, #970), checked by the wall.
        if (
            wall.element_b is not None
            and wall.element_b not in self.elements
            and wall.element_b not in self.nodes
        ):
            raise ValueError(f"Unknown element_b/node_b '{wall.element_b}' for wall '{wall.id}'.")

        if wall.film is not None:
            provider = self.film_provider(wall.film)
            if provider is None:
                raise ValueError(
                    f"Wall '{wall.id}': film {wall.film!r} is not a film provider "
                    "(a FilmCoolingElement, SlotFilmElement, EffusionPlateElement or EffusionLiner)."
                )
            error = getattr(provider, "film_provider_error", lambda: None)()
            if error is not None:
                raise ValueError(f"Wall '{wall.id}': film {wall.film!r}: {error}.")
            from .components import film_sources  # deferred, as the imports above

            side = wall.film_side_id
            own = {provider.id}
            for src in film_sources(provider):
                own |= {src.element_id, src.supply_node}
            if side in own:
                raise ValueError(
                    f"Wall '{wall.id}': the film covers the GAS side downstream of "
                    f"the holes; {side!r} is the provider or its coolant supply (#386)."
                )
            if getattr(provider, "film_htc", "none") != "none" and not float(wall.film_end) > 0.0:
                raise ValueError(
                    f"Wall '{wall.id}': film {wall.film!r} carries a film-side h "
                    f"({provider.film_htc!r}), which is infinite at the slot lip; "
                    "give the wall a span (film_end > 0) (#962)."
                )

        self.walls[wall.id] = wall

    def film_provider(self, film_id: str) -> Any:
        """The film provider ``film_id`` names: an element with a film span
        average, or one registered in ``film_providers``; None otherwise."""
        el = self.elements.get(film_id)
        if el is not None and hasattr(el, "film_span_average"):
            return el
        return self.film_providers.get(film_id)

    def get_upstream_elements(self, node_id: str) -> list[NetworkElement]:
        """Return all elements feeding into the specified node."""
        if node_id not in self.nodes:
            raise KeyError(f"Node '{node_id}' not found.")
        return [self.elements[eid] for eid in self._upstream_of_node[node_id]]

    def get_downstream_elements(self, node_id: str) -> list[NetworkElement]:
        """Return all elements flowing out of the specified node."""
        if node_id not in self.nodes:
            raise KeyError(f"Node '{node_id}' not found.")
        return [self.elements[eid] for eid in self._downstream_of_node[node_id]]

    def resolve_all_topology(self) -> None:
        """
        Trigger the auto-discovery phase where all components inspect the graph
        to extract their needed operational parameters (e.g. Orifice discovering
        upstream pipe diameter).
        """
        for node in self.nodes.values():
            node.resolve_topology(self)
        for element in self.elements.values():
            element.resolve_topology(self)
        # Second pass: elements that inherit from other elements (e.g. TeeJunction.F_C from a
        # channel whose D was itself inherited) get a second chance after the first pass resolved
        # all direct values.
        for element in self.elements.values():
            element.resolve_topology(self)

    def _reachable_from(self, start_id: str, through_closed: bool = True) -> set[str]:
        """BFS over undirected adjacency to find all node IDs reachable from start_id.

        ``through_closed=False`` does not cross a closed restriction (a zero
        bore, #1024): what lies beyond one is reached only by another path."""
        from .components import _is_closed_restriction

        visited: set[str] = set()
        queue = [start_id]
        while queue:
            nid = queue.pop()
            if nid in visited:
                continue
            visited.add(nid)
            for eid in self._upstream_of_node[nid] + self._downstream_of_node[nid]:
                elem = self.elements[eid]
                if not through_closed and _is_closed_restriction(elem):
                    continue
                if nid in elem.all_sink_nodes():
                    neighbours = elem.all_source_nodes()
                else:
                    neighbours = elem.all_sink_nodes()
                for nb in neighbours:
                    if nb not in visited:
                        queue.append(nb)
        return visited

    def _flow_independent_loop(self) -> tuple[list[str], tuple[str, str] | None] | None:
        """One closed loop made only of elements whose total-pressure relation
        does not depend on their own flow (``_flow_independent_dpt``; #988 for
        lossless connections), or None.

        Every ``PressureBoundary`` fixes its Pt, so all of them count as one
        shared node (#1003): a path of such elements from one boundary to
        another closes a loop through the fixed pressures. Returns the loop's
        element ids and, for such a path, the two boundaries it joins (the
        ids then run in order from the first to the second); for a loop among
        the network's own nodes, or through a single boundary, None.

        Union-find over the nodes: such an edge whose ends are already joined
        by such edges closes a loop, and the tree path between its ends
        completes it."""
        from .components import PressureBoundary

        shared = "<pressure boundaries>"

        def key(n: str) -> str:
            return shared if isinstance(self.nodes.get(n), PressureBoundary) else n

        parent: dict[str, str] = {}

        def root(n: str) -> str:
            parent.setdefault(n, n)
            while parent[n] != n:
                parent[n] = parent[parent[n]]
                n = parent[n]
            return n

        tree: dict[str, list[tuple[str, str]]] = {}
        for eid, elem in self.elements.items():
            if _flow_independent_dpt(elem) is None:
                continue
            a, b = key(elem.from_node), key(elem.to_node)
            ra, rb = root(a), root(b)
            if ra != rb:
                parent[ra] = rb
                tree.setdefault(a, []).append((b, eid))
                tree.setdefault(b, []).append((a, eid))
                continue
            # a and b are already joined: the tree path from a to b plus eid.
            # a == b only for an element between two boundaries (add_element
            # refuses an element whose ends are one node).
            came: dict[str, tuple[str, str]] = {a: (a, "")}
            queue = deque([a])
            while queue:
                n = queue.popleft()
                if n == b:
                    break
                for nb, via in tree.get(n, []):
                    if nb not in came:
                        came[nb] = (n, via)
                        queue.append(nb)
            nodes, edges, n = [b], [], b
            while n != a:
                n, via = came[n]
                nodes.append(n)
                edges.append(via)
            nodes.reverse()
            edges.reverse()
            # edges[i] joins nodes[i] and nodes[i + 1]; eid closes nodes[-1]
            # back to nodes[0].
            loop = [*edges, eid]
            if shared not in nodes:
                return loop, None
            i = nodes.index(shared)
            path = loop[i:] + loop[:i]
            first, last = self.elements[path[0]], self.elements[path[-1]]
            if len(path) == 1:
                ends = (first.from_node, first.to_node)
            else:
                ends = (
                    next(m for m in (first.from_node, first.to_node) if key(m) == shared),
                    next(m for m in (last.from_node, last.to_node) if key(m) == shared),
                )
            if ends[0] == ends[1]:
                # Out of and back into one boundary: a loop of its own.
                return loop, None
            order = list(self.nodes)
            if order.index(ends[0]) > order.index(ends[1]):
                # Name the boundaries in the order they were added.
                path.reverse()
                ends = (ends[1], ends[0])
            return path, ends
        return None

    def validate(self) -> None:
        """
        Validate the network topology.
        - Must have at least one node and one element.
        - Must have at least one PressureBoundary.
        - Must contain at least one element with losses (to prevent infinite flow).
        - No closed loop made only of lossless connections (#988), or more
          generally of elements whose Pt relation does not depend on their
          own flow (``_flow_independent_dpt``): no solution, or an
          undetermined split. Every PressureBoundary counts as one shared
          node, so a path of such elements between two boundaries is one
          such loop (#1003).
        - No isolated subgraphs (all nodes reachable from first node).
        - Every interior node must be connected to at least one element.
        """
        from .components import PressureBoundary

        if not self.nodes:
            raise ValueError("FlowNetwork validation failed: the network contains no nodes.")

        if len(self.elements) == 0:
            raise ValueError("FlowNetwork validation failed: the network contains no elements.")

        has_pressure_bc = any(isinstance(node, PressureBoundary) for node in self.nodes.values())
        if not has_pressure_bc:
            raise ValueError(
                "FlowNetwork validation failed: the network must contain at least "
                "one PressureBoundary to serve as a reference."
            )

        from .components import LosslessConnectionElement, PressureLossElement

        all_lossless = all(isinstance(e, LosslessConnectionElement) for e in self.elements.values())
        if all_lossless:
            raise ValueError(
                "Network contains only lossless connection elements. "
                "There must be at least one pressure drop element to solve flow."
            )

        # A closed loop of lossless connections (two in parallel, or a ring)
        # holds Pt equal all round it, so any flow circulating around it meets
        # every row: the steady problem has a one-parameter family of roots
        # and the solver returns whichever its start lies nearest (#988).
        # The same holds for any loop of elements whose Pt relation does not
        # depend on their own flow (a fraction loss without m_ref, a fixed
        # pressure ratio): the product of their Pt ratios round the loop is
        # fixed, so either it is not 1 and there is no root, or it is and the
        # circulation is free. A fraction loss's only "root" at a mismatch is
        # its zero-flow blend, where the forward and reverse forms cancel.
        found = self._flow_independent_loop()
        if found is not None:
            loop, ends = found
            members = [self.elements[eid] for eid in loop]
            if ends is not None:
                # A path between two boundaries closes the loop through their
                # fixed Pt (#1003): the GUI's plain edges to two outlets.
                described = ", ".join(
                    f"{e.id!r} ({e.from_node} -> {e.to_node}, {_flow_independent_dpt(e)})"
                    for e in members
                )
                raise ValueError(
                    f"FlowNetwork validation failed: the elements {described} join the "
                    f"pressure boundaries {ends[0]!r} and {ends[1]!r} by a path on which "
                    "every total-pressure relation is fixed independently of the flow. Each "
                    "boundary fixes its Pt, so the pressure ratio along the path is fixed at "
                    "both ends: unless it matches the boundaries' Pt ratio exactly there is no "
                    "solution, and if it does the flow along the path is not determined (an "
                    "outlet split can take any value, one outlet can even run backwards). Put "
                    "a flow-dependent loss (a duct, an orifice, a head loss) on each outlet "
                    "path, or merge the boundaries into one."
                )
            if all(isinstance(e, LosslessConnectionElement) for e in members):
                described = ", ".join(f"{e.id!r} ({e.from_node} -> {e.to_node})" for e in members)
                raise ValueError(
                    f"FlowNetwork validation failed: the lossless connections {described} "
                    "form a closed loop. Total pressure is equal all around it, so the flow "
                    "circulating around it is not determined: any value is a root, and the "
                    "solve would return an arbitrary one. Remove one of the connections, or "
                    "give one path a loss element (a real duct has one)."
                )
            described = ", ".join(
                f"{e.id!r} ({e.from_node} -> {e.to_node}, {_flow_independent_dpt(e)})"
                for e in members
            )
            raise ValueError(
                f"FlowNetwork validation failed: the elements {described} form a closed loop "
                "in which every total-pressure relation is fixed independently of the flow. "
                "Around such a loop the pressure ratios either fail to close (no solution: "
                "the solve can only sit where a fraction loss blends to zero flow, or fail) "
                "or close exactly (the split between the paths is not determined). Put a "
                "flow-dependent loss (an orifice, a channel, a head loss) on one path, or "
                "give each fraction loss on the loop an m_ref (set it wherever the flow may "
                "run low or reverse), or remove one of the elements."
            )

        # Per-element specific validation
        from .components import declared_row_kinds

        for _, elem in self.elements.items():
            elem.validate()
            # RESIDUAL-ROW-SCALING: the units of every residual row are
            # declared, where the residual is defined (raises TypeError).
            declared_row_kinds(elem)

        # A momentum chamber narrower than a pipe it adjoins (#1023): a
        # contraction or expansion with no loss. Warned, not refused; someone
        # may mean it.
        from .components import MomentumChamberNode

        for node in self.nodes.values():
            if isinstance(node, MomentumChamberNode):
                for message in node.pipe_area_warnings(self):
                    warnings.warn(message, UserWarning, stacklevel=3)

        # A wall exchanges heat between two convective surfaces. A side without
        # one (a PlenumNode, an orifice; a combustor or chamber with no surface
        # area) gave no heat and no error: the wall was silently inert. Checked
        # here, after topology, since a surface's area can be resolved from the
        # geometry.
        if self.thermal_coupling_enabled:
            for wid, wall in self.walls.items():
                for side in wall.endpoints:
                    obj = self.elements.get(side) or self.nodes.get(side)
                    if obj is not None and not obj.has_convective_surface:
                        raise ValueError(
                            f"FlowNetwork validation failed: wall {wid!r} side {side!r} "
                            f"({type(obj).__name__}) has no convective surface, so the "
                            "wall would exchange no heat. Wall a channel, an impingement "
                            "target, or a loss element, combustor or momentum chamber "
                            "with a surface area, or remove the wall."
                        )

        # PressureLossElement ambiguity: both endpoints have has_theta=True and no
        # explicit theta_source override.  Requires user disambiguation.
        for eid, elem in self.elements.items():
            if isinstance(elem, PressureLossElement) and elem._both_endpoints_theta:
                raise ValueError(
                    f"FlowNetwork validation failed: PressureLossElement '{eid}' has "
                    f"has_theta=True nodes on both endpoints "
                    f"('{elem.from_node}' and '{elem.to_node}'). "
                    f"Pass theta_source='<node_id>' to disambiguate."
                )

        # Isolated subgraph check: all nodes must be reachable from *some* PressureBoundary.
        # This prevents islands that have elements but no reference pressure.
        reachable_from_bc = set()
        for nid, node in self.nodes.items():
            if isinstance(node, PressureBoundary):
                reachable_from_bc.update(self._reachable_from(nid))

        isolated = set(self.nodes) - reachable_from_bc
        if isolated:
            raise ValueError(
                f"FlowNetwork validation failed: isolated node(s) detected: {sorted(isolated)}. "
                "Every node must be reachable from a PressureBoundary."
            )

        # A closed restriction (a zero bore, #1024) is a wall: its row is m = 0
        # and reads no pressure. A part of the network that closed elements
        # alone join to the rest keeps solving where it holds a
        # PressureBoundary of its own; one they cut off from every boundary
        # has no pressure level (its stagnant gas could be at any), so the
        # problem has no unique root: refused, naming the closed elements.
        open_from_bc: set[str] = set()
        for nid, node in self.nodes.items():
            if isinstance(node, PressureBoundary):
                open_from_bc.update(self._reachable_from(nid, through_closed=False))
        cut_off = [n for n in self.nodes if n not in open_from_bc]
        if cut_off:
            from .components import _is_closed_restriction

            closing = sorted(
                eid
                for eid, e in self.elements.items()
                if _is_closed_restriction(e) and (e.from_node in cut_off) != (e.to_node in cut_off)
            )
            raise ValueError(
                f"FlowNetwork validation failed: node(s) {cut_off} reach a "
                f"PressureBoundary only through closed restriction(s) {closing} (zero "
                "area: a wall, no flow). Nothing sets the pressure of the gas trapped "
                "beyond them, so the network has no unique solution. Open the "
                "restriction (give it an area), remove the part it closes off, or give "
                "that part a PressureBoundary of its own."
            )

        # Every interior node must have at least 2 connected elements to satisfy continuity (in and out)
        # Boundary nodes only need at least 1.
        for node_id, node in self.nodes.items():
            upstream = self._upstream_of_node[node_id]
            downstream = self._downstream_of_node[node_id]
            connected = upstream + downstream
            is_boundary = node.is_boundary

            if not connected:
                raise ValueError(
                    f"FlowNetwork validation failed: node '{node_id}' has no connected elements."
                )

            if not is_boundary and len(connected) < 2:
                raise ValueError(
                    f"FlowNetwork validation failed: interior node '{node_id}' has only {len(connected)} "
                    "connected element(s). Interior nodes must have at least 2 connections to satisfy "
                    "mass conservation."
                )

            # WallNode is a dead end - no flow may exit it.
            if type(node).__name__ == "WallNode" and downstream:
                raise ValueError(
                    f"FlowNetwork validation failed: WallNode '{node_id}' has downstream "
                    "connections. A wall is a closed end - connect only one element leading "
                    "into the wall, not out of it."
                )

            # Every interior node needs at least one element delivering flow in and one taking
            # flow out. All-upstream or all-downstream means mass cannot be conserved --
            # the most common cause is a tee junction wired with the wrong port as inlet.
            if not is_boundary and (not upstream or not downstream):
                direction = "upstream" if not upstream else "downstream"
                raise ValueError(
                    f"FlowNetwork validation failed: interior node '{node_id}' has no {direction} "
                    "connections -- flow has nowhere to go. This often means a tee junction port "
                    "is wired backwards (e.g. C on a branching tee should be the inlet)."
                )

            # MomentumChamberNode topology guard (#174). MCN's scalar closure
            # Pt = P + 0.5 rho v^2 only models a single bulk velocity, so it
            # cannot represent merging or splitting streams in the chamber
            # itself. >1 incoming or >1 outgoing edges puts MCN outside its
            # valid envelope; use MultiPortChamberBase (momentum-CV
            # junction with one MCN per port) for those topologies.
            #
            # A MERGE is the exception (#471): with ``main_inlet`` declared the
            # chamber closes several inflows by axial momentum (main face
            # offset, see MomentumChamberNode). A split stays refused -- one
            # momentum equation cannot set two outlet pressures.
            is_merge = getattr(node, "main_inlet", None) is not None
            if type(node).__name__ == "MomentumChamberNode" and (
                (len(upstream) > 1 and not is_merge) or len(downstream) > 1
            ):
                direction = "outgoing" if len(downstream) > 1 else "incoming"
                edges = downstream if len(downstream) > 1 else upstream
                hint = (
                    "A split needs a loss closure: use a tee / MultiPortChamberBase junction."
                    if direction == "outgoing"
                    else "Declare main_inlet to make it a merge chamber: the "
                    "other inflows become side streams closed by axial momentum."
                )
                raise ValueError(
                    f"FlowNetwork validation failed: MomentumChamberNode "
                    f"'{node_id}' has multiple {direction} edges ({sorted(edges)}). "
                    f"{hint}"
                )
            node_validate = getattr(node, "validate", None)
            if callable(node_validate):
                node_validate()
            # EnergyBoundary(Tt_target=...) refusals (#792)
            from .components import validate_energy_boundaries

            validate_energy_boundaries(node)

    def to_dict(self) -> dict:
        """
        Serialize the entire network topology and physics settings into a declarative dictionary suitable for JSON.
        Uses inspect.signature to dynamically capture constructor arguments from component instances.
        """
        nodes_data = {}
        for nid, node in self.nodes.items():
            sig = inspect.signature(node.__init__)
            # Extract arguments identically matching the __init__ signature keys
            args = {
                k: getattr(node, k, getattr(node, f"_{k}", None))
                for k in sig.parameters
                if k != "self" and k != "args" and k != "kwargs"
            }
            _restore_init_kwargs(node, args)
            nodes_data[nid] = {"type": type(node).__name__, "kwargs": args}
            # Attached after construction, so not among the arguments (#1099).
            ebs = getattr(node, "energy_boundaries", None)
            if ebs:
                nodes_data[nid]["energy_boundaries"] = [
                    {"id": eb.id, "Q": eb.Q, "fraction": eb.fraction, "Tt_target": eb.Tt_target}
                    for eb in ebs
                ]

        elements_data = {}
        for eid, element in self.elements.items():
            sig = inspect.signature(element.__init__)
            args = {
                k: getattr(element, k, getattr(element, f"_{k}", None))
                for k in sig.parameters
                if k != "self" and k != "args" and k != "kwargs"
            }
            _restore_init_kwargs(element, args)
            elements_data[eid] = {"type": type(element).__name__, "kwargs": args}

        walls_data = {}
        for wid, wall in self.walls.items():
            sig = inspect.signature(wall.__init__)
            args = {
                k: getattr(wall, k, getattr(wall, f"_{k}", None))
                for k in sig.parameters
                if k != "self" and k != "args" and k != "kwargs"
            }
            _restore_init_kwargs(wall, args)
            walls_data[wid] = {"type": type(wall).__name__, "kwargs": args}

        return {
            "nodes": nodes_data,
            "elements": elements_data,
            "walls": walls_data,
            "thermal_coupling_enabled": self.thermal_coupling_enabled,
        }

    @classmethod
    def from_dict(cls, data: dict) -> Self:
        """
        Deserialize a flow network from a declarative dictionary generated by to_dict().
        """
        from .components import EnergyBoundary

        network = cls()

        # Restore thermal coupling toggle
        network.thermal_coupling_enabled = data.get("thermal_coupling_enabled", True)

        for _nid, n_data in data.get("nodes", {}).items():
            NodeClass = _component_class(n_data["type"])
            node = NodeClass(**n_data["kwargs"])
            for eb_kwargs in n_data.get("energy_boundaries", ()):
                node.add_energy_boundary(EnergyBoundary(**eb_kwargs))
            network.add_node(node)

        for _eid, e_data in data.get("elements", {}).items():
            ElemClass = _component_class(e_data["type"])
            network.add_element(ElemClass(**e_data["kwargs"]))

        for _wid, w_data in data.get("walls", {}).items():
            WallClass = _component_class(w_data["type"])
            network.add_wall(WallClass(**w_data["kwargs"]))

        return network


def _restore_init_kwargs(obj: Any, args: dict) -> None:
    """Replace serialised arguments whose attribute is not the constructor's value.

    ``to_dict`` reads each ``__init__`` argument back from the attribute of the
    same name. A class whose attribute differs from the argument (rewritten
    by ``resolve_topology``, say, or derived) declares
    ``_init_kwarg_attrs = {argument: attribute name or function of the
    object}``; the maps of its base classes apply too, the class's own
    winning (#1099).
    """
    aliases: dict = {}
    for klass in reversed(type(obj).__mro__):
        aliases.update(vars(klass).get("_init_kwarg_attrs", {}))
    for k, source in aliases.items():
        if k in args:
            args[k] = source(obj) if callable(source) else getattr(obj, source)


def _component_class(name: str) -> type:
    """The network class ``to_dict`` recorded as ``type``: in ``components``,
    the package's exports (the junction elements live in their own modules)
    or ``ejector_element`` (not exported) (#1099)."""
    import importlib

    for mod in ("components", "", "ejector_element"):
        module = importlib.import_module(f"{__package__}.{mod}" if mod else __package__)
        cls = getattr(module, name, None)
        if isinstance(cls, type):
            return cls
    raise ValueError(f"FlowNetwork.from_dict: unknown component type {name!r}")


def _flow_independent_dpt(element: NetworkElement) -> str | None:
    """Why this element's total-pressure relation does not depend on its own
    flow, or None if it does (or cannot be told).

    These are the edges of ``FlowNetwork._flow_independent_loop``: a closed
    loop made only of them has no solution or an undetermined split.

    - A pass-through (#1023): an ``AreaChangeElement`` with equal faces, an
      ``OrificeElement`` whose bore equals its pipe (#788): Pt_out = Pt_in.
      Known once the topology is resolved.
    - ``LosslessConnectionElement``, ``VortexElement`` (a free vortex
      conserves Pt, #661) and a flow splitter's reference outlet (its valve
      open): Pt_out = Pt_in.
    - ``CompressorElement`` / ``TurbineElement``: Pt_out = PR^(+-1) Pt_in at a
      fixed pressure ratio (#893).
    - A ``PressureLossElement`` with ``ConstantFractionLoss`` or
      ``LinearThetaFractionLoss`` and ``m_ref`` None: Pt_out = (1 - xi) Pt_in
      at every flow in either direction (``k theta + xi0`` for the linear
      form: theta is a temperature ratio at its source node, not this
      element's flow). With m_ref set the loss fades to zero at zero flow, so
      it is flow-dependent, unless its xi is zero at every flow.

    Not included: head losses (they scale with the flow's dynamic head), user
    callables (their dependence is not known before a solve), metered
    splitter outlets (their row fixes the flow and leaves Pt free) and every
    other element, whose drop depends on its flow.
    """
    from .components import (
        AreaChangeElement,
        LosslessConnectionElement,
        OrificeElement,
        PressureLossElement,
        SplitterOutletElement,
        TurbomachineElement,
        VortexElement,
    )
    from .pressure_loss import ConstantFractionLoss, LinearThetaFractionLoss

    if isinstance(element, LosslessConnectionElement):
        return "lossless"
    if isinstance(element, (AreaChangeElement, OrificeElement)) and element._pass_through:
        return "pass-through (area equal to its neighbour's), Pt conserved"
    if isinstance(element, VortexElement):
        return "free vortex, Pt conserved"
    if isinstance(element, SplitterOutletElement) and element.reference:
        return "splitter reference outlet, Pt conserved"
    if isinstance(element, TurbomachineElement):
        return f"{type(element).__name__} at fixed pressure ratio {element.pressure_ratio:g}"
    if isinstance(element, PressureLossElement):
        c = element.correlation
        if isinstance(c, ConstantFractionLoss):
            zero = c.xi == 0.0
        elif isinstance(c, LinearThetaFractionLoss):
            zero = c.k == 0.0 and c.xi0 == 0.0
        else:
            return None
        if zero:
            return f"{type(c).__name__} with zero loss"
        if c.m_ref is None:
            return f"{type(c).__name__} without m_ref"
    return None
