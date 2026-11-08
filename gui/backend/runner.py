"""Standalone network runner for programmatic and parametric use.

Loads a saved GUI network JSON, applies boundary-condition overrides
keyed by node label, solves, and returns a NetworkResult that exposes
a full-detail DataFrame, raw state access, and a sweep helper for
parametric studies.

Override key format: ``"<label>.<attribute>"``, where ``<label>`` is the
node's GUI label and ``<attribute>`` is the data field to override
(e.g. ``"air_inlet.m_dot"``, ``"fuel_inlet.Tt"``, ``"outlet.Pt"``).

Typical notebook usage::

    from gui.backend.runner import NetworkRunner

    runner = NetworkRunner.from_file("my_model.json")
    result = runner.solve({"air_inlet.m_dot": 1.2})
    print(result.success, result.final_norm)
    df = result.to_dataframe()

    import pandas as pd
    params = pd.DataFrame({"air_inlet.m_dot": [0.8, 1.0, 1.2]})
    sweep_df = runner.sweep(params, metrics=["channel_1.m_dot"])
    sweep_df["outputs"]["channel_1", "m_dot"]   # one row per case (#784)

The sweep's table, and ``load_results`` / ``append_case`` for exported CSVs
and single results, are in ``gui.backend.results``.
"""

import json
import math
import warnings
from pathlib import Path
from typing import Any

import pandas as pd

import combaero as cb

from .graph_builder import build_network_from_schema
from .results import input_cells
from .schemas import (
    DEFAULT_SOLVE_TIMEOUT,
    ElementResult,
    NetworkGraphSchema,
    NodeResult,
    StateResult,
)
from .second_law_view import element_fields
from .solver_setup import make_solver
from .units import label_with_unit
from .validity_view import element_fields as validity_fields

# ---------------------------------------------------------------------------
# Shared result-building helper (used by both the HTTP API and NetworkRunner)
# ---------------------------------------------------------------------------


def _schema_maps(
    schema: NetworkGraphSchema,
) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    """Build (id_to_label, id_to_kind, id_to_method) dicts from a schema.

    Shared by NetworkRunner.__init__ and the HTTP export endpoint so label
    resolution is consistent in both paths.
    """
    id_to_label: dict[str, str] = {}
    id_to_kind: dict[str, str] = {}
    id_to_method: dict[str, str] = {}
    for n in schema.nodes:
        lbl = n.data.get("label") or ""
        id_to_label[n.id] = lbl
        id_to_kind[n.id] = n.type
        if n.type == "combustor":
            id_to_method[n.id] = n.data.get("method", "complete")
    for e in schema.edges:
        lbl = (e.data or {}).get("label") or ""
        id_to_label[e.id] = lbl
        id_to_kind[e.id] = (e.data or {}).get("type") or "flow"
        auto_id = f"__auto_link__{e.source}__{e.target}"
        id_to_label[auto_id] = lbl
        id_to_kind[auto_id] = "lossless_connection"
    return id_to_label, id_to_kind, id_to_method


# Builder-internal ids (#783): what each part of an array or liner is, by its
# suffix, per parent kind (impingement_array.py, effusion_liner.py).
_ARRAY_ROLES: dict[str, dict[str, str]] = {
    "impingement_array": {"p": "plate row", "x": "crossflow segment", "c": "crossflow node"},
    "effusion_liner": {"b": "duct station", "s": "duct segment", "p": "panel"},
}


def _synthetic_meta(
    item_id: str, schema: NetworkGraphSchema, id_to_label: dict[str, str]
) -> tuple[str, str, str] | None:
    """(kind, label, parent) of an id the builder made, not one the user drew
    (#783): the junction between two drawn elements, a junction's port
    chamber, an auto link, or a row / panel / segment of an array. None for a
    drawn id or an unknown pattern. The parent is the drawn node or edge it
    belongs to; the label is derived from the parent's."""

    def name(i: str) -> str:
        return id_to_label.get(i) or i

    def edge_between(a: str, b: str) -> str:
        return next((e.id for e in schema.edges if e.source == a and e.target == b), "")

    if item_id.startswith("__junction__"):
        a, _, b = item_id[len("__junction__") :].partition("__")
        return "junction", f"junction {name(a)} -> {name(b)}", edge_between(a, b)
    if item_id.startswith("__auto_link__"):
        a, _, b = item_id[len("__auto_link__") :].partition("__")
        edge = edge_between(a, b)
        label = id_to_label.get(edge) or f"link {name(a)} -> {name(b)}"
        return "lossless_connection", label, edge
    if item_id.startswith("__tee_jct__"):
        rest = item_id[len("__tee_jct__") :]
        drawn = sorted((n.id for n in schema.nodes if rest.startswith(n.id + "_")), key=len)
        if drawn:
            jct = drawn[-1]
            port = rest[len(jct) + 1 :]
            return "port_chamber", f"{name(jct)} {port} port", jct
        return None
    parent, sep, part = item_id.rpartition("__")
    if sep and part[:1] and part[1:].isdigit():
        kind = next((n.type for n in schema.nodes if n.id == parent), None)
        role = _ARRAY_ROLES.get(kind or "", {}).get(part[:1])
        if role:
            return role.replace(" ", "_"), f"{name(parent)} {role} {part[1:]}", parent
    return None


def _build_result_objects(
    result: dict,
    net,
) -> tuple[dict[str, NodeResult], dict[str, ElementResult], dict[str, dict]]:
    """Extract typed result objects from the raw solver output dict.

    Factored out of the HTTP solve handler so that both the FastAPI layer
    and the standalone NetworkRunner share identical result construction.
    """
    success = bool(result.get("__success__", False))

    node_results: dict[str, NodeResult] = {}
    for node_id in net.nodes:
        prefix = f"{node_id}."
        node_vals = {k[len(prefix) :]: v for k, v in result.items() if k.startswith(prefix)}

        y_vals: list[float] = []
        x_vals: list[float] = []
        i = 0
        while f"Y[{i}]" in node_vals:
            y_vals.append(float(node_vals.pop(f"Y[{i}]")))
            if f"X[{i}]" in node_vals:
                x_vals.append(float(node_vals.pop(f"X[{i}]")))
            i += 1

        # A pressure boundary emits only Pt/Tt: it is a reservoir, so its
        # static state IS its total state (its rho, h, ... are evaluated there).
        # The old fallback was a made-up 101325 Pa / 300 K, which the probe and
        # the CSV export showed beside the real Pt.
        Pt = node_vals.pop("Pt", None)
        Tt = node_vals.pop("Tt", None)
        state_res = StateResult(
            T=float(node_vals.pop("T", math.nan if Tt is None else Tt)),
            P=float(node_vals.pop("P", math.nan if Pt is None else Pt)),
            Pt=Pt,
            Tt=Tt,
            Y=y_vals,
            X=x_vals or None,
            **node_vals,
        )
        node_results[node_id] = NodeResult(state=state_res, success=success)

    # A heat source (#792) reports the heat it applied -- the SOLVED Q in a
    # Tt_target mode, the set one otherwise, NaN when undetermined (no
    # through-flow) -- beside the state of the node it heats, whose Tt is
    # that node's TOTAL temperature (T its static one, #782).
    for hs_id, (target, eb) in getattr(net, "_gui_heat_sources", {}).items():
        node_res = node_results.get(target)
        if node_res is None:
            continue
        extra = {"Q": float(result.get(f"{hs_id}.Q", math.nan))}
        if eb.is_target:
            extra["Tt_target"] = float(eb.Tt_target)
        node_results[hs_id] = NodeResult(
            state=node_res.state, success=success, heated_node=target, **extra
        )

    element_results: dict[str, ElementResult] = {}
    elem_diags = result.get("__element_diag__", {})
    for elem_id in net.elements:
        diag = elem_diags.get(elem_id, {}).copy()
        # Junction elements have no own `.m_dot` unknown (flows live on the
        # connecting elements). Report a representative through-flow from their
        # diagnostics: the tee's common leg (`m_dot_com`), or the ejector's
        # total outlet flow (`m_dot_outlet`).
        m_dot = float(
            result.get(
                f"{elem_id}.m_dot",
                diag.get("m_dot_com", diag.get("m_dot_outlet", 0.0)),
            )
        )
        diag.pop("m_dot", None)
        diag.update(element_fields((result.get("__second_law__") or {}).get(elem_id)))
        diag.update(validity_fields((result.get("__incompressible_validity__") or {}).get(elem_id)))
        element_results[elem_id] = ElementResult(m_dot=m_dot, success=success, **diag)

    # A tee port's hidden sharp step (#536) is reported on the tee the user
    # drew: the step's loss coefficient, its total-pressure drop and the
    # channel-to-port area ratio, per port.
    for step_id, (tee_id, port, ratio) in getattr(net, "_gui_tee_steps", {}).items():
        step = elem_diags.get(step_id)
        tee_res = element_results.get(tee_id)
        if not step or tee_res is None:
            continue
        extra = {
            f"step_{port}_zeta": float(step.get("zeta", 0.0)),
            f"step_{port}_dPt": float(step.get("dPt", 0.0)),
            f"step_{port}_connected_area_ratio": float(ratio),
        }
        element_results[tee_id] = ElementResult(**{**tee_res.model_dump(), **extra})

    def _safe_float(v: Any) -> Any:
        try:
            if isinstance(v, (list, tuple, dict)):
                return v
            return float(v)
        except (TypeError, ValueError):
            return v

    edge_results: dict[str, dict] = {}
    for edge_id in net.walls:
        edge_keys = {
            k.split(".", 1)[1]: _safe_float(v)
            for k, v in result.items()
            if k.startswith(f"{edge_id}.")
        }
        if edge_keys:
            edge_results[edge_id] = edge_keys

    # Impingement arrays (#465) are one GUI node but many network elements:
    # report the array under its own id, and each wall to it as one edge.
    for arr_id, arr in getattr(net, "_gui_arrays", {}).items():
        if not elem_diags:
            continue
        summary = arr.summarize(elem_diags)
        m_total = summary.pop("m_dot")
        element_results[arr_id] = ElementResult(m_dot=m_total, success=success, **summary)
        if not hasattr(arr, "plate_id"):
            continue  # an effusion liner owns its walls: no wall edge to aggregate
        plates = {arr.plate_id(i) for i in range(1, int(arr.n_rows) + 1)}
        rows: dict[str, list[dict]] = {}
        for wid, wall in net.walls.items():
            if wid in edge_results and {wall.element_a, wall.element_b} & plates:
                # ImpingementArray.add_wall names row walls f"{edge_id}__r{i}".
                rows.setdefault(wid.rpartition("__r")[0], []).append(edge_results[wid])
        for base, per_row in rows.items():
            edge_results[base] = _aggregate_array_wall(per_row)
        # The per-row walls are internal ids (no drawn edge); their values live
        # on in the aggregate's rows_* lists. Kept, they doubled the array's
        # heat in any sum over the export's wall rows.
        for wid in [w for w in edge_results if w.rpartition("__r")[0] in rows]:
            del edge_results[wid]

    return node_results, element_results, edge_results


def _aggregate_array_wall(per_row: list[dict]) -> dict:
    """One GUI wall spread over an array's rows: the total heat, and the
    profile of the HOTTEST row (the one a designer checks the wall against).
    The per-row values stay available as ``rows_Q`` and ``rows_T_hot``."""
    hottest = max(per_row, key=lambda r: float(r.get("T_hot", 0.0)))
    out = dict(hottest)
    out["Q"] = float(sum(float(r.get("Q", 0.0)) for r in per_row))
    out["rows_Q"] = [float(r.get("Q", 0.0)) for r in per_row]
    out["rows_T_hot"] = [float(r.get("T_hot", 0.0)) for r in per_row]
    return out


# ---------------------------------------------------------------------------
# NetworkResult
# ---------------------------------------------------------------------------


class NetworkResult:
    """Result of a single network solve.

    Carries the full solved state for programmatic access, DataFrame export,
    and use as a warm start for follow-on operations such as boundary swaps.

    Attributes:
        success: True when the solve both converged AND passed every element's
            physical-consistency check. It is NOT "converged within the
            residual tolerance" -- a root that a junction rejects as
            unphysical reports ``success=False`` with a small ``final_norm``,
            which read as a contradiction until ``converged`` and
            ``consistent`` were reported separately.
        converged: True when the root finder reached the residual tolerance,
            regardless of what the consistency checks then said.
        consistent: True/False from the elements' own physical-consistency
            checks, or None when nothing checked -- either because the solve
            never converged, or because no element in this network has a
            verifier. None is NOT a pass.
        inconsistent_elements: ids of the elements that rejected the solution.
        outcome: machine-readable reason the solve ended
            (``combaero.network.SolveOutcome``). Prefer this over matching on
            ``message``, whose wording comes partly from SciPy.
        message: Human-readable solver status string.
        final_norm: Euclidean norm of the residual vector at convergence.
        worst_residuals: the rows carrying the residual, largest first.
    """

    def __init__(
        self,
        raw: dict,
        node_results: dict[str, NodeResult],
        element_results: dict[str, ElementResult],
        edge_results: dict[str, dict],
        net: Any,
        schema: NetworkGraphSchema,
        id_to_label: dict[str, str],
        id_to_kind: dict[str, str],
        id_to_method: dict[str, str],
    ) -> None:
        self.success: bool = bool(raw.get("__success__", False))
        self.converged: bool = bool(raw.get("__converged__", False))
        self.consistent: bool | None = raw.get("__consistent__")
        self.inconsistent_elements: list[str] = list(raw.get("__inconsistent_elements__") or [])
        self.outcome: str = str(raw.get("__outcome__", "not_converged"))
        self.message: str = str(raw.get("__message__", ""))
        self.final_norm: float | None = raw.get("__final_norm__")
        self.worst_residuals: list[dict] = list(raw.get("__worst_residuals__") or [])
        self._x_solution: Any = raw.get("__x_solution__")
        unk_names: list[str] = raw.get("__unknown_names__", [])
        x_sol: list[float] = raw.get("__x_solution__") or []
        self._solved_state: dict[str, float] = dict(zip(unk_names, x_sol, strict=False))
        self._raw = raw
        self._node_results = node_results
        self._element_results = element_results
        self._edge_results = edge_results
        self._net = net
        self._schema = schema
        self._id_to_label = id_to_label
        self._id_to_kind = id_to_kind
        self._id_to_method = id_to_method
        # How the case was solved, for the results table's attrs (#784): not
        # inputs, since they do not change a converged answer. The GUI routes
        # solve with the schema's settings; NetworkRunner.solve overwrites
        # these with the arguments it was given.
        settings = schema.solver_settings
        self._solver_options: dict[str, Any] = {
            "method": settings.method,
            "init_strategy": settings.init_strategy,
            "timeout": settings.timeout,
            "stall_handover": settings.stall_handover,
            "incompressible_mach_limit": settings.incompressible_mach_limit,
        }
        self._warm_start: Any = None

    def get(self, key: str) -> float:
        """Get a scalar result by solver key or ``'<label>.<quantity>'``.

        Tries the key directly first (raw ``'<id>.<quantity>'`` form), then
        resolves the left-hand part as a node/element label so callers can
        write ``result.get("combustor.Tt")`` instead of a UUID key.
        """
        if key in self._raw:
            return float(self._raw[key])
        if "." in key:
            label, quantity = key.split(".", 1)
            try:
                node_id = self._resolve_label(label)
                raw_key = f"{node_id}.{quantity}"
                if raw_key in self._raw:
                    return float(self._raw[raw_key])
            except KeyError:
                pass
        raise KeyError(
            f"Result key '{key}' not found. "
            "Keys follow the format '<id>.<quantity>' or '<label>.<quantity>'."
        )

    def node_state(self, label: str) -> dict:
        """Return the full solved thermodynamic state for a node by label or id.

        The returned dict contains T, P, Pt, Tt, rho, mach, Y, X and any
        additional quantities exposed by that node type.
        """
        node_id = self._resolve_label(label)
        res = self._node_results.get(node_id)
        if res is None:
            raise KeyError(f"No node with label or id '{label}' in results.")
        return res.state.model_dump()

    def _resolve_label(self, label: str) -> str:
        """Return the node/element id for a given label, or the label itself if it is an id."""
        for nid, lbl in self._id_to_label.items():
            if lbl == label:
                return nid
        if label in self._node_results or label in self._element_results:
            return label
        raise KeyError(f"No node or element with label or id '{label}'.")

    def to_dataframe(self) -> pd.DataFrame:
        """Full-detail DataFrame equivalent to the GUI CSV export.

        Returns one row per network entity (nodes, flow elements, thermal walls)
        with all available diagnostic quantities.  Column names carry unit
        annotations matching the GUI export (e.g. ``'P [Pa]'``, ``'T [K]'``).
        Every row carries the solve's status (``success``, ``converged``,
        ``consistent``, ``outcome``, ``final_norm``): check ``success`` before
        trusting the numbers.

        After the result columns come the case's inputs (#784, 0.8.0):
        ``input:<key> [<unit>]`` on each entity's row, ``initial_guess:<key>``,
        and a ``network`` row (id ``_network``) with the network-wide
        settings and the ``solver:<key>`` options. An edge with inputs and no
        result row gets an ``edge`` row. ``gui.backend.results.load_results``
        reads this frame, or its CSV, into the one-row-per-case table.
        """
        species_labels: list[str] = list(cb.species.names)
        data: list[dict] = []

        for node_id, res in self._node_results.items():
            state_data = res.state.model_dump()
            Y = state_data.pop("Y", [])
            X = state_data.pop("X", None)
            kind = self._id_to_kind.get(node_id, "")
            row: dict = {
                "type": "node",
                "kind": kind,
                "is_boundary": bool(getattr(self._net.nodes.get(node_id), "is_boundary", False)),
                "combustion_method": self._id_to_method.get(node_id, ""),
                "id": node_id,
                "label": self._id_to_label.get(node_id, ""),
                **state_data,
            }
            for i, y_val in enumerate(Y):
                name = species_labels[i] if i < len(species_labels) else str(i)
                row[f"Y[{name}]"] = y_val
            if X:
                for i, x_val in enumerate(X):
                    name = species_labels[i] if i < len(species_labels) else str(i)
                    row[f"X[{name}]"] = x_val
            data.append(row)

        for elem_id, res in self._element_results.items():
            elem_data = res.model_dump()
            data.append(
                {
                    "type": "element",
                    "kind": self._id_to_kind.get(elem_id, ""),
                    "is_boundary": False,
                    "combustion_method": "",
                    "id": elem_id,
                    "label": self._id_to_label.get(elem_id, ""),
                    **elem_data,
                }
            )

        for wall_id, res in self._edge_results.items():
            data.append(
                {
                    "type": "thermal_wall",
                    "kind": "thermal_wall",
                    "is_boundary": False,
                    "combustion_method": "",
                    "id": wall_id,
                    "label": self._id_to_label.get(wall_id, ""),
                    **res,
                }
            )

        # The case's inputs (#784), after every result column: so an exported
        # CSV reads back into the full results table (gui.backend.results).
        cells, input_rows = input_cells(
            self._schema, [r["id"] for r in data], self._solver_options, self._warm_start
        )
        data.extend(input_rows)

        # Builder-internal rows (#783) say what they are and which drawn item
        # they belong to; a drawn row's parent is empty.
        for row in data:
            meta = _synthetic_meta(str(row["id"]), self._schema, self._id_to_label)
            if meta is not None:
                kind, label, parent = meta
                row["kind"] = row.get("kind") or kind
                row["label"] = row.get("label") or label
                row["parent"] = parent
            else:
                row.setdefault("parent", "")

        df = pd.DataFrame(data)
        if "parent" in df.columns and "label" in df.columns:
            cols = list(df.columns)
            cols.insert(cols.index("label") + 1, cols.pop(cols.index("parent")))
            df = df[cols]
        # The solve's own status on EVERY row, beside the label: a timed-out or
        # rejected solve exports numbers that look valid, and a filtered or
        # pivoted frame must still say which solve its rows came from. Element
        # rows used to carry `success` alone and node rows nothing.
        status = {
            "success": self.success,
            "converged": self.converged,
            "consistent": self.consistent,
            "outcome": self.outcome,
            "final_norm": self.final_norm,
        }
        df = df.drop(columns=[c for c in status if c in df.columns])
        at = list(df.columns).index("label") + 1 if "label" in df.columns else 0
        for offset, (col, value) in enumerate(status.items()):
            df.insert(at + offset, col, [value] * len(df))
        df = df.rename(columns={col: label_with_unit(col) for col in df.columns})
        inputs = pd.DataFrame([cells.get(rid, {}) for rid in df["id"]], index=df.index)
        return pd.concat([df, inputs], axis=1)

    def swap_boundary(self, label: str, new_type: str) -> "NetworkRunner":
        """Return a new NetworkRunner with one boundary node retyped.

        The required field for ``new_type`` is auto-populated from this solved
        state so that the first re-solve reproduces the same physical operating
        point:

        * ``mass_boundary`` -> ``pressure_boundary``: copies solved ``Pt``.
        * ``pressure_boundary`` -> ``mass_boundary``: sums ``m_dot`` of all
          adjacent elements.

        ``Tt``, ``composition``, ``label``, and ``initial_guess`` are
        preserved unchanged.

        Args:
            label: Node label or raw node id.
            new_type: ``"mass_boundary"`` or ``"pressure_boundary"``.

        Returns:
            A new :class:`NetworkRunner` ready for ``.solve()`` or ``.sweep()``.

        Raises:
            ValueError: If ``new_type`` is not supported, the node is already
                that type, or the existing type cannot be swapped.
            KeyError: If ``label`` does not resolve to a node.
        """
        _SUPPORTED = {"mass_boundary", "pressure_boundary"}
        if new_type not in _SUPPORTED:
            raise ValueError(f"swap_boundary only supports {_SUPPORTED!r}, got '{new_type}'.")
        node_id = self._resolve_label(label)
        schema = self._schema.model_copy(deep=True)
        node = next((n for n in schema.nodes if n.id == node_id), None)
        if node is None:
            raise KeyError(f"Node '{label}' not found in schema.")
        old_type = node.type
        if old_type == new_type:
            raise ValueError(f"Node '{label}' is already '{new_type}'.")
        if old_type not in _SUPPORTED:
            raise ValueError(
                f"Node '{label}' has type '{old_type}', which cannot be swapped. "
                f"Only {_SUPPORTED!r} are supported."
            )
        state = self._node_results[node_id].state
        if new_type == "pressure_boundary":
            if state.Pt is None:
                raise ValueError(
                    f"Solved state for '{label}' has no Pt -- cannot swap to pressure_boundary."
                )
            node.data.pop("m_dot", None)
            node.data["Pt"] = float(state.Pt)
        else:
            node.data.pop("Pt", None)
            node.data["m_dot"] = self._adjacent_m_dot(node_id)
        node.type = new_type

        # Seed initial_guess on internal (non-pressure-boundary) nodes so that
        # the first solve after the swap starts from a physically valid point.
        # PressureBoundary nodes must NOT have P/Pt/T/Tt in initial_guess --
        # _get_node_state() prefers initial_guess over node.Pt/Tt, so seeding
        # them would prevent BC overrides (e.g. pressure sweeps) from taking effect.
        for n in schema.nodes:
            if n.type == "pressure_boundary":
                continue
            res = self._node_results.get(n.id)
            if res is None:
                continue
            s = res.state
            guess: dict[str, float] = {}
            if s.T is not None:
                guess["T"] = float(s.T)
            if s.P is not None:
                guess["P"] = float(s.P)
            if s.Pt is not None:
                guess["Pt"] = float(s.Pt)
            if s.Tt is not None:
                guess["Tt"] = float(s.Tt)
            if guess:
                n.data["initial_guess"] = guess

        for e in schema.edges:
            eres = self._element_results.get(e.id)
            if eres is None:
                continue
            if e.data is None:
                e.data = {}
            e.data["initial_guess"] = {"m_dot": float(eres.m_dot)}

        return NetworkRunner(schema, _warm_state=self._solved_state)

    def _adjacent_m_dot(self, node_id: str) -> float:
        """Sum m_dot of all elements directly connected to a boundary node."""
        downstream = self._net.get_downstream_elements(node_id)
        upstream = self._net.get_upstream_elements(node_id)
        connected = downstream or upstream
        if not connected:
            raise ValueError(
                f"Node '{node_id}' has no connected elements -- cannot determine m_dot."
            )
        return sum(
            float(self._element_results[e.id].m_dot)
            for e in connected
            if e.id in self._element_results
        )


# ---------------------------------------------------------------------------
# NetworkRunner
# ---------------------------------------------------------------------------


class NetworkRunner:
    """Programmatic driver for a saved combaero-gui network.

    Load a network JSON file saved from the GUI, apply boundary-condition
    overrides by node label, and run parametric sweeps without starting
    the web server.

    Override key format: ``"<label>.<attribute>"`` where ``<label>`` is the
    node's GUI label and ``<attribute>`` is the data field name
    (e.g. ``"air_inlet.m_dot"``, ``"fuel_inlet.Tt"``, ``"outlet.Pt"``).
    The label falls back to the node id when no label is set.

    Example::

        runner = NetworkRunner.from_file("network.json")

        # Single solve with BC overrides
        result = runner.solve({"air_inlet.m_dot": 1.2, "fuel_inlet.Tt": 850.0})
        print(result.success, result.final_norm)

        # Parametric sweep
        import pandas as pd
        params = pd.DataFrame({
            "air_inlet.m_dot": [0.8, 1.0, 1.2],
            "fuel_inlet.m_dot": [0.03, 0.04, 0.05],
        })
        sweep_df = runner.sweep(params, metrics=["combustor.Tt", "hot_channel.m_dot"])
        T_exit = sweep_df["outputs"]["combustor", "Tt"]
    """

    def __init__(
        self,
        schema: NetworkGraphSchema,
        _warm_state: dict[str, float] | None = None,
    ) -> None:
        self._schema = schema
        self._id_to_label, self._id_to_kind, self._id_to_method = _schema_maps(schema)
        self._label_to_node_id: dict[str, str] = {
            lbl: nid for nid, lbl in self._id_to_label.items() if lbl
        }
        self._warm_state = _warm_state or {}

    @classmethod
    def from_file(cls, path: str | Path) -> "NetworkRunner":
        """Load a GUI-saved network JSON file and return a runner."""
        with Path(path).open() as fh:
            return cls.from_dict(json.load(fh))

    @classmethod
    def from_dict(cls, d: dict) -> "NetworkRunner":
        """Construct a runner from an already-loaded dict."""
        return cls(NetworkGraphSchema.model_validate(d))

    # ------------------------------------------------------------------
    # Solve
    # ------------------------------------------------------------------

    def solve(
        self,
        overrides: dict[str, Any] | None = None,
        method: str = "hybr",
        init_strategy: str = "default",
        timeout: float | None = DEFAULT_SOLVE_TIMEOUT,
        _x0: Any = None,
    ) -> NetworkResult:
        """Solve the network, optionally patching boundary conditions.

        Each call is stateless: a deep copy of the schema is made before
        applying overrides, so repeated calls (e.g. inside an optimisation
        loop) do not interfere with each other.

        Args:
            overrides: Dict of ``"<label>.<attribute>": value`` pairs.
                Examples: ``{"air_inlet.m_dot": 1.2}``,
                ``{"outlet.Pt": 200000.0}``.
            method: scipy root-finding method passed to NetworkSolver.
            init_strategy: Solver initialisation strategy
                (``"default"``, ``"analytical_pt_prop"``,
                ``"outletref_warmstart"``, ``"homotopy"``).
            timeout: Per-solve soft timeout in seconds (default
                ``DEFAULT_SOLVE_TIMEOUT``, 30 s, as in the GUI; ``None``: no limit).
            _x0: Internal -- solver solution vector from a prior solve used
                as the initial guess (passed by ``sweep`` for continuation).

        Returns:
            NetworkResult with solved state and DataFrame export.
        """
        import numpy as np

        schema = self._schema.model_copy(deep=True)
        if overrides:
            self._apply_overrides(schema, overrides)

        net = build_network_from_schema(schema)
        solver = make_solver(net, schema)
        x0_array: np.ndarray | None = None
        if _x0 is not None:
            x0_array = np.asarray(_x0)
        elif self._warm_state:
            default_x0 = solver._build_x0()
            x0_array = np.array(
                [
                    self._warm_state.get(name, default_x0[i])
                    for i, name in enumerate(solver.unknown_names)
                ],
                dtype=float,
            )
        raw = solver.solve(
            method=method,
            init_strategy=init_strategy,
            timeout=timeout,
            x0=x0_array,
            stall_handover=schema.solver_settings.stall_handover,
        )
        node_results, element_results, edge_results = _build_result_objects(raw, net)

        result = NetworkResult(
            raw=raw,
            node_results=node_results,
            element_results=element_results,
            edge_results=edge_results,
            net=net,
            schema=schema,
            id_to_label=dict(self._id_to_label),
            id_to_kind=dict(self._id_to_kind),
            id_to_method=dict(self._id_to_method),
        )
        result._solver_options.update(
            {"method": method, "init_strategy": init_strategy, "timeout": timeout}
        )
        if _x0 is not None:
            result._warm_start = "x0"
        elif self._warm_state:
            result._warm_start = "solved_state"
        return result

    # ------------------------------------------------------------------
    # Sweep
    # ------------------------------------------------------------------

    def sweep(
        self,
        params: pd.DataFrame,
        metrics: list[str] | None = None,
        method: str = "hybr",
        init_strategy: str = "default",
        timeout: float | None = DEFAULT_SOLVE_TIMEOUT,
    ) -> pd.DataFrame:
        """Run a parametric sweep over rows of ``params``.

        Each row of ``params`` is a set of BC overrides (column names are
        override keys in ``<label>.<attribute>`` format) and one case of the
        result.

        Each row first tries the previous converged row's solution as its
        initial guess; when that does not converge, the row is solved cold,
        as ``solve()`` would. A row's ``warm_start`` names the case it started
        from, if any.

        Args:
            params: DataFrame whose columns are override keys and whose rows
                are individual parameter combinations. Its index labels the
                cases (``0, 1, 2...`` for a default index; set it to names
                for readable rows).
            metrics: Optional list of ``"<label>.<quantity>"`` keys (or
                ``"<id>.<quantity>"``) to keep in ``outputs``, e.g.
                ``["combustor.Tt", "hot_channel.m_dot"]``. ``None`` keeps every
                quantity of every entity. A key the solve did not produce is a
                NaN column.
            method: scipy root-finding method forwarded to the solver.
            init_strategy: Solver initialisation strategy.
            timeout: Per-solve timeout in seconds (default
                ``DEFAULT_SOLVE_TIMEOUT``, 30 s, as in the GUI; ``None``: no limit).

        Returns:
            The results table of ``gui.backend.results`` (#784): one row per
            case, ``(kind, entity, quantity)`` columns. ``df["inputs"]`` is
            each case's full network inputs (the overridden values included),
            ``df["outputs"]`` its results with the solve status under
            ``_solve``; units, solver options, initial guesses and the warm
            start are in ``df.attrs``. A row whose solve raised has its
            inputs, ``_solve.success`` False and ``_solve.outcome``
            ``"exception"``. Until 0.8.0 this returned a flat frame (metric
            columns, or ``to_dataframe()`` stacked per case).
        """
        from .results import _build_case, _case_from, _frame, _schema_frame

        cases: list[Any] = []
        _prev_x: Any = None
        _prev_name: Any = None
        for idx, param_row in params.iterrows():
            overrides = {col: param_row[col] for col in params.columns}
            # Continuation (#1054): the previous row's solution seeds this
            # one. The seed carries absolute pressures, so a jump in pressure
            # level can put it far from this row's root; the solver skips its
            # retries from a given x0, so that fails fast, and the row is
            # then solved cold, exactly as solve() would, with every retry.
            result = None
            if _prev_x is not None:
                # A failed attempt is not this row's answer: its warnings
                # (did not converge) would only be noise.
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    try:
                        warm = self.solve(
                            overrides=overrides,
                            method=method,
                            init_strategy=init_strategy,
                            timeout=timeout,
                            _x0=_prev_x,
                        )
                    except Exception:
                        warm = None
                if warm is not None and warm.success:
                    result = warm
            from_prev = result is not None
            try:
                if result is None:
                    result = self.solve(
                        overrides=overrides,
                        method=method,
                        init_strategy=init_strategy,
                        timeout=timeout,
                    )
            except Exception as exc:
                warnings.warn(
                    f"NetworkRunner.sweep: solve raised an exception at index {idx}: {exc}",
                    stacklevel=2,
                )
                schema = self._schema.model_copy(deep=True)
                try:
                    self._apply_overrides(schema, overrides)
                except (KeyError, ValueError):
                    schema = self._schema
                options = {
                    "method": method,
                    "init_strategy": init_strategy,
                    "timeout": timeout,
                    "stall_handover": schema.solver_settings.stall_handover,
                    "incompressible_mach_limit": schema.solver_settings.incompressible_mach_limit,
                }
                case = _build_case(idx, _schema_frame(schema, options))
                case.values[("outputs", "_solve", "success")] = False
                case.values[("outputs", "_solve", "outcome")] = "exception"
                for m in metrics or []:
                    entity, _, q = m.partition(".")
                    case.values[("outputs", entity, q)] = math.nan
                cases.append(case)
                continue
            case = _case_from(result, idx, metrics=metrics)
            case.warm_start = _prev_name if from_prev else result._warm_start
            cases.append(case)
            if result.success:
                _prev_x = result._x_solution
                _prev_name = idx

        return _frame(cases)

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _apply_overrides(
        self,
        schema: NetworkGraphSchema,
        overrides: dict[str, Any],
    ) -> None:
        """Patch schema node data in-place with BC override values."""
        for key, value in overrides.items():
            if "." not in key:
                raise ValueError(
                    f"Override key '{key}' must be '<label>.<attribute>' (e.g. 'air_inlet.m_dot')."
                )
            label, attr = key.split(".", 1)
            node_id = self._label_to_node_id.get(label)
            if node_id is None:
                if any(n.id == label for n in schema.nodes):
                    node_id = label
                else:
                    raise KeyError(f"No node with label or id '{label}' found in this network.")
            for node in schema.nodes:
                if node.id == node_id:
                    node.data[attr] = value
                    break
