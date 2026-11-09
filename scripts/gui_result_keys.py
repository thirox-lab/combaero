"""The result keys the GUI backend emits, per palette kind.

The Probe node, the Inspector's Live Telemetry, the node cards and the
Display Fields panel read solve results by key and look each key up in the
frontend's hand-kept ``QUANTITY_CATALOGUE`` (label, display unit, format).
Nothing tied that catalogue to what the backend actually emits, so keys were
renamed or added on the element side while the catalogue kept the old ones.

This script solves one small network per palette kind (the defaults the
palette drops, plus the variants that switch on optional diagnostics: a
thermal wall, a tee per junction model and direction, an orifice in a pipe,
...) through the HTTP ``/solve`` route and records every numeric key
(null included: an undefined value) the dropped node's result carries. The output,
``gui/frontend/src/utils/resultKeys.generated.json``, is

* read by the frontend: the Probe offers a kind's keys before the first solve,
  and a key with no catalogue entry still shows its SI unit;
* checked by ``python/tests/test_gui_result_keys_sync.py`` (the file is what
  this script writes now, and every key has a unit in ``gui/backend/units.py``);
* checked by ``gui/frontend/src/utils/quantities.test.ts`` (every catalogue
  entry is emitted by some kind, in a unit consistent with the backend's).

Usage (repo root)::

    uv run python scripts/gui_result_keys.py          # rewrite the JSON
    uv run python scripts/gui_result_keys.py --check  # exit 1 when stale
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from collections.abc import Iterator
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "gui" / "frontend" / "src" / "utils" / "resultKeys.generated.json"
# The palette networks and the HTTP client live with the route tests.
for _p in (REPO, REPO / "python" / "tests"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

#: Kinds that are not solved elements or nodes (display-only).
_DISPLAY_ONLY = {"probe"}
#: The pseudo-kind under which thermal-wall (edge) result keys are listed.
THERMAL_WALL = "thermal_wall"


def _node(nid: str, ntype: str, data: dict | None = None) -> dict:
    return {"id": nid, "type": ntype, "position": {"x": 0.0, "y": 0.0}, "data": data or {}}


def _edge(
    eid: str,
    src: str,
    tgt: str,
    sh: str | None = None,
    th: str | None = None,
    data: dict | None = None,
) -> dict:
    e: dict = {"id": eid, "source": src, "target": tgt, "data": data or {}}
    if sh:
        e["sourceHandle"] = sh
    if th:
        e["targetHandle"] = th
    return e


def _with(graph: dict, node_id: str, **data: object) -> dict:
    for n in graph["nodes"]:
        if n["id"] == node_id:
            n["data"] = {**n["data"], **data}
    return graph


def _merge_tee(data: dict) -> dict:
    """Two inlets on straight and branch, one outlet on common."""
    pb = "pressure_boundary"
    hi, lo = {"Pt": 1.05e5, "Tt": 300.0}, {"Pt": 1.0e5, "Tt": 300.0}
    ch = {"L": 0.5, "D": 0.05}
    return {
        "nodes": [
            _node("i1", pb, hi),
            _node("i2", pb, hi),
            _node("out", pb, lo),
            _node("c1", "channel", ch),
            _node("c2", "channel", ch),
            _node("c3", "channel", ch),
            _node("x", "mpce_tee", {"flow_direction": "merge", **data}),
        ],
        "edges": [
            _edge("e1", "i1", "c1"),
            _edge("e2", "c1", "x", th="port-straight-target"),
            _edge("e3", "i2", "c2"),
            _edge("e4", "c2", "x", th="port-branch-target"),
            _edge("e5", "x", "c3", sh="port-common-source"),
            _edge("e6", "c3", "out"),
        ],
        "solver_settings": {"timeout": 20.0},
    }


def _heated_pair() -> dict:
    """Two channels sharing a thermal wall: edge results, heated-duct flags."""
    pb = "pressure_boundary"
    wall = {
        "type": "thermal",
        "thickness": 0.003,
        "conductivity": 20.0,
        "area": 0.05,
        "layers": [{"thickness": 0.003, "conductivity": 20.0, "material": "generic"}],
    }
    return {
        "nodes": [
            _node("hin", pb, {"Pt": 1.05e5, "Tt": 900.0}),
            _node("hout", pb, {"Pt": 1.0e5, "Tt": 900.0}),
            _node("cin", pb, {"Pt": 1.05e5, "Tt": 300.0}),
            _node("cout", pb, {"Pt": 1.0e5, "Tt": 300.0}),
            _node("x", "channel", {"L": 0.5, "D": 0.05}),
            _node("c", "channel", {"L": 0.5, "D": 0.05}),
        ],
        "edges": [
            _edge("e1", "hin", "x"),
            _edge("e2", "x", "hout"),
            _edge("e3", "cin", "c"),
            _edge("e4", "c", "cout"),
            _edge("w", "x", "c", data=wall),
        ],
        "solver_settings": {"timeout": 20.0},
    }


def _radiating_pair(side_a: str, side_b: str) -> dict:
    """The heated pair with radiation on both faces (#793, #1071): each
    face's flux, rate and surface temperature, and a gas face's emissivity."""
    graph = _heated_pair()
    graph["edges"][-1]["data"].update(radiation_a={"kind": side_a}, radiation_b={"kind": side_b})
    return graph


def _walled_loss(graph: dict) -> dict:
    """A discrete loss walled to a warm channel: its convective diagnostics.

    The GUI's loss has no wetted-area field; the wall gives it one (#1022),
    so only a walled loss reports h, Nu and T_aw.
    """
    pb = "pressure_boundary"
    graph["nodes"][0] = _node("in", "mass_boundary", {"m_dot": 0.05, "Tt": 300.0})
    graph["nodes"] += [
        _node("win", pb, {"Pt": 1.05e5, "Tt": 500.0}),
        _node("wout", pb, {"Pt": 1.0e5, "Tt": 500.0}),
        _node("cw", "channel", {"L": 0.5, "D": 0.05}),
    ]
    wall = {"type": "thermal", "thickness": 0.003, "conductivity": 20.0, "area": 0.01}
    graph["edges"] += [
        _edge("w1", "win", "cw"),
        _edge("w2", "cw", "wout"),
        _edge("w", "x", "cw", data=wall),
    ]
    return graph


def _orifice_in_pipe(data: dict) -> dict:
    """An orifice with an approach pipe: beta > 0, the ISO 5167 taps."""
    pb = "pressure_boundary"
    return {
        "nodes": [
            _node("in", pb, {"Pt": 1.05e5, "Tt": 300.0}),
            _node("out", pb, {"Pt": 1.0e5, "Tt": 300.0}),
            _node("c1", "channel", {"L": 0.5, "D": 0.1}),
            _node("c2", "channel", {"L": 0.5, "D": 0.1}),
            _node("x", "orifice", {"diameter": 0.05, **data}),
        ],
        "edges": [
            _edge("e1", "in", "c1"),
            _edge("e2", "c1", "x"),
            _edge("e3", "x", "c2"),
            _edge("e4", "c2", "out"),
        ],
        "solver_settings": {"timeout": 20.0},
    }


def networks() -> Iterator[tuple[str, str, dict]]:
    """(kind, variant, graph): the node to read is always ``x``."""
    import test_gui_effusion_plate as eff
    import test_gui_http_routes as http
    import test_gui_impingement_array as imp

    for kind in http._palette_kinds():
        if kind in _DISPLAY_ONLY:
            continue
        graph = http._palette_network(kind)
        if kind == "discrete_loss":
            # A constant-fraction loss between two pressure boundaries fixes
            # Pt_out / Pt_in and leaves the flow undetermined: feed a flow.
            graph["nodes"][0] = _node("in", "mass_boundary", {"m_dot": 0.05, "Tt": 300.0})
        if kind == "vortex":
            # A free vortex conserves Pt, so two pressure boundaries leave no
            # solution: feed a flow, at a shaft speed (the palette's null is
            # the global setting, unset here, i.e. no swirl).
            _with(graph, "x", shaft_speed_rpm=3000.0)
            graph["nodes"][0] = _node("in", "mass_boundary", {"m_dot": 0.05, "Tt": 300.0})
        yield kind, "palette", graph

    yield "channel", "heated (thermal wall)", _heated_pair()
    yield "channel", "incompressible", _with(_heated_pair(), "x", regime="incompressible")
    yield "channel", "heated, gas radiation on A, casing on B", _radiating_pair("gas", "casing")
    yield "channel", "heated, casing radiation on A, gas on B", _radiating_pair("casing", "gas")
    yield "orifice", "in a pipe", _orifice_in_pipe({})
    for model in ("nozzle", "jobson1955", "iso5167"):
        yield "orifice", f"in a pipe, {model}", _orifice_in_pipe({"compressible_model": model})
    for corr in ("Lichtarowicz", "IdelchikThick", "McGreehanSchotsch", "fixed"):
        # A 0.1 m plate: l/d = 2, inside every discharge-hole correlation.
        yield (
            "orifice",
            f"in a pipe, {corr}",
            _orifice_in_pipe({"correlation": corr, "plate_thickness": 0.1}),
        )
    for model in ("mynard", "constant_k", "ideal_momentum"):
        yield (
            "mpce_tee",
            f"branch, {model}",
            _with(http._palette_network("mpce_tee"), "x", junction_model=model),
        )
        yield "mpce_tee", f"merge, {model}", _merge_tee({"junction_model": model})
    yield (
        "ejector",
        "kracik_dvorak",
        _with(http._palette_network("ejector"), "x", anchor_closure="kracik_dvorak"),
    )
    jet_pump = http._palette_network("ejector")
    jet_pump["solver_settings"]["global_regime"] = "incompressible"
    yield "ejector", "incompressible (jet pump)", jet_pump
    contraction = http._palette_network("area_change")
    for n in contraction["nodes"]:
        if n["id"] == "c1":
            n["data"] = {"L": 0.5, "D": 0.08}
        if n["id"] == "c2":
            n["data"] = {"L": 0.5, "D": 0.05}
    yield "area_change", "contraction", contraction
    yield (
        "area_change",
        "conical",
        _with(http._palette_network("area_change"), "x", model_type="conical", length=0.05),
    )
    for corr in ("constant_head", "linear_theta_fraction", "linear_theta_head"):
        graph = http._palette_network("discrete_loss")
        graph["nodes"][0] = _node("in", "mass_boundary", {"m_dot": 0.05, "Tt": 300.0})
        yield "discrete_loss", corr, _with(graph, "x", correlation_type=corr)
    yield "discrete_loss", "thermal wall", _walled_loss(http._palette_network("discrete_loss"))
    yield (
        "heat_source",
        "outlet Tt (total)",
        _with(http._palette_network("heat_source"), "x", mode="Tt_target", Tt_target=320.0),
    )
    yield (
        "heat_source",
        "fraction",
        _with(http._palette_network("heat_source"), "x", mode="fraction", fraction=0.01),
    )
    yield (
        "combustor",
        "equilibrium",
        _with(http._palette_network("combustor"), "x", method="equilibrium"),
    )
    yield (
        "effusion_liner",
        "baldauf_sellers film",
        _with(http._palette_network("effusion_liner"), "x", gas_film="baldauf_sellers"),
    )

    def _renamed(graph: dict, old: str) -> dict:
        for n in graph["nodes"]:
            if n["id"] == old:
                n["id"] = "x"
        for e in graph["edges"]:
            e["source"] = "x" if e["source"] == old else e["source"]
            e["target"] = "x" if e["target"] == old else e["target"]
        graph["solver_settings"] = {"timeout": 20.0}
        return graph

    yield "effusion_plate", "on a gas-side chamber", _renamed(eff._liner(), "eff")
    yield (
        "effusion_plate",
        "on a gas-side chamber, baldauf_sellers film",
        _renamed(eff._liner({"gas_film": "baldauf_sellers"}), "eff"),
    )
    yield (
        "effusion_plate",
        "on a gas-side chamber, flat-plate gas baseline",
        _renamed(eff._liner({"gas_baseline": "flat_plate", "gas_start_length": 0.02}), "eff"),
    )
    yield (
        "effusion_plate",
        "on a gas-side chamber, gas radiation",
        _renamed(eff._liner({"radiation": {"kind": "gas"}}), "eff"),
    )
    yield (
        "impingement_array",
        "bypass from a boundary",
        _renamed(imp._bypass_schema("boundary"), "arr"),
    )
    yield "impingement_array", "thermal wall", _renamed(imp._schema(with_wall=True), "arr")


def _numeric(v: object) -> bool:
    """A number, or null: /solve writes an undefined (NaN) or unbounded
    (inf) value as null, and the KEY is still the result's -- the list is the
    schema, not the values one solve happened to produce."""
    return v is None or isinstance(v, (int, float)) and not isinstance(v, bool)


def _numeric_keys(res: dict) -> set[str]:
    keys = set()
    for k, v in res.items():
        if k == "state" and isinstance(v, dict):
            keys |= {f"state.{sk}" for sk, sv in v.items() if _numeric(sv)}
        elif _numeric(v):
            keys.add(k)
    return keys


def collect(verbose: bool = False) -> dict:
    """Solve every network and return the JSON document."""
    import test_gui_http_routes as http

    from gui.backend.units import unit_of

    kinds: dict[str, set[str]] = {}
    variants: dict[str, list[str]] = {}
    for kind, variant, graph in networks():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = http._post("/solve", graph)
        if r.status_code != 200:
            raise RuntimeError(f"{kind} ({variant}): HTTP {r.status_code}: {r.text[:300]}")
        body = r.json()
        if not body["success"]:
            # A failed solve may stop before every diagnostic is computed; the
            # key list must not depend on where an iteration gave up.
            raise RuntimeError(f"{kind} ({variant}) did not converge: {body['message'][:300]}")
        res = body["element_results"].get("x") or body["node_results"].get("x")
        if res is None:
            raise RuntimeError(f"{kind} ({variant}): no result for the dropped node")
        kinds.setdefault(kind, set()).update(_numeric_keys(res))
        for edge_res in body.get("edge_results", {}).values():
            kinds.setdefault(THERMAL_WALL, set()).update(_numeric_keys(edge_res))
        variants.setdefault(kind, []).append(variant)
        if verbose:
            print(f"{kind:>20} | {variant:<45} success={body['success']}", file=sys.stderr)

    every = sorted(set().union(*kinds.values()))
    units = {}
    for key in every:
        raw = key.split(".", 1)[-1]
        # None: no unit in gui/backend/units.py (the sync test names it).
        units[raw] = unit_of(raw)
    return {
        "_generated_by": "uv run python scripts/gui_result_keys.py (do not edit by hand)",
        "kinds": {k: sorted(v) for k, v in sorted(kinds.items())},
        "variants": dict(sorted(variants.items())),
        "units": dict(sorted(units.items())),
    }


def render(doc: dict) -> str:
    return json.dumps(doc, indent="\t", ensure_ascii=True) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true", help="exit 1 when the file is stale")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()
    text = render(collect(verbose=args.verbose))
    if args.check:
        if OUT.read_text(encoding="utf-8") != text:
            print(f"{OUT.relative_to(REPO)} is stale: rerun scripts/gui_result_keys.py")
            return 1
        return 0
    OUT.write_text(text, encoding="utf-8")
    print(f"wrote {OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
