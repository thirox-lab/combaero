"""GUI hole nodes carry their Cd as ONE ``hole_cd`` object (#1098, PR 5).

``hole_cd`` is the backend's HoleCdData: the library's HoleCd field for field,
passed straight through by graph_builder. A file saved before it (up to
v0.11.1) carries ``correlation``, ``Cd`` and an impingement array's
``cd_crossflow`` loose on the node; they are lifted into ``hole_cd`` on load
and build the same network. ``data/gui_v0_11_1_hole_cd.json`` is such a file.
"""

from __future__ import annotations

import copy
import dataclasses
import json
import pathlib
import typing
import warnings

import numpy as np
import pytest
from pydantic import ValidationError

from combaero.network import HoleCd
from combaero.network.hole_cd import (
    CD_INCLINATION_MODELS,
    HOLE_CD_CORRELATIONS,
    HOLE_CD_EXIT_CROSSFLOW_MODELS,
    HOLE_CD_INLET_CROSSFLOW_MODELS,
)
from gui.backend import schemas
from gui.backend.export import FILE_VERSION
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import HoleCdData, NetworkGraphSchema, OrificeData
from gui.backend.solver_setup import make_solver

DATA = pathlib.Path(__file__).parent / "data"
LEGACY = DATA / "gui_v0_11_1_hole_cd.json"
EXAMPLES = pathlib.Path(__file__).resolve().parents[2] / "gui" / "examples"

HOLE_NODES = {
    "orifice": schemas.OrificeData,
    "effusion_plate": schemas.EffusionPlateData,
    "film_cooling": schemas.FilmCoolingData,
    "slot_film": schemas.SlotFilmData,
    "impingement_array": schemas.ImpingementArrayData,
    "effusion_liner": schemas.EffusionLinerData,
}


def _literal(name: str) -> set:
    """The values a HoleCdData field's Literal allows (None aside)."""
    ann = HoleCdData.model_fields[name].annotation
    if typing.get_origin(ann) is typing.Literal:
        return set(typing.get_args(ann))
    return {v for arg in typing.get_args(ann) for v in typing.get_args(arg)}


def test_hole_cd_data_mirrors_the_library_spec() -> None:
    assert set(HoleCdData.model_fields) == {f.name for f in dataclasses.fields(HoleCd)}
    assert _literal("correlation") == set(HOLE_CD_CORRELATIONS)
    assert _literal("inclination") == set(CD_INCLINATION_MODELS)
    assert _literal("inlet_crossflow") == set(HOLE_CD_INLET_CROSSFLOW_MODELS)
    assert _literal("exit_crossflow") == set(HOLE_CD_EXIT_CROSSFLOW_MODELS)
    assert HoleCd(**HoleCdData().model_dump()) == HoleCd()


def test_hole_cd_data_refuses_what_hole_cd_refuses() -> None:
    with pytest.raises(ValidationError, match="computes its own Cd"):
        HoleCdData(correlation="IdelchikThick", Cd=0.6)
    with pytest.raises(ValidationError, match="finite"):
        HoleCdData(correlation="fixed", Cd=float("nan"))
    with pytest.raises(ValidationError, match="Miller"):
        HoleCdData(correlation="Miller")


@pytest.mark.parametrize("ntype", sorted(HOLE_NODES))
def test_every_hole_node_carries_one_hole_cd(ntype: str) -> None:
    cls = HOLE_NODES[ntype]
    assert cls.model_fields["hole_cd"].annotation is HoleCdData
    assert cls().hole_cd == HoleCdData()
    for loose in ("correlation", "Cd", "cd_crossflow", "cd_inclination", "cd_exit_crossflow"):
        assert loose not in cls.model_fields, (ntype, loose)


def _legacy() -> dict:
    return json.loads(LEGACY.read_text())


def test_a_v0_11_1_file_lifts_its_loose_fields_into_hole_cd() -> None:
    nodes = {n.id: n for n in NetworkGraphSchema(**_legacy()).nodes}
    got = {
        nid: HOLE_NODES[n.type](**n.data).hole_cd.model_dump(exclude_defaults=True)
        for nid, n in nodes.items()
        if n.type in HOLE_NODES
    }
    assert got == {
        "inlet": {"correlation": "fixed", "Cd": 0.8},
        # A correlation's Cd was never read (#872): dropped.
        "film": {"correlation": "IdelchikThick"},
        "slot": {"correlation": "fixed", "Cd": 0.65},
        "plate": {"correlation": "Lichtarowicz"},
        "arr": {
            "correlation": "fixed",
            "Cd": 0.79,
            "exit_crossflow": "florschuetz_metzger_su_1982",
        },
    }


def test_the_lift_fills_the_v0_11_1_defaults_of_absent_keys() -> None:
    # A fixed impingement Cd without a correlation key: 'fixed' was the default.
    assert schemas.ImpingementArrayData(Cd=0.83).hole_cd.model_dump(exclude_defaults=True) == {
        "correlation": "fixed",
        "Cd": 0.83,
    }
    # 'fixed' without a Cd key: the schema's 0.6 then.
    assert OrificeData(correlation="fixed").hole_cd.Cd == 0.6
    # A metering orifice's unused Cd.
    assert OrificeData(Cd=0.6).hole_cd.model_dump(exclude_defaults=True) == {
        "correlation": "ReaderHarrisGallagher"
    }


def test_loose_fields_beside_hole_cd_are_refused() -> None:
    with pytest.raises(ValidationError, match="together"):
        OrificeData(correlation="fixed", hole_cd={"correlation": "fixed", "Cd": 0.7})


def _migrated(doc: dict) -> dict:
    doc = copy.deepcopy(doc)
    for n in doc["nodes"]:
        cls = HOLE_NODES.get(n["type"])
        if cls is not None:
            n["data"] = cls.lift_loose_hole_cd(n["data"])
    return doc


def _fingerprint(doc: dict) -> tuple[bytes, bytes]:
    schema = NetworkGraphSchema(**doc)
    solver = make_solver(build_network_from_schema(schema), schema)
    F, J = solver._residuals_and_jacobian(solver._build_x0())
    return np.asarray(F).tobytes(), J.toarray().tobytes()


def test_a_legacy_file_builds_the_same_network_as_its_migrated_twin() -> None:
    legacy = _legacy()
    new = _migrated(legacy)
    assert all("correlation" not in n["data"] for n in new["nodes"])
    assert _fingerprint(legacy) == _fingerprint(new)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = NetworkRunner.from_dict(legacy).solve(timeout=120.0)
        b = NetworkRunner.from_dict(new).solve(timeout=120.0)
    assert a.success and b.success, (a.message, b.message)
    assert a._raw["__x_solution__"] == pytest.approx(b._raw["__x_solution__"], rel=0, abs=0)


def test_the_builder_passes_hole_cd_through() -> None:
    doc = _migrated(_legacy())
    for n in doc["nodes"]:
        if n["id"] == "plate":
            n["data"]["angle_deg"] = 30.0
            n["data"]["hole_cd"] = {"correlation": "IdelchikThick", "inclination": "none"}
        if n["id"] == "arr":
            n["data"]["hole_cd"] = {"correlation": "fixed", "Cd": 0.7, "exit_crossflow": "none"}
    net = build_network_from_schema(NetworkGraphSchema(**doc))
    assert net.elements["plate"].hole_cd.inclination == "none"
    assert net.elements["plate"]._hole_cd.inclination_deg == 90.0
    row = net.elements["arr__p1"]
    assert (row.hole_cd.correlation, row.Cd, row.hole_cd.exit_crossflow) == ("fixed", 0.7, "none")


def test_the_shipped_examples_are_in_the_hole_cd_shape() -> None:
    for f in sorted(EXAMPLES.glob("*.json")):
        doc = json.loads(f.read_text())
        assert doc["version"] == FILE_VERSION, f.name
        for n in doc["nodes"]:
            if n["type"] in HOLE_NODES:
                assert not {"correlation", "Cd", "cd_crossflow"} & set(n["data"]), f.name
