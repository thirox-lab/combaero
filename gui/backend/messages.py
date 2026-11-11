"""User-facing wording for solver errors and messages (#321 item 1).

The network elements raise with their internal ids: a tee or ejector port is
the auto-inserted node ``__tee_jct__<node id>_<port>``, a sharp step the
builder put between a tee port and a channel of another bore is
``__tee_step__<node id>_<port>`` (its far node adds ``__node``), and every
element is its React Flow id. The GUI shows neither -- a tee's ports are labelled c, s, b
and an ejector's p, s, o -- so the message is rewritten here, at the HTTP
boundary, and the elements keep raising what a script user needs.
"""

from __future__ import annotations

import re

from .schemas import SPLITTER_MAX_OUTLETS, NetworkGraphSchema

#: handle port name -> the letter the node draws next to that port
PORT_TAG = {
    "common": "c",
    "straight": "s",
    "branch": "b",
    "primary": "p",
    "secondary": "s",
    "outlet": "o",
    "coolantin": "coolant in",
    "coolantout": "coolant out",
    "discharge": "discharge",
    # A flow splitter's (#903) inlet and outlets, drawn unlabelled.
    "in": "in",
    **{f"out{i}": f"outlet {i + 1}" for i in range(SPLITTER_MAX_OUTLETS)},
}

#: node type -> the name the palette uses, where the type string differs
KIND_NAME = {
    "mpce_tee": "Tee",
    "momentum_chamber": "Momentum Chamber",
    "mass_boundary": "Mass Flow Boundary",
    "fuel_boundary": "Fuel Boundary",
    "compressor": "Compressor",
    "flow_splitter": "Flow Splitter",
    "turbine": "Turbine",
    "pressure_boundary": "Pressure Boundary",
}

#: solver class name -> what the Tee Inspector's "Junction Model" select calls it
MODEL_NAME = {
    "IdealMomentumTeeElement": 'the "Ideal momentum (parameter-free)" junction model',
    "MultiPortChamberElement": 'Tee (junction model "Momentum CV (Mynard)")',
    "ConstantKTeeElement": 'Tee (junction model "Constant K (handbook)")',
}

#: solver attribute wording -> the Inspector field it is set by
FIELD_NAME = (
    ("(flow_direction='merge')", '(Flow Direction "Merge")'),
    ("(flow_direction='branch')", '(Flow Direction "Branch")'),
    ("change flow_direction", "change the tee's Flow Direction"),
)

_PORT_ID = re.compile(r"__tee_(jct|step)__(.+?)_(" + "|".join(PORT_TAG) + r")(?:__node)?\b")

HINTS = (
    (
        "is declared as an",
        "Hint: this usually means the junction's flow direction (merge / branch) "
        "does not match how it is wired.",
    ),
    (
        "must be a MomentumChamberNode",
        "Hint: connect this port through an element (for example a channel), "
        "not straight to a boundary or plenum.",
    ),
    (
        # The Inspector shows the hole-angle field for IdelchikThick only.
        "inclination='idelchik_3_2') is Idelchik",
        'Hint: the Hole Cd section\'s Hole Angle Term "idelchik_3_2" needs the '
        'Discharge Model "Deep Hole in Wall (Idelchik 4-18a)"; select it to see the '
        "field and set the term to Auto or None, or keep that Discharge Model.",
    ),
)


def node_name(schema: NetworkGraphSchema, node_id: str) -> str:
    """The user's label for a node, with its kind; the kind alone if unlabelled."""
    node = next((n for n in schema.nodes if n.id == node_id), None)
    if node is None:
        return node_id
    kind = KIND_NAME.get(node.type, node.type.replace("_", " ").title())
    label = (node.data or {}).get("label")
    return f"{label} ({kind})" if label else kind


def humanize(message: str, schema: NetworkGraphSchema) -> str:
    """Rewrite internal ids in a solver message into what the canvas shows."""
    if not message:
        return message
    splitters = {n.id for n in schema.nodes if n.type == "flow_splitter"}
    out = _PORT_ID.sub(
        lambda m: (
            # A flow splitter's ports are drawn unlabelled: "<node> outlet 2".
            f"{node_name(schema, m[2])} {PORT_TAG[m[3]]}"
            if m[2] in splitters
            else ("the area step at " if m[1] == "step" else "")
            + f"{node_name(schema, m[2])} port {m[3]} ({PORT_TAG[m[3]]})"
        ),
        message,
    )
    for node in schema.nodes:
        quoted = f"'{node.id}'"
        if quoted in out:
            out = out.replace(quoted, f"'{node.id}' [{node_name(schema, node.id)}]")
    # The declared-direction hint (#640) names solver classes and attributes;
    # a GUI user picks those in the Tee Inspector under other names.
    for cls, name in MODEL_NAME.items():
        out = out.replace(cls, name)
    for raw, shown in FIELD_NAME:
        out = out.replace(raw, shown)
    hints = [hint for trigger, hint in HINTS if trigger in out]
    return "\n".join([out, *hints])
