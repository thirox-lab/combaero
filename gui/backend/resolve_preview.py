"""Resolve a network's geometry inheritance WITHOUT solving it (#559).

The solver's `resolve_all_topology` is the authoritative geometry-inheritance
pass: each element/node fills its diameter / area / F0 / F1 / port areas from its
neighbours, or refuses (raises) when the value IS the physics and nothing can be
inferred. This module runs that pass ahead of a solve and reports, per element
and per geometry-bearing node:

  * the resolved value of each geometry field,
  * its provenance (`"user"` / `"inherited from <id>"` / `"channel '<id>'"` /
    `"default"`), read from the `_area_source` / `_F0_source` / ... attributes
    the elements record (B1a/B1b), and
  * whether it is unresolved -- the element would refuse, so the network is not
    yet solvable, and why.

The GUI reads this instead of re-deriving inheritance in its own graph walk, so
the placeholder it shows and the value the solver uses are the same number, and a
"you must set this" is surfaced before Solve rather than as a bare failure.

Read-only to the caller: it resolves the given network in place (as a solve
would) but does not solve or mutate solver state.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import combaero as cb
from combaero.network.components import (
    AreaChangeElement,
    ChannelElement,
    CombustorNode,
    MomentumChamberNode,
    MultiPortChamberBase,
    OrificeElement,
    PressureLossElement,
    TeeJunctionElement,
)
from combaero.network.ejector_element import EjectorElement

if TYPE_CHECKING:
    from combaero.network.graph import FlowNetwork


@dataclass
class FieldResolution:
    """One geometry field's resolved value and where it came from."""

    value: float | None
    #: "user" | "inherited from <id>" | "channel '<id>'" | "default" | "" (unset)
    source: str


@dataclass
class ItemResolution:
    """Resolution of one element or geometry-bearing node."""

    #: "resolved" (every field has a value) or "unresolved" (it would refuse).
    status: str
    fields: dict[str, FieldResolution] = field(default_factory=dict)
    #: When unresolved, the refusal message (which field to set, and where it
    #: searched). Empty when resolved.
    reason: str = ""
    #: A warning the solve will issue about this item's geometry, which does
    #: not stop it (#1023: a momentum chamber narrower than an adjoining
    #: pipe). Empty when there is none.
    warning: str = ""
    #: A resolved geometry the solver will refuse at validate() (#1024):
    #: overlapping holes, a zero-area duct or node, a Dh no section of its
    #: area has. Empty when there is nothing to refuse.
    refusal: str = ""


def _s(obj: Any, name: str) -> str:
    """A source attribute, defaulting to "" when absent."""
    return str(getattr(obj, name, "") or "")


def _extract(item: Any) -> dict[str, FieldResolution]:
    """The geometry fields and provenance of one element or node, by type.

    Only covers the items that inherit/derive geometry; anything else returns an
    empty dict and is reported as trivially resolved.
    """
    # Builder-resolved geometry (a GUI tee's F_C / F_branch are resolved in
    # graph_builder and handed to the element as explicit port areas, so the
    # element itself would report "user"). The builder records the real source.
    preview_fields = getattr(item, "_preview_fields", None)
    if preview_fields:
        return {
            name: FieldResolution(float(value) if value is not None else None, str(source))
            for name, (value, source) in preview_fields.items()
        }
    # Nodes.
    if isinstance(item, (MomentumChamberNode, CombustorNode)):
        return {
            "area": FieldResolution(_f(item, "area"), _s(item, "_area_source")),
            "Dh": FieldResolution(_f(item, "Dh"), _s(item, "_Dh_source")),
        }
    # Elements. Order matters: Effusion/Impingement plates subclass Orifice but
    # take their geometry from construction, so they are not inheriting items.
    if isinstance(item, ChannelElement):
        return {
            "diameter": FieldResolution(_f(item, "diameter"), _s(item, "_area_source")),
            "Dh": FieldResolution(_f(item, "Dh"), _s(item, "_Dh_source")),
        }
    if isinstance(item, AreaChangeElement):
        # Equal faces are no area change (#1023): pass_through = 1, zero
        # loss. An explicit face contradicting its pipe is refused, so it
        # arrives here as unresolved.
        pass_through = bool(item._pass_through)
        return {
            "F0": FieldResolution(_f(item, "F0"), _s(item, "_F0_source")),
            "F1": FieldResolution(_f(item, "F1"), _s(item, "_F1_source")),
            "pass_through": FieldResolution(
                1.0 if pass_through else 0.0,
                "F0 = F1: no area change" if pass_through else "",
            ),
        }
    if isinstance(item, TeeJunctionElement):
        return {"F_C": FieldResolution(_f(item, "F_C"), _s(item, "_area_source"))}
    if isinstance(item, EjectorElement):
        # Throat and nozzle exit are internal geometry (always the user's or the
        # schema default); only the mixing duct A_3 inherits, from the outlet pipe.
        return {
            "throat_area": FieldResolution(_f(item, "throat_area"), "user"),
            "nozzle_exit_area": FieldResolution(_f(item, "nozzle_exit_area"), "user"),
            "mixing_area": FieldResolution(_f(item, "mixing_area"), _s(item, "_mixing_source")),
        }
    if isinstance(item, MultiPortChamberBase):
        sources = getattr(item, "_port_area_sources", None) or [""] * item.N
        areas = getattr(item, "port_areas", None) or [None] * item.N
        return {
            f"port_{i}_area": FieldResolution(
                float(areas[i]) if areas[i] is not None else None, str(sources[i] or "")
            )
            for i in range(item.N)
        }
    if isinstance(item, PressureLossElement):
        return {"area": FieldResolution(_f(item, "area"), _s(item, "_area_source"))}
    if type(item) is OrificeElement:  # not the plate subclasses
        # The bore never inherits (it is the restriction); the upstream pipe's
        # diameter does, and sets beta = d / D for the velocity-of-approach
        # factor. beta = 0 means no correction: no single upstream channel, or
        # in_pipe=False (a lumped Cd*A restriction). A bore equal to the pipe
        # is no orifice (pass_through = 1, beta = 1, zero loss); a larger one
        # is refused, so it arrives here as unresolved (#788). The same for a
        # bore against an adjacent momentum chamber's area (#1016); the
        # field's source says which. Reversed flow
        # approaches through the to-node's single pipe (#554). E comes from the
        # C++ factor the kernels use, never re-derived in the frontend.
        beta = _f(item, "beta") or 0.0
        beta_rev = _f(item, "beta_reverse") or 0.0
        return {
            **_pass_through_field(item),
            "diameter": FieldResolution(_f(item, "diameter"), _s(item, "_area_source")),
            "pipe_diameter": FieldResolution(
                _f(item, "upstream_diameter"), _s(item, "_upstream_source")
            ),
            "beta": FieldResolution(beta, "derived from bore / pipe D"),
            "E": FieldResolution(float(cb.velocity_of_approach_factor(beta)), "derived from beta"),
            "pipe_diameter_reverse": FieldResolution(
                _f(item, "_reverse_pipe_diameter"), _s(item, "_reverse_pipe_source")
            ),
            "beta_reverse": FieldResolution(beta_rev, "derived from bore / pipe D"),
            "E_reverse": FieldResolution(
                float(cb.velocity_of_approach_factor(beta_rev)), "derived from beta"
            ),
            # The compressible model the element runs (#861): 1 = jobson1955,
            # 0 = anything else; selected at build unless the user set one.
            "jobson1955": FieldResolution(
                1.0 if getattr(item, "compressible_model", "nozzle") == "jobson1955" else 0.0,
                "user"
                if getattr(item, "_compressible_model_explicit", True)
                else "selected at build",
            ),
        }
    if isinstance(item, OrificeElement):
        # The plates (effusion, film, slot, impingement) take their geometry
        # from construction; only the verdict against an adjacent momentum
        # chamber is resolved (#1016). A wider one is refused: unresolved.
        return _pass_through_field(item)
    return {}


def _pass_through_field(item: Any) -> dict[str, FieldResolution]:
    """1 when the orifice is no orifice (#788 bore = pipe, #1016 bore area =
    an adjacent momentum chamber's), with the element's own reason."""
    pass_through = bool(getattr(item, "_pass_through", False))
    return {
        "pass_through": FieldResolution(
            1.0 if pass_through else 0.0,
            _s(item, "_pass_through_reason") if pass_through else "",
        )
    }


def _f(obj: Any, name: str) -> float | None:
    v = getattr(obj, name, None)
    return float(v) if v is not None else None


def resolve_preview(net: FlowNetwork) -> dict[str, ItemResolution]:
    """Resolve ``net``'s geometry inheritance and report it per item.

    Runs the solver's own resolve (``FlowNetwork._resolve_to_fixed_point``,
    the passes ``resolve_all_topology`` repeats until the geometry settles)
    but catches each item's refusal (``ValueError``) instead of stopping at
    the first, so a partially-specified network still yields a full report.
    The geometry it reports is the one a solve uses: a node that inherits
    from an element whose own geometry is inherited settles only on a later
    pass (#1148). Keys are element and node ids (including the builder's
    synthetic ids; the caller maps them to canvas nodes).
    """
    refusals: dict[str, str] = {}
    net._resolve_to_fixed_point(refusals)

    report: dict[str, ItemResolution] = {}
    for item_id, item in {**net.nodes, **net.elements}.items():
        reason = refusals.get(item_id, "")
        fields = _extract(item)
        # The warning FlowNetwork.validate issues (#1023), from the same
        # method, on the resolved geometry.
        warning = (
            " ".join(item.pipe_area_warnings(net)) if isinstance(item, MomentumChamberNode) else ""
        )
        if isinstance(item, OrificeElement) and item.closed:
            # A zero bore (#1024): closed, a wall in the flow path, no flow.
            fields["closed"] = FieldResolution(1.0, "zero area: closed, no flow")
        if reason:
            report[item_id] = ItemResolution(
                status="unresolved", fields=fields, reason=reason, warning=warning
            )
        else:
            report[item_id] = ItemResolution(
                status="resolved",
                fields=fields,
                warning=warning,
                refusal=_geometry_refusal(item),
            )
    return report


def _geometry_refusal(item: Any) -> str:
    """The element's or node's own geometry refusal (#1024), the message its
    validate() raises, or "": checked on the resolved values, before Solve."""
    check = getattr(item, "_geometry_error", None)
    if not callable(check):
        return ""
    try:
        return str(check() or "")
    except (ValueError, TypeError, ArithmeticError):
        return ""
