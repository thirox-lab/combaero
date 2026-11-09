import math as _math

from pydantic import ValidationError

from combaero.network import (
    AreaChangeElement,
    ChannelElement,
    CombustorNode,
    CompressorElement,
    ConstantFractionLoss,
    ConstantHeadLoss,
    ConvectiveSurface,
    EffusionLiner,
    EffusionPlateElement,
    EnergyBoundary,
    FilmCoolingElement,
    FlowNetwork,
    FlowSplitter,
    FuelBoundary,
    ImpingementArray,
    ImpingementModel,
    LinearThetaFractionLoss,
    LinearThetaHeadLoss,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    MultiPortChamberBase,
    OrificeElement,
    PinFinModel,
    PlenumNode,
    PortDischargeElement,
    PressureBoundary,
    PressureLossElement,
    RibbedModel,
    SingleJetImpingementModel,
    SlotFilmElement,
    SmoothModel,
    ThermalWall,
    TurbineElement,
    VortexElement,
    WallLayer,
    WallNode,
)
from combaero.network.components import DEFAULT_CHANNEL_DIAMETER
from combaero.network.ejector_element import EjectorElement
from combaero.network.ideal_tee_element import IdealMomentumTeeElement
from combaero.network.mpce_element import ConstantKTeeElement, MultiPortChamberElement

from .schemas import (
    AreaChangeData,
    ChannelData,
    CombustorData,
    CompositionData,
    ConstantFractionLossData,
    ConstantHeadLossData,
    DiscreteLossData,
    EffusionLinerData,
    EffusionPlateData,
    EjectorData,
    FilmCoolingData,
    FlowSplitterData,
    FuelBoundaryData,
    HeatSourceData,
    ImpingementArrayData,
    ImpingementModelData,
    LinearThetaFractionLossData,
    LinearThetaHeadLossData,
    MassBoundaryData,
    MomentumChamberData,
    MPCETeeData,
    NetworkGraphSchema,
    OrificeData,
    PinFinModelData,
    PlenumData,
    PressureBoundaryData,
    RibbedModelData,
    SingleJetImpingementModelData,
    SlotFilmData,
    SmoothModelData,
    ThermalWallData,
    TurbomachineData,
    VortexData,
    WallData,
    removed_surface_message,
)

#: |channel area / port area - 1| above which a tee port and its directly
#: connected channel are joined through a sharp step (#536). Below it the bores
#: count as the same size, so round-off in a typed diameter inserts no step.
#: The same bound as the junction's unmodelled-step diagnostic (#539), so a
#: step the builder leaves out is never one that diagnostic reports.
TEE_PORT_STEP_TOL = MultiPortChamberBase.PORT_AREA_MISMATCH_TOL

#: The impingement array's bypass handle (#467): an initial crossflow into the
#: gap's upstream end. Must match ImpingementArrayNode.tsx.
ARRAY_BYPASS_HANDLE = "port-bypass-target"


#: Canvas kinds the second pass builds into elements. Any other kind that is
#: not a physical node (first pass) or a probe is refused: skipping it would
#: solve the network without it.
ELEMENT_KINDS = frozenset(
    {
        "channel",
        "orifice",
        "effusion_plate",
        "film_cooling",
        "slot_film",
        "impingement_array",
        "area_change",
        "lossless_connection",
        "discrete_loss",
        "vortex",
        "mpce_tee",
        "ejector",
        "effusion_liner",
        "compressor",
        "turbine",
        "flow_splitter",
    }
)


#: A flow splitter's (#903) per-outlet results, reported on its canvas node as
#: ``rows_<key>`` lists in port order (out0, out1, ...).
SPLITTER_ROW_KEYS = (
    "m_dot",
    "fraction",
    "fraction_achieved",
    "valve_dPt",
    "valve_CdA",
    "valve_Cd",
    "valve_zeta",
    "valve_needs_pump",
)


class GuiFlowSplitter:
    """A ``FlowSplitter`` as one canvas node (#903): ``summarize`` gives the
    runner its result, scalars (the reference in use and the one started
    from) plus a ``rows_<key>`` list per outlet in port order, None where an outlet has no such value (the open reference has no
    metering; Cd and zeta need a reference area)."""

    def __init__(self, splitter: FlowSplitter) -> None:
        self.splitter = splitter

    def summarize(self, element_diag: dict[str, dict]) -> dict:
        s = self.splitter.summarize(element_diag)
        outs = list(self.splitter.outlets)
        rows = [s["outlets"][o] for o in outs]
        out: dict = {
            "m_dot": float(s["m_dot"]),
            "n_outlets": float(len(outs)),
            # The reference the solve ended with, and the one it started from
            # (they differ when it was moved, #903). Counted from 1, as the
            # card and Inspector number the outlets (port out<i> is i + 1).
            "reference_outlet": float(outs.index(s["reference"]) + 1),
            "reference_start": float(outs.index(s["reference_start"]) + 1),
            "reference_swapped": float(s["reference"] != s["reference_start"]),
            "n_needs_pump": float(len(s["needs_pump"])),
        }
        for key in SPLITTER_ROW_KEYS:
            out[f"rows_{key}"] = [float(r[key]) if r.get(key) is not None else None for r in rows]
        # Always reported, NaN when no outlet row carries one (no flow to
        # split): a result's keys follow the network, not the operating point.
        dpt = [float(r["valve_dPt"]) for r in rows if "valve_dPt" in r]
        out["valve_dPt_max"] = max(dpt) if dpt else _math.nan
        err = [
            abs(float(r["fraction_achieved"]) - float(r["fraction"]))
            for r in rows
            if "fraction_achieved" in r
        ]
        out["fraction_error_max"] = max(err) if err else _math.nan
        return out


#: The heat source canvas kind (#792) and the kinds it may heat.
HEAT_SOURCE_KIND = "heat_source"
HEAT_SOURCE_TARGETS = frozenset({"plenum", "momentum_chamber", "combustor"})


def _parse(model, data: dict | None, owner_id: str, what: str = "Node"):
    """Validate one node's (or wall's) data, naming it in the error.

    Pydantic's own message names the model ("ChannelData"), not the node,
    and carries a documentation link per error; the user needs which node,
    which field and why.
    """
    try:
        return model(**(data or {}))
    except ValidationError as e:
        parts = []
        for err in e.errors():
            loc = ".".join(str(p) for p in err["loc"]) or "data"
            parts.append(f"{loc}: {err['msg'].removeprefix('Value error, ')}")
        raise ValueError(f"{what} '{owner_id}': " + "; ".join(parts)) from e


def _fixed_cd(
    data: OrificeData | EffusionPlateData | FilmCoolingData | SlotFilmData | ImpingementArrayData,
) -> float | None:
    """The Cd an orifice-family element is built with (#872): the field's
    value for the 'fixed' correlation, None for any other, which computes its
    own. The node keeps a Cd field whatever its correlation (the Inspector
    shows it for 'fixed' only), and the element refuses a Cd given with a
    correlation; a file saved with both loads, and its Cd is not used, as it
    never was."""
    return data.Cd if data.correlation == "fixed" else None


def _cd_inclination(data: EffusionPlateData | FilmCoolingData) -> str:
    """The opt-in hole-angle Cd term an inclined cooling-hole element is built
    with: the field's value with 'IdelchikThick', the only correlation it
    composes with, and "none" otherwise. Like _fixed_cd, the node keeps the
    field whatever its correlation (the Inspector shows it for IdelchikThick
    only), so switching the correlation does not leave an opt-in the element
    would refuse."""
    return data.cd_inclination if data.correlation == "IdelchikThick" else "none"


def _checked_edges(schema: NetworkGraphSchema) -> NetworkGraphSchema:
    """Refuse edges naming a node not on the canvas; drop edges to probes.

    A probe is display-only, but its ``probe-in`` handle takes a flow edge;
    counted, that edge was the element's second outlet. An edge to a missing
    node was ignored, or blamed on its element as "not fully connected".
    """
    ids = {n.id for n in schema.nodes}
    for edge in schema.edges:
        for end in (edge.source, edge.target):
            if end not in ids:
                raise ValueError(
                    f"Edge '{edge.id}' connects '{edge.source}' to '{edge.target}', "
                    f"but '{end}' is not on the canvas."
                )
    # Probes are display-only; a heat source's one edge names the node it
    # heats (_heat_source_links), not a flow path.
    dropped = {n.id for n in schema.nodes if n.type in ("probe", HEAT_SOURCE_KIND)}
    if not dropped:
        return schema
    edges = [e for e in schema.edges if e.source not in dropped and e.target not in dropped]
    return schema.model_copy(update={"edges": edges})


def _heat_source_links(schema: NetworkGraphSchema) -> dict[str, str]:
    """{heat source id: the canvas node it heats} (#792).

    A heat source is wired by exactly ONE edge, to a plenum, momentum chamber
    or combustor. Anything else is refused by name: an unwired source would
    solve the network without its heat, an element has no energy balance of
    its own to take it, and two nodes would leave its heat split by nothing.
    """
    kinds = {n.id: n.type for n in schema.nodes}
    links: dict[str, str] = {}
    for n in schema.nodes:
        if n.type != HEAT_SOURCE_KIND:
            continue
        ends = sorted(
            {
                e.target if e.source == n.id else e.source
                for e in schema.edges
                if n.id in (e.source, e.target)
            }
        )
        if not ends:
            raise ValueError(
                f"Heat source '{n.id}' is not wired: connect it by one edge to the "
                "plenum, momentum chamber or combustor it heats."
            )
        if len(ends) > 1:
            raise ValueError(
                f"Heat source '{n.id}' is wired to {len(ends)} nodes ({', '.join(ends)}): "
                "a heat source heats one node. Use one heat source per node."
            )
        target = ends[0]
        if kinds[target] not in HEAT_SOURCE_TARGETS:
            raise ValueError(
                f"Heat source '{n.id}' is wired to '{target}', a {kinds[target].replace('_', ' ')}: "
                "a heat source heats a plenum, momentum chamber or combustor (a node whose "
                "outlet temperature comes from an energy balance), not an element or a boundary."
            )
        links[n.id] = target
    return links


def _heat_source_boundary(node_id: str, data: HeatSourceData) -> EnergyBoundary:
    """The EnergyBoundary a heat source puts on its node, by its mode."""
    if data.mode == "Tt_target":
        return EnergyBoundary(node_id, Tt_target=data.Tt_target)
    if data.mode == "fraction":
        return EnergyBoundary(node_id, fraction=data.fraction)
    return EnergyBoundary(node_id, Q=data.Q)


def build_pressure_loss(data, area: float = 0.1):
    """Converts PressureLossData schema to a pressure-loss callable."""
    if data is None:
        return None
    t = getattr(data, "type", None)
    if isinstance(data, ConstantFractionLossData) or t == "constant_fraction":
        return ConstantFractionLoss(xi=data.xi)
    if isinstance(data, LinearThetaFractionLossData) or t == "linear_theta_fraction":
        return LinearThetaFractionLoss(k=data.k, xi0=data.xi0)
    if isinstance(data, ConstantHeadLossData) or t == "constant_head":
        return ConstantHeadLoss(zeta=data.zeta, area=area)
    if isinstance(data, LinearThetaHeadLossData) or t == "linear_theta_head":
        return LinearThetaHeadLoss(k=data.k, zeta0=data.zeta0, area=area)
    return None


def map_surface_model(data):
    """Maps UI SurfaceModelData to combaero.network models."""
    if isinstance(data, SmoothModelData) or data.type == "smooth":
        return SmoothModel()
    if isinstance(data, RibbedModelData) or data.type == "ribbed":
        # getattr with the schema's own defaults: a network saved before a
        # field existed should load with that field defaulted, not raise an
        # AttributeError the user cannot act on.
        d = RibbedModelData()
        return RibbedModel(
            e_D=getattr(data, "e_D", d.e_D),
            p_e=getattr(data, "p_e", d.p_e),
            alpha_deg=getattr(data, "alpha_deg", d.alpha_deg),
            W_H=getattr(data, "W_H", d.W_H),
            n_ribbed_walls=getattr(data, "n_ribbed_walls", d.n_ribbed_walls),
            smooth_wall_Nu_multiplier=getattr(
                data, "smooth_wall_Nu_multiplier", d.smooth_wall_Nu_multiplier
            ),
        )
    if isinstance(data, ImpingementModelData) or data.type == "impingement":
        d = ImpingementModelData()
        return ImpingementModel(
            d_jet=getattr(data, "d_jet", d.d_jet),
            xn_d=getattr(data, "xn_d", d.xn_d),
            yn_d=getattr(data, "yn_d", d.yn_d),
            z_d=getattr(data, "z_d", d.z_d),
            row=getattr(data, "row", d.row),
            C_D=getattr(data, "C_D", d.C_D),
        )
    if isinstance(data, SingleJetImpingementModelData) or data.type == "single_jet_impingement":
        import combaero as cb

        d = SingleJetImpingementModelData()
        bc_str = getattr(data, "bc", d.bc)
        bc = (
            cb.ImpingementThermalBC.ConstantWallTemperature
            if bc_str == "constant_wall_temperature"
            else cb.ImpingementThermalBC.ConstantHeatFlux
        )
        return SingleJetImpingementModel(
            bc=bc,
            d_jet=getattr(data, "d_jet", d.d_jet),
            L_D=getattr(data, "L_D", d.L_D),
            R_D=getattr(data, "R_D", d.R_D),
        )
    if isinstance(data, PinFinModelData) or data.type == "pin_fin":
        import combaero as cb

        d = PinFinModelData()
        D = getattr(data, "pin_diameter", d.pin_diameter)
        H = getattr(data, "channel_height", d.channel_height)
        return PinFinModel(
            pin_diameter=D,
            H_D=H / D if D else d.channel_height / d.pin_diameter,
            S_D=getattr(data, "S_D", d.S_D),
            X_D=getattr(data, "X_D", d.X_D),
            N_rows=getattr(data, "N_rows", d.N_rows),
            arrangement=(
                cb.PinArrangement.Staggered
                if getattr(data, "is_staggered", d.is_staggered)
                else cb.PinArrangement.Inline
            ),
            k_pin=getattr(data, "k_pin", d.k_pin),
        )

    # Dimpled was removed in 0.7.0. Reject rather than silently
    # falling back to smooth: a saved network asking for one of these and
    # getting an unenhanced surface back is a wrong answer presented as a
    # right one.
    raise ValueError(removed_surface_message(getattr(data, "type", "unknown")))


def resolve_composition(
    comp: CompositionData, T: float, P: float, solver_settings=None
) -> list[float]:
    """
    Converts UI CompositionData into a mass fraction list (Y)
    using combaero.species utilities.
    """
    import combaero as cb

    if comp.source == "dry_air":
        Y = cb.species.dry_air_mass()
    elif comp.source == "humid_air":
        s = solver_settings
        ambient_T = (
            comp.ambient_T
            if comp.ambient_T is not None
            else (s.ambient_T if s is not None and s.ambient_T is not None else 288.15)
        )
        ambient_P = (
            comp.ambient_P
            if comp.ambient_P is not None
            else (s.ambient_P if s is not None and s.ambient_P is not None else 101325.0)
        )
        ambient_RH = (
            comp.relative_humidity
            if comp.relative_humidity is not None
            else (s.ambient_RH if s is not None and s.ambient_RH is not None else 0.6)
        )
        Y = cb.species.humid_air_mass(ambient_T, ambient_P, ambient_RH)
    elif comp.source == "fuel":
        # Default fuel to CH4 for now, or use a specific fuel if provided
        Y = cb.species.pure_species("CH4")
        Y = cb.species.to_mass(Y)
    elif comp.source == "custom" and comp.custom_fractions:
        vec = cb.species.from_mapping(comp.custom_fractions)
        Y = cb.species.to_mass(vec) if comp.mode == "mole" else vec
    else:
        Y = cb.species.dry_air_mass()

    return [float(y) for y in Y]


def _expand_initial_guess(short_guess: dict, prefix: str) -> dict:
    """Translate short GUI guess keys (``P``, ``T``, ``m_dot``) into the fully
    qualified solver unknown names (e.g. ``"plenum1.P"``). ``P`` and ``T`` are
    broadcast to both static and total variants because the solver exposes
    them as separate unknowns on some node types; the solver silently ignores
    entries that do not correspond to an actual unknown.
    Any key already containing a ``.`` is passed through unchanged so that
    hand-authored qualified guesses still work. Explicit short keys
    (e.g. ``Pt``) take precedence over broadcast values from ``P``/``T``.
    """
    if not short_guess:
        return {}
    result: dict = {}
    # Pass 1 -- broadcast P and T to both static and total variants.
    for key, val in short_guess.items():
        if "." in key:
            continue
        if key == "P":
            result[f"{prefix}.P"] = val
            result[f"{prefix}.Pt"] = val
        elif key == "T":
            result[f"{prefix}.T"] = val
            result[f"{prefix}.Tt"] = val
    # Pass 2 -- explicit short keys (Pt, Tt, m_dot, X[i], etc.)
    # are applied after the broadcast so they win on conflict.
    for key, val in short_guess.items():
        if "." in key:
            result[key] = val
            continue
        if key in ("P", "T"):
            continue  # already handled in pass 1
        result[f"{prefix}.{key}"] = val
    return result


def _guess_from_prior_result(node_data: dict, prefix: str) -> dict:
    """Build an initial-guess dict from a prior solve result stored in node_data.

    The GUI frontend writes the last solve result into node_data["result"].
    When data.initial_guess is empty (no explicit user override), this
    provides a warm-start from the previously found state so that repeated
    solves on a difficult network benefit from accumulated Newton progress.
    Works for both pressure-node results (state.P/T/Pt/Tt) and element
    results (m_dot stored at the top level).

    Only CONVERGED prior results are used: the best iterate of a failed
    solve is stored too, and seeding from it parks Newton at the stall
    point -- one failed solve (e.g. a timeout) would otherwise poison
    every subsequent run of a network that converges fine from cold.
    """
    prior = node_data.get("result")
    if not prior or not prior.get("success"):
        return {}
    short: dict[str, float] = {}
    state = prior.get("state")
    if state:
        if (p := state.get("P")) is not None:
            short["P"] = float(p)
        if (t := state.get("T")) is not None:
            short["T"] = float(t)
        if (pt := state.get("Pt")) is not None:
            short["Pt"] = float(pt)
        if (tt := state.get("Tt")) is not None:
            short["Tt"] = float(tt)
    if (m_dot := prior.get("m_dot")) is not None:
        short["m_dot"] = float(m_dot)
    return _expand_initial_guess(short, prefix)


def _build_discrete_loss_correlation(
    corr_type: str,
    xi: float,
    k: float,
    xi0: float,
    zeta: float,
    zeta0: float,
    area: float,
    m_ref: float | None = None,
):
    """Build a pressure loss correlation object from discrete loss params."""
    if corr_type == "constant_fraction":
        return ConstantFractionLoss(xi=xi, m_ref=m_ref)
    if corr_type == "linear_theta_fraction":
        return LinearThetaFractionLoss(k=k, xi0=xi0, m_ref=m_ref)
    if corr_type == "constant_head":
        return ConstantHeadLoss(zeta=zeta, area=area)
    if corr_type == "linear_theta_head":
        return LinearThetaHeadLoss(k=k, zeta0=zeta0, area=area)
    return ConstantFractionLoss(xi=xi, m_ref=m_ref)


def _channel_D_at_tee_port(
    schema_edges,
    schema_nodes,
    tee_id: str,
    port: str,
) -> float | None:
    """Return the inner diameter [m] of a channel directly connected to a tee port.

    The user-facing schema connects channels straight to tee handles (the
    port-MCN is auto-inserted at network-build time, not at schema level).
    Used by ``mpce_tee`` dispatch to inherit ``F_C`` / ``F_branch``
    from the connected channel geometry when the
    user has not explicitly set them. Returns None when no channel is
    directly connected (e.g., plenum, boundary, multi-hop).
    """
    hit = _channel_at_tee_port(schema_edges, schema_nodes, tee_id, port)
    return None if hit is None else hit[1]


def _channel_at_tee_port(
    schema_edges,
    schema_nodes,
    tee_id: str,
    port: str,
):
    """The edge joining a tee port straight to a channel, and that channel's D [m].

    None when no channel with a positive D is directly connected.
    """
    src_handle = f"port-{port}-source"
    tgt_handle = f"port-{port}-target"
    for edge in schema_edges:
        if edge.data and edge.data.get("type") == "thermal":
            continue
        neighbour_id: str | None = None
        if edge.source == tee_id and edge.sourceHandle == src_handle:
            neighbour_id = edge.target
        elif edge.target == tee_id and edge.targetHandle == tgt_handle:
            neighbour_id = edge.source
        if neighbour_id is None:
            continue
        neighbour = next((n for n in schema_nodes if n.id == neighbour_id), None)
        if neighbour is None or neighbour.type != "channel":
            continue
        d = neighbour.data.get("D") if isinstance(neighbour.data, dict) else None
        if d is None:
            continue
        try:
            d_val = float(d)
        except (TypeError, ValueError):
            continue
        if d_val > 0.0:
            return edge, d_val
    return None


def _channel_id_at_tee_port(edge, tee_id: str) -> str:
    """The channel at the far end of a tee-port edge."""
    return edge.target if edge.source == tee_id else edge.source


#: Common-arm area [m^2] of a tee with no channel on its common or straight
#: port and no user F_C (reported as "default"; the GUI mirrors it).
DEFAULT_TEE_F_C = 0.01


def _resolve_tee_areas_with_sources(
    schema_edges,
    schema_nodes,
    tee_id: str,
    F_C_explicit: float | None,
    F_branch_explicit: float | None,
    psi_fallback: float,
) -> tuple[float, float, float, str, str]:
    """Resolve (F_C, F_branch, psi) for a 3-port tee, with each area's source.

    Precedence:
      1. Explicit user-set value (non-null) wins.            -> "user"
      2. Channel D at the connected port (auto-inherited).   -> "channel '<id>'"
      3. Fallback: F_C=0.01 m^2                              -> "default";
         F_branch derived from psi_fallback                  -> "derived from F_C / psi"

    The tee's areas are resolved HERE, in the builder, and handed to the
    element as explicit port areas -- so the element itself would report them
    all as "user". The sources returned here are the real provenance, attached
    to the element for ``resolve_preview`` (#559).

    Returns (F_C, F_branch, psi, F_C_source, F_branch_source).
    """
    F_C = F_C_explicit
    F_C_source = "user" if F_C is not None else ""
    if F_C is None:
        # Try common port first, then straight (both are F_C in Bassett/Idelchik
        # geometry where F_s = F_c by definition).
        for port in ("common", "straight"):
            hit = _channel_at_tee_port(schema_edges, schema_nodes, tee_id, port)
            if hit is not None:
                F_C = _math.pi / 4.0 * hit[1] * hit[1]
                F_C_source = f"channel '{_channel_id_at_tee_port(hit[0], tee_id)}'"
                break
    if F_C is None:
        F_C = DEFAULT_TEE_F_C
        F_C_source = "default"

    F_branch = F_branch_explicit
    F_branch_source = "user" if F_branch is not None else ""
    if F_branch is None:
        hit = _channel_at_tee_port(schema_edges, schema_nodes, tee_id, "branch")
        if hit is not None:
            F_branch = _math.pi / 4.0 * hit[1] * hit[1]
            F_branch_source = f"channel '{_channel_id_at_tee_port(hit[0], tee_id)}'"

    if F_branch is not None and F_branch > 0.0:
        psi = F_C / F_branch
    else:
        # No branch channel; use stored psi to derive F_branch
        psi = max(float(psi_fallback), 1e-12)
        F_branch = F_C / psi
        F_branch_source = "derived from F_C / psi"

    return F_C, F_branch, psi, F_C_source, F_branch_source


def _resolve_tee_areas(
    schema_edges,
    schema_nodes,
    tee_id: str,
    F_C_explicit: float | None,
    F_branch_explicit: float | None,
    psi_fallback: float,
) -> tuple[float, float, float]:
    """Resolve (F_C, F_branch, psi) for a 3-port tee, applying inheritance.

    See ``_resolve_tee_areas_with_sources`` for the precedence; this returns
    just the three resolved floats.
    """
    F_C, F_branch, psi, _, _ = _resolve_tee_areas_with_sources(
        schema_edges, schema_nodes, tee_id, F_C_explicit, F_branch_explicit, psi_fallback
    )
    return F_C, F_branch, psi


def _port_from_handle(handle: str | None) -> str | None:
    """Extract port name from a tee handle string.

    'port-common-source' -> 'common', 'port-branch-target' -> 'branch', etc.
    Returns None for non-tee handles.
    """
    if (
        handle
        and handle.startswith("port-")
        and (handle.endswith("-source") or handle.endswith("-target"))
    ):
        return handle.split("-")[1]
    return None


def _find_junction_port(
    edges: list,
    junction_id: str,
    port: str,
    nodes_map: dict,
    tee_port_node_map: dict | None = None,
    wall_edge_remap: dict | None = None,
    all_ports: tuple[str, ...] = ("common", "straight", "branch"),
) -> str:
    """Return the physical-node ID connected to a named junction port.

    Shared by the tee (ports common/straight/branch) and the ejector (ports
    primary/secondary/outlet); ``all_ports`` only shapes the error message.
    Checks both edge directions: junction as source (sourceHandle ==
    port-X-source) and junction as target (targetHandle == port-X-target).
    Falls back to *tee_port_node_map* for auto-created junction nodes when an
    element is connected directly to a port without an explicit plenum.
    """
    src_handle = f"port-{port}-source"
    tgt_handle = f"port-{port}-target"
    for edge in edges:
        if edge.data and edge.data.get("type") == "thermal":
            continue
        if edge.source == junction_id and edge.sourceHandle == src_handle:
            nid = edge.target
        elif edge.target == junction_id and edge.targetHandle == tgt_handle:
            nid = edge.source
        else:
            continue
        # An auto-inserted port node wins: it is also created for a plenum or
        # boundary wired straight to the port (build pass 1b).
        if tee_port_node_map and (junction_id, port) in tee_port_node_map:
            return tee_port_node_map[(junction_id, port)]
        if nid in nodes_map:
            if wall_edge_remap and edge.id in wall_edge_remap:
                return wall_edge_remap[edge.id]
            return nid
        raise ValueError(
            f"Junction '{junction_id}' port '{port}' connects to '{nid}', which is "
            "not a physical node (plenum or boundary). All junction ports must "
            "connect directly to plena or boundaries."
        )
    raise ValueError(
        f"Junction '{junction_id}' port '{port}' is not connected. All ports "
        f"({', '.join(all_ports)}) must be wired to a node."
    )


def _wetted_area_from_wall(
    net: FlowNetwork, side: str, side_area: float | None, contact_area: float | None
) -> None:
    """Give a walled combustor, momentum chamber or discrete loss its wetted area.

    The GUI has no wetted-area field on these components. Their surface used
    to be the flow area (passed here when typed, written at resolve when
    inherited); a wetted area is not a flow area, and neither the library nor
    the GUI derives one from the other any more (#1022). A wall drawn to the
    component names the area that exchanges heat -- its side's face, else the
    wall's area, which the coupling uses either way -- so that is the
    surface's area. An unwalled one has no convective surface.
    """
    obj = net.nodes.get(side) or net.elements.get(side)
    if not isinstance(obj, CombustorNode | MomentumChamberNode | PressureLossElement):
        return
    area = side_area if side_area is not None else contact_area
    if obj.surface.area == 0.0 and area:
        obj.surface.area = float(area)


def build_network_from_schema(schema: NetworkGraphSchema) -> FlowNetwork:
    checked = _checked_edges(schema)  # refuses an edge to a missing node first
    heat_links = _heat_source_links(schema)
    schema = checked
    net = FlowNetwork()
    nodes_map = {}  # ID -> Physical Node
    plenum_node_ids: set[str] = set()
    element_nodes = []  # List of node schemas that are actually elements

    # 1. First Pass: Create Physical Nodes
    for node_schema in schema.nodes:
        node_id = node_schema.id
        node_data = node_schema.data
        node_type = node_schema.type

        if node_type == "plenum":
            data = _parse(PlenumData, node_data, node_id)
            node = PlenumNode(node_id)
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            if not node.initial_guess:
                node.initial_guess = _guess_from_prior_result(node_data, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
            plenum_node_ids.add(node_id)
        elif node_type == "mass_boundary":
            data = _parse(MassBoundaryData, node_data, node_id)
            Y = resolve_composition(data.composition, data.Tt, 101325.0, schema.solver_settings)
            node = MassFlowBoundary(node_id, m_dot=data.m_dot, Tt=data.Tt, Y=Y)
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            if not node.initial_guess:
                node.initial_guess = _guess_from_prior_result(node_data, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "fuel_boundary":
            data = _parse(FuelBoundaryData, node_data, node_id)
            if not data.combustor_id:
                raise ValueError(
                    f"Fuel boundary '{node_id}': no combustor to control. Connect it to "
                    "one, or drag its link handle onto the combustor it sets."
                )
            Y = resolve_composition(data.composition, data.Tt, 101325.0, schema.solver_settings)
            node = FuelBoundary(
                node_id,
                combustor=data.combustor_id,
                Tt_target=data.Tt_target,
                branch=data.branch,
                Tt=data.Tt,
                Y=Y,
                m_dot_guess=data.m_dot_guess,
            )
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "pressure_boundary":
            data = _parse(PressureBoundaryData, node_data, node_id)
            Y = resolve_composition(data.composition, data.Tt, data.Pt, schema.solver_settings)
            node = PressureBoundary(node_id, Pt=data.Pt, Tt=data.Tt, Y=Y)
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "combustor":
            data = _parse(CombustorData, node_data, node_id)
            node = CombustorNode(
                node_id,
                method=data.method,
                area=data.area,
                Dh=data.Dh,
                surface=ConvectiveSurface(
                    # No wetted-area field: a wall drawn to it gives one (#1022).
                    area=0.0,
                    model=map_surface_model(data.surface),
                    Nu_multiplier=data.Nu_multiplier
                    * (schema.solver_settings.Nu_multiplier or 1.0),
                    f_multiplier=data.f_multiplier * (schema.solver_settings.f_multiplier or 1.0),
                ),
            )
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            if not node.initial_guess:
                node.initial_guess = _guess_from_prior_result(node_data, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "momentum_chamber":
            data = _parse(MomentumChamberData, node_data, node_id)
            node = MomentumChamberNode(
                node_id,
                area=data.area,
                Dh=data.Dh,
                surface=ConvectiveSurface(
                    # No wetted-area field: a wall drawn to it gives one (#1022).
                    area=0.0,
                    model=map_surface_model(data.surface),
                    Nu_multiplier=data.Nu_multiplier
                    * (schema.solver_settings.Nu_multiplier or 1.0),
                    f_multiplier=data.f_multiplier * (schema.solver_settings.f_multiplier or 1.0),
                ),
            )
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            if not node.initial_guess:
                node.initial_guess = _guess_from_prior_result(node_data, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "probe":
            pass  # Display-only diagnostic nodes -- not part of the solver graph
        elif node_type == HEAT_SOURCE_KIND:
            pass  # An EnergyBoundary on the node it is wired to (below)
        elif node_type == "wall":
            data = _parse(WallData, node_data, node_id)
            node = WallNode(node_id)
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type in ELEMENT_KINDS:
            element_nodes.append(node_schema)
        else:
            raise ValueError(
                f"Node '{node_id}' is of unknown kind {node_type!r}; this version "
                "cannot build it, and solving without it would be a different network."
            )

    # Heat sources (#792): each an EnergyBoundary, named by the source's id, on
    # the node it is wired to. Two Tt_target sources on one node are refused
    # by the node itself (one outlet temperature).
    heat_sources: dict[str, tuple[str, EnergyBoundary]] = {}
    hs_data = {n.id: n.data for n in schema.nodes if n.type == HEAT_SOURCE_KIND}
    for hs_id, target in heat_links.items():
        data = _parse(HeatSourceData, hs_data[hs_id], hs_id)
        eb = _heat_source_boundary(hs_id, data)
        try:
            nodes_map[target].add_energy_boundary(eb)
        except ValueError as exc:
            raise ValueError(f"Heat source '{hs_id}': {exc}") from exc
        heat_sources[hs_id] = (target, eb)
    net._gui_heat_sources = heat_sources

    element_ids = {node.id for node in element_nodes}
    # Impingement arrays (#465) expand into a row chain. Their ends fan out to
    # N plates (supply) or collect the last crossflow segment (exit), which a
    # MomentumChamberNode junction cannot represent, so element-to-array links
    # get a PlenumNode instead -- the jet supply plenum the source assumes.
    # A link on an array's bypass handle gets a PlenumNode too: the lead-in
    # accelerates the bypass from that node's total pressure, which a feed
    # channel keeps into a plenum, and an upstream array's last gap segment
    # continues through it with its momentum (#467).
    array_ids = {node.id for node in element_nodes if node.type == "impingement_array"}
    gui_arrays: dict[str, ImpingementArray] = {}
    # Junction elements with named handle ports (common/straight/branch for
    # the tee; primary/secondary/outlet for the ejector). Both get the same
    # per-port MomentumChamberNode auto-insertion below.
    junction_ids = {
        node.id
        for node in element_nodes
        if node.type in ("mpce_tee", "ejector", "effusion_liner", "flow_splitter")
    }
    # An effusion liner's ports are a reservoir, a dump and a discharge that
    # N panels merge into: plenum semantics. A momentum chamber there would
    # refuse the N inflows and mis-state the reservoir (#471).
    liner_ids = {node.id for node in element_nodes if node.type == "effusion_liner"}
    # A flow splitter's (#903) inlet is the node its N outlets draw from, and
    # each outlet a valve discharging into a node at rest: plenum ports too. A
    # momentum chamber has one bulk velocity and refuses a split.
    splitter_ids = {node.id for node in element_nodes if node.type == "flow_splitter"}
    plenum_port_ids = liner_ids | splitter_ids

    # 1b. Create implicit junction nodes for element -> element visual links.
    # For direct element<->element connections we auto-insert a MomentumChamberNode.
    # For direct element<->tee-port connections we auto-insert a PlenumNode so that
    # users don't need to manually place a plenum on every tee arm.
    edge_junction_map: dict[tuple[str, str], str] = {}
    _ejector_seeds: list[tuple[EjectorElement, tuple[str, str, str]]] = []
    tee_port_node_map: dict[tuple[str, str], str] = {}  # (junction_id, port) -> node_id
    # (from_node, to_node) lossless links joining a directly-wired plenum or
    # boundary to the port node inserted for it.
    port_node_links: list[tuple[str, str]] = []
    for edge in schema.edges:
        if edge.data and edge.data.get("type") == "thermal":
            continue
        src_elem = edge.source in element_ids
        tgt_elem = edge.target in element_ids
        src_tee = edge.source in junction_ids
        tgt_tee = edge.target in junction_ids

        if src_elem and tgt_elem and not src_tee and not tgt_tee:
            # Standard element-to-element: auto-insert MomentumChamberNode
            junction_id = f"__junction__{edge.source}__{edge.target}"
            if junction_id not in nodes_map:
                to_array = edge.source in array_ids or edge.target in array_ids
                junction_node = (
                    PlenumNode(junction_id) if to_array else MomentumChamberNode(junction_id)
                )
                net.add_node(junction_node)
                nodes_map[junction_id] = junction_node
            edge_junction_map[(edge.source, edge.target)] = junction_id

        elif src_tee and tgt_elem:
            # Junction port -> element (e.g. tee.common-source -> channel): auto-insert
            # MomentumChamberNode so that Pt != P_static and velocity is accounted for.
            port = _port_from_handle(edge.sourceHandle)
            if port and (edge.source, port) not in tee_port_node_map:
                jid = f"__tee_jct__{edge.source}_{port}"
                if jid not in nodes_map:
                    jn = (
                        PlenumNode(jid)
                        if edge.source in plenum_port_ids
                        else MomentumChamberNode(jid)
                    )
                    net.add_node(jn)
                    nodes_map[jid] = jn
                tee_port_node_map[(edge.source, port)] = jid
                edge_junction_map[(edge.source, edge.target)] = jid

        elif src_elem and tgt_tee:
            # Element -> junction port (e.g. channel -> tee.straight-target): auto-insert MomentumChamberNode
            port = _port_from_handle(edge.targetHandle)
            if port and (edge.target, port) not in tee_port_node_map:
                jid = f"__tee_jct__{edge.target}_{port}"
                if jid not in nodes_map:
                    jn = (
                        PlenumNode(jid)
                        if edge.target in plenum_port_ids
                        else MomentumChamberNode(jid)
                    )
                    net.add_node(jn)
                    nodes_map[jid] = jn
                tee_port_node_map[(edge.target, port)] = jid
                edge_junction_map[(edge.source, edge.target)] = jid

        elif (src_tee and edge.target in nodes_map) or (tgt_tee and edge.source in nodes_map):
            # Junction port <-> plenum / boundary wired directly. The junction
            # needs a MomentumChamberNode on each port (MultiPortChamberBase
            # refuses anything else), so insert one and join it to the node with
            # a lossless link, as a channel on the port would get (#559 audit).
            # A momentum chamber the user placed is used as-is; a liner's ports
            # are plenums by design.
            jct, other = (edge.source, edge.target) if src_tee else (edge.target, edge.source)
            handle = edge.sourceHandle if src_tee else edge.targetHandle
            port = _port_from_handle(handle)
            if jct in splitter_ids:
                # A plenum or boundary is used as it is. A momentum chamber
                # at the inlet would carry the split itself: it feeds an
                # inserted plenum through a lossless link instead.
                if (
                    port == "in"
                    and isinstance(nodes_map[other], MomentumChamberNode)
                    and (jct, port) not in tee_port_node_map
                ):
                    jid = f"__tee_jct__{jct}_{port}"
                    if jid not in nodes_map:
                        jn = PlenumNode(jid)
                        net.add_node(jn)
                        nodes_map[jid] = jn
                    tee_port_node_map[(jct, port)] = jid
                    port_node_links.append((other, jid))
            elif (
                port
                and jct not in liner_ids
                and not isinstance(nodes_map[other], MomentumChamberNode)
                and (jct, port) not in tee_port_node_map
            ):
                jid = f"__tee_jct__{jct}_{port}"
                if jid not in nodes_map:
                    jn = MomentumChamberNode(jid)
                    net.add_node(jn)
                    nodes_map[jid] = jn
                tee_port_node_map[(jct, port)] = jid
                # Flow direction follows the canvas edge: node -> port (inlet)
                # or port -> node (outlet).
                port_node_links.append((other, jid) if tgt_tee else (jid, other))

    # 1c. A tee port whose channel has a different area meets it through a
    # sharp step (#536). The tee's ports keep the closures' geometry (F_S =
    # F_C), so a straight channel narrower than the main pipe used to join a
    # port of the main area with nothing in between: the ideal tee, which
    # couples on static pressure, then raised the channel's Pt by its dynamic
    # head, and Mynard, which couples on Pt, contracted losslessly. Bolted
    # together, two bores of different size are a sharp contraction or
    # expansion, so that is what is inserted: a hidden AreaChangeElement
    # (Weisbach / Borda-Carnot) and a momentum chamber on the channel side.
    tee_steps: dict[str, tuple[str, str, float]] = {}
    for tee_node in element_nodes:
        if tee_node.type != "mpce_tee":
            continue
        _td = _parse(MPCETeeData, tee_node.data, tee_node.id)
        _F_C, _F_B, _ = _resolve_tee_areas(
            schema.edges, schema.nodes, tee_node.id, _td.F_C, _td.F_branch, _td.psi
        )
        for port, A_port in (("common", _F_C), ("straight", _F_C), ("branch", _F_B)):
            hit = _channel_at_tee_port(schema.edges, schema.nodes, tee_node.id, port)
            if hit is None:
                continue
            edge, D = hit
            A_ch = _math.pi / 4.0 * D * D
            if abs(A_ch / A_port - 1.0) <= TEE_PORT_STEP_TOL:
                continue
            key = (edge.source, edge.target)
            port_node = edge_junction_map.get(key)
            if port_node is None:
                continue
            step_id = f"__tee_step__{tee_node.id}_{port}"
            side_id = f"{step_id}__node"
            side = MomentumChamberNode(side_id)
            net.add_node(side)
            nodes_map[side_id] = side
            if edge.source == tee_node.id:  # drawn tee -> channel
                step = AreaChangeElement(
                    step_id, port_node, side_id, F0=A_port, F1=A_ch, model_type="sharp"
                )
            else:
                step = AreaChangeElement(
                    step_id, side_id, port_node, F0=A_ch, F1=A_port, model_type="sharp"
                )
            step.regime = schema.solver_settings.global_regime
            net.add_element(step)
            edge_junction_map[key] = side_id
            tee_steps[step_id] = (tee_node.id, port, A_ch / A_port)
    net._gui_tee_steps = tee_steps

    # Pre-scan: for each edge whose target is a WallNode, the 2nd+ connection
    # gets its own clone so the solver enforces m_dot=0 per branch individually
    # (not just sum(m_dot)=0, which allows non-zero flows to cancel).
    # This covers both regular elements and tee junction port connections.
    _wall_edge_remap: dict[str, str] = {}  # edge_id -> effective target node_id
    _wall_target_count: dict[str, int] = {}
    for _edge in schema.edges:
        if _edge.data and _edge.data.get("type") == "thermal":
            continue
        _tgt = _edge.target
        if isinstance(nodes_map.get(_tgt), WallNode):
            _wall_target_count[_tgt] = _wall_target_count.get(_tgt, 0) + 1
            if _wall_target_count[_tgt] > 1:
                _n = _wall_target_count[_tgt]
                _clone_id = f"{_tgt}_w{_n}"
                _clone = WallNode(_clone_id)
                net.add_node(_clone)
                nodes_map[_clone_id] = _clone
                _wall_edge_remap[_edge.id] = _clone_id

    # 2. Second Pass: Create Elements
    for elem_node_schema in element_nodes:
        elem_id = elem_node_schema.id
        elem_data = elem_node_schema.data
        elem_type = elem_node_schema.type

        # 3-port element (Mynard-based). The flow_direction field constrains
        # which side is upstream.
        if elem_type == "mpce_tee":
            _md = _parse(MPCETeeData, elem_data, elem_id)
            # Resolve areas with inheritance precedence (see _resolve_tee_areas).
            _F_C, _A_branch, _psi, _F_C_src, _A_branch_src = _resolve_tee_areas_with_sources(
                schema.edges,
                schema.nodes,
                elem_id,
                F_C_explicit=_md.F_C,
                F_branch_explicit=_md.F_branch,
                psi_fallback=_md.psi,
            )

            _common = _find_junction_port(
                schema.edges, elem_id, "common", nodes_map, tee_port_node_map, _wall_edge_remap
            )
            _straight = _find_junction_port(
                schema.edges, elem_id, "straight", nodes_map, tee_port_node_map, _wall_edge_remap
            )
            _branch = _find_junction_port(
                schema.edges, elem_id, "branch", nodes_map, tee_port_node_map, _wall_edge_remap
            )

            # Branch (separating): single inlet at common, two outlets.
            # Merge (joining): two inlets at straight+branch, single outlet at common.
            if _md.flow_direction == "branch":
                _inlets = [_common]
                _outlets = [_straight, _branch]
                _inlet_angles = [0.0]
                _outlet_angles = [0.0, _md.theta_deg]
                _port_areas = [_F_C, _F_C, _A_branch]
            else:  # "merge"
                _inlets = [_straight, _branch]
                _outlets = [_common]
                _inlet_angles = [0.0, _md.theta_deg]
                _outlet_angles = [0.0]
                _port_areas = [_F_C, _A_branch, _F_C]

            # strict=False: strict mode raises whenever a Newton ITERATE
            # (not the converged answer) explores a wrong-sign port
            # m_dot, killing GUI solves that would otherwise converge
            # ("Unexpected error during residual evaluation"). Soft mode
            # lets the solver explore; the post-solve
            # verify_solution_consistent guard in NetworkSolver rejects
            # any wrong-direction artifact root at convergence.
            if _md.junction_model == "constant_k":
                # Simplest-model tier: fixed handbook K per non-common
                # port, referenced to the common-leg dynamic head. Port
                # order is inlets first: branch = [common, straight,
                # branch]; merge = [straight, branch, common].
                # A K left blank comes from Bassett's analytical pair at the
                # design split and the resolved geometry.
                _K_str_d, _K_bra_d = ConstantKTeeElement.bassett_design_K(
                    _md.flow_direction, _md.design_branch_share, _psi, _md.theta_deg
                )
                _K_str = _md.K_straight if _md.K_straight is not None else _K_str_d
                _K_bra = _md.K_branch if _md.K_branch is not None else _K_bra_d
                # A pair that creates total pressure at its own design split
                # converges and is then rejected by the second-law check; refuse
                # it here with the cause instead (#779). No clamp: a clamped
                # zero is still not the user's K.
                _net = ConstantKTeeElement.design_net_loss(_md.design_branch_share, _K_str, _K_bra)
                if _net < 0.0:
                    _blank = [
                        name
                        for name, typed in (
                            ("K_straight", _md.K_straight),
                            ("K_branch", _md.K_branch),
                        )
                        if typed is None
                    ]
                    _origin = (
                        f" ({' and '.join(_blank)} from Bassett 2001's analytical pair)"
                        if _blank
                        else ""
                    )
                    raise ValueError(
                        f"Tee '{elem_id}': constant-K pair K_straight = {_K_str:.3f}, "
                        f"K_branch = {_K_bra:.3f}{_origin} loses net total pressure "
                        f"{_net:.3f} < 0 at design branch share {_md.design_branch_share:g}: "
                        "the tee would create total pressure and the solve would be "
                        "rejected by the second-law check. Enter K_straight and K_branch "
                        "explicitly (for example measured values), move the design share, "
                        "or use the Mynard junction model."
                    )
                if _md.flow_direction == "branch":
                    _k_ports = {1: _K_str, 2: _K_bra}
                else:
                    _k_ports = {0: _K_str, 1: _K_bra}
                _mpce = ConstantKTeeElement(
                    id=elem_id,
                    inlet_nodes=_inlets,
                    outlet_nodes=_outlets,
                    inlet_angles_deg=_inlet_angles,
                    outlet_angles_deg=_outlet_angles,
                    port_areas=_port_areas,
                    flow_direction=_md.flow_direction,
                    strict=False,
                    K_ports=_k_ports,
                )
            elif _md.junction_model == "ideal_momentum":
                # Parameter-free baseline (#517): no K, no declared
                # direction -- the inlet/outlet split above only wires the
                # ports and sets the angle convention. loss_factor (#526) is
                # the user's knob; 1.0 is the closure itself.
                _mpce = IdealMomentumTeeElement(
                    id=elem_id,
                    inlet_nodes=_inlets,
                    outlet_nodes=_outlets,
                    inlet_angles_deg=_inlet_angles,
                    outlet_angles_deg=_outlet_angles,
                    port_areas=_port_areas,
                    loss_factor=_md.loss_factor,
                )
            else:
                _mpce = MultiPortChamberElement(
                    id=elem_id,
                    inlet_nodes=_inlets,
                    outlet_nodes=_outlets,
                    inlet_angles_deg=_inlet_angles,
                    outlet_angles_deg=_outlet_angles,
                    port_areas=_port_areas,
                    flow_direction=_md.flow_direction,
                    strict=False,
                    joining_etransfer_alpha=_md.joining_etransfer_alpha,
                )
            _guess = dict(_md.initial_guess)
            if _mpce.JUNCTION_UNKNOWN == "Pt_jct" and "P_jct" in _guess:
                # A file saved before #773 named the Mynard / constant-K
                # junction TOTAL "P_jct"; it is "Pt_jct" now. Migrated
                # silently, as the frontend's loader does.
                _guess.setdefault("Pt_jct", _guess["P_jct"])
                del _guess["P_jct"]
            _mpce.initial_guess = _expand_initial_guess(_guess, elem_id)
            if not _mpce.initial_guess:
                _mpce.initial_guess = _guess_from_prior_result(elem_data, elem_id)
            # The element sees these areas as explicit port areas ("user"); the
            # builder knows where they really came from. resolve_preview reports
            # these in place of the per-port sources (#559).
            _mpce._preview_fields = {
                "F_C": (_F_C, _F_C_src),
                "F_branch": (_A_branch, _A_branch_src),
            }
            net.add_element(_mpce)
            continue

        # Ejector: 3-port MultiPortChamberBase with fixed roles (primary +
        # secondary inlets, one outlet). Reuses the same named-port resolution
        # and per-port MomentumChamberNode auto-insertion as the tee.
        if elem_type == "ejector":
            _ed = _parse(EjectorData, elem_data, elem_id)
            _ports = ("primary", "secondary", "outlet")
            _primary = _find_junction_port(
                schema.edges,
                elem_id,
                "primary",
                nodes_map,
                tee_port_node_map,
                _wall_edge_remap,
                all_ports=_ports,
            )
            _secondary = _find_junction_port(
                schema.edges,
                elem_id,
                "secondary",
                nodes_map,
                tee_port_node_map,
                _wall_edge_remap,
                all_ports=_ports,
            )
            _outlet = _find_junction_port(
                schema.edges,
                elem_id,
                "outlet",
                nodes_map,
                tee_port_node_map,
                _wall_edge_remap,
                all_ports=_ports,
            )
            _ejector = EjectorElement(
                id=elem_id,
                primary_node=_primary,
                secondary_node=_secondary,
                outlet_node=_outlet,
                throat_area=_ed.throat_area,
                nozzle_exit_area=_ed.nozzle_exit_area,
                mixing_area=_ed.mixing_area,
                recovery_efficiency=_ed.recovery_efficiency,
                anchor_closure=_ed.anchor_closure,
                # The network's regime (#588): incompressible is Sanger's
                # jet pump, lossless apart from mixing.
                regime=schema.solver_settings.global_regime,
                k_nozzle=_ed.k_nozzle,
                k_suction=_ed.k_suction,
                k_throat=_ed.k_throat or 0.0,
                k_diffuser=_ed.k_diffuser or 0.0,
            )
            _ejector.initial_guess = _expand_initial_guess(_ed.initial_guess, elem_id)
            if not _ejector.initial_guess:
                _ejector.initial_guess = _guess_from_prior_result(elem_data, elem_id)
            net.add_element(_ejector)
            # Warm-started once the whole network exists (end of the build):
            # the flow seed walks each port's element chain, and an inherited
            # A_3 needs the outlet pipe, which may be built after the ejector.
            _ejector_seeds.append((_ejector, (_primary, _secondary, _outlet)))
            continue

        # Effusion liner (#471): three named ports, expanded into N stations.
        if elem_type == "effusion_liner":
            _ld = _parse(EffusionLinerData, elem_data, elem_id)
            _lports = ("coolantin", "coolantout", "discharge")
            _ends = [
                _find_junction_port(
                    schema.edges,
                    elem_id,
                    port,
                    nodes_map,
                    tee_port_node_map,
                    _wall_edge_remap,
                    all_ports=_lports,
                )
                for port in _lports
            ]
            liner = EffusionLiner(
                elem_id,
                n_segments=_ld.n_segments,
                length=_ld.length,
                width=_ld.width,
                duct_height=_ld.duct_height,
                hole_diameter=_ld.hole_diameter,
                wall_thickness=_ld.wall_thickness,
                pitch_x=_ld.pitch_x,
                pitch_y=_ld.pitch_y,
                angle_deg=_ld.angle_deg,
                entry_K=_ld.entry_K,
                roughness=_ld.roughness,
                wall_conductivity=_ld.wall_conductivity,
                gas_film=_ld.gas_film,
                gas_augmentation=_ld.gas_augmentation,
                turbulence_intensity=_ld.turbulence_intensity,
                gas_heat_flux=_ld.gas_heat_flux,
                film_direction=_ld.film_direction,
                material=_ld.material or None,
            )
            liner.add_to(net, _ends[0], _ends[1], _ends[2])
            gui_arrays[elem_id] = liner
            continue

        # Flow splitter (#903): one inlet and N outlet ports, expanded into one
        # SplitterOutletElement per outlet.
        if elem_type == "flow_splitter":
            _sd = _parse(FlowSplitterData, elem_data, elem_id)
            _n = len(_sd.fractions)
            _sports = ("in", *(f"out{i}" for i in range(_n)))
            _ends = [
                _find_junction_port(
                    schema.edges,
                    elem_id,
                    port,
                    nodes_map,
                    tee_port_node_map,
                    _wall_edge_remap,
                    all_ports=_sports,
                )
                for port in _sports
            ]
            _inlet, _outs = _ends[0], _ends[1:]
            _dup = sorted({o for o in _outs if _outs.count(o) > 1})
            if _dup:
                raise ValueError(
                    f"Flow splitter '{elem_id}': outlets {', '.join(_dup)} are wired to "
                    "the same node; give each outlet its own node."
                )
            _ref = _sd.reference
            try:
                splitter = FlowSplitter(
                    elem_id,
                    _inlet,
                    dict(zip(_outs, _sd.fractions, strict=True)),
                    reference=_outs[_ref] if _ref is not None else None,
                    areas={
                        o: a for o, a in zip(_outs, _sd.areas or [], strict=False) if a is not None
                    },
                    auto_reference=_sd.auto_reference,
                )
            except ValueError as exc:
                raise ValueError(
                    str(exc).replace(f"FlowSplitter {elem_id!r}", f"Flow splitter '{elem_id}'")
                ) from exc
            splitter.add_to(net)
            gui_arrays[elem_id] = GuiFlowSplitter(splitter)
            continue

        source_id = None
        target_id = None
        bypass_id = None
        incoming_count = 0
        outgoing_count = 0

        for edge in schema.edges:
            if edge.data and edge.data.get("type") == "thermal":
                continue
            if (
                edge.target == elem_id
                and elem_type == "impingement_array"
                and edge.targetHandle == ARRAY_BYPASS_HANDLE
            ):
                # Initial (bypass) crossflow into the gap's upstream end
                # (#467): declared by wiring this handle. Unwired, the end is
                # closed (#465).
                if bypass_id is not None:
                    raise ValueError(f"Impingement array '{elem_id}' has two bypass links.")
                if edge.source in nodes_map:
                    bypass_id = edge.source
                elif edge.source in element_ids:
                    bypass_id = edge_junction_map.get((edge.source, elem_id))
                continue
            if edge.target == elem_id:
                incoming_count += 1
                if incoming_count > 1:
                    raise ValueError(f"Element '{elem_id}' has multiple upstream links.")
                if edge.source in nodes_map:
                    source_id = edge.source
                elif edge.source in element_ids:
                    source_id = edge_junction_map.get((edge.source, elem_id))
            if edge.source == elem_id:
                outgoing_count += 1
                if outgoing_count > 1:
                    raise ValueError(f"Element '{elem_id}' has multiple downstream links.")
                if edge.target in nodes_map:
                    target_id = _wall_edge_remap.get(edge.id, edge.target)
                elif edge.target in element_ids:
                    target_id = edge_junction_map.get((elem_id, edge.target))

        if not source_id or not target_id:
            raise ValueError(f"Element '{elem_id}' is not fully connected.")

        if elem_type == "channel":
            data = _parse(ChannelData, elem_data, elem_id)
            # The convective area is left to the element: it depends on the
            # surface model (#462) and, when D is inherited, on the resolved
            # diameter. ChannelElement.resolve_topology fills it.
            conv_area = 0.0
            elem = ChannelElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                length=data.L,
                diameter=data.D,
                Dh=data.Dh,
                # Unset D inherits; with nothing to inherit from, a marked
                # default keeps a freshly dropped chain solvable (#570).
                default_diameter=DEFAULT_CHANNEL_DIAMETER,
                roughness=data.roughness,
                friction_model=data.friction_model,
                surface=ConvectiveSurface(
                    area=conv_area,
                    model=map_surface_model(data.surface),
                    Nu_multiplier=data.Nu_multiplier
                    * (schema.solver_settings.Nu_multiplier or 1.0),
                    f_multiplier=data.f_multiplier * (schema.solver_settings.f_multiplier or 1.0),
                ),
            )
            regime = (
                data.regime if data.regime != "default" else schema.solver_settings.global_regime
            )
            elem.regime = regime
            elem.initial_guess = _expand_initial_guess(data.initial_guess, elem_id)
            if not elem.initial_guess:
                elem.initial_guess = _guess_from_prior_result(elem_data, elem_id)
            net.add_element(elem)
        elif elem_type == "orifice":
            data = _parse(OrificeData, elem_data, elem_id)
            elem = OrificeElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                diameter=data.diameter,
                Cd=_fixed_cd(data),
                correlation=data.correlation,
                plate_thickness=data.plate_thickness,
                edge_radius=data.edge_radius,
                bevel_depth=data.bevel_depth,
                tap=data.tap,
                in_pipe=data.in_pipe,
                # "default": the element selects at build (#861).
                compressible_model=(
                    None if data.compressible_model == "default" else data.compressible_model
                ),
                # Cv is Jobson's alone; a value left from an earlier pick of
                # that model must not refuse the others.
                Cv=data.Cv if data.compressible_model == "jobson1955" else 1.0,
            )
            regime = (
                data.regime if data.regime != "default" else schema.solver_settings.global_regime
            )
            elem.regime = regime
            elem.initial_guess = _expand_initial_guess(data.initial_guess, elem_id)
            if not elem.initial_guess:
                elem.initial_guess = _guess_from_prior_result(elem_data, elem_id)
            net.add_element(elem)
        elif elem_type == "effusion_plate":
            data = _parse(EffusionPlateData, elem_data, elem_id)
            elem = EffusionPlateElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                hole_diameter=data.hole_diameter,
                wall_thickness=data.wall_thickness,
                pitch_x=data.pitch_x,
                pitch_y=data.pitch_y,
                panel_length=data.panel_length,
                panel_width=data.panel_width,
                angle_deg=data.angle_deg,
                correlation=data.correlation,
                Cd=_fixed_cd(data),
                internal_Nu_multiplier=data.internal_Nu_multiplier,
                wall_conductivity=data.wall_conductivity,
                gas_film=data.gas_film,
                gas_augmentation=data.gas_augmentation,
                turbulence_intensity=data.turbulence_intensity,
                gas_heat_flux=data.gas_heat_flux,
                material=data.material or None,
                gas_baseline=data.gas_baseline,
                gas_start_length=(
                    data.gas_start_length if data.gas_baseline == "flat_plate" else None
                ),
                cd_inclination=_cd_inclination(data),
            )
            net.add_element(elem)
        elif elem_type == "film_cooling":
            data = _parse(FilmCoolingData, elem_data, elem_id)
            shaped = data.hole_type == "shaped"
            elem = FilmCoolingElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                hole_diameter=data.hole_diameter,
                holes_per_row=data.holes_per_row,
                lateral_pitch=data.lateral_pitch,
                wall_thickness=data.wall_thickness,
                n_rows=data.n_rows,
                row_pitch=data.row_pitch if data.n_rows > 1 else None,
                angle_deg=data.angle_deg,
                hole_type=data.hole_type,
                area_ratio=data.area_ratio if shaped else None,
                breakout_width=data.breakout_width if shaped else None,
                turbulence_intensity=data.turbulence_intensity,
                correlation=data.correlation,
                Cd=_fixed_cd(data),
                cd_inclination=_cd_inclination(data),
            )
            net.add_element(elem)
        elif elem_type == "slot_film":
            data = _parse(SlotFilmData, elem_data, elem_id)
            elem = SlotFilmElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                slot_height=data.slot_height,
                span=data.span,
                lip_thickness=data.lip_thickness,
                turbulence_intensity=data.turbulence_intensity,
                C_m_multiplier=data.C_m_multiplier,
                correlation=data.correlation,
                Cd=_fixed_cd(data),
                film_htc=data.film_htc,
            )
            net.add_element(elem)
        elif elem_type == "impingement_array":
            data = _parse(ImpingementArrayData, elem_data, elem_id)
            arr = ImpingementArray(
                elem_id,
                n_rows=data.n_rows,
                d_jet=data.d_jet,
                xn_d=data.xn_d,
                yn_d=data.yn_d,
                z_d=data.z_d,
                span=data.span,
                plate_thickness=data.plate_thickness,
                pattern=data.pattern,
                correlation=data.correlation,
                Cd=_fixed_cd(data),
                cd_crossflow=data.cd_crossflow,
                Nu_multiplier=data.Nu_multiplier * (schema.solver_settings.Nu_multiplier or 1.0),
                roughness=data.roughness,
                bypass_length=data.bypass_length,
                bypass_entry_K=data.bypass_entry_K,
            )
            arr.add_to(net, source_id, target_id, bypass=bypass_id)
            gui_arrays[elem_id] = arr
        elif elem_type == "area_change":
            data = _parse(AreaChangeData, elem_data, elem_id)
            elem = AreaChangeElement(
                id=elem_id,
                from_node=source_id,
                to_node=target_id,
                F0=data.F0,
                F1=data.F1,
                model_type=data.model_type,
                length=data.length,
                D_h=data.D_h,
                # The network's regime (#588): the closure's M -> 0 limit when
                # incompressible (#731).
                regime=schema.solver_settings.global_regime,
            )
            elem.initial_guess = _expand_initial_guess(data.initial_guess, elem_id)
            if not elem.initial_guess:
                elem.initial_guess = _guess_from_prior_result(elem_data, elem_id)
            net.add_element(elem)
        elif elem_type == "lossless_connection":
            elem = LosslessConnectionElement(elem_id, from_node=source_id, to_node=target_id)
            net.add_element(elem)
        elif elem_type in ("compressor", "turbine"):
            data = _parse(TurbomachineData, elem_data, elem_id, "Element")
            cls = CompressorElement if elem_type == "compressor" else TurbineElement
            eta = (
                {"eta_poly": data.efficiency}
                if data.efficiency_type == "polytropic"
                else {"eta_is": data.efficiency}
            )
            elem = cls(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                pressure_ratio=data.pressure_ratio,
                **eta,
            )
            elem.initial_guess = _expand_initial_guess(data.initial_guess, elem_id)
            if not elem.initial_guess:
                elem.initial_guess = _guess_from_prior_result(elem_data, elem_id)
            net.add_element(elem)
        elif elem_type == "discrete_loss":
            data = _parse(DiscreteLossData, elem_data, elem_id)
            # area=None signals "inherit from the topology in resolve_topology".
            # The node-area lookup that used to live here is now part of the
            # library's shared inference path, so GUI and headless networks
            # resolve the same area (issue #262).
            area: float | None = data.area
            correlation = _build_discrete_loss_correlation(
                data.correlation_type,
                data.xi,
                data.k,
                data.xi0,
                data.zeta,
                data.zeta0,
                area if area and area > 0 else 0.1,
                m_ref=data.m_ref if data.m_ref and data.m_ref > 0 else None,
            )
            theta_source = data.theta_source if data.theta_source not in (None, "none") else None
            surface = ConvectiveSurface(
                # No wetted-area field: a wall drawn to it gives one (#1022).
                area=0.0,
                model=map_surface_model(data.surface),
                Nu_multiplier=data.Nu_multiplier,
                f_multiplier=data.f_multiplier,
            )
            elem = PressureLossElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                correlation=correlation,
                theta_source=theta_source,
                area=area,
                surface=surface,
            )
            elem.initial_guess = _expand_initial_guess(data.initial_guess, elem_id)
            if not elem.initial_guess:
                elem.initial_guess = _guess_from_prior_result(elem_data, elem_id)
            net.add_element(elem)
        elif elem_type == "vortex":
            data = _parse(VortexData, elem_data, elem_id)
            shaft_speed = (
                data.shaft_speed_rpm
                if data.shaft_speed_rpm is not None
                else (schema.solver_settings.shaft_speed_rpm or 0.0)
            )
            elem = VortexElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                r_out=data.r_out,
                r_in=data.r_in,
                shaft_speed_rpm=float(shaft_speed),
                v_theta_in=data.v_theta_in,
            )
            elem.initial_guess = _expand_initial_guess(data.initial_guess, elem_id)
            if not elem.initial_guess:
                elem.initial_guess = _guess_from_prior_result(elem_data, elem_id)
            net.add_element(elem)

    # 3. Third Pass: Auto-connections
    for edge in schema.edges:
        if (
            edge.source in nodes_map
            and edge.target in nodes_map
            and not (edge.data or {}).get("type")
        ):
            auto_id = f"__auto_link__{edge.source}__{edge.target}"
            if auto_id not in net.elements:
                elem = LosslessConnectionElement(
                    auto_id, from_node=edge.source, to_node=edge.target
                )
                net.add_element(elem)

    for src, dst in port_node_links:
        link_id = f"__auto_link__{src}__{dst}"
        if link_id not in net.elements:
            # A port chamber discharging into a pressure boundary loses its
            # jet head there and, choked, leaves under-expanded (#637): a
            # Pt-equal link made the reservoir recover the whole head, which
            # has no root once a critical-mode ejector's outlet is sonic.
            if isinstance(nodes_map.get(dst), PressureBoundary):
                link = PortDischargeElement(link_id, from_node=src, to_node=dst)
            else:
                link = LosslessConnectionElement(link_id, from_node=src, to_node=dst)
            net.add_element(link)

    # Ejector warm-starts: port pressures, P_jct and the port-chain flows from
    # the ejector's closed-form physics (only where the user gave no guess).
    # The seeds invert CHOKED flow: compressible ejectors only. The
    # incompressible jet pump has no choke and starts from the solver's own
    # propagation (#588).
    # The seed lives on the element (EjectorElement.seed_warm_start) and runs
    # again at every resolve_topology; seeding here too keeps the guesses
    # visible on the built network.
    for _ej, _ej_ports in _ejector_seeds:
        _ej.seed_warm_start(net)

    # Result assembly reports each array under its GUI id (runner.py).
    net._gui_arrays = gui_arrays

    # Merge chambers (#471): inflows on a momentum chamber's "side-target"
    # handle are side streams; the one inflow on its flow handle is the main
    # inlet, which sets the chamber's axis and area. Declared, never inferred.
    for node_schema in schema.nodes:
        if node_schema.type != "momentum_chamber":
            continue
        nid = node_schema.id
        inflows = [
            e
            for e in schema.edges
            if e.target == nid and not (e.data and e.data.get("type") == "thermal")
        ]
        side = [e for e in inflows if e.targetHandle == "side-target"]
        if not side:
            continue
        main = [e for e in inflows if e.targetHandle != "side-target"]
        if len(main) != 1:
            raise ValueError(
                f"Momentum chamber '{nid}' has side streams and {len(main)} main "
                "inlets; it needs exactly one main inlet on its flow handle."
            )
        src = main[0].source
        if src in splitter_ids:
            # a splitter outlet wired straight to the chamber (#903)
            nodes_map[nid].main_inlet = f"{src}__{nid}"
            continue
        nodes_map[nid].main_inlet = src if src in net.elements else f"__auto_link__{src}__{nid}"

    # 4. Fourth Pass: Thermal Walls
    for edge in schema.edges:
        if edge.data and edge.data.get("type") == "thermal":
            if edge.source in plenum_node_ids or edge.target in plenum_node_ids:
                continue
            wall_id = edge.id
            data = _parse(ThermalWallData, edge.data, wall_id, what="Thermal wall")

            # Build layers list
            layers = []
            if data.layers:
                # New multi-layer format
                layers = [
                    WallLayer(
                        thickness=L.thickness, conductivity=L.conductivity, material=L.material
                    )
                    for L in data.layers
                ]
            else:
                # Legacy scalar format
                t = data.thickness if data.thickness is not None else 0.003
                k = data.conductivity if data.conductivity is not None else 20.0
                layers = [WallLayer(thickness=t, conductivity=k)]

            if edge.source in gui_arrays or edge.target in gui_arrays:
                # One wall to an array is one wall per row, each on that
                # row's target footprint; the edge's own area does not apply.
                array_is_a = edge.source in gui_arrays
                arr = gui_arrays[edge.source if array_is_a else edge.target]
                if edge.source in gui_arrays and edge.target in gui_arrays:
                    raise ValueError(f"Thermal wall '{wall_id}' joins two impingement arrays.")
                arr.add_wall(
                    net,
                    wall_id,
                    hot_element=edge.target if array_is_a else edge.source,
                    layers=layers,
                    R_fouling=data.R_fouling,
                    array_is_side_a=array_is_a,
                )
                continue
            _wetted_area_from_wall(net, edge.source, data.area_a, data.area)
            _wetted_area_from_wall(net, edge.target, data.area_b, data.area)
            net.add_wall(
                ThermalWall(
                    id=wall_id,
                    element_a=edge.source,
                    element_b=edge.target,
                    layers=layers,
                    # Per-side areas (#887), when both are given, replace the one
                    # area; ThermalWall refuses one without the other, by name.
                    contact_area=(
                        data.area if data.area_a is None and data.area_b is None else None
                    ),
                    area_a=data.area_a,
                    area_b=data.area_b,
                    wall_area=data.wall_area,
                    R_fouling=data.R_fouling,
                    orientation=data.orientation,
                    film=data.film or None,
                    film_side=data.film_side,
                    film_start=data.film_start,
                    film_end=data.film_end,
                )
            )
            net.thermal_coupling_enabled = True

    # Every momentum chamber -- drawn, or inserted at a junction port, an edge
    # junction or a merge -- runs the network's regime (#588): an
    # incompressible network keeps no sonic limit hidden in its chambers.
    for node in net.nodes.values():
        if isinstance(node, MomentumChamberNode):
            node.regime = schema.solver_settings.global_regime

    return net
