# GUI Technical Reference

This document provides the technical specification for the JSON schemas and data flow used by the CombAero Network Designer. This is intended for developers extending the GUI and agents performing automated network generation or inspection.

## Data Flow Architecture

The GUI follows a decoupled architecture:
1.  **Frontend (React Flow)**: Manages the visual graph, nodes, and edges.
2.  **State Persistence**: The graph is serialized to a JSON format (see `NetworkGraphSchema` below).
3.  **Backend (FastAPI)**: Receives the JSON, maps React Flow nodes/edges to C++ physics objects, and invokes the `NetworkSolver`.
4.  **Results Propagation**: The solver returns a `SolveResponse` containing the converged state of every node and element.

## Network Schema

The primary data structure for saving/loading networks is the `NetworkGraphSchema`.

### Root Structure
```json
{
  "nodes": [],           // List of ReactFlowNode objects
  "edges": [],           // List of ReactFlowEdge objects
  "solver_settings": {}  // SolverSettings object
}
```

### Node Data Schemas
Each node in the `nodes` list carries a `data` object corresponding to its type.

#### Pressure Boundary (`pressure_boundary`)
Used for inlets/outlets with fixed total conditions.
- `Pt`: Total Pressure [Pa] (Default: 101325.0)
- `Tt`: Total Temperature [K] (Default: 300.0)
- `composition`: Composition object (Mode: mole/mass, Source: dry_air/fuel/custom)

#### Plenum Node (`plenum`)
Internal mixing node; its pressure is solved (P and Pt are the unknowns, and
a plenum's gas is at rest, so Pt = P). Temperature is not an unknown: it is
propagated from the inflows at every evaluation.
- `initial_guess`: `{"P": float}` (Optional; seeds both P and Pt). The
  Inspector offers P only; a `T` key in a hand-written file has nothing to seed.

#### Combustor Node (`combustor`)
Internal node with heat release logic.
- `method`: "complete" or "equilibrium"
- `area`: Cross-sectional area [m²]

#### Heat Source (`heat_source`)
A heat source or sink (#792): an `EnergyBoundary` on the plenum, momentum
chamber or combustor that its ONE edge (handle `heat-out`, edge
`data.type = "heat"`) points at. Not a flow path. The backend refuses an
unwired source, an edge to an element or boundary, and edges to two nodes.
- `mode`: `"Q"` (default), `"fraction"` or `"Tt_target"`
- `Q`: Heat flow [W], > 0 heats the fluid (Default: 1000.0)
- `fraction`: Fraction of the inflow's sensible enthalpy [-] (Default: 0.01)
- `Tt_target`: The node's outlet TOTAL temperature [K]; the heat is then
  solved for (Default: 600.0)

Its result (`node_results[<id>]`) carries `Q` [W], the heat it applied in
every mode (NaN when no through-flow determines it), `Tt_target` in that mode,
`heated_node`, and the heated node's `state`.

#### Fuel Boundary (`fuel_boundary`)
A fuel controller (#892): a `FuelBoundary` whose fuel flow is SOLVED so that
the combustor named by `combustor_id` is at `Tt_target`. It only supplies
(one outlet handle, `flow-source`). The canvas draws a dashed yellow link to
its combustor, derived from `combustor_id` like a probe's from `target_id`
(never stored as an edge); drag from its link handle (`probe-out`) onto a
combustor to choose one, or connect the fuel straight to a combustor, which
links it when it names none. Changing the combustor is a structural edit.
The backend refuses a fuel boundary with no `combustor_id`.
- `combustor_id`: the combustor it controls (required)
- `Tt_target`: Combustor exit total temperature [K] (Default: 1500.0)
- `branch`: `"lean"` (default) or `"rich"`
- `Tt`: Fuel total temperature [K] (Default: 300.0)
- `composition`: as a boundary's (Default: source `fuel`, CH4)
- `m_dot_guess`: optional start [kg/s]; otherwise sized on the target

Its result `state` carries `m_dot`, `phi`, `Tt_target`, `Tt_without_fuel`
and `rich_branch`.

### Element Data Schemas
Elements are defined by edges in the graph but carry their own physical parameters.

#### Orifice Element (`orifice`)
- `diameter`: Orifice diameter [m]
- `hole_cd`: the discharge coefficient (see Hole Cd below); its
  `correlation` is "ReaderHarrisGallagher" (the default), "IdelchikThick",
  "IdelchikBeveled", "IdelchikRounded", "McGreehanSchotsch", "Lichtarowicz",
  or "fixed" with `Cd`.
  - "Stolz" was removed in 0.9.0 (#876). A saved file with it still opens in the editor, but the backend refuses it: `/solve` returns HTTP 400 naming the orifice and "ReaderHarrisGallagher". It is not migrated silently, since that would change the result.
- `regime`: "default", "incompressible", or "compressible".

#### Hole Cd (`hole_cd`)
Every node built as a hole element -- `orifice`, `effusion_plate`,
`film_cooling`, `slot_film`, `impingement_array`, `effusion_liner` -- carries
its discharge coefficient as one object, `hole_cd`: the library's `HoleCd`
field for field (`HoleCdData` in `gui/backend/schemas.py`), which the
builder passes through unchanged. Fields: `correlation`, `Cd` (with
`"fixed"` only), `inclination`, `inlet_crossflow`, `exit_crossflow`. `null` or
`"auto"` is the node type's own default (`HoleCd`'s docstring lists them);
a choice the element has no such term for is refused when the network is
built. The Inspector edits it as one "Hole Cd" section offering only each
node type's choices (`gui/frontend/src/utils/holeCd.ts`). The exit-crossflow
term's geometry, `exit_crossflow_region` and
`exit_crossflow_orientation_deg`, stays on the film and effusion-plate node.

File version 1.3. A file saved before it (up to v0.11.1) carries
`correlation`, `Cd` and, on an impingement array, `cd_crossflow` loose on the
node: they are lifted into `hole_cd` on load (backend and frontend alike) as
that build read them -- the correlation (its default where the key is
absent), the Cd for `"fixed"` only, `cd_crossflow` as `exit_crossflow` -- and
solve to the same result. A node with both the loose fields and `hole_cd` is
refused.

#### Compressor and Turbine (`compressor`, `turbine`)
A `CompressorElement` / `TurbineElement` at a fixed total-to-total pressure
ratio (#893): Pt_out = PR x Pt_in (compressor) or Pt_in / PR (turbine); the
network sets the flow. One inlet (`flow-target`) and one outlet
(`flow-source`); no surface, so no thermal handles.
- `pressure_ratio`: >= 1 (Default: 4 compressor, 3.8 turbine: a compressor and a turbine dropped together leave 5 % of Pt for the combustor and ducts between them)
- `efficiency_type`: `"polytropic"` (default) or `"isentropic"`
- `efficiency`: in (0, 1] (Default: 0.88 compressor, 0.9 turbine)

Its result carries `power` [W] > 0 (absorbed / delivered), `dh` [J/kg] signed
into the flow, `Tt_in`, `Tt_out`, `Tt_out_isentropic`, `eta_poly`, `eta_is`,
`pressure_ratio`. A reversed-flow root is refused; a machine between two
pressure boundaries is refused (nothing sets its flow).

#### Flow Splitter (`flow_splitter`)
A `FlowSplitter` (#903): one inlet (`port-in-target`) to N outlets
(`port-out<i>-source`, i = 0 .. N-1, top to bottom), each a mass-flow
controller at its fraction of the inlet flow. The card grows with its
outlets; the Inspector edits the fractions and adds or removes the last
outlet (with its connection). No surface, so no thermal handles.
- `fractions`: one per outlet, each in (0, 1); at most one `null`, the
  remainder (1 minus the others). Without a remainder they must sum to 1
  within 1e-9; the solve is refused otherwise, never renormalised.
  (Default: `[0.5, null]`.)
- `reference`: the index of the outlet whose valve is open (Pt_out = Pt_in)
  to START from, `null` for the remainder (else the first guess).
- `auto_reference`: `true` (default): when the start is not the path needing
  the most pressure (a metered outlet would need a pump), the solve moves
  the reference there and solves again; the result's `reference_outlet` is
  the one in use, `reference_start` the start, `reference_swapped` 1 when
  they differ (both counted from 1, as the card numbers the outlets), and
  the solver message says so. `false` keeps the start, and
  such a solve is refused with a message naming the outlet to use.
- `areas`: optional valve reference area [m^2] per outlet, for `valve_Cd`
  and `valve_zeta`.

A port wired to an element meets it at an inserted plenum, as do an
effusion liner's; a plenum or boundary wired straight to a port is used as
it is; a momentum chamber at the inlet feeds an inserted plenum through a
lossless link (it refuses a split). Each outlet carries the inlet's Tt and
composition. The result carries `m_dot` (the throughput), `n_outlets`,
`reference_outlet`, `reference_start`, `reference_swapped`, `n_needs_pump`, `valve_dPt_max`, `fraction_error_max`
and per outlet the lists `rows_m_dot`, `rows_fraction`,
`rows_fraction_achieved`, `rows_valve_dPt`, `rows_valve_CdA`, `rows_valve_Cd`,
`rows_valve_zeta`, `rows_valve_needs_pump` (`null` where an outlet has no
value: the open reference has no metering). `network_to_gui` exports a
script's splitter as one such node.

#### Channel Element (`channel`)
- `L`: Length [m]
- `D`: Hydraulic diameter [m]
- `roughness`: Surface roughness [m]
- `surface`: Surface model (smooth, ribbed, dimpled, pin_fin, impingement).

#### Ejector Element (`ejector`)
Critical-mode supersonic ejector (3 ports: `primary` + `secondary` inlets,
one `outlet`; handle ids `port-{primary,secondary,outlet}-{target,source}`).
- `throat_area`: Primary nozzle throat area A_t [m^2].
- `nozzle_exit_area`: Primary nozzle exit area A_p1 [m^2] (must exceed throat).
- `mixing_area`: Constant-area mixing-section area A_3 [m^2] (must exceed nozzle exit).
- `recovery_efficiency`: Multiplies the lossless mixed stagnation pressure to
  give the critical back pressure P_c* (default 1.0 = no artificial loss).

Diagnostics returned on solve include `entrainment_ratio`,
`P_c_star` (critical back pressure), and `gamma`.

#### Radiation (`radiation`, `radiation_a`, `radiation_b`)
Opt-in radiation (#1071) on an `effusion_plate` or `effusion_liner` hot face
(`radiation`, kind `"none"` or `"gas"`) and on each side of a thermal wall
edge (`radiation_a` on the source, `radiation_b` on the target, also
`"casing"`). Absent or `"none"` builds no radiation. `"gas"` reads
`beam_length` [m], `emissivity_wall`, `luminosity` and the optional overrides
`T_gas` [K], `fuel_air_ratio`, `pressure` [Pa]; `"casing"` reads
`area_ratio`, `emissivity_wall`, `emissivity_casing` and `T_casing` [K]
(`gui/backend/schemas.py` holds the defaults). Results: `q_rad_hot` and
`eps_g` on the plate; `q_rad_a`, `Q_rad_a`, `T_surface_a`, `eps_g_a` (and
`_b`) on the wall edge.

## Solver Settings
```json
{
  "global_regime": "incompressible" | "compressible",
  "init_strategy": "default" | "analytical_pt_prop" | "outletref_warmstart" | "homotopy" | "continuation",
  "method": "hybr" | "lm" | "broyden1" | ...,
  "timeout": 30.0
}
```

## Solve Results Schema
After a solve, the backend returns a `SolveResponse`.

- `success`: boolean
- `message`: Error message if failed.
- `node_results`: Map of node IDs to results.
  - `state`: `{"P": float, "T": float, "Pt": float, "Tt": float, "Y": list[float], ...}`
- `element_results`: Map of element IDs to results.
  - `m_dot`: Mass flow rate [kg/s]

## CSV Export
`POST /export` returns `NetworkResult.to_dataframe()` of the DISPLAYED result
as CSV (ASCII, UTF-8 without BOM, #783). `/solve` keeps its last few solves
keyed by what a solve depends on: the nodes and edges without their
display-only data (`result`, `label`, `rotation`, `label_offset`,
`selected`), probes left out, and the solver settings. An export of the same
network, at any layout and with any labels, writes that solve with the
current labels and does not solve again. Only where none matches (the canvas
changed since the solve, or the backend restarted) does it solve; the header
`X-Combaero-Export-Source` says `displayed` or `solved`, and the GUI shows a
notice for `solved`. The format is one row per entity:

| Columns | Rows | Content |
|---|---|---|
| `type`, `kind`, `is_boundary`, `combustion_method`, `id`, `label`, `parent` | all | what the row is |
| `success`, `converged`, `consistent`, `outcome`, `final_norm` | all | the solve's status |
| `<key> [<unit>]` (`T [K]`, `m_dot [kg/s]`, `Y[N2] [-]`, ...) | result rows | results (#774) |
| `input:<key> [<unit>]` (`input:L [m]`, `input:composition.source`) | entity rows | the case's inputs, nested fields dotted (#784) |
| `initial_guess:<key>` | entity rows | the entity's initial guess |
| `solver:<key>` (`solver:method`, `solver:timeout`, ...) | network row | how the case was solved |

Result rows are `node`, `element` and `thermal_wall`, as before. Since 0.8.0
the input columns follow every result column, and two row types are added:
one `network` row (id `_network`) with the network-wide settings
(`input:global_regime`, ambient conditions, multipliers) and the solver
options, and an `edge` row for an edge with inputs and no result row. A flow
edge's inputs go on its auto-linked element row. Every numeric input carries
its unit (`INPUT_UNIT_MAP` in `gui/backend/units.py`); text inputs are bare.
A text cell starting with `= + - @`, a tab or a carriage return gets a
leading `'`.

**Rows the builder made** (#783) are kept: they are real solved state. The
junction between two drawn elements (`__junction__a__b`), a junction's port
chamber (`__tee_jct__x_port`), an auto link (`__auto_link__a__b`), and an
array's rows and segments (`x__p1`, `x__c2`, ...) each get a `kind`
(`junction`, `port_chamber`, `lossless_connection`, `plate_row`,
`crossflow_segment`, `panel`, ...), a label derived from their parent's
(`junction c1 -> x`, `Jets plate row 1`) and a `parent`: the drawn node or
edge they belong to. A drawn row's `parent` is empty.

**List cells** (`rows_*`, `T_interface`): the DataFrame holds real lists; the
CSV writes each as a strict JSON array, a non-finite value as `null`.
`pd.read_csv(path, converters={col: json.loads})` reads them, and
`load_results` reads them back (and the Python-repr lists of an export made
before #783).

`gui.backend.results.load_results("combaero_results.csv")` reads the file back
into the one-row-per-case results table, inputs included; see
[API_PYTHON.md](API_PYTHON.md#results-table-one-row-per-case).
