import axios from "axios";
import {
	addEdge,
	applyEdgeChanges,
	applyNodeChanges,
	type Connection,
	type Edge,
	type EdgeChange,
	MarkerType,
	type Node,
	type NodeChange,
} from "reactflow";
import { create } from "zustand";
import { exportResults, isBackendReachable, resolvePreview } from "../api";
import { discoverFields, REFUSED_ROOT_KEY } from "../utils/diagnostics";
import { repairEdgeHandles } from "../utils/edgeHandles";
import { type ExampleInfo, parseExampleInfo } from "../utils/exampleInfo";
import { HEAT_SOURCE_HANDLE, heatSourceEdgeProps } from "../utils/heatSource";
import { preReleaseHoleCdKeys } from "../utils/holeCd";
import type { PreviewItems } from "../utils/inheritance";
import { migrateNodeData, migrateSolverSettings } from "../utils/legacyKeys";
import {
	migrateNodePositions,
	NETWORK_FILE_VERSION,
} from "../utils/nodeOrigin";
import { describeNonFinite, findNonFinite } from "../utils/nonFinite";
import {
	DISPLAY_ONLY_KINDS,
	isProbeTarget,
	isStructuralNodeChange,
	LINK_FIELD,
	PROBE_HANDLE,
	PROBE_LINK_HANDLE,
	probeOfLink,
} from "../utils/probeLink";
import { solveErrorMessage } from "../utils/solveError";
import { dropSurfacelessWalls } from "../utils/thermalSides";

// Stored per-node/edge results warm-start the next solve (backend
// graph_builder). A structural edit (add/remove node or edge, new
// connection) changes the system those results solved, so they must be
// dropped: a converged 5-tee solution seeding a 6-tee network poisons
// the first solve after the topology edit. Position/selection changes
// keep results, and so does anything done to a probe (#856): the backend
// drops probes before the solve, so a probe never changes the system.
export const stripResults = <T extends Node | Edge>(items: T[]): T[] =>
	items.map((item) => {
		if (!item.data?.result) return item;
		const { result: _result, ...rest } = item.data;
		return { ...item, data: rest };
	});

/** The (invisible) handle on a probe node a connection can be drawn to. */
export { PROBE_HANDLE };

const isStructuralEdgeChange = (changes: EdgeChange[]): boolean =>
	changes.some((c) => c.type === "add" || c.type === "remove");

/** The probe whose derived link (#856) an edge change is about, or null. */
const linkProbe = (c: EdgeChange): string | null =>
	"id" in c ? probeOfLink(c.id) : null;

export interface RFState {
	nodes: Node[];
	edges: Edge[];
	solveResults: any | null;
	onNodesChange: (changes: NodeChange[]) => void;
	onEdgesChange: (changes: EdgeChange[]) => void;
	onConnect: (connection: Connection) => void;
	/** Point a probe at a node or thermal wall, or (null) at nothing. Keeps
	 *  every result: a probe is display-only. An invalid target is ignored. */
	setProbeTarget: (probeId: string, targetId: string | null) => void;
	/** Ids of the selected probe links (utils/probeLink.ts). The links are
	 *  derived from target_id, so their selection lives here, not on edges. */
	selectedProbeLinks: string[];
	setNodes: (nodes: Node[]) => void;
	setEdges: (edges: Edge[]) => void;
	updateNodeData: (
		nodeId: string,
		data: any,
		preserveResults?: boolean,
	) => void;
	updateEdgeData: (
		edgeId: string,
		data: any,
		preserveResults?: boolean,
	) => void;
	setSolveResults: (results: any) => void;
	// Authoritative geometry-inheritance view from /resolve_preview (#559),
	// keyed by canvas node id. The Inspector reads this instead of re-deriving
	// inheritance itself. refreshPreview re-fetches it for the current graph.
	previewItems: PreviewItems;
	refreshPreview: () => Promise<void>;
	speciesMetadata: { names: string[]; molar_masses: number[] } | null;
	fetchSpeciesMetadata: () => Promise<void>;
	displaySettings: string[];
	/** Every field a solve has offered so far: only a field not in it is
	 *  auto-added to displaySettings, so an unchecked field stays unchecked. */
	knownFields: string[];
	setDisplaySettings: (settings: string[]) => void;
	/** Node-card rows the user chose, per palette kind (#775); a kind not in
	 *  it shows its curated default (utils/cardRows.ts). */
	cardRows: Partial<Record<string, string[]>>;
	/** Set a kind's card rows; null restores the curated default. */
	setCardRows: (kind: string, keys: string[] | null) => void;
	solverSettings: {
		global_regime: "incompressible" | "compressible";
		init_strategy:
			| "default"
			| "analytical_pt_prop"
			| "outletref_warmstart"
			| "homotopy"
			| "continuation";
		method: string;
		timeout: number | null;
		/** Hand a stalled cold attempt's budget to the next rung (#864). */
		stall_handover: boolean;
		shaft_speed_rpm: number | null;
		Nu_multiplier: number | null;
		f_multiplier: number | null;
		/** Mach above which an incompressible element is flagged (#588); null = 0.2. */
		incompressible_mach_limit: number | null;
		ambient_T: number | null;
		ambient_P: number | null;
		ambient_RH: number | null;
	};
	updateSolverSettings: (settings: Partial<RFState["solverSettings"]>) => void;
	unitPreferences: {
		pressure: "Pa" | "kPa" | "MPa";
		length: "m" | "mm";
		area: "m²" | "mm²";
		rotSpeed: "rpm" | "Hz";
	};
	setPressureUnit: (unit: "Pa" | "kPa" | "MPa") => void;
	setLengthUnit: (unit: "m" | "mm") => void;
	setAreaUnit: (unit: "m²" | "mm²") => void;
	setRotSpeedUnit: (unit: "rpm" | "Hz") => void;
	saveNetwork: (filename?: string) => void | Promise<void>;
	loadNetwork: (data: any) => void;
	// What a load had to change in the file, shown once above the canvas;
	// also why a load or a save was refused (#1038).
	loadNotice: string | null;
	dismissLoadNotice: () => void;
	// A worked example's description (#793), from the file's `example` key:
	// shown above the canvas until dismissed or another file replaces it.
	exampleInfo: ExampleInfo | null;
	dismissExampleInfo: () => void;
	// When the canvas came back from the autosave (#596), for the note above
	// the canvas. Any load or a Clear replaces that canvas, so both reset it.
	restoredAt: string | null;
	setRestoredAt: (when: string | null) => void;
	// Clear (#596): an empty canvas with no results and no restore note.
	clearCanvas: () => void;
	dismissSolveStatus: () => void;
	isExporting: boolean;
	exportNetworkResults: () => Promise<void>;
}

const useStore = create<RFState>((set, get) => ({
	nodes: [],
	edges: [],
	solveResults: null,
	loadNotice: null,
	exampleInfo: null,
	restoredAt: null,
	isExporting: false,
	previewItems: {},
	selectedProbeLinks: [],

	refreshPreview: async () => {
		const { nodes, edges, solverSettings } = get();
		if (nodes.length === 0) {
			set({ previewItems: {} });
			return;
		}
		try {
			const res = await resolvePreview({
				nodes,
				edges,
				solver_settings: solverSettings,
			});
			set({ previewItems: res.items });
		} catch {
			// Advisory only: on a backend error (or backend down) keep the last
			// good preview; the solve path surfaces real errors.
		}
	},

	onNodesChange: (changes: NodeChange[]) => {
		const before = get().nodes;
		const nodes = applyNodeChanges(changes, before);
		if (isStructuralNodeChange(changes, before)) {
			set({
				nodes: stripResults(nodes),
				edges: stripResults(get().edges),
				solveResults: null,
			});
		} else {
			// A move, a selection, a measured size, or a probe added or removed:
			// the solved system is unchanged, so its results and the solve's
			// status (ledger, second-law and Mach flags) stay. Clicking a node
			// used to drop the status, and with it the Inspector's flags.
			set({ nodes });
		}
	},

	onEdgesChange: (changes: EdgeChange[]) => {
		// A probe's link (#856) is derived from its target_id, not stored:
		// its selection is kept here, and deleting it clears the target.
		const linkChanges = changes.filter((c) => linkProbe(c) !== null);
		if (linkChanges.length > 0) {
			let selectedLinks = get().selectedProbeLinks;
			for (const c of linkChanges) {
				if (c.type === "select") {
					selectedLinks = selectedLinks.filter((id) => id !== c.id);
					if (c.selected) selectedLinks = [...selectedLinks, c.id];
				} else if (c.type === "remove") {
					selectedLinks = selectedLinks.filter((id) => id !== c.id);
					get().setProbeTarget(linkProbe(c) as string, null);
				}
			}
			set({ selectedProbeLinks: selectedLinks });
		}
		const flowChanges = changes.filter((c) => linkProbe(c) === null);
		if (flowChanges.length === 0) return;
		const edges = applyEdgeChanges(flowChanges, get().edges);
		if (isStructuralEdgeChange(flowChanges)) {
			set({
				edges: stripResults(edges),
				nodes: stripResults(get().nodes),
				solveResults: null,
			});
		} else {
			set({ edges });
		}
	},

	setProbeTarget: (probeId: string, targetId: string | null) => {
		const { nodes, edges } = get();
		if (targetId !== null && !isProbeTarget(probeId, targetId, nodes, edges))
			return;
		const owner = nodes.find((n) => n.id === probeId);
		const field = LINK_FIELD[owner?.type ?? ""];
		if (!owner || !field) return;
		const next = nodes.map((n) => {
			if (n.id !== probeId) return n;
			if (targetId === null) {
				const { [field]: _old, ...rest } = n.data ?? {};
				return { ...n, data: rest };
			}
			return { ...n, data: { ...n.data, [field]: targetId } };
		});
		if (DISPLAY_ONLY_KINDS.has(owner.type ?? "")) {
			set({ nodes: next });
			return;
		}
		// A fuel boundary's combustor (#892) is solved: a structural edit.
		set({
			nodes: stripResults(next),
			edges: stripResults(edges),
			solveResults: null,
		});
	},

	onConnect: (connection: Connection) => {
		// A connection drawn to a probe's handle points the probe at the node
		// it starts from. It is not a flow edge: drawn as one it changed
		// nothing the probe showed, and (a structural edit) cleared every
		// result in the network.
		if (connection.targetHandle === PROBE_HANDLE) {
			if (connection.target && connection.source)
				get().setProbeTarget(connection.target, connection.source);
			return;
		}
		// The same from the probe's link handle (#856), should a drag from it
		// end on another node's handle. NetworkCanvas attaches the probe on
		// drop (onConnectEnd), so this is the fallback.
		if (connection.sourceHandle === PROBE_LINK_HANDLE) {
			if (connection.source && connection.target)
				get().setProbeTarget(connection.source, connection.target);
			return;
		}
		// A heat source's edge (#792) names the node it heats: no flow.
		if (connection.sourceHandle === HEAT_SOURCE_HANDLE) {
			const edge: Edge = {
				...connection,
				id: `edge_${Date.now()}`,
				...heatSourceEdgeProps(),
			} as Edge;
			set({
				edges: stripResults(addEdge(edge, get().edges)),
				nodes: stripResults(get().nodes),
				solveResults: null,
			});
			return;
		}
		const isThermal =
			connection.sourceHandle?.includes("thermal") ||
			connection.targetHandle?.includes("thermal");

		const edge: Edge = {
			...connection,
			id: `edge_${Date.now()}`,
			type: isThermal ? "thermal" : "default",
			markerEnd: isThermal
				? {
						type: MarkerType.ArrowClosed,
						color: "rgba(255, 152, 0, 0.45)",
					}
				: undefined,
			data: isThermal
				? {
						type: "thermal",
						thickness: 0.003,
						conductivity: 20.0,
						area: 0.05,
						layers: [
							{ thickness: 0.003, conductivity: 20.0, material: "generic" },
						],
						probe_mode: "preset",
						probe_preset: { type: "hot", index: 0 },
					}
				: {},
		} as Edge;

		// A fuel boundary (#892) connected straight to a combustor controls
		// that one unless it already names another.
		const nodes = get().nodes;
		const src = nodes.find((n) => n.id === connection.source);
		const tgt = nodes.find((n) => n.id === connection.target);
		const autoLink =
			src?.type === "fuel_boundary" &&
			!src.data?.combustor_id &&
			tgt?.type === "combustor";
		set({
			edges: stripResults(addEdge(edge, get().edges)),
			nodes: stripResults(
				autoLink
					? nodes.map((n) =>
							n.id === src.id
								? { ...n, data: { ...n.data, combustor_id: tgt.id } }
								: n,
						)
					: nodes,
			),
			solveResults: null,
		});
	},

	setNodes: (nodes: Node[]) => {
		set({ nodes });
	},

	setEdges: (edges: Edge[]) => {
		set({ edges });
	},

	updateNodeData: (nodeId: string, data: any, preserveResults = false) => {
		if (preserveResults) {
			set({
				nodes: get().nodes.map((node) => {
					if (node.id === nodeId) {
						return { ...node, data: { ...node.data, ...data } };
					}
					return node;
				}),
			});
		} else {
			set({
				nodes: get().nodes.map((node) => {
					// Clear stale result on every node when any edit invalidates the solve
					const { result: _r, ...rest } = node.data ?? {};
					if (node.id === nodeId) {
						return { ...node, data: { ...rest, ...data } };
					}
					return { ...node, data: rest };
				}),
				edges: get().edges.map((edge) => {
					const { result: _r, ...rest } = edge.data ?? {};
					return { ...edge, data: rest };
				}),
				solveResults: null,
			});
		}
	},
	updateEdgeData: (edgeId: string, data: any, preserveResults = false) => {
		if (preserveResults) {
			set({
				edges: get().edges.map((edge) => {
					if (edge.id === edgeId) {
						return { ...edge, data: { ...edge.data, ...data } };
					}
					return edge;
				}),
			});
		} else {
			set({
				nodes: get().nodes.map((node) => {
					const { result: _r, ...rest } = node.data ?? {};
					return { ...node, data: rest };
				}),
				edges: get().edges.map((edge) => {
					const { result: _r, ...rest } = edge.data ?? {};
					if (edge.id === edgeId) {
						return { ...edge, data: { ...rest, ...data } };
					}
					return { ...edge, data: rest };
				}),
				solveResults: null,
			});
		}
	},

	setSolveResults: (results: any) => {
		// null clears the status banner only (start of a solve, dismiss,
		// export): the stored per-node results stay.
		set({ solveResults: results });
		if (!results) return;

		// Auto-promote every NEWLY discovered scalar field into displaySettings
		// so the Live Telemetry shows it without manual opt-in (node cards
		// show their own per-kind rows, utils/cardRows.ts). A field
		// seen before is left as the user set it: one they unchecked stays
		// unchecked at the next solve (it used to be re-added every solve).
		const discovered = discoverFields(results).map((f) => f.key);
		const known = new Set(get().knownFields);
		const fresh = discovered.filter((k) => !known.has(k));
		if (fresh.length > 0) {
			const shown = new Set(get().displaySettings);
			set({
				knownFields: [...known, ...fresh],
				displaySettings: [
					...get().displaySettings,
					...fresh.filter((k) => !shown.has(k)),
				],
			});
		}

		// Inject each result into its node / edge. A node or edge the response
		// has no result for loses its old one: after a failed request (no
		// results at all) or for an item the solve did not cover, a kept result
		// would read as current in the probe, the telemetry and the cards.
		// A refused root (converged, outcome "inconsistent") is told from a
		// solve that did not converge on each result, so the telemetry and the
		// probes keep the difference after the status banner is dismissed.
		const refused = results.success === false && results.converged === true;
		const withResult = <T extends Node | Edge>(item: T, result: any): T => {
			if (result && refused && result.success === false)
				result = { ...result, [REFUSED_ROOT_KEY]: true };
			if (result) return { ...item, data: { ...item.data, result } };
			if (!item.data?.result) return item;
			const { result: _old, ...rest } = item.data;
			return { ...item, data: rest };
		};
		set({
			nodes: get().nodes.map((node) =>
				withResult(
					node,
					results.node_results?.[node.id] ?? results.element_results?.[node.id],
				),
			),
			edges: get().edges.map((edge) => {
				const r = results.edge_results?.[edge.id];
				// Wall results carry no success flag of their own: they share the
				// solve's, so a probe on a wall can flag a non-converged solve.
				return withResult(
					edge,
					r && r.success === undefined
						? { ...r, success: results.success !== false }
						: r,
				);
			}),
		});
	},

	speciesMetadata: null,
	fetchSpeciesMetadata: async () => {
		try {
			const res = await fetch("/metadata/species");
			const data = await res.json();
			set({ speciesMetadata: data });
		} catch (error) {
			console.error("Failed to fetch species metadata:", error);
		}
	},

	displaySettings: ["Tt", "P", "m_dot", "rho"],
	knownFields: ["Tt", "P", "m_dot", "rho"],
	setDisplaySettings: (settings: string[]) => {
		set({ displaySettings: settings });
	},
	cardRows: {},
	setCardRows: (kind: string, keys: string[] | null) => {
		const { [kind]: _old, ...rest } = get().cardRows;
		set({ cardRows: keys === null ? rest : { ...rest, [kind]: keys } });
	},
	solverSettings: {
		// Mirrors combaero.network.components.DEFAULT_REGIME; kept equal by
		// python/tests/test_regime_default_sync.py (#585).
		global_regime: "compressible",
		init_strategy: "default",
		method: "hybr",
		// Mirrors gui/backend/schemas.py DEFAULT_SOLVE_TIMEOUT; kept equal by
		// python/tests/test_solve_timeout_default_sync.py.
		timeout: 30,
		stall_handover: true,
		shaft_speed_rpm: null,
		Nu_multiplier: null,
		f_multiplier: null,
		incompressible_mach_limit: null,
		ambient_T: null,
		ambient_P: null,
		ambient_RH: null,
	},
	updateSolverSettings: (settings: Partial<RFState["solverSettings"]>) => {
		set({
			solverSettings: { ...get().solverSettings, ...settings },
		});
	},

	unitPreferences: {
		pressure: "kPa",
		length: "m",
		area: "m²",
		rotSpeed: "rpm",
	},
	setPressureUnit: (unit: "Pa" | "kPa" | "MPa") => {
		set((state) => ({
			unitPreferences: { ...state.unitPreferences, pressure: unit },
		}));
	},
	setLengthUnit: (unit: "m" | "mm") => {
		set((state) => ({
			unitPreferences: { ...state.unitPreferences, length: unit },
		}));
	},
	setAreaUnit: (unit: "m²" | "mm²") => {
		set((state) => ({
			unitPreferences: { ...state.unitPreferences, area: unit },
		}));
	},
	setRotSpeedUnit: (unit: "rpm" | "Hz") => {
		set((state) => ({
			unitPreferences: { ...state.unitPreferences, rotSpeed: unit },
		}));
	},

	saveNetwork: async (filename?: string) => {
		const { nodes, edges, solverSettings } = get();
		// JSON writes Infinity and NaN as null, which the backend refuses for
		// a number: never write a file that would not load and solve (#1038).
		const nonFinite = findNonFinite({ nodes, edges, solverSettings });
		if (nonFinite.length > 0) {
			set({
				loadNotice: `Not saved: JSON cannot hold a non-finite number, and it would be written as null (${describeNonFinite(nonFinite)}). Set a finite value and save again.`,
			});
			return;
		}
		const data = {
			nodes,
			edges,
			solverSettings,
			version: NETWORK_FILE_VERSION,
		};
		const content = JSON.stringify(data, null, 2);

		// Try File System Access API (Native "Save As" Dialog)
		if ("showSaveFilePicker" in window) {
			try {
				const handle = await (window as any).showSaveFilePicker({
					suggestedName: filename
						? `${filename.trim()}.json`
						: `combaero_network_${new Date().toISOString().split("T")[0]}.json`,
					types: [
						{
							description: "CombAero Network JSON",
							accept: { "application/json": [".json"] },
						},
					],
				});
				const writable = await handle.createWritable();
				await writable.write(content);
				await writable.close();
				return;
			} catch (err) {
				// AbortError is perfectly normal (user clicked Cancel)
				if ((err as Error).name === "AbortError") return;
				console.warn(
					"Native Save Picker failed, falling back to legacy download",
					err,
				);
			}
		} else {
			console.info(
				"File System Access API not supported in this browser. Using legacy download mode.",
			);
		}

		// Fallback: the classic <a download> method
		const blob = new Blob([content], {
			type: "application/json",
		});
		const url = URL.createObjectURL(blob);
		const link = document.createElement("a");
		link.href = url;

		// Smarter filename for sorting in legacy mode (adds timestamp if collision likely)
		const timestamp = new Date()
			.toLocaleTimeString("en-GB", {
				hour: "2-digit",
				minute: "2-digit",
				second: "2-digit",
			})
			.replace(/:/g, "-");

		link.download = filename
			? `${filename.trim()}.json`
			: `combaero_network_${new Date().toISOString().split("T")[0]}_${timestamp}.json`;
		document.body.appendChild(link);
		link.click();
		document.body.removeChild(link);
		URL.revokeObjectURL(url);
	},

	loadNetwork: (data: any) => {
		if (data.nodes && data.edges) {
			// An overflowing literal (1e999) parses as Infinity; carried into
			// the canvas it would be saved as null and refused (#1038). The
			// canvas stays as it was.
			const nonFinite = findNonFinite(data);
			if (nonFinite.length > 0) {
				set({
					loadNotice: `Not loaded: the file has a non-finite number (${describeNonFinite(nonFinite)}). Every numeric field must be finite; edit the file and load it again.`,
				});
				return;
			}
			// Keys renamed since the file was saved (#773: omega_rpm).
			const loadedSettings = migrateSolverSettings(
				data.solverSettings || get().solverSettings,
			);
			// Legacy migration: the inlet-referenced incompressible warm start
			// was removed in 0.8.0; its successor is the outlet-referenced
			// variant (tracks the compressible root 6-7x closer). The backend
			// maps the old value the same way for headless workflows.
			if (loadedSettings.init_strategy === "incompressible_warmstart") {
				loadedSettings.init_strategy = "outletref_warmstart";
			}
			// Files saved before #864 have no stall_handover: the default (on).
			if (typeof loadedSettings.stall_handover !== "boolean") {
				loadedSettings.stall_handover = true;
			}
			// A wall on a node with no convective surface (#720): the solver
			// refuses it and the node has no handle to draw it from.
			const { edges: wallChecked, dropped } = dropSurfacelessWalls(
				data.nodes,
				data.edges,
			);
			// A flow edge naming a handle its node lacks is drawn disconnected
			// though it solves (#901): move it, or say so.
			const {
				edges: keptEdges,
				repaired,
				unresolved,
			} = repairEdgeHandles(data.nodes, wallChecked);
			if (repaired.length > 0 || unresolved.length > 0) {
				console.warn("loadNetwork: edge handles", { repaired, unresolved });
			}
			const preRelease = preReleaseHoleCdKeys(data.nodes);
			const notices = [
				preRelease.length > 0
					? `This file has hole-Cd fields of a pre-release development build, which are no longer read: the solver refuses them. Move each into the node's hole_cd in the file and remove it (${preRelease.join("; ")}).`
					: null,
				dropped.length > 0
					? `Removed ${dropped.length} thermal wall(s) on an orifice, area change, plenum or hole element (${dropped.join(", ")}): it has no surface to exchange heat. Wall a channel, loss element or chamber instead.`
					: null,
				repaired.length > 0
					? `Reconnected ${repaired.length} flow edge end(s) that named a handle their node does not have (${repaired.join("; ")}).`
					: null,
				unresolved.length > 0
					? `${unresolved.length} flow edge end(s) name a port their node does not have and are drawn disconnected; reconnect them by hand (${unresolved.join("; ")}).`
					: null,
			].filter((m): m is string => m !== null);
			set({
				// Up to file version 1.1 positions were top-left corners (#865).
				nodes: migrateNodePositions(data.nodes, data.version).map((n: any) => ({
					...n,
					data: { ...migrateNodeData(n.data, n.type), result: undefined },
				})),
				edges: keptEdges.map((e: any) => ({
					...e,
					data: e.data ? { ...e.data, result: undefined } : e.data,
				})),
				solverSettings: loadedSettings,
				solveResults: null,
				selectedProbeLinks: [],
				restoredAt: null,
				exampleInfo: parseExampleInfo(data.example),
				loadNotice: notices.length > 0 ? notices.join(" ") : null,
			});
		}
	},

	dismissSolveStatus: () => {
		set({ solveResults: null });
	},

	dismissLoadNotice: () => {
		set({ loadNotice: null });
	},

	dismissExampleInfo: () => {
		set({ exampleInfo: null });
	},

	setRestoredAt: (when) => {
		set({ restoredAt: when });
	},

	clearCanvas: () => {
		set({
			nodes: [],
			edges: [],
			solveResults: null,
			selectedProbeLinks: [],
			restoredAt: null,
			loadNotice: null,
			exampleInfo: null,
		});
	},

	exportNetworkResults: async () => {
		const { nodes, edges, solverSettings } = get();
		// The displayed results stay: clearing them wiped the status badge, the
		// diagnostics and the second-law flags, and nothing put them back.
		set({ isExporting: true });

		try {
			const { blob, success, outcome, source } = await exportResults({
				nodes,
				edges,
				solver_settings: solverSettings,
			});
			// The export writes the displayed result (#783). Where none matched
			// the canvas (edited since the solve, or the backend restarted) it
			// solved, and the CSV is that solve: say so. A failed or rejected
			// solve still writes a CSV (every row carries success = False).
			const notes: string[] = [];
			if (source === "solved")
				notes.push(
					"Exported from a fresh solve: no displayed result matched this canvas (it changed since the last solve, or the backend restarted).",
				);
			if (!success)
				notes.push(
					`The exported solve did not succeed (${outcome || "unknown outcome"}): the CSV's rows carry success = False. Check the solve before using these numbers.`,
				);
			if (notes.length > 0) set({ loadNotice: notes.join(" ") });

			// 1. Try Native "Save As" Dialog (Chrome/Edge)
			if ("showSaveFilePicker" in window) {
				try {
					const handle = await (window as any).showSaveFilePicker({
						suggestedName: `combaero_results_${new Date().toISOString().split("T")[0]}.csv`,
						types: [
							{
								description: "CSV File",
								accept: { "text/csv": [".csv"] },
							},
						],
					});
					const writable = await handle.createWritable();
					await writable.write(blob);
					await writable.close();
					return;
				} catch (err) {
					if ((err as Error).name === "AbortError") return;
					console.warn(
						"Native Save Picker failed, falling back to legacy",
						err,
					);
				}
			}

			// 2. Fallback to classic link download
			const url = window.URL.createObjectURL(blob);
			const a = document.createElement("a");
			a.href = url;
			a.download = `combaero_results_${new Date().toISOString().split("T")[0]}.csv`;
			document.body.appendChild(a);
			a.click();
			document.body.removeChild(a);
			window.URL.revokeObjectURL(url);
		} catch (error) {
			console.error("Export failed:", error);
			const backendReachable = axios.isAxiosError(error)
				? await isBackendReachable()
				: true;
			set({
				solveResults: {
					success: false,
					...solveErrorMessage(error, backendReachable, "Export"),
				},
			});
		} finally {
			set({ isExporting: false });
		}
	},
}));

export default useStore;
