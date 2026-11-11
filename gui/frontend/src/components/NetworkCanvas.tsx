import ReactFlow, {
	Background,
	type Connection,
	Controls,
	MiniMap,
	type OnConnectStartParams,
	useReactFlow,
} from "reactflow";
import "reactflow/dist/style.css";
import { useCallback, useMemo, useRef } from "react";
import { useCopyPaste } from "../hooks/useCopyPaste";
import useStore, { stripResults } from "../store/useStore";
import { defaultFilmSide, FILM_HANDLE } from "../utils/filmLink";
import { heatSourceConnectionError } from "../utils/heatSource";
import { NODE_ORIGIN, SNAP_GRID } from "../utils/nodeOrigin";
import {
	DISPLAY_ONLY_KINDS,
	dropTargetId,
	PROBE_LINK_HANDLE,
	probeLinkEdges,
} from "../utils/probeLink";
import { edgeTypes, nodeTypes } from "./flowTypes";
import SolverStatusBadge from "./SolverStatusBadge";

const NetworkCanvas = () => {
	const {
		nodes,
		edges,
		setNodes,
		setEdges,
		onNodesChange,
		onEdgesChange,
		onConnect,
		updateEdgeData,
		setProbeTarget,
		selectedProbeLinks,
	} = useStore();
	useCopyPaste();
	// Each probe's dotted link to its target (#856), derived from target_id:
	// drawn by React Flow but never stored, saved or sent to the backend.
	const shownEdges = useMemo(
		() => [...edges, ...probeLinkEdges(nodes, edges, selectedProbeLinks)],
		[nodes, edges, selectedProbeLinks],
	);
	const { screenToFlowPosition } = useReactFlow();

	const onDragOver = useCallback((event: React.DragEvent) => {
		event.preventDefault();
		event.dataTransfer.dropEffect = "move";
	}, []);

	const onDrop = useCallback(
		(event: React.DragEvent) => {
			event.preventDefault();

			const type = event.dataTransfer.getData("application/reactflow");

			if (typeof type === "undefined" || !type) {
				return;
			}

			// Map client drop position to flow coordinates based on current zoom/pan panning
			const position = screenToFlowPosition({
				x: event.clientX,
				y: event.clientY,
			});

			const id = `node_${Date.now()}`;
			let data = { id };

			if (type === "mass_boundary") {
				data = {
					...data,
					m_dot: 1.0,
					Tt: 300,
					composition: { source: "humid_air", mode: "mole" },
				} as any;
			} else if (type === "fuel_boundary") {
				// #892: CH4 to 1500 K, lean. The combustor is chosen by
				// connecting the fuel to it or dragging its link handle there.
				data = {
					...data,
					Tt_target: 1500,
					branch: "lean",
					Tt: 300,
					composition: { source: "fuel", mode: "mole" },
				} as any;
			} else if (type === "pressure_boundary") {
				data = {
					...data,
					Pt: 101325,
					Tt: 300,
					composition: { source: "humid_air", mode: "mole" },
				} as any;
			} else if (type === "channel") {
				data = {
					...data,
					L: 1.0,
					// Unset: inherits from its neighbours, else the marked 0.1 m
					// default (#570).
					D: null,
					roughness: 1e-5,
					Nu_multiplier: 1.0,
					f_multiplier: 1.0,
				} as any;
			} else if (type === "compressor") {
				// #893: a fixed total-to-total pressure ratio; the network sets
				// the flow. Polytropic efficiency by default.
				data = {
					...data,
					pressure_ratio: 4,
					efficiency_type: "polytropic",
					efficiency: 0.88,
				} as any;
			} else if (type === "turbine") {
				data = {
					...data,
					pressure_ratio: 3.8,
					efficiency_type: "polytropic",
					efficiency: 0.9,
				} as any;
			} else if (type === "flow_splitter") {
				// #903: two outlets, half each; the second takes the rest and
				// is the reference (its valve open).
				data = { ...data, fractions: [0.5, null], reference: null } as any;
			} else if (type === "orifice") {
				data = {
					...data,
					diameter: 0.08,
					hole_cd: { correlation: "ReaderHarrisGallagher" },
				} as any;
			} else if (type === "effusion_liner") {
				// Andrews plate C holes and pitch over a 30 mm backside duct
				// (#471). U1/Vi follows the duct drop against the hole drive.
				data = {
					...data,
					n_segments: 4,
					length: 0.152,
					width: 0.152,
					duct_height: 0.03,
					hole_diameter: 3.27e-3,
					wall_thickness: 6.3e-3,
					pitch_x: 15.24e-3,
					pitch_y: 15.24e-3,
					angle_deg: 90.0,
					entry_K: 0.5,
					wall_conductivity: 20.0,
					gas_film: "none",
					gas_augmentation: 1.0,
					turbulence_intensity: 0.05,
					gas_heat_flux: 0.0,
					radiation: { kind: "none" },
					hole_cd: { correlation: "McGreehanSchotsch" },
				} as any;
			} else if (type === "film_cooling") {
				// Colban Eq. (19) at Gritsch configuration B, the geometry it is
				// scored on (#386): a working example, not invented typical values.
				data = {
					...data,
					hole_type: "shaped",
					hole_diameter: 1.0e-3,
					holes_per_row: 20,
					lateral_pitch: 6.0e-3,
					n_rows: 1,
					wall_thickness: 3.0e-3,
					angle_deg: 30.0,
					area_ratio: 3.5,
					breakout_width: 2.94e-3,
					turbulence_intensity: 0.04,
					hole_cd: { correlation: "IdelchikThick" },
				} as any;
			} else if (type === "slot_film") {
				// Marek & Tacina's D-7958 rig, on which the slot film is scored
				// (#714): a working example, not invented typical values.
				data = {
					...data,
					slot_height: 6.3e-3,
					span: 0.37,
					lip_thickness: 1.0e-3,
					turbulence_intensity: 0.14,
					C_m_multiplier: 1.0,
					hole_cd: { correlation: "IdelchikThick" },
				} as any;
			} else if (type === "effusion_plate") {
				// Andrews et al.'s plate C, the geometry the closure is scored
				// on (#471): a working example, not invented typical values.
				data = {
					...data,
					hole_diameter: 3.27e-3,
					wall_thickness: 6.3e-3,
					pitch_x: 15.24e-3,
					pitch_y: 15.24e-3,
					panel_length: 0.152,
					panel_width: 0.152,
					angle_deg: 90.0,
					hole_cd: { correlation: "IdelchikThick" },
					wall_conductivity: 20.0,
					gas_film: "none",
					gas_augmentation: 1.0,
					gas_baseline: "chamber",
					turbulence_intensity: 0.05,
					gas_heat_flux: 0.0,
					radiation: { kind: "none" },
				} as any;
			} else if (type === "impingement_array") {
				// A Florschuetz test geometry, so the default is a working
				// example inside the correlation's range (#465).
				data = {
					...data,
					n_rows: 10,
					d_jet: 0.00254,
					xn_d: 5.0,
					yn_d: 4.0,
					z_d: 2.0,
					span: 0.122,
					plate_thickness: 0.00254,
					pattern: "inline",
					hole_cd: { correlation: "fixed", Cd: 0.79 },
					Nu_multiplier: 1.0,
				} as any;
			} else if (type === "discrete_loss") {
				data = {
					...data,
					correlation_type: "constant_fraction",
					xi: 0.03,
					zeta: 1.0,
					k: 0.001,
					xi0: 0.02,
					zeta0: 1.0,
					area: null,
					theta_source: null,
				} as any;
			} else if (type === "combustor") {
				data = {
					...data,
					method: "complete",
					// area null: the chamber inherits it from a neighbour (or the 0.1 m^2
					// default), shown greyed; a value here would count as user-set.
					area: null,
					Nu_multiplier: 1.0,
				} as any;
			} else if (type === "momentum_chamber") {
				data = {
					...data,
					area: null,
					Nu_multiplier: 1.0,
					f_multiplier: 1.0,
				} as any;
			} else if (type === "mpce_tee") {
				// F_C / F_branch / psi intentionally omitted so the backend
				// inherits geometry from connected channels (per schema default
				// `None = inherit`).
				data = {
					...data,
					flow_direction: "branch",
					theta_deg: 90,
				} as any;
			} else if (type === "ejector") {
				data = {
					...data,
					throat_area: 3.14e-5,
					nozzle_exit_area: 1.0e-4,
					// null: inherit A_3 from the outlet pipe (or the 8e-4 m^2 default).
					mixing_area: null,
					recovery_efficiency: 1.0,
				} as any;
			} else if (type === "heat_source") {
				data = {
					...data,
					// A working example (#792): 1 kW into the node it is wired to.
					mode: "Q",
					Q: 1000.0,
					fraction: 0.01,
					Tt_target: 600.0,
				} as any;
			} else if (type === "vortex") {
				data = {
					...data,
					r_out: 0.1,
					r_in: 0.05,
					shaft_speed_rpm: null,
					v_theta_in: null,
				} as any;
			}

			const newNode = {
				id,
				type,
				position,
				data,
			};

			// Adding a node is a structural edit: stored results solved a
			// different network and must not warm-start the next solve. A
			// probe is not (#856): it keeps them, and shows its target's value
			// as soon as it is pointed at one.
			setNodes(
				(DISPLAY_ONLY_KINDS.has(type) ? nodes : stripResults(nodes)).concat(
					newNode,
				),
			);
		},
		[nodes, setNodes, screenToFlowPosition],
	);

	// A film (#386): a drag from the holes' film handle that ends on a thermal
	// wall's label box names the holes as that wall's film and selects the
	// wall, so the link and the Inspector's Film section show at once.
	const filmDrag = useRef<string | null>(null);
	// A probe's link (#856): a drag from the probe's link handle that ends on
	// a node or a thermal wall (its path or label box) points the probe at it.
	const probeDrag = useRef<string | null>(null);
	const onConnectStart = useCallback(
		(_: unknown, { nodeId, handleId }: OnConnectStartParams) => {
			filmDrag.current = handleId === FILM_HANDLE ? nodeId : null;
			probeDrag.current = handleId === PROBE_LINK_HANDLE ? nodeId : null;
		},
		[],
	);
	const onConnectEnd = useCallback(
		(event: MouseEvent | TouchEvent) => {
			const holes = filmDrag.current;
			const probe = probeDrag.current;
			filmDrag.current = null;
			probeDrag.current = null;
			const pt = "changedTouches" in event ? event.changedTouches[0] : event;
			if (probe) {
				const targetId = dropTargetId(
					document.elementFromPoint(pt.clientX, pt.clientY),
				);
				// setProbeTarget ignores the probe itself, another probe and
				// a flow edge (no results to read).
				if (targetId) setProbeTarget(probe, targetId);
				return;
			}
			if (!holes) return;
			const wallId = document
				.elementFromPoint(pt.clientX, pt.clientY)
				?.closest("[data-film-wall]")
				?.getAttribute("data-film-wall");
			const wall = edges.find((e) => e.id === wallId);
			if (!wall) return;
			updateEdgeData(wall.id, {
				film: holes,
				film_side: defaultFilmSide(wall, holes, edges),
			});
			setNodes(
				useStore.getState().nodes.map((n) => ({ ...n, selected: false })),
			);
			setEdges(
				useStore
					.getState()
					.edges.map((e) => ({ ...e, selected: e.id === wall.id })),
			);
		},
		[edges, updateEdgeData, setNodes, setEdges, setProbeTarget],
	);
	const isValidConnection = useCallback(
		(c: Connection) =>
			c.sourceHandle !== FILM_HANDLE &&
			c.targetHandle !== FILM_HANDLE &&
			// The probe's link attaches on drop (onConnectEnd), never as an edge.
			c.sourceHandle !== PROBE_LINK_HANDLE &&
			heatSourceConnectionError(c, nodes, edges) === null,
		[nodes, edges],
	);

	return (
		<main
			className="flex-grow h-full bg-stone-50 relative"
			onDragOver={onDragOver}
			onDrop={onDrop}
		>
			<SolverStatusBadge />
			<ReactFlow
				nodes={nodes}
				edges={shownEdges}
				onNodesChange={onNodesChange}
				onEdgesChange={onEdgesChange}
				onConnect={onConnect}
				onConnectStart={onConnectStart}
				onConnectEnd={onConnectEnd}
				isValidConnection={isValidConnection}
				nodeTypes={nodeTypes}
				edgeTypes={edgeTypes}
				nodeOrigin={NODE_ORIGIN}
				snapToGrid
				snapGrid={SNAP_GRID}
				fitView
				fitViewOptions={{ padding: 0.5, maxZoom: 1.0 }}
			>
				<Background />
				<Controls />
				<MiniMap />
			</ReactFlow>
		</main>
	);
};

export default NetworkCanvas;
