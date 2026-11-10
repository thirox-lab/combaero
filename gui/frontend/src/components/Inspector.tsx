import { Activity, Crosshair, RotateCw } from "lucide-react";
import React from "react";
import useStore from "../store/useStore";
import { isCardKind } from "../utils/cardRows";
import {
	probeKeyOptions,
	resolveField,
	unacceptedNote,
	unacceptedResult,
} from "../utils/diagnostics";
import { defaultFilmSide } from "../utils/filmLink";
import { solvedFlowAsMref } from "../utils/fractionLoss";
import {
	formatHeat,
	formatKelvin,
	HEAT_SOURCE_DEFAULTS,
	HEAT_SOURCE_MODES,
	type HeatSourceMode,
	heatedNodeId,
} from "../utils/heatSource";
import { warningsFor } from "../utils/incompressible";
import {
	inheritedField,
	inheritedNote,
	orificePipeNote,
	previewWarning,
	unresolvedReason,
} from "../utils/inheritance";
import { junctionGuessKey } from "../utils/legacyKeys";
import { nodeTypeLabel } from "../utils/nodeUtils";
import { withPortArm } from "../utils/portArms";
import {
	formatSignificant,
	QUANTITY_CATALOGUE,
	quantityMeta,
	THERMAL_WALL_KIND,
} from "../utils/quantities";
import { describeElement, failuresFor } from "../utils/secondLaw";
import {
	resolvedFractions,
	SPLITTER_MAX_OUTLETS,
	splitterFractions,
	splitterReference,
	validateSplitter,
} from "../utils/splitter";
import { bassettDesignK, designNetLoss } from "../utils/teeK";
import { validateNetwork } from "../utils/validation";
import {
	designArrangement,
	orientationGlyph,
	orientationOf,
	sideRun,
	solvedArrangement,
} from "../utils/wallArrangement";
import AreaInput from "./AreaInput";
import { CardRowPin, CardRowsReset } from "./CardRowControls";
import CompositionEditor from "./CompositionEditor";
import GeometryFields from "./GeometryFields";
import LengthInput from "./LengthInput";
import NumericInput from "./NumericInput";
import PlateMaterialSelect from "./PlateMaterialSelect";
import RadiationEditor from "./RadiationEditor";
import SurfaceEnhancementInspector from "./SurfaceEnhancementInspector";
import UnitInput from "./UnitInput";
import WallLayersEditor from "./WallLayersEditor";

// The canvas types built as one OrificeElement (or a plate subclass) under
// their own id, whose bore the solver checks against an adjacent momentum
// chamber (#1016).
const ORIFICE_FAMILY_TYPES = new Set([
	"orifice",
	"effusion_plate",
	"film_cooling",
	"slot_film",
]);

const Inspector = () => {
	const {
		nodes,
		edges,
		updateNodeData,
		updateEdgeData,
		speciesMetadata,
		unitPreferences,
		setRotSpeedUnit,
		isExporting,
		exportNetworkResults,
		displaySettings,
		solverSettings,
		solveResults,
		previewItems,
	} = useStore();

	const validationErrors = validateNetwork(nodes);

	// Get selected node
	const selectedNode = nodes.find((n) => n.selected);
	const selectedEdge = edges.find((e) => e.selected);

	if (!selectedNode && !selectedEdge) {
		return (
			<aside className="w-80 border-l bg-white p-4 flex flex-col gap-4 italic text-gray-400">
				Select a node or edge to edit properties
			</aside>
		);
	}

	// ── Probe Inspector ──────────────────────────────────────────
	if (selectedNode?.type === "probe") {
		const probeData = selectedNode.data;
		const targetElement =
			nodes.find((n) => n.id === probeData.target_id) ||
			edges.find((e) => e.id === probeData.target_id);
		const result = targetElement?.data?.result;
		const targetKind =
			(targetElement as any)?.data?.type === "thermal"
				? THERMAL_WALL_KIND
				: targetElement?.type;

		// The keys this target's result carries (finite numbers only), or,
		// before a solve, the keys a solve gives its kind
		// (resultKeys.generated.json, from the backend).
		const availableKeys = probeKeyOptions(result, targetKind).map((key) => {
			const meta = quantityMeta(key);
			const raw = key.startsWith("state.") ? key.slice(6) : key;
			const label = withPortArm(
				meta.label,
				key,
				targetElement?.type,
				targetElement?.data,
			);
			const unit = meta.unit === "—" ? "" : ` [${meta.unit}]`;
			return {
				key,
				label: meta.catalogued ? `${label}${unit} (${raw})` : `${label}${unit}`,
			};
		});
		availableKeys.sort((a, b) => a.label.localeCompare(b.label));
		// A slot key the target does not carry (picked for another target, or
		// before a solve) stays selectable, marked, rather than vanishing.
		for (const slotKey of [probeData.slot1_key, probeData.slot2_key]) {
			if (slotKey && !availableKeys.some((o) => o.key === slotKey)) {
				availableKeys.push({
					key: slotKey,
					label: `${quantityMeta(slotKey).label} (not on this target)`,
				});
			}
		}

		const slotSelect = (
			slot: "slot1_key" | "slot2_key",
			current: string | undefined,
		) => (
			<select
				className="p-1.5 border rounded text-xs bg-white w-full"
				value={current ?? ""}
				onChange={(e) =>
					selectedNode &&
					updateNodeData(
						selectedNode.id,
						{ [slot]: e.target.value || undefined },
						true,
					)
				}
			>
				<option value="">— select quantity —</option>
				{availableKeys.map((item) => (
					<option key={item.key} value={item.key}>
						{item.label}
					</option>
				))}
			</select>
		);

		return (
			<aside
				key={selectedNode.id}
				className="w-80 border-l bg-white p-4 flex flex-col gap-4 overflow-y-auto"
			>
				<h2 className="text-lg font-bold border-b pb-2 flex items-center gap-2 text-blue-600">
					<Crosshair size={18} />
					<span>Probe</span>
				</h2>

				{/* Label */}
				<div className="flex flex-col gap-1">
					<label className="text-xs font-bold text-gray-500 uppercase">
						Label
					</label>
					<input
						type="text"
						value={probeData.label ?? "Probe"}
						onChange={(e) =>
							updateNodeData(selectedNode.id, { label: e.target.value }, true)
						}
						className="p-1.5 border rounded text-xs w-full focus:ring-1 focus:ring-blue-500 outline-none"
					/>
				</div>

				{/* Target */}
				<div className="flex flex-col gap-1">
					<label className="text-xs font-bold text-gray-500 uppercase">
						Target Element
					</label>
					<select
						className="p-1.5 border rounded text-xs bg-white w-full focus:ring-1 focus:ring-blue-500 outline-none"
						value={probeData.target_id ?? ""}
						// Display-only: re-targeting a probe changes nothing the
						// solve depends on, so it keeps the solution it shows
						// (it used to clear every result in the network).
						onChange={(e) =>
							updateNodeData(
								selectedNode.id,
								{ target_id: e.target.value || undefined },
								true,
							)
						}
					>
						<optgroup label="Network Nodes">
							<option value="">— select target —</option>
							{nodes
								.filter((n) => n.id !== selectedNode?.id && n.type !== "probe")
								.map((n) => (
									<option key={n.id} value={n.id}>
										{n.data?.label || nodeTypeLabel(n.type) || n.id}
									</option>
								))}
						</optgroup>
						<optgroup label="Thermal Walls">
							{edges
								.filter((e) => e.data?.type === "thermal")
								.map((e) => (
									<option key={e.id} value={e.id}>
										{e.data?.label || `Thermal Wall (${e.id})`}
									</option>
								))}
						</optgroup>
					</select>
					{probeData.target_id && !targetElement && (
						<span className="text-amber-500 text-[10px]">
							⚠ Target not found
						</span>
					)}
					<span className="text-stone-400 text-[10px]">
						Or drag from the probe's right-hand dot onto a node or thermal wall.
						Select the dotted link and press Delete to detach it.
					</span>
				</div>

				{/* Slot selectors */}
				<div className="flex flex-col gap-3 border-t pt-3">
					<div className="text-[10px] font-bold text-stone-400 uppercase">
						Display Slots
					</div>
					{!probeData.target_id && (
						<p className="text-[10px] text-stone-400 italic">
							Select a target first to see available quantities.
						</p>
					)}
					{probeData.target_id && (
						<>
							<div className="flex flex-col gap-1">
								<label className="text-[10px] text-stone-500">Slot 1</label>
								{slotSelect("slot1_key", probeData.slot1_key)}
							</div>
							<div className="flex flex-col gap-1">
								<label className="text-[10px] text-stone-500">Slot 2</label>
								{slotSelect("slot2_key", probeData.slot2_key)}
							</div>
						</>
					)}
					{probeData.target_id && targetElement && !result && (
						<p className="text-[10px] text-stone-400 italic">
							Not solved: the list is what a solve gives this kind. Run Solve to
							see values.
						</p>
					)}
					{(() => {
						const why = unacceptedResult(result, solveResults);
						return (
							why && (
								<p className="text-[10px] text-red-500">
									{unacceptedNote(why)}
								</p>
							)
						);
					})()}
				</div>
			</aside>
		);
	}

	if (selectedNode) {
		return (
			// key: remount the whole subtree when the selection changes.
			// React Flow prevents default on node mousedown, so a focused
			// input never blurs when the user clicks another node; without
			// the remount React rewires the still-focused input (holding the
			// typed value) to the NEW node, and the eventual blur commits the
			// old node's edit onto it (observed as "diameter edit applied to
			// 2 nodes").
			<aside
				key={selectedNode.id}
				className="w-80 border-l bg-white p-4 flex flex-col gap-4 overflow-y-auto"
			>
				{validationErrors.length > 0 && (
					<div className="bg-amber-50 border border-amber-200 p-2 rounded text-[10px] text-amber-800 flex flex-col gap-1">
						<div className="font-bold uppercase">Network Warnings</div>
						<ul className="list-disc pl-3">
							{validationErrors.map((err, i) => (
								<li key={i}>{err}</li>
							))}
						</ul>
					</div>
				)}
				<h2 className="text-lg font-bold border-b pb-2 capitalize flex justify-between items-center">
					<span>{nodeTypeLabel(selectedNode.type) || "unknown"} Node</span>
					<div className="flex gap-2">
						<button
							type="button"
							onClick={() => {
								const currentRotation = selectedNode.data.rotation || 0;
								updateNodeData(
									selectedNode.id,
									{ rotation: (currentRotation + 90) % 360 },
									true,
								);
							}}
							className="flex items-center gap-1 text-[10px] bg-blue-50 hover:bg-blue-100 text-blue-600 px-2 py-1 rounded border border-blue-200 transition-colors uppercase font-bold"
							title="Rotate node 90° clockwise"
						>
							<RotateCw size={12} />
							Rotate
						</button>
						<button
							type="button"
							onClick={exportNetworkResults}
							disabled={isExporting}
							className="text-[10px] bg-stone-100 hover:bg-stone-200 px-2 py-1 rounded border border-stone-300 transition-colors uppercase font-bold disabled:opacity-50 min-w-[80px]"
						>
							{isExporting ? "Exporting..." : "Export CSV"}
						</button>
					</div>
				</h2>

				{(() => {
					// #1024: the solver's own geometry verdict from the resolve pass.
					// A geometry it refuses at Solve (overlapping holes, a zero-area
					// duct or node, a Dh no section of its area has), and a closed
					// restriction (zero bore: a wall, no flow).
					const item = previewItems?.[selectedNode.id];
					const refusal = item?.refusal;
					const closed = item?.fields?.closed?.value === 1;
					if (!refusal && !closed) return null;
					return (
						<div className="mb-2 flex flex-col gap-1 text-[10px]">
							{refusal && (
								<p className="rounded border border-amber-200 bg-amber-50 p-2 text-amber-800">
									The solver refuses this geometry: {refusal}
								</p>
							)}
							{closed && (
								<p className="rounded border border-stone-200 bg-stone-50 p-2 text-stone-600">
									Closed: zero area, a wall in the flow path. No flow passes;
									the pressures either side are independent.
								</p>
							)}
						</div>
					);
				})()}

				<div className="flex flex-col gap-2 mb-4">
					<label
						htmlFor={`node_label_${selectedNode.id}`}
						className="text-xs font-bold text-stone-500 uppercase"
					>
						Label / Name
					</label>
					<input
						id={`node_label_${selectedNode.id}`}
						type="text"
						value={selectedNode.data.label || ""}
						onChange={(e) =>
							updateNodeData(selectedNode.id, { label: e.target.value }, true)
						}
						placeholder="e.g. Custom Label"
						className="p-2 border rounded border-stone-200 outline-none focus:ring-2 focus:ring-blue-100"
					/>
				</div>

				<div className="flex flex-col gap-2 opacity-60">
					<label
						htmlFor={`node_id_${selectedNode.id}`}
						className="text-xs font-bold text-gray-500 uppercase"
					>
						Node ID
					</label>
					<input
						id={`node_id_${selectedNode.id}`}
						type="text"
						value={selectedNode.id}
						className="p-1 border rounded bg-gray-50 cursor-not-allowed text-xs"
						disabled
					/>
				</div>
				{/* Second law (#545): this node's control volume -- or a hidden
				    element it owns, such as a tee's port step -- destroys entropy:
				    a passive element creating total pressure. Reported only. */}
				{failuresFor(selectedNode.id, solveResults?.second_law_failures)
					.length > 0 && (
					<div className="rounded border border-amber-300 bg-amber-50 px-2 py-1.5 text-[10px] leading-snug text-amber-800">
						<div className="font-bold uppercase mb-0.5">
							Second law: total pressure created
						</div>
						{failuresFor(
							selectedNode.id,
							solveResults?.second_law_failures,
						).map((f) => (
							<div key={f.element_id} className="font-mono">
								{describeElement(f)}: σ/(R ṁ) = {f.sigma_norm.toExponential(2)}
							</div>
						))}
						<div className="italic text-amber-700 mt-0.5">
							In the solved state total pressure rises across{" "}
							{failuresFor(selectedNode.id, solveResults?.second_law_failures)
								.length > 1
								? "these"
								: "it"}{" "}
							without work or cooling to pay for it, which no passive element
							can do. Usually a modelling gap (an area step or loss left out);
							the solve is not rejected for it.
						</div>
					</div>
				)}
				{/* Rayleigh omission (#551): a wall-heated compressible channel whose
				    momentum leaves out the heating's total-pressure change -- since
				    #672 only where no exit state can be formed. */}
				{selectedNode.data.result?.rayleigh_omitted === 1 && (
					<div className="rounded border border-amber-300 bg-amber-50 px-2 py-1.5 text-[10px] leading-snug text-amber-800">
						<div className="font-bold uppercase mb-0.5">
							Heating not in the momentum
						</div>
						<div className="font-mono">
							omitted ΔPt ≈{" "}
							{(
								(selectedNode.data.result.rayleigh_omitted_frac ?? 0) * 100
							).toFixed(1)}
							% of inlet Pt
						</div>
						<div className="italic text-amber-700 mt-0.5">
							No flow enters this duct's receiving node at this point, so its
							exit state could not be formed and the duct is solved as adiabatic
							friction flow: its wall heat raises the outlet temperature but not
							its momentum. Heating costs total pressure (cooling gains it) as
							Ma²: this is the first-order Rayleigh estimate. The solve is not
							rejected.
						</div>
					</div>
				)}
				{/* Heat-transfer area beyond the friction-wetted area (#551). */}
				{selectedNode.data.result?.heat_area_exceeds_wetted === 1 && (
					<div className="rounded border border-amber-300 bg-amber-50 px-2 py-1.5 text-[10px] leading-snug text-amber-800">
						<div className="font-bold uppercase mb-0.5">
							Heat area beyond the wetted wall
						</div>
						<div className="font-mono">
							heat area ={" "}
							{(selectedNode.data.result.heat_area_ratio ?? 0).toFixed(2)} × the
							duct's wetted area
						</div>
						<div className="italic text-amber-700 mt-0.5">
							Heat flows over surface the duct's friction never sees, so heat
							transfer and pressure loss disagree. In a compressible duct a
							cooling gain can then outgrow the friction and the solve may fail.
							Use the duct's own wetted area, or a ribbed or pin-fin surface
							that carries its own friction.
						</div>
					</div>
				)}
				{/* Incompressible validity (#588): an element on an incompressible
				    closure running above the Mach limit. Reported only. */}
				{warningsFor(selectedNode.id, solveResults?.incompressible_warnings)
					.length > 0 && (
					<div className="rounded border border-amber-300 bg-amber-50 px-2 py-1.5 text-[10px] leading-snug text-amber-800">
						<div className="font-bold uppercase mb-0.5">
							Outside incompressible validity
						</div>
						{warningsFor(
							selectedNode.id,
							solveResults?.incompressible_warnings,
						).map((f) => (
							<div key={f.element_id} className="font-mono">
								{f.element_id === selectedNode.id
									? "this element"
									: f.element_id}
								: Ma = {f.mach.toFixed(3)} &gt; {f.limit}
							</div>
						))}
						<div className="italic text-amber-700 mt-0.5">
							An incompressible closure here drifts from the compressible result
							as Ma² (about 4% on a channel at Ma 0.2, 13% at 0.4). The solve is
							not rejected; set the regime to compressible where it matters.
						</div>
					</div>
				)}
				{selectedNode.type === "plenum" && (
					<InitialGuessEditor node={selectedNode} />
				)}

				{selectedNode.type === "mass_boundary" && (
					<>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`m_dot_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Mass Flow (kg/s)
							</label>
							<NumericInput
								id={`m_dot_${selectedNode.id}`}
								value={selectedNode.data.m_dot ?? 1.0}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										m_dot: val,
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`Tt_mb_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Tt (K)
							</label>
							<NumericInput
								id={`Tt_mb_${selectedNode.id}`}
								value={selectedNode.data.Tt ?? 300}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										Tt: val,
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<CompositionEditor
							key={selectedNode.id}
							nodeId={selectedNode.id}
							data={selectedNode.data}
						/>
						<InitialGuessEditor node={selectedNode} />
					</>
				)}

				{selectedNode.type === "fuel_boundary" &&
					(() => {
						// #892: the fuel flow is solved, set by the target Tt.
						const d = selectedNode.data;
						const combustors = nodes.filter((n) => n.type === "combustor");
						const res = d.result?.state;
						return (
							<div className="flex flex-col gap-3">
								<div className="text-[10px] text-stone-500 leading-snug">
									The fuel flow is solved: it puts the combustor&apos;s exit at
									the target Tt. A target has a lean and a rich root; the one
									asked for is the answer, the other is refused. A target not
									above the combustor&apos;s Tt without this fuel, or above the
									stoichiometric peak, is refused with the reason.
								</div>
								<div className="flex flex-col gap-2">
									<label
										htmlFor={`fb_comb_${selectedNode.id}`}
										className="text-xs font-bold text-gray-500 uppercase"
									>
										Combustor
									</label>
									<select
										id={`fb_comb_${selectedNode.id}`}
										value={d.combustor_id ?? ""}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												combustor_id: e.target.value || undefined,
											})
										}
										className="p-2 border rounded text-sm"
									>
										<option value="">(none)</option>
										{combustors.map((c) => (
											<option key={c.id} value={c.id}>
												{c.data?.label || c.id}
											</option>
										))}
									</select>
								</div>
								<div className="flex flex-col gap-2">
									<label
										htmlFor={`fb_Ttt_${selectedNode.id}`}
										className="text-xs font-bold text-gray-500 uppercase"
									>
										Target Tt (K)
									</label>
									<NumericInput
										id={`fb_Ttt_${selectedNode.id}`}
										value={d.Tt_target ?? 1500}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { Tt_target: val })
										}
										className="p-2 border rounded"
									/>
								</div>
								<div className="flex flex-col gap-2">
									<label
										htmlFor={`fb_branch_${selectedNode.id}`}
										className="text-xs font-bold text-gray-500 uppercase"
									>
										Branch
									</label>
									<select
										id={`fb_branch_${selectedNode.id}`}
										value={d.branch ?? "lean"}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												branch: e.target.value,
											})
										}
										className="p-2 border rounded text-sm"
									>
										<option value="lean">Lean (phi &lt; 1)</option>
										<option value="rich">
											Rich (phi &gt; 1; RQL, reformers)
										</option>
									</select>
								</div>
								<div className="flex flex-col gap-2">
									<label
										htmlFor={`fb_Tt_${selectedNode.id}`}
										className="text-xs font-bold text-gray-500 uppercase"
									>
										Fuel Tt (K)
									</label>
									<NumericInput
										id={`fb_Tt_${selectedNode.id}`}
										value={d.Tt ?? 300}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { Tt: val })
										}
										className="p-2 border rounded"
									/>
								</div>
								<div className="flex flex-col gap-2">
									<label
										htmlFor={`fb_guess_${selectedNode.id}`}
										className="text-xs font-bold text-gray-500 uppercase"
									>
										Fuel flow guess (kg/s, optional)
									</label>
									<NumericInput
										id={`fb_guess_${selectedNode.id}`}
										value={d.m_dot_guess ?? null}
										placeholder="from the target"
										onChange={(val) =>
											updateNodeData(selectedNode.id, { m_dot_guess: val })
										}
										onClear={() =>
											updateNodeData(selectedNode.id, {
												m_dot_guess: undefined,
											})
										}
										className="p-2 border rounded"
									/>
								</div>
								{res && typeof res.m_dot === "number" && (
									<div className="text-[11px] font-mono text-stone-600 leading-snug">
										m_dot {res.m_dot.toPrecision(5)} kg/s
										{typeof res.phi === "number"
											? `, phi ${res.phi.toFixed(4)}`
											: ""}
										{typeof res.Tt_without_fuel === "number"
											? `; ${res.Tt_without_fuel.toFixed(1)} K without it`
											: ""}
									</div>
								)}
								<CompositionEditor
									key={selectedNode.id}
									nodeId={selectedNode.id}
									data={selectedNode.data}
								/>
								<InitialGuessEditor node={selectedNode} />
							</div>
						);
					})()}

				{selectedNode.type === "flow_splitter" &&
					(() => {
						// #903: one inlet to N outlets at given fractions.
						const d = selectedNode.data;
						const fractions = splitterFractions(d);
						const n = fractions.length;
						const shown = resolvedFractions(fractions);
						const ref = splitterReference(d);
						const areas: (number | null)[] =
							d.areas ?? fractions.map(() => null);
						const set = (patch: Record<string, unknown>) =>
							updateNodeData(selectedNode.id, patch);
						const setFraction = (i: number, val: number) =>
							set({ fractions: fractions.map((f, j) => (j === i ? val : f)) });
						const makeRest = (i: number) =>
							// The old remainder keeps the share it had.
							set({
								fractions: fractions.map((f, j) =>
									j === i ? null : f == null ? shown[j] : f,
								),
							});
						const addOutlet = () => {
							if (n >= SPLITTER_MAX_OUTLETS) return;
							const hasRest = fractions.includes(null);
							set({
								fractions: [...fractions, hasRest ? 0.1 : null],
								...(d.areas ? { areas: [...areas, null] } : {}),
							});
						};
						const removeOutlet = () => {
							if (n <= 2) return;
							const last = n - 1;
							const { edges, setEdges } = useStore.getState();
							setEdges(
								edges.filter(
									(e) =>
										!(
											e.source === selectedNode.id &&
											e.sourceHandle === `port-out${last}-source`
										),
								),
							);
							set({
								fractions: fractions.slice(0, last),
								...(d.reference === last ? { reference: null } : {}),
								...(d.areas ? { areas: areas.slice(0, last) } : {}),
							});
						};
						const r = d.result;
						const errors = validateSplitter(d);
						const auto = d.auto_reference ?? true;
						// The outlet whose valve the solve left open (1-based in the
						// result), and the one it started from.
						const inUse =
							typeof r?.reference_outlet === "number"
								? r.reference_outlet - 1
								: null;
						const started =
							typeof r?.reference_start === "number"
								? r.reference_start - 1
								: null;
						const { unitPreferences } = useStore.getState();
						const mm2 = unitPreferences.area === "mm²";
						const aScale = mm2 ? 1e6 : 1;
						return (
							<div className="flex flex-col gap-3">
								<div className="text-[10px] text-stone-500 leading-snug">
									Each outlet takes its fraction of the inlet flow, a mass-flow
									controller: its valve takes whatever drop the path downstream
									leaves (shown after a solve). Fractions sum to 1; "Rest" takes
									what the others leave. One valve stays open: it has to be on
									the path that needs the most pressure, and the solve moves it
									there when it starts elsewhere.
								</div>
								<table className="text-[11px]">
									<thead>
										<tr className="text-gray-500 text-[10px] uppercase">
											<th className="text-left">#</th>
											<th className="text-left">Fraction</th>
											<th>Rest</th>
											<th title="The outlet whose valve is open; the solve starts here">
												Open
											</th>
											<th className="text-left">
												Valve area {mm2 ? "mm²" : "m²"}
											</th>
										</tr>
									</thead>
									<tbody>
										{fractions.map((f, i) => (
											<tr key={`out_${i + 1}`}>
												<td>{i + 1}</td>
												<td>
													{f == null ? (
														<span className="font-mono text-gray-500">
															{shown[i].toPrecision(4)}
														</span>
													) : (
														<NumericInput
															id={`fs_f_${selectedNode.id}_${i}`}
															value={f}
															min={0}
															max={1}
															onChange={(val) => setFraction(i, val)}
															className="p-1 border rounded w-20"
														/>
													)}
												</td>
												<td className="text-center">
													<input
														type="radio"
														name={`fs_rest_${selectedNode.id}`}
														aria-label={`Outlet ${i + 1} takes the rest`}
														checked={f == null}
														onChange={() => makeRest(i)}
													/>
												</td>
												<td className="text-center">
													<input
														type="radio"
														name={`fs_ref_${selectedNode.id}`}
														aria-label={`Outlet ${i + 1} starts with its valve open`}
														checked={ref === i}
														onChange={() => set({ reference: i })}
													/>
												</td>
												<td>
													<NumericInput
														id={`fs_a_${selectedNode.id}_${i}`}
														value={areas[i] == null ? null : areas[i] * aScale}
														min={0}
														placeholder="CdA only"
														onChange={(val) =>
															set({
																areas: areas.map((a, j) =>
																	j === i ? val / aScale : a,
																),
															})
														}
														onClear={() =>
															set({
																areas: areas.map((a, j) =>
																	j === i ? null : a,
																),
															})
														}
														className="p-1 border rounded w-20"
													/>
												</td>
											</tr>
										))}
									</tbody>
								</table>
								<div className="flex gap-2">
									<button
										type="button"
										className="px-2 py-1 border rounded text-xs bg-white hover:bg-stone-50"
										onClick={addOutlet}
										disabled={n >= SPLITTER_MAX_OUTLETS}
									>
										+ Outlet
									</button>
									<button
										type="button"
										className="px-2 py-1 border rounded text-xs bg-white hover:bg-stone-50"
										onClick={removeOutlet}
										disabled={n <= 2}
										title="Removes the last outlet and its connection"
									>
										- Last outlet
									</button>
								</div>
								<label className="flex items-center gap-2 text-[11px] text-stone-600">
									<input
										type="checkbox"
										checked={auto}
										onChange={(e) => set({ auto_reference: e.target.checked })}
									/>
									Move the open valve to the path needing the most pressure
								</label>
								{errors.length > 0 && (
									<div className="text-[10px] text-red-600 leading-snug">
										{errors.join("; ")}
									</div>
								)}
								{inUse != null && started != null && inUse !== started && (
									<div className="text-[10px] text-amber-700 bg-amber-50 border border-amber-200 rounded p-2 leading-snug flex flex-col gap-1">
										<span>
											The solve moved the open valve from outlet {started + 1}{" "}
											to outlet {inUse + 1}: its path needs the most pressure.
										</span>
										<button
											type="button"
											className="self-start px-2 py-0.5 border rounded bg-white hover:bg-stone-50"
											onClick={() => set({ reference: inUse })}
										>
											Start from outlet {inUse + 1}
										</button>
									</div>
								)}
								{r && Array.isArray(r.rows_m_dot) && (
									<div className="flex flex-col gap-1">
										<label className="text-[10px] font-bold text-blue-500 uppercase tracking-wider">
											Per Outlet ({Number(r.m_dot).toPrecision(4)} kg/s in)
										</label>
										<table className="text-[10px] font-mono">
											<thead>
												<tr className="text-gray-500">
													<th className="text-left">#</th>
													<th className="text-right">m g/s</th>
													<th className="text-right">frac</th>
													<th className="text-right">dPt kPa</th>
													<th className="text-right">CdA mm2</th>
													<th className="text-right">zeta</th>
												</tr>
											</thead>
											<tbody>
												{(r.rows_m_dot as number[]).map((m, i) => {
													const cda = r.rows_valve_CdA?.[i];
													const zeta = r.rows_valve_zeta?.[i];
													return (
														<tr key={`res_${i + 1}`}>
															<td>{i + 1}</td>
															<td className="text-right">
																{(1000 * m).toFixed(3)}
															</td>
															<td className="text-right">
																{Number(r.rows_fraction_achieved?.[i]).toFixed(
																	4,
																)}
															</td>
															<td className="text-right">
																{i === inUse
																	? "open"
																	: (
																			Number(r.rows_valve_dPt?.[i]) / 1000
																		).toFixed(3)}
															</td>
															<td className="text-right">
																{typeof cda === "number"
																	? (1e6 * cda).toPrecision(4)
																	: "-"}
															</td>
															<td className="text-right">
																{typeof zeta === "number"
																	? zeta.toPrecision(4)
																	: "-"}
															</td>
														</tr>
													);
												})}
											</tbody>
										</table>
										{r.n_needs_pump > 0 && (
											<div className="text-[10px] text-red-600">
												An outlet needs more pressure than the open one gets:
												open the valve on the most demanding path instead.
											</div>
										)}
									</div>
								)}
							</div>
						);
					})()}

				{(selectedNode.type === "compressor" ||
					selectedNode.type === "turbine") &&
					(() => {
						// #893: a fixed total-to-total pressure ratio.
						const d = selectedNode.data;
						const isTurbine = selectedNode.type === "turbine";
						const res = d.result;
						return (
							<div className="flex flex-col gap-3">
								<div className="text-[10px] text-stone-500 leading-snug">
									{isTurbine
										? "Pt_out = Pt_in / PR. The stream leaves at the expansion's outlet Tt; the network sets the flow."
										: "Pt_out = PR x Pt_in. The stream leaves at the compression's outlet Tt; the network sets the flow."}{" "}
									A reversed flow (surge) is refused.
								</div>
								<div className="flex flex-col gap-2">
									<label
										htmlFor={`tm_pr_${selectedNode.id}`}
										className="text-xs font-bold text-gray-500 uppercase"
									>
										Pressure ratio (-)
									</label>
									<NumericInput
										id={`tm_pr_${selectedNode.id}`}
										value={d.pressure_ratio ?? (isTurbine ? 3.8 : 4)}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { pressure_ratio: val })
										}
										className="p-2 border rounded"
									/>
								</div>
								<div className="flex flex-col gap-2">
									<label
										htmlFor={`tm_et_${selectedNode.id}`}
										className="text-xs font-bold text-gray-500 uppercase"
									>
										Efficiency
									</label>
									<select
										id={`tm_et_${selectedNode.id}`}
										value={d.efficiency_type ?? "polytropic"}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												efficiency_type: e.target.value,
											})
										}
										className="p-2 border rounded text-sm"
									>
										<option value="polytropic">Polytropic</option>
										<option value="isentropic">Isentropic</option>
									</select>
									<NumericInput
										id={`tm_eta_${selectedNode.id}`}
										value={d.efficiency ?? (isTurbine ? 0.9 : 0.88)}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { efficiency: val })
										}
										className="p-2 border rounded"
									/>
								</div>
								{res && typeof res.power === "number" && (
									<div className="text-[11px] font-mono text-stone-600 leading-snug">
										power {(res.power / 1000).toPrecision(5)} kW; Tt{" "}
										{Number(res.Tt_in).toFixed(1)} {"->"}{" "}
										{Number(res.Tt_out).toFixed(1)} K; eta_poly{" "}
										{Number(res.eta_poly).toFixed(4)}, eta_is{" "}
										{Number(res.eta_is).toFixed(4)}
									</div>
								)}
								<InitialGuessEditor node={selectedNode} />
							</div>
						);
					})()}

				{selectedNode.type === "pressure_boundary" && (
					<>
						<UnitInput
							label="Pt"
							value={selectedNode.data.Pt ?? 101325}
							onChange={(val) => updateNodeData(selectedNode.id, { Pt: val })}
						/>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`Tt_pb_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Tt (K)
							</label>
							<NumericInput
								id={`Tt_pb_${selectedNode.id}`}
								value={selectedNode.data.Tt ?? 300}
								onChange={(val) =>
									updateNodeData(selectedNode.id, {
										Tt: val,
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<CompositionEditor
							key={selectedNode.id}
							nodeId={selectedNode.id}
							data={selectedNode.data}
						/>
						<InitialGuessEditor node={selectedNode} />
					</>
				)}

				{selectedNode.type === "wall" && (
					<div className="text-xs text-stone-400 italic px-1">
						Closed end — zero mass flow. No configurable parameters.
					</div>
				)}

				{selectedNode.type === "heat_source" &&
					(() => {
						// #792: an EnergyBoundary on the node its one edge points at.
						const d = selectedNode.data;
						const mode: HeatSourceMode = d.mode ?? HEAT_SOURCE_DEFAULTS.mode;
						const heated = heatedNodeId(selectedNode.id, edges);
						const heatedNode = nodes.find((n) => n.id === heated);
						const res = d.result;
						return (
							<div className="flex flex-col gap-3">
								<div className="text-[10px] text-stone-500 leading-snug">
									Heats{" "}
									<span className="font-mono font-bold">
										{heatedNode
											? heatedNode.data.label || heatedNode.id
											: "nothing yet"}
									</span>
									{heatedNode ? ` (${nodeTypeLabel(heatedNode.type)})` : ""}. Q
									&gt; 0 heats the fluid, Q &lt; 0 cools it.
								</div>
								<div className="flex flex-col gap-2">
									<label
										htmlFor={`hs_mode_${selectedNode.id}`}
										className="text-xs font-bold text-gray-500 uppercase"
									>
										Set
									</label>
									<select
										id={`hs_mode_${selectedNode.id}`}
										value={mode}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												mode: e.target.value as HeatSourceMode,
											})
										}
										className="p-2 border rounded text-sm"
									>
										{HEAT_SOURCE_MODES.map((m) => (
											<option key={m.value} value={m.value}>
												{m.label}
											</option>
										))}
									</select>
								</div>
								{mode === "Q" && (
									<div className="flex flex-col gap-2">
										<label
											htmlFor={`hs_Q_${selectedNode.id}`}
											className="text-xs font-bold text-gray-500 uppercase"
										>
											Q (W)
										</label>
										<NumericInput
											id={`hs_Q_${selectedNode.id}`}
											value={d.Q ?? HEAT_SOURCE_DEFAULTS.Q}
											onChange={(val) =>
												updateNodeData(selectedNode.id, { Q: val })
											}
											className="p-2 border rounded"
										/>
									</div>
								)}
								{mode === "fraction" && (
									<div className="flex flex-col gap-2">
										<label
											htmlFor={`hs_f_${selectedNode.id}`}
											className="text-xs font-bold text-gray-500 uppercase"
										>
											Fraction of sensible enthalpy (-)
										</label>
										<NumericInput
											id={`hs_f_${selectedNode.id}`}
											value={d.fraction ?? HEAT_SOURCE_DEFAULTS.fraction}
											onChange={(val) =>
												updateNodeData(selectedNode.id, { fraction: val })
											}
											className="p-2 border rounded"
										/>
										<span className="text-[9px] text-stone-400">
											Of h(T) - h(298.15 K) of the inflow (of the products on a
											combustor). -0.05 removes 5%.
										</span>
									</div>
								)}
								{mode === "Tt_target" && (
									<div className="flex flex-col gap-2">
										<label
											htmlFor={`hs_T_${selectedNode.id}`}
											className="text-xs font-bold text-gray-500 uppercase"
										>
											Outlet Tt (total) (K)
										</label>
										<NumericInput
											id={`hs_T_${selectedNode.id}`}
											value={d.Tt_target ?? HEAT_SOURCE_DEFAULTS.Tt_target}
											onChange={(val) =>
												updateNodeData(selectedNode.id, { Tt_target: val })
											}
											className="p-2 border rounded"
										/>
										<span className="text-[9px] text-stone-400">
											The node's TOTAL temperature after the heat; the heat flow
											is solved for. On a combustor: the exit Tt, reached by
											heat taken out after the reaction (the fuel flow is
											unchanged).
										</span>
									</div>
								)}
								{res && (
									<div className="grid grid-cols-2 gap-x-2 gap-y-2 p-2 bg-rose-50 border border-rose-200 rounded">
										<div className="flex flex-col">
											<span className="text-stone-400 text-[9px] font-bold">
												{mode === "Tt_target" ? "Q (solved)" : "Q (applied)"}
											</span>
											<span className="font-mono text-xs font-bold">
												{formatHeat(res.Q)}
											</span>
										</div>
										<div className="flex flex-col">
											<span className="text-stone-400 text-[9px] font-bold">
												Tt (total) of {heatedNode?.data.label || heated}
											</span>
											<span className="font-mono text-xs font-bold">
												{formatKelvin(res.state?.Tt)}
											</span>
										</div>
										{!Number.isFinite(res.Q) && (
											<div className="col-span-2 text-[10px] text-red-600 leading-snug">
												No heat determined: without flow through the node no
												heat sets its outlet temperature.
											</div>
										)}
									</div>
								)}
							</div>
						);
					})()}

				{selectedNode.type === "channel" &&
					(() => {
						// Inheritance comes from the solver's own resolve pass
						// (/resolve_preview, #559), not a frontend graph walk, so the
						// value shown is the value the solve will use. The solver never
						// inherits a channel's Dh (unset = circular), so no inherited Dh.
						const geomD = inheritedField(
							previewItems,
							selectedNode.id,
							"diameter",
						);
						return (
							<>
								<LengthInput
									id={`L_${selectedNode.id}`}
									label="Length Flow Path"
									value={selectedNode.data.L ?? 1.0}
									onChange={(val) =>
										updateNodeData(selectedNode.id, { L: val })
									}
								/>
								<GeometryFields
									idPrefix={`ch_${selectedNode.id}`}
									title="Channel cross-section"
									stores="D"
									userValue={selectedNode.data.D}
									inherited={geomD}
									onChange={(v) => updateNodeData(selectedNode.id, { D: v })}
									dh={{
										userValue: selectedNode.data.Dh,
										onChange: (v) => updateNodeData(selectedNode.id, { Dh: v }),
									}}
									missingHint="Nothing to inherit a diameter from — set D or A for this channel."
									defaultHint="Nothing to inherit a diameter from — the 0.1 m default is used. Connect a sized neighbour, or untick Inherit to set D or A."
								/>
								<div className="flex flex-col gap-2">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Friction Model
									</label>
									<select
										className="p-2 border rounded bg-white text-xs"
										value={selectedNode.data.friction_model || "haaland"}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												friction_model: e.target.value,
											})
										}
									>
										<option value="haaland">
											Haaland (default, ~2–3% vs Colebrook)
										</option>
										<option value="serghides">
											Serghides (&lt;0.3% vs Colebrook)
										</option>
										<option value="colebrook">
											Colebrook-White (reference)
										</option>
										<option value="petukhov">
											Petukhov (smooth pipe, pairs with Petukhov Nu)
										</option>
									</select>
								</div>
								<LengthInput
									id={`roughness_${selectedNode.id}`}
									label="Surface Roughness"
									value={selectedNode.data.roughness ?? 1e-5}
									onChange={(val) =>
										updateNodeData(selectedNode.id, { roughness: val })
									}
									disabled={selectedNode.data.friction_model === "petukhov"}
								/>
								<div className="grid grid-cols-2 gap-3">
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											Nu Multiplier
										</label>
										<NumericInput
											value={selectedNode.data.Nu_multiplier ?? 1.0}
											onChange={(val) =>
												updateNodeData(selectedNode.id, { Nu_multiplier: val })
											}
											className="p-2 border rounded"
										/>
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.Nu_multiplier ?? 1}`}
										</span>
									</div>
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											F Multiplier
										</label>
										<NumericInput
											value={selectedNode.data.f_multiplier ?? 1.0}
											onChange={(val) =>
												updateNodeData(selectedNode.id, { f_multiplier: val })
											}
											className="p-2 border rounded"
										/>
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.f_multiplier ?? 1}`}
										</span>
									</div>
								</div>
								<div className="flex flex-col gap-2">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Regime
									</label>
									<select
										className="p-2 border rounded bg-white text-xs"
										value={selectedNode.data.regime || "default"}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												regime: e.target.value,
											})
										}
									>
										<option value="default">Default (Global)</option>
										<option value="incompressible">
											Forced Incompressible
										</option>
										<option value="compressible">Forced Compressible</option>
									</select>
								</div>
								<SurfaceEnhancementInspector
									surface={selectedNode.data.surface || { type: "smooth" }}
									onChange={(surface) =>
										updateNodeData(selectedNode.id, { surface })
									}
								/>
								<InitialGuessEditor node={selectedNode} />
							</>
						);
					})()}

				{selectedNode.type != null &&
					ORIFICE_FAMILY_TYPES.has(selectedNode.type) &&
					(() => {
						// #1016: the bore's own area (a plate's total hole area)
						// against an adjacent momentum chamber, the solver's verdict
						// from the resolve pass: equal is no orifice, larger refused.
						const pt = inheritedField(
							previewItems,
							selectedNode.id,
							"pass_through",
						);
						const chamberPassThrough =
							pt?.value === 1 && pt.source.includes("chamber");
						const reason = unresolvedReason(previewItems, selectedNode.id);
						const chamberRefusal =
							reason?.includes("momentum chamber") === true;
						if (!chamberPassThrough && !chamberRefusal) return null;
						return (
							<p className="text-[9px] text-amber-700">
								{chamberPassThrough
									? `${pt.source.charAt(0).toUpperCase()}${pt.source.slice(1)} — zero loss, the flow passes unchanged.`
									: `The solver refuses this geometry: ${reason}`}
							</p>
						);
					})()}

				{selectedNode.type === "orifice" &&
					(() => {
						// Bore and pipe diameter from the solver's resolve pass
						// (/resolve_preview, #559). The bore never inherits (it is the
						// restriction; unset = the 0.08 m default); the single upstream
						// channel's D -- inherited or not -- sets beta = d / D.
						const geomBore = inheritedField(
							previewItems,
							selectedNode.id,
							"diameter",
						);
						const geomPipe = inheritedField(
							previewItems,
							selectedNode.id,
							"pipe_diameter",
						);
						const boreDiameter: number =
							selectedNode.data.diameter ?? geomBore?.value ?? 0.08;
						const pipeD =
							geomPipe?.value != null && geomPipe.value > 0
								? geomPipe.value
								: null;
						// beta and E come from the solver's resolve pass, E from the
						// same C++ factor the kernels use; beta = 0 means no
						// correction. Reversed flow approaches through the to-node's
						// single pipe (#554).
						const field = (name: string) =>
							inheritedField(previewItems, selectedNode.id, name);
						const positive = (v: number | null | undefined) =>
							v != null && v > 0 ? v : null;
						const beta = positive(field("beta")?.value);
						const E = beta != null ? (field("E")?.value ?? null) : null;
						const geomPipeRev = field("pipe_diameter_reverse");
						const pipeDRev = positive(geomPipeRev?.value);
						const betaRev = positive(field("beta_reverse")?.value);
						const ERev =
							betaRev != null ? (field("E_reverse")?.value ?? null) : null;
						const res = selectedNode.data.result ?? {};
						const q = (key: string) => {
							const v = res[key];
							const spec = QUANTITY_CATALOGUE[key];
							return typeof v === "number" && Number.isFinite(v) && spec
								? `${spec.format(v)} ${spec.unit === "—" ? "" : spec.unit}`.trim()
								: null;
						};
						const metering =
							(selectedNode.data.correlation || "ReaderHarrisGallagher") ===
							"ReaderHarrisGallagher";
						const lu = unitPreferences.length === "mm" ? 1000 : 1;
						// #788: in_pipe = false is a lumped Cd*A restriction (no pipe
						// diameter, beta = 0). In a pipe, a bore equal to it is no
						// orifice (the solver's pass_through verdict); a larger one is
						// refused, and the resolve pass reports it unresolved.
						const inPipe = selectedNode.data.in_pipe !== false;
						// The pipe's verdict only: a chamber's (#1016) has its own note.
						const passThrough =
							field("pass_through")?.value === 1 &&
							!field("pass_through")?.source.includes("chamber");
						// The solver's refusal from the resolve pass, whichever pipe
						// caused it: the downstream one sets beta in reverse (#1039).
						const pipeNote = orificePipeNote({
							preview: previewItems,
							nodeId: selectedNode.id,
							inPipe,
							pipeD,
							passThrough,
							boreDiameter,
						});
						return (
							<>
								<GeometryFields
									idPrefix={`ori_${selectedNode.id}`}
									title="Orifice bore"
									stores="D"
									userValue={selectedNode.data.diameter}
									inherited={geomBore ?? { value: 0.08, source: "default" }}
									inheritable={false}
									onChange={(v) =>
										updateNodeData(selectedNode.id, { diameter: v })
									}
									defaultHint="No bore set — the 0.08 m default bore is used. Set D or A."
								/>
								<div className="rounded bg-stone-50 border border-stone-200 px-3 py-2 text-[10px] text-stone-600 flex flex-col gap-0.5">
									<div className="flex justify-between">
										<span className="text-stone-400">Pipe D (upstream)</span>
										<span className="font-mono">
											{pipeD != null
												? `${(pipeD * lu).toFixed(unitPreferences.length === "mm" ? 1 : 4)} ${unitPreferences.length}`
												: "—"}
										</span>
									</div>
									{pipeD != null && geomPipe?.source && (
										<span className="text-[9px] italic text-stone-400">
											{inheritedNote(geomPipe.source)}.
										</span>
									)}
									{/* Only against a pipe: a chamber pass-through (#1016) also sets
									    beta = 1, with no pipe and the "beta = 0" hint below. */}
									{beta != null && E != null && pipeD != null && (
										<>
											<div className="flex justify-between">
												<span className="text-stone-400">β = d / D</span>
												<span className="font-mono font-semibold">
													{beta.toFixed(4)}
												</span>
											</div>
											<div className="flex justify-between">
												<span className="text-stone-400">
													E (velocity of approach)
												</span>
												<span className="font-mono">{E.toFixed(4)}</span>
											</div>
										</>
									)}
									{pipeD == null && inPipe && (
										<span className="text-[9px] italic text-stone-400">
											No single upstream channel: no velocity-of-approach
											correction (β = 0).
										</span>
									)}
									{!inPipe && (
										<span className="text-[9px] italic text-stone-400">
											No pipe diameter: a lumped Cd·A restriction, a hole in a
											wall (β = 0 both ways, no recovery).
										</span>
									)}
									{betaRev != null && pipeDRev != null && (
										<>
											<div className="flex justify-between mt-1">
												<span className="text-stone-400">
													Pipe D (reverse flow)
												</span>
												<span className="font-mono">
													{`${(pipeDRev * lu).toFixed(unitPreferences.length === "mm" ? 1 : 4)} ${unitPreferences.length}`}
												</span>
											</div>
											{geomPipeRev?.source && (
												<span className="text-[9px] italic text-stone-400">
													{inheritedNote(geomPipeRev.source)}.
												</span>
											)}
											<div className="flex justify-between">
												<span className="text-stone-400">β, E (reverse)</span>
												<span className="font-mono">
													{betaRev.toFixed(4)}, {ERev?.toFixed(4) ?? "—"}
												</span>
											</div>
										</>
									)}
								</div>
								<label className="flex items-center gap-2 text-xs text-gray-600 cursor-pointer select-none">
									<input
										type="checkbox"
										checked={inPipe}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												in_pipe: e.target.checked,
											})
										}
										className="rounded"
									/>
									In a pipe (the adjoining channel sets β)
								</label>
								{pipeNote?.kind === "refused" && (
									<p className="text-[9px] text-amber-700">
										The solver refuses this geometry: {pipeNote.reason}
									</p>
								)}
								{pipeNote?.kind === "equal" && (
									<p className="text-[9px] text-amber-700">
										Bore = pipe D: no orifice — zero loss, the flow passes
										unchanged. Untick "In a pipe" for a lumped Cd·A restriction.
									</p>
								)}
								{pipeNote?.kind === "larger" && (
									<p className="text-[9px] text-amber-700">
										Bore larger than the pipe it sits in: the solver refuses
										this geometry. Untick "In a pipe" for a lumped Cd·A
										restriction.
									</p>
								)}

								<div className="flex flex-col gap-2">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Discharge Model (Cd)
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={
											selectedNode.data.correlation || "ReaderHarrisGallagher"
										}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												correlation: e.target.value,
											})
										}
									>
										<option value="ReaderHarrisGallagher">
											Reader-Harris/Gallagher (Sharp)
										</option>
										<option value="IdelchikThick">
											Deep Hole in Wall (Idelchik 4-18a)
										</option>
										<option value="IdelchikBeveled">
											Beveled Hole in Wall (Idelchik 4-18b)
										</option>
										<option value="IdelchikRounded">
											Rounded Hole in Wall (Idelchik 4-18c)
										</option>
										<option value="McGreehanSchotsch">
											Cooling Hole (M-S 1988, plenum-fed)
										</option>
										<option value="Lichtarowicz">
											Long Orifice, l/d 2-10 (Lichtarowicz 1965)
										</option>
										<option value="fixed">Manual / Fixed Value</option>
									</select>
								</div>

								{/* Conditional Inputs based on correlation */}
								{(selectedNode.data.correlation === "IdelchikThick" ||
									selectedNode.data.correlation === "McGreehanSchotsch" ||
									selectedNode.data.correlation === "Lichtarowicz") && (
									<LengthInput
										id={`plate_thickness_${selectedNode.id}`}
										label="Hole Length (L)"
										value={selectedNode.data.plate_thickness ?? 0.0}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { plate_thickness: val })
										}
									/>
								)}

								{(selectedNode.data.correlation === "IdelchikRounded" ||
									selectedNode.data.correlation === "McGreehanSchotsch") && (
									<LengthInput
										id={`edge_radius_${selectedNode.id}`}
										label="Inlet Edge Radius (r)"
										value={selectedNode.data.edge_radius ?? 0.0}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { edge_radius: val })
										}
									/>
								)}

								{selectedNode.data.correlation === "IdelchikBeveled" && (
									<LengthInput
										id={`bevel_depth_${selectedNode.id}`}
										label="Bevel Depth (l)"
										value={selectedNode.data.bevel_depth ?? 0.0}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { bevel_depth: val })
										}
									/>
								)}

								{metering && (
									<div className="flex flex-col gap-1">
										<label
											htmlFor={`tap_${selectedNode.id}`}
											className="text-xs font-bold text-gray-500 uppercase"
										>
											Pressure Taps (reported)
										</label>
										<select
											id={`tap_${selectedNode.id}`}
											className="p-2 border rounded bg-white text-xs border-stone-200"
											value={selectedNode.data.tap || "Flange"}
											onChange={(e) =>
												updateNodeData(selectedNode.id, { tap: e.target.value })
											}
										>
											<option value="Flange">Flange (25.4 mm)</option>
											<option value="Corner">Corner</option>
											<option value="D_D2">D and D/2</option>
										</select>
										<span className="text-[9px] italic text-stone-400">
											Picks the differential reported as dp at taps. The flow
											does not depend on it: the network evaluates the plate at
											D and D/2 taps.
										</span>
									</div>
								)}

								{/* Solved approach state, tap differential and permanent loss */}
								{(q("dP_tap") || q("P_approach")) && (
									<div className="rounded bg-stone-50 border border-stone-200 px-3 py-2 text-[10px] text-stone-600 flex flex-col gap-0.5">
										{(
											[
												["P_approach", "p1 (approach pipe static)"],
												[
													"dP_tap",
													`dP at ${String(res.tap ?? "Flange").replace("D_D2", "D-D/2")} taps`,
												],
												["dP_permanent", "dP permanent"],
												["loss_ratio", "Permanent / tap dP"],
												["E", "E (along the flow)"],
												["Cd_pipe_ratio", "Hole Cd factor in pipe"],
											] as const
										).map(([key, label]) => {
											const v = q(key);
											// The hole factor only means something off 1 (#645).
											const hidden =
												key === "Cd_pipe_ratio" && res.Cd_pipe_ratio === 1;
											return v == null || hidden ? null : (
												<div key={key} className="flex justify-between">
													<span className="text-stone-400">{label}</span>
													<span className="font-mono">{v}</span>
												</div>
											);
										})}
										{res.pressure_recovery === 0 && q("dP_tap") && (
											<span className="text-[9px] italic text-stone-400">
												No pipe recovery downstream: the jet's head is lost.
											</span>
										)}
									</div>
								)}
								{res.recovery_in_range === 0 && (
									<p className="text-[9px] text-amber-700">
										Tap pressure below 0.75 p1: outside ISO 5167's
										incompressible permanent-loss relation.
									</p>
								)}

								{/* Calculated Result (for correlation models) */}
								{selectedNode.data.correlation !== "fixed" && (
									<div className="flex flex-col gap-1 mb-2 bg-blue-50/50 p-2 rounded border border-blue-100/50">
										<label className="text-[10px] font-bold text-blue-500 uppercase tracking-wider">
											Calculated Cd
										</label>
										<div className="text-sm font-mono font-bold text-blue-700">
											{selectedNode.data.result?.Cd?.toFixed(4) ||
												"— (Pending Solve)"}
										</div>
									</div>
								)}

								{/* Manual Entry (only for Fixed model) */}
								{selectedNode.data.correlation === "fixed" && (
									<div className="flex flex-col gap-1 mt-1">
										<label
											htmlFor={`Cd_${selectedNode.id}`}
											className="text-xs font-bold text-gray-500 uppercase"
										>
											Fixed Cd Value
										</label>
										<NumericInput
											id={`Cd_${selectedNode.id}`}
											value={selectedNode.data.Cd ?? 0.6}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													Cd: val,
												})
											}
											className="p-1.5 h-8 text-sm border rounded bg-white"
											placeholder="0.6"
										/>
									</div>
								)}

								<div className="flex flex-col gap-2 mt-2">
									<label className="text-xs font-bold text-gray-500 uppercase tracking-wide">
										Flow Regime
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={selectedNode.data.regime || "default"}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												regime: e.target.value,
											})
										}
									>
										<option value="default">Default (Global)</option>
										<option value="incompressible">
											Forced Incompressible
										</option>
										<option value="compressible">Forced Compressible</option>
									</select>
								</div>

								{/* Compressible expansion model (#854); the default is
								    selected at build from correlation + geometry (#861) */}
								<div className="flex flex-col gap-1 mt-2">
									<label
										htmlFor={`compressible_model_${selectedNode.id}`}
										className="text-xs font-bold text-gray-500 uppercase tracking-wide"
									>
										Compressible Model
									</label>
									<select
										id={`compressible_model_${selectedNode.id}`}
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={selectedNode.data.compressible_model || "default"}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												compressible_model: e.target.value,
											})
										}
									>
										<option value="default">Default (selected at build)</option>
										<option value="nozzle">Isentropic nozzle</option>
										<option value="jobson1955">
											Sharp-edged jet (Jobson 1955)
										</option>
										<option
											value="iso5167"
											disabled={
												(selectedNode.data.correlation ||
													"ReaderHarrisGallagher") !== "ReaderHarrisGallagher"
											}
										>
											ISO 5167-2 expansibility (RHG only)
										</option>
									</select>
									<span className="text-[9px] italic text-stone-400">
										Compressible regime only. Jobson: the jet widens as the
										pressure ratio falls, through choke; ISO: valid for p2/p1 ≥
										0.75. Default: Jobson for an incompressible-referenced Cd
										(RHG, Idelchik, Lichtarowicz, McGreehan-Schotsch) on a
										sharp-edged bore at β ≤ 0.75; the nozzle for fixed Cd, a
										rounded or beveled edge or a wider bore.
									</span>
									{(selectedNode.data.compressible_model || "default") ===
										"default" &&
										field("jobson1955") != null && (
											<span className="text-[9px] text-stone-500">
												Selected:{" "}
												{field("jobson1955")?.value === 1
													? "Jobson 1955"
													: "isentropic nozzle"}
											</span>
										)}
									{selectedNode.data.compressible_model === "iso5167" &&
										(selectedNode.data.correlation ||
											"ReaderHarrisGallagher") !== "ReaderHarrisGallagher" && (
											<p className="text-[9px] text-amber-700">
												ISO 5167-2 expansibility needs the
												Reader-Harris/Gallagher Cd: the solver refuses this
												pairing.
											</p>
										)}
									{selectedNode.data.result?.model_in_range === 0 && (
										<p className="text-[9px] text-amber-700">
											Outside the compressible model's stated range.
										</p>
									)}
								</div>
								{selectedNode.data.compressible_model === "jobson1955" && (
									<div className="flex flex-col gap-1 mt-1">
										<label
											htmlFor={`Cv_${selectedNode.id}`}
											className="text-xs font-bold text-gray-500 uppercase"
										>
											Velocity Coefficient Cv
										</label>
										<NumericInput
											id={`Cv_${selectedNode.id}`}
											value={selectedNode.data.Cv ?? 1.0}
											onChange={(val) =>
												updateNodeData(selectedNode.id, { Cv: val })
											}
											className="p-1.5 h-8 text-sm border rounded bg-white"
											placeholder="1.0"
										/>
										<span className="text-[9px] italic text-stone-400">
											Cd = Cv Ci: the contraction Ci is Cd / Cv (Jobson used
											0.982).
										</span>
									</div>
								)}
								<InitialGuessEditor node={selectedNode} />
							</>
						);
					})()}

				{selectedNode.type === "effusion_liner" &&
					(() => {
						const d = selectedNode.data;
						const set = (patch: Record<string, unknown>) =>
							updateNodeData(selectedNode.id, patch);
						const num = (key: string, label: string, fallback: number) => (
							<div className="flex flex-col gap-1">
								<label
									htmlFor={`${key}_${selectedNode.id}`}
									className="text-xs font-bold text-gray-500 uppercase"
								>
									{label}
								</label>
								<NumericInput
									id={`${key}_${selectedNode.id}`}
									value={d[key] ?? fallback}
									onChange={(val) => set({ [key]: val })}
									className="p-1.5 h-8 text-sm border rounded bg-white"
									placeholder={String(fallback)}
								/>
							</div>
						);
						const len = (key: string, label: string, fallback: number) => (
							<LengthInput
								id={`${key}_${selectedNode.id}`}
								label={label}
								value={d[key] ?? fallback}
								onChange={(val) => set({ [key]: val })}
							/>
						);
						const r = d.result;
						const rows: number[] | undefined = r?.rows_m_dot_bleed;
						return (
							<>
								<div className="text-[10px] text-gray-500 bg-stone-50 p-2 rounded border border-stone-200">
									Backside coolant duct (C in to C out) whose wall bleeds
									through effusion holes into the gas (G: connect to a chamber's
									"s" handle). Each station is a Bassett bleed and a panel whose
									Cd sees the duct crossflow (McGreehan-Schotsch U1/Vi).
									Coolant-side heat transfer is still the plenum-fed Andrews
									model: provisional (#785).
								</div>
								{num("n_segments", "Stations", 4)}
								{len("length", "Liner Length (along duct)", 0.152)}
								{len("width", "Liner Width", 0.152)}
								{len("duct_height", "Backside Duct Height", 0.03)}
								{len("hole_diameter", "Hole Diameter", 3.27e-3)}
								{len("wall_thickness", "Wall Thickness", 6.3e-3)}
								{len("pitch_x", "Streamwise Pitch", 15.24e-3)}
								{len("pitch_y", "Spanwise Pitch", 15.24e-3)}
								{num("angle_deg", "Hole Angle to Wall [deg]", 90)}
								{num("entry_K", "Duct Entry Loss K", 0.5)}
								{num("wall_conductivity", "Wall Conductivity [W/m/K]", 20)}
								<PlateMaterialSelect
									value={d.material}
									onChange={(material) => set({ material })}
								/>
								<div className="text-[10px] font-bold text-gray-500 uppercase mt-2">
									Gas Side
								</div>
								{num("gas_augmentation", "Gas h Augmentation (user)", 1.0)}
								<div className="flex flex-col gap-1">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Film
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={d.gas_film || "none"}
										onChange={(e) => set({ gas_film: e.target.value })}
									>
										<option value="none">None (validated closure)</option>
										<option value="baldauf_sellers">
											Baldauf + Sellers (not supported by Andrews data)
										</option>
									</select>
								</div>
								{num("gas_heat_flux", "Imposed Heat Flux [W/m²] (plenum)", 0)}
								<div className="flex flex-col gap-1">
									<label
										htmlFor={`film_direction_${selectedNode.id}`}
										className="text-xs font-bold text-gray-500 uppercase"
									>
										Gas Flow Past the Panels
									</label>
									<select
										id={`film_direction_${selectedNode.id}`}
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={d.film_direction || "with_coolant"}
										onChange={(e) => set({ film_direction: e.target.value })}
									>
										<option value="with_coolant">
											With the coolant (station 1 upstream)
										</option>
										<option value="against_coolant">Against the coolant</option>
									</select>
									<p className="text-[10px] text-gray-500 leading-snug">
										Orders the panels for the film this liner lays on a wall
										downstream (drag its "f" handle onto that wall).
									</p>
								</div>
								<RadiationEditor
									id={`radiation_${selectedNode.id}`}
									label="Gas Radiation on the Panels (opt-in)"
									value={d.radiation}
									onChange={(radiation) => set({ radiation })}
								/>
								{rows && rows.length > 0 && (
									<div className="flex flex-col gap-1 mt-2">
										<label className="text-[10px] font-bold text-blue-500 uppercase tracking-wider">
											Per Station (bleed {(100 * r.bleed_fraction).toFixed(0)}%)
										</label>
										<table className="text-[10px] font-mono">
											<thead>
												<tr className="text-gray-500">
													<th className="text-left">#</th>
													<th className="text-right">ṁ g/s</th>
													<th className="text-right">U1/Vi</th>
													<th className="text-right">Cd</th>
													<th className="text-right">P kPa</th>
													<th className="text-right">T wall</th>
												</tr>
											</thead>
											<tbody>
												{rows.map((m, i) => (
													<tr key={`st_${i + 1}`}>
														<td>{i + 1}</td>
														<td className="text-right">
															{(1000 * m).toFixed(2)}
														</td>
														<td className="text-right">
															{r.rows_U1_over_Vi[i].toFixed(3)}
														</td>
														<td className="text-right">
															{r.rows_Cd[i].toFixed(3)}
														</td>
														<td className="text-right">
															{(r.rows_P_backside[i] / 1000).toFixed(2)}
														</td>
														<td className="text-right">
															{r.rows_T_wall_hot
																? r.rows_T_wall_hot[i].toFixed(0)
																: "—"}
														</td>
													</tr>
												))}
											</tbody>
										</table>
										{r.crossflow_cd_degraded > 0 ? (
											<div className="text-[10px] text-amber-600">
												U1/Vi above 0.35: the crossflow Cd degrades badly there
												(Rohde-scored). It follows the duct's pressure drop
												against the hole drive, not the duct size.
											</div>
										) : r.crossflow_cd_beyond_8pct > 0 ? (
											<div className="text-[10px] text-amber-600">
												U1/Vi above 0.2: crossflow Cd beyond its ~8% band.
											</div>
										) : null}
										{r.station_psi_beyond_bassett > 0 && (
											<div className="text-[10px] text-gray-500">
												Duct/hole area ratio {r.station_psi_min.toFixed(1)}{" "}
												beyond Bassett's measured 3 (bleed station extrapolated
												in psi).
											</div>
										)}
										{r.is_ingesting > 0 && (
											<div className="text-[10px] text-red-600">
												At least one station ingests hot gas.
											</div>
										)}
										{r.wall_above_limit === 1 && (
											<div className="text-[10px] text-red-600">
												A panel's wall is above its material's limit (
												{Number(r.T_wall_limit).toFixed(0)} K): it fails.
											</div>
										)}
									</div>
								)}
							</>
						);
					})()}

				{selectedNode.type === "effusion_plate" &&
					(() => {
						const d = selectedNode.data;
						const set = (patch: Record<string, unknown>) =>
							updateNodeData(selectedNode.id, patch);
						const num = (key: string, label: string, fallback: number) => (
							<div className="flex flex-col gap-1">
								<label
									htmlFor={`${key}_${selectedNode.id}`}
									className="text-xs font-bold text-gray-500 uppercase"
								>
									{label}
								</label>
								<NumericInput
									id={`${key}_${selectedNode.id}`}
									value={d[key] ?? fallback}
									onChange={(val) => set({ [key]: val })}
									className="p-1.5 h-8 text-sm border rounded bg-white"
									placeholder={String(fallback)}
								/>
							</div>
						);
						const len = (key: string, label: string, fallback: number) => (
							<LengthInput
								id={`${key}_${selectedNode.id}`}
								label={label}
								value={d[key] ?? fallback}
								onChange={(val) => set({ [key]: val })}
							/>
						);
						const r = d.result;
						return (
							<>
								<div className="text-[10px] text-gray-500 bg-stone-50 p-2 rounded border border-stone-200">
									Plenum-fed plate that owns its wall. Discharge into a momentum
									chamber: the chamber's flow and surface give the gas side.
									Discharge into a plenum: the gas side is the imposed heat flux
									below. Connect it to a chamber's "s" handle to add its jets'
									momentum.
								</div>
								{len("hole_diameter", "Hole Diameter", 3.27e-3)}
								{len("wall_thickness", "Wall Thickness", 6.3e-3)}
								{len("pitch_x", "Streamwise Pitch", 15.24e-3)}
								{len("pitch_y", "Spanwise Pitch", 15.24e-3)}
								{len("panel_length", "Panel Length", 0.152)}
								{len("panel_width", "Panel Width", 0.152)}
								{num("angle_deg", "Hole Angle to Wall [deg]", 90)}
								<div className="flex flex-col gap-1">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Hole Discharge Model (Cd)
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={d.correlation || "IdelchikThick"}
										onChange={(e) => set({ correlation: e.target.value })}
									>
										<option value="IdelchikThick">
											Deep Hole in Wall (Idelchik 4-18a)
										</option>
										<option value="Lichtarowicz">
											Long Orifice, l/d 2-10 (Lichtarowicz 1965)
										</option>
										<option value="McGreehanSchotsch">
											Cooling Hole (M-S 1988, plenum-fed)
										</option>
										<option value="fixed">Manual / Fixed Value</option>
									</select>
								</div>
								{d.correlation === "fixed" && num("Cd", "Fixed Cd Value", 0.6)}
								{(d.correlation || "IdelchikThick") === "IdelchikThick" && (
									<div className="flex flex-col gap-1">
										<label className="text-xs font-bold text-gray-500 uppercase">
											Hole Angle Cd Term
										</label>
										<select
											className="p-2 border rounded bg-white text-xs border-stone-200"
											value={d.cd_inclination || "auto"}
											onChange={(e) => set({ cd_inclination: e.target.value })}
										>
											<option value="auto">
												Default (inclined inlet, Idelchik 3-2)
											</option>
											<option value="none">
												None (angle via hole length only)
											</option>
											<option value="idelchik_3_2">
												Inclined inlet (Idelchik 3-2, forced)
											</option>
										</select>
									</div>
								)}
								{num("wall_conductivity", "Wall Conductivity [W/m/K]", 20)}
								<PlateMaterialSelect
									value={d.material}
									onChange={(material) => set({ material })}
								/>
								<div className="text-[10px] font-bold text-gray-500 uppercase mt-2">
									Gas Side
								</div>
								<div className="flex flex-col gap-1">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Gas h Baseline (h₀)
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={d.gas_baseline || "chamber"}
										onChange={(e) =>
											set(
												e.target.value === "flat_plate"
													? {
															gas_baseline: "flat_plate",
															gas_start_length: d.gas_start_length ?? 0.0,
														}
													: { gas_baseline: "chamber" },
											)
										}
									>
										<option value="chamber">
											Chamber surface (fully developed)
										</option>
										<option value="flat_plate">
											Flat plate, boundary layer restarted upstream
										</option>
									</select>
								</div>
								{d.gas_baseline === "flat_plate" &&
									len(
										"gas_start_length",
										"BL Start Upstream of Panel (bleed, lip, step)",
										0.0,
									)}
								{num("gas_augmentation", "Gas h Augmentation (user)", 1.0)}
								{d.result?.h_gas_unblown_model && (
									<p className="text-[9px] italic text-gray-400 -mt-1">
										Multiplies the unblown h₀ ={" "}
										{Number(d.result.h_gas_unblown).toFixed(1)} W/m²K:{" "}
										{d.result.h_gas_unblown_model} at Re{" "}
										{Number(d.result.h_gas_unblown_Re).toExponential(2)}, Pr{" "}
										{Number(d.result.h_gas_unblown_Pr).toFixed(3)}. An
										augmentation from the literature only carries over against
										the same baseline.
									</p>
								)}
								<div className="flex flex-col gap-1">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Film
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={d.gas_film || "none"}
										onChange={(e) => set({ gas_film: e.target.value })}
									>
										<option value="none">None (validated closure)</option>
										<option value="baldauf_sellers">
											Baldauf + Sellers (not supported by Andrews data)
										</option>
									</select>
								</div>
								{d.gas_film === "baldauf_sellers" &&
									num("turbulence_intensity", "Turbulence Intensity", 0.05)}
								{num("gas_heat_flux", "Imposed Heat Flux [W/m²] (plenum)", 0)}
								<RadiationEditor
									id={`radiation_${selectedNode.id}`}
									label="Gas Radiation on the Hot Face (opt-in)"
									value={d.radiation}
									onChange={(radiation) => set({ radiation })}
								/>
								{r && r.T_wall_hot != null && (
									<div className="flex flex-col gap-0.5 mt-2 bg-blue-50/50 p-2 rounded border border-blue-100/50 text-[10px] font-mono">
										<span>
											Gas side:{" "}
											{r.gas_side_chamber ? "chamber (flow)" : "imposed flux"}
										</span>
										<span>
											T wall hot/cold: {r.T_wall_hot.toFixed(0)} /{" "}
											{r.T_wall_cold.toFixed(0)} K
										</span>
										<span>q: {(r.q_wall / 1000).toFixed(1)} kW/m²</span>
										{r.gas_side_chamber ? (
											<>
												<span>
													T gas: {r.T_gas.toFixed(0)} K (
													{r.T_gas_from_main_inlet ? "main inlet" : "chamber"})
												</span>
												<span>η overall: {r.eta_overall.toFixed(3)}</span>
												<span>h gas: {r.h_gas.toFixed(0)} W/m²/K</span>
												<span>M: {r.blowing_ratio.toFixed(2)}</span>
												{r.q_rad_hot != null &&
													Number.isFinite(r.q_rad_hot) && (
														<span>
															q radiation: {(r.q_rad_hot / 1000).toFixed(2)}{" "}
															kW/m² (ε gas {Number(r.eps_g).toFixed(3)})
														</span>
													)}
											</>
										) : null}
										{r.film_extrapolated > 0 && (
											<span className="text-amber-600">
												Film outside Baldauf's range
											</span>
										)}
										{r.coolant_crossflow_ignored > 0 && (
											<span className="text-amber-600">
												Supply has crossflow: not modelled (plenum-fed plate)
											</span>
										)}
										{r.wall_above_limit === 1 && (
											<span className="text-red-600">
												Wall above its material's limit (
												{Number(r.T_wall_limit).toFixed(0)} K): it fails
												{r.Re_hole === 0 ? " (no coolant: adiabatic)" : ""}.
											</span>
										)}
									</div>
								)}
								{r && r.no_cooling_path === 1 && (
									<div className="text-[10px] text-red-600 mt-2">
										No cooling path: wall T undetermined (heats until failure).
										No coolant flows and the gas side has no flow.
									</div>
								)}
							</>
						);
					})()}

				{selectedNode.type === "film_cooling" &&
					(() => {
						const d = selectedNode.data;
						const set = (patch: Record<string, unknown>) =>
							updateNodeData(selectedNode.id, patch);
						const num = (key: string, label: string, fallback: number) => (
							<div className="flex flex-col gap-1">
								<label
									htmlFor={`${key}_${selectedNode.id}`}
									className="text-xs font-bold text-gray-500 uppercase"
								>
									{label}
								</label>
								<NumericInput
									id={`${key}_${selectedNode.id}`}
									value={d[key] ?? fallback}
									onChange={(val) => set({ [key]: val })}
									className="p-1.5 h-8 text-sm border rounded bg-white"
									placeholder={String(fallback)}
								/>
							</div>
						);
						const len = (key: string, label: string, fallback: number) => (
							<LengthInput
								id={`${key}_${selectedNode.id}`}
								label={label}
								value={d[key] ?? fallback}
								onChange={(val) => set({ [key]: val })}
							/>
						);
						const shaped = (d.hole_type ?? "shaped") === "shaped";
						return (
							<>
								<div className="text-[10px] text-gray-500 bg-stone-50 p-2 rounded border border-stone-200">
									Film-cooling holes: coolant from the inlet into the gas. They
									own no wall. Select this node as the Film of a thermal wall on
									the gas duct downstream; that side's driving temperature then
									follows the film.
								</div>
								<div className="flex flex-col gap-1">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Hole Type
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={d.hole_type || "shaped"}
										onChange={(e) => set({ hole_type: e.target.value })}
									>
										<option value="shaped">Shaped (Colban 2011)</option>
										<option value="cylindrical">
											Cylindrical (Baldauf 2002)
										</option>
									</select>
								</div>
								{len(
									"hole_diameter",
									shaped ? "Metering Diameter" : "Hole Diameter",
									1.0e-3,
								)}
								{num("holes_per_row", "Holes per Row", 20)}
								{len("lateral_pitch", "Lateral Pitch", 6.0e-3)}
								{num("n_rows", "Rows", 1)}
								{(d.n_rows ?? 1) > 1 && len("row_pitch", "Row Pitch", 0.01)}
								{len("wall_thickness", "Wall Thickness", 3.0e-3)}
								{num("angle_deg", "Hole Angle to Surface [deg]", 30)}
								{shaped ? (
									<>
										{num("area_ratio", "Area Ratio A_exit/A_inlet", 3.5)}
										{len("breakout_width", "Breakout Width t", 2.94e-3)}
									</>
								) : (
									num("turbulence_intensity", "Turbulence Intensity", 0.04)
								)}
								<div className="flex flex-col gap-1">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Hole Discharge Model (Cd)
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={d.correlation || "IdelchikThick"}
										onChange={(e) => set({ correlation: e.target.value })}
									>
										<option value="IdelchikThick">
											Deep Hole in Wall (Idelchik 4-18a)
										</option>
										<option value="Lichtarowicz">
											Long Orifice, l/d 2-10 (Lichtarowicz 1965)
										</option>
										<option value="McGreehanSchotsch">
											Cooling Hole (M-S 1988, plenum-fed)
										</option>
										<option value="fixed">Manual / Fixed Value</option>
									</select>
								</div>
								{d.correlation === "fixed" && num("Cd", "Fixed Cd Value", 0.6)}
								{(d.correlation || "IdelchikThick") === "IdelchikThick" && (
									<div className="flex flex-col gap-1">
										<label className="text-xs font-bold text-gray-500 uppercase">
											Hole Angle Cd Term
										</label>
										<select
											className="p-2 border rounded bg-white text-xs border-stone-200"
											value={d.cd_inclination || "auto"}
											onChange={(e) => set({ cd_inclination: e.target.value })}
										>
											<option value="auto">
												Default (Idelchik 3-2 if cylindrical)
											</option>
											<option value="none">
												None (angle via hole length only)
											</option>
											<option value="idelchik_3_2">
												Inclined inlet (Idelchik 3-2, forced)
											</option>
										</select>
									</div>
								)}
							</>
						);
					})()}

				{selectedNode.type === "slot_film" &&
					(() => {
						const d = selectedNode.data;
						const set = (patch: Record<string, unknown>) =>
							updateNodeData(selectedNode.id, patch);
						const num = (key: string, label: string, fallback: number) => (
							<div className="flex flex-col gap-1">
								<label
									htmlFor={`${key}_${selectedNode.id}`}
									className="text-xs font-bold text-gray-500 uppercase"
								>
									{label}
								</label>
								<NumericInput
									id={`${key}_${selectedNode.id}`}
									value={d[key] ?? fallback}
									onChange={(val) => set({ [key]: val })}
									className="p-1.5 h-8 text-sm border rounded bg-white"
									placeholder={String(fallback)}
								/>
							</div>
						);
						const len = (key: string, label: string, fallback: number) => (
							<LengthInput
								id={`${key}_${selectedNode.id}`}
								label={label}
								value={d[key] ?? fallback}
								onChange={(val) => set({ [key]: val })}
							/>
						);
						return (
							<>
								<div className="text-[10px] text-gray-500 bg-stone-50 p-2 rounded border border-stone-200">
									A 2-D film slot (Marek &amp; Tacina 1975): coolant from the
									inlet into the gas. It owns no wall. Select this node as the
									Film of a thermal wall on the gas duct downstream; the span is
									measured from the slot lip. Defaults: the D-7958 rig.
								</div>
								{len("slot_height", "Slot Height S", 6.3e-3)}
								{len("span", "Slot Span", 0.37)}
								{len("lip_thickness", "Lip Thickness", 1.0e-3)}
								{num(
									"turbulence_intensity",
									"Free-Stream Turbulence Tu [-] (0.14 = 14%)",
									0.14,
								)}
								{num("C_m_multiplier", "C_m Multiplier (rig tuning)", 1.0)}
								<div className="flex flex-col gap-1">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Slot Discharge Model (Cd)
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={d.correlation || "IdelchikThick"}
										onChange={(e) => set({ correlation: e.target.value })}
									>
										<option value="IdelchikThick">
											Deep Hole in Wall (Idelchik 4-18a, on Dh)
										</option>
										<option value="Lichtarowicz">
											Long Orifice (Lichtarowicz 1965, on Dh)
										</option>
										<option value="McGreehanSchotsch">
											Cooling Hole (M-S 1988, on Dh)
										</option>
										<option value="fixed">Manual / Fixed Value</option>
									</select>
								</div>
								{d.correlation === "fixed" && num("Cd", "Fixed Cd Value", 0.6)}
								<div className="flex flex-col gap-1">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Film-Side h on the Filmed Wall
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={d.film_htc || "none"}
										onChange={(e) => set({ film_htc: e.target.value })}
									>
										<option value="none">None (the duct's own h)</option>
										<option value="seban1960">
											Seban 1960 (slot flow; recommended)
										</option>
										<option value="ballal_lefebvre1973">
											Ballal-Lefebvre wall jet (book Eq. 8.44, m &gt; 1.3)
										</option>
									</select>
								</div>
							</>
						);
					})()}

				{selectedNode.type === "impingement_array" &&
					(() => {
						const d = selectedNode.data;
						const set = (patch: Record<string, unknown>) =>
							updateNodeData(selectedNode.id, patch);
						const dJet: number = d.d_jet ?? 0.00254;
						const ynD: number = d.yn_d ?? 4.0;
						const span: number = d.span ?? 0.122;
						const holesPerRow = Math.max(1, Math.round(span / (ynD * dJet)));
						const ratio = (key: string, label: string, fallback: number) => (
							<div className="flex flex-col gap-1">
								<label
									htmlFor={`${key}_${selectedNode.id}`}
									className="text-xs font-bold text-gray-500 uppercase"
								>
									{label}
								</label>
								<NumericInput
									id={`${key}_${selectedNode.id}`}
									value={d[key] ?? fallback}
									onChange={(val) => set({ [key]: val })}
									className="p-1.5 h-8 text-sm border rounded bg-white"
									placeholder={String(fallback)}
								/>
							</div>
						);
						const rows: number[] | undefined = d.result?.rows_Nu;
						const hasBypass = edges.some(
							(e) =>
								e.target === selectedNode.id &&
								e.targetHandle === "port-bypass-target",
						);
						return (
							<>
								<div className="text-[10px] text-gray-500 bg-stone-50 p-2 rounded border border-stone-200">
									Jet plate fed from the upstream node, impinging on a target;
									the spent air leaves down the gap to the downstream node
									(Florschuetz, Truman and Metzger 1981). Gc/Gj comes from the
									solved crossflow. Spent air leaving through the target (#468)
									is a different element.
								</div>
								<div
									className={`text-[10px] p-2 rounded border ${
										hasBypass
											? "text-purple-700 bg-purple-50 border-purple-200"
											: "text-gray-500 bg-stone-50 border-stone-200"
									}`}
								>
									{hasBypass ? (
										<>
											<b>Bypass wired (b):</b> an initial crossflow enters the
											gap ahead of row 1 (#467). Flow split checked against
											Florschuetz, Metzger and Su (1982, NASA CR-3630). From
											another array's exit the gap continues with its momentum.
											h is referenced to the jet supply: valid for a bypass at
											the jet temperature only, flagged otherwise.
										</>
									) : (
										<>
											<b>Upstream end closed.</b> Wire the "b" handle for an
											initial (bypass) crossflow (#467).
										</>
									)}
								</div>
								{ratio("n_rows", "Rows", 10)}
								<LengthInput
									id={`d_jet_${selectedNode.id}`}
									label="Hole Diameter (d)"
									value={dJet}
									onChange={(val) => set({ d_jet: val })}
								/>
								{ratio("xn_d", "Streamwise Pitch xn/d", 5.0)}
								{ratio("yn_d", "Spanwise Pitch yn/d", 4.0)}
								{ratio("z_d", "Gap z/d", 2.0)}
								<LengthInput
									id={`span_${selectedNode.id}`}
									label="Span"
									value={span}
									onChange={(val) => set({ span: val })}
								/>
								<div className="text-[9px] text-gray-400 -mt-1">
									{holesPerRow} holes per row
								</div>
								<LengthInput
									id={`plate_thickness_${selectedNode.id}`}
									label="Jet Plate Thickness"
									value={d.plate_thickness ?? 0.00254}
									onChange={(val) => set({ plate_thickness: val })}
								/>
								<div className="flex flex-col gap-1">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Hole Pattern
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={d.pattern || "inline"}
										onChange={(e) => set({ pattern: e.target.value })}
									>
										<option value="inline">Inline</option>
										<option value="staggered">Staggered</option>
									</select>
								</div>
								<div className="flex flex-col gap-1">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Hole Discharge Model (Cd)
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={d.correlation || "fixed"}
										onChange={(e) => set({ correlation: e.target.value })}
									>
										<option value="fixed">Fixed (Florschuetz 0.79)</option>
										<option value="IdelchikThick">
											Deep Hole in Wall (Idelchik 4-18a)
										</option>
										<option value="Lichtarowicz">
											Long Orifice, l/d 2-10 (Lichtarowicz 1965)
										</option>
										<option value="McGreehanSchotsch">
											Cooling Hole (M-S 1988, plenum-fed)
										</option>
									</select>
								</div>
								{(d.correlation || "fixed") === "fixed" &&
									ratio("Cd", "Fixed Cd Value", 0.79)}
								<div className="flex flex-col gap-1">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Crossflow on the Holes (Cd)
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={d.cd_crossflow || "rogers_hersh_1976"}
										onChange={(e) => set({ cd_crossflow: e.target.value })}
									>
										<option value="rogers_hersh_1976">
											Grazing flow (Rogers &amp; Hersh 1976)
										</option>
										<option value="florschuetz_metzger_su_1982">
											CR-3630 Table 4.1 (tested plates only)
										</option>
										<option value="none">None (constant Cd)</option>
									</select>
									<div className="text-[9px] text-gray-400">
										A strong crossflow bends the jets over and lowers the Cd.
										Below Gc/Gj ~1.5 the hole's own Cd holds.
									</div>
								</div>
								{ratio("Nu_multiplier", "Nu Multiplier (rig matching)", 1.0)}
								{hasBypass && (
									<>
										<LengthInput
											id={`bypass_length_${selectedNode.id}`}
											label="Bypass Lead-in to Row 1"
											value={d.bypass_length ?? 0.5 * (d.xn_d ?? 5.0) * dJet}
											onChange={(val) => set({ bypass_length: val })}
										/>
										{ratio("bypass_entry_K", "Bypass Entry K", 0.0)}
										<div className="text-[9px] text-gray-400 -mt-1">
											0 for a bypass arriving along a duct, 0.5 for an abrupt
											sharp-edged entry from a large plenum. Unused when the gap
											continues another array's.
										</div>
									</>
								)}
								{rows && rows.length > 0 && (
									<div className="flex flex-col gap-1 mt-2">
										<label className="text-[10px] font-bold text-blue-500 uppercase tracking-wider">
											Per Row
										</label>
										<table className="text-[10px] font-mono">
											<thead>
												<tr className="text-gray-500">
													<th className="text-left">row</th>
													<th className="text-right">Re_j</th>
													<th className="text-right">Gc/Gj</th>
													<th className="text-right">Nu</th>
													<th className="text-right">h</th>
												</tr>
											</thead>
											<tbody>
												{rows.map((nu, i) => (
													<tr key={`row_${i + 1}`}>
														<td>{i + 1}</td>
														<td className="text-right">
															{d.result.rows_Re_j[i].toFixed(0)}
														</td>
														<td className="text-right">
															{d.result.rows_Gc_Gj[i].toFixed(3)}
														</td>
														<td className="text-right">{nu.toFixed(1)}</td>
														<td className="text-right">
															{d.result.rows_htc[i].toFixed(0)}
														</td>
													</tr>
												))}
											</tbody>
										</table>
										{d.result.surface_extrapolated > 0 && (
											<div className="text-[10px] text-amber-600">
												Outside Florschuetz's range on at least one row.
											</div>
										)}
									</div>
								)}
							</>
						);
					})()}

				{selectedNode.type === "area_change" &&
					(() => {
						// Face areas come from the solver's resolve pass
						// (/resolve_preview, #559): F0 from the upstream side, F1 from
						// the downstream side, exactly as the solve will use them.
						const geomF0 = inheritedField(previewItems, selectedNode.id, "F0");
						const geomF1 = inheritedField(previewItems, selectedNode.id, "F1");
						const effectiveF0: number | undefined =
							selectedNode.data.F0 ?? geomF0?.value ?? undefined;
						const effectiveF1: number | undefined =
							selectedNode.data.F1 ?? geomF1?.value ?? undefined;
						return (
							<>
								<div className="flex flex-col gap-2">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Model Type
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={selectedNode.data.model_type || "sharp"}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												model_type: e.target.value,
											})
										}
									>
										<option value="sharp">Sharp-Edged (Sudden)</option>
										<option value="conical">Conical (Gradual)</option>
									</select>
								</div>

								<GeometryFields
									idPrefix={`F0_${selectedNode.id}`}
									title="Upstream face (F0)"
									stores="A"
									userValue={selectedNode.data.F0}
									inherited={geomF0}
									onChange={(v) => updateNodeData(selectedNode.id, { F0: v })}
									missingHint="Nothing upstream to inherit F0 from — set D or A."
								/>
								<GeometryFields
									idPrefix={`F1_${selectedNode.id}`}
									title="Downstream face (F1)"
									stores="A"
									userValue={selectedNode.data.F1}
									inherited={geomF1}
									onChange={(v) => updateNodeData(selectedNode.id, { F1: v })}
									missingHint="Nothing downstream to inherit F1 from — set D or A."
								/>

								{(() => {
									// Only once both faces are known (set, or resolved by the
									// solver); an unresolved face has no ratio to show.
									if (effectiveF0 == null || effectiveF1 == null) return null;
									const ratio =
										effectiveF0 > 0 ? effectiveF1 / effectiveF0 : null;
									const label =
										ratio == null
											? null
											: ratio > 1.005
												? "expansion"
												: ratio < 0.995
													? "contraction"
													: "straight";
									return ratio != null ? (
										<div className="rounded bg-stone-50 border border-stone-200 px-3 py-2 text-[10px] text-stone-600 flex flex-col gap-0.5">
											<div className="flex justify-between">
												<span className="font-medium text-stone-500 uppercase tracking-wide text-[8px]">
													Area ratio
												</span>
												<span
													className={
														label === "expansion"
															? "text-blue-600 font-semibold"
															: label === "contraction"
																? "text-amber-600 font-semibold"
																: "text-stone-500"
													}
												>
													{label}
												</span>
											</div>
											<div className="flex justify-between mt-0.5">
												<span className="text-stone-400">F₁/F₀</span>
												<span className="font-mono font-semibold">
													{ratio.toFixed(4)}
												</span>
											</div>
										</div>
									) : null;
								})()}

								{(() => {
									// #1023, from the solver's resolve pass: equal faces are no
									// area change (zero loss); an explicit face that contradicts
									// its pipe is refused.
									const passThrough =
										inheritedField(
											previewItems,
											selectedNode.id,
											"pass_through",
										)?.value === 1;
									const reason = unresolvedReason(
										previewItems,
										selectedNode.id,
									);
									const faceRefusal = reason?.includes("contradicts") === true;
									if (!passThrough && !faceRefusal) return null;
									return (
										<p className="text-[9px] text-amber-700">
											{passThrough
												? "F₀ = F₁: no area change — zero loss, the flow passes unchanged."
												: `The solver refuses this geometry: ${reason}`}
										</p>
									);
								})()}

								{selectedNode.data.model_type === "conical" && (
									<LengthInput
										id={`length_${selectedNode.id}`}
										label="Axial Length"
										value={selectedNode.data.length ?? 0.1}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { length: val })
										}
									/>
								)}

								{selectedNode.data.model_type !== "conical" && (
									<div className="flex flex-col gap-1 mt-1">
										<LengthInput
											id={`D_h_${selectedNode.id}`}
											label="Hydraulic Diameter (D_h)"
											value={selectedNode.data.D_h}
											placeholder="Circular (Auto)"
											onChange={(val) =>
												updateNodeData(selectedNode.id, { D_h: val })
											}
											onClear={() =>
												updateNodeData(selectedNode.id, { D_h: undefined })
											}
										/>
										<p className="text-[9px] text-gray-400 italic">
											{effectiveF0 != null && effectiveF1 != null ? (
												<>
													Defaults to circular diameter:{" "}
													{(
														Math.sqrt(
															(4 * Math.min(effectiveF0, effectiveF1)) /
																Math.PI,
														) * (unitPreferences.length === "mm" ? 1000 : 1)
													).toFixed(2)}{" "}
													{unitPreferences.length}
												</>
											) : (
												"Defaults to the circular diameter of the smaller face."
											)}
										</p>
									</div>
								)}

								<InitialGuessEditor node={selectedNode} />
							</>
						);
					})()}

				{selectedNode.type === "mpce_tee" &&
					(() => {
						const flowDirection = selectedNode.data.flow_direction ?? "branch";
						const isMerge = flowDirection === "merge";

						// The tee's areas are resolved by the backend builder
						// (_resolve_tee_areas) and reported by /resolve_preview (#559)
						// with their real source: a channel connected directly to the
						// port, the 0.01 m^2 default, or F_branch derived from F_C / psi.
						const userSetFC = selectedNode.data.F_C != null;
						const userSetFBranch = selectedNode.data.F_branch != null;
						const geomFC = inheritedField(previewItems, selectedNode.id, "F_C");
						const geomFBranch = inheritedField(
							previewItems,
							selectedNode.id,
							"F_branch",
						);
						const inheritedFC = userSetFC
							? undefined
							: (geomFC?.value ?? undefined);
						const inheritedFBranch = userSetFBranch
							? undefined
							: (geomFBranch?.value ?? undefined);

						// A tee always resolves (it has fallbacks), so these are always
						// numbers in the solve. Until the preview arrives, show the
						// builder's own fallbacks (F_C 0.01 m^2, F_branch = F_C / psi).
						const effectiveFC: number = userSetFC
							? (selectedNode.data.F_C as number)
							: (inheritedFC ?? 0.01);
						const effectiveFBranch: number = userSetFBranch
							? (selectedNode.data.F_branch as number)
							: (inheritedFBranch ??
								effectiveFC /
									Math.max((selectedNode.data.psi as number) ?? 1.0, 1e-12));
						const psi =
							effectiveFBranch > 0 ? effectiveFC / effectiveFBranch : null;
						// Constant-K: a blank K shows (and solves with) Bassett's pair
						// at the design share and these resolved areas.
						const designShare: number =
							selectedNode.data.design_branch_share ?? 0.5;
						const [autoKStraight, autoKBranch] = bassettDesignK(
							isMerge ? "merge" : "branch",
							designShare,
							psi ?? 1,
							selectedNode.data.theta_deg ?? 90,
						);
						const netLoss = designNetLoss(
							designShare,
							selectedNode.data.K_straight ?? autoKStraight,
							selectedNode.data.K_branch ?? autoKBranch,
						);
						const dFromF = (f: number) => Math.sqrt((4 * f) / Math.PI);
						const lu = unitPreferences.length === "mm" ? 1000 : 1;
						const lu_label = unitPreferences.length;

						const result = selectedNode.data.result;

						return (
							<div className="flex flex-col gap-4">
								{/* Model explainer: the canvas calls this element just
								    "Tee"; the underlying model name lives here. */}
								<div className="text-[9px] leading-snug text-stone-500 italic">
									Momentum-based multi-port chamber junction (MPCE): a
									momentum-control-volume tee model (Mynard closure) for
									splitting ("branch") and joining ("merge") flow.
								</div>
								{/* Port legend (MPCE indigo / rose theme) */}
								<div className="rounded border border-stone-200 bg-stone-50 px-3 py-2 text-[9px] leading-snug text-stone-600">
									<div className="font-bold uppercase text-stone-400 mb-1">
										Port Guide
									</div>
									<div className="grid grid-cols-3 gap-1 text-center font-mono mb-1.5">
										<div>
											<span
												className="font-extrabold"
												style={{ color: isMerge ? "#3b82f6" : "#f43f5e" }}
											>
												S
											</span>{" "}
											straight
										</div>
										<div>
											<span
												className="font-extrabold"
												style={{ color: isMerge ? "#f43f5e" : "#3b82f6" }}
											>
												C
											</span>{" "}
											common
										</div>
										<div>
											<span
												className="font-extrabold"
												style={{ color: isMerge ? "#3b82f6" : "#f43f5e" }}
											>
												B
											</span>{" "}
											branch
										</div>
									</div>
									<div className="text-stone-500">
										{isMerge ? (
											<>
												<span style={{ color: "#3b82f6" }}>S, B</span> →{" "}
												<span style={{ color: "#f43f5e" }}>C</span>
												{"  "}(two inlets, one outlet at C)
											</>
										) : (
											<>
												<span style={{ color: "#3b82f6" }}>C</span> →{" "}
												<span style={{ color: "#f43f5e" }}>S, B</span>
												{"  "}(one inlet at C, two outlets)
											</>
										)}
									</div>
								</div>

								{/* Flow direction (constrained — the Mynard residual asserts) */}
								<div className="flex flex-col gap-2">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Flow Direction
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={flowDirection}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												flow_direction: e.target.value,
											})
										}
									>
										<option value="branch">Branch (1 inlet, 2 outlets)</option>
										<option value="merge">Merge (2 inlets, 1 outlet)</option>
									</select>
									<p className="text-[9px] text-gray-400 italic">
										{(selectedNode.data.junction_model ?? "mynard") ===
										"ideal_momentum"
											? "Wiring only: the ideal momentum closure takes each port's role from its flow, so any port may run reversed."
											: "Constrained: the residual asserts at solve time and errors if the observed flow direction disagrees with this declaration."}
									</p>
								</div>

								{/* Junction closure model */}
								<div className="flex flex-col gap-2">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Junction Model
									</label>
									<select
										className="p-2 border rounded bg-white text-xs border-stone-200"
										value={selectedNode.data.junction_model ?? "mynard"}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												junction_model: e.target.value,
											})
										}
									>
										<option value="mynard">Momentum CV (Mynard)</option>
										<option value="constant_k">Constant K (handbook)</option>
										<option value="ideal_momentum">
											Ideal momentum (parameter-free)
										</option>
									</select>
									{selectedNode.data.junction_model === "constant_k" && (
										<p className="text-[9px] text-gray-400 italic">
											Constant K: fixed loss coefficients referenced to the
											common-leg dynamic head (Idelchik convention); the Mynard
											closure is disabled and K does not vary with the flow
											split.
										</p>
									)}
									{selectedNode.data.junction_model === "ideal_momentum" && (
										<p className="text-[9px] text-gray-400 italic">
											Ideal momentum: a baseline with nothing tuned. Suppliers
											share one junction static pressure; each receiving port
											gets the momentum it can draw from the supplied streams,
											the misaligned part is lost (merge = Blaisdell &amp;
											Manson, split = vector Borda-Carnot). Expect higher losses
											than measured on dividing tees (Bassett: K straight +0.26,
											K branch +0.23). Valid for subsonic supplying ports only.
										</p>
									)}
								</div>

								{selectedNode.data.junction_model === "ideal_momentum" && (
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											Loss Factor [-]
										</label>
										<NumericInput
											value={selectedNode.data.loss_factor ?? 1.0}
											onChange={(val) =>
												updateNodeData(selectedNode.id, { loss_factor: val })
											}
											className="p-2 border rounded"
											placeholder="1.0"
											min={0.25}
										/>
										<p className="text-[9px] text-gray-400 italic">
											Scales each receiving port&apos;s loss against the mixed
											supplied stream; the mixing of several suppliers is not
											scaled. 1.0 is the parameter-free model. A knob for
											matching your rig, not a validated setting: nothing here
											recommends a value. Minimum 0.25, below which the solve
											loses robustness and a nearly lossless split creates flow
											work.
										</p>
									</div>
								)}

								{(selectedNode.data.junction_model ?? "mynard") ===
									"constant_k" && (
									<>
										<div className="flex flex-col gap-2">
											<label className="text-xs font-bold text-gray-500 uppercase">
												Design Branch Share [-]
											</label>
											<NumericInput
												value={designShare}
												onChange={(val) =>
													updateNodeData(selectedNode.id, {
														design_branch_share: val,
													})
												}
												className="p-2 border rounded"
												placeholder="0.5"
												min={0}
												max={1}
											/>
											<p className="text-[9px] text-gray-400 italic">
												Branch share of the common-leg flow the K values are
												meant for. A blank K comes from Bassett 2001&apos;s
												analytical pair at this share, the arm areas and the
												angle; a constant K is right only near its design split.
											</p>
										</div>
										<div className="flex flex-col gap-2">
											<label className="text-xs font-bold text-gray-500 uppercase">
												K Straight [-]
											</label>
											<NumericInput
												value={selectedNode.data.K_straight ?? null}
												placeholder={autoKStraight.toFixed(3)}
												onChange={(val) =>
													updateNodeData(selectedNode.id, { K_straight: val })
												}
												onClear={() =>
													updateNodeData(selectedNode.id, { K_straight: null })
												}
												className="p-2 border rounded"
											/>
										</div>
										<div className="flex flex-col gap-2">
											<label className="text-xs font-bold text-gray-500 uppercase">
												K Branch [-]
											</label>
											<NumericInput
												value={selectedNode.data.K_branch ?? null}
												placeholder={autoKBranch.toFixed(3)}
												onChange={(val) =>
													updateNodeData(selectedNode.id, { K_branch: val })
												}
												onClear={() =>
													updateNodeData(selectedNode.id, { K_branch: null })
												}
												className="p-2 border rounded"
											/>
											<p className="text-[9px] text-gray-400 italic">
												Pt_arm = Pt_jct -/+ K * q_common (separating/joining).
												Blank = from the design share (greyed value); type to
												override, clear to return to it.
											</p>
										</div>
										{netLoss < 0 && (
											<p className="text-[10px] leading-snug text-red-700 bg-red-50 border border-red-200 rounded px-2 py-1">
												This pair loses net total pressure {netLoss.toFixed(3)}{" "}
												&lt; 0 at the design share: the tee would create total
												pressure, so the solve refuses it. Enter both K values
												(for example measured ones), move the design share, or
												use the Mynard model. Bassett&apos;s analytical pair
												runs low for joining tees at shallow angles.
											</p>
										)}
									</>
								)}

								{/* Branch angle */}
								<div className="flex flex-col gap-2">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Branch Angle (deg)
									</label>
									<NumericInput
										value={selectedNode.data.theta_deg ?? 90}
										onChange={(val) =>
											updateNodeData(selectedNode.id, { theta_deg: val })
										}
										className="p-2 border rounded"
										placeholder="90"
										min={-90}
										max={90}
									/>
									<p className="text-[9px] text-gray-400 italic">
										Branch arm angle measured from the main axis, 0 &lt; |theta|
										&lt;= 90 deg; a negative angle is the mirror image.
									</p>
								</div>

								{/* ── Arm areas (3-arm with Mynard labels) ── */}

								{/* Common arm (F_C) */}
								<GeometryFields
									idPrefix={`F_C_mpce_${selectedNode.id}`}
									title="Common arm (F_C)"
									stores="A"
									userValue={selectedNode.data.F_C}
									inherited={geomFC ?? { value: 0.01, source: "default" }}
									onChange={(v) => updateNodeData(selectedNode.id, { F_C: v })}
									defaultHint="No channel connected directly to the common or straight port — the 0.01 m² default is used. Untick Inherit to set F_C."
								/>

								{/* Straight arm (read-only = F_C — Mynard joining type-6 sets F_S = F_C) */}
								<div className="flex flex-col gap-1">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Straight Arm Area (F_S)
									</label>
									<div className="rounded bg-stone-50 border border-stone-200 px-3 py-2 text-[10px] text-stone-500 flex justify-between items-center">
										<span className="font-mono">
											{effectiveFC.toExponential(3)} m²{"  "}(D ={" "}
											{(dFromF(effectiveFC) * lu).toFixed(2)} {lu_label})
										</span>
										<span className="italic text-[9px] text-stone-400">
											= F_C
										</span>
									</div>
									<p className="text-[9px] text-stone-400 italic">
										Every junction model takes F_S = F_C (Mynard's
										pseudosupplier geometry requires it). A channel on any port
										whose bore differs by more than 5% meets the port through a
										sharp contraction or expansion, inserted automatically and
										reported under Port steps.
									</p>
								</div>

								{/* Branch arm (F_branch) */}
								<GeometryFields
									idPrefix={`F_branch_mpce_${selectedNode.id}`}
									title="Branch arm (F_B)"
									stores="A"
									userValue={selectedNode.data.F_branch}
									inherited={
										geomFBranch ?? {
											value:
												inheritedFBranch ??
												effectiveFC /
													Math.max(
														(selectedNode.data.psi as number) ?? 1.0,
														1e-12,
													),
											source: "derived from F_C / psi",
										}
									}
									onChange={(v) =>
										updateNodeData(selectedNode.id, { F_branch: v })
									}
								/>

								{/* Derived ratio summary */}
								{psi != null && (
									<div className="rounded bg-stone-50 border border-stone-200 px-3 py-2 text-[10px] text-stone-600 flex flex-col gap-0.5">
										<div className="flex justify-between">
											<span className="font-medium text-stone-500 uppercase tracking-wide text-[8px]">
												Area ratio
											</span>
											<span
												className={
													psi > 1.005
														? "text-amber-600 font-semibold"
														: psi < 0.995
															? "text-blue-600 font-semibold"
															: "text-stone-500"
												}
											>
												{psi > 1.005
													? "C larger than B"
													: psi < 0.995
														? "B larger than C"
														: "equal areas"}
											</span>
										</div>
										<div className="flex justify-between mt-0.5">
											<span className="text-stone-400">psi = F_C / F_B</span>
											<span className="font-mono font-semibold">
												{psi.toFixed(4)}
											</span>
										</div>
									</div>
								)}

								{/* Initial guess panel */}
								<InitialGuessEditor node={selectedNode} />

								{/* Post-solve diagnostics — emitted by MPCEv2Element.diagnostics */}
								{result && (
									<div className="mt-4 p-3 bg-stone-50 border border-stone-200 rounded">
										<h3 className="text-xs font-bold text-stone-500 uppercase mb-3">
											{isMerge ? "Joining" : "Separating"} Flow Diagnostics
										</h3>
										<div className="grid grid-cols-2 gap-x-2 gap-y-3">
											{/* Junction total (Pt_jct) and static (P_jct) pressure (#773) */}
											{(
												[
													["Pt_jct", "Pt junction"],
													["P_jct", "P junction (static)"],
												] as const
											).map(([key, label]) =>
												result[key] != null ? (
													<div key={key} className="flex flex-col">
														<span className="text-stone-400 text-[9px] font-bold uppercase">
															{label}
														</span>
														<span className="font-mono text-xs font-bold">
															{Number(result[key]).toFixed(1)}{" "}
															<span className="text-[9px] font-normal text-stone-400">
																Pa
															</span>
														</span>
													</div>
												) : null,
											)}
											{/* Per-port mdots (port_0_m_dot, etc.). Indices depend on
											    flow direction:
											      branch: 0=com, 1=straight, 2=branch
											      merge:  0=straight, 1=branch, 2=com */}
											<div className="flex flex-col">
												<span className="text-stone-400 text-[9px] font-bold uppercase">
													ṁ common
												</span>
												<span className="font-mono text-xs font-bold">
													{Math.abs(
														(isMerge
															? result.port_2_m_dot
															: result.port_0_m_dot) ?? 0,
													).toFixed(4)}{" "}
													<span className="text-[9px] font-normal text-stone-400">
														kg/s
													</span>
												</span>
											</div>
											<div className="flex flex-col">
												<span className="text-stone-400 text-[9px] font-bold uppercase">
													ṁ straight
												</span>
												<span className="font-mono text-xs font-bold">
													{Math.abs(
														(isMerge
															? result.port_0_m_dot
															: result.port_1_m_dot) ?? 0,
													).toFixed(4)}{" "}
													<span className="text-[9px] font-normal text-stone-400">
														kg/s
													</span>
												</span>
											</div>
											<div className="flex flex-col">
												<span className="text-stone-400 text-[9px] font-bold uppercase">
													ṁ branch
												</span>
												<span className="font-mono text-xs font-bold">
													{Math.abs(
														(isMerge
															? result.port_1_m_dot
															: result.port_2_m_dot) ?? 0,
													).toFixed(4)}{" "}
													<span className="text-[9px] font-normal text-stone-400">
														kg/s
													</span>
												</span>
											</div>
											<div className="flex flex-col">
												<span className="text-stone-400 text-[9px] font-bold uppercase">
													ṁ branch / ṁ common
												</span>
												<span className="font-mono text-xs font-bold">
													{(result.mass_flow_ratio ?? 0).toFixed(3)}{" "}
													<span className="text-[9px] font-normal text-stone-400">
														kg/kg
													</span>
												</span>
											</div>
											{/* K coefficients — names differ by direction */}
											{isMerge ? (
												<>
													<div className="flex flex-col">
														<span className="text-stone-400 text-[9px] font-bold uppercase">
															K11 (str → com)
														</span>
														<span className="font-mono text-xs font-bold">
															{(result.K11 ?? 0).toFixed(3)}
														</span>
													</div>
													<div className="flex flex-col">
														<span className="text-stone-400 text-[9px] font-bold uppercase">
															K12 (bra → com)
														</span>
														<span className="font-mono text-xs font-bold">
															{(result.K12 ?? 0).toFixed(3)}
														</span>
													</div>
												</>
											) : (
												<>
													<div className="flex flex-col">
														<span className="text-stone-400 text-[9px] font-bold uppercase">
															K straight
														</span>
														<span className="font-mono text-xs font-bold">
															{(result.K_straight ?? 0).toFixed(3)}
														</span>
													</div>
													<div className="flex flex-col">
														<span className="text-stone-400 text-[9px] font-bold uppercase">
															K branch
														</span>
														<span className="font-mono text-xs font-bold">
															{(result.K_branch ?? 0).toFixed(3)}
														</span>
													</div>
												</>
											)}
										</div>
										{/* Port area vs the connected channel (#536) */}
										{(result.port_area_mismatch ?? 0) > 0 && (
											<div className="mt-3 text-[10px] text-amber-600">
												Area step not modelled at{" "}
												{(isMerge
													? ["straight", "branch", "common"]
													: ["common", "straight", "branch"]
												)
													.map((name, i) => ({
														name,
														r: result[`port_${i}_connected_area_ratio`] as
															| number
															| undefined,
														a: result[`port_${i}_area`] as number | undefined,
														o: result[`port_${i}_outer_area`] as
															| number
															| undefined,
													}))
													.filter(
														(p) => p.r != null && Math.abs(p.r - 1) > 0.05,
													)
													.map(
														(p) =>
															`${p.name} (tee D ${(dFromF(p.a ?? 0) * lu).toFixed(1)}, channel D ${(dFromF(p.o ?? 0) * lu).toFixed(1)} ${lu_label})`,
													)
													.join(", ")}
												. The junction model uses the tee's area; the step to
												the channel carries no contraction or expansion loss.
												The ideal momentum model can then report a
												total-pressure gain and the solve is rejected; size the
												channel to the tee, or set F_C / F_branch.
											</div>
										)}
										{/* Sharp steps the solver inserted where a channel's bore
										    differs from the port's (#536). */}
										{(["common", "straight", "branch"] as const).some(
											(p) => result[`step_${p}_zeta`] != null,
										) && (
											<div className="mt-2 rounded border border-amber-200 bg-amber-50 px-2 py-1.5 text-[9px] leading-snug text-amber-800">
												<div className="font-bold uppercase mb-0.5">
													Port steps
												</div>
												{(["common", "straight", "branch"] as const)
													.filter((p) => result[`step_${p}_zeta`] != null)
													.map((p) => (
														<div key={p} className="font-mono">
															{p}: A_channel / A_port{" "}
															{(
																result[
																	`step_${p}_connected_area_ratio`
																] as number
															).toFixed(2)}
															, ζ{" "}
															{(result[`step_${p}_zeta`] as number).toFixed(3)},
															ΔPt{" "}
															{(result[`step_${p}_dPt`] as number).toFixed(0)}{" "}
															Pa
														</div>
													))}
												<div className="italic text-amber-700 mt-0.5">
													The connected channel's bore differs from the port's,
													so the two meet through a sharp contraction or
													expansion (Weisbach / Borda-Carnot) inserted at the
													port.
												</div>
											</div>
										)}
									</div>
								)}
							</div>
						);
					})()}

				{selectedNode.type === "vortex" && (
					<div className="flex flex-col gap-4">
						<div className="flex flex-col gap-1">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Rotational Speed
							</label>
							<div className="flex items-center gap-1">
								<NumericInput
									id={`shaft_speed_rpm_${selectedNode.id}`}
									className="p-1 border rounded text-xs bg-white h-7 outline-none focus:ring-1 focus:ring-stone-200 flex-1"
									value={
										selectedNode.data.shaft_speed_rpm != null
											? selectedNode.data.shaft_speed_rpm *
												(unitPreferences.rotSpeed === "Hz" ? 1 / 60 : 1)
											: null
									}
									onChange={(val) =>
										updateNodeData(selectedNode.id, {
											shaft_speed_rpm:
												val * (unitPreferences.rotSpeed === "Hz" ? 60 : 1),
										})
									}
									onClear={() =>
										updateNodeData(selectedNode.id, { shaft_speed_rpm: null })
									}
									min={0}
									placeholder={`global: ${(
										(solverSettings.shaft_speed_rpm ?? 0) *
											(unitPreferences.rotSpeed === "Hz" ? 1 / 60 : 1)
									).toFixed(unitPreferences.rotSpeed === "Hz" ? 3 : 0)}`}
								/>
								<select
									value={unitPreferences.rotSpeed}
									onChange={(e) =>
										setRotSpeedUnit(e.target.value as "rpm" | "Hz")
									}
									className="w-14 h-[30px] border rounded text-[9px] bg-white outline-none shrink-0"
								>
									<option value="rpm">rpm</option>
									<option value="Hz">Hz</option>
								</select>
							</div>
						</div>
						{(selectedNode.data.r_c != null || selectedNode.data.n != null) && (
							<div className="rounded border border-amber-300 bg-amber-50 px-2 py-1.5 text-[10px] leading-snug text-amber-800 flex flex-col gap-1">
								<div className="font-bold">Vortex model changed (#661)</div>
								<div>
									Saved with the old model, which added a Vatistas pressure rise
									to the total pressure without any work input. It is now a free
									vortex: total pressure and temperature are conserved, and only
									the static pressure redistributes with radius. Results differ
									from the saved ones. Core radius and shape n are no longer
									used.
								</div>
								<button
									type="button"
									className="self-start rounded border border-amber-300 bg-white px-2 py-0.5 text-[10px] hover:bg-amber-100"
									onClick={() =>
										updateNodeData(
											selectedNode.id,
											{ r_c: undefined, n: undefined },
											true,
										)
									}
								>
									Dismiss
								</button>
							</div>
						)}
						<LengthInput
							id={`r_in_${selectedNode.id}`}
							label="Inlet Radius r_in"
							value={selectedNode.data.r_in ?? 0.05}
							onChange={(val) => updateNodeData(selectedNode.id, { r_in: val })}
						/>
						<LengthInput
							id={`r_out_${selectedNode.id}`}
							label="Outlet Radius r_out"
							value={selectedNode.data.r_out ?? 0.1}
							onChange={(val) =>
								updateNodeData(selectedNode.id, { r_out: val })
							}
						/>
						<div className="flex flex-col gap-1">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Inlet Swirl V<sub>θ</sub> at r_in [m/s]
							</label>
							<NumericInput
								id={`v_theta_in_${selectedNode.id}`}
								className="p-1 border rounded text-xs bg-white h-7 outline-none focus:ring-1 focus:ring-stone-200"
								value={selectedNode.data.v_theta_in ?? null}
								onChange={(val) =>
									updateNodeData(selectedNode.id, { v_theta_in: val })
								}
								onClear={() =>
									updateNodeData(selectedNode.id, { v_theta_in: null })
								}
								placeholder={`solid body ω·r_in: ${(
									(((selectedNode.data.shaft_speed_rpm ??
										solverSettings.shaft_speed_rpm ??
										0) *
										2 *
										Math.PI) /
										60) *
										(selectedNode.data.r_in ?? 0.05)
								).toFixed(1)}`}
							/>
							<span className="text-[9px] text-stone-400">
								Free vortex: r·V<sub>θ</sub> conserved, total pressure and
								temperature unchanged.
							</span>
						</div>
						{selectedNode.data.result && (
							<div className="rounded border border-stone-200 bg-stone-50 px-3 py-2 text-xs">
								<div className="font-bold uppercase text-stone-400 mb-1 text-[9px]">
									Result
								</div>
								<div className="flex justify-between">
									<span className="text-stone-500">
										V<sub>θ</sub> in → out
									</span>
									<span className="font-mono font-bold">
										{(selectedNode.data.result.v_theta_in ?? 0).toFixed(1)} →{" "}
										{(selectedNode.data.result.v_theta_out ?? 0).toFixed(1)} m/s
									</span>
								</div>
								{selectedNode.data.result.dP_swirl_static != null && (
									<div className="flex justify-between">
										<span className="text-stone-500">ΔP static (swirl)</span>
										<span className="font-mono font-bold">
											{(
												selectedNode.data.result.dP_swirl_static / 1000
											).toFixed(2)}{" "}
											kPa
										</span>
									</div>
								)}
								<div className="flex justify-between">
									<span className="text-stone-500">ṁ</span>
									<span className="font-mono font-bold">
										{(selectedNode.data.result.m_dot ?? 0).toFixed(4)} kg/s
									</span>
								</div>
							</div>
						)}
					</div>
				)}

				{selectedNode.type === "ejector" && (
					<div className="flex flex-col gap-4">
						{/* Throat and nozzle exit are internal geometry and never
						    inherit; the mixing duct A_3 is the same duct as the outlet
						    pipe, so it inherits from it (or the pipe from it) (#559). */}
						<GeometryFields
							idPrefix={`ej_t_${selectedNode.id}`}
							title="Nozzle throat (A_t)"
							stores="A"
							userValue={selectedNode.data.throat_area}
							inherited={{ value: 3.14e-5, source: "default" }}
							inheritable={false}
							onChange={(v) =>
								updateNodeData(selectedNode.id, { throat_area: v })
							}
						/>
						<GeometryFields
							idPrefix={`ej_p1_${selectedNode.id}`}
							title="Nozzle exit (A_p1)"
							stores="A"
							userValue={selectedNode.data.nozzle_exit_area}
							inherited={{ value: 1.0e-4, source: "default" }}
							inheritable={false}
							onChange={(v) =>
								updateNodeData(selectedNode.id, { nozzle_exit_area: v })
							}
						/>
						<GeometryFields
							idPrefix={`ej_3_${selectedNode.id}`}
							title="Mixing duct / outlet (A_3)"
							stores="A"
							userValue={selectedNode.data.mixing_area}
							inherited={
								inheritedField(
									previewItems,
									selectedNode.id,
									"mixing_area",
								) ?? {
									value: 8.0e-4,
									source: "default",
								}
							}
							onChange={(v) =>
								updateNodeData(selectedNode.id, { mixing_area: v })
							}
							defaultHint="No outlet pipe with a known diameter — the 8e-4 m² default is used. Connect a pipe, or untick Inherit to set A_3."
						/>
						{(() => {
							// #1023: a mixing area wider than the outlet pipe is refused
							// (a wider pipe is a diffuser exit and stays accepted).
							const reason = unresolvedReason(previewItems, selectedNode.id);
							return reason?.includes("outlet pipe") === true ? (
								<p className="text-[9px] text-amber-700">
									{`The solver refuses this geometry: ${reason}`}
								</p>
							) : null;
						})()}
						<div className="flex flex-col gap-1">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Recovery Efficiency
							</label>
							<NumericInput
								id={`recovery_efficiency_${selectedNode.id}`}
								className="p-1 border rounded text-xs bg-white h-7 outline-none focus:ring-1 focus:ring-stone-200"
								value={selectedNode.data.recovery_efficiency ?? 1.0}
								onChange={(val) =>
									updateNodeData(selectedNode.id, { recovery_efficiency: val })
								}
								min={0}
								max={1}
								placeholder="1.0"
							/>
							<span className="text-[9px] text-stone-400">
								Multiplies mixed stagnation P to give P_c*; 1.0 = no loss
							</span>
						</div>
						{solverSettings.global_regime !== "incompressible" && (
							<div className="flex flex-col gap-1">
								<label
									htmlFor={`anchor_closure_${selectedNode.id}`}
									className="text-xs font-bold text-gray-500 uppercase"
								>
									Back-Pressure Anchors
								</label>
								<select
									id={`anchor_closure_${selectedNode.id}`}
									className="p-1 border rounded bg-white text-xs h-7"
									value={selectedNode.data.anchor_closure || "huang1999"}
									onChange={(e) =>
										updateNodeData(selectedNode.id, {
											anchor_closure: e.target.value,
										})
									}
								>
									<option value="huang1999">
										Huang 1999 + sudden-expansion dead head (default)
									</option>
									<option value="kracik_dvorak">
										Kracik &amp; Dvorak (lossless mixing)
									</option>
								</select>
								<span className="text-[9px] text-stone-400">
									P_c* and the dead head P_b0 that end the plateau and the
									droop. Huang's pair scores 2-3x closer on measured air rigs;
									K&amp;D's lossless pair is kept for comparison.
								</span>
							</div>
						)}
						{solverSettings.global_regime === "incompressible" && (
							<div className="flex flex-col gap-1 rounded border border-stone-200 p-2">
								<div className="text-[10px] font-bold text-gray-500 uppercase">
									Jet-pump losses (incompressible)
								</div>
								<span className="text-[9px] text-stone-400">
									Sanger&apos;s K on each dynamic head. Empty nozzle / suction =
									the compressible nozzle efficiencies (K = 1/η − 1), so both
									regimes agree at low Mach.
								</span>
								{(
									[
										["k_nozzle", "K nozzle", "0.053"],
										["k_suction", "K suction", "0.176"],
										["k_throat", "K throat", "0"],
										["k_diffuser", "K diffuser", "0"],
									] as const
								).map(([key, label, ph]) => (
									<div key={key} className="flex items-center gap-2">
										<label
											htmlFor={`${key}_${selectedNode.id}`}
											className="text-[10px] text-gray-500 w-20 shrink-0"
										>
											{label}
										</label>
										<NumericInput
											id={`${key}_${selectedNode.id}`}
											className="p-1 border rounded text-xs bg-white h-7 outline-none focus:ring-1 focus:ring-stone-200 flex-1 min-w-0"
											value={selectedNode.data[key] ?? null}
											onChange={(val) =>
												updateNodeData(selectedNode.id, { [key]: val })
											}
											onClear={() =>
												updateNodeData(selectedNode.id, { [key]: null })
											}
											min={0}
											placeholder={ph}
										/>
									</div>
								))}
							</div>
						)}
						{selectedNode.data.result && (
							<div className="rounded border border-stone-200 bg-stone-50 px-3 py-2 text-xs">
								<div className="font-bold uppercase text-stone-400 mb-1 text-[9px]">
									Result
								</div>
								<div className="flex justify-between">
									<span className="text-stone-500">Entrainment ratio</span>
									<span className="font-mono font-bold">
										{(selectedNode.data.result.entrainment_ratio ?? 0).toFixed(
											4,
										)}
									</span>
								</div>
								{selectedNode.data.result.P_c_star != null && (
									<div className="flex justify-between">
										<span className="text-stone-500">P_c*</span>
										<span className="font-mono font-bold">
											{(selectedNode.data.result.P_c_star / 1000).toFixed(2)}{" "}
											kPa
										</span>
									</div>
								)}
								{selectedNode.data.result.head_ratio_N != null && (
									<div className="flex justify-between">
										<span className="text-stone-500">Head ratio N</span>
										<span className="font-mono font-bold">
											{selectedNode.data.result.head_ratio_N.toFixed(4)}
										</span>
									</div>
								)}
								{selectedNode.data.result.efficiency != null && (
									<div className="flex justify-between">
										<span className="text-stone-500">Efficiency M·N</span>
										<span className="font-mono font-bold">
											{(selectedNode.data.result.efficiency * 100).toFixed(1)} %
										</span>
									</div>
								)}
								{selectedNode.data.result.mach_jet != null && (
									<div className="flex justify-between">
										<span className="text-stone-500">Jet Mach</span>
										<span className="font-mono font-bold">
											{selectedNode.data.result.mach_jet.toFixed(3)}
										</span>
									</div>
								)}
							</div>
						)}
					</div>
				)}

				{selectedNode.type === "discrete_loss" &&
					(() => {
						// Upstream neighbour, used below for the theta auto-default.
						const upstreamEdge = edges.find(
							(e) => e.target === selectedNode.id && !e.data?.type,
						);
						const upstreamNode = upstreamEdge
							? nodes.find((n) => n.id === upstreamEdge.source)
							: undefined;

						// The reference area comes from the solver's resolve pass
						// (/resolve_preview, #559): nearest channel, then any area-bearing
						// neighbour, else the 0.1 m^2 chamber default.
						const userSetArea =
							selectedNode.data.area !== null &&
							selectedNode.data.area !== undefined;
						const geomArea = inheritedField(
							previewItems,
							selectedNode.id,
							"area",
						);
						const inheritedArea = userSetArea
							? undefined
							: (geomArea?.value ?? undefined);

						// Compute available has_theta nodes (combustors)
						const thetaNodes = nodes.filter(
							(n) => n.type === "combustor" && n.id !== selectedNode.id,
						);

						// Determine auto-default: first downstream combustor, then upstream
						const downstreamEdge = edges.find(
							(e) => e.source === selectedNode.id && !e.data?.type,
						);
						const downstreamNode = downstreamEdge
							? nodes.find((n) => n.id === downstreamEdge.target)
							: undefined;
						const autoDefault =
							downstreamNode?.type === "combustor"
								? downstreamNode.id
								: upstreamNode?.type === "combustor"
									? upstreamNode.id
									: null;

						const thetaSource: string | null =
							selectedNode.data.theta_source ?? null;
						const hasThetaSource =
							thetaSource !== null && thetaSource !== "none";
						const corrType: string =
							selectedNode.data.correlation_type ?? "constant_fraction";
						const isThetaCorr =
							corrType === "linear_theta_fraction" ||
							corrType === "linear_theta_head";

						const areaWarning =
							(userSetArea ? selectedNode.data.area : inheritedArea) != null &&
							(userSetArea ? selectedNode.data.area : inheritedArea) <= 0;

						return (
							<div className="flex flex-col gap-4">
								{/* Correlation type */}
								<div className="flex flex-col gap-2">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Correlation Type
									</label>
									<select
										className="p-2 border rounded bg-white text-sm"
										value={corrType}
										onChange={(e) => {
											const t = e.target.value;
											const needsTheta =
												t === "linear_theta_fraction" ||
												t === "linear_theta_head";
											if (
												needsTheta &&
												!hasThetaSource &&
												thetaNodes.length > 0
											) {
												updateNodeData(selectedNode.id, {
													correlation_type: t,
													theta_source: autoDefault ?? thetaNodes[0].id,
												});
											} else if (needsTheta && thetaNodes.length === 0) {
												return;
											} else {
												updateNodeData(selectedNode.id, {
													correlation_type: t,
												});
											}
										}}
									>
										<option value="constant_fraction">
											Constant fraction — dP/P = ξ
										</option>
										<option value="constant_head">
											Constant head — dP = ζ · q
										</option>
										{thetaNodes.length > 0 || hasThetaSource ? (
											<option value="linear_theta_fraction">
												Linear Θ fraction — dP/P = k·Θ + ξ₀
											</option>
										) : null}
										{thetaNodes.length > 0 || hasThetaSource ? (
											<option value="linear_theta_head">
												Linear Θ head — dP = (k·Θ + ζ₀) · q
											</option>
										) : null}
									</select>
									{isThetaCorr && !hasThetaSource && (
										<p className="text-[9px] text-amber-500 italic">
											Select a theta source below to enable this correlation.
										</p>
									)}
								</div>

								{/* Flow Area — always shown, inherited from upstream by default */}
								<GeometryFields
									idPrefix={`dl_${selectedNode.id}`}
									title="Flow cross-section"
									stores="A"
									userValue={selectedNode.data.area}
									inherited={geomArea ?? { value: 0.1, source: "default" }}
									onChange={(v) => updateNodeData(selectedNode.id, { area: v })}
									defaultHint="No neighbour with a known area — the 0.1 m² default is used. Connect a channel, or untick Inherit to set the area."
								/>
								{areaWarning && (
									<div className="text-[10px] text-amber-600 bg-amber-50 border border-amber-200 rounded px-2 py-1">
										⚠ Area ≤ 0 — head loss will be undefined.
									</div>
								)}

								{/* Parameters */}
								{(corrType === "constant_fraction" ||
									corrType === "linear_theta_fraction") && (
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											{corrType === "constant_fraction"
												? "ξ — loss fraction [−]"
												: "ξ₀ — cold-flow base fraction [−]"}
										</label>
										<NumericInput
											value={
												corrType === "constant_fraction"
													? (selectedNode.data.xi ?? 0.03)
													: (selectedNode.data.xi0 ?? 0.02)
											}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													[corrType === "constant_fraction" ? "xi" : "xi0"]:
														val,
												})
											}
											min={0}
											className="p-2 border rounded"
										/>
										{corrType === "linear_theta_fraction" && (
											<>
												<label className="text-xs font-bold text-gray-500 uppercase">
													k — Theta sensitivity [−/−]
												</label>
												<NumericInput
													value={selectedNode.data.k ?? 0.001}
													onChange={(val) =>
														updateNodeData(selectedNode.id, { k: val })
													}
													min={0}
													className="p-2 border rounded"
												/>
												<p className="text-[9px] text-gray-400 italic">
													P_out = P_in · (1 − k·Θ − ξ₀). Θ = T_burned/T_in − 1
												</p>
											</>
										)}
										{corrType === "constant_fraction" && (
											<p className="text-[9px] text-gray-400 italic">
												P_out = P_in · (1 − ξ). e.g. 0.03 = 3%
											</p>
										)}
										{/* Design flow (#550): a fraction loss quoted at its
										    design point scales with the dynamic head below it,
										    so the flow can run low or reverse. */}
										<label
											htmlFor={`m_ref_${selectedNode.id}`}
											className="text-xs font-bold text-gray-500 uppercase"
										>
											Design flow ṁ_ref [kg/s]
										</label>
										<NumericInput
											id={`m_ref_${selectedNode.id}`}
											value={selectedNode.data.m_ref ?? null}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													m_ref: val > 0 ? val : null,
												})
											}
											onClear={() =>
												updateNodeData(selectedNode.id, { m_ref: null })
											}
											min={0}
											format={formatSignificant}
											placeholder="none: constant at every flow"
											className="p-2 border rounded"
										/>
										{/* Opt-in (#1004): copy the last converged flow. At
										    ṁ = ṁ_ref the loss is exactly ξ, so the point stays. */}
										{(() => {
											const solved = solvedFlowAsMref(selectedNode.data);
											return (
												<button
													type="button"
													disabled={solved === null}
													onClick={() => {
														if (solved !== null) {
															updateNodeData(selectedNode.id, {
																m_ref: solved,
															});
														}
													}}
													title={
														solved === null
															? "Needs a converged solve in which this loss carries flow (and is not below a set ṁ_ref)."
															: "At ṁ = ṁ_ref the loss is exactly the fraction: the solved point is unchanged, and below it the loss now fades with the flow."
													}
													className="self-start px-2 py-0.5 border rounded bg-white text-[10px] hover:bg-stone-50 disabled:opacity-50 disabled:hover:bg-white"
												>
													Use the solved flow as ṁ_ref
													{solved !== null
														? ` (${solved.toPrecision(4)} kg/s)`
														: ""}
												</button>
											);
										})()}
										{selectedNode.data.result?.fraction_unanchored === 1 && (
											<div className="text-[10px] text-amber-600 bg-amber-50 border border-amber-200 rounded px-2 py-1">
												This loss sits at ~zero flow (|ṁ| below 1e-4 kg/s, the
												blend of its forward and reverse forms) in a network
												that is driven: with no ṁ_ref it holds ΔPt = ξ·Pt at
												every flow but zero, and the network gives it less. The
												loss law is not met at this point. Set ṁ_ref, its design
												flow: type it, or solve at a point where it flows and
												use the button above.
											</div>
										)}
										<p className="text-[9px] text-gray-400 italic">
											At and above ṁ_ref the loss is exactly the fraction; below
											it the loss scales with the dynamic head and vanishes at
											zero flow. Set it if the flow may run low or reverse: a
											constant fraction holds a pressure difference at zero
											flow, and the solver cannot carry a reversing flow through
											it.
										</p>
									</div>
								)}

								{(corrType === "constant_head" ||
									corrType === "linear_theta_head") && (
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											{corrType === "constant_head"
												? "ζ — Euler loss coeff. [−]"
												: "ζ₀ — cold-flow Euler coeff. [−]"}
										</label>
										<NumericInput
											value={
												corrType === "constant_head"
													? (selectedNode.data.zeta ?? 1.0)
													: (selectedNode.data.zeta0 ?? 1.0)
											}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													[corrType === "constant_head" ? "zeta" : "zeta0"]:
														val,
												})
											}
											min={0}
											className="p-2 border rounded"
										/>
										{corrType === "linear_theta_head" && (
											<>
												<label className="text-xs font-bold text-gray-500 uppercase">
													k — Theta sensitivity [−/−]
												</label>
												<NumericInput
													value={selectedNode.data.k ?? 0.001}
													onChange={(val) =>
														updateNodeData(selectedNode.id, { k: val })
													}
													min={0}
													className="p-2 border rounded"
												/>
												<p className="text-[9px] text-gray-400 italic">
													dP = (k·Θ + ζ₀) · ½ρv². Θ = T_burned/T_in − 1
												</p>
											</>
										)}
									</div>
								)}

								{/* Theta Source */}
								<div className="flex flex-col gap-2 border-t pt-3">
									<label className="text-xs font-bold text-gray-500 uppercase">
										Theta Source
									</label>
									<select
										className="p-2 border rounded bg-white text-sm"
										value={thetaSource ?? ""}
										onChange={(e) => {
											const val = e.target.value || null;
											updateNodeData(selectedNode.id, { theta_source: val });
										}}
									>
										<option value="">
											Auto{autoDefault ? ` (→ ${autoDefault})` : " (cold-flow)"}
										</option>
										<option value="none">None (cold-flow)</option>
										{thetaNodes.map((n) => (
											<option key={n.id} value={n.id}>
												{n.data?.label
													? `${n.data.label} (${n.id})`
													: `Combustor: ${n.id}`}
											</option>
										))}
									</select>
									{thetaNodes.length === 0 && (
										<p className="text-[9px] text-gray-400 italic">
											No combustor nodes in network. Only ξ/ζ correlations
											available.
										</p>
									)}
								</div>

								<div className="grid grid-cols-2 gap-3 border-t pt-3">
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											Nu Multiplier
										</label>
										<NumericInput
											value={selectedNode.data.Nu_multiplier ?? 1.0}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													Nu_multiplier: val,
												})
											}
											className="p-2 border rounded"
										/>
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.Nu_multiplier ?? 1}`}
										</span>
									</div>
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											F Multiplier
										</label>
										<NumericInput
											value={selectedNode.data.f_multiplier ?? 1.0}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													f_multiplier: val,
												})
											}
											className="p-2 border rounded"
										/>
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.f_multiplier ?? 1}`}
										</span>
									</div>
								</div>

								<div className="border-t pt-3">
									<SurfaceEnhancementInspector
										surface={selectedNode.data.surface || { type: "smooth" }}
										onChange={(surface) =>
											updateNodeData(selectedNode.id, { surface })
										}
									/>
								</div>

								<InitialGuessEditor node={selectedNode} />
							</div>
						);
					})()}

				{selectedNode.type === "combustor" &&
					(() => {
						// Area and Dh come from the solver's resolve pass
						// (/resolve_preview, #559): Dh from the nearest upstream channel,
						// area derived from it (or inferred from a neighbour), else the
						// 0.1 m^2 chamber default.
						const geomAreaComb = inheritedField(
							previewItems,
							selectedNode.id,
							"area",
						);
						const geomDhComb = inheritedField(
							previewItems,
							selectedNode.id,
							"Dh",
						);
						return (
							<div className="flex flex-col gap-4">
								<div className="flex flex-col gap-2">
									<label
										htmlFor={`method_${selectedNode.id}`}
										className="text-xs font-bold text-gray-500 uppercase"
									>
										Combustion Method
									</label>
									<select
										id={`method_${selectedNode.id}`}
										className="p-2 border rounded bg-white text-sm"
										value={selectedNode.data.method || "complete"}
										onChange={(e) =>
											updateNodeData(selectedNode.id, {
												method: e.target.value,
											})
										}
									>
										<option value="complete">Complete (CO2 + H2O)</option>
										<option value="equilibrium">
											Partial equilibrium (CO/H2 shift)
										</option>
									</select>
									<p className="text-[10px] text-gray-500 leading-snug">
										{(selectedNode.data.method || "complete") === "complete"
											? "All fuel burns to CO2 and H2O. Fast; the default."
											: "Complete combustion, then reforming and water-gas shift among CO, H2, CO2 and H2O. No dissociation (no O, OH or NO): the same as Complete when lean, cooler when rich. Usually closer to a real combustor than full equilibrium, which it rarely reaches."}
									</p>
								</div>

								<GeometryFields
									idPrefix={`comb_${selectedNode.id}`}
									title="Combustor cross-section"
									stores="A"
									userValue={selectedNode.data.area}
									inherited={geomAreaComb ?? { value: 0.1, source: "default" }}
									onChange={(v) => updateNodeData(selectedNode.id, { area: v })}
									dh={{
										userValue: selectedNode.data.Dh,
										inherited: geomDhComb,
										onChange: (v) => updateNodeData(selectedNode.id, { Dh: v }),
									}}
									defaultHint="Nothing upstream to size the combustor from — the 0.1 m² default is used. Connect a channel, or untick Inherit to set the area or Dh."
								/>
								<div className="grid grid-cols-2 gap-3">
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											Nu Multiplier
										</label>
										<NumericInput
											id={`Nu_multiplier_comb_${selectedNode.id}`}
											value={selectedNode.data.Nu_multiplier ?? 1.0}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													Nu_multiplier: val,
												})
											}
											className="p-2 border rounded"
										/>
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.Nu_multiplier ?? 1}`}
										</span>
									</div>
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											F Multiplier
										</label>
										<NumericInput
											id={`f_multiplier_comb_${selectedNode.id}`}
											value={selectedNode.data.f_multiplier ?? 1.0}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													f_multiplier: val,
												})
											}
											className="p-2 border rounded"
										/>
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.f_multiplier ?? 1}`}
										</span>
									</div>
								</div>
								<SurfaceEnhancementInspector
									surface={selectedNode.data.surface || { type: "smooth" }}
									onChange={(surface) =>
										updateNodeData(selectedNode.id, { surface })
									}
								/>

								<InitialGuessEditor node={selectedNode} />
							</div>
						);
					})()}

				{selectedNode.type === "momentum_chamber" &&
					(() => {
						// Area and Dh come from the solver's resolve pass
						// (/resolve_preview, #559): a merge chamber takes its main
						// inlet's area; otherwise Dh from the nearest upstream channel and
						// the area from it (or a neighbour), else the 0.1 m^2 default.
						const geomAreaMom = inheritedField(
							previewItems,
							selectedNode.id,
							"area",
						);
						const geomDhMom = inheritedField(
							previewItems,
							selectedNode.id,
							"Dh",
						);
						return (
							<div className="flex flex-col gap-4">
								<GeometryFields
									idPrefix={`mom_${selectedNode.id}`}
									title="Chamber cross-section"
									stores="A"
									userValue={selectedNode.data.area}
									inherited={geomAreaMom ?? { value: 0.1, source: "default" }}
									onChange={(v) => updateNodeData(selectedNode.id, { area: v })}
									dh={{
										userValue: selectedNode.data.Dh,
										inherited: geomDhMom,
										onChange: (v) => updateNodeData(selectedNode.id, { Dh: v }),
									}}
									defaultHint="Nothing upstream to size the chamber from — the 0.1 m² default is used. Connect a channel, or untick Inherit to set the area or Dh."
								/>
								{(() => {
									// #1023: narrower than an adjoining pipe, a contraction or
									// expansion with no loss. The solve warns; it does not refuse.
									const warning = previewWarning(previewItems, selectedNode.id);
									return warning ? (
										<p className="text-[9px] text-amber-700">{warning}</p>
									) : null;
								})()}
								<div className="grid grid-cols-2 gap-3">
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											Nu Multiplier
										</label>
										<NumericInput
											id={`Nu_multiplier_mom_${selectedNode.id}`}
											value={selectedNode.data.Nu_multiplier ?? 1.0}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													Nu_multiplier: val,
												})
											}
											className="p-2 border rounded"
										/>
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.Nu_multiplier ?? 1}`}
										</span>
									</div>
									<div className="flex flex-col gap-2">
										<label className="text-xs font-bold text-gray-500 uppercase">
											F Multiplier
										</label>
										<NumericInput
											id={`f_multiplier_mom_${selectedNode.id}`}
											value={selectedNode.data.f_multiplier ?? 1.0}
											onChange={(val) =>
												updateNodeData(selectedNode.id, {
													f_multiplier: val,
												})
											}
											className="p-2 border rounded"
										/>
										<span className="text-[9px] text-stone-400">
											{`global: ×${solverSettings.f_multiplier ?? 1}`}
										</span>
									</div>
								</div>
								<SurfaceEnhancementInspector
									surface={selectedNode.data.surface || { type: "smooth" }}
									onChange={(surface) =>
										updateNodeData(selectedNode.id, { surface })
									}
								/>

								<InitialGuessEditor node={selectedNode} />
							</div>
						);
					})()}

				{selectedNode.data.result && (
					<div className="mt-4 p-3 bg-stone-50 border border-stone-200 rounded">
						<div className="flex justify-between items-center mb-3">
							<h3 className="text-xs font-bold text-stone-500 uppercase flex items-center gap-2">
								<Activity size={14} />
								Live Telemetry
							</h3>
							{isCardKind(selectedNode.type) && (
								<CardRowsReset kind={selectedNode.type} />
							)}
						</div>
						{(() => {
							const why = unacceptedResult(
								selectedNode.data.result,
								solveResults,
							);
							return (
								why && (
									<p className="text-[10px] text-red-500 mb-2">
										{unacceptedNote(why)}
									</p>
								)
							);
						})()}

						{displaySettings.length === 0 ? (
							<p className="text-[10px] text-stone-400 italic text-center py-2">
								Select fields in the Display Fields panel to show diagnostics
								here.
							</p>
						) : (
							<div className="flex flex-col gap-3">
								<div className="grid grid-cols-2 gap-x-2 gap-y-3">
									{displaySettings.map((key) => {
										const val = resolveField(selectedNode.data.result, key);
										if (val === undefined) return null;
										const meta = quantityMeta(key);
										// A junction's port_<i>_* keys carry the arm letter
										// the node card draws (utils/portArms.ts).
										const label = withPortArm(
											meta.label,
											key,
											selectedNode.type,
											selectedNode.data,
										);
										return (
											<div key={key} className="flex flex-col">
												<span className="flex items-center gap-1">
													<span
														className="text-stone-400 text-[9px] font-bold"
														title={key}
													>
														{label}
													</span>
													{isCardKind(selectedNode.type) && (
														<CardRowPin
															kind={selectedNode.type}
															fieldKey={key}
														/>
													)}
												</span>
												<span className="font-mono text-xs font-bold whitespace-nowrap">
													{meta.format(val)}{" "}
													<span className="text-[9px] font-normal text-stone-400 ml-0.5">
														{meta.unit === "—" ? "" : meta.unit}
													</span>
												</span>
											</div>
										);
									})}
								</div>

								{/* Composition — not a scalar field, always shown when present */}
								{selectedNode.data.result.state?.X && (
									<div className="border-t border-stone-100 pt-3">
										<div className="text-[10px] font-bold text-stone-400 uppercase mb-2">
											Composition (Mole %)
										</div>
										<div className="max-h-40 overflow-y-auto pr-1 flex flex-col gap-1">
											{selectedNode.data.result.state.X.map(
												(x: number, i: number) => ({
													x,
													name: speciesMetadata?.names[i] || `Species ${i}`,
												}),
											)
												.filter((item: any) => item.x > 1e-6)
												.map((item: any, i: number) => (
													<div
														key={i}
														className="flex justify-between font-mono text-[9px] bg-white p-1 rounded border border-stone-50"
													>
														<span className="text-stone-500 text-ellipsis overflow-hidden shrink pr-2">
															{item.name}
														</span>
														<span className="font-bold flex-shrink-0">
															{(item.x * 100).toFixed(4)}%
														</span>
													</div>
												))}
										</div>
									</div>
								)}
							</div>
						)}
					</div>
				)}
			</aside>
		);
	}

	// Edge Inspector
	if (selectedEdge && selectedEdge.data?.type === "thermal") {
		const edgeQ: number | undefined = selectedEdge.data.result?.Q;
		const flowsForward = edgeQ === undefined || edgeQ >= 0;

		let probeTemp: number | null = null;
		let probeLabel = "Probe Temp";
		if (selectedEdge.data.result?.T_interface) {
			const targetX_manual = selectedEdge.data.probe_depth ?? 0;
			const tInt = selectedEdge.data.result.T_interface as number[];

			// Fallback layers for legacy/initial elements
			const layers = selectedEdge.data.layers || [
				{
					thickness: selectedEdge.data.thickness || 0.003,
					conductivity: selectedEdge.data.conductivity || 20.0,
				},
			];

			let targetX = targetX_manual;
			let hotSideShortcut = false;
			if (
				selectedEdge.data.probe_mode === "preset" &&
				selectedEdge.data.probe_preset
			) {
				const { type, index } = selectedEdge.data.probe_preset;
				// Ensure index is valid for current layers
				const safeIndex = Math.min(index, layers.length - 1);
				const layer = layers[safeIndex];

				let runningX = 0;
				for (let i = 0; i < safeIndex; i++) {
					runningX += layers[i].thickness;
				}

				const L = layer?.thickness || 0;
				if (type === "hot") {
					targetX = flowsForward ? runningX : runningX + L;
					probeLabel = `L${safeIndex + 1} ${flowsForward ? "Hot" : "Cold"} Side`;
				} else if (type === "avg") {
					targetX = runningX + L / 2;
					probeLabel = `L${safeIndex + 1} Average`;
				} else {
					targetX = flowsForward ? runningX + L : runningX;
					probeLabel = `L${safeIndex + 1} ${flowsForward ? "Cold" : "Hot"} Side`;
				}
			} else if (!selectedEdge.data.probe_mode) {
				hotSideShortcut = true;
				probeLabel = "Hot Side";
			} else {
				const scale = unitPreferences.length === "mm" ? 1000 : 1;
				probeLabel = `Custom Depth d=${(targetX * scale).toFixed(unitPreferences.length === "mm" ? 1 : 3)}${unitPreferences.length}`;
			}

			if (tInt.length === layers.length + 1) {
				if (hotSideShortcut) {
					probeTemp = flowsForward ? tInt[0] : tInt[tInt.length - 1];
				} else if (targetX <= 0) {
					probeTemp = tInt[0];
				} else {
					let found = false;
					let iterX = 0;
					for (let i = 0; i < layers.length; i++) {
						const t = layers[i].thickness;
						const nextX = iterX + t;
						if (targetX <= nextX + 1e-9) {
							const frac = t > 0 ? (targetX - iterX) / t : 0;
							probeTemp = tInt[i] + frac * (tInt[i + 1] - tInt[i]);
							found = true;
							break;
						}
						iterX = nextX;
					}
					if (!found) probeTemp = tInt[tInt.length - 1];
				}
			} else {
				probeLabel = "(solve needed)";
				probeTemp = null;
			}
		}

		// Calculate derived HTCs for diagnostics
		const layers = selectedEdge.data.layers || [];
		const R_wall =
			layers.reduce(
				(acc: number, l: any) =>
					acc + (l.thickness || 0) / (l.conductivity || 1),
				0,
			) + (selectedEdge.data.R_fouling || 0);

		const h_wall = R_wall > 0 ? 1 / R_wall : Infinity;
		const h_a = selectedEdge.data.result?.h_a ?? Infinity;
		const h_b = selectedEdge.data.result?.h_b ?? Infinity;

		const h_total_inv =
			(h_a > 0 ? 1 / h_a : 0) + R_wall + (h_b > 0 ? 1 / h_b : 0);
		const h_total = h_total_inv > 0 ? 1 / h_total_inv : 0;

		return (
			<aside
				key={selectedEdge.id}
				className="w-80 border-l bg-white p-4 flex flex-col gap-4 overflow-y-auto"
			>
				<h2 className="text-lg font-bold border-b pb-2 uppercase text-orange-600">
					Thermal Wall
				</h2>

				{/* Wall Name Label */}
				<div className="flex flex-col gap-1">
					<label className="text-xs font-bold text-gray-500 uppercase">
						Wall Name
					</label>
					<input
						type="text"
						value={selectedEdge.data.label ?? ""}
						onChange={(e) =>
							updateEdgeData(selectedEdge.id, { label: e.target.value }, true)
						}
						placeholder="e.g. Custom Wall"
						className="p-1.5 border rounded text-xs w-full focus:ring-1 focus:ring-orange-500 outline-none"
					/>
				</div>

				{/* Show label on graph */}
				<label className="flex items-center gap-2 text-xs text-gray-600 cursor-pointer select-none">
					<input
						type="checkbox"
						checked={!!selectedEdge.data?.show_label}
						onChange={(e) =>
							updateEdgeData(
								selectedEdge.id,
								{ show_label: e.target.checked },
								true,
							)
						}
						className="rounded"
					/>
					Show label in graph
				</label>
				{/* Per-side areas (#887): a wall whose two faces differ (a hole
				    in a thick wall, a fin, a curved wall). */}
				<label className="flex items-center gap-2 text-xs text-gray-600 cursor-pointer select-none">
					<input
						type="checkbox"
						checked={selectedEdge.data.area_a != null}
						onChange={(e) => {
							const A = selectedEdge.data.area ?? 0.05;
							updateEdgeData(
								selectedEdge.id,
								e.target.checked
									? { area_a: A, area_b: A }
									: { area_a: undefined, area_b: undefined },
							);
						}}
						className="rounded"
					/>
					Per-side areas (faces differ)
				</label>
				{selectedEdge.data.area_a != null ? (
					<>
						<AreaInput
							id={`area_a_edge_${selectedEdge.id}`}
							label="Area, side A (source)"
							value={selectedEdge.data.area_a}
							onChange={(val) =>
								updateEdgeData(selectedEdge.id, { area_a: val })
							}
						/>
						<AreaInput
							id={`area_b_edge_${selectedEdge.id}`}
							label="Area, side B (target)"
							value={selectedEdge.data.area_b ?? selectedEdge.data.area_a}
							onChange={(val) =>
								updateEdgeData(selectedEdge.id, { area_b: val })
							}
						/>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`wall_area_${selectedEdge.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Conduction area
							</label>
							<select
								id={`wall_area_${selectedEdge.id}`}
								value={selectedEdge.data.wall_area ?? "log_mean"}
								onChange={(e) =>
									updateEdgeData(selectedEdge.id, {
										wall_area: e.target.value,
									})
								}
								className="p-2 border rounded text-sm"
							>
								<option value="log_mean">
									Log mean (curved wall; exact for a cylinder)
								</option>
								<option value="slab_a">Slab on side A</option>
								<option value="slab_b">Slab on side B</option>
							</select>
						</div>
					</>
				) : (
					<AreaInput
						id={`area_edge_${selectedEdge.id}`}
						label="Contact Area"
						value={selectedEdge.data.area ?? 0.05}
						onChange={(val) => updateEdgeData(selectedEdge.id, { area: val })}
					/>
				)}

				<div className="flex flex-col gap-2">
					<label className="text-xs font-bold text-gray-500 uppercase">
						Fouling Resistance (m²K/W)
					</label>
					<NumericInput
						value={selectedEdge.data.R_fouling ?? 0}
						onChange={(val) =>
							updateEdgeData(selectedEdge.id, { R_fouling: val })
						}
						className="p-1.5 border rounded text-sm bg-stone-50 focus:bg-white"
					/>
				</div>

				{(() => {
					const arrayTypes = ["impingement_array", "effusion_liner"];
					const onArray = [selectedEdge.source, selectedEdge.target].some(
						(nid) =>
							arrayTypes.includes(nodes.find((n) => n.id === nid)?.type ?? ""),
					);
					const orientation = orientationOf(selectedEdge.data);
					const runs = [selectedEdge.source, selectedEdge.target].map((nid) =>
						sideRun(nid, nodes, edges),
					);
					const runText = (r: (typeof runs)[number]) =>
						`${r.from ?? "?"} → ${r.to ?? "?"}`;
					return (
						<div className="flex flex-col gap-1 border-t pt-4">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Orientation
							</label>
							{onArray ? (
								<p className="text-[10px] text-gray-500 leading-snug">
									An array's row walls exchange with the rows they face; the
									wall has no orientation to declare.
								</p>
							) : (
								<>
									<select
										aria-label="Wall orientation"
										value={orientation}
										onChange={(e) =>
											updateEdgeData(selectedEdge.id, {
												orientation: e.target.value,
											})
										}
										className="p-1.5 border rounded text-xs bg-white w-full border-stone-200 outline-none focus:ring-1 focus:ring-orange-500"
									>
										<option value="aligned">
											{orientationGlyph("aligned")} Aligned: both sides run the
											same way along the wall
										</option>
										<option value="opposed">
											{orientationGlyph("opposed")} Opposed: the sides run
											opposite ways
										</option>
									</select>
									<div className="text-[10px] text-gray-500 font-mono leading-snug">
										<div>
											A {runs[0].label}: {runText(runs[0])}
										</div>
										<div>
											B {runs[1].label}: {runText(runs[1])}
										</div>
									</div>
									<p className="text-[10px] text-gray-500 leading-snug">
										Forward flow on both sides runs{" "}
										{designArrangement(orientation)} flow. A reversed side flips
										it; the solve follows the flow signs.
									</p>
								</>
							)}
						</div>
					);
				})()}

				{(() => {
					// A film covering one side (#386): the holes upstream on the
					// gas path, named here; that side's T_aw becomes
					// T_aw - eta (T_aw - T_coolant).
					const providers = nodes.filter(
						(n) =>
							n.type === "film_cooling" ||
							n.type === "slot_film" ||
							n.type === "effusion_plate" ||
							n.type === "effusion_liner",
					);
					const d = selectedEdge.data;
					const r = d.result;
					const sideLabel = (nid: string) => {
						const n = nodes.find((m) => m.id === nid);
						return n?.data?.label || nid;
					};
					return (
						<div className="flex flex-col gap-2 border-t pt-4">
							<label className="text-xs font-bold text-gray-500 uppercase">
								Film
							</label>
							<select
								aria-label="Wall film"
								value={d.film ?? ""}
								onChange={(e) => {
									const film = e.target.value || null;
									updateEdgeData(selectedEdge.id, {
										film,
										...(film && {
											film_side: defaultFilmSide(selectedEdge, film, edges),
										}),
									});
								}}
								className="p-1.5 border rounded text-xs bg-white w-full border-stone-200 outline-none focus:ring-1 focus:ring-orange-500"
							>
								<option value="">None</option>
								{providers.map((n) => (
									<option key={n.id} value={n.id}>
										{n.data?.label || n.id} ({nodeTypeLabel(n.type ?? "")})
									</option>
								))}
							</select>
							{providers.length === 0 && (
								<p className="text-[10px] text-gray-500 leading-snug">
									Add Film Cooling Holes (or an Effusion Plate with a panel
									length, or an Effusion Liner) upstream on the gas path to
									cover this wall. Dragging from their "f" handle onto this
									wall's label sets it too.
								</p>
							)}
							{d.film && (
								<>
									<select
										aria-label="Film side"
										value={d.film_side ?? "a"}
										onChange={(e) =>
											updateEdgeData(selectedEdge.id, {
												film_side: e.target.value,
											})
										}
										className="p-1.5 border rounded text-xs bg-white w-full border-stone-200 outline-none focus:ring-1 focus:ring-orange-500"
									>
										<option value="a">
											Covers {sideLabel(selectedEdge.source)}
										</option>
										<option value="b">
											Covers {sideLabel(selectedEdge.target)}
										</option>
									</select>
									<LengthInput
										id={`film_start_${selectedEdge.id}`}
										label="Wall Starts After Last Row"
										value={d.film_start ?? 0}
										onChange={(val) =>
											updateEdgeData(selectedEdge.id, { film_start: val })
										}
									/>
									<LengthInput
										id={`film_end_${selectedEdge.id}`}
										label="Wall Ends After Last Row"
										value={d.film_end ?? d.film_start ?? 0}
										onChange={(val) =>
											updateEdgeData(selectedEdge.id, { film_end: val })
										}
									/>
									<p className="text-[10px] text-gray-500 leading-snug">
										The covered side must be the gas duct downstream of the
										holes. Only its driving temperature changes; its h and
										stream stay its own. Not scored as a coupled wall.
									</p>
									{r?.eta_film != null && (
										<div className="flex flex-col gap-0.5 bg-blue-50/50 p-2 rounded border border-blue-100/50 text-[10px] font-mono">
											<span>η film: {Number(r.eta_film).toFixed(3)}</span>
											<span>
												T_aw: {Number(r.T_aw_unfilmed).toFixed(0)} →{" "}
												{Number(
													d.film_side === "b" ? r.T_aw_b : r.T_aw_a,
												).toFixed(0)}{" "}
												K (coolant {Number(r.T_coolant_film).toFixed(0)} K)
											</span>
											<span>
												M: {Number(r.film_blowing_ratio).toFixed(2)} · DR:{" "}
												{Number(r.film_density_ratio).toFixed(2)}
											</span>
											{r.film_extrapolated > 0 && (
												<span className="text-amber-600">
													Outside the film correlation's fitted range
												</span>
											)}
										</div>
									)}
								</>
							)}
						</div>
					);
				})()}

				{(() => {
					// Radiation on each face (#793, #1071): side A is the edge's
					// source, side B its target, as for the per-side areas.
					const d = selectedEdge.data;
					const r = d.result;
					const sideLabel = (nid: string) => {
						const n = nodes.find((m) => m.id === nid);
						return n?.data?.label || nid;
					};
					const face = (side: "a" | "b") => {
						const q = r?.[`q_rad_${side}`];
						if (q == null || !Number.isFinite(q)) return null;
						const eps = r[`eps_g_${side}`];
						return (
							<span>
								{side.toUpperCase()}: q {(q / 1000).toFixed(2)} kW/m² · T
								surface {Number(r[`T_surface_${side}`]).toFixed(0)} K
								{eps != null && ` · ε gas ${Number(eps).toFixed(3)}`}
							</span>
						);
					};
					return (
						<div className="flex flex-col gap-3 border-t pt-4">
							<RadiationEditor
								id={`radiation_a_${selectedEdge.id}`}
								label={`Radiation, side A (${sideLabel(selectedEdge.source)})`}
								value={d.radiation_a}
								onChange={(radiation_a) =>
									updateEdgeData(selectedEdge.id, { radiation_a })
								}
								allowCasing
							/>
							<RadiationEditor
								id={`radiation_b_${selectedEdge.id}`}
								label={`Radiation, side B (${sideLabel(selectedEdge.target)})`}
								value={d.radiation_b}
								onChange={(radiation_b) =>
									updateEdgeData(selectedEdge.id, { radiation_b })
								}
								allowCasing
							/>
							{(face("a") || face("b")) && (
								<div className="flex flex-col gap-0.5 bg-blue-50/50 p-2 rounded border border-blue-100/50 text-[10px] font-mono">
									{face("a")}
									{face("b")}
								</div>
							)}
						</div>
					);
				})()}

				<div className="border-t pt-4">
					<WallLayersEditor
						layers={
							selectedEdge.data.layers || [
								{
									thickness: selectedEdge.data.thickness ?? 0.003,
									conductivity: selectedEdge.data.conductivity ?? 20.0,
									material: "generic",
								},
							]
						}
						onChange={(layers) => updateEdgeData(selectedEdge.id, { layers })}
					/>
				</div>

				<div className="flex flex-col gap-3 pt-4 border-t border-stone-100">
					<div className="flex flex-col gap-1">
						<label className="text-xs font-bold text-gray-500 uppercase">
							Probe Location
						</label>
						<select
							className="p-1.5 border rounded text-xs bg-white w-full border-stone-200 outline-none focus:ring-1 focus:ring-orange-500"
							value={
								selectedEdge.data.probe_mode === "preset" &&
								selectedEdge.data.probe_preset
									? `preset:${selectedEdge.data.probe_preset.type}:${selectedEdge.data.probe_preset.index}`
									: "custom"
							}
							onChange={(e) => {
								const val = e.target.value;
								if (val === "custom") {
									updateEdgeData(
										selectedEdge.id,
										{ probe_mode: "custom" },
										true,
									);
								} else {
									const [_, type, idx] = val.split(":");
									updateEdgeData(
										selectedEdge.id,
										{
											probe_mode: "preset",
											probe_preset: { type, index: Number.parseInt(idx, 10) },
										},
										true,
									);
								}
							}}
						>
							<option value="custom">Custom Position (d=x)</option>
							{/* Predefined locations for Layer 1 */}
							<option value="preset:hot:0">
								L1 {flowsForward ? "Hot" : "Cold"} Side
							</option>
							<option value="preset:avg:0">L1 Average</option>
							<option value="preset:cold:0">
								L1 {flowsForward ? "Cold" : "Hot"} Side
							</option>
							{/* Dynamic locations for additional layers */}
							{(selectedEdge.data.layers || [])
								.slice(1)
								.map((_: any, i: number) => (
									<React.Fragment key={i + 1}>
										<option value={`preset:avg:${i + 1}`}>
											L{i + 2} Average
										</option>
										<option value={`preset:cold:${i + 1}`}>
											L{i + 2} {flowsForward ? "Cold" : "Hot"} Side
										</option>
									</React.Fragment>
								))}
						</select>
					</div>

					{selectedEdge.data.probe_mode === "custom" && (
						<LengthInput
							id={`probe_depth_${selectedEdge.id}`}
							label="Custom Depth d"
							value={selectedEdge.data.probe_depth ?? 0}
							onChange={(val) => {
								const totalThickness =
									(selectedEdge.data.layers || []).reduce(
										(acc: number, l: any) => acc + (l.thickness || 0),
										0,
									) || 0.003;
								const clamped = Math.min(Math.max(0, val), totalThickness);
								updateEdgeData(selectedEdge.id, { probe_depth: clamped }, true);
							}}
							placeholder="0.00"
						/>
					)}
				</div>

				{selectedEdge.data.result && (
					<div className="mt-4 border-t pt-4 flex flex-col gap-3">
						<h3 className="text-xs font-bold text-orange-600 uppercase">
							Solve Results
						</h3>
						<div className="bg-stone-50 rounded p-3 flex flex-col gap-2 border border-stone-100 shadow-sm">
							<div className="flex justify-between items-center">
								<span className="text-[10px] text-stone-500 font-bold uppercase">
									Heat Flow (Q)
								</span>
								<span className="font-mono text-xs font-black text-orange-700">
									{selectedEdge.data.result.Q?.toLocaleString()} W
								</span>
							</div>
							{solvedArrangement(selectedEdge.data.result) && (
								<div className="flex justify-between items-center">
									<span className="text-[10px] text-stone-500 font-bold uppercase">
										Arrangement
									</span>
									<span className="font-mono text-xs font-black">
										{solvedArrangement(selectedEdge.data.result)}
										{selectedEdge.data.result.NTU != null &&
											` · NTU ${Number(selectedEdge.data.result.NTU).toPrecision(3)}`}
										{selectedEdge.data.result.effectiveness != null &&
											` · ε ${Number(selectedEdge.data.result.effectiveness).toFixed(3)}`}
									</span>
								</div>
							)}
							<div className="flex justify-between items-center">
								<span className="text-[10px] text-stone-500 font-bold uppercase">
									Hot Side Wall Temp
								</span>
								<span className="font-mono text-xs font-black">
									{(() => {
										const tInt = selectedEdge.data.result.T_interface as
											| number[]
											| undefined;
										if (tInt && tInt.length >= 2) {
											return flowsForward
												? tInt[0].toFixed(1)
												: tInt[tInt.length - 1].toFixed(1);
										}
										return selectedEdge.data.result.T_hot?.toFixed(1);
									})()} K
								</span>
							</div>
							{selectedEdge.data.result.material_in_range === 0 && (
								<div className="text-[10px] text-amber-700 bg-amber-50 border border-amber-200 rounded px-2 py-1">
									A layer's k is outside its material's sourced range: it
									follows the k(T) envelope (continued from the data,
									unsourced).
								</div>
							)}
							<div className="grid grid-cols-4 gap-1 pt-2 border-t border-stone-100">
								<div className="flex flex-col">
									<span className="text-[7px] text-stone-400 uppercase">
										h (Side A)
									</span>
									<span className="text-[10px] font-mono font-bold">
										{h_a < 1e6 ? h_a.toFixed(1) : "∞"}
									</span>
								</div>
								<div className="flex flex-col">
									<span className="text-[7px] text-stone-400 uppercase">
										h (Wall)
									</span>
									<span className="text-[10px] font-mono font-bold text-orange-600">
										{h_wall < 1e6 ? h_wall.toFixed(1) : "∞"}
									</span>
								</div>
								<div className="flex flex-col">
									<span className="text-[7px] text-stone-400 uppercase">
										h (Side B)
									</span>
									<span className="text-[10px] font-mono font-bold">
										{h_b < 1e6 ? h_b.toFixed(1) : "∞"}
									</span>
								</div>
								<div className="flex flex-col">
									<span className="text-[7px] text-stone-400 uppercase">
										h (Total)
									</span>
									<span className="text-[10px] font-mono font-bold text-blue-600">
										{h_total.toFixed(1)}
									</span>
								</div>
							</div>
							{probeTemp !== null && (
								<div className="flex justify-between items-center pt-2 mt-2 border-t border-stone-200 bg-orange-50/50 p-2 rounded border">
									<span className="text-[10px] text-orange-600 font-bold uppercase tracking-wider">
										{probeLabel}
									</span>
									<span className="font-mono text-sm font-black text-orange-600">
										{probeTemp.toFixed(1)} K
									</span>
								</div>
							)}
						</div>
					</div>
				)}
			</aside>
		);
	}

	if (selectedEdge) {
		return (
			<aside
				key={selectedEdge.id}
				className="w-80 border-l bg-white p-4 flex flex-col gap-4 overflow-y-auto"
			>
				<h2 className="text-lg font-bold border-b pb-2 uppercase text-slate-600">
					Connection
				</h2>

				{/* Custom label */}
				<div className="flex flex-col gap-1">
					<label className="text-xs font-bold text-gray-500 uppercase">
						Custom Label (optional)
					</label>
					<input
						type="text"
						value={selectedEdge.data?.label ?? ""}
						onChange={(e) =>
							updateEdgeData(
								selectedEdge.id,
								{ label: e.target.value || undefined },
								true,
							)
						}
						placeholder="e.g. Main Feed"
						className="p-1.5 border rounded text-xs w-full focus:ring-1 focus:ring-slate-400 outline-none"
					/>
				</div>

				{/* Show label on graph */}
				<label className="flex items-center gap-2 text-xs text-gray-600 cursor-pointer select-none">
					<input
						type="checkbox"
						checked={!!selectedEdge.data?.show_label}
						onChange={(e) =>
							updateEdgeData(
								selectedEdge.id,
								{ show_label: e.target.checked },
								true,
							)
						}
						className="rounded"
					/>
					Show label in graph
				</label>

				{/* Auto ID for reference */}
				<div className="flex flex-col gap-1 mt-2 border-t pt-3">
					<span className="text-[10px] font-bold text-stone-400 uppercase">
						Internal ID
					</span>
					<span className="font-mono text-[10px] text-stone-400 break-all">
						{selectedEdge.id}
					</span>
				</div>
			</aside>
		);
	}

	return null;
};

const InitialGuessEditor = ({ node }: { node: any }) => {
	const { updateNodeData } = useStore();
	const initialGuess = node.data.initial_guess || {};

	const updateGuess = (key: string, val: number | undefined) => {
		const newGuess = { ...initialGuess };
		if (val === undefined || Number.isNaN(val)) {
			delete newGuess[key];
		} else {
			newGuess[key] = val;
		}
		updateNodeData(node.id, { initial_guess: newGuess });
	};

	// Type-aware guess field selection. Elements carry only a solvable m_dot
	// unknown; nodes carry P / Pt. Temperature is never an unknown: it is
	// propagated from the inflows each evaluation, so a T guess has nothing
	// to seed. Pressure-boundary nodes are fully specified (no unknowns) and
	// are handled by returning an empty field list below.
	const FIELDS_BY_TYPE: Record<string, { key: string; unit: string }[]> = {
		// Volume nodes: P and Pt are the unknowns ("P" seeds both).
		plenum: [{ key: "P", unit: "Pa" }],
		combustor: [{ key: "P", unit: "Pa" }],
		momentum_chamber: [{ key: "P", unit: "Pa" }],
		// Boundaries: Pt is solved (mass-flow inlet) or fixed (pressure
		// inlet). Short key "P" is broadcast to both by the backend; we show
		// it only where at least one of the two is actually an unknown.
		mass_boundary: [{ key: "P", unit: "Pa" }],
		fuel_boundary: [{ key: "P", unit: "Pa" }],
		// Elements (edges in React Flow but still nodes here): m_dot only.
		channel: [{ key: "m_dot", unit: "kg/s" }],
		orifice: [{ key: "m_dot", unit: "kg/s" }],
		compressor: [{ key: "m_dot", unit: "kg/s" }],
		turbine: [{ key: "m_dot", unit: "kg/s" }],
		discrete_loss: [{ key: "m_dot", unit: "kg/s" }],
		area_change: [{ key: "m_dot", unit: "kg/s" }],
		// MPCE tee: the only element-owned unknown is the junction pressure,
		// named by what it is (#773): the TOTAL Pt_jct for the Mynard and
		// constant-K models, the static P_jct for the ideal momentum tee
		// (see junctionGuessKey below). Warm-starting it is the most useful
		// nudge for the slower MPCE solve.
		mpce_tee: [
			{ key: junctionGuessKey(node.data?.junction_model), unit: "Pa" },
		],
		// Ejector: element-owned unknown is P_jct (repurposed as the critical
		// back pressure P_c*); warm-starting it helps this stiff solve.
		ejector: [{ key: "P_jct", unit: "Pa" }],
		vortex: [{ key: "m_dot", unit: "kg/s" }],
	};
	const guessFields = FIELDS_BY_TYPE[node.type as string] ?? [];
	if (guessFields.length === 0) {
		return null;
	}

	return (
		<div className="mt-4 border-t pt-4">
			<h3 className="text-xs font-bold text-stone-500 uppercase mb-2">
				Initial Guess Overrides
			</h3>
			<div className="flex flex-col gap-2">
				{guessFields.map(({ key, unit }) => (
					<div key={key} className="flex items-center justify-between gap-2">
						<span className="text-[10px] font-mono text-gray-500">
							{key}{" "}
							<span className="text-[9px] text-gray-400 font-normal">
								[{unit}]
							</span>
						</span>
						<NumericInput
							placeholder="Auto"
							value={initialGuess[key]}
							onChange={(val) => updateGuess(key, val || undefined)}
							className="p-1 border rounded text-xs w-24"
						/>
					</div>
				))}
			</div>
		</div>
	);
};

export default Inspector;
