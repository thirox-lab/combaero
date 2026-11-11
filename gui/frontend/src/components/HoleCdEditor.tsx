import type React from "react";
import {
	autoInletCrossflow,
	correlationLabel,
	exitCrossflowLabel,
	HOLE_CD_RULES,
	inclinationLabel,
	inletCrossflowLabel,
	resolvedCorrelation,
	resolvedExitCrossflow,
	withCorrelation,
} from "../utils/holeCd";
import type {
	HoleCdCorrelation,
	HoleCdData,
	HoleCdExitCrossflow,
	HoleCdInclination,
	HoleCdInletCrossflow,
	HoleCdNodeType,
} from "./flowTypes";
import NumericInput from "./NumericInput";

// A hole node's discharge coefficient (#1098): the node's `hole_cd`, the
// library's HoleCd, edited as one section. Each select offers only what this
// node type builds (utils/holeCd.ts), and its "auto" names what it is here.

interface HoleCdEditorProps {
	id: string;
	type: HoleCdNodeType;
	value: HoleCdData | null | undefined;
	onChange: (next: HoleCdData) => void;
	// The node fields "auto" inclination reads.
	geom?: { hole_type?: string; angle_deg?: number };
	// Shown under the exit-crossflow select while Rowbury et al. is chosen:
	// the term's geometry, which is the node's (region, orientation).
	rowburyFields?: React.ReactNode;
}

const labelCls = "text-xs font-bold text-gray-500 uppercase";
const selectCls = "p-2 border rounded bg-white text-xs border-stone-200";

const HoleCdEditor: React.FC<HoleCdEditorProps> = ({
	id,
	type,
	value,
	onChange,
	geom = {},
	rowburyFields,
}) => {
	const hc: HoleCdData = value ?? {};
	const rules = HOLE_CD_RULES[type];
	const correlation = resolvedCorrelation(type, hc);
	const set = (patch: Partial<HoleCdData>) => onChange({ ...hc, ...patch });
	const showInclination =
		rules.inclination.length > 0 && correlation === "IdelchikThick";
	const inlet = hc.inlet_crossflow ?? "auto";
	// The liner's duct-fed panels need McGreehan-Schotsch's crossflow term;
	// another correlation is refused unless the duct is taken as a plenum.
	const ductRefused =
		rules.inletCrossflow.length > 0 &&
		correlation !== "McGreehanSchotsch" &&
		(inlet === "mcgreehan_schotsch" ||
			(inlet === "auto" && autoInletCrossflow(type) === "mcgreehan_schotsch"));
	return (
		<div className="flex flex-col gap-2 rounded border border-stone-200 p-2">
			<span className="text-xs font-bold uppercase text-gray-500">Hole Cd</span>
			<div className="flex flex-col gap-1">
				<label htmlFor={`${id}_correlation`} className={labelCls}>
					Discharge Model
				</label>
				<select
					id={`${id}_correlation`}
					className={selectCls}
					value={hc.correlation ?? ""}
					onChange={(e) =>
						onChange(
							withCorrelation(
								type,
								hc,
								(e.target.value || null) as HoleCdCorrelation | null,
							),
						)
					}
				>
					<option value="">
						{`Default: ${correlationLabel(type, rules.defaultCorrelation)}`}
					</option>
					{rules.correlations.map((c) => (
						<option key={c} value={c}>
							{correlationLabel(type, c)}
						</option>
					))}
				</select>
			</div>
			{correlation === "fixed" && (
				<div className="flex flex-col gap-1">
					<label htmlFor={`${id}_Cd`} className={labelCls}>
						Fixed Cd Value
					</label>
					<NumericInput
						id={`${id}_Cd`}
						value={hc.Cd ?? rules.fixedCd}
						onChange={(val) => set({ Cd: val })}
						className="p-1.5 h-8 text-sm border rounded bg-white"
						placeholder={String(rules.fixedCd)}
						min={0}
					/>
				</div>
			)}
			{showInclination && (
				<div className="flex flex-col gap-1">
					<label htmlFor={`${id}_inclination`} className={labelCls}>
						Hole Angle Term
					</label>
					<select
						id={`${id}_inclination`}
						className={selectCls}
						value={hc.inclination ?? "auto"}
						onChange={(e) =>
							set({ inclination: e.target.value as HoleCdInclination })
						}
					>
						{rules.inclination.map((v) => (
							<option key={v} value={v}>
								{inclinationLabel(v, type, hc, geom)}
							</option>
						))}
					</select>
				</div>
			)}
			{rules.inletCrossflow.length > 0 && (
				<div className="flex flex-col gap-1">
					<label htmlFor={`${id}_inlet_crossflow`} className={labelCls}>
						Supply Crossflow
					</label>
					<select
						id={`${id}_inlet_crossflow`}
						className={selectCls}
						value={inlet}
						onChange={(e) =>
							set({
								inlet_crossflow: e.target.value as HoleCdInletCrossflow,
							})
						}
					>
						{rules.inletCrossflow.map((v) => (
							<option
								key={v}
								value={v}
								disabled={
									v === "mcgreehan_schotsch" &&
									correlation !== "McGreehanSchotsch"
								}
							>
								{inletCrossflowLabel(v, type)}
							</option>
						))}
					</select>
					{ductRefused && (
						<p className="text-[9px] text-amber-700">
							The duct crossflow term is McGreehan-Schotsch's own: with this
							correlation the solver refuses the duct-fed panels. Pick None to
							take the duct as a plenum.
						</p>
					)}
				</div>
			)}
			{rules.exitCrossflow.length > 0 && (
				<div className="flex flex-col gap-1">
					<label htmlFor={`${id}_exit_crossflow`} className={labelCls}>
						{type === "impingement_array"
							? "Crossflow on the Holes"
							: "Exit Crossflow"}
					</label>
					<select
						id={`${id}_exit_crossflow`}
						className={selectCls}
						value={hc.exit_crossflow ?? "auto"}
						onChange={(e) =>
							set({ exit_crossflow: e.target.value as HoleCdExitCrossflow })
						}
					>
						{rules.exitCrossflow.map((v) => (
							<option key={v} value={v}>
								{exitCrossflowLabel(v, type)}
							</option>
						))}
					</select>
					<div className="text-[9px] text-gray-400">
						{type === "impingement_array"
							? "A strong crossflow bends the jets over and lowers the Cd. Below Gc/Gj ~1.5 the hole's own Cd holds."
							: "Rowbury et al.: the mainstream of a flowing discharge chamber lowers the Cd; the solver refuses it into a plenum."}
					</div>
					{resolvedExitCrossflow(type, hc) === "rowbury_2001" && rowburyFields}
				</div>
			)}
		</div>
	);
};

export default HoleCdEditor;
