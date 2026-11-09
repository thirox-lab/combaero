import type React from "react";
import LengthInput from "./LengthInput";
import NumericInput from "./NumericInput";
import UnitInput from "./UnitInput";

// Radiation on one face (#1071): the backend's GasRadiationData (a plate's or
// liner's hot face) or WallRadiationData (a thermal wall side), which build
// the library's GasRadiation / CasingRadiation. Kind "none" builds nothing.
export type RadiationKind = "none" | "gas" | "casing";

export interface RadiationData {
	kind?: RadiationKind;
	beam_length?: number;
	emissivity_wall?: number;
	luminosity?: number;
	T_gas?: number | null;
	fuel_air_ratio?: number | null;
	pressure?: number | null;
	area_ratio?: number;
	emissivity_casing?: number;
	T_casing?: number | null;
}

// The backend's defaults (gui/backend/schemas.py): the library's, with a
// working-example beam length (a 50 mm gas layer) and Lefebvre's
// tubo-annular A_w/A_c of 0.8.
export const RADIATION_DEFAULTS = {
	beam_length: 0.05,
	emissivity_wall: 0.7,
	luminosity: 1.0,
	area_ratio: 0.8,
	emissivity_casing: 0.4,
} as const;

interface RadiationEditorProps {
	id: string;
	label: string;
	value: RadiationData | null | undefined;
	onChange: (next: RadiationData) => void;
	// A plate's hot face sees the gas only; a wall side may see a casing.
	allowCasing?: boolean;
}

const labelCls = "text-xs font-bold text-gray-500 uppercase";
const inputCls = "p-1.5 h-8 text-sm border rounded bg-white";

const RadiationEditor: React.FC<RadiationEditorProps> = ({
	id,
	label,
	value,
	onChange,
	allowCasing = false,
}) => {
	const v: RadiationData = value ?? {};
	const kind: RadiationKind = v.kind ?? "none";
	const set = (patch: Partial<RadiationData>) => onChange({ ...v, ...patch });
	const num = (
		key: keyof RadiationData,
		text: string,
		fallback: number,
		min?: number,
		max?: number,
	) => (
		<div className="flex flex-col gap-1">
			<label htmlFor={`${id}_${key}`} className={labelCls}>
				{text}
			</label>
			<NumericInput
				id={`${id}_${key}`}
				value={(v[key] as number | undefined) ?? fallback}
				onChange={(val) => set({ [key]: val })}
				className={inputCls}
				placeholder={String(fallback)}
				min={min}
				max={max}
			/>
		</div>
	);
	// An optional override: blank is the face's own stream value.
	const optional = (key: keyof RadiationData, text: string) => (
		<div className="flex flex-col gap-1">
			<label htmlFor={`${id}_${key}`} className={labelCls}>
				{text}
			</label>
			<NumericInput
				id={`${id}_${key}`}
				value={(v[key] as number | null | undefined) ?? null}
				onChange={(val) => set({ [key]: val })}
				onClear={() => set({ [key]: null })}
				className={inputCls}
				placeholder="from the stream"
				min={0}
			/>
		</div>
	);
	return (
		<div className="flex flex-col gap-2">
			<label htmlFor={`${id}_kind`} className={labelCls}>
				{label}
			</label>
			<select
				id={`${id}_kind`}
				aria-label={label}
				className="p-2 border rounded bg-white text-xs border-stone-200"
				value={kind}
				onChange={(e) => {
					const next = e.target.value as RadiationKind;
					// The fields stay whatever the kind; a first choice fills
					// the defaults so the node carries what it builds.
					onChange(
						next === "none"
							? { ...v, kind: next }
							: { ...RADIATION_DEFAULTS, ...v, kind: next },
					);
				}}
			>
				<option value="none">None</option>
				<option value="gas">Gas radiation (Lefebvre Eq. 8.7)</option>
				{allowCasing && (
					<option value="casing">Wall to casing (Lefebvre Eq. 8.17)</option>
				)}
			</select>
			{kind === "gas" && (
				<>
					<LengthInput
						id={`${id}_beam_length`}
						label="Mean Beam Length (3.4 V/A)"
						value={v.beam_length ?? RADIATION_DEFAULTS.beam_length}
						onChange={(val) => set({ beam_length: val })}
					/>
					{num(
						"emissivity_wall",
						"Wall Emissivity [-]",
						RADIATION_DEFAULTS.emissivity_wall,
						0,
						1,
					)}
					{num(
						"luminosity",
						"Luminosity Factor [-] (1: non-luminous)",
						RADIATION_DEFAULTS.luminosity,
						0,
					)}
					{optional("T_gas", "Gas Temperature Override [K]")}
					{optional("fuel_air_ratio", "Fuel-Air Ratio Override [-]")}
					<UnitInput
						id={`${id}_pressure`}
						label="Gas Pressure Override"
						value={v.pressure ?? null}
						onChange={(val) => set({ pressure: val })}
						onClear={() => set({ pressure: null })}
					/>
					<p className="text-[10px] text-gray-500 leading-snug">
						Blank overrides take the face's stream: its Tt, its fuel-air ratio
						by elements and its static pressure.
					</p>
				</>
			)}
			{kind === "casing" && (
				<>
					{num(
						"area_ratio",
						"Wall / Casing Area A_w/A_c [-]",
						RADIATION_DEFAULTS.area_ratio,
						0,
					)}
					{num(
						"emissivity_wall",
						"Wall Emissivity [-]",
						RADIATION_DEFAULTS.emissivity_wall,
						0,
						1,
					)}
					{num(
						"emissivity_casing",
						"Casing Emissivity [-]",
						RADIATION_DEFAULTS.emissivity_casing,
						0,
						1,
					)}
					{optional("T_casing", "Casing Temperature Override [K]")}
					<p className="text-[10px] text-gray-500 leading-snug">
						The heat goes into this side's stream; a blank casing temperature is
						the stream's Tt.
					</p>
				</>
			)}
		</div>
	);
};

export default RadiationEditor;
