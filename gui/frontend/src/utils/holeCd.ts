import type {
	HoleCdCorrelation,
	HoleCdData,
	HoleCdExitCrossflow,
	HoleCdInclination,
	HoleCdInletCrossflow,
	HoleCdNodeType,
} from "../components/flowTypes";

// Which hole-Cd choices each node type takes, and what its null / "auto"
// resolves to: the library's rules (python/combaero/network/hole_cd.py,
// HoleCd's docstring), kept here so the Inspector offers only what the
// element builds and names what "auto" is for this node.

interface HoleCdRules {
	// Offered in this order; null resolves to `defaultCorrelation`.
	correlations: HoleCdCorrelation[];
	defaultCorrelation: HoleCdCorrelation;
	// The Cd a switch to 'fixed' fills in (the element's default where it
	// has one: Florschuetz's 0.79 on an impingement array).
	fixedCd: number;
	// The choices of each term on this node type; a type with no such term
	// takes only "auto"/"none" and the Inspector does not show the select.
	inclination: HoleCdInclination[];
	inletCrossflow: HoleCdInletCrossflow[];
	exitCrossflow: HoleCdExitCrossflow[];
}

const DISCHARGE: HoleCdCorrelation[] = [
	"IdelchikThick",
	"Lichtarowicz",
	"McGreehanSchotsch",
	"fixed",
];
const INCLINED: HoleCdInclination[] = ["auto", "none", "idelchik_3_2"];
const ROWBURY: HoleCdExitCrossflow[] = ["auto", "none", "rowbury_2001"];

export const HOLE_CD_RULES: Record<HoleCdNodeType, HoleCdRules> = {
	orifice: {
		correlations: [
			"ReaderHarrisGallagher",
			"IdelchikThick",
			"IdelchikBeveled",
			"IdelchikRounded",
			"McGreehanSchotsch",
			"Lichtarowicz",
			"fixed",
		],
		defaultCorrelation: "ReaderHarrisGallagher",
		fixedCd: 0.6,
		inclination: [],
		inletCrossflow: [],
		exitCrossflow: [],
	},
	effusion_plate: {
		correlations: DISCHARGE,
		defaultCorrelation: "IdelchikThick",
		fixedCd: 0.6,
		inclination: INCLINED,
		inletCrossflow: [],
		exitCrossflow: ROWBURY,
	},
	film_cooling: {
		correlations: DISCHARGE,
		defaultCorrelation: "IdelchikThick",
		fixedCd: 0.6,
		inclination: INCLINED,
		inletCrossflow: [],
		exitCrossflow: ROWBURY,
	},
	slot_film: {
		correlations: DISCHARGE,
		defaultCorrelation: "IdelchikThick",
		fixedCd: 0.6,
		inclination: [],
		inletCrossflow: [],
		exitCrossflow: [],
	},
	impingement_array: {
		correlations: [
			"fixed",
			"IdelchikThick",
			"Lichtarowicz",
			"McGreehanSchotsch",
		],
		defaultCorrelation: "fixed",
		fixedCd: 0.79,
		inclination: [],
		inletCrossflow: [],
		exitCrossflow: [
			"auto",
			"rogers_hersh_1976",
			"florschuetz_metzger_su_1982",
			"none",
		],
	},
	effusion_liner: {
		correlations: [
			"McGreehanSchotsch",
			"IdelchikThick",
			"Lichtarowicz",
			"fixed",
		],
		defaultCorrelation: "McGreehanSchotsch",
		fixedCd: 0.6,
		inclination: INCLINED,
		inletCrossflow: ["auto", "none", "mcgreehan_schotsch"],
		exitCrossflow: ROWBURY,
	},
};

export const isHoleCdNode = (type?: string): type is HoleCdNodeType =>
	type != null && type in HOLE_CD_RULES;

export const CORRELATION_LABELS: Record<HoleCdCorrelation, string> = {
	ReaderHarrisGallagher: "Reader-Harris/Gallagher (Sharp)",
	IdelchikThick: "Deep Hole in Wall (Idelchik 4-18a)",
	IdelchikBeveled: "Beveled Hole in Wall (Idelchik 4-18b)",
	IdelchikRounded: "Rounded Hole in Wall (Idelchik 4-18c)",
	McGreehanSchotsch: "Cooling Hole (M-S 1988, plenum-fed)",
	Lichtarowicz: "Long Orifice, l/d 2-10 (Lichtarowicz 1965)",
	fixed: "Manual / Fixed Value",
};

// A node type's own wording where the hole differs: a slot's correlations
// see its hydraulic diameter, the liner's panels are duct-fed, and an
// impingement array's fixed value is Florschuetz's.
const CORRELATION_LABEL_OVERRIDES: Partial<
	Record<HoleCdNodeType, Partial<Record<HoleCdCorrelation, string>>>
> = {
	slot_film: {
		IdelchikThick: "Deep Hole in Wall (Idelchik 4-18a, on Dh)",
		Lichtarowicz: "Long Orifice (Lichtarowicz 1965, on Dh)",
		McGreehanSchotsch: "Cooling Hole (M-S 1988, on Dh)",
	},
	effusion_liner: {
		McGreehanSchotsch: "Cooling Hole (M-S 1988, duct-fed)",
	},
	impingement_array: { fixed: "Fixed value (Florschuetz default 0.79)" },
};

export const correlationLabel = (
	type: HoleCdNodeType,
	c: HoleCdCorrelation,
): string => CORRELATION_LABEL_OVERRIDES[type]?.[c] ?? CORRELATION_LABELS[c];

export const resolvedCorrelation = (
	type: HoleCdNodeType,
	hc: HoleCdData | null | undefined,
): HoleCdCorrelation =>
	hc?.correlation ?? HOLE_CD_RULES[type].defaultCorrelation;

// The node fields "auto" inclination reads besides the spec.
interface HoleGeometry {
	hole_type?: string;
	angle_deg?: number;
}

// The node type's own defaults of those fields (gui/backend/schemas.py).
const DEFAULT_ANGLE: Partial<Record<HoleCdNodeType, number>> = {
	film_cooling: 30,
	effusion_plate: 90,
	effusion_liner: 90,
};

/** What "auto" inclination is on this node: Idelchik 3-2 for an inclined
 *  cylindrical hole on IdelchikThick, none otherwise (_resolve_cd_inclination). */
export const autoInclination = (
	type: HoleCdNodeType,
	hc: HoleCdData | null | undefined,
	geom: HoleGeometry,
): "none" | "idelchik_3_2" => {
	if (HOLE_CD_RULES[type].inclination.length === 0) return "none";
	const holeType =
		type === "film_cooling" ? (geom.hole_type ?? "shaped") : "cylindrical";
	const angle = geom.angle_deg ?? DEFAULT_ANGLE[type] ?? 90;
	return resolvedCorrelation(type, hc) === "IdelchikThick" &&
		holeType === "cylindrical" &&
		angle < 90
		? "idelchik_3_2"
		: "none";
};

/** Why "auto" inclination is off, for its label. */
const inclinationOffReason = (
	type: HoleCdNodeType,
	hc: HoleCdData | null | undefined,
	geom: HoleGeometry,
): string => {
	if (resolvedCorrelation(type, hc) !== "IdelchikThick")
		return "IdelchikThick only";
	if (type === "film_cooling" && (geom.hole_type ?? "shaped") === "shaped")
		return "shaped hole";
	return "normal hole";
};

export const inclinationLabel = (
	v: HoleCdInclination,
	type: HoleCdNodeType,
	hc: HoleCdData | null | undefined,
	geom: HoleGeometry,
): string => {
	if (v === "auto")
		return autoInclination(type, hc, geom) === "idelchik_3_2"
			? "Auto: inclined inlet (Idelchik 3-2)"
			: `Auto: none (${inclinationOffReason(type, hc, geom)})`;
	return v === "none"
		? "None (angle via hole length only)"
		: "Inclined inlet (Idelchik 3-2, forced)";
};

/** What "auto" supply crossflow is: the backside duct's U1/Vi on the liner's
 *  duct-fed panels, none (plenum-fed) everywhere else. */
export const autoInletCrossflow = (
	type: HoleCdNodeType,
): "none" | "mcgreehan_schotsch" =>
	type === "effusion_liner" ? "mcgreehan_schotsch" : "none";

export const inletCrossflowLabel = (
	v: HoleCdInletCrossflow,
	type: HoleCdNodeType,
): string => {
	if (v === "auto")
		return autoInletCrossflow(type) === "mcgreehan_schotsch"
			? "Auto: duct crossflow U1/Vi (M-S 1988)"
			: "Auto: none (plenum-fed)";
	return v === "none"
		? "None (duct taken as a plenum)"
		: "Duct crossflow U1/Vi (M-S 1988)";
};

/** What "auto" exit crossflow is: the gap's crossflow law on an impingement
 *  array (#800), none on every other hole (Rowbury is opt-in, #1074). */
export const autoExitCrossflow = (
	type: HoleCdNodeType,
): "none" | "rogers_hersh_1976" =>
	type === "impingement_array" ? "rogers_hersh_1976" : "none";

const EXIT_LABELS: Record<Exclude<HoleCdExitCrossflow, "auto">, string> = {
	none: "None (constant Cd)",
	rowbury_2001: "Mainstream past the exit (Rowbury et al. 2001)",
	rogers_hersh_1976: "Grazing flow (Rogers & Hersh 1976)",
	florschuetz_metzger_su_1982: "CR-3630 Table 4.1 (tested plates only)",
};

export const exitCrossflowLabel = (
	v: HoleCdExitCrossflow,
	type: HoleCdNodeType,
): string => {
	if (v !== "auto") return EXIT_LABELS[v];
	return autoExitCrossflow(type) === "rogers_hersh_1976"
		? "Auto: grazing flow (Rogers & Hersh 1976)"
		: "Auto: none";
};

export const resolvedExitCrossflow = (
	type: HoleCdNodeType,
	hc: HoleCdData | null | undefined,
): Exclude<HoleCdExitCrossflow, "auto"> => {
	const v = hc?.exit_crossflow ?? "auto";
	return v === "auto" ? autoExitCrossflow(type) : v;
};

/** The Cd a 'fixed' hole_cd WITHOUT one builds with: the element's own
 *  default where it has one (Florschuetz's 0.79 on an impingement array),
 *  null where the solver refuses it (HoleCd: None = the element default). */
export const fixedCdDefault = (type: HoleCdNodeType): number | null =>
	type === "impingement_array" ? HOLE_CD_RULES[type].fixedCd : null;

/** A correlation change (null = the node type's default). The library
 *  refuses a forced term that does not compose with the new correlation, so
 *  leaving IdelchikThick resets a forced "idelchik_3_2" and leaving
 *  McGreehan-Schotsch a forced "mcgreehan_schotsch" to "auto" (#1094): no
 *  hidden choice makes the build fail. The Cd is the 'fixed' correlation's
 *  only (HoleCd refuses it with any other): switching to 'fixed' keeps the
 *  node's Cd, else restores `remembered` (the user's last fixed Cd on this
 *  node, held by the editor, not in the file), else fills the type's
 *  default; switching away drops it. */
export const withCorrelation = (
	type: HoleCdNodeType,
	hc: HoleCdData | null | undefined,
	correlation: HoleCdCorrelation | null,
	remembered?: number | null,
): HoleCdData => {
	const { Cd, ...rest } = hc ?? {};
	const next: HoleCdData = { ...rest, correlation };
	const resolved = correlation ?? HOLE_CD_RULES[type].defaultCorrelation;
	if (resolved !== "IdelchikThick" && next.inclination === "idelchik_3_2")
		next.inclination = "auto";
	if (
		resolved !== "McGreehanSchotsch" &&
		next.inlet_crossflow === "mcgreehan_schotsch"
	)
		next.inlet_crossflow = "auto";
	if (resolved === "fixed")
		next.Cd = Cd ?? remembered ?? HOLE_CD_RULES[type].fixedCd;
	return next;
};

// The node types whose v0.11.1 files carry the Cd choices loose on the node,
// with the defaults that build read: (correlation, Cd). Only the impingement
// array had `cd_crossflow`; on any other type it was ignored, and still is.
const LEGACY_DEFAULTS: Partial<
	Record<HoleCdNodeType, [HoleCdCorrelation, number]>
> = {
	orifice: ["ReaderHarrisGallagher", 0.6],
	effusion_plate: ["IdelchikThick", 0.6],
	film_cooling: ["IdelchikThick", 0.6],
	slot_film: ["IdelchikThick", 0.6],
	impingement_array: ["fixed", 0.79],
};
const legacyKeys = (type: HoleCdNodeType): string[] =>
	type === "impingement_array"
		? ["correlation", "Cd", "cd_crossflow"]
		: ["correlation", "Cd"];

// Node keys of the unreleased dev builds (#1074 to #1098): not lifted, not
// ignored -- the backend refuses them, and the load notice says where each
// goes (an ignored "none" would silently switch the term back to "auto").
export const PRE_RELEASE_HOLE_CD_KEYS: Record<string, string> = {
	cd_inclination: "hole_cd.inclination",
	cd_exit_crossflow: "hole_cd.exit_crossflow",
};

/** "node.key -> hole_cd.field" for every pre-release hole-Cd key in a file. */
export const preReleaseHoleCdKeys = (nodes: any[]): string[] =>
	(nodes ?? []).flatMap((n: any) =>
		isHoleCdNode(n?.type) && n.data && typeof n.data === "object"
			? Object.keys(PRE_RELEASE_HOLE_CD_KEYS)
					.filter((k) => k in n.data)
					.map((k) => `${n.id}.${k} -> ${PRE_RELEASE_HOLE_CD_KEYS[k]}`)
			: [],
	);

/** A file saved before `hole_cd` (up to v0.11.1): its loose `correlation`,
 *  `Cd` and (impingement array) `cd_crossflow` lifted into `hole_cd` as the
 *  backend's HoleNodeData does, so the Inspector shows what is built. A node
 *  with both is left alone: the backend refuses it, naming the fields. */
export const liftLooseHoleCd = (data: any, type?: string): any => {
	if (!data || typeof data !== "object" || !isHoleCdNode(type)) return data;
	const legacy = LEGACY_DEFAULTS[type];
	if (!legacy) return data;
	const keys = legacyKeys(type);
	if (!keys.some((k) => k in data)) return data;
	if (data.hole_cd != null) return data;
	const rest: any = { ...data };
	for (const k of keys) delete rest[k];
	const corr: HoleCdCorrelation = data.correlation || legacy[0];
	const hole_cd: HoleCdData = { correlation: corr };
	if (corr === "fixed") hole_cd.Cd = data.Cd ?? legacy[1];
	if (keys.includes("cd_crossflow") && data.cd_crossflow != null)
		hole_cd.exit_crossflow = data.cd_crossflow;
	return { ...rest, hole_cd };
};
