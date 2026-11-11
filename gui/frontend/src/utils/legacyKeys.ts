// Renamed input keys in saved networks (#773). A file written before a
// rename still loads: the old key moves to the new one unless the new one is
// already set. The backend accepts the old names too (pydantic aliases and
// the junction guess migration in gui/backend), so a headless payload needs
// no migration.

import { liftLooseHoleCd } from "./holeCd";

/** omega_rpm became shaft_speed_rpm: omega is 2 pi f, not a speed in rpm. */
const RENAMED: Record<string, string> = { omega_rpm: "shaft_speed_rpm" };

const renameKey = (obj: any, oldKey: string, newKey: string): any => {
	if (!obj || typeof obj !== "object" || !(oldKey in obj)) return obj;
	const { [oldKey]: value, ...rest } = obj;
	return newKey in rest && rest[newKey] != null
		? rest
		: { ...rest, [newKey]: value };
};

const renameKeys = (obj: any): any => {
	let out = obj;
	for (const [oldKey, newKey] of Object.entries(RENAMED)) {
		out = renameKey(out, oldKey, newKey);
	}
	return out;
};

/** The tee's owned unknown, named by what it holds (#773): the junction
 *  TOTAL pressure for the Mynard and constant-K models, the static for the
 *  ideal momentum tee. */
export const junctionGuessKey = (junctionModel: unknown): "P_jct" | "Pt_jct" =>
	junctionModel === "ideal_momentum" ? "P_jct" : "Pt_jct";

/** Node data under its current key names: the vortex's shaft speed, a
 *  Mynard / constant-K tee's junction guess (P_jct before #773), and a hole
 *  node's loose Cd fields in one `hole_cd` (file version 1.3). */
export const migrateNodeData = (data: any, type?: string): any => {
	let out = liftLooseHoleCd(renameKeys(data), type);
	if (
		type === "mpce_tee" &&
		out?.initial_guess &&
		junctionGuessKey(out.junction_model) === "Pt_jct"
	) {
		out = {
			...out,
			initial_guess: renameKey(out.initial_guess, "P_jct", "Pt_jct"),
		};
	}
	return out;
};

/** Solver settings (the global shaft speed) under their current key names. */
export const migrateSolverSettings = (settings: any): any =>
	renameKeys(settings);
