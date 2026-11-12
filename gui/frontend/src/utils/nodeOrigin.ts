import type { Node, XYPosition } from "reactflow";
import { SPLITTER_WIDTH, splitterFractions, splitterHeight } from "./splitter";

/**
 * A node's `position` is its card CENTRE (#865): cards of different heights
 * placed at one y have their left/right flow handles level, so a main path
 * drawn on one row runs straight. React Flow subtracts the origin when it
 * draws; `positionAbsolute` stays the centre, so geometry built from it must
 * go through `nodeTopLeft`.
 */
export const NODE_ORIGIN: [number, number] = [0.5, 0.5];

/** Canvas snap grid [px]; centres snap to it, so cards dragged near one row
 *  land on it. */
export const SNAP_GRID: [number, number] = [10, 10];

/** Saved-file version. 1.2: `position` is the card centre; up to 1.1 it was
 *  the top-left corner. 1.3: a hole node's Cd choices are one `hole_cd`
 *  object; an older file's loose fields are lifted on load, by key, not by
 *  version (legacyKeys.ts). */
export const NETWORK_FILE_VERSION = "1.3";

/** The first version whose `position` is the card centre (#865). */
export const CENTRE_POSITION_VERSION = "1.2";

// Card sizes [px] as rendered (measured headlessly, #865), for files that
// carry no measured size. The lossless connection and the probe size to their
// content; theirs are the sizes with a short label.
const CARD_SIZE: Record<string, [number, number]> = {
	lossless_connection: [120, 36],
	probe: [140, 96],
	plenum: [160, 72],
	mass_boundary: [144, 56],
	fuel_boundary: [150, 64],
	pressure_boundary: [128, 56],
	channel: [140, 56],
	orifice: [140, 48],
	compressor: [140, 56],
	turbine: [140, 56],
	impingement_array: [172, 48],
	effusion_plate: [140, 48],
	film_cooling: [140, 48],
	slot_film: [140, 48],
	effusion_liner: [150, 56],
	combustor: [160, 72],
	momentum_chamber: [160, 72],
	discrete_loss: [144, 58],
	area_change: [140, 48],
	mpce_tee: [140, 56],
	ejector: [140, 56],
	vortex: [140, 56],
	wall: [128, 56],
	heat_source: [150, 60],
};
const FALLBACK_SIZE: [number, number] = [140, 48];

export function cardSize(n: Node): [number, number] {
	if (n.type === "flow_splitter") {
		// It grows with its outlets (#903).
		const rows = splitterFractions(n.data ?? {}).length;
		return [n.width ?? SPLITTER_WIDTH, n.height ?? splitterHeight(rows)];
	}
	const [w, h] = CARD_SIZE[n.type ?? ""] ?? FALLBACK_SIZE;
	return [n.width ?? w, n.height ?? h];
}

/** The card's top-left corner in flow coordinates. */
export function nodeTopLeft(n: Node): XYPosition {
	const p = n.positionAbsolute ?? n.position;
	const [w, h] = cardSize(n);
	return { x: p.x - w * NODE_ORIGIN[0], y: p.y - h * NODE_ORIGIN[1] };
}

function versionBefore(version: unknown, than: string): boolean {
	const parse = (v: string) =>
		v.split(".").map((s) => Number.parseInt(s, 10) || 0);
	const [a, b] = [
		parse(typeof version === "string" ? version : "0"),
		parse(than),
	];
	for (let i = 0; i < Math.max(a.length, b.length); i++) {
		if ((a[i] ?? 0) !== (b[i] ?? 0)) return (a[i] ?? 0) < (b[i] ?? 0);
	}
	return false;
}

/** Positions of a file saved before 1.2 (top-left corners) as card centres,
 *  so the network lands where it was drawn. Newer files pass through. */
export function migrateNodePositions(nodes: Node[], version: unknown): Node[] {
	if (!versionBefore(version, CENTRE_POSITION_VERSION)) return nodes;
	return nodes.map((n) => {
		const [w, h] = cardSize(n);
		const at = (p: XYPosition | undefined) =>
			p && { x: p.x + w * NODE_ORIGIN[0], y: p.y + h * NODE_ORIGIN[1] };
		return {
			...n,
			position: at(n.position) as XYPosition,
			...(n.positionAbsolute
				? { positionAbsolute: at(n.positionAbsolute) }
				: {}),
		};
	});
}
