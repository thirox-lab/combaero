import type { Node } from "reactflow";
import { describe, expect, it } from "vitest";
import {
	CENTRE_POSITION_VERSION,
	cardSize,
	migrateNodePositions,
	NETWORK_FILE_VERSION,
	nodeTopLeft,
} from "./nodeOrigin";

const node = (type: string, x: number, y: number, extra = {}): Node =>
	({ id: type, type, position: { x, y }, data: {}, ...extra }) as Node;

describe("node position = card centre (#865)", () => {
	it("writes files as 1.3 (hole_cd); positions are centres since 1.2", () => {
		expect(NETWORK_FILE_VERSION).toBe("1.3");
		expect(CENTRE_POSITION_VERSION).toBe("1.2");
	});

	it("takes the measured size over the per-kind table", () => {
		expect(cardSize(node("orifice", 0, 0))).toEqual([140, 48]);
		expect(cardSize(node("orifice", 0, 0, { width: 150, height: 50 }))).toEqual(
			[150, 50],
		);
		expect(cardSize(node("no_such_kind", 0, 0))).toEqual([140, 48]);
	});

	it("puts the top-left corner half a card from the centre", () => {
		expect(nodeTopLeft(node("channel", 270, 0))).toEqual({ x: 200, y: -28 });
		const placed = node("channel", 0, 0, {
			positionAbsolute: { x: 70, y: 28 },
		});
		expect(nodeTopLeft(placed)).toEqual({ x: 0, y: 0 });
	});

	it("moves a 1.1 file's corners to centres, so it lands where it was drawn", () => {
		const [ch, ori] = migrateNodePositions(
			[
				node("channel", 200, 0),
				node("orifice", 400, 0, { width: 140, height: 48 }),
			],
			"1.1",
		);
		expect(ch.position).toEqual({ x: 270, y: 28 });
		expect(ori.position).toEqual({ x: 470, y: 24 });
		expect(nodeTopLeft(ch)).toEqual({ x: 200, y: 0 });
		expect(nodeTopLeft(ori)).toEqual({ x: 400, y: 0 });
	});

	it("treats a file with no version as old, and leaves 1.2 and later alone", () => {
		const n = [node("channel", 200, 0)];
		expect(migrateNodePositions(n, undefined)[0].position).toEqual({
			x: 270,
			y: 28,
		});
		expect(migrateNodePositions(n, "1.2")).toBe(n);
		expect(migrateNodePositions(n, "1.3")).toBe(n);
		expect(migrateNodePositions(n, "1.10")).toBe(n);
	});
});
