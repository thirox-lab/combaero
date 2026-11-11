/** @vitest-environment jsdom */
import type { Edge, Node } from "reactflow";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { create } from "zustand";
import {
	AUTOSAVE_KEY,
	type AutosaveSnapshot,
	clearAutosave,
	readAutosave,
	startAutosave,
	writeAutosave,
} from "./autosave";

// An in-memory Storage, so the tests do not depend on what jsdom offers.
const memoryStorage = (): Storage => {
	const m = new Map<string, string>();
	return {
		get length() {
			return m.size;
		},
		clear: () => m.clear(),
		getItem: (k: string) => (m.has(k) ? (m.get(k) as string) : null),
		key: (i: number) => Array.from(m.keys())[i] ?? null,
		removeItem: (k: string) => {
			m.delete(k);
		},
		setItem: (k: string, v: string) => {
			m.set(k, String(v));
		},
	};
};

const node = (id: string, result?: object): Node => ({
	id,
	type: "plenum",
	position: { x: 0, y: 0 },
	data: result ? { id, result } : { id },
});

const snap = (nodes: Node[], edges: Edge[] = []): AutosaveSnapshot => ({
	nodes,
	edges,
	solverSettings: { global_regime: "compressible" },
});

describe("autosave storage (#596)", () => {
	let storage: Storage;
	const get = () => storage;

	beforeEach(() => {
		storage = memoryStorage();
	});

	it("round-trips the network without its solve results", () => {
		const t = new Date("2026-10-18T14:32:00Z");
		expect(writeAutosave(snap([node("a", { P: 1e5 })]), get, t)).toBe(true);
		const rec = readAutosave(get);
		expect(rec).not.toBeNull();
		expect(rec?.nodes[0].id).toBe("a");
		expect(rec?.nodes[0].data.result).toBeUndefined();
		expect(rec?.solverSettings).toEqual({ global_regime: "compressible" });
		expect(rec?.version).toBe("1.3");
		expect(rec?.savedAt).toBe(t.toISOString());
	});

	it("an empty network removes the record, so a clear stays cleared", () => {
		writeAutosave(snap([node("a")]), get);
		writeAutosave(snap([]), get);
		expect(storage.getItem(AUTOSAVE_KEY)).toBeNull();
		expect(readAutosave(get)).toBeNull();
	});

	it("an unreadable record is ignored, not thrown", () => {
		storage.setItem(AUTOSAVE_KEY, "{not json");
		expect(readAutosave(get)).toBeNull();
		storage.setItem(AUTOSAVE_KEY, JSON.stringify({ nodes: "x", edges: [] }));
		expect(readAutosave(get)).toBeNull();
	});

	it("blocked storage (a SecurityError on access) never breaks the editor", () => {
		const blocked = () => null;
		expect(writeAutosave(snap([node("a")]), blocked)).toBe(false);
		expect(readAutosave(blocked)).toBeNull();
		expect(() => clearAutosave(blocked)).not.toThrow();
	});

	it("a full storage (QuotaExceededError) returns false", () => {
		const full = memoryStorage();
		full.setItem = () => {
			throw new DOMException("full", "QuotaExceededError");
		};
		expect(writeAutosave(snap([node("a")]), () => full)).toBe(false);
	});

	it("clearAutosave deletes the record", () => {
		writeAutosave(snap([node("a")]), get);
		clearAutosave(get);
		expect(readAutosave(get)).toBeNull();
	});
});

describe("startAutosave", () => {
	let storage: Storage;
	const get = () => storage;
	const makeStore = () =>
		create<AutosaveSnapshot & { other: number }>(() => ({
			nodes: [],
			edges: [],
			solverSettings: {},
			other: 0,
		}));

	beforeEach(() => {
		storage = memoryStorage();
		vi.useFakeTimers();
	});
	afterEach(() => {
		vi.useRealTimers();
	});

	it("writes once, a debounce after the last edit", () => {
		const store = makeStore();
		const setItem = vi.spyOn(storage, "setItem");
		const stop = startAutosave(store, get, 1000);
		store.setState({ nodes: [node("a")] });
		vi.advanceTimersByTime(500);
		store.setState({ nodes: [node("a"), node("b")] });
		vi.advanceTimersByTime(999);
		expect(setItem).not.toHaveBeenCalled();
		vi.advanceTimersByTime(1);
		expect(setItem).toHaveBeenCalledTimes(1);
		expect(readAutosave(get)?.nodes).toHaveLength(2);
		stop();
	});

	it("ignores state the record does not hold", () => {
		const store = makeStore();
		const setItem = vi.spyOn(storage, "setItem");
		const stop = startAutosave(store, get, 1000);
		store.setState({ other: 1 });
		vi.advanceTimersByTime(5000);
		expect(setItem).not.toHaveBeenCalled();
		stop();
	});

	it("flushes at once on pagehide (Safari and mobile skip beforeunload)", () => {
		const store = makeStore();
		const stop = startAutosave(store, get, 1000);
		store.setState({ nodes: [node("a")] });
		window.dispatchEvent(new Event("pagehide"));
		expect(readAutosave(get)?.nodes).toHaveLength(1);
		stop();
	});

	it("flushes at once when the tab is hidden", () => {
		const store = makeStore();
		const stop = startAutosave(store, get, 1000);
		store.setState({ nodes: [node("a")] });
		Object.defineProperty(document, "visibilityState", {
			configurable: true,
			get: () => "hidden",
		});
		document.dispatchEvent(new Event("visibilitychange"));
		expect(readAutosave(get)?.nodes).toHaveLength(1);
		Object.defineProperty(document, "visibilityState", {
			configurable: true,
			get: () => "visible",
		});
		stop();
	});

	it("stop flushes a pending edit and unsubscribes", () => {
		const store = makeStore();
		const stop = startAutosave(store, get, 1000);
		store.setState({ nodes: [node("a")] });
		stop();
		expect(readAutosave(get)?.nodes).toHaveLength(1);
		store.setState({ nodes: [node("a"), node("b")] });
		vi.advanceTimersByTime(5000);
		expect(readAutosave(get)?.nodes).toHaveLength(1);
	});
});
