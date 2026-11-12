import { describe, expect, it } from "vitest";
import { nodeTypes } from "../components/flowTypes";
import {
	expectedResultKeys,
	formatGeneric,
	formatMassFlow,
	formatSignificant,
	QUANTITY_CATALOGUE,
	quantityMeta,
	resultKinds,
	siUnit,
	THERMAL_WALL_KIND,
} from "./quantities";

// resultKeys.generated.json is what the backend emits per palette kind
// (scripts/gui_result_keys.py; python/tests/test_gui_result_keys_sync.py
// keeps it current). These tests tie the hand-kept catalogue to it.

const rawKey = (k: string) => (k.startsWith("state.") ? k.slice(6) : k);
const emitted = new Set(
	resultKinds().flatMap((kind) => expectedResultKeys(kind).map(rawKey)),
);

// Catalogue entries the generator's networks do not reach, each with why it
// is kept. Anything else never emitted is a stale entry: delete it.
const NOT_REACHED: Record<string, string> = {
	rayleigh_omitted:
		"heated compressible duct with no derivable exit state (#551); since #672 every GUI topology derives one",
	rayleigh_omitted_frac: "as rayleigh_omitted",
};

// A display unit may differ from the backend's SI unit only by a display
// scaling (the format converts) or by notation.
const DISPLAY_OF_SI: Record<string, string[]> = {
	"-": ["—", "%", "% Pt", "kg/kg"],
	Pa: ["bar"],
	"W/m²": ["kW/m²"],
	"J/kg/K": ["J/kg·K"],
	"W/m/K": ["W/m-K"],
	Pa·s: ["Pa·s"],
};

describe("QUANTITY_CATALOGUE against what the backend emits", () => {
	it("every catalogue entry is emitted by at least one kind", () => {
		const stale = Object.keys(QUANTITY_CATALOGUE).filter(
			(k) => !emitted.has(k) && !(k in NOT_REACHED),
		);
		expect(stale).toEqual([]);
	});

	it("an entry kept as not reached is really not emitted", () => {
		// Once a network reaches it, it is checked like any other entry.
		for (const k of Object.keys(NOT_REACHED)) {
			expect(emitted.has(k)).toBe(false);
			expect(QUANTITY_CATALOGUE[k]).toBeDefined();
		}
	});

	it("every catalogue unit agrees with the backend's SI unit", () => {
		const wrong = Object.entries(QUANTITY_CATALOGUE)
			.filter(([k]) => emitted.has(k))
			.filter(([k, cat]) => {
				const si = siUnit(k);
				return (
					si === undefined ||
					(cat.unit !== si && !(DISPLAY_OF_SI[si] ?? []).includes(cat.unit))
				);
			})
			.map(([k, cat]) => `${k}: catalogue ${cat.unit}, backend ${siUnit(k)}`);
		expect(wrong).toEqual([]);
	});

	it("every emitted key has an SI unit", () => {
		expect([...emitted].filter((k) => siUnit(k) === undefined)).toEqual([]);
	});
});

describe("the generated key list covers the palette", () => {
	it("every solvable palette kind and the thermal wall have keys", () => {
		const kinds = Object.keys(nodeTypes).filter((k) => k !== "probe");
		for (const kind of [...kinds, THERMAL_WALL_KIND]) {
			expect(expectedResultKeys(kind).length, kind).toBeGreaterThan(0);
		}
		expect(expectedResultKeys("probe")).toEqual([]);
		expect(expectedResultKeys(undefined)).toEqual([]);
	});

	it("node kinds report their state under state.*, elements at the top", () => {
		expect(expectedResultKeys("plenum")).toContain("state.P");
		expect(expectedResultKeys("plenum")).not.toContain("P");
		expect(expectedResultKeys("channel")).toContain("P_in");
		expect(expectedResultKeys("channel")).not.toContain("state.P");
	});
});

describe("quantityMeta", () => {
	it("a catalogued key uses its entry, with or without state.", () => {
		expect(quantityMeta("P").label).toBe("P");
		expect(quantityMeta("P").unit).toBe("bar");
		expect(quantityMeta("P").format(1e5)).toBe("1.000");
		expect(quantityMeta("state.P").unit).toBe("bar");
		expect(quantityMeta("P").catalogued).toBe(true);
	});

	it("an uncatalogued key keeps its name and shows its SI unit", () => {
		const meta = quantityMeta("port_0_P");
		expect(meta.catalogued).toBe(false);
		expect(meta.label).toBe("port_0_P");
		expect(meta.unit).toBe("Pa");
		expect(quantityMeta("NTU").unit).toBe("—");
	});

	it("a deprecated alias shows its new key's meta, marked", () => {
		const m = quantityMeta("U1_over_Vi_max");
		const n = quantityMeta("Cd_inlet_crossflow_U1_over_Vi_max");
		expect(m.label).toBe(
			`${n.label} (deprecated, use Cd_inlet_crossflow_U1_over_Vi_max)`,
		);
		expect(m.unit).toBe(n.unit);
		expect(m.catalogued).toBe(true);
	});

	it("a key nobody knows still renders", () => {
		const meta = quantityMeta("no_such_key");
		expect(meta.unit).toBe("—");
		expect(meta.format(1.5)).toBe("1.5");
	});
});

describe("formatGeneric", () => {
	it("is readable across magnitudes", () => {
		expect(formatGeneric(0)).toBe("0");
		expect(formatGeneric(101325)).toBe("101330");
		expect(formatGeneric(0.123456)).toBe("0.12346");
		expect(formatGeneric(2.5e-5)).toBe("2.500e-5");
		expect(formatGeneric(3.2e7)).toBe("3.200e+7");
		expect(formatGeneric(-1)).toBe("-1");
	});
});

// Components that read a result key directly (not through the catalogue).
// A key renamed on the backend leaves such a read silently undefined, as
// AreaChangeNode's `result.ratio` was after the area_ratio rename.
const SOURCES = import.meta.glob(["../components/**/*.tsx"], {
	query: "?raw",
	import: "default",
	eager: true,
}) as Record<string, string>;

// Read directly but not scalars: lists, strings and the state block.
const NON_SCALAR_READS = new Set([
	"state",
	"X",
	"Y",
	"T_interface",
	"rows_Re_j",
	"rows_Nu",
	"rows_htc",
	"rows_Gc_Gj",
	"rows_valve_dPt", // flow splitter (#903): the drop at each outlet
	"h_gas_unblown_model",
]);

describe("direct result reads in the components", () => {
	it("every key a component reads is one the backend emits", () => {
		const unknown: string[] = [];
		for (const [file, text] of Object.entries(SOURCES)) {
			for (const m of text.matchAll(
				/result\??\.(?:state\??\.)?([A-Za-z_][A-Za-z0-9_]*)/g,
			)) {
				const key = m[1];
				if (
					!emitted.has(key) &&
					!NON_SCALAR_READS.has(key) &&
					!(key in NOT_REACHED)
				) {
					unknown.push(`${file}: result.${key}`);
				}
			}
		}
		expect(unknown).toEqual([]);
	});
});

describe("formatMassFlow", () => {
	it("keeps four decimals for an ordinary flow", () => {
		expect(formatMassFlow(0.2284676518958348)).toBe("0.2285");
		expect(formatMassFlow(12.5)).toBe("12.5000");
		expect(formatMassFlow(0.001)).toBe("0.0010");
		expect(formatMassFlow(-0.05)).toBe("-0.0500");
		expect(formatMassFlow(0)).toBe("0.0000");
	});

	it("goes exponential below 1e-3 kg/s instead of rounding to a flow", () => {
		// A blocked fraction loss (#1006): 5.5e-5 kg/s read as "0.0001".
		expect(formatMassFlow(5.5e-5)).toBe("5.50e-5");
		expect(formatMassFlow(-5.5e-5)).toBe("-5.50e-5");
		expect(formatMassFlow(4e-6)).toBe("4.00e-6");
		expect(formatMassFlow(9.99e-4)).toBe("9.99e-4");
	});

	it("shows float residue as zero, not as a flow", () => {
		// An unfed effusion plate's card read "-4.51e-27 kg/s".
		expect(formatMassFlow(-4.51e-27)).toBe("0.0000");
		expect(formatMassFlow(9e-13)).toBe("0.0000");
		expect(formatMassFlow(-0)).toBe("0.0000");
		expect(formatMassFlow(2e-12)).toBe("2.00e-12");
		expect(formatMassFlow(5.53e-5)).toBe("5.53e-5");
	});

	it("is the format of every kg/s catalogue entry", () => {
		const massFlows = Object.entries(QUANTITY_CATALOGUE).filter(
			([, cat]) => cat.unit === "kg/s",
		);
		expect(massFlows.length).toBeGreaterThan(5);
		for (const [k, cat] of massFlows) {
			expect(cat.format(5.5e-5), k).toBe("5.50e-5");
			expect(cat.format(0.25), k).toBe("0.2500");
		}
	});
});

describe("plate-limit and fraction-loss flags have labels", () => {
	it("shows a name, not the raw key", () => {
		expect(quantityMeta("T_wall_limit").label).toBe("T wall limit");
		expect(quantityMeta("T_wall_limit").unit).toBe("K");
		expect(quantityMeta("T_wall_limit").format(1574.15)).toBe("1574");
		for (const k of [
			"wall_above_limit",
			"no_cooling_path",
			"fraction_unanchored",
		]) {
			const meta = quantityMeta(k);
			expect(meta.catalogued, k).toBe(true);
			expect(meta.label, k).not.toBe(k);
			expect(meta.format(1), k).toBe("yes");
			expect(meta.format(0), k).toBe("no");
		}
	});
});

describe("formatSignificant", () => {
	it("shows a solved m_ref at 6 significant figures", () => {
		// The "Use the solved flow as m_ref" button stores |m_dot| exactly.
		expect(formatSignificant(0.2284676518958348)).toBe("0.228468");
		expect(formatSignificant(0.25)).toBe("0.25");
		expect(formatSignificant(1234567.8)).toBe("1234570");
		expect(formatSignificant(5.5e-5)).toBe("0.000055");
		expect(formatSignificant(0.2284676518958348, 4)).toBe("0.2285");
	});
});
