import { describe, expect, it } from "vitest";
import {
	autoInclination,
	exitCrossflowLabel,
	HOLE_CD_RULES,
	inclinationLabel,
	liftLooseHoleCd,
	resolvedCorrelation,
	withCorrelation,
} from "./holeCd";
import { migrateNodeData } from "./legacyKeys";

describe("hole_cd: one object per hole node (#1098)", () => {
	it("resolves a null correlation to the node type's default", () => {
		expect(resolvedCorrelation("orifice", {})).toBe("ReaderHarrisGallagher");
		expect(resolvedCorrelation("film_cooling", undefined)).toBe(
			"IdelchikThick",
		);
		expect(resolvedCorrelation("impingement_array", {})).toBe("fixed");
		expect(resolvedCorrelation("effusion_liner", {})).toBe("McGreehanSchotsch");
	});

	it("offers only the terms each node type builds", () => {
		expect(HOLE_CD_RULES.orifice.inclination).toEqual([]);
		expect(HOLE_CD_RULES.slot_film.exitCrossflow).toEqual([]);
		expect(HOLE_CD_RULES.film_cooling.exitCrossflow).toContain("rowbury_2001");
		expect(HOLE_CD_RULES.impingement_array.exitCrossflow).not.toContain(
			"rowbury_2001",
		);
		expect(HOLE_CD_RULES.effusion_liner.inletCrossflow).toContain(
			"mcgreehan_schotsch",
		);
		expect(HOLE_CD_RULES.effusion_plate.inletCrossflow).toEqual([]);
	});

	it("says what auto is on this node", () => {
		const cyl = { hole_type: "cylindrical", angle_deg: 30 };
		expect(autoInclination("film_cooling", {}, cyl)).toBe("idelchik_3_2");
		expect(autoInclination("film_cooling", {}, { angle_deg: 30 })).toBe("none");
		expect(inclinationLabel("auto", "film_cooling", {}, {})).toBe(
			"Auto: none (shaped hole)",
		);
		expect(autoInclination("effusion_plate", {}, {})).toBe("none");
		expect(autoInclination("effusion_plate", {}, { angle_deg: 30 })).toBe(
			"idelchik_3_2",
		);
		expect(
			autoInclination(
				"effusion_plate",
				{ correlation: "Lichtarowicz" },
				{ angle_deg: 30 },
			),
		).toBe("none");
		expect(exitCrossflowLabel("auto", "impingement_array")).toContain(
			"Rogers & Hersh",
		);
		expect(exitCrossflowLabel("auto", "film_cooling")).toBe("Auto: none");
	});

	it("leaving IdelchikThick resets a forced inclination (#1094)", () => {
		const hc = { correlation: "IdelchikThick", inclination: "idelchik_3_2" };
		expect(withCorrelation("film_cooling", hc as any, "Lichtarowicz")).toEqual({
			correlation: "Lichtarowicz",
			inclination: "auto",
		});
		const kept = { correlation: "IdelchikThick", inclination: "none" };
		expect(
			withCorrelation("film_cooling", kept as any, "Lichtarowicz"),
		).toEqual({ correlation: "Lichtarowicz", inclination: "none" });
	});

	it("leaving McGreehan-Schotsch resets a forced duct crossflow", () => {
		const hc = {
			correlation: "McGreehanSchotsch",
			inlet_crossflow: "mcgreehan_schotsch",
		};
		expect(
			withCorrelation("effusion_liner", hc as any, "IdelchikThick"),
		).toEqual({ correlation: "IdelchikThick", inlet_crossflow: "auto" });
	});

	it("carries a Cd with 'fixed' only", () => {
		expect(withCorrelation("orifice", {}, "fixed")).toEqual({
			correlation: "fixed",
			Cd: 0.6,
		});
		expect(withCorrelation("impingement_array", {}, null)).toEqual({
			correlation: null,
			Cd: 0.79,
		});
		expect(
			withCorrelation("orifice", { correlation: "fixed", Cd: 0.7 }, "fixed"),
		).toEqual({ correlation: "fixed", Cd: 0.7 });
		expect(
			withCorrelation(
				"orifice",
				{ correlation: "fixed", Cd: 0.7 },
				"IdelchikThick",
			),
		).toEqual({ correlation: "IdelchikThick" });
	});
});

describe("a v0.11.1 file's loose Cd fields are lifted into hole_cd", () => {
	it("keeps the correlation and the fixed Cd, drops an unused Cd", () => {
		expect(
			liftLooseHoleCd(
				{ diameter: 0.05, correlation: "IdelchikThick", Cd: 0.6 },
				"orifice",
			),
		).toEqual({ diameter: 0.05, hole_cd: { correlation: "IdelchikThick" } });
		expect(
			liftLooseHoleCd({ correlation: "fixed", Cd: 0.7 }, "slot_film"),
		).toEqual({ hole_cd: { correlation: "fixed", Cd: 0.7 } });
		expect(liftLooseHoleCd({ Cd: 0.6 }, "orifice")).toEqual({
			hole_cd: { correlation: "ReaderHarrisGallagher" },
		});
	});

	it("an impingement array's cd_crossflow becomes exit_crossflow", () => {
		expect(
			migrateNodeData(
				{ correlation: "fixed", Cd: 0.79, cd_crossflow: "none", n_rows: 4 },
				"impingement_array",
			),
		).toEqual({
			n_rows: 4,
			hole_cd: { correlation: "fixed", Cd: 0.79, exit_crossflow: "none" },
		});
		expect(migrateNodeData({}, "impingement_array")).toEqual({});
	});

	it("leaves new files and other node types alone", () => {
		const d = { hole_cd: { correlation: "Lichtarowicz" } };
		expect(liftLooseHoleCd(d, "film_cooling")).toBe(d);
		const loss = { correlation_type: "constant_head", zeta: 1 };
		expect(liftLooseHoleCd(loss, "discrete_loss")).toBe(loss);
		expect(liftLooseHoleCd({ correlation: "x" }, "effusion_liner")).toEqual({
			correlation: "x",
		});
	});
});
