"""Score the supply-crossflow Cd against Kaddah (1964), Vol. 1 Tables 11-30.

WHAT IS COMPARED: ACCURACY, CROSS-SOURCE. Kaddah's rig (Cranfield, 1964) is
not McGreehan and Schotsch's, nor Rohde's, nor NACA Lewis's, so a miss is
the model's limitation, never a fidelity bug. Nothing is fitted or tuned.

THE DATA ARE THE THESIS'S OWN PRINTED TABLES (data/kaddah1964/), kept as
printed in the thesis's units. Everything below that turns them into the
correlation's coordinates is code, not data:

* K is the velocity head ratio, (P1 - pj)/(P1 - p1) -- Rohde's VHR and TN
  3663's Phi -- and Kaddah's C is referenced to the approach TOTAL pressure
  (his Vj is from P1 - pj). conventions.py declares it Cd_total and the
  registered Cd_total -> Cd_static conversion, sqrt(K/(K-1)), is applied.
  The relative error is the same on either basis (one factor multiplies
  both sides), so it is reported once and the basis is stated.
* U1/Vi, TWO DEFINITIONS, always both and always labelled, never one
  silently. Incompressibly Vi/V1 = sqrt(K - 1) (Vi from the duct static):
    - "approach": U1 = Kaddah's approach V1  ->  U1/Vi = 1/sqrt(K - 1)
    - "mean": combaero's (m_a + m_b)/(2 rho A), the mean of the approaching
      and departing duct flow (crossflow_velocity_ratio)  ->
      U1/Vi = (1 - alpha/2)/sqrt(K - 1)
  Kaddah's bleed runs 2.7-100%, so the two differ by up to 2x; in Rohde
  and TN 3663 the bleed was small and they coincided.

TWO MODES, as orifice_runner:

* "eq17": McGreehan-Schotsch Eq. (17) on a baseline anchored to each
  table's own no-cross-flow block, the way the source validates Eq. (17)
  on Rohde ("a set baseline point at U1/Vi = 0"). The no-cross-flow rows
  still carry a small U1/Vi (the duct feeds the holes), so the baseline is
  the mean over that block of the cd_base that makes Eq. (17) return each
  row's measured Cd, under the same U1 definition. Only the 80 and 130 ft/s
  rows are scored; the no-cross-flow rows are the anchor.
* "chain": the shipped chain, mcgreehan_schotsch_1988_cd(Re, r/d = 0,
  L/d = t/d, U1/Vi), predicting the baseline too. Circular holes only: the
  chain needs a bore, and no slot is in its range. Every row is scored.

EXCLUSIONS AND READ CONFIDENCE come from data/kaddah1964/cell_flags.csv and
tables.yaml, never chosen here: a row whose scoring inputs (C, K, alpha, Lj,
Lth) carry `printed_inconsistent`, and a table marked
`excluded_from_scoring`, are excluded by default; a row whose scoring inputs
carry `sequence_judged` is scored but reported separately.

The jet angle phi is transcribed but not scored here: combaero has no
jet-angle closure yet (#813).
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Literal

import combaero as cb
import yaml

from validation.cooling import conventions

DATA = Path(__file__).parent / "data" / "kaddah1964"
SOURCE = "kaddah1964"
SET_EQ17 = "mcgreehan_schotsch_1988_crossflow_cd"
SET_CHAIN = "mcgreehan_schotsch_1988_cd"

U1Def = Literal["approach", "mean"]
Mode = Literal["eq17", "chain"]
U1_DEFS: tuple[U1Def, ...] = ("approach", "mean")
MODES: tuple[Mode, ...] = ("eq17", "chain")

# Columns the Cd score reads (directly, or through C = Lj/Lth).
CD_INPUTS = ("C", "K", "alpha_%", "Lj_lb/hr", "Lth_lb/hr")

# Unit conversions of the printed values (definitions, not data).
LB_PER_FT3_TO_KG_PER_M3 = 16.018463
FT_TO_M = 0.3048
IN_TO_M = 0.0254
IN_HG_TO_PA = 3386.389

# Velocity-head-ratio bands, as orifice_runner reports Rohde: the
# total-to-static factor is 1.41 at K 2 and 1.01 at K 50, so one pooled
# number would mix a near-exact comparison with a conversion-dominated one.
K_BANDS = ((1.0, 2.0), (2.0, 10.0), (10.0, math.inf))


@dataclass(frozen=True)
class Row:
    table: int
    code: str
    block: str
    row: int
    printed: dict[str, str]
    flags: dict[str, str] = field(default_factory=dict)  # column -> flag

    def f(self, col: str) -> float:
        return float(self.printed[col])

    @property
    def alpha(self) -> float:
        return self.f("alpha_%") / 100.0

    @property
    def excluded(self) -> bool:
        return any(self.flags.get(c) == "printed_inconsistent" for c in CD_INPUTS)

    @property
    def sequence_judged(self) -> bool:
        return any(self.flags.get(c) == "sequence_judged" for c in CD_INPUTS)


@dataclass(frozen=True)
class Table:
    meta: dict
    rows: tuple[Row, ...]

    @property
    def code(self) -> str:
        return self.meta["code"]

    @property
    def envelope(self) -> str:
        return self.meta["liner_envelope"]

    @property
    def excluded(self) -> bool:
        return bool(self.meta.get("excluded_from_scoring", False))


@lru_cache(maxsize=1)
def load_meta() -> dict:
    return yaml.safe_load((DATA / "tables.yaml").read_text())


@lru_cache(maxsize=1)
def load_flags() -> dict[tuple[int, str, int], dict[str, str]]:
    out: dict[tuple[int, str, int], dict[str, str]] = {}
    with open(DATA / "cell_flags.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            key = (int(r["table"]), r["block"], int(r["row"]))
            out.setdefault(key, {})[r["column"]] = r["flag"]
    return out


def _read_rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(line for line in fh if not line.startswith("#")))


@lru_cache(maxsize=1)
def load_tables() -> tuple[Table, ...]:
    meta = load_meta()
    flags = load_flags()
    out = []
    for tm in meta["tables"]:
        rows = []
        for r in _read_rows(DATA / tm["file"]):
            key = (tm["table"], r["block"], int(r["row"]))
            rows.append(Row(tm["table"], tm["code"], r["block"], int(r["row"]), r, flags.get(key, {})))
        out.append(Table(tm, tuple(rows)))
    return tuple(out)


# ---------------------------------------------------------------------------
# Kaddah's coordinates -> the correlation's.
# ---------------------------------------------------------------------------


def _source_like():
    from types import SimpleNamespace

    return SimpleNamespace(label=SOURCE, y_axis="Cd", source=SimpleNamespace(name=SOURCE))


def static_factor(k: float) -> float:
    """Cd_total -> Cd_static, the REGISTERED conversion (conventions.py)."""
    res = conventions.resolve(_source_like(), SET_CHAIN)
    if res.status != conventions.STATUS_CONVERTED:
        raise ValueError(f"{SOURCE}: {res.detail}")
    conversion = conventions.CONVERSIONS[("Cd_total", "Cd_static")]
    assert conversion.factor is not None
    return conversion.factor(k)


def crossflow_ratio(row: Row, u1: U1Def) -> float:
    """U1/Vi under the named U1 definition (see the module docstring)."""
    base = 1.0 / math.sqrt(row.f("K") - 1.0)
    if u1 == "approach":
        return base
    return (1.0 - row.alpha / 2.0) * base


def measured_cd_static(row: Row) -> float:
    return row.f("C") * static_factor(row.f("K"))


def hole_reynolds(row: Row, d_in: float) -> float:
    """rhoj Vj d / mu at the approach static temperature (p1 = 33 in Hg)."""
    p1 = load_meta()["rig"]["approach_static_in_Hg_abs"] * IN_HG_TO_PA
    rho1 = row.f("rho1_lb/cu.ft") * LB_PER_FT3_TO_KG_PER_M3
    x_air = list(cb.species.dry_air())
    t1 = p1 / (rho1 * cb.specific_gas_constant(x_air))
    mu = cb.viscosity(t1, p1, x_air)
    rhoj = row.f("rhoj_lb/cu.ft") * LB_PER_FT3_TO_KG_PER_M3
    return rhoj * row.f("Vj_ft/sec") * FT_TO_M * d_in * IN_TO_M / mu


def _eq17_baseline(table: Table, u1: U1Def) -> float:
    """Mean over the no-cross-flow block of the cd_base for which Eq. (17)
    returns each row's measured static Cd at that row's own U1/Vi."""
    vals = []
    for r in table.rows:
        if r.block != "noxflow" or r.excluded:
            continue
        x = crossflow_ratio(r, u1)
        target = measured_cd_static(r)
        lo, hi = 0.05, 1.5
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if cb.mcgreehan_schotsch_1988_crossflow_cd(mid, x) < target:
                lo = mid
            else:
                hi = mid
        vals.append(0.5 * (lo + hi))
    if not vals:
        raise ValueError(f"table {table.meta['table']}: no usable no-cross-flow row")
    return sum(vals) / len(vals)


@dataclass(frozen=True)
class Record:
    row: Row
    envelope: str
    u1: U1Def
    mode: Mode
    x: float  # U1/Vi
    k: float
    measured: float  # Cd_static
    predicted: float | None
    reason: str | None = None

    @property
    def rel_error(self) -> float | None:
        if self.predicted is None:
            return None
        return self.predicted / self.measured - 1.0

    @property
    def k_band(self) -> str:
        for lo, hi in K_BANDS:
            if lo <= self.k < hi:
                return f"K {lo:g}-{hi:g}"
        return "K ?"


def run_table(table: Table, u1: U1Def, mode: Mode) -> list[Record]:
    out: list[Record] = []
    if table.excluded:
        return out
    base = _eq17_baseline(table, u1) if mode == "eq17" else None
    d_in = table.meta["dimensions_in"].get("D_in")
    for r in table.rows:
        if r.excluded:
            continue
        if mode == "eq17" and r.block == "noxflow":
            continue
        x = crossflow_ratio(r, u1)
        meas = measured_cd_static(r)
        if mode == "chain" and d_in is None:
            out.append(Record(r, table.envelope, u1, mode, x, r.f("K"), meas, None,
                              "chain needs a bore diameter; slots are outside its range"))
            continue
        if mode == "eq17":
            pred = cb.mcgreehan_schotsch_1988_crossflow_cd(base, x)
        else:
            re = hole_reynolds(r, d_in)
            pred = cb.mcgreehan_schotsch_1988_cd(re, 0.0, table.meta["t_over_d"], x)
        out.append(Record(r, table.envelope, u1, mode, x, r.f("K"), meas, pred))
    return out


def run_all(u1: U1Def, mode: Mode) -> list[Record]:
    out: list[Record] = []
    for t in load_tables():
        out.extend(run_table(t, u1, mode))
    return out


@dataclass(frozen=True)
class Score:
    n: int
    bias: float
    mae: float
    rms: float


def score(records: list[Record]) -> Score:
    errs = [r.rel_error for r in records if r.rel_error is not None]
    n = len(errs)
    if n == 0:
        return Score(0, math.nan, math.nan, math.nan)
    return Score(n, sum(errs) / n, sum(abs(e) for e in errs) / n, math.sqrt(sum(e * e for e in errs) / n))


def report() -> str:
    """Plain-text report: ACCURACY, CROSS-SOURCE, every split kept apart."""
    lines = [
        "Kaddah (1964) Tables 11-30 vs McGreehan-Schotsch -- ACCURACY, CROSS-SOURCE",
        "Cd on the duct-STATIC basis (Kaddah's total-referenced C x sqrt(K/(K-1)));",
        "relative errors are identical on the total basis.",
        "",
    ]
    for mode in MODES:
        for u1 in U1_DEFS:
            recs = run_all(u1, mode)
            lines.append(f"[{mode}]  U1 = {'Kaddah approach V1' if u1 == 'approach' else 'combaero mean of approaching and departing flow'}")
            for env in ("inside", "marginal", "outside"):
                for seq in (False, True):
                    sel = [r for r in recs if r.envelope == env and r.row.sequence_judged == seq and r.predicted is not None]
                    if not sel:
                        continue
                    tag = "sequence-judged inputs" if seq else "clean inputs"
                    s = score(sel)
                    lines.append(f"  {env:8s} {tag:22s} n={s.n:3d} bias={s.bias:+.1%} MAE={s.mae:.1%}")
                    if not seq:
                        for lo, hi in K_BANDS:
                            b = [r for r in sel if lo <= r.k < hi]
                            if b:
                                sb = score(b)
                                lines.append(f"      K {lo:g}-{hi:g}: n={sb.n:3d} bias={sb.bias:+.1%} MAE={sb.mae:.1%}")
            lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Jet angle phi: REPORTED, NOT SCORED. combaero has no jet-angle closure yet
# (#813); these are the two printed equations a future one would start from,
# evaluated on Kaddah's own printed columns. Nothing here enters the solver.
# ---------------------------------------------------------------------------

# Lefebvre & Ballal plot NREC's Eq. 97 against Kaddah's phi with Cd_inf 0.625
# (their Fig. 4.5); the alternative anchors Cd_inf on each table's own
# no-cross-flow mean C.
NREC_CD_INF_LEFEBVRE = 0.625


def phi_report1373_eq9(row: Row) -> float:
    """Report 1373 Eq. 9, cos phi = V_A/V_j: the jet keeps the approach
    velocity. V_A = Kaddah's V1, V_j = his printed (ideal) Vj."""
    return math.degrees(math.acos(min(1.0, row.f("V1_ft/sec") / row.f("Vj_ft/sec"))))


def phi_nrec_eq97(row: Row, cd_inf: float) -> float:
    """NREC CR-72374 Eq. 97 (= Lefebvre & Ballal Eq. 4.17), sin^2 phi =
    Cd/Cd_inf, with Kaddah's printed C."""
    return math.degrees(math.asin(math.sqrt(min(1.0, row.f("C") / cd_inf))))


def jet_angle_report() -> str:
    lines = ["Jet angle phi, crossflow rows (80 and 130 ft/s) -- REPORTED, NOT SCORED (#813)",
             "error = predicted - measured, degrees"]
    tables = load_tables()
    for env in ("inside", "marginal", "outside"):
        res: dict[str, list[float]] = {"eq9": [], "eq97_0.625": [], "eq97_own": []}
        for t in tables:
            if t.envelope != env or t.excluded:
                continue
            own = [r.f("C") for r in t.rows if r.block == "noxflow" and not r.excluded]
            cd_inf = sum(own) / len(own)
            for r in t.rows:
                if r.block == "noxflow" or r.flags.get("phi_deg") in ("printed_inconsistent", "sequence_judged"):
                    continue
                m = r.f("phi_deg")
                res["eq9"].append(phi_report1373_eq9(r) - m)
                res["eq97_0.625"].append(phi_nrec_eq97(r, NREC_CD_INF_LEFEBVRE) - m)
                res["eq97_own"].append(phi_nrec_eq97(r, cd_inf) - m)
        for k, v in res.items():
            if v:
                n = len(v)
                lines.append(f"  {env:8s} {k:11s} n={n:3d} bias={sum(v) / n:+5.1f} MAE={sum(abs(e) for e in v) / n:4.1f}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(report())
    print(jet_angle_report())
