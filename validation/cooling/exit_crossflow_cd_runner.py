"""The opt-in exit-crossflow Cd term (Rowbury, Oldfield & Lock 2001, #1074).

THREE QUESTIONS, NEVER ONE NUMBER (docs/VALIDATION_POLICY.md). Nothing is
fitted or tuned here, and every series stays `scores: null`: this runner
REPORTS, it does not score into the scorecard.

1. FIDELITY -- does combaero's term reproduce the AUTHORS' OWN numbers?
   The CHTT data the coefficients were fitted on (Figs. 4-8) are printed
   only as small plots, and the engine-vane data are proprietary, so
   fidelity rests on the one case the authors COMPUTED and printed with
   reproducible inputs: their Fig. 9 predictions for Gritsch et al.'s 30 deg
   cylindrical hole at Ma_m 0.3, 0.6, 1.2 (suction-surface set, theta 30,
   alpha 0, Cd,no = Lichtarowicz Eq. 1 at l/d 6 = 0.776, which is their flat
   "Predict Mm = 0" line). The predicted points were READ BY EYE off a
   600 dpi render of the printed figure (FIG9_PREDICTED, +-0.005 in Cd,
   +-0.005 in PR). They live in this runner's source only, NOT in
   validation/cooling/data/ (no reviewed digitisation, no metadata), so
   they are a reproduction check, not a dataset. Fidelity therefore means:
   the implementation reproduces the authors' generalised curves to the
   reading accuracy. It is judged on Eq. 8 AS PRINTED (`as_printed`); the
   shipped term also fades Eq. 8's offset beyond the paper's fitted k
   (k_perp > 3.3, include/crossflow_cd.h), our extrapolation, which the
   authors' Ma_m 0.3 markers above PR ~1.5 already reach; both are reported.

2. ACCURACY -- CROSS-SOURCE. Gritsch et al. 1998 (Karlsruhe, digitised,
   validation/cooling/data/gritsch1998, Fig. 3: Ma_c 0, Ma_m 0.3/0.6/1.2) and
   Hay et al. 1983 (Bath, validation/cooling/data/hay1983, Mc 0, Minf 0.31:
   Figs. 6, 8, 12). Neither lab is Rowbury's; Gritsch's 30 deg hole is the
   geometry the authors themselves predict in Fig. 9, but the coefficients
   come from the CHTT, so a miss is the model's limitation, never a fidelity
   bug. Cd,no is combaero's own (IdelchikThick, the inclination term off
   and on), so the before/after isolates the exit term on what a user gets.

3. TUNING: none. `Nu_multiplier`-style factors do not exist for this term
   and none is introduced.

THE POINT SOLVE MIRRORS THE ELEMENT. FilmCoolingElement with
cd_exit_crossflow='rowbury_2001' solves m = Cd A_eff-nozzle flow with
Cd = rowbury(Cd0(Re(m)), k(m)), k = G_h^2 f / (rho_h rho_m u_m^2), f =
sin(theta) at alpha 0, rho_h the isentropic exit-plane density
(rowbury_2001_exit_density; OrificeElement._exit_crossflow). The element reads
rho_m u_m^2 from a MomentumChamberNode, which is subsonic, and Gritsch's
Ma_m 1.2 is not; so each point here is solved on the papers' own ideal flow
(isentropic, coolant total to mainstream static, compressible to choking)
with rho_m u_m^2 = gamma p_m Ma_m^2, and the same C++ calls the element
makes. tests (test_exit_crossflow_cd.py) pin that the element in a network
gives this point solve's Cd at its solved state. A second reason to solve
points: a hole into a moving chamber blows only while its supply Pt exceeds
the chamber's Pt (#1087), so the series' low-PR points inside p_m < p0c <
Pt_m (Ma_m 0.3 below PR 1.064, 0.6 below 1.276, 1.2 all of them) are not
reachable through the element in a network at all.

RIG BOUNDARY CONDITIONS (assumed, as in hole_cd_inclination_runner): exit
static 1 atm, air at 290 K (Gritsch) / 293 K (Hay), Tt_c = Tt_m.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cache

import combaero as cb
from validation.cooling.hole_cd_inclination_runner import GRITSCH, HAY, P_INF, ideal_flow
from validation.cooling.schema import SeriesMetadata, load_dataset, load_points

_X_AIR = list(cb.species.dry_air())
THICK = cb.DischargeCdCorrelation.Idelchik1966Thick
SUCTION = cb.Rowbury2001Set.CylindricalSuctionSurface

# Lichtarowicz et al. (1965) Eq. 1 as Rowbury print it (Eq. 1): the authors'
# own Cd,no for Gritsch's hole, l/d 6, Re_h > 2e4.
LICHTAROWICZ_CDU_L_OVER_D_6 = 0.827 - 0.0085 * 6.0

# Rowbury et al. Fig. 9 "Predict" markers, (Ma_m, PR, Cd), read by eye off a
# 600 dpi render of p. 263 (axis 1380 px per unit PR, 1594 px per unit Cd);
# +-0.005 in both. NOT a reviewed digitisation (see the module docstring).
FIG9_PREDICTED = (
    (0.3, 1.016, 0.643),
    (0.3, 1.028, 0.691),
    (0.3, 1.059, 0.739),
    (0.3, 1.090, 0.755),
    (0.3, 1.133, 0.766),
    (0.3, 1.197, 0.776),
    (0.3, 1.348, 0.786),
    (0.3, 1.552, 0.791),
    (0.6, 1.048, 0.558),
    (0.6, 1.068, 0.612),
    (0.6, 1.081, 0.639),
    (0.6, 1.117, 0.677),
    (0.6, 1.153, 0.702),
    (0.6, 1.206, 0.719),
    (0.6, 1.267, 0.736),
    (0.6, 1.356, 0.748),
    (1.2, 1.111, 0.429),
    (1.2, 1.147, 0.491),
    (1.2, 1.164, 0.517),
    (1.2, 1.248, 0.590),
    (1.2, 1.325, 0.630),
    (1.2, 1.403, 0.655),
    (1.2, 1.489, 0.673),
    (1.2, 1.578, 0.686),
    (1.2, 1.704, 0.701),
    (1.2, 1.829, 0.712),
)

GRITSCH_SERIES = ("fig3_Mac0p0_Mam0p3", "fig3_Mac0p0_Mam0p6", "fig3_Mac0p0_Mam1p2")
HAY_SERIES = ("fig6_open_circles", "fig8_open_circles", "fig12_open_circles")
# Hay's "parallel" and "perpendicular" name the COOLANT duct; at the exit
# every hole leans downstream with the mainstream (alpha = 0).
HAY_MINF = 0.31


def _gamma(T: float, p: float) -> float:
    return cb.complete_state(T, p, _X_AIR).thermo.gamma


@dataclass(frozen=True)
class Rig:
    d: float
    n: int
    l_over_d: float
    Tt: float
    theta: float


def _cd0(rig: Rig, re: float, inclination: bool) -> float:
    h = cb.DischargeHoleGeometry(d=rig.d, L=rig.l_over_d * rig.d)
    if inclination and rig.theta < 90.0:
        h.inclination_deg = rig.theta
    return float(cb.discharge_cd(THICK, h, cb.DischargeHoleState(Re=re)))


@cache
def solve_point(
    rig: Rig,
    pr: float,
    mach_m: float,
    exit_term: bool,
    inclination: bool,
    cd_no: float | None = None,
    rho_h_supply_tt: bool = False,
    as_printed: bool = False,
) -> tuple[float, float, bool]:
    """(Cd, k_perp, held) at one operating point.

    cd_no fixes the no-crossflow Cd (the authors' Lichtarowicz constant)
    instead of combaero's IdelchikThick at the iterate's Re.
    rho_h_supply_tt takes rho(Tt_c, p_m) instead of the isentropic exit
    density (the first version's choice, a sensitivity). as_printed
    evaluates Eq. 8 as the paper prints it (the k floor kept, no fade
    beyond the fitted k)."""
    m_ideal = ideal_flow(pr, rig.Tt, rig.d, rig.n)
    area = rig.n * math.pi * rig.d**2 / 4.0
    g_ideal = m_ideal / area
    pt = pr * P_INF
    mu = cb.transport_state(rig.Tt, pt, _X_AIR).mu
    if rho_h_supply_tt:
        rho_h = float(cb.density(rig.Tt, P_INF, _X_AIR))
    else:
        rho_h = float(cb.rowbury_2001_exit_density(pt, rig.Tt, P_INF, _X_AIR).rho)
    q_m = _gamma(rig.Tt, P_INF) * P_INF * mach_m**2

    def residual(cd: float) -> tuple[float, float, bool]:
        re = 4.0 * cd * m_ideal / rig.n / (math.pi * rig.d * mu)
        c0 = cd_no if cd_no is not None else _cd0(rig, re, inclination)
        if not exit_term or q_m <= 0.0:
            return cd - c0, math.nan, False
        k = cb.rowbury_2001_k_perp(cd * g_ideal, rho_h, q_m, rig.theta, 0.0).k
        if as_printed:
            c = cb.rowbury_2001_coefficients(SUCTION, rig.theta, 0.0)
            kf = cb.rowbury_2001_k_floor(SUCTION, rig.theta, 0.0)
            kk = max(k, kf)
            return cd - (1.0 / c0**2 + c.A * kk**-c.B - c.C) ** -0.5, k, k < kf
        r = cb.rowbury_2001_exit_crossflow_cd(c0, k, SUCTION, rig.theta, 0.0)
        return cd - r.Cd, k, bool(r.held)

    lo, hi = 1e-6, 1.0
    f_lo = residual(lo)[0]
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f = residual(mid)[0]
        if (f < 0.0) == (f_lo < 0.0):
            lo, f_lo = mid, f
        else:
            hi = mid
        if hi - lo < 1e-12:
            break
    cd = 0.5 * (lo + hi)
    _, k, held = residual(cd)
    return cd, k, held


# --- 1. fidelity ----------------------------------------------------------

GRITSCH_RIG = Rig(GRITSCH["d"], GRITSCH["n"], GRITSCH["l_over_d"], GRITSCH["Tt"], 30.0)


def fidelity_rows(
    rho_h_supply_tt: bool = False, as_printed: bool = True
) -> list[tuple[float, float, float, float]]:
    """(Ma_m, PR, read Cd, combaero Cd) on the authors' Fig. 9 predictions."""
    out = []
    for mach, pr, cd_read in FIG9_PREDICTED:
        cd, _, _ = solve_point(
            GRITSCH_RIG,
            pr,
            mach,
            True,
            False,
            cd_no=LICHTAROWICZ_CDU_L_OVER_D_6,
            rho_h_supply_tt=rho_h_supply_tt,
            as_printed=as_printed,
        )
        out.append((mach, pr, cd_read, cd))
    return out


# --- 2. accuracy ----------------------------------------------------------


@dataclass(frozen=True)
class Row:
    label: str
    pr: float
    measured: float
    cd: dict[tuple[bool, bool], float]  # (inclination, exit_term) -> Cd
    k_perp: float
    held: bool

    def err(self, inclination: bool, exit_term: bool) -> float:
        return self.cd[(inclination, exit_term)] / self.measured - 1.0


def _series(ds: list[SeriesMetadata], source: str, stems: tuple[str, ...]) -> list[SeriesMetadata]:
    by = {s.path.stem: s for s in ds if s.source.name == source}
    return [by[s] for s in stems]


def rows_for(series: SeriesMetadata) -> list[Row]:
    g = series.geometry or {}
    if series.source.name == "gritsch1998":
        rig, mach = GRITSCH_RIG, float(g["Ma_m"])
    else:
        rig = Rig(HAY["d"], HAY["n"], HAY["l_over_d"], HAY["Tt"], float(series.alpha_deg))
        mach = float(g["M_inf"])
    out = []
    for p in load_points(series):
        cds = {}
        k, held = math.nan, False
        for inc in (False, True):
            for ex in (False, True):
                cd, kk, hh = solve_point(rig, p.x, mach, ex, inc)
                cds[(inc, ex)] = cd
                if ex and not inc:
                    k, held = kk, hh
        out.append(Row(series.label, p.x, p.y, cds, k, held))
    return out


def accuracy(ds: list[SeriesMetadata] | None = None) -> list[tuple[str, list[Row]]]:
    ds = ds if ds is not None else load_dataset()
    sel = _series(ds, "gritsch1998", GRITSCH_SERIES) + _series(ds, "hay1983", HAY_SERIES)
    return [(s.label, rows_for(s)) for s in sel]


def summary(rows: list[Row], inclination: bool, exit_term: bool) -> tuple[float, float]:
    e = [r.err(inclination, exit_term) for r in rows]
    return sum(e) / len(e), sum(abs(x) for x in e) / len(e)


# --- 3. the DLR effusion rig (Behrendt et al. 2008), reported only -------


def behrendt_effect(ds: list[SeriesMetadata] | None = None) -> list[dict[str, float | str]]:
    """The term's effect on each Behrendt 2008 operating point, through the
    network (effusion_network_runner, chamber baseline): Cd, coolant flow and
    eta_overall off -> on. Reported, never scored: the term was not built
    for, or judged on, this rig."""
    from validation.cooling import effusion_network_runner as en

    ds = ds if ds is not None else load_dataset()
    out = []
    seen = set()
    for s in ds:
        if s.source.name != "behrendt2008" or s.y_axis != "eta_overall":
            continue
        key = tuple(sorted((s.geometry or {}).items()))
        if key in seen:
            continue
        seen.add(key)
        off = en.predict(s, "chamber")
        try:
            on = en.predict(s, "chamber", cd_exit_crossflow="rowbury_2001")
        except ValueError as exc:
            # A 20 deg hole is below the suction set's A = 0 crossing.
            out.append({"series": s.label, "refused": str(exc).split(": ", 2)[-1]})
            continue
        out.append(
            {
                "series": s.label,
                "Cd_off": off["Cd"],
                "Cd_on": on["Cd"],
                "m_ratio": on["m_dot"] / off["m_dot"],
                "k_perp": on["exit_crossflow_k_perp"],
                "eta_off": off["eta_overall"],
                "eta_on": on["eta_overall"],
            }
        )
    return out


def _main() -> None:
    print("1. FIDELITY: the authors' Fig. 9 predictions (eye-read, +-0.005), Cd,no 0.776")
    variants = (
        (False, True, "Eq. 8 as printed, rho_h isentropic exit (the element's)"),
        (False, False, "as shipped: offset faded beyond k_perp 3.3"),
        (True, True, "Eq. 8 as printed, rho_h = rho(Tt_c, p_m) (sensitivity)"),
    )
    for supply, printed, tag in variants:
        rows = fidelity_rows(supply, printed)
        iso = supply or not printed
        print(f"  {tag}")
        for mach in (0.3, 0.6, 1.2):
            sel = [r for r in rows if r[0] == mach]
            d = [r[3] - r[2] for r in sel]
            print(
                f"    Ma_m {mach:.1f}: n {len(sel)} mean dCd {sum(d) / len(d):+.4f} "
                f"max |dCd| {max(abs(x) for x in d):.4f}"
            )
            if not iso:
                for _, pr, read, cd in sel:
                    print(f"      PR {pr:.3f} read {read:.3f} combaero {cd:.3f}")
    print("\n2. ACCURACY, CROSS-SOURCE (Rowbury coefficients vs Karlsruhe / Bath data)")
    print(
        f"{'series':36s} {'n':>3s} {'k_perp':>11s}  "
        "bias(no-incl: off -> on) |e|        bias(incl: off -> on) |e|"
    )
    for label, rows in accuracy():
        ks = [r.k_perp for r in rows]
        line = f"{label:36s} {len(rows):3d} {min(ks):5.3f}-{max(ks):5.2f}"
        for inc in (False, True):
            b0, a0 = summary(rows, inc, False)
            b1, a1 = summary(rows, inc, True)
            line += f"  {b0:+6.1%} -> {b1:+6.1%} ({a0:5.1%} -> {a1:5.1%})"
        held = sum(r.held for r in rows)
        print(line + (f"  held {held}" if held else ""))
    print("\n3. BEHRENDT 2008 (DLR effusion), reported only: Cd, flow, eta off -> on")
    for r in behrendt_effect():
        if "refused" in r:
            print(f"  {r['series']:44s} refused: {r['refused']}")
            continue
        print(
            f"  {r['series']:44s} k_perp {r['k_perp']:6.2f} Cd {r['Cd_off']:.4f} -> "
            f"{r['Cd_on']:.4f} m x{r['m_ratio']:.4f} eta {r['eta_off']:.4f} -> {r['eta_on']:.4f}"
        )


if __name__ == "__main__":
    _main()
