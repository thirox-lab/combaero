"""The opt-in hole-inclination Cd term against Hay 1983 and Gritsch 1998.

WHAT IS COMPARED: ACCURACY, CROSS-SOURCE. The term is Idelchik's (diagram
3-2 composed into Eq. 4-20's inlet coefficient of the thick-walled hole, see
orifice::idelchik::inclination_* in include/orifice.h); neither paper is
Idelchik's, so a miss is the model's limitation, never a fidelity bug.
Nothing is fitted or tuned here, and every series stays `scores: null`: this
runner REPORTS, it does not score into the scorecard.

THE NETWORK, NOT THE FORMULA. Each operating point is a rig built as a
network and solved:

    PressureBoundary(coolant, Pt = PR p_inf) -> FilmCoolingElement
        (cylindrical, IdelchikThick, cd_inclination off / on)
        -> PressureBoundary(exit, p_inf)

and the model's Cd is its solved flow over the papers' ideal flow
(conventions.py `Cd_coolant_total_exit_static`: isentropic from the coolant
TOTAL pressure to the exit STATIC pressure, compressible up to choking, on
the inlet area pi D^2/4). Both boundaries are plena, so the element's own
reference (upstream total = static of a still plenum, to the downstream
plenum) is the papers' pair exactly and no conversion is needed; the
element's nozzle chokes where the ideal flow does.

IN SCOPE: the closure is plenum to plenum, so only NO-COOLANT-CROSSFLOW data
(Mc = Ma_c = 0) is scored:

* gritsch1998 fig3_Mac0p0_Mam0p0: 30 deg cylindrical, Ma_c 0, Ma_m 0 -- the
  one clean plenum-to-plenum series on the 30 deg hole.
* hay1983 Mc = 0 series at Minf ~0.31 (fig6 90 deg, fig8 30 deg parallel,
  fig12 30 deg perpendicular): no coolant crossflow but a mainstream
  crossflow on the EXIT, which no closure here carries. Reported as such.
* hay1983's printed carpet values at Mc = 0, Minf ~0.01-0.02 (Figs. 9, 10,
  13, 14; extractions/hay_1983_crossflow_cd.md, PR >= 1.6 plateau): the
  cleanest plenum-to-plenum numbers Hay give, read at PR 1.6 and 2.0.

OUT OF SCOPE (reported, not scored): every coolant-crossflow series (hay1983
Mc > 0, gritsch1998 Ma_c > 0), the mainstream-crossflow-only Gritsch series
(Ma_m > 0), shaped holes (Gritsch Figs. 7-14), and hay1983 Fig. 15 (its
abscissa is a velocity ratio).

RIG BOUNDARY CONDITIONS, NOT MODEL CHOICES (all reported, none tuned): the
exit static pressure is not printed by either paper; 1 atm is assumed (it
only moves the hole Reynolds number, and the Thick Cd is flat in Re above
1e4). Room-temperature air, 293 K (Hay) and 290 K (Gritsch). Hole length L
is taken along the axis (L/D 6 in both papers; Hay do not say which).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import combaero as cb
from combaero.network import FilmCoolingElement, FlowNetwork, NetworkSolver, PressureBoundary

from validation.cooling.schema import SeriesMetadata, load_dataset, load_points

P_INF = 101325.0  # [Pa] exit static, assumed (not printed)
_X_AIR = list(cb.species.dry_air())
_Y_AIR = list(cb.mole_to_mass(_X_AIR))
INCLINATION = ("none", "idelchik_3_2")

# hay1983: in-scope digitised series (Mc = 0) and the out-of-scope rest.
HAY_IN_SCOPE = ("fig6_open_circles", "fig8_open_circles", "fig12_open_circles")
GRITSCH_IN_SCOPE = ("fig3_Mac0p0_Mam0p0",)

# Hay et al.'s printed carpet Cd at Mc = 0, Minf ~0.01-0.02 (PR >= 1.6
# plateau), transcribed in extractions/hay_1983_crossflow_cd.md:
# (figure, rig, alpha_deg, Cd).
HAY_CARPET_NO_CROSSFLOW = (
    ("9", "parallel", 90.0, 0.81),
    ("10", "parallel", 30.0, 0.73),
    ("13", "perpendicular", 90.0, 0.78),
    ("14", "perpendicular", 30.0, 0.71),
)
CARPET_PR = (1.6, 2.0)

# Rig geometry. hay1983: D 3.5 mm, a row of four holes at S/D 2.5, L/D 6.
# gritsch1998: one D 10 mm hole, L/D 6 (pitch irrelevant to one hole).
HAY = {"d": 3.5e-3, "n": 4, "pitch_over_d": 2.5, "l_over_d": 6.0, "Tt": 293.0}
GRITSCH = {"d": 10.0e-3, "n": 1, "pitch_over_d": 6.0, "l_over_d": 6.0, "Tt": 290.0}


def ideal_flow(pr: float, Tt: float, d: float, n: int) -> float:
    """Isentropic flow from Pt = pr P_INF to P_INF, compressible up to
    choking, on n inlet areas (the papers' ideal flow)."""
    pt = pr * P_INF
    st = cb.complete_state(Tt, pt, _X_AIR).thermo
    k = st.gamma
    R = st.cp - st.cv
    p2 = max(P_INF, pt * (2.0 / (k + 1.0)) ** (k / (k - 1.0)))
    r = p2 / pt
    flux = pt * r ** ((k + 1.0) / (2.0 * k)) * math.sqrt(
        2.0 * k / ((k - 1.0) * R * Tt) * (r ** (-(k - 1.0) / k) - 1.0)
    )
    return flux * n * math.pi * d * d / 4.0


@lru_cache(maxsize=None)
def model_cd(
    pr: float, alpha: float, d: float, n: int, pitch_over_d: float, l_over_d: float, Tt: float,
    inclination: str,
) -> tuple[float, float, float, float]:
    """(Cd as the papers define it, the element's Cd, Re_hole, Cd_in_range)."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("cool", Pt=pr * P_INF, Tt=Tt, Y=_Y_AIR))
    g.add_node(PressureBoundary("exit", Pt=P_INF, Tt=Tt, Y=_Y_AIR))
    g.add_element(
        FilmCoolingElement(
            "holes",
            "cool",
            "exit",
            hole_diameter=d,
            holes_per_row=n,
            lateral_pitch=pitch_over_d * d,
            wall_thickness=l_over_d * d * math.sin(math.radians(alpha)),
            angle_deg=alpha,
            hole_type="cylindrical",
            correlation="IdelchikThick",
            cd_inclination=inclination,
        )
    )
    r = NetworkSolver(g).solve()
    if not r.get("__success__"):
        raise RuntimeError(f"hole Cd network did not solve: {r.get('__message__')}")
    diag = r["__element_diag__"]["holes"]
    m = abs(float(r["holes.m_dot"]))
    mu = cb.transport_state(Tt, pr * P_INF, _X_AIR).mu
    re = 4.0 * m / n / (math.pi * d * mu)
    return m / ideal_flow(pr, Tt, d, n), float(diag["Cd"]), re, float(diag["Cd_in_range"])


@dataclass(frozen=True)
class Row:
    label: str
    alpha: float
    pr: float
    measured: float
    off: float
    on: float
    re: float
    in_range_on: bool

    def err(self, which: str) -> float:
        return (self.off if which == "none" else self.on) / self.measured - 1.0


def _rig(series: SeriesMetadata) -> dict:
    return HAY if series.source.name == "hay1983" else GRITSCH


def rows_for(series: SeriesMetadata) -> list[Row]:
    rig = _rig(series)
    alpha = float(series.alpha_deg)
    out = []
    for p in load_points(series):
        args = (p.x, alpha, rig["d"], rig["n"], rig["pitch_over_d"], rig["l_over_d"], rig["Tt"])
        off, _, re, _ = model_cd(*args, "none")
        on, _, _, rng = model_cd(*args, "idelchik_3_2")
        out.append(Row(series.label, alpha, p.x, p.y, off, on, re, bool(rng)))
    return out


def carpet_rows() -> list[Row]:
    out = []
    for fig, rig_name, alpha, cd in HAY_CARPET_NO_CROSSFLOW:
        for pr in CARPET_PR:
            args = (pr, alpha, HAY["d"], HAY["n"], HAY["pitch_over_d"], HAY["l_over_d"], HAY["Tt"])
            off, _, re, _ = model_cd(*args, "none")
            on, _, _, rng = model_cd(*args, "idelchik_3_2")
            out.append(Row(f"hay1983/carpet_fig{fig}_{rig_name}", alpha, pr, cd, off, on, re, bool(rng)))
    return out


def in_scope(dataset: list[SeriesMetadata]) -> list[SeriesMetadata]:
    keep = {("hay1983", s) for s in HAY_IN_SCOPE} | {("gritsch1998", s) for s in GRITSCH_IN_SCOPE}
    return [s for s in dataset if (s.source.name, s.path.stem) in keep]


def out_of_scope(dataset: list[SeriesMetadata]) -> list[tuple[str, str]]:
    """Every other hay1983 / gritsch1998 series, with why it is not scored."""
    out = []
    for s in dataset:
        if s.source.name not in ("hay1983", "gritsch1998") or s in in_scope(dataset):
            continue
        g = s.geometry or {}
        if s.x_axis != "pc0_over_pinf":
            why = f"abscissa {s.x_axis}"
        elif s.source.name == "hay1983" and "Mc" not in g:
            why = "unassigned mark (series unknown)"
        elif g.get("hole_type", "cylindrical") != "cylindrical":
            why = "shaped hole"
        elif float(g.get("Mc", g.get("Ma_c", 0.0)) or 0.0) > 0.0:
            why = "coolant (inlet) crossflow"
        else:
            why = "mainstream (exit) crossflow"
        out.append((s.label, why))
    return out


def _interp(pts: list[tuple[float, float]], x: float) -> float | None:
    for (x0, y0), (x1, y1) in zip(pts, pts[1:], strict=False):
        if x0 <= x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return None


def angle_ratio(dataset: list[SeriesMetadata], prs: tuple[float, ...]) -> list[dict]:
    """Cd(30)/Cd(90) on Hay's parallel rig, Mc 0, Minf .31 (fig8 / fig6):
    measured (linear interpolation between marks) against off / on."""
    by = {s.path.stem: s for s in dataset if s.source.name == "hay1983"}
    m30 = [(p.x, p.y) for p in load_points(by["fig8_open_circles"])]
    m90 = [(p.x, p.y) for p in load_points(by["fig6_open_circles"])]
    out = []
    for pr in prs:
        a, b = _interp(m30, pr), _interp(m90, pr)
        rig = (HAY["d"], HAY["n"], HAY["pitch_over_d"], HAY["l_over_d"], HAY["Tt"])
        row = {"PR": pr, "measured": None if a is None or b is None else a / b}
        for inc in INCLINATION:
            row[inc] = model_cd(pr, 30.0, *rig, inc)[0] / model_cd(pr, 90.0, *rig, inc)[0]
        out.append(row)
    return out


def _summary(rows: list[Row], which: str) -> tuple[float, float]:
    e = [r.err(which) for r in rows]
    return sum(e) / len(e), sum(abs(x) for x in e) / len(e)


def _main() -> None:
    ds = load_dataset()
    print("Cd, coolant total -> exit static -- ACCURACY, CROSS-SOURCE (Idelchik vs Hay / Gritsch)")
    print(f"{'series':42s} {'alpha':>5s} {'n':>3s} {'Re':>13s} "
          f"{'off bias':>8s} {'off |e|':>7s} {'on bias':>8s} {'on |e|':>7s}")
    groups: list[tuple[str, list[Row]]] = [(s.label, rows_for(s)) for s in in_scope(ds)]
    carpet = carpet_rows()
    for key in sorted({r.label for r in carpet}):
        groups.append((key, [r for r in carpet if r.label == key]))
    for label, rows in groups:
        ob, oa = _summary(rows, "none")
        nb, na = _summary(rows, "idelchik_3_2")
        res = [r.re for r in rows]
        print(f"{label:42s} {rows[0].alpha:5.0f} {len(rows):3d} {min(res):6.0f}-{max(res):6.0f} "
              f"{ob:+8.1%} {oa:7.1%} {nb:+8.1%} {na:7.1%}")
        for r in rows:
            print(f"    PR {r.pr:.3f} meas {r.measured:.3f} off {r.off:.3f} on {r.on:.3f} "
                  f"in_range(on) {int(r.in_range_on)}")
    print("\nCd(30)/Cd(90), Hay parallel rig, Mc 0, Minf .31 (fig8 / fig6)")
    for row in angle_ratio(ds, (1.1, 1.2, 1.3, 1.5, 1.7, 1.9)):
        meas = "  -  " if row["measured"] is None else f"{row['measured']:.3f}"
        print(f"  PR {row['PR']:.1f} meas {meas} off {row['none']:.3f} on {row['idelchik_3_2']:.3f}")
    print("Cd(30)/Cd(90), Hay carpets Mc 0, Minf ~0: parallel 0.73/0.81 = "
          f"{0.73 / 0.81:.3f}, perpendicular 0.71/0.78 = {0.71 / 0.78:.3f}")
    print("\nOut of scope (not scored):")
    for label, why in out_of_scope(ds):
        print(f"  {label:44s} {why}")


if __name__ == "__main__":
    _main()
