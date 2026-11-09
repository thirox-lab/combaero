"""Score EffusionPlateElement's network eta_overall against Behrendt et al. (2008).

GT2008-51191: total (overall) cooling effectiveness of effusion liner samples
behind a kerosene combustor at 5 and 10 bar (DLR Cologne), eta_tot(x/d) on a
conducting Nimonic 75 wall. Data and record: validation/cooling/data/
behrendt2008/ and extractions/behrendt_2008_effusion_eta_tot.md.

WHAT A MISS MEANS: ACCURACY, CROSS-SOURCE. None of the element's closures
comes from this paper or lab -- the coolant side is Andrews 86-GT-225 (the
scoring set, `andrews_1986_effusion_internal`), the gas side the chamber's
surface correlation or the CR-2648 flat plate. Nothing here is fitted:
`gas_augmentation`, `internal_Nu_multiplier` and the chamber's
`Nu_multiplier` stay 1.0 for every scored number, and no inversion is ever
applied (docs/VALIDATION_POLICY.md).

THE NETWORK, NOT THE FORMULA. Unlike `effusion_overall_runner` (Andrews'
closure called directly), this builds the rig as a network and reads the
plate's own `wall_heat_transfer` diagnostics, so it scores what a user
gets:

    MassFlowBoundary(gas) -> Lossless -> MomentumChamberNode(A_duct) -> Lossless -> exit
    PressureBoundary(coolant) -> EffusionPlateElement -> chamber

BOTH GAS-SIDE BASELINES ARE SCORED, EACH ON ITS OWN ROW (reporting group),
never pooled and never ranked:

* `chamber` (the default): the chamber's own SmoothModel(gnielinski), fully
  developed, on Dh = sqrt(4 A/pi) = 319 mm.
* `flat_plate`: a turbulent boundary layer restarted at the rig's bleed,
  21 d upstream of the sample's leading edge (pp. 1803-1804), averaged over
  the scored window. Its Re_x (2-5e4) is below the CR-2786 data span
  (1.5e5), so those rows are out-of-domain BY THE KERNEL'S OWN FLAG.

THE SCORED WINDOW IS THE MEASURED ONE. The element is homogenised: it
predicts one eta for its panel. The data are a profile, and the plateau the
record defines is x/d 18-50, which all nine curves span (E20 stops at 50).
So the panel is that window -- 32 d long, starting 18 d into the sample --
and every vertex inside it is scored against the window's prediction;
vertices beyond x/d 50 belong to another window and are returned unscored.
That placement is geometry (where the measurement is), not a choice made on
the score: the whole-sample panel (x/d 0-143) is reported beside it by
`sensitivity`, and for the chamber baseline the window does not matter (the
coefficient has no x).

RIG BOUNDARY CONDITIONS, NOT MODEL CHOICES (all reported, none tuned):

* Duct area 80e3 mm^2 (ICAS 2006). The ICAS paper's y = 40 mm symmetry plane
  suggests an 80 mm deep duct, which at that area is 1 m wide: possibly a
  typo for 8e3 mm^2 (80 x 100 mm, Dh 89 mm). Only the CHAMBER baseline sees
  it (Dh enters its correlation); the flat plate takes the measured u_h,
  which the solved chamber velocity reproduces either way.
* Hot gas: no kerosene species, so lean n-heptane/air complete-combustion
  products at the ICAS operating point's 1750 K flame from 650 K preheat; T_h
  is the paper's NEAR-WALL gas at x/d -10 (Table 2), as conventions.py's
  eta_overall entry requires, via the chamber's main inlet.
* Coolant: dry air at T_c, supply total pressure p_h (1 + dp_rc).
* Wall: Nimonic 75, k at 750 K from the materials database.

GAS RADIATION (#1071) IS SCORED ON ITS OWN ROWS, beside each baseline without
it, never pooled: `radiation=GasRadiation` on the plate, Lefebvre & Ballal Eq.
8.7/8.11 at L = 1 (the authors' spectra showed no soot radiation, p. 1805).
Its inputs are rig boundary conditions too:

* Beam length 3.4 V/A of the duct cross-section (Eq. 8.10), 80 x 1000 mm
  with A_DUCT; the 80 x 100 mm reading in `sensitivity`.
* Wall emissivity 0.9: the hot face was coated "with a high emissivity
  colour" (p. 1803), no value printed; 0.7 (Nimonic, the spec's default) in
  `sensitivity`.
* Gas: T_h, the chamber's static P, the hot gas's composition (q from its
  elements) -- the spec's stream defaults.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import combaero as cb
from combaero.network import (
    EffusionPlateElement,
    GasRadiation,
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    PressureBoundary,
)
from combaero.network import NetworkSolver

from validation.cooling.schema import SeriesMetadata, load_points

SOURCE = "behrendt2008"
SET = "andrews_1986_effusion_internal"
BASELINES = ("chamber", "flat_plate")
RADIATION = (False, True)

# The rig (pp. 1803-1804 and ICAS 2006-152).
A_DUCT = 80.0e3 * 1e-6  # [m^2]
DH_ALT = 4.0 * (0.080 * 0.100) / (2.0 * (0.080 + 0.100))  # the 80 x 100 mm reading
SAMPLE_LENGTH_OVER_D = 143.0  # l/d, l = 100 mm
SAMPLE_WIDTH = 0.080  # b = 80 mm
BLEED_OVER_D = 21.0  # BL bleed at dx/d -21: where the hot boundary layer restarts
WINDOW_OVER_D = (18.0, 50.0)  # the record's plateau window
T_WALL_REF = 750.0  # [K] for the plate's k
T_PREHEAT = 650.0  # [K] ICAS 2006 combustor air
T_FLAME = 1750.0  # [K] ICAS 2006 global flame temperature
DUCT_DEPTH = 0.080  # [m] the y = 40 mm symmetry plane
EPS_WALL = 0.9  # "high emissivity colour" on the hot face, no value printed
EPS_WALL_ALT = 0.7  # Nimonic, GasRadiation's default


def beam_length(width: float = A_DUCT / DUCT_DEPTH) -> float:
    """3.4 V/A of the duct cross-section per unit length (Eq. 8.10) [m]."""
    return cb.mean_beam_length(DUCT_DEPTH * width, 2.0 * (DUCT_DEPTH + width))

_Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


@lru_cache(maxsize=1)
def hot_gas() -> tuple[tuple[float, ...], float]:
    """(Y, phi): lean n-heptane/air complete-combustion products at T_FLAME.

    A rig boundary condition. Kerosene (H/C ~1.9) has no species here;
    n-heptane (2.29) is the nearest liquid surrogate, and at phi ~0.47 the
    products' transport barely depends on which.
    """
    x_air = list(cb.species.dry_air())
    i_f = cb.species.species_index("NC7H16")
    lo, hi = 0.1, 1.0
    for _ in range(80):
        phi = 0.5 * (lo + hi)
        n_f = phi * 0.2095 / 11.0  # mol fuel per mol air: C7H16 + 11 O2
        xr = [x / (1.0 + n_f) for x in x_air]
        xr[i_f] += n_f / (1.0 + n_f)
        st = cb.complete_combustion(T_PREHEAT, xr, 5.0e5)
        if st.T > T_FLAME:
            hi = phi
        else:
            lo = phi
    return tuple(cb.mole_to_mass(list(st.X))), phi


@dataclass(frozen=True)
class Record:
    """One extracted vertex, evaluated (one per baseline)."""

    series: SeriesMetadata
    x: float  # x/d
    measured: float
    predicted: float | None
    extrapolated: bool
    reason: str | None = None
    group: str | None = None

    @property
    def rel_error(self) -> float | None:
        if self.predicted is None or self.measured == 0.0:
            return None
        return self.predicted / self.measured - 1.0

    @property
    def within_uncertainty(self) -> bool:
        band = self.series.uncertainty
        err = self.rel_error
        if band is None or err is None:
            return False
        return abs(err) <= band


def owns(series: SeriesMetadata) -> bool:
    """Every behrendt2008 series; the Fig. 8 heat fluxes come back unscored."""
    return series.source.name == SOURCE


def _group(baseline: str, radiation: bool = False) -> str:
    return f"eta_tot, {baseline} gas side" + (" + gas radiation" if radiation else "")


@lru_cache(maxsize=None)
def _solve(
    angle: float,
    d: float,
    t_over_d: float,
    sx_over_d: float,
    sz_over_d: float,
    p_h: float,
    T_h: float,
    u_h: float,
    T_c: float,
    dp: float,
    baseline: str,
    x_start_over_d: float,
    length_over_d: float,
    dh: float | None,
    radiation: tuple[float, float] | None = None,
) -> dict[str, float]:
    """`radiation` = (beam length, wall emissivity), or None for none."""
    y_gas, _ = hot_gas()
    y_gas = list(y_gas)
    x_gas = list(cb.mass_to_mole(y_gas))
    m_gas = cb.density(T_h, p_h, x_gas) * u_h * A_DUCT
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("gas_in", m_dot=m_gas, Tt=T_h, Y=y_gas))
    ex = PressureBoundary("exit")
    ex.Pt, ex.Tt, ex.Y = p_h, T_h, y_gas
    g.add_node(ex)
    cool = PressureBoundary("cool")
    cool.Pt, cool.Tt, cool.Y = p_h * (1.0 + dp), T_c, _Y_AIR
    g.add_node(cool)
    g.add_node(MomentumChamberNode("ch", area=A_DUCT, Dh=dh, main_inlet="feed"))
    g.add_element(LosslessConnectionElement("feed", "gas_in", "ch"))
    g.add_element(LosslessConnectionElement("tail", "ch", "exit"))
    kw: dict = {}
    if baseline == "flat_plate":
        kw = {"gas_baseline": "flat_plate", "gas_start_length": x_start_over_d * d}
    if radiation is not None:
        kw["radiation"] = GasRadiation(beam_length=radiation[0], emissivity_wall=radiation[1])
    g.add_element(
        EffusionPlateElement(
            "eff",
            "cool",
            "ch",
            hole_diameter=d,
            wall_thickness=t_over_d * d,
            pitch_x=sx_over_d * d,
            pitch_y=sz_over_d * d,
            panel_length=length_over_d * d,
            panel_width=SAMPLE_WIDTH,
            angle_deg=angle,
            wall_conductivity=cb.get_material_conductivity("nimonic75", T_WALL_REF),
            material="nimonic75",
            **kw,
        )
    )
    r = NetworkSolver(g).solve()
    if not r.get("__success__"):
        raise RuntimeError(f"behrendt2008 network did not solve: {r.get('__message__')}")
    return dict(r["__element_diag__"]["eff"])


def predict(
    series: SeriesMetadata,
    baseline: str,
    *,
    bleed_over_d: float = BLEED_OVER_D,
    window_over_d: tuple[float, float] = WINDOW_OVER_D,
    dh: float | None = None,
    radiation: bool = False,
    beam: float | None = None,
    eps_wall: float = EPS_WALL,
) -> dict[str, float]:
    """The plate's diagnostics for one series' operating point.

    `bleed_over_d` and `window_over_d` set the flat plate's x range: the
    boundary layer starts `bleed_over_d` upstream of the sample's leading
    edge and the panel covers the window. `radiation` adds gas radiation at
    `beam` (default `beam_length()`) and `eps_wall`. Defaults are the rig's.
    """
    g = series.geometry or {}
    lo, hi = window_over_d
    return _solve(
        float(g["angle_deg"]),
        float(g["d_mm"]) * 1e-3,
        float(g["t_over_d"]),
        float(g["sx_over_d"]),
        float(g["sz_over_d"]),
        float(g["p_h_bar"]) * 1e5,
        float(g["T_h_K"]),
        float(g["u_h_m_s"]),
        float(g["T_c_K"]),
        float(g["dp_rc_pct"]) * 1e-2,
        baseline,
        bleed_over_d + lo,
        hi - lo,
        dh,
        (beam if beam is not None else beam_length(), eps_wall) if radiation else None,
    )


def run_series(series: SeriesMetadata) -> list[Record]:
    points = load_points(series)
    if series.scores is None:
        return [
            Record(series, p.x, p.y, None, False, "not scored by any set") for p in points
        ]
    if series.y_axis != "eta_overall" or series.x_axis != "x_over_d":
        return [
            Record(series, p.x, p.y, None, False, f"y_axis={series.y_axis} not implemented")
            for p in points
        ]
    lo, hi = WINDOW_OVER_D
    out: list[Record] = []
    for baseline in BASELINES:
        for rad in RADIATION:
            d = predict(series, baseline, radiation=rad)
            extrapolated = bool(d.get("h_gas_unblown_extrapolated", 0.0))
            group = _group(baseline, rad)
            for p in points:
                if lo <= p.x <= hi:
                    out.append(
                        Record(series, p.x, p.y, d["eta_overall"], extrapolated, None, group)
                    )
                else:
                    out.append(
                        Record(
                            series,
                            p.x,
                            p.y,
                            None,
                            extrapolated,
                            "outside the scored window x/d 18-50",
                            group,
                        )
                    )
    return out


def run_all(dataset: list[SeriesMetadata]) -> list[Record]:
    out: list[Record] = []
    for series in dataset:
        out.extend(run_series(series))
    return out


def plateau(series: SeriesMetadata, window_over_d: tuple[float, float] = WINDOW_OVER_D) -> float:
    """Trapezoid mean of the measured eta over the window (the record's plateau)."""
    pts = [p for p in load_points(series) if window_over_d[0] <= p.x <= window_over_d[1]]
    area = sum(0.5 * (a.y + b.y) * (b.x - a.x) for a, b in zip(pts, pts[1:], strict=False))
    return area / (pts[-1].x - pts[0].x)


def scored(dataset: list[SeriesMetadata]) -> list[SeriesMetadata]:
    return [s for s in dataset if owns(s) and s.scores is not None]


# The record's one-parameter trends: (name, series stem from, series stem to).
TRENDS = (
    ("angle 20 -> 30 deg (TC1 dp 3% Tc 540 K)", "fig9_eta_tot_E20_dp3", "fig6_eta_tot_E30_Tc540"),
    ("T_c 470 -> 580 K (E30 TC1 dp 3%)", "fig6_eta_tot_E30_Tc470", "fig6_eta_tot_E30_Tc580"),
    ("dp 2 -> 4% (E20 TC1 Tc 540 K)", "fig9_eta_tot_E20_dp2", "fig9_eta_tot_E20_dp4"),
    ("p_h 5 -> 10 bar (E30 dp 3%)", "fig6_eta_tot_E30_Tc470", "fig10_eta_tot_E30_10bar_dp3"),
    ("p_h 5 -> 10 bar (E30 dp ~2%)", "fig10_eta_tot_E30_5bar_dp2", "fig10_eta_tot_E30_10bar_dp2"),
)


def trends(dataset: list[SeriesMetadata]) -> list[dict[str, float | str]]:
    """Measured plateau differences against each baseline's predicted ones."""
    by_stem = {s.path.stem: s for s in scored(dataset)}
    rows: list[dict[str, float | str]] = []
    for name, a, b in TRENDS:
        sa, sb = by_stem[a], by_stem[b]
        row: dict[str, float | str] = {"trend": name, "measured": plateau(sb) - plateau(sa)}
        for baseline in BASELINES:
            for rad in RADIATION:
                key = baseline + ("+rad" if rad else "")
                row[key] = (
                    predict(sb, baseline, radiation=rad)["eta_overall"]
                    - predict(sa, baseline, radiation=rad)["eta_overall"]
                )
        rows.append(row)
    return rows


def sensitivity(series: SeriesMetadata) -> dict[str, float]:
    """eta for the inputs that are boundary conditions, REPORTED, never chosen.

    * the flat plate's start: the bleed (default), half and twice as far, and
      a boundary layer from the sample's own leading edge;
    * the whole-sample panel (x/d 0-143) instead of the scored window;
    * the chamber at the 80 x 100 mm duct reading (Dh 89 mm), which the flat
      plate does not see;
    * gas radiation on the flat plate at the 80 x 100 mm beam length and at
      wall emissivity 0.7.
    """
    out = {
        "flat_plate": predict(series, "flat_plate")["eta_overall"],
        "flat_plate_bleed_half": predict(series, "flat_plate", bleed_over_d=0.5 * BLEED_OVER_D)[
            "eta_overall"
        ],
        "flat_plate_bleed_double": predict(
            series, "flat_plate", bleed_over_d=2.0 * BLEED_OVER_D
        )["eta_overall"],
        "flat_plate_from_sample_le": predict(series, "flat_plate", bleed_over_d=0.0)[
            "eta_overall"
        ],
        "flat_plate_whole_sample": predict(
            series, "flat_plate", window_over_d=(0.0, SAMPLE_LENGTH_OVER_D)
        )["eta_overall"],
        "chamber": predict(series, "chamber")["eta_overall"],
        "chamber_dh_alt": predict(series, "chamber", dh=DH_ALT)["eta_overall"],
        "flat_plate_dh_alt": predict(series, "flat_plate", dh=DH_ALT)["eta_overall"],
        "flat_plate_rad": predict(series, "flat_plate", radiation=True)["eta_overall"],
        "flat_plate_rad_beam_alt": predict(
            series, "flat_plate", radiation=True, beam=beam_length(0.100)
        )["eta_overall"],
        "flat_plate_rad_eps_alt": predict(
            series, "flat_plate", radiation=True, eps_wall=EPS_WALL_ALT
        )["eta_overall"],
    }
    return out


def _main() -> None:
    from validation.cooling.schema import load_dataset

    ds = load_dataset()
    print("eta_tot plateau x/d 18-50, gas_augmentation 1.0 -- ACCURACY, cross-source")
    print(f"{'series':34s} {'meas':>6s} {'chamber':>8s} {'err':>7s} {'flat':>7s} {'err':>7s} "
          f"{'ch+rad':>7s} {'err':>7s} {'fp+rad':>7s} {'err':>7s} {'q_rad':>6s} {'eps_g':>6s}")
    for s in scored(ds):
        m = plateau(s)
        dc, df = predict(s, "chamber"), predict(s, "flat_plate")
        rc, rf = predict(s, "chamber", radiation=True), predict(s, "flat_plate", radiation=True)
        print(
            f"{s.path.stem:34s} {m:6.3f} {dc['eta_overall']:8.3f} "
            f"{dc['eta_overall'] / m - 1:+7.1%} {df['eta_overall']:7.3f} "
            f"{df['eta_overall'] / m - 1:+7.1%} {rc['eta_overall']:7.3f} "
            f"{rc['eta_overall'] / m - 1:+7.1%} {rf['eta_overall']:7.3f} "
            f"{rf['eta_overall'] / m - 1:+7.1%} {rf['q_rad_hot'] * 1e-3:6.1f} {rf['eps_g']:6.3f}"
        )
    for row in trends(ds):
        print(
            f"{row['trend']:42s} meas {row['measured']:+.4f} "
            + " ".join(f"{k} {v:+.4f}" for k, v in row.items() if k not in ("trend", "measured"))
        )
    for s in scored(ds):
        print(s.path.stem, {k: round(v, 4) for k, v in sensitivity(s).items()})


if __name__ == "__main__":
    _main()
