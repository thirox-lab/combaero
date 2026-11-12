"""Redraw a digitised figure from the committed data.

A numeric check tells you a series is self-consistent. It does not tell
you the cloud has the shape the page shows, that a class sits where your
eye says it should, or that two series were swapped in a way the
arithmetic happens to tolerate. Putting the data back on axes and holding
it next to the scan catches that in a second.

This reads the SAME metadata the verifier and scorecard use -- axis type,
tick span and plot limits all come from each series' figure card -- so a
plot cannot drift from what the checks believe.

    uv run python -m validation.cooling.plot --figure 4.46
    uv run python -m validation.cooling.plot --list

Needs the `examples` extra for matplotlib:  uv pip install -e ".[examples]"
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from validation.cooling.schema import SeriesMetadata, load_dataset, load_points

# Distinct markers so ten overlapping symbol classes stay separable -- the
# same problem the source figure has, which is why classes are hard to
# attribute in the first place.
MARKERS = ["o", "s", "^", "v", "D", "<", ">", "P", "X", "*", "h", "p"]

PANEL_ORDER = {"upper": 0, "single": 1, "lower": 2}

# Redrawn figures are committed. They are namespaced BY SOURCE, as data
# and cards are: a figure number alone is not unique -- two books can each
# have a figure 4.46 -- and a collision would silently overwrite one
# source's plot with another's.
PLOT_DIR = Path(__file__).parent / "plots"


def plot_path(source: str, figure: str) -> Path:
    """Where a source's redrawn figure lives."""
    return PLOT_DIR / source / f"fig{figure}_redrawn.png"


def resolve_figure(dataset, wanted: str) -> tuple[str, str]:
    """Resolve a figure reference to (source, figure).

    Accepts a bare figure number while it is unambiguous, and the
    qualified `source/figure` form always. A bare number that matches more
    than one source is an error naming the alternatives rather than a
    silent pick.
    """
    if "/" in wanted:
        source, figure = wanted.split("/", 1)
        hits = [s for s in dataset if s.source.name == source and s.figure == figure]
        if not hits:
            raise SystemExit(f"no series for {wanted}")
        return source, figure

    sources = sorted({s.source.name for s in dataset if s.figure == wanted})
    if not sources:
        raise SystemExit(f"no series for figure {wanted}")
    if len(sources) > 1:
        alts = ", ".join(f"{src}/{wanted}" for src in sources)
        raise SystemExit(
            f"figure {wanted} is ambiguous across sources; use one of: {alts}"
        )
    return sources[0], wanted

# What the figures actually print on their axes, rather than the terse
# keys the metadata uses.
AXIS_NAMES = {
    "e_plus": "e+ = (e/D) Re (f/2)^0.5",
    "G": "G",
    "G_bar": "G_bar",
    "R_normalised": "R / (P/e/10)^0.35",
    "f_ratio": "f / f_0",
    "Nu_ratio": "Nu_r / Nu_0",
    # Han 2013 Figs. 4.88/4.90 (Chandra et al.): the book defines neither
    # baseline, so these are not "Nu_ratio" and carry no convention.
    "Nu_ratio_smooth_wall": "Nu_s / Nu_0, smooth walls (baseline not printed)",
    "Nu_ratio_ribbed_wall": "Nu_r / Nu_0, ribbed walls (baseline not printed)",
    # Kays and Perkins, NOT Dittus-Boelter: 0.873 of it at Pr 0.72.
    "Nu_ratio_KP": "Nu / Nu_inf, Nu_inf = 0.0176 Re^0.8 (Kays & Perkins)",
    "Nu_ratio_DB": "Nu / Nu0, Nu0 = 0.023 Re^0.8 Pr^0.4 (Dittus-Boelter)",
    # Naik et al. (2014): Fanning over Petukhov's Fanning, Eqs. 4-5.
    "f_ratio_petukhov": "f / f0, Fanning, f0 = 2 (2.236 ln Re - 4.639)^-2",
    # Naik et al. (2014) Fig. 18 as printed -- NOT its own Eq. 6 for the
    # paper's own points (3 Nu/f); a cross-check, never scored.
    "performance_index_printed": "performance index PI as plotted",
    "G_coolant": "G [kg/s/m2], coolant per unit plate area",
    # The two h conventions spelled out, because a bare "h" on an axis is
    # exactly how 86-GT-225's h_m and 88-GT-290's h get confused.
    "h_internal": "h on plate area A = X2 - piD2/4 [W/m2K]",
    "h_hole_length": "h_m on hole area piDL [W/m2K]",
    "h_imp_eff_wall": "imp./effusion wall h on A = X2 - piD2/4 [W/m2K]",
    "h_impingement_target": "impingement target-plate h [W/m2K]",
    "eta_overall": "overall cooling effectiveness (Tg-Tw)/(Tg-Tc)",
    "x_over_d": "axial coordinate x/d",
    "dp_rc_percent": "cooling air pressure drop dp_r,c [%]",
    # Behrendt et al. (2008) Fig. 8: biased by the sensor, never scored.
    "q_hot_side_W_cm2": "hot-gas-side heat flux q'' [W/cm2]",
    "pressure_ratio": "pressure ratio r = p / p0",
    # Hay et al. (1983): the INVERSE sense of Jobson's r, and on the coolant
    # total / mainstream static pair -- hence its own name.
    "pc0_over_pinf": "pressure ratio pc+ / p_inf (coolant total / mainstream static)",
    "uj_over_uinf": "velocity ratio u_j / u_inf",
    "Ma_c": "coolant-channel (internal crossflow) Mach number Ma_c",
    "K_mass_flow": "mass-flow coefficient K = m / (A sqrt(p0 rho0))",
}


def _curve(spec: dict, xs: list[float]) -> list[float]:
    """Evaluate a printed_curve spec across xs."""
    if "constant" in spec:
        return [float(spec["constant"])] * len(xs)
    if "power_law" in spec:
        pl = spec["power_law"]
        return [float(pl["C"]) * x ** float(pl["n"]) for x in xs]
    if "polynomial" in spec:
        poly = spec["polynomial"]
        scale = float(poly.get("variable_scale", 1.0))
        return [
            sum(float(c) * (x / scale) ** i for i, c in enumerate(poly["coeffs"]))
            for x in xs
        ]
    raise ValueError(f"unrecognised printed_curve: {sorted(spec)}")


def _label(s: SeriesMetadata) -> str:
    geom = s.geometry or {}
    parts = []
    if "e_D" in geom:
        parts.append(f"e/D {geom['e_D']}")
    if "p_e" in geom:
        parts.append(f"P/e {geom['p_e']}")
    if "W_H" in geom:
        parts.append(f"W/H {geom['W_H']}")
    base = "  ".join(parts) if parts else s.series[:46]
    if s.class_confidence == "disputed":
        base += "  [DISPUTED]"
    return base


def plot_figure(
    figure: str,
    out: Path,
    dataset=None,
    show_frames: bool = False,
    source: str | None = None,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dataset = dataset or load_dataset()
    chosen = [
        s
        for s in dataset
        if s.figure == figure and (source is None or s.source.name == source)
    ]
    if not chosen:
        raise SystemExit(f"no series for figure {figure}")

    panels: dict[str, list[SeriesMetadata]] = defaultdict(list)
    for s in chosen:
        panels[s.panel or "single"].append(s)
    names = sorted(panels, key=lambda p: PANEL_ORDER.get(p, 1))

    fig, axes = plt.subplots(
        len(names), 1, figsize=(9.5, 5.2 * len(names)), squeeze=False
    )

    for ax, panel in zip(axes[:, 0], names, strict=False):
        marker_i = 0
        card0 = None
        frame_y: list[float] = []
        for s in sorted(panels[panel], key=lambda s: (s.kind, s.path.name)):
            card = s.verification or {}
            card0 = card0 or card
            pts = load_points(s)
            xs = [p.x for p in pts]
            ys = [p.y for p in pts]

            if s.kind == "frame":
                # Off by default. A frame sits outside the labelled ticks,
                # so drawing it forces the axis open and squashes the data
                # into a strip -- and its job (measuring distortion) is
                # already done numerically by the verifier.
                if not show_frames:
                    continue
                ax.plot(xs, ys, ":", color="0.45", lw=1.2,
                        label="panel frame (distortion standard)")
                frame_y.extend(ys)
                continue

            if s.kind == "correlation":
                ax.plot(xs, ys, "--", lw=1.4, alpha=0.9,
                        label=f"drawn: {s.series[:44]}")
                printed = card.get("printed_curve")
                if printed:
                    grid = sorted(xs)
                    ax.plot(grid, _curve(printed, grid), "-", lw=1.1,
                            color="k", alpha=0.55, label="printed equation")
                continue

            ax.plot(xs, ys, MARKERS[marker_i % len(MARKERS)], ms=6,
                    mfc="none", mew=1.3, label=_label(s))
            marker_i += 1

        card0 = card0 or {}
        ax.set_xscale(card0.get("x_axis_type", "log"))
        ax.set_yscale(card0.get("y_axis_type", "log"))
        for setter, key in ((ax.set_xlim, "plot_x_limits"),
                            (ax.set_ylim, "plot_y_limits")):
            lim = card0.get(key)
            if lim:
                lo, hi = lim
                if key == "plot_y_limits" and frame_y:
                    lo = min(lo, min(frame_y) * 0.97)
                    hi = max(hi, max(frame_y) * 1.03)
                setter(lo, hi)
        real = [s for s in panels[panel] if s.kind != "frame"]
        sample = (real or panels[panel])[0]
        ax.set_xlabel(AXIS_NAMES.get(sample.x_axis, sample.x_axis))
        ylabels = {s.y_axis for s in real}
        ax.set_ylabel(
            " or ".join(AXIS_NAMES.get(y, y) for y in sorted(ylabels))
            if ylabels
            else sample.y_axis
        )
        ax.set_title(
            f"{chosen[0].source.name}  figure {figure} -- {panel} panel",
            fontsize=10,
        )
        ax.grid(True, which="both", alpha=0.25, lw=0.5)
        handles, labels = ax.get_legend_handles_labels()
        seen: dict[str, object] = {}
        for h, lab in zip(handles, labels, strict=False):
            seen.setdefault(lab, h)
        ax.legend(seen.values(), seen.keys(), fontsize=6.5, loc="best",
                  ncol=2, framealpha=0.9)

    # Source citation travels with the image, always. Wrapped, because a
    # citation clipped at the page edge is not a citation.
    import textwrap

    lines = ["Redrawn from digitised measurements."]
    lines += textwrap.wrap(f"Source: {chosen[0].source.citation}", 110)
    if chosen[0].after:
        lines += textwrap.wrap(f"After: {chosen[0].after}", 110)
    fig.text(0.5, 0.004, "\n".join(lines), ha="center", va="bottom",
             fontsize=6, color="0.35", linespacing=1.5)

    fig.tight_layout(rect=(0, 0.013 * len(lines) + 0.01, 1, 1))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--figure",
        help='which figure: "4.46", or "han2012/4.46" when ambiguous',
    )
    ap.add_argument("--out", type=Path, help="output image path")
    ap.add_argument("--list", action="store_true", help="list known figures")
    ap.add_argument("--all", action="store_true",
                    help="redraw every figure in the dataset")
    ap.add_argument("--frames", action="store_true",
                    help="draw panel frames too (opens the axes wide)")
    args = ap.parse_args()

    dataset = load_dataset()

    if args.all:
        pairs = sorted({(s.source.name, s.figure) for s in dataset if s.figure})
        for source, figure in pairs:
            out = plot_path(source, figure)
            print(
                f"wrote {plot_figure(figure, out, dataset, args.frames, source)}"
            )
        return

    if args.list or not args.figure:
        seen: dict[tuple[str, str], set] = defaultdict(set)
        for s in dataset:
            if s.figure:
                seen[(s.source.name, s.figure)].add(s.panel)
        for source, figure in sorted(seen):
            n = sum(
                1
                for s in dataset
                if s.figure == figure and s.source.name == source
            )
            ref = f"{source}/{figure}"
            print(f"  {ref:<20} {n:>3} series, panels: "
                  f"{', '.join(sorted(seen[(source, figure)]))}")
        return

    source, figure = resolve_figure(dataset, args.figure)
    out = args.out or plot_path(source, figure)
    print(f"wrote {plot_figure(figure, out, dataset, args.frames, source)}")


if __name__ == "__main__":
    main()
