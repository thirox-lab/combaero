"""Score a ribbed ChannelElement's pressure drop against measured rib friction (#1106).

`runner.py` scores `evaluate_rib`, the correlation. This scores the ELEMENT:
the drop a ribbed `ChannelElement` puts in its residual row, read back as the
passage-average Fanning friction the sources measure,

    fbar = dP / {4 (L/D_h) rho V^2 / 2}

(Han and Park 1988, p. 187). A correlation that is right can still be applied
wrongly by the element -- in the wrong friction convention, or for the wrong
number of ribbed walls -- and only a score through the element sees that.

THE MEASURED SIDE IS CONVERTED BY DEFINITION ONLY. Han's R series plot the
roughness function of the FOUR-SIDED friction f (Eq. 7); Han formed it from
his measured two-wall fbar by area weighting (Eq. 6), so the measured fbar is
recovered by running the same definitions backwards:

    f    = 2 / (R - 2.5 ln((2e/D)(2W/(W+H))) - 2.5)^2       (Eq. 7 inverted)
    Re   = e+ / ((e/D) sqrt(f/2))                          (e+ definition)
    fbar = (f W/H + f_s) / (W/H + 1),  f_s = 0.046 Re^-0.2 (Eq. 6 inverted)

R-versus-alpha series carry no e+, so their Re is ARBITRARY_RE (R is
Re-independent; only f_s moves, weakly). Taslim and Spring's series are fbar
already; Naik's f/f0 goes through the registered Petukhov conversion.

BASES ARE KEPT APART (docs/VALIDATION_POLICY.md). Fidelity is a set on its
own paper's data, n = 2; accuracy is everything else, labelled cross-source.
One-wall series (Taslim and Spring, open markers) are accuracy for the
area-weighted n = 1 rule, which no Han paper measured. Nothing here tunes:
`f_multiplier` stays 1.0.

    uv run python -m validation.cooling.rib_element_friction_runner
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass

import combaero as cb
from combaero import _solver_tools
from combaero.network.components import (
    ChannelElement,
    ConvectiveSurface,
    NetworkMixtureState,
    RibbedModel,
)
from validation.cooling import conventions, runner, schema, scorecard

X_AIR = cb.species.dry_air()
Y_AIR = cb.mole_to_mass(X_AIR)
T_RIG, P_RIG = 300.0, 1.5e5
LENGTH, DIAMETER = 1.0, 0.05

R_AXES = ("R", "R_normalised", "R_normalised_angled")
F_AXES = ("f_fanning_passage",)
RATIO_AXES = ("f_ratio_petukhov",)

# Taslim's one-side series are scored by the set their two-side twins use.
ONE_WALL_SET = "han_1988_orthogonal"


@dataclass
class Row:
    source: str
    series: str
    set_name: str
    basis: str
    n_walls: int
    re: float
    measured: float  # passage fbar, Fanning
    predicted: float  # element fbar, Fanning


def _geometry(series: schema.SeriesMetadata, rib_set) -> cb.RibGeometry:
    geom = runner._probe_geometry(rib_set)
    g = series.geometry or {}
    geom.e_D = float(g.get("e_D", geom.e_D))
    geom.p_e = float(g.get("p_e", geom.p_e))
    geom.W_H = float(g.get("W_H", geom.W_H))
    if series.alpha_deg is not None:
        geom.alpha_deg = float(series.alpha_deg)
    elif "alpha_deg" in g:
        geom.alpha_deg = float(g["alpha_deg"])
    return geom


def _f_s(re: float) -> float:
    return runner.SMOOTH_F_COEF * re**runner.SMOOTH_F_EXP


def _four_sided_from_R(r_abs: float, geom: cb.RibGeometry) -> float:
    group = (2.0 * geom.e_D) * (2.0 * geom.W_H / (geom.W_H + 1.0))
    bracket = r_abs - 2.5 * math.log(group) - 2.5
    return 2.0 / (bracket * bracket)


def element_fbar(rib_set, geom: cb.RibGeometry, n_walls: int, re: float, f_mult: float = 1.0) -> float:
    """The element's residual drop at Re, as a passage Fanning fbar."""
    el = ChannelElement(
        id="ch",
        from_node="A",
        to_node="B",
        length=LENGTH,
        diameter=DIAMETER,
        roughness=0.0,
        surface=ConvectiveSurface(
            area=0.1,
            model=RibbedModel(
                correlation_set=rib_set,
                e_D=geom.e_D,
                p_e=geom.p_e,
                alpha_deg=geom.alpha_deg,
                W_H=geom.W_H,
                n_ribbed_walls=n_walls,
            ),
        ),
    )
    el.surface.f_multiplier = f_mult
    rho, _, _ = _solver_tools.density_and_jacobians(T_RIG, P_RIG, X_AIR)
    mu = cb.complete_state(T_RIG, P_RIG, X_AIR).transport.mu
    area = el.area
    m_dot = re * area * mu / DIAMETER
    a = NetworkMixtureState(P=P_RIG, Pt=P_RIG, Tt=T_RIG, m_dot=m_dot, Y=Y_AIR)
    b = NetworkMixtureState(P=0.9 * P_RIG, Pt=0.9 * P_RIG, Tt=T_RIG, m_dot=m_dot, Y=Y_AIR)
    res, _ = el.residuals(a, b)
    dP = a.Pt - b.Pt - res[0]
    return dP * 2.0 * rho * area * area / (4.0 * (LENGTH / DIAMETER) * m_dot * m_dot)


def _rows_for(series: schema.SeriesMetadata, set_name: str, n_walls: int, f_mult: float) -> list[Row]:
    rib_set = runner.SETS[set_name]()
    geom = _geometry(series, rib_set)
    basis = scorecard.basis_of(series, set_name) if n_walls == 2 else "accuracy"
    to_fbar = None
    if series.y_axis in RATIO_AXES:
        resolution = conventions.resolve(series, set_name)
        if resolution.status != conventions.STATUS_CONVERTED or resolution.conversion is None:
            return []
        to_fbar = resolution.conversion.factor
    out: list[Row] = []
    for p in schema.load_points(series):
        if series.y_axis in R_AXES:
            if series.x_axis == "alpha_deg":
                geom.alpha_deg = p.x
            r_abs = p.y
            if series.y_axis != "R":
                r_abs = p.y / runner._normalised_R(rib_set, geom, 1.0)
            f4 = _four_sided_from_R(r_abs, geom)
            if series.x_axis == "e_plus":
                re = p.x / (geom.e_D * math.sqrt(f4 / 2.0))
            else:
                re = runner.ARBITRARY_RE
            measured = (f4 * geom.W_H + _f_s(re)) / (geom.W_H + 1.0)
        elif series.y_axis in F_AXES:
            re, measured = p.x, p.y
        else:
            re, measured = p.x, p.y * to_fbar(p.x)
        predicted = element_fbar(rib_set, geom, n_walls, re, f_mult)
        out.append(
            Row(series.source.name, series.path.name, set_name, basis, n_walls, re, measured, predicted)
        )
    return out


def run(f_mult: float = 1.0) -> list[Row]:
    rows: list[Row] = []
    for series in schema.load_dataset():
        if series.kind != "measured":
            continue
        if series.y_axis not in R_AXES + F_AXES + RATIO_AXES:
            continue
        if series.x_axis not in ("e_plus", "alpha_deg", "Re_Dh"):
            continue
        n_walls = int((series.geometry or {}).get("turbulated_walls", 2))
        set_name = series.scores or (ONE_WALL_SET if n_walls == 1 else None)
        if set_name is None:
            continue
        rows.extend(_rows_for(series, set_name, n_walls, f_mult))
    return rows


def _stats(rs: list[Row]) -> tuple[int, float, float]:
    errs = [r.predicted / r.measured - 1.0 for r in rs]
    return len(errs), statistics.fmean(errs), statistics.fmean(abs(e) for e in errs)


def render(rows: list[Row]) -> list[str]:
    lines = [
        "Ribbed ChannelElement friction, element fbar vs measured passage fbar (Fanning)",
        f"{'basis':<9} {'n':>1} {'source':<18} {'series':<52} {'pts':>4} {'bias':>8} {'MAE':>7}",
    ]
    keys = sorted({(r.basis, r.n_walls, r.source, r.series, r.set_name) for r in rows})
    for basis, n, src, ser, _ in keys:
        sub = [r for r in rows if (r.basis, r.n_walls, r.source, r.series) == (basis, n, src, ser)]
        k, bias, mae = _stats(sub)
        lines.append(f"{basis:<9} {n:>1} {src:<18} {ser:<52} {k:>4} {bias:>+8.1%} {mae:>7.1%}")
    lines.append("")
    for basis in ("fidelity", "accuracy"):
        for n in (1, 2, 4):
            sub = [r for r in rows if r.basis == basis and r.n_walls == n]
            if sub:
                k, bias, mae = _stats(sub)
                label = "cross-source" if basis == "accuracy" else "own data"
                lines.append(f"POOLED {basis} ({label}) n={n}: {k} pts, bias {bias:+.1%}, MAE {mae:.1%}")
    return lines


def main() -> None:
    print("\n".join(render(run())))


if __name__ == "__main__":
    main()
