"""Evaluate an implemented correlation set against the digitised series.

The runner never computes G from a set's coefficients itself. It calls
cb.evaluate_rib and bisects the Reynolds number until the chain's own e+
lands on the digitised abscissa, so the whole chain -- friction factor,
e+, then G -- is what gets scored. Reimplementing the formula here would
score the model against a second copy of itself, which is the failure
mode #333 exists to prevent.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import combaero as cb

from validation.cooling.conventions import STATUS_CONVERTED, STATUS_DIRECT
from validation.cooling.schema import Point, SeriesMetadata, load_points

# Han states G_bar (the ribbed/smooth average) as 1.2 G. Confirmed as item
# 10 of validation/cooling/extractions/han_ribbed.md, from the label on
# Figure 4.47, and independently by Figure 4.46's printed pair 3.7 and 4.5
# (ratio 1.216).
#
# IT IS A FIT, AND HAN MADE IT A CONSTANT. That is the point: the ratio
# is measurably not constant, and 1.2 is the single number Han settled on
# anyway. Pooling every closed-form measurement available -- 91-GT-3's
# table 3, Lau's table 2, and CR-3837's per-run Nu(R)/Nu(AV) split, 16
# configurations across three rigs at e+ = 300:
#
#     min 1.096   max 1.413   mean 1.284   spread 29%
#
# and it drifts with e+ within a configuration too (91-GT-3's 90 deg rib
# runs 1.29 down to 1.13 over e+ 150-1000).
#
# A generalised replacement was considered and REJECTED on the evidence:
# 69% of that variance is BETWEEN RIGS, not within them (per-source means
# 1.160 / 1.323 / 1.327, internal spreads 12% / 3% / 14%). No angle or
# shape term can reach a rig offset, so a "better" correlation would be
# fitting rig identity.
#
# ALL THREE RIGS ARE TEXAS A&M, and that cuts both ways. 91-GT-3 and
# CR-3837 are Han's own laboratory; lau1990 is Lau's, in the same
# department, and Lau sat on the committees for the rest. So:
#
#   * the "no correlation" conclusion gets STRONGER. If three rigs within
#     one laboratory tradition differ by 69% of the variance between
#     them, the offset is rig-level -- geometry, instrumentation, data
#     reduction -- and not an artefact of comparing institutions. There
#     is even less for an angle or shape term to reach.
#   * the "MAE 8.3%" gets WEAKER. This is not a sample of the world's
#     rigs, so the 29% spread is a LOWER BOUND on the true one and Han's
#     1.2 may score worse against a non-TAMU rig than it does here.
#
# Recorded because "16 configurations across three rigs" reads as more
# independent than it is. See #403 item 2 on the same monoculture.
#
# Against the pooled population Han's 1.2 scores
# bias -6.1%, MAE 8.3% -- a defensible compromise for one number, and
# close to what a least-squares fit over the same population would give.
#
# So it is applied AS PUBLISHED, at every rib angle. The error that
# produces is ACCURACY -- the model's limitation, reported and labelled --
# not a harness defect. See docs/VALIDATION_POLICY.md. An earlier version
# refused anything off 90 deg; that withheld a number Han does publish.
G_BAR_OVER_G = 1.2

# Smooth-channel references for the performance-curve view (Nu_ratio
# against f_ratio, figure 4.53), as printed by Han, Zhang and Lee (1991)
# JHT 113, 590, Eqs. (2) and (4) -- the paper whose continuous-rib data
# figure 4.53 reprints for comparison (extractions/han_zhang_lee_1991_jht.md):
#
#     f_0  = 0.046 Re^-0.2               Blasius form, smooth circular tube
#     Nu_0 = 0.023 Re^0.8 Pr^0.4         McAdams / Dittus-Boelter
#
# Superseded CR-4015's four-sided 0.079 Re^-0.25 (Han's lab, 1986), which
# #440 used before this paper was on disk; the two differ by under 3% over
# Re 15k-80k and the score by 0.3%. Written out here rather than calling
# nusselt_dittus_boelter, which warns below Re 10,000 on every bisection
# probe. The same f_0 is the f_s in Han's four-sided decomposition, so it
# also converts a measured passage friction (_channel_fbar).
SMOOTH_F_COEF, SMOOTH_F_EXP = 0.046, -0.2
SMOOTH_NU_COEF, SMOOTH_NU_RE_EXP, SMOOTH_NU_PR_EXP = 0.023, 0.8, 0.4

# Bracket for the Re bisection. Deliberately far wider than any set's
# stated validity: a series may sit outside it, and reporting that as
# extrapolated is the point rather than something to avoid.
# What the Re_Dh path predicts: passage friction (Taslim & Spring), the same
# friction as a ratio over Petukhov (Naik et al. 2014), and the ribbed-wall
# Nu ratio on a smooth-tube reference -- Kays and Perkins (CR-4396, #452) or
# Dittus-Boelter itself (Naik et al. 2014).
RE_PATH_Y_AXES = ("f_fanning_passage", "f_ratio_petukhov", "Nu_ratio_KP", "Nu_ratio_DB")
# The Nu-ratio axes, each with the ONE resolution status its name allows: a
# Kays-Perkins axis must convert, a Dittus-Boelter axis must not. A series
# whose axis and declared convention disagree is refused, so a KP series
# mislabelled DB cannot be scored on the raw basis (and vice versa).
NU_RATIO_AXIS_STATUS = {"Nu_ratio_KP": STATUS_CONVERTED, "Nu_ratio_DB": STATUS_DIRECT}

RE_LO = 1.0e3
RE_HI = 1.0e8
BISECT_STEPS = 200


@dataclass(frozen=True)
class Record:
    """One digitised point, evaluated."""

    series: SeriesMetadata
    x: float  # e+
    measured: float
    predicted: float | None  # None when the set has no binding
    extrapolated: bool
    re_used: float | None
    reason: str | None = None  # why there is no prediction
    # The set that made `predicted`. A series may be scored by several
    # (`also_scores`); each set's records stay separate all the way to the
    # report, so nothing averages or picks between them.
    scored_by: str | None = None

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


SETS = {
    "han_1988_orthogonal": cb.han_1988_orthogonal,
    "rallabandi_2009_high_re": cb.rallabandi_2009_high_re,
    "han_park_1988_angled": cb.han_park_1988_angled,
    "han_1989_narrow_channel": cb.han_1989_narrow_channel,
}

# Han, Zhang and Lee (1991) Table 2: one set per tested configuration, keyed
# by the name the C++ factory gives it (han_zhang_lee_1991_60vee, ...).
HZL_1991_CONFIGS = (
    ("90", "Transverse", 90.0),
    ("60par", "Parallel", 60.0), ("60crs", "Crossed", 60.0),
    ("60vee", "V", 60.0), ("60lam", "Lambda", 60.0),
    ("45par", "Parallel", 45.0), ("45crs", "Crossed", 45.0),
    ("45vee", "V", 45.0), ("45lam", "Lambda", 45.0),
)
for _tag, _shape, _alpha in HZL_1991_CONFIGS:
    SETS[f"han_zhang_lee_1991_{_tag}"] = (
        lambda shape=_shape, alpha=_alpha: cb.han_zhang_lee_1991(
            getattr(cb.RibShape, shape), alpha
        )
    )

# Series metadata names a rib shape as a lowercase string; v_broken has no
# RibShape because no implemented set models broken ribs.
SERIES_RIB_SHAPE = {
    "transverse": "Transverse",
    "parallel": "Parallel",
    "crossed": "Crossed",
    "v": "V",
    "lambda": "Lambda",
}


def _binding_reason(
    rib_set: "cb.RibCorrelationSet", series: SeriesMetadata
) -> str | None:
    """Why this set cannot answer for this series, or None if it can.

    A set outside its NUMERIC validity is still answering the question --
    that is an extrapolation, and the scorecard reports it as one. A set
    asked about a CONFIGURATION it has no binding for is not: scoring
    han_1988_orthogonal, whose valid_alpha is [90, 90] and whose G carries
    no alpha term, against 45 deg rib data produced a 26.7% bias that reads
    as model error when it is really the wrong correlation entirely. That
    is the shape of misleading metric this harness exists to prevent, so
    the configuration mismatch is refused rather than scored.
    """
    set_shape = rib_set.shape
    if set_shape != cb.RibShape.Unspecified:
        declared = series.rib_shape
        if declared is None:
            # Transverse and parallel sets predate shape metadata and score
            # series that never recorded one; a crossed, V or Lambda set must
            # never assume it.
            if set_shape not in (cb.RibShape.Transverse, cb.RibShape.Parallel):
                return (
                    f"series records no rib shape; {rib_set.name} binds "
                    f"{set_shape.name} ribs only"
                )
        elif SERIES_RIB_SHAPE.get(declared) != set_shape.name:
            return f"{declared} ribs; {rib_set.name} binds {set_shape.name} ribs only"

    alpha = series.alpha_deg
    rng = rib_set.valid_alpha
    if alpha is not None and rng.hi >= rng.lo:
        if not (rng.lo <= alpha <= rng.hi):
            return (
                f"{alpha:g} deg ribs; {rib_set.name} binds "
                f"{rng.lo:g}-{rng.hi:g} deg only"
            )
    return None


def _probe_geometry(rib_set: "cb.RibCorrelationSet") -> cb.RibGeometry:
    """A geometry at the centre of the set's stated validity.

    Only legitimate when the set's G carries no geometry dependence -- see
    _assert_geometry_free. The probe then fixes only WHICH Reynolds number
    reaches a given e+, not the G reported there.
    """
    geom = cb.RibGeometry()
    geom.e_D = _mid(rib_set.valid_eD, 0.0625)
    geom.p_e = _mid(rib_set.valid_pe, 10.0)
    geom.W_H = _mid(rib_set.valid_WH, 1.0)
    geom.alpha_deg = 90.0
    return geom


def _mid(rng: "cb.RibRange", fallback: float) -> float:
    # Same rule as RibRange::bounded(): hi == lo is one tested value.
    return 0.5 * (rng.lo + rng.hi) if rng.hi > 0.0 and rng.hi >= rng.lo else fallback


def _gbar_reason(series: SeriesMetadata) -> str | None:
    """Why a G_bar series cannot be converted, or None if it can.

    Nothing is refused on rib angle any more. `G_bar = 1.2 G` is what Han
    publishes, so applying it at every angle is the faithful
    implementation; what that costs is measured and reported as accuracy
    rather than hidden by declining to answer. See G_BAR_OVER_G.

    Kept as a hook because a future G_bar series may be unconvertible for
    a reason that is NOT rib angle -- a different averaging basis, say.
    Lau's Gbar is a four-wall average like Han's (#392), but another
    source's need not be.
    """
    del series  # no refusal applies today; see the docstring
    return None


def _binds_geometry(rib_set: "cb.RibCorrelationSet") -> bool:
    """True when the set's G actually depends on the rig geometry.

    A shape switch (SquareVsRectangular, NarrowChannelAlphaSwitch) binds G
    to geometry through dedicated fields, not through G_eD/G_pe/G_WH/G_alpha
    -- han_park_1988_angled happens to also set a nonzero G_WH exponent, so
    the term check alone has caught it so far, but that was never the real
    reason. Checked explicitly so a future shape-only set doesn't slip
    through silently.
    """
    return rib_set.G_shape_model != cb.GShapeModel.Fixed or any(
        t.exponent != 0.0
        for t in (rib_set.G_eD, rib_set.G_pe, rib_set.G_WH, rib_set.G_alpha)
    )


def _assert_geometry_free(rib_set: "cb.RibCorrelationSet") -> None:
    """Refuse to score a geometry-dependent set against geometry-less data.

    Applies only to series with NO recorded rig geometry -- the Han and
    Zhang (1992) figures, whose rig the extracted text never states. For a
    set whose G has zero exponents on e/D, P/e, W/H and alpha this costs
    nothing. For any other set it would mean substituting a default and
    reporting the result as a measurement, so the runner stops instead.

    Figure 4.46 states a geometry per class in its legend, so those series
    are exempt and reach the disputed-label check above instead.
    """
    terms = {
        "e/D": rib_set.G_eD,
        "P/e": rib_set.G_pe,
        "W/H": rib_set.G_WH,
        "alpha": rib_set.G_alpha,
    }
    bound = [name for name, t in terms.items() if t.exponent != 0.0]
    if rib_set.G_shape_model != cb.GShapeModel.Fixed:
        # Same reasoning as _binds_geometry: a shape switch binds G to
        # geometry through its own dedicated fields, not through the terms
        # above.
        bound.append("G_shape_model")
    if bound:
        raise ValueError(
            f"{rib_set.name} binds G to {', '.join(bound)}, but the digitised "
            "series carries no rig geometry. Supply the geometry in "
            "metadata.yaml from the primary paper before scoring this set."
        )


def _normalised_R(
    rib_set: "cb.RibCorrelationSet", geom: cb.RibGeometry, raw_R: float
) -> float:
    """Raw R divided back down to what a figure's y-axis actually plots.

    Both R shapes in this project print their correlation as R divided by
    its own normalisers (Han 1988: R/(p/e/10)^0.35; Han and Park 1988:
    R/[(p/e/10)^0.35 (W/H)^m]) -- never raw R. Recomputing the divisor from
    the SAME exposed fields evaluate_rib used keeps this from silently
    drifting out of step with a future change to either set's shape.
    """
    pe_term = rib_set.R_pe
    pe_factor = (
        (geom.p_e / pe_term.reference) ** pe_term.exponent
        if pe_term.reference and pe_term.exponent
        else 1.0
    )
    WH_factor = 1.0
    if rib_set.R_alpha_shape == cb.RAlphaShape.QuadraticAlpha:
        at_90 = abs(geom.alpha_deg - 90.0) < 1e-9
        m = (
            rib_set.R_quad_WH_exponent_at_90
            if at_90
            else rib_set.R_quad_WH_exponent_off_90
        )
        if m:
            W_H = geom.W_H
            if rib_set.R_quad_WH_cap > 0.0:
                W_H = min(W_H, rib_set.R_quad_WH_cap)
            WH_factor = W_H**m
    divisor = pe_factor * WH_factor
    return raw_R / divisor if divisor else raw_R


def _at_eplus(
    rib_set: "cb.RibCorrelationSet", geom: cb.RibGeometry, target: float
) -> tuple[float, float, float, bool, float, float] | None:
    """G, absolute R, normalised R and G_bar at a target e+, via the real chain.

    The lower panel of Figure 4.46 plots R/(P/e/10)^0.35, which is what
    the set's C_R is defined as, so the normalised value is returned
    alongside G rather than being recomputed by the caller. Absolute R is
    returned too: sources that TABULATE the roughness functions print R
    itself, not the figure's normalised ordinate.
    """
    lo, hi = RE_LO, RE_HI
    if cb.evaluate_rib(rib_set, geom, lo).e_plus > target:
        return None
    if cb.evaluate_rib(rib_set, geom, hi).e_plus < target:
        return None
    for _ in range(BISECT_STEPS):
        mid = 0.5 * (lo + hi)
        if cb.evaluate_rib(rib_set, geom, mid).e_plus < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-9 * max(1.0, hi):
            break
    re = 0.5 * (lo + hi)
    res = cb.evaluate_rib(rib_set, geom, re)
    pe_term = rib_set.R_pe
    norm = (geom.p_e / pe_term.reference) ** pe_term.exponent if pe_term.reference else 1.0
    r_norm = res.R / norm if norm else res.R
    # The source's own printed G_bar when the set carries one; Han's 1.2
    # otherwise (see G_BAR_OVER_G).
    g_bar = res.G_bar if res.has_G_bar else res.G * G_BAR_OVER_G
    return res.G, res.R, r_norm, res.extrapolated, re, g_bar


def _r_at_alpha(
    rib_set: "cb.RibCorrelationSet", geom: cb.RibGeometry, alpha_deg: float
) -> float:
    """R at a given rib angle. No Re bisection needed.

    Every set in this project has R carrying no e+ term (see
    rib_correlation.h's header comment: "R CARRIES NO e+ TERM by
    construction"), so R does not depend on Reynolds number at all. Any Re
    gives the same R -- ARBITRARY_RE below is not a choice that matters,
    only a value evaluate_rib needs to run.
    """
    geom.alpha_deg = alpha_deg
    return cb.evaluate_rib(rib_set, geom, ARBITRARY_RE).R


ARBITRARY_RE = 30000.0


def _channel_fbar(f_four_sided: float, w_h: float, re: float) -> tuple[float, float]:
    """(fbar, f_s): the MEASURED channel friction from the four-sided one.

    Han measures `fbar`, the passage-average Fanning friction of a channel
    with two opposite ribbed walls, and converts it to the four-sided
    equivalent `f = fbar + (H/W)(fbar - f_s)` before forming R (Han 1988;
    Han, Zhang and Lee 1991 Eq. 6 at W/H = 1). evaluate_rib returns that
    four-sided `f`, so comparing it with a measured passage friction needs
    the inverse, `fbar = (f W/H + f_s) / (W/H + 1)`, with the same `f_s` the
    decomposition used.
    """
    f_s = SMOOTH_F_COEF * re**SMOOTH_F_EXP
    return (f_four_sided * w_h + f_s) / (w_h + 1.0), f_s


def _performance_ratios(
    rib_set: "cb.RibCorrelationSet", geom: cb.RibGeometry, re: float
) -> tuple[float, float, bool]:
    """(f_ratio, Nu_ratio, extrapolated) at one Reynolds number.

    f_ratio is CHANNEL-AVERAGE friction over smooth: Han's measured `fbar`
    for two opposite ribbed walls, not the four-sided equivalent `f` that
    evaluate_rib returns. Inverting Han's `f = fbar + (H/W)(fbar - f_s)`:
    `fbar = (f W/H + f_s) / (W/H + 1)`. Taking `f/f_s` instead is a 1.8x
    error at figure 4.53's geometry -- the same trap CR-3837's `Rbar` set
    (han_park_lei_1984_cr3837.md, "The conversion").

    Nu_ratio is the RIBBED-side Nusselt number over Dittus-Boelter, which is
    what St_r already is.
    """
    res = cb.evaluate_rib(rib_set, geom, re)
    fbar, f_s = _channel_fbar(res.f, geom.W_H, re)
    pr = rib_set.valid_Pr
    nu_r = res.St_r * re * pr
    nu_s = SMOOTH_NU_COEF * re**SMOOTH_NU_RE_EXP * pr**SMOOTH_NU_PR_EXP
    return fbar / f_s, nu_r / nu_s, res.extrapolated


def _at_f_ratio(
    rib_set: "cb.RibCorrelationSet", geom: cb.RibGeometry, target: float
) -> tuple[float, bool, float] | None:
    """Nu_ratio where the set's own f_ratio reaches `target`, via Re.

    The performance curve plots both ratios against each other with Re
    eliminated, and the digitised points carry no Re. So Re is recovered
    the way the e+ path recovers it from e+: bisect the real chain until
    its f_ratio lands on the abscissa. A friction error therefore shifts
    WHERE on the curve a point is compared, which is how a reader of the
    figure compares curves too. f_ratio rises monotonically with Re for
    this family (f is Re-independent, f_s falls).
    """
    lo, hi = RE_LO, RE_HI
    if not (
        _performance_ratios(rib_set, geom, lo)[0]
        < target
        < _performance_ratios(rib_set, geom, hi)[0]
    ):
        return None
    for _ in range(BISECT_STEPS):
        mid = (lo * hi) ** 0.5
        if _performance_ratios(rib_set, geom, mid)[0] < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-9 * hi:
            break
    re = (lo * hi) ** 0.5
    _, nu_ratio, extrapolated = _performance_ratios(rib_set, geom, re)
    return nu_ratio, extrapolated, re


def run_series(
    series: SeriesMetadata,
    points: list[Point] | None = None,
    set_name: str | None = None,
) -> list[Record]:
    """Evaluate one series. Unscored series yield records with no prediction.

    `points` overrides the series' own CSV. Its only use is asking the set
    what it predicts at abscissae the digitisation could not pick -- see
    `recovery.py` -- so the recovered bounds are checked against the same
    chain as everything else rather than a second implementation of it.

    `set_name` scores with another set than `series.scores` -- one of its
    `also_scores`, or a candidate `eligible_sets` is testing.
    """
    if points is None:
        points = load_points(series)
    set_name = set_name or series.scores
    return [
        dataclasses.replace(r, scored_by=set_name)
        for r in _run_series(series, points, set_name)
    ]


def _run_series(
    series: SeriesMetadata, points: list[Point], set_name: str | None
) -> list[Record]:
    if set_name is None:
        return [
            Record(series, p.x, p.y, None, False, None, "not scored by any set")
            for p in points
        ]
    if series.x_axis not in ("e_plus", "alpha_deg", "f_ratio", "Re_Dh"):
        return [
            Record(series, p.x, p.y, None, False, None, f"x axis is {series.x_axis}")
            for p in points
        ]
    if series.x_axis == "Re_Dh":
        refusal = None
        if series.y_axis not in RE_PATH_Y_AXES:
            refusal = (
                f"the Re_Dh path only predicts {', '.join(RE_PATH_Y_AXES)}, "
                f"not {series.y_axis}"
            )
        elif not series.geometry:
            refusal = "friction depends on e/D and W/H; no rig geometry recorded"
        elif series.geometry.get("turbulated_walls", 2) != 2:
            # Han's decomposition assumes two OPPOSITE ribbed walls; a
            # one-side-turbulated passage is a different configuration.
            refusal = "one-side-turbulated passage; the set's f decomposition is for two ribbed walls"
        if refusal is not None:
            return [Record(series, p.x, p.y, None, False, None, refusal) for p in points]
    if series.x_axis == "f_ratio" and series.y_axis != "Nu_ratio":
        return [
            Record(
                series, p.x, p.y, None, False, None,
                f"the f_ratio path only predicts Nu_ratio, not {series.y_axis}",
            )
            for p in points
        ]
    if series.x_axis == "f_ratio" and not series.geometry:
        # Unlike G on the e+ path, both performance ratios depend on e/D
        # through f for EVERY set, so there is no geometry-free case.
        return [
            Record(
                series, p.x, p.y, None, False, None,
                "performance ratios depend on e/D through f; no rig geometry recorded",
            )
            for p in points
        ]
    if series.x_axis == "alpha_deg" and not series.y_axis.startswith("R"):
        # The alpha-indexed path below only computes R (R, R_normalised,
        # R_normalised_angled -- every R-family name in this dataset starts
        # with "R"), the only quantity independent of e+/Re in every set
        # implemented so far.
        return [
            Record(
                series, p.x, p.y, None, False, None,
                f"alpha-indexed scoring is only implemented for an R "
                f"quantity, not {series.y_axis}",
            )
            for p in points
        ]

    rib_set = SETS[set_name]()

    gbar_block = _gbar_reason(series)
    if gbar_block is not None:
        return [
            Record(series, p.x, p.y, None, False, None, gbar_block)
            for p in points
        ]

    if series.class_confidence == "disputed" and _binds_geometry(rib_set):
        # A disputed label is usable while nothing depends on it. The
        # moment a set binds geometry, the label IS the input, and a
        # disputed input must not be fed in silently.
        return [
            Record(
                series, p.x, p.y, None, False, None,
                f"class label disputed and {rib_set.name} binds geometry",
            )
            for p in points
        ]

    if not series.geometry:
        # Only series with no recorded rig geometry need this. The figure
        # 4.46 legend states one per class, so those are exempt.
        _assert_geometry_free(rib_set)

    unsupported = _binding_reason(rib_set, series)
    if unsupported is not None:
        return [
            Record(series, p.x, p.y, None, False, None, unsupported) for p in points
        ]

    geom = _probe_geometry(rib_set)
    if series.geometry:
        # The figure 4.46 legend states each class's geometry, so use it.
        # It matters for R, whose plotted ordinate is normalised by P/e.
        geom.e_D = float(series.geometry.get("e_D", geom.e_D))
        geom.p_e = float(series.geometry.get("p_e", geom.p_e))
        geom.W_H = float(series.geometry.get("W_H", geom.W_H))
    # Rib angle, which _probe_geometry defaults to 90. Omitting it scored
    # every angled series as though its ribs were transverse: figure 4.51's
    # 45 and 60 deg classes were evaluated at 90 against
    # han_park_1988_angled, the one set whose whole subject is rib angle.
    # `alpha_deg` is a top-level field, not part of the geometry mapping,
    # which is how it was missed when e_D/p_e/W_H were wired through.
    if series.alpha_deg is not None:
        geom.alpha_deg = float(series.alpha_deg)

    if series.x_axis == "Re_Dh" and series.y_axis in NU_RATIO_AXIS_STATUS:
        return _nu_ratio_vs_re(series, rib_set, set_name, geom, points)

    if series.x_axis == "Re_Dh":
        # Measured passage friction against Re directly: no e+ bisection and
        # no R conversion on the measured side. Comparing in f, not R, keeps
        # the law-of-the-wall formalism out of the arbitration -- which
        # matters, since that formalism is part of what is in question.
        # A ratio over a smooth baseline (Naik's f/f0) is first multiplied
        # back by that baseline, through the conversion conventions.py
        # registers -- never a constant written here.
        to_fbar = None
        if series.y_axis == "f_ratio_petukhov":
            from validation.cooling import conventions

            resolution = conventions.resolve(series, set_name)
            if (
                resolution.status != conventions.STATUS_CONVERTED
                or resolution.conversion is None
            ):
                return [
                    Record(
                        series, p.x, p.y, None, False, None,
                        f"f_ratio_petukhov needs a registered conversion to the "
                        f"set's basis; got {resolution.status}: {resolution.detail}",
                    )
                    for p in points
                ]
            to_fbar = resolution.conversion.factor
        records = []
        for p in points:
            res = cb.evaluate_rib(rib_set, geom, p.x)
            fbar, _ = _channel_fbar(res.f, geom.W_H, p.x)
            measured = p.y * to_fbar(p.x) if to_fbar is not None else p.y
            records.append(Record(series, p.x, measured, fbar, res.extrapolated, p.x))
        return records

    if series.x_axis == "f_ratio":
        records = []
        for p in points:
            found = _at_f_ratio(rib_set, geom, p.x)
            if found is None:
                records.append(
                    Record(series, p.x, p.y, None, True, None, "f ratio unreachable")
                )
                continue
            nu_ratio, extrapolated, re = found
            records.append(Record(series, p.x, p.y, nu_ratio, extrapolated, re))
        return records

    if series.x_axis == "alpha_deg":
        # R vs alpha (figure 4.47's own axis): no e+/Re bisection needed --
        # R does not depend on Re for any set in this project. alpha itself
        # is the point's x-value, not a fixed series field, so extrapolation
        # is read per point from evaluate_rib's own flag rather than the
        # series-level _binding_reason check above (which only applies when
        # alpha is fixed for the whole series).
        records = []
        for pt in points:
            geom.alpha_deg = pt.x
            res = cb.evaluate_rib(rib_set, geom, ARBITRARY_RE)
            predicted = res.R
            if series.y_axis in ("R_normalised", "R_normalised_angled"):
                # The figure plots R divided by its own normalisers, not raw
                # R -- dividing them back out here, from the SAME fields
                # evaluate_rib used, rather than duplicating Eq. 4.17's
                # formula. Mirrors what _at_eplus does for han_1988's
                # R_pe division below.
                predicted = _normalised_R(rib_set, geom, res.R)
            records.append(
                Record(series, pt.x, pt.y, predicted, res.extrapolated, ARBITRARY_RE)
            )
        return records

    records: list[Record] = []
    for p in points:
        found = _at_eplus(rib_set, geom, p.x)
        if found is None:
            records.append(
                Record(series, p.x, p.y, None, True, None, "e+ unreachable")
            )
            continue
        g, r_abs, r_norm, extrapolated, re, g_bar = found
        # Explicit per quantity. This was a catch-all `else: predicted = g`,
        # so any y_axis that was not G_bar or R_normalised -- an absolute R
        # among them -- was silently scored against G. Nothing hit it only
        # because every absolute-R series carried `scores: null`.
        if series.y_axis == "G":
            predicted = g
        elif series.y_axis == "G_bar":
            predicted = g_bar
        elif series.y_axis == "R_normalised":
            predicted = r_norm
        elif series.y_axis == "R":
            predicted = r_abs
        else:
            records.append(
                Record(series, p.x, p.y, None, extrapolated, re,
                       f"no prediction for y_axis {series.y_axis} on the e+ path")
            )
            continue
        records.append(Record(series, p.x, p.y, predicted, extrapolated, re))
    return records


def _nu_ratio_vs_re(
    series: SeriesMetadata,
    rib_set: "cb.RibCorrelationSet",
    set_name: str,
    geom: cb.RibGeometry,
    points: list[Point],
) -> list[Record]:
    """Ribbed-wall Nu ratio against Re, on the Dittus-Boelter basis.

    The prediction is the chain's own ribbed-wall Nusselt number, St_r Re Pr
    at the set's Pr, over 0.023 Re^0.8 Pr^0.4 -- what _performance_ratios
    calls Nu_ratio. The MEASURED side is on whatever smooth reference its
    source used (CR-4396: Kays and Perkins), so it goes through the
    conversion conventions.py registers, never a constant written here. The
    axis name fixes which resolution is acceptable (NU_RATIO_AXIS_STATUS): a
    Kays-Perkins axis that resolved direct, or a Dittus-Boelter axis that
    needed a conversion, contradicts its own declaration and is refused, so
    no conversion is ever skipped or applied twice.
    """
    from validation.cooling import conventions

    resolution = conventions.resolve(series, set_name)
    wanted = NU_RATIO_AXIS_STATUS[series.y_axis]
    ok = resolution.status == wanted and (
        wanted == conventions.STATUS_DIRECT or resolution.conversion is not None
    )
    if not ok:
        need = (
            "a registered conversion to the set's basis"
            if wanted == conventions.STATUS_CONVERTED
            else "the set's own Dittus-Boelter basis, unconverted"
        )
        return [
            Record(
                series, p.x, p.y, None, False, None,
                f"{series.y_axis} needs {need}; "
                f"got {resolution.status}: {resolution.detail}",
            )
            for p in points
        ]
    to_db = (
        resolution.conversion.factor()
        if wanted == conventions.STATUS_CONVERTED
        else 1.0
    )
    pr = rib_set.valid_Pr
    records = []
    for p in points:
        res = cb.evaluate_rib(rib_set, geom, p.x)
        nu_db = SMOOTH_NU_COEF * p.x**SMOOTH_NU_RE_EXP * pr**SMOOTH_NU_PR_EXP
        predicted = res.St_r * p.x * pr / nu_db
        records.append(
            Record(series, p.x, p.y * to_db, predicted, res.extrapolated, p.x)
        )
    return records


def sets_for(series: SeriesMetadata) -> tuple[str, ...]:
    """Every set a series is scored by: `scores`, then `also_scores`."""
    if series.scores is None:
        return ()
    return (series.scores, *series.also_scores)


def eligible_sets(series: SeriesMetadata) -> tuple[str, ...]:
    """Every rib set this harness CAN score the series with.

    Derived by rule, never chosen: the set binds the series' rib shape and
    angle, its output MATCHES the series' measured quantity (directly or by
    a registered conversion -- an undeclared convention is scorable for the
    series' own set but qualifies no other), and the runner returns at
    least one prediction. A set validated at one rib angle does not qualify
    for a series that does not say its angle: figure 4.47 pools 30-90 deg,
    and a 90 deg set would score it as though every rib were transverse.

    Every scored rib series must list exactly these, so the choice of which
    sets to show cannot be made after seeing their scores.
    """
    from validation.cooling import conventions

    accepted = (conventions.STATUS_DIRECT, conventions.STATUS_CONVERTED)
    angle_known = series.alpha_deg is not None or series.x_axis == "alpha_deg"
    out = []
    for name in SETS:
        if name != series.scores and conventions.resolve(series, name).status not in accepted:
            continue
        if not conventions.resolve(series, name).scorable:
            continue
        alpha = SETS[name]().valid_alpha
        if not angle_known and alpha.hi > 0.0 and alpha.hi == alpha.lo:
            continue
        try:
            records = run_series(series, set_name=name)
        except ValueError:
            # _assert_geometry_free: the set binds geometry the series
            # does not record -- not scorable, so not eligible.
            continue
        if any(r.predicted is not None for r in records):
            out.append(name)
    return tuple(sorted(out))


def run_all(dataset: list[SeriesMetadata]) -> list[Record]:
    out: list[Record] = []
    for series in dataset:
        for name in sets_for(series) or (None,):
            out.extend(run_series(series, set_name=name))
    return out
