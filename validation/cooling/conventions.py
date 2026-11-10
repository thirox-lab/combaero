"""What each measured quantity is REFERENCED TO, declared rather than implied.

THE FAILURE THIS EXISTS TO CATCH. Two numbers for the same physical thing
under different conventions are different numbers, and every coordinate is
in range, so no bounds check can see it. #389 calls it a category error
rather than a disagreement, and this project has now hit it three times:

  * **Rohde's `Cd` is referenced to duct TOTAL pressure; McGreehan and
    Schotsch's to STATIC.** Related by `sqrt(VHR/(VHR-1))` -- 1.33 at
    VHR 2. `orifice_runner` applies it correctly and always has, but
    nothing DECLARED that a conversion was needed, so nothing would have
    noticed a future series arriving the other way round.
  * **Andrews 86-GT-225's `h_m` is over the HOLE LENGTH; 88-GT-290's `h`
    is on the PLATE AREA.** Same author, same group, same figure number.
    Taking one for the other scores -60% instead of -10%; it was caught
    by reading both nomenclatures by hand (#430).
  * **A gas-side `h` is meaningless without its unblown baseline.** The
    augmentation Andrews' Fig. 10 requires (plate B, with Baldauf's film)
    is 2.7-4.7 against a fully developed Dittus-Boelter and 1.5-2.7 against
    the same duct entry-corrected -- a factor of 1.76 from the baseline
    alone (#431; table in extractions/andrews_effusion_overall_eta.md, #476).

THE SHAPE OF THE FIX. A measured series declares what its y quantity is
referenced to; a correlation set declares what it produces; scoring either
matches them, applies a REGISTERED conversion, or refuses. It never
silently compares two different quantities.

WHY THE DECLARATION IS A TABLE AND NOT 103 YAML EDITS. The convention is a
property of (source, quantity), not of the individual series -- every
Rohde `Cd` is total-referenced because that is how Rohde's rig was
instrumented. One auditable table beats the same string copied 103 times
and drifting. `SERIES_OVERRIDE` exists for the exception, and
`test_every_scored_series_declares_a_convention` makes a new source fail
loudly rather than default.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

# ---------------------------------------------------------------------------
# The conventions themselves. The string is the declaration; the comment is
# what a reader needs to know to tell two of them apart.
# ---------------------------------------------------------------------------

CONVENTION_NOTES: dict[str, str] = {
    # Discharge coefficient
    "Cd_total": "Cd referenced to duct TOTAL pressure (Rohde 1969)",
    "Cd_static": "Cd referenced to duct STATIC pressure (McGreehan-Schotsch)",
    "Cd_coolant_total_exit_static": (
        "Cd of a crossflow-fed hole with ideal flow from the COOLANT "
        "STAGNATION pressure to the MAINSTREAM (exit) STATIC pressure, "
        "compressible up to choking, against PR = pc+/p_inf (Hay, Lampard & "
        "Benmansour 1983, p. 245). Upstream total like Cd_total, but with a "
        "second crossflow on the exit side and no VHR, so Cd_total's "
        "conversion does not apply. Area basis: Gritsch, Schulz & Wittig "
        "(1998, p. 559) print pi D^2/4 on the INLET (metering) diameter for "
        "every hole shape, so their shaped-hole Cd can exceed 1; Hay et al. "
        "do not print theirs (plain holes, so inlet = exit)"
    ),
    # Compressible orifice mass flow (#854)
    "K_reservoir_total": (
        "mass-flow coefficient K = m_dot / (A sqrt(p0 rho0)) on the GEOMETRIC "
        "orifice area and the upstream RESERVOIR (stagnation) state, against "
        "r = p/p0 (Jobson 1955). Not a Cd: it carries the isentropic flux "
        "K_N(r) as well as the contraction"
    ),
    # Rib roughness functions
    "G_ribbed_wall": "heat-transfer roughness function on the RIBBED wall only",
    "G_four_wall": "G_bar, averaged over two ribbed and two smooth walls",
    "R_ribbed_wall": "friction roughness function on the ribbed wall",
    "R_normalised_pe": "R divided by (P/e/10)^0.35, Han's figure 4.46 ordinate",
    "Nu_turbulated_wall": (
        "absolute Nusselt number on the TURBULATED wall, on D_H and the "
        "projected wall area (Taslim & Spring 1987)"
    ),
    "f_channel_two_ribbed": (
        "passage-average FANNING friction of a channel with two opposite "
        "ribbed walls -- Han's measured fbar, NOT the four-sided f_r = "
        "fbar + (H/W)(fbar - f_s) his R is built on"
    ),
    "f_ratio_passage_petukhov": (
        "passage-average FANNING friction over Petukhov's smooth-tube FANNING "
        "f0 = 2 [2.236 ln Re - 4.639]^-2, both on the passage D_h (Naik et al. "
        "2014, Eqs. 4-5). Not Han's 0.046 Re^-0.2 baseline: multiplying by "
        "f0(Re) recovers the passage fbar exactly, so the conversion is "
        "definition"
    ),
    "R_normalised_angled": "R normalised for the angled-rib quadratic form",
    "Nu_ratio_ribbed_DB_vs_fbar": (
        "RIBBED-side Nu over Dittus-Boelter 0.023 Re^0.8 Pr^0.4, plotted "
        "against CHANNEL-AVERAGE f over 0.046 Re^-0.2 (Han, Zhang & Lee 1991, "
        "Eqs. 2/4). Not the four-sided f: 1.8x apart at e/D 0.0625. Nu is "
        "averaged over X/D 0-20 including the entrance region"
    ),
    "Nu_ratio_ribbed_DB": (
        "RIBBED-wall Nu over Dittus-Boelter 0.023 Re^0.8 Pr^0.4 at the same "
        "Re, against Re itself (no friction axis)"
    ),
    "Nu_ratio_ribbed_KP": (
        "RIBBED-wall Nu over Kays and Perkins' smooth-tube constant-wall-"
        "temperature Nu_inf = 0.0176 Re^0.8, a Colburn form with Pr fixed "
        "at 0.72 (NASA CR-4396 Vol. II, p. 18). NOT Dittus-Boelter: 0.873 "
        "of it"
    ),
    # Jet impingement
    "Nu_over_Nu1": "row Nusselt divided by the first-row value, same correlation",
    # Film and effusion
    "eta_adiabatic": (
        "ADIABATIC film effectiveness, (Tg-Taw)/(Tg-Tc). Measured on an "
        "insulated wall, so it carries no heat transfer coefficient"
    ),
    "eta_overall": (
        "OVERALL cooling effectiveness, (Tg-Tw)/(Tg-Tc). Contains the "
        "internal convection and the gas-side coefficient as well as the "
        "film -- NOT the same quantity as eta_adiabatic and not convertible "
        "to it without the whole resistance network"
    ),
    "h_hole_internal_area": "h referenced to the hole internal surface, pi D L",
    "h_plate_approach_area": "h referenced to the plate area, X^2 - pi D^2/4",
    "h_target_plate_area": "h referenced to the impingement target plate's plan area",
}

# ---------------------------------------------------------------------------
# What each SOURCE publishes, per quantity. Keyed (source, y_axis).
# ---------------------------------------------------------------------------

SOURCE_PUBLISHES: dict[tuple[str, str], str] = {
    ("rohde1969", "Cd"): "Cd_total",
    # Kaddah's C = Lj/Lth with Vj from P1 - pj (approach TOTAL to jet static),
    # against K = (P1 - pj)/(P1 - p1), Rohde's VHR. Read by kaddah_runner.
    ("kaddah1964", "Cd"): "Cd_total",
    ("mcgreehan_schotsch1988", "Cd"): "Cd_static",
    # Unscored today (no angle or two-sided crossflow term in combaero);
    # declared so a future set must meet the basis explicitly.
    ("hay1983", "Cd"): "Cd_coolant_total_exit_static",
    # Same pressure pair, printed isentropic ideal flow, inlet-area basis.
    ("gritsch1998", "Cd"): "Cd_coolant_total_exit_static",
    ("jobson1955", "K_mass_flow"): "K_reservoir_total",
    ("han2012", "G"): "G_ribbed_wall",
    ("han_park_lei1984", "G"): "G_ribbed_wall",
    ("lau1990", "G"): "G_ribbed_wall",
    ("han2012", "G_bar"): "G_four_wall",
    ("lau1990", "G_bar"): "G_four_wall",
    ("han_park_lei1984", "R"): "R_ribbed_wall",
    # Taslim & Spring (1987) Fig. 11, f = dP gc / (2 (L/D_H) rho V^2) over
    # the passage -- the same definition as Han's Eq. (1) (#403 item 2).
    ("taslim_spring1987", "f_fanning_passage"): "f_channel_two_ribbed",
    ("taslim_spring1987", "Nu_turbulated"): "Nu_turbulated_wall",
    ("lau1990", "R"): "R_ribbed_wall",
    # Figure 4.48's R panel (#402) plots raw R, unlike figures 4.46/4.47's
    # normalised ordinate -- same quantity as han_park_lei1984/lau1990's "R".
    ("han2012", "R"): "R_ribbed_wall",
    ("han2012", "R_normalised"): "R_normalised_pe",
    # Figure 4.53, Han and Zhang (1992). Their own normalisation is not on
    # disk; the 1991 paper it reprints from prints it (#435).
    ("han2012", "Nu_ratio"): "Nu_ratio_ribbed_DB_vs_fbar",
    ("han2012", "R_normalised_angled"): "R_normalised_angled",
    ("florschuetz1981", "Nu_over_Nu1"): "Nu_over_Nu1",
    ("andrei2014", "eta_adiabatic"): "eta_adiabatic",
    ("murray2018", "eta_adiabatic"): "eta_adiabatic",
    ("colban2011", "eta_adiabatic"): "eta_adiabatic",
    ("andrews1988", "eta_overall"): "eta_overall",
    # Reacting 5-10 bar rig (GT2008-51191). T_g is the NEAR-WALL gas at
    # x/d -10, y/d 4, not a bulk temperature; effusion_network_runner feeds
    # it as the chamber's main-inlet stream, which the plate's T_gas is.
    ("behrendt2008", "eta_overall"): "eta_overall",
    # The pair that cost -60% when confused. Same lab, same figure number.
    ("andrews1986", "h_hole_length"): "h_hole_internal_area",
    ("andrews1988", "h_internal"): "h_plate_approach_area",
    # Same Leeds rig and area as andrews1988, but the WHOLE impingement/
    # effusion wall: its own y_axis so the effusion-internal runner, which
    # claims h_internal, cannot score it as an effusion plate alone.
    ("aldabagh1990", "h_imp_eff_wall"): "h_plate_approach_area",
    ("aldabagh1990", "h_impingement_target"): "h_target_plate_area",
    # NASA CR-4396 Vol. II, Figs. 6.1/6.2 (#452): Kays and Perkins, not DB.
    ("johnson1993", "Nu_ratio_KP"): "Nu_ratio_ribbed_KP",
    # Naik et al. 2014 Fig. 16 (Alstom wedge TE passage): Nu/Nu0 with Nu0 the
    # printed Dittus-Boelter (Eq. 3), averaged between the ribs (rib surfaces
    # excluded). Declared ribbed-wall although the paper does not say which
    # walls Fig. 16 averages: the two ribbed walls are ~80% of the perimeter
    # and a by-eye read of Fig. 11c puts a whole-perimeter mean ~2-5% below
    # theirs -- inside the 9.5% band. See extractions/naik_2014_te_wedge_ribs.md.
    ("naik2014", "Nu_ratio_DB"): "Nu_ratio_ribbed_DB",
    # Fig. 17: Fanning passage f (Eq. 4) over Petukhov's Fanning f0 (Eq. 5).
    ("naik2014", "f_ratio_petukhov"): "f_ratio_passage_petukhov",
}

# A series whose convention differs from the rest of its source's. None yet;
# the hook exists so the exception does not force the table to be abandoned.
SERIES_OVERRIDE: dict[str, str] = {}

# ---------------------------------------------------------------------------
# What each CORRELATION SET produces.
# ---------------------------------------------------------------------------

SET_PRODUCES: dict[str, str] = {
    "han_1988_orthogonal": "G_ribbed_wall",
    "han_park_1988_angled": "G_ribbed_wall",
    "han_1989_narrow_channel": "G_ribbed_wall",
    "rallabandi_2009_high_re": "G_ribbed_wall",
    "florschuetz_1981_inline": "Nu_over_Nu1",
    "mcgreehan_schotsch_1988_cd": "Cd_static",
    "mcgreehan_schotsch_1988_crossflow_cd": "Cd_static",
    "jobson_1955_mass_flow": "K_reservoir_total",
    "baldauf_2002_sellers": "eta_adiabatic",
    "colban_2011_shaped": "eta_adiabatic",
    # One set, two quantities: it produces an internal h, and the overall
    # eta runner builds a network result from it. Declared on the quantity
    # the runner actually compares, which `resolve` is told per series.
    "andrews_1986_effusion_internal": "h_hole_internal_area",
}

# A set that legitimately produces more than one quantity declares the
# others here, so `resolve` can accept any of them rather than needing a
# conversion that is really a different output.
SET_ALSO_PRODUCES: dict[str, tuple[str, ...]] = {
    # The rib sets produce R as well as G, in several normalisations the
    # runner selects between; and G_bar through Han's published 1.2.
    # Nu_ratio_ribbed_DB is St_r Re Pr over Dittus-Boelter at the set's own
    # Pr -- the ribbed-wall Nu the chain already computes.
    "han_1988_orthogonal": (
        "R_ribbed_wall", "R_normalised_pe", "G_four_wall", "f_channel_two_ribbed",
        "Nu_ratio_ribbed_DB_vs_fbar", "Nu_ratio_ribbed_DB",
    ),
    # f_channel_two_ribbed: Han and Park (1988) MEASURE the passage-average
    # fbar of a channel with two opposite ribbed walls (p. 186), build R on
    # the four-sided f = fbar + (H/W)(fbar - f(FD)) (Eq. 6, f(FD) = 0.046
    # Re^-0.2, Eq. 8), and state that fbar "can be predicted by combining
    # equations (7)-(9)" (p. 193). runner._channel_fbar is that inverse.
    "han_park_1988_angled": (
        "R_ribbed_wall", "R_normalised_pe", "R_normalised_angled",
        "G_four_wall", "f_channel_two_ribbed", "Nu_ratio_ribbed_DB",
    ),
    "rallabandi_2009_high_re": (
        "R_ribbed_wall", "R_normalised_pe", "Nu_ratio_ribbed_DB",
    ),
    # Fig. 4.48's R panel plots raw R (no (p/e/10)^0.35 normalisation unlike
    # figures 4.46/4.47), so this set needs no R_normalised_* entry.
    "han_1989_narrow_channel": ("R_ribbed_wall", "Nu_ratio_ribbed_DB"),
    # The effusion set's overall-eta closure is a network result built from
    # the same internal h -- an OUTPUT of the set, not a conversion of it.
    "andrews_1986_effusion_internal": (
        "h_plate_approach_area", "eta_overall",
    ),
}

# Taslim & Spring (1987) ratio sets (#444): turbulated-wall Nu and the
# passage-average Fanning f, both the paper's own quantities.
for _ar, _ed in ((0.5, 0.125), (0.5, 0.250), (1.0, 0.083), (1.0, 0.167),
                 (3.5, 0.053), (3.5, 0.107), (3.5, 0.161)):
    SET_PRODUCES[f"taslim_spring_1987_ar{_ar:.1f}_eD{_ed:.3f}"] = "Nu_turbulated_wall"
    SET_ALSO_PRODUCES[f"taslim_spring_1987_ar{_ar:.1f}_eD{_ed:.3f}"] = ("f_channel_two_ribbed",)

# Han, Zhang and Lee (1991) Table 2: R and G on Han's basis (four-sided f_r,
# ribbed-side St), plus the source's OWN printed G_bar; the 60 deg V set also
# answers figure 4.53's performance curve through the f_ratio path.
for _tag in ("90", "60par", "60crs", "60vee", "60lam", "45par", "45crs", "45vee", "45lam"):
    SET_PRODUCES[f"han_zhang_lee_1991_{_tag}"] = "G_ribbed_wall"
    SET_ALSO_PRODUCES[f"han_zhang_lee_1991_{_tag}"] = (
        "R_ribbed_wall", "G_four_wall", "Nu_ratio_ribbed_DB_vs_fbar",
        "Nu_ratio_ribbed_DB",
    )


# Kays and Perkins (1973) smooth tube at constant wall temperature, as
# CR-4396 Vol. II prints it: Nu_inf = 0.0176 Re^0.8, "with a Prandtl number
# equal to 0.72" (p. 18).
KAYS_PERKINS_COEF = 0.0176
KAYS_PERKINS_PR = 0.72

# Petukhov's smooth-tube friction as Naik et al. (2014) print it, Eq. 5, in
# FANNING form: f0 = 2 [2.236 ln Re - 4.639]^-2. Petukhov's Darcy form
# (0.790 ln Re - 1.64)^-2 over 4 is the same curve to 0.2% (2.236/sqrt(2) =
# 1.5811, 4.639/sqrt(2) = 3.2803); the paper's own constants are kept because
# they define its ordinate.
NAIK_PETUKHOV_LN_COEF = 2.236
NAIK_PETUKHOV_OFFSET = 4.639


def naik_petukhov_fanning(re: float) -> float:
    """Naik et al. (2014) Eq. 5: Petukhov's smooth-tube FANNING friction."""
    return 2.0 / (NAIK_PETUKHOV_LN_COEF * math.log(re) - NAIK_PETUKHOV_OFFSET) ** 2


@dataclass(frozen=True)
class Conversion:
    """A registered way to turn one convention into another.

    `factor` takes whatever the runner can supply as context and returns
    the multiplier. Registered conversions are PURE DEFINITION -- geometry
    or an algebraic identity -- never a fitted correction. A conversion
    that needed fitting would be a model, and would belong in the library
    with its own evidence.
    """

    derivation: str
    applied_by: str
    factor: Callable[..., float] | None = None


CONVERSIONS: dict[tuple[str, str], Conversion] = {
    # Rohde's Cd on total -> the static basis McGreehan-Schotsch produces.
    # The factor diverges as VHR -> 1, which is why those scores are
    # reported per VHR band and never pooled.
    ("Cd_total", "Cd_static"): Conversion(
        derivation="Cd_static = Cd_total * sqrt(VHR/(VHR-1)); 1.33 at VHR 2, "
        "1.01 at VHR 50. See extractions/rohde_1969_orifice.md.",
        applied_by="validation.cooling.orifice_runner.vhr_to_static_cd",
        factor=lambda vhr: math.sqrt(vhr / (vhr - 1.0)),
    ),
    # Andrews' two h bases. Pure geometry: A_h/A with A_h = pi D L and
    # A = X^2 - pi D^2/4. Between 1.4 and 9.8 over his plates.
    ("h_hole_internal_area", "h_plate_approach_area"): Conversion(
        derivation="h_plate = h_hole * (pi D L)/(X^2 - pi D^2/4). Pure "
        "geometry, no correlation. See extractions/"
        "andrews_effusion_internal_h.md.",
        applied_by="validation.cooling.effusion_internal_runner.predict",
        factor=None,  # needs the series geometry; the runner owns it
    ),
    # CR-4396's smooth-tube reference -> Dittus-Boelter. Both are Re^0.8, so
    # the ratio is a constant; the Pr is the source's own declared 0.72,
    # because 0.0176 IS Kays and Perkins evaluated there. Multiplying a
    # measured Nu/Nu_KP by it gives Nu/Nu_DB. It carries the assumption every
    # Nu-ratio comparison does: the ribbed and smooth Nu share Pr^0.4, so a
    # ratio at Pr 0.72 meets the set's at its own 0.70 (1.1% if it did not).
    ("Nu_ratio_ribbed_KP", "Nu_ratio_ribbed_DB"): Conversion(
        derivation="Nu/Nu_DB = (Nu/Nu_KP) * 0.0176 / (0.023 * 0.72^0.4) = "
        "0.873 (Nu/Nu_KP). Both references are Re^0.8; Pr 0.72 is CR-4396's "
        "own (Vol. II p. 18). Definition, not a fit. See "
        "extractions/johnson_1993_cr4396.md.",
        applied_by="validation.cooling.runner.run_series",
        factor=lambda: KAYS_PERKINS_COEF / (0.023 * KAYS_PERKINS_PR**0.4),
    ),
    # Naik's f/f0 -> the passage-average Fanning fbar Han's sets are compared
    # in (_channel_fbar). Both Fanning, both on the passage D_h and Re: f0(Re)
    # is the paper's own printed denominator, so multiplying it back is
    # definition. Which SETS may use it is a separate question --
    # SET_ALSO_PRODUCES declares f_channel_two_ribbed for han_1988_orthogonal
    # and han_park_1988_angled only; any other set refuses it (unregistered).
    ("f_ratio_passage_petukhov", "f_channel_two_ribbed"): Conversion(
        derivation="fbar = (f/f0) * f0(Re), f0 = 2 [2.236 ln Re - 4.639]^-2 "
        "(Naik et al. 2014 Eq. 5, Fanning). Definition, not a fit. See "
        "extractions/naik_2014_te_wedge_ribs.md.",
        applied_by="validation.cooling.runner.run_series",
        factor=naik_petukhov_fanning,
    ),
}

# Pairs that must NEVER be converted, with why. Listing them is the point:
# an unregistered pair is refused by default, but these are the ones a
# reader might assume are convertible.
INCOMPATIBLE: dict[tuple[str, str], str] = {
    ("eta_adiabatic", "eta_overall"): (
        "An adiabatic effectiveness is measured on a wall passing no heat, "
        "so it contains no coefficient. An overall effectiveness contains "
        "the internal convection AND the gas-side coefficient. Going from "
        "one to the other needs the whole resistance network, which is a "
        "model and not a conversion -- and an overall-eta correlation used "
        "as a closure would double-count the convection the network "
        "computes (#387)."
    ),
    ("G_ribbed_wall", "G_four_wall"): (
        "Han's G_bar = 1.2 G is a published CONSTANT fitted at 90 degrees, "
        "not an identity: the measured ratio spans 1.096 to 1.413 across "
        "16 configurations, with 69% of the variance between rigs (#401). "
        "It is applied as the source's own correlation, with its cost "
        "measured -- not registered here as though it were definitional."
    ),
}

STATUS_DIRECT = "direct"
STATUS_CONVERTED = "converted"
STATUS_UNDECLARED = "undeclared"
STATUS_INCOMPATIBLE = "incompatible"
STATUS_UNREGISTERED = "unregistered"
# A series no set scores. Distinct from `undeclared`, which means a
# scored series whose convention the table does not yet carry -- the
# count that must not grow.
STATUS_NOT_SCORED = "not-scored"


@dataclass(frozen=True)
class Resolution:
    status: str
    detail: str
    conversion: Conversion | None = None

    @property
    def scorable(self) -> bool:
        """Whether a comparison may proceed.

        `undeclared` scores: the table is incomplete by construction while
        sources are being added, and refusing would hide working results
        behind bookkeeping. It is COUNTED instead, and
        `test_the_undeclared_count_does_not_grow` is what keeps that from
        becoming permanent.
        """
        return self.status in (
            STATUS_DIRECT,
            STATUS_CONVERTED,
            STATUS_UNDECLARED,
            STATUS_NOT_SCORED,
        )


def series_convention(series) -> str | None:
    """What this series' y quantity is referenced to, or None if undeclared."""
    override = SERIES_OVERRIDE.get(series.label)
    if override:
        return override
    return SOURCE_PUBLISHES.get((series.source.name, series.y_axis))


def resolve(series, set_name: str | None) -> Resolution:
    """Whether `series` and `set_name` measure the same quantity.

    Returns a Resolution rather than raising, so the scorecard can REPORT
    a refusal beside the rows it refused -- a silent drop is the failure
    mode this is replacing.
    """
    if not set_name:
        return Resolution(STATUS_NOT_SCORED, "not scored by any set")

    measured = series_convention(series)
    produced = SET_PRODUCES.get(set_name)
    if measured is None or produced is None:
        missing = []
        if measured is None:
            missing.append(f"series ({series.source.name}, {series.y_axis})")
        if produced is None:
            missing.append(f"set ({set_name})")
        return Resolution(
            STATUS_UNDECLARED, "no convention declared for " + " and ".join(missing)
        )

    if measured == produced or measured in SET_ALSO_PRODUCES.get(set_name or "", ()):
        return Resolution(STATUS_DIRECT, measured)

    why = INCOMPATIBLE.get((produced, measured)) or INCOMPATIBLE.get(
        (measured, produced)
    )
    if why:
        return Resolution(
            STATUS_INCOMPATIBLE, f"{produced} vs {measured}: {why}"
        )

    # To the set's primary quantity, or to any other it also produces: a rib
    # set's Nu ratio is not its primary G, but it is still its output.
    for target in (produced, *SET_ALSO_PRODUCES.get(set_name, ())):
        conversion = CONVERSIONS.get((measured, target))
        if conversion is not None:
            return Resolution(
                STATUS_CONVERTED, f"{measured} -> {target}", conversion
            )

    return Resolution(
        STATUS_UNREGISTERED,
        f"{produced} vs {measured}: no registered conversion. Scoring two "
        "different quantities is a category error, not a disagreement -- "
        "register the conversion with its derivation, or declare the pair "
        "INCOMPATIBLE.",
    )
