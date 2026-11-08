# Extraction: Hay, Lampard and Benmansour (1983), hole Cd with crossflow on both sides

**Status: DATA ONLY -- digitised 2026-11-08, user-confirmed on scan overlays
("all figs okay").** No correlation is extracted and none is implemented: the
paper prints none. The data live in
`validation/cooling/data/hay1983/` (25 series, 174 marks, all `scores: null`).

## For the reviewer

- **The two suspected carpet misprints** (below) are transcribed as printed and
  not used anywhere. Settling them needs the Benmansour thesis.
- **One Fig. 6 mark has no series** (`fig6_unassigned_open_square`, filed
  `class_confidence: disputed`); one Fig. 12 filled square is doubtful
  (`fig12_filled_squares`, `provisional`).
- **Whether `L` in L/D is the inclined hole length** is not stated. It matters
  only if a model built on this data maps a wall thickness onto L/D.

### What I still need sent

- Benmansour, S. (1983), M.Phil. thesis, Nottingham: the tabulated runs behind
  every figure. Not held (dropped from wanted, 2026-10-28); Hay & Lampard 1996
  Fig. 4 reprints the inlet-crossflow points.

## Source, pinned

> Hay, N., Lampard, D. and Benmansour, S. (1983). Effect of Crossflows on the
> Discharge Coefficient of Film Cooling Holes. *ASME J. Engineering for Power*
> **105**(2), 243-248. Paper No. 82-GT-147.
> `sources/orifices/hay_lampard_benmansour_1983_JEP_105_243_crossflow_cd.pdf`
> (gitignored, copyrighted). Register entry in `validation/sources.yaml`.

A scan, not vector: each page is one 1-bit CCITT image at 300 dpi with an OCR
layer on top. Worked on the native bitmaps (`pdfimages -png`), losslessly.

## Extracted items

| # | item | value | state |
|---|---|---|---|
| 1 | Hole | D = 3.5 mm, four holes in one row, S/D = 2.5; flat plate between two channels (Fig. 3), cross-shaped plate for perpendicular flows (Fig. 4). Number of holes shown negligible | confirmed |
| 2 | L/D | 2 and 6 tested in the parallel rig; **every published figure is L/D 6**; "for the range of L/D of interest (2-10) it had hardly any effect on Cd" (p. 244). Whether L is the inclined length: not stated | confirmed / missing |
| 3 | Angle | alpha = hole axis to the mainstream: 30, 60, 90 deg; 60 deg fell between and was dropped. Parallel rig: the 30 deg hole leans downstream with both streams. Perpendicular rig: coolant at 90 deg to the hole's lean plane | confirmed |
| 4 | Inlet edge | not stated for the test plate; production holes are "round holes with sharp edges" (p. 244) -- read as sharp | suspect |
| 5 | Range | room-temperature air; hole Re 1e4-1e5 (not printed per point); Mc 0-0.4, Minf 0-0.5; PR 1-2, "just beyond the sonic condition within the hole" | confirmed |
| 6 | **Cd definition** | actual / ideal flow, the ideal flow "based on the coolant stagnation pressure and mainstream static pressure" (p. 245); "this takes into account the dynamic head of the flow approaching the hole". Abscissa PR = pc+/p_inf on the same pair. Ideal-flow formula (compressible isentropic implied) and area basis not printed | confirmed / missing |
| 7 | **Uncertainty (Table 3)** | a measurement band. Minf 0.3: 8-3.5 % at 1.01 < PR < 1.2, 3.5-1.5 % at 1.2-1.5, 1.5 % above 1.5; one-sided crossflow ~6 / 2.5 / 1.5 % | confirmed |
| 8 | Fig. 11 caption | "at Minf = .3" is a misprint: its legend holds Mc ~0.29 and varies Minf (the perpendicular twin of Fig. 5) | confirmed |

The convention (item 6) is registered in `validation/cooling/conventions.py` as
`Cd_coolant_total_exit_static`, distinct from Rohde's `Cd_total` (no exit
crossflow, VHR abscissa) and McGreehan-Schotsch's `Cd_static`; no conversion is
registered, so a future set producing either must refuse rather than score.

Item 7 is PR-dependent, so the schema's scalar `uncertainty` is null on every
series; the band is recorded here and in the metadata header. It is never to be
replaced by a number derived from a model's error on these series.

## Figures digitised

| Fig. | flows | alpha | held | varied | series | marks | flagged |
|---|---|---|---|---|---|---|---|
| 5 | parallel | 90 | Mc 0.29 | Minf 0, .12, .31, .52 | 4 | 26 | 3 |
| 6 | parallel | 90 | Minf 0.31 | Mc 0, .09, .19, .29, .40 | 5 + 1 unassigned | 36 | 10 |
| 8 | parallel | 30 | Minf 0.30 | Mc 0, .10, .19, .29 | 4 | 32 | 4 |
| 11 | perpendicular | 90 | Mc ~0.29 | Minf 0, .11, .33, .52 | 4 | 28 | 5 |
| 12 | perpendicular | 30 | Minf 0.31 | Mc 0, .10, .20, .30, .38 | 5 | 33 | 10 |
| 15 | parallel | 30 | Mc 0 | Minf .30, .52 (x = u_j/u_inf) | 2 | 19 | 1 |

Legend values are recorded per series (e.g. Fig. 6's pointers are Minf .30,
the others .31). Flagged marks -- merged with another symbol or a curve, or in
the overplotted PR 1.9-2.1 clusters of Figs. 6 and 12 -- are kept and listed in
each series' `cross_check`.

**Not digitised.** Fig. 7 (90 deg, Minf 0, against Rohde, log abscissa): its
axes calibrate, but the present-work symbols are printed over Rohde's hatched
markers and could not be separated; Rohde is held and digitised
(`data/rohde1969/`). Fig. 15's "Prediction" curve is the authors' model, not
data. The grey halftone curve bands of Figs. 5, 6, 11, 12 are curves.

## Digitisation method and calibration evidence

Automated candidate detection on the 300 dpi bitmap, every candidate reviewed
and assigned by hand, missed marks picked manually; then a scan overlay
(scan + calibrated frame + points + numbered suspects) for the user, confirmed
2026-11-08. Scan and overlays stay out of git.

Axes were fitted as long dark lines (robust fit, rms 0.4-0.6 px) and ticks as
short perpendicular runs, one linear map per axis; every tick re-predicted:

| Fig. | x max abs residual | y max abs residual | ticks not found (x / y) | y-axis tilt | x-axis tilt |
|---|---|---|---|---|---|
| 5 | 0.0018 | 0.0012 | 0 / 0 | -0.35 deg | -0.05 deg |
| 6 | 0.0016 | 0.0025 | 0 / 0 | -0.85 | -0.20 |
| 8 | 0.0007 | 0.0016 | 0 / 0 | -0.95 | +0.38 |
| 11 | 0.0018 | 0.0015 | 2 / 2 | -0.70 | +0.11 |
| 12 | 0.0010 | 0.0024 | 7 / 2 | -0.73 | +0.27 |
| 15 | 0.0022 | 0.0049 | 0 / 0 | -0.62 | +0.12 |

Both axes are **linear** on every figure (the log/linear trap was checked; only
Fig. 7's abscissa is log). The tick files are committed beside the data
(`fig<N>_xticks.csv`, `_yticks.csv`, `_corners.csv`), not listed under `files:`.

**Drafting slant, not removed.** The y axis leans -0.35 to -0.95 deg in every
panel while the x axis tilts between -0.20 and +0.38 deg; a page rotation would
tilt both equally, so this is a drafting shear. The map is oblique (constant-PR
lines parallel to the drawn y axis). The committed `x` is that reading; `x_perp`
(constant PR perpendicular to the x axis) is kept as the sensitivity, +0.002 to
+0.008 in PR, largest at the top of Figs. 6/8. Cd is unaffected. No correction
is applied to the files.

**Third channel.** Figs. 5 and 6 both plot the same run (crosses, Minf .31,
Mc .29). Digitised independently in the two panels, all 7 marks agree to
|dPR| <= 0.010 and |dCd| <= 0.009 (mean ~0.005).

`uv run python -m validation.cooling.verify` checks every series' span and
count against its figure card.

## Printed carpet values (PR >= 1.6 plateau), Figs. 9, 10, 13, 14

Printed NUMBERS, transcribed, not digitised points -- not filed as series
(the carpets' axes were not calibrated, so no figure card could check them).
Cd as printed; (Mc, Minf) of each label read from its dot by eye, +-0.02.

| Fig. | Minf (read) | Mc (read) : Cd (printed) |
|---|---|---|
| 9 parallel 90 | 0.50-0.52 | 0.00 : 0.82; 0.09 : 0.78; 0.29 : 0.76; 0.40 : 0.74 |
| 9 | 0.29-0.30 | 0.00 : 0.80; 0.09 : 0.79; 0.18 : 0.78; 0.29 : 0.77; 0.39 : 0.76 |
| 9 | 0.12 | 0.00 : 0.79; 0.11 : 0.78; 0.19 : 0.78; 0.29 : **0.83** |
| 9 | 0.01 | 0.00 : 0.81; 0.29 : **0.83** |
| 10 parallel 30 | 0.52 | 0.00 : 0.70; 0.12 : 0.76; 0.20 : 0.83; 0.29 : 0.93 |
| 10 | 0.30 | 0.00 : 0.72; 0.10 : 0.74; 0.20 : 0.84; 0.30 : 0.95 |
| 10 | 0.12 | 0.00 : 0.73; 0.10 : 0.82; 0.20 : 0.85; 0.29 : 0.94; 0.39 : 0.97 |
| 10 | 0.02 | 0.00 : 0.73; 0.30 : 0.95 |
| 13 perpendicular 90 | 0.52 | 0.00 : 0.78; 0.11 : 0.77; 0.21 : 0.75; 0.29 : 0.74 |
| 13 | 0.33 | 0.00 : 0.80; 0.10 : 0.82; 0.20 : 0.80; 0.30 : 0.76; 0.41 : 0.75 |
| 13 | 0.12 | 0.00 : 0.79; 0.10 : 0.82; 0.22 : 0.78; 0.30 : 0.76 |
| 13 | 0.01 | 0.00 : 0.78; 0.10 : 0.81; 0.23 : 0.75; 0.30 : 0.74 |
| 14 perpendicular 30 | 0.52-0.54 | 0.00 : 0.70; 0.12 : 0.69; 0.18 : 0.68; 0.27 : 0.67 |
| 14 | 0.32 | 0.00 : 0.71; 0.09 : 0.72; 0.18 : 0.70; 0.27 : 0.70; 0.36 : 0.68 |
| 14 | 0.13 | 0.00 : 0.72; 0.09 : 0.73; 0.18 : 0.69; 0.27 : 0.68 |
| 14 | 0.01 | 0.00 : 0.71; 0.09 : 0.74; 0.18 : 0.72; 0.27 : **0.79** |

**Suspected misprints (bold), left as printed.** Fig. 9 prints 0.83 at Mc ~0.3
for Minf ~0 and ~0.12, although its own contour lines (and Fig. 6) fall with
Mc. Fig. 14 prints 0.79 at Mc ~0.3, Minf ~0 against 0.68-0.72 everywhere
around it. Check against the thesis before use. The two Mc = 0, Minf ~0
corners (parallel vs perpendicular rigs, nominally the same flow) differ by
0.03 (90 deg) and 0.02 (30 deg) -- a rough rig-to-rig scatter.

## What the data can test

- **A hole-angle term in Cd** (30 vs 90 deg at equal L/D) and **the direction
  of the coolant crossflow relative to the hole's lean** (parallel vs
  perpendicular rig). No combaero closure carries either, so nothing is scored.
- Against any combaero closure this is **cross-source ACCURACY** and must be
  labelled so in the report (docs/VALIDATION_POLICY.md). **Fidelity** applies
  only to a model built from this paper, and none exists.
- Table 3 is the band such a score is judged against, PR by PR. Matching this
  rig is a user's `Cd`/multiplier choice, never the harness's.

## First look (a look, not a score)

Linear interpolation between neighbouring digitised marks, no extrapolation;
points are scattered, so read +-0.03-0.05 at low PR. Nothing below is pinned
by a test or used to tune anything.

- **30 / 90 deg, parallel, no coolant crossflow, Minf 0.3** (Fig. 8 / Fig. 6):
  0.85-0.86 at PR 1.1-1.2, 0.90 at 1.3-1.5, 0.93 at 1.8; the plateau carpets
  give 0.90 with no crossflow at all. The 30 deg hole discharges ~10 % less.
- **Coolant crossflow co-aligned with the 30 deg lean** (parallel, Minf .3):
  the ratio flips -- Mc .29: 1.25-1.31 at PR 1.3-1.5; carpets 1.14 (Minf 0) to
  1.23 (Minf .3). Cd(Mc .29)/Cd(Mc 0) at 30 deg is +25-33 % (Table 4: "about
  30 % in all ranges").
- **90 deg, parallel** (Fig. 6): coolant crossflow lowers Cd at low PR (Mc .19:
  0.66 / 0.76 of the Mc 0 value at PR 1.1 / 1.2), ~1.0 by PR 1.8 -- Table 4's
  -15 to -20 % at low PR and ~-6 % above 1.5, in sign and size.
- **Perpendicular, 30 deg** (Fig. 12 / Fig. 11, Mc .30, Minf .31): 30/90 =
  0.88-0.94 -- crossflow across the lean plane does not help the angled hole.
- **Mainstream crossflow** (Mc .29 held): parallel 90 deg Minf .52 vs 0,
  0.85-0.86 at PR 1.2-1.3; perpendicular 90 deg shows no reduction within
  scatter.
- **Cross-source, against combaero's plain-hole closures** (sharp inlet, L/d
  6, no angle term; Idelchik1966Thick / McGreehanSchotsch1988 / Lichtarowicz1965
  give 0.77-0.81 at Re 1e4-1e5): Hay's no-crossflow plateau is 0.78-0.81 at
  90 deg (within -5 to +4 %) and 0.71-0.73 at 30 deg (combaero 5-14 % high).
  Caveats: Hay's Cd is coolant-total to mainstream-static and compressible on
  the plateau, combaero's a plenum-to-plenum closure; they coincide only when
  both crossflows are ~0; hole Re per point is not printed; the inlet edge is
  inferred.

## What it does NOT cover

- L/D other than 6 (L/D 2 was measured but not published as data); Fig. 7.
- Hole Reynolds number per point, inlet-edge geometry, the ideal-flow formula.
- Hot or density-ratio-varied flow (all isothermal room-temperature air).
- Thin plates (impingement, L/D ~1) or the Gc/Gj ~125 stagnant-row regime.
- Shaped holes, compound angles, or rows other than four holes at S/D 2.5.

## Review log

| date | who | entry |
|---|---|---|
| 2026-11-08 | Claude | Digitised Figs. 5, 6, 8, 11, 12, 15 (174 marks, 33 flagged) and transcribed the carpets of Figs. 9, 10, 13, 14; scan overlays and per-figure review pages sent to the user. |
| 2026-11-08 | user | **CONFIRMED** -- "all figs okay". Promoted to `validation/cooling/data/hay1983/`, unscored. |
