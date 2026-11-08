# Extraction: Gritsch, Schulz and Wittig (1998), Cd of cylindrical and shaped 30 deg holes

**Status: DATA ONLY -- digitised 2026-11-08, user-confirmed over three review
rounds.** The paper prints no correlation, and none is extracted or
implemented. The data are in `validation/cooling/data/gritsch1998/`: 63 series
and 500 marks (465 confirmed, 35 provisional), all `scores: null`.

## For the reviewer

- **35 provisional marks** are filed apart in `_provisional` series. Most sit
  where several curves collapse onto one another (Figs. 7, 9, 12). They are
  kept because the collapse is data (user, 2026-11-08). A coincident point
  appears with the same coordinates in every series that passes through it.
- **The shaped-hole gap below** comes from a first look, not a score. Closing
  it needs a shaped-hole Cd closure. Neither this paper nor anything held
  supplies one, and the authors' 1997 ISABE correlation is not held.

### What I still need sent

- Gritsch, Schulz and Wittig (1997), ISABE: the correlation fitted to these
  data. Wanted only if a shaped-hole closure is to be built from it.

## Source, pinned

> Gritsch, M., Schulz, A. and Wittig, S. (1998). Discharge Coefficient
> Measurements of Film-Cooling Holes With Expanded Exits. *ASME J.
> Turbomachinery* **120**(3), 557-563. Paper No. 97-GT-165.
> `sources/orifices/gritsch_schulz_wittig_1998_JT_120_557_expanded_exit_hole_cd.pdf`
> (gitignored, copyrighted). The register entry is in `validation/sources.yaml`.

The source is a 1-bit CCITT raster at 300 dpi. It has no vector data, and the
figures were read from the bitmaps.

## Extracted items

| # | item | value | state |
|---|---|---|---|
| 1 | Holes | Single scaled-up holes. alpha 30 deg, sharp-edged, smooth, D = 10 mm, L/D = 6 for all three shapes | confirmed |
| 2 | Fan-shaped | Cylindrical inlet 2D long, then 14 deg lateral expansion to a 30 mm exit width. A_exit/A_inlet = 3.0 (areas normal to the hole axis) | confirmed |
| 3 | Laidback fan-shaped | As item 2, plus a 15 deg forward expansion (layback) to a 40 mm exit length. Area ratio 3.1 | confirmed |
| 4 | Rig | Coolant channel 20 x 60 mm, COFLOWING with the main flow. Main channel 90 x 41 mm (32 mm supersonic). Tt_c ~290 K, Tt_c/Tt_m = 1. Re_D up to 2.5e5 (internal) and 2.1e5 (main). Tu_c 1 %, Tu_m < 1.5 % | confirmed |
| 5 | Range | Ma_c 0 / 0.3 / 0.6; Ma_m 0 / 0.3 / 0.6 / 1.2; PR = p_tc/p_m from 1 to 2 | confirmed |
| 6 | **Cd definition** (p. 559, printed) | Cd = m_dot / m_ideal, with `m_ideal = p_tc (p_m/p_tc)^((k+1)/2k) sqrt(2k/((k-1) R T_tc) ((p_tc/p_m)^((k-1)/k) - 1)) pi D^2/4`. This is isentropic 1-D flow from the coolant TOTAL pressure to the mainstream STATIC pressure | confirmed |
| 7 | **Area reference** | pi D^2/4 on the INLET (metering) diameter for all three holes, so a shaped-hole Cd can exceed 1 | confirmed |
| 8 | Pressure taps | p_tc and T_tc are taken 1D upstream of the hole inlet, 2D off the channel centreline. p_m is taken on the side wall 10D upstream of the exit (subsonic) or on the top wall opposite the exit (Ma_m 1.2) | confirmed |
| 9 | **Uncertainty** | Kline-McClintock: flow by orifice 2.5 %. Cd below 2 % in most cases, up to 4.5 % at low PR / low flow | confirmed |
| 10 | Figs. 6/10/14 PR levels | cylindrical 1.11/1.24/1.42/1.63; fan 1.09/1.25/1.41/1.62; laidback 1.10/1.25/1.40/1.60 | confirmed |

**Convention.** Items 6-7 use the same pressure pair as Hay et al. (1983). Both
sources are therefore declared under `Cd_coolant_total_exit_static` in
`validation/cooling/conventions.py`, and its note now records the inlet-area
basis. A combaero `Cd` on the metering area matches this basis. A Cd on any
other area does not.

**Uncertainty.** Item 9 depends on PR, so the schema's scalar `uncertainty` is
null on every series, as for hay1983. The band is recorded here and in the
metadata header. It must never be replaced by a number derived from a model's
error on these series.

## Figures digitised

| Fig. | hole | x | held | curves | marks (confirmed / provisional) |
|---|---|---|---|---|---|
| 3 | cylindrical | PR | Ma_c 0.0 | Ma_m 0, .3, .6, 1.2 | 49 / 0 |
| 4 | cylindrical | PR | Ma_c 0.3 | same | 39 / 0 |
| 5 | cylindrical | PR | Ma_c 0.6 | same | 37 / 0 |
| 6 | cylindrical | Ma_c | Ma_m 0 | PR 1.11, 1.24, 1.42, 1.63 | 44 / 0 |
| 7 | fan-shaped | PR | Ma_c 0.0 | Ma_m 0, .3, .6, 1.2 | 38 / 13 |
| 8 | fan-shaped | PR | Ma_c 0.3 | same | 41 / 4 |
| 9 | fan-shaped | PR | Ma_c 0.6 | same | 27 / 9 |
| 10 | fan-shaped | Ma_c | Ma_m 0 | PR 1.09, 1.25, 1.41, 1.62 | 47 / 0 |
| 11 | laidback | PR | Ma_c 0.0 | Ma_m 0, .3, .6, 1.2 | 34 / 0 |
| 12 | laidback | PR | Ma_c 0.3 | same | 39 / 7 |
| 13 | laidback | PR | Ma_c 0.6 | same | 27 / 2 |
| 14 | laidback | Ma_c | Ma_m 0 | PR 1.10, 1.25, 1.40, 1.60 | 43 / 0 |

Legend symbols: square, circle, triangle and diamond, in the order of each
legend. Figs. 15-18 regroup the same 36 PR curves by Ma_m and were not
digitised separately.

Fig. 8 has two symbols that the first pass left out. A is a square peak near
PR 1.07, Cd 1.15. B is an overlapping square and circle near PR 1.25. The user
confirmed both in a final review. A is filed as a provisional square. B is a
confirmed square and a confirmed circle.

## Digitisation method and calibration evidence

1. Symbols were detected by template matching on the 300 dpi bitmap. The
   markers follow gnuplot conventions: the point is the template origin, and
   for a triangle it lies at 0.69 of the height below the apex.
2. The user reviewed scan overlays in three rounds on 2026-11-08 and confirmed
   every figure. Each decision is in the CSV's `note` column. Every point the
   user did not single out was declared correct and is `confirmed`.
3. The scan and the overlays stay out of git.

Each panel was calibrated on its printed dotted grid lines, because the 1 px
ticks are partly lost. The fit uses two half-bands per axis and an affine map,
which absorbs the 0.1-0.2 deg page skew. All axes are linear.

| Fig. | grid residual max x / y | frame edge through the fit vs printed |
|---|---|---|
| 3 | 0.0013 / 0.0017 | left 0.9987 vs 1 |
| 4 | 0.0006 / 0.0007 | left 1.0007, right 1.9964, top 0.8993 |
| 5 | 0.0011 / 0.0004 | left 0.9989, top 0.8999, bottom 0.1012 |
| 7 | 0.0010 / 0.0008 | left 0.9988, top 0.9004, bottom 0.0998 |
| 8 | 0.0006 / 0.0009 | left 1.0002, right 1.9960 |
| 9 | 0.0006 / 0.0010 | left 1.0002, right 1.9965 |
| 11 | 0.0010 / 0.0016 | left 0.9991, bottom 0.1997 |
| 12 | 0.0007 / 0.0014 | left 1.0001, right 1.9971 |
| 13 | 0.0009 / 0.0007 | left 1.0008, right 1.9970, top 1.0999, bottom 0.5001 |
| 6 | 0.0004 / 0.0002 | left 0.0006, right 0.5992 |
| 10 | 0.0003 / 0.0005 | left 0.0009, right 0.5993, bottom 0.8003 |
| 14 | 0.0002 / 0.0004 | left 0.0007, right 0.5991, top 1.1505, bottom 0.7503 |

- **Edge agreement.** All edges agree with the printed numbers to 1-3 px.
- **Right frame of the PR panels.** In every PR panel the right frame reads
  0.003-0.004 below 2.0. The grid fits are clean, so the offset is in the
  drawn frame, not the scale.
- **Point precision.** +-1 px is +-0.001-0.002 in Cd.
- **Calibration file.** The per-panel fit is
  `data/gritsch1998/calibration_grid.json`. It holds the grid-line pixel
  positions, coefficients, residuals and edge checks. It is committed but not
  listed under `files:`.
- **Figure cards.** In Figs. 9 and 12 the lowest marks lie below the last
  labelled grid line. Their cards therefore bound y with the frame edges read
  through the fit (`y_limits`) instead of the tick span.

`uv run python -m validation.cooling.verify` checks the span and the point
count of every series against its figure card.

## What the data can test

- **A shaped-hole Cd.** The data give fan-shaped and laidback against
  cylindrical holes on the same rig and the same metering area. They also give
  its dependence on coolant and mainstream crossflow up to Ma_m 1.2.
- **Scoring.** combaero has no closure for any of this, so nothing is scored.
  A future closure scored on these data is **cross-source ACCURACY** and must
  be labelled so (docs/VALIDATION_POLICY.md). **Fidelity** would apply only to
  a model built from this paper, such as the authors' 1997 correlation.
- **Tuning.** The data are never a tuner. Matching a rig is the user's
  multiplier choice.

### The gap, stated

`FilmCoolingElement` with hole_type "shaped" applies the plain-hole Cd on the
METERING area. The area basis matches this paper's, but the plain-hole value
does not carry the shaped-hole gain. On these data that under-predicts
shaped-hole flow:

- by about **4-11 %** without crossflow (the fan/cyl and laidback/cyl ratios
  at Ma_c 0, Ma_m 0-0.3);
- by about **20-75 %** with coolant crossflow (Ma_c 0.3-0.6 at low and middle
  PR);
- by more at Ma_m 1.2 near PR 1.2-1.4, where the cylindrical hole's Cd
  collapses.

These figures come from the ratios below, not from a score.

## First look: shaped / cylindrical ratios (a look, not a score)

The values are linearly interpolated between neighbouring digitised marks,
with no extrapolation. Provisional marks are included. Each ratio is between
two holes of the same rig at nominally the same condition. Figs. 6/10/14 pair
legend PRs that differ by up to 0.03. Nothing here is pinned by a test or used
to tune anything.

Summary:

- **No crossflow (Ma_c 0, Ma_m 0).** fan/cyl is 1.04-1.11 and laidback/cyl is
  1.07-1.08. The shaped holes run nearly flat in PR (0.80-0.82). The
  cylindrical hole rises from 0.74 to 0.78.
- **Coolant crossflow (Ma_c 0.3, Ma_m 0).** fan/cyl is 1.43 at PR 1.1 and
  1.19 at PR 1.3. The shaped-hole Cd exceeds 1 at low PR.
- **Ma_c 0.6.** fan/cyl and laidback/cyl are 1.55-1.75 at PR 1.3 and fall to
  1.09-1.19 by PR 1.6-1.8.
- **Figs. 6/10/14 at PR ~1.11.** The ratios grow with Ma_c from 1.12 at
  Ma_c 0 to 1.45 at Ma_c 0.3, because the cylindrical Cd peaks at Ma_c ~0.2
  and then falls while the shaped holes keep rising.
- **Layback.** Its effect on Cd is small: laidback/fan is within a few per
  cent almost everywhere, as the paper states.

## What it does NOT cover

- **Geometry.** Other angles, L/D, area ratios or compound angles; any
  inlet-edge shape other than sharp.
- **Supply flow.** Counter-flowing or perpendicular coolant supply. The supply
  here always coflows.
- **Density ratio.** Tt_c/Tt_m 0.54-1 was tested with no effect at matched
  Ma_m, but no data are plotted.
- **Rows.** Multi-hole rows (single holes only).
- **Figs. 15-18.** These replot the same curves.

## Review log

| date | who | entry |
|---|---|---|
| 2026-11-08 | Claude | Digitised Figs. 3-14. First-pass overlays with numbered suspects were sent for review. |
| 2026-11-08 | user | Review rounds 1-3: suspects resolved, missed marks added, collapsed-curve points kept as provisional. The Fig. 8 A/B symbols were confirmed last (A provisional square; B square and circle, confirmed). **CONFIRMED.** Promoted to `validation/cooling/data/gritsch1998/`, unscored. |

## Appendix: ratio tables (ratios.md, first look)


**Ma_c = 0.0, Ma_m = 0.0** (points used: cyl 11, fan 13, laidback 7)

| PR | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 1.1 | 0.737 | 0.816 | 0.794 | 1.107 | 1.077 |
| 1.2 | 0.745 | 0.821 | 0.799 | 1.101 | 1.072 |
| 1.3 |   -   | 0.816 |   -   |   -   |   -   |
| 1.4 |   -   | 0.817 |   -   |   -   |   -   |
| 1.6 | 0.769 | 0.820 |   -   | 1.066 |   -   |
| 1.8 | 0.784 | 0.814 |   -   | 1.038 |   -   |

**Ma_c = 0.0, Ma_m = 0.3** (points used: cyl 12, fan 12, laidback 7)

| PR | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 1.1 | 0.707 | 0.793 | 0.787 | 1.122 | 1.113 |
| 1.2 | 0.735 | 0.809 | 0.797 | 1.100 | 1.085 |
| 1.3 | 0.744 | 0.811 | 0.800 | 1.090 | 1.075 |
| 1.4 | 0.755 | 0.817 | 0.803 | 1.081 | 1.064 |
| 1.6 | 0.770 | 0.819 | 0.815 | 1.063 | 1.058 |
| 1.8 | 0.784 | 0.814 |   -   | 1.039 |   -   |

**Ma_c = 0.0, Ma_m = 0.6** (points used: cyl 12, fan 13, laidback 7)

| PR | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 1.1 | 0.554 | 0.713 | 0.732 | 1.287 | 1.321 |
| 1.2 | 0.664 | 0.776 | 0.773 | 1.169 | 1.164 |
| 1.3 | 0.699 | 0.796 | 0.787 | 1.138 | 1.127 |
| 1.4 | 0.727 | 0.808 |   -   | 1.112 |   -   |
| 1.6 | 0.755 | 0.815 |   -   | 1.079 |   -   |
| 1.8 |   -   | 0.814 |   -   |   -   |   -   |

**Ma_c = 0.0, Ma_m = 1.2** (points used: cyl 14, fan 13, laidback 13)

| PR | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 1.1 |   -   | 0.421 | 0.467 |   -   |   -   |
| 1.2 | 0.477 | 0.633 | 0.615 | 1.327 | 1.290 |
| 1.3 | 0.581 | 0.715 | 0.688 | 1.231 | 1.184 |
| 1.4 | 0.642 | 0.742 | 0.725 | 1.156 | 1.130 |
| 1.6 | 0.709 | 0.798 | 0.778 | 1.125 | 1.096 |
| 1.8 | 0.753 | 0.815 | 0.808 | 1.082 | 1.074 |

**Ma_c = 0.3, Ma_m = 0.0** (points used: cyl 8, fan 11, laidback 12)

| PR | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 1.1 | 0.790 | 1.126 | 1.117 | 1.426 | 1.414 |
| 1.2 | 0.842 | 1.053 | 1.025 | 1.250 | 1.217 |
| 1.3 | 0.855 | 1.015 | 0.981 | 1.188 | 1.147 |
| 1.4 |   -   | 0.973 | 0.956 |   -   |   -   |
| 1.6 |   -   |   -   |   -   |   -   |   -   |
| 1.8 | 0.868 |   -   | 0.918 |   -   | 1.057 |

**Ma_c = 0.3, Ma_m = 0.3** (points used: cyl 9, fan 10, laidback 10)

| PR | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 1.1 | 0.716 | 1.088 | 1.091 | 1.519 | 1.523 |
| 1.2 | 0.830 | 1.045 | 1.018 | 1.259 | 1.226 |
| 1.3 | 0.845 | 1.013 | 0.986 | 1.198 | 1.167 |
| 1.4 | 0.856 | 0.970 | 0.955 | 1.133 | 1.116 |
| 1.6 |   -   |   -   |   -   |   -   |   -   |
| 1.8 |   -   | 0.905 |   -   |   -   |   -   |

**Ma_c = 0.3, Ma_m = 0.6** (points used: cyl 10, fan 9, laidback 10)

| PR | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 1.1 |   -   | 0.912 | 1.019 |   -   |   -   |
| 1.2 | 0.746 | 0.995 | 0.993 | 1.333 | 1.332 |
| 1.3 | 0.801 | 0.990 | 0.965 | 1.236 | 1.205 |
| 1.4 | 0.827 |   -   | 0.951 |   -   | 1.150 |
| 1.6 | 0.846 |   -   |   -   |   -   |   -   |
| 1.8 |   -   |   -   |   -   |   -   |   -   |

**Ma_c = 0.3, Ma_m = 1.2** (points used: cyl 12, fan 15, laidback 14)

| PR | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 1.1 |   -   |   -   |   -   |   -   |   -   |
| 1.2 | 0.427 | 0.830 | 0.808 | 1.946 | 1.893 |
| 1.3 | 0.618 | 0.897 | 0.872 | 1.452 | 1.412 |
| 1.4 | 0.706 | 0.924 | 0.899 | 1.309 | 1.273 |
| 1.6 | 0.795 | 0.921 | 0.919 | 1.159 | 1.156 |
| 1.8 | 0.825 | 0.911 | 0.919 | 1.104 | 1.113 |

**Ma_c = 0.6, Ma_m = 0.0** (points used: cyl 12, fan 9, laidback 9)

| PR | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 1.1 |   -   |   -   |   -   |   -   |   -   |
| 1.2 |   -   |   -   |   -   |   -   |   -   |
| 1.3 | 0.622 | 0.988 | 0.963 | 1.589 | 1.549 |
| 1.4 | 0.738 | 0.964 | 0.964 | 1.305 | 1.306 |
| 1.6 | 0.800 | 0.914 | 0.929 | 1.142 | 1.160 |
| 1.8 | 0.820 | 0.897 | 0.902 | 1.094 | 1.099 |

**Ma_c = 0.6, Ma_m = 0.3** (points used: cyl 6, fan 7, laidback 5)

| PR | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 1.1 |   -   |   -   |   -   |   -   |   -   |
| 1.2 |   -   |   -   |   -   |   -   |   -   |
| 1.3 | 0.565 | 0.971 | 0.991 | 1.718 | 1.754 |
| 1.4 | 0.726 | 0.956 | 0.963 | 1.316 | 1.326 |
| 1.6 | 0.795 | 0.914 |   -   | 1.150 |   -   |
| 1.8 |   -   |   -   |   -   |   -   |   -   |

**Ma_c = 0.6, Ma_m = 0.6** (points used: cyl 7, fan 9, laidback 8)

| PR | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 1.1 |   -   |   -   |   -   |   -   |   -   |
| 1.2 |   -   |   -   |   -   |   -   |   -   |
| 1.3 |   -   | 0.535 | 0.973 |   -   |   -   |
| 1.4 | 0.653 | 0.947 | 0.950 | 1.450 | 1.454 |
| 1.6 | 0.768 | 0.914 |   -   | 1.191 |   -   |
| 1.8 | 0.808 |   -   |   -   |   -   |   -   |

**Ma_c = 0.6, Ma_m = 1.2** (points used: cyl 12, fan 11, laidback 7)

| PR | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 1.1 |   -   |   -   |   -   |   -   |   -   |
| 1.2 |   -   |   -   |   -   |   -   |   -   |
| 1.3 |   -   |   -   |   -   |   -   |   -   |
| 1.4 | 0.242 | 0.592 | 0.616 | 2.447 | 2.546 |
| 1.6 | 0.633 | 0.891 | 0.893 | 1.407 | 1.411 |
| 1.8 | 0.743 | 0.888 | 0.896 | 1.195 | 1.206 |

### Figs. 6 / 10 / 14 (Ma_m = 0): Cd vs Ma_c at fixed PR


**PR ~ 1.11** (cyl 1.11 / fan 1.09 / laidback 1.10; points 13/15/8)

| Ma_c | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 0.0 | 0.733 | 0.819 |   -   | 1.118 |   -   |
| 0.1 | 0.826 | 0.965 | 0.935 | 1.168 | 1.132 |
| 0.2 | 0.848 | 1.098 | 1.043 | 1.295 | 1.229 |
| 0.3 | 0.781 | 1.136 | 1.123 | 1.454 | 1.438 |
| 0.4 |   -   |   -   | 0.967 |   -   |   -   |
| 0.5 |   -   |   -   |   -   |   -   |   -   |

**PR ~ 1.24** (cyl 1.24 / fan 1.25 / laidback 1.25; points 10/12/12)

| Ma_c | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 0.0 |   -   | 0.824 |   -   |   -   |   -   |
| 0.1 | 0.810 | 0.914 | 0.866 | 1.128 | 1.069 |
| 0.2 |   -   | 0.991 | 0.952 |   -   |   -   |
| 0.3 | 0.844 | 1.041 | 1.000 | 1.234 | 1.186 |
| 0.4 | 0.799 | 1.057 | 1.034 | 1.323 | 1.295 |
| 0.5 | 0.695 | 1.022 | 1.028 | 1.470 | 1.480 |

**PR ~ 1.42** (cyl 1.42 / fan 1.41 / laidback 1.40; points 8/10/10)

| Ma_c | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 0.0 | 0.757 | 0.810 |   -   | 1.070 |   -   |
| 0.1 |   -   | 0.902 |   -   |   -   |   -   |
| 0.2 |   -   | 0.946 |   -   |   -   |   -   |
| 0.3 | 0.853 | 0.967 | 0.958 | 1.133 | 1.123 |
| 0.4 | 0.842 | 0.979 | 0.980 | 1.164 | 1.165 |
| 0.5 | 0.803 | 0.977 | 0.979 | 1.217 | 1.220 |

**PR ~ 1.63** (cyl 1.63 / fan 1.62 / laidback 1.60; points 13/10/13)

| Ma_c | Cd cyl | Cd fan | Cd laidback | fan/cyl | laidback/cyl |
|---|---|---|---|---|---|
| 0.0 | 0.770 |   -   | 0.812 |   -   | 1.055 |
| 0.1 | 0.826 |   -   | 0.859 |   -   | 1.040 |
| 0.2 | 0.854 | 0.901 | 0.912 | 1.055 | 1.068 |
| 0.3 | 0.863 | 0.920 | 0.936 | 1.065 | 1.084 |
| 0.4 | 0.858 | 0.925 | 0.937 | 1.078 | 1.091 |
| 0.5 | 0.838 | 0.916 | 0.936 | 1.093 | 1.117 |
