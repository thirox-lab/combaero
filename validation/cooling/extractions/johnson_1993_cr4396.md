# NASA CR-4396 Vol. II (Johnson, Wagner and Steuber, 1993) -- stationary trips

**Status: digitised, user-checked and scored, 2026-10-22 (#452).**
Source: `sources/heat_transfer/rib/CR4396/19940011518.pdf` (gitignored), contract
NAS3-23691, Pratt & Whitney / UTRC. Data: `data/johnson1993/`.

## Why this source matters

Every other rib heat-transfer source in this dataset traces to Texas A&M
(Han's or Lau's laboratory; #403 item 2). This is an industrial rig with
cast-style trips, so it is an **independent-lab accuracy check** on the Han
sets. It is never fidelity, and the harness does not tune to it.

**Circularity gate.** CR-4396 is not among the data any shipped rib set was
fitted on. Its Fig. 10.1 *compares* test 301 with Han et al. (1986). That is
a comparison, not pooling.

## What the report prints

| item | value | where |
|---|---|---|
| passage | square, W/H = 1 (H/w = 1.0) | p. 21 |
| trips | e/d = 0.1, P/e = 10, leading + trailing walls only, **staggered** | p. 9, 21 |
| trip shape | curvilinear (rounded, cast-style) | p. 9, 21 |
| skewed model | 45 deg, both walls parallel in projection (Fig. 3.5) | p. 11, 21 |
| thermal BC | all four walls heated, constant wall temperature, inlet (Tw-Tb)/Tw = 0.13 | p. 18 |
| properties | film temperature (Tw+Tb)/2 | p. 18 |
| area | projected; true area 1.11x projected (normal), 1.15x (skewed) | p. 18 |
| reference | Nu_inf = 0.0176 Re^0.8, Kays and Perkins (1973), constant wall T, "Prandtl number equal to 0.72" | p. 18 |
| uncertainty | Kline-McClintock: +-6% at the inlet to +-30% at the exit (stationary); repeats agree within about half that | p. 18-19 |

**Fig. 3.5's arrangement.** The leading-wall trips (dashed) and the
trailing-wall trips (solid) run the same way in the cross-sectional view,
interleaved. That is Han's "parallel" arrangement with a stagger, not a
"crossed" one.

## The conversion, and why it is not optional

The figures plot Nu/Nu_inf on **Kays and Perkins**. The shipped sets' Nu
ratio is the ribbed-wall Nu over **Dittus-Boelter** 0.023 Re^0.8 Pr^0.4.
Both references scale as Re^0.8, so the conversion is one constant:

    Nu/Nu_DB = (Nu/Nu_KP) * 0.0176 / (0.023 * 0.72^0.4) = 0.873 (Nu/Nu_KP)

The Pr in it is the source's own 0.72, because 0.0176 *is* Kays and Perkins
at 0.72. This is definition, not a fit. It is registered in
`conventions.py` as `("Nu_ratio_ribbed_KP", "Nu_ratio_ribbed_DB")`.

**Skipping it hides the result.** On the raw Kays-Perkins basis the Han sets
read -2.2% (normal) and +3.3% (skewed), which looks like agreement. On a
common basis they over-predict by 12-18% (below). This is the category
error #389 names, and the same size as the ~12% gap #444 criterion 6
warned about.

A second-order assumption sits under every Nu-ratio comparison: the ribbed
and smooth Nu share the same Pr^0.4. Here the ratio is taken at Pr 0.72
against a set evaluated at its own 0.70. If they did not share it, the
difference would be 1.1%.

## Digitisation

400 dpi render of PDF pages 48-49. Scope is the first leg, segments B-D,
both walls, four Re: 48 points.

- **Abscissa.** Each symbol's nominal Re from the legend: square 12,500
  (test 302), circle 25,000 (301), triangle 50,000 (303), diamond 75,000
  (304). The figure's own x is S/d, which identifies the segment.
- **Calibration.** The y axis is log, calibrated per panel on the printed
  ticks. The worst tick residual is 1.5 px, or 0.4% in Nu/Nu_inf.
- **Squares** come from the flat top edge (+19.5 px to centre). **Diamonds**
  come from the bottom tip (-30 px to centre). Both are cross-checked
  against an outline-template match.
- **The diamond offset was measured, not assumed.** One plot diamond has a
  clean upper outline (Fig. 6.2, leading wall, segment R). Its top tip is
  row 1597 and its widest row is 1627, so the centre is tip + 30, and the
  lower half mirrors the upper within 1 px. 1 px is 0.26% in Nu/Nu_inf.
  The legend diamonds are drawn about 15% larger than the plot's and are
  not used.
- **Circles and triangles** come from the template match only. A
  triangle's position is its bounding-box centre. If the plotter centred
  triangles on the centroid, they read about 1.8% high.
- **User review.** Nine points were flagged and checked by the user on a
  scan overlay, all confirmed: two template/anchor disagreements and seven
  weak template matches, mostly Fig. 6.2 triangles inside other outlines.
- **Overplotting.** Four symbols share a cluster about 60 px tall, so the
  score is limited by the read. The edge-anchored half (24 points) is the
  best-located.

## Score (cross-source accuracy, on the Dittus-Boelter basis)

Every eligible rib set scores every series, each on its own row
(`also_scores`, docs/VALIDATION_POLICY.md). Eligibility is derived by rule
in `runner.eligible_sets` and pinned by
`test_cr4396_is_scored_by_every_eligible_set_and_no_fewer`, so the list was
not chosen after seeing the numbers. Sets are listed by name, not by score.

| figure | set | n | in its box | bias (all) | bias (in box) |
|---|---|---|---|---|---|
| 6.1 normal | `han_1988_orthogonal` | 24 | 0 | +12.1% | - |
| 6.1 normal | `han_zhang_lee_1991_90` | 24 | 0 | +5.2% | - |
| 6.2 skewed 45 deg | `han_1989_narrow_channel` | 24 | 0 | +7.1% | - |
| 6.2 skewed 45 deg | `han_park_1988_angled` | 24 | 0 | +18.4% | - |
| 6.2 skewed 45 deg | `han_zhang_lee_1991_45par` | 24 | 0 | -1.6% | - |
| 6.2 skewed 45 deg | `rallabandi_2009_high_re` | 24 | 12 | +5.0% | -0.6% |

Only Rallabandi's box (e/D 0.1-0.18, 45 deg, Re >= 30k) covers any of the
rig; its 12 in-box points are Re 50k and 75k (MAE 3.9%). Every other set
is outside on e/D (0.1 against 0.047-0.078, or HZL's single tested
0.0625). The spread across sets is the number a designer inherits by
picking one: 7 points on the normal trips, 20 on the skewed.

**Not independent of each other.** All six sets are Han's laboratory at
Texas A&M (#403 item 2). That they bracket CR-4396 is not six confirmations;
CR-4396 is the only independent lab here.

For the two original sets, by wall: 6.1 leading +16.4%, trailing +7.7%;
6.2 leading +17.0%, trailing +19.7%. By Re the bias is flat: 6.1 runs +9.8
to +15.0%, 6.2 runs +15.2 to +20.3%. So the miss is a level offset, not a
Re trend.

## Not modelled

Staggered rather than in-line trips, rounded trip profiles, constant wall
temperature, heating on all four walls, and the extra trip area (1.11x and
1.15x the projected area) are all left out. Any of them could carry part of
the offset. None is fitted for.
