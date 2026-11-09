# Naik, Retzko, Gritsch and Sedlov (2014), GT2014-25841 -- wedge trailing-edge passage ribs

**Status: digitised, user-approved on scan overlays (2026-11-09), committed and
scored as cross-source accuracy.**
Source: `sources/heat_transfer/rib/Naik_Retzko_Gritsch_Sedlov_2014_GT2014-25841_TE_wedge_turbulators.pdf`
(gitignored; sha256 `04e8ae90...125af`, see `validation/sources.yaml`). Alstom,
Baden and Moscow. Data: `data/naik2014/`.

## For the reviewer

1. **Which walls Fig. 16 averages is not stated.** It is declared ribbed-wall
   (`Nu_ratio_ribbed_DB`). The ribbed PS and SS walls are about 80% of the
   Fig. 11c perimeter. A by-eye read of the measured line average there gives
   the LE wall ~1.5 Nu0, the TE apex ~2.2 and the ribbed walls ~2.25, so a
   whole-perimeter mean would sit ~2-5% below the ribbed-wall mean. If
   Fig. 16 is the perimeter mean, every prediction reads that much high
   against a ribbed-wall truth. That is inside the 9.5% band but not nothing.
   It is recorded, not corrected. A digitisation of Fig. 11c would tighten
   it, but would not settle the definition.
2. **The friction is not scored, by the existing rule.** Only
   `han_1988_orthogonal` declares the passage-average friction
   (`f_channel_two_ribbed`), and it binds transverse 90 deg ribs. So
   `eligible_sets` is empty for Fig. 17. Declaring it for
   `han_park_1988_angled` is a one-line change to `SET_ALSO_PRODUCES`. It
   would also make that set eligible for all 8 two-sided Taslim & Spring
   friction series, which would gain rows. What it would read on Naik is
   below ("Friction, off the harness"). This is your call; the commit does
   not make it.
3. **Uncertainty is the stated 9.5% / 10%, not the drawn +-15%.** Fig. 16
   draws +-15% bars on every point (read: 14.9-15.0%). The text gives HTC
   9.5% and f 10% (Moffat, p. 5) and does not say what the bars are. The
   stated figure is used.

Nothing further is needed from the source.

## What the paper prints

| item | value | where |
|---|---|---|
| passage | wedge (triangular) trailing-edge passage, AR = W/H = 6.64, H at the thick LE end | p. 2, Table 1 |
| ribbed walls | pressure and suction side (the two wide walls); LE wall smooth | Figs. 2-4 |
| angled-tapered | 45 deg, PS/SS staggered p/2; e / local H constant; P/e 15 avg (8 at LE end); e/Dh 0.089 avg (0.163 max) | p. 3, Table 1 |
| broken | 45 deg discrete, two rows; LE half P/e 15, e/Dh 0.134; TE half P/e 30, e/Dh 0.066, half width | p. 3, Table 1 |
| chevron | small 45 deg chevrons, two rows, same P/e and e/Dh split | p. 3, Table 1 |
| rib cross-section | not stated | -- |
| Re | 22k, 60k, 145k, on the passage Dh (Dh, W, H not given numerically) | Table 1 |
| method | transient TLC on Perspex, over the last 3-4 pitches before the outlet; 12 thermocouples, 4 pressure taps | p. 4 |
| fully developed | angled: equal line averages at two stations | p. 7 |
| averaging | rib-covered area excluded | pp. 6-7 |
| Nu0 | 0.023 Re^0.8 Pr^0.4 (Eq. 3); Pr not printed (air, 0.70 assumed) | p. 7 |
| f | dp / (4 (dL/dh) (1/2) (m/A)^2/rho), FANNING, passage average (Eq. 4) | p. 7 |
| f0 | 2 [2.236 ln Re - 4.639]^-2, Petukhov in Fanning form (Eq. 5) | p. 7 |
| PI | (Nu/Nu0) / (f/f0)^(1/3) (Eq. 6) | p. 10 |
| uncertainty | HTC 9.5%, f 10% (mass flow 2%, properties 2%) | p. 5 |

**Eq. 5 is Fanning.** 2.236/sqrt(2) = 1.5811 and 4.639/sqrt(2) = 3.2803, so
Eq. 5 is Petukhov's Darcy (0.790 ln Re - 1.64)^-2 divided by 4, to 0.2%
(`test_naik_petukhov_is_fanning`). Eq. 4 is Fanning too (the 4 dL/dh). A
Darcy reading of either would be 4x off.

**Rib arrangement: parallel.** The angled ribs set up a pair of
counter-rotating vortices. The near-wall flow runs TE to LE on both walls,
symmetric (Fig. 10). That is Han's "parallel" arrangement with a stagger, as
in CR-4396. The paper compares itself to Wright et al.'s "Parallel Angled".
So `rib_shape: parallel`.

## Digitisation

First pass 2026-11-09. Scan overlays were reviewed and approved by the user
("naik digitization looks solid").

- **Fig. 16 is a raster.** y is calibrated on the 7 gridlines 1.0-4.0, worst
  residual 0.0021 in Nu/Nu0. Markers are located by pixel extent. At 60k and
  145k the angled circle and the chevron triangle overplot. Three of the five
  marks there are located from the visible edge plus the error-bar cap
  midpoint, which agree to 1 px. **1 of the 3 angled-tapered marks is
  clean** (22k).
- **Figs. 17 and 18 are vector.** Each marker is its own path object with
  its own colour. The centre is the path bounding box, mapped by a fit on
  every labelled tick: y residual <= 0.0086 in f/f0 and <= 0.00075 in PI,
  x residual <= 47 in Re. Nothing can hide a mark. In Fig. 17, Rhee 45 deg V
  and Wright W coincide at (1e4, 9.0); both were read from separate paths.
- **x is the nominal Re** of Table 1, or of the reprints' legends.

| series | nominal Re | read Re (F16 / F17 / F18) | Nu/Nu0 | f/f0 | PI (F18) | F16 clean |
|---|---|---|---|---|---|---|
| angled_tapered | 22,000 | 22,120 / 21,980 / 21,997 | 2.7050 | 7.2237 | 1.1253 | yes |
| angled_tapered | 60,000 | 59,943 / 60,015 / 59,971 | 2.1485 | 7.9032 | 0.8168 | no |
| angled_tapered | 145,000 | 145,043 / 144,993 / 145,026 | 2.1427 | 10.8889 | 0.5900 | no |
| broken | 22,000 | as above | 2.5820 | 11.1188 | 0.6972 | yes |
| broken | 60,000 | | 1.8791 | 12.3914 | 0.4555 | yes |
| broken | 145,000 | | 1.9201 | 15.8749 | 0.3639 | yes |
| chevron | 22,000 | -- / 21,980 / -- | (none) | 14.7267 | (none) | -- |
| chevron | 60,000 | | 2.1017 | 18.2581 | 0.3458 | no |
| chevron | 145,000 | | 2.2013 | 23.3111 | 0.2831 | no |

The chevron has **no 22k mark in Figs. 16 or 18**: no marker and no third
error bar. Fig. 17 has one. This is a gap in the figure, not in the read.
The reprints' 24 marks (Rhee 45/60 deg V at 10/20/30k; Wright parallel and
discrete W at 10/20/40k) are all clean, and in the CSVs.

## Fig. 18 is mis-computed in the paper

Fig. 18 is the check that could have failed. It is a third figure, read
independently, and Eq. 6 ties it to Figs. 16 and 17.

| group | n | PI against Eq. 6 | PI against 3 (Nu/Nu0)/(f/f0) |
|---|---|---|---|
| reprints (Rhee 2003, Wright 2004) | 12 | 0.13% worst | 14-62% off |
| the paper's own (angled, broken, chevron) | 8 | 24% to 172% low | 0.3% worst; factor 2.998-3.009 |

So the digitisation of all three figures is confirmed, to 0.13% on the
reprints and 0.3% on the paper's own points. **The paper's own Fig. 18
column is 3 (Nu/Nu0)/(f/f0), not Eq. 6.** That is what A/(B^1/3) evaluates
to when `^` binds before `/`: A / ((B^1)/3) = 3A/B. The figure's own curves,
and the abstract's claim that the angled ribs give the best PI, are on that
quantity. The ranking survives on Eq. 6: angled 1.40/1.08/0.97, broken
1.16/0.81/0.76, chevron -/0.80/0.77. The magnitudes do not. Fig. 18 is
committed as this cross-check and is **never scored**. Pinned by
`test_naik_fig18_own_points_are_3_nu_over_f` and
`test_naik_fig18_reprints_follow_eq6`.

## What is committed

- **Fig. 16 angled-tapered Nu/Nu0**: scored (below).
- **Fig. 16 broken and chevron**: unscored context. No shipped set models
  discrete ribs in two staggered rows of different height and pitch
  (`rib_shape: discrete_angled`, `discrete_chevron`). Han, Zhang and Lee's V
  ribs span the wall.
- **Fig. 17, all series**: unscored (reviewer item 2 for the angled case).
- **Fig. 18, all series**: the cross-check above. Never scored.
- **Reprints**: `reprint_*`, with the primary named in `after`, geometry
  null, unscored. They are second-hand, and Naik gives neither their e/Dh
  nor their P/e, so the Re_Dh path refuses them. They carry the Fig. 18
  proof.

## Score: cross-source accuracy, out of domain, an upper bound

Every eligible set scores the series, each on its own row. Eligibility is
derived by rule (`runner.eligible_sets`): 45 deg parallel ribs, a Nu ratio
already on the Dittus-Boelter basis (resolves direct), a prediction
returned. It is pinned by
`test_naik2014_nu_is_scored_by_every_eligible_set_as_accuracy`. Sets are
listed by name, not by score. All four are Texas A&M (Han), so their
spread is not four confirmations. Naik is the independent lab.

| set | 22k | 60k | 145k | bias | within 9.5% |
|---|---|---|---|---|---|
| `han_1989_narrow_channel` | 3.853 (+42.4%) | 3.342 (+55.5%) | 2.943 (+37.3%) | +45.1% | 0 of 3 |
| `han_park_1988_angled` | 2.268 (-16.2%) | 1.962 (-8.7%) | 1.726 (-19.5%) | -14.8% | 1 of 3 |
| `han_zhang_lee_1991_45par` | 3.604 (+33.2%) | 3.089 (+43.8%) | 2.693 (+25.7%) | +34.2% | 0 of 3 |
| `rallabandi_2009_high_re` | 3.525 (+30.3%) | 2.857 (+33.0%) | 2.367 (+10.5%) | +24.6% | 0 of 3 |
| measured Nu/Nu0 | 2.705 | 2.149 | 2.143 | | |

**Every point is outside every set's box**: W/H 6.64 against fitted 1-4
(han_park), 0.25-1 (han_1989) or 1 (HZL, Rallabandi); e/Dh 0.089 against
<= 0.078 (Rallabandi: >= 0.10); Re 145k above 60k/90k. It is also a wedge,
and the ribs are tapered. **Sampling is `unknown`, and 1 of the 3 marks is
clean**, so these rows are an upper bound on the error, not an estimate.

**What drives `han_park_1988_angled`'s number.** Its G is 2.24 (W/H)^0.1
(e+)^0.35 for any rectangular channel. Its R takes W/H capped at 2 (Eq.
4.17). Off the harness, its Nu ratio at e/Dh 0.089 and P/e 15 reads -16/-7/
-18% at W/H 2, -15/-7/-18% at W/H 4 and -16/-9/-20% at W/H 6.64. The miss
is the same at W/H 2 as at 6.64. Above the square-to-rectangular switch,
W/H moves the prediction by under 2%, so the Nu miss is not a W/H
extrapolation effect. At W/H 1 (the square branch) the same inputs read
+8/+20/+8%. The measured ratio is flat from 60k to 145k (2.15 to 2.14); every
set falls with Re.

**Friction, off the harness** (reviewer item 2; not a score). The measured
fbar is (f/f0) f0(Re) = 0.0460 / 0.0397 / 0.0453 at 22k / 60k / 145k. If
`han_park_1988_angled` declared the passage friction, its fbar at W/H 6.64
would be 0.0653 / 0.0651 / 0.0650, which is +42 / +64 / +44% (mean +50%).
Unlike Nu, this grows with W/H through Han's four-sided decomposition
`fbar = (f W/H + f_s)/(W/H + 1)`, which weights the ribbed walls by W/H. The
same set and inputs read -11 / +2 / -11% at W/H 2 and +22 / +41 / +23% at
W/H 4. With e/Dh held at its in-box 0.078, W/H 6.64 still reads +26 / +45 /
+27%.

## What it can and cannot settle

- **Can**: how far the shipped rectangular-channel rib sets sit from one
  real high-AR wedge passage, measured by an independent lab, on heat
  transfer. It also gives evidence for an envelope decision on W/H > 4.
  For Nu the evidence is that W/H above 2 barely moves the prediction, and
  the miss is a level. For friction (off the harness) it is that the
  decomposition drives a growing over-prediction.
- **Cannot**: fidelity. There is no correlation of Naik's. It cannot
  separate W/H from the wedge shape, the taper or e/Dh: one geometry varies
  all four at once. It cannot tune anything. With 3 points, 1 of them
  clean, it cannot measure a Re trend. It cannot validate broken or chevron
  ribs, which no shipped set models.

## Not modelled

The wedge section (Dh, secondary-flow cells and the smooth LE wall), tapered
ribs, the inter-rib-only average (rib surfaces excluded, where Han's copper
segments include the rib heat), TLC on Perspex against Han's heated copper,
and an unstated rib cross-section. Any of them could carry part of the
offset. None is fitted for.

## Review log

- 2026-11-09: screened. Registered in `validation/sources.yaml` (PR #1078).
- 2026-11-09: first pass. Fig. 16 raster calibration and marker extents;
  Figs. 17-18 vector paths. Cross-check found the Fig. 18 error (reprints on
  Eq. 6 to 0.13%; own points on 3 Nu/f to 0.3%).
- 2026-11-09: user reviewed the scan overlays and approved them ("naik
  digitization looks solid").
- 2026-11-09: promoted to `data/naik2014/` and scored. Fig. 16 angled-tapered
  goes against the four eligible sets, all accuracy and out of domain.
  Friction is unscored by rule, and the decision is left to the reviewer.
