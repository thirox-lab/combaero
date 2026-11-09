# Rib friction and heat transfer across aspect ratios -- Taslim & Spring 1987

Taslim, M.E. and Spring, S.D. (1987). "Friction Factors and Heat Transfer
Coefficients in Turbulated Cooling Passages of Different Aspect Ratios,
Part I: Experimental Results." AIAA-87-2009, 23rd Joint Propulsion
Conference, San Diego. DOI 10.2514/6.1987-2009.
`sources/heat_transfer/rib/Taslim_Spring_1987_AIAA-87-2009_turbulated_passages_aspect_ratio.pdf`

## Why this source

**It is the first rib source in this dataset that is not Texas A&M.**
Northeastern (Taslim) with General Electric (Spring). `han2012`,
`han_park_lei1984` and `lau1990` all trace to Han's or Lau's laboratory,
and #403 item 2 is blocked for exactly that reason: no TAMU source can
arbitrate a TAMU disagreement. #401's "16 configurations across three
rigs" has the same limitation, now recorded.

It also covers **Han's narrow branch**: Taslim's AR 3.5 is Han's
W/H = 0.286, inside the `1/4 <= W/H < 1/2` range of Eq. 4.19 that #402
will implement from Han, Ou, Park & Lei (1989). That set would otherwise
ship with TAMU-only validation.

## Scope

Square-sectioned rig, turbulated sidewall width `B` = 3.0 in fixed, the
non-turbulated wall height `A` varied. **90 degree transverse, in-line**
turbulators on two opposite walls, and separately on one wall.
`e` = .25/.50/.75 in with `e/s` = 0.10 throughout, so `P/e` = 10.
Re 20,000 to 190,000.

## THE ASPECT RATIO IS INVERTED RELATIVE TO HAN

Taslim's `AR = A/B` is HEIGHT/WIDTH; Han's `W/H` is WIDTH/HEIGHT.

| Taslim AR | `D_H` in | Han `W/H` | `e/D_H` tested |
|---|---|---|---|
| 0.5 | 2.000 | 2.0 | .125 .250 |
| 1.0 | 3.000 | 1.0 | .083 .167 .250 |
| 3.5 | 4.667 | **0.286** | .053 .107 .161 |

Confusing them silently relabels a narrow channel as a wide one. Both are
recorded in every `geometry` block.

**`e/D_H` is COUPLED to AR**, because `e` is fixed at three sizes and
`D_H` follows from the aspect ratio. Low blockage therefore occurs only at
high AR -- which is why **no Taslim configuration reaches Han's
0.047-0.078 band at a square channel**: AR 1.0 bottoms out at 0.083. That
is a property of the rig, not of this paper, and it means every
independent arbitration of #403 item 2 will carry an extrapolation cost.

## Verification -- five figures, each cross-checking another

### Calibration

| figure | x | y | note |
|---|---|---|---|
| 4 | 0.78% | 0.23% | best of the set |
| 5 | 1.06% | 0.20% | |
| 9 | 0.60% | 0.57% | |
| 11 | 1.20% | 0.20% | log-log, x anchored on ticks 1 and 20 |
| 12 | 1.80% | 0.28% | linear |

Worst tick error as a fraction of the axis span. Every x axis shows the
same bow -- near zero at the anchors, peaking mid-range -- so it is the
scan, not the digitiser. Harmless on fig11, where only `y` is read
because the paper's own result is that `f` is Re-independent.

**One of these numbers was my error, not the data's.** The fig4 x check
first read 17.4% because I compared against fig5's tick list; fig4's axis
has **no `5` tick**. Corrected it is 0.78%, the best in the set.

### The smooth-duct curve is a free check, and it passes twice

Figures 4 and 5 each draw `Nu_s = 0.023 Re^0.8 Pr^0.4`. Digitised and
compared against that formula at Pr = 0.70:

| | bias | MAE | slope |
|---|---|---|---|
| fig4 | **+0.30%** | 0.64% | 0.816 |
| fig5 | -1.77% | 1.77% | 0.802 |

Expected slope 0.8. This confirms the ordinate is `Nu x 1e-3` without
being told, and validates both calibrations against no correlation at all.

### `Nu_T ~ Re^0.6` recovered from the marks

    fig5 series   0.602  0.647  0.580  0.606
    fig4 series   0.518  0.616  0.594
    fig4's own drawn fit line   0.613   (printed as Re^0.6)

**The 0.518 is real and is not one bad mark**: leave-one-out over the
`e/D` 0.053 series gives 0.511 to 0.528. It is also the series most buried
in fig4's symbol clusters. Recorded as a disagreement with the paper's
stated exponent, not corrected -- the author's own fit line covers all
three series jointly at 0.613, so the paper treats 0.6 as applying here.

### Figures 4+5 reconstruct Figure 9 to 5.4%

Figure 9's ordinate is `(Nu_T/Nu_s)(Re_DH/Re_ref)^0.2`, `Re_ref = 1e4` --
the 0.2 cancels the residual Re dependence, since `Nu_T ~ Re^0.6` and
`Nu_s ~ Re^0.8`. Rebuilding it from the fig4 and fig5 marks reproduces all
seven of fig9's points at **5.4% MAE, -0.4% bias**. The naive `Nu_T/Nu_s`
reads 42% low, so the normalisation is confirmed rather than assumed.

**The exponent is 0.2 and the figure appears to say 2.** The superscript's
decimal point does not survive the scan. The body text states ".2" twice
and explains why. Exponent 2 would put the ordinate at 35 to 523 on an
axis that stops at 5. Same class as the lost decimal in Han's
`(e/D)^0.014`, which was worth 22% -- see `han_ribbed_high_re.md` item 31.

### Figure 12 misdraws one marker, and it is the most valuable one

Six of fig12's seven marks land on their declared `e/D` to within the x
bow (+0.002 to +0.004). The seventh, **AR 1.0 at `e/D` 0.083**, is drawn
at x = 0.1306 -- off by +0.047, thirteen times the bow -- and its `f`
reads 0.0746 against fig11's **0.0468** from eight marks. Wrong in BOTH
coordinates. Confirmed by rendering at 450 dpi: there is no marker
anywhere near `e/D` 0.083 on that curve.

Fig12's own **fit curve** at 0.083 gives 0.0437, within 7% of fig11's
value. So the curve is right and the symbol is misplaced.

This matters more than the other six because `e/D` 0.083 at a square
channel is the closest any independent measurement gets to Han's
0.047-0.078.

### Use fig11 over fig12 wherever they overlap

fig11 gives 4 to 11 marks per series on a two-decade log axis; fig12 gives
one mark per configuration on a linear 0-0.4 axis that cannot resolve
small `f`. At AR 3.5 the fig12 CURVE disagrees with fig11 by up to 41% --
which is 0.0044 absolute, about 1% of fig12's full scale. The marks agree
far better, 1.2% to 13%.

**What fig12 cannot show:** AR 0.5 at `e/D` 0.250 has `f` = 0.5165
(fig11), above fig12's 0.4 frame top. Its curve exits the top at x = 0.233
and the marker is genuinely absent -- independently confirming the
off-scale reading.

### Sealed marker counts

Following the Rohde Fig. 10 precedent, fig4's per-series counts were
written to `SEALED_fig4_marker_counts.txt` before the manual pass and not
quoted in the brief. Sealed **6/7/6**; the manual pass found **7/7/7**.
Both undercounts fall in the three symbol clusters (Re ~ 5.4, 14.4, 17.5)
the seal had flagged LOW confidence. The known failure mode of an
automated centroid read in dense overlap, quantified rather than assumed.

A second seal (`SEALED_fig12_predicts_fig11.txt`) held fig12's curve
predictions for the seven two-side series before fig11 was digitised; they
came in at 4-7% on the AR 0.5 and 1.0 curves and 7-41% on AR 3.5, which is
how the fig12 resolution limit was found.

**One sealed claim was wrong, and it was mine.** The fig4 seal asserted
`.161 > .107 > .053` at both ends of the Re range and said a crossing
would indicate a mis-assigned symbol. The data crosses: at low Re the
order is `.161 > .053 > .107`. That is what the figure shows, and it is
physical -- conclusion 3 says Nu is only *slightly* sensitive to `e/D` at
AR 3.5, so near-coincident series can cross in the scatter. The seal
assumed monotonicity it had no right to.

## What is here

| AR = Han W/H | `e/D` | `f` | `Nu` |
|---|---|---|---|
| 0.5 = 2.0 | .125, .250 | yes | yes |
| 1.0 = 1.0 | **.083**, .167 | yes | yes |
| 1.0 = 1.0 | .250 | yes | **not measured** |
| 3.5 = **0.286** | .053, .107, .161 | yes | yes |

Plus the one-side-turbulated counterpart of all eight friction series.

**Seven configurations carry both; the eighth carries only friction.**
AR 1.0 at `e/D` 0.250 appears in no Nusselt figure -- not 4, not 5, and
not 9, whose AR 1.0 fit curve stops at `e/D` 0.215. That is the study, not
a gap in the digitisation.

## Square-channel friction SCORED (2026-10-02, #403 item 2)

The three **AR 1.0, two-side-turbulated** friction series (Fig. 11) now
score against `han_1988_orthogonal`. At AR 1.0 the aspect-ratio term inside
the logarithm below is exactly 1, and the comparison is made in `f`, not
`R`: Taslim's passage-average Fanning `f` is Han's measured `fbar` (same
definition, Han's Eq. 1), so the set's four-sided `f` is converted with
Han's own `fbar = (f W/H + f_s)/(W/H + 1)` and nothing on the measured side
passes through the law-of-the-wall formalism. Fig. 11's x axis is printed
as `Reynolds No. (x10^-4)`; every Fig. 11 series now carries
`x_scale: 1.0e4`.

| `e/D_H` | n | Han `fbar` vs measured | Taslim `R` on Han's basis |
|---|---|---|---|
| 0.083 | 8 | **-14.9%** | **2.77** (Han 3.20) |
| 0.167 | 8 | -24.9% | 2.73 |
| 0.250 | 7 | -49.1% | 2.50 |

All rows are out of domain: every `e/D` here is above Han's 0.047-0.078
(0.083 by 6%). `R` is flat from 0.083 to 0.167 and falls at 0.250 -- Han's
`e/D`-independence holds at moderate blockage and fails at high, as
Rallabandi found. **What it settles:** Lau's `R` sits ~12% ABOVE Han's,
Taslim's ~13% BELOW, so the independent lab does not side with either --
see `lau_1990_v_ribs.md` L2. The rib cross-section shape is not stated in
the paper, and nor is a measurement uncertainty.

## Every aspect ratio SCORED (2026-10-22)

The AR 0.5 and 3.5 two-side friction series had stayed unscored under the
note below, which #443 had already made stale: the comparison is made in
`f`, so the aspect-ratio term inside Han's `R` never touches the measured
side. W/H enters only through the set's own four-sided `f` (evaluated at the
rig's W/H) and Han's `fbar` inverse, the same path as AR 1.0. By the
CR-4396 rule they are scored by every eligible set -- here only
`han_1988_orthogonal`. Leaving them out kept only the square channel, which
reads best.

| Taslim AR (Han W/H) | `e/D_H` | n | Han `fbar` vs measured |
|---|---|---|---|
| 0.5 (2.0) | 0.250 | 4 | -11.7% |
| 0.5 (2.0) | 0.125 | 11 | +12.6% |
| 3.5 (0.286) | 0.053 | 8 | -34.6% |
| 3.5 (0.286) | 0.107 | 7 | -34.3% |
| 3.5 (0.286) | 0.161 | 8 | -43.7% |

All out of domain: AR 0.5's `e/D` is above Han's band; AR 3.5's W/H 0.286 is
below Han's 1-4. AR 3.5 is the only independent-lab check of Han's friction
in a narrow channel, and it reads ~38% low. Accuracy, cross-source.

Still unscored: the one-side-turbulated series (refused by the runner: Han's
decomposition is for two ribbed walls). The Nusselt series are scored by
this source's own ratio sets (#444).

## The original UNSCORED note (superseded for friction)

This source publishes a passage-average **Fanning** friction factor
(`f = dP D_H g_c / 2 L rho V^2`) and a turbulated-surface Nusselt number.
combaero's rib sets produce `R` and `G`, the law-of-the-wall roughness
functions. Getting from one to the other needs this family's own
decomposition -- Han, Ou, Park & Lei (1989) Eqs. (5)-(9) -- whose `R`
carries an aspect-ratio term INSIDE the logarithm:

    R = (2/f)^1/2 + 2.5 ln{(2e/D)[2H/(W+H)]} + 2.5

That is **not the same quantity** as the square-channel `R` the shipped
sets produce. Scoring it without declaring the conversion is the category
error #389 named and #433 built the machinery to refuse. The runner work
belongs with #402 and #403; the data is committed now because it is
verified, and because the two-side friction set is what those issues have
been waiting for.

## Not implemented

**One-side-turbulated.** Eight series are committed and no combaero set
models the configuration. Kept because the paper's conclusion 4 -- that
one-side Nu approximately equals large-AR two-side Nu -- is a structural
claim worth having data for if a one-sided element ever exists.
