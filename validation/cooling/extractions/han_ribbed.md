# Extraction: Han rib correlations

**Status: items 1-35 CONFIRMED (2026-09-10). Items 36-46 ADDED later and
UNCONFIRMED.**

The originally reviewed items stand. A later pass added the narrow-aspect-ratio
correlation (Eq. 4.19, Figure 4.48), which extends coverage to `W/H` below 1
and contradicts nothing already confirmed. Those items await review.

Every item was checked against the source or resolved by derivation, and the
reviewer has signed off in the review log. This document is released for
implementation under #334.

Modelling decisions accepted with it -- **D1** (treat `R` as constant
in `e+`), **D2** (use `12.31`), **D5** (extend `RibCorrelationSet` rather
than give Eq. 4.17/4.18 its own type), and **D6** (fit `R` for `W/H < 1` to
Fig. 4.48a's drawn line, since no printed equation exists). **D3 was
withdrawn**: it recorded a departure from the printed four-sided friction
factor that turned out not to exist. Nothing in the implementation
diverges from the source, except D6's R, which is explicitly NOT an
extraction and is labelled `Fitted`.

Eq. 4.17 (R vs alpha) and Eq. 4.18 (G, square/rectangular switch)
IMPLEMENTED 2026-09-20 as `combaero.han_park_1988_angled()`.

Eq. 4.19 (G, narrow-channel alpha switch; R from D6's fit) IMPLEMENTED
2026-10-01 as `combaero.han_1989_narrow_channel()` (#402).

Changes after this point go in the review log, not silently into the tables. If
implementation surfaces something the extraction got wrong, that reopens the
status rather than being patched in the code.

Tracked by #334, under #339.

---

## For the reviewer

Six items need judgement against the book. Everything else is either confirmed
by two independent channels or marked missing. Work through these; the detail
for each is in the tables below under the same number.

| # | what to check | why it matters |
|---|---|---|
| **7** | Fig. 4.46 lower panel: is the drawn curve meant to be flat? Its label reads `= 3.2`, a constant, but measurement gives a ~4.5% rise across the range, ~100x the scan skew. Is there surrounding text stating an `e+` range over which `R` is taken as constant? | decides whether an implementation is closed form or needs an inner solve with its own derivative -- the difference between a cheap and an expensive #334 |
| **6** | The four-sided friction factor. Is it printed as `f = f_bar + (H/W)(f_bar - f_s) * f_bar`? The transcription is settled; the question is whether the book itself is right, since as printed the correction is negligible. | 1.18x to 3.77x on `f`, depending on `W/H` |
| **10** | Which of `G = 3.7(e+)^0.28` and `G_bar = 4.5(e+)^0.28` is the ribbed wall and which the channel average? | a 22% error in `St` if swapped |
| **23** | Fig. 4.46 and Eq. 4.18 give different `G`. Are they meant for different configurations -- Fig. 4.46 for 90-degree orthogonal ribs, 4.17/4.18 for broad-aspect ducts with angled ribs? What applies to a 90-degree rib in a broad-aspect duct, which fits both? | up to 22% in `G`, and it decides which correlation an implementation selects |
| **18** | Fig. 4.47, angled ribs | the angled-rib branch cannot be implemented without it, and it is where 4.17/4.18 belong |
| -- | Does the book carry a **worked example** for ribs? | would resolve item 6 and give the whole chain an end-to-end check |

### Flag list -- checkable against the printed page

Every quantity appearing in both the page images and the text extraction
(`sources/heat_transfer/han/han_ribs.md`) was compared. These agree exactly and
need no checking: `u+`, `T+`, `e+`, Eq. 4.15, Eq. 4.16, Eq. 4.18, the `m`/`n`
rules, the `W/H` cap, the validity ranges, `R = 3.2 (P/e/10)^0.35`, and
`G = 3.7 (e+)^0.28`.

**Open: none.** Every flagged item has been checked against the book or
resolved by derivation. The document is ready for a status decision.

**Resolved:**

| # | outcome |
|---|---|
| 6 | Han prints `f = f_bar + (H/W)(f_bar - f_s)`, which the area-weighted assumption derives exactly. The apparent trailing `* f_bar` was a **typesetting artefact in our reading** -- a full stop followed by the next sentence's subject. Three channels misread it; the algebraic cross-check caught it. See D3, withdrawn |
| 7 | `R` is independent of `e+`; valid `e+ >= 50`, 6% for 95% of data. The measured 4.5% curve slope was artefact |
| 10 | `G_bar = 1.2 G`, **stated on Fig. 4.47**. Not a Prandtl normalisation -- an earlier resolution here was wrong and is recorded as such |
| 15-17 | `m` rule, `W/H` cap, Eq. 4.17 validity (p. 377) |
| 18 | Figure 4.47 read; items 28-34 |
| 19-22 | Eq. 4.18 and its `m`/`n` exponents |
| 23 | Fig. 4.46 is Han (1988) JHT 110, 321; Fig. 4.47 is Han and Park (1988) IJHMT 31(1), 183 -- different papers, different configurations |
| 24 | Eq. 4.17 **does** carry `(W/H)^m`. The image was right; the text extraction dropped it |
| 25 | `G ~ Pr^0.57` from Webb, Eckert and Goldstein Eq. (2). A **cross-source decision**, not something already in Han's figure |
| 26 | Eq. 4.14: both wall laws are `2.5 ln(y/e) + term`, with `R(e+)` for `u+` and `G(e+, Pr)` for `T+` |
| 27 | Citation corrected: Webb, R.L. and Eckert, E.R.G. (1972), "Application of rough surfaces to heat exchanger design", *IJHMT* **15(9)**, 1647-1658. The repo says 15(8) |
| 27b | Webb and Eckert (1972) derived the `1/3` exponent. The **name** "thermal performance factor", and the factor as commonly defined, comes later -- Gee, D.L. and Webb, R.L. (1980), "Forced convection heat transfer in helically rib-roughened tubes", *IJHMT* **23(8)**, 1127-1136 |
| 35 | `12.3` (p. 376) vs `12.31` (Fig. 4.47). Resolved by decision D2: use `12.31` |

### Item 25: the Prandtl gap, and the trap in closing it

**Confirmed by the text: the correlations are presented for `Pr ~ 0.7`.**

`G` is written `G(e+, Pr)` in Eq. 4.14 and Eq. 4.16 -- a function of Prandtl
number by definition. Both correlations that let you *evaluate* it,
`G = 3.7 (e+)^0.28` and Eq. 4.18, carry no `Pr` term, and the text ties them to
`Pr ~ 0.7`.

**Preferred candidate: Webb, Eckert and Goldstein (1972).**

> Webb, R.L., Eckert, E.R.G. and Goldstein, R.J. (1972). Generalized heat
> transfer and friction correlations for tubes with repeated-rib roughness.
> *International Journal of Heat and Mass Transfer*, 15(1), 180-184.
> doi:10.1016/0017-9310(72)90179-2

Two structural reasons it fits better than Dipprey and Sabersky:

- **Same roughness family.** It treats *repeated-rib* roughness, which is Han's
  geometry class. Dipprey and Sabersky fitted sand grains, so borrowing from
  them would be an extrapolation across roughness types.
- **Same formalism.** "Generalized ... for tubes with repeated-rib roughness"
  is the `R(e+)` / `G(e+, Pr)` wall-function framework Han builds on, which
  raises the possibility that Han's `Pr ~ 0.7` fit is a special case of it
  rather than something needing a correction bolted on.

**But the book does not cite it.** Section 4.2.3 references only Nikuradse
(1950) and Dipprey and Sabersky (1963). So adopting it is a **modelling
decision made on merit**, not an extraction, and it belongs in the decisions
section with its own justification once the paper is in hand. It must not be
recorded in a way that reads as something Han said.

**It still needs the paper.** A Prandtl exponent recalled rather than read is
the failure this process exists to prevent.

**Where Dipprey and Sabersky (1963) enters.** The book cites it at Eq. 4.14, as
the conceptual origin of the wall-law formulation, alongside Nikuradse (1950) --
not as a correlation Han evaluates. But its title is *"Heat and momentum
transfer in smooth and rough tubes at various Prandtl number"*, which makes it
the natural place to look for a `Pr`-explicit form of `G`.

> Dipprey, D.F. and Sabersky, R.H. (1963). Heat and momentum transfer in smooth
> and rough tubes at various Prandtl number. *International Journal of Heat and
> Mass Transfer*, 6, 329-353.

**Two cautions before that is pursued.**

*It needs the paper.* A `Pr` exponent recalled rather than read is exactly the
failure this process exists to prevent -- see the note on `C = 0.38` becoming
`~40` under "Bearing on the removed implementation".

*Having the paper is not sufficient.* Dipprey and Sabersky fitted **sand-grain
roughness**; Han fitted **ribs**. Carrying a Prandtl exponent from one to
extend the other is not extraction -- it is a cross-source construction, and
cross-source construction is what produced the correlations removed in #332. If
it is done, it belongs in the modelling decisions section with its own
justification and its own error estimate, not in the extracted items.

**A third possible outcome, if Webb, Eckert and Goldstein carries a
`Pr`-explicit `G` for repeated ribs:** adopt it as a labelled decision, and
check that it reproduces Han's `G = 3.7 (e+)^0.28` at `Pr = 0.703` within
Han's stated 8%. That check is the thing that would make the adoption
defensible rather than merely plausible -- and it can fail.

**The two outcomes if it does not:**

1. Han's `G` gets `Pr ~ 0.7` as a **declared validity condition**, with a
   warning outside it. For a `ChannelElement` running combustion products at
   high temperature that will fire regularly, which is a real limit on the
   re-add's usefulness -- but a stated one.
2. A `Pr` generalisation is adopted as a **labelled modelling decision**,
   traceable to its own source, with the extrapolation from sand grains to ribs
   named as the assumption it is.

What is not available is quietly evaluating a `Pr`-dependent function with no
`Pr` in it.

### Superseded framing

`G` is written `G(e+, Pr)` in Eq. 4.14 and Eq. 4.16 -- it is a function of
Prandtl number by definition. But both correlations that would let you
*evaluate* it, `G = 3.7 (e+)^0.28` and Eq. 4.18, carry no `Pr` term at all, and
the text attaches the first to `Pr = 0.703` specifically.

A `ChannelElement` computes `Pr` from the local mixture and temperature, so it
will routinely sit away from 0.703 -- combustion products and high temperatures
both move it. Using either correlation there is extrapolation in a variable the
correlation does not expose, with no stated behaviour.

Either the book gives a Pr-generalised form elsewhere, or the implementation
must declare `Pr = 0.703` as a validity condition and warn outside it. That is
a decision for the re-add, and it needs settling before implementation rather
than after.

### What I still need sent

- **Figure 4.47** and its surrounding text
- any **worked example** using these correlations
- the text around Figure 4.46 that distinguishes `G` from `G_bar` (item 10)

### How to record a correction

Put it in the review log at the bottom rather than editing the tables. A
correction is evidence about the extraction process, not only about the number
-- item 7 has now been claimed three ways, and that history is worth more than
any one of the three claims.

---

## Source

> Han, J., Dutta, S. and Ekkad, S. (2012). *Gas Turbine Heat Transfer and
> Cooling Technology*. 2nd Edition. CRC Press.

Secondary source. The rib correlations are **after Han, J.C., ASME J. Heat
Transfer, 110, 321, 1988**, per the caption of Figure 4.46. Figure 4.46 also
carries data attributed to Han (1984).

Page images are kept in `sources/scans/`, which is gitignored -- they are
copyrighted and must not be committed. This document is the tracked artifact.

| item | file in `sources/scans/` |
|---|---|
| Eq. 4.15, 4.16, e+, four-sided f | `han_2012_rib_eq_4.15_4.16.png` |
| Fig. 4.46, Eq. 4.17 | `han_2012_rib_fig_4.46_eq_4.17.png` |
| page 377 text: `m` rules, validity, Eq. 4.18 exponents | supplied as text |
| Eq. 4.18 | `han_2012_rib_eq_4.18.png` |

---

## Extracted items

Each line below is what a reviewer should compare against the book, one at a
time. State is one of **confirmed**, **suspect**, or **missing**.

### Definitions

Eq. 4.14 gives the wall laws these all sit on, after Nikuradse (1950) and
Dipprey and Sabersky (1963):

```
u+ = 2.5 ln(y/e) + R(e+)
T+ = 2.5 ln(y/e) + G(e+, Pr)
```

`R` and `G` are the dimensionless velocity and temperature **at the rib tip**,
`y = e`. This is why `G` has `R`'s structure with temperature in place of
velocity, and it is the external anchor that made cross-check 1 below
non-circular.


| # | item | as extracted | state |
|---|---|---|---|
| 1 | dimensionless average velocity | `u+ = (2/f)^(1/2)` | confirmed |
| 2 | dimensionless average temperature | `T+ = (f/2)^(1/2) / St` | confirmed |
| 3 | roughness Reynolds number | `e+ = (e/D) * Re * (f/2)^(1/2)` | confirmed |
| 4 | Eq. 4.15, roughness function | `R(e+) = (2/f)^(1/2) + 2.5 ln( (2e/D) * (2W/(W+H)) ) + 2.5` | confirmed |
| 5 | Eq. 4.16, heat-transfer roughness function | `G(e+,Pr) = R(e+) + [ (f/(2*St_r)) - 1 ] / (f/2)^(1/2)` | confirmed |
| 6 | four-sided ribbed channel friction factor | `f = f_bar + (H/W) * (f_bar - f_s) * f_bar` | **suspect** |

**Item 6 is the one to check most carefully.** Two independent readings -- a
visual read and macOS Vision OCR -- both report a trailing multiplication by
`f_bar` after `(f_bar - f_s)`, and a high-resolution crop shows a centred dot
followed by an overbarred `f`. So the *transcription* is settled. What is not
settled is whether it is correct **in the book**: as printed the sidewall
correction is negligible, which is odd for something given its own equation.

| H/W | as printed | without the trailing `* f_bar` | ratio |
|---|---|---|---|
| 0.25 | 0.020075 | 0.023750 | 1.18x |
| 1.0 | 0.020300 | 0.035000 | 1.72x |
| 4.0 | 0.021200 | 0.080000 | 3.77x |

(evaluated at `f_bar = 0.020`, `f_s = 0.005`)

Resolving this needs a worked example in the book, or Han (1988) itself. It
must **not** be resolved by judging which looks more physical -- that is the
move that produced the defects this rebuild exists to undo.

### Correlations, 90 degree orthogonal ribs (Figure 4.46)

| # | item | as extracted | state |
|---|---|---|---|
| 7 | roughness function | `R = 3.2 * (P/e/10)^0.35`, **independent of `e+`**, for `e+ >= 50` | confirmed by the section text |
| 8 | heat-transfer roughness function, solid line | `G = 3.7 * (e+)^0.28` | confirmed |
| 9 | heat-transfer roughness function, dashed line | `G_bar = 4.5 * (e+)^0.28` | confirmed |
| 10 | what `G_bar = 4.5 (e+)^0.28` is | the **Prandtl-normalised** roughness function, `G_bar = G * Pr^(-0.57)`, in Webb, Eckert and Goldstein's notation | **resolved -- see below** |
| 11 | figure validity, this study | `alpha = 90 deg`, `10,000 <= Re <= 60,000` | confirmed |
| 12 | figure validity, Han (1984) data | `alpha = 90 deg`, `8,000 <= Re <= 80,000` | confirmed |
| 13 | plotted e+ range | roughly 40 to 1000 | confirmed |

**Item 7 is now settled by the running text**, which states it directly:

> The correlation of friction roughness function R is `R = 3.2*((P/e)/10)^0.35`
> for `e+` greater than or equal to 50. The equation correlates 95% of the
> experimental data within 6% deviation. **Note that R is independent of the
> roughness Reynolds number e+.** This implies that the average friction factor
> is independent of roughness Reynolds number.

So `R` does not depend on `e+`, the prediction chain below is closed form, and
the ~4.5% slope measured off the drawn curve is drafting or marker
contamination rather than physics. Two further facts come with it: the
correlation holds for **`e+ >= 50`**, and it carries a stated accuracy of
**6% for 95% of the data**.

The measurement that prompted the doubt is kept below, because a reading that
was challenged, measured, and then confirmed by the text is a better record
than a reading that was simply right.

**Superseded discussion.** The label is printed as
a constant, `R/(P/e/10)^0.35 = 3.2`, joined to the plotted curve by a leader
line. But the lower panel's y-axis is **logarithmic**, and the drawn curve is
not quite flat.

Measured from the page image, using the panel frame as a control for scan skew:

| structure | fitted slope | rise across 1000 px |
|---|---|---|
| upper panel bottom frame (control) | -0.000124 px/px | 0.12 px |
| the R curve | -0.011800 px/px | 11.80 px |

The curve's slope is ~100x the frame's, so it is **not** scan skew. Across the
plotted range that is roughly 10.6 px, or about **4.5%** on a log axis of
~554 px/decade -- equivalently `R` proportional to `(e+)^0.02`.

So the two readings disagree:

- the **printed equation** says `R` is a constant, which makes Eq. 4.15 invert
  directly for `f` and the whole chain closed form
- the **drawn curve** rises ~4.5% across the range, which would make `R` weakly
  dependent on `e+`, and the chain implicit again

4.5% is small -- propagated through Eq. 4.15 it moves `f` by about 4% -- but
the distinction decides whether an implementation needs an inner iteration.
A reviewer with the book should judge whether the drawn curve is meant to be
horizontal (drafting, or my measurement picking up the marker trend) or
genuinely sloped.

**This corrects an earlier reading of mine.** I recorded the panel as a
horizontal line from visual inspection, and asserted on #334 that the chain is
therefore closed form. The measurement above was prompted by a reviewer
challenging that reading, and it does not support the confident version.

Geometries carried in the figure legend, usable as harness cases:

| set | e/D | P/e | W/H |
|---|---|---|---|
| this study | 0.047 | 10 | 1 |
| this study | 0.047 | 20 | 1 |
| this study | 0.047 | 10 | 2 |
| this study | 0.047 | 20 | 2 |
| this study | 0.078 | 10 | 4 |
| this study | 0.078 | 20 | 4 |
| Han (1984) | 0.063 | 10 | 1 |
| Han (1984) | 0.063 | 20 | 1 |
| Han (1984) | 0.042 | 10 | 1 |
| Han (1984) | 0.021 | 10 | 1 |

### Item 10: `G_bar = 1.2 G`, stated on Fig. 4.47

**This corrects an earlier resolution in this document.** Figure 4.47 labels its
dashed line directly:

```
G_bar = 1.2 G
```

That is the answer. `G_bar` is not a Prandtl normalisation.

Physically it is consistent: `G` is inversely proportional to Stanton number,
so `G_bar = 1.2 G` describes the *lower*-heat-transfer quantity -- a four-wall
channel average dragged down by the smooth walls, against `G` for the ribbed
wall alone.

**The withdrawn argument, kept because the way it failed is instructive.**
Webb, Eckert and Goldstein's Eq. (2) defines `g_bar = G Pr^(-0.57)`, and at
`Pr = 0.703` that factor is 1.2225. Applied to Han's `G = 3.7 (e+)^0.28` it
predicts 4.5231 against the printed 4.5 -- 0.51%. That looked like a
confirmation across two independently read sources.

It was a coincidence. The two candidate explanations of Fig. 4.46's ratio:

| explanation | factor | implies `G_bar` | vs printed 4.5 |
|---|---|---|---|
| `G_bar = 1.2 G`, **stated on Fig. 4.47** | 1.2000 | 4.440 | 1.33% |
| `G_bar = G Pr^(-0.57)`, inferred | 1.2225 | 4.523 | 0.51% |

The inference fits the digits better and is still wrong, because the other is
printed. **A stated relationship beats a better-fitting coincidence** -- and a
numerical agreement at half a percent is not evidence of mechanism when the
mechanism is written down two pages away.

Note also that this restores the reading item 10 started with -- ribbed wall
versus channel average -- which was abandoned on the strength of the
coincidence and then confirmed by the figure.

### Correlation with rib angle (Eq. 4.17)

| # | item | as extracted | state |
|---|---|---|---|
| 14 | Eq. 4.17 | `R / [ (P/e/10)^0.35 * (W/H)^m ] = 12.3 - 27.07*(alpha/90) + 17.86*(alpha/90)^2` | confirmed |
| 15 | the exponent `m` in Eq. 4.17 | `m = 0` for `alpha = 90 deg`; `m = 0.35` for `alpha < 90 deg` | confirmed (p. 377) |
| 16 | aspect-ratio cap | if `W/H > 2`, set `W/H = 2` | confirmed (p. 377) |
| 17 | Eq. 4.17 validity | `P/e = 10-20`, `e/D = 0.047-0.078`, `alpha = 90-30 deg`, `W/H = 1-4`, `Re = 10,000-60,000` | confirmed (p. 377) |
| 18 | Figure 4.47, angled ribs | referenced in the text, not extracted | **missing** |

**Item 15 carries its own corroboration.** `m = 0` at `alpha = 90 deg` means
`(W/H)^m = 1` for *any* aspect ratio. Figure 4.46 is all-90-degree data and its
legend spans `W/H = 1, 2 and 4` -- all of which fall on the single line at 3.2.
The rule and the figure support each other without either being derived from
the other.

### Correlation for G with e+ (Eq. 4.18, page 377)

| # | item | as extracted | state |
|---|---|---|---|
| 19 | Eq. 4.18 | `G = 2.24 * (W/H)^0.1 * (alpha/90)^m * (P/e/10)^n * (e+)^0.35` | confirmed |
| 20 | exponents, square channel | `m = 0.35`, `n = 0.1` | confirmed (p. 377) |
| 21 | exponents, rectangular channel | `m = n = 0` | confirmed (p. 377) |
| 22 | stated consequence | rib angle `alpha` and rib spacing `P/e` are not significant for `G` in a rectangular channel | confirmed (p. 377) |

**Symbol collision, worth a comment in any implementation.** `m` in Eq. 4.17
(0 or 0.35, selected by rib angle) and `m` in Eq. 4.18 (0.35 or 0, selected by
channel shape) are different quantities that happen to share a letter and even
share the value 0.35. Items 15 and 20 must not be conflated.

Item 19 shows what items 20 and 21 attach to: `m` is the exponent of
`(alpha/90)` and `n` the exponent of `(P/e/10)`. Both vanish for a rectangular
channel, which is the stated consequence in item 22. Note that `(W/H)^0.1` sits
**outside** the switch and therefore applies always.

### Conflict: two G correlations that disagree (item 23)

| # | item | state |
|---|---|---|
| 23 | Fig. 4.46's `G` and Eq. 4.18's `G` disagree | **resolved -- different configurations** |

At `alpha = 90 deg`, `P/e = 10`, over the plotted range:

| `e+` | Fig. 4.46 | Eq. 4.18, `W/H`=1 | ratio |
|---|---|---|---|
| 40 | 10.39 | 8.15 | 0.78 |
| 100 | 13.43 | 11.23 | 0.84 |
| 400 | 19.81 | 18.24 | 0.92 |
| 1000 | 25.60 | 25.13 | 0.98 |

The `e+` exponents differ (0.28 against 0.35), and Eq. 4.18's `(W/H)^0.1` term
lies outside the `m`/`n` switch, so it cannot collapse `W/H` onto a single line
-- yet Fig. 4.46 shows `W/H` = 1, 2 and 4 doing exactly that.

**Resolved by the section text.** They are correlations for **different
configurations**, and the text frames them separately: Figure 4.46 is Han
(1988) for two-sided orthogonal 90-degree ribs, while Figure 4.47 covers
broad-aspect ratio rectangular ducts with angled ribs. They are not obliged to
agree, and the earlier concern about the R/G asymmetry falls away -- the R
comparison being close at 90 degrees is a coincidence of the fits, not evidence
that the two describe the same thing.

Original reasoning, kept for the record: Fig. 4.46 is Han (1988) for 90-degree orthogonal ribs, while
Eqs. 4.17 and 4.18 belong to Figure 4.47, introduced in the text as
"broad-aspect ratio rectangular ducts with angled ribs". On that reading they
are not obliged to agree.

What makes it worth a reviewer's judgement is the **asymmetry**: at
`alpha = 90 deg` the two `R` correlations agree to 3.44%, while the two `G`
correlations differ by up to 22%. Two correlations for genuinely different
configurations would not be expected to agree closely on one function and
poorly on the other.

An implementation must know which correlation applies where, and what happens
at the boundary -- a 90-degree rib in a broad-aspect duct satisfies both
descriptions.

---

### Figure 4.47 (item 18) and what it settles

Read from `sources/heat_transfer/han/han_fig4.47_pp377.png`.

| # | item | as extracted | state |
|---|---|---|---|
| 28 | Fig. 4.47 left panel, R correlation | `R/[(P/e/10)^0.35 (W/H)^m] = 12.31 - 27.07(alpha/90) + 17.86(alpha/90)^2` | confirmed |
| 29 | Fig. 4.47 right panel, G correlation | `G = 2.24 (W/H)^0.1 (alpha/90)^m (p/e/10)^n (e+)^0.35` | confirmed, matches Eq. 4.18 |
| 30 | dashed line | `G_bar = 1.2 G` | confirmed -- resolves item 10 |
| 31 | exponent switch | square channel `m = 0.35, n = 0.1`; **rectangular channels I and II** `m = 0, n = 0` | confirmed, matches items 20-21 |
| 32 | validity box, both panels | `P/e = 10-20`, `e/D = 0.047-0.078`, `alpha = 90-30 deg`, `W/H = 1-4`, `Re = 10,000-60,000` | confirmed, matches item 17 |
| 33 | `e+` axis range | 50 to 1000 | consistent with the `e+ >= 50` floor |
| 34 | source of Fig. 4.47 | Han, J.C. and Park, J.S., *Int. J. Heat Mass Transfer*, 31(1), 183, 1988 | confirmed |

**Item 34 independently confirms item 23.** Figure 4.46 is credited to Han,
J.C., *ASME J. Heat Transfer*, 110, 321, 1988; Figure 4.47 to Han and Park,
*IJHMT*, 31(1), 183, 1988. Two different papers, two different configurations.
They are not obliged to agree, which is what item 23 concluded on weaker
grounds.

### Discrepancy inside the book (item 35)

| # | item | state |
|---|---|---|
| 35 | Eq. 4.17's leading constant | **12.3** in the printed equation on p. 376, **12.31** in Figure 4.47 | needs review |

Small -- at `alpha = 90 deg` it moves `R/(P/e/10)^0.35` from 3.09 to 3.10, about
0.3% -- but it is a disagreement between two places in the same book rather
than a reading problem, and an implementation has to pick one. `12.31` is the
likelier intent, being the more precise of the two, but that is an inference.

### Transcription discrepancy between two of our own records (item 24)

| # | item | state |
|---|---|---|
| 24 | Eq. 4.17's denominator | **needs a third look** |

Two transcriptions of Eq. 4.17 disagree:

| source | denominator |
|---|---|
| page image `han_2012_rib_fig_4.46_eq_4.17.png` | `(P/e/10)^0.35 * (W/H)^m` |
| text extraction, `sources/heat_transfer/han/han_ribs.md` line 125 | `(P/e/10)^0.35` only |

The image version is the more likely correct one, because the sentence
immediately following the equation reads "where `m = 0` for `alpha = 90` and
`m = 0.35` for `alpha < 90`" -- which is meaningless unless `m` appears in the
equation. Recorded rather than assumed: it needs one look at the printed page.

This is what two extraction channels are for. The text extraction carried
things the images did not -- Eq. 4.14, the `e+ >= 50` floor, the stated
deviation bands -- and the image carried a factor the text dropped.

## Narrow-aspect-ratio channels (Eq. 4.19, Figure 4.48) -- ADDED 2026-09-13

Extends coverage to `W/H < 1`. The confirmed items cover `W/H = 1-4`; this
covers `1/4 < W/H < 1`, so together they span `1/4` to `4`.

after: Han, J.C. et al., *Int. J. Heat Mass Transfer*, **31**(1), 183, 1988,
with the `G` correlation attributed in the text to Han et al. (1989).

| # | item | as extracted | state |
|---|---|---|---|
| 36 | Eq. 4.19 form | `G = C (e+)^n`, at `Pr = 0.7` | confirmed by text |
| 37 | `1/2 < W/H < 1` | `n = 0.35`; `C = 2.24` at `alpha = 90 deg`, `C = 1.80` for `30 < alpha < 90` | confirmed by text |
| 38 | `1/4 < W/H < 1/2` | `n = 0.35 (W/H)^0.44`; `C = 2.24 (W/H)^-0.76` at `alpha = 90 deg`, `C = 1.80 (W/H)^-0.76` for `30 < alpha < 90` | confirmed by text |
| 39 | `R` behaviour, Fig. 4.48 | `R` increases with decreasing `W/H`; all three aspect ratios coincide at `alpha = 30 deg` | confirmed by text |
| 40 | `R` correlation for narrow channels | plotted in Fig. 4.48a, no closed form in the extracted text | **resolved by D6: fitted, not extracted** |

### Cross-check: Eq. 4.18 and Eq. 4.19 meet exactly at `W/H = 1`

Different papers, adjoining ranges. At their shared boundary, `alpha = 90 deg`,
`P/e = 10`:

| `e+` | Eq. 4.18 | Eq. 4.19 | difference |
|---|---|---|---|
| 50 | 8.8082 | 8.8082 | 0 |
| 100 | 11.2266 | 11.2266 | 0 |
| 300 | 16.4908 | 16.4908 | 0 |
| 1000 | 25.1332 | 25.1332 | 0 |

Exact to every digit: Eq. 4.18's square-channel branch reduces to
`2.24 (e+)^0.35`, which is Eq. 4.19's upper branch. A check that could have
failed and did not.

### Discontinuity at the `W/H = 1/2` branch boundary (item 41)

| # | item | state |
|---|---|---|
| 41 | Eq. 4.19's two branches disagree at `W/H = 1/2` | **RESOLVED, see below** |

Evaluated at `W/H = 1/2`, `alpha = 90 deg`, from either side:

| `e+` | upper branch | lower branch | ratio |
|---|---|---|---|
| 100 | 11.2266 | 12.4459 | 1.109 |
| 1000 | 25.1332 | 22.5435 | 0.897 |

A **10% step** at the boundary, in opposite directions at the two ends of the
range. Unlike the `W/H = 1` junction above, these branches do not meet.

The text gives open intervals, so `W/H = 1/2` exactly may be intended to fall
in neither. But an implementation has to choose, and a 10% discontinuity inside
the validated range will appear in any solve that crosses it.

**Resolved by digitising Figure 4.48's own `W/H = 2/4` correlation line**
(`validation/cooling/data/han2012/fig4.48_G_correlation_WH0.5.csv`). Free fit
gives `G = 2.458 (e+)^0.3302`, against:

| branch | form at `W/H=0.5`, `alpha=90` | match to the digitised line |
|---|---|---|
| wide (`1/2<W/H<1`) | `C=2.24`, `n=0.35` | 10% in `C`, 6% in `n` |
| narrow (`1/4<W/H<1/2`) | `C=3.79`, `n=0.258` | 69% in `C`, 22% in `n` |

The wide-branch form is unambiguously closer, on both parameters. Han's own
plotted correlation at the boundary uses the wide branch. The discontinuity
itself is not removed -- the two forms still disagree by 10% approaching `1/2`
from below -- but which side the boundary belongs to is now settled by the
source's own figure rather than left to an implementer's guess.

### Eq. 4.17 outside its stated 30-90 degree range -- not a validity claim

Han states Eq. 4.17 for `alpha = 30-90 deg`. Investigated the raw quadratic's
behaviour beyond that purely to characterise it numerically -- e.g. what a
solver would see if an iterate briefly overshot the range -- not to extend the
correlation's claimed validity.

At `W/H=1`, `m=0`:

| `alpha` | `R` (raw, unclamped) |
|---|---|
| 0 | 12.31 |
| 30 | 5.27 |
| 68.2 (vertex) | 2.05 |
| 90 | 3.10 |
| 120 | 7.97 |
| 150 | 16.80 |
| 180 | 29.61 |

The quadratic has no real roots (discriminant `-146.6`), so it stays positive
and finite for any `alpha` -- no blow-up, no sign change. But its vertex sits
at `alpha=68.2 deg`, not `90`, so the raw formula has no reflective symmetry
about `90` and simply keeps climbing past it: `R(180) = 29.6`, higher than
`R(0) = 12.3`, a number with no connection to anything measured.

**Reflection-folding gives a bounded, physically motivated alternative.**
Mapping `alpha_eff = 180 - alpha` for `alpha > 90` (equivalently `abs`-folding
into `[0, 90]`) exactly reproduces the fitted 0-90 curve by construction, so it
never exceeds the range's own extremes (`R` between `2.05` and `12.31`):

| `alpha` | raw (unclamped) | folded (`alpha_eff = 180 - alpha`) |
|---|---|---|
| 120 | 7.97 | 2.20 (= `R(60)`) |
| 150 | 16.80 | 5.27 (= `R(30)`) |
| 180 | 29.61 | 12.31 (= `R(0)`) |

This is the same convention `RibCorrelationSet::symmetric` already applies
elsewhere (true for 90 deg orthogonal and parallel angled ribs, where
reversing the flow is a mirror operation) -- reflection-folding the *angle*
input is that same idea applied consistently, not a new one.

**Not implemented.** No guard exists in `evaluate_rib` for out-of-range
`alpha`; this is a numerical characterisation for whoever adds one, recording
that reflection-fold is the bounded, convention-consistent choice and that
leaving the raw polynomial unclamped is not dangerous (no negative or complex
output) but is unbounded and physically unmotivated beyond the fitted range.

**FALSIFIED BY DATA (2026-09-24, #385).** Reflection-folding is bounded and
convention-consistent, and it is also **wrong as physics**. Lau, Kukreja and
McMillin (1990) measured both sides of the fold -- V-45 and V-60 against their
reversals V-135 and V-120 -- and they are not the same:

| `e+` | `G(V-60)` | `G(V-120)` | | `G(V-45)` | `G(V-135)` | |
|---|---|---|---|---|---|---|
| 100 | 8.16 | 10.03 | **+22.9%** | 9.33 | 10.55 | **+13.1%** |
| 400 | 14.19 | 16.34 | +15.2% | 15.26 | 17.43 | +14.2% |

Higher `G` is lower heat transfer, so reversed arrays transfer less, and the
paper's own conclusion 4 states it independently: "Reversing the 45 and 60
V-shaped rib arrays lowers the heat transfer from the channel walls and
increases the channel pressure drop."

So a guard **must not** be written on the reflection-fold basis. The bounded
behaviour it offers is real; the physical claim underneath it is not. If a
guard is ever needed past 90 deg, it should refuse or extrapolate explicitly
rather than pretend `alpha` and `180 - alpha` are equivalent. See
`lau_1990_v_ribs.md`.

### Prandtl number

Item 36 states `Pr = 0.7`, consistent with item 25. The narrow-channel
correlation carries no `Pr` term either.

### Figure 4.48 digitised: both panels, plus a labelled `[Ref. 3]` curve

`validation/cooling/data/han2012/fig4.48_*.csv`, 16 series, 111 points.
Attribution by marker fill (solid = `W/H=1/4`, open = `W/H=2/4`), cross-checked
against the physics: `R` must rise as `W/H` falls, and the solid points do
measure higher `R` than the open points at every `alpha`.

A third curve, labelled `W/H = 2/2 [Ref. 3]` on the figure, carries no scatter
of its own -- only a drawn line, with error bars at a few points and no raw
markers. `[Ref. 3]` is not resolved to a specific citation; what is measurable
is that this curve's `R` matches Eq. 4.17 (`W/H=1`, `m=0`) to 4.16% RMS, and its
`G` matches Eq. 4.18's square-channel form (`C=2.24`, `n=0.35`) to 3.70% RMS --
a third confirmation of both equations, after Figures 4.46 and 4.47.

**A cross-paper gap surfaced by the same comparison.** Eq. 4.17 evaluated at
`alpha=90` gives `R=3.10`, against Han (1988)'s own direct `R=3.2` (`-3.1%`).
Eq. 4.18's square-channel form gives `G=2.24(e+)^0.35`, against
`han_1988_orthogonal`'s `G=3.7(e+)^0.28` -- diverging from `-20.9%` at `e+=70`
to `-11.4%` at `e+=404`. Both are small, real disagreements between two
published forms for nominally the same case, not a defect in either. **Nothing
in Figure 4.48 is scored against `han_1988_orthogonal`** for three independent
reasons: the new `W/H<1` data falls outside its valid range; the `[Ref. 3]`
curve is a drawn line (never a source of truth, per the harness rule); and
scoring the overlapping points would report the cross-paper gap as model
error. See the metadata header for the full reasoning.

**Error bars gave a third independent measurement of Han's accuracy band.**
Four reconstructed centers on the `R` panel: half-widths 6.5%, 8.0%, 6.2%,
8.2%, mean 7.2%. This is the first time in this dataset the uncertainty comes
from the source's own printed error bars rather than from scatter or a stated
blanket percentage.

## High-performance ribs (section 4.2.4) -- ADDED 2026-09-15

Angled-rib derivatives: V-shaped, broken, wedge- and delta-shaped. Recorded for
scope, because the section contains **no closed-form correlations** -- it is
comparative, presenting `G` and `R` against `e+` as plots.

| # | item | as extracted | state |
|---|---|---|---|
| 42 | V-shaped ribs, Figs. 4.49-4.51 | after Han, J.C. et al., *ASME J. Heat Transfer*, **113**, 590, 1991 | confirmed |
| 43 | broken ribs, Figs. 4.52-4.54 | after Han, J.C. and Zhang, Y.M., *IJHMT*, **35**(2), 513, 1992 | confirmed |
| 44 | wedge- and delta-shaped ribs, Figs. 4.55-4.58 | after Han, J.C. et al., *Enhanced Heat Transfer*, **1**(1), 37, 1993 | confirmed |
| 45 | closed-form correlations for any of these | **none given** -- performance is reported as plots and as ratios in the text | **missing** |
| 46 | best performers, per the text | 60 deg V-shaped broken rib; backward-aligned delta ribs at 3-4x heat transfer for 7-9x pressure drop | confirmed |

**Implication for #334.** These configurations cannot be implemented as
correlations from this source. What the section supports is a *comparative*
statement -- that a 60 deg V-broken rib outperforms a 90 deg continuous rib --
not a predictive one. Implementing them would need the three primary papers,
and is out of scope for the current issue.

### Cross-check: Figure 4.54 independently confirms items 8 and 10

Figure 4.54 (Han and Zhang 1992) plots `G_bar` against `e+` for broken ribs,
with the 90 deg continuous rib as a reference line. The section text states:

> The heat-transfer roughness function for the 90 deg rib is the same as the
> previous correlation that was developed for the 90 deg rib by Han (1988).

Testing that against the digitised line, using **only** quantities from the
confirmed items -- `G = 3.7 (e+)^0.28` (item 8) and `G_bar = 1.2 G` (item 10):

| comparison | mean pred/data | RMS |
|---|---|---|
| `G_bar = 4.44 (e+)^0.28` vs Fig. 4.54 line | 0.967 | **3.5%** |
| same vs Fig. 4.54 data points | 0.949 | 5.4% |

Free fit to the digitised line gives `G_bar = 4.201 (e+)^0.2952` against the
confirmed `4.44 (e+)^0.28`.

**This corroborates item 10 from a third independent direction**, which matters
because item 10 is the one this document got wrong and then right: it was first
assumed, then "resolved" via a Prandtl coincidence that fitted to 0.51% and was
wrong, then settled by a label printed on Figure 4.47. A different paper, a
different figure and a different rib family now reproduce it to 3.5%.

### Digitised reference point for a high-performance rib

Not a correlation, and not to be implemented as one. Recorded because it
quantifies the section's qualitative claim:

```
90 deg continuous (Fig. 4.54 line)   G_bar = 4.201 (e+)^0.2952
45 deg V-shaped broken (data)        G_bar = 3.242 (e+)^0.2882
at e+ = 500:  26.3 against 19.4, a ratio of 0.739
```

Lower `G` means better heat transfer, so the 45 deg V-broken rib sits about 26%
better on this measure -- consistent with the text, which calls the broken
V-shaped ribs the best performers of the group.

## Modelling decisions

Choices made in implementing what was extracted. Kept apart from the extracted
items above, because those are facts about a page and these are not -- a
reviewer should be able to disagree with one of these without doubting the
transcription.

### D7. The element's drop: Han's area weighting run forwards, in Darcy form (#1106)

**Decision.** A ribbed channel's drop uses
`f_D = 4 [phi f + (1 - phi) 0.046 Re^-0.2]` (C++ `rib_channel_friction`),
`phi` the ribbed perimeter fraction (`W/(W+H)` for two walls, half for one,
1 for four).

**Why.** Two definitions in Han and Park (1988). Fanning: p. 187 defines the
measured `fbar = dP / {4 (dL/D)[G^2/(2 rho g_c)]}` against the Fanning Blasius
`f(FD) = 0.046 Re^-0.2` (Eq. 8). Four-sided: Eq. 6 forms `f` from `fbar` by area
weighting. The element previously used `f (L/D) rho V^2/2` -- Fanning as Darcy
and four-sided for any wall count -- and scored -59.1% on Han's own two-wall
data through the element (`rib_element_friction_runner`), -0.9% / 5.9% MAE
after. The correlation-level harness never saw it: it inverted Eq. 6 itself.

**Alternatives.** Keeping `f` and only multiplying by 4: 1.8x high at W/H 1
for two walls. A per-wall-count constant: no source.

**Outside the measured configuration.** No Han paper measured one or four
ribbed walls, so both set `extrapolated`. Evidence for the weighting there is
cross-source: Taslim and Spring's one-side series agree model-free
(`fbar_1 = (fbar_2 + f_s)/2`) within -12%..+11% at e/D <= 0.107 and miss by
+50% to +275% at e/D >= 0.125 (blockage beyond every Han set); Chandra and
Cook's printed two- and four-wall drop ratios (Han 2013 p. 414) agree within
3%. Pinned by `python/tests/test_rib_element_friction.py` and
`test_ribbed_channel_element.py`.

### D5. Extend `RibCorrelationSet` rather than give Eq. 4.17/4.18 a separate type

**Decision, 2026-09-20.** `han_park_1988_angled()` needed two shapes the
existing `RibCorrelationSet` schema could not express: R as a QUADRATIC in
`alpha` rather than a power law (Eq. 4.17), and G's `alpha`/`p_e` exponents
SWITCHING on channel shape (Eq. 4.18's square-vs-rectangular split) rather
than being fixed. Both are genuine discontinuities the source states, not
numerical artefacts.

Two ways to add this: extend the shared struct with new, opt-in fields
(`RAlphaShape`/`GShapeModel` discriminators plus the quadratic/switch
fields), or give the angled correlation its own parallel type and
`evaluate_rib`-like function.

**Chose extension.** `RibbedModel.correlation_set` (the Python/GUI surface)
already accepts any `RibCorrelationSet` and needs nothing else to use a new
one -- the whole point of "users pick and define correlations, with
provenance and validity visible" (this issue's own framing). A parallel
type would break that substitutability and need dispatch logic wherever a
set is consumed. The cost is a larger shared struct; the two new
discriminators default to the plain power-law behaviour, so
`han_1988_orthogonal` and `rallabandi_2009_high_re` are provably unaffected
-- their full pre-existing test suites pass unmodified, byte for byte,
after this change.

**What is NOT done: no smoothing across either switch.** At `alpha = 90`
the `(W/H)^m` jump reaches 62.5% at `W/H = 4`; at `W/H = 1` the
square/rectangular switch reaches 27% at low `alpha`. Neither source states
a smooth transition, so none is invented -- consistent with the earlier
decision to record, not guess, Eq. 4.17's behaviour outside 30-90 degrees.
A caller whose solver traverses either boundary exactly needs a guard this
evaluator does not supply; the jump sizes are pinned by
`HanParkRAndGSwitchesAreGenuineDiscontinuities` in `test_rib_correlation.cpp`
so a future "helpful" smoothing change cannot land silently.

### D6. Fit `R` for `W/H < 1` to Fig. 4.48a's drawn line -- not an extraction

**Decision, 2026-10-01 (#402).** Item 40: Han et al. (1989)/Fig. 4.48 give a
closed-form `G` (Eq. 4.19) but print no equation for `R` below `W/H = 1` --
only a drawn correlation line. `han_1989_narrow_channel()` needed a value
for `R` anyway, since it feeds directly into `f` through the wall-law
inversion every other set in this file also uses.

**What was fitted, and against what.** A quadratic in `alpha/90`, the same
functional shape Eq. 4.17 already uses for angled ribs -- one triple of
coefficients per `W/H` sub-band (`1/2 <= W/H < 1` and `1/4 < W/H < 1/2`,
the same split item 41 already established), least-squares fit to the two
already-digitised drawn-line series,
`fig4.48_R_correlation_WH0.5.csv`/`_WH0.25.csv`. RMS against the line
itself: 0.99% (wide band), 0.62% (narrow band) -- both far inside Han's
usual 6% scatter, so the quadratic shape is not a stretch for this family.

**Why this is not an extraction.** Every other coefficient in this file is
transcribed from a printed equation and checked against the source. These
two quadratics were never printed anywhere -- they are this project's own
choice of functional form, fitted to a figure we also digitised. That is
exactly what `RibProvenance::Fitted` exists to distinguish from
`Extracted`, and it is the first time this project has actually populated
it (previously defined, unused by any of the three earlier sets).

**The set-level `provenance` is `Fitted`, even though `G` alone is
`Extracted`.** `RibCorrelationSet::provenance` is one enum per set, and
splitting it into `provenance_R`/`provenance_G` would touch the Python
bindings, the GUI surface and every existing test for a distinction that
is explained in one sentence in the factory's own doc comment. The more
cautious of the two values was kept rather than the more flattering one --
a reader who only checks `provenance` should not be told the whole set is
`Extracted` when a third of its claim (the friction side) is not.

**Accuracy is `Measured`, against the RAW scatter, separately from the
fit.** The coefficients above were fitted to the drawn LINE; `accuracy_R`
(pooled RMS 5.4%, bias -0.9%/+0.4% per band) is measured against the 16
raw digitised points in `fig4.48_R_WH0.25.csv`/`_WH0.5.csv` instead -- the
same fit/score separation `rallabandi_2009_high_re` already observes.
Never a judging band (`AccuracyProvenance::Measured`), same as every other
measured figure in this file.

**`G`'s accuracy is also `Measured`, at `alpha = 90` deg uniformly, and
that is a checked simplification, not an unresolved gap.** The raw `G`
scatter (`fig4.48_G_WH0.25.csv`/`_WH0.5.csv`) pools multiple rib angles
onto one e+-vs-G curve with no per-point angle recorded, so Eq. 4.19's
alpha switch (`C=2.24` at 90 deg, `1.80` off-axis -- a nominal 24% jump)
cannot be resolved per point. Scored at `alpha=90` uniformly regardless
(pooled RMS 5.4%, bias +1.9% over 21 points), consistent with how
`G_BAR_OVER_G` in `runner.py` already applies a published constant
everywhere and reports the resulting error as accuracy rather than
refusing. Checked rather than assumed: a free power-law fit to the pooled
scatter ALONE, ignoring the text's alpha switch entirely, already lands
within 3.9-5.0% RMS -- so the angles are not separable on this axis in the
first place, and scoring at one fixed angle is not hiding a 24%-scale
defect.

**What this does NOT do.** It does not extend to any `W/H` or `alpha`
Han's own figure does not cover, and it does not smooth the genuine
discontinuities Eq. 4.19 states (the alpha switch; the `W/H` sub-band
boundary at 1/2). Both are pinned by
`HanNarrowRAndGSwitchesAreGenuineDiscontinuities` in
`test_rib_correlation.cpp`, the same discipline D5 established.

### D2. Use `12.31`, not `12.3`, in Eq. 4.17

**Decision.** Where the printed equation on p. 376 gives `12.3` and Figure 4.47
gives `12.31`, use `12.31`.

**Why.** The figure carries the more precise value, and a rounded `12.3` in
running text is the likelier abbreviation of a fitted `12.31` than the reverse.
The difference is 0.3% at `alpha = 90 deg` -- 3.09 against 3.10 -- so nothing
turns on it numerically. Recorded so the choice is visible rather than
arbitrary.

### D3. WITHDRAWN -- there was no decision to make

**Superseded.** This section previously recorded a decision to implement the
four-sided friction factor as derived, departing from a printed form that
carried a trailing `* f_bar`. There is no such departure. The book prints

```
f = f_bar + (H/W) * (f_bar - f_s)
```

which is exactly what Han's area-weighted assumption derives. The apparent
trailing factor was **a typesetting artefact in our reading**: the dot is a
full stop ending the equation, and the `f_bar` after it is the *subject of the
next sentence* -- "`f_bar` is the average friction factor in a channel with two
opposite ribbed walls" -- sitting on the same line.

**Why this episode is worth keeping.** Three independent extraction channels
all reported the trailing factor: a visual read, macOS Vision OCR, and a
separate text extraction. All three were reading the same real glyphs; all
three misattributed them. Agreement between channels does not establish
meaning, only that they saw the same marks.

What caught it was the **algebraic cross-check**. Deriving `f` from Han's
stated area-weighting gave the form without the trailing factor, disagreeing
with every reading channel. That disagreement was the signal, and it was right.

The tell was also present and misread. The prose immediately after the equation
was flagged in this document as an orphaned subject -- "is the average friction
factor", with nothing in front of it -- and attributed to that region of the
page extracting badly. The cause was the opposite: the subject was not lost, it
had been absorbed into the equation line. The anomaly and the spurious factor
were one artefact, not two.

**Net effect on implementation:** none, beyond simplification. The derived form
and the printed form are the same, so nothing departs from the source and no
worked example in the book can disagree with us on this point.

### D1. Treat `R` as independent of `e+`, following Han

**Decision.** Implement `R = 3.2 (P/e/10)^0.35` as a constant in `e+`, as the
text states, rather than fitting the slope visible in the drawn curve.

**Why it is defensible, not merely convenient.** The drawn curve does rise --
measured at ~4.5% across the plotted range, about 100x the scan skew of the
panel frame, so it is real. But Han's stated accuracy for this correlation is
**6% for 95% of the data**. The neglected slope is *smaller than the scatter
the correlation already admits*, so treating `R` as constant is consistent with
its own accuracy claim rather than in tension with it.

**What it buys.** Eq. 4.15 inverts directly for `f`, making the whole
prediction chain closed form. Keeping the slope would make `R` depend on `e+`,
`e+` depends on `f`, and `f` would then need an inner solve carrying its own
analytic derivative to satisfy the Solver (f, J) rule -- a substantial cost for
a correction below the noise floor.

**What to watch.** If the harness ever shows systematic bias in `f` that grows
with `e+`, this decision is the first place to look. The effect is bounded at
roughly 4.5% and one-signed, so it would appear as a drift rather than scatter.

## Cross-checks performed

These are checks that could have failed and did not. They test the extraction,
not the physics.

**1. Eq. 4.15 and Eq. 4.16 combine into the analogy structure.**
The numerator of 4.16 is `T+ - u+` in disguise, and `u+` is the leading term of
`R`, so `u+` cancels exactly:

```
G = T+ + 2.5 ln( (2e/D) * (2W/(W+H)) ) + 2.5
```

which is `R` with temperature substituted for velocity -- what the
heat-momentum analogy requires. The expectation is external to the page, so
this is not circular. Verified numerically at four operating points, agreeing
to 3.55e-15.

Falsified: five plausible misreadings of 4.16 move `G` by 3.65 to 25.00 at a
reference point where the correct reading gives exactly 0.

**2. Eq. 4.17 reduces to Figure 4.46 at alpha = 90 degrees.**

```
Eq. 4.17 at alpha=90 : 12.3 - 27.07 + 17.86 = 3.0900
Fig. 4.46            :                        3.2
agreement            : 3.44%
```

`(W/H)^m` drops out at `W/H = 1`, so this holds without knowing `m` (item 15).
Two separately extracted equations from different parts of the book agreeing at
their common point.

**3. The prediction chain lands the figure's own legend geometries inside the
figure's own e+ range.** Eleven of twelve fall within 40 to 1000; the twelfth
(e/D=0.078, W/H=4, Re=60,000) reaches 1122, 12% beyond the axis, at the extreme
corner of the plotted set.

**4. Two independent extraction channels.** Visual reading and macOS Vision OCR
(`scripts/ocr_page.swift`). Both report `12.3`, `27.07`, `17.86`, `4.5`, `3.7`,
`0.28` and `3.2` identically. OCR garbles the typeset fractions -- it read
`(3)"/2` for `(2/f)^(1/2)` -- so it is a disagreement detector, not a source of
truth. Its per-line confidence self-flags: the one badly garbled region
reported 0.30 where everything else reported 1.00.

---

## Derived prediction chain

Follows from items 4, 5, 7 and 8, **and only if item 7 is a true constant.**
If `R` carries the ~4.5% `e+` dependence the drawn curve suggests, the first
two steps become implicit and need an inner solve with its own analytic
derivative. The chain below assumes the printed equation.

```
R    = 3.2 * (P/e/10)^0.35                                  item 7
f    = 2 / [ R - 2.5 ln( (2e/D)*(2W/(W+H)) ) - 2.5 ]^2      invert item 4
e+   = (e/D) * Re * (f/2)^(1/2)                             item 3
G    = 3.7 * (e+)^0.28                                      item 8
St_r = f / ( 2 * [ 1 + (G - R) * (f/2)^(1/2) ] )            invert item 5
```

History of this item, kept because the oscillation is itself evidence about the
process. It was first claimed on #334 that the correlation is implicit, from
reasoning about the structure rather than reading the figure. Reading the
figure suggested a horizontal line and a closed-form chain. Measuring the
figure, after a reviewer objected, showed a small but real slope. It is
currently **open**, and the printed equation is the better authority until a
reviewer rules otherwise.

Worked values from the chain, for review against the book if it carries an
example:

| e/D | P/e | W/H | Re | f | e+ | G | St_r |
|---|---|---|---|---|---|---|---|
| 0.047 | 10 | 1 | 10,000 | 0.04576 | 71.1 | 12.21 | 0.00968 |
| 0.047 | 10 | 1 | 60,000 | 0.04576 | 426.6 | 20.16 | 0.00642 |
| 0.078 | 10 | 4 | 10,000 | 0.11503 | 187.1 | 16.01 | 0.01413 |

Note `f` is independent of `Re` in this chain, which follows from item 7 and is
worth a reviewer's attention.

---

## Validation target

Figure 4.46 carries **no data table**, so the scatter would have to be
digitised. But the correlation lines are printed on the plot as equations
(items 7, 8, 9), which splits the job in two:

**The correlation line is the target for correctness, and it needs no
digitising.** An implementation of Han's correlation must reproduce
`R/(P/e/10)^0.35 = 3.2` and `G = 3.7 (e+)^0.28` exactly, because those *are*
the correlation. Reproducing them is a closed-form check.

**The tolerance is stated in the text, so the scatter need not be digitised
after all.** Han gives it directly:

| correlation | stated accuracy | validity |
|---|---|---|
| `R = 3.2 (P/e/10)^0.35` | 95% of data within **6%** | `e+ >= 50` |
| `G = 3.7 (e+)^0.28` | 95% of data within **8%** | `e+ >= 50`, `Pr = 0.703` |

That is a better tolerance than a digitised band would have been: it is the
author's own figure for his own fit, rather than our reading of his plot. The
harness should use 6% and 8%, and treat them as the multiplicative bands they
are on a log axis.

**The scatter, if ever digitised, is the target for tolerance.** How far the experimental points
fall from the printed line is the accuracy of Han's correlation against his own
measurements -- which is what a harness tolerance should be set from, rather
than a number chosen for convenience. Digitising the scatter is therefore worth
doing, but for the band, not for the values.

Both panels use **logarithmic y-axes**, so a band read off them is
multiplicative. A tolerance expressed as a ratio, not as an absolute offset.

Keeping these apart matters. Our code is not obliged to reproduce Han's
measurements; it is obliged to reproduce Han's correlation, and to be honest
about how well that correlation matched reality. Conflating the two is how a
tolerance ends up wide enough to admit a 4-5x error, which is what
`assert 1.0 < multiplier < 10.0` did in the removed implementation.

Not yet done: digitising the scatter of either panel. Until then the harness
has a correctness target but no defensible tolerance.

---

## Bearing on the removed implementation

Not part of the extraction; recorded because it explains what was removed
in #332 and should inform the re-add.

Han gives `R = 3.2 * (P/e/10)^0.35`. The deleted `rib_enhancement_factor` was:

```
3.5 * (e/D)^0.35 * (P/e/10)^-0.15 * f_alpha
```

The constant is close (3.5 against 3.2), the `(P/e/10)` group is exactly Han's,
and the exponent `0.35` is present -- but attached to `(e/D)` rather than to
`(P/e/10)`, and the `-0.15` has no counterpart. It reads as a garbled
recollection of `R`.

More seriously: **`R` is a roughness function, inversely related to the
friction factor.** It is not a Nusselt ratio. The removed code returned it as
`Nu_rib / Nu_smooth`, which is a category error rather than a constant error,
and is consistent with the friction side ending up 4-5x low.

---

## Validity, against what the removed code claimed

Item 17 gives Han's operating range. The implementation removed in #332
validated a **wider** band, so it accepted extrapolation silently and threw
only outside its own inflated limits.

| parameter | Han (item 17) | removed code | effect |
|---|---|---|---|
| `e/D` | 0.047-0.078 | 0.02-0.1 | accepted down to less than half the source's lower bound |
| `P/e` | 10-20 | 5-20 | accepted half the source's lower bound |
| `alpha` | 30-90 deg | 30-90 deg | matches |
| `W/H` | 1-4 | not validated | absent |
| `Re` | 10,000-60,000 | not validated | absent |

A re-add should validate item 17's band, and apply item 16's `W/H` cap rather
than extrapolating past it.

## Open items

**Swept 2026-09-26 against the document's own Resolved table and the
dataset. Almost nothing was still open; the list had simply not been
updated as items closed.** That is the fourth time a superseded record in
this project was read as current (see also item 10 here, the `lau1990`
metadata in #392, and the orphaned docstring in `components.py`).

Closed in the sweep, with where each was already resolved:

- **6** four-sided friction factor -- resolved by derivation, Resolved
  table above. The apparent trailing `* f_bar` was a typesetting artefact
  in our reading; D3 withdrawn. It had stayed on this list saying
  "correctness suspect" after the correctness question was answered.
- **18** Figure 4.47 -- read; items 28-34 confirmed, and six series are
  digitised from it in `data/han2012/`. It had stayed here as "missing".
- **23** -- was recorded TWICE with different content: resolved above as
  "Fig. 4.46 and Fig. 4.47 are different papers and different
  configurations", and open here as "they disagree on `G` by up to 22%".
  Both are true and they are not the same statement. Superseded by the
  finding below, which is evidence rather than an internal comparison.
- ~~digitise the Figure 4.46 scatter to set the harness tolerance~~ --
  **do not.** The digitised subset over-represents the tails about 5x
  (#393): 26% of picked marks fall outside a band the source says holds
  95% of its data, because only spatially isolated marks can be picked. A
  tolerance set from it would be far too wide. CR-3837's tabulated runs
  give the true spread instead.

### Recorded as a known accuracy limit, not an open question

**Eq. 4.18 reads low on `G`.** NASA CR-3837's 33 tabulated runs score
against it with a NEGATIVE bias at every rib angle, -0.8% to -19.9%
(#398). That is out-of-sample: CR-3837 is `e/D` = 0.063, a blockage Han
and Park never measured. It is the same direction as the old item 23
internal disagreement, so an independent dataset now takes the figure's
side.

**A second dataset agrees (2026-10-01/02).** Han, Zhang and Lee (1991)'s
parallel classes, scored at that paper's stated rig (square, `e/D`
0.0625, `P/e` 10), read low against Eq. 4.18: 45 deg -17.2% and 60 deg
-9.3% bias in domain (`G`), 60 deg -8.2% (`G_bar`). Until 2026-10-02 they
were scored at the runner's probe `P/e` 15, which reads Eq. 4.18 ~4% high
through its `(P/e/10)^0.1` and had been hiding part of the gap (see
`han_zhang_lee_1991_jht.md`).

`han_park_1988_angled` ships Eq. 4.18, so this is a property of the model
a user should know, not a defect to fix. Per
`docs/VALIDATION_POLICY.md` it is an ACCURACY result and the harness
reports it as one. **Not reconciled, and deliberately not corrected** --
"improving" a published correlation against one dataset is a new feature
with its own evidence bar, not a bug fix.

### Genuinely still open

**Nothing.** The last item -- whether the book carries a worked example --
was answered on 2026-09-26: **it does not.**

Searched all 865 pages of the text layer: no `Example N` or `Solution:`
block anywhere, no problems or exercises section, no walk-through
language (`consider a`, `as an example`, `let us`, `suppose`) in the rib
chapter, and nothing that substitutes numbers into Eq. 4.17 or 4.18. The
book states correlations and plots data; it never works one through.

**The end-to-end check it would have given now exists anyway, and is
better.** NASA CR-3837's appendix tabulates `Re`, the friction chain, `e+`,
`R` and `Gbar` per run for 33 runs (#398), which exercises the same chain
against MEASURED numbers rather than against the author's own arithmetic.
A worked example would only have confirmed that Han computes his own
equations consistently; the tabulated runs test whether the equations
match a rig.

So item 6 stands as resolved by derivation, and the independent
confirmation it was waiting for came from a different and stronger place.

---

## Review log

| date | reviewer | outcome |
|---|---|---|
| 2026-10-02 | Claude | **Figure 4.51 scored at its stated rig** (`e/D` 0.0625, `P/e` 10, square; Han, Zhang and Lee 1991 p. 590) instead of the runner's probe `P/e` 15, which had read Eq. 4.18 ~4% high. `han_park_1988_angled` in-domain bias: 60 deg parallel -6.8% -> -9.3%, 45 deg -13.8% -> -17.2%. The known low reading, larger; 90 deg unchanged. |
| 2026-10-01 | reviewer + Claude | **#403 item 1 resolved against the primary paper**, Han, Zhang and Lee (1991) JHT 113, 590 (`han_zhang_lee_1991_jht.md`). Its Table 2 prints `G` and `G_bar` for all nine configurations; every confirmed figure 4.51 class lands on its own panel's fit within 2.8%, while BOTH panels of the 45 deg parallel and crossed classes land on the printed `G`. So the `G_bar` panel carried the `G` marks twice: `G_45par`/`G_45crs` confirmed, `Gbar_45par`/`Gbar_45crs` kept disputed and unscored. Printed `G_bar/G` is 1.07-1.26 across all nine, never below 1. Also corrected: 91-GT-3 had been recorded as this paper's conference version; it is the heat-flux-ratio companion, with the same `R` and a `G` up to 23% different. |
| 2026-10-01 | Claude | **Figure 4.53 (Han and Zhang 1992) scored on its own axes (#435).** The 90 deg continuous series goes through a new `f_ratio` path in `runner.py` against `han_1988_orthogonal`: MAE 7.9%, bias -2.1% over 5 points, ACCURACY. The issue said "no source needed"; that was not quite true -- the rig is unstated, so e/D 0.0625 / P/e 10 / square are borrowed from 91-GT-3, and the smooth references (Blasius `0.079 Re^-0.25`, Dittus-Boelter) are the lab's as printed in NASA CR-4015 Eqs. (4)/(5). `f_ratio` is the channel-average `fbar`, not the four-sided `f` (1.8x apart); a test falsifies the four-sided reading. The 60 deg V-shaped series stays refused pending #434. |
| 2026-10-01 | Claude | **Eq. 4.19 implemented** as `combaero.han_1989_narrow_channel()` (#402), settling item 40 via decision D6. `G` is text-extracted (items 36-38); `R` has no printed equation below `W/H=1` and is this project's own quadratic-in-alpha fit to Fig. 4.48a's drawn line (0.99%/0.62% RMS against the two sub-bands' lines), scored separately against the raw digitised scatter (`accuracy_R` 5.4% RMS over 16 points). `G`'s raw scatter pools multiple rib angles with no per-point angle recorded, so it is scored at `alpha=90` uniformly; checked rather than assumed that this does not hide Eq. 4.19's own ~24% alpha switch, since a free fit to the pooled points alone already lands within 3.9-5.0% RMS (`accuracy_G` 5.4% RMS over 21 points). Wired 37 previously-unscored Fig. 4.48 points (`fig4.48_{R,G}_WH{0.25,0.5}`); the `*_WH1_ref3` drawn-line and error-bar series stay deliberately unscored, unaffected -- the issue's "58 points" estimate did not account for that, and the real unlock is 37. New schema shapes (`RAlphaShape::QuadraticAlphaTwoBand`, `GShapeModel::NarrowChannelAlphaSwitch`) follow D5's precedent: additive fields, `han_1988_orthogonal`/`rallabandi_2009_high_re`/`han_park_1988_angled`'s test suites pass unmodified. Writing the Stanton-derivative test caught a real bug before it shipped: the generic `dG/dRe` formula used the shape's unscaled base `e+` exponent rather than the `(W/H)`-corrected one actually applied in the narrow sub-band, silently wrong by as much as 40% there until fixed. |
| 2026-09-26 | reviewer + Claude | **Han's rib implementation wrapped up; open-items list swept.** Items 6, 18 and 23 were all already resolved in this document's own Resolved table while still listed as open -- item 23 twice, with different content in each place. Nothing was reopened by the sweep. Two decisions recorded rather than left hanging. **(1) Eq. 4.18 reads low on `G`**, confirmed out-of-sample by CR-3837 at every rib angle (-0.8% to -19.9%); recorded as a known accuracy limit of the published correlation, not a defect, since "improving" a paper against one dataset is a new feature with its own evidence bar. **(2) No generalised `G_bar/G`.** Pooling 91-GT-3's table 3 (extracted for this), Lau's table 2 and CR-3837's per-run split gives 16 configurations over three rigs spanning 1.096-1.413, and **69% of that variance is BETWEEN rigs** -- no angle or shape term can reach a rig offset, so a "better" correlation would fit rig identity. Han's 1.2 scores MAE 8.3% against the pooled population and is applied as published at every angle (#401). Pinned by a test that fails if the between-rig share drops, which is the condition for revisiting. |
| 2026-09-26 | reviewer + Claude | **`han_park_1988_angled` was scoring every angled series at 90 degrees.** Found while wiring NASA CR-3837 in as a tabulated primary (see `han_park_lei_1984_cr3837.md`): on the `e+` path `run_series` copied `e_D`, `p_e` and `W_H` out of `series.geometry` but never `alpha_deg`, which is a top-level field rather than part of that mapping, so `_probe_geometry`'s default of 90 survived. Every angled series was therefore scored as though its ribs were transverse -- against the one set whose entire subject is rib angle. **This invalidates the "RMS 8.5%" cross-check against figure 4.51's parallel-rib classes recorded in the 2026-09-20 entry below**; that number was computed at 90 degrees for 45 and 60 degree data. Re-scored after the fix, `fig4.51_G_60par` reports MAE 6.8%, bias -6.8%, and the set as a whole moves from bias +4.4% to +1.5%. Fixed with a test that checks the prediction against `evaluate_rib` at the declared angle rather than against a stored number, falsified by reverting the fix. The defect is the shape this project exists to catch: nothing failed, the numbers merely meant something other than what they said. |
| 2026-09-26 | reviewer + Claude | **The `e/D` 0.047 / `P/e` 10 / `W/H` 2 class dispute is closed, from behind the figure rather than by re-reading it.** Its two panels paired on only 2 of 4 marks, and the digitised data could not say which symbol they belonged to. Resolved with **NASA CR-4015** (= AVSCOM 86-C-25), Han, Park and Ibrahim, contract NAS3-24227, September 1986 -- the report Fig. 4.46's "This study" classes come from, traced through the acknowledgement in Han and Park (1988), which names NAS3-24227 where CR-3837 carries NAG3-311. Its appendix 7.3 tabulates every run with a self-identifying header (`E/D`, `P/E`, `ALPHA`, `HYD DIA`), so the class is found BY LABEL and the overplotted cluster is never resolved. Five runs on pp.141-145 at `Re` 10111, 18869, 32117, 60898, 64193 give `e+` 80.7, 150.5, 256.2, 485.8, 512.1, and **every mark in both panels lands on one of them**. Two recorded claims were wrong: the class is not "absent near `e+` ~ 80" (run 22 sits at 80.7 and the `R` panel picked it up to 0.1%), and its predicted ceiling of 479 was too low (the true maximum is 512.1, so the 546.6 mark is 6.7% out, not outside the class). The non-pairing was overplotting: the panels share `e+`, so a mark visible in one proves the run exists in both and silence in the other only means buried -- see #393. Reviewer independently confirmed the one clipped digit (`Re` 32117) by mass-flow ratio against run 24. |
| 2026-09-20 | Claude | **Eq. 4.17/4.18 implemented** as `combaero.han_park_1988_angled()`, following the schema extension recorded as decision D5 (quadratic-in-alpha R, square/rectangular G switch, both opt-in and provably not affecting `han_1988_orthogonal` or `rallabandi_2009_high_re`). Scored against this document's own source, figure 4.47: R at RMS 10.5% over 39 points, G at RMS 8.8% over 115 points (representative geometry, the cloud carries no legend). Independently cross-checked against figure 4.51's parallel-rib classes (a different paper, Han et al. 1991) at RMS 8.5%; the other seven rib shapes on that figure (crossed, V, lambda) are deliberately not scored, since this equation cannot represent rib shape beyond angle and doing so would repeat item 23's lesson. Along the way, corrected two more stale `scores` fields: figure 4.51's 45/60 deg parallel series had been scored against `han_1988_orthogonal`, silently refused every time by the alpha guard (valid range 90-90) since before `han_park_1988_angled` existed to score them correctly. |
| 2026-09-20 | reviewer + Claude | **Item 41 resolved, item 40 given data.** Figure 4.48 digitised: both panels, `[Ref. 3]` curve, and its error bars (16 series, 111 points). The `W/H=2/4` correlation at the branch boundary matches the WIDE branch of Eq. 4.19 (10%/6%), not the narrow one (69%/22%), settling which form Han's own figure uses at `W/H=1/2`. Along the way, corrected an initial misreading of panel (a)'s y-axis as log -- the half-integer minor ticks (1.5, 2.5, 3.5, 4.5) are only possible on a linear axis, confirmed on the printed page and independently by the digitised frame. One error-bar end (alpha=30, top panel) was initially obscured; the rescanned value agreed with the earlier partial read to 0.7%. |
| 2026-09-18 | reviewer + Claude | **Item 7 closed with evidence, not deference.** D1 stands, but the reason improves: the `R` slope is in the DRAWING, not the scan. Both panel frames were digitised as distortion standards -- a frame is horizontal by construction, so its fitted slope is the panel's distortion and nothing else, and unlike a printed equation it does not presume the draftsman drew the equation faithfully. Upper frame `+0.00078`, lower frame `-0.00036`, neither monotonic, i.e. picking noise. The drawn `R` line rises at `+0.00857` -- **23x the distortion of the panel it sits in, and in the opposite direction**. No page skew lifts `R` by 2.7% while leaving its own frame flat to 0.04%. Independently, `G_bar` in the upper panel reproduces its printed `4.5 (e+)^0.28` to 0.25% RMS off the same scan. Filed as `validation/cooling/data/han2012/fig4.46_frame_*.csv`. |
| 2026-09-10 | extracted by Claude | UNCONFIRMED -- submitted for review |
| 2026-09-10 | reviewer | **CONFIRMED.** All items checked against the book or resolved by derivation; no open flags. D1 and D2 accepted, D3 withdrawn. Released for implementation under #334. |
| 2026-09-10 | reviewer | re-read the layout: the trailing `* f_bar` is a full stop plus the next sentence's subject on the same line. Han prints `f = f_bar + (H/W)(f_bar - f_s)`, matching the derivation. **D3 withdrawn** -- no departure from the source exists. |
| 2026-09-10 | reviewer | closed the last flags. 24 and 26 confirmed against the book. **Item 6 resolved by derivation**: supplied Han's area-weighted assumption, from which `f = f_bar + (H/W)(f_bar - f_s)` follows exactly -- the printed trailing `* f_bar` is a typo (decision D3). Item 35 decided in favour of the figure's `12.31` (D2). Item 27b: the name comes from Gee & Webb (1980), not Webb & Eckert (1972). **No open items remain.** |
| 2026-09-10 | reviewer | supplied the Webb, Eckert & Goldstein paper. **Items 10 and 25 resolved together**: their Eq. (2) gives `g_bar = G Pr^(-0.57)`, and `3.7/0.703^0.57 = 4.5231` reproduces Fig. 4.46's dashed `4.5` to 0.51%. The dashed line is the Prandtl-normalised `G`, not a four-wall average, and the Pr scaling is already inside Han's figure rather than borrowed. |
| 2026-09-10 | reviewer | proposed Webb, Eckert & Goldstein (1972), IJHMT 15(1), 180-184 as the Prandtl basis. Better fit than Dipprey & Sabersky: same roughness family (repeated ribs) and same `R`/`G` formalism. The book does not cite it, so adopting it is a modelling decision, not an extraction. Opens item 27, a citation collision with the Webb & Eckert 1972 already cited for `thermal_performance_factor`. |
| 2026-09-10 | reviewer | item 25: confirmed the text presents the correlations for `Pr ~ 0.7`. Raised Dipprey and Sabersky (1963) as a possible source of a `Pr`-explicit form; its role in the book is at Eq. 4.14, as the conceptual origin of the wall laws. Not yet consulted. |
| 2026-09-10 | reviewer | accepted D1: follow Han in treating `R` as constant in `e+`. The figure does show a slope, but it is smaller than Han's own 6% band and following the text keeps the chain closed form. |
| 2026-09-10 | extracted by Claude | systematic comparison of both channels. Five quantities agree exactly. Flag list raised: items 24 (`(W/H)^m` in Eq. 4.17), 6, 10, 25 (`Pr = 0.703` only) and 26 (Eq. 4.14 has one channel). |
| 2026-09-10 | reviewer | section text supplied (`sources/heat_transfer/han/han_ribs.md`). **Item 7 resolved**: the text states R is independent of `e+`, valid `e+ >= 50`, 6% deviation for 95% of data -- the measured 4.5% slope was artefact. **Item 23 resolved**: Fig. 4.46 and Fig. 4.47 cover different configurations. **Item 10 partly resolved.** Adds Eq. 4.14, the stated tolerances, and opens item 24 -- the text extraction dropped `(W/H)^m` from Eq. 4.17 where the image has it. |
| 2026-09-10 | reviewer | Eq. 4.18 image supplied: item 19 confirmed, items 20-22 now interpretable. Opened item 23 -- Eq. 4.18 and Fig. 4.46 disagree on `G` by up to 22%, with different `e+` exponents (0.35 vs 0.28). |
| 2026-09-10 | reviewer | page 377 supplied: resolves the `m` rule (item 15), adds the `W/H` cap (16), the validity range (17) and the Eq. 4.18 exponents (20-22). Eq. 4.18 itself still missing -- plain text drops equation images. |
| 2026-09-10 | reviewer | item 7 challenged: y-axis is logarithmic and the curve is not flat. Measured against the panel frame as a skew control: curve slope is ~100x the frame's, about 4.5% rise across the range. Item 7 reopened; the closed-form claim now conditional. |

Change **Status** at the top of this file when reviewed, and record corrections
here rather than silently editing the tables above -- a correction is evidence
about the extraction process, not just about the number.
